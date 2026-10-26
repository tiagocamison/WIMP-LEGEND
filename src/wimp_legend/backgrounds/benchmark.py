"""BG-BENCHMARK-001: replaceable poster components, not a LEGEND prediction.

Activities precede all detector efficiencies. The flat component is already
observed. Response ports are supplied by the caller; this module owns no detector
configuration. Nuclear/solar inputs are versioned, separately sourced snapshots.
"""
from dataclasses import dataclass
from importlib.resources import files
from math import fsum
from types import MappingProxyType
from typing import Mapping
import hashlib,json
import numpy as np
from numpy.polynomial.legendre import leggauss

from .._analysis import scalar,text
from ..constants import SECONDS_PER_DAY
from ..energy import EnergyDomain
from .beta import AllowedBetaShape,TabulatedBetaShape,M_E_KEV
from .model import CompositeBackground,FlatBackground
from .routing import TruthBackground
from .radioactive import (BetaBranch,CaptureShell,ElectronCaptureBranch,DecayChain,
                           mean_activity,probability)
from .neutrinos import FluxQuadrature,SolarSurvivalBenchmark,SolarScatteringSpectrum


_INPUT_MANIFEST_SHA256 = "a2b295083c73caf7d3d542ff24911967d0bf712609a5bd3f54ba53bf0683af20"


def benchmark_inputs():
    """Read and check the versioned source files; return independent dictionaries."""
    root=files('wimp_legend.backgrounds').joinpath('data','benchmark_v1')
    manifest_bytes=root.joinpath('manifest.json').read_bytes()
    if hashlib.sha256(manifest_bytes).hexdigest()!=_INPUT_MANIFEST_SHA256:
        raise ValueError('background source manifest checksum mismatch')
    manifest=json.loads(manifest_bytes)
    for name,digest in manifest['files'].items():
        if hashlib.sha256(root.joinpath(name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'background source checksum mismatch: {name}')
    nuclear=json.loads(root.joinpath('nuclear.json').read_text())
    solar=json.loads(root.joinpath('solar.json').read_text())
    # Original source rows are retained beside this portable central spectrum.
    raw=np.array([[float(v) for v in line.split(',')] for line in
                  root.joinpath('ar39.csv').read_text().splitlines()[1:]])
    ar=TabulatedBetaShape(tuple(raw[:,0]),tuple(raw[:,1]),
        root.joinpath('ar39_provenance.json').read_text())
    if ar.q_keV!=nuclear['decays']['Ar39']['endpoint_keV']:
        raise ValueError('Ar39 table does not cover the complete sourced endpoint')
    return nuclear,solar,ar,manifest


@dataclass(frozen=True)
class BenchmarkConfig:
    label: str
    activities_per_kg_day: tuple[tuple[str,float],...]
    flat_rate_per_kg_day_keV_obs: float
    ar39_activity_Bq_per_kg_Ar: float
    argon_mass_per_kg_Ge: float
    lar_transfer_probability: float
    lar_veto_survival_probability: float
    cooldown_day: float
    run_duration_day: float
    ga68_initialization: str
    ga68_initial_activity_per_kg_day: float | None
    zn65_gamma_escape_probability: float
    ga68_gamma_escape_probability: float
    annihilation_escape_probability: float
    solar_flux_components: tuple[str,...]
    nu_target_A: float
    nu_target_mass_GeV: float
    nu_targets_per_kg: float
    solar_sin2_theta12: float
    solar_sin2_theta13: float
    solar_delta_m21_eV2: float
    solar_rho_Ye_g_cm3: float
    provenance: str

    def __post_init__(self):
        for name in ('label','provenance'):object.__setattr__(self,name,text(getattr(self,name),name))
        pairs=tuple((text(k,'isotope'),scalar(v,'activity_per_kg_day')) for k,v in self.activities_per_kg_day)
        if len(pairs)!=4 or {k for k,v in pairs}!={'H3','Ge68','Zn65','Fe55'}:
            raise ValueError('supply H3/Ge68/Zn65/Fe55 exactly once, in decays/(kg Ge day)')
        object.__setattr__(self,'activities_per_kg_day',pairs)
        positive={'run_duration_day','nu_target_A','nu_target_mass_GeV','nu_targets_per_kg','solar_delta_m21_eV2'}
        probs={'lar_transfer_probability','lar_veto_survival_probability','zn65_gamma_escape_probability',
               'ga68_gamma_escape_probability','annihilation_escape_probability','solar_sin2_theta12','solar_sin2_theta13'}
        for name in ('flat_rate_per_kg_day_keV_obs','ar39_activity_Bq_per_kg_Ar','argon_mass_per_kg_Ge',
                     'lar_transfer_probability','lar_veto_survival_probability','cooldown_day','run_duration_day',
                     'zn65_gamma_escape_probability','ga68_gamma_escape_probability','annihilation_escape_probability',
                     'nu_target_A','nu_target_mass_GeV','nu_targets_per_kg','solar_sin2_theta12',
                     'solar_sin2_theta13','solar_delta_m21_eV2','solar_rho_Ye_g_cm3'):
            value=probability(getattr(self,name),name) if name in probs else scalar(getattr(self,name),name,positive=name in positive)
            object.__setattr__(self,name,value)
        if self.nu_target_A<32:raise ValueError('Ge neutrino benchmark requires A>=32')
        if self.ga68_initialization not in ('zero','secular_equilibrium','explicit'):
            raise ValueError('unknown Ga68 initialization')
        if self.ga68_initialization=='explicit':
            object.__setattr__(self,'ga68_initial_activity_per_kg_day',scalar(self.ga68_initial_activity_per_kg_day,'Ga68 initial activity'))
        elif self.ga68_initial_activity_per_kg_day is not None:
            raise ValueError('an explicit Ga68 initial activity requires explicit initialization')
        names=tuple(self.solar_flux_components)
        if len(set(names))!=len(names) or any(n not in ('pp','B8') for n in names):
            raise ValueError('only pinned pp/B8 spectra supported; replace the component for other fluxes')
        object.__setattr__(self,'solar_flux_components',names)


@dataclass(frozen=True)
class BenchmarkModel:
    components: Mapping
    config: BenchmarkConfig
    provenance: str

    def __post_init__(self):
        components=dict(self.components)
        # Existing composite performs the standard port validation.
        CompositeBackground(tuple(components.values()))
        object.__setattr__(self,'components',MappingProxyType(components))
        object.__setattr__(self,'provenance',text(self.provenance,'provenance'))

    @property
    def composite(self):return CompositeBackground(tuple(self.components.values()))

    def expected_counts(self, observed_bin_edges_keV, exposure_kg_day):
        return self.composite.expected_counts(observed_bin_edges_keV,exposure_kg_day)


def solar_fluxes(solar, names, order=4):
    """Gauss-Legendre quadrature of the pinned piecewise-linear full spectra."""
    x,w=leggauss(order);result=[]
    for record in solar['components']:
        if record['name'] not in names:continue
        e=np.asarray(record['energies_MeV']);p=np.asarray(record['density_per_MeV'])
        nodes=(e[:-1,None]+e[1:,None])/2+np.diff(e)[:,None]*x/2
        density=p[:-1,None]+(p[1:]-p[:-1])[:,None]*(nodes-e[:-1,None])/np.diff(e)[:,None]
        weights=record['total_flux_cm2_s']*np.diff(e)[:,None]*w*density/2
        result.append(FluxQuadrature(tuple(nodes.ravel()),tuple(weights.ravel()),json.dumps(record,sort_keys=True)))
    return tuple(result)


def _decay_background(record, activity, response, *, gamma_escape, annihilation_escape):
    raw_total=fsum(r['fraction'] for r in record['ec_branches']+record['beta_plus_branches'])
    if not 0.999<raw_total<1.001:raise ValueError('decay branches are not a complete rounded scheme')
    provenance=json.dumps(dict(source=record,mean_activity_per_kg_Ge_day=activity,branch_normalization_factor=1/raw_total,
        transport='all-or-nothing nuclear gamma/annihilation containment; no partial Compton deposit model'),sort_keys=True)
    lines=[];beta=[]
    for row in record['ec_branches']:
        shells=tuple(CaptureShell(name,energy,p) for name,energy,p in zip(
            ('K','L','M+'),record['vacancy_energies_keV'],(row['PK'],row['PL'],row['PMplus'])))
        branch=ElectronCaptureBranch(row['fraction']/raw_total,shells,row['excitation_keV'],provenance)
        lines.extend(branch.lines(activity,gamma_escape_probability=gamma_escape))
    daughter_z={'Cu':29,'Zn':30}[record['daughter']] if record['beta_plus_branches'] else 0
    A=65 if record['daughter']=='Cu' else 68
    for row in record['beta_plus_branches']:
        shape=AllowedBetaShape(row['endpoint_keV'],daughter_z,A,beta_minus=False)
        gammas=((0.,1.),) if row['excitation_keV']==0 else ((0.,gamma_escape),(row['excitation_keV'],1-gamma_escape))
        for energy,gweight in gammas:
            for annihilation,aweight in ((0.,annihilation_escape),(2*M_E_KEV,1-annihilation_escape)):
                weight=row['fraction']/raw_total*gweight*aweight
                if weight:beta.append(BetaBranch(shape,row['endpoint_keV'],weight,energy+annihilation,provenance))
    def spectrum(e):return activity*sum((np.asarray(branch(e)) for branch in beta),np.zeros_like(np.asarray(e),dtype=float))
    if beta:
        interval=(min(b.deposited_offset_keV for b in beta),max(b.deposited_offset_keV+b.endpoint_keV for b in beta))
        points=tuple(sorted({b.deposited_offset_keV for b in beta}|{b.deposited_offset_keV+b.endpoint_keV for b in beta}))
        return TruthBackground(EnergyDomain.ER,response,spectrum,interval,tuple(lines),points,provenance)
    return TruthBackground(EnergyDomain.ER,response,lines=tuple(lines),provenance=provenance)


def build_benchmark(config: BenchmarkConfig, *, er_response, nr_response, component_responses=None):
    """Assemble existing routed components; optional response overrides by name.

    Override keys: H3,Ge68,Ga68,Zn65,Fe55,Ar39,solar_CEvNS,solar_ES.
    The flat term is observed and receives no response/efficiency a second time.
    """
    if not isinstance(config,BenchmarkConfig):raise TypeError('BenchmarkConfig required')
    overrides=dict(component_responses or {})
    allowed={'H3','Ge68','Ga68','Zn65','Fe55','Ar39','solar_CEvNS','solar_ES'}
    if set(overrides)-allowed:raise ValueError('unknown background response override')
    def response(name,domain):return overrides.get(name,nr_response if domain is EnergyDomain.NR else er_response)
    nuclear,solar,ar39,manifest=benchmark_inputs();D=nuclear['decays'];initial=dict(config.activities_per_kg_day)
    time=dict(start_day=config.cooldown_day,duration_day=config.run_duration_day)
    averages={name:mean_activity(value,D[name]['half_life_day'],**time) for name,value in initial.items()}
    ga0=(0 if config.ga68_initialization=='zero' else initial['Ge68']) if config.ga68_initialization!='explicit' else config.ga68_initial_activity_per_kg_day
    chain=DecayChain(D['Ge68']['half_life_day'],D['Ga68']['half_life_day'])
    averages['Ge68'],averages['Ga68']=chain.mean_activities(initial['Ge68'],ga0,**time)
    components={}
    h3=AllowedBetaShape(D['H3']['endpoint_keV'],2,3)
    components['H3']=TruthBackground(EnergyDomain.ER,response('H3',EnergyDomain.ER),
        lambda e:averages['H3']*np.asarray(h3(e)),(0,h3.q_keV),provenance=json.dumps(D['H3']))
    for name in ('Ge68','Ga68','Zn65','Fe55'):
        escape=config.zn65_gamma_escape_probability if name=='Zn65' else config.ga68_gamma_escape_probability
        components[name]=_decay_background(D[name],averages[name],response(name,EnergyDomain.ER),
            gamma_escape=escape,annihilation_escape=config.annihilation_escape_probability)
    ar_initial=config.ar39_activity_Bq_per_kg_Ar*SECONDS_PER_DAY*config.argon_mass_per_kg_Ge
    ar_average=mean_activity(ar_initial,D['Ar39']['half_life_day'],**time)
    # Ge transfer and veto selection are explicitly separate from activity.
    ar_intensity=ar_average*config.lar_transfer_probability*config.lar_veto_survival_probability
    components['Ar39']=TruthBackground(EnergyDomain.ER,response('Ar39',EnergyDomain.ER),
        lambda e:ar_intensity*np.asarray(ar39(e)),(0,ar39.q_keV),
        spectrum_breakpoints_keV=ar39.energies_keV,
        provenance=json.dumps(dict(spectrum_source=json.loads(ar39.source_id),
            physical_mean_decays_per_kg_Ge_day_before_transfer=ar_average,
            transfer_probability=config.lar_transfer_probability,veto_survival_probability=config.lar_veto_survival_probability,
            transport_assumption='conditional Ge energy spectrum has the emitted beta shape: toy identity transfer, not transport MC'),sort_keys=True))
    fluxes=solar_fluxes(solar,config.solar_flux_components)
    if fluxes:
        weak=nuclear['weak_inputs']
        survival=SolarSurvivalBenchmark(config.solar_sin2_theta12,config.solar_sin2_theta13,
            config.solar_delta_m21_eV2,config.solar_rho_Ye_g_cm3,weak['G_F_GeV_minus2'],weak['avogadro_per_mol'])
        for process,domain in (('CEvNS',EnergyDomain.NR),('ES',EnergyDomain.ER)):
            spec=SolarScatteringSpectrum(process,fluxes,config.nu_target_A,32,config.nu_target_mass_GeV,
                config.nu_targets_per_kg,weak['G_F_GeV_minus2'],weak['sin2_theta_W'],weak['m_e_GeV'],survival,
                tuple((b,int(n)) for b,n in nuclear['ge_electron_shells_keV_occupancy']),
                solar['coverage']+'; effective Ge target and adiabatic density benchmark; binding stepping for ES')
            maximum=max(max(f.energies_MeV) for f in fluxes)*1e-3
            upper=2*maximum*maximum/config.nu_target_mass_GeV*1e6 if process=='CEvNS' else 2*maximum*maximum/(weak['m_e_GeV']+2*maximum)*1e6
            # pp NR recoils occupy only a few eV: declare their endpoint and
            # subintervals so generic quadrature cannot miss the narrow support.
            endpoints=[2*(max(f.energies_MeV)*1e-3)**2/config.nu_target_mass_GeV*1e6
                       if process=='CEvNS' else 2*(max(f.energies_MeV)*1e-3)**2/(weak['m_e_GeV']+2*max(f.energies_MeV)*1e-3)*1e6
                       for f in fluxes]
            points=tuple(sorted({value*fraction for value in endpoints
                                for fraction in (.001,.01,.05,.1,.25,.5,1.)}))
            # Every fixed flux node has its own recoil cutoff in both channels.
            # Expose all of them, including CEvNS derivative kinks, so default
            # public composite counts do not rely on the plotting-only fold.
            recoils=[]
            for flux in fluxes:
                e=np.asarray(flux.energies_MeV)*1e-3
                recoils.extend(2*e*e/config.nu_target_mass_GeV*1e6 if process=='CEvNS'
                               else 2*e*e/(weak['m_e_GeV']+2*e)*1e6)
            points=tuple(sorted(set(points+tuple(recoils))))
            if process=='ES':points+=tuple(b for b,n in nuclear['ge_electron_shells_keV_occupancy'])
            name='solar_'+process
            components[name]=TruthBackground(domain,response(name,domain),spec,(0,upper),
                spectrum_breakpoints_keV=points,provenance=spec.provenance)
    components['flat_observed']=FlatBackground(config.flat_rate_per_kg_day_keV_obs)
    metadata=json.dumps(dict(status='BENCHMARK, not validated LEGEND-1000 prediction',source_manifest=manifest,
        source_activities_at_reference_per_kg_day=initial,mean_intrinsic_activities_per_kg_day=averages,
        mean_ar39_decays_per_kg_Ge_day_before_transfer=ar_average,ar39_intensity_after_transfer_and_veto=ar_intensity,
        solar_coverage=solar['coverage']),sort_keys=True)
    return BenchmarkModel(components,config,metadata)
