"""Small orchestration layer over accepted WIMP-LEGEND physics and statistics.

Raw recoil rates come from the existing generic rate engine / SIO1 adapter.
PCHIP interpolates controlled *truth rates*, never coarse observed-bin counts.
Every observed bin is integrated again against the detector probabilities.
"""
from dataclasses import asdict
from functools import partial
from pathlib import Path
from math import sqrt
import csv
import importlib.util
import json
import numpy as np
from scipy.interpolate import PchipInterpolator
from numpy.polynomial.legendre import leggauss
from wimp_legend.backgrounds import TruthBackground, CompositeBackground
from wimp_legend.backgrounds.routing import MonoenergeticLine
from wimp_legend.energy import EnergyDomain
from wimp_legend.halo.shm import StandardHaloModel
from wimp_legend.interactions import WilsonCoefficients
from wimp_legend.kinematics import maximum_recoil_energy_keV
from wimp_legend.rates import NREFTRateConfig, differential_rate_per_kg_day_keV
from wimp_legend.si import SIO1Signal
from wimp_legend.statistics import PoissonCLsCounting, BinnedPoissonProfile
from wimp_legend.sensitivity import AnalysisResult
from wimp_legend.targets import load_germanium_afh2014, germanium_afh2014_momentum
from .configuration import resolve_profile, validate_config
from .folding import prune_gaussian_bins
from .models import target, isotope_weights, responses, build_background
from .storage import ROOT, NumericalCache, source_identity, digest, file_hash, canonical


def _tool(name):
    spec=importlib.util.spec_from_file_location('legend_poster_'+name,
                ROOT/'tools/background_benchmark'/f'{name}.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

_fold = _tool('folding').quadrature_counts


def observed_edges(config):
    pieces=[]
    for lo,hi,step in config['numerical']['bin_segments']:
        # Explicit endpoint insertion avoids arange drift and ragged last bins.
        pieces.extend(np.r_[np.arange(lo,hi,step),hi])
    for record in config['background']['scenarios'].values():
        name=record.get('detector',config['detector']['default'])
        pieces.append(config['detector']['presets'][name]['threshold_keV_obs'])
    return np.unique(np.round(pieces,12))


def _sigma_bound(component):
    response=component.response
    hi=component.spectrum_energy_range_keV[1]
    mean=hi if component.energy_domain is EnergyDomain.ER else response.quenching.mean_energy_keV(hi)
    return float(response.resolution.sigma_keV(mean))


def fold_component(component,bins,numerical,*,signal=False,refined=False):
    prefix='signal_fold' if signal else 'background'
    step=numerical[prefix+'_step_keV']/(2 if refined else 1)
    order=min(32,numerical[prefix+'_order']+(2 if refined else 0))
    return _fold(prune_gaussian_bins(component,numerical['gaussian_tail_sigma']),bins,1.,step_keV=step,order=order,
                 tail_sigma=numerical['gaussian_tail_sigma'],sigma_bound=_sigma_bound)


def convergence(coarse,fine,numerical,kind):
    coarse=np.asarray(coarse);fine=np.asarray(fine)
    if np.any(fine<0) or np.any(~np.isfinite(fine)):raise ValueError('invalid integrated count')
    total=abs(float(coarse.sum()-fine.sum()))/max(float(fine.sum()),1e-300)
    active=fine>float(fine.max(initial=0))*numerical['significant_bin_fraction']
    bins=float(np.max(abs((coarse[active]-fine[active])/fine[active]))) if np.any(active) else 0.
    return dict(total_relative=total,significant_bin_relative=bins,
                passed=bool(total<=numerical[kind+'_total_rtol'] and bins<=numerical[kind+'_bin_rtol']))


def _check_truth_support(component,config,preset):
    if isinstance(component,CompositeBackground):
        for c in component.components:_check_truth_support(c,config,preset)
    elif isinstance(component,TruthBackground):
        d=config['detector']['presets'][preset]
        key='NR_energy_range_keV' if component.energy_domain is EnergyDomain.NR else 'ER_energy_range_keV'
        hi=d[key][1]
        if component.spectrum is not None and component.spectrum_energy_range_keV[1]>hi:
            raise ValueError('increase detector '+key+' to cover source support')
        if any(line.energy_keV>hi for line in component.lines):raise ValueError('line outside assigned truth range')


def background_counts(config,scenario,bins,cache):
    record=config['background']['scenarios'][scenario]
    model=build_background(config,scenario);n=config['numerical']
    inputs=dict(background=config['background']['base'],scenario=record,target=config['target'],
        solar_target_mode=config['background']['solar_target_mode'],detector=config['detector'],
        bins=bins.tolist(),numerical={k:v for k,v in n.items() if k.startswith('background') or k in ('gaussian_tail_sigma','significant_bin_fraction')})
    diagnostics=[];result={}
    for name,component in model.components.items():
        option=record['components'].get(name,{})
        scale=option.get('scale',1.) if option.get('enabled',True) else 0.
        preset=option.get('response',record.get('detector',config['detector']['default']))
        _check_truth_support(component,config,preset)
        def compute():
            if scale==0:return dict(coarse=np.zeros(len(bins)-1),fine=np.zeros(len(bins)-1))
            return dict(coarse=fold_component(component,bins,n)*scale,
                        fine=fold_component(component,bins,n,refined=True)*scale)
        component_inputs=dict(inputs)
        if name.startswith('solar_'):
            combined=dict(config['background']['base']);combined.update(record['overrides'])
            component_inputs=dict(solar={k:v for k,v in combined.items() if k.startswith('solar_')},
                target=config['target']['presets'][record['target']],isotopes=config['target']['isotopes'],
                mass_conversion=config['target']['mass_conversion'],mode=config['background']['solar_target_mode'],
                detector=config['detector']['presets'][preset],option=option,
                bins=bins.tolist(),numerical=inputs['numerical'])
        arrays=cache.arrays('background:'+name,component_inputs,compute)
        check=convergence(arrays['coarse'],arrays['fine'],n,'background')
        if not check['passed']:raise ValueError(f'background convergence failed: {scenario}/{name}: {check}')
        diagnostics.append(dict(scenario=scenario,component=name,kind='background',**check))
        result[name]=arrays['fine']
    result['Total']=sum(result.values(),np.zeros(len(bins)-1))
    return result,diagnostics,json.loads(model.provenance)


def signal_model(config,isotope,mass,operator):
    s=config['signal'];h=s['halo']
    halo=StandardHaloModel.from_km_s(**{k:h[k] for k in ('v0_km_s','vesc_km_s','v_lab_km_s')})
    dataset=load_germanium_afh2014(isotope,momentum=germanium_afh2014_momentum(isotope.A))
    weights=isotope_weights(config,next(iter(config['target']['presets'])))
    common=dict(m_chi_GeV=mass,dataset=dataset,rho_chi_GeV_cm3=h['rho_chi_GeV_cm3'],
        targets_per_kg=weights[isotope.A]['pure_targets_per_kg'],j_chi=s['j_chi'],
        config=NREFTRateConfig(s['m_N_reference_GeV']),speed_pdf=halo.lab_speed_pdf,
        beta_max=halo.max_lab_speed,breakpoints=halo.integration_breakpoints,
        epsrel=config['numerical']['rate_epsrel'],epsabs=0.)
    if operator==1:
        spectrum=SIO1Signal(sigma_p_ref_cm2=s['sigma_p_ref_cm2'],m_proton_GeV=s['m_proton_GeV'],
                    fn_over_fp=s['fn_over_fp'],proton_sign=s['proton_sign'],**common)
    else:
        cp=s['proton_sign']*s['coefficient_ref_GeV_minus2']
        coefficients=WilsonCoefficients.from_proton_neutron({operator:(cp,s['fn_over_fp']*cp)})
        spectrum=partial(differential_rate_per_kg_day_keV,coefficients=coefficients,**common)
    endpoint=float(maximum_recoil_energy_keV(mass,isotope.mass_GeV,halo.max_lab_speed))
    return spectrum,endpoint


def pure_signal_counts(config,isotope,mass,operator,preset,bins,cache):
    n=config['numerical'];_,response=responses(config,preset)
    spectrum,endpoint=signal_model(config,isotope,mass,operator)
    d=config['detector']['presets'][preset]
    if endpoint>d['NR_energy_range_keV'][1]:raise ValueError('signal endpoint outside configured NR range')
    s=config['signal'];inputs=dict(isotope=config['target']['isotopes'][str(isotope.A)],
        mass_conversion=config['target']['mass_conversion'],signal={k:v for k,v in s.items() if k not in ('masses_GeV','operators')},
        mass=mass,operator=operator,detector=d,bins=bins.tolist(),numerical=n)
    def calculate(points):
        grid_inputs=dict(inputs,grid_points=points)
        truth_inputs={k:inputs[k] for k in ('isotope','mass_conversion','signal','mass','operator')}
        truth_inputs.update(grid_points=points,rate_epsrel=n['rate_epsrel'],grid_min_keV=n['signal_grid_min_keV'])
        def evaluate():
            grid=np.unique(np.r_[0.,endpoint*np.linspace(0,1,points)**2,
                min(n['signal_grid_min_keV'],endpoint),endpoint])
            values=np.asarray(spectrum(grid),dtype=float)
            if np.any(values<0):raise ValueError('negative core nuclear rate; do not clip source diffraction-node artifacts')
            return dict(energy=grid,rate=values)
        truth=cache.arrays('truth-rate',truth_inputs,evaluate)
        interp=PchipInterpolator(truth['energy'],truth['rate'],extrapolate=False)
        def kernel(energy):
            e=np.asarray(energy);out=np.zeros_like(e,dtype=float);active=(e>=0)&(e<=endpoint)
            out[active]=interp(e[active]);return out
        component=TruthBackground(EnergyDomain.NR,response,kernel,(0.,endpoint),
            spectrum_breakpoints_keV=tuple(truth['energy']),provenance='Interpolated accepted core rate; numerically refined, no observed-count interpolation')
        def integrate():return dict(counts=fold_component(component,bins,n,signal=True,refined=True),
                                    coarse=fold_component(component,bins,n,signal=True))
        folded=cache.arrays('folded-signal',grid_inputs,integrate)
        fold_check=convergence(folded['coarse'],folded['counts'],n,'signal')
        if not fold_check['passed']:raise ValueError('signal folding quadrature failed convergence; refine fold step/order')
        return folded['counts'],fold_check
    points=n['signal_grid_points'];coarse,_=calculate(points)
    check=None
    for refinement in range(n['signal_refinements']):
        points=2*points-1;fine,fold_check=calculate(points)
        check=convergence(coarse,fine,n,'signal')
        if check['passed']:
            return fine,dict(kind='signal',isotope=isotope.A,operator=operator,mass_GeV=mass,
                             grid_points=points,fold_total_relative=fold_check['total_relative'],
                             fold_bin_relative=fold_check['significant_bin_relative'],**check)
        coarse=fine
    raise ValueError(f'signal interpolation convergence failed: Ge{isotope.A}/O{operator}/{mass}GeV: {check}; increase grid points/refinements')


def historical_counts(config,bins,cache):
    h=config['background']['historical'];path=ROOT/h['notebook']
    if file_hash(path)!=h['sha256']:raise ValueError('historical notebook checksum mismatch')
    def compute():
        function,_=_tool('legacy').load_legacy(path)
        x,w=leggauss(16)
        nodes=((bins[:-1,None]+bins[1:,None])/2+np.diff(bins)[:,None]*x/2).ravel()
        weights=(np.diff(bins)[:,None]*w/2).ravel()
        values=function(nodes,target=h['target'],argon=h['argon'])
        return {k:np.sum((np.asarray(v)*weights).reshape(len(bins)-1,-1),axis=1)/config['analysis']['days_per_year']
                for k,v in values.items()}
    return cache.arrays('historical-observed-background',dict(historical=h,bins=bins.tolist(),days_per_year=config['analysis']['days_per_year']),compute)


def roi_slice(signal,bins,threshold,fraction):
    first=int(np.searchsorted(bins,round(threshold,12)))
    if first>=len(signal):raise ValueError('threshold outside bins')
    cumulative=np.cumsum(signal[first:]);total=float(cumulative[-1])
    if total<=0:return slice(first,len(signal)),0.
    stop=first+int(np.searchsorted(cumulative,total*fraction))+1
    return slice(first,stop),float(signal[first:stop].sum()/total)


def count_limit(signal,background,bins,exposure,threshold,config,method):
    region,fraction=roi_slice(signal,bins,threshold,config['analysis']['roi']['signal_fraction'])
    if method=='binned_asimov':
        region=slice(int(np.searchsorted(bins,round(threshold,12))),len(signal));fraction=1. if signal[region].sum()>0 else 0.
    if method=='poisson_cls_counting':
        s=np.array([signal[region].sum()*exposure]);b=np.array([background[region].sum()*exposure])
        out_edges=[float(bins[region.start]),float(bins[region.stop])]
        backend=PoissonCLsCounting()
    else:
        s=signal[region]*exposure;b=background[region]*exposure
        out_edges=bins[region.start:region.stop+1];backend=BinnedPoissonProfile()
    result=backend.expected_upper_limit(s,b,confidence_level=config['analysis']['confidence_level'])
    analysis=AnalysisResult(tuple(s),tuple(b),tuple(out_edges),exposure,(0.,config['detector']['presets'][config['detector']['default']]['NR_energy_range_keV'][1]),result)
    return analysis,dict(roi_lower_keV=float(out_edges[0]),roi_upper_keV=float(out_edges[-1]),
        retained_signal_fraction=fraction,reference_signal_events=float(sum(s)),background_events=float(sum(b)),
        mu_upper=result.mu_upper,status=result.status,confidence_level=result.confidence_level,
        expected_kind=result.kind,pseudo_observed_count=float(sum(result.reference_counts)),
        statistical_notes=json.dumps(result.notes),statistical_construction=result.construction,
        statistical_approximation=result.approximation)


def run_analysis(config,*,profile=None,cache_dir=None,output_dir=None):
    from threadpoolctl import threadpool_limits
    c=resolve_profile(config,profile)
    with threadpool_limits(limits=c['numerical']['blas_threads']):
        return _execute_analysis(c,cache_dir=cache_dir,output_dir=output_dir)


def _execute_analysis(c,*,cache_dir=None,output_dir=None):
    identity=source_identity();bins=observed_edges(c)
    cache=NumericalCache(cache_dir or Path(__file__).parent/'cache',identity['physics_fingerprint'])
    numerical_id=digest(dict(config={k:v for k,v in c.items() if k!='plotting'},physics=identity['physics_fingerprint']))
    run_id=digest(dict(numerical_id=numerical_id,config=c,git_commit=identity['git_commit'],plot_code=file_hash(Path(__file__).with_name('plotting.py'))))
    path=Path(output_dir or Path(__file__).parent/'runs')/run_id[:20];path.mkdir(parents=True,exist_ok=True)
    backgrounds={};signals={};diagnostics=[];rows=[];assumptions={};weights={}
    pure={}
    for scenario in c['analysis']['scenarios']:
        record=c['background']['scenarios'][scenario];target_id=record['target']
        preset=record.get('detector',c['detector']['default']);d=c['detector']['presets'][preset]
        background,checks,provenance=background_counts(c,scenario,bins,cache)
        diagnostics.extend(checks);assumptions[scenario]=provenance
        for name,values in background.items():backgrounds[scenario+'__'+name]=values
        weights[scenario]=isotope_weights(c,target_id)
        comparisons={'revised':background['Total']}
        if c['background']['historical']['enabled'] and scenario==c['background']['historical']['pair_with']:
            old=historical_counts(c,bins,cache)
            for name,values in old.items():backgrounds[scenario+'__historical__'+name]=values
            comparisons['historical']=old['Total']
        for mass in c['signal']['masses_GeV']:
            for operator in c['signal']['operators']:
                signal=np.zeros(len(bins)-1)
                for isotope,fraction in target(c,target_id):
                    weight=weights[scenario][isotope.A]['bulk_mass_fraction']
                    if weight==0:continue
                    key=(isotope.A,mass,operator,preset)
                    if key not in pure:
                        values,check=pure_signal_counts(c,isotope,mass,operator,preset,bins,cache)
                        pure[key]=values;diagnostics.append(check)
                        signals[f'pure_Ge{isotope.A}__O{operator}__m{mass:g}__{preset}']=values
                    signal+=weight*pure[key]
                signals[f'{scenario}__O{operator}__m{mass:g}']=signal
                for exposure in c['analysis']['exposures_kg_day']:
                    for label,b in comparisons.items():
                        for method in c['analysis']['methods']:
                            _,result=count_limit(signal,b,bins,exposure,d['threshold_keV_obs'],c,method)
                            mu=result['mu_upper']
                            if operator==1:
                                # Use the same public SI adapter conversion as the core.
                                isotope=next(i for i,w in target(c,target_id) if w>0)
                                adapter,_=signal_model(c,isotope,mass,operator)
                                limit=adapter.cross_section_limit_cm2(mu);unit='cm^2'
                            else:
                                limit=sqrt(mu)*abs(c['signal']['coefficient_ref_GeV_minus2']);unit='GeV^-2'
                            rows.append(dict(scenario=scenario,target=target_id,background=label,
                                mass_GeV=mass,operator=operator,exposure_kg_day=exposure,
                                method=method,limit_value=limit,limit_unit=unit,
                                convention_id=c['signal']['convention_id'],**result))
    np.savez_compressed(path/'backgrounds.npz',observed_edges_keV=bins,**backgrounds)
    np.savez_compressed(path/'signals.npz',observed_edges_keV=bins,**signals)
    for name,data in [('limits.csv',rows),('convergence.json',diagnostics)]:
        if name.endswith('.csv'):
            with (path/name).open('w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        else:(path/name).write_text(json.dumps(data,indent=2)+'\n')
    manifest=dict(schema_version=1,run_id=run_id,numerical_id=numerical_id,resolved_configuration=c,
        source=identity,isotope_weights=weights,
        nuclear_masses_GeV={str(i.A):i.mass_GeV for i,_ in target(c,next(iter(c['target']['presets'])))},
        background_provenance=assumptions,
        cache=dict(hits=cache.hits,misses=cache.misses,keys=sorted(set(cache.used))),
        numerical_notes=['All observed bins reintegrated; PCHIP used only for core truth rates.',
            'Gaussian continuum/source and batch-bin pruning at configured sigma bound; combined omitted probability <=2 Phi(-tail_sigma) per event; exact unpruned line probabilities.',
            'Per-kg normalization uses neutral atomic mass; kinematics uses derived bare nuclear masses.',
            'Solar CEvNS/ES summed over configured isotope nuclei when isotope_mixture selected.',
            '95% ROI conditional on signal inside threshold-to-upper analysis window; not full unobserved signal.',
            'Historical backgrounds already observed; no additional folding/efficiency; frozen historical response.'],
        output_sha256={name:file_hash(path/name) for name in ('backgrounds.npz','signals.npz','limits.csv','convergence.json')})
    (path/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    return load_run(path)


def load_run(path):
    path=Path(path);manifest=json.loads((path/'manifest.json').read_text())
    for name,expected in manifest['output_sha256'].items():
        if file_hash(path/name)!=expected:raise ValueError('run output checksum mismatch: '+name)
    arrays={}
    for name in ('backgrounds','signals'):
        with np.load(path/(name+'.npz'),allow_pickle=False) as f:arrays[name]={k:f[k].copy() for k in f.files}
    with (path/'limits.csv').open() as f:
        rows=list(csv.DictReader(f))
    for row in rows:
        for key in ('mass_GeV','exposure_kg_day','limit_value','roi_lower_keV','roi_upper_keV','mu_upper','reference_signal_events','background_events','retained_signal_fraction'):row[key]=float(row[key])
        row['operator']=int(row['operator'])
    return dict(path=path,manifest=manifest,limits=rows,convergence=json.loads((path/'convergence.json').read_text()),**arrays)
