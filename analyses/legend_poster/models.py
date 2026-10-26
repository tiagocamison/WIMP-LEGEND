"""Build existing target/halo/detector/background APIs from resolved input data."""
from dataclasses import replace
import json
import numpy as np
from wimp_legend.constants import C_CM_S, GEV_TO_J
from wimp_legend.targets import Isotope, SourceReference, TargetComposition
from wimp_legend.detector import (DetectorResponse, ElectronRecoilResponse,
    LindhardQuenching, GaussianResolution, ConstantEfficiency)
from wimp_legend.backgrounds import BenchmarkConfig, build_benchmark, CompositeBackground

COMPONENTS = ('H3','Ge68','Ga68','Zn65','Fe55','Ar39','solar_CEvNS','solar_ES','flat_observed')


def target(config, preset):
    data = config['target']
    isotopes = []
    for a, row in data['isotopes'].items():
        mass_source = SourceReference(row['source']['reference'], row['source']['url'],
            'implementation', notes=str(data['mass_conversion']))
        spin_source = SourceReference('AFH PRC89 065501 supplement; NIST Ge atomic handbook',
            'ground-state Ge isotope spins', 'primary_literature')
        isotopes.append(Isotope('Ge', int(a), 32, row['mass_GeV'], row['spin'],
            {'identity':spin_source,'mass_GeV':mass_source,'spin':spin_source}))
    record = data['presets'][preset]
    source = SourceReference(record['provenance'], 'master configuration target preset '+preset,
                             'synthetic' if record['status']=='hypothetical' else 'primary_literature')
    return TargetComposition(tuple((i, record['fractions'][str(i.A)]) for i in isotopes),
                             record['fraction_kind'], source)


def isotope_weights(config, preset):
    """Number, nuclear-mass and neutral-atomic bulk-mass weights, explicitly distinct.

    Mass fractions supplied to TargetComposition always mean nuclear-ensemble
    fractions. Per-kg detector counts use neutral atomic masses from the same
    source, so electron mass is included in the bulk mass but not kinematics.
    """
    comp = target(config, preset)
    numbers = {i.A:x for i,x in comp.to_number_fractions()}
    nuclear = {i.A:w for i,w in comp.to_mass_fractions()}
    conversion = config['target']['mass_conversion']['u_GeV']
    atoms = {int(a):v['atomic_mass_u']*conversion for a,v in config['target']['isotopes'].items()}
    average = sum(numbers[a]*atoms[a] for a in atoms)
    kg_GeV = (C_CM_S/100)**2/GEV_TO_J
    return {a:dict(number_fraction=numbers[a], nuclear_mass_fraction=nuclear[a],
                   bulk_mass_fraction=numbers[a]*atoms[a]/average,
                   targets_per_kg=numbers[a]*kg_GeV/average,
                   pure_targets_per_kg=kg_GeV/atoms[a]) for a in numbers}


def responses(config, preset=None):
    name = preset or config['detector']['default']
    d = config['detector']['presets'][name]
    def sigma(energy):
        return np.sqrt(d['noise_sigma_keV']**2+d['fano_factor']*d['pair_energy_keV']*np.asarray(energy))
    kwargs = dict(resolution=GaussianResolution(sigma),
                  observed_efficiency=ConstantEfficiency(d['efficiency']))
    return (ElectronRecoilResponse(**kwargs),
            DetectorResponse(LindhardQuenching(d['Z'],d['lindhard_k']),**kwargs))


def build_background(config, scenario):
    record = config['background']['scenarios'][scenario]
    preset = record.get('detector',config['detector']['default'])
    er,nr = responses(config,preset)
    comp = target(config,record['target'])
    weights = isotope_weights(config,record['target'])
    selected = [(i,f) for i,f in comp.to_number_fractions() if f>0]
    nt = sum(w['targets_per_kg'] for w in weights.values())
    kwargs = dict(config['background']['base']);kwargs.update(record['overrides'])
    kwargs.update(label=scenario,nu_target_A=sum(i.A*f for i,f in selected),
        nu_target_mass_GeV=sum(i.mass_GeV*f for i,f in selected),nu_targets_per_kg=nt)
    overrides = {k:responses(config,v['response'])[1 if k=='solar_CEvNS' else 0]
                 for k,v in record['components'].items() if 'response' in v}
    model = build_benchmark(BenchmarkConfig(**kwargs),er_response=er,nr_response=nr,
                            component_responses=overrides)
    components = dict(model.components)
    # Delegate isotope-specific solar kernels/endpoints to the existing factory,
    # not a new weak-scattering calculation or effective isotope surrogate.
    if config['background']['solar_target_mode']=='isotope_mixture':
        solar = {k:[] for k in ('solar_CEvNS',) if k in components}
        # ES has the same Z/shell/kernel for every Ge isotope; the factory's
        # total atom count is exactly the mixture sum, not an effective-mass
        # approximation to ES. Only CEvNS needs separate isotope kernels.
        for isotope, _ in selected:
            isotope_cfg = replace(model.config,nu_target_A=float(isotope.A),
                nu_target_mass_GeV=isotope.mass_GeV,
                nu_targets_per_kg=weights[isotope.A]['targets_per_kg'])
            part = build_benchmark(isotope_cfg,er_response=er,nr_response=nr,
                                   component_responses=overrides)
            for k in solar:solar[k].append(part.components[k])
        for k,v in solar.items():components[k]=CompositeBackground(tuple(v))
    provenance=json.loads(model.provenance)
    provenance['workbench_solar_target']=dict(mode=config['background']['solar_target_mode'],
        isotope_weights=weights,
        note='CEvNS explicit isotope kernels; ES uses exact total Ge atom count. Factory effective fields are bookkeeping only in mixture mode.')
    return replace(model,components=components,provenance=json.dumps(provenance,sort_keys=True))
