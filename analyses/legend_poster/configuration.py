"""The saved master JSON is the sole working physical configuration."""
from copy import deepcopy
import json
from pathlib import Path
import math
from wimp_legend.interactions import NREFT_CONVENTION_ID
from .models import target, responses, build_background, COMPONENTS

DEFAULT_CONFIG = Path(__file__).with_name('config.json')


def _positive(value, name, zero=False):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
        raise ValueError(name+' must be a finite real number')
    if value<0 or (not zero and value==0):raise ValueError(name+' outside allowed domain')


def validate_config(config):
    c=deepcopy(config)
    try:
        if c['schema_version']!=1:raise ValueError('unsupported config schema')
        if set(c['target']['isotopes'])!={'70','72','73','74','76'}:
            raise ValueError('supply all five Ge isotope property records')
        if c['target']['kg_normalization']!='neutral_atomic_detector_mass':
            raise ValueError('explicit neutral atomic kg normalization required')
        if c['target']['mass_fraction_definition']!='bare_nuclear_ensemble (TargetComposition convention)':
            raise ValueError('mass-fraction semantics must match TargetComposition')
        conv=c['target']['mass_conversion']
        for v in ('u_GeV','electron_mass_GeV','electron_binding_eV'):_positive(conv[v],v)
        for a,row in c['target']['isotopes'].items():
            expected=row['atomic_mass_u']*conv['u_GeV']-32*conv['electron_mass_GeV']+conv['electron_binding_eV']*1e-9
            if not math.isclose(expected,row['mass_GeV'],rel_tol=1e-13):
                raise ValueError('mass conversion inconsistent for Ge'+a)
            if not row['source']['reference'] or not row['source']['url']:raise ValueError('isotope source required')
        for name,record in c['target']['presets'].items():
            if set(record['fractions'])!=set(c['target']['isotopes']):raise ValueError('five explicit fractions required')
            if not record['provenance'] or record['status'] not in ('hypothetical','literature_snapshot','measured'):
                raise ValueError('composition status/provenance required')
            target(c,name)
        for name,d in c['detector']['presets'].items():
            for k in ('threshold_keV_obs','noise_sigma_keV','fano_factor','pair_energy_keV','lindhard_k'):_positive(d[k],k)
            _positive(d['efficiency'],'efficiency',zero=True)
            if d['efficiency']>1 or d['Z']!=32:raise ValueError('invalid Ge detector')
            for k in ('NR_energy_range_keV','ER_energy_range_keV'):
                lo,hi=d[k];_positive(lo,k,zero=True);_positive(hi,k)
                if lo!=0 or hi<=lo:raise ValueError('truth ranges must start at zero')
            responses(c,name)
        s=c['signal']
        if s.get('si_rate_engine','generic') not in ('generic','factorized_shm'):
            raise ValueError('si_rate_engine must be generic or factorized_shm')
        if s['convention_id']!=NREFT_CONVENTION_ID or s['coefficient_basis']!='proton_neutron' or s['coefficient_unit']!='GeV^-2':
            raise ValueError('unsupported Wilson convention; no automatic conversion')
        if not s['operators'] or any(type(op)!=int or op not in (1,8,11) for op in s['operators']):
            raise ValueError('workbench supports single O1/O8/O11 only')
        if len(set(s['operators']))!=len(s['operators']):raise ValueError('duplicate operator')
        for k in ('sigma_p_ref_cm2','coefficient_ref_GeV_minus2','m_proton_GeV','m_N_reference_GeV'):_positive(s[k],k)
        _positive(s['j_chi'],'j_chi',zero=True)
        if 2*s['j_chi']!=int(2*s['j_chi']):raise ValueError('spin must be a half integer')
        if type(s['proton_sign'])!=int or s['proton_sign'] not in (-1,1):raise ValueError('explicit proton sign required')
        if isinstance(s['fn_over_fp'],bool) or not math.isfinite(s['fn_over_fp']):raise ValueError('invalid neutron coupling ratio')
        for m in s['masses_GeV']:_positive(m,'WIMP mass')
        if not s['masses_GeV'] or len(set(s['masses_GeV']))!=len(s['masses_GeV']):raise ValueError('nonempty distinct mass grid required')
        from wimp_legend.halo.shm import StandardHaloModel
        h=s['halo']
        for k in ('v0_km_s','vesc_km_s','v_lab_km_s'):_positive(h[k],k,zero=k=='v_lab_km_s')
        StandardHaloModel.from_km_s(**{k:h[k] for k in ('v0_km_s','vesc_km_s','v_lab_km_s')})
        _positive(h['rho_chi_GeV_cm3'],'halo density')
        a=c['analysis']
        if isinstance(a['confidence_level'],bool) or isinstance(a['roi']['signal_fraction'],bool):raise ValueError('CL/ROI cannot be boolean')
        if not 0<a['confidence_level']<1 or not 0<a['roi']['signal_fraction']<=1:raise ValueError('invalid CL or ROI fraction')
        if a['binned_roi']!='full_analysis_window':raise ValueError('binned ROI must be the full analysis window')
        if a['roi']['method']!='threshold_anchored':raise ValueError('unsupported ROI definition')
        if a['roi']['normalization']!='signal inside configured observed analysis window':raise ValueError('ROI normalization must be explicit')
        if not a['methods'] or not set(a['methods'])<= {'poisson_cls_counting','binned_asimov'}:raise ValueError('unknown statistical method')
        if 'binned_asimov' in a['methods'] and a['confidence_level']<=.5:raise ValueError('asymptotic CL must exceed .5')
        for x in a['exposures_kg_day']:_positive(x,'exposure_kg_day')
        if not a['exposures_kg_day']:raise ValueError('exposure list is empty')
        _positive(a['days_per_year'],'days_per_year')
        n=c['numerical'];previous=None
        _positive(n['gaussian_tail_probability_max'],'Gaussian tail probability bound')
        if math.erfc(n['gaussian_tail_sigma']/math.sqrt(2))>n['gaussian_tail_probability_max']:
            raise ValueError('Gaussian tail truncation exceeds requested probability bound')
        for lo,hi,step in n['bin_segments']:
            _positive(lo,'bin lower',zero=True);_positive(hi,'bin upper');_positive(step,'bin step')
            if hi<=lo or (previous is not None and lo!=previous):raise ValueError('bin segments must be increasing and contiguous')
            previous=hi
        for key in ('background_step_keV','signal_fold_step_keV','signal_grid_min_keV','gaussian_tail_sigma'):_positive(n[key],key)
        for key in ('rate_epsrel','signal_total_rtol','signal_bin_rtol','background_total_rtol','background_bin_rtol','significant_bin_fraction'):
            _positive(n[key],key)
            if n[key]>=1:raise ValueError(key+' must be <1')
        for key in ('background_order','signal_fold_order','signal_grid_points','signal_refinements','blas_threads'):
            if type(n[key])!=int or n[key]<1:raise ValueError('positive integer '+key+' required')
        if n['signal_grid_points']<8 or not 2<=n['background_order']<=32 or not 2<=n['signal_fold_order']<=32:
            raise ValueError('unsupported quadrature/grid size')
        for scenario,record in c['background']['scenarios'].items():
            target(c,record['target'])
            if set(record['components'])-set(COMPONENTS):raise ValueError('unknown component override')
            for k,v in record['components'].items():
                if set(v)-{'enabled','scale','response'}:raise ValueError('unknown component option')
                if 'enabled' in v and type(v['enabled'])!=bool:raise ValueError('enabled must be bool')
                if 'scale' in v:_positive(v['scale'],'component scale',zero=True)
                if k=='flat_observed' and 'response' in v:raise ValueError('observed residual cannot be refolded')
                if 'response' in v:responses(c,v['response'])
            build_background(c,scenario)
            d=c['detector']['presets'][record.get('detector',c['detector']['default'])]
            if not n['bin_segments'][0][0]<=d['threshold_keV_obs']<n['bin_segments'][-1][1]:raise ValueError('threshold outside plotted domain')
        if c['background']['solar_target_mode'] not in ('isotope_mixture','effective_isotope'):
            raise ValueError('solar target approximation must be explicit')
        if not set(a['scenarios'])<=set(c['background']['scenarios']) or not a['scenarios']:raise ValueError('unknown/empty scenario selection')
        if c['background']['historical']['pair_with'] not in c['background']['scenarios']:raise ValueError('invalid historical pair')
        if c['execution']['default_profile'] not in c['execution']['profiles']:raise ValueError('unknown default profile')
        allowed=set(n)|{'masses_GeV','scenarios','operators','label'}
        for profile,values in c['execution']['profiles'].items():
            if set(values)-allowed:raise ValueError('unknown profile field')
        style=c['plotting']
        for name in ('energy_range_keV','threshold_zoom_keV','mass_range_GeV','sigma_range_cm2'):
            lo,hi=style[name];_positive(lo,name,zero=name.endswith('keV'));_positive(hi,name)
            if hi<=lo:raise ValueError('invalid plotting range '+name)
        if type(style['dpi'])!=int or style['dpi']<=0:raise ValueError('positive plotting dpi required')
        if style['exposure_comparison_scenario'] not in c['background']['scenarios']:raise ValueError('unknown exposure-plot scenario')
        for index in style['limit_exposure_indices']:
            if type(index)!=int or not -len(a['exposures_kg_day'])<=index<len(a['exposures_kg_day']):
                raise ValueError('plot exposure index outside exposure list')
        available=json.loads((DEFAULT_CONFIG.parent/'data/si_limits/manifest.json').read_text())
        if set(style['comparisons'])-{v['id'] for v in available['curves']}:raise ValueError('unknown comparison dataset')
        json.dumps(c,allow_nan=False)
    except (KeyError,TypeError,OverflowError) as exc:
        raise ValueError('invalid/incomplete master configuration: '+str(exc)) from exc
    return c


def load_config(path=DEFAULT_CONFIG):
    return validate_config(json.loads(Path(path).read_text(encoding='utf-8')))


def save_config(config,path=DEFAULT_CONFIG):
    c=validate_config(config);path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(c,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    temporary.replace(path)
    return path


def resolve_profile(config,profile=None):
    c=validate_config(config);name=profile or c['execution']['default_profile']
    if name not in c['execution']['profiles']:raise ValueError('unknown execution profile '+name)
    values=c['execution']['profiles'][name]
    for k,v in values.items():
        if k in ('masses_GeV','operators'):c['signal'][k]=deepcopy(v)
        elif k=='scenarios':c['analysis'][k]=deepcopy(v)
        elif k!='label':c['numerical'][k]=deepcopy(v)
    c['execution']['resolved_profile']=name
    c['execution']['resolved_label']=values.get('label',name)
    return validate_config(c)
