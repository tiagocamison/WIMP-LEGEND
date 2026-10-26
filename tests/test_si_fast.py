"""Independent adaptive generic NREFT oracles for the optional SI reduction."""
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import sys
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from wimp_legend.halo.shm import StandardHaloModel
from wimp_legend.halo.integrals import velocity_moment
from wimp_legend.halo.shm_eta import mean_inverse_speed_shm
from wimp_legend.interactions import WilsonCoefficients,CANONICAL_OPERATORS
from wimp_legend.si_fast import is_o1_only,differential_rate_shm_per_kg_day_keV
from analyses.legend_poster import load_config,resolve_profile,validate_config,load_run
from analyses.legend_poster.models import target,isotope_weights
from analyses.legend_poster.analysis import signal_model,pure_signal_counts,observed_edges,count_limit
from analyses.legend_poster.storage import NumericalCache,source_identity


@pytest.fixture(scope='module')
def config():return load_config()


@pytest.mark.parametrize('lab',[0.,1e-8,266.,544.,700.])
def test_eta_against_adaptive_moment_and_support(lab):
    h=StandardHaloModel.from_km_s(v_lab_km_s=lab)
    a=np.unique(np.r_[0.,np.linspace(0,h.max_lab_speed,25),h.integration_breakpoints,
                       h.max_lab_speed*1.1])
    a=a[(a<h.max_lab_speed*.999)|(a>=h.max_lab_speed)]
    reference=velocity_moment(h.lab_speed_pdf,h.max_lab_speed,a,n=-1,
                              breakpoints=h.integration_breakpoints,epsrel=1e-10)
    actual=mean_inverse_speed_shm(h,a)
    np.testing.assert_allclose(actual,reference,rtol=2e-9,atol=reference.max()*2e-14)
    assert np.all(actual>=0) and np.all(np.diff(actual)<=reference.max()*2e-14)
    assert mean_inverse_speed_shm(h,h.max_lab_speed)==0
    assert mean_inverse_speed_shm(h,np.inf)==0
    assert mean_inverse_speed_shm(h,-1)==mean_inverse_speed_shm(h,0)
    assert mean_inverse_speed_shm(h,np.empty((0,3))).shape==(0,3)
    with pytest.raises(ValueError):mean_inverse_speed_shm(h,np.nan)


def test_small_escape_sphere_eta():
    h=StandardHaloModel(v0=.001,vesc=.00001,v_lab=0.)
    a=np.linspace(0,h.max_lab_speed,15)
    expected=velocity_moment(h.lab_speed_pdf,h.max_lab_speed,a,n=-1,epsrel=1e-10)
    np.testing.assert_allclose(mean_inverse_speed_shm(h,a),expected,rtol=1e-12,atol=0)


@pytest.mark.parametrize('lab',[1e-8,266.,700.])
def test_eta_tail_against_decimal_angular_integral(lab):
    # Independent 70-digit exp subtraction in the angularly integrated PDF,
    # not the production erf/expm1 branches or double-precision adaptive oracle.
    from decimal import Decimal as D,localcontext
    h=StandardHaloModel.from_km_s(v_lab_km_s=lab)
    nodes,weights=np.polynomial.legendre.leggauss(32)
    with localcontext() as context:
        context.prec=70
        v0,escape,boost=map(D.from_float,(h.v0,h.vesc,h.v_lab))
        factor=D.from_float(float(np.pi))*v0*v0/(D.from_float(h.normalization_3d)*boost)
        for fraction in (1e-4,1e-7):
            a=h.max_lab_speed*(1-fraction)
            points=[D.from_float(a)]
            kink=abs(escape-boost)
            if points[0]<kink<escape+boost:points.append(kink)
            points.append(escape+boost)
            integral=D(0)
            for lo,hi in zip(points[:-1],points[1:]):
                for node,weight in zip(nodes,weights):
                    v=(lo+hi)/2+(hi-lo)*D.from_float(float(node))/2
                    upper=min(v+boost,escape)
                    kernel=(-((v-boost)/v0)**2).exp()-(-(upper/v0)**2).exp()
                    integral+=D.from_float(float(weight))*(hi-lo)/2*factor*kernel
            assert mean_inverse_speed_shm(h,a)==pytest.approx(float(integral),rel=1e-7,abs=1e-25)


def arguments(config,A,mass,cp=1e-5,cn=1e-5):
    isotope=next(i for i,_ in target(config,'natural') if i.A==A)
    model,end=signal_model(config,isotope,mass,1)
    args=model._rate_arguments()
    for name in ('speed_pdf','beta_max','breakpoints'):args.pop(name)
    args['coefficients']=WilsonCoefficients.from_proton_neutron({1:(cp,cn)})
    h=StandardHaloModel.from_km_s(**{k:config['signal']['halo'][k] for k in ('v0_km_s','vesc_km_s','v_lab_km_s')})
    return args,h,end


@pytest.mark.parametrize('A',[70,72,73,74,76])
@pytest.mark.parametrize('mass',[3.,10.,100.])
@pytest.mark.parametrize('cp,cn',[(1e-5,1e-5),(1e-5,.4e-5),(1e-5,-1e-5),
                               (0.,1e-5),(1e-5,0.),(0.,0.),(-1e-5,-.4e-5)])
def test_all_isotopes_masses_and_signed_couplings(config,A,mass,cp,cn):
    args,h,end=arguments(config,A,mass,cp,cn)
    e=end*np.array([[0.,.001,.25],[.999999,1.,1.01]])
    reference=differential_rate_shm_per_kg_day_keV(e,halo=h,**args)
    actual=differential_rate_shm_per_kg_day_keV(e,halo=h,use_fast_si=True,**args)
    np.testing.assert_allclose(actual,reference,rtol=3e-9,atol=np.max(abs(reference))*2e-13)
    # The inverse endpoint can land one ulp below vmax; do not change the
    # generic accessibility predicate to force a rounded energy to zero.
    assert actual[1,2]==0
    assert abs(actual[1,1])<=np.max(abs(actual))*1e-25
    assert isinstance(differential_rate_shm_per_kg_day_keV(0.,halo=h,use_fast_si=True,**args),float)
    assert differential_rate_shm_per_kg_day_keV(np.empty((0,2)),halo=h,use_fast_si=True,**args).shape==(0,2)


@pytest.mark.parametrize('op',[op for op in CANONICAL_OPERATORS if op!=1])
@pytest.mark.parametrize('tau',[0,1])
@pytest.mark.parametrize('coefficient',[np.nextafter(0.,1.),-1e-5])
def test_every_non_si_entry_forces_generic(config,monkeypatch,op,tau,coefficient):
    import wimp_legend.si_fast as module
    args,h,_=arguments(config,73,10.)
    pair=[0.,0.];pair[tau]=coefficient
    args['coefficients']=WilsonCoefficients.from_isospin({1:(1e-5,0.),op:pair})
    assert not is_o1_only(args['coefficients'])
    seen=[]
    def generic(e,**kw):seen.append(kw);return 123.
    def forbidden(*a,**kw):raise AssertionError('SI eta was selected for a mixed operator')
    monkeypatch.setattr(module,'differential_rate_per_kg_day_keV',generic)
    monkeypatch.setattr(module,'mean_inverse_speed_shm',forbidden)
    assert module.differential_rate_shm_per_kg_day_keV(.1,halo=h,use_fast_si=True,**args)==123.
    assert seen[0]['coefficients'] is args['coefficients']


@pytest.mark.parametrize('op',[3,8,11])
def test_mixed_operator_actual_fallback_is_identical(config,op):
    args,h,_=arguments(config,73,10.)
    args['coefficients']=WilsonCoefficients.from_proton_neutron({1:(1e-5,.4e-5),op:(1e-3,-2e-3)})
    e=np.array([0.,.1,1.])
    np.testing.assert_array_equal(differential_rate_shm_per_kg_day_keV(e,halo=h,**args),
        differential_rate_shm_per_kg_day_keV(e,halo=h,use_fast_si=True,**args))


def test_opt_in_validation_and_dataset_gates(config):
    args,h,_=arguments(config,76,10.)
    bad=replace(args['dataset'],metadata=replace(args['dataset'].metadata,convention_id='other'))
    with pytest.raises(ValueError):differential_rate_shm_per_kg_day_keV(0,halo=h,use_fast_si=True,**dict(args,dataset=bad))
    entries=dict(args['dataset'].responses);entries.pop(next(k for k in entries if k.channel!='M'))
    bad=replace(args['dataset'],responses=entries)
    with pytest.raises(KeyError):differential_rate_shm_per_kg_day_keV(0,halo=h,use_fast_si=True,**dict(args,dataset=bad))
    for e in (True,-1,np.nan,1j):
        with pytest.raises(ValueError):differential_rate_shm_per_kg_day_keV(e,halo=h,use_fast_si=True,**args)
    with pytest.raises(ValueError):differential_rate_shm_per_kg_day_keV(0,halo=h,use_fast_si=1,**args)
    model,_=signal_model(config,args['dataset'].isotope,10.,1)
    assert not model.use_fast_si
    with pytest.raises(ValueError):replace(model,use_fast_si=True,halo=h,speed_pdf=lambda b:0.)
    c=deepcopy(config);c['signal']['si_rate_engine']='unknown'
    with pytest.raises(ValueError):validate_config(c)


def test_disabled_dispatch_and_ordered_signed_synthetic_response(config,monkeypatch):
    import wimp_legend.si_fast as module
    from test_rates import synthetic
    args,h,_=arguments(config,76,10.,2e-5,-4e-5)
    # Deliberately distinguish 01 from 10 to detect accidental symmetrization.
    args['dataset']=synthetic({('M',0,0):(3.,),('M',0,1):(-.2,),
                               ('M',1,0):(.7,),('M',1,1):(2.,)})
    e=np.array([0.,.1,1.])
    reference=module.differential_rate_shm_per_kg_day_keV(e,halo=h,**args)
    actual=module.differential_rate_shm_per_kg_day_keV(e,halo=h,use_fast_si=True,**args)
    np.testing.assert_allclose(actual,reference,rtol=3e-12,atol=0)
    seen=[]
    def spy(e,**kw):seen.append(kw);return 42.
    monkeypatch.setattr(module,'differential_rate_per_kg_day_keV',spy)
    assert module.differential_rate_shm_per_kg_day_keV(0,halo=h,**args)==42.
    assert seen and is_o1_only(args['coefficients'])


def test_si_option_does_not_invalidate_non_si_cache(config,tmp_path,monkeypatch):
    import analyses.legend_poster.analysis as module
    c=resolve_profile(config,'smoke');isotope=next(i for i,_ in target(c,'natural') if i.A==76)
    cache=NumericalCache(tmp_path,source_identity()['physics_fingerprint']);bins=observed_edges(c)
    a,_=pure_signal_counts(c,isotope,6.,8,'benchmark',bins,cache);misses=cache.misses
    c['signal']['si_rate_engine']='factorized_shm'
    def forbidden(*a,**kw):raise AssertionError('unchanged non-SI rate was recalculated')
    monkeypatch.setattr(module,'differential_rate_per_kg_day_keV',forbidden)
    b,_=pure_signal_counts(c,isotope,6.,8,'benchmark',bins,cache)
    np.testing.assert_array_equal(a,b);assert cache.misses==misses


def test_isospin_violating_and_pure_neutron_folded_limits(config):
    from scipy.interpolate import PchipInterpolator
    from wimp_legend.backgrounds import TruthBackground
    from wimp_legend.energy import EnergyDomain
    from analyses.legend_poster.models import responses
    from analyses.legend_poster.analysis import fold_component
    _,response=responses(config)
    bins=np.array([.1,.2,.5,1.,2.,4.])
    for cp,cn in [(1e-5,-.7e-5),(0.,1e-5),(1e-5,0.),(0.,0.)]:
        args,h,end=arguments(config,73,10.,cp,cn)
        grid=end*np.linspace(0,1,32)**2
        counts=[]
        for enabled in (False,True):
            truth=differential_rate_shm_per_kg_day_keV(grid,halo=h,use_fast_si=enabled,**args)
            source=TruthBackground(EnergyDomain.NR,response,PchipInterpolator(grid,truth),
                (0.,end),spectrum_breakpoints_keV=tuple(grid))
            counts.append(fold_component(source,bins,config['numerical']))
        np.testing.assert_allclose(counts[1],counts[0],rtol=3e-9,atol=max(counts[0])*2e-13)
        # A synthetic observed-flat background isolates rate/fold/statistical
        # equivalence; it is not a new physical LEGEND background component.
        background=.01*np.diff(bins)
        a,_=count_limit(counts[0],background,bins,365.25,.1,config,'poisson_cls_counting')
        b,_=count_limit(counts[1],background,bins,365.25,.1,config,'poisson_cls_counting')
        assert b.limit.mu_upper==pytest.approx(a.limit.mu_upper,rel=3e-9)


def test_folded_cache_ignores_exposure_background_and_threshold(config,tmp_path,monkeypatch):
    from wimp_legend.si import SIO1Signal
    c=resolve_profile(config,'smoke');c['signal']['si_rate_engine']='factorized_shm'
    isotope=next(i for i,_ in target(c,'natural') if i.A==76);bins=observed_edges(c)
    cache=NumericalCache(tmp_path,source_identity()['physics_fingerprint'])
    a,_=pure_signal_counts(c,isotope,6.,1,'benchmark',bins,cache)
    misses=cache.misses
    def forbidden(*a,**kw):raise AssertionError('truth rate was needlessly recalculated')
    monkeypatch.setattr(SIO1Signal,'__call__',forbidden)
    c['detector']['presets']['benchmark']['threshold_keV_obs']=.2
    c['analysis']['exposures_kg_day']=[99.]
    c['background']['base']['lar_transfer_probability']*=2
    c['numerical']['background_order']+=2
    np.testing.assert_array_equal(observed_edges(c),bins)
    b,_=pure_signal_counts(c,isotope,6.,1,'benchmark',bins,cache)
    np.testing.assert_array_equal(a,b);assert cache.misses==misses
    c['detector']['presets']['benchmark']['noise_sigma_keV']*=2
    b,_=pure_signal_counts(c,isotope,6.,1,'benchmark',bins,cache)
    assert cache.misses>misses and not np.array_equal(a,b)
