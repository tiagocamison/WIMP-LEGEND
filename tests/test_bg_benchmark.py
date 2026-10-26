"""Focused benchmark checks against mathematical and evaluated source references."""
from dataclasses import FrozenInstanceError,replace
from pathlib import Path
import json,math,sys
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import ndtr
from wimp_legend._analysis import IntegrationPolicy
from wimp_legend.backgrounds import (AllowedBetaShape,TabulatedBetaShape,BetaBranch,
    CaptureShell,ElectronCaptureBranch,DecayChain,mean_activity,BenchmarkConfig,
    build_benchmark,TruthBackground,CompositeBackground,BinnedBackground)
from wimp_legend.backgrounds.benchmark import benchmark_inputs,solar_fluxes
from wimp_legend.backgrounds.neutrinos import FluxQuadrature,SolarSurvivalBenchmark,SolarScatteringSpectrum
from wimp_legend.detector import (ElectronRecoilResponse,DetectorResponse,NoSmearing,
    ConstantQuenching,ConstantEfficiency,GaussianResolution)
from wimp_legend.energy import EnergyDomain

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools/background_benchmark'))
from folding import quadrature_counts
from legacy import load_legacy


@pytest.fixture
def cfg():
    record=json.loads((ROOT/'examples/bg_benchmark_config.json').read_text())['scenarios']['enriched_UAr']
    return BenchmarkConfig(**record)


def responses(efficiency=1.):
    kwargs=dict(observed_efficiency=ConstantEfficiency(efficiency),policy=IntegrationPolicy(epsrel=1e-5,epsabs=1e-12))
    return ElectronRecoilResponse(NoSmearing(),**kwargs),DetectorResponse(ConstantQuenching(.2),NoSmearing(),**kwargs)


@pytest.mark.parametrize('minus',[True,False])
@pytest.mark.parametrize('Q,Z,A',[(18.564,2,3),(565,19,39),(100,0,10)])
def test_allowed_shape_full_support_normalization_and_owned_state(minus,Q,Z,A):
    shape=AllowedBetaShape(Q,Z,A,beta_minus=minus)
    assert quad(shape,0,Q,epsabs=2e-10)[0]==pytest.approx(1,rel=3e-10)
    e=np.array([-.1,0,.1,Q/2,Q,Q+1]);before=e.copy();out=shape(e)
    assert out.shape==e.shape and np.all(out>=0)
    assert out[0]==out[-1]==out[-2]==0
    assert out[1]>0 if minus and Z else out[1]==0
    assert isinstance(shape(.1),float);np.testing.assert_array_equal(e,before)
    with pytest.raises(FrozenInstanceError):shape.q_keV=5
    # No cached_property mutation through the public object's __dict__.
    assert shape._normalization>0


def test_allowed_zero_charge_independent_phase_space_oracle():
    Q=30.;me=510.99895
    reference=lambda e:math.sqrt(e*(e+2*me))*(e+me)*(Q-e)**2
    norm=quad(reference,0,Q,epsabs=.01)[0]
    shape=AllowedBetaShape(Q,0,10)
    for e in (.001,.1,10,29):assert shape(e)==pytest.approx(reference(e)/norm,rel=2e-12)


def test_tabulation_full_domain_normalization_no_caller_mutation():
    e=np.array([0.,1.,2.]);r=np.array([1.,2.,0.])
    shape=TabulatedBetaShape(e,r,'synthetic linear density')
    e[1]=.1;r[1]=100
    assert shape.q_keV==2
    assert shape.unnormalized_dndE==(1.,2.,0.)
    assert shape.source_integral==2.5
    assert quad(shape,0,2,points=[1])[0]==pytest.approx(1)
    np.testing.assert_array_equal(shape([-1,3]),[0,0])
    assert shape(1)==pytest.approx(2/2.5)
    for bad in ([False,1,2],[0,np.inf,2]):
        with pytest.raises(ValueError):TabulatedBetaShape(bad,[1,1,0],'test')


def test_ar39_first_forbidden_unique_frozen_evaluated_input():
    data,solar,ar,manifest=benchmark_inputs()
    assert isinstance(ar,TabulatedBetaShape)
    assert '1st-forbidden-unique' in ar.source_id
    assert '22b8691dd3ea34e6b69fcb93167ab872673ffd5b' in ar.source_id
    e=np.array(ar.energies_keV);r=ar(e)
    assert e[0]==0 and e[-1]==565 and r[-1]==0
    assert np.sum(np.diff(e)*(r[:-1]+r[1:])/2)==pytest.approx(1,rel=3e-15)
    # ENSDF EAV=218.8(21) keV and BetaShape report218(2) keV.
    moment=np.sum(np.diff(e)*(e[:-1]*r[:-1]+e[1:]*r[1:])/2)
    assert moment==pytest.approx(218.8,abs=2.1)
    # Low-energy density is finite; legacy no-Coulomb allowed spectrum vanishes.
    assert ar(0)>0.002
    assert 'default numerical screening' in ar.source_id
    h3=AllowedBetaShape(data['decays']['H3']['endpoint_keV'],2,3)
    assert quad(lambda t:t*h3(t),0,h3.q_keV)[0]==pytest.approx(5.68,abs=.02)


@pytest.mark.parametrize('name',('Ge68','Ga68','Zn65','Fe55'))
def test_evaluated_ec_shell_probability_denominators_and_total_cascades(name):
    data,*_=benchmark_inputs();rec=data['decays'][name]
    for row in rec['ec_branches']:
        probabilities=[row['PK'],row['PL'],row['PMplus']]
        assert math.fsum(probabilities)==pytest.approx(1,abs=2e-16)
        assert probabilities[2]>0
        assert rec['vacancy_energies_keV'][0]>rec['vacancy_energies_keV'][1]>rec['vacancy_energies_keV'][2]>0
    if name=='Fe55':assert rec['vacancy_energies_keV'][0]==6.539  # Mn total K vacancy, not5.90keV X-ray
    if name=='Ge68':assert rec['ec_branches'][0]['PK']==.8639


@pytest.mark.parametrize('escape',(0.,.3,1.))
def test_ec_gamma_splitting_preserves_parent_decay_branch(escape):
    branch=ElectronCaptureBranch(.4,(CaptureShell('K',10,.8),CaptureShell('L',1,.15),CaptureShell('M+',.1,.05)),1000,'synthetic')
    lines=branch.lines(5,gamma_escape_probability=escape)
    assert len(lines)==6
    assert sum(l.rate_per_kg_day for l in lines)==pytest.approx(2.)
    assert sum(l.rate_per_kg_day for l in lines if l.energy_keV<20)==pytest.approx(2*escape)
    assert {l.energy_domain for l in lines}=={EnergyDomain.ER}
    with pytest.raises(ValueError):ElectronCaptureBranch(1,(CaptureShell('K',10,.8),),0,'incomplete')
    with pytest.raises(ValueError):branch.lines(5,gamma_escape_probability=1.1)


@pytest.mark.parametrize('times',( (10.,1.),(1.,10.),(1.,1.),(1.,1.+1e-12) ))
def test_bateman_initial_and_integrated_activities(times):
    hp,hd=times;chain=DecayChain(hp,hd,.7);Ap,Ad=3.,2.
    assert chain.activities_at(Ap,Ad,0)==(Ap,Ad)
    mean=chain.mean_activities(Ap,Ad,start_day=2,duration_day=5)
    # Independent time integration of analytic two-member solution, stable at equality.
    lp,ld=math.log(2)/hp,math.log(2)/hd
    def daughter(t):
        if abs(lp-ld)<1e-10:return math.exp(-ld*t)*(Ad+.7*ld*Ap*t)
        return Ad*math.exp(-ld*t)+.7*ld*Ap*(math.exp(-lp*t)-math.exp(-ld*t))/(ld-lp)
    assert mean[0]==pytest.approx(quad(lambda t:Ap*math.exp(-lp*t),2,7)[0]/5,rel=1e-13)
    assert mean[1]==pytest.approx(quad(daughter,2,7)[0]/5,rel=2e-12)
    assert mean_activity(Ap,hp,start_day=2,duration_day=5)==pytest.approx(mean[0],rel=1e-13)


def test_ge68_ga68_growth_and_equilibrium():
    data,*_=benchmark_inputs();D=data['decays'];chain=DecayChain(D['Ge68']['half_life_day'],D['Ga68']['half_life_day'])
    Ap,Ad=chain.activities_at(10,0,1)
    assert Ad/Ap==pytest.approx(1.,rel=.001)
    Ap,Ad=chain.activities_at(10,0,10)
    assert Ad/Ap==pytest.approx(1/(1-chain.daughter_half_life_day/chain.parent_half_life_day),rel=1e-12)
    assert chain.activities_at(10,0,0)[1]==0
    with pytest.raises(ValueError):DecayChain(0,1)
    with pytest.raises(ValueError):chain.mean_activities(1,0,start_day=0,duration_day=0)


def test_fixed_source_flux_normalization_and_survival():
    _,solar,_,_=benchmark_inputs()
    fluxes=solar_fluxes(solar,('pp','B8'))
    assert math.fsum(fluxes[0].weights_cm2_s)==pytest.approx(5.98e10,rel=3e-15)
    assert math.fsum(fluxes[1].weights_cm2_s)==pytest.approx(5.46e6,rel=3e-15)
    P=SolarSurvivalBenchmark(.307,.022,7.53e-5,65,1.1663787e-5,6.02214076e23)
    e=np.array([0,.01,.1,1,10,100]);v=P(e)
    assert np.all((v>=0)&(v<=1)) and np.all(np.diff(v)<0)
    vacuum=(1-.022)**2*(1-2*.307*(1-.307))+.022**2
    assert v[0]==pytest.approx(vacuum,rel=2e-15)
    assert v[-1]<v[1]-.1
    vacuum_model=replace(P,rho_Ye_g_cm3=0)
    np.testing.assert_allclose(vacuum_model(e),np.full(e.shape,vacuum),rtol=2e-15)


def test_solar_es_absolute_integral_and_cevens_zero_momentum():
    GF=1.1663787e-5;me=.00051099895;sw=.23863;E=.001;flux=1e10
    f=FluxQuadrature((1.,),(flux,),'synthetic monoenergetic neutrinos')
    spec=SolarScatteringSpectrum('ES',(f,),76,32,70.8,1,GF,sw,me,lambda e:np.ones_like(e),(),'synthetic free-electron check')
    tmax=2*E*E/(me+2*E)
    gl=.5+sw;gr=sw
    primitive=gl*gl*tmax+gr*gr*(tmax-tmax*tmax/E+tmax**3/(3*E*E))-gl*gr*me*tmax*tmax/(2*E*E)
    # Independent SI hbar*c conversion, not the production constant.
    hc=6.62607015e-34*29979245800/(2*math.pi*1.602176634e-10)
    expected=2*GF*GF*me/math.pi*primitive*hc*hc*flux*86400*32
    assert quad(spec,0,tmax*1e6,epsabs=1e-20)[0]==pytest.approx(expected,rel=2e-12)
    assert spec(tmax*1e6+1)==0
    nr=replace(spec,process='CEvNS')
    Q=44-(1-4*sw)*32
    assert nr(0)==pytest.approx(GF*GF*70.8*Q*Q/(4*math.pi)*hc*hc*1e-6*flux*86400,rel=2e-15)
    assert nr(2*E*E/70.8*1e6)==pytest.approx(0,abs=1e-20)


def test_binding_step_and_factory_energy_domains(cfg):
    er,nr=responses();m=build_benchmark(cfg,er_response=er,nr_response=nr)
    assert isinstance(m.composite,CompositeBackground)
    assert m.components['solar_CEvNS'].energy_domain is EnergyDomain.NR
    for key in ('H3','Ge68','Ga68','Zn65','Fe55','Ar39','solar_ES'):
        assert m.components[key].energy_domain is EnergyDomain.ER
        assert m.components[key].response is er
    es=m.components['solar_ES'].spectrum
    ratio=es(np.array([.02919,.02921]));assert ratio[1]>ratio[0]*2
    for ch in ('Ge68','Ga68','Zn65','Fe55'):
        assert all(l.energy_domain is EnergyDomain.ER for l in m.components[ch].lines)
    with pytest.raises(ValueError):build_benchmark(cfg,er_response=er,nr_response=nr,component_responses={'Ge68':nr})
    with pytest.raises(ValueError):build_benchmark(cfg,er_response=er,nr_response=nr,component_responses={'typo':er})
    with pytest.raises(TypeError):m.components['new']=None
    with pytest.raises(FrozenInstanceError):m.config.cooldown_day=1


def test_activity_transfer_veto_and_response_efficiency_separation(cfg):
    er,nr=responses(.4);m=build_benchmark(cfg,er_response=er,nr_response=nr)
    er2,_=responses(.8);changed=build_benchmark(cfg,er_response=er2,nr_response=nr)
    bins=[0,20]
    # EC branch intensities precede response. Efficiency enters once only.
    key='Ge68';one=m.components[key].expected_counts(bins,2)
    two=changed.components[key].expected_counts(bins,2)
    assert two==pytest.approx(2*one)
    assert m.components['flat_observed'].expected_counts(bins,2)==pytest.approx(changed.components['flat_observed'].expected_counts(bins,2))
    half=build_benchmark(replace(cfg,lar_transfer_probability=cfg.lar_transfer_probability/2),er_response=er,nr_response=nr)
    assert half.components['Ar39'].spectrum(1)==pytest.approx(m.components['Ar39'].spectrum(1)/2)
    assert json.loads(half.provenance)['mean_ar39_decays_per_kg_Ge_day_before_transfer']==json.loads(m.provenance)['mean_ar39_decays_per_kg_Ge_day_before_transfer']
    no_veto=build_benchmark(replace(cfg,lar_veto_survival_probability=0),er_response=er,nr_response=nr)
    assert no_veto.components['Ar39'].spectrum(1)==0
    np.testing.assert_array_equal(m.components['Ge68'].expected_counts(bins,0),[0])
    np.testing.assert_allclose(m.components['Ge68'].expected_counts(bins,10),5*one)


def test_nuclear_gamma_containment_and_parent_normalization(cfg):
    er,nr=responses();full=build_benchmark(cfg,er_response=er,nr_response=nr)
    escaped=build_benchmark(replace(cfg,zn65_gamma_escape_probability=1),er_response=er,nr_response=nr)
    a=full.components['Zn65'];b=escaped.components['Zn65']
    assert sum(l.rate_per_kg_day for l in a.lines)==pytest.approx(sum(l.rate_per_kg_day for l in b.lines))
    low=lambda c:sum(l.rate_per_kg_day for l in c.lines if l.energy_keV<15)
    assert low(b)/low(a)==pytest.approx((.4835+.5023)/.4835)
    # Beta+ annihilation is contained by default: continuum support >=1022keV.
    assert a.spectrum_energy_range_keV[0]==pytest.approx(2*510.99895)
    ratio=sum(l.rate_per_kg_day for l in full.components['Ga68'].lines)/sum(l.rate_per_kg_day for l in full.components['Ge68'].lines)
    assert .10<ratio<.12  # Ge/Ga activities linked; only ~11.1% of Ga decays are EC.


def test_plot_quadrature_independent_integral_and_lines():
    response=ElectronRecoilResponse(GaussianResolution(.1),truth_efficiency=ConstantEfficiency(.5),observed_efficiency=ConstantEfficiency(.7))
    source=TruthBackground(EnergyDomain.ER,response,lambda e:np.ones_like(e),(0,1))
    bins=np.array([-.1,.2,.6,1.2]);actual=quadrature_counts(source,bins,3,step_keV=.08,order=8)
    phi=lambda z:math.exp(-z*z/2)/math.sqrt(2*math.pi)
    G=lambda z:z*ndtr(z)+phi(z)
    cdf_int=lambda edge:.1*(G(edge/.1)-G((edge-1)/.1))
    expected=3*.5*.7*np.diff([cdf_int(e) for e in bins])
    np.testing.assert_allclose(actual,expected,rtol=2e-13,atol=1e-15)
    clipped=quadrature_counts(source,bins,3,step_keV=.04,order=6,tail_sigma=12,sigma_bound=lambda c:.1)
    np.testing.assert_allclose(clipped,expected,rtol=2e-13)
    # Observed binned spectra bypass response entirely.
    observed=BinnedBackground(bins,[1,2,3],'per_unit_exposure')
    np.testing.assert_array_equal(quadrature_counts(observed,bins,3),[3,6,9])


def test_legacy_notebook_reference_is_not_modified():
    source=ROOT/'WIMP_O1_spectrum.ipynb';before=source.read_bytes()
    legacy,digest=load_legacy(source)
    cfg=json.loads((ROOT/'examples/bg_benchmark_config.json').read_text())
    assert digest==cfg['legacy_notebook_sha256']
    values=legacy(np.array([.1,1.,5.9,10.37]),target='enrGe',argon='UAr')
    assert set(values)>={'68Ge','65Zn','55Fe','39Ar','CEvNS','nu_e_ES','Total'}
    assert np.isfinite(values['Total']).all()
    np.testing.assert_allclose(values['Total'],sum(v for k,v in values.items() if k!='Total'),rtol=2e-15)
    assert source.read_bytes()==before


def test_pp_coherent_total_flux_and_narrow_recoil_support(cfg,monkeypatch):
    import wimp_legend.backgrounds.neutrinos as module
    er,nr=responses();model=build_benchmark(cfg,er_response=er,nr_response=nr)
    component=model.components['solar_CEvNS'];spec=component.spectrum
    pp=spec.fluxes[0];spec=replace(spec,fluxes=(pp,))
    endpoint=2*(max(pp.energies_MeV)*1e-3)**2/spec.m_target_GeV*1e6
    assert endpoint<.01
    assert any(0<p<endpoint for p in component.spectrum_breakpoints_keV)
    monkeypatch.setattr(module,'helm_squared',lambda q,A:np.ones_like(q))
    points=sorted({2*(e*1e-3)**2/spec.m_target_GeV*1e6 for e in pp.energies_MeV})
    got=quad(spec,0,endpoint,points=points[:-1],epsabs=1e-12,epsrel=1e-9,limit=1000)[0]
    Q=spec.A-spec.Z-(1-4*spec.sin2_theta_W)*spec.Z
    hc=6.62607015e-34*29979245800/(2*math.pi*1.602176634e-10)
    moment=math.fsum(w*(e*1e-3)**2 for w,e in zip(pp.weights_cm2_s,pp.energies_MeV))
    expected=spec.G_F_GeV_minus2**2*Q**2/(4*math.pi)*moment*hc*hc*86400*spec.targets_per_kg
    assert got==pytest.approx(expected,rel=2e-12)


def test_composite_counts_exposure_and_replacement(cfg):
    er,nr=responses(.8);model=build_benchmark(cfg,er_response=er,nr_response=nr)
    bins=np.array([.1,1,5,15]);one=model.expected_counts(bins,1.)
    assert np.all(one>0)
    np.testing.assert_allclose(model.expected_counts(bins,2.),2*one,rtol=1e-12)
    np.testing.assert_allclose(one,sum((c.expected_counts(bins,1) for c in model.components.values()),np.zeros(3)),rtol=1e-12)
    replacement=BinnedBackground(bins,[1,2,3],'per_unit_exposure',provenance='future measured/MC observed template')
    changed=dict(model.components);changed['Ar39']=replacement
    result=CompositeBackground(tuple(changed.values())).expected_counts(bins,1)
    expected=one-model.components['Ar39'].expected_counts(bins,1)+np.array([1,2,3])
    np.testing.assert_allclose(result,expected,rtol=1e-12)


@pytest.mark.parametrize('change',[
 {'lar_transfer_probability':-1},{'lar_veto_survival_probability':2},
 {'run_duration_day':0},{'cooldown_day':-1},{'ga68_initialization':'implicit'},
 {'ga68_initial_activity_per_kg_day':1},{'solar_flux_components':['CNO']},
 {'activities_per_kg_day':[('H3',1),('Ge68',2),('Ge68',2),('Fe55',1)]},
])
def test_invalid_config_is_rejected(cfg,change):
    with pytest.raises(ValueError):replace(cfg,**change)


def test_config_input_ownership(cfg):
    pairs=[list(p) for p in cfg.activities_per_kg_day]
    fluxes=['pp','B8'];owned=replace(cfg,activities_per_kg_day=pairs,solar_flux_components=fluxes)
    pairs[0][1]=999;fluxes.clear()
    assert dict(owned.activities_per_kg_day)['H3']==dict(cfg.activities_per_kg_day)['H3']
    assert owned.solar_flux_components==('pp','B8')


@pytest.mark.parametrize('escape',(0.,.3,1.))
def test_all_parent_branches_are_accounted_once(cfg,escape):
    er,nr=responses();model=build_benchmark(replace(cfg,zn65_gamma_escape_probability=escape,
        ga68_gamma_escape_probability=escape,annihilation_escape_probability=escape),er_response=er,nr_response=nr)
    activities=json.loads(model.provenance)['mean_intrinsic_activities_per_kg_day']
    for name in ('Ge68','Ga68','Zn65','Fe55'):
        component=model.components[name]
        integrated=sum(line.rate_per_kg_day for line in component.lines)
        if component.spectrum is not None:
            integrated+=quad(component.spectrum,*component.spectrum_energy_range_keV,
                points=[v for v in component.spectrum_breakpoints_keV if component.spectrum_energy_range_keV[0]<v<component.spectrum_energy_range_keV[1]],
                epsabs=1e-12,epsrel=1e-10,limit=100)[0]
        assert integrated==pytest.approx(activities[name],rel=2e-10)


@pytest.mark.parametrize('filename',('manifest.json','nuclear.json','ar39.csv'))
def test_background_input_integrity_rejects_tampering(tmp_path,monkeypatch,filename):
    import shutil
    from wimp_legend.backgrounds import benchmark as module
    root=Path(module.__file__).parent/'data/benchmark_v1'
    shutil.copytree(root,tmp_path/'data/benchmark_v1')
    changed=tmp_path/'data/benchmark_v1'/filename
    with changed.open('a') as f:f.write(' ')
    monkeypatch.setattr(module,'files',lambda package:tmp_path)
    with pytest.raises(ValueError,match='checksum mismatch'):module.benchmark_inputs()
