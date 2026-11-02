"""Equivalent optional Gaussian batch backend, including reference fallbacks."""
from copy import deepcopy
from pathlib import Path
import sys
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from analyses.legend_poster import load_config,validate_config
from analyses.legend_poster.analysis import fold_component
from analyses.legend_poster.folding import VectorizedGaussianBins,LocalGaussianBins,prepare_folding
from analyses.legend_poster.models import responses
from wimp_legend.detector import (DetectorResponse,ElectronRecoilResponse,GaussianResolution,
    ConstantEfficiency,ConstantQuenching,LindhardQuenching,NoSmearing)
from wimp_legend.backgrounds import TruthBackground,MonoenergeticLine,CompositeBackground
from wimp_legend.energy import EnergyDomain

@pytest.mark.parametrize('domain',[EnergyDomain.ER,EnergyDomain.NR])
@pytest.mark.parametrize('sigma',[.003,.05,lambda e:.01+np.sqrt(e)*.02])
@pytest.mark.parametrize('epsilon',[None,ConstantEfficiency(.73)])
def test_gaussian_batch_matches_reference(domain,sigma,epsilon):
    resolution=GaussianResolution(sigma)
    response=(ElectronRecoilResponse(resolution,observed_efficiency=epsilon) if domain is EnergyDomain.ER
        else DetectorResponse(LindhardQuenching(32,.157),resolution,observed_efficiency=epsilon))
    bins=np.unique(np.r_[np.linspace(-.2,.2,73),np.linspace(.21,3,61)])
    e=np.linspace(0,5,270).reshape(3,3,30)
    original=e.copy();edges=bins.copy()
    actual=VectorizedGaussianBins(response,12).bin_probabilities(e,bins)
    expected=LocalGaussianBins(response,12).bin_probabilities(e,bins)
    np.testing.assert_array_equal(actual,expected)
    np.testing.assert_array_equal(e,original);np.testing.assert_array_equal(bins,edges)
    assert actual.shape==e.shape+(len(bins)-1,)
    assert np.all(actual>=0) and np.all(actual.sum(axis=-1)<=1+1e-14)
    assert VectorizedGaussianBins(response,12).bin_probabilities(.1,bins).shape==(len(bins)-1,)
    assert VectorizedGaussianBins(response,12).bin_probabilities(np.array([]),bins).shape==(0,len(bins)-1)

@pytest.mark.parametrize('domain',[EnergyDomain.ER,EnergyDomain.NR])
def test_weighted_efficiency_falls_back_exactly(domain):
    efficiency=lambda e:.4+.1*np.asarray(e)
    resolution=GaussianResolution(.03)
    response=(ElectronRecoilResponse(resolution,observed_efficiency=efficiency) if domain is EnergyDomain.ER
        else DetectorResponse(ConstantQuenching(.2),resolution,observed_efficiency=efficiency))
    fast=VectorizedGaussianBins(response,12)
    assert not fast._vectorized
    e=np.array([.1,.5]);bins=np.array([0,.1,.3])
    np.testing.assert_array_equal(fast.bin_probabilities(e,bins),LocalGaussianBins(response,12).bin_probabilities(e,bins))

@pytest.mark.parametrize('domain',[EnergyDomain.ER,EnergyDomain.NR])
@pytest.mark.parametrize('noise',[.003,.05])
def test_continuum_and_lines_use_same_quadrature(domain,noise):
    c=load_config();c['detector']['presets']['benchmark']['noise_sigma_keV']=noise
    er,nr=responses(c);response=er if domain is EnergyDomain.ER else nr
    # A nonconstant truth efficiency belongs outside the Gaussian kernel.
    from dataclasses import replace
    response=replace(response,truth_efficiency=lambda e:.5+.1*np.asarray(e)/(1+np.asarray(e)))
    component=TruthBackground(domain,response,lambda e:np.exp(-np.asarray(e)),(0.,3.),
        lines=(MonoenergeticLine(.13,.7,domain),MonoenergeticLine(1.2,.2,domain)),
        spectrum_breakpoints_keV=(.2,1.))
    bins=np.unique(np.r_[np.linspace(.03,1,99),1.2,1.3,2.])
    n=c['numerical'];reference=fold_component(component,bins,n)
    fast=dict(n,folding_backend='vectorized')
    actual=fold_component(component,bins,fast)
    np.testing.assert_allclose(actual,reference,rtol=2e-14,atol=max(reference)*1e-15)
    proxy=prepare_folding(component,12,'vectorized')
    np.testing.assert_array_equal(proxy.response.expected_line_counts(.13,.7,bins,2),response.expected_line_counts(.13,.7,bins,2))
    composite=CompositeBackground((component,component))
    np.testing.assert_allclose(fold_component(composite,bins,fast),2*actual,rtol=2e-14,atol=0)


def test_backend_default_validation_and_unsupported_resolution():
    c=load_config();assert c['numerical']['folding_backend']=='reference'
    del c['numerical']['folding_backend'];validate_config(c)
    c['numerical']['folding_backend']='invented'
    with pytest.raises(ValueError,match='folding_backend'):validate_config(c)
    response=ElectronRecoilResponse(NoSmearing())
    cpt=TruthBackground(EnergyDomain.ER,response,lines=(MonoenergeticLine(.13,1,EnergyDomain.ER),))
    assert prepare_folding(cpt,12,'vectorized') is cpt
    with pytest.raises(ValueError):prepare_folding(cpt,12,'invalid')


def test_zero_width_matches_pruned_reference():
    response=ElectronRecoilResponse(GaussianResolution(lambda e:0. if e<.5 else .01),observed_efficiency=ConstantEfficiency(.6))
    e=np.array([.123,.234,.7]);bins=np.linspace(0,1,101)
    np.testing.assert_array_equal(VectorizedGaussianBins(response,12).bin_probabilities(e,bins),LocalGaussianBins(response,12).bin_probabilities(e,bins))

@pytest.mark.parametrize('mass',[3.,10.,100.])
@pytest.mark.parametrize('A',[70,72,73,74,76])
def test_physical_isotope_folds_and_limits(mass,A,tmp_path):
    from analyses.legend_poster import resolve_profile
    from analyses.legend_poster.models import target
    from analyses.legend_poster.analysis import pure_signal_counts,observed_edges,count_limit
    from analyses.legend_poster.storage import NumericalCache
    c=resolve_profile(load_config(),'standard');c['signal']['si_rate_engine']='factorized_shm'
    isotope=next(i for i,_ in target(c,'natural') if i.A==A)
    bins=observed_edges(c);cache=NumericalCache(tmp_path,'test-source')
    reference,refcheck=pure_signal_counts(c,isotope,mass,1,'benchmark',bins,cache)
    c['numerical']['folding_backend']='vectorized'
    actual,check=pure_signal_counts(c,isotope,mass,1,'benchmark',bins,cache)
    assert check['passed'] and refcheck['passed']
    np.testing.assert_allclose(actual,reference,rtol=3e-8,atol=max(reference)*1e-13)
    background=.01*np.diff(bins)
    results=[]
    for values in (reference,actual):
        _,record=count_limit(values,background,bins,3652.5,.1,c,'poisson_cls_counting');results.append(record)
    assert results[0]['roi_upper_keV']==results[1]['roi_upper_keV']
    assert results[1]['mu_upper']==pytest.approx(results[0]['mu_upper'],rel=3e-8,abs=0)


def test_custom_kernel_and_observed_data_keep_their_own_routes():
    from wimp_legend.detector import KernelResponse
    from wimp_legend.backgrounds import FlatBackground
    from wimp_legend._analysis import IntegrationPolicy
    response=KernelResponse(lambda e,bins:np.array([.2,.4]),includes_efficiency=True,
        policy=IntegrationPolicy(epsabs=1e-10))
    cpt=TruthBackground(EnergyDomain.NR,response,lambda e:2.,(0.,1.))
    n=load_config()['numerical'];n['folding_backend']='vectorized';bins=np.array([0,.1,1.])
    np.testing.assert_array_equal(fold_component(cpt,bins,n),cpt.expected_counts(bins,1.))
    observed=FlatBackground(.7)
    np.testing.assert_array_equal(fold_component(observed,bins,n),observed.expected_counts(bins,1.))
