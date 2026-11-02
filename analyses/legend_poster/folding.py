"""Bounded Gaussian-bin pruning while delegating probabilities to core APIs.

This is a numerical optimization, not a new detector model. Outside bins are
at least tail_sigma standard deviations from every truth mean in this batch.
Point lines retain the original exact, unpruned core probability calculation.
"""
from dataclasses import replace
import numpy as np
from wimp_legend.backgrounds import TruthBackground,CompositeBackground


class LocalGaussianBins:
    def __init__(self,response,tail_sigma):
        self.response=response;self.tail_sigma=tail_sigma
        self.input_energy_domain=response.input_energy_domain
        self.resolution=response.resolution
        self.truth_efficiency=response.truth_efficiency
        self.observed_efficiency=response.observed_efficiency
        if hasattr(response,'quenching'):self.quenching=response.quenching

    def bin_probabilities(self,energy,bins):
        from wimp_legend.energy import EnergyDomain
        e=np.asarray(energy)
        mean=e if self.input_energy_domain is EnergyDomain.ER else self.quenching.mean_energy_keV(e)
        sigma=np.asarray(self.resolution.sigma_at(mean))
        lower=float(np.min(mean-self.tail_sigma*sigma));upper=float(np.max(mean+self.tail_sigma*sigma))
        first=max(0,int(np.searchsorted(bins,lower,side='right')-1))
        last=min(len(bins)-1,int(np.searchsorted(bins,upper,side='left')))
        result=np.zeros(e.shape+(len(bins)-1,))
        if last>first:
            result[...,first:last]=self.response.bin_probabilities(e,np.asarray(bins)[first:last+1])
        return result

    def expected_line_counts(self,*args,**kwargs):
        return self.response.expected_line_counts(*args,**kwargs)

    def expected_counts(self,*args,**kwargs):
        return self.response.expected_counts(*args,**kwargs)


def prune_gaussian_bins(component,tail_sigma):
    if isinstance(component,CompositeBackground):
        return CompositeBackground(tuple(prune_gaussian_bins(c,tail_sigma) for c in component.components))
    if isinstance(component,TruthBackground):
        return replace(component,response=LocalGaussianBins(component.response,tail_sigma))
    return component


class VectorizedGaussianBins(LocalGaussianBins):
    """Same local Gaussian kernel, evaluated in bounded NumPy batches.

    Only the concrete built-in ER/NR responses with Gaussian resolution and
    absent/constant observed efficiency qualify. All other response models
    retain their own probability implementation. Truth efficiency is applied
    by the unchanged quadrature, and point lines use the original response.
    """
    def __init__(self,response,tail_sigma):
        from wimp_legend.detector import DetectorResponse, ElectronRecoilResponse
        from wimp_legend.detector.resolution import GaussianResolution
        from wimp_legend.detector.efficiency import ConstantEfficiency
        super().__init__(response,tail_sigma)
        self._vectorized=(type(response) in (DetectorResponse,ElectronRecoilResponse)
            and type(response.resolution) is GaussianResolution
            and (response.observed_efficiency is None
                 or type(response.observed_efficiency) is ConstantEfficiency))

    def bin_probabilities(self,energy,bins):
        from scipy.special import ndtr
        from wimp_legend._analysis import array,edges,output
        from wimp_legend.energy import EnergyDomain
        if not self._vectorized:
            return super().bin_probabilities(energy,bins)
        e=array(energy,'truth energy_keV');bins=edges(bins)
        mean=e if self.input_energy_domain is EnergyDomain.ER else array(
            self.quenching.mean_energy_keV(e),'E_ee_keV')
        if mean.shape!=e.shape:raise ValueError('quenching must preserve energy shape')
        from .models import NoiseFanoSigma
        sigma_function=self.resolution.sigma_keV
        sigmas=(array(sigma_function(mean),'sigma_keV') if type(sigma_function) is NoiseFanoSigma
                else np.asarray(self.resolution.sigma_at(mean)))
        # Preserve the reference wrapper's deterministic endpoint rules in
        # degenerate batches; do not reimplement zero-width migration here.
        if np.any(sigmas==0):return super().bin_probabilities(e,bins)
        result=np.zeros((mean.size,len(bins)-1))
        if not mean.size:return output(result.reshape(mean.shape+(len(bins)-1,)))
        means=mean.ravel();sigmas=sigmas.ravel()
        epsilon=1. if self.observed_efficiency is None else self.observed_efficiency.epsilon
        # Bound temporary storage even for direct array calls outside quadrature.
        # Slice the same batch-union bin window as the reference wrapper.
        lower=float(np.min(means-self.tail_sigma*sigmas));upper=float(np.max(means+self.tail_sigma*sigmas))
        first=max(0,int(np.searchsorted(bins,lower,side='right')-1))
        last=min(len(bins)-1,int(np.searchsorted(bins,upper,side='left')))
        for start in range(0,len(means),128):
            end=min(start+128,len(means));m=means[start:end];s=sigmas[start:end]
            positive=s>0
            if last>first and np.any(positive):
                with np.errstate(over='raise',invalid='raise',divide='raise'):
                    z=(bins[None,first:last+1]-m[positive,None])/s[positive,None]
                a,b=z[:,:-1],z[:,1:]
                # Survival probabilities on the positive side match the core
                # exactly; do not subtract two CDFs rounded close to one.
                side=a>=0
                hi=ndtr(np.where(side,-a,b));lo=ndtr(np.where(side,-b,a))
                result[start+np.flatnonzero(positive),first:last]=(hi-lo)*epsilon
        return output(result.reshape(mean.shape+(len(bins)-1,)))


def prepare_folding(component,tail_sigma,backend='reference'):
    """Optional numerical backend; no changes to component physics or routing."""
    if backend not in ('reference','vectorized'):raise ValueError('unknown folding backend')
    if isinstance(component,CompositeBackground):
        return CompositeBackground(tuple(prepare_folding(c,tail_sigma,backend) for c in component.components))
    if isinstance(component,TruthBackground):
        wrapper=VectorizedGaussianBins if backend=='vectorized' else LocalGaussianBins
        # Replace only the built-in analytic responses. Other ports may lack
        # Gaussian attributes; leave their existing integration API intact.
        from wimp_legend.detector import DetectorResponse,ElectronRecoilResponse
        from wimp_legend.detector.resolution import GaussianResolution
        if (type(component.response) in (DetectorResponse,ElectronRecoilResponse)
                and type(component.response.resolution) is GaussianResolution):
            return replace(component,response=wrapper(component.response,tail_sigma))
        if backend=='reference':return prune_gaussian_bins(component,tail_sigma)
    return component
