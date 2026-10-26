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
