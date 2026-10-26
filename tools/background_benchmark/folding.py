"""Numerical plotting helper; production response/count interfaces stay intact.

Fixed Gauss-Legendre continuum integration, point lines through the existing
exact bin-probability API. Optional Gaussian tail bounds are EXPLICIT numerical
approximations: a supplied global sigma upper bound gives omitted mass at most
2*Phi(-tail_sigma) per event. The plotting caller supplies the bound from its
monotonic noise+Fano model. Other templates use their existing expected_counts.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from wimp_legend._analysis import edges,scalar,output
from wimp_legend.energy import EnergyDomain
from wimp_legend.backgrounds import CompositeBackground,TruthBackground


def quadrature_counts(background,bins,exposure,*,step_keV=.05,order=6,tail_sigma=12.,sigma_bound=None):
    bins=edges(bins);exposure=scalar(exposure,'exposure')
    step=scalar(step_keV,'step_keV',positive=True);tail=scalar(tail_sigma,'tail_sigma',positive=True)
    if isinstance(order,bool) or not isinstance(order,int) or not 2<=order<=32:raise ValueError('quadrature order must be 2..32')
    if exposure==0:return np.zeros(len(bins)-1)
    if isinstance(background,CompositeBackground):
        return output(sum((quadrature_counts(c,bins,exposure,step_keV=step,order=order,
                        tail_sigma=tail,sigma_bound=sigma_bound) for c in background.components),np.zeros(len(bins)-1)))
    if not isinstance(background,TruthBackground):return background.expected_counts(bins,exposure)
    response=background.response;result=np.zeros(len(bins)-1)
    for line in background.lines:
        result+=response.expected_line_counts(line.energy_keV,line.rate_per_kg_day,bins,exposure)
    if background.spectrum is None:return output(result)
    lo,hi=background.spectrum_energy_range_keV
    if not callable(getattr(response,'bin_probabilities',None)):
        # A custom response/MC port remains independently replaceable.
        return background.expected_counts(bins,exposure)
    mean=(lambda e:float(e)) if background.energy_domain is EnergyDomain.ER else (
        lambda e:float(response.quenching.mean_energy_keV(e)))
    if sigma_bound is not None:
        bound=scalar(sigma_bound(background),'global Gaussian sigma bound')
        min_mean,max_mean=mean(lo),mean(hi)
        left,right=bins[0]-tail*bound,bins[-1]+tail*bound
        if min_mean>right or max_mean<left:return output(result)
        if min_mean<left:lo=brentq(lambda e:mean(e)-left,lo,hi)
        if max_mean>right:hi=brentq(lambda e:mean(e)-right,lo,hi)
    if lo>=hi:return output(result)
    # Bin crossings make the no-smearing limit and thresholds well resolved.
    crossings=[]
    for edge in bins:
        if mean(lo)<edge<mean(hi):crossings.append(brentq(lambda e:mean(e)-edge,lo,hi))
    count=int(np.ceil((hi-lo)/step))
    if count>100000:raise ValueError('plot quadrature exceeds 100000 truth intervals; increase step_keV')
    points=np.unique(np.r_[np.linspace(lo,hi,count+1),crossings,
        [v for v in background.spectrum_breakpoints_keV if lo<v<hi]])
    x,w=leggauss(order)
    nodes=((points[:-1,None]+points[1:,None])/2+np.diff(points)[:,None]*x/2).ravel()
    weights=(np.diff(points)[:,None]*w/2).ravel()
    for start in range(0,len(nodes),128):
        e=nodes[start:start+128];weight=weights[start:start+128]
        rate=np.asarray(background.spectrum(e))
        if rate.shape!=e.shape or np.any(~np.isfinite(rate)) or np.any(rate<0):raise ValueError('invalid vector truth spectrum')
        probabilities=np.asarray(response.bin_probabilities(e,bins))
        truth_efficiency=getattr(response,'truth_efficiency',None)
        efficiency=np.ones_like(e) if truth_efficiency is None else np.asarray(truth_efficiency(e))
        if efficiency.shape!=e.shape or np.any(~np.isfinite(efficiency)) or np.any((efficiency<0)|(efficiency>1)):
            raise ValueError('invalid truth efficiency')
        result+=exposure*((weight*rate*efficiency)@probabilities)
    return output(result)
