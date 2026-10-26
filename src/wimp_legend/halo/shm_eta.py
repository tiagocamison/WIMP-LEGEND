"""Vectorized SHM mean inverse speed, in beta=v/c units.

Integrating the accepted shm.lab_speed_pdf analytically gives the usual erf
branches. Zero boost and the upper tail use expm1 to avoid cancellation.
Reference: Freese, Lisanti & Savage, Rev.Mod.Phys.85 (2013) 1561,
arXiv:1209.3339v3 Appendix B Eq.(B4), hard truncation (their beta=0).
Narrow exponential integrals are evaluated by 16-node Gaussian quadrature;
there is no approximate small-boost replacement or physical cutoff.
"""
import numpy as np
from scipy.special import erfc, gammainc
from numpy.polynomial.legendre import leggauss
from .shm import StandardHaloModel
from .integrals import _thresholds


def mean_inverse_speed_shm(halo, v_min_c):
    """eta=integral f1(beta)/beta d beta, preserving scalar/array shape.

    Threshold policy matches velocity_moment: negative clamps to zero,
    +infinity is outside support, NaN is rejected. No normalization in km/s.
    In the final 0.01*v0 of support, a stable 16-node Gaussian integral of
    the exponential difference replaces a cancellation-prone erf subtraction.
    """
    if not isinstance(halo, StandardHaloModel):
        raise TypeError('halo must be a StandardHaloModel')
    a=_thresholds(v_min_c);out=np.zeros_like(a)
    active=a<halo.max_lab_speed
    if not np.any(active):return float(out) if out.ndim==0 else out
    x=a[active]/halo.v0;z=halo.z;y=halo.v_lab/halo.v0
    norm=float(gammainc(1.5,z*z))
    if y==0:
        value=(np.exp(-x*x)*(-np.expm1(-(z-x)*(z+x)))/norm)*(2/(np.sqrt(np.pi)*halo.v0))
    else:
        value=np.empty_like(x)
        full=x<abs(z-y)
        if y>=z:
            value[full]=1/halo.v_lab
        else:
            u=x[full]
            if y<1e-3 or z<.1:
                nodes,weights=leggauss(16)
                t=u[:,None]+y*nodes
                integrand=np.exp(-t*t)*(-np.expm1(-(z-t)*(z+t)))
                value[full]=((integrand@weights)/norm)/(np.sqrt(np.pi)*halo.v0)
            else:
                value[full]=(erfc(u-y)-erfc(u+y)-4*y/np.sqrt(np.pi)*np.exp(-z*z))/(2*norm*halo.v_lab)
        partial=~full;u=x[partial]
        # Form delta in dimensional speeds: x+y-z can erase the last ulp of
        # a threshold next to the endpoint, changing a physical tiny tail to 0.
        delta=(halo.max_lab_speed-a[active][partial])/halo.v0
        tail=(delta<.01)|(z<.1)
        branch=np.empty_like(u)
        branch[~tail]=(erfc(u[~tail]-y)-erfc(z)-2*delta[~tail]/np.sqrt(np.pi)*np.exp(-z*z))/(2*norm*halo.v_lab)
        if np.any(tail):
            nodes,weights=leggauss(16)
            t=delta[tail,None]*(nodes+1)/2
            integrand=np.exp(-(z-t)**2)*(-np.expm1(-t*(2*z-t)))
            branch[tail]=(integrand@weights)*delta[tail]/(2*np.sqrt(np.pi)*norm*halo.v_lab)
        value[partial]=branch
    if np.any(~np.isfinite(value)) or np.any(value<0):
        raise ValueError('SHM eta is not finite/nonnegative at this numerical scale')
    out[active]=value
    return float(out) if out.ndim==0 else out
