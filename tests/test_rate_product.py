"""Independent second audit. Synthetic nuclei only; no repository test imports."""
from dataclasses import replace
from decimal import Decimal as D, localcontext
from math import pi, sqrt, exp, erf
from pathlib import Path
import ast
import numpy as np
import pytest
from scipy.integrate import quad, IntegrationWarning
from wimp_legend.targets import (SourceReference, Isotope, TargetComposition, ResponseMetadata,
    MomentumVariable, ResponseKey, PolynomialResponse, NuclearResponseDataset)
from wimp_legend.interactions import WilsonCoefficients, particle_response
from wimp_legend.kinematics import momentum_transfer_GeV, minimum_speed_c
from wimp_legend.halo import StandardHaloModel, build_velocity_integral
from wimp_legend.rates import (NREFTRateConfig, transition_probability as prob,
    differential_cross_section_GeV_minus3 as cross, integrate_speed_flux as flux,
    differential_rate_per_kg_day_keV as rate)
from wimp_legend import constants

ID='anand2014_prc89_065501_half_isospin_v1'
CH=('M','Sigma_double_prime','Sigma_prime','Phi_double_prime','Phi_double_prime_M',
    'Phi_tilde_prime','Delta','Delta_Sigma_prime')
KEYS=[(k,t,u) for k in CH for t in (0,1) for u in (0,1)]
SOURCE=SourceReference('Independent synthetic audit', 'audit-only', 'synthetic')

def fixture(entries=None, J=1.5, momentum=None):
    iso=Isotope('X', 31, 12, 23., J, dict.fromkeys(('identity','mass_GeV','spin'),SOURCE))
    meta=ResponseMetadata('synthetic basis',('0','1'),'half isospin','raw W',
        'dimensionless',momentum or MomentumVariable('q2_GeV2'),SOURCE,ID)
    responses={ResponseKey(*k):PolynomialResponse((0.,),0,None,SOURCE) for k in KEYS}
    for k,v in (entries or {}).items():
        responses[ResponseKey(*k)]=PolynomialResponse(tuple(v),0,None,SOURCE)
    return NuclearResponseDataset('independent synthetic fixture',iso,meta,responses)

def pdf(L):
    return lambda b: 3*(np.asarray(b)/L)**2/L

def si():
    with localcontext() as ctx:
        ctx.prec=75
        p=D('3.141592653589793238462643383279502884197169399375105820974944592307816406286')
        hc=D('6.62607015e-34')*D('299792458')*100/(2*p*D('1.602176634e-19')*D('1e9'))
        return float(hc),float(hc*hc),float(hc*hc*D('1e-6')*D('29979245800')*86400)

def kwargs(ds=None,raw=None,mref=1.7,L=.006):
    return dict(dataset=fixture({('M',0,0):(1.3,)}) if ds is None else ds,
        coefficients=WilsonCoefficients({1:(3e-6,0)} if raw is None else raw),
        config=NREFTRateConfig(mref),j_chi=1.,m_chi_GeV=17.,
        rho_chi_GeV_cm3=.43,targets_per_kg=2.8e25,speed_pdf=pdf(L),beta_max=L)

def test_RATE_AUDIT_001_representable_final_rate_erased():
    # Intentionally red falsifier: an intermediate area conversion underflows
    # before a large, ordinary nuclei/kg factor restores the final result.
    k=kwargs(raw={1:(1e-150,0)})
    with localcontext() as ctx:
        ctx.prec=90
        p=D('3.141592653589793238462643383279502884197169399375105820974944592307816406286')
        hc=D('6.62607015e-34')*D('29979245800')/(2*p*D('1.602176634e-10'))
        # At E=0, J=3/2, P=pi*w*g² and integral=mT*w*g²*3/(4L).
        expected=float(hc**2*D('1e-6')*D('29979245800')*86400*D('.43')/17*
            D('2.8e25')*23*D('1.3')*D('1e-150')**2*3/(4*D('.006')))
    assert expected>0 and np.isfinite(expected)
    assert rate(0.,**k)==pytest.approx(expected,rel=1e-12,abs=0)


@pytest.mark.parametrize('g', [1e-150, 3e-6, .02])
def test_rate_product_scaling_and_sign(g):
    k=kwargs(raw={1:(g,0)})
    baseline=rate(0.,**k)
    for factor in [.25, 4.]:
        scaled=dict(k, targets_per_kg=k['targets_per_kg']*factor)
        assert rate(0.,**scaled)==pytest.approx(baseline*factor,rel=2e-15,abs=0)
    k['dataset']=fixture({('M',0,0):(-1.3,)})
    assert rate(0.,**k)==-baseline

@pytest.mark.parametrize('factors,divisor', [
    ((1e-300,1e-100,1e300),1.), ((1e300,1e300,1e-300),1.),
    ((1e300,1e300),1e300), ((1e-300,1e-300),1e-300),
    ((-2.,3.,1e-320),1.), ((0.,1e308,1e308),1.),
])
def test_scaled_product_decimal(factors,divisor):
    from wimp_legend.rates.recoil import _scaled_product
    with localcontext() as ctx:
        ctx.prec=100
        expected=D(1)
        for f in factors:expected*=D(f)
        expected=float(expected/D(divisor))
    assert _scaled_product(*factors,divisor=divisor)==pytest.approx(expected,rel=5e-16,abs=0)

@pytest.mark.parametrize('factors', [(1e308,1e308),(1e-300,1e-300)])
def test_scaled_product_true_range_failure(factors):
    from wimp_legend.rates.recoil import _scaled_product
    with pytest.raises(ValueError,match='range'):_scaled_product(*factors)
