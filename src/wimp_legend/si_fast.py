"""Opt-in O1-only SHM factorization; the generic NREFT engine is unchanged.

From the accepted recoil assembly (AFH PRC89 Eq.40,50):
dR/dE = Nt*rho/mchi * c*day*(hbar*c)^2*1e-6
         * 2*mT/(2*J+1) * sum_(tau,u) c1^tau*c1^u*W_M^(tau,u) * eta.
Wilson coefficients are dimensionful GeV^-2 in the frozen half-isospin basis.
Ordered entries/signs are preserved; no Helm model or coherent A^2 repair.
"""
import numpy as np
from .constants import C_CM_S,GEV_MINUS2_TO_CM2,KEV_TO_GEV,SECONDS_PER_DAY
from .interactions import WilsonCoefficients
from .halo.shm import StandardHaloModel
from .halo.shm_eta import mean_inverse_speed_shm
from .kinematics import minimum_speed_c,momentum_transfer_GeV
from .rates import differential_rate_per_kg_day_keV
from .rates.recoil import _array,_output,_scaled_product


def is_o1_only(coefficients):
    """Exact stored-coefficient predicate, including all-zero coefficients.

    No tolerance, squared-amplitude test, spin/target selection rule or coupling
    cancellation can authorize this optimization for a nonzero non-O1 entry.
    """
    if not isinstance(coefficients,WilsonCoefficients):
        raise TypeError('coefficients must be WilsonCoefficients')
    return all(c==0 for op,pair in coefficients.values.items() if op!=1 for c in pair)


def differential_rate_shm_per_kg_day_keV(E_nr_keV,*,halo,use_fast_si=False,**arguments):
    """Dispatch to factorized SI iff opted in AND every non-O1 entry is zero.

    Other arguments are exactly the generic single-isotope rate arguments
    except speed_pdf/beta_max/breakpoints, which are owned by the explicit SHM.
    Generic fallback handles any non-SI coefficient, however small. Both modes
    enforce the same dataset/convention/config validation, including all 32 W
    entries. Inaccessible recoil energies do not evaluate nuclear responses.
    Signed nuclear contractions remain signed; no abs/clipping/renormalization.
    """
    if type(use_fast_si)!=bool:raise ValueError('use_fast_si must be bool')
    if not isinstance(halo,StandardHaloModel):raise TypeError('halo must be a StandardHaloModel')
    generic=dict(arguments,speed_pdf=halo.lab_speed_pdf,beta_max=halo.max_lab_speed,
                 breakpoints=halo.integration_breakpoints)
    # Do not silently accept a second, contradictory PDF/support declaration.
    if any(k in arguments for k in ('speed_pdf','beta_max','breakpoints')):
        raise ValueError('SHM dispatcher owns speed_pdf, beta_max and breakpoints')
    selected=use_fast_si and is_o1_only(arguments['coefficients'])
    if not selected:return differential_rate_per_kg_day_keV(E_nr_keV,**generic)
    differential_rate_per_kg_day_keV(np.empty(0),**generic)
    energy=_array(E_nr_keV,'E_nr_keV');dataset=arguments['dataset']
    mass=arguments['m_chi_GeV'];minimum=minimum_speed_c(energy,mass,dataset.isotope.mass_GeV)
    active=minimum<halo.max_lab_speed;result=np.zeros_like(energy)
    if np.any(active):
        q=momentum_transfer_GeV(energy[active],dataset.isotope.mass_GeV);q2=q*q
        if np.any((q>0)&(q2==0)):raise ValueError('q2_GeV2 underflowed; unsupported numerical scale')
        contraction=np.zeros_like(q2);c=arguments['coefficients']
        with np.errstate(over='raise',invalid='raise',divide='raise'):
            for t in (0,1):
                for u in (0,1):
                    nuclear=dataset.evaluate(channel='M',tau=t,tau_prime=u,q2_GeV2=q2)
                    contraction+=c.get(1,t)*c.get(1,u)*nuclear
        eta=mean_inverse_speed_shm(halo,minimum[active])
        values=np.empty_like(eta)
        for i in np.ndindex(eta.shape):
            values[i]=_scaled_product(contraction[i],eta[i],2.,dataset.isotope.mass_GeV,
                1/(2*dataset.isotope.spin+1),
                arguments['targets_per_kg'],arguments['rho_chi_GeV_cm3'],C_CM_S,SECONDS_PER_DAY,
                GEV_MINUS2_TO_CM2,KEV_TO_GEV,divisor=mass)
        result[active]=values
    return _output(result)
