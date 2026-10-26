"""Solar scattering benchmarks with explicit flux quadratures and weak inputs.

CEvNS: leading nonrelativistic Freedman cross section (PRD 9, 1389, 1974).
ES: tree-level neutrino-electron scattering, Vogel & Engel PRD 39, 3378 (1989).
No NREFT nuclear response, halo, composition registry or detector is imported.
"""
from dataclasses import dataclass
from typing import Callable
import numpy as np
from scipy.special import spherical_jn
from .._analysis import array, output, scalar, text
from ..constants import GEV_MINUS2_TO_CM2, SECONDS_PER_DAY, HBAR_C_GEV_CM


@dataclass(frozen=True)
class FluxQuadrature:
    """Point fluxes or quadrature weights, each in cm^-2 s^-1 (not a density)."""
    energies_MeV: tuple[float,...]
    weights_cm2_s: tuple[float,...]
    provenance: str

    def __post_init__(self):
        e=array(self.energies_MeV,'energies_MeV');w=array(self.weights_cm2_s,'weights_cm2_s')
        if e.ndim!=1 or not e.size or w.shape!=e.shape or np.any(e<=0) or not np.any(w>0):
            raise ValueError('positive neutrino energies and matching nonzero flux weights required')
        object.__setattr__(self,'energies_MeV',tuple(e));object.__setattr__(self,'weights_cm2_s',tuple(w))
        object.__setattr__(self,'provenance',text(self.provenance,'provenance'))


@dataclass(frozen=True)
class SolarSurvivalBenchmark:
    """Adiabatic, decohered MSW at one effective production density.

    sin² theta12/13 and dm² are explicit analysis inputs, as is rho*Ye in
    g/cm³. This neglects production-region averaging, Earth effects and solar
    density uncertainties. PDG neutrino-mixing review, adiabatic MSW formula.
    """
    sin2_theta12: float
    sin2_theta13: float
    delta_m21_eV2: float
    rho_Ye_g_cm3: float
    G_F_GeV_minus2: float
    avogadro_per_mol: float

    def __post_init__(self):
        for name in self.__dataclass_fields__:
            object.__setattr__(self,name,scalar(getattr(self,name),name,positive=name in ('delta_m21_eV2','G_F_GeV_minus2','avogadro_per_mol')))
        if self.sin2_theta12>1 or self.sin2_theta13>1:
            raise ValueError('mixing sin² must be in [0,1]')

    def __call__(self, E_nu_MeV):
        E=array(E_nu_MeV,'E_nu_MeV')*1e6
        s12=self.sin2_theta12;s13=self.sin2_theta13
        ne=self.rho_Ye_g_cm3*self.avogadro_per_mol
        potential=np.sqrt(2)*self.G_F_GeV_minus2*1e-18*ne*(HBAR_C_GEV_CM*1e9)**3
        cosine=1-2*s12
        matter=cosine-2*E*potential*(1-s13)/self.delta_m21_eV2
        denominator=np.sqrt(matter*matter+4*s12*(1-s12))
        cm=np.divide(matter,denominator,out=np.zeros_like(matter),where=denominator>0)
        return output((1-s13)**2*(.5+.5*cosine*cm)+s13**2)


def helm_squared(q_GeV, A):
    """Lewin & Smith Astropart.Phys.6 (1996) 87, Helm benchmark prescription.

    c=1.23*A^(1/3)-0.60 fm, a=0.52 fm, s=0.9 fm. Nuclear weak distribution
    approximated by this charge-style form factor; not an AFH/NREFT dataset.
    """
    q=array(q_GeV,'q_GeV');A=scalar(A,'A',positive=True)
    c=1.23*A**(1/3)-.60;a=.52;s=.9
    r2=c*c+7*np.pi*np.pi*a*a/3-5*s*s
    if r2<=0:raise ValueError('Helm prescription requires positive effective radius')
    qfm=q/(HBAR_C_GEV_CM*1e13);x=qfm*np.sqrt(r2)
    j=np.ones_like(x);nz=x!=0;j[nz]=3*spherical_jn(1,x[nz])/x[nz]
    return output(j*j*np.exp(-(qfm*s)**2))


@dataclass(frozen=True)
class SolarScatteringSpectrum:
    """Truth spectrum, events/(kg day keV), before detector response.

    process='CEvNS' means NR; 'ES' means ER. A may be an explicitly documented
    effective-mass-number surrogate, never an implicit physical isotope dataset.
    electron_shells contains (binding_keV, occupancy); empty means free Z
    electrons. The optional stepping approximation uses source shell energies.
    """
    process: str
    fluxes: tuple[FluxQuadrature,...]
    A: float
    Z: int
    m_target_GeV: float
    targets_per_kg: float
    G_F_GeV_minus2: float
    sin2_theta_W: float
    m_e_GeV: float
    survival_probability: Callable
    electron_shells: tuple[tuple[float,int],...]
    provenance: str

    def __post_init__(self):
        if self.process not in ('CEvNS','ES'):raise ValueError('process must be CEvNS or ES')
        fluxes=tuple(self.fluxes)
        if not fluxes or any(not isinstance(x,FluxQuadrature) for x in fluxes):raise TypeError('FluxQuadrature inputs required')
        object.__setattr__(self,'fluxes',fluxes)
        for name in ('A','m_target_GeV','targets_per_kg','G_F_GeV_minus2','m_e_GeV'):
            object.__setattr__(self,name,scalar(getattr(self,name),name,positive=True))
        if isinstance(self.Z,bool) or not isinstance(self.Z,int) or not 0<self.Z<=self.A:raise ValueError('invalid Z')
        sw=scalar(self.sin2_theta_W,'sin2_theta_W')
        if sw>1:raise ValueError('weak mixing sin² must be <=1')
        object.__setattr__(self,'sin2_theta_W',sw)
        if not callable(self.survival_probability):raise TypeError('explicit survival model required')
        shells=[]
        for binding,occupancy in self.electron_shells:
            binding=scalar(binding,'electron binding')
            if isinstance(occupancy,bool) or not isinstance(occupancy,int) or occupancy<=0:raise ValueError('positive integer shell occupancy required')
            shells.append((binding,occupancy))
        if shells and sum(n for _,n in shells)!=self.Z:raise ValueError('electron occupancies must sum to Z')
        object.__setattr__(self,'electron_shells',tuple(shells))
        object.__setattr__(self,'provenance',text(self.provenance,'provenance'))

    def __call__(self, recoil_keV):
        energies=array(recoil_keV,'recoil_keV');shape=energies.shape
        T=energies.reshape(-1,1)*1e-6
        total=np.zeros(len(T))
        for flux in self.fluxes:
            E=np.asarray(flux.energies_MeV)[None,:]*1e-3
            if self.process=='CEvNS':
                weak_charge=self.A-self.Z-(1-4*self.sin2_theta_W)*self.Z
                bracket=np.maximum(0.,1-self.m_target_GeV*T/(2*E*E))
                form=np.asarray(helm_squared(np.sqrt(2*self.m_target_GeV*T[:,0]),self.A))[:,None]
                ds=(self.G_F_GeV_minus2**2*self.m_target_GeV/(4*np.pi)*weak_charge**2*bracket*form)
            else:
                me=self.m_e_GeV;sw=self.sin2_theta_W;gr=sw
                Pee=np.asarray(self.survival_probability(np.asarray(flux.energies_MeV)))
                if Pee.shape!=(E.shape[1],) or np.any(~np.isfinite(Pee)) or np.any((Pee<0)|(Pee>1)):
                    raise ValueError('survival model must return matching probabilities')
                def kernel(gl):return gl*gl+gr*gr*(1-T/E)**2-gl*gr*me*T/(E*E)
                terms=Pee*kernel(.5+sw)+(1-Pee)*kernel(-.5+sw)
                active=(T<=2*E*E/(me+2*E))
                electrons=np.full(len(T),float(self.Z))
                if self.electron_shells:
                    electrons=np.zeros(len(T))
                    for binding,n in self.electron_shells:electrons+=n*(energies.reshape(-1)>=binding)
                ds=2*self.G_F_GeV_minus2**2*me/np.pi*np.maximum(0,terms)*active*electrons[:,None]
            total+=(ds@np.asarray(flux.weights_cm2_s))*GEV_MINUS2_TO_CM2*1e-6
        return output((total*self.targets_per_kg*SECONDS_PER_DAY).reshape(shape))
