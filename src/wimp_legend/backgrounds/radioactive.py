"""Small radioactive source primitives, independent of detector and statistics.

Activities are decays/(kg day) of the explicitly named source mass. Branch
fractions refer to parent decays; shell probabilities are conditional on EC.
No geometric transfer, veto or detector efficiency is a nuclear branching ratio.
"""
from dataclasses import dataclass
from math import fsum, isclose, log
from typing import Callable
import numpy as np
from scipy.linalg import expm

from .._analysis import array, output, scalar, text
from ..energy import EnergyDomain
from .routing import MonoenergeticLine


def probability(value, name):
    value = scalar(value, name)
    if value > 1:
        raise ValueError(f'{name} must lie in [0,1]')
    return value


@dataclass(frozen=True)
class BetaBranch:
    shape: Callable
    endpoint_keV: float
    fraction_per_decay: float
    deposited_offset_keV: float = 0.
    provenance: str = 'user-supplied benchmark'

    def __post_init__(self):
        if not callable(self.shape):
            raise TypeError('beta shape must be callable')
        for name in ('endpoint_keV','deposited_offset_keV'):
            object.__setattr__(self,name,scalar(getattr(self,name),name,positive=name=='endpoint_keV'))
        object.__setattr__(self,'fraction_per_decay',probability(self.fraction_per_decay,'fraction_per_decay'))
        object.__setattr__(self,'provenance',text(self.provenance,'provenance'))

    def __call__(self, E_er_keV):
        energy = array(E_er_keV,'E_er_keV',nonnegative=False)
        return output(self.fraction_per_decay*np.asarray(self.shape(energy-self.deposited_offset_keV)))


@dataclass(frozen=True)
class CaptureShell:
    label: str
    vacancy_energy_keV: float
    probability_given_ec: float

    def __post_init__(self):
        object.__setattr__(self,'label',text(self.label,'label'))
        object.__setattr__(self,'vacancy_energy_keV',scalar(self.vacancy_energy_keV,'vacancy_energy_keV'))
        object.__setattr__(self,'probability_given_ec',probability(self.probability_given_ec,'probability_given_ec'))


@dataclass(frozen=True)
class ElectronCaptureBranch:
    fraction_per_decay: float
    shells: tuple[CaptureShell, ...]
    nuclear_excitation_keV: float
    provenance: str

    def __post_init__(self):
        shells=tuple(self.shells)
        if not shells or any(not isinstance(s,CaptureShell) for s in shells):
            raise TypeError('capture shells must be nonempty CaptureShell records')
        if len({s.label for s in shells}) != len(shells):
            raise ValueError('duplicate capture shell')
        if not isclose(fsum(s.probability_given_ec for s in shells),1,rel_tol=0,abs_tol=1e-12):
            raise ValueError('conditional shell probabilities must sum to one; no silent normalization')
        object.__setattr__(self,'shells',shells)
        object.__setattr__(self,'fraction_per_decay',probability(self.fraction_per_decay,'fraction_per_decay'))
        object.__setattr__(self,'nuclear_excitation_keV',scalar(self.nuclear_excitation_keV,'nuclear_excitation_keV'))
        object.__setattr__(self,'provenance',text(self.provenance,'provenance'))

    def lines(self, activity_per_kg_day, *, gamma_escape_probability):
        """All-or-nothing nuclear gamma containment; atomic cascade contained.

        Escape is a geometry BENCHMARK probability, independent of veto/efficiency. Partial gamma
        deposits and atomic escape require a future measured/MC template.
        M+ shell labels explicitly aggregate M and higher captures at M1.
        """
        activity=scalar(activity_per_kg_day,'activity_per_kg_day')
        escape=probability(gamma_escape_probability,'gamma_escape_probability')
        result=[]
        outcomes=((0.,1.),) if self.nuclear_excitation_keV == 0 else (
            (0.,escape),(self.nuclear_excitation_keV,1-escape))
        for shell in self.shells:
            for energy,weight in outcomes:
                result.append(MonoenergeticLine(shell.vacancy_energy_keV+energy,
                    activity*self.fraction_per_decay*shell.probability_given_ec*weight,EnergyDomain.ER))
        return tuple(result)


@dataclass(frozen=True)
class DecayChain:
    """Closed two-member Bateman chain; exact mean activities over a run.

    Initial activities are specified at t=0, before cooling. No replenishment,
    plate-out, removal or activity-to-acceptance conversion is inferred. An
    augmented matrix exponential handles equal/near-equal half-lives safely.
    """
    parent_half_life_day: float
    daughter_half_life_day: float
    daughter_production_branch: float = 1.

    def __post_init__(self):
        for name in ('parent_half_life_day','daughter_half_life_day'):
            object.__setattr__(self,name,scalar(getattr(self,name),name,positive=True))
        object.__setattr__(self,'daughter_production_branch',probability(
            self.daughter_production_branch,'daughter_production_branch'))

    def _matrix(self):
        lp=log(2)/self.parent_half_life_day;ld=log(2)/self.daughter_half_life_day
        return np.array([[-lp,0,0,0],[self.daughter_production_branch*ld,-ld,0,0],
                         [1,0,0,0],[0,1,0,0]],dtype=float)

    def activities_at(self, parent_activity_per_kg_day, daughter_activity_per_kg_day, time_day):
        initial=np.array([scalar(parent_activity_per_kg_day,'parent activity'),
                          scalar(daughter_activity_per_kg_day,'daughter activity'),0,0])
        time=scalar(time_day,'time_day')
        return tuple(output((expm(self._matrix()*time)@initial)[:2]))

    def mean_activities(self, parent_activity_per_kg_day, daughter_activity_per_kg_day, *,
                        start_day, duration_day):
        start=scalar(start_day,'start_day');duration=scalar(duration_day,'duration_day',positive=True)
        parent,daughter=self.activities_at(parent_activity_per_kg_day,daughter_activity_per_kg_day,start)
        integrated=expm(self._matrix()*duration)@np.array([parent,daughter,0,0])
        return tuple(output(integrated[2:]/duration))


def mean_activity(activity_per_kg_day, half_life_day, *, start_day, duration_day):
    activity=scalar(activity_per_kg_day,'activity_per_kg_day')
    half_life=scalar(half_life_day,'half_life_day',positive=True)
    start=scalar(start_day,'start_day');duration=scalar(duration_day,'duration_day',positive=True)
    decay=log(2)/half_life
    return float(activity*np.exp(-decay*start)*(-np.expm1(-decay*duration))/(decay*duration))
