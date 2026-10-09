"""Pure efficiency functions; application stage is owned by DetectorResponse."""
from dataclasses import dataclass
from typing import Callable, Protocol
import numpy as np
from .._analysis import array, elementwise, output, scalar


class Efficiency(Protocol):
    def __call__(self, energy_keV):
        """Return a probability in [0,1], scalar/array matching the input."""
        ...


@dataclass(frozen=True)
class ConstantEfficiency:
    epsilon: float

    def __post_init__(self):
        value = scalar(self.epsilon, 'epsilon')
        if value > 1:
            raise ValueError('epsilon must be <= 1')
        object.__setattr__(self, 'epsilon', value)

    def __call__(self, energy_keV):
        # Negative reconstructed energies are legitimate Gaussian outcomes.
        energy = array(energy_keV, 'energy_keV', nonnegative=False)
        return output(np.full_like(energy, self.epsilon))


@dataclass(frozen=True)
class CallableEfficiency:
    callback: Callable

    def __post_init__(self):
        if not callable(self.callback):
            raise TypeError('efficiency callback must be callable')

    def __call__(self, energy_keV):
        return elementwise(self.callback, energy_keV, 'efficiency', probability=True,
                           nonnegative_input=False)
