"""Electron-recoil response: mean E_ee=E_er, with no NR quenching input."""
from dataclasses import dataclass, field

from .._analysis import IntegrationPolicy
from ..energy import EnergyDomain
from .efficiency import Efficiency
from .quenching import ConstantQuenching
from .resolution import ResolutionModel
from .response import DetectorResponse


@dataclass(frozen=True)
class ElectronRecoilResponse:
    """ER keV -> observed keV. Resolution sigma is evaluated at E_er_keV.

    Reuse the validated detector fold with an identity energy mapping, never an
    NR quenching curve. Efficiencies are ER-truth/observed probabilities. No
    empirical electron nonlinearity or energy calibration is assumed/supplied.
    """
    resolution: ResolutionModel
    truth_efficiency: Efficiency | None = None
    observed_efficiency: Efficiency | None = None
    policy: IntegrationPolicy = IntegrationPolicy()
    observed_breakpoints: tuple[float, ...] = ()
    _response: DetectorResponse = field(init=False, repr=False)

    def __post_init__(self):
        response = DetectorResponse(ConstantQuenching(1.), self.resolution,
            self.truth_efficiency, self.observed_efficiency, self.policy,
            self.observed_breakpoints)
        object.__setattr__(self, 'observed_breakpoints', response.observed_breakpoints)
        object.__setattr__(self, '_response', response)

    @property
    def input_energy_domain(self):
        return EnergyDomain.ER

    def bin_probabilities(self, E_er_keV, observed_bin_edges_keV):
        return self._response.bin_probabilities(E_er_keV, observed_bin_edges_keV)

    def expected_counts(self, truth_spectrum, observed_bin_edges_keV, exposure_kg_day,
                        truth_energy_range_keV, *, truth_breakpoints=()):
        """truth_spectrum(E_er_keV) has units events/(kg day keV_er)."""
        return self._response.expected_counts(truth_spectrum, observed_bin_edges_keV,
            exposure_kg_day, truth_energy_range_keV, truth_breakpoints=truth_breakpoints)

    def expected_line_counts(self, E_er_keV, rate_per_kg_day,
                             observed_bin_edges_keV, exposure_kg_day):
        return self._response.expected_line_counts(E_er_keV, rate_per_kg_day,
                                                   observed_bin_edges_keV, exposure_kg_day)
