"""Generic reference recoil assembly; explicit physics inputs and unit boundaries."""

from .recoil import (
    NREFTRateConfig, transition_probability, differential_cross_section_GeV_minus3,
    integrate_speed_flux, differential_rate_per_kg_day_keV,
)

__all__ = [
    'NREFTRateConfig', 'transition_probability', 'differential_cross_section_GeV_minus3',
    'integrate_speed_flux', 'differential_rate_per_kg_day_keV',
]
