"""Explicit energy domains for response routing; all energies are in keV."""
from enum import Enum


class EnergyDomain(Enum):
    NR = 'keV_nr'
    ER = 'keV_er'
    OBSERVED = 'keV_obs'
