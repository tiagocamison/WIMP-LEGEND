"""Exact unit conversions; no halo, target, or interaction parameters."""

from typing import Final

C_KM_S: Final[float] = 299_792.458
KEV_TO_GEV: Final[float] = 1e-6

# Exact SI definitions: BIPM, SI Brochure 9th edition, Sections 2 and 4/Table 8.
# https://www.bipm.org/en/measurement-units/si-defining-constants
# https://www.bipm.org/en/publications/si-brochure
# h and e are exact; pi and the derived float values are rounded numerically.
from math import pi

PLANCK_J_S: Final[float] = 6.62607015e-34
ELEMENTARY_CHARGE_C: Final[float] = 1.602176634e-19
GEV_TO_J: Final[float] = 1e9 * ELEMENTARY_CHARGE_C
C_CM_S: Final[float] = 1e5 * C_KM_S
SECONDS_PER_DAY: Final[float] = 86400.0
HBAR_C_GEV_CM: Final[float] = PLANCK_J_S * C_CM_S / (2 * pi * GEV_TO_J)
GEV_MINUS2_TO_CM2: Final[float] = HBAR_C_GEV_CM**2
