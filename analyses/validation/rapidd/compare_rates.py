"""Absolute elastic O1/O11 isotope rates from the pinned RAPIDD C calculator."""

import argparse
import ctypes as ct
import json
from itertools import product
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np

from compare_kernels import GE_ISOTOPES, shared_fixture
from compare_responses import PIN
from wimp_legend.constants import (C_CM_S, C_KM_S, GEV_MINUS2_TO_CM2,
                                   GEV_TO_J, KEV_TO_GEV, SECONDS_PER_DAY)
from wimp_legend.halo.shm import StandardHaloModel
from wimp_legend.halo.shm_eta import mean_inverse_speed_shm
from wimp_legend.interactions import WilsonCoefficients
from wimp_legend.kinematics import maximum_recoil_energy_keV, minimum_speed_c
from wimp_legend.rates import NREFTRateConfig, differential_rate_per_kg_day_keV
from wimp_legend.targets import Isotope, SourceReference

OPERATORS = (1, 11)
ISOTOPES = (73, 76)
MASSES_GEV = (5., 40.)
RECOIL_FRACTIONS = (*np.geomspace(.001, .95, 32), .99, 1., 1.01)
HALO_KM_S = (220., 544., 232.)  # v0, vesc, fixed lab speed
RHO_GEV_CM3 = .3
COUPLINGS_GEV_MINUS2 = (1e-6, 1e-6)
M_N_FORM_FACTOR_GEV = .938
VEV_GEV = 246.22  # phys_consts.h at the pinned revision
KG_IN_GEV = (1000*C_KM_S)**2/GEV_TO_J
# This stage uses only the inverse-speed term in difRateGen.c.
RAPIDD_ETA_CONVERSION = 1.698603e14
WIMP_ETA_CONVERSION = (
    KG_IN_GEV * GEV_MINUS2_TO_CM2 * KEV_TO_GEV * C_CM_S * SECONDS_PER_DAY * C_KM_S
)


class _RateParameters(ct.Structure):
    """Mirror the pinned difcros_params layout from difRateGen.h."""
    _fields_ = [
        ('target', ct.c_char_p), ('Er', ct.c_double), ('mchi', ct.c_double),
        ('jchi', ct.c_double), ('v', ct.c_double), ('Nucleon', ct.c_char_p),
        ('exposure', ct.c_double), ('delta', ct.c_double),
    ]


def _bind(library, name, result, *arguments):
    function = getattr(library, name)
    function.restype = result
    function.argtypes = arguments
    return function


def _build_library(checkout, build):
    cmake = shutil.which('cmake')
    if cmake is None:
        raise RuntimeError('Live absolute rates require CMake, a C compiler, and GSL development files.')
    for args in (
        ('-S', str(checkout/'lib'), '-B', str(build)),
        ('--build', str(build), '--target', 'RAPIDD', '--parallel', '2'),
    ):
        completed = subprocess.run([cmake, *args], capture_output=True, text=True, timeout=120)
        if completed.returncode:
            raise RuntimeError(f'RAPIDD build failed:\n{completed.stdout}\n{completed.stderr}')
    libraries = [path for name in ('libRAPIDD.so', 'libRAPIDD.dylib')
                 if (path := build/name).is_file()]
    if len(libraries) != 1:
        raise RuntimeError('Expected one RAPIDD shared library from the Linux/macOS build.')
    return ct.CDLL(str(libraries[0]))


def _calculate(library, isotopes, masses, operators, fractions, halo, cp, cn, rho):
    # Bind the C entry points directly. All state belongs to this fresh shared library.
    ff = _bind(library, 'FF_pn_Fitz', ct.c_double, ct.c_int, *([ct.c_char_p]*4), ct.c_double)
    mass_of = _bind(library, 'approx_mass_nucleus', ct.c_double, ct.c_double, ct.c_double)
    set_halo = _bind(library, 'fv_v', ct.c_int, ct.c_char_p, *([ct.c_double]*5),
                     ct.c_int, *([ct.c_double]*5))
    eta = _bind(library, 'halo', ct.c_double, ct.c_double, ct.c_int)
    reset = _bind(library, 'set_coeffs', None)
    set_coefficient = _bind(library, 'set_any_Ncoeff', None, ct.c_double, ct.c_int, ct.c_char_p)
    rate = _bind(library, 'total_difrate_isotope_dEr', ct.c_double,
                 ct.c_int, ct.c_int, ct.c_double, ct.POINTER(_RateParameters))
    v0, vesc, lab = halo.v0*C_KM_S, halo.vesc*C_KM_S, halo.v_lab*C_KM_S
    # Fill the same in-memory spline grid used by the C rate function.
    # The SHM path ignores beta, k, and the five LMC arguments.
    if set_halo(b'SHM', vesc, v0, 0., lab, 0., 2, 0., 0., 0., 0., 0.) != 0:
        raise RuntimeError('RAPIDD could not initialize the SHM table')

    source = SourceReference('RAPIDD '+PIN, 'Pointwise pn_Fitz form factors', 'implementation')
    config = NREFTRateConfig(M_N_FORM_FACTOR_GEV)
    cases = []
    for A, mchi, op in product(isotopes, masses, operators):
        kinematic_mass = mass_of(A, 32)
        # effFormFact.c uses its own rounded nucleon masses for q², unlike the rate routine.
        form_factor_mass = 32*.938+(A-32)*.940
        isotope = Isotope('Ge', A, 32, kinematic_mass, 4.5 if A == 73 else 0.,
                          {name: source for name in ('identity', 'mass_GeV', 'spin')})
        endpoint = float(maximum_recoil_energy_keV(mchi, kinematic_mass, halo.max_lab_speed))
        energies = fractions * endpoint
        vmin = minimum_speed_c(energies, mchi, kinematic_mass)
        coeff = WilsonCoefficients.from_proton_neutron({op: (cp, cn)})
        reset()
        set_coefficient(cp*VEV_GEV**2, op, b'p')
        set_coefficient(cn*VEV_GEV**2, op, b'n')
        params = _RateParameters(b'Ge', 0., mchi, .5, 0., b'pn_Fitz', 1., 0.)
        native, reference = [], []
        for energy in energies:
            params.Er = energy
            native.append(rate(A, 32, rho, ct.byref(params)))
            fixture = shared_fixture(ff, isotope, energy, source, map_phi=False)
            reference.append(differential_rate_per_kg_day_keV(
                energy, m_chi_GeV=mchi, rho_chi_GeV_cm3=rho,
                targets_per_kg=KG_IN_GEV/kinematic_mass, dataset=fixture, coefficients=coeff,
                j_chi=.5, config=config, speed_pdf=halo.lab_speed_pdf,
                beta_max=halo.max_lab_speed, breakpoints=halo.integration_breakpoints,
            ))
        # For these two operators, only O11 has an explicit q² multiplying W_M.
        # Keep the same halo kinematics and match that one factor to RAPIDD.
        q2_factor = form_factor_mass/kinematic_mass if op == 11 else 1.
        cases.append(dict(
            A=A, mchi_GeV=mchi, operator=op, endpoint_keV=endpoint,
            rate_mass_GeV=kinematic_mass, form_factor_mass_GeV=form_factor_mass,
            q2_factor=q2_factor, energy_keV=energies.tolist(),
            rapidd_rate=native, wimp_rate_same_mass=reference,
            wimp_rate_matched=(np.asarray(reference)*q2_factor).tolist(),
            rapidd_eta_per_km_s=[eta(v*C_KM_S, 0) for v in vmin],
            wimp_eta_per_km_s=(mean_inverse_speed_shm(halo, vmin)/C_KM_S).tolist(),
        ))
    return cases


def run_absolute_rates(checkout, *, isotopes=ISOTOPES, masses_GeV=MASSES_GEV,
                       operators=OPERATORS, recoil_fractions=RECOIL_FRACTIONS,
                       couplings=COUPLINGS_GEV_MINUS2, halo_km_s=HALO_KM_S,
                       rho_GeV_cm3=RHO_GEV_CM3):
    """Build pristine RAPIDD and return matched, absolute true recoil spectra."""
    checkout = Path(checkout).resolve()
    revision = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'],
                                       text=True).strip()
    status = subprocess.check_output(['git', '-C', str(checkout), 'status', '--porcelain'],
                                     text=True).strip()
    if revision != PIN or status:
        raise ValueError(f'Require a clean RAPIDD checkout at {PIN}.')
    isotopes, operators, masses = tuple(isotopes), tuple(operators), tuple(masses_GeV)
    fractions = np.asarray(tuple(recoil_fractions), dtype=float)
    cp, cn = map(float, couplings)
    v0, vesc, lab = map(float, halo_km_s)
    if (not isotopes or not set(isotopes) <= set(GE_ISOTOPES)
            or not operators or not set(operators) <= set(OPERATORS)):
        raise ValueError('Select stable Ge isotopes and O1/O11.')
    if not masses or not np.all(np.isfinite(masses)) or min(masses) <= 0:
        raise ValueError('Supply finite positive WIMP masses.')
    if (fractions.ndim != 1 or not fractions.size or not np.all(np.isfinite(fractions))
            or np.any(fractions <= 0) or np.any(np.diff(fractions) <= 0)):
        raise ValueError('Recoil fractions must be finite, positive, and strictly increasing.')
    if not np.all(np.isfinite([cp, cn, rho_GeV_cm3])) or rho_GeV_cm3 < 0:
        raise ValueError('Require finite couplings and finite nonnegative density.')
    if not 0 < lab < vesc:
        raise ValueError('The pinned RAPIDD SHM formula requires 0 < v_lab < vesc.')
    halo = StandardHaloModel.from_km_s(v0_km_s=v0, vesc_km_s=vesc, v_lab_km_s=lab)
    with tempfile.TemporaryDirectory(prefix='rapidd-rates-') as tmp:
        cases = _calculate(_build_library(checkout, Path(tmp)), isotopes, masses,
                           operators, fractions, halo, cp, cn, rho_GeV_cm3)
    return dict(schema_version=1, rapidd_commit=revision, basis='pn_Fitz', j_chi=.5,
                halo_km_s=dict(v0=v0, vesc=vesc, v_lab=lab), rho_GeV_cm3=rho_GeV_cm3,
                cp_GeV_minus2=cp, cn_GeV_minus2=cn, higgs_vev_GeV=VEV_GEV,
                units='events/(kg day keV_nr)',
                unit_prefactor_ratio=RAPIDD_ETA_CONVERSION/WIMP_ETA_CONVERSION,
                cases=cases)


def rate_residuals(case, *, relative_floor=1e-6):
    """Return energies, absolute residuals and ratios away from vanishing rates."""
    energy = np.asarray(case['energy_keV'])
    rapidd = np.asarray(case['rapidd_rate'])
    wimp = np.asarray(case['wimp_rate_matched'])
    eta = np.asarray(case['wimp_eta_per_km_s'])
    absolute = rapidd - wimp
    relative = np.full(energy.shape, np.nan)
    eta_relative = relative.copy()
    significant = wimp > relative_floor*np.max(wimp)
    np.divide(absolute, wimp, where=significant, out=relative)
    np.divide(np.asarray(case['rapidd_eta_per_km_s']) - eta, eta,
              where=significant & (eta > 0), out=eta_relative)
    return energy, absolute, relative, eta_relative


def plot_absolute_rates(report):
    """Plot the rates and ratios from any live or saved run."""
    import matplotlib.pyplot as plt

    figures = []
    cases = report['cases']
    unit_ratio = report['unit_prefactor_ratio']
    for operator in sorted({case['operator'] for case in cases}):
        group = [case for case in cases if case['operator'] == operator]
        fig, axes = plt.subplots(2, len(group), squeeze=False, sharex='col',
                                 figsize=(3.3*len(group), 5), layout='constrained')
        for case, (top, bottom) in zip(group, axes.T, strict=True):
            energy, _, residual, _ = rate_residuals(case)
            for name, label, style in (
                ('rapidd_rate', 'RAPIDD', '-'), ('wimp_rate_matched', 'WIMP-LEGEND', '--')
            ):
                rate = np.asarray(case[name])
                visible = (energy < case['endpoint_keV']) & (rate > 0)
                top.plot(energy, np.where(visible, rate, np.nan), style, label=label)
                if np.any(visible):
                    top.set_yscale('log')
            top.set_title(f"Ge-{case['A']}, {case['mchi_GeV']:g} GeV")
            bottom.axhline(0., color='0.5', linewidth=.8, linestyle='--')
            bottom.plot(energy, 1e6*((1 + residual)/unit_ratio - 1), '.-')
            bottom.set_xlabel(r'$E_R$ [keV]')
        axes[0, 0].set_ylabel(r'events/(kg day keV$_{\rm nr}$)')
        axes[1, 0].set_ylabel('Residual after unit ratio [ppm]')
        axes[0, 0].legend()
        fig.suptitle(f'O{operator}, shared RAPIDD nuclear inputs')
        figures.append(fig)
    return figures


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rapidd_checkout', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = run_absolute_rates(args.rapidd_checkout)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
