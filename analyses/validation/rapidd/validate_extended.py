"""Compare all elastic Ge rates, signed interference and isotope sums to RAPIDD."""

import argparse
import ctypes as ct
from dataclasses import replace
from itertools import product
import json
from math import pi
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from compare_kernels import COUPLINGS, GE_ISOTOPES, shared_fixture
from compare_rates import (_bind, _build_library, _RateParameters, HALO_KM_S,
                           RHO_GEV_CM3, VEV_GEV, RAPIDD_ETA_CONVERSION,
                           WIMP_ETA_CONVERSION, M_N_FORM_FACTOR_GEV, KG_IN_GEV)
from compare_responses import PIN
from wimp_legend.constants import C_KM_S
from wimp_legend.halo import StandardHaloModel, velocity_moment
from wimp_legend.halo.shm_eta import mean_inverse_speed_shm
from wimp_legend.interactions import CANONICAL_OPERATORS, WilsonCoefficients
from wimp_legend.kinematics import maximum_recoil_energy_keV, minimum_speed_c
from wimp_legend.rates import (NREFTRateConfig, differential_rate_per_kg_day_keV,
                               transition_probability)
from wimp_legend.targets import (Isotope, SourceReference, TargetComposition,
                                 load_germanium_afh2014, germanium_afh2014_momentum)

RAPIDD_XI_CONVERSION = 1.889947e3
GE_NUMBER_FRACTIONS = (.2123, .2766, .0773, .3594, .0774)
MIXED_PAIRS = ((1, 3), (4, 5), (4, 6), (8, 9), (11, 12), (11, 15), (12, 15))
ZERO_CONTROLS = ((1, 4), (5, 8))
MIXED_COUPLINGS = (((1., 0.), (1., 0.)), ((0., 1.), (0., 1.)),
                   ((1., 0.), (0., 1.)), ((0., 1.), (1., 0.)))
AMPLITUDE = 1e-6
MASSES = (5., 40.)
RECOIL_FRACTIONS = (.05, .25, .65, .95, 1.01)
BASIS = b'pn_Fitz'
CONFIG = NREFTRateConfig(M_N_FORM_FACTOR_GEV)
SOURCE = SourceReference('RAPIDD '+PIN, 'Shared Ge responses and source constants', 'implementation')


class Rapidd:
    """Read the unmodified C rate and signed form-factor contributions."""
    def __init__(self, library, halo_km_s):
        self.library = library
        self.ff = _bind(library, 'FF_pn_Fitz', ct.c_double, ct.c_int,
                        *([ct.c_char_p]*4), ct.c_double)
        self.mass = _bind(library, 'approx_mass_nucleus', ct.c_double, ct.c_double, ct.c_double)
        self.moment = _bind(library, 'halo', ct.c_double, ct.c_double, ct.c_int)
        self.reset = _bind(library, 'set_coeffs', None)
        self.set = _bind(library, 'set_any_Ncoeff', None, ct.c_double, ct.c_int, ct.c_char_p)
        self.total = _bind(library, 'total_difrate_isotope_dEr', ct.c_double,
                           ct.c_int, ct.c_int, ct.c_double, ct.POINTER(_RateParameters))
        self.signed = _bind(library, 'difrate_isotope_dEr', ct.c_double, ct.c_int, ct.c_int,
                            ct.c_double, ct.POINTER(_RateParameters), ct.c_int, ct.c_int)
        self.ge = _bind(library, 'difrate_dER_2', ct.c_double, ct.c_double,
                        ct.POINTER(_RateParameters), ct.c_double, ct.c_char_p)
        self.form = tuple(_bind(library, f'FormFact_v{order}', ct.c_double, ct.c_char_p,
                                *([ct.c_int]*4), ct.c_char_p, ct.c_char_p, *([ct.c_double]*3))
                          for order in (0, 2))
        v0, vesc, lab = halo_km_s
        setup = _bind(library, 'fv_v', ct.c_int, ct.c_char_p, *([ct.c_double]*5),
                      ct.c_int, *([ct.c_double]*5))
        if setup(b'SHM', vesc, v0, 0., lab, 0., 2, 0., 0., 0., 0., 0.) != 0:
            raise RuntimeError('RAPIDD halo initialization failed.')

    def coefficients(self, couplings):
        self.reset()
        for op, (cp, cn) in couplings.items():
            self.set(cp*VEV_GEV**2, op, b'p')
            self.set(cn*VEV_GEV**2, op, b'n')

    def isotope(self, A):
        return Isotope('Ge', A, 32, self.mass(A, 32), 4.5 if A == 73 else 0.,
                       {key: SOURCE for key in ('identity', 'mass_GeV', 'spin')})

    def components(self, A, params, couplings, orders):
        """F0 and F2 using the C entries in the order used by its public rate."""
        return tuple(sum(cn*cm*form(BASIS, A, 32, i, j, n, m, params.Er, params.mchi, .5)
                         for i, j in orders
                         for n, cn in zip((b'p', b'n'), couplings[i])
                         for m, cm in zip((b'p', b'n'), couplings[j]))
                     for form in self.form)


def _parameters(energy, mass):
    return _RateParameters(b'Ge', energy, mass, .5, 0., BASIS, 1., 0.)


def _q2(energy, A):
    # The form-factor mass differs from the mass in RAPIDD's v_min.
    return 2*(32*.938+(A-32)*.940)*energy*1e-6


def _kernel(fixture, q2, coeff, velocity):
    return transition_probability(fixture, coeff, q2_GeV2=q2, V=velocity, j_chi=.5, config=CONFIG)


def _wimp_rate(fixture, q2, coeff, vmin, eta, xi, mass):
    # Eq. (40) is affine in V. Integrate the constant and V terms separately.
    constant, at_one = _kernel(fixture, q2, coeff, (0., 1.))
    return RHO_GEV_CM3/(2*pi*mass)*WIMP_ETA_CONVERSION*(
        constant*eta+(at_one-constant)*(xi-vmin**2*eta)/C_KM_S**2)


def _rescaled_c(terms, eta, xi, vmin, mass):
    constant, velocity = terms
    return RHO_GEV_CM3/(2*pi*mass)*WIMP_ETA_CONVERSION*(
        constant*eta+velocity*(xi-vmin**2*eta)/C_KM_S**2)


def _stock_c_from_terms(terms, eta, xi, mass):
    constant, velocity = terms
    return RHO_GEV_CM3/(2*pi*mass)*(
        constant*eta*RAPIDD_ETA_CONVERSION+velocity*xi*RAPIDD_XI_CONVERSION)


def _moments(rpd, halo, isotope, energies, mass):
    beta = minimum_speed_c(energies, mass, isotope.mass_GeV)
    speed = beta*C_KM_S
    eta = mean_inverse_speed_shm(halo, beta)/C_KM_S
    xi = velocity_moment(halo.lab_speed_pdf, halo.max_lab_speed, beta, n=1,
                         breakpoints=halo.integration_breakpoints)*C_KM_S
    return zip(speed, eta, xi, (rpd.moment(v, 0) for v in speed),
               (rpd.moment(v, 2) for v in speed), strict=True)


def _metric(left, right):
    return abs(left-right)/max(abs(left), abs(right)) if left or right else 0.


def _summary(records, operator_columns):
    result = []
    groups = sorted({tuple(row[key] for key in operator_columns) for row in records})
    for group in groups:
        rows = [row for row in records if tuple(row[key] for key in operator_columns) == group]
        worst = max(rows, key=lambda row: _metric(row['rapidd_aligned'], row['wimp']))
        raw_worst = max(rows, key=lambda row: _metric(row['rapidd_raw'], row['wimp']))
        summary = dict(zip(operator_columns, group), cases=len(rows),
                           max_raw_gap=_metric(raw_worst['rapidd_raw'], raw_worst['wimp']),
                           max_aligned_gap=_metric(worst['rapidd_aligned'], worst['wimp']),
                           max_c_reconstruction_gap=max(r['c_gap'] for r in rows),
                           worst=worst, worst_raw=raw_worst)
        if operator_columns == ('operator',):
            summary['negative_signed_count'] = sum(r['rapidd_signed'] < 0 for r in rows)
        if group == (8, 9):
            summary['max_reverse_order_gap'] = max(
                _metric(r['rapidd_reversed_aligned'], r['wimp']) for r in rows)
        result.append(summary)
    return result


def _isotope_checks(rpd, halo):
    singles, mixed = [], []
    for A, mass in product(GE_ISOTOPES, MASSES):
        isotope = rpd.isotope(A)
        afh = load_germanium_afh2014(isotope, momentum=germanium_afh2014_momentum(A))
        endpoint = maximum_recoil_energy_keV(mass, isotope.mass_GeV, halo.max_lab_speed)
        energies = endpoint*np.asarray(RECOIL_FRACTIONS)
        for energy, (vmin, eta, xi, c_eta, c_xi) in zip(
                energies, _moments(rpd, halo, isotope, energies, mass), strict=True):
            fixture = shared_fixture(rpd.ff, isotope, energy, SOURCE, map_phi=True)
            q2 = _q2(energy, A)
            params = _parameters(energy, mass)
            point = dict(A=A, mchi_GeV=mass, recoil_keV=energy)
            for op, pair in product(CANONICAL_OPERATORS, COUPLINGS):
                coupling = {op: tuple(AMPLITUDE*c for c in pair)}
                coeff = WilsonCoefficients.from_proton_neutron(coupling)
                rpd.coefficients(coupling)
                raw = rpd.total(A, 32, RHO_GEV_CM3, ct.byref(params))
                signed = rpd.signed(A, 32, RHO_GEV_CM3, ct.byref(params), op, op)
                terms = rpd.components(A, params, coupling, ((op, op),))
                row = dict(point, operator=op, couplings=pair, rapidd_raw=raw,
                                    rapidd_signed=signed,
                                    rapidd_aligned=_rescaled_c(terms, c_eta, c_xi, vmin, mass),
                                    wimp=_wimp_rate(fixture, q2, coeff, vmin, eta, xi, mass),
                                    c_gap=max(_metric(raw, abs(signed)),
                                              _metric(signed, _stock_c_from_terms(
                                                  terms, c_eta, c_xi, mass))))
                if signed < 0:
                    row['wimp_afh_at_negative'] = _wimp_rate(afh, q2, coeff, vmin, eta, xi, mass)
                singles.append(row)
            for (i, j), (first, second) in product(MIXED_PAIRS+ZERO_CONTROLS, MIXED_COUPLINGS):
                plus = {i: tuple(AMPLITUDE*c for c in first), j: tuple(AMPLITUDE*c for c in second)}
                minus = plus.copy()  # Flip the second operator in both rate calculators.
                minus[j] = tuple(-c for c in plus[j])
                rpd.coefficients(plus)
                total_plus = rpd.total(A, 32, RHO_GEV_CM3, ct.byref(params))
                orders = ((i, j), (j, i)) if (i, j) == (4, 6) else ((i, j),)
                signed = sum(rpd.signed(A, 32, RHO_GEV_CM3, ct.byref(params), a, b)
                             for a, b in orders)
                terms = rpd.components(A, params, plus, orders)
                reverse = _rescaled_c(rpd.components(A, params, plus, ((j, i),)),
                                      c_eta, c_xi, vmin, mass)
                rpd.coefficients(minus)
                total_minus = rpd.total(A, 32, RHO_GEV_CM3, ct.byref(params))
                w_plus, w_minus = (_wimp_rate(fixture, q2,
                    WilsonCoefficients.from_proton_neutron(c), vmin, eta, xi, mass)
                    for c in (plus, minus))
                raw = (total_plus-total_minus)/2
                scale = max(total_plus, total_minus)
                reconstruction = _stock_c_from_terms(terms, c_eta, c_xi, mass)
                mixed.append(dict(point, first=i, second=j, couplings=(first, second),
                                  rapidd_raw=raw, rapidd_signed=signed,
                                  rapidd_aligned=_rescaled_c(terms, c_eta, c_xi, vmin, mass),
                                  rapidd_reversed_aligned=reverse, wimp=(w_plus-w_minus)/2,
                                  c_gap=max(abs(raw-signed), abs(signed-reconstruction))/scale
                                  if scale else 0.))
    return singles, mixed


def _mixture_checks(rpd, halo):
    isotopes = tuple(map(rpd.isotope, GE_ISOTOPES))
    comp = TargetComposition.from_reported(
        tuple(zip(isotopes, GE_NUMBER_FRACTIONS)), 'number', SOURCE, normalize=True)
    normalized = np.array([weight for _, weight in comp.to_number_fractions()])
    mass_weights = np.array([w for _, w in comp.to_mass_fractions()])
    energies = (.1, .5, 1., 3., 10., 30.)
    native, table_impact = [], []
    for mass, op in product(MASSES, CANONICAL_OPERATORS):
        coupling = {op: (AMPLITUDE, AMPLITUDE)}
        coeff = WilsonCoefficients.from_proton_neutron(coupling)
        rpd.coefficients(coupling)
        spectra = []
        for isotope, velocities in zip(isotopes, (_moments(rpd, halo, iso, energies, mass)
                                                  for iso in isotopes), strict=True):
            A = isotope.A
            afh = load_germanium_afh2014(isotope, momentum=germanium_afh2014_momentum(A))
            for energy, (vmin, eta, xi, _, _) in zip(energies, velocities, strict=True):
                params = _parameters(energy, mass)
                raw = rpd.total(A, 32, RHO_GEV_CM3, ct.byref(params))
                fixture = shared_fixture(rpd.ff, isotope, energy, SOURCE, map_phi=True)
                q2 = _q2(energy, A)
                wimp_c = _wimp_rate(fixture, q2, coeff, vmin, eta, xi, mass)
                wimp_afh = _wimp_rate(afh, q2, coeff, vmin, eta, xi, mass)
                # Compare tables at a common q² and halo. Keep only the candidate
                # response from RAPIDD to isolate its effect on O12 and O13.
                if A == 73 and op in (12, 13):
                    hybrid = replace(fixture, responses={
                        key: value if key.channel == 'Phi_tilde_prime' else
                        replace(value, coefficients=(afh.evaluate(channel=key.channel,
                            tau=key.tau, tau_prime=key.tau_prime, q2_GeV2=q2),))
                        for key, value in fixture.responses.items()
                    })
                    table_impact.append(dict(operator=op, mchi_GeV=mass, recoil_keV=energy,
                        afh=wimp_afh, rapidd_table=wimp_c,
                        candidate_only=_wimp_rate(hybrid, q2, coeff, vmin, eta, xi, mass)))
                spectra.append((raw, wimp_c, wimp_afh))
        rates = np.asarray(spectra).reshape(len(isotopes), len(energies), 3)
        for energy, per_iso in zip(energies, rates.transpose(1, 0, 2), strict=True):
            params = _parameters(energy, mass)
            native.append(dict(operator=op, mchi_GeV=mass, recoil_keV=energy,
                rapidd_ge=rpd.ge(RHO_GEV_CM3, ct.byref(params), energy, b'None'),
                rapidd_number_sum=np.dot(GE_NUMBER_FRACTIONS, per_iso[:, 0]),
                rapidd_normalized_number_sum=np.dot(normalized, per_iso[:, 0]),
                rapidd_mass_sum=np.dot(mass_weights, per_iso[:, 0]),
                wimp_matched_number_sum=np.dot(GE_NUMBER_FRACTIONS, per_iso[:, 1]),
                wimp_matched_mass_sum=np.dot(mass_weights, per_iso[:, 1]),
                wimp_afh_mass_sum=np.dot(mass_weights, per_iso[:, 2])))
    return native, table_impact, normalized.tolist(), mass_weights.tolist()


def _check_moment_factorization(rpd, halo):
    # Compare the moment reduction with the existing speed-quadrature path
    # using a common q². The RAPIDD two-mass diagnostic is a separate operation.
    isotope = rpd.isotope(73)
    energy, mass = 3., 40.
    fixture = shared_fixture(rpd.ff, isotope, energy, SOURCE, map_phi=True)
    q2 = 2*isotope.mass_GeV*energy*1e-6
    vmin, eta, xi, _, _ = next(_moments(rpd, halo, isotope, (energy,), mass))
    records = []
    for operators in ((8,), (13,), (8, 9), (12, 15)):
        coeff = WilsonCoefficients.from_proton_neutron(
            {op: (AMPLITUDE, -.7*AMPLITUDE) for op in operators})
        factored = _wimp_rate(fixture, q2, coeff, vmin, eta, xi, mass)
        direct = differential_rate_per_kg_day_keV(
            energy, m_chi_GeV=mass, rho_chi_GeV_cm3=RHO_GEV_CM3,
            targets_per_kg=KG_IN_GEV/isotope.mass_GeV, dataset=fixture, coefficients=coeff,
            j_chi=.5, config=CONFIG, speed_pdf=halo.lab_speed_pdf,
            beta_max=halo.max_lab_speed, breakpoints=halo.integration_breakpoints)
        records.append(dict(operators=operators, factored=factored, quadrature=direct))
    return records


def run_extended(checkout):
    """Run the signed-kernel, isotope, detector-folding and inference checks."""
    from validate_folding import run_folding_checks

    checkout = Path(checkout).resolve()
    revision = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'],
                                       text=True).strip()
    if revision != PIN or subprocess.check_output(
            ['git', '-C', str(checkout), 'status', '--porcelain'], text=True).strip():
        raise ValueError(f'Require clean RAPIDD checkout {PIN}.')
    v0, vesc, lab = HALO_KM_S
    halo = StandardHaloModel.from_km_s(v0_km_s=v0, vesc_km_s=vesc, v_lab_km_s=lab)
    with tempfile.TemporaryDirectory(prefix='rapidd-extended-') as temp:
        temp = Path(temp)
        rpd = Rapidd(_build_library(checkout, temp), HALO_KM_S)
        singles, mixed = _isotope_checks(rpd, halo)
        mixtures, candidate, number_weights, mass_weights = _mixture_checks(rpd, halo)
        moment_checks = _check_moment_factorization(rpd, halo)
        folding = run_folding_checks(rpd, checkout, temp)
    return dict(schema_version=1, rapidd_commit=revision, halo_km_s=HALO_KM_S,
                units='events/(kg day keV_nr)', rho_GeV_cm3=RHO_GEV_CM3, j_chi=.5,
                coupling_amplitude_GeV_minus2=AMPLITUDE, masses_GeV=MASSES,
                recoil_endpoint_fractions=RECOIL_FRACTIONS, diagonal_couplings=COUPLINGS,
                interference_couplings=MIXED_COUPLINGS,
                singles=_summary(singles, ('operator',)), mixed=_summary(mixed, ('first', 'second')),
                mixtures=mixtures, phi_candidate=candidate, moment_checks=moment_checks,
                folding=folding,
                isotope_A=GE_ISOTOPES, raw_ge_weights=GE_NUMBER_FRACTIONS,
                normalized_number_weights=number_weights, mass_weights=mass_weights)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rapidd_checkout', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run_extended(args.rapidd_checkout), indent=2,
                                      allow_nan=False)+'\n')
