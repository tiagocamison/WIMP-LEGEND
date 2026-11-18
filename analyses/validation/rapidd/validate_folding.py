"""Matched Ge bin counts and fixed-background inference for the RAPIDD benchmark."""

import ctypes as ct
from functools import lru_cache
import runpy

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr, ndtri

from compare_rates import _bind, _RateParameters, HALO_KM_S, RHO_GEV_CM3
from wimp_legend.detector import (ConstantEfficiency, ConstantQuenching, DetectorResponse,
                                  GaussianResolution, NoSmearing)
from wimp_legend.statistics import BinnedPoissonProfile


def run_folding_checks(rpd, checkout, temporary):
    """Use the same stock C Ge spectrum for count and likelihood comparisons."""
    library = rpd.library
    params = _RateParameters(b'Ge', 0., 5., .5, 0., b'pn_Fitz', 1., 0.)
    rpd.coefficients({1: (1e-6, 1e-6)})
    integrate_c = _bind(library, 'counts_bin_noread', ct.c_double, ct.c_double,
                        ct.POINTER(_RateParameters), *([ct.c_double]*3), ct.c_char_p)
    smear_c = _bind(library, 'counts_bin_smeared', ct.c_double, ct.c_char_p, ct.c_double,
                    ct.POINTER(_RateParameters), *([ct.c_double]*3), ct.c_char_p, ct.c_double)
    write_halo = _bind(library, 'write_fv_v', ct.c_int, ct.c_char_p, ct.c_char_p,
                       *([ct.c_double]*5), ct.c_int, *([ct.c_double]*5))
    exposure = 10.
    rho = RHO_GEV_CM3

    @lru_cache(maxsize=4096)
    def spectrum(energy):
        return rpd.ge(rho, ct.byref(params), energy, b'None')

    truth_edges = (.25, .5, 1., 2., 3.)
    native = [integrate_c(rho, ct.byref(params), exposure, lo, hi, b'None')
              for lo, hi in zip(truth_edges[:-1], truth_edges[1:])]
    direct_wimp = DetectorResponse(ConstantQuenching(1.), NoSmearing()).expected_counts(
        spectrum, truth_edges, exposure, (truth_edges[0], truth_edges[-1]))

    # Exercise the claimed native smearing path at deliberately separated widths.
    # The branch under comparison loads its on-disk table inside this function.
    v0, vesc, lab = HALO_KM_S
    halo_file = str(temporary/'halo.dat').encode()
    if write_halo(halo_file, b'SHM', vesc, v0, 0., lab, 0., 2, 0., 0., 0., 0., 0.) != 0:
        raise RuntimeError('Cannot write the RAPIDD halo table.')
    sigmas = (.01, .3)
    c_smeared = [[smear_c(halo_file, rho, ct.byref(params), exposure, lo, hi, b'None', sigma)
                  for lo, hi in zip(truth_edges[:-1], truth_edges[1:])]
                 for sigma in sigmas]

    # Test a controlled Ge response using one shared truth spectrum and an
    # independent CDF integral. It is a benchmark, with no fitted DL/TL inputs.
    quenching, sigma, efficiency, floor = .2, .035, .9, .02
    observed = (.10, .15, .25, .45, .8)
    interval = (.05, 4.5)
    cutoff = floor/quenching
    spectrum.cache_clear()  # C halo state was just read from disk.
    response = DetectorResponse(ConstantQuenching(quenching), GaussianResolution(sigma),
                                truth_efficiency=lambda energy: float(energy >= cutoff),
                                observed_efficiency=ConstantEfficiency(efficiency))
    folded_wimp = response.expected_counts(spectrum, observed, exposure, interval,
                                          truth_breakpoints=(cutoff,))
    folded_reference = [
        exposure*efficiency*quad(lambda energy: spectrum(energy)*(
            ndtr((hi-quenching*energy)/sigma)-ndtr((lo-quenching*energy)/sigma)),
            cutoff, interval[-1], epsrel=1e-7, epsabs=1e-9, limit=200)[0]
        for lo, hi in zip(observed[:-1], observed[1:])
    ]

    # RAPIDD's statistic is the likelihood ratio to mu=0. An Asimov
    # background-only observation puts the WIMP-LEGEND fit at that same point.
    statistic = runpy.run_path(str(checkout/'rapidd/stats.py'))['test_statistic']
    background = np.array([20., 15., 10., 5.])
    model = BinnedPoissonProfile()
    signal = np.asarray(folded_reference)
    strengths = (0., .01, .1, 1.)
    q_rapidd = [statistic(mu, signal, background, background) for mu in strengths]
    q_wimp = [model.q_mu(mu, signal, background, background) for mu in strengths]
    threshold = float(ndtri(.90)**2)

    def limit_at(q_threshold):
        high = 1.
        while statistic(high, signal, background, background) < q_threshold:
            high *= 2
        return brentq(lambda mu: statistic(mu, signal, background, background)-q_threshold, 0., high)

    return dict(exposure_kg_day=exposure, truth_bins_keV=truth_edges, stock_c_bins=native,
                wimp_unsmeared=direct_wimp.tolist(),
                stock_smearing_sigma_keV=sigmas, stock_smearing_counts=c_smeared,
                observed_bins_keV=observed, truth_interval_keV=interval,
                quenching=quenching, resolution_keV=sigma, observed_efficiency=efficiency,
                pre_smearing_floor_keVee=floor, gaussian_reference=folded_reference,
                gaussian_wimp=folded_wimp.tolist(), background=background.tolist(),
                strengths=strengths, rapidd_q=q_rapidd, wimp_q=q_wimp,
                matched_q_threshold=threshold, mu_rapidd_matched=limit_at(threshold),
                mu_rapidd_default_threshold=limit_at(2.706),
                mu_wimp=model.expected_upper_limit(signal, background, confidence_level=.90).mu_upper)
