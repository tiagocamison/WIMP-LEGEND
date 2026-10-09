"""Run: PYTHONPATH=src python examples/generic_analysis.py

Entirely synthetic; no Ge response, detector calibration, or experimental data.
The supplied truth callable could instead adapt the existing generic rate API.
"""
from wimp_legend.detector import (
    DetectorResponse, ConstantQuenching, GaussianResolution, ConstantEfficiency,
)
from wimp_legend.backgrounds import FlatBackground
from wimp_legend.statistics import PoissonCLsCounting, BinnedPoissonProfile
from wimp_legend.sensitivity import SensitivityAnalysis, scan_parameter


def truth(E_nr_keV):
    # Synthetic rate in events/(kg day keV_nr), used only within explicit [0,5].
    return 2.0


def main():
    response = DetectorResponse(
        ConstantQuenching(.2), GaussianResolution(.1),
        truth_efficiency=ConstantEfficiency(.4),
        observed_efficiency=ConstantEfficiency(.6),
    )
    background = FlatBackground(3.)  # synthetic events/(kg day keV_obs)
    bins, exposure, truth_range = [0., .2, .6, 1.2], 100., (0., 5.)
    for backend in (PoissonCLsCounting(), BinnedPoissonProfile()):
        analysis = SensitivityAnalysis(response, background, backend)
        result = analysis.expected_limit(truth, bins, exposure, truth_range)
        print(result.limit.method, result.limit.kind)
        print('signal counts:', result.signal_counts)
        print('background counts:', result.background_counts)
        print('mu upper:', result.limit.mu_upper, result.limit.approximation)
        print('status:', result.limit.status, result.limit.notes)

    scan = scan_parameter(
        lambda theta: [5*(theta*theta-1)**2],
        parameter_domain=(-3.,3.), parameter_grid=[-3.,-2.,-1.,0.,1.,2.,3.],
        background_counts=[0.], observed_counts=[0], backend=PoissonCLsCounting(),
    )
    print('nonlinear accepted sample ranges:', scan.accepted_sample_ranges)
    print(scan.notes)


if __name__ == '__main__':
    main()
