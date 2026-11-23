import numpy as np
import pytest

from wimp_legend.statistics import (
    BinnedPoissonProfile, GaussianPenalty, NuisanceParameter,
    ProfiledBinnedPoisson, RootPolicy,
)


def test_unconstrained_conditional_and_global_fits():
    model = ProfiledBinnedPoisson(
        lambda mu, t: [mu+t['b'], t['b']],
        [NuisanceParameter('b', 5, 2, (0, None))],
    )
    fixed = model.fit([10, 4], mu=0)
    best = model.fit([10, 4])
    assert fixed.converged and best.converged
    assert fixed.theta_hat['b'] == pytest.approx(7)
    assert best.mu_hat == pytest.approx(6)
    assert best.theta_hat['b'] == pytest.approx(4)
    assert best.deviance == pytest.approx(0, abs=1e-12)
    assert model.q_mu(3, [10, 4]) == 0


def test_penalty_in_deviance_units_including_correlations():
    cov = np.array([[4., 1.], [1., 2.]])
    penalty = GaussianPenalty({'a': 0, 'b': 0}, cov)
    cov[:] = 0
    model = ProfiledBinnedPoisson(lambda mu, t: [mu, t['a']],
                                 [NuisanceParameter('a', 1, bounds=(0, None)),
                                  NuisanceParameter('b')], penalty)
    fit = model.fit([3, 2], mu=3)
    # Conditional on a, b=a/4, reducing the penalty to a²/4.
    # Minimize 2(a-2 log a)+a²/4: a²+4a-8=0.
    a = np.sqrt(12)-2
    assert fit.converged
    assert fit.theta_hat == pytest.approx({'a': a, 'b': a/4}, abs=1e-6)
    assert penalty({'a': 2, 'b': .5}) == pytest.approx(1)
    with pytest.raises(ValueError, match='positive definite'):
        GaussianPenalty({'a': 0, 'b': 0}, [[1, 1], [1, 1]])


def test_independent_gaussian_and_custom_penalty_at_upper_bound():
    assert GaussianPenalty.independent({'x': 2}, {'x': .5})({'x': 3}) == 4
    model = ProfiledBinnedPoisson(lambda mu, t: [mu+t['x']],
                                  [NuisanceParameter('x', .3, bounds=(0, 1))],
                                  lambda t: abs(t['x']-.5)**3)
    fit = model.fit([5], mu=0)
    assert fit.converged and fit.theta_hat['x'] == pytest.approx(1)
    assert fit.deviance == pytest.approx(2*(-4+5*np.log(5))+.125)


def test_unused_nuisance_and_fixed_parameter():
    model = ProfiledBinnedPoisson(lambda mu, t: [mu+t['fixed']],
                                  [NuisanceParameter('unused'),
                                   NuisanceParameter('fixed', 2, bounds=(2, 2))])
    fit = model.fit([4])
    assert fit.converged and fit.locally_flat == ('unused',)
    assert fit.theta_hat['fixed'] == 2
    assert model.observed_upper_limit([4]).mu_upper == pytest.approx(
        BinnedPoissonProfile().observed_upper_limit([1], [2], [4]).mu_upper)


@pytest.mark.parametrize('observed', ([9, 3], [1, 0], [0, 0]))
def test_nuisance_free_observed_parity(observed):
    s, b = np.array([2., 1.]), np.array([2., 1.])
    generic = ProfiledBinnedPoisson(lambda mu, t: mu*s+b)
    known = BinnedPoissonProfile()
    assert generic.observed_upper_limit(observed).mu_upper == pytest.approx(
        known.observed_upper_limit(s, b, observed).mu_upper, abs=1e-7)
    assert generic.q_mu(4, observed) == pytest.approx(known.q_mu(4, s, b, observed))


def test_asimov_fractional_counts_and_constrained_limit():
    s, b = np.array([1., 2.]), np.array([3.5, 1.25])
    known = BinnedPoissonProfile().expected_upper_limit(s, b)
    plain = ProfiledBinnedPoisson(lambda mu, t: mu*s+b)
    result = plain.expected_upper_limit()
    assert result.reference_counts == tuple(b)
    assert result.mu_upper == pytest.approx(known.mu_upper)
    model = ProfiledBinnedPoisson(lambda mu, t: mu*s+b*t['scale'],
                                  [NuisanceParameter('scale', 1, .8, (0, None))],
                                  GaussianPenalty.independent({'scale': 1}, {'scale': .5}))
    assert model.expected_upper_limit().mu_upper > result.mu_upper
    with pytest.raises(ValueError, match='integers'):
        model.fit(b)


def test_impossible_data_and_no_crossing():
    impossible = ProfiledBinnedPoisson(lambda mu, t: [0])
    assert not impossible.fit([1], mu=1).converged
    with pytest.raises(RuntimeError, match='fit failed'):
        impossible.observed_upper_limit([1])
    flat = ProfiledBinnedPoisson(lambda mu, t: [2], root_policy=RootPolicy(max_bracket_steps=4))
    with pytest.raises(RuntimeError, match='bracket'):
        flat.observed_upper_limit([2])


def test_nonconvergence_is_public_and_not_used_for_limits():
    model = ProfiledBinnedPoisson(lambda mu, t: [mu + t['b'], t['b']],
                                  [NuisanceParameter('b', 1, bounds=(0, None))], maxiter=1)
    assert not model.fit([10, 4]).converged
    with pytest.raises(RuntimeError, match='fit failed'):
        model.observed_upper_limit([10, 4])


def test_signal_absorbed_by_nuisance_has_no_finite_upper_crossing():
    model = ProfiledBinnedPoisson(lambda mu, t: [np.exp(mu+t['offset'])],
                                  [NuisanceParameter('offset')],
                                  root_policy=RootPolicy(max_bracket_steps=5))
    with pytest.raises(RuntimeError, match='bracket'):
        model.observed_upper_limit([1])


def test_bad_callback_shape_and_invalid_domain():
    wrong = ProfiledBinnedPoisson(lambda mu, t: [1, 2])
    with pytest.raises(ValueError, match='identical shapes'):
        wrong.fit([1], mu=0)
    excluded = ProfiledBinnedPoisson(lambda mu, t: None if mu < 1 else [mu])
    assert not excluded.fit([1], mu=0).converged
    assert excluded.fit([1]).mu_hat == pytest.approx(1)
