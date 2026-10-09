from dataclasses import FrozenInstanceError
import numpy as np
import pytest
from wimp_legend.backgrounds import (ZeroBackground, FlatBackground, CallableBackground,
                                    BinnedBackground, CompositeBackground)


def test_flat_zero_and_composite_analytic_counts():
    bins=np.array([-.2,.1,.7,2.])
    np.testing.assert_array_equal(ZeroBackground().expected_counts(bins,9),[0,0,0])
    flat=FlatBackground(3)
    np.testing.assert_allclose(flat.expected_counts(bins,9),27*np.diff(bins),rtol=1e-15)
    components=[flat,FlatBackground(2)]
    c=CompositeBackground(components);components.clear()
    np.testing.assert_allclose(c.expected_counts(bins,9),45*np.diff(bins),rtol=1e-15)
    assert len(c.components)==2
    assert np.all(c.expected_counts(bins,0)==0)
    assert np.all(CompositeBackground([]).expected_counts(bins,1)==0)


def test_callable_rate_integral_scalar_array_breakpoint():
    c=CallableBackground(lambda e:1 if e<.5 else 3,breakpoints_keV_obs=[.5])
    np.testing.assert_allclose(c.expected_counts([0,1,2],4),[8,12],rtol=1e-14)
    np.testing.assert_array_equal(c.differential_rate([0,1]),[1,3])
    assert c.differential_rate(.2)==1
    assert c.differential_rate([]).shape==(0,)
    polynomial=CallableBackground(lambda e:2+e*e)
    bins=np.array([0,.3,1.7])
    np.testing.assert_allclose(polynomial.expected_counts(bins,2),2*np.diff(2*bins+bins**3/3),rtol=1e-14)


def test_binned_modes_edges_exposure_and_ownership():
    bins=np.array([0.,1.,3.]);values=np.array([2.,4.])
    c=BinnedBackground(bins,values,'per_unit_exposure')
    f=BinnedBackground(bins,values,'fixed_expected_counts',5,provenance='synthetic fixture')
    bins[:]=9;values[:]=9
    np.testing.assert_array_equal(c.expected_counts([0,1,3],5),[10,20])
    np.testing.assert_array_equal(f.expected_counts([0,1,3],5),[2,4])
    with pytest.raises(ValueError):f.expected_counts([0,1,3],6)
    with pytest.raises(ValueError):c.expected_counts([0,2,3],5)
    with pytest.raises(FrozenInstanceError):c.normalization='fixed_expected_counts'
    with pytest.raises(ValueError):BinnedBackground([0,1],[2],'fixed_expected_counts')
    with pytest.raises(ValueError):BinnedBackground([0,1],[2],'per_unit_exposure',5)
    with pytest.raises(ValueError):BinnedBackground([0,1],[2],'fixed_expected_counts',0)


@pytest.mark.parametrize('bad',[-1,np.nan,np.inf,True,1j,'1'])
def test_bad_rates_and_exposures(bad):
    with pytest.raises(ValueError):FlatBackground(bad)
    with pytest.raises(ValueError):ZeroBackground().expected_counts([0,1],bad)
    with pytest.raises(ValueError):CallableBackground(lambda e:bad).expected_counts([0,1],1)


@pytest.mark.parametrize('values',[[1,2],[-1],[np.nan],[True],[]])
def test_bad_binned_values(values):
    with pytest.raises(ValueError):BinnedBackground([0,1],values,'per_unit_exposure')


def test_background_provenance_and_callback_failures():
    with pytest.raises(ValueError):CallableBackground(lambda e:1,provenance=np.array('mutable'))
    with pytest.raises(ValueError):CallableBackground(lambda e:[1,2]).expected_counts([0,1],1)
    class WrongShape:
        def expected_counts(self,*args):return [1,2]
    with pytest.raises(ValueError):CompositeBackground([WrongShape()]).expected_counts([0,1],1)
