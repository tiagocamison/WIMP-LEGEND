"""Focused permanent regressions for CODE_AUDIT_001, synthetic inputs only."""
from dataclasses import replace
from decimal import Decimal, localcontext

import numpy as np
import pytest
from scipy.integrate import quad

from wimp_legend.targets import (
    Isotope, SourceReference, TargetComposition, MomentumVariable, ResponseKey,
    ResponseMetadata, NuclearResponseDataset,
)
from wimp_legend.interactions import WilsonCoefficients, particle_response
from wimp_legend.kinematics import momentum_transfer_GeV, minimum_speed_c
from wimp_legend.halo import StandardHaloModel, build_velocity_integral


def records():
    source = SourceReference('synthetic audit regression', 'test_audit_regressions.py', 'synthetic')
    iso = Isotope('X', 20, 4, 8, 0, dict.fromkeys(('identity', 'mass_GeV', 'spin'), source))
    meta = ResponseMetadata('synthetic', ('0', '1'), 'synthetic labels', 'raw',
                            'dimensionless', MomentumVariable('q2_GeV2'), source, 'test-v1 ')
    return source, iso, meta


@pytest.mark.parametrize('field', ['kind', 'reference', 'locator', 'notes', 'transformation'])
def test_audit001_source_text_rejects_mutable_arrays(field):
    s, _, _ = records()
    if field == 'transformation':
        s = SourceReference('derived', 'conversion', 'implementation', derived_from=s,
                            transformation='number_to_mass', input_fractions=(.5, .5))
    with pytest.raises(ValueError):
        replace(s, **{field: np.array(getattr(s, field))})


@pytest.mark.parametrize('field', ['basis', 'isospin_definition', 'normalization',
                                   'response_units', 'convention_id'])
def test_audit001_metadata_text_rejects_mutable_arrays(field):
    _, _, meta = records()
    with pytest.raises(ValueError):
        replace(meta, **{field: np.array(getattr(meta, field))})


def test_audit001_momentum_and_fraction_tags():
    s, iso, _ = records()
    with pytest.raises(ValueError):
        MomentumVariable(np.array('y'), 4, s)
    pairs = [(iso, .5), (replace(iso, A=21, mass_GeV=16), .5)]
    for normalize in (False, True):
        with pytest.raises(ValueError):
            TargetComposition.from_reported(pairs, np.array('number'), s, normalize=normalize)
    with pytest.raises(ValueError):
        TargetComposition(pairs, np.array('mass'), s)
    c = TargetComposition(pairs, 'number', s)
    assert [f for _, f in c.to_mass_fractions()] == pytest.approx([1/3, 2/3])


def test_audit001_other_text_and_exact_owned_ids():
    s, iso, meta = records()
    for call in [lambda: ResponseKey(np.array('M'), 0, 0),
                 lambda: replace(iso, symbol=np.array('X')),
                 lambda: replace(s, upstream=[np.array('citation')]),
                 lambda: replace(meta, isospin_labels=[np.array('0'), '1']),
                 lambda: NuclearResponseDataset(np.array('id'), iso, meta, {})]:
        with pytest.raises(ValueError):
            call()
    for field in ('basis', 'isospin_definition', 'normalization', 'response_units', 'convention_id'):
        owned = replace(meta, **{field: np.str_(getattr(meta, field))})
        assert type(getattr(owned, field)) is str
        assert getattr(owned, field) == getattr(meta, field)
    assert meta.convention_id == 'test-v1 '


@pytest.mark.parametrize('offset', [1e-12, -1e-12, 0])
def test_audit002_factored_sigma_high_precision(offset):
    q = .001
    raw = {4: (-q + offset, -q - 2e-12), 6: (1, 1)}
    c = WilsonCoefficients(raw)
    for t, u in ((0, 0), (0, 1), (1, 0), (1, 1)):
        with localcontext() as ctx:
            ctx.prec = 80
            expected = float((Decimal(raw[4][t])+Decimal(q)) *
                             (Decimal(raw[4][u])+Decimal(q)) / 16)
        got = particle_response('Sigma_double_prime', c, Q=q, V=0, j_chi=.5, tau=t, tau_prime=u)
        assert got == pytest.approx(expected, rel=2e-15, abs=0)
        if t == u:
            assert got >= 0
        elif offset > 0:
            assert got < 0


@pytest.mark.parametrize('energy,mass', [(1e-300, 1e308), (1e-200, 1e-200), (1, 50), (1e308, 1e308)])
def test_audit003_scaled_momentum_and_speed(energy, mass):
    with localcontext() as ctx:
        ctx.prec = 80
        q = (2 * Decimal(mass) * Decimal(energy) * Decimal('1e-6')).sqrt()
        # Equal WIMP/target masses imply 2 mu = mass.
        speed = q / Decimal(mass)
    assert momentum_transfer_GeV(energy, mass) == pytest.approx(float(q), rel=7e-16, abs=0)
    assert minimum_speed_c(energy, mass, mass) == pytest.approx(float(speed), rel=7e-16, abs=0)


def test_audit003_nonrepresentable_speed_rejected():
    with pytest.raises((ValueError, FloatingPointError)):
        minimum_speed_c(1, 1e-320, 1)


@pytest.mark.parametrize('boost', [0, 1e-12, 1, 2, 3])
def test_audit004_dimensionless_pdf_scaling_and_normalization(boost):
    tiny = StandardHaloModel(1e-100, 2e-100, boost*1e-100)
    ordinary = StandardHaloModel(.001, .002, boost*.001)
    x = np.linspace(0, 2+boost, 151, endpoint=False)
    for name in ('lab_speed_pdf', 'galactic_speed_pdf'):
        np.testing.assert_allclose(getattr(tiny, name)(x*1e-100)*1e-100,
                                   getattr(ordinary, name)(x*.001)*.001,
                                   rtol=2e-12, atol=1e-15)
    points = [v/tiny.max_lab_speed for v in tiny.integration_breakpoints]
    norm = quad(lambda x: float(tiny.lab_speed_pdf(x*tiny.max_lab_speed))*tiny.max_lab_speed,
                0, 1, points=points, epsabs=1e-12)[0]
    assert norm == pytest.approx(1, rel=1e-11)


@pytest.mark.parametrize('length', [1e-160, 1e-100, .003])
@pytest.mark.parametrize('power', [-1, 0])
def test_audit005_dimensionless_table(length, power):
    pdf = lambda v: 3*(np.asarray(v)/length)**2/length
    table = build_velocity_integral(pdf, length, power=power, n_points=200,
                                    breakpoints=[.31*length])
    a = np.array([[0., .2], [.5, 1.]])
    expected = 3*length**power/(power+3)*(1-a**(power+3))
    np.testing.assert_allclose(table(length*a), expected, rtol=1e-4, atol=0)
    assert table(np.empty((0, 2))).shape == (0, 2)
    assert table(np.inf) == 0


def test_audit005_unrepresentable_table_fails():
    length = .003
    with pytest.raises((ValueError, OverflowError, FloatingPointError)):
        build_velocity_integral(lambda v: 1e308*np.minimum((v/(.1*length))**2, 1),
                                length, n_points=100)


@pytest.mark.parametrize('value', [np.nextafter(0., 1.), 1e-320, 1e-200, 1., 1e308])
def test_audit006_equal_and_opposite_range(value):
    value = float(value)
    for sign in (1, -1):
        assert WilsonCoefficients.from_proton_neutron({1: (value, sign*value)}).values[1] == (
            (value, 0.) if sign == 1 else (0., value))
