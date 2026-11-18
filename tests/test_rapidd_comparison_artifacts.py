"""Check the recorded independent C comparison, not physical-data acceptance."""
import json
from pathlib import Path
import pytest
from wimp_legend.interactions import CANONICAL_OPERATORS

ROOT = Path(__file__).resolve().parents[1]/'analyses/validation/rapidd'
RESPONSES = json.loads((ROOT/'response_comparison.json').read_text())
KERNELS = json.loads((ROOT/'kernel_comparison.json').read_text())


def test_pinned_comparison_completeness_and_scope():
    assert RESPONSES['rapidd_commit'] == KERNELS['rapidd_commit'] == 'c36b290abfc8f68dc4e53445bc18c185f4b466b0'
    assert len(RESPONSES['records']) == 560
    assert len(KERNELS['records']) == 180
    assert RESPONSES['unsupported'] == ['Phi_tilde_prime']
    assert {r['A'] for r in RESPONSES['records']} == {70,72,73,74,76}
    keys = {(r['A'],r['channel'],r['n1'],r['n2'],r['y']) for r in RESPONSES['records']}
    assert len(keys) == 560


@pytest.mark.parametrize('operator',(1,11))
def test_shared_nuclear_input_o1_o11_kernel_agreement(operator):
    rows = [r for r in KERNELS['records'] if r['operator']==operator]
    assert len(rows) == 60
    for row in rows:
        assert row['rapidd_raw_velocity_kernel'] == pytest.approx(row['ours'],rel=5e-12)


def test_o8_velocity_discrepancy_isolated_from_nuclear_data():
    for row in KERNELS['records']:
        if row['operator'] != 8:
            continue
        assert row['shifted_velocity_diagnostic'] == pytest.approx(row['ours'],rel=5e-12)
        assert row['raw_relative_difference'] > 0
        if row['A'] != 73:  # even isotopes: Delta=0, isolate M term
            mass = 32*.938+(row['A']-32)*.940
            mu = mass*row['mchi_GeV']/(mass+row['mchi_GeV'])
            vmin2 = mass*row['ER_keV']*1e-6/(2*mu**2)
            v2 = (row['v_km_s']/299792.458)**2
            assert row['raw_relative_difference'] == pytest.approx(vmin2/(v2-vmin2),rel=5e-11)


def test_source_particle_number_deficit_not_rescaled():
    for A in (70,72,73,74,76):
        rows=[r for r in RESPONSES['records'] if r['A']==A and r['channel']=='M' and r['y']==0]
        assert sum(r['rapidd'] for r in rows) == pytest.approx(A*A,rel=1e-15)
        assert 0 < sum(r['afh_source'] for r in rows) < A*A


def test_all_operator_snapshot_checks_zero_and_nonzero_kernels():
    reports = json.loads((ROOT/'all_operator_comparison.json').read_text())
    kernels = reports['kernels']
    assert kernels['rapidd_commit'] == KERNELS['rapidd_commit']
    rows = [dict(zip(kernels['columns'], row, strict=True)) for row in kernels['rows']]
    assert len(rows) == 840
    assert {row['operator'] for row in rows} == set(CANONICAL_OPERATORS)
    grid = {
        (A, energy, cp, cn) for A in (70, 72, 73, 74, 76)
        for energy in (1., 10., 30.) for cp, cn in ((1., 1.), (1., 0.), (1., -.7), (0., 1.))
    }
    for op in CANONICAL_OPERATORS:
        group = [row for row in rows if row['operator'] == op]
        assert len(group) == len(grid)
        assert {(r['A'], r['ER_keV'], r['cp'], r['cn']) for r in group} == grid
    for row in rows:
        if row['ours'] == 0:
            assert row['raw_relative_difference'] is None
            assert row['shifted_relative_difference'] is None
            assert row['rapidd_raw_velocity_kernel'] == row['shifted_velocity_diagnostic'] == 0
        else:
            assert row['shifted_velocity_diagnostic'] == pytest.approx(row['ours'], rel=5e-12, abs=0)

    responses = reports['responses']
    assert responses['rapidd_commit'] == kernels['rapidd_commit']
    assert responses['phi_mapping'] == kernels['phi_mapping']
    assert len(responses['rows']) == 640
    rows = (dict(zip(responses['columns'], row, strict=True)) for row in responses['rows'])
    phi_nn = next(r for r in rows if r['channel'] == 'Phi_tilde_prime'
                  and r['A'] == 73 and r['n1'] == r['n2'] == 1 and r['y'] == 0)
    assert phi_nn['rapidd']/phi_nn['afh_source'] == pytest.approx(5., rel=1e-4)
