"""Pinned AFH PRC supplement Ge responses; no abundances or nuclear masses.

Anand, Fitzpatrick, Haxton, PRC 89 (2014) 065501, supplement SHA-256
943a403e216633840887bbfcbec55c56e0c5813b86335b49868acb89e74a49b9.
See docs/GE_NUCLEAR_DATA.md for extraction, ordering, precision and limitations.
Only explicit, sourced ``y`` momentum mappings are accepted. This importer
performs no runtime convention conversion; the versioned data is canonical.
"""

import hashlib
from importlib.resources import files
import json
from math import sqrt

from .isotopes import Isotope
from .nuclear_responses import (
    MomentumVariable, NuclearResponseDataset, PolynomialResponse,
    ResponseKey, ResponseMetadata,
)
from .provenance import SourceReference, _integer

_VERSION = "afh2014_ge_obdme_v1"
_MANIFEST_SHA256 = "b17df65326c24a46a4effef6f2b5881304c9c0f03792417c2635e921a6a76f33"
_CONVENTION = "anand2014_prc89_065501_half_isospin_v1"
_CHANNELS = ("M", "Sigma_double_prime", "Sigma_prime", "Phi_double_prime",
             "Phi_tilde_prime", "Delta", "Phi_double_prime_M", "Delta_Sigma_prime")


def _read(A):
    A = _integer(A, "A")
    if A not in (70, 72, 73, 74, 76):
        raise ValueError("AFH Ge data supports A=70,72,73,74,76 only")
    root = files("wimp_legend.targets").joinpath("data", "ge_afh2014_v1")
    manifest_bytes = root.joinpath("manifest.json").read_bytes()
    if hashlib.sha256(manifest_bytes).hexdigest() != _MANIFEST_SHA256:
        raise ValueError("AFH Ge manifest checksum mismatch")
    manifest = json.loads(manifest_bytes)
    entry = manifest["isotopes"][str(A)]
    raw = root.joinpath(entry["file"]).read_bytes()
    if hashlib.sha256(raw).hexdigest() != entry["sha256"]:
        raise ValueError("AFH Ge dataset checksum mismatch")
    data = json.loads(raw)
    if (data["schema_version"] != 1 or data["dataset_version"] != _VERSION
            or data["convention_id"] != _CONVENTION):
        raise ValueError("incompatible AFH Ge dataset schema/version/convention")
    return data, entry["sha256"]


def germanium_afh2014_momentum(A: int) -> MomentumVariable:
    """Explicitly select the supplement's oscillator/unit prescription for A.

    b_fm=sqrt(41.467/(45*A**(-1/3)-25*A**(-2/3))), SetIsotope line 1479;
    b_GeV_inv=5.0677*b_fm, Femtometer line 17. These are source numerical
    inputs, not universal constants or empirical validity bounds. Callers may
    instead construct a sourced MomentumVariable('y', ...) themselves.
    """
    data, digest = _read(A)
    A = data["isotope"]["A"]
    b_fm = sqrt(41.467/(45*A**(-1/3)-25*A**(-2/3)))
    source = SourceReference(
        reference="Anand, Fitzpatrick, Haxton, PRC 89, 065501 (2014), supplement",
        locator="SetIsotope line 1479; Femtometer line 17",
        kind="primary_literature",
        notes=f"Explicit supplement b(A) prescription; b_fm={b_fm!r}; "
              f"rounded source fm=5.0677 GeV^-1; dataset_sha256={digest}; "
              f"supplement_sha256={data['primary_source']['supplement_sha256']}",
    )
    return MomentumVariable("y", b_fm*5.0677, source)


def load_germanium_afh2014(isotope: Isotope, *, momentum: MomentumVariable) -> NuclearResponseDataset:
    """Load all 32 canonical W entries for a supplied Ge ground-state isotope.

    ``isotope`` retains the caller's explicitly sourced mass. Identity and spin
    must match the supplement. ``momentum`` must explicitly select a sourced y
    mapping (e.g. germanium_afh2014_momentum(isotope.A)); no hidden b(A) default.
    The analytic mathematical domain is y>=0, unbounded. No physical validated
    range is known. Retains the source particle-number deficit, never A²-fixes.
    """
    if not isinstance(isotope, Isotope):
        raise TypeError("isotope must be an Isotope with sourced mass and spin")
    if isotope.symbol != "Ge" or isotope.Z != 32:
        raise ValueError("AFH Ge loader requires germanium, Z=32")
    data, digest = _read(isotope.A)
    expected_spin = 4.5 if isotope.A == 73 else 0.0
    if isotope.spin != expected_spin:
        raise ValueError("isotope spin disagrees with the AFH ground state")
    if not isinstance(momentum, MomentumVariable) or momentum.variable != "y":
        raise ValueError("AFH polynomials require an explicit sourced y momentum mapping")
    provenance = {k: v for k, v in data.items()
                  if k not in ("responses", "source_pn_polynomials", "density_rows")}
    provenance["data_sha256"] = digest
    source = SourceReference(
        reference="Anand, Fitzpatrick, Haxton, PRC 89, 065501 (2014), DOI:10.1103/PhysRevC.89.065501",
        locator=f"pfGe{isotope.A}; dmformfactor-prc.m lines 803–1031; analytic export {_VERSION}/ge{isotope.A}.json",
        kind="implementation",
        upstream=("Anand et al., PRC 89 (2014) 065501, DOI:10.1103/PhysRevC.89.065501, pinned primary supplement",
                  "Fitzpatrick et al., JCAP 02 (2013) 004, arXiv:1203.3542v3, Appendix A.3"),
        notes=json.dumps(provenance, sort_keys=True),
    )
    responses = {}
    for entry in data["responses"]:
        key = ResponseKey(entry["channel"], entry["tau"], entry["tau_prime"])
        if key in responses:
            raise ValueError("duplicate response in AFH export")
        coefficient_source = SourceReference(
            reference=source.reference, kind="implementation", upstream=source.upstream,
            locator=f"{source.locator}; source response {entry['source_response_index']}, "
                    f"ordered source isospin {entry['source_isospin_indices']}; "
                    f"canonical {key.channel}[{key.tau},{key.tau_prime}], ascending powers",
            notes=f"data_sha256={digest}; {entry['verification_status']}; "
                  "exact written-decimal OBDMEs; exp(-2y); factor 4; transpose mixed source channels 7/8; "
                  "no particle-number rescaling; mathematical domain y>=0; physical validated range unknown",
        )
        responses[key] = PolynomialResponse(tuple(map(float, entry["coefficients"])),
                                             2.0, None, coefficient_source)
    expected = {ResponseKey(ch, t, u) for ch in _CHANNELS for t in (0, 1) for u in (0, 1)}
    if set(responses) != expected:
        raise ValueError("incomplete AFH Ge response dataset")
    metadata = ResponseMetadata(
        basis="Anand–Fitzpatrick–Haxton one-body Galilean nuclear responses, PRC 89 (2014) Eq. (41)",
        isospin_labels=("0", "1"),
        isospin_definition="t0=1, t1=tau3; proton +1, neutron -1; c0=(cp+cn)/2, c1=(cp-cn)/2",
        normalization="raw reduced-element bilinears W; no spin average, 4pi, external q2/m_N2, "
                      "interference multiplicity, Wilson coefficients or abundance absorbed; "
                      "4pi/(2J+1) belongs to the canonical contraction",
        response_units="dimensionless", momentum=momentum, source=source, convention_id=_CONVENTION,
    )
    return NuclearResponseDataset(f"Ge{isotope.A}_{_VERSION}_b{momentum.b_GeV_inv.hex()}",
                                  isotope, metadata, responses)
