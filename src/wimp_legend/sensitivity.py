"""Thin count-level analysis orchestration, with no signal-physics assumptions."""
from dataclasses import dataclass
from collections.abc import Mapping
from typing import Callable
import numpy as np
from ._analysis import array, counts, edges, scalar, text
from .backgrounds import BackgroundModel
from .detector import CountResponse
from .statistics import LimitResult, StatisticsBackend
from .statistics.base import confidence, data


@dataclass(frozen=True)
class AnalysisResult:
    signal_counts: tuple[float, ...]
    background_counts: tuple[float, ...]
    observed_bin_edges_keV: tuple[float, ...]
    exposure_kg_day: float
    truth_energy_range_keV: tuple[float, float]
    limit: LimitResult

    def __post_init__(self):
        s,b,_ = data(self.signal_counts,self.background_counts)
        bins = edges(self.observed_bin_edges_keV)
        interval = edges(self.truth_energy_range_keV, 'truth_energy_range_keV', nonnegative=True)
        if len(bins)-1 != len(s) or len(interval) != 2:
            raise ValueError('analysis result has incompatible edges/range/counts')
        if not isinstance(self.limit,LimitResult):
            raise TypeError('backend must return a LimitResult with construction metadata')
        object.__setattr__(self,'signal_counts',tuple(s))
        object.__setattr__(self,'background_counts',tuple(b))
        object.__setattr__(self,'observed_bin_edges_keV',tuple(bins))
        object.__setattr__(self,'truth_energy_range_keV',tuple(interval))
        object.__setattr__(self,'exposure_kg_day',scalar(self.exposure_kg_day,'exposure_kg_day'))


@dataclass(frozen=True)
class GridResult:
    template_id: str
    exposure_kg_day: float
    analysis: AnalysisResult

    def __post_init__(self):
        object.__setattr__(self,'template_id',text(self.template_id,'template_id'))
        object.__setattr__(self,'exposure_kg_day',scalar(self.exposure_kg_day,'exposure_kg_day'))
        if not isinstance(self.analysis,AnalysisResult):
            raise TypeError('analysis must be an AnalysisResult')


@dataclass(frozen=True)
class SensitivityAnalysis:
    response: CountResponse
    background: BackgroundModel
    backend: StatisticsBackend

    def __post_init__(self):
        for model,methods in [(self.response,('expected_counts',)),
                              (self.background,('expected_counts',)),
                              (self.backend,('observed_upper_limit','expected_upper_limit','acceptance_margin'))]:
            if any(not callable(getattr(model,method,None)) for method in methods):
                raise TypeError('analysis component does not implement its protocol')

    def _counts(self,truth_spectrum,bins,exposure,interval,truth_breakpoints):
        observed = edges(bins)
        exposure = scalar(exposure,'exposure_kg_day')
        interval = edges(interval,'truth_energy_range_keV',nonnegative=True)
        if len(interval) != 2:
            raise ValueError('truth_energy_range_keV needs two endpoints')
        s = self.response.expected_counts(truth_spectrum,observed.copy(),exposure,
                                         interval.copy(),truth_breakpoints=truth_breakpoints)
        b = self.background.expected_counts(observed.copy(),exposure)
        s,b,_ = data(s,b)
        if len(s) != len(observed)-1:
            raise ValueError('response/background count shape does not match bins')
        return s,b,observed,exposure,interval

    def observed_limit(self,truth_spectrum,observed_bin_edges_keV,exposure_kg_day,
                       truth_energy_range_keV,observed_counts,*,confidence_level=.90,truth_breakpoints=()):
        s,b,bins,exposure,interval = self._counts(truth_spectrum,observed_bin_edges_keV,
            exposure_kg_day,truth_energy_range_keV,truth_breakpoints)
        _,_,n = data(s,b,observed_counts)
        limit = self.backend.observed_upper_limit(s.copy(),b.copy(),n.copy(),confidence_level=confidence_level)
        return AnalysisResult(tuple(s),tuple(b),tuple(bins),exposure,tuple(interval),limit)

    def expected_limit(self,truth_spectrum,observed_bin_edges_keV,exposure_kg_day,
                       truth_energy_range_keV,*,confidence_level=.90,truth_breakpoints=()):
        s,b,bins,exposure,interval = self._counts(truth_spectrum,observed_bin_edges_keV,
            exposure_kg_day,truth_energy_range_keV,truth_breakpoints)
        limit = self.backend.expected_upper_limit(s.copy(),b.copy(),confidence_level=confidence_level)
        return AnalysisResult(tuple(s),tuple(b),tuple(bins),exposure,tuple(interval),limit)

    def expected_grid(self,truth_spectra,observed_bin_edges_keV,exposures_kg_day,
                      truth_energy_range_keV,*,confidence_level=.90,truth_breakpoints=()):
        """Named supplied templates crossed with supplied exposures; no theory scan."""
        if not isinstance(truth_spectra,Mapping):
            raise TypeError('truth_spectra must map template IDs to callables')
        templates = tuple((text(name,'template_id'),spectrum) for name,spectrum in truth_spectra.items())
        exposures = array(exposures_kg_day,'exposures_kg_day')
        if exposures.ndim != 1:
            raise ValueError('exposures_kg_day must be a vector')
        # Own iterator-valued breakpoints once, so every grid point sees them.
        points = tuple(truth_breakpoints)
        return tuple(GridResult(name,float(exposure),self.expected_limit(
            spectrum,observed_bin_edges_keV,float(exposure),truth_energy_range_keV,
            confidence_level=confidence_level,truth_breakpoints=points))
            for name,spectrum in templates for exposure in exposures)


@dataclass(frozen=True)
class ParameterScan:
    parameter_domain: tuple[float, float]
    parameters: tuple[float, ...]
    acceptance_margins: tuple[float, ...]
    accepted: tuple[bool, ...]
    confidence_level: float
    method: str
    construction: str
    approximation: str
    status: str = 'completed_explicit_grid'
    notes: tuple[str, ...] = (
        'Acceptance is evaluated only at supplied grid points; between-point crossings are unresolved',
        'Accepted sample ranges are not continuous confidence intervals or a unique upper limit',
        'Construction and coverage are defined by the supplied backend acceptance_margin',
    )

    def __post_init__(self):
        domain = edges(self.parameter_domain,'parameter_domain',nonnegative=False)
        grid = edges(self.parameters,'parameters',nonnegative=False)
        margins = array(self.acceptance_margins,'acceptance_margins',nonnegative=False)
        if len(domain) != 2 or grid[0] != domain[0] or grid[-1] != domain[-1] or margins.shape != grid.shape:
            raise ValueError('scan arrays/domain are incompatible')
        accepted = tuple(self.accepted)
        if any(type(x) is not bool for x in accepted) or accepted != tuple(bool(x>=0) for x in margins):
            raise ValueError('accepted flags must agree with margin signs')
        object.__setattr__(self,'parameter_domain',tuple(domain))
        object.__setattr__(self,'parameters',tuple(grid))
        object.__setattr__(self,'acceptance_margins',tuple(margins))
        object.__setattr__(self,'accepted',accepted)
        object.__setattr__(self,'confidence_level',confidence(self.confidence_level))
        for name in ('method','construction','approximation','status'):
            object.__setattr__(self,name,text(getattr(self,name),name))
        object.__setattr__(self,'notes',tuple(text(x,'note') for x in self.notes))

    @property
    def accepted_sample_ranges(self):
        ranges, start = [],None
        for i,accepted in enumerate(self.accepted):
            if accepted and start is None:
                start = i
            if start is not None and (not accepted or i == len(self.accepted)-1):
                end = i if accepted else i-1
                ranges.append((self.parameters[start],self.parameters[end]))
                start = None
        return tuple(ranges)


def scan_parameter(signal_counts: Callable, *, parameter_domain, parameter_grid,
                   background_counts, observed_counts, backend: StatisticsBackend,
                   confidence_level=.90):
    """Explicit nonlinear count-model scan, without assuming monotonicity.

    Each supplied signal_counts(theta) is tested at strength 1 by the backend.
    For BinnedPoissonProfile its alternative fit is along that point's fixed
    signal shape, not a profile fit over theta. Choose an external backend if
    the intended nonlinear model needs another construction.
    """
    if not callable(signal_counts) or any(not callable(getattr(backend,name,None))
            for name in ('acceptance_margin','acceptance_metadata')):
        raise TypeError('signal_counts/backend must implement callable count/test contracts')
    domain = edges(parameter_domain,'parameter_domain',nonnegative=False)
    grid = edges(parameter_grid,'parameter_grid',nonnegative=False)
    if len(domain) != 2 or grid[0] != domain[0] or grid[-1] != domain[-1]:
        raise ValueError('parameter_grid must span the complete explicit parameter_domain')
    b,n = counts(background_counts,'background counts'),counts(observed_counts,'observed counts',observed=True)
    cl = confidence(confidence_level)
    values = []
    for parameter in grid:
        s = counts(signal_counts(float(parameter)),'nonlinear signal counts')
        s,b,_ = data(s,b,n)
        margin = array(backend.acceptance_margin(s.copy(),b.copy(),n.copy(),confidence_level=cl),
                       'acceptance margin',nonnegative=False)
        if margin.ndim != 0:
            raise ValueError('backend acceptance margin must be a finite scalar')
        values.append(float(margin))
    metadata = tuple(backend.acceptance_metadata())
    if len(metadata) != 3:
        raise ValueError('backend acceptance metadata needs method, construction, approximation')
    return ParameterScan(tuple(domain),tuple(grid),tuple(values),tuple(x >= 0 for x in values),cl,*metadata)
