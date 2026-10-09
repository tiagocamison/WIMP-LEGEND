"""Count-only replaceable statistical backends; exact counting vs asymptotic shape."""
from .base import StatisticsBackend, LimitResult, RootPolicy
from .counting import PoissonCLsCounting
from .binned import BinnedPoissonProfile, FitResult
__all__ = ['StatisticsBackend','LimitResult','RootPolicy','PoissonCLsCounting',
           'BinnedPoissonProfile','FitResult']
