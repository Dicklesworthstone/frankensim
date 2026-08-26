"""
FrankenSim Python SDK - Structured Exceptions.

Exit codes and error conditions from the Rust engine are mapped to explicit,
actionable Python exception classes.
"""

from __future__ import annotations
from typing import List, Optional
from .models import Diagnostic, ExitCode


class FrankenSimError(Exception):
    """Base class for all FrankenSim exceptions."""
    def __init__(
        self,
        message: str,
        exit_code: int = ExitCode.REFUSED,
        diagnostics: Optional[List[Diagnostic]] = None,
    ):
        super().__init__(message)
        self.message = message
        self.exit_code = exit_code
        self.diagnostics = diagnostics or []

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(exit_code={self.exit_code}, message={self.message!r})"


class UsageError(FrankenSimError):
    """Command arguments or flags do not match the documented grammar (exit 2)."""
    def __init__(self, message: str, diagnostics: Optional[List[Diagnostic]] = None):
        super().__init__(message, exit_code=ExitCode.USAGE, diagnostics=diagnostics)


class InputError(FrankenSimError):
    """Input file could not be read, decoded, or admitted within resource caps (exit 3)."""
    def __init__(self, message: str, diagnostics: Optional[List[Diagnostic]] = None):
        super().__init__(message, exit_code=ExitCode.INPUT, diagnostics=diagnostics)


class RefusalError(FrankenSimError):
    """Input was read but refused by schema, semantic checks, or physics contracts (exit 4)."""
    def __init__(self, message: str, diagnostics: Optional[List[Diagnostic]] = None):
        super().__init__(message, exit_code=ExitCode.REFUSED, diagnostics=diagnostics)


class UnavailableError(FrankenSimError):
    """Requested command or stage is unavailable pending producing capability (exit 5)."""
    def __init__(self, message: str, diagnostics: Optional[List[Diagnostic]] = None):
        super().__init__(message, exit_code=ExitCode.UNAVAILABLE, diagnostics=diagnostics)


class BudgetExceededError(FrankenSimError):
    """Solver stopped honestly at a budget limit with a durable prefix (exit 6)."""
    def __init__(self, message: str, diagnostics: Optional[List[Diagnostic]] = None):
        super().__init__(message, exit_code=ExitCode.BUDGET, diagnostics=diagnostics)


class CancellationError(FrankenSimError):
    """Execution was cancelled via cancellation gate or SIGINT (exit 130)."""
    def __init__(self, message: str, diagnostics: Optional[List[Diagnostic]] = None):
        super().__init__(message, exit_code=ExitCode.CANCELLED, diagnostics=diagnostics)


class FrankenSimTimeoutError(FrankenSimError):
    """The CLIENT-imposed wall-clock timeout killed the engine (SIGKILL).

    Distinct from :class:`BudgetExceededError`, which is the ENGINE's honest
    exit-6 budget stop with a durable prefix. A killed child produced no
    durable prefix, so resuming from receipts would be a lie.
    """

    def __init__(self, message: str, diagnostics: Optional[List[Diagnostic]] = None):
        super().__init__(message, exit_code=ExitCode.CANCELLED, diagnostics=diagnostics)
