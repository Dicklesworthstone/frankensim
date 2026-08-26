"""
FrankenSim Python SDK.

A certified, zero-FFI Python interface for FrankenSim computational geometry,
physics simulation, and optimization.
"""

from .client import FrankenSimClient
from .exceptions import (
    FrankenSimError,
    FrankenSimTimeoutError,
    UsageError,
    InputError,
    RefusalError,
    UnavailableError,
    BudgetExceededError,
    CancellationError,
)
from .models import (
    CompareResult,
    Diagnostic,
    EngineeringReport,
    EvidenceColor,
    EvidenceColorKind,
    ExitCode,
    ImportResult,
    PackageAuditResult,
    QoiDiffItem,
    QoiItem,
    RunOutcome,
    SolveOutcome,
    SolveStageReceipt,
    ValidationFinding,
    ValidationResult,
)

__all__ = [
    "FrankenSimClient",
    "FrankenSimError",
    "UsageError",
    "InputError",
    "RefusalError",
    "UnavailableError",
    "BudgetExceededError",
    "CancellationError",
    "CompareResult",
    "Diagnostic",
    "EngineeringReport",
    "EvidenceColor",
    "EvidenceColorKind",
    "ExitCode",
    "ImportResult",
    "PackageAuditResult",
    "QoiDiffItem",
    "QoiItem",
    "RunOutcome",
    "SolveOutcome",
    "SolveStageReceipt",
    "ValidationFinding",
    "ValidationResult",
]

__version__ = "0.0.1"
