"""
FrankenSim Python SDK - Strongly Typed Models and Data Structures.

Under ADPT-2026-07 (decision record fs-govern::adapter_policy), Python surfaces
are zero-FFI out-of-process client wrappers over the certified frankensim binary.
All scientific claims, error budgets, evidence colors, and units are preserved
bit-accurately without loss or laundering.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple


class ExitCode(int, Enum):
    """Stable exit classes from the FrankenSim CLI contract."""
    SUCCESS = 0
    USAGE = 2
    INPUT = 3
    REFUSED = 4
    UNAVAILABLE = 5
    BUDGET = 6
    CANCELLED = 130


class EvidenceColorKind(str, Enum):
    """Epistemic evidence classifications."""
    VERIFIED = "verified"
    VALIDATED = "validated"
    ESTIMATED = "estimated"
    WAIVED = "waived"


@dataclass(frozen=True)
class EvidenceColor:
    """Rigorous evidence classification with quantitative bounds."""
    kind: EvidenceColorKind
    bounds: Optional[Tuple[float, float]] = None
    regime: Optional[str] = None
    estimator: Optional[str] = None
    dispersion: Optional[float] = None

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> EvidenceColor:
        if "verified" in data or data.get("kind") == "verified":
            lo = float(data.get("lo", data.get("verified", {}).get("lo", 0.0)))
            hi = float(data.get("hi", data.get("verified", {}).get("hi", 0.0)))
            return cls(kind=EvidenceColorKind.VERIFIED, bounds=(lo, hi))
        elif "validated" in data or data.get("kind") == "validated":
            regime = str(data.get("regime", data.get("validated", {}).get("regime", "")))
            return cls(kind=EvidenceColorKind.VALIDATED, regime=regime)
        elif "waived" in data or data.get("kind") == "waived":
            return cls(kind=EvidenceColorKind.WAIVED)
        else:
            estimator = str(data.get("estimator", data.get("estimated", {}).get("estimator", "unspecified")))
            dispersion = float(data.get("dispersion", data.get("estimated", {}).get("dispersion", 0.0)))
            return cls(kind=EvidenceColorKind.ESTIMATED, estimator=estimator, dispersion=dispersion)


@dataclass(frozen=True)
class Diagnostic:
    """Structured diagnostic record emitted by FrankenSim."""
    schema: str
    command: str
    severity: str
    code: str
    message: str
    fix: str
    subject: Optional[str] = None
    hint: Optional[str] = None

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> Diagnostic:
        return cls(
            schema=data.get("schema", "frankensim.cli.diagnostic.v1"),
            command=data.get("command", ""),
            severity=data.get("severity", "error"),
            code=data.get("code", ""),
            message=data.get("message", ""),
            fix=data.get("fix", ""),
            subject=data.get("subject"),
            hint=data.get("hint"),
        )


@dataclass(frozen=True)
class ValidationFinding:
    """Finding from project structural validation."""
    code: str
    what: str
    fix: str


@dataclass(frozen=True)
class ValidationResult:
    """Outcome of project schema and semantic admission."""
    status: str
    command: str
    project: str
    project_hash: Optional[str]
    fsim_version: Optional[str]
    finding_count: int
    findings: List[ValidationFinding] = field(default_factory=list)
    authority: str = "structural-project-admission"
    no_claim: str = ""
    exit_code: int = 0
    diagnostics: List[Diagnostic] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        return self.status == "ok" and self.finding_count == 0


@dataclass(frozen=True)
class ImportResult:
    """Outcome of geometric mesh/STEP import."""
    status: str
    command: str
    project: str
    project_hash: Optional[str]
    ledger: str
    op_id: int
    summary_artifact: Optional[str]
    artifact_count: int
    assignment_table: str
    authority: str
    no_claim: str
    exit_code: int = 0
    diagnostics: List[Diagnostic] = field(default_factory=list)


@dataclass(frozen=True)
class SolveStageReceipt:
    """Durable receipt for a completed solve stage."""
    stage: str
    wall_s: float
    receipt: str


@dataclass(frozen=True)
class SolveOutcome:
    """Result of solving a simulation project."""
    status: str
    command: str
    subject: str
    run_id: str
    stages_completed: int
    run_receipt: Optional[str] = None
    stages: List[SolveStageReceipt] = field(default_factory=list)
    exit_code: int = 0
    diagnostics: List[Diagnostic] = field(default_factory=list)

    @property
    def is_completed(self) -> bool:
        return self.status == "completed" and self.exit_code == 0


@dataclass(frozen=True)
class QoiItem:
    """Quantity of Interest report record."""
    name: str
    description: str
    nominal_value: float
    unit: str
    color: EvidenceColor
    discretization_error: float
    parameter_uncertainty: float
    surrogate_error: float
    total_uncertainty_budget: float
    source_root: str


@dataclass(frozen=True)
class EngineeringReport:
    """Complete HTML and JSON twin engineering report."""
    run_id: str
    project_name: str
    content_hash: str
    html_path: str
    json_path: str
    qois: List[QoiItem] = field(default_factory=list)
    exit_code: int = 0
    diagnostics: List[Diagnostic] = field(default_factory=list)


@dataclass(frozen=True)
class PackageAuditResult:
    """Outcome of offline evidence package generation and audit."""
    status: str
    command: str
    package_path: str
    content_hash: str
    receipt_count: int
    audit_verdict: str
    exit_code: int = 0
    diagnostics: List[Diagnostic] = field(default_factory=list)

    @property
    def is_verified(self) -> bool:
        return self.status == "ok" and self.audit_verdict == "admissible"


@dataclass(frozen=True)
class RunOutcome:
    """Unified result of one-command `run` workflow."""
    status: str
    command: str
    subject: str
    run_id: str
    report_file: str
    package_file: str
    exit_code: int = 0
    diagnostics: List[Diagnostic] = field(default_factory=list)
    verdict: str = ""
    stage: Optional[str] = None


@dataclass(frozen=True)
class QoiDiffItem:
    """Semantic comparison record for a single QoI."""
    name: str
    unit_left: str
    unit_right: str
    nominal_left: float
    nominal_right: float
    delta: float
    rel_delta: float
    color_left: str
    color_right: str
    color_evolution: str
    classification: str


@dataclass(frozen=True)
class CompareResult:
    """Outcome of semantic run comparison."""
    status: str
    command: str
    left_run: str
    right_run: str
    summary: str
    qoi_count: int
    qoi_diffs: List[QoiDiffItem] = field(default_factory=list)
    #: True when any retained receipt row differs between the two runs.
    changed: bool = False
    #: Canonical project hashes of the two runs; a design change is a
    #: different project hash, so they may differ.
    project_hash_left: str = ""
    project_hash_right: str = ""
    same_project: bool = True
    authority: str = "projection-of-retained-receipts"
    no_claim: str = ""
    exit_code: int = 0
    diagnostics: List[Diagnostic] = field(default_factory=list)
