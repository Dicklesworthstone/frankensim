"""
FrankenSim Python SDK - Client Implementation.

Provides an intuitive, strongly-typed Python API over the frankensim binary.
Operates strictly out-of-process via stdin/stdout JSON lines in accordance
with ADPT-2026-07.
"""

from __future__ import annotations
import json
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any, Dict, List, Optional, Tuple, Union

from .exceptions import (
    FrankenSimError,
    UsageError,
    InputError,
    RefusalError,
    UnavailableError,
    BudgetExceededError,
    FrankenSimTimeoutError,
    CancellationError,
)
from .models import (
    Diagnostic,
    EngineeringReport,
    EvidenceColor,
    ExitCode,
    ImportResult,
    PackageAuditResult,
    QoiItem,
    RunOutcome,
    SolveOutcome,
    SolveStageReceipt,
    ValidationFinding,
    ValidationResult,
)


class FrankenSimClient:
    """Out-of-process client for the FrankenSim simulation and optimization engine."""

    def __init__(
        self,
        binary_path: Optional[Union[str, Path]] = None,
        default_timeout_s: float = 120.0,
    ):
        self.binary_path = self._resolve_binary(binary_path)
        self.default_timeout_s = default_timeout_s

    @staticmethod
    def _resolve_binary(explicit: Optional[Union[str, Path]]) -> str:
        if explicit:
            p = Path(explicit)
            if p.is_file() and os.access(p, os.X_OK):
                return str(p.resolve())
            raise InputError(f"Specified frankensim binary not found or not executable: {explicit}")

        env_bin = os.environ.get("FRANKENSIM_BIN")
        if env_bin:
            p = Path(env_bin)
            if p.is_file() and os.access(p, os.X_OK):
                return str(p.resolve())

        # Check CARGO_TARGET_DIR
        cargo_target = os.environ.get("CARGO_TARGET_DIR")
        if cargo_target:
            p_rel = Path(cargo_target) / "release" / "frankensim"
            if p_rel.is_file() and os.access(p_rel, os.X_OK):
                return str(p_rel.resolve())
            p_dbg = Path(cargo_target) / "debug" / "frankensim"
            if p_dbg.is_file() and os.access(p_dbg, os.X_OK):
                return str(p_dbg.resolve())

        # Check in target debug / release directories relative to repo root
        current = Path(__file__).resolve()
        for parent in current.parents:
            for cand in [
                parent / "target" / "release" / "frankensim",
                parent / "target" / "debug" / "frankensim",
                Path("/Volumes/USB_NVME/cargo-target/release/frankensim"),
                Path("/Volumes/USB_NVME/cargo-target/debug/frankensim"),
            ]:
                if cand.is_file() and os.access(cand, os.X_OK):
                    return str(cand.resolve())

        # Check in system PATH
        path_bin = shutil.which("frankensim")
        if path_bin:
            return path_bin

        # Fallback default executable name for PATH lookup
        return "frankensim"

    def _execute(
        self,
        subcommand_args: List[str],
        timeout_s: Optional[float] = None,
    ) -> Tuple[int, Dict[str, Any], List[Diagnostic], str, str]:
        """Execute frankensim binary with --json and parse output."""
        cmd = [self.binary_path, "--json"] + subcommand_args
        timeout = timeout_s if timeout_s is not None else self.default_timeout_s

        try:
            proc = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            # A client-side wall-clock kill is NOT the engine's honest exit-6
            # budget stop (no durable prefix exists); conflating them would
            # make callers resume from receipts that were never written.
            raise FrankenSimTimeoutError(
                f"Client timeout killed execution after {timeout}s: {' '.join(cmd)}"
            ) from exc
        except FileNotFoundError as exc:
            raise InputError(
                f"Could not find or invoke executable at `{self.binary_path}`. "
                "Ensure frankensim is built or set FRANKENSIM_BIN."
            ) from exc

        stdout_raw = proc.stdout
        stderr_raw = proc.stderr

        # Parse diagnostics from stderr JSON lines
        diagnostics: List[Diagnostic] = []
        for line in stderr_raw.splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                data = json.loads(line)
                if isinstance(data, dict) and data.get("schema") == "frankensim.cli.diagnostic.v1":
                    diagnostics.append(Diagnostic.from_dict(data))
            except json.JSONDecodeError:
                pass

        # Parse final result record from stdout JSON
        result_data: Dict[str, Any] = {}
        try:
            full_data = json.loads(stdout_raw.strip())
            if isinstance(full_data, dict):
                result_data = full_data
        except (json.JSONDecodeError, ValueError):
            for line in reversed(stdout_raw.splitlines()):
                line = line.strip()
                if not line:
                    continue
                try:
                    data = json.loads(line)
                    if isinstance(data, dict) and "schema" in data:
                        result_data = data
                        break
                except json.JSONDecodeError:
                    pass

        return proc.returncode, result_data, diagnostics, stdout_raw, stderr_raw

    def _raise_for_exit_code(
        self,
        exit_code: int,
        diagnostics: List[Diagnostic],
        raw_stderr: str,
    ) -> None:
        """Map non-zero exit codes to typed Python exceptions."""
        msg = "; ".join(d.message for d in diagnostics) if diagnostics else raw_stderr.strip()
        if not msg:
            msg = f"frankensim process exited with code {exit_code}"

        if exit_code == ExitCode.USAGE:
            raise UsageError(msg, diagnostics)
        elif exit_code == ExitCode.INPUT:
            raise InputError(msg, diagnostics)
        elif exit_code == ExitCode.REFUSED:
            raise RefusalError(msg, diagnostics)
        elif exit_code == ExitCode.UNAVAILABLE:
            raise UnavailableError(msg, diagnostics)
        elif exit_code == ExitCode.BUDGET:
            raise BudgetExceededError(msg, diagnostics)
        elif exit_code == ExitCode.CANCELLED:
            raise CancellationError(msg, diagnostics)
        elif exit_code != ExitCode.SUCCESS:
            raise FrankenSimError(msg, exit_code=exit_code, diagnostics=diagnostics)

    def validate(
        self,
        project_path: Union[str, Path],
        strict: bool = True,
        timeout_s: Optional[float] = None,
    ) -> ValidationResult:
        """Validate a .fsim or .json simulation project."""
        p = Path(project_path)
        exit_code, data, diagnostics, _, stderr_raw = self._execute(
            ["validate", str(p)],
            timeout_s=timeout_s,
        )

        findings = [
            ValidationFinding(code=d.code, what=d.message, fix=d.fix)
            for d in diagnostics
        ]

        result = ValidationResult(
            status=data.get("status", "unknown"),
            command="validate",
            project=data.get("project", str(p)),
            project_hash=data.get("project_hash"),
            fsim_version=data.get("fsim_version"),
            finding_count=int(data.get("finding_count", len(findings))),
            findings=findings,
            authority=data.get("authority", "structural-project-admission"),
            no_claim=data.get("no_claim", ""),
            exit_code=exit_code,
            diagnostics=diagnostics,
        )

        if strict and exit_code != ExitCode.SUCCESS:
            self._raise_for_exit_code(exit_code, diagnostics, stderr_raw)

        return result

    def import_mesh(
        self,
        project_path: Union[str, Path],
        source_path: Union[str, Path],
        ledger_path: Union[str, Path],
        unit: str = "m",
        max_hole_edges: Optional[int] = None,
        step_root: Optional[int] = None,
        target_h: Optional[float] = None,
        strict: bool = True,
        timeout_s: Optional[float] = None,
    ) -> ImportResult:
        """Import CAD / STL geometry into durable ledger."""
        args = [
            "import",
            str(project_path),
            str(source_path),
            str(ledger_path),
            "--unit",
            unit,
        ]

        step_given = step_root is not None
        target_given = target_h is not None
        if step_given != target_given:
            raise UsageError(
                "import_mesh: step_root and target_h must be provided together"
            )
        if step_given and max_hole_edges is not None:
            raise UsageError(
                "import_mesh: max_hole_edges is a mesh-only option and cannot "
                "be combined with the STEP import path (step_root/target_h)"
            )

        if step_given:
            args.extend(["--step-root", str(step_root), "--target-h", str(target_h)])
        elif max_hole_edges is not None:
            args.extend(["--max-hole-edges", str(max_hole_edges)])
        else:
            args.extend(["--max-hole-edges", "0"])

        exit_code, data, diagnostics, _, stderr_raw = self._execute(args, timeout_s=timeout_s)

        result = ImportResult(
            status=data.get("status", "unknown"),
            command="import",
            project=data.get("project", str(project_path)),
            project_hash=data.get("project_hash"),
            ledger=data.get("ledger", str(ledger_path)),
            op_id=int(data.get("op_id", 0)),
            summary_artifact=data.get("summary_artifact"),
            artifact_count=int(data.get("artifact_count", 0)),
            assignment_table=data.get("assignment_table", ""),
            authority=data.get("authority", "retained-import-and-assignment-evidence"),
            no_claim=data.get("no_claim", ""),
            exit_code=exit_code,
            diagnostics=diagnostics,
        )

        if strict and exit_code != ExitCode.SUCCESS:
            self._raise_for_exit_code(exit_code, diagnostics, stderr_raw)

        return result

    def solve(
        self,
        project_path: Union[str, Path],
        ledger_path: Union[str, Path],
        materials: Optional[List[Union[str, Path]]] = None,
        interfaces: Optional[List[Union[str, Path]]] = None,
        strict: bool = True,
        timeout_s: Optional[float] = None,
    ) -> SolveOutcome:
        """Execute solve orchestration for a project and durable ledger."""
        args = ["solve", str(project_path), str(ledger_path)]
        if materials:
            for mat in materials:
                args.extend(["--materials", str(mat)])
        if interfaces:
            for iface in interfaces:
                args.extend(["--interfaces", str(iface)])

        exit_code, data, diagnostics, _, stderr_raw = self._execute(args, timeout_s=timeout_s)

        result = SolveOutcome(
            status=data.get("status", "unknown"),
            command="solve",
            subject=data.get("subject", str(project_path)),
            run_id=data.get("run", ""),
            stages_completed=int(data.get("stages_completed", 0)),
            run_receipt=data.get("run_receipt"),
            exit_code=exit_code,
            diagnostics=diagnostics,
        )

        if strict and exit_code != ExitCode.SUCCESS:
            self._raise_for_exit_code(exit_code, diagnostics, stderr_raw)

        return result

    def report(
        self,
        run_id: str,
        ledger_path: Optional[Union[str, Path]] = None,
        strict: bool = True,
        timeout_s: Optional[float] = None,
    ) -> EngineeringReport:
        """Generate deterministic HTML engineering report and JSON twin."""
        args = ["report", run_id]
        if ledger_path is not None:
            args.append(str(ledger_path))

        exit_code, data, diagnostics, _, stderr_raw = self._execute(args, timeout_s=timeout_s)

        html_file = data.get("html_report", "")
        json_file = data.get("json_twin", "")

        result = EngineeringReport(
            run_id=run_id,
            project_name=data.get("project_name", ""),
            content_hash=data.get("content_hash", ""),
            html_path=html_file,
            json_path=json_file,
            exit_code=exit_code,
            diagnostics=diagnostics,
        )

        if strict and exit_code != ExitCode.SUCCESS:
            self._raise_for_exit_code(exit_code, diagnostics, stderr_raw)

        return result

    def package(
        self,
        run_id: str,
        ledger_path: Optional[Union[str, Path]] = None,
        strict: bool = True,
        timeout_s: Optional[float] = None,
    ) -> PackageAuditResult:
        """Generate and independently verify an evidence package archive."""
        args = ["package", run_id]
        if ledger_path is not None:
            args.append(str(ledger_path))

        exit_code, data, diagnostics, _, stderr_raw = self._execute(args, timeout_s=timeout_s)

        result = PackageAuditResult(
            status=data.get("status", "unknown"),
            command="package",
            package_path=data.get("package", ""),
            merkle_root=data.get("merkle_root", ""),
            claim_count=int(data.get("claim_count", 0)),
            checker=data.get("checker", "unspecified"),
            exit_code=exit_code,
            diagnostics=diagnostics,
            authority=data.get("authority", ""),
        )

        if strict and exit_code != ExitCode.SUCCESS:
            self._raise_for_exit_code(exit_code, diagnostics, stderr_raw)

        return result

    def run(
        self,
        project_path: Union[str, Path],
        ledger_path: Union[str, Path],
        materials: Optional[List[Union[str, Path]]] = None,
        interfaces: Optional[List[Union[str, Path]]] = None,
        strict: bool = True,
        timeout_s: Optional[float] = None,
    ) -> RunOutcome:
        """Execute the one-command convenience run workflow."""
        args = ["run", str(project_path), str(ledger_path)]
        if materials:
            for mat in materials:
                args.extend(["--materials", str(mat)])
        if interfaces:
            for iface in interfaces:
                args.extend(["--interfaces", str(iface)])

        exit_code, data, diagnostics, _, stderr_raw = self._execute(args, timeout_s=timeout_s)
        result = RunOutcome(
            status=data.get("status", "unknown"),
            command=data.get("command", "run"),
            subject=data.get("subject", str(project_path)),
            run_id=data.get("run", ""),
            # Empty means the engine minted no artifact (e.g. an unavailable
            # stage under strict=False): never fabricate success-shaped paths.
            report_file=data.get("report_json", ""),
            package_file=data.get("package", ""),
            verdict=data.get("verdict", ""),
            stage=data.get("stage"),
            exit_code=exit_code,
            diagnostics=diagnostics,
        )

        if strict and exit_code != ExitCode.SUCCESS:
            self._raise_for_exit_code(exit_code, diagnostics, stderr_raw)

        return result

    def compare(
        self,
        left_run: str,
        right_run: str,
        ledger_path: Optional[Union[str, Path]] = None,
        strict: bool = True,
        timeout_s: Optional[float] = None,
    ) -> CompareResult:
        """Semantically compare two retained simulation runs or packages."""
        from .models import CompareResult, QoiDiffItem

        args = ["compare", left_run, right_run]
        if ledger_path is not None:
            args.append(str(ledger_path))

        exit_code, data, diagnostics, _, stderr_raw = self._execute(args, timeout_s=timeout_s)

        qoi_diffs = []
        for q in data.get("qoi_diffs", []):
            qoi_diffs.append(
                QoiDiffItem(
                    name=q.get("name", ""),
                    unit_left=q.get("unit_left", ""),
                    unit_right=q.get("unit_right", ""),
                    nominal_left=float(q.get("nominal_left", 0.0)),
                    nominal_right=float(q.get("nominal_right", 0.0)),
                    delta=float(q.get("delta", 0.0)),
                    rel_delta=float(q.get("rel_delta", 0.0)),
                    color_left=q.get("color_left", "unspecified"),
                    color_right=q.get("color_right", "unspecified"),
                    color_evolution=q.get("color_evolution", "same"),
                    classification=q.get("classification", "same"),
                )
            )

        result = CompareResult(
            status=data.get("status", "unknown"),
            command="compare",
            left_run=data.get("left_run", left_run),
            right_run=data.get("right_run", right_run),
            summary=data.get("summary", ""),
            qoi_count=int(data.get("qoi_count", len(qoi_diffs))),
            qoi_diffs=qoi_diffs,
            changed=bool(data.get("changed", False)),
            project_hash_left=data.get("project_hash_left", ""),
            project_hash_right=data.get("project_hash_right", ""),
            same_project=bool(data.get("same_project", True)),
            authority=data.get("authority", "none"),
            no_claim=data.get("no_claim", ""),
            exit_code=exit_code,
            diagnostics=diagnostics,
        )

        if strict and exit_code != ExitCode.SUCCESS:
            self._raise_for_exit_code(exit_code, diagnostics, stderr_raw)

        return result
