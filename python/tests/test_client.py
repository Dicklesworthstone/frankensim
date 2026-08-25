"""
Unit and integration tests for the FrankenSim Python SDK.
"""

import json
from pathlib import Path
import tempfile
import unittest

from frankensim import (
    FrankenSimClient,
    EvidenceColor,
    EvidenceColorKind,
    ExitCode,
    UsageError,
    InputError,
    RefusalError,
    UnavailableError,
    BudgetExceededError,
)


class TestFrankenSimModels(unittest.TestCase):
    def test_evidence_color_parsing(self):
        v = EvidenceColor.from_dict({"verified": {"lo": 300.0, "hi": 350.0}})
        self.assertEqual(v.kind, EvidenceColorKind.VERIFIED)
        self.assertEqual(v.bounds, (300.0, 350.0))

        val = EvidenceColor.from_dict({"validated": {"regime": "turbulent"}})
        self.assertEqual(val.kind, EvidenceColorKind.VALIDATED)
        self.assertEqual(val.regime, "turbulent")

        est = EvidenceColor.from_dict({"estimated": {"estimator": "gci", "dispersion": 0.05}})
        self.assertEqual(est.kind, EvidenceColorKind.ESTIMATED)
        self.assertEqual(est.dispersion, 0.05)

        w = EvidenceColor.from_dict({"waived": {}})
        self.assertEqual(w.kind, EvidenceColorKind.WAIVED)


class TestFrankenSimExceptions(unittest.TestCase):
    def test_exception_exit_codes(self):
        self.assertEqual(UsageError("bad args").exit_code, ExitCode.USAGE)
        self.assertEqual(InputError("missing file").exit_code, ExitCode.INPUT)
        self.assertEqual(RefusalError("invalid physics").exit_code, ExitCode.REFUSED)
        self.assertEqual(UnavailableError("bead open").exit_code, ExitCode.UNAVAILABLE)
        self.assertEqual(BudgetExceededError("timeout").exit_code, ExitCode.BUDGET)


class TestFrankenSimClientIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Locate project files relative to repository root
        repo_root = Path(__file__).resolve().parents[2]
        cls.repo_root = repo_root
        cls.fsim = repo_root / "examples" / "heatsink-fan" / "heatsink-fan.fsim"
        cls.stl = repo_root / "examples" / "heatsink-fan" / "heatsink.stl"
        cls.pack = repo_root / "data" / "reference-project" / "aa6061.fsmcdpk"

        # Build frankensim binary if needed or locate
        cls.client = FrankenSimClient()

    def test_validate_success(self):
        if not self.fsim.exists():
            self.skipTest("heatsink-fan.fsim not found")
        res = self.client.validate(self.fsim, strict=False)
        self.assertEqual(res.exit_code, ExitCode.SUCCESS)
        self.assertTrue(res.is_valid)
        self.assertEqual(res.finding_count, 0)

    def test_validate_nonexistent_fails_with_input_error(self):
        with self.assertRaises(InputError):
            self.client.validate("/nonexistent/fake_project.fsim", strict=True)

    def test_import_and_solve_workflow(self):
        if not (self.fsim.exists() and self.stl.exists() and self.pack.exists()):
            self.skipTest("Test fixtures missing")

        with tempfile.TemporaryDirectory() as tmpdir:
            ledger = Path(tmpdir) / "python_test_ledger.db"

            # 1. Import
            imp = self.client.import_mesh(
                project_path=self.fsim,
                source_path=self.stl,
                ledger_path=ledger,
                unit="m",
                max_hole_edges=0,
                strict=True,
            )
            self.assertEqual(imp.exit_code, ExitCode.SUCCESS)
            self.assertGreater(imp.artifact_count, 0)
            self.assertEqual(imp.op_id, 1)

            # 2. Solve (stops at conduction gap with exit code UNAVAILABLE = 5)
            solve_res = self.client.solve(
                project_path=self.fsim,
                ledger_path=ledger,
                materials=[self.pack],
                strict=False,
            )
            self.assertEqual(solve_res.exit_code, ExitCode.UNAVAILABLE)
            self.assertTrue(any("conduction" in d.message or "conduction" in d.fix for d in solve_res.diagnostics))

            # 3. Compare runs
            comp = self.client.compare("baseline_run", "candidate_run", ledger_path=ledger, strict=True)
            self.assertEqual(comp.status, "changed")
            self.assertEqual(comp.exit_code, ExitCode.SUCCESS)
            self.assertGreaterEqual(comp.qoi_count, 2)
            self.assertTrue(any(q.name == "junction_maximum" for q in comp.qoi_diffs))


if __name__ == "__main__":
    unittest.main()
