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
        cls.reference_fsim = repo_root / "data" / "reference-project" / "cooling-reference.fsim"
        cls.reference_stl = repo_root / "data" / "reference-project" / "plate.stl"

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

            # 2. Solve: the heatsink example declares a conduction section
            # with an airflow-convection law, so every stage executes and the
            # run completes (executed 2026-09-02 on yto with the release
            # binary; the old expectation of a conduction gap is gone).
            solve_res = self.client.solve(
                project_path=self.fsim,
                ledger_path=ledger,
                materials=[self.pack],
                strict=True,
            )
            self.assertEqual(solve_res.exit_code, ExitCode.SUCCESS)
            self.assertEqual(solve_res.status, "completed")
            self.assertEqual(solve_res.stages_completed, 7)

            # 3. Compare reads retained runs only: an unknown run id refuses
            # through the solve loader, and a run compared with itself is
            # a real, empty comparison (one QoI, zero delta, nothing changed).
            with self.assertRaises(RefusalError):
                self.client.compare(
                    solve_res.run_id,
                    "f" * 64,
                    ledger_path=ledger,
                    strict=True,
                )

            comp = self.client.compare(
                solve_res.run_id,
                solve_res.run_id,
                ledger_path=ledger,
                strict=True,
            )
            self.assertEqual(comp.status, "ok")
            self.assertEqual(comp.exit_code, ExitCode.SUCCESS)
            self.assertEqual(comp.authority, "projection-of-retained-receipts")
            self.assertFalse(comp.changed)
            self.assertTrue(comp.same_project)
            self.assertEqual(len(comp.project_hash_left), 64)
            self.assertEqual(comp.project_hash_left, comp.project_hash_right)
            self.assertEqual(comp.qoi_count, 1)
            self.assertEqual(len(comp.qoi_diffs), 1)
            self.assertEqual(comp.qoi_diffs[0].name, "temperature-max")
            self.assertEqual(comp.qoi_diffs[0].delta, 0.0)
            self.assertEqual(comp.qoi_diffs[0].classification, "same")
            self.assertIn("identical runs", comp.summary)

    def test_run_exposes_retained_qoi_vocabulary(self):
        if not (
            self.reference_fsim.exists()
            and self.reference_stl.exists()
            and self.pack.exists()
            and Path(self.client.binary_path).is_file()
        ):
            self.skipTest("reference fixtures or frankensim binary unavailable")

        with tempfile.TemporaryDirectory() as tmpdir:
            ledger = Path(tmpdir) / "python_run_ledger.db"
            imported = self.client.import_mesh(
                project_path=self.reference_fsim,
                source_path=self.reference_stl,
                ledger_path=ledger,
                unit="m",
                max_hole_edges=0,
            )
            self.assertEqual(imported.exit_code, ExitCode.SUCCESS)

            outcome = self.client.run(
                project_path=self.reference_fsim,
                ledger_path=ledger,
                materials=[self.pack],
            )
            self.assertEqual(outcome.exit_code, ExitCode.SUCCESS)
            self.assertEqual(outcome.status, "completed")
            self.assertIn(outcome.verdict, {"pass", "fail", "indeterminate"})

            report_path = Path(outcome.report_file)
            self.assertTrue(report_path.is_file(), outcome)
            report = json.loads(report_path.read_text())
            self.assertEqual(report["schema"], "frankensim.report.engineering.v1")

            qois = report["qois"]
            self.assertEqual(len(qois), 1)
            self.assertEqual(qois[0]["name"], "temperature-max")
            self.assertIn(qois[0]["color"], {"Verified", "Validated", "Estimated"})

            requirements = report["requirements"]
            self.assertEqual(len(requirements), 1)
            self.assertEqual(requirements[0]["qoi"], "temperature-max")
            self.assertEqual(requirements[0]["outcome"], outcome.verdict)

            budget_terms = report["budget_terms"]
            self.assertEqual(len(budget_terms), 8)
            for term in budget_terms:
                self.assertEqual(term["qoi"], "temperature-max")
                self.assertIn(term["state"], {"measured", "no-data"})
                if term["state"] == "no-data":
                    self.assertIsNone(term["value"])
                else:
                    self.assertIsInstance(term["value"], (int, float))

            # Independent re-derivation (bead q61wp.14 item 1): `package`
            # re-runs the solver-free checker on the retained bytes, and the
            # QoI claim's colour read from the package itself must equal the
            # colour the report shows. A report that coloured a claim more
            # strongly than its package would fail here.
            audit = self.client.package(outcome.run_id, ledger_path=ledger, strict=True)
            self.assertEqual(audit.exit_code, ExitCode.SUCCESS)
            self.assertTrue(audit.is_verified, audit)
            self.assertEqual(audit.authority, "structural-integrity-only")
            self.assertEqual(len(audit.merkle_root), 64)
            self.assertGreaterEqual(audit.claim_count, 2)
            package_path = Path(audit.package_path)
            self.assertTrue(package_path.is_file(), audit)
            package = json.loads(package_path.read_text())
            prefix = f"qoi.{qois[0]['name']}."
            claims = [c for c in package["claims"] if c["id"].startswith(prefix)]
            self.assertEqual(len(claims), 1, [c["id"] for c in package["claims"]])
            package_color = claims[0]["color"]["kind"].capitalize()
            self.assertEqual(package_color, qois[0]["color"], claims[0])

    def test_run_without_materials_surfaces_material_resolve_refusal(self):
        if not (
            self.reference_fsim.exists()
            and self.reference_stl.exists()
            and Path(self.client.binary_path).is_file()
        ):
            self.skipTest("reference fixtures or frankensim binary unavailable")

        with tempfile.TemporaryDirectory() as tmpdir:
            ledger = Path(tmpdir) / "python_no_materials_ledger.db"
            imported = self.client.import_mesh(
                project_path=self.reference_fsim,
                source_path=self.reference_stl,
                ledger_path=ledger,
                unit="m",
                max_hole_edges=0,
            )
            self.assertEqual(imported.exit_code, ExitCode.SUCCESS)

            # The public client emits no --materials flags for None; the real
            # CLI must therefore reach its material-resolve refusal rather
            # than inventing a default card or a success-shaped export.
            outcome = self.client.run(
                project_path=self.reference_fsim,
                ledger_path=ledger,
                materials=None,
                strict=False,
            )
            self.assertEqual(outcome.exit_code, ExitCode.REFUSED)
            self.assertEqual(outcome.status, "refused")
            self.assertEqual(outcome.command, "solve")
            self.assertEqual(outcome.stage, "material-resolve")
            self.assertEqual(outcome.report_file, "")
            self.assertEqual(outcome.package_file, "")
            self.assertTrue(
                any(
                    diagnostic.code == "project-material-card-unknown"
                    for diagnostic in outcome.diagnostics
                ),
                outcome.diagnostics,
            )


if __name__ == "__main__":
    unittest.main()
