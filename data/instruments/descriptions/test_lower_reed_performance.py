#!/usr/bin/env python3
"""G0/G3 regressions for physical inputs at the existing wind-description boundary.

Run from the repository root with:
    python3 -B data/instruments/descriptions/test_lower_reed_performance.py
These tests exercise declaration lowering, not acoustic or specimen validation.
"""
import importlib.util
import math
from pathlib import Path
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location(
    "lower_reed_performance", Path(__file__).with_name("lower_reed_performance.py")
)
lower = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lower)


def description():
    """Explicit manufactured description; no source material defaults."""
    return {
        "audio.rate_hz": "48000",
        "audio.duration_s": "0.1",
        "audio.full_scale_pa": "20",
        "reed.rest_opening_m": "0.0004",
        "reed.width_m": "0.013",
        "reed.stiffness_per_area_pa_m": "8660000",
        "reed.effective_area_m2": "0.0000761",
        "reed.mass_per_area_kg_m2": "0.05",
        "reed.damping_rate_1_s": "3000",
        "reed.operating_blowing_pa": "1800",
        "reed.tip_opening_cap_m": "0.00115",
        "ambient.temperature_k": "300.15",
        "ambient.pressure_pa": "101325",
        "ambient.relative_humidity": "0.7",
        "observer.record": "baffled_outlet 0 0 0.5 4 16 1000",
        "observer.law": "stokes-kirchhoff",
        "observer.termination": "flanged",
        "compile_limits": "100000 1024",
        "gesture.control_rate_hz": "1000",
        "gesture.attack_start_s": "0",
        "gesture.attack_s": "0.01",
        "gesture.release_start_s": "0.08",
        "gesture.release_s": "0.01",
    }


SEGMENTS = [
    "cylinder 0.0073 0.2",
    "hole 0.002 0.004 0.0073 1",
    "cylinder 0.0073 0.37",
]


def line_starting(source, prefix):
    return next(line for line in source.splitlines() if line.startswith(prefix))


class ReedPhysicsInputs(unittest.TestCase):
    def test_pressure_sweep_changes_drive_without_rebuilding_the_reed(self):
        desc = description()
        baseline = lower.performance(desc, SEGMENTS, 1800)
        expected = line_starting(baseline, "reed ")
        for pressure in (0, 600, 1800, 3600, 7200):
            with self.subTest(pressure_pa=pressure):
                source = lower.performance(desc, SEGMENTS, pressure)
                self.assertEqual(line_starting(source, "reed "), expected)
                self.assertEqual(source.split("schedule\n")[0], baseline.split("schedule\n")[0])
                self.assertIn(f"pressure-pa\t{lower.exp(pressure)}\n", source)
        opening, _, closing, mass, stiffness, damping = map(float, expected.split()[1:])
        self.assertEqual(opening, 0.0004)
        self.assertEqual(closing, 8660000 * opening)
        self.assertAlmostEqual(mass, 0.05 * 0.0000761, delta=1e-20)
        self.assertAlmostEqual(stiffness, 8660000 * 0.0000761, delta=1e-12)
        self.assertAlmostEqual(damping, 3000 / (2 * math.sqrt(stiffness / mass)), delta=1e-15)

    def test_explicit_opening_change_changes_mechanical_closing_pressure(self):
        desc = description()
        before = list(map(float, lower.reed_record(desc).split()))
        desc["reed.rest_opening_m"] = "0.0008"
        after = list(map(float, lower.reed_record(desc).split()))
        self.assertEqual(after[0], 2 * before[0])
        self.assertEqual(after[2], 2 * before[2])
        self.assertEqual(after[3:], before[3:])
        desc["reed.rest_opening_m"] = "0.002"
        with self.assertRaisesRegex(SystemExit, "exceeds tip opening cap"):
            lower.reed_record(desc)

    def test_literal_mechanical_parameters_stay_literal(self):
        desc = description()
        del desc["reed.stiffness_per_area_pa_m"]
        desc.update({
            "reed.closing_pressure_pa": "6000", "reed.mass_kg": "0.00001",
            "reed.stiffness_n_m": "100", "reed.damping_ratio": "0.35",
        })
        self.assertEqual(lower.reed_record(desc), "0.0004 0.013 6000 0.00001 100 0.35")


class IndependentExteriorGas(unittest.TestCase):
    def test_mixture_exterior_reaches_open_hole_performance_without_changing_interior(self):
        desc = description()
        desc.update({
            "exterior.temperature_k": "293.15", "exterior.pressure_pa": "101325",
            "exterior.composition": "mole 3 N2 0.7 O2 0.2 He 0.1",
        })
        source = lower.performance(desc, SEGMENTS, 1800)
        self.assertIn("ambient 300.15 101325 0.7\n", source)
        self.assertIn("exterior mixture 293.15 101325 mole 3 N2 0.7 O2 0.2 He 0.1\n", source)
        self.assertIn("segments 3\n" + "\n".join(SEGMENTS) + "\n", source)
        self.assertLess(source.index("baffled_outlet "), source.index("exterior "))
        self.assertLess(source.index("exterior "), source.index("observer_law "))
        self.assertLess(source.index("observer_law "), source.index("termination "))

    def test_all_gas_forms_use_the_same_declaration_for_either_region(self):
        cases = [
            ({"temperature_k": "293.15", "pressure_pa": "101325", "relative_humidity": "0.5"},
             "293.15 101325 0.5"),
            ({"temperature_k": "303.15", "pressure_pa": "90000", "composition": "mole 2 N2 0.8 CO2 0.2"},
             "mixture 303.15 90000 mole 2 N2 0.8 CO2 0.2"),
            ({"gas_pack": "declared-gas"}, "pack declared-gas"),
        ]
        for fields, expected in cases:
            for prefix in ("ambient", "exterior"):
                with self.subTest(prefix=prefix, fields=fields):
                    desc = {f"{prefix}.{key}": value for key, value in fields.items()}
                    with patch.object(lower.os.path, "isfile", return_value=True):
                        self.assertEqual(lower.gas_record(desc, prefix), expected)

    def test_pack_authority_is_preserved_and_missing_pack_refuses(self):
        desc = description()
        desc.update({"exterior.gas_pack": "declared-gas", "exterior.composition": "UNKNOWN"})
        with patch.object(lower.os.path, "isfile", return_value=True) as exists:
            source = lower.performance(desc, SEGMENTS, 1800)
        self.assertIn("exterior pack declared-gas\n", source)
        exists.assert_called_once_with("data/matdb/seed-v1/declared-gas/manifest.tsv")
        with patch.object(lower.os.path, "isfile", return_value=False):
            with self.assertRaisesRegex(SystemExit, "missing gas pack manifest"):
                lower.performance(desc, SEGMENTS, 1800)

    def test_absent_exterior_is_omitted_and_partial_exterior_never_defaults(self):
        desc = description()
        self.assertNotIn("\nexterior ", lower.performance(desc, SEGMENTS, 1800))
        desc["exterior.relative_humidity"] = "0.5"
        with self.assertRaisesRegex(SystemExit, "missing exterior.temperature_k"):
            lower.performance(desc, SEGMENTS, 1800)
        desc["exterior.temperature_k"] = "UNKNOWN"
        with self.assertRaisesRegex(SystemExit, "exterior.temperature_k is UNKNOWN"):
            lower.performance(desc, SEGMENTS, 1800)


if __name__ == "__main__":
    unittest.main()
