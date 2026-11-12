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


class SpecimenInputsWithoutClarinetNames(unittest.TestCase):
    """Manufactured rows: diameters printed at holes, other hole quantities, no mouthpiece."""

    DIMS = {
        ("bore", "station-000"): {"axial-position": 0.0, "diameter": 0.002},
        ("bore", "station-500"): {"axial-position": 0.5, "diameter": 0.016},
        ("bore", "at-hole-01"): {"diameter": 0.006},
        ("bore", "at-hole-02"): {"diameter": 0.006},
        ("tone-holes", "hole-01"): {"axial-position": 0.2, "height": 0.003, "diameter": 0.004},
        ("tone-holes", "hole-02"): {"axial-position": 0.2, "height": 0.003, "diameter": 0.005},
        ("tone-holes", "hole-03"): {"axial-position": 0.3, "height": 0.004, "diameter": 0.006},
    }

    def test_diameter_printed_at_a_hole_takes_that_holes_position_once(self):
        self.assertEqual(lower.bore_knots(self.DIMS), [(0.0, 0.001), (0.2, 0.003), (0.5, 0.008)])

    def test_hole_quantities_are_named_and_holes_are_omitted_by_key(self):
        desc = {"tonehole.height_quantity": "height", "tonehole.diameter_quantity": "diameter", "tonehole.omit": "1,2"}
        self.assertEqual(lower.holes_of(desc, self.DIMS), {3: (0.3, 0.004, 0.003)})

    def test_an_omitted_hole_stays_in_when_the_fingering_opens_it(self):
        desc = {"tonehole.height_quantity": "height", "tonehole.diameter_quantity": "diameter", "tonehole.omit": "1,2"}
        self.assertEqual(sorted(lower.holes_of(desc, self.DIMS, {2})), [2, 3])

    def test_an_opened_hole_with_no_printed_height_takes_the_declared_one_or_refuses(self):
        dims = {("tone-holes", "hole-01"): {"axial-position": 0.1, "diameter": 0.0003}}
        desc = {"tonehole.height_quantity": "height", "tonehole.diameter_quantity": "diameter", "tonehole.omit": "1"}
        self.assertEqual(lower.holes_of(desc, dims), {})
        with self.assertRaisesRegex(SystemExit, "hole-01 is missing height"):
            lower.holes_of(desc, dims, {1})
        desc["tonehole.height_unprinted_m"] = "0.006"
        self.assertEqual(lower.holes_of(desc, dims, {1}), {1: (0.1, 0.006, 0.00015)})

    def test_without_the_keys_the_clarinet_quantities_are_required(self):
        with self.assertRaisesRegex(SystemExit, "hole-0[123] is missing chimney-height"):
            lower.holes_of({}, self.DIMS)
        dims = {("tone-holes", "hole-01"): {"axial-position": 0.2, "chimney-height": 0.003, "min-diameter": 0.004}}
        self.assertEqual(lower.holes_of({}, dims), {1: (0.2, 0.003, 0.002)})

    def test_declared_absent_mouthpiece_adds_no_segments(self):
        self.assertEqual(lower.mouthpiece_segments({"mouthpiece.part": "-"}, self.DIMS), [])
        with self.assertRaisesRegex(SystemExit, "barrel_pull_m needs a mouthpiece"):
            lower.mouthpiece_segments({"mouthpiece.part": "-", "mouthpiece.barrel_pull_m": "0.001"}, self.DIMS)


class OptionalCardRows(unittest.TestCase):
    """The optional rows, the per-note replacements and the extra fingerings, on manufactured values."""

    def test_rows_appear_only_with_their_keys_and_a_per_note_dash_removes_one(self):
        desc = description()
        plain = lower.performance(desc, SEGMENTS, 1800, "x")
        for prefix in ("slit_length", "gas_field", "open_hole_interaction", "tonehole_lengths", "constriction_noise", "open_hole_orifice"):
            self.assertNotIn(f"\n{prefix} ", plain)
        desc.update({
            "reed.lay_stiffness": "5e6", "reed.lay_exponent": "2.0", "reed.lay_onset_m": "1e-4", "reed.lay_loss_s_per_m": "0.0",
            "reed.channel_length_m": "0.0014", "tonehole.interaction_law": "tmmi", "tonehole.interaction_epsilon": "2",
            "tonehole.lengths": "lefebvre-scavone-2012", "reed.constriction_noise_seed": "7",
            "reed.constriction_noise_duct_radius_m": "0.0069", "tonehole.orifice_holes": "open",
            "at.y.tonehole.lengths": "-", "at.y.reed.channel_length_m": "0.025",
        })
        full = lower.performance(desc, SEGMENTS, 1800, "x")
        self.assertTrue(line_starting(full, "reed ").endswith(" 5e6 2.0 1e-4 0.0"))
        self.assertIn(f"\n{line_starting(full, 'reed ')}\nslit_length 0.0014\nbaffled_outlet ", full)
        tail = full.split("cylinder 0.0073 0.37\n", 1)[1].split("compile_limits", 1)[0]
        self.assertEqual(tail, "open_hole_interaction tmmi 2\ntonehole_lengths lefebvre-scavone-2012\n"
                               "constriction_noise 7 0.0069\nopen_hole_orifice 1\n")
        other = lower.performance(desc, SEGMENTS, 1800, "y")
        self.assertNotIn("\ntonehole_lengths", other)
        self.assertIn("\nslit_length 0.025\n", other)

    def test_an_orifice_row_needs_an_open_hole_and_a_lay_needs_all_four_fields(self):
        desc = description()
        desc["tonehole.orifice_holes"] = "open"
        with self.assertRaisesRegex(SystemExit, "opens no hole"):
            lower.performance(desc, ["cylinder 0.0073 0.2", "hole 0.002 0.004 0.0073 0", "cylinder 0.0073 0.37"], 1800, "x")
        desc = description()
        desc["reed.lay_stiffness"] = "5e6"
        with self.assertRaisesRegex(SystemExit, "missing reed.lay_exponent"):
            lower.reed_record(desc, "x")

    def test_gas_field_counts_its_stations_and_refuses_a_short_station(self):
        desc = description()
        desc.update({
            "bore.gas_field_mode": "continuous", "bore.gas_field_steps": "8", "bore.gas_field_basis": "mole",
            "bore.gas_field_species": "N2 CO2", "bore.gas_station.01": "0 300 1 0", "bore.gas_station.02": "0.57 310 0.5 0.5",
        })
        self.assertIn("\ngas_field continuous 8 2 mole 2 N2 CO2\ngas_station 0 300 1 0\ngas_station 0.57 310 0.5 0.5\n",
                      lower.performance(desc, SEGMENTS, 1800, "x"))
        desc["bore.gas_station.02"] = "0.57 310 0.5"
        with self.assertRaisesRegex(SystemExit, "one fraction per species"):
            lower.performance(desc, SEGMENTS, 1800, "x")

    def test_per_note_keys_must_name_a_declared_fingering_and_a_per_note_key(self):
        chart = [{"fingering": "e3", "written": "E3"}, {"fingering": "a3", "written": "A3"}]
        rows = lower.extra_fingerings({"fingering.gs3": "G#3 F#3 first 1,3,4"}, chart)
        self.assertEqual([r["fingering"] for r in rows], ["e3", "gs3", "a3"])
        self.assertEqual(rows[1]["open_holes"], "1,3,4")
        with self.assertRaisesRegex(SystemExit, "repeats a fingering"):
            lower.extra_fingerings({"fingering.e3": "E3 D3 first -"}, chart)
        with self.assertRaisesRegex(SystemExit, "names fingering zz"):
            lower.check_per_note_keys({"at.zz.tonehole.lengths": "-"}, {"e3"})
        with self.assertRaisesRegex(SystemExit, "is not one of the per-note keys"):
            lower.check_per_note_keys({"at.e3.observer.law": "none"}, {"e3"})


class ClarinetCardsAreGenerated(unittest.TestCase):
    """The committed clarinet cards and their index are what the description lowers to, byte for byte."""

    def test_the_description_lowers_to_every_committed_card(self):
        import os
        import subprocess
        import sys
        import tempfile
        root = Path(__file__).resolve().parents[3]
        committed = root / "data/instruments/descriptions/clarinet/performances"
        with tempfile.TemporaryDirectory() as out:
            subprocess.run(
                [sys.executable, "-I", str(Path(__file__).with_name("lower_reed_performance.py")),
                 "data/instruments/descriptions/clarinet/description.tsv", out],
                cwd=root, check=True, capture_output=True,
            )
            written = sorted(os.listdir(out))
            self.assertEqual(written, sorted(p.name for p in committed.iterdir()))
            for name in written:
                self.assertEqual((Path(out) / name).read_bytes(), (committed / name).read_bytes(), name)


if __name__ == "__main__":
    unittest.main()
