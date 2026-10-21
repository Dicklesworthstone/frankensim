"""The cello bow-kinematics derivation, checked on a rig whose answer is known by construction.

Run: python3 -m unittest data/instruments/descriptions/players/test_derive_cello_bow_kinematics.py -v   (from the repository root)
"""
import csv
import math
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import derive_cello_bow_kinematics as D  # noqa: E402


def rig(n: int = 50, fs: float = 240.0, v: float = 0.25, beta: float = 0.1, tilt_deg: float = 0.0, lift: float = 0.0):
    """Nut at the origin, string along +y 0.7 m long, bow axis along +x (square to the string) with the frog 0.1 m before the string, so
    the closest approach is 0.1 m along the bow; the frog then moves away from the string at v along the bow's own axis, so the
    contact point slides toward the tip at v. Lifted `lift` m above the string plane. The bow crosses the string at y = (1 - beta) * 0.7."""
    rows = []
    ang = math.radians(tilt_deg)
    for i in range(n):
        t = i / fs
        shift = -v * t  # the frog moves by -v t along the bow axis, the contact point moves by +v t along it
        ux, uy = math.cos(ang), math.sin(ang)
        frog = (-0.1 + shift * ux, (1 - beta) * 0.7 + shift * uy, lift)
        tip = (frog[0] + 0.69 * ux, frog[1] + 0.69 * uy, lift)
        row = [t] + list(tip) + list(frog) + [0.0] * 9 + [0.0, 0.0, 0.0] + [0.0, 0.7, 0.0]
        rows.append(row)
    return rows


class Geometry(unittest.TestCase):
    def test_a_bow_square_to_the_string_reads_position_velocity_distance_and_skew_exactly(self):
        out = D.derive(rig(v=0.25, beta=0.1, lift=0.02))
        for fr, s, vel, beta, sep, skew in out:
            self.assertAlmostEqual(vel, 0.25, places=9)
            self.assertAlmostEqual(beta, 0.1, places=9)
            self.assertAlmostEqual(sep, 0.02, places=9)
            self.assertAlmostEqual(skew, 90.0, places=7)
        self.assertAlmostEqual(out[0][1], 0.1, places=9)  # the crossing is 0.1 m from the frog at t = 0
        self.assertAlmostEqual(out[10][1], 0.1 + 0.25 * 10 / 240.0, places=9)

    def test_the_velocity_is_relative_to_the_string_line_so_the_instruments_motion_cancels(self):
        rows = rig(v=0.0, beta=0.1)
        for i, r in enumerate(rows):  # slide the whole instrument (nut and bridge) along the bow direction; the bow stays put
            drift = 0.3 * (i / 240.0)
            for k in (16, 19):  # MNut x, MBridge x
                r[k] += drift
        out = D.derive(rows)
        self.assertAlmostEqual(out[5][2], 0.3, places=9)  # the string moves toward the tip relative to the bow at the instrument's speed: relative motion, as it must

    def test_a_bow_slanted_across_the_string_is_not_read_as_square(self):
        out = D.derive(rig(tilt_deg=15.0))
        self.assertAlmostEqual(out[0][5], 75.0, places=6)  # 90 - 15 degrees between the lines

    def test_a_bow_parallel_to_the_string_is_refused_not_guessed(self):
        rows = rig()
        for r in rows:
            r[1], r[2], r[3] = r[4], r[5] + 0.69, r[6]  # tip straight along +y from the frog
        with self.assertRaises(ValueError):
            D.derive(rows)

    def test_the_gradient_is_central_inside_and_one_sided_at_the_ends(self):
        t = [0.0, 1.0, 3.0]
        self.assertEqual(D.gradient([0.0, 1.0, 5.0], t), [1.0, 5.0 / 3.0, 2.0])


class Files(unittest.TestCase):
    def test_a_file_with_another_header_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "x.csv"
            p.write_text("a,b,c\n1,2,3\n")
            with self.assertRaises(SystemExit):
                D.read(p)

    def test_round_trip_writes_the_schema_line_columns_and_one_row_per_frame(self):
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = Path(tmp) / "in.csv", Path(tmp) / "out.tsv"
            with src.open("w", newline="") as f:
                w = csv.writer(f)
                w.writerow(D.EXPECTED_HEADER)
                w.writerows(rig(n=30))
            D.write(D.derive(D.read(src)), dst)
            lines = dst.read_text().splitlines()
            self.assertEqual((lines[0], lines[1].split("\t")), (D.SCHEMA, D.COLUMNS))
            self.assertEqual(len(lines), 2 + 30)
            self.assertEqual(lines[2].split("\t")[0], "0")


if __name__ == "__main__":
    unittest.main()
