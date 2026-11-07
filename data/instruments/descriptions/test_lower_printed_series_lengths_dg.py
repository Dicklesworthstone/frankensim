#!/usr/bin/env python3
"""G0/G3 regressions for the printed series length t_a lowered as a bore shortening on a DG card.

Run from the repository root with:
    python3 -B data/instruments/descriptions/test_lower_printed_series_lengths_dg.py
The reference numbers are the hand-computed values in fs-duct/src/hole_lengths/tests.rs (a hole of radius 4 mm on a 7.5 mm bore,
chimney 6 mm): Eq. 33 t_a(open) = -3.29996986190055e-4 m and Eq. 34 t_a(closed) = -3.29666958445790e-4 m. The geometry tests use a
manufactured DG card whose shortening is worked by hand. These tests exercise the lowering, not acoustic or specimen validation.
"""
import importlib.util
import math
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location(
    "lower_printed_series_lengths_dg", Path(__file__).with_name("lower_printed_series_lengths_dg.py")
)
lower = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lower)

B, A, T = 4.0e-3, 7.5e-3, 6.0e-3
TA_OPEN = -3.29996986190055e-4
TA_CLOSED = -3.29666958445790e-4


def live_card(open_flags):
    rows = "\n".join(f"hole {B!r} {T!r} {A!r} {int(flag)}" for flag in open_flags)
    return f"frankensim-reed-performance-v1\nreed 1 2 3\n{rows}\ncompile_limits 1 2\n"


def dg_card(segment_lengths, faces, gas, elements=None, xs=None, closed=(3,)):
    """A DG card: cylinders of the given lengths, two chimneys on the faces (z from the reed), gas stations at z."""
    total = sum(segment_lengths)
    k = elements if elements is not None else len(segment_lengths) + len([z for z in gas if 1e-9 < z < total - 1e-9])
    lines = ["frankensim-material-modal-dg-valve-v1", "audio 96000 96000 100 20",
             f"dg_bore {len(segment_lengths)} baffled {k} 1 4 50 2000 0.01 1e-11 1000 2 min_poles"]
    lines += [f"dg_segment cylinder {A!r} {length!r}" for length in segment_lengths]
    lines += [f"gas_station {z!r} 300.0 1.0 0.1" for z in gas]
    lines.append(f"dg_side_openings {len(faces)}")
    for n, face in enumerate(faces):
        x = xs[n] if xs else total - face
        lines.append(f"dg_side_chimney {n + 2} {n + 3} {x!r} 0 {B!r} {T!r} 1 1")
        lines.append("dg_chimney_gas pack x")
    lines.append(f"dg_terminal_channels {len(closed)} 0")
    lines += [f"dg_terminal_closed {p}" for p in closed]
    lines.append("dg_radiation x")
    return "\n".join(lines) + "\n"


SEGMENTS = [0.10, 0.05, 0.05, 0.10]  # faces at 0.10 (hole 1) and 0.20 (hole 2); total 0.30
GAS = [0.0, 0.075, 0.30]  # one station inside segment 1: 5 spans, 5 elements


class PrintedSeriesLengths(unittest.TestCase):
    def test_the_equations_are_the_hand_computed_values(self):
        self.assertAlmostEqual(lower.series_length_m(B, T, A, True), TA_OPEN, delta=1e-12)
        self.assertAlmostEqual(lower.series_length_m(B, T, A, False), TA_CLOSED, delta=1e-12)
        self.assertLess(TA_OPEN, 0.0)
        self.assertLess(TA_CLOSED, 0.0)

    def test_the_bore_is_shortened_at_each_face_and_everything_follows_by_hand(self):
        text, rows, total, new_total = lower.lower(live_card([True, False]), dg_card(SEGMENTS, [0.10, 0.20], GAS))
        d1, d2 = -TA_OPEN, -TA_CLOSED
        self.assertAlmostEqual(total - new_total, d1 + d2, delta=1e-12)
        lines = text.split("\n")
        seg = [float(l.split()[3]) for l in lines if l.startswith("dg_segment ")]
        for got, want in zip(seg, [0.10 - d1, 0.05, 0.05 - d2, 0.10]):
            self.assertAlmostEqual(got, want, delta=1e-12)
        gas = [float(l.split()[1]) for l in lines if l.startswith("gas_station ")]
        self.assertEqual(gas[0], 0.0)
        self.assertAlmostEqual(gas[1], 0.75 * (0.10 - d1), delta=1e-12)  # the same relative place inside segment 1
        self.assertAlmostEqual(gas[2], 0.30 - d1 - d2, delta=1e-12)  # the engine needs the last station at the bore length
        x = [float(l.split()[3]) for l in lines if l.startswith("dg_side_chimney ")]
        self.assertAlmostEqual(x[0], 0.20 - d2, delta=1e-12)  # a hole upstream of a shortening moves toward the bell by it
        self.assertAlmostEqual(x[1], 0.10, delta=1e-12)  # nothing downstream of hole 2 was shortened
        # the face station (bore length minus X_M) is the new segment boundary
        bounds = [0.0]
        for length in seg:
            bounds.append(bounds[-1] + length)
        self.assertAlmostEqual(new_total - x[0], bounds[1], delta=1e-12)
        self.assertAlmostEqual(new_total - x[1], bounds[3], delta=1e-12)

    def test_only_segments_stations_and_hole_positions_change(self):
        before = dg_card(SEGMENTS, [0.10, 0.20], GAS)
        text, *_ = lower.lower(live_card([True, False]), before)
        old, new = before.split("\n"), text.split("\n")
        self.assertEqual(len(old), len(new))
        changed = {o.split()[0] for o, n in zip(old, new) if o != n}
        self.assertEqual(changed, {"dg_segment", "gas_station", "dg_side_chimney"})
        for o, n in zip(old, new):
            if o.startswith("dg_side_chimney "):
                a, b = o.split(" "), n.split(" ")
                self.assertEqual([a[i] for i in (0, 1, 2, 4, 5, 6, 7, 8)], [b[i] for i in (0, 1, 2, 4, 5, 6, 7, 8)])

    def test_a_card_that_breaks_an_invariant_is_refused_by_name(self):
        live = live_card([True, False])
        cases = {
            "not a segment boundary": dg_card(SEGMENTS, [0.10, 0.20], GAS, xs=[0.2005, 0.10]),
            "spans for": dg_card(SEGMENTS, [0.10, 0.20], GAS, elements=7),
            "dg_terminal_closed list says the opposite": dg_card(SEGMENTS, [0.10, 0.20], GAS, closed=(2,)),
            "more than half": dg_card([0.10, 0.05, 0.0002, 0.1498], [0.10, 0.1502], GAS, xs=[0.2, 0.1498]),
            "sits on a segment boundary": dg_card(SEGMENTS, [0.10, 0.20], [0.0, 0.10, 0.30]),
        }
        for message, card in cases.items():
            with self.assertRaises(SystemExit) as caught:
                lower.lower(live, card)
            self.assertIn(message, str(caught.exception.code), message)
        with self.assertRaises(SystemExit) as caught:
            lower.lower(live_card([True]), dg_card(SEGMENTS, [0.10, 0.20], GAS))
        self.assertIn("1 holes and the DG card 2 chimneys", str(caught.exception.code))


if __name__ == "__main__":
    unittest.main()
