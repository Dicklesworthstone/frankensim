#!/usr/bin/env python3
"""G0/G3 regressions for the printed tone-hole lengths on a DG chimney card.

Run from the repository root with:
    python3 -B data/instruments/descriptions/test_lower_printed_hole_lengths_dg.py
The reference numbers are the hand-computed values in fs-duct/src/hole_lengths/tests.rs (a hole of radius 4 mm on a 7.5 mm bore,
chimney 6 mm, d = 8/15), evaluated there from the Lefebvre and Scavone (2012) equations in plain power form. These tests exercise
the lowering, not acoustic or specimen validation.
"""
import importlib.util
import math
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location(
    "lower_printed_hole_lengths_dg", Path(__file__).with_name("lower_printed_hole_lengths_dg.py")
)
lower = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lower)

B, A, T = 4.0e-3, 7.5e-3, 6.0e-3
HAND_INNER = 0.539207619423868 * B       # Eq. 31 times b, G = 1
HAND_MATCHING = 2.75040711111111e-4      # Eq. 8
HAND_RADIATION = 2.57754131520697e-3     # Eq. 11


def live_card(open_flags):
    rows = "\n".join(f"hole {B!r} {T!r} {A!r} {int(flag)}" for flag in open_flags)
    return f"frankensim-reed-performance-v1\nreed 1 2 3\n{rows}\ncompile_limits 1 2\n"


def dg_card(open_flags):
    chimneys = "\n".join(f"dg_side_chimney {k + 2} {k + 5} 0.1 0 {B!r} {T!r} 1 1\ndg_chimney_gas pack x" for k in range(len(open_flags)))
    closed = [k + 2 for k, flag in enumerate(open_flags) if not flag]
    block = f"dg_terminal_channels {len(closed)} 0\n" + "".join(f"dg_terminal_closed {p}\n" for p in closed)
    return f"frankensim-material-modal-dg-valve-v1\ndg_side_openings {len(open_flags)}\n{chimneys}\n{block}dg_radiation x\n"


class PrintedLengths(unittest.TestCase):
    def test_the_equations_are_the_hand_computed_values(self):
        d = B / A
        self.assertTrue(math.isclose(lower.inner_over_b(d) * B, HAND_INNER, rel_tol=1e-12))
        self.assertTrue(math.isclose(lower.matching_m(B, d), HAND_MATCHING, rel_tol=1e-12))
        self.assertTrue(math.isclose(lower.radiation_m(B, A, T), HAND_RADIATION, rel_tol=1e-12))
        self.assertTrue(math.isclose(
            lower.open_extra_m(B, A, T),
            HAND_INNER + HAND_MATCHING + (HAND_RADIATION - 0.8216 * B), rel_tol=1e-12))
        self.assertTrue(math.isclose(lower.sealed_extra_m(B, A), HAND_MATCHING, rel_tol=1e-12))

    def test_only_the_chimney_lengths_change_and_by_state(self):
        flags = [False, True, False, True]
        dg = dg_card(flags)
        out, rows = lower.lower(live_card(flags), dg)
        before, after = dg.split("\n"), out.split("\n")
        self.assertEqual(len(before), len(after))
        changed = [k for k, (x, y) in enumerate(zip(before, after)) if x != y]
        self.assertEqual([before[k].split(" ")[0] for k in changed], ["dg_side_chimney"] * 4)
        for k, (x, y) in zip(changed, [(before[k], after[k]) for k in changed]):
            fx, fy = x.split(" "), y.split(" ")
            self.assertEqual(fx[:6] + fx[7:], fy[:6] + fy[7:])  # every field but LENGTH
        for (hole, port, is_open, b, a, t, extra, new), line in zip(rows, [after[k] for k in changed]):
            self.assertEqual(float(line.split(" ")[6]), new)
            self.assertEqual(new, t + extra)
            self.assertEqual(extra, lower.open_extra_m(b, a, t) if is_open else lower.sealed_extra_m(b, a))
        self.assertGreater(rows[1][6], rows[0][6])  # an open hole grows by more than a closed one

    def test_a_card_that_disagrees_with_the_live_holes_is_refused_by_name(self):
        flags = [False, True]
        live, dg = live_card(flags), dg_card(flags)
        for bad_live, bad_dg, word in [
            (live_card([False, True, True]), dg, "3 holes"),
            (live_card([True, True]), dg, "dg_terminal_closed"),
            (live.replace(repr(B), "0.0041", 1), dg, "does not match"),
            (live_card(flags).replace(repr(A), repr(B), 1), dg, "inside the bore radius"),
            (live, dg.replace("frankensim-material-modal-dg-valve-v1", "x"), "first line"),
            (live.replace(" 0\n", " 2\n", 1), dg, "malformed"),
        ]:
            with self.assertRaises(SystemExit) as caught:
                lower.lower(bad_live, bad_dg)
            self.assertIn(word, str(caught.exception))


if __name__ == "__main__":
    unittest.main()
