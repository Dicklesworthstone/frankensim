"""The corpus-mode chain of classify_violin_modes.py leaves out a mode whose bridge kinetic share exceeds one half, and keeps one at or under it.

Run with: uv run --python 3.12 --with numpy python -m unittest discover -s crates/fs-couple/examples/shell_body_modes -p 'test_classify_violin_modes.py' -v
The mode rows below are written by hand: the rule chain and the share parser are exercised, not a body run.
"""
import unittest

import classify_violin_modes as C


def row(rank, f, bridge, neck=0.0, dv=0.0, odd=0.0, waist=0.0, sign="same"):
    return dict(rank=rank, f=f, neck=neck, dv=dv, tw=0.0, bw=0.0, q=50.0, qrad=100.0, qtot=33.0, odd=odd, waist=waist, sign=sign, bridge=bridge)


def a0():
    return row(1, 243.0, 0.04, neck=0.95, dv=0.16, odd=0.25, waist=0.22)


def cbr_like(rank, f, bridge):
    """Odd in x and concentrated at the waist: what the CBR rule asks for."""
    return row(rank, f, bridge, odd=0.99, waist=0.59)


class BridgeShareThreshold(unittest.TestCase):
    def test_threshold_is_one_half(self):
        self.assertEqual(C.BRIDGE_SHARE_MAX, 0.5)

    def test_a_share_just_above_one_half_leaves_the_chain(self):
        found = C.identify([a0(), cbr_like(2, 400.0, 0.5001), cbr_like(3, 476.0, 0.4999)])
        self.assertTrue(C.is_bridge_mode(cbr_like(2, 400.0, 0.5001)))
        self.assertEqual(found["cbr"]["rank"], 3)

    def test_a_share_just_below_one_half_stays_in_the_chain(self):
        found = C.identify([a0(), cbr_like(2, 400.0, 0.4999), cbr_like(3, 476.0, 0.10)])
        self.assertFalse(C.is_bridge_mode(cbr_like(2, 400.0, 0.4999)))
        self.assertEqual(found["cbr"]["rank"], 2)

    def test_exactly_one_half_does_not_exceed_it(self):
        self.assertFalse(C.is_bridge_mode(cbr_like(2, 400.0, 0.5)))

    def test_a_bridge_mode_cannot_be_b1_minus(self):
        # even, plates with opposite signs, small volume change: the B1- rule's shape, carried by the bridge
        b1m_shape = dict(odd=0.06, sign="opposite", dv=0.007)
        rows = [a0(), cbr_like(2, 476.0, 0.0), row(3, 695.6, 0.88, **b1m_shape), row(4, 720.0, 0.02, **b1m_shape)]
        found = C.identify(rows)
        self.assertEqual(found["b1m"]["rank"], 4)
        rows[2]["bridge"] = 0.12
        self.assertEqual(C.identify(rows)["b1m"]["rank"], 3)

    def test_without_shares_nothing_is_left_out(self):
        rows = [a0(), cbr_like(2, 400.0, None)]
        self.assertFalse(C.is_bridge_mode(rows[1]))
        self.assertEqual(C.identify(rows)["cbr"]["rank"], 2)


class ShareBlockParser(unittest.TestCase):
    def test_reads_the_bridge_column_by_name(self):
        text = [
            "open sound hole (shell modes + collective air spring + sound-hole plug):",
            "kinetic energy share of the shell parts by group, open sound hole (lumped nodal mass):",
            "  mode   f_Hz       back       bass     bridge        rib  soundpost        top",
            "     1   243.23      0.440      0.050      0.040      0.020      0.000      0.450",
            "     6   695.65      0.070      0.000      0.880      0.010      0.000      0.040",
            "",
        ]
        self.assertEqual(C.bridge_shares(text), {1: 0.04, 6: 0.88})

    def test_an_output_without_the_block_gives_none(self):
        self.assertIsNone(C.bridge_shares(["open sound hole (shell modes + collective air spring + sound-hole plug):", ""]))


if __name__ == "__main__":
    unittest.main()
