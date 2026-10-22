"""identity.py (the current Q2 and the UNREGISTERED candidates A, B, C) on synthetic instruments whose life is known by construction.

uv run --with numpy --with matplotlib python3 -I test_identity.py -v
Every test prints what it built, the set it was judged against and the verdict with each measure's margin, so a failure log
says why.
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import decompose as dc  # noqa: E402
import freeze  # noqa: E402
import identity as idn  # noqa: E402
from decompose import SR  # noqa: E402

FORMANTS_A = ((600.0, 12.0, 0.5), (2500.0, 9.0, 0.4))  # (centre Hz, height dB, width in octaves): instrument A's resonances
FORMANTS_B = ((350.0, 10.0, 0.4), (1400.0, 12.0, 0.3), (4000.0, 6.0, 0.5))  # instrument B's


def played(midi, dyn=0.0, seed=0, dur=2.4, jitter_cents=None, shimmer_db=1.0, noise_db=None, life=True, formants=FORMANTS_A):
    """One note of a synthetic instrument that lives: a fixed formant envelope (so its notes share an identity), a brighter tilt
    at louder dynamics (`dyn` -1..1), a 40 ms attack, f0 wandering as low-passed noise (`jitter_cents` rms, default 6 c, more at
    pp), each partial's level wandering on its own (`shimmer_db` rms), and breath noise `noise_db` under the note (default
    -28 - 4 dyn dB: a soft note is breathier) shaped by the same envelope. life=False: no wandering, no breath (a periodic render)."""
    rng = np.random.default_rng(seed)
    noise_db = -28.0 - 4.0 * dyn if noise_db is None else noise_db
    jitter_cents = 6.0 * (1.0 - 0.3 * dyn) if jitter_cents is None else jitter_cents
    n = int(dur * SR)
    t = np.arange(n) / SR
    f0 = 440.0 * 2 ** ((midi - 69) / 12)

    def wander(rms, hz):
        w = rng.standard_normal(n)
        W = np.fft.rfft(w)
        W[np.fft.rfftfreq(n, 1 / SR) > hz] = 0
        w = np.fft.irfft(W, n)
        return rms * w / np.sqrt(np.mean(w * w))

    cents = wander(jitter_cents, 12.0) if life else np.zeros(n)
    phase = 2 * np.pi * np.cumsum(f0 * 2 ** (cents / 1200)) / SR

    def envelope_db(f):
        peaks = sum(h * np.exp(-0.5 * (np.log2(f / c) / w) ** 2) for c, h, w in formants)
        return peaks - 9 * (1 - 0.3 * dyn) * np.log2(np.maximum(f, 100.0) / 300.0)

    x = np.zeros(n)
    for k in range(1, int(14000 // f0) + 1):
        lv = 10 ** ((envelope_db(k * f0) + (wander(shimmer_db, 8.0) if life else 0.0)) / 20)
        x += lv * np.sin(k * phase + 0.4 * k)
    a = np.clip(t / 0.04, 0, 1) * np.clip((dur - t) / 0.05, 0, 1)
    x *= a
    if life:
        X = np.fft.rfft(rng.standard_normal(n))
        f = np.fft.rfftfreq(n, 1 / SR)
        X *= 10 ** (envelope_db(np.maximum(f, 50.0)) / 20)
        b = np.fft.irfft(X, n) * a
        x += b * np.sqrt(np.mean(x * x) / np.mean(b * b)) * 10 ** (noise_db / 20)
    lead = np.zeros(int(0.2 * SR))
    return 0.25 * np.concatenate([lead, x]) / np.abs(x).max()


def feats(y, midi):
    """A signal graded as a candidate is: the render floor, the alignment, the grader's Note, then identity.features."""
    f0 = 440.0 * 2 ** ((midi - 69) / 12)
    z, info = dc.align(y, f0, dc.render_floor(dc.highpass(y, f0), f0))
    return idn.features(dc.Note(z, f0, "sustained", False, "synthetic", info, max_cents=400.0))


def a_set(seed0, **kw):
    """3 pitches x 3 dynamics, as a library holds a real instrument's notes around a candidate's pitch."""
    return [(f"m{m}-d{d:+.0f}", feats(played(m, d, seed=seed0 + 10 * m + int(d), **kw), m)) for m in (67, 69, 71) for d in (-1.0, 0.0, 1.0)]


def shift(y, cents):
    """`y` played `cents` higher: band-limited (FFT) resampling to len(y) / 2^(cents/1200) samples. A linear-interpolation
    resampler is not a pitch shift here: its interpolation error lifts the steep top of a note's spectrum by several dB."""
    n = len(y)
    m = int(round(n / 2 ** (cents / 1200)))
    return np.fft.irfft(np.fft.rfft(y)[: m // 2 + 1], m) * (m / n)


def show(name, j):
    ms = "  ".join(f'{m} {v["nearest_db"]:.2f}/{v["ruler_db"]:.2f} ({v["margin_db"]:+.2f})' for m, v in j["measures"].items()
                   if v.get("margin_db") is not None)
    print(f"\n  {name}: passes {j['recognisable']}, margin {j['margin_db']} dB over a set of {j['set_size']}; nearest/ruler (margin) dB: {ms}")



class SyntheticInstrumentTests(unittest.TestCase):
    """Instrument A, whose life is mostly in its partials' wandering, and instrument B, which only differs in its formants."""

    @classmethod
    def setUpClass(cls):
        cls.members = a_set(5)
        cls.others = {"B": a_set(7005, formants=FORMANTS_B)}
        cls.cand = played(69, 0.3, seed=999)
        print(f"\n  set: {len(cls.members)} living notes of A (G4, A4, B4 x pp, mf, ff; 4-8 c f0 wander, 1 dB partial wander, "
              f"breath at -24..-32 dB); {len(cls.others['B'])} of B (other formants)")

    def judge(self, name, y, midi=69, measures=idn.LIFE_MEASURES):
        j = idn.judge(feats(y, midi), self.members, measures)
        show(name, j)
        return j

    def test_q2_alive_is_the_life_rule_on_sustained_notes_and_makes_no_claim_on_struck_ones(self):
        """The registered statistic (pexfp.357): candidate B on a sustained note; no claim at all on a decaying one."""
        f, frozen = feats(self.cand, 69), feats(freeze.split(self.cand, 440.0).frozen, 69)
        alive, still, struck = idn.q2_alive(f, self.members, "sustained"), idn.q2_alive(frozen, self.members, "sustained"), \
            idn.q2_alive(f, self.members, "decaying")
        print(f"\n  {idn.STATISTIC}: living A4 {alive['recognisable']} ({alive['margin_db']:+.2f} dB), frozen {still['recognisable']} "
              f"({still['margin_db']:+.2f} dB), as a decaying note {struck['recognisable']} ({struck['reason']})")
        self.assertEqual(alive["statistic"], idn.STATISTIC)
        self.assertTrue(alive["recognisable"])
        self.assertEqual(alive["margin_db"], idn.judge(f, self.members, idn.LIFE_MEASURES)["margin_db"])
        self.assertFalse(still["recognisable"])
        self.assertIsNone(struck["recognisable"])
        self.assertIn("no claim", struck["reason"])

    def test_life_passes_a_living_note_of_the_instrument(self):
        self.assertTrue(self.judge("B: another living A4 of A", self.cand)["recognisable"])

    def test_life_fails_the_same_note_frozen(self):
        s = freeze.split(self.cand, 440.0)
        j = self.judge("B: the same A4 frozen (freeze.split)", s.frozen)
        self.assertFalse(j["recognisable"])
        self.assertLess(j["measures"]["flux_level_db"]["margin_db"], -10.0)  # not a near miss: its stillness

    def test_life_and_set_membership_fail_a_periodic_render_of_the_instrument(self):
        y = played(69, 0.3, seed=999, life=False)
        self.assertFalse(self.judge("B: a periodic render of the A4 (no wander, no breath)", y)["recognisable"])
        self.assertFalse(self.judge("A: the same", y, measures=idn.SET_MEASURES)["recognisable"])

    def test_life_fails_the_four_tone_nulls(self):
        y = self.cand[int(0.19 * SR):]
        for name, z in freeze.null_signals(y, 440.0).items():
            if name != "frozen":
                self.assertFalse(self.judge(f"B: null {name}", z)["recognisable"], name)

    def test_the_current_definition_passes_a_frozen_note_that_life_fails(self):
        """A defect the redefinition is for: the bead's Q2 compares spectra only, and a frozen note keeps its recording's spectrum
        wherever the residual does not shape it. (Freezing also drops the residual: a breathy note's top bands go with it, and
        then the spectrum fails too. This note of A has its wander and no breath.)"""
        played_dry = played(69, 0.3, seed=999, noise_db=-90.0)
        s = freeze.split(played_dry, 440.0)
        ref, frozen = feats(played_dry, 69), feats(s.frozen, 69)
        cur = idn.bead_q2(frozen, ref, self.others["B"])
        print(f"\n  current Q2 on a breathless A4 frozen (reference: the same A4 unfrozen, others: B's notes): {cur}")
        self.assertTrue(cur["recognisable"])
        j = idn.judge(frozen, self.members, idn.LIFE_MEASURES)
        show("B: the same frozen A4", j)
        self.assertFalse(j["recognisable"])

    def test_identity_by_the_nearest_set_tells_a_from_b_and_life_does_not(self):
        own = idn.nearest_set(feats(self.cand, 69), self.members, self.others)
        b = feats(played(69, 0.3, seed=998, formants=FORMANTS_B), 69)
        other = idn.nearest_set(b, self.members, self.others)
        print(f"\n  C: A's A4 against A's set and B's: {own}\n  C: B's A4 against A's set and B's: {other}")
        self.assertTrue(own["own_nearest"])
        self.assertFalse(other["own_nearest"])
        self.assertEqual(other["nearest_other"], "B")
        self.assertTrue(idn.judge(b, self.members, idn.LIFE_MEASURES)["recognisable"], "life alone does not tell A from B")

    def test_a_twenty_cent_shift_moves_no_margin_by_a_db(self):
        fa, fb = feats(self.cand, 69), feats(shift(self.cand, 20.0), 69)
        d = {}
        for name, ms in (("A", idn.SET_MEASURES), ("B", idn.LIFE_MEASURES)):
            a, b = idn.judge(fa, self.members, ms), idn.judge(fb, self.members, ms)
            show(f"{name}: living A4", a)
            show(f"{name}: living A4 shifted +20 c", b)
            d.update({f"{name}:{m}": abs(a["measures"][m]["margin_db"] - b["measures"][m]["margin_db"]) for m in ms})
        ia, ib = idn.nearest_set(fa, self.members, self.others), idn.nearest_set(fb, self.members, self.others)
        d["C:" + idn.IDENTITY_MEASURE] = abs(ia["margin_db"] - ib["margin_db"])
        ca, cb = idn.bead_q2(fa, self.members[4][1], self.others["B"]), idn.bead_q2(fb, self.members[4][1], self.others["B"])
        d["current:spec_db"] = abs(ca["margin_db"] - cb["margin_db"])
        print("  |margin change| dB: " + "  ".join(f"{m} {v:.3f}" for m, v in d.items())
              + " (a verdict this near the edge of a synthetic set may flip: its flux ruler is a few tenths of a dB, tighter than any real set's)")
        self.assertLess(max(d.values()), 1.0)

    def test_a_set_too_small_makes_no_claim(self):
        j = idn.judge(feats(self.cand, 69), self.members[: idn.MIN_SET - 1], idn.LIFE_MEASURES)
        print(f"\n  set of {idn.MIN_SET - 1}: passes {j['recognisable']} ({j.get('reason')})")
        self.assertIsNone(j["recognisable"])

    def test_the_residual_is_seen_on_a_frozen_note(self):
        """Residual present vs absent. On the frozen note the residual is the only movement left, so life must see it come and go.
        (Removed from a note whose partials still wander, it barely moves the flux: the flux is an RMS over the cells that are in the
        map, and a partial's own band carries its share of the residual. And with it, a frozen note may pass: stationary noise moves
        a spectrum too.)"""
        s = freeze.split(self.cand, 440.0)
        bare, noisy = self.judge("B: A4 frozen, residual removed", s.frozen), self.judge("B: A4 frozen, residual put back", s.frozen + s.residual)
        rise = noisy["measures"]["flux_level_db"]["margin_db"] - bare["measures"]["flux_level_db"]["margin_db"]
        print(f"  residual {s.stats(self.cand)['residual_db']:.1f} dB re note; it raises the flux margin by {rise:+.1f} dB")
        self.assertFalse(bare["recognisable"])
        self.assertGreater(rise, 10.0)


if __name__ == "__main__":
    unittest.main()
