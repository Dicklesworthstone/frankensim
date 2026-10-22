"""chain.py on synthetic floors and signals whose answers are fixed by construction (the session test is test_chain_session.py).

uv run --with numpy --with matplotlib python3 -I test_chain.py -v      (ffmpeg on PATH: the chain's resampler is convert_refs.py's)
Each test prints what it built and what it measured.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chain  # noqa: E402
from chain import FS_REC, Chain, ChainError, Stretch  # noqa: E402
from instrgrade import SR  # noqa: E402


def level_db(x):
    return 10 * math.log10(float(np.mean(np.asarray(x) ** 2)) + 1e-30)


def floor_noise(seconds, db, seed, rumble=False):
    """Stationary floor at `db` dB re full-scale power: white, or (rumble) mostly below 100 Hz, as an air handler's."""
    rng = np.random.default_rng(seed)
    n = int(seconds * FS_REC)
    x = rng.standard_normal(n)
    if rumble:
        X = np.fft.rfft(x)
        f = np.fft.rfftfreq(n, 1 / FS_REC)
        X *= 1 / np.sqrt(1 + (f / 60.0) ** 4)
        x = np.fft.irfft(X, n) + 0.05 * rng.standard_normal(n)
    return x * 10 ** (db / 20) / np.sqrt(np.mean(x * x))


def synthetic_chain(seconds=20.0, db=-80.0, seed=1, bits=16):
    x = floor_noise(seconds, db, seed)
    return Chain(id="synthetic@test", page="test", record={}, stretches=[Stretch("amb", "0" * 64, 0, len(x), "ambient")],
                 flags=[dict(flag="synthetic")], bits=bits, sources=[dict(path="amb", sha256="0" * 64, role="ambient")], _audio={"amb": x})


class QuantizeTests(unittest.TestCase):
    def test_chain_quantization_bit_exact(self):
        """G0: the declared rule, round half to even on the 2^-(bits-1) grid, no dither; clips counted, never wrapped."""
        lsb = 1 / 32768
        x = np.array([0.5, 1.5, -0.5, -1.5, 2.5, 0.49999, -2.5, 32767.4, 32767.5, -32768.6]) * lsb
        q, clips = chain.quantize(x, 16)
        print(f"\n  16 bit: {x / lsb} LSB -> {q.tolist()}, clips {clips}")
        self.assertEqual(q.tolist(), [0, 2, 0, -2, 2, 0, -2, 32767, 32767, -32768])
        self.assertEqual(clips, 2)  # 32767.5 rounds to 32768 and -32768.6 to -32769: both out of range, both counted
        q24, c24 = chain.quantize(np.array([0.5, 1.5, -1.0]) / 8388608, 24)
        self.assertEqual(q24.tolist(), [0, 2, -1])
        self.assertEqual(c24, 0)


class GapRuleTests(unittest.TestCase):
    def estimate(self, x):
        clean = chain.exclude_events([Stretch("g", "0", 0, len(x), "gap")], {"g": x})
        kept = np.concatenate([x[s.start:s.end] for s in clean])
        return level_db(kept), len(kept) / len(x)

    def test_chain_gap_floor_excludes_events(self):
        """G3: breath and key events injected into a stationary floor are excluded, so the floor estimate is unbiased within
        0.5 dB, for a white floor and for an air-handler rumble whose 50 ms levels scatter by dBs."""
        for rumble in (False, True):
            floor = floor_noise(30.0, -80.0, seed=3, rumble=rumble)
            x = floor.copy()
            rng = np.random.default_rng(4)
            t = np.arange(len(x)) / FS_REC
            for c in np.arange(1.5, 29.0, 3.0):  # a breath every 3 s: 0.4 s of noise 12 dB over the floor, raised-cosine edges
                w = np.clip(1 - np.abs(t - c) / 0.2, 0, 1)
                x += np.sqrt(0.5 - 0.5 * np.cos(np.pi * w)) * rng.standard_normal(len(x)) * 10 ** (-68 / 20)
            for c in np.arange(0.7, 29.0, 1.7):  # a key click every 1.7 s: 3 ms, 25 dB over the floor
                i = int(c * FS_REC)
                x[i:i + int(0.003 * FS_REC)] += rng.standard_normal(int(0.003 * FS_REC)) * 10 ** (-55 / 20)
            truth = level_db(floor)
            raw = level_db(x)
            est, kept = self.estimate(x)
            print(f"\n  {'rumble' if rumble else 'white'} floor {truth:.2f} dB; with events {raw:.2f} dB; after the gap rule {est:.2f} dB "
                  f"({100 * kept:.0f} % of the gap kept)")
            self.assertGreater(abs(raw - truth), 0.5, "the events bias an unfiltered estimate: the test discriminates")
            self.assertLess(abs(est - truth), 0.5)
            clean_est, clean_kept = self.estimate(floor)
            print(f"  the same floor without events: {clean_est:.2f} dB, {100 * clean_kept:.0f} % kept")
            self.assertLess(abs(clean_est - truth), 0.5)  # the rule does not shave a stationary floor's own fluctuation

    def test_gaps_leave_each_note_its_lead_and_tail(self):
        g = chain.gaps(10 * FS_REC, [(2.0, 3.0), (5.0, 6.0)])
        s = [(a / FS_REC, b / FS_REC) for a, b in g]
        print(f"\n  gaps of a 10 s file with notes at 2-3 s and 5-6 s: {s}")
        self.assertEqual([(round(a, 3), round(b, 3)) for a, b in s], [(0.0, 1.85), (3.5, 4.85), (6.5, 10.0)])


class FloorTests(unittest.TestCase):
    def test_chain_deterministic(self):
        """G5: the same (description, chain id, block) keys give the same samples, bit for bit; another description, others."""
        ch = synthetic_chain()
        a, pa = ch.floor(int(3.2 * FS_REC), "clarinet D3 mf render abc")
        b, pb = ch.floor(int(3.2 * FS_REC), "clarinet D3 mf render abc")
        c, pc = ch.floor(int(3.2 * FS_REC), "clarinet D3 mf render abd")
        print(f"\n  block starts: {[p['start'] for p in pa]}  /  another description: {[p['start'] for p in pc]}")
        self.assertEqual(hashlib.sha256(a.tobytes()).hexdigest(), hashlib.sha256(b.tobytes()).hexdigest())
        self.assertEqual(pa, pb)
        self.assertNotEqual([p["start"] for p in pa], [p["start"] for p in pc])

    def test_crossfaded_blocks_keep_the_floor_level(self):
        ch = synthetic_chain()
        fl, picks = ch.floor(10 * FS_REC, "level")
        src = level_db(ch._audio["amb"])
        lv = chain.frame_levels(fl, int(0.05 * FS_REC))
        print(f"\n  source floor {src:.2f} dB, assembled {level_db(fl):.2f} dB over {len(picks)} blocks; 50 ms levels "
              f"{lv.min():.2f}..{lv.max():.2f} dB")
        self.assertLess(abs(level_db(fl) - src), 0.2)
        self.assertLess(lv.max() - lv.min(), 3.0)  # no join is a hole or a bump

class ResampleTests(unittest.TestCase):
    def test_chain_resample_declared_filter(self):
        """G2: the declared filter is convert_refs.py's (ffmpeg aresample, swr defaults). Declared spec: within the grader's band
        (to 16 kHz) a 48 -> 44.1 -> 48 kHz round trip keeps every tone's level within 0.1 dB; content between 22.05 and 24 kHz
        can alias only above 20.1 kHz, outside that band, so its attenuation is measured and reported, not required."""
        t = np.arange(int(1.0 * SR)) / SR
        worst = 0.0
        for f in (100.0, 1000.0, 5000.0, 10000.0, 14000.0, 16000.0):
            x = 0.25 * np.sin(2 * np.pi * f * t)
            y = chain.resample(chain.resample(x, SR, FS_REC), FS_REC, SR)
            n = min(len(x), len(y))
            d = level_db(y[2000:n - 2000]) - level_db(x[2000:n - 2000])
            worst = max(worst, abs(d))
            print(f"\n  {f:7.0f} Hz: round-trip level change {d:+.4f} dB", end="")
        x = 0.25 * np.sin(2 * np.pi * 23000.0 * t)
        y = chain.resample(x, SR, FS_REC)
        Y = np.abs(np.fft.rfft(y[1000:-1000] * np.hanning(len(y) - 2000)))
        fy = np.fft.rfftfreq(len(y) - 2000, 1 / FS_REC)
        alias = 20 * math.log10(Y.max() / (np.abs(np.fft.rfft(x[1000:-1000] * np.hanning(len(x) - 2000))).max()) + 1e-30)
        print(f"\n  23 kHz tone after 48 -> 44.1 kHz: strongest line {fy[np.argmax(Y)]:.0f} Hz at {alias:.1f} dB (reported)")
        self.assertLess(worst, 0.1)


def aiff_bytes(x, rate=FS_REC):
    """`x` as a 16-bit mono AIFF at `rate`, written by ffmpeg to a file (a pipe leaves the header's sizes unset)."""
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "x.aiff"
        chain._ffmpeg(["-f", "f64le", "-ar", str(rate), "-ac", "1", "-i", "pipe:0", "-c:a", "pcm_s16be", "-y", str(p)],
                      np.asarray(x, dtype="<f8").tobytes())
        return p.read_bytes()


class RefusalTests(unittest.TestCase):
    def manifest(self, d, ambient_bytes=None, sha=None):
        amb = Path(d) / "refs" / "Amb.ambient.silence.aiff"
        amb.parent.mkdir(parents=True)
        if ambient_bytes is not None:
            amb.write_bytes(ambient_bytes)
            sha = sha or hashlib.sha256(ambient_bytes).hexdigest()
        m = dict(chains=[dict(id="c1", era="pre-2011", mic="KM 84", format="16-bit/44.1 kHz mono", signal_chain=["Mackie mixer"],
                              low_cutoff_status="undocumented")],
                 pages=[dict(page="P.html", chain="c1", ambient="amb1")],
                 files=[dict(id="amb1", path="refs/Amb.ambient.silence.aiff", page="P.html", kind="ambient-silence", gradable=True,
                             channels="mono", sha256=sha or "f" * 64)])
        p = Path(d) / "references.json"
        p.write_text(json.dumps(m))
        return p

    def test_chain_refuses_missing_ambient(self):
        """G0: a page whose ambient file is not there, or not the pinned bytes, is refused by name; so is an unknown session."""
        with tempfile.TemporaryDirectory() as d:
            p = self.manifest(d)
            with self.assertRaises(ChainError) as e:
                Chain.from_references(p, "P.html", Path(d), [])
            print(f"\n  absent: {e.exception}")
            self.assertIn("IRPCHAIN_AMBIENT_MISSING", str(e.exception))
            with self.assertRaises(ChainError) as e:
                Chain.from_references(p, "Q.html", Path(d), [])
            self.assertIn("IRPCHAIN_UNKNOWN_SESSION", str(e.exception))
        with tempfile.TemporaryDirectory() as d:
            p = self.manifest(d, ambient_bytes=b"not the pinned bytes", sha="e" * 64)
            with self.assertRaises(ChainError) as e:
                Chain.from_references(p, "P.html", Path(d), [])
            print(f"  altered: {e.exception}")
            self.assertIn("IRPCHAIN_HASH_MISMATCH", str(e.exception))

    def test_an_ambient_file_builds_a_flagged_chain_and_one_too_short_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p = self.manifest(d, ambient_bytes=aiff_bytes(floor_noise(3.0, -75.0, 6)))
            ch = Chain.from_references(p, "P.html", Path(d), [])
            flags = sorted(f["flag"] for f in ch.flags)
            print(f"\n  chain {ch.id}: {ch.bits} bit, {len(ch.stretches)} clean stretch(es), flags {flags}")
            self.assertEqual(ch.bits, 16)
            self.assertIn("mixer_eq_undocumented", flags)
            self.assertIn("low_cutoff_undocumented", flags)
            self.assertIn("capsule_not_applied", flags)
            self.assertNotIn("ambient_absent_on_page", flags)
        with tempfile.TemporaryDirectory() as d:
            p = self.manifest(d, ambient_bytes=aiff_bytes(floor_noise(0.4, -75.0, 7)))
            with self.assertRaises(ChainError) as e:
                Chain.from_references(p, "P.html", Path(d), [])
            print(f"  0.4 s of ambient: {e.exception}")
            self.assertIn("IRPCHAIN_FLOOR_TOO_SHORT", str(e.exception))

    def test_a_session_recorded_at_96_khz_keeps_its_own_rate(self):
        """The 2011 guitar .mono files are 96 kHz although their chain record says 44.1: the probed rate is used and flagged."""
        with tempfile.TemporaryDirectory() as d:
            p = self.manifest(d, ambient_bytes=aiff_bytes(floor_noise(3.0, -75.0, 8), rate=96000))
            ch = Chain.from_references(p, "P.html", Path(d), [])
            sig = 0.3 * np.sin(2 * np.pi * 440 * np.arange(2 * SR) / SR)
            y, r = ch.apply(sig, "a 440 Hz tone at 96 kHz")
            n = min(len(sig), len(y))
            added = level_db(y[1000:n - 1000] - sig[1000:n - 1000])
            # a white floor at 96 kHz spans 0-48 kHz; the grader's 48 kHz keeps the half below 24 kHz
            expect = -75.0 - 10 * math.log10(96000 / SR)
            print(f"\n  chain {ch.id}: {ch.bits} bit at {ch.fs} Hz, flags {sorted(f['flag'] for f in ch.flags)}; floor -75 dB over 0-48 kHz, "
                  f"{expect:.2f} dB below 24 kHz; added {added:.2f} dB")
            self.assertEqual(ch.fs, 96000)
            self.assertEqual(r["fs_hz"], 96000)
            self.assertIn("rate_differs_from_record", [f["flag"] for f in ch.flags])
            self.assertLess(abs(added - expect), 0.5)

    def test_a_render_too_loud_for_the_session_is_refused(self):
        ch = synthetic_chain()
        with self.assertRaises(ChainError) as e:
            ch.apply(1.2 * np.sin(2 * np.pi * 440 * np.arange(SR) / SR), "loud")
        print(f"\n  {e.exception}")
        self.assertIn("IRPCHAIN_CLIPS", str(e.exception))

    def test_an_applied_chain_adds_the_floor_at_its_recorded_level_and_says_so(self):
        ch = synthetic_chain(db=-70.0)
        sig = 0.3 * np.sin(2 * np.pi * 440 * np.arange(2 * SR) / SR)
        y, r = ch.apply(sig, "a 440 Hz tone")
        n = min(len(sig), len(y))
        resid = y[1000:n - 1000] - sig[1000:n - 1000]
        print(f"\n  floor -70.0 dB; what the chain added to the tone: {level_db(resid):.2f} dB; receipt keys {sorted(r)}")
        self.assertLess(abs(level_db(resid) - (-70.0)), 0.5)
        for k in ("chain_id", "ambient_sha256", "flags", "bits", "fs_hz", "floor_db_per_band", "gap_rule", "blocks"):
            self.assertIn(k, r)


if __name__ == "__main__":
    unittest.main()
