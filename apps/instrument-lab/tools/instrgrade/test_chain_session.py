"""chain.py against a real UIowa session's fetched files (the synthetic tests are test_chain.py).

  FRANKENSIM_REFERENCE_DIR   the directory fetch-references.sh filled (UIowa AIFFs, never committed)
  IRP_LIBRARY                decompose's library JSON (file, start, end of every segmented note), for the gaps
  IRP_REFS_WAV48             convert_refs.py's output, for the resampler check
  FRANKENSIM_REFERENCES_MANIFEST   optional; default data/instruments/references.json of this checkout

Without the files these tests FAIL with IRPCHAIN_AMBIENT_MISSING: they are never skipped.
uv run --with numpy --with matplotlib python3 -I test_chain_session.py -v
"""
from __future__ import annotations

import math
import os
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chain  # noqa: E402
from chain import FS_REC, Chain  # noqa: E402
from instrgrade import load  # noqa: E402

SESSION = "MISBbclarinet.html"  # the one session whose page names an ambient file (BbClar.ambient.silence.aiff)
BANDS = 20.0 * 2 ** (np.arange(0, 31) / 3)  # third octaves from 20 Hz to 20 kHz
TOL_DB = 1.0


def need(var):
    v = os.environ.get(var)
    if not v or not Path(v).exists():
        raise AssertionError(f"IRPCHAIN_AMBIENT_MISSING: set {var} to the fetched UIowa files (fetch-references.sh); this test is not skipped")
    return Path(v)


def manifest():
    return Path(os.environ.get("FRANKENSIM_REFERENCES_MANIFEST", Path(__file__).resolve().parents[4] / "data/instruments/references.json"))


def band_psd(x, fs, n=8192):
    """Per third-octave band (BANDS): (mean power in dB over Hann frames, the 95 % half-width of that mean from the frames' scatter)."""
    w = np.hanning(n)
    frames = np.lib.stride_tricks.sliding_window_view(x, n)[:: n // 2]
    P = np.abs(np.fft.rfft(frames * w, axis=1)) ** 2 / (w ** 2).sum()
    f = np.fft.rfftfreq(n, 1 / fs)
    mean, half = [], []
    for fc in BANDS:
        sel = (f >= fc * 2 ** (-1 / 6)) & (f < fc * 2 ** (1 / 6))
        pb = P[:, sel].sum(axis=1)  # each frame's band power
        m = pb.mean()
        # frames overlap by half: about len/2 independent ones
        se = pb.std(ddof=1) / math.sqrt(max(1.0, len(pb) / 2))
        mean.append(10 * math.log10(m + 1e-30))
        half.append(10 * math.log10((m + 1.96 * se) / m))
    return np.array(mean), np.array(half)


class SessionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.refdir = need("FRANKENSIM_REFERENCE_DIR")
        lib = need("IRP_LIBRARY")
        import json
        cls.ch = Chain.from_references(manifest(), SESSION, cls.refdir, json.loads(lib.read_text()))

    def test_chain_floor_matches_measured_floor(self):
        """G2: the emulated floor's third-octave PSD matches the session's measured ambient PSD within 1 dB per band, 20 Hz to
        20 kHz, each band printed with the 95 % interval of the difference from the frames' scatter in both."""
        amb = next(s for s in self.ch.sources if s["role"] == "ambient")
        ambient = self.ch._audio[amb["path"]]
        fl, _ = self.ch.floor(20 * self.ch.fs, "test_chain_floor_matches_measured_floor")
        q, _ = chain.quantize(fl, self.ch.bits)  # what a graded file holds: the floor as recorded, on the 16-bit grid
        ma, ha = band_psd(ambient, self.ch.fs)
        me, he = band_psd(q / 32768.0, self.ch.fs)
        d = me - ma
        ci = np.sqrt(ha ** 2 + he ** 2)
        print(f"\n  session {self.ch.id} (floor_source {self.ch.floor_source}): ambient {len(ambient) / self.ch.fs:.2f} s, emulated 20 s "
              f"from {len(self.ch.stretches)} clean stretch(es)")
        bad = []
        for fc, a, e, dd, c in zip(BANDS, ma, me, d, ci):
            ok = abs(dd) <= TOL_DB
            print(f"  {fc:8.1f} Hz  ambient {a:7.2f}  emulated {e:7.2f}  diff {dd:+6.2f} dB  (95 % +-{c:.2f})  {'ok' if ok else 'MISMATCH'}")
            if not ok:
                bad.append(float(round(fc, 1)))
        print(f"  bands off by more than {TOL_DB} dB: {bad}")
        self.assertEqual(bad, [])

    def test_a_session_with_an_ambient_file_takes_its_floor_from_it_alone(self):
        """The ruling (chain.py docstring): where a page lists an ambient file, it is the floor and no gap is read; its receipt says
        so and carries no gap_floor_unvalidated flag."""
        roles = sorted({s["role"] for s in self.ch.sources})
        flags = [f["flag"] for f in self.ch.flags]
        print(f"\n  {self.ch.id}: sources {roles}, floor_source {self.ch.floor_source}, flags {flags}")
        self.assertEqual(roles, ["ambient"])
        self.assertEqual(self.ch.floor_source, "ambient")
        self.assertNotIn("gap_floor_unvalidated", flags)

    def test_the_resampler_is_convert_refs(self):
        """The chain's way back to 48 kHz is the conversion every reference went through: one of the session's graded files,
        decoded and resampled by chain.resample_pcm, equals convert_refs.py's WAV of it sample for sample."""
        wavdir = need("IRP_REFS_WAV48")
        import json
        src = next(f for f in json.loads(manifest().read_text())["files"] if Path(f["path"]).stem == "BbClar.ff.C4B4")
        x, bits, fs = chain.read_pcm(self.refdir, src)
        q, _ = chain.quantize(x, bits)
        mine = chain.resample_pcm(q, bits, fs, 48000)
        theirs = load(wavdir / (Path(src["path"]).stem.replace(" ", "") + ".wav"))
        n = min(len(mine), len(theirs))
        dmax = float(np.max(np.abs(mine[:n] - theirs[:n])))
        print(f"\n  {src['path']}: {len(mine)} vs {len(theirs)} samples, largest difference {dmax:.3e} ({chain.resampler_id()})")
        self.assertEqual(len(mine), len(theirs))
        self.assertEqual(dmax, 0.0)


if __name__ == "__main__":
    unittest.main()
