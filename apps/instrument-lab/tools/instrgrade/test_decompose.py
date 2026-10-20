"""decompose.py on synthetic notes with known answers.

uv run --with numpy --with matplotlib python3 -I test_decompose.py -v
Every expected value is fixed by how the signal is built (a delay, a detuning, a missing partial), not by a recording.
"""
from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import decompose as dc
from decompose import SR, Note

F0 = 220.0


def tone(f0=F0, dur=2.6, attack=0.04, lead=0.3, delay=0.0, odd_only=False, cents=0.0, noise_db=None, tilt=1.0, seed=1):
    """Sustained harmonic tone with a raised-cosine attack, `lead` seconds of leading silence, optional noise floor."""
    n = int(dur * SR)
    t = np.arange(n) / SR
    f = f0 * 2 ** (cents / 1200)
    x = np.zeros(n)
    for k in range(1, 25):
        if odd_only and k % 2 == 0:
            continue
        if k * f > 15000:
            break
        x += np.sin(2 * np.pi * k * f * t + 0.3 * k) / k ** tilt
    env = np.clip((t - lead - delay) / attack, 0, 1)
    x *= 0.5 - 0.5 * np.cos(np.pi * env)
    x *= 0.2 / np.abs(x).max()
    if noise_db is not None:
        x += np.random.default_rng(seed).standard_normal(n) * 10 ** (noise_db / 20)
    return x


def note_of(x, f0=F0, kind="sustained"):
    floor = dc.floor_of(dc.highpass(x, f0))
    y, info = dc.align(x, f0, floor)
    return Note(y, f0, kind, False, "synthetic", info)


class AlignmentTests(unittest.TestCase):
    def test_delayed_copy_aligns_to_the_millisecond(self):
        a, b = note_of(tone()), note_of(tone(delay=0.137, noise_db=-80))
        self.assertLess(abs(a.info["cut_onset_s"] - b.info["cut_onset_s"]) * 1000, dc.ALIGN_TOL_MS)
        self.assertAlmostEqual(a.info["cut_onset_s"], dc.PREROLL_S, delta=0.0015)
        self.assertGreater(abs(a.info["onset_s"] - b.info["onset_s"]), 0.1)  # the cut really moved them

    def test_silence_refuses(self):
        with self.assertRaises(dc.OnsetMisaligned):
            dc.onset_times(np.zeros(SR), -120.0)


class GapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = note_of(tone())

    def gaps(self, x):
        return dc.compare(note_of(x), self.ref)

    def test_identical_note_has_zero_gaps(self):
        g = dc.compare(self.ref, self.ref)
        for m in ("spec_db", "stft_lsd_db", "cwt_lsd_db", "cepstral_env_db", "partial_freq_cents", "f0_cents", "attack_log2",
                  "centroid_st", "flux_db", "odd_even_db"):
            self.assertAlmostEqual(g[m][0], 0.0, places=6, msg=m)

    def test_level_and_delay_do_not_register(self):
        g = self.gaps(0.25 * tone(delay=0.2, noise_db=-90))
        self.assertLess(g["stft_lsd_db"][0], 2.0)
        self.assertLess(g["attack_log2"][0], 0.1)
        self.assertLess(g["spec_db"][0], 1.0)

    def test_detuning_is_read_in_cents(self):
        g = self.gaps(tone(cents=25.0))
        self.assertAlmostEqual(g["f0_cents"][1], 25.0, delta=1.0)
        self.assertLess(g["partial_freq_cents"][0], 3.0)  # a pure detune moves f0, not the partial ratios

    def test_attack_ratio_is_log2(self):
        g = self.gaps(tone(attack=0.16))  # 4x slower
        self.assertAlmostEqual(g["attack_log2"][0], 2.0, delta=0.35)

    def test_odd_only_spectrum_is_far_in_odd_even_and_harmonic_levels(self):
        g = self.gaps(tone(odd_only=True))
        self.assertGreater(g["odd_even_db"][0], 6.0)
        self.assertGreater(g["harm_db"][0], 3.0)

    def test_darker_tone_has_lower_centroid(self):
        g = self.gaps(tone(tilt=2.0))
        self.assertLess(g["centroid_st"][1], -1.0)
        self.assertGreater(g["cepstral_env_db"][0], 1.0)

    def test_stretched_partials_show_in_partial_frequency(self):
        n = int(2.6 * SR)
        t = np.arange(n) / SR
        B = 4e-4
        x = sum(np.sin(2 * np.pi * k * F0 * math.sqrt(1 + B * k * k) * t) / k for k in range(1, 20))
        env = np.clip((t - 0.3) / 0.04, 0, 1)
        x = 0.2 * x * (0.5 - 0.5 * np.cos(np.pi * env)) / np.abs(x).max()
        g = self.gaps(x)
        self.assertGreater(g["inharmonicity_B"][0], 1e-4)
        self.assertGreater(g["partial_freq_cents"][0], 3.0)


class ContractTests(unittest.TestCase):
    """write_pair against P8's validator (tools/build_instruments.py) on a synthetic library."""

    @classmethod
    def setUpClass(cls):
        import json
        cls.tmp = Path(tempfile.mkdtemp())
        wav = cls.tmp / "wav"
        wav.mkdir()
        # a UIowa-style file: three notes, each preceded by silence, plus a quieter second dynamic of the middle note
        notes = [(57, F0), (58, F0 * 2 ** (1 / 12)), (59, F0 * 2 ** (2 / 12))]

        def scale(scale_gain):
            out = [np.zeros(SR)]
            for _, f in notes:
                out += [tone(f0=f, dur=2.0, lead=0.0, attack=0.04) * scale_gain, np.zeros(SR)]
            return np.concatenate(out)
        from instrgrade import write_wav
        write_wav(wav / "Syn.mf.A3B3.wav", scale(1.0))
        write_wav(wav / "Syn.ff.A3B3.wav", 2.0 * scale(0.5))
        lib = []
        for fname, dyn, g in (("Syn.mf.A3B3.wav", "mf", 1.0), ("Syn.ff.A3B3.wav", "ff", 1.0)):
            for k, (m, f) in enumerate(notes):
                t0 = 1.0 + 3.0 * k
                lib.append(dict(file=fname, start=t0 - 0.005, end=t0 + 2.2, note=m, name=dc.name_of(m), kind="sustained", inst="syn", art="",
                                dyn=dyn, sul="", mallet="", range="A3B3", desc={}))
        (cls.tmp / "lib.json").write_text(json.dumps(lib))
        (cls.tmp / "refs.json").write_text(json.dumps({"files": [
            dict(path="x/Syn.mf.A3B3.aiff", page="p.html", chain="c", id="syn-mf"), dict(path="x/Syn.ff.A3B3.aiff", page="p.html", chain="c", id="syn-ff")]}))
        cls.refs = dc.Refs(cls.tmp / "lib.json", wav, cls.tmp / "refs.json")
        cls.render = cls.tmp / "render.wav"
        write_wav(cls.render, 0.3 * tone(f0=F0 * 2 ** (1 / 12), dur=2.4, lead=0.2, noise_db=-85, cents=8))

    def validate(self, root, family, note):
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
        import build_instruments as bi
        import json
        base = Path(root) / family / note
        pj = json.loads((base / "pair.json").read_text())
        problems = bi.validate_pair(pj, family, note, base, {"syn-mf": {}, "syn-ff": {}}) + bi.check_audio(base, pj)
        self.assertEqual(problems, [])
        return pj

    def test_reference_only_pair_meets_the_contract(self):
        root = self.tmp / "ref"
        ref_e = self.refs.find("syn", 58, "mf")[0]
        y, info = self.refs.aligned(ref_e)
        pj = dc.write_pair(root, "syn", "Bb3", ref_e, self.refs, y, info)
        self.validate(root, "syn", "Bb3")
        self.assertIsNone(pj["ours"])
        self.assertTrue(any(m["ruler"] for m in pj["metrics"]))
        self.assertTrue(all(m["gap"] is None for m in pj["metrics"]))

    def test_render_pair_meets_the_contract_and_reads_the_detune(self):
        root = self.tmp / "pair"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="abc", observer="test", recipe="r.sh")
        pj = dc.run_pair(self.refs, c, root)
        self.validate(root, "syn", "Bb3")
        by = {m["id"]: m for m in pj["metrics"]}
        self.assertAlmostEqual(by["f0_cents"]["gap"], 8.0, delta=1.5)
        self.assertAlmostEqual(pj["alignment"]["level_offset_db"], 20 * math.log10(1 / 0.3), delta=1.0)
        self.assertEqual(by["f0_cents"]["gap_unit"], "cents")
        self.assertNotIn("neighbour", " ".join(r["label"] for r in by["f0_cents"]["ruler"]))  # pitch-absolute: same note only
        self.assertTrue(any("neighbouring" in r["label"] for r in by["attack_log2"]["ruler"]))

    def test_render_pair_refuses_misaligned_onsets(self):
        ref_e = self.refs.find("syn", 58, "mf")[0]
        y, info = self.refs.aligned(ref_e)
        raw = tone(f0=F0 * 2 ** (1 / 12), dur=2.4, lead=0.2)
        cy, ci = dc.align(raw, F0 * 2 ** (1 / 12), dc.floor_of(raw))
        ci = dict(ci, cut_onset_s=ci["cut_onset_s"] + 0.01)  # the render's onset sits 10 ms late in its own cut
        meta = dict(wav="r.wav", git_sha="x", observer="o", sample_rate=48000, recipe="r.sh")
        with self.assertRaises(dc.OnsetMisaligned):
            dc.write_pair(self.tmp / "bad", "syn", "Bb3", ref_e, self.refs, y, info, cy, ci, meta)

    def test_json_has_no_nan(self):
        import json
        s = json.dumps(dc.clean({"a": float("nan"), "b": np.float32(1.5), "c": [np.inf]}), allow_nan=False)
        self.assertEqual(json.loads(s), {"a": None, "b": 1.5, "c": [None]})


if __name__ == "__main__":
    unittest.main()
