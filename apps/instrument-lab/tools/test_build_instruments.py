"""What the #instruments cards are allowed to claim: refusals, the within/outside verdict and the badge.

Run with: uv run --with numpy --with pillow python -m unittest discover -s apps/instrument-lab/tools -p 'test_build_instruments.py' -v
The pairs below are throwaway inputs to the generator, written to a temp directory; none of it reaches the page.
"""
import copy
import json
import math
import tempfile
import unittest
import wave
from pathlib import Path
from unittest import mock

import numpy as np
from PIL import Image

import build_instruments as B

REFS = json.loads(B.REFERENCES.read_text())
REF_ID = next(f["id"] for f in REFS["files"] if f["instrument"] == "Bb clarinet" and f["kind"] != "ambient-silence")


def write_wav(path: Path, peak: float, seconds: float = 1.0) -> None:
    t = np.arange(int(48000 * seconds)) / 48000
    x = (peak * np.sin(2 * math.pi * 220 * t) * 32767).astype("<i2")
    with wave.open(str(path), "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(48000)
        w.writeframes(x.tobytes())


def write_png(path: Path) -> None:
    Image.fromarray(np.full((40, 80, 3), 90, dtype=np.uint8)).save(path)


def fig(with_ours: bool) -> dict:
    return {"real": "s_real.png", "ours": "s_ours.png" if with_ours else None, "diff": "s_diff.png" if with_ours else None,
            "f_lo_hz": 40, "f_hi_hz": 16000, "span_db": 90, "diff_range_db": 24}


def track(n: int, frames: int = 4) -> dict:
    return {"n": n, "f_hz": [220.0 * n * (1 + 0.001 * k) for k in range(frames)], "level_db": [-20.0 - 6 * n] * frames}


def metric(id_: str, gap, ruler) -> dict:
    return {"id": id_, "label": id_, "unit": "dB", "real": None, "ours": None, "gap": gap, "ruler": [{"label": f"r{i}", "value": v} for i, v in enumerate(ruler)], "note": ""}


def make_pair(root: Path, with_ours: bool = True, ours_peak: float = 0.5) -> dict:
    d = root / "clarinet" / "D3"
    d.mkdir(parents=True)
    write_wav(d / "real.wav", 0.5)
    for name in ("s_real", "s_ours", "s_diff"):
        write_png(d / f"{name}.png")
    if with_ours:
        write_wav(d / "ours.wav", ours_peak)
    pj = {
        "schema": B.PAIR_SCHEMA, "family": "clarinet", "note": "D3", "dynamic": "mf", "duration_s": 1.0,
        "real": {"reference_id": REF_ID, "file": "BbClar.mf.C3B3.aiff", "segment_s": [0.0, 1.0]},
        "ours": {"wav": "/x/ours.wav", "git_sha": "abc123", "observer": "baffled_outlet 0 0 0.5", "sample_rate": 48000, "recipe": "data/instruments/descriptions/clarinet/render.sh"} if with_ours else None,
        "audio": {"real": "real.wav", "ours": "ours.wav" if with_ours else None, "level_match": "equal RMS over the first 250 ms"},
        "alignment": None,
        "figures": {"spectrogram": fig(with_ours), "scalogram": fig(with_ours)},
        "partials": {"t_s": [0, 1, 2, 3], "f0_hz": {"real": 220.0, "ours": 220.5 if with_ours else None},
                     "real": [track(1), track(2)], "ours": [track(1), track(2)] if with_ours else None},
        "metrics": [metric("a", 1.0, [2.0, 3.0]), metric("b", 5.0, [2.0, 3.0]), metric("c", 1.0, [])],
    }
    if not with_ours:
        for m in pj["metrics"]:
            m["gap"] = None
    (d / "pair.json").write_text(json.dumps(pj))
    return pj


class Verdicts(unittest.TestCase):
    def test_within_means_no_larger_than_the_largest_real_vs_real_distance(self):
        m = B.metric_view(metric("x", 3.0, [1.0, 3.0]), True)
        self.assertEqual((m["verdict"], m["ruler_max"], m["ratio"]), ("within", 3.0, 1.0))
        self.assertEqual(B.metric_view(metric("x", 3.01, [1.0, 3.0]), True)["verdict"], "outside")

    def test_no_ruler_makes_no_claim(self):
        self.assertEqual(B.metric_view(metric("x", 9.0, []), True)["verdict"], "no-ruler")

    def test_a_pair_without_a_render_has_no_gap_even_if_one_is_supplied(self):
        self.assertEqual(B.metric_view(metric("x", 9.0, [1.0]), False)["verdict"], "no-gap")


class Badges(unittest.TestCase):
    def pair(self, verdicts):
        ms = [B.metric_view(metric(f"m{i}", g, [2.0] if g is not None else []), True) for i, g in enumerate(verdicts)]
        return {"note": "D3", "ours": {"x": 1}, "metrics": ms}

    def test_no_pairs_says_no_render_yet(self):
        j = B.judge(None, [])
        self.assertEqual((j["badge"], j["label"]), ("none", "No render yet"))

    def test_blocked_family_says_why(self):
        j = B.judge({"blocked": "needs the tone-hole operator"}, [])
        self.assertEqual((j["badge"], j["line"]), ("blocked", "needs the tone-hole operator"))

    def test_reference_only_pairs_are_not_a_render(self):
        ref_only = {"note": "D3", "ours": None, "metrics": []}
        self.assertEqual(B.judge(None, [ref_only])["badge"], "none")

    def test_close_mixed_off(self):
        self.assertEqual(B.judge(None, [self.pair([1.0, 2.0])])["badge"], "close")
        self.assertEqual(B.judge(None, [self.pair([1.0, 5.0])])["badge"], "mixed")
        self.assertEqual(B.judge(None, [self.pair([5.0, 6.0, 1.0])])["badge"], "off")

    def test_rendered_but_nothing_to_compare_against_is_not_graded(self):
        self.assertEqual(B.judge(None, [self.pair([None])])["badge"], "ungraded")

    def test_the_sentence_names_the_worst_measure_by_its_computed_ratio(self):
        j = B.judge(None, [self.pair([1.0, 8.0])])
        self.assertIn("1 of 2", j["line"])
        self.assertIn("4.0×", j["line"])


class Refusals(unittest.TestCase):
    def check(self, mutate=None, **kw):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pj = make_pair(root, **kw)
            if mutate:
                mutate(pj)
            base = root / "clarinet" / "D3"
            by_id = {f["id"]: f for f in REFS["files"]}
            problems = B.validate_pair(pj, "clarinet", "D3", base, by_id)
            return problems or B.check_audio(base, pj)

    def test_a_good_pair_passes(self):
        self.assertEqual(self.check(), [])
        self.assertEqual(self.check(with_ours=False), [])

    def test_a_clip_at_full_scale_is_refused(self):
        self.assertTrue(any("full scale" in p for p in self.check(ours_peak=1.0)))

    def test_clip_length_must_match_duration(self):
        self.assertTrue(any("lasts" in p for p in self.check(lambda pj: pj.update(duration_s=2.0))))

    def test_unknown_reference_is_refused(self):
        self.assertTrue(any("references.json" in p for p in self.check(lambda pj: pj["real"].update(reference_id="nope"))))

    def test_a_ruler_value_must_be_a_number(self):
        self.assertTrue(any("ruler entry" in p for p in self.check(lambda pj: pj["metrics"][0]["ruler"].append({"label": "x", "value": None}))))

    def test_track_length_must_match_the_frame_times(self):
        self.assertTrue(any("must have 4 entries" in p for p in self.check(lambda pj: pj["partials"]["real"][0].update(f_hz=[1.0]))))

    def test_a_missing_render_cannot_carry_ours_files(self):
        self.assertTrue(any("must be null" in p for p in self.check(lambda pj: pj["audio"].update(ours="real.wav"), with_ours=False)))

    def test_family_must_match_its_directory(self):
        self.assertTrue(any("is not its directory" in p for p in self.check(lambda pj: pj.update(family="oboe"))))


class Partials(unittest.TestCase):
    def test_cents_are_relative_to_the_real_recordings_median_so_a_detuned_render_shows_up(self):
        real = {"n": 1, "f_hz": [100.0, 100.0, 100.0], "level_db": [-10, -10, -10]}
        ours = {"n": 1, "f_hz": [101.0, 101.0, None], "level_db": [-12, -12, None]}
        v = B.partials_view({"t_s": [0, 1, 2], "f0_hz": {"real": 100.0}, "real": [real], "ours": [ours]})
        t = v["tracks"][0]
        self.assertEqual(t["real"]["cents"], [0.0, 0.0, 0.0])
        self.assertAlmostEqual(t["ours"]["cents"][0], 1200 * math.log2(1.01), places=2)
        self.assertIsNone(t["ours"]["cents"][2])

    def test_only_the_first_eight_partials_are_kept(self):
        tracks = [track(n) for n in range(1, 13)]
        v = B.partials_view({"t_s": [0, 1, 2, 3], "f0_hz": {}, "real": tracks, "ours": None})
        self.assertEqual([t["n"] for t in v["tracks"]], list(range(1, 9)))


class EndToEnd(unittest.TestCase):
    def test_main_writes_cards_for_every_family_and_a_graded_clarinet(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "analysis", Path(tmp) / "out"
            make_pair(root)
            with mock.patch("sys.argv", ["build_instruments.py", "--analysis", str(root), "--out", str(out)]):
                self.assertEqual(B.main(), 0)
            index = json.loads((out / "index.json").read_text())
            self.assertEqual(len(index["families"]), len(B.FAMILIES) + len([p for p in B.DESCRIPTIONS.glob("*/status.json") if p.parent.name not in {f[0] for f in B.FAMILIES}]))
            by_id = {f["id"]: f for f in index["families"]}
            self.assertEqual(by_id["oboe"]["badge"], "none")
            self.assertEqual(by_id["clarinet"]["badge"], "mixed")
            card = json.loads((out / "clarinet" / "card.json").read_text())
            p = card["pairs"][0]
            for path in (p["audio"]["real"], p["audio"]["ours"], p["figures"]["spectrogram"]["diff"], p["figures"]["scalogram"]["real"]):
                self.assertTrue((out / path).is_file(), path)
            self.assertEqual(p["real"]["reference_id"], REF_ID)

    def test_a_refused_pair_is_listed_and_the_exit_code_says_so(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "analysis", Path(tmp) / "out"
            make_pair(root, ours_peak=1.0)
            with mock.patch("sys.argv", ["build_instruments.py", "--analysis", str(root), "--out", str(out)]):
                self.assertEqual(B.main(), 2)
            card = json.loads((out / "clarinet" / "card.json").read_text())
            self.assertEqual(card["pairs"], [])
            self.assertEqual(card["refused"][0]["note"], "D3")
            self.assertEqual(card["status"]["badge"], "none")


if __name__ == "__main__":
    unittest.main()
