"""What the #instruments cards are allowed to claim: refusals, the within/outside verdict and the badge.

Run with: uv run --with numpy --with pillow python -m unittest discover -s apps/instrument-lab/tools -p 'test_build_instruments.py' -v
The pairs below are throwaway inputs to the generator, written to a temp directory; none of it reaches the page.
"""
import collections
import contextlib
import copy
import io
import json
import math
import shutil
import subprocess
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


def metric(id_: str, gap, ruler, strict=None) -> dict:
    """strict None: the pair predates the strict ruler. Otherwise, like decompose.py: ruler_strict plus both inside_* verdicts."""
    m = {"id": id_, "label": id_, "unit": "dB", "real": None, "ours": None, "gap": gap, "ruler": [{"label": f"r{i}", "value": v} for i, v in enumerate(ruler)], "note": ""}
    if strict is not None:
        m.update(ruler_strict=[{"label": f"s{i}", "value": v} for i, v in enumerate(strict)], inside_loose=B.inside(gap, ruler), inside_strict=B.inside(gap, strict))
    return m


# A hand-made set of measures, one of each case the strict ruler exists to tell apart (gap, loose ruler, strict ruler):
STRICT_CASES = [
    metric("both_inside", 1.0, [2.0, 4.0], [1.5]),        # inside both spreads
    metric("loose_only", 3.0, [2.0, 4.0], [1.0, 2.0]),    # inside the neighbouring-note spread, OUTSIDE the same-note one
    metric("both_outside", 9.0, [2.0, 4.0], [1.5]),       # outside both
    metric("no_same_note", 2.0, [3.0], []),               # no other recording of this note: no strict claim
    metric("no_claim", 1.0, [], []),                      # no ruler of either kind
]


def make_pair(root: Path, with_ours: bool = True, ours_peak: float = 0.5, metrics: list | None = None) -> dict:
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
        "metrics": metrics if metrics is not None else [metric("a", 1.0, [2.0, 3.0]), metric("b", 5.0, [2.0, 3.0]), metric("c", 1.0, [])],
    }
    if not with_ours:
        for m in pj["metrics"]:
            m["gap"] = None
    (d / "pair.json").write_text(json.dumps(pj))
    return pj


class StrictRuler(unittest.TestCase):
    """The two rulers, from a hand-made pair.json (STRICT_CASES): loose = this and neighbouring notes, strict = the same note only."""

    def views(self) -> dict:
        return {v["id"]: v for v in (B.metric_view(m, True) for m in STRICT_CASES)}

    def test_each_measure_carries_both_verdicts_and_a_loose_only_measure_is_told_apart(self):
        v = self.views()
        self.assertEqual((v["both_inside"]["inside_loose"], v["both_inside"]["inside_strict"]), (True, True))
        self.assertEqual((v["loose_only"]["inside_loose"], v["loose_only"]["inside_strict"]), (True, False))
        self.assertEqual((v["loose_only"]["verdict"], v["loose_only"]["verdict_strict"], v["loose_only"]["ratio_strict"]), ("within", "outside", 1.5))
        self.assertEqual((v["both_outside"]["verdict"], v["both_outside"]["verdict_strict"]), ("outside", "outside"))
        self.assertEqual((v["no_same_note"]["inside_loose"], v["no_same_note"]["inside_strict"], v["no_same_note"]["verdict_strict"]), (True, None, "no-ruler"))
        self.assertEqual((v["no_claim"]["inside_loose"], v["no_claim"]["inside_strict"]), (None, None))
        self.assertEqual((v["loose_only"]["ruler_strict_max"], v["loose_only"]["ruler_strict_min"], v["loose_only"]["ruler_max"]), (2.0, 1.0, 4.0))

    def test_a_pair_without_the_strict_fields_shows_the_loose_verdict_only(self):
        v = B.metric_view(metric("old", 3.0, [2.0, 4.0]), True)
        self.assertEqual((v["verdict"], v["ruler_strict"], v["verdict_strict"], v["inside_strict"]), ("within", None, None, None))

    def test_the_pair_counts_both_verdicts(self):
        want = {"measures": 5, "loose": {"inside": 3, "outside": 1, "no_ruler": 1}, "strict": {"inside": 1, "outside": 2, "no_ruler": 2}}
        self.assertEqual(B.pair_verdicts(list(self.views().values())), want)
        self.assertIsNone(B.pair_verdicts([B.metric_view(metric("old", 3.0, [2.0]), True)])["strict"])

    def pair_dict(self, metrics: list, verdicts=None) -> dict:
        return {"note": "D3", "real": {}, "ours": {"x": 1}, "metrics": [B.metric_view(m, True) for m in metrics], "verdicts": verdicts}

    def test_the_badge_is_the_weaker_of_the_two_and_the_sentence_gives_both_counts(self):
        # Loose says nothing is outside (close); strict has one measure outside (1 of 2 inside: mixed): the badge is mixed.
        j = B.judge(None, [self.pair_dict([metric("a", 1.0, [3.0], [2.0]), metric("b", 2.5, [3.0], [2.0]), metric("c", 1.0, [3.0], [])])])
        self.assertEqual(j["badge"], "mixed")
        self.assertEqual((j["counts"]["within"], j["counts"]["outside"], j["counts"]["strict_within"], j["counts"]["strict_outside"], j["counts"]["strict_unruled"]), (3, 0, 1, 1, 1))
        self.assertIn("inside the spread of real recordings of this and neighbouring notes, 3 of 3 measures", j["line"])
        self.assertIn("Inside the spread of real recordings of the same note, 1 of 2 measures (1 has no same-note comparison)", j["line"])
        self.assertIn("The badge is the weaker of the two.", j["line"])
        self.assertIn("Largest gap against the same-note spread: b, 1.2× it (D3)", j["line"])
        # Strict inside everywhere it can speak, loose too: close. Strict outside nearly everywhere: off, even though loose is mixed.
        self.assertEqual(B.judge(None, [self.pair_dict([metric("a", 1.0, [3.0], [2.0])])])["badge"], "close")
        off = B.judge(None, [self.pair_dict([metric("a", 1.0, [3.0], [0.5]), metric("b", 5.0, [9.0], [0.5]), metric("c", 2.0, [3.0], [0.5]), metric("d", 9.5, [9.0], [0.5])])])
        self.assertEqual((off["counts"]["within"], off["counts"]["outside"], off["badge"]), (3, 1, "off"))

    def test_when_no_note_has_a_same_note_recording_the_badge_is_the_loose_one_and_says_so(self):
        j = B.judge(None, [self.pair_dict([metric("a", 1.0, [3.0], []), metric("b", 5.0, [3.0], [])])])
        self.assertEqual((j["badge"], j["counts"]["strict_notes"]), ("mixed", 1))
        self.assertIn("No same-note comparison could be made", j["line"])
        self.assertNotIn("The badge is the weaker", j["line"])

    def test_a_family_of_older_pairs_says_the_same_note_comparison_is_not_computed(self):
        j = B.judge(None, [self.pair_dict([metric("a", 1.0, [3.0]), metric("b", 5.0, [3.0])])])
        self.assertEqual((j["badge"], j["counts"]["strict_notes"]), ("mixed", 0))
        self.assertIn("The same-note comparison has not been computed for these notes.", j["line"])

    def test_a_pair_whose_own_verdicts_disagree_with_its_numbers_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pj = make_pair(root, metrics=[dict(m) for m in STRICT_CASES])
            base = root / "clarinet" / "D3"
            check = lambda p: B.validate_pair(p, "clarinet", "D3", base, {REF_ID: next(f for f in REFS["files"] if f["id"] == REF_ID)})  # noqa: E731
            pj["verdicts"] = B.pair_verdicts([B.metric_view(m, True) for m in pj["metrics"]])
            self.assertEqual(check(pj), [], "the pair decompose.py would write is accepted")
            wrong = copy.deepcopy(pj)
            wrong["metrics"][1]["inside_strict"] = True  # the loose-only measure claimed inside the same-note spread
            self.assertTrue(any("inside_strict" in p and "loose_only" in p for p in check(wrong)), check(wrong))
            wrong = copy.deepcopy(pj)
            wrong["verdicts"]["strict"]["inside"] = 4
            self.assertTrue(any("verdicts.strict" in p for p in check(wrong)), check(wrong))
            wrong = copy.deepcopy(pj)
            wrong["metrics"][0]["ruler_strict"] = [{"label": "s0", "value": -1.0}]
            self.assertTrue(any("ruler_strict" in p for p in check(wrong)), check(wrong))
            old = copy.deepcopy(pj)
            for m in old["metrics"]:
                for key in ("ruler_strict", "inside_loose", "inside_strict"):
                    m.pop(key)
            old["verdicts"] = None
            self.assertEqual(check(old), [], "a pair that predates the strict ruler is still accepted")

    def test_the_card_carries_both_rulers_and_the_pairs_verdicts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "analysis", Path(tmp) / "out"
            pj = make_pair(root, metrics=[dict(m) for m in STRICT_CASES])
            pj["verdicts"] = B.pair_verdicts([B.metric_view(m, True) for m in pj["metrics"]])
            (root / "clarinet" / "D3" / "pair.json").write_text(json.dumps(pj))
            with mock.patch("sys.argv", ["build_instruments.py", "--analysis", str(root), "--out", str(out)]):
                self.assertEqual(B.main(), 0)
            card = json.loads((out / "clarinet" / "card.json").read_text())
            p = card["pairs"][0]
            self.assertEqual(p["verdicts"], {"measures": 5, "loose": {"inside": 3, "outside": 1, "no_ruler": 1}, "strict": {"inside": 1, "outside": 2, "no_ruler": 2}})
            self.assertEqual([m["inside_strict"] for m in p["metrics"]], [True, False, False, None, None])
            self.assertEqual(card["status"]["badge"], "off", "strict: 1 of 3 inside, the weaker of the two")
            self.assertEqual((card["status"]["counts"]["strict_within"], card["status"]["counts"]["strict_outside"], card["status"]["counts"]["strict_unruled"]), (1, 2, 2))


class Verdicts(unittest.TestCase):
    def test_within_means_no_larger_than_the_largest_real_vs_real_distance(self):
        m = B.metric_view(metric("x", 3.0, [1.0, 3.0]), True)
        self.assertEqual((m["verdict"], m["ruler_max"], m["ratio"]), ("within", 3.0, 1.0))
        self.assertEqual(B.metric_view(metric("x", 3.01, [1.0, 3.0]), True)["verdict"], "outside")

    def test_tiny_values_keep_their_digits(self):
        # Real clarinet inharmonicity rulers are about 3e-6 to 9e-6: rounding to 4 decimals published them as 0.
        m = B.metric_view(metric("B", 4.0e-6, [9.2181e-06, 3.1258e-06]), True)
        self.assertEqual((m["ruler_min"], m["ruler_max"], m["gap"]), (3.1258e-06, 9.2181e-06, 4.0e-06))
        self.assertEqual(m["verdict"], "within")

    def test_no_ruler_makes_no_claim(self):
        self.assertEqual(B.metric_view(metric("x", 9.0, []), True)["verdict"], "no-ruler")

    def test_a_pair_without_a_render_has_no_gap_even_if_one_is_supplied(self):
        self.assertEqual(B.metric_view(metric("x", 9.0, [1.0]), False)["verdict"], "no-gap")


class Badges(unittest.TestCase):
    def pair(self, verdicts):
        ms = [B.metric_view(metric(f"m{i}", g, [2.0] if g is not None else []), True) for i, g in enumerate(verdicts)]
        return {"note": "D3", "real": {}, "ours": {"x": 1}, "metrics": ms}

    def test_no_pairs_says_no_render_yet(self):
        j = B.judge(None, [])
        self.assertEqual((j["badge"], j["label"]), ("none", "No render yet"))

    def test_blocked_family_says_why(self):
        j = B.judge({"blocked": "needs the tone-hole operator"}, [])
        self.assertEqual((j["badge"], j["line"]), ("blocked", "needs the tone-hole operator"))

    def test_reference_only_pairs_are_not_a_render(self):
        ref_only = {"note": "D3", "real": {}, "ours": None, "metrics": []}
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

    def test_a_render_none_of_whose_partials_matched_is_still_shown(self):
        self.assertEqual(self.check(lambda pj: pj["partials"].update(ours=None)), [])

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


def write_refusal(root: Path, family: str, note: str, **extra) -> None:
    d = root / family / note
    d.mkdir(parents=True, exist_ok=True)
    (d / "pair.json").write_text(json.dumps({"schema": B.PAIR_SCHEMA, "ours": None, "refused": {"reason": "reference measures -113 cents from the named note", "render_sha": "99f76b610", "time": "2026-10-20T23:40:00Z"}, **extra}))


VOCAL_SOURCE = {"id": "vocalset-m10-long-forte-a", "name": "VocalSet", "credit": "TEST CREDIT: Wilkins et al. 2018, Zenodo 1193957", "licence": "CC BY 4.0",
                "licence_url": "https://creativecommons.org/licenses/by/4.0/", "url": "https://zenodo.org/record/1193957"}


class Credits(unittest.TestCase):
    def check(self, mutate=None):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pj = make_pair(root)
            if mutate:
                mutate(pj)
            return B.validate_pair(pj, "clarinet", "D3", root / "clarinet" / "D3", {f["id"]: f for f in REFS["files"]})

    def external(self, pj, **over):
        pj["real"] = {"reference_id": "vocalset-m10-long-forte-a", "file": "m10_long_forte_a.wav", "segment_s": None, "source": {**VOCAL_SOURCE, **over}}

    def test_a_library_recording_is_credited_to_the_electronic_music_studios(self):
        ref = next(f for f in REFS["files"] if f["id"] == REF_ID)
        lic = {l["id"]: l for l in REFS["licences"]}
        s = B.real_source({"reference_id": REF_ID}, ref, lic)
        self.assertEqual(s["credit"], "University of Iowa Electronic Music Studios")
        self.assertFalse(s["external"])
        self.assertTrue(s["licence_url"].startswith("https://theremin.music.uiowa.edu"))
        self.assertIn("without restrictions", s["licence"])

    def test_every_library_licence_has_a_credit(self):
        self.assertTrue({f["licence"] for f in REFS["files"]} <= set(B.CREDITS))

    def test_a_library_entry_without_a_licence_is_a_problem_never_an_exception(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pj = make_pair(root)
            base = root / "clarinet" / "D3"
            for entry in ({"id": REF_ID, "url": "https://x/y.aiff"}, {"id": REF_ID, "licence": None, "url": "https://x/y.aiff"}, {"id": REF_ID, "licence": "nobody-knows", "url": "https://x/y.aiff"}):
                problems = B.validate_pair(pj, "clarinet", "D3", base, {REF_ID: entry})
                self.assertTrue(any("no credit line in CREDITS" in p for p in problems), entry)

    def test_a_library_entry_without_a_url_is_accepted_and_has_no_link(self):
        # P2's contract tests use synthetic references with a licence and no url: still a valid, creditable recording.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pj = make_pair(root)
            base = root / "clarinet" / "D3"
            entry = {"licence": "uiowa-mis"}
            self.assertEqual(B.validate_pair(pj, "clarinet", "D3", base, {REF_ID: entry}), [])
            pair = B.build_pair(pj, base, root / "out" / "clarinet" / "D3", {REF_ID: entry}, {}, {l["id"]: l for l in REFS["licences"]})
            self.assertIsNone(pair["real"]["url"])
            self.assertEqual(pair["real"]["source"]["credit"], "University of Iowa Electronic Music Studios")

    def test_a_minimal_library_entry_that_validates_also_builds(self):
        # P2's contract tests give the generator synthetic references with only an id, a licence and a url.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pj = make_pair(root)
            base = root / "clarinet" / "D3"
            entry = {"id": REF_ID, "licence": "uiowa-mis", "url": "https://x/y.aiff"}
            self.assertEqual(B.validate_pair(pj, "clarinet", "D3", base, {REF_ID: entry}), [])
            pair = B.build_pair(pj, base, root / "out" / "clarinet" / "D3", {REF_ID: entry}, {}, {l["id"]: l for l in REFS["licences"]})
            self.assertEqual(pair["real"]["source"]["credit"], "University of Iowa Electronic Music Studios")
            self.assertIsNone(pair["real"]["sha256"])

    def test_a_recording_outside_the_library_is_accepted_with_its_own_source(self):
        self.assertEqual(self.check(lambda pj: self.external(pj)), [])

    def test_attribution_is_mandatory_for_an_outside_recording(self):
        for key in ("credit", "licence", "name", "url", "id"):
            self.assertTrue(any(f"real.source.{key}" in p for p in self.check(lambda pj: self.external(pj, **{key: "  "}))), key)

    def test_an_unknown_recording_with_no_source_is_refused(self):
        self.assertTrue(any("real.source" in p for p in self.check(lambda pj: pj["real"].update(reference_id="nope"))))

    def test_the_card_carries_the_source_and_never_claims_a_library_reference_for_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out, desc = Path(tmp) / "analysis", Path(tmp) / "out", Path(tmp) / "desc"
            pj = make_pair(root, with_ours=True)
            pj["family"] = "voice"
            self.external(pj)
            d = root / "clarinet" / "D3"
            (root / "voice").mkdir()
            d.rename(root / "voice" / "D3")
            (root / "voice" / "D3" / "pair.json").write_text(json.dumps({**pj, "note": "D3"}))
            (desc / "voice").mkdir(parents=True)
            (desc / "voice" / "status.json").write_text(json.dumps({"schema": B.STATUS_SCHEMA, "family": "voice", "title": "Voice", "reference_instrument": "adult male voice (not in the library)", "summary": "", "physics_present": [], "physics_missing": [], "blocked": None, "beads": []}))
            with mock.patch("sys.argv", ["build_instruments.py", "--analysis", str(root), "--out", str(out), "--descriptions", str(desc)]):
                self.assertEqual(B.main(), 0)
            card = json.loads((out / "voice" / "card.json").read_text())
            self.assertIsNone(card["reference"], "a reference that is not in the library must not be described as library recordings")
            src = card["pairs"][0]["real"]["source"]
            self.assertEqual((src["credit"], src["external"], src["licence"]), (VOCAL_SOURCE["credit"], True, "CC BY 4.0"))
            index = json.loads((out / "index.json").read_text())
            self.assertIsNone(next(f for f in index["families"] if f["id"] == "voice")["reference"])


def make_listening(root: Path, family: str = "pan-pipes", note: str = "C5") -> dict:
    """A render with no reference recording, in the shape P2 and P8 agreed: real null, no_reference.reason, our audio and figures only."""
    d = root / family / note
    d.mkdir(parents=True)
    write_wav(d / "ours.wav", 0.5)
    for name in ("s_ours", "c_ours"):
        write_png(d / f"{name}.png")
    one = {"real": None, "ours": None, "diff": None, "f_lo_hz": 40, "f_hi_hz": 16000, "span_db": 90}
    pj = {
        "schema": B.PAIR_SCHEMA, "family": family, "note": note, "dynamic": "mf", "duration_s": 1.0,
        "real": None, "no_reference": {"reason": "No licensed recording of this instrument is in the reference library."},
        "ours": {"wav": "/x/ours.wav", "git_sha": "abc123", "observer": "baffled_outlet 0 0 0.5", "sample_rate": 48000, "recipe": "data/instruments/descriptions/pan-pipes/render.sh"},
        "audio": {"real": None, "ours": "ours.wav", "level_match": "peak at -6 dBFS; there is no recording to match"},
        "alignment": None,
        "figures": {"spectrogram": dict(one, ours="s_ours.png"), "scalogram": dict(one, ours="c_ours.png")},
        "partials": None, "metrics": [],
    }
    (d / "pair.json").write_text(json.dumps(pj))
    return pj


class Listening(unittest.TestCase):
    def check(self, pj: dict, root: Path, family: str = "pan-pipes", note: str = "C5") -> list[str]:
        return B.validate_pair(pj, family, note, root / family / note, {})

    def test_a_render_with_no_reference_is_valid_and_has_nothing_to_grade(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pj = make_listening(root)
            self.assertEqual(self.check(pj, root), [])
            self.assertEqual(B.check_audio(root / "pan-pipes" / "C5", pj), [])

    def test_anything_that_needs_a_recording_is_refused_on_such_a_note(self):
        cases = {
            "audio.real must be null": lambda p: p["audio"].update(real="real.wav"),
            "figures.spectrogram.real must be null": lambda p: p["figures"]["spectrogram"].update(real="s_ours.png"),
            "figures.scalogram.diff must be null": lambda p: p["figures"]["scalogram"].update(diff="c_ours.png"),
            "metrics must be []": lambda p: p.update(metrics=[metric("a", 1.0, [2.0])]),
            "partials must be null": lambda p: p.update(partials={"t_s": [0, 1], "f0_hz": {}, "real": None, "ours": None}),
            "alignment must be null": lambda p: p.update(alignment={"onset_real_s": 0.1, "onset_ours_s": 0.1, "highpass_hz": 50, "level_offset_db": 0}),
            "no_reference must be {reason}": lambda p: p.update(no_reference={"reason": "  "}),
            "needs a render": lambda p: p.update(ours=None),
            "figures.spectrogram.ours names no file": lambda p: p["figures"]["spectrogram"].update(ours="missing.png"),
        }
        for expect, mutate in cases.items():
            with self.subTest(expect), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                pj = make_listening(root)
                mutate(pj)
                self.assertTrue(any(expect in p for p in self.check(pj, root)), self.check(pj, root))

    def test_real_null_without_no_reference_is_not_a_listening_note(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pj = make_listening(root)
            del pj["no_reference"]
            self.assertTrue(any("real is missing" in p for p in self.check(pj, root)))

    def test_the_card_offers_our_audio_and_pictures_only_and_the_badge_says_there_is_nothing_to_grade_against(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "analysis", Path(tmp) / "out"
            make_listening(root)
            with mock.patch("sys.argv", ["build_instruments.py", "--analysis", str(root), "--out", str(out), "--descriptions", str(root / "none"), "--audio-ext", "mp4"]):
                self.assertEqual(B.main(), 0)
            card = json.loads((out / "pan-pipes" / "card.json").read_text())
            p = card["pairs"][0]
            self.assertEqual((card["status"]["badge"], card["status"]["label"]), ("listening", "Rendered, no reference to grade against"))
            self.assertEqual((p["real"], p["audio"]["real"], p["metrics"], p["partials"]), (None, None, [], None))
            self.assertEqual(p["no_reference"], "No licensed recording of this instrument is in the reference library.")
            self.assertEqual((card["status"]["counts"]["within"], card["status"]["counts"]["outside"], card["status"]["counts"]["notes_no_reference"]), (0, 0, 1))
            self.assertEqual(sorted(f.name for f in (out / "pan-pipes" / "C5").iterdir()), ["ours.mp4", "scalogram-ours.webp", "spectrogram-ours.webp"])
            self.assertEqual(p["figures"]["spectrogram"]["real"], None)
            index = json.loads((out / "index.json").read_text())
            self.assertEqual(next(f for f in index["families"] if f["id"] == "pan-pipes")["badge"], "listening")

    def test_next_to_graded_notes_it_does_not_change_the_families_badge(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "analysis", Path(tmp) / "out"
            make_pair(root)
            make_listening(root, "clarinet", "E3")
            with mock.patch("sys.argv", ["build_instruments.py", "--analysis", str(root), "--out", str(out), "--descriptions", str(root / "none")]):
                self.assertEqual(B.main(), 0)
            card = json.loads((out / "clarinet" / "card.json").read_text())
            self.assertEqual(card["status"]["badge"], "mixed")
            self.assertEqual([p["note"] for p in card["pairs"]], ["D3", "E3"])
            self.assertEqual((card["status"]["counts"]["within"], card["status"]["counts"]["outside"]), (1, 1), "the graded note alone decides the counts")


class GraderRefusals(unittest.TestCase):
    def run_main(self, root: Path, out: Path) -> int:
        with mock.patch("sys.argv", ["build_instruments.py", "--analysis", str(root), "--out", str(out), "--descriptions", str(root / "none")]):
            return B.main()

    def test_a_refusal_record_is_listed_with_its_reason_and_counts_in_no_badge(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "analysis", Path(tmp) / "out"
            make_pair(root)
            write_refusal(root, "clarinet", "E5")
            self.assertEqual(self.run_main(root, out), 0, "a refusal by the grader is data, not a build failure")
            card = json.loads((out / "clarinet" / "card.json").read_text())
            self.assertEqual([p["note"] for p in card["pairs"]], ["D3"])
            (r,) = card["refused"]
            self.assertEqual((r["note"], r["kind"], r["render_sha"]), ("E5", "grader", "99f76b610"))
            self.assertIn("-113 cents", r["reason"])
            self.assertEqual(card["status"]["counts"]["notes_refused_by_grader"], 1)
            with tempfile.TemporaryDirectory() as tmp2:
                root2, out2 = Path(tmp2) / "analysis", Path(tmp2) / "out"
                make_pair(root2)
                self.run_main(root2, out2)
                without = json.loads((out2 / "clarinet" / "card.json").read_text())["status"]
            keep = lambda s: {k: s[k] for k in ("badge", "label", "line")} | {k: s["counts"][k] for k in ("within", "outside", "unruled", "notes_with_render")}  # noqa: E731
            self.assertEqual(keep(card["status"]), keep(without), "the refused note must not move the badge, the sentence or any count")

    def test_a_family_whose_only_render_was_refused_reads_rendered_not_graded(self):
        j = B.judge(None, [], [{"note": "A1", "reason": "reference is +54 cents off", "render_sha": "945988669d92"}])
        self.assertEqual(j["badge"], "ungraded")
        self.assertIn("reference is +54 cents off", j["line"])
        self.assertIn("nothing is graded", j["line"])

    def test_a_refused_note_without_a_render_does_not_claim_one(self):
        j = B.judge(None, [], [{"note": "A1", "reason": "no signal", "render_sha": None}])
        self.assertEqual(j["badge"], "none")

    def test_a_malformed_refusal_is_a_site_refusal_and_fails_the_build(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "analysis", Path(tmp) / "out"
            d = root / "clarinet" / "E5"
            d.mkdir(parents=True)
            (d / "pair.json").write_text(json.dumps({"refused": {"reason": ""}, "ours": None}))
            self.assertEqual(self.run_main(root, out), 2)
            card = json.loads((out / "clarinet" / "card.json").read_text())
            self.assertEqual(card["refused"][0]["kind"], "site")

    def test_a_refusal_in_the_wrong_directory_is_refused(self):
        with self.assertRaises(B.Refused):
            B.grader_refusal({"family": "oboe", "refused": {"reason": "x"}}, "clarinet", "E5")


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
            # No oboe has been rendered here, so its badge is whatever its own status card says: blocked once a family lane
            # writes a `blocked` reason, none while it has no status at all.
            oboe = B.DESCRIPTIONS / "oboe" / "status.json"
            blocked = oboe.is_file() and bool(json.loads(oboe.read_text()).get("blocked"))
            self.assertEqual(by_id["oboe"]["badge"], "blocked" if blocked else "none")
            self.assertEqual(by_id["clarinet"]["badge"], "mixed")
            card = json.loads((out / "clarinet" / "card.json").read_text())
            p = card["pairs"][0]
            for path in (p["audio"]["real"], p["audio"]["ours"], p["figures"]["spectrogram"]["diff"], p["figures"]["scalogram"]["real"]):
                self.assertTrue((out / path).is_file(), path)
            self.assertEqual(p["real"]["reference_id"], REF_ID)

    def test_audio_ext_names_the_clips_and_the_card_agrees_with_the_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "analysis"
            make_pair(root)
            for ext, argv in (("m4a", []), ("mp4", ["--audio-ext", "mp4"])):
                out = Path(tmp) / f"out-{ext}"
                with mock.patch("sys.argv", ["build_instruments.py", "--analysis", str(root), "--out", str(out), *argv]):
                    self.assertEqual(B.main(), 0)
                p = json.loads((out / "clarinet" / "card.json").read_text())["pairs"][0]
                for side in ("real", "ours"):
                    self.assertRegex(p["audio"][side], rf"/{side}\.{ext}$")
                    self.assertTrue((out / p["audio"][side]).is_file())
                    codec = subprocess.run(["ffprobe", "-v", "error", "-select_streams", "a:0", "-show_entries", "stream=codec_name", "-of", "csv=p=0",
                                            str(out / p["audio"][side])], capture_output=True, text=True).stdout.strip()
                    self.assertEqual(codec, "aac")
                other = "mp4" if ext == "m4a" else "m4a"
                self.assertEqual(list(out.rglob(f"*.{other}")), [], f"no .{other} clip beside the .{ext} ones")

    def test_only_shows_the_named_notes_and_lists_the_rest_as_withheld_with_no_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "analysis"
            make_pair(root)
            shutil.copytree(root / "clarinet" / "D3", root / "clarinet" / "E3")
            e3 = root / "clarinet" / "E3" / "pair.json"
            e3.write_text(json.dumps(dict(json.loads(e3.read_text()), note="E3")))
            everything, some = Path(tmp) / "all", Path(tmp) / "some"
            with mock.patch("sys.argv", ["build_instruments.py", "--analysis", str(root), "--out", str(everything)]):
                self.assertEqual(B.main(), 0)
            card = json.loads((everything / "clarinet" / "card.json").read_text())
            self.assertEqual(([p["note"] for p in card["pairs"]], card["withheld"]), (["D3", "E3"], []))
            err = io.StringIO()
            with mock.patch("sys.argv", ["build_instruments.py", "--analysis", str(root), "--out", str(some), "--only", "clarinet=D3", "--only", "clarinet=G9"]), contextlib.redirect_stderr(err):
                self.assertEqual(B.main(), 0)
            card = json.loads((some / "clarinet" / "card.json").read_text())
            self.assertEqual(([p["note"] for p in card["pairs"]], card["withheld"]), (["D3"], ["E3"]))
            self.assertFalse((some / "clarinet" / "E3").exists(), "a withheld note ships no files")
            self.assertEqual(card["status"]["counts"]["notes_with_render"], 1, "the badge covers the notes shown")
            self.assertIn("clarinet=G9", err.getvalue(), "a note named in --only that does not exist is reported")

    def test_recording_only_notes_are_left_out_of_a_family_with_no_render_and_only_then(self):
        def build(root: Path, out: Path, *flag: str) -> dict:
            with mock.patch("sys.argv", ["build_instruments.py", "--analysis", str(root), "--out", str(out), *flag]):
                self.assertEqual(B.main(), 0)
            return json.loads((out / "clarinet" / "card.json").read_text())

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "analysis"
            make_pair(root, with_ours=False)  # a recording-only note, and nothing rendered in the family
            card = build(root, Path(tmp) / "all")
            self.assertEqual(([p["note"] for p in card["pairs"]], card["withheld"]), (["D3"], []))
            card = build(root, Path(tmp) / "lean", "--recordings-only-if-rendered")
            self.assertEqual(([p["note"] for p in card["pairs"]], card["withheld"]), ([], ["D3"]))
            self.assertFalse((Path(tmp) / "lean" / "clarinet" / "D3").exists(), "a withheld note ships no files")
            # Once the family has a render, its recording-only notes are shown again.
            e3 = root / "clarinet" / "E3"
            shutil.copytree(root / "clarinet" / "D3", e3)
            make_listening(root, "clarinet", "F3")
            e3_pair = json.loads((e3 / "pair.json").read_text())
            (e3 / "pair.json").write_text(json.dumps(dict(e3_pair, note="E3")))
            card = build(root, Path(tmp) / "rendered", "--recordings-only-if-rendered")
            self.assertEqual(([p["note"] for p in card["pairs"]], card["withheld"]), (["F3", "D3", "E3"], []))

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


def tally(t: tuple) -> dict:
    return {"inside": t[0], "outside": t[1], "no_ruler": t[2]}


def rec(slot: str, note: str, measures: int, loose: tuple, strict: tuple, out_loose: list, out_strict: list, retuned: bool = False) -> dict:
    """A row of P2's recording_list: one real recording graded like a render."""
    return {"slot": slot, "note": note, "recording": slot, "reference_file": "ref.aiff", "retuned_by_authors": retuned, "measures": measures,
            "loose": tally(loose), "strict": tally(strict), "outside_loose": out_loose, "outside_strict": out_strict}


def fraction(t: dict) -> dict:
    n = t["inside"] + t["outside"]
    return {**t, "inside_fraction": None if n == 0 else round(t["inside"] / n, 4)}


# Two graded real recordings of D3 (loose 5+6 inside, 3+2 outside; strict 4+6 inside, 3+1 outside), the way P2's file lists them.
ROWS_D3 = [rec("D3-s1", "D3", 10, (5, 3, 2), (4, 3, 3), ["both_outside", "loose_only", "no_claim"], ["both_outside", "loose_only", "both_inside"]),
           rec("D3-s2", "D3", 8, (6, 2, 0), (6, 1, 1), ["both_outside", "loose_only"], ["both_outside"], retuned=True)]


def p2_family(rows=None, refused=(), live=None, **over) -> dict:
    """A family entry of cross-library-baseline.json that agrees with itself, built the way P2's counting script builds it, so a test
    breaks exactly one thing at a time."""
    rows = copy.deepcopy(ROWS_D3 if rows is None else rows)
    live = CrossLibraryBaseline.LIVE if live is None else live
    add = lambda ts: {k: sum(t[k] for t in ts) for k in ("inside", "outside", "no_ruler")}  # noqa: E731
    nothing = lambda r, kind: r[kind]["outside"] == 0 and r[kind]["inside"] > 0  # noqa: E731  (a recording nothing of which was judged is not "nothing outside")
    real = {kind: add([r[kind] for r in rows]) for kind in ("loose", "strict")}
    have = [r for r in rows if r["note"] in live]
    ours = {kind: add([live[r["note"]][kind] for r in have]) for kind in ("loose", "strict")} if have else {}
    out = {kind: collections.Counter(m for r in rows for m in r[f"outside_{kind}"]) for kind in ("loose", "strict")}
    fam = {
        "source": {"name": "TEST LIBRARY, v1", "url": "https://example.org/record/0", "licence": "CC BY 4.0", "licence_url": "https://creativecommons.org/licenses/by/4.0/", "credit": "TEST CREDIT: A, B and C"},
        "independence": "other_library", "compared_with": "TEST reference (a different room, microphone and player)", "grader": ["decompose-1/instrgrade-2"],
        "recordings": len(rows), "refused": [dict(x) for x in refused], "retuned_by_authors": sum(r["retuned_by_authors"] for r in rows),
        "real_recordings": {kind: fraction(real[kind]) for kind in real},
        "recordings_with_nothing_outside": {"loose": sum(nothing(r, "loose") for r in rows), "strict": sum(nothing(r, "strict") for r in rows), "both": sum(nothing(r, "loose") and nothing(r, "strict") for r in rows)},
        "our_renders_same_notes": {**{kind: fraction(ours[kind]) for kind in ours}, "render_pairs": len(have)},
        "ruler_rows_per_measure": {"loose": {"median": 5, "min": 0, "max": 6}, "strict": {"median": 2, "min": 0, "max": 3}},
        "by_measure": {mid: {kind: {"outside": out[kind].get(mid, 0), "ruled": len(rows)} for kind in out} for mid in ("both_inside", "both_outside", "loose_only", "no_claim")},
        "recording_list": rows,
    }
    fam.update(over)
    return fam


class CrossLibraryBaseline(unittest.TestCase):
    """What a real recording from another library earns against our rulers, beside the render's own share (P2's file)."""
    LIVE = {"D3": {"measures": 5, "loose": {"inside": 3, "outside": 1, "no_ruler": 1}, "strict": {"inside": 1, "outside": 2, "no_ruler": 2}}}
    LABELS = {"both_outside": "Both outside", "loose_only": "Loose only"}

    def entry(self, family=None, live=None):
        return B.baseline_entry("clarinet", p2_family() if family is None else family, self.LIVE if live is None else live, self.LABELS)

    def refusal(self, family, live=None) -> str:
        with self.assertRaises(B.Refused) as cm:
            self.entry(family, live)
        return "; ".join(cm.exception.problems)

    def test_the_real_side_is_summed_from_the_recordings_and_the_render_side_comes_from_the_live_pairs(self):
        e, disagreements = self.entry(p2_family(refused=[{"slot": "D3-s3", "reason": "IRPDECOMP_ONSET_MISALIGNED: cut onset 17.5 ms"}]))
        self.assertEqual(disagreements, [])
        self.assertEqual(e["real"]["loose"], {"inside": 11, "outside": 5, "no_ruler": 2, "pct": 69})
        self.assertEqual(e["real"]["strict"], {"inside": 10, "outside": 4, "no_ruler": 4, "pct": 71})
        # one render pair per recording, repeats counted: the render's D3 verdicts twice
        self.assertEqual(e["render"]["loose"], {"inside": 6, "outside": 2, "no_ruler": 2, "pct": 75})
        self.assertEqual(e["render"]["strict"], {"inside": 2, "outside": 4, "no_ruler": 4, "pct": 33})
        self.assertEqual((e["render"]["recordings"], e["render"]["missing_notes"]), (2, []))
        self.assertEqual(e["recordings"], {"graded": 2, "refused": [{"slot": "D3-s3", "reason": "IRPDECOMP_ONSET_MISALIGNED: cut onset 17.5 ms"}], "retuned": 1, "notes": ["D3"], "nothing_outside": {"loose": 0, "strict": 0, "both": 0}})
        self.assertEqual((e["source"]["name"], e["source"]["licence"], e["independence"], e["grader"]), ("TEST LIBRARY, v1", "CC BY 4.0", "other_library", ["decompose-1/instrgrade-2"]))
        self.assertEqual(e["ruler_rows"]["loose"]["median"], 5)

    def test_a_recording_none_of_whose_measures_had_a_ruler_was_not_judged_and_is_not_nothing_outside(self):
        rows = ROWS_D3 + [rec("D3-s9", "D3", 6, (0, 0, 6), (0, 0, 6), [], [])]  # like the voice take that had no ruler left
        e, _ = self.entry(p2_family(rows))
        self.assertEqual(e["recordings"]["nothing_outside"], {"loose": 0, "strict": 0, "both": 0})
        self.assertEqual(e["real"]["loose"]["no_ruler"], 8)
        self.assertIn("recordings_with_nothing_outside says", self.refusal(p2_family(rows, recordings_with_nothing_outside={"loose": 1, "strict": 1, "both": 1})))

    def test_the_measures_are_listed_worst_first_with_the_pages_own_labels(self):
        e, _ = self.entry()
        self.assertEqual([(m["id"], m["label"], m["loose"]["outside"], m["loose"]["ruled"]) for m in e["by_measure"]],
                         [("both_outside", "Both outside", 2, 2), ("loose_only", "Loose only", 2, 2), ("no_claim", "no_claim", 1, 2), ("both_inside", "both_inside", 0, 2)])
        self.assertEqual(e["by_measure"][2]["strict"]["outside"], 0)

    def test_the_percent_is_of_the_measures_that_had_a_ruler_and_rounds_half_up(self):
        self.assertEqual(B.share({"inside": 1, "outside": 1, "no_ruler": 9})["pct"], 50)
        self.assertEqual(B.share({"inside": 1, "outside": 7, "no_ruler": 0})["pct"], 13)  # 12.5
        self.assertEqual(B.share({"inside": 125, "outside": 6, "no_ruler": 0})["pct"], 95)
        self.assertIsNone(B.share({"inside": 0, "outside": 0, "no_ruler": 3})["pct"])  # nothing had a ruler: no percent, never 0%
        self.assertIsNone(B.share(None))

    def test_a_recording_whose_note_has_no_graded_render_is_listed_and_not_counted_on_the_render_side(self):
        rows = ROWS_D3 + [rec("E3-s1", "E3", 4, (4, 0, 0), (4, 0, 0), [], [])]
        e, disagreements = self.entry(p2_family(rows))
        self.assertEqual(disagreements, [])
        self.assertEqual((e["render"]["recordings"], e["render"]["missing_notes"]), (2, ["E3"]))
        self.assertEqual(e["real"]["loose"]["inside"], 15)  # the real side still counts every graded recording
        self.assertEqual((e["recordings"]["notes"], e["recordings"]["nothing_outside"]), (["D3", "E3"], {"loose": 1, "strict": 1, "both": 1}))
        e, _ = self.entry(p2_family(rows, live={}), live={})
        self.assertEqual((e["render"]["loose"], e["render"]["strict"], e["render"]["missing_notes"]), (None, None, ["D3", "E3"]))

    def test_a_strict_share_is_not_computed_when_a_live_pair_lacks_the_strict_verdicts(self):
        old = dict(self.LIVE["D3"], strict=None)
        e, _ = self.entry(live={"D3": old})
        self.assertIsNone(e["render"]["strict"])
        self.assertEqual(e["render"]["loose"]["pct"], 75)

    def test_a_file_that_disagrees_with_itself_is_refused(self):
        def broken(mutate):
            fam = p2_family()
            mutate(fam)
            return self.refusal(fam)

        self.assertIn("real_recordings.loose says", broken(lambda f: f["real_recordings"]["loose"].update(inside=125, outside=6, no_ruler=4)))
        self.assertIn("real_recordings.strict says", broken(lambda f: f["real_recordings"]["strict"].update(outside=5)))
        self.assertIn("inside_fraction", broken(lambda f: f["real_recordings"]["loose"].update(inside_fraction=0.95)))
        self.assertIn("recordings says 9", broken(lambda f: f.update(recordings=9)))
        self.assertIn("retuned_by_authors says 5", broken(lambda f: f.update(retuned_by_authors=5)))
        self.assertIn("recordings_with_nothing_outside says", broken(lambda f: f["recordings_with_nothing_outside"].update(loose=2)))
        self.assertIn("by_measure.both_outside.loose.outside is 1 but the rows list it outside 2 times", broken(lambda f: f["by_measure"]["both_outside"]["loose"].update(outside=1)))
        self.assertIn("needs outside <= ruled", broken(lambda f: f["by_measure"]["both_outside"]["strict"].update(ruled=0)))
        self.assertIn("by_measure.loose_only.loose needs outside <= ruled", broken(lambda f: f["by_measure"].pop("loose_only")))  # a measure the rows list as missed, absent from the table

    def test_a_recording_row_whose_numbers_do_not_add_up_is_refused(self):
        def one(mutate):
            row = copy.deepcopy(ROWS_D3[0])
            mutate(row)
            return self.refusal(p2_family([row, ROWS_D3[1]]))

        self.assertIn("is not its 10 measures", one(lambda r: r["loose"].update(no_ruler=1)))
        self.assertIn("strict: inside 4 + outside 3 + no_ruler 4 is not its 10 measures", one(lambda r: r["strict"].update(no_ruler=4)))
        self.assertIn("outside_loose lists 3 measures but loose.outside is 2", one(lambda r: r["loose"].update(inside=6, outside=2)))
        self.assertIn("whole number", one(lambda r: r.update(measures=10.5)))
        self.assertIn("whole numbers", one(lambda r: r["loose"].update(inside=True)))
        self.assertIn("appears twice", self.refusal(p2_family([ROWS_D3[0], ROWS_D3[0]])))
        self.assertIn("needs a slot and a note", one(lambda r: r.update(slot="")))

    def test_the_source_the_independence_and_the_rulers_must_be_stated(self):
        self.assertIn("source needs a name", self.refusal(p2_family(source={"name": "x", "url": "https://x", "licence": "CC BY 4.0", "credit": " "})))
        self.assertIn("independence must be one of", self.refusal(p2_family(independence="independent")))
        self.assertIn("compared_with", self.refusal(p2_family(compared_with="")))
        self.assertIn("recording_list must be a non-empty list", self.refusal(p2_family([], recordings=0)))
        self.assertIn("ruler_rows_per_measure.strict", self.refusal(p2_family(ruler_rows_per_measure={"loose": {"median": 5, "min": 0, "max": 6}})))
        self.assertIn("refused must list", self.refusal(p2_family(refused=[{"slot": "x"}])))

    def test_a_render_figure_the_file_states_wrongly_is_reported_and_the_recomputed_one_is_shown(self):
        fam = p2_family()
        fam["our_renders_same_notes"]["loose"] = fraction(tally((60, 2, 2)))
        e, disagreements = self.entry(fam)
        self.assertEqual(e["render"]["loose"]["pct"], 75)  # the live pairs' own figure
        self.assertEqual(len(disagreements), 1)
        self.assertIn("clarinet loose: the file says the renders were (60, 2, 2)", disagreements[0])
        fam = p2_family()
        fam["our_renders_same_notes"]["render_pairs"] = 3
        self.assertIn("3 render pairs", self.entry(fam)[1][0])

    def test_the_render_side_is_recomputed_from_the_raw_pair_never_trusted(self):
        pj = {"ours": {"wav": "x"}, "real": {"reference_id": REF_ID}, "metrics": STRICT_CASES}
        self.assertEqual(B.live_verdicts(pj), self.LIVE["D3"])
        self.assertIsNone(B.live_verdicts(dict(pj, ours=None)))                     # a recording-only note
        self.assertIsNone(B.live_verdicts(dict(pj, refused={"reason": "x"})))       # a grader refusal
        self.assertIsNone(B.live_verdicts(dict(pj, metrics=[])))
        self.assertIsNone(B.live_verdicts(dict(pj, metrics=[{"id": "broken"}])))    # a pair this generator cannot read makes no claim

    def build(self, root: Path, out: Path, *flag: str, families=None, raw=None, **top):
        doc = {"schema": B.BASELINE_SCHEMA, "generated_utc": "2026-10-21T00:00:00Z", "method": ["Graded like a render."], "caveats": ["A small sample."], "families": families, **top}
        if families is not None:
            (root / B.BASELINE_FILE).write_text(json.dumps(doc))
        if raw is not None:
            (root / B.BASELINE_FILE).write_text(raw)
        err = io.StringIO()
        with mock.patch("sys.argv", ["build_instruments.py", "--analysis", str(root), "--out", str(out), *flag]), contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
            code = B.main()
        return code, err.getvalue(), lambda fam: json.loads((out / fam / "card.json").read_text())

    def test_the_card_carries_the_baseline_and_a_family_with_no_entry_carries_none(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "analysis"
            make_pair(root, metrics=STRICT_CASES)
            code, err, card = self.build(root, Path(tmp) / "out", families={"clarinet": p2_family()})
            self.assertEqual(code, 0, err)
            b = card("clarinet")["baseline"]
            self.assertTrue(b["shown"])
            self.assertEqual((b["real"]["loose"]["pct"], b["render"]["loose"]["pct"], b["render"]["strict"]["pct"]), (69, 75, 33))
            self.assertEqual((b["method"], b["caveats"], b["generated_utc"]), (["Graded like a render."], ["A small sample."], "2026-10-21T00:00:00Z"))
            self.assertEqual(b["by_measure"][0]["label"], "both_outside")  # the label of a measure the family's pairs know
            self.assertIsNone(card("flute")["baseline"])
            self.assertIsNone(card("oboe")["baseline"])

    def test_a_note_the_build_withholds_still_counts_on_the_render_side_of_the_baseline(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "analysis"
            make_pair(root, metrics=STRICT_CASES)
            code, err, card = self.build(root, Path(tmp) / "out", "--only", "clarinet=Z9", families={"clarinet": p2_family()})
            self.assertEqual(card("clarinet")["withheld"], ["D3"])
            self.assertEqual(card("clarinet")["baseline"]["render"]["loose"]["pct"], 75)

    def test_a_refused_family_baseline_says_why_and_the_exit_code_says_so(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "analysis"
            make_pair(root, metrics=STRICT_CASES)
            bad = p2_family(recordings=9)
            code, err, card = self.build(root, Path(tmp) / "out", families={"clarinet": bad, "ghost": p2_family()})
            self.assertEqual(code, 2)
            self.assertFalse(card("clarinet")["baseline"]["shown"])
            self.assertIn("recordings says 9", card("clarinet")["baseline"]["why"])
            self.assertIn("REFUSED baseline for clarinet", err)
            self.assertIn("'ghost', which is not on the page", err)

    def test_a_stale_render_figure_still_builds_with_the_live_figure_and_a_warning_and_exit_2(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "analysis"
            make_pair(root, metrics=STRICT_CASES)
            fam = p2_family()
            fam["our_renders_same_notes"]["strict"] = fraction(tally((9, 9, 0)))
            code, err, card = self.build(root, Path(tmp) / "out", families={"clarinet": fam})
            self.assertEqual(code, 2)
            self.assertIn("WARNING baseline out of date", err)
            self.assertEqual(card("clarinet")["baseline"]["render"]["strict"]["pct"], 33)

    def test_no_baseline_file_means_no_baseline_and_an_unreadable_one_stops_the_build(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "analysis"
            make_pair(root, metrics=STRICT_CASES)
            code, err, card = self.build(root, Path(tmp) / "out")
            self.assertEqual(code, 0)
            self.assertIsNone(card("clarinet")["baseline"])
            for n, raw in enumerate(('{"schema": "frankensim.instrument-cross-library-baseline.v1", "fam', json.dumps({"schema": "wrong", "families": {}, "method": ["x"], "caveats": ["y"]}),
                                     json.dumps({"schema": B.BASELINE_SCHEMA, "families": {}, "caveats": ["y"]}), json.dumps({"schema": B.BASELINE_SCHEMA, "families": {}, "method": ["x"], "caveats": []}))):
                with self.assertRaises(SystemExit):
                    self.build(root, Path(tmp) / f"out-{n}", raw=raw)


if __name__ == "__main__":
    unittest.main()
