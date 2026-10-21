"""What the #instruments cards are allowed to claim: refusals, the within/outside verdict and the badge.

Run with: uv run --with numpy --with pillow python -m unittest discover -s apps/instrument-lab/tools -p 'test_build_instruments.py' -v
The pairs below are throwaway inputs to the generator, written to a temp directory; none of it reaches the page.
"""
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
