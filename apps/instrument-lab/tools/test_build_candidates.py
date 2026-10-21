"""What the reference-candidates data is allowed to claim: a record that fails a check is refused, a licence is never guessed, nothing is
said to have been listened to, and every clip plays at the same loudness.

Run with: uv run --with numpy python -m unittest discover -s apps/instrument-lab/tools -p 'test_build_candidates.py' -v
The clips below are throwaway sine tones written to a temp directory; none of it reaches the page.
"""
import contextlib
import copy
import hashlib
import io
import json
import math
import tempfile
import unittest
import wave
from pathlib import Path
from unittest import mock

import numpy as np

import build_candidates as C

STANDARD = "YouTube standard licence (all rights reserved)"
CC = "Creative Commons Attribution license (reuse allowed)"


def write_wav(path: Path, peak: float = 0.3, seconds: float = 1.0, sr: int = 44100) -> None:
    t = np.arange(int(sr * seconds)) / sr
    x = (peak * np.sin(2 * math.pi * 220 * t) * 32767).astype("<i2")
    with wave.open(str(path), "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(sr)
        w.writeframes(x.tobytes())


def record(root: Path, family: str = "violin", cid: str = "violin-yt-aaa", peak: float = 0.3, seconds: float = 1.0, licence: str = STANDARD, source_type: str = "youtube",
           usable=("G3",), facts=None, **over) -> dict:
    """A candidate as P2's README describes it, with its clip and facts.json written under <root>/<family>/<id>/."""
    d = root / family / cid
    d.mkdir(parents=True, exist_ok=True)
    write_wav(d / "clip.wav", peak, seconds)
    (d / "facts.json").write_text(json.dumps(facts if facts is not None else {"peak_dbfs": -10.5, "rms_dbfs": -23.0, "bandwidth_99pct_hz": 3000, "pitch": {"median_hz": 196.0, "nearest_note": "G3", "cents_from_note": 1.0, "fraction_within_30_cents_of_median": 0.9}}))
    rec = {
        "id": cid, "family": family, "instrument": "TEST violin open strings", "notes": ["G3", "D4"], "dynamic": "not stated", "articulation": "arco, open strings",
        "source_type": source_type, "url": "https://www.youtube.com/watch?v=aaa", "title": "TEST Violin Open Strings", "channel_or_author": "A Channel",
        "licence": licence, "licence_evidence": "yt-dlp license field of the watch page", "start_s": 12.0, "end_s": 12.0 + seconds,
        "conditions": "STATED BY THE SOURCE: a description. MEASURED by script on the clip: peak -10.5 dBFS. NOT LISTENED TO.",
        "sample_rate_hz": 44100, "reason": "two steady notes found by a pitch tracker",
        "files": {"clip": f"{cid}/clip.wav", "facts": f"{cid}/facts.json"},
        "sha256_clip": hashlib.sha256((d / "clip.wav").read_bytes()).hexdigest(), "usable_for_render_notes": list(usable),
    }
    rec.update(over)
    return rec


def write_doc(root: Path, family: str, cands: list, name: str = "candidates-external-test.json", looked=None, none_found=None, **over) -> None:
    doc = {"schema": C.SCHEMA, "family": family, "made_by": "test", "made_utc": "2026-10-21T00:00:00Z",
           "looked": looked if looked is not None else [{"source": "YouTube", "query_or_url": "violin open strings", "result": "8 results, 2 usable"}],
           "none_found": none_found, "candidates": cands, **over}
    (root / family).mkdir(parents=True, exist_ok=True)
    (root / family / name).write_text(json.dumps(doc))


def build(root: Path, out: Path, *flags: str):
    err = io.StringIO()
    with mock.patch("sys.argv", ["build_candidates.py", "--candidates", str(root), "--out", str(out), *flags]), contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
        code = C.main()
    return code, err.getvalue(), json.loads((out / "index.json").read_text())


def family(index: dict, fid: str) -> dict:
    return next(f for f in index["families"] if f["id"] == fid)


class Loudness(unittest.TestCase):
    def test_a_quiet_clip_is_brought_to_the_target_rms(self):
        x = 0.01 * np.sin(2 * math.pi * 220 * np.arange(44100) / 44100)
        rms = 20 * math.log10(math.sqrt(float(np.mean(x ** 2))))
        self.assertAlmostEqual(rms + C.playback_gain_db(x), C.TARGET_RMS_DBFS, places=1)

    def test_a_clip_that_would_peak_too_high_is_held_back_by_its_peak(self):
        x = np.zeros(44100)
        x[100] = 0.9  # a click: tiny RMS, big peak
        self.assertAlmostEqual(20 * math.log10(0.9) + C.playback_gain_db(x), C.PEAK_LIMIT_DBFS, places=1)

    def test_a_silent_clip_is_refused(self):
        with self.assertRaises(C.Refused):
            C.playback_gain_db(np.zeros(100))


class Coverage(unittest.TestCase):
    """usable_for_render_notes is optional and the hunters wrote it four ways: all four are read, none is refused for its shape."""

    def test_a_list_of_notes_is_the_graded_notes(self):
        self.assertEqual(C.coverage(["B3", " D4 ", ""])["graded"], ["B3", "D4"])
        self.assertEqual(C.coverage([])["graded"], [])
        self.assertEqual(C.coverage(None), {"graded": [], "reference_only": [], "other": [], "text": None})

    def test_a_map_from_note_to_category_is_split_by_category(self):
        c = C.coverage({"A3": "rendered_and_graded", "B5": "reference_only", "X": "listening_only_or_refused", "Y": "something else"})
        self.assertEqual((c["graded"], c["reference_only"], c["other"]), (["A3"], ["B5"], ["X", "Y"]))

    def test_a_sentence_is_kept_and_its_graded_notes_are_read_from_it(self):
        c = C.coverage("graded notes covered: G2, D3, A3; reference-only notes cover B5")
        self.assertEqual((c["graded"], c["text"]), (["G2", "D3", "A3"], "graded notes covered: G2, D3, A3; reference-only notes cover B5"))
        self.assertEqual(C.coverage("graded notes covered: C2, G2 (partly)")["graded"], ["C2", "G2"])
        self.assertEqual(C.coverage("graded notes covered: Bb4 and Ab3")["graded"], ["Bb4", "Ab3"])
        none = C.coverage("none of the notes in families.json for this family (extra pitches)")
        self.assertEqual((none["graded"], none["text"].startswith("none of the notes")), ([], True))
        self.assertEqual(C.coverage("reference-only notes covered: B5")["graded"], [], "a sentence that does not say graded is not read as graded")

    def test_any_other_shape_is_refused(self):
        for bad in (7, 1.5, True, ["G3", 4], {"G3": 1}, {1: "reference_only"}):
            self.assertIsNone(C.coverage(bad), bad)


class Licence(unittest.TestCase):
    def test_creative_commons_is_read_from_the_start_of_the_licence_never_searched_for(self):
        self.assertTrue(C.is_cc(CC))
        self.assertTrue(C.is_cc("CC BY-SA 4.0"))
        self.assertTrue(C.is_cc("CC0"))
        self.assertFalse(C.is_cc(STANDARD))
        self.assertFalse(C.is_cc("Standard YouTube licence (not Creative Commons): for the owner's listening review only"), "the old wording says 'not Creative Commons'")
        self.assertFalse(C.is_cc("Public domain"))


class Records(unittest.TestCase):
    def check(self, mutate=None, **kw):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rec = record(root, **kw)
            if mutate:
                mutate(rec, root / "violin" / "violin-yt-aaa")
            return C.validate(rec, "violin", root / "violin")

    def problems(self, mutate=None, **kw) -> str:
        return "; ".join(self.check(mutate, **kw)[0])

    def test_a_good_record_passes_and_its_facts_come_with_it(self):
        problems, facts = self.check()
        self.assertEqual(problems, [])
        self.assertEqual(facts["pitch"]["nearest_note"], "G3")

    def test_the_folder_is_the_family_and_the_id_is_filesystem_safe(self):
        self.assertIn("family 'cello'", self.problems(lambda r, d: r.update(family="cello")))
        for bad in ("a/b", "has space", "", "../up"):
            self.assertIn("id must be a filesystem-safe", self.problems(lambda r, d, bad=bad: r.update(id=bad)), bad)

    def test_a_licence_and_where_it_was_read_must_be_stated_never_guessed(self):
        for key in ("licence", "licence_evidence"):
            for bad in ("", "   ", None):
                self.assertIn(f"{key} must be stated", self.problems(lambda r, d, key=key, bad=bad: r.update({key: bad})), (key, bad))

    def test_nobody_may_claim_to_have_listened(self):
        self.assertIn("NOT LISTENED TO", self.problems(lambda r, d: r.update(conditions="Studio, close mic, sounds warm.")))
        self.assertIn("NOT LISTENED TO", self.problems(lambda r, d: r.update(conditions="")))
        self.assertEqual(self.problems(lambda r, d: r.update(conditions="Measured: peak -3 dBFS. Not listened to.")), "")

    def test_the_clip_must_be_the_file_the_record_describes(self):
        self.assertIn("sha256_clip", self.problems(lambda r, d: r.update(sha256_clip="0" * 64)))
        self.assertIn("lasts", self.problems(lambda r, d: r.update(end_s=r["start_s"] + 5.0)))
        self.assertIn("sample_rate_hz says 48000", self.problems(lambda r, d: r.update(sample_rate_hz=48000)))
        self.assertIn("is missing", self.problems(lambda r, d: r["files"].update(clip="nope/clip.wav")))
        self.assertIn("files.clip", self.problems(lambda r, d: r.update(files={})))
        self.assertIn("files.facts", self.problems(lambda r, d: r["files"].update(facts="nope/facts.json")))

    def test_the_other_fields_are_checked_too(self):
        self.assertIn("source_type must be one of", self.problems(lambda r, d: r.update(source_type="rumour")))
        self.assertIn("url must be an http(s) address", self.problems(lambda r, d: r.update(url="not a url")))
        self.assertIn("start_s < end_s", self.problems(lambda r, d: r.update(start_s=5, end_s=5)))
        self.assertIn("notes must be", self.problems(lambda r, d: r.update(notes=[])))
        self.assertIn("notes must be", self.problems(lambda r, d: r.update(notes="G3")))
        self.assertIn("sample_rate_hz must be a whole number", self.problems(lambda r, d: r.update(sample_rate_hz=44100.5)))
        self.assertIn("usable_for_render_notes", self.problems(lambda r, d: r.update(usable_for_render_notes=7)))
        self.assertIn("usable_for_render_notes", self.problems(lambda r, d: r.update(usable_for_render_notes=["G3", 4])))
        self.assertIn("usable_for_render_notes", self.problems(lambda r, d: r.update(usable_for_render_notes={"G3": 1})))
        self.assertIn("title must be a non-empty string", self.problems(lambda r, d: r.update(title=" ")))
        self.assertEqual(self.problems(lambda r, d: r.pop("usable_for_render_notes")), "")  # optional


class Build(unittest.TestCase):
    def test_a_candidate_is_built_with_its_clip_gain_licence_class_and_facts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "c", Path(tmp) / "out"
            write_doc(root, "violin", [record(root, cid="violin-yt-aaa"), record(root, cid="violin-yt-bbb", licence=CC)])
            code, err, index = build(root, out)
            self.assertEqual(code, 0, err)
            by_id = {c["id"]: c for c in family(index, "violin")["candidates"]}
            self.assertFalse(by_id["violin-yt-aaa"]["creative_commons"])
            self.assertTrue(by_id["violin-yt-bbb"]["creative_commons"])
            c = by_id["violin-yt-aaa"]
            self.assertTrue((out / c["audio"]).is_file())
            self.assertTrue(c["audio"].endswith(".mp4"))
            self.assertEqual(c["playback_gain_db"], C.playback_gain_db(C.decode(root / "violin" / "violin-yt-aaa" / "clip.wav")[0]))
            self.assertAlmostEqual(c["clip_duration_s"], 1.0, places=2)
            self.assertEqual(c["facts"]["pitch"]["nearest_note"], "G3")
            self.assertEqual((c["render_coverage"]["graded"], c["source_type"], c["licence_evidence"]), (["G3"], "youtube", "yt-dlp license field of the watch page"))
            self.assertEqual(family(index, "violin")["looked"][0]["query_or_url"], "violin open strings")
            self.assertEqual(family(index, "violin")["made_by"], ["test"])

    def test_until_a_family_is_merged_the_library_and_external_files_are_read_together_and_a_merged_file_wins(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "c", Path(tmp) / "out"
            write_doc(root, "violin", [record(root, cid="violin-lib-1", source_type="open-library")], name="candidates-library.json",
                      looked=[{"source": "IRCAM TinySOL", "query_or_url": "Zenodo 3685367", "result": "violin ordinario mf"}])
            write_doc(root, "violin", [record(root, cid="violin-yt-1")], name="candidates-external-p8.json")
            write_doc(root, "violin", [record(root, cid="violin-yt-2")], name="candidates-external-other.json")
            code, err, index = build(root, out)
            f = family(index, "violin")
            self.assertEqual(sorted(c["id"] for c in f["candidates"]), ["violin-lib-1", "violin-yt-1", "violin-yt-2"])
            self.assertEqual(sorted(x["source"] for x in f["looked"]), ["IRCAM TinySOL", "YouTube", "YouTube"])
            write_doc(root, "violin", [record(root, cid="violin-yt-1")], name="candidates.json")  # P2's merge
            code, err, index = build(root, Path(tmp) / "out2")
            self.assertEqual([c["id"] for c in family(index, "violin")["candidates"]], ["violin-yt-1"], "candidates.json is the family's file once it exists")

    def test_a_record_that_fails_a_check_is_listed_and_the_exit_code_says_so(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "c", Path(tmp) / "out"
            bad = record(root, cid="violin-yt-bad")
            bad["sha256_clip"] = "0" * 64
            write_doc(root, "violin", [record(root, cid="violin-yt-aaa"), bad, record(root, cid="violin-yt-aaa")])
            code, err, index = build(root, out)
            self.assertEqual(code, 2)
            f = family(index, "violin")
            self.assertEqual([c["id"] for c in f["candidates"]], ["violin-yt-aaa"])
            self.assertEqual(sorted((r["id"], "sha256_clip" in r["problems"][0] or "twice" in r["problems"][0]) for r in f["refused"]), [("violin-yt-aaa", True), ("violin-yt-bad", True)])
            self.assertIn("REFUSED violin/violin-yt-bad", err)
            self.assertFalse((out / "violin" / "violin-yt-bad.mp4").exists())

    def test_a_file_that_does_not_say_where_it_looked_is_refused_whole(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "c", Path(tmp) / "out"
            write_doc(root, "violin", [record(root)], looked=[])
            code, err, index = build(root, out)
            f = family(index, "violin")
            self.assertEqual((code, f["candidates"]), (2, []))
            self.assertIn("looked must list every source and query", f["refused"][0]["problems"][0])
            (root / "cello").mkdir()
            (root / "cello" / "candidates.json").write_text("{not json")
            code, err, index = build(root, Path(tmp) / "out2")
            self.assertIn("is not JSON", family(index, "cello")["refused"][0]["problems"][0])
            write_doc(root, "viola", [record(root, family="viola", cid="viola-yt-x")], name="candidates.json")
            doc = json.loads((root / "viola" / "candidates.json").read_text())
            doc["family"] = "cello"
            (root / "viola" / "candidates.json").write_text(json.dumps(doc))
            code, err, index = build(root, Path(tmp) / "out3")
            self.assertIn("family 'cello' is not 'viola'", family(index, "viola")["refused"][0]["problems"][0])

    def test_none_found_says_why_and_a_family_with_nothing_recorded_says_so(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "c", Path(tmp) / "out"
            write_doc(root, "duduk", [], none_found={"reason": "only melodies", "looked_summary": "youtube: 8 results, none a single note"},
                      looked=[{"source": "YouTube", "query_or_url": "duduk long note", "result": "8 results, melodies"}])
            code, err, index = build(root, out)
            self.assertEqual(code, 0, err)
            d = family(index, "duduk")
            self.assertEqual((d["candidates"], d["none_found"]["reason"], d["none_found"]["looked_summary"]), ([], "only melodies", "youtube: 8 results, none a single note"))
            m = family(index, "marimba")
            self.assertEqual((m["candidates"], m["none_found"], m["looked"]), ([], None, []))
            write_doc(root, "tuba", [], none_found=None)
            code, err, index = build(root, Path(tmp) / "out2")
            self.assertEqual(code, 2)
            self.assertIn("no candidates and no none_found.reason", family(index, "tuba")["refused"][0]["problems"][0])

    def test_only_the_best_ranked_few_of_a_family_are_built_and_the_rest_are_listed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "c", Path(tmp) / "out"
            lib = record(root, cid="violin-lib-g3", source_type="open-library", licence="CC BY 4.0", usable=("G3", "D4"))
            mid = record(root, cid="violin-yt-mid", usable=("G3",))
            low = record(root, cid="violin-yt-low", usable=(), facts={"pitch": {"fraction_within_30_cents_of_median": 0.1}})
            write_doc(root, "violin", [low, mid, lib])
            code, err, index = build(root, out, "--max-per-family", "2")
            f = family(index, "violin")
            self.assertEqual(([c["id"] for c in f["candidates"]], f["omitted"]), (["violin-lib-g3", "violin-yt-mid"], ["violin-yt-low"]))
            self.assertFalse((out / "violin" / "violin-yt-low.mp4").exists(), "an omitted candidate gets no file")

    def test_the_rank_counts_coverage_of_the_graded_notes_first(self):
        base = {"source_type": "youtube", "licence": STANDARD, "usable_for_render_notes": []}
        covers = {**base, "usable_for_render_notes": ["G3", "D4"]}
        library = {**base, "source_type": "open-library"}
        cc = {**base, "licence": CC}
        self.assertGreater(C.rank(covers, None), C.rank(library, None))
        self.assertGreater(C.rank(library, None), C.rank(base, None))
        self.assertGreater(C.rank(cc, None), C.rank(base, None))
        steady = {"pitch": {"fraction_within_30_cents_of_median": 1.0}}
        self.assertGreater(C.rank(base, steady), C.rank(base, {"pitch": {"fraction_within_30_cents_of_median": 0.2}}))
        self.assertEqual(C.rank(base, {"steady_notes_by_pitch_tracker": [{"clip_t0": 0, "clip_t1": 40}]}), C.rank(base, {"steady_notes_by_pitch_tracker": [{"clip_t0": 0, "clip_t1": 12}]}), "capped")

    def test_the_notes_our_renders_are_graded_on_are_passed_to_the_family(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "c", Path(tmp) / "out"
            write_doc(root, "violin", [record(root)])
            (root / "families.json").write_text(json.dumps({"families": {"violin": {"rendered_and_graded": ["A4", "G3"], "reference_only": ["B3"], "listening_only_or_refused": []}}}))
            code, err, index = build(root, out)
            self.assertEqual(family(index, "violin")["render_notes"]["rendered_and_graded"], ["A4", "G3"])
            self.assertIsNone(family(index, "cello")["render_notes"])

    def test_an_unknown_family_folder_still_gets_an_entry_and_scratch_folders_do_not(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "c", Path(tmp) / "out"
            write_doc(root, "theremin", [record(root, family="theremin", cid="theremin-yt-aaa")])
            (root / "_scratch").mkdir()
            code, err, index = build(root, out)
            self.assertEqual(len(family(index, "theremin")["candidates"]), 1)
            self.assertNotIn("_scratch", [f["id"] for f in index["families"]])

    def test_the_families_are_the_owners_list_in_his_order(self):
        self.assertEqual([f[0] for f in C.FAMILIES][:5], ["violin", "viola", "cello", "double-bass", "guitar"])
        self.assertEqual(len(C.FAMILIES), 22)


if __name__ == "__main__":
    unittest.main()
