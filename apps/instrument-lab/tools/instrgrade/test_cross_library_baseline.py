from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cross_library_baseline as cl  # noqa: E402

TINYSOL = cl.SOURCES["tinysol"]["prefix"] + " Strings/Violin/ordinario/"


def metric(mid, loose, strict, n_loose=3, n_strict=1):
    """One measure: inside_* True/False/None (None = no ruler), a gap unless the measure has none."""
    return {"id": mid, "gap": 1.0, "inside_loose": loose, "inside_strict": strict, "ruler": [0] * n_loose, "ruler_strict": [0] * n_strict}


def counts(metrics):
    out = {}
    for kind in ("loose", "strict"):
        v = [m["inside_" + kind] for m in metrics if m["gap"] is not None]
        out[kind] = {"inside": v.count(True), "outside": v.count(False), "no_ruler": v.count(None)}
    return {"measures": sum(m["gap"] is not None for m in metrics), **out}


def pair(family, note, metrics, wav="Vn-ord-A4-mf-2c-N.wav", observer=None, ref="Violin.arco.mf.A4.aif"):
    return {
        "family": family, "note": note, "metrics": metrics, "verdicts": counts(metrics),
        "real": {"file": ref}, "ours": {"wav": "/x/" + wav, "observer": observer or TINYSOL + wav},
        "decompose": {"version": "decompose-1", "grader": "instrgrade-2"},
    }


def write(root, family, slot, doc):
    d = Path(root) / family / slot
    d.mkdir(parents=True)
    (d / "pair.json").write_text(json.dumps(doc))


class BaselineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.certify, self.live = Path(self.tmp.name) / "certify", Path(self.tmp.name) / "live"
        a = [metric("a", True, True), metric("b", False, True), metric("c", True, None), metric("d", None, None)]
        b = [metric("a", True, True), metric("b", True, True), metric("c", True, True), metric("d", True, False)]
        write(self.certify, "violin", "A4-mf-s2", pair("violin", "A4-mf-s2", a))
        write(self.certify, "violin", "A4-mf-s3", pair("violin", "A4-mf-s3", b, wav="Vn-ord-A4-mf-3c-T16u.wav"))
        write(self.certify, "violin", "E5-mf-s1", {"refused": {"reason": "IRPDECOMP_ONSET_MISALIGNED: cut onset 17.5 ms"}})
        write(self.live, "violin", "A4", pair("violin", "A4", [metric("a", True, True), metric("b", False, False), metric("c", False, None), metric("d", None, None)]))
        write(self.live, "cello", "C2", pair("cello", "C2", [metric("a", True, True)]))
        self.doc = cl.build(self.certify, self.live)
        self.violin = self.doc["families"]["violin"]

    def tearDown(self):
        self.tmp.cleanup()

    def test_pooled_counts_leave_out_measures_without_a_ruler(self):
        real = self.violin["real_recordings"]
        # loose: first recording 2 in, 1 out, 1 no ruler; second 4 in. strict: first 2 in, 0 out, 2 no ruler; second 3 in, 1 out.
        self.assertEqual(real["loose"], {"inside": 6, "outside": 1, "no_ruler": 1, "inside_fraction": round(6 / 7, 4)})
        self.assertEqual(real["strict"], {"inside": 5, "outside": 1, "no_ruler": 2, "inside_fraction": round(5 / 6, 4)})

    def test_a_refusal_is_listed_and_not_counted(self):
        self.assertEqual(self.violin["recordings"], 2)
        self.assertEqual([r["slot"] for r in self.violin["refused"]], ["E5-mf-s1"])

    def test_our_renders_are_the_live_pair_of_the_same_note_once_per_recording(self):
        ours = self.violin["our_renders_same_notes"]
        self.assertEqual(ours["render_pairs"], 2)
        # the live A4 pair has loose 1 in / 2 out and strict 1 in / 1 out (one measure each without a verdict); counted for both recordings
        self.assertEqual((ours["loose"]["inside"], ours["loose"]["outside"]), (2, 4))
        self.assertEqual((ours["strict"]["inside"], ours["strict"]["outside"]), (2, 2))

    def test_page_rule_inside_needs_zero_outside_and_the_weaker_ruler_speaks(self):
        self.assertEqual(cl.badge(5, 0), "Inside the real spread")
        self.assertEqual(cl.badge(5, 1), "Partly inside the real spread")
        self.assertEqual(cl.badge(1, 1), "Partly inside the real spread")
        self.assertEqual(cl.badge(1, 2), "Outside the real spread")
        self.assertIsNone(cl.badge(0, 0))
        self.assertEqual(cl.weaker("Inside the real spread", "Outside the real spread"), "Outside the real spread")
        self.assertEqual(cl.weaker(None, "Partly inside the real spread"), "Partly inside the real spread")
        self.assertEqual(self.violin["real_recordings"]["badge_by_page_rule"], "Partly inside the real spread")

    def test_recordings_with_nothing_outside_are_counted_per_ruler(self):
        # first recording: loose has one outside; second: loose clean, strict has one outside
        self.assertEqual(self.violin["recordings_with_nothing_outside"], {"loose": 1, "strict": 1, "both": 0})

    def test_by_measure_counts_misses_and_the_recordings_that_have_a_ruler(self):
        bm = self.violin["by_measure"]
        self.assertEqual(bm["b"]["loose"], {"outside": 1, "ruled": 2})
        self.assertEqual(bm["d"]["loose"], {"outside": 0, "ruled": 1})  # recording one has no ruler there
        self.assertEqual(bm["d"]["strict"], {"outside": 1, "ruled": 1})
        self.assertEqual(self.doc["all_recordings_by_measure"]["b"]["loose"], {"outside": 1, "ruled": 2})

    def test_retuned_recordings_are_marked_from_the_file_name(self):
        flags = {r["slot"]: r["retuned_by_authors"] for r in self.violin["recording_list"]}
        self.assertEqual(flags, {"A4-mf-s2": False, "A4-mf-s3": True})
        self.assertEqual(self.violin["retuned_by_authors"], 1)

    def test_families_on_the_page_with_no_recording_are_listed_as_without_a_baseline(self):
        self.assertEqual(self.doc["families_without_baseline"], ["cello"])

    def test_an_observer_text_that_names_no_known_source_stops_the_export(self):
        write(self.certify, "viola", "A4-mf-s1", pair("viola", "A4-mf-s1", [metric("a", True, True)], observer="a render"))
        with self.assertRaises(SystemExit):
            cl.build(self.certify, self.live)

    def test_voice_is_labelled_as_another_singer_not_another_library(self):
        obs = cl.SOURCES["vocalset"]["prefix"] + " m6_long_straight_a.wav, held out of the ruler"
        write(self.certify, "voice", "C3-a", pair("voice", "C3-a", [metric("a", True, True)], wav="m6.wav", observer=obs))
        write(self.live, "voice", "C3-a", pair("voice", "C3-a", [metric("a", False, False)]))
        got = cl.build(self.certify, self.live)["families"]["voice"]
        self.assertEqual(got["independence"], "other_singer_same_library")
        self.assertFalse(got["recording_list"][0]["retuned_by_authors"])
        self.assertEqual(got["our_renders_same_notes"]["loose"]["outside"], 1)
        self.assertEqual(self.violin["independence"], "other_library")


if __name__ == "__main__":
    unittest.main()
