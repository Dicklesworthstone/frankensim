"""What the blind A/B data is allowed to do: a pair that could give its answer away, or that has nothing to tell apart, is refused; the clips carry
opaque names; and the two sounds and the reference play at one loudness.

Run with: uv run --with numpy python -m unittest discover -s apps/instrument-lab/tools -p 'test_build_blind_ab.py' -v
The clips below are throwaway sine tones written to a temp directory; none of it reaches the page.
"""
import contextlib
import copy
import io
import json
import math
import tempfile
import unittest
import wave
from pathlib import Path
from unittest import mock

import numpy as np

import build_blind_ab as B
import build_candidates as C


def write_wav(path: Path, peak: float = 0.3, freq: float = 220.0, seconds: float = 1.0, sr: int = 44100) -> None:
    t = np.arange(int(sr * seconds)) / sr
    x = (peak * np.sin(2 * math.pi * freq * t) * 32767).astype("<i2")
    with wave.open(str(path), "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(sr)
        w.writeframes(x.tobytes())


def rec(pid: str = "p9-bars-1", **over) -> dict:
    r = {"schema": B.SCHEMA, "id": pid, "set": "P9 bars", "title": "Xylophone C5, pair 1", "question": "Which is closer to the real xylophone C5?",
         "change": "silicone vs boxwood mallet, everything else identical", "made_by": "P9",
         "a": {"label": "silicone mallet", "file": "a.wav", "commit": "0123456789ab", "recipe": "cargo test bar_a"},
         "b": {"label": "boxwood mallet", "file": "b.wav", "commit": "fedcba987654", "recipe": "cargo test bar_b"},
         "reference": {"label": "real xylophone C5", "file": "real.wav", "source": "UIowa MIS", "licence": "UIowa MIS terms"}}
    r.update(over)
    return r


def pair_dir(root: Path, pid: str = "p9-bars-1", peaks=(0.3, 0.1, 0.5), freqs=(220.0, 233.0, 220.0), seconds=(1.0, 1.0, 1.0), **over) -> Path:
    d = root / pid
    d.mkdir(parents=True)
    for name, peak, freq, sec in zip(("a.wav", "b.wav", "real.wav"), peaks, freqs, seconds):
        write_wav(d / name, peak=peak, freq=freq, seconds=sec)
    (d / "pair.json").write_text(json.dumps(rec(pid, **over)))
    return d


def build(pairs: Path | None, out: Path, *flags: str):
    err = io.StringIO()
    argv = ["build_blind_ab.py", "--out", str(out), *flags] + (["--pairs", str(pairs)] if pairs else [])
    with mock.patch("sys.argv", argv), contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
        code = B.main()
    return code, err.getvalue(), json.loads((out / "index.json").read_text())


def level(path: Path) -> float:
    x, _ = C.decode(path)
    return B.rms_db(x)


class Validate(unittest.TestCase):
    def problems(self, mutate=None, **over) -> str:
        with tempfile.TemporaryDirectory() as tmp:
            d = pair_dir(Path(tmp), **over)
            r = json.loads((d / "pair.json").read_text())
            if mutate:
                mutate(r, d)
            return "; ".join(B.validate(r, d))

    def test_a_good_pair_passes(self):
        self.assertEqual(self.problems(), "")

    def test_the_title_and_question_are_shown_before_the_answer_so_they_may_not_name_either_render(self):
        self.assertIn("contains the label 'silicone mallet'", self.problems(title="Xylophone, silicone mallet"))
        self.assertIn("contains the label 'boxwood mallet'", self.problems(question="Is the BOXWOOD mallet closer?"))
        self.assertIn("contains the commit of a", self.problems(title="Pair 0123456789ab"))

    def test_the_two_renders_must_be_tellable_apart_afterwards(self):
        self.assertIn("same", self.problems(lambda r, d: r["b"].update(label="Silicone Mallet")))

    def test_fields_are_checked(self):
        self.assertIn("schema must be", self.problems(schema="v0"))
        self.assertIn("not the folder's name", self.problems(id="other"))
        self.assertIn("id must be lower case", self.problems(lambda r, d: r.update(id="Bad Id")))
        for k in ("set", "title", "question", "change", "made_by"):
            self.assertIn(f"{k} must be a non-empty string", self.problems(lambda r, d, k=k: r.update({k: " "})), k)
        self.assertIn("a.commit must be a git sha", self.problems(lambda r, d: r["a"].update(commit="main")))
        self.assertIn("b.recipe must be", self.problems(lambda r, d: r["b"].update(recipe="")))
        self.assertIn("a.file must be a relative path", self.problems(lambda r, d: r["a"].update(file="../a.wav")))
        self.assertIn("a.file nope.wav is missing", self.problems(lambda r, d: r["a"].update(file="nope.wav")))
        self.assertIn("reference.source must be", self.problems(lambda r, d: r["reference"].update(source="")))
        self.assertIn("reference.file must be a relative path", self.problems(lambda r, d: r["reference"].update(file="/etc/passwd")))
        self.assertIn("reference must be an object", self.problems(lambda r, d: r.update(reference=None)))
        self.assertIn("level_match must be", self.problems(level_match="loud"))
        self.assertIn("placeholder must be true or false", self.problems(placeholder="yes"))


class Build(unittest.TestCase):
    def test_the_clips_get_opaque_names_and_the_data_names_which_is_which(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "pairs", Path(tmp) / "out"
            root.mkdir()
            pair_dir(root)
            code, err, index = build(root, out)
            self.assertEqual((code, index["refused"]), (0, []), err)
            p = index["pairs"][0]
            for slot in ("a", "b", "reference"):
                self.assertRegex(p[slot]["audio"], r"^p9-bars-1/[0-9a-f]{10}\.mp4$")
                self.assertTrue((out / p[slot]["audio"]).is_file())
            self.assertEqual((p["a"]["label"], p["b"]["commit"], p["reference"]["label"]), ("silicone mallet", "fedcba987654", "real xylophone C5"))
            self.assertEqual(sorted(f.name for f in (out / "p9-bars-1").iterdir() if f.is_file() and f.name.startswith(".")), [], "no temporary file is left")
            self.assertEqual(len({p[s]["audio"] for s in ("a", "b", "reference")}), 3)

    def test_the_two_sounds_and_the_reference_play_at_one_loudness(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "pairs", Path(tmp) / "out"
            root.mkdir()
            pair_dir(root, peaks=(0.3, 0.05, 0.6))
            _, _, index = build(root, out)
            p = index["pairs"][0]
            self.assertEqual(p["target_rms_dbfs"], C.TARGET_RMS_DBFS)
            levels = [level(out / p[s]["audio"]) for s in ("a", "b", "reference")]
            for v in levels:
                self.assertAlmostEqual(v, C.TARGET_RMS_DBFS, delta=0.6)
            self.assertGreater(p["b"]["gain_db"], p["a"]["gain_db"], "the quiet clip is brought up more")

    def test_a_peaky_clip_lowers_the_common_target_instead_of_clipping_or_unmatching_the_others(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "pairs", Path(tmp) / "out"
            root.mkdir()
            d = pair_dir(root)
            x = np.zeros(44100)
            x[::4410] = 0.9  # a click train: tiny RMS next to its peak, so reaching -26 dBFS RMS would put the peak far above -1 dBFS
            with wave.open(str(d / "a.wav"), "wb") as w:
                w.setnchannels(1); w.setsampwidth(2); w.setframerate(44100)
                w.writeframes((x * 32767).astype("<i2").tobytes())
            _, _, index = build(root, out)
            p = index["pairs"][0]
            self.assertLess(p["target_rms_dbfs"], C.TARGET_RMS_DBFS)
            self.assertLessEqual(p["a"]["peak_dbfs"] + p["a"]["gain_db"], C.PEAK_LIMIT_DBFS + 0.01)
            for s in ("a", "b", "reference"):
                self.assertAlmostEqual(p[s]["rms_dbfs"] + p[s]["gain_db"], p["target_rms_dbfs"], places=1)

    def test_level_match_none_plays_the_clips_as_rendered(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "pairs", Path(tmp) / "out"
            root.mkdir()
            pair_dir(root, level_match="none")
            _, _, index = build(root, out)
            p = index["pairs"][0]
            self.assertEqual((p["level_match"], p["target_rms_dbfs"], p["a"]["gain_db"], p["b"]["gain_db"], p["reference"]["gain_db"]), ("none", None, 0.0, 0.0, 0.0))

    def test_a_long_clip_plays_its_first_part_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "pairs", Path(tmp) / "out"
            root.mkdir()
            pair_dir(root, seconds=(40.0, 1.0, 1.0))
            _, _, index = build(root, out)
            p = index["pairs"][0]
            self.assertEqual((p["a"]["playback_s"], p["a"]["source_duration_s"]), (C.MAX_PLAY_S, 40.0))
            self.assertEqual(p["b"]["playback_s"], 1.0)

    def test_pairs_that_cannot_be_judged_are_refused_listed_and_the_rest_still_build(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "pairs", Path(tmp) / "out"
            root.mkdir()
            pair_dir(root, "ok-pair")
            same = pair_dir(root, "same-bytes")
            (same / "b.wav").write_bytes((same / "a.wav").read_bytes())
            pair_dir(root, "leaks", title="Xylophone with the silicone mallet")
            short = pair_dir(root, "too-short")
            write_wav(short / "b.wav", seconds=0.05)
            silent = pair_dir(root, "silent")
            write_wav(silent / "a.wav", peak=0.0)
            (root / "broken").mkdir()
            (root / "broken" / "pair.json").write_text("{not json")
            code, err, index = build(root, out)
            self.assertEqual(code, 2)
            self.assertEqual([p["id"] for p in index["pairs"]], ["ok-pair"])
            why = {r["id"]: "; ".join(r["problems"]) for r in index["refused"]}
            self.assertEqual(sorted(why), ["broken", "leaks", "same-bytes", "silent", "too-short"])
            self.assertIn("same bytes", why["same-bytes"])
            self.assertIn("contains the label", why["leaks"])
            self.assertIn("under 0.2 s", why["too-short"])
            self.assertIn("silent", why["silent"])
            self.assertIn("cannot be read", why["broken"])
            self.assertFalse((out / "leaks").exists() and any((out / "leaks").iterdir()), "a refused pair gets no files")

    def test_ids_are_unique_and_underscore_folders_are_scratch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "pairs", Path(tmp) / "out"
            root.mkdir()
            pair_dir(root, "one")
            pair_dir(root, "_scratch")
            _, _, index = build(root, out)
            self.assertEqual([p["id"] for p in index["pairs"]], ["one"])

    def test_pairs_sort_real_before_placeholder_then_by_set_and_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "pairs", Path(tmp) / "out"
            root.mkdir()
            pair_dir(root, "z-last", set="P9 bars")
            pair_dir(root, "a-first", set="P9 bars")
            pair_dir(root, "m-flute", set="P10 flute")
            pair_dir(root, "ph", placeholder=True)
            _, _, index = build(root, out)
            self.assertEqual([p["id"] for p in index["pairs"]], ["m-flute", "a-first", "z-last", "ph"])


class Placeholder(unittest.TestCase):
    def analysis(self, root: Path) -> Path:
        for note, freq in (("A4", 440.0), ("Bb4", 466.16)):
            d = root / "clarinet" / note
            d.mkdir(parents=True)
            write_wav(d / "ours.wav", peak=0.2, freq=freq)
            write_wav(d / "real.wav", peak=0.4, freq=440.0)
            (d / "pair.json").write_text(json.dumps({"ours": {"git_sha": "596fe876f", "recipe": "render.sh", "observer": "baffled_outlet 0 0 0.5"}}))
        return root

    def test_a_placeholder_pair_is_made_of_two_existing_renders_and_says_so(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out"
            code, _, index = build(None, out, "--placeholder-from", str(self.analysis(Path(tmp) / "analysis")))
            self.assertEqual((code, index["refused"]), (0, []))
            p = index["pairs"][0]
            self.assertTrue(p["placeholder"])
            self.assertEqual((p["id"], p["a"]["commit"], p["a"]["label"], p["b"]["label"]), (B.PLACEHOLDER_ID, "596fe876f", "our clarinet A4 render", "our clarinet Bb4 render"))
            self.assertIn("NOT A REAL TEST", p["change"])
            self.assertEqual(p["set"], "Page trial")
            for s in ("a", "b", "reference"):
                self.assertTrue((out / p[s]["audio"]).is_file())

    def test_no_clarinet_renders_means_no_placeholder_not_a_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out"
            (Path(tmp) / "empty").mkdir()
            code, err, index = build(None, out, "--placeholder-from", str(Path(tmp) / "empty"))
            self.assertEqual((code, index["pairs"]), (0, []))
            self.assertIn("no placeholder", err)

    def test_placeholder_and_real_pairs_build_together(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = Path(tmp) / "pairs", Path(tmp) / "out"
            root.mkdir()
            pair_dir(root, "real-one")
            _, _, index = build(root, out, "--placeholder-from", str(self.analysis(Path(tmp) / "analysis")))
            self.assertEqual([(p["id"], p["placeholder"]) for p in index["pairs"]], [("real-one", False), (B.PLACEHOLDER_ID, True)])


if __name__ == "__main__":
    unittest.main()
