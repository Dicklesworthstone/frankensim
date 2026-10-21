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


def note_of(x, f0=F0, kind="sustained", max_cents=50.0):
    floor = dc.render_floor(dc.highpass(x, f0), f0)
    y, info = dc.align(x, f0, floor)
    return Note(y, f0, kind, False, "synthetic", info, max_cents=max_cents)


class AlignmentTests(unittest.TestCase):
    def test_delayed_copy_aligns_to_the_millisecond(self):
        a, b = note_of(tone()), note_of(tone(delay=0.137, noise_db=-80))
        self.assertLess(abs(a.info["cut_onset_s"] - b.info["cut_onset_s"]) * 1000, dc.ALIGN_TOL_MS)
        self.assertAlmostEqual(a.info["cut_onset_s"], dc.preroll_s(F0), delta=0.0015)
        self.assertGreater(abs(a.info["onset_s"] - b.info["onset_s"]), 0.1)  # the cut really moved them

    def test_legato_reattack_is_found_through_the_previous_notes_tail(self):
        n = int(2.0 * SR)
        t = np.arange(n) / SR
        first = np.where(t < 0.5, np.exp(-(t - 0.0) * 6), 0) * np.sin(2 * np.pi * 233 * t)  # decaying tail, still audible at 0.5 s
        env = np.clip((t - 0.5) / 0.03, 0, 1)
        second = 3.0 * env * np.sin(2 * np.pi * 247 * t)
        x = 0.1 * (first + second) / 3.0
        # the search starts 20 ms before the new note, inside the tail
        ta, tb = dc.onset_times(dc.highpass(x, 247.0), -100.0, anchor_s=0.52, f0_hz=247.0)
        self.assertAlmostEqual(max(ta, tb), 0.5, delta=0.01)

    def test_slow_ramp_from_rest_with_no_quiet_frames_aligns(self):
        # a bowed render: starts from silence at sample 0, ramps over 0.4 s, never has a quiet 20 ms frame afterwards
        n = note_of(tone(attack=0.4, lead=0.0, dur=3.0))
        self.assertAlmostEqual(n.info["cut_onset_s"], dc.preroll_s(F0), delta=0.0015)
        again = note_of(n.y)  # cutting at the onset and detecting again finds the same crossing
        self.assertAlmostEqual(again.info["cut_onset_s"], n.info["cut_onset_s"], delta=0.0015)

    def test_note_that_begins_on_a_loud_tail_is_refused_not_guessed(self):
        n = int(2.0 * SR)
        t = np.arange(n) / SR
        x = 0.1 * (np.sin(2 * np.pi * 233 * t) * np.where(t < 0.5, 1.0, 0) + 1.0 * np.clip((t - 0.5) / 0.03, 0, 1) * np.sin(2 * np.pi * 247 * t))
        with self.assertRaisesRegex(dc.OnsetMisaligned, "NO_QUIET_BEFORE_ONSET"):
            dc.onset_times(dc.highpass(x, 247.0), -100.0, anchor_s=0.52, f0_hz=247.0)

    def test_anchor_inside_a_slow_rise_finds_the_start_of_the_rise(self):
        x = tone(attack=0.2, lead=0.5)
        y = dc.highpass(x, F0)
        ta, tb = dc.onset_times(y, -100.0, anchor_s=0.9, f0_hz=F0)  # the anchor lags well inside the ramp
        first = dc.onset_times(y, -100.0, f0_hz=F0)
        self.assertAlmostEqual(max(ta, tb), max(first), delta=0.002)

    def test_a_recording_that_begins_at_its_own_onset_aligns_only_from_the_file_start(self):
        n = int(2.0 * SR)
        tt = np.arange(n) / SR
        x = 0.2 * np.exp(-tt * 2.0) * np.sin(2 * np.pi * 196 * tt)  # a plucked string: loudest at sample 0
        floor = -90.0
        cut, info = dc.align(x, 196.0, floor, anchor_s=0.02, at_start=True)
        self.assertAlmostEqual(info["cut_onset_s"], dc.preroll_s(196.0), delta=0.0015)
        self.assertLess(abs(info["onset_s"]), 0.01)  # in the time of the file: its first sample
        with self.assertRaisesRegex(dc.OnsetMisaligned, "NO_QUIET_BEFORE_ONSET"):
            dc.align(x, 196.0, floor, anchor_s=0.02, at_start=False)  # a window cut from the middle of a file: legato
        render, rinfo = dc.align(x, 196.0, floor)  # a render that starts at sample 0 is the same note
        self.assertAlmostEqual(rinfo["cut_onset_s"], dc.preroll_s(196.0), delta=0.0015)

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


class PitchTests(unittest.TestCase):
    def weak_fundamental(self, f0, cents=0.0):
        n = int(2.6 * SR)
        t = np.arange(n) / SR
        f = f0 * 2 ** (cents / 1200)
        amps = {1: 0.35, 2: 1.0, 3: 0.2, 4: 0.12, 5: 0.06}  # the second partial dominates, as in a bowed bass
        x = sum(a * np.sin(2 * np.pi * k * f * t + k) for k, a in amps.items())
        return 0.2 * x / np.abs(x).max() * (0.5 - 0.5 * np.cos(np.pi * np.clip((t - 0.3) / 0.05, 0, 1)))

    def test_comb_finds_the_fundamental_when_the_second_partial_dominates(self):
        for cents in (0.0, 4.0, -30.0, 120.0):
            x = self.weak_fundamental(55.0, cents)
            y = dc.highpass(x, 55.0)
            cut, _ = dc.align(y, 55.0, dc.render_floor(y, 55.0))
            f = dc.f0_comb(cut, 55.0)
            self.assertAlmostEqual(1200 * math.log2(f / (55.0 * 2 ** (cents / 1200))), 0.0, delta=1.0, msg=str(cents))

    def test_note_pitch_gap_is_the_true_detune_not_an_estimator_artifact(self):
        ref = note_of(self.weak_fundamental(55.0), 55.0)
        for cents in (4.0, 150.0):
            ours = note_of(self.weak_fundamental(55.0, cents), 55.0, max_cents=400.0)
            self.assertAlmostEqual(dc.compare(ours, ref)["f0_cents"][1], cents, delta=2.0)

    def test_a_tone_with_nothing_at_its_fundamental_is_refused(self):
        n = int(2.6 * SR)
        tt = np.arange(n) / SR
        env = 0.5 - 0.5 * np.cos(np.pi * np.clip((tt - 0.3) / 0.05, 0, 1))
        for f in (1558.2, 2047.3):  # the flute lock-on tones against a D4 / G4 name
            x = 0.2 * env * (np.sin(2 * np.pi * f * tt) + 0.2 * np.sin(2 * np.pi * 3 * f * tt))
            with self.assertRaisesRegex(dc.GraderRefusalError, "NOT_THE_NAMED_NOTE|PITCH_MISMATCH"):
                note_of(x, 293.66 if f < 2000 else 392.0, max_cents=400.0)

    def test_a_stretched_harmonic_note_passes_the_fundamental_check(self):
        n = int(2.6 * SR)
        tt = np.arange(n) / SR
        x = sum(np.sin(2 * np.pi * k * 110.0 * math.sqrt(1 + 6e-4 * k * k) * tt) / k for k in range(1, 16))
        x = 0.2 * x * (0.5 - 0.5 * np.cos(np.pi * np.clip((tt - 0.3) / 0.05, 0, 1))) / np.abs(x).max()
        note_of(x, 110.0)

    def test_a_reference_more_than_fifty_cents_off_its_name_is_refused(self):
        with self.assertRaisesRegex(dc.GraderRefusalError, "PITCH_MISMATCH"):
            note_of(self.weak_fundamental(55.0, 90.0), 55.0)


# the keys the page generator reads from a references.json entry: id, licence (a CREDITS key) and the url of the recording
REFERENCE_ENTRIES = {i: dict(id=i, licence="uiowa-mis", url=f"https://theremin.music.uiowa.edu/{i}.aiff") for i in ("syn-mf", "syn-ff")}


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
        problems = bi.validate_pair(pj, family, note, base, REFERENCE_ENTRIES) + bi.check_audio(base, pj)
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
        cy, ci = dc.align(raw, F0 * 2 ** (1 / 12), dc.render_floor(raw, F0 * 2 ** (1 / 12)))
        ci = dict(ci, cut_onset_s=ci["cut_onset_s"] + 0.01)  # the render's onset sits 10 ms late in its own cut
        meta = dict(wav="r.wav", git_sha="x", observer="o", sample_rate=48000, recipe="r.sh")
        with self.assertRaises(dc.OnsetMisaligned):
            dc.write_pair(self.tmp / "bad", "syn", "Bb3", ref_e, self.refs, y, info, cy, ci, meta)

    def test_reference_only_run_never_replaces_a_graded_pair(self):
        root = self.tmp / "keep"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="abc", observer="test", recipe="r.sh")
        dc.run_pair(self.refs, c, root)
        self.assertTrue(dc.has_render(root, "syn", "Bb3"))
        self.assertFalse(dc.has_render(root, "syn", "A3"))
        dc._REFS = self.refs
        res = dc._ref_job((dict(family="syn", art="", dyn="mf", notes=["Bb3"], review=False), root))
        self.assertIn("untouched", res["kept"])
        self.assertIsNotNone(__import__("json").loads((root / "syn" / "Bb3" / "pair.json").read_text())["ours"])

    def test_reference_only_run_leaves_a_refusal_record_of_a_render_alone(self):
        root = self.tmp / "keeprefusal"
        dc.refuse_pair(root, "syn", "Bb3", "IRPDECOMP_NOT_THE_NAMED_NOTE: test", render_sha="zz", render_wav="r.wav", dyn="mf")
        self.assertTrue(dc.has_render(root, "syn", "Bb3"))
        dc.refuse_pair(root, "syn", "A3", "stale pair superseded")  # no render_sha: it names no render
        self.assertFalse(dc.has_render(root, "syn", "A3"))

    def test_a_restarted_reference_run_skips_what_it_already_regenerated(self):
        import time
        root = self.tmp / "resume"
        dc._REFS = self.refs
        ref_e = self.refs.find("syn", 58, "mf")[0]
        y, info = self.refs.aligned(ref_e)
        dc.write_pair(root, "syn", "Bb3", ref_e, self.refs, y, info)
        job = dict(family="syn", art="", dyn="mf", notes=["Bb3"], review=False)
        self.assertEqual(dc._ref_job((dict(job, since=time.time() - 60), root))["kept"], "already regenerated since --since")
        res = dc._ref_job((dict(job, since=time.time() + 60), root))  # an older pair than --since is regenerated
        self.assertNotIn("kept", res)
        self.assertEqual(len(list((root / "syn" / "Bb3" / "superseded").iterdir())), 1)

    def test_reference_only_run_falls_back_to_the_next_verifying_note(self):
        root = self.tmp / "fallback"
        dc._REFS = self.refs
        res = dc._ref_job((dict(family="syn", art="", dyn="mf", notes=["C9", "Bb3"], review=False), root))
        self.assertTrue(res["ok"])
        self.assertEqual(res["note"], "Bb3")
        self.assertEqual([r["note"] for r in res["refused"]], ["C9"])

    def test_refusal_supersedes_by_moving_and_never_deletes(self):
        import json
        root = self.tmp / "refuse"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="abc", observer="test", recipe="r.sh")
        dc.run_pair(self.refs, c, root)
        before = {p.name for p in (root / "syn" / "Bb3").iterdir()}
        self.assertIn("ours.wav", before)
        dc.refuse_pair(root, "syn", "Bb3", "IRPDECOMP_PITCH_MISMATCH: test", "abc", self.render, "mf")
        d = root / "syn" / "Bb3"
        rec = json.loads((d / "pair.json").read_text())
        self.assertIsNone(rec["ours"])
        self.assertEqual(rec["refused"]["render_sha"], "abc")
        self.assertIn("PITCH_MISMATCH", rec["refused"]["reason"])
        stamps = list((d / "superseded").iterdir())
        self.assertEqual(len(stamps), 1)
        self.assertEqual({p.name for p in stamps[0].iterdir()}, before)  # everything moved, including the old pair.json
        self.assertEqual({p.name for p in d.iterdir()}, {"pair.json", "superseded"})
        dc.refuse_pair(root, "syn", "Bb3", "again")  # a second refusal moves only the first record, keeping the first stamp
        self.assertEqual(len(list((d / "superseded").iterdir())), 2)

    def test_a_regrade_moves_the_earlier_files_aside_instead_of_overwriting(self):
        root = self.tmp / "regrade"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="abc", observer="test", recipe="r.sh")
        dc.run_pair(self.refs, c, root)
        first = (root / "syn" / "Bb3" / "ours.wav").read_bytes()
        dc.run_pair(self.refs, dict(c, git_sha="def"), root)
        d = root / "syn" / "Bb3"
        stamps = list((d / "superseded").iterdir())
        self.assertEqual(len(stamps), 1)
        self.assertEqual((stamps[0] / "ours.wav").read_bytes(), first)
        self.assertTrue((stamps[0] / "pair.json").exists())
        self.assertEqual(__import__("json").loads((d / "pair.json").read_text())["ours"]["git_sha"], "def")

    def test_a_recording_outside_references_json_carries_its_own_source(self):
        import json
        refs = dc.Refs(self.tmp / "lib.json", self.tmp / "wav", self.tmp / "refs.json")
        refs.external = dict(name="Set", url="https://example.org/set", licence="CC BY 4.0", licence_url="https://example.org/by", credit="Someone")
        root = self.tmp / "external"
        ref_e = refs.find("syn", 58, "mf")[0]
        y, info = refs.aligned(ref_e)
        pj = dc.write_pair(root, "voice", "Bb3", ref_e, refs, y, info)
        self.assertEqual(pj["real"]["source"]["credit"], "Someone")
        self.assertEqual(pj["real"]["source"]["id"], "syn-mf")
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
        import build_instruments as bi
        base = root / "voice" / "Bb3"
        self.assertEqual(bi.validate_pair(json.loads((base / "pair.json").read_text()), "voice", "Bb3", base, {}), [])

    def test_a_render_that_locks_onto_another_tone_is_refused_for_that_not_for_its_onset(self):
        import json
        root = self.tmp / "lock"
        n = int(2.4 * SR)
        tt = np.arange(n) / SR
        x = 0.2 * np.sin(2 * np.pi * 1915.8 * tt) * np.clip(tt / 0.01, 0, 1)  # sounds at once, on the wrong tone
        from instrgrade import write_wav
        wav = self.tmp / "lock.wav"
        write_wav(wav, x)
        res = dc._batch_job((dict(family="syn", note="Bb3", dyn="mf", wav=str(wav), git_sha="q", observer="o", recipe="r.sh"), root))
        self.assertFalse(res["ok"])
        self.assertIn("NOT_THE_NAMED_NOTE", res["error"])
        self.assertIn("NOT_THE_NAMED_NOTE", json.loads((root / "syn" / "Bb3" / "pair.json").read_text())["refused"]["reason"])

    def test_diff_reports_what_moved_between_two_gradings_of_a_note(self):
        from instrgrade import write_wav
        root = self.tmp / "diffroot"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="old1", observer="o", recipe="r.sh")
        dc.run_pair(self.refs, c, root)
        far = self.tmp / "far.wav"
        write_wav(far, 0.3 * tone(f0=F0 * 2 ** (1 / 12), dur=2.4, lead=0.2, noise_db=-85, cents=60))
        dc.run_pair(self.refs, dict(c, wav=str(far), git_sha="new2"), root)
        d = root / "syn" / "Bb3"
        old = __import__("json").loads(next((d / "superseded").glob("*/pair.json")).read_text())
        new = __import__("json").loads((d / "pair.json").read_text())
        rows, tot = dc.diff_pairs(old, new)
        f0 = next(r for r in rows if r["id"] == "f0_cents")
        self.assertAlmostEqual(f0["gap_old"], 8.0, delta=1.5)
        self.assertAlmostEqual(f0["gap_new"], 60.0, delta=2.0)
        self.assertGreater(f0["gap_new"], f0["gap_old"])
        self.assertNotEqual(old["ours"]["git_sha"], new["ours"]["git_sha"])
        self.assertEqual(sum(tot["old"]), sum(tot["new"]))  # the same measures are judged both times

    def test_a_band_limited_render_is_graded_only_inside_its_band(self):
        import json
        from instrgrade import write_wav
        root = self.tmp / "band"
        f0 = F0 * 2 ** (1 / 12)
        full = 0.3 * tone(f0=f0, dur=2.4, lead=0.2, noise_db=-85)
        limited = self.tmp / "limited.wav"
        write_wav(limited, dc.lowpass(full, 3500.0))  # a render that models nothing above 3.5 kHz
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(limited), git_sha="b", observer="o", recipe="r.sh")
        wide = dc.run_pair(self.refs, c, root / "wide")
        narrow = dc.run_pair(self.refs, dict(c, band_hz=3900.0), root / "narrow")
        gap = lambda pj, m: next(x["gap"] for x in pj["metrics"] if x["id"] == m)  # noqa: E731
        self.assertLess(gap(narrow, "stft_lsd_db"), 0.6 * gap(wide, "stft_lsd_db"))
        self.assertEqual(narrow["figures"]["spectrogram"]["f_hi_hz"], 3900.0)
        self.assertEqual(narrow["figures"]["scalogram"]["f_hi_hz"], 3900.0)
        self.assertEqual(narrow["alignment"]["band_hz"], 3900.0)
        self.assertIn("low-passed at 3900 Hz", narrow["audio"]["level_match"])
        self.assertIsNone(wide["alignment"]["band_hz"])

    def test_a_failed_regrade_leaves_the_notes_earlier_pair_in_place(self):
        import json
        from unittest import mock
        root = self.tmp / "atomic"
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.render), git_sha="first", observer="o", recipe="r.sh")
        dc.run_pair(self.refs, c, root)
        d = root / "syn" / "Bb3"
        before = (d / "pair.json").read_bytes()
        with mock.patch.object(self.refs, "ruler", side_effect=RuntimeError("boom")):
            with self.assertRaises(RuntimeError):
                dc.run_pair(self.refs, dict(c, git_sha="second"), root)
        self.assertEqual((d / "pair.json").read_bytes(), before)  # the old grading is still there, untouched, mid-failure
        aborted = list((d / "superseded").glob("aborted-*"))  # the failed attempt's partial files are kept aside, not deleted
        self.assertEqual(len(aborted), 1)
        self.assertTrue(any(aborted[0].iterdir()))
        dc.run_pair(self.refs, dict(c, git_sha="third"), root)
        self.assertEqual(json.loads((d / "pair.json").read_text())["ours"]["git_sha"], "third")
        self.assertFalse([p for p in d.iterdir() if p.name.startswith(".staging")])  # nothing left staged

    def test_a_refused_render_batch_job_writes_the_record(self):
        import json
        root = self.tmp / "refbatch"
        dc._REFS = self.refs
        c = dict(family="syn", note="Bb3", dyn="mf", wav=str(self.tmp / "nonexistent.wav"), git_sha="zz", observer="o", recipe="r.sh")
        res = dc._batch_job((c, root))
        self.assertFalse(res["ok"])
        self.assertEqual(json.loads((root / "syn" / "Bb3" / "pair.json").read_text())["refused"]["render_sha"], "zz")

    def test_json_has_no_nan(self):
        import json
        s = json.dumps(dc.clean({"a": float("nan"), "b": np.float32(1.5), "c": [np.inf]}), allow_nan=False)
        self.assertEqual(json.loads(s), {"a": None, "b": 1.5, "c": [None]})


if __name__ == "__main__":
    unittest.main()
