"""refscreen.py on synthetic recordings whose onset, pitch, noise and decay are known by construction.

uv run --with numpy python3 -I test_refscreen.py -v      (needs ffmpeg and ffprobe on PATH)
Each source is written as a 16-bit WAV at 44.1 kHz, the way a library ships it, and goes through the same ffmpeg decode as a
real source. Each test prints what it built and what the screen measured, so a failure log says why.
"""
from __future__ import annotations

import json
import struct
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import refscreen as rs  # noqa: E402

SR_SRC = 44100


def write_pcm16(path, x, sr=SR_SRC):
    """[n] or [n, channels] floats in [-1, 1] as a 16-bit PCM WAV."""
    x = np.atleast_2d(np.asarray(x, np.float64).T).T
    q = np.clip(np.round(x * 32767), -32768, 32767).astype("<i2")
    data = q.tobytes()
    ch = x.shape[1]
    Path(path).write_bytes(b"RIFF" + struct.pack("<I", 36 + len(data)) + b"WAVEfmt " +
                           struct.pack("<IHHIIHH", 16, 1, ch, sr, sr * 2 * ch, 2 * ch, 16) + b"data" + struct.pack("<I", len(data)) + data)


def note(f0=440.0, lead=0.5, dur=5.0, decay_db_s=25.0, amp=0.3, noise_db=-75.0, seed=3, kmax=12, attack_s=0.004):
    """lead seconds of noise, then a struck harmonic note (1/k partials, the k-th decaying k^0.5 times faster) in that noise."""
    rng = np.random.default_rng(seed)
    n, n0 = int((lead + dur) * SR_SRC), int(lead * SR_SRC)
    t = np.arange(n - n0) / SR_SRC
    y = np.zeros(n)
    for k in range(1, kmax + 1):
        if k * f0 > 0.45 * SR_SRC:
            break
        y[n0:] += np.sin(2 * np.pi * k * f0 * t + 0.3 * k) / k * 10 ** (-decay_db_s * k ** 0.5 * t / 20)
    y[n0:] *= np.clip(t / attack_s, 0, 1)
    y *= amp / np.max(np.abs(y))
    return y + 10 ** (noise_db / 20) * rng.standard_normal(n)


class Screen(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def run_cut(self, x, name="A4", kind="decaying", **kw):
        src, out = self.d / "src.wav", self.d / "clip.wav"
        write_pcm16(src, x)
        row = rs.cut(src, out, name, kind, **kw)
        print(json.dumps({k: row.get(k) for k in ("note", "status", "failed", "trim_start_s", "trim_end_s", "snr_db", "cents_comb",
                                                   "cents_yin", "reattack_db", "tail_db", "clipped_runs")}))
        return row, src, out

    def test_clean_struck_note_is_valid_and_cut_around_its_onset(self):
        row, src, out = self.run_cut(np.stack([note(), 0.8 * note()], axis=1))
        self.assertEqual(row["status"], "VALID", row["failed"])
        # the clean note's loudest 20 ms over the -75 dB white noise it was built in, read from the clip's 60 ms lead-in and a 20 ms
        # frame grid of its own: within 1.5 dB
        clean = note(noise_db=-400.0)[int(0.5 * SR_SRC):]
        hop = int(0.02 * SR_SRC)
        true_snr = 10 * np.log10(np.max(np.mean(clean[: len(clean) // hop * hop].reshape(-1, hop) ** 2, axis=1))) + 75.0
        print(f"true peak-to-noise {true_snr:.2f} dB, screened {row['snr_db']:.2f} dB")
        self.assertLess(abs(row["snr_db"] - true_snr), 1.5)
        self.assertEqual(row["noise_method"], "lead-in")
        self.assertAlmostEqual(row["trim_start_s"], 0.5 - rs.LEAD_S, delta=0.004)
        self.assertGreater(row["trim_end_s"], 2.0)  # the note decays for seconds before it reaches its floor
        self.assertEqual(row["channels"], 2)
        self.assertEqual(row["sample_rate_hz"], SR_SRC)
        self.assertEqual(row["clip_sha256"], rs.sha256(out))
        x, sr = rs.decode(src)
        y, _ = rs.decode(out)
        i0 = int(round(row["trim_start_s"] * sr))
        np.testing.assert_array_equal(y, x[i0:i0 + len(y)])  # the clip IS the source's samples, unchanged

    def test_wrong_pitch_fails(self):
        row, *_ = self.run_cut(note(f0=440.0 * 2 ** (2 / 12)))  # a B4 named A4
        self.assertEqual(row["status"], "INVALID")
        self.assertEqual([f.split(":")[0] for f in row["failed"]], ["pitch"])

    def test_octave_below_is_caught_by_yin_not_by_the_harmonic_sum(self):
        row, *_ = self.run_cut(note(f0=220.0, kmax=24))  # an A3 named A4: every partial of A4 is a partial of A3
        self.assertLessEqual(abs(row["cents_comb"]), rs.PITCH_TOL_CENTS)
        self.assertAlmostEqual(row["cents_yin"], -1200.0, delta=20.0)
        self.assertEqual([f.split(":")[0] for f in row["failed"]], ["pitch"])

    def held(self, under_db):
        """A held A4 (2.5 s at a constant level, 50 ms release) between 0.3 s of white noise either side, the noise `under_db`
        under the note's loudest 20 ms (measured on the clean note): (signal, the noise's level under that peak as built). The
        note fills most of the recording, as a single-note reference does, so the grader's steady level (the median of frames
        within 40 dB of the loudest) is the note's; with more noise than note the median falls in the noise and even the grader's
        onset refuses the recording."""
        n0, n1, n = int(0.3 * SR_SRC), int(2.8 * SR_SRC), int(3.1 * SR_SRC)
        t = np.arange(n) / SR_SRC
        env = np.clip((t - 0.3) / 0.03, 0, 1) * np.where(t < 2.8, 1.0, np.exp(-(t - 2.8) / 0.05))
        clean = 0.2 * env * sum(np.sin(2 * np.pi * k * 440.0 * t + 0.3 * k) / k for k in range(1, 9))
        hop = int(0.02 * SR_SRC)
        peak = 10 * np.log10(np.max(np.mean(clean[n0:n1][: (n1 - n0) // hop * hop].reshape(-1, hop) ** 2, axis=1)))
        return clean + 10 ** ((peak - under_db) / 20) * np.random.default_rng(13).standard_normal(n), under_db

    def test_noise_inside_the_criterion_b_margin_is_not_measurable(self):
        # 12 dB under the note: inside STEADY_DROP_DB (20 dB), so no lead-in frame lies 20 dB under the note's steady level, and
        # the stationary tail lies within TAIL_UNDER_PEAK_DB (20 dB) of the peak: there is nothing quiet to read the noise in
        x, under = self.held(12.0)
        self.assertLess(under, min(rs.dc.STEADY_DROP_DB, rs.TAIL_UNDER_PEAK_DB))
        row, *_ = self.run_cut(x, kind="sustained")
        self.assertEqual(row["noise_method"], "none")
        self.assertIsNone(row["snr_db"])
        self.assertIn("snr", [f.split(":")[0] for f in row["failed"]])

    def test_noise_outside_the_margin_is_measured_and_fails_snr(self):
        # 30 dB under the note: outside the 20 dB margin, so the lead-in measures it, within 1.5 dB, and 30 dB is under SNR_MIN_DB
        x, under = self.held(30.0)
        self.assertGreater(under, rs.dc.STEADY_DROP_DB)
        self.assertLess(under, rs.SNR_MIN_DB)
        row, *_ = self.run_cut(x, kind="sustained")
        print(f"noise built {under:.2f} dB under the peak, screened snr {row['snr_db']}")
        self.assertEqual(row["noise_method"], "lead-in")
        self.assertLess(abs(row["snr_db"] - under), 1.5)
        self.assertEqual([f.split(":")[0] for f in row["failed"]], ["snr"])

    def test_second_strike_fails_single(self):
        x = note()
        i = int((0.5 + 0.8) * SR_SRC)
        x[i:] += note(lead=0.0, dur=len(x[i:]) / SR_SRC, seed=4)[: len(x) - i]
        # a strike 3.4 s in, past the grader's 3 s window, on a note that decays slowly enough to be still sounding then
        j = int((0.5 + 3.4) * SR_SRC)
        y = note(decay_db_s=6.0, dur=8.0)
        y[j:] += note(lead=0.0, dur=len(y[j:]) / SR_SRC + 0.01, seed=6, decay_db_s=6.0)[: len(y) - j]
        late, *_ = self.run_cut(y)
        self.assertIn("single", [f.split(":")[0] for f in late["failed"]])
        row, *_ = self.run_cut(x)
        self.assertIn("single", [f.split(":")[0] for f in row["failed"]])

    def test_beating_unison_strings_are_one_note(self):
        # two strings 1.5 Hz apart: the sum's envelope drops into a notch every 0.67 s and climbs back, 30 dB deep and more
        n0, n = int(0.4 * SR_SRC), int(6.4 * SR_SRC)
        t = np.arange(n - n0) / SR_SRC
        y = np.zeros(n)
        for k in range(1, 9):
            for f in (830.6, 832.1):
                y[n0:] += np.sin(2 * np.pi * k * f * t + 0.4 * k) / k * 10 ** (-8.0 * k ** 0.5 * t / 20)
        y[n0:] *= np.clip(t / 0.003, 0, 1)
        y = 0.3 * y / np.max(np.abs(y)) + 1e-4 * np.random.default_rng(11).standard_normal(n)
        row, *_ = self.run_cut(y, name="Ab5")
        self.assertNotIn("single", [f.split(":")[0] for f in row["failed"]])
        self.assertLess(row["reattack_db"], rs.REATTACK_DB)

    def test_clipped_recording_fails_clip(self):
        row, *_ = self.run_cut(np.clip(note(amp=1.6), -1.0, 1.0))
        self.assertIn("clip", [f.split(":")[0] for f in row["failed"]])
        self.assertGreater(row["clipped_runs"], 0)

    def test_note_cut_off_before_its_decay_fails_decay(self):
        row, *_ = self.run_cut(note(decay_db_s=3.0), until=1.5)
        self.assertEqual([f.split(":")[0] for f in row["failed"]], ["decay"])

    def test_unpitched_strike_needs_no_pitch(self):
        rng = np.random.default_rng(9)
        n0, n = int(0.4 * SR_SRC), int(2.4 * SR_SRC)
        t = np.arange(n - n0) / SR_SRC
        burst = 0.4 * rng.standard_normal(n - n0) * np.exp(-t / 0.12)
        # band-limited to 18 kHz as a real source is: a white burst carries full energy at 22.05 kHz, and the grader's FFT
        # resample to 48 kHz (instrgrade.load, refscreen.to_sr) spreads that ahead of the onset (it read 50 dB for 71)
        B = np.fft.rfft(burst)
        B[np.fft.rfftfreq(len(burst), 1 / SR_SRC) > 18000.0] = 0.0
        burst = np.fft.irfft(B, len(burst))
        x = 1e-4 * rng.standard_normal(n)
        x[n0:] += burst
        row, *_ = self.run_cut(x, name=None)
        self.assertEqual(row["status"], "VALID", row["failed"])
        self.assertNotIn("cents_comb", row)
        hop = int(0.02 * SR_SRC)
        true_snr = 10 * np.log10(np.max(np.mean(burst[: len(burst) // hop * hop].reshape(-1, hop) ** 2, axis=1))) + 80.0
        print(f"true peak-to-noise {true_snr:.2f} dB, screened {row['snr_db']:.2f} dB")
        self.assertLess(abs(row["snr_db"] - true_snr), 1.5)

    def test_hint_cuts_the_second_of_two_separated_notes(self):
        a = note(lead=0.3, dur=2.0, decay_db_s=40.0)
        b = note(f0=440.0 * 2 ** (3 / 12), lead=0.7, dur=3.0, seed=5)
        row, *_ = self.run_cut(np.concatenate([a, b]), name="C5", hint=2.3 + 0.7 + 0.3)
        self.assertEqual(row["status"], "VALID", row["failed"])
        self.assertAlmostEqual(row["trim_start_s"], 2.3 + 0.7 - rs.LEAD_S, delta=0.004)

    def test_sustained_note_with_a_release(self):
        n0, n1, n = int(0.4 * SR_SRC), int(2.4 * SR_SRC), int(3.4 * SR_SRC)
        t = np.arange(n) / SR_SRC
        env = np.clip((t - 0.4) / 0.06, 0, 1) * np.where(t < 2.4, 1.0, np.exp(-(t - 2.4) / 0.05))
        vib = 2 * np.pi * np.cumsum(330.0 * 2 ** (12 * np.sin(2 * np.pi * 5.5 * t) / 1200)) / SR_SRC
        x = 0.2 * env * sum(np.sin(k * vib) / k for k in range(1, 10)) + 3e-5 * np.random.default_rng(2).standard_normal(n)
        row, *_ = self.run_cut(x, name="E4", kind="sustained")
        self.assertEqual(row["status"], "VALID", row["failed"])
        self.assertGreater(row["trim_end_s"], 2.4)
        self.assertLess(row["trim_end_s"], 3.0)
        del n0, n1

    def test_library_sample_cut_at_its_attack_that_fades_out_has_no_measurable_noise(self):
        n = int(3.0 * SR_SRC)
        t = np.arange(n) / SR_SRC
        x = 0.2 * np.clip(t / 0.08, 0, 1) * sum(np.sin(2 * np.pi * k * 523.25 * t) / k for k in range(1, 8))
        x += 3e-3 * np.random.default_rng(8).standard_normal(n)  # breath and room at -50 dB: under the note, never alone
        x *= np.clip((3.0 - t) / 0.3, 0, 1)  # the library's fade to digital silence
        x[-int(0.05 * SR_SRC):] = 0.0
        row, *_ = self.run_cut(x, name="C5", kind="sustained")
        self.assertEqual(row["noise_method"], "none")
        self.assertIsNone(row["snr_db"])
        self.assertEqual([f.split(":")[0] for f in row["failed"]], ["snr"])

    def test_library_sample_cut_at_its_attack_that_decays_into_noise_reads_the_tail(self):
        row, *_ = self.run_cut(note(lead=0.0, dur=4.0, decay_db_s=30.0, noise_db=-80.0))
        self.assertEqual(row["noise_method"], "tail")
        self.assertEqual(row["status"], "VALID", row["failed"])
        self.assertGreater(row["snr_db"], 50.0)

    def test_loop_sample_cut_mid_sustain_is_not_its_own_noise(self):
        n = int(2.0 * SR_SRC)
        t = np.arange(n) / SR_SRC
        x = 0.4 * np.clip(t / 0.05, 0, 1) * sum(np.sin(2 * np.pi * k * 698.46 * t) / k for k in range(1, 8))
        x += 1e-3 * np.random.default_rng(12).standard_normal(n)  # breath under the note; the file stops while the note sounds
        row, *_ = self.run_cut(x, name="F5", kind="sustained")
        self.assertEqual(row["noise_method"], "none")
        self.assertEqual(sorted(f.split(":")[0] for f in row["failed"]), ["decay", "snr"])

    def test_verify_recomputes_each_clips_hash_and_verdict(self):
        row, _, out = self.run_cut(note())
        (self.d / "piano").mkdir()
        clip = self.d / "piano" / "piano-test-A4.wav"
        out.rename(clip)
        entry = dict(id="piano-test-A4", family="piano", note="A4", kind="decaying", clip_sha256=row["clip_sha256"], status=row["status"])
        # an unpitched clip carries no note at all, as the manifest writes it
        rng = np.random.default_rng(9)
        hit = 1e-4 * rng.standard_normal(int(2.0 * SR_SRC))
        hit[int(0.4 * SR_SRC):] += 0.4 * rng.standard_normal(int(1.6 * SR_SRC)) * np.exp(-np.arange(int(1.6 * SR_SRC)) / SR_SRC / 0.05)
        write_pcm16(self.d / "hit.wav", hit)
        (self.d / "kit").mkdir()
        drum = rs.cut(self.d / "hit.wav", self.d / "kit" / "kit-test-hit.wav", None, "decaying")
        unpitched = dict(id="kit-test-hit", family="kit", piece="snare", kind="decaying", clip_sha256=drum["clip_sha256"], status=drum["status"])
        man = self.d / "manifest.json"
        man.write_text(json.dumps(dict(files=[dict(id="src", members=[dict(clips=[entry])]), dict(id="hit", clips=[unpitched])])))
        args = rs.argparse.Namespace(manifest=str(man), root=str(self.d))
        self.assertEqual(rs.cmd_verify(args), 0)
        man.write_text(json.dumps(dict(files=[dict(id="src", clips=[dict(entry, status="INVALID")])])))
        self.assertEqual(rs.cmd_verify(args), 1)  # a verdict the clip does not earn
        b = bytearray(clip.read_bytes())
        b[-8:] = bytes(8)
        clip.write_bytes(bytes(b))
        man.write_text(json.dumps(dict(files=[dict(id="src", clips=[entry])])))
        self.assertEqual(rs.cmd_verify(args), 1)  # a clip that is not the one pinned


class Coverage(unittest.TestCase):
    def test_counts_valid_pitches_per_dynamic(self):
        rows = [dict(family="recorder", note=n, dynamic=d, status=s, failed=[] if s == "VALID" else ["pitch: x"])
                for n, d, s in [("C5", "mf", "VALID"), ("D5", "mf", "VALID"), ("E5", "mf", "VALID"), ("F5", "pp", "VALID"),
                                ("G5", "mf", "INVALID")]]
        rows.append(dict(family="kit", note=None, piece="snare", dynamic="v3", status="VALID", failed=[]))
        cov = rs.coverage(dict(rows=rows))
        print(json.dumps(cov))
        self.assertEqual(cov["recorder"]["by_dynamic"], {"mf": ["C5", "D5", "E5"], "pp": ["F5"]})
        self.assertTrue(cov["recorder"]["done_when"])
        self.assertEqual(cov["recorder"]["failed_checks"], {"pitch": 1})
        self.assertEqual(cov["kit"]["by_dynamic"], {"v3": ["snare"]})
        self.assertFalse(cov["kit"]["done_when"])


if __name__ == "__main__":
    unittest.main()
