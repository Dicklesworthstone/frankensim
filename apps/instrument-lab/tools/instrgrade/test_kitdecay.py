"""kitdecay.py on synthetic signals whose band decays are known by construction.

uv run --with numpy python3 -I test_kitdecay.py -v
"""
from __future__ import annotations

import math
import struct
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kitdecay as kd  # noqa: E402

SR = 48000


def tone(fc, t60, seconds=3.0, start=0.1):
    t = np.arange(int(seconds * SR)) / SR
    return np.where(t >= start, np.sin(2 * math.pi * fc * (t - start)) * 10 ** (-3 * (t - start) / t60), 0.0)


def floors():
    impulse = np.zeros(3 * SR)
    impulse[SR // 2] = 1.0
    return [r["decay_ms"] for r in kd.decays(impulse, SR)]


class KnownDecays(unittest.TestCase):
    def test_an_exponential_tone_reads_a_third_of_its_t60_in_its_own_band(self):
        f = floors()
        for k, fc in enumerate(kd.OCT):
            for t60 in (0.15, 0.5, 2.0):
                row = kd.decays(tone(fc, t60), SR, f)[k]
                want = 1000.0 * t60 / 3.0
                print(f"{fc} Hz, T60 {t60} s: want {want:.0f} ms, got {row['mark']}{row['decay_ms']}")
                # one 10 ms frame, plus the band's own impulse resolution
                self.assertEqual(row["mark"], "")
                self.assertLessEqual(abs(row["decay_ms"] - want), 10.0 + f[k])

    def test_the_impulse_resolution_floor_is_at_most_two_frames(self):
        self.assertTrue(all(v <= 20.0 for v in floors()), floors())

    def test_a_band_that_never_falls_twenty_db_is_a_lower_bound(self):
        row = kd.decays(tone(500, 30.0, seconds=1.0), SR)[3]
        self.assertEqual(row["mark"], ">")
        # measured to the end of the analysed window, EDGE_S before the clip's end
        end_ms = 1000.0 * (1.0 - kd.EDGE_S)
        self.assertAlmostEqual(row["decay_ms"], end_ms - 1000.0 * row["peak_s"], delta=10.0)

    def test_the_analysis_stops_before_the_band_filters_response_to_the_cut(self):
        # The first run's tom-12 still rang at its last sample, and the band filter spread that
        # step into its final frames. A tom-12-like steady spectrum stopped dead shows the
        # same: a strong 125 Hz tone, every other band 25 dB down. The uncut 63 Hz band rises
        # at the very end; inside the analysed window no band does.
        t = np.arange(2 * SR) / SR
        x = np.sin(2 * math.pi * 125 * t) + 10 ** (-25 / 20) * sum(
            np.sin(2 * math.pi * fc * t) for fc in kd.OCT if fc != 125)

        def rise(e):
            return float(e[-3:].max() - np.median(e))

        self.assertGreater(rise(kd.band_energy_db(x, SR)[0]), 4.0)
        for fc, e in zip(kd.OCT, kd.analysed_band_energy_db(x, SR)):
            self.assertLess(rise(e), 2.0, fc)

    def test_a_recording_band_within_26_db_of_its_floor_is_not_measurable(self):
        noisy = tone(1000, 2.0, start=0.0) + 0.5 * np.random.default_rng(1).standard_normal(3 * SR)
        row = kd.decays(noisy, SR, recording=True)[4]
        print(row)
        self.assertEqual(row["mark"], "nm")
        self.assertIsNone(row["decay_ms"])
        # A render is noiseless: the same test is not applied to it.
        self.assertNotEqual(kd.decays(noisy, SR)[4]["mark"], "nm")


class Wav(unittest.TestCase):
    def write(self, path, payload, tag, ch, bits, extensible=False):
        if extensible:
            fmt = struct.pack("<HHIIHHHHI", 0xFFFE, ch, SR, SR * ch * bits // 8, ch * bits // 8, bits, 22, bits, 0)
            fmt += struct.pack("<H", tag) + b"\x00\x00\x00\x00\x10\x00\x80\x00\x00\xaa\x00\x38\x9b\x71"
        else:
            fmt = struct.pack("<HHIIHH", tag, ch, SR, SR * ch * bits // 8, ch * bits // 8, bits)
        body = b"WAVE" + b"fmt " + struct.pack("<I", len(fmt)) + fmt + b"data" + struct.pack("<I", len(payload)) + payload
        Path(path).write_bytes(b"RIFF" + struct.pack("<I", len(body)) + body)

    def test_float_pcm16_pcm24_and_extensible_stereo_read_as_the_channel_mean(self):
        left = np.array([0.5, -0.25, 0.125, 0.0])
        right = np.array([0.25, 0.25, -0.125, 0.5])
        inter = np.stack([left, right], axis=1)
        mean = (left + right) / 2
        with tempfile.TemporaryDirectory() as d:
            cases = {
                "f32.wav": (inter.astype("<f4").tobytes(), 3, 32, False),
                "f32x.wav": (inter.astype("<f4").tobytes(), 3, 32, True),
                "i16.wav": (np.round(inter * 32768).astype("<i2").tobytes(), 1, 16, False),
            }
            q24 = np.round(inter * (1 << 23)).astype(np.int64).ravel()
            cases["i24.wav"] = (b"".join(int(v & 0xFFFFFF).to_bytes(3, "little") for v in q24), 1, 24, False)
            for name, (payload, tag, bits, ext) in cases.items():
                p = Path(d) / name
                self.write(p, payload, tag, 2, bits, ext)
                x, sr = kd.read_wav(p)
                self.assertEqual(sr, SR, name)
                np.testing.assert_allclose(x, mean, atol=1e-6, err_msg=name)


if __name__ == "__main__":
    unittest.main()
