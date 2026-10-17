#!/usr/bin/env python3
"""Unit and integration tests for the instrgrade -> fs-psycho subprocess bridge (M7.5-W1)."""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from instrgrade import SR, write_wav  # noqa: E402
from psycho import (  # noqa: E402
    PSYCHO_METRICS,
    compare_psycho_blocks,
    evaluate_psycho_signal,
    psycho_bin_sha256,
    resolve_psycho_bin,
)

REPO_ROOT = Path(__file__).resolve().parents[4]
DEFAULT_BIN = REPO_ROOT / "target" / "debug" / "examples" / "psycho_descriptors"


def get_test_bin() -> Path:
    env_bin = resolve_psycho_bin()
    if env_bin is not None:
        return env_bin
    if DEFAULT_BIN.is_file() and os.access(DEFAULT_BIN, os.X_OK):
        return DEFAULT_BIN
    raise unittest.SkipTest(
        f"psycho_descriptors binary not found at {DEFAULT_BIN}; build with "
        "`cargo build -p fs-psycho --example psycho_descriptors`"
    )


class TestPsychoBridge(unittest.TestCase):
    def test_iso_1khz_40db_one_sone(self) -> None:
        """G2: 1 kHz 40 dB SPL pure tone at 48 kHz free-field yields 1.0 +/- 0.15 sone."""
        bin_path = get_test_bin()
        n = SR  # 1.0 s at 48 kHz
        t = np.arange(n, dtype=np.float64) / SR
        # 40 dB SPL -> p_rms = 20 uPa * 10^(40/20) = 2e-3 Pa
        p_rms_40 = 20e-6 * (10.0 ** (40.0 / 20.0))
        x_40 = (math.sqrt(2.0) * p_rms_40 * np.sin(2.0 * math.pi * 1000.0 * t)).astype(np.float32)

        res_40 = evaluate_psycho_signal(
            x_40,
            pa_per_full_scale=1.0,
            calibration_source="sidecar",
            field="free",
            bin_override=bin_path,
        )
        self.assertEqual(res_40["schema"], "fs-psycho-descriptors/1")
        self.assertFalse(res_40["level_unverified"])
        sone_40 = res_40["loudness_sone"]
        self.assertIsNotNone(sone_40)
        self.assertAlmostEqual(sone_40, 1.0, delta=0.15)

        # 60 dB SPL -> p_rms = 0.02 Pa; ISO 532-1 Annex B.3 tone = 4.019 sone, DIN sharpness ~ 1.0 acum
        p_rms_60 = 20e-6 * (10.0 ** (60.0 / 20.0))
        x_60 = (math.sqrt(2.0) * p_rms_60 * np.sin(2.0 * math.pi * 1000.0 * t)).astype(np.float32)
        res_60 = evaluate_psycho_signal(
            x_60,
            pa_per_full_scale=1.0,
            calibration_source="sidecar",
            field="free",
            bin_override=bin_path,
        )
        self.assertAlmostEqual(res_60["loudness_sone"], 4.019, delta=0.05)
        self.assertGreaterEqual(res_60["sharpness_acum"], 0.7)
        self.assertLessEqual(res_60["sharpness_acum"], 1.3)

    def test_am_tone_roughness_peak_near_70hz(self) -> None:
        """G2: 100% AM 1 kHz 60 dB SPL sweep peaks at 70 Hz and matches MoSQITo reference."""
        bin_path = get_test_bin()
        n = 8192
        t = np.arange(n, dtype=np.float64) / SR
        p_rms = 20e-6 * (10.0 ** (60.0 / 20.0))
        amp = math.sqrt(2.0) * p_rms

        # Pinned reference values from crates/fs-psycho/tests/psycho.rs::roughness_am_tone_peaks_near_70_hz
        ref_curve = [
            (20.0, 0.197_327_796_978_835_65),
            (40.0, 0.660_982_441_631_880_9),
            (55.0, 0.949_126_590_618_143_1),
            (70.0, 1.044_791_270_606_281_9),
            (90.0, 0.940_671_646_615_856_8),
            (120.0, 0.641_113_651_282_113_5),
            (160.0, 0.335_153_389_958_774_8),
        ]
        observed: list[tuple[float, float]] = []
        for fmod, expected_asper in ref_curve:
            sig = (amp * (1.0 + np.sin(2.0 * math.pi * fmod * t)) * np.sin(2.0 * math.pi * 1000.0 * t)).astype(
                np.float32
            )
            res = evaluate_psycho_signal(
                sig,
                pa_per_full_scale=1.0,
                calibration_source="sidecar",
                field="free",
                bin_override=bin_path,
            )
            r_val = res["roughness_asper"]
            self.assertIsNotNone(r_val, f"roughness_asper refused at {fmod} Hz: {res.get('refusals')}")
            rel_err = abs(r_val - expected_asper) / expected_asper
            self.assertLess(rel_err, 1e-5, f"fmod={fmod} Hz: got {r_val}, expected {expected_asper}")
            observed.append((fmod, float(r_val)))

        peak_fmod, peak_asper = max(observed, key=lambda kv: kv[1])
        self.assertGreaterEqual(peak_fmod, 55.0)
        self.assertLessEqual(peak_fmod, 90.0)
        self.assertEqual(peak_fmod, 70.0)
        self.assertAlmostEqual(peak_asper, 1.04479127, places=5)

    def test_missing_binary_is_typed_refusal(self) -> None:
        """G0: Missing FS_PSYCHO_BIN records typed refusal and None metrics (never 0 or NaN)."""
        sig = np.zeros(SR, dtype=np.float32)
        res = evaluate_psycho_signal(
            sig,
            pa_per_full_scale=1.0,
            calibration_source="sidecar",
            field="free",
            bin_override="/nonexistent/path/to/psycho_descriptors",
        )
        for m in PSYCHO_METRICS:
            self.assertIsNone(res[m], f"metric {m} must be None on missing binary, got {res[m]!r}")
            self.assertIn(m, res["refusals"])
            self.assertTrue(res["refusals"][m].startswith("FS_PSYCHO_BIN_MISSING:"))

    def test_reference_marked_level_unverified(self) -> None:
        """G0: Reference with assumed full-scale marks level-dependent comparisons level_unverified."""
        bin_path = get_test_bin()
        t = np.arange(SR, dtype=np.float64) / SR
        env = np.minimum(1.0, t / 0.05) * np.exp(-1.5 * t)
        cand_sig = (0.25 * env * np.sin(2.0 * math.pi * 440.0 * t)).astype(np.float32)
        ref_sig = (0.20 * env * np.sin(2.0 * math.pi * 440.0 * t)).astype(np.float32)
        ruler_sig = (0.30 * env * np.sin(2.0 * math.pi * 440.0 * t)).astype(np.float32)

        ours = evaluate_psycho_signal(
            cand_sig,
            pa_per_full_scale=2.0,
            calibration_source="sidecar",
            field="free",
            bin_override=bin_path,
        )
        real = evaluate_psycho_signal(
            ref_sig,
            pa_per_full_scale=1.0,
            calibration_source="assumed",
            field="free",
            bin_override=bin_path,
        )
        ruler = evaluate_psycho_signal(
            ruler_sig,
            pa_per_full_scale=1.0,
            calibration_source="assumed",
            field="free",
            bin_override=bin_path,
        )
        self.assertFalse(ours["level_unverified"])
        self.assertTrue(real["level_unverified"])

        cmp_block = compare_psycho_blocks(ours, real, ruler_blocks=[ruler])
        self.assertTrue(cmp_block["level_unverified"])
        self.assertTrue(cmp_block["metrics"]["loudness_sone"]["level_unverified"])
        self.assertTrue(cmp_block["metrics"]["roughness_asper"]["level_unverified"])
        self.assertFalse(cmp_block["metrics"]["log_attack_time"]["level_unverified"])
        self.assertIsNotNone(cmp_block["metrics"]["loudness_sone"]["inside_real_spread"])

        # End-to-end check through grade_run.py
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            wavdir = root / "real" / "wav48"
            gradedir = root / "grade"
            outdir = root / "out"
            wavdir.mkdir(parents=True)
            gradedir.mkdir(parents=True)

            write_wav(wavdir / "ref.wav", ref_sig)
            cand_wav = root / "cand.wav"
            write_wav(cand_wav, cand_sig)

            from instrgrade import describe, hz_of, midi_of

            f0 = hz_of(midi_of("A4"))
            ref_desc = describe(ref_sig, "sustained", f0)
            lib = [
                {
                    "inst": "flute",
                    "art": "vib",
                    "dyn": "mf",
                    "kind": "sustained",
                    "name": "A4",
                    "note": midi_of("A4"),
                    "f0": f0,
                    "file": "ref.wav",
                    "start": 0.0,
                    "end": 1.0,
                    "desc": ref_desc,
                }
            ]
            (gradedir / "library.json").write_text(json.dumps(lib))
            cands = [
                {
                    "id": "flute_a4_psycho",
                    "label": "Flute A4 psycho test",
                    "wav": str(cand_wav),
                    "inst": "flute",
                    "note": "A4",
                    "dyn": "mf",
                    "art": "vib",
                    "kind": "sustained",
                    "pa_per_full_scale": 2.0,
                    "calibration_source": "sidecar",
                }
            ]
            cands_path = root / "candidates.json"
            cands_path.write_text(json.dumps(cands))

            # The temporary scratch holds grade/library.json and real/wav48; the
            # Reproduce rows' inherited LIBFILE and IRP_UIOWA_WAV48 must not leak.
            env_vars = dict(os.environ, FS_PSYCHO_BIN=str(bin_path), LIBFILE="library.json",
                            IRP_UIOWA_WAV48=str(root / "real" / "wav48"))
            grade_run_py = Path(__file__).resolve().parent / "grade_run.py"
            proc = subprocess.run(
                [sys.executable, "-I", str(grade_run_py), str(root), str(cands_path), str(outdir)],
                capture_output=True,
                text=True,
                env=env_vars,
                timeout=60,
            )
            self.assertEqual(proc.returncode, 0, msg=proc.stderr or proc.stdout)
            results = json.loads((outdir / "results.json").read_text())
            self.assertEqual(len(results), 1)
            self.assertEqual(results[0]["psycho_bin_sha256"], psycho_bin_sha256(bin_path))
            self.assertTrue(results[0]["psycho"]["level_unverified"])
            self.assertIsNotNone(results[0]["psycho"]["ours"]["loudness_sone"])
            receipts = [
                json.loads(line)
                for line in (outdir / "receipt.jsonl").read_text().splitlines()
                if line.strip()
            ]
            cand_receipt = next(r for r in receipts if r["event"] == "candidate")
            self.assertEqual(cand_receipt["psycho_bin_sha256"], psycho_bin_sha256(bin_path))
            self.assertTrue(cand_receipt["psycho"]["level_unverified"])


if __name__ == "__main__":
    unittest.main()
