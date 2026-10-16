"""G0 guards and G2 audit reproduction (dB bands are grader-only).

uv run --with numpy --with matplotlib python3 -I test_instrgrade.py -v
Reproduce requires IRP_GRADE_SCRATCH/grade/library.json, IRP_UIOWA_WAV48,
and IRP_CONVERSION_CMD. The clarinet case additionally requires
IRP_CLARINET_BORE_WAV, IRP_CLARINET_RENDER_SHA and IRP_CLARINET_RENDER_COMMAND.
Missing external inputs fail explicitly; no skipped test counts as evidence.
The performance fixture transcribes Appendix C's bore/pressure with the
existing reed fixture parameters. Its historical attack timing is unretained;
only a successful native reproduction can establish the audit correspondence.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

TOOLS = Path(__file__).resolve().parent
ROOT = TOOLS.parents[3]
sys.path.insert(0, str(TOOLS))
import instrgrade as ig
import radiate


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def required_env(name):
    value = os.environ.get(name)
    if not value:
        raise AssertionError(f"IRPGRADE_INPUT_MISSING: set {name}")
    return value


def required_file(path):
    path = Path(path)
    if not path.is_file():
        raise AssertionError(f"IRPGRADE_INPUT_MISSING: {path}")
    return path


def aligned(path):
    x = ig.highpass(ig.load(path))
    return x[ig.onset(x):]


def run_script(name, *args, env=None):
    return subprocess.run([sys.executable, "-I", str(TOOLS / name), *map(str, args)],
                          cwd=ROOT, capture_output=True, text=True, env=env, timeout=60)


class NoAbsolutePaths(unittest.TestCase):
    def test_no_absolute_paths(self):
        """G0: zero machine-local path literals in the committed tool set."""
        tracked = subprocess.run(["git", "ls-files", "--cached", "--others", "--exclude-standard", "--",
                                  str(TOOLS.relative_to(ROOT))], cwd=ROOT, check=True,
                                 capture_output=True, text=True, timeout=10).stdout.splitlines()
        self.assertTrue(tracked)
        forbidden = [b"/" + part + b"/" for part in (b"Users", b"home", b"Volumes")]
        forbidden += [b"Down" + b"loads", b"scratch" + b"pad"]
        for name in tracked:
            content = (ROOT / name).read_bytes()
            for token in forbidden:
                self.assertNotIn(token, content, name)


class Identity(unittest.TestCase):
    def test_pack_descriptor_code_is_unchanged(self):
        """G0: the version constant is the sole change to the pinned grader."""
        original = (TOOLS / "instrgrade.py").read_bytes().replace(b'GRADER_VERSION = "instrgrade-1"\n', b"", 1)
        self.assertEqual(hashlib.sha256(original).hexdigest(),
                         "880ef40539b19bbde71b8b600dd9590957fca0f89499ed33f5154b14346b29d9")
        self.assertEqual(sha256(TOOLS / "uiowa-links.json"),
                         "46b50e792be079a2b118a3adeb632ef719d2b785f3750a6b71ea2ed14960e74d")

    def test_distance_self_is_zero(self):
        """G0: spec self-distance 0; three distance components symmetric to 1e-12 dB."""
        t = np.arange(2 * ig.SR) / ig.SR
        x = sum(np.sin(2 * np.pi * 220 * k * t) / k for k in range(1, 13))
        y = sum(np.sin(2 * np.pi * 220 * k * t) / k ** 2 for k in range(1, 13))
        a, b = ig.describe(x, f0_hint=220), ig.describe(y, f0_hint=220)
        self.assertEqual(ig.distance(a, a)["spec_db"], 0.0)
        ab, ba = ig.distance(a, b), ig.distance(b, a)
        for key in ("spec_db", "env_db", "harm_db"):
            self.assertAlmostEqual(ab[key], ba[key], delta=1e-12)


class Guards(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Retain tiny test outputs; no shared-tree cleanup or external audio copies.
        cls.work = Path(tempfile.mkdtemp(prefix="instrgrade-test-"))

    def test_radiated_candidate_refused(self):
        """G0: flag or suffix refuses with exit 2 before opening audio or references."""
        for index, (wav, flag) in enumerate((("candidate.wav", True), ("candidate-radiated.wav", False))):
            with self.subTest(wav=wav, radiated=flag):
                candidates = self.work / f"guard-{index}.json"
                candidates.write_text(json.dumps([dict(id="guard", wav=wav, radiated=flag)]))
                output = self.work / f"refused-{index}"
                result = run_script("grade_run.py", self.work, candidates, output)
                self.assertEqual(result.returncode, 2, result.stderr)
                refusal = json.loads(result.stderr.strip().splitlines()[-1])
                self.assertEqual(refusal["error"], "IRPGRADE_RADIATED_CANDIDATE_REFUSED")
                self.assertIn("DIAGNOSTIC-ONLY", refusal["message"])
                self.assertFalse(output.exists())

    def test_results_and_receipt_record_pass_and_failed_reproduction(self):
        """G0: real CLI computes identity, binds hashes, and fails a 1 dB mismatch."""
        scratch = self.work / "receipt"
        wavdir = scratch / "real" / "wav48"
        wavdir.mkdir(parents=True)
        (scratch / "grade").mkdir()
        t = np.arange(2 * ig.SR) / ig.SR
        wav = wavdir / "Synthetic.mf.A3.wav"
        ig.write_wav(wav, sum(np.sin(2 * np.pi * 220 * k * t) / (20 * k) for k in range(1, 13)))
        desc = ig.describe(aligned(wav), f0_hint=220)
        library = [dict(file=wav.name, start=0, end=2, note=57, name="A3", kind="sustained",
                        inst="synthetic", dyn="mf", art="", desc=desc)]
        libpath = scratch / "grade" / "library.json"
        libpath.write_text(json.dumps(library))
        env = dict(os.environ, IRP_UIOWA_WAV48=str(wavdir), LIBFILE="library.json")
        for index, expected in enumerate((0.0, 1.0)):
            candidates = scratch / f"candidates-{index}.json"
            candidates.write_text(json.dumps([dict(id="identity", label="synthetic identity control",
                wav=str(wav), inst="synthetic", note="A3", dyn="mf", radiated=False,
                expected_db=expected, tolerance_db=1e-12, conversion_cmd="synthetic float WAV; no conversion")]))
            output = scratch / f"output-{index}"
            result = run_script("grade_run.py", scratch, candidates, output, env=env)
            self.assertEqual(result.returncode, index, result.stderr)
            row = json.loads((output / "results.json").read_text())[0]
            receipt, terminal = map(json.loads, (output / "receipt.jsonl").read_text().splitlines())
            self.assertEqual(row["grader_version"], "instrgrade-1")
            self.assertEqual(receipt["candidate_wav_sha256"], sha256(wav))
            self.assertEqual(receipt["library_sha256"], sha256(libpath))
            self.assertEqual(receipt["tool_sha256"]["instrgrade.py"], sha256(TOOLS / "instrgrade.py"))
            self.assertEqual(receipt["reference_file"], wav.name)
            self.assertEqual(receipt["reference_note"], "A3")
            self.assertEqual(receipt["reference_dynamic"], "mf")
            self.assertEqual(receipt["reference_wav_sha256"], sha256(wav))
            self.assertAlmostEqual(receipt["spec_db"], 0.0, delta=1e-12)
            for field in ("env_db", "harm_db", "rulers", "expected_db", "tolerance_db", "conversion_cmd",
                          "elapsed_s", "host", "git_sha", "command"):
                self.assertIn(field, receipt)
            self.assertEqual(terminal["event"], "terminal")
            self.assertEqual(terminal["verdict"], "pass" if index == 0 else "fail")
        diagnostic_dir = scratch / "diagnostic"
        diagnostic = run_script("radiate.py", libpath, wav, diagnostic_dir,
                                "--instrument", "synthetic", "--note", "A3",
                                "--segments", '[["cyl", 0.0073, 0.57]]', "--flanged")
        self.assertEqual(diagnostic.returncode, 0, diagnostic.stderr)
        self.assertEqual(diagnostic.stdout.splitlines()[0], "DIAGNOSTIC-ONLY: never a graded render")
        diagnostic_row = json.loads(diagnostic.stdout.splitlines()[-1])
        self.assertIs(diagnostic_row["radiated"], True)
        self.assertIs(diagnostic_row["diagnostic_only"], True)
        self.assertEqual(diagnostic_row["grader_version"], "instrgrade-1")
        self.assertTrue((diagnostic_dir / "Synthetic.mf.A3-radiated.wav").is_file())


class Reproduce(unittest.TestCase):
    def inputs(self, note):
        scratch = Path(required_env("IRP_GRADE_SCRATCH"))
        wavdir = Path(required_env("IRP_UIOWA_WAV48"))
        required_env("IRP_CONVERSION_CMD")
        libpath = required_file(scratch / "grade" / os.environ.get("LIBFILE", "library.json"))
        library = json.loads(libpath.read_text())
        targets = ig.find(library, "clarinet", ig.midi_of(note), dyn="mf")
        self.assertTrue(targets, f"IRPGRADE_REFERENCE_MISSING: clarinet {note} mf")
        target = targets[0]
        required_file(wavdir / target["file"])
        live = ig.describe(ig.ref_audio(target, wavdir), f0_hint=ig.hz_of(ig.midi_of(note)))
        self.assertLessEqual(ig.distance(live, target["desc"])["spec_db"], 0.1,
                             "IRPGRADE_REFERENCE_DRIFT: converted audio differs from the cached descriptor")
        return scratch, wavdir, target, library

    def test_reed_sustained_row(self):
        """G2: original tracked WAV vs E3 mf; |delta| <= 0.1 dB, no Rust build."""
        scratch, wavdir, target, library = self.inputs("E3")
        wav = required_file(ROOT / "data" / "listening" / "reed-sustained.wav")
        desc = ig.describe(aligned(wav), f0_hint=ig.hz_of(ig.midi_of("E3")))
        observed = ig.distance(desc, target["desc"])["spec_db"]
        self.assertAlmostEqual(observed, 23.5407, delta=0.1)
        print(json.dumps(dict(grader_version=ig.GRADER_VERSION, test=self.id(), spec_db=observed,
                              expected_db=23.5407, tolerance_db=0.1, candidate_wav_sha256=sha256(wav))))

    def test_clarinet_bore_row(self):
        """G2: native D3 mf row 28.6 +/-0.1; separate radiation diagnostic 12.6 +/-0.1."""
        scratch, wavdir, target, library = self.inputs("D3")
        wav = required_file(required_env("IRP_CLARINET_BORE_WAV"))
        render_sha = required_env("IRP_CLARINET_RENDER_SHA")
        self.assertRegex(render_sha, r"^[0-9a-f]{40}$", "IRPGRADE_RENDER_SHA_INVALID")
        required_env("IRP_CLARINET_RENDER_COMMAND")
        x = aligned(wav)
        hint = ig.hz_of(ig.midi_of("D3"))
        raw = ig.describe(x, f0_hint=hint)
        grade = ig.distance(raw, target["desc"])["spec_db"]
        self.assertAlmostEqual(grade, 28.6, delta=0.1)
        # The historical screen used 146.8 Hz as its hint; retain it exactly.
        y, _, _ = radiate.apply(x, [("cyl", 0.0073, 0.57)], True)
        diagnostic = ig.distance(ig.describe(y, f0_hint=146.8), target["desc"])["spec_db"]
        self.assertAlmostEqual(diagnostic, 12.6, delta=0.1)
        print(json.dumps(dict(grader_version=ig.GRADER_VERSION, test=self.id(), spec_db=grade,
                              diagnostic_only_spec_db=diagnostic, diagnostic_only=True,
                              expected_db=28.6, diagnostic_expected_db=12.6, tolerance_db=0.1,
                              render_git_sha=render_sha, candidate_wav_sha256=sha256(wav))))


if __name__ == "__main__":
    unittest.main()
