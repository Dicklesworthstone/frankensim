"""Source-to-native launch tests; physical correctness is owned by the Rust tests.

python3 -m unittest discover -s apps/instrument-lab/tools -p test_render_material_strings.py
"""
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import render_material_strings as R

REPO = Path(__file__).resolve().parents[3]
EXAMPLES = REPO / "crates/fs-couple/examples/material_strings"
CARDS = REPO / "data/instruments/descriptions/violin"


class MaterialPreparation(unittest.TestCase):
    def test_every_open_note_retains_the_complete_material_bank_and_only_rebinds_the_bow(self):
        original = (EXAMPLES / "cannone-dominant-g3.material-strings").read_text()
        original_header = original.split("\nschedule\n")[0].splitlines()
        for index, note in enumerate(("G3", "D4", "A4", "E5")):
            card = (CARDS / f"arco-open-{note}.performance").read_text()
            text = R.bow_schedule(original, card, index, 0.1)
            header, schedule = text.split("\nschedule\n")
            # The material density, stiffness, loss, torsion, geometry, tension,
            # body, cavity and numerical rung remain original data for all notes.
            self.assertEqual([r for r in header.splitlines() if not r.startswith("audio ")],
                             [r for r in original_header if not r.startswith("audio ")])
            self.assertIn("audio 144000 432000 1 128", header)
            self.assertIn(f"track\tbow\tbow\t{index}", schedule)
            commands = [r.split("\t") for r in schedule.splitlines() if r.startswith(("initial\t", "event\t"))]
            old = [r.split("\t") for r in card.splitlines() if r.startswith(("initial\t", "event\t"))]
            self.assertEqual([r[:-1] for r in commands], [r[:-1] for r in old])
            self.assertTrue(all(r[-1] == "7.707571801566581e-1" for r in commands))

    def test_assembly_dimensions_come_from_the_native_conversion(self):
        original = (EXAMPLES / "cannone-dominant-g3.material-strings").read_text()
        report = "body-length-m 0.4\nbody-width-m 0.2\nreference-thickness-m 0.003\nbridge-snap-m 0 0.001\ncavity-volume-m3 0.002\ncavity-neck-radius-m 0.03\ncavity-neck-length-m 0.004\n"
        changed = R.assembly_geometry(original, report)
        self.assertIn("body violin-spruce-norway-viala2020-shell.fsmatpk 0.4 0.2 0.003 30 1", changed)
        self.assertIn("cavity physical-neck 0.002 0.03 0.004 flanged flanged 1", changed)
        self.assertEqual(changed.split("\nschedule\n")[1], original.split("\nschedule\n")[1])
        with self.assertRaises(KeyError):
            R.assembly_geometry(original, "body-length-m 0.4\n")

    def test_pack_replacement_is_bytes_only_and_refuses_unconsumed_or_duplicate_names(self):
        root = Path(tempfile.mkdtemp())
        packs = root / "packs"
        packs.mkdir()
        (packs / "wire.fsmatpk").write_bytes(b"original pack bytes")
        replacement = root / "replacement.fsmatpk"
        replacement.write_bytes(b"replacement pack with its own source identity")
        output = root / "variant"
        output.mkdir()
        R.stage_packs(["string wire.fsmatpk 0.3"], packs, [f"wire.fsmatpk={replacement}"], output)
        self.assertEqual((output / "wire.fsmatpk").read_bytes(), replacement.read_bytes())
        self.assertEqual((packs / "wire.fsmatpk").read_bytes(), b"original pack bytes")
        for bindings in ([f"unknown.fsmatpk={replacement}"], [f"wire.fsmatpk={replacement}"] * 2):
            with self.assertRaises(ValueError):
                R.stage_packs(["string wire.fsmatpk 0.3"], packs, bindings, output)

    def test_native_launch_consumes_the_substituted_gas_pack_and_original_tension(self):
        root = Path(tempfile.mkdtemp())
        # This executable is an I/O witness, not a mock physical model. It checks
        # the actual staged files received by the subprocess and records them.
        executable = root / "material_strings"
        executable.write_text("#!/usr/bin/env python3\nimport json, sys\nfrom pathlib import Path\np=Path(sys.argv[1])\nt=Path(sys.argv[sys.argv.index('--contact-thermal')+1])\nprint(json.dumps({'performance':p.read_text(),'wire':(p.parent/'wire.fsmatpk').read_text(),'args':sys.argv[3:],'thermal':t.read_bytes().hex()}))\nPath(sys.argv[2]).write_bytes(b'render witness')\n")
        executable.chmod(0o755)
        original = (EXAMPLES / "pluck.material-strings").read_text()
        source = root / "source.material-strings"
        source.write_text(original)
        (root / "wire.fsmatpk").write_bytes(b"original-wire")
        (root / "body.fsmatpk").write_bytes(b"original-body")
        other = root / "other.fsmatpk"
        other.write_bytes(b"different-density-and-stiffness")
        gas = "mixture 293.15 50662.5 mole 3 N2 0.7 O2 0.2 He 0.1"
        thermal = root / "source.contact-thermal"
        thermal.write_bytes(b"frankensim-contact-thermal-v1\r\nboundary adiabatic\r\n")
        args = ["--renderer", str(executable), "--performance", str(source), "--out", str(root / "run"),
                "--gas-record", gas, "--pack", f"wire.fsmatpk={other}", "--full-scale-pa", "5",
                "--contact-thermal", str(thermal)]
        self.assertEqual(R.main(args), 0)
        witness = json.loads((root / "run/render.log").read_text())
        self.assertIn("gas " + gas, witness["performance"])
        self.assertEqual(witness["wire"], "different-density-and-stiffness")
        self.assertEqual(witness["thermal"], thermal.read_bytes().hex())
        self.assertEqual((root / "run/contact-thermal.profile").read_bytes(), thermal.read_bytes())
        self.assertEqual(witness["args"], ["--contact-thermal", str(root / "run/contact-thermal.profile")])
        self.assertEqual(R.row(witness["performance"].splitlines(), "string"), R.row(original.splitlines(), "string"))
        self.assertEqual(source.read_text(), original)
        with self.assertRaises(FileExistsError):
            R.main(args)


if __name__ == "__main__":
    unittest.main()
