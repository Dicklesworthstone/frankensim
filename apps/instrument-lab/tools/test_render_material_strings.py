"""Source-to-native launch tests; physical correctness is owned by the Rust tests.

python3 -m unittest discover -s apps/instrument-lab/tools -p test_render_material_strings.py
"""
import json
import base64
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


    def test_material_stop_model_sources_and_explicit_native_paths_reach_the_complete_experiment(self):
        root = Path(tempfile.mkdtemp())
        original = (EXAMPLES / "pluck.material-strings").read_text()
        source = root / "source.material-strings"
        source.write_text(original)
        sources = {"wire.fsmatpk": b"original-wire", "body.fsmatpk": b"original-body",
                   "contacts/receiver.fsmatpk": b"receiver-bulk\x00\xff",
                   "contacts/normal.fsintpk": b"ordered-interface\x00\xff",
                   "contacts/normal.fsmodelpk": b"original-normal-model\x00\xff"}
        for name, data in sources.items():
            path = root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        stop = root / "fingers.stops"
        stop.write_bytes(b"frankensim-string-stops-v1\r\n"
                         b"material-stop finger 0 0.25 positive fixed -1e-9 contacts/receiver.fsmatpk "
                         b"contacts/normal.fsintpk contacts/normal.fsmodelpk "
                         b"5000 5000 0.1 0.00005 0.00005 0.01 0.02 0.03 10 1e-8\r\n")
        for index, explicit in enumerate([False, True]):
            profile = stop.read_bytes()
            args = []
            expected = dict(sources)
            if explicit:
                # The native loader permits arbitrary declared paths; suffixes
                # discover conventional sources but cannot restrict an explicit binding.
                profile = profile.replace(b"contacts/normal.fsmodelpk", b"contacts/normal.model-data")
                replacement = root / "replacement.pack"
                replacement.write_bytes(b"substituted-source-model\x00\xff")
                args = ["--pack", f"contacts/normal.model-data={replacement}"]
                del expected["contacts/normal.fsmodelpk"]
                expected["contacts/normal.model-data"] = replacement.read_bytes()
            chosen = root / f"fingers-{index}.stops"
            chosen.write_bytes(profile)
            output = root / f"prepared-{index}"
            self.assertEqual(R.main(["--prepare-only", "--performance", str(source), "--out", str(output),
                                     "--stop-profile", str(chosen), *args]), 0)
            packet = json.loads((output / "experiment.native.json").read_text())
            attached = {file["path"]: base64.b64decode(file["bytes"]) for file in packet["attachments"]}
            self.assertEqual(packet["stringProfiles"]["stopPath"], "stops.profile")
            self.assertEqual(attached.pop("stops.profile"), profile)
            self.assertEqual(attached, expected)
            self.assertEqual(base64.b64decode(packet["source"]), original.encode())

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

    def test_prepare_only_exports_complete_profiles_and_original_gas_sources_without_an_engine(self):
        root = Path(tempfile.mkdtemp())
        original = (EXAMPLES / "pluck.material-strings").read_text()
        source = root / "source.material-strings"
        source.write_text(original)
        packs = {"wire.fsmatpk": b"original-wire\x00\xff", "body.fsmatpk": b"original-body",
                 "contact.fsmatpk": b"original contact material"}
        for name, data in packs.items():
            (root / name).write_bytes(data)
        body, stop, thermal = (root / name for name in ("body.source", "stop.source", "heat.source"))
        body.write_bytes(b"original body geometry\r\n# \xc2\xb5\r\n")
        stop.write_bytes(b"original moving bow profile contact.fsmatpk\r\n")
        thermal.write_bytes(b"original contact thermal history\r\n")
        gas_files = []
        for i in range(6):
            path = root / f"gas-source-{i}"
            path.write_bytes(bytes([0, 255, i, 13, 10]))
            gas_files.append(path)
        gas = "thermochemical 307 82000 mole 2 oxygen-source 0.3 inert-source 0.7"
        output = root / "prepared"
        args = ["--prepare-only", "--performance", str(source), "--out", str(output),
                "--body-profile", str(body), "--stop-profile", str(stop), "--contact-thermal", str(thermal),
                "--gas-record", gas, "--gas-component", *map(str, gas_files[:3]),
                "--gas-component", *map(str, gas_files[3:]), "--gas-water", "none", "--gas-phase", "frozen"]
        with mock.patch.object(R.subprocess, "run", side_effect=AssertionError("preparation must not render audio")):
            self.assertEqual(R.main(args), 0)
        packet = json.loads((output / "experiment.native.json").read_text())
        self.assertEqual(packet["format"], "frankensim-native-wind-input-v1")
        self.assertEqual(packet["renderer"], "strings")
        self.assertFalse(packet["decimate"])
        prepared = base64.b64decode(packet["source"])
        self.assertEqual(prepared, (output / "performance.material-strings").read_bytes())
        self.assertIn("gas " + gas, prepared.decode())
        self.assertEqual(R.row(prepared.decode().splitlines(), "string"), R.row(original.splitlines(), "string"))
        files = {file["path"]: base64.b64decode(file["bytes"]) for file in packet["attachments"]}
        self.assertEqual(packet["stringProfiles"], {"bodyPath": "body.profile", "stopPath": "stops.profile",
                                                  "contactThermalPath": "contact-thermal.profile"})
        for name, data in packs.items():
            self.assertEqual(files[name], data)
        for name, path in (("body.profile", body), ("stops.profile", stop), ("contact-thermal.profile", thermal)):
            self.assertEqual(files[name], path.read_bytes())
        self.assertEqual(packet["gasSources"]["phase"], "frozen")
        self.assertIsNone(packet["gasSources"]["waterSpecies"])
        actual_gases = [files[name] for component in packet["gasSources"]["components"] for name in component.values()]
        self.assertEqual(actual_gases, [path.read_bytes() for path in gas_files])
        self.assertFalse((output / "render.wav").exists())
        self.assertEqual(source.read_text(), original)

    def test_violin_material_bank_export_keeps_original_sections_and_explicit_decimation(self):
        root = Path(tempfile.mkdtemp())
        source = EXAMPLES / "cannone-dominant-g3.material-strings"
        output = root / "violin"
        self.assertEqual(R.main(["--prepare-only", "--performance", str(source), "--out", str(output),
                                 "--pack-root", str(EXAMPLES / "packs"), "--decimate"]), 0)
        packet = json.loads((output / "experiment.native.json").read_text())
        self.assertTrue(packet["decimate"])
        # This tests the committed four-string materials, not full specimen geometry.
        self.assertEqual(base64.b64decode(packet["source"]), source.read_bytes())
        for file in packet["attachments"]:
            self.assertEqual(base64.b64decode(file["bytes"]), (EXAMPLES / "packs" / file["path"]).read_bytes())
        self.assertGreaterEqual(len(packet["attachments"]), 5)

    def test_lab_export_requires_explicit_source_policies_and_never_overwrites_an_experiment(self):
        root = Path(tempfile.mkdtemp())
        source = root / "performance.material-strings"
        source.write_bytes(b"frankensim-material-strings-v3\r\n# \xc2\xb5\r\n")
        path = R.write_lab_experiment(source, root, [], {}, None, True)
        self.assertEqual(base64.b64decode(json.loads(path.read_text())["source"]), source.read_bytes())
        with self.assertRaises(FileExistsError):
            R.write_lab_experiment(source, root, [], {}, None, False)
        for policy in (["--gas-component", "species", "model", "transport"],
                       ["--gas-water", "none", "--gas-phase", "frozen"]):
            with self.assertRaises(SystemExit):
                R.main(["--prepare-only", "--performance", str(source), "--out", str(root / "refused"), *policy])
        self.assertFalse((root / "refused").exists())

    def test_fiber_contact_and_native_decimation_keep_the_original_mechanical_clock(self):
        root = Path(tempfile.mkdtemp())
        executable = root / "material_strings"
        executable.write_text("""#!/usr/bin/env python3
import json, sys
from pathlib import Path
p = Path(sys.argv[1])
stop = Path(sys.argv[sys.argv.index('--stop-profile') + 1])
print(json.dumps({'performance': p.read_text(), 'stop': stop.read_bytes().hex(),
    'packs': {name: (p.parent / name).read_bytes().hex()
              for name in ('wire.fsmatpk', 'body.fsmatpk', 'fiber.fsmatpk')},
    'args': sys.argv[3:]}))
Path(sys.argv[2]).write_bytes(b'render witness')
""")
        executable.chmod(0o755)
        # This is the existing manufactured contact fixture, not measured hair
        # or a claim of realistic sound. The subprocess witnesses transport only.
        original = (EXAMPLES / "fiber-contact.material-strings").read_text().replace(
            "audio 48000 4096 ", "audio 144000 12288 ")
        source = root / "source.material-strings"
        source.write_text(original)
        stop_bytes = (EXAMPLES / "fiber-contact.stops").read_bytes().replace(b"\n", b"\r\n")
        stops = root / "contacts.stops"
        stops.write_bytes(stop_bytes)
        packs = {name: ("original " + name).encode()
                 for name in ("wire.fsmatpk", "body.fsmatpk", "fiber.fsmatpk")}
        for name, data in packs.items():
            (root / name).write_bytes(data)
        self.assertEqual(R.main(["--renderer", str(executable), "--performance", str(source),
                                 "--stop-profile", str(stops), "--decimate", "--out", str(root / "run")]), 0)
        witness = json.loads((root / "run/render.log").read_text())
        self.assertEqual(witness["performance"], original)
        self.assertEqual(witness["stop"], stop_bytes.hex())
        self.assertEqual((root / "run/stops.profile").read_bytes(), stop_bytes)
        self.assertEqual(witness["packs"], {name: data.hex() for name, data in packs.items()})
        self.assertEqual(witness["args"], ["--stop-profile", str(root / "run/stops.profile"), "--decimate"])
        self.assertEqual(source.read_text(), original)
        self.assertEqual(stops.read_bytes(), stop_bytes)


if __name__ == "__main__":
    unittest.main()
