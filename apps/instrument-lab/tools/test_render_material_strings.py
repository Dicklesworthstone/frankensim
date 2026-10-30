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

    def test_independent_cavity_gas_preserves_mechanics_and_original_record_order(self):
        original = (EXAMPLES / "pluck.material-strings").read_text()
        exterior = "mixture 301 101325 mass 3 N2 0.6 O2 0.3 CO2 0.1"
        cavity = "thermochemical 290 101325 mole 2 helium-source 0.8 oxygen-source 0.2"
        gas_line = next(line for line in original.splitlines(keepends=True) if line.startswith("gas "))
        inherited = R.gas_records(original, exterior, None)
        self.assertNotIn("cavity_gas", inherited)
        inserted = R.gas_records(original, None, cavity)
        self.assertEqual(inserted, original.replace(gas_line, gas_line + "cavity_gas " + cavity + "\n", 1))
        independent = R.gas_records(inserted, exterior, None)
        self.assertIn("cavity_gas " + cavity + "\n", independent)
        self.assertEqual(R.gas_records(inserted, None, None), inserted)
        replacement = "mixture 288 101325 mole 1 He 1"
        replaced = R.gas_records(independent, None, replacement)
        self.assertEqual(replaced.count("cavity_gas "), 1)
        self.assertIn("gas " + exterior + "\ncavity_gas " + replacement + "\n", replaced)
        self.assertEqual([line for line in replaced.splitlines() if line.split()[0] not in ("gas", "cavity_gas")],
                         [line for line in original.splitlines() if line.split()[0] != "gas"])
        # The original molecular declaration follows the inserted medium;
        # its law and budgets must never be recalculated or replaced here.
        law = "cavity_molecular 32 1024 original-law-data\n"
        molecular = original.replace(gas_line, gas_line + law, 1)
        self.assertIn(gas_line + "cavity_gas " + cavity + "\n" + law,
                      R.gas_records(molecular, None, cavity))

    def test_gas_record_overrides_refuse_multiple_or_ambiguous_records(self):
        original = (EXAMPLES / "pluck.material-strings").read_text()
        for record in ("", "   ", "mixture\nstate 1 1", "mixture\rstate 1 1"):
            for exterior, cavity in ((record, None), (None, record)):
                with self.assertRaisesRegex(ValueError, "one nonempty record"):
                    R.gas_records(original, exterior, cavity)
        valid = R.gas_records(original, None, "mixture 293.15 101325 mole 1 He 1")
        for source in (valid.replace("body body.fsmatpk", "cavity_gas mixture 293.15 101325 mole 1 He 1\nbody body.fsmatpk"),
                       original.replace("strings 1", "cavity_gas mixture 293.15 101325 mole 1 He 1\nstrings 1")):
            with self.assertRaisesRegex(ValueError, "single record immediately after gas"):
                R.gas_records(source, None, "mixture 293.15 101325 mole 1 N2 1")

    def test_prepare_only_normalizes_crlf_to_the_native_canonical_schedule(self):
        root = Path(tempfile.mkdtemp())
        source = root / "source.material-strings"
        original = (EXAMPLES / "pluck.material-strings").read_bytes().replace(b"\n", b"\r\n")
        source.write_bytes(original)
        for name in ("wire.fsmatpk", "body.fsmatpk"):
            (root / name).write_bytes(b"original material source")
        output = root / "prepared"
        self.assertEqual(R.main(["--prepare-only", "--performance", str(source), "--out", str(output)]), 0)
        canonical = original.replace(b"\r\n", b"\n")
        self.assertEqual((output / "performance.material-strings").read_bytes(), canonical)
        packet = json.loads((output / "experiment.native.json").read_text())
        self.assertEqual(base64.b64decode(packet["source"]), canonical)
        self.assertIn(b"\nschedule\n", canonical)
        self.assertEqual(source.read_bytes(), original)

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


    def test_source_manifests_reach_measured_packets_and_native_render_without_rewriting_setup(self):
        root = Path(tempfile.mkdtemp()).resolve()
        original = (EXAMPLES / "measured-cello-g.material-strings").read_text()
        # A native source path has no required filename extension.
        original = original.replace("cello-g-whole.fsmatpk", "strings/whole.source")
        source = root / "source.material-strings"
        source.write_text(original)
        (root / "strings").mkdir()
        stale = root / "strings/whole.source"
        stale.write_bytes(b"stale pack must never win an explicit compilation")
        body = root / "replacement body.fsmatpk"
        body.write_bytes(b"original body source bytes\x00\xff")
        source_dir = root / "source tables (original)"
        source_dir.mkdir()
        manifest = source_dir / "whole string.manifest.tsv"
        source_pack = REPO / "data/matdb/seed-v1/cello-g-strings-whole-lampis-2024"
        manifest.write_bytes((source_pack / "manifest.tsv").read_bytes()
                             .replace(b"citation\t", "citation\tµ ".encode()).replace(b"\n", b"\r\n"))
        table = source_dir / "properties.tsv"
        table.write_bytes((source_pack / "properties.tsv").read_bytes().replace(b"\n", b"\r\n"))
        before = (manifest.read_bytes(), table.read_bytes())
        compiler, renderer = root / "source compiler", root / "native renderer"
        emitted = b"compiler transport witness\x00\xff"
        calls = []

        def native_call(command, **kwargs):
            calls.append(command)
            if command[0] == str(compiler):
                source_args = command[-2:] if portable else []
                if portable:
                    self.assertEqual(source_args[0], "--source-inputs")
                    self.assertEqual(Path(source_args[1]).parent.parent, output)
                self.assertEqual(command, [str(compiler), "matdb-pack", "--manifest",
                                           str(manifest), "--out", str(output / "strings/whole.source"),
                                           *(["--material-state", *identity] if identity else []), *source_args])
                self.assertEqual(kwargs, {"check": True})
                self.assertEqual((manifest.read_bytes(), table.read_bytes()), before)
                # This witness is only an executable boundary. The real Rust
                # compiler owns manifest decoding and normalized pack semantics.
                Path(command[command.index("--out") + 1]).write_bytes(emitted)
                if portable:
                    Path(source_args[1]).write_text(json.dumps({
                        "format": "frankensim-material-source-inputs-v1",
                        "manifest": before[0].decode("utf-8"),
                        "sources": [{"path": "properties.tsv", "text": before[1].decode("utf-8")}]}))
            else:
                self.assertEqual(command, [str(renderer), str(output / "performance.material-strings"),
                                           str(output / "render.wav")])
                self.assertEqual(Path(command[1]).read_text(), original)
                self.assertEqual((output / "strings/whole.source").read_bytes(), emitted)
                self.assertEqual((output / "body.fsmatpk").read_bytes(), body.read_bytes())
                Path(command[2]).write_bytes(b"native renderer transport witness")
            return mock.Mock(returncode=0)

        with mock.patch.object(R.shutil, "which", side_effect=lambda name: name), \
             mock.patch.object(R.subprocess, "run", side_effect=native_call):
            card_identity = ["declared whole-string", "solid", "source specimen", "4294967295"]
            for name, mode, identity in (
                    ("prepared", ["--prepare-only"], []),
                    ("rendered", ["--renderer", str(renderer), "--lab-experiment"], []),
                    ("card-prepared", ["--prepare-only"], card_identity),
                    ("card-rendered", ["--renderer", str(renderer), "--lab-experiment"], card_identity),
                    ("native-only", ["--renderer", str(renderer)], [])):
                output = root / name
                portable = "--lab-experiment" in mode or "--prepare-only" in mode
                self.assertEqual(R.main(["--performance", str(source), "--out", str(output),
                    "--pack-compiler", str(compiler),
                    "--pack-manifest", f"strings/whole.source={manifest}",
                    "--pack", f"body.fsmatpk={body}",
                    *(["--pack-material-state", "strings/whole.source", *identity] if identity else []),
                    *mode]), 0)
                if not portable:
                    self.assertFalse((output / "experiment.native.json").exists())
                    continue
                packet = json.loads((output / "experiment.native.json").read_text())
                self.assertEqual(base64.b64decode(packet["source"]), original.encode())
                self.assertEqual({item["path"]: base64.b64decode(item["bytes"])
                                  for item in packet["attachments"]},
                                 {"body.fsmatpk": body.read_bytes(),
                                  "material-sources-0/manifest.tsv": before[0],
                                  "material-sources-0/properties.tsv": before[1]})
                binding = {"packPath": "strings/whole.source",
                           "manifestPath": "material-sources-0/manifest.tsv"}
                if identity:
                    binding["materialState"] = dict(zip(("chemistry", "phase", "process", "revision"),
                                                        [*identity[:3], int(identity[3])]))
                self.assertEqual(packet["materialManifests"], [binding])
                self.assertEqual(packet["renderer"], "strings")
        self.assertEqual(len(calls), 8)  # five compiles, three renders; native-only calls retain their old arguments
        self.assertEqual((manifest.read_bytes(), table.read_bytes()), before)
        self.assertEqual(source.read_text(), original)
        self.assertEqual(stale.read_bytes(), b"stale pack must never win an explicit compilation")

    def test_editable_material_export_preserves_nested_source_paths_and_avoids_existing_names(self):
        root = Path(tempfile.mkdtemp())
        source = root / "performance.material-strings"
        source.write_bytes(b"frankensim-material-strings-v3\r\n")
        original = {"wire.fsmatpk": b"native compiled pack", "material-sources-0": b"other original pack"}
        for name, data in original.items():
            (root / name).write_bytes(data)
        inputs = {"wire.fsmatpk": {"manifest": "original µ manifest\r\n", "sources": [
            {"path": "manifest.tsv", "text": "original source named manifest.tsv\r\n"},
            {"path": "tables/wire.tsv", "text": "original µ source table\r\n"}]}}
        packet = json.loads(R.write_lab_experiment(source, root, list(original), {}, None, False,
                                                   material_inputs=inputs).read_text())
        files = {item["path"]: base64.b64decode(item["bytes"]) for item in packet["attachments"]}
        prefix = "material-sources-0-input/"
        self.assertEqual(files, {"material-sources-0": original["material-sources-0"],
            prefix + "manifest.tsv": b"original source named manifest.tsv\r\n",
            prefix + "tables/wire.tsv": "original µ source table\r\n".encode(),
            prefix + "manifest-input.tsv": "original µ manifest\r\n".encode()})
        self.assertEqual(packet["materialManifests"], [{"packPath": "wire.fsmatpk",
            "manifestPath": prefix + "manifest-input.tsv"}])
        self.assertEqual((root / "wire.fsmatpk").read_bytes(), original["wire.fsmatpk"])

    def test_editable_material_export_counts_generated_packs_and_refuses_unportable_inputs(self):
        root = Path(tempfile.mkdtemp())
        source = root / "performance.material-strings"
        source.write_bytes(b"frankensim-material-strings-v3\n")
        (root / "wire.fsmatpk").write_bytes(b"native compiled pack")
        input_for = lambda path: {"wire.fsmatpk": {"manifest": "original manifest",
            "sources": [{"path": path, "text": "original source"}]}}
        for path in ("../outside.tsv", "/outside.tsv", "with space.tsv", "nested/./source.tsv"):
            with self.assertRaises(ValueError):
                R.write_lab_experiment(source, root, ["wire.fsmatpk"], {}, None, False,
                                       material_inputs=input_for(path))
        large = b"x" * (8 * 1024 * 1024)
        # Actual byte budgets, without writing large fixture files to disk.
        small_source = source.read_bytes()
        names = ["wire.fsmatpk", "body.fsmatpk", "bow.fsmatpk", "contact.fsmatpk"]
        with mock.patch.object(R.Path, "read_bytes",
                               side_effect=lambda path=None: large) as reader:
            reader.side_effect = [*(large for _ in names), small_source]
            with self.assertRaisesRegex(ValueError, "byte limits"):
                R.write_lab_experiment(source, root, names, {}, None, False,
                                       material_inputs=input_for("properties.tsv"))
        many = {"wire.fsmatpk": {"manifest": "original manifest", "sources":
            [{"path": f"table-{i}.tsv", "text": "source"} for i in range(31)]}}
        with self.assertRaisesRegex(ValueError, "source-file"):
            R.write_lab_experiment(source, root, ["wire.fsmatpk"], {}, None, False, material_inputs=many)
        self.assertFalse((root / "experiment.native.json").exists())

    def test_manifest_bindings_refuse_conflicts_missing_sources_and_failed_compilation(self):
        root = Path(tempfile.mkdtemp())
        (root / "wire.fsmatpk").write_bytes(b"stale wire")
        output = root / "new"
        output.mkdir()
        manifest = root / "source.tsv"
        manifest.write_bytes(b"original manifest")
        binding = f"wire.fsmatpk={manifest}"
        cases = [
            (["string wire.fsmatpk"], [binding], [binding]),
            (["string wire.fsmatpk"], [], [binding, binding]),
            (["string wire.fsmatpk"], [], [f"unused.fsmatpk={manifest}"]),
            (["string ../wire.fsmatpk"], [], [f"../wire.fsmatpk={manifest}"]),
            (["string wire.fsmatpk"], [], ["wire.fsmatpk="]),
            (["string wire.fsmatpk", "body missing.fsmatpk"], [], [binding]),
        ]
        with mock.patch.object(R.subprocess, "run", side_effect=AssertionError("compiler before binding admission")):
            for texts, packs, manifests in cases:
                with self.assertRaises(ValueError):
                    R.stage_packs(texts, root, packs, output, manifests, "xtask")
        self.assertEqual(list(output.iterdir()), [])
        failure = R.subprocess.CalledProcessError(1, ["xtask", "matdb-pack"])
        with mock.patch.object(R.subprocess, "run", side_effect=failure):
            with self.assertRaises(R.subprocess.CalledProcessError):
                R.stage_packs(["string wire.fsmatpk"], root, [], output, [binding], "xtask")
        self.assertFalse((output / "wire.fsmatpk").exists())
        self.assertEqual((root / "wire.fsmatpk").read_bytes(), b"stale wire")
        with mock.patch.object(R.subprocess, "run", return_value=mock.Mock(returncode=0)):
            with self.assertRaisesRegex(ValueError, "did not create a nonempty pack"):
                R.stage_packs(["string wire.fsmatpk"], root, [], output, [binding], "xtask")

    def test_material_state_bindings_refuse_ignored_or_ambiguous_identities_before_io(self):
        identity = ["wire.fsmatpk", "Cu", "solid", "annealed", "0"]
        manifests = ["wire.fsmatpk=source.tsv"]
        self.assertEqual(R.material_state_bindings([identity], manifests),
                         {"wire.fsmatpk": identity[1:]})
        for bindings, sources in [
            ([identity], []),
            ([identity, identity], manifests),
            ([["other.pack", *identity[1:]]], manifests),
            ([[*identity[:-1], "-1"]], manifests),
            ([[*identity[:-1], "4294967296"]], manifests),
            ([[*identity[:-1], "01"]], manifests),
            ([["wire.fsmatpk", "", "solid", "annealed", "0"]], manifests),
            ([identity[:-1]], manifests),
        ]:
            with self.assertRaises(ValueError):
                R.material_state_bindings(bindings, sources)
        for options in (
                ["--pack-material-state", *identity],
                ["--pack-compiler", "xtask", "--pack-manifest", manifests[0],
                 "--pack-material-state", *identity, "--pack-material-state", *identity],
                ["--pack-compiler", "xtask", "--pack-manifest", manifests[0],
                 "--pack-material-state", *identity[:-1], "4294967296"]):
            with mock.patch.object(R.Path, "read_text", side_effect=AssertionError("source I/O before refusal")), \
                 mock.patch.object(R.Path, "mkdir", side_effect=AssertionError("directory before refusal")), \
                 mock.patch.object(R.subprocess, "run", side_effect=AssertionError("compiler before refusal")):
                with self.assertRaises(SystemExit):
                    R.main(["--prepare-only", "--performance", "unread.material-strings",
                            "--out", "uncreated", *options])

    def test_manifest_compiler_options_refuse_before_source_io(self):
        for options in (["--pack-manifest", "wire.fsmatpk=source.tsv"],
                        ["--pack-compiler", "missing-xtask"],
                        ["--pack-compiler", "missing-xtask", "--pack-manifest", "wire.fsmatpk=source.tsv"]):
            with mock.patch.object(R.shutil, "which", return_value=None), \
                 mock.patch.object(R.Path, "read_text", side_effect=AssertionError("source I/O before refusal")), \
                 mock.patch.object(R.Path, "mkdir", side_effect=AssertionError("directory before refusal")), \
                 mock.patch.object(R.subprocess, "run", side_effect=AssertionError("compiler before refusal")):
                with self.assertRaises(SystemExit):
                    R.main(["--prepare-only", "--performance", "unread.material-strings",
                            "--out", "uncreated", *options])

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
        cavity = "mixture 300 50662.5 mole 3 N2 0.6 O2 0.2 CO2 0.2"
        thermal = root / "source.contact-thermal"
        thermal.write_bytes(b"frankensim-contact-thermal-v1\r\nboundary adiabatic\r\n")
        body = root / "source.body"
        body.write_bytes(b"original coated shell profile\r\n")
        args = ["--renderer", str(executable), "--performance", str(source), "--out", str(root / "run"),
                "--gas-record", gas, "--cavity-gas-record", cavity,
                "--pack", f"wire.fsmatpk={other}", "--full-scale-pa", "5",
                "--contact-thermal", str(thermal), "--body-profile", str(body),
                "--coating-strain-guard", "536870912", "1000000000000"]
        self.assertEqual(R.main(args), 0)
        witness = json.loads((root / "run/render.log").read_text())
        self.assertIn("gas " + gas, witness["performance"])
        self.assertIn("gas " + gas + "\ncavity_gas " + cavity + "\n", witness["performance"])
        self.assertEqual(witness["wire"], "different-density-and-stiffness")
        self.assertEqual(witness["thermal"], thermal.read_bytes().hex())
        self.assertEqual((root / "run/contact-thermal.profile").read_bytes(), thermal.read_bytes())
        self.assertEqual(witness["args"], ["--body-profile", str(root / "run/body.profile"),
            "--contact-thermal", str(root / "run/contact-thermal.profile"),
            "--coating-strain-guard", "536870912", "1000000000000"])
        self.assertEqual(R.row(witness["performance"].splitlines(), "string")[1], R.row(original.splitlines(), "string")[1])
        self.assertEqual(R.row(witness["performance"].splitlines(), "state")[1], R.row(original.splitlines(), "state")[1])
        self.assertEqual(source.read_text(), original)
        with self.assertRaises(FileExistsError):
            R.main(args)

    def test_portable_coating_guard_refuses_invalid_or_inexact_budgets_before_io(self):
        root = Path(tempfile.mkdtemp())
        # JavaScript must retain the native budgets exactly through a saved packet.
        for index, budget in enumerate((["0", "2"], ["1", "9007199254740992"])):
            output = root / f"refused-{index}"
            with self.assertRaises(SystemExit):
                R.main(["--prepare-only", "--performance", str(root / "unread"),
                        "--out", str(output), "--body-profile", str(root / "unread-body"),
                        "--coating-strain-guard", *budget])
            self.assertFalse(output.exists())

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
        cavity = "thermochemical 290 82000 mole 2 oxygen-source 0.8 inert-source 0.2"
        output = root / "prepared"
        args = ["--prepare-only", "--performance", str(source), "--out", str(output),
                "--body-profile", str(body), "--stop-profile", str(stop), "--contact-thermal", str(thermal),
                "--coating-strain-guard", "536870912", "1000000000000",
                "--gas-record", gas, "--cavity-gas-record", cavity,
                "--gas-component", *map(str, gas_files[:3]),
                "--gas-component", *map(str, gas_files[3:]), "--gas-water", "none", "--gas-phase", "frozen"]
        with mock.patch.object(R.subprocess, "run", side_effect=AssertionError("preparation must not render audio")):
            self.assertEqual(R.main(args), 0)
        packet = json.loads((output / "experiment.native.json").read_text())
        self.assertEqual(packet["format"], "frankensim-native-wind-input-v1")
        self.assertEqual(packet["renderer"], "strings")
        self.assertFalse(packet["decimate"])
        self.assertEqual(packet["coatingStrainGuard"], {"maximumBytes": 536870912, "maximumWork": 1000000000000})
        prepared = base64.b64decode(packet["source"])
        self.assertEqual(prepared, (output / "performance.material-strings").read_bytes())
        self.assertIn("gas " + gas, prepared.decode())
        self.assertIn("gas " + gas + "\ncavity_gas " + cavity + "\n", prepared.decode())
        self.assertEqual(R.row(prepared.decode().splitlines(), "string")[1], R.row(original.splitlines(), "string")[1])
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


class LoadedResponses(unittest.TestCase):
    def test_coherent_arguments_preserve_physical_amplitude_phase_and_declaration_order(self):
        sweep = ["101", "2000", "32", "coherent"]
        forces = [["1", "0", "-2.5e-3"], ["0", "0.25", "0"], ["1", "0", "2.5e-3"]]
        expected = [*sweep,
                    "--body-response-force", "1", "0", "-2.5e-3",
                    "--body-response-force", "0", "0.25", "0",
                    "--body-response-force", "1", "0", "2.5e-3",
                    "--body-response-strain", "1048576", "1000000"]
        self.assertEqual(R.loaded_response_arguments(sweep, forces, [1048576, 1000000]), expected)
        self.assertEqual(forces[0], ["1", "0", "-2.5e-3"])
        self.assertEqual(R.loaded_response_arguments(None, [], None), [])
        for bridge in ("0", "all"):
            query = ["101", "2000", "32", bridge]
            self.assertEqual(R.loaded_response_arguments(query, [], None), query)
        for query, loads, strain in [
            (None, [["0", "1", "0"]], None),
            (None, [], [1, 2]),
            (sweep, [], None),
            (["101", "2000", "32", "all"], forces, None),
            (sweep, forces * 22, None),
            (sweep, forces, [0, 2]),
        ]:
            with self.assertRaises(ValueError):
                R.loaded_response_arguments(query, loads, strain)

    def test_harmonic_playback_and_export_conflicts_refuse_before_source_io(self):
        base = ["--renderer", "native-strings", "--performance", "unread.material-strings",
                "--out", "uncreated", "--body-profile", "unread.body",
                "--loaded-response-only", "101", "2000", "8", "coherent",
                "--body-response-force", "0", "-1e-3", "2e-3"]
        for conflict in (["--decimate"], ["--contact-thermal", "unread.thermal"],
                         ["--stop-profile", "unread.stops"],
                         ["--coating-strain-guard", "1", "2"], ["--prepare-only"],
                         ["--lab-experiment"], ["--analysis", "uncreated-analysis"]):
            with mock.patch.object(R.Path, "read_text", side_effect=AssertionError("source I/O before refusal")), \
                 mock.patch.object(R.Path, "read_bytes", side_effect=AssertionError("source I/O before refusal")), \
                 mock.patch.object(R.Path, "mkdir", side_effect=AssertionError("directory before refusal")), \
                 mock.patch.object(R.subprocess, "run", side_effect=AssertionError("native call before refusal")):
                with self.assertRaises(SystemExit):
                    R.main([*base, *conflict])

    def test_loaded_native_launch_stages_original_sources_and_emits_no_wav(self):
        root = Path(tempfile.mkdtemp())
        executable = root / "material_strings"
        # This is an I/O witness only. No harmonic or acoustic value is fabricated.
        executable.write_text("""#!/usr/bin/env python3
import json, sys
from pathlib import Path
assert sys.argv[2] == '--loaded-body-response-only'
source = Path(sys.argv[1])
body = Path(sys.argv[sys.argv.index('--body-profile') + 1])
print(json.dumps({'performance': source.read_text(), 'body': body.read_bytes().hex(),
                  'wire': (source.parent / 'wire.fsmatpk').read_bytes().hex(),
                  'args': sys.argv[2:]}))
Path(sys.argv[3]).write_text('transport witness only\\n')
""")
        executable.chmod(0o755)
        original = (EXAMPLES / "pluck.material-strings").read_text()
        source = root / "source.material-strings"
        source.write_text(original)
        (root / "wire.fsmatpk").write_bytes(b"original wire")
        (root / "body.fsmatpk").write_bytes(b"original body")
        replacement = root / "other.fsmatpk"
        replacement.write_bytes(b"replacement wire with its own source identity")
        profile = root / "source.body"
        profile.write_bytes(b"original shell profile\r\n")
        output = root / "loaded"
        sweep = ["101", "2000", "8", "coherent"]
        args = ["--renderer", str(executable), "--performance", str(source), "--out", str(output),
                "--body-profile", str(profile), "--pack", f"wire.fsmatpk={replacement}",
                "--loaded-response-only", *sweep,
                "--body-response-force", "0", "-1e-3", "2e-3",
                "--body-response-force", "0", "0.0005", "0",
                "--body-response-strain", "1048576", "1000000000"]
        self.assertEqual(R.main(args), 0)
        witness = json.loads((output / "loaded-response.log").read_text())
        self.assertEqual(witness["performance"], original)
        self.assertEqual(witness["body"], profile.read_bytes().hex())
        self.assertEqual(witness["wire"], replacement.read_bytes().hex())
        self.assertEqual(witness["args"], [
            "--loaded-body-response-only", str(output / "loaded-response.tsv"), *sweep,
            "--body-response-force", "0", "-1e-3", "2e-3",
            "--body-response-force", "0", "0.0005", "0",
            "--body-response-strain", "1048576", "1000000000",
            "--body-profile", str(output / "body.profile")])
        self.assertEqual((output / "loaded-response.tsv").read_text(), "transport witness only\n")
        self.assertFalse((output / "render.wav").exists())
        self.assertFalse((output / "experiment.native.json").exists())
        self.assertEqual(source.read_text(), original)
        with self.assertRaises(FileExistsError):
            R.main(args)


if __name__ == "__main__":
    unittest.main()
