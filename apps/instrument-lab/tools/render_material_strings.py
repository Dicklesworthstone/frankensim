#!/usr/bin/env python3
"""Render or analyze an existing material-string assembly through its native solver.

All mechanics, material admission, gas resolution and PCM encoding stay in the
native material_strings executable. This launcher stages its original inputs;
it neither estimates missing properties nor tunes a string after substitution.
Explicit --pack-manifest bindings invoke the existing Rust matdb-pack compiler;
its output enters the native solve, while portable Lab experiments retain the
original editable manifests and source tables for recompilation.
Contact profiles reach that same native solve. Optional native decimation changes
only the observed audio clock, retaining the authored mechanical resolution.
The output directory must be new. Python's standard library suffices for native
rendering and harmonic analysis; --analysis additionally uses the existing
instrgrade dependencies for listening exports.
"""
from __future__ import annotations

import argparse
import base64
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from decimal import Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent


def row(lines: list[str], key: str) -> tuple[int, list[str]]:
    matches = [(i, line.split()) for i, line in enumerate(lines) if line.startswith(key + " ")]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one {key} record")
    return matches[0]


def gas_records(text: str, exterior: str | None, cavity: str | None) -> str:
    """Replace only explicitly selected native media, without evaluating either.

    An absent cavity record inherits the exterior in the native solver. An
    existing one remains independent, including when only the exterior changes.
    The launcher supplies canonical LF text required by the native schedule.
    Keep the original material state, molecular law and mechanics.
    """
    for option, value in (("--gas-record", exterior), ("--cavity-gas-record", cavity)):
        if value is not None and (not value.strip() or "\n" in value or "\r" in value):
            raise ValueError(f"{option} must be one nonempty record")
    if exterior is None and cavity is None:
        return text
    header, schedule = text.split("\nschedule\n", 1)
    lines = header.splitlines()
    keys = [line.split()[0] if line.split() else "" for line in lines]
    gas_rows = [i for i, key in enumerate(keys) if key == "gas"]
    cavity_rows = [i for i, key in enumerate(keys) if key == "cavity_gas"]
    if len(gas_rows) != 1:
        raise ValueError("expected exactly one gas record")
    gas_row = gas_rows[0]
    if len(cavity_rows) > 1 or (cavity_rows and cavity_rows[0] != gas_row + 1):
        raise ValueError("cavity_gas must be a single record immediately after gas")

    if exterior is not None:
        lines[gas_row] = "gas " + exterior
    if cavity is not None:
        if cavity_rows:
            index = cavity_rows[0]
            lines[index] = "cavity_gas " + cavity
        else:
            lines.insert(gas_row + 1, "cavity_gas " + cavity)
    return "\n".join(lines) + "\nschedule\n" + schedule


def rust_e(value: float) -> str:
    """Canonical gesture f64 exponent spelling, as used by the existing note generator."""
    if not math.isfinite(value):
        raise ValueError("nonfinite gesture value")
    if value == 0:
        return "0e0"
    sign, digits, exponent = Decimal(repr(value)).as_tuple()
    digits = list(digits)
    while len(digits) > 1 and digits[-1] == 0:
        digits.pop()
        exponent += 1
    mantissa = str(digits[0]) + ("." + "".join(map(str, digits[1:])) if len(digits) > 1 else "")
    return ("-" if sign else "") + f"{mantissa}e{exponent + len(digits) - 1}"


def bow_schedule(text: str, card: str, string: int, distance_from_bridge: float) -> str:
    """Use only a card's physical bow commands/horizon, never its old string or body.

    The requested distance is a fraction of the speaking span between the nut
    and the material string's bridge. No length, tension, mass or loss is changed.
    A fixed bow station is explicit; moving-station source cards refuse.
    """
    header, _ = text.split("\nschedule\n", 1)
    old_header, schedule = card.split("\nschedule\n", 1)
    lines = header.splitlines()
    strings = [line.split() for line in lines if line.startswith("string ")]
    if not 0 <= string < len(strings) or not 0 < distance_from_bridge < 1:
        raise ValueError("bow string index or bridge-relative station is outside the supplied geometry")
    fields = strings[string]
    # v1 has a scalar tension; v2/v3 explicitly select tension or extension.
    bridge_index = 7 if lines[0].endswith("-v1") else 9 if fields[4] == "extension" else 8
    at = rust_e(float(fields[bridge_index]) * (1 - distance_from_bridge))
    records = [line.split("\t") for line in schedule.splitlines()]
    if [r for r in records if r[0] == "tracks"] != [["tracks", "1"]]:
        raise ValueError("a bow card must contain exactly one canonical bow track")
    positions = set()
    for r in records:
        if r[0] == "track":
            if len(r) != 4 or r[2] != "bow":
                raise ValueError("source card has a non-bow track")
            r[3] = str(string)
        elif r[0] in ("initial", "event"):
            kind = 1 if r[0] == "initial" else 3
            if len(r) != kind + 4 or r[kind] != "bow":
                raise ValueError("source card has a non-bow command")
            positions.add(float(r[-1]))
            r[-1] = at
    if len(positions) != 1:
        raise ValueError("moving-station bow cards need an explicit material-string performance")
    i, audio = row(lines, "audio")
    _, old_audio = row(old_header.splitlines(), "audio")
    samples, remainder = divmod(int(old_audio[2]) * int(audio[1]), int(old_audio[1]))
    if remainder:
        raise ValueError("source bow duration is not an integer number of target-clock samples")
    audio[2] = str(samples)
    lines[i] = " ".join(audio)
    return "\n".join(lines) + "\nschedule\n" + "\n".join("\t".join(r) for r in records) + "\n"


def assembly_geometry(text: str, report: str) -> str:
    """Bind the native shell converter's dimensions, with no parallel geometry calculation."""
    values = {r[0]: r[1] for line in report.splitlines() if len(r := line.split()) == 2}
    header, schedule = text.split("\nschedule\n", 1)
    lines = header.splitlines()
    i, body = row(lines, "body")
    body[2:5] = [values[k] for k in ("body-length-m", "body-width-m", "reference-thickness-m")]
    lines[i] = " ".join(body)
    i, cavity = row(lines, "cavity")
    if cavity[1] != "none":
        cavity[2:5] = [values[k] for k in ("cavity-volume-m3", "cavity-neck-radius-m", "cavity-neck-length-m")]
        lines[i] = " ".join(cavity)
    return "\n".join(lines) + "\nschedule\n" + schedule


def material_state_bindings(bindings: list[list[str]], manifests: list[str]) -> dict[str, list[str]]:
    """Explicit original card identities; never infer chemistry or a contact pair.

    The Rust compiler owns material-card admission. Here we only prevent an
    ignored, ambiguous or malformed CLI binding before staging any sources.
    """
    targets = {binding.partition("=")[0] for binding in manifests}
    states = {}
    for binding in bindings:
        if len(binding) != 5:
            raise ValueError("--pack-material-state needs PACK CHEMISTRY PHASE PROCESS REVISION")
        name, *identity = binding
        if name not in targets or name in states:
            raise ValueError("--pack-material-state must name one distinct --pack-manifest target")
        if any(not part.strip() or part.startswith("--")
               or any(ord(char) < 32 or ord(char) == 127 for char in part)
               for part in identity[:3]):
            raise ValueError("material chemistry, phase and process must be nonempty identity fields")
        revision = identity[3]
        if not re.fullmatch(r"0|[1-9][0-9]*", revision) or int(revision) > 4_294_967_295:
            raise ValueError("material revision must be a canonical unsigned 32-bit integer")
        states[name] = identity
    return states


def stage_packs(texts: list[str], root: Path, replacements: list[str], output: Path,
                manifests: list[str] | None = None, compiler: str | None = None,
                material_states: list[list[str]] | None = None,
                source_inputs: dict | None = None) -> list[str]:
    # These are the existing binary pack formats, copied byte for byte.
    # Rust remains the material decoder and property/domain authority.
    states = material_state_bindings(material_states or [], manifests or [])
    suffixes = (".fsmatpk", ".fsintpk", ".fsmodelpk")
    available = {str(p.relative_to(root)): p for suffix in suffixes for p in root.rglob("*" + suffix)}
    tokens = {token for text in texts for line in text.splitlines() if not line.startswith("#") for token in line.split()}
    if manifests and compiler is None:
        raise ValueError("--pack-manifest requires --pack-compiler")
    bound = set()
    compiled = {}
    for option, bindings, sources in (("--pack", replacements, available),
                                      ("--pack-manifest", manifests or [], compiled)):
        for binding in bindings:
            name, sep, source = binding.partition("=")
            # Explicit bindings may use any native-referenced filename. Rust
            # compiles and decodes the actual bulk/interface/model source.
            if not sep or not name or not source or name not in tokens or name in bound:
                raise ValueError(f"{option} must name one distinct referenced path: {binding!r}")
            path = Path(name)
            if path.is_absolute() or any(p in ("", ".", "..") for p in name.split("/")):
                raise ValueError("material pack paths must stay relative to the performance directory")
            sources[name] = Path(source).expanduser().resolve()
            bound.add(name)
    required = {t for t in tokens if t.endswith(suffixes)} | bound
    for name in sorted(required):
        if name not in available and name not in compiled:
            raise ValueError(f"missing material pack {name}; supply --pack-root, --pack {name}=PATH or --pack-manifest {name}=MANIFEST")
    input_dir = Path(tempfile.mkdtemp(prefix="material-inputs-", dir=output)) if compiled and source_inputs is not None else None
    for index, name in enumerate(sorted(required)):
        destination = output / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        if name in compiled:
            if destination.exists():
                raise FileExistsError(destination)
            # The existing compiler resolves sibling source tables, source
            # identity, dimensions, uncertainty and validity from the manifest.
            # A failed compilation stops here; never use a stale available pack.
            state_args = ["--material-state", *states[name]] if name in states else []
            input_path = input_dir / f"{index}.json" if input_dir else None
            input_args = ["--source-inputs", str(input_path)] if input_path else []
            subprocess.run([compiler, "matdb-pack", "--manifest", str(compiled[name]),
                            "--out", str(destination), *state_args, *input_args], check=True)
            if not destination.is_file() or destination.stat().st_size == 0:
                raise ValueError(f"material compiler did not create a nonempty pack: {name}")
            if input_path is not None:
                with input_path.open("rb") as stream:
                    raw = stream.read(48 * 1024 * 1024 + 1)
                if len(raw) > 48 * 1024 * 1024:
                    raise ValueError("material source-input export exceeds 48 MiB")
                inputs = json.loads(raw)
                if not isinstance(inputs, dict) or inputs.get("format") != "frankensim-material-source-inputs-v1":
                    raise ValueError("compiler did not return original material source inputs")
                if name in states:
                    chemistry, phase, process, revision = states[name]
                    inputs["materialState"] = {"chemistry": chemistry, "phase": phase,
                                               "process": process, "revision": int(revision)}
                source_inputs[name] = inputs
        else:
            with destination.open("xb") as stream, available[name].open("rb") as source:
                shutil.copyfileobj(source, stream)
    return sorted(required)


def write_lab_experiment(source: Path, output: Path, attachments: list[str], profiles: dict,
                         gas_sources: dict | None, decimate: bool,
                         coating_strain_guard: list[int] | None = None,
                         material_inputs: dict | None = None) -> Path:
    """Export the existing native HTTP packet; material/physics decoding stays in Rust.

    The lab's shared packet reader checks its transport bounds when opening this
    file. Prepared source bytes and all referenced source paths remain intact.
    """
    encoded = lambda data: base64.b64encode(data).decode("ascii")
    material_inputs = material_inputs or {}
    original = {name: (output / name).read_bytes() for name in attachments}
    if len(original) != len(attachments) or not material_inputs.keys() <= original.keys():
        raise ValueError("portable material inputs must bind distinct staged pack paths")
    names = {"input.performance", "output.wav", "output.provenance.json"}
    def admit_path(name):
        if (not isinstance(name, str) or not name or len(name.encode("utf-8")) > 512
                or name.startswith("-") or "\\" in name or ":" in name
                or any(ord(char) <= 32 or ord(char) == 127 for char in name)
                or any(not part or part.startswith(".") for part in name.split("/"))
                or any(name == prior or name.startswith(prior + "/")
                       or prior.startswith(name + "/") for prior in names)):
            raise ValueError("portable source paths must be confined, distinct and free of whitespace")
        names.add(name)
        return name
    for name in original:
        admit_path(name)
    files = {name: data for name, data in original.items() if name not in material_inputs}
    bindings = []
    for index, (target, inputs) in enumerate(material_inputs.items()):
        # Only the canonical compiler discovers dependencies and interprets the
        # manifest. Relative source paths and exact UTF-8 bytes remain unchanged.
        prefix = f"material-sources-{index}"
        while any(prefix == name or prefix.startswith(name + "/")
                  or name.startswith(prefix + "/") for name in names):
            prefix += "-input"
        if not isinstance(inputs, dict):
            raise ValueError("invalid canonical material source-input export")
        source_files = inputs.get("sources")
        if not isinstance(inputs.get("manifest"), str) or not isinstance(source_files, list) or not source_files:
            raise ValueError("incomplete canonical material source-input export")
        for item in source_files:
            if not isinstance(item, dict) or not isinstance(item.get("text"), str):
                raise ValueError("invalid canonical material source file")
            relative = item.get("path")
            if not isinstance(relative, str) or relative.startswith("/"):
                raise ValueError("invalid canonical material source path")
            files[admit_path(prefix + "/" + relative)] = item["text"].encode("utf-8")
        manifest_name = prefix + "/manifest.tsv"
        while any(manifest_name == name or manifest_name.startswith(name + "/")
                  or name.startswith(manifest_name + "/") for name in names):
            manifest_name = manifest_name.removesuffix(".tsv") + "-input.tsv"
        files[admit_path(manifest_name)] = inputs["manifest"].encode("utf-8")
        bindings.append({"packPath": target, "manifestPath": manifest_name,
                         **({"materialState": inputs["materialState"]} if "materialState" in inputs else {})})
    sizes = [len(data) for data in files.values()] + [len(original[name]) for name in material_inputs]
    source_bytes = source.read_bytes()
    if (len(files) + len(bindings) > 32 or any(not size or size > 8 * 1024 * 1024 for size in sizes)
            or sum(sizes) > 32 * 1024 * 1024 or not source_bytes or len(source_bytes) > 1024 * 1024):
        raise ValueError("portable experiment exceeds the existing source-file or byte limits")
    packet = {"format": "frankensim-native-wind-input-v1", "renderer": "strings",
              "source": encoded(source_bytes), "decimate": decimate,
              "attachments": [{"path": name, "bytes": encoded(data)} for name, data in files.items()]}
    if bindings:
        packet["materialManifests"] = bindings
    if profiles:
        packet["stringProfiles"] = profiles
    if gas_sources is not None:
        packet["gasSources"] = gas_sources
    if coating_strain_guard is not None:
        packet["coatingStrainGuard"] = {"maximumBytes": coating_strain_guard[0],
                                        "maximumWork": coating_strain_guard[1]}
    path = output / "experiment.native.json"
    serialized = json.dumps(packet, separators=(",", ":"), ensure_ascii=False)
    if len(serialized.encode("utf-8")) > 4 * ((33 * 1024 * 1024 + 2) // 3) + 65536:
        raise ValueError("portable experiment exceeds the existing complete request limit")
    with path.open("x", encoding="utf-8") as stream:
        stream.write(serialized)
    return path


def loaded_response_arguments(sweep: list[str] | None, forces: list[list[str]],
                              strain: list[int] | None) -> list[str]:
    """Transport native physical inputs without evaluating or normalizing them.

    The native parser and original source own numerical bounds, actual bridge
    existence, force summation, source bands and response units. This preflight
    refuses options that would otherwise be discarded before staging any files.
    """
    if sweep is None:
        if forces or strain is not None:
            raise ValueError("body-response force/strain options require --loaded-response-only")
        return []
    coherent = sweep[3] == "coherent"
    if coherent != bool(forces):
        raise ValueError("select coherent with at least one --body-response-force; unit/all selections take no force records")
    if len(forces) > 64:
        raise ValueError("at most 64 body-response force records are supported")
    if strain is not None and min(strain) <= 0:
        raise ValueError("--body-response-strain needs positive byte and work allowances")
    arguments = list(sweep)
    for force in forces:
        arguments += ["--body-response-force", *force]
    if strain is not None:
        arguments += ["--body-response-strain", *map(str, strain)]
    return arguments


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    # argparse otherwise mistakes a signed scientific phasor such as -1e-3
    # for an option. Numerical/source admission still belongs to the native CLI.
    parser._negative_number_matcher = re.compile(r"^-(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$")

    parser.add_argument("--renderer", help="path to material_strings; required for rendering, harmonic analysis or shell assembly")
    parser.add_argument("--performance", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True, help="new directory for staged sources, native output and log")
    parser.add_argument("--pack-root", type=Path)
    parser.add_argument("--pack", action="append", default=[], metavar="REFERENCED_PATH=SOURCE_PACK")
    parser.add_argument("--pack-manifest", action="append", default=[], metavar="REFERENCED_PATH=MANIFEST",
                        help="compile this original source manifest into its referenced pack with the existing Rust material compiler")
    parser.add_argument("--pack-compiler", help="path to the existing xtask executable; required with --pack-manifest")
    parser.add_argument("--pack-material-state", action="append", nargs=5, default=[],
                        metavar=("PACK", "CHEMISTRY", "PHASE", "PROCESS", "REVISION"),
                        help="bind an explicit original material-card identity to one --pack-manifest target; omit for the original property/interface/model/species envelope")
    parser.add_argument("--body-profile", type=Path)
    parser.add_argument("--stop-profile", type=Path, help="original finger/fret/bow-surface contact profile")
    parser.add_argument("--contact-thermal", type=Path, help="original finite contact-heat profile; thermal properties come from the same assembly packs")
    parser.add_argument("--shell-table", type=Path)
    parser.add_argument("--assembly-options", type=Path)
    parser.add_argument("--gas-record", help="exact exterior/string gas record after 'gas'; the cavity inherits it only when no cavity_gas is declared")
    parser.add_argument("--cavity-gas-record", help="exact independent interior gas record after 'cavity_gas'; requires a native cavity and the same static pressure as the exterior")
    parser.add_argument("--gas-component", action="append", nargs=3, default=[], metavar=("SPECIES", "NASA9", "TRANSPORT"))
    parser.add_argument("--gas-water")
    parser.add_argument("--gas-phase", choices=("frozen", "liquid-water"))
    parser.add_argument("--full-scale-pa", type=float, help="explicit PCM scale; use the same value for every comparison")
    parser.add_argument("--bow-card", type=Path, help="reuse only this legacy card's canonical bow commands and duration")
    parser.add_argument("--string-index", type=int)
    parser.add_argument("--bow-distance-from-bridge", type=float)
    parser.add_argument("--decimate", action="store_true", help="use the native 48 kHz output decimator without changing the mechanical clock")
    parser.add_argument("--coating-strain-guard", nargs=2, type=int, metavar=("MAX_BYTES", "MAX_WORK"),
                        help="check all coating films at every native mechanical endpoint with explicit resources, retained in lab experiments")
    parser.add_argument("--loaded-response-only", nargs=4, metavar=("MIN_HZ", "MAX_HZ", "POINTS", "BRIDGE"),
                        help="save loaded-response.tsv through native released harmonic analysis; BRIDGE is an index, all or coherent")
    parser.add_argument("--body-response-force", action="append", nargs=3, default=[],
                        metavar=("BRIDGE", "REAL_N", "IMAG_N"),
                        help="physical peak force phasor for coherent loaded analysis; repeat to add actuators")
    parser.add_argument("--body-response-strain", nargs=2, type=int, metavar=("MAX_BYTES", "MAX_WORK"),
                        help="include native coating-strain bounds in the loaded harmonic response")
    parser.add_argument("--lab-experiment", action="store_true", help="also save experiment.native.json for Open experiment in Instrument Lab")
    parser.add_argument("--prepare-only", action="store_true", help="stage sources and save the lab experiment without rendering audio")
    parser.add_argument("--analysis", type=Path, help="also write standard calibrated listening pairs for build_instruments.py")
    parser.add_argument("--family")
    parser.add_argument("--note")
    parser.add_argument("--variant", default="material", help="analysis slot suffix; keeps original graded notes intact")
    parser.add_argument("--git-sha", help="actual renderer revision, required when exporting to the lab")
    parser.add_argument("--observer", help="actual receiver description, required when exporting to the lab")
    parser.add_argument("--recipe", default=str(Path(__file__).resolve()))
    args = parser.parse_args(argv)
    if bool(args.pack_manifest) != bool(args.pack_compiler):
        parser.error("--pack-manifest and --pack-compiler must be supplied together")
    try:
        material_state_bindings(args.pack_material_state, args.pack_manifest)
    except ValueError as error:
        parser.error(str(error))
    compiler = shutil.which(args.pack_compiler) if args.pack_compiler else None
    if args.pack_compiler and compiler is None:
        parser.error("--pack-compiler must name an available xtask executable")
    compiler = str(Path(compiler).resolve()) if compiler else None
    try:
        response_arguments = loaded_response_arguments(
            args.loaded_response_only, args.body_response_force, args.body_response_strain)
    except ValueError as error:
        parser.error(str(error))
    if response_arguments:
        if args.prepare_only or args.lab_experiment or args.analysis:
            parser.error("--loaded-response-only produces a native response, not an audio/lab-experiment export")
        if args.decimate or args.contact_thermal or args.stop_profile or args.coating_strain_guard:
            parser.error("loaded harmonic analysis cannot select playback decimation, nonlinear stops, contact heat or the transient coating guard")
        if not (args.body_profile or args.shell_table):
            parser.error("--loaded-response-only requires a shell body profile or shell table")
    renderer = shutil.which(args.renderer) if args.renderer else None
    if (args.renderer and renderer is None) or (renderer is None and (not args.prepare_only or args.shell_table)):
        parser.error("--renderer must name an available material_strings executable")
    renderer = str(Path(renderer).resolve()) if renderer else None
    export_experiment = args.lab_experiment or args.prepare_only
    if args.coating_strain_guard:
        if min(args.coating_strain_guard) <= 0:
            parser.error("--coating-strain-guard needs positive byte and work allowances")
        if export_experiment and max(args.coating_strain_guard) > 9_007_199_254_740_991:
            parser.error("portable coating strain budgets must be safe integers no larger than 9007199254740991")
        if not (args.body_profile or args.shell_table):
            parser.error("--coating-strain-guard requires a shell body profile or shell table")
    if args.prepare_only and args.analysis:
        parser.error("--prepare-only cannot export listening analysis without an audio render")
    if export_experiment and (args.gas_component or args.gas_water or args.gas_phase):
        if not (args.gas_component and args.gas_water and args.gas_phase):
            parser.error("a lab experiment with gas sources needs --gas-component, --gas-water and --gas-phase")
    if bool(args.shell_table) != bool(args.assembly_options) or (args.shell_table and args.body_profile):
        parser.error("supply --body-profile or both --shell-table and --assembly-options")
    if args.bow_card and (args.string_index is None or args.bow_distance_from_bridge is None):
        parser.error("--bow-card requires --string-index and --bow-distance-from-bridge")
    if args.analysis and not all((args.family, args.note, args.git_sha, args.observer)):
        parser.error("--analysis requires --family, --note, --git-sha and --observer")
    for label in (args.family, args.note, args.variant):
        if label and (Path(label).name != label or label in (".", "..")):
            parser.error("family, note and variant must be single path components")
    text = args.performance.read_text()
    if args.bow_card:
        text = bow_schedule(text, args.bow_card.read_text(), args.string_index, args.bow_distance_from_bridge)
    try:
        text = gas_records(text, args.gas_record, args.cavity_gas_record)
    except ValueError as error:
        parser.error(str(error))
    header, schedule = text.split("\nschedule\n", 1)
    lines = header.splitlines()
    if args.full_scale_pa is not None:
        if not math.isfinite(args.full_scale_pa) or args.full_scale_pa <= 0:
            parser.error("--full-scale-pa must be finite positive pascals")
        i, audio = row(lines, "audio")
        audio[3] = repr(args.full_scale_pa)
        lines[i] = " ".join(audio)
    text = "\n".join(lines) + "\nschedule\n" + schedule
    thermal_profile = args.contact_thermal.read_bytes() if args.contact_thermal else None
    stop_profile = args.stop_profile.read_bytes() if args.stop_profile else None
    output = args.out.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=False)
    profile_bytes = args.body_profile.read_bytes() if args.body_profile else b""
    profile = profile_bytes.decode("utf-8")
    profiles = {}
    attachments = []
    if args.shell_table:
        result = subprocess.run([renderer, "assemble-shell", str(args.shell_table.resolve()),
                                 str(args.assembly_options.resolve()), str(output / "body.profile")],
                                check=True, capture_output=True, text=True)
        (output / "assembly.log").write_text(result.stdout + result.stderr)
        profile = (output / "body.profile").read_bytes().decode("utf-8")
        text = assembly_geometry(text, result.stdout)
    elif profile:
        (output / "body.profile").write_bytes(profile_bytes)
    if profile:
        profiles["bodyPath"] = "body.profile"
        attachments.append("body.profile")
    if stop_profile is not None:
        (output / "stops.profile").write_bytes(stop_profile)
        profiles["stopPath"] = "stops.profile"
        attachments.append("stops.profile")
    if thermal_profile is not None:
        (output / "contact-thermal.profile").write_bytes(thermal_profile)
        profiles["contactThermalPath"] = "contact-thermal.profile"
        attachments.append("contact-thermal.profile")
    # Thermal cells refer to the SAME assembly packs. Do not introduce an
    # independent thermal pack loader or bypass the native card-identity check.
    material_inputs = {} if export_experiment else None
    attachments += stage_packs([text, profile, (stop_profile or b"").decode("utf-8")],
                              args.pack_root or args.performance.parent, args.pack, output,
                              args.pack_manifest, compiler, args.pack_material_state, material_inputs)
    source = output / "performance.material-strings"
    source.write_text(text)
    gas_components = []
    for index, paths in enumerate(args.gas_component):
        component = {}
        for key, filename, path in zip(("speciesPath", "nasa9Path", "transportPath"),
                                       ("species.pack", "nasa9.pack", "transport.pack"), paths):
            name = f"gas/{index}/{filename}"
            destination = output / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            with destination.open("xb") as stream, Path(path).expanduser().open("rb") as original:
                shutil.copyfileobj(original, stream)
            component[key] = name
            attachments.append(name)
        gas_components.append(component)
    if export_experiment:
        gas_sources = {"components": gas_components, "waterSpecies": None if args.gas_water == "none" else args.gas_water,
                       "phase": args.gas_phase} if gas_components else None
        packet = write_lab_experiment(source, output, attachments, profiles, gas_sources, args.decimate,
                                      args.coating_strain_guard, material_inputs)
        print(packet)
    if args.prepare_only:
        return 0
    wav = output / "render.wav"
    response = output / "loaded-response.tsv"
    command = [renderer, str(source)]
    if response_arguments:
        command += ["--loaded-body-response-only", str(response), *response_arguments]
    else:
        command.append(str(wav))
    if profile:
        command += ["--body-profile", str(output / "body.profile")]
    if stop_profile is not None:
        command += ["--stop-profile", str(output / "stops.profile")]
    if thermal_profile is not None:
        command += ["--contact-thermal", str(output / "contact-thermal.profile")]
    if args.decimate:
        command += ["--decimate"]
    if args.coating_strain_guard:
        command += ["--coating-strain-guard", *map(str, args.coating_strain_guard)]
    for component in gas_components:
        command += ["--gas-component", *(str(output / name) for name in component.values())]
    for key in ("water", "phase"):
        if value := getattr(args, "gas_" + key):
            command += ["--gas-" + key, value]
    log_name = "loaded-response.log" if response_arguments else "render.log"
    with (output / log_name).open("x") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    print(response if response_arguments else wav)
    if args.analysis:
        lines = text.split("\nschedule\n", 1)[0].splitlines()
        _, audio = row(lines, "audio")
        gas = " ".join(row(lines, "gas")[1][1:])
        cavity = [line for line in lines if line.split()[:1] == ["cavity_gas"]]
        atmosphere = "; gas " + gas + ("; " + cavity[0] if cavity else "; cavity_gas inherited")
        subprocess.run([sys.executable, str(HERE / "instrgrade" / "decompose.py"), "listen",
                        "--render", str(wav), "--family", args.family, "--note", args.note,
                        "--slot", args.note + "-" + args.variant, "--git-sha", args.git_sha,
                        "--observer", args.observer + atmosphere, "--recipe", args.recipe,
                        "--reason", "This material/atmosphere configuration has no co-characterized reference recording; this is an ungraded physics render.",
                        "--full-scale-pa", audio[3], "--out", str(args.analysis)], check=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error)) from None
