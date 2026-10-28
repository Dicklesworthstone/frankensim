#!/usr/bin/env python3
"""Render an existing material-string assembly and optionally add it to Instrument Lab.

All mechanics, material admission, gas resolution and PCM encoding stay in the
native material_strings executable. This launcher stages its original inputs;
it neither estimates missing properties nor tunes a string after substitution.
Contact profiles reach that same native solve. Optional native decimation changes
only the observed audio clock, retaining the authored mechanical resolution.
The output directory must be new. Python's standard library suffices for native
rendering; --analysis additionally uses the existing instrgrade dependencies.
"""
from __future__ import annotations

import argparse
import base64
import json
import math
import shutil
import subprocess
import sys
from decimal import Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent


def row(lines: list[str], key: str) -> tuple[int, list[str]]:
    matches = [(i, line.split()) for i, line in enumerate(lines) if line.startswith(key + " ")]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one {key} record")
    return matches[0]


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


def stage_packs(texts: list[str], root: Path, replacements: list[str], output: Path) -> list[str]:
    # These are the existing binary pack formats, copied byte for byte.
    # Rust remains the material decoder and property/domain authority.
    available = {str(p.relative_to(root)): p for suffix in ("*.fsmatpk", "*.fsintpk") for p in root.rglob(suffix)}
    tokens = {token for text in texts for line in text.splitlines() if not line.startswith("#") for token in line.split()}
    replaced = set()
    for binding in replacements:
        name, sep, source = binding.partition("=")
        if not sep or not name or name not in tokens or name in replaced or not name.endswith((".fsmatpk", ".fsintpk")):
            raise ValueError(f"--pack must name one distinct referenced path: {binding!r}")
        path = Path(name)
        if path.is_absolute() or any(p in (".", "..") for p in name.split("/")):
            raise ValueError("material pack paths must stay relative to the performance directory")
        available[name] = Path(source).expanduser().resolve()
        replaced.add(name)
    required = {t for t in tokens if t.endswith((".fsmatpk", ".fsintpk"))} | replaced
    for name in sorted(required):
        if name not in available:
            raise ValueError(f"missing material pack {name}; supply --pack-root or --pack {name}=PATH")
        destination = output / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("xb") as stream, available[name].open("rb") as source:
            shutil.copyfileobj(source, stream)
    return sorted(required)


def write_lab_experiment(source: Path, output: Path, attachments: list[str], profiles: dict,
                         gas_sources: dict | None, decimate: bool) -> Path:
    """Export the existing native HTTP packet; material/physics decoding stays in Rust.

    The lab's shared packet reader checks its transport bounds when opening this
    file. Prepared source bytes and all referenced source paths remain intact.
    """
    encoded = lambda data: base64.b64encode(data).decode("ascii")
    packet = {"format": "frankensim-native-wind-input-v1", "renderer": "strings",
              "source": encoded(source.read_bytes()), "decimate": decimate,
              "attachments": [{"path": name, "bytes": encoded((output / name).read_bytes())}
                              for name in attachments]}
    if profiles:
        packet["stringProfiles"] = profiles
    if gas_sources is not None:
        packet["gasSources"] = gas_sources
    path = output / "experiment.native.json"
    with path.open("x", encoding="utf-8") as stream:
        json.dump(packet, stream, separators=(",", ":"), ensure_ascii=False)
    return path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--renderer", help="path to material_strings; required for audio rendering or shell assembly")
    parser.add_argument("--performance", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True, help="new directory for the staged sources, WAV and native log")
    parser.add_argument("--pack-root", type=Path)
    parser.add_argument("--pack", action="append", default=[], metavar="REFERENCED_PATH=SOURCE_PACK")
    parser.add_argument("--body-profile", type=Path)
    parser.add_argument("--stop-profile", type=Path, help="original finger/fret/bow-surface contact profile")
    parser.add_argument("--contact-thermal", type=Path, help="original finite contact-heat profile; thermal properties come from the same assembly packs")
    parser.add_argument("--shell-table", type=Path)
    parser.add_argument("--assembly-options", type=Path)
    parser.add_argument("--gas-record", help="exact shared gas record after 'gas', e.g. 'mixture 293.15 101325 mole 2 N2 0.8 O2 0.2'")
    parser.add_argument("--gas-component", action="append", nargs=3, default=[], metavar=("SPECIES", "NASA9", "TRANSPORT"))
    parser.add_argument("--gas-water")
    parser.add_argument("--gas-phase", choices=("frozen", "liquid-water"))
    parser.add_argument("--full-scale-pa", type=float, help="explicit PCM scale; use the same value for every comparison")
    parser.add_argument("--bow-card", type=Path, help="reuse only this legacy card's canonical bow commands and duration")
    parser.add_argument("--string-index", type=int)
    parser.add_argument("--bow-distance-from-bridge", type=float)
    parser.add_argument("--decimate", action="store_true", help="use the native 48 kHz output decimator without changing the mechanical clock")
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
    renderer = shutil.which(args.renderer) if args.renderer else None
    if (args.renderer and renderer is None) or (renderer is None and (not args.prepare_only or args.shell_table)):
        parser.error("--renderer must name an available material_strings executable")
    renderer = str(Path(renderer).resolve()) if renderer else None
    export_experiment = args.lab_experiment or args.prepare_only
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
    header, schedule = text.split("\nschedule\n", 1)
    lines = header.splitlines()
    if args.gas_record is not None:
        if not args.gas_record.strip() or "\n" in args.gas_record or "\r" in args.gas_record:
            parser.error("--gas-record must be one nonempty record")
        lines[row(lines, "gas")[0]] = "gas " + args.gas_record
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
    attachments += stage_packs([text, profile, (stop_profile or b"").decode("utf-8")],
                              args.pack_root or args.performance.parent, args.pack, output)
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
        packet = write_lab_experiment(source, output, attachments, profiles, gas_sources, args.decimate)
        print(packet)
    if args.prepare_only:
        return 0
    wav = output / "render.wav"
    command = [renderer, str(source), str(wav)]
    if profile:
        command += ["--body-profile", str(output / "body.profile")]
    if stop_profile is not None:
        command += ["--stop-profile", str(output / "stops.profile")]
    if thermal_profile is not None:
        command += ["--contact-thermal", str(output / "contact-thermal.profile")]
    if args.decimate:
        command += ["--decimate"]
    for component in gas_components:
        command += ["--gas-component", *(str(output / name) for name in component.values())]
    for key in ("water", "phase"):
        if value := getattr(args, "gas_" + key):
            command += ["--gas-" + key, value]
    with (output / "render.log").open("x") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    print(wav)
    if args.analysis:
        _, audio = row(text.split("\nschedule\n", 1)[0].splitlines(), "audio")
        gas = " ".join(row(text.split("\nschedule\n", 1)[0].splitlines(), "gas")[1][1:])
        subprocess.run([sys.executable, str(HERE / "instrgrade" / "decompose.py"), "listen",
                        "--render", str(wav), "--family", args.family, "--note", args.note,
                        "--slot", args.note + "-" + args.variant, "--git-sha", args.git_sha,
                        "--observer", args.observer + "; gas " + gas, "--recipe", args.recipe,
                        "--reason", "This material/atmosphere configuration has no co-characterized reference recording; this is an ungraded physics render.",
                        "--full-scale-pa", audio[3], "--out", str(args.analysis)], check=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error)) from None
