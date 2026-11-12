#!/usr/bin/env python3
"""Lower a wind description plus a specimen table set to frankensim-reed-performance-v1 files.

usage: python3 -I lower_reed_performance.py DESCRIPTION.tsv OUTDIR [FINGERING,...]

The description (frankensim.wind-description.v1) names a specimen directory
(frankensim.specimen-* tables), a fingering chart, a lumped reed, a physical
observer and a player-pressure table. Nothing is defaulted: a missing key, a
missing specimen row or a hole wider than its bore stops the run with the name
of the field. Every bore station and chimney becomes one record; the bore is
split at each chimney position and is never decimated.

The authored reed geometry and constitutive response stay fixed when the player
pressure changes. Changing embouchure requires changing its explicit physical
description; a pressure table must not silently resize the reed opening.

Optional exterior.* gas keys use exactly the ambient.* representations: gas_pack,
temperature_k/pressure_pa/composition, or temperature_k/pressure_pa/relative_humidity.
The native gas owner resolves both declarations and checks their compatibility.
Omitting exterior.* retains the native same-medium behavior; it never invents
room conditions. A supplied exterior reaches both opening loads and the physical
microphone, including descriptions with open tone holes.

A specimen may print other hole quantities: optional tonehole.height_quantity and
tonehole.diameter_quantity name them (chimney-height and min-diameter when absent)
and tonehole.omit lists holes left out unless a fingering opens them, when a hole
that prints no height takes tonehole.height_unprinted_m. A bore row at-hole-NN is a diameter at that
hole's position. mouthpiece.part - declares a reed with no mouthpiece, whose bore
stations start at the reed tip.

Optional card rows, each written only when its key is present (values are copied as written):
reed.lay_stiffness with reed.lay_exponent, reed.lay_onset_m and reed.lay_loss_s_per_m (the reed
row's four lay fields), reed.channel_length_m (slit_length), mouthpiece.barrel_pull_m (one cylinder
at the mouthpiece exit radius before the bore), bore.gas_field_mode with bore.gas_field_steps,
bore.gas_field_basis, bore.gas_field_species and the bore.gas_station.NN rows (the gas_field block),
tonehole.interaction_law with tonehole.interaction_epsilon (open_hole_interaction), tonehole.lengths,
reed.constriction_noise_seed with reed.constriction_noise_duct_radius_m, and tonehole.orifice_holes
`open` (open_hole_orifice N, N = the fingering's open-hole count).
A row `at.<fingering>.<key>` replaces one of PER_NOTE_KEYS for that fingering; the value `-` there
leaves the row out. A row `fingering.<name>` WRITTEN CONCERT REGISTER OPEN_HOLES adds a fingering the
chart does not print, placed before the first chart row of higher written pitch.
"""
import csv
import math
import os
import sys
from decimal import Decimal

CM_H2O_PA = 98.0665
PITCH_CLASS = {"C": 0, "D": 2, "E": 4, "F": 5, "G": 7, "A": 9, "B": 11}
FLAT_NAME = ["C", "Db", "D", "Eb", "E", "F", "Gb", "G", "Ab", "A", "Bb", "B"]
PER_NOTE_KEYS = ("reed.lay_stiffness", "reed.channel_length_m", "tonehole.lengths", "tonehole.orifice_holes")


def refuse(what):
    sys.exit(f"REFUSED: {what}")


def option(desc, key, fingering):
    """An optional row's value for one fingering: `at.<fingering>.<key>` replaces KEY, and `-` means no row."""
    value = desc.get(f"at.{fingering}.{key}", desc.get(key))
    if value is None or value == "-":
        return None
    if value in ("", "UNKNOWN"):
        refuse(f"description key {key} is empty or UNKNOWN and cannot be lowered")
    return value


def check_per_note_keys(desc, fingerings):
    for key in desc:
        if key.startswith("at."):
            name, _, rest = key[3:].partition(".")
            if name not in fingerings:
                refuse(f"{key} names fingering {name}, which neither the chart nor a fingering.* row declares")
            if rest not in PER_NOTE_KEYS:
                refuse(f"{key}: {rest} is not one of the per-note keys {', '.join(PER_NOTE_KEYS)}")


def extra_fingerings(desc, chart):
    """Chart rows plus each `fingering.<name>` row, inserted before the first chart row of higher written pitch."""
    rows = list(chart)
    for key in sorted(k for k in desc if k.startswith("fingering.")):
        fields = desc[key].split()
        if len(fields) != 4:
            refuse(f"{key} needs WRITTEN CONCERT REGISTER OPEN_HOLES")
        row = dict(zip(("written", "concert", "register", "open_holes"), fields), fingering=key[len("fingering."):])
        if any(r["fingering"] == row["fingering"] for r in rows):
            refuse(f"{key} repeats a fingering the chart already prints")
        at = next((i for i, r in enumerate(rows) if midi(r["written"]) > midi(row["written"])), len(rows))
        rows.insert(at, row)
    return rows


def table(path, schema):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or rows[0] != [schema]:
        refuse(f"{path}: first line must be {schema}")
    header = rows[1]
    return [dict(zip(header, row)) for row in rows[2:] if row]


def read_description(path):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or rows[0] != ["frankensim.wind-description.v1"]:
        refuse(f"{path}: first line must be frankensim.wind-description.v1")
    keys = {}
    for row in rows[2:]:
        if row:
            keys[row[0]] = row[1]
    return keys


def need(desc, key):
    if key not in desc or desc[key] == "":
        refuse(f"description is missing {key}")
    if desc[key] == "UNKNOWN":
        refuse(f"description key {key} is UNKNOWN and cannot be lowered")
    return desc[key]


def midi(name):
    letter, rest = name[0], name[1:]
    shift = 0
    while rest and rest[0] in "#b":
        shift += 1 if rest[0] == "#" else -1
        rest = rest[1:]
    return 12 * (int(rest) + 1) + PITCH_CLASS[letter] + shift


def uiowa_name(value):
    return f"{FLAT_NAME[value % 12]}{value // 12 - 1}"


def num(value):
    """Plain decimal with at most 10 places, no exponent."""
    text = format(round(value, 10), ".10f").rstrip("0").rstrip(".")
    return text if text not in ("", "-0") else "0"


def exp(value):
    """Shortest round-trip digits in the canonical LowerExp spelling (3.7e3, 2e-2, 0e0)."""
    if value == 0:
        return "0e0"
    sign, digits, exponent = Decimal(repr(float(value))).as_tuple()
    digits = list(digits)
    while len(digits) > 1 and digits[-1] == 0:
        digits.pop()
        exponent += 1
    lead = str(digits[0]) + ("." + "".join(map(str, digits[1:])) if len(digits) > 1 else "")
    return ("-" if sign else "") + lead + "e" + str(exponent + len(digits) - 1)


def dimensions(specimen_dir):
    found = {}
    for row in table(os.path.join(specimen_dir, "dimensions.tsv"), "frankensim.specimen-dimensions.v1"):
        found.setdefault((row["part"], row["feature"]), {})[row["quantity"]] = float(row["value_si"])
    return found


def bore_knots(dims):
    stations = {
        (q["axial-position"], q["diameter"] / 2.0)
        for (part, feature), q in dims.items()
        if part == "bore" and feature.startswith("station-")
    }
    # A diameter printed at a hole ("at-hole-NN") takes its position from that hole; two holes may share a station.
    stations |= {
        (dims[("tone-holes", feature[3:])]["axial-position"], q["diameter"] / 2.0)
        for (part, feature), q in dims.items()
        if part == "bore" and feature.startswith("at-hole-")
    }
    stations = sorted(stations)
    if len(stations) < 2:
        refuse("specimen has fewer than two bore stations")
    return stations


def holes_of(desc, dims, open_holes=()):
    """Each hole's position, height and radius; a specimen may name other quantities and omit unusable holes by key."""
    height = desc.get("tonehole.height_quantity", "chimney-height")
    diameter = desc.get("tonehole.diameter_quantity", "min-diameter")
    omitted = {int(h) for h in desc.get("tonehole.omit", "").split(",") if h}
    holes = {}
    for (part, feature), q in dims.items():
        if part == "tone-holes" and feature.startswith("hole-") and (int(feature[5:]) not in omitted or int(feature[5:]) in open_holes):
            if int(feature[5:]) in omitted and height not in q and "tonehole.height_unprinted_m" in desc:
                q = {**q, height: float(desc["tonehole.height_unprinted_m"])}
            for quantity in ("axial-position", height, diameter):
                if quantity not in q:
                    refuse(f"{feature} is missing {quantity}")
            holes[int(feature[5:])] = (q["axial-position"], q[height], q[diameter] / 2.0)
    return holes


def radius_at(knots, x):
    for (xa, ra), (xb, rb) in zip(knots, knots[1:]):
        if xa <= x <= xb:
            return ra + (rb - ra) * (x - xa) / (xb - xa)
    refuse(f"position {x} m lies outside the bore stations")


def mouthpiece_segments(desc, dims):
    part = need(desc, "mouthpiece.part")
    if part == "-":  # a double reed has no mouthpiece: the bore stations start at the reed tip
        if "mouthpiece.barrel_pull_m" in desc:
            refuse("mouthpiece.barrel_pull_m needs a mouthpiece")
        return []
    near_name, exit_name = need(desc, "mouthpiece.exit_cone_features").split(",")
    volume_feature = need(desc, "mouthpiece.volume_feature")
    near, exit_ = dims[(part, near_name)], dims[(part, exit_name)]
    volume = dims[(part, volume_feature)]
    cone_length = exit_["axial-position"] - near["axial-position"]
    r_in, r_out = near["diameter"] / 2.0, exit_["diameter"] / 2.0
    cone_volume = math.pi / 3.0 * cone_length * (r_in**2 + r_in * r_out + r_out**2)
    interior_length = volume["total-length"] - cone_length
    interior_volume = volume["internal-volume"] - cone_volume
    if interior_length <= 0 or interior_volume <= 0:
        refuse("mouthpiece total length or volume does not exceed its exit cone")
    l_eq = interior_volume / (math.pi * r_in**2)
    records = [
        f"cylinder {num(r_in)} {num(l_eq)}",
        f"cone {num(r_in)} {num(r_out)} {num(cone_length)}",
    ]
    if "mouthpiece.barrel_pull_m" in desc:  # the gap a pulled barrel opens, at the mouthpiece's exit radius
        records.append(f"cylinder {num(r_out)} {need(desc, 'mouthpiece.barrel_pull_m')}")
    return records


def axial_record(xa, ra, xb, rb):
    length = xb - xa
    if length <= 0:
        refuse(f"non-positive axial length between {xa} and {xb}")
    if abs(ra - rb) <= 1e-12:
        return f"cylinder {num(ra)} {num(length)}"
    return f"cone {num(ra)} {num(rb)} {num(length)}"


def duct_records(desc, dims, open_holes):
    knots = bore_knots(dims)
    holes = holes_of(desc, dims, open_holes)
    points = {x: (r, []) for x, r in knots}
    for number, (x, height, radius) in sorted(holes.items()):
        near = [k for k in points if abs(k - x) < 1e-9]
        key = near[0] if near else x
        if not near:
            points[key] = (radius_at(knots, x), [])
        points[key][1].append((number, height, radius))
    records = mouthpiece_segments(desc, dims)
    ordered = sorted(points.items())
    for index, (x, (bore_radius, here)) in enumerate(ordered):
        for number, height, radius in here:
            if radius >= bore_radius:
                refuse(f"hole-{number:02d} radius {radius} is not below its bore radius {bore_radius}")
            sigma = 1 if number in open_holes else 0
            records.append(f"hole {num(radius)} {num(height)} {num(bore_radius)} {sigma}")
        if index + 1 < len(ordered):
            nx, (nr, _) = ordered[index + 1]
            records.append(axial_record(x, bore_radius, nx, nr))
    return records


def pressure_curve(desc, root):
    table_rel = need(desc, "player.pressure_table")
    if table_rel == "-":
        return None
    group = need(desc, "player.instrument_group")
    level = need(desc, "player.level")
    path = os.path.join(root, table_rel)
    sums = {}
    for row in table(path, "frankensim.player-blowing-pressure.v1"):
        if row["instrument"] == group and row["level"] == level:
            sums.setdefault(midi(row["pitch"]), []).append(float(row["mean_cm_h2o"]))
    if not sums:
        refuse(f"no {level} rows for {group} in {path}")
    return sorted((m, sum(v) / len(v) * CM_H2O_PA) for m, v in sums.items())


def pressure_for(desc, curve, note, reach, step):
    if curve is None:
        value = float(need(desc, "reed.operating_blowing_pa"))
        return round(value / step) * step, "fixed-operating-point"
    lo, hi = curve[0][0], curve[-1][0]
    if note < lo - reach or note > hi + reach:
        refuse(f"{uiowa_name(note)} lies more than {reach} semitones outside the published pressure range")
    if note <= lo:
        value, rule = curve[0][1], "held-nearest" if note < lo else "tabulated"
    elif note >= hi:
        value, rule = curve[-1][1], "held-nearest" if note > hi else "tabulated"
    else:
        value, rule = None, "interpolated"
        for (ma, pa), (mb, pb) in zip(curve, curve[1:]):
            if ma <= note <= mb:
                value = pa + (pb - pa) * (note - ma) / (mb - ma)
                rule = "tabulated" if note in (ma, mb) else "interpolated"
    return round(value / step) * step, rule


def gesture(desc, pressure):
    attack = float(need(desc, "gesture.attack_s"))
    start = float(need(desc, "gesture.attack_start_s"))
    if "tongue.contact_duration_s" in desc:
        start = round(start - float(need(desc, "tongue.contact_duration_s")), 10)
        if start < 0:
            refuse("tongue.contact_duration_s exceeds gesture.attack_start_s")
    release_at = float(need(desc, "gesture.release_start_s"))
    release = float(need(desc, "gesture.release_s"))
    return (
        "frankensim-gesture-schedule-v1\n"
        f"control_rate_hz\t{need(desc, 'gesture.control_rate_hz')}\n"
        "tracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t2\n"
        f"event\t{exp(start)}\t{exp(attack)}\tpressure-pa\t{exp(pressure)}\n"
        f"event\t{exp(release_at)}\t{exp(release)}\tpressure-pa\t0e0\n"
    )


def tongue_record(desc):
    if "tongue.contact_area_m2" not in desc:
        return ""
    area_t = float(need(desc, "tongue.contact_area_m2"))
    area_r = float(need(desc, "reed.effective_area_m2"))
    force_n = float(need(desc, "tongue.contact_force_n"))
    release_start = float(need(desc, "gesture.attack_start_s"))
    contact_duration = float(need(desc, "tongue.contact_duration_s"))
    unload_start = round(release_start - contact_duration, 10)
    unload_s = float(need(desc, "gesture.attack_s"))
    release_s = float(need(desc, "tongue.release_s"))
    if unload_start < 0 or unload_start + unload_s > release_start:
        refuse("tongue contact window does not enclose the pre-release mouth-pressure ramp")
    return f"tongue_contact {num(area_t)} {num(area_r)} {num(force_n)} {num(release_start)} {num(release_s)}\n"


def reed_record(desc, fingering=None):
    """The reed row: the record below plus, with reed.lay_stiffness, the four lay fields as written."""
    record = reed_values(desc)
    if option(desc, "reed.lay_stiffness", fingering) is None:
        return record
    lay = [option(desc, "reed.lay_stiffness", fingering)] + [
        need(desc, key) for key in ("reed.lay_exponent", "reed.lay_onset_m", "reed.lay_loss_s_per_m")
    ]
    return f"{record} {' '.join(lay)}"


def reed_values(desc):
    """The `reed` header values from the authored geometry and constitutive response.

    Per-area mode (reed.stiffness_per_area_pa_m present): at the published operating point
    p_m = reed.operating_blowing_pa (Chatziioannou and van Walstijn 2009 Table 1), the lip-reduced
    equilibrium rest opening is held at y_m = h0 = 0.40 mm across all notes with closing pressure
    P_M = K_a y_m = 3464 Pa, stiffness k = K_a S_r = 659.026 N/m, mass m = m_a S_r = 3.805e-6 kg,
    and damping ratio zeta = g / (2 omega_r) with omega_r = sqrt(K_a / m_a) = sqrt(k / m).
    """
    height = need(desc, "reed.rest_opening_m")
    width = need(desc, "reed.width_m")
    if "reed.stiffness_per_area_pa_m" not in desc:
        keys = ("reed.closing_pressure_pa", "reed.mass_kg", "reed.stiffness_n_m", "reed.damping_ratio")
        return " ".join([height, width] + [need(desc, key) for key in keys])
    h0 = float(height)
    area = float(need(desc, "reed.effective_area_m2"))
    k_a = float(need(desc, "reed.stiffness_per_area_pa_m"))
    m_a = float(need(desc, "reed.mass_per_area_kg_m2"))
    g = float(need(desc, "reed.damping_rate_1_s"))
    closing = k_a * h0
    if "reed.tip_opening_cap_m" in desc and h0 > float(desc["reed.tip_opening_cap_m"]):
        refuse(f"equilibrium reed opening {h0} m exceeds tip opening cap {desc['reed.tip_opening_cap_m']} m")
    stiffness = k_a * area
    mass = m_a * area
    zeta = g / (2.0 * math.sqrt(k_a / m_a))
    return f"{num(h0)} {width} {num(closing)} {mass!r} {stiffness!r} {zeta!r}"


def gas_record(desc, prefix):
    """Serialize one physical gas declaration without deriving or replacing its properties.

    Both the bore and the exterior go through the same native gas resolver.
    A pack remains the authority when a description also retains its source
    coordinates/composition as explanatory columns, as the clarinet does.
    """
    if f"{prefix}.gas_pack" in desc:
        pack_id = need(desc, f"{prefix}.gas_pack")
        manifest = os.path.join("data", "matdb", "seed-v1", pack_id, "manifest.tsv")
        if not os.path.isfile(manifest):
            refuse(f"missing gas pack manifest {manifest}")
        return f"pack {pack_id}"
    if f"{prefix}.composition" in desc:
        return (
            f"mixture {need(desc, f'{prefix}.temperature_k')} {need(desc, f'{prefix}.pressure_pa')} "
            f"{need(desc, f'{prefix}.composition')}"
        )
    return " ".join(
        need(desc, f"{prefix}.{key}")
        for key in ("temperature_k", "pressure_pa", "relative_humidity")
    )


def gas_field_rows(desc):
    """The gas_field block: its header from the bore.gas_field_* keys, then the bore.gas_station.NN rows as written."""
    if "bore.gas_field_mode" not in desc:
        return ""
    species = need(desc, "bore.gas_field_species").split()
    stations = [need(desc, key) for key in sorted(k for k in desc if k.startswith("bore.gas_station."))]
    if len(stations) < 2:
        refuse("bore.gas_field_mode needs at least two bore.gas_station.NN rows")
    for station in stations:
        if len(station.split()) != 2 + len(species):
            refuse(f"gas station '{station}' needs X_M T_K and one fraction per species ({len(species)})")
    head = (
        f"gas_field {need(desc, 'bore.gas_field_mode')} {need(desc, 'bore.gas_field_steps')} {len(stations)} "
        f"{need(desc, 'bore.gas_field_basis')} {len(species)} {' '.join(species)}"
    )
    return "\n".join([head] + [f"gas_station {s}" for s in stations]) + "\n"


def tail_rows(desc, records, fingering):
    """The rows between the bore and compile_limits, in the order the reed parser reads them."""
    rows = gas_field_rows(desc)
    if "tonehole.interaction_law" in desc:
        rows += f"open_hole_interaction {need(desc, 'tonehole.interaction_law')} {need(desc, 'tonehole.interaction_epsilon')}\n"
    lengths = option(desc, "tonehole.lengths", fingering)
    if lengths is not None:
        rows += f"tonehole_lengths {lengths}\n"
    if "reed.constriction_noise_seed" in desc:
        rows += (
            f"constriction_noise {need(desc, 'reed.constriction_noise_seed')} "
            f"{need(desc, 'reed.constriction_noise_duct_radius_m')}\n"
        )
    orifice = option(desc, "tonehole.orifice_holes", fingering)
    if orifice is not None:
        if orifice != "open":
            refuse(f"tonehole.orifice_holes is '{orifice}'; the lowerer takes 'open' (every open hole) or '-'")
        count = sum(1 for r in records if r.startswith("hole ") and r.endswith(" 1"))
        if count == 0:
            refuse(f"tonehole.orifice_holes open on fingering {fingering}, which opens no hole")
        rows += f"open_hole_orifice {count}\n"
    return rows


def performance(desc, records, pressure, fingering=None):
    rate = int(need(desc, "audio.rate_hz"))
    samples = round(float(need(desc, "audio.duration_s")) * rate)
    reed = reed_record(desc, fingering)
    channel = option(desc, "reed.channel_length_m", fingering)
    if channel is not None:
        reed += f"\nslit_length {channel}"
    ambient = gas_record(desc, "ambient")
    exterior = (
        f"exterior {gas_record(desc, 'exterior')}\n"
        if any(key.startswith("exterior.") for key in desc) else ""
    )
    observer_law = f"observer_law {need(desc, 'observer.law')}\n" if "observer.law" in desc else ""
    head = (
        "frankensim-reed-performance-v1\n"
        f"audio {rate} {samples} {need(desc, 'audio.full_scale_pa')}\n"
        f"ambient {ambient}\nreed {reed}\n{tongue_record(desc)}{need(desc, 'observer.record')}\n"
        f"{exterior}{observer_law}termination {need(desc, 'observer.termination')}\nsegments {len(records)}\n"
    )
    return (
        head + "\n".join(records) + "\n" + tail_rows(desc, records, fingering)
        + f"compile_limits {need(desc, 'compile_limits')}\nschedule\n" + gesture(desc, pressure)
    )


def main(argv):
    if len(argv) not in (3, 4):
        sys.exit(__doc__)
    desc_path, out_dir = argv[1], argv[2]
    wanted = set(argv[3].split(",")) if len(argv) == 4 else None
    root = os.getcwd()
    desc = read_description(desc_path)
    specimen = os.path.join(root, need(desc, "specimen_dir"))
    dims = dimensions(specimen)
    chart = table(os.path.join(root, need(desc, "fingering_chart")), "frankensim.specimen-fingerings.v1")
    chart = extra_fingerings(desc, chart)
    check_per_note_keys(desc, {row["fingering"] for row in chart})
    curve = pressure_curve(desc, root)
    reach = int(need(desc, "player.extrapolation_semitones_max"))
    step = float(need(desc, "player.pressure_round_pa"))
    level = need(desc, "player.level")
    prefix = need(desc, "output.prefix")
    os.makedirs(out_dir, exist_ok=True)
    index = ["fingering\twritten\tconcert\tregister\topen_holes\tpressure_pa\tpressure_rule\tsegments\ttone_holes\tfile"]
    for row in chart:
        if wanted is not None and row["fingering"] not in wanted:
            continue
        open_holes = set() if row["open_holes"] == "-" else {int(h) for h in row["open_holes"].split(",")}
        note = midi(row["concert"])
        pressure, rule = pressure_for(desc, curve, note, reach, step)
        records = duct_records(desc, dims, open_holes)
        name = f"{prefix}.{level}.{uiowa_name(note)}.performance"
        with open(os.path.join(out_dir, name), "w", encoding="utf-8", newline="\n") as handle:
            handle.write(performance(desc, records, pressure, row["fingering"]))
        holes = sum(1 for r in records if r.startswith("hole "))
        index.append(
            f"{row['fingering']}\t{row['written']}\t{uiowa_name(note)}\t{row['register']}\t"
            f"{row['open_holes']}\t{num(pressure)}\t{rule}\t{len(records)}\t{holes}\t{name}"
        )
    with open(os.path.join(out_dir, "index.tsv"), "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(index) + "\n")
    print(f"wrote {len(index) - 1} performances to {out_dir}")


if __name__ == "__main__":
    main(sys.argv)
