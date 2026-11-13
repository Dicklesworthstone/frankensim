//! `music_render instrument`: a wind description lowered in Rust to the reed card `wind` renders
//! (bead rc-1109 .12.1). The description is `frankensim.wind-description.v1` with its specimen
//! tables and fingering chart, as in `data/instruments/descriptions/<family>/description.tsv`.
//!
//! ```text
//! music_render instrument DESCRIPTION.tsv FINGERING OUT.wav [wind options]
//! music_render instrument DESCRIPTION.tsv --lower OUTDIR
//! ```
//!
//! The first form lowers one fingering and writes its card beside OUT.wav as OUT.performance.
//! That file is refused if it already exists. It then renders the card through `music_render wind`
//! with the same options, so the description and the card give the same WAV. The second form writes every
//! fingering's card and `index.tsv`, the files `data/instruments/descriptions/lower_reed_performance.py`
//! writes, byte for byte; that script stays as the cross-check until this parity is recorded.
//!
//! Paths inside the description are relative to the working directory, the repository root. Nothing
//! is defaulted: a missing key, specimen row or chart column refuses with its name.
//! The per-note rows (`at.<fingering>.<key>`, `-` for no row) and the extra fingerings
//! (`fingering.<name>`) are documented in the Python lowerer's header.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const DESCRIPTION_SCHEMA: &str = "frankensim.wind-description.v1";
const CM_H2O_PA: f64 = 98.0665;
const FLAT_NAME: [&str; 12] = [
    "C", "Db", "D", "Eb", "E", "F", "Gb", "G", "Ab", "A", "Bb", "B",
];
const PER_NOTE_KEYS: [&str; 4] = [
    "reed.lay_stiffness",
    "reed.channel_length_m",
    "tonehole.lengths",
    "tonehole.orifice_holes",
];
const INDEX_HEADER: &str = "fingering\twritten\tconcert\tregister\topen_holes\tpressure_pa\tpressure_rule\tsegments\ttone_holes\tfile";

type Lowered<T> = Result<T, String>;

/// The rows of a tab-separated file whose first line is exactly `schema`. A blank line is no row.
fn tsv_rows(text: &str, origin: &str, schema: &str) -> Lowered<Vec<Vec<String>>> {
    let rows: Vec<Vec<String>> = text
        .lines()
        .map(|line| {
            if line.is_empty() {
                Vec::new()
            } else {
                line.split('\t').map(str::to_owned).collect()
            }
        })
        .collect();
    if rows.first().map(|row| row.len() == 1 && row[0] == schema) != Some(true) {
        return Err(format!("{origin}: first line must be {schema}"));
    }
    if let Some(field) = rows.iter().flatten().find(|field| field.starts_with('"')) {
        return Err(format!(
            "{origin}: quoted field {field:?}; the lowerer reads unquoted tab-separated text"
        ));
    }
    Ok(rows)
}

fn read_text(path: &Path) -> Lowered<String> {
    std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))
}

/// A specimen or chart table: the header row names each later row's fields.
fn table(path: &Path, schema: &str) -> Lowered<Vec<BTreeMap<String, String>>> {
    let origin = path.display().to_string();
    let rows = tsv_rows(&read_text(path)?, &origin, schema)?;
    let header = rows
        .get(1)
        .ok_or_else(|| format!("{origin}: no header row"))?;
    Ok(rows[2..]
        .iter()
        .filter(|row| !row.is_empty())
        .map(|row| header.iter().cloned().zip(row.iter().cloned()).collect())
        .collect())
}

fn field<'a>(row: &'a BTreeMap<String, String>, name: &str, origin: &str) -> Lowered<&'a str> {
    row.get(name)
        .map(String::as_str)
        .ok_or_else(|| format!("{origin}: a row has no {name} field"))
}

fn parse_f64(text: &str, what: &str) -> Lowered<f64> {
    text.parse::<f64>()
        .map_err(|_| format!("{what} is not a number: {text:?}"))
}

/// The description's keys; a later row with the same key replaces an earlier one.
pub(super) struct Description {
    keys: BTreeMap<String, String>,
}

impl Description {
    pub(super) fn parse(text: &str, origin: &str) -> Lowered<Self> {
        let rows = tsv_rows(text, origin, DESCRIPTION_SCHEMA)?;
        let mut keys = BTreeMap::new();
        for row in rows.iter().skip(2).filter(|row| !row.is_empty()) {
            let value = row
                .get(1)
                .ok_or_else(|| format!("{origin}: key {} has no value column", row[0]))?;
            keys.insert(row[0].clone(), value.clone());
        }
        Ok(Self { keys })
    }

    fn get(&self, key: &str) -> Option<&str> {
        self.keys.get(key).map(String::as_str)
    }

    fn has(&self, key: &str) -> bool {
        self.keys.contains_key(key)
    }

    fn has_prefix(&self, prefix: &str) -> bool {
        self.keys.keys().any(|key| key.starts_with(prefix))
    }

    fn need(&self, key: &str) -> Lowered<&str> {
        match self.get(key) {
            None | Some("") => Err(format!("description is missing {key}")),
            Some("UNKNOWN") => Err(format!(
                "description key {key} is UNKNOWN and cannot be lowered"
            )),
            Some(value) => Ok(value),
        }
    }

    fn number(&self, key: &str) -> Lowered<f64> {
        parse_f64(self.need(key)?, key)
    }

    /// An optional row's value for one fingering: `at.<fingering>.<key>` replaces KEY, `-` means no row.
    fn option(&self, key: &str, fingering: &str) -> Lowered<Option<&str>> {
        let value = self
            .get(&format!("at.{fingering}.{key}"))
            .or_else(|| self.get(key));
        match value {
            None | Some("-") => Ok(None),
            Some("") | Some("UNKNOWN") => Err(format!(
                "description key {key} is empty or UNKNOWN and cannot be lowered"
            )),
            Some(value) => Ok(Some(value)),
        }
    }
}

fn midi(name: &str) -> Lowered<i32> {
    let mut chars = name.chars();
    let pitch = match chars.next() {
        Some('C') => 0,
        Some('D') => 2,
        Some('E') => 4,
        Some('F') => 5,
        Some('G') => 7,
        Some('A') => 9,
        Some('B') => 11,
        _ => {
            return Err(format!(
                "note name {name:?} does not start with a letter A to G"
            ));
        }
    };
    let mut rest = chars.as_str();
    let mut shift = 0;
    while let Some(sign) = rest.chars().next().filter(|c| *c == '#' || *c == 'b') {
        shift += if sign == '#' { 1 } else { -1 };
        rest = &rest[1..];
    }
    let octave: i32 = rest
        .parse()
        .map_err(|_| format!("note name {name:?} has no octave number"))?;
    Ok(12 * (octave + 1) + pitch + shift)
}

fn uiowa_name(value: i32) -> String {
    format!(
        "{}{}",
        FLAT_NAME[value.rem_euclid(12) as usize],
        value.div_euclid(12) - 1
    )
}

/// Plain decimal with at most 10 places, no exponent; "0" for zero and for a rounded "-0".
fn num(value: f64) -> String {
    let text = format!("{value:.10}");
    let text = text.trim_end_matches('0').trim_end_matches('.');
    if text.is_empty() || text == "-0" {
        "0".into()
    } else {
        text.into()
    }
}

/// Shortest round-trip digits in the LowerExp spelling (3.7e3, 2e-2, 0e0).
fn exp(value: f64) -> String {
    if value == 0.0 {
        "0e0".into()
    } else {
        format!("{value:e}")
    }
}

/// Python's `repr` of a float: shortest round-trip digits, fixed for decimal exponents -4 to 15,
/// otherwise `d.ddde-06` with a sign and at least two exponent digits.
fn py_repr(value: f64) -> Lowered<String> {
    if !value.is_finite() {
        return Err(format!("reed value {value} is not finite"));
    }
    if value == 0.0 {
        return Ok(if value.is_sign_negative() {
            "-0.0"
        } else {
            "0.0"
        }
        .into());
    }
    let sci = format!("{:e}", value.abs());
    let (mantissa, exponent) = sci.split_once('e').ok_or("LowerExp has no exponent")?;
    let exponent: i32 = exponent
        .parse()
        .map_err(|_| "LowerExp exponent is not an integer")?;
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    let body = if (-4..16).contains(&exponent) {
        if exponent >= 0 {
            let whole = exponent as usize + 1;
            if digits.len() <= whole {
                format!("{digits}{}.0", "0".repeat(whole - digits.len()))
            } else {
                format!("{}.{}", &digits[..whole], &digits[whole..])
            }
        } else {
            format!("0.{}{digits}", "0".repeat((-exponent - 1) as usize))
        }
    } else {
        format!(
            "{mantissa}e{}{:02}",
            if exponent < 0 { '-' } else { '+' },
            exponent.abs()
        )
    };
    Ok(format!("{}{body}", if value < 0.0 { "-" } else { "" }))
}

/// Python's `round(value, 10)`: the double nearest the value rounded to ten decimals.
fn round10(value: f64) -> f64 {
    format!("{value:.10}").parse().unwrap_or(value)
}

type Dimensions = BTreeMap<(String, String), BTreeMap<String, f64>>;

fn dimensions(specimen: &Path) -> Lowered<Dimensions> {
    let path = specimen.join("dimensions.tsv");
    let origin = path.display().to_string();
    let mut found = Dimensions::new();
    for row in table(&path, "frankensim.specimen-dimensions.v1")? {
        let key = (
            field(&row, "part", &origin)?.to_owned(),
            field(&row, "feature", &origin)?.to_owned(),
        );
        let value = parse_f64(
            field(&row, "value_si", &origin)?,
            &format!("{origin} {} {}", key.0, key.1),
        )?;
        found
            .entry(key)
            .or_default()
            .insert(field(&row, "quantity", &origin)?.to_owned(), value);
    }
    Ok(found)
}

fn quantity(q: &BTreeMap<String, f64>, name: &str, what: &str) -> Lowered<f64> {
    q.get(name)
        .copied()
        .ok_or_else(|| format!("{what} is missing {name}"))
}

fn bore_knots(dims: &Dimensions) -> Lowered<Vec<(f64, f64)>> {
    let mut stations = Vec::new();
    for ((part, feature), q) in dims {
        if part == "bore" && feature.starts_with("station-") {
            stations.push((
                quantity(q, "axial-position", feature)?,
                quantity(q, "diameter", feature)? / 2.0,
            ));
        }
        // A diameter printed at a hole ("at-hole-NN") takes its position from that hole.
        if part == "bore" && feature.starts_with("at-hole-") {
            let hole = &feature[3..];
            let position = dims
                .get(&("tone-holes".to_owned(), hole.to_owned()))
                .ok_or_else(|| {
                    format!("{feature} names {hole}, which the specimen does not print")
                })?;
            stations.push((
                quantity(position, "axial-position", hole)?,
                quantity(q, "diameter", feature)? / 2.0,
            ));
        }
    }
    stations.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    stations.dedup();
    if stations.len() < 2 {
        return Err("specimen has fewer than two bore stations".into());
    }
    Ok(stations)
}

fn hole_number(feature: &str) -> Lowered<i64> {
    feature[5..]
        .parse()
        .map_err(|_| format!("{feature} is not hole-NN"))
}

/// Each hole's position, height and radius; omitted holes stay out unless the fingering opens them.
fn holes_of(
    desc: &Description,
    dims: &Dimensions,
    open: &BTreeSet<i64>,
) -> Lowered<BTreeMap<i64, (f64, f64, f64)>> {
    let height = desc
        .get("tonehole.height_quantity")
        .unwrap_or("chimney-height");
    let diameter = desc
        .get("tonehole.diameter_quantity")
        .unwrap_or("min-diameter");
    let mut omitted = BTreeSet::new();
    for hole in desc
        .get("tonehole.omit")
        .unwrap_or("")
        .split(',')
        .filter(|h| !h.is_empty())
    {
        omitted.insert(
            hole.parse::<i64>()
                .map_err(|_| format!("tonehole.omit lists {hole:?}, not a hole number"))?,
        );
    }
    let mut holes = BTreeMap::new();
    for ((part, feature), q) in dims {
        if part != "tone-holes" || !feature.starts_with("hole-") {
            continue;
        }
        let number = hole_number(feature)?;
        if omitted.contains(&number) && !open.contains(&number) {
            continue;
        }
        let mut q = q.clone();
        if omitted.contains(&number) && !q.contains_key(height) {
            if let Some(unprinted) = desc.get("tonehole.height_unprinted_m") {
                q.insert(
                    height.to_owned(),
                    parse_f64(unprinted, "tonehole.height_unprinted_m")?,
                );
            }
        }
        let position = quantity(&q, "axial-position", feature)?;
        let hole_height = quantity(&q, height, feature)?;
        let hole_diameter = quantity(&q, diameter, feature)?;
        holes.insert(number, (position, hole_height, hole_diameter / 2.0));
    }
    Ok(holes)
}

fn radius_at(knots: &[(f64, f64)], x: f64) -> Lowered<f64> {
    for pair in knots.windows(2) {
        let ((xa, ra), (xb, rb)) = (pair[0], pair[1]);
        if xa <= x && x <= xb {
            return Ok(ra + (rb - ra) * (x - xa) / (xb - xa));
        }
    }
    Err(format!("position {x} m lies outside the bore stations"))
}

fn mouthpiece_segments(desc: &Description, dims: &Dimensions) -> Lowered<Vec<String>> {
    let part = desc.need("mouthpiece.part")?;
    if part == "-" {
        // A double reed has no mouthpiece: the bore stations start at the reed tip.
        if desc.has("mouthpiece.barrel_pull_m") {
            return Err("mouthpiece.barrel_pull_m needs a mouthpiece".into());
        }
        return Ok(Vec::new());
    }
    let features = desc.need("mouthpiece.exit_cone_features")?;
    let (near_name, exit_name) = features
        .split_once(',')
        .filter(|(_, exit)| !exit.contains(','))
        .ok_or("mouthpiece.exit_cone_features needs two features separated by one comma")?;
    let volume_feature = desc.need("mouthpiece.volume_feature")?;
    let feature = |name: &str| {
        dims.get(&(part.to_owned(), name.to_owned()))
            .ok_or_else(|| format!("specimen has no {part} {name} row"))
    };
    let (near, exit, volume) = (
        feature(near_name)?,
        feature(exit_name)?,
        feature(volume_feature)?,
    );
    let cone_length =
        quantity(exit, "axial-position", exit_name)? - quantity(near, "axial-position", near_name)?;
    let r_in = quantity(near, "diameter", near_name)? / 2.0;
    let r_out = quantity(exit, "diameter", exit_name)? / 2.0;
    let cone_volume =
        std::f64::consts::PI / 3.0 * cone_length * (r_in * r_in + r_in * r_out + r_out * r_out);
    let interior_length = quantity(volume, "total-length", volume_feature)? - cone_length;
    let interior_volume = quantity(volume, "internal-volume", volume_feature)? - cone_volume;
    if interior_length <= 0.0 || interior_volume <= 0.0 {
        return Err("mouthpiece total length or volume does not exceed its exit cone".into());
    }
    let l_eq = interior_volume / (std::f64::consts::PI * (r_in * r_in));
    let mut records = vec![
        format!("cylinder {} {}", num(r_in), num(l_eq)),
        format!("cone {} {} {}", num(r_in), num(r_out), num(cone_length)),
    ];
    if desc.has("mouthpiece.barrel_pull_m") {
        // The gap a pulled barrel opens, at the mouthpiece's exit radius.
        records.push(format!(
            "cylinder {} {}",
            num(r_out),
            desc.need("mouthpiece.barrel_pull_m")?
        ));
    }
    Ok(records)
}

fn axial_record(xa: f64, ra: f64, xb: f64, rb: f64) -> Lowered<String> {
    let length = xb - xa;
    if length <= 0.0 {
        return Err(format!("non-positive axial length between {xa} and {xb}"));
    }
    Ok(if (ra - rb).abs() <= 1e-12 {
        format!("cylinder {} {}", num(ra), num(length))
    } else {
        format!("cone {} {} {}", num(ra), num(rb), num(length))
    })
}

/// Every bore station and chimney as one record; the bore is split at each chimney, never decimated.
fn duct_records(
    desc: &Description,
    dims: &Dimensions,
    open: &BTreeSet<i64>,
) -> Lowered<Vec<String>> {
    let knots = bore_knots(dims)?;
    let holes = holes_of(desc, dims, open)?;
    // Points in insertion order, as the Python lowerer's dict: a hole joins the first point within 1 nm.
    let mut points: Vec<(f64, f64, Vec<(i64, f64, f64)>)> = Vec::new();
    for &(x, r) in &knots {
        match points.iter_mut().find(|point| point.0 == x) {
            Some(point) => *point = (x, r, Vec::new()),
            None => points.push((x, r, Vec::new())),
        }
    }
    for (&number, &(x, height, radius)) in &holes {
        let at = match points.iter().position(|point| (point.0 - x).abs() < 1e-9) {
            Some(at) => at,
            None => {
                points.push((x, radius_at(&knots, x)?, Vec::new()));
                points.len() - 1
            }
        };
        points[at].2.push((number, height, radius));
    }
    points.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut records = mouthpiece_segments(desc, dims)?;
    for (index, (x, bore_radius, here)) in points.iter().enumerate() {
        for &(number, height, radius) in here {
            if radius >= *bore_radius {
                return Err(format!(
                    "hole-{number:02} radius {radius} is not below its bore radius {bore_radius}"
                ));
            }
            let sigma = u8::from(open.contains(&number));
            records.push(format!(
                "hole {} {} {} {sigma}",
                num(radius),
                num(height),
                num(*bore_radius)
            ));
        }
        if let Some((next_x, next_r, _)) = points.get(index + 1) {
            records.push(axial_record(*x, *bore_radius, *next_x, *next_r)?);
        }
    }
    Ok(records)
}

fn pressure_curve(desc: &Description, root: &Path) -> Lowered<Option<Vec<(i32, f64)>>> {
    let table_rel = desc.need("player.pressure_table")?;
    if table_rel == "-" {
        return Ok(None);
    }
    let group = desc.need("player.instrument_group")?;
    let level = desc.need("player.level")?;
    let path = root.join(table_rel);
    let origin = path.display().to_string();
    let mut sums: BTreeMap<i32, Vec<f64>> = BTreeMap::new();
    for row in table(&path, "frankensim.player-blowing-pressure.v1")? {
        if field(&row, "instrument", &origin)? == group && field(&row, "level", &origin)? == level {
            let pressure = parse_f64(
                field(&row, "mean_cm_h2o", &origin)?,
                &format!("{origin} mean_cm_h2o"),
            )?;
            sums.entry(midi(field(&row, "pitch", &origin)?)?)
                .or_default()
                .push(pressure);
        }
    }
    if sums.is_empty() {
        return Err(format!("no {level} rows for {group} in {origin}"));
    }
    Ok(Some(
        sums.into_iter()
            .map(|(note, values)| {
                (
                    note,
                    values.iter().fold(0.0, |sum, v| sum + v) / values.len() as f64 * CM_H2O_PA,
                )
            })
            .collect(),
    ))
}

fn pressure_for(
    desc: &Description,
    curve: Option<&[(i32, f64)]>,
    note: i32,
    reach: i32,
    step: f64,
) -> Lowered<(f64, &'static str)> {
    let Some(curve) = curve else {
        let value = desc.number("reed.operating_blowing_pa")?;
        return Ok((
            (value / step).round_ties_even() * step,
            "fixed-operating-point",
        ));
    };
    let (lo, hi) = (curve[0].0, curve[curve.len() - 1].0);
    if note < lo - reach || note > hi + reach {
        return Err(format!(
            "{} lies more than {reach} semitones outside the published pressure range",
            uiowa_name(note)
        ));
    }
    let (value, rule) = if note <= lo {
        (
            curve[0].1,
            if note < lo {
                "held-nearest"
            } else {
                "tabulated"
            },
        )
    } else if note >= hi {
        (
            curve[curve.len() - 1].1,
            if note > hi {
                "held-nearest"
            } else {
                "tabulated"
            },
        )
    } else {
        // The last window that holds the note decides, as in the Python lowerer.
        let mut found = (f64::NAN, "interpolated");
        for pair in curve.windows(2) {
            let ((ma, pa), (mb, pb)) = (pair[0], pair[1]);
            if ma <= note && note <= mb {
                let value = pa + (pb - pa) * f64::from(note - ma) / f64::from(mb - ma);
                found = (
                    value,
                    if note == ma || note == mb {
                        "tabulated"
                    } else {
                        "interpolated"
                    },
                );
            }
        }
        found
    };
    Ok(((value / step).round_ties_even() * step, rule))
}

fn gesture(desc: &Description, pressure: f64) -> Lowered<String> {
    let attack = desc.number("gesture.attack_s")?;
    let mut start = desc.number("gesture.attack_start_s")?;
    if desc.has("tongue.contact_duration_s") {
        start = round10(start - desc.number("tongue.contact_duration_s")?);
        if start < 0.0 {
            return Err("tongue.contact_duration_s exceeds gesture.attack_start_s".into());
        }
    }
    let release_at = desc.number("gesture.release_start_s")?;
    let release = desc.number("gesture.release_s")?;
    Ok(format!(
        "frankensim-gesture-schedule-v1\ncontrol_rate_hz\t{}\n\
         tracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t0e0\nevents\t2\n\
         event\t{}\t{}\tpressure-pa\t{}\nevent\t{}\t{}\tpressure-pa\t0e0\n",
        desc.need("gesture.control_rate_hz")?,
        exp(start),
        exp(attack),
        exp(pressure),
        exp(release_at),
        exp(release),
    ))
}

fn tongue_record(desc: &Description) -> Lowered<String> {
    if !desc.has("tongue.contact_area_m2") {
        return Ok(String::new());
    }
    let area_t = desc.number("tongue.contact_area_m2")?;
    let area_r = desc.number("reed.effective_area_m2")?;
    let force_n = desc.number("tongue.contact_force_n")?;
    let release_start = desc.number("gesture.attack_start_s")?;
    let contact_duration = desc.number("tongue.contact_duration_s")?;
    let unload_start = round10(release_start - contact_duration);
    let unload_s = desc.number("gesture.attack_s")?;
    let release_s = desc.number("tongue.release_s")?;
    if unload_start < 0.0 || unload_start + unload_s > release_start {
        return Err(
            "tongue contact window does not enclose the pre-release mouth-pressure ramp".into(),
        );
    }
    Ok(format!(
        "tongue_contact {} {} {} {} {}\n",
        num(area_t),
        num(area_r),
        num(force_n),
        num(release_start),
        num(release_s)
    ))
}

/// The `reed` values: literal mechanical parameters, or the per-area set of Chatziioannou and van
/// Walstijn (closing pressure K_a h0, mass m_a S_r, stiffness K_a S_r, damping g / (2 sqrt(K_a/m_a))).
fn reed_values(desc: &Description) -> Lowered<String> {
    let height = desc.need("reed.rest_opening_m")?;
    let width = desc.need("reed.width_m")?;
    if !desc.has("reed.stiffness_per_area_pa_m") {
        let mut fields = vec![height, width];
        for key in [
            "reed.closing_pressure_pa",
            "reed.mass_kg",
            "reed.stiffness_n_m",
            "reed.damping_ratio",
        ] {
            fields.push(desc.need(key)?);
        }
        return Ok(fields.join(" "));
    }
    let h0 = parse_f64(height, "reed.rest_opening_m")?;
    let area = desc.number("reed.effective_area_m2")?;
    let k_a = desc.number("reed.stiffness_per_area_pa_m")?;
    let m_a = desc.number("reed.mass_per_area_kg_m2")?;
    let g = desc.number("reed.damping_rate_1_s")?;
    let closing = k_a * h0;
    if let Some(cap) = desc.get("reed.tip_opening_cap_m") {
        if h0 > parse_f64(cap, "reed.tip_opening_cap_m")? {
            return Err(format!(
                "equilibrium reed opening {} m exceeds tip opening cap {cap} m",
                py_repr(h0)?
            ));
        }
    }
    let stiffness = k_a * area;
    let mass = m_a * area;
    let zeta = g / (2.0 * (k_a / m_a).sqrt());
    Ok(format!(
        "{} {width} {} {} {} {}",
        num(h0),
        num(closing),
        py_repr(mass)?,
        py_repr(stiffness)?,
        py_repr(zeta)?
    ))
}

fn reed_record(desc: &Description, fingering: &str) -> Lowered<String> {
    let record = reed_values(desc)?;
    let Some(stiffness) = desc.option("reed.lay_stiffness", fingering)? else {
        return Ok(record);
    };
    let mut lay = vec![stiffness];
    for key in [
        "reed.lay_exponent",
        "reed.lay_onset_m",
        "reed.lay_loss_s_per_m",
    ] {
        lay.push(desc.need(key)?);
    }
    Ok(format!("{record} {}", lay.join(" ")))
}

/// One gas declaration as written: a pack (its manifest must exist), a mixture, or T, p and RH.
fn gas_record(desc: &Description, prefix: &str, root: &Path) -> Lowered<String> {
    let pack_key = format!("{prefix}.gas_pack");
    if desc.has(&pack_key) {
        let pack = desc.need(&pack_key)?;
        let manifest = root
            .join("data/matdb/seed-v1")
            .join(pack)
            .join("manifest.tsv");
        if !manifest.is_file() {
            return Err(format!(
                "missing gas pack manifest data/matdb/seed-v1/{pack}/manifest.tsv"
            ));
        }
        return Ok(format!("pack {pack}"));
    }
    let field = |name: &str| desc.need(&format!("{prefix}.{name}"));
    if desc.has(&format!("{prefix}.composition")) {
        return Ok(format!(
            "mixture {} {} {}",
            field("temperature_k")?,
            field("pressure_pa")?,
            field("composition")?
        ));
    }
    Ok(format!(
        "{} {} {}",
        field("temperature_k")?,
        field("pressure_pa")?,
        field("relative_humidity")?
    ))
}

fn gas_field_rows(desc: &Description) -> Lowered<String> {
    if !desc.has("bore.gas_field_mode") {
        return Ok(String::new());
    }
    let species: Vec<&str> = desc
        .need("bore.gas_field_species")?
        .split_whitespace()
        .collect();
    let stations: Vec<&str> = desc
        .keys
        .iter()
        .filter(|(key, _)| key.starts_with("bore.gas_station."))
        .map(|(key, _)| desc.need(key))
        .collect::<Lowered<_>>()?;
    if stations.len() < 2 {
        return Err("bore.gas_field_mode needs at least two bore.gas_station.NN rows".into());
    }
    if let Some(short) = stations
        .iter()
        .find(|s| s.split_whitespace().count() != 2 + species.len())
    {
        return Err(format!(
            "gas station '{short}' needs X_M T_K and one fraction per species ({})",
            species.len()
        ));
    }
    let mut rows = format!(
        "gas_field {} {} {} {} {} {}\n",
        desc.need("bore.gas_field_mode")?,
        desc.need("bore.gas_field_steps")?,
        stations.len(),
        desc.need("bore.gas_field_basis")?,
        species.len(),
        species.join(" ")
    );
    for station in stations {
        rows.push_str(&format!("gas_station {station}\n"));
    }
    Ok(rows)
}

/// The rows between the bore and compile_limits, in the order the reed parser reads them.
fn tail_rows(desc: &Description, records: &[String], fingering: &str) -> Lowered<String> {
    let mut rows = gas_field_rows(desc)?;
    if desc.has("tonehole.interaction_law") {
        rows.push_str(&format!(
            "open_hole_interaction {} {}\n",
            desc.need("tonehole.interaction_law")?,
            desc.need("tonehole.interaction_epsilon")?
        ));
    }
    if let Some(lengths) = desc.option("tonehole.lengths", fingering)? {
        rows.push_str(&format!("tonehole_lengths {lengths}\n"));
    }
    if desc.has("reed.constriction_noise_seed") {
        rows.push_str(&format!(
            "constriction_noise {} {}\n",
            desc.need("reed.constriction_noise_seed")?,
            desc.need("reed.constriction_noise_duct_radius_m")?
        ));
    }
    if let Some(orifice) = desc.option("tonehole.orifice_holes", fingering)? {
        if orifice != "open" {
            return Err(format!(
                "tonehole.orifice_holes is '{orifice}'; the lowerer takes 'open' (every open hole) or '-'"
            ));
        }
        let count = records
            .iter()
            .filter(|r| r.starts_with("hole ") && r.ends_with(" 1"))
            .count();
        if count == 0 {
            return Err(format!(
                "tonehole.orifice_holes open on fingering {fingering}, which opens no hole"
            ));
        }
        rows.push_str(&format!("open_hole_orifice {count}\n"));
    }
    Ok(rows)
}

fn performance(
    desc: &Description,
    records: &[String],
    pressure: f64,
    fingering: &str,
    root: &Path,
) -> Lowered<String> {
    let rate: u32 = desc
        .need("audio.rate_hz")?
        .parse()
        .map_err(|_| "audio.rate_hz is not a whole number")?;
    let samples = (desc.number("audio.duration_s")? * f64::from(rate)).round_ties_even();
    if !(0.0..=9.0e15).contains(&samples) {
        return Err(format!("audio.duration_s gives {samples} samples"));
    }
    let samples = samples as u64;
    let mut reed = reed_record(desc, fingering)?;
    if let Some(channel) = desc.option("reed.channel_length_m", fingering)? {
        reed.push_str(&format!("\nslit_length {channel}"));
    }
    let exterior = if desc.has_prefix("exterior.") {
        format!("exterior {}\n", gas_record(desc, "exterior", root)?)
    } else {
        String::new()
    };
    let observer_law = if desc.has("observer.law") {
        format!("observer_law {}\n", desc.need("observer.law")?)
    } else {
        String::new()
    };
    Ok(format!(
        "frankensim-reed-performance-v1\naudio {rate} {samples} {}\nambient {}\nreed {reed}\n{}{}\n\
         {exterior}{observer_law}termination {}\nsegments {}\n{}\n{}compile_limits {}\nschedule\n{}",
        desc.need("audio.full_scale_pa")?,
        gas_record(desc, "ambient", root)?,
        tongue_record(desc)?,
        desc.need("observer.record")?,
        desc.need("observer.termination")?,
        records.len(),
        records.join("\n"),
        tail_rows(desc, records, fingering)?,
        desc.need("compile_limits")?,
        gesture(desc, pressure)?,
    ))
}

struct Fingering {
    name: String,
    written: String,
    concert: String,
    register: String,
    open_holes: String,
}

/// The chart's rows plus each `fingering.<name>` row, placed before the first row of higher written pitch.
fn fingerings(desc: &Description, root: &Path) -> Lowered<Vec<Fingering>> {
    let path = root.join(desc.need("fingering_chart")?);
    let origin = path.display().to_string();
    let mut rows = Vec::new();
    for row in table(&path, "frankensim.specimen-fingerings.v1")? {
        let get = |name: &str| field(&row, name, &origin).map(str::to_owned);
        rows.push(Fingering {
            name: get("fingering")?,
            written: get("written")?,
            concert: get("concert")?,
            register: get("register")?,
            open_holes: get("open_holes")?,
        });
    }
    for (key, value) in desc
        .keys
        .iter()
        .filter(|(key, _)| key.starts_with("fingering."))
    {
        let fields: Vec<&str> = value.split_whitespace().collect();
        let &[written, concert, register, open_holes] = fields.as_slice() else {
            return Err(format!("{key} needs WRITTEN CONCERT REGISTER OPEN_HOLES"));
        };
        let name = &key["fingering.".len()..];
        if rows.iter().any(|row: &Fingering| row.name == name) {
            return Err(format!(
                "{key} repeats a fingering the chart already prints"
            ));
        }
        let pitch = midi(written)?;
        let mut at = rows.len();
        for (index, row) in rows.iter().enumerate() {
            if midi(&row.written)? > pitch {
                at = index;
                break;
            }
        }
        rows.insert(
            at,
            Fingering {
                name: name.to_owned(),
                written: written.to_owned(),
                concert: concert.to_owned(),
                register: register.to_owned(),
                open_holes: open_holes.to_owned(),
            },
        );
    }
    for key in desc.keys.keys().filter(|key| key.starts_with("at.")) {
        let (name, rest) = key[3..].split_once('.').unwrap_or((&key[3..], ""));
        if !rows.iter().any(|row| row.name == name) {
            return Err(format!(
                "{key} names fingering {name}, which neither the chart nor a fingering.* row declares"
            ));
        }
        if !PER_NOTE_KEYS.contains(&rest) {
            return Err(format!(
                "{key}: {rest} is not one of the per-note keys {}",
                PER_NOTE_KEYS.join(", ")
            ));
        }
    }
    Ok(rows)
}

/// The cards (file name, text) of every fingering, or of the one named, and their index.tsv text.
pub(super) fn lower(
    desc: &Description,
    root: &Path,
    wanted: Option<&str>,
) -> Lowered<(Vec<(String, String)>, String)> {
    let specimen = root.join(desc.need("specimen_dir")?);
    let dims = dimensions(&specimen)?;
    let chart = fingerings(desc, root)?;
    let curve = pressure_curve(desc, root)?;
    let reach: i32 = desc
        .need("player.extrapolation_semitones_max")?
        .parse()
        .map_err(|_| "player.extrapolation_semitones_max is not a whole number")?;
    let step = desc.number("player.pressure_round_pa")?;
    let level = desc.need("player.level")?;
    let prefix = desc.need("output.prefix")?;
    let mut cards: Vec<(String, String)> = Vec::new();
    let mut index = vec![INDEX_HEADER.to_owned()];
    for row in chart
        .iter()
        .filter(|row| wanted.is_none_or(|name| name == row.name))
    {
        let mut open = BTreeSet::new();
        if row.open_holes != "-" {
            for hole in row.open_holes.split(',') {
                open.insert(
                    hole.parse::<i64>()
                        .map_err(|_| format!("fingering {} opens {hole:?}", row.name))?,
                );
            }
        }
        let note = midi(&row.concert)?;
        let (pressure, rule) = pressure_for(desc, curve.as_deref(), note, reach, step)?;
        let records = duct_records(desc, &dims, &open)?;
        let name = format!("{prefix}.{level}.{}.performance", uiowa_name(note));
        if cards.iter().any(|(existing, _)| *existing == name) {
            return Err(format!(
                "fingering {} lowers to {name}, which another fingering already wrote",
                row.name
            ));
        }
        let holes = records.iter().filter(|r| r.starts_with("hole ")).count();
        index.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{rule}\t{}\t{holes}\t{name}",
            row.name,
            row.written,
            uiowa_name(note),
            row.register,
            row.open_holes,
            num(pressure),
            records.len()
        ));
        cards.push((
            name,
            performance(desc, &records, pressure, &row.name, root)?,
        ));
    }
    if let Some(name) = wanted.filter(|_| cards.is_empty()) {
        return Err(format!(
            "no fingering {name} in the chart or the description's fingering.* rows"
        ));
    }
    Ok((cards, index.join("\n") + "\n"))
}

fn write_new(path: &Path, bytes: &[u8]) -> Lowered<()> {
    use std::io::Write;
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|e| format!("refusing to write {}: {e}", path.display()))?;
    file.write_all(bytes)
        .map_err(|e| format!("{}: {e}", path.display()))
}

pub(super) fn run(args: &[String], gas_options: &super::thermochemical::Options) -> Lowered<()> {
    let usage = "usage: music_render instrument DESCRIPTION.tsv FINGERING OUT.wav [wind options] | DESCRIPTION.tsv --lower OUTDIR";
    let (Some(description), Some(second), Some(third)) = (args.first(), args.get(1), args.get(2))
    else {
        return Err(usage.into());
    };
    let root = std::env::current_dir().map_err(|e| format!("working directory: {e}"))?;
    let desc = Description::parse(&read_text(Path::new(description))?, description)?;
    if second == "--lower" {
        if args.len() != 3 {
            return Err(usage.into());
        }
        let (cards, index) = lower(&desc, &root, None)?;
        let dir = Path::new(third);
        std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        for (name, text) in &cards {
            write_new(&dir.join(name), text.as_bytes())?;
        }
        write_new(&dir.join("index.tsv"), index.as_bytes())?;
        println!(
            "{{\"suite\":\"music-render\",\"verdict\":\"lowered\",\"performances\":{},\"dir\":{}}}",
            cards.len(),
            super::json::json_string(third)
        );
        return Ok(());
    }
    let (cards, _) = lower(&desc, &root, Some(second))?;
    let card: PathBuf = Path::new(third).with_extension("performance");
    if card == Path::new(third) {
        return Err(format!(
            "OUT {third} must not end in .performance: the lowered card is written there"
        ));
    }
    write_new(&card, cards[0].1.as_bytes())?;
    let card = card
        .to_str()
        .ok_or("the card path must be UTF-8")?
        .to_owned();
    let wind_args: Vec<String> = std::iter::once(card)
        .chain(args[2..].iter().cloned())
        .collect();
    super::wind_input::run_with_gas_options(&wind_args, gas_options)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn root() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
    }

    fn clarinet_text() -> String {
        read_text(&root().join("data/instruments/descriptions/clarinet/description.tsv")).unwrap()
    }

    #[test]
    fn the_clarinet_description_lowers_to_every_committed_card_and_the_index_byte_for_byte() {
        let desc = Description::parse(&clarinet_text(), "clarinet").unwrap();
        let (cards, index) = lower(&desc, &root(), None).unwrap();
        let committed = root().join("data/instruments/descriptions/clarinet/performances");
        let mut names: Vec<String> = std::fs::read_dir(&committed)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().into_string().unwrap())
            .filter(|name| name.ends_with(".performance"))
            .collect();
        names.sort();
        let mut lowered: Vec<&String> = cards.iter().map(|(name, _)| name).collect();
        lowered.sort();
        assert_eq!(lowered, names.iter().collect::<Vec<_>>());
        assert_eq!(names.len(), 24);
        for (name, text) in &cards {
            let expected = std::fs::read_to_string(committed.join(name)).unwrap();
            assert!(*text == expected, "{name} differs from the committed card");
        }
        assert_eq!(
            index,
            std::fs::read_to_string(committed.join("index.tsv")).unwrap()
        );
    }

    #[test]
    fn one_fingering_lowers_alone_to_its_card_and_an_unknown_one_is_refused() {
        let desc = Description::parse(&clarinet_text(), "clarinet").unwrap();
        let (cards, index) = lower(&desc, &root(), Some("gs3")).unwrap();
        assert_eq!(cards.len(), 1);
        assert_eq!(cards[0].0, "BbClar.mf.Gb3.performance");
        let committed = root()
            .join("data/instruments/descriptions/clarinet/performances/BbClar.mf.Gb3.performance");
        assert!(cards[0].1 == std::fs::read_to_string(committed).unwrap());
        assert_eq!(index.lines().count(), 2);
        let refused = lower(&desc, &root(), Some("zz9")).err().unwrap();
        assert!(refused.contains("no fingering zz9"), "{refused}");
    }

    #[test]
    fn a_missing_key_and_a_stray_per_note_row_are_refused_by_name() {
        let without_rate: String = clarinet_text()
            .lines()
            .filter(|line| !line.starts_with("audio.rate_hz\t"))
            .map(|l| format!("{l}\n"))
            .collect();
        let desc = Description::parse(&without_rate, "clarinet").unwrap();
        let refused = lower(&desc, &root(), Some("e3")).err().unwrap();
        assert_eq!(refused, "description is missing audio.rate_hz");
        let stray = format!(
            "{}at.zz.tonehole.lengths\t-\t-\tauthored\t-\n",
            clarinet_text()
        );
        let refused = lower(
            &Description::parse(&stray, "clarinet").unwrap(),
            &root(),
            None,
        )
        .err()
        .unwrap();
        assert!(refused.contains("names fingering zz"), "{refused}");
        let wrong = format!(
            "{}at.e3.observer.law\tnone\t-\tauthored\t-\n",
            clarinet_text()
        );
        let refused = lower(
            &Description::parse(&wrong, "clarinet").unwrap(),
            &root(),
            None,
        )
        .err()
        .unwrap();
        assert!(
            refused.contains("is not one of the per-note keys"),
            "{refused}"
        );
        let unknown = clarinet_text().replace(
            "reed.channel_length_m\t0.0014\t",
            "reed.channel_length_m\tUNKNOWN\t",
        );
        let refused = lower(
            &Description::parse(&unknown, "clarinet").unwrap(),
            &root(),
            Some("e3"),
        )
        .err()
        .unwrap();
        assert!(
            refused.contains("reed.channel_length_m is empty or UNKNOWN"),
            "{refused}"
        );
    }

    #[test]
    fn numbers_are_spelled_as_the_python_lowerer_spells_them() {
        assert_eq!(
            py_repr(3.8050000000000004e-06).unwrap(),
            "3.8050000000000004e-06"
        );
        assert_eq!(py_repr(659.0260000000001).unwrap(), "659.0260000000001");
        assert_eq!(py_repr(0.11397702454754344).unwrap(), "0.11397702454754344");
        assert_eq!(py_repr(1e16).unwrap(), "1e+16");
        assert_eq!(py_repr(1e15).unwrap(), "1000000000000000.0");
        assert_eq!(py_repr(1e-4).unwrap(), "0.0001");
        assert_eq!(py_repr(1e-5).unwrap(), "1e-05");
        assert_eq!(py_repr(-2.5).unwrap(), "-2.5");
        assert_eq!(num(0.1 + 0.2), "0.3");
        assert_eq!(num(-1e-12), "0");
        assert_eq!(num(12.0), "12");
        assert_eq!(exp(1800.0), "1.8e3");
        assert_eq!(exp(0.02), "2e-2");
        assert_eq!(exp(-0.0), "0e0");
        assert_eq!(round10(0.15 - 0.035), 0.115);
        assert_eq!(uiowa_name(midi("F#5").unwrap()), "Gb5");
        assert_eq!(uiowa_name(midi("D#3").unwrap()), "Eb3");
    }

    #[test]
    fn a_pressure_table_interpolates_holds_one_step_and_refuses_beyond() {
        let desc = Description::parse(&clarinet_text(), "clarinet").unwrap();
        let curve = [(60, 1000.0), (64, 2000.0), (67, 2600.0)];
        assert_eq!(
            pressure_for(&desc, Some(&curve), 62, 1, 1.0).unwrap(),
            (1500.0, "interpolated")
        );
        assert_eq!(
            pressure_for(&desc, Some(&curve), 64, 1, 1.0).unwrap(),
            (2000.0, "tabulated")
        );
        assert_eq!(
            pressure_for(&desc, Some(&curve), 59, 1, 1.0).unwrap(),
            (1000.0, "held-nearest")
        );
        assert_eq!(
            pressure_for(&desc, Some(&curve), 68, 1, 1.0).unwrap(),
            (2600.0, "held-nearest")
        );
        let refused = pressure_for(&desc, Some(&curve), 57, 1, 1.0).err().unwrap();
        assert!(
            refused.contains("A3 lies more than 1 semitones outside"),
            "{refused}"
        );
        assert_eq!(
            pressure_for(&desc, None, 62, 1, 1.0).unwrap(),
            (1800.0, "fixed-operating-point")
        );
    }
}
