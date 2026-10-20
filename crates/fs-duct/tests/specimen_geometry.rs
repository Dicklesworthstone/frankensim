//! M8-W3 specimen geometry consumed by the TMM duct solver (IRP bead
//! `frankensim-irp-v5-instrument-realism-pexfp.197`, partial phase).
//!
//! Specimens: the published Buffet Crampon Prestige of Taillard & Sanchez
//! (CFA 2016: bore, 24 chimneys, Vandoren M30 facts) and the Moers &
//! Kergomard 2011 Yamaha Y250 lattice, whose normal fingering set is mapped
//! onto the Buffet by hole index. Every geometric value comes from
//! `data/instruments/geometry/<id>/`; this file authors none.
//!
//! Gauntlet: G0 completeness and refusals (cap and cap+1, deleted rows,
//! units, frames), resolution ladders on the measured stations and on the
//! frequency grid, the 100-cent gross-error guard, and a falsifier proving
//! the ladder observable can fail. Each test prints `IRPDATA {json}` lines
//! for `scripts/e2e/instrument_data.sh`.
//!
//! The one closure: the mouthpiece upstream of its measured exit cone is a
//! cylinder holding the rest of the measured 11.2 cm3 (`u-mouthpiece-interior`
//! names the unknown chamber). No-claim: the specimen is not an R13, and
//! undercuts and pads are unknown, so the peaks below are a gross-error guard
//! and a data-resolution check, never a pitch or same-specimen claim.

use fs_duct::{
    Duct, DuctError, Fingering, FingeringTable, HoleState, LossModel, Segment, Termination,
    input_impedance,
};
use fs_material::gas::{GasSpec, GasState};
use fs_vvreg::instrument_data::{
    ArtifactRef, Decision, INSTRUMENT_DATA_LEDGER_LOCATOR, parse_instrument_data_ledger,
};
use fs_vvreg::specimen::{
    DimensionRow, Dimensions, Fingerings, KeyedTable, SPECIMEN_GEOMETRY_ROOT, SiUnit, Sigma,
    Specimen, SpecimenDescription, SpecimenError, SpecimenErrorKind, SpecimenTable,
};
use std::f64::consts::PI;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

const BEAD: &str = "frankensim-irp-v5-instrument-realism-pexfp.197";
const BUFFET: &str = "clarinet-buffet-prestige-1981-taillard2016";
const YAMAHA: &str = "clarinet-yamaha-y250-moers2011";
/// Boehm B-flat clarinet family: 24 chimneys, hole 24 the register hole,
/// hole 1 farthest from the mouthpiece.
const HOLES: u16 = 24;
/// Family plausibility cap on any radius: catches millimetres written as metres.
const RADIUS_CAP_M: f64 = 0.05;
/// 1/5 of X-Cl-R0's 10-cent band.
const STATION_LADDER_CENTS: f64 = 2.0;
const GRID_LADDER_CENTS: f64 = 0.2;
const GROSS_ERROR_CENTS: f64 = 100.0;
/// X-Cl-R0's band: a 5 mm hole move must exceed it.
const FALSIFIER_CENTS: f64 = 10.0;
const HOLE_SHIFT_M: f64 = 0.005;
/// Lattice correspondence between the two instruments, relative to hole 24.
const CORRESPONDENCE_M: f64 = 0.005;
const CORRESPONDENCE_RADIUS_M: f64 = 0.0005;
/// Coarse search step and its upper search limit for impedance peaks.
const COARSE_HZ: f64 = 2.0;
const SEARCH_HI_HZ: f64 = 4000.0;
const LADDER_FINGERINGS: [&str; 4] = ["e3", "a3", "d4", "g5"];

// ------------------------------------------------------------------ loading

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_owned()
}
fn table_path(spec: &str, table: SpecimenTable) -> PathBuf {
    Path::new(SPECIMEN_GEOMETRY_ROOT)
        .join(spec)
        .join(table.file_name())
}
const TABLES: [SpecimenTable; 6] = [
    SpecimenTable::Specimen,
    SpecimenTable::Frames,
    SpecimenTable::Dimensions,
    SpecimenTable::Parts,
    SpecimenTable::Unknowns,
    SpecimenTable::Fingerings,
];
/// Load a specimen directory, letting a test edit any table's text first.
fn load_with(
    spec: &str,
    edit: impl Fn(SpecimenTable, String) -> String,
) -> Result<SpecimenDescription, SpecimenError> {
    let text = |table| {
        let path = root().join(table_path(spec, table));
        fs::read_to_string(&path)
            .map(|t| edit(table, t))
            .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
    };
    SpecimenDescription::try_new(
        Specimen::parse(&text(SpecimenTable::Specimen))?,
        KeyedTable::parse(SpecimenTable::Frames, &text(SpecimenTable::Frames))?,
        Dimensions::parse(&text(SpecimenTable::Dimensions))?,
        KeyedTable::parse(SpecimenTable::Parts, &text(SpecimenTable::Parts))?,
        KeyedTable::parse(SpecimenTable::Unknowns, &text(SpecimenTable::Unknowns))?,
        Some(Fingerings::parse(&text(SpecimenTable::Fingerings))?),
    )
}
fn load(spec: &str) -> SpecimenDescription {
    load_with(spec, |_, t| t).unwrap_or_else(|e| panic!("{spec}: {e}"))
}
/// The text with the first line equal to `line` removed.
fn without_line(text: &str, line: &str) -> String {
    let (mut out, mut removed) = (String::new(), false);
    for l in text.lines() {
        if !removed && l == line {
            removed = true;
            continue;
        }
        out.push_str(l);
        out.push('\n');
    }
    assert!(removed, "line not found: {line}");
    out
}

// ------------------------------------------------------------------ IRPDATA

fn json(value: &str) -> String {
    let mut out = String::from("\"");
    for c in value.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if c.is_control() => out.push_str(&format!("\\u{:04x}", u32::from(c))),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}
fn digests(spec: &str, desc: &SpecimenDescription) -> Vec<String> {
    let mut hashes = vec![
        (SpecimenTable::Specimen, desc.specimen.digest()),
        (SpecimenTable::Frames, desc.frames.digest()),
        (SpecimenTable::Dimensions, desc.dimensions.digest()),
        (SpecimenTable::Parts, desc.parts.digest()),
        (SpecimenTable::Unknowns, desc.unknowns.digest()),
    ];
    if let Some(f) = &desc.fingerings {
        hashes.push((SpecimenTable::Fingerings, f.digest()));
    }
    hashes
        .into_iter()
        .map(|(t, h)| {
            format!(
                "{{\"path\":{},\"blake3\":{}}}",
                json(&table_path(spec, t).to_string_lossy()),
                json(&h.to_hex())
            )
        })
        .collect()
}
struct Report<'a> {
    test: &'a str,
    sources: Vec<String>,
    artifacts: Vec<String>,
    rows: usize,
    refusals: usize,
    ladder: Vec<String>,
    tolerance: f64,
    observed: f64,
    start: Instant,
}
impl Report<'_> {
    fn emit(&self, verdict: &str) {
        println!(
            "IRPDATA {{\"bead\":{},\"test\":{},\"source_ids\":[{}],\"artifacts\":[{}],\"rows_checked\":{},\
             \"refusals_exercised\":{},\"ladder\":[{}],\"tolerance\":{},\"observed\":{},\"verdict\":{},\"elapsed_ms\":{}}}",
            json(BEAD),
            json(self.test),
            self.sources
                .iter()
                .map(|s| json(s))
                .collect::<Vec<_>>()
                .join(","),
            self.artifacts.join(","),
            self.rows,
            self.refusals,
            self.ladder.join(","),
            self.tolerance,
            self.observed,
            json(verdict),
            self.start.elapsed().as_millis()
        );
    }
}
fn report<'a>(test: &'a str, spec: &str, desc: &SpecimenDescription) -> Report<'a> {
    Report {
        test,
        sources: desc.specimen.source_ids().to_vec(),
        artifacts: digests(spec, desc),
        rows: desc.dimensions.rows().count(),
        refusals: 0,
        ladder: Vec::new(),
        tolerance: 0.0,
        observed: 0.0,
        start: Instant::now(),
    }
}

// ------------------------------------------------------------------ builder

/// Why the descriptions could not become a duct; never a silent default.
#[derive(Debug)]
enum BuildError {
    Spec(SpecimenError),
    Unit { what: String },
    Cap { what: String, value: f64 },
    Geometry(String),
    Duct(DuctError),
}
impl From<SpecimenError> for BuildError {
    fn from(e: SpecimenError) -> Self {
        Self::Spec(e)
    }
}
impl From<DuctError> for BuildError {
    fn from(e: DuctError) -> Self {
        Self::Duct(e)
    }
}

fn name(row: &DimensionRow) -> String {
    format!("{}/{}/{}", row.part, row.feature, row.quantity)
}
fn length(row: &DimensionRow) -> Result<f64, BuildError> {
    if row.unit != SiUnit::Metre {
        return Err(BuildError::Unit { what: name(row) });
    }
    Ok(row.value_si)
}
fn admit_radius(what: String, r: f64) -> Result<f64, BuildError> {
    if r > 0.0 && r <= RADIUS_CAP_M {
        Ok(r)
    } else {
        Err(BuildError::Cap { what, value: r })
    }
}
fn radius(row: &DimensionRow) -> Result<f64, BuildError> {
    admit_radius(name(row), length(row)? / 2.0)
}
/// Feature numbers `prefix-NN`, required to be exactly 1..=N.
fn numbered(d: &Dimensions, part: &str, prefix: &str) -> Result<Vec<u16>, BuildError> {
    let mut numbers: Vec<u16> = d
        .features(part)
        .iter()
        .map(|f| {
            f.strip_prefix(prefix)
                .and_then(|n| n.strip_prefix('-'))
                .and_then(|n| n.parse().ok())
                .ok_or_else(|| BuildError::Geometry(format!("{part}: unexpected feature {f}")))
        })
        .collect::<Result<_, _>>()?;
    numbers.sort_unstable();
    if numbers
        .iter()
        .enumerate()
        .any(|(i, n)| usize::from(*n) != i + 1)
    {
        return Err(BuildError::Geometry(format!(
            "{part}: {prefix} numbering has a gap"
        )));
    }
    Ok(numbers)
}

#[derive(Debug)]
struct Build {
    template: Duct,
    /// Hole numbers in chain (mouthpiece-first) order.
    chain_holes: Vec<u16>,
    closure_radius: f64,
}

/// The Boehm-clarinet description as a TMM chain: mouthpiece closure and
/// measured exit cone, every `stride`-th measured station (always the last)
/// as cones, and the 24 chimneys at their measured positions on the
/// interpolated bore. `shift` moves one hole along the axis (falsifier).
fn clarinet_duct(
    desc: &SpecimenDescription,
    stride: usize,
    shift: Option<(u16, f64)>,
) -> Result<Build, BuildError> {
    let d = &desc.dimensions;
    let total = length(d.require("mouthpiece", "vandoren-m30", "total-length")?)?;
    let volume_row = d.require("mouthpiece", "vandoren-m30", "internal-volume")?;
    if volume_row.unit != SiUnit::CubicMetre {
        return Err(BuildError::Unit {
            what: name(volume_row),
        });
    }
    let x45 = length(d.require("mouthpiece", "exit-station-45mm", "axial-position")?)?;
    let r45 = radius(d.require("mouthpiece", "exit-station-45mm", "diameter")?)?;
    let x0 = length(d.require("mouthpiece", "exit-station-0mm", "axial-position")?)?;
    let r0 = radius(d.require("mouthpiece", "exit-station-0mm", "diameter")?)?;
    let cone = x0 - x45;
    let cone_volume = PI * cone * (r45 * r45 + r45 * r0 + r0 * r0) / 3.0;
    let (rest, rest_volume) = (total - cone, volume_row.value_si - cone_volume);
    if !(cone > 0.0 && rest > 0.0 && rest_volume > 0.0) {
        return Err(BuildError::Geometry(
            "mouthpiece cone exceeds its measured length or volume".into(),
        ));
    }
    let closure_radius = admit_radius(
        "mouthpiece closure".into(),
        (rest_volume / (PI * rest)).sqrt(),
    )?;
    let mut segments = vec![
        Segment::Cylinder {
            radius: closure_radius,
            length: rest,
        },
        Segment::Cone {
            inlet_radius: r45,
            outlet_radius: r0,
            length: cone,
        },
    ];

    let mut stations = Vec::new();
    for n in numbered(d, "bore", "station")? {
        let f = format!("station-{n:02}");
        stations.push((
            length(d.require("bore", &f, "axial-position")?)?,
            radius(d.require("bore", &f, "diameter")?)?,
        ));
    }
    if stations.len() < 2 || stations.windows(2).any(|w| w[1].0 <= w[0].0) {
        return Err(BuildError::Geometry(
            "bore stations must strictly increase".into(),
        ));
    }
    if stations[0].0 != x0 {
        return Err(BuildError::Geometry(
            "bore does not start at the mouthpiece exit".into(),
        ));
    }
    let last = stations.len() - 1;
    let kept: Vec<(f64, f64)> = (0..=last)
        .filter(|i| i % stride == 0 || *i == last)
        .map(|i| stations[i])
        .collect();

    let numbers = numbered(d, "tone-holes", "hole")?;
    if numbers.len() != usize::from(HOLES) {
        return Err(BuildError::Geometry(format!(
            "{} holes, the family has {HOLES}",
            numbers.len()
        )));
    }
    let mut holes = Vec::new();
    for n in numbers {
        let f = format!("hole-{n:02}");
        let moved = shift.filter(|(h, _)| *h == n).map_or(0.0, |(_, dx)| dx);
        holes.push((
            length(d.require("tone-holes", &f, "axial-position")?)? + moved,
            radius(d.require("tone-holes", &f, "min-diameter")?)?,
            length(d.require("tone-holes", &f, "chimney-height")?)?,
            n,
        ));
    }
    if shift.is_none() && holes.windows(2).any(|w| w[1].0 >= w[0].0) {
        return Err(BuildError::Geometry(
            "hole-01 must be farthest from the mouthpiece".into(),
        ));
    }
    holes.sort_by(|a, b| a.0.total_cmp(&b.0));
    if holes.first().is_some_and(|h| h.0 < kept[0].0)
        || holes.last().is_some_and(|h| h.0 >= kept[kept.len() - 1].0)
    {
        return Err(BuildError::Geometry("a hole lies outside the bore".into()));
    }
    let mut next = holes.iter().peekable();
    let mut chain_holes = Vec::new();
    for w in kept.windows(2) {
        let ((xa, ra), (xb, rb)) = (w[0], w[1]);
        let (mut cursor, mut r_cursor) = (xa, ra);
        while let Some(&&(x, r_hole, h, n)) = next.peek() {
            if x >= xb {
                break;
            }
            let r_at = ra + (rb - ra) * (x - xa) / (xb - xa);
            if x > cursor {
                segments.push(Segment::Cone {
                    inlet_radius: r_cursor,
                    outlet_radius: r_at,
                    length: x - cursor,
                });
            }
            segments.push(Segment::ToneHole {
                hole_radius: r_hole,
                chimney_height: h,
                bore_radius: r_at,
                state: HoleState::Closed,
            });
            chain_holes.push(n);
            (cursor, r_cursor) = (x, r_at);
            next.next();
        }
        segments.push(Segment::Cone {
            inlet_radius: r_cursor,
            outlet_radius: rb,
            length: xb - cursor,
        });
    }
    for s in &segments {
        s.validate()?;
    }
    Ok(Build {
        template: Duct { segments },
        chain_holes,
        closure_radius,
    })
}

fn fingering_table(build: &Build, chart: &Fingerings) -> Result<FingeringTable, BuildError> {
    let mut fingerings = Vec::new();
    for row in chart.rows() {
        if let Some(o) = row
            .open_holes
            .iter()
            .find(|o| !build.chain_holes.contains(&o.hole))
        {
            return Err(BuildError::Geometry(format!(
                "{}: hole {} is not in the duct",
                row.fingering, o.hole
            )));
        }
        let holes = build
            .chain_holes
            .iter()
            .map(|n| match row.open_holes.iter().find(|o| o.hole == *n) {
                Some(o) if o.fraction >= 1.0 => HoleState::Open,
                Some(o) => HoleState::Vent(o.fraction),
                None => HoleState::Closed,
            })
            .collect();
        fingerings.push(Fingering {
            label: row.fingering.clone(),
            holes,
        });
    }
    Ok(FingeringTable::try_new(build.template.clone(), fingerings)?)
}

// ------------------------------------------------------------------ acoustics

/// Declared test conditions: dry air at 20 °C and 1 atm (the source states
/// no measurement environment, `u-environment`).
fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}
fn z_abs(duct: &Duct, air: &GasState, hz: f64) -> f64 {
    input_impedance(
        duct,
        air,
        2.0 * PI * hz,
        LossModel::Bessel,
        Termination::FlangedOpen,
    )
    .unwrap_or_else(|e| panic!("impedance at {hz} Hz: {e:?}"))
    .impedance
    .abs()
}
/// The first `count` |Z_in| maxima: a `coarse` search from 60 Hz, then a
/// 41-point local grid of step `coarse/20` and a parabolic vertex.
fn peaks(duct: &Duct, count: usize, coarse: f64) -> Vec<f64> {
    let air = air();
    let mut found = Vec::new();
    let mut f = 60.0 + coarse;
    let (mut a, mut b) = (z_abs(duct, &air, 60.0), z_abs(duct, &air, f));
    while found.len() < count && f < SEARCH_HI_HZ {
        let c = z_abs(duct, &air, f + coarse);
        if b > a && b > c {
            let h = coarse / 20.0;
            let grid: Vec<(f64, f64)> = (0..=40)
                .map(|i| {
                    let hz = f - coarse + h * f64::from(i);
                    (hz, z_abs(duct, &air, hz))
                })
                .collect();
            let i = (1..40)
                .max_by(|&i, &j| grid[i].1.total_cmp(&grid[j].1))
                .unwrap();
            let (ya, yb, yc) = (grid[i - 1].1, grid[i].1, grid[i + 1].1);
            found.push(grid[i].0 + 0.5 * h * (ya - yc) / (ya - 2.0 * yb + yc));
        }
        (a, b, f) = (b, c, f + coarse);
    }
    assert_eq!(
        found.len(),
        count,
        "only {} of {count} peaks below {SEARCH_HI_HZ} Hz",
        found.len()
    );
    found
}
fn cents(f: f64, reference: f64) -> f64 {
    1200.0 * (f / reference).log2()
}
/// Equal-tempered A4 = 440 Hz frequency of a note such as `A#3`.
fn nominal_hz(note: &str) -> f64 {
    const NAMES: [&str; 12] = [
        "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B",
    ];
    let (pc, octave) = note.split_at(note.len() - 1);
    let pc = NAMES
        .iter()
        .position(|n| *n == pc)
        .unwrap_or_else(|| panic!("note {note}"));
    let midi = 12 * (octave.parse::<i32>().expect("octave") + 1) + pc as i32;
    440.0 * 2f64.powf(f64::from(midi - 69) / 12.0)
}
/// Which impedance resonance sounds: the first in the first register, the
/// twelfth (second peak) with the register hole open. Other labels refuse.
fn sounding_peak(register: &str) -> usize {
    match register {
        "first" => 0,
        "second" => 1,
        other => panic!("register {other:?} has no declared sounding resonance"),
    }
}
fn label_duct(table: &FingeringTable, label: &str) -> Duct {
    table
        .duct(label)
        .unwrap_or_else(|e| panic!("{label}: {e:?}"))
}

// ------------------------------------------------------------------ tests

#[test]
fn buffet_tables_parse_complete_with_uncertainty() {
    let desc = load(BUFFET);
    let mut r = report(
        "buffet_tables_parse_complete_with_uncertainty",
        BUFFET,
        &desc,
    );
    let d = &desc.dimensions;
    let stations = numbered(d, "bore", "station").unwrap();
    assert_eq!(
        stations.len(),
        64,
        "Taillard 2016 Tableau 3 has 64 Buffet stations"
    );
    for n in &stations {
        for q in ["axial-position", "diameter"] {
            d.require("bore", &format!("station-{n:02}"), q).unwrap();
        }
    }
    assert_eq!(
        numbered(d, "tone-holes", "hole").unwrap().len(),
        usize::from(HOLES)
    );
    for n in 1..=HOLES {
        for q in ["axial-position", "min-diameter", "chimney-height"] {
            d.require("tone-holes", &format!("hole-{n:02}"), q).unwrap();
        }
    }
    let stated = d
        .rows()
        .filter(|row| matches!(row.sigma, Sigma::Value(_)))
        .count();
    let unstated = d.rows().filter(|row| row.sigma == Sigma::Unstated).count();
    assert_eq!(
        stated + unstated,
        d.rows().count(),
        "every row carries a sigma or `unstated`"
    );
    // The micrometer precision is stated for the bore outside the last 20 cm.
    let bell_from = d
        .require("bore", "station-64", "axial-position")
        .unwrap()
        .value_si
        - 0.2;
    for n in &stations {
        let f = format!("station-{n:02}");
        let x = d.require("bore", &f, "axial-position").unwrap().value_si;
        let sigma = d.require("bore", &f, "diameter").unwrap().sigma;
        assert_eq!(
            matches!(sigma, Sigma::Value(_)),
            x < bell_from,
            "{f}: sigma {sigma:?} at x = {x}"
        );
    }
    for id in [
        "u-model-identity",
        "u-undercut",
        "u-pad-geometry",
        "u-key-mechanics",
        "u-joint-lengths",
        "u-hole-azimuth",
        "u-hole-sigma",
        "u-bore-bell",
        "u-mouthpiece-interior",
        "u-reed",
        "u-materials",
        "u-fingering-nonnormal",
        "u-caliper-crosscheck",
        "u-player-anatomy",
    ] {
        assert!(
            desc.unknowns.row(id).is_some(),
            "absent feature must be a named unknown: {id}"
        );
    }
    let chart = desc.fingerings.as_ref().unwrap();
    assert_eq!(
        chart.rows().len(),
        21,
        "Moers 2011 normal set: 12 first- and 9 second-register notes"
    );
    for row in chart.rows() {
        let (written, concert) = (nominal_hz(&row.written), nominal_hz(&row.concert));
        assert!(
            (cents(written, concert) - 200.0).abs() < 1e-9,
            "{}: B-flat transposition",
            row.fingering
        );
        let register_open = row.open_holes.iter().any(|o| o.hole == HOLES);
        assert_eq!(
            register_open,
            row.register == "second",
            "{}: register hole state",
            row.fingering
        );
    }
    assert!(
        desc.specimen
            .get("target_relation")
            .unwrap()
            .contains("not the R13")
    );
    r.rows += chart.rows().len() + desc.unknowns.keys().count();
    r.observed = unstated as f64;
    r.emit("pass");
    println!("buffet sigma rows: {stated} stated, {unstated} unstated");
}

#[test]
fn buffet_units_and_frames_refuse_mismatches() {
    let desc = load(BUFFET);
    let mut r = report("buffet_units_and_frames_refuse_mismatches", BUFFET, &desc);
    r.tolerance = RADIUS_CAP_M;
    assert!(admit_radius("cap".into(), RADIUS_CAP_M).is_ok());
    assert!(matches!(
        admit_radius("cap+1".into(), RADIUS_CAP_M.next_up()),
        Err(BuildError::Cap { .. })
    ));
    // 14.9 mm written as 14.9 m: the cap names the row.
    let row = d_line(&desc, "bore", "station-10", "diameter");
    let mm_as_m = load_with(BUFFET, |t, text| {
        if t == SpecimenTable::Dimensions {
            text.replace(&row, &row.replace("\t0.0149\t", "\t14.9\t"))
        } else {
            text
        }
    })
    .unwrap();
    match clarinet_duct(&mm_as_m, 1, None) {
        Err(BuildError::Cap { what, value }) => {
            assert_eq!(what, "bore/station-10/diameter");
            assert!((value - 7.45).abs() < 1e-12);
        }
        other => panic!("mm-as-m must refuse at the cap: {other:?}"),
    }
    // A length declared in the wrong unit refuses by name.
    let as_volume = load_with(BUFFET, |t, text| {
        if t == SpecimenTable::Dimensions {
            text.replace(&row, &row.replace("\tm\t", "\tm3\t"))
        } else {
            text
        }
    })
    .unwrap();
    assert!(
        matches!(clarinet_duct(&as_volume, 1, None), Err(BuildError::Unit { what }) if what == "bore/station-10/diameter")
    );
    // A dimensioned part without a frame refuses.
    let frame = desc.frames.row("tone-holes").unwrap().join("\t");
    let no_frame = load_with(BUFFET, |t, text| {
        if t == SpecimenTable::Frames {
            without_line(&text, &frame)
        } else {
            text
        }
    });
    assert_eq!(
        no_frame.unwrap_err().kind,
        SpecimenErrorKind::Inconsistent("part tone-holes has no frame".to_owned())
    );
    r.refusals = 4;
    r.emit("pass");
}

/// The exact `dimensions.tsv` line of one value.
fn d_line(desc: &SpecimenDescription, part: &str, feature: &str, quantity: &str) -> String {
    let text =
        fs::read_to_string(root().join(table_path(BUFFET, SpecimenTable::Dimensions))).unwrap();
    let row = desc.dimensions.require(part, feature, quantity).unwrap();
    text.lines()
        .find(|l| {
            l.starts_with(&format!(
                "{}\t{}\t{}\t",
                row.part, row.feature, row.quantity
            ))
        })
        .unwrap()
        .to_owned()
}

#[test]
fn buffet_caliper_cross_checks_are_recorded_unavailable() {
    let desc = load(BUFFET);
    let mut r = report(
        "buffet_caliper_cross_checks_are_recorded_unavailable",
        BUFFET,
        &desc,
    );
    // No independent caliper measurement of this specimen exists: the gap is
    // named and blocks closure, rather than a check being skipped silently.
    let blocks = desc
        .unknowns
        .field("u-caliper-crosscheck", "blocks")
        .unwrap();
    assert!(blocks.contains("M8-W3a"), "{blocks}");
    // Source-internal consistency only (one source, so not independent): the
    // stated barrel length ends at the station where the bore steps down,
    // and the stated barrel bore matches the stations inside it.
    let d = &desc.dimensions;
    let barrel_len = d
        .require("barrel", "original-barrel", "useful-length")
        .unwrap()
        .value_si;
    let barrel_d = d
        .require("barrel", "original-barrel", "bore-diameter")
        .unwrap()
        .value_si;
    let mut inside = 0;
    for n in numbered(d, "bore", "station").unwrap() {
        let f = format!("station-{n:02}");
        let x = d.require("bore", &f, "axial-position").unwrap().value_si;
        let dia = d.require("bore", &f, "diameter").unwrap().value_si;
        if x <= barrel_len + 1e-12 {
            inside += 1;
            assert!(
                (dia - barrel_d).abs() <= 0.02e-3 + 1e-12,
                "{f}: {dia} vs barrel {barrel_d}"
            );
        } else if inside > 0 && usize::from(n) == inside + 1 {
            assert!(
                dia < barrel_d - 0.1e-3,
                "{f}: the bore steps down after the barrel"
            );
        }
    }
    assert_eq!(inside, 6, "the barrel spans stations 1-6");
    r.rows = inside;
    r.emit("pass-unavailable-recorded");
}

#[test]
fn buffet_builds_solver_object_without_defaults() {
    let desc = load(BUFFET);
    let mut r = report(
        "buffet_builds_solver_object_without_defaults",
        BUFFET,
        &desc,
    );
    let build = clarinet_duct(&desc, 1, None).unwrap();
    let table = fingering_table(&build, desc.fingerings.as_ref().unwrap()).unwrap();
    assert_eq!(table.labels().len(), 21);
    assert_eq!(
        build.chain_holes.first(),
        Some(&HOLES),
        "the register hole is nearest the mouthpiece"
    );
    println!(
        "mouthpiece closure radius {:.4} mm (u-mouthpiece-interior)",
        build.closure_radius * 1e3
    );
    // Deleting any one row the solver object consumes refuses, naming it.
    let mut refused = 0;
    for row in desc.dimensions.rows().filter(|row| row.part != "barrel") {
        let line = d_line(&desc, &row.part, &row.feature, &row.quantity);
        let edited = load_with(BUFFET, |t, s| {
            if t == SpecimenTable::Dimensions {
                without_line(&s, &line)
            } else {
                s
            }
        })
        .unwrap();
        match clarinet_duct(&edited, 1, None) {
            Err(BuildError::Spec(e)) => {
                assert_eq!(
                    e.kind,
                    SpecimenErrorKind::Missing(format!("dimension {}", name(row)))
                );
                refused += 1;
            }
            other => panic!("deleting {} must refuse by name: {other:?}", name(row)),
        }
    }
    assert_eq!(refused, 64 * 2 + 24 * 3 + 6);
    // A chart hole the duct does not carry refuses too.
    let bad_chart = load_with(BUFFET, |t, s| {
        if t == SpecimenTable::Fingerings {
            s.replacen("\t1,3,5\t", "\t1,3,5,25\t", 1)
        } else {
            s
        }
    })
    .unwrap();
    assert!(fingering_table(&build, bad_chart.fingerings.as_ref().unwrap()).is_err());
    r.refusals = refused + 1;
    r.emit("pass");
}

/// The M8-W3a closure gate, kept exactly as specified and RED on this source:
/// measured 2026-10-20 (csd -> hz2), the all-closed `e3` peaks move up to
/// 5.46 c between every-2nd and every published station (peak 1: +6.6 c then
/// +5.1 c, not converging), because the sparse bell stations dominate when no
/// hole is open; open-hole fingerings' low peaks agree within 0.7 c. The
/// published table therefore limits bell-dominated peaks above 2 c
/// (`u-bore-bell`). Run with `--ignored`; never widen the band to pass.
#[test]
#[ignore = "red closure gate: published Buffet bell stations limit the all-closed peaks to ~5.5 c (2 c required); pexfp.197 stays open"]
fn buffet_station_ladder_converges() {
    let desc = load(BUFFET);
    let mut r = report("buffet_station_ladder_converges", BUFFET, &desc);
    r.tolerance = STATION_LADDER_CENTS;
    let chart = desc.fingerings.as_ref().unwrap();
    let rungs: Vec<(usize, FingeringTable)> = [4, 2, 1]
        .into_iter()
        .map(|stride| {
            let build = clarinet_duct(&desc, stride, None).unwrap();
            (stride, fingering_table(&build, chart).unwrap())
        })
        .collect();
    let mut worst = 0.0f64;
    for label in LADDER_FINGERINGS {
        let found: Vec<Vec<f64>> = rungs
            .iter()
            .map(|(_, t)| peaks(&label_duct(t, label), 6, COARSE_HZ))
            .collect();
        for k in 0..6 {
            let (coarse, mid, fine) = (found[0][k], found[1][k], found[2][k]);
            let last = cents(mid, fine);
            worst = worst.max(last.abs());
            r.ladder.push(format!(
                "{{\"fingering\":{},\"peak\":{},\"stride4_hz\":{coarse:.4},\"stride2_hz\":{mid:.4},\"stride1_hz\":{fine:.4},\
                 \"c_4_to_2\":{:.3},\"c_2_to_1\":{last:.3}}}",
                json(label),
                k + 1,
                cents(coarse, mid)
            ));
        }
    }
    r.observed = worst;
    let verdict = if worst <= STATION_LADDER_CENTS {
        "pass"
    } else {
        "fail"
    };
    r.emit(verdict);
    assert!(
        worst <= STATION_LADDER_CENTS,
        "station ladder: {worst:.3} c between every-2nd and every station exceeds {STATION_LADDER_CENTS} c"
    );
}

#[test]
fn buffet_frequency_grid_ladder_converges() {
    let desc = load(BUFFET);
    let mut r = report("buffet_frequency_grid_ladder_converges", BUFFET, &desc);
    r.tolerance = GRID_LADDER_CENTS;
    let build = clarinet_duct(&desc, 1, None).unwrap();
    let table = fingering_table(&build, desc.fingerings.as_ref().unwrap()).unwrap();
    let duct = label_duct(&table, "e3");
    let found: Vec<Vec<f64>> = [2.0 * COARSE_HZ, COARSE_HZ, 0.5 * COARSE_HZ]
        .into_iter()
        .map(|step| peaks(&duct, 6, step))
        .collect();
    let mut worst = 0.0f64;
    for k in 0..6 {
        let last = cents(found[1][k], found[2][k]);
        worst = worst.max(last.abs());
        r.ladder.push(format!(
            "{{\"peak\":{},\"df4_hz\":{:.5},\"df2_hz\":{:.5},\"df1_hz\":{:.5},\"c_2_to_1\":{last:.4}}}",
            k + 1,
            found[0][k],
            found[1][k],
            found[2][k]
        ));
    }
    r.observed = worst;
    r.emit(if worst <= GRID_LADDER_CENTS {
        "pass"
    } else {
        "fail"
    });
    assert!(
        worst <= GRID_LADDER_CENTS,
        "grid ladder {worst:.4} c > {GRID_LADDER_CENTS} c"
    );
}

#[test]
fn buffet_first_peak_gross_error_guard() {
    let desc = load(BUFFET);
    let mut r = report("buffet_first_peak_gross_error_guard", BUFFET, &desc);
    r.tolerance = GROSS_ERROR_CENTS;
    let build = clarinet_duct(&desc, 1, None).unwrap();
    let chart = desc.fingerings.as_ref().unwrap();
    let table = fingering_table(&build, chart).unwrap();
    let mut worst = 0.0f64;
    let mut failures = Vec::new();
    for row in chart.rows() {
        let k = sounding_peak(&row.register);
        let f = peaks(&label_duct(&table, &row.fingering), k + 1, COARSE_HZ)[k];
        let nominal = nominal_hz(&row.concert);
        let c = cents(f, nominal);
        worst = worst.max(c.abs());
        r.ladder.push(format!(
            "{{\"fingering\":{},\"concert\":{},\"peak\":{},\"model_hz\":{f:.3},\"nominal_hz\":{nominal:.3},\"cents\":{c:.1}}}",
            json(&row.fingering),
            json(&row.concert),
            k + 1
        ));
        if c.abs() > GROSS_ERROR_CENTS {
            failures.push(format!("{} {c:+.1} c", row.fingering));
        }
    }
    r.observed = worst;
    r.emit(if failures.is_empty() { "pass" } else { "fail" });
    assert!(
        failures.is_empty(),
        "probable hole-state or unit errors (investigate the source, never edit geometry): {failures:?}"
    );
}

#[test]
fn falsifier_buffet_hole_shift_moves_the_peak() {
    let desc = load(BUFFET);
    let mut r = report("falsifier_buffet_hole_shift_moves_the_peak", BUFFET, &desc);
    r.tolerance = FALSIFIER_CENTS;
    let chart = desc.fingerings.as_ref().unwrap();
    // a3 opens holes 1, 3 and 5; hole 5 is its first open hole.
    let base = fingering_table(&clarinet_duct(&desc, 1, None).unwrap(), chart).unwrap();
    let moved = fingering_table(
        &clarinet_duct(&desc, 1, Some((5, HOLE_SHIFT_M))).unwrap(),
        chart,
    )
    .unwrap();
    let f0 = peaks(&label_duct(&base, "a3"), 1, COARSE_HZ)[0];
    let f1 = peaks(&label_duct(&moved, "a3"), 1, COARSE_HZ)[0];
    let shift = cents(f1, f0);
    r.observed = shift.abs();
    r.ladder.push(format!(
        "{{\"base_hz\":{f0:.4},\"moved_hz\":{f1:.4},\"cents\":{shift:.2}}}"
    ));
    r.emit(if shift.abs() > FALSIFIER_CENTS {
        "pass"
    } else {
        "fail"
    });
    assert!(
        shift < -FALSIFIER_CENTS,
        "moving hole 5 by 5 mm toward the bell must flatten a3 by > 10 c: {shift:.2}"
    );
}

#[test]
fn buffet_hole_index_corresponds_to_the_fingering_source() {
    let (buffet, yamaha) = (load(BUFFET), load(YAMAHA));
    let mut r = report(
        "buffet_hole_index_corresponds_to_the_fingering_source",
        BUFFET,
        &buffet,
    );
    r.tolerance = CORRESPONDENCE_M;
    let x = |desc: &SpecimenDescription, n: u16| {
        desc.dimensions
            .require("tone-holes", &format!("hole-{n:02}"), "axial-position")
            .unwrap()
            .value_si
    };
    let radius_b = |n: u16| {
        buffet
            .dimensions
            .require("tone-holes", &format!("hole-{n:02}"), "min-diameter")
            .unwrap()
            .value_si
            / 2.0
    };
    let radius_y = |n: u16| {
        yamaha
            .dimensions
            .require("tone-holes", &format!("hole-{n:02}"), "hole-radius")
            .unwrap()
            .value_si
    };
    // Positions relative to each instrument's register hole remove the two
    // different axis origins (barrel entry vs reed tip).
    let rel = |desc: &SpecimenDescription, n: u16| x(desc, n) - x(desc, HOLES);
    let mut worst: f64 = 0.0;
    for n in 1..=HOLES {
        let dx = rel(&buffet, n) - rel(&yamaha, n);
        let dr = radius_b(n) - radius_y(n);
        worst = worst.max(dx.abs());
        r.ladder.push(format!(
            "{{\"hole\":{n},\"dx_mm\":{:.2},\"dr_mm\":{:.2}}}",
            dx * 1e3,
            dr * 1e3
        ));
        assert!(dx.abs() <= CORRESPONDENCE_M, "hole {n}: {:.2} mm", dx * 1e3);
        assert!(
            dr.abs() <= CORRESPONDENCE_RADIUS_M,
            "hole {n}: radius {:.2} mm",
            dr * 1e3
        );
    }
    // The check can fail: an off-by-one index breaks it.
    let shifted = (1..HOLES)
        .map(|n| (rel(&buffet, n) - rel(&yamaha, n + 1)).abs())
        .fold(0.0, f64::max);
    assert!(
        shifted > CORRESPONDENCE_M,
        "an index shift must break the correspondence ({shifted})"
    );
    r.observed = worst;
    r.refusals = 1;
    r.sources
        .extend(yamaha.specimen.source_ids().iter().cloned());
    r.artifacts.extend(digests(YAMAHA, &yamaha));
    r.emit("pass");
}

#[test]
fn buffet_ledger_rows_cover_every_source() {
    let text = fs::read_to_string(root().join(INSTRUMENT_DATA_LEDGER_LOCATOR)).unwrap();
    let ledger = parse_instrument_data_ledger(&text).expect("committed ledger parses");
    let mut checked = 0;
    for spec in [BUFFET, YAMAHA] {
        let desc = load(spec);
        let mut cited: Vec<&str> = desc
            .specimen
            .source_ids()
            .iter()
            .map(String::as_str)
            .collect();
        cited.extend(desc.dimensions.rows().map(|r| r.source_id.as_str()));
        cited.extend(
            desc.fingerings
                .iter()
                .flat_map(|f| f.rows().iter().map(|r| r.source_id.as_str())),
        );
        cited.sort_unstable();
        cited.dedup();
        for source in cited {
            let row = ledger
                .by_source(source)
                .unwrap_or_else(|| panic!("{spec}: no ledger row for {source}"));
            assert!(
                matches!(row.decision, Decision::Admit | Decision::AdmitFactsOnly),
                "{source}: {:?}",
                row.decision
            );
            assert!(
                row.artifacts
                    .contains(&ArtifactRef::Geometry(spec.to_owned())),
                "{source} must name geometry:{spec}"
            );
            checked += 1;
        }
        let mut r = report("buffet_ledger_rows_cover_every_source", spec, &desc);
        r.rows = checked;
        r.emit("pass");
    }
    assert_eq!(checked, 3, "Buffet cites two sources, the Yamaha one");
}
