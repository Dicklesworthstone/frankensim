//! M2.5-W4a free aluminium bar specimen (IRP bead
//! `frankensim-irp-v5-instrument-realism-pexfp.105`, G3 fixture): the clean
//! 2024-T4 beams of NASA TM X-2852 (Howell 1974, public domain), transcribed
//! from Tables I and III.
//!
//! This is a data-form and consistency check, never a physical claim. It
//! checks that the five specimen tables parse and agree with each other, that
//! the printed masses and thicknesses imply one alloy density, and that the
//! free-decay table carries both air and vacuum rows where the report prints
//! them, with the air rows no less damped. `unknowns.tsv` must keep the
//! temper substitution and the radiation-plus-viscous air loss on record.

use fs_vvreg::specimen::{
    Dimensions, KeyedTable, ReferenceState, SPECIMEN_GEOMETRY_ROOT, Specimen, SpecimenDescription,
    SpecimenTable,
};
use std::fs;
use std::path::{Path, PathBuf};

const SPEC: &str = "aluminum-beams-howell1974";
const DECAY: &str = "data/instruments/descriptions/bars/aluminum-beam-free-decay-howell1974.tsv";
const BEAMS: [&str; 5] = ["beam-40", "beam-45", "beam-46", "beam-55", "beam-56"];

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .expect("workspace root")
        .to_path_buf()
}

fn read(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn load() -> SpecimenDescription {
    let dir = root().join(SPECIMEN_GEOMETRY_ROOT).join(SPEC);
    let text = |table: SpecimenTable| read(&dir.join(table.file_name()));
    SpecimenDescription::try_new(
        Specimen::parse(&text(SpecimenTable::Specimen)).expect("specimen.tsv"),
        KeyedTable::parse(SpecimenTable::Frames, &text(SpecimenTable::Frames)).expect("frames"),
        Dimensions::parse(&text(SpecimenTable::Dimensions)).expect("dimensions"),
        KeyedTable::parse(SpecimenTable::Parts, &text(SpecimenTable::Parts)).expect("parts"),
        KeyedTable::parse(SpecimenTable::Unknowns, &text(SpecimenTable::Unknowns))
            .expect("unknowns"),
        None,
    )
    .expect("tables agree with each other")
}

fn value(d: &Dimensions, part: &str, quantity: &str) -> f64 {
    d.get(part, "prism", quantity)
        .unwrap_or_else(|| panic!("missing {part}/prism/{quantity}"))
        .value_si
}

/// `scalar` rows of a description-source table, by row name.
fn decay_rows() -> Vec<(String, f64)> {
    let text = read(&root().join(DECAY));
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some("frankensim.description-source.v1"));
    lines
        .filter(|l| l.starts_with("scalar\t"))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert_eq!(f.len(), 7, "{l}");
            assert_eq!(f[2], "howell-1974-strain-gage-clean-beams");
            assert_eq!(f[6], "constant");
            (f[1].to_owned(), f[4].parse().expect("finite scalar"))
        })
        .collect()
}

#[test]
fn howell_beam_tables_parse_and_agree() {
    let d = load();
    assert_eq!(d.specimen.id(), SPEC);
    assert_eq!(d.specimen.reference_state(), ReferenceState::Unloaded);
    assert!(
        d.specimen
            .source_ids()
            .iter()
            .any(|s| s == "nasa-howell-1974-aluminum-beams")
    );
    for beam in BEAMS {
        assert!(d.parts.row(beam).is_some(), "{beam} parts row");
        assert!(d.frames.row(beam).is_some(), "{beam} frame");
    }
    // One alloy: every beam's class mass over its measured volume gives the
    // same density to within the spread of nominal length, width and class
    // mass (2.79 to 2.85 g/cm3 as printed).
    let densities: Vec<f64> = BEAMS
        .iter()
        .map(|beam| {
            value(&d.dimensions, beam, "mass")
                / (value(&d.dimensions, beam, "length")
                    * value(&d.dimensions, beam, "width")
                    * value(&d.dimensions, beam, "thickness"))
        })
        .collect();
    let (lo, hi) = densities
        .iter()
        .fold((f64::INFINITY, 0.0_f64), |(lo, hi), &r| {
            (lo.min(r), hi.max(r))
        });
    println!(
        "{{\"case\":\"howell-beam-implied-density\",\"observed_kg_m3\":{densities:?},\"spread\":{}}}",
        hi / lo - 1.0
    );
    assert!(
        lo > 2650.0 && hi < 2900.0 && hi / lo < 1.03,
        "{densities:?}"
    );
    let unknowns: Vec<&str> = d.unknowns.keys().collect();
    for id in [
        "u-elastic-constants",
        "u-air-loss-split",
        "u-gage-and-threads",
    ] {
        assert!(unknowns.contains(&id), "{id} must stay recorded");
    }
}

#[test]
fn howell_free_decay_rows_pair_air_with_vacuum() {
    let rows = decay_rows();
    assert_eq!(rows.len(), 18);
    let get = |name: &str| {
        rows.iter()
            .find(|(n, _)| n == name)
            .unwrap_or_else(|| panic!("missing {name}"))
            .1
    };
    // Air adds loss: each printed vacuum decrement is below its air row.
    for (beam, vacuum) in [
        ("beam-40", "3.0e-4atm"),
        ("beam-45", "5.6e-4atm"),
        ("beam-46", "5.3e-4atm"),
        ("beam-56", "5.6e-4atm"),
    ] {
        let air = get(&format!("{beam}-log-decrement-1atm"));
        let low = get(&format!("{beam}-log-decrement-{vacuum}"));
        assert!(low < air && low > 0.0, "{beam}: {low} vs {air}");
        let f_air = get(&format!("{beam}-frequency-1atm"));
        let f_low = get(&format!("{beam}-frequency-{vacuum}"));
        assert!((f_air - f_low).abs() / f_air < 0.01, "{beam}");
    }
    // Specimen 55 is printed at 1 atm only.
    assert!(get("beam-55-log-decrement-1atm") > 0.0);
    assert!(
        !rows
            .iter()
            .any(|(n, _)| n.starts_with("beam-55-log-decrement-5"))
    );
}
