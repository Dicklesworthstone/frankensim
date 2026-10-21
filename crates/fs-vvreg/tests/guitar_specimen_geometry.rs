//! M8-W3d classical guitar geometry specimen (IRP bead
//! `frankensim-irp-v5-instrument-realism-pexfp.203`): the measured plan of a
//! 1964 Madrid classical-type guitar by Zimmer and Mores (Zenodo
//! 10.5281/zenodo.11383935, CC BY 4.0), transcribed as printed text.
//!
//! This is a data-form and consistency check, never a physical claim. It
//! checks that the six specimen tables parse and agree with each other, that
//! printed values which must agree with each other do, that the point maps
//! are well-formed and lie inside the body they were read from, and that the
//! identity discrepancy and the not-the-Iowa-guitar boundary stay recorded in
//! `unknowns.tsv`. The specimen is NOT the UIowa Raimundo 118.

use fs_vvreg::specimen::{
    Dimensions, KeyedTable, SPECIMEN_GEOMETRY_ROOT, Specimen, SpecimenDescription, SpecimenTable,
};
use std::fs;
use std::path::{Path, PathBuf};

const SPEC: &str = "guitar-esteso-1964-zimmer2024";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .expect("workspace root")
        .to_path_buf()
}

fn read(file: &str) -> String {
    let path = root().join(SPECIMEN_GEOMETRY_ROOT).join(SPEC).join(file);
    fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn load() -> SpecimenDescription {
    let text = |table: SpecimenTable| read(table.file_name());
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

/// Rows of a point table: schema line and header checked, fields split.
fn points(file: &str, header: &str, columns: usize) -> Vec<Vec<String>> {
    let text = read(file);
    let mut lines = text.lines();
    assert!(
        lines
            .next()
            .is_some_and(|l| l.starts_with("frankensim.specimen-"))
    );
    assert_eq!(lines.next(), Some(header), "{file} header");
    lines
        .map(|l| {
            let row: Vec<String> = l.split('\t').map(str::to_owned).collect();
            assert_eq!(row.len(), columns, "{file}: {l}");
            row
        })
        .collect()
}

fn value(d: &Dimensions, part: &str, feature: &str, quantity: &str) -> f64 {
    d.get(part, feature, quantity)
        .unwrap_or_else(|| panic!("missing {part}/{feature}/{quantity}"))
        .value_si
}

#[test]
fn esteso_tables_parse_and_agree() {
    let d = load();
    assert_eq!(d.specimen.id(), SPEC);
    assert_eq!(
        d.specimen.source_ids(),
        [
            "zimmer-mores-2024-esteso1964-plan",
            "zimmer-mores-2024-esteso1964-plan-positions"
        ]
    );
    assert!(d.dimensions.rows().count() >= 60);
    for part in d.dimensions.parts() {
        assert!(d.frames.row(part).is_some(), "frame for {part}");
        assert!(d.parts.row(part).is_some(), "parts row for {part}");
    }
    let spec = d.specimen.get("target_relation").expect("target_relation");
    assert!(spec.contains("NOT the UIowa Raimundo 118"));
    for id in ["u-specimen-identity", "u-not-iowa-guitar", "u-materials"] {
        assert!(d.unknowns.row(id).is_some(), "unknown {id} stays recorded");
    }
}

#[test]
fn printed_values_that_must_agree_do() {
    let dims = load().dimensions;
    // Ring width is printed 17 and the two radii R46 and R63.
    let ring = value(&dims, "top", "soundhole-annulus", "outer-radius")
        - value(&dims, "top", "soundhole-annulus", "inner-radius");
    assert!((ring - value(&dims, "top", "soundhole-annulus", "width")).abs() < 1e-12);
    // The hole lies inside the inner ring, the lower bout is wider than the upper bout and the
    // waist narrower than both.
    assert!(
        value(&dims, "top", "soundhole", "radius")
            < value(&dims, "top", "soundhole-annulus", "inner-radius")
    );
    let (lower, upper, waist) = (
        value(&dims, "top", "outline", "lower-bout-width"),
        value(&dims, "top", "outline", "upper-bout-width"),
        value(&dims, "top", "outline", "waist-width"),
    );
    assert!(lower > upper && upper > waist);
    // Back bars are not perpendicular to the axis: both printed sides exist and differ, and every
    // gap is smaller than the printed overall length.
    for pair in ["brace-1-to-2", "brace-2-to-3", "brace-3-to-4"] {
        let treble = value(&dims, "back-braces", pair, "spacing-treble-side");
        let bass = value(&dims, "back-braces", pair, "spacing-bass-side");
        assert!(treble > 0.0 && bass > 0.0 && (treble - bass).abs() > 1e-4);
        assert!(treble < value(&dims, "back", "outline", "overall-length"));
    }
    // Printed angles in radians equal their printed degrees.
    let angle = |feature: &str, quantity: &str| value(&dims, "top-braces", feature, quantity);
    assert!(
        (angle("layout-angle-treble-fan", "angle-to-transverse-line") - 83.0_f64.to_radians())
            .abs()
            < 1e-4
    );
    assert!(
        (angle("layout-angle-bass-fan", "angle-to-transverse-line") - 82.0_f64.to_radians()).abs()
            < 1e-4
    );
    assert!((angle("layout-angle-diagonal", "angle-to-axis") - 69.0_f64.to_radians()).abs() < 1e-4);
    assert!(
        (value(&dims, "neck", "headstock-break", "angle") - 18.0_f64.to_radians()).abs() < 1e-4
    );
}

#[test]
fn thickness_point_maps_are_well_formed_and_inside_their_bodies() {
    let header = "point\tx_mm\ty_mm\tthickness_mm";
    let top = points("top-thickness.tsv", header, 4);
    let back = points("back-thickness.tsv", header, 4);
    let side = points("side-thickness.tsv", header, 4);
    assert_eq!((top.len(), back.len(), side.len()), (54, 55, 9));
    for (rows, x_max, y_max) in [
        (&top, 489.4, 183.6),
        (&back, 495.0, 184.0),
        (&side, 489.4, 120.0),
    ] {
        for row in rows {
            let [x, y, t] = [&row[1], &row[2], &row[3]].map(|s| s.parse::<f64>().expect("number"));
            assert!((0.0..=x_max).contains(&x), "x {x}");
            assert!(y.abs() <= y_max, "y {y}");
            assert!((1.3..=3.0).contains(&t), "thickness {t}");
        }
    }
    // Plates are thicker than ribs at every printed point.
    let min_plate = top
        .iter()
        .chain(&back)
        .map(|r| r[3].parse::<f64>().expect("number"))
        .fold(f64::INFINITY, f64::min);
    let max_rib = side
        .iter()
        .map(|r| r[3].parse::<f64>().expect("number"))
        .fold(0.0, f64::max);
    assert!(min_plate > max_rib);
    let ids: Vec<&str> = top.iter().map(|r| r[0].as_str()).collect();
    assert_eq!(ids.first(), Some(&"T-01"));
    assert_eq!(
        ids.len(),
        ids.iter().collect::<std::collections::BTreeSet<_>>().len()
    );
}

#[test]
fn vector_reads_agree_with_the_printed_numbers() {
    let dims = load().dimensions;
    let header = "path\tindex\tx_mm\ty_mm";
    let outline = points("top-outline-vector.tsv", header, 4);
    let of = |name: &str| -> Vec<(f64, f64)> {
        outline
            .iter()
            .filter(|r| r[0] == name)
            .map(|r| (r[2].parse().expect("x"), r[3].parse().expect("y")))
            .collect()
    };
    let extent = |pts: &[(f64, f64)], axis: usize| {
        let v = |p: &(f64, f64)| if axis == 0 { p.0 } else { p.1 };
        let max = pts.iter().map(v).fold(f64::NEG_INFINITY, f64::max);
        let min = pts.iter().map(v).fold(f64::INFINITY, f64::min);
        max - min
    };
    // The plot is 1:1: the vector outline height equals the printed lower-bout width.
    let top = of("outline");
    assert_eq!(top.len(), 486);
    let printed = 1000.0 * value(&dims, "top", "outline", "lower-bout-width");
    assert!(
        (extent(&top, 1) - printed).abs() < 0.05,
        "{}",
        extent(&top, 1)
    );
    // The sound hole circle read from the vectors has the printed radius.
    let hole = of("soundhole-circle");
    let radius = 1000.0 * value(&dims, "top", "soundhole", "radius");
    assert!((extent(&hole, 0) / 2.0 - radius).abs() < 0.1);
    assert!((extent(&hole, 1) / 2.0 - radius).abs() < 0.1);
    // Fret lines: the nut line and ten fret lines coincide with equal temperament for 650 mm to 0.1 mm;
    // other lines (a short sixth line, a string tick line, a body-edge line) are 5 mm or more away.
    let frets = points(
        "fret-lines-vector.tsv",
        "line\tx_pt\tmm_from_nut_line\twidth_mm\tet_mm_for_650_scale\tdifference_mm",
        6,
    );
    let diff: Vec<f64> = frets.iter().map(|r| r[5].parse().expect("diff")).collect();
    assert_eq!(diff.iter().filter(|d| d.abs() < 0.1).count(), 11);
    assert!(diff.iter().any(|d| d.abs() > 5.0));
    let back = points("back-outline-vector.tsv", header, 4);
    assert_eq!(back.len(), 570);
}

#[test]
fn brace_strips_read_from_vectors_match_the_printed_widths_and_heights() {
    let dims = load().dimensions;
    let strips = points(
        "top-braces-vector.tsv",
        "strip\tfeature\tx0_mm\ty0_mm\tx1_mm\ty1_mm\tlength_mm\tangle_deg\twidth_mm\tprinted_width_match_mm\tannotated_height_mm\tnote",
        12,
    );
    assert_eq!(strips.len(), 12);
    let mut features = std::collections::BTreeSet::new();
    for row in &strips {
        let width: f64 = row[8].parse().expect("width");
        // Every vector strip width equals one of the printed widths to 0.05 mm.
        let printed: f64 = row[9].parse().expect("a printed width matched");
        assert!(
            (width - printed).abs() < 0.05,
            "{}: {width} vs {printed}",
            row[0]
        );
        let height: f64 = row[10].parse().expect("height");
        let feature = row[1].as_str();
        features.insert(feature.to_owned());
        // The dimension table carries the same printed numbers for the same feature.
        let w = 1000.0 * value(&dims, "top-braces", feature, "width");
        let h = 1000.0 * value(&dims, "top-braces", feature, "annotated-height");
        assert!(
            (w - printed).abs() < 1e-6,
            "{feature}: width {w} vs {printed}"
        );
        assert!(
            (h - height).abs() < 1e-6,
            "{feature}: height {h} vs {height}"
        );
    }
    // Eleven bars; the slanted bar is read as two collinear runs, so twelve strips.
    assert_eq!(features.len(), 11);
    let lines = points(
        "top-dashed-lines-vector.tsv",
        "line\tx0_mm\ty0_mm\tx1_mm\ty1_mm\tlength_mm\tangle_deg\tmax_deviation_mm\tdashes",
        9,
    );
    assert_eq!(lines.len(), 35);
}

#[test]
fn signature_modes_and_annotations_are_well_formed() {
    let modes = points(
        "signature-modes.tsv",
        "mode\tfrequency_hz\tdescription\tsource_id",
        4,
    );
    let names: Vec<&str> = modes.iter().map(|r| r[0].as_str()).collect();
    assert_eq!(names, ["A0", "T1", "T2", "T3", "T4", "A1"]);
    let hz: Vec<f64> = modes.iter().map(|r| r[1].parse().expect("hz")).collect();
    assert_eq!(hz, [98.0, 200.0, 240.0, 277.0, 317.0, 408.0]);
    let annotations = points(
        "printed-annotations.tsv",
        "view\ttext\tvertical\theight_label\tx_pt\ty_pt\tx_mm\ty_mm",
        8,
    );
    assert!(annotations.len() >= 100);
    assert!(annotations.iter().any(|r| r[1] == "367,02"));
}
