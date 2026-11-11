//! The duduk specimen of Maugeais and Dalmont 2024 (Acta Acustica 8, 60, CC BY 4.0), as the plain-bore card
//! `data/instruments/geometry/duduk-a-simon-maugeais-dalmont2024`.
//!
//! A data-form and consistency check, never a physical claim: the six tables parse under the repository's own
//! specimen parser and agree with each other, the hole positions are the printed ones of Figure 1 plus the exposed
//! reed (0.088 m), the bore is the plain 12 mm duct of 0.442 m the card states, the chart closes the first k holes of
//! Table 1, and the open items stay recorded in `unknowns.tsv`. No committed test loaded this specimen before; the
//! first load found two vocabulary refusals (reference_state, route) that this test now guards.

use fs_vvreg::specimen::{
    Dimensions, Fingerings, KeyedTable, ReferenceState, SPECIMEN_GEOMETRY_ROOT, Specimen, SpecimenDescription,
    SpecimenTable,
};
use std::fs;
use std::path::{Path, PathBuf};

const SPEC: &str = "duduk-a-simon-maugeais-dalmont2024";
/// Hole positions printed in Figure 1, from the body's top plane [mm].
const PRINTED_MM: [f64; 10] = [45.0, 57.3, 79.2, 111.2, 143.2, 170.7, 201.2, 231.7, 260.1, 297.6];

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).parent().and_then(Path::parent).expect("workspace root").to_path_buf()
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
        KeyedTable::parse(SpecimenTable::Unknowns, &text(SpecimenTable::Unknowns)).expect("unknowns"),
        Some(Fingerings::parse(&text(SpecimenTable::Fingerings)).expect("fingerings")),
    )
    .expect("the six tables agree with each other")
}

#[test]
fn the_six_tables_parse_agree_and_carry_the_card_s_geometry() {
    let d = load();
    assert_eq!(d.specimen.id(), SPEC);
    assert_eq!(d.specimen.reference_state(), ReferenceState::Assembled);
    assert_eq!(
        d.specimen.source_ids(),
        ["maugeais-dalmont-2024".to_owned(), "maugeais-dalmont-2024-figure-reads".to_owned()]
    );
    let value = |part: &str, feature: &str, quantity: &str| {
        d.dimensions.get(part, feature, quantity).unwrap_or_else(|| panic!("row {part}/{feature}/{quantity}")).value_si
    };
    // the plain bore: 12 mm over 0.442 m from the reed tip (0.103 m reed region + 0.339 m below it)
    assert_eq!(value("bore", "station-000", "diameter"), 0.012);
    assert_eq!(value("bore", "station-001", "diameter"), 0.012);
    let length = value("bore", "station-001", "axial-position") - value("bore", "station-000", "axial-position");
    assert!((length - 0.442).abs() < 1e-12, "bore length {length}");
    // ten holes: the printed positions plus the exposed reed, all 9.0 mm with the 5 mm wall
    for (i, printed) in PRINTED_MM.iter().enumerate() {
        let hole = format!("hole-{:02}", i + 1);
        let position = value("tone-holes", &hole, "axial-position");
        assert!((position - (0.088 + printed / 1000.0)).abs() < 1e-12, "{hole} at {position}");
        assert_eq!(value("tone-holes", &hole, "hole-diameter"), 0.009);
        assert_eq!(value("tone-holes", &hole, "chimney-height"), 0.005);
    }
    // Table 1: the chart lists the first k holes closed, k = 10 down to 0, so row i opens holes k+1..=10
    let chart = d.fingerings.as_ref().expect("fingerings").rows();
    assert_eq!(chart.len(), 11);
    for (i, row) in chart.iter().enumerate() {
        let closed = 10 - u16::try_from(i).expect("row index");
        let expected: Vec<u16> = ((closed + 1)..=10).collect();
        let open: Vec<u16> = row.open_holes.iter().map(|o| o.hole).collect();
        assert_eq!(open, expected, "{}", row.fingering);
        assert_eq!(open.len(), i, "{}", row.fingering);
        assert!(row.open_holes.iter().all(|o| o.fraction == 1.0));
    }
    // the open items stay recorded
    for id in ["u-undercut", "u-reed-face-area", "u-reed-interior", "u-fitted-reed"] {
        assert!(d.unknowns.row(id).is_some(), "unknown {id} is gone");
    }
}
