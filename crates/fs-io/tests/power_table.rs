//! Quarantine battery for tabular component-power ingestion (bead f85xj.11.7).
//!
//! The headline case is the milliwatt hazard: the SAME digits under a `mW`
//! declaration and a `W` declaration must produce wattages a factor of 1000
//! apart, with the declared unit and the applied factor both recorded. That is
//! the error the bead exists to prevent, and it is silent in every system that
//! infers the unit from a column header.
//!
//! Every other test is a refusal. The house rule from the PLY face-list model
//! is: assert the variant, assert the position exactly, and assert the message
//! names the offending value. Here "position" is BOTH the physical source line
//! and the data-row ordinal, because a spreadsheet export is edited by a human
//! in a text editor.

use fs_io::power_table::{POWER_TABLE_SEMANTICS_VERSION, PowerUnit, parse_power_table};

const GOOD: &str = "\
#power-unit: W
#declared-total: 4.2
component,power,half_width,confidence
die,2.5,0.1,0.95
regulator,1.7,,
";

fn expect_refusal(input: &str) -> fs_io::power_table::PowerTableError {
    parse_power_table(input).expect_err("input must refuse")
}

/// The happy path, and the seam the scenario layer consumes.
#[test]
fn a_well_formed_table_is_admitted_with_its_provenance() {
    let table = parse_power_table(GOOD).expect("well-formed table");
    assert_eq!(table.rows().len(), 2);

    let die = &table.rows()[0];
    assert_eq!(die.name(), "die");
    assert!((die.watts() - 2.5).abs() < 1e-12);
    assert_eq!(die.half_width_w(), Some(0.1));
    assert_eq!(die.confidence(), Some(0.95));
    // The physical line, not the record ordinal: header is line 3, so row 1 is 4.
    assert_eq!(die.line(), 4);

    let regulator = &table.rows()[1];
    assert_eq!(regulator.name(), "regulator");
    assert_eq!(
        regulator.half_width_w(),
        None,
        "an absent interval stays absent"
    );
    assert_eq!(regulator.confidence(), None);
    assert_eq!(regulator.line(), 5);

    // The parser carries names and watts. It resolves NOTHING to geometry —
    // binding a name to a vertex set is the scenario layer's job, and a parser
    // that did it would need a mesh.
    let receipt = table.receipt();
    assert_eq!(receipt.semantics(), POWER_TABLE_SEMANTICS_VERSION);
    assert_eq!(receipt.declared_unit(), PowerUnit::Watts);
    assert!((receipt.watts_factor() - 1.0).abs() < 1e-12);
    assert_eq!(receipt.rows(), 2);
    assert_eq!(receipt.declared_total_w(), Some(4.2));

    let json = receipt.to_json();
    assert!(json.contains("\"kind\":\"power-table-receipt\""), "{json}");
    assert!(json.contains("\"declared_unit\":\"W\""), "{json}");
    assert!(json.contains("\"trust\":\"quarantined\""), "{json}");
    assert!(
        json.contains("\"source_hash\""),
        "provenance binds the exact bytes: {json}"
    );
}

/// THE hazard: identical digits, different declared unit, 1000x apart — and
/// the receipt records which one was applied.
#[test]
fn the_declared_unit_changes_the_wattage_by_exactly_its_factor() {
    let body = "component,power\ndie,1200\n";
    let in_watts = parse_power_table(&format!("#power-unit: W\n{body}")).expect("watts");
    let in_milliwatts = parse_power_table(&format!("#power-unit: mW\n{body}")).expect("milliwatts");
    let in_kilowatts = parse_power_table(&format!("#power-unit: kW\n{body}")).expect("kilowatts");

    assert!((in_watts.rows()[0].watts() - 1200.0).abs() < 1e-9);
    assert!((in_milliwatts.rows()[0].watts() - 1.2).abs() < 1e-9);
    assert!((in_kilowatts.rows()[0].watts() - 1_200_000.0).abs() < 1e-3);

    // A thousandfold difference from one directive line. This is exactly the
    // error that is silent when a reader assumes watts from a column header.
    assert!((in_watts.rows()[0].watts() / in_milliwatts.rows()[0].watts() - 1000.0).abs() < 1e-6);

    // The applied factor is auditable, not invisible.
    assert_eq!(
        in_milliwatts.receipt().declared_unit(),
        PowerUnit::Milliwatts
    );
    assert!((in_milliwatts.receipt().watts_factor() - 1e-3).abs() < 1e-18);
    assert!(
        in_milliwatts
            .receipt()
            .to_json()
            .contains("\"declared_unit\":\"mW\"")
    );

    // Half widths convert with the same factor as the values they bound.
    let with_interval = parse_power_table(
        "#power-unit: mW\ncomponent,power,half_width,confidence\ndie,1200,50,0.9\n",
    )
    .expect("interval table");
    assert!((with_interval.rows()[0].half_width_w().expect("half width") - 0.05).abs() < 1e-12);
}

/// An undeclared or unrecognized unit refuses rather than assuming watts.
#[test]
fn an_undeclared_or_unknown_unit_refuses() {
    let missing = expect_refusal("component,power\ndie,2.5\n");
    assert_eq!(missing.field, "power-unit");
    assert_eq!(missing.line, 0, "a whole-document refusal has no line");
    assert!(
        missing.reason.contains("assume watts"),
        "the refusal must say why: {}",
        missing.reason
    );
    assert!(missing.reason.contains("W, mW, kW"), "{}", missing.reason);

    let unknown = expect_refusal("#power-unit: horsepower\ncomponent,power\ndie,2.5\n");
    assert_eq!(unknown.field, "power-unit");
    assert_eq!(unknown.line, 1, "the directive line is named");
    assert!(
        unknown.reason.contains("horsepower"),
        "must name the offender: {}",
        unknown.reason
    );

    // Declaring it twice is ambiguous, not a last-one-wins.
    let twice = expect_refusal("#power-unit: W\n#power-unit: mW\ncomponent,power\ndie,1\n");
    assert_eq!((twice.field, twice.line), ("power-unit", 2));
}

/// A malformed cell names BOTH the physical line and the data row, and quotes
/// the offending text.
#[test]
fn a_malformed_cell_names_the_line_the_row_and_the_offender() {
    let input = "\
#power-unit: W

# a comment that still advances the line counter
component,power
die,2.5
regulator,abc
";
    let error = expect_refusal(input);
    // Header is line 4; regulator is data row 2, hence line 6.
    assert_eq!(
        (error.line, error.row),
        (6, Some(2)),
        "comments and blanks advance the line counter: {error:?}"
    );
    assert!(
        error.reason.contains("abc"),
        "must name the offender: {}",
        error.reason
    );
    assert!(error.reason.contains("power"), "{}", error.reason);
    // The rendered message carries both coordinates.
    let rendered = error.to_string();
    assert!(rendered.contains("line 6"), "{rendered}");
    assert!(rendered.contains("data row 2"), "{rendered}");
}

/// Missing required columns, empty names, and negative power all refuse.
#[test]
fn structural_faults_refuse_by_name() {
    let missing_column = expect_refusal("#power-unit: W\ncomponent\ndie\n");
    assert_eq!(missing_column.field, "header");
    assert!(
        missing_column.reason.contains("power"),
        "{}",
        missing_column.reason
    );

    let empty_name = expect_refusal("#power-unit: W\ncomponent,power\n,2.5\n");
    assert_eq!(empty_name.line, 3);
    assert!(
        empty_name.reason.contains("empty") || empty_name.reason.contains("required"),
        "{}",
        empty_name.reason
    );

    // Negative dissipation is outside the declared column range.
    let negative = expect_refusal("#power-unit: W\ncomponent,power\ndie,-2.5\n");
    assert_eq!((negative.line, negative.row), (3, Some(1)));
    assert!(negative.reason.contains("-2.5"), "{}", negative.reason);

    // A unit directive with no table at all.
    let headerless = expect_refusal("#power-unit: W\n");
    assert_eq!(headerless.field, "body");

    // A header with no rows is not an empty-but-valid table.
    let no_rows = expect_refusal("#power-unit: W\ncomponent,power\n");
    assert_eq!(no_rows.field, "body");
    assert!(
        no_rows.reason.contains("no component rows"),
        "{}",
        no_rows.reason
    );
}

/// A duplicate component name refuses and names BOTH occurrences.
#[test]
fn a_duplicate_component_refuses_naming_the_first_occurrence() {
    let error = expect_refusal("#power-unit: W\ncomponent,power\ndie,2.5\ndie,1.0\n");
    assert_eq!((error.line, error.row), (4, Some(2)));
    assert_eq!(error.field, "component");
    assert!(error.reason.contains("die"), "{}", error.reason);
    assert!(
        error.reason.contains("line 3"),
        "the first occurrence is named so the fix is obvious: {}",
        error.reason
    );
}

/// An interval without a confidence states no uncertainty, so it refuses
/// rather than being silently promoted to a bare half width.
#[test]
fn an_interval_without_a_confidence_refuses() {
    let error =
        expect_refusal("#power-unit: W\ncomponent,power,half_width,confidence\ndie,2.5,0.1,\n");
    assert_eq!((error.line, error.row), (3, Some(1)));
    assert_eq!(error.field, "confidence");
    assert!(
        error.reason.contains("states no uncertainty"),
        "{}",
        error.reason
    );

    // A confidence outside (0, 1) is refused by the column range.
    let out_of_range =
        expect_refusal("#power-unit: W\ncomponent,power,half_width,confidence\ndie,2.5,0.1,1.5\n");
    assert_eq!(out_of_range.line, 3);
    assert!(
        out_of_range.reason.contains("1.5"),
        "{}",
        out_of_range.reason
    );
}

/// ADVERSARIAL: a declared total that disagrees with its rows is carried
/// through VERBATIM. The parser must not correct it, because the authoritative
/// balance refusal lives in `PowerMap` and two implementations would drift.
#[test]
fn a_disagreeing_declared_total_is_carried_through_not_corrected() {
    let table = parse_power_table(
        "#power-unit: W\n#declared-total: 100\ncomponent,power\ndie,2.5\nregulator,1.5\n",
    )
    .expect("a disagreeing total is not a PARSE fault");

    assert_eq!(table.declared_total_w(), Some(100.0));
    assert!((table.row_total_w() - 4.0).abs() < 1e-12);
    assert_ne!(
        table.declared_total_w(),
        Some(table.row_total_w()),
        "the parser must not reconcile the two; PowerMap owns that refusal"
    );
    // The receipt preserves the declared total so the downstream refusal can
    // quote what the file actually claimed.
    assert_eq!(table.receipt().declared_total_w(), Some(100.0));

    // A malformed total is a parse fault, though.
    let bad_total = expect_refusal("#power-unit: W\n#declared-total: lots\ncomponent,power\nd,1\n");
    assert_eq!((bad_total.field, bad_total.line), ("declared-total", 2));
    assert!(bad_total.reason.contains("lots"), "{}", bad_total.reason);

    let negative_total =
        expect_refusal("#power-unit: W\n#declared-total: -5\ncomponent,power\nd,1\n");
    assert_eq!(negative_total.field, "declared-total");
    assert!(
        negative_total.reason.contains("-5"),
        "{}",
        negative_total.reason
    );
}

/// The declared total converts with the same factor as the rows, so a mW file
/// stays internally consistent.
#[test]
fn the_declared_total_converts_with_its_rows() {
    let table = parse_power_table(
        "#power-unit: mW\n#declared-total: 4000\ncomponent,power\na,1000\nb,3000\n",
    )
    .expect("mW table");
    assert!((table.declared_total_w().expect("total") - 4.0).abs() < 1e-12);
    assert!((table.row_total_w() - 4.0).abs() < 1e-12);
}

/// Ingestion is deterministic and binds the exact bytes.
#[test]
fn ingestion_is_deterministic_and_byte_bound() {
    let first = parse_power_table(GOOD).expect("parse");
    let second = parse_power_table(GOOD).expect("parse");
    assert_eq!(first, second);
    assert_eq!(first.receipt().to_json(), second.receipt().to_json());

    // A single changed digit is a different source, and the receipt says so.
    let altered = parse_power_table(&GOOD.replace("2.5", "2.6")).expect("parse");
    assert_ne!(
        first.receipt().source_hash(),
        altered.receipt().source_hash(),
        "the receipt binds the exact bytes, not the parsed values"
    );
}
