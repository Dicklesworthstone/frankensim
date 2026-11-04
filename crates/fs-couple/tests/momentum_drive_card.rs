//! The card row `momentum_drive KAPPA` (MB 16590), parsed on the committed D4 flute card at the coarse
//! forced-jet row. KAPPA is required and only the two printed weights are accepted: `1`, Thwaites and
//! Fletcher (1983, J. Acoust. Soc. Am. 74, 400) Eq. (1) p. 401, and `2.6`, their Sec. II fit p. 402.
//! The row needs `jet_response` before it. A card without the row has no momentum drive.

use fs_couple::flue_loop::momentum_drive::{EQ1_KAPPA, KAPPA_REFUSED, SEC2_KAPPA};
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_exec::CancelGate;

const CARD: &str =
    include_str!("../../../data/instruments/descriptions/flute/flute-d4.performance");
const JET: &str = "jet 90 1e-6 3 0\n";
const RESPONSE: &str = "jet_response forced-bickley 0.25 40\n";

fn parse_with(tail: &str) -> Result<FlueFilePerformance, String> {
    assert_eq!(
        CARD.matches(JET).count(),
        1,
        "the D4 jet row changed; update this test"
    );
    let text = CARD.replace(JET, &format!("{JET}{tail}"));
    FlueFilePerformance::from_bytes(text.as_bytes(), 256, &CancelGate::new_clock_free())
        .map_err(|e| format!("{e:?}"))
}

fn parse(row: &str) -> Result<FlueFilePerformance, String> {
    parse_with(&format!("{RESPONSE}{row}\n"))
}

#[test]
fn both_printed_weights_attach_and_a_card_without_the_row_has_none() {
    let none = parse_with(RESPONSE).expect("parses");
    assert!(none.renderer().voice().forced_jet().is_some());
    assert!(none.renderer().voice().momentum_drive().is_none());
    for (row, kappa) in [
        ("momentum_drive 1", EQ1_KAPPA),
        ("momentum_drive 1.0", EQ1_KAPPA),
        ("momentum_drive 1e0", EQ1_KAPPA),
        ("momentum_drive 2.6", SEC2_KAPPA),
        ("momentum_drive 2.60", SEC2_KAPPA),
    ] {
        let card = parse(row).unwrap_or_else(|e| panic!("'{row}' must parse: {e}"));
        let voice = card.renderer().voice();
        assert!(voice.forced_jet().is_some());
        let drive = voice.momentum_drive().expect("attached");
        assert_eq!(drive.kappa().to_bits(), kappa.to_bits(), "{row}");
        assert_eq!(
            drive.bore_area_m2().to_bits(),
            card.geometry().bore_area_m2.to_bits()
        );
    }
    assert_eq!((EQ1_KAPPA, SEC2_KAPPA), (1.0, 2.6));
}

#[test]
fn kappa_is_required_and_any_other_value_is_refused_naming_both_printed_values() {
    for row in [
        "momentum_drive",
        "momentum_drive 0",
        "momentum_drive -1",
        "momentum_drive 0.5",
        "momentum_drive 2",
        "momentum_drive 3",
        "momentum_drive 4",
        "momentum_drive 2.6000001",
        "momentum_drive 1.0000001",
        "momentum_drive nan",
        "momentum_drive inf",
        "momentum_drive one",
    ] {
        let err = parse(row)
            .err()
            .unwrap_or_else(|| panic!("'{row}' must be refused"));
        assert!(err.contains(KAPPA_REFUSED), "'{row}' refused with: {err}");
        println!("'{row}' refused: {err}");
    }
    assert!(KAPPA_REFUSED.contains("1 (") && KAPPA_REFUSED.contains("2.6 ("));
}

#[test]
fn malformed_rows_and_a_wrong_place_are_refused() {
    // Trailing fields.
    assert!(parse("momentum_drive 1 extra").is_err());
    assert!(parse("momentum_drive 1 2.6").is_err());
    // It acts on the forced jet: no `jet_response` row, no momentum drive.
    let err = parse_with("momentum_drive 1\n").err().expect("refused");
    assert!(
        err.contains("requires the jet_response row before it"),
        "{err}"
    );
    // Before `jet_response` is the wrong place.
    let err = parse_with(&format!("momentum_drive 1\n{RESPONSE}"))
        .err()
        .expect("refused");
    assert!(
        err.contains("requires the jet_response row before it"),
        "{err}"
    );
    // Twice.
    assert!(parse("momentum_drive 1\nmomentum_drive 1").is_err());
}

#[test]
fn the_row_changes_the_loop_not_the_card() {
    let a = parse("momentum_drive 1").expect("parses");
    let b = parse("momentum_drive 2.6").expect("parses");
    let none = parse_with(RESPONSE).expect("parses");
    // Same card, same forced-jet kernel: only the added source differs.
    for card in [&a, &b] {
        let (x, y) = (
            card.renderer().voice().forced_jet().unwrap(),
            none.renderer().voice().forced_jet().unwrap(),
        );
        assert_eq!(
            x.step_response().values().len(),
            y.step_response().values().len()
        );
        for (p, q) in x
            .step_response()
            .values()
            .iter()
            .zip(y.step_response().values())
        {
            assert_eq!(p.to_bits(), q.to_bits());
        }
    }
    assert!(
        b.renderer().voice().momentum_drive().unwrap().kappa()
            > a.renderer().voice().momentum_drive().unwrap().kappa()
    );
}
