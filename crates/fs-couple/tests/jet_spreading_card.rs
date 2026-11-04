//! The card row `jet_spreading LAW` (MB 97), parsed on the committed D4 flute card at the coarse
//! forced-jet row. LAW is required and only `bickley` (the laminar Bickley jet matched at the exit)
//! and `fletcher-1976` (an organ-flue Pitot measurement carried to the card's jet) are accepted.
//! The row needs `jet_response` before it, and goes after `momentum_drive` when both are present.
//! A card without the row has the parallel base and the unspread cut.

use fs_couple::flue_loop::forced_jet::{BAND_HALF_HEIGHT, JetSpreadingLaw};
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_exec::CancelGate;
use fs_math::det;

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

#[test]
fn both_laws_attach_and_a_card_without_the_row_is_unspread() {
    let none = parse_with(RESPONSE).expect("parses");
    let voice = none.renderer().voice();
    let kernel = voice.forced_jet().expect("forced jet");
    assert_eq!(kernel.spreading(), JetSpreadingLaw::None);
    assert!(kernel.station_spread().is_none());
    assert_eq!(
        kernel.step_response().resolution().half_height,
        BAND_HALF_HEIGHT
    );
    let unspread = voice.cut_half_width_m();
    for (row, law) in [
        ("jet_spreading bickley", JetSpreadingLaw::Bickley),
        ("jet_spreading fletcher-1976", JetSpreadingLaw::Fletcher1976),
    ] {
        let card = parse_with(&format!("{RESPONSE}{row}\n"))
            .unwrap_or_else(|e| panic!("'{row}' must parse: {e}"));
        let voice = card.renderer().voice();
        let kernel = voice.forced_jet().expect("forced jet");
        assert_eq!(kernel.spreading(), law);
        let (beta, flux) = kernel.station_spread().expect("spread");
        assert!(beta > 1.0, "{row}: beta {beta}");
        assert_eq!(flux, det::sqrt(beta));
        assert_eq!(voice.cut_half_width_m(), unspread * beta);
        println!(
            "{row}: station {:.4} b0, b(X)/b0 {beta:.4}, u_c(X)/U0 {:.4}, band {} b0",
            kernel.step_response().station(),
            1.0 / flux,
            kernel.step_response().resolution().half_height
        );
    }
}

/// The voice provenance and the render sidecar print `source()`; Fletcher's law must say it is
/// an organ-flue measurement carried over to the card's jet, not measured on it (MB 125, MB 1349).
#[test]
fn each_law_names_where_it_comes_from() {
    let fletcher = JetSpreadingLaw::Fletcher1976.source();
    assert!(
        fletcher.starts_with(
            "an organ-flue Pitot measurement (E4 pipe, flue 0.25-1.0 mm) carried to this jet, not measured on it: "
        ),
        "{fletcher}"
    );
    assert!(fletcher.contains("Fletcher 1976 eqs. 4-6, phi 5 deg (Table I)"));
    assert!(JetSpreadingLaw::Bickley.source().starts_with("the laminar Bickley jet"));
}

#[test]
fn the_row_follows_the_momentum_drive_and_needs_the_forced_jet() {
    let both = parse_with(&format!(
        "{RESPONSE}momentum_drive 1\njet_spreading fletcher-1976\n"
    ))
    .expect("momentum then spreading parses");
    let voice = both.renderer().voice();
    assert!(voice.momentum_drive().is_some());
    assert_eq!(
        voice.forced_jet().unwrap().spreading(),
        JetSpreadingLaw::Fletcher1976
    );
    for (tail, text) in [
        (
            format!("{RESPONSE}jet_spreading fletcher-1976\nmomentum_drive 1\n"),
            "momentum_drive goes before jet_spreading",
        ),
        (
            "jet_spreading bickley\n".to_owned(),
            "jet_spreading requires the jet_response row before it",
        ),
        (
            format!("{RESPONSE}jet_spreading\n"),
            "jet_spreading requires bickley or fletcher-1976",
        ),
        (
            format!("{RESPONSE}jet_spreading tan-phi\n"),
            "jet_spreading requires bickley or fletcher-1976",
        ),
        (
            format!("{RESPONSE}jet_spreading bickley 2\n"),
            "unexpected trailing fields",
        ),
    ] {
        match parse_with(&tail) {
            Ok(_) => panic!("accepted:\n{tail}"),
            Err(got) => assert!(
                got.contains(text),
                "{tail}: refused as `{got}`, expected `{text}`"
            ),
        }
    }
}
