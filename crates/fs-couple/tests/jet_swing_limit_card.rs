//! The card rows `jet_swing_limit fletcher-1976` and `jet_cut_width station` (MB 1784, 1816), parsed
//! on the committed D4 flute card at the coarse forced-jet row. The limit needs a named
//! `jet_spreading` law, because it takes the spread width at the lip. The width needs `jet_response`.
//! Without the rows the voice has neither.

use fs_couple::flue_loop::jet_swing_limit::{LAW, NAME, OURS, PRINTED};
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_exec::CancelGate;

const CARD: &str =
    include_str!("../../../data/instruments/descriptions/flute/flute-d4.performance");
const JET: &str = "jet 90 1e-6 3 0\n";
const RESPONSE: &str = "jet_response forced-bickley 0.25 40\n";
const SPREAD: &str = "jet_spreading fletcher-1976\n";

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
fn without_the_rows_the_voice_has_neither() {
    let card = parse_with(&format!("{RESPONSE}{SPREAD}")).expect("parses");
    let voice = card.renderer().voice();
    assert!(voice.jet_swing_limit().is_none());
    assert!(!voice.jet_cut_width_station());
    let (beta, _) = voice.forced_jet().unwrap().station_spread().unwrap();
    assert_eq!(
        voice.cut_half_width_m(),
        0.5 * card.geometry().flue_height_m * beta
    );
}

#[test]
fn the_limit_takes_the_stations_width_and_the_cut_up() {
    let card = parse_with(&format!(
        "{RESPONSE}{SPREAD}jet_swing_limit fletcher-1976\n"
    ))
    .expect("parses");
    let voice = card.renderer().voice();
    let kernel = voice.forced_jet().unwrap();
    let (beta, _) = kernel.station_spread().unwrap();
    let limit = voice.jet_swing_limit().expect("limit");
    assert_eq!(limit.spread_width_m(), beta * kernel.half_width_m());
    assert_eq!(limit.cut_up_m(), card.geometry().cut_up_m);
    assert_eq!(limit.gain(), 1.0);
    // The limit moves no scale of the cut.
    assert_eq!(
        voice.cut_half_width_m(),
        0.5 * card.geometry().flue_height_m * beta
    );
}

#[test]
fn the_cut_width_is_the_jets_own_with_or_without_a_spreading_law() {
    for (tail, spread) in [
        (RESPONSE.to_owned(), false),
        (format!("{RESPONSE}{SPREAD}"), true),
    ] {
        let card = parse_with(&format!("{tail}jet_cut_width station\n")).expect("parses");
        let voice = card.renderer().voice();
        let kernel = voice.forced_jet().unwrap();
        let beta = kernel.station_spread().map_or(1.0, |(beta, _)| beta);
        assert_eq!(kernel.station_spread().is_some(), spread);
        assert!(voice.jet_cut_width_station());
        assert_eq!(voice.cut_half_width_m(), kernel.half_width_m() * beta);
    }
}

#[test]
fn both_rows_attach_together() {
    let card = parse_with(&format!(
        "{RESPONSE}{SPREAD}jet_swing_limit fletcher-1976\njet_cut_width station\n"
    ))
    .expect("parses");
    let voice = card.renderer().voice();
    assert!(voice.jet_swing_limit().is_some());
    assert!(voice.jet_cut_width_station());
}

/// The row and the page must say the print gives an envelope and the time form is ours (MB 1784).
#[test]
fn the_row_says_what_is_printed_and_what_is_ours() {
    assert_eq!(NAME, "fletcher-1976");
    assert!(LAW.contains("x_t = X acosh(W w/v) / acosh(y w/v)"), "{LAW}");
    assert!(PRINTED.contains("his model choice"), "{PRINTED}");
    assert!(PRINTED.contains("has very little effect"), "{PRINTED}");
    assert!(PRINTED.contains("de la Cuadra and Fabre 2004"), "{PRINTED}");
    assert!(PRINTED.contains("Verge 1995"), "{PRINTED}");
    assert!(
        OURS.starts_with("the time form is ours, not printed"),
        "{OURS}"
    );
    let card = parse_with(&format!(
        "{RESPONSE}{SPREAD}jet_swing_limit fletcher-1976\n"
    ))
    .expect("parses");
    let provenance = &card.renderer().voice().island().provenance;
    assert!(provenance.contains(PRINTED), "{provenance}");
    assert!(provenance.contains(OURS), "{provenance}");
}

#[test]
fn the_rows_refuse_what_they_cannot_read() {
    for (tail, text) in [
        (
            format!("{RESPONSE}jet_swing_limit fletcher-1976\n"),
            "jet_swing_limit requires the jet_spreading row before it",
        ),
        (
            "jet_swing_limit fletcher-1976\n".to_owned(),
            "jet_swing_limit requires the jet_spreading row before it",
        ),
        (
            format!("{RESPONSE}{SPREAD}jet_swing_limit\n"),
            "jet_swing_limit requires fletcher-1976",
        ),
        (
            format!("{RESPONSE}{SPREAD}jet_swing_limit 15\n"),
            "jet_swing_limit requires fletcher-1976",
        ),
        (
            format!("{RESPONSE}{SPREAD}jet_swing_limit fletcher-1976 15\n"),
            "unexpected trailing fields",
        ),
        (
            "jet_cut_width station\n".to_owned(),
            "jet_cut_width requires the jet_response row before it",
        ),
        (
            format!("{RESPONSE}jet_cut_width\n"),
            "jet_cut_width requires station",
        ),
        (
            format!("{RESPONSE}jet_cut_width half-height\n"),
            "jet_cut_width requires station",
        ),
        (
            format!("{RESPONSE}jet_cut_width station 1.5\n"),
            "unexpected trailing fields",
        ),
        (
            format!("{RESPONSE}{SPREAD}jet_cut_width station\njet_swing_limit fletcher-1976\n"),
            "jet_swing_limit goes before jet_cut_width",
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
