//! The card row `jet_forcing_area open-strip` (MB 1089), parsed on the committed D4 flute card at
//! the coarse forced-jet row. The row follows a `window embouchure` row and needs `jet_response`.
//! A card without it forces the jet over flue width by cut-up.

use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_exec::CancelGate;

const CARD: &str =
    include_str!("../../../data/instruments/descriptions/flute/flute-d4.performance");
const JET: &str = "jet 90 1e-6 3 0\n";
const RESPONSE: &str = "jet_response forced-bickley 0.25 40\n";
const WINDOW: &str =
    "window embouchure 16 0.0098740 0.0118488 0.0058387 lip one-sided face 0.0018\n";
const ROW: &str = "jet_forcing_area open-strip\n";
const ROW_SOURCE: &str = "jet_forcing_area open-strip source\n";

fn parse(response: &str, window: &str) -> Result<FlueFilePerformance, String> {
    assert_eq!(CARD.matches(JET).count(), 1, "the D4 jet row changed; update this test");
    assert_eq!(CARD.matches(WINDOW).count(), 1, "the D4 window row changed; update this test");
    let text = CARD
        .replace(JET, &format!("{JET}{response}"))
        .replace(WINDOW, window);
    FlueFilePerformance::from_bytes(text.as_bytes(), 256, &CancelGate::new_clock_free())
        .map_err(|e| format!("{e:?}"))
}

#[test]
fn the_row_names_the_open_strip_and_a_card_without_it_keeps_the_flue_rows_area() {
    let without = parse(RESPONSE, WINDOW).expect("parses");
    assert_eq!(without.jet_forcing_area_m2(), None);
    assert_eq!(
        without.renderer().voice().jet_forcing_area_m2(),
        without.geometry().mouth_area_m2()
    );

    let with = parse(RESPONSE, &format!("{WINDOW}{ROW}")).expect("parses");
    let strip = 0.0118488 * 0.0058387;
    assert_eq!(with.jet_forcing_area_m2(), Some(strip));
    assert_eq!(with.renderer().voice().jet_forcing_area_m2(), strip);
    // Only the jet's velocity moves: the series mass stays where the window put it.
    assert_eq!(
        with.renderer().voice().mouth_inertance_pa_s2_m3(),
        without.renderer().voice().mouth_inertance_pa_s2_m3()
    );
}

#[test]
fn the_source_term_moves_the_source_coefficient_to_the_strip_and_the_bare_row_does_not() {
    let without = parse(RESPONSE, WINDOW).expect("parses");
    let bare = parse(RESPONSE, &format!("{WINDOW}{ROW}")).expect("parses");
    let with = parse(RESPONSE, &format!("{WINDOW}{ROW_SOURCE}")).expect("parses");
    assert!(!without.jet_source_on_strip());
    assert!(!bare.jet_source_on_strip());
    assert!(with.jet_source_on_strip());

    let plain = without.renderer().voice().jet_source_coefficient();
    assert_eq!(bare.renderer().voice().jet_source_coefficient().to_bits(), plain.to_bits());

    let strip = 0.0118488 * 0.0058387;
    let s_m = without.geometry().mouth_area_m2();
    let moved = with.renderer().voice().jet_source_coefficient();
    assert!((moved / plain - (s_m / strip).sqrt()).abs() < 1e-14);
    // The velocity's area stays the strip, and the series mass stays where the window put it.
    assert_eq!(with.jet_forcing_area_m2(), Some(strip));
    assert_eq!(with.renderer().voice().jet_forcing_area_m2(), strip);
    assert_eq!(
        with.renderer().voice().mouth_inertance_pa_s2_m3(),
        without.renderer().voice().mouth_inertance_pa_s2_m3()
    );
}

#[test]
fn the_row_refuses_what_it_cannot_read() {
    for (response, window, text) in [
        (RESPONSE, format!("{WINDOW}jet_forcing_area\n"), "jet_forcing_area requires open-strip"),
        (RESPONSE, format!("{WINDOW}jet_forcing_area riser\n"), "jet_forcing_area requires open-strip"),
        (RESPONSE, format!("{WINDOW}jet_forcing_area open-strip 2\n"), "takes only the term source"),
        (RESPONSE, format!("{WINDOW}jet_forcing_area open-strip contraction\n"), "takes only the term source"),
        (RESPONSE, format!("{WINDOW}jet_forcing_area open-strip source 2\n"), "unexpected trailing fields"),
        ("", format!("{WINDOW}{ROW_SOURCE}"), "jet forcing area requires a cold voice with the forced jet"),
        (
            RESPONSE,
            ROW.to_owned(),
            "jet_forcing_area open-strip requires a window embouchure row before it",
        ),
        (
            RESPONSE,
            format!("window wall-cut 16\n{ROW}"),
            "jet_forcing_area open-strip requires a window embouchure row before it",
        ),
        ("", format!("{WINDOW}{ROW}"), "jet forcing area requires a cold voice with the forced jet"),
    ] {
        match parse(response, &window) {
            Ok(_) => panic!("accepted:\n{response}{window}"),
            Err(got) => assert!(got.contains(text), "{window}: refused as `{got}`, expected `{text}`"),
        }
    }
}
