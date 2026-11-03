//! The committed card `flute-g4-extrap-flue-jet-noise-re700.performance` is the card the posted flute G4 was rendered from, byte for
//! byte (md5 3b448d0d7b69a5aac36c296888cc6a9e, 5346 bytes). The card format has no comment grammar (every record is read in order,
//! and a `#` line is refused: `flue input line 2: clock`), so the card's provenance is stated here and in ledger row
//! `delacuadra-2008-flute-controls`:
//!
//! Posted G4 (4/14, deploy 194). Origin: flute-g4.performance as committed at 2b6fcc760 (estimated flue: width 5.94 mm, height
//! 0.90 mm, cut-up 7.08 mm; 379 Pa hold, Re_h 1494), before a7afc67bd replaced it with the measured flue (width 5.89, height 0.98,
//! cut-up 6.25 mm; graded 3/14), plus ONE added line, the jet_noise row: Verge 1995's measured flue-jet noise fraction
//! (1.3e-3 rho U_j^2/2 over St_h 0.045-0.115, slit 1 mm), carried to this flue as an extrapolation. RE_ON 700 is Verge 1995
//! p. 196's simulation recipe onset; the first version of this card used 1490 (SunnySnow 46d937e3).
//!
//! The test parses the card with the reader and checks what the posted render depends on: the `flue` record in the reader's field
//! order (`flue_width_m, flue_height_m, cut_up_m, labium_offset_m, jet_angle_rad, mouth_loss_rel`), the 379 Pa blowing-pressure
//! hold, and the `jet_noise` record with its onset Reynolds number 700. It is a card with no other opt-in jet row.

use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_exec::CancelGate;
use fs_scenario::gesture::GestureValue;

const CARD: &str = include_str!(
    "../../../data/instruments/descriptions/flute/flute-g4-extrap-flue-jet-noise-re700.performance"
);

fn rows(prefix: &str) -> Vec<&'static str> {
    CARD.lines()
        .filter(|line| line.starts_with(prefix))
        .collect()
}

#[test]
fn the_posted_g4_card_has_its_estimated_flue_its_379_pa_hold_and_the_jet_noise_gate_at_re_700() {
    assert_eq!(
        CARD.len(),
        5346,
        "the card is the posted bytes; do not edit it, add a card"
    );
    assert_eq!(rows("flue "), ["flue 0.00594 0.00090 0.00708 0.00029 0 0"]);
    assert_eq!(rows("jet_noise "), ["jet_noise 14 1.3e-3 0.045 0.115 700"]);
    for absent in [
        "jet_response ",
        "jet_transfer ",
        "momentum_drive ",
        "mouth_loss ",
    ] {
        assert!(
            rows(absent).is_empty(),
            "the card has no `{absent}` row: {:?}",
            rows(absent)
        );
    }
    assert_eq!(rows("windway "), ["windway none"]);

    let card = FlueFilePerformance::from_bytes(CARD.as_bytes(), 256, &CancelGate::new_clock_free())
        .expect("the posted G4 card parses");
    let geometry = card.geometry();
    for (name, got, want) in [
        ("flue width", geometry.flue_width_m, 5.94e-3),
        ("flue height", geometry.flue_height_m, 0.90e-3),
        ("cut-up", geometry.cut_up_m, 7.08e-3),
        ("labium offset", geometry.labium_offset_m, 0.29e-3),
        ("jet angle", geometry.jet_angle_rad, 0.0),
        ("mouth loss", geometry.mouth_loss_rel, 0.0),
    ] {
        assert_eq!(got.to_bits(), f64::to_bits(want), "{name}");
    }

    let voice = card.renderer().voice();
    let noise = voice.jet_noise().expect("the jet_noise row attaches");
    assert_eq!(noise.seed, 14);
    for (name, got, want) in [
        ("gain", noise.gain, 1.3e-3),
        ("St_h low", noise.strouhal_lo, 0.045),
        ("St_h high", noise.strouhal_hi, 0.115),
        ("onset Re_h", noise.reynolds_onset, 700.0),
    ] {
        assert_eq!(got.to_bits(), f64::to_bits(want), "{name}");
    }
    assert!(
        voice.forced_jet().is_none(),
        "the card keeps the island jet"
    );

    let tracks = card.renderer().schedule().tracks();
    assert_eq!(tracks.len(), 1);
    assert_eq!(tracks[0].id, "blow");
    let pressures: Vec<f64> = tracks[0]
        .events
        .iter()
        .map(|event| match event.value {
            GestureValue::PressurePa(pa) => pa,
            ref other => panic!("a blowing-pressure event, not {other:?}"),
        })
        .collect();
    assert_eq!(pressures, [379.0, 0.0], "a 379 Pa hold, then the release");
    assert_eq!(
        tracks[0].events[0].transition_s.to_bits(),
        5e-2f64.to_bits()
    );
}
