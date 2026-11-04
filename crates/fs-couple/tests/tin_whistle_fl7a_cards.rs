//! The four committed cards `tin-whistle-<note>-fletcher-spreading-borrowed-offset-re700.performance` are the cards the posted whistle set (D5, G5, B5, D6; deploy 234) was rendered
//! from, byte for byte. The card format has no comment grammar (a `#` line is refused: `flue input line 2: clock`), so their provenance is stated here and in ledger row
//! `tin-whistle-fl7a-posted-cards`:
//!
//! Each card is the committed `tin-whistle-<note>.performance` plus `mouth_loss contraction 0.6 0.3 Auvray2014-TableI-alpha_vc`, `observer_gas`, `inner 0.0061000 matching series` on the
//! open chimneys (G5 holes 4 to 6, B5 holes 2 to 6, D6 hole 1), `jet_spreading fletcher-1976` (a TRIAL law: Fletcher's flue 0.25 to 1.0 mm, this one 1.48 mm), `jet_noise 14 1.3e-3 0.045 0.115 700`
//! (Verge 1995 extrapolation, onset Re_h 700) and the labium offset 0.000203 m in the fourth `flue` field. That offset is BORROWED from the recorder (0.1374 h = the recorder card's mid-face
//! ratio 0.1455 / 1.059 mm), not measured on any whistle; it was 0 (UNKNOWN).
//!
//! The test parses each card with the reader and checks what the posted render depends on: the exact length, the `flue` record in the reader's field order (`flue_width_m, flue_height_m,
//! cut_up_m, labium_offset_m, jet_angle_rad, mouth_loss_rel`), the `jet_spreading` law the forced jet attaches, the `jet_noise` record with its onset Reynolds number 700, the blowing-pressure
//! hold and its release, and that no `momentum_drive` or `jet_transfer` row is present.

use fs_couple::flue_loop::forced_jet::JetSpreadingLaw;
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_exec::CancelGate;
use fs_scenario::gesture::GestureValue;

const D5: &str = include_str!(
    "../../../data/instruments/descriptions/tin-whistle/tin-whistle-d5-fletcher-spreading-borrowed-offset-re700.performance"
);
const G5: &str = include_str!(
    "../../../data/instruments/descriptions/tin-whistle/tin-whistle-g5-fletcher-spreading-borrowed-offset-re700.performance"
);
const B5: &str = include_str!(
    "../../../data/instruments/descriptions/tin-whistle/tin-whistle-b5-fletcher-spreading-borrowed-offset-re700.performance"
);
const D6: &str = include_str!(
    "../../../data/instruments/descriptions/tin-whistle/tin-whistle-d6-fletcher-spreading-borrowed-offset-re700.performance"
);

/// (note, card, posted byte length, blowing-pressure hold in Pa)
const CARDS: [(&str, &str, usize, f64); 4] = [
    ("D5", D5, 2112, 93.447),
    ("G5", G5, 2280, 151.84),
    ("B5", B5, 2390, 198.32),
    ("D6", D6, 2167, 292.44),
];

fn rows(card: &'static str, prefix: &str) -> Vec<&'static str> {
    card.lines()
        .filter(|line| line.starts_with(prefix))
        .collect()
}

#[test]
fn the_posted_whistle_cards_keep_their_bytes_rows_and_borrowed_offset() {
    for (note, card_text, bytes, hold) in CARDS {
        assert_eq!(
            card_text.len(),
            bytes,
            "{note}: the card is the posted bytes; do not edit it, add a card"
        );
        assert_eq!(
            rows(card_text, "flue "),
            ["flue 0.00800 0.00148 0.00530 0.000203 0 0"],
            "{note}"
        );
        assert_eq!(
            rows(card_text, "jet_spreading "),
            ["jet_spreading fletcher-1976"],
            "{note}"
        );
        assert_eq!(
            rows(card_text, "jet_noise "),
            ["jet_noise 14 1.3e-3 0.045 0.115 700"],
            "{note}"
        );
        assert_eq!(
            rows(card_text, "mouth_loss "),
            ["mouth_loss contraction 0.6 0.3 Auvray2014-TableI-alpha_vc"],
            "{note}"
        );
        for absent in ["momentum_drive ", "jet_transfer "] {
            assert!(
                rows(card_text, absent).is_empty(),
                "{note}: the card has no `{absent}` row"
            );
        }

        let card = FlueFilePerformance::from_bytes(
            card_text.as_bytes(),
            256,
            &CancelGate::new_clock_free(),
        )
        .unwrap_or_else(|e| panic!("{note}: the posted card parses: {e:?}"));
        let geometry = card.geometry();
        for (name, got, want) in [
            ("flue width", geometry.flue_width_m, 8.0e-3),
            ("flue height", geometry.flue_height_m, 1.48e-3),
            ("cut-up", geometry.cut_up_m, 5.30e-3),
            ("labium offset", geometry.labium_offset_m, 0.203e-3),
            ("jet angle", geometry.jet_angle_rad, 0.0),
            ("mouth loss", geometry.mouth_loss_rel, 0.0),
        ] {
            assert_eq!(got.to_bits(), f64::to_bits(want), "{note}: {name}");
        }

        let voice = card.renderer().voice();
        let kernel = voice.forced_jet().expect("the forced jet attaches");
        assert_eq!(kernel.spreading(), JetSpreadingLaw::Fletcher1976, "{note}");
        let noise = voice.jet_noise().expect("the jet_noise row attaches");
        assert_eq!(noise.seed, 14, "{note}");
        for (name, got, want) in [
            ("gain", noise.gain, 1.3e-3),
            ("St_h low", noise.strouhal_lo, 0.045),
            ("St_h high", noise.strouhal_hi, 0.115),
            ("onset Re_h", noise.reynolds_onset, 700.0),
        ] {
            assert_eq!(got.to_bits(), f64::to_bits(want), "{note}: {name}");
        }

        let tracks = card.renderer().schedule().tracks();
        assert_eq!(tracks.len(), 1, "{note}");
        assert_eq!(tracks[0].id, "blow", "{note}");
        let pressures: Vec<f64> = tracks[0]
            .events
            .iter()
            .map(|event| match event.value {
                GestureValue::PressurePa(pa) => pa,
                ref other => panic!("{note}: a blowing-pressure event, not {other:?}"),
            })
            .collect();
        assert_eq!(pressures, [hold, 0.0], "{note}: the hold, then the release");
        assert_eq!(
            tracks[0].events[0].transition_s.to_bits(),
            1.25e-2f64.to_bits(),
            "{note}: the 12.5 ms onset ramp"
        );
    }
}
