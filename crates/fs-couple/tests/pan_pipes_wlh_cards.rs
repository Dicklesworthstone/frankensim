//! The four committed cards `pan-pipes-<note>-forced-jet-spreading-half-width-offset-re700.performance` (D3, G3, D4, B4) are the cards the posted pan-pipe set (deploy 243, the `wlh` rung)
//! was rendered from, byte for byte. The card format has no comment grammar (a `#` line is refused: `flue input line 2: clock`), so their provenance is stated here and in ledger row
//! `pan-pipes-wlh-posted-cards`:
//!
//! Each card is the committed `pan-pipes-<note>.performance` plus: the forced Bickley jet (`jet_response forced-bickley 0.1 90`) with the momentum drive at weight 1 (`momentum_drive 1`,
//! Thwaites and Fletcher 1983, Eq. (1), as printed) and `jet_spreading fletcher-1976` (a TRIAL law: Fletcher's flues were 0.25 to 1.0 mm high, these are 1.00 to 2.00 mm), `jet_noise 14 1.3e-3 0.045 0.115 700`
//! (Verge 1995's flue-jet noise fraction carried to these flues, onset Re_h 700 from his simulation recipe, p. 196), `mouth_loss contraction 0.6 0.3 Auvray2014-TableI-alpha_vc` (a recorder mouth), the computed
//! baffled radiation, the Stokes-Kirchhoff observer gas, and the labium offset, fourth `flue` field, set to HALF THE SPREAD CUT HALF-WIDTH (offset/b 0.50 on all four, the stopped-pipe offset Fletcher 2005
//! Eqs. (4) to (6) gives as about 0.5 b), with the flue's sixth field (mouth loss) 0 because the contraction row replaces it.
//!
//! The test parses each card with the reader and checks what the posted render depends on: the exact length, the `flue` record in the reader's field order (`flue_width_m, flue_height_m, cut_up_m,
//! labium_offset_m, jet_angle_rad, mouth_loss_rel`), the order of the jet rows, the forced jet with its spreading law, the momentum drive at weight 1 on the bore's area, the contraction row, the
//! `jet_noise` record with its onset Reynolds number 700, and the blowing-pressure hold with its 50 ms onset ramp and its release.
//!
//! Cost: the reader builds a renderer for every card it parses, which takes minutes in a release build (175.4 s for the four on ts2), so one test pins all four cards by their bytes and rows (the length, the
//! `flue` record, the jet rows, the contraction row and the two blowing-pressure event rows) and a second test parses ONE card, B4 (the shortest tube), and checks what the reader makes of it.

use fs_couple::flue_loop::forced_jet::JetSpreadingLaw;
use fs_couple::flue_loop::momentum_drive::EQ1_KAPPA;
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_exec::CancelGate;
use fs_scenario::gesture::GestureValue;

const D3: &str = include_str!(
    "../../../data/instruments/descriptions/pan-pipes/pan-pipes-d3-forced-jet-spreading-half-width-offset-re700.performance"
);
const G3: &str = include_str!(
    "../../../data/instruments/descriptions/pan-pipes/pan-pipes-g3-forced-jet-spreading-half-width-offset-re700.performance"
);
const D4: &str = include_str!(
    "../../../data/instruments/descriptions/pan-pipes/pan-pipes-d4-forced-jet-spreading-half-width-offset-re700.performance"
);
const B4: &str = include_str!(
    "../../../data/instruments/descriptions/pan-pipes/pan-pipes-b4-forced-jet-spreading-half-width-offset-re700.performance"
);

/// (note, card, posted byte length, `flue` row, the blowing-pressure hold as written in the card's first event row,
/// [width, height, cut-up, labium offset] in m)
const CARDS: [(&str, &str, usize, &str, &str, [f64; 4]); 4] = [
    (
        "D3",
        D3,
        2359,
        "flue 0.02033 0.00200 0.010053 0.0010863 0 0",
        "1.1747e2",
        [0.02033, 0.00200, 0.010053, 0.0010863],
    ),
    (
        "G3",
        G3,
        2357,
        "flue 0.01794 0.00176 0.008943 0.0009616 0 0",
        "1.568e2",
        [0.01794, 0.00176, 0.008943, 0.0009616],
    ),
    (
        "D4",
        D4,
        2358,
        "flue 0.01506 0.00143 0.007635 0.0008028 0 0",
        "2.3493e2",
        [0.01506, 0.00143, 0.007635, 0.0008028],
    ),
    (
        "B4",
        B4,
        2358,
        "flue 0.01201 0.00100 0.006300 0.0006175 0 0",
        "3.9511e2",
        [0.01201, 0.00100, 0.006300, 0.0006175],
    ),
];

/// The one card the reader is run on (a renderer per card costs minutes).
const PARSED: &str = "B4";

fn rows(card: &'static str, prefix: &str) -> Vec<&'static str> {
    card.lines()
        .filter(|line| line.starts_with(prefix))
        .collect()
}

#[test]
fn the_posted_pan_pipe_cards_keep_their_bytes_and_rows() {
    for (note, card_text, bytes, flue_row, hold, _) in CARDS {
        assert_eq!(
            card_text.len(),
            bytes,
            "{note}: the card is the posted bytes; do not edit it, add a card"
        );
        assert_eq!(rows(card_text, "flue "), [flue_row], "{note}");
        assert_eq!(
            rows(card_text, "jet_response "),
            ["jet_response forced-bickley 0.1 90"],
            "{note}"
        );
        assert_eq!(
            rows(card_text, "momentum_drive "),
            ["momentum_drive 1"],
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
        assert!(
            rows(card_text, "jet_transfer ").is_empty(),
            "{note}: the card has no `jet_transfer` row"
        );
        assert_eq!(
            rows(card_text, "event\t"),
            [
                format!("event\t0e0\t5e-2\tpressure-pa\t{hold}"),
                "event\t2.2e0\t1.05e-2\tpressure-pa\t0e0".to_string()
            ],
            "{note}: the 50 ms onset ramp to the hold, then the 10.5 ms release"
        );
    }
}

#[test]
fn the_reader_makes_the_forced_jet_the_momentum_drive_and_the_noise_of_one_posted_card() {
    let (note, card_text, _, _, hold, flue) = CARDS
        .into_iter()
        .find(|(note, ..)| *note == PARSED)
        .expect("the parsed card is one of the four");
    let card =
        FlueFilePerformance::from_bytes(card_text.as_bytes(), 256, &CancelGate::new_clock_free())
            .unwrap_or_else(|e| panic!("{note}: the posted card parses: {e:?}"));
    let geometry = card.geometry();
    for (name, got, want) in [
        ("flue width", geometry.flue_width_m, flue[0]),
        ("flue height", geometry.flue_height_m, flue[1]),
        ("cut-up", geometry.cut_up_m, flue[2]),
        ("labium offset", geometry.labium_offset_m, flue[3]),
        ("jet angle", geometry.jet_angle_rad, 0.0),
        ("mouth loss", geometry.mouth_loss_rel, 0.0),
    ] {
        assert_eq!(got.to_bits(), f64::to_bits(want), "{note}: {name}");
    }

    let voice = card.renderer().voice();
    let kernel = voice.forced_jet().expect("the forced jet attaches");
    assert_eq!(kernel.spreading(), JetSpreadingLaw::Fletcher1976, "{note}");
    let drive = voice
        .momentum_drive()
        .expect("the momentum drive is attached");
    assert_eq!(drive.kappa().to_bits(), EQ1_KAPPA.to_bits(), "{note}");
    assert_eq!(
        drive.bore_area_m2().to_bits(),
        card.geometry().bore_area_m2.to_bits(),
        "{note}"
    );
    assert!(
        voice.mouth_contraction().is_some(),
        "{note}: the contraction row is attached"
    );
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
    let want_hold: f64 = hold.parse().expect("the hold is a number");
    assert_eq!(
        pressures,
        [want_hold, 0.0],
        "{note}: the hold, then the release"
    );
    assert_eq!(
        tracks[0].events[0].transition_s.to_bits(),
        5.0e-2f64.to_bits(),
        "{note}: the 50 ms onset ramp"
    );
    assert_eq!(
        tracks[0].events[1].transition_s.to_bits(),
        1.05e-2f64.to_bits(),
        "{note}: the 10.5 ms release"
    );
}
