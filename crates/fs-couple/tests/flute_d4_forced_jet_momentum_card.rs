//! The committed card `flute-d4-forced-jet-momentum-jet-noise-re700.performance` is the card the posted flute D4 (6/14, deploy 194) was
//! rendered from, byte for byte (md5 0f672e973794c8d3c88ffad8b22605c5, 5974 bytes). It is the measured-onset card
//! `flute-d4-onset-noise-contraction.performance` (md5 38b1bb2a: the committed flute-d4.performance plus `mouth_loss contraction`,
//! `jet_noise` at RE_ON 1700 and the 12 thesis Fig. 3.4 onset segments) plus two added rows, `jet_response forced-bickley 0.1 90` and
//! `momentum_drive 1` (Thwaites and Fletcher 1983, Eq. (1), as printed), and RE_ON 700 in place of 1700 (Verge 1995 p. 196's simulation
//! recipe onset; his 390 Pa is the measurement floor the 1700 came from). The card format has no comment grammar (a `#` line is refused),
//! so its source line is this header and the ledger row `delacuadra-2008-flute-controls-figure-reads`.
//!
//! The test parses the card with the reader and checks what the posted render depends on: the order of the three jet rows, the forced jet
//! with the momentum drive at weight 1 on the inlet section's area, the contraction row, the `jet_noise` record with onset Reynolds
//! number 700, and the 12 onset segments of the blowing-pressure schedule up to the 332 Pa hold and the release.

use fs_couple::flue_loop::momentum_drive::EQ1_KAPPA;
use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_exec::CancelGate;
use fs_scenario::gesture::GestureValue;

const CARD: &str = include_str!(
    "../../../data/instruments/descriptions/flute/flute-d4-forced-jet-momentum-jet-noise-re700.performance"
);

fn rows(prefix: &str) -> Vec<&'static str> {
    CARD.lines()
        .filter(|line| line.starts_with(prefix))
        .collect()
}

#[test]
fn the_posted_d4_card_has_the_forced_jet_the_momentum_drive_the_onset_and_the_gate_at_re_700() {
    assert_eq!(
        CARD.len(),
        5974,
        "the card is the posted bytes; do not edit it, add a card"
    );
    assert_eq!(rows("flue "), ["flue 0.00618 0.00113 0.00750 0.00029 0 0"]);
    assert_eq!(
        rows("jet_response "),
        ["jet_response forced-bickley 0.1 90"]
    );
    assert_eq!(rows("momentum_drive "), ["momentum_drive 1"]);
    assert_eq!(rows("jet_noise "), ["jet_noise 14 1.3e-3 0.045 0.115 700"]);
    assert_eq!(rows("mouth_loss ").len(), 1);
    let lines: Vec<&str> = CARD.lines().collect();
    let jet = lines
        .iter()
        .position(|line| *line == "jet 90 1e-6 3 0")
        .expect("the jet row");
    assert_eq!(
        &lines[jet + 1..jet + 3],
        ["jet_response forced-bickley 0.1 90", "momentum_drive 1"],
        "the momentum row follows jet_response directly"
    );

    let card = FlueFilePerformance::from_bytes(CARD.as_bytes(), 256, &CancelGate::new_clock_free())
        .expect("the posted D4 card parses");
    let voice = card.renderer().voice();
    assert!(
        voice.forced_jet().is_some(),
        "the forced Bickley jet is attached"
    );
    let drive = voice
        .momentum_drive()
        .expect("the momentum drive is attached");
    assert_eq!(drive.kappa().to_bits(), EQ1_KAPPA.to_bits());
    assert_eq!(
        drive.bore_area_m2().to_bits(),
        card.geometry().bore_area_m2.to_bits()
    );
    assert!(
        voice.mouth_contraction().is_some(),
        "the contraction row is attached"
    );
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
    assert_eq!(
        pressures,
        [
            60.8, 109.5, 101.4, 113.5, 141.9, 174.3, 198.6, 202.7, 218.9, 247.3, 291.9, 332.0, 0.0
        ],
        "12 onset segments to the 332 Pa hold, then the release"
    );
}
