//! The wired corpus gate (music bead `frankensim-music-v8-root-3ez8g.1.1`,
//! third DONE-WHEN): the Ernoult 2021 fingering ladder driven END-TO-END
//! from the registered fs-vvreg acoustic corpus rows — geometry, measured
//! peaks, and acceptance envelopes all come FROM the corpus, so this test
//! proves the loop data -> corpus row -> instrument model -> gate closes.
//!
//! This deliberately duplicates none of the inline `fs-duct` ladder test:
//! that one pins the crate's own conformance surface with local constants;
//! THIS one proves the corpus registration is consumable — the duct is
//! REBUILT from the corpus rows' context ranges, and the cents envelopes
//! are read from the rows' acceptance records. If the corpus geometry or
//! envelope drifts from what the model needs, this test fails while the
//! inline one still passes — exactly the seam it exists to watch.
//!
//! Per-fingering JSON-lines logging (suite/case/rows/verdict) mirrors the
//! repo's structured-test-log convention so a reviewer can re-derive the
//! verdict from the output alone.

use fs_duct::{Duct, HoleState, LossModel, Segment, Termination, impedance_peaks, impedance_sweep};
use fs_material::gas::{GasSpec, GasState};
use fs_vvreg::acoustic::{AcousticAcceptance, AcousticCase, acoustic_cases};

/// Pull a context value (a point range) out of a corpus row by axis name.
fn context_value(case: &AcousticCase, name: &str) -> f64 {
    let ctx = case
        .context
        .iter()
        .find(|ctx| ctx.name == name)
        .unwrap_or_else(|| panic!("{}: missing context axis {name}", case.id));
    assert!(
        (ctx.hi - ctx.lo).abs() < f64::EPSILON,
        "{}: {name} must be a point range to rebuild exact geometry",
        case.id
    );
    ctx.lo
}

/// Rebuild the four-hole cylinder FROM a corpus row's context geometry.
fn duct_from_corpus(case: &AcousticCase, states: [HoleState; 4]) -> Duct {
    let bore = context_value(case, "bore-radius-m");
    let total_length = context_value(case, "bore-length-m");
    let mut segments = Vec::new();
    let mut cursor = 0.0;
    for (index, state) in states.iter().enumerate() {
        let position = context_value(case, &format!("hole-{}-position-m", index + 1));
        let radius = context_value(case, &format!("hole-{}-radius-m", index + 1));
        let chimney = context_value(case, &format!("hole-{}-chimney-m", index + 1));
        segments.push(Segment::Cylinder {
            radius: bore,
            length: position - cursor,
        });
        segments.push(Segment::ToneHole {
            hole_radius: radius,
            chimney_height: chimney,
            bore_radius: bore,
            state: *state,
        });
        cursor = position;
    }
    segments.push(Segment::Cylinder {
        radius: bore,
        length: total_length - cursor,
    });
    Duct { segments }
}

fn first_peak_hz(duct: &Duct, state: &GasState) -> f64 {
    let sweep = impedance_sweep(
        duct,
        state,
        2.0 * core::f64::consts::PI * 150.0,
        2.0 * core::f64::consts::PI * 1000.0,
        12_000,
        LossModel::WideTube,
        Termination::UnflangedOpen,
    )
    .expect("sweep");
    let peaks = impedance_peaks(&sweep);
    assert!(!peaks.is_empty(), "no impedance peak in [150, 1000] Hz");
    sweep[peaks[0]].omega / (2.0 * core::f64::consts::PI)
}

#[test]
fn ernoult_ladder_driven_from_the_registered_corpus_rows() {
    use HoleState::{Closed as X, Open as O};
    use core::fmt::Write as _;
    let fingerings: [(&str, [HoleState; 4]); 5] = [
        ("acoustic-ernoult-2021-xxxx", [X, X, X, X]),
        ("acoustic-ernoult-2021-xxxo", [X, X, X, O]),
        ("acoustic-ernoult-2021-xxox", [X, X, O, X]),
        ("acoustic-ernoult-2021-xoxx", [X, O, X, X]),
        ("acoustic-ernoult-2021-oxxx", [O, X, X, X]),
    ];
    let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air");
    let mut rows = String::new();
    let mut previous = 0.0f64;
    for (index, (case_id, states)) in fingerings.iter().enumerate() {
        let case = acoustic_cases()
            .iter()
            .find(|case| case.id == *case_id)
            .unwrap_or_else(|| panic!("corpus row {case_id} must be registered"));
        // The acceptance envelope comes FROM the row, never from a local
        // constant — the corpus is the authority this gate consumes.
        let AcousticAcceptance::CentsEnvelope { outer, inner } = case.acceptance else {
            panic!("{case_id}: fingering rows carry cents envelopes");
        };
        let duct = duct_from_corpus(case, *states);
        let model_hz = first_peak_hz(&duct, &air);
        let measured_hz = case.reference_value_si;
        let cents = 1200.0 * (model_hz / measured_hz).ln() / core::f64::consts::LN_2;
        assert!(
            cents.abs() < outer,
            "{case_id}: {model_hz:.1} Hz vs measured {measured_hz} Hz = {cents:+.1} cents \
             (outer envelope {outer})"
        );
        assert!(
            cents.abs() < inner,
            "{case_id}: post-Dalmont band {inner} cents exceeded ({cents:+.1})"
        );
        assert!(
            model_hz > previous,
            "fingering ladder must rise monotonically: {model_hz:.1} after {previous:.1}"
        );
        previous = model_hz;
        write!(
            rows,
            "{}{{\"case\":\"{case_id}\",\"model_hz\":{model_hz:.1},\"measured_hz\":{measured_hz},\
             \"cents\":{cents:+.1},\"outer\":{outer},\"inner\":{inner}}}",
            if index == 0 { "" } else { "," }
        )
        .expect("write");
    }
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"ernoult-corpus-gate\",\"source\":\"fs-vvreg \
         acoustic corpus (bead 3ez8g.1.1)\",\"rows\":[{rows}],\"verdict\":\"pass\"}}"
    );
}

#[test]
fn corpus_geometry_matches_the_papers_table() {
    // Cross-check: the corpus context must carry exactly the Table 1
    // geometry the inline conformance test encodes locally. If either
    // drifts, ONE of the two tests fails and names the seam.
    let case = acoustic_cases()
        .iter()
        .find(|case| case.id == "acoustic-ernoult-2021-xxxx")
        .expect("registered row");
    let expected: [(&str, f64); 14] = [
        ("bore-radius-m", 2.0e-3),
        ("bore-length-m", 0.2875),
        ("hole-1-position-m", 0.100),
        ("hole-1-radius-m", 1.5e-3),
        ("hole-1-chimney-m", 1.7e-3),
        ("hole-2-position-m", 0.130),
        ("hole-2-radius-m", 1.75e-3),
        ("hole-2-chimney-m", 1.3e-3),
        ("hole-3-position-m", 0.180),
        ("hole-3-radius-m", 1.75e-3),
        ("hole-3-chimney-m", 1.5e-3),
        ("hole-4-position-m", 0.240),
        ("hole-4-radius-m", 1.25e-3),
        ("hole-4-chimney-m", 1.4e-3),
    ];
    for (name, value) in expected {
        assert!(
            (context_value(case, name) - value).abs() < 1.0e-12,
            "context {name} drifted from Table 1"
        );
    }
}

#[test]
fn falsifier_a_detuned_corpus_geometry_fails_the_gate() {
    // The gate must actually discriminate: rebuild the all-closed duct
    // with hole 1 moved 15 mm toward the mouthpiece and confirm the
    // first peak leaves the corpus row's INNER envelope. A gate that
    // cannot fail on a broken input is not a gate.
    use HoleState::Closed as X;
    let case = acoustic_cases()
        .iter()
        .find(|case| case.id == "acoustic-ernoult-2021-oxxx")
        .expect("registered row");
    let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air");
    let AcousticAcceptance::CentsEnvelope { inner, .. } = case.acceptance else {
        panic!("cents envelope expected");
    };
    // Open hole 1 (the oxxx fingering) but at a WRONG position: the open
    // hole dominates the effective length, so mis-placing it detunes hard.
    let bore = 2.0e-3;
    let wrong = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: bore,
                length: 0.085, // paper: 0.100 — 15 mm early
            },
            Segment::ToneHole {
                hole_radius: 1.5e-3,
                chimney_height: 1.7e-3,
                bore_radius: bore,
                state: HoleState::Open,
            },
            Segment::Cylinder {
                radius: bore,
                length: 0.030,
            },
            Segment::ToneHole {
                hole_radius: 1.75e-3,
                chimney_height: 1.3e-3,
                bore_radius: bore,
                state: X,
            },
            Segment::Cylinder {
                radius: bore,
                length: 0.050,
            },
            Segment::ToneHole {
                hole_radius: 1.75e-3,
                chimney_height: 1.5e-3,
                bore_radius: bore,
                state: X,
            },
            Segment::Cylinder {
                radius: bore,
                length: 0.060,
            },
            Segment::ToneHole {
                hole_radius: 1.25e-3,
                chimney_height: 1.4e-3,
                bore_radius: bore,
                state: X,
            },
            Segment::Cylinder {
                radius: bore,
                length: 0.0475,
            },
        ],
    };
    let model_hz = first_peak_hz(&wrong, &air);
    let cents = 1200.0 * (model_hz / case.reference_value_si).ln() / core::f64::consts::LN_2;
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"ernoult-corpus-falsifier\",\"model_hz\":{model_hz:.1},\
         \"cents\":{cents:+.1},\"inner\":{inner},\"verdict\":\"discriminates\"}}"
    );
    assert!(
        cents.abs() > inner,
        "a 15 mm hole-position error must leave the {inner}-cent band (got {cents:+.1}); \
         if it does not, the gate cannot discriminate geometry errors"
    );
}
