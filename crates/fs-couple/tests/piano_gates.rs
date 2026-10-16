//! Piano claim gates (music bead `frankensim-music-v8-root-3ez8g.5.3`):
//! the simulate-vs-authority board cross-check, the out-of-band felt
//! honesty refusal, and the gate-summary enumeration — the review that
//! converts the piano vertical's executed physics into registry
//! evidence. The tilt-vs-COUPON gate stays NAMED (the felt coupon's
//! retention is blocked — absence row `acoustic-absent-felt-coupon`);
//! the tilt-vs-velocity MECHANISM gates live in the felt island and
//! pv-001.

use fs_couple::piano_vertical::{HammerLaw, PedalState, PianoStringSpec, PianoVertical};
use fs_modal::SliceOptions;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh, PlateSection, assemble, modes};

fn base_spec() -> PianoStringSpec {
    PianoStringSpec {
        f0_hz: 220.0,
        b_inharmonicity: 3.5e-4,
        detune_cents: 0.0,
        n_modes: 12,
        damping_ratio: 4.0e-4,
    }
}

#[test]
fn pg_001_board_from_the_fs_plate_authority() {
    // The vertical's board modes come FROM the fs-plate modal authority
    // (an isotropic spruce-scale panel; the orthotropic card upgrade
    // rides the same seam), and the Weinreich aftersound still emerges
    // on the authority-derived board — simulate-vs-authority, executed.
    let section = PlateSection::isotropic(10.0e9, 0.35, 0.009, 420.0).expect("spruce-ish");
    let mesh = PlateMesh::rectangle(1.0, 0.6, 12, 8);
    let boundary = PlateMesh::rectangle_boundary(12, 8);
    let model = assemble(
        &mesh,
        &section,
        &boundary,
        &[],
        &AssemblyOptions {
            pretension: 0.0,
            support: EdgeSupport::SimplySupported,
        },
    )
    .expect("assemble");
    let tau = core::f64::consts::TAU;
    let window = ((tau * 40.0f64).powi(2), (tau * 500.0f64).powi(2));
    let report = modes(&model, window, &SliceOptions::default()).expect("modes");
    assert!(
        report.modes.len() >= 4,
        "expected at least four board modes in 40..500 Hz, got {}",
        report.modes.len()
    );
    let board: Vec<(f64, f64, f64)> = report
        .modes
        .iter()
        .take(4)
        .enumerate()
        .map(|(i, m)| {
            (
                m.lambda.sqrt() / tau,
                0.02 + 0.005 * i as f64, // authored loss (no measured board Q)
                1.0 / (1.0 + i as f64 * 0.4),
            )
        })
        .collect();
    let f_first = board[0].0;
    let mut pv = PianoVertical::new_with_board(
        base_spec(),
        1.2,
        HammerLaw::Felt,
        PedalState {
            sustain: true,
            una_corda: false,
        },
        0.03,
        &board,
    )
    .expect("vertical with the authority board");
    pv.strike(2.0);
    let mut ensemble = Vec::new();
    for k in 0..48_000usize {
        let _ = pv.step();
        if k % 240 == 0 {
            let e = pv.string_energies();
            ensemble.push((e[0] + e[1] + e[2]).max(1e-300));
        }
    }
    let rate = |w0: usize, w1: usize| -> f64 {
        (ensemble[w0] / ensemble[w1]).ln() / ((w1 - w0) as f64 * 240.0 / 24_000.0)
    };
    let early = rate(10, 40);
    let late = rate(120, 190);
    let two_stage = early > 1.6 * late;
    assert!(
        two_stage,
        "the aftersound must survive the authority board (early {early:.2} late {late:.2})"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"pg-001-board-authority\",\"verdict\":\"pass\",\
         \"board_modes_hz\":[{:.1},{:.1},{:.1},{:.1}],\"first_mode_hz\":{f_first:.1},\
         \"early_decay\":{early:.2},\"late_decay\":{late:.2}}}",
        board[0].0, board[1].0, board[2].0, board[3].0
    );
}

#[test]
fn pg_002_hammer_band_out_of_band_refuses() {
    // Out-of-band honesty (scope 5): a felt Prony chain certified over
    // a band that does NOT cover the hammer spectrum REFUSES by name at
    // the hammer's own frequencies — never extrapolates.
    use fs_material::visco::{FractionalZener, ViscoError, lower_to_prony};
    let fz = FractionalZener::new(4.0e5, 1.2e6, 0.35, 1.0e-3).expect("FZ");
    // A band deliberately BELOW the hammer spectrum's top octave.
    let lowered = lower_to_prony(&fz, 50.0, 1_000.0, 12, 0.05).expect("narrow band certifies");
    let tau = core::f64::consts::TAU;
    // The 87zbd strike's contact time (~1.5-3 ms) puts real content at
    // 2-5 kHz: querying there must refuse.
    let in_band = lowered.modulus_checked(tau * 500.0);
    let out_hammer = lowered.modulus_checked(tau * 3_000.0);
    assert!(in_band.is_ok(), "mid-band evaluates");
    assert!(
        matches!(out_hammer, Err(ViscoError::OutOfBand { .. })),
        "a hammer spectrum exceeding the certified band must refuse, got {out_hammer:?}"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"pg-002-out-of-band\",\"verdict\":\"pass\",\
         \"band_hz\":[50,1000],\"refused_at_hz\":3000}}"
    );
}

#[test]
fn pg_003_piano_gate_summary_enumerates_every_row() {
    // The committed piano gate summary must name every hammer-felt row
    // (the piano filling's registry spelling) and carry the honest
    // statuses; the air-column image is a STRUCTURAL NON-ROW for this
    // filling (asserted absent — it was never legal, so it is not a
    // refused row either).
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf();
    let registry = std::fs::read_to_string(root.join("instrument-claims.json")).expect("registry");
    let summary = std::fs::read_to_string(root.join("data/claims/piano-gate-summary.tsv"))
        .expect("committed piano gate summary");
    assert!(summary.starts_with("# frankensim-piano-gate-summary-v1"));
    let mut images = Vec::new();
    let mut cursor = 0usize;
    while let Some(hit) = registry[cursor..].find("\"filling\": \"hammer-felt\"") {
        let at = cursor + hit;
        let tail = &registry[at..];
        let tag = "\"image\": \"";
        let img_at = tail.find(tag).expect("image") + tag.len();
        let img_end = tail[img_at..].find('"').expect("end");
        images.push(tail[img_at..img_at + img_end].to_string());
        cursor = at + tag.len();
    }
    assert!(images.len() >= 4, "the piano menu: {images:?}");
    for image in &images {
        assert!(summary.contains(image.as_str()), "summary omits {image}");
    }
    // The structural non-row: no air-column image for the piano filling
    // anywhere in the registry.
    assert!(
        !images
            .iter()
            .any(|i| i.contains("air") || i.contains("column")),
        "an air-column image on the piano filling was never legal"
    );
    // The linear-spring debug row exists with its narrowed claim.
    assert!(
        summary.contains("linear-spring\tdebug-control"),
        "the linear-spring debug row must be enumerated"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"pg-003-gate-summary\",\"verdict\":\"pass\",\
         \"rows\":{}}}",
        images.len()
    );
}
