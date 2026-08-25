//! The 3-D slot-jet Re sweep campaign (music bead
//! `frankensim-music-v8-root-3ez8g.10.1`).
//!
//! Drives the staged-but-never-executed 3-D ladder on the D3Q19
//! central-moment operator at fixed geometry, classifies every rung
//! through the shared measurement pipeline, archives one JSONL row per
//! rung into `receipts/slot-jet-3d-re-sweep.jsonl`, and prints the
//! regime map. HONEST OUTCOMES (both are wins): a demonstrated
//! broadband rung with its Re/resolution boundary, OR a quantified
//! refusal extending the recorded 2-D tonal-lock refusal. Neither
//! outcome mints an experimental claim; [`crate::SCOPE_STATEMENT`]
//! travels in every run.
//!
//! STATUS (2026-08-25): type-check PENDING — this file was written
//! against source-verified APIs but never compiled: the host hit
//! critical memory pressure and the RCH fleet refused all slots at
//! once. First runner with a fleet slot: `cargo check -p fs-aeroac
//! --test slot_jet_3d_sweep` BEFORE running anything; on failure,
//! revert the introducing commit rather than patching blind.
//!
//! Run with `cargo test -p fs-aeroac --test slot_jet_3d_sweep --release
//! -- --ignored --nocapture`. Heavy: each rung settles then records;
//! budget wall time accordingly and let the checkpoint machinery
//! (chunked runner) absorb worker job walls if routed through RCH.

use fs_aeroac::slot_jet_3d::{SlotJet3dConfig, classify_rung, run_slot_jet_3d};
use fs_lbm::d3q19::CollisionModel3;

/// Fixed geometry for the whole ladder (the clean actuator is the
/// collision rate at fixed `u_jet`, per the executed ramp protocol).
fn base_config(second_order_rate: f64) -> SlotJet3dConfig {
    SlotJet3dConfig {
        nx: 80,
        ny: 40,
        nz: 12,
        slot_half: 2.5,
        u_jet: 0.04,
        collision: CollisionModel3::CentralMoment {
            second_order_rate,
            higher_order_rate: second_order_rate + 0.1,
        },
        nozzle_thickness: 1,
        edge_distance: 8,
        plate_length: 16,
        fringe_width: 12,
        fringe_sigma: 0.4,
        seed_amplitude: 0.05,
        steps_settle: 4_000,
        steps_record: 4_096,
    }
}

/// Ladder from the sub-tonal floor to the central-moment stability
/// edge. Higher rate = lower nu = higher jet Re.
const LADDER: [f64; 6] = [1.60, 1.75, 1.85, 1.92, 1.96, 1.98];

#[test]
#[ignore = "heavy: full 3-D Re sweep campaign (bead 3ez8g.10.1)"]
fn re_sweep_campaign() {
    let receipt_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/receipts/slot-jet-3d-re-sweep.jsonl");
    if let Some(parent) = receipt_path.parent() {
        std::fs::create_dir_all(parent).expect("receipt dir");
    }

    let mut jsonl = String::from(
        "{\"schema\":\"fs-aeroac.slot-jet-3d.re-sweep/v1\",\"scope\":\"campaign header; per-rung rows follow\"}\n",
    );
    let mut regime_map = String::from(
        "runr\tRe\trate\tflatness\ttonal\tStrouhal\tpeak_bin\tprominence\tforce_rms\tamplitude_qualified\tmach_max\tflux_imbalance\n",
    );

    for rate in LADDER {
        let cfg = base_config(rate);
        let started = std::time::Instant::now();
        let run = run_slot_jet_3d(&cfg).expect("rung executes");
        let rung = classify_rung(&run, &cfg).expect("classification succeeds");
        jsonl.push_str(&rung.to_jsonl());
        jsonl.push('\n');
        println!(
            "rung\t{:.3}\tRe={:.1}\tflatness={:.3e}\ttonal={}\tSt={:.4}\tprominence={:.2e}\trms={:.3e}\tqualified={}\tmach={:.4}\tseconds\t{:.1}",
            rate,
            rung.reynolds,
            rung.flatness,
            rung.tonal,
            rung.strouhal,
            rung.prominence,
            rung.force_rms,
            rung.amplitude_qualified,
            run.diagnostics.mach_max_lattice,
            started.elapsed().as_secs_f64()
        );
        regime_map.push_str(&format!(
            "{:.3}\t{:.1}\t{:.3e}\t{}\t{:.4}\t{}\t{:.2e}\t{:.3e}\t{}\t{:.4}\t{:.3}\n",
            rate,
            rung.reynolds,
            rung.flatness,
            rung.tonal,
            rung.strouhal,
            rung.peak_bin,
            rung.prominence,
            rung.force_rms,
            rung.amplitude_qualified,
            run.diagnostics.mach_max_lattice,
            rung.flux_imbalance
        ));
    }

    std::fs::write(&receipt_path, &jsonl).expect("archive sweep receipts");
    println!(
        "\nREGIME MAP (rate/Re/flatness/tonal/St/bin/prominence/rms/qualified/mach/imbalance):\n{regime_map}"
    );
    println!("archived: {}", receipt_path.display());
}
