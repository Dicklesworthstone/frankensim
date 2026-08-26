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
//! Run with `cargo test -p fs-aeroac --test slot_jet_3d_sweep --release
//! -- --ignored --nocapture`. Heavy: each rung settles then records;
//! budget wall time accordingly and let the checkpoint machinery
//! (chunked runner) absorb worker job walls if routed through RCH.

use fs_aeroac::slot_jet_3d::{
    SlotJet3dConfig, SweepProgress, classify_rung, run_slot_jet_3d_chunked,
};
use fs_lbm::d3q19::CollisionModel3;
use std::fmt::Write as _;

/// Fixed geometry for the whole ladder (the clean actuator is the
/// collision rate at fixed `u_jet`, per the executed ramp protocol).
fn base_config(second_order_rate: f64) -> SlotJet3dConfig {
    SlotJet3dConfig {
        nx: 96,
        ny: 48,
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
        steps_settle: 8_000,
        steps_record: 8_192,
    }
}

/// Ladder from the sub-tonal floor to the central-moment stability
/// edge. Higher rate = lower nu = higher jet Re.
const LADDER: [f64; 6] = [1.60, 1.75, 1.85, 1.92, 1.96, 1.98];

/// Optional single-rung filter for distributed execution:
/// `SWEEP_ONLY_RUNG=<index into LADDER>` runs just that rung and
/// archives to a per-rung receipt file, so independent machines can
/// take different rungs concurrently. Unset = the full serial ladder
/// (the default for one dedicated host).
fn only_rung() -> Option<usize> {
    std::env::var("SWEEP_ONLY_RUNG")
        .ok()
        .and_then(|v| v.parse().ok())
}

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

    let only = only_rung();
    let mut filtered_rate = Option::<f64>::None;
    for (rung_idx, rate) in LADDER.iter().enumerate() {
        if let Some(o) = only {
            if rung_idx != o {
                continue;
            }
            filtered_rate = Some(*rate);
        }
        let cfg = base_config(*rate);
        let started = std::time::Instant::now();
        // Chunked execution with an atomic checkpoint: repeated
        // invocations (or post-crash reruns) resume bit-identically.
        // Checkpoints live INSIDE the synced tree (gitignored) so
        // repeated RCH invocations resume on any worker; temp_dir is
        // worker-local and loses state across routing changes.
        let ckpt = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../target-sweep-ckpt")
            .join(format!("{rate:.2}"));
        let run = loop {
            match run_slot_jet_3d_chunked(&cfg, &ckpt, 1_024).expect("chunk executes") {
                SweepProgress::Complete(run) => break *run,
                SweepProgress::Partial { steps_done } => {
                    println!("  rate {rate:.2}: {steps_done} steps done");
                }
            }
        };
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
        let _ = writeln!(
            regime_map,
            "{:.3}\t{:.1}\t{:.3e}\t{}\t{:.4}\t{}\t{:.2e}\t{:.3e}\t{}\t{:.4}\t{:.3}",
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
        );
    }

    // Box-sensitivity octave (DONE-WHEN): double the spanwise extent
    // on the highest-Re rung — 2-D-ness is exactly what we are testing
    // our way out of, so a too-thin box must be ruled out explicitly.
    // Skipped in single-rung distributed mode unless that rung IS the
    // highest-Re one (its owner also owns the octave).
    if only.is_some_and(|o| o + 1 != LADDER.len()) {
        println!("octave skipped in single-rung mode");
        return;
    }
    let mut cfg_hi = base_config(*LADDER.last().expect("non-empty ladder"));
    cfg_hi.nz *= 2;
    let ckpt_hi =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target-sweep-ckpt/octave");
    let run_hi = loop {
        match run_slot_jet_3d_chunked(&cfg_hi, &ckpt_hi, 1_024).expect("chunk executes") {
            SweepProgress::Complete(run) => break *run,
            SweepProgress::Partial { steps_done } => {
                println!("  octave: {steps_done} steps done");
            }
        }
    };
    let cls_hi = classify_rung(&run_hi, &cfg_hi).expect("octave classification");
    jsonl.push_str(&cls_hi.to_jsonl());
    jsonl.push('\n');
    println!(
        "octave\tnz={}\tRe={:.1}\tflatness={:.3e}\ttonal={}\tSt={:.4}",
        cfg_hi.nz, cls_hi.reynolds, cls_hi.flatness, cls_hi.tonal, cls_hi.strouhal
    );

    if let Some(rate) = filtered_rate {
        let rung_path =
            receipt_path.with_file_name(format!("slot-jet-3d-re-sweep-rung{rate:.2}.jsonl"));
        std::fs::write(&rung_path, &jsonl).expect("archive per-rung receipts");
        println!("archived: {}", rung_path.display());
    } else {
        std::fs::write(&receipt_path, &jsonl).expect("archive sweep receipts");
    }
    println!(
        "\nREGIME MAP (rate/Re/flatness/tonal/St/bin/prominence/rms/qualified/mach/imbalance):\n{regime_map}"
    );
    println!("archived: {}", receipt_path.display());
}
