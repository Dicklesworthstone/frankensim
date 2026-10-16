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
    RE_LADDER, SlotJet3dConfig, SweepProgress, classify_rung, classify_rung_parity_filtered,
    ladder_config, load_completed_run, run_slot_jet_3d_chunked,
};
use fs_lbm::d3q19::CollisionModel3;
use std::fmt::Write as _;

/// Fixed geometry for the whole ladder (the clean actuator is the
/// collision rate at fixed `u_jet`, per the executed ramp protocol).
fn base_config(second_order_rate: f64) -> SlotJet3dConfig {
    // The rig definition lives in the library so the card minter's
    // provenance fingerprint and this driver cannot drift apart.
    ladder_config(second_order_rate)
}

/// Ladder from the sub-tonal floor to the central-moment stability
/// edge. Higher rate = lower nu = higher jet Re.
const LADDER: [f64; 6] = RE_LADDER;

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
        // HONEST OUTCOMES (bead law): a rung that destabilizes inside the
        // rig's own containment (typed NonFinite refusal from the density
        // guard) is a RESULT — the stability edge of this operator family —
        // and is recorded as a typed refusal row rather than aborting the
        // campaign. Executed: rate 1.92 with the clamped nuisance rate 1.99
        // blew up during settle (rho -0.244) on 2026-09-01 while 1.96 and
        // 1.98 stayed stable — a nonmonotonic edge worth mapping, not hiding.
        let outcome = loop {
            match run_slot_jet_3d_chunked(&cfg, &ckpt, 1_024) {
                Ok(SweepProgress::Complete(run)) => break Ok(*run),
                Ok(SweepProgress::Partial { steps_done }) => {
                    println!("  rate {rate:.2}: {steps_done} steps done");
                }
                Err(refusal) => break Err(refusal),
            }
        };
        let run = match outcome {
            Ok(run) => run,
            // Only the rig's own containment refusal is a RESULT. Any
            // other error (checkpoint I/O on a full disk, a fingerprint
            // mismatch, config refusal) is infrastructure and must abort
            // loudly rather than be archived as a physics refusal row —
            // executed hazard 2026-09-01: yto's root filesystem hit 100%
            // while three rungs were 800 steps from their terminal write.
            Err(refusal @ fs_aeroac::AeroacError::NonFinite { .. }) => {
                let CollisionModel3::CentralMoment {
                    second_order_rate,
                    higher_order_rate,
                } = cfg.collision
                else {
                    unreachable!("base_config pins the central-moment operator");
                };
                let row = format!(
                    "{{\"schema\":\"fs-aeroac.slot-jet-3d.rung-refusal/v1\",\
                     \"second_order_rate\":{second_order_rate},\
                     \"higher_order_rate\":{higher_order_rate},\
                     \"refusal\":\"{refusal}\"}}"
                );
                jsonl.push_str(&row);
                jsonl.push('\n');
                println!("rung\t{rate:.3}\tREFUSED\t{refusal}");
                let _ = writeln!(
                    regime_map,
                    "{rate:.3}\tREFUSED: {refusal}\t-\t-\t-\t-\t-\t-\t-\t-\t-"
                );
                continue;
            }
            Err(other) => panic!(
                "rate {rate:.2}: infrastructure error, not a rung result (nothing archived): {other}"
            ),
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
    // highest-Re one (its owner also owns the octave). The archive
    // still happens: the pre-fix early return here silently discarded
    // every non-octave single-rung receipt (executed on 2026-09-01 —
    // the recorded lesson behind the shared archive path below).
    if only.is_some_and(|o| o + 1 != LADDER.len()) {
        println!("octave skipped in single-rung mode");
        archive_receipts(&receipt_path, filtered_rate, &jsonl, &regime_map);
        return;
    }
    let mut cfg_hi = base_config(*LADDER.last().expect("non-empty ladder"));
    cfg_hi.nz *= 2;
    let ckpt_hi =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target-sweep-ckpt/octave");
    let octave_outcome = loop {
        match run_slot_jet_3d_chunked(&cfg_hi, &ckpt_hi, 1_024) {
            Ok(SweepProgress::Complete(run)) => break Ok(*run),
            Ok(SweepProgress::Partial { steps_done }) => {
                println!("  octave: {steps_done} steps done");
            }
            Err(refusal) => break Err(refusal),
        }
    };
    match octave_outcome {
        Ok(run_hi) => {
            let cls_hi = classify_rung(&run_hi, &cfg_hi).expect("octave classification");
            jsonl.push_str(&cls_hi.to_jsonl());
            jsonl.push('\n');
            println!(
                "octave\tnz={}\tRe={:.1}\tflatness={:.3e}\ttonal={}\tSt={:.4}",
                cfg_hi.nz, cls_hi.reynolds, cls_hi.flatness, cls_hi.tonal, cls_hi.strouhal
            );
        }
        Err(other) if !matches!(other, fs_aeroac::AeroacError::NonFinite { .. }) => {
            panic!("octave: infrastructure error, not a rung result (nothing archived): {other}")
        }
        Err(refusal) => {
            let row = format!(
                "{{\"schema\":\"fs-aeroac.slot-jet-3d.rung-refusal/v1\",\
                 \"octave\":true,\"nz\":{},\"refusal\":\"{refusal}\"}}",
                cfg_hi.nz
            );
            jsonl.push_str(&row);
            jsonl.push('\n');
            println!("octave\tnz={}\tREFUSED\t{refusal}", cfg_hi.nz);
            let _ = writeln!(
                regime_map,
                "octave\tREFUSED: {refusal}\t-\t-\t-\t-\t-\t-\t-\t-\t-"
            );
        }
    }

    archive_receipts(&receipt_path, filtered_rate, &jsonl, &regime_map);
}

/// Post-hoc re-classification of every completed rung under the
/// disclosed parity filter (bead law: logging sufficient to re-derive
/// the classification). Terminal checkpoints retain the raw force
/// record and `run_slot_jet_3d_chunked` returns a Complete run
/// idempotently, so this reads physics that already executed and
/// runs no lattice step. Archives one receipt file with, per rung, the
/// raw row, the parity-filtered row, and the Nyquist-edge verdict of
/// the raw peak (the artifact class the 2-D battery pins at
/// `bin < n/2 - 8`). Rungs without a Complete checkpoint are skipped
/// and named.
#[test]
#[ignore = "post-hoc analysis of executed sweep checkpoints (bead 3ez8g.10.1)"]
fn re_sweep_parity_reclassify() {
    let receipt_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/receipts/slot-jet-3d-re-sweep-parity.jsonl");
    let mut jsonl = String::from(
        "{\"schema\":\"fs-aeroac.slot-jet-3d.re-sweep-parity/v1\",\"scope\":\"per completed rung: raw row, parity-filtered row (pair-averaged record, bin width doubled), and the raw peak's Nyquist-edge verdict; no lattice step executed\"}\n",
    );
    let mut any = false;
    for rate in LADDER {
        let cfg = base_config(rate);
        let ckpt = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../target-sweep-ckpt")
            .join(format!("{rate:.2}"));
        // Load-only: a checkpoint a live rung still owns is never
        // stepped or rewritten by the analyser.
        let run = match load_completed_run(&cfg, &ckpt) {
            Ok(Some(run)) => run,
            Ok(None) => {
                println!("rate {rate:.2}: no terminal checkpoint; skipped");
                continue;
            }
            Err(refusal) => {
                println!("rate {rate:.2}: checkpoint refused ({refusal}); skipped");
                continue;
            }
        };
        let raw = classify_rung(&run, &cfg).expect("raw classification");
        let filtered = classify_rung_parity_filtered(&run, &cfg).expect("filtered classification");
        let nyquist_bins = run.diagnostics.record_len / 2;
        let raw_edge = raw.peak_bin + 8 >= nyquist_bins;
        // A filtered peak sitting exactly on the pipeline's low-bin
        // guard (n'/8) means the admitted spectrum is monotone from the
        // guard down: drift, not a tone. Recorded, not asserted.
        let filtered_guard = filtered.peak_bin == run.diagnostics.record_len / 2 / 8;
        let _ = writeln!(
            jsonl,
            "{{\"schema\":\"fs-aeroac.slot-jet-3d.rung-parity/v1\",\"second_order_rate\":{},\"raw_peak_at_nyquist_edge\":{raw_edge},\"filtered_peak_at_guard_floor\":{filtered_guard},\"raw\":{},\"filtered\":{}}}",
            rate,
            raw.to_jsonl(),
            filtered.to_jsonl()
        );
        println!(
            "rate {rate:.2}\tRe={:.1}\traw: bin {} of {} (edge={raw_edge}) flatness={:.3e} tonal={}\tfiltered: St={:.4} bin {} (guard_floor={filtered_guard}) flatness={:.3e} tonal={} prominence={:.2e}",
            raw.reynolds,
            raw.peak_bin,
            nyquist_bins,
            raw.flatness,
            raw.tonal,
            filtered.strouhal,
            filtered.peak_bin,
            filtered.flatness,
            filtered.tonal,
            filtered.prominence
        );
        any = true;
    }
    assert!(any, "no completed rung checkpoint to re-classify");
    std::fs::write(&receipt_path, jsonl).expect("archive parity receipts");
    println!("archived: {}", receipt_path.display());
}

/// Shared terminal archive: the per-rung file in single-rung mode, the
/// full campaign file otherwise, plus the printed regime map.
fn archive_receipts(
    receipt_path: &std::path::Path,
    filtered_rate: Option<f64>,
    jsonl: &str,
    regime_map: &str,
) {
    if let Some(rate) = filtered_rate {
        let rung_path =
            receipt_path.with_file_name(format!("slot-jet-3d-re-sweep-rung{rate:.2}.jsonl"));
        std::fs::write(&rung_path, jsonl).expect("archive per-rung receipts");
        println!("archived: {}", rung_path.display());
    } else {
        std::fs::write(receipt_path, jsonl).expect("archive sweep receipts");
        println!("archived: {}", receipt_path.display());
    }
    println!(
        "\nREGIME MAP (rate/Re/flatness/tonal/St/bin/prominence/rms/qualified/mach/imbalance):\n{regime_map}"
    );
}
