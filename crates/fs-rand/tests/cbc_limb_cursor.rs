//! Limb-cursor ratchets for CBC exact arithmetic (bead frankensim-epic-
//! bedrock-6ys.20.6).
//!
//! Every exact primitive is a resumable microprogram with a persistent
//! limb cursor; the tile contract carries an explicit limb-microstep
//! block; and poll spacing is bounded by admitted mutation cells rather
//! than by an unbounded big-integer operation. These tests prove
//! partition-invariant results, typed tamper refusal, persistent tile
//! envelopes across resumed calls, and deterministic replay.

#![deny(unsafe_code)]
use fs_rand::cbc::{CbcBudget, CbcExecutionMode, CbcProblem};
use fs_rand::cbc_exec::{CbcControl, CbcExecError, CbcRunStatus, CbcTileShape, DEFAULT_LIMB_BLOCK};
use fs_rand::cbc_limb::{LimbCursor, StepOutcome, factor_limbs_u32, step_add_multiply};

fn admitted(n: u32, dimension: usize, mode: CbcExecutionMode) -> fs_rand::cbc::CbcAdmission {
    let problem = CbcProblem::new(n, dimension).expect("structural fixture");
    problem
        .admit_for(mode, CbcBudget::UNBOUNDED)
        .expect("fixture admits under an unbounded budget")
}

fn tile(candidate_block: u32, point_block: u32, limb_block: u32) -> CbcTileShape {
    CbcTileShape::with_limbs(candidate_block, point_block, limb_block).expect("nonzero test tile")
}

fn keep_going() -> CbcControl {
    CbcControl::Continue
}

fn drive(executor: &mut fs_rand::cbc_exec::CbcExecutor, shape: CbcTileShape) -> CbcRunStatus {
    executor
        .run(&mut keep_going, shape, u128::MAX)
        .expect("unbounded allowance cannot exhaust")
}

#[test]
fn clc_001_results_are_invariant_across_candidate_point_and_limb_partitions() {
    // One-shot goldens must not depend on how finely execution is sliced:
    // identical lattice prefix, exact work totals, and certificate bytes
    // across every candidate/point/limb partition, including the extreme
    // one-cell-per-poll granularity.
    let golden_partition = tile(1, 1, 1);
    let golden = {
        let mut executor =
            fs_rand::cbc_exec::CbcExecutor::new(admitted(16, 4, CbcExecutionMode::Certified))
                .expect("golden executor admits");
        executor.enable_certificates().expect("certified mode");
        assert_eq!(
            drive(&mut executor, golden_partition),
            CbcRunStatus::Completed
        );
        (
            executor.prefix().to_vec(),
            executor.work_spent(),
            executor.certificates().to_vec(),
        )
    };

    for &(c, p, l) in &[
        (2, 3, 2),
        (4, 8, 64),
        (16, 16, 1),
        (3, 7, DEFAULT_LIMB_BLOCK),
    ] {
        let mut executor =
            fs_rand::cbc_exec::CbcExecutor::new(admitted(16, 4, CbcExecutionMode::Certified))
                .expect("partition executor admits");
        executor.enable_certificates().expect("certified mode");
        assert_eq!(drive(&mut executor, tile(c, p, l)), CbcRunStatus::Completed);
        assert_eq!(executor.prefix(), golden.0.as_slice(), "prefix identity");
        assert_eq!(executor.work_spent(), golden.1, "work totals");
        assert_eq!(
            executor.certificates(),
            golden.2.as_slice(),
            "certificate bytes"
        );
    }
}

#[test]
fn clc_002_poll_frequency_scales_with_limb_granularity() {
    fn polls_under(limb_block: u32) -> u32 {
        let mut executor =
            fs_rand::cbc_exec::CbcExecutor::new(admitted(8, 3, CbcExecutionMode::Construction))
                .expect("executor admits");
        let mut polls = 0_u32;
        let status = executor
            .run(
                &mut || {
                    polls += 1;
                    CbcControl::Continue
                },
                tile(4, 2, limb_block),
                u128::MAX,
            )
            .expect("unbounded run");
        assert_eq!(status, CbcRunStatus::Completed);
        polls
    }

    let fine = polls_under(1);
    let coarse = polls_under(DEFAULT_LIMB_BLOCK);
    assert!(
        fine > coarse * 2,
        "one-cell granularity must interleave far more observations \
         (fine={fine}, coarse={coarse})"
    );
}

#[test]
fn clc_003_cancellation_at_microstep_boundaries_never_half_commits() {
    // Golden one-shot reference.
    let golden = {
        let mut executor =
            fs_rand::cbc_exec::CbcExecutor::new(admitted(8, 3, CbcExecutionMode::Certified))
                .expect("golden executor admits");
        executor.enable_certificates().expect("certified mode");
        assert_eq!(drive(&mut executor, tile(4, 2, 1)), CbcRunStatus::Completed);
        (
            executor.prefix().to_vec(),
            executor.work_spent(),
            executor.certificates().to_vec(),
        )
    };

    // Cancel at EVERY poll with one-cell tiles: pauses land inside exact
    // operations. Public prefixes stay whole-component throughout and the
    // resumed result replays the golden byte-for-byte.
    let mut executor =
        fs_rand::cbc_exec::CbcExecutor::new(admitted(8, 3, CbcExecutionMode::Certified))
            .expect("adversarial executor admits");
    executor.enable_certificates().expect("certified mode");
    let mut guard = 0_u32;
    loop {
        let mut first = true;
        let status = executor
            .run(
                &mut || {
                    if first {
                        first = false;
                        CbcControl::Cancel
                    } else {
                        CbcControl::Continue
                    }
                },
                tile(4, 2, 1),
                u128::MAX,
            )
            .expect("adversarial run cannot refuse storage");
        match status {
            CbcRunStatus::Completed => break,
            CbcRunStatus::Cancelled(_) => {
                // Whole-component visibility: the public prefix only ever
                // contains committed generator components.
                assert!(
                    executor.prefix().len() <= 3,
                    "prefix never escapes the admitted dimension"
                );
            }
            CbcRunStatus::AllowanceExhausted(_) => {}
        }
        guard += 1;
        assert!(guard < 10_000, "adversarial stepping failed to converge");
    }
    assert_eq!(executor.prefix(), golden.0.as_slice());
    assert_eq!(executor.work_spent(), golden.1);
    assert_eq!(executor.certificates().to_vec(), golden.2);
    let _ = guard;
}

#[test]
fn clc_004_progress_receipt_exposes_positions_without_buffers() {
    let mut executor =
        fs_rand::cbc_exec::CbcExecutor::new(admitted(8, 3, CbcExecutionMode::Construction))
            .expect("executor admits");

    // Before any scanning there is no pending microprogram.
    assert!(executor.micro_progress().is_none());

    // Drive into the scan phase and pause mid-operation using one-cell
    // tiles: the receipt must expose a pending cursor position.
    let mut saw_pending = false;
    loop {
        let mut first = true;
        let status = executor
            .run(
                &mut || {
                    if first {
                        first = false;
                        CbcControl::Cancel
                    } else {
                        CbcControl::Continue
                    }
                },
                tile(4, 2, 1),
                u128::MAX,
            )
            .expect("progress probing run");
        if let Some(progress) = executor.micro_progress() {
            if let Some(pending) = progress.pending {
                saw_pending = true;
                assert!(
                    pending.dst_pos < 4096,
                    "cursor positions are logical, never buffer identities"
                );
            }
            assert!(
                progress.scan_tile_candidates <= 4,
                "persistent counter respects the tile envelope"
            );
        }
        match status {
            CbcRunStatus::Completed => break,
            _ => {}
        }
    }
    assert!(saw_pending, "at least one pause observed a pending cursor");
    assert!(
        executor.micro_progress().is_none(),
        "completed runs park cleanly"
    );
}

#[test]
fn clc_005_tampered_or_stale_cursors_are_refused_without_mutation() {
    let mut dst: Vec<u32> = vec![7];
    let (words, flen) = factor_limbs_u32(3);

    // Version skew: a cursor from a foreign encoding generation refuses.
    let stale = LimbCursor {
        version: fs_rand::cbc_limb::LIMB_CURSOR_VERSION + 1,
        ..LimbCursor::begin_add_multiply(1, flen, dst.len(), 1 + flen + 1, usize::MAX)
    };
    let mut probe = stale;
    assert_eq!(
        step_add_multiply(&mut dst, &[1], &words, flen, &mut probe, 8),
        StepOutcome::Refused
    );

    // Out-of-bounds source position refuses rather than clamping.
    let beyond = LimbCursor {
        src_pos: 9,
        ..LimbCursor::begin_add_multiply(1, flen, dst.len(), 1 + flen + 1, usize::MAX)
    };
    let mut probe = beyond;
    assert_eq!(
        step_add_multiply(&mut dst, &[1], &words, flen, &mut probe, 8),
        StepOutcome::Refused
    );
    assert_eq!(dst, vec![7], "refusals never mutate the destination");
}

#[test]
fn clc_006_tile_admission_rejects_zero_limb_blocks() {
    assert!(CbcTileShape::new(1, 1).is_ok());
    assert!(CbcTileShape::with_limbs(1, 1, 1).is_ok());
    assert_eq!(
        CbcTileShape::with_limbs(1, 1, 0),
        Err(CbcExecError::InvalidLimbBlock { limb_block: 0 })
    );
}
