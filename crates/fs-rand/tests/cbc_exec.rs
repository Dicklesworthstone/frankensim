//! G0/G3/G4-class battery for the tiled exact-CBC executor (bead 6ys.20,
//! execution tranche): byte-equivalence with the synchronous authority,
//! tile-shape and pause/resume invariance, request-drain-finalize
//! cancellation with no half-committed component, allowance exhaustion at
//! named boundaries, zero-allowance no-work, and schedule conformance
//! against the admission estimate.

use fs_rand::cbc::{CbcBudget, CbcProblem};
use fs_rand::cbc_exec::{
    CBC_EXECUTOR_SCHEMA_VERSION, CbcBoundary, CbcControl, CbcExecError, CbcExecutor, CbcRunStatus,
    CbcTileShape,
};
use fs_rand::qmc::Lattice;

const CASES: [(u32, usize); 5] = [(5, 3), (8, 3), (16, 4), (127, 4), (257, 5)];

fn admitted_executor(n: u32, dim: usize) -> CbcExecutor {
    let problem = CbcProblem::new(n, dim).expect("battery cases are structurally valid");
    let admission = problem
        .admit(CbcBudget::UNBOUNDED)
        .expect("battery cases admit under the unbounded budget");
    CbcExecutor::new(admission)
}

fn run_to_completion(executor: &mut CbcExecutor, tile: CbcTileShape) {
    let mut keep_going = || CbcControl::Continue;
    let status = executor
        .run(&mut keep_going, tile, u128::MAX)
        .expect("an admitted run must not breach its own schedule");
    assert_eq!(status, CbcRunStatus::Completed);
}

#[test]
fn cbx_001_executor_matches_synchronous_authority() {
    for (n, dim) in CASES {
        let reference = Lattice::cbc(n, dim);
        let mut executor = admitted_executor(n, dim);
        run_to_completion(
            &mut executor,
            CbcTileShape::new(3, 7).expect("static tile shape"),
        );
        let lattice = executor
            .into_lattice()
            .expect("completed construction yields the lattice");
        assert_eq!(
            lattice.z, reference.z,
            "n={n} dim={dim}: executor and Lattice::cbc disagree"
        );
        assert_eq!(lattice.n, reference.n);
    }
    // The independently pinned n=257 KAT from the exact-CBC repair.
    let mut executor = admitted_executor(257, 5);
    run_to_completion(
        &mut executor,
        CbcTileShape::new(64, 64).expect("static tile shape"),
    );
    assert_eq!(
        executor
            .into_lattice()
            .expect("completed construction yields the lattice")
            .z,
        vec![1, 71, 56, 106, 21],
        "the pinned n=257 generator vector moved"
    );
}

#[test]
fn cbx_002_tile_shape_never_changes_bytes_or_debits() {
    for (n, dim) in [(8, 3), (127, 4)] {
        let mut reference_z: Option<(Vec<u32>, u128)> = None;
        for (candidate_block, point_block) in [(1, 1), (3, 7), (64, 64), (1024, 1024)] {
            let mut executor = admitted_executor(n, dim);
            run_to_completion(
                &mut executor,
                CbcTileShape::new(candidate_block, point_block).expect("nonzero tile shape"),
            );
            let spent = executor.work_spent();
            let z = executor
                .into_lattice()
                .expect("completed construction yields the lattice")
                .z;
            match &reference_z {
                None => reference_z = Some((z, spent)),
                Some((reference, reference_spent)) => {
                    assert_eq!(
                        &z, reference,
                        "n={n} dim={dim} tile ({candidate_block},{point_block}) changed bytes"
                    );
                    assert_eq!(
                        spent, *reference_spent,
                        "n={n} dim={dim} tile ({candidate_block},{point_block}) changed debits"
                    );
                }
            }
        }
    }
}

#[test]
fn cbx_003_cancellation_drains_and_never_half_commits() {
    let reference = Lattice::cbc(16, 4);
    for cancel_after in [0_u32, 1, 2, 5, 11, 23] {
        let mut executor = admitted_executor(16, 4);
        let tile = CbcTileShape::new(2, 4).expect("static tile shape");
        let mut polls = 0_u32;
        let mut poll = move || {
            polls += 1;
            if polls > cancel_after {
                CbcControl::Cancel
            } else {
                CbcControl::Continue
            }
        };
        let status = executor
            .run(&mut poll, tile, u128::MAX)
            .expect("cancellation is a status, not an error");
        match status {
            CbcRunStatus::Cancelled(_) => {
                // The committed prefix is a byte-exact prefix of the final
                // vector: no half-committed component exists.
                let prefix = executor.prefix().to_vec();
                assert!(
                    prefix.len() <= reference.z.len(),
                    "prefix cannot exceed the dimension"
                );
                assert_eq!(
                    prefix.as_slice(),
                    &reference.z[..prefix.len()],
                    "cancel_after={cancel_after}: committed prefix diverged"
                );
                // Resume with no further cancellation: identical final bytes.
                let mut keep_going = || CbcControl::Continue;
                let resumed = executor
                    .run(&mut keep_going, tile, u128::MAX)
                    .expect("resume after cancellation");
                assert_eq!(resumed, CbcRunStatus::Completed);
                assert_eq!(
                    executor
                        .into_lattice()
                        .expect("completed construction yields the lattice")
                        .z,
                    reference.z,
                    "cancel_after={cancel_after}: resume diverged from one-shot"
                );
            }
            CbcRunStatus::Completed => {
                assert!(
                    cancel_after > 0,
                    "an immediate cancel must not complete the construction"
                );
            }
            CbcRunStatus::AllowanceExhausted(_) => {
                panic!("an unbounded allowance cannot exhaust")
            }
        }
    }
}

#[test]
fn cbx_004_allowance_pause_resume_equals_one_shot() {
    let reference = Lattice::cbc(127, 3);
    let mut one_shot = admitted_executor(127, 3);
    run_to_completion(
        &mut one_shot,
        CbcTileShape::new(5, 32).expect("static tile shape"),
    );
    let one_shot_spent = one_shot.work_spent();

    let mut executor = admitted_executor(127, 3);
    let tile = CbcTileShape::new(5, 32).expect("static tile shape");
    let mut keep_going = || CbcControl::Continue;
    let mut runs = 0_u32;
    loop {
        runs += 1;
        assert!(runs < 1_000_000, "allowance loop failed to converge");
        match executor
            .run(&mut keep_going, tile, one_shot_spent / 97 + 1)
            .expect("sliced runs must not breach the schedule")
        {
            CbcRunStatus::Completed => break,
            CbcRunStatus::AllowanceExhausted(boundary) => {
                assert!(
                    matches!(
                        boundary,
                        CbcBoundary::PointBlock | CbcBoundary::CandidateBlock | CbcBoundary::Prefix
                    ),
                    "in-flight exhaustion must name a real tile boundary"
                );
            }
            CbcRunStatus::Cancelled(_) => panic!("nothing requested cancellation"),
        }
    }
    assert!(runs > 1, "the sliced allowance must actually pause");
    assert_eq!(
        executor.work_spent(),
        one_shot_spent,
        "pause/resume changed the debit total"
    );
    assert_eq!(
        executor
            .into_lattice()
            .expect("completed construction yields the lattice")
            .z,
        reference.z,
        "pause/resume changed the constructed bytes"
    );
}

#[test]
fn cbx_005_zero_allowance_invokes_no_work() {
    let mut executor = admitted_executor(16, 4);
    let mut poll_count = 0_u32;
    let mut counting_poll = || {
        poll_count += 1;
        CbcControl::Continue
    };
    let status = executor
        .run(
            &mut counting_poll,
            CbcTileShape::new(2, 4).expect("static tile shape"),
            0,
        )
        .expect("a zero allowance is a status, not an error");
    assert_eq!(status, CbcRunStatus::AllowanceExhausted(CbcBoundary::Entry));
    assert_eq!(executor.work_spent(), 0, "zero allowance performed work");
    assert!(
        executor.prefix().is_empty(),
        "zero allowance committed state"
    );
}

#[test]
fn cbx_006_debits_stay_within_the_admission_estimate() {
    for (n, dim) in CASES {
        let problem = CbcProblem::new(n, dim).expect("battery cases are structurally valid");
        let estimate = problem.estimate().expect("battery cases estimate");
        let mut executor = admitted_executor(n, dim);
        run_to_completion(
            &mut executor,
            CbcTileShape::new(1, 1).expect("static tile shape"),
        );
        assert!(
            executor.work_spent() <= estimate.work_units(),
            "n={n} dim={dim}: executor debits {} exceed the admitted schedule {}",
            executor.work_spent(),
            estimate.work_units(),
        );
        assert!(
            executor.work_spent() > 0,
            "n={n} dim={dim}: a completed run must have debited work"
        );
    }
}

#[test]
fn cbx_007_structural_refusals_are_typed() {
    assert_eq!(
        CbcTileShape::new(0, 4),
        Err(CbcExecError::InvalidTileShape {
            candidate_block: 0,
            point_block: 4
        })
    );
    assert_eq!(
        CbcTileShape::new(4, 0),
        Err(CbcExecError::InvalidTileShape {
            candidate_block: 4,
            point_block: 0
        })
    );

    let mut executor = admitted_executor(5, 2);
    run_to_completion(
        &mut executor,
        CbcTileShape::new(1, 1).expect("static tile shape"),
    );
    let mut keep_going = || CbcControl::Continue;
    assert_eq!(
        executor.run(
            &mut keep_going,
            CbcTileShape::new(1, 1).expect("static tile shape"),
            u128::MAX
        ),
        Err(CbcExecError::AlreadyComplete),
        "running a completed executor must refuse"
    );

    assert_eq!(CBC_EXECUTOR_SCHEMA_VERSION, 1);
}
