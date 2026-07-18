//! G3/G5 full-study replay for production world-fork steering (7tv.21.45).
//!
//! The retained identity covers the causal seed, complete resumable state, and
//! ordered steering lineage. This fixture independently recomputes every final
//! objective, replays a trunk plus two diverged branches exactly, and proves
//! that recomputed identities refuse isolated seed and lineage mutations.
//!
//! This target does not claim convergence, optimizer superiority, all seeds or
//! configurations, cross-ISA equality, cancellation, persistence,
//! authenticated admission, or performance.

#![deny(unsafe_code)]

use fs_dfo::steer::{
    STEERED_STUDY_IDENTITY_DOMAIN, STEERED_STUDY_IDENTITY_SCHEMA_VERSION, SteeredStudy,
};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, EventKind, Severity};

const SUITE: &str = "fs-dfo/steer-study-replay";
const CASE: &str = "nested-world-fork-complete-identity";
const INPUT_SEED: u64 = 0x5EED_2145;
const DIMENSION: usize = 3;
const POPULATION: usize = 12;
const OBJECTIVES: usize = 2;
const BOUNDS: (f64, f64) = (-2.0, 2.0);
const TRUNK_GENERATIONS: u64 = 2;
const REFORK_GENERATIONS: u64 = 2;
const TAIL_GENERATIONS: u64 = 10;
const LEFT_WEIGHTS: [f64; OBJECTIVES] = [0.9, 0.1];
const LEFT_REFORK_WEIGHTS: [f64; OBJECTIVES] = [0.75, 0.25];
const RIGHT_WEIGHTS: [f64; OBJECTIVES] = [0.1, 0.9];

#[derive(Debug, Clone)]
struct StudyFrame {
    fixture: ReplayIdentity,
    parent_at_fork: ReplayIdentity,
    trunk: SteeredStudy,
    left: SteeredStudy,
    right: SteeredStudy,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct IdentityMismatch {
    expected: u64,
    found: u64,
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixture cardinality fits u64")
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

/// Two independently evaluable objectives with optima at opposite corners.
fn objectives(decision: &[f64]) -> Vec<f64> {
    let toward_positive = decision
        .iter()
        .map(|value| {
            let residual = value - 1.0;
            residual * residual
        })
        .sum();
    let toward_negative = decision
        .iter()
        .map(|value| {
            let residual = value + 1.0;
            residual * residual
        })
        .sum();
    vec![toward_positive, toward_negative]
}

/// Loop-form oracle kept separate from the iterator-form production callback.
fn objective_oracle(decision: &[f64]) -> [f64; OBJECTIVES] {
    let mut toward_positive = 0.0;
    let mut toward_negative = 0.0;
    for &value in decision {
        let positive_residual = value - 1.0;
        let negative_residual = value + 1.0;
        toward_positive += positive_residual * positive_residual;
        toward_negative += negative_residual * negative_residual;
    }
    [toward_positive, toward_negative]
}

fn fixture_identity() -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-dfo-steered-study-fixture-v1")
        .str("algorithm", "fs_dfo::SteeredStudy")
        .str("objective", "opposed-shifted-spheres-v1")
        .str("units", "dimensionless")
        .u64("dimension", usize_u64(DIMENSION))
        .u64("population", usize_u64(POPULATION))
        .u64("objectives", usize_u64(OBJECTIVES))
        .f64_bits("lower-bound", BOUNDS.0)
        .f64_bits("upper-bound", BOUNDS.1)
        .u64("trunk-generations", TRUNK_GENERATIONS)
        .u64("refork-generations", REFORK_GENERATIONS)
        .u64("tail-generations", TAIL_GENERATIONS)
        .u64("input-seed", INPUT_SEED)
        .str("study-identity-domain", STEERED_STUDY_IDENTITY_DOMAIN)
        .u64(
            "study-identity-schema-version",
            u64::from(STEERED_STUDY_IDENTITY_SCHEMA_VERSION),
        )
        .str("fs-dfo-version", fs_dfo::VERSION)
        .str("fs-obs-version", fs_obs::VERSION);
    for &weight in &LEFT_WEIGHTS {
        builder = builder.f64_bits("left-weight", weight);
    }
    for &weight in &LEFT_REFORK_WEIGHTS {
        builder = builder.f64_bits("left-refork-weight", weight);
    }
    for &weight in &RIGHT_WEIGHTS {
        builder = builder.f64_bits("right-weight", weight);
    }
    builder.finish()
}

fn run_frame() -> StudyFrame {
    let mut objective = objectives;
    let mut trunk = SteeredStudy::start(
        &mut objective,
        DIMENSION,
        BOUNDS,
        POPULATION,
        OBJECTIVES,
        INPUT_SEED,
    );
    trunk.advance(&mut objective, BOUNDS, TRUNK_GENERATIONS);
    let parent_at_fork = trunk.replay_identity();

    let mut left = trunk.fork(LEFT_WEIGHTS.to_vec());
    let mut right = trunk.fork(RIGHT_WEIGHTS.to_vec());
    assert_eq!(
        trunk.replay_identity().canonical_bytes(),
        parent_at_fork.canonical_bytes(),
        "forking must leave the parent identity unchanged"
    );

    left.advance(&mut objective, BOUNDS, REFORK_GENERATIONS);
    left = left.fork(LEFT_REFORK_WEIGHTS.to_vec());
    left.advance(&mut objective, BOUNDS, TAIL_GENERATIONS);
    right.advance(
        &mut objective,
        BOUNDS,
        REFORK_GENERATIONS + TAIL_GENERATIONS,
    );
    trunk.advance(
        &mut objective,
        BOUNDS,
        REFORK_GENERATIONS + TAIL_GENERATIONS,
    );

    StudyFrame {
        fixture: fixture_identity(),
        parent_at_fork,
        trunk,
        left,
        right,
    }
}

fn frame_identity(frame: &StudyFrame) -> ReplayIdentity {
    IdentityBuilder::new("fs-dfo-steered-study-frame-v1")
        .child("fixture", &frame.fixture)
        .child("parent-at-fork", &frame.parent_at_fork)
        .child("trunk", &frame.trunk.replay_identity())
        .child("left", &frame.left.replay_identity())
        .child("right", &frame.right.replay_identity())
        .finish()
}

fn assert_exact_study_replay(label: &str, original: &SteeredStudy, replay: &SteeredStudy) {
    assert_eq!(original.seed, replay.seed, "{label}: causal seed");
    assert_eq!(
        state_bits(original),
        state_bits(replay),
        "{label}: complete public state"
    );
    assert_eq!(original.lineage, replay.lineage, "{label}: ordered lineage");
    assert_eq!(
        original.replay_identity().canonical_bytes(),
        replay.replay_identity().canonical_bytes(),
        "{label}: complete canonical replay identity"
    );
}

fn state_bits(study: &SteeredStudy) -> Vec<u64> {
    let mut state = vec![
        study.state.stream_index,
        usize_u64(study.state.population.len()),
    ];
    for individual in &study.state.population {
        state.push(usize_u64(individual.x.len()));
        state.extend(bits(&individual.x));
        state.push(usize_u64(individual.f.len()));
        state.extend(bits(&individual.f));
    }
    state.push(usize_u64(study.state.weights.len()));
    state.extend(bits(&study.state.weights));
    state
}

fn assert_population_oracle(study: &SteeredStudy) {
    assert_eq!(study.state.population.len(), POPULATION);
    assert_eq!(study.state.weights.len(), OBJECTIVES);
    for individual in &study.state.population {
        assert_eq!(individual.x.len(), DIMENSION);
        assert!(
            individual
                .x
                .iter()
                .all(|value| value.is_finite() && (BOUNDS.0..=BOUNDS.1).contains(value))
        );
        assert_eq!(
            bits(&individual.f),
            bits(&objective_oracle(&individual.x)),
            "every retained objective must match the independent analytic oracle"
        );
    }
}

fn admit_identity(
    expected: &ReplayIdentity,
    candidate: &SteeredStudy,
) -> Result<(), IdentityMismatch> {
    let found = candidate.replay_identity();
    if found.canonical_bytes() == expected.canonical_bytes() {
        Ok(())
    } else {
        Err(IdentityMismatch {
            expected: expected.root(),
            found: found.root(),
        })
    }
}

fn next_generation_trace(study: &SteeredStudy) -> Vec<Vec<u64>> {
    let mut study = study.clone();
    let mut trace = Vec::with_capacity(POPULATION);
    let mut objective = |decision: &[f64]| {
        trace.push(bits(decision));
        objectives(decision)
    };
    study.advance(&mut objective, BOUNDS, 1);
    trace
}

fn emit_receipt(
    frame: &ReplayIdentity,
    seed_error: IdentityMismatch,
    lineage_error: IdentityMismatch,
) {
    let mut emitter = Emitter::new(SUITE, CASE);
    let receipt = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "steered-study-replay-receipt".to_string(),
            json: format!(
                "{{\"input_seed\":{INPUT_SEED},\"frame\":\"{}\",\
                 \"seed_mutant\":\"{:016x}\",\"lineage_mutant\":\"{:016x}\"}}",
                frame.hex(),
                seed_error.found,
                lineage_error.found
            ),
        },
        None,
    );
    let receipt_line = receipt.to_jsonl();
    fs_obs::validate_line(&receipt_line).expect("replay receipt must use the fs-obs wire schema");
    println!("{receipt_line}");

    let verdict = emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: CASE.to_string(),
            pass: true,
            detail: "trunk and nested forks replayed exactly; objective oracle passed; causal seed and ordered-lineage mutations were refused"
                .to_string(),
            seed: INPUT_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&verdict).expect("replay verdict must be replayable");
    let verdict_line = verdict.to_jsonl();
    fs_obs::validate_line(&verdict_line).expect("replay verdict must use the fs-obs wire schema");
    println!("{verdict_line}");
}

#[test]
fn steered_study_replays_and_refuses_seed_and_lineage_mutations() {
    let original = run_frame();
    let replay = run_frame();
    assert_eq!(
        original.fixture.canonical_bytes(),
        replay.fixture.canonical_bytes(),
        "fixture identity must replay exactly"
    );
    assert_eq!(
        original.parent_at_fork.canonical_bytes(),
        replay.parent_at_fork.canonical_bytes(),
        "parent-at-fork canonical identity must replay exactly"
    );
    assert_exact_study_replay("trunk", &original.trunk, &replay.trunk);
    assert_exact_study_replay("left", &original.left, &replay.left);
    assert_exact_study_replay("right", &original.right, &replay.right);

    let original_identity = frame_identity(&original);
    let replay_identity = frame_identity(&replay);
    assert_eq!(
        original_identity.canonical_bytes(),
        replay_identity.canonical_bytes(),
        "the complete trunk/fork frame must replay exactly"
    );

    for study in [&original.trunk, &original.left, &original.right] {
        assert_population_oracle(study);
        assert_eq!(
            study.fingerprint(),
            study.replay_identity().root(),
            "the compatibility fingerprint must project the complete identity"
        );
    }
    assert!(original.trunk.lineage.is_empty());
    assert_eq!(original.left.lineage.len(), 2);
    assert_eq!(original.right.lineage.len(), 1);
    assert_eq!(original.left.lineage[0].at_generation, TRUNK_GENERATIONS);
    assert_eq!(
        original.left.lineage[1].at_generation,
        TRUNK_GENERATIONS + REFORK_GENERATIONS
    );
    assert_eq!(bits(&original.left.lineage[0].from), bits(&[0.5, 0.5]));
    assert_eq!(bits(&original.left.lineage[0].to), bits(&LEFT_WEIGHTS));
    assert_eq!(bits(&original.left.lineage[1].from), bits(&LEFT_WEIGHTS));
    assert_eq!(
        bits(&original.left.lineage[1].to),
        bits(&LEFT_REFORK_WEIGHTS)
    );
    assert_eq!(
        bits(&original.left.state.weights),
        bits(&LEFT_REFORK_WEIGHTS)
    );
    assert_ne!(
        original.trunk.replay_identity().canonical_bytes(),
        original.left.replay_identity().canonical_bytes()
    );
    assert_ne!(
        original.left.replay_identity().canonical_bytes(),
        original.right.replay_identity().canonical_bytes()
    );
    assert_ne!(
        state_bits(&original.left),
        state_bits(&original.right),
        "opposite steering histories must produce genuinely diverged states"
    );

    let reference = original.left.replay_identity();
    let mut seed_mutant = original.left.clone();
    seed_mutant.seed ^= 1;
    assert_eq!(state_bits(&seed_mutant), state_bits(&original.left));
    assert_eq!(seed_mutant.lineage, original.left.lineage);
    let seed_error = admit_identity(&reference, &seed_mutant)
        .expect_err("changing only the causal seed must move the complete identity");
    assert_ne!(seed_error.expected, seed_error.found);
    assert_ne!(
        next_generation_trace(&seed_mutant),
        next_generation_trace(&original.left),
        "the bound seed is causal: it changes the next generated candidates"
    );

    let mut lineage_mutant = original.left.clone();
    lineage_mutant.lineage.swap(0, 1);
    assert_eq!(lineage_mutant.seed, original.left.seed);
    assert_eq!(state_bits(&lineage_mutant), state_bits(&original.left));
    let lineage_error = admit_identity(&reference, &lineage_mutant)
        .expect_err("reordering only the retained lineage must move the complete identity");
    assert_ne!(lineage_error.expected, lineage_error.found);

    emit_receipt(&original_identity, seed_error, lineage_error);
}
