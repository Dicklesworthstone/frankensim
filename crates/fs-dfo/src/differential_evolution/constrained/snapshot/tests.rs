use super::*;
use super::super::{ConstrainedEvaluation, DifferentialEvolutionConfig, DifferentialEvolutionRunError};
use super::super::super::{DifferentialEvolutionState, DifferentialEvolutionStop};
use fs_alloc::{ArenaConfig, ArenaPool, OperationMemoryLease};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey, TilePool};
use std::convert::Infallible;
use std::sync::atomic::{AtomicUsize, Ordering};

const MODEL: [u8; 32] = [49; 32];
const LIMITS: DifferentialEvolutionRestoreLimits = DifferentialEvolutionRestoreLimits {
    max_bytes: 1 << 20, max_dimension: 2, max_population_entries: 16,
};
fn config(limit: usize) -> DifferentialEvolutionConfig {
    DifferentialEvolutionConfig { population_size: 8, differential_weight: 0.8,
        crossover_probability: 0.9, seed: 17, max_evaluations: limit,
        max_population_entries: 16, target: None }
}
fn with_cx<T>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> T) -> T {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&gate, &Cx::new(
        &gate, arena, StreamKey { seed: 17, kernel_id: 55, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn state(cx: &Cx<'_>, limit: usize) -> ConstrainedDifferentialEvolutionState {
    ConstrainedDifferentialEvolutionState::new(cx, &[[0.0, 1.0]; 2], Some(&[0.0; 2]),
        config(limit), 2, 0.0).unwrap()
}
fn evaluation(x: &[f64]) -> Result<ConstrainedEvaluation, Infallible> {
    Ok(ConstrainedEvaluation { objective: x[0]*x[0] + x[1]*x[1],
        inequalities: vec![0.6-x[0], x[0]+x[1]-1.2] })
}
fn save(cx: &Cx<'_>, state: &ConstrainedDifferentialEvolutionState) -> Vec<u8> {
    state.checkpoint_bytes(cx, MODEL, LIMITS.max_bytes).unwrap()
}
fn restore(cx: &Cx<'_>, bytes: &[u8]) -> ConstrainedDifferentialEvolutionState {
    ConstrainedDifferentialEvolutionState::from_checkpoint_bytes(cx, bytes, MODEL, 2, 0.0, LIMITS).unwrap()
}
fn rehash(mut bytes: Vec<u8>) -> Vec<u8> {
    let end = bytes.len()-32;
    let digest = fs_blake3::hash_domain(DOMAIN, &bytes[..end]);
    bytes[end..].copy_from_slice(digest.as_bytes());
    bytes
}

#[test]
fn g5_every_early_checkpoint_preserves_constrained_continuation_and_budget() {
    with_cx(|_, cx| {
        for cut in 0..=33 {
            let mut original = state(cx, cut);
            original.advance(cx, usize::MAX, &mut evaluation).unwrap();
            let bytes = save(cx, &original);
            let mut resumed = restore(cx, &bytes);
            assert_eq!(save(cx, &resumed), bytes);
            assert_eq!(resumed.evaluations(), cut);
            let exhausted = resumed.advance(cx, 1, &mut |_| -> Result<ConstrainedEvaluation, Infallible> {
                panic!("restore replenished the spent budget")
            }).unwrap();
            assert_eq!(exhausted.stop, DifferentialEvolutionStop::EvaluationLimit);
            original.extend_evaluation_budget(129).unwrap();
            resumed.extend_evaluation_budget(129).unwrap();
            let expected = original.advance(cx, usize::MAX, &mut evaluation).unwrap();
            let mut actual = resumed.advance(cx, 7, &mut evaluation).unwrap();
            while actual.stop == DifferentialEvolutionStop::Paused {
                actual = resumed.advance(cx, 7, &mut evaluation).unwrap();
            }
            assert_eq!(actual, expected);
            assert_eq!(save(cx, &resumed), save(cx, &original));
        }
    });
}

#[test]
fn g4_restart_retries_one_failed_parallel_member_not_seven_successful_siblings() {
    let bytes = with_cx(|_, cx| {
        let mut search = state(cx, 32);
        let pool = TilePool::for_host(4, 17);
        let lease = OperationMemoryLease::bounded(64 << 20);
        assert!(search.advance_parallel(cx, &pool, &lease, 8, 8, &|_, x| {
            if x == [0.0; 2] { Err("physical refusal") }
            else { Ok(evaluation(x).unwrap()) }
        }).is_err());
        assert_eq!(search.evaluations(), 8);
        assert!(search.population().is_empty());
        save(cx, &search)
    });
    with_cx(|_, cx| {
        let mut resumed = restore(cx, &bytes);
        let calls = AtomicUsize::new(0);
        let report = resumed.advance(cx, 1, &mut |x| {
            assert_eq!(x, [0.0; 2]);
            calls.fetch_add(1, Ordering::Relaxed);
            evaluation(x)
        }).unwrap();
        assert_eq!(calls.load(Ordering::Relaxed), 1);
        assert_eq!(report.evaluations, 9);
        assert_eq!(report.initialized_members, 8);
        let mut reference = state(cx, 8);
        reference.advance(cx, 8, &mut evaluation).unwrap();
        assert_eq!(resumed.population(), reference.population());
        assert_eq!(resumed.best(), reference.best());
    });
}

#[test]
fn g4_cancelled_completed_callback_survives_drop_and_restore() {
    let search = with_cx(|gate, cx| {
        let mut search = state(cx, 32);
        let result = search.advance(cx, 1, &mut |x| {
            gate.request();
            evaluation(x)
        });
        assert!(matches!(result, Err(DifferentialEvolutionRunError::Solver(Error::Cancelled))));
        assert_eq!(search.evaluations(), 1);
        assert!(search.best().is_none());
        assert!(search.checkpoint_bytes(cx, MODEL, LIMITS.max_bytes).is_err());
        search
    });
    let bytes = with_cx(|_, cx| save(cx, &search));
    drop(search);
    with_cx(|_, cx| {
        let mut resumed = restore(cx, &bytes);
        let report = resumed.advance(cx, 0, &mut |_| -> Result<ConstrainedEvaluation, Infallible> {
            panic!("completed physical callback repeated after restart")
        }).unwrap();
        assert_eq!(report.evaluations, 1);
        assert_eq!(report.initialized_members, 1);
        assert!(!report.best.unwrap().value.is_feasible());
    });
}

#[test]
fn g0_changed_model_policy_malformed_score_and_bounded_input_refuse() {
    with_cx(|gate, cx| {
        let mut search = state(cx, 32);
        search.advance(cx, 11, &mut evaluation).unwrap();
        let bytes = save(cx, &search);
        for (model, count, tolerance) in [([50;32],2,0.0), (MODEL,3,0.0),
            (MODEL,2,1e-6), (MODEL,0,0.0), (MODEL,2,f64::NAN)] {
            assert!(ConstrainedDifferentialEvolutionState::from_checkpoint_bytes(
                cx, &bytes, model, count, tolerance, LIMITS).is_err());
        }
        for limits in [DifferentialEvolutionRestoreLimits { max_bytes: bytes.len()-1, ..LIMITS },
            DifferentialEvolutionRestoreLimits { max_dimension: 1, ..LIMITS },
            DifferentialEvolutionRestoreLimits { max_population_entries: 15, ..LIMITS }] {
            assert!(ConstrainedDifferentialEvolutionState::from_checkpoint_bytes(
                cx, &bytes, MODEL, 2, 0.0, limits).is_err());
        }
        for n in [0, 31, 64, bytes.len()-1] {
            assert!(ConstrainedDifferentialEvolutionState::from_checkpoint_bytes(
                cx, &bytes[..n], MODEL, 2, 0.0, LIMITS).is_err());
        }
        let end = bytes.len()-32;
        for invalid_violation in [-1.0_f64, f64::INFINITY, f64::NAN] {
            let mut bad = bytes.clone();
            bad[end-8..end].copy_from_slice(&invalid_violation.to_bits().to_le_bytes());
            assert!(ConstrainedDifferentialEvolutionState::from_checkpoint_bytes(
                cx, &rehash(bad), MODEL, 2, 0.0, LIMITS).is_err());
        }
        let mut bad = bytes.clone(); bad[8..12].copy_from_slice(&99_u32.to_le_bytes());
        assert!(ConstrainedDifferentialEvolutionState::from_checkpoint_bytes(
            cx, &rehash(bad), MODEL, 2, 0.0, LIMITS).is_err());
        let mut trailing = bytes[..end].to_vec(); trailing.push(7); trailing.extend_from_slice(&[0;32]);
        assert!(ConstrainedDifferentialEvolutionState::from_checkpoint_bytes(
            cx, &rehash(trailing), MODEL, 2, 0.0, LIMITS).is_err());
        let mut corrupt = bytes.clone(); corrupt[end] ^= 1;
        assert!(ConstrainedDifferentialEvolutionState::from_checkpoint_bytes(
            cx, &corrupt, MODEL, 2, 0.0, LIMITS).is_err());
        assert!(search.checkpoint_bytes(cx, MODEL, bytes.len()-1).is_err());
        gate.request();
        assert!(ConstrainedDifferentialEvolutionState::from_checkpoint_bytes(
            cx, &bytes, MODEL, 2, 0.0, LIMITS).is_err());
    });
}

#[test]
fn g3_feasible_ties_preserve_the_earliest_best_even_after_its_parent_is_replaced() {
    with_cx(|_, cx| {
        let mut search = ConstrainedDifferentialEvolutionState::new(cx,
            &[[0.0, 1.0];2], Some(&[0.0;2]), config(32), 2, 0.5).unwrap();
        search.advance(cx, 8, &mut |_| Ok::<_, Infallible>(ConstrainedEvaluation {
            objective: 1.0, inequalities: vec![0.1, 0.0] })).unwrap();
        let earliest = search.best().unwrap().clone();
        search.advance(cx, 8, &mut |_| Ok::<_, Infallible>(ConstrainedEvaluation {
            objective: 1.0, inequalities: vec![0.4, 0.0] })).unwrap();
        assert!(search.population().iter().all(|p| p.value.max_scaled_violation() == 0.4));
        assert_eq!(search.best(), Some(&earliest));
        let bytes = save(cx, &search);
        let resumed = ConstrainedDifferentialEvolutionState::from_checkpoint_bytes(
            cx, &bytes, MODEL, 2, 0.5, LIMITS).unwrap();
        assert_eq!(resumed.best(), Some(&earliest));
        assert_eq!(save(cx, &resumed), bytes);
    });
}

#[test]
fn g0_restored_infeasible_objective_target_is_not_reported_as_reached() {
    with_cx(|_, cx| {
        let mut cfg = config(9); cfg.target = Some(0.0);
        let mut search = ConstrainedDifferentialEvolutionState::new(cx,
            &[[0.0, 1.0];2], None, cfg, 2, 0.0).unwrap();
        let mut evaluate = |_: &[f64]| Ok::<_, Infallible>(ConstrainedEvaluation {
            objective: -1.0, inequalities: vec![1.0, -1.0] });
        search.advance(cx, 3, &mut evaluate).unwrap();
        let mut resumed = restore(cx, &save(cx, &search));
        let report = resumed.advance(cx, usize::MAX, &mut evaluate).unwrap();
        assert_eq!(report.stop, DifferentialEvolutionStop::EvaluationLimit);
        assert!(!report.best.unwrap().value.is_feasible());
    });
}

#[test]
fn g5_scalar_v1_wire_format_is_unchanged_and_score_formats_cannot_be_confused() {
    with_cx(|_, cx| {
        // Independently encoded from the pre-refactor v1 field grammar. No
        // callback or shared encoder builds the expected payload.
        let hex = concat!(
            "4653444543503031010000000100000001de0000313131313131313131313131",
            "3131313131313131313131313131313131313131020000000000000008000000",
            "000000009a9999999999e93fcdccccccccccec3f110000000000000028000000",
            "0000000010000000000000000000000000000000000000000000000000000000",
            "0000000000000000000000f03f0000000000000000000000000000f03f010000",
            "0000000000800000000000000000000000000000000000000000000000000000",
            "0000000000000000"
        );
        let mut expected: Vec<u8> = (0..hex.len()).step_by(2)
            .map(|i| u8::from_str_radix(&hex[i..i+2], 16).unwrap()).collect();
        let digest = fs_blake3::hash_domain("fs-dfo/differential-evolution-checkpoint-v1", &expected);
        expected.extend_from_slice(digest.as_bytes());
        let scalar = DifferentialEvolutionState::new(cx, &[[0.0, 1.0];2],
            Some(&[-0.0, 0.0]), config(40)).unwrap();
        let scalar_bytes = scalar.checkpoint_bytes(cx, MODEL, LIMITS.max_bytes).unwrap();
        assert_eq!(scalar_bytes, expected);
        assert!(ConstrainedDifferentialEvolutionState::from_checkpoint_bytes(
            cx, &scalar_bytes, MODEL, 2, 0.0, LIMITS).is_err());
        let constrained_bytes = save(cx, &state(cx, 40));
        assert!(DifferentialEvolutionState::from_checkpoint_bytes(
            cx, &constrained_bytes, MODEL, LIMITS).is_err());
    });
}
