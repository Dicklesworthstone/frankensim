use super::*;
use super::super::{DifferentialEvolutionStop, DifferentialEvolutionRunError};
use fs_alloc::{ArenaConfig, ArenaPool, OperationMemoryLease};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey, TilePool};
use std::convert::Infallible;
use std::sync::atomic::{AtomicUsize, Ordering};

const MODEL: [u8; 32] = [19; 32];
const LIMITS: DifferentialEvolutionRestoreLimits = DifferentialEvolutionRestoreLimits {
    max_bytes: 1 << 20, max_dimension: 3, max_population_entries: 24,
};
fn with_cx<T>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> T) -> T {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&gate, &Cx::new(
        &gate, arena, StreamKey { seed: 17, kernel_id: 55, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn state(cx: &Cx<'_>, limit: usize) -> DifferentialEvolutionState {
    DifferentialEvolutionState::new(cx, &[[-5.0, 5.0]; 3], Some(&[-0.0, 0.0, 0.0]),
        DifferentialEvolutionConfig { population_size: 8, differential_weight: 0.8,
            crossover_probability: 0.9, seed: 17, max_evaluations: limit,
            max_population_entries: 24, target: None }).unwrap()
}
fn objective(x: &[f64]) -> f64 { x.iter().map(|v| (v - 0.25).powi(2)).sum() }
fn score(x: &[f64]) -> Result<f64, Infallible> { Ok(objective(x)) }
fn save(cx: &Cx<'_>, state: &DifferentialEvolutionState) -> Vec<u8> {
    state.checkpoint_bytes(cx, MODEL, LIMITS.max_bytes).unwrap()
}
fn restore(cx: &Cx<'_>, bytes: &[u8]) -> DifferentialEvolutionState {
    DifferentialEvolutionState::from_checkpoint_bytes(cx, bytes, MODEL, LIMITS).unwrap()
}
fn rehash(mut bytes: Vec<u8>) -> Vec<u8> {
    let end = bytes.len() - DIGEST_BYTES;
    let digest = fs_blake3::hash_domain(DOMAIN, &bytes[..end]);
    bytes[end..].copy_from_slice(digest.as_bytes());
    bytes
}

#[test]
fn g5_round_trip_at_every_early_budget_preserves_exact_continuation() {
    with_cx(|_, cx| {
        for cut in 0..=33 {
            let mut original = state(cx, cut);
            original.advance(cx, usize::MAX, &mut score).unwrap();
            let bytes = save(cx, &original);
            let mut resumed = restore(cx, &bytes);
            assert_eq!(save(cx, &resumed), bytes);
            assert_eq!(resumed.evaluations(), cut);
            assert_eq!(resumed.advance(cx, 0, &mut score).unwrap().stop,
                DifferentialEvolutionStop::EvaluationLimit);
            original.extend_evaluation_budget(113).unwrap();
            resumed.extend_evaluation_budget(113).unwrap();
            let expected = original.advance(cx, usize::MAX, &mut score).unwrap();
            let mut actual = resumed.advance(cx, 7, &mut score).unwrap();
            while actual.stop == DifferentialEvolutionStop::Paused {
                actual = resumed.advance(cx, 7, &mut score).unwrap();
            }
            assert_eq!(actual, expected);
            assert_eq!(save(cx, &resumed), save(cx, &original));
        }
    });
}

#[test]
fn g4_saved_failed_wave_reuses_successful_siblings_in_a_new_search() {
    let bytes = with_cx(|_, cx| {
        let mut search = state(cx, 32);
        let pool = TilePool::for_host(4, 17);
        let lease = OperationMemoryLease::bounded(64 << 20);
        assert!(search.advance_parallel(cx, &pool, &lease, 8, 8, &|_, x| {
            if x == [0.0; 3] { Err("physical refusal") } else { Ok(objective(x)) }
        }).is_err());
        assert_eq!(search.evaluations(), 8);
        assert!(search.population().is_empty());
        save(cx, &search)
    });
    with_cx(|_, cx| {
        let mut resumed = restore(cx, &bytes);
        let calls = AtomicUsize::new(0);
        let report = resumed.advance(cx, 1, &mut |x| {
            assert_eq!(x, [0.0; 3]);
            calls.fetch_add(1, Ordering::Relaxed);
            score(x)
        }).unwrap();
        assert_eq!(calls.load(Ordering::Relaxed), 1);
        assert_eq!(report.evaluations, 9);
        assert_eq!(report.initialized_members, 8);
        let mut reference = state(cx, 8);
        reference.advance(cx, 8, &mut score).unwrap();
        assert_eq!(resumed.population(), reference.population());
        assert_eq!(resumed.best(), reference.best());
    });
}

#[test]
fn g4_snapshot_keeps_a_completed_callback_across_cancellation_and_process_boundary() {
    let search = with_cx(|gate, cx| {
        let mut search = state(cx, 32);
        let result = search.advance(cx, 1, &mut |x| {
            gate.request();
            score(x)
        });
        assert!(matches!(result, Err(DifferentialEvolutionRunError::Solver(Error::Cancelled))));
        assert_eq!(search.evaluations(), 1);
        assert!(search.population().is_empty());
        assert!(search.checkpoint_bytes(cx, MODEL, LIMITS.max_bytes).is_err());
        search
    });
    let bytes = with_cx(|_, cx| save(cx, &search));
    drop(search);
    with_cx(|_, cx| {
        let mut resumed = restore(cx, &bytes);
        let result = resumed.advance(cx, 0, &mut |_| -> Result<f64, Infallible> {
            panic!("completed physical result must not be recomputed")
        }).unwrap();
        assert_eq!(result.evaluations, 1);
        assert_eq!(result.initialized_members, 1);
        assert_eq!(resumed.population()[0].coordinates[0].to_bits(), (-0.0_f64).to_bits());
    });
}

#[test]
fn g0_corruption_model_version_bounds_and_trailing_data_refuse() {
    with_cx(|_, cx| {
        let mut search = state(cx, 32);
        search.advance(cx, 11, &mut score).unwrap();
        let bytes = save(cx, &search);
        assert!(DifferentialEvolutionState::from_checkpoint_bytes(cx, &bytes, [20; 32], LIMITS).is_err());
        for limits in [DifferentialEvolutionRestoreLimits { max_bytes: bytes.len()-1, ..LIMITS },
            DifferentialEvolutionRestoreLimits { max_dimension: 2, ..LIMITS },
            DifferentialEvolutionRestoreLimits { max_population_entries: 23, ..LIMITS }]
        { assert!(DifferentialEvolutionState::from_checkpoint_bytes(cx, &bytes, MODEL, limits).is_err()); }
        for n in [0, 31, 64, bytes.len()-1] {
            assert!(DifferentialEvolutionState::from_checkpoint_bytes(cx, &bytes[..n], MODEL, LIMITS).is_err());
        }
        for offset in [0, 8, 12, 16, 52, bytes.len()-1] {
            let mut bad = bytes.clone(); bad[offset] ^= 0x80;
            assert!(DifferentialEvolutionState::from_checkpoint_bytes(cx, &bad, MODEL, LIMITS).is_err());
        }
        let mut version = bytes.clone(); version[8..12].copy_from_slice(&99_u32.to_le_bytes());
        assert!(DifferentialEvolutionState::from_checkpoint_bytes(cx, &rehash(version), MODEL, LIMITS).is_err());
        let mut trailing = bytes[..bytes.len()-DIGEST_BYTES].to_vec();
        trailing.push(7); trailing.extend_from_slice(&[0; DIGEST_BYTES]);
        assert!(DifferentialEvolutionState::from_checkpoint_bytes(cx, &rehash(trailing), MODEL, LIMITS).is_err());
        assert!(search.checkpoint_bytes(cx, MODEL, bytes.len()-1).is_err());
        assert_eq!(save(cx, &search), bytes);
    });
}

#[test]
fn g0_inconsistent_counters_or_foreign_cached_candidates_do_not_restore() {
    with_cx(|_, cx| {
        let mut search = state(cx, 32);
        search.advance(cx, 9, &mut score).unwrap();
        search.evaluations = 0;
        assert!(DifferentialEvolutionState::from_checkpoint_bytes(cx, &save(cx, &search), MODEL, LIMITS).is_err());
        search.evaluations = 9;
        search.pending = Some(Pending { coordinates: vec![4.0; 3], value: None });
        assert!(DifferentialEvolutionState::from_checkpoint_bytes(cx, &save(cx, &search), MODEL, LIMITS).is_err());
    });
}
