use super::*;
use std::{cell::Cell, convert::Infallible};

fn config(n: usize, population: usize, budget: usize) -> DifferentialEvolutionConfig {
    DifferentialEvolutionConfig { population_size: population, differential_weight: 0.8,
        crossover_probability: 0.9, seed: 17, max_evaluations: budget,
        max_population_entries: n*population, target: None }
}
fn state(bounds: &[[f64; 2]], initial: Option<&[f64]>, cfg: DifferentialEvolutionConfig)
    -> DifferentialEvolutionState
{ DifferentialEvolutionState::new_with(bounds, initial, cfg, &mut || true).unwrap() }
fn sphere(x: &[f64]) -> f64 { x.iter().map(|x| x*x).sum() }
fn run(state: &mut DifferentialEvolutionState, count: usize, trace: &mut Vec<Vec<f64>>) -> DifferentialEvolutionReport {
    state.advance_with(count, &mut |x| { trace.push(x.to_vec()); Ok::<_, Infallible>(sphere(x)) }, &mut || true).unwrap()
}

#[test]
fn g0_every_budget_boundary_and_pause_partition_preserves_candidate_trace() {
    let bounds = [[-5.0, 5.0]; 3];
    let mut reference = state(&bounds, None, config(3, 8, 128));
    let mut expected = Vec::new();
    let report = run(&mut reference, 128, &mut expected);
    assert_eq!(report.stop, DifferentialEvolutionStop::EvaluationLimit);
    for width in [1, 3, 7, 17, 127] {
        let mut replay = state(&bounds, None, config(3, 8, 128));
        let mut trace = Vec::new();
        while replay.evaluations() < 128 { run(&mut replay, width, &mut trace); }
        assert_eq!(trace, expected);
        assert_eq!(replay.population(), reference.population());
        assert_eq!(replay.best(), reference.best());
    }
    let mut replay = state(&bounds, None, config(3, 8, 0));
    let mut trace = Vec::new();
    assert!(run(&mut replay, 100, &mut trace).best.is_none());
    for limit in 1..=128 {
        replay.extend_evaluation_budget(limit).unwrap();
        let report = run(&mut replay, 1000, &mut trace);
        assert_eq!(report.evaluations, limit);
        assert_eq!(trace, expected[..limit]);
    }
    assert_eq!(replay.population(), reference.population());
    assert!(replay.extend_evaluation_budget(128).is_err());
}

#[test]
fn g0_donors_are_distinct_and_exclude_target_even_for_minimum_population() {
    for count in 4..24 {
        for target in 0..count {
            for generation in 0..8 {
                let mut rng = CandidateRandom { key: StreamKey { seed: 17,
                    kernel: DIFFERENTIAL_EVOLUTION_STREAM_KERNEL, tile: target as u32 },
                    next: generation*512, end: (generation+1)*512 };
                let d = donors(&mut rng, count, target, &mut || true).unwrap();
                assert!(d.iter().all(|&i| i<count && i!=target));
                assert!(d[0]!=d[1] && d[0]!=d[2] && d[1]!=d[2]);
            }
        }
    }
    let mut rng = CandidateRandom { key: StreamKey { seed: 0, kernel: 0, tile: 0 }, next: 9, end: 9 };
    assert_eq!(rng.uniform(), Err(DifferentialEvolutionError::RandomBudget));
}

#[test]
fn g3_frozen_generations_publish_atomically_and_neutral_trials_keep_exploring() {
    let mut de = state(&[[0.0, 1.0]; 2], Some(&[0.0, 0.0]), config(2, 8, 24));
    let mut flat = |_: &[f64]| Ok::<_, Infallible>(1.0);
    de.advance_with(8, &mut flat, &mut || true).unwrap();
    let parents = de.population().to_vec();
    for completed in 1..8 {
        let report = de.advance_with(1, &mut flat, &mut || true).unwrap();
        assert_eq!(de.population(), parents);
        assert_eq!(report.pending_members, completed);
        assert_eq!(report.generations, 0);
    }
    let report = de.advance_with(1, &mut flat, &mut || true).unwrap();
    assert_eq!(report.generations, 1);
    assert_eq!(report.pending_members, 0);
    assert_ne!(de.population(), parents);
    assert_eq!(de.best().unwrap().coordinates, [0.0, 0.0]);
    let report = de.advance_with(8, &mut flat, &mut || true).unwrap();
    assert_eq!(report.stop, DifferentialEvolutionStop::EvaluationLimit);
    assert_eq!(report.best.unwrap().value, 1.0);
}

#[test]
fn g4_completed_callback_survives_cancellation_without_repeating_physics() {
    let permit = Cell::new(true);
    let calls = Cell::new(0);
    let mut de = state(&[[0.0, 1.0]], Some(&[0.0]), config(1, 4, 20));
    let result = de.advance_with(1, &mut |x| {
        calls.set(calls.get()+1);
        permit.set(false);
        Ok::<_, Infallible>(sphere(x))
    }, &mut || permit.get());
    assert!(matches!(result, Err(DifferentialEvolutionRunError::Solver(DifferentialEvolutionError::Cancelled))));
    assert!(de.population().is_empty());
    assert_eq!(de.evaluations(), 1);
    permit.set(true);
    let report = de.advance_with(0, &mut |_| -> Result<f64, Infallible> {
        panic!("completed callback must not be repeated")
    }, &mut || permit.get()).unwrap();
    assert_eq!(report.initialized_members, 1);
    assert_eq!(report.best.unwrap().value, 0.0);
    assert_eq!(calls.get(), 1);
}

#[test]
fn g0_failed_and_nonfinite_callbacks_spend_attempts_without_fabricated_scores() {
    let mut de = state(&[[-2.0, 2.0]], None, config(1, 4, 10));
    let mut attempted = Vec::new();
    let result = de.advance_with(10, &mut |x| { attempted.push(x.to_vec()); Err::<f64,_>("physical refusal") }, &mut || true);
    assert!(matches!(result, Err(DifferentialEvolutionRunError::Evaluation("physical refusal"))));
    assert_eq!(de.evaluations(), 1);
    assert!(de.best().is_none());
    let result = de.advance_with(10, &mut |x| { attempted.push(x.to_vec()); Ok::<_, Infallible>(f64::NAN) }, &mut || true);
    assert!(matches!(result, Err(DifferentialEvolutionRunError::Solver(DifferentialEvolutionError::NonFinite(_)))));
    assert_eq!(de.evaluations(), 2);
    let report = de.advance_with(1, &mut |x| { attempted.push(x.to_vec()); Ok::<_, Infallible>(sphere(x)) }, &mut || true).unwrap();
    assert_eq!(report.evaluations, 3);
    assert_eq!(report.initialized_members, 1);
    assert!(attempted.windows(2).all(|p| p[0]==p[1]));
}

#[test]
fn g1_population_search_reaches_multimodal_and_flat_region_targets() {
    let mut cfg = config(2, 40, 20_000);
    cfg.target = Some(1e-7);
    let mut de = state(&[[-5.0, 5.0]; 2], Some(&[0.5, 0.5]), cfg);
    let result = de.advance_with(20_000, &mut |x| Ok::<_, Infallible>(x.iter()
        .map(|x| x*x + 10.0*(1.0-(2.0*std::f64::consts::PI*x).cos())).sum()), &mut || true).unwrap();
    assert_eq!(result.stop, DifferentialEvolutionStop::TargetReached);
    assert!(result.best.unwrap().value <= 1e-7);
    let mut cfg = config(1, 12, 3000);
    cfg.target = Some(1e-12);
    let mut de = state(&[[0.0, 1.0]], Some(&[0.0]), cfg);
    let result = de.advance_with(3000, &mut |x| Ok::<_, Infallible>(((x[0]-0.65).max(0.0)-0.15).powi(2)), &mut || true).unwrap();
    assert_eq!(result.stop, DifferentialEvolutionStop::TargetReached);
    assert!((result.best.unwrap().coordinates[0]-0.8).abs() <= 1e-6);
}

#[test]
fn g0_bounds_forced_crossover_and_invalid_admission_are_explicit() {
    let bounds = [[-3.0, -1.0], [0.0, 1e-200], [1e100, 2e100]];
    let mut cfg = config(3, 8, 80); cfg.crossover_probability = 0.0; cfg.differential_weight = 2.0;
    let mut de = state(&bounds, None, cfg.clone());
    de.advance_with(8, &mut |_| Ok::<_, Infallible>(1.0), &mut || true).unwrap();
    let parents = de.population().to_vec();
    let mut trials = Vec::new();
    de.advance_with(8, &mut |x| {
        assert!(x.iter().zip(bounds).all(|(&v,[lo,hi])| v.is_finite() && v>=lo && v<=hi));
        trials.push(x.to_vec()); Ok::<_, Infallible>(1.0)
    }, &mut || true).unwrap();
    assert!(trials.iter().zip(&parents).all(|(x,p)| x.iter().zip(&p.coordinates).filter(|(a,b)| a!=b).count()==1));
    for bounds in [vec![], vec![[1.0,1.0]], vec![[f64::NEG_INFINITY,1.0]], vec![[-f64::MAX,f64::MAX]]] {
        assert!(DifferentialEvolutionState::new_with(&bounds, None, cfg.clone(), &mut || true).is_err());
    }
    for weight in [0.0, -1.0, 2.1, f64::NAN] {
        cfg.differential_weight=weight;
        assert!(DifferentialEvolutionState::new_with(&[[0.0,1.0]], None, cfg.clone(), &mut || true).is_err());
    }
    assert!(DifferentialEvolutionState::new_with(&[[0.0,1.0]], Some(&[2.0]), config(1,4,10), &mut || true).is_err());
    assert!(DifferentialEvolutionState::new_with(&[[0.0,1.0]], None, config(1,3,10), &mut || true).is_err());
    let mut small=config(1,4,10);small.max_population_entries=3;
    assert!(DifferentialEvolutionState::new_with(&[[0.0,1.0]], None, small, &mut || true).is_err());
    assert!(matches!(DifferentialEvolutionState::new_with(&[[0.0,1.0]], None, config(1,4,10), &mut || false), Err(DifferentialEvolutionError::Cancelled)));
}
