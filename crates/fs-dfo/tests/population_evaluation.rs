use fs_alloc::OperationMemoryLease;
use fs_dfo::{
    PopulationCandidate, PopulationEvaluator, PopulationLimits, PopulationProvenance,
    PopulationPublishError, PopulationPublisher, PopulationRefusal,
};
use fs_exec::{Budget, CancelGate, PoolConfig, RunId, TilePool};

const POOL_SEED: u64 = 0x7A24_1600;

fn candidates() -> Vec<PopulationCandidate> {
    vec![
        PopulationCandidate {
            identity: 40,
            decision: vec![4.0],
        },
        PopulationCandidate {
            identity: 10,
            decision: vec![1.0],
        },
        PopulationCandidate {
            identity: 30,
            decision: vec![3.0],
        },
        PopulationCandidate {
            identity: 20,
            decision: vec![2.0],
        },
        PopulationCandidate {
            identity: 50,
            decision: vec![5.0],
        },
    ]
}

fn provenance(generation: u64, run: u64, population: usize, width: usize) -> PopulationProvenance {
    PopulationProvenance {
        schema_version: 1,
        generation,
        run: RunId(run),
        individuals: population,
        objective_dimension: 2,
        tiles: if population == 0 {
            0
        } else {
            population.div_ceil(width) as u64
        },
    }
}

fn limits() -> PopulationLimits {
    PopulationLimits {
        max_individuals: 8,
        max_work_units: 8,
        max_output_bytes: 8 * 24,
    }
}

fn budget(work: u64) -> Budget {
    Budget::new().with_cost_quota(work).with_poll_quota(64)
}

fn lease() -> OperationMemoryLease {
    OperationMemoryLease::bounded(1 << 20)
}

#[test]
fn tiles_preserve_semantic_order_and_replay_across_worker_counts() {
    let population = candidates();
    let run = provenance(0, 9, population.len(), 2);
    let objective = |x: &[f64]| vec![x[0], x[0] * x[0]];
    let one_lease = lease();
    let one = PopulationEvaluator::new(&TilePool::new(PoolConfig::for_host(1, POOL_SEED)), 2)
        .evaluate(
            &population,
            2,
            limits(),
            run,
            &CancelGate::new(),
            budget(population.len() as u64),
            &one_lease,
            objective,
        )
        .expect("one worker evaluates");
    let many_lease = lease();
    let many = PopulationEvaluator::new(&TilePool::new(PoolConfig::for_host(4, POOL_SEED)), 2)
        .evaluate(
            &population,
            2,
            limits(),
            run,
            &CancelGate::new(),
            budget(population.len() as u64),
            &many_lease,
            objective,
        )
        .expect("many workers evaluates");
    assert_eq!(one.generation(), many.generation());
    assert_eq!(one.provenance(), many.provenance());
    assert_eq!(
        one.evaluations().collect::<Vec<_>>(),
        many.evaluations().collect::<Vec<_>>(),
        "G5: tile placement cannot alter semantic population rows"
    );
    assert_ne!(
        one.identity_root(),
        many.identity_root(),
        "exact generation identity is bound to the placement-specific completion receipt"
    );
    assert_eq!(
        one.evaluations()
            .map(|row| row.identity())
            .collect::<Vec<_>>(),
        vec![40, 10, 30, 20, 50]
    );
    assert_eq!(
        one.evaluations()
            .map(|row| row.stream_key())
            .collect::<Vec<_>>(),
        many.evaluations()
            .map(|row| row.stream_key())
            .collect::<Vec<_>>(),
        "G5: actual TilePool streams, not worker placement, bind provenance"
    );
}

#[test]
fn mid_tile_cancellation_is_reported_as_executor_error() {
    let population = candidates();
    let gate = CancelGate::new_clock_free();
    let evaluation_lease = lease();
    let calls = std::sync::atomic::AtomicUsize::new(0);
    let result = PopulationEvaluator::new(&TilePool::new(PoolConfig::for_host(1, 3)), 2).evaluate(
        &population,
        2,
        limits(),
        provenance(0, 3, population.len(), 2),
        &gate,
        budget(population.len() as u64),
        &evaluation_lease,
        |x| {
            if calls.fetch_add(1, std::sync::atomic::Ordering::Relaxed) == 0 {
                gate.request();
            }
            vec![x[0], x[0]]
        },
    );
    assert!(
        matches!(result, Err(PopulationPublishError::Executor(_))),
        "G4 cancellation is an executor outcome"
    );
    assert_eq!(calls.load(std::sync::atomic::Ordering::Relaxed), 1);
}

#[test]
fn preflight_refuses_duplicate_identity_before_callback_or_publication() {
    let mut population = candidates();
    population[4].identity = 40;
    let calls = std::sync::atomic::AtomicUsize::new(0);
    let evaluation_lease = lease();
    let result = PopulationEvaluator::new(&TilePool::new(PoolConfig::for_host(2, 4)), 2).evaluate(
        &population,
        2,
        limits(),
        provenance(0, 4, population.len(), 2),
        &CancelGate::new_clock_free(),
        budget(population.len() as u64),
        &evaluation_lease,
        |_| {
            calls.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            vec![0.0, 0.0]
        },
    );
    assert!(matches!(
        result,
        Err(PopulationPublishError::Refused(
            PopulationRefusal::DuplicateIdentity { identity: 40 }
        ))
    ));
    assert_eq!(calls.load(std::sync::atomic::Ordering::Relaxed), 0);
}

#[test]
fn checkpoint_resume_publishes_each_generation_once() {
    let population = candidates();
    let pool = TilePool::new(PoolConfig::for_host(2, 5));
    let evaluator = PopulationEvaluator::new(&pool, 2);
    let publisher = PopulationPublisher::new(10);
    let first_lease = lease();
    let first = evaluator
        .evaluate(
            &population,
            2,
            limits(),
            provenance(0, 5, population.len(), 2),
            &CancelGate::new(),
            budget(population.len() as u64),
            &first_lease,
            |x| vec![x[0], -x[0]],
        )
        .unwrap();
    publisher.publish(first).unwrap();
    let resumed = PopulationPublisher::from_checkpoint(publisher.checkpoint(), 10).unwrap();
    let second_lease = lease();
    let second = evaluator
        .evaluate(
            &population,
            2,
            limits(),
            provenance(1, 6, population.len(), 2),
            &CancelGate::new(),
            budget(population.len() as u64),
            &second_lease,
            |x| vec![x[0], -x[0]],
        )
        .unwrap();
    resumed.publish(second).unwrap();
    assert_eq!(resumed.checkpoint().committed_work_units(), 10);
    assert_eq!(resumed.checkpoint().committed().unwrap().generation(), 1);
}

#[test]
fn publication_uses_sealed_completion_not_a_late_live_gate() {
    let population = candidates();
    let pool = TilePool::new(PoolConfig::for_host(1, 6));
    let evaluator = PopulationEvaluator::new(&pool, 2);
    let gate = CancelGate::new_clock_free();
    let evaluation_lease = lease();
    let generation = evaluator
        .evaluate(
            &population,
            2,
            limits(),
            provenance(0, 6, population.len(), 2),
            &gate,
            budget(population.len() as u64),
            &evaluation_lease,
            |x| vec![x[0], -x[0]],
        )
        .unwrap();
    let publisher = PopulationPublisher::new(population.len() as u64);
    let generation_identity = generation.identity_root();
    assert_eq!(generation.completion_witness().verify(), Ok(()));
    assert!(!generation.completion_witness().cancellation_requested());
    gate.request();
    publisher.publish(generation).unwrap();
    assert_eq!(
        publisher.checkpoint().committed_identity_root(),
        Some(generation_identity),
        "publication is bound to executor completion, not a later mutable gate"
    );
}

#[test]
fn explicit_cost_budget_refuses_before_callback_or_lease_admission() {
    let population = candidates();
    let calls = std::sync::atomic::AtomicUsize::new(0);
    let evaluation_lease = OperationMemoryLease::bounded(0);
    let result = PopulationEvaluator::new(&TilePool::new(PoolConfig::for_host(1, 7)), 2).evaluate(
        &population,
        2,
        limits(),
        provenance(0, 7, population.len(), 2),
        &CancelGate::new_clock_free(),
        budget((population.len() - 1) as u64),
        &evaluation_lease,
        |_| {
            calls.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            vec![0.0, 0.0]
        },
    );
    assert!(matches!(
        result,
        Err(PopulationPublishError::Refused(
            PopulationRefusal::BudgetWorkLimit {
                requested: 5,
                maximum: 4,
            }
        ))
    ));
    assert_eq!(calls.load(std::sync::atomic::Ordering::Relaxed), 0);
    assert_eq!(evaluation_lease.receipt().requested_bytes, 0);
}

#[test]
fn tampered_provenance_is_refused_before_callback() {
    let population = candidates();
    let calls = std::sync::atomic::AtomicUsize::new(0);
    let evaluation_lease = lease();
    let mut forged = provenance(0, 8, population.len(), 2);
    forged.tiles += 1;
    let result = PopulationEvaluator::new(&TilePool::new(PoolConfig::for_host(1, 8)), 2).evaluate(
        &population,
        2,
        limits(),
        forged,
        &CancelGate::new_clock_free(),
        budget(population.len() as u64),
        &evaluation_lease,
        |_| {
            calls.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            vec![0.0, 0.0]
        },
    );
    assert!(matches!(
        result,
        Err(PopulationPublishError::Refused(
            PopulationRefusal::ProvenanceMismatch
        ))
    ));
    assert_eq!(calls.load(std::sync::atomic::Ordering::Relaxed), 0);
}

#[test]
fn wrong_dimension_and_nonfinite_callbacks_never_mint_a_generation() {
    let population = candidates();
    let pool = TilePool::new(PoolConfig::for_host(1, 11));
    let evaluator = PopulationEvaluator::new(&pool, 2);
    let wrong_dimension_lease = lease();
    assert!(matches!(
        evaluator.evaluate(
            &population,
            2,
            limits(),
            provenance(0, 11, population.len(), 2),
            &CancelGate::new_clock_free(),
            budget(population.len() as u64),
            &wrong_dimension_lease,
            |_| vec![0.0],
        ),
        Err(PopulationPublishError::ObjectiveDimension {
            expected: 2,
            actual: 1,
            ..
        })
    ));
    let nonfinite_lease = lease();
    assert!(matches!(
        evaluator.evaluate(
            &population,
            2,
            limits(),
            provenance(0, 12, population.len(), 2),
            &CancelGate::new_clock_free(),
            budget(population.len() as u64),
            &nonfinite_lease,
            |_| vec![f64::NAN, 0.0],
        ),
        Err(PopulationPublishError::NonFiniteObjective { .. })
    ));
}

#[test]
fn output_limit_refuses_before_callback_or_lease_admission() {
    let population = candidates();
    let calls = std::sync::atomic::AtomicUsize::new(0);
    let zero_lease = OperationMemoryLease::bounded(0);
    let result = PopulationEvaluator::new(&TilePool::new(PoolConfig::for_host(1, 13)), 2).evaluate(
        &population,
        2,
        PopulationLimits {
            max_output_bytes: 1,
            ..limits()
        },
        provenance(0, 13, population.len(), 2),
        &CancelGate::new_clock_free(),
        budget(population.len() as u64),
        &zero_lease,
        |_| {
            calls.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            vec![0.0, 0.0]
        },
    );
    assert!(matches!(
        result,
        Err(PopulationPublishError::Refused(
            PopulationRefusal::OutputLimit {
                requested: 120,
                maximum: 1,
            }
        ))
    ));
    assert_eq!(calls.load(std::sync::atomic::Ordering::Relaxed), 0);
    assert_eq!(zero_lease.receipt().requested_bytes, 0);
}

#[test]
fn cumulative_publish_and_resume_caps_reject_excess_work() {
    let population = candidates();
    let pool = TilePool::new(PoolConfig::for_host(1, 14));
    let evaluator = PopulationEvaluator::new(&pool, 2);
    let first_lease = lease();
    let first = evaluator
        .evaluate(
            &population,
            2,
            limits(),
            provenance(0, 14, population.len(), 2),
            &CancelGate::new_clock_free(),
            budget(population.len() as u64),
            &first_lease,
            |x| vec![x[0], -x[0]],
        )
        .unwrap();
    let publisher = PopulationPublisher::new(population.len() as u64);
    publisher.publish(first).unwrap();
    assert!(matches!(
        PopulationPublisher::from_checkpoint(publisher.checkpoint(), 4),
        Err(PopulationPublishError::CumulativeWorkLimit {
            requested: 5,
            maximum: 4,
        })
    ));
    let second_lease = lease();
    let second = evaluator
        .evaluate(
            &population,
            2,
            limits(),
            provenance(1, 15, population.len(), 2),
            &CancelGate::new_clock_free(),
            budget(population.len() as u64),
            &second_lease,
            |x| vec![x[0], -x[0]],
        )
        .unwrap();
    assert!(matches!(
        publisher.publish(second),
        Err(PopulationPublishError::CumulativeWorkLimit {
            requested: 10,
            maximum: 5,
        })
    ));
}
