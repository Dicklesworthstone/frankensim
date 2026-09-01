use fs_dfo::{
    PopulationCandidate, PopulationEvaluator, PopulationLimits, PopulationProvenance,
    PopulationPublishError, PopulationPublisher, PopulationRefusal,
};
use fs_exec::{CancelGate, PoolConfig, RunId, TilePool};

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
        seed: 0x7A24_1600,
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

#[test]
fn tiles_preserve_semantic_order_and_replay_across_worker_counts() {
    let population = candidates();
    let run = provenance(0, 9, population.len(), 2);
    let objective = |x: &[f64]| vec![x[0], x[0] * x[0]];
    let one = PopulationEvaluator::new(&TilePool::new(PoolConfig::for_host(1, run.seed)), 2)
        .evaluate(&population, 2, limits(), run, &CancelGate::new(), objective)
        .expect("one worker evaluates");
    let many = PopulationEvaluator::new(&TilePool::new(PoolConfig::for_host(4, run.seed)), 2)
        .evaluate(&population, 2, limits(), run, &CancelGate::new(), objective)
        .expect("many workers evaluates");
    assert_eq!(
        one, many,
        "G5: tile placement cannot alter accepted generation bits"
    );
    assert_eq!(
        one.evaluations
            .iter()
            .map(|row| row.identity)
            .collect::<Vec<_>>(),
        vec![40, 10, 30, 20, 50]
    );
}

#[test]
fn cancellation_drains_and_never_publishes_a_partial_generation() {
    let population = candidates();
    let gate = CancelGate::new_clock_free();
    gate.request();
    let publisher = PopulationPublisher::new();
    let result = PopulationEvaluator::new(&TilePool::new(PoolConfig::for_host(2, 3)), 1).evaluate(
        &population,
        2,
        limits(),
        provenance(0, 3, population.len(), 1),
        &gate,
        |x| vec![x[0], x[0]],
    );
    assert!(
        matches!(result, Err(PopulationPublishError::Executor(_))),
        "G4 cancellation is an executor outcome"
    );
    assert_eq!(
        publisher.checkpoint().committed,
        None,
        "no partial population can mint a generation"
    );
    assert_eq!(publisher.checkpoint().committed_work_units, 0);
}

#[test]
fn preflight_refuses_duplicate_identity_before_callback_or_publication() {
    let mut population = candidates();
    population[4].identity = 40;
    let calls = std::sync::atomic::AtomicUsize::new(0);
    let result = PopulationEvaluator::new(&TilePool::new(PoolConfig::for_host(2, 4)), 2).evaluate(
        &population,
        2,
        limits(),
        provenance(0, 4, population.len(), 2),
        &CancelGate::new_clock_free(),
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
    let evaluator = PopulationEvaluator::new(&TilePool::new(PoolConfig::for_host(2, 5)), 2);
    let publisher = PopulationPublisher::new();
    let first = evaluator
        .evaluate(
            &population,
            2,
            limits(),
            provenance(0, 5, population.len(), 2),
            &CancelGate::new(),
            |x| vec![x[0], -x[0]],
        )
        .unwrap();
    publisher.publish(first).unwrap();
    let resumed = PopulationPublisher::from_checkpoint(publisher.checkpoint());
    let second = evaluator
        .evaluate(
            &population,
            2,
            limits(),
            provenance(1, 6, population.len(), 2),
            &CancelGate::new(),
            |x| vec![x[0], -x[0]],
        )
        .unwrap();
    resumed.publish(second).unwrap();
    assert_eq!(resumed.checkpoint().committed_work_units, 10);
    assert_eq!(resumed.checkpoint().committed.unwrap().generation, 1);
}
