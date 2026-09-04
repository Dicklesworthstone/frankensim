//! The walking objective must pay for walking.
//!
//! Before v069 it did not. Measured at the curriculum seed, contact-schedule
//! matching was 79% of the shaping score and speed tracking 7.8%, with no
//! forward-progress reward at all — only a penalty for going backwards. A
//! 300-generation study showed exactly what that buys: the objective improved
//! from 1.32 to -11.94 while the distance walked FELL from 0.329 m to 0.260 m.
//! The optimizer was working perfectly; it had been asked to keep time rather
//! than to travel.
//!
//! This test fails if that ever comes back. It runs a real search against the
//! real owner and asserts the physical outcome, because an objective that
//! improves while the robot walks less far is the failure worth catching.
#![cfg(test)]

use fs_cmaes_viz_wasm::g1_walking::{
    g1_walking_curriculum_mean, G1Challenge, G1Task, G1WalkingConfig, G1WalkingEvaluator,
};
use fs_dfo::cma_family::{CmaConfig, CmaFamily, CmaOptimizer};

const POPULATION: usize = 16;
const GENERATIONS: usize = 80;

#[test]
fn a_short_search_makes_the_robot_walk_measurably_further() {
    let config = G1WalkingConfig {
        task: G1Task::Walking,
        challenge: G1Challenge::TerrainAndPush,
        ..G1WalkingConfig::default()
    };
    let evaluator = G1WalkingEvaluator::new(config.clone()).expect("evaluator");
    let mean = g1_walking_curriculum_mean();
    let seed = evaluator.evaluate(&mean).expect("seed rollout");

    let mut optimizer = CmaOptimizer::new(CmaConfig {
        family: CmaFamily::LmMa,
        mean: mean.to_vec(),
        sigma: 0.0005,
        max_evaluations: POPULATION * (GENERATIONS + 2),
        seed: 0x4731_5050,
        population_size: Some(POPULATION),
        memory: Some(12),
    })
    .expect("cma admission");

    let mut best_objective = seed.objective;
    let mut best_policy = mean.to_vec();
    for _ in 0..GENERATIONS {
        let ask = optimizer.ask().expect("ask");
        // Evaluated across threads: this is a real search against the real
        // owner, and a serial pass makes the guard slow enough that people
        // start skipping it.
        let candidates = ask.candidates();
        let threads = std::thread::available_parallelism()
            .map(|value| value.get())
            .unwrap_or(4);
        let chunk = candidates.len().div_ceil(threads);
        let mut objectives = vec![f64::INFINITY; candidates.len()];
        std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for (index, slice) in candidates.chunks(chunk).enumerate() {
                let config = config.clone();
                handles.push((
                    index,
                    scope.spawn(move || {
                        let local = G1WalkingEvaluator::new(config).expect("evaluator");
                        slice
                            .iter()
                            .map(|candidate| {
                                local
                                    .evaluate(candidate)
                                    .map(|receipt| receipt.objective)
                                    .unwrap_or(f64::INFINITY)
                            })
                            .collect::<Vec<_>>()
                    }),
                ));
            }
            for (index, handle) in handles {
                for (offset, value) in handle.join().expect("worker").into_iter().enumerate() {
                    objectives[index * chunk + offset] = value;
                }
            }
        });
        for (candidate, objective) in candidates.iter().zip(objectives.iter()) {
            if *objective < best_objective {
                best_objective = *objective;
                best_policy = candidate.clone();
            }
        }
        optimizer.tell(&ask, &objectives).expect("tell");
    }

    let learned = evaluator.evaluate(&best_policy).expect("learned rollout");

    // The objective improved...
    assert!(
        learned.objective < seed.objective,
        "search did not improve the objective: {} -> {}",
        seed.objective,
        learned.objective
    );
    // ...and the robot actually walks further for it. This is the assertion
    // that the pre-v069 weighting failed.
    assert!(
        learned.distance_m > seed.distance_m * 1.25,
        "objective improved but the robot did not walk further: {:.3} m -> {:.3} m",
        seed.distance_m,
        learned.distance_m
    );
    // Survival stays lexicographically primary: walking further must not be
    // bought by falling over sooner.
    assert_eq!(
        learned.completed_steps, seed.completed_steps,
        "the further-walking policy did not survive the full horizon"
    );
    // And it must not be bought by skating. Absolute slip is the wrong test —
    // covering twice the ground can cost a little more sliding in total and
    // still be a better walk — so the invariant is slip PER METRE TRAVELLED,
    // which is what sliding the feet forward to fake progress would ruin.
    let seed_slip_per_meter = seed.slip_integral / seed.distance_m;
    let learned_slip_per_meter = learned.slip_integral / learned.distance_m;
    assert!(
        learned_slip_per_meter < seed_slip_per_meter,
        "distance was bought with foot slip: {seed_slip_per_meter:.2} -> {learned_slip_per_meter:.2} slip per metre"
    );
}
