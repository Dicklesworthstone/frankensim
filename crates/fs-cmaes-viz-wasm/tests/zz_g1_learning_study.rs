//! Diagnostic: does the walking search actually make the robot walk?
//!
//! Runs the real owner under real CMA and reports the PHYSICAL result of the
//! best policy as generations pass, so a falling objective can be checked
//! against distance and survival instead of taken on trust.
#![cfg(test)]

use fs_cmaes_viz_wasm::g1_walking::{
    g1_walking_curriculum_mean, G1Challenge, G1Task, G1WalkingConfig, G1WalkingEvaluator,
};
use fs_dfo::cma_family::{CmaConfig, CmaFamily, CmaOptimizer};
use std::time::Instant;

struct Study {
    label: &'static str,
    duration_s: f64,
    target_speed: f64,
    sigma: f64,
    generations: usize,
}

fn evaluate_batch(config: &G1WalkingConfig, candidates: &[Vec<f64>]) -> Vec<f64> {
    let threads = std::thread::available_parallelism().map(|v| v.get()).unwrap_or(4);
    let chunk = candidates.len().div_ceil(threads);
    let mut objectives = vec![0.0; candidates.len()];
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for (index, slice) in candidates.chunks(chunk).enumerate() {
            let config = config.clone();
            handles.push((index, scope.spawn(move || {
                let evaluator = G1WalkingEvaluator::new(config).expect("evaluator");
                slice
                    .iter()
                    .map(|c| evaluator.evaluate(c).map(|r| r.objective).unwrap_or(f64::INFINITY))
                    .collect::<Vec<_>>()
            })));
        }
        for (index, handle) in handles {
            let part = handle.join().expect("worker");
            for (offset, value) in part.into_iter().enumerate() {
                objectives[index * chunk + offset] = value;
            }
        }
    });
    objectives
}

#[test]
#[ignore = "diagnostic study; run explicitly"]
fn report_whether_the_search_makes_the_robot_walk() {
    let studies = [
        Study { label: "shipped (1.5s, 0.65, 5e-4)", duration_s: 1.5, target_speed: 0.65, sigma: 0.0005, generations: 200 },

    ];

    for study in studies {
        let config = G1WalkingConfig {
            task: G1Task::Walking,
            challenge: G1Challenge::TerrainAndPush,
            duration_s: study.duration_s,
            target_forward_speed_m_per_s: study.target_speed,
            ..G1WalkingConfig::default()
        };
        let probe = G1WalkingEvaluator::new(config.clone()).expect("evaluator");
        let mean = g1_walking_curriculum_mean();
        let seed_receipt = probe.evaluate(&mean).expect("seed");
        let population = 16;
        let mut optimizer = CmaOptimizer::new(CmaConfig {
            family: CmaFamily::LmMa,
            mean: mean.to_vec(),
            sigma: study.sigma,
            max_evaluations: population * (study.generations + 2),
            seed: 0x4731_5050,
            population_size: Some(population),
            memory: Some(12),
        })
        .expect("cma");

        let started = Instant::now();
        let mut best_objective = seed_receipt.objective;
        let mut best_policy = mean.to_vec();
        let seed_dur = seed_receipt.completed_steps as f64 * config.step_s;
        eprintln!(
            "\nSTUDY {} \n  SEED     obj={:>9.3} dist={:.3} m  speed={:.3} m/s  steps={:>4}  | contact={:.4} slip={:.4} clearance={:.4} single={:.2}s double={:.2}s flight={:.2}s",
            study.label,
            seed_receipt.objective,
            seed_receipt.distance_m,
            seed_receipt.distance_m / seed_dur,
            seed_receipt.completed_steps,
            seed_receipt.contact_schedule_mismatch_integral / seed_dur,
            seed_receipt.slip_integral / seed_dur,
            seed_receipt.swing_clearance_error_integral / seed_dur,
            seed_receipt.single_support_s,
            seed_receipt.double_support_s,
            seed_receipt.flight_s,
        );
        for generation in 1..=study.generations {
            let ask = optimizer.ask().expect("ask");
            let candidates: Vec<Vec<f64>> = ask.candidates().iter().map(|c| c.to_vec()).collect();
            let objectives = evaluate_batch(&config, &candidates);
            for (candidate, objective) in candidates.iter().zip(objectives.iter()) {
                if *objective < best_objective {
                    best_objective = *objective;
                    best_policy = candidate.clone();
                }
            }
            optimizer.tell(&ask, &objectives).expect("tell");
            if generation % 100 == 0 || generation == study.generations {
                let receipt = probe.evaluate(&best_policy).expect("best");
                let dur = receipt.completed_steps as f64 * config.step_s;
                eprintln!(
                    "  gen {:>4} obj={:>9.3} dist={:.3} m  speed={:.3} m/s  steps={:>4}  | contact={:.4} slip={:.4} clearance={:.4} single={:.2}s double={:.2}s flight={:.2}s  ({:.0}s)",
                    generation,
                    receipt.objective,
                    receipt.distance_m,
                    receipt.distance_m / dur,
                    receipt.completed_steps,
                    receipt.contact_schedule_mismatch_integral / dur,
                    receipt.slip_integral / dur,
                    receipt.swing_clearance_error_integral / dur,
                    receipt.single_support_s,
                    receipt.double_support_s,
                    receipt.flight_s,
                    started.elapsed().as_secs_f64()
                );
            }
        }
    }
}
