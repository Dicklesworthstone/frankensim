//! Can the arm learn to do its task with the work surface declared solid?
//!
//! The curriculum policies were learned before the counter was a body, so they
//! sweep mid-arm links through it. This asks whether a search that STARTS from
//! them can find policies that keep the task and respect the surface.
#![cfg(test)]

use fs_cmaes_viz_wasm::manipulation::{
    manipulation_curriculum_mean, ManipulationConfig, ManipulationEvaluator, ManipulationTask,
    ObstacleBox, ARM_POLICY_DIMENSION,
};
use fs_dfo::cma_family::{CmaConfig, CmaFamily, CmaOptimizer};
use fs_ga::Vec3;
use fs_scene::BodyRole;

const POPULATION: usize = 24;

/// The counter the browser draws, in the owner frame.
fn counter(support_height_m: f64) -> ObstacleBox {
    ObstacleBox {
        center_m: Vec3::new(-0.85, 0.0, support_height_m - 0.045),
        half_extents_m: Vec3::new(0.7, 0.825, 0.045),
        yaw_rad: 0.0,
        role: BodyRole::Support,
    }
}

fn config_for(task: ManipulationTask, support_height_m: Option<f64>) -> ManipulationConfig {
    ManipulationConfig {
        task,
        obstacles: support_height_m.map(|h| vec![counter(h)]).unwrap_or_default(),
        ..Default::default()
    }
}

fn evaluate_batch(config: &ManipulationConfig, candidates: &[Vec<f64>]) -> Vec<f64> {
    let threads = std::thread::available_parallelism().map(|v| v.get()).unwrap_or(4);
    let chunk = candidates.len().div_ceil(threads);
    let mut objectives = vec![f64::INFINITY; candidates.len()];
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for (index, slice) in candidates.chunks(chunk).enumerate() {
            let config = config.clone();
            handles.push((index, scope.spawn(move || {
                let evaluator = ManipulationEvaluator::new(config).expect("evaluator");
                slice.iter()
                    .map(|c| evaluator.evaluate(c).map(|r| r.objective).unwrap_or(f64::INFINITY))
                    .collect::<Vec<_>>()
            })));
        }
        for (index, handle) in handles {
            for (offset, value) in handle.join().expect("worker").into_iter().enumerate() {
                objectives[index * chunk + offset] = value;
            }
        }
    });
    objectives
}

#[test]
#[ignore = "study; run explicitly"]
fn can_the_arm_learn_to_respect_its_own_table() {
    for task in [
        ManipulationTask::KitchenMug,
        ManipulationTask::LivingRoomRemote,
        ManipulationTask::BackyardTrowel,
    ] {
        let bare = ManipulationEvaluator::new(config_for(task, None)).expect("evaluator");
        let support = bare.scene().support_height_m;
        let mean = manipulation_curriculum_mean(task);
        let before = bare.evaluate(&mean).expect("baseline");

        let config = config_for(task, Some(support));
        let declared = ManipulationEvaluator::new(config.clone()).expect("evaluator");
        let blocked = declared.evaluate(&mean).expect("declared baseline");
        eprintln!(
            "\nSTUDY {task:?} support={support:.4}\n  seed undeclared: placed={} obj={:.3}\n  seed declared:   placed={} obj={:.3}",
            before.placed, before.objective, blocked.placed, blocked.objective
        );

        for &sigma in &[0.02_f64, 0.08] {
            let mut optimizer = CmaOptimizer::new(CmaConfig {
                family: CmaFamily::Full,
                mean: mean.to_vec(),
                sigma,
                max_evaluations: POPULATION * 402,
                seed: 0x5150_4152,
                population_size: Some(POPULATION),
                memory: None,
            })
            .expect("cma");
            let mut best = blocked.objective;
            let mut best_policy: Vec<f64> = mean.to_vec();
            for generation in 1..=400 {
                let ask = optimizer.ask().expect("ask");
                let candidates: Vec<Vec<f64>> =
                    ask.candidates().iter().map(|c| c.to_vec()).collect();
                let objectives = evaluate_batch(&config, &candidates);
                for (candidate, objective) in candidates.iter().zip(objectives.iter()) {
                    if *objective < best {
                        best = *objective;
                        best_policy = candidate.clone();
                    }
                }
                optimizer.tell(&ask, &objectives).expect("tell");
                if generation % 200 == 0 {
                    let r = declared.evaluate(&best_policy).expect("best");
                    eprintln!(
                        "  sigma={sigma} gen {generation}: placed={} obj={:.3} err={:.4} lift={:.3} dim={}",
                        r.placed, r.objective, r.final_object_error_m, r.maximum_lift_m,
                        ARM_POLICY_DIMENSION
                    );
                }
            }
        }
    }
}
