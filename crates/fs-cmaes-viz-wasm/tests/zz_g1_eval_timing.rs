#![cfg(test)]
use fs_cmaes_viz_wasm::g1_walking::{
    g1_walking_curriculum_mean, G1Challenge, G1Task, G1WalkingConfig, G1WalkingEvaluator,
};
use std::time::Instant;

#[test]
fn report_g1_evaluation_cost() {
    let evaluator = G1WalkingEvaluator::new(G1WalkingConfig {
        task: G1Task::Walking,
        challenge: G1Challenge::TerrainAndPush,
        ..G1WalkingConfig::default()
    })
    .unwrap();
    let mean = g1_walking_curriculum_mean();
    let started = Instant::now();
    for _ in 0..10 {
        let _ = evaluator.evaluate(&mean).unwrap();
    }
    let receipt = evaluator.evaluate(&mean).unwrap();
    eprintln!("PROBE per_eval_ms={:.1}", started.elapsed().as_secs_f64() * 100.0);
    eprintln!(
        "PROBE objective={:.4} distance={:.4} steps={} speed_err={:.4}",
        receipt.objective, receipt.distance_m, receipt.completed_steps, receipt.speed_error_integral
    );
}
