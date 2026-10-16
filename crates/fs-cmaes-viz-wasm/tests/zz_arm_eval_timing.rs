#![cfg(test)]
use fs_cmaes_viz_wasm::manipulation::{
    manipulation_curriculum_mean, ManipulationConfig, ManipulationEvaluator, ManipulationTask,
};
use std::time::Instant;

#[test]
fn report_evaluation_cost() {
    let task = ManipulationTask::LivingRoomRemote;
    let evaluator =
        ManipulationEvaluator::new(ManipulationConfig { task, ..Default::default() }).unwrap();
    let mean = manipulation_curriculum_mean(task);
    let started = Instant::now();
    for _ in 0..10 {
        let _ = evaluator.evaluate(&mean).unwrap();
    }
    eprintln!("PROBE per_evaluation_ms={:.1}", started.elapsed().as_secs_f64() * 100.0);
    eprintln!("PROBE threads={}", std::thread::available_parallelism().map(|v| v.get()).unwrap_or(1));
}
