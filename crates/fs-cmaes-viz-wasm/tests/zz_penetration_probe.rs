//! Diagnostic: report the exact penetration the guard publishes for the
//! mid-body wall case, so the no-tunnelling bound is pinned to a measured
//! number rather than an assumed one.
#![cfg(test)]

use fs_cmaes_viz_wasm::g1_walking::{
    G1Challenge, G1Task, G1TerminationReason, G1WalkingConfig, G1WalkingEvaluator, ObstacleBox,
};
use fs_scene::BodyRole;

#[test]
fn report_mid_body_wall_penetration() {
    let wall = ObstacleBox {
        center_m: [0.20, 0.0, 0.75],
        half_extents_m: [0.03, 1.5, 1.0],
        yaw_rad: 0.0,
        role: BodyRole::KeepOut,
    };
    let cfg = G1WalkingConfig {
        task: G1Task::Walking,
        challenge: G1Challenge::Flat,
        obstacles: vec![wall],
        ..G1WalkingConfig::default()
    };
    let evaluator = G1WalkingEvaluator::new(cfg).expect("evaluator");
    let curriculum = fs_cmaes_viz_wasm::g1_walking::g1_walking_curriculum_mean();
    let receipt = evaluator.evaluate(&curriculum).expect("rollout runs");
    eprintln!(
        "PROBE reason={:?} steps={} penetration={:.6}",
        receipt.termination_reason, receipt.completed_steps, receipt.maximum_body_penetration_m
    );
    assert!(matches!(
        receipt.termination_reason,
        G1TerminationReason::BodyObstacle
    ));
}
