//! Diagnostic: find the support skin the owner's own link envelopes require.
//!
//! Lowering the declared surface by d is equivalent to widening the skin by d,
//! so sweeping the surface height measures the skin against the exact geometry
//! the collision path uses rather than against a sphere approximation.
#![cfg(test)]

use fs_cmaes_viz_wasm::manipulation::{
    manipulation_curriculum_mean, ManipulationConfig, ManipulationEvaluator, ManipulationTask,
    ObstacleBox, ARM_SUPPORT_SURFACE_SKIN_M,
};
use fs_ga::Vec3;
use fs_scene::BodyRole;

#[test]
fn report_required_support_skin() {
    for task in [
        ManipulationTask::KitchenMug,
        ManipulationTask::LivingRoomRemote,
        ManipulationTask::BackyardTrowel,
    ] {
        let probe = ManipulationEvaluator::new(ManipulationConfig {
            task,
            ..ManipulationConfig::default()
        })
        .expect("evaluator");
        let support = probe.scene().support_height_m;
        let mean = manipulation_curriculum_mean(task);
        let mut required = f64::NAN;
        for step in 0..40 {
            let drop = step as f64 * 0.0025;
            let evaluator = ManipulationEvaluator::new(ManipulationConfig {
                task,
                obstacles: vec![ObstacleBox {
                    center_m: Vec3::new(-0.85, 0.0, support - 0.045 - drop),
                    half_extents_m: Vec3::new(0.7, 0.825, 0.045),
                    yaw_rad: 0.0,
                    role: BodyRole::Support,
                }],
                ..ManipulationConfig::default()
            })
            .expect("evaluator");
            let receipt = evaluator.evaluate(&mean).expect("rollout");
            if receipt.placed {
                required = ARM_SUPPORT_SURFACE_SKIN_M + drop;
                break;
            }
        }
        eprintln!("PROBE task={task:?} required_skin={required:.3} m");
    }
}
