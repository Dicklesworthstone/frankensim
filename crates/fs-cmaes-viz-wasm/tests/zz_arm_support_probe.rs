//! Diagnostic: measure how far the arm's own links reach into the work
//! surface the browser draws, so the support skin is a measured number.
//!
//! The surface is a finite slab beside the arm, NOT a half-space: the arm is
//! floor-mounted at the owner origin and its base legitimately sits a third of
//! a metre below the table top, outside the slab's footprint.
#![cfg(test)]

use fs_cmaes_viz_wasm::manipulation::{
    manipulation_curriculum_mean, ManipulationConfig, ManipulationEvaluator, ManipulationTask,
};

const RADII: [f64; 7] = [0.082, 0.072, 0.068, 0.064, 0.058, 0.052, 0.046];

// The browser's `armCounterSlabObstacle`, in the owner frame: a 1.4 x 1.65 m
// slab 0.09 m thick whose top face is the admitted support height.
const SLAB_CENTER_XY: [f64; 2] = [-0.85, 0.0];
const SLAB_HALF: [f64; 3] = [0.7, 0.825, 0.045];

#[test]
fn report_link_penetration_into_the_work_surface() {
    for task in [
        ManipulationTask::KitchenMug,
        ManipulationTask::LivingRoomRemote,
        ManipulationTask::BackyardTrowel,
    ] {
        let evaluator = ManipulationEvaluator::new(ManipulationConfig {
            task,
            ..ManipulationConfig::default()
        })
        .expect("evaluator");
        let receipt = evaluator
            .trace(&manipulation_curriculum_mean(task))
            .expect("trace");
        let support = evaluator.scene().support_height_m;
        // Sink-limit body: the slab's footprint, from `skin` below the top
        // face downward through the whole column beneath the table. This is
        // what a support surface actually forbids -- going below it -- and it
        // works for a slab thinner than the skin, which an all-face shrink
        // would collapse to nothing.
        const SKIN: f64 = 0.05;
        const COLUMN_DEPTH: f64 = 1.0;
        let limit_top = support - SKIN;
        let center = [
            SLAB_CENTER_XY[0],
            SLAB_CENTER_XY[1],
            limit_top - 0.5 * COLUMN_DEPTH,
        ];
        let half = [SLAB_HALF[0], SLAB_HALF[1], 0.5 * COLUMN_DEPTH];
        let mut deepest = 0.0_f64;
        let mut deepest_link = usize::MAX;
        for sample in &receipt.trace {
            for (index, pose) in sample.link_pose.iter().enumerate() {
                let radius = RADII[index.min(RADII.len() - 1)];
                let depth = fs_scene::sphere_box_overlap_depth(
                    &[pose[0], pose[1], pose[2]],
                    radius,
                    &center,
                    &half,
                    1.0,
                    0.0,
                );
                if depth > deepest {
                    deepest = depth;
                    deepest_link = index;
                }
            }
        }
        eprintln!(
            "PROBE task={task:?} support={support:.4} deepest_below_sink_limit={deepest:.4} m (link {deepest_link})"
        );
    }
}
