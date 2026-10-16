//! Diagnostic: is the arm blocked by the surface under the object, or by the
//! rest of the tabletop its links sweep across to get there?
#![cfg(test)]

use fs_cmaes_viz_wasm::manipulation::{
    manipulation_curriculum_mean, ManipulationConfig, ManipulationEvaluator, ManipulationTask,
    ObstacleBox,
};
use fs_ga::Vec3;
use fs_scene::BodyRole;

#[test]
fn report_which_part_of_the_surface_blocks_the_task() {
    let task = ManipulationTask::LivingRoomRemote;
    let probe = ManipulationEvaluator::new(ManipulationConfig { task, ..Default::default() })
        .expect("evaluator");
    let scene = probe.scene();
    let support = scene.support_height_m;
    let start = scene.initial_object_position_m;
    let goal = scene.goal_object_position_m;
    let mean = manipulation_curriculum_mean(task);

    let run = |label: &str, center: Vec3, half: Vec3| {
        let evaluator = ManipulationEvaluator::new(ManipulationConfig {
            task,
            obstacles: vec![ObstacleBox {
                center_m: center,
                half_extents_m: half,
                yaw_rad: 0.0,
                role: BodyRole::Support,
            }],
            ..Default::default()
        })
        .expect("evaluator");
        let r = evaluator.evaluate(&mean).expect("rollout");
        eprintln!("PROBE {label}: placed={} lift={:.3}", r.placed, r.maximum_lift_m);
    };

    eprintln!(
        "PROBE support={support:.4} start=({:.2},{:.2}) goal=({:.2},{:.2})",
        start.x, start.y, goal.x, goal.y
    );
    let z = Vec3::new(0.0, 0.0, support - 0.045);
    run("full tabletop", Vec3::new(-0.85, 0.0, z.z), Vec3::new(0.7, 0.825, 0.045));
    // Which x-band of the top conflicts? Each is a 0.1 m slice of the tabletop.
    for step in 0..14 {
        let x_lo = -1.55 + step as f64 * 0.1;
        let center_x = x_lo + 0.05;
        run(
            &format!("band x=[{:.2},{:.2}]", x_lo, x_lo + 0.1),
            Vec3::new(center_x, 0.0, z.z),
            Vec3::new(0.05, 0.825, 0.045),
        );
    }
}
