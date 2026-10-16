//! What work surface can this arm actually work at?
//!
//! Declaring the counter as drawn refuses all three tasks. Either the policies
//! are wrong or the furniture is: a 1.65 m deep top whose near edge is 150 mm
//! from a floor-mounted arm's base axis is not a workcell anyone would build.
//! This searches table geometries for one that holds both object stations and
//! leaves the arm room to work.
#![cfg(test)]

use fs_cmaes_viz_wasm::manipulation::{
    manipulation_curriculum_mean, ManipulationConfig, ManipulationEvaluator, ManipulationTask,
    ObstacleBox,
};
use fs_ga::Vec3;
use fs_scene::BodyRole;

const TASKS: [ManipulationTask; 3] = [
    ManipulationTask::KitchenMug,
    ManipulationTask::LivingRoomRemote,
    ManipulationTask::BackyardTrowel,
];

#[test]
fn report_a_table_the_arm_can_work_at() {
    // Where the objects actually are: the top has to hold these.
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    let mut y_abs = 0.0_f64;
    for task in TASKS {
        let e = ManipulationEvaluator::new(ManipulationConfig { task, ..Default::default() })
            .expect("evaluator");
        let s = e.scene();
        for p in [s.initial_object_position_m, s.goal_object_position_m] {
            x_lo = x_lo.min(p.x);
            x_hi = x_hi.max(p.x);
            y_abs = y_abs.max(p.y.abs());
        }
        eprintln!(
            "PROBE {task:?} start=({:.2},{:.2}) goal=({:.2},{:.2}) support={:.4}",
            s.initial_object_position_m.x, s.initial_object_position_m.y,
            s.goal_object_position_m.x, s.goal_object_position_m.y, s.support_height_m
        );
    }
    eprintln!("PROBE stations span x=[{x_lo:.2},{x_hi:.2}] |y|<={y_abs:.2}");

    // What each task does with NO table declared, so a failure below can be
    // attributed to the table rather than inherited from the policy.
    for task in TASKS {
        let e = ManipulationEvaluator::new(ManipulationConfig { task, ..Default::default() })
            .expect("evaluator");
        let r = e.evaluate(&manipulation_curriculum_mean(task)).expect("rollout");
        eprintln!(
            "PROBE UNDECLARED {task:?}: placed={} err={:.4} lift={:.3} grasped={}",
            r.placed, r.final_object_error_m, r.maximum_lift_m, r.ever_grasped
        );
    }

    // Candidate tops, each with a 60 mm margin around the stations at minimum.
    let candidates: [(f64, f64, f64); 6] = [
        // (near edge x, far edge x, y half-extent)
        (-0.15, -1.55, 0.825), // as drawn today
        (-0.35, -1.45, 0.55),
        (-0.35, -1.25, 0.45),
        (-0.40, -1.20, 0.45),
        (-0.40, -1.10, 0.40),
        (-0.45, -1.05, 0.40),
    ];

    for (near, far, y_half) in candidates {
        let center_x = 0.5 * (near + far);
        let half_x = 0.5 * (near - far).abs();
        let holds_stations = center_x + half_x >= x_hi - 1e-9
            && center_x - half_x <= x_lo + 1e-9
            && y_half >= y_abs;
        let mut verdict = Vec::new();
        for task in TASKS {
            let probe = ManipulationEvaluator::new(ManipulationConfig { task, ..Default::default() })
                .expect("evaluator");
            let support = probe.scene().support_height_m;
            let evaluator = ManipulationEvaluator::new(ManipulationConfig {
                task,
                obstacles: vec![ObstacleBox {
                    center_m: Vec3::new(center_x, 0.0, support - 0.045),
                    half_extents_m: Vec3::new(half_x, y_half, 0.045),
                    yaw_rad: 0.0,
                    role: BodyRole::Support,
                }],
                ..Default::default()
            })
            .expect("evaluator");
            let r = evaluator
                .evaluate(&manipulation_curriculum_mean(task))
                .expect("rollout");
            verdict.push(r.placed);
        }
        eprintln!(
            "PROBE top x=[{far:.2},{near:.2}] y=+-{y_half:.2} holds_stations={holds_stations} placed={verdict:?}"
        );
    }
}
