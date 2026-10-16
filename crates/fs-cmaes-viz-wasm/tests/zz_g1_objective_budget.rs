//! Diagnostic: where does the walking objective's weight actually go?
#![cfg(test)]
use fs_cmaes_viz_wasm::g1_walking::{
    g1_walking_curriculum_mean, G1Challenge, G1Task, G1WalkingConfig, G1WalkingEvaluator,
};

#[test]
fn report_objective_term_budget() {
    let config = G1WalkingConfig {
        task: G1Task::Walking,
        challenge: G1Challenge::TerrainAndPush,
        ..G1WalkingConfig::default()
    };
    let evaluator = G1WalkingEvaluator::new(config.clone()).unwrap();
    let r = evaluator.evaluate(&g1_walking_curriculum_mean()).unwrap();
    let dur = r.completed_steps as f64 * config.step_s;
    let scale: f64 = config.target_forward_speed_m_per_s.max(0.25);
    let target_distance = (config.target_forward_speed_m_per_s * dur).max(0.10);
    let speed_tracking = r.speed_error_integral / (scale * scale * dur);
    let terms: [(&str, f64); 9] = [
        ("speed tracking      x20", 20.0 * speed_tracking),
        ("contact mismatch   x180", 180.0 * r.contact_schedule_mismatch_integral / dur),
        ("swing clearance    x100", 100.0 * r.swing_clearance_error_integral / dur),
        ("posture             x12", 12.0 * r.posture_integral / dur),
        ("lateral             x12", 12.0 * r.lateral_error_integral / dur),
        ("heading             x10", 10.0 * r.heading_error_integral / dur),
        ("joint limit          x8", 8.0 * r.joint_limit_integral / dur),
        ("impact               x6", 6.0 * r.impact_integral / dur),
        ("backward            x20", 20.0 * r.backward_distance_m / target_distance),
    ];
    let total: f64 = terms.iter().map(|(_, v)| v).sum();
    eprintln!("PROBE distance={:.3} m of {:.3} target; completed {} steps", r.distance_m, target_distance, r.completed_steps);
    for (label, value) in terms {
        eprintln!("PROBE {label}: {value:>9.2}  ({:>5.1}% of accounted)", 100.0 * value / total);
    }
    eprintln!("PROBE accounted total={total:.2}");
}
