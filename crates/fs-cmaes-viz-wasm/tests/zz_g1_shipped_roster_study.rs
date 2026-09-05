//! Does the walking search improve with the roster the SITE actually declares?
//!
//! The 0.6.20 regression came from tuning a weight in a bare configuration:
//! no obstacles, where the page declares 48 house bodies. The weight that won
//! the bare sweep walked LESS far in the browser. This study carries the real
//! roster so a tuning question can be asked honestly and cheaply.
#![cfg(test)]

use fs_cmaes_viz_wasm::g1_walking::{
    g1_walking_curriculum_mean, G1Challenge, G1Task, G1WalkingConfig, G1WalkingEvaluator,
    ObstacleBox,
};
use fs_dfo::cma_family::{CmaConfig, CmaFamily, CmaOptimizer};
use fs_scene::BodyRole;

/// The roster the browser sends, dumped from G1_KERNEL_OBSTACLES.
fn house_roster() -> Vec<ObstacleBox> {
    vec![
    ObstacleBox { center_m: [1.7500, 1.0500, -0.1500], half_extents_m: [4.2000, 5.2000, 0.1500], yaw_rad: 0.0000, role: BodyRole::Support }, // house floor
    ObstacleBox { center_m: [0.5500, -1.3500, 0.3750], half_extents_m: [0.5500, 0.3000, 0.3750], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // library table
    ObstacleBox { center_m: [1.2500, -0.2500, 0.9000], half_extents_m: [0.2000, 0.2000, 0.9000], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // coat-rack
    ObstacleBox { center_m: [1.3500, -0.7500, 0.6000], half_extents_m: [0.3000, 0.1500, 0.6000], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // bookshelf-small
    ObstacleBox { center_m: [1.9500, 0.8500, 0.4250], half_extents_m: [0.6000, 0.2000, 0.4250], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // console-hall
    ObstacleBox { center_m: [-0.5375, 1.5500, 1.2500], half_extents_m: [1.2125, 0.0500, 1.2500], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // wall 7, segment 1
    ObstacleBox { center_m: [-1.7500, -0.4500, 0.6000], half_extents_m: [0.2500, 0.2500, 0.6000], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // plant-floor
    ObstacleBox { center_m: [-1.4500, -1.2500, 0.8000], half_extents_m: [0.1750, 0.1750, 0.8000], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // table-lamp
    ObstacleBox { center_m: [-1.5500, -1.3500, 0.2750], half_extents_m: [0.2500, 0.2500, 0.2750], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // lamp-table
    ObstacleBox { center_m: [1.9500, 0.6500, 0.7500], half_extents_m: [0.3000, 0.0150, 0.7500], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // mirror-hall
    ObstacleBox { center_m: [1.9500, -1.3500, 1.2500], half_extents_m: [0.8000, 0.0500, 1.2500], yaw_rad: -1.5708, role: BodyRole::KeepOut }, // wall 5, segment 1
    ObstacleBox { center_m: [-1.7500, 1.1500, 0.0150], half_extents_m: [0.0900, 0.1250, 0.0150], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // book-bedroom
    ObstacleBox { center_m: [1.4500, 2.1500, 1.2500], half_extents_m: [0.6000, 0.0500, 1.2500], yaw_rad: 1.5708, role: BodyRole::KeepOut }, // wall 8, segment 1
    ObstacleBox { center_m: [-1.5500, -2.3500, 0.6500], half_extents_m: [0.2000, 0.7000, 0.6500], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // fireplace
    ObstacleBox { center_m: [-0.8500, 3.1500, 0.2750], half_extents_m: [0.7000, 1.0000, 0.2750], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // bed-frame
    ObstacleBox { center_m: [-2.2500, 0.2500, 1.2500], half_extents_m: [5.5000, 0.0600, 1.2500], yaw_rad: -1.5708, role: BodyRole::KeepOut }, // wall 1, segment 1
    ObstacleBox { center_m: [-0.8500, 3.1500, 0.3000], half_extents_m: [0.7000, 0.9500, 0.3000], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // bed-master
    ObstacleBox { center_m: [2.6625, 1.5500, 1.2500], half_extents_m: [1.0375, 0.0500, 1.2500], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // wall 7, segment 2
    ObstacleBox { center_m: [-0.8500, 3.1500, 0.1250], half_extents_m: [0.7000, 1.0000, 0.1250], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // mattress
    ObstacleBox { center_m: [2.1500, 1.2500, 0.7500], half_extents_m: [0.3000, 0.0150, 0.7500], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // mirror
    ObstacleBox { center_m: [2.1500, 1.7500, 0.4250], half_extents_m: [0.4000, 0.2750, 0.4250], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // vanity
    ObstacleBox { center_m: [2.4500, -0.2500, 0.3000], half_extents_m: [0.1500, 0.1500, 0.3000], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // umbrella-stand
    ObstacleBox { center_m: [-0.6500, -2.3500, 0.1250], half_extents_m: [0.2250, 0.1000, 0.1250], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // cat-fireplace
    ObstacleBox { center_m: [2.5625, -0.4500, 1.2500], half_extents_m: [0.2125, 0.0500, 1.2500], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // wall 6, segment 1
    ObstacleBox { center_m: [2.7500, 1.6500, 0.4750], half_extents_m: [0.6000, 0.2500, 0.4750], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // dresser-2
    ObstacleBox { center_m: [2.7500, -1.3500, 0.4750], half_extents_m: [0.2250, 0.2500, 0.4750], yaw_rad: -3.1416, role: BodyRole::KeepOut }, // dining-chair-5
    ObstacleBox { center_m: [2.7500, -1.7500, 0.4750], half_extents_m: [0.2250, 0.2500, 0.4750], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // dining-chair-3
    ObstacleBox { center_m: [-1.5500, -3.3500, 0.4500], half_extents_m: [0.4750, 0.4500, 0.4500], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // armchair-1
    ObstacleBox { center_m: [0.8500, 3.3500, 0.4750], half_extents_m: [0.5500, 0.2500, 0.4750], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // dresser
    ObstacleBox { center_m: [3.2500, 0.1000, 0.0100], half_extents_m: [0.1350, 0.1350, 0.0100], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // plate-1
    ObstacleBox { center_m: [0.5500, -3.6500, 0.4000], half_extents_m: [0.9500, 0.4250, 0.4000], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // sofa
    ObstacleBox { center_m: [1.4500, 3.4500, 1.0000], half_extents_m: [0.7000, 0.3000, 1.0000], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // wardrobe-1
    ObstacleBox { center_m: [3.7500, 0.5500, 0.4600], half_extents_m: [0.5000, 0.3000, 0.4600], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // kitchen-island
    ObstacleBox { center_m: [3.3500, -2.3500, 0.3800], half_extents_m: [0.7000, 0.4500, 0.3800], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // dining-table
    ObstacleBox { center_m: [3.2500, 0.3000, 0.0500], half_extents_m: [0.0450, 0.0450, 0.0500], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // mug-1
    ObstacleBox { center_m: [3.4500, 0.1000, 0.0100], half_extents_m: [0.1350, 0.1350, 0.0100], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // plate-2
    ObstacleBox { center_m: [0.2500, 3.6500, 0.3000], half_extents_m: [0.1750, 0.1750, 0.3000], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // lamp-bedside-2
    ObstacleBox { center_m: [3.4500, 0.3000, 0.0650], half_extents_m: [0.0350, 0.0350, 0.0650], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // glass-1
    ObstacleBox { center_m: [0.1500, 3.7500, 0.2750], half_extents_m: [0.2500, 0.2500, 0.2750], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // nightstand-2
    ObstacleBox { center_m: [2.8500, 2.5500, 0.0250], half_extents_m: [0.3000, 0.0750, 0.0250], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // towel-rack
    ObstacleBox { center_m: [2.7500, -2.9500, 0.4750], half_extents_m: [0.2250, 0.2500, 0.4750], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // dining-chair-1
    ObstacleBox { center_m: [4.3375, -0.4500, 1.2500], half_extents_m: [0.6125, 0.0500, 1.2500], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // wall 6, segment 2
    ObstacleBox { center_m: [0.8500, -4.1500, 0.4500], half_extents_m: [0.4750, 0.4500, 0.4500], yaw_rad: 0.3000, role: BodyRole::KeepOut }, // armchair-2
    ObstacleBox { center_m: [4.5000, -0.7000, 0.4500], half_extents_m: [0.7500, 0.2500, 0.4500], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // sideboard
    ObstacleBox { center_m: [1.4500, 3.8000, 1.2500], half_extents_m: [0.2500, 0.0500, 1.2500], yaw_rad: 1.5708, role: BodyRole::KeepOut }, // wall 8, segment 2
    ObstacleBox { center_m: [-1.7500, 3.6500, 0.3000], half_extents_m: [0.1750, 0.1750, 0.3000], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // lamp-bedside-1
    ObstacleBox { center_m: [1.9500, -3.7500, 1.2500], half_extents_m: [0.4000, 0.0500, 1.2500], yaw_rad: -1.5708, role: BodyRole::KeepOut }, // wall 5, segment 2
    ObstacleBox { center_m: [-1.8500, 3.7500, 0.2750], half_extents_m: [0.2500, 0.2500, 0.2750], yaw_rad: 0.0000, role: BodyRole::KeepOut }, // nightstand-1
    ]
}

fn evaluate_batch(config: &G1WalkingConfig, candidates: &[Vec<f64>]) -> Vec<f64> {
    let threads = std::thread::available_parallelism().map(|v| v.get()).unwrap_or(4);
    let chunk = candidates.len().div_ceil(threads);
    let mut objectives = vec![f64::INFINITY; candidates.len()];
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for (index, slice) in candidates.chunks(chunk).enumerate() {
            let config = config.clone();
            handles.push((index, scope.spawn(move || {
                let evaluator = G1WalkingEvaluator::new(config).expect("evaluator");
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
fn report_walking_progress_with_the_shipped_roster() {
    let config = G1WalkingConfig {
        task: G1Task::Walking,
        challenge: G1Challenge::TerrainAndPush,
        obstacles: house_roster(),
        ..G1WalkingConfig::default()
    };
    let probe = G1WalkingEvaluator::new(config.clone()).expect("evaluator");
    let mean = g1_walking_curriculum_mean();
    let seed = probe.evaluate(&mean).expect("seed");
    eprintln!(
        "STUDY roster={} seed obj={:.3} dist={:.3} steps={}",
        config.obstacles.len(), seed.objective, seed.distance_m, seed.completed_steps
    );
    let population = 16;
    let generations = 320;
    let mut optimizer = CmaOptimizer::new(CmaConfig {
        family: CmaFamily::LmMa,
        mean: mean.to_vec(),
        sigma: 0.0005,
        max_evaluations: population * (generations + 2),
        seed: 0x4731_5050,
        population_size: Some(population),
        memory: Some(12),
    })
    .expect("cma");
    let mut best = seed.objective;
    let mut best_policy = mean.to_vec();
    for generation in 1..=generations {
        let ask = optimizer.ask().expect("ask");
        let candidates: Vec<Vec<f64>> = ask.candidates().iter().map(|c| c.to_vec()).collect();
        let objectives = evaluate_batch(&config, &candidates);
        for (candidate, objective) in candidates.iter().zip(objectives.iter()) {
            if *objective < best { best = *objective; best_policy = candidate.clone(); }
        }
        optimizer.tell(&ask, &objectives).expect("tell");
        if generation % 64 == 0 {
            let r = probe.evaluate(&best_policy).expect("best");
            eprintln!(
                "  gen {:>4} obj={:>9.3} dist={:.3} m steps={:>4}",
                generation, r.objective, r.distance_m, r.completed_steps
            );
        }
    }
}
