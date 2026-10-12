//! Synthetic fan power, temperature and sound limits with unequal job durations.
//! No threads/sleeps and no measured hardware or timing-performance claim.
//! cargo run -p fs-bo --example asynchronous_constrained_cooling
use fs_bo::{Kernel, Matern};
use fs_bo::constrained::{AsyncConstrainedStudy, ConstrainedBoConfig,
    ConstrainedObservation, OutcomeConstraint};
use fs_bo::noisy::{NoisyBoConfig, NoisyObservation};
use fs_bo::study::NoisyRequest;
use std::collections::BTreeMap;

fn kernel(signal: f64) -> Kernel {
    Kernel { family: Matern::FiveHalves, signal, lengthscales: vec![0.3] }
}
fn simulate(request: &NoisyRequest) -> ConstrainedObservation {
    let mut rng = fs_rand::StreamKey { seed: 2026, kernel: 0x4143_4F4F,
        tile: u32::try_from(request.evaluation).expect("18 jobs") }.stream();
    let mut noisy = |mean: f64, amplitude: f64| NoisyObservation {
        value: mean + amplitude * (2.0 * rng.next_f64() - 1.0),
        noise_variance: amplitude * amplitude / 3.0,
    };
    let speed = request.x[0];
    ConstrainedObservation { objective: noisy(20.0 * speed.powi(3), 0.3),
        constraints: vec![noisy(350.0 - 80.0 * speed, 1.5), noisy(45.0 + 35.0 * speed, 1.0)] }
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let config = ConstrainedBoConfig {
        search: NoisyBoConfig { bounds: (0.0, 1.0), kernel: kernel(100.0), prior_mean: 10.0,
            q: 2, mc_samples: 128, acq_starts: 2, acq_evals: 64, seed: 2026 },
        objective_learning: None,
        constraints: vec![
            OutcomeConstraint { kernel: kernel(400.0), prior_mean: 320.0, upper_bound: 310.0, learning: None },
            OutcomeConstraint { kernel: kernel(100.0), prior_mean: 55.0, upper_bound: 70.0, learning: None },
        ], reference: 25.0, recommendation_probability: 0.9,
    };
    let mut study = AsyncConstrainedStudy::new(2026, 1, 6, 18, 3, &config)?;
    let mut running: BTreeMap<usize, (u64, NoisyRequest)> = BTreeMap::new();
    let mut now = 0u64;
    let mut completed = 0;
    while !study.is_complete() {
        // ask is not a lease: this map prevents duplicate physical dispatch.
        for request in study.ask() {
            running.entry(request.evaluation).or_insert_with(|| {
                let duration = if request.evaluation == 0 { 100 } else { 1 + request.evaluation as u64 % 3 };
                (now + duration, request)
            });
        }
        let (&id, &(time, _)) = running.iter().min_by_key(|(id, (time, _))| (*time, **id))
            .ok_or("no outstanding job before finalization")?;
        now = time;
        let (_, request) = running.remove(&id).expect("selected job");
        study.tell(&request, simulate(&request))?;
        completed += 1;
        let issued = study.refill()?;
        println!("event={completed} logical_time={now} received={id} new_jobs={issued} unanswered={}", study.in_flight());
    }
    assert_eq!(completed, 18);
    let report = study.model_report().expect("final modeled outcomes");
    if let Some(best) = &report.recommendation {
        let id = report.observed_evaluations[best.observation_index];
        println!("job={id} speed={:.6} power_mean_w={:.6} temperature_mean_k={:.6} sound_mean_db={:.6} modeled_joint_feasibility={:.6}",
            study.jobs()[id].request.x[0], best.objective_mean, best.constraint_means[0],
            best.constraint_means[1], best.joint_feasibility);
    } else {
        println!("No evaluated design meets the modeled joint-feasibility threshold.");
    }
    println!("Synthetic model estimates; physical evaluations are not guaranteed feasible.");
    Ok(())
}
