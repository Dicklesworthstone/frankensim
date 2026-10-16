//! Synthetic fan-power design with noisy temperature and sound constraints.
//! These declared algebraic responses are NOT a validated cooling/acoustic model.
//! cargo run -p fs-bo --example constrained_cooling
use fs_bo::{Kernel, Matern};
use fs_bo::constrained::{ConstrainedBoConfig, ConstrainedObservation, OutcomeConstraint, minimize_constrained};
use fs_bo::hyper::HeteroFitConfig;
use fs_bo::learning::NoisyLearningConfig;
use fs_bo::noisy::{NoisyBoConfig, NoisyObservation};

fn kernel(signal: f64) -> Kernel {
    Kernel { family: Matern::FiveHalves, signal, lengthscales: vec![0.3] }
}
fn learning(signal_bounds: (f64, f64)) -> Option<NoisyLearningConfig> {
    Some(NoisyLearningConfig { refit_every: 2, fit: HeteroFitConfig {
        lengthscale_bounds: vec![(0.08, 2.0)], signal_bounds, starts: 2,
        max_iterations: 12, max_evaluations: 24, gradient_tolerance: 1e-6, seed: 2026,
    }})
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let c = ConstrainedBoConfig {
        search: NoisyBoConfig { bounds: (0.0, 1.0), kernel: kernel(100.0), prior_mean: 10.0,
            q: 2, mc_samples: 64, acq_starts: 1, acq_evals: 32, seed: 2026 },
        objective_learning: learning((1.0, 400.0)),
        constraints: vec![
            OutcomeConstraint { kernel: kernel(400.0), prior_mean: 320.0,
                upper_bound: 310.0, learning: learning((10.0, 1600.0)) },
            OutcomeConstraint { kernel: kernel(100.0), prior_mean: 55.0,
                upper_bound: 70.0, learning: None },
        ],
        reference: 25.0,
        recommendation_probability: 0.9,
    };
    let mut random = fs_rand::StreamKey { seed: 2026, kernel: 0x434F_4F4C, tile: 0 }.stream();
    let mut noisy = |mean: f64, amplitude: f64| NoisyObservation {
        value: mean + amplitude * (2.0 * random.next_f64() - 1.0),
        noise_variance: amplitude * amplitude / 3.0,
    };
    let r = minimize_constrained(&mut |x| {
        let speed = x[0];
        ConstrainedObservation {
            objective: noisy(20.0 * speed.powi(3), 0.3),
            constraints: vec![noisy(350.0 - 80.0 * speed, 1.5), noisy(45.0 + 35.0 * speed, 1.0)],
        }
    }, 1, 6, 4, &c)?;
    for (stage, recommendation) in r.recommendations.iter().enumerate() {
        match recommendation {
            Some(best) => println!(
                "batch={stage} speed={:.6} power_mean_w={:.6} temperature_mean_k={:.6} sound_mean_db={:.6} modeled_joint_feasibility={:.6}",
                r.x[best.observation_index][0], best.objective_mean, best.constraint_means[0],
                best.constraint_means[1], best.joint_feasibility),
            None => println!("batch={stage} no evaluated design meets the modeled feasibility threshold"),
        }
    }
    println!("evaluations={} likelihood_probes={} posterior_only_fits={}",
        r.x.len(), r.likelihood_evaluations, r.posterior_only_fits);
    println!("Synthetic model estimates; candidate evaluations are not guaranteed feasible.");
    Ok(())
}
