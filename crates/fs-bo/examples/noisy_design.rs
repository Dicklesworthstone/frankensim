//! Declared heteroscedastic objective noise, q-NEI, and posterior design
//! recommendations. Synthetic numerical example, not a physical model.
//!
//! cargo run -p fs-bo --example noisy_design
//! cargo run -p fs-bo --example noisy_design -- --learn-kernel

use fs_bo::hyper::HeteroFitConfig;
use fs_bo::learning::{NoisyLearningConfig, minimize_noisy_with_learning};
use fs_bo::noisy::{NoisyBoConfig, NoisyObservation, minimize_noisy};
use fs_bo::{Kernel, Matern};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut arguments = std::env::args().skip(1);
    let learn = match arguments.next().as_deref() {
        None => false,
        Some("--learn-kernel") => true,
        _ => return Err("usage: noisy_design [--learn-kernel]".into()),
    };
    if arguments.next().is_some() {
        return Err("usage: noisy_design [--learn-kernel]".into());
    }
    let config = NoisyBoConfig {
        bounds: (0.0, 1.0),
        kernel: Kernel {
            family: Matern::FiveHalves,
            signal: 0.04,
            lengthscales: vec![0.25],
        },
        prior_mean: 0.1,
        q: 2,
        mc_samples: 128,
        acq_starts: 2,
        acq_evals: 64,
        seed: 2026,
    };
    let mut noise = fs_rand::StreamKey {
        seed: 2026,
        kernel: 0x4F42_5345,
        tile: 0,
    }
    .stream();
    let mut objective = |x: &[f64]| {
        let amplitude = 0.01 + 0.03 * x[0];
        NoisyObservation {
            value: (x[0] - 0.37).powi(2) + amplitude * (2.0 * noise.next_f64() - 1.0),
            // Variance of uniform [-amplitude, amplitude] measurement noise.
            noise_variance: amplitude * amplitude / 3.0,
        }
    };
    let result = if learn {
        let learning = NoisyLearningConfig {
            refit_every: 2,
            fit: HeteroFitConfig {
                lengthscale_bounds: vec![(0.03, 2.0)],
                signal_bounds: (1e-4, 1.0),
                starts: 3, max_iterations: 40, max_evaluations: 100,
                gradient_tolerance: 1e-7, seed: 2026,
            },
        };
        let learned = minimize_noisy_with_learning(&mut objective, 1, 6, 6, &config, &learning)?;
        for fit in &learned.fits {
            println!("refit_after_batches={} observations={} lengthscale={:.8} signal_variance={:.8} lml_gain={:.8} likelihood_probes={}",
                fit.after_batches, fit.observation_count, fit.kernel.lengthscales[0],
                fit.kernel.signal, fit.lml - fit.initial_lml, fit.evaluations);
        }
        println!("likelihood_probes={} posterior_only_fits={}",
            learned.likelihood_evaluations, learned.posterior_only_fits);
        learned.report
    } else {
        minimize_noisy(&mut objective, 1, 6, 6, &config)
    };
    for (batch, incumbent) in result.incumbent_trace.iter().enumerate() {
        println!(
            "batch={batch} x={:.8} posterior_mean={:.8} posterior_variance={:.8}",
            result.x[incumbent.observation_index][0], incumbent.mean, incumbent.variance,
        );
    }
    println!("evaluations={} (model estimates, not certificates)", result.x.len());
    Ok(())
}
