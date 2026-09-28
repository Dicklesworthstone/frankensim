//! Declared heteroscedastic objective noise, fixed-prior q-NEI, and posterior
//! design recommendations. Synthetic numerical example, not a physical model.
//!
//! cargo run -p fs-bo --example noisy_design

use fs_bo::noisy::{NoisyBoConfig, NoisyObservation, minimize_noisy};
use fs_bo::{Kernel, Matern};

fn main() {
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
    let result = minimize_noisy(&mut objective, 1, 6, 6, &config);
    for (batch, incumbent) in result.incumbent_trace.iter().enumerate() {
        println!(
            "batch={batch} x={:.8} posterior_mean={:.8} posterior_variance={:.8}",
            result.x[incumbent.observation_index][0], incumbent.mean, incumbent.variance,
        );
    }
    println!("evaluations={} (model estimates, not certificates)", result.x.len());
}
