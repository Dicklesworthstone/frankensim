//! Caller-owned noisy simulations delivered out of order, with an in-memory
//! partial-batch checkpoint. Synthetic example, not experimental validation.
//!
//! cargo run -p fs-bo --example external_noisy_design

use fs_bo::{Kernel, Matern};
use fs_bo::noisy::{NoisyBoConfig, NoisyObservation};
use fs_bo::study::{NoisyRequest, NoisyStudy, TellDisposition};

fn simulate(request: &NoisyRequest) -> NoisyObservation {
    // In a real dispatcher this key travels with the remote physics job.
    // Never draw shared sequential randomness in arrival order.
    let mut rng = fs_rand::StreamKey { seed: 2026, kernel: 0x4558_544E,
        tile: u32::try_from(request.evaluation).expect("this example has 18 jobs") }.stream();
    let amplitude = 0.01 + 0.03 * request.x[0];
    NoisyObservation {
        value: (request.x[0] - 0.37).powi(2) + amplitude * (2.0 * rng.next_f64() - 1.0),
        noise_variance: amplitude * amplitude / 3.0,
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let config = NoisyBoConfig { bounds: (0.0, 1.0),
        kernel: Kernel { family: Matern::FiveHalves, signal: 0.04, lengthscales: vec![0.25] },
        prior_mean: 0.1, q: 2, mc_samples: 128, acq_starts: 2, acq_evals: 64, seed: 2026 };
    let mut study = NoisyStudy::new(2026, 1, 6, 6, &config);
    let first = study.ask()[0].clone();
    let observation = simulate(&first);
    study.tell(&first, observation)?;
    let checkpoint = study.clone();
    study = checkpoint;
    assert_eq!(study.tell(&first, observation)?, TellDisposition::Replayed);
    let mut physical_calls = 1;
    while !study.is_complete() {
        // A real job system may return this order, any other order, or partial
        // arrivals between polls. Do not dispatch an already in-flight ID twice.
        let mut jobs = study.ask();
        jobs.reverse();
        for job in jobs {
            let result = simulate(&job);
            physical_calls += 1;
            study.tell(&job, result)?;
        }
        study.advance()?;
    }
    assert_eq!(physical_calls, 18);
    let report = study.report();
    let best = report.incumbent_trace.last().expect("completed initial design");
    println!("physical_calls={physical_calls} x={:.8} posterior_mean={:.8} posterior_variance={:.8}",
        report.x[best.observation_index][0], best.mean, best.variance);
    println!("model_fit_attempts={} acquisition_batch_attempts={}",
        study.work().model_fit_attempts, study.work().acquisition_batch_attempts);
    println!("Synthetic model estimates; no optimizer-quality or physical-validation claim.");
    Ok(())
}
