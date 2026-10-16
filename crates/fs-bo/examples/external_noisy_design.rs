//! Caller-owned noisy simulations, learned kernels and process checkpoints.
//! Synthetic example, not experimental validation. Resume only checkpoints
//! produced by this example under the same numerical build/profile.
//!
//! cargo run -p fs-bo --example external_noisy_design
//! cargo run -p fs-bo --example external_noisy_design -- --learn-kernel --checkpoint study.bin
//! cargo run -p fs-bo --example external_noisy_design -- --resume study.bin

use std::io::{Read, Write};
use fs_bo::{Kernel, Matern};
use fs_bo::hyper::HeteroFitConfig;
use fs_bo::learning::NoisyLearningConfig;
use fs_bo::noisy::{NoisyBoConfig, NoisyObservation};
use fs_bo::study::{CheckpointLimits, NoisyRequest, NoisyStudy};

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
    let mut args = std::env::args().skip(1);
    let mut learn = false;
    let mut save: Option<String> = None;
    let mut resume: Option<String> = None;
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--learn-kernel" if !learn => learn = true,
            "--checkpoint" if save.is_none() && resume.is_none() => {
                save = Some(args.next().ok_or("--checkpoint needs a fresh path")?);
            }
            "--resume" if save.is_none() && resume.is_none() => {
                resume = Some(args.next().ok_or("--resume needs a checkpoint path")?);
            }
            _ => return Err("usage: external_noisy_design [--learn-kernel] [--checkpoint PATH] | --resume PATH".into()),
        }
    }
    if learn && resume.is_some() { return Err("a checkpoint already owns its learning policy".into()); }
    // Refuse an existing output before spending any simulation work.
    let mut checkpoint_file = if let Some(path) = &save {
        Some(std::fs::OpenOptions::new().write(true).create_new(true).open(path)?)
    } else { None };
    let limits = CheckpointLimits { max_bytes: 1 << 20, max_evaluations: 18 };
    let mut study = if let Some(path) = &resume {
        // Bound the file read as well as the decoder; never allocate an entire
        // oversized file before consulting the checkpoint limits.
        let file = std::fs::File::open(path)?;
        let mut bytes = Vec::new();
        file.take(limits.max_bytes as u64 + 1).read_to_end(&mut bytes)?;
        NoisyStudy::from_checkpoint(&bytes, limits)?
    } else {
        let config = NoisyBoConfig { bounds: (0.0, 1.0),
            kernel: Kernel { family: Matern::FiveHalves, signal: 0.04, lengthscales: vec![0.25] },
            prior_mean: 0.1, q: 2, mc_samples: 128, acq_starts: 2, acq_evals: 64, seed: 2026 };
        if learn {
            let policy = NoisyLearningConfig { refit_every: 2, fit: HeteroFitConfig {
                lengthscale_bounds: vec![(0.03, 2.0)], signal_bounds: (1e-4, 1.0),
                starts: 3, max_iterations: 40, max_evaluations: 100,
                gradient_tolerance: 1e-7, seed: 2026,
            } };
            NoisyStudy::new_with_learning(2026, 1, 6, 6, &config, &policy)?
        } else { NoisyStudy::new(2026, 1, 6, 6, &config) }
    };
    let retained = study.report().x.len() + study.pending().iter().filter(|p| p.observation.is_some()).count();
    let mut physical_calls = 0;
    while !study.is_complete() {
        let mut jobs = study.ask(); jobs.reverse();
        for job in jobs {
            study.tell(&job, simulate(&job))?;
            physical_calls += 1;
            if let (Some(path), Some(file)) = (&save, &mut checkpoint_file) {
                // Stop with a modeled prefix AND a partial adaptive batch, so
                // the file exercises accepted kernels and in-flight results.
                if !study.report().x.is_empty() {
                    file.write_all(&study.checkpoint_bytes())?; file.sync_all()?;
                    println!("checkpoint={path} retained_results={} committed_refits={}",
                        retained + physical_calls, study.fit_records().len());
                    return Ok(());
                }
            }
        }
        study.advance()?;
        // Same roundtrip as file restoration; no model work or callbacks.
        study = NoisyStudy::from_checkpoint(&study.checkpoint_bytes(), limits)?;
    }
    assert_eq!(retained + physical_calls, 18);
    let report = study.report();
    let best = report.incumbent_trace.last().expect("completed initial design");
    println!("retained_results={retained} physical_calls_this_process={physical_calls} x={:.8} posterior_mean={:.8} posterior_variance={:.8}",
        report.x[best.observation_index][0], best.mean, best.variance);
    println!("model_fit_attempts={} acquisition_batch_attempts={} committed_refits={}",
        study.work().model_fit_attempts, study.work().acquisition_batch_attempts, study.fit_records().len());
    println!("Synthetic model estimates; no optimizer-quality or physical-validation claim.");
    Ok(())
}
