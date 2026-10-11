//! Simulated external workers with unequal completion times. No real threads,
//! sleeps or measured performance claim; all objective noise is keyed by job.
//! cargo run -p fs-bo --example asynchronous_noisy_design
//! cargo run -p fs-bo --example asynchronous_noisy_design -- --cached-jitter 1e-8

use fs_bo::async_study::AsyncNoisyStudy;
use fs_bo::noisy::{NoisyBoConfig, NoisyObservation};
use fs_bo::study::NoisyRequest;
use fs_bo::{Kernel, Matern};
use std::collections::BTreeMap;

fn simulate(r: &NoisyRequest) -> NoisyObservation {
    let mut noise = fs_rand::StreamKey { seed: 2026, kernel: 0x4153_594E,
        tile: u32::try_from(r.evaluation).expect("18 jobs") }.stream();
    let amplitude = 0.01 + 0.03 * r.x[0];
    NoisyObservation { value: (r.x[0] - 0.37).powi(2) + amplitude * (2.0 * noise.next_f64() - 1.0),
        noise_variance: amplitude * amplitude / 3.0 }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let arguments: Vec<_> = std::env::args().skip(1).collect();
    let jitter = match arguments.as_slice() {
        [] => None,
        [flag, value] if flag == "--cached-jitter" => {
            let value: f64 = value.parse()?;
            if !value.is_finite() || value < 0.0 {
                return Err("cached posterior jitter must be finite and nonnegative".into());
            }
            Some(value)
        }
        _ => return Err("usage: asynchronous_noisy_design [--cached-jitter variance]".into()),
    };
    if let Some(jitter) = jitter { println!("cached_posterior_jitter={jitter:e}"); }
    let c = NoisyBoConfig { bounds: (0.0, 1.0),
        kernel: Kernel { family: Matern::FiveHalves, signal: 0.04, lengthscales: vec![0.25] },
        prior_mean: 0.1, q: 2, mc_samples: 128, acq_starts: 2, acq_evals: 64, seed: 2026 };
    let mut s = AsyncNoisyStudy::new(2026, 1, 6, 18, 4, &c)?;
    let mut running: BTreeMap<usize, (u64, NoisyRequest)> = BTreeMap::new();
    let mut now = 0u64;
    let mut completed = 0;
    while !s.is_complete() {
        // ask is a view, NOT a claim/lease operation. This map prevents duplicate
        // dispatch while a real remote worker would still be processing an ID.
        for r in s.ask() {
            running.entry(r.evaluation).or_insert_with(|| {
                let duration = if r.evaluation == 0 { 100 } else { 1 + (r.evaluation as u64 % 3) };
                (now + duration, r)
            });
        }
        let (&id, &(time, _)) = running.iter().min_by_key(|(id, (time, _))| (*time, **id))
            .ok_or("no outstanding jobs before completion")?;
        now = time;
        let (_, request) = running.remove(&id).expect("selected job");
        s.tell(&request, simulate(&request))?;
        completed += 1;
        let replacements = match jitter {
            Some(jitter) => s.refill_cached(jitter)?,
            None => s.refill()?,
        };
        println!("event={completed} logical_time={now} received={id} new_jobs={replacements} unanswered={}", s.in_flight());
    }
    assert_eq!(completed, 18);
    let best = s.model_report().expect("final posterior");
    println!("recommended_job={} x={:.8} posterior_mean={:.8} posterior_variance={:.8}",
        best.recommended_evaluation, s.jobs()[best.recommended_evaluation].request.x[0], best.mean, best.variance);
    println!("Synthetic logical-time illustration; model estimates, not experimental validation.");
    Ok(())
}
