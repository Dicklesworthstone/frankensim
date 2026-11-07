//! Estimate a cooling rate and shared thermometer offset from one fixed batch.
//! The readings are synthetic. The derivative comes from the declared cooling
//! law; the result is local Gaussian inference, not experimental validation.
use fs_assimilate::nonlinear::iterated::{
    IndependentMeasurements, IteratedSettings, ModelEvaluation, assimilate_iterated,
};
use fs_assimilate::{AssimError, Belief};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};

fn main() -> Result<(), AssimError> {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    let clock = VirtualClock::new();
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 6418,
                kernel_id: 18,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        )
        .with_time_source(&clock);
        let times = [1.0_f64, 2.0, 4.0, 6.0];
        let noise = [0.10, -0.08, 0.06, -0.04];
        let readings = times
            .iter()
            .zip(noise)
            .map(|(&time, error)| 293.15 + 60.0 * (-0.25 * time).exp() + 1.5 + error)
            .collect();
        let instruments = (0..times.len())
            .map(|i| format!("synthetic-reading-{i}"))
            .collect();
        // A single uncertain bias is represented in the physical state. Only
        // the remaining per-reading noise is declared independent.
        let measurements = IndependentMeasurements::new(
            readings,
            vec![0.04; times.len()],
            instruments,
            "newton-cooling-known-ambient-and-initial-temperature",
        )?;
        let prior = Belief::diagonal(vec![0.6, 0.0], &[0.25, 4.0], &cx)?;
        let result = assimilate_iterated(
            &prior,
            &measurements,
            IteratedSettings {
                max_iterations: 32,
                max_model_evaluations: 128,
                max_backtracks: 24,
                gradient_tolerance: 1e-7,
                step_tolerance: 1e-12,
            },
            |state, _| {
                if state[0] <= 0.0 {
                    return Ok(ModelEvaluation::OutsideDomain);
                }
                let decay: Vec<_> = times
                    .iter()
                    .map(|&t| 60.0 * (-state[0] * t).exp())
                    .collect();
                Ok(ModelEvaluation::Value {
                    predictions: decay.iter().map(|&d| 293.15 + d + state[1]).collect(),
                    jacobian: decay
                        .iter()
                        .zip(times)
                        .map(|(&d, t)| vec![-t * d, 1.0])
                        .collect(),
                })
            },
            &cx,
        )?;
        println!(
            "Synthetic cooling-rate and shared-offset calibration; independent residual noise."
        );
        println!(
            "outcome={:?} rate={:.8}/s offset={:.8}K evaluations={} rejected={}",
            result.outcome(),
            result.state()[0],
            result.state()[1],
            result.model_evaluations(),
            result.rejected_trials()
        );
        for (i, row) in result.trace().iter().enumerate() {
            println!(
                "iteration={i} objective={:.10} gradient={:.3e} step={:.6}",
                row.objective, row.gradient_norm, row.step_fraction
            );
        }
        if let Some(posterior) = result.posterior() {
            println!(
                "local variances: rate={:.8}/s^2 offset={:.8}K^2",
                posterior.variance(0)?,
                posterior.variance(1)?
            );
        }
        println!("receipt={}", result.receipt_id());
        Ok(())
    })
}
