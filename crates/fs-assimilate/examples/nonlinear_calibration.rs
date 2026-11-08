//! Estimate a cooling rate and shared thermometer offset from one fixed batch.
//! Compare explicit offset estimation with marginalizing the same uncertain
//! offset into a full observation covariance; both retain its shared error.
//! The readings are synthetic. The derivative comes from the declared cooling
//! law; the result is local Gaussian inference, not experimental validation.
use fs_assimilate::nonlinear::iterated::{
    CorrelatedMeasurements, IndependentMeasurements, IteratedSettings, ModelEvaluation,
    assimilate_correlated_iterated, assimilate_iterated,
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
        let readings: Vec<f64> = times
            .iter()
            .zip(noise)
            .map(|(&time, error)| 293.15 + 60.0 * (-0.25 * time).exp() + 1.5 + error)
            .collect();
        let instruments: Vec<String> = (0..times.len())
            .map(|i| format!("synthetic-reading-{i}"))
            .collect();
        // A single uncertain bias is represented in the physical state. Only
        // the remaining per-reading noise is declared independent.
        let measurements = IndependentMeasurements::new(
            readings.clone(),
            vec![0.04; times.len()],
            instruments.clone(),
            "newton-cooling-known-ambient-and-initial-temperature",
        )?;
        let prior = Belief::diagonal(vec![0.6, 0.0], &[0.25, 4.0], &cx)?;
        let settings = IteratedSettings {
            max_iterations: 32,
            max_model_evaluations: 128,
            max_backtracks: 24,
            gradient_tolerance: 1e-7,
            step_tolerance: 1e-12,
        };
        let result = assimilate_iterated(
            &prior,
            &measurements,
            settings,
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

        // Integrating out the same N(0,4) common offset yields R=0.04 I+4 11'.
        // Its variance belongs in every off-diagonal as well as the diagonals.
        let mut covariance = vec![vec![4.0; times.len()]; times.len()];
        for (i, row) in covariance.iter_mut().enumerate() {
            row[i] += 0.04;
        }
        let correlated = CorrelatedMeasurements::new(
            readings,
            covariance,
            instruments,
            "newton-cooling-marginalized-shared-offset",
            &cx,
        )?;
        let marginalized = assimilate_correlated_iterated(
            &Belief::scalar(0.6, 0.25)?,
            &correlated,
            settings,
            |state, _| {
                if state[0] <= 0.0 {
                    return Ok(ModelEvaluation::OutsideDomain);
                }
                let decay: Vec<_> = times
                    .iter()
                    .map(|&t| 60.0 * (-state[0] * t).exp())
                    .collect();
                Ok(ModelEvaluation::Value {
                    predictions: decay.iter().map(|&d| 293.15 + d).collect(),
                    jacobian: decay
                        .iter()
                        .zip(times)
                        .map(|(&d, t)| vec![-t * d])
                        .collect(),
                })
            },
            &cx,
        )?;
        println!(
            "Marginalized shared offset: outcome={:?} rate={:.8}/s evaluations={}",
            marginalized.outcome(),
            marginalized.state()[0],
            marginalized.model_evaluations(),
        );
        if let Some(posterior) = marginalized.posterior() {
            println!("local rate variance={:.8}/s^2", posterior.variance(0)?);
        }
        println!("receipt={}", marginalized.receipt_id());
        Ok(())
    })
}
