//! Run an explicitly seeded synthetic study through production sensing APIs.
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_oed_e2e::plan::acquisition::simulation::{GaussianMeasurementSimulator, SimulationTruth};
use fs_oed_e2e::plan::acquisition::simulation::study::SimulatedStudy;
use fs_oed_e2e::plan::measurement::{MeasurementAlternative, MeasurementError};
use fs_oed_e2e::{AcquisitionCost, ObjectiveValue};
use fs_qty::QtyAny;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let usage = "usage: simulated_sensorforge SEED MAX_READINGS (0..=128)";
    let seed: u64 = args.next().ok_or(usage)?.parse()?;
    let cap: usize = args.next().ok_or(usage)?.parse()?;
    if cap > 128 || args.next().is_some() { return Err(usage.into()); }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| -> Result<(), MeasurementError> {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed, kernel_id: 20, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic).with_time_source(&clock);
        let mut alternatives = Vec::new();
        for (name, mean) in [("A", 10.0), ("B", 10.2)] {
            alternatives.push(MeasurementAlternative::new(name, ObjectiveValue::dimensionless(mean)?,
                QtyAny::dimensionless(1.0), QtyAny::dimensionless(0.25), AcquisitionCost::relative(1.0)?)?);
        }
        let provider = GaussianMeasurementSimulator::new(seed, vec![
            SimulationTruth::new("A", ObjectiveValue::dimensionless(10.1)?)?,
            SimulationTruth::new("B", ObjectiveValue::dimensionless(10.2)?)?,
        ], &cx)?;
        let mut run = SimulatedStudy::new(alternatives, ObjectiveValue::dimensionless(0.025)?, cap, provider, &cx)?;
        println!("synthetic dimensionless model; seed={seed} max_readings={cap} source={}", run.source_identity());
        // One extra recommendation reports the actual stop after the reading cap.
        let progress = run.advance(cap + 1, &cx);
        // Retain and report accepted data even when a later operation refuses.
        for update in run.study().updates() {
            println!("{} {} reading={:.9} posterior_mean={:.9} posterior_variance={:.9}",
                update.request().ordinal(), update.request().alternative(), update.outcome().value().value(),
                update.posterior_mean().value(), update.posterior_variance().value);
        }
        let progress = progress?;
        println!("stop={:?} best={} estimated_evpi={:.9} synthetic_realized_loss={:.9}",
            progress.stop(), run.study().best_alternative(), run.study().evpi().value(), run.realized_loss()?.value());
        println!("Gaussian-model estimates only; no calibrated risk or optional-stopping certificate.");
        Ok(())
    });
    assert!(pool.stats().quiescent());
    result?;
    Ok(())
}
