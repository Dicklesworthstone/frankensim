//! Synthetic thermocouple tracking with uncertain cooling rate and shared bias.
//! The nonlinear Newton-cooling step is an exact discrete model here; inference
//! is a local Gaussian approximation, not experimental validation.
use fs_assimilate::factor::{FactorBelief, forecast::NonlinearDynamics};
use fs_assimilate::groups::RowId;
use fs_assimilate::groups::jacobian::{JacobianBatch, SharedNoiseModel};
use fs_assimilate::groups::stream::{SharedNoiseStream, design::SharedNoiseProbe};
use fs_assimilate::{AssimError, Belief};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};

fn main() -> Result<(), AssimError> {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    let clock = VirtualClock::new();
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate, arena, StreamKey { seed: 6417, kernel_id: 17, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic,
        ).with_time_source(&clock);
        // Physical state: temperature in kelvin, cooling rate in 1/second.
        // The thermometer's additive kelvin bias is retained separately.
        let prior = Belief::diagonal(vec![350.0, 0.25], &[4.0, 0.01], &cx)?;
        let calibration = SharedNoiseModel::new(
            "synthetic-thermometer", vec!["temperature-bias".into()], vec![vec![1.0]], &cx,
        )?;
        let mut stream = SharedNoiseStream::new(&prior, &calibration, 16, &cx)?;
        let dynamics = NonlinearDynamics::new(
            FactorBelief::diagonal(vec![0.0, 0.0], &[0.01, 0.0], &cx)?,
            vec![0.001, 0.00001], "newton-cooling-ambient-293.15K-step-1s", &cx,
        )?;
        println!("Synthetic readings; local Gaussian state estimation with a retained shared bias.");
        for (ordinal, reading) in [338.2, 328.5, 320.8, 314.8, 310.1].into_iter().enumerate() {
            let from = ordinal as f64;
            stream.forecast_physical_nonlinear(&dynamics, from, from + 1.0, |state, _| {
                if state[1] < 0.0 {
                    return Err(AssimError::InvalidForecast { reason: "cooling rate must be nonnegative" });
                }
                Ok(vec![293.15 + (state[0] - 293.15) * (-state[1]).exp(), state[1]])
            }, &cx)?;
            let sensor = SharedNoiseProbe::new(
                &calibration, RowId::try_new(format!("thermocouple-{ordinal}"))?,
                vec![1.0, 0.0], vec![1.0], 0.09,
            )?;
            let batch = JacobianBatch::new(
                calibration.clone(),
                vec![sensor.record(&calibration, reading, "synthetic-thermocouple")?], &cx,
            )?;
            stream.assimilate(&batch, &cx)?;
            let physical = stream.physical_belief(&cx)?;
            println!(
                "t={:.0}s temperature={:.5}K variance={:.6}K^2 cooling-rate={:.6}/s bias={:.5}K",
                from + 1.0, physical.mean()[0], physical.variance(0)?,
                physical.mean()[1], stream.joint().mean()[2],
            );
        }
        Ok(())
    })
}
