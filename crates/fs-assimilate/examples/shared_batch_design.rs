//! Joint sensor planning followed by a separate acquired-reading update.
//! The small unitless Gaussian fixture is synthetic, not calibrated metrology.
//! Optional arguments provide values for rows `a` and `c`; selection never
//! receives these values. With no arguments, explicitly labeled synthetic
//! readings demonstrate the same public acquisition API.
use fs_assimilate::groups::{RowId, jacobian::SharedNoiseModel, stream::SharedNoiseStream};
use fs_assimilate::groups::stream::design::SharedNoiseProbe;
use fs_assimilate::groups::stream::design::batch::{CostedProbe, ProbeBatchBudget, ProbeReading};
use fs_assimilate::Belief;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).take(3).collect();
    let (values, label) = match args.as_slice() {
        [] => ([3.0, 1.0], "synthetic demonstration readings"),
        [a, c] => ([a.parse::<f64>()?, c.parse::<f64>()?], "caller-supplied readings"),
        _ => return Err("usage: shared_batch_design [reading_a reading_c]".into()),
    };
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| -> Result<(), Box<dyn std::error::Error>> {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 6420, kernel_id: 20, tile: 0, iteration: 0 },
            Budget::INFINITE.with_cost_quota(10_000_000).with_poll_quota(10_000),
            ExecMode::Deterministic).with_time_source(&clock);
        let model = SharedNoiseModel::new("zero-reference-error", vec!["offset".into()], vec![vec![0.0]], &cx)?;
        let prior = Belief::diagonal(vec![0.0, 0.0], &[2.0, 1.0], &cx)?;
        let mut stream = SharedNoiseStream::new(&prior, &model, 8, &cx)?;
        let mut candidates = Vec::new();
        for (row, response) in [("a", vec![1.0, 0.0]), ("b", vec![1.0, 0.0]), ("c", vec![0.0, 1.0])] {
            candidates.push(CostedProbe::new(SharedNoiseProbe::new(&model, RowId::try_new(row)?, response, vec![0.0], 0.01)?, 1)?);
        }
        // Unit of acquisition cost: one equal-effort reading. Numerical compute
        // limits above are separate, and no real-world units are inferred.
        let target = [1.0, 1.0];
        let plan = stream.plan_probe_batch(&target, &candidates, ProbeBatchBudget::new(2, 2)?, &cx)?;
        println!("selection stop: {:?}; acquisition cost: {}", plan.stop(), plan.total_cost());
        for step in plan.steps() {
            println!("row {}: target variance {:.12} -> {:.12}", step.candidate().probe().row().as_str(),
                step.prior_target_variance(), step.posterior_target_variance());
        }
        plan.check_current(&stream)?;
        // Acquisition is intentionally separate and supplied in reverse order.
        // No reading or oracle truth was passed to plan_probe_batch.
        let readings = vec![
            ProbeReading::new(RowId::try_new("c")?, values[1], "c-instrument")?,
            ProbeReading::new(RowId::try_new("a")?, values[0], "a-instrument")?,
        ];
        println!("input: {label}");
        let update = plan.assimilate_readings(&mut stream, &model, &readings, &cx)?;
        let posterior = stream.physical_belief(&cx)?;
        let actual_variance: f64 = posterior.covariance().iter().flatten().sum();
        println!("accepted rows: {}; posterior mean: {:?}", update.accepted_indices().len(), posterior.mean());
        println!("predicted target variance: {:.12}; actual: {:.12}", plan.posterior_target_variance(), actual_variance);
        if (actual_variance - plan.posterior_target_variance()).abs() > 1.0e-10 {
            return Err("synthetic forecast/acquisition comparison failed".into());
        }
        Ok(())
    });
    assert!(pool.stats().quiescent());
    result
}
