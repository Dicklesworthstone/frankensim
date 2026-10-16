//! A synthetic common offset/tilt model; the numbers are illustrative, not a
//! calibrated instrument specification. Run with:
//! cargo run -p fs-assimilate --example shared_calibration
use fs_assimilate::groups::RowId;
use fs_assimilate::groups::jacobian::{JacobianBatch, JacobianRecord, SharedNoiseModel};
use fs_assimilate::{AssimError, Belief, point_sensor};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};

fn main() -> Result<(), AssimError> {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    let clock = VirtualClock::new();
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate, arena,
            StreamKey { seed: 6416, kernel_id: 16, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic,
        ).with_time_source(&clock);
        let model = SharedNoiseModel::new(
            "illustrative-offset-tilt-fit",
            vec!["offset".into(), "tilt".into()],
            vec![vec![0.25, 0.0625], vec![0.0625, 0.125]],
            &cx,
        )?;
        let mut rows = Vec::new();
        for (i, (position, reading)) in [(-1.0, 1.8), (0.0, 2.0), (1.0, 2.1)].into_iter().enumerate() {
            rows.push(JacobianRecord::new(
                point_sensor(0, 1, reading, 0.0625, "scan-instrument")?,
                RowId::try_new(format!("sample-{i}"))?,
                &model,
                vec![1.0, position],
            )?);
        }
        let batch = JacobianBatch::new(model, rows, &cx)?;
        let result = batch.assimilate(&Belief::scalar(0.0, 10.0)?, &cx)?;
        let posterior = result.robust().posterior().expect("all readings available");
        println!("total observation covariance: {:?}", batch.covariance());
        println!("posterior mean: {:.8}", posterior.mean()[0]);
        println!("posterior variance: {:.8}", posterior.variance(0)?);
        println!("model-bound result: {}", result.content_id());
        Ok(())
    })
}
