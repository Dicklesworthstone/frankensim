//! Synthetic streaming readings followed by outcome-free sensor selection.
use fs_assimilate::groups::RowId;
use fs_assimilate::groups::jacobian::{JacobianBatch, SharedNoiseModel};
use fs_assimilate::groups::stream::SharedNoiseStream;
use fs_assimilate::groups::stream::design::SharedNoiseProbe;
use fs_assimilate::{AssimError, Belief};
use fs_assimilate::factor::{FactorBelief, forecast::LinearDynamics};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};

fn main() -> Result<(), AssimError> {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    let clock = VirtualClock::new();
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 6416, kernel_id: 16, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic).with_time_source(&clock);
        let model = SharedNoiseModel::new("fixed-calibration", vec!["offset".into()], vec![vec![4.0]], &cx)?;
        let mut stream = SharedNoiseStream::new(&Belief::scalar(0.0, 10.0)?, &model, 64, &cx)?;
        // These are synthetic acquired values, not oracle truth used by planning.
        for (index, value) in [1.1, 0.9, 1.2, 1.0, 0.8, 1.1, 1.0, 0.9].into_iter().enumerate() {
            let sensor = SharedNoiseProbe::new(&model, RowId::try_new(format!("sample-{index}"))?, vec![1.0], vec![1.0], 0.25)?;
            let batch = JacobianBatch::new(model.clone(), vec![sensor.record(&model, value, "scan-instrument")?], &cx)?;
            stream.assimilate(&batch, &cx)?;
        }
        let candidates = [
            SharedNoiseProbe::new(&model, RowId::try_new("repeat-direction")?, vec![1.0], vec![1.0], 0.25)?,
            SharedNoiseProbe::new(&model, RowId::try_new("opposite-direction")?, vec![1.0], vec![-1.0], 0.25)?,
        ];
        let scores = stream.rank_probes(&[1.0], &candidates, &cx)?;
        for score in &scores {
            println!("{}: predicted target variance {:.8}, reduction {:.8}", score.row().as_str(), score.posterior_target_variance(), score.variance_reduction());
        }
        let chosen = candidates.iter().find(|probe| probe.row() == scores[0].row()).expect("ranked probe exists");
        // Acquisition is a separate operation and explicitly supplies a reading.
        let batch = JacobianBatch::new(model.clone(), vec![chosen.record(&model, 1.05, "next-instrument")?], &cx)?;
        stream.assimilate(&batch, &cx)?;
        let physical = stream.physical_belief(&cx)?;
        println!("physical mean {:.8}; variance {:.8}; accepted rows {}", physical.mean()[0], physical.variance(0)?, stream.accepted_rows());
        let replay = stream.assimilate(&batch, &cx)?;
        assert!(replay.accepted_indices().is_empty());
        println!("exact replay ignored; retained history {}", stream.history_id());
        // A complete synthetic one-second transition, with independent process
        // variance. The learned calibration state and cross-covariance survive.
        let dynamics = LinearDynamics::new(vec![vec![0.9]], vec![0.2],
            FactorBelief::diagonal(vec![0.0], &[0.05], &cx)?, "synthetic-one-second", &cx)?;
        stream.forecast_physical(&dynamics, 0.0, 1.0, &cx)?;
        let forecast = stream.physical_belief(&cx)?;
        println!("forecast t=1: physical mean {:.8}; variance {:.8}; calibration mean {:.8}",
            forecast.mean()[0], forecast.variance(0)?, stream.joint().mean()[1]);
        Ok(())
    })
}
