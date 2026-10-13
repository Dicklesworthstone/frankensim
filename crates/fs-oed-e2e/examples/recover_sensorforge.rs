//! Separate-process acquisition/recovery using immutable, create-new files.
//! The unitless two-design fixture has NO hidden truth or generated outcomes.
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{self, Read, Write};
#[cfg(unix)]
use std::path::Path;

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_oed_e2e::plan::acquisition::{AcquiredMeasurement, MeasurementProvider};
use fs_oed_e2e::plan::measurement::{MeasurementAlternative, MeasurementDecision, MeasurementError,
    MeasurementOutcome, MeasurementRequest, MeasurementStudy};
use fs_oed_e2e::plan::measurement::checkpoint::MAX_MEASUREMENT_CHECKPOINT_BYTES;
use fs_oed_e2e::{AcquisitionCost, ObjectiveValue};
use fs_qty::QtyAny;

type Result<T> = std::result::Result<T, Box<dyn Error>>;

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    let clock = VirtualClock::new();
    pool.scope(|arena| run(&args, &Cx::new(&gate, arena,
        StreamKey { seed: 920, kernel_id: 20, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic).with_time_source(&clock)))
}

fn run(args: &[String], cx: &Cx<'_>) -> Result<()> {
    match args {
        [command, output, cap] if command == "new" => {
            let alternatives = [("A", 0.0, 4.0), ("B", 0.1, 9.0)].into_iter()
                .map(|(name, mean, variance)| MeasurementAlternative::new(name,
                    ObjectiveValue::dimensionless(mean).expect("fixed finite fixture"),
                    QtyAny::dimensionless(variance), QtyAny::dimensionless(1.0),
                    AcquisitionCost::relative(1.0).expect("fixed positive cost")))
                .collect::<std::result::Result<Vec<_>, _>>()?;
            let study = MeasurementStudy::new(alternatives, ObjectiveValue::dimensionless(0.0)?, cap.parse()?, cx)?;
            write_new(output, &study.to_checkpoint(cx)?)?;
            println!("Saved synthetic unitless priors; no readings have been acquired.");
            describe(&study, cx)?;
        }
        [command, input] if command == "next" => {
            describe(&MeasurementStudy::from_checkpoint(&read_bounded(input)?, cx)?, cx)?;
        }
        [command, input, output, alternative, value, instrument] if command == "record" => {
            let study = MeasurementStudy::from_checkpoint(&read_bounded(input)?, cx)?;
            let request = next(&study, cx)?;
            if request.alternative() != alternative {
                return Err(io::Error::new(io::ErrorKind::InvalidInput, "reading names a different alternative than the current request").into());
            }
            let mut provider = SuppliedReading { value: ObjectiveValue::dimensionless(value.parse()?)?, instrument };
            let acquired = study.acquire(&request, &mut provider, cx)?;
            write_new(output, &acquired.to_checkpoint(cx)?)?;
            println!("Saved unapplied reading for {} ({}); source study unchanged.", request.alternative(), request.identity());
        }
        [command, input, pending, output] if command == "apply" => {
            let mut study = MeasurementStudy::from_checkpoint(&read_bounded(input)?, cx)?;
            let acquired = AcquiredMeasurement::from_checkpoint(&study, &read_bounded(pending)?, cx)?;
            acquired.apply(&mut study, cx)?;
            // Output is a NEW immutable generation. A failed/truncated write
            // cannot replace the source study or the acquired reading.
            write_new(output, &study.to_checkpoint(cx)?)?;
            println!("Saved accepted generation; no provider was called. Source files retained.");
            describe(&study, cx)?;
        }
        _ => return Err(io::Error::new(io::ErrorKind::InvalidInput,
            "usage: recover_sensorforge new STATE MAX_READINGS | next STATE | record STATE PENDING ALTERNATIVE VALUE INSTRUMENT | apply STATE PENDING NEW_STATE").into()),
    }
    Ok(())
}

fn next(study: &MeasurementStudy, cx: &Cx<'_>) -> Result<MeasurementRequest> {
    match study.recommend(cx)? {
        MeasurementDecision::Acquire(request) => Ok(request),
        MeasurementDecision::Stop(stop) => Err(io::Error::new(io::ErrorKind::InvalidInput, format!("study stopped: {stop:?}")).into()),
    }
}
fn describe(study: &MeasurementStudy, cx: &Cx<'_>) -> Result<()> {
    println!("accepted={} best={} modeled_evpi={} history={}", study.updates().len(),
        study.best_alternative(), study.evpi().value(), study.history_id());
    match study.recommend(cx)? {
        MeasurementDecision::Acquire(request) => println!("next={} ordinal={} request={}", request.alternative(), request.ordinal(), request.identity()),
        MeasurementDecision::Stop(stop) => println!("stop={stop:?}"),
    }
    Ok(())
}
struct SuppliedReading<'a> { value: ObjectiveValue, instrument: &'a str }
impl MeasurementProvider for SuppliedReading<'_> {
    type Error = MeasurementError;
    fn read(&mut self, request: &MeasurementRequest, _: &Cx<'_>) -> std::result::Result<MeasurementOutcome, Self::Error> {
        request.reading(self.value, self.instrument)
    }
}
fn read_bounded(path: &str) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    File::open(path)?.take((MAX_MEASUREMENT_CHECKPOINT_BYTES + 1) as u64).read_to_end(&mut bytes)?;
    if bytes.len() > MAX_MEASUREMENT_CHECKPOINT_BYTES {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "checkpoint exceeds byte limit").into());
    }
    Ok(bytes)
}
fn write_new(path: &str, bytes: &[u8]) -> Result<()> {
    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;
    file.write_all(bytes)?;
    file.sync_all()?;
    #[cfg(unix)]
    {
        let parent = Path::new(path).parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or_else(|| Path::new("."));
        File::open(parent)?.sync_all()?;
    }
    Ok(())
}
