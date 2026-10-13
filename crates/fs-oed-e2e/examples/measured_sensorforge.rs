//! Read actual outcomes from stdin, only after selecting the requested design.
//! Example model is explicitly synthetic and dimensionless; no truth is stored.
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_oed_e2e::plan::acquisition::{AcquisitionError, MeasurementProvider};
use fs_oed_e2e::plan::measurement::*;
use fs_oed_e2e::{AcquisitionCost, ObjectiveValue};
use fs_qty::QtyAny;
use std::io::{self, BufRead, Read, Write};

struct Input<R>(R);
impl<R: BufRead> MeasurementProvider for Input<R> {
    type Error = io::Error;
    fn read(&mut self, request: &MeasurementRequest, cx: &Cx<'_>) -> io::Result<MeasurementOutcome> {
        cx.checkpoint().map_err(|_| io::Error::other("acquisition cancelled"))?;
        eprintln!("measure {} (request {}, variance after {:.8}); enter: alternative value instrument",
            request.alternative(), request.identity(), request.predicted_posterior_variance().value);
        io::stderr().flush()?;
        // Bound line storage. This example uses blocking stdin: real device
        // providers must implement their own cancellable I/O and timeouts.
        let mut bytes = Vec::new();
        let count = self.0.by_ref().take(1025).read_until(b'\n', &mut bytes)?;
        if count == 0 { return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "no reading supplied")); }
        if count > 1024 { return Err(io::Error::new(io::ErrorKind::InvalidData, "reading exceeds 1024 bytes")); }
        let line = std::str::from_utf8(&bytes).map_err(|_| io::Error::other("reading is not UTF-8"))?;
        let mut fields = line.split_whitespace();
        let name = fields.next().ok_or_else(|| io::Error::other("missing alternative"))?;
        let value: f64 = fields.next().ok_or_else(|| io::Error::other("missing value"))?.parse().map_err(|_| io::Error::other("invalid numeric reading"))?;
        let instrument = fields.next().ok_or_else(|| io::Error::other("missing instrument"))?;
        if name != request.alternative() || fields.next().is_some() { return Err(io::Error::other("reading must name the requested alternative and have exactly three fields")); }
        request.reading(ObjectiveValue::dimensionless(value).map_err(io::Error::other)?, instrument).map_err(io::Error::other)
    }
}

fn run(cx: &Cx<'_>) -> Result<(), Box<dyn std::error::Error>> {
    let objective = |x| ObjectiveValue::dimensionless(x);
    let alternatives = ["A", "B"].into_iter().map(|name| {
        MeasurementAlternative::new(name, objective(10.0)?, QtyAny::dimensionless(1.0),
            QtyAny::dimensionless(0.01), AcquisitionCost::relative(1.0)?)
    }).collect::<Result<Vec<_>, MeasurementError>>()?;
    let mut study = MeasurementStudy::new(alternatives, objective(0.03)?, 4, cx)?;
    let input = io::stdin(); let mut provider = Input(input.lock());
    println!("Synthetic dimensionless two-design model; values come only from stdin.");
    loop {
        let request = match study.recommend(cx)? {
            MeasurementDecision::Stop(reason) => {
                println!("stop={reason:?} chosen={} evpi={:.8} accepted={} history={}", study.best_alternative(), study.evpi().value(), study.updates().len(), study.history_id());
                break;
            }
            MeasurementDecision::Acquire(request) => request,
        };
        let acquired = match study.acquire(&request, &mut provider, cx) {
            Ok(acquired) => acquired,
            Err(error) => {
                if let AcquisitionError::After { acquired, .. } = &error {
                    eprintln!("unapplied reading retained: {} {} {}", acquired.request().alternative(), acquired.outcome().value().value(), acquired.outcome().instrument());
                }
                eprintln!("accepted={} history={} (study unchanged by failed acquisition)", study.updates().len(), study.history_id());
                return Err(error.into());
            }
        };
        match acquired.apply(&mut study, cx) {
            Ok(update) => println!("accepted {} mean={:.8} variance={:.8} evpi={:.8}", update.request().alternative(), update.posterior_mean().value(), update.posterior_variance().value, update.evpi().value()),
            Err(error) => {
                eprintln!("unapplied reading: {} {} {}; accepted={} history={}", acquired.request().alternative(), acquired.outcome().value().value(), acquired.outcome().instrument(), study.updates().len(), study.history_id());
                return Err(error.into());
            }
        }
    }
    Ok(())
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    pool.scope(|arena| run(&Cx::new(&gate, arena,
        StreamKey { seed: 20, kernel_id: 4, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
