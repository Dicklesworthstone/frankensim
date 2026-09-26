//! Native project uncertainty: every observation is a completed ordinary
//! geometry-import/solve/QoI pipeline. The statistical executor owns sampling
//! and recovery; this adapter binds its observations to retained child runs.

use std::fmt::{self, Write as _};
use std::io::Read as _;
use std::path::Path;
use std::time::Instant;

use fs_blake3::{ContentHash, hash_bytes};
use fs_exec::CancelGate;
use fs_ledger::{EdgeRole, FiveExplicits, Ledger, LedgerError, OpOutcome};
use fs_package::{Claim, EvidencePackage, Provenance};
use fs_uq::{
    CorrelationModel, ParameterUncertainty, PropagationMethod, UqExecution, UqPlan, UqStatus,
};

use super::STUDY_RUN_RECEIPT_SCHEMA;
use crate::json_read::JsonValue as J;
use crate::{CommandOutput, Diagnostic, OutputMode, exit, push_json_string, refusal};

#[path = "study_uncertainty/legacy.rs"]
mod legacy;
#[path = "study_uncertainty/model.rs"]
mod model;
use model::{Model, Sample};

const DRIVER: &str = "native-cooling-uncertainty-v1";
const RECEIPT_KIND: &str = "study-run-receipt";
const REPORT_SCHEMA: &str = "frankensim.cli.native-uncertainty-result.v1";
const MAX_ARTIFACT_BYTES: u64 = 16 * 1024 * 1024;
const NO_CLAIM: &str = "Empirical propagation through the declared native numerical cooling model with explicitly independent uniform inputs. Fixed-count means, spread, quantiles and pass fractions are descriptive estimates, not confidence intervals or optional-stopping decisions. Child engineering uncertainty budgets and verdicts remain unchanged. No continuum, physical-model, experimental-validation or safety-signoff claim; source/card tolerances are not probability distributions. Interrupted or refused samples are never replaced, clipped or skipped.";

type Result<T> = std::result::Result<T, Failure>;
#[derive(Debug)]
struct Failure {
    code: &'static str,
    message: String,
    exit: u8,
}
impl fmt::Display for Failure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.code, self.message)
    }
}
impl From<LedgerError> for Failure {
    fn from(error: LedgerError) -> Self {
        fail("cli-uncertainty-ledger", error.to_string())
    }
}
fn fail(code: &'static str, message: impl Into<String>) -> Failure {
    Failure {
        code,
        message: message.into(),
        exit: exit::REFUSED,
    }
}
fn quoted(value: &str) -> String {
    let mut text = String::new();
    push_json_string(&mut text, value);
    text
}
fn optional(value: Option<f64>) -> String {
    value.map_or_else(|| "null".into(), |v| v.to_string())
}
fn output_error(command: &'static str, mode: OutputMode, e: Failure) -> CommandOutput {
    refusal(
        mode,
        e.exit,
        &Diagnostic::new(
            command,
            e.code,
            e.message,
            "inspect the named native study input or retained run; each sample must complete the ordinary cooling solve",
        ),
        None,
    )
}
fn artifact(ledger: &Ledger, hash: ContentHash, kind: &str, cap: u64) -> Result<Vec<u8>> {
    let info = ledger
        .artifact_info(&hash)?
        .ok_or_else(|| fail("cli-uncertainty-artifact", "missing retained artifact"))?;
    if info.kind != kind {
        return Err(fail(
            "cli-uncertainty-artifact",
            format!("expected {kind}, found {}", info.kind),
        ));
    }
    let bytes = ledger
        .get_artifact_bounded(&hash, cap)?
        .ok_or_else(|| fail("cli-uncertainty-artifact", "missing retained bytes"))?;
    if hash_bytes(&bytes) != hash {
        return Err(fail(
            "cli-uncertainty-artifact",
            "retained content hash differs",
        ));
    }
    Ok(bytes)
}
fn parse(bytes: &[u8]) -> Result<J> {
    let text =
        std::str::from_utf8(bytes).map_err(|e| fail("cli-uncertainty-receipt", e.to_string()))?;
    J::parse(text).map_err(|e| fail("cli-uncertainty-receipt", e.to_string()))
}
fn integer(value: &J, key: &str) -> Result<usize> {
    value
        .get(key)
        .and_then(J::number_raw)
        .and_then(|v| v.parse().ok())
        .ok_or_else(|| fail("cli-uncertainty-receipt", format!("invalid {key}")))
}
fn hash_field(value: &J, key: &str) -> Result<ContentHash> {
    value
        .str_field(key)
        .and_then(ContentHash::from_hex)
        .ok_or_else(|| fail("cli-uncertainty-receipt", format!("missing {key} hash")))
}
fn linked(ledger: &Ledger, value: &J, key: &str, kind: &str) -> Result<Vec<u8>> {
    artifact(ledger, hash_field(value, key)?, kind, MAX_ARTIFACT_BYTES)
}
fn budget(text: Option<&str>) -> Result<Option<usize>> {
    text.map(|text| {
        text.parse::<usize>()
            .ok()
            .filter(|n| *n <= fs_project::uncertainty::MAX_SAMPLES)
            .ok_or_else(|| {
                fail(
                    "cli-uncertainty-budget",
                    "--budget must be an additional sample count in 0..=256",
                )
            })
    })
    .transpose()
}
fn plan(model: &Model) -> UqPlan {
    let study = model.bound.study();
    let mut plan = UqPlan::new(study.qoi(), PropagationMethod::MonteCarlo, study.samples())
        .with_correlation(CorrelationModel::Independent)
        .with_compliance_threshold(model.bound.threshold_k());
    plan.seed = study.seed();
    for p in study.parameters() {
        plan = plan.with_parameter(ParameterUncertainty::uniform(
            &p.name,
            p.low,
            p.high,
            p.target.unit(),
        ));
    }
    plan
}
fn rows_json(rows: &[Sample]) -> String {
    let mut text = String::from("[");
    for (ordinal, row) in rows.iter().enumerate() {
        if ordinal > 0 {
            text.push(',');
        }
        let parameters = row
            .parameters
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(",");
        let _ = write!(
            text,
            "{{\"ordinal\":{ordinal},\"parameters\":[{parameters}],\"run\":{},\"project_hash\":{},\"qoi_receipt\":{},\"value_k\":{}}}",
            quoted(&row.run),
            quoted(&row.project_hash),
            quoted(&row.qoi_receipt.to_hex()),
            row.value_k
        );
    }
    text.push(']');
    text
}
fn read_rows(bytes: &[u8]) -> Result<Vec<Sample>> {
    let value = parse(bytes)?;
    let rows = value
        .as_array()
        .filter(|r| r.len() <= fs_project::uncertainty::MAX_SAMPLES)
        .ok_or_else(|| {
            fail(
                "cli-uncertainty-observations",
                "expected a bounded observation array",
            )
        })?;
    rows.iter()
        .enumerate()
        .map(|(ordinal, row)| {
            if integer(row, "ordinal")? != ordinal {
                return Err(fail(
                    "cli-uncertainty-observations",
                    "sample ordinals are not contiguous",
                ));
            }
            let parameters = row
                .get("parameters")
                .and_then(J::as_array)
                .filter(|p| p.len() <= 32)
                .ok_or_else(|| fail("cli-uncertainty-observations", "invalid sample parameters"))?
                .iter()
                .map(|p| {
                    p.as_f64().filter(|v| v.is_finite()).ok_or_else(|| {
                        fail("cli-uncertainty-observations", "nonfinite sample parameter")
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            Ok(Sample {
                parameters,
                run: row
                    .str_field("run")
                    .ok_or_else(|| fail("cli-uncertainty-observations", "missing child run"))?
                    .into(),
                project_hash: row
                    .str_field("project_hash")
                    .ok_or_else(|| {
                        fail("cli-uncertainty-observations", "missing child project hash")
                    })?
                    .into(),
                qoi_receipt: hash_field(row, "qoi_receipt")?,
                value_k: row
                    .f64_field("value_k")
                    .filter(|v| v.is_finite())
                    .ok_or_else(|| fail("cli-uncertainty-observations", "nonfinite sample QoI"))?,
            })
        })
        .collect()
}

struct Outcome {
    pointer: String,
    receipt: String,
    status: String,
}
fn render(out: Outcome, mode: OutputMode) -> CommandOutput {
    let exit_code = match out.status.as_str() {
        "completed" => exit::SUCCESS,
        "refused" => exit::REFUSED,
        "cancelled" => exit::CANCELLED,
        _ => exit::BUDGET,
    };
    let stdout = match mode {
        OutputMode::Json => format!(
            "{{\"command\":\"study\",\"status\":{},\"run_id\":{},\"run\":{},\"receipt\":{}}}\n",
            quoted(&out.status),
            quoted(&out.pointer),
            quoted(&out.pointer),
            out.receipt
        ),
        OutputMode::Text => format!(
            "command=study\nstatus={}\nrun={}\nauthority=estimated-native-model-sampling\n",
            out.status, out.pointer
        ),
    };
    CommandOutput {
        exit_code,
        stdout,
        stderr: String::new(),
    }
}
fn op_ir(id: ContentHash, n: usize) -> String {
    format!(
        "{{\"driver\":{DRIVER:?},\"model\":{},\"samples_completed\":{n}}}",
        quoted(&id.to_hex())
    )
}

#[allow(clippy::too_many_arguments)]
fn persist(
    model: &Model,
    ledger: &Ledger,
    execution: &UqExecution,
    rows: &[Sample],
    status: &str,
    termination: &str,
    used_wall: f64,
    predecessor: Option<ContentHash>,
) -> Result<Outcome> {
    let n = rows.len();
    if n != execution.observations().len()
        || rows
            .iter()
            .zip(execution.observations())
            .any(|(r, v)| r.value_k.to_bits() != v.to_bits())
    {
        return Err(fail(
            "cli-uncertainty-observations",
            "retained child runs disagree with statistical observations",
        ));
    }
    let id = model.identity();
    let report = execution.report();
    let samples = rows_json(rows);
    let statistics = if status == "completed" {
        let quantiles = report.percentiles.map_or_else(
            || "null".into(),
            |v| format!("[{},{},{}]", v[0], v[1], v[2]),
        );
        format!(
            "{{\"mean_k\":{},\"std_dev_k\":{},\"sampling_standard_error_k\":{},\"empirical_probability_of_compliance\":{},\"empirical_range_k\":[{},{}],\"quantiles_p05_p50_p95_k\":{quantiles}}}",
            optional(report.mean),
            optional(report.std_dev),
            optional(report.std_dev.map(|_| report.sampling_error)),
            optional(report.probability_of_compliance),
            report.interval_bounds[0],
            report.interval_bounds[1]
        )
    } else {
        "null".into()
    };
    let failure = report
        .rejection_reason
        .as_deref()
        .map_or_else(|| "null".into(), quoted);
    let summary = format!(
        "{{\"schema\":{REPORT_SCHEMA:?},\"driver\":{DRIVER:?},\"study_id\":{},\"status\":{},\"termination\":{},\"qoi\":\"temperature-max\",\"unit\":\"K\",\"authority\":\"Estimated\",\"method\":\"monte-carlo\",\"correlation\":\"independent\",\"seed\":{},\"samples_evaluated\":{n},\"samples_planned\":{},\"evaluations_attempted\":{},\"temperature_limit_k\":{},\"observations\":{samples},\"statistics\":{statistics},\"failure\":{failure},\"no_claim\":{}}}",
        quoted(&id.to_hex()),
        quoted(status),
        quoted(termination),
        quoted(&model.bound.study().seed().to_string()),
        model.bound.study().samples(),
        execution.evaluations_attempted(),
        model.bound.threshold_k(),
        quoted(NO_CLAIM)
    );
    let mut table = String::new();
    for (i, row) in rows.iter().enumerate() {
        let _ = write!(
            table,
            "<tr><td>{}</td><td>{}</td><td><code>{}</code></td></tr>",
            i + 1,
            row.value_k,
            row.run
        );
    }
    let html = format!(
        "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Native cooling uncertainty</title><body><h1>Native cooling uncertainty</h1><p>Status: {status}; {n}/{} completed samples. Estimated.</p><p>Mean: {} K; sample standard deviation: {} K; empirical pass fraction: {}.</p><p>{NO_CLAIM}</p><table><tr><th>Sample</th><th>Maximum temperature (K)</th><th>Retained solve</th></tr>{table}</table></body></html>",
        model.bound.study().samples(),
        if status == "completed" {
            optional(report.mean)
        } else {
            "unavailable".into()
        },
        if status == "completed" {
            optional(report.std_dev)
        } else {
            "unavailable".into()
        },
        if status == "completed" {
            optional(report.probability_of_compliance)
        } else {
            "unavailable".into()
        }
    );
    let mut package = EvidencePackage::new(Provenance::new(
        format!("fs-cli/{}+{DRIVER}", env!("CARGO_PKG_VERSION")),
        id.to_hex(),
    ));
    if status == "completed" {
        package = package.with_claim(Claim::estimated("cooling.uncertainty.sample-mean",
            format!("{} K across {n} completed native solves; descriptive Monte Carlo standard error {} K. Result {}. {NO_CLAIM}",
                optional(report.mean), report.sampling_error, hash_bytes(summary.as_bytes()).to_hex()),
            "fixed-count-native-monte-carlo-descriptive-standard-error", report.sampling_error));
    }
    let package = package
        .to_json()
        .map_err(|e| fail("cli-uncertainty-package", e.to_string()))?;
    let checkpoint = if report.status == UqStatus::Refused {
        None
    } else {
        Some(
            execution
                .checkpoint(id)
                .map_err(|e| fail("cli-uncertainty-checkpoint", e.to_string()))?,
        )
    };
    let seed = model.bound.study().seed().to_le_bytes();
    let versions = format!(
        "{{\"driver\":{DRIVER:?},\"solve_driver\":{}}}",
        crate::SOLVE_DRIVER_VERSION
    );
    let budgets = format!(
        "{{\"samples\":{},\"wall_s\":{},\"consumed_wall_s\":{used_wall}}}",
        model.bound.study().samples(),
        model.bound.study().wall_seconds()
    );
    ledger.begin()?;
    let result = (|| {
        let op = ledger.begin_op(Some(id.as_bytes()), &op_ir(id, n), &FiveExplicits {
            seed: &seed, versions: &versions, budget: &budgets,
            capability: "{\"ops\":[\"native-project-uncertainty\"],\"physical_ops\":\"inherited-from-retained-base-project\"}",
        }, 0)?;
        let model_hash = model.retain(ledger, op)?;
        if model_hash != id {
            return Err(fail(
                "cli-uncertainty-model",
                "retained model identity differs",
            ));
        }
        if let Some(previous) = predecessor {
            ledger.link(op, &previous, EdgeRole::In)?;
        }
        let mut child_receipts = std::collections::BTreeSet::new();
        for row in rows {
            if child_receipts.insert(row.qoi_receipt) {
                ledger.link(op, &row.qoi_receipt, EdgeRole::In)?;
            }
        }
        let mut refs = String::new();
        for (key, kind, bytes) in [
            (
                "observations",
                "native-uncertainty-observations",
                samples.as_bytes(),
            ),
            ("report_json", "study-report-json", summary.as_bytes()),
            ("report_html", "study-report-html", html.as_bytes()),
            ("package", "study-package", package.as_bytes()),
        ] {
            let a = ledger.put_artifact(kind, bytes, None)?;
            ledger.link(op, &a.hash, EdgeRole::Out)?;
            let _ = write!(refs, ",{key:?}:{}", quoted(&a.hash.to_hex()));
        }
        let checkpoint_hash = if let Some(bytes) = checkpoint {
            let a = ledger.put_artifact("native-uncertainty-checkpoint", &bytes, None)?;
            ledger.link(op, &a.hash, EdgeRole::Out)?;
            quoted(&a.hash.to_hex())
        } else {
            "null".into()
        };
        let receipt = format!(
            "{{\"schema\":{STUDY_RUN_RECEIPT_SCHEMA:?},\"driver\":{DRIVER:?},\"study_id\":{},\"model\":{},\"status\":{},\"termination\":{},\"samples_completed\":{n},\"samples_planned\":{},\"consumed_wall_s\":{used_wall},\"failure\":{failure},\"checkpoint\":{checkpoint_hash},\"predecessor\":{}{refs}}}",
            quoted(&id.to_hex()),
            quoted(&id.to_hex()),
            quoted(status),
            quoted(termination),
            model.bound.study().samples(),
            predecessor.map_or_else(|| "null".into(), |h| quoted(&h.to_hex()))
        );
        let a = ledger.put_artifact(RECEIPT_KIND, receipt.as_bytes(), None)?;
        ledger.link(op, &a.hash, EdgeRole::Out)?;
        if ledger.artifact_output_seal(&a.hash)?.is_none() {
            ledger.seal_artifact_output(&a.hash, op)?;
        }
        ledger.finish_op(op, OpOutcome::Ok, None, 1)?;
        Ok(Outcome {
            pointer: format!("study-{}", a.hash.to_hex()),
            receipt,
            status: status.into(),
        })
    })();
    match result {
        Ok(out) => match ledger.commit() {
            Ok(()) => Ok(out),
            Err(e) => {
                ledger.rollback()?;
                Err(e.into())
            }
        },
        Err(e) => {
            ledger.rollback()?;
            Err(e)
        }
    }
}

struct Loaded {
    hash: ContentHash,
    value: J,
    bytes: Vec<u8>,
}
fn load(ledger: &Ledger, pointer: &str) -> Result<Loaded> {
    let hash = pointer
        .strip_prefix("study-")
        .and_then(ContentHash::from_hex)
        .ok_or_else(|| {
            fail(
                "cli-uncertainty-run",
                "expected study- followed by a receipt hash",
            )
        })?;
    let bytes = artifact(ledger, hash, RECEIPT_KIND, MAX_ARTIFACT_BYTES)?;
    let value = parse(&bytes)?;
    if value.str_field("schema") != Some(STUDY_RUN_RECEIPT_SCHEMA)
        || value.str_field("driver") != Some(DRIVER)
    {
        return Err(fail(
            "cli-uncertainty-run",
            "unsupported native uncertainty receipt",
        ));
    }
    let id = hash_field(&value, "study_id")?;
    let n = integer(&value, "samples_completed")?;
    let op_id = ledger
        .artifact_output_seal(&hash)?
        .ok_or_else(|| fail("cli-uncertainty-run", "unsealed study receipt"))?;
    let op = ledger
        .op(op_id)?
        .ok_or_else(|| fail("cli-uncertainty-run", "missing study producer"))?;
    if hash_field(&value, "model")? != id
        || n > fs_project::uncertainty::MAX_SAMPLES
        || op.session.as_deref() != Some(id.as_bytes().as_slice())
        || op.ir != op_ir(id, n)
        || op.outcome.as_deref() != Some("ok")
        || !ledger.edge_exists(op_id, &hash, EdgeRole::Out)?
    {
        return Err(fail(
            "cli-uncertainty-run",
            "receipt differs from its completed producing operation",
        ));
    }
    for (key, kind, role) in [
        ("model", model::MANIFEST_KIND, EdgeRole::In),
        (
            "observations",
            "native-uncertainty-observations",
            EdgeRole::Out,
        ),
        ("report_json", "study-report-json", EdgeRole::Out),
        ("report_html", "study-report-html", EdgeRole::Out),
        ("package", "study-package", EdgeRole::Out),
    ] {
        let h = hash_field(&value, key)?;
        if !ledger.edge_exists(op_id, &h, role)? {
            return Err(fail(
                "cli-uncertainty-run",
                format!("missing {key} lineage"),
            ));
        }
        artifact(ledger, h, kind, MAX_ARTIFACT_BYTES)?;
    }
    if value.str_field("status") != Some("refused") {
        let h = hash_field(&value, "checkpoint")?;
        if !ledger.edge_exists(op_id, &h, EdgeRole::Out)? {
            return Err(fail("cli-uncertainty-run", "missing checkpoint lineage"));
        }
        artifact(
            ledger,
            h,
            "native-uncertainty-checkpoint",
            MAX_ARTIFACT_BYTES,
        )?;
    }
    Ok(Loaded { hash, value, bytes })
}

fn drive(
    model: &Model,
    ledger: &Ledger,
    cap: Option<usize>,
    gate: &CancelGate,
    started: Instant,
    prior: Option<&Loaded>,
) -> Result<Outcome> {
    let plan = plan(model);
    let (mut execution, mut rows, used_before) = if let Some(old) = prior {
        if hash_field(&old.value, "model")? != model.identity() {
            return Err(fail("cli-uncertainty-resume", "retained model changed"));
        }
        let checkpoint = linked(
            ledger,
            &old.value,
            "checkpoint",
            "native-uncertainty-checkpoint",
        )?;
        let execution = UqExecution::restore(&plan, model.identity(), &checkpoint)
            .map_err(|e| fail("cli-uncertainty-resume", e.to_string()))?;
        let rows = read_rows(&linked(
            ledger,
            &old.value,
            "observations",
            "native-uncertainty-observations",
        )?)?;
        if rows.len() != execution.observations().len()
            || rows.len() != integer(&old.value, "samples_completed")?
        {
            return Err(fail(
                "cli-uncertainty-resume",
                "checkpoint and child-run counts differ",
            ));
        }
        // Replay only the cheap sampler, checking each actual child against its
        // exact addressed parameters; the physical prefix is never re-solved.
        let mut sampler = UqExecution::new(&plan).map_err(|e| fail("cli-uncertainty-plan", e))?;
        for row in &rows {
            model.verify_sample(ledger, row)?;
            let result = sampler.advance(
                1,
                || false,
                |parameters| {
                    if parameters
                        .iter()
                        .map(|v| v.to_bits())
                        .ne(row.parameters.iter().map(|v| v.to_bits()))
                    {
                        Err("retained parameters differ from the seeded sample ordinal")
                    } else {
                        Ok(row.value_k)
                    }
                },
            );
            if result.status == UqStatus::Refused {
                return Err(fail(
                    "cli-uncertainty-resume",
                    result.rejection_reason.unwrap_or_default(),
                ));
            }
        }
        if rows
            .iter()
            .zip(execution.observations())
            .any(|(r, v)| r.value_k.to_bits() != v.to_bits())
        {
            return Err(fail(
                "cli-uncertainty-resume",
                "checkpoint values differ from child QoIs",
            ));
        }
        let used = old
            .value
            .f64_field("consumed_wall_s")
            .filter(|v| v.is_finite() && *v >= 0.0)
            .ok_or_else(|| fail("cli-uncertainty-resume", "invalid retained wall charge"))?;
        (execution, rows, used)
    } else {
        (
            UqExecution::new(&plan).map_err(|e| fail("cli-uncertainty-plan", e))?,
            Vec::new(),
            0.0,
        )
    };
    let start_count = rows.len();
    let mut predecessor = prior.map(|old| old.hash);
    let mut interrupted = false;
    loop {
        let used = used_before + started.elapsed().as_secs_f64();
        let report = execution.report();
        let (status, termination) = if report.status == UqStatus::Refused {
            ("refused", "child-refused")
        } else if report.status == UqStatus::Complete {
            ("completed", "fixed-sample-count")
        } else if used >= model.bound.study().wall_seconds() {
            ("budget-exhausted", "wall-budget")
        } else if gate.is_requested() {
            ("cancelled", "cancelled")
        } else if interrupted {
            ("budget-exhausted", "child-budget")
        } else if cap.is_some_and(|cap| rows.len() - start_count >= cap) {
            ("budget-exhausted", "invocation-sample-budget")
        } else {
            ("running", "sampling")
        };
        let out = persist(
            model,
            ledger,
            &execution,
            &rows,
            status,
            termination,
            used,
            predecessor,
        )?;
        if status != "running" {
            return Ok(out);
        }
        predecessor = ContentHash::from_hex(out.pointer.trim_start_matches("study-"));
        let remaining =
            (model.bound.study().wall_seconds() - used_before - started.elapsed().as_secs_f64())
                .max(0.0);
        let mut accepted = None;
        execution.advance_interruptible(
            1,
            || gate.is_requested(),
            |parameters| {
                let sample = model.sample(ledger, gate, parameters, remaining)?;
                let value = sample.as_ref().map(|s| s.value_k);
                accepted = sample;
                Ok::<_, Failure>(value)
            },
        );
        if let Some(sample) = accepted {
            rows.push(sample);
        } else {
            interrupted = true;
        }
    }
}

pub(crate) fn looks_like(path: &Path) -> bool {
    let mut text = String::new();
    if std::fs::File::open(path)
        .and_then(|f| {
            f.take(fs_project::uncertainty::MAX_SOURCE_BYTES as u64 + 1)
                .read_to_string(&mut text)
        })
        .is_err()
    {
        return false;
    }
    // This is only routing. The full typed parser owns all admission.
    text.lines()
        .map(|line| line.split(';').next().unwrap_or(""))
        .collect::<String>()
        .trim_start()
        .starts_with("(fsim-uncertainty-study")
}
pub(crate) fn owns_run(pointer: &str, path: &Path) -> bool {
    if legacy::owns(pointer) {
        return true;
    }
    if !path.is_file() {
        return false;
    }
    let Some(hash) = pointer
        .strip_prefix("study-")
        .and_then(ContentHash::from_hex)
    else {
        return false;
    };
    let Some(path) = path.to_str() else {
        return false;
    };
    let Ok(ledger) = Ledger::open(path) else {
        return false;
    };
    artifact(&ledger, hash, RECEIPT_KIND, MAX_ARTIFACT_BYTES)
        .and_then(|bytes| parse(&bytes))
        .is_ok_and(|value| value.str_field("driver") == Some(DRIVER))
}
pub(crate) fn study_path(
    path: &Path,
    ledger_path: &Path,
    override_text: Option<&str>,
    mode: OutputMode,
) -> CommandOutput {
    let result = (|| {
        let cap = budget(override_text)?;
        let model = Model::load(path)?;
        let ledger = Ledger::open(
            ledger_path
                .to_str()
                .ok_or_else(|| fail("cli-uncertainty-ledger", "ledger path is not UTF-8"))?,
        )?;
        drive(
            &model,
            &ledger,
            cap,
            &CancelGate::new_clock_free(),
            Instant::now(),
            None,
        )
    })();
    match result {
        Ok(out) => render(out, mode),
        Err(e) => output_error("study", mode, e),
    }
}
pub(crate) fn resume_path(
    pointer: &str,
    path: &Path,
    override_text: Option<&str>,
    mode: OutputMode,
) -> CommandOutput {
    if legacy::owns(pointer) {
        return legacy::resume(mode);
    }
    let result = (|| {
        let cap = budget(override_text)?;
        if !path.is_file() {
            return Err(fail(
                "cli-uncertainty-ledger",
                "resume requires an existing ledger",
            ));
        }
        let ledger = Ledger::open(
            path.to_str()
                .ok_or_else(|| fail("cli-uncertainty-ledger", "ledger path is not UTF-8"))?,
        )?;
        let started = Instant::now();
        let old = load(&ledger, pointer)?;
        let status = old
            .value
            .str_field("status")
            .ok_or_else(|| fail("cli-uncertainty-run", "missing status"))?;
        if matches!(status, "completed" | "refused") {
            return Ok(Outcome {
                pointer: pointer.into(),
                receipt: String::from_utf8(old.bytes)
                    .map_err(|e| fail("cli-uncertainty-run", e.to_string()))?,
                status: status.into(),
            });
        }
        let model = Model::restore(&ledger, hash_field(&old.value, "model")?)?;
        drive(
            &model,
            &ledger,
            cap,
            &CancelGate::new_clock_free(),
            started,
            Some(&old),
        )
    })();
    match result {
        Ok(out) => render(out, mode),
        Err(e) => output_error("study", mode, e),
    }
}
pub(crate) fn export(
    command: &'static str,
    pointer: &str,
    path: Option<&Path>,
    mode: OutputMode,
) -> CommandOutput {
    if legacy::owns(pointer) {
        return legacy::export(command, pointer, path, mode);
    }
    let result = (|| {
        let path = path.filter(|p| p.is_file()).ok_or_else(|| {
            fail(
                "cli-uncertainty-ledger",
                "export requires an existing ledger",
            )
        })?;
        let ledger = Ledger::open(
            path.to_str()
                .ok_or_else(|| fail("cli-uncertainty-ledger", "ledger path is not UTF-8"))?,
        )?;
        let old = load(&ledger, pointer)?;
        let fields: &[(&str, &str, &str)] = if command == "package" {
            &[("package", "study-package", "fspkg")]
        } else {
            &[
                ("report_json", "study-report-json", "json"),
                ("report_html", "study-report-html", "html"),
            ]
        };
        let mut paths = String::new();
        let mut text_paths = String::new();
        for &(key, kind, extension) in fields {
            let bytes = linked(&ledger, &old.value, key, kind)?;
            if key == "package" {
                let package = EvidencePackage::from_json(
                    std::str::from_utf8(&bytes)
                        .map_err(|e| fail("cli-uncertainty-package", e.to_string()))?,
                )
                .map_err(|e| fail("cli-uncertainty-package", e.to_string()))?;
                if !fs_checker::check(&package).passed() {
                    return Err(fail(
                        "cli-uncertainty-package",
                        "retained package failed its structural check",
                    ));
                }
            }
            let dest = path
                .parent()
                .unwrap_or_else(|| Path::new("."))
                .join(format!("{pointer}.{extension}"));
            crate::report::write_retained(&dest, &bytes)
                .map_err(|e| fail("cli-uncertainty-export", e))?;
            let _ = write!(paths, ",{key:?}:{}", quoted(&dest.to_string_lossy()));
            let _ = writeln!(
                text_paths,
                "{key}={}",
                crate::escape_text(&dest.to_string_lossy())
            );
        }
        let status = old.value.str_field("status").unwrap_or("unknown");
        Ok(CommandOutput {
            exit_code: exit::SUCCESS,
            stderr: String::new(),
            stdout: match mode {
                OutputMode::Json => format!(
                    "{{\"command\":{command:?},\"status\":\"ok\",\"run\":{},\"study_status\":{},\"authority\":\"projection-of-retained-estimates\",\"verification\":\"sealed-evidence\"{paths}}}\n",
                    quoted(pointer),
                    quoted(status)
                ),
                OutputMode::Text => format!(
                    "command={command}\nstatus=ok\nrun={pointer}\nstudy_status={status}\n{text_paths}"
                ),
            },
        })
    })();
    match result {
        Ok(out) => out,
        Err(e) => output_error(command, mode, e),
    }
}
