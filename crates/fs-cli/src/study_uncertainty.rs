//! Native probability studies use the ordinary import/solve/report producers.
//! A finite, sealed Kelvin QoI is the only successful sample. Probability laws
//! are explicit assumptions; empirical statistics never replace a solve verdict.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::Instant;

use fs_blake3::{ContentHash, hash_bytes, hash_domain};
use fs_exec::CancelGate;
use fs_ledger::{EdgeRole, FiveExplicits, Ledger, OpOutcome};
use fs_project::uncertainty::{BoundStudy, UncertaintyStudy};
use fs_uq::{CorrelationModel, ParameterUncertainty, PropagationMethod, UqExecution, UqPlan, UqStatus};

use crate::json_read::JsonValue;
use crate::{CardPackKind, CardPackSet, CommandOutput, Diagnostic, GeometryImportLimits,
    OutputMode, RawCardPack, RawGeometryLibrary, SolveRunStatus, exit, push_json_string, refusal};

const DRIVER: &str = "native-cooling-uncertainty-v1";
const KIND: &str = "native-uncertainty-receipt";
const PREFIX: &str = "study-uq-";
const MAX_INPUT_BYTES: u64 = 64 * 1024 * 1024;
const MAX_RECEIPT_BYTES: u64 = 16 * 1024 * 1024;
const NO_CLAIM: &str = "Estimated fixed-count Monte Carlo of explicit independent uniform inputs. Numerical threshold frequency is not an engineering compliance verdict. No certified output enclosure, optional-stopping confidence bound, model validation or guaranteed memory bound. Total wall time is checked at import/solve phase boundaries; individual native kernels retain their own budgets and cancellation.";
type Result<T> = std::result::Result<T, String>;

fn quote(text: &str) -> String {
    let mut result = String::new();
    push_json_string(&mut result, text);
    result
}
fn diagnostic(command: &'static str, mode: OutputMode, message: impl Into<String>) -> CommandOutput {
    refusal(mode, exit::REFUSED, &Diagnostic::new(command, "cli-native-uncertainty",
        message, "inspect the explicit study inputs or retained native solve refusal; failed samples cannot be skipped"), None)
}
fn bounded(path: &Path, cap: u64) -> Result<Vec<u8>> {
    let file = std::fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let metadata = file.metadata().map_err(|e| e.to_string())?;
    if !metadata.is_file() || metadata.len() > cap {
        return Err(format!("{} must be a regular file within {cap} bytes", path.display()));
    }
    let mut bytes = Vec::new();
    file.take(cap.saturating_add(1)).read_to_end(&mut bytes).map_err(|e| e.to_string())?;
    if bytes.len() as u64 > cap { return Err("input grew beyond its allowance".into()); }
    Ok(bytes)
}
fn text(bytes: &[u8]) -> Result<&str> {
    std::str::from_utf8(bytes).map_err(|e| e.to_string())
}
fn allowance(value: Option<&str>, default: usize) -> Result<usize> {
    match value {
        None => Ok(default),
        Some(value) => value.parse::<usize>().ok().filter(|n| (1..=256).contains(n))
            .ok_or_else(|| "--budget must be an integer in 1..=256".into()),
    }
}
pub(crate) fn looks_like(path: &Path) -> bool {
    bounded(path, fs_project::uncertainty::MAX_SOURCE_BYTES as u64).ok()
        .and_then(|bytes| String::from_utf8(bytes).ok())
        .is_some_and(|source| source.contains("(fsim-uncertainty-study"))
}
pub(crate) fn owns_run(pointer: &str) -> bool { pointer.starts_with(PREFIX) }

struct Blob {
    kind: &'static str,
    bytes: Vec<u8>,
}
struct Inputs {
    bound: BoundStudy,
    raw: RawGeometryLibrary,
    cards: CardPackSet,
    limits: GeometryImportLimits,
    blobs: Vec<Blob>,
    manifest: String,
    identity: ContentHash,
}
fn plan(bound: &BoundStudy) -> UqPlan {
    UqPlan {
        target_qoi: bound.study().qoi().into(),
        compliance_threshold: Some(bound.threshold_k()),
        parameters: bound.study().parameters().iter().map(|p|
            ParameterUncertainty::uniform(&p.name, p.low, p.high, p.target.unit())).collect(),
        correlation: CorrelationModel::Independent,
        method: PropagationMethod::MonteCarlo,
        budget_max_samples: bound.study().samples(),
        seed: bound.study().seed(),
    }
}
fn manifest(blobs: &[Blob]) -> String {
    let entries = blobs.iter().map(|b| format!("{{\"kind\":{},\"hash\":{}}}",
        quote(b.kind), quote(&hash_bytes(&b.bytes).to_hex()))).collect::<Vec<_>>().join(",");
    format!("{{\"driver\":{DRIVER:?},\"solver_version\":{},\"constellation\":{},\"inputs\":[{entries}]}}",
        crate::SOLVE_DRIVER_VERSION, quote(&hash_bytes(include_bytes!("../../../constellation.lock")).to_hex()))
}
fn assemble(blobs: Vec<Blob>) -> Result<Inputs> {
    if blobs.len() < 2 || blobs[0].kind != "native-uq-source" || blobs[1].kind != "native-uq-project" {
        return Err("missing native source and base project".into());
    }
    let study = UncertaintyStudy::parse(text(&blobs[0].bytes)?).map_err(|e| e.detail)?;
    let decoded = fs_project::parse_sexpr_migrating(text(&blobs[1].bytes)?).map_err(|e| e.detail)?.decoded;
    let bound = study.bind(&decoded.spec).map_err(|e| e.detail)?;
    let declaration = bound.study();
    if blobs.len() != 2 + declaration.geometry().len() + declaration.materials().len() + declaration.interfaces().len() {
        return Err("retained source count disagrees with the declaration".into());
    }
    let memory = bound.base().budgets.as_ref().ok_or("missing native memory budget")?.memory_bytes;
    let cap = MAX_INPUT_BYTES.min(memory / 4);
    let total = blobs.iter().try_fold(0_u64, |sum, blob| sum.checked_add(blob.bytes.len() as u64))
        .ok_or("source length overflow")?;
    if total > cap { return Err(format!("native study sources exceed the {cap}-byte admission allowance")); }
    let mut limits = GeometryImportLimits::DEFAULT;
    limits.max_source_bytes = limits.max_source_bytes.min(cap as usize);
    limits.max_total_source_bytes = limits.max_total_source_bytes.min(cap as usize);
    let mut raw = RawGeometryLibrary::new();
    let mut cursor = 2;
    for source in declaration.geometry() {
        let blob = &blobs[cursor];
        if blob.kind != "native-uq-mesh" { return Err("wrong retained geometry kind".into()); }
        let artifact = bound.base().geometry.as_ref().ok_or("missing native geometry")?.iter()
            .find(|row| row.role == source.role).ok_or("missing native geometry role")?;
        raw.insert_mesh(artifact, source.path.clone(), blob.bytes.clone(), source.unit.clone(),
            source.max_hole_edges, Vec::new());
        cursor += 1;
    }
    let mut packs = Vec::new();
    for (paths, kind, expected) in [
        (declaration.materials(), CardPackKind::Material, "native-uq-material"),
        (declaration.interfaces(), CardPackKind::Interface, "native-uq-interface"),
    ] {
        for path in paths {
            let blob = &blobs[cursor];
            if blob.kind != expected { return Err("wrong retained card kind".into()); }
            packs.push(RawCardPack { kind, source: path.clone(), bytes: blob.bytes.clone(), expect: None });
            cursor += 1;
        }
    }
    let cards = CardPackSet::admit(packs).map_err(|e| format!("{}: {}", e.code, e.what))?;
    let manifest = manifest(&blobs);
    let identity = hash_domain("org.frankensim.cli.native-uq-model.v1", manifest.as_bytes());
    Ok(Inputs { bound, raw, cards, limits, blobs, manifest, identity })
}
fn admit(path: &Path) -> Result<Inputs> {
    let source = bounded(path, fs_project::uncertainty::MAX_SOURCE_BYTES as u64)?;
    let study = UncertaintyStudy::parse(text(&source)?).map_err(|e| e.detail)?;
    let root = path.parent().unwrap_or_else(|| Path::new("."));
    let decoded = crate::read_project_for_solve(&root.join(study.project_path()), OutputMode::Json)
        .map_err(|e| e.stderr)?;
    let memory = decoded.spec.budgets.as_ref().ok_or("missing memory budget")?.memory_bytes;
    let cap = MAX_INPUT_BYTES.min(memory / 4);
    let mut blobs = vec![
        Blob { kind: "native-uq-source", bytes: study.canonical().as_bytes().to_vec() },
        Blob { kind: "native-uq-project", bytes: fs_project::print_sexpr(&decoded.spec).into_bytes() },
    ];
    // Resolve every path relative to the study, and read it ONCE. Later samples
    // and resumes consume the retained bytes, never a mutable path alias.
    let mut remaining = cap.checked_sub(blobs.iter().map(|b| b.bytes.len() as u64).sum())
        .ok_or("source metadata exceeds memory admission")?;
    for (path, kind) in study.geometry().iter().map(|g| (&g.path, "native-uq-mesh"))
        .chain(study.materials().iter().map(|p| (p, "native-uq-material")))
        .chain(study.interfaces().iter().map(|p| (p, "native-uq-interface"))) {
        let bytes = bounded(&root.join(path), remaining)?;
        remaining -= bytes.len() as u64;
        blobs.push(Blob { kind, bytes });
    }
    assemble(blobs)
}

fn transaction<T>(ledger: &Ledger, work: impl FnOnce() -> Result<T>) -> Result<T> {
    if ledger.in_transaction() { return Err("native UQ requires its own transaction".into()); }
    ledger.begin().map_err(|e| e.to_string())?;
    let result = work().and_then(|value| ledger.commit().map(|()| value).map_err(|e| e.to_string()));
    match result {
        Ok(value) => Ok(value),
        Err(error) => match ledger.rollback() {
            Ok(()) => Err(error),
            Err(rollback) => Err(format!("{error}; rollback also failed: {rollback}")),
        },
    }
}
fn retain_inputs(ledger: &Ledger, input: &Inputs) -> Result<ContentHash> {
    transaction(ledger, || {
        for blob in &input.blobs { ledger.put_artifact(blob.kind, &blob.bytes, None).map_err(|e| e.to_string())?; }
        Ok(ledger.put_artifact("native-uq-manifest", input.manifest.as_bytes(), None)
            .map_err(|e| e.to_string())?.hash)
    })
}
fn qoi_value(value: &JsonValue) -> Result<f64> {
    let rows = value.get("qoi").and_then(JsonValue::as_array).ok_or("missing QoI rows")?;
    let mut matching = rows.iter().filter(|r| r.str_field("name") == Some("temperature-max"));
    let row = matching.next().ok_or("missing temperature-max")?;
    if matching.next().is_some() || row.str_field("unit") != Some("K") {
        return Err("ambiguous temperature-max or non-Kelvin unit".into());
    }
    row.f64_field("value").filter(|v| v.is_finite() && *v >= 0.0)
        .ok_or_else(|| "temperature-max must be a finite nonnegative Kelvin observation".into())
}
fn evaluate(input: &Inputs, ledger: &Ledger, parameters: &[f64], gate: &CancelGate,
    started: Instant, wall: f64, ordinal: usize) -> Result<Option<(f64, String)>> {
    if gate.is_requested() || started.elapsed().as_secs_f64() >= wall { return Ok(None); }
    let project = input.bound.sample_project(parameters).map_err(|e| e.detail)?;
    let decoded = fs_project::parse_sexpr_migrating(&fs_project::print_sexpr(&project))
        .map_err(|e| e.detail)?.decoded;
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let imported = pool.scope(|arena| {
        let cx = fs_exec::Cx::new(gate, arena, fs_exec::StreamKey {
            seed: project.seeds.as_ref().map_or(0, |s| s.root),
            kernel_id: 0x66_73_75_71, tile: 0, iteration: 0,
        }, fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
        crate::import_project_geometry(&project, &input.raw, ledger, input.limits, &cx)
    });
    if gate.is_requested() { return Ok(None); }
    imported.map_err(|e| format!("native import {}: {}", e.code, e.what))?;
    if started.elapsed().as_secs_f64() >= wall { return Ok(None); }
    let sample_start = Instant::now();
    let mut clock = || sample_start.elapsed().as_secs_f64();
    let mut progress = Vec::new();
    let solved = crate::run_solve(ledger, gate, &mut clock, &decoded, &input.cards, &mut progress);
    if gate.is_requested() { return Ok(None); }
    let solved = solved.map_err(|e| format!("native solve {}: {}", e.code, e.what))?;
    if !matches!(solved.status, SolveRunStatus::Completed) { return Ok(None); }
    // Use the ordinary export attestation, not an unsealed intermediate field.
    let loaded = crate::report::load_export_from("study", &solved.run, ledger,
        PathBuf::from("."), OutputMode::Json).map_err(|e| e.stderr)?;
    let mut stages = loaded.export.stages.iter().filter(|(stage, _, _)| *stage == "qoi");
    let (_, _, receipt) = stages.next().ok_or("completed solve has no QoI receipt")?;
    if stages.next().is_some() { return Err("completed solve has duplicate QoI receipts".into()); }
    let hash = ContentHash::from_hex(receipt).ok_or("invalid QoI receipt hash")?;
    let bytes = ledger.get_artifact_bounded(&hash, MAX_RECEIPT_BYTES).map_err(|e| e.to_string())?
        .ok_or("missing retained QoI receipt")?;
    let qoi = JsonValue::parse(text(&bytes)?).map_err(|e| e.to_string())?;
    if qoi.str_field("run") != Some(solved.run.as_str())
        || qoi.str_field("project_hash") != Some(loaded.export.project_hash.as_str()) {
        return Err("QoI receipt is not bound to this sample solve".into());
    }
    let value = qoi_value(&qoi)?;
    let parameters = parameters.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",");
    let row = format!("{{\"ordinal\":{ordinal},\"parameters\":[{parameters}],\"value_kelvin\":{value},\"run\":{},\"project_hash\":{},\"qoi_receipt\":{}}}",
        quote(&solved.run), quote(&loaded.export.project_hash), quote(receipt));
    // Once a solve completed, retain its paid observation even if the wall
    // allowance elapsed during that indivisible phase. Stop before the next one.
    Ok(Some((value, row)))
}
fn number(value: Option<f64>) -> String {
    value.filter(|v| v.is_finite()).map_or_else(|| "null".into(), |v| v.to_string())
}
fn statistics(execution: &UqExecution) -> String {
    let report = execution.report();
    let percentiles = report.percentiles.map_or_else(|| "null".into(), |p|
        format!("[{},{},{}]", p[0], p[1], p[2]));
    format!("{{\"mean_kelvin\":{},\"std_dev_kelvin\":{},\"standard_error_kelvin\":{},\"percentiles_05_50_95_kelvin\":{percentiles},\"empirical_threshold_frequency\":{},\"failure\":{}}}",
        number(report.mean), number(report.std_dev),
        number(report.std_dev.map(|_| report.sampling_error)),
        number(report.probability_of_compliance), report.rejection_reason.as_deref().map_or_else(|| "null".into(), quote))
}
struct Outcome { pointer: String, receipt: String, status: String }
fn retain(ledger: &Ledger, input: &Inputs, manifest: ContentHash, execution: &UqExecution,
    rows: &[String], wall: f64, status: &str, predecessor: Option<ContentHash>) -> Result<Outcome> {
    let checkpoint = if execution.report().status == UqStatus::Refused { None }
        else { Some(execution.checkpoint(input.identity).map_err(|e| e.to_string())?) };
    transaction(ledger, || {
        let seed = input.bound.study().seed().to_le_bytes();
        let versions = format!("{{\"driver\":{DRIVER:?},\"solver\":{}}}", crate::SOLVE_DRIVER_VERSION);
        let budget = format!("{{\"samples\":{},\"wall_s\":{},\"consumed_wall_s\":{wall}}}",
            input.bound.study().samples(), input.bound.study().wall_seconds());
        let ir = format!("{{\"driver\":{DRIVER:?},\"model\":{},\"accepted\":{}}}",
            quote(&input.identity.to_hex()), execution.observations().len());
        let op = ledger.begin_op(Some(input.identity.as_bytes()), &ir, &FiveExplicits {
            seed: &seed, versions: &versions, budget: &budget,
            capability: "{\"ops\":[\"uncertainty.native-cooling\"],\"authority\":\"estimated\"}",
        }, 0).map_err(|e| e.to_string())?;
        ledger.link(op, &manifest, EdgeRole::In).map_err(|e| e.to_string())?;
        if let Some(previous) = predecessor { ledger.link(op, &previous, EdgeRole::In).map_err(|e| e.to_string())?; }
        let checkpoint_hash = if let Some(bytes) = checkpoint {
            let artifact = ledger.put_artifact("native-uq-checkpoint", &bytes, None).map_err(|e| e.to_string())?;
            ledger.link(op, &artifact.hash, EdgeRole::Out).map_err(|e| e.to_string())?;
            quote(&artifact.hash.to_hex())
        } else { "null".into() };
        let receipt = format!("{{\"schema\":{},\"driver\":{DRIVER:?},\"model_identity\":{},\"manifest\":{},\"status\":{},\"samples_accepted\":{},\"samples_attempted\":{},\"target_samples\":{},\"consumed_wall_s\":{wall},\"checkpoint\":{checkpoint_hash},\"predecessor\":{},\"observations\":[{}],\"statistics\":{},\"authority\":\"Estimated\",\"no_claim\":{}}}",
            quote(super::STUDY_RUN_RECEIPT_SCHEMA), quote(&input.identity.to_hex()), quote(&manifest.to_hex()), quote(status),
            execution.observations().len(), execution.evaluations_attempted(), input.bound.study().samples(),
            predecessor.map_or_else(|| "null".into(), |h| quote(&h.to_hex())), rows.join(","), statistics(execution), quote(NO_CLAIM));
        let artifact = ledger.put_artifact(KIND, receipt.as_bytes(), None).map_err(|e| e.to_string())?;
        ledger.link(op, &artifact.hash, EdgeRole::Out).map_err(|e| e.to_string())?;
        if ledger.artifact_output_seal(&artifact.hash).map_err(|e| e.to_string())?.is_none() {
            ledger.seal_artifact_output(&artifact.hash, op).map_err(|e| e.to_string())?;
        }
        ledger.finish_op(op, OpOutcome::Ok, None, 1).map_err(|e| e.to_string())?;
        Ok(Outcome { pointer: format!("{PREFIX}{}", artifact.hash.to_hex()), receipt, status: status.into() })
    })
}
fn render(command: &'static str, outcome: Outcome, mode: OutputMode) -> CommandOutput {
    let exit_code = match outcome.status.as_str() {
        "complete" => exit::SUCCESS, "cancelled" => exit::CANCELLED,
        "refused" => exit::REFUSED, _ => exit::BUDGET,
    };
    let stdout = match mode {
        OutputMode::Json => format!("{{\"command\":{},\"status\":{},\"run\":{},\"run_id\":{},\"receipt\":{}}}\n",
            quote(command), quote(&outcome.status), quote(&outcome.pointer), quote(&outcome.pointer), outcome.receipt),
        OutputMode::Text => format!("command={command}\nstatus={}\nrun={}\nauthority=Estimated\nreceipt={}\n",
            outcome.status, outcome.pointer, outcome.receipt),
    };
    CommandOutput { exit_code, stdout, stderr: String::new() }
}
fn drive(input: &Inputs, ledger: &Ledger, cap: usize, mut execution: UqExecution,
    mut rows: Vec<String>, used_wall: f64, mut predecessor: Option<ContentHash>) -> Result<Outcome> {
    let manifest = retain_inputs(ledger, input)?;
    let started = Instant::now();
    let remaining = (input.bound.study().wall_seconds() - used_wall).max(0.0);
    let gate = CancelGate::new();
    let mut last = None;
    for _ in 0..cap {
        if execution.report().status == UqStatus::Complete { break; }
        if started.elapsed().as_secs_f64() >= remaining { break; }
        let ordinal = execution.observations().len();
        let mut row = None;
        let report = execution.advance_interruptible(1, || gate.is_requested(), |parameters| {
            evaluate(input, ledger, parameters, &gate, started, remaining, ordinal).map(|result|
                result.map(|(value, receipt)| { row = Some(receipt); value }))
        });
        if let Some(row) = row { rows.push(row); }
        let status = if report.status == UqStatus::Cancelled && !gate.is_requested() {
            "budget-truncated"
        } else { report.status.label() };
        let outcome = retain(ledger, input, manifest, &execution, &rows,
            used_wall + started.elapsed().as_secs_f64(), status, predecessor)?;
        predecessor = outcome.pointer.strip_prefix(PREFIX).and_then(ContentHash::from_hex);
        last = Some(outcome);
        if report.status != UqStatus::BudgetTruncated { break; }
    }
    match last {
        Some(outcome) => Ok(outcome),
        None => {
            let report = execution.report();
            let status = if report.status == UqStatus::Complete { "complete" }
                else if gate.is_requested() { "cancelled" } else { "budget-truncated" };
            retain(ledger, input, manifest, &execution, &rows,
                used_wall + started.elapsed().as_secs_f64(), status, predecessor)
        }
    }
}
pub(crate) fn study_path(path: &Path, ledger_path: &Path, override_text: Option<&str>, mode: OutputMode) -> CommandOutput {
    let result: Result<Outcome> = (|| {
        let input = admit(path)?;
        let cap = allowance(override_text, input.bound.study().samples())?;
        let execution = UqExecution::new(&plan(&input.bound)).map_err(str::to_string)?;
        let ledger = Ledger::open(ledger_path.to_str().ok_or("ledger path is not UTF-8")?).map_err(|e| e.to_string())?;
        drive(&input, &ledger, cap, execution, Vec::new(), 0.0, None)
    })();
    match result { Ok(outcome) => render("study", outcome, mode), Err(e) => diagnostic("study", mode, e) }
}
pub(crate) fn resume_path(_pointer: &str, _path: &Path, _budget: Option<&str>, mode: OutputMode) -> CommandOutput {
    diagnostic("study", mode, "native UQ checkpoints are retained; this driver revision does not yet admit resumption")
}
pub(crate) fn export(command: &'static str, pointer: &str, path: Option<&Path>, mode: OutputMode) -> CommandOutput {
    let result: Result<Outcome> = (|| {
        let (ledger, _) = crate::report::open_export_ledger(command, pointer, path, mode).map_err(|e| e.stderr)?;
        let hash = pointer.strip_prefix(PREFIX).and_then(ContentHash::from_hex).ok_or("invalid native study pointer")?;
        let info = ledger.artifact_info(&hash).map_err(|e| e.to_string())?.ok_or("missing study receipt")?;
        if info.kind != KIND || ledger.artifact_output_seal(&hash).map_err(|e| e.to_string())?.is_none() {
            return Err("study receipt is not a sealed native UQ result".into());
        }
        let bytes = ledger.get_artifact_bounded(&hash, MAX_RECEIPT_BYTES).map_err(|e| e.to_string())?.ok_or("missing study bytes")?;
        let receipt = text(&bytes)?.to_string();
        let parsed = JsonValue::parse(&receipt).map_err(|e| e.to_string())?;
        if parsed.str_field("driver") != Some(DRIVER) { return Err("wrong study driver".into()); }
        let status = parsed.str_field("status").ok_or("missing study status")?.to_string();
        Ok(Outcome { pointer: pointer.into(), receipt, status })
    })();
    match result { Ok(outcome) => render(command, outcome, mode), Err(e) => diagnostic(command, mode, e) }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn quantity_requires_one_finite_kelvin_observation() {
        for source in [r#"{"qoi":[]}"#, r#"{"qoi":[{"name":"temperature-max","unit":"C","value":20}]}"#,
            r#"{"qoi":[{"name":"temperature-max","unit":"K","value":-1}]}"#,
            r#"{"qoi":[{"name":"temperature-max","unit":"K","value":300},{"name":"temperature-max","unit":"K","value":300}]}"#] {
            assert!(qoi_value(&JsonValue::parse(source).unwrap()).is_err());
        }
        let parsed = JsonValue::parse(r#"{"qoi":[{"name":"temperature-max","unit":"K","value":310.25}]}"#).unwrap();
        assert_eq!(qoi_value(&parsed).unwrap(), 310.25);
    }
    #[test]
    fn sample_override_never_changes_the_lifetime_plan() {
        assert_eq!(allowance(None, 10).unwrap(), 10);
        assert_eq!(allowance(Some("1"), 10).unwrap(), 1);
        for invalid in ["0", "257", "-1", "NaN", "1.5"] { assert!(allowance(Some(invalid), 10).is_err()); }
    }
    #[test]
    fn missing_dispersion_is_null_not_zero_certainty() {
        let mut p = UqPlan::new("temperature-max", PropagationMethod::MonteCarlo, 2)
            .with_parameter(ParameterUncertainty::uniform("power", 1.0, 2.0, "W"));
        p.correlation = CorrelationModel::Independent;
        let mut execution = UqExecution::new(&p).unwrap();
        execution.advance(1, || false, |_| Ok::<_, String>(300.0));
        let value = JsonValue::parse(&statistics(&execution)).unwrap();
        assert_eq!(value.f64_field("mean_kelvin"), Some(300.0));
        assert!(value.f64_field("standard_error_kelvin").is_none());
        execution.advance(1, || false, |_| Err::<f64, _>("native solve refused"));
        let value = JsonValue::parse(&statistics(&execution)).unwrap();
        assert!(value.f64_field("mean_kelvin").is_none());
        assert!(value.str_field("failure").unwrap().contains("native solve refused"));
    }
}
