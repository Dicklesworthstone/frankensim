//! Executable normalized thermal studies. The numerical producer owns accepted
//! iterates; this membrane admits intent, meters work and retains real results.
//! Export/checker integrity does not turn an estimated PDE result into a bound.

use std::fmt::Write as _;
use std::io::Read;
use std::path::Path;
use std::time::Instant;

use fs_blake3::{ContentHash, hash_bytes, hash_domain};
use fs_exec::CancelGate;
use fs_ledger::{EdgeRole, FiveExplicits, Ledger, LedgerError, OpOutcome};
use fs_marquee::study::{PlateWithHoles, StudyConfig, StudyRunner};
use fs_package::{Claim, EvidencePackage, Provenance};
use fs_project::study::{StudySpec, parse_study_strict, print_study_sexpr};
use fs_session::{CapabilityToken, Charge, Enforcement, Governor, SessionId};

use crate::json_read::JsonValue;
use crate::{
    CommandOutput, Diagnostic, MAX_PROJECT_BYTES, OutputMode, exit, push_json_string, refusal,
};

pub const STUDY_RUN_RECEIPT_SCHEMA: &str = "frankensim.cli.study-run-receipt.v1";
const RECEIPT_KIND: &str = "study-run-receipt";
const MAX_ARTIFACT_BYTES: u64 = 16 * 1024 * 1024;
const DRIVER: &str = "normalized-thermal-study-v1";
const NO_CLAIM: &str = "Normalized scalar Poisson compliance; fixed centers and circular hole radii. Estimated DWR and algebraic terms, no guaranteed error bound, elasticity, free-boundary topology, KKT or optimality claim. Memory is an admission estimate, not measured peak RSS. Cancellation is checked between bounded iterations, not inside CutFEM solves. The checker proves package structure only.";

type Result<T> = std::result::Result<T, Failure>;
#[derive(Debug)]
struct Failure {
    code: &'static str,
    message: String,
    exit: u8,
}
impl From<LedgerError> for Failure {
    fn from(e: LedgerError) -> Self {
        fail("cli-study-ledger", e.to_string())
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
fn failure_output(command: &'static str, mode: OutputMode, e: Failure) -> CommandOutput {
    refusal(
        mode,
        e.exit,
        &Diagnostic::new(
            command,
            e.code,
            e.message,
            "inspect the named declaration or ledger failure; use examples/marquee/thermal-2d.fsim for the supported study",
        ),
        None,
    )
}

struct Admitted {
    spec: StudySpec,
    source: String,
    id: ContentHash,
    config: StudyConfig,
    wall_s: f64,
    memory_bytes: u64,
}
fn admit(source: &str, json: bool) -> Result<Admitted> {
    let spec = parse_study_strict(source, json).map_err(|e| fail(e.code, e.detail))?;
    if let Some(v) = spec.validate().first() {
        return Err(fail(v.code, &v.what));
    }
    if spec.versions.as_ref().expect("validated versions").schema != fs_project::STUDY_FSIM_VERSION
    {
        return Err(fail(
            "cli-study-version",
            "versions.schema must match the admitted study syntax version",
        ));
    }
    let physics = spec.physics.as_ref().expect("validated physics");
    let domain = spec.domain.as_ref().expect("validated domain");
    let scenario = spec.scenario.as_ref().expect("validated scenario");
    let objective = spec.objective.as_ref().expect("validated objective");
    let opt = spec.optimizer.as_ref().expect("validated optimizer");
    let budgets = spec.budgets.as_ref().expect("validated budgets");
    if physics.physics_type != "thermal-poisson-2d-normalized" {
        return Err(fail(
            "cli-study-physics-unavailable",
            "only thermal-poisson-2d-normalized is implemented; elasticity/topology remains q61wp.16",
        ));
    }
    if domain.domain_type != "sdf-plate-with-holes"
        || domain.bounds != ([0.0, 0.0], [1.0, 1.0])
        || scenario.fixed_boundary != "all-boundaries-zero"
        || scenario.load_region != "domain-unit-source"
        || objective.objective_type != "compliance"
        || objective.sense != "minimize"
        || objective.unit != "1"
        || spec.units.as_ref().expect("units").storage != "normalized"
        || opt.optimizer_type != "projected-gradient"
    {
        return Err(fail(
            "cli-study-model",
            "declare the normalized unit plate, unit domain source, zero temperature on all boundaries, dimensionless compliance and projected-gradient optimizer",
        ));
    }
    let caps = spec.capabilities.as_ref().expect("capabilities");
    if caps.len() != 3
        || ![
            "optimization.thermal-radii",
            "geometry.sdf",
            "physics.cutfem",
        ]
        .iter()
        .all(|cap| caps.iter().any(|c| c == cap))
    {
        return Err(fail(
            "cli-study-capability",
            "required capabilities are optimization.thermal-radii, geometry.sdf and physics.cutfem",
        ));
    }
    let metadata = spec.metadata.as_ref().expect("metadata");
    if metadata.decision_gate != fs_project::DecisionGate::ScopingEstimate
        || metadata.consequence != fs_project::ConsequenceClass::Advisory
    {
        return Err(fail(
            "cli-study-decision",
            "this estimated study supports advisory scoping only",
        ));
    }
    let wall = budgets
        .wall_time
        .as_ref()
        .ok_or_else(|| fail("cli-study-budget", "wall time must be explicit"))?;
    let memory_bytes = budgets
        .memory_bytes
        .ok_or_else(|| fail("cli-study-budget", "memory must be explicit"))?;
    let max_iterations = budgets
        .max_iterations
        .ok_or_else(|| fail("cli-study-budget", "max-iterations must be explicit"))?;
    if !wall.value.is_finite()
        || wall.value <= 0.0
        || wall.dims != fs_qty::Dims([0, 0, 1, 0, 0, 0])
        || max_iterations == 0
        || max_iterations > 256
        || opt.steps > 256
        || domain.initial_holes.len() > 32
        || physics.mesh_level > 5
        || memory_bytes < 128 * 1024 * 1024
    {
        return Err(fail(
            "cli-study-budget",
            "supported envelope: positive wall seconds, 128 MiB memory admission, at most 256 iterations, 32 holes and mesh level 5",
        ));
    }
    let config = StudyConfig {
        level: physics.mesh_level,
        steps: opt.steps,
        step_size: opt.step_size,
        area_target: spec
            .constraints
            .as_ref()
            .expect("constraints")
            .volume_fraction,
        r_min: opt.r_min,
        r_max: opt.r_max,
    };
    let canonical = print_study_sexpr(&spec);
    let identity = format!(
        "{DRIVER}\n{}\n{canonical}",
        hash_bytes(include_bytes!("../../../constellation.lock")).to_hex()
    );
    let id = hash_domain(
        "org.frankensim.cli.normalized-thermal-study.v1",
        identity.as_bytes(),
    );
    let wall_s = wall.value;
    // Admission includes the actual projection, so no invalid layout creates a ledger.
    StudyRunner::new(design(&spec), config.clone())
        .map_err(|e| fail("cli-study-geometry", e.to_string()))?;
    Ok(Admitted {
        spec,
        source: canonical,
        id,
        config,
        wall_s,
        memory_bytes,
    })
}
fn design(spec: &StudySpec) -> PlateWithHoles {
    let holes = &spec.domain.as_ref().expect("admitted domain").initial_holes;
    PlateWithHoles {
        centers: holes.iter().map(|h| h.center).collect(),
        radii: holes.iter().map(|h| h.radius).collect(),
    }
}
fn budget(value: Option<&str>) -> Result<Option<usize>> {
    value
        .map(|v| {
            v.parse::<usize>()
                .ok()
                .filter(|n| *n > 0 && *n <= 256)
                .ok_or_else(|| {
                    fail(
                        "cli-study-budget-override",
                        "--budget must be an integer in 1..=256",
                    )
                })
        })
        .transpose()
}
fn ir(id: ContentHash, ordinal: usize) -> String {
    format!(
        "{{\"driver\":{DRIVER:?},\"study_id\":\"{}\",\"ordinal\":{ordinal},\"units\":\"1\"}}",
        id.to_hex()
    )
}

#[derive(Debug)]
struct Outcome {
    pointer: String,
    receipt: String,
    status: &'static str,
}
fn render_output(out: Outcome, mode: OutputMode) -> CommandOutput {
    let exit_code = match out.status {
        "completed" => exit::SUCCESS,
        "cancelled" => exit::CANCELLED,
        _ => exit::BUDGET,
    };
    let stdout = match mode {
        OutputMode::Json => format!(
            "{{\"command\":\"study\",\"status\":{0:?},\"run_id\":{1:?},\"run\":{1:?},\"receipt\":{2}}}\n",
            out.status, out.pointer, out.receipt
        ),
        OutputMode::Text => format!(
            "command=study\nstatus={}\nrun={}\nauthority=estimated\n",
            out.status, out.pointer
        ),
    };
    CommandOutput {
        exit_code,
        stdout,
        stderr: String::new(),
    }
}

pub(crate) fn study_path(
    path: &Path,
    ledger_path: &Path,
    override_text: Option<&str>,
    mode: OutputMode,
) -> CommandOutput {
    let result = (|| {
        let cap = budget(override_text)?;
        let json = match path.extension().and_then(|s| s.to_str()) {
            Some("fsim") => false,
            Some("json") => true,
            _ => {
                return Err(fail(
                    "cli-study-format",
                    "study inputs require .fsim or .json",
                ));
            }
        };
        let mut bytes = Vec::new();
        std::fs::File::open(path)
            .and_then(|f| f.take(MAX_PROJECT_BYTES + 1).read_to_end(&mut bytes))
            .map_err(|e| Failure {
                code: "cli-study-read",
                message: e.to_string(),
                exit: exit::INPUT,
            })?;
        if bytes.len() as u64 > MAX_PROJECT_BYTES {
            return Err(Failure {
                code: "cli-study-size",
                message: "study exceeds 16 MiB".into(),
                exit: exit::INPUT,
            });
        }
        let text = std::str::from_utf8(&bytes).map_err(|e| Failure {
            code: "cli-study-utf8",
            message: e.to_string(),
            exit: exit::INPUT,
        })?;
        let admitted = admit(text, json)?;
        let ledger = Ledger::open(
            ledger_path
                .to_str()
                .ok_or_else(|| fail("cli-study-ledger-path", "ledger path is not UTF-8"))?,
        )?;
        let start = Instant::now();
        drive(
            &admitted,
            &ledger,
            cap,
            &CancelGate::new(),
            &|| start.elapsed().as_secs_f64(),
            None,
        )
    })();
    match result {
        Ok(out) => render_output(out, mode),
        Err(e) => failure_output("study", mode, e),
    }
}

fn drive(
    a: &Admitted,
    ledger: &Ledger,
    cap: Option<usize>,
    gate: &CancelGate,
    clock: &dyn Fn() -> f64,
    prior: Option<&Loaded>,
) -> Result<Outcome> {
    if gate.is_requested() {
        return Err(Failure {
            code: "cli-study-cancelled",
            message: "cancelled before publication".into(),
            exit: exit::CANCELLED,
        });
    }
    if ledger.in_transaction() {
        return Err(fail(
            "cli-study-transaction",
            "study requires its own ledger transaction",
        ));
    }
    let mut runner = StudyRunner::new(design(&a.spec), a.config.clone())
        .map_err(|e| fail("cli-study-geometry", e.to_string()))?;
    let mut used_wall = 0.0;
    let mut predecessor = None;
    if let Some(old) = prior {
        if old.value.str_field("study_id") != Some(a.id.to_hex().as_str()) {
            return Err(fail(
                "cli-study-resume-identity",
                "retained study identity changed",
            ));
        }
        let count = integer(&old.value, "iterations_completed")?;
        if count > a.config.steps {
            return Err(fail(
                "cli-study-resume-identity",
                "retained iteration count exceeds requested steps",
            ));
        }
        used_wall = old
            .value
            .f64_field("consumed_wall_s")
            .filter(|v| v.is_finite() && *v >= 0.0)
            .ok_or_else(|| fail("cli-study-resume-budget", "invalid retained wall charge"))?;
        let mut replay_clock = 0.0;
        for _ in 0..count {
            if gate.is_requested() {
                return Err(Failure {
                    code: "cli-study-cancelled",
                    message: "cancelled while replaying retained prefix; no publication".into(),
                    exit: exit::CANCELLED,
                });
            }
            let now = clock();
            if !now.is_finite() || now < replay_clock {
                return Err(fail(
                    "cli-study-clock",
                    "replay clock must be finite and monotonic",
                ));
            }
            replay_clock = now;
            if now >= a.wall_s - used_wall {
                return Err(Failure {
                    code: "cli-study-resume-budget",
                    message: format!(
                        "wall budget exhausted during prefix replay; retained checkpoint: study-{}",
                        old.hash.to_hex()
                    ),
                    exit: exit::BUDGET,
                });
            }
            runner
                .advance()
                .map_err(|e| fail("cli-study-replay", e.to_string()))?;
        }
        if old.value.str_field("trace_hash") != Some(runner.report().trace_hash.as_str()) {
            return Err(fail(
                "cli-study-replay",
                "replayed prefix does not reproduce the retained trace",
            ));
        }
        predecessor = Some(old.hash);
    }
    let governor = Governor::new();
    let mut session_bytes = [0; 8];
    session_bytes.copy_from_slice(&a.id.as_bytes()[..8]);
    let session = SessionId(u64::from_le_bytes(session_bytes));
    let token = CapabilityToken {
        session,
        ops: vec!["optimization.thermal-radii".into()],
        core_s: a.wall_s,
        mem_bytes: a.memory_bytes,
        wall_s: a.wall_s,
        cores: 1,
        ledger_scope: a.id.to_hex(),
    };
    let open = governor
        .session_open_id(session, "study-open")
        .map_err(|e| fail("cli-study-session", e.to_string()))?;
    governor
        .open_session_declared(open, token)
        .map_err(|e| fail("cli-study-session", e.to_string()))?;
    let charge = |ordinal: usize, seconds: f64| -> Result<bool> {
        let key = governor
            .meter_report_id(session, &format!("study-meter-{ordinal}"))
            .map_err(|e| fail("cli-study-session", e.to_string()))?;
        let receipt = governor
            .charge(
                key,
                Charge {
                    core_s: seconds,
                    wall_s: seconds,
                    mem_peak_bytes: 0,
                },
            )
            .map_err(|e| fail("cli-study-session", e.to_string()))?;
        Ok(!matches!(receipt.enforcement(), Enforcement::Ok))
    };
    let start_count = runner.iterations().len();
    let limit = a
        .spec
        .budgets
        .as_ref()
        .and_then(|b| b.max_iterations)
        .expect("admitted budget")
        .min(a.config.steps);
    let mut previous_clock = 0.0;
    let mut exhausted = charge(0, used_wall)?;
    loop {
        let now = clock();
        if !now.is_finite() || now < previous_clock {
            return Err(fail(
                "cli-study-clock",
                "study clock must be finite and monotonic",
            ));
        }
        let delta = now - previous_clock;
        previous_clock = now;
        let next_wall = used_wall + delta;
        if !next_wall.is_finite() {
            return Err(fail(
                "cli-study-clock",
                "accumulated wall charge overflowed",
            ));
        }
        used_wall = next_wall;
        exhausted |= charge(runner.iterations().len() + 1, delta)?;
        let n = runner.iterations().len();
        let status = if gate.is_requested() {
            "cancelled"
        } else if exhausted {
            "budget-exhausted"
        } else if n == a.config.steps {
            "completed"
        } else if n >= limit || cap.is_some_and(|c| n - start_count >= c) {
            "budget-exhausted"
        } else {
            "running"
        };
        let out = persist(a, ledger, &runner, status, used_wall, predecessor)?;
        predecessor = ContentHash::from_hex(out.pointer.trim_start_matches("study-"));
        if status != "running" {
            return Ok(out);
        }
        runner.advance().map_err(|e| {
            fail(
                "cli-study-solve",
                format!("{e}; last durable checkpoint: {}", out.pointer),
            )
        })?;
    }
}

fn persist(
    a: &Admitted,
    ledger: &Ledger,
    runner: &StudyRunner,
    status: &'static str,
    wall_s: f64,
    predecessor: Option<ContentHash>,
) -> Result<Outcome> {
    let report = runner.report();
    let n = report.iterations.len();
    let mut rows = format!(
        "{{\"schema\":\"thermal-study-iterations-v1\",\"study_id\":\"{}\",\"iterations\":[",
        a.id.to_hex()
    );
    for (index, row) in report.iterations.iter().enumerate() {
        if index > 0 {
            rows.push(',');
        }
        rows.push_str(row.jsonl_row().trim());
    }
    rows.push_str("]}");
    let design_json = format!(
        "{{\"model\":\"unit-plate-circular-holes\",\"centers\":{:?},\"radii\":{:?},\"area\":{:.17e}}}",
        report.design.centers,
        report.design.radii,
        report.design.area()
    );
    let mut package = EvidencePackage::new(Provenance::new(
        format!("fs-cli/{};{DRIVER}", env!("CARGO_PKG_VERSION")),
        hash_bytes(include_bytes!("../../../constellation.lock")).to_hex(),
    ));
    let final_compliance = report.iterations.last().map_or("null".to_string(), |r| {
        format!("{:.17e}", r.accepted_compliance)
    });
    if let Some(last) = report.iterations.last() {
        package = package.with_claim(Claim::estimated("thermal.accepted-compliance",
            format!("Accepted normalized compliance {} after {} transitions; source {}; trace {}; status {}. {}", final_compliance, n, a.id.to_hex(), report.trace_hash, status, NO_CLAIM),
            "fs-marquee/dwr-plus-algebraic", last.accepted_cert_dwr + last.accepted_cert_algebraic));
    }
    let package_json = package
        .to_json()
        .map_err(|e| fail("cli-study-package", e.to_string()))?;
    if !fs_checker::check(&package).passed() {
        return Err(fail(
            "cli-study-package",
            "checker refused the produced package",
        ));
    }
    let summary = format!(
        "{{\"study_id\":\"{}\",\"status\":{status:?},\"iterations_completed\":{n},\"target_iterations\":{},\"final_compliance\":{final_compliance},\"final_area\":{:.17e},\"trace_hash\":\"{}\",\"authority\":\"Estimated\",\"no_claim\":{}}}",
        a.id.to_hex(),
        a.config.steps,
        report.design.area(),
        report.trace_hash,
        quoted(NO_CLAIM)
    );
    let mut table = String::new();
    let mut points = String::new();
    let scale = report
        .iterations
        .first()
        .map_or(1.0, |r| r.compliance.abs().max(f64::MIN_POSITIVE));
    for (i, r) in report.iterations.iter().enumerate() {
        let _ = writeln!(
            table,
            "<tr><td>{}</td><td>{:.8e}</td><td>{:.5e}</td><td>{:.5e}</td><td>{}</td><td>Estimated</td></tr>",
            i + 1,
            r.accepted_compliance,
            r.accepted_cert_dwr,
            r.accepted_cert_algebraic,
            r.backtracks
        );
        let _ = write!(
            points,
            "{},{} ",
            20.0 + 560.0 * i as f64 / n.max(1) as f64,
            180.0 - 160.0 * r.accepted_compliance / scale
        );
    }
    let html = format!(
        "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Thermal radius study</title><body><h1>Normalized thermal radius study</h1><p>Status: {status}. {n}/{} iterations. Estimated.</p><p>{NO_CLAIM}</p><p>Accepted compliance: {final_compliance}; material area: {:.8}; target: {:.8}.</p><svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 600 200\" role=\"img\" aria-label=\"Accepted normalized compliance by iteration\"><polyline fill=\"none\" stroke=\"#2563eb\" points=\"{points}\"/></svg><table><tr><th>Iteration</th><th>Accepted compliance</th><th>DWR estimate</th><th>Algebraic estimate</th><th>Backtracks</th><th>Authority</th></tr>{table}</table><p>Trace: {}</p></body></html>",
        a.config.steps,
        report.design.area(),
        a.config.area_target,
        report.trace_hash
    );
    let versions = format!(
        "{{\"driver\":{DRIVER:?},\"crate\":{:?},\"constellation_lock\":\"{}\"}}",
        env!("CARGO_PKG_VERSION"),
        hash_bytes(include_bytes!("../../../constellation.lock")).to_hex()
    );
    let budget_json = format!(
        "{{\"wall_s\":{},\"memory_bytes\":{},\"max_iterations\":{},\"consumed_wall_s\":{wall_s}}}",
        a.wall_s,
        a.memory_bytes,
        a.spec
            .budgets
            .as_ref()
            .and_then(|b| b.max_iterations)
            .expect("budget")
    );
    let seed = a.spec.seeds.as_ref().expect("seed").root.to_le_bytes();
    ledger.begin()?;
    let result = (|| -> Result<Outcome> {
        let op = ledger.begin_op(Some(a.id.as_bytes()), &ir(a.id, n), &FiveExplicits { seed: &seed, versions: &versions, budget: &budget_json, capability: "{\"ops\":[\"optimization.thermal-radii\",\"geometry.sdf\",\"physics.cutfem\"]}" }, 0)?;
        if let Some(previous) = predecessor {
            ledger.link(op, &previous, EdgeRole::In)?;
        }
        let source = ledger.put_artifact("study-source", a.source.as_bytes(), None)?;
        ledger.link(op, &source.hash, EdgeRole::In)?;
        let mut refs = String::new();
        for (name, kind, bytes) in [
            ("iterations", "study-iterations", rows.as_bytes()),
            ("design", "study-design", design_json.as_bytes()),
            ("report_html", "study-report-html", html.as_bytes()),
            ("report_json", "study-report-json", summary.as_bytes()),
            ("package", "study-package", package_json.as_bytes()),
        ] {
            let artifact = ledger.put_artifact(kind, bytes, None)?;
            ledger.link(op, &artifact.hash, EdgeRole::Out)?;
            let _ = write!(refs, ",{name:?}:\"{}\"", artifact.hash.to_hex());
        }
        let receipt = format!(
            "{{\"schema\":{STUDY_RUN_RECEIPT_SCHEMA:?},\"study_id\":\"{}\",\"status\":{status:?},\"source\":\"{}\",\"iterations_completed\":{n},\"target_iterations\":{},\"trace_hash\":\"{}\",\"consumed_wall_s\":{wall_s},\"predecessor\":{},\"stages\":[\"study-admit\",{},\"study-report\",\"study-package\"]{refs}}}",
            a.id.to_hex(),
            source.hash.to_hex(),
            a.config.steps,
            report.trace_hash,
            predecessor.map_or("null".into(), |h| quoted(&h.to_hex())),
            if n == 0 {
                "\"study-optimize-not-started\""
            } else {
                "\"study-optimize\""
            }
        );
        let artifact = ledger.put_artifact(RECEIPT_KIND, receipt.as_bytes(), None)?;
        ledger.link(op, &artifact.hash, EdgeRole::Out)?;
        if ledger.artifact_output_seal(&artifact.hash)?.is_none() {
            ledger.seal_artifact_output(&artifact.hash, op)?;
        }
        ledger.finish_op(op, OpOutcome::Ok, None, 1)?;
        Ok(Outcome {
            pointer: format!("study-{}", artifact.hash.to_hex()),
            receipt,
            status,
        })
    })();
    match result {
        Ok(out) => {
            if let Err(e) = ledger.commit() {
                return match ledger.rollback() {
                    Ok(()) => Err(e.into()),
                    Err(r) => Err(fail(
                        "cli-study-ledger",
                        format!("commit failed: {e}; rollback also failed: {r}"),
                    )),
                };
            }
            Ok(out)
        }
        Err(e) => {
            let rollback = ledger.rollback();
            match rollback {
                Ok(()) => Err(e),
                Err(r) => Err(fail(
                    "cli-study-ledger",
                    format!("{}; rollback also failed: {r}", e.message),
                )),
            }
        }
    }
}

struct Loaded {
    hash: ContentHash,
    value: JsonValue,
}
fn integer(value: &JsonValue, key: &str) -> Result<usize> {
    value
        .get(key)
        .and_then(JsonValue::number_raw)
        .and_then(|s| s.parse::<usize>().ok())
        .ok_or_else(|| fail("cli-study-receipt", format!("invalid {key}")))
}
fn artifact(ledger: &Ledger, hash: ContentHash, kind: &str) -> Result<Vec<u8>> {
    let info = ledger
        .artifact_info(&hash)?
        .ok_or_else(|| fail("cli-study-artifact", "missing study artifact"))?;
    if info.kind != kind {
        return Err(fail(
            "cli-study-artifact",
            format!("expected {kind}, found {}", info.kind),
        ));
    }
    ledger
        .get_artifact_bounded(&hash, MAX_ARTIFACT_BYTES)?
        .ok_or_else(|| fail("cli-study-artifact", "missing study bytes"))
}
fn linked(ledger: &Ledger, value: &JsonValue, key: &str, kind: &str) -> Result<Vec<u8>> {
    let hash = value
        .str_field(key)
        .and_then(ContentHash::from_hex)
        .ok_or_else(|| fail("cli-study-receipt", format!("missing {key} hash")))?;
    artifact(ledger, hash, kind)
}
fn load(ledger: &Ledger, pointer: &str) -> Result<Loaded> {
    let hash = pointer
        .strip_prefix("study-")
        .and_then(ContentHash::from_hex)
        .ok_or_else(|| {
            fail(
                "cli-study-run-id",
                "expected study- followed by a 64-digit receipt hash",
            )
        })?;
    let bytes = artifact(ledger, hash, RECEIPT_KIND)?;
    let text = std::str::from_utf8(&bytes).map_err(|e| fail("cli-study-receipt", e.to_string()))?;
    let value = JsonValue::parse(text).map_err(|e| fail("cli-study-receipt", e.to_string()))?;
    if value.str_field("schema") != Some(STUDY_RUN_RECEIPT_SCHEMA) {
        return Err(fail("cli-study-receipt", "unsupported study receipt"));
    }
    let study_id = value
        .str_field("study_id")
        .and_then(ContentHash::from_hex)
        .ok_or_else(|| fail("cli-study-receipt", "missing study identity"))?;
    let ordinal = integer(&value, "iterations_completed")?;
    let producer = ledger
        .artifact_output_seal(&hash)?
        .ok_or_else(|| fail("cli-study-receipt", "unsealed study receipt"))?;
    let op = ledger
        .op(producer)?
        .ok_or_else(|| fail("cli-study-receipt", "missing producing operation"))?;
    if op.session.as_deref() != Some(study_id.as_bytes().as_slice())
        || op.ir != ir(study_id, ordinal)
        || op.outcome.as_deref() != Some("ok")
        || !ledger.edge_exists(producer, &hash, EdgeRole::Out)?
    {
        return Err(fail(
            "cli-study-receipt",
            "receipt is not bound to the completed study operation",
        ));
    }
    for (key, kind, role) in [
        ("source", "study-source", EdgeRole::In),
        ("iterations", "study-iterations", EdgeRole::Out),
        ("design", "study-design", EdgeRole::Out),
        ("report_html", "study-report-html", EdgeRole::Out),
        ("report_json", "study-report-json", EdgeRole::Out),
        ("package", "study-package", EdgeRole::Out),
    ] {
        let h = value
            .str_field(key)
            .and_then(ContentHash::from_hex)
            .ok_or_else(|| fail("cli-study-receipt", format!("missing {key}")))?;
        if !ledger.edge_exists(producer, &h, role)? {
            return Err(fail("cli-study-receipt", format!("missing {key} lineage")));
        }
        artifact(ledger, h, kind)?;
    }
    Ok(Loaded { hash, value })
}

pub(crate) fn resume_path(
    pointer: &str,
    path: &Path,
    override_text: Option<&str>,
    mode: OutputMode,
) -> CommandOutput {
    let result = (|| {
        let cap = budget(override_text)?;
        if !path.is_file() {
            return Err(fail(
                "cli-study-ledger-missing",
                "resume requires an existing ledger",
            ));
        }
        let ledger = Ledger::open(
            path.to_str()
                .ok_or_else(|| fail("cli-study-ledger-path", "ledger path is not UTF-8"))?,
        )?;
        let old = load(&ledger, pointer)?;
        let source = linked(&ledger, &old.value, "source", "study-source")?;
        let text =
            std::str::from_utf8(&source).map_err(|e| fail("cli-study-receipt", e.to_string()))?;
        let a = admit(text, false)?;
        let start = Instant::now();
        drive(
            &a,
            &ledger,
            cap,
            &CancelGate::new(),
            &|| start.elapsed().as_secs_f64(),
            Some(&old),
        )
    })();
    match result {
        Ok(out) => render_output(out, mode),
        Err(e) => failure_output("study", mode, e),
    }
}

pub(crate) fn export(
    command: &'static str,
    pointer: &str,
    path: Option<&Path>,
    mode: OutputMode,
) -> CommandOutput {
    let result = (|| -> Result<String> {
        let path = path.filter(|p| p.is_file()).ok_or_else(|| {
            fail(
                "cli-study-ledger-missing",
                "study export requires an existing ledger operand",
            )
        })?;
        let ledger = Ledger::open(
            path.to_str()
                .ok_or_else(|| fail("cli-study-ledger-path", "ledger path is not UTF-8"))?,
        )?;
        let loaded = load(&ledger, pointer)?;
        let dir = path.parent().unwrap_or_else(|| Path::new("."));
        let mut paths = String::new();
        let fields: &[(&str, &str, &str)] = if command == "package" {
            &[("package", "study-package", "fspkg")]
        } else {
            &[
                ("report_html", "study-report-html", "html"),
                ("report_json", "study-report-json", "json"),
            ]
        };
        for &(key, kind, extension) in fields {
            let bytes = linked(&ledger, &loaded.value, key, kind)?;
            if key == "package" {
                let text = std::str::from_utf8(&bytes)
                    .map_err(|e| fail("cli-study-package", e.to_string()))?;
                let package = EvidencePackage::from_json(text)
                    .map_err(|e| fail("cli-study-package", e.to_string()))?;
                if !fs_checker::check(&package).passed() {
                    return Err(fail(
                        "cli-study-package",
                        "retained package failed structural verification",
                    ));
                }
            }
            let dest = dir.join(format!("{pointer}.{extension}"));
            crate::report::write_retained(&dest, &bytes)
                .map_err(|e| fail("cli-study-export", e))?;
            let _ = write!(paths, ",{key:?}:{}", quoted(&dest.to_string_lossy()));
        }
        Ok(format!(
            "{{\"command\":{command:?},\"status\":\"ok\",\"run\":{pointer:?},\"study_status\":{},\"authority\":\"projection-of-retained-estimates\",\"verification\":\"sealed-evidence\"{}{paths}}}\n",
            quoted(loaded.value.str_field("status").unwrap_or("unknown")),
            if command == "package" {
                ",\"checker\":\"pass\",\"checker_authority\":\"structural-integrity-only\""
            } else {
                ""
            }
        ))
    })();
    match result {
        Ok(mut stdout) => {
            if matches!(mode, OutputMode::Text) {
                let value = JsonValue::parse(&stdout).expect("generated export JSON");
                stdout = format!(
                    "command={command}\nstatus=ok\nrun={pointer}\nauthority=projection-of-retained-estimates\n"
                );
                for field in [
                    "study_status",
                    "report_html",
                    "report_json",
                    "package",
                    "checker",
                    "checker_authority",
                ] {
                    if let Some(value) = value.str_field(field) {
                        let _ = writeln!(stdout, "{field}={}", crate::escape_text(value));
                    }
                }
            }
            CommandOutput {
                exit_code: exit::SUCCESS,
                stdout,
                stderr: String::new(),
            }
        }
        Err(e) => failure_output(command, mode, e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;
    use std::sync::atomic::{AtomicUsize, Ordering};

    const FIXTURE: &str = include_str!("../../../examples/marquee/thermal-2d.fsim");
    fn short_source() -> String {
        FIXTURE
            .replace(":steps 8", ":steps 2")
            .replace(":mesh-level 4", ":mesh-level 3")
    }
    fn scratch() -> std::path::PathBuf {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "fs-cli-study-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&path).expect("scratch directory");
        path
    }
    fn value(out: &CommandOutput) -> JsonValue {
        JsonValue::parse(&out.stdout)
            .unwrap_or_else(|e| panic!("{e}: {} / {}", out.stdout, out.stderr))
    }
    fn cli(args: &[&str]) -> CommandOutput {
        crate::run(&args.iter().map(|s| s.to_string()).collect::<Vec<_>>())
    }

    #[test]
    fn g0_executable_study_refuses_lost_declarations_and_wrong_units() {
        let admitted = admit(FIXTURE, false).expect("complete thermal declaration");
        let json = fs_project::study::print_study_json(&admitted.spec);
        assert_eq!(admit(&json, true).expect("JSON study").id, admitted.id);
        for (from, to, code) in [
            (":schema 1", ":schema 2", "cli-study-version"),
            (
                ":unit \"1\"",
                ":unit \"J\"",
                "study-objective-dimension-mismatch",
            ),
            (
                ":bounds ((0 0) (1 1))",
                ":bounds ((0 0) (2 1))",
                "study-noncanonical-declaration",
            ),
            (
                ":memory 1073741824 B",
                ":memory -1 B",
                "study-noncanonical-declaration",
            ),
            (
                ":wall-time 600 s",
                ":wall-time 600 kg",
                "study-noncanonical-declaration",
            ),
            (
                ":storage \"normalized\"",
                ":storage \"SI\"",
                "cli-study-model",
            ),
        ] {
            let bad = FIXTURE.replace(from, to);
            let e = match admit(&bad, false) {
                Err(e) => e,
                Ok(_) => panic!("admitted {to}"),
            };
            assert_eq!(e.code, code, "{}", e.message);
        }
        let bad = FIXTURE.replace("(0.7 0.5)", "(0.31 0.5)");
        assert!(matches!(
            admit(&bad, false),
            Err(Failure {
                code: "cli-study-geometry",
                ..
            })
        ));
    }

    #[test]
    fn g3_cli_budget_resume_report_and_package_use_real_iterates() {
        let dir = scratch();
        let source = dir.join("study.fsim");
        let db = dir.join("ledger.db");
        std::fs::write(&source, short_source()).expect("fixture");
        let partial = cli(&[
            "--json",
            "study",
            source.to_str().unwrap(),
            db.to_str().unwrap(),
            "--budget",
            "1",
        ]);
        assert_eq!(partial.exit_code, exit::BUDGET, "{}", partial.stderr);
        let first = value(&partial);
        assert_eq!(
            first
                .path(&["receipt", "iterations_completed"])
                .and_then(JsonValue::as_f64),
            Some(1.0)
        );
        let pointer = first.str_field("run_id").unwrap();
        let resumed = cli(&["--json", "study", "--resume", pointer, db.to_str().unwrap()]);
        assert_eq!(resumed.exit_code, exit::SUCCESS, "{}", resumed.stderr);
        let last = value(&resumed);
        let pointer = last.str_field("run_id").unwrap();
        let report = cli(&["--json", "report", pointer, db.to_str().unwrap()]);
        assert_eq!(report.exit_code, exit::SUCCESS, "{}", report.stderr);
        let report_json = value(&report);
        let summary =
            std::fs::read_to_string(report_json.str_field("report_json").unwrap()).unwrap();
        let summary = JsonValue::parse(&summary).unwrap();
        assert!(summary.f64_field("final_compliance").unwrap() > 0.0);
        assert!((summary.f64_field("final_area").unwrap() - 0.853).abs() < 1e-12);
        let html_path = report_json.str_field("report_html").unwrap();
        let html = std::fs::read_to_string(html_path).unwrap();
        assert!(html.contains("<polyline") && html.contains("DWR estimate"));
        let exported = cli(&["--json", "package", pointer, db.to_str().unwrap()]);
        assert_eq!(exported.exit_code, exit::SUCCESS, "{}", exported.stderr);
        let exported_json = value(&exported);
        let bytes = std::fs::read_to_string(exported_json.str_field("package").unwrap()).unwrap();
        let package = EvidencePackage::from_json(&bytes).unwrap();
        assert!(fs_checker::check(&package).passed());
        assert_eq!(package.declared_claims_unverified().len(), 1);
        let full_db = dir.join("full.db");
        let full = cli(&[
            "--json",
            "study",
            source.to_str().unwrap(),
            full_db.to_str().unwrap(),
        ]);
        assert_eq!(full.exit_code, exit::SUCCESS, "{}", full.stderr);
        assert_eq!(
            value(&full).path(&["receipt", "trace_hash"]),
            last.path(&["receipt", "trace_hash"])
        );
        let ledger = Ledger::open(db.to_str().unwrap()).unwrap();
        assert!(ledger.lint().unwrap().is_clean());
        std::fs::write(html_path, "user-owned conflicting bytes").unwrap();
        let conflict = cli(&["--json", "report", pointer, db.to_str().unwrap()]);
        assert_eq!(conflict.exit_code, exit::REFUSED);
        assert_eq!(
            std::fs::read_to_string(html_path).unwrap(),
            "user-owned conflicting bytes"
        );
    }

    #[test]
    fn g4_cancel_budget_and_failed_transaction_do_not_publish_success() {
        let a = admit(&short_source(), false).unwrap();
        let ledger = Ledger::open(":memory:").unwrap();
        let gate = CancelGate::new_clock_free();
        gate.request();
        assert_eq!(
            drive(&a, &ledger, None, &gate, &|| 0.0, None)
                .unwrap_err()
                .exit,
            exit::CANCELLED
        );
        assert_eq!(ledger.table_count("ops").unwrap(), 0);
        let gate = CancelGate::new_clock_free();
        let ticks = Cell::new(0);
        let clock = || {
            let n = ticks.get();
            ticks.set(n + 1);
            if n == 1 {
                gate.request();
            }
            n as f64
        };
        let cancelled = drive(&a, &ledger, None, &gate, &clock, None).unwrap();
        assert_eq!(cancelled.status, "cancelled");
        assert_eq!(
            integer(
                &load(&ledger, &cancelled.pointer).unwrap().value,
                "iterations_completed"
            )
            .unwrap(),
            1
        );
        let exhausted = drive(
            &a,
            &ledger,
            None,
            &CancelGate::new_clock_free(),
            &|| 601.0,
            None,
        )
        .unwrap();
        assert_eq!(exhausted.status, "budget-exhausted");
        let checkpoint = load(&ledger, &cancelled.pointer).unwrap();
        let before_replay = ledger.table_count("ops").unwrap();
        let e = drive(
            &a,
            &ledger,
            None,
            &CancelGate::new_clock_free(),
            &|| 601.0,
            Some(&checkpoint),
        )
        .unwrap_err();
        assert_eq!(e.code, "cli-study-resume-budget");
        assert_eq!(e.exit, exit::BUDGET);
        assert!(e.message.contains(&cancelled.pointer));
        assert_eq!(ledger.table_count("ops").unwrap(), before_replay);
        let before = ledger.table_count("ops").unwrap();
        ledger.begin().unwrap();
        let e = drive(
            &a,
            &ledger,
            None,
            &CancelGate::new_clock_free(),
            &|| 0.0,
            None,
        )
        .unwrap_err();
        assert_eq!(e.code, "cli-study-transaction");
        assert_eq!(ledger.table_count("ops").unwrap(), before);
        ledger.rollback().unwrap();
        let bad = scratch().join("missing-parent").join("ledger.db");
        let source = scratch().join("study.fsim");
        std::fs::write(&source, short_source()).unwrap();
        let refused = cli(&[
            "--json",
            "study",
            source.to_str().unwrap(),
            bad.to_str().unwrap(),
        ]);
        assert_ne!(refused.exit_code, exit::SUCCESS);
    }
}
