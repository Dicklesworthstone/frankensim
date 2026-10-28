//! Native delivered-power sizing through the ordinary retained cooling solver.
//! Every decision uses a completed physical candidate. The one-coordinate
//! bracket assumes nondecreasing temperature with delivered power; neither a
//! small bracket nor an accepted sampled peak proves global/continuum safety.

use super::*;
use std::fmt::Write as _;
use fs_project::uncertainty::{DesignPolicy, UniformParameter};

pub(super) const DRIVER: &str = "native-cooling-power-design-v1";
const OBSERVATIONS: &str = "native-power-design-observations";
const PROJECT: &str = "solve-project-source";
const SCOPE: &str = "Maximize one region's delivered watts over the declared interval under the explicit nondecreasing-temperature-with-power assumption. Every feasibility decision uses a completed native solve and the unchanged requirement limit minus margin. Transient decisions use the accepted fine-grid sampled peak including the original initial state; final native QoI is retained separately. The selected project is the highest-power observed feasible candidate. No global monotonicity proof, continuous-time or continuum temperature bound, probability statement, physical validation or compliance signoff is supplied.";

#[derive(Debug, Clone, PartialEq)]
struct Trial {
    sample: Sample,
    temperature_k: f64,
    time_s: Option<f64>,
    state: Option<usize>,
    conduction: ContentHash,
}

impl Trial {
    fn power(&self) -> f64 { self.sample.parameters[0] }
    fn json(&self) -> String {
        format!(
            "{{\"power_w\":{},\"temperature_k\":{},\"final_temperature_k\":{},\"time_s\":{},\"state\":{},\"run\":{},\"project_hash\":{},\"qoi_receipt\":{},\"conduction_receipt\":{}}}",
            self.power(), self.temperature_k, self.sample.value_k, optional(self.time_s),
            self.state.map_or_else(|| "null".into(), |v| v.to_string()),
            quoted(&self.sample.run), quoted(&self.sample.project_hash),
            quoted(&self.sample.qoi_receipt.to_hex()), quoted(&self.conduction.to_hex())
        )
    }
}

fn invalid(message: impl Into<String>) -> Failure {
    fail("cli-native-power-design", message)
}

fn transient(model: &Model) -> bool {
    model.bound.base().cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .is_some_and(|c| c.transient.is_some())
}

/// The child run has already passed the native sealed-QoI/geometry checks.
/// Bind the forward sampled peak to that SAME run and requirement region.
fn observe(model: &Model, ledger: &Ledger, sample: Sample) -> Result<Trial> {
    if sample.parameters.len() != 1 { return Err(invalid("design needs one physical coordinate")); }
    let sealed = crate::solve::load_completed_run(ledger, &sample.run)
        .map_err(|e| invalid(format!("{}: {}", e.code, e.what)))?;
    let conduction = sealed.stages.iter().find(|row| row.0 == "conduction")
        .and_then(|row| ContentHash::from_hex(&row.2))
        .ok_or_else(|| invalid("completed design candidate has no conduction receipt"))?;
    let bytes = artifact(ledger, conduction, "solve-stage-receipt", 4 * 1024 * 1024)?;
    let receipt = parse(&bytes)?;
    if receipt.str_field("run") != Some(sample.run.as_str()) {
        return Err(invalid("design observation belongs to another native run"));
    }
    let (temperature_k, time_s, state) = if transient(model) {
        let time = receipt.get("transient").ok_or_else(|| invalid("missing native trajectory"))?;
        let peak = time.get("sampled_peak").ok_or_else(|| invalid("native forward trajectory has no sampled-peak observation"))?;
        let region = &model.bound.base().requirements.as_ref().expect("bound requirement")[0].region;
        let steps = integer(time, "fine_steps")?;
        let samples = integer(peak, "samples")?;
        let state = integer(peak, "state")?;
        let horizon = time.f64_field("final_time_s").filter(|v| v.is_finite() && *v > 0.0)
            .ok_or_else(|| invalid("invalid retained trajectory horizon"))?;
        let at = peak.f64_field("time_s").filter(|v| v.is_finite() && *v >= 0.0 && *v <= horizon)
            .ok_or_else(|| invalid("invalid retained peak time"))?;
        let value = peak.f64_field("value_k").filter(|v| v.is_finite() && *v >= 0.0)
            .ok_or_else(|| invalid("invalid retained peak temperature"))?;
        if time.str_field("status") != Some("completed")
            || peak.str_field("region") != Some(region.as_str())
            || peak.get("includes_initial") != Some(&J::Bool(true))
            || steps.checked_add(1) != Some(samples) || state >= samples
            || (state == 0) != (at == 0.0) || value < sample.value_k
        {
            return Err(invalid("sampled peak does not cover this candidate's full fine trajectory and requirement region"));
        }
        (value, Some(at), Some(state))
    } else {
        (sample.value_k, None, None)
    };
    Ok(Trial { sample, temperature_k, time_s, state, conduction })
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Next {
    Evaluate(f64),
    Complete(&'static str),
    Refuse(&'static str),
}

fn selected(rows: &[Trial], limit: f64) -> Option<&Trial> {
    rows.iter().filter(|row| row.temperature_k <= limit)
        .max_by(|a, b| a.power().total_cmp(&b.power()))
}

fn bracket(rows: &[Trial], limit: f64) -> Option<(&Trial, &Trial)> {
    let passing = selected(rows, limit)?;
    let failing = rows.iter().filter(|row| row.temperature_k > limit && row.power() > passing.power())
        .min_by(|a, b| a.power().total_cmp(&b.power()))?;
    Some((passing, failing))
}

fn next(control: &UniformParameter, policy: DesignPolicy, limit: f64, rows: &[Trial]) -> Next {
    if rows.is_empty() { return Next::Evaluate(control.low); }
    if rows[0].temperature_k > limit { return Next::Refuse("minimum-power-infeasible"); }
    // This detects a violated assumption in measured candidates, without
    // claiming that an unobserved nonmonotone branch has been excluded.
    for a in rows {
        for b in rows {
            if a.power() < b.power() && a.temperature_k - b.temperature_k > policy.temperature_tolerance_k {
                return Next::Refuse("nonmonotone-response");
            }
        }
    }
    if rows.len() == 1 { return Next::Evaluate(control.high); }
    if rows[1].temperature_k <= limit { return Next::Complete("maximum-feasible"); }
    let Some((passing, failing)) = bracket(rows, limit) else {
        return Next::Refuse("missing-feasible-bracket");
    };
    let width = failing.power() - passing.power();
    if width <= policy.power_tolerance_w && limit - passing.temperature_k <= policy.temperature_tolerance_k {
        return Next::Complete("target-bracketed");
    }
    let candidate = passing.power() + 0.5 * width;
    if !(candidate > passing.power() && candidate < failing.power()) {
        return Next::Refuse("power-resolution-exhausted");
    }
    Next::Evaluate(candidate)
}

fn ir(id: ContentHash, n: usize) -> String {
    format!("{{\"driver\":{DRIVER:?},\"model\":{},\"evaluations_completed\":{n}}}", quoted(&id.to_hex()))
}

#[allow(clippy::too_many_arguments)]
fn persist(
    model: &Model, ledger: &Ledger, rows: &[Trial], status: &str, termination: &str,
    failure: Option<&str>, used_wall: f64, predecessor: Option<ContentHash>,
) -> Result<Outcome> {
    let study = model.bound.study();
    let policy = study.design().ok_or_else(|| invalid("missing deterministic design policy"))?;
    let control = &study.parameters()[0];
    let limit = model.bound.threshold_k();
    let chosen = selected(rows, limit);
    if status == "completed" && chosen.is_none() {
        return Err(invalid("a completed design must retain an evaluated feasible project"));
    }
    let id = model.identity();
    let n = rows.len();
    let history = rows.iter().map(Trial::json).collect::<Vec<_>>().join(",");
    let chosen_json = chosen.map_or_else(|| "null".into(), Trial::json);
    let width = bracket(rows, limit).map(|(a, b)| b.power() - a.power());
    let failure_json = failure.map_or_else(|| "null".into(), quoted);
    let observation = if transient(model) { "fine-grid-sampled-peak-including-initial" } else { "steady" };
    let region = &model.bound.base().requirements.as_ref().expect("bound requirement")[0].region;
    let summary = format!(
        "{{\"schema\":\"frankensim.cli.native-power-design-result.v1\",\"driver\":{DRIVER:?},\"study_id\":{},\"status\":{},\"termination\":{},\"objective\":\"maximize-delivered-regional-power\",\"power_region\":{},\"temperature_region\":{},\"temperature_observation\":{observation:?},\"temperature_limit_k\":{limit},\"power_tolerance_w\":{},\"temperature_tolerance_k\":{},\"evaluations_completed\":{n},\"evaluations_planned\":{},\"selected\":{chosen_json},\"bracket_width_w\":{},\"history\":[{history}],\"failure\":{failure_json},\"authority\":\"Estimated\",\"no_claim\":{}}}",
        quoted(&id.to_hex()), quoted(status), quoted(termination), quoted(&control.entity),
        quoted(region), policy.power_tolerance_w, policy.temperature_tolerance_k, study.samples(),
        optional(width), quoted(SCOPE)
    );
    let chosen_html = chosen.map_or_else(|| "No evaluated feasible design yet.".into(), |row|
        format!("Selected delivered power: {} W; evaluated temperature: {} K; final native QoI: {} K; child run: {}.",
            row.power(), row.temperature_k, row.sample.value_k, row.sample.run));
    let mut table = String::new();
    for row in rows {
        let _ = write!(table, "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
            row.power(), row.temperature_k, row.sample.value_k, row.sample.run);
    }
    let html = format!("<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Native power design</title><body><h1>Native power design</h1><p>Status: {status}; termination: {termination}; {n} completed candidates.</p><p>{chosen_html}</p><p>Temperature limit after margin: {limit} K. Observation: {observation}.</p><table><tr><th>Delivered W</th><th>Decision K</th><th>Final QoI K</th><th>Native run</th></tr>{table}</table><p>{SCOPE}</p></body></html>");
    // A retained feasible numerical result is useful, but no certified optimum
    // or physical-temperature error bound has been established by this search.
    let package = EvidencePackage::new(Provenance::new(
        format!("fs-cli/{}+{DRIVER}", env!("CARGO_PKG_VERSION")), id.to_hex()))
        .to_json().map_err(|e| invalid(e.to_string()))?;
    let samples = rows_json(&rows.iter().map(|row| row.sample.clone()).collect::<Vec<_>>());
    let observations = format!("{{\"driver\":{DRIVER:?},\"evaluations\":{samples}}}");
    let project = chosen.map(|row| model.project(&[row.power()])).transpose()?;
    if project.as_ref().zip(chosen).is_some_and(|(project, row)|
        project.hash().to_hex() != row.sample.project_hash) {
        return Err(invalid("selected export differs from the evaluated native project"));
    }
    let seed = study.seed().to_le_bytes();
    let versions = format!("{{\"driver\":{DRIVER:?},\"solve_driver\":{}}}", crate::SOLVE_DRIVER_VERSION);
    let budgets = format!("{{\"evaluations\":{},\"wall_s\":{},\"consumed_wall_s\":{used_wall}}}",
        study.samples(), study.wall_seconds());
    ledger.begin()?;
    let result = (|| {
        let op = ledger.begin_op(Some(id.as_bytes()), &ir(id, n), &FiveExplicits {
            seed: &seed, versions: &versions, budget: &budgets,
            capability: "{\"ops\":[\"native-delivered-power-design\"],\"physical_ops\":\"inherited-from-retained-base-project\"}",
        }, 0)?;
        if model.retain(ledger, op)? != id { return Err(invalid("retained native model identity differs")); }
        if let Some(previous) = predecessor { ledger.link(op, &previous, EdgeRole::In)?; }
        let mut children = std::collections::BTreeSet::new();
        for row in rows {
            for hash in [row.sample.qoi_receipt, row.conduction] {
                if children.insert(hash) { ledger.link(op, &hash, EdgeRole::In)?; }
            }
        }
        let mut refs = String::new();
        for (key, kind, bytes) in [
            ("observations", OBSERVATIONS, observations.as_bytes()),
            ("report_json", "study-report-json", summary.as_bytes()),
            ("report_html", "study-report-html", html.as_bytes()),
            ("package", "study-package", package.as_bytes()),
        ] {
            let a = ledger.put_artifact(kind, bytes, None)?;
            ledger.link(op, &a.hash, EdgeRole::Out)?;
            let _ = write!(refs, ",{key:?}:{}", quoted(&a.hash.to_hex()));
        }
        let selected_project = if let Some(project) = &project {
            let a = ledger.put_artifact(PROJECT, project.canonical.as_bytes(), None)?;
            ledger.link(op, &a.hash, EdgeRole::Out)?;
            quoted(&a.hash.to_hex())
        } else { "null".into() };
        let receipt = format!(
            "{{\"schema\":{STUDY_RUN_RECEIPT_SCHEMA:?},\"driver\":{DRIVER:?},\"study_id\":{},\"model\":{},\"status\":{},\"termination\":{},\"evaluations_completed\":{n},\"evaluations_planned\":{},\"consumed_wall_s\":{used_wall},\"selected_project\":{selected_project},\"failure\":{failure_json},\"predecessor\":{}{refs}}}",
            quoted(&id.to_hex()), quoted(&id.to_hex()), quoted(status), quoted(termination),
            study.samples(), predecessor.map_or_else(|| "null".into(), |h| quoted(&h.to_hex()))
        );
        let a = ledger.put_artifact(RECEIPT_KIND, receipt.as_bytes(), None)?;
        ledger.link(op, &a.hash, EdgeRole::Out)?;
        if ledger.artifact_output_seal(&a.hash)?.is_none() { ledger.seal_artifact_output(&a.hash, op)?; }
        ledger.finish_op(op, OpOutcome::Ok, None, 1)?;
        Ok(Outcome { pointer: format!("study-{}", a.hash.to_hex()), receipt,
            status: status.into(), authority: "estimated-native-design" })
    })();
    match result {
        Ok(out) => match ledger.commit() {
            Ok(()) => Ok(out),
            Err(error) => { ledger.rollback()?; Err(error.into()) }
        },
        Err(error) => { ledger.rollback()?; Err(error) }
    }
}

pub(super) fn validate_loaded(ledger: &Ledger, hash: ContentHash, value: J, bytes: Vec<u8>) -> Result<Loaded> {
    let id = hash_field(&value, "study_id")?;
    let n = integer(&value, "evaluations_completed")?;
    let op_id = ledger.artifact_output_seal(&hash)?.ok_or_else(|| invalid("unsealed design receipt"))?;
    let op = ledger.op(op_id)?.ok_or_else(|| invalid("missing design producer"))?;
    if value.str_field("schema") != Some(STUDY_RUN_RECEIPT_SCHEMA)
        || hash_field(&value, "model")? != id || n > fs_project::uncertainty::MAX_SAMPLES
        || op.session.as_deref() != Some(id.as_bytes().as_slice())
        || op.ir != ir(id, n) || op.outcome.as_deref() != Some("ok")
        || !ledger.edge_exists(op_id, &hash, EdgeRole::Out)?
    {
        return Err(invalid("design receipt differs from its completed producing operation"));
    }
    for (key, kind, role) in [
        ("model", model::MANIFEST_KIND, EdgeRole::In),
        ("observations", OBSERVATIONS, EdgeRole::Out),
        ("report_json", "study-report-json", EdgeRole::Out),
        ("report_html", "study-report-html", EdgeRole::Out),
        ("package", "study-package", EdgeRole::Out),
    ] {
        let h = hash_field(&value, key)?;
        if !ledger.edge_exists(op_id, &h, role)? { return Err(invalid(format!("missing design {key} lineage"))); }
        artifact(ledger, h, kind, MAX_ARTIFACT_BYTES)?;
    }
    if value.str_field("selected_project").is_some() {
        let h = hash_field(&value, "selected_project")?;
        if !ledger.edge_exists(op_id, &h, EdgeRole::Out)? { return Err(invalid("missing selected-project lineage")); }
        artifact(ledger, h, PROJECT, crate::MAX_PROJECT_BYTES)?;
    }
    Ok(Loaded { hash, value, bytes })
}

pub(super) fn drive(
    model: &Model, ledger: &Ledger, cap: Option<usize>, gate: &CancelGate,
    started: Instant, prior: Option<&Loaded>,
) -> Result<Outcome> {
    let study = model.bound.study();
    let policy = study.design().ok_or_else(|| invalid("missing deterministic design policy"))?;
    let control = &study.parameters()[0];
    let limit = model.bound.threshold_k();
    let mut rows = Vec::new();
    let used_before = if let Some(old) = prior {
        if hash_field(&old.value, "model")? != model.identity() {
            return Err(invalid("retained design model changed"));
        }
        let observations = parse(&linked(ledger, &old.value, "observations", OBSERVATIONS)?)?;
        if observations.str_field("driver") != Some(DRIVER) {
            return Err(invalid("retained observations belong to another study producer"));
        }
        let samples = read_rows_value(observations.get("evaluations")
            .ok_or_else(|| invalid("missing retained design evaluations"))?)?;
        if samples.len() != integer(&old.value, "evaluations_completed")? || samples.len() > study.samples() {
            return Err(invalid("retained design evaluation count differs"));
        }
        for sample in samples {
            let Next::Evaluate(expected) = next(control, policy, limit, &rows) else {
                return Err(invalid("retained design continues beyond a terminal search"));
            };
            if sample.parameters.len() != 1 || sample.parameters[0].to_bits() != expected.to_bits() {
                return Err(invalid("retained design coordinate differs from the deterministic bracket"));
            }
            model.verify_sample(ledger, &sample)?;
            rows.push(observe(model, ledger, sample)?);
        }
        old.value.f64_field("consumed_wall_s").filter(|v| v.is_finite() && *v >= 0.0)
            .ok_or_else(|| invalid("invalid retained design wall charge"))?
    } else { 0.0 };
    let start_count = rows.len();
    let mut predecessor = prior.map(|old| old.hash);
    let mut interrupted = false;
    let mut failure = None::<String>;
    loop {
        let used = used_before + started.elapsed().as_secs_f64();
        let decision = next(control, policy, limit, &rows);
        let (status, termination) = if failure.is_some() {
            ("refused", "child-refused")
        } else if let Next::Complete(reason) = decision {
            ("completed", reason)
        } else if let Next::Refuse(reason) = decision {
            ("refused", reason)
        } else if rows.len() >= study.samples() {
            ("budget-exhausted", "lifetime-evaluation-budget")
        } else if used >= study.wall_seconds() {
            ("budget-exhausted", "wall-budget")
        } else if gate.is_requested() {
            ("cancelled", "cancelled")
        } else if interrupted {
            ("budget-exhausted", "child-budget")
        } else if cap.is_some_and(|cap| rows.len() - start_count >= cap) {
            ("budget-exhausted", "invocation-evaluation-budget")
        } else { ("running", "bracketing") };
        let out = persist(model, ledger, &rows, status, termination, failure.as_deref(), used, predecessor)?;
        if status != "running" { return Ok(out); }
        predecessor = ContentHash::from_hex(out.pointer.trim_start_matches("study-"));
        let Next::Evaluate(power) = decision else { return Err(invalid("missing next design coordinate")); };
        let remaining = (study.wall_seconds() - used_before - started.elapsed().as_secs_f64()).max(0.0);
        match model.sample(ledger, gate, &[power], remaining) {
            Ok(Some(sample)) => match observe(model, ledger, sample) {
                Ok(trial) => rows.push(trial),
                Err(error) => failure = Some(error.to_string()),
            },
            Ok(None) => interrupted = true,
            Err(error) => failure = Some(error.to_string()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_project::uncertainty::Target;
    fn control() -> UniformParameter {
        UniformParameter { name: "delivered-power".into(), target: Target::Power,
            entity: "chip".into(), low: 0.0, high: 16.0 }
    }
    fn policy() -> DesignPolicy {
        DesignPolicy { power_tolerance_w: 0.01, temperature_tolerance_k: 0.01 }
    }
    fn trial(power: f64, temperature: f64) -> Trial {
        Trial { sample: Sample { run: String::new(), project_hash: String::new(),
            qoi_receipt: hash_bytes(b"qoi"), value_k: temperature, parameters: vec![power] },
            temperature_k: temperature, time_s: None, state: None, conduction: hash_bytes(b"conduction") }
    }
    #[test]
    fn g1_power_bracket_keeps_a_measured_feasible_endpoint() {
        let mut rows = Vec::new();
        for _ in 0..32 {
            match next(&control(), policy(), 305.0, &rows) {
                Next::Evaluate(power) => rows.push(trial(power, 300.0 + 2.0 * power)),
                Next::Complete(reason) => {
                    assert_eq!(reason, "target-bracketed");
                    let chosen = selected(&rows, 305.0).unwrap();
                    assert!(chosen.temperature_k <= 305.0);
                    assert!((chosen.power() - 2.5).abs() <= 0.01);
                    return;
                },
                Next::Refuse(reason) => panic!("{reason}"),
            }
        }
        panic!("affine physical response failed to bracket");
    }
    #[test]
    fn g0_infeasible_and_nonmonotone_responses_never_complete() {
        assert_eq!(next(&control(), policy(), 305.0, &[trial(0.0, 306.0)]),
            Next::Refuse("minimum-power-infeasible"));
        assert_eq!(next(&control(), policy(), 305.0, &[trial(0.0, 304.0), trial(16.0, 303.0)]),
            Next::Refuse("nonmonotone-response"));
        assert_eq!(next(&control(), policy(), 305.0, &[trial(0.0, 300.0), trial(16.0, 304.0)]),
            Next::Complete("maximum-feasible"));
    }
}
