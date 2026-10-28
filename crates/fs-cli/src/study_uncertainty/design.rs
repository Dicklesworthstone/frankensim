//! Native power and fan sizing through the ordinary retained cooling solver.
//! Every decision uses a completed physical candidate. The one-coordinate
//! bracket uses an explicitly declared monotone temperature response; neither a
//! small bracket nor an accepted sampled peak proves global/continuum safety.

use super::*;
use std::fmt::Write as _;
use fs_project::uncertainty::{DesignControl, DesignPolicy, UniformParameter};

mod cases;

pub(super) const DRIVER: &str = "native-cooling-power-design-v1";
pub(super) const FAN_DRIVER: &str = "native-cooling-fan-design-v1";
const POWER_CASE_DRIVER: &str = "native-cooling-power-case-design-v1";
const FAN_CASE_DRIVER: &str = "native-cooling-fan-case-design-v1";
const OBSERVATIONS: &str = "native-power-design-observations";
const FAN_OBSERVATIONS: &str = "native-fan-design-observations";
const PROJECT: &str = "solve-project-source";
const SCOPE: &str = "Maximize one region's delivered watts over the declared interval under the explicit nondecreasing-temperature-with-power assumption. Every feasibility decision uses a completed native solve and the unchanged requirement limit minus margin. Transient decisions use the accepted fine-grid sampled peak including the original initial state; final native QoI is retained separately. The selected project is the highest-power observed feasible candidate. No global monotonicity proof, continuous-time or continuum temperature bound, probability statement, physical validation or compliance signoff is supplied.";
const FAN_SCOPE: &str = "Minimize one named fan bank's absolute speed ratio over its declared interval under the explicit nonincreasing-temperature-with-fan-speed assumption. Every candidate recomputes the native hydraulic operating point and conjugate cooling response with the original fan curves, other banks, power, materials and geometry. Every feasibility decision uses a completed native solve and the unchanged requirement limit minus margin. Transient decisions use the accepted fine-grid sampled peak including the original initial state; final native QoI is retained separately. The selected project is the lowest-speed observed feasible candidate. No global monotonicity proof, minimum electrical power or acoustic noise, continuous-time or continuum temperature bound, probability statement, physical validation or compliance signoff is supplied.";

pub(super) fn is_driver(value: Option<&str>) -> bool {
    matches!(value, Some(DRIVER | FAN_DRIVER | POWER_CASE_DRIVER | FAN_CASE_DRIVER))
}

fn fan(policy: DesignPolicy) -> bool {
    matches!(policy.control, DesignControl::FanSpeedRatio { .. })
}

fn driver(policy: DesignPolicy) -> &'static str {
    if fan(policy) { FAN_DRIVER } else { DRIVER }
}

fn study_driver(study: &fs_project::uncertainty::UncertaintyStudy, policy: DesignPolicy) -> &'static str {
    if study.operating_cases().is_empty() { driver(policy) }
    else if fan(policy) { FAN_CASE_DRIVER } else { POWER_CASE_DRIVER }
}

fn observations_kind(driver: &str) -> &'static str {
    match driver {
        FAN_CASE_DRIVER => "native-fan-case-design-observations",
        POWER_CASE_DRIVER => "native-power-case-design-observations",
        FAN_DRIVER => FAN_OBSERVATIONS,
        _ => OBSERVATIONS,
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Trial {
    sample: Sample,
    temperature_k: f64,
    time_s: Option<f64>,
    state: Option<usize>,
    conduction: ContentHash,
}

impl Trial {
    fn coordinate(&self) -> f64 { self.sample.parameters[0] }
    fn json(&self, policy: DesignPolicy) -> String {
        let coordinate = if fan(policy) { "fan_speed_ratio" } else { "power_w" };
        format!(
            "{{{coordinate:?}:{},\"temperature_k\":{},\"final_temperature_k\":{},\"time_s\":{},\"state\":{},\"run\":{},\"project_hash\":{},\"qoi_receipt\":{},\"conduction_receipt\":{}}}",
            self.coordinate(), self.temperature_k, self.sample.value_k, optional(self.time_s),
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
    if sample.parameters.len() != model.bound.study().parameters().len() {
        return Err(invalid("design observation does not match its declared physical coordinates"));
    }
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

fn selected(rows: &[Trial], policy: DesignPolicy, limit: f64) -> Option<&Trial> {
    rows.iter().filter(|row| row.temperature_k <= limit)
        .max_by(|a, b| if fan(policy) { b.coordinate().total_cmp(&a.coordinate()) }
            else { a.coordinate().total_cmp(&b.coordinate()) })
}

fn bracket(rows: &[Trial], policy: DesignPolicy, limit: f64) -> Option<(&Trial, &Trial)> {
    let passing = selected(rows, policy, limit)?;
    let failing = rows.iter().filter(|row| row.temperature_k > limit &&
        if fan(policy) { row.coordinate() < passing.coordinate() }
        else { row.coordinate() > passing.coordinate() })
        .min_by(|a, b| if fan(policy) { b.coordinate().total_cmp(&a.coordinate()) }
            else { a.coordinate().total_cmp(&b.coordinate()) })?;
    Some((passing, failing))
}

fn next(control: &UniformParameter, policy: DesignPolicy, limit: f64, rows: &[Trial]) -> Next {
    let minimize = fan(policy);
    if rows.is_empty() { return Next::Evaluate(if minimize { control.high } else { control.low }); }
    if rows[0].temperature_k > limit {
        return Next::Refuse(if minimize { "maximum-fan-speed-infeasible" } else { "minimum-power-infeasible" });
    }
    // This detects a violated assumption in measured candidates, without
    // claiming that an unobserved nonmonotone branch has been excluded.
    for a in rows {
        for b in rows {
            let increase = if minimize { b.temperature_k - a.temperature_k }
                else { a.temperature_k - b.temperature_k };
            if a.coordinate() < b.coordinate() && increase > policy.temperature_tolerance_k {
                return Next::Refuse("nonmonotone-response");
            }
        }
    }
    if rows.len() == 1 { return Next::Evaluate(if minimize { control.low } else { control.high }); }
    if rows[1].temperature_k <= limit {
        return Next::Complete(if minimize { "minimum-feasible" } else { "maximum-feasible" });
    }
    let Some((passing, failing)) = bracket(rows, policy, limit) else {
        return Next::Refuse("missing-feasible-bracket");
    };
    let (low, high) = if minimize { (failing.coordinate(), passing.coordinate()) }
        else { (passing.coordinate(), failing.coordinate()) };
    let width = high - low;
    if width <= policy.control.tolerance() && limit - passing.temperature_k <= policy.temperature_tolerance_k {
        return Next::Complete("target-bracketed");
    }
    let candidate = low + 0.5 * width;
    if !(candidate > low && candidate < high) {
        return Next::Refuse(if minimize { "fan-speed-resolution-exhausted" } else { "power-resolution-exhausted" });
    }
    Next::Evaluate(candidate)
}

fn ir(driver: &str, id: ContentHash, n: usize) -> String {
    format!("{{\"driver\":{driver:?},\"model\":{},\"evaluations_completed\":{n}}}", quoted(&id.to_hex()))
}

#[allow(clippy::too_many_arguments)]
fn persist(
    model: &Model, ledger: &Ledger, rows: &[Trial], status: &str, termination: &str,
    failure: Option<&str>, used_wall: f64, predecessor: Option<ContentHash>,
) -> Result<Outcome> {
    let study = model.bound.study();
    let policy = study.design().ok_or_else(|| invalid("missing deterministic design policy"))?;
    let fan = fan(policy);
    let driver = study_driver(study, policy);
    let scope = cases::scope(study, if fan { FAN_SCOPE } else { SCOPE });
    let observations_kind = observations_kind(driver);
    let control = &study.parameters()[0];
    let limit = model.bound.threshold_k();
    let candidate_rows = cases::completed(rows, cases::count(study));
    let chosen = selected(&candidate_rows, policy, limit);
    if status == "completed" && chosen.is_none() {
        return Err(invalid("a completed design must retain an evaluated feasible project"));
    }
    let id = model.identity();
    let n = rows.len();
    let history = candidate_rows.iter().map(|row| row.json(policy)).collect::<Vec<_>>().join(",");
    let chosen_json = chosen.map_or_else(|| "null".into(), |row| row.json(policy));
    let width = bracket(&candidate_rows, policy, limit).map(|(a, b)| (b.coordinate() - a.coordinate()).abs());
    let case_json = cases::json(study, rows, policy, chosen);
    let failure_json = failure.map_or_else(|| "null".into(), quoted);
    let observation = if transient(model) { "fine-grid-sampled-peak-including-initial" } else { "steady" };
    let region = &model.bound.base().requirements.as_ref().expect("bound requirement")[0].region;
    let schema = match (fan, study.operating_cases().is_empty()) {
        (true, true) => "frankensim.cli.native-fan-design-result.v1",
        (false, true) => "frankensim.cli.native-power-design-result.v1",
        (true, false) => "frankensim.cli.native-fan-case-design-result.v1",
        (false, false) => "frankensim.cli.native-power-case-design-result.v1",
    };
    let objective = if fan { "minimize-fan-speed-ratio" } else { "maximize-delivered-regional-power" };
    let entity = if fan { "fan_bank" } else { "power_region" };
    let tolerance = if fan { "fan_speed_ratio_tolerance" } else { "power_tolerance_w" };
    let bracket_width = if fan { "bracket_width_ratio" } else { "bracket_width_w" };
    let summary = format!(
        "{{\"schema\":{schema:?},\"driver\":{driver:?},\"study_id\":{},\"status\":{},\"termination\":{},\"objective\":{objective:?},{entity:?}:{},\"temperature_region\":{},\"temperature_observation\":{observation:?},\"temperature_limit_k\":{limit},{tolerance:?}:{},\"temperature_tolerance_k\":{},\"evaluations_completed\":{n},\"evaluations_planned\":{},\"selected\":{chosen_json},{bracket_width:?}:{},\"history\":[{history}],\"failure\":{failure_json}{case_json},\"authority\":\"Estimated\",\"no_claim\":{}}}",
        quoted(&id.to_hex()), quoted(status), quoted(termination), quoted(&control.entity),
        quoted(region), policy.control.tolerance(), policy.temperature_tolerance_k, study.samples(),
        optional(width), quoted(&scope)
    );
    let label = if fan { "fan speed ratio" } else { "delivered power" };
    let unit = if fan { " (absolute ratio)" } else { " W" };
    let title = if fan { "Native fan design" } else { "Native power design" };
    let column = if fan { "Fan speed ratio" } else { "Delivered W" };
    let chosen_html = chosen.map_or_else(|| "No evaluated feasible design yet.".into(), |row|
        format!("Selected {label}: {}{unit}; evaluated temperature: {} K; final native QoI: {} K; child run: {}.",
            row.coordinate(), row.temperature_k, row.sample.value_k, row.sample.run));
    let mut table = String::new();
    for row in candidate_rows.iter() {
        let _ = write!(table, "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
            row.coordinate(), row.temperature_k, row.sample.value_k, row.sample.run);
    }
    let candidate_count = candidate_rows.len();
    let case_html = cases::html(study, rows);
    let html = format!("<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>{title}</title><body><h1>{title}</h1><p>Status: {status}; termination: {termination}; {candidate_count} completed candidates.</p><p>{chosen_html}</p><p>Temperature limit after margin: {limit} K. Observation: {observation}.</p><table><tr><th>{column}</th><th>Decision K</th><th>Final QoI K</th><th>Native run</th></tr>{table}</table><p>{scope}</p>{case_html}</body></html>");
    // A retained feasible numerical result is useful, but no certified optimum
    // or physical-temperature error bound has been established by this search.
    let package = EvidencePackage::new(Provenance::new(
        format!("fs-cli/{}+{driver}", env!("CARGO_PKG_VERSION")), id.to_hex()))
        .to_json().map_err(|e| invalid(e.to_string()))?;
    let samples = rows_json(&rows.iter().map(|row| row.sample.clone()).collect::<Vec<_>>());
    let observations = format!("{{\"driver\":{driver:?},\"evaluations\":{samples}}}");
    let project = chosen.map(|row| model.project(&row.sample.parameters)).transpose()?;
    if project.as_ref().zip(chosen).is_some_and(|(project, row)|
        project.hash().to_hex() != row.sample.project_hash) {
        return Err(invalid("selected export differs from the evaluated native project"));
    }
    let seed = study.seed().to_le_bytes();
    let versions = format!("{{\"driver\":{driver:?},\"solve_driver\":{}}}", crate::SOLVE_DRIVER_VERSION);
    let budgets = format!("{{\"evaluations\":{},\"wall_s\":{},\"consumed_wall_s\":{used_wall}}}",
        study.samples(), study.wall_seconds());
    ledger.begin()?;
    let result = (|| {
        let op = ledger.begin_op(Some(id.as_bytes()), &ir(driver, id, n), &FiveExplicits {
            seed: &seed, versions: &versions, budget: &budgets,
            capability: if fan {
                "{\"ops\":[\"native-fan-speed-design\"],\"physical_ops\":\"inherited-from-retained-base-project\"}"
            } else {
                "{\"ops\":[\"native-delivered-power-design\"],\"physical_ops\":\"inherited-from-retained-base-project\"}"
            },
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
            ("observations", observations_kind, observations.as_bytes()),
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
            "{{\"schema\":{STUDY_RUN_RECEIPT_SCHEMA:?},\"driver\":{driver:?},\"study_id\":{},\"model\":{},\"status\":{},\"termination\":{},\"evaluations_completed\":{n},\"evaluations_planned\":{},\"consumed_wall_s\":{used_wall},\"selected_project\":{selected_project},\"failure\":{failure_json},\"predecessor\":{}{refs}}}",
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
    let driver = value.str_field("driver").filter(|value| is_driver(Some(value)))
        .ok_or_else(|| invalid("unknown native design producer"))?;
    let observations = observations_kind(driver);
    let id = hash_field(&value, "study_id")?;
    let n = integer(&value, "evaluations_completed")?;
    let op_id = ledger.artifact_output_seal(&hash)?.ok_or_else(|| invalid("unsealed design receipt"))?;
    let op = ledger.op(op_id)?.ok_or_else(|| invalid("missing design producer"))?;
    if value.str_field("schema") != Some(STUDY_RUN_RECEIPT_SCHEMA)
        || hash_field(&value, "model")? != id || n > fs_project::uncertainty::MAX_SAMPLES
        || op.session.as_deref() != Some(id.as_bytes().as_slice())
        || op.ir != ir(driver, id, n) || op.outcome.as_deref() != Some("ok")
        || !ledger.edge_exists(op_id, &hash, EdgeRole::Out)?
    {
        return Err(invalid("design receipt differs from its completed producing operation"));
    }
    for (key, kind, role) in [
        ("model", model::MANIFEST_KIND, EdgeRole::In),
        ("observations", observations, EdgeRole::Out),
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
    let driver = study_driver(study, policy);
    let observations = observations_kind(driver);
    let control = &study.parameters()[0];
    let limit = model.bound.threshold_k();
    let mut rows = Vec::new();
    let used_before = if let Some(old) = prior {
        if hash_field(&old.value, "model")? != model.identity()
            || old.value.str_field("driver") != Some(driver) {
            return Err(invalid("retained design model or objective changed"));
        }
        let observations = parse(&linked(ledger, &old.value, "observations", observations)?)?;
        if observations.str_field("driver") != Some(driver) {
            return Err(invalid("retained observations belong to another study producer"));
        }
        let samples = read_rows_value(observations.get("evaluations")
            .ok_or_else(|| invalid("missing retained design evaluations"))?)?;
        if samples.len() != integer(&old.value, "evaluations_completed")? || samples.len() > study.samples() {
            return Err(invalid("retained design evaluation count differs"));
        }
        for sample in samples {
            let decision = {
                let candidates = cases::completed(&rows, cases::count(study));
                next(control, policy, limit, &candidates)
            };
            let Next::Evaluate(expected) = decision else {
                return Err(invalid("retained design continues beyond a terminal search"));
            };
            let expected = cases::point(study, expected, rows.len());
            if sample.parameters.len() != expected.len()
                || sample.parameters.iter().zip(&expected).any(|(a, b)| a.to_bits() != b.to_bits())
            {
                return Err(invalid("retained design coordinate or operating case differs from the deterministic search"));
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
        let decision = {
            let candidates = cases::completed(&rows, cases::count(study));
            next(control, policy, limit, &candidates)
        };
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
        let Next::Evaluate(coordinate) = decision else { return Err(invalid("missing next design coordinate")); };
        let remaining = (study.wall_seconds() - used_before - started.elapsed().as_secs_f64()).max(0.0);
        let values = cases::point(study, coordinate, rows.len());
        match model.sample(ledger, gate, &values, remaining) {
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
        DesignPolicy { control: DesignControl::DeliveredPower { tolerance_w: 0.01 }, temperature_tolerance_k: 0.01 }
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
                    let chosen = selected(&rows, policy(), 305.0).unwrap();
                    assert!(chosen.temperature_k <= 305.0);
                    assert!((chosen.coordinate() - 2.5).abs() <= 0.01);
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

    fn fan_control() -> UniformParameter {
        UniformParameter { name: "fan-speed".into(), target: Target::FanSpeedRatio,
            entity: "bank".into(), low: 0.5, high: 1.5 }
    }

    fn fan_policy() -> DesignPolicy {
        DesignPolicy { control: DesignControl::FanSpeedRatio { tolerance: 0.001 }, temperature_tolerance_k: 0.01 }
    }

    #[test]
    fn g1_fan_bracket_minimizes_speed_from_a_measured_feasible_baseline() {
        let mut rows = Vec::new();
        assert_eq!(next(&fan_control(), fan_policy(), 305.0, &rows), Next::Evaluate(1.5));
        for _ in 0..32 {
            match next(&fan_control(), fan_policy(), 305.0, &rows) {
                Next::Evaluate(speed) => {
                    // Nonlinear cooling with an independently known threshold.
                    rows.push(trial(speed, 300.0 + 4.0 / (speed * speed)));
                    let chosen = selected(&rows, fan_policy(), 305.0).unwrap();
                    assert!(chosen.temperature_k <= 305.0);
                },
                Next::Complete(reason) => {
                    assert_eq!(reason, "target-bracketed");
                    let chosen = selected(&rows, fan_policy(), 305.0).unwrap();
                    assert!((chosen.coordinate() - 0.8_f64.sqrt()).abs() <= 0.001);
                    assert!(chosen.coordinate() >= 0.8_f64.sqrt());
                    let (passing, failing) = bracket(&rows, fan_policy(), 305.0).unwrap();
                    assert!(passing.coordinate() - failing.coordinate() <= 0.001);
                    assert!(failing.temperature_k > 305.0);
                    return;
                },
                Next::Refuse(reason) => panic!("{reason}"),
            }
        }
        panic!("nonlinear cooling response failed to bracket");
    }

    #[test]
    fn g0_fan_boundary_and_reversed_monotonicity_outcomes_are_explicit() {
        assert_eq!(next(&fan_control(), fan_policy(), 305.0, &[trial(1.5, 306.0)]),
            Next::Refuse("maximum-fan-speed-infeasible"));
        assert_eq!(next(&fan_control(), fan_policy(), 305.0, &[trial(1.5, 300.0)]),
            Next::Evaluate(0.5));
        let feasible = [trial(1.5, 300.0), trial(0.5, 305.0)];
        assert_eq!(next(&fan_control(), fan_policy(), 305.0, &feasible),
            Next::Complete("minimum-feasible"));
        assert_eq!(selected(&feasible, fan_policy(), 305.0).unwrap().coordinate(), 0.5);
        assert_eq!(next(&fan_control(), fan_policy(), 305.0,
            &[trial(1.5, 304.0), trial(0.5, 303.0)]), Next::Refuse("nonmonotone-response"));
    }

    #[test]
    fn g0_fan_history_uses_speed_units_and_a_distinct_producer() {
        let row = trial(0.8, 304.0);
        let fan_json = J::parse(&row.json(fan_policy())).unwrap();
        assert_eq!(fan_json.f64_field("fan_speed_ratio"), Some(0.8));
        assert!(fan_json.get("power_w").is_none());
        assert_eq!(driver(fan_policy()), FAN_DRIVER);
        assert_eq!(observations_kind(FAN_DRIVER), FAN_OBSERVATIONS);
        assert_ne!(ir(FAN_DRIVER, hash_bytes(b"model"), 1), ir(DRIVER, hash_bytes(b"model"), 1));
        assert!(is_driver(Some(FAN_DRIVER)));
        assert!(is_driver(Some(DRIVER)));
        assert!(!is_driver(Some("other-design")));
    }
}
