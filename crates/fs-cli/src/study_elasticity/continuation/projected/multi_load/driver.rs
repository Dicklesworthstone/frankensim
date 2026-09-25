//! Native orchestration of the existing simultaneous independent-load optimizer.
use super::*;

fn stopped(gate: &CancelGate, start: Instant, consumed: f64, spec: &ElasticitySpec) -> Option<&'static str> {
    stop_status(gate.is_requested(), consumed + start.elapsed().as_secs_f64(), spec.wall_s)
}

fn retain(spec: &ElasticitySpec, ledger: &Ledger, owner: &MultiLoadProjectedOptimizer,
    report: &OptimizeReport, status: &'static str, wall: f64, predecessor: Option<ContentHash>,
    evidence: &mut Evidence) -> Result<Outcome>
{
    let retained = evidence.projected.as_mut().expect("constrained history");
    if restoration::updates(&retained.baseline, &retained.accepted, &retained.policy)
        != owner.restoration_updates()
        || (owner.is_restoring_stress() && status == "completed")
    { return Err(malformed("restoration progress cannot be published as feasible completion")); }
    retained.family.as_mut()
        .expect("independent load history").capture(owner);
    persist(spec, ledger, owner.geometry(), report, status, wall, predecessor, evidence)
}

pub(in super::super) fn drive(spec: &ElasticitySpec, ledger: &Ledger, cap: Option<usize>,
    gate: &CancelGate, prior: Option<&Loaded>) -> Result<Outcome>
{
    drive_observed(spec, ledger, cap, gate, prior, |_| {})
}

pub(super) fn drive_observed(spec: &ElasticitySpec, ledger: &Ledger, cap: Option<usize>,
    gate: &CancelGate, prior: Option<&Loaded>, mut observe: impl FnMut(MultiLoadProjectedStage)) -> Result<Outcome>
{
    if ledger.in_transaction() { return Err(malformed("multi-load study requires its own ledger transaction")); }
    let policy = stress_controls(spec)?;
    let family = policy.family.as_ref().ok_or_else(|| malformed("missing independent load family"))?;
    let cases = family.cases(spec)?;
    let start = Instant::now();
    let mut evidence = Evidence { producer: producer_identity()?, updates: 0, legacy_replayed: 0,
        projected: None, volume: None };
    let mut report = OptimizeReport::default();
    let mut consumed = 0.0;
    let mut predecessor = prior.map(|old| old.hash);
    let mut last = None;
    let mut owner = if let Some(old) = prior {
        if old.value.str_field("study_id") != Some(spec.id.to_hex().as_str()) {
            return Err(malformed("retained multi-load study identity changed"));
        }
        let binding = old.value.get("continuation").ok_or_else(|| malformed("missing continuation binding"))?;
        if integer(binding, "version")? != 1
            || binding.str_field("producer") != Some(evidence.producer.to_hex().as_str())
        { return Err(malformed("multi-load resume requires the identical executable")); }
        let design = document(&linked(ledger, &old.value, "design", "study-design")?)?;
        let iterations = document(&linked(ledger, &old.value, "iterations", "study-iterations")?)?;
        let (phi, decoded) = decode(spec, &old.value, &design, &iterations)?;
        report = decoded;
        let mut retained = ConstraintEvidence::read(binding.get("constraints")
            .ok_or_else(|| malformed("missing multi-load constraints"))?, &report,
            spec.projected.as_ref().expect("admitted stress policy"))?;
        if retained.current().snapshot != snapshot(&phi) { return Err(malformed("multi-load field and stress differ")); }
        let expected_repairs = restoration::updates(&retained.baseline, &retained.accepted, policy);
        if old.value.str_field("status") == Some("completed") && !restoration::feasible(retained.current(), policy) {
            return Err(malformed("completed multi-load receipt is not stress feasible"));
        }
        let history = retained.family.as_mut().ok_or_else(|| malformed("missing multi-load checkpoint"))?;
        if cases_json(&history.cases) != cases_json(&cases) {
            return Err(malformed("retained primary or additional load differs from the native source"));
        }
        let status = match old.value.str_field("status") {
            Some("running") => "running", Some("completed") => "completed",
            Some("cancelled") => "cancelled", Some("budget-exhausted") => "budget-exhausted",
            Some("no-feasible-descent") => "no-feasible-descent",
            _ => return Err(malformed("unknown multi-load study terminal")),
        };
        consumed = old.value.f64_field("consumed_wall_s").filter(|value| value.is_finite() && *value >= 0.0)
            .ok_or_else(|| malformed("invalid retained multi-load wall charge"))?;
        last = Some(Outcome { pointer: format!("study-{}", old.hash.to_hex()), receipt: old.bytes.clone(), status });
        // These are sealed, read-only terminals. Do not rerun physics or consume
        // recovery work merely to return an already completed/exhausted receipt.
        if status == "completed" || status == "no-feasible-descent"
            || (status == "budget-exhausted"
                && (family.max_solves - history.solves < cases.len() || report.rows.len() == spec.steps))
        { return last.ok_or_else(|| malformed("missing terminal receipt")); }
        if let Some(status) = stopped(gate, start, consumed, spec) {
            return Err(constraints_stop(status, last.as_ref()));
        }
        let prepared = regions::prepare(spec, &policy.regions, |_| {
            match stopped(gate, start, consumed, spec) {
                Some(status) => ControlFlow::Break(status), None => ControlFlow::Continue(()),
            }
        })?;
        let prepared = match prepared {
            ControlFlow::Continue(prepared) => prepared,
            ControlFlow::Break(status) => return Err(constraints_stop(status, last.as_ref())),
        };
        let mut remaining = family.max_recovery_solves - history.recovery_solves;
        if remaining < 2 * cases.len() {
            return Err(retained_error(Failure { code: "cli-study-multi-load-recovery-budget",
                message: "recovery needs two complete independent load families; declared lifetime recovery allowance is exhausted".into(),
                exit: exit::BUDGET }, last.as_ref()));
        }
        // The library restores its original baseline, current solutions, global
        // ordinal, multiplier and spent study solves. No prefix update is run.
        // This owner currently exposes synchronous baseline/recovery solves;
        // native reports state that boundary rather than claiming preemption.
        let replay = MultiLoadProjectedOptimizer::restore_checkpoint(&history.checkpoint, &mut remaining)
            .map_err(|error| malformed(&error.to_string()));
        history.recovery_solves = family.max_recovery_solves - remaining;
        let checked = replay.and_then(|restored| {
            if restored.load_cases().iter().any(|case| case.edge() != DesignBoxEdge::Right) {
                return Err(malformed("checkpoint changed a native load edge"));
            }
            bind(&restored, spec, policy, family, &prepared.fixed_nodes)?;
            let stress = restored.current_stress().ok_or_else(|| malformed("recovered family has no stress"))?;
            let restored_cases = case_states(stress);
            let expected_cases = history.accepted.last().unwrap_or(&history.baseline);
            if restored.next_iteration() != report.rows.len() || restored.solves_started() != history.solves
                || restored.restoration_updates() != expected_repairs
                || restored.geometry().nodes().iter().zip(phi.nodes()).any(|(a, b)| a.to_bits() != b.to_bits())
                || restored.geometry().n() != phi.n()
                || states_json(&restored_cases) != states_json(expected_cases)
                || states_json(&case_states(restored.baseline_stress().ok_or_else(|| malformed("missing original baseline stress"))?))
                    != states_json(&history.baseline)
                || restored.search_multiplier().to_bits() != report.ell.last().copied().unwrap_or(spec.ell0).to_bits()
            { return Err(malformed("multi-load checkpoint disagrees with the retained native history")); }
            Ok(restored)
        });
        evidence.projected = Some(retained);
        match checked {
            Ok(restored) => restored,
            Err(error) => {
                // Charge attempted recovery even on a later replay refusal.
                // Re-retain the OLD geometry/history, never the failed result.
                // This running receipt is a recovery pointer, not a success
                // returned to the user; the command still returns the error.
                let charged = persist(spec, ledger, &phi, &report, "running",
                    consumed + start.elapsed().as_secs_f64(), predecessor, &evidence)
                    .map_err(|write| retained_error(write, last.as_ref()))?;
                return Err(retained_error(error, Some(&charged)));
            }
        }
    } else {
        let prepared = regions::prepare(spec, &policy.regions, |_| {
            match stopped(gate, start, 0.0, spec) {
                Some(status) => ControlFlow::Break(status), None => ControlFlow::Continue(()),
            }
        })?;
        let prepared = match prepared {
            ControlFlow::Continue(prepared) => prepared,
            ControlFlow::Break(status) => return Err(constraints_stop(status, None)),
        };
        if let Some(status) = stopped(gate, start, 0.0, spec) { return Err(constraints_stop(status, None)); }
        let owner = MultiLoadProjectedOptimizer::new(prepared.geometry, &cases, settings(spec, spec.steps),
            family.aggregate, prepared.fixed_nodes, policy.area, family.controls(policy))
            .and_then(|owner| match family.restoration_reduction {
                Some(reduction) => owner.with_stress_restoration(policy.stress, reduction),
                None => owner.with_sampled_stress_limit(policy.stress),
            })
            .map_err(|error| malformed(&error.to_string()))?;
        let baseline = case_states(owner.baseline_stress().ok_or_else(|| malformed("missing baseline family stress"))?);
        let combined = summary(&baseline, &cases, family.aggregate)?;
        if combined.compliance.to_bits() != owner.baseline().objective.to_bits() {
            return Err(malformed("native family summary disagrees with the numerical objective"));
        }
        evidence.projected = Some(ConstraintEvidence {
            policy: policy.clone(), baseline: combined, accepted: Vec::new(), attempts: Vec::new(), refusals: Vec::new(),
            family: Some(History { cases, baseline, accepted: Vec::new(), checkpoint: Vec::new(),
                solves: owner.solves_started(), recovery_solves: 0 }),
        });
        owner
    };
    let target = spec.steps.min(report.rows.len().saturating_add(cap.unwrap_or(spec.steps - report.rows.len())));
    // Publish admitted recovery work before starting another candidate, charging
    // it once along this receipt chain. Initial expired admission publishes none.
    let stop = stopped(gate, start, consumed, spec);
    if last.is_none() && stop.is_some() { return Err(constraints_stop(stop.unwrap(), None)); }
    let accepted = retain(spec, ledger, &owner, &report, stop.unwrap_or("running"),
        consumed + start.elapsed().as_secs_f64(), predecessor, &mut evidence)
        .map_err(|error| retained_error(error, last.as_ref()))?;
    if stop.is_some() { return Ok(accepted); }
    predecessor = accepted.pointer.strip_prefix("study-").and_then(ContentHash::from_hex);
    last = Some(accepted);
    loop {
        let status = stopped(gate, start, consumed, spec).or_else(||
            (owner.next_iteration() == target).then_some(
                if owner.next_iteration() == spec.steps && !owner.is_restoring_stress() {
                    "completed"
                } else { "budget-exhausted" }));
        if let Some(status) = status {
            return retain(spec, ledger, &owner, &report, status, consumed + start.elapsed().as_secs_f64(),
                predecessor, &mut evidence).map_err(|error| retained_error(error, last.as_ref()));
        }
        let progress = owner.advance_one_polling(policy.search.poll_iters, |stage| {
            observe(stage);
            match stopped(gate, start, consumed, spec) {
                Some(status) => ControlFlow::Break(status), None => ControlFlow::Continue(()),
            }
        }).map_err(|error| retained_error(malformed(&error.to_string()), last.as_ref()))?;
        let progress = match progress {
            ControlFlow::Continue(progress) => progress,
            ControlFlow::Break(status) => return retain(spec, ledger, &owner, &report, status,
                consumed + start.elapsed().as_secs_f64(), predecessor, &mut evidence)
                .map_err(|error| retained_error(error, last.as_ref())),
        };
        let retained = evidence.projected.as_mut().expect("admitted multi-load history");
        match progress {
            MultiLoadProjectedProgress::Accepted(step) => {
                let measurements = case_states(owner.current_stress().ok_or_else(|| malformed("accepted family has no stress"))?);
                let current = summary(&measurements, owner.load_cases(), family.aggregate)?;
                let repairing = restoration::transition(retained.current(), &current, policy)
                    .map_err(|error| retained_error(error, last.as_ref()))?;
                if step.restoration != repairing || step.iteration != report.rows.len() || step.state.snapshot != current.snapshot
                    || step.state.objective.to_bits() != current.compliance.to_bits()
                    || step.state.volume.to_bits() != current.volume.to_bits()
                { return Err(retained_error(malformed("accepted multi-load metrics disagree"), last.as_ref())); }
                let ell = owner.search_multiplier();
                let drift = step.proposal_audit.interface_drift_h;
                let pad = step.proposal_load_pad_nodes;
                report.rows.push(format!("{{\"iter\":{},\"compliance\":{:.17e},\"volume\":{:.17e},\"ell\":{ell:.17e},\"drift_h\":{drift:.17e},\"load_pad_nodes\":{pad},\"snapshot\":\"{:#018x}\"}}",
                    step.iteration, current.compliance, current.volume, current.snapshot));
                report.compliance.push(current.compliance); report.volume.push(current.volume);
                report.ell.push(ell); report.snapshots.push(current.snapshot); report.load_pad_nodes.push(pad);
                retained.accepted.push(current); retained.attempts.push(step.attempts.len()); retained.refusals.clear();
                retained.family.as_mut().expect("family history").accepted.push(measurements);
                evidence.updates += 1;
            }
            other => {
                let (status, attempts) = match other {
                    MultiLoadProjectedProgress::NoDescent(attempts) => ("no-feasible-descent", attempts),
                    MultiLoadProjectedProgress::SolveBudget(attempts) => ("budget-exhausted", attempts),
                    _ => return Err(malformed("multi-load optimizer completed before its native target")),
                };
                retained.refusals = attempts.into_iter().map(|attempt| attempt.refusal
                    .unwrap_or_else(|| "candidate family was not admitted".into())).collect();
                return retain(spec, ledger, &owner, &report, status, consumed + start.elapsed().as_secs_f64(),
                    predecessor, &mut evidence).map_err(|error| retained_error(error, last.as_ref()));
            }
        }
        let accepted = retain(spec, ledger, &owner, &report, "running", consumed + start.elapsed().as_secs_f64(),
            predecessor, &mut evidence).map_err(|error| retained_error(error, last.as_ref()))?;
        predecessor = accepted.pointer.strip_prefix("study-").and_then(ContentHash::from_hex);
        last = Some(accepted);
    }
}
