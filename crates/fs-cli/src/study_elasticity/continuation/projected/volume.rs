//! Hard-area-only native elasticity studies over the existing projected owner.
//! No invented stress limit, new optimizer, or alternate checkpoint format.
use super::*;
use fs_topols::EvaluatedFinalState;
use fs_topols::projected::ProjectedSetupStage;

#[path = "volume/mesh.rs"]
mod mesh;

#[derive(Debug, Clone)]
pub(crate) struct Controls {
    area: VolumeProjectionSettings,
    search: ProjectedSettings,
    regions: Vec<DesignRegion>,
    resolution: Option<mesh::ResolutionPolicy>,
}

impl Controls {
    pub(crate) fn parse(fields: &[Node], target: f64) -> Result<Self> {
        let real = |key| super::super::super::number(field(fields, key)?, key);
        let count = |key| integer_node(field(fields, key)?, key);
        let policy = Self {
            area: VolumeProjectionSettings {
                target, tolerance: real("area-tolerance-m2")?,
                max_shift: real("max-projection-shift")?,
                max_evaluations: count("max-area-evaluations")?,
            },
            search: ProjectedSettings {
                max_candidates: count("max-candidates")?, contraction: real("contraction")?,
                min_relative_improvement: real("min-relative-improvement")?,
                poll_iters: count("cg-poll-iters")?,
            },
            regions: regions::parse_regions(fields)?,
            resolution: mesh::parse(fields)?,
        };
        let a = policy.area;
        let s = policy.search;
        if !(a.tolerance > 0.0 && a.tolerance <= 0.01 && a.max_shift > 0.0)
            || !(3..=64).contains(&a.max_evaluations)
            || !(1..=16).contains(&s.max_candidates)
            || !(s.contraction > 0.0 && s.contraction < 1.0)
            || !(0.0..1.0).contains(&s.min_relative_improvement)
            || !(1..=60_000).contains(&s.poll_iters)
        { return Err(malformed("projected-volume requires area tolerance in (0,0.01], positive shift, 3..=64 area evaluations, 1..=16 candidates, contraction (0,1), relative decrease [0,1), and 1..=60000 CG polling iterations")); }
        // The parent canonical comparison refuses unknown, duplicate, reordered
        // and stress/load-family fields instead of silently ignoring them.
        Ok(policy)
    }

    pub(crate) fn canonical(&self, out: &mut String) {
        let _ = writeln!(out, "    :constraint-mode projected-volume");
        let _ = writeln!(out, "    :area-tolerance-m2 {}", canonical_float(self.area.tolerance));
        let _ = writeln!(out, "    :max-projection-shift {}", canonical_float(self.area.max_shift));
        let _ = writeln!(out, "    :max-area-evaluations {}", self.area.max_evaluations);
        let _ = writeln!(out, "    :max-candidates {}", self.search.max_candidates);
        let _ = writeln!(out, "    :contraction {}", canonical_float(self.search.contraction));
        let _ = writeln!(out, "    :min-relative-improvement {}", canonical_float(self.search.min_relative_improvement));
        if self.regions.is_empty() && self.resolution.is_none() {
            let _ = writeln!(out, "    :cg-poll-iters {})", self.search.poll_iters);
        } else {
            let _ = writeln!(out, "    :cg-poll-iters {}", self.search.poll_iters);
            if let Some(policy) = self.resolution { mesh::canonical(policy, out); }
            if self.regions.is_empty() { let _ = writeln!(out, "  )"); }
            else { regions::canonical(&self.regions, out); }
        }
    }
}

/// The measured projection of an independently evaluated design, not a stress
/// evaluation and not a replacement PDE report.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Measured {
    pub(crate) compliance: f64,
    pub(crate) volume: f64,
    snapshot: u64,
}
impl From<EvaluatedFinalState> for Measured {
    fn from(value: EvaluatedFinalState) -> Self {
        Self { compliance: value.compliance, volume: value.volume, snapshot: value.snapshot }
    }
}
impl Measured {
    fn json(self) -> String {
        format!("{{\"compliance_j\":{:.17e},\"area_m2\":{:.17e},\"snapshot\":\"{:#018x}\"}}",
            self.compliance, self.volume, self.snapshot)
    }
    fn read(value: &JsonValue, policy: &Controls) -> Result<Self> {
        let measured = Self {
            compliance: number(value, "compliance_j")?.0,
            volume: number(value, "area_m2")?.0,
            snapshot: hexadecimal(value.str_field("snapshot")
                .ok_or_else(|| malformed("missing volume-only state snapshot"))?, true)?,
        };
        if measured.compliance < 0.0 || measured.volume <= 0.0
            || (measured.volume - policy.area.target).abs() > policy.area.tolerance {
            return Err(malformed("retained projected-volume state violates material feasibility"));
        }
        Ok(measured)
    }
    fn same(self, other: Self) -> bool {
        self.snapshot == other.snapshot && self.compliance.to_bits() == other.compliance.to_bits()
            && self.volume.to_bits() == other.volume.to_bits()
    }
}

/// History is bounded by the native 32-update admission. Baseline projection is
/// NOT counted as compliance improvement from the overfilled input geometry.
pub(crate) struct VolumeEvidence {
    policy: Controls,
    baseline: Measured,
    accepted: Vec<Measured>,
    attempts: Vec<usize>,
    refusals: Vec<String>,
    mesh: Option<mesh::LastCheck>,
}
impl VolumeEvidence {
    pub(crate) fn current(&self) -> Measured {
        self.accepted.last().copied().unwrap_or(self.baseline)
    }
    pub(crate) fn json(&self) -> String {
        let reduction = if self.baseline.compliance > 0.0 {
            format!("{:.17e}", 1.0 - self.current().compliance / self.baseline.compliance)
        } else { "null".into() };
        format!(concat!("{{\"mode\":\"projected-volume-v1\",\"baseline_scope\":\"same-material-same-load\",",
            "\"area_target_m2\":{:.17e},\"area_tolerance_m2\":{:.17e},",
            "\"min_relative_improvement\":{:.17e},\"baseline\":{},\"accepted\":[{}],",
            "\"candidate_counts\":[{}],\"terminal_refusals\":[{}],",
            "\"relative_reduction\":{},\"area_constraint_satisfied\":true,",
            "\"stress_evaluation\":\"not-requested\"{}}}"),
            self.policy.area.target, self.policy.area.tolerance, self.policy.search.min_relative_improvement,
            self.baseline.json(), self.accepted.iter().map(|state| state.json()).collect::<Vec<_>>().join(","),
            self.attempts.iter().map(usize::to_string).collect::<Vec<_>>().join(","),
            self.refusals.iter().map(|s| quoted(s)).collect::<Vec<_>>().join(","), reduction,
            regions::json_field(&self.policy.regions) + &mesh::field(self.policy.resolution, self.mesh.as_ref()))
    }
    pub(crate) fn html(&self) -> String {
        format!("<p>Hard material area: {:.8e} m² ± {:.8e} m². Independently solved same-material, same-load baseline: {:.8e} J; current compliance: {:.8e} J. Baseline area projection is feasibility preparation, not an optimization improvement. No stress limit or stress evaluation was requested. Prescribed material/void regions: {}. The area is numerical cut quadrature; iteration completion is not convergence or optimality.</p>",
            self.policy.area.target, self.policy.area.tolerance, self.baseline.compliance,
            self.current().compliance, self.policy.regions.len())
            + &mesh::html(self.policy.resolution, self.mesh.as_ref())
    }
    fn read(value: &JsonValue, report: &OptimizeReport, policy: &Controls) -> Result<Self> {
        if value.str_field("mode") != Some("projected-volume-v1")
            || value.str_field("baseline_scope") != Some("same-material-same-load")
            || value.str_field("stress_evaluation") != Some("not-requested")
            || value.get("area_constraint_satisfied") != Some(&JsonValue::Bool(true)) {
            return Err(malformed("missing volume-only constraint identity"));
        }
        for (key, expected) in [("area_target_m2", policy.area.target),
            ("area_tolerance_m2", policy.area.tolerance),
            ("min_relative_improvement", policy.search.min_relative_improvement)] {
            if number(value, key)?.0.to_bits() != expected.to_bits() {
                return Err(malformed("retained volume-only policy differs from the canonical study"));
            }
        }
        regions::check_retained(value, &policy.regions)?;
        let array = |key| value.get(key).and_then(JsonValue::as_array)
            .ok_or_else(|| malformed("missing volume-only history array"));
        let baseline = Measured::read(value.get("baseline")
            .ok_or_else(|| malformed("missing same-material baseline"))?, policy)?;
        let rows = array("accepted")?;
        let attempts = array("candidate_counts")?;
        let refusals = array("terminal_refusals")?;
        if rows.len() != report.rows.len() || attempts.len() != rows.len()
            || refusals.len() > policy.search.max_candidates {
            return Err(malformed("projected-volume history lengths disagree"));
        }
        let mut accepted = Vec::with_capacity(rows.len());
        let mut counts = Vec::with_capacity(rows.len());
        let mut previous = baseline;
        for (index, (row, count)) in rows.iter().zip(attempts).enumerate() {
            let state = Measured::read(row, policy)?;
            let count = count.number_raw().and_then(|raw| raw.parse::<usize>().ok())
                .ok_or_else(|| malformed("invalid volume-only candidate count"))?;
            if !(1..=policy.search.max_candidates).contains(&count)
                || state.snapshot != report.snapshots[index]
                || state.compliance.to_bits() != report.compliance[index].to_bits()
                || state.volume.to_bits() != report.volume[index].to_bits()
                || !(state.compliance < previous.compliance * (1.0 - policy.search.min_relative_improvement)) {
                return Err(malformed("projected-volume history is not a decreasing same-area trajectory"));
            }
            accepted.push(state); counts.push(count); previous = state;
        }
        let refusals = refusals.iter().map(|v| v.as_str().map(str::to_string)
            .ok_or_else(|| malformed("invalid projected-volume refusal"))).collect::<Result<Vec<_>>>()?;
        let current = accepted.last().copied().unwrap_or(baseline);
        let previous = accepted.get(accepted.len().saturating_sub(2)).copied()
            .filter(|_| accepted.len() >= 2).unwrap_or(baseline);
        let mesh = mesh::read(value, policy, previous, current, accepted.len())?;
        let evidence = Self { policy: policy.clone(), baseline, accepted, attempts: counts, refusals, mesh };
        // Recompute rather than trusting the retained headline improvement.
        let expected = document(evidence.json().as_bytes())?;
        if value.get("relative_reduction") != expected.get("relative_reduction") {
            return Err(malformed("retained improvement differs from its same-area baseline"));
        }
        Ok(evidence)
    }
}

#[derive(Debug, Clone, Copy)]
enum VolumeStage {
    Regions(DesignRegionStage),
    Setup(ProjectedSetupStage),
    Update(ProjectedStage),
    Mesh(mesh::MeshCheckStage),
}
fn stopped(gate: &CancelGate, start: Instant, consumed: f64, spec: &ElasticitySpec) -> Option<&'static str> {
    stop_status(gate.is_requested(), consumed + start.elapsed().as_secs_f64(), spec.wall_s)
}

pub(in super::super) fn drive(spec: &ElasticitySpec, ledger: &Ledger, cap: Option<usize>,
    gate: &CancelGate, prior: Option<&Loaded>) -> Result<Outcome> {
    drive_observed(spec, ledger, cap, gate, prior, |_| {})
}

fn drive_observed(spec: &ElasticitySpec, ledger: &Ledger, cap: Option<usize>,
    gate: &CancelGate, prior: Option<&Loaded>, mut observe: impl FnMut(VolumeStage)) -> Result<Outcome> {
    if ledger.in_transaction() { return Err(malformed("projected-volume requires its own ledger transaction")); }
    let Some(ProjectedControls::Volume(policy)) = spec.projected.as_ref()
        else { return Err(malformed("missing explicit projected-volume policy")); };
    if let Some(policy) = policy.resolution {
        policy.validate(settings(spec, spec.steps).level).map_err(|e| malformed(&e.to_string()))?;
    }
    let start = Instant::now();
    let mut evidence = Evidence { producer: producer_identity()?, updates: 0, legacy_replayed: 0,
        projected: None, volume: None };
    let mut report = OptimizeReport::default();
    let mut consumed = 0.0;
    let mut predecessor = prior.map(|old| old.hash);
    let mut last = None;
    let mut state = if let Some(old) = prior {
        if old.value.str_field("study_id") != Some(spec.id.to_hex().as_str()) {
            return Err(malformed("retained projected-volume study identity changed"));
        }
        let binding = old.value.get("continuation").ok_or_else(|| malformed("missing continuation binding"))?;
        if integer(binding, "version")? != 1
            || binding.str_field("producer") != Some(evidence.producer.to_hex().as_str()) {
            return Err(malformed("projected-volume continuation requires the identical executable"));
        }
        let design = document(&linked(ledger, &old.value, "design", "study-design")?)?;
        let rows = document(&linked(ledger, &old.value, "iterations", "study-iterations")?)?;
        let (phi, decoded) = decode(spec, &old.value, &design, &rows)?;
        report = decoded;
        let retained = VolumeEvidence::read(binding.get("constraints")
            .ok_or_else(|| malformed("missing retained area constraints"))?, &report, policy)?;
        if retained.current().snapshot != snapshot(&phi) {
            return Err(malformed("projected-volume field and evaluated state disagree"));
        }
        if retained.mesh.as_ref().is_some_and(|check|
            check.baseline.rungs[0].level != settings(spec, spec.steps).level) {
            return Err(malformed("retained mesh-check levels differ from the study"));
        }
        let status = match old.value.str_field("status") {
            Some("running") => "running", Some("completed") => "completed",
            Some("cancelled") => "cancelled", Some("budget-exhausted") => "budget-exhausted",
            Some("no-feasible-descent") => "no-feasible-descent",
            Some("mesh-unresolved") if retained.mesh.as_ref().is_some_and(|c| c.outcome == "baseline-unresolved") => "mesh-unresolved",
            _ => return Err(malformed("unknown projected-volume terminal")),
        };
        consumed = old.value.f64_field("consumed_wall_s").filter(|v| v.is_finite() && *v >= 0.0)
            .ok_or_else(|| malformed("invalid retained projected-volume wall charge"))?;
        evidence.volume = Some(retained);
        last = Some(Outcome { pointer: format!("study-{}", old.hash.to_hex()), receipt: old.bytes.clone(), status });
        if matches!(status, "completed" | "no-feasible-descent" | "mesh-unresolved") {
            return last.ok_or_else(|| malformed("terminal volume-only state has no receipt"));
        }
        if let Some(status) = stopped(gate, start, consumed, spec) {
            return Err(constraints_stop(status, last.as_ref()));
        }
        let prepared = regions::prepare(spec, &policy.regions, |stage| {
            observe(VolumeStage::Regions(stage));
            match stopped(gate, start, consumed, spec) {
                Some(status) => ControlFlow::Break(status), None => ControlFlow::Continue(()),
            }
        })?;
        let prepared = match prepared {
            ControlFlow::Continue(prepared) => prepared,
            ControlFlow::Break(status) => return persist(spec, ledger, &phi, &report, status,
                consumed + start.elapsed().as_secs_f64(), predecessor, &evidence),
        };
        let checkpoint = OptimizeCheckpoint::restore(phi, fixture(spec), settings(spec, spec.steps),
            report.rows.len(), report.ell.last().copied().unwrap_or(spec.ell0))
            .map_err(|error| malformed(&error.to_string()))?;
        let restored = ProjectedOptimizer::from_checkpoint_controlled(&checkpoint, prepared.fixed_nodes,
            policy.area, policy.search, |stage| {
                observe(VolumeStage::Setup(stage));
                match stopped(gate, start, consumed, spec) {
                    Some(status) => ControlFlow::Break(status), None => ControlFlow::Continue(()),
                }
            });
        match restored {
            Ok(ControlFlow::Continue(restored)) => {
                if !Measured::from(restored.current()).same(evidence.volume.as_ref().expect("retained area state").current()) {
                    return Err(malformed("independent projected-volume replay differs from the retained endpoint"));
                }
                restored
            }
            Ok(ControlFlow::Break(status)) => return persist(spec, ledger, checkpoint.geometry(), &report,
                status, consumed + start.elapsed().as_secs_f64(), predecessor, &evidence),
            Err(error) => {
                let charged = persist(spec, ledger, checkpoint.geometry(), &report, "running",
                    consumed + start.elapsed().as_secs_f64(), predecessor, &evidence)
                    .map_err(|error| retained_error(error, last.as_ref()))?;
                return Err(retained_error(malformed(&error.to_string()), Some(&charged)));
            }
        }
    } else {
        let prepared = regions::prepare(spec, &policy.regions, |stage| {
            observe(VolumeStage::Regions(stage));
            match stopped(gate, start, 0.0, spec) {
                Some(status) => ControlFlow::Break(status), None => ControlFlow::Continue(()),
            }
        })?;
        let prepared = match prepared {
            ControlFlow::Continue(prepared) => prepared,
            ControlFlow::Break(status) => return Err(constraints_stop(status, None)),
        };
        let state = ProjectedOptimizer::new_controlled(&prepared.geometry, fixture(spec), settings(spec, spec.steps),
            prepared.fixed_nodes, policy.area, policy.search, |stage| {
                observe(VolumeStage::Setup(stage));
                match stopped(gate, start, 0.0, spec) {
                    Some(status) => ControlFlow::Break(status), None => ControlFlow::Continue(()),
                }
            }).map_err(|error| malformed(&error.to_string()))?;
        let state = match state {
            ControlFlow::Continue(state) => state,
            ControlFlow::Break(status) => return Err(constraints_stop(status, None)),
        };
        evidence.volume = Some(VolumeEvidence { policy: policy.clone(), baseline: state.current().into(),
            accepted: Vec::new(), attempts: Vec::new(), refusals: Vec::new(), mesh: None });
        state
    };
    let target = spec.steps.min(report.rows.len().saturating_add(cap.unwrap_or(spec.steps - report.rows.len())));
    let stop = stopped(gate, start, consumed, spec);
    if last.is_none() && stop.is_some() { return Err(constraints_stop(stop.unwrap(), None)); }
    // Retain the feasible baseline and charge recovery before another update.
    let admitted = persist(spec, ledger, state.checkpoint().geometry(), &report, stop.unwrap_or("running"),
        consumed + start.elapsed().as_secs_f64(), predecessor, &evidence)
        .map_err(|error| retained_error(error, last.as_ref()))?;
    if stop.is_some() { return Ok(admitted); }
    predecessor = admitted.pointer.strip_prefix("study-").and_then(ContentHash::from_hex);
    last = Some(admitted);
    loop {
        let status = stopped(gate, start, consumed, spec).or_else(||
            (state.checkpoint().next_iteration() == target).then_some(
                if state.checkpoint().is_complete() { "completed" } else { "budget-exhausted" }));
        if let Some(status) = status {
            return persist(spec, ledger, state.checkpoint().geometry(), &report, status,
                consumed + start.elapsed().as_secs_f64(), predecessor, &evidence)
                .map_err(|error| retained_error(error, last.as_ref()));
        }
        let progress = if let Some(resolution) = policy.resolution {
            let checked = state.advance_one_resolution_controlled(resolution, |stage| {
                observe(VolumeStage::Mesh(stage));
                match stopped(gate, start, consumed, spec) {
                    Some(status) => ControlFlow::Break(status), None => ControlFlow::Continue(()),
                }
            }).map_err(|error| retained_error(malformed(&error.to_string()), last.as_ref()))?;
            match checked {
                ControlFlow::Break(status) => ControlFlow::Break(status),
                ControlFlow::Continue(checked) => {
                    let (progress, check) = mesh::LastCheck::capture(checked)?;
                    let retained = evidence.volume.as_mut().expect("admitted volume-only state");
                    if let Some(reason) = &check.reason { retained.refusals = vec![reason.clone()]; }
                    retained.mesh = Some(check);
                    ControlFlow::Continue(progress)
                }
            }
        } else {
            match state.advance_one_controlled(|stage| {
                observe(VolumeStage::Update(stage));
                match stopped(gate, start, consumed, spec) {
                    Some(status) => ControlFlow::Break(status), None => ControlFlow::Continue(()),
                }
            }).map_err(|error| retained_error(malformed(&error.to_string()), last.as_ref()))? {
                ControlFlow::Break(status) => ControlFlow::Break(status),
                ControlFlow::Continue(progress) => ControlFlow::Continue(Some(progress)),
            }
        };
        let progress = match progress {
            ControlFlow::Continue(Some(progress)) => progress,
            stopped => {
                let status = match stopped {
                    ControlFlow::Break(status) => status,
                    _ => "mesh-unresolved",
                };
                return persist(spec, ledger, state.checkpoint().geometry(), &report,
                    status, consumed + start.elapsed().as_secs_f64(), predecessor, &evidence)
                    .map_err(|error| retained_error(error, last.as_ref()));
            }
        };
        let retained = evidence.volume.as_mut().expect("admitted volume-only state");
        match progress {
            ProjectedProgress::Accepted(step) => {
                let current = Measured::from(state.current());
                if step.iteration != report.rows.len() || !Measured::from(step.state).same(current)
                    || !Measured::from(step.previous).same(retained.current()) {
                    return Err(retained_error(malformed("accepted projected-volume mechanics disagree"), last.as_ref()));
                }
                let drift = step.proposal.audits.first().ok_or_else(|| malformed("missing proposal audit"))?.interface_drift_h;
                let pad = step.proposal.load_pad_nodes.first().copied().ok_or_else(|| malformed("missing proposal load-pad count"))?;
                let ell = state.checkpoint().ell();
                report.rows.push(format!("{{\"iter\":{},\"compliance\":{:.17e},\"volume\":{:.17e},\"ell\":{ell:.17e},\"drift_h\":{drift:.17e},\"load_pad_nodes\":{pad},\"snapshot\":\"{:#018x}\"}}",
                    step.iteration, current.compliance, current.volume, current.snapshot));
                report.compliance.push(current.compliance); report.volume.push(current.volume);
                report.ell.push(ell); report.snapshots.push(current.snapshot); report.load_pad_nodes.push(pad);
                retained.accepted.push(current); retained.attempts.push(step.attempts.len()); retained.refusals.clear();
                evidence.updates += 1;
            }
            ProjectedProgress::NoDescent(attempts) => {
                retained.refusals = attempts.into_iter().map(|attempt| attempt.refusal
                    .unwrap_or_else(|| "candidate was not admitted".into())).collect();
                return persist(spec, ledger, state.checkpoint().geometry(), &report, "no-feasible-descent",
                    consumed + start.elapsed().as_secs_f64(), predecessor, &evidence)
                    .map_err(|error| retained_error(error, last.as_ref()));
            }
            ProjectedProgress::IterationLimit => return Err(malformed("projected-volume completed before target")),
        }
        let accepted = persist(spec, ledger, state.checkpoint().geometry(), &report, "running",
            consumed + start.elapsed().as_secs_f64(), predecessor, &evidence)
            .map_err(|error| retained_error(error, last.as_ref()))?;
        predecessor = accepted.pointer.strip_prefix("study-").and_then(ContentHash::from_hex);
        last = Some(accepted);
    }
}

#[cfg(test)]
#[path = "volume/tests.rs"]
mod tests;
