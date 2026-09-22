//! Hard-area and sampled-stress admission in the existing native study/ledger.
//! Numerical work belongs to fs-topols; this module only binds declarations,
//! retained accepted state, stopping and reports to that existing controller.
use super::*;
use std::fmt::Write as _;
use fs_topols::projected::{ProjectedOptimizer, ProjectedProgress, ProjectedSettings, ProjectedStage};
use fs_topols::projected_stress::ProjectedStressSetupStage;
use fs_topols::volume::VolumeProjectionSettings;
use fs_topols::{ProjectedStressOptimizer, SampledStressEvaluation, SampledStressLimit};
use fs_topols::design_regions::{DesignRegion, DesignRegionStage};

#[path = "projected/regions.rs"]
mod regions;

#[path = "projected/multi_load.rs"]
mod multi_load;
use multi_load::restoration;

pub(crate) const PROJECTED_SCOPE: &str = "2-D plane-strain CutFEM with numerical material-area equality and a declared deterministic sampled von Mises limit. Explicit restoration may retain overstressed designs only while reducing measured worst stress excess; these are not stress-feasible results. Compliance improvement starts at the first stress-feasible design under identical loads. Every accepted state is independently re-solved and durably retained. Candidate CG and stress-cell boundaries are cancellable; the load-family report states startup/recovery limitations. Assembly, area quadrature and ledger I/O remain indivisible. Iteration completion is not convergence. Drift and nucleation diagnostics describe proposals, not projected geometry. No physical validation, continuous stress/volume certificate, stress-adjoint/KKT/global optimum, 3-D result or guaranteed discretization-error bound is claimed.";

#[derive(Debug, Clone)]
pub(crate) struct Controls {
    area: VolumeProjectionSettings,
    search: ProjectedSettings,
    stress: SampledStressLimit,
    regions: Vec<DesignRegion>,
    family: Option<multi_load::LoadFamily>,
}

pub(crate) fn parse_controls(fields: &[Node], target: f64) -> Result<Option<Controls>> {
    let mode = fields.windows(2).find(|pair|
        matches!(&pair[0].kind, NodeKind::Keyword(key) if key == "constraint-mode"));
    let Some(mode) = mode else { return Ok(None) };
    if !matches!(&mode[1].kind, NodeKind::Symbol(value) if value == "projected-stress") {
        return Err(malformed("constraint-mode must be projected-stress"));
    }
    let real = |key| super::super::number(field(fields, key)?, key);
    let count = |key| integer_node(field(fields, key)?, key);
    let controls = Controls {
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
        stress: SampledStressLimit::new(real("sampled-stress-limit-pa")?, real("stress-tolerance-pa")?)
            .map_err(|error| malformed(&error.to_string()))?,
        regions: regions::parse_regions(fields)?,
        family: multi_load::parse(fields)?,
    };
    let a = controls.area;
    let s = controls.search;
    if !(a.tolerance > 0.0 && a.tolerance <= 0.01 && a.max_shift > 0.0)
        || !(3..=64).contains(&a.max_evaluations)
        || !(1..=16).contains(&s.max_candidates)
        || !(s.contraction > 0.0 && s.contraction < 1.0)
        || !(0.0..1.0).contains(&s.min_relative_improvement)
        || !(1..=60_000).contains(&s.poll_iters)
    {
        return Err(malformed("projected-stress requires area tolerance in (0,0.01], positive shift, 3..=64 area evaluations, 1..=16 candidates, contraction (0,1), relative decrease [0,1), and 1..=60000 CG polling iterations"));
    }
    Ok(Some(controls))
}

impl Controls {
    pub(crate) fn canonical(&self, out: &mut String) {
        let _ = writeln!(out, "    :constraint-mode projected-stress");
        for (key, value) in [("area-tolerance-m2", self.area.tolerance),
            ("max-projection-shift", self.area.max_shift)]
        { let _ = writeln!(out, "    :{key} {}", canonical_float(value)); }
        let _ = writeln!(out, "    :max-area-evaluations {}", self.area.max_evaluations);
        let _ = writeln!(out, "    :max-candidates {}", self.search.max_candidates);
        let _ = writeln!(out, "    :contraction {}", canonical_float(self.search.contraction));
        let _ = writeln!(out, "    :min-relative-improvement {}", canonical_float(self.search.min_relative_improvement));
        let _ = writeln!(out, "    :cg-poll-iters {}", self.search.poll_iters);
        let _ = writeln!(out, "    :sampled-stress-limit-pa {}", canonical_float(self.stress.max_von_mises));
        if self.regions.is_empty() && self.family.is_none() {
            let _ = writeln!(out, "    :stress-tolerance-pa {})", canonical_float(self.stress.absolute_tolerance));
        } else {
            let _ = writeln!(out, "    :stress-tolerance-pa {}", canonical_float(self.stress.absolute_tolerance));
            if let Some(family) = &self.family { family.canonical(out); }
            if self.regions.is_empty() {
                let _ = writeln!(out, "  )");
            } else {
                regions::canonical(&self.regions, out);
            }
        }
    }
}

fn stress_json(s: &SampledStressEvaluation) -> String {
    format!("{{\"compliance_j\":{:.17e},\"area_m2\":{:.17e},\"sampled_von_mises_pa\":{:.17e},\"max_x\":{:.17e},\"max_y\":{:.17e},\"sample_count\":{},\"snapshot\":\"{:#018x}\"}}",
        s.compliance, s.volume, s.sampled_max_von_mises, s.max_location[0], s.max_location[1], s.sample_count, s.snapshot)
}

fn read_stress(value: &JsonValue, policy: &Controls) -> Result<SampledStressEvaluation> {
    let real = |key| number(value, key).map(|(number, _)| number);
    let state = SampledStressEvaluation {
        compliance: real("compliance_j")?, volume: real("area_m2")?,
        sampled_max_von_mises: real("sampled_von_mises_pa")?,
        max_location: [real("max_x")?, real("max_y")?],
        sample_count: integer(value, "sample_count")?,
        snapshot: hexadecimal(value.str_field("snapshot").ok_or_else(|| malformed("missing stress snapshot"))?, true)?,
    };
    if state.compliance < 0.0 || state.volume <= 0.0 || state.sample_count == 0
        || state.sample_count > 100_000_000
        || state.max_location.iter().any(|v| !(0.0..=1.0).contains(v))
        || (state.volume - policy.area.target).abs() > policy.area.tolerance
        || state.sampled_max_von_mises < 0.0
        || (restoration::reduction(policy).is_none()
            && state.sampled_max_von_mises > policy.stress.admitted_max())
    { return Err(malformed("retained projected stress state violates the declared constraints")); }
    Ok(state)
}

fn same(a: &SampledStressEvaluation, b: &SampledStressEvaluation) -> bool {
    a.snapshot == b.snapshot && a.sample_count == b.sample_count
        && [a.compliance, a.volume, a.sampled_max_von_mises, a.max_location[0], a.max_location[1]]
            .iter().zip([b.compliance, b.volume, b.sampled_max_von_mises, b.max_location[0], b.max_location[1]])
            .all(|(a, b)| a.to_bits() == b.to_bits())
}

/// Small accepted-state history, bounded by the native 32-update envelope.
/// Detailed proposal attempts remain owned by fs-topols; final refusal reasons
/// are retained here rather than disguised as convergence or a completed study.
pub(super) struct ConstraintEvidence {
    policy: Controls,
    baseline: SampledStressEvaluation,
    accepted: Vec<SampledStressEvaluation>,
    attempts: Vec<usize>,
    refusals: Vec<String>,
    family: Option<multi_load::History>,
}

impl ConstraintEvidence {
    pub(super) fn current(&self) -> &SampledStressEvaluation {
        self.accepted.last().unwrap_or(&self.baseline)
    }

    pub(super) fn html(&self) -> String {
        let baseline_label = if restoration::reduction(&self.policy).is_some() {
            "Area-feasible input compliance"
        } else { "Feasible baseline compliance" };
        let mut html = format!("<p>Hard material area: {:.8e} m² ± {:.8e} m². Sampled stress limit: {:.8e} Pa + {:.8e} Pa allowance. {baseline_label}: {:.8e} J. Current sampled maximum: {:.8e} Pa. Stress is sample-scoped, not a continuous-domain bound.</p>",
            self.policy.area.target, self.policy.area.tolerance,
            self.policy.stress.max_von_mises, self.policy.stress.absolute_tolerance,
            self.baseline.compliance, self.current().sampled_max_von_mises);
        html.push_str(&restoration::html(&self.baseline, &self.accepted, &self.policy));
        if !self.policy.regions.is_empty() {
            let _ = write!(html, "<p>{} protected material/void regions are imposed on every intersected cell through its corner nodes. Coverage may extend by less than one cell per side. The phi margin is a field-value margin, not a certified physical clearance or wall thickness.</p>", self.policy.regions.len());
        }
        if let (Some(family), Some(history)) = (&self.policy.family, &self.family) {
            html.push_str(&family.html(history));
        }
        html
    }

    pub(super) fn json(&self) -> String {
        let reduction = restoration::relative_reduction(&self.baseline, &self.accepted, &self.policy);
        format!(concat!("{{\"mode\":\"projected-stress-v1\",\"baseline_scope\":{},",
            "\"area_target_m2\":{:.17e},\"area_tolerance_m2\":{:.17e},",
            "\"stress_limit_pa\":{:.17e},\"stress_tolerance_pa\":{:.17e},",
            "\"baseline\":{},\"accepted\":[{}],\"candidate_counts\":[{}],",
            "\"terminal_refusals\":[{}],\"relative_reduction\":{}{}{}{}}}"),
            quoted(restoration::baseline_scope(&self.policy)),
            self.policy.area.target, self.policy.area.tolerance,
            self.policy.stress.max_von_mises, self.policy.stress.absolute_tolerance,
            stress_json(&self.baseline), self.accepted.iter().map(stress_json).collect::<Vec<_>>().join(","),
            self.attempts.iter().map(usize::to_string).collect::<Vec<_>>().join(","),
            self.refusals.iter().map(|v| quoted(v)).collect::<Vec<_>>().join(","), reduction,
            regions::json_field(&self.policy.regions),
            self.family.as_ref().zip(self.policy.family.as_ref())
                .map_or_else(String::new, |(history, family)| history.json_field(family)),
            restoration::json_field(&self.baseline, &self.accepted, &self.policy))
    }

    fn read(value: &JsonValue, report: &OptimizeReport, policy: &Controls) -> Result<Self> {
        if value.str_field("mode") != Some("projected-stress-v1")
            || value.str_field("baseline_scope") != Some(restoration::baseline_scope(policy))
        { return Err(malformed("missing projected-stress baseline identity")); }
        regions::check_retained(value, &policy.regions)?;
        for (key, expected) in [("area_target_m2", policy.area.target),
            ("area_tolerance_m2", policy.area.tolerance),
            ("stress_limit_pa", policy.stress.max_von_mises),
            ("stress_tolerance_pa", policy.stress.absolute_tolerance)]
        {
            if number(value, key)?.0.to_bits() != expected.to_bits() {
                return Err(malformed("retained constraints differ from the canonical study"));
            }
        }
        let array = |key| value.get(key).and_then(JsonValue::as_array)
            .ok_or_else(|| malformed("missing constrained history array"));
        let baseline = read_stress(value.get("baseline").ok_or_else(|| malformed("missing feasible baseline"))?, policy)?;
        let accepted_values = array("accepted")?;
        let attempts = array("candidate_counts")?;
        let refused = array("terminal_refusals")?;
        if accepted_values.len() != report.rows.len() || attempts.len() != report.rows.len()
            || refused.len() > policy.search.max_candidates
        { return Err(malformed("constraint history lengths disagree with the accepted trajectory")); }
        let mut accepted = Vec::with_capacity(accepted_values.len());
        let mut counts = Vec::with_capacity(attempts.len());
        for (i, (state, count)) in accepted_values.iter().zip(attempts).enumerate() {
            let state = read_stress(state, policy)?;
            let count: usize = count.number_raw().and_then(|v| v.parse().ok())
                .ok_or_else(|| malformed("invalid candidate count"))?;
            if !(1..=policy.search.max_candidates).contains(&count)
                || state.snapshot != report.snapshots[i]
                || state.compliance.to_bits() != report.compliance[i].to_bits()
                || state.volume.to_bits() != report.volume[i].to_bits()
            { return Err(malformed("accepted constraint history does not match the decreasing trajectory")); }
            restoration::transition(accepted.last().unwrap_or(&baseline), &state, policy)?;
            accepted.push(state);
            counts.push(count);
        }
        let refusals = refused.iter().map(|v| v.as_str().map(str::to_string)
            .ok_or_else(|| malformed("invalid candidate refusal"))).collect::<Result<Vec<_>>>()?;
        let family = multi_load::History::read(value, &baseline, &accepted, policy)?;
        restoration::check_retained(value, &baseline, &accepted, policy)?;
        Ok(Self { policy: policy.clone(), baseline, accepted, attempts: counts, refusals, family })
    }
}

#[derive(Debug, Clone, Copy)]
enum Stage {
    Regions(DesignRegionStage),
    Setup(ProjectedStressSetupStage),
    Update(ProjectedStage),
}

fn constraints_stop(status: &'static str, last: Option<&Outcome>) -> Failure {
    retained_error(Failure {
        code: "cli-study-elasticity-constraint-stop",
        message: format!("{status} before complete constrained-state admission; no new feasible state published"),
        exit: if status == "cancelled" { exit::CANCELLED } else { exit::BUDGET },
    }, last)
}

pub(super) fn drive(spec: &ElasticitySpec, ledger: &Ledger, cap: Option<usize>,
    gate: &CancelGate, prior: Option<&Loaded>) -> Result<Outcome> {
    if spec.projected.as_ref().is_some_and(|policy| policy.family.is_some()) {
        return multi_load::drive(spec, ledger, cap, gate, prior);
    }
    drive_observed(spec, ledger, cap, gate, prior, |_| {})
}

fn drive_observed(spec: &ElasticitySpec, ledger: &Ledger, cap: Option<usize>,
    gate: &CancelGate, prior: Option<&Loaded>, mut observe: impl FnMut(Stage)) -> Result<Outcome> {
    if ledger.in_transaction() { return Err(malformed("constrained study requires its own ledger transaction")); }
    let policy = spec.projected.as_ref().ok_or_else(|| malformed("missing projected controls"))?;
    let start = Instant::now();
    let mut evidence = Evidence { producer: producer_identity()?, updates: 0, legacy_replayed: 0, projected: None };
    let mut predecessor = prior.map(|old| old.hash);
    let mut last = None;
    let mut consumed = 0.0;
    let mut report = OptimizeReport::default();
    let mut state = if let Some(old) = prior {
        if old.value.str_field("study_id") != Some(spec.id.to_hex().as_str()) {
            return Err(malformed("retained study identity changed"));
        }
        let binding = old.value.get("continuation").ok_or_else(|| malformed("constrained resume needs accepted-state continuation"))?;
        if integer(binding, "version")? != 1
            || binding.str_field("producer") != Some(evidence.producer.to_hex().as_str())
        { return Err(malformed("constrained resume requires the identical executable")); }
        let design = document(&linked(ledger, &old.value, "design", "study-design")?)?;
        let iterations = document(&linked(ledger, &old.value, "iterations", "study-iterations")?)?;
        let (phi, decoded) = decode(spec, &old.value, &design, &iterations)?;
        report = decoded;
        let retained = ConstraintEvidence::read(binding.get("constraints")
            .ok_or_else(|| malformed("missing retained constrained history"))?, &report, policy)?;
        if snapshot(&phi) != retained.current().snapshot {
            return Err(malformed("retained stress and geometry differ"));
        }
        let status = match old.value.str_field("status") {
            Some("running") => "running", Some("completed") => "completed",
            Some("cancelled") => "cancelled", Some("budget-exhausted") => "budget-exhausted",
            Some("no-feasible-descent") => "no-feasible-descent",
            _ => return Err(malformed("unknown constrained study terminal")),
        };
        last = Some(Outcome { pointer: format!("study-{}", old.hash.to_hex()), receipt: old.bytes.clone(), status });
        consumed = old.value.f64_field("consumed_wall_s").filter(|v| v.is_finite() && *v >= 0.0)
            .ok_or_else(|| malformed("invalid retained wall charge"))?;
        // Rebuild the ORIGINAL prescriptions. Re-authoring the saved endpoint
        // would conceal a changed fixed value or silently repair a violation.
        let prepared = regions::prepare(spec, &policy.regions, |stage| {
            observe(Stage::Regions(stage));
            match stop_status(gate.is_requested(), consumed + start.elapsed().as_secs_f64(), spec.wall_s) {
                Some(status) => ControlFlow::Break(status), None => ControlFlow::Continue(()),
            }
        }).map_err(|error| retained_error(error, last.as_ref()))?;
        let prepared = match prepared {
            ControlFlow::Continue(prepared) => prepared,
            ControlFlow::Break(status) => return Err(constraints_stop(status, last.as_ref())),
        };
        let checkpoint = OptimizeCheckpoint::restore(phi, fixture(spec), settings(spec, spec.steps),
            report.rows.len(), report.ell.last().copied().unwrap_or(spec.ell0))
            .map_err(|error| malformed(&error.to_string()))?;
        let restored = ProjectedStressOptimizer::from_checkpoint_controlled(&checkpoint, prepared.fixed_nodes,
            policy.area, policy.search, policy.stress, |stage| {
                observe(Stage::Setup(stage));
                match stop_status(gate.is_requested(), consumed + start.elapsed().as_secs_f64(), spec.wall_s) {
                    Some(status) => ControlFlow::Break(status), None => ControlFlow::Continue(()),
                }
            }).map_err(|error| retained_error(malformed(&error.to_string()), last.as_ref()))?;
        let restored = match restored {
            ControlFlow::Continue(restored) => restored,
            ControlFlow::Break(status) => return Err(constraints_stop(status, last.as_ref())),
        };
        if !same(restored.current(), retained.current()) {
            return Err(malformed("independent constrained replay differs from the retained endpoint"));
        }
        evidence.projected = Some(retained);
        if status == "completed" || status == "no-feasible-descent" {
            return last.ok_or_else(|| malformed("terminal constrained study has no receipt"));
        }
        restored
    } else {
        let prepared = regions::prepare(spec, &policy.regions, |stage| {
            observe(Stage::Regions(stage));
            match stop_status(gate.is_requested(), start.elapsed().as_secs_f64(), spec.wall_s) {
                Some(status) => ControlFlow::Break(status), None => ControlFlow::Continue(()),
            }
        })?;
        let prepared = match prepared {
            ControlFlow::Continue(prepared) => prepared,
            ControlFlow::Break(status) => return Err(constraints_stop(status, None)),
        };
        let area = ProjectedOptimizer::new_controlled(&prepared.geometry, fixture(spec), settings(spec, spec.steps),
            prepared.fixed_nodes, policy.area, policy.search, |stage| {
                observe(Stage::Setup(ProjectedStressSetupStage::Area(stage)));
                match stop_status(gate.is_requested(), start.elapsed().as_secs_f64(), spec.wall_s) {
                    Some(status) => ControlFlow::Break(status), None => ControlFlow::Continue(()),
                }
            }).map_err(|error| malformed(&error.to_string()))?;
        let area = match area { ControlFlow::Continue(area) => area,
            ControlFlow::Break(status) => return Err(constraints_stop(status, None)) };
        let state = ProjectedStressOptimizer::new_controlled(&area, policy.stress, |stage| {
            observe(Stage::Setup(stage));
            match stop_status(gate.is_requested(), start.elapsed().as_secs_f64(), spec.wall_s) {
                Some(status) => ControlFlow::Break(status), None => ControlFlow::Continue(()),
            }
        }).map_err(|error| malformed(&error.to_string()))?;
        let state = match state { ControlFlow::Continue(state) => state,
            ControlFlow::Break(status) => return Err(constraints_stop(status, None)) };
        evidence.projected = Some(ConstraintEvidence { policy: policy.clone(), baseline: state.current().clone(),
            accepted: Vec::new(), attempts: Vec::new(), refusals: Vec::new(), family: None });
        state
    };
    let target = spec.steps.min(report.rows.len().saturating_add(cap.unwrap_or(spec.steps - report.rows.len())));
    if last.is_none() {
        if let Some(status) = stop_status(gate.is_requested(), start.elapsed().as_secs_f64(), spec.wall_s) {
            return Err(constraints_stop(status, None));
        }
        let initial = persist(spec, ledger, state.checkpoint().geometry(), &report, "running",
            consumed + start.elapsed().as_secs_f64(), predecessor, &evidence)?;
        predecessor = initial.pointer.strip_prefix("study-").and_then(ContentHash::from_hex);
        last = Some(initial);
    }
    loop {
        let status = stop_status(gate.is_requested(), consumed + start.elapsed().as_secs_f64(), spec.wall_s)
            .or_else(|| (state.checkpoint().next_iteration() == target)
                .then_some(if state.checkpoint().is_complete() { "completed" } else { "budget-exhausted" }));
        if let Some(status) = status {
            return persist(spec, ledger, state.checkpoint().geometry(), &report, status,
                consumed + start.elapsed().as_secs_f64(), predecessor, &evidence)
                .map_err(|error| retained_error(error, last.as_ref()));
        }
        let update = state.advance_one_controlled(|stage| {
            observe(Stage::Update(stage));
            match stop_status(gate.is_requested(), consumed + start.elapsed().as_secs_f64(), spec.wall_s) {
                Some(status) => ControlFlow::Break(status), None => ControlFlow::Continue(()),
            }
        }).map_err(|error| retained_error(malformed(&error.to_string()), last.as_ref()))?;
        let update = match update {
            ControlFlow::Continue(update) => update,
            ControlFlow::Break(status) => return persist(spec, ledger, state.checkpoint().geometry(), &report,
                status, consumed + start.elapsed().as_secs_f64(), predecessor, &evidence)
                .map_err(|error| retained_error(error, last.as_ref())),
        };
        let retained = evidence.projected.as_mut().expect("admitted constrained state");
        match update.progress {
            ProjectedProgress::Accepted(step) => {
                let current = state.current();
                if step.iteration != report.rows.len() || step.state.snapshot != current.snapshot
                    || step.state.compliance.to_bits() != current.compliance.to_bits()
                    || step.state.volume.to_bits() != current.volume.to_bits()
                { return Err(retained_error(malformed("accepted mechanics and stress are inconsistent"), last.as_ref())); }
                let audit = step.proposal.audits.first().ok_or_else(|| malformed("proposal audit is missing"))?;
                let drift = audit.interface_drift_h;
                let pad = step.proposal.load_pad_nodes.first().copied().ok_or_else(|| malformed("proposal load-pad count is missing"))?;
                let ell = state.checkpoint().ell();
                report.rows.push(format!("{{\"iter\":{},\"compliance\":{:.17e},\"volume\":{:.17e},\"ell\":{ell:.17e},\"drift_h\":{drift:.17e},\"load_pad_nodes\":{pad},\"snapshot\":\"{:#018x}\"}}",
                    step.iteration, current.compliance, current.volume, current.snapshot));
                report.compliance.push(current.compliance);
                report.volume.push(current.volume);
                report.ell.push(ell);
                report.snapshots.push(current.snapshot);
                report.load_pad_nodes.push(pad);
                retained.accepted.push(current.clone());
                retained.attempts.push(step.attempts.len());
                retained.refusals.clear();
                evidence.updates += 1;
            }
            ProjectedProgress::NoDescent(attempts) => {
                retained.refusals = attempts.into_iter().map(|attempt| attempt.refusal
                    .unwrap_or_else(|| "candidate was not admitted".into())).collect();
                return persist(spec, ledger, state.checkpoint().geometry(), &report, "no-feasible-descent",
                    consumed + start.elapsed().as_secs_f64(), predecessor, &evidence)
                    .map_err(|error| retained_error(error, last.as_ref()));
            }
            ProjectedProgress::IterationLimit => return Err(malformed("constrained optimizer completed before target")),
        }
        let accepted = persist(spec, ledger, state.checkpoint().geometry(), &report, "running",
            consumed + start.elapsed().as_secs_f64(), predecessor, &evidence)
            .map_err(|error| retained_error(error, last.as_ref()))?;
        predecessor = accepted.pointer.strip_prefix("study-").and_then(ContentHash::from_hex);
        last = Some(accepted);
    }
}

#[cfg(test)]
#[path = "projected/tests.rs"]
mod tests;
