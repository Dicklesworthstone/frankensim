//! Whole-load-case selection using the completed physical observation Jacobian.
use super::{Result, CancelGate, EquilibriumDesignFile, ObservationEvaluation, quoted};
use std::fmt::Write as _;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::EquilibriumDesign;
use fs_uq::experimental_design::{ExperimentDesignConfig, ExperimentDesignStop, ExperimentGroup, select_experiments};

#[derive(Clone, Copy)]
pub(super) struct Options { count: usize, ridge: f64, factorizations: usize }
impl Options {
    pub(super) fn new(count: usize, ridge: f64, factorizations: usize) -> Result<Self> {
        let options = Self { count, ridge, factorizations };
        options.config(1).admit(64, 1)?;
        Ok(options)
    }
    fn config(self, parameters: usize) -> ExperimentDesignConfig {
        ExperimentDesignConfig { parameters, maximum_selected: self.count,
            ridge_precision: self.ridge, maximum_factorizations: self.factorizations }
    }
    pub(super) fn admit(self, problem: &EquilibriumDesign) -> Result<()> {
        let rows = problem.load_cases().iter().try_fold(0usize, |n, case| n.checked_add(case.targets.len()))
            .ok_or("observation count overflow")?;
        self.config(problem.variables().len()).admit(problem.load_cases().len(), rows)?;
        Ok(())
    }
}

pub(super) fn run(loaded: &EquilibriumDesignFile, data: &ObservationEvaluation, options: Options,
    gate: &CancelGate) -> Result<String>
{
    let problem = loaded.problem();
    options.admit(problem)?;
    let mut groups: Vec<ExperimentGroup> = problem.load_cases().iter().map(|case|
        ExperimentGroup { name: case.name.clone(), rows: Vec::with_capacity(case.targets.len()) }).collect();
    for row in data.rows() {
        if gate.is_requested() { return Err("experiment selection cancelled".into()); }
        groups.get_mut(row.case).ok_or("unknown observation case")?.rows.push(row.residual_gradient.clone());
    }
    let report = select_experiments(&groups, options.config(data.point().len()), || gate.is_requested())?;
    let stop = match report.stop {
        ExperimentDesignStop::SelectionLimit => "selection-limit",
        ExperimentDesignStop::EvaluationLimit => "factorization-budget",
        ExperimentDesignStop::NoResolvedGain => "no-resolved-gain",
    };
    let mut out = format!("{{\"method\":\"greedy-whole-case-logdet\",\"stop_reason\":{},\"requested_cases\":{},\"selected_count\":{},\"ridge_precision\":{:.17e},\"maximum_factorizations\":{},\"factorizations\":{},\"log_determinant_gain\":{:.17e},\"criterion\":\"log det(I + sum selected J_case^T J_case / ridge_precision)\",\"scope\":\"Local weighted-Jacobian heuristic at the declared point, not a global optimal subset, calibrated-noise model, posterior or rank certificate. The top-level information describes ALL candidate cases, not the selected subset. No physical constraint is waived or newly certified.\",\"candidate_physics_reused\":true,\"selected_cases\":[",
        quoted(stop), options.count, report.selected.len(), options.ridge,
        options.factorizations, report.factorizations, report.log_determinant_gain);
    for (ordinal, step) in report.steps.iter().enumerate() {
        if ordinal != 0 { out.push(','); }
        write!(&mut out, "{{\"case_index\":{},\"case\":{},\"marginal_gain\":{:.17e},\"cumulative_gain\":{:.17e},\"observation_indices\":[",
            step.group, quoted(&groups[step.group].name), step.marginal_gain, step.log_determinant_gain).expect("String write");
        let mut first = true;
        for (index, _) in data.rows().iter().enumerate().filter(|(_, r)| r.case == step.group) {
            if !first { out.push(','); } first = false;
            write!(&mut out, "{index}").expect("String write");
        }
        out.push_str("]}");
    }
    out.push_str("]}");
    if gate.is_requested() { return Err("experiment selection cancelled".into()); }
    Ok(out)
}
