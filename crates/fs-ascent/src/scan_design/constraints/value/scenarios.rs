//! One authored heating design evaluated against a finite set of physical cases.
//!
//! Cases may change their meshes, material charts, initial states, boundaries,
//! observations and requirements, but not the authored source/control map. All
//! inequalities survive in case order; averaging a loss never averages away a
//! violation. This is finite-case design, not a probability-of-failure bound or
//! a certificate over an unsampled uncertainty set.

use super::{ValueConstrainedScanStudy, ScanConstraintValueBatch, ScanConstraintValues,
    ScanPowerValueEvaluation, ScanObservation, ScanObservationSpec, ScanSensitivityError,
    ScanningEnthalpyRun, MovingHeatSource};
use fs_exec::Cx;
use std::{collections::BTreeSet, fmt};

/// A caller-declared case. Weights are positive preferences, not inferred probabilities.
pub struct ScanScenario<F, O = ScanObservation> {
    pub name: String,
    pub weight: f64,
    pub study: ValueConstrainedScanStudy<F, O>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScanScenarioObjective {
    /// Sum(weight * original loss) / sum(weight), using overflow-safe scaling.
    WeightedMean,
    /// Largest original case loss. Weights do not alter this policy.
    WorstCase,
}

/// Aggregate shape limits. Each physical producer retains its own solve budgets.
/// Vertices count retained final fields, not bytes or allocator/workspace memory.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanScenarioBudget {
    pub max_cases: usize,
    pub max_constraint_rows: usize,
    pub max_retained_vertices: usize,
    pub max_name_bytes: usize,
}

#[derive(Debug)]
pub enum ScanScenarioError {
    Invalid(&'static str),
    Cancelled,
    /// Original physical error and its zero-based declared case index.
    Physical { case: usize, source: ScanSensitivityError },
}
impl fmt::Display for ScanScenarioError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(message) => f.write_str(message),
            Self::Cancelled => f.write_str("thermal case evaluation cancelled"),
            Self::Physical { case, source } => write!(f, "thermal case {case}: {source}"),
        }
    }
}
impl std::error::Error for ScanScenarioError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Physical { source, .. } => Some(source), _ => None }
    }
}

#[derive(Debug)]
pub struct ScanScenarioCaseEvaluation {
    pub name: String,
    pub weight: f64,
    pub physical: ScanPowerValueEvaluation,
    pub constraints: ScanConstraintValues,
}

#[derive(Debug)]
pub struct ScanScenarioEvaluation {
    pub objective: f64,
    pub policy: ScanScenarioObjective,
    /// Earliest case attaining the largest original fitting loss.
    pub worst_objective_case: usize,
    pub cases: Vec<ScanScenarioCaseEvaluation>,
    /// EVERY signed inequality, concatenated without reweighting or rescaling.
    pub residuals: Vec<f64>,
    /// Case i owns residuals[offsets[i]..offsets[i+1]]. Includes final sentinel.
    pub constraint_offsets: Vec<usize>,
}

/// Immutable case set, common design coordinates and original physical budgets.
/// No case may be removed, reweighted or replaced during an optimization run.
pub struct ScanScenarioStudy<F, O = ScanObservation> {
    cases: Vec<ValueConstrainedScanStudy<F, O>>,
    names: Vec<String>,
    weights: Vec<f64>,
    rows: Vec<usize>,
    policy: ScanScenarioObjective,
    budget: ScanScenarioBudget,
    constraint_count: usize,
}

pub struct ScanScenarioValueBatch<'a, F, O = ScanObservation> {
    cases: Vec<ScanConstraintValueBatch<'a, F, O>>,
    names: &'a [String],
    weights: &'a [f64],
    rows: &'a [usize],
    policy: ScanScenarioObjective,
    budget: ScanScenarioBudget,
}

impl<F, O: ScanObservationSpec> ScanScenarioStudy<F, O> {
    /// Admit without invoking any physical factory. Identical coordinate counts
    /// alone are insufficient: compare sources, selected entries, unit scales,
    /// coordinate maps and the source-validation instant as well.
    pub fn new(cx: &Cx<'_>, cases: Vec<ScanScenario<F, O>>, policy: ScanScenarioObjective,
        budget: ScanScenarioBudget) -> Result<Self, ScanScenarioError>
    {
        check(cx)?;
        if cases.is_empty() || cases.len() > budget.max_cases {
            return Err(invalid("thermal cases must be nonempty and fit their case allowance"));
        }
        let first = &cases[0].study.study;
        let mut names_seen = BTreeSet::new();
        let (mut name_bytes, mut constraint_count) = (0usize, 0usize);
        for case in &cases {
            check(cx)?;
            name_bytes = name_bytes.checked_add(case.name.len()).filter(|&n| n <= budget.max_name_bytes)
                .ok_or_else(|| invalid("thermal case names exceed their byte allowance"))?;
            if case.name.trim().is_empty() || !names_seen.insert(case.name.as_str())
                || !case.weight.is_finite() || case.weight <= 0.0 {
                return Err(invalid("thermal cases need unique nonempty names and positive finite weights"));
            }
            let other = &case.study.study;
            if first.sources != other.sources || first.controls != other.controls
                || first.geometry_controls != other.geometry_controls
                || first.power_coordinates != other.power_coordinates
                || first.objective.source_reference_time_s.to_bits() != other.objective.source_reference_time_s.to_bits()
                || !same_words(&first.initial, &other.initial) {
                return Err(invalid("thermal cases must use the same authored sources and complete design-coordinate map"));
            }
            constraint_count = constraint_count.checked_add(case.study.constraint_count())
                .filter(|&n| n <= budget.max_constraint_rows)
                .ok_or_else(|| invalid("thermal case inequalities exceed their aggregate allowance"))?;
        }
        drop(names_seen);
        let mut studies = reserve(cases.len())?;
        let mut names = reserve(cases.len())?;
        let mut weights = reserve(cases.len())?;
        let mut rows = reserve(cases.len())?;
        for case in cases {
            check(cx)?;
            rows.push(case.study.constraint_count());
            studies.push(case.study);
            names.push(case.name);
            weights.push(case.weight);
        }
        // Reject unrepresentable relative weights rather than silently dropping a case.
        let maximum = weights.iter().copied().fold(0.0_f64, f64::max);
        if weights.iter().any(|w| w / maximum == 0.0) {
            return Err(invalid("thermal case weight ratio underflows"));
        }
        Ok(Self { cases: studies, names, weights, rows, policy, budget, constraint_count })
    }

    pub fn initial_coordinates(&self) -> &[f64] { self.cases[0].initial_coordinates() }
    pub fn case_names(&self) -> &[String] { &self.names }
    pub fn weights(&self) -> &[f64] { &self.weights }
    pub const fn policy(&self) -> ScanScenarioObjective { self.policy }
    pub const fn constraint_count(&self) -> usize { self.constraint_count }
    pub fn evaluations(&self) -> Vec<usize> { self.cases.iter().map(|s| s.evaluations()).collect() }
    pub fn evaluation_limits(&self) -> Vec<usize> { self.cases.iter().map(|s| s.evaluation_limit()).collect() }
    /// Number of complete candidates affordable by EVERY case. A failed late
    /// case leaves earlier work charged, so individual counters may differ.
    pub fn remaining_evaluations(&self) -> usize {
        self.cases.iter().map(|s| s.evaluation_limit().saturating_sub(s.evaluations())).min().unwrap_or(0)
    }

    /// Validate all totals before changing any original allowance. Prior probes
    /// and partially failed multi-case evaluations are never refunded.
    pub fn extend_evaluation_budgets(&mut self, totals: &[usize]) -> Result<(), ScanScenarioError> {
        if totals.len() != self.cases.len()
            || totals.iter().zip(&self.cases).any(|(&n,s)| n < s.evaluation_limit())
            || totals.iter().zip(&self.cases).all(|(&n,s)| n == s.evaluation_limit()) {
            return Err(invalid("case allowances may only increase and must name every case"));
        }
        for (case, (study, &total)) in self.cases.iter_mut().zip(totals).enumerate() {
            if total > study.evaluation_limit() {
                study.extend_evaluation_budget(total).map_err(|source| ScanScenarioError::Physical { case, source })?;
            }
        }
        Ok(())
    }

    /// One common physical source declaration. Export is not a feasibility claim.
    pub fn sources_at(&self, cx: &Cx<'_>, point: &[f64]) -> Result<Vec<MovingHeatSource>, ScanScenarioError> {
        self.cases[0].sources_at(cx, point).map_err(|source| ScanScenarioError::Physical { case: 0, source })
    }

    /// Run every case through its original producer. No partial aggregate or
    /// feasibility score escapes failure. Completed work remains charged; a
    /// retry recomputes this candidate's cases, not just its failed suffix.
    pub fn evaluate_at<'h, 'm, 'c>(&mut self, cx: &Cx<'_>, point: &[f64])
        -> Result<ScanScenarioEvaluation, ScanScenarioError>
    where F: FnMut(&Cx<'_>, Vec<MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
    {
        check(cx)?;
        if self.remaining_evaluations() == 0 { return Err(invalid("a thermal case has exhausted its physical allowance")); }
        collect(cx, &self.names, &self.weights, &self.rows, self.policy, self.budget,
            |case| self.cases[case].evaluate_at(cx, point)
                .map_err(|source| ScanScenarioError::Physical { case, source }))
    }

    /// Borrow all original studies for concurrent CANDIDATES, evaluating cases
    /// in declaration order within each callback. Per-case unused reservations
    /// settle on drop, including failure/unwinding. No fresh worker budgets.
    pub fn value_batch(&mut self, cx: &Cx<'_>, allowance: usize)
        -> Result<ScanScenarioValueBatch<'_, F, O>, ScanScenarioError>
    {
        check(cx)?;
        if allowance > self.remaining_evaluations() { return Err(invalid("thermal case batch exceeds a remaining physical allowance")); }
        let mut cases = reserve(self.cases.len())?;
        for (case, study) in self.cases.iter_mut().enumerate() {
            cases.push(study.value_batch(cx, allowance)
                .map_err(|source| ScanScenarioError::Physical { case, source })?);
        }
        Ok(ScanScenarioValueBatch { cases, names: &self.names, weights: &self.weights,
            rows: &self.rows, policy: self.policy, budget: self.budget })
    }
}

impl<F, O: ScanObservationSpec> ScanScenarioValueBatch<'_, F, O> {
    pub fn evaluate_at<'h, 'm, 'c>(&self, cx: &Cx<'_>, point: &[f64])
        -> Result<ScanScenarioEvaluation, ScanScenarioError>
    where F: Fn(&Cx<'_>, Vec<MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
    {
        collect(cx, self.names, self.weights, self.rows, self.policy, self.budget,
            |case| self.cases[case].evaluate_at(cx, point)
                .map_err(|source| ScanScenarioError::Physical { case, source }))
    }
}

#[allow(clippy::too_many_arguments)]
fn collect(cx: &Cx<'_>, names: &[String], weights: &[f64], rows: &[usize],
    policy: ScanScenarioObjective, budget: ScanScenarioBudget,
    mut evaluate: impl FnMut(usize) -> Result<(ScanPowerValueEvaluation, ScanConstraintValues), ScanScenarioError>)
    -> Result<ScanScenarioEvaluation, ScanScenarioError>
{
    check(cx)?;
    let row_count = rows.iter().try_fold(0usize, |a,&b| a.checked_add(b))
        .ok_or_else(|| invalid("thermal case row count overflow"))?;
    let mut cases = reserve(names.len())?;
    let mut residuals = reserve(row_count)?;
    let mut offsets = reserve(names.len().checked_add(1).ok_or_else(|| invalid("case count overflow"))?)?;
    let mut losses = reserve(names.len())?;
    let mut vertices = 0usize;
    offsets.push(0);
    for case in 0..names.len() {
        check(cx)?;
        let (physical, constraints) = evaluate(case)?;
        vertices = vertices.checked_add(physical.physical.final_step.solution().temperature.len())
            .filter(|&n| n <= budget.max_retained_vertices)
            .ok_or_else(|| invalid("thermal case final fields exceed their retained-vertex allowance"))?;
        if constraints.residuals.len() != rows[case] { return Err(invalid("thermal case changed its constraint layout")); }
        for (i, &row) in constraints.residuals.iter().enumerate() {
            if i % 256 == 0 { check(cx)?; }
            if !row.is_finite() { return Err(invalid("thermal case returned a nonfinite inequality")); }
            residuals.push(row);
        }
        offsets.push(residuals.len());
        losses.push(physical.physical.value);
        cases.push(ScanScenarioCaseEvaluation { name: names[case].clone(), weight: weights[case], physical, constraints });
    }
    let (objective, worst_objective_case) = aggregate(&losses, weights, policy)?;
    check(cx)?;
    Ok(ScanScenarioEvaluation { objective, policy, worst_objective_case, cases,
        residuals, constraint_offsets: offsets })
}

// Scale both the weights and the nonnegative losses. Each numerator term is
// <= its denominator term, preventing a finite weighted mean from overflowing.
fn aggregate(losses: &[f64], weights: &[f64], policy: ScanScenarioObjective)
    -> Result<(f64, usize), ScanScenarioError>
{
    if losses.is_empty() || losses.len() != weights.len()
        || losses.iter().any(|v| !v.is_finite() || *v < 0.0)
        || weights.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        return Err(invalid("case aggregation requires finite nonnegative losses and positive weights"));
    }
    let mut worst = 0;
    for i in 1..losses.len() { if losses[i] > losses[worst] { worst = i; } }
    let maximum = losses[worst];
    if policy == ScanScenarioObjective::WorstCase || maximum == 0.0 { return Ok((maximum, worst)); }
    let weight_scale = weights.iter().copied().fold(0.0_f64, f64::max);
    let (mut numerator, mut denominator) = (0.0, 0.0);
    for (&loss, &weight) in losses.iter().zip(weights) {
        let w = weight / weight_scale;
        if w == 0.0 { return Err(invalid("thermal case weight ratio underflows")); }
        numerator += w * (loss / maximum);
        denominator += w;
    }
    let value = maximum * (numerator / denominator);
    if !value.is_finite() { return Err(invalid("thermal case objective is not representable")); }
    Ok((value, worst))
}
fn same_words(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(a,b)| a.to_bits() == b.to_bits())
}
fn reserve<T>(n: usize) -> Result<Vec<T>, ScanScenarioError> {
    let mut values = Vec::new();
    values.try_reserve_exact(n).map_err(|_| invalid("thermal case allocation refused"))?;
    Ok(values)
}
fn invalid(message: &'static str) -> ScanScenarioError { ScanScenarioError::Invalid(message) }
fn check(cx: &Cx<'_>) -> Result<(), ScanScenarioError> { cx.checkpoint().map_err(|_| ScanScenarioError::Cancelled) }

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g0_scaled_means_preserve_extreme_finite_losses_and_weights() {
        assert_eq!(aggregate(&[f64::MAX, f64::MAX], &[f64::MAX, f64::MAX],
            ScanScenarioObjective::WeightedMean).unwrap(), (f64::MAX, 0));
        let (value, worst) = aggregate(&[0., 12.], &[1., 3.], ScanScenarioObjective::WeightedMean).unwrap();
        assert_eq!((value, worst), (9., 1));
        assert_eq!(aggregate(&[2., 2.], &[1., 1.], ScanScenarioObjective::WorstCase).unwrap(), (2., 0));
        assert_eq!(aggregate(&[0., 0.], &[1., 2.], ScanScenarioObjective::WeightedMean).unwrap(), (0., 0));
    }
    #[test]
    fn g0_bad_losses_and_vanishing_case_weights_never_become_valid_scores() {
        for losses in [vec![], vec![f64::NAN], vec![-1.], vec![f64::INFINITY]] {
            assert!(aggregate(&losses, &[1.], ScanScenarioObjective::WeightedMean).is_err());
        }
        assert!(aggregate(&[1., 2.], &[f64::from_bits(1), f64::MAX], ScanScenarioObjective::WeightedMean).is_err());
        assert!(aggregate(&[1.], &[0.], ScanScenarioObjective::WeightedMean).is_err());
    }
}
