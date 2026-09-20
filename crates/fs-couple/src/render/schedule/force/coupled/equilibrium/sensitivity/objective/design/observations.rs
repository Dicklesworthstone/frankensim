//! Observation-wise calibration derivatives over the original equilibrium owner.
//! An aggregate objective gradient is not an observation Jacobian: it vanishes
//! at an exact fit even when each observation remains sensitive to parameters.
use super::*;
use super::constraints::{ConstraintSense, ResponseQuantity};

/// Bounded output and caller-owned cumulative observation-adjoint allowance.
/// All targets count, including targets with zero objective weight. Counts are
/// attempts, not flops. A failed query never refunds attempted adjoints.
#[derive(Clone, Debug)]
pub struct ObservationControl {
    maximum_rows: usize,
    maximum_adjoints: usize,
    adjoints: usize,
}
impl ObservationControl {
    /// Each query admits at most 1024 rows and 128 decision columns. The
    /// original design control independently limits evaluations and cases.
    pub fn new(maximum_rows: usize, maximum_adjoints: usize) -> Result<Self, DesignError> {
        if !(1..=1024).contains(&maximum_rows) {
            return Err(bad("observation row cap must be in 1..=1024"));
        }
        Ok(Self { maximum_rows, maximum_adjoints, adjoints: 0 })
    }
    /// Cumulative attempted observation adjoints, including returned failures.
    #[must_use]
    pub const fn adjoints_attempted(&self) -> usize { self.adjoints }
    /// Raise the allowance without resetting work or changing the output cap.
    pub fn extend(&mut self, maximum_adjoints: usize) -> Result<(), DesignError> {
        if maximum_adjoints < self.maximum_adjoints {
            return Err(bad("observation adjoint allowance may only increase"));
        }
        self.maximum_adjoints = maximum_adjoints;
        Ok(())
    }
}

/// One physical observation and one weighted least-squares residual row.
#[derive(Clone, Debug, PartialEq)]
pub struct ObservationRow {
    /// Original load-case index.
    pub case: usize,
    /// Original target index inside that case, even for zero-weight targets.
    pub target: usize,
    /// Physical displacement [m].
    pub value_m: f64,
    /// d(displacement)/dx, in metres per dimensionless decision coordinate.
    /// Includes shared-field sums and case-local physical load derivatives.
    pub gradient_m: Vec<f64>,
    /// sqrt(weight)*(displacement-target_m)/scale_m; zero for zero weight.
    pub residual: f64,
    /// Derivative of `residual` in the same declared decision coordinates.
    pub residual_gradient: Vec<f64>,
    /// Recomputed full-coordinate residual from the original adjoint solve.
    pub adjoint_relative_residual: f64,
}

/// Complete case-major/target-major rows for one immutable physical candidate.
/// This does not assess response constraints or change optimizer state. The
/// private container prevents caller mutation invalidating derived operations.
#[derive(Clone, Debug, PartialEq)]
pub struct ObservationEvaluation {
    point: Vec<f64>,
    physical_parameters: Vec<f64>,
    rows: Vec<ObservationRow>,
    equilibria: Vec<EquilibriumLinearizationReport>,
    unassessed_response_constraints: usize,
    objective: f64,
    gradient: Vec<f64>,
}
impl ObservationEvaluation {
    /// Complete dimensionless decision point, not an optimizer's next trial.
    #[must_use] pub fn point(&self) -> &[f64] { &self.point }
    /// Affine-decoded physical parameters in variable declaration order.
    #[must_use] pub fn physical_parameters(&self) -> &[f64] { &self.physical_parameters }
    /// Complete observation rows; none are suppressed by fitting weights.
    #[must_use] pub fn rows(&self) -> &[ObservationRow] { &self.rows }
    /// Original force-balance/activity reports in load-case order.
    #[must_use] pub fn equilibria(&self) -> &[EquilibriumLinearizationReport] { &self.equilibria }
    /// Requirements not assessed by this calibration query (not presumed met).
    #[must_use] pub const fn unassessed_response_constraints(&self) -> usize { self.unassessed_response_constraints }
    /// 0.5*sum(residual^2), reconstructed from the weighted observation rows.
    /// Algebraically the original objective; a different reduction can round differently.
    #[must_use] pub const fn objective(&self) -> f64 { self.objective }
    /// J^T r, not a claim about parameter uniqueness or experimental noise.
    #[must_use] pub fn gradient(&self) -> &[f64] { &self.gradient }

    /// J^T J v without a square matrix or further physics evaluations.
    /// This is the Gauss-Newton action, NOT the full Hessian away from zero
    /// residual: sum(r_i * Hessian(r_i)) is intentionally not supplied.
    pub fn gauss_newton_product(&self, direction: &[f64], gate: &CancelGate)
        -> Result<Vec<f64>, DesignError>
    {
        checkpoint(gate)?;
        if direction.len() != self.point.len() || direction.iter().any(|v| !v.is_finite()) {
            return Err(bad("Gauss-Newton direction must match every finite decision coordinate"));
        }
        let mut product = vec![0.0; direction.len()];
        for row in &self.rows {
            checkpoint(gate)?;
            let projection = dot(&row.residual_gradient, direction).map_err(|e| case_error(row.case, e))?;
            for (value, derivative) in product.iter_mut().zip(&row.residual_gradient) {
                *value = number(*value + number(derivative * projection)?)?;
            }
        }
        checkpoint(gate)?;
        Ok(product)
    }
}

impl EquilibriumDesign {
    /// Differentiate every declared displacement observation, including an
    /// exact fit or zero-weight target. One primal and one tangent preparation
    /// per load case; one existing response adjoint per observation. No finite
    /// differences, objective-gradient division, perturbation solves or new
    /// mechanics. Contact-margin and derivative work gates remain unchanged.
    ///
    /// Admit the whole row/adjoint family before physics. `DesignControl`
    /// retains its original case/evaluation accounting; `ObservationControl`
    /// additionally counts attempted row adjoints, including failures. A late
    /// error returns no partial matrix and leaves both ledgers charged.
    /// Weighted residuals encode fitting priorities, NOT inferred noise precision.
    pub fn evaluate_observations(&self, point: &[f64], control: &mut DesignControl,
        observations: &mut ObservationControl, gate: &CancelGate)
        -> Result<ObservationEvaluation, DesignError>
    {
        checkpoint(gate)?;
        let count = self.cases.iter().try_fold(0usize, |n, case| n.checked_add(case.targets.len()))
            .ok_or_else(|| bad("observation row count overflow"))?;
        if count > observations.maximum_rows
            || count > observations.maximum_adjoints.saturating_sub(observations.adjoints) {
            return Err(DesignError::Budget { what: "complete observation/adjoint family" });
        }
        let prepared = self.prepare(point, control, gate)?;
        let mut result = ObservationEvaluation { point: point.to_vec(), physical_parameters: prepared.parameters.clone(),
            rows: Vec::with_capacity(count), equilibria: Vec::with_capacity(prepared.cases.len()),
            unassessed_response_constraints: self.constraints.len(), objective: 0.0, gradient: vec![0.0; point.len()] };
        for (index, case) in prepared.cases.iter().enumerate() {
            let solved = self.solve_case(&prepared, index, control, gate)?;
            let linear = EquilibriumLinearization::new(&solved.network, &solved.external,
                &prepared.contacts, self.budget.sensitivity, gate).map_err(|e| case_error(index, e))?;
            for (target_index, target) in case.targets.iter().enumerate() {
                checkpoint(gate)?;
                // This is a temporary linear response query, not a new design
                // constraint. Unit normalization preserves the raw dy/dx even
                // when the objective weight is zero or the fitting residual vanishes.
                let response = ResponseConstraint { name: String::new(), case: index,
                    quantity: ResponseQuantity::Displacement(target.attachment.clone()),
                    sense: ConstraintSense::Equal, bound: 0.0, scale: 1.0 };
                observations.adjoints += 1;
                let raw = constraints::evaluate(&response, &linear, &prepared.contacts, case, &self.variables, gate)
                    .map_err(|e| case_error(index, e))?;
                let mut residual = 0.0;
                let mut residual_gradient = vec![0.0; point.len()];
                if target.weight > 0.0 {
                    let root = target.weight.sqrt();
                    residual = number(root * number(number(raw.value - target.target_m)? / target.scale_m)?)?;
                    for (value, derivative) in residual_gradient.iter_mut().zip(&raw.gradient) {
                        *value = number(root * number(derivative / target.scale_m)?)?;
                    }
                }
                result.objective = number(result.objective + number(0.5 * residual * residual)?)?;
                for (g, derivative) in result.gradient.iter_mut().zip(&residual_gradient) {
                    *g = number(*g + number(derivative * residual)?)?;
                }
                result.rows.push(ObservationRow { case: index, target: target_index, value_m: raw.value,
                    gradient_m: raw.gradient, residual, residual_gradient,
                    adjoint_relative_residual: raw.adjoint_relative_residual });
            }
            result.equilibria.push(linear.report());
        }
        checkpoint(gate)?;
        Ok(result)
    }
}
