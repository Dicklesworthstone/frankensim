//! Exact quadratic response-loss differences under prescribed boundary motion.
//! The fine adjoint linearizes the loss at the solved fine field. Its Taylor
//! remainder is -1/2 sum_j w_j ((q_j^T(u_f-Pu_c))/s_j)^2, NOT zero. Localize
//! that remainder as well as DWR, so an exactly fitted fine target does not
//! erase the refinement signal merely because its fine adjoint vanishes.
use super::*;
use fs_cutfem::elastic3::adaptive::{AdaptiveElasticity3, enrichment::CellResidual3};

/// One mesh-independent physical experiment; g is imposed on retained supports.
/// The RHS is f+b_g(current scales), not an external force or actuator work.
#[derive(Clone, Copy)]
pub struct EquilibriumLoad3<'a> {
    pub external: ReferenceLoad3<'a>,
    pub prescribed: Option<&'a dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
}
impl EquilibriumLoad3<'_> {
    /// Reuse the actual reference-force and current-density lifting integrators.
    /// No field solve or model mutation; returned failures expose no partial RHS.
    pub fn assemble(self, op: &AdaptiveElasticity3, mut checkpoint: impl FnMut() -> ControlFlow<()>)
        -> Result<Vec<f64>, GoalError3> {
        poll(&mut checkpoint)?;
        if self.prescribed.is_some() && op.embedded_dirichlet_penalty().is_none() {
            return Err(GoalError3::Invalid("prescribed motion requires embedded support"));
        }
        let mut rhs = op.reference_load(self.external, &mut checkpoint)?;
        if let Some(g) = self.prescribed {
            let lifting = op.prescribed_displacement_load(g, &mut checkpoint)?;
            for (i, (r, b)) in rhs.iter_mut().zip(lifting).enumerate() {
                if i % 192 == 0 { poll(&mut checkpoint)?; }
                *r += b;
                if !r.is_finite() { return Err(GoalError3::Invalid("equilibrium RHS overflow")); }
            }
        }
        poll(&mut checkpoint)?; Ok(rhs)
    }
    pub(super) fn residuals(self, op: &AdaptiveElasticity3, u: &[f64], w: &[f64],
        checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<Vec<CellResidual3>, GoalError3> {
        let mut rows = op.cell_reference_residuals(u, w, self.external, &mut *checkpoint)?;
        if let Some(g) = self.prescribed {
            let lifting = op.prescribed_displacement_scale_work(g, w, &mut *checkpoint)?;
            for ((row, work), scale) in rows.iter_mut().zip(lifting).zip(op.scales()) {
                poll(checkpoint)?;
                row.load += scale * work;
                if !row.load.is_finite() || !row.residual().is_finite() {
                    return Err(GoalError3::Invalid("motion residual overflow"));
                }
            }
        }
        poll(checkpoint)?; Ok(rows)
    }
}

/// Reference volume/surface observation. Target and scale use response units;
/// functional densities are reintegrated on each grid, never nodally transferred.
#[derive(Clone, Copy)]
pub struct ResponseObservation3<'a> {
    pub functional: ReferenceLoad3<'a>,
    pub target: f64,
    pub scale: f64,
    pub weight: f64,
}
/// One complete loss evaluation and its aggregate adjoint RHS, without a solve.
#[derive(Debug, Clone)]
pub struct ResponseLinearization3 {
    pub value: f64,
    pub responses: Vec<f64>,
    pub adjoint_rhs: Vec<f64>,
}
/// Evaluate sum w/2 ((q^T u-target)/scale)^2 and its exact derivative in u.
/// This is also the RHS producer for the externally controlled fine adjoint.
/// A zero-weight observation is measured but never contributes arithmetic to
/// the objective. max_observations bounds the whole supplied observation slice.
pub fn linearize_responses3(op: &AdaptiveElasticity3, observations: &[ResponseObservation3<'_>],
    u: &[f64], max_observations: usize, mut checkpoint: impl FnMut() -> ControlFlow<()>)
    -> Result<ResponseLinearization3, GoalError3> {
    poll(&mut checkpoint)?;
    if observations.is_empty() || observations.len() > max_observations
        || u.len() != 3*op.nodes().len() || u.iter().any(|v| !v.is_finite())
        || u.iter().enumerate().any(|(i,v)| op.fixed()[i/3] && *v != 0.0)
        || observations.iter().any(|o| !o.target.is_finite() || !o.scale.is_finite() || o.scale <= 0.0
            || !o.weight.is_finite() || o.weight < 0.0) {
        return Err(GoalError3::Invalid("invalid response observations, field, or budget"));
    }
    let mut out = ResponseLinearization3 { value: 0.0, responses: Vec::with_capacity(observations.len()), adjoint_rhs: vec![0.0; u.len()] };
    for observation in observations {
        poll(&mut checkpoint)?;
        let q = op.reference_load(observation.functional, &mut checkpoint)?;
        let value = dot(&q,u);
        if !value.is_finite() { return Err(GoalError3::Invalid("response observation overflow")); }
        out.responses.push(value);
        if observation.weight == 0.0 { continue; }
        let residual = (value-observation.target)/observation.scale;
        let coefficient = observation.weight*residual/observation.scale;
        out.value += 0.5*observation.weight*residual*residual;
        if !coefficient.is_finite() || !out.value.is_finite() { return Err(GoalError3::Invalid("response loss overflow")); }
        for (i,(r,q)) in out.adjoint_rhs.iter_mut().zip(q).enumerate() {
            if i % 192 == 0 { poll(&mut checkpoint)?; }
            *r += coefficient*q;
            if !r.is_finite() { return Err(GoalError3::Invalid("response adjoint RHS overflow")); }
        }
    }
    poll(&mut checkpoint)?; Ok(out)
}
#[derive(Debug, Clone, Copy, Default)]
pub struct ResponseGoalCell3 {
    pub linear_dwr: f64,
    pub quadratic_remainder: f64,
    /// Absolute fine-cell DWR plus absolute fine-cell curvature contributions.
    /// Do not take an absolute value only after summing cancelling children.
    pub marking_mass: f64,
}
/// Actual quadratic objective difference, distinct from linearized goal work.
#[derive(Debug, Clone)]
pub struct ResponseGoalEstimate3 {
    /// Its values/goal_transfer are derivative work, NOT the quadratic objective.
    pub linearized: GoalEstimate3,
    pub coarse_value: f64,
    pub fine_value: f64,
    /// J_f(Pu_c)-J_c(u_c), including observation quadrature differences.
    pub objective_transfer: f64,
    /// Exact negative quadratic remainder about the fine solved field.
    pub quadratic_remainder: f64,
    pub identity_relative_defect: f64,
    pub coarse_responses: Vec<f64>,
    pub fine_responses: Vec<f64>,
    pub cells: BTreeMap<Octant3, ResponseGoalCell3>,
}
impl ResponseGoalEstimate3 {
    pub fn correction(&self) -> f64 {
        self.linearized.dwr + self.linearized.coarse_space + self.objective_transfer
            + self.quadratic_remainder + self.linearized.algebraic
    }
    pub fn mark(&self, theta: f64, max_marks: usize, checkpoint: impl FnMut() -> ControlFlow<()>)
        -> Result<GoalMarking3, GoalError3> {
        dorfler3(&self.cells.iter().map(|(&c,r)| (c,r.marking_mass)).collect(), theta, max_marks, checkpoint)
    }
}
/// Revalidate all four fields and reconstruct the actual fitting-loss change.
/// Coarse/fine adjoints must use their OWN aggregate loss derivatives, not a
/// compliance load or a coarse-only linearization. No physical solve occurs.
/// Fine scales must be exactly inherited; g and all observation laws must stay
/// unchanged between grids. Error beyond this enriched space, load quadrature
/// error, patch equivalence and continuum accuracy are NOT certified.
#[allow(clippy::too_many_arguments)]
pub fn estimate_response_goal3(transfer: &AdaptiveTransfer3<'_>, load: EquilibriumLoad3<'_>,
    observations: &[ResponseObservation3<'_>], fields: GoalFields3<'_>, options: GoalOptions3,
    max_observations: usize, mut checkpoint: impl FnMut() -> ControlFlow<()>)
    -> Result<ResponseGoalEstimate3, GoalError3> {
    admit_transfer3(transfer, options, &mut checkpoint)?;
    let coarse = linearize_responses3(transfer.coarse(), observations, fields.coarse_primal, max_observations, &mut checkpoint)?;
    let fine = linearize_responses3(transfer.fine(), observations, fields.fine_primal, max_observations, &mut checkpoint)?;
    let v = transfer.prolongate(fields.coarse_primal, &mut checkpoint)?;
    let projected = linearize_responses3(transfer.fine(), observations, &v, max_observations, &mut checkpoint)?;
    let linearized = estimate_vectors3(transfer, load, &coarse.adjoint_rhs, &fine.adjoint_rhs, fields, options, &mut checkpoint)?;
    let mut cells: BTreeMap<_,_> = linearized.cells.iter().map(|(&c,r)| (c,ResponseGoalCell3 {
        linear_dwr: r.signed_dwr, quadratic_remainder: 0.0, marking_mass: r.marking_mass,
    })).collect();
    let error: Vec<_> = fields.fine_primal.iter().zip(&v).map(|(u,v)| u-v).collect();
    let zero = vec![0.0; fields.fine_primal.len()];
    let mut remainder = 0.0;
    for (index, observation) in observations.iter().enumerate() {
        poll(&mut checkpoint)?;
        if observation.weight == 0.0 { continue; }
        let delta = (fine.responses[index]-projected.responses[index])/observation.scale;
        remainder -= 0.5*observation.weight*delta*delta;
        let coefficient = -0.5*observation.weight*delta/observation.scale;
        // At u=0 the same residual owner supplies precisely q_cell(error).
        let local = transfer.fine().cell_reference_residuals(&zero, &error, observation.functional, &mut checkpoint)?;
        for (row,&parent) in local.iter().zip(transfer.parents()) {
            poll(&mut checkpoint)?;
            let part = coefficient*row.load;
            let cell = cells.get_mut(&transfer.coarse().leaves()[parent]).expect("admitted transfer parent");
            cell.quadratic_remainder += part; cell.marking_mass += part.abs();
        }
    }
    let mut out = ResponseGoalEstimate3 { linearized, coarse_value: coarse.value, fine_value: fine.value,
        objective_transfer: projected.value-coarse.value, quadratic_remainder: remainder,
        identity_relative_defect: 0.0, coarse_responses: coarse.responses, fine_responses: fine.responses, cells };
    let values = [out.coarse_value,out.fine_value,out.linearized.dwr,out.linearized.coarse_space,
        out.objective_transfer,out.quadratic_remainder,out.linearized.algebraic];
    if values.iter().any(|v| !v.is_finite()) || !out.correction().is_finite()
        || out.cells.values().any(|r| !r.quadratic_remainder.is_finite() || !r.marking_mass.is_finite()) {
        return Err(GoalError3::Identity { relative_defect: f64::INFINITY });
    }
    let scale = values.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
    if scale > 0.0 {
        let sum: f64 = values[2..].iter().map(|v|v/scale).sum();
        let localized: f64 = out.cells.values().map(|r|r.quadratic_remainder/scale).sum();
        out.identity_relative_defect = ((out.fine_value/scale-out.coarse_value/scale)-sum).abs()
            .max((localized-out.quadratic_remainder/scale).abs());
    }
    if out.identity_relative_defect > options.identity_tolerance { return Err(GoalError3::Identity { relative_defect: out.identity_relative_defect }); }
    poll(&mut checkpoint)?; Ok(out)
}
