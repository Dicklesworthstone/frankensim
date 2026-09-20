//! Two-level, goal-weighted weak residuals on locally refined 3-D CutFEM.
//!
//! With v=P*u_c, t=P*z_c and A_f*z_f=q_f, the hierarchical term is
//! R_f(z_f-t). Its localization uses the actual bulk and ghost operators.
//! The coarse-space term R_f(t) is retained: different cut quadrature and
//! stabilization invalidate an automatic Galerkin-orthogonality assumption.
//! Goal transfer and finite primal/adjoint solve defects are also separate.
//! Their sum reconstructs J_f(u_f)-J_c(u_c), a TWO-LEVEL difference, not a
//! certified continuum error, a saturation theorem, or an optimality verdict.
use std::collections::BTreeMap;
use std::ops::ControlFlow;
use fs_cutfem::elastic3::ElasticityError3;
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octant3;

/// Externally solved homogeneous master fields, checked against actual operators.
#[derive(Clone, Copy)]
pub struct GoalFields3<'a> {
    /// Coarse primal for the supplied reference load.
    pub coarse_primal: &'a [f64],
    /// Enriched primal with the SAME inherited physical stiffness distribution.
    pub fine_primal: &'a [f64],
    /// Coarse adjoint for the supplied linear volume/surface goal.
    pub coarse_adjoint: &'a [f64],
    /// Enriched adjoint; computing it only on the coarse space is insufficient.
    pub fine_adjoint: &'a [f64],
}
impl<'a> GoalFields3<'a> {
    /// Compliance is self-adjoint: no separate adjoint solve is required.
    #[must_use]
    pub const fn compliance(coarse: &'a [f64], fine: &'a [f64]) -> Self {
        Self { coarse_primal: coarse, fine_primal: fine, coarse_adjoint: coarse, fine_adjoint: fine }
    }
}
/// Numerical admission gates; neither tolerance is a continuum-error bound.
#[derive(Debug, Clone, Copy)]
pub struct GoalOptions3 {
    /// Maximum recomputed relative Euclidean residual of each supplied field.
    pub residual_tolerance: f64,
    /// Maximum scaled defect of the complete two-level algebraic identity.
    pub identity_tolerance: f64,
}
impl Default for GoalOptions3 {
    fn default() -> Self { Self { residual_tolerance: 1e-8, identity_tolerance: 1e-8 } }
}
/// No partial estimate is returned on a refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum GoalError3 {
    /// Operator/field/callback refusal, including cancellation.
    Physics(ElasticityError3),
    /// Invalid controls, changed material distribution, or absent enrichment.
    Invalid(&'static str),
    /// A stale, partial or inaccurate field cannot supply a goal estimate.
    FieldResidual { field: &'static str, value: f64 },
    /// Nonfinite arithmetic or an inconsistent numerical identity.
    Identity { relative_defect: f64 },
}
impl From<ElasticityError3> for GoalError3 {
    fn from(e: ElasticityError3) -> Self { Self::Physics(e) }
}
impl std::fmt::Display for GoalError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { write!(f, "3-D goal estimate refused: {self:?}") }
}
impl std::error::Error for GoalError3 {}
fn poll(c: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), GoalError3> {
    if c().is_break() { Err(ElasticityError3::Cancelled.into()) } else { Ok(()) }
}
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(a, b)| a * b).sum() }

/// Hierarchical residual contributions grouped onto one original coarse leaf.
#[derive(Debug, Clone, Copy, Default)]
pub struct GoalCell3 {
    /// Signed sum; cancellations remain visible in the global correction.
    pub signed_dwr: f64,
    /// Sum of absolute FINE-cell contributions, preventing child cancellation
    /// from hiding a coarse region during refinement marking.
    pub marking_mass: f64,
}
/// Numerical two-level evidence. No field in this report certifies a bound.
#[derive(Debug, Clone)]
pub struct GoalEstimate3 {
    /// Goal value on the coarse solved field.
    pub coarse_value: f64,
    /// Goal value on the enriched solved field, with inherited physical scales.
    pub fine_value: f64,
    /// Signed hierarchical correction sum R_f(z_f-Pz_c).
    pub dwr: f64,
    /// R_f(Pz_c), including quadrature/stabilization inconsistency and coarse
    /// algebraic error. It is NOT silently assigned zero by orthogonality.
    pub coarse_space: f64,
    /// q_f^T P u_c - q_c^T u_c: the goal discretization/transfer difference.
    pub goal_transfer: f64,
    /// -z_f^T r_primal + r_adjoint^T(u_f-Pu_c), from finite enriched solves.
    pub algebraic: f64,
    /// Scaled discrepancy of the complete reconstructed goal difference.
    pub identity_relative_defect: f64,
    /// Coarse primal, fine primal, coarse adjoint, fine adjoint true residuals.
    pub field_residuals: [f64; 4],
    /// Deterministic coarse-leaf localization for marking.
    pub cells: BTreeMap<Octant3, GoalCell3>,
}
impl GoalEstimate3 {
    /// Signed two-level correction, not its absolute continuum-error bound.
    #[must_use]
    pub fn correction(&self) -> f64 { self.dwr + self.coarse_space + self.goal_transfer + self.algebraic }
    /// Dörfler marking from cancellation-safe absolute local masses.
    pub fn mark(&self, theta: f64, max_marks: usize, checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<GoalMarking3, GoalError3> {
        dorfler3(&self.cells.iter().map(|(&k, v)| (k, v.marking_mass)).collect(), theta, max_marks, checkpoint)
    }
}

/// Evaluate a linear VOLUME goal `integral goal_density dot u` on two grids.
/// This body-only entry point uses the same mixed-load estimator below.
/// For compliance pass body for both and use `GoalFields3::compliance`.
pub fn estimate_goal3(transfer: &AdaptiveTransfer3<'_>, body: &dyn Fn([f64; 3]) -> [f64; 3],
    goal_density: &dyn Fn([f64; 3]) -> [f64; 3], fields: GoalFields3<'_>, options: GoalOptions3,
    checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<GoalEstimate3, GoalError3> {
    estimate_reference_goal3(transfer, ReferenceLoad3::body(body), ReferenceLoad3::body(goal_density),
        fields, options, checkpoint)
}

/// Evaluate a linear volume/surface goal for a mixed reference load. Both laws
/// are independently integrated on each grid using its retained bulk and surface
/// rules; nodal loads are NOT prolonged. Pressure is inward-positive and fixed
/// in the reference configuration. Observation density may differ from loading,
/// requiring its own actual coarse and enriched adjoint fields.
///
/// No solves occur here: the caller owns their budgets. Every supplied field is
/// rechecked against its actual load and operator. Exact parent-inherited
/// PHYSICAL stiffness is required. The decomposition retains surface quadrature
/// and goal-transfer differences rather than assuming Galerkin orthogonality.
/// Missing surface rules refuse, never mean zero traction. Error beyond the
/// enriched space and the numerical surface/volume quadrature remain unbounded.
pub fn estimate_reference_goal3(transfer: &AdaptiveTransfer3<'_>, load: ReferenceLoad3<'_>,
    goal: ReferenceLoad3<'_>, fields: GoalFields3<'_>, options: GoalOptions3,
    mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<GoalEstimate3, GoalError3> {
    poll(&mut checkpoint)?;
    if ![options.residual_tolerance, options.identity_tolerance].iter().all(|v| v.is_finite() && *v > 0.0 && *v < 1.0) {
        return Err(GoalError3::Invalid("invalid numerical goal gates"));
    }
    let (coarse, fine) = (transfer.coarse(), transfer.fine());
    if !fine.leaves().iter().zip(transfer.parents()).any(|(f, &c)| f.level() > coarse.leaves()[c].level()) {
        return Err(GoalError3::Invalid("goal estimate requires an enriched active space"));
    }
    for (&scale, &parent) in fine.scales().iter().zip(transfer.parents()) {
        if scale.to_bits() != coarse.scales()[parent].to_bits() {
            return Err(GoalError3::Invalid("enriched physical stiffness is not parent-inherited"));
        }
    }
    let bc = coarse.reference_load(load, &mut checkpoint)?;
    let bf = fine.reference_load(load, &mut checkpoint)?;
    let qc = coarse.reference_load(goal, &mut checkpoint)?;
    let qf = fine.reference_load(goal, &mut checkpoint)?;
    let mut residuals = [0.0; 4];
    for (i, (name, op, x, b)) in [
        ("coarse-primal", coarse, fields.coarse_primal, &bc),
        ("fine-primal", fine, fields.fine_primal, &bf),
        ("coarse-adjoint", coarse, fields.coarse_adjoint, &qc),
        ("fine-adjoint", fine, fields.fine_adjoint, &qf),
    ].into_iter().enumerate() {
        residuals[i] = op.field_residual(x, b, &mut checkpoint)?;
        if residuals[i] > options.residual_tolerance {
            return Err(GoalError3::FieldResidual { field: name, value: residuals[i] });
        }
    }
    let v = transfer.prolongate(fields.coarse_primal, &mut checkpoint)?;
    let t = transfer.prolongate(fields.coarse_adjoint, &mut checkpoint)?;
    let weight: Vec<_> = fields.fine_adjoint.iter().zip(&t).map(|(z, t)| z - t).collect();
    let error: Vec<_> = fields.fine_primal.iter().zip(&v).map(|(u, v)| u - v).collect();
    let local = fine.cell_reference_residuals(&v, &weight, load, &mut checkpoint)?;
    let consistency = fine.cell_reference_residuals(&v, &t, load, &mut checkpoint)?;
    let primal_defect = fine.cell_reference_residuals(fields.fine_primal, fields.fine_adjoint, load, &mut checkpoint)?;
    let adjoint_defect = fine.cell_reference_residuals(fields.fine_adjoint, &error, goal, &mut checkpoint)?;
    let mut cells: BTreeMap<_, GoalCell3> = coarse.leaves().iter().map(|&c| (c, GoalCell3::default())).collect();
    for (term, &parent) in local.iter().zip(transfer.parents()) {
        poll(&mut checkpoint)?;
        let cell = cells.get_mut(&coarse.leaves()[parent]).expect("admitted parent");
        cell.signed_dwr += term.residual();
        cell.marking_mass += term.residual().abs();
    }
    let mut report = GoalEstimate3 {
        coarse_value: dot(&qc, fields.coarse_primal), fine_value: dot(&qf, fields.fine_primal),
        dwr: local.iter().map(|t| t.residual()).sum(),
        coarse_space: consistency.iter().map(|t| t.residual()).sum(),
        goal_transfer: dot(&qf, &v) - dot(&qc, fields.coarse_primal),
        algebraic: -primal_defect.iter().map(|t| t.residual()).sum::<f64>() + adjoint_defect.iter().map(|t| t.residual()).sum::<f64>(),
        identity_relative_defect: 0.0, field_residuals: residuals, cells,
    };
    let numbers = [report.coarse_value, report.fine_value, report.dwr, report.coarse_space, report.goal_transfer, report.algebraic];
    if !numbers.iter().all(|v| v.is_finite()) || report.cells.values().any(|v| !v.signed_dwr.is_finite() || !v.marking_mass.is_finite()) {
        return Err(GoalError3::Identity { relative_defect: f64::INFINITY });
    }
    let scale = numbers.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if scale > 0.0 {
        report.identity_relative_defect = ((report.fine_value / scale - report.coarse_value / scale)
            - (report.dwr / scale + report.coarse_space / scale + report.goal_transfer / scale + report.algebraic / scale)).abs();
    }
    if !report.correction().is_finite() || report.identity_relative_defect > options.identity_tolerance {
        return Err(GoalError3::Identity { relative_defect: report.identity_relative_defect });
    }
    poll(&mut checkpoint)?; Ok(report)
}

/// Bounded marking result; zero signal and an insufficient mark budget are not convergence.
#[derive(Debug, Clone)]
pub struct GoalMarking3 {
    /// Stable descending-mass/ascending-key prefix.
    pub marked: Vec<Octant3>,
    /// Fraction of total numerical mass carried by the marked prefix.
    pub achieved_fraction: f64,
    /// False for zero signal or a cap that prevents reaching theta.
    pub target_met: bool,
}
/// Same Dörfler ordering as the 2-D path, with explicit admission and a mark cap.
/// Normalize before summing so finite extreme-scale indicators do not overflow.
pub fn dorfler3(masses: &BTreeMap<Octant3, f64>, theta: f64, max_marks: usize,
    mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<GoalMarking3, GoalError3> {
    poll(&mut checkpoint)?;
    if !theta.is_finite() || theta <= 0.0 || theta > 1.0
        || masses.values().any(|v| !v.is_finite() || *v < 0.0) { return Err(GoalError3::Invalid("invalid marking policy/mass")); }
    let scale = masses.values().copied().fold(0.0_f64, f64::max);
    let mut result = GoalMarking3 { marked: Vec::new(), achieved_fraction: 0.0, target_met: false };
    if scale == 0.0 { poll(&mut checkpoint)?; return Ok(result); }
    let total: f64 = masses.values().map(|v| v / scale).sum();
    let mut selected = 0.0;
    for (key, mass) in crate::mark::indicator_order(masses) {
        poll(&mut checkpoint)?;
        if result.marked.len() == max_marks || selected / total >= theta || mass == 0.0 { break; }
        selected += mass / scale; result.marked.push(key);
    }
    result.achieved_fraction = (selected / total).min(1.0);
    result.target_met = result.achieved_fraction >= theta;
    poll(&mut checkpoint)?; Ok(result)
}
