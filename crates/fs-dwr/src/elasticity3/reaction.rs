//! Two-grid goals containing the actual embedded-support reaction.
//! A reaction is AFFINE in displacement: R=q(s)^T u+c(s,g). Both q and c
//! change with the grid's retained Nitsche rule, even with inherited scales.
//! Never discard c_f-c_c or use the imposed primal motion as adjoint data.
use super::*;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;

/// A common reference observation, optionally plus a support-force/moment mode.
/// The reaction mode and imposed motion must be density-independent pure laws.
/// Reaction sign and support selection are exactly those of embedded_reaction.
#[derive(Clone, Copy)]
pub struct AffineMotionGoal3<'a> {
    pub displacement: ReferenceLoad3<'a>,
    pub reaction_mode: Option<&'a PrescribedMotion3<'a>>,
}
/// Independently assembled on ONE grid at its current physical stiffness.
/// Reassemble after changing scales, geometry or prescribed motion.
#[derive(Debug, Clone)]
pub struct AffineGoalForm3 {
    pub gradient: Vec<f64>,
    pub offset: f64,
}
/// Integrate q and c using existing load and reaction owners, without a solve.
/// The reaction is evaluated at zero displacement to retain its offset without
/// subtracting two large evaluated reaction/gradient contractions.
pub fn assemble_affine_motion_goal3(op: &AdaptiveElasticity3, goal: AffineMotionGoal3<'_>,
    prescribed: Option<&PrescribedMotion3<'_>>, mut checkpoint: impl FnMut() -> ControlFlow<()>)
    -> Result<AffineGoalForm3, GoalError3> {
    poll(&mut checkpoint)?;
    let mut gradient = op.reference_load(goal.displacement, &mut checkpoint)?;
    let mut offset = 0.0;
    if let Some(mode) = goal.reaction_mode {
        let zero = vec![0.0; gradient.len()];
        let reaction = op.embedded_reaction(&zero, prescribed, mode, &mut checkpoint)?;
        offset = reaction.value;
        for (i, (q, r)) in gradient.iter_mut().zip(reaction.displacement_gradient).enumerate() {
            if i % 192 == 0 { poll(&mut checkpoint)?; }
            *q += r;
            if !q.is_finite() { return Err(GoalError3::Invalid("affine reaction gradient overflow")); }
        }
    }
    if !offset.is_finite() { return Err(GoalError3::Invalid("affine reaction offset overflow")); }
    poll(&mut checkpoint)?;
    Ok(AffineGoalForm3 { gradient, offset })
}

/// Actual affine goal difference with the boundary offset reported separately.
/// `linearized` retains the existing residual, transfer and algebraic terms.
/// Its goal values are derivative work, not the affine reaction values below.
#[derive(Debug, Clone)]
pub struct AffineGoalEstimate3 {
    pub linearized: GoalEstimate3,
    pub coarse_value: f64,
    pub fine_value: f64,
    pub coarse_offset: f64,
    pub fine_offset: f64,
    pub offset_transfer: f64,
    pub identity_relative_defect: f64,
}
impl AffineGoalEstimate3 {
    pub fn correction(&self) -> f64 { self.linearized.correction() + self.offset_transfer }
    /// Mark from real hierarchical residual contributions. The separate offset
    /// change is NOT fabricated into an adjoint cell indicator. Empty marks in
    /// the presence of a transfer discrepancy do not establish accuracy.
    pub fn mark(&self, theta: f64, max_marks: usize, checkpoint: impl FnMut() -> ControlFlow<()>)
        -> Result<GoalMarking3, GoalError3> {
        self.linearized.mark(theta, max_marks, checkpoint)
    }
}

/// Revalidate all four fields and reconstruct the full two-grid affine goal.
/// Primal fields use external loading plus imposed motion; adjoints use ONLY
/// the assembled goal gradient, with homogeneous boundary motion. No new
/// residual localization, solver, quadrature, or continuum certificate.
/// Exact parent-inherited physical scales and genuine enrichment are required.
pub fn estimate_affine_motion_goal3(transfer: &AdaptiveTransfer3<'_>, load: EquilibriumLoad3<'_>,
    goal: AffineMotionGoal3<'_>, fields: GoalFields3<'_>, options: GoalOptions3,
    mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<AffineGoalEstimate3, GoalError3> {
    admit_transfer3(transfer, options, &mut checkpoint)?;
    let coarse = assemble_affine_motion_goal3(transfer.coarse(), goal, load.prescribed, &mut checkpoint)?;
    let fine = assemble_affine_motion_goal3(transfer.fine(), goal, load.prescribed, &mut checkpoint)?;
    let linearized = estimate_vectors3(transfer, load, &coarse.gradient, &fine.gradient,
        fields, options, &mut checkpoint)?;
    let mut report = AffineGoalEstimate3 {
        coarse_value: linearized.coarse_value + coarse.offset,
        fine_value: linearized.fine_value + fine.offset,
        coarse_offset: coarse.offset, fine_offset: fine.offset,
        offset_transfer: fine.offset - coarse.offset,
        linearized, identity_relative_defect: 0.0,
    };
    let values = [report.coarse_value, report.fine_value, report.coarse_offset,
        report.fine_offset, report.offset_transfer, report.linearized.correction(), report.correction()];
    if values.iter().any(|v| !v.is_finite()) {
        return Err(GoalError3::Invalid("affine goal difference overflow"));
    }
    // The offset and derivative work can be individually large and cancel.
    // Scale by the actual identity terms, not a claimed reaction-error bound.
    let scale = values.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if scale > 0.0 {
        report.identity_relative_defect = ((report.fine_value/scale - report.coarse_value/scale)
            - (report.linearized.correction()/scale + report.offset_transfer/scale)).abs();
    }
    if !report.identity_relative_defect.is_finite() || report.identity_relative_defect > options.identity_tolerance {
        return Err(GoalError3::Identity { relative_defect: report.identity_relative_defect });
    }
    poll(&mut checkpoint)?;
    Ok(report)
}
