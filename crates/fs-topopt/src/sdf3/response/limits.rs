//! Explicit response inequalities on the same cut-cell equilibrium operator.
//!
//! Each signed residual c_i is nonpositive when its physical limit holds. The
//! scalar constraint is the UNNORMALIZED log-sum-exp upper approximation
//! G = max(c) + log(sum(exp(beta*(c_i-max(c)))))/beta. In exact arithmetic,
//! max(c) <= G <= max(c)+log(m)/beta. Do not subtract log(m): the resulting
//! averaged aggregate could hide a violated limit. This is a numerical design
//! constraint, not a continuum guarantee; individual residuals remain visible.
//!
//! One current-density preparation serves all primal and aggregate-adjoint
//! solves. Reaction state AND direct material derivatives, prescribed-motion
//! lifting, ghost/Nitsche stiffness and the original density pullback all enter.
use super::*;
use fs_cutfem::elastic3::adaptive::dirichlet::PrescribedMotion3;

mod physical;
pub mod binary;
pub mod design;
pub mod refinement;
pub mod strength;

/// Physical scalar whose value is bounded. Observation laws remain fixed in
/// the reference configuration and must use the current grid's DOF ordering.
#[derive(Clone, Copy)]
pub enum BoundedResponse3<'a> {
    /// q^T u; assemble q using the existing reference body/surface integrators.
    Displacement(&'a [f64]),
    /// Embedded numerical support traction tested by a force or moment mode.
    /// This is not an extra strongly clamped box-node reaction or actuator work.
    Reaction(&'a PrescribedMotion3<'a>),
}
/// Signed inequality direction; no automatic absolute value or sign reversal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResponseLimitSense3 { AtMost, AtLeast }
impl ResponseLimitSense3 {
    fn sign(self) -> f64 { match self { Self::AtMost => 1.0, Self::AtLeast => -1.0 } }
}
/// One upper or lower response limit. To constrain a window, supply both sides.
/// Scale has the response's physical units and is independent of its bound;
/// zero and negative bounds are valid. There are no averaging/priority weights.
#[derive(Clone, Copy)]
pub struct ResponseLimit3<'a> {
    pub response: BoundedResponse3<'a>,
    pub sense: ResponseLimitSense3,
    pub bound: f64,
    pub scale: f64,
}
/// One independent fixed-load/prescribed-motion experiment and its inequalities.
#[derive(Clone, Copy)]
pub struct ResponseLimitCase3<'a> {
    pub force: &'a [f64],
    pub prescribed: Option<&'a PrescribedMotion3<'a>>,
    pub limits: &'a [ResponseLimit3<'a>],
}
#[derive(Debug, Clone, Copy)]
pub struct ResponseLimitOptions3 {
    pub max_cases: usize,
    /// Total scalar limits, including both sides of any window.
    pub max_limits: usize,
    /// Positive dimensionless inverse smoothing width. Larger values reduce
    /// conservatism but can make the optimization less well conditioned.
    pub sharpness: f64,
}
impl Default for ResponseLimitOptions3 {
    fn default() -> Self { Self { max_cases: 32, max_limits: 256, sharpness: 40.0 } }
}
/// Complete family at ONE raw design. Fields/adjoints differentiate the signed
/// smooth constraint, not the squared fitting objective or external-load work.
#[derive(Debug, Clone)]
pub struct ResponseLimitEvaluation3 {
    pub rho: Vec<f64>,
    pub projected_rho: Vec<f64>,
    pub volume_fraction: f64,
    pub volume_gradient: Vec<f64>,
    /// Physical material inventory, populated by a mass-objective design.
    /// Ordinary response-only evaluation does not request this extra pullback.
    pub mass_kg: Option<f64>,
    /// Exact raw-design derivative of `mass_kg`, including filter/projection.
    pub mass_gradient: Option<Vec<f64>>,
    pub constraint: f64,
    pub constraint_gradient: Vec<f64>,
    /// Maximum ORIGINAL signed normalized residual, without smoothing bias.
    pub worst_residual: f64,
    /// Same experiment/limit ordering as the input, no discarded zero weights.
    pub measurements: Vec<Vec<f64>>,
    pub residuals: Vec<Vec<f64>>,
    pub displacements: Vec<Vec<f64>>,
    pub adjoints: Vec<Vec<f64>>,
    pub work: SolveWork,
    pub(super) scales: Vec<f64>,
}
impl ResponseLimitEvaluation3 {
    /// Both the aggregate and original numerical limits must pass. A caller's
    /// tolerance is dimensionless and does not confer continuum certification.
    #[must_use]
    pub fn feasible(&self, tolerance: f64) -> bool {
        tolerance.is_finite() && tolerance >= 0.0 && self.constraint.is_finite()
            && self.worst_residual.is_finite() && self.constraint <= tolerance
            && self.worst_residual <= tolerance
    }
}
fn admit_limits<O: AdaptiveSdf3Elasticity>(study: &CutDensityStudy3<O>,
    cases: &[ResponseLimitCase3<'_>], options: ResponseLimitOptions3) -> Result<(), ResponseError3> {
    let invalid = ResponseError3::Invalid;
    if cases.is_empty() || cases.len() > options.max_cases || !options.sharpness.is_finite()
        || options.sharpness <= 0.0 || !(1.0/options.sharpness).is_finite() {
        return Err(invalid("invalid response limit family or smoothing policy"));
    }
    let n = study.operator.n(); let mut count = 0usize;
    for case in cases {
        count = count.checked_add(case.limits.len()).filter(|c| *c <= options.max_limits)
            .ok_or(invalid("response limit allowance exhausted"))?;
        if case.limits.is_empty() || case.force.len() != n || case.force.iter().any(|v| !v.is_finite()) {
            return Err(invalid("empty limits or invalid response-limit force"));
        }
        if case.prescribed.is_some() && study.operator.adaptive().embedded_dirichlet_penalty().is_none() {
            return Err(invalid("prescribed motion requires embedded supports"));
        }
        for limit in case.limits {
            if !limit.bound.is_finite() || !limit.scale.is_finite() || limit.scale <= 0.0
                || !(1.0/limit.scale).is_finite() {
                return Err(invalid("response bounds/scales must be finite with positive invertible scales"));
            }
            match limit.response {
                BoundedResponse3::Displacement(q) if q.len() != n || q.iter().any(|v| !v.is_finite()) =>
                    return Err(invalid("invalid bounded displacement observation")),
                BoundedResponse3::Reaction(_) if study.operator.adaptive().embedded_dirichlet_penalty().is_none() =>
                    return Err(invalid("bounded reaction requires embedded supports")),
                _ => {},
            }
        }
    }
    Ok(())
}

fn aggregate(residuals: &[f64], beta: f64) -> Result<(f64, f64, Vec<f64>), ResponseError3> {
    if residuals.is_empty() || residuals.iter().any(|v| !v.is_finite())
        || !beta.is_finite() || beta <= 0.0 {
        return Err(ResponseError3::Invalid("invalid smooth-maximum input"));
    }
    let max = residuals.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut weights: Vec<f64> = residuals.iter().map(|c| {
        let x = beta*(c-max);
        // Negative infinite exponents mean negligible terms, not overflow of
        // the final finite log-sum. The maximal term is exactly exp(0)=1.
        if x == f64::NEG_INFINITY { 0.0 } else { fs_math::det::exp(x) }
    }).collect();
    let sum = finite(weights.iter().sum())?;
    if sum < 1.0 || weights.iter().any(|w| !w.is_finite() || *w < 0.0) {
        return Err(ResponseError3::Invalid("invalid smooth-maximum weights"));
    }
    let value = finite(max+fs_math::det::ln(sum)/beta)?;
    for w in &mut weights { *w /= sum; }
    Ok((value, max, weights))
}

impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Evaluate all declared limits using real independent primal solves and
    /// one aggregate adjoint per experiment. The mathematical upper aggregate
    /// is conservative, but individual finite-arithmetic residuals are ALSO
    /// retained and checked by `feasible`; no certificate is inferred from it.
    ///
    /// Geometry, force vectors, observations and prescribed motion must be
    /// density independent. Reaction functionals themselves retain their exact
    /// density dependence. The existing filter/projection honors protected
    /// physical regions. No matrix of per-limit density gradients is stored.
    /// Reaction traces are reintegrated in the adjoint pass to bound storage.
    /// Every solve/trace shares `control`. BOTH success and returned errors
    /// restore incoming scales. Pure callback panics are not caught.
    pub fn evaluate_response_limits(&mut self, rho: &[f64], cases: &[ResponseLimitCase3<'_>],
        options: ResponseLimitOptions3, control: &mut SolveControl<'_>)
        -> Result<ResponseLimitEvaluation3, ResponseError3> {
        self.response_limits_with_mass(rho, cases, options, false, control)
    }

    // The minimum-material owner requests physical inventory from the SAME
    // prepared design. No second forward filter or elasticity solve is needed.
    fn response_limits_with_mass(&mut self, rho: &[f64], cases: &[ResponseLimitCase3<'_>],
        options: ResponseLimitOptions3, include_mass: bool, control: &mut SolveControl<'_>)
        -> Result<ResponseLimitEvaluation3, ResponseError3> {
        control.checkpoint("response-limits-start")?;
        admit_limits(self, cases, options)?;
        if rho.len() != self.cells() || rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r)) {
            return Err(ResponseError3::Invalid("one raw response-limit density in [0,1] per cell required"));
        }
        let design = self.design(rho, control)?;
        let mass = if include_mass { Some(self.mass_prepared(&design, control)?) } else { None };
        let previous = self.operator.scales().to_vec();
        self.operator.set_scales(&design.scales)?;
        let result = (|| {
            let prepared = self.operator.prepare_elasticity(control)?;
            let mut evaluation = self.response_limits_prepared(rho, cases, options, &design, &prepared, control)?;
            if let Some((value, gradient)) = mass {
                evaluation.mass_kg = Some(value);
                evaluation.mass_gradient = Some(gradient);
            }
            Ok(evaluation)
        })();
        self.operator.set_scales(&previous).expect("previous admitted scales remain valid");
        result
    }

    // Shared with the joint response/stress owner. The caller has admitted the
    // complete family and owns the current-density transaction. An immutable
    // borrow keeps this material and its prepared action fixed through the
    // component solves; there is no public cached-field admission shortcut.
    #[allow(clippy::too_many_arguments)]
    fn response_limits_prepared(&self, rho: &[f64], cases: &[ResponseLimitCase3<'_>],
        options: ResponseLimitOptions3, design: &crate::sdf3::Design3,
        prepared: &impl fs_sparse::precond::Precond, control: &mut SolveControl<'_>)
        -> Result<ResponseLimitEvaluation3, ResponseError3> {
            let op = self.operator.adaptive();
            let physical::LimitPrimals3 { displacements, measurements, residuals, all } =
                self.response_limit_primals(cases, prepared, control)?;
            control.checkpoint("response-limits-aggregate")?;
            let (constraint, worst_residual, weights) = aggregate(&all, options.sharpness)?;
            let mut local = vec![0.0; self.cells()];
            let mut adjoints = Vec::with_capacity(cases.len());
            let mut offset = 0;
            for (case, u) in cases.iter().zip(&displacements) {
                let mut rhs = vec![0.0; op.n()];
                for limit in case.limits {
                    control.checkpoint("response-limits-pullback")?;
                    let coefficient = finite(weights[offset]*limit.sense.sign()/limit.scale)?;
                    offset += 1;
                    if coefficient == 0.0 { continue; }
                    match limit.response {
                        BoundedResponse3::Displacement(q) => {
                            for (i, (r, q)) in rhs.iter_mut().zip(q).enumerate() {
                                if !op.fixed()[i/3] { *r = finite(*r+coefficient*q)?; }
                            }
                        }
                        BoundedResponse3::Reaction(mode) => {
                            let reaction = op.embedded_reaction(u, case.prescribed, mode, || poll(control))?;
                            for (r, q) in rhs.iter_mut().zip(reaction.displacement_gradient) {
                                *r = finite(*r+coefficient*q)?;
                            }
                            for (r, direct) in local.iter_mut().zip(reaction.scale_gradient) {
                                *r = finite(*r+coefficient*direct)?;
                            }
                        }
                    }
                }
                let z = checked_solve_preconditioned(&self.operator, prepared, &rhs, 1e-11,
                    "response-limits-adjoint", control)?;
                let stiffness = op.scale_bilinear_forms(&z, u, || poll(control))?;
                let lifting = match case.prescribed {
                    Some(g) => op.prescribed_displacement_scale_work(g, &z, || poll(control))?,
                    None => vec![0.0; self.cells()],
                };
                for ((l, k), b) in local.iter_mut().zip(stiffness).zip(lifting) { *l = finite(*l+b-k)?; }
                adjoints.push(z);
            }
            let p = self.params;
            for ((l, r), slope) in local.iter_mut().zip(&design.projected).zip(&design.slope) {
                let power = if p.penal == 1.0 { 1.0 } else if *r == 0.0 { 0.0 }
                    else { fs_math::det::pow(*r, p.penal-1.0) };
                *l = finite(*l*(1.0-p.e_min)*p.penal*power*slope)?;
            }
            let constraint_gradient = self.pullback(&local, control)?;
            let dv: Vec<_> = self.mass.iter().zip(&design.slope).map(|(m,s)| m*s).collect();
            let volume_gradient = self.pullback(&dv, control)?;
            control.checkpoint("response-limits-publish")?;
            Ok(ResponseLimitEvaluation3 { rho: rho.to_vec(), projected_rho: design.projected.clone(),
                scales: design.scales.clone(), volume_fraction: design.volume, volume_gradient,
                mass_kg: None, mass_gradient: None,
                constraint, constraint_gradient, worst_residual, measurements, residuals,
                displacements, adjoints, work: control.work() })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g0_upper_aggregate_cannot_average_away_a_violated_limit() {
        let (g, m, w) = aggregate(&[0.1, -100.0, -100.0], 40.0).unwrap();
        assert!(g >= m && m == 0.1); assert_eq!(w[0], 1.0);
        let (g, m, w) = aggregate(&[0.0; 5], 20.0).unwrap();
        assert!(g > m); assert!((g-fs_math::det::ln(5.0)/20.0).abs() < 1e-15);
        assert!(w.iter().all(|w| (*w-0.2).abs() < 1e-15));
        assert_eq!(aggregate(&[-0.7], 4.0).unwrap(), (-0.7, -0.7, vec![1.0]));
    }
    #[test]
    fn g0_shifted_aggregation_handles_wide_ranges_and_refuses_invalid_data() {
        let a = aggregate(&[1e200, -1e200], 1e200).unwrap();
        assert_eq!(a.0, 1e200); assert_eq!(a.2, vec![1.0, 0.0]);
        assert!(aggregate(&[], 1.0).is_err());
        assert!(aggregate(&[f64::NAN], 1.0).is_err());
        assert!(aggregate(&[1.0], 0.0).is_err());
    }
}
