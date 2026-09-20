//! Displacement-response objectives with density-dependent prescribed motion.
//!
//! K(s) u = f + b_g(s), J = alpha V + sum w/2 ((q^T u-target)/scale)^2.
//! One adjoint per load case solves K z = dJ/du. Its cell-scale derivative is
//! z^T db_g/ds - z^T dK/ds u, followed by the ORIGINAL SIMP/projection/filter
//! pullback. Neither differentiating CG nor treating b_g as fixed is correct.
//! These are numerical discrete responses, not reaction/actuator work or bounds.
use super::*;
use crate::SolveWork;
use fs_cutfem::elastic3::ElasticityError3;
use std::ops::ControlFlow;

/// One fixed linear observation of independent displacement coordinates.
/// Assemble `q` with the operator's body/reference-load integrators to observe
/// a volume/surface integral. The coefficients and target must describe the
/// current grid; no nearest-node or automatic inter-grid observation transfer.
#[derive(Clone, Copy)]
pub struct ResponseTarget3<'a> {
    pub q: &'a [f64],
    pub target: f64,
    /// Positive physical response scale. Misfit is (q^T u-target)/scale.
    pub scale: f64,
    /// Nonnegative weight; no automatic normalization.
    pub weight: f64,
}
/// Independent experiment on one fixed reference geometry and support patch.
#[derive(Clone, Copy)]
pub struct ResponseCase3<'a> {
    /// Fixed external body/traction load, NOT a prescribed-motion lifting.
    pub force: &'a [f64],
    /// Pure reference displacement law on the operator's embedded support.
    /// None requests homogeneous supports and works with legacy box clamps.
    pub prescribed: Option<&'a dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
    pub targets: &'a [ResponseTarget3<'a>],
}
#[derive(Debug, Clone, Copy)]
pub struct ResponseOptions3 {
    pub max_cases: usize,
    /// Total observations across all cases, not a per-case allowance.
    pub max_observations: usize,
    /// Nonnegative coefficient of numerical projected volume in the objective.
    pub volume_weight: f64,
}
impl Default for ResponseOptions3 {
    fn default() -> Self { Self { max_cases: 32, max_observations: 256, volume_weight: 0.0 } }
}
#[derive(Debug, Clone, PartialEq)]
pub enum ResponseError3 {
    Invalid(&'static str),
    Evaluation(EvaluationStop),
    Physics(ElasticityError3),
}
impl From<EvaluationStop> for ResponseError3 {
    fn from(e: EvaluationStop) -> Self { Self::Evaluation(e) }
}
impl From<ElasticityError3> for ResponseError3 {
    fn from(e: ElasticityError3) -> Self {
        if matches!(e, ElasticityError3::Cancelled) { Self::Evaluation(EvaluationStop::Cancelled) }
        else { Self::Physics(e) }
    }
}
impl std::fmt::Display for ResponseError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { write!(f, "response design refused: {self:?}") }
}
impl std::error::Error for ResponseError3 {}

/// Complete experiment family at one explicitly retained design. Adjoint
/// fields differentiate the summed response misfit, not augmented-RHS work.
#[derive(Debug, Clone)]
pub struct ResponseEvaluation3 {
    pub rho: Vec<f64>,
    pub projected_rho: Vec<f64>,
    pub objective: f64,
    pub gradient: Vec<f64>,
    pub volume_fraction: f64,
    pub volume_gradient: Vec<f64>,
    pub responses: Vec<Vec<f64>>,
    pub displacements: Vec<Vec<f64>>,
    pub adjoints: Vec<Vec<f64>>,
    /// f^T u, separate from both the fitting objective and (f+b_g)^T u.
    pub external_work: Vec<f64>,
    /// Includes variational lifting work; NOT physical actuator work.
    pub augmented_rhs_work: Vec<f64>,
    pub work: SolveWork,
    // Numeric state for the owning optimizer, never supplied by a caller.
    pub(super) scales: Vec<f64>,
}
fn poll(control: &mut SolveControl<'_>) -> ControlFlow<()> {
    if control.checkpoint("sdf3-response-physics").is_ok() { ControlFlow::Continue(()) }
    else { ControlFlow::Break(()) }
}
fn finite(x: f64) -> Result<f64, ResponseError3> {
    if x.is_finite() { Ok(x) } else { Err(ResponseError3::Invalid("nonfinite response arithmetic")) }
}
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)| a*b).sum() }

pub(super) fn admit<O: AdaptiveSdf3Elasticity>(study: &CutDensityStudy3<O>, cases: &[ResponseCase3<'_>],
    options: ResponseOptions3) -> Result<(), ResponseError3> {
    let invalid = ResponseError3::Invalid;
    if cases.is_empty() || cases.len() > options.max_cases || !options.volume_weight.is_finite() || options.volume_weight < 0.0 {
        return Err(invalid("invalid response family, weight or case budget"));
    }
    let n = study.operator.n(); let mut count = 0usize; let mut positive = options.volume_weight > 0.0;
    for case in cases {
        count = count.checked_add(case.targets.len()).filter(|v| *v <= options.max_observations)
            .ok_or(invalid("observation budget exhausted"))?;
        if case.force.len() != n || case.force.iter().any(|x| !x.is_finite()) || case.targets.is_empty() {
            return Err(invalid("invalid external load or empty observation family"));
        }
        if case.prescribed.is_some() && study.operator.adaptive().embedded_dirichlet_penalty().is_none() {
            return Err(invalid("prescribed motion requires embedded Dirichlet support"));
        }
        for target in case.targets {
            if target.q.len() != n || target.q.iter().any(|x| !x.is_finite()) || !target.target.is_finite()
                || !target.scale.is_finite() || target.scale <= 0.0 || !target.weight.is_finite() || target.weight < 0.0 {
                return Err(invalid("invalid observation shape, target, scale or weight"));
            }
            positive |= target.weight > 0.0;
        }
    }
    if !positive { return Err(invalid("objective requires a positive observation or volume weight")); }
    Ok(())
}
impl<O: AdaptiveSdf3Elasticity> CutDensityStudy3<O> {
    /// Evaluate response fitting, with one shared preparation and one primal/
    /// aggregate-adjoint pair per independent case. All physics/filter/setup
    /// work consumes the existing control, including work discarded on failure.
    ///
    /// This evaluation is observational: the incoming operator scales are
    /// restored on BOTH success and every returned error. Returned fields carry
    /// their own raw/projected design. The response optimizer installs only its
    /// accepted evaluation; a caller cannot mistake the last rejected trial for
    /// the study's material state. Callback panics are not caught here.
    ///
    /// q and f must be density-independent. g may vary between cases but not
    /// with density. Geometry/support/load rebuilding across grids is explicit.
    pub fn evaluate_responses(&mut self, rho: &[f64], cases: &[ResponseCase3<'_>], options: ResponseOptions3,
        control: &mut SolveControl<'_>) -> Result<ResponseEvaluation3, ResponseError3> {
        control.checkpoint("sdf3-response-start")?;
        admit(self, cases, options)?;
        if rho.len() != self.cells() || rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r)) {
            return Err(ResponseError3::Invalid("one raw density in [0,1] per cell required"));
        }
        let design = self.design(rho, control)?;
        let previous = self.operator.scales().to_vec();
        self.operator.set_scales(&design.scales)?;
        let result = (|| {
            let prepared = self.operator.prepare_elasticity(control)?;
            let op = self.operator.adaptive();
            let mut objective = 0.0;
            let mut local = vec![0.0; self.cells()];
            let mut responses = Vec::with_capacity(cases.len());
            let mut displacements = Vec::with_capacity(cases.len());
            let mut adjoints = Vec::with_capacity(cases.len());
            let mut external_work = Vec::with_capacity(cases.len());
            let mut augmented_rhs_work = Vec::with_capacity(cases.len());
            for case in cases {
                control.checkpoint("sdf3-response-case")?;
                let external: Vec<_> = case.force.iter().enumerate().map(|(i,&v)| if op.fixed()[i/3] { 0.0 } else { v }).collect();
                let mut rhs = external.clone();
                if let Some(g) = case.prescribed {
                    let lifting = op.prescribed_displacement_load(g, || poll(control))?;
                    for (r, b) in rhs.iter_mut().zip(lifting) { *r = finite(*r+b)?; }
                }
                let u = checked_solve_preconditioned(&self.operator, &prepared, &rhs, 1e-11, "sdf3-response-primal", control)?;
                let mut observed = Vec::with_capacity(case.targets.len());
                let mut adjoint_rhs = vec![0.0; op.n()];
                for target in case.targets {
                    control.checkpoint("sdf3-response-observation")?;
                    let value = finite(dot(target.q, &u))?;
                    observed.push(value);
                    if target.weight == 0.0 { continue; }
                    let error = finite((value-target.target)/target.scale)?;
                    objective = finite(objective + 0.5*target.weight*error*error)?;
                    let coefficient = finite((target.weight*error)/target.scale)?;
                    for (i, (b, &q)) in adjoint_rhs.iter_mut().zip(target.q).enumerate() {
                        if !op.fixed()[i/3] { *b = finite(*b+coefficient*q)?; }
                    }
                }
                let z = checked_solve_preconditioned(&self.operator, &prepared, &adjoint_rhs, 1e-11, "sdf3-response-adjoint", control)?;
                let stiffness = op.scale_bilinear_forms(&z, &u, || poll(control))?;
                let lifting = match case.prescribed {
                    Some(g) => op.prescribed_displacement_scale_work(g, &z, || poll(control))?,
                    None => vec![0.0; self.cells()],
                };
                for ((s, k), b) in local.iter_mut().zip(stiffness).zip(lifting) { *s = finite(*s+b-k)?; }
                external_work.push(finite(dot(&external, &u))?);
                augmented_rhs_work.push(finite(dot(&rhs, &u))?);
                displacements.push(u); adjoints.push(z); responses.push(observed);
            }
            let p = self.params;
            for ((value, &r), &slope) in local.iter_mut().zip(&design.projected).zip(&design.slope) {
                let power = if p.penal == 1.0 { 1.0 } else if r == 0.0 { 0.0 } else { fs_math::det::pow(r,p.penal-1.0) };
                *value = finite(*value*(1.0-p.e_min)*p.penal*power*slope)?;
            }
            let mut gradient = self.pullback(&local, control)?;
            let dv: Vec<_> = self.mass.iter().zip(&design.slope).map(|(m,s)| m*s).collect();
            let volume_gradient = self.pullback(&dv, control)?;
            objective = finite(objective+options.volume_weight*design.volume)?;
            for (g,v) in gradient.iter_mut().zip(&volume_gradient) { *g = finite(*g+options.volume_weight*v)?; }
            control.checkpoint("sdf3-response-publish")?;
            Ok(ResponseEvaluation3 { rho: rho.to_vec(), projected_rho: design.projected, scales: design.scales,
                objective, gradient, volume_fraction: design.volume, volume_gradient, responses,
                displacements, adjoints, external_work, augmented_rhs_work, work: control.work() })
        })();
        self.operator.set_scales(&previous).expect("previous admitted scales remain valid");
        result
    }
}
