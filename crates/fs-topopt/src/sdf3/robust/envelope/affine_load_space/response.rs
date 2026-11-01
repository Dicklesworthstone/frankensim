//! Exact support-function response limits over the same continuous load ball.
//!
//! For y(a)=y0+v^T a and ||a||<=1, each signed inequality has worst residual
//! [sign*(y0-bound)+||v||]/scale. Each response has its OWN adverse direction;
//! the worst-compliance direction is not a surrogate for serviceability.
//! Native embedded reaction state AND direct material terms are differentiated.
use super::grid::{finite, force_work, normalize, relative_difference};
use super::{AffineLoadSpace3, AffineLoadSpaceEvaluation3, AffineLoadSpaceGrid3, RobustOptions3};
use crate::sdf3::response::ResponseError3;
use crate::sdf3::response::limits::{
    BoundedResponse3, ResponseLimit3, ResponseLimitSense3, aggregate,
};
use crate::sdf3::{
    CutDensityStudy3, Design3, Sdf3Elasticity, checked_solve_preconditioned, failure,
};
use crate::{EvaluationStop, SolveControl, SolveWork};
use fs_cutfem::elastic3::ElasticityError3;
use fs_sparse::precond::Precond;
use std::borrow::Cow;
use std::ops::ControlFlow;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AffineResponseOptions3 {
    /// Total original upper/lower inequalities, including both sides of windows.
    pub max_limits: usize,
    /// Positive inverse width of the UNNORMALIZED smooth upper maximum.
    pub sharpness: f64,
}
impl Default for AffineResponseOptions3 {
    fn default() -> Self {
        Self {
            max_limits: 256,
            sharpness: 40.0,
        }
    }
}

/// The same admitted observation types, units and signed-bound conventions as
/// finite-case response limits. Every inequality applies to the ENTIRE ball.
/// Observation vectors and reaction modes are pure, density-independent
/// reference laws in the current native displacement/support coordinates.
#[derive(Clone, Copy)]
pub struct AffineResponseLimits3<'a> {
    pub limits: &'a [ResponseLimit3<'a>],
    pub options: AffineResponseOptions3,
}

#[derive(Debug, Clone)]
pub struct AffineResponseWitness3 {
    pub baseline: f64,
    pub basis_responses: Vec<f64>,
    /// Euclidean norm of the basis responses. Zero has no smooth norm gradient.
    pub radius: f64,
    pub coefficients: Vec<f64>,
    /// Closed-form adverse physical extremum: baseline + sign*radius.
    pub worst_value: f64,
    /// Same observation re-evaluated on the actual combined native field.
    pub witness_value: f64,
    pub bound: f64,
    pub sense: ResponseLimitSense3,
    pub scale: f64,
    pub residual: f64,
    pub witness_residual: f64,
    /// Numerical support/witness identity, scaled by baseline/radius/extrema.
    /// This is neither an interval enclosure nor near-zero relative accuracy.
    pub relative_defect: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AffineResponseAdjointLayout3 {
    /// Two-grid assessment requested values only; no adjoint was attempted.
    NotRequested,
    /// Some zero-radius response prevents a smooth aggregate gradient.
    Unavailable,
    /// One weighted signed adjoint in original observation order.
    Observations,
    /// Aggregate adjoints in [offset, basis0, ...] order. Only offset has B^T z.
    AffineFields,
}

#[derive(Debug, Clone)]
pub struct AffineResponseEvaluation3 {
    pub options: AffineResponseOptions3,
    pub constraint: f64,
    pub gradient: Option<Vec<f64>>,
    /// Worst original support-function OR actual witness residual.
    pub worst_residual: f64,
    pub witnesses: Vec<AffineResponseWitness3>,
    pub adjoint_layout: AffineResponseAdjointLayout3,
    /// At most min(limit_count, basis_count+1) adjoints, under the stated layout.
    pub adjoints: Vec<Vec<f64>>,
    pub work: SolveWork,
}
impl AffineResponseEvaluation3 {
    /// Require the conservative aggregate AND every original physical bound.
    /// No optimizer tolerance excuses a physical-cap excess. This is discrete
    /// linear-model evidence, not a continuum or floating-point enclosure.
    #[must_use]
    pub fn physically_feasible(&self) -> bool {
        !self.witnesses.is_empty()
            && self.constraint.is_finite()
            && self.constraint <= 0.0
            && self.worst_residual.is_finite()
            && self.worst_residual <= 0.0
            && self.witnesses.iter().all(|w| {
                [
                    w.worst_value,
                    w.witness_value,
                    w.bound,
                    w.residual,
                    w.witness_residual,
                ]
                .iter()
                .all(|v| v.is_finite())
                    && w.residual <= 0.0
                    && w.witness_residual <= 0.0
                    && match w.sense {
                        ResponseLimitSense3::AtMost => {
                            w.worst_value <= w.bound && w.witness_value <= w.bound
                        }
                        ResponseLimitSense3::AtLeast => {
                            w.worst_value >= w.bound && w.witness_value >= w.bound
                        }
                    }
            })
    }
}

pub(super) fn admit<O: Sdf3Elasticity>(
    study: &CutDensityStudy3<O>,
    family: AffineResponseLimits3<'_>,
) -> Result<(), EvaluationStop> {
    admit_operator(&study.operator, family)
}

pub(super) fn admit_policy(
    options: AffineResponseOptions3,
    count: usize,
) -> Result<(), EvaluationStop> {
    if count == 0
        || count > options.max_limits
        || !options.sharpness.is_finite()
        || options.sharpness <= 0.0
        || !(1.0 / options.sharpness).is_finite()
    {
        return Err(failure("sdf3-affine-response-policy"));
    }
    Ok(())
}

pub(super) fn admit_bound(bound: f64, scale: f64) -> Result<(), EvaluationStop> {
    if !bound.is_finite() || !scale.is_finite() || scale <= 0.0 || !(1.0 / scale).is_finite() {
        return Err(failure("sdf3-affine-response-bound"));
    }
    Ok(())
}

fn admit_operator(
    operator: &impl Sdf3Elasticity,
    family: AffineResponseLimits3<'_>,
) -> Result<(), EvaluationStop> {
    admit_policy(family.options, family.limits.len())?;
    for limit in family.limits {
        admit_bound(limit.bound, limit.scale)?;
        match limit.response {
            BoundedResponse3::Displacement(q) => {
                if q.len() != operator.n() || q.iter().any(|v| !v.is_finite()) {
                    return Err(failure("sdf3-affine-response-observation"));
                }
            }
            BoundedResponse3::Reaction(_) if operator.embedded_dirichlet_penalty().is_none() => {
                return Err(failure("sdf3-affine-response-embedded-support"));
            }
            BoundedResponse3::Reaction(_) => {}
        }
    }
    Ok(())
}

fn controlled<T>(
    control: &mut SolveControl<'_>,
    action: impl FnOnce(&mut dyn FnMut() -> ControlFlow<()>) -> Result<T, ElasticityError3>,
) -> Result<T, EvaluationStop> {
    control.checkpoint("sdf3-affine-response-physics")?;
    let mut stopped = None;
    let result = action(
        &mut || match control.checkpoint("sdf3-affine-response-physics") {
            Ok(()) => ControlFlow::Continue(()),
            Err(error) => {
                stopped = Some(error);
                ControlFlow::Break(())
            }
        },
    );
    if let Some(error) = stopped {
        return Err(error);
    }
    result.map_err(|_| failure("sdf3-affine-response-physics"))
}

struct Trace<'a> {
    value: f64,
    gradient: Cow<'a, [f64]>,
    direct: Option<Vec<f64>>,
}

fn trace<'a>(
    operator: &impl Sdf3Elasticity,
    response: BoundedResponse3<'a>,
    u: &[f64],
    control: &mut SolveControl<'_>,
) -> Result<Trace<'a>, EvaluationStop> {
    match response {
        BoundedResponse3::Displacement(q) => Ok(Trace {
            value: force_work(q, u, operator.fixed(), control)?,
            gradient: Cow::Borrowed(q),
            direct: None,
        }),
        BoundedResponse3::Reaction(mode) => {
            let native = controlled(control, |poll| {
                operator.embedded_reaction(u, None, mode, poll)
            })?;
            Ok(Trace {
                value: native.value,
                gradient: Cow::Owned(native.displacement_gradient),
                direct: Some(native.scale_gradient),
            })
        }
    }
}

fn norm(values: &[f64]) -> Result<f64, EvaluationStop> {
    let scale = values.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return Ok(0.0);
    }
    let sum: f64 = values.iter().map(|v| (v / scale) * (v / scale)).sum();
    finite(
        scale * fs_math::det::sqrt(sum),
        "sdf3-affine-response-radius",
    )
}

fn displacement(
    grid: &AffineLoadSpaceGrid3,
    coefficients: &[f64],
    control: &mut SolveControl<'_>,
) -> Result<Vec<f64>, EvaluationStop> {
    let mut u = grid.offset_displacement.clone();
    for (i, value) in u.iter_mut().enumerate() {
        if i % 256 == 0 {
            control.checkpoint("sdf3-affine-response-witness")?;
        }
        for (a, column) in coefficients.iter().zip(&grid.basis_displacements) {
            *value += a * column[i];
        }
        finite(*value, "sdf3-affine-response-witness")?;
    }
    Ok(u)
}

fn product(a: f64, b: f64) -> Result<f64, EvaluationStop> {
    let value = a * b;
    if !value.is_finite() || (a != 0.0 && b != 0.0 && value == 0.0) {
        Err(failure("sdf3-affine-response-gradient-range"))
    } else {
        Ok(value)
    }
}

fn measure(
    operator: &impl Sdf3Elasticity,
    grid: &AffineLoadSpaceGrid3,
    family: AffineResponseLimits3<'_>,
    control: &mut SolveControl<'_>,
) -> Result<Vec<AffineResponseWitness3>, EvaluationStop> {
    let mut witnesses = Vec::with_capacity(family.limits.len());
    for limit in family.limits {
        control.checkpoint("sdf3-affine-response-measure")?;
        let baseline = trace(operator, limit.response, &grid.offset_displacement, control)?;
        let mut basis_responses = Vec::with_capacity(grid.basis_displacements.len());
        // Homogeneous support makes the native reaction exactly linear in u;
        // reuse its state row instead of re-integrating every basis trace.
        for u in &grid.basis_displacements {
            basis_responses.push(force_work(
                &baseline.gradient,
                u,
                operator.fixed(),
                control,
            )?);
        }
        let radius = norm(&basis_responses)?;
        let coefficients = if radius == 0.0 {
            vec![0.0; basis_responses.len()]
        } else {
            basis_responses
                .iter()
                .map(|v| {
                    normalize(
                        limit.sense.sign() * v,
                        radius,
                        "sdf3-affine-response-direction",
                    )
                })
                .collect::<Result<Vec<_>, _>>()?
        };
        let worst_value = finite(
            baseline.value + limit.sense.sign() * radius,
            "sdf3-affine-response-extremum",
        )?;
        let u = displacement(grid, &coefficients, control)?;
        let witness_value = trace(operator, limit.response, &u, control)?.value;
        let scale = baseline
            .value
            .abs()
            .max(radius)
            .max(worst_value.abs())
            .max(witness_value.abs());
        let relative_defect = relative_difference(worst_value, witness_value, scale);
        if !relative_defect.is_finite() || relative_defect > 1e-9 {
            return Err(failure("sdf3-affine-response-identity"));
        }
        // Preserve a small positive radius when baseline==bound even if adding
        // it to a large baseline rounds the displayed physical extremum away.
        let residual = normalize(
            finite(
                limit.sense.sign()
                    * finite(
                        baseline.value - limit.bound,
                        "sdf3-affine-response-residual",
                    )?
                    + radius,
                "sdf3-affine-response-residual",
            )?,
            limit.scale,
            "sdf3-affine-response-residual",
        )?;
        let witness_residual = normalize(
            limit.sense.sign()
                * finite(witness_value - limit.bound, "sdf3-affine-response-residual")?,
            limit.scale,
            "sdf3-affine-response-residual",
        )?;
        witnesses.push(AffineResponseWitness3 {
            baseline: baseline.value,
            basis_responses,
            radius,
            coefficients,
            worst_value,
            witness_value,
            bound: limit.bound,
            sense: limit.sense,
            scale: limit.scale,
            residual,
            witness_residual,
            relative_defect,
        });
    }
    Ok(witnesses)
}

fn add_direct(
    local: &mut [f64],
    direct: Option<Vec<f64>>,
    alpha: f64,
    control: &mut SolveControl<'_>,
) -> Result<(), EvaluationStop> {
    if let Some(direct) = direct {
        for (i, (value, d)) in local.iter_mut().zip(direct).enumerate() {
            if i % 256 == 0 {
                control.checkpoint("sdf3-affine-response-direct")?;
            }
            *value = finite(*value + product(alpha, d)?, "sdf3-affine-response-direct")?;
        }
    }
    Ok(())
}

fn subtract_stiffness<O: Sdf3Elasticity>(
    study: &CutDensityStudy3<O>,
    z: &[f64],
    u: &[f64],
    local: &mut [f64],
    control: &mut SolveControl<'_>,
) -> Result<(), EvaluationStop> {
    let stiffness = controlled(control, |poll| {
        study.operator.scale_bilinear_forms(z, u, poll)
    })?;
    for (i, (value, k)) in local.iter_mut().zip(stiffness).enumerate() {
        if i % 256 == 0 {
            control.checkpoint("sdf3-affine-response-stiffness")?;
        }
        *value = finite(*value - k, "sdf3-affine-response-stiffness")?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn gradient<O: Sdf3Elasticity>(
    study: &CutDensityStudy3<O>,
    design: &Design3,
    prepared: &impl Precond,
    grid: &AffineLoadSpaceGrid3,
    family: AffineResponseLimits3<'_>,
    witnesses: &[AffineResponseWitness3],
    weights: &[f64],
    control: &mut SolveControl<'_>,
) -> Result<(Vec<f64>, Vec<Vec<f64>>, AffineResponseAdjointLayout3), EvaluationStop> {
    let n = study.operator.n();
    let fields = grid.basis_displacements.len() + 1;
    let mut local = vec![0.0; study.cells()];
    let mut occupancy = vec![0.0; study.cells()];
    let mut adjoints = Vec::with_capacity(family.limits.len().min(fields));
    let layout = if family.limits.len() <= fields {
        // Few observations: each has a distinct frozen worst field and one
        // adjoint. Its force derivative contains the SAME single body map B.
        for ((limit, witness), &weight) in family.limits.iter().zip(witnesses).zip(weights) {
            control.checkpoint("sdf3-affine-response-pullback")?;
            let alpha = normalize(
                weight * limit.sense.sign(),
                limit.scale,
                "sdf3-affine-response-weight",
            )?;
            let u = displacement(grid, &witness.coefficients, control)?;
            let observation = trace(&study.operator, limit.response, &u, control)?;
            add_direct(&mut local, observation.direct, alpha, control)?;
            let mut rhs = vec![0.0; n];
            for (i, (&q, value)) in observation.gradient.iter().zip(&mut rhs).enumerate() {
                if i % 256 == 0 {
                    control.checkpoint("sdf3-affine-response-adjoint-load")?;
                }
                if !study.operator.fixed()[i / 3] {
                    *value = product(alpha, q)?;
                }
            }
            let z = checked_solve_preconditioned(
                &study.operator,
                prepared,
                &rhs,
                1e-11,
                "sdf3-affine-response-adjoint",
                control,
            )?;
            subtract_stiffness(study, &z, &u, &mut local, control)?;
            study.add_material_body_gradient(&z, 1.0, &mut occupancy, control)?;
            adjoints.push(z);
        }
        AffineResponseAdjointLayout3::Observations
    } else {
        // Many observations: factor the aggregate by the existing affine
        // primal family. This uses m+1 adjoints, not one solve per observation.
        let mut rhs = vec![vec![0.0; n]; fields];
        for ((limit, witness), &weight) in family.limits.iter().zip(witnesses).zip(weights) {
            control.checkpoint("sdf3-affine-response-pullback")?;
            let alpha = normalize(
                weight * limit.sense.sign(),
                limit.scale,
                "sdf3-affine-response-weight",
            )?;
            let u = displacement(grid, &witness.coefficients, control)?;
            let observation = trace(&study.operator, limit.response, &u, control)?;
            add_direct(&mut local, observation.direct, alpha, control)?;
            for (i, &q) in observation.gradient.iter().enumerate() {
                if i % 256 == 0 {
                    control.checkpoint("sdf3-affine-response-adjoint-load")?;
                }
                if study.operator.fixed()[i / 3] {
                    continue;
                }
                let base = product(alpha, q)?;
                rhs[0][i] = finite(rhs[0][i] + base, "sdf3-affine-response-adjoint-load")?;
                for (j, &a) in witness.coefficients.iter().enumerate() {
                    rhs[j + 1][i] = finite(
                        rhs[j + 1][i] + product(base, a)?,
                        "sdf3-affine-response-adjoint-load",
                    )?;
                }
            }
        }
        for (j, rhs) in rhs.into_iter().enumerate() {
            let z = checked_solve_preconditioned(
                &study.operator,
                prepared,
                &rhs,
                1e-11,
                "sdf3-affine-response-adjoint",
                control,
            )?;
            let u = if j == 0 {
                &grid.offset_displacement
            } else {
                &grid.basis_displacements[j - 1]
            };
            subtract_stiffness(study, &z, u, &mut local, control)?;
            if j == 0 {
                study.add_material_body_gradient(&z, 1.0, &mut occupancy, control)?;
            }
            adjoints.push(z);
        }
        AffineResponseAdjointLayout3::AffineFields
    };
    let gradient = study.physical_gradient(design, &local, &occupancy, control)?;
    Ok((gradient, adjoints, layout))
}

fn measure_evaluation(
    operator: &impl Sdf3Elasticity,
    grid: &AffineLoadSpaceGrid3,
    family: AffineResponseLimits3<'_>,
    control: &mut SolveControl<'_>,
) -> Result<(AffineResponseEvaluation3, Vec<f64>), EvaluationStop> {
    control.checkpoint("sdf3-affine-response-assess")?;
    admit_operator(operator, family)?;
    if grid.basis_displacements.is_empty()
        || grid.basis_displacements.len() > fs_ascent::quadratic_ball::MAX_QUADRATIC_BALL_DIMENSION
    {
        return Err(failure("sdf3-affine-response-field-shape"));
    }
    for field in std::iter::once(&grid.offset_displacement).chain(&grid.basis_displacements) {
        if field.len() != operator.n() {
            return Err(failure("sdf3-affine-response-field-shape"));
        }
        for (i, value) in field.iter().enumerate() {
            if i % 256 == 0 {
                control.checkpoint("sdf3-affine-response-field-shape")?;
            }
            if !value.is_finite() {
                return Err(failure("sdf3-affine-response-field-shape"));
            }
        }
    }
    let witnesses = measure(operator, grid, family, control)?;
    let residuals: Vec<_> = witnesses.iter().map(|w| w.residual).collect();
    let (constraint, mut worst_residual, weights) = aggregate(&residuals, family.options.sharpness)
        .map_err(|error| match error {
            ResponseError3::Evaluation(stop) => stop,
            _ => failure("sdf3-affine-response-aggregate"),
        })?;
    for witness in &witnesses {
        worst_residual = worst_residual.max(witness.witness_residual);
    }
    Ok((
        AffineResponseEvaluation3 {
            options: family.options,
            constraint,
            gradient: None,
            worst_residual,
            witnesses,
            adjoint_layout: AffineResponseAdjointLayout3::NotRequested,
            adjoints: Vec::new(),
            work: control.work(),
        },
        weights,
    ))
}

/// Same-transaction native grid assessment without filtering or solving any
/// additional primal/adjoint. The caller owns the actual material installation
/// and fresh affine fields; this private seam never admits public field caches.
/// Zero-radius witnesses remain valid because no smooth gradient is requested.
pub(super) fn measure_affine_responses(
    operator: &impl Sdf3Elasticity,
    grid: &AffineLoadSpaceGrid3,
    family: AffineResponseLimits3<'_>,
    control: &mut SolveControl<'_>,
) -> Result<AffineResponseEvaluation3, EvaluationStop> {
    measure_evaluation(operator, grid, family, control).map(|(report, _)| report)
}

fn evaluate<O: Sdf3Elasticity>(
    study: &CutDensityStudy3<O>,
    design: &Design3,
    prepared: &impl Precond,
    grid: &AffineLoadSpaceGrid3,
    family: AffineResponseLimits3<'_>,
    control: &mut SolveControl<'_>,
) -> Result<AffineResponseEvaluation3, EvaluationStop> {
    let (mut report, weights) = measure_evaluation(&study.operator, grid, family, control)?;
    if report.witnesses.iter().all(|w| w.radius > 0.0) {
        let (gradient, adjoints, layout) = gradient(
            study,
            design,
            prepared,
            grid,
            family,
            &report.witnesses,
            &weights,
            control,
        )?;
        report.gradient = Some(gradient);
        report.adjoints = adjoints;
        report.adjoint_layout = layout;
    } else {
        report.adjoint_layout = AffineResponseAdjointLayout3::Unavailable;
    }
    report.work = control.work();
    Ok(report)
}

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Joint compliance and response assessment from ONE native affine primal
    /// family and ONE eroded preparation. All original response bounds cover
    /// the entire load ball. The exact discrete material derivative uses native
    /// stiffness contractions, reaction direct terms and B^T z, followed by the
    /// original physical-region/projection/filter chain. No sampled directions,
    /// finite-difference optimization or gradients through Krylov are used.
    /// A zero response radius retains valid values but no smooth gradient.
    /// All returned outcomes preserve incoming scales and cumulative work.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_robust_affine_response_limits(
        &mut self,
        rho: &[f64],
        dead_offset: &[f64],
        basis: &[&[f64]],
        robustness: RobustOptions3,
        policy: &AffineLoadSpace3,
        family: AffineResponseLimits3<'_>,
        control: &mut SolveControl<'_>,
    ) -> Result<AffineLoadSpaceEvaluation3, EvaluationStop> {
        self.affine_load_space_response_with_mass(
            rho,
            dead_offset,
            basis,
            robustness,
            policy,
            family,
            false,
            control,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn affine_load_space_response_with_mass(
        &mut self,
        rho: &[f64],
        dead_offset: &[f64],
        basis: &[&[f64]],
        robustness: RobustOptions3,
        policy: &AffineLoadSpace3,
        family: AffineResponseLimits3<'_>,
        include_mass: bool,
        control: &mut SolveControl<'_>,
    ) -> Result<AffineLoadSpaceEvaluation3, EvaluationStop> {
        control.checkpoint("sdf3-affine-response-start")?;
        admit(self, family)?;
        let (mut result, mut responses) = self.affine_load_space_extended(
            rho,
            dead_offset,
            basis,
            robustness,
            policy,
            include_mass,
            control,
            |study, design, prepared, grid, control| {
                evaluate(study, design, prepared, grid, family, control)
            },
        )?;
        control.checkpoint("sdf3-affine-response-publish")?;
        responses.work = control.work();
        result.work = control.work();
        result.responses = Some(responses);
        Ok(result)
    }
}
