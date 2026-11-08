//! Constrained allocation of element Young's-modulus multipliers.
//!
//! Minimize the caller's positive weighted modulus cost under a lower bound
//! on the first positive proportional buckling factor. Every trial re-solves
//! the reference equilibrium and uses its TOTAL buckling derivative, including
//! prestress redistribution. The optimizer is the existing projected AL owner.
//!
//! This is fixed-geometry elastic-stiffness allocation. A modulus multiplier
//! is not a physical mass, density, thickness or section-area law. P1/Q1,
//! plane stress/strain and Standard/B-bar follow the native solid problem.
//! The original first-mode entry requires a separated mode. The explicit
//! aggregate entry constrains a conservative smooth KS lower envelope of a
//! complete mode packet, including repeated modes inside that packet. Both
//! consume the same native optimizer state and total equilibrium derivatives.
use fs_ascent::projected_al::{
    ProjectedAlError, ProjectedAlOptions, ProjectedAlReport, ProjectedAlSample, ProjectedAlState,
    ProjectedAlStop, ProjectedAlWork,
};
use fs_solid::SolidError;
use fs_solid::linear::LinearProblem;
use fs_solid::stability::sensitivity::aggregate::{
    ElementScaleAggregateSensitivity, element_scale_aggregate_sensitivities_cancellable,
};
use fs_solid::stability::sensitivity::{
    ElementScaleSensitivity, element_scale_sensitivities_cancellable,
};
use std::ops::ControlFlow;

/// Physical buckling requirement and existing numerical optimizer allowances.
#[derive(Debug, Clone, Copy)]
pub struct BucklingDesignOptions {
    /// Positive dimensionless multiplier of the problem's reference loading.
    pub minimum_load_factor: f64,
    /// Native eigen-iteration allowance per complete forward/adjoint evaluation.
    /// Native static-solve, residual, separation and factor-size limits remain.
    pub spectral_steps: usize,
    pub optimizer: ProjectedAlOptions,
}
impl Default for BucklingDesignOptions {
    fn default() -> Self {
        Self {
            minimum_load_factor: 1.0,
            spectral_steps: 600,
            optimizer: ProjectedAlOptions::default(),
        }
    }
}

/// One fully evaluated design. Its cost, derivative and native physics belong
/// to the same immutable element scales and original loading problem.
#[derive(Debug, Clone)]
pub struct BucklingDesignPoint<Physics = ElementScaleSensitivity> {
    physics: Physics,
    sample: ProjectedAlSample,
    meets_floor: bool,
}
impl<Physics> BucklingDesignPoint<Physics> {
    /// Complete equilibrium-dependent mode, total derivatives and residuals.
    #[must_use]
    pub fn physics(&self) -> &Physics {
        &self.physics
    }
    /// Cost and normalized constraint with their exact discrete derivatives.
    /// The inequality is 1 - selected_factor/minimum_load_factor <= 0.
    /// The aggregate entry selects its conservative KS lower envelope.
    #[must_use]
    pub fn sample(&self) -> &ProjectedAlSample {
        &self.sample
    }
    /// True only when the selected critical factor or KS envelope meets the floor;
    /// optimizer feasibility tolerance cannot excuse a physical shortfall.
    #[must_use]
    pub fn meets_load_floor(&self) -> bool {
        self.meets_floor
    }
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}

fn evaluate(
    problem: &LinearProblem<'_>,
    weights: &[f64],
    scales: &[f64],
    options: BucklingDesignOptions,
    cancelled: &impl Fn() -> bool,
) -> Result<BucklingDesignPoint, SolidError> {
    let physics = element_scale_sensitivities_cancellable(
        problem,
        scales,
        1,
        options.spectral_steps,
        cancelled,
    )?;
    let factor = *physics
        .buckling()
        .loads
        .first()
        .ok_or_else(|| invalid("buckling design requires a positive separated mode"))?;
    let sample = design_sample(weights, scales, factor, &physics.gradients()[0], options)?;
    Ok(BucklingDesignPoint {
        meets_floor: factor >= options.minimum_load_factor,
        physics,
        sample,
    })
}

// A KS envelope can be negative even when all individual critical factors are
// positive. Such a finite, infeasible initial point must remain optimizable.
fn design_sample(
    weights: &[f64],
    scales: &[f64],
    factor: f64,
    derivative: &[f64],
    options: BucklingDesignOptions,
) -> Result<ProjectedAlSample, SolidError> {
    let normalized = factor / options.minimum_load_factor;
    if !normalized.is_finite() || (factor != 0.0 && normalized == 0.0) {
        return Err(invalid(
            "buckling design load normalization is not representable",
        ));
    }
    let mut constraint_gradient = Vec::with_capacity(weights.len());
    for &value in derivative {
        let normalized = -value / options.minimum_load_factor;
        if !normalized.is_finite() || (value != 0.0 && normalized == 0.0) {
            return Err(invalid(
                "buckling design derivative normalization is not representable",
            ));
        }
        constraint_gradient.push(normalized);
    }
    let cost: f64 = weights.iter().zip(scales).map(|(w, s)| w * s).sum();
    if !cost.is_finite() || cost <= 0.0 {
        return Err(invalid("weighted modulus cost is not representable"));
    }
    Ok(ProjectedAlSample {
        objective: cost,
        gradient: weights.to_vec(),
        constraint: 1.0 - normalized,
        constraint_gradient,
    })
}

/// Resumable first-mode buckling-constrained modulus allocation.
///
/// Accepted steps need not be physically feasible or lower the original cost.
/// best_feasible() separately retains the least-cost accepted point meeting
/// the actual load floor. Its native physics is retained, not reconstructed
/// from an optimizer scalar or a frozen initial stress field.
///
/// Geometry, Poisson ratio and the reference force/prescribed-displacement
/// functions must remain unchanged for this session. No mesh, load, material
/// or support is mutated. An unresolved spectral crossing, native solve error
/// or cancellation cannot install a partial trial. Already completed steps
/// remain available if a later trial stops. Numerical KKT diagnostics concern
/// this discrete smooth branch; they do not certify continuum stability.
pub struct BucklingDesignStudy<'problem, 'mesh> {
    problem: &'problem LinearProblem<'mesh>,
    weights: Vec<f64>,
    options: BucklingDesignOptions,
    engine: DesignState<ElementScaleSensitivity>,
}

impl<'problem, 'mesh> BucklingDesignStudy<'problem, 'mesh> {
    /// Bind a positive scalar weight and finite positive box to each original
    /// element. Weights declare both the resource cost and its numerical scale;
    /// no physical mass or volume interpretation is supplied.
    ///
    /// # Errors
    /// Invalid inputs, absent/clustered positive modes and native failures refuse.
    pub fn new(
        problem: &'problem LinearProblem<'mesh>,
        initial: &[f64],
        lower: &[f64],
        upper: &[f64],
        weights: &[f64],
        options: BucklingDesignOptions,
    ) -> Result<Self, ProjectedAlError<SolidError>> {
        Self::new_cancellable(problem, initial, lower, upper, weights, options, &|| false)
    }

    /// Cancellable construction using the original native forward/adjoint owner.
    ///
    /// # Errors
    /// As new(), with cancellation before or during any native numerical phase.
    #[allow(clippy::too_many_arguments)]
    pub fn new_cancellable(
        problem: &'problem LinearProblem<'mesh>,
        initial: &[f64],
        lower: &[f64],
        upper: &[f64],
        weights: &[f64],
        options: BucklingDesignOptions,
        cancelled: &impl Fn() -> bool,
    ) -> Result<Self, ProjectedAlError<SolidError>> {
        validate_design(problem, initial, lower, upper, weights, options)?;
        let engine = initialize(
            initial,
            lower,
            upper,
            options.optimizer,
            &mut |scales| evaluate(problem, weights, scales, options, cancelled),
            cancelled,
        )?;
        Ok(Self {
            problem,
            weights: weights.to_vec(),
            options,
            engine,
        })
    }

    /// The last accepted original element multipliers.
    #[must_use]
    pub fn element_scales(&self) -> &[f64] {
        self.engine.state.point()
    }
    /// Physics and full derivative at exactly the current optimizer point.
    #[must_use]
    pub fn accepted(&self) -> &BucklingDesignPoint {
        &self.engine.accepted
    }
    /// Lowest-cost accepted design satisfying the ORIGINAL positive load floor.
    #[must_use]
    pub fn best_feasible(&self) -> Option<&BucklingDesignPoint> {
        self.engine.best_feasible.as_ref()
    }
    /// Cumulative work includes rejected and interrupted physics attempts.
    #[must_use]
    pub fn work(&self) -> ProjectedAlWork {
        self.engine.state.work()
    }

    /// Advance at most this many accepted steps using the same optimizer state.
    ///
    /// # Errors
    /// Native failures and unresolved crossings stop without installing that trial.
    pub fn run(
        &mut self,
        additional_steps: usize,
    ) -> Result<ProjectedAlReport, ProjectedAlError<SolidError>> {
        self.run_cancellable(additional_steps, &|| false)
    }

    /// Advance with cancellation through native equilibrium, buckling and adjoints.
    /// A zero-step call reports cached state without additional physics.
    ///
    /// # Errors
    /// As run(); cancellation and native failures remain explicit stops.
    pub fn run_cancellable(
        &mut self,
        additional_steps: usize,
        cancelled: &impl Fn() -> bool,
    ) -> Result<ProjectedAlReport, ProjectedAlError<SolidError>> {
        let (problem, weights, options) = (self.problem, &self.weights, self.options);
        self.engine.advance(
            additional_steps,
            cancelled,
            ElementScaleSensitivity::element_scales,
            &mut |scales| evaluate(problem, weights, scales, options, cancelled),
        )
    }
}

fn validate_design(
    problem: &LinearProblem<'_>,
    initial: &[f64],
    lower: &[f64],
    upper: &[f64],
    weights: &[f64],
    options: BucklingDesignOptions,
) -> Result<(), ProjectedAlError<SolidError>> {
    options.optimizer.validate(initial.len())?;
    let n = problem.mesh.elems.len();
    if n == 0
        || initial.len() != n
        || lower.len() != n
        || upper.len() != n
        || weights.len() != n
        || options.spectral_steps == 0
        || !options.minimum_load_factor.is_finite()
        || options.minimum_load_factor <= 0.0
    {
        return Err(ProjectedAlError::Invalid(
            "invalid buckling design shape or load policy",
        ));
    }
    let mut maximum_cost = 0.0;
    for i in 0..n {
        let (lo, hi, s, w) = (lower[i], upper[i], initial[i], weights[i]);
        if !lo.is_finite()
            || !hi.is_finite()
            || !s.is_finite()
            || lo <= 0.0
            || lo > hi
            || s < lo
            || s > hi
            || !w.is_finite()
            || w <= 0.0
            || w * lo == 0.0
            || !(w * hi).is_finite()
        {
            return Err(ProjectedAlError::Invalid(
                "positive finite modulus box and cost weights required",
            ));
        }
        maximum_cost += w * hi;
    }
    if !maximum_cost.is_finite() {
        return Err(ProjectedAlError::Invalid(
            "maximum weighted modulus cost overflows",
        ));
    }
    Ok(())
}

struct DesignState<Physics> {
    state: ProjectedAlState,
    accepted: BucklingDesignPoint<Physics>,
    best_feasible: Option<BucklingDesignPoint<Physics>>,
}

fn initialize<Physics: Clone>(
    initial: &[f64],
    lower: &[f64],
    upper: &[f64],
    optimizer: ProjectedAlOptions,
    evaluate: &mut impl FnMut(&[f64]) -> Result<BucklingDesignPoint<Physics>, SolidError>,
    cancelled: &impl Fn() -> bool,
) -> Result<DesignState<Physics>, ProjectedAlError<SolidError>> {
    let mut initial_physics = None;
    let state = ProjectedAlState::try_new(
        initial,
        lower,
        upper,
        optimizer,
        &mut |scales| {
            let point = evaluate(scales)?;
            let sample = point.sample.clone();
            initial_physics = Some(point);
            Ok(Some(sample))
        },
        |_| {
            if cancelled() {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        },
    )?;
    let accepted = initial_physics.expect("accepted initial sample has matching native physics");
    let best_feasible = accepted.meets_floor.then(|| accepted.clone());
    Ok(DesignState {
        state,
        accepted,
        best_feasible,
    })
}

impl<Physics: Clone> DesignState<Physics> {
    fn advance(
        &mut self,
        additional_steps: usize,
        cancelled: &impl Fn() -> bool,
        native_scales: fn(&Physics) -> &[f64],
        evaluate: &mut impl FnMut(&[f64]) -> Result<BucklingDesignPoint<Physics>, SolidError>,
    ) -> Result<ProjectedAlReport, ProjectedAlError<SolidError>> {
        let mut remaining = additional_steps;
        loop {
            let before = self.state.work().iterations;
            let mut last = None;
            let result = self.state.try_run(
                remaining.min(1),
                &mut |scales| {
                    // No error becomes an unavailable or frozen-gradient sample.
                    let point = evaluate(scales)?;
                    let sample = point.sample.clone();
                    last = Some(point);
                    Ok(Some(sample))
                },
                |_| {
                    if cancelled() {
                        ControlFlow::Break(())
                    } else {
                        ControlFlow::Continue(())
                    }
                },
            );
            let completed = self.state.work().iterations - before;
            if completed > 0 {
                let accepted = last.expect("accepted step has matching native buckling physics");
                assert_eq!(native_scales(&accepted.physics), self.state.point());
                if accepted.meets_floor
                    && self
                        .best_feasible
                        .as_ref()
                        .is_none_or(|best| accepted.sample.objective < best.sample.objective)
                {
                    self.best_feasible = Some(accepted.clone());
                }
                self.accepted = accepted;
                remaining -= completed;
            }
            let report = result?;
            if remaining == 0 || completed == 0 || report.stop != ProjectedAlStop::IterationLimit {
                return Ok(report);
            }
        }
    }
}

/// Fixed complete-mode policy for a smooth conservative buckling constraint.
///
/// For m retained loads, min(loads) - ln(m)/sharpness <= KS <= min(loads).
/// Increasing sharpness reduces this conservative bias and sharpens the
/// optimization landscape. Keep the policy fixed throughout a study.
#[derive(Debug, Clone, Copy)]
pub struct AggregateBucklingDesignOptions {
    /// Original physical load floor, spectral allowance and native AL policy.
    pub design: BucklingDesignOptions,
    /// Number of lowest positive modes. The packet must contain entire clusters.
    pub modes: usize,
    /// Finite positive inverse load-factor scale in the KS exponential.
    pub sharpness: f64,
}
impl Default for AggregateBucklingDesignOptions {
    fn default() -> Self {
        Self {
            design: BucklingDesignOptions::default(),
            modes: 2,
            sharpness: 32.0,
        }
    }
}

/// A conservative complete-packet design and its single-adjoint total gradient.
pub type AggregateBucklingDesignPoint = BucklingDesignPoint<ElementScaleAggregateSensitivity>;

fn evaluate_aggregate(
    problem: &LinearProblem<'_>,
    weights: &[f64],
    scales: &[f64],
    options: AggregateBucklingDesignOptions,
    cancelled: &impl Fn() -> bool,
) -> Result<AggregateBucklingDesignPoint, SolidError> {
    let physics = element_scale_aggregate_sensitivities_cancellable(
        problem,
        scales,
        options.modes,
        options.design.spectral_steps,
        options.sharpness,
        cancelled,
    )?;
    let sample = design_sample(
        weights,
        scales,
        physics.value(),
        physics.gradient(),
        options.design,
    )?;
    Ok(BucklingDesignPoint {
        meets_floor: physics.value() >= options.design.minimum_load_factor,
        physics,
        sample,
    })
}

/// Resumable modulus allocation through repeated or crossing buckling modes.
///
/// The inequality is 1 - KS/minimum_load_factor <= 0. Since KS is a lower
/// bound on the smallest physical critical load, satisfying this inequality
/// conservatively satisfies the physical requirement. Finite negative KS
/// values remain valid infeasible starts; they are never physical eigenvalues.
///
/// Repeated modes strictly inside the fixed packet are smooth and admissible.
/// A cluster cut at its outer boundary or a shortfall of positive modes refuses
/// the trial. Increase the packet size in a new study if that boundary matters;
/// the objective and its derivatives never switch smoothing policy mid-run.
///
/// Every evaluation re-equilibrates the continuum and uses one aggregate
/// adjoint, independent of packet size. Original P1/Q1, Standard/B-bar and
/// prescribed-displacement physics remain with fs-solid. Accepted state,
/// cancellation, work accounting and best-feasible retention use the same
/// projected AL owner as BucklingDesignStudy. Geometry and loading functions
/// must remain fixed throughout a study.
pub struct AggregateBucklingDesignStudy<'problem, 'mesh> {
    problem: &'problem LinearProblem<'mesh>,
    weights: Vec<f64>,
    options: AggregateBucklingDesignOptions,
    engine: DesignState<ElementScaleAggregateSensitivity>,
}

impl<'problem, 'mesh> AggregateBucklingDesignStudy<'problem, 'mesh> {
    /// Bind a finite positive modulus box and cost to each continuum element.
    ///
    /// # Errors
    /// Invalid policies, incomplete mode packets and native numerical failures
    /// refuse construction without publishing a partially evaluated design.
    pub fn new(
        problem: &'problem LinearProblem<'mesh>,
        initial: &[f64],
        lower: &[f64],
        upper: &[f64],
        weights: &[f64],
        options: AggregateBucklingDesignOptions,
    ) -> Result<Self, ProjectedAlError<SolidError>> {
        Self::new_cancellable(problem, initial, lower, upper, weights, options, &|| false)
    }

    /// Construct using the caller's cancellation callback through native solves.
    ///
    /// # Errors
    /// As new(), plus cancellation before or during native numerical work.
    #[allow(clippy::too_many_arguments)]
    pub fn new_cancellable(
        problem: &'problem LinearProblem<'mesh>,
        initial: &[f64],
        lower: &[f64],
        upper: &[f64],
        weights: &[f64],
        options: AggregateBucklingDesignOptions,
        cancelled: &impl Fn() -> bool,
    ) -> Result<Self, ProjectedAlError<SolidError>> {
        validate_design(problem, initial, lower, upper, weights, options.design)?;
        if options.modes == 0 || !options.sharpness.is_finite() || options.sharpness <= 0.0 {
            return Err(ProjectedAlError::Invalid(
                "positive complete packet size and finite positive KS sharpness required",
            ));
        }
        let engine = initialize(
            initial,
            lower,
            upper,
            options.design.optimizer,
            &mut |scales| evaluate_aggregate(problem, weights, scales, options, cancelled),
            cancelled,
        )?;
        Ok(Self {
            problem,
            weights: weights.to_vec(),
            options,
            engine,
        })
    }

    /// Original element modulus multipliers at the last accepted step.
    #[must_use]
    pub fn element_scales(&self) -> &[f64] {
        self.engine.state.point()
    }

    /// Paired physical modes, KS value, total gradient and optimizer sample.
    #[must_use]
    pub fn accepted(&self) -> &AggregateBucklingDesignPoint {
        &self.engine.accepted
    }

    /// Least-cost accepted design whose KS value strictly meets the load floor.
    ///
    /// A physical first mode above the floor alone does not qualify if the
    /// selected conservative envelope still violates the constraint.
    #[must_use]
    pub fn best_feasible(&self) -> Option<&AggregateBucklingDesignPoint> {
        self.engine.best_feasible.as_ref()
    }

    /// Cumulative optimizer work, including rejected and interrupted evaluations.
    #[must_use]
    pub fn work(&self) -> ProjectedAlWork {
        self.engine.state.work()
    }

    /// Advance at most this many accepted steps with a fixed spectral policy.
    ///
    /// # Errors
    /// Native failures and unresolved packet boundaries preserve the last fully
    /// accepted point; no frozen or partial trial is installed.
    pub fn run(
        &mut self,
        additional_steps: usize,
    ) -> Result<ProjectedAlReport, ProjectedAlError<SolidError>> {
        self.run_cancellable(additional_steps, &|| false)
    }

    /// Advance with cancellation through equilibrium, spectral and adjoint work.
    /// A zero-step call reports cached state without a new physics evaluation.
    ///
    /// # Errors
    /// As run(), plus cancellation from the native numerical owners.
    pub fn run_cancellable(
        &mut self,
        additional_steps: usize,
        cancelled: &impl Fn() -> bool,
    ) -> Result<ProjectedAlReport, ProjectedAlError<SolidError>> {
        let (problem, weights, options) = (self.problem, &self.weights, self.options);
        self.engine.advance(
            additional_steps,
            cancelled,
            ElementScaleAggregateSensitivity::element_scales,
            &mut |scales| evaluate_aggregate(problem, weights, scales, options, cancelled),
        )
    }
}
