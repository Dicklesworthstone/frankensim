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
//! The first mode must remain separated: unresolved crossings STOP the smooth
//! design session, without frozen prestress or a substituted branch gradient.
use fs_ascent::projected_al::{
    ProjectedAlError, ProjectedAlOptions, ProjectedAlReport, ProjectedAlSample, ProjectedAlState,
    ProjectedAlStop, ProjectedAlWork,
};
use fs_solid::SolidError;
use fs_solid::linear::LinearProblem;
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
pub struct BucklingDesignPoint {
    physics: ElementScaleSensitivity,
    sample: ProjectedAlSample,
    meets_floor: bool,
}
impl BucklingDesignPoint {
    /// Complete equilibrium-dependent mode, total derivatives and residuals.
    #[must_use]
    pub fn physics(&self) -> &ElementScaleSensitivity {
        &self.physics
    }
    /// Cost and normalized constraint with their exact discrete derivatives.
    /// The inequality is 1 - critical_factor/minimum_load_factor <= 0.
    #[must_use]
    pub fn sample(&self) -> &ProjectedAlSample {
        &self.sample
    }
    /// True only when the returned critical factor meets the original floor;
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
    let normalized = factor / options.minimum_load_factor;
    if !normalized.is_finite() || normalized <= 0.0 {
        return Err(invalid(
            "buckling design load normalization is not representable",
        ));
    }
    let mut constraint_gradient = Vec::with_capacity(weights.len());
    for &derivative in &physics.gradients()[0] {
        let normalized = -derivative / options.minimum_load_factor;
        if !normalized.is_finite() || (derivative != 0.0 && normalized == 0.0) {
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
    Ok(BucklingDesignPoint {
        meets_floor: factor >= options.minimum_load_factor,
        physics,
        sample: ProjectedAlSample {
            objective: cost,
            gradient: weights.to_vec(),
            constraint: 1.0 - normalized,
            constraint_gradient,
        },
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
    state: ProjectedAlState,
    accepted: BucklingDesignPoint,
    best_feasible: Option<BucklingDesignPoint>,
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
        let mut initial_physics = None;
        let state = ProjectedAlState::try_new(
            initial,
            lower,
            upper,
            options.optimizer,
            &mut |scales| {
                let point = evaluate(problem, weights, scales, options, cancelled)?;
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
        let accepted =
            initial_physics.expect("accepted initial sample has matching native physics");
        let best_feasible = accepted.meets_floor.then(|| accepted.clone());
        Ok(Self {
            problem,
            weights: weights.to_vec(),
            options,
            state,
            accepted,
            best_feasible,
        })
    }

    /// The last accepted original element multipliers.
    #[must_use]
    pub fn element_scales(&self) -> &[f64] {
        self.state.point()
    }
    /// Physics and full derivative at exactly the current optimizer point.
    #[must_use]
    pub fn accepted(&self) -> &BucklingDesignPoint {
        &self.accepted
    }
    /// Lowest-cost accepted design satisfying the ORIGINAL positive load floor.
    #[must_use]
    pub fn best_feasible(&self) -> Option<&BucklingDesignPoint> {
        self.best_feasible.as_ref()
    }
    /// Cumulative work includes rejected and interrupted physics attempts.
    #[must_use]
    pub fn work(&self) -> ProjectedAlWork {
        self.state.work()
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
        let mut remaining = additional_steps;
        loop {
            let before = self.state.work().iterations;
            let mut last = None;
            let result = self.state.try_run(
                remaining.min(1),
                &mut |scales| {
                    // No error becomes an unavailable or frozen-gradient sample.
                    let point =
                        evaluate(self.problem, &self.weights, scales, self.options, cancelled)?;
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
                assert_eq!(accepted.physics.element_scales(), self.state.point());
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
