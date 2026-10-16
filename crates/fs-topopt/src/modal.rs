//! Matrix-free structural modes and eigenfrequency design sensitivities.
//!
//! This is the L4 domain adapter to L1 fs-spectral, not another eigensolver.
//! The existing dense eigenfrequency API remains an independent small-model
//! reference. No dense fallback, mass lumping, or synthetic convergence status.

use fs_solver::LinearOp;
use fs_sparse::precond::Precond;
use fs_spectral::generalized::{GeneralizedError, GeneralizedOp, GeneralizedOptions, GeneralizedState};

use crate::control::{EvaluationStop, SolveControl, SolveWork};
use crate::eigenfreq::{mass_interp, mass_interp_derivative, smooth_min};
use crate::elasticity::DensityElasticity;
use crate::filter::heaviside_derivative;
use crate::pipeline::DesignPipeline;

/// Domain failures and spectral numerical refusals, including original budgets.
pub type MatrixFreeEigenError = GeneralizedError<EvaluationStop>;

/// Explicit settings for the matrix-free modal calculation.
#[derive(Debug, Clone, Copy)]
pub struct MatrixFreeEigenOptions {
    /// Number of low modes requested; a repeated cluster should be kept whole.
    pub count: usize,
    /// Additional search vectors beyond count (limited by the free dimension).
    pub oversampling: usize,
    /// Maximum block inverse steps; zero permits only an already-converged seed.
    pub max_iterations: usize,
    /// Recomputed relative pencil-residual threshold.
    pub relative_tolerance: f64,
    /// CG relative tolerance, normally much tighter than the modal tolerance.
    /// Per-solve and total iteration limits come from the shared SolveControl.
    pub linear_tolerance: f64,
    /// Explicit, deterministic search-block initialization seed.
    pub seed: u64,
}

impl Default for MatrixFreeEigenOptions {
    fn default() -> Self {
        Self { count: 3, oversampling: 4, max_iterations: 128,
            relative_tolerance: 1e-8, linear_tolerance: 1e-11, seed: 42 }
    }
}

/// Full, residual-qualified modal fields for the original mesh.
#[derive(Debug, Clone)]
pub struct MatrixFreeEigenReport {
    /// Eigenvalues lambda = omega^2, in ascending order.
    pub values: Vec<f64>,
    /// Consistent-mass-normalized full displacement vectors. Fixed entries are
    /// exactly zero; constraint identity rows never participate as fake modes.
    pub modes: Vec<Vec<f64>>,
    /// ||K phi - lambda M phi|| / (||K phi|| + |lambda| ||M phi||).
    pub relative_residuals: Vec<f64>,
    /// Accepted block inverse steps.
    pub iterations: usize,
    /// Cumulative actual work in the caller's shared control, not just this call.
    pub work: SolveWork,
}

/// Model and aggregation settings, separate from eigensolver tolerances.
#[derive(Debug, Clone, Copy)]
pub struct EigenfrequencyObjectiveOptions {
    /// Modal solver settings, including the number of aggregated modes.
    pub solver: MatrixFreeEigenOptions,
    /// Positive smooth-min sharpness, in inverse eigenvalue units.
    pub aggregation_beta: f64,
    /// Positive full-material mass density. With SI geometry and elasticity,
    /// this is kg/m^3 and the modal eigenvalues have units s^-2.
    pub reference_density: f64,
}

impl Default for EigenfrequencyObjectiveOptions {
    fn default() -> Self {
        Self { solver: MatrixFreeEigenOptions::default(), aggregation_beta: 1.0, reference_density: 1.0 }
    }
}

/// Smooth low-mode objective and its raw-density derivative at the SAME design.
#[derive(Debug, Clone)]
pub struct MatrixFreeEigenObjective {
    /// Unnormalized smooth minimum, a lower bound on the computed low values.
    /// It may be negative; it is an optimization objective, not a squared mode.
    pub aggregate: f64,
    /// Derivative through consistent mass, SIMP, projection, and filter transpose.
    pub gradient: Vec<f64>,
    /// Actual modal fields and solve diagnostics used for this derivative.
    pub modal: MatrixFreeEigenReport,
}

struct Diagonal(Vec<f64>);
impl Precond for Diagonal {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for ((&r, z), &d) in r.iter().zip(z).zip(&self.0) { *z = r / d; }
    }
}

struct ElasticPencil<'a, 'b, 'c> {
    elasticity: &'a DensityElasticity,
    mass: &'a [f64],
    free: Vec<usize>,
    diagonal: Diagonal,
    tolerance: f64,
    control: &'b mut SolveControl<'c>,
}

impl ElasticPencil<'_, '_, '_> {
    fn scatter(&self, reduced: &[f64]) -> Vec<f64> {
        let mut full = vec![0.0; self.elasticity.n()];
        for (&dof, &value) in self.free.iter().zip(reduced) { full[dof] = value; }
        full
    }
    fn gather(&self, full: &[f64], reduced: &mut [f64]) {
        for (&dof, value) in self.free.iter().zip(reduced) { *value = full[dof]; }
    }
}

impl GeneralizedOp for ElasticPencil<'_, '_, '_> {
    type Error = EvaluationStop;
    fn dim(&self) -> usize { self.free.len() }
    fn apply_a(&mut self, x: &[f64], y: &mut [f64]) -> Result<(), Self::Error> {
        self.control.checkpoint("eigen-stiffness")?;
        let full = self.scatter(x);
        let mut output = vec![0.0; self.elasticity.n()];
        self.elasticity.apply(&full, &mut output);
        self.gather(&output, y);
        self.control.checkpoint("eigen-stiffness")
    }
    fn apply_b(&mut self, x: &[f64], y: &mut [f64]) -> Result<(), Self::Error> {
        self.control.checkpoint("eigen-mass")?;
        let full = self.scatter(x);
        let mut output = vec![0.0; self.elasticity.n()];
        self.elasticity.apply_mass(self.mass, &full, &mut output);
        self.gather(&output, y);
        self.control.checkpoint("eigen-mass")
    }
    fn solve_a(&mut self, rhs: &[f64], y: &mut [f64]) -> Result<(), Self::Error> {
        let full_rhs = self.scatter(rhs);
        let solution = self.control.solve_preconditioned(self.elasticity, &self.diagonal,
            &full_rhs, self.tolerance, usize::MAX, "eigen-elasticity")?;
        self.gather(&solution, y);
        Ok(())
    }
    fn checkpoint(&mut self) -> Result<(), Self::Error> { self.control.checkpoint("eigen") }
}

/// Compute low elastic modes using the original consistent element mass.
/// `mass_densities` contains one physical mass weight per cell. Zero weights
/// are allowed only when every free node still has strictly positive mass.
/// Supports must eliminate rigid modes; symmetry/SPD are model preconditions.
///
/// All returned modes satisfy the explicitly recomputed residual threshold.
/// This is not a certified proof of spectral completeness or continuum error.
pub fn controlled_matrix_free_eigenpairs(
    elasticity: &DensityElasticity,
    mass_densities: &[f64],
    options: MatrixFreeEigenOptions,
    control: &mut SolveControl<'_>,
) -> Result<MatrixFreeEigenReport, MatrixFreeEigenError> {
    if mass_densities.len() != elasticity.cells()
        || mass_densities.iter().any(|r| !r.is_finite() || *r < 0.0) {
        return Err(GeneralizedError::InvalidInput("one finite nonnegative mass density per cell is required"));
    }
    if elasticity.moduli.len() != elasticity.cells()
        || elasticity.moduli.iter().any(|e| !e.is_finite() || *e <= 0.0) {
        return Err(GeneralizedError::InvalidInput("stiffness weights must be finite, positive, and match the cells"));
    }
    if !options.linear_tolerance.is_finite() || options.linear_tolerance <= 0.0
        || options.linear_tolerance >= 1.0 {
        return Err(GeneralizedError::InvalidInput("linear tolerance must lie in (0, 1)"));
    }
    control.checkpoint("eigen-setup").map_err(GeneralizedError::Operator)?;
    let free: Vec<usize> = elasticity.free().iter().enumerate()
        .filter_map(|(d, is_free)| is_free.then_some(d)).collect();
    if options.count == 0 || options.count > free.len() {
        return Err(GeneralizedError::InvalidInput("mode count must lie in 1..=free dofs"));
    }
    let block_size = options.count.checked_add(options.oversampling)
        .ok_or(GeneralizedError::InvalidInput("search-block width overflows"))?.min(free.len());
    let diagonal = elasticity.stiffness_diagonal();
    let mass_diagonal = elasticity.mass_diagonal(mass_densities);
    if free.iter().any(|&d| !diagonal[d].is_finite() || diagonal[d] <= 0.0
        || !mass_diagonal[d].is_finite() || mass_diagonal[d] <= 0.0) {
        return Err(GeneralizedError::InvalidInput("every free dof requires positive finite stiffness and consistent mass"));
    }
    let mut op = ElasticPencil { elasticity, mass: mass_densities, free,
        diagonal: Diagonal(diagonal), tolerance: options.linear_tolerance, control };
    let settings = GeneralizedOptions { count: options.count, block_size,
        relative_tolerance: options.relative_tolerance, seed: options.seed };
    let mut state = GeneralizedState::new(&mut op, settings, None)?;
    let report = state.run(&mut op, options.max_iterations)?;
    let modes = report.vectors.iter().map(|v| op.scatter(v)).collect();
    op.control.checkpoint("eigen-publish").map_err(GeneralizedError::Operator)?;
    Ok(MatrixFreeEigenReport { values: report.values, modes,
        relative_residuals: report.relative_residuals, iterations: report.iterations,
        work: op.control.work() })
}

/// Evaluate a frequency objective and its exact discrete-design chain rule
/// using residual-qualified matrix-free modes. Includes the derivative of
/// physical mass density (including the low-density rho^6 interpolation).
/// The whole repeated cluster must be included for a differentiable aggregate;
/// a fixed truncated cluster is not smooth at a crossing with an omitted mode.
///
/// Filter/mesh compatibility retains DesignPipeline's constructor contracts.
/// On a returned error, restore the previous elasticity moduli and do not
/// publish a partial gradient. Work spent on the failed evaluation remains
/// charged to the shared control. On success moduli match the evaluated design.
pub fn controlled_matrix_free_eigenfrequency_objective(
    pipeline: &DesignPipeline,
    elasticity: &mut DensityElasticity,
    rho: &[f64],
    options: EigenfrequencyObjectiveOptions,
    control: &mut SolveControl<'_>,
) -> Result<MatrixFreeEigenObjective, MatrixFreeEigenError> {
    if rho.len() != elasticity.cells() || rho.is_empty()
        || rho.iter().any(|r| !r.is_finite() || !(0.0..=1.0).contains(r)) {
        return Err(GeneralizedError::InvalidInput("raw densities must match the cells and lie in [0, 1]"));
    }
    if !options.reference_density.is_finite() || options.reference_density <= 0.0
        || !options.aggregation_beta.is_finite() || options.aggregation_beta <= 0.0 {
        return Err(GeneralizedError::InvalidInput("reference density and aggregation beta must be finite and positive"));
    }
    let previous_moduli = elasticity.moduli.clone();
    let result = (|| {
        let (filtered, projected, moduli) = pipeline.try_forward(rho, control)
            .map_err(GeneralizedError::Operator)?;
        elasticity.moduli = moduli;
        let mass: Vec<f64> = projected.iter()
            .map(|&r| options.reference_density * mass_interp(r)).collect();
        let mut modal = controlled_matrix_free_eigenpairs(elasticity, &mass, options.solver, control)?;
        let (aggregate, weights) = smooth_min(&modal.values, options.aggregation_beta);
        let p = &pipeline.params;
        let mut local = vec![0.0; rho.len()];
        for ((&lambda, mode), &weight) in modal.values.iter().zip(&modal.modes).zip(&weights) {
            control.checkpoint("eigen-sensitivity").map_err(GeneralizedError::Operator)?;
            let strain = elasticity.cell_energies(mode);
            let kinetic = elasticity.cell_kinetic(mode);
            for cell in 0..rho.len() {
                let rb = projected[cell].clamp(0.0, 1.0);
                let stiffness_slope = (1.0 - p.e_min) * p.penal
                    * fs_math::det::pow(rb.max(1e-12), p.penal - 1.0);
                let mass_slope = options.reference_density * mass_interp_derivative(rb);
                let derivative = stiffness_slope.mul_add(strain[cell], -lambda * mass_slope * kinetic[cell]);
                local[cell] = weight.mul_add(derivative, local[cell]);
            }
        }
        for (derivative, &r) in local.iter_mut().zip(&filtered) {
            *derivative *= heaviside_derivative(r, p.beta, p.eta);
        }
        if !aggregate.is_finite() || local.iter().any(|d| !d.is_finite()) {
            return Err(GeneralizedError::Breakdown("non-finite frequency objective or sensitivity"));
        }
        let gradient = pipeline.filter.try_apply_transpose(&local, control)
            .map_err(GeneralizedError::Operator)?;
        if gradient.iter().any(|d| !d.is_finite()) {
            return Err(GeneralizedError::Breakdown("non-finite frequency design gradient"));
        }
        control.checkpoint("eigen-objective-publish").map_err(GeneralizedError::Operator)?;
        modal.work = control.work();
        Ok(MatrixFreeEigenObjective { aggregate, gradient, modal })
    })();
    if result.is_err() { elasticity.moduli = previous_moduli; }
    result
}
