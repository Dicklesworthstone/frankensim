//! Material-aware free vibrations and eigenfrequency design on the native 3-D grid.
//!
//! The stiffness, consistent physical mass and filter retain their original
//! owners. This is an adapter to the shared matrix-free generalized eigensolver,
//! not a second modal solver or a lumped-mass approximation. Strong clamp rows
//! are removed from the spectral space; adaptive hanging constraints already
//! belong to both native operators.

pub mod design;
pub mod refinement;
pub mod refit;

use fs_cutfem::elastic3::mass::MaterialMass3;
use fs_sparse::precond::Precond;
use fs_spectral::generalized::{
    GeneralizedError, GeneralizedOp, GeneralizedOptions, GeneralizedState,
};

use super::operator::native_controlled;
use super::{CutDensityStudy3, Sdf3Elasticity, checked_solve_preconditioned, failure};
use crate::eigenfreq::{mass_interp, mass_interp_derivative, smooth_min};
use crate::modal::{
    MatrixFreeEigenError, MatrixFreeEigenObjective, MatrixFreeEigenOptions, MatrixFreeEigenReport,
};
use crate::{EvaluationStop, SolveControl};

/// Solver controls and smooth-minimum scale for a 3-D frequency objective.
#[derive(Debug, Clone, Copy)]
pub struct Sdf3EigenfrequencyOptions3 {
    /// Low modes to aggregate; a further mode is solved when needed to check
    /// that the selected cluster does not end at an unresolved crossing.
    pub solver: MatrixFreeEigenOptions,
    /// Positive sharpness in inverse eigenvalue units (squared seconds in SI).
    pub aggregation_beta: f64,
}

impl Default for Sdf3EigenfrequencyOptions3 {
    fn default() -> Self {
        Self {
            solver: MatrixFreeEigenOptions::default(),
            aggregation_beta: 1.0,
        }
    }
}

/// Frequency and physical inventory at one exact filtered/projected design.
#[derive(Debug, Clone)]
pub struct Sdf3EigenfrequencyEvaluation3 {
    pub objective: MatrixFreeEigenObjective,
    pub projected_rho: Vec<f64>,
    pub volume_fraction: f64,
    pub volume_gradient: Vec<f64>,
    /// Linear physical occupancy times native density, in kg for SI geometry.
    /// This inventory does not use the low-density modal mass interpolation.
    pub physical_mass_kg: f64,
    pub physical_mass_gradient: Vec<f64>,
    /// First excluded numerical mode, when a truncated cluster was used.
    /// Its gap is an admission check, not a spectral completeness certificate.
    pub next_eigenvalue: Option<f64>,
}

const MASS_SWITCH: &str =
    "sdf3 eigenfrequency mass interpolation is nondifferentiable at rho=0.1";
const CLUSTER_CROSSING: &str =
    "sdf3 eigenfrequency selected cluster has no resolved outer gap";

/// Trial optimizers may backtrack over these known nonsmooth points without
/// discarding cancellation, linear budgets or a genuine spectral breakdown.
pub(super) fn nondifferentiable(error: &MatrixFreeEigenError) -> bool {
    match error {
        GeneralizedError::InvalidInput(message) | GeneralizedError::Breakdown(message) => {
            *message == MASS_SWITCH || *message == CLUSTER_CROSSING
        }
        _ => false,
    }
}

struct ElasticPencil<'a, 'b, 'c, O: Sdf3Elasticity, P: Precond> {
    elasticity: &'a O,
    mass: &'a MaterialMass3,
    occupancy: &'a [f64],
    free: Vec<usize>,
    prepared: P,
    full: Vec<f64>,
    output: Vec<f64>,
    tolerance: f64,
    control: &'b mut SolveControl<'c>,
}

impl<O: Sdf3Elasticity, P: Precond> ElasticPencil<'_, '_, '_, O, P> {
    fn scatter(&mut self, reduced: &[f64]) {
        // Fixed entries stay zero for the lifetime of this scratch vector.
        for (&dof, &value) in self.free.iter().zip(reduced) {
            self.full[dof] = value;
        }
    }

    fn gather(&self, full: &[f64], reduced: &mut [f64]) {
        for (&dof, value) in self.free.iter().zip(reduced) {
            *value = full[dof];
        }
    }
}

impl<O: Sdf3Elasticity, P: Precond> GeneralizedOp for ElasticPencil<'_, '_, '_, O, P> {
    type Error = EvaluationStop;

    fn dim(&self) -> usize {
        self.free.len()
    }

    fn apply_a(&mut self, x: &[f64], y: &mut [f64]) -> Result<(), Self::Error> {
        self.control.checkpoint("sdf3-eigen-stiffness")?;
        self.scatter(x);
        self.elasticity.apply(&self.full, &mut self.output);
        self.gather(&self.output, y);
        self.control.checkpoint("sdf3-eigen-stiffness")
    }

    fn apply_b(&mut self, x: &[f64], y: &mut [f64]) -> Result<(), Self::Error> {
        self.scatter(x);
        native_controlled(self.control, "sdf3-eigen-mass", |poll| {
            self.mass
                .apply_into(self.occupancy, &self.full, &mut self.output, poll)
        })?;
        self.gather(&self.output, y);
        Ok(())
    }

    fn solve_a(&mut self, rhs: &[f64], y: &mut [f64]) -> Result<(), Self::Error> {
        self.scatter(rhs);
        let solution = checked_solve_preconditioned(
            self.elasticity,
            &self.prepared,
            &self.full,
            self.tolerance,
            "sdf3-eigen-elasticity",
            self.control,
        )?;
        self.gather(&solution, y);
        Ok(())
    }

    fn checkpoint(&mut self) -> Result<(), Self::Error> {
        self.control.checkpoint("sdf3-eigen")
    }
}

/// Low modes of the current native stiffness and an explicitly weighted
/// material mass. The prepared mass is reusable across density iterations and
/// is validated against this exact geometry, material and constraint owner.
///
/// `occupancy` contains the mass multiplier per active cell; no stiffness
/// scale, ersatz density or mass lumping is substituted. Every free DOF must
/// retain positive mass. Supports must eliminate rigid motion and the native
/// stiffness must be SPD on that space. Each inner solve passes the shared
/// true-residual gate and consumes the caller's original cumulative budget.
///
/// These are unforced vibrations about the undeformed reference state.
/// External forces, gravity and prescribed nonzero support motion do not
/// enter this pencil; prestress, geometric stiffness and damping are absent.
/// Mode normalization is arbitrary physical amplitude, so constitutive strain
/// ceilings are not imposed on mass-normalized eigenvectors. Residuals and
/// Ritz values are numerical estimates, not continuum or completeness bounds.
pub fn controlled_sdf3_eigenpairs<O: Sdf3Elasticity>(
    elasticity: &O,
    mass: &MaterialMass3,
    occupancy: &[f64],
    options: MatrixFreeEigenOptions,
    control: &mut SolveControl<'_>,
) -> Result<MatrixFreeEigenReport, MatrixFreeEigenError> {
    control
        .checkpoint("sdf3-eigen-setup")
        .map_err(GeneralizedError::Operator)?;
    elasticity.validate_material_mass(mass).map_err(|_| {
        GeneralizedError::InvalidInput("sdf3 modal mass does not match the native operator")
    })?;
    if occupancy.len() != elasticity.cells()
        || occupancy
            .iter()
            .any(|r| !r.is_finite() || !(0.0..=1.0).contains(r))
    {
        return Err(GeneralizedError::InvalidInput(
            "one finite mass occupancy in [0,1] per active cell is required",
        ));
    }
    if elasticity.scales().len() != elasticity.cells()
        || elasticity
            .scales()
            .iter()
            .any(|e| !e.is_finite() || *e <= 0.0)
    {
        return Err(GeneralizedError::InvalidInput(
            "stiffness weights must be finite and positive",
        ));
    }
    if elasticity.fixed().len().checked_mul(3) != Some(elasticity.n()) {
        return Err(GeneralizedError::InvalidInput(
            "native node constraints do not match the displacement space",
        ));
    }
    if !options.linear_tolerance.is_finite()
        || options.linear_tolerance <= 0.0
        || options.linear_tolerance >= 1.0
    {
        return Err(GeneralizedError::InvalidInput(
            "linear tolerance must lie in (0,1)",
        ));
    }
    let free: Vec<usize> = elasticity
        .fixed()
        .iter()
        .enumerate()
        .filter(|(_, fixed)| !**fixed)
        .flat_map(|(node, _)| [3 * node, 3 * node + 1, 3 * node + 2])
        .collect();
    if options.count == 0 || options.count > free.len() {
        return Err(GeneralizedError::InvalidInput(
            "mode count must lie in 1..=free dofs",
        ));
    }
    let block_size = options
        .count
        .checked_add(options.oversampling)
        .ok_or(GeneralizedError::InvalidInput(
            "search-block width overflows",
        ))?
        .min(free.len());
    if block_size > 64 {
        return Err(GeneralizedError::InvalidInput(
            "generalized search block exceeds 64 modes",
        ));
    }
    let diagonal = native_controlled(control, "sdf3-eigen-mass-admission", |poll| {
        mass.diagonal(occupancy, poll)
    })
    .map_err(GeneralizedError::Operator)?;
    if free
        .iter()
        .any(|&d| !diagonal[d].is_finite() || diagonal[d] <= 0.0)
    {
        return Err(GeneralizedError::InvalidInput(
            "every free dof requires positive finite consistent mass",
        ));
    }
    let prepared = elasticity
        .prepare_elasticity(control)
        .map_err(GeneralizedError::Operator)?;
    let mut op = ElasticPencil {
        elasticity,
        mass,
        occupancy,
        free,
        prepared,
        full: vec![0.0; elasticity.n()],
        output: vec![0.0; elasticity.n()],
        tolerance: options.linear_tolerance,
        control,
    };
    let settings = GeneralizedOptions {
        count: options.count,
        block_size,
        relative_tolerance: options.relative_tolerance,
        seed: options.seed,
    };
    let mut state = GeneralizedState::new(&mut op, settings, None)?;
    let report = state.run(&mut op, options.max_iterations)?;
    let mut modes = Vec::with_capacity(report.vectors.len());
    for reduced in &report.vectors {
        op.control
            .checkpoint("sdf3-eigen-fields")
            .map_err(GeneralizedError::Operator)?;
        let mut full = vec![0.0; elasticity.n()];
        for (&dof, &value) in op.free.iter().zip(reduced) {
            full[dof] = value;
        }
        modes.push(full);
    }
    op.control
        .checkpoint("sdf3-eigen-publish")
        .map_err(GeneralizedError::Operator)?;
    Ok(MatrixFreeEigenReport {
        values: report.values,
        modes,
        relative_residuals: report.relative_residuals,
        iterations: report.iterations,
        work: op.control.work(),
    })
}

/// One selected cluster and the same outer-gap admission on every native grid.
/// This helper evaluates modes only; assessment does not request a design
/// gradient at a transferred physical field or impose a strain amplitude.
pub(super) fn selected_modes<O: Sdf3Elasticity>(
    elasticity: &O,
    mass: &MaterialMass3,
    occupancy: &[f64],
    options: Sdf3EigenfrequencyOptions3,
    control: &mut SolveControl<'_>,
) -> Result<(MatrixFreeEigenReport, Option<f64>), MatrixFreeEigenError> {
    let free = elasticity
        .fixed()
        .iter()
        .filter(|fixed| !**fixed)
        .count()
        .checked_mul(3)
        .ok_or(GeneralizedError::InvalidInput("free dimension overflows"))?;
    let count = options.solver.count;
    if count == 0 || count > free {
        return Err(GeneralizedError::InvalidInput(
            "mode count must lie in 1..=free dofs",
        ));
    }
    let guarded = count < free;
    let mut solver = options.solver;
    if guarded {
        if count >= 64 {
            return Err(GeneralizedError::InvalidInput(
                "selected cluster requires an additional guard mode within the 64-mode cap",
            ));
        }
        solver.count += 1;
        solver.oversampling = solver.oversampling.saturating_sub(1);
    }
    let mut modal = controlled_sdf3_eigenpairs(elasticity, mass, occupancy, solver, control)?;
    let next = if guarded {
        let last = modal.values[count - 1];
        let next = modal.values[count];
        let scale = last.abs().max(next.abs());
        let gap = next - last;
        let tolerance = (64.0 * f64::EPSILON + 2.0 * solver.relative_tolerance) * scale;
        if !gap.is_finite() || gap <= tolerance {
            return Err(GeneralizedError::Breakdown(CLUSTER_CROSSING));
        }
        modal.values.truncate(count);
        modal.modes.truncate(count);
        modal.relative_residuals.truncate(count);
        Some(next)
    } else {
        None
    };
    Ok((modal, next))
}

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Source-material modes, smooth minimum, physical inventories and exact
    /// discrete design derivatives through the SAME filter/projection/SIMP map.
    /// Native mass uses the shared low-density `mass_interp`; kg inventory uses
    /// linear occupancy. Protected solid/void regions retain their actual
    /// physical map and zero local design slope in both mechanisms.
    ///
    /// The selected modal cluster must be separated from its first excluded
    /// mode by the declared numerical tolerance. One extra mode checks this
    /// boundary, consuming one oversampling slot where available. Unresolved
    /// crossings and the variable-cell mass-interpolation switch withhold the
    /// gradient. This is not a certified proof of spectral completeness.
    ///
    /// Success and every returned error restore incoming stiffness scales.
    /// The prepared native mass remains reusable. All consumed work, including
    /// filter/inventory pullbacks, remains charged to the original control.
    pub fn evaluate_eigenfrequency(
        &mut self,
        mass: &MaterialMass3,
        rho: &[f64],
        options: Sdf3EigenfrequencyOptions3,
        control: &mut SolveControl<'_>,
    ) -> Result<Sdf3EigenfrequencyEvaluation3, MatrixFreeEigenError> {
        control
            .checkpoint("sdf3-eigen-objective")
            .map_err(GeneralizedError::Operator)?;
        if rho.len() != self.cells()
            || rho.is_empty()
            || rho
                .iter()
                .any(|r| !r.is_finite() || !(0.0..=1.0).contains(r))
        {
            return Err(GeneralizedError::InvalidInput(
                "raw densities must match the cells and lie in [0,1]",
            ));
        }
        if !options.aggregation_beta.is_finite() || options.aggregation_beta <= 0.0 {
            return Err(GeneralizedError::InvalidInput(
                "aggregation beta must be finite and positive",
            ));
        }
        self.operator.validate_material_mass(mass).map_err(|_| {
            GeneralizedError::InvalidInput("sdf3 modal mass does not match the native operator")
        })?;
        let design = self
            .design(rho, control)
            .map_err(GeneralizedError::Operator)?;
        if design
            .projected
            .iter()
            .zip(&design.slope)
            .any(|(&r, &slope)| r == 0.1 && slope != 0.0)
        {
            return Err(GeneralizedError::InvalidInput(MASS_SWITCH));
        }
        let occupancy: Vec<f64> = design.projected.iter().map(|&r| mass_interp(r)).collect();
        let previous = self.operator.scales().to_vec();
        self.operator
            .set_scales(&design.scales)
            .map_err(|_| GeneralizedError::Operator(failure("sdf3-eigen-scales")))?;
        let result = (|| {
            let (mut modal, next_eigenvalue) =
                selected_modes(&self.operator, mass, &occupancy, options, control)?;
            let (aggregate, weights) = smooth_min(&modal.values, options.aggregation_beta);
            let mut stiffness = vec![0.0; self.cells()];
            let mut inertia = vec![0.0; self.cells()];
            for ((&lambda, mode), &weight) in
                modal.values.iter().zip(&modal.modes).zip(&weights)
            {
                let strain =
                    native_controlled(control, "sdf3-eigen-stiffness-derivative", |poll| {
                        self.operator.scale_bilinear_forms(mode, mode, poll)
                    })
                    .map_err(GeneralizedError::Operator)?;
                let kinetic =
                    native_controlled(control, "sdf3-eigen-mass-derivative", |poll| {
                        mass.occupancy_pullback(mode, mode, poll)
                    })
                    .map_err(GeneralizedError::Operator)?;
                for cell in 0..self.cells() {
                    if cell % 64 == 0 {
                        control
                            .checkpoint("sdf3-eigen-gradient")
                            .map_err(GeneralizedError::Operator)?;
                    }
                    stiffness[cell] = weight.mul_add(strain[cell], stiffness[cell]);
                    let derivative = -lambda
                        * mass_interp_derivative(design.projected[cell])
                        * kinetic[cell];
                    inertia[cell] = weight.mul_add(derivative, inertia[cell]);
                }
            }
            if !aggregate.is_finite() {
                return Err(GeneralizedError::Breakdown(
                    "nonfinite sdf3 eigenfrequency aggregate",
                ));
            }
            let gradient = self
                .physical_gradient(&design, &stiffness, &inertia, control)
                .map_err(GeneralizedError::Operator)?;
            let local_volume: Vec<f64> = self
                .mass
                .iter()
                .zip(&design.slope)
                .map(|(m, slope)| m * slope)
                .collect();
            let volume_gradient = self
                .pullback(&local_volume, control)
                .map_err(GeneralizedError::Operator)?;
            let (physical_mass_kg, physical_mass_gradient) = self
                .mass_prepared(&design, control)
                .map_err(GeneralizedError::Operator)?;
            if !design.volume.is_finite() || volume_gradient.iter().any(|v| !v.is_finite()) {
                return Err(GeneralizedError::Breakdown(
                    "nonfinite sdf3 modal volume or gradient",
                ));
            }
            control
                .checkpoint("sdf3-eigen-objective-publish")
                .map_err(GeneralizedError::Operator)?;
            modal.work = control.work();
            Ok(Sdf3EigenfrequencyEvaluation3 {
                objective: MatrixFreeEigenObjective {
                    aggregate,
                    gradient,
                    modal,
                },
                projected_rho: design.projected,
                volume_fraction: design.volume,
                volume_gradient,
                physical_mass_kg,
                physical_mass_gradient,
                next_eigenvalue,
            })
        })();
        self.operator
            .set_scales(&previous)
            .expect("previous scales valid");
        result
    }
}
