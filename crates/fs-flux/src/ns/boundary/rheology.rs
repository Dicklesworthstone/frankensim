//! Checked generalized-Newtonian Stokes flow with open traction boundaries.

use super::{FluxParams, FluxSystem, StokesBoundary, StokesBoundaryError, StokesSolution};
use crate::ns::FluxRheologyError;

/// Controls for a frozen-load, steady constitutive solve.
#[derive(Clone, Copy, Debug)]
pub struct StokesRheologyControls {
    /// Maximum number of actual linearized flow solves.
    pub maximum_iterations: usize,
    /// Required relative velocity correction AND both nonlinear residuals.
    pub relative_tolerance: f64,
    /// Fraction of each correction used after the first admissible flow.
    /// Must lie in (0,1]. Under-relaxation can stabilize shear-thickening laws;
    /// it is not a convergence guarantee for an arbitrary constitutive law.
    pub relaxation: f64,
}

impl Default for StokesRheologyControls {
    fn default() -> Self {
        Self {
            maximum_iterations: 60,
            relative_tolerance: 1e-8,
            relaxation: 1.0,
        }
    }
}

/// Accepted steady flow at its own shear-dependent viscosity.
pub struct StokesRheologySolution {
    /// Pressure retains the admitted boundary datum. The field's two residuals
    /// are independently reassembled NONLINEAR endpoint residuals, tested
    /// against the nonlinear tolerance rather than the inner linear tolerance.
    pub solution: StokesSolution,
    /// Actual linearized flow solves, bounded by maximum_iterations.
    pub iterations: usize,
    /// Last relative velocity correction BEFORE under-relaxation.
    pub relative_update: f64,
    /// Larger of the finite original and row-scaled endpoint residuals.
    pub nonlinear_residual: f64,
}

impl<'m> FluxSystem<'m> {
    /// Solve -div(2 nu(gamma) D(u)) + grad(p) = f with mixed velocity and
    /// TOTAL-traction boundaries, where gamma = sqrt(2 D:D).
    ///
    /// The law receives (cell, physical point, current shear rate) and returns
    /// positive finite kinematic viscosity. Temperature and other captured
    /// material context must remain fixed throughout this solve; params.nu is
    /// unused. The supplied initial velocity is never changed, and its actual
    /// rates initialize the law without probing an assumed zero-shear state.
    ///
    /// Boundary samples are admitted once and reused by every nonlinear
    /// iterate. The first solved flow establishes compatible velocity moments
    /// and continuity. Subsequent corrections may be under-relaxed; relaxation
    /// never scales the reported correction or the physical residual.
    /// Acceptance requires both a small correction and the original AND
    /// row-scaled residuals reassembled at the returned velocity. Reassembly is
    /// deferred until that correction is small or the budget is exhausted.
    ///
    /// This is steady creeping flow: it adds neither convection/backflow
    /// policy, thermal evolution nor viscoelastic memory.
    ///
    /// # Errors
    /// Refuses invalid controls/state, inadmissible boundaries or coefficients,
    /// unconverged linear/nonlinear solves and unrepresentable arithmetic.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_stokes_mixed_rheology(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: StokesRheologyControls,
    ) -> Result<StokesRheologySolution, StokesBoundaryError> {
        self.solve_stokes_mixed_rheology_checked(
            params,
            law,
            force,
            boundaries,
            initial_velocity,
            controls,
            &|| true,
        )
    }

    // A material owner keeps its typed source error outside the scalar law.
    // Check it after each assembly, before linear work or final publication.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve_stokes_mixed_rheology_checked(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: StokesRheologyControls,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<StokesRheologySolution, StokesBoundaryError> {
        if initial_velocity.len() != self.n_u
            || initial_velocity.iter().any(|x| !x.is_finite())
            || controls.maximum_iterations == 0
            || !controls.relative_tolerance.is_finite()
            || controls.relative_tolerance <= 0.0
            || !controls.relaxation.is_finite()
            || controls.relaxation <= 0.0
            || controls.relaxation > 1.0
        {
            return Err(StokesBoundaryError::InvalidInput(
                "finite velocity dofs, positive tolerance/budget and relaxation in (0,1] required",
            ));
        }
        let admitted = self.admit_stokes_problem(params, boundaries)?;
        let assemble = |velocity: &[f64]| {
            let rates = self.shear_rates(velocity);
            if rates.iter().any(|rate| !rate.is_finite()) {
                return Err(StokesBoundaryError::InvalidInput("nonfinite shear rate"));
            }
            self.assemble_stokes_admitted_viscosity(
                params,
                &|cell, point| law(cell, point, rates[cell]),
                force,
                &admitted,
            )
        };
        let mut current = vec![0.0; self.n];
        current[..self.n_u].copy_from_slice(initial_velocity);
        for iteration in 1..=controls.maximum_iterations {
            let assembled = assemble(&current);
            if !source_is_valid() {
                return Err(StokesBoundaryError::Rheology(FluxRheologyError::Invalid(
                    "constitutive source refused",
                )));
            }
            let linear = assembled?.solve()?;
            let mut next = linear.field;
            let mut difference = 0.0_f64;
            let mut magnitude = 0.0_f64;
            let mut old_magnitude = 0.0_f64;
            for i in 0..self.n_u {
                difference = difference.hypot(next.x[i] - current[i]);
                magnitude = magnitude.hypot(next.x[i]);
                old_magnitude = old_magnitude.hypot(current[i]);
            }
            let update = difference / magnitude.max(old_magnitude).max(f64::MIN_POSITIVE);
            if !update.is_finite() {
                return Err(StokesBoundaryError::InvalidInput(
                    "unrepresentable constitutive correction",
                ));
            }
            if iteration > 1 && controls.relaxation < 1.0 {
                for (value, old) in next.x.iter_mut().zip(&current) {
                    *value = (1.0 - controls.relaxation) * old + controls.relaxation * *value;
                }
            }
            if next.x.iter().any(|value| !value.is_finite()) {
                return Err(StokesBoundaryError::InvalidInput("nonfinite relaxed flow"));
            }
            if update <= controls.relative_tolerance || iteration == controls.maximum_iterations {
                let endpoint = assemble(&next.x);
                if !source_is_valid() {
                    return Err(StokesBoundaryError::Rheology(FluxRheologyError::Invalid(
                        "constitutive source refused",
                    )));
                }
                let endpoint = endpoint?;
                let mut residual = vec![0.0; self.n];
                let (original, scaled) = crate::ns::sparse::residuals(
                    &endpoint.matrix,
                    &endpoint.rhs,
                    &next.x,
                    &mut residual,
                );
                let defect = if original.is_finite() && scaled.is_finite() {
                    original.max(scaled)
                } else {
                    f64::INFINITY
                };
                if update <= controls.relative_tolerance && defect <= controls.relative_tolerance {
                    next.rel_residual = original;
                    next.scaled_rel_residual = scaled;
                    return Ok(StokesRheologySolution {
                        solution: endpoint.report_solution(next)?,
                        iterations: iteration,
                        relative_update: update,
                        nonlinear_residual: defect,
                    });
                }
                if iteration == controls.maximum_iterations {
                    return Err(StokesBoundaryError::Rheology(
                        FluxRheologyError::NonlinearConvergence {
                            iterations: iteration,
                            relative_update: update,
                            residual: defect,
                        },
                    ));
                }
            }
            current = next.x;
        }
        unreachable!("a nonzero budget accepts or refuses its final candidate")
    }
}
