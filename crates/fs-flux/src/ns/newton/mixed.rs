//! Mixed-boundary Navier--Stokes on the original consistent Newton owner.

use super::{FluxNewtonControls, FluxNewtonError};
use crate::ns::boundary::{NavierStokesBoundary, OpenBoundaryPolicy, StokesSolution};
use crate::ns::linearization::{
    FluxRheologyLinearization, FluxViscosityTangent, rheology_residual,
};
use crate::ns::{FluxParams, FluxSystem};

/// Accepted mixed flow, retaining its pressure datum and physical reactions.
pub struct FluxMixedNewtonStep {
    /// Original field, component pressure datum and endpoint reactions.
    pub solution: StokesSolution,
    /// Completed original Newton correction solves.
    pub iterations: usize,
    /// Rejected line-search candidates across all corrections.
    pub backtracks: usize,
    /// Last accepted relative velocity correction.
    pub relative_update: f64,
    /// Maximum original, row-scaled and pressure-origin endpoint defect.
    pub nonlinear_residual: f64,
    /// Actual stress policy on each original open edge, including reservoirs.
    pub open_boundary_policies: Vec<(usize, OpenBoundaryPolicy)>,
}

impl FluxSystem<'_> {
    /// Differentiate the original mixed-boundary Navier--Stokes residual.
    ///
    /// The full state includes velocity and pressure. Zero inverse timestep
    /// selects steady flow; a positive value requires a fixed mass history.
    /// Constitutive, transported-velocity and advecting-velocity derivatives
    /// share the original operator. Natural traction, dissipative backflow,
    /// reservoir inflow and prescribed-normal slip keep their distinct fluxes.
    /// The active nonnegative branch is a semismooth selection at zero flux.
    ///
    /// # Errors
    /// Refuses invalid state/history, unsupported steady boundary partitions
    /// and the original boundary or constitutive admission failures.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_navier_stokes_mixed_rheology(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        state: &[f64],
        inverse_timestep: f64,
        mass_history: Option<&[f64]>,
    ) -> Result<FluxRheologyLinearization, FluxNewtonError> {
        let anchors = vec![inverse_timestep > 0.0; self.pressure_roots.len()];
        let admitted = self
            .admit_navier_stokes_boundaries_with_velocity_anchors(params, boundaries, &anchors)
            .map_err(FluxNewtonError::Boundary)?;
        self.linearize_rheology_admitted(
            params,
            law,
            force,
            &|_| [0.0; 2],
            state,
            true,
            inverse_timestep,
            mass_history,
            Some(&admitted),
        )
        .map_err(FluxNewtonError::Linearization)
    }

    /// Solve steady mixed-boundary flow with the original damped Newton loop.
    ///
    /// No time mass is introduced. Each component needs velocity or Navier-slip
    /// support. Stress, convection and physical boundary reactions use the same
    /// endpoint viscosity and velocity. Natural traction is never silently
    /// changed to a dissipative policy, and a reservoir supplies only incoming
    /// momentum, not a prescribed normal flow rate.
    ///
    /// The law and external data must be pure and fixed throughout the solve.
    /// `params.nu` is unused. The native work limits, pressure-origin residual
    /// checks and constitutive/semismooth qualifications remain authoritative;
    /// this does not establish global flow uniqueness or backflow stability.
    ///
    /// # Errors
    /// Preserves boundary, tangent, linear-solve, line-search and nonlinear
    /// convergence refusals. No candidate or mutated initial state is published.
    pub fn solve_navier_stokes_mixed_rheology_newton(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: FluxNewtonControls,
    ) -> Result<FluxMixedNewtonStep, FluxNewtonError> {
        self.mixed_rheology_newton_endpoint(
            params,
            law,
            force,
            boundaries,
            initial_velocity,
            0.0,
            None,
            controls,
        )
    }

    /// Fully implicit BDF1 with the same admitted mixed-boundary Newton owner.
    ///
    /// Boundary and external fields describe the new time. Finite physical
    /// mass anchors velocity even on an all-open component. Previous velocity
    /// is immutable, and accepted reactions contain the actual endpoint flux.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_navier_stokes_mixed_rheology_newton_step(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: FluxNewtonControls,
    ) -> Result<FluxMixedNewtonStep, FluxNewtonError> {
        self.bdf2_navier_stokes_mixed_rheology_newton_variable_step(
            params,
            law,
            force,
            boundaries,
            previous_velocity,
            None,
            timestep_s,
            controls,
        )
    }

    /// Fully implicit mixed-boundary BDF2 with actual unequal physical intervals.
    ///
    /// `older_history` pairs the older velocity with its interval to the previous
    /// velocity. The original history owner admits new/previous ratios in
    /// `(0, 2]`; `None` is BDF1 startup. Only the mass load is combined:
    /// the actual previous velocity remains the nonlinear initial guess.
    /// Reactions and dissipation are instantaneous endpoint rates, not a
    /// two-history BDF energy balance or integrated thermal deposition.
    ///
    /// # Errors
    /// Refuses invalid histories/intervals and the same physical/numerical
    /// failures as the steady owner. Neither history is modified.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_navier_stokes_mixed_rheology_newton_variable_step(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        older_history: Option<(&[f64], f64)>,
        timestep_s: f64,
        controls: FluxNewtonControls,
    ) -> Result<FluxMixedNewtonStep, FluxNewtonError> {
        if previous_velocity.len() != self.n_u
            || previous_velocity.iter().any(|value| !value.is_finite())
            || !timestep_s.is_finite()
            || timestep_s <= 0.0
            || !timestep_s.recip().is_finite()
            || timestep_s.recip() <= 0.0
        {
            return Err(FluxNewtonError::Invalid(
                "finite original velocity and representable positive timestep required",
            ));
        }
        let (history, effective_timestep) = self
            .mixed_bdf2_history(previous_velocity, older_history, timestep_s)
            .map_err(FluxNewtonError::History)?;
        self.mixed_rheology_newton_endpoint(
            params,
            law,
            force,
            boundaries,
            previous_velocity,
            effective_timestep.recip(),
            Some(history.as_deref().unwrap_or(previous_velocity)),
            controls,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn mixed_rheology_newton_endpoint(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> FluxViscosityTangent,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        initial_velocity: &[f64],
        inverse_timestep: f64,
        mass_history: Option<&[f64]>,
        controls: FluxNewtonControls,
    ) -> Result<FluxMixedNewtonStep, FluxNewtonError> {
        let anchors = vec![inverse_timestep > 0.0; self.pressure_roots.len()];
        let admitted = self
            .admit_navier_stokes_boundaries_with_velocity_anchors(params, boundaries, &anchors)
            .map_err(FluxNewtonError::Boundary)?;
        let mut step = self.rheology_newton_solve_admitted_checked(
            params,
            law,
            force,
            &|_| [0.0; 2],
            initial_velocity,
            true,
            inverse_timestep,
            mass_history,
            controls,
            &|| true,
            Some(&admitted),
        )?;
        let rates = self.shear_rates(&step.solution.x);
        // Reuse the ORIGINAL reaction and pressure-datum owner, including mass
        // history and open advective transport, at this same accepted endpoint.
        let endpoint = self
            .assemble_stokes_admitted_operator(
                params,
                &|cell, point| law(cell, point, rates[cell]).nu,
                force,
                &admitted,
                None,
                inverse_timestep,
                mass_history,
                Some(&step.solution.x[..self.n_u]),
            )
            .map_err(FluxNewtonError::Boundary)?;
        let (_, original, scaled) = rheology_residual(
            self,
            endpoint.matrix(),
            endpoint.rhs(),
            &step.solution.x,
            true,
        );
        step.solution.rel_residual = original;
        step.solution.scaled_rel_residual = scaled;
        if !step.solution.converged(controls.residual_tolerance) {
            return Err(FluxNewtonError::NonConvergence {
                iterations: step.iterations,
                residual: original.max(scaled),
            });
        }
        let policies = (0..self.mesh.edges.len())
            .filter_map(|edge| admitted.open_policy(edge).map(|policy| (edge, policy)))
            .collect();
        Ok(FluxMixedNewtonStep {
            solution: endpoint
                .report_solution(step.solution)
                .map_err(FluxNewtonError::Boundary)?,
            iterations: step.iterations,
            backtracks: step.backtracks,
            relative_update: step.relative_update,
            nonlinear_residual: original.max(scaled),
            open_boundary_policies: policies,
        })
    }
}
