//! Checked generalized-Newtonian Stokes flow with open traction boundaries.

use super::flow_control::{FlowControl, checked_port_flows, pressure_rhs};
use super::{
    FluxParams, FluxSystem, StokesAssembly, StokesBoundary, StokesBoundaryError,
    StokesFlowControlSolution, StokesSolution,
};
use crate::ns::FluxRheologyError;

/// Controls for frozen-load mixed constitutive or convective Picard solves.
#[derive(Clone, Copy, Debug)]
pub struct StokesRheologyControls {
    /// Maximum number of Picard flow updates. With integral flow control, each
    /// update solves one baseline and one pressure response per controlled port.
    pub maximum_iterations: usize,
    /// Required relative velocity correction AND both nonlinear residuals.
    pub relative_tolerance: f64,
    /// Fraction of each correction used after the first admissible flow.
    /// Must lie in (0,1]. Under-relaxation can stabilize Picard corrections;
    /// it is not a convergence guarantee for arbitrary laws or flow data.
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
    /// Completed Picard flow updates, bounded by maximum_iterations.
    /// Flow control uses one baseline plus one pressure solve per port per update.
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
        self.solve_stokes_mixed_rheology_with_control(
            params,
            law,
            force,
            boundaries,
            initial_velocity,
            controls,
            None,
            source_is_valid,
        )
        .map(|result| result.flow)
    }

    /// Solve shear-dependent creeping flow with prescribed integral port flows.
    ///
    /// Each Picard update resolves the pressure increments at the supplied
    /// traction ports using the current viscosity. Neither port velocity
    /// profiles nor a frozen viscosity are imposed. Pressure offsets relax
    /// with the velocity and pressure field. Final acceptance re-evaluates
    /// viscosity at the returned shear and checks both nonlinear equation
    /// residuals plus every requested flow against its original tolerance.
    ///
    /// Law, force and traction have the same units as
    /// [`Self::solve_stokes_mixed_rheology`]. Targets are signed outward m²/s
    /// per depth; pressure increments are divided by the same reference
    /// density as the equations. Ports obey [`StokesAssembly::prepare_flow_control`]
    /// ownership and pressure-datum requirements. Boundary samples are fixed
    /// throughout the solve. Each update costs one baseline solve plus one
    /// pressure-response solve per port, within the declared Picard budget.
    ///
    /// # Errors
    /// Preserves native material, boundary, port and linear/nonlinear
    /// convergence checks. Refusal never modifies the initial velocity or
    /// publishes a partially converged field or pressure.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_stokes_mixed_rheology_flow_control(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        ports: &[&[usize]],
        outward_volume_rates_m2_s: &[f64],
        initial_velocity: &[f64],
        controls: StokesRheologyControls,
        absolute_tolerance_m2_s: f64,
    ) -> Result<StokesFlowControlSolution<StokesRheologySolution>, StokesBoundaryError> {
        self.solve_stokes_mixed_rheology_flow_control_checked(
            params,
            law,
            force,
            boundaries,
            ports,
            outward_volume_rates_m2_s,
            initial_velocity,
            controls,
            absolute_tolerance_m2_s,
            &|| true,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve_stokes_mixed_rheology_flow_control_checked(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        ports: &[&[usize]],
        outward_volume_rates_m2_s: &[f64],
        initial_velocity: &[f64],
        controls: StokesRheologyControls,
        absolute_tolerance_m2_s: f64,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<StokesFlowControlSolution<StokesRheologySolution>, StokesBoundaryError> {
        self.solve_stokes_mixed_rheology_with_control(
            params,
            law,
            force,
            boundaries,
            initial_velocity,
            controls,
            Some(FlowControl {
                ports,
                targets: outward_volume_rates_m2_s,
                absolute_tolerance_m2_s,
            }),
            source_is_valid,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_stokes_mixed_rheology_with_control(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: StokesRheologyControls,
        flow_control: Option<FlowControl<'_>>,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<StokesFlowControlSolution<StokesRheologySolution>, StokesBoundaryError> {
        if let Some(control) = flow_control {
            control.admit()?;
        }
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
        self.solve_mixed_picard_with_flow_control(
            initial_velocity,
            controls,
            0.0,
            &assemble,
            source_is_valid,
            flow_control,
        )
    }

    // Shared bounded iteration and endpoint acceptance for constitutive Stokes
    // and convective mixed-boundary flow. Admission precedes this private loop.
    pub(super) fn solve_mixed_picard<'a>(
        &'a self,
        initial_velocity: &[f64],
        controls: StokesRheologyControls,
        assemble: &dyn Fn(&[f64]) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError>,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<StokesRheologySolution, StokesBoundaryError> {
        self.solve_mixed_picard_with_velocity_tolerance(
            initial_velocity,
            controls,
            0.0,
            assemble,
            source_is_valid,
        )
    }

    // The absolute criterion measures the MAXIMUM unrelaxed correction of
    // the original affine BDM normal velocity traces [m/s]. Their endpoint
    // values are m0 +/- 6*m1, so the exact supremum is |dm0| + 6*|dm1|.
    // This does not grow with edge count or mix pressure into a velocity norm.
    // Zero retains the original relative-only convergence contract.
    pub(super) fn solve_mixed_picard_with_velocity_tolerance<'a>(
        &'a self,
        initial_velocity: &[f64],
        controls: StokesRheologyControls,
        absolute_velocity_tolerance_m_s: f64,
        assemble: &dyn Fn(&[f64]) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError>,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<StokesRheologySolution, StokesBoundaryError> {
        self.solve_mixed_picard_with_flow_control(
            initial_velocity,
            controls,
            absolute_velocity_tolerance_m_s,
            assemble,
            source_is_valid,
            None,
        )
        .map(|result| result.flow)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn solve_mixed_picard_with_flow_control<'a>(
        &'a self,
        initial_velocity: &[f64],
        controls: StokesRheologyControls,
        absolute_velocity_tolerance_m_s: f64,
        assemble: &dyn Fn(&[f64]) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError>,
        source_is_valid: &dyn Fn() -> bool,
        flow_control: Option<FlowControl<'_>>,
    ) -> Result<StokesFlowControlSolution<StokesRheologySolution>, StokesBoundaryError> {
        if !absolute_velocity_tolerance_m_s.is_finite() || absolute_velocity_tolerance_m_s < 0.0 {
            return Err(StokesBoundaryError::InvalidInput(
                "finite nonnegative absolute velocity correction tolerance required",
            ));
        }
        let mut current = vec![0.0; self.n];
        current[..self.n_u].copy_from_slice(initial_velocity);
        let mut current_pressures = Vec::new();
        for iteration in 1..=controls.maximum_iterations {
            let assembled = assemble(&current);
            if !source_is_valid() {
                return Err(StokesBoundaryError::Rheology(FluxRheologyError::Invalid(
                    "constitutive source refused",
                )));
            }
            let (linear, mut pressures) = if let Some(control) = flow_control {
                let controlled = assembled?
                    .prepare_flow_control(control.ports)?
                    .solve(control.targets, control.absolute_tolerance_m2_s)?;
                (controlled.flow, controlled.pressure_offsets_m2_s2)
            } else {
                (assembled?.solve()?, Vec::new())
            };
            let mut next = linear.field;
            let mut difference = 0.0_f64;
            let mut absolute_update_m_s = 0.0_f64;
            let mut magnitude = 0.0_f64;
            let mut old_magnitude = 0.0_f64;
            for i in 0..self.n_u {
                difference = difference.hypot(next.x[i] - current[i]);
                if i % 2 == 0 {
                    absolute_update_m_s = absolute_update_m_s.max(
                        (next.x[i] - current[i]).abs()
                            + 6.0 * (next.x[i + 1] - current[i + 1]).abs(),
                    );
                }
                magnitude = magnitude.hypot(next.x[i]);
                old_magnitude = old_magnitude.hypot(current[i]);
            }
            let update = difference / magnitude.max(old_magnitude).max(f64::MIN_POSITIVE);
            if !update.is_finite() {
                return Err(StokesBoundaryError::InvalidInput(
                    "unrepresentable constitutive correction",
                ));
            }
            let correction_converged = update <= controls.relative_tolerance
                || absolute_update_m_s <= absolute_velocity_tolerance_m_s;
            if iteration > 1 && controls.relaxation < 1.0 {
                for (value, old) in next.x.iter_mut().zip(&current) {
                    *value = (1.0 - controls.relaxation) * old + controls.relaxation * *value;
                }
                // These are additional pressure unknowns of the same coupled
                // equations. Relaxing only the field would mismatch its load.
                for (value, old) in pressures.iter_mut().zip(&current_pressures) {
                    *value = (1.0 - controls.relaxation) * old + controls.relaxation * *value;
                }
            }
            if next.x.iter().any(|value| !value.is_finite()) {
                return Err(StokesBoundaryError::InvalidInput("nonfinite relaxed flow"));
            }
            if pressures.iter().any(|value| !value.is_finite()) {
                return Err(StokesBoundaryError::InvalidInput(
                    "nonfinite relaxed port pressure",
                ));
            }
            if correction_converged || iteration == controls.maximum_iterations {
                let endpoint = assemble(&next.x);
                if !source_is_valid() {
                    return Err(StokesBoundaryError::Rheology(FluxRheologyError::Invalid(
                        "constitutive source refused",
                    )));
                }
                let mut endpoint = endpoint?;
                let (flows, maximum_flow_residual_m2_s) = if let Some(control) = flow_control {
                    endpoint.rhs = pressure_rhs(&endpoint, control.ports, &pressures)?;
                    checked_port_flows(
                        &endpoint,
                        control.ports,
                        control.targets,
                        control.absolute_tolerance_m2_s,
                        &next.x,
                    )?
                } else {
                    (Vec::new(), 0.0)
                };
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
                if correction_converged && defect <= controls.relative_tolerance {
                    next.rel_residual = original;
                    next.scaled_rel_residual = scaled;
                    let mut solution = endpoint.report_solution(next)?;
                    if let Some(control) = flow_control {
                        solution
                            .reactions
                            .add_port_pressures(control.ports, &pressures)?;
                    }
                    return Ok(StokesFlowControlSolution {
                        flow: StokesRheologySolution {
                            solution,
                            iterations: iteration,
                            relative_update: update,
                            nonlinear_residual: defect,
                        },
                        pressure_offsets_m2_s2: pressures,
                        outward_volume_rates_m2_s: flows,
                        maximum_flow_residual_m2_s,
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
            current_pressures = pressures;
        }
        unreachable!("a nonzero budget accepts or refuses its final candidate")
    }
}
