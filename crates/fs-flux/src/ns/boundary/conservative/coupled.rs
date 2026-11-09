//! Endpoint-advector coupling of conservative density and native BDM momentum.
//!
//! Each nonlinear trial starts from the SAME immutable old density, inventories
//! and velocity. Only the volume advector changes. The accepted transport uses
//! the returned endpoint velocity, and its original compatible mass flux enters
//! the reassembled momentum equation and physical observations.

use std::cell::RefCell;

use super::{
    DensityMomentumControls, DensityMomentumStep, DensityTransportCell, FluxSystem,
    MaterialFluxNumerics, MaterialNavierStokesBoundary, Result, invalid, observe_step, poll,
};
use crate::ns::boundary::flow_control::FlowControl;
use crate::ns::boundary::{StokesBoundaryError, StokesFlowControlSolution, StokesRheologyControls};

/// Finite external RCR storage coupled to the same transported-density endpoint.
#[path = "coupled/windkessel.rs"]
pub mod windkessel;

/// Joint density/inventory/flow endpoint at its own transport velocity.
pub struct CoupledDensityMomentumStep<const N: usize> {
    /// Physical inventories, pressure, boundary reactions and native velocity.
    /// Both reported flow residuals belong to the reassembled coupled endpoint.
    pub step: DensityMomentumStep<N>,
    /// Completed native Picard momentum solves.
    pub iterations: usize,
    /// Last velocity correction BEFORE under-relaxation.
    pub relative_update: f64,
    /// Maximum original and row-scaled coupled endpoint momentum residual.
    pub nonlinear_residual: f64,
}

/// Coupled transported-density endpoint with native integral volume-flow ports.
pub struct CoupledDensityMomentumFlowControlStep<const N: usize> {
    /// Accepted transport, physical momentum and nonlinear endpoint evidence.
    pub step: CoupledDensityMomentumStep<N>,
    /// Uniform pressure added to each port's original traction as `-p*n` [Pa].
    /// The same increments enter the returned pressure and boundary reactions.
    pub pressure_offsets_pa: Vec<f64>,
    /// Actual signed outward volume rate per unit depth [m²/s] at each port.
    /// Multiply by the declared depth for m³/s; this is not a mass-flow rate.
    pub outward_volume_rates_m2_s: Vec<f64>,
    /// Largest absolute endpoint error in the requested volume rates [m²/s].
    pub maximum_flow_residual_m2_s: f64,
}

impl<'m> FluxSystem<'m> {
    /// Couple the endpoint volume velocity to density and native BDF1 momentum.
    ///
    /// The original conservative owner advances old density and every passive
    /// inventory with the current Picard velocity, then assembles momentum with
    /// that SAME transport's final density, old inertia and compatible mass flux.
    /// Force and dynamic viscosity are evaluated at these candidate cells.
    /// Every trial starts again from the immutable old state; rejected nonlinear
    /// iterations are not additional timesteps and never consume inventories.
    ///
    /// The native Picard owner checks the UNRELAXED velocity correction and both
    /// original/row-scaled momentum residuals reassembled at the returned field.
    /// The accepted density, mass flux, force, viscosity, pressure and reactions
    /// all come from that last endpoint assembly, not the preceding trial.
    /// All input callbacks must be pure throughout the solve.
    ///
    /// This is first-order endpoint-ADVECTOR coupling: donor densities/specifics
    /// still originate in the old state and use this FluxSystem's configured
    /// spatial reconstruction.
    /// It is not a fully implicit density update. Positivity retains the original
    /// outgoing-volume CFL cap, enforced at EVERY trial, and the stated local
    /// divergence tolerance. The caller reduces the timestep after a refusal;
    /// no advector projection, flux clipping or unreported substep occurs.
    ///
    /// The previous velocity is the initial guess. Its divergence must already
    /// meet transport admission, while endpoint boundary values may change.
    /// Picard relaxation preserves incompressibility between admitted iterates.
    /// Boundary donor states must cover every incoming stream encountered,
    /// including any flow reversal. Missing donors refuse the complete step.
    ///
    /// No phase-volume source, compressible EOS, heat diffusion, density-source
    /// inference, new rheology law or full coupled adjoint is supplied. Numerical
    /// kinetic losses remain separate from physical heating.
    ///
    /// Storage and synchronous linear-solve limitations follow the original
    /// owners. Each trial discards the preceding transport before preparing its
    /// replacement; no iteration history is accumulated. Cancellation is polled
    /// around native solves and at transport/observation cell and edge boundaries.
    /// A refused/cancelled solve publishes no partial state.
    #[allow(clippy::too_many_arguments)]
    pub fn density_momentum_coupled_bdf1_step<const N: usize>(
        &self,
        numerics: MaterialFluxNumerics,
        dynamic_viscosity: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> f64,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensityMomentumControls,
        nonlinear: StokesRheologyControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CoupledDensityMomentumStep<N>> {
        self.coupled_density_momentum_bdf1(
            numerics,
            dynamic_viscosity,
            force_density,
            boundaries,
            previous_velocity,
            cells,
            inflow,
            timestep_s,
            controls,
            nonlinear,
            None,
            cancelled,
        )
        .map(|solved| solved.flow)
    }

    /// Couple transported density and native BDF1 momentum to integral port flow.
    ///
    /// The port pressures are unknowns of the SAME nonlinear endpoint equations.
    /// Each frozen trial resolves them through the native pressure responses;
    /// their endpoint traction, momentum residual and boundary work are retained.
    /// No velocity profile is imposed at the controlled ports. Only the signed
    /// outward integral is constrained. The original traction remains present,
    /// including its tangential part; the unknown adds uniform `-p*n`.
    ///
    /// Each nonempty port must contain disjoint original natural-traction edges
    /// in one fluid component. Every controlled component needs an uncontrolled
    /// natural-traction edge to retain its pressure datum. Both `Open` and
    /// `OpenReservoir` retain their authored momentum/backflow behavior.
    /// The caller must supply density/inventory donors for every inflow reached
    /// by any trial, including reversal at an uncontrolled outlet.
    ///
    /// Targets and their absolute tolerance are volume rates per unit depth
    /// [m²/s]. Divide an authored m³/s target by `controls.depth_m` first.
    /// Mass flow follows the original transported donors and cannot be inferred
    /// by multiplying these rates by the conditioning reference density.
    /// Returned pressure increments are physical Pa; the reference density only
    /// normalizes the common equations. Physical boundary work uses that same
    /// normalization and the declared depth, exactly as for the retained flow.
    ///
    /// All endpoint-advector, reconstruction, CFL, cancellation and atomic
    /// refusal semantics of [`Self::density_momentum_coupled_bdf1_step`] apply.
    /// Pressure responses are rebuilt when candidate density or momentum changes;
    /// no response from an earlier frozen trial is reused.
    #[allow(clippy::too_many_arguments)]
    pub fn density_momentum_coupled_bdf1_flow_control_step<const N: usize>(
        &self,
        numerics: MaterialFluxNumerics,
        dynamic_viscosity: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> f64,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        ports: &[&[usize]],
        outward_volume_rates_m2_s: &[f64],
        absolute_tolerance_m2_s: f64,
        previous_velocity: &[f64],
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensityMomentumControls,
        nonlinear: StokesRheologyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoupledDensityMomentumFlowControlStep<N>> {
        let solved = self.coupled_density_momentum_bdf1(
            numerics,
            dynamic_viscosity,
            force_density,
            boundaries,
            previous_velocity,
            cells,
            inflow,
            timestep_s,
            controls,
            nonlinear,
            Some(FlowControl::targets(
                ports,
                outward_volume_rates_m2_s,
                absolute_tolerance_m2_s,
            )),
            &mut cancelled,
        )?;
        let mut pressure_offsets_pa = solved.pressure_offsets_m2_s2;
        for pressure in &mut pressure_offsets_pa {
            poll(&mut cancelled)?;
            let physical = *pressure * controls.reference_density_kg_m3;
            if !physical.is_finite() || (physical == 0.0 && *pressure != 0.0) {
                return Err(invalid("unrepresentable physical port pressure"));
            }
            *pressure = physical;
        }
        poll(&mut cancelled)?;
        Ok(CoupledDensityMomentumFlowControlStep {
            step: solved.flow,
            pressure_offsets_pa,
            outward_volume_rates_m2_s: solved.outward_volume_rates_m2_s,
            maximum_flow_residual_m2_s: solved.maximum_flow_residual_m2_s,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn coupled_density_momentum_bdf1<const N: usize>(
        &self,
        numerics: MaterialFluxNumerics,
        dynamic_viscosity: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> f64,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensityMomentumControls,
        nonlinear: StokesRheologyControls,
        flow_control: Option<FlowControl<'_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StokesFlowControlSolution<CoupledDensityMomentumStep<N>>> {
        poll(&mut cancelled)?;
        self.admit_mixed_bdf1_history(previous_velocity, timestep_s)?;
        if nonlinear.maximum_iterations == 0
            || !nonlinear.relative_tolerance.is_finite()
            || nonlinear.relative_tolerance <= 0.0
            || !nonlinear.relaxation.is_finite()
            || nonlinear.relaxation <= 0.0
            || nonlinear.relaxation > 1.0
        {
            return Err(invalid(
                "positive nonlinear budget/tolerance and relaxation in (0,1] required",
            ));
        }
        let cancelled = RefCell::new(cancelled);
        let failure = RefCell::new(None);
        let retained = RefCell::new(None);
        let assemble = |velocity: &[f64]| {
            let candidate = (|| -> Result<_> {
                poll(&mut *cancelled.borrow_mut())?;
                // An old candidate owns no accepted state. Drop it before the
                // replacement's transport allocation, even when this trial fails.
                drop(retained.borrow_mut().take());
                let mut prepared = self.prepare_density_momentum_bdf1(
                    numerics,
                    force_density,
                    boundaries,
                    &velocity[..self.n_u],
                    previous_velocity,
                    cells,
                    inflow,
                    timestep_s,
                    controls,
                    &mut *cancelled.borrow_mut(),
                )?;
                let weights = core::mem::take(&mut prepared.current_weights);
                let momentum = prepared.assemble(
                    &|cell, point| dynamic_viscosity(cell, point, &prepared.cells[cell]),
                    weights,
                )?;
                poll(&mut *cancelled.borrow_mut())?;
                *retained.borrow_mut() = Some(prepared);
                Ok(momentum)
            })();
            candidate.map_err(|error| {
                *failure.borrow_mut() = Some(error);
                StokesBoundaryError::InvalidInput("coupled density/momentum assembly refused")
            })
        };
        let solved = self.solve_mixed_picard_with_flow_control(
            previous_velocity,
            nonlinear,
            0.0,
            &assemble,
            &|| failure.borrow().is_none(),
            flow_control,
        );
        if let Some(error) = failure.into_inner() {
            return Err(error);
        }
        let solved = solved?;
        let prepared = retained
            .into_inner()
            .ok_or_else(|| invalid("accepted coupled endpoint transport required"))?;
        let step = observe_step(
            self,
            solved.flow.solution,
            &prepared.cells,
            prepared.transport,
            cells,
            previous_velocity,
            &prepared.conservative.faces,
            controls,
            &mut *cancelled.borrow_mut(),
        )?;
        poll(&mut *cancelled.borrow_mut())?;
        Ok(StokesFlowControlSolution {
            flow: CoupledDensityMomentumStep {
                step,
                iterations: solved.flow.iterations,
                relative_update: solved.flow.relative_update,
                nonlinear_residual: solved.flow.nonlinear_residual,
            },
            pressure_offsets_m2_s2: solved.pressure_offsets_m2_s2,
            outward_volume_rates_m2_s: solved.outward_volume_rates_m2_s,
            maximum_flow_residual_m2_s: solved.maximum_flow_residual_m2_s,
        })
    }
}
