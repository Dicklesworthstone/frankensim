//! Temperature-dependent source viscosity on the existing fixed Stokes owners.
//!
//! Each interval samples the original finite thermal cells BEFORE its flow
//! solve, resolves one source material at their volume-weighted temperature,
//! and holds that material during the existing implicit motion/heat step.
//! This is a first-order material split, not a converged thermoviscous solve.
//! A caller-authored temperature-spread bound limits the homogeneous-material
//! approximation. Density must remain exactly the originally admitted density;
//! neither thermal expansion nor variable-mass momentum is introduced.
//!
//! The original caloric laws, contacts and thermal observers remain their own
//! owners. New source heat capacity/conductivity do not replace those laws or
//! conductances. The caller must resolve the same physical material and refine
//! both motion and thermal timesteps to assess the split. Geometry remains at
//! the reference pose; advection and pressure work are absent.
//!
//! Ordinary motion/heat methods continue using the construction material.
//! Call the feedback method on EVERY interval that should sample temperature.
//! Each receipt retains the actual resolved source, including its query and
//! evidence, rather than replacing it with an anonymous transport tuple.

use fs_flux::MaterialFluxError;
use fs_material::MaterialError;
use fs_material::fluid_transport::{
    FluidTransport, FluidTransportState, ThermalFluidTransport, ThermalFluidTransportState,
};
use fs_thermal::enthalpy::{DepositionObserver, DepositionSimulation, StepConfig};

use super::heating::{
    GuidedStokesHeatError, GuidedStokesHeatStep, GuidedStokesObservedHeatError, StokesHeatBinding,
    StokesHeatFailure, StokesSlipHeatRoute, StokesThermalCells, unobserved_error,
};
use super::{
    GuidedStokesControls, GuidedStokesMotion, MaterialStokesBoundary, nonnegative, poll, positive,
};

/// Original fluid-cell temperatures at the shared START of a motion interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StokesTemperatureSample {
    /// Volume-weighted temperature [K] supplied to the source resolver.
    pub temperature_k: f64,
    /// Smallest mapped fluid-cell temperature [K].
    pub minimum_temperature_k: f64,
    /// Largest mapped fluid-cell temperature [K].
    pub maximum_temperature_k: f64,
}

/// Accepted original motion/thermal receipts and the material that drove them.
pub struct StokesTemperatureFeedbackStep<M, H> {
    /// The original resolved material, retaining its source and query authority.
    pub material: M,
    /// Checked thermal view of that material; exactly the flow transport tuple.
    pub transport: ThermalFluidTransportState,
    /// Thermal state at which the source was evaluated before this interval.
    pub sample: StokesTemperatureSample,
    /// Existing guided or planar heat receipt; its heat is already in its audit.
    pub heat: H,
}

impl<F: FluidTransport + ?Sized> GuidedStokesMotion<'_, '_, '_, F> {
    /// Refresh source viscosity from the original thermal state, then advance
    /// the existing guided body/flow and finite heat owners once.
    ///
    /// `resolve` is a pure material query at the supplied absolute temperature.
    /// Its returned source must have coherent thermal/momentum views, that
    /// same temperature (up to roundoff), and exactly the original density.
    /// The mapped cells' temperature range must not exceed the finite,
    /// nonnegative `maximum_temperature_spread_k`. This bounds a homogeneous
    /// viscosity approximation, not the error of a nonlinear within-step solve.
    ///
    /// `routes` has the complete ownership requirements of
    /// [`Self::step_with_slip_heating`]; use an empty slice without Navier slip.
    /// Material sampling precedes flow, and the accepted endpoint's physical
    /// bulk/slip heat reaches the existing thermal timeline. The original
    /// frozen binding is retained for ordinary steps; call this method again
    /// to refresh the material on the next interval.
    ///
    /// # Errors
    /// Original map, source, density, temperature, solver, heat or cancellation
    /// refusals preserve both accepted owners, including all clocks.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_temperature_feedback<M: ThermalFluidTransport>(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        cells: StokesThermalCells<'_>,
        routes: &[StokesSlipHeatRoute],
        maximum_temperature_spread_k: f64,
        resolve: impl FnOnce(f64) -> Result<M, MaterialError>,
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: f64,
        controls: GuidedStokesControls,
        thermal_config: StepConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StokesTemperatureFeedbackStep<M, GuidedStokesHeatStep>, GuidedStokesHeatError> {
        self.step_with_temperature_feedback_observed(
            thermal, &mut (), cells, routes, maximum_temperature_spread_k, resolve,
            timestep_s, force_density, exterior, pressure_offsets_pa, external_force_n,
            controls, thermal_config, cancelled,
        )
        .map_err(unobserved_error)
    }

    /// Also retain another physical owner on the common finite thermal clock.
    ///
    /// The original observed heat path forwards every observer hook and keeps
    /// physical viscous/slip heat separate from the additional source. Observer
    /// cloning must isolate tentative history as required by
    /// `DepositionObserver`. The material resolver and source queries must be
    /// pure; externally shared side effects cannot be rolled back.
    ///
    /// # Errors
    /// Every refusal preserves motion, thermal state and observer history.
    /// The original observer's refusal type is retained.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_temperature_feedback_observed<M: ThermalFluidTransport, O: DepositionObserver>(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        observer: &mut O,
        cells: StokesThermalCells<'_>,
        routes: &[StokesSlipHeatRoute],
        maximum_temperature_spread_k: f64,
        resolve: impl FnOnce(f64) -> Result<M, MaterialError>,
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: f64,
        controls: GuidedStokesControls,
        thermal_config: StepConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StokesTemperatureFeedbackStep<M, GuidedStokesHeatStep>,
        GuidedStokesObservedHeatError<O::Error>> {
        poll(&mut cancelled)?;
        let binding = StokesHeatBinding {
            mesh: self.system.mesh,
            density_kg_m3: self.fluid.transport().density_kg_m3(),
            depth_m: self.config.depth_m,
            time_s: self.state.time_s,
        };
        binding.check_thermal_cells(thermal, cells, &mut cancelled)?;
        binding.check_slip_heat_routes(
            thermal, cells, &self.edges, exterior, routes, &mut cancelled,
        )?;
        let (material, transport, sample) = resolve_material(
            &binding, thermal, cells, maximum_temperature_spread_k, resolve, &mut cancelled,
        )?;
        let (state, heat) = {
            let fluid = self.fluid.with_material(&material).map_err(StokesHeatFailure::from)?;
            check_transport(&transport, fluid.transport(), binding.density_kg_m3)?;
            // Continue the SAME accepted native history and reference pose.
            // Only the prepared material binding changes in this candidate.
            let mut candidate = GuidedStokesMotion {
                system: self.system,
                fluid,
                properties: self.properties,
                reference_pose: self.reference_pose,
                config: self.config,
                edges: self.edges.clone(),
                state: self.state.clone(),
            };
            let heat = candidate.step_with_slip_heating_observed(
                thermal, observer, cells, routes, timestep_s,
                force_density, exterior, pressure_offsets_pa, external_force_n,
                controls, thermal_config, &mut cancelled,
            )?;
            (candidate.state, heat)
        };
        // The existing observed step has completed its final cancellation
        // gate and committed the thermal/observer candidates. Publish the same
        // native motion history immediately; nothing fallible follows it.
        self.state = state;
        Ok(StokesTemperatureFeedbackStep { material, transport, sample, heat })
    }
}

/// Shared guided/planar material admission. The caller first admits the
/// complete cell map and slip routes against the unchanged physical owners.
pub(super) fn resolve_material<M: ThermalFluidTransport>(
    binding: &StokesHeatBinding<'_>,
    thermal: &DepositionSimulation<'_>,
    cells: StokesThermalCells<'_>,
    maximum_temperature_spread_k: f64,
    resolve: impl FnOnce(f64) -> Result<M, MaterialError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(M, ThermalFluidTransportState, StokesTemperatureSample), StokesHeatFailure> {
    if !nonnegative(maximum_temperature_spread_k) {
        return Err(StokesHeatFailure::Binding(
            "finite nonnegative homogeneous temperature-spread bound required",
        ));
    }
    let total_area: f64 = binding.mesh.areas.iter().sum();
    if !positive(total_area) || cells.triangle_cells.is_empty() {
        return Err(StokesHeatFailure::Binding("representable nonempty fluid volume required"));
    }
    let reference = thermal.temperature_k(cells.triangle_cells[0])?;
    let mut minimum = reference;
    let mut maximum = reference;
    let mut delta = 0.0;
    for (index, (&cell, &area)) in cells.triangle_cells.iter().zip(&binding.mesh.areas).enumerate() {
        if index % 256 == 0 && cancelled() {
            return Err(StokesHeatFailure::Cancelled);
        }
        let temperature = thermal.temperature_k(cell)?;
        minimum = minimum.min(temperature);
        maximum = maximum.max(temperature);
        // Normalize first, avoiding overflow of temperature times volume.
        delta += (area / total_area) * (temperature - reference);
    }
    let temperature = reference + delta;
    if !positive(temperature) || maximum - minimum > maximum_temperature_spread_k {
        return Err(StokesHeatFailure::Binding(
            "fluid-cell temperatures exceed the homogeneous material approximation",
        ));
    }
    if cancelled() {
        return Err(StokesHeatFailure::Cancelled);
    }
    let material = resolve(temperature)
        .map_err(MaterialFluxError::from)
        .map_err(StokesHeatFailure::from)?;
    let transport = material.thermal_transport_state()
        .map_err(MaterialFluxError::from)
        .map_err(StokesHeatFailure::from)?;
    let evaluated = transport.temperature_k();
    if (evaluated - temperature).abs()
        > 64.0 * f64::EPSILON * evaluated.abs().max(temperature.abs())
    {
        return Err(StokesHeatFailure::Binding(
            "resolved fluid source must describe the sampled absolute temperature",
        ));
    }
    Ok((material, transport, StokesTemperatureSample {
        temperature_k: temperature,
        minimum_temperature_k: minimum,
        maximum_temperature_k: maximum,
    }))
}

pub(super) fn check_transport(
    thermal: &ThermalFluidTransportState,
    momentum: &FluidTransportState,
    reference_density_kg_m3: f64,
) -> Result<(), StokesHeatFailure> {
    if thermal.transport() != *momentum {
        return Err(StokesHeatFailure::Binding(
            "source thermal and momentum transport views must agree exactly",
        ));
    }
    if momentum.density_kg_m3() != reference_density_kg_m3 {
        return Err(StokesHeatFailure::Binding(
            "fixed-mass Stokes feedback requires exactly unchanged fluid density",
        ));
    }
    Ok(())
}
