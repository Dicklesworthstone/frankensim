//! Temperature-refreshed material on the original planar motion and heat path.

use fs_la::LuWorkspace;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, ThermalFluidTransport};
use fs_thermal::enthalpy::{DepositionObserver, DepositionSimulation, StepConfig};

use super::{PlanarStokesControls, PlanarStokesError, PlanarStokesMotion, poll};
use super::super::feedback::{
    StokesTemperatureFeedbackStep, check_transport, resolve_material,
};
use super::super::heating::{
    PlanarStokesHeatError, PlanarStokesHeatStep, PlanarStokesObservedHeatError,
    StokesHeatBinding, StokesHeatFailure, StokesSlipHeatRoute, StokesThermalCells,
    unobserved_planar_error,
};
use super::super::MaterialStokesBoundary;

impl<F: FluidTransport + ?Sized> PlanarStokesMotion<'_, '_, '_, F> {
    /// Resolve viscosity at the original fluid cells' current temperature, then
    /// advance native planar translation, rotation and finite heat once.
    ///
    /// This shares guided temperature/source admission and the original planar
    /// force/torque solve. The finite nonnegative temperature-spread bound limits
    /// the homogeneous material approximation. The returned source must describe
    /// that temperature, agree across thermal and momentum transport views, and
    /// retain exactly the original density and thus the existing cell masses.
    ///
    /// `resolve` is a pure query of the same physical material. Its original
    /// source is retained in the receipt. `routes` must completely partition
    /// exterior Navier-slip losses; supply an empty slice when no slip is present.
    /// Thermal calorics, conductances and contacts remain their original owners.
    ///
    /// Material is sampled at the start and held for this interval. Refine both
    /// motion and thermal steps: this is a first-order split without thermal
    /// advection, a moving mesh, or converged within-step viscosity feedback.
    /// Ordinary methods retain the construction binding; call this method on
    /// each interval that should use the evolved temperature.
    ///
    /// # Errors
    /// Material, motion, thermal or cancellation refusals leave all accepted
    /// physical histories unchanged.
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
        external_force_n: [f64; 2],
        external_torque_n_m: f64,
        controls: PlanarStokesControls,
        thermal_config: StepConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StokesTemperatureFeedbackStep<M, PlanarStokesHeatStep>, PlanarStokesHeatError> {
        self.step_with_temperature_feedback_observed(
            thermal, &mut (), cells, routes, maximum_temperature_spread_k, resolve,
            timestep_s, force_density, exterior, pressure_offsets_pa, external_force_n,
            external_torque_n_m, controls, thermal_config, cancelled,
        )
        .map_err(unobserved_planar_error)
    }

    /// Retain another physical observer on the same proposed thermal timeline.
    ///
    /// All native source, transition, conductance, interval and timestep hooks
    /// remain active. Observer cloning must isolate tentative history, and
    /// material callbacks must be pure. The flow/material candidate and both
    /// thermal/observer candidates remain private through the final cancellation
    /// check; only then are the three accepted histories published together.
    ///
    /// # Errors
    /// The original observer refusal type is preserved. Every late refusal
    /// preserves motion, thermal cells, events, sources and observer history.
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
        external_force_n: [f64; 2],
        external_torque_n_m: f64,
        controls: PlanarStokesControls,
        thermal_config: StepConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StokesTemperatureFeedbackStep<M, PlanarStokesHeatStep>,
        PlanarStokesObservedHeatError<O::Error>> {
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
            let mut candidate = PlanarStokesMotion {
                system: self.system,
                fluid,
                properties: self.properties,
                masses: self.masses,
                reference_pose: self.reference_pose,
                config: self.config,
                edges: self.edges.clone(),
                schur: LuWorkspace::new(3).map_err(PlanarStokesError::from)?,
                state: self.state.clone(),
            };
            let heat = candidate.step_with_slip_heating_observed(
                thermal, observer, cells, routes, timestep_s,
                force_density, exterior, pressure_offsets_pa, external_force_n,
                external_torque_n_m, controls, thermal_config, &mut cancelled,
            )?;
            (candidate.state, heat)
        };
        // Native heat acceptance is the final fallible boundary. Complete the
        // motion publication immediately; no extra field copies or late poll
        // may split this same accepted physical instant.
        self.state = state;
        Ok(StokesTemperatureFeedbackStep { material, transport, sample, heat })
    }
}
