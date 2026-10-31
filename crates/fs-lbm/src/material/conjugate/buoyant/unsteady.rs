//! Source-material natural and mixed convection on the native transient clock.

use super::{MaterialBuoyancyError, MaterialBuoyancySample, poll, samples};
use crate::conjugate::unsteady::{MaterialCoefficients, MaterialHeatStep, compact_owners, march_conjugate_with_material};
use crate::conjugate::{
    ConjugateMarch, EnergyConfig, RadiationConfig, SimpleConfig, SolidMaterial, ThermalSetup,
    TransientRecord, Turbulence, UnsteadyConfig, Voxel, VoxelDomain,
};
use crate::material::conjugate::MaterialChtFluid;
use crate::material::transport::boussinesq::BoussinesqConfig;
use fs_exec::CancelGate;
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
};

mod solid;
pub use solid::{MaterialSolidBuoyancyMarch, SolidThermalConfig, SolidThermalRecord};

/// Native transient solution with the exact final applied and endpoint sources.
pub struct MaterialBuoyancyMarch<F> {
    /// Fixed reference inertia and continuity. Sensible storage uses reference
    /// Cp unless an explicit source-caloric coupling was selected.
    pub reference: FluidTransportState,
    /// Original time histories, conservative fluxes, temperatures and heat ledger.
    pub march: ConjugateMarch,
    /// Source snapshots used for the last momentum and energy step.
    pub applied: Vec<MaterialBuoyancySample<F>>,
    /// Original sources freshly admitted at the returned temperature field.
    pub endpoint: Vec<MaterialBuoyancySample<F>>,
    /// Time of the temperatures used for the last explicit material refresh.
    pub applied_time_s: f64,
    /// Time of the returned thermal state and endpoint source snapshots.
    pub endpoint_time_s: f64,
}

impl<F: FluidTransport + ?Sized> MaterialChtFluid<'_, F> {
    /// March material-driven natural or mixed convection on the existing
    /// SIMPLEC and backward-Euler thermal owners, retaining native velocity,
    /// pressure and all momentum time levels between steps.
    ///
    /// At each accepted temperature, one original source query supplies density
    /// for `(rho(T)-rho0)*gravity` and viscosity/conductivity for both equations.
    /// The same admitted snapshot drives both operators during the next step.
    /// Density and Cp departures are bounded by `config` at the initial state
    /// and EVERY endpoint, including the last. Inertia, divergence-free volume
    /// flux, heat-capacity advection and sensible storage retain reference rho/Cp.
    /// Source identity, pressure and composition remain the callback's context.
    ///
    /// Material feedback is explicit and FIRST ORDER in time even with BDF2
    /// momentum. No endpoint fixed-point convergence or second-order coupled
    /// accuracy is implied. The caller's time-step ladder determines accuracy.
    /// Only two source epochs are retained; coefficient buffers are reused.
    ///
    /// Solids retain their explicitly authored/frozen heat capacity and
    /// grid-aligned conductivity. Existing contacts, sources, inlet schedules,
    /// fixed-temperature reservoirs and converged endpoint radiation retain
    /// their original assembly and accounting. Source-dependent solid caloric
    /// and conductivity response uses `march_natural_convection_with_material_solids`.
    /// Variable-density fluid enthalpy storage requires its own owner.
    ///
    /// # Errors
    /// Refuses SGS/LVEL, independently authored eddy conductivity or solid
    /// conductivity tables, invalid time/solver inputs, and unsupported source
    /// states. The entry has no separate scalar-beta input. A refused endpoint,
    /// failed solve or cancellation returns no partial result. Source callbacks
    /// must be repeatable equilibrium queries and must not commit physical history.
    #[allow(clippy::too_many_arguments)]
    pub fn march_natural_convection<S: FluidTransport>(
        &self,
        domain: &VoxelDomain,
        solids: &[SolidMaterial],
        flow_config: &SimpleConfig,
        thermal: &ThermalSetup,
        energy_config: &EnergyConfig,
        initial_temperature_k: &[f64],
        config: BoussinesqConfig<3>,
        unsteady: &UnsteadyConfig,
        inlet_schedule: impl Fn(f64) -> f64,
        power_schedule: impl Fn(f64) -> f64,
        radiation: Option<&RadiationConfig>,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<S, MaterialError>,
    ) -> Result<MaterialBuoyancyMarch<S>, MaterialBuoyancyError> {
        self.march_natural_convection_inner(
            domain, solids, flow_config, thermal, energy_config, initial_temperature_k,
            config, unsteady, inlet_schedule, power_schedule, radiation, None, gate, resolve,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn march_natural_convection_inner<S: FluidTransport>(
        &self,
        domain: &VoxelDomain,
        solids: &[SolidMaterial],
        flow_config: &SimpleConfig,
        thermal: &ThermalSetup,
        energy_config: &EnergyConfig,
        initial_temperature_k: &[f64],
        config: BoussinesqConfig<3>,
        unsteady: &UnsteadyConfig,
        inlet_schedule: impl Fn(f64) -> f64,
        power_schedule: impl Fn(f64) -> f64,
        radiation: Option<&RadiationConfig>,
        material_heat: Option<&mut dyn FnMut(&MaterialHeatStep<'_>, &mut MaterialCoefficients)
            -> Result<(Vec<f64>, TransientRecord), MaterialBuoyancyError>>,
        gate: &CancelGate,
        mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<S, MaterialError>,
    ) -> Result<MaterialBuoyancyMarch<S>, MaterialBuoyancyError> {
        poll(gate)?;
        config.validate().map_err(|_| {
            MaterialBuoyancyError::Invalid(
                "gravity and reference density/heat-capacity departure bounds must be admitted",
            )
        })?;
        if initial_temperature_k.len() != domain.cell_count() || domain.fluid_count() == 0 {
            return Err(MaterialBuoyancyError::Invalid(
                "one initial temperature per cell and at least one fluid cell required",
            ));
        }
        if flow_config.turbulence != Turbulence::Laminar
            || !thermal.eddy_conductivity_w_m_k.is_empty()
            || solids
                .iter()
                .any(|solid| !solid.conductivity_table.is_empty())
        {
            return Err(MaterialBuoyancyError::Invalid(
                "material transient feedback requires laminar flow without a second conductivity owner",
            ));
        }
        let compact = compact_owners(domain, &thermal.compact_components, gate)?;
        let mut fields = MaterialCoefficients {
            viscosity_pa_s: vec![self.transport.dynamic_viscosity_pa_s(); domain.cell_count()],
            conductivity_w_m_k: vec![
                [self.transport.thermal_conductivity_w_m_k(); 3];
                domain.cell_count()
            ],
            force_n_m3: vec![[0.0; 3]; domain.cell_count()],
        };
        for (cell, values) in fields.conductivity_w_m_k.iter_mut().enumerate() {
            poll(gate)?;
            if compact.get(cell).is_some_and(|&owner| owner != usize::MAX) {
                // This numeric field slot is inactive: the authored junction
                // links replace every bulk face of a collapsed package.
                continue;
            }
            if let Voxel::Solid(index) = domain.voxel_at(cell) {
                *values = solids
                    .get(usize::from(index))
                    .ok_or(MaterialBuoyancyError::Invalid(
                        "solid voxel refers to an absent material",
                    ))?
                    .axis_conductivity();
            }
        }
        drop(compact);
        let mut applied = Vec::new();
        let mut endpoint = Vec::new();
        let mut refresh = |temperature: &[f64], fields: &mut MaterialCoefficients| {
            let next = samples(
                domain,
                &self.transport,
                temperature,
                config,
                gate,
                &mut resolve,
            )?;
            for sample in &next {
                let transport = sample.material.reference_storage_transport();
                fields.viscosity_pa_s[sample.cell] = transport.dynamic_viscosity_pa_s();
                fields.conductivity_w_m_k[sample.cell] =
                    [transport.thermal_conductivity_w_m_k(); 3];
                fields.force_n_m3[sample.cell] = sample.force_n_m3;
            }
            applied = core::mem::replace(&mut endpoint, next);
            Ok::<(), MaterialBuoyancyError>(())
        };
        let march = march_conjugate_with_material(
            domain,
            &self.properties,
            solids,
            thermal,
            flow_config,
            unsteady,
            initial_temperature_k,
            inlet_schedule,
            power_schedule,
            energy_config,
            radiation,
            &mut fields,
            &mut refresh,
            material_heat,
            gate,
        )?;
        poll(gate)?;
        Ok(MaterialBuoyancyMarch {
            reference: self.transport,
            march,
            applied,
            endpoint,
            applied_time_s: (unsteady.steps - 1) as f64 * unsteady.time_step_s,
            endpoint_time_s: unsteady.steps as f64 * unsteady.time_step_s,
        })
    }
}
