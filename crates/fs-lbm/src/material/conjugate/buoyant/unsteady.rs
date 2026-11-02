//! Source-material natural and mixed convection on the native transient clock.

use super::{MaterialBuoyancyError, MaterialBuoyancySample, poll, samples};
#[cfg(feature = "enthalpy-coupling")]
use crate::conjugate::unsteady::MaterialViscousHeat;
use crate::conjugate::unsteady::{
    MaterialCoefficients, MaterialEndpointCoupling, MaterialHeatStep, compact_owners,
    march_conjugate_with_material,
};
pub use crate::conjugate::unsteady::{
    MaterialEndpointConfig, MaterialEndpointMarch, MaterialEndpointStepReport,
};
use crate::conjugate::{
    ConjugateMarch, EnergyConfig, FvSgsConfig, RadiationConfig, SimpleConfig, SolidMaterial,
    ThermalSetup, TransientRecord, Turbulence, UnsteadyConfig, Voxel, VoxelDomain,
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
    /// Physical time of the final applied source epoch. Endpoint-coupled
    /// entries use the returned time; explicit entries use the preceding time.
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
            domain,
            solids,
            flow_config,
            thermal,
            energy_config,
            initial_temperature_k,
            config,
            unsteady,
            None,
            #[cfg(feature = "enthalpy-coupling")]
            None,
            inlet_schedule,
            power_schedule,
            radiation,
            None,
            None,
            gate,
            resolve,
        )
    }

    /// Close source viscosity, conductivity and density-derived buoyancy at
    /// each native thermal endpoint, retaining the original momentum history.
    ///
    /// Every bounded outer pass solves the same BE/BDF2 momentum interval and
    /// backward-Euler heat interval. Previous/older velocity levels and the old
    /// accepted temperature stay fixed throughout the nonlinear iteration;
    /// physical records, averages and time levels advance once after acceptance.
    /// Source queries at the actual unrelaxed candidate must close temperature,
    /// mu/k and every force component independently. Relaxation changes only
    /// the next source-temperature iterate. A sharp admissible material change
    /// cannot pass solely because its temperature increment is small.
    ///
    /// This removes the explicit source lag but coupled time accuracy remains
    /// FIRST ORDER because heat uses backward Euler. Reference rho/Cp still
    /// owns mass, inertia, advection and storage. Use the source-enthalpy implicit
    /// entry when storage and heat transport require the original h(T) law.
    /// Original solids, contacts, compact capacities, schedules and endpoint
    /// radiation retain their owners. Every discarded momentum/energy/source
    /// trial consumes the whole-march budgets in `endpoint_controls`.
    ///
    /// # Errors
    /// Unsupported source state, native solve failure, exhausted work, failure
    /// of any unrelaxed closure test or cancellation returns no partial march.
    /// This laminar entry does not select SGS, wall closures or viscous heating.
    #[allow(clippy::too_many_arguments)]
    pub fn march_natural_convection_implicit<S: FluidTransport>(
        &self,
        domain: &VoxelDomain,
        solids: &[SolidMaterial],
        flow_config: &SimpleConfig,
        thermal: &ThermalSetup,
        energy_config: &EnergyConfig,
        initial_temperature_k: &[f64],
        config: BoussinesqConfig<3>,
        unsteady: &UnsteadyConfig,
        endpoint_controls: MaterialEndpointConfig,
        inlet_schedule: impl Fn(f64) -> f64,
        power_schedule: impl Fn(f64) -> f64,
        radiation: Option<&RadiationConfig>,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<S, MaterialError>,
    ) -> Result<MaterialEndpointMarch<MaterialBuoyancyMarch<S>>, MaterialBuoyancyError> {
        let mut endpoint = MaterialEndpointCoupling::new(endpoint_controls)?;
        let solution = self.march_natural_convection_inner(
            domain,
            solids,
            flow_config,
            thermal,
            energy_config,
            initial_temperature_k,
            config,
            unsteady,
            None,
            #[cfg(feature = "enthalpy-coupling")]
            None,
            inlet_schedule,
            power_schedule,
            radiation,
            None,
            Some(&mut endpoint),
            gate,
            resolve,
        )?;
        poll(gate)?;
        Ok(endpoint.finish(solution))
    }

    /// March source-dependent molecular viscosity with explicitly selected
    /// transient SGS transport on the original momentum and thermal clock.
    ///
    /// The source supplies molecular mu(T), k(T) and buoyancy at each accepted
    /// thermal endpoint. Every nonlinear momentum sweep recomputes nu_t from
    /// the resolved velocity and checks its ratio to the LOCAL molecular nu.
    /// Their sum supplies the full symmetric stress; material jumps retain
    /// harmonic edge transmission. No eddy coefficient is folded back into
    /// the source state or compounded on a subsequent refresh. No-slip wall
    /// traction still uses the molecular half-cell viscosity.
    ///
    /// Projected nu_t adds one gradient heat flux with the original reference
    /// rho*Cp and Pr_t=0.9, separately from source molecular conductivity.
    /// Storage, advection and the density/Cp departure bounds are those of
    /// `march_natural_convection`. The source feedback is first order in time;
    /// refinement must resolve both material gradients and the wall layer.
    /// No molecular or modeled dissipation is relabelled as a heating source.
    /// This is an exploratory algebraic closure, not a turbulent-DNS claim.
    ///
    /// # Errors
    /// The original source/solver refusals, invalid SGS controls, acceptance
    /// use, a local eddy/molecular ratio above its explicit cap, or LVEL/another
    /// authored eddy-conductivity field. Failure returns no partial march.
    #[allow(clippy::too_many_arguments)]
    pub fn march_natural_convection_with_sgs<S: FluidTransport>(
        &self,
        domain: &VoxelDomain,
        solids: &[SolidMaterial],
        flow_config: &SimpleConfig,
        thermal: &ThermalSetup,
        energy_config: &EnergyConfig,
        initial_temperature_k: &[f64],
        config: BoussinesqConfig<3>,
        unsteady: &UnsteadyConfig,
        sgs: FvSgsConfig,
        inlet_schedule: impl Fn(f64) -> f64,
        power_schedule: impl Fn(f64) -> f64,
        radiation: Option<&RadiationConfig>,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<S, MaterialError>,
    ) -> Result<MaterialBuoyancyMarch<S>, MaterialBuoyancyError> {
        self.march_natural_convection_inner(
            domain,
            solids,
            flow_config,
            thermal,
            energy_config,
            initial_temperature_k,
            config,
            unsteady,
            Some(sgs),
            #[cfg(feature = "enthalpy-coupling")]
            None,
            inlet_schedule,
            power_schedule,
            radiation,
            None,
            None,
            gate,
            resolve,
        )
    }

    /// Close material/buoyancy feedback with explicitly selected native SGS.
    ///
    /// Every material trial solves the original BE/BDF2 momentum equations
    /// from the same accepted time levels. Its selected native SGS field is
    /// recomputed from the resolved velocity, including the final pressure
    /// projection, with the original local eddy/molecular viscosity cap.
    /// That current trial supplies `rho0*Cp0*nu_t/0.9` to the native heat
    /// operator exactly once, separately from the original molecular k(T).
    /// Neither old SGS conductivity nor a compounded material coefficient is
    /// carried into a new trial. Unrelaxed endpoint temperature, mu/k and
    /// density-derived force must independently meet `endpoint_controls`.
    ///
    /// Histories, averages and physical records advance once on acceptance;
    /// all discarded native/source/heat work consumes the same explicit
    /// budgets as `march_natural_convection_implicit`. Thermal accuracy is
    /// first order. This algebraic SGS composition is exploratory: it does
    /// not establish resolved turbulent spectra or a thermal wall layer.
    /// No molecular or modeled dissipation is deposited as heat.
    ///
    /// # Errors
    /// The original implicit/source/native refusals, inadmissible SGS or
    /// acceptance use, a local eddy/molecular cap violation, equilibrium LVEL,
    /// or another authored eddy-conductivity field. No partial march is returned.
    #[allow(clippy::too_many_arguments)]
    pub fn march_natural_convection_with_sgs_implicit<S: FluidTransport>(
        &self,
        domain: &VoxelDomain,
        solids: &[SolidMaterial],
        flow_config: &SimpleConfig,
        thermal: &ThermalSetup,
        energy_config: &EnergyConfig,
        initial_temperature_k: &[f64],
        config: BoussinesqConfig<3>,
        unsteady: &UnsteadyConfig,
        sgs: FvSgsConfig,
        endpoint_controls: MaterialEndpointConfig,
        inlet_schedule: impl Fn(f64) -> f64,
        power_schedule: impl Fn(f64) -> f64,
        radiation: Option<&RadiationConfig>,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<S, MaterialError>,
    ) -> Result<MaterialEndpointMarch<MaterialBuoyancyMarch<S>>, MaterialBuoyancyError> {
        let mut endpoint = MaterialEndpointCoupling::new(endpoint_controls)?;
        let solution = self.march_natural_convection_inner(
            domain,
            solids,
            flow_config,
            thermal,
            energy_config,
            initial_temperature_k,
            config,
            unsteady,
            Some(sgs),
            #[cfg(feature = "enthalpy-coupling")]
            None,
            inlet_schedule,
            power_schedule,
            radiation,
            None,
            Some(&mut endpoint),
            gate,
            resolve,
        )?;
        poll(gate)?;
        Ok(endpoint.finish(solution))
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
        sgs: Option<FvSgsConfig>,
        #[cfg(feature = "enthalpy-coupling")] viscous_heat: Option<&mut MaterialViscousHeat>,
        inlet_schedule: impl Fn(f64) -> f64,
        power_schedule: impl Fn(f64) -> f64,
        radiation: Option<&RadiationConfig>,
        material_heat: Option<
            &mut dyn FnMut(
                &MaterialHeatStep<'_>,
                &mut MaterialCoefficients,
            )
                -> Result<(Vec<f64>, TransientRecord), MaterialBuoyancyError>,
        >,
        endpoint_coupling: Option<&mut MaterialEndpointCoupling>,
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
        let implicit = endpoint_coupling.is_some();
        let march = march_conjugate_with_material(
            domain,
            &self.properties,
            solids,
            thermal,
            flow_config,
            unsteady,
            sgs,
            #[cfg(feature = "enthalpy-coupling")]
            viscous_heat,
            initial_temperature_k,
            inlet_schedule,
            power_schedule,
            energy_config,
            radiation,
            &mut fields,
            &mut refresh,
            material_heat,
            endpoint_coupling,
            gate,
        )?;
        poll(gate)?;
        Ok(MaterialBuoyancyMarch {
            reference: self.transport,
            march,
            applied,
            endpoint,
            applied_time_s: (unsteady.steps - usize::from(!implicit)) as f64 * unsteady.time_step_s,
            endpoint_time_s: unsteady.steps as f64 * unsteady.time_step_s,
        })
    }
}
