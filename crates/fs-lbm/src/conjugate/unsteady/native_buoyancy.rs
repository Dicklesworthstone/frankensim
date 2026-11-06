//! Implicit native Boussinesq feedback on the existing momentum/heat march.
//!
//! The carrier's conductivity field participates only in endpoint closure. Heat
//! retains its integrated solid k(T), contact, wall and compact-storage laws.

use super::{
    Boussinesq, ChtError, ConjugateMarch, EnergyConfig, FluidProperties, FvSgsConfig, MaterialCoefficients,
    MaterialEndpointConfig, MaterialEndpointCoupling, MaterialEndpointMarch, MaterialHeatStep,
    RadiationConfig, SimpleConfig, SolidMaterial, ThermalSetup, Turbulence, UnsteadyConfig,
    VoxelDomain, finite, finite_positive, march_conjugate_with_material, poll,
};
use crate::conjugate::Voxel;
use fs_exec::CancelGate;

/// March laminar native flow and heat with buoyancy closed at each endpoint.
///
/// Unlike the lagged [`super::march_conjugate`] entry, each physical interval
/// iterates the original momentum and energy equations until the temperature
/// and Boussinesq force agree within `endpoint` tolerances. Previous and older
/// velocity levels and the old temperature remain fixed during these trials;
/// only accepted steps enter histories and averages. Energy remains backward
/// Euler; the selected native BE/BDF2 momentum history is unchanged.
///
/// Fluid properties and the signed thermal expansion coefficient retain the
/// ordinary constant-property Boussinesq interpretation. Native solid k(T)
/// tables, integrated interface and reservoir laws, compact storage, radiation
/// and the selected temperature convection scheme keep their original owners.
/// The coefficient refresh count includes native force and cell-conductivity
/// evaluations; integrated heat remains in its original inner solve.
///
/// # Errors
///
/// The original march's input, convergence and cancellation refusals; an
/// equilibrium turbulence closure or authored eddy-conductivity field; invalid
/// endpoint controls; or exhaustion of a nonlinear/whole-march work limit.
/// An unclosed candidate never becomes an accepted physical time level.
#[allow(clippy::too_many_arguments)]
pub fn march_conjugate_implicit_buoyancy(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    config: &SimpleConfig,
    unsteady: &UnsteadyConfig,
    buoyancy: &Boussinesq,
    endpoint: MaterialEndpointConfig,
    initial_temperature: &[f64],
    inlet_schedule: impl Fn(f64) -> f64,
    power_schedule: impl Fn(f64) -> f64,
    energy: &EnergyConfig,
    radiation: Option<&RadiationConfig>,
    gate: &CancelGate,
) -> Result<MaterialEndpointMarch<ConjugateMarch>, ChtError> {
    march_with_sgs(
        domain, fluid, solids, setup, config, unsteady, None, buoyancy,
        endpoint, initial_temperature, inlet_schedule, power_schedule,
        energy, radiation, gate,
    )
}

/// Converge native buoyancy and heat with explicit transient SGS transport.
///
/// Each nonlinear trial solves the original momentum equations with the
/// selected Smagorinsky, WALE or Vreman model. Its projected velocity supplies
/// one heat-diffusion field `rho * cp * nu_t / Pr_t`, using the independently
/// declared [`SimpleConfig::turbulent_prandtl`]. That field is replaced on every
/// trial; it is never accumulated into molecular conductivity. Wall and
/// solid-interface heat transfer retain their native molecular laws.
///
/// Previous/older momentum states and the old temperature stay fixed through
/// the coupling trials as in [`march_conjugate_implicit_buoyancy`]. Native
/// nonlinear solid conductivity, compact storage and radiation remain inside
/// the same heat solve. A zero model coefficient exactly recovers the
/// laminar implicit entry. This is an explicit exploratory SGS closure, with
/// its original viscosity cap and use policy; it adds no SGS dissipation heat
/// or thermal wall model and establishes no resolved turbulent spectrum.
///
/// # Errors
///
/// The refusals of [`march_conjugate_implicit_buoyancy`] and the original SGS
/// controls, including an incompatible closure, viscosity cap or use policy.
#[allow(clippy::too_many_arguments)]
pub fn march_conjugate_with_sgs_implicit_buoyancy(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    config: &SimpleConfig,
    unsteady: &UnsteadyConfig,
    sgs: FvSgsConfig,
    buoyancy: &Boussinesq,
    endpoint: MaterialEndpointConfig,
    initial_temperature: &[f64],
    inlet_schedule: impl Fn(f64) -> f64,
    power_schedule: impl Fn(f64) -> f64,
    energy: &EnergyConfig,
    radiation: Option<&RadiationConfig>,
    gate: &CancelGate,
) -> Result<MaterialEndpointMarch<ConjugateMarch>, ChtError> {
    march_with_sgs(
        domain, fluid, solids, setup, config, unsteady, Some(sgs), buoyancy,
        endpoint, initial_temperature, inlet_schedule, power_schedule,
        energy, radiation, gate,
    )
}

#[allow(clippy::too_many_arguments)]
fn march_with_sgs(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    config: &SimpleConfig,
    unsteady: &UnsteadyConfig,
    sgs: Option<FvSgsConfig>,
    buoyancy: &Boussinesq,
    endpoint: MaterialEndpointConfig,
    initial_temperature: &[f64],
    inlet_schedule: impl Fn(f64) -> f64,
    power_schedule: impl Fn(f64) -> f64,
    energy: &EnergyConfig,
    radiation: Option<&RadiationConfig>,
    gate: &CancelGate,
) -> Result<MaterialEndpointMarch<ConjugateMarch>, ChtError> {
    poll(gate)?;
    if let Some(sgs) = sgs {
        sgs.validate()?;
    }
    fluid.validate()?;
    domain.check_materials(solids.len())?;
    if config.turbulence != Turbulence::Laminar || !setup.eddy_conductivity_w_m_k.is_empty() {
        return Err(ChtError::InvalidInput {
            field: "unsteady.implicit_buoyancy",
            reason: "native implicit buoyancy excludes equilibrium turbulence and authored eddy conductivity".into(),
        });
    }
    finite("buoyancy.expansion_per_k", buoyancy.expansion_per_k)?;
    finite(
        "buoyancy.reference_temperature_k",
        buoyancy.reference_temperature_k,
    )?;
    for gravity in buoyancy.gravity_m_s2 {
        finite("buoyancy.gravity_m_s2", gravity)?;
    }
    let viscosity = fluid.density_kg_m3 * fluid.kinematic_viscosity_m2_s;
    finite_positive("fluid.dynamic_viscosity_pa_s", viscosity)?;
    let cells = domain.cell_count();
    let mut fields = MaterialCoefficients {
        viscosity_pa_s: vec![viscosity; cells],
        conductivity_w_m_k: vec![[fluid.conductivity_w_m_k; 3]; cells],
        force_n_m3: vec![[0.0; 3]; cells],
    };
    let mut refresh = |temperature: &[f64], fields: &mut MaterialCoefficients| {
        refresh_native(domain, fluid, solids, buoyancy, temperature, fields, gate)
    };
    // Never pass the cell-centre bookkeeping field as the heat operator:
    // doing so would replace the native integrated k(T) face and wall laws.
    let mut heat =
        |interval: &MaterialHeatStep<'_>, _: &mut MaterialCoefficients| interval.solve_native();
    let mut feedback = MaterialEndpointCoupling::new(endpoint)?;
    let solution = march_conjugate_with_material::<ChtError>(
        domain,
        fluid,
        solids,
        setup,
        config,
        unsteady,
        sgs,
        #[cfg(feature = "enthalpy-coupling")]
        None,
        initial_temperature,
        inlet_schedule,
        power_schedule,
        energy,
        radiation,
        &mut fields,
        &mut refresh,
        Some(&mut heat),
        Some(&mut feedback),
        gate,
    )?;
    Ok(feedback.finish(solution))
}

fn refresh_native(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    buoyancy: &Boussinesq,
    temperature: &[f64],
    fields: &mut MaterialCoefficients,
    gate: &CancelGate,
) -> Result<(), ChtError> {
    for (cell, &temperature) in temperature.iter().enumerate() {
        if cell.is_multiple_of(256) {
            poll(gate)?;
        }
        finite("unsteady.initial_temperature", temperature)?;
        fields.force_n_m3[cell] = [0.0; 3];
        fields.conductivity_w_m_k[cell] = match domain.voxel_at(cell) {
            Voxel::Fluid => {
                if buoyancy.expansion_per_k != 0.0
                    && buoyancy.gravity_m_s2.iter().any(|&g| g != 0.0)
                {
                    let scale = -fluid.density_kg_m3
                        * buoyancy.expansion_per_k
                        * (temperature - buoyancy.reference_temperature_k);
                    fields.force_n_m3[cell] = buoyancy.gravity_m_s2.map(|gravity| scale * gravity);
                    for force in fields.force_n_m3[cell] {
                        finite("buoyancy.force_n_m3", force)?;
                    }
                }
                [fluid.conductivity_w_m_k; 3]
            }
            Voxel::Solid(material) => {
                let solid = &solids[usize::from(material)];
                let factor = solid.conductivity_factor(temperature);
                let axes = solid.axis_conductivity().map(|k| k * factor);
                for k in axes {
                    finite_positive("solid.endpoint_conductivity_w_m_k", k)?;
                }
                axes
            }
        };
    }
    Ok(())
}
