//! Declared chemistry coupled to the original moving-liquid thermal timeline.
//!
//! Formation energy stays in the existing gas inventory. A chemical step at
//! fixed total internal energy changes composition, temperature and pressure;
//! the SAME source-backed gas then exchanges heat with native VOF liquid and
//! finite walls. Original liquid material laws update viscosity, conductivity,
//! buoyancy and any attached electrical response before native momentum evolves.
//! This is a first-order reaction/thermal/flow split over one physical interval.

use crate::cavity_phs::piston::reactive::{
    GasReaction, GasReactionError, GasReactionSolve, GasReactionStep,
};
use fs_exec::CancelGate;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceCellHeat, FreeSurfaceEnthalpyError, MaterialFreeSurfaceEnthalpy3,
    conduction::{
        FreeSurfaceBodyContact, FreeSurfaceBodyInterface, FreeSurfaceConductionControls,
        FreeSurfaceConductionFace, FreeSurfaceGasFace, FreeSurfaceThermalBody,
    },
    viscous::FreeSurfaceViscousConductionReport,
};
use fs_material::thermochemical_gas::energy::{IsochoricGasParcel, parcel::ClosedGasParcel};

/// Original chemical and thermal work limits plus a whole-step admission.
#[derive(Clone, Copy, Debug)]
pub struct ReactiveLiquidGasSolve {
    /// Original reaction budget, including the one shared caloric inversion policy.
    pub reaction: GasReactionSolve,
    /// Original local material conduction/subcycling and energy budget.
    pub thermal: FreeSurfaceConductionControls,
    /// Both chemical intermediate and final gas departures from initial T [K].
    pub maximum_gas_temperature_increment_k: f64,
    /// Actual gas + liquid + wall storage minus authored/electrical/viscous heat [J].
    pub total_energy_tolerance_j: f64,
    /// Deposit original molecular stress power once; numerical losses are never heat.
    pub with_viscous_heating: bool,
}

/// Source, numerical and cancellation failures preserve every original state.
#[derive(Debug)]
pub enum ReactiveLiquidGasError {
    Reaction(GasReactionError),
    Liquid(FreeSurfaceEnthalpyError),
    Invalid(&'static str),
    Cancelled,
}
impl core::fmt::Display for ReactiveLiquidGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reacting gas/moving liquid refused: {self:?}")
    }
}
impl std::error::Error for ReactiveLiquidGasError {}
impl From<GasReactionError> for ReactiveLiquidGasError {
    fn from(error: GasReactionError) -> Self {
        if error == GasReactionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Reaction(error)
        }
    }
}
impl From<FreeSurfaceEnthalpyError> for ReactiveLiquidGasError {
    fn from(error: FreeSurfaceEnthalpyError) -> Self {
        Self::Liquid(error)
    }
}

/// Independent complete-interval inventory check and unchanged owner reports.
#[derive(Clone, Debug)]
pub struct ReactiveLiquidGasReport {
    /// The native flow timestep, shared by chemistry and heat, not added together.
    pub duration_s: f64,
    /// Original adiabatic reaction. Chemical formation energy is NOT a heat source.
    pub reaction: GasReactionStep,
    /// Original finite gas/solid/liquid conduction and native VOF response.
    pub thermal: FreeSurfaceViscousConductionReport,
    /// Actual gas internal-energy change from BEFORE chemistry [J].
    pub gas_internal_energy_change_j: f64,
    /// Actual liquid enthalpy change after native mass/enthalpy transport [J].
    pub liquid_enthalpy_change_j: f64,
    /// Actual sum of finite solid wall enthalpy changes [J].
    pub wall_enthalpy_change_j: f64,
    /// External authored/electrical heat and physical molecular dissipation [J].
    /// Gas/liquid and solid/liquid exchange are internal and do not enter this sum.
    pub supplied_thermal_energy_j: f64,
    /// Sum of actual storage changes minus supplied_thermal_energy_j [J].
    pub combined_energy_residual_j: f64,
}

/// Complete candidate; no accepted caller-owned state is changed by its creation.
pub struct ReactiveLiquidGasCandidate<'s, 'm, 'b> {
    liquid: MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: IsochoricGasParcel,
    bodies: Vec<FreeSurfaceThermalBody<'b>>,
    report: ReactiveLiquidGasReport,
}
impl<'s, 'm, 'b> ReactiveLiquidGasCandidate<'s, 'm, 'b> {
    pub const fn liquid(&self) -> &MaterialFreeSurfaceEnthalpy3<'s, 'm> {
        &self.liquid
    }
    pub const fn gas(&self) -> &IsochoricGasParcel {
        &self.gas
    }
    pub fn bodies(&self) -> &[FreeSurfaceThermalBody<'b>] {
        &self.bodies
    }
    pub const fn report(&self) -> &ReactiveLiquidGasReport {
        &self.report
    }
    /// Publish the accepted original owners together, with their source bindings intact.
    pub fn into_parts(
        self,
    ) -> (
        MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        IsochoricGasParcel,
        Vec<FreeSurfaceThermalBody<'b>>,
        ReactiveLiquidGasReport,
    ) {
        (self.liquid, self.gas, self.bodies, self.report)
    }
}

fn poll(gate: &CancelGate) -> Result<(), ReactiveLiquidGasError> {
    if gate.is_requested() {
        Err(ReactiveLiquidGasError::Cancelled)
    } else {
        Ok(())
    }
}

impl GasReaction {
    /// React a finite atmosphere, exchange heat, and evolve the SAME moving liquid.
    ///
    /// The interval comes from the liquid's retained dx/dt mapping. Chemistry
    /// uses the accepted gas inventory and original NASA-9/transport owners;
    /// thermal exchange starts from its reacted state and cools that same gas.
    /// Reusing the candidate on the next call therefore closes temperature and
    /// composition feedback into rates without resetting to a reference mixture.
    ///
    /// Existing body/interface, gas-film, electrical, viscosity and buoyancy
    /// models remain attached. All source and energy gates must pass before any
    /// state can publish. Gas, native VOF populations/carry, caloric histories,
    /// contact history and finite walls are cloned for that one transaction.
    ///
    /// The gas is homogeneous at fixed mass/volume. Film geometry must describe
    /// actually exposed liquid faces at the start of this interval. This adds
    /// no gas momentum/traction, flame or shock, evaporation, phase conversion,
    /// moving solid, surface-energy or whole-system mechanical energy theorem.
    /// The final account is thermal storage with physical molecular heating;
    /// native VOF numerical kinetic losses never become heat. The chemical
    /// entropy gate remains local to the isolated chemical stage. Refine the
    /// complete timestep to assess this explicit first-order split.
    #[allow(clippy::too_many_arguments)]
    pub fn liquid_candidate<'s, 'm, 'b>(
        &self,
        liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        gas: &IsochoricGasParcel,
        bodies: &[FreeSurfaceThermalBody<'b>],
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        contacts: &[FreeSurfaceBodyContact],
        interfaces: &[FreeSurfaceBodyInterface<'_>],
        exposed: &[FreeSurfaceGasFace],
        controls: ReactiveLiquidGasSolve,
        gate: &CancelGate,
    ) -> Result<ReactiveLiquidGasCandidate<'s, 'm, 'b>, ReactiveLiquidGasError> {
        poll(gate)?;
        if !controls.maximum_gas_temperature_increment_k.is_finite()
            || controls.maximum_gas_temperature_increment_k <= 0.0
            || !controls.total_energy_tolerance_j.is_finite()
            || controls.total_energy_tolerance_j < 0.0
            || exposed.is_empty()
            || exposed.len() > controls.thermal.maximum_faces
        {
            return Err(ReactiveLiquidGasError::Invalid(
                "bounded exposed patches and finite gas temperature/combined energy limits required",
            ));
        }
        let duration_s = liquid.flow().scaling().timestep_s();
        let initial_temperature = gas.state().state().temperature;
        let initial_liquid_h = liquid.total_enthalpy_j()?;
        let initial_wall_h: f64 = bodies.iter().map(FreeSurfaceThermalBody::enthalpy_j).sum();
        let reaction = self.parcel_candidate(
            &ClosedGasParcel::from_isochoric(gas.clone()),
            duration_s,
            controls.reaction,
            || gate.is_requested(),
        )?;
        if (reaction.gas().state().state().temperature - initial_temperature).abs()
            > controls.maximum_gas_temperature_increment_k
        {
            return Err(ReactiveLiquidGasError::Invalid(
                "chemical intermediate exceeds the gas temperature-increment allowance",
            ));
        }
        let (reacted, reaction) = reaction.into_parts();
        let mut next_gas = reacted.into_isochoric();
        let mut next_liquid = liquid.clone();
        let mut next_bodies = bodies.to_vec();
        let thermal = next_liquid.step_with_finite_gas_environment(
            faces,
            heat,
            controls.thermal,
            &mut next_bodies,
            contacts,
            interfaces,
            &mut next_gas,
            exposed,
            controls.reaction.gas,
            controls.with_viscous_heating,
            gate,
        )?;
        let gas_internal_energy_change_j = gas.mass_kg()
            * (next_gas.state().specific_internal_energy_j_kg()
                - gas.state().specific_internal_energy_j_kg());
        let liquid_enthalpy_change_j = next_liquid.total_enthalpy_j()? - initial_liquid_h;
        let wall_enthalpy_change_j = next_bodies
            .iter()
            .map(FreeSurfaceThermalBody::enthalpy_j)
            .sum::<f64>()
            - initial_wall_h;
        let environment =
            thermal
                .conduction
                .environment
                .as_ref()
                .ok_or(ReactiveLiquidGasError::Invalid(
                    "finite atmosphere lost its native heat account",
                ))?;
        let flow = &thermal.conduction.flow;
        let supplied_thermal_energy_j =
            flow.external_heat_j - environment.heat_to_liquid_j + flow.viscous_heat_j;
        let combined_energy_residual_j =
            gas_internal_energy_change_j + liquid_enthalpy_change_j + wall_enthalpy_change_j
                - supplied_thermal_energy_j;
        if [
            gas_internal_energy_change_j,
            liquid_enthalpy_change_j,
            wall_enthalpy_change_j,
            supplied_thermal_energy_j,
            combined_energy_residual_j,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || combined_energy_residual_j.abs() > controls.total_energy_tolerance_j
            || next_gas.mass_kg().to_bits() != gas.mass_kg().to_bits()
            || next_gas.volume_m3().to_bits() != gas.volume_m3().to_bits()
            || (next_gas.state().state().temperature - initial_temperature).abs()
                > controls.maximum_gas_temperature_increment_k
        {
            return Err(ReactiveLiquidGasError::Invalid(
                "reaction/thermal/native-flow interval does not close its original inventories or bounds",
            ));
        }
        poll(gate)?;
        Ok(ReactiveLiquidGasCandidate {
            liquid: next_liquid,
            gas: next_gas,
            bodies: next_bodies,
            report: ReactiveLiquidGasReport {
                duration_s,
                reaction,
                thermal,
                gas_internal_energy_change_j,
                liquid_enthalpy_change_j,
                wall_enthalpy_change_j,
                supplied_thermal_energy_j,
                combined_energy_residual_j,
            },
        })
    }
}
