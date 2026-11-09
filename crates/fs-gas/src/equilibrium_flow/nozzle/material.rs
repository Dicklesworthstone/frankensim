//! Original material publication of a reacting converging-nozzle endpoint.
//!
//! The shared SP search establishes the physical branch. The existing material
//! HP path publishes the same source-relative enthalpy drop with original
//! species cards, transport and feed mass. Its numerical inverse cannot replace
//! nozzle acceptance: actual material energy, density and flow, chemical entropy
//! and reacting sonic residual are checked again before the stream is returned.

use super::{
    EquilibriumNozzleControl, EquilibriumNozzleError, EquilibriumNozzleFlow,
    equilibrium_nozzle_flow,
};
use fs_matdb::QueryPoint;
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::equilibrium::pressure::{
        EquilibriumGasEnthalpyPath, GasEquilibriumEnthalpyError, GasEquilibriumEnthalpySolve,
        GasEquilibriumEnthalpyState,
    },
};
use fs_math::{dd::Dd, det};
use fs_qty::{Amount, Area, Energy, MassFlowRate, Pressure};
use fs_thermochem::{SpeciesId, equilibrium::composition::IdealGasEquilibriumV1};

/// Native flow controls plus independent original-material publication limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialEquilibriumNozzleControl {
    pub nozzle: EquilibriumNozzleControl,
    /// Original HP work/source bounds and material mass/H/freezing tolerances.
    /// Its Gibbs controls must equal the nozzle SP controls exactly.
    pub material: GasEquilibriumEnthalpySolve,
    /// Maximum absolute material/native density and mass-flux ratio minus one.
    /// In (0,0.01); this numerical endpoint allowance is not a discharge loss.
    pub relative_endpoint_tolerance: f64,
}

/// Actual material flow compared with its original reacting SP branch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialEquilibriumNozzleBalance {
    /// Direct material Delta H divided by exactly retained feed mass [J/kg].
    pub specific_enthalpy_change_j_kg: f64,
    /// Actual material Delta h plus actual v^2/2 [J/kg].
    pub stagnation_enthalpy_residual_j_kg: f64,
    /// Chemical entropy change of the actual published HP endpoint [J/K].
    pub entropy_residual_j_k: f64,
    /// Actual material density divided by the accepted native SP density minus one.
    pub relative_density_residual: f64,
    /// Actual material mass flux divided by the native SP flux minus one.
    /// Exactly zero for the original zero-flow endpoint.
    pub relative_mass_flux_residual: f64,
    /// Actual reacting c^2 minus actual material v^2; absent at zero flow.
    pub sonic_residual_m2_s2: Option<f64>,
}

/// Both original thermodynamic and material endpoints, without a second EOS.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialEquilibriumNozzleFlow {
    /// Native SP search, source evidence and branch/pressure-work accounts.
    pub isentropic: EquilibriumNozzleFlow,
    /// Source-bearing material reference, with all original transport selections.
    pub inlet: GasEquilibriumEnthalpyState,
    /// Actual accepted material products and their chemical-equilibrium response.
    pub outlet: GasEquilibriumEnthalpyState,
    /// Actual density times velocity times the authored effective area.
    pub mass_flow: MassFlowRate,
    pub mass_flux_kg_m2_s: f64,
    pub velocity_m_s: f64,
    /// Actual velocity divided by the published chemical-equilibrium sound speed.
    /// Frozen material sound speed remains separate and is never used for choking.
    pub mach_number: f64,
    pub balance: MaterialEquilibriumNozzleBalance,
    /// One original material HP reference plus at most one final HP query.
    pub material_pressure_evaluations: usize,
    /// Total actual Gibbs work across SP search and both material stages.
    pub composition_iterations: usize,
    /// Total actual temperature iterations across SP and HP queries.
    pub temperature_iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MaterialEquilibriumNozzleError {
    Invalid(&'static str),
    Nozzle(EquilibriumNozzleError),
    Material(GasEquilibriumEnthalpyError),
    /// Actual material/chemical quantities fail the original flow tolerances.
    Endpoint(MaterialEquilibriumNozzleBalance),
    Cancelled,
}
impl core::fmt::Display for MaterialEquilibriumNozzleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reacting material converging nozzle: {self:?}")
    }
}
impl std::error::Error for MaterialEquilibriumNozzleError {}
impl From<EquilibriumNozzleError> for MaterialEquilibriumNozzleError {
    fn from(error: EquilibriumNozzleError) -> Self {
        if matches!(error, EquilibriumNozzleError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Nozzle(error)
        }
    }
}
impl From<GasEquilibriumEnthalpyError> for MaterialEquilibriumNozzleError {
    fn from(error: GasEquilibriumEnthalpyError) -> Self {
        if matches!(error, GasEquilibriumEnthalpyError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Material(error)
        }
    }
}
type Error = MaterialEquilibriumNozzleError;

/// Resolve pressure-driven flow into original-owner material product parcels.
///
/// The context's original typed T,p define the stagnation equilibrium. The bank
/// must contain every eligible species with the exact original thermochemical
/// sources. Actual products retain the source domain and neutral-transport
/// admission of the existing material HP path. No reservoir is depleted.
///
/// The SP search is bounded by `nozzle`; one HP preparation and at most one HP
/// query additionally inherit `material`'s original work limits. The final HP
/// chemical entropy must meet the original SP tolerance, and its actual reacting
/// sonic residual must meet the original nozzle criterion. Tighter HP tolerances
/// may be necessary; none of these acceptance gates is relaxed automatically.
pub fn equilibrium_material_nozzle_flow(
    bank: &ThermochemicalGasModel,
    chemistry: &IdealGasEquilibriumV1,
    context: &QueryPoint,
    feed: &[(SpeciesId, Amount)],
    back_pressure: Pressure,
    effective_area: Area,
    controls: MaterialEquilibriumNozzleControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MaterialEquilibriumNozzleFlow, Error> {
    poll(&mut cancelled)?;
    controls.nozzle.validate()?;
    if controls.material.equilibrium.equilibrium != controls.nozzle.isentropic.equilibrium
        || !controls.relative_endpoint_tolerance.is_finite()
        || !(0.0..0.01).contains(&controls.relative_endpoint_tolerance)
        || controls.relative_endpoint_tolerance == 0.0
    {
        return Err(Error::Invalid(
            "identical original Gibbs policies and a positive bounded endpoint allowance required",
        ));
    }
    let mut path = EquilibriumGasEnthalpyPath::new(
        bank,
        chemistry,
        context,
        feed,
        controls.material,
        &mut cancelled,
    )?;
    let inlet = path.reference_state().clone();
    let initial = inlet.equilibrium().equilibrium();
    let isentropic = equilibrium_nozzle_flow(
        chemistry,
        initial.temperature(),
        initial.pressure(),
        feed,
        back_pressure,
        effective_area,
        controls.nozzle,
        &mut cancelled,
    )?;
    let mass = inlet.balance().mass_kg;
    if inlet.equilibrium() != &isentropic.inlet || mass != isentropic.feed_mass_kg {
        return Err(Error::Invalid(
            "material HP and nozzle SP must retain the identical original equilibrium inlet",
        ));
    }
    let moving = isentropic.velocity_m_s > 0.0;
    let target = isentropic.sample.specific_enthalpy_change_j_kg * mass;
    if !target.is_finite() || (moving && target >= 0.0) {
        return Err(Error::Invalid(
            "material nozzle enthalpy target is not representable",
        ));
    }
    let outlet = if moving {
        path.sample(
            isentropic.sample.pressure,
            Energy::new(target),
            &mut cancelled,
        )?
        .clone()
    } else {
        inlet.clone()
    };
    let sample = outlet.native_sample();
    if sample.pressure != isentropic.sample.pressure || outlet.balance().mass_kg != mass {
        return Err(Error::Invalid(
            "material outlet lost the accepted pressure or original feed mass",
        ));
    }
    let dh = outlet.balance().enthalpy_change_j / mass;
    let v2 = -2.0 * dh;
    if !v2.is_finite() || v2 < 0.0 || (moving && v2 == 0.0) {
        return Err(Error::Invalid(
            "material expansion enthalpy does not resolve positive kinetic energy",
        ));
    }
    let velocity = det::sqrt(v2);
    let density = outlet.parcel().state().state().density;
    let flux = density * velocity;
    let rate = effective_area.value() * flux;
    let sound = sample.response.sound_speed_m_s;
    let mach = velocity / sound;
    let sonic = sound * sound - v2;
    let residual = (Dd::from_f64(dh)
        + Dd::from_f64(0.5) * Dd::from_f64(velocity) * Dd::from_f64(velocity))
    .to_f64();
    let balance = MaterialEquilibriumNozzleBalance {
        specific_enthalpy_change_j_kg: dh,
        stagnation_enthalpy_residual_j_kg: residual,
        entropy_residual_j_k: sample.entropy_change_j_k,
        relative_density_residual: density / isentropic.sample.response.density_kg_m3 - 1.0,
        relative_mass_flux_residual: if moving {
            flux / isentropic.mass_flux_kg_m2_s - 1.0
        } else {
            0.0
        },
        sonic_residual_m2_s2: moving.then_some(sonic),
    };
    let energy_tolerance = controls.nozzle.absolute_energy_tolerance_j_kg
        + controls.nozzle.relative_energy_tolerance * dh.abs();
    let sonic_tolerance = controls.nozzle.absolute_sonic_tolerance_m2_s2
        + controls.nozzle.relative_sonic_tolerance * (sound * sound).max(v2);
    if !density.is_finite()
        || density <= 0.0
        || !sound.is_finite()
        || sound <= 0.0
        || !flux.is_finite()
        || !rate.is_finite()
        || !mach.is_finite()
        || (moving && (flux == 0.0 || rate == 0.0))
        || !residual.is_finite()
        || !energy_tolerance.is_finite()
        || residual.abs() > energy_tolerance
        || !balance.entropy_residual_j_k.is_finite()
        || balance.entropy_residual_j_k.abs() > controls.nozzle.isentropic.entropy_tolerance_j_k
        || !balance.relative_density_residual.is_finite()
        || balance.relative_density_residual.abs() > controls.relative_endpoint_tolerance
        || !balance.relative_mass_flux_residual.is_finite()
        || balance.relative_mass_flux_residual.abs() > controls.relative_endpoint_tolerance
        || !sonic.is_finite()
        || !sonic_tolerance.is_finite()
        || (isentropic.choked && sonic.abs() > sonic_tolerance)
        || (moving && !isentropic.choked && sonic < 0.0)
    {
        return Err(Error::Endpoint(balance));
    }
    let composition_iterations = isentropic.composition_iterations
        + inlet.equilibrium().equilibrium().composition_iterations()
        + if moving {
            sample.composition_iterations
        } else {
            0
        };
    let temperature_iterations = isentropic.temperature_iterations
        + if moving {
            sample.temperature_iterations
        } else {
            0
        };
    poll(&mut cancelled)?;
    Ok(MaterialEquilibriumNozzleFlow {
        isentropic,
        inlet,
        outlet,
        mass_flow: MassFlowRate::new(rate),
        mass_flux_kg_m2_s: flux,
        velocity_m_s: velocity,
        mach_number: mach,
        balance,
        material_pressure_evaluations: if moving { 2 } else { 1 },
        composition_iterations,
        temperature_iterations,
    })
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
