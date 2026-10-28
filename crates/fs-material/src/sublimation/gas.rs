//! Nonisothermal molecular exchange between solids and the original gas.
//!
//! Both temperatures evolve. The shared molecular owner supplies outgoing and
//! incoming mass fluxes at their respective temperatures; the original NASA-9
//! species supplies the energy of each stream. Maxwellian molecules crossing
//! a surface carry h(T)-R_specific*T/2, including the original energy datum.
//! The flux-weighted translational mean is 2RT rather than the bulk 3RT/2:
//! Nathanson (2021), sections 1 and 3,
//! <https://doi.org/10.1007/978-3-030-63963-1_27>.
//!
//! Assumes equilibrium internal modes, energy/angle-independent equal sticking
//! probabilities, specular reflection of nonsticking molecules, and a local
//! free-molecular vapor. Carrier/nonsticking molecules supply no additional
//! thermal accommodation. There is no gas film, Knudsen-layer flow, jet energy,
//! external heat, solid pressure work, or moving-mesh update in this operation.
//!
//! Rates are explicit, frozen at the initial temperatures. BOTH gross streams
//! are inventory-bounded, including when their net mass is zero. Refine the
//! physical timestep to resolve evolving temperatures and backpressure.

/// Joint exchange of heterogeneous surfaces against one evolving gas snapshot.
pub mod chamber;
/// Coupled original heat histories and one evolving thermochemical gas.
#[path = "gas/timeline.rs"]
pub mod timeline;

use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{
        GasEnergyError,
        parcel::{
            ClosedGasParcel,
            transfer::species::boundary::{
                GasSpeciesBoundaryCandidate, GasSpeciesBoundarySolve, SpeciesBoundaryTransfer,
            },
        },
    },
};
use fs_thermal::sublimation::{
    SolidSurfaceState, SublimationError, SublimationLaw, SublimationProperties, VaporReservoir,
    molecular::molecular_mass_flux,
};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};

/// Explicit caller binding; neither chemical identity nor datum follows from a name.
#[derive(Clone, Debug, PartialEq)]
pub struct GasSurfaceBinding {
    /// Canonical gas species explicitly associated with the sublimation law.
    pub species: SpeciesId,
    pub accommodation_coefficient: f64,
    /// Added to the original solid h to express storage in the gas species datum.
    /// This is a fixed reference offset, not extra latent heat or supplied heat.
    pub solid_enthalpy_offset_j_kg: f64,
}

/// Existing gas reconstruction plus explicit phase-coupling admission budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasSurfaceSolve {
    pub gas: GasSpeciesBoundarySolve,
    /// Each gross withdrawal / its own initial solid or vapor inventory, in (0,1).
    /// This controls counterflow stiffness even at zero net mass.
    pub maximum_gross_inventory_fraction: f64,
    /// Absolute agreement required for h_v(T)-(h_s+offset)=L(T), at BOTH endpoints.
    /// A reference offset cannot repair incompatible temperature dependence.
    pub latent_compatibility_tolerance_j_kg: f64,
}

/// Original-owner refusal or a failed physical compatibility/balance check.
#[derive(Clone, Debug, PartialEq)]
pub enum GasSurfaceError {
    Invalid(&'static str),
    /// Frozen gross rates exceed an original inventory allowance. A caller may
    /// retry the entire interval at or below this conservative duration.
    InventoryLimit {
        maximum_duration_s: f64,
    },
    Surface(SublimationError),
    Gas(GasEnergyError),
    CaloricMismatch {
        temperature_k: f64,
        residual_j_kg: f64,
        tolerance_j_kg: f64,
    },
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for GasSurfaceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonisothermal solid/gas exchange: {self:?}")
    }
}
impl std::error::Error for GasSurfaceError {}
impl From<SublimationError> for GasSurfaceError {
    fn from(value: SublimationError) -> Self {
        Self::Surface(value)
    }
}
impl From<GasEnergyError> for GasSurfaceError {
    fn from(value: GasEnergyError) -> Self {
        Self::Gas(value)
    }
}

/// Accounts for the SAME accepted molecular streams, in coherent SI units.
#[derive(Clone, Debug, PartialEq)]
pub struct GasSurfaceStep {
    pub outgoing_mass_kg: f64,
    pub incoming_mass_kg: f64,
    /// Positive sublimation, negative deposition.
    pub transferred_mass_kg: f64,
    pub surface_temperature_k: [f64; 2],
    pub gas_temperature_k: [f64; 2],
    pub vapor_partial_pressure_pa: [f64; 2],
    pub outgoing_specific_energy_j_kg: f64,
    /// None when no incoming molecules require this caloric query.
    pub incoming_specific_energy_j_kg: Option<f64>,
    /// Complete signed energy delivered to the gas; already includes phase energy.
    pub energy_into_gas_j: f64,
    /// Actual solid energy change after applying the fixed datum offset.
    pub solid_energy_change_j: f64,
    pub gas_energy_change_j: f64,
    pub energy_residual_j: f64,
    pub mass_residual_kg: f64,
    /// Positive thickness loss, negative frost growth, using endpoint density.
    pub recession_m: f64,
}

/// Publish the solid and gas together after the surrounding coupling accepts.
#[derive(Clone, Debug, PartialEq)]
pub struct GasSurfaceCandidate<'catalog> {
    solid: SolidSurfaceState,
    gas: GasSpeciesBoundaryCandidate<'catalog>,
    step: GasSurfaceStep,
}
impl<'catalog> GasSurfaceCandidate<'catalog> {
    #[must_use]
    pub const fn solid(&self) -> SolidSurfaceState {
        self.solid
    }
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.gas.gas()
    }
    #[must_use]
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.gas.catalog()
    }
    #[must_use]
    pub const fn gas_exchange(&self) -> &GasSpeciesBoundaryCandidate<'catalog> {
        &self.gas
    }
    #[must_use]
    pub const fn step(&self) -> &GasSurfaceStep {
        &self.step
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        SolidSurfaceState,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasSurfaceStep,
    ) {
        let (gas, catalog, _) = self.gas.into_parts();
        (self.solid, gas, catalog, self.step)
    }
}

fn invalid(message: &'static str) -> GasSurfaceError {
    GasSurfaceError::Invalid(message)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasSurfaceError> {
    if cancelled() {
        Err(SublimationError::Cancelled.into())
    } else {
        Ok(())
    }
}

// Admission, not a constitutive law: the supplied original chart still owns
// all values, domains and temperature dependence.
fn phase_at<L: SublimationLaw + ?Sized>(
    law: &L,
    temperature: f64,
    molar_mass: f64,
) -> Result<SublimationProperties, GasSurfaceError> {
    let [lower, upper] = law.temperature_range_k();
    if !positive(lower) || !positive(upper) || lower > upper {
        return Err(invalid(
            "invalid solid/vapor constitutive temperature range",
        ));
    }
    if !temperature.is_finite() || temperature < lower || temperature > upper {
        return Err(SublimationError::TemperatureOutOfRange {
            temperature_k: temperature,
            minimum_k: lower,
            maximum_k: upper,
        }
        .into());
    }
    let properties = law.properties(temperature)?;
    if !positive(properties.latent_heat_j_kg)
        || !positive(properties.solid_density_kg_m3)
        || !positive(properties.molar_mass_kg_mol)
        || !nonnegative(properties.saturation_pressure_pa)
        || !nonnegative(properties.saturation_pressure_derivative_pa_k)
        || !properties.latent_heat_derivative_j_kg_k.is_finite()
        || (properties.molar_mass_kg_mol / molar_mass - 1.0).abs() > 64.0 * f64::EPSILON
    {
        return Err(invalid(
            "inadmissible phase data or disagreement with the original gas species molar mass",
        ));
    }
    Ok(properties)
}
fn compatibility(
    vapor_enthalpy: f64,
    solid_enthalpy: f64,
    offset: f64,
    latent: f64,
    temperature: f64,
    tolerance: f64,
) -> Result<(), GasSurfaceError> {
    let residual = (vapor_enthalpy - (solid_enthalpy + offset)) - latent;
    if !residual.is_finite() || residual.abs() > tolerance {
        return Err(GasSurfaceError::CaloricMismatch {
            temperature_k: temperature,
            residual_j_kg: residual,
            tolerance_j_kg: tolerance,
        });
    }
    Ok(())
}

#[derive(Clone, Copy)]
struct GasSpeciesSnapshot {
    molar_mass: f64,
    specific_gas_constant: f64,
    temperature_k: f64,
    partial_pressure_pa: f64,
    vapor_mass_kg: f64,
    incoming_specific_energy_j_kg: Option<f64>,
}

#[derive(Clone, Copy)]
struct SurfaceSnapshot {
    solid: SolidSurfaceState,
    gas: GasSpeciesSnapshot,
    temperature_k: f64,
    outgoing_mass_kg: f64,
    incoming_mass_kg: f64,
    transferred_mass_kg: f64,
    outgoing_specific_energy_j_kg: f64,
    incoming_specific_energy_j_kg: Option<f64>,
}

#[derive(Clone, Copy)]
struct PreparedSurface {
    solid: SolidSurfaceState,
    surface_temperature_k: [f64; 2],
    outgoing_mass_kg: f64,
    incoming_mass_kg: f64,
    transferred_mass_kg: f64,
    outgoing_specific_energy_j_kg: f64,
    incoming_specific_energy_j_kg: Option<f64>,
    energy_into_gas_j: f64,
    solid_energy_change_j: f64,
    energy_residual_j: f64,
    recession_m: f64,
}

fn validate_controls(controls: GasSurfaceSolve, duration_s: f64) -> Result<(), GasSurfaceError> {
    controls.gas.gas.validate()?;
    if !nonnegative(duration_s)
        || !positive(controls.maximum_gross_inventory_fraction)
        || controls.maximum_gross_inventory_fraction >= 1.0
        || !positive(controls.latent_compatibility_tolerance_j_kg)
        || !positive(controls.gas.mass_tolerance_kg)
        || !positive(controls.gas.energy_tolerance_j)
    {
        return Err(invalid(
            "finite duration and positive bounded solve controls required",
        ));
    }
    Ok(())
}

// Cache only evaluated numerical views of the ORIGINAL old gas. A chamber
// queries each participating species once, independently of surface count.
fn species_snapshot(
    gas: &ClosedGasParcel,
    catalog: &ThermochemicalGasModel,
    id: &SpeciesId,
) -> Result<GasSpeciesSnapshot, GasSurfaceError> {
    let species = catalog
        .thermodynamics()
        .components()
        .iter()
        .find(|species| species.species() == id)
        .ok_or(invalid(
            "sublimating species needs its complete original gas source",
        ))?;
    let molar_mass = species.molar_mass().value();
    let specific_gas_constant = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / molar_mass;
    if !positive(molar_mass) || !positive(specific_gas_constant) {
        return Err(invalid(
            "original gas species has an unrepresentable molar mass or gas constant",
        ));
    }
    let temperature_k = gas.state().state().temperature;
    let partial_pressure_pa = gas.species_partial_pressure_pa(id)?;
    let vapor_mass_kg = gas.species_mass_kg(id)?;
    let incoming_specific_energy_j_kg = if vapor_mass_kg == 0.0 {
        None
    } else {
        let energy = catalog.species_specific_enthalpy_j_kg(id, temperature_k)?
            - 0.5 * specific_gas_constant * temperature_k;
        if !energy.is_finite() {
            return Err(invalid("incoming molecular energy is unrepresentable"));
        }
        Some(energy)
    };
    Ok(GasSpeciesSnapshot {
        molar_mass,
        specific_gas_constant,
        temperature_k,
        partial_pressure_pa,
        vapor_mass_kg,
        incoming_specific_energy_j_kg,
    })
}

// No endpoint or gas reconstruction here: all surfaces must first see the
// same old gas and jointly admit their gross inventory budget.
fn snapshot_surface<L, F>(
    law: &L,
    binding: &GasSurfaceBinding,
    solid: SolidSurfaceState,
    gas: GasSpeciesSnapshot,
    catalog: &ThermochemicalGasModel,
    temperature_from_enthalpy: &F,
    duration_s: f64,
    controls: GasSurfaceSolve,
) -> Result<SurfaceSnapshot, GasSurfaceError>
where
    L: SublimationLaw + ?Sized,
    F: Fn(f64) -> Result<f64, SublimationError> + ?Sized,
{
    if !positive(solid.mass_kg)
        || !solid.specific_enthalpy_j_kg.is_finite()
        || !nonnegative(solid.area_m2)
        || !binding.solid_enthalpy_offset_j_kg.is_finite()
        || !binding.accommodation_coefficient.is_finite()
        || !(0.0..=1.0).contains(&binding.accommodation_coefficient)
    {
        return Err(invalid(
            "finite surface state, datum and admissible sticking probability required",
        ));
    }
    let temperature_k = temperature_from_enthalpy(solid.specific_enthalpy_j_kg)?;
    let phase = phase_at(law, temperature_k, gas.molar_mass)?;
    let vapor_enthalpy = catalog.species_specific_enthalpy_j_kg(&binding.species, temperature_k)?;
    compatibility(
        vapor_enthalpy,
        solid.specific_enthalpy_j_kg,
        binding.solid_enthalpy_offset_j_kg,
        phase.latent_heat_j_kg,
        temperature_k,
        controls.latent_compatibility_tolerance_j_kg,
    )?;
    let flux = molecular_mass_flux(
        law,
        temperature_k,
        binding.accommodation_coefficient,
        VaporReservoir::Prescribed {
            partial_pressure_pa: gas.partial_pressure_pa,
            temperature_k: gas.temperature_k,
        },
    )?;
    let exposure = solid.area_m2 * duration_s;
    let outgoing_mass_kg = flux.outgoing_kg_m2_s * exposure;
    let incoming_mass_kg = flux.incoming_kg_m2_s * exposure;
    let transferred_mass_kg = flux.net_kg_m2_s * exposure;
    let outgoing_specific_energy_j_kg =
        vapor_enthalpy - 0.5 * gas.specific_gas_constant * temperature_k;
    if !nonnegative(exposure)
        || !nonnegative(outgoing_mass_kg)
        || !nonnegative(incoming_mass_kg)
        || !transferred_mass_kg.is_finite()
        || !outgoing_specific_energy_j_kg.is_finite()
        || (gas.vapor_mass_kg == 0.0 && incoming_mass_kg != 0.0)
    {
        return Err(invalid(
            "molecular gross streams or energy are unrepresentable",
        ));
    }
    let incoming_specific_energy_j_kg = if incoming_mass_kg == 0.0 {
        None
    } else {
        gas.incoming_specific_energy_j_kg
    };
    Ok(SurfaceSnapshot {
        solid,
        gas,
        temperature_k,
        outgoing_mass_kg,
        incoming_mass_kg,
        transferred_mass_kg,
        outgoing_specific_energy_j_kg,
        incoming_specific_energy_j_kg,
    })
}

// Called only for finite, nonnegative gross transfers. Rounding the admitted
// duration downward makes the returned bound conservative; zero means no
// positive duration was representable under this original inventory budget.
fn inventory_limit(
    duration_s: f64,
    gross_mass_kg: f64,
    inventory_kg: f64,
    maximum_fraction: f64,
) -> Option<f64> {
    if gross_mass_kg == 0.0
        || (inventory_kg > 0.0 && gross_mass_kg / inventory_kg <= maximum_fraction)
    {
        return None;
    }
    let limit = duration_s * ((inventory_kg * maximum_fraction) / gross_mass_kg);
    Some(if limit > 0.0 {
        limit.next_down().max(0.0)
    } else {
        0.0
    })
}

fn prepare_surface<L, F>(
    law: &L,
    binding: &GasSurfaceBinding,
    snapshot: SurfaceSnapshot,
    catalog: &ThermochemicalGasModel,
    temperature_from_enthalpy: &F,
    controls: GasSurfaceSolve,
) -> Result<PreparedSurface, GasSurfaceError>
where
    L: SublimationLaw + ?Sized,
    F: Fn(f64) -> Result<f64, SublimationError> + ?Sized,
{
    let SurfaceSnapshot {
        solid,
        gas,
        temperature_k,
        outgoing_mass_kg,
        incoming_mass_kg,
        transferred_mass_kg,
        outgoing_specific_energy_j_kg,
        incoming_specific_energy_j_kg,
    } = snapshot;
    let difference =
        incoming_specific_energy_j_kg.map_or(0.0, |energy| outgoing_specific_energy_j_kg - energy);
    let energy_into_gas_j =
        transferred_mass_kg * outgoing_specific_energy_j_kg + incoming_mass_kg * difference;
    let next_mass = solid.mass_kg - transferred_mass_kg;
    if !energy_into_gas_j.is_finite()
        || !positive(next_mass)
        || (transferred_mass_kg != 0.0 && next_mass == solid.mass_kg)
    {
        return Err(invalid(
            "molecular solid mass or energy increment is unrepresentable",
        ));
    }
    let offset = binding.solid_enthalpy_offset_j_kg;
    let specific_change = (transferred_mass_kg
        * (solid.specific_enthalpy_j_kg + offset - outgoing_specific_energy_j_kg)
        - incoming_mass_kg * difference)
        / next_mass;
    let next_enthalpy = solid.specific_enthalpy_j_kg + specific_change;
    if !specific_change.is_finite()
        || !next_enthalpy.is_finite()
        || (specific_change != 0.0 && next_enthalpy == solid.specific_enthalpy_j_kg)
    {
        return Err(invalid(
            "molecular solid caloric increment is unrepresentable",
        ));
    }
    let final_temperature = temperature_from_enthalpy(next_enthalpy)?;
    let final_phase = phase_at(law, final_temperature, gas.molar_mass)?;
    let final_vapor_enthalpy =
        catalog.species_specific_enthalpy_j_kg(&binding.species, final_temperature)?;
    compatibility(
        final_vapor_enthalpy,
        next_enthalpy,
        offset,
        final_phase.latent_heat_j_kg,
        final_temperature,
        controls.latent_compatibility_tolerance_j_kg,
    )?;
    let solid_energy_change_j = next_mass * (next_enthalpy - solid.specific_enthalpy_j_kg)
        + (next_mass - solid.mass_kg) * (solid.specific_enthalpy_j_kg + offset);
    let energy_residual_j = solid_energy_change_j + energy_into_gas_j;
    if !solid_energy_change_j.is_finite()
        || !energy_residual_j.is_finite()
        || energy_residual_j.abs() > controls.gas.energy_tolerance_j
    {
        return Err(GasSurfaceError::EnergyBalance {
            residual_j: energy_residual_j,
            tolerance_j: controls.gas.energy_tolerance_j,
        });
    }
    let recession_m = if solid.area_m2 == 0.0 {
        0.0
    } else {
        transferred_mass_kg / final_phase.solid_density_kg_m3 / solid.area_m2
    };
    if !recession_m.is_finite() {
        return Err(invalid("molecular surface recession is unrepresentable"));
    }
    Ok(PreparedSurface {
        solid: SolidSurfaceState {
            mass_kg: next_mass,
            specific_enthalpy_j_kg: next_enthalpy,
            area_m2: solid.area_m2,
        },
        surface_temperature_k: [temperature_k, final_temperature],
        outgoing_mass_kg,
        incoming_mass_kg,
        transferred_mass_kg,
        outgoing_specific_energy_j_kg,
        incoming_specific_energy_j_kg,
        energy_into_gas_j,
        solid_energy_change_j,
        energy_residual_j,
        recession_m,
    })
}

/// One explicit conservative molecular interval against an evolving original gas.
///
/// The callback inverts the original solid specific-enthalpy chart. The binding
/// must align that chart with the original species h, with the declared latent
/// difference at both endpoints. No constant-cp gas or second energy inversion
/// is used. The original parcel reconstructs all accepted gas properties once.
///
/// All gross rates use the initial state. A step that exceeds either gross
/// inventory allowance returns a conservative duration bound instead of clipping
/// a stream. Exact extinction/vacuum and unsupported final states refuse.
/// Both endpoints are returned together; no accepted input is ever modified.
#[allow(clippy::too_many_arguments)]
pub fn advance_surface_with_gas<'catalog, L: SublimationLaw + ?Sized>(
    law: &L,
    binding: &GasSurfaceBinding,
    solid: SolidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    temperature_from_enthalpy: impl Fn(f64) -> Result<f64, SublimationError>,
    duration_s: f64,
    controls: GasSurfaceSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasSurfaceCandidate<'catalog>, GasSurfaceError> {
    poll(&mut cancelled)?;
    validate_controls(controls, duration_s)?;
    let gas_snapshot = species_snapshot(gas, catalog, &binding.species)?;
    let snapshot = snapshot_surface(
        law,
        binding,
        solid,
        gas_snapshot,
        catalog,
        &temperature_from_enthalpy,
        duration_s,
        controls,
    )?;
    let limit = [
        inventory_limit(
            duration_s,
            snapshot.outgoing_mass_kg,
            solid.mass_kg,
            controls.maximum_gross_inventory_fraction,
        ),
        inventory_limit(
            duration_s,
            snapshot.incoming_mass_kg,
            gas_snapshot.vapor_mass_kg,
            controls.maximum_gross_inventory_fraction,
        ),
    ]
    .into_iter()
    .flatten()
    .reduce(f64::min);
    if let Some(maximum_duration_s) = limit {
        return Err(GasSurfaceError::InventoryLimit { maximum_duration_s });
    }
    poll(&mut cancelled)?;
    let prepared = prepare_surface(
        law,
        binding,
        snapshot,
        catalog,
        &temperature_from_enthalpy,
        controls,
    )?;
    poll(&mut cancelled)?;
    let row = SpeciesBoundaryTransfer {
        species: binding.species.clone(),
        mass_into_gas_kg: prepared.transferred_mass_kg,
        energy_into_gas_j: prepared.energy_into_gas_j,
    };
    let endpoint = gas.exchange_species_boundary(catalog, &[row], controls.gas, &mut cancelled)?;
    let gas_change = endpoint.step().stored_energy_change_j;
    let residual = prepared.solid_energy_change_j + gas_change;
    if !residual.is_finite() || residual.abs() > controls.gas.energy_tolerance_j {
        return Err(GasSurfaceError::EnergyBalance {
            residual_j: residual,
            tolerance_j: controls.gas.energy_tolerance_j,
        });
    }
    let mass_residual =
        (prepared.solid.mass_kg - solid.mass_kg) + (endpoint.gas().mass_kg() - gas.mass_kg());
    if !mass_residual.is_finite() || mass_residual.abs() > controls.gas.mass_tolerance_kg {
        return Err(invalid("molecular endpoints failed total mass accounting"));
    }
    let final_pressure = endpoint
        .gas()
        .species_partial_pressure_pa(&binding.species)?;
    let final_gas_temperature = endpoint.gas().state().state().temperature;
    poll(&mut cancelled)?;
    Ok(GasSurfaceCandidate {
        solid: prepared.solid,
        gas: endpoint,
        step: GasSurfaceStep {
            outgoing_mass_kg: prepared.outgoing_mass_kg,
            incoming_mass_kg: prepared.incoming_mass_kg,
            transferred_mass_kg: prepared.transferred_mass_kg,
            surface_temperature_k: prepared.surface_temperature_k,
            gas_temperature_k: [gas_snapshot.temperature_k, final_gas_temperature],
            vapor_partial_pressure_pa: [gas_snapshot.partial_pressure_pa, final_pressure],
            outgoing_specific_energy_j_kg: prepared.outgoing_specific_energy_j_kg,
            incoming_specific_energy_j_kg: prepared.incoming_specific_energy_j_kg,
            energy_into_gas_j: prepared.energy_into_gas_j,
            solid_energy_change_j: prepared.solid_energy_change_j,
            gas_energy_change_j: gas_change,
            energy_residual_j: residual,
            mass_residual_kg: mass_residual,
            recession_m: prepared.recession_m,
        },
    })
}

#[cfg(test)]
mod tests;
