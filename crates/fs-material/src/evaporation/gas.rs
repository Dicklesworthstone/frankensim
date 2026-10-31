//! Continuum liquid evaporation/condensation into the original finite gas.
//!
//! The shared Stefan-film law determines signed species mass. The original
//! liquid enthalpy chart and NASA-9 gas owner determine the two energy changes;
//! no constant-cp gas, humidity clipping, or second equation of state is used.
//! A fixed, explicit caloric-datum offset must reconcile the phase and gas
//! sources. An admitted source mismatch is retained in the energy transfer,
//! rather than silently losing mass times that mismatch.
//!
//! This is a first-order, frozen-rate interval. The film uses the initial gas
//! temperature, pressure and composition and a caller-supplied mass-transfer
//! coefficient. Incoming vapor carries its upstream gas enthalpy; outgoing
//! vapor carries its surface enthalpy. The optional [`heat`] path first solves
//! reciprocal sensible film heat on the same original calorics. A step that
//! overshoots equilibrium, consumes too much inventory or
//! changes temperature too far refuses with a smaller duration to retry.
//!
//! Gas volume and thermal reference geometry are fixed. The conserved thermal
//! account is liquid enthalpy plus gas internal energy; this is not a moving
//! liquid boundary, a rigid two-phase vessel, a boiling or a spray model.

/// Simultaneous retained liquid films sharing one finite gas inventory.
pub mod chamber;
/// Implicit reciprocal sensible film heat, separately from vapor-carried energy.
pub mod heat;
/// Original heat transport and finite humidity feedback on one physical clock.
pub mod timeline;
/// Fully-wetting wall films with dry onset and complete finite-inventory dry-out.
pub mod wall_film;
pub use heat::{
    GasEvaporationHeatCandidate, GasFilmHeatSolve, GasFilmHeatStep,
    advance_surface_with_gas_and_heat, advance_surface_with_gas_and_heat_mass_selection,
};

use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{
        GasEnergyError,
        parcel::{
            ClosedGasParcel,
            transfer::species::boundary::{
                GasSpeciesBoundaryCandidate, GasSpeciesBoundarySolve,
                SpeciesBoundaryRelativeTransfer,
            },
        },
    },
};
use fs_thermal::evaporation::{
    EvaporationError, EvaporationLaw, GasFilm, LiquidSurfaceState, properties, surface_flux,
};
use fs_thermochem::SpeciesId;

/// Neither volatile identity nor its energy reference follows from a name.
#[derive(Clone, Debug, PartialEq)]
pub struct GasLiquidBinding {
    pub species: SpeciesId,
    /// Added to the original liquid h to express storage in the gas datum.
    pub liquid_enthalpy_offset_j_kg: f64,
}

/// Explicit physical-step and original gas reconstruction budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasEvaporationSolve {
    pub gas: GasSpeciesBoundarySolve,
    /// Maximum withdrawal divided by initial donor mass, strictly in (0,1).
    pub maximum_inventory_fraction: f64,
    /// Maximum absolute departure of either surface or gas temperature [K].
    pub maximum_temperature_increment_k: f64,
    /// Admitted |h_v - (h_l + offset) - L| at both surface endpoints [J/kg].
    pub latent_compatibility_tolerance_j_kg: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum GasEvaporationError {
    Invalid(&'static str),
    Surface(EvaporationError),
    Gas(GasEnergyError),
    /// Retry from the same original states; no state was published.
    StepSizeLimit {
        maximum_duration_s: f64,
        reason: &'static str,
    },
    CaloricMismatch {
        temperature_k: f64,
        residual_j_kg: f64,
        tolerance_j_kg: f64,
    },
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
    /// The admitted original gas/liquid branches do not contain the heat endpoint.
    HeatOutsideBracket {
        residual_temperature_k: [f64; 2],
    },
    /// The implicit sensible heat equation could not meet its temperature budget.
    HeatNotConverged {
        iterations: usize,
        residual_temperature_k: f64,
    },
    /// The admitted source domains do not bracket the joint heat-energy root.
    HeatEnergyOutsideBracket {
        residual_energy_j: [f64; 2],
    },
    /// Bounded scalar gas-temperature work was exhausted without a heat endpoint.
    HeatEnergyNotConverged {
        iterations: usize,
        residual_energy_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for GasEvaporationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-gas liquid evaporation: {self:?}")
    }
}
impl std::error::Error for GasEvaporationError {}
impl From<EvaporationError> for GasEvaporationError {
    fn from(value: EvaporationError) -> Self {
        Self::Surface(value)
    }
}
impl From<GasEnergyError> for GasEvaporationError {
    fn from(value: GasEnergyError) -> Self {
        Self::Gas(value)
    }
}

/// All accounts refer to the same admitted signed species transfer.
#[derive(Clone, Debug, PartialEq)]
pub struct GasEvaporationStep {
    pub duration_s: f64,
    /// Original film flux integrated over area and duration, before the native
    /// inventory representation selects its nearby admitted increment.
    pub requested_transferred_mass_kg: f64,
    /// Positive liquid-to-gas evaporation; negative condensation.
    pub transferred_mass_kg: f64,
    pub surface_temperature_k: [f64; 2],
    pub gas_temperature_k: [f64; 2],
    pub vapor_partial_pressure_pa: [f64; 2],
    /// Phase-law latent portion only, already included in carried vapor energy.
    pub latent_energy_j: f64,
    /// Complete signed carried vapor energy, including the original gas datum.
    pub energy_into_gas_j: f64,
    /// Actual caloric increment plus intended transferred-mass reference energy.
    /// Represented inventory rounding is separately bounded by the mass budget.
    pub liquid_enthalpy_change_j: f64,
    /// Original gas owner's source-relative storage change, with the same
    /// separately budgeted mass-reference rounding convention.
    pub gas_internal_energy_change_j: f64,
    /// Formation references cancel before evaluating this actual thermal defect.
    pub energy_residual_j: f64,
    /// Actual liquid mass change plus actual volatile gas mass change.
    pub mass_residual_kg: f64,
    /// Positive thickness loss; geometry remains caller-owned.
    pub recession_m: f64,
}

/// Proposed original gas and liquid states, accepted or discarded together.
#[derive(Clone, Debug, PartialEq)]
pub struct GasEvaporationCandidate<'catalog> {
    liquid: LiquidSurfaceState,
    gas: GasSpeciesBoundaryCandidate<'catalog>,
    step: GasEvaporationStep,
}
impl<'catalog> GasEvaporationCandidate<'catalog> {
    pub const fn liquid(&self) -> LiquidSurfaceState {
        self.liquid
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        self.gas.gas()
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.gas.catalog()
    }
    pub const fn step(&self) -> &GasEvaporationStep {
        &self.step
    }
    pub const fn gas_exchange(&self) -> &GasSpeciesBoundaryCandidate<'catalog> {
        &self.gas
    }
    pub fn into_parts(
        self,
    ) -> (
        LiquidSurfaceState,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        GasEvaporationStep,
    ) {
        let (gas, catalog, _) = self.gas.into_parts();
        (self.liquid, gas, catalog, self.step)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasEvaporationError> {
    if cancelled() {
        Err(GasEvaporationError::Cancelled)
    } else {
        Ok(())
    }
}
fn smaller_step(duration_s: f64, reason: &'static str) -> GasEvaporationError {
    GasEvaporationError::StepSizeLimit {
        maximum_duration_s: (0.5 * duration_s).max(0.0),
        reason,
    }
}
fn compatibility(
    vapor: f64,
    liquid: f64,
    offset: f64,
    latent: f64,
    temperature: f64,
    tolerance: f64,
) -> Result<f64, GasEvaporationError> {
    let residual = (vapor - (liquid + offset)) - latent;
    if !residual.is_finite() || residual.abs() > tolerance {
        return Err(GasEvaporationError::CaloricMismatch {
            temperature_k: temperature,
            residual_j_kg: residual,
            tolerance_j_kg: tolerance,
        });
    }
    Ok(residual)
}

pub(super) fn validate_controls(
    controls: GasEvaporationSolve,
    duration_s: f64,
    mass_transfer_coefficient_m_s: f64,
) -> Result<(), GasEvaporationError> {
    controls.gas.gas.validate()?;
    if !nonnegative(duration_s)
        || !nonnegative(mass_transfer_coefficient_m_s)
        || !positive(controls.maximum_inventory_fraction)
        || controls.maximum_inventory_fraction >= 1.0
        || !positive(controls.maximum_temperature_increment_k)
        || !nonnegative(controls.latent_compatibility_tolerance_j_kg)
        || !positive(controls.gas.mass_tolerance_kg)
        || !positive(controls.gas.energy_tolerance_j)
    {
        return Err(GasEvaporationError::Invalid(
            "finite film/duration and positive bounded solve controls required",
        ));
    }
    Ok(())
}

/// Transfer one bounded interval through the shared continuum film law.
///
/// The callback inverts the original, admitted pure-liquid enthalpy branch.
/// The retained catalog permits a currently absent vapor to enter with its
/// original species, transport and caloric sources. Carrier species amounts,
/// auxiliary query coordinates and gas volume survive reconstruction.
///
/// Refine physical time for the frozen-film/rate approximation. The temperature
/// and equilibrium-crossing checks are admission limits, not error estimates.
/// Every failure, including cancellation, leaves both caller-owned inputs intact.
#[allow(clippy::too_many_arguments)]
pub fn advance_surface_with_gas<'catalog, L: EvaporationLaw + ?Sized>(
    law: &L,
    binding: &GasLiquidBinding,
    liquid: LiquidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    temperature_from_enthalpy: impl Fn(f64) -> Result<f64, EvaporationError>,
    mass_transfer_coefficient_m_s: f64,
    duration_s: f64,
    controls: GasEvaporationSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasEvaporationCandidate<'catalog>, GasEvaporationError> {
    advance_surface_with_gas_mass_selection(
        law,
        binding,
        liquid,
        gas,
        catalog,
        temperature_from_enthalpy,
        mass_transfer_coefficient_m_s,
        duration_s,
        controls,
        Ok,
        cancelled,
    )
}

/// Original phase/gas source admission shared by free and wall-bound films.
struct GasLiquidSource {
    phase: fs_thermal::evaporation::EvaporationProperties,
    molecular_mass: f64,
    alignment: f64,
}

fn admit_liquid_source<L: EvaporationLaw + ?Sized>(
    law: &L,
    binding: &GasLiquidBinding,
    catalog: &ThermochemicalGasModel,
    surface_t: f64,
    liquid_h: f64,
    controls: GasEvaporationSolve,
) -> Result<GasLiquidSource, GasEvaporationError> {
    use GasEvaporationError as Error;
    if !liquid_h.is_finite() || !binding.liquid_enthalpy_offset_j_kg.is_finite() {
        return Err(Error::Invalid(
            "finite original liquid caloric datum required",
        ));
    }
    let phase = properties(law, surface_t)?;
    let original_species = catalog
        .thermodynamics()
        .components()
        .iter()
        .find(|source| source.species() == &binding.species)
        .ok_or(Error::Invalid(
            "evaporating species requires its complete original gas source",
        ))?;
    let molecular_mass = original_species.molar_mass().value();
    if !positive(molecular_mass)
        || (phase.molar_mass_kg_mol / molecular_mass - 1.0).abs() > 64.0 * f64::EPSILON
    {
        return Err(Error::Invalid(
            "phase and original gas species molar masses disagree",
        ));
    }
    let vapor_h = catalog.species_specific_enthalpy_j_kg(&binding.species, surface_t)?;
    let alignment = compatibility(
        vapor_h,
        liquid_h,
        binding.liquid_enthalpy_offset_j_kg,
        phase.latent_heat_j_kg,
        surface_t,
        controls.latent_compatibility_tolerance_j_kg,
    )?;
    Ok(GasLiquidSource {
        phase,
        molecular_mass,
        alignment,
    })
}

/// Prepared entirely against one original gas state. The chamber and legacy
/// single-surface paths share this rate, liquid-caloric and source-alignment owner.
struct PreparedGasLiquidSurface {
    original: LiquidSurfaceState,
    liquid: LiquidSurfaceState,
    requested_dm: f64,
    dm: f64,
    surface_temperature_k: [f64; 2],
    gas_temperature_k: f64,
    vapor_partial_pressure_pa: f64,
    old_vapor_mass_kg: f64,
    relative_energy_into_gas_j: f64,
    liquid_energy_residual_j: f64,
    liquid_mass_residual_kg: f64,
    liquid_enthalpy_change_j: f64,
    latent_energy_j: f64,
    recession_m: f64,
}

#[allow(clippy::too_many_arguments)]
fn prepare_surface<L: EvaporationLaw + ?Sized>(
    law: &L,
    binding: &GasLiquidBinding,
    liquid: LiquidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &ThermochemicalGasModel,
    temperature_from_enthalpy: impl Fn(f64) -> Result<f64, EvaporationError>,
    mass_transfer_coefficient_m_s: f64,
    duration_s: f64,
    controls: GasEvaporationSolve,
    select_mass: impl FnOnce(f64) -> Result<f64, GasEvaporationError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PreparedGasLiquidSurface, GasEvaporationError> {
    use GasEvaporationError as Error;
    poll(&mut cancelled)?;
    validate_controls(controls, duration_s, mass_transfer_coefficient_m_s)?;
    if !positive(liquid.mass_kg)
        || !nonnegative(liquid.area_m2)
        || !liquid.specific_enthalpy_j_kg.is_finite()
        || !binding.liquid_enthalpy_offset_j_kg.is_finite()
    {
        return Err(Error::Invalid(
            "finite liquid/film data and positive bounded solve controls required",
        ));
    }
    let surface_t = temperature_from_enthalpy(liquid.specific_enthalpy_j_kg)?;
    let gas_t = gas.state().state().temperature;
    let partial_pressure = gas.species_partial_pressure_pa(&binding.species)?;
    let film = GasFilm {
        total_pressure_pa: gas.state().state().pressure,
        vapor_partial_pressure_pa: partial_pressure,
        temperature_k: gas_t,
        mass_transfer_coefficient_m_s,
    };
    let flux = surface_flux(law, surface_t, film)?;
    let source = admit_liquid_source(
        law,
        binding,
        catalog,
        surface_t,
        liquid.specific_enthalpy_j_kg,
        controls,
    )?;
    let phase = source.phase;
    let molecular_mass = source.molecular_mass;
    let alignment = source.alignment;
    let exposure = liquid.area_m2 * duration_s;
    let requested_dm = exposure * flux.mass_flux_kg_m2_s;
    if !nonnegative(exposure)
        || !requested_dm.is_finite()
        || (duration_s > 0.0 && liquid.area_m2 > 0.0 && exposure == 0.0)
        || (exposure > 0.0 && flux.mass_flux_kg_m2_s != 0.0 && requested_dm == 0.0)
    {
        return Err(Error::Invalid(
            "exposure or nonzero transferred mass is unrepresentable",
        ));
    }
    let dm = select_mass(requested_dm)?;
    if !dm.is_finite()
        || (requested_dm == 0.0 && dm != 0.0)
        || (requested_dm != 0.0
            && (dm == 0.0 || dm.is_sign_negative() != requested_dm.is_sign_negative()))
        || (dm - requested_dm).abs() > controls.gas.mass_tolerance_kg
    {
        return Err(Error::Invalid(
            "native phase mass selection changes or loses the physical transfer",
        ));
    }
    let old_vapor_mass = gas.species_mass_kg(&binding.species)?;
    let donor_mass = if dm < 0.0 {
        old_vapor_mass
    } else {
        liquid.mass_kg
    };
    let withdrawal = dm.abs().max(requested_dm.abs());
    if dm != 0.0
        && (donor_mass <= 0.0 || withdrawal / donor_mass > controls.maximum_inventory_fraction)
    {
        let duration =
            duration_s * ((controls.maximum_inventory_fraction * donor_mass) / withdrawal);
        return Err(Error::StepSizeLimit {
            maximum_duration_s: duration.next_down().max(0.0),
            reason: "initial liquid or vapor inventory withdrawal limit",
        });
    }
    let next_mass = liquid.mass_kg - dm;
    let liquid_mass_residual = (next_mass - liquid.mass_kg) + dm;
    if !positive(next_mass)
        || (dm != 0.0 && next_mass == liquid.mass_kg)
        || !liquid_mass_residual.is_finite()
        || liquid_mass_residual.abs() > controls.gas.mass_tolerance_kg
    {
        return Err(Error::Invalid("liquid inventory change is unrepresentable"));
    }
    // Use the original source increment: formation constants cannot erase
    // sensible heat when either source's absolute h has a large datum.
    let vapor_h_difference = if dm == 0.0 {
        0.0
    } else {
        catalog.species_specific_enthalpy_change_j_kg(&binding.species, gas_t, surface_t)?
    };
    let carried_gap =
        phase.latent_heat_j_kg + alignment - if dm < 0.0 { vapor_h_difference } else { 0.0 };
    let liquid_thermal_target = -dm * carried_gap;
    let delta_h = liquid_thermal_target / next_mass;
    let next_h = liquid.specific_enthalpy_j_kg + delta_h;
    if !liquid_thermal_target.is_finite()
        || !delta_h.is_finite()
        || !next_h.is_finite()
        || (liquid_thermal_target != 0.0 && delta_h == 0.0)
        || (delta_h != 0.0 && next_h == liquid.specific_enthalpy_j_kg)
    {
        return Err(Error::Invalid(
            "liquid caloric increment is unrepresentable",
        ));
    }
    let final_surface_t = temperature_from_enthalpy(next_h)?;
    let final_phase = properties(law, final_surface_t)?;
    compatibility(
        catalog.species_specific_enthalpy_j_kg(&binding.species, final_surface_t)?,
        next_h,
        binding.liquid_enthalpy_offset_j_kg,
        final_phase.latent_heat_j_kg,
        final_surface_t,
        controls.latent_compatibility_tolerance_j_kg,
    )?;
    if (final_phase.molar_mass_kg_mol / molecular_mass - 1.0).abs() > 64.0 * f64::EPSILON {
        return Err(Error::Invalid(
            "volatile molar mass must remain the same original species",
        ));
    }
    if (final_surface_t - surface_t).abs() > controls.maximum_temperature_increment_k {
        return Err(smaller_step(duration_s, "surface temperature departure"));
    }
    let liquid_thermal_change = next_mass * (next_h - liquid.specific_enthalpy_j_kg);
    let liquid_residual = liquid_thermal_change - liquid_thermal_target;
    let relative_energy = if dm > 0.0 {
        dm * vapor_h_difference
    } else {
        0.0
    };
    let liquid_enthalpy_change = liquid_thermal_change
        - dm * (liquid.specific_enthalpy_j_kg + binding.liquid_enthalpy_offset_j_kg);
    let recession = if liquid.area_m2 == 0.0 {
        0.0
    } else {
        dm / final_phase.liquid_density_kg_m3 / liquid.area_m2
    };
    let latent_energy = dm * phase.latent_heat_j_kg;
    if ![
        liquid_residual,
        liquid_enthalpy_change,
        recession,
        latent_energy,
        relative_energy,
    ]
    .into_iter()
    .all(f64::is_finite)
    {
        return Err(Error::Invalid(
            "prepared liquid/gas accounts are unrepresentable",
        ));
    }
    Ok(PreparedGasLiquidSurface {
        original: liquid,
        liquid: LiquidSurfaceState {
            mass_kg: next_mass,
            specific_enthalpy_j_kg: next_h,
            area_m2: liquid.area_m2,
        },
        requested_dm,
        dm,
        surface_temperature_k: [surface_t, final_surface_t],
        gas_temperature_k: gas_t,
        vapor_partial_pressure_pa: partial_pressure,
        old_vapor_mass_kg: old_vapor_mass,
        relative_energy_into_gas_j: relative_energy,
        liquid_energy_residual_j: liquid_residual,
        liquid_mass_residual_kg: liquid_mass_residual,
        liquid_enthalpy_change_j: liquid_enthalpy_change,
        latent_energy_j: latent_energy,
        recession_m: recession,
    })
}

/// Every surface checks the same accepted gas endpoint; local preparation is
/// insufficient to certify a joint finite-atmosphere step.
#[allow(clippy::too_many_arguments)]
fn check_surface_endpoint<L: EvaporationLaw + ?Sized>(
    prepared: &PreparedGasLiquidSurface,
    law: &L,
    binding: &GasLiquidBinding,
    gas: &ClosedGasParcel,
    mass_transfer_coefficient_m_s: f64,
    duration_s: f64,
    controls: GasEvaporationSolve,
) -> Result<f64, GasEvaporationError> {
    check_interface_endpoint(
        law,
        binding,
        prepared.surface_temperature_k[1],
        prepared.gas_temperature_k,
        prepared.dm,
        gas,
        mass_transfer_coefficient_m_s,
        duration_s,
        controls,
    )
}

#[allow(clippy::too_many_arguments)]
fn check_interface_endpoint<L: EvaporationLaw + ?Sized>(
    law: &L,
    binding: &GasLiquidBinding,
    surface_temperature_k: f64,
    initial_gas_temperature_k: f64,
    transferred_mass_kg: f64,
    gas: &ClosedGasParcel,
    mass_transfer_coefficient_m_s: f64,
    duration_s: f64,
    controls: GasEvaporationSolve,
) -> Result<f64, GasEvaporationError> {
    let final_gas_t = gas.state().state().temperature;
    let final_pressure = gas.species_partial_pressure_pa(&binding.species)?;
    if (final_gas_t - initial_gas_temperature_k).abs() > controls.maximum_temperature_increment_k {
        return Err(smaller_step(duration_s, "gas temperature departure"));
    }
    let final_flux = surface_flux(
        law,
        surface_temperature_k,
        GasFilm {
            total_pressure_pa: gas.state().state().pressure,
            vapor_partial_pressure_pa: final_pressure,
            temperature_k: final_gas_t,
            mass_transfer_coefficient_m_s,
        },
    )?;
    if (transferred_mass_kg > 0.0 && final_flux.mass_flux_kg_m2_s < 0.0)
        || (transferred_mass_kg < 0.0 && final_flux.mass_flux_kg_m2_s > 0.0)
    {
        return Err(smaller_step(
            duration_s,
            "finite-atmosphere equilibrium overshoot",
        ));
    }
    Ok(final_pressure)
}

/// The same original film and caloric solve with a represented mass selector.
///
/// Native volume/mass owners may select the nearest increment they can store.
/// The callback receives signed liquid-to-gas mass [kg] and runs BEFORE either
/// caloric update. It cannot reverse a transfer, swallow a nonzero transfer,
/// create one from zero or change it beyond the original absolute mass budget.
/// Inventory, temperature and equilibrium limits apply to the selected amount.
/// This does not permit clipping a physically inadmissible step to inventory.
#[allow(clippy::too_many_arguments)]
pub fn advance_surface_with_gas_mass_selection<'catalog, L: EvaporationLaw + ?Sized>(
    law: &L,
    binding: &GasLiquidBinding,
    liquid: LiquidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    temperature_from_enthalpy: impl Fn(f64) -> Result<f64, EvaporationError>,
    mass_transfer_coefficient_m_s: f64,
    duration_s: f64,
    controls: GasEvaporationSolve,
    select_mass: impl FnOnce(f64) -> Result<f64, GasEvaporationError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasEvaporationCandidate<'catalog>, GasEvaporationError> {
    let prepared = prepare_surface(
        law,
        binding,
        liquid,
        gas,
        catalog,
        temperature_from_enthalpy,
        mass_transfer_coefficient_m_s,
        duration_s,
        controls,
        select_mass,
        &mut cancelled,
    )?;
    poll(&mut cancelled)?;
    let endpoint = gas.exchange_species_boundary_relative(
        catalog,
        &[SpeciesBoundaryRelativeTransfer {
            species: binding.species.clone(),
            mass_into_gas_kg: prepared.dm,
            relative_energy_into_gas_j: prepared.relative_energy_into_gas_j,
        }],
        controls.gas,
        &mut cancelled,
    )?;
    let final_pressure = check_surface_endpoint(
        &prepared,
        law,
        binding,
        endpoint.gas(),
        mass_transfer_coefficient_m_s,
        duration_s,
        controls,
    )?;
    let residual = prepared.liquid_energy_residual_j + endpoint.step().energy_residual_j;
    let mass_residual = (prepared.liquid.mass_kg - liquid.mass_kg)
        + (endpoint.gas().species_mass_kg(&binding.species)? - prepared.old_vapor_mass_kg);
    if !residual.is_finite()
        || !mass_residual.is_finite()
        || mass_residual.abs() > controls.gas.mass_tolerance_kg
    {
        return Err(GasEvaporationError::Invalid(
            "accepted liquid/gas accounts are unrepresentable or fail mass balance",
        ));
    }
    if residual.abs() > controls.gas.energy_tolerance_j {
        return Err(GasEvaporationError::EnergyBalance {
            residual_j: residual,
            tolerance_j: controls.gas.energy_tolerance_j,
        });
    }
    poll(&mut cancelled)?;
    let step = GasEvaporationStep {
        duration_s,
        requested_transferred_mass_kg: prepared.requested_dm,
        transferred_mass_kg: prepared.dm,
        surface_temperature_k: prepared.surface_temperature_k,
        gas_temperature_k: [
            prepared.gas_temperature_k,
            endpoint.gas().state().state().temperature,
        ],
        vapor_partial_pressure_pa: [prepared.vapor_partial_pressure_pa, final_pressure],
        latent_energy_j: prepared.latent_energy_j,
        energy_into_gas_j: endpoint.step().energy_into_gas_j,
        liquid_enthalpy_change_j: prepared.liquid_enthalpy_change_j,
        gas_internal_energy_change_j: endpoint.step().stored_energy_change_j,
        energy_residual_j: residual,
        mass_residual_kg: mass_residual,
        recession_m: prepared.recession_m,
    };
    Ok(GasEvaporationCandidate {
        liquid: prepared.liquid,
        gas: endpoint,
        step,
    })
}

#[cfg(test)]
mod tests;
