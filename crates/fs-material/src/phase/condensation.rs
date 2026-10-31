//! Finite moist-gas saturation adjustment with a retained liquid inventory.
//!
//! The original NASA-9 species own gas caloric storage, EOS and transport;
//! IAPWS SR1-86(1992) owns saturation and liquid caloric storage. A caller-chosen
//! reference temperature aligns the IAPWS vapor enthalpy to the explicitly
//! identified NASA water species ONCE. Subtracting unaligned NASA formation
//! enthalpy and IAPWS liquid enthalpy would invent a large phase-change heat.
//!
//! At prescribed T and fixed GAS volume, water partial pressure is limited by
//! p_sat(T). Excess vapor becomes liquid, and existing liquid can re-evaporate
//! only up to its finite inventory. Every noncondensing species is retained.
//! Gas and liquid share one temperature; a thermal owner must solve their
//! combined internal-energy change against supplied heat. Latent energy is
//! already included and must never be added as a second source.
//!
//! This is a near-atmospheric ideal-mixture approximation over 273.16--323.15 K
//! and 20--200 kPa total gas pressure. Liquid properties are the saturated-liquid
//! approximation; liquid volume does not displace the prescribed gas volume.
//! No fugacity/enhancement/Poynting/Henry correction, ice, nucleation rate,
//! interface transport, liquid motion or measured error bound is claimed.
//! See IAPWS SR1-86(1992), <https://iapws.org/relguide/Supp-sat.html>, and the
//! distinct humid-air equilibrium treatment in IAPWS G8-10/G11-15.

use crate::MaterialError;
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{GasEnergyError, IsochoricGasParcel, parcel::ClosedGasParcel},
};
use fs_qty::Temperature;
use fs_thermochem::water::{SaturatedWaterState, WaterSaturationError};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};
use std::collections::BTreeMap;
use std::sync::Arc;

/// Lower temperature admitted without a separate ice owner [K].
pub const MINIMUM_TEMPERATURE_K: f64 = 273.16;
/// Upper temperature of this dilute moist-gas approximation [K].
pub const MAXIMUM_TEMPERATURE_K: f64 = 323.15;
/// Lower total gas pressure admitted by this near-atmospheric model [Pa].
pub const MINIMUM_PRESSURE_PA: f64 = 20_000.0;
/// Upper total gas pressure admitted without a nonideal humid-air model [Pa].
pub const MAXIMUM_PRESSURE_PA: f64 = 200_000.0;

/// A failed proposal never changes the original gas or liquid inventories.
#[derive(Clone, Debug, PartialEq)]
pub enum WaterCondensationError {
    /// Invalid input, incompatible sources or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Original material/source refusal.
    Material(MaterialError),
    /// Original conservative gas/species owner refused the endpoint.
    Gas(GasEnergyError),
    /// Original saturation source refused the temperature.
    Water(WaterSaturationError),
    /// Caller cancellation observed before returning a complete proposal.
    Cancelled,
}
impl core::fmt::Display for WaterCondensationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite water condensation: {self:?}")
    }
}
impl std::error::Error for WaterCondensationError {}
impl From<MaterialError> for WaterCondensationError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<GasEnergyError> for WaterCondensationError {
    fn from(error: GasEnergyError) -> Self {
        Self::Gas(error)
    }
}
impl From<WaterSaturationError> for WaterCondensationError {
    fn from(error: WaterSaturationError) -> Self {
        Self::Water(error)
    }
}
fn bad(message: &'static str) -> WaterCondensationError {
    WaterCondensationError::Invalid(message)
}
fn finite(value: f64) -> Result<f64, WaterCondensationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(bad("phase arithmetic is unrepresentable"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), WaterCondensationError> {
    if cancelled() {
        Err(WaterCondensationError::Cancelled)
    } else {
        Ok(())
    }
}
fn saturation(temperature_k: f64) -> Result<SaturatedWaterState, WaterCondensationError> {
    if !(MINIMUM_TEMPERATURE_K..=MAXIMUM_TEMPERATURE_K).contains(&temperature_k) {
        return Err(bad("liquid-water condensation requires 273.16..=323.15 K"));
    }
    Ok(SaturatedWaterState::at_temperature(Temperature::new(
        temperature_k,
    ))?)
}
fn admit_pressure(pressure_pa: f64) -> Result<(), WaterCondensationError> {
    if !(MINIMUM_PRESSURE_PA..=MAXIMUM_PRESSURE_PA).contains(&pressure_pa) {
        return Err(bad(
            "ideal moist-gas condensation requires 20..=200 kPa total pressure",
        ));
    }
    Ok(())
}

/// Immutable source catalog and one explicit, constant gas/liquid caloric datum.
/// The complete catalog survives depletion of water from the active gas.
#[derive(Clone, Debug, PartialEq)]
pub struct WaterCondensationModel {
    source: Arc<ThermochemicalGasModel>,
    water: SpeciesId,
    reference_temperature_k: f64,
    reference_vapor_enthalpy_j_kg: f64,
    liquid_energy_offset_j_kg: f64,
    water_gas_constant_j_kg_k: f64,
}
impl WaterCondensationModel {
    /// Bind the original sources; water identity is explicit, never name-inferred.
    /// At the chosen Tref, the original IAPWS vaporization enthalpy is retained.
    /// Away from Tref the NASA vapor and IAPWS liquid heat capacities determine
    /// the gap. Their distinct model forms are not silently forced to coincide.
    pub fn new(
        source: ThermochemicalGasModel,
        water: SpeciesId,
        reference_temperature_k: f64,
    ) -> Result<Self, WaterCondensationError> {
        let reference = saturation(reference_temperature_k)?;
        if source.water_species() != Some(&water) {
            return Err(bad(
                "condensation needs the source's explicitly identified water species",
            ));
        }
        let species = source
            .thermodynamics()
            .components()
            .iter()
            .find(|species| species.species() == &water)
            .ok_or(bad("water source is absent"))?;
        // IAPWS ordinary water and common NASA/JANAF tables round the molar
        // mass slightly differently. Permit 10 ppm, never a different molecule
        // relabelled as water or an unqualified isotopologue substitution.
        if (species.molar_mass().value() / 0.018_015_268 - 1.0).abs() > 1.0e-5 {
            return Err(bad(
                "IAPWS ordinary water requires a compatible H2O molar mass",
            ));
        }
        let water_gas_constant_j_kg_k =
            UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / species.molar_mass().value();
        let reference_vapor_enthalpy_j_kg = reference.vapor_specific_enthalpy().value();
        let liquid_energy_offset_j_kg = finite(
            source.species_specific_enthalpy_j_kg(&water, reference_temperature_k)?
                - reference_vapor_enthalpy_j_kg,
        )?;
        Ok(Self {
            source: Arc::new(source),
            water,
            reference_temperature_k,
            reference_vapor_enthalpy_j_kg,
            liquid_energy_offset_j_kg,
            water_gas_constant_j_kg_k,
        })
    }
    /// Original full NASA-9/transport catalog, including depleted species.
    pub fn source(&self) -> &ThermochemicalGasModel {
        &self.source
    }
    /// Explicit canonical water identity.
    pub const fn water_species(&self) -> &SpeciesId {
        &self.water
    }
    /// Chosen common vapor-enthalpy reference temperature [K].
    pub const fn reference_temperature_k(&self) -> f64 {
        self.reference_temperature_k
    }
    /// Add this CONSTANT to original IAPWS liquid h AND u to share the gas datum.
    /// Consumers using native IAPWS storage must account for it on mass transfer.
    pub const fn liquid_energy_offset_j_kg(&self) -> f64 {
        self.liquid_energy_offset_j_kg
    }
    /// Latent enthalpy in the aligned model [J/kg], without subtracting two
    /// potentially large formation references. Not an additional heat source.
    pub fn latent_enthalpy_j_kg(&self, temperature_k: f64) -> Result<f64, WaterCondensationError> {
        let liquid = saturation(temperature_k)?;
        let gap = finite(
            self.reference_vapor_enthalpy_j_kg
                + self.source.species_specific_enthalpy_change_j_kg(
                    &self.water,
                    self.reference_temperature_k,
                    temperature_k,
                )?
                - liquid.liquid_specific_enthalpy().value(),
        )?;
        if gap <= 0.0 {
            return Err(bad(
                "aligned water sources have no positive vaporization gap",
            ));
        }
        Ok(gap)
    }
    fn latent_internal_energy_j_kg(
        &self,
        liquid: SaturatedWaterState,
    ) -> Result<f64, WaterCondensationError> {
        let gap = finite(
            self.latent_enthalpy_j_kg(liquid.temperature().value())?
                - self.water_gas_constant_j_kg_k * liquid.temperature().value()
                + liquid.pressure().value() / liquid.liquid_density().value(),
        )?;
        if gap <= 0.0 {
            return Err(bad(
                "aligned water sources have no positive latent internal energy",
            ));
        }
        Ok(gap)
    }
}

/// One gas and its finite, isothermal liquid-water store at fixed GAS volume.
/// Canonical species masses are retained between trials to prevent composition
/// renormalization from accumulating carrier-inventory drift.
#[derive(Clone, Debug, PartialEq)]
pub struct WetGasParcel {
    model: WaterCondensationModel,
    gas: IsochoricGasParcel,
    species_mass_kg: BTreeMap<SpeciesId, f64>,
    total_water_mass_kg: f64,
    liquid_mass_kg: f64,
    liquid_saturation: SaturatedWaterState,
}
impl WetGasParcel {
    /// Retain an original gas and nonnegative liquid inventory at the same T.
    /// Initial supersaturation or dry gas alongside liquid is allowed: the next
    /// endpoint request resolves that nonequilibrium inventory without deleting it.
    pub fn new(
        model: WaterCondensationModel,
        gas: IsochoricGasParcel,
        liquid_mass_kg: f64,
    ) -> Result<Self, WaterCondensationError> {
        if !liquid_mass_kg.is_finite() || liquid_mass_kg < 0.0 {
            return Err(bad("finite nonnegative liquid inventory required"));
        }
        let liquid_saturation = saturation(gas.state().state().temperature)?;
        admit_pressure(gas.state().state().pressure)?;
        let closed = ClosedGasParcel::from_isochoric(gas.clone());
        let mut species_mass_kg = BTreeMap::new();
        for species in model.source.thermodynamics().components() {
            species_mass_kg.insert(
                species.species().clone(),
                closed.species_mass_kg(species.species())?,
            );
        }
        if !species_mass_kg
            .iter()
            .any(|(id, mass)| id != &model.water && *mass > 0.0)
        {
            return Err(bad("mixed-gas condensation needs a noncondensing carrier"));
        }
        // Reuse catalog/source/pack admission. This trial is not published and
        // therefore cannot replace the original inlet's represented state.
        closed.at_species_inventory_temperature(
            &model.source,
            &species_mass_kg,
            gas.state().state().temperature,
        )?;
        let vapor_mass_kg = species_mass_kg[&model.water];
        let total_water_mass_kg = finite(vapor_mass_kg + liquid_mass_kg)?;
        if (vapor_mass_kg > 0.0 && total_water_mass_kg == liquid_mass_kg)
            || (liquid_mass_kg > 0.0 && total_water_mass_kg == vapor_mass_kg)
        {
            return Err(bad(
                "positive water inventory is lost in its represented total",
            ));
        }
        model.latent_internal_energy_j_kg(liquid_saturation)?;
        let parcel = Self {
            model,
            gas,
            species_mass_kg,
            total_water_mass_kg,
            liquid_mass_kg,
            liquid_saturation,
        };
        finite(parcel.total_internal_energy_j())?;
        Ok(parcel)
    }
    /// Same original gas owner for acoustic, transport and heat consumers.
    pub const fn gas(&self) -> &IsochoricGasParcel {
        &self.gas
    }
    /// Retained full sources and common caloric reference.
    pub const fn model(&self) -> &WaterCondensationModel {
        &self.model
    }
    /// Common gas/film temperature [K].
    pub fn temperature_k(&self) -> f64 {
        self.gas.state().state().temperature
    }
    /// Actual stored liquid water [kg]; zero after complete evaporation.
    pub const fn liquid_mass_kg(&self) -> f64 {
        self.liquid_mass_kg
    }
    /// Conserved vapor plus liquid inventory [kg].
    pub const fn total_water_mass_kg(&self) -> f64 {
        self.total_water_mass_kg
    }
    /// Canonical retained vapor mass [kg].
    pub fn vapor_mass_kg(&self) -> f64 {
        self.species_mass_kg[&self.model.water]
    }
    /// Canonical gas species inventories [kg], including an exact absent-water zero.
    pub const fn species_masses_kg(&self) -> &BTreeMap<SpeciesId, f64> {
        &self.species_mass_kg
    }
    /// Original unshifted IAPWS liquid caloric and saturation source.
    pub const fn liquid_saturation(&self) -> SaturatedWaterState {
        self.liquid_saturation
    }
    /// Absolute gas plus aligned liquid internal energy [J]. Use the candidate's
    /// direct increment for thermal balances; subtracting these totals loses
    /// small changes beside arbitrary chemical formation datums.
    pub fn total_internal_energy_j(&self) -> f64 {
        self.gas.mass_kg() * self.gas.state().specific_internal_energy_j_kg()
            + self.liquid_mass_kg
                * (self
                    .liquid_saturation
                    .liquid_specific_internal_energy()
                    .value()
                    + self.model.liquid_energy_offset_j_kg)
    }

    /// Propose equilibrium at T, preserving gas volume and every original species.
    ///
    /// At fixed gas volume the saturation capacity is p_sat*V/(R_water*T),
    /// independent of carrier partial pressure. Using total pressure here would
    /// incorrectly suppress condensation. All retained liquid is available for
    /// evaporation, but no external water reservoir is assumed. Transport and
    /// acoustics are resolved through the original gas owner at the endpoint.
    /// No accepted caller state changes on cancellation or any refusal.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterCondensationStep, WaterCondensationError> {
        poll(&mut cancelled)?;
        let liquid = saturation(temperature_k)?;
        let capacity = finite(
            (liquid.pressure().value() / temperature_k) * self.gas.volume_m3()
                / self.model.water_gas_constant_j_kg_k,
        )?;
        if capacity <= 0.0 {
            return Err(bad("saturation inventory is unrepresentable"));
        }
        let vapor_mass = self.total_water_mass_kg.min(capacity);
        let liquid_mass = self.total_water_mass_kg - vapor_mass;
        let condensed_mass = self.vapor_mass_kg() - vapor_mass;
        let liquid_mass_change = liquid_mass - self.liquid_mass_kg;
        // A large film must not hide a lost transfer behind a film-sized mass
        // tolerance while the energy ledger credits the requested latent heat.
        let phase_roundoff =
            256.0 * f64::EPSILON * (self.vapor_mass_kg() + vapor_mass + condensed_mass.abs());
        if (condensed_mass != 0.0 && liquid_mass_change == 0.0)
            || (liquid_mass_change - condensed_mass).abs() > phase_roundoff
        {
            return Err(bad(
                "phase transfer is not representable in the liquid inventory",
            ));
        }
        let mut masses = self.species_mass_kg.clone();
        masses.insert(self.model.water.clone(), vapor_mass);
        let old = ClosedGasParcel::from_isochoric(self.gas.clone());
        let next =
            old.at_species_inventory_temperature(&self.model.source, &masses, temperature_k)?;
        admit_pressure(next.state().state().pressure)?;
        let mut maximum_species_mass_residual_kg: f64 = 0.0;
        for (id, target) in &masses {
            poll(&mut cancelled)?;
            maximum_species_mass_residual_kg =
                maximum_species_mass_residual_kg.max((next.species_mass_kg(id)? - target).abs());
        }
        let mass_residual_kg = finite((next.mass_kg() - self.gas.mass_kg()) + liquid_mass_change)?;
        let roundoff_mass = 256.0 * f64::EPSILON * (self.gas.mass_kg() + next.mass_kg());
        if mass_residual_kg.abs().max(maximum_species_mass_residual_kg) > roundoff_mass {
            return Err(bad(
                "phase endpoint failed original species/total mass conservation",
            ));
        }
        let gas_sensible = self.gas.mass_kg()
            * self
                .gas
                .model()
                .thermodynamics()
                .specific_internal_energy_change(
                    Temperature::new(self.temperature_k()),
                    Temperature::new(temperature_k),
                )
                .map_err(|error| {
                    GasEnergyError::Material(MaterialError::Parameters {
                        what: format!("condensation gas caloric increment: {error}"),
                    })
                })?
                .value();
        let liquid_sensible = self.liquid_mass_kg
            * (liquid.liquid_specific_internal_energy().value()
                - self
                    .liquid_saturation
                    .liquid_specific_internal_energy()
                    .value());
        let gap_u = self.model.latent_internal_energy_j_kg(liquid)?;
        let total_internal_energy_change_j =
            finite(gas_sensible + liquid_sensible - condensed_mass * gap_u)?;
        let water_u = self
            .model
            .source
            .species_specific_enthalpy_j_kg(&self.model.water, temperature_k)?
            - self.model.water_gas_constant_j_kg_k * temperature_k;
        let gas_internal_energy_change_j = finite(gas_sensible - condensed_mass * water_u)?;
        let liquid_internal_energy_change_j = finite(
            liquid_sensible
                + condensed_mass
                    * (liquid.liquid_specific_internal_energy().value()
                        + self.model.liquid_energy_offset_j_kg),
        )?;
        let liquid_enthalpy_change_j = finite(
            self.liquid_mass_kg
                * (liquid.liquid_specific_enthalpy().value()
                    - self.liquid_saturation.liquid_specific_enthalpy().value())
                + condensed_mass
                    * (liquid.liquid_specific_enthalpy().value()
                        + self.model.liquid_energy_offset_j_kg),
        )?;
        let saturation_residual_pa =
            next.species_partial_pressure_pa(&self.model.water)? - liquid.pressure().value();
        if saturation_residual_pa > 256.0 * f64::EPSILON * liquid.pressure().value() {
            return Err(bad(
                "accepted gas exceeds its original water saturation pressure",
            ));
        }
        poll(&mut cancelled)?;
        let parcel = Self {
            model: self.model.clone(),
            gas: next.into_isochoric(),
            species_mass_kg: masses,
            total_water_mass_kg: self.total_water_mass_kg,
            liquid_mass_kg: liquid_mass,
            liquid_saturation: liquid,
        };
        finite(parcel.total_internal_energy_j())?;
        Ok(WaterCondensationStep {
            parcel,
            condensed_mass_kg: condensed_mass,
            gas_internal_energy_change_j,
            liquid_internal_energy_change_j,
            liquid_enthalpy_change_j,
            total_internal_energy_change_j,
            latent_enthalpy_j: finite(
                condensed_mass * self.model.latent_enthalpy_j_kg(temperature_k)?,
            )?,
            mass_residual_kg,
            maximum_species_mass_residual_kg,
            saturation_residual_pa,
        })
    }
}

/// Complete, immutable gas/film proposal and independently inspectable balances.
#[derive(Clone, Debug, PartialEq)]
pub struct WaterCondensationStep {
    /// Proposed complete state; publish atomically with the thermal partner.
    pub parcel: WetGasParcel,
    /// Positive vapor -> liquid, negative finite-film evaporation [kg].
    pub condensed_mass_kg: f64,
    /// Original NASA gas internal-energy change, including species datum [J].
    pub gas_internal_energy_change_j: f64,
    /// Aligned liquid internal-energy change, including transferred mass [J].
    pub liquid_internal_energy_change_j: f64,
    /// Aligned liquid enthalpy change for an explicit enthalpy consumer [J].
    /// This is a different coordinate from rigid-volume internal energy.
    pub liquid_enthalpy_change_j: f64,
    /// Reference-invariant combined Delta U; required net heat at fixed volume [J].
    /// Uses sensible increments and a datum-cancelled latent gap. This is the
    /// thermal root's target, not a subtraction of absolute formation energies.
    pub total_internal_energy_change_j: f64,
    /// Diagnostic dm*L_h(T) [J], ALREADY contained in the caloric account.
    pub latent_enthalpy_j: f64,
    /// Actual gas plus liquid mass change [kg].
    pub mass_residual_kg: f64,
    /// Maximum emitted-gas versus canonical species-mass difference [kg].
    pub maximum_species_mass_residual_kg: f64,
    /// Actual water partial pressure minus source saturation pressure [Pa].
    /// Negative for an unsaturated state after exhausting the liquid inventory.
    pub saturation_residual_pa: f64,
}

#[cfg(test)]
mod tests;
