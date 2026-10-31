//! The same finite water partition at prescribed pressure rather than volume.

use super::*;
use crate::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_matdb::QueryPoint;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Amount, QuantitySpec};
use fs_thermochem::{Composition, CompositionBasis};

/// Source-resolved equilibrium of an explicitly supplied finite water inventory.
/// A steady transport owner may use a declared reference-duration packet and
/// divide its masses and joules by that duration. This object is not a stored
/// wall film, a flow-rate law, or an isochoric parcel.
#[derive(Clone, Debug, PartialEq)]
pub struct WaterIsobaricPartition {
    /// Actual source-bound gas at the caller's pressure and temperature.
    pub gas: ResolvedThermochemicalGasState,
    /// Full catalog-order mole composition, including exact absent-species zeros.
    pub mole_composition: Composition,
    /// Canonical gas masses, with every original noncondensing species retained.
    pub gas_species_masses_kg: BTreeMap<SpeciesId, f64>,
    /// Sum of the canonical gas inventory [kg].
    pub gas_mass_kg: f64,
    /// Actual remaining liquid inventory [kg].
    pub liquid_mass_kg: f64,
    /// Positive vapor-to-liquid transfer, negative finite-film evaporation [kg].
    pub condensed_mass_kg: f64,
    /// Original saturated-liquid h plus the model's constant NASA datum [J/kg].
    pub liquid_specific_enthalpy_j_kg: f64,
    /// Reference-independent vapor-minus-liquid specific enthalpy [J/kg].
    pub latent_enthalpy_j_kg: f64,
    /// Phase contribution to total H at THIS fixed T: -dm*L [J].
    /// The thermal owner adds the original species/liquid sensible increments;
    /// it must not add this energy twice when transporting liquid enthalpy.
    pub isothermal_phase_enthalpy_change_j: f64,
    /// Actual total gas-plus-liquid mass change [kg].
    pub mass_residual_kg: f64,
    /// Largest emitted-gas versus canonical species-inventory difference [kg].
    pub maximum_species_mass_residual_kg: f64,
    /// Actual water partial pressure minus the original saturation pressure [Pa].
    pub saturation_residual_pa: f64,
}

impl WaterCondensationModel {
    /// Partition finite water at an explicitly prescribed common pressure and T.
    ///
    /// Carrier amounts remain fixed. Ideal-mixture equilibrium gives
    /// n_water,sat = n_carrier * p_sat/(p-p_sat). Retained liquid can evaporate
    /// only up to its inventory; excess vapor condenses. This pressure constraint
    /// is different from the fixed-gas-volume `WetGasParcel` constraint.
    ///
    /// The original IAPWS source and the same 273.16..=323.15 K, 20..=200 kPa
    /// approximation apply. The original gas owner resolves EOS and transport;
    /// there is no pressure-drop, stream-energy or kinetics model here. Callers
    /// supply source conditions without temperature or pressure and close their
    /// own enthalpy/work account before publishing any physical state.
    #[allow(clippy::too_many_arguments)]
    pub fn partition_at_pressure_temperature(
        &self,
        species_masses_kg: &BTreeMap<SpeciesId, f64>,
        liquid_mass_kg: f64,
        temperature_k: f64,
        pressure_pa: f64,
        source_conditions: &QueryPoint,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterIsobaricPartition, WaterCondensationError> {
        poll(&mut cancelled)?;
        let liquid = saturation(temperature_k)?;
        admit_pressure(pressure_pa)?;
        if !liquid_mass_kg.is_finite()
            || liquid_mass_kg < 0.0
            || source_conditions.axes().contains_key("temperature")
            || source_conditions.axes().contains_key("pressure")
            || species_masses_kg
                .values()
                .any(|mass| !mass.is_finite() || *mass < 0.0)
            || species_masses_kg.keys().any(|id| {
                !self
                    .source
                    .thermodynamics()
                    .components()
                    .iter()
                    .any(|s| s.species() == id)
            })
        {
            return Err(bad(
                "isobaric condensation needs original finite inventories and unowned source coordinates",
            ));
        }
        let original_water = species_masses_kg.get(&self.water).copied().unwrap_or(0.0);
        let total_water = finite(original_water + liquid_mass_kg)?;
        if (original_water > 0.0 && total_water == liquid_mass_kg)
            || (liquid_mass_kg > 0.0 && total_water == original_water)
        {
            return Err(bad(
                "isobaric total water cannot represent both phase inventories",
            ));
        }
        let mut carrier_moles = 0.0;
        let mut masses = BTreeMap::new();
        for species in self.source.thermodynamics().components() {
            poll(&mut cancelled)?;
            let mass = species_masses_kg
                .get(species.species())
                .copied()
                .unwrap_or(0.0);
            masses.insert(species.species().clone(), mass);
            if species.species() != &self.water {
                let amount = finite(mass / species.molar_mass().value())?;
                if mass > 0.0 && amount == 0.0 {
                    return Err(bad("isobaric carrier amount is unrepresentable"));
                }
                carrier_moles = finite(carrier_moles + amount)?;
            }
        }
        if carrier_moles <= 0.0 || pressure_pa <= liquid.pressure().value() {
            return Err(bad(
                "isobaric liquid-water partition needs a noncondensing carrier below boiling",
            ));
        }
        let water_molar_mass = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / self.water_gas_constant_j_kg_k;
        let capacity = finite(
            carrier_moles
                * (liquid.pressure().value() / (pressure_pa - liquid.pressure().value()))
                * water_molar_mass,
        )?;
        if capacity <= 0.0 {
            return Err(bad("isobaric saturation inventory is unrepresentable"));
        }
        let vapor = total_water.min(capacity);
        let remaining_liquid = total_water - vapor;
        let condensed = original_water - vapor;
        let stored_liquid_change = remaining_liquid - liquid_mass_kg;
        let transfer_roundoff =
            256.0 * f64::EPSILON * original_water.max(vapor).max(condensed.abs());
        if (condensed != 0.0 && stored_liquid_change == 0.0)
            || (stored_liquid_change - condensed).abs() > transfer_roundoff
        {
            return Err(bad(
                "isobaric phase transfer cannot be represented in liquid storage",
            ));
        }
        masses.insert(self.water.clone(), vapor);
        let gas_mass = finite(masses.values().sum())?;
        let original_gas_mass = finite(species_masses_kg.values().sum())?;
        let mass_residual = finite((gas_mass - original_gas_mass) + stored_liquid_change)?;
        let mass_roundoff = 256.0 * f64::EPSILON * gas_mass.max(original_gas_mass);
        if gas_mass <= 0.0 || mass_residual.abs() > mass_roundoff {
            return Err(bad(
                "isobaric phase endpoint failed total mass conservation",
            ));
        }
        let amounts: Vec<_> = self
            .source
            .thermodynamics()
            .components()
            .iter()
            .map(|species| Amount::new(masses[species.species()] / species.molar_mass().value()))
            .collect();
        let model = self.source.with_mole_amounts(&amounts)?;
        let active = model.thermodynamics();
        let full_fractions: Vec<_> = self
            .source
            .thermodynamics()
            .components()
            .iter()
            .map(|species| {
                active
                    .components()
                    .iter()
                    .position(|s| s.species() == species.species())
                    .map_or(0.0, |index| active.mole_composition().fractions()[index])
            })
            .collect();
        let mole_composition = Composition::new(CompositionBasis::MoleFraction, full_fractions)
            .map_err(|_| bad("isobaric canonical mole composition is unrepresentable"))?;
        let active_molar_masses: Vec<_> =
            active.components().iter().map(|s| s.molar_mass()).collect();
        let emitted = active
            .mole_composition()
            .to_mass_fractions(&active_molar_masses)
            .map_err(|_| bad("isobaric emitted species inventory is unrepresentable"))?;
        let mut maximum_species_mass_residual_kg: f64 = 0.0;
        for (species, &fraction) in active.components().iter().zip(emitted.fractions()) {
            maximum_species_mass_residual_kg = maximum_species_mass_residual_kg
                .max((gas_mass * fraction - masses[species.species()]).abs());
        }
        if maximum_species_mass_residual_kg > mass_roundoff {
            return Err(bad("isobaric source gas failed species conservation"));
        }
        let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
        let point = source_conditions
            .clone()
            .with_quantity(
                "temperature",
                typed(QuantityKind::AbsoluteTemperature),
                temperature_k,
            )
            .and_then(|point| {
                point.with_quantity("pressure", typed(QuantityKind::Pressure), pressure_pa)
            })
            .map_err(|_| bad("invalid isobaric material coordinates"))?;
        let gas = if model.water_species().is_some() {
            model.evaluate_with_iapws_liquid_water_limit(&point)?
        } else {
            model.evaluate(&point)?
        };
        let saturation_residual_pa =
            gas.state().water_mole_fraction * pressure_pa - liquid.pressure().value();
        let latent_enthalpy_j_kg = self.latent_enthalpy_j_kg(temperature_k)?;
        let liquid_specific_enthalpy_j_kg =
            finite(liquid.liquid_specific_enthalpy().value() + self.liquid_energy_offset_j_kg)?;
        poll(&mut cancelled)?;
        Ok(WaterIsobaricPartition {
            gas,
            mole_composition,
            gas_species_masses_kg: masses,
            gas_mass_kg: gas_mass,
            liquid_mass_kg: remaining_liquid,
            condensed_mass_kg: condensed,
            liquid_specific_enthalpy_j_kg,
            latent_enthalpy_j_kg,
            isothermal_phase_enthalpy_change_j: finite(-condensed * latent_enthalpy_j_kg)?,
            mass_residual_kg: mass_residual,
            maximum_species_mass_residual_kg,
            saturation_residual_pa,
        })
    }
}
