//! Conserved original-source noncondensable gas inside a thermal cavity.
//!
//! One immutable original parcel supplies its exact mass, composition, NASA-9
//! calorics, transport claims and auxiliary query conditions. The enclosing
//! bubble owns the current temperature and volume; this material view has no
//! independently accepted state or clock. Trial samples evaluate only the
//! original heat capacities and ideal-gas EOS, with no allocations or caloric
//! inversion on successful paths. Energy changes integrate the original Cv,
//! so arbitrary species formation references cannot change the dynamics.
//!
//! Full endpoint qualification re-evaluates the original source at the dry
//! component's partial pressure. Its transport and sound speed describe that
//! component, NOT a water-vapor mixture. A consumer must separately admit gas
//! inertia, continuum behavior, ideal-mixture equilibrium, internal thermal
//! and species mixing, and negligible dissolution of the conserved inventory.
//! In particular, dry-gas conductivity is not a wet-mixture transport law.

use crate::cavitation::invalid;
use crate::thermochemical_gas::ResolvedThermochemicalGasState;
use crate::thermochemical_gas::energy::parcel::ClosedGasParcel;
use crate::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_qty::Temperature;

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn source_error(stage: &str, error: impl core::fmt::Display) -> GasEnergyError {
    GasEnergyError::Material(invalid(&format!("conserved bubble gas {stage}: {error}")))
}

/// Immutable finite dry-gas inventory and its original caloric source domain.
///
/// The retained parcel is the reference inventory, not another evolved
/// thermal state. Its exact accepted mass never changes with T or bubble size.
#[derive(Clone, Debug, PartialEq)]
pub struct ConservedBubbleGas {
    reference: ClosedGasParcel,
    controls: GasEnergySolve,
    gas_constant_j_kg_k: f64,
}

/// Allocation-free original dry-component state for one bubble T and volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubbleGasSample {
    pub temperature_k: f64,
    pub volume_m3: f64,
    pub mass_kg: f64,
    pub density_kg_m3: f64,
    /// Dry-component partial pressure, not total wet-mixture pressure.
    pub pressure_pa: f64,
    pub gas_constant_j_kg_k: f64,
    pub cp_j_kg_k: f64,
    pub cv_j_kg_k: f64,
    /// Actual extensive dry-gas capacity m_g Cv(T), in J/K.
    pub heat_capacity_j_k: f64,
    pub gamma: f64,
    /// Frozen-composition sound speed of the DRY component alone.
    pub sound_speed_m_s: f64,
}

impl ConservedBubbleGas {
    /// Retain the exact original inventory and admit one continuous Cv branch.
    ///
    /// Water must be explicitly excluded from the source's noncondensable
    /// composition. Other species' gas-phase stability and negligible liquid
    /// solubility remain declared applicability conditions. The original
    /// caloric controls admit source continuity and positive capacity over the
    /// whole temperature bracket; no inverse is performed by trial samples.
    pub fn new(
        reference: ClosedGasParcel,
        controls: GasEnergySolve,
    ) -> Result<Self, GasEnergyError> {
        controls.validate()?;
        if reference.model().water_species().is_some() {
            return Err(GasEnergyError::Invalid(
                "conserved noncondensable bubble gas must exclude the declared water species",
            ));
        }
        let gas_constant = reference.state().state().specific_gas_constant;
        if !positive(reference.mass_kg())
            || !positive(reference.volume_m3())
            || !positive(gas_constant)
        {
            return Err(GasEnergyError::Invalid(
                "positive original dry-gas inventory, volume and gas constant required",
            ));
        }
        reference
            .model()
            .prepare_caloric_branch(controls, || false)?;
        let result = Self {
            reference,
            controls,
            gas_constant_j_kg_k: gas_constant,
        };
        result.sample(
            result.reference.state().state().temperature,
            result.reference.volume_m3(),
        )?;
        Ok(result)
    }

    /// Original immutable parcel, source identities and auxiliary conditions.
    pub const fn reference_parcel(&self) -> &ClosedGasParcel {
        &self.reference
    }

    /// Original continuous caloric admission bracket and numerical policy.
    pub const fn controls(&self) -> GasEnergySolve {
        self.controls
    }

    /// Exact conserved mass of the original noncondensable inventory [kg].
    pub const fn mass_kg(&self) -> f64 {
        self.reference.mass_kg()
    }

    /// Original frozen-composition specific gas constant [J/(kg K)].
    pub const fn gas_constant_j_kg_k(&self) -> f64 {
        self.gas_constant_j_kg_k
    }

    fn admit_temperature(&self, temperature_k: f64) -> Result<(), GasEnergyError> {
        if !temperature_k.is_finite()
            || temperature_k < self.controls.minimum_temperature_k
            || temperature_k > self.controls.maximum_temperature_k
        {
            return Err(GasEnergyError::Invalid(
                "dry-gas temperature outside the admitted original caloric branch",
            ));
        }
        Ok(())
    }

    /// Sample original Cp/Cv and the conserved-mass EOS without source receipts.
    ///
    /// Successful calls allocate nothing. Full source/transport qualification
    /// remains necessary before accepting an enclosing bubble endpoint.
    pub fn sample(
        &self,
        temperature_k: f64,
        volume_m3: f64,
    ) -> Result<BubbleGasSample, GasEnergyError> {
        self.admit_temperature(temperature_k)?;
        if !positive(volume_m3) {
            return Err(GasEnergyError::Invalid(
                "positive finite dry-gas volume required",
            ));
        }
        let capacity = self
            .reference
            .model()
            .thermodynamics()
            .evaluate_heat_capacities(Temperature::new(temperature_k))
            .map_err(|error| source_error("heat capacity", error))?;
        let cp = capacity.mass_specific_cp().value();
        let cv = capacity.mass_specific_cv().value();
        let gamma = cp / cv;
        let mass = self.mass_kg();
        let density = mass / volume_m3;
        let rt = self.gas_constant_j_kg_k * temperature_k;
        let pressure = density * rt;
        let extensive_capacity = mass * cv;
        let sound_speed = fs_math::det::sqrt(gamma * rt);
        if ![cp, cv, density, pressure, extensive_capacity, sound_speed]
            .into_iter()
            .all(positive)
            || !gamma.is_finite()
            || gamma <= 1.0
            || gamma > 5.0 / 3.0 + 1.0e-9
        {
            return Err(GasEnergyError::Invalid(
                "represented original gas capacity, EOS and frozen sound speed required",
            ));
        }
        Ok(BubbleGasSample {
            temperature_k,
            volume_m3,
            mass_kg: mass,
            density_kg_m3: density,
            pressure_pa: pressure,
            gas_constant_j_kg_k: self.gas_constant_j_kg_k,
            cp_j_kg_k: cp,
            cv_j_kg_k: cv,
            heat_capacity_j_k: extensive_capacity,
            gamma,
            sound_speed_m_s: sound_speed,
        })
    }

    /// Original m_g integral(Cv dT) [J], independent of formation references.
    ///
    /// Both endpoints and the connecting path belong to the same previously
    /// admitted continuous source branch. The enclosing physical owner must
    /// balance this energy against heat, vapor energy and total pressure work.
    pub fn energy_change_j(
        &self,
        initial_temperature_k: f64,
        final_temperature_k: f64,
    ) -> Result<f64, GasEnergyError> {
        self.admit_temperature(initial_temperature_k)?;
        self.admit_temperature(final_temperature_k)?;
        let change = self.mass_kg()
            * self
                .reference
                .model()
                .thermodynamics()
                .specific_internal_energy_change(
                    Temperature::new(initial_temperature_k),
                    Temperature::new(final_temperature_k),
                )
                .map_err(|error| source_error("energy increment", error))?
                .value();
        if !change.is_finite()
            || (initial_temperature_k != final_temperature_k && change == 0.0)
        {
            return Err(GasEnergyError::Invalid(
                "represented original dry-gas caloric increment required",
            ));
        }
        Ok(change)
    }

    /// Fully qualify an endpoint through the original material and query owner.
    ///
    /// Only T and dry partial p change in the retained typed query. The exact
    /// fixed mass remains owned here; resolving a component state never
    /// recreates its inventory from rounded pressure or total bubble volume.
    /// No state is modified on success or refusal.
    pub fn qualify(
        &self,
        temperature_k: f64,
        volume_m3: f64,
    ) -> Result<ResolvedThermochemicalGasState, GasEnergyError> {
        let sample = self.sample(temperature_k, volume_m3)?;
        let point = self
            .reference
            .point()
            .clone()
            .with("temperature", temperature_k)
            .and_then(|point| point.with("pressure", sample.pressure_pa))
            .map_err(|error| source_error("endpoint query", error))?;
        let resolved = self.reference.model().evaluate(&point)?;
        let relative = resolved.state().density / sample.density_kg_m3 - 1.0;
        if !relative.is_finite() || relative.abs() > 64.0 * f64::EPSILON {
            return Err(GasEnergyError::Invalid(
                "original dry-gas endpoint does not retain the conserved density",
            ));
        }
        Ok(resolved)
    }
}

#[cfg(test)]
#[path = "gas/tests.rs"]
mod tests;
