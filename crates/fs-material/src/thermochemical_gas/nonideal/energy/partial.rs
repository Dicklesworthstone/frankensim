//! Original NASA species enthalpies plus their mixture-dependent PR departures.
//!
//! A selective stream carries partial enthalpy at its actual upstream state,
//! not pure-species or mixture-average real-gas enthalpy. This constitutive
//! query supplies no species rate, separation work or entropy-production law.

use super::{GasEnergyError, ThermochemicalPengRobinsonModel};
use crate::thermochemical_gas::failure;
use fs_qty::{Density, Temperature};
use fs_thermochem::eos::peng_robinson::{
    MAX_PENG_ROBINSON_COMPONENTS_V1, partial_caloric::PengRobinsonPartialCaloricEvaluationV1,
};

/// Partial specific enthalpies [J/kg of species] in the original canonical
/// species order. The borrowed model retains species, NASA and PR identity.
/// Fixed-capacity storage and the underlying EOS query allocate no memory.
#[derive(Debug)]
pub struct RealGasPartialCalorics<'a> {
    model: &'a ThermochemicalPengRobinsonModel,
    departure: PengRobinsonPartialCaloricEvaluationV1<'a>,
    enthalpies_j_kg: [f64; MAX_PENG_ROBINSON_COMPONENTS_V1],
}

impl RealGasPartialCalorics<'_> {
    pub const fn model(&self) -> &ThermochemicalPengRobinsonModel {
        self.model
    }
    pub fn partial_specific_enthalpies_j_kg(&self) -> &[f64] {
        &self.enthalpies_j_kg[..self.departure.components().len()]
    }
    /// Original fixed-T,p partial molar volumes and enthalpy departures.
    pub const fn departure(&self) -> &PengRobinsonPartialCaloricEvaluationV1<'_> {
        &self.departure
    }
}

impl ThermochemicalPengRobinsonModel {
    /// Resolve original partial species enthalpies at one explicit density/T.
    /// The original homogeneous mechanical and complete-Cv checks apply. NASA
    /// source support is checked; transport is checked when resolving an actual
    /// material endpoint. No phase selection or composition perturbation occurs.
    pub fn partial_calorics_at_density(
        &self,
        density_kg_m3: f64,
        temperature_k: f64,
    ) -> Result<RealGasPartialCalorics<'_>, GasEnergyError> {
        let temperature = Temperature::new(temperature_k);
        let departure = self
            .mechanics()
            .evaluate_density_partial_calorics(temperature, Density::new(density_kg_m3))
            .map_err(|error| failure("real-gas partial calorics", error))?;
        self.caloric_from_eos(departure.state())?;
        let mut enthalpies_j_kg = [0.0; MAX_PENG_ROBINSON_COMPONENTS_V1];
        for (i, (source, residual)) in self
            .source()
            .thermodynamics()
            .components()
            .iter()
            .zip(departure.components())
            .enumerate()
        {
            let ideal = source
                .evaluate_properties(temperature)
                .map_err(|error| failure("real-gas partial ideal enthalpy", error))?;
            let h = (ideal.h().value() + residual.residual_enthalpy_j_mol())
                / source.molar_mass().value();
            if !h.is_finite() {
                return Err(GasEnergyError::Invalid(
                    "partial species enthalpy is unrepresentable",
                ));
            }
            enthalpies_j_kg[i] = h;
        }
        Ok(RealGasPartialCalorics {
            model: self,
            departure,
            enthalpies_j_kg,
        })
    }
}
