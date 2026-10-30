//! Explicit Peng--Robinson binding over the existing thermochemical species bank.
//!
//! The ideal-gas `ThermochemicalGasModel::evaluate` remains unchanged. This
//! module exposes nonideal mechanics separately so consumers cannot combine PR
//! density with ideal-gas sound speed or caloric departures by accident.

use fs_matdb::QueryPoint;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::eos::peng_robinson::{
    PengRobinsonComponentV1, PengRobinsonEvaluationV1, PengRobinsonMixtureV1,
    flash::{PengRobinsonFlashOptionsV1, PengRobinsonFlashV1},
};
use fs_thermochem::mixture::FrozenIdealGasMixtureEvaluationV1;

use super::{ThermochemicalGasModel, coordinates, failure, invalid};
use crate::MaterialError;

#[derive(Clone, Debug, PartialEq)]
pub struct PengRobinsonSpeciesBinding {
    pub species: SpeciesId,
    pub critical_temperature: Temperature,
    pub critical_pressure: Pressure,
    pub acentric_factor: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ThermochemicalPengRobinsonModel {
    source: ThermochemicalGasModel,
    mechanics: PengRobinsonMixtureV1,
}

impl ThermochemicalPengRobinsonModel {
    pub fn new(
        source: ThermochemicalGasModel,
        mut bindings: Vec<PengRobinsonSpeciesBinding>,
        binary_interactions: Vec<f64>,
    ) -> Result<Self, MaterialError> {
        let components = source.thermodynamics().components();
        if bindings.len() != components.len() {
            return Err(invalid(
                "real-gas binding requires exactly one critical-property record per active species",
            ));
        }
        bindings.sort_by(|a, b| a.species.cmp(&b.species));
        if bindings.iter().zip(components).any(
            |(binding, component)| &binding.species != component.species()
        ) {
            return Err(invalid(
                "real-gas critical properties contain duplicate, missing or foreign species",
            ));
        }
        let pr_components = bindings.iter().zip(components).map(
            |(binding, component)| PengRobinsonComponentV1 {
                critical_temperature: binding.critical_temperature,
                critical_pressure: binding.critical_pressure,
                acentric_factor: binding.acentric_factor,
                molar_mass: component.molar_mass(),
            }
        ).collect();
        let mechanics = PengRobinsonMixtureV1::new(
            pr_components,
            source.thermodynamics().mole_composition().fractions().to_vec(),
            binary_interactions,
        ).map_err(|error| failure("Peng-Robinson binding", error))?;
        Ok(Self { source, mechanics })
    }

    #[must_use]
    pub const fn source(&self) -> &ThermochemicalGasModel { &self.source }

    #[must_use]
    pub const fn mechanics(&self) -> &PengRobinsonMixtureV1 { &self.mechanics }

    pub fn evaluate_mechanics(
        &self,
        point: &QueryPoint,
    ) -> Result<ResolvedThermochemicalPengRobinsonState, MaterialError> {
        let (temperature, pressure) = coordinates(point)?;
        let thermodynamics = self.source.thermodynamics().evaluate(temperature, pressure)
            .map_err(|error| failure("nonideal source thermodynamics", error))?;
        let mechanics = self.mechanics.evaluate_pt(temperature, pressure)
            .map_err(|error| failure("Peng-Robinson mechanics", error))?;
        Ok(ResolvedThermochemicalPengRobinsonState {
            temperature, pressure, thermodynamics, mechanics,
        })
    }

    pub fn two_phase_flash(
        &self,
        point: &QueryPoint,
        options: PengRobinsonFlashOptionsV1,
    ) -> Result<PengRobinsonFlashV1, MaterialError> {
        let (temperature, pressure) = coordinates(point)?;
        self.source.thermodynamics().evaluate(temperature, pressure)
            .map_err(|error| failure("nonideal source thermodynamics", error))?;
        self.mechanics.two_phase_flash_pt(temperature, pressure, options)
            .map_err(|error| failure("Peng-Robinson flash", error))
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermochemicalPengRobinsonState {
    temperature: Temperature,
    pressure: Pressure,
    thermodynamics: FrozenIdealGasMixtureEvaluationV1,
    mechanics: PengRobinsonEvaluationV1,
}
impl ResolvedThermochemicalPengRobinsonState {
    #[must_use] pub const fn temperature(&self) -> Temperature { self.temperature }
    #[must_use] pub const fn pressure(&self) -> Pressure { self.pressure }
    #[must_use] pub const fn thermodynamics(&self) -> &FrozenIdealGasMixtureEvaluationV1 {
        &self.thermodynamics
    }
    #[must_use] pub const fn mechanics(&self) -> &PengRobinsonEvaluationV1 { &self.mechanics }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermochemical_gas::tests::{model, point};

    fn bindings() -> Vec<PengRobinsonSpeciesBinding> {
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B").unwrap(),
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A").unwrap(),
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
        ]
    }

    #[test]
    fn canonical_species_bind_to_original_composition_and_molar_mass() {
        let source = model(vec![0.75, 0.25], 0.0);
        let bound = ThermochemicalPengRobinsonModel::new(
            source.clone(), bindings(), vec![]
        ).unwrap();
        let state = bound.evaluate_mechanics(&point(300.0, 5.0e6)).unwrap();
        assert_eq!(
            state.mechanics().molar_mass().value(),
            state.thermodynamics().properties().molar_mass().value()
        );
        assert_ne!(
            state.mechanics().phase_roots()[0].compressibility_factor().to_bits(),
            1.0f64.to_bits()
        );
        assert_eq!(bound.source(), &source);
    }

    #[test]
    fn wrapped_mechanics_match_direct_pr_owner_exactly() {
        let bound = ThermochemicalPengRobinsonModel::new(
            model(vec![0.75, 0.25], 0.0),
            bindings(),
            vec![0.0, 0.03, 0.03, 0.0],
        ).unwrap();
        let query = point(280.0, 8.0e6);
        let wrapped = bound.evaluate_mechanics(&query).unwrap();
        let direct = bound.mechanics().evaluate_pt(
            wrapped.temperature(), wrapped.pressure()
        ).unwrap();
        assert_eq!(*wrapped.mechanics(), direct);
    }

    #[test]
    fn foreign_species_refuses_instead_of_reassigning_critical_data() {
        let mut bad = bindings();
        bad[0].species = SpeciesId::new("FOREIGN").unwrap();
        assert!(ThermochemicalPengRobinsonModel::new(
            model(vec![0.75, 0.25], 0.0), bad, vec![]
        ).is_err());
    }
}
