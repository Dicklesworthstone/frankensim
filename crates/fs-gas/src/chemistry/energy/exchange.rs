//! Shared collision inputs and the original neutral thermal moment.

use super::*;

/// Original collision laws and explicit physical domains for reactive thermal
/// evolution. Numerical time integration remains the reacting solver's owner.
/// A supplied Coulomb model binds the initial populations exactly; each later
/// evaluation explicitly refreshes it at the actual source-bound state.
#[derive(Clone, Debug, PartialEq)]
pub struct ElectronHeavyCollisionSources {
    /// One original elastic momentum source for every populated neutral target
    /// and every neutral target that the reaction can create.
    pub momentum_sources: Vec<ElectronCrossSectionV1>,
    /// Original Coulomb screening, domain limits and initial population receipt.
    pub coulomb: Option<CoulombModelV1>,
    /// Charged mole-fraction limit, in (0,1]; one requires a Coulomb source.
    pub maximum_charged_mole_fraction: f64,
    /// Joint recoil and target/electron thermal-speed bound, in (0,1).
    pub maximum_target_thermal_speed_ratio: f64,
    /// Original-source Maxwellian integration controls.
    pub rate_controls: MaxwellianRateControlsV1,
}

#[derive(Clone, Copy, Debug, Default)]
pub(in crate::chemistry) struct ExchangeConductance {
    pub neutral: f64,
    pub coulomb: f64,
}
impl ExchangeConductance {
    pub fn total(self) -> f64 {
        self.neutral + self.coulomb
    }
}

/// Source matching is done once, outside the nonlinear/time integration loops.
pub(in crate::chemistry) struct PreparedExchange<'a> {
    original: &'a FrozenIonizedGasStateV1,
    inputs: &'a ElectronHeavyCollisionSources,
    targets: [Option<&'a ElectronCrossSectionV1>; MAX_EQUILIBRIUM_SPECIES_V1],
    minimum_electron_heavy_ratio: f64,
}
impl<'a> PreparedExchange<'a> {
    pub fn new(
        initial: &'a FrozenIonizedGasStateV1,
        inputs: &'a ElectronHeavyCollisionSources,
        created_neutral: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ElectronHeavyEnergyError> {
        poll(cancelled)?;
        if !positive(inputs.maximum_charged_mole_fraction)
            || inputs.maximum_charged_mole_fraction > 1.0
            || (inputs.coulomb.is_none() && inputs.maximum_charged_mole_fraction == 1.0)
            || !positive(inputs.maximum_target_thermal_speed_ratio)
            || inputs.maximum_target_thermal_speed_ratio >= 1.0
            || inputs
                .coulomb
                .as_ref()
                .is_some_and(|m| !m.is_bound_to_populations(initial))
        {
            return Err(ElectronHeavyEnergyError::Invalid(
                "collision domain or original Coulomb population binding",
            ));
        }
        let bank = initial.thermodynamics();
        let species = bank.elemental_matrix().species();
        if inputs.momentum_sources.len() > species.len() {
            return Err(ElectronHeavyEnergyError::Invalid(
                "too many neutral collision sources",
            ));
        }
        let mut prepared = Self {
            original: initial,
            inputs,
            targets: [None; MAX_EQUILIBRIUM_SPECIES_V1],
            minimum_electron_heavy_ratio: 0.0,
        };
        for source in &inputs.momentum_sources {
            poll(cancelled)?;
            let index = species
                .binary_search(source.target_species())
                .map_err(|_| {
                    ElectronHeavyEnergyError::Invalid("collision target absent from source bank")
                })?;
            if source.kind() != ElectronCollisionKindV1::ElasticMomentumTransfer
                || bank.charge_vector().get(index) != Some(0)
                || prepared.targets[index].is_some()
            {
                return Err(ElectronHeavyEnergyError::Invalid(
                    "one elastic momentum source per neutral target",
                ));
            }
            prepared.targets[index] = Some(source);
        }
        let has_electrons = initial
            .electron_index()
            .is_some_and(|i| initial.amounts()[i].value() > 0.0);
        if has_electrons {
            let electron = initial
                .electron_index()
                .expect("present structural electrons");
            let mass = bank.models()[electron].molar_mass().value() / AVOGADRO_PER_MOL;
            if ((mass - ELECTRON_MASS_KG) / ELECTRON_MASS_KG).abs() > 16.0 * f64::EPSILON {
                return Err(ElectronHeavyEnergyError::Invalid(
                    "physical collision electron mass",
                ));
            }
            for (i, amount) in initial.amounts().iter().enumerate() {
                if bank.charge_vector().get(i) != Some(0)
                    || (amount.value() == 0.0 && i != created_neutral)
                {
                    continue;
                }
                if prepared.targets[i].is_none() {
                    return Err(ElectronHeavyEnergyError::MissingSource(species[i].clone()));
                }
                let mass = bank.models()[i].molar_mass().value() / AVOGADRO_PER_MOL;
                let ratio =
                    (ELECTRON_MASS_KG / mass) / inputs.maximum_target_thermal_speed_ratio.powi(2);
                if !ratio.is_finite() || ratio > 1.0 {
                    return Err(ElectronHeavyEnergyError::OutsideDomain);
                }
                prepared.minimum_electron_heavy_ratio =
                    prepared.minimum_electron_heavy_ratio.max(ratio);
            }
        }
        prepared.coupling(initial, cancelled)?;
        Ok(prepared)
    }

    pub const fn minimum_electron_heavy_ratio(&self) -> f64 {
        self.minimum_electron_heavy_ratio
    }

    pub fn coupling(
        &self,
        state: &FrozenIonizedGasStateV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ExchangeConductance, ElectronHeavyEnergyError> {
        if state.thermodynamics() != self.original.thermodynamics() {
            return Err(ElectronHeavyEnergyError::Invalid(
                "collision source bank changed",
            ));
        }
        let mut total = 0.0;
        let mut charged = 0.0;
        for (i, amount) in state.amounts().iter().enumerate() {
            total += amount.value();
            if state.thermodynamics().charge_vector().get(i) != Some(0) {
                charged += amount.value();
            }
        }
        if !positive(total) || charged / total > self.inputs.maximum_charged_mole_fraction {
            return Err(ElectronHeavyEnergyError::OutsideDomain);
        }
        let coulomb = if let Some(source) = &self.inputs.coulomb {
            source
                .evaluate_at_state(state, &mut *cancelled)
                .map_err(coulomb_error)?
                .electron_heavy_thermal_conductance_w_m3_k()
                * state.volume().value()
        } else {
            0.0
        };
        let neutral = neutral_conductance(
            state,
            &self.targets,
            state.heavy_temperature().value(),
            state.electron_temperature().value(),
            self.inputs.maximum_target_thermal_speed_ratio,
            self.inputs.rate_controls,
            cancelled,
        )?;
        if !coulomb.is_finite() || coulomb < 0.0 || !(neutral + coulomb).is_finite() {
            return Err(ElectronHeavyEnergyError::Arithmetic(
                "electron-heavy conductance",
            ));
        }
        Ok(ExchangeConductance { neutral, coulomb })
    }
}

/// W/K for the parcel. Both frozen and reacting owners use this same leading
/// massive-target Maxwellian moment; no second neutral relaxation law exists.
pub(super) fn neutral_conductance(
    state: &FrozenIonizedGasStateV1,
    sources: &[Option<&ElectronCrossSectionV1>; MAX_EQUILIBRIUM_SPECIES_V1],
    th: f64,
    te: f64,
    maximum_target_ratio: f64,
    rate_controls: MaxwellianRateControlsV1,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ElectronHeavyEnergyError> {
    let Some(electron) = state.electron_index() else {
        return Ok(0.0);
    };
    let electron_amount = state.amounts()[electron].value();
    if electron_amount == 0.0 {
        return Ok(0.0);
    }
    let mut frequency = 0.0;
    let mut correction = 0.0;
    for (i, amount) in state.amounts().iter().enumerate() {
        if amount.value() == 0.0 || state.thermodynamics().charge_vector().get(i) != Some(0) {
            continue;
        }
        poll(cancelled)?;
        let source = sources[i].ok_or_else(|| {
            ElectronHeavyEnergyError::MissingSource(
                state.thermodynamics().models()[i].species().clone(),
            )
        })?;
        let mass = state.thermodynamics().models()[i].molar_mass().value() / AVOGADRO_PER_MOL;
        let ratio_squared = (th / te) * (ELECTRON_MASS_KG / mass);
        if fs_math::det::sqrt(ELECTRON_MASS_KG / mass) > maximum_target_ratio
            || !ratio_squared.is_finite()
            || ratio_squared > maximum_target_ratio.powi(2)
        {
            return Err(ElectronHeavyEnergyError::OutsideDomain);
        }
        let evaluation = source
            .evaluate_properties(
                Temperature::new(te),
                MaxwellianRateWeightingV1::WeakDriftMomentum,
                rate_controls,
                &mut *cancelled,
            )
            .map_err(|error| match error {
                ElectronCollisionErrorV1::Cancelled => ElectronHeavyEnergyError::Cancelled,
                error => ElectronHeavyEnergyError::Collision {
                    species: source.target_species().clone(),
                    source: error,
                },
            })?;
        let coefficient =
            evaluation
                .coefficient_m3_s()
                .ok_or(ElectronHeavyEnergyError::Arithmetic(
                    "unrepresentable elastic momentum coefficient",
                ))?;
        let density = state
            .number_density_per_m3(i)
            .expect("canonical source target");
        let increment = (density * coefficient) * (ELECTRON_MASS_KG / mass) - correction;
        let next = frequency + increment;
        correction = (next - frequency) - increment;
        frequency = next;
    }
    // R=kB*NA; molar amount avoids an unnecessary total particle count.
    let conductance = (3.0 * R * electron_amount) * frequency;
    if !conductance.is_finite() || conductance < 0.0 {
        return Err(ElectronHeavyEnergyError::Arithmetic(
            "neutral thermal conductance",
        ));
    }
    Ok(conductance)
}
