//! Effective carrier drag, retaining either original material cards or the
//! electron momentum-transfer cross sections evaluated at the actual state.

use super::{
    AVOGADRO_PER_MOL, FrozenIonizedGasStateV1, IonizedGasError, MOMENTUM_RELAXATION_RATE_DIMS,
    SpeciesId, poll, positive,
};
use crate::state_point::ResolvedMaterialStatePoint;
use fs_matdb::PropertyKey;
use fs_math::det::sqrt;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Pressure, QuantitySpec, Temperature};
use fs_thermochem::electron_collision::{
    ElectronCollisionEvaluationV1, ElectronCollisionKindV1, MaxwellianRateWeightingV1,
};

/// Complete original evidence for one effective momentum-relaxation rate.
#[derive(Clone, Debug, PartialEq)]
pub enum CarrierMomentumSource {
    /// Caller-selected effective drag against its authored background.
    Material {
        material: ResolvedMaterialStatePoint,
        property: PropertyKey,
    },
    /// Weak-drift shifted-Maxwellian electron drag against every present neutral
    /// in this exact frozen state. This is a fluid moment closure, not an exact
    /// solution of the energy-dependent Lorentz-Boltzmann conductivity problem.
    ElectronNeutral {
        state: FrozenIonizedGasStateV1,
        contributions: Vec<ElectronNeutralMomentumContribution>,
        /// Joint bound on sqrt(me/M_target) and target/electron thermal speed.
        maximum_target_thermal_speed_ratio: f64,
    },
}

/// One retained neutral target's contribution to electron momentum relaxation.
#[derive(Clone, Debug, PartialEq)]
pub struct ElectronNeutralMomentumContribution {
    evaluation: ElectronCollisionEvaluationV1,
    number_density_per_m3: f64,
    rate_per_s: f64,
    target_thermal_speed_ratio: f64,
}
impl ElectronNeutralMomentumContribution {
    pub const fn evaluation(&self) -> &ElectronCollisionEvaluationV1 {
        &self.evaluation
    }
    pub const fn number_density_per_m3(&self) -> f64 {
        self.number_density_per_m3
    }
    pub const fn rate_per_s(&self) -> f64 {
        self.rate_per_s
    }
    /// sqrt((Th/m_target)/(Te/m_e)); the neutral thermal speed relative to electrons.
    pub const fn target_thermal_speed_ratio(&self) -> f64 {
        self.target_thermal_speed_ratio
    }
}

/// One positive effective carrier/background momentum-relaxation rate [1/s].
/// Material applicability includes its authored mixture/background; a card is
/// not inferred from species names. Cross-section-derived rates bind the full
/// frozen population as well as both actual temperatures.
#[derive(Clone, Debug, PartialEq)]
pub struct CarrierMomentumRelaxation {
    species: SpeciesId,
    source: CarrierMomentumSource,
    heavy_temperature_k: f64,
    electron_temperature_k: f64,
    pressure_pa: f64,
    rate_per_s: f64,
}
impl CarrierMomentumRelaxation {
    /// Select a positive static rate from a source-resolved bundle with typed
    /// `T` [K] and `p` [Pa]. This common-temperature declaration cannot be used
    /// when the electron and heavy species have different temperatures.
    pub fn new(
        species: SpeciesId,
        material: ResolvedMaterialStatePoint,
        property: PropertyKey,
    ) -> Result<Self, IonizedGasError> {
        Self::from_material(species, material, property, ["T", "T"])
    }

    /// Select a rate explicitly resolved at typed `T_h`, `T_e` [K] and `p` [Pa].
    /// Both temperatures must match the actual source-bound state, including
    /// for a rate whose authored value is constant throughout that domain.
    pub fn new_two_temperature(
        species: SpeciesId,
        material: ResolvedMaterialStatePoint,
        property: PropertyKey,
    ) -> Result<Self, IonizedGasError> {
        Self::from_material(species, material, property, ["T_h", "T_e"])
    }

    fn from_material(
        species: SpeciesId,
        material: ResolvedMaterialStatePoint,
        property: PropertyKey,
        temperature_axes: [&str; 2],
    ) -> Result<Self, IonizedGasError> {
        let invalid = IonizedGasError::Invalid;
        let source = material.property_by_key(&property).ok_or(invalid(
            "explicit carrier momentum-relaxation source is absent",
        ))?;
        if property.quantity() != QuantitySpec::dimensional(MOMENTUM_RELAXATION_RATE_DIMS)
            || !source.is_static_at_fixed_external_conditions()
            || !positive(source.value_si())
        {
            return Err(invalid(
                "require a positive static momentum-relaxation source",
            ));
        }
        let rate_per_s = source.value_si();
        let coordinate = |axis, dims, kind| {
            let quantity = material.axis_quantities().get(axis);
            if !quantity.is_some_and(|q| {
                *q == QuantitySpec::dimensional(dims)
                    || *q == QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
            }) {
                return Err(invalid(
                    "collision coordinate has the wrong quantity meaning",
                ));
            }
            material
                .state_coordinate(axis)
                .filter(|value| positive(*value))
                .ok_or(invalid(
                    "collision state requires positive actual temperatures and pressure",
                ))
        };
        let heavy_temperature_k = coordinate(
            temperature_axes[0],
            Temperature::DIMS,
            QuantityKind::AbsoluteTemperature,
        )?;
        let electron_temperature_k = coordinate(
            temperature_axes[1],
            Temperature::DIMS,
            QuantityKind::AbsoluteTemperature,
        )?;
        let pressure_pa = coordinate("p", Pressure::DIMS, QuantityKind::Pressure)?;
        Ok(Self {
            species,
            source: CarrierMomentumSource::Material { material, property },
            heavy_temperature_k,
            electron_temperature_k,
            pressure_pa,
            rate_per_s,
        })
    }

    /// Derive electron drag from source-retained momentum cross sections and
    /// the actual density of every present neutral. The rate owner supplies
    /// `(2/3)<sigma_mt v E/(kTe)>`; an event rate `<sigma v>` is refused.
    ///
    /// The stationary-heavy approximation requires an explicit positive joint
    /// bound r below one on both sqrt(me/m_target) and
    /// sqrt((Th/m_target)/(Te/me)). Thus cold targets cannot evade the massive
    /// recoil-target requirement. These are approximation-domain bounds, not
    /// certified errors. Every target must pass both.
    /// This excludes Coulomb/electron-ion drag, non-Maxwellian distributions and
    /// magnetic response; those effects cannot be supplied as extra neutral rows.
    pub fn from_electron_neutral_cross_sections(
        state: &FrozenIonizedGasStateV1,
        mut evaluations: Vec<ElectronCollisionEvaluationV1>,
        maximum_target_thermal_speed_ratio: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IonizedGasError> {
        let invalid = IonizedGasError::Invalid;
        poll(&mut cancelled)?;
        if !positive(maximum_target_thermal_speed_ratio)
            || maximum_target_thermal_speed_ratio >= 1.0
            || evaluations.len() > state.amounts().len()
        {
            return Err(invalid(
                "bounded stationary-target speed ratio and source count required",
            ));
        }
        let electron = state
            .electron_index()
            .ok_or(invalid("structural electron model required"))?;
        if !state.number_density_per_m3(electron).is_some_and(positive) {
            return Err(invalid(
                "electron drag requires a present electron population",
            ));
        }
        let chemistry = state.thermodynamics();
        let models = chemistry.models();
        let electron_mass = models[electron].molar_mass().value() / AVOGADRO_PER_MOL;
        let th = state.heavy_temperature().value();
        let te = state.electron_temperature().value();
        evaluations.sort_by(|a, b| a.target_species().cmp(b.target_species()));
        if evaluations
            .windows(2)
            .any(|p| p[0].target_species() == p[1].target_species())
        {
            return Err(invalid("duplicate neutral momentum-transfer target"));
        }
        let expected = evaluations.len();
        let mut evaluations = evaluations.into_iter();
        let mut contributions = Vec::with_capacity(expected);
        let mut rate_per_s = 0.0;
        for (i, model) in models.iter().enumerate() {
            poll(&mut cancelled)?;
            let density = state
                .number_density_per_m3(i)
                .ok_or(invalid("missing target density"))?;
            if chemistry.charge_vector().get(i) != Some(0) || density == 0.0 {
                continue;
            }
            let evaluation = evaluations
                .next()
                .filter(|e| e.target_species() == model.species())
                .ok_or(invalid(
                    "one momentum cross section is required for every present neutral",
                ))?;
            if evaluation.weighting() != MaxwellianRateWeightingV1::WeakDriftMomentum
                || evaluation.source().kind() != ElectronCollisionKindV1::ElasticMomentumTransfer
                || evaluation.electron_temperature() != state.electron_temperature()
                || (electron_mass / evaluation.electron_mass_kg() - 1.0).abs() > 16.0 * f64::EPSILON
            {
                return Err(invalid(
                    "momentum weighting, electron mass and actual Te must agree",
                ));
            }
            let coefficient = evaluation
                .coefficient_m3_s()
                .filter(|v| positive(*v))
                .ok_or(invalid(
                    "positive representable momentum rate coefficient required",
                ))?;
            let target_mass = model.molar_mass().value() / AVOGADRO_PER_MOL;
            let mass_ratio_root = sqrt(electron_mass / target_mass);
            let speed_ratio = sqrt(th / te) * mass_ratio_root;
            let contribution = density * coefficient;
            if !mass_ratio_root.is_finite()
                || mass_ratio_root > maximum_target_thermal_speed_ratio
                || !speed_ratio.is_finite()
                || speed_ratio > maximum_target_thermal_speed_ratio
                || !positive(contribution)
            {
                return Err(invalid(
                    "neutral target violates stationary-background or rate domain",
                ));
            }
            rate_per_s += contribution;
            contributions.push(ElectronNeutralMomentumContribution {
                evaluation,
                number_density_per_m3: density,
                rate_per_s: contribution,
                target_thermal_speed_ratio: speed_ratio,
            });
        }
        if contributions.len() != expected || !positive(rate_per_s) {
            return Err(invalid(
                "complete nonempty neutral target coverage and finite drag required",
            ));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            species: models[electron].species().clone(),
            source: CarrierMomentumSource::ElectronNeutral {
                state: state.clone(),
                contributions,
                maximum_target_thermal_speed_ratio,
            },
            heavy_temperature_k: th,
            electron_temperature_k: te,
            pressure_pa: state.pressure().value(),
            rate_per_s,
        })
    }

    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }
    /// Original material source or complete evaluated cross-section evidence.
    pub const fn source(&self) -> &CarrierMomentumSource {
        &self.source
    }
    pub const fn rate_per_s(&self) -> f64 {
        self.rate_per_s
    }
    pub const fn heavy_temperature_k(&self) -> f64 {
        self.heavy_temperature_k
    }
    pub const fn electron_temperature_k(&self) -> f64 {
        self.electron_temperature_k
    }
    pub const fn pressure_pa(&self) -> f64 {
        self.pressure_pa
    }

    pub(super) fn matches_state(&self, state: &FrozenIonizedGasStateV1, pressure_pa: f64) -> bool {
        self.heavy_temperature_k == state.heavy_temperature().value()
            && self.electron_temperature_k == state.electron_temperature().value()
            && match &self.source {
                CarrierMomentumSource::Material { .. } => self.pressure_pa == pressure_pa,
                CarrierMomentumSource::ElectronNeutral {
                    state: original, ..
                } => original == state,
            }
    }
}
