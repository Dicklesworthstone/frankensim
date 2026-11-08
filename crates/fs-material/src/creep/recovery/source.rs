//! Original source binding for primary creep with static resistance recovery.
//!
//! Recovery has its own reference rate, activation energy and reference
//! temperature. None is inferred from Norton flow, plasticity or a melting
//! point. All thirteen original claims are pinned together, and every knot of
//! the admitted temperature interval is checked before an interval is prepared.
//! The two Arrhenius laws own their separate temperature dependences; the
//! reference coefficients and elastic storage stay fixed over this interval.
//! Explicit thermoelastic construction instead admits sourced affine E/nu through
//! the original thermal potential, without changing either kinetic memory law.

use std::sync::Arc;

use super::{CreepRecoveryParameters, RecoveringCreepState, RecoveringNortonCreep};
use crate::creep::source::{
    CREEP_ACTIVATION_ENERGY_DIMS, CREEP_RATE_DIMS, CreepSourceError, CreepSourceFamily,
    NortonCreepSource, NortonSourceQuery, norton_strain_hardening_requirements,
};
use crate::creep::{CreepError, CreepIncrement, CreepLaw, CreepTolerances, NortonCreep};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::ModelCard;
use fs_matdb::{ClaimSelection, MaterialCard};
use fs_qty::{QuantitySpec, Temperature};

/// Independent static resistance-recovery rate at its reference temperature [1/s].
/// Zero explicitly selects no recovery; an absent claim is never treated as zero.
pub const CREEP_RECOVERY_REFERENCE_RATE_PROPERTY: &str = "creep_recovery_reference_rate";
/// Independent molar activation energy of resistance recovery [J/mol].
pub const CREEP_RECOVERY_ACTIVATION_ENERGY_PROPERTY: &str = "creep_recovery_activation_energy";
/// Absolute temperature of the recovery reference rate [K].
pub const CREEP_RECOVERY_REFERENCE_TEMPERATURE_PROPERTY: &str =
    "creep_recovery_reference_temperature";

/// Complete thirteen-claim selection: primary Norton flow and independent
/// Arrhenius resistance recovery. This is model input, not fitted-data evidence.
pub fn recovering_norton_requirements()
-> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
    let mut requirements = norton_strain_hardening_requirements()?;
    for (name, dims, domain) in [
        (
            CREEP_RECOVERY_REFERENCE_RATE_PROPERTY,
            CREEP_RATE_DIMS,
            ScalarAdmissibility::NonNegative,
        ),
        (
            CREEP_RECOVERY_ACTIVATION_ENERGY_PROPERTY,
            CREEP_ACTIVATION_ENERGY_DIMS,
            ScalarAdmissibility::NonNegative,
        ),
        (
            CREEP_RECOVERY_REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        ),
    ] {
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, domain)?);
    }
    Ok(requirements)
}

#[derive(Debug)]
struct Binding<'a> {
    card: &'a MaterialCard,
    creep: NortonCreepSource<'a>,
    reference: ResolvedMaterialStatePoint,
    recovery: CreepRecoveryParameters,
}

/// Immutable source family shared by cells, physical steps and thermal trials.
/// It owns no specimen history, so resolving an endpoint cannot age a specimen.
#[derive(Clone, Debug)]
pub struct RecoveringNortonSource<'a> {
    binding: Arc<Binding<'a>>,
}

impl<'a> RecoveringNortonSource<'a> {
    /// Admit the complete source interval. All original claim identities and
    /// fixed external coordinates are preserved when temperature changes.
    pub fn new(
        card: &'a MaterialCard,
        query: NortonSourceQuery<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CreepSourceError> {
        Self::prepare(card, query, false, cancelled)
    }

    // Complete smooth elasticity and expansion admission belongs to the public
    // thermal wrapper. This intermediate source refuses bare continuation.
    pub(crate) fn prepare_thermoelastic(
        card: &'a MaterialCard,
        query: NortonSourceQuery<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CreepSourceError> {
        Self::prepare(card, query, true, cancelled)
    }

    fn prepare(
        card: &'a MaterialCard,
        query: NortonSourceQuery<'_>,
        thermal_elasticity: bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CreepSourceError> {
        poll(&mut cancelled)?;
        let reference = resolve_material_state_point(
            card,
            query.lower,
            &recovering_norton_requirements()?,
            query.selection,
        )?;
        let selection = MaterialPropertySelection::PinnedByProperty(
            norton_strain_hardening_requirements()?
                .into_iter()
                .map(|requirement| {
                    let property = reference
                        .property(requirement.name())
                        .expect("combined source contains every primary Norton coefficient");
                    (
                        requirement.name().into(),
                        property.answer().receipt.selected,
                    )
                })
                .collect(),
        );
        let creep_query = NortonSourceQuery {
            lower: query.lower,
            upper: query.upper,
            temperature_axis: query.temperature_axis,
            selection,
            strain_limit: query.strain_limit,
        };
        let creep = if thermal_elasticity {
            NortonCreepSource::prepare_thermoelastic(card, creep_query, true, &mut cancelled)?
        } else {
            NortonCreepSource::new_with_strain_hardening(card, creep_query, &mut cancelled)?
        };
        for name in [
            CREEP_RECOVERY_REFERENCE_RATE_PROPERTY,
            CREEP_RECOVERY_ACTIVATION_ENERGY_PROPERTY,
            CREEP_RECOVERY_REFERENCE_TEMPERATURE_PROPERTY,
        ] {
            poll(&mut cancelled)?;
            let property = reference
                .property(name)
                .expect("combined source contains every recovery coefficient");
            if !property.has_only_fixed_external_conditions(Some(query.temperature_axis), false) {
                return Err(CreepSourceError::Invalid(
                    "recovery coefficients cannot freeze evolving stress, strain, rate or frequency axes",
                ));
            }
            let selected = property.answer().receipt.selected;
            let claim = card
                .claims()
                .claim(selected)
                .ok_or(CreepSourceError::Invalid(
                    "selected recovery coefficient is absent from the original card",
                ))?;
            let line = card
                .claims()
                .query_line_typed(
                    &claim.key,
                    query.temperature_axis,
                    query.lower,
                    query.upper,
                    ClaimSelection::Pinned(selected),
                )
                .map_err(|error| CreepSourceError::Support {
                    property: name.into(),
                    reason: error.to_string(),
                })?;
            for (index, (_, answer)) in line.knots.iter().enumerate() {
                if index % 256 == 0 {
                    poll(&mut cancelled)?;
                }
                if answer.evidence.value.value != property.value_si() {
                    return Err(CreepSourceError::Invalid(
                        "recovery reference coefficients must be constant over the complete thermal span",
                    ));
                }
            }
        }
        let value = |name| {
            reference
                .property(name)
                .expect("required recovery coefficient")
                .value_si()
        };
        let recovery = CreepRecoveryParameters {
            reference_rate_s: value(CREEP_RECOVERY_REFERENCE_RATE_PROPERTY),
            activation_energy_j_mol: value(CREEP_RECOVERY_ACTIVATION_ENERGY_PROPERTY),
            reference_temperature_k: value(CREEP_RECOVERY_REFERENCE_TEMPERATURE_PROPERTY),
        };
        let source = Self {
            binding: Arc::new(Binding {
                card,
                creep,
                reference,
                recovery,
            }),
        };
        for temperature_k in source.creep_source().parameters().temperature_range_k {
            if thermal_elasticity {
                source.resolve_temperature(temperature_k, 1.0, true)?;
            } else {
                RecoveringNortonCreep::new(
                    source.creep_source().prepared_law(temperature_k, 1.0)?,
                    recovery,
                )?;
            }
        }
        poll(&mut cancelled)?;
        Ok(source)
    }

    /// Every original flow, resistance and elastic source receipt.
    pub fn reference(&self) -> &ResolvedMaterialStatePoint {
        &self.binding.reference
    }

    pub fn creep_source(&self) -> &NortonCreepSource<'a> {
        &self.binding.creep
    }

    pub fn recovery_parameters(&self) -> &CreepRecoveryParameters {
        &self.binding.recovery
    }

    /// Prepare one physical interval, refreshing the full pinned claim bundle
    /// once. Trial stress evaluations use the prepared native scalar return and
    /// do not revisit material tables or perform a second nonlinear solve.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<ResolvedRecoveringNorton, CreepSourceError> {
        self.resolve_temperature(temperature_k, duration_s, false)
    }

    fn resolve_temperature(
        &self,
        temperature_k: f64,
        duration_s: f64,
        with_thermal_storage: bool,
    ) -> Result<ResolvedRecoveringNorton, CreepSourceError> {
        // Ordinary source calls retain their native-before-receipt validation.
        let fixed_law = if with_thermal_storage {
            None
        } else {
            Some(RecoveringNortonCreep::new(
                self.binding.creep.prepared_law(temperature_k, duration_s)?,
                self.binding.recovery,
            )?)
        };
        let point = self
            .reference()
            .query()
            .with(self.binding.creep.temperature_axis(), temperature_k)
            .map_err(|error| CreepSourceError::Support {
                property: self.binding.creep.temperature_axis().into(),
                reason: error.to_string(),
            })?;
        let material = self.reference().at_query_point(self.binding.card, &point)?;
        let law = match fixed_law {
            Some(law) => law,
            None => RecoveringNortonCreep::new(
                self.binding.creep.prepared_thermoelastic_law(
                    &material,
                    temperature_k,
                    duration_s,
                )?,
                self.binding.recovery,
            )?,
        };
        Ok(ResolvedRecoveringNorton {
            material: Arc::new(material),
            law,
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CreepSourceError> {
    if cancelled() {
        Err(CreepSourceError::Cancelled)
    } else {
        Ok(())
    }
}

/// Native recovery law and its complete, original source bundle.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedRecoveringNorton {
    material: Arc<ResolvedMaterialStatePoint>,
    law: RecoveringNortonCreep,
}

impl ResolvedRecoveringNorton {
    pub const fn law(&self) -> &RecoveringNortonCreep {
        &self.law
    }
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
}

impl SmallStrainLaw for ResolvedRecoveringNorton {
    type State = RecoveringCreepState;
    fn initial_state(&self) -> Self::State {
        self.law.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        self.law.stress(strain, state)
    }
    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        self.law.tangent(strain, state)
    }
    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        self.law.update_state(strain, state)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.law.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.law.card()
    }
}

impl CreepLaw for ResolvedRecoveringNorton {
    fn kinetics(&self) -> &NortonCreep {
        self.law.kinetics()
    }
    fn validate_creep_history(
        &self,
        strain: &Voigt,
        state: &Self::State,
    ) -> Result<(), CreepError> {
        self.law.validate_creep_history(strain, state)
    }
    fn creep_increment(
        &self,
        old_strain: &Voigt,
        old: &Self::State,
        strain: &Voigt,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<Self::State>, CreepError> {
        self.law
            .creep_increment(old_strain, old, strain, tolerances)
    }
}

impl<'a> CreepSourceFamily<'a> for RecoveringNortonSource<'a> {
    type Material = ResolvedRecoveringNorton;
    fn creep_source(&self) -> &NortonCreepSource<'a> {
        RecoveringNortonSource::creep_source(self)
    }
    fn at_temperature(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<Self::Material, CreepSourceError> {
        RecoveringNortonSource::at_temperature(self, temperature_k, duration_s)
    }
    fn at_thermoelastic_temperature(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<Self::Material, CreepSourceError> {
        self.resolve_temperature(temperature_k, duration_s, true)
    }
}
