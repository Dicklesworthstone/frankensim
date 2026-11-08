//! Source-bound finite primary creep with independent thermal recovery.
//!
//! The three exact selectors below add `creep_reference_recovery_rate_s` [1/s]
//! and `creep_recovery_activation_temperature_k` [K] to the corresponding
//! primary Norton, primary Norton/J2, or primary Norton/J2/Voce thermal model.
//! Law version 1 and state schema 3 retain the original inverse inelastic
//! metric, permanent plastic/creep paths, and an independent resistance. A
//! zero recovery coefficient still selects that complete state schema.
//!
//! Recovery has its own Arrhenius activation, independent of creep viscosity.
//! It changes kinetic resistance, never the permanent metric or stored energy.
//! Kirchhoff stress, original-volume energy, logarithmic expansion, thermal
//! derivatives and all source bounds retain the parent family's conventions.

use super::{
    K, ResolvedAffineThermalFiniteJ2, THERMAL_FINITE_PRIMARY_NORTON_J2_LAW_ID,
    THERMAL_FINITE_PRIMARY_NORTON_J2_VOCE_LAW_ID, THERMAL_FINITE_PRIMARY_NORTON_LAW_ID, parameter,
    parse_law, refusal,
};
use crate::MaterialError;
use crate::plastic::finite::creep::recovery::FiniteRecoveringCreepState;
use crate::plastic::finite::thermal::{
    AffineThermalRecoveringCreep, FiniteCreepThermalRecovery, FiniteJ2ThermalError,
    FiniteJ2ThermalIncrement, FiniteThermalInelasticMaterial, FrozenThermalRecoveringCreep,
};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, QueryPoint};
use fs_qty::Dims;

/// Sole finite primary creep with a separately recovering kinetic resistance.
pub const THERMAL_FINITE_RECOVERING_NORTON_LAW_ID: &str =
    "finite-recovering-primary-norton-hencky-kirchhoff-thermal";
/// Simultaneous finite recovering creep and energetic linear J2 hardening.
pub const THERMAL_FINITE_RECOVERING_NORTON_J2_LAW_ID: &str =
    "finite-recovering-primary-norton-j2-hencky-kirchhoff-thermal";
/// Simultaneous finite recovering creep and energetic Voce J2 hardening.
pub const THERMAL_FINITE_RECOVERING_NORTON_J2_VOCE_LAW_ID: &str =
    "finite-recovering-primary-norton-j2-voce-hencky-kirchhoff-thermal";

const RECOVERY_RATE: &str = "creep_reference_recovery_rate_s";
const RECOVERY_ACTIVATION: &str = "creep_recovery_activation_temperature_k";

/// Exact original card, reference context and complete recovering finite law.
/// No ordinary finite-J2 adapter can silently discard the resistance history.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedAffineThermalRecoveringCreep {
    source: ResolvedAffineThermalFiniteJ2,
    law: AffineThermalRecoveringCreep,
}

impl ResolvedAffineThermalRecoveringCreep {
    pub fn try_from_model_card_at(
        card: &ConstitutiveModelCard,
        context: &QueryPoint,
        temperature_axis: &str,
    ) -> Result<Self, MaterialError> {
        card.validate()
            .map_err(|error| refusal(error.to_string()))?;
        let (base_selector, base_schema) = match card.law.0.as_str() {
            THERMAL_FINITE_RECOVERING_NORTON_LAW_ID => (THERMAL_FINITE_PRIMARY_NORTON_LAW_ID, 1),
            THERMAL_FINITE_RECOVERING_NORTON_J2_LAW_ID => {
                (THERMAL_FINITE_PRIMARY_NORTON_J2_LAW_ID, 2)
            }
            THERMAL_FINITE_RECOVERING_NORTON_J2_VOCE_LAW_ID => {
                (THERMAL_FINITE_PRIMARY_NORTON_J2_VOCE_LAW_ID, 2)
            }
            _ => {
                return Err(refusal(
                    "an exact finite thermal recovering Norton selector is required",
                ));
            }
        };
        if card.law_version != 1
            || card.state_schema_version != 3
            || card.initial_state != InitialStatePolicy::ZeroInternalState
        {
            return Err(refusal(
                "finite thermal recovery requires law v1, complete state v3 and virgin internal history",
            ));
        }
        let recovery = FiniteCreepThermalRecovery {
            reference_recovery_rate_s: parameter(card, RECOVERY_RATE, Dims([0, 0, -1, 0, 0, 0]))?,
            activation_temperature_k: parameter(card, RECOVERY_ACTIVATION, K)?,
        };
        // This private parameter projection reuses the exact native parser. It
        // never becomes a source receipt: admission below retains the ORIGINAL
        // card, law identity, state schema, hashes and applicability domain.
        let mut base_card = card.clone();
        base_card.law = LawId(base_selector.into());
        base_card.state_schema_version = base_schema;
        base_card.parameters.remove(RECOVERY_RATE);
        base_card.parameters.remove(RECOVERY_ACTIVATION);
        let base = parse_law(&base_card)?;
        let law = AffineThermalRecoveringCreep::new(base.clone(), recovery)
            .map_err(|error| refusal(format!("inadmissible finite thermal recovery: {error}")))?;
        let source = ResolvedAffineThermalFiniteJ2::admit_family(
            card,
            context,
            temperature_axis,
            base,
            None,
        )?;
        Ok(Self { source, law })
    }

    pub const fn model_card(&self) -> &ConstitutiveModelCard {
        self.source.model_card()
    }
    pub const fn context(&self) -> &QueryPoint {
        self.source.context()
    }
    pub fn temperature_axis(&self) -> &str {
        self.source.temperature_axis()
    }
    pub const fn law(&self) -> &AffineThermalRecoveringCreep {
        &self.law
    }
}

impl FiniteThermalInelasticMaterial for ResolvedAffineThermalRecoveringCreep {
    type State = FiniteRecoveringCreepState;
    type Frozen = FrozenThermalRecoveringCreep;

    fn initial_state(&self) -> Self::State {
        self.law.initial_state()
    }
    fn reference_temperature_k(&self) -> f64 {
        self.law.reference_temperature_k()
    }
    fn temperature_range_k(&self) -> [f64; 2] {
        self.law.temperature_range_k()
    }
    fn minimum_mechanical_heat_capacity_j_m3_k(&self) -> f64 {
        self.law.minimum_mechanical_heat_capacity_j_m3_k()
    }
    fn at_temperature_for_inspection(
        &self,
        temperature_k: f64,
    ) -> Result<Self::Frozen, FiniteJ2ThermalError> {
        self.law.at_temperature_for_inspection(temperature_k)
    }
    fn at_temperature_for_step(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<Self::Frozen, FiniteJ2ThermalError> {
        self.law.at_temperature_for_step(temperature_k, duration_s)
    }
    fn checked_increment_for_step(
        &self,
        old_f: &[f64; 9],
        old: &Self::State,
        old_temperature_k: f64,
        new_f: &[f64; 9],
        new_temperature_k: f64,
        yield_tolerance_pa: f64,
        duration_s: f64,
    ) -> Result<FiniteJ2ThermalIncrement<Self::State>, FiniteJ2ThermalError> {
        self.law.checked_increment_for_step(
            old_f,
            old,
            old_temperature_k,
            new_f,
            new_temperature_k,
            yield_tolerance_pa,
            duration_s,
        )
    }
}
