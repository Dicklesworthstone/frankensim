//! Reaction-specific bath collision laws on the original evolving gas inventory.
//!
//! The source card supplies nonnegative efficiencies epsilon_i and a positive
//! reference concentration c0 [mol/m3]. Actual moles/volume and the complete
//! current mixture give q = sum_i(epsilon_i c_i)/c0, including inert species.
//! Every catalog species has an explicit property binding; no default collider,
//! gas-name rule, duplicate species registry or caller-supplied pressure exists.
//!
//! Third-body kinetics multiply the original forward AND reverse rates by q.
//! Lindemann falloff treats the base scale as the high-pressure limit k_inf,
//! resolves a normalized low-pressure scale k0 from the same card, and applies
//! Pr/(1+Pr), where Pr = (k0/k_inf) q. Both scales have units mol/(m3 s) because
//! species activities and q are dimensionless. An imported concentration-based
//! low-pressure coefficient must include the explicit c0 and activity-basis
//! conversion. Explicit Troe3/Troe4 regimes multiply this same Lindemann factor
//! by the canonical Troe broadening function, using the actual gas temperature.
//! Three parameters omit the T2 term; four parameters require its source claim.
//! This bounded rung admits alpha in [0,1], positive characteristic temperatures
//! and 0 < Fcent <= 1. It does not clamp coefficients or supply SRI, chemically
//! activated kinetics or validated combustion-mechanism data.
//!
//! Reference conventions:
//! <https://www.cantera.org/stable/reference/kinetics/reaction-rates.html>
//! <https://www.cantera.org/stable/reference/kinetics/rate-constants.html>
//!
//! All coefficients, including currently absent species' efficiencies, use
//! original pinned material claims at the actual gas state. Their values and
//! sources stay in the original material bundle; the native kinetic receipt
//! records the common multiplier. Zero collision concentration gives a real
//! zero multiplier; a positive quantity lost to underflow refuses instead.
use super::arrhenius::MOLAR_CONCENTRATION_DIMS;
use super::{
    ACTIVITY_RATE_SCALE_DIMS, ClosedGasParcel, GasReactionError, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, invalid, poll, positive,
};
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_qty::{Dims, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::equilibrium::MAX_REACTION_SPECIES_V1;
use std::collections::BTreeSet;

/// Explicit positive normalization c0 in mol/m3.
pub const REFERENCE_CONCENTRATION_PROPERTY: &str = "reaction_collider_reference_concentration";
/// Normalized low-pressure activity scale k0 in mol/(m3 s), used by falloff laws.
/// This can be a source scalar, curve or surface at the actual state point.
pub const LOW_PRESSURE_SCALE_PROPERTY: &str = "reaction_low_pressure_rate_scale";
/// Dimensionless Troe weight alpha in [0,1].
pub const TROE_ALPHA_PROPERTY: &str = "reaction_troe_alpha";
/// Troe T3 (T***) in kelvin, in (1-alpha)*exp(-T/T3).
pub const TROE_T3_PROPERTY: &str = "reaction_troe_temperature_3";
/// Troe T1 (T*) in kelvin, in alpha*exp(-T/T1).
pub const TROE_T1_PROPERTY: &str = "reaction_troe_temperature_1";
/// Troe T2 (T**) in kelvin, in exp(-T2/T); required only for Troe4.
pub const TROE_T2_PROPERTY: &str = "reaction_troe_temperature_2";

/// A source property for ONE species' reaction-specific collision efficiency.
/// Multiple explicitly named species may share one property.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ColliderEfficiency {
    pub species: SpeciesId,
    pub property: String,
}

/// Caller-declared kinetic regime; neither is inferred from species names.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ColliderRegime {
    /// Multiply both existing directional rates by normalized collider density.
    ThirdBody,
    /// Interpolate low/high-pressure limits with the Lindemann Pr/(1+Pr) law.
    Lindemann,
    /// Troe broadening with explicitly sourced alpha, T3 and T1; no T2 term.
    Troe3,
    /// Troe broadening with explicitly sourced alpha, T3, T1 and T2.
    Troe4,
}

/// Immutable species/property bindings for a single reaction channel.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MaterialColliderLaw {
    regime: ColliderRegime,
    efficiencies: Vec<ColliderEfficiency>,
}
impl MaterialColliderLaw {
    /// Admit explicit bounded bindings, canonically ordered by species identity.
    /// Complete source-catalog coverage is checked when binding the reaction.
    pub fn new(
        regime: ColliderRegime,
        mut efficiencies: Vec<ColliderEfficiency>,
    ) -> Result<Self, GasReactionError> {
        if efficiencies.is_empty() || efficiencies.len() > MAX_REACTION_SPECIES_V1 {
            return Err(invalid(
                "collider law requires a bounded nonempty species binding",
            ));
        }
        efficiencies.sort_by(|a, b| a.species.cmp(&b.species));
        if efficiencies
            .windows(2)
            .any(|pair| pair[0].species == pair[1].species)
        {
            return Err(invalid("collider law repeats a species identity"));
        }
        let law = Self {
            regime,
            efficiencies,
        };
        let required = law.requirements()?;
        if required.len() > fs_material::state_point::MAX_MATERIAL_STATE_PROPERTIES {
            return Err(invalid(
                "collider coefficient requirements exceed the shared material bound",
            ));
        }
        Ok(law)
    }

    #[must_use]
    pub const fn regime(&self) -> ColliderRegime {
        self.regime
    }
    /// Original source bindings, in canonical species order.
    #[must_use]
    pub fn efficiencies(&self) -> &[ColliderEfficiency] {
        &self.efficiencies
    }

    pub(super) fn validate_catalog(
        &self,
        catalog: &ThermochemicalGasModel,
    ) -> Result<(), GasReactionError> {
        let species = catalog.thermodynamics().components();
        if species.len() != self.efficiencies.len()
            || species
                .iter()
                .zip(&self.efficiencies)
                .any(|(s, binding)| s.species() != &binding.species)
        {
            return Err(invalid(
                "collider bindings must cover every catalog species exactly once",
            ));
        }
        Ok(())
    }

    pub(super) fn requirements(&self) -> Result<Vec<ScalarPropertyRequirement>, GasReactionError> {
        let mut required = vec![ScalarPropertyRequirement::try_new(
            REFERENCE_CONCENTRATION_PROPERTY,
            MOLAR_CONCENTRATION_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?];
        if self.regime != ColliderRegime::ThirdBody {
            required.push(ScalarPropertyRequirement::try_new(
                LOW_PRESSURE_SCALE_PROPERTY,
                ACTIVITY_RATE_SCALE_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?);
        }
        if matches!(self.regime, ColliderRegime::Troe3 | ColliderRegime::Troe4) {
            required.push(ScalarPropertyRequirement::try_new(
                TROE_ALPHA_PROPERTY,
                Dims([0; 6]),
                ScalarAdmissibility::NonNegative,
            )?);
            for name in [TROE_T3_PROPERTY, TROE_T1_PROPERTY] {
                required.push(ScalarPropertyRequirement::try_new(
                    name,
                    Temperature::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                )?);
            }
            if self.regime == ColliderRegime::Troe4 {
                required.push(ScalarPropertyRequirement::try_new(
                    TROE_T2_PROPERTY,
                    Temperature::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                )?);
            }
        }
        let mut names = BTreeSet::new();
        for binding in &self.efficiencies {
            if binding.property == REFERENCE_CONCENTRATION_PROPERTY
                || binding.property == LOW_PRESSURE_SCALE_PROPERTY
                || required.iter().any(|r| r.name() == binding.property)
            {
                return Err(invalid(
                    "an efficiency property cannot also be a collider scale",
                ));
            }
            names.insert(binding.property.as_str());
        }
        for name in names {
            required.push(ScalarPropertyRequirement::try_new(
                name,
                Dims([0; 6]),
                ScalarAdmissibility::NonNegative,
            )?);
        }
        Ok(required)
    }

    pub(super) fn multiplier(
        &self,
        resolved: &ResolvedMaterialStatePoint,
        gas: &ClosedGasParcel,
        base_scale: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, GasReactionError> {
        poll(cancelled)?;
        let property = |name| {
            resolved
                .property(name)
                .map(|p| p.value_si())
                .ok_or_else(|| invalid("resolved source omitted a required collider coefficient"))
        };
        // Selected source validity and the physical Troe coefficient domain
        // still apply when the complete effective collider population is zero.
        let log_center = self.troe_log_center(resolved, gas.state().state().temperature)?;
        let Some(log_q) = super::third_body::normalized_concentration_log(
            &self.efficiencies,
            resolved,
            gas,
            cancelled,
        )?
        else {
            return Ok(0.0);
        };
        let factor = match self.regime {
            ColliderRegime::ThirdBody => fs_math::det::exp(log_q),
            ColliderRegime::Lindemann | ColliderRegime::Troe3 | ColliderRegime::Troe4 => {
                let low = property(LOW_PRESSURE_SCALE_PROPERTY)?;
                if !positive(base_scale) || !positive(low) {
                    return Err(invalid(
                        "positive low/high-pressure kinetic scales required",
                    ));
                }
                let log_reduced = fs_math::det::ln(low) - fs_math::det::ln(base_scale) + log_q;
                if !log_reduced.is_finite() {
                    return Err(invalid("falloff reduced pressure is unrepresentable"));
                }
                if let Some(log_center) = log_center {
                    troe_multiplier(log_reduced, log_center)?
                } else {
                    // Avoid overflow in Pr or 1+Pr at a perfectly finite high-
                    // pressure limit. A lost positive low-pressure rate still refuses.
                    if log_reduced >= 0.0 {
                        1.0 / (1.0 + fs_math::det::exp(-log_reduced))
                    } else {
                        let reduced = fs_math::det::exp(log_reduced);
                        reduced / (1.0 + reduced)
                    }
                }
            }
        };
        if !positive(factor) {
            return Err(invalid(
                "positive collider multiplier underflows or overflows",
            ));
        }
        poll(cancelled)?;
        Ok(factor)
    }

    fn troe_log_center(
        &self,
        resolved: &ResolvedMaterialStatePoint,
        temperature: f64,
    ) -> Result<Option<f64>, GasReactionError> {
        if !matches!(self.regime, ColliderRegime::Troe3 | ColliderRegime::Troe4) {
            return Ok(None);
        }
        let property = |name| {
            resolved
                .property(name)
                .map(|p| p.value_si())
                .ok_or_else(|| invalid("resolved source omitted a required Troe coefficient"))
        };
        let alpha = property(TROE_ALPHA_PROPERTY)?;
        if !alpha.is_finite() || !(0.0..=1.0).contains(&alpha) || !positive(temperature) {
            return Err(invalid(
                "Troe requires alpha in [0,1] and positive gas temperature",
            ));
        }
        let mut terms = [f64::NEG_INFINITY; 3];
        for (i, (weight, name)) in [(1.0 - alpha, TROE_T3_PROPERTY), (alpha, TROE_T1_PROPERTY)]
            .into_iter()
            .enumerate()
        {
            let characteristic = property(name)?;
            if !positive(characteristic) {
                return Err(invalid("Troe characteristic temperatures must be positive"));
            }
            if weight == 0.0 {
                continue;
            }
            let ratio = temperature / characteristic;
            if !positive(ratio) {
                return Err(invalid("Troe temperature ratio underflows or overflows"));
            }
            terms[i] = fs_math::det::ln(weight) - ratio;
        }
        if self.regime == ColliderRegime::Troe4 {
            let characteristic = property(TROE_T2_PROPERTY)?;
            let ratio = characteristic / temperature;
            if !positive(characteristic) || !positive(ratio) {
                return Err(invalid("Troe T2 temperature ratio underflows or overflows"));
            }
            terms[2] = -ratio;
        }
        // Keep even an exponentially small positive center in logarithmic
        // form. A final representable rate must not depend on exp(-T/Ti)
        // fitting separately in a direct floating-point value.
        let maximum = terms.into_iter().fold(f64::NEG_INFINITY, f64::max);
        if !maximum.is_finite() {
            return Err(invalid("Troe center has no representable logarithm"));
        }
        let sum: f64 = terms
            .into_iter()
            .filter(|v| v.is_finite())
            .map(|v| fs_math::det::exp(v - maximum))
            .sum();
        let log_center = maximum + fs_math::det::ln(sum);
        if !log_center.is_finite() || log_center > 0.0 {
            return Err(invalid("Troe broadening requires 0 < Fcent <= 1"));
        }
        Ok(Some(log_center))
    }
}

/// Canonical Troe formula, evaluated without constructing Pr, Fcent or x/y.
fn troe_multiplier(log_reduced: f64, log_center: f64) -> Result<f64, GasReactionError> {
    let center10 = log_center / core::f64::consts::LN_10;
    let c = -0.4 - 0.67 * center10;
    let n = 0.75 - 1.27 * center10;
    let x = log_reduced / core::f64::consts::LN_10 + c;
    let y = n - 0.14 * x;
    if ![c, n, x, y].iter().all(|v| v.is_finite()) || (x == 0.0 && y == 0.0) {
        return Err(invalid("Troe broadening coordinates are unrepresentable"));
    }
    // 1/(1+(x/y)^2) = y^2/(x^2+y^2), including the continuous y=0
    // value zero. Scaling by the larger coordinate avoids either square
    // overflowing and never divides by the vanishing denominator in f1.
    let weight = if x.abs() <= y.abs() {
        let ratio = x / y;
        1.0 / (1.0 + ratio * ratio)
    } else {
        let ratio = y / x;
        let square = ratio * ratio;
        square / (1.0 + square)
    };
    let log_lindemann = if log_reduced >= 0.0 {
        -fs_math::det::ln(1.0 + fs_math::det::exp(-log_reduced))
    } else {
        log_reduced - fs_math::det::ln(1.0 + fs_math::det::exp(log_reduced))
    };
    let log_factor = log_lindemann + log_center * weight;
    if !log_factor.is_finite() {
        return Err(invalid("Troe multiplier logarithm is unrepresentable"));
    }
    Ok(fs_math::det::exp(log_factor))
}
