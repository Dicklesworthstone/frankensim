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
//! conversion. This is Lindemann (broadening factor one), not Troe/SRI falloff,
//! chemically activated kinetics or validated combustion-mechanism data.
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
use crate::cavity_phs::piston::reactive::mechanism::Sum;
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_qty::Dims;
use fs_thermochem::SpeciesId;
use fs_thermochem::equilibrium::MAX_REACTION_SPECIES_V1;
use std::collections::BTreeSet;

/// Explicit positive normalization c0 in mol/m3.
pub const REFERENCE_CONCENTRATION_PROPERTY: &str = "reaction_collider_reference_concentration";
/// Normalized low-pressure activity scale k0 in mol/(m3 s), used by Lindemann.
/// This can be a source scalar, curve or surface at the actual state point.
pub const LOW_PRESSURE_SCALE_PROPERTY: &str = "reaction_low_pressure_rate_scale";

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
        if self.regime == ColliderRegime::Lindemann {
            required.push(ScalarPropertyRequirement::try_new(
                LOW_PRESSURE_SCALE_PROPERTY,
                ACTIVITY_RATE_SCALE_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?);
        }
        let mut names = BTreeSet::new();
        for binding in &self.efficiencies {
            if binding.property == REFERENCE_CONCENTRATION_PROPERTY
                || binding.property == LOW_PRESSURE_SCALE_PROPERTY
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
        let c0 = property(REFERENCE_CONCENTRATION_PROPERTY)?;
        let model = gas.model().thermodynamics();
        let mut weighted = Sum::default();
        for (species, fraction) in model
            .components()
            .iter()
            .zip(model.mole_composition().fractions())
        {
            poll(cancelled)?;
            let index = self
                .efficiencies
                .binary_search_by(|b| b.species.cmp(species.species()))
                .map_err(|_| invalid("live gas has no declared collider efficiency"))?;
            let efficiency = property(&self.efficiencies[index].property)?;
            let contribution = efficiency * fraction;
            if !contribution.is_finite()
                || (efficiency > 0.0 && *fraction > 0.0 && contribution == 0.0)
            {
                return Err(invalid("positive collider contribution is unrepresentable"));
            }
            weighted.add(contribution)?;
        }
        // A real absence of effective colliders is distinct from numerical
        // underflow; the native kinetic owner preserves actual activities.
        if weighted.value == 0.0 {
            poll(cancelled)?;
            return Ok(0.0);
        }
        let molar_mass = gas
            .state()
            .thermodynamics()
            .properties()
            .molar_mass()
            .value();
        let concentration = (gas.mass_kg() / molar_mass) / gas.volume_m3();
        if !positive(concentration) || !positive(weighted.value) || !positive(c0) {
            return Err(invalid(
                "positive represented collider concentration required",
            ));
        }
        let log_q = fs_math::det::ln(weighted.value) + fs_math::det::ln(concentration)
            - fs_math::det::ln(c0);
        if !log_q.is_finite() {
            return Err(invalid(
                "collider concentration ratio logarithm is unrepresentable",
            ));
        }
        let factor = match self.regime {
            ColliderRegime::ThirdBody => fs_math::det::exp(log_q),
            ColliderRegime::Lindemann => {
                let low = property(LOW_PRESSURE_SCALE_PROPERTY)?;
                if !positive(base_scale) || !positive(low) {
                    return Err(invalid(
                        "positive low/high-pressure kinetic scales required",
                    ));
                }
                let log_reduced = fs_math::det::ln(low) - fs_math::det::ln(base_scale) + log_q;
                if !log_reduced.is_finite() {
                    return Err(invalid("Lindemann reduced pressure is unrepresentable"));
                }
                // Avoid overflow in Pr or 1+Pr at a perfectly finite high-
                // pressure limit. A lost positive low-pressure rate still refuses.
                if log_reduced >= 0.0 {
                    1.0 / (1.0 + fs_math::det::exp(-log_reduced))
                } else {
                    let reduced = fs_math::det::exp(log_reduced);
                    reduced / (1.0 + reduced)
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
}
