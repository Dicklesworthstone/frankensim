//! A source-bound K(T) path consumed by the shared crystallization law.
//!
//! Preparation retains original observations, applicability and source receipts.
//! It neither fits rates nor converts a coefficient with units s^-n into K.
//! The supplied key must explicitly name the Nakamura rate in inverse seconds;
//! Avrami exponent, attainable crystallinity and melting limit stay explicit.

use super::{
    CrystallizationCaloricPair, CrystallizationError, CrystallizationLaw, CrystallizationRate,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, PropertyLineSlice, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};

/// The complete supported temperature line for a single original kinetic claim.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedCrystallizationRate {
    material_card_identity: ContentHash,
    identity: ContentHash,
    source: PropertyLineSlice,
}

impl PreparedCrystallizationRate {
    /// Query every source breakpoint once, fixing pressure, composition, process
    /// and all other coordinates. Gaps, source switching, unknown temperature
    /// support and negative/nonfinite rates refuse through the original owner.
    /// A positive constant is admissible only on an explicitly bounded interval.
    pub fn new(
        card: &MaterialCard,
        rate_key: &PropertyKey,
        temperature_axis: &str,
        lower: &QueryPoint,
        upper: &QueryPoint,
        selection: ClaimSelection,
    ) -> Result<Self, CrystallizationError> {
        let invalid = CrystallizationError::Invalid;
        if rate_key.quantity() != QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])) {
            return Err(invalid(
                "Nakamura K must be a rate in inverse seconds, not the Avrami s^-n coefficient",
            ));
        }
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        for point in [lower, upper] {
            if !point
                .axes()
                .get(temperature_axis)
                .is_some_and(|t| t.is_finite() && *t > 0.0)
                || !point
                    .axis_quantities()
                    .get(temperature_axis)
                    .is_some_and(|q| {
                        *q == absolute || *q == QuantitySpec::dimensional(Temperature::DIMS)
                    })
            {
                return Err(invalid(
                    "kinetic source temperatures must be positive absolute kelvin",
                ));
            }
        }
        if lower.axes()[temperature_axis] >= upper.axes()[temperature_axis] {
            return Err(invalid(
                "kinetic source temperature interval must have positive width",
            ));
        }
        let source = card
            .claims()
            .query_line_typed(rate_key, temperature_axis, lower, upper, selection)
            .map_err(|_| {
                invalid("kinetic source lacks continuous support at the supplied typed conditions")
            })?;
        let claim = card
            .claims()
            .claim(source.knots[0].1.receipt.selected)
            .ok_or(invalid("kinetic source claim is unavailable"))?;
        if !claim.validity.bounds().contains_key(temperature_axis) {
            return Err(invalid(
                "kinetic source must explicitly bound its temperature applicability",
            ));
        }
        if source
            .knots
            .iter()
            .any(|(_, a)| !a.evidence.value.value.is_finite() || a.evidence.value.value < 0.0)
        {
            return Err(invalid(
                "Nakamura source rates must be finite and nonnegative",
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-material.crystallization-rate.v1");
        hash.update(card.content_hash().as_bytes());
        for (t, answer) in &source.knots {
            hash.update(&t.to_bits().to_le_bytes());
            hash.update(answer.receipt.content_hash().as_bytes());
        }
        Ok(Self {
            material_card_identity: card.content_hash(),
            identity: hash.finalize(),
            source,
        })
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub const fn source(&self) -> &PropertyLineSlice {
        &self.source
    }
    #[must_use]
    pub fn temperature_range_k(&self) -> [f64; 2] {
        [
            self.source.knots[0].0,
            self.source.knots[self.source.knots.len() - 1].0,
        ]
    }
}

impl CrystallizationRate for PreparedCrystallizationRate {
    fn rate_s_inv(&self, temperature_k: f64) -> Result<f64, CrystallizationError> {
        let [minimum_k, maximum_k] = self.temperature_range_k();
        if !temperature_k.is_finite() || temperature_k < minimum_k || temperature_k > maximum_k {
            return Err(CrystallizationError::TemperatureOutOfRange {
                temperature_k,
                minimum_k,
                maximum_k,
            });
        }
        let knots = &self.source.knots;
        let upper = knots.partition_point(|(t, _)| *t < temperature_k);
        if knots[upper].0 == temperature_k {
            return Ok(knots[upper].1.evidence.value.value);
        }
        let (a, va) = &knots[upper - 1];
        let (b, vb) = &knots[upper];
        let weight = (temperature_k - a) / (b - a);
        let rate = (1.0 - weight) * va.evidence.value.value + weight * vb.evidence.value.value;
        if !rate.is_finite() || rate < 0.0 {
            return Err(CrystallizationError::Invalid(
                "source kinetic interpolation is not representable",
            ));
        }
        Ok(rate)
    }
}

impl CrystallizationCaloricPair<'_> {
    /// Prepare one law without caller-authored source hashes. The same material
    /// identity and complete rate support must cover both original calorics.
    /// Selected source knots and receipts become part of every kinetic state ID.
    pub fn law_with_sourced_rate<'b>(
        &'b self,
        rate: &'b PreparedCrystallizationRate,
        avrami_exponent: f64,
        maximum_crystalline_mass_fraction: f64,
        melting_temperature_k: f64,
    ) -> Result<CrystallizationLaw<'b>, CrystallizationError> {
        let [low, high] = self.temperature_range_k();
        let [rate_low, rate_high] = rate.temperature_range_k();
        if rate.material_card_identity != self.parent().curve().material_card_identity()
            || low < rate_low
            || high > rate_high
        {
            return Err(CrystallizationError::Invalid(
                "kinetic source must match the caloric material and cover its domain",
            ));
        }
        self.law(
            rate,
            rate.identity(),
            avrami_exponent,
            maximum_crystalline_mass_fraction,
            melting_temperature_k,
        )
    }
}
