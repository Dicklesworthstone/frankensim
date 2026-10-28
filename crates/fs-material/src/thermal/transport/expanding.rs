//! Local thermal transport composed with constant-mass free-expansion density.
//!
//! The retained expansion model owns geometry, reference normalization and its
//! physical assumptions. Heat capacity and directional scalar conductivity stay
//! original selected source models. No density-versus-temperature measurement
//! is manufactured, and derived density or transport products never interpolate.
//! The caller selects a compatible caloric constraint and conductivity direction;
//! this composition does not add enthalpy, mechanical work, phase change, mass
//! evolution or uncertainty propagation.

use super::{
    ThermalTransportError, ThermalTransportSample, checked_sample,
    expansion::PreparedFreeExpansionDensity, positive, source_line, source_value,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, PropertyLineSlice, QueryPoint};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

/// Heat-capacity/conductivity grid only; density retains its separately bounded grid.
pub const MAX_EXPANDING_THERMAL_KNOTS: usize = 2 * fs_matdb::MAX_PROPERTY_LINE_KNOTS;

/// Selected thermal sources on exactly the existing density model's source path.
#[derive(Clone, Copy, Debug)]
pub struct ExpandingThermalSourceQuery<'a> {
    /// Specific heat for the declared caloric constraint [J/(kg K)].
    pub specific_heat: &'a PropertyKey,
    /// Static scalar conductivity in the caller-selected heat-flow direction [W/(m K)].
    pub conductivity: &'a PropertyKey,
    pub temperature_axis: &'a str,
    pub lower: &'a QueryPoint,
    pub upper: &'a QueryPoint,
    /// Specific heat, then conductivity; source independence is not inferred.
    pub selections: [ClaimSelection; 2],
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct ThermalKnot {
    temperature: f64,
    values: [f64; 2],
}

/// Immutable composition retaining all original density and thermal evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedExpandingThermalTransport {
    density: PreparedFreeExpansionDensity,
    sources: [PropertyLineSlice; 2],
    knots: Vec<ThermalKnot>,
    identity: ContentHash,
}

impl PreparedExpandingThermalTransport {
    /// Compose only the original material card, exact typed interval and fixed
    /// query conditions. Each thermal source must support the whole density path
    /// under its explicit selection; no trimming, extrapolation or claim stitching.
    /// Original uncertainty and dependency identities survive without covariance
    /// assumptions. Re-prepare when any source, card or fixed condition changes.
    pub fn try_new(
        card: &MaterialCard,
        density: PreparedFreeExpansionDensity,
        query: ExpandingThermalSourceQuery<'_>,
    ) -> Result<Self, ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        if density.material_card_identity() != card.content_hash()
            || density.temperature_axis() != query.temperature_axis
            || density.query_endpoints() != [query.lower, query.upper]
        {
            return Err(invalid(
                "expanding thermal transport requires the density model's exact card, axis and query endpoints",
            ));
        }
        let keys = [query.specific_heat, query.conductivity];
        let conductivity = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        ));
        if keys[0].quantity() != QuantitySpec::dimensional(super::SPECIFIC_HEAT)
            || (keys[1].quantity() != conductivity
                && keys[1].quantity() != QuantitySpec::dimensional(super::CONDUCTIVITY))
            || keys.iter().any(|key| {
                key.hardness_test().is_some()
                    || key.elastic_component().is_some()
                    || key.strain_component().is_some()
                    || key.stress_component().is_some()
            })
        {
            return Err(invalid(
                "expanding thermal transport requires mass-specific heat and scalar static conductivity",
            ));
        }
        let resolve = |index| {
            let line = source_line(
                card,
                keys[index],
                query.temperature_axis,
                query.lower,
                query.upper,
                query.selections[index],
            )?;
            if line
                .knots
                .iter()
                .any(|(_, answer)| !positive(answer.evidence.value.value))
            {
                return Err(invalid(
                    "specific heat and conductivity must remain finite and positive",
                ));
            }
            Ok(line)
        };
        let sources = [resolve(0)?, resolve(1)?];
        let mut temperatures: Vec<_> = sources
            .iter()
            .flat_map(|source| source.knots.iter().map(|(t, _)| *t))
            .collect();
        temperatures.sort_by(f64::total_cmp);
        temperatures.dedup();
        if temperatures.len() > MAX_EXPANDING_THERMAL_KNOTS {
            return Err(invalid(
                "expanding thermal source union exceeds its breakpoint budget",
            ));
        }
        let mut identity =
            DomainHasher::new("org.frankensim.fs-material.prepared-expanding-thermal-transport.v1");
        identity.update(density.identity().as_bytes());
        for source in &sources {
            identity.update(&(source.knots.len() as u64).to_le_bytes());
            for (_, answer) in &source.knots {
                identity.update(answer.receipt.content_hash().as_bytes());
            }
        }
        let knots = temperatures
            .into_iter()
            .map(|temperature| ThermalKnot {
                temperature,
                values: std::array::from_fn(|i| source_value(&sources[i], temperature)),
            })
            .collect();
        let prepared = Self {
            density,
            sources,
            knots,
            identity: identity.finalize(),
        };
        // Check every input breakpoint, including density-only knots. Keep the
        // original density model intact rather than resampling its derived rho.
        for t in prepared
            .temperatures_k()
            .chain(prepared.density.temperatures_k())
        {
            prepared.at_temperature(t)?;
        }
        Ok(prepared)
    }

    /// Evaluate the original expansion model and the two affine thermal inputs.
    /// At most two binary searches; no allocation, hashing, source query or
    /// mutation. Density is never replaced by a linear endpoint approximation.
    /// Preparation checks all source knots; extreme interior products or ratios
    /// can still be unrepresentable and refuse without changing this model.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ThermalTransportSample, ThermalTransportError> {
        let rho = self.density.at_temperature(temperature_k)?.density_kg_m3();
        let index = self
            .knots
            .partition_point(|knot| knot.temperature < temperature_k);
        let upper = self.knots[index];
        let values = if upper.temperature == temperature_k {
            upper.values
        } else {
            let lower = self.knots[index - 1];
            let weight =
                (temperature_k - lower.temperature) / (upper.temperature - lower.temperature);
            std::array::from_fn(|i| lower.values[i] + weight * (upper.values[i] - lower.values[i]))
        };
        checked_sample(self.identity, temperature_k, [rho, values[0], values[1]])
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub const fn material_card_identity(&self) -> ContentHash {
        self.density.material_card_identity()
    }
    /// Original density geometry, reference answer, source receipts and assumptions.
    #[must_use]
    pub const fn density_model(&self) -> &PreparedFreeExpansionDensity {
        &self.density
    }
    /// Original specific-heat and conductivity source lines, in that order.
    #[must_use]
    pub const fn source_lines(&self) -> &[PropertyLineSlice; 2] {
        &self.sources
    }
    #[must_use]
    pub fn temperature_axis(&self) -> &str {
        self.density.temperature_axis()
    }
    #[must_use]
    pub fn query_endpoints(&self) -> [&QueryPoint; 2] {
        self.density.query_endpoints()
    }
    #[must_use]
    pub fn temperature_span_k(&self) -> [f64; 2] {
        self.density.temperature_span_k()
    }
    /// Cp/k breakpoint union only; the retained density model has its own grid.
    pub fn temperatures_k(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.knots.iter().map(|knot| knot.temperature)
    }
}
