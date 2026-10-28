//! Local thermal transport from volumetric storage and conductivity.
//!
//! Transient-plane-source and related measurements can supply C_vol and k
//! without independent density or mass-specific heat. Keep that source basis:
//! this adapter never decomposes C_vol, mints a density, or reinterprets it as
//! thermodynamic constant-volume mass-specific heat Cv. The caller selects the
//! source's caloric constraint and scalar conductivity direction explicitly.

/// Continuously supported volumetric storage and transport paths.
#[path = "volumetric/prepared.rs"]
mod prepared;
pub use prepared::{
    PreparedVolumetricThermalTransport, VolumetricThermalPathQuery, VolumetricThermalSample,
};

use super::{ThermalTransportError, admit_temperature_point, source_line, transport_values};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialAnswer, MaterialCard, PropertyKey, QueryPoint};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

/// Volumetric heat capacity [J/(m3 K)], distinct from specific heat [J/(kg K)].
pub const VOLUMETRIC_HEAT_CAPACITY_DIMS: Dims = Dims([-1, 1, -2, -1, 0, 0]);

/// Independent source choices at one exact, fully specified condition.
#[derive(Clone, Copy, Debug)]
pub struct VolumetricThermalSourceQuery<'a> {
    pub volumetric_heat_capacity: &'a PropertyKey,
    pub conductivity: &'a PropertyKey,
    pub temperature_axis: &'a str,
    pub point: &'a QueryPoint,
    /// Volumetric heat capacity, then conductivity.
    pub selections: [ClaimSelection; 2],
}

/// Original two-property snapshot and its algebraic transport consequences.
/// Single-temperature support remains single-temperature support. This value
/// does not authorize a finite thermal excursion or supply a caloric integral.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedVolumetricThermalTransport {
    material_card: ContentHash,
    point: QueryPoint,
    temperature_axis: String,
    sources: [MaterialAnswer; 2],
    // C_vol, k, alpha, effusivity; no invented rho or mass-specific c.
    values: [f64; 4],
    identity: ContentHash,
}

impl ResolvedVolumetricThermalTransport {
    /// Resolve both exact source keys through the canonical evaluator and the
    /// same static thermal admission as PreparedThermalTransport. Independent
    /// selection does not imply independent observations or known covariance.
    pub fn try_new(
        card: &MaterialCard,
        query: VolumetricThermalSourceQuery<'_>,
    ) -> Result<Self, ThermalTransportError> {
        let keys = [query.volumetric_heat_capacity, query.conductivity];
        admit_keys(keys)?;
        admit_temperature_point(query.point, query.temperature_axis)?;
        let resolve = |index| -> Result<MaterialAnswer, ThermalTransportError> {
            let mut line = source_line(
                card,
                keys[index],
                query.temperature_axis,
                query.point,
                query.point,
                query.selections[index],
            )?;
            Ok(line
                .knots
                .pop()
                .expect("canonical point line is nonempty")
                .1)
        };
        let sources = [resolve(0)?, resolve(1)?];
        let capacity = sources[0].evidence.value.value;
        let conductivity = sources[1].evidence.value.value;
        let [diffusivity, effusivity] = transport_values(capacity, conductivity)?;
        let material_card = card.content_hash();
        let mut identity = DomainHasher::new(
            "org.frankensim.fs-material.resolved-volumetric-thermal-transport.v1",
        );
        identity.update(material_card.as_bytes());
        identity.update(&(query.temperature_axis.len() as u64).to_le_bytes());
        identity.update(query.temperature_axis.as_bytes());
        for answer in &sources {
            identity.update(answer.receipt.content_hash().as_bytes());
        }
        Ok(Self {
            material_card,
            point: query.point.clone(),
            temperature_axis: query.temperature_axis.to_owned(),
            sources,
            values: [capacity, conductivity, diffusivity, effusivity],
            identity: identity.finalize(),
        })
    }

    /// Re-query the same card and exact selected claims at a new condition.
    /// Both source quantities must support the update; failure leaves this
    /// snapshot unchanged. Coordinate names and quantity descriptors stay fixed.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, ThermalTransportError> {
        if card.content_hash() != self.material_card
            || point.axis_quantities() != self.point.axis_quantities()
            || !point.axes().keys().eq(self.point.axes().keys())
        {
            return Err(ThermalTransportError::Invalid(
                "volumetric thermal refresh requires the original card and coordinate schema",
            ));
        }
        let keys = self.sources.each_ref().map(|answer| {
            &card
                .claims()
                .claim(answer.receipt.selected)
                .expect("unchanged card retains the selected source")
                .key
        });
        Self::try_new(
            card,
            VolumetricThermalSourceQuery {
                volumetric_heat_capacity: keys[0],
                conductivity: keys[1],
                temperature_axis: &self.temperature_axis,
                point,
                selections: self
                    .sources
                    .each_ref()
                    .map(|answer| ClaimSelection::Pinned(answer.receipt.selected)),
            },
        )
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub const fn material_card_identity(&self) -> ContentHash {
        self.material_card
    }
    #[must_use]
    pub const fn query_point(&self) -> &QueryPoint {
        &self.point
    }
    /// Original source answers in capacity/conductivity order, with uncertainty.
    #[must_use]
    pub const fn sources(&self) -> &[MaterialAnswer; 2] {
        &self.sources
    }
    #[must_use]
    pub const fn volumetric_heat_capacity_j_m3_k(&self) -> f64 {
        self.values[0]
    }
    #[must_use]
    pub const fn conductivity_w_m_k(&self) -> f64 {
        self.values[1]
    }
    #[must_use]
    pub const fn diffusivity_m2_s(&self) -> f64 {
        self.values[2]
    }
    #[must_use]
    pub const fn effusivity_j_m2_k_sqrt_s(&self) -> f64 {
        self.values[3]
    }
}

fn admit_keys(keys: [&PropertyKey; 2]) -> Result<(), ThermalTransportError> {
    let thermal = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::ThermalConductivity,
        ValueForm::Static,
    ));
    if keys[0].quantity() != QuantitySpec::dimensional(VOLUMETRIC_HEAT_CAPACITY_DIMS)
        || (keys[1].quantity() != thermal
            && keys[1].quantity() != QuantitySpec::dimensional(super::CONDUCTIVITY))
        || keys.iter().any(|key| {
            key.hardness_test().is_some()
                || key.elastic_component().is_some()
                || key.strain_component().is_some()
                || key.stress_component().is_some()
        })
    {
        return Err(ThermalTransportError::Invalid(
            "volumetric storage and scalar static conductivity require their exact quantity roles",
        ));
    }
    Ok(())
}
