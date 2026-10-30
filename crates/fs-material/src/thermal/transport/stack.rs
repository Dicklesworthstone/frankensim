//! Retain each material's conductivity source on an ordered planar stack.
//!
//! The caller declares the scalar chart aligned with each layer's heat-flow
//! direction. A source's unknown orientation stays unknown; no isotropic material,
//! matched specimen, contact law, caloric law or moisture trajectory is inferred.

/// Spatially resolved layer storage on the original conservative thermal owner.
#[path = "stack/transient.rs"]
pub mod transient;

use super::{CONDUCTIVITY, ThermalTransportError, admit_temperature_point, positive, source_line};
use crate::thermal::conduction::{
    ConductivityCurve, LayerConductionConfig, LayerStackConduction, LayerStackResponse,
    MAX_CONDUCTION_LAYERS,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, PropertyLineSlice, QueryPoint};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::DepositionError;

/// One layer's geometry and original complete fixed-context source selection.
#[derive(Clone, Copy, Debug)]
pub struct LayerStackSourceQuery<'a> {
    /// Stable layer ID; repeated material cards may have distinct layer IDs.
    pub id: u64,
    pub card: &'a MaterialCard,
    pub conductivity: &'a PropertyKey,
    pub thickness_m: f64,
    pub temperature_axis: &'a str,
    pub lower: &'a QueryPoint,
    pub upper: &'a QueryPoint,
    pub selection: ClaimSelection,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SourcedConductionLayer {
    id: u64,
    material_card: ContentHash,
    thickness: f64,
    temperature_axis: String,
    endpoints: [QueryPoint; 2],
    source: PropertyLineSlice,
}

impl SourcedConductionLayer {
    #[must_use]
    pub const fn id(&self) -> u64 {
        self.id
    }
    #[must_use]
    pub const fn material_card_identity(&self) -> ContentHash {
        self.material_card
    }
    #[must_use]
    pub const fn thickness_m(&self) -> f64 {
        self.thickness
    }
    #[must_use]
    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }
    #[must_use]
    pub const fn query_endpoints(&self) -> &[QueryPoint; 2] {
        &self.endpoints
    }
    #[must_use]
    pub const fn source_line(&self) -> &PropertyLineSlice {
        &self.source
    }
}

/// Immutable material and geometry binding for repeated multilayer flux queries.
///
/// Original source lines are retained separately; they are never resampled onto
/// a common grid or merged into an invented effective conductivity. The common
/// support is the intersection of the selected intervals. Layer order, thickness,
/// area, contact resistance and original receipts all contribute to identity.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedLayerStackConduction {
    layers: Vec<SourcedConductionLayer>,
    curves: Vec<ConductivityCurve>,
    thicknesses: Vec<f64>,
    contacts: Vec<f64>,
    area: f64,
    span: [f64; 2],
    identity: ContentHash,
}

impl PreparedLayerStackConduction {
    pub fn try_new(
        queries: &[LayerStackSourceQuery<'_>],
        contact_resistances_m2_k_w: &[f64],
        area_m2: f64,
    ) -> Result<Self, ThermalTransportError> {
        let n = queries.len();
        if n == 0
            || n > MAX_CONDUCTION_LAYERS
            || contact_resistances_m2_k_w.len() != n - 1
            || contact_resistances_m2_k_w
                .iter()
                .any(|r| !r.is_finite() || *r < 0.0)
            || !positive(area_m2)
        {
            return Err(invalid(
                "sourced stack requires 1..=64 layers, positive area and N-1 nonnegative contacts",
            ));
        }
        if queries
            .iter()
            .enumerate()
            .any(|(i, query)| queries[..i].iter().any(|q| q.id == query.id))
        {
            return Err(invalid("sourced stack layer identities must be distinct"));
        }
        let mut layers = Vec::with_capacity(n);
        let mut curves = Vec::with_capacity(n);
        let mut thicknesses = Vec::with_capacity(n);
        let mut span = [0.0_f64, f64::INFINITY];
        for query in queries {
            admit(query)?;
            let source = source_line(
                query.card,
                query.conductivity,
                query.temperature_axis,
                query.lower,
                query.upper,
                query.selection,
            )?;
            let curve = ConductivityCurve::new(
                source
                    .knots
                    .iter()
                    .map(|(t, a)| [*t, a.evidence.value.value])
                    .collect(),
            )
            .map_err(|_| {
                invalid("layer conductivity must be positive, finite and continuously supported")
            })?;
            span[0] = span[0].max(query.lower.axes()[query.temperature_axis]);
            span[1] = span[1].min(query.upper.axes()[query.temperature_axis]);
            layers.push(SourcedConductionLayer {
                id: query.id,
                material_card: query.card.content_hash(),
                thickness: query.thickness_m,
                temperature_axis: query.temperature_axis.to_owned(),
                endpoints: [query.lower.clone(), query.upper.clone()],
                source,
            });
            curves.push(curve);
            thicknesses.push(query.thickness_m);
        }
        if span[0] >= span[1] {
            return Err(invalid(
                "sourced layers need a positive-width common temperature interval",
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-material.prepared-layer-stack.v1");
        hash.update(&(n as u64).to_le_bytes());
        hash.update(&area_m2.to_bits().to_le_bytes());
        for layer in &layers {
            hash.update(&layer.id.to_le_bytes());
            hash.update(layer.material_card.as_bytes());
            hash.update(&layer.thickness.to_bits().to_le_bytes());
            hash.update(&(layer.temperature_axis.len() as u64).to_le_bytes());
            hash.update(layer.temperature_axis.as_bytes());
            hash.update(&(layer.source.knots.len() as u64).to_le_bytes());
            for (_, answer) in &layer.source.knots {
                hash.update(answer.receipt.content_hash().as_bytes());
            }
        }
        for resistance in contact_resistances_m2_k_w {
            hash.update(&resistance.to_bits().to_le_bytes());
        }
        Ok(Self {
            layers,
            curves,
            thicknesses,
            contacts: contact_resistances_m2_k_w.to_vec(),
            area: area_m2,
            span,
            identity: hash.finalize(),
        })
    }

    /// Feed the returned secant conductance to the existing thermal link owner;
    /// endpoint derivatives are separately available for a Newton assembly.
    /// Evaluation allocates nothing and neither queries nor hashes source data.
    pub fn evaluate(
        &self,
        temperatures_k: [f64; 2],
        config: LayerConductionConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<LayerStackResponse, DepositionError> {
        LayerStackConduction::new(&self.curves, &self.thicknesses, &self.contacts, self.area)?
            .evaluate(temperatures_k, config, checkpoint)
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub fn layers(&self) -> &[SourcedConductionLayer] {
        &self.layers
    }
    #[must_use]
    pub const fn temperature_span_k(&self) -> [f64; 2] {
        self.span
    }
    #[must_use]
    pub fn contact_resistances_m2_k_w(&self) -> &[f64] {
        &self.contacts
    }
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area
    }
}

fn admit(query: &LayerStackSourceQuery<'_>) -> Result<(), ThermalTransportError> {
    let static_k = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::ThermalConductivity,
        ValueForm::Static,
    ));
    if !positive(query.thickness_m)
        || (query.conductivity.quantity() != static_k
            && query.conductivity.quantity() != QuantitySpec::dimensional(CONDUCTIVITY))
        || query.conductivity.hardness_test().is_some()
        || query.conductivity.elastic_component().is_some()
        || query.conductivity.strain_component().is_some()
        || query.conductivity.stress_component().is_some()
    {
        return Err(invalid(
            "a layer needs positive thickness and scalar static thermal conductivity",
        ));
    }
    for point in [query.lower, query.upper] {
        admit_temperature_point(point, query.temperature_axis)?;
    }
    if query.lower.axis_quantities() != query.upper.axis_quantities()
        || !query.lower.axes().keys().eq(query.upper.axes().keys())
        || query.lower.axes().iter().any(|(axis, value)| {
            axis != query.temperature_axis && query.upper.axes()[axis] != *value
        })
        || query.lower.axes()[query.temperature_axis] >= query.upper.axes()[query.temperature_axis]
    {
        return Err(invalid(
            "each layer path must retain its complete fixed source context",
        ));
    }
    Ok(())
}

fn invalid(reason: &'static str) -> ThermalTransportError {
    ThermalTransportError::Invalid(reason)
}
