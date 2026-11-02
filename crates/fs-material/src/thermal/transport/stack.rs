//! Retain each material's conductivity source on an ordered planar or radial stack.
//!
//! The caller declares the scalar chart aligned with each layer's heat-flow
//! direction. A source's unknown orientation stays unknown; no isotropic material,
//! matched specimen, contact law, caloric law or moisture trajectory is inferred.

/// Spatially resolved layer storage on the original conservative thermal owner.
#[path = "stack/transient.rs"]
pub mod transient;

/// Source-bound phase layers on the original enthalpy timeline.
#[path = "stack/phase_transient.rs"]
pub mod phase_transient;

#[path = "stack/geometry.rs"]
mod geometry;
pub use geometry::LayerStackGeometry;

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
    resistance_lengths: Vec<f64>,
    contacts: Vec<f64>,
    reference_contacts: Vec<f64>,
    geometry: LayerStackGeometry,
    layer_bounds: Vec<[f64; 2]>,
    face_areas: Vec<f64>,
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
        Self::try_new_geometric(
            queries,
            contact_resistances_m2_k_w,
            LayerStackGeometry::Planar { area_m2 },
        )
    }

    /// Prepare the same source-integrated stack in fixed physical geometry.
    /// Radial legs use A_reference integral(dr/A); contacts retain their actual
    /// interface areas. The original physical thicknesses remain on `layers`.
    /// The reference area is the inner face, never a substitute storage area.
    pub fn try_new_geometric(
        queries: &[LayerStackSourceQuery<'_>],
        contact_resistances_m2_k_w: &[f64],
        geometry: LayerStackGeometry,
    ) -> Result<Self, ThermalTransportError> {
        let area_m2 = geometry.reference_area_m2().map_err(invalid)?;
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
        let mut resistance_lengths = Vec::with_capacity(n);
        let mut layer_bounds = Vec::with_capacity(n);
        let mut face_areas = Vec::with_capacity(n + 1);
        face_areas.push(area_m2);
        let mut coordinate = geometry.start_m();
        let mut span = [0.0_f64, f64::INFINITY];
        for query in queries {
            admit(query)?;
            let end = coordinate + query.thickness_m;
            if !end.is_finite() || end <= coordinate {
                return Err(invalid("layer position is not representable"));
            }
            let bounds = [coordinate, end];
            let distance = if geometry.is_planar() {
                query.thickness_m
            } else {
                geometry.resistance_length_m(bounds).map_err(invalid)?
            };
            resistance_lengths.push(distance);
            layer_bounds.push(bounds);
            face_areas.push(geometry.face_area_m2(end).map_err(invalid)?);
            coordinate = end;
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
        }
        if span[0] >= span[1] {
            return Err(invalid(
                "sourced layers need a positive-width common temperature interval",
            ));
        }
        let reference_contacts = contact_resistances_m2_k_w
            .iter()
            .enumerate()
            .map(|(i, &r)| {
                geometry
                    .reference_contact_resistance(r, layer_bounds[i][1])
                    .map_err(invalid)
            })
            .collect::<Result<Vec<_>, _>>()?;
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
        // Preserve planar identities; radial coordinates must not alias a
        // planar stack or another shell with the same reference surface area.
        match geometry {
            LayerStackGeometry::Planar { .. } => {}
            LayerStackGeometry::Cylindrical {
                inner_radius_m,
                length_m,
            } => {
                hash.update(b"cylindrical");
                hash.update(&inner_radius_m.to_bits().to_le_bytes());
                hash.update(&length_m.to_bits().to_le_bytes());
            }
            LayerStackGeometry::Spherical { inner_radius_m } => {
                hash.update(b"spherical");
                hash.update(&inner_radius_m.to_bits().to_le_bytes());
            }
        }
        Ok(Self {
            layers,
            curves,
            resistance_lengths,
            contacts: contact_resistances_m2_k_w.to_vec(),
            reference_contacts,
            geometry,
            layer_bounds,
            face_areas,
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
        LayerStackConduction::new(
            &self.curves,
            &self.resistance_lengths,
            &self.reference_contacts,
            self.area,
        )?
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
    pub const fn geometry(&self) -> LayerStackGeometry {
        self.geometry
    }
    /// Physical axial positions or radii of each original layer's two faces.
    #[must_use]
    pub fn layer_bounds_m(&self) -> &[[f64; 2]] {
        &self.layer_bounds
    }
    /// Actual areas of the N+1 ordered layer faces, including both boundaries.
    #[must_use]
    pub fn face_areas_m2(&self) -> &[f64] {
        &self.face_areas
    }
    /// Common planar area, or the radial INNER reference face area. Radial
    /// storage and boundary/contact exchange must use their actual geometry.
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
