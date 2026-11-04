//! Re-evaluate one spatial rod description from actual material-card data.
//!
//! Reference geometry, source selection and support conditions are independent.
//! Each refresh evaluates the original elastic sources at the new local
//! temperature and integrates the SAME alpha claim from the original
//! stress-free temperature. It does not add successive endpoint-alpha
//! increments or select another fit on heating.

use fs_exec::CancelGate;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer,
    LayerIsotropicParameters, MAX_CONCENTRIC_LAYERS};
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::state_point::{
    DENSITY_PROPERTY, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    YOUNG_MODULUS_PROPERTY, integrate_isotropic_thermal_expansion, resolve_material_state_point,
};
use fs_qty::{Density, Pressure};

use super::{SpatialMaterialRod, SpatialRodError, checkpoint};
use crate::rod::MAX_MATERIAL_ROD_SEGMENTS;
use crate::rod::material_field::reference::RodReferenceGeometry;

/// The same pinned rod description lowered to complete layer thermodynamics.
mod thermoelastic;

/// One spatial layer, its immutable card and complete stress-free state.
#[derive(Clone, Debug)]
pub struct RodLayerSource<'a> {
    /// Outer radius in metres; layers remain in spatial order.
    pub outer_radius_m: f64,
    /// Actual material source, not a name used to select a preset.
    pub card: &'a MaterialCard,
    /// Complete reference coordinates, including positive absolute `T` in K.
    /// Nonthermal coordinates and their quantity conventions remain unchanged.
    pub reference_point: QueryPoint,
    /// Selection for the mechanical source requirements chosen by the constructor.
    /// Explicit isotropic layers retain their own density and independent pair.
    pub mechanical_selection: MaterialPropertySelection,
    /// Selection for the canonical linear-expansion coefficient alone.
    pub expansion_selection: MaterialPropertySelection,
}

/// One equal-reference-length segment of the rod.
#[derive(Clone, Debug)]
pub struct RodSegmentSource<'a> {
    /// Explicit central bore radius [m], zero for a solid core.
    pub inner_radius_m: f64,
    /// Explicit shear correction of the existing concentric-section model.
    pub shear_correction: f64,
    /// Spatially ordered material assignments; never sorted by card identity.
    pub layers: Vec<RodLayerSource<'a>>,
}

/// Reusable material description with source claims pinned at the reference.
///
/// This is not a registry, heat solver, phase selector or dynamic history. Each
/// returned rod is a quasistatic snapshot on the SAME frozen reference geometry.
/// Current density remains a source property; closed-specimen reference mass is
/// computed once from the reference density and is never recomputed on heating.
#[derive(Clone, Debug)]
pub struct RodMaterialSource<'a> {
    segments: Vec<RodSegmentSource<'a>>,
    reference: SpatialMaterialRod,
    linear_strain_limit: f64,
    reference_temperatures_k: Vec<f64>,
    reference_mass_kg: f64,
}

enum ElasticSources<'b> {
    NamedShear(&'b str),
    Isotropic(&'b [Vec<LayerIsotropicParameters>]),
}

impl<'a> RodMaterialSource<'a> {
    /// Admit canonical density/Young's modulus and one named shear modulus,
    /// then pin every selected claim. E and G define the original rod response;
    /// no additional three-dimensional isotropy assumption is introduced.
    /// The existing section and rod owners check geometry and static support.
    ///
    /// # Errors
    /// Invalid counts/geometry, missing or incompatible properties, ambiguous
    /// selections, unsupported reference conditions and cancellation refuse.
    pub fn try_new(
        length_m: f64,
        segments: Vec<RodSegmentSource<'a>>,
        shear_property: &str,
        linear_strain_limit: f64,
        gate: &CancelGate,
    ) -> Result<Self, SpatialRodError> {
        Self::build(length_m, segments, ElasticSources::NamedShear(shear_property),
            linear_strain_limit, gate)
    }

    /// Admit one explicit density source and isotropic elastic pair per layer.
    /// Parameter groups follow segment/layer order exactly; source names and
    /// E/nu, E/G or K/G choices may differ between layers. The existing material
    /// projector derives E/G without adding fictitious claims to the cards.
    /// Heating preserves these roles, original claims and reference inventory.
    ///
    /// # Errors
    /// Refuses incomplete parameter groups, unsupported or inadmissible source
    /// data, invalid geometry, reference conditions or cancellation.
    pub fn try_new_with_isotropic_parameters(
        length_m: f64,
        segments: Vec<RodSegmentSource<'a>>,
        parameters: &[Vec<LayerIsotropicParameters>],
        linear_strain_limit: f64,
        gate: &CancelGate,
    ) -> Result<Self, SpatialRodError> {
        Self::build(length_m, segments, ElasticSources::Isotropic(parameters),
            linear_strain_limit, gate)
    }

    fn build(
        length_m: f64,
        mut segments: Vec<RodSegmentSource<'a>>,
        sources: ElasticSources<'_>,
        linear_strain_limit: f64,
        gate: &CancelGate,
    ) -> Result<Self, SpatialRodError> {
        checkpoint(gate)?;
        if !length_m.is_finite() || length_m <= 0.0
            || segments.is_empty() || segments.len() > MAX_MATERIAL_ROD_SEGMENTS
            || segments.iter().any(|s| s.layers.is_empty() || s.layers.len() > MAX_CONCENTRIC_LAYERS)
            || !linear_strain_limit.is_finite() || linear_strain_limit <= 0.0 || linear_strain_limit >= 1.0
        {
            return Err(SpatialRodError::Invalid("bounded nonempty rod/layers and a linear strain limit in (0,1) are required"));
        }
        match &sources {
            ElasticSources::NamedShear(shear_property) => {
                if *shear_property == DENSITY_PROPERTY || *shear_property == YOUNG_MODULUS_PROPERTY
                    || *shear_property == LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY
                {
                    return Err(SpatialRodError::Invalid("shear modulus must have its own property key"));
                }
            }
            ElasticSources::Isotropic(parameters) => {
                if parameters.len() != segments.len() || parameters.iter().zip(&segments)
                    .any(|(bindings, segment)| bindings.len() != segment.layers.len())
                {
                    return Err(SpatialRodError::Invalid("one isotropic parameterization per segment/layer is required"));
                }
            }
        }
        let mut reference_sections = Vec::with_capacity(segments.len());
        let mut reference_temperatures_k = Vec::new();
        for (segment_index, segment) in segments.iter_mut().enumerate() {
            let mut layers = Vec::with_capacity(segment.layers.len());
            let mut expansions = Vec::with_capacity(segment.layers.len());
            for (layer_index, layer) in segment.layers.iter_mut().enumerate() {
                checkpoint(gate)?;
                let requirements = match &sources {
                    ElasticSources::NamedShear(shear_property) => [
                        (DENSITY_PROPERTY, Density::DIMS),
                        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS),
                        (*shear_property, Pressure::DIMS),
                    ].into_iter().map(|(name, dims)| {
                        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                    }).collect::<Result<Vec<_>, _>>()?,
                    ElasticSources::Isotropic(parameters) => {
                        let binding = &parameters[segment_index][layer_index];
                        std::iter::once(&binding.density).chain(binding.parameters.requirements())
                            .cloned().collect()
                    }
                };
                let material = resolve_material_state_point(layer.card, &layer.reference_point,
                    &requirements, layer.mechanical_selection.clone())?;
                let expansion = integrate_isotropic_thermal_expansion(layer.card,
                    &layer.reference_point, &layer.reference_point, layer.expansion_selection.clone())?;
                let temperature = *layer.reference_point.axes().get("T")
                    .ok_or(SpatialRodError::Invalid("reference state needs positive absolute T"))?;
                reference_temperatures_k.push(temperature);
                layer.mechanical_selection = MaterialPropertySelection::PinnedByProperty(
                    material.properties().iter().map(|p| (p.requirement().name().to_owned(),
                        p.answer().receipt.selected)).collect());
                layer.expansion_selection = MaterialPropertySelection::PinnedByProperty(vec![(
                    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY.to_owned(), expansion.selected_claim(),
                )]);
                layers.push(ConcentricMaterialLayer { outer_radius_m: layer.outer_radius_m, material });
                expansions.push(expansion);
            }
            let elastic = match &sources {
                ElasticSources::NamedShear(_) =>
                    ConcentricElasticSection::try_new(segment.inner_radius_m, layers)?,
                ElasticSources::Isotropic(parameters) =>
                    ConcentricElasticSection::try_new_with_isotropic_parameters(
                        segment.inner_radius_m, layers, &parameters[segment_index])?,
            };
            // Admit the common stress-free reference before accepting the source.
            ConcentricThermoelasticSection::try_new(elastic.clone(), expansions)?;
            reference_sections.push(match &sources {
                ElasticSources::NamedShear(shear_property) =>
                    elastic.with_shear_modulus(shear_property, segment.shear_correction)?,
                ElasticSources::Isotropic(_) => elastic.with_isotropic_shear(segment.shear_correction)?,
            });
        }
        let reference = SpatialMaterialRod::try_straight(length_m, reference_sections)?;
        let reference_mass_kg = reference.reference_mass().total_mass_kg();
        if !reference_mass_kg.is_finite() || reference_mass_kg <= 0.0 {
            return Err(SpatialRodError::Invalid("reference rod mass is unrepresentable"));
        }
        checkpoint(gate)?;
        Ok(Self { segments, reference,
            linear_strain_limit, reference_temperatures_k, reference_mass_kg })
    }

    /// Refresh ALL material coefficients at local temperatures in segment/layer
    /// order. Thermal paths always start at the retained reference, not the last
    /// observation. The first unsupported layer refuses the whole candidate.
    ///
    /// Supports and loads are applied to the returned rod through existing APIs.
    /// Its initial pose is the retained reference; no branch tracking or dynamic remap is implied.
    pub fn resolve_at(
        &self, temperatures_k: &[f64], gate: &CancelGate,
    ) -> Result<SpatialMaterialRod, SpatialRodError> {
        checkpoint(gate)?;
        if temperatures_k.len() != self.reference_temperatures_k.len()
            || temperatures_k.iter().any(|t| !t.is_finite() || *t <= 0.0)
        {
            return Err(SpatialRodError::Invalid("one finite positive temperature per source layer is required"));
        }
        let mut mechanical = Vec::with_capacity(self.segments.len());
        let mut thermal = Vec::with_capacity(self.segments.len());
        let mut temperatures = temperatures_k.iter();
        for (segment, reference) in self.segments.iter().zip(self.reference.sections()) {
            let mut layers = Vec::with_capacity(segment.layers.len());
            let mut expansions = Vec::with_capacity(segment.layers.len());
            for (layer, original) in segment.layers.iter().zip(reference.elastic().layers()) {
                checkpoint(gate)?;
                let temperature = *temperatures.next().expect("complete temperature count admitted");
                let current = point_at_temperature(&layer.reference_point, temperature)?;
                let material = original.material.at_query_point(layer.card, &current)?;
                expansions.push(integrate_isotropic_thermal_expansion(layer.card,
                    &layer.reference_point, &current, layer.expansion_selection.clone())?);
                layers.push(ConcentricMaterialLayer { outer_radius_m: layer.outer_radius_m, material });
            }
            let elastic = reference.elastic().with_material_layers(layers)?;
            thermal.push(ConcentricThermoelasticSection::try_new(elastic.clone(), expansions)?);
            mechanical.push(reference.with_material_section(elastic)?);
        }
        // Retain the original per-segment inventory and its density receipts.
        // Current rho(T) must not redefine weight on frozen reference geometry.
        let mut rod = SpatialMaterialRod::try_straight_with_mass(
            self.reference.reference_length_m(), mechanical, Some(self.reference.reference_mass.clone()),
        )?;
        rod.reference_geometry = self.reference.reference_geometry.clone();
        rod.rod.positions.clone_from(&self.reference.rod.positions);
        rod.rod.frames.clone_from(&self.reference.rod.frames);
        rod.with_thermal_sections(thermal, self.linear_strain_limit, gate)
    }

    /// Bind an unloaded curved shape to this already admitted source description.
    /// Material-coordinate length and section count must match. Source pins,
    /// reference mass and temperature paths are unchanged; heating cannot erase
    /// intrinsic curvature or mistake it for mechanical strain.
    /// This makes a new static description, not a live rest-shape update.
    pub fn with_reference_geometry(
        &self, geometry: RodReferenceGeometry, gate: &CancelGate,
    ) -> Result<Self, SpatialRodError> {
        checkpoint(gate)?;
        if geometry.length_m().to_bits() != self.reference.reference_length_m().to_bits() {
            return Err(SpatialRodError::Invalid("source and reference geometry must have the same material-coordinate length"));
        }
        let mut candidate = self.clone();
        candidate.reference = SpatialMaterialRod::try_from_reference(
            geometry, self.reference.sections.clone())?;
        candidate.reference.reference_mass = self.reference.reference_mass.clone();
        checkpoint(gate)?;
        Ok(candidate)
    }

    /// Original reference material bundles and undeformed geometry.
    #[must_use]
    pub const fn reference_rod(&self) -> &SpatialMaterialRod { &self.reference }
    /// Pinned source descriptions in segment/layer order; cards remain borrowed.
    #[must_use]
    pub fn segments(&self) -> &[RodSegmentSource<'a>] { &self.segments }
    /// Flattened reference temperatures in the order expected by `resolve_at`.
    #[must_use]
    pub fn reference_temperatures_k(&self) -> &[f64] { &self.reference_temperatures_k }
    /// Fixed reference mass [kg], not current density times unchanged volume.
    #[must_use]
    pub const fn reference_mass_kg(&self) -> f64 { self.reference_mass_kg }
    /// Original small-strain envelope, shared by force and caloric consumers.
    #[must_use]
    pub const fn linear_strain_limit(&self) -> f64 { self.linear_strain_limit }
}

fn point_at_temperature(reference: &QueryPoint, temperature: f64) -> Result<QueryPoint, SpatialRodError> {
    let mut point = QueryPoint::new();
    for (axis, &value) in reference.axes() {
        let value = if axis == "T" { temperature } else { value };
        point = match reference.axis_quantities().get(axis) {
            Some(quantity) => point.with_quantity(axis, *quantity, value),
            None => point.with(axis, value),
        }.map_err(|_| SpatialRodError::Invalid("could not retain source state coordinates"))?;
    }
    Ok(point)
}
