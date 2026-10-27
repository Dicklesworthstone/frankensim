//! Re-evaluate one spatial rod description from actual material-card data.
//!
//! Reference geometry, source selection and support conditions are independent.
//! Each refresh evaluates E/rho/G at the new local temperature and integrates
//! the SAME alpha claim from the original stress-free temperature. It does not
//! add successive endpoint-alpha increments or select another fit on heating.

use fs_exec::CancelGate;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer, MAX_CONCENTRIC_LAYERS};
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
    /// Selection for exactly density, Young's modulus and the named shear modulus.
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
    requirements: Vec<ScalarPropertyRequirement>,
    shear_property: String,
    linear_strain_limit: f64,
    reference_temperatures_k: Vec<f64>,
    reference_mass_kg: f64,
}

impl<'a> RodMaterialSource<'a> {
    /// Admit geometry and reference sources, then pin every selected claim.
    /// The existing section and rod owners check geometry and static support.
    ///
    /// # Errors
    /// Invalid counts/geometry, missing or incompatible properties, ambiguous
    /// selections, unsupported reference conditions and cancellation refuse.
    pub fn try_new(
        length_m: f64,
        mut segments: Vec<RodSegmentSource<'a>>,
        shear_property: &str,
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
        if shear_property == DENSITY_PROPERTY || shear_property == YOUNG_MODULUS_PROPERTY
            || shear_property == LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY
        {
            return Err(SpatialRodError::Invalid("shear modulus must have its own property key"));
        }
        let requirements = [
            (DENSITY_PROPERTY, Density::DIMS),
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS),
            (shear_property, Pressure::DIMS),
        ].into_iter().map(|(name, dims)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
        }).collect::<Result<Vec<_>, _>>()?;
        let mut reference_sections = Vec::with_capacity(segments.len());
        let mut reference_temperatures_k = Vec::new();
        for segment in &mut segments {
            let mut layers = Vec::with_capacity(segment.layers.len());
            let mut expansions = Vec::with_capacity(segment.layers.len());
            for layer in &mut segment.layers {
                checkpoint(gate)?;
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
            let elastic = ConcentricElasticSection::try_new(segment.inner_radius_m, layers)?;
            // Admit the common stress-free reference before accepting the source.
            ConcentricThermoelasticSection::try_new(elastic.clone(), expansions)?;
            reference_sections.push(elastic.with_shear_modulus(shear_property, segment.shear_correction)?);
        }
        let reference = SpatialMaterialRod::try_straight(length_m, reference_sections)?;
        let reference_mass_kg = reference.reference_mass().total_mass_kg();
        if !reference_mass_kg.is_finite() || reference_mass_kg <= 0.0 {
            return Err(SpatialRodError::Invalid("reference rod mass is unrepresentable"));
        }
        checkpoint(gate)?;
        Ok(Self { segments, reference, requirements, shear_property: shear_property.into(),
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
        for segment in &self.segments {
            let mut layers = Vec::with_capacity(segment.layers.len());
            let mut expansions = Vec::with_capacity(segment.layers.len());
            for layer in &segment.layers {
                checkpoint(gate)?;
                let temperature = *temperatures.next().expect("complete temperature count admitted");
                let current = point_at_temperature(&layer.reference_point, temperature)?;
                let material = resolve_material_state_point(layer.card, &current,
                    &self.requirements, layer.mechanical_selection.clone())?;
                expansions.push(integrate_isotropic_thermal_expansion(layer.card,
                    &layer.reference_point, &current, layer.expansion_selection.clone())?);
                layers.push(ConcentricMaterialLayer { outer_radius_m: layer.outer_radius_m, material });
            }
            let elastic = ConcentricElasticSection::try_new(segment.inner_radius_m, layers)?;
            thermal.push(ConcentricThermoelasticSection::try_new(elastic.clone(), expansions)?);
            mechanical.push(elastic.with_shear_modulus(&self.shear_property, segment.shear_correction)?);
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
