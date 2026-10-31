//! Shared source admission and read-only physical inventory for spatial consumers.

use super::*;

/// The same typed source contract for both slab and Cartesian finite volumes.
pub(super) fn bind_source(
    card: &MaterialCard,
    point: &QueryPoint,
    diffusivities: &[&str],
    selection: MaterialPropertySelection,
    maximum_moisture: f64,
) -> Result<(ResolvedMaterialStatePoint, ResolvedMaterialStatePoint, f64), MoistureError> {
    if point.axis_quantities().get("temperature")
        != Some(&quantity(QuantityKind::AbsoluteTemperature))
        || point.axis_quantities().get("moisture_content")
            != Some(&quantity(QuantityKind::Moisture(MoistureBasis::DryMass)))
    {
        return Err(MoistureError::Invalid(
            "absolute temperature and dry-mass moisture coordinates required",
        ));
    }
    let temperature = point.axes().get("temperature").copied().unwrap_or(f64::NAN);
    let reference_content = point
        .axes()
        .get("moisture_content")
        .copied()
        .unwrap_or(f64::NAN);
    if point.axes().contains_key("T")
        || point.axis_quantities().iter().any(|(axis, q)| {
            (axis != "temperature" && *q == quantity(QuantityKind::AbsoluteTemperature))
                || (axis != "moisture_content"
                    && *q == quantity(QuantityKind::Moisture(MoistureBasis::DryMass)))
        })
    {
        return Err(MoistureError::Invalid(
            "duplicate physical temperature/moisture coordinates are unsupported",
        ));
    }
    if !positive(temperature)
        || !reference_content.is_finite()
        || reference_content < 0.0
        || reference_content > maximum_moisture
        || diffusivities.is_empty()
        || diffusivities.len() > 3
        || diffusivities.contains(&"density")
    {
        return Err(MoistureError::Invalid(
            "positive temperature, supported reference moisture and distinct density/transport roles required",
        ));
    }
    let mut requirements = Vec::new();
    for (i, &name) in diffusivities.iter().enumerate() {
        if !diffusivities[..i].contains(&name) {
            requirements.push(ScalarPropertyRequirement::try_new(
                name,
                MOISTURE_DIFFUSIVITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )?);
        }
    }
    let mut all = vec![ScalarPropertyRequirement::try_new(
        "density",
        Density::DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )?];
    all.extend(requirements.iter().cloned());
    let reference = resolve_material_state_point(card, point, &all, selection)?;
    if reference
        .properties()
        .iter()
        .any(|p| !p.is_static_at_fixed_external_conditions())
    {
        return Err(MoistureError::Invalid(
            "transport and reference density need static environmental sources",
        ));
    }
    let pins = requirements
        .iter()
        .map(|r| {
            (
                r.name().to_owned(),
                reference
                    .property(r.name())
                    .expect("required diffusivity")
                    .answer()
                    .receipt
                    .selected,
            )
        })
        .collect();
    let transport = resolve_material_state_point(
        card,
        point,
        &requirements,
        MaterialPropertySelection::PinnedByProperty(pins),
    )?;
    let dry_density = reference
        .property("density")
        .expect("required density")
        .value_si()
        / (1.0 + reference_content);
    if !positive(dry_density) {
        return Err(MoistureError::Invalid(
            "dry reference density is unrepresentable",
        ));
    }
    Ok((reference, transport, dry_density))
}

/// Geometry-independent view of an accepted equal-reference-volume inventory.
/// The mechanical owner supplies actual cell geometry, never a surrogate slab.
#[derive(Clone, Copy)]
pub(super) struct MoistureField<'a> {
    pub card: &'a MaterialCard,
    pub reference: &'a ResolvedMaterialStatePoint,
    pub moisture: &'a [f64],
    pub dry_mass_per_cell_kg: f64,
    pub time_s: f64,
    pub identity: ContentHash,
}
impl MoistureField<'_> {
    pub fn point_at(&self, moisture: f64) -> Result<QueryPoint, MoistureError> {
        self.reference
            .query()
            .with("moisture_content", moisture)
            .map_err(|_| MoistureError::Invalid("unrepresentable dry-basis moisture coordinate"))
    }
    pub fn wet_mass_kg(&self) -> f64 {
        self.dry_mass_per_cell_kg * self.moisture.len() as f64
            + self
                .moisture
                .iter()
                .map(|m| m * self.dry_mass_per_cell_kg)
                .sum::<f64>()
    }
}
