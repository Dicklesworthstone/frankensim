//! Shared source admission and read-only physical inventory for spatial consumers.

use super::*;

/// The same typed source contract for both slab and Cartesian finite volumes.
pub(super) fn bind_source(
    card: &MaterialCard,
    point: &QueryPoint,
    diffusivities: &[&str],
    selection: MaterialPropertySelection,
    maximum_moisture: f64,
    allow_degenerate: bool,
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
                ScalarAdmissibility::NonNegative,
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
    if !allow_degenerate {
        admit_transport_state(&transport)?;
    }
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

/// A zero value is an impermeable direction only when the source is independent
/// of evolving moisture. A degenerate D(M)=0 at a dry front is not an infinite
/// material resistance; freezing its harmonic face value would prevent wetting.
fn admit_transport_state(state: &ResolvedMaterialStatePoint) -> Result<(), MoistureError> {
    if state.properties().iter().any(|property| {
        property.value_si() == 0.0
            && property
                .source_axes()
                .iter()
                .any(|(axis, _)| axis == "moisture_content")
    }) {
        return Err(MoistureError::Invalid(
            "zero diffusivity requires a moisture-independent impermeable source",
        ));
    }
    Ok(())
}

pub(super) fn resolve_transport_state(
    card: &MaterialCard,
    binding: &ResolvedMaterialStatePoint,
    point: &QueryPoint,
) -> Result<ResolvedMaterialStatePoint, MoistureError> {
    let state = binding.at_query_point(card, point)?;
    admit_transport_state(&state)?;
    Ok(state)
}

/// A prescribed uniform-temperature event changes the transport environment,
/// never the reference density datum or either conserved inventory. Reuse the
/// originally selected claims, including on a return to an earlier temperature.
/// Bind at an occupied state: the original reference moisture need not remain
/// in the current temperature cut of the transport source's support.
pub(super) fn at_prescribed_temperature(
    card: &MaterialCard,
    binding: &ResolvedMaterialStatePoint,
    temperature_kelvin: f64,
    moisture: &[f64],
    allow_degenerate: bool,
    cx: &Cx<'_>,
) -> Result<ResolvedMaterialStatePoint, MoistureError> {
    poll(cx)?;
    if !positive(temperature_kelvin) {
        return Err(MoistureError::Invalid(
            "prescribed absolute temperature must be finite and positive",
        ));
    }
    let point = binding
        .query()
        .with("temperature", temperature_kelvin)
        .and_then(|point| point.with("moisture_content", moisture[0]))
        .map_err(|_| MoistureError::Invalid("unrepresentable prescribed temperature context"))?;
    let transport = binding.at_query_point(card, &point)?;
    if !allow_degenerate {
        admit_transport_state(&transport)?;
    }
    for &m in &moisture[1..] {
        poll(cx)?;
        let point = point
            .clone()
            .with("moisture_content", m)
            .map_err(|_| MoistureError::Invalid("unrepresentable dry-basis moisture coordinate"))?;
        let state = transport.at_query_point(card, &point)?;
        if !allow_degenerate {
            admit_transport_state(&state)?;
        }
    }
    poll(cx)?;
    Ok(transport)
}

/// Immutable dry-skeleton measure of a native finite-volume field.
/// Curved/nonuniform cells retain their individual physical inventories.
#[derive(Clone, Copy)]
pub(super) enum DryMasses<'a> {
    Uniform(f64),
    Cells(&'a [f64]),
}

/// Geometry-independent view of an accepted physical inventory.
/// The mechanical owner supplies actual cell geometry, never a surrogate slab.
#[derive(Clone, Copy)]
pub(super) struct MoistureField<'a> {
    pub card: &'a MaterialCard,
    pub reference: &'a ResolvedMaterialStatePoint,
    pub transport: &'a ResolvedMaterialStatePoint,
    pub moisture: &'a [f64],
    pub dry_masses: DryMasses<'a>,
    pub time_s: f64,
    pub identity: ContentHash,
}
impl MoistureField<'_> {
    pub fn dry_mass_at(&self, cell: usize) -> f64 {
        match self.dry_masses {
            DryMasses::Uniform(mass) => mass,
            DryMasses::Cells(masses) => masses[cell],
        }
    }
    pub fn dry_mass_kg(&self) -> f64 {
        match self.dry_masses {
            DryMasses::Uniform(mass) => mass * self.moisture.len() as f64,
            DryMasses::Cells(masses) => masses.iter().sum(),
        }
    }
    pub fn same_dry_masses(&self, other: &Self) -> bool {
        match (self.dry_masses, other.dry_masses) {
            (DryMasses::Uniform(a), DryMasses::Uniform(b)) => a.to_bits() == b.to_bits(),
            (DryMasses::Cells(a), DryMasses::Cells(b)) => {
                a.len() == b.len() && a.iter().zip(b).all(|(a, b)| a.to_bits() == b.to_bits())
            }
            _ => false,
        }
    }
    pub fn temperature_kelvin(&self) -> f64 {
        self.transport
            .state_coordinate("temperature")
            .expect("admitted absolute temperature")
    }
    pub fn point_at(&self, moisture: f64) -> Result<QueryPoint, MoistureError> {
        self.transport
            .query()
            .with("moisture_content", moisture)
            .map_err(|_| MoistureError::Invalid("unrepresentable dry-basis moisture coordinate"))
    }
    pub fn wet_mass_kg(&self) -> f64 {
        self.dry_mass_kg() + self.water_mass_kg()
    }
    pub fn water_mass_kg(&self) -> f64 {
        self.moisture
            .iter()
            .enumerate()
            .map(|(i, m)| m * self.dry_mass_at(i))
            .sum()
    }
}
