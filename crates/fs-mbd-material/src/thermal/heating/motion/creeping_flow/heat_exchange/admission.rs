//! Source/context and boundary ownership checks outside the transport hot loop.
use super::*;
use fs_material::state_point::{DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_PROPERTY};

pub(super) fn positive(value: f64) -> Result<f64, HeatedMotionError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(HeatedMotionError::Invalid(
            "sphere thermal coefficient must be positive and representable",
        ))
    }
}

pub(super) fn admit_solid<F: FluidTransport + ?Sized>(
    creeping: &CreepingExpansion<'_, '_, '_, F>,
    exchange: SphereHeatExchange<'_>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), HeatedMotionError> {
    let limits = exchange.limits;
    if [
        limits.maximum_peclet,
        limits.maximum_rayleigh,
        limits.maximum_biot,
        limits.maximum_diffusion_ratio,
        limits.maximum_relative_temperature_contrast,
    ]
    .into_iter()
    .any(|limit| !limit.is_finite() || limit <= 0.0 || limit > 0.1)
    {
        return Err(HeatedMotionError::Invalid(
            "sphere conduction limits must each be in (0,0.1]",
        ));
    }
    let source = creeping.motion.expansion.reference.regions()[0].material();
    if exchange.solid.material_card_identity() != source.card_identity() {
        return Err(HeatedMotionError::Invalid(
            "sphere conductivity and original solid material card disagree",
        ));
    }
    for endpoint in exchange.solid.query_endpoints() {
        checkpoint()?;
        if endpoint.axis_quantities() != source.axis_quantities()
            || endpoint.axes().len() != source.query_point().len()
            || source.query_point().iter().any(|(axis, value)| {
                endpoint
                    .axes()
                    .get(axis)
                    .is_none_or(|actual| axis != "T" && actual.to_bits() != value.to_bits())
            })
        {
            return Err(HeatedMotionError::Invalid(
                "sphere conductivity must retain the original solid source coordinates and non-temperature context",
            ));
        }
    }
    let region = &creeping.motion.expansion.regions[0];
    for (line, name) in exchange.solid.source_lines()[..2]
        .iter()
        .zip([DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_PROPERTY])
    {
        let selected = region
            .caloric
            .receipts()
            .iter()
            .find(|receipt| receipt.property == name)
            .ok_or(HeatedMotionError::Invalid(
                "sphere thermal transport lacks its original caloric source",
            ))?
            .selected;
        for (_, answer) in &line.knots {
            checkpoint()?;
            if answer.receipt.property != name || answer.receipt.selected != selected {
                return Err(HeatedMotionError::Invalid(
                    "sphere conductivity preparation must retain the same density and Cp sources as caloric storage",
                ));
            }
        }
    }
    Ok(())
}

pub(super) fn admit_boundary(
    thermal: &DepositionSimulation<'_>,
    id: CellId,
    reference_radius: f64,
) -> Result<(), HeatedMotionError> {
    let cell = thermal.cells().get(&id).ok_or(HeatedMotionError::Invalid(
        "inactive spherical thermal cell",
    ))?;
    let area = positive(4.0 * std::f64::consts::PI * reference_radius * reference_radius)?;
    if thermal.cells().len() != 1
        || !thermal.links().is_empty()
        || cell.boundary.htc_w_m2_k != 0.0
        || cell.boundary.reservoir.is_some()
        || (cell.geometry.surface_area_m2 / area - 1.0).abs() > 1e-10
        || (thermal.exposed_area_m2(id)? / area - 1.0).abs() > 1e-10
    {
        return Err(HeatedMotionError::Invalid(
            "spherical exterior heat needs its complete exposed area without another Robin, reservoir or contact owner",
        ));
    }
    Ok(())
}
