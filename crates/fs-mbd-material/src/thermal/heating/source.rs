//! Validate storage and mechanics against the same reference material context.
use super::*;
use crate::ResolvedRegion;
use fs_material::state_point::{DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_PROPERTY};

pub(super) fn admit(original: &ResolvedRegion, region: &ThermalRegion<'_>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), HeatingError> {
    let source = original.material();
    if region.card.content_hash() != source.card_identity()
        || region.card.id() != source.material()
        || region.caloric.curve().material_card_identity() != source.card_identity()
        || region.geometry.volume_m3.to_bits() != original.geometry().volume_m3().to_bits() {
        return Err(HeatingError::Binding("caloric card, reference material and volume must agree"));
    }
    let density_claim = source.property(DENSITY_PROPERTY)
        .expect("reference body resolved density").answer().receipt.selected;
    for receipt in region.caloric.receipts() {
        checkpoint()?;
        if !matches!(receipt.property.as_str(), DENSITY_PROPERTY | SPECIFIC_HEAT_CAPACITY_PROPERTY)
            || (receipt.property == DENSITY_PROPERTY && receipt.selected != density_claim)
            || &receipt.axis_quantities != source.axis_quantities()
            || receipt.query_point.len() != source.query_point().len()
            || !receipt.query_point.iter().zip(source.query_point()).all(|((a,x),(b,y))| {
                a == b && (a == "T" || x.to_bits() == y.to_bits())
            }) {
            return Err(HeatingError::Binding("caloric claims must retain the original density source and non-temperature context"));
        }
    }
    let state = region.caloric.curve()
        .state_at_specific_enthalpy(region.initial_specific_enthalpy_j_kg)?;
    let reference_t = source.query_point().iter().find(|(axis,_)| axis == "T")
        .map(|(_,t)| *t).ok_or(HeatingError::Binding("reference temperature is required"))?;
    // Only roundoff from the existing h/T chart is admitted, not a thermal
    // mismatch or a measurement uncertainty silently converted to tolerance.
    let tolerance = 64.0 * f64::EPSILON * reference_t.abs().max(state.temperature_k().abs());
    if !state.is_fully_solid() || (state.temperature_k() - reference_t).abs() > tolerance {
        return Err(HeatingError::Binding("initial caloric state must be solid at the reference temperature"));
    }
    Ok(())
}
