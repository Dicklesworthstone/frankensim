//! Chain the native normalized-resistance adjoint through the selected card.
//! Calibration coordinates are physical pressure/thickness/gap/torque, NOT
//! resistance multipliers. No additional thermal solve or probe is performed.
use super::*;
use fs_matdb::{ClaimId, ClaimSelection, SelectionPolicy};

pub(super) fn convert(model: &Model, sample: &Sample, report: &J, coefficients: &mut [f64]) -> Result<()> {
    let parameters = model.bound.study().parameters();
    if !parameters.iter().any(|p| p.target.contact_axis().is_some() && p.low != p.high) { return Ok(()); }
    if coefficients.len() != parameters.len()
        || !matches!(report.str_field("output"), Some(output::CONTACT | output::COMBINED)) {
        return Err(invalid("joint-state calibration requires the complete contact adjoint report"));
    }
    // The marginal means generally differ from the base project. Reconstruct
    // exactly the calibration point, never query a card at the old base state.
    let project = model.project(&sample.parameters)?;
    let library = model.cards.library();
    let requirements = fs_project::BindingRequirements::thermal_steady_v1();
    let resolution = fs_project::resolve_bindings(&project.spec, &library, &requirements);
    if !resolution.admissible() { return Err(invalid("joint-state calibration cannot reproduce its card bindings")); }
    let rows = array(report,"parameters",256)?;
    for (i, parameter) in parameters.iter().enumerate() {
        let Some(axis) = parameter.target.contact_axis() else { continue; };
        if parameter.low == parameter.high { continue; }
        let binding = project.spec.interface_cards.as_deref().unwrap_or(&[]).iter()
            .find(|b| b.interface == parameter.entity).ok_or_else(|| invalid("missing sampled interface"))?;
        let resolved = resolution.bindings.iter().find(|b|
            matches!(&b.target, fs_project::BindingTarget::Interface(name) if name == &parameter.entity))
            .ok_or_else(|| invalid("missing resolved calibration interface"))?;
        let property = resolved.properties.iter().find(|p|p.property==fs_project::CONTACT_RESISTANCE_PROPERTY)
            .ok_or_else(|| invalid("missing calibrated resistance property"))?;
        let card = library.interface(&binding.card).ok_or_else(|| invalid("missing calibrated interface card"))?;
        let selection = match binding.claim.as_deref() {
            Some(pin) => ClaimSelection::Pinned(ClaimId(ContentHash::from_hex(pin)
                .ok_or_else(|| invalid("invalid interface claim pin"))?)),
            None => ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        };
        let response = fs_project::interface_state::resistance_sensitivity(card.claims(), &binding.state,
            &requirements.temperature_axis, resolved.range_lo, resolved.range_hi, selection).map_err(project_error)?;
        if response.axis != axis || response.unit != parameter.target.unit()
            || response.nominal.to_bits() != sample.parameters[i].to_bits()
            || response.claim.0.to_hex() != property.selected_claim
            || response.resistance.to_bits() != property.value_lo.to_bits()
            || property.value_lo.to_bits() != property.value_hi.to_bits()
        { return Err(invalid("joint-state derivative does not reproduce its nominal resistance and source")); }
        let row = rows.iter().find(|r|r.str_field("target")==Some("contact-resistance-multiplier")
            && r.str_field("entity")==Some(parameter.entity.as_str()))
            .ok_or_else(|| invalid("missing normalized contact-resistance derivative"))?;
        if row.str_field("interface_card") != Some(binding.card.as_str())
            || row.f64_field("reference_value") != Some(1.0)
            || row.get("mapped") != Some(&J::Bool(false))
        { return Err(invalid("joint-state calibration needs an unmapped card-backed normalized contact control")); }
        // dT/dx = dT/ds|_{R -> sR,s=1} * (dR/dx)/R.
        coefficients[i] *= response.derivative / response.resistance;
        if !coefficients[i].is_finite() { return Err(invalid("joint-state thermal derivative is nonfinite")); }
    }
    Ok(())
}

#[cfg(test)]
#[path = "contact/tests.rs"]
mod tests;
