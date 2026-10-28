//! Chain the native normalized-resistance adjoint through the selected card.
//! Calibration coordinates are physical pressure/thickness/gap/torque, NOT
//! resistance multipliers. No additional thermal solve or probe is performed.
use super::*;
use fs_matdb::{ClaimId, ClaimSelection, SelectionPolicy};

/// Both steady and transient producers differentiate a unit multiplier on
/// the original contact resistance. Their report layouts differ; the source
/// chain rule below is the same and must use the sampled physical state.
fn control_row<'a>(report: &'a J, entity: &str) -> Result<&'a J> {
    let (rows, identity) = match report.str_field("schema") {
        Some("frankensim.cli.nominal-adjoint.v1") => {
            if !matches!(report.str_field("output"), Some(output::CONTACT | output::COMBINED)) {
                return Err(invalid("joint-state calibration requires the complete contact adjoint report"));
            }
            (array(report, "parameters", 256)?, "entity")
        }
        Some("fs-cli-native-transient-adjoint-v1" | "fs-cli-native-enthalpy-adjoint-v1") => {
            let peak = report.str_field("output") == Some("temperature-sampled-peak-material-adjoint");
            if !peak && report.str_field("output") != Some("temperature-final-material-adjoint") {
                return Err(invalid("joint-state trajectory calibration requires the selected material adjoint"));
            }
            let controls = report.get("material_controls")
                .ok_or_else(|| invalid("transient contact calibration has no material controls"))?;
            let steps = report.get("state_index").and_then(J::number_raw)
                .and_then(|value| value.parse::<usize>().ok())
                .filter(|steps| peak || *steps > 0)
                .ok_or_else(|| invalid("transient contact calibration has no selected state index"))?;
            if controls.f64_field("control_steps") != Some(steps as f64)
                || controls.str_field("parameter_scope") != Some(
                    "constant-through-selected-trajectory-independent-regional-effective-conductivity-and-named-matching-resistance-multipliers")
                || controls.str_field("units") != Some("K")
            {
                return Err(invalid("contact controls do not span the same selected trajectory in kelvin"));
            }
            (array(controls, "contact_interfaces", 64)?, "interface")
        }
        _ => return Err(invalid("unsupported normalized contact adjoint producer")),
    };
    let mut matches = rows.iter().filter(|row| {
        row.str_field(identity) == Some(entity)
            && (identity == "interface"
                || row.str_field("target") == Some("contact-resistance-multiplier"))
    });
    let row = matches.next()
        .ok_or_else(|| invalid("missing normalized contact-resistance derivative"))?;
    if matches.next().is_some() {
        return Err(invalid("ambiguous normalized contact-resistance derivative"));
    }
    if row.f64_field("reference_value") != Some(1.0)
        || row.get("mapped") != Some(&J::Bool(false))
        || (identity == "interface" && (row.f64_field("log_reference_value") != Some(0.0)
            || row.f64_field("face_pairs").is_none_or(|value|
                !value.is_finite() || value <= 0.0 || value.fract() != 0.0)))
    {
        return Err(invalid("joint-state calibration needs an unmapped card-backed unit resistance control"));
    }
    Ok(row)
}

/// The caller has already bound the full receipt to the selected physical state.
/// Return the normalized derivative; `convert` applies dR/dx at that state.
pub(super) fn transient_coefficient(report: &J, entity: &str) -> Result<f64> {
    control_row(report, entity)?.f64_field("dtemperature_dresistance_multiplier_k")
        .filter(|value| value.is_finite())
        .ok_or_else(|| invalid("transient contact control has no finite resistance derivative"))
}

pub(super) fn convert(model: &Model, sample: &Sample, report: &J, coefficients: &mut [f64]) -> Result<()> {
    let parameters = model.bound.study().parameters();
    if !parameters.iter().any(|p| p.target.contact_axis().is_some() && p.low != p.high) { return Ok(()); }
    if coefficients.len() != parameters.len() {
        return Err(invalid("joint-state calibration coefficient layout differs"));
    }
    // The marginal means generally differ from the base project. Reconstruct
    // exactly the calibration point, never query a card at the old base state.
    let project = model.project(&sample.parameters)?;
    let library = model.cards.library();
    let requirements = fs_project::BindingRequirements::thermal_steady_v1();
    let resolution = fs_project::resolve_bindings(&project.spec, &library, &requirements);
    if !resolution.admissible() { return Err(invalid("joint-state calibration cannot reproduce its card bindings")); }
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
        let row = control_row(report, &parameter.entity)?;
        if row.str_field("interface_card") != Some(binding.card.as_str()) {
            return Err(invalid("joint-state calibration contact card differs from its sampled source"));
        }
        // dT/dx = dT/ds|_{R -> sR,s=1} * (dR/dx)/R.
        coefficients[i] *= response.derivative / response.resistance;
        if !coefficients[i].is_finite() { return Err(invalid("joint-state thermal derivative is nonfinite")); }
    }
    Ok(())
}

#[cfg(test)]
#[path = "contact/tests.rs"]
mod tests;
