//! Source-preserving controls of the complete capacity/enclosure transpose.
//! The numerical boundary, reflected tangent and remap belong to fs-conduction.

use fs_conduction::radiation::GrayDiffuseEnclosure;
use crate::solve::radiation::LoweredRadiation;

use super::{SolveRefusal, bad};
use super::super::super::{json_string, number};

pub(super) fn fragment(
    lowering: &LoweredRadiation,
    model: &GrayDiffuseEnclosure,
    emissivities: &[f64],
    steps: usize,
    radiosity_iterations: usize,
    worst_radiosity_residual: f64,
) -> Result<String, SolveRefusal> {
    if model.surfaces().len() != lowering.surface_count()
        || model.surfaces().len() != emissivities.len()
    {
        return Err(bad("capacity enclosure controls lost their declared surface layout"));
    }
    let mut rows = Vec::with_capacity(emissivities.len());
    for (index, ((name, target), surface)) in lowering.surface_labels()
        .zip(model.surfaces()).enumerate()
    {
        if target != surface.name() {
            return Err(bad("capacity enclosure control order differs from its view factors"));
        }
        let emissivity = surface.emissivity();
        let directions = if emissivity.value() == 1.0 {
            "nonpositive-at-emissivity-one"
        } else {
            "two-sided-within-zero-exclusive-one-inclusive"
        };
        rows.push(format!(
            "{{\"name\":{},\"target\":{},\"card\":{},\"material_state\":{},\"emissivity_receipt\":{},\"emissivity\":{},\"query_temperature_k\":{},\"dtemperature_demissivity_k\":{},\"emissivity_feasible_directions\":{}}}",
            json_string(name), json_string(target),
            json_string(&emissivity.card_identity().to_hex()),
            json_string(emissivity.material_state()),
            json_string(&emissivity.receipt().content_hash().to_hex()),
            number(emissivity.value())?, number(emissivity.temperature_k())?,
            number(emissivities[index])?, json_string(directions),
        ));
    }
    Ok(format!(
        ",\"radiation\":{{\"method\":\"complete-closed-gray-diffuse-capacity-tangent\",\"patch_count\":{},\"control_steps\":{},\"radiosity_adjoint_krylov_iterations\":{},\"maximum_radiosity_dual_relative_residual\":{},\"parameters\":\"absolute-per-surface-emissivity-throughout-selected-trajectory\",\"view_factor_identity\":{},\"surfaces\":[{}]}}",
        emissivities.len(), steps, radiosity_iterations, number(worst_radiosity_residual)?,
        json_string(&model.view_factors().identity().to_hex()), rows.join(","),
    ))
}
