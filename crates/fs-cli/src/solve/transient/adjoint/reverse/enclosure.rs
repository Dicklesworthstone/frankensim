//! Source-preserving controls of the complete capacity/enclosure transpose.
//! The numerical boundary, reflected tangent and remap belong to fs-conduction.

use crate::solve::radiation::LoweredRadiation;

use super::{EnclosureModel, SolveRefusal, bad};
use super::super::super::{json_string, number};

pub(super) fn fragment(
    lowering: &LoweredRadiation,
    model: &EnclosureModel,
    emissivities: &[f64],
    ambient_derivative: f64,
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
    let (method, parameters, factors, bath) = match model {
        EnclosureModel::Closed(model) => (
            "complete-closed-gray-diffuse-capacity-tangent",
            "absolute-per-surface-emissivity-throughout-selected-trajectory",
            model.view_factors(), String::new()),
        EnclosureModel::Open(model, ambient) => (
            "complete-open-gray-diffuse-capacity-tangent",
            "absolute-per-surface-emissivity-and-shared-black-environment-temperature-throughout-selected-trajectory",
            model.view_factors(), format!(
                ",\"black_environment\":{{\"temperature_k\":{},\"dtemperature_dambient_temperature_k_per_k\":{},\"feasible_directions\":\"two-sided-within-positive-kelvin\"}}",
                number(*ambient)?, number(ambient_derivative)?)),
    };
    Ok(format!(
        ",\"radiation\":{{\"method\":\"{method}\",\"patch_count\":{},\"control_steps\":{},\"radiosity_adjoint_krylov_iterations\":{},\"maximum_radiosity_dual_relative_residual\":{},\"parameters\":\"{parameters}\",\"view_factor_identity\":{}{bath},\"surfaces\":[{}]}}",
        emissivities.len(), steps, radiosity_iterations, number(worst_radiosity_residual)?,
        json_string(&factors.identity().to_hex()), rows.join(","),
    ))
}
