//! Report the actual nodal emissivities used by the material radiation owner.
//! The reference query selects a source; it is not an emission coefficient.

use super::{LoweredRadiation, PointwiseRadiationReport, SolveRefusal, bad, json_string, num, source};

pub(super) fn receipt(
    report: &PointwiseRadiationReport,
    lowering: &LoweredRadiation,
    convective_out_w: f64,
) -> Result<String, SolveRefusal> {
    if lowering.is_enclosure() || report.patches.len() != lowering.surface_count() {
        return Err(bad("local material emissivity report differs from its declared surfaces"));
    }
    let rows = report.patches.iter().map(|row| {
        let source = source(lowering, row.patch.region())?;
        if row.material_states.is_empty() {
            return Err(bad("local material emissivity cannot use a frozen source snapshot"));
        }
        let mut minimum = f64::INFINITY;
        let mut maximum = f64::NEG_INFINITY;
        let states = row.material_states.iter().map(|state| {
            let emissivity = &state.emissivity;
            if emissivity.card_identity() != source.emissivity.card_identity()
                || emissivity.receipt().selected != source.emissivity.receipt().selected
                || state.phase.material_card_identity() != emissivity.card_identity()
                || state.material.card_identity() != emissivity.card_identity()
            {
                return Err(bad("local emissivity report lost its original material and selected source"));
            }
            minimum = minimum.min(emissivity.value());
            maximum = maximum.max(emissivity.value());
            Ok(format!(
                concat!("{{\"vertex\":{},\"specific_enthalpy_j_kg\":{},\"temperature_k\":{},",
                    "\"liquid_mass_fraction\":{},\"emissivity\":{},\"phase_curve_identity\":{},",
                    "\"material_state_identity\":{},\"selected_claim\":{},\"usage_receipt\":{}}}"),
                state.vertex, num(state.phase.specific_enthalpy_j_kg())?,
                num(state.phase.temperature_k())?, num(state.phase.liquid_mass_fraction())?,
                num(emissivity.value())?, json_string(&state.phase.phase_curve_identity().to_hex()),
                json_string(&state.material.identity().to_hex()),
                json_string(&emissivity.receipt().selected.0.to_hex()),
                json_string(&emissivity.receipt().content_hash().to_hex()),
            ))
        }).collect::<Result<Vec<_>, SolveRefusal>>()?;
        Ok(format!(
            concat!("{{\"name\":{},\"target\":{},\"card\":{},\"material_state\":{},",
                "\"selected_claim\":{},\"reference_emissivity_receipt\":{},",
                "\"reference_emissivity\":{},\"query_temperature_k\":{},",
                "\"ambient_temperature_k\":{},\"area_m2\":{},\"mean_temperature_k\":{},",
                "\"minimum_temperature_k\":{},\"maximum_temperature_k\":{},",
                "\"minimum_emissivity\":{},\"maximum_emissivity\":{},",
                "\"outward_heat_w\":{},\"material_states\":[{}]}}"),
            json_string(&source.name), json_string(&source.target),
            json_string(&source.emissivity.card_identity().to_hex()),
            json_string(source.emissivity.material_state()),
            json_string(&source.emissivity.receipt().selected.0.to_hex()),
            json_string(&source.emissivity.receipt().content_hash().to_hex()),
            num(source.emissivity.value())?, num(source.emissivity.temperature_k())?,
            num(row.patch.ambient_temperature_k())?, num(row.area_m2)?,
            num(row.mean_surface_temperature_k)?, num(row.minimum_surface_temperature_k)?,
            num(row.maximum_surface_temperature_k)?, num(minimum)?, num(maximum)?,
            num(row.outward_heat_w)?, states.join(","),
        ))
    }).collect::<Result<Vec<_>, SolveRefusal>>()?;
    Ok(format!(
        concat!("{{\"model\":\"pointwise-p1-material-gray-surface-to-black-reservoir\",",
            "\"spatial_model\":\"pointwise\",\"emissivity_response\":\"local-temperature\",",
            "\"surfaces\":[{}],\"convective_out_w\":{},\"radiative_out_w\":{},",
            "\"nonlinear_radiative_out_w\":{},\"solver\":\"enthalpy-newton-fgmres\",",
            "\"authority\":\"Estimated\",\"no_claim\":{}}}"),
        rows.join(","), num(convective_out_w)?, num(report.outward_heat_w)?,
        num(report.outward_heat_w)?, json_string(SCOPE),
    ))
}

const SCOPE: &str = "The original selected scalar or linear temperature-curve emissivity is queried at every surface node's actual caloric state. Nodal temperature and emissivity are interpolated on the same P1 trace; positive degree-six surface quadrature integrates their product with fourth-power emission. The complete emissivity-temperature feedback enters the original enthalpy Newton operator and accepted physical energy audit. The reference query only selects and identifies the retained source. Every radiating node must use that same material card and exact caloric owner; unsupported source coordinates, source switching and extrapolation refuse. Each surface sees its declared black reservoir with unit view factor. Original convection alone heats transported air. No separate coating, spectral or participating-medium model, moving visibility, phase-dependent emissivity, continuum bound or experimental validation is inferred.";
