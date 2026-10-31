//! Actual local-emission reports, never relabeled as a patch-mean fixed point.
use fs_conduction::{PointwiseEnclosureReport, PointwiseRadiationReport};
use super::{LoweredRadiation, SolveRefusal, SurfaceSource, conduction_error, json_string, num};

pub(super) enum PointwiseEndpoint {
    Ambient(PointwiseRadiationReport),
    Enclosure(PointwiseEnclosureReport),
}

fn bad(what: &str) -> SolveRefusal {
    conduction_error(
        "cli-solve-conduction-pointwise-radiation", what,
        "retain the declared spatial radiation model and its actual endpoint report",
    )
}

impl PointwiseEndpoint {
    pub(super) fn new(
        ambient: Option<PointwiseRadiationReport>,
        enclosure: Option<PointwiseEnclosureReport>,
    ) -> Result<Self, SolveRefusal> {
        match (ambient, enclosure) {
            (Some(report), None) => Ok(Self::Ambient(report)),
            (None, Some(report)) => Ok(Self::Enclosure(report)),
            _ => Err(bad("a pointwise endpoint must retain exactly its declared radiation report")),
        }
    }

    pub(super) fn receipt(
        &self, lowering: &LoweredRadiation, convective_out_w: f64,
    ) -> Result<String, SolveRefusal> {
        match self {
            Self::Ambient(report) => receipt(Some(report), None, lowering, convective_out_w),
            Self::Enclosure(report) => receipt(None, Some(report), lowering, convective_out_w),
        }
    }
}

fn source<'a>(lowering: &'a LoweredRadiation, target: &str) -> Result<&'a SurfaceSource, SolveRefusal> {
    lowering.surfaces.iter().find(|row| row.target == target)
        .ok_or_else(|| bad("pointwise radiation report has an undeclared surface"))
}

fn source_fields(source: &SurfaceSource) -> Result<String, SolveRefusal> {
    Ok(format!(
        "\"name\":{},\"target\":{},\"card\":{},\"material_state\":{},\"emissivity_receipt\":{},\"emissivity\":{},\"query_temperature_k\":{}",
        json_string(&source.name), json_string(&source.target),
        json_string(&source.emissivity.card_identity().to_hex()),
        json_string(source.emissivity.material_state()),
        json_string(&source.emissivity.receipt().content_hash().to_hex()),
        num(source.emissivity.value())?, num(source.emissivity.temperature_k())?,
    ))
}

pub(in crate::solve) fn receipt(
    ambient: Option<&PointwiseRadiationReport>,
    enclosure: Option<&PointwiseEnclosureReport>,
    lowering: &LoweredRadiation,
    convective_out_w: f64,
) -> Result<String, SolveRefusal> {
    if !lowering.is_pointwise() {
        return Err(bad("pointwise radiation report has a patch-mean declaration"));
    }
    let (model, heat, rows, extra) = match (ambient, enclosure) {
        (Some(report), None) if !lowering.is_enclosure() => {
            let rows = report.patches.iter().map(|row| {
                let source = source(lowering, row.patch.region())?;
                Ok(format!(
                    "{{{},\"ambient_temperature_k\":{},\"area_m2\":{},\"mean_temperature_k\":{},\"minimum_temperature_k\":{},\"maximum_temperature_k\":{},\"outward_heat_w\":{}}}",
                    source_fields(source)?, num(row.patch.ambient_temperature_k())?,
                    num(row.area_m2)?, num(row.mean_surface_temperature_k)?,
                    num(row.minimum_surface_temperature_k)?, num(row.maximum_surface_temperature_k)?,
                    num(row.outward_heat_w)?,
                ))
            }).collect::<Result<Vec<_>, SolveRefusal>>()?;
            ("pointwise-p1-gray-surface-to-black-reservoir", report.outward_heat_w, rows, String::new())
        }
        (None, Some(report)) if lowering.is_enclosure() => {
            if lowering.open_ambient_temperature_k() != report.ambient.as_ref().map(|bath| bath.temperature_k) {
                return Err(bad("pointwise enclosure report lost its declared black environment"));
            }
            let rows = report.patches.iter().map(|row| {
                let source = source(lowering, &row.name)?;
                Ok(format!(
                    "{{{},\"area_m2\":{},\"mean_temperature_k\":{},\"minimum_temperature_k\":{},\"maximum_temperature_k\":{},\"mean_blackbody_emission_w_m2\":{},\"radiosity_w_m2\":{},\"irradiation_w_m2\":{},\"outward_heat_w\":{}}}",
                    source_fields(source)?, num(row.area_m2)?, num(row.mean_temperature_k)?,
                    num(row.minimum_temperature_k)?, num(row.maximum_temperature_k)?,
                    num(row.mean_blackbody_emission_w_m2)?, num(row.radiosity_w_m2)?,
                    num(row.irradiation_w_m2)?, num(row.outward_heat_w)?,
                ))
            }).collect::<Result<Vec<_>, SolveRefusal>>()?;
            let mut extra = format!(
                ",\"view_factor_identity\":{},\"radiosity_residual_w_m2\":{},\"absolute_patch_heat_w\":{},\"enclosure_energy_closure_w\":{},\"radiosity_heat_tolerance_w\":{},\"radiosity_relative_heat_tolerance\":{}",
                json_string(&report.view_factor_identity.to_hex()),
                num(report.linear_residual_max_w_m2)?, num(report.absolute_patch_heat_w)?,
                num(report.system_energy_closure_w)?, num(lowering.config.balance_tolerance_w)?,
                num(lowering.config.balance_relative_tolerance)?,
            );
            let model = if let Some(bath) = &report.ambient {
                extra.push_str(&format!(
                    ",\"ambient_temperature_k\":{},\"ambient_effective_area_m2\":{},\"ambient_absorbed_w\":{},\"ambient_radiosity_w_m2\":{},\"ambient_irradiation_w_m2\":{}",
                    num(bath.temperature_k)?, num(bath.effective_area_m2)?,
                    num(-bath.outward_heat_w)?, num(bath.radiosity_w_m2)?,
                    num(bath.irradiation_w_m2)?,
                ));
                "pointwise-p1-open-gray-diffuse-enclosure"
            } else {
                "pointwise-p1-closed-gray-diffuse-enclosure"
            };
            (model, report.outward_heat_w, rows, extra)
        }
        _ => return Err(bad("pointwise radiation report differs from its declared model")),
    };
    if rows.len() != lowering.surface_count() {
        return Err(bad("pointwise radiation omitted a declared surface"));
    }
    Ok(format!(
        "{{\"model\":\"{model}\",\"spatial_model\":\"pointwise\",\"surfaces\":[{}],\"convective_out_w\":{},\"radiative_out_w\":{},\"nonlinear_radiative_out_w\":{}{extra},\"solver\":\"enthalpy-newton-fgmres\",\"authority\":\"Estimated\",\"no_claim\":{}}}",
        rows.join(","), num(convective_out_w)?, num(heat)?, num(heat)?,
        json_string(SCOPE),
    ))
}

const SCOPE: &str = "Local P1 fourth-power emission and its analytic tangent enter the same total-enthalpy Newton solve as conduction, contact and latent storage. Positive degree-five surface quadrature retains the actual nodal temperature field. Fixed emissivity comes from the declared source query and every surface node must remain in that claim's temperature support. Ambient surfaces see their declared black reservoir with unit view factor; enclosures retain patch-uniform irradiation from the declared reciprocal view factors. Open-enclosure row deficits see one explicit black environment whose independently evaluated power remains outside the solid degrees of freedom and air convection. Physical surface heat enters solid storage once; enclosure_energy_closure_w includes the bath when present. Original convection alone heats air. Reported mean temperature is an observation, not the emission input. Supported closed/ambient trajectory adjoints separately bind the complete spatial tangent and reflected transpose; native open-bath adjoint controls are not yet exposed. No moving visibility, participating medium, emissivity-temperature derivative, continuum bound or experimental validation is inferred.";
