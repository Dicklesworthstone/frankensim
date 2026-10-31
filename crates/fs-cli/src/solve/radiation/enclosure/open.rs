//! Native open-cavity exchange retains the black environment as an external
//! heat account. It never inserts that reservoir into the solid trace list.
use super::*;
use fs_conduction::radiation::{OpenDiffuseEnclosure, solve_with_open_gray_diffuse_enclosure};

pub(in crate::solve::radiation) struct OpenResult {
    model: OpenDiffuseEnclosure,
    report: CoupledRadiationReport,
    decomposition_residual_w: f64,
}

pub(in crate::solve::radiation) fn bind(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    lowering: &LoweredRadiation,
) -> Result<OpenDiffuseEnclosure, ConductionError> {
    let declaration = &lowering
        .enclosure
        .as_ref()
        .ok_or_else(|| physical("missing open-enclosure declaration"))?
        .declaration;
    if declaration.ambient_temperature.is_none() {
        return Err(physical(
            "open exchange requires its explicit black environment",
        ));
    }
    OpenDiffuseEnclosure::new(
        cx,
        bind_surfaces(cx, problem, lowering)?,
        declaration.view_factors.clone(),
        evidence(declaration),
        ViewFactorTolerance {
            row_sum_abs: declaration.closure_tolerance,
            reciprocity_rel: declaration.reciprocity_relative_tolerance,
        },
    )
}

pub(super) fn solve(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    config: SolveConfig,
    lowering: &LoweredRadiation,
) -> Result<SolidSolution, ConductionError> {
    let declaration = &lowering
        .enclosure
        .as_ref()
        .ok_or_else(|| physical("missing open-enclosure declaration"))?
        .declaration;
    let ambient_k = declaration
        .ambient_temperature
        .ok_or_else(|| physical("missing black environment temperature"))?
        .value;
    let model = bind(cx, problem, lowering)?;
    let mut solved = solve_with_open_gray_diffuse_enclosure(
        cx,
        problem,
        interfaces,
        &model,
        ambient_k,
        config,
        CoupledRadiationConfig {
            max_iterations: declaration.max_iterations as usize,
            surface_temperature_atol_k: declaration.temperature_tolerance.value,
            surface_temperature_rtol: 0.0,
            relaxation: declaration.relaxation,
        },
    )?;
    let n = model.surfaces().len();
    let applied = &solved.radiation.radiosity;
    let actual = &solved.radiation.equilibrium.nonlinear_radiosity;
    let applied_out: f64 = applied.net_outward_heat_w[..n].iter().sum();
    let decomposition = solved.conduction.report.energy.robin_out_w
        + solved.conduction.report.energy.neumann_out_w
        - solved.convective_out_w
        - solved.prescribed_neumann_out_w
        - applied_out;
    let tolerance = declaration.heat_tolerance.value;
    let largest_area = model
        .view_factors()
        .areas_m2()
        .iter()
        .copied()
        .fold(0.0_f64, f64::max);
    if !tolerance.is_finite()
        || tolerance <= 0.0
        || !applied_out.is_finite()
        || !decomposition.is_finite()
        || decomposition.abs() > tolerance
        || !solved
            .radiation
            .equilibrium
            .max_patch_heat_mismatch_w
            .is_finite()
        || solved.radiation.equilibrium.max_patch_heat_mismatch_w > tolerance
        || [applied, actual].iter().any(|report| {
            !report.enclosure_energy_closure_w.is_finite()
                || report.enclosure_energy_closure_w.abs() > tolerance
                || !report.linear_residual_max_w_m2.is_finite()
                || report.linear_residual_max_w_m2 * largest_area > tolerance
        })
    {
        return Err(physical(
            "open exchange failed its declared heat, radiosity, solid/bath closure or convection decomposition allowance",
        ));
    }
    let energy = &mut solved.conduction.report.energy;
    energy.neumann_out_w = solved.prescribed_neumann_out_w;
    energy.robin_out_w = solved.convective_out_w + applied_out;
    energy.closure_w =
        energy.source_w + energy.dirichlet_in_w - energy.neumann_out_w - energy.robin_out_w;
    energy.scale_w = energy
        .source_w
        .abs()
        .max(energy.dirichlet_in_w.abs())
        .max(energy.neumann_out_w.abs())
        .max(energy.robin_out_w.abs())
        .max(f64::MIN_POSITIVE);
    poll(cx, solved.radiation.iterations)?;
    Ok(SolidSolution {
        conduction: solved.conduction,
        convective_fluxes: Some(solved.convective_robin_fluxes),
        convective_out: Some(solved.convective_out_w),
        radiation: None,
        pointwise: None,
        combined_boundary: Some(solved.combined_boundary),
        enclosure: Some(EnclosureResult::Open(OpenResult {
            model,
            report: solved.radiation,
            decomposition_residual_w: decomposition,
        })),
    })
}

impl OpenResult {
    pub(super) fn receipt(
        &self,
        lowering: &LoweredRadiation,
        convective_out_w: f64,
    ) -> Result<String, SolveRefusal> {
        let declaration = &lowering
            .enclosure
            .as_ref()
            .ok_or_else(|| bad("open exchange lost its declaration"))?
            .declaration;
        let ambient_k = declaration
            .ambient_temperature
            .ok_or_else(|| bad("open exchange lost its black environment"))?
            .value;
        let n = self.model.surfaces().len();
        let actual = &self.report.equilibrium.nonlinear_radiosity;
        let applied = &self.report.radiosity;
        let mut rows = Vec::with_capacity(n);
        for (i, (surface, source)) in self
            .model
            .surfaces()
            .iter()
            .zip(&lowering.surfaces)
            .enumerate()
        {
            rows.push(format!(
                "{{\"name\":{},\"target\":{},\"card\":{},\"material_state\":{},\"emissivity_receipt\":{},\"emissivity\":{},\"query_temperature_k\":{},\"area_m2\":{},\"mean_temperature_k\":{},\"ambient_view_factor\":{},\"applied_heat_w\":{},\"nonlinear_heat_w\":{},\"radiosity_w_m2\":{},\"irradiation_w_m2\":{}}}",
                json_string(&source.name), json_string(&source.target),
                json_string(&source.emissivity.card_identity().to_hex()),
                json_string(source.emissivity.material_state()),
                json_string(&source.emissivity.receipt().content_hash().to_hex()),
                num(source.emissivity.value())?, num(source.emissivity.temperature_k())?,
                num(surface.area_m2())?, num(actual.surface_temperatures_k[i])?,
                num(self.model.view_factors().factors()[i][n])?,
                num(applied.net_outward_heat_w[i])?, num(actual.net_outward_heat_w[i])?,
                num(actual.radiosity_w_m2[i])?, num(actual.irradiation_w_m2[i])?,
            ));
        }
        let factors = self.model.view_factors().factors()[..n]
            .iter()
            .map(|row| {
                row[..n]
                    .iter()
                    .map(|v| num(*v))
                    .collect::<Result<Vec<_>, _>>()
                    .map(|values| format!("[{}]", values.join(",")))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let origin = match &declaration.evidence {
            EnclosureViewFactorEvidence::Analytic { geometry } => format!(
                "{{\"kind\":\"analytic\",\"geometry\":{}}}",
                json_string(geometry)
            ),
            EnclosureViewFactorEvidence::ExternalQmc {
                seed,
                samples,
                generator,
            } => format!(
                "{{\"kind\":\"external-qmc\",\"seed\":{seed},\"samples\":{samples},\"generator\":{}}}",
                json_string(generator)
            ),
        };
        Ok(format!(
            "{{\"model\":\"open-gray-diffuse-enclosure\",\"surfaces\":[{}],\"view_factors\":[{}],\"view_factor_identity\":{},\"view_factor_evidence\":{origin},\"ambient_temperature_k\":{},\"ambient_effective_area_m2\":{},\"ambient_absorbed_w\":{},\"ambient_radiosity_w_m2\":{},\"ambient_irradiation_w_m2\":{},\"enclosure_energy_closure_w\":{},\"convective_out_w\":{},\"radiative_out_w\":{},\"nonlinear_radiative_out_w\":{},\"decomposition_residual_w\":{},\"solid_solves\":{},\"surface_iterations\":{},\"max_heat_mismatch_w\":{},\"nonlinear_residual_w\":{},\"nonlinear_residual_threshold_w\":{},\"radiosity_residual_w_m2\":{},\"authority\":\"Estimated\",\"no_claim\":{}}}",
            rows.join(","),
            factors.join(","),
            json_string(&self.model.view_factors().identity().to_hex()),
            num(ambient_k)?,
            num(self.model.view_factors().areas_m2()[n])?,
            num(-actual.net_outward_heat_w[n])?,
            num(actual.radiosity_w_m2[n])?,
            num(actual.irradiation_w_m2[n])?,
            num(actual.enclosure_energy_closure_w)?,
            num(convective_out_w)?,
            num(applied.net_outward_heat_w[..n].iter().sum())?,
            num(actual.net_outward_heat_w[..n].iter().sum())?,
            num(self.decomposition_residual_w)?,
            self.report.iterations + 1,
            self.report.iterations,
            num(self.report.equilibrium.max_patch_heat_mismatch_w)?,
            num(self.report.equilibrium.residual_w)?,
            num(self.report.equilibrium.residual_threshold_w)?,
            num(actual.linear_residual_max_w_m2)?,
            json_string(SCOPE),
        ))
    }
}

const SCOPE: &str = "Opaque diffuse gray physical surfaces exchange reflected radiation and see a single prescribed black environment through their exact row deficits. Actual mesh areas admit reciprocity; supplied factors are not geometrically certified. The black environment has no solid degree of freedom. Its power is independently evaluated and stays outside air convection. This steady path uses patch-mean fourth-power emission, frozen source-backed emissivity and the shared complete nonlinear residual gate; initialization needs a conductive or convective thermal anchor. Use explicit pointwise total-enthalpy storage for vacuum cooling histories. Native open-bath adjoint controls, moving visibility, participating media, spectral or directional surfaces, continuum error bounds and experimental validation are not claimed.";
