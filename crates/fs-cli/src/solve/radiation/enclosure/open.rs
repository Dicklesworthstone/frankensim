//! Native open-cavity exchange retains the black environment as an external
//! heat account. It never inserts that reservoir into the solid trace list.
use super::*;
use fs_conduction::radiation::{OpenDiffuseEnclosure, solve_with_open_gray_diffuse_enclosure};

pub(in crate::solve::radiation) struct OpenResult {
    model: OpenDiffuseEnclosure,
    ambient_temperature_k: f64,
    report: EnclosureReport,
    decomposition_residual_w: f64,
    contact_work: Option<ContactEnclosureWork>,
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
    let applied_out: f64 = applied.net_outward_heat_w[..n].iter().sum();
    let decomposition = solved.conduction.report.energy.robin_out_w
        + solved.conduction.report.energy.neumann_out_w
        - solved.convective_out_w
        - solved.prescribed_neumann_out_w
        - applied_out;
    check_steady(&model, &solved.radiation, declaration.heat_tolerance.value, decomposition)?;
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
        radiation_thermal_responses: None,
        combined_boundary: Some(solved.combined_boundary),
        enclosure: Some(EnclosureResult::Open(OpenResult {
            model,
            ambient_temperature_k: ambient_k,
            report: EnclosureReport::Steady(solved.radiation),
            decomposition_residual_w: decomposition,
            contact_work: None,
        })),
    })
}

fn check_steady(
    model: &OpenDiffuseEnclosure, radiation: &CoupledRadiationReport,
    tolerance: f64, decomposition: f64,
) -> Result<(), ConductionError> {
    let n = model.surfaces().len();
    let applied = &radiation.radiosity;
    let actual = &radiation.equilibrium.nonlinear_radiosity;
    let largest_area = model.view_factors().areas_m2().iter().copied().fold(0.0_f64, f64::max);
    if !tolerance.is_finite() || tolerance <= 0.0
        || !applied.net_outward_heat_w[..n].iter().sum::<f64>().is_finite()
        || !decomposition.is_finite() || decomposition.abs() > tolerance
        || !radiation.equilibrium.max_patch_heat_mismatch_w.is_finite()
        || radiation.equilibrium.max_patch_heat_mismatch_w > tolerance
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
    Ok(())
}

impl SolidSolution {
    /// Keep the contact owner's complete physical endpoint. The black bath is
    /// the final radiosity row, never part of the solid's outward heat sum.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::solve) fn from_contact_open_enclosure(
        conduction: ConductionSolution, model: OpenDiffuseEnclosure, ambient_temperature_k: f64,
        radiation: CoupledRadiationReport, applied_boundary: ThermalBoundary,
        convective_fluxes: Vec<RobinFlux>, convective_out_w: f64,
        work: ContactEnclosureWork, heat_tolerance_w: f64,
    ) -> Result<Self, ConductionError> {
        let actual = &radiation.equilibrium.nonlinear_radiosity;
        let decomposition = conduction.report.energy.robin_out_w - convective_out_w
            - actual.net_outward_heat_w[..model.surfaces().len()].iter().sum::<f64>();
        check_steady(&model, &radiation, heat_tolerance_w, decomposition)?;
        Ok(Self { conduction, convective_fluxes: Some(convective_fluxes),
            convective_out: Some(convective_out_w), radiation: None, pointwise: None,
            radiation_thermal_responses: Some(work.thermal_responses),
            combined_boundary: Some(applied_boundary),
            enclosure: Some(EnclosureResult::Open(OpenResult { model, ambient_temperature_k,
                report: EnclosureReport::Steady(radiation), decomposition_residual_w: decomposition, contact_work: Some(work),
            })),
        })
    }

    /// Retain the original capacity endpoint without re-solving or substituting
    /// a steady balance. Storage and prescribed reactions stay with its time owner.
    pub(in crate::solve) fn from_open_enclosure_endpoint(
        conduction: ConductionSolution, model: OpenDiffuseEnclosure,
        ambient_temperature_k: f64, endpoint: super::super::EnclosureEndpoint,
    ) -> Self {
        Self {
            conduction, convective_fluxes: Some(endpoint.convective_robin_fluxes),
            convective_out: Some(endpoint.convective_out_w), radiation: None, pointwise: None,
            radiation_thermal_responses: None, combined_boundary: Some(endpoint.combined_boundary),
            enclosure: Some(EnclosureResult::Open(OpenResult {
                model, ambient_temperature_k,
                report: EnclosureReport::Transient {
                    radiation: endpoint.report, time_integration: endpoint.time_integration,
                    dt_s: endpoint.dt_s, residual_j: endpoint.physical_residual_norm_j,
                    threshold_j: endpoint.physical_residual_tolerance_j,
                    energy_residual_j: endpoint.physical_energy_residual_j,
                },
                decomposition_residual_w: endpoint.decomposition_residual_w, contact_work: None,
            })),
        }
    }
}

impl OpenResult {
    pub(super) fn model(&self) -> (&OpenDiffuseEnclosure, f64) {
        (&self.model, self.ambient_temperature_k)
    }

    pub(super) fn receipt(
        &self,
        lowering: &LoweredRadiation,
        convective_out_w: f64,
        thermal_responses: Option<usize>,
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
        let (actual, applied_heat, mut detail, scope) = match &self.report {
            EnclosureReport::Steady(report) => (
                &report.equilibrium.nonlinear_radiosity, report.radiosity.net_outward_heat_w.as_slice(),
                format!("\"solid_solves\":{},\"surface_iterations\":{},\"max_heat_mismatch_w\":{},\"nonlinear_residual_w\":{},\"nonlinear_residual_threshold_w\":{}",
                    thermal_responses.unwrap_or(report.iterations + 1), report.iterations,
                    num(report.equilibrium.max_patch_heat_mismatch_w)?, num(report.equilibrium.residual_w)?,
                    num(report.equilibrium.residual_threshold_w)?),
                if self.contact_work.is_some() { CONTACT_SCOPE } else { SCOPE },
            ),
            EnclosureReport::Transient {
                radiation, time_integration, dt_s, residual_j, threshold_j, energy_residual_j,
            } => (
                &radiation.radiosity, radiation.applied_heat_w.as_slice(),
                format!("\"time_integration\":{},\"time_step_s\":{},\"solid_solves\":{},\"radiation_trials\":{},\"solid_iterations\":{},\"krylov_iterations\":{},\"max_heat_mismatch_w\":{},\"nonlinear_residual_j\":{},\"nonlinear_residual_threshold_j\":{},\"nonlinear_residual_w\":{},\"nonlinear_residual_threshold_w\":{},\"physical_energy_residual_j\":{}",
                    json_string(time_integration), num(*dt_s)?, thermal_responses.unwrap_or(radiation.iterations),
                    radiation.iterations, radiation.solid_iterations, radiation.krylov_iterations,
                    num(radiation.max_heat_mismatch_w)?, num(*residual_j)?, num(*threshold_j)?,
                    num(*residual_j / *dt_s)?, num(*threshold_j / *dt_s)?, num(*energy_residual_j)?),
                CAPACITY_SCOPE,
            ),
        };
        if actual.net_outward_heat_w.len() != n + 1 || applied_heat.len() < n {
            return Err(bad("open enclosure receipt lost its original surface and black-reservoir rows"));
        }
        if let Some(work) = self.contact_work {
            let EnclosureReport::Steady(report) = &self.report else {
                return Err(bad("steady contact work cannot replace capacity trajectory counts"));
            };
            detail.push_str(&format!(
                ",\"radiation_responses\":{},\"solid_iterations\":{},\"krylov_iterations\":{},\"residual_evaluations\":{},\"jacobian_assemblies\":{}",
                report.iterations + 1, work.nonlinear_iterations, work.krylov_iterations,
                work.residual_evaluations, work.jacobian_assemblies));
        }
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
                num(applied_heat[i])?, num(actual.net_outward_heat_w[i])?,
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
            "{{\"model\":\"open-gray-diffuse-enclosure\",\"surfaces\":[{}],\"view_factors\":[{}],\"view_factor_identity\":{},\"view_factor_evidence\":{origin},\"ambient_temperature_k\":{},\"ambient_effective_area_m2\":{},\"ambient_absorbed_w\":{},\"ambient_radiosity_w_m2\":{},\"ambient_irradiation_w_m2\":{},\"enclosure_energy_closure_w\":{},\"convective_out_w\":{},\"radiative_out_w\":{},\"nonlinear_radiative_out_w\":{},\"decomposition_residual_w\":{},{detail},\"radiosity_residual_w_m2\":{},\"authority\":\"Estimated\",\"no_claim\":{}}}",
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
            num(applied_heat[..n].iter().sum())?,
            num(actual.net_outward_heat_w[..n].iter().sum())?,
            num(self.decomposition_residual_w)?,
            num(actual.linear_residual_max_w_m2)?,
            json_string(scope),
        ))
    }
}

const SCOPE: &str = "Opaque diffuse gray physical surfaces exchange reflected radiation and see a single prescribed black environment through their exact row deficits. Actual mesh areas admit reciprocity; supplied factors are not geometrically certified. The black environment has no solid degree of freedom. Its power is independently evaluated and stays outside air convection. This steady path uses patch-mean fourth-power emission, frozen source-backed emissivity and the shared complete nonlinear residual gate; initialization needs a conductive or convective thermal anchor. Use explicit pointwise total-enthalpy storage for vacuum cooling histories. Requested steady nominal derivatives include complete reflected exchange, each physical emissivity and one shared bath temperature through the same thermal dual. Moving visibility, participating media, spectral or directional surfaces, continuum error bounds and experimental validation are not claimed.";
const CONTACT_SCOPE: &str = "Live source-backed contact and heterogeneous conduction share the original open gray-diffuse patch-mean radiosity iteration. The unchanged final field retains freshly queried contact and physical radiation residuals, reservoir reactions and watt gates. The explicit black environment has no solid degree of freedom; its independently evaluated heat stays outside the physical-surface outward sum and transported-air convection. Source emissivity queries, view factors and joint state remain fixed; supplied factors are not geometrically certified. All nested contact and numerical responses are counted. Requested steady selected-node derivatives compose the complete contact/material and reflected-radiosity tangent with original natural/air feedback; source, contact-law, emissivity and shared bath controls contract that same dual. No arbitrary-law convergence or uniqueness, pointwise emission, continuum bound or experimental validation is supplied.";
const CAPACITY_SCOPE: &str = "Original volumetric capacity, bulk conductivity and any live source contacts advance through the shared backward-Euler endpoint and open area-mean radiosity iteration. Every trial retains the same physical old field. Fresh physical-surface heat enters the nodal residual, prescribed-node storage reaction and solid energy balance; the independent black-reservoir row has no solid capacity and remains outside the solid outward heat sum. Original convection alone heats transported air. Fixed emissivity sources, geometry and view factors retain their admitted authority and domains. All nested thermal and numerical responses are counted. This composition supplies no open-capacity radiation adjoint, pointwise emission, arbitrary-law convergence or uniqueness, continuum or continuous-time error certificate, or experimental validation.";
