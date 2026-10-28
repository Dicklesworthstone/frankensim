//! Native closed-enclosure exchange through the shared conduction/radiosity owner.
//! The declared matrix order is preserved; actual mesh areas admit reciprocity.
//! Only original convection reaches air transport. Internal radiation is separate.

use std::collections::BTreeSet;

use fs_conduction::radiation::{
    CoupledRadiationConfig, CoupledRadiationReport, GrayDiffuseEnclosure, RadiationSurface,
    ViewFactorEvidence, ViewFactorMatrix, ViewFactorTolerance,
    solve_with_convective_gray_diffuse_enclosure,
};
use fs_conduction::{ConductionError, ConductionProblem, SolveConfig, ThermalInterfaces};
use fs_exec::Cx;
use fs_conduction::transient::enthalpy::enclosure::EnthalpyEnclosureReport;
use fs_project::spec::{ConductionSetup, ThermalBoundaryCondition};
use fs_project::{ConductionEnclosure, EnclosureViewFactorEvidence, ProjectSpec};

use super::{
    AmbientRadiationConfig, CardPackSet, LoweredRadiation, SolidSolution, SolveRefusal,
    SurfaceSource, conduction_error, json_string, num, resolve_emissivity,
};

pub(super) struct LoweredEnclosure {
    declaration: ConductionEnclosure,
}

pub(super) struct EnclosureResult {
    model: GrayDiffuseEnclosure,
    report: EnclosureReport,
    decomposition_residual_w: f64,
}

enum EnclosureReport {
    Steady(CoupledRadiationReport),
    Transient {
        radiation: EnthalpyEnclosureReport,
        time_integration: &'static str,
        dt_s: f64,
        residual_j: f64,
        threshold_j: f64,
        energy_residual_j: f64,
    },
}

fn bad(what: impl Into<String>) -> SolveRefusal {
    conduction_error(
        "cli-solve-conduction-enclosure",
        what,
        "declare sourced emissivities, explicit ordered view factors, uniform convection or prescribed-flux traces, and explicit capacity or total-enthalpy storage for a trajectory",
    )
}

fn physical(what: impl Into<String>) -> ConductionError {
    ConductionError::Radiation {
        surface: "native-enclosure".into(),
        what: what.into(),
        fix: "inspect the exact patch partition and reciprocal view factors; tighten solid/surface tolerances or increase the declared coupling budget".into(),
    }
}

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "native-enclosure-binding",
        at,
    })
}

/// Unsupported numerical claims refuse before geometry or a physical solve.
pub(super) fn lower(
    spec: &ProjectSpec,
    setup: &ConductionSetup,
    cards: &CardPackSet,
) -> Result<LoweredRadiation, SolveRefusal> {
    let declaration = setup
        .enclosure_radiation
        .as_ref()
        .ok_or_else(|| bad("missing closed-enclosure declaration"))?;
    if setup.radiation.is_some() {
        return Err(bad(
            "ambient-reservoir and closed-enclosure radiation cannot own the same native model",
        ));
    }
    if spec
        .solver
        .as_ref()
        .is_some_and(|solver| solver.fidelity == "adaptive")
        || spec.outputs.iter().flatten().any(|output| {
            output.name == "temperature-volume-mean-bound"
        })
    {
        return Err(bad(
            "enclosure adaptive goals and continuum mean bounds require their own complete enclosure error estimator; a frozen radiative load is not that model",
        ));
    }
    if setup.transient.is_some() {
        // The trajectory owner admits final/selected-prefix goals and rejects
        // steady objectives before any physical solve. It also owns the
        // storage, differentiated natural-card and time-grid contracts.
        super::super::transient::admit(spec)?;
    } else {
        let mut adjoint = false;
        for output in spec.outputs.iter().flatten().filter(|output| output.name.contains("adjoint")) {
            if !matches!(output.name.as_str(),
                "temperature-max-adjoint" | "temperature-max-contact-adjoint"
                    | "temperature-max-boundary-adjoint" | "temperature-max-contact-boundary-adjoint"
                    | "temperature-volume-mean-adjoint")
            {
                return Err(bad("the requested enclosure adjoint is not a supported native steady objective"));
            }
            adjoint = true;
        }
        if adjoint {
            super::super::nominal_adjoint::requested(spec)?;
        }
    }
    if !(2..=64).contains(&declaration.surfaces.len()) {
        return Err(bad("a native enclosure needs 2..=64 ordered surfaces"));
    }
    let mut targets = BTreeSet::new();
    let mut names = BTreeSet::new();
    let mut surfaces = Vec::with_capacity(declaration.surfaces.len());
    for surface in &declaration.surfaces {
        if !names.insert(surface.name.as_str()) || !targets.insert(surface.target.as_str()) {
            return Err(bad("enclosure names and physical targets must be unique"));
        }
        let law = setup
            .boundaries
            .iter()
            .find(|boundary| boundary.target == surface.target)
            .ok_or_else(|| {
                bad(format!(
                    "enclosure target {} has no declared boundary",
                    surface.target
                ))
            })?;
        if !matches!(
            law.condition,
            ThermalBoundaryCondition::Convection { .. }
                | ThermalBoundaryCondition::AirflowConvection { .. }
                | ThermalBoundaryCondition::NaturalConvection { .. }
                | ThermalBoundaryCondition::HeatFlux { .. }
        ) {
            return Err(bad(
                "each native enclosure patch must retain its original convection or prescribed outward flux; zero flux declares vacuum",
            ));
        }
        surfaces.push(SurfaceSource {
            name: surface.name.clone(),
            target: surface.target.clone(),
            emissivity: resolve_emissivity(
                &surface.name,
                &surface.card,
                surface.claim.as_deref(),
                surface.query_temperature.value,
                cards,
            )?,
        });
    }
    Ok(LoweredRadiation {
        patches: Vec::new(),
        config: AmbientRadiationConfig {
            max_iterations: declaration.max_iterations as usize,
            temperature_tolerance_k: declaration.temperature_tolerance.value,
            balance_tolerance_w: declaration.heat_tolerance.value,
            balance_relative_tolerance: 0.0,
            relaxation: declaration.relaxation,
        },
        surfaces,
        enclosure: Some(LoweredEnclosure {
            declaration: declaration.clone(),
        }),
    })
}

fn evidence(declaration: &ConductionEnclosure) -> ViewFactorEvidence {
    match &declaration.evidence {
        EnclosureViewFactorEvidence::Analytic { geometry } => ViewFactorEvidence::Analytic {
            geometry: geometry.clone(),
        },
        EnclosureViewFactorEvidence::ExternalQmc {
            seed,
            samples,
            generator,
        } => ViewFactorEvidence::ExternalQmc {
            seed: *seed,
            samples: *samples,
            generator: generator.clone(),
        },
    }
}

/// Bind by the existing boundary's exact topological faces, never by a second
/// geometric selector or a user-supplied area. Rebinding happens on every mesh.
pub(super) fn bind(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    lowering: &LoweredRadiation,
) -> Result<GrayDiffuseEnclosure, ConductionError> {
    let declaration = &lowering
        .enclosure
        .as_ref()
        .ok_or_else(|| physical("missing enclosure lowering"))?
        .declaration;
    let mut surfaces = Vec::with_capacity(lowering.surfaces.len());
    for (ordinal, source) in lowering.surfaces.iter().enumerate() {
        poll(cx, ordinal)?;
        let region = problem
            .boundary
            .region_names()
            .iter()
            .position(|name| name == &source.target)
            .ok_or_else(|| {
                physical(format!(
                    "enclosure target {} did not bind a boundary",
                    source.target
                ))
            })?;
        let mut faces = BTreeSet::new();
        for (slot, face) in problem.mesh.boundary().iter().enumerate() {
            if slot % 256 == 0 {
                poll(cx, slot)?;
            }
            if problem.boundary.region_for(slot) == Some(region) {
                faces.insert(face.vertices);
            }
        }
        // The numerical owner keys its exact physical trace by target. The
        // distinct authored surface name remains in the retained source rows.
        surfaces.push(RadiationSurface::new(
            problem.mesh,
            &source.target,
            |face| faces.contains(&face.vertices),
            source.emissivity.clone(),
        )?);
    }
    let areas = surfaces.iter().map(RadiationSurface::area_m2).collect();
    let factors = ViewFactorMatrix::admit(
        areas,
        declaration.view_factors.clone(),
        evidence(declaration),
        ViewFactorTolerance {
            row_sum_abs: declaration.closure_tolerance,
            reciprocity_rel: declaration.reciprocity_relative_tolerance,
        },
    )?;
    poll(cx, surfaces.len())?;
    GrayDiffuseEnclosure::new(surfaces, factors)
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
        .ok_or_else(|| physical("missing enclosure lowering"))?
        .declaration;
    let model = bind(cx, problem, lowering)?;
    let mut solved = solve_with_convective_gray_diffuse_enclosure(
        cx,
        problem,
        interfaces,
        &model,
        config,
        CoupledRadiationConfig {
            max_iterations: declaration.max_iterations as usize,
            surface_temperature_atol_k: declaration.temperature_tolerance.value,
            surface_temperature_rtol: 0.0,
            relaxation: declaration.relaxation,
        },
    )?;
    poll(cx, solved.radiation.iterations)?;
    let tolerance = declaration.heat_tolerance.value;
    let applied = &solved.radiation.radiosity;
    let actual = &solved.radiation.equilibrium.nonlinear_radiosity;
    let largest_area = model
        .surfaces()
        .iter()
        .map(RadiationSurface::area_m2)
        .fold(0.0_f64, f64::max);
    let decomposition = solved.conduction.report.energy.robin_out_w
        + solved.conduction.report.energy.neumann_out_w
        - solved.convective_out_w - solved.prescribed_neumann_out_w
        - applied.enclosure_energy_closure_w;
    // The owner already gates the full nonlinear free residual. Additionally
    // enforce the project's explicit patch-watt and closed-energy allowances.
    if !tolerance.is_finite()
        || tolerance <= 0.0
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
            "the converged enclosure failed its declared patch heat, radiosity, closed-energy or convection decomposition gate",
        ));
    }
    // Native common energy always exposes authored Neumann heat separately.
    // Its existing combined convection/radiation channel includes the signed
    // internal enclosure closure; applied per-trace evidence remains retained.
    let energy = &mut solved.conduction.report.energy;
    energy.neumann_out_w = solved.prescribed_neumann_out_w;
    energy.robin_out_w = solved.convective_out_w + applied.enclosure_energy_closure_w;
    energy.closure_w = energy.source_w + energy.dirichlet_in_w
        - energy.neumann_out_w - energy.robin_out_w;
    energy.scale_w = energy.source_w.abs().max(energy.dirichlet_in_w.abs())
        .max(energy.neumann_out_w.abs()).max(energy.robin_out_w.abs()).max(f64::MIN_POSITIVE);
    Ok(SolidSolution {
        conduction: solved.conduction,
        convective_fluxes: Some(solved.convective_robin_fluxes),
        convective_out: Some(solved.convective_out_w),
        radiation: None,
        combined_boundary: Some(solved.combined_boundary),
        enclosure: Some(EnclosureResult {
            model,
            report: EnclosureReport::Steady(solved.radiation),
            decomposition_residual_w: decomposition,
        }),
    })
}

impl EnclosureResult {
    pub(super) fn model(&self) -> &GrayDiffuseEnclosure {
        &self.model
    }

    pub(super) fn from_transient(
        model: GrayDiffuseEnclosure,
        endpoint: &super::EnclosureEndpoint,
    ) -> Self {
        Self {
            model,
            report: EnclosureReport::Transient {
                radiation: endpoint.report.clone(),
                time_integration: endpoint.time_integration,
                dt_s: endpoint.dt_s,
                residual_j: endpoint.physical_residual_norm_j,
                threshold_j: endpoint.physical_residual_tolerance_j,
                energy_residual_j: endpoint.physical_energy_residual_j,
            },
            decomposition_residual_w: endpoint.decomposition_residual_w,
        }
    }

    pub(super) fn receipt(
        &self,
        lowering: &LoweredRadiation,
        convective_out_w: f64,
    ) -> Result<String, SolveRefusal> {
        let declaration = &lowering
            .enclosure
            .as_ref()
            .ok_or_else(|| bad("enclosure receipt lost its declaration"))?
            .declaration;
        let (actual, applied_heat, applied_closure, detail) = match &self.report {
            EnclosureReport::Steady(report) => (
                &report.equilibrium.nonlinear_radiosity,
                report.radiosity.net_outward_heat_w.as_slice(),
                report.radiosity.enclosure_energy_closure_w,
                format!(
                    "\"solid_solves\":{},\"surface_iterations\":{},\"equilibrium_checks\":{},\"max_temperature_change_k\":{},\"max_heat_mismatch_w\":{},\"nonlinear_residual_w\":{},\"nonlinear_residual_threshold_w\":{}",
                    report.iterations + 1, report.iterations, report.equilibrium_checks,
                    num(report.surface_update_history_k.last().copied().unwrap_or(0.0))?,
                    num(report.equilibrium.max_patch_heat_mismatch_w)?,
                    num(report.equilibrium.residual_w)?,
                    num(report.equilibrium.residual_threshold_w)?,
                ),
            ),
            EnclosureReport::Transient {
                radiation, time_integration, dt_s, residual_j, threshold_j, energy_residual_j,
            } => (
                &radiation.radiosity,
                radiation.applied_heat_w.as_slice(),
                radiation.applied_heat_w.iter().sum(),
                format!(
                    "\"time_integration\":{},\"time_step_s\":{},\"solid_solves\":{},\"solid_iterations\":{},\"krylov_iterations\":{},\"max_heat_mismatch_w\":{},\"nonlinear_residual_j\":{},\"nonlinear_residual_threshold_j\":{},\"nonlinear_residual_w\":{},\"nonlinear_residual_threshold_w\":{},\"physical_energy_residual_j\":{}",
                    json_string(time_integration), num(*dt_s)?, radiation.iterations, radiation.solid_iterations,
                    radiation.krylov_iterations, num(radiation.max_heat_mismatch_w)?,
                    num(*residual_j)?, num(*threshold_j)?,
                    num(*residual_j / *dt_s)?, num(*threshold_j / *dt_s)?,
                    num(*energy_residual_j)?,
                ),
            ),
        };
        let mut rows = Vec::with_capacity(self.model.surfaces().len());
        for (i, (surface, source)) in self
            .model
            .surfaces()
            .iter()
            .zip(&lowering.surfaces)
            .enumerate()
        {
            rows.push(format!(
                "{{\"name\":{},\"target\":{},\"card\":{},\"material_state\":{},\"emissivity_receipt\":{},\"emissivity\":{},\"query_temperature_k\":{},\"area_m2\":{},\"mean_temperature_k\":{},\"applied_heat_w\":{},\"nonlinear_heat_w\":{},\"radiosity_w_m2\":{},\"irradiation_w_m2\":{}}}",
                json_string(&source.name), json_string(&source.target),
                json_string(&source.emissivity.card_identity().to_hex()),
                json_string(source.emissivity.material_state()),
                json_string(&source.emissivity.receipt().content_hash().to_hex()),
                num(source.emissivity.value())?, num(source.emissivity.temperature_k())?,
                num(surface.area_m2())?, num(actual.surface_temperatures_k[i])?,
                num(applied_heat[i])?, num(actual.net_outward_heat_w[i])?,
                num(actual.radiosity_w_m2[i])?, num(actual.irradiation_w_m2[i])?,
            ));
        }
        let factors = self
            .model
            .view_factors()
            .factors()
            .iter()
            .map(|row| {
                row.iter()
                    .map(|v| num(*v))
                    .collect::<Result<Vec<_>, _>>()
                    .map(|values| format!("[{}]", values.join(",")))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let evidence = match &declaration.evidence {
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
            "{{\"model\":\"closed-gray-diffuse-enclosure\",\"surfaces\":[{}],\"view_factors\":[{}],\"view_factor_identity\":{},\"view_factor_evidence\":{},\"max_reciprocity_residual\":{},{detail},\"radiosity_residual_w_m2\":{},\"convective_out_w\":{},\"radiative_out_w\":{},\"nonlinear_radiative_out_w\":{},\"decomposition_residual_w\":{},\"controls\":{{\"max_iterations\":{},\"temperature_tolerance_k\":{},\"heat_tolerance_w\":{},\"relaxation\":{},\"closure_tolerance\":{},\"reciprocity_relative_tolerance\":{}}},\"authority\":\"Estimated\",\"no_claim\":{}}}",
            rows.join(","),
            factors.join(","),
            json_string(&self.model.view_factors().identity().to_hex()),
            evidence,
            num(self.model.view_factors().max_reciprocity_residual())?,
            num(actual.linear_residual_max_w_m2)?,
            num(convective_out_w)?,
            num(applied_closure)?,
            num(actual.enclosure_energy_closure_w)?,
            num(self.decomposition_residual_w)?,
            declaration.max_iterations,
            num(declaration.temperature_tolerance.value)?,
            num(declaration.heat_tolerance.value)?,
            num(declaration.relaxation)?,
            num(declaration.closure_tolerance)?,
            num(declaration.reciprocity_relative_tolerance)?,
            json_string(SCOPE),
        ))
    }
}

const SCOPE: &str = "Card-backed opaque gray-diffuse closed-enclosure exchange between exact native P1 boundary traces. User-declared view factors are checked for row closure and actual-area reciprocity on every mesh, not generated or geometrically validated. The shared owner checks the full fresh nonlinear residual; original convection alone exchanges heat with reservoirs or air, prescribed Neumann heat remains external, and mutual radiation is internal. A zero-flux patch declares vacuum without an artificial convection coefficient. Steady initialization still needs a physical conductive thermal anchor; an isolated steady enclosure cannot set its own absolute temperature. Emissivity is frozen at its selected material-card query and retains its source validity. Surface names and matrix order are explicit. Capacity and total-enthalpy endpoints keep the same accepted temperature or enthalpy history in every radiation, natural-convection and air trial and report complete storage/heat residuals in joules; their separate trajectory receipt retains both time grids and all trial work. This receipt describes the accepted primal endpoint; requested steady or final/sampled-peak trajectory derivatives have their own complete-transpose receipt. No continuum error bound, shadowing reconstruction, participating medium, temperature-dependent emissivity or experimental validation is claimed. Solid-solves counts work in this retained enclosure endpoint; nested natural/air iteration counts remain in their own receipts.";
