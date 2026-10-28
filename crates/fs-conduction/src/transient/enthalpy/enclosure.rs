//! Closed gray-diffuse exchange on existing natural enthalpy boundaries.
//!
//! Each patch emits at its P1 area mean and applies a uniform outward flux.
//! The shifted Robin reference `T_air-q_rad/htc` is algebraic, never an air
//! temperature. On a prescribed-flux or insulated trace, radiation adds to the
//! original Neumann load without introducing convection. Every trial solves
//! the same immutable specific-enthalpy
//! history. Acceptance checks the complete nonlinear residual in joules and
//! the independent energy balance; a closed enclosure transfers heat internally.

/// Physical accepted-endpoint tangent, including reflected radiation.
pub mod adjoint;

use fs_exec::Cx;
use fs_solver::{NewtonKrylovState, norm2};

use super::{
    EnthalpyBackwardEuler, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution, StepContext,
    finite, poll,
};
use crate::{
    AmbientRadiationConfig, ConductionError, ConductionProblem, RobinFlux, ScalarField, ThermalBc,
    ThermalBoundary, ThermalInterfaces,
    assemble::ASSEMBLY_TILE,
    radiation::{GrayDiffuseEnclosure, RadiosityReport},
};

/// Accepted exchange and actual work over all enclosure trials.
#[derive(Debug, Clone)]
pub struct EnthalpyEnclosureReport {
    /// Fresh radiosity evaluated at the accepted chart-resolved temperatures.
    pub radiosity: RadiosityReport,
    /// Signed patch watts actually supplied to the accepted inner solve.
    pub applied_heat_w: Vec<f64>,
    /// Maximum applied-versus-fresh patch heat difference, W.
    pub max_heat_mismatch_w: f64,
    /// Number of inner solid solves; none advances history independently.
    pub iterations: usize,
    /// Sum of actual Newton attempts across all inner solves.
    pub solid_iterations: usize,
    /// Sum of actual Newton linear iterations across all inner solves.
    pub krylov_iterations: usize,
}

/// One accepted physical step with separately identified convection and exchange.
#[derive(Debug, Clone)]
pub struct EnthalpyEnclosureStepSolution {
    /// Accepted h/phase state; radiating Robin and Neumann rows include the
    /// applied enclosure flux. Use the separated accounts below externally.
    pub conduction: EnthalpyStepSolution,
    /// Exact shifted boundary used by the accepted inner solve.
    pub combined_boundary: ThermalBoundary,
    /// Original convection only, in the original boundary declaration order.
    pub convective_robin_fluxes: Vec<RobinFlux>,
    /// Original convection heat, W outward; this alone is passed to air.
    pub convective_out_w: f64,
    /// Original prescribed outward flux integrated over all Neumann faces, W.
    /// Excludes internal enclosure exchange, including on vacuum boundaries.
    pub prescribed_neumann_out_w: f64,
    /// Reflected radiation and all actual inner work.
    pub radiation: EnthalpyEnclosureReport,
    /// Complete nonlinear residual at the accepted h field, J.
    pub physical_residual_norm_j: f64,
    /// Original physical initial-residual target, J.
    pub physical_residual_tolerance_j: f64,
    /// Storage minus net input using fresh enclosure exchange, J.
    pub physical_energy_residual_j: f64,
    /// Inward bath reaction using fresh enclosure exchange, W, distinct from
    /// the applied-boundary reaction in `conduction`.
    pub physical_dirichlet_in_w: f64,
}

struct Patch {
    region: usize,
    law: PatchLaw,
}

enum PatchLaw {
    Convection { htc: f64, reference: f64 },
    PrescribedFlux { outward_flux: f64 },
}

impl EnthalpyBackwardEuler<'_, '_> {
    /// Advance closed enclosure exchange on existing natural boundary traces.
    /// Each named radiation surface must exactly match one positive uniform
    /// Robin or uniform Neumann region. Zero Neumann flux declares vacuum;
    /// nonzero authored flux remains additive and is never replaced by exchange.
    /// The existing mesh, source, contact, phase and nonlinear conductivity
    /// owners remain authoritative. At most 64 patches are admitted.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_with_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        step_config: EnthalpyStepConfig,
        enclosure: &GrayDiffuseEnclosure,
        radiation_config: AmbientRadiationConfig,
    ) -> Result<EnthalpyEnclosureStepSolution, EnthalpyError> {
        self.admit_step(cx, problem, old_h, dt_s, step_config)?;
        radiation_config.validate()?;
        radiation_config
            .max_iterations
            .checked_mul(step_config.max_newton_iterations)
            .and_then(|v| v.checked_mul(step_config.newton.linear_restart))
            .and_then(|v| v.checked_mul(step_config.newton.max_linear_cycles))
            .ok_or(EnthalpyError::InvalidInput(
                "enclosure Newton work overflow",
            ))?;
        let context = StepContext {
            storage: self,
            cx,
            problem,
            interfaces,
            old: old_h,
            dt: dt_s,
            pointwise_radiation: None,
        };
        let old_temperature = context.temperatures(old_h)?;
        let patches = bind(cx, problem, enclosure, &old_temperature)?;
        let prescribed_neumann_out_w = prescribed_neumann_heat(cx, problem)?;
        let target = initial_target(&context, step_config, enclosure, &patches, radiation_config)?;
        let mut driving = means(cx, problem, enclosure, &old_temperature)?;
        let mut solid_iterations = 0_usize;
        let mut krylov_iterations = 0_usize;
        let mut last_change = 0.0_f64;
        let mut last_mismatch = 0.0_f64;
        for iteration in 0..radiation_config.max_iterations {
            poll(cx, iteration)?;
            let applied = radiosity(cx, enclosure, &driving, radiation_config)?;
            let boundary = shifted_boundary(problem, &patches, &applied)?;
            let trial_context = StepContext {
                problem: ConductionProblem {
                    boundary: &boundary,
                    ..problem
                },
                ..context
            };
            let mut initial = vec![0.0; old_h.len()];
            trial_context.residual(&trial_context.initial_values()?, &mut initial)?;
            let norm = finite(norm2(&initial))?;
            let mut inner = step_config;
            inner.newton.absolute_tolerance = target * 0.25;
            if target != 0.0 || norm != 0.0 {
                inner.newton.relative_tolerance = inner
                    .newton
                    .relative_tolerance
                    .min((target * 0.25) / norm.max(f64::MIN_POSITIVE));
            }
            if inner.newton.relative_tolerance <= 0.0 {
                return Err(EnthalpyError::InvalidInput(
                    "enclosure residual tolerance is not representable",
                ));
            }
            let solid = self.advance(cx, trial_context.problem, interfaces, old_h, dt_s, inner)?;
            solid_iterations = solid_iterations
                .checked_add(solid.newton.iterations)
                .ok_or(EnthalpyError::InvalidInput(
                    "enclosure Newton work overflow",
                ))?;
            for row in &solid.newton.history {
                krylov_iterations = krylov_iterations.checked_add(row.linear_iterations).ok_or(
                    EnthalpyError::InvalidInput("enclosure Krylov work overflow"),
                )?;
            }
            let actual = means(cx, problem, enclosure, &solid.temperature)?;
            let fresh = radiosity(cx, enclosure, &actual, radiation_config)?;
            let (convection, convective_out) = convection(cx, problem, &solid)?;
            let mut heat_ok = true;
            last_mismatch = 0.0;
            for (index, patch) in patches.iter().enumerate() {
                poll(cx, index)?;
                let actual_applied = match patch.law {
                    PatchLaw::Convection { .. } => {
                        let name = &problem.boundary.region_names()[patch.region];
                        let full = solid.robin_fluxes.iter().find(|row| row.region == *name)
                            .ok_or(EnthalpyError::InvalidInput(
                                "enclosure inner solve omitted a convective patch",
                            ))?;
                        let air = convection.iter().find(|row| row.region == *name)
                            .ok_or(EnthalpyError::InvalidInput(
                                "enclosure convection omitted a patch",
                            ))?;
                        finite(full.heat_rate_w - air.heat_rate_w)?
                    }
                    PatchLaw::PrescribedFlux { outward_flux } => {
                        let ThermalBc::Neumann { outward_flux: ScalarField::Uniform(combined) }
                            = boundary.conditions()[patch.region]
                        else {
                            return Err(EnthalpyError::InvalidInput(
                                "enclosure changed prescribed-flux boundary ownership",
                            ));
                        };
                        finite((combined - outward_flux) * enclosure.surfaces()[index].area_m2())?
                    }
                };
                if finite(actual_applied - applied.net_outward_heat_w[index])?.abs()
                    > heat_tolerance(
                        radiation_config,
                        actual_applied,
                        applied.net_outward_heat_w[index],
                    )?
                {
                    return Err(EnthalpyError::InvalidInput(
                        "enclosure boundary heat does not separate into authored load and radiation",
                    ));
                }
                let mismatch = finite(actual_applied - fresh.net_outward_heat_w[index])?.abs();
                last_mismatch = last_mismatch.max(mismatch);
                heat_ok &= mismatch
                    <= heat_tolerance(
                        radiation_config,
                        actual_applied,
                        fresh.net_outward_heat_w[index],
                    )?;
            }
            last_change = actual
                .iter()
                .zip(&driving)
                .try_fold(0.0_f64, |m, (&a, &b)| {
                    Ok::<_, EnthalpyError>(m.max(finite(a - b)?.abs()))
                })?;
            let physical_boundary = shifted_boundary(problem, &patches, &fresh)?;
            let physical = StepContext {
                problem: ConductionProblem {
                    boundary: &physical_boundary,
                    ..problem
                },
                ..context
            };
            let mut residual = vec![0.0; old_h.len()];
            physical.residual(&solid.specific_enthalpy_j_kg, &mut residual)?;
            let residual_norm = finite(norm2(&residual))?;
            let physical_dirichlet_in_w =
                physical.physical_prescribed_heat_input_w(&solid.specific_enthalpy_j_kg)?;
            let energy = finite(
                solid.stored_energy_change_j
                    + dt_s
                        * (prescribed_neumann_out_w + convective_out + fresh.enclosure_energy_closure_w
                            - solid.source_w - physical_dirichlet_in_w),
            )?;
            if heat_ok
                && last_change <= radiation_config.temperature_tolerance_k
                && residual_norm <= target
                && energy.abs() <= step_config.energy_tolerance_j
            {
                poll(cx, iteration)?;
                return Ok(EnthalpyEnclosureStepSolution {
                    conduction: solid,
                    combined_boundary: boundary,
                    convective_robin_fluxes: convection,
                    convective_out_w: convective_out,
                    prescribed_neumann_out_w,
                    radiation: EnthalpyEnclosureReport {
                        radiosity: fresh,
                        applied_heat_w: applied.net_outward_heat_w,
                        max_heat_mismatch_w: last_mismatch,
                        iterations: iteration + 1,
                        solid_iterations,
                        krylov_iterations,
                    },
                    physical_residual_norm_j: residual_norm,
                    physical_residual_tolerance_j: target,
                    physical_energy_residual_j: energy,
                    physical_dirichlet_in_w,
                });
            }
            for (old, new) in driving.iter_mut().zip(actual) {
                *old = finite(
                    (1.0 - radiation_config.relaxation) * *old + radiation_config.relaxation * new,
                )?;
            }
        }
        Err(ConductionError::Radiation {surface:"enthalpy-enclosure".into(),
            what:format!("{} solid solves exhausted; final temperature change {last_change} K and patch mismatch {last_mismatch} W; physical residual or energy may also remain unresolved",radiation_config.max_iterations),
            fix:"provide adequate explicit enclosure/Newton budgets or reduce the timestep".into()}.into())
    }
}

fn bind(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    enclosure: &GrayDiffuseEnclosure,
    temperature: &[f64],
) -> Result<Vec<Patch>, EnthalpyError> {
    if enclosure.surfaces().is_empty() || enclosure.surfaces().len() > 64 {
        return Err(EnthalpyError::InvalidInput(
            "enthalpy enclosure requires 1..=64 patches",
        ));
    }
    let mut patches = Vec::with_capacity(enclosure.surfaces().len());
    for surface in enclosure.surfaces() {
        poll(cx, patches.len())?;
        // Public trace query checks the exact mesh identity before slot access.
        surface.mean_temperature(problem.mesh, temperature)?;
        let region = problem
            .boundary
            .region_names()
            .iter()
            .position(|name| name == surface.name())
            .ok_or(EnthalpyError::InvalidInput(
                "unknown enthalpy enclosure region",
            ))?;
        let law = match problem.boundary.conditions()[region] {
            ThermalBc::Robin {
                htc: ScalarField::Uniform(htc),
                t_ref: ScalarField::Uniform(reference),
            } if htc > 0.0 => PatchLaw::Convection { htc, reference },
            ThermalBc::Neumann { outward_flux: ScalarField::Uniform(outward_flux) } => {
                PatchLaw::PrescribedFlux { outward_flux }
            }
            _ => {
                return Err(EnthalpyError::InvalidInput(
                    "enclosure requires uniform Robin or Neumann patches",
                ));
            }
        };
        let mut count = 0_usize;
        for slot in 0..problem.mesh.boundary().len() {
            if slot % ASSEMBLY_TILE == 0 {
                poll(cx, slot)?;
            }
            if problem.boundary.region_for(slot) == Some(region) {
                count += 1;
            }
        }
        if count != surface.face_slots().len()
            || surface
                .face_slots()
                .iter()
                .any(|&slot| problem.boundary.region_for(slot) != Some(region))
        {
            return Err(EnthalpyError::InvalidInput(
                "enclosure patch must exactly cover its named natural region",
            ));
        }
        patches.push(Patch { region, law });
    }
    Ok(patches)
}

fn means(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    enclosure: &GrayDiffuseEnclosure,
    temperature: &[f64],
) -> Result<Vec<f64>, EnthalpyError> {
    enclosure
        .surfaces()
        .iter()
        .enumerate()
        .map(|(i, surface)| {
            poll(cx, i)?;
            Ok(surface.mean_temperature(problem.mesh, temperature)?)
        })
        .collect()
}

fn heat_tolerance(config: AmbientRadiationConfig, a: f64, b: f64) -> Result<f64, EnthalpyError> {
    finite(
        config
            .balance_relative_tolerance
            .mul_add(a.abs().max(b.abs()), config.balance_tolerance_w),
    )
}

fn radiosity(
    cx: &Cx<'_>,
    enclosure: &GrayDiffuseEnclosure,
    temperature: &[f64],
    config: AmbientRadiationConfig,
) -> Result<RadiosityReport, EnthalpyError> {
    poll(cx, 0)?;
    let report = enclosure.solve(cx, temperature)?;
    for &value in report.radiosity_w_m2.iter().chain(&report.irradiation_w_m2) {
        if !value.is_finite() || value < 0.0 {
            return Err(EnthalpyError::InvalidInput(
                "nonphysical enclosure radiosity",
            ));
        }
    }
    for &value in report
        .net_outward_flux_w_m2
        .iter()
        .chain(&report.net_outward_heat_w)
    {
        finite(value)?;
    }
    let area = enclosure
        .surfaces()
        .iter()
        .map(|surface| surface.area_m2())
        .fold(0.0_f64, f64::max);
    let tolerance = heat_tolerance(config, report.enclosure_energy_scale_w, 0.0)?;
    if finite(report.linear_residual_max_w_m2 * area)? > tolerance
        || finite(report.enclosure_energy_closure_w)?.abs() > tolerance
    {
        return Err(EnthalpyError::InvalidInput(
            "enclosure radiosity equations or closed energy sum exceed their allowance",
        ));
    }
    Ok(report)
}

fn shifted_boundary(
    problem: ConductionProblem<'_>,
    patches: &[Patch],
    radiosity: &RadiosityReport,
) -> Result<ThermalBoundary, EnthalpyError> {
    let mut robin = Vec::new();
    let mut neumann = Vec::new();
    for (patch, &flux) in patches.iter().zip(&radiosity.net_outward_flux_w_m2) {
        match patch.law {
            PatchLaw::Convection { htc, reference } => {
                robin.push((patch.region, htc, finite(reference - flux / htc)?));
            }
            PatchLaw::PrescribedFlux { .. } => neumann.push((patch.region, flux)),
        }
    }
    Ok(problem.boundary.with_uniform_robin_replacements(&robin)?
        .with_uniform_neumann_additions(&neumann)?)
}

// Integrate the authored boundary independently. A mixed Robin/Neumann
// enclosure generally transfers nonzero net radiative heat through either
// subset, so the inner solve's Neumann watts cannot be treated as external.
fn prescribed_neumann_heat(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
) -> Result<f64, EnthalpyError> {
    let mut heat = 0.0;
    for (slot, face) in problem.mesh.boundary().iter().enumerate() {
        if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
        if let Some(ThermalBc::Neumann { outward_flux }) = problem.boundary.condition_for(slot) {
            let mean = face.vertices.iter().map(|&v| outward_flux.at(v as usize)).sum::<f64>() / 3.0;
            heat = finite(heat + face.area * mean)?;
        }
    }
    poll(cx, problem.mesh.boundary().len())?;
    Ok(heat)
}

fn convection(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    solid: &EnthalpyStepSolution,
) -> Result<(Vec<RobinFlux>, f64), EnthalpyError> {
    let mut total = 0.0;
    let mut rows = Vec::with_capacity(solid.robin_fluxes.len());
    for (ordinal, row) in solid.robin_fluxes.iter().enumerate() {
        poll(cx, ordinal)?;
        let region = problem
            .boundary
            .region_names()
            .iter()
            .position(|name| *name == row.region)
            .ok_or(EnthalpyError::InvalidInput(
                "unknown original convection region",
            ))?;
        let mut original = row.clone();
        match &problem.boundary.conditions()[region] {
            ThermalBc::Robin {
                htc: ScalarField::Uniform(h),
                t_ref: ScalarField::Uniform(reference),
            } => {
                original.mean_htc_w_per_m2_k = *h;
                original.mean_reference_temperature_k = *reference;
                original.heat_rate_w =
                    finite(h * row.area_m2 * (row.mean_wall_temperature_k - reference))?;
            }
            // Only radiating rows are changed. Unrelated nonuniform Robin
            // rows retain their exact inner-solve integrals.
            ThermalBc::Robin { .. } => {}
            _ => {
                return Err(EnthalpyError::InvalidInput(
                    "enclosure report contains a non-Robin row",
                ));
            }
        }
        total = finite(total + original.heat_rate_w)?;
        rows.push(original);
    }
    Ok((rows, total))
}

fn initial_target(
    context: &StepContext<'_, '_, '_, '_>,
    config: EnthalpyStepConfig,
    enclosure: &GrayDiffuseEnclosure,
    patches: &[Patch],
    radiation_config: AmbientRadiationConfig,
) -> Result<f64, EnthalpyError> {
    let initial_h = context.initial_values()?;
    let temperature = context.temperatures(&initial_h)?;
    let means = means(context.cx, context.problem, enclosure, &temperature)?;
    let initial = radiosity(context.cx, enclosure, &means, radiation_config)?;
    let boundary = shifted_boundary(context.problem, patches, &initial)?;
    let physical = StepContext {
        problem: ConductionProblem {
            boundary: &boundary,
            ..context.problem
        },
        storage: context.storage,
        cx: context.cx,
        interfaces: context.interfaces,
        old: context.old,
        dt: context.dt,
        pointwise_radiation: None,
    };
    let stage = physical.stage(&initial_h)?;
    let initial = NewtonKrylovState::new(&stage, initial_h, config.newton)
        .map_err(|error| stage.take_failure().unwrap_or(EnthalpyError::Newton(error)))?;
    finite(
        config
            .newton
            .absolute_tolerance
            .max(config.newton.relative_tolerance * initial.residual_norm().max(f64::MIN_POSITIVE)),
    )
}
