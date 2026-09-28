//! Ambient gray radiation on a named, uniform convective P1 trace.
//!
//! This preserves the cooling lab's area-mean patch model: its integrated
//! radiative loss is epsilon sigma A (T_mean^4 - T_ambient^4). A positive
//! secant coefficient multiplies the pointwise Robin trace during each solid
//! solve. It is neither an integral of pointwise T(x)^4 nor uniform patch
//! flux. The reservoir is black, isothermal and occupies view factor one;
//! there is no enclosure, occlusion or participating medium.
//!
//! The returned conduction report retains the combined Robin operator that
//! was actually solved. Separate ORIGINAL convection rows carry heat to air.
//! Fixed contact, heterogeneous and nonlinear conductivity stay in every
//! inner solve. No derivative of this radiative fixed point is provided.

use super::{
    STEFAN_BOLTZMANN_W_M2_K4, SurfaceEmissivity, radiation_error, require_temperature,
    run_conduction,
};
use crate::{
    ConductionError, ConductionProblem, ConductionSolution, DofMap, InitialGuess, RobinFlux,
    ScalarField, SolveConfig, ThermalBc, ThermalBoundary, ThermalInterfaces,
};
use fs_exec::Cx;
use std::collections::BTreeSet;

/// One immutable material-card-backed ambient patch on an existing Robin row.
#[derive(Debug, Clone, PartialEq)]
pub struct AmbientRadiationPatch {
    region: String,
    emissivity: SurfaceEmissivity,
    ambient_temperature_k: f64,
}

impl AmbientRadiationPatch {
    /// Bind the material declaration and reservoir temperature. The named
    /// Robin trace is admitted against the actual problem at solve time.
    pub fn new(
        region: impl Into<String>,
        emissivity: SurfaceEmissivity,
        ambient_temperature_k: f64,
    ) -> Result<Self, ConductionError> {
        let region = region.into();
        if region.trim().is_empty() {
            return Err(radiation_error(
                "<unnamed>",
                "ambient patch name is blank",
                "name the existing convective boundary region",
            ));
        }
        // The reservoir need not lie within the SURFACE material's validity
        // interval. Only actual and driving patch temperatures query that law.
        require_temperature(&region, "radiative ambient", ambient_temperature_k)?;
        let patch = Self {
            region,
            emissivity,
            ambient_temperature_k,
        };
        patch.secant_coefficient_w_m2_k(patch.emissivity.temperature_k())?;
        Ok(patch)
    }

    /// Named existing uniform Robin region.
    #[must_use]
    pub fn region(&self) -> &str {
        &self.region
    }

    /// Surface material declaration and exact property-use receipt.
    #[must_use]
    pub const fn emissivity(&self) -> &SurfaceEmissivity {
        &self.emissivity
    }

    /// Black reservoir absolute temperature, K.
    #[must_use]
    pub const fn ambient_temperature_k(&self) -> f64 {
        self.ambient_temperature_k
    }

    /// Positive, stable secant of epsilon sigma (T^4 - T_ambient^4).
    /// The factorization remains defined at T = T_ambient and avoids
    /// subtracting nearly equal fourth powers.
    pub fn secant_coefficient_w_m2_k(
        &self,
        surface_temperature_k: f64,
    ) -> Result<f64, ConductionError> {
        self.emissivity.validate_temperature(
            &self.region,
            "ambient-radiation surface temperature",
            surface_temperature_k,
        )?;
        let t = surface_temperature_k;
        let a = self.ambient_temperature_k;
        let coefficient =
            self.emissivity.value() * STEFAN_BOLTZMANN_W_M2_K4 * (t + a) * t.mul_add(t, a * a);
        if coefficient.is_finite() && coefficient > 0.0 {
            Ok(coefficient)
        } else {
            Err(radiation_error(
                &self.region,
                "radiative secant is not finite and positive",
                "use a supported finite temperature and emissivity range",
            ))
        }
    }

    /// Nonlinear patch heat-flux density, W/m2; negative heats the solid.
    pub fn heat_flux_w_m2(&self, surface_temperature_k: f64) -> Result<f64, ConductionError> {
        finite(
            self.secant_coefficient_w_m2_k(surface_temperature_k)?
                * (surface_temperature_k - self.ambient_temperature_k),
            "radiative heat flux",
        )
    }
}

/// Bounded outer iteration, with independent raw temperature and heat gates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AmbientRadiationConfig {
    /// Maximum solid solves, including the first radiative secant iterate.
    pub max_iterations: usize,
    /// Maximum unrelaxed patch-temperature difference, K.
    pub temperature_tolerance_k: f64,
    /// Absolute heat mismatch allowance per patch and for the heat split, W.
    pub balance_tolerance_w: f64,
    /// Relative allowance against the larger applied/nonlinear patch heat.
    pub balance_relative_tolerance: f64,
    /// Relaxation of the next driving temperature, in (0, 1].
    pub relaxation: f64,
}

impl Default for AmbientRadiationConfig {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            temperature_tolerance_k: 1e-9,
            balance_tolerance_w: 1e-8,
            balance_relative_tolerance: 1e-10,
            relaxation: 0.5,
        }
    }
}

impl AmbientRadiationConfig {
    fn validate(self) -> Result<(), ConductionError> {
        if self.max_iterations == 0
            || !(self.temperature_tolerance_k.is_finite() && self.temperature_tolerance_k > 0.0)
            || !(self.balance_tolerance_w.is_finite() && self.balance_tolerance_w >= 0.0)
            || !(self.balance_relative_tolerance.is_finite()
                && self.balance_relative_tolerance >= 0.0)
            || self.balance_tolerance_w == 0.0 && self.balance_relative_tolerance == 0.0
            || !(self.relaxation.is_finite() && self.relaxation > 0.0 && self.relaxation <= 1.0)
        {
            return Err(radiation_error(
                "ambient-coupling",
                "inadmissible iteration controls",
                "declare positive iterations and temperature tolerance, a positive watt or relative allowance, and relaxation in (0, 1]",
            ));
        }
        Ok(())
    }

    fn heat_tolerance(self, a: f64, b: f64) -> Result<f64, ConductionError> {
        finite(
            self.balance_relative_tolerance
                .mul_add(a.abs().max(b.abs()), self.balance_tolerance_w),
            "radiative heat tolerance",
        )
    }
}

/// Heat and constitutive evidence for one accepted patch, in input order.
#[derive(Debug, Clone, PartialEq)]
pub struct AmbientRadiationPatchReport {
    /// Exact physical declaration and its emissivity receipt.
    pub patch: AmbientRadiationPatch,
    /// Integrated trace area, m2.
    pub area_m2: f64,
    /// Solved area-mean surface temperature, K.
    pub mean_surface_temperature_k: f64,
    /// Temperature whose secant was actually assembled, K.
    pub driving_temperature_k: f64,
    /// Radiative secant coefficient actually assembled, W/(m2 K).
    pub applied_coefficient_w_m2_k: f64,
    /// Radiative heat carried by the assembled Robin row, W outward.
    pub applied_heat_w: f64,
    /// epsilon sigma A (T_mean^4 - T_ambient^4), W outward.
    pub nonlinear_heat_w: f64,
    /// Absolute difference of nonlinear and applied heat, W.
    pub heat_mismatch_w: f64,
    /// Independent convergence threshold applied to that difference, W.
    pub heat_tolerance_w: f64,
}

/// Accepted outer fixed point; this is an Estimated discrete patch model.
#[derive(Debug, Clone, PartialEq)]
pub struct AmbientRadiationReport {
    /// Number of solid solves actually performed.
    pub iterations: usize,
    /// Sum of nonlinear solid iterations over all outer iterations.
    pub solid_iterations: usize,
    /// Sum of actual inner Krylov work over all outer iterations.
    pub krylov_iterations: usize,
    /// Maximum raw surface-temperature update on the accepted iterate, K.
    pub max_temperature_change_k: f64,
    /// Maximum patch applied/nonlinear heat mismatch, W.
    pub max_heat_mismatch_w: f64,
    /// Exact bounded controls used for this solve.
    pub config: AmbientRadiationConfig,
    /// Each patch in declaration order, with constitutive evidence.
    pub patches: Vec<AmbientRadiationPatchReport>,
    /// Sum of applied radiative heat, W outward; includes hot reservoirs.
    pub applied_radiation_out_w: f64,
    /// Sum of the accepted nonlinear patch heat, W outward.
    pub nonlinear_radiation_out_w: f64,
    /// Combined Robin total minus convection minus applied radiation, W.
    pub decomposition_residual_w: f64,
}

/// Solved combined boundary physics and separate ORIGINAL convection records.
#[derive(Debug, Clone, PartialEq)]
pub struct AmbientRadiationSolution {
    /// Actual solved field/report. Its Robin rows and energy include both
    /// convection and applied radiation; do not feed them directly into air.
    pub conduction: ConductionSolution,
    /// Exact combined partition/coefficient/reference that produced the field.
    pub combined_boundary: ThermalBoundary,
    /// Exact initial state, nonlinear stop and linear controls of the final
    /// solid solve; a derivative consumer can reconstruct the same binding.
    /// This does not itself include a radiative derivative.
    pub final_solve_config: SolveConfig,
    /// Heat rates on the ORIGINAL convection laws, in declaration order.
    pub convective_robin_fluxes: Vec<RobinFlux>,
    /// Sum of original convective heat rates only, W outward.
    pub convective_out_w: f64,
    /// Nonlinear radiative heat, convergence and cumulative work evidence.
    pub radiation: AmbientRadiationReport,
}

struct BoundPatch<'a> {
    patch: &'a AmbientRadiationPatch,
    region: usize,
    convective_coefficient: f64,
    convective_reference: f64,
}

/// Solve the existing conduction/contact problem with ambient patch radiation
/// in addition to its named uniform Robin laws. Every iteration uses the
/// original heterogeneous/nonlinear conductivity and finite contact operator.
///
/// Patches must be unique nonempty existing uniform Robin traces. Failed
/// material validity, incomplete convergence, exhausted work and cancellation
/// return errors, never a frozen-coefficient success. This function does not
/// provide radiative derivatives or continuum/uncertainty bounds.
pub fn solve_with_ambient_radiation(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    patches: &[AmbientRadiationPatch],
    conduction_config: SolveConfig,
    config: AmbientRadiationConfig,
) -> Result<AmbientRadiationSolution, ConductionError> {
    poll(cx, 0)?;
    config.validate()?;
    if patches.is_empty() {
        return Err(radiation_error(
            "ambient-coupling",
            "no radiation patches were declared",
            "use the ordinary conduction solve when radiation is absent",
        ));
    }
    let mut bound = Vec::with_capacity(patches.len());
    let mut names = BTreeSet::new();
    for patch in patches {
        poll(cx, bound.len())?;
        if !names.insert(patch.region()) {
            return Err(radiation_error(
                patch.region(),
                "duplicate ambient patch",
                "declare one radiation owner per boundary region",
            ));
        }
        let region = problem
            .boundary
            .region_names()
            .iter()
            .position(|name| name == patch.region())
            .ok_or_else(|| {
                radiation_error(
                    patch.region(),
                    "no matching boundary region",
                    "name an existing uniform convective boundary",
                )
            })?;
        let ThermalBc::Robin {
            htc: ScalarField::Uniform(h),
            t_ref: ScalarField::Uniform(reference),
        } = &problem.boundary.conditions()[region]
        else {
            return Err(radiation_error(
                patch.region(),
                "ambient patch is not a uniform Robin boundary",
                "bind the patch to uniform convection; prescribed, flux and nodal Robin rows are unsupported",
            ));
        };
        require_temperature(patch.region(), "convective reference", *reference)?;
        let mut area = 0.0;
        for (slot, face) in problem.mesh.boundary().iter().enumerate() {
            poll(cx, slot)?;
            if problem.boundary.region_for(slot) == Some(region) {
                area = finite(area + face.area, "radiative patch area")?;
            }
        }
        if area <= 0.0 {
            return Err(radiation_error(
                patch.region(),
                "ambient patch owns no boundary area",
                "bind radiation to a nonempty exterior trace",
            ));
        }
        bound.push(BoundPatch {
            patch,
            region,
            convective_coefficient: *h,
            convective_reference: *reference,
        });
    }
    let dofs = DofMap::new(problem.boundary, problem.mesh.vertex_count())?;
    let mut driving: Vec<f64> = patches
        .iter()
        .map(|patch| patch.emissivity.temperature_k())
        .collect();
    let mut next_config = conduction_config.clone();
    let mut solid_iterations = 0_usize;
    let mut krylov_iterations = 0_usize;
    let mut last_change = f64::INFINITY;
    let mut last_mismatch = f64::INFINITY;
    for iteration in 0..config.max_iterations {
        poll(cx, iteration)?;
        let mut applied = Vec::with_capacity(bound.len());
        let mut replacements = Vec::with_capacity(bound.len());
        for (patch, &temperature) in bound.iter().zip(&driving) {
            let h_rad = patch.patch.secant_coefficient_w_m2_k(temperature)?;
            let h = finite(patch.convective_coefficient + h_rad, "combined coefficient")?;
            let reference = finite(
                (patch.convective_coefficient / h) * patch.convective_reference
                    + (h_rad / h) * patch.patch.ambient_temperature_k,
                "combined reference",
            )?;
            replacements.push((patch.region, h, reference));
            applied.push(h_rad);
        }
        let combined_boundary = problem
            .boundary
            .with_uniform_robin_replacements(&replacements)?;
        let combined_problem = ConductionProblem {
            boundary: &combined_boundary,
            ..problem
        };
        let final_solve_config = next_config.clone();
        let conduction = run_conduction(cx, combined_problem, interfaces, next_config)?;
        poll(cx, iteration)?;
        if conduction.report.final_residual > conduction.report.residual_threshold {
            return Err(radiation_error(
                "ambient-coupling",
                "inner solid stopped without satisfying its residual gate",
                "tighten the solid stop rule and disable premature step-only convergence",
            ));
        }
        solid_iterations = count(solid_iterations, conduction.report.iterations)?;
        for solve in &conduction.report.linear {
            krylov_iterations = count(krylov_iterations, solve.iterations)?;
        }
        let mut convective_robin_fluxes = conduction.report.robin_fluxes.clone();
        let mut rows = Vec::with_capacity(bound.len());
        let mut temperature_converged = true;
        let mut heat_converged = true;
        last_change = 0.0;
        last_mismatch = 0.0;
        let mut applied_radiation_out_w = 0.0;
        let mut nonlinear_radiation_out_w = 0.0;
        for (index, patch) in bound.iter().enumerate() {
            poll(cx, index)?;
            let flux = convective_robin_fluxes
                .iter_mut()
                .find(|flux| flux.region == patch.patch.region)
                .ok_or_else(|| {
                    radiation_error(
                        patch.patch.region(),
                        "solved radiation trace is missing",
                        "report the boundary binding defect",
                    )
                })?;
            let mean = flux.mean_wall_temperature_k;
            let nonlinear_heat_w = finite(
                patch.patch.heat_flux_w_m2(mean)? * flux.area_m2,
                "nonlinear radiative heat",
            )?;
            let applied_heat_w = finite(
                applied[index] * flux.area_m2 * (mean - patch.patch.ambient_temperature_k),
                "applied radiative heat",
            )?;
            let convective_heat = finite(
                patch.convective_coefficient * flux.area_m2 * (mean - patch.convective_reference),
                "convective heat",
            )?;
            let split_error = finite(
                flux.heat_rate_w - convective_heat - applied_heat_w,
                "patch heat decomposition",
            )?;
            if split_error.abs()
                > config.heat_tolerance(
                    flux.heat_rate_w,
                    convective_heat.abs() + applied_heat_w.abs(),
                )?
            {
                return Err(radiation_error(
                    patch.patch.region(),
                    format!("combined Robin heat split differs by {split_error} W"),
                    "use resolvable heat tolerances or report the boundary decomposition defect",
                ));
            }
            let change = finite(mean - driving[index], "patch temperature update")?.abs();
            let mismatch =
                finite(nonlinear_heat_w - applied_heat_w, "radiative heat mismatch")?.abs();
            let tolerance = config.heat_tolerance(nonlinear_heat_w, applied_heat_w)?;
            last_change = last_change.max(change);
            last_mismatch = last_mismatch.max(mismatch);
            temperature_converged &= change <= config.temperature_tolerance_k;
            heat_converged &= mismatch <= tolerance;
            applied_radiation_out_w = finite(
                applied_radiation_out_w + applied_heat_w,
                "applied radiative total",
            )?;
            nonlinear_radiation_out_w = finite(
                nonlinear_radiation_out_w + nonlinear_heat_w,
                "nonlinear radiative total",
            )?;
            rows.push(AmbientRadiationPatchReport {
                patch: patch.patch.clone(),
                area_m2: flux.area_m2,
                mean_surface_temperature_k: mean,
                driving_temperature_k: driving[index],
                applied_coefficient_w_m2_k: applied[index],
                applied_heat_w,
                nonlinear_heat_w,
                heat_mismatch_w: mismatch,
                heat_tolerance_w: tolerance,
            });
            flux.mean_htc_w_per_m2_k = patch.convective_coefficient;
            flux.mean_reference_temperature_k = patch.convective_reference;
            flux.heat_rate_w = convective_heat;
        }
        let mut convective_out_w = 0.0;
        for flux in &convective_robin_fluxes {
            convective_out_w = finite(convective_out_w + flux.heat_rate_w, "convective total")?;
        }
        let decomposition_residual_w = finite(
            conduction.report.energy.robin_out_w - convective_out_w - applied_radiation_out_w,
            "boundary heat decomposition",
        )?;
        if decomposition_residual_w.abs()
            > config.heat_tolerance(
                conduction.report.energy.robin_out_w,
                convective_out_w.abs() + applied_radiation_out_w.abs(),
            )?
        {
            return Err(radiation_error(
                "ambient-coupling",
                "whole-boundary heat decomposition failed",
                "use resolvable heat tolerances or report the boundary decomposition defect",
            ));
        }
        if temperature_converged && heat_converged {
            poll(cx, iteration)?;
            return Ok(AmbientRadiationSolution {
                conduction,
                combined_boundary,
                final_solve_config,
                convective_robin_fluxes,
                convective_out_w,
                radiation: AmbientRadiationReport {
                    iterations: iteration + 1,
                    solid_iterations,
                    krylov_iterations,
                    max_temperature_change_k: last_change,
                    max_heat_mismatch_w: last_mismatch,
                    config,
                    patches: rows,
                    applied_radiation_out_w,
                    nonlinear_radiation_out_w,
                    decomposition_residual_w,
                },
            });
        }
        for (value, row) in driving.iter_mut().zip(&rows) {
            *value = finite(
                config.relaxation.mul_add(
                    row.mean_surface_temperature_k,
                    (1.0 - config.relaxation) * *value,
                ),
                "relaxed patch temperature",
            )?;
        }
        next_config = conduction_config.clone();
        next_config.initial = InitialGuess::Free(dofs.gather(&conduction.temperature));
    }
    Err(ConductionError::AmbientRadiationNotConverged {
        iterations: config.max_iterations,
        temperature_change_k: last_change,
        temperature_tolerance_k: config.temperature_tolerance_k,
        heat_mismatch_w: last_mismatch,
    })
}

fn finite(value: f64, what: &str) -> Result<f64, ConductionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(radiation_error(
            "ambient-coupling",
            format!("non-finite {what}"),
            "use a supported finite temperature and heat range",
        ))
    }
}

fn count(total: usize, increment: usize) -> Result<usize, ConductionError> {
    total.checked_add(increment).ok_or_else(|| {
        radiation_error(
            "ambient-coupling",
            "cumulative solver work overflow",
            "reduce the declared solve budgets",
        )
    })
}

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "ambient-radiation",
        at,
    })
}
