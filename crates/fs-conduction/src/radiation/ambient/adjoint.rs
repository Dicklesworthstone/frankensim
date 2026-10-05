//! Total derivatives of the SAME mean-temperature secant Robin model as the
//! ambient-radiation producer. In particular this is not uniform patch flux
//! and is not the integral of pointwise T(x)^4 over a non-isothermal wall.

use super::{AmbientRadiationPatch, ConductionError, ConductionProblem, Cx,
    ScalarField, ThermalBc, ThermalInterfaces, finite, poll, radiation_error,
    STEFAN_BOLTZMANN_W_M2_K4};
use crate::LinearConfig;
use crate::adjoint::{RobinGradient, RobinResponse};
use std::collections::BTreeSet;

/// Discrete nodal-temperature sensitivities through ambient radiation.
/// Every ordering is explicit; no physical or uncertainty certification is
/// inferred from the residual-checked numerical derivative.
#[derive(Debug, Clone, PartialEq)]
pub struct AmbientRadiationGradient {
    /// Original convective Robin regions, in the caller's selected order.
    pub regions: Vec<String>,
    /// Derivatives with respect to ORIGINAL convective references, ln(h_c),
    /// and assembled nodal loads. These include radiative state feedback;
    /// they are not partials of the artificially combined Robin row.
    pub convection: RobinGradient,
    /// Radiating regions, in the input patch order.
    pub radiation_regions: Vec<String>,
    /// Derivative with respect to each black-reservoir temperature, per K.
    /// Includes both the reservoir reference and its secant-coefficient term.
    pub reservoir_temperatures: Vec<f64>,
    /// Derivative per unit absolute hemispherical emissivity. The selected
    /// card value is a control; its query temperature, identity and validity
    /// domain are fixed. At parameter limits only feasible directions apply.
    pub emissivities: Vec<f64>,
}

fn invalid(what: &str) -> ConductionError {
    radiation_error("ambient-adjoint", what,
        "use the retained field, original uniform convection rows, exact radiation patches, and explicit derivative budgets")
}

impl AmbientRadiationPatch {
    /// Analytic secant partials [d h_rad / d wall mean, d h_rad / d reservoir]
    /// in W/(m^2 K^2). Emissivity remains the fixed card-query value.
    ///
    /// Differentiates epsilon*sigma*(S+A)*(S^2+A^2), including its finite
    /// value at S=A. The wall must lie inside the surface-material validity
    /// interval; no two-sided derivative is supplied at an interval endpoint.
    pub fn secant_partials_w_m2_k2(&self, wall_mean_k: f64)
        -> Result<[f64; 2], ConductionError>
    {
        self.secant_coefficient_w_m2_k(wall_mean_k)?;
        self.secant_coefficient_w_m2_k(wall_mean_k.next_down())?;
        self.secant_coefficient_w_m2_k(wall_mean_k.next_up())?;
        let s = wall_mean_k;
        let a = self.ambient_temperature_k();
        let scale = self.emissivity().value() * STEFAN_BOLTZMANN_W_M2_K4;
        Ok([
            finite(scale * (3.0*s*s + 2.0*s*a + a*a), "radiative wall partial")?,
            finite(scale * (s*s + 2.0*s*a + 3.0*a*a), "radiative reservoir partial")?,
        ])
    }
}

/// Differentiate a nodal-temperature goal on an accepted ambient-radiation field.
///
/// `problem.boundary` is the ORIGINAL convective partition, not the combined
/// boundary from a radiation solve. `regions` selects uniform Robin controls
/// and MUST include every radiation patch and every mean-dependent convection
/// law. Its `convective_mean_slopes_w_m2_k2` entries are dh_c/d(mean wall T);
/// use zero for prescribed coefficients. The caller owns and binds these
/// convection laws to the supplied field (for example a natural-convection
/// fixed point). Other Robin coefficients are fixed.
///
/// Rebuild the exact physical secant at the ACTUAL wall mean, not the final
/// driving-temperature approximation. The production residual must pass the
/// caller's linear tolerance on this unchanged field. An outer solve that is
/// adequate for visualization can therefore refuse a tighter derivative:
/// tighten its temperature/watt tolerances rather than freeze the radiation
/// or run a hidden replacement primal. Smooth k(T) and matching contacts stay
/// in the same Jacobian. Complete low-rank feedback and its true transpose use
/// the existing bounded FGMRES path; no dense trace matrix is assembled.
///
/// The returned convection entries are independent physical controls. A
/// control also changing h_c must compose `convection.log_htc*d(log h_c)`;
/// its indirect effect through wall temperature is already in the dual.
/// This does not differentiate geometry, material laws, card selection, air
/// advection/recirculation, enclosure radiosity, or the location of a maximum.
/// `max_feedback_entries` caps two nodal factor vectors per selected region,
/// not total memory/RSS. At most 64 regions/patches are accepted.
#[allow(clippy::too_many_arguments)]
pub fn pullback_ambient_radiation(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    linear: LinearConfig,
    temperature: &[f64],
    regions: &[&str],
    convective_mean_slopes_w_m2_k2: &[f64],
    patches: &[AmbientRadiationPatch],
    nodal_weights: &[f64],
    max_feedback_entries: usize,
) -> Result<AmbientRadiationGradient, ConductionError> {
    poll(cx, 0)?;
    let n = problem.mesh.vertex_count();
    if regions.len() > 64 || patches.is_empty() || patches.len() > 64
        || n.checked_mul(2).and_then(|v| v.checked_mul(regions.len()))
            .is_none_or(|v| v > max_feedback_entries)
    {
        return Err(invalid("ambient adjoint exceeds the factor-entry/64-region limit or has no patches"));
    }
    if temperature.len() != n || nodal_weights.len() != n
        || convective_mean_slopes_w_m2_k2.len() != regions.len()
    {
        return Err(invalid("ambient adjoint field, goal or convection-slope length differs from its binding"));
    }
    for (i, &value) in temperature.iter().chain(nodal_weights)
        .chain(convective_mean_slopes_w_m2_k2).enumerate()
    {
        if i % 512 == 0 { poll(cx, i)?; }
        finite(value, "ambient adjoint input")?;
    }
    let mut seen = BTreeSet::new();
    for &name in regions {
        if !seen.insert(name) { return Err(invalid("duplicate ambient-adjoint Robin region")); }
    }
    seen.clear();
    for patch in patches {
        if !seen.insert(patch.region()) || !regions.contains(&patch.region()) {
            return Err(invalid("every radiation patch must have exactly one selected Robin region"));
        }
    }

    let mut replacements = Vec::with_capacity(regions.len());
    let mut h_slopes = Vec::with_capacity(regions.len());
    let mut reference_slopes = Vec::with_capacity(regions.len());
    // (h_c, T_c, h_total, T_combined). Original controls must not be lost when
    // the producer folds two independent reservoirs into one numerical row.
    let mut physical = Vec::with_capacity(regions.len());
    let mut radiative = vec![None; patches.len()];
    for (i, &name) in regions.iter().enumerate() {
        poll(cx, i)?;
        let region = problem.boundary.region_names().iter().position(|r| r == name)
            .ok_or_else(|| invalid("unknown ambient-adjoint Robin region"))?;
        let ThermalBc::Robin { htc: ScalarField::Uniform(hc), t_ref: ScalarField::Uniform(tc) }
            = &problem.boundary.conditions()[region]
        else { return Err(invalid("ambient adjoints require selected uniform Robin rows")); };
        let hc_slope = convective_mean_slopes_w_m2_k2[i];
        if let Some((j, patch)) = patches.iter().enumerate().find(|(_, p)| p.region() == name) {
            let mut area = 0.0;
            let mut integral = 0.0;
            for (slot, face) in problem.mesh.boundary().iter().enumerate() {
                if slot % 512 == 0 { poll(cx, slot)?; }
                if problem.boundary.region_for(slot) != Some(region) { continue; }
                area = finite(area + face.area, "ambient-adjoint patch area")?;
                for &v in &face.vertices {
                    integral = finite(integral + (face.area/3.0)*temperature[v as usize],
                        "ambient-adjoint wall integral")?;
                }
            }
            if area <= 0.0 { return Err(invalid("radiating region has no positive trace area")); }
            let mean = finite(integral/area, "ambient-adjoint wall mean")?;
            let hr = patch.secant_coefficient_w_m2_k(mean)?;
            let [hr_slope, hr_ambient] = patch.secant_partials_w_m2_k2(mean)?;
            let ambient = patch.ambient_temperature_k();
            let h = finite(hc + hr, "combined adjoint coefficient")?;
            if h <= 0.0 { return Err(invalid("combined adjoint coefficient is not positive")); }
            let reference = finite((hc/h)*tc + (hr/h)*ambient, "combined adjoint reference")?;
            replacements.push((region, h, reference));
            h_slopes.push(finite(hc_slope + hr_slope, "combined coefficient slope")?);
            reference_slopes.push(finite((hc_slope/h)*(tc-reference)
                + (hr_slope/h)*(ambient-reference), "combined reference slope")?);
            physical.push((*hc, *tc, h, reference));
            radiative[j] = Some((i, hr, hr_ambient));
        } else {
            // Preserve a fixed or mean-dependent convection-only row exactly.
            replacements.push((region, *hc, *tc));
            h_slopes.push(hc_slope);
            reference_slopes.push(0.0);
            physical.push((*hc, *tc, *hc, *tc));
        }
    }
    let combined = problem.boundary.with_uniform_robin_replacements(&replacements)?;
    let combined_problem = ConductionProblem { boundary: &combined, ..problem };
    let mut gradient = RobinResponse::pullback_mean_robin(cx, combined_problem, interfaces,
        linear, temperature, regions, &h_slopes, &reference_slopes, nodal_weights,
        max_feedback_entries)?;
    let mut reservoirs = Vec::with_capacity(patches.len());
    let mut emissivities = Vec::with_capacity(patches.len());
    for (j, patch) in patches.iter().enumerate() {
        poll(cx, j)?;
        let (i, hr, hr_ambient) = radiative[j]
            .ok_or_else(|| invalid("radiation derivative is missing its bound trace"))?;
        let (hc, tc, h, reference) = physical[i];
        let g_ref = gradient.references[i];
        let g_log = gradient.log_htc[i];
        let reservoir = patch.ambient_temperature_k();
        let coefficient_bar = finite(g_log + g_ref*(reservoir-reference),
            "radiative coefficient contraction")?;
        reservoirs.push(finite(g_ref*(hr/h) + coefficient_bar*(hr_ambient/h),
            "radiative reservoir derivative")?);
        emissivities.push(finite(coefficient_bar*(hr/h)/patch.emissivity().value(),
            "emissivity derivative")?);
        gradient.references[i] = finite(g_ref*(hc/h), "original convection reference derivative")?;
        gradient.log_htc[i] = finite((hc/h)*(g_log + g_ref*(tc-reference)),
            "original convection coefficient derivative")?;
    }
    poll(cx, gradient.iterations)?;
    Ok(AmbientRadiationGradient {
        regions: regions.iter().map(|name| (*name).to_owned()).collect(),
        convection: gradient,
        radiation_regions: patches.iter().map(|patch| patch.region().to_owned()).collect(),
        reservoir_temperatures: reservoirs,
        emissivities,
    })
}
