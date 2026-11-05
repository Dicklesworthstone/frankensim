//! Total derivatives of the SAME mean-temperature secant Robin model as the
//! ambient-radiation producer. In particular this is not uniform patch flux
//! and is not the integral of pointwise T(x)^4 over a non-isothermal wall.

use super::{AmbientRadiationPatch, ConductionError, ConductionProblem, Cx,
    ScalarField, ThermalBc, ThermalInterfaces, finite, poll, radiation_error,
    STEFAN_BOLTZMANN_W_M2_K4};
use crate::{LinearConfig, ThermalBoundary};
use crate::adjoint::{RobinGradient, RobinResponse};
use std::collections::BTreeSet;

/// Discrete nodal-temperature sensitivities through ambient radiation.
/// Every ordering is explicit; no physical or uncertainty certification is
/// inferred from the residual-checked numerical derivative.
#[derive(Debug, Clone, PartialEq)]
pub struct AmbientRadiationGradient {
    /// Selected thermal regions, in the caller's order: original convection
    /// controls and explicitly zero-flux radiation-only traces.
    pub regions: Vec<String>,
    /// Derivatives with respect to ORIGINAL convective references, ln(h_c),
    /// and assembled nodal loads. These include radiative state feedback;
    /// they are not partials of the artificially combined Robin row. Reference
    /// and log-HTC entries for radiation-only traces are exactly zero; they do
    /// not create a convection law or define the logarithm of a zero HTC.
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

/// Exact unchanged-field ambient-radiation boundary and its mean feedback.
///
/// Separating boundary binding from the transposed solve lets a complete live
/// contact tangent share the same radiation law and original-control remapping.
/// This value does not solve or certify the thermal endpoint. Its consumer must
/// validate the physical residual, solve the complete transpose, then remap the
/// combined Robin controls through this same binding.
pub struct AmbientRadiationLinearization {
    boundary: ThermalBoundary,
    regions: Vec<String>,
    h_slopes: Vec<f64>,
    reference_slopes: Vec<f64>,
    reference_feedback: Option<Vec<f64>>,
    physical: Vec<(f64, f64, f64, f64)>,
    radiative: Vec<Option<(usize, f64, f64)>>,
    patches: Vec<AmbientRadiationPatch>,
    vertices: usize,
}

fn invalid(what: &str) -> ConductionError {
    radiation_error("ambient-adjoint", what,
        "use the retained field, original uniform Robin or zero-flux radiating rows, exact radiation patches, and explicit derivative budgets")
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
/// `problem.boundary` is the ORIGINAL physical partition, not the combined
/// boundary from a radiation solve. `regions` selects uniform Robin controls
/// and uniform zero-flux traces carrying an explicit radiation patch,
/// and MUST include every radiation patch and every mean-dependent convection
/// law. Its `convective_mean_slopes_w_m2_k2` entries are dh_c/d(mean wall T);
/// use zero for prescribed coefficients and radiation-only traces. The caller
/// owns and binds these convection laws to the supplied field (for example a natural-convection
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
    pullback(cx, problem, interfaces, linear, temperature, regions,
        convective_mean_slopes_w_m2_k2, None, patches, nodal_weights, max_feedback_entries)
}

impl AmbientRadiationGradient {
    /// Differentiate coupled convection references, optionally with radiation.
    ///
    /// `reference_feedback` is the row-major m-by-m partial of each ORIGINAL
    /// convective reference with respect to all selected wall means, in K/K.
    /// The caller evaluates these references at the retained field before
    /// calling (for example by marching the actual fixed-flow air paths).
    /// Radiation supplies its OWN secant and weighted-reference derivatives;
    /// only h_c/(h_c+h_rad) of the convective-reference feedback enters the
    /// combined Robin reference. Adding the unweighted air matrix is wrong.
    /// Empty patches admit nonlinear solid/air coupling without radiation.
    /// A radiation-only trace has no convective reference: its feedback ROW
    /// and its convective mean slope must be zero. Its wall mean may still
    /// drive another trace's declared convection reference through a column.
    ///
    /// Returned convection references are partials of independent additive
    /// offsets of the convective law through the full state Jacobian. Air
    /// inlet controls must chain these through the inlet-to-reference map.
    /// A physical air h or flow control must ALSO differentiate that map;
    /// `log_htc` alone does not supply such a derivative. No hydraulics,
    /// temperature-dependent air properties, or radiative heat into air is
    /// inferred. Smooth k(T), matching contact and both state feedbacks share
    /// one transposed solve. Four nodal factor vectors per region must fit
    /// `max_feedback_entries`; at most 64 regions and 4096 coefficients.
    #[allow(clippy::too_many_arguments)]
    pub fn pullback_with_reference_feedback(
        cx: &Cx<'_>, problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
        linear: LinearConfig, temperature: &[f64], regions: &[&str],
        convective_mean_slopes_w_m2_k2: &[f64], reference_feedback: &[f64],
        patches: &[AmbientRadiationPatch], nodal_weights: &[f64], max_feedback_entries: usize,
    ) -> Result<Self, ConductionError> {
        pullback(cx, problem, interfaces, linear, temperature, regions,
            convective_mean_slopes_w_m2_k2, Some(reference_feedback), patches,
            nodal_weights, max_feedback_entries)
    }
}

#[allow(clippy::too_many_arguments)]
fn pullback(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
    linear: LinearConfig, temperature: &[f64], regions: &[&str],
    convective_mean_slopes_w_m2_k2: &[f64], reference_feedback: Option<&[f64]>,
    patches: &[AmbientRadiationPatch], nodal_weights: &[f64], max_feedback_entries: usize,
) -> Result<AmbientRadiationGradient, ConductionError> {
    let bound = AmbientRadiationLinearization::bind(cx, problem, temperature, regions,
        convective_mean_slopes_w_m2_k2, reference_feedback, patches, max_feedback_entries)?;
    let combined_problem = ConductionProblem { boundary: bound.boundary(), ..problem };
    let gradient = if let Some(feedback) = bound.reference_feedback() {
        RobinResponse::pullback_mean_robin_with_reference_feedback(cx, combined_problem,
            interfaces, linear, temperature, regions, bound.coefficient_mean_slopes(),
            bound.local_reference_slopes(), feedback, nodal_weights, max_feedback_entries)?
    } else {
        RobinResponse::pullback_mean_robin(cx, combined_problem, interfaces,
            linear, temperature, regions, bound.coefficient_mean_slopes(),
            bound.local_reference_slopes(), nodal_weights, max_feedback_entries)?
    };
    bound.remap_gradient(cx, gradient, linear)
}

impl AmbientRadiationLinearization {
    /// Rebuild the physical area-mean secant boundary at `temperature`.
    ///
    /// The boundary is the original convection/zero-flux partition; the selected
    /// regions must contain every radiating and mean-dependent convection trace.
    /// Slopes are dh_c/dmean(T), and optional row-major reference feedback is the
    /// derivative of the ORIGINAL convective references. Radiation contributes
    /// its own coefficient/reference slopes and scales each feedback row by
    /// h_c/(h_c+h_rad). Zero-flux radiating traces have no convective controls.
    ///
    /// The original two- or four-factor entry allowance and 64-trace cap apply.
    /// Smoothness/domain checks remain attached to the selected emissivity card.
    /// No primal field, contact law, material law or source selection is changed.
    #[allow(clippy::too_many_arguments)]
    pub fn bind(
        cx: &Cx<'_>, problem: ConductionProblem<'_>, temperature: &[f64],
        regions: &[&str], convective_mean_slopes_w_m2_k2: &[f64],
        reference_feedback: Option<&[f64]>, patches: &[AmbientRadiationPatch],
        max_feedback_entries: usize,
    ) -> Result<Self, ConductionError> {
    poll(cx, 0)?;
    let n = problem.mesh.vertex_count();
    let factors = if reference_feedback.is_some() { 4 } else { 2 };
    if regions.len() > 64 || (patches.is_empty() && reference_feedback.is_none()) || patches.len() > 64
        || n.checked_mul(factors).and_then(|v| v.checked_mul(regions.len()))
            .is_none_or(|v| v > max_feedback_entries)
    {
        return Err(invalid("ambient adjoint exceeds the factor-entry/64-region limit or has no patches"));
    }
    if temperature.len() != n
        || convective_mean_slopes_w_m2_k2.len() != regions.len()
    {
        return Err(invalid("ambient adjoint field or convection-slope length differs from its binding"));
    }
    for (i, &value) in temperature.iter()
        .chain(convective_mean_slopes_w_m2_k2).enumerate()
    {
        if i % 512 == 0 { poll(cx, i)?; }
        finite(value, "ambient adjoint input")?;
    }
    if let Some(feedback) = reference_feedback {
        if feedback.len() != regions.len()*regions.len() {
            return Err(invalid("coupled convection reference matrix has the wrong region count"));
        }
        for (i, &value) in feedback.iter().enumerate() {
            if i % 512 == 0 { poll(cx, i)?; }
            finite(value, "convective reference feedback")?;
        }
    }
    let mut seen = BTreeSet::new();
    for &name in regions {
        if !seen.insert(name) { return Err(invalid("duplicate ambient-adjoint thermal region")); }
    }
    seen.clear();
    for patch in patches {
        if !seen.insert(patch.region()) || !regions.contains(&patch.region()) {
            return Err(invalid("every radiation patch must have exactly one selected thermal region"));
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
            .ok_or_else(|| invalid("unknown ambient-adjoint thermal region"))?;
        let patch = patches.iter().enumerate().find(|(_, p)| p.region() == name);
        let hc_slope = convective_mean_slopes_w_m2_k2[i];
        let (hc, tc) = match &problem.boundary.conditions()[region] {
            ThermalBc::Robin { htc: ScalarField::Uniform(hc), t_ref: ScalarField::Uniform(tc) } => (*hc, *tc),
            ThermalBc::Neumann { outward_flux: ScalarField::Uniform(flux) } if *flux == 0.0 => {
                let (_, radiation) = patch.ok_or_else(||
                    invalid("a selected zero-flux trace requires its own radiation patch"))?;
                if hc_slope != 0.0 || reference_feedback.is_some_and(|matrix|
                    matrix[i*regions.len()..(i+1)*regions.len()].iter().any(|&v| v != 0.0)) {
                    return Err(invalid("a radiation-only trace has no convective slope or reference feedback"));
                }
                // The reservoir is only an algebraically inert placeholder
                // for T_c: every convection contribution is multiplied by 0.
                (0.0, radiation.ambient_temperature_k())
            }
            _ => return Err(invalid("ambient adjoints require uniform Robin or explicitly zero-flux radiating rows")),
        };
        if let Some((j, patch)) = patch {
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
            physical.push((hc, tc, h, reference));
            radiative[j] = Some((i, hr, hr_ambient));
        } else {
            // Preserve a fixed or mean-dependent convection-only row exactly.
            replacements.push((region, hc, tc));
            h_slopes.push(hc_slope);
            reference_slopes.push(0.0);
            physical.push((hc, tc, hc, tc));
        }
    }
    let boundary = problem.boundary.with_radiative_robin_replacements(&replacements)?;
    let reference_feedback = if let Some(feedback) = reference_feedback {
        let m = regions.len();
        let mut combined_feedback = Vec::with_capacity(feedback.len());
        for (i, &(hc, _, h, _)) in physical.iter().enumerate() {
            poll(cx, i)?;
            for &entry in &feedback[i*m..(i+1)*m] {
                combined_feedback.push(finite((hc/h)*entry, "combined convective feedback")?);
            }
        }
        Some(combined_feedback)
    } else {
        None
    };
    poll(cx, n)?;
    Ok(Self { boundary, regions: regions.iter().map(|name| (*name).to_owned()).collect(),
        h_slopes, reference_slopes, reference_feedback, physical, radiative,
        patches: patches.to_vec(), vertices: n })
    }

    /// Fresh physical secant boundary for the complete thermal/contact response.
    #[must_use]
    pub const fn boundary(&self) -> &ThermalBoundary { &self.boundary }
    /// Selected thermal trace order, including any radiation-only traces.
    #[must_use]
    pub fn regions(&self) -> &[String] { &self.regions }
    /// Derivatives of the combined coefficients with respect to local wall mean.
    #[must_use]
    pub fn coefficient_mean_slopes(&self) -> &[f64] { &self.h_slopes }
    /// Local weighted-reference slopes, excluding additional reference feedback.
    #[must_use]
    pub fn local_reference_slopes(&self) -> &[f64] { &self.reference_slopes }
    /// Additional reference feedback with each original convection row weighted.
    #[must_use]
    pub fn reference_feedback(&self) -> Option<&[f64]> { self.reference_feedback.as_deref() }

    /// Recover original convection, emissivity and bath controls from ONE dual.
    ///
    /// The supplied gradient must differentiate this binding's combined Robin
    /// rows through the complete state operator, in `regions()` order. Length,
    /// finiteness and the requested dual residual are checked; this is not an
    /// independent proof that an arbitrary caller's dual belongs to that solve.
    pub fn remap_gradient(
        &self, cx: &Cx<'_>, mut gradient: RobinGradient, linear: LinearConfig,
    ) -> Result<AmbientRadiationGradient, ConductionError> {
    poll(cx, 0)?;
    if gradient.nodal_load.len() != self.vertices
        || gradient.references.len() != self.regions.len()
        || gradient.log_htc.len() != self.regions.len()
        || !linear.tolerance.is_finite() || linear.tolerance <= 0.0
        || !gradient.relative_residual.is_finite()
        || gradient.relative_residual < 0.0
        || gradient.relative_residual > linear.tolerance
    {
        return Err(invalid("ambient gradient must retain the bound nodal/trace layout and meet its dual residual target"));
    }
    for (i, &value) in gradient.nodal_load.iter().chain(&gradient.references)
        .chain(&gradient.log_htc).enumerate()
    {
        if i % 512 == 0 { poll(cx, i)?; }
        finite(value, "ambient gradient input")?;
    }
    let mut reservoirs = Vec::with_capacity(self.patches.len());
    let mut emissivities = Vec::with_capacity(self.patches.len());
    for (j, patch) in self.patches.iter().enumerate() {
        poll(cx, j)?;
        let (i, hr, hr_ambient) = self.radiative[j]
            .ok_or_else(|| invalid("radiation derivative is missing its bound trace"))?;
        let (hc, tc, h, reference) = self.physical[i];
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
        regions: self.regions.clone(),
        convection: gradient,
        radiation_regions: self.patches.iter().map(|patch| patch.region().to_owned()).collect(),
        reservoir_temperatures: reservoirs,
        emissivities,
    })
    }
}
