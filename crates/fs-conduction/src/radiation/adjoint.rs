//! Complete steady FEM transpose for reflected gray-diffuse enclosure exchange.
//!
//! The original convection and uniform enclosure flux share one trace:
//! h(T - T_c) + q = h(T - (T_c - q/h)). This algebraic reference is rebuilt
//! at the unchanged accepted field. Its complete cross-patch derivative is
//! composed with the production material/contact Jacobian, not frozen.

use std::collections::BTreeSet;

use fs_exec::Cx;

use super::{
    EnclosureRadiosityLinearization, GrayDiffuseEnclosure, radiation_error,
    radiation_mesh_identity_with_poll, require_temperature,
};
use crate::{
    ConductionError, ConductionProblem, LinearConfig, ScalarField, ThermalBc, ThermalInterfaces,
    adjoint::{RobinGradient, RobinResponse},
};

/// Discrete local derivatives through the complete reflected-radiation law.
#[derive(Debug, Clone, PartialEq)]
pub struct EnclosureRadiationGradient {
    /// Selected original convective regions, in caller order.
    pub regions: Vec<String>,
    /// Original convective reference, ln(h), and assembled nodal-load partials.
    /// These use the complete thermal dual, including reflected exchange.
    pub convection: RobinGradient,
    /// Exact enclosure region identities, in view-factor row order.
    pub radiation_regions: Vec<String>,
    /// Derivative per unit absolute emissivity, in view-factor row order.
    /// Geometry, card selection and query temperature remain fixed. At eps=1,
    /// only feasible parameter directions have physical meaning.
    pub emissivities: Vec<f64>,
}

fn invalid(what: impl Into<String>) -> ConductionError {
    radiation_error(
        "steady-enclosure-adjoint",
        what,
        "retain the exact enclosure, original uniform positive convection, unchanged field, and explicit derivative budgets",
    )
}

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "steady-enclosure-adjoint",
        at,
    })
}

fn checked(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("nonfinite enclosure derivative arithmetic"))
    }
}

/// Differentiate a retained steady field with state-independent convection.
///
/// Every enclosure surface must exactly cover its same-named original
/// uniform positive Robin region, and be included in regions. Other selected
/// regions provide independent convection controls. Unselected convection is
/// fixed. The full smooth heterogeneous material and matching-contact tangent
/// is retained; emissivity and view factors do not depend on temperature.
///
/// Fresh radiosity must pass heat_tolerance_w and the unchanged field must
/// pass the original linear residual tolerance. No replacement primal solve
/// or derivative of the outer iteration is performed. Four nodal vectors per
/// selected region must fit max_feedback_entries; at most 64 are admitted.
/// The result is a discrete local derivative, not an inverse certificate,
/// continuum error bound, geometric derivative or uncertainty bound.
#[allow(clippy::too_many_arguments)]
pub fn pullback_gray_diffuse_enclosure(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    linear: LinearConfig,
    temperature: &[f64],
    regions: &[&str],
    enclosure: &GrayDiffuseEnclosure,
    heat_tolerance_w: f64,
    nodal_weights: &[f64],
    max_feedback_entries: usize,
) -> Result<EnclosureRadiationGradient, ConductionError> {
    if regions.len() > 64 {
        return Err(invalid("enclosure derivative exceeds 64 selected regions"));
    }
    let m = regions.len();
    EnclosureRadiationGradient::pullback_with_reference_feedback(
        cx, problem, interfaces, linear, temperature, regions, &vec![0.0; m],
        &vec![0.0; m * m], enclosure, heat_tolerance_w, nodal_weights, max_feedback_entries,
    )
}

impl EnclosureRadiationGradient {
    /// Include caller-bound smooth convection and cross-reference feedback.
    ///
    /// convective_mean_slopes_w_m2_k2 contains dh_i/d(mean_i T).
    /// reference_feedback is the row-major partial d(T_c_i)/d(mean_j T),
    /// including its diagonal, in regions order. The caller must evaluate
    /// every original h_i and T_c_i at the retained field before calling.
    /// Only original convection belongs in an air law; enclosure radiation
    /// stays internal. This entry point does not itself bind air transport,
    /// natural convection, hydraulic controls or material-card selection.
    ///
    /// For a radiating region, the shifted reference derivative is
    /// dT_c_i/dmean_j - (dq_i/dmean_j)/h_i, with the additional local
    /// q_i*(dh_i/dmean_i)/h_i^2 term. The returned ln(h) partial restores
    /// the original physical reference. Supplying a shifted reference as
    /// the original boundary would count radiation twice and is invalid.
    #[allow(clippy::too_many_arguments)]
    pub fn pullback_with_reference_feedback(
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        linear: LinearConfig,
        temperature: &[f64],
        regions: &[&str],
        convective_mean_slopes_w_m2_k2: &[f64],
        reference_feedback: &[f64],
        enclosure: &GrayDiffuseEnclosure,
        heat_tolerance_w: f64,
        nodal_weights: &[f64],
        max_feedback_entries: usize,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        let n = problem.mesh.vertex_count();
        let m = regions.len();
        let p = enclosure.surfaces().len();
        if !(2..=64).contains(&p) || m > 64 || p > m
            || n.checked_mul(4).and_then(|v| v.checked_mul(m))
                .is_none_or(|v| v > max_feedback_entries)
        {
            return Err(invalid("enclosure derivative exceeds its factor-entry or 64-region budget"));
        }
        if temperature.len() != n || nodal_weights.len() != n
            || convective_mean_slopes_w_m2_k2.len() != m
            || reference_feedback.len() != m * m
            || !(linear.tolerance.is_finite() && linear.tolerance > 0.0 && linear.tolerance < 1.0)
            || linear.max_iterations == 0 || linear.restart == 0
        {
            return Err(invalid("enclosure derivative needs complete fields, slopes and explicit linear limits"));
        }
        for (i, &v) in temperature.iter().chain(nodal_weights)
            .chain(convective_mean_slopes_w_m2_k2).chain(reference_feedback).enumerate()
        {
            if i % 256 == 0 { poll(cx, i)?; }
            checked(v)?;
        }
        let mut selected = Vec::with_capacity(m);
        let mut physical = Vec::with_capacity(m);
        let mut seen = BTreeSet::new();
        for &name in regions {
            poll(cx, selected.len())?;
            if !seen.insert(name) { return Err(invalid("duplicate selected convection region")); }
            let region = problem.boundary.region_names().iter().position(|r| r == name)
                .ok_or_else(|| invalid("selected convection region is absent"))?;
            let ThermalBc::Robin {
                htc: ScalarField::Uniform(h),
                t_ref: ScalarField::Uniform(reference),
            } = problem.boundary.conditions()[region] else {
                return Err(invalid("enclosure derivatives require original uniform Robin controls"));
            };
            if !h.is_finite() || h <= 0.0 {
                return Err(invalid("enclosure derivatives require finite positive convection"));
            }
            require_temperature(name, "original convective reference", reference)?;
            selected.push(region);
            physical.push((h, reference));
        }

        // One cancellable mesh identity and one coverage scan, rather than
        // hashing the complete volume independently for every surface.
        let mesh_identity = radiation_mesh_identity_with_poll(problem.mesh, |at| poll(cx, at))?;
        let mut counts = vec![0usize; problem.boundary.region_names().len()];
        for slot in 0..problem.mesh.boundary().len() {
            if slot % 256 == 0 { poll(cx, slot)?; }
            if let Some(region) = problem.boundary.region_for(slot) { counts[region] += 1; }
        }
        let mut ports = Vec::with_capacity(p);
        let mut means = Vec::with_capacity(p);
        for surface in enclosure.surfaces() {
            poll(cx, ports.len())?;
            if surface.mesh_identity != mesh_identity {
                return Err(invalid("enclosure surface belongs to a different exact mesh"));
            }
            let port = regions.iter().position(|&name| name == surface.name())
                .ok_or_else(|| invalid("every enclosure patch needs one selected same-named convection region"))?;
            let region = selected[port];
            if counts[region] != surface.face_slots().len() {
                return Err(invalid("enclosure surface must exactly cover its convection region"));
            }
            let mut integral = 0.0;
            for (i, &slot) in surface.face_slots().iter().enumerate() {
                if i % 256 == 0 { poll(cx, i)?; }
                if problem.boundary.region_for(slot) != Some(region) {
                    return Err(invalid("enclosure surface must exactly cover its convection region"));
                }
                let face = &problem.mesh.boundary()[slot];
                let mut sum = 0.0;
                for &v in &face.vertices { sum = checked(sum + temperature[v as usize])?; }
                let mean = checked(sum / 3.0)?;
                require_temperature(surface.name(), "surface nodal average", mean)?;
                integral = checked(face.area.mul_add(mean, integral))?;
            }
            let mean = checked(integral / surface.area_m2())?;
            if let Some((lower, upper)) = surface.emissivity().temperature_validity_k()
                && (mean <= lower || mean >= upper)
            {
                return Err(invalid("enclosure derivative requires interior emissivity temperature validity"));
            }
            ports.push(port);
            means.push(mean);
        }

        let fresh = enclosure.solve(cx, &means)?;
        let radiosity = EnclosureRadiosityLinearization::new(
            cx, enclosure, &means, heat_tolerance_w,
        )?;
        let mut feedback = reference_feedback.to_vec();
        let mut local_reference_slopes = vec![0.0; m];
        let mut replacements = Vec::with_capacity(p);
        for (i, &port) in ports.iter().enumerate() {
            poll(cx, i)?;
            let (h, reference) = physical[port];
            let flux = checked(fresh.net_outward_flux_w_m2[i])?;
            replacements.push((selected[port], h, checked(reference - flux / h)?));
            local_reference_slopes[port] =
                checked((flux / h) * (convective_mean_slopes_w_m2_k2[port] / h))?;

            // A flux-basis pullback supplies row i of the full nonlocal
            // patch Jacobian through the shared radiosity transpose owner.
            let mut weights = vec![0.0; p];
            weights[i] = 1.0;
            let partials = radiosity.pullback(cx, &weights, linear.tolerance, linear.max_iterations)?;
            for (j, &other) in ports.iter().enumerate() {
                feedback[port * m + other] =
                    checked(feedback[port * m + other] - partials.temperatures[j] / h)?;
            }
        }
        let boundary = problem.boundary.with_uniform_robin_replacements(&replacements)?;
        let combined_problem = ConductionProblem { boundary: &boundary, ..problem };
        let mut convection = RobinResponse::pullback_mean_robin_with_reference_feedback(
            cx, combined_problem, interfaces, linear, temperature, regions,
            convective_mean_slopes_w_m2_k2, &local_reference_slopes, &feedback,
            nodal_weights, max_feedback_entries,
        )?;

        // Restore physical convection controls before returning. The mean
        // feedback is already in the dual and must not be applied a second time.
        let mut flux_weights = vec![0.0; p];
        for (i, &port) in ports.iter().enumerate() {
            poll(cx, i)?;
            let h = physical[port].0;
            let flux = fresh.net_outward_flux_w_m2[i];
            flux_weights[i] = checked(-convection.references[port] / h)?;
            convection.log_htc[port] =
                checked(convection.log_htc[port] + convection.references[port] * (flux / h))?;
        }
        let radiation = radiosity.pullback(
            cx, &flux_weights, linear.tolerance, linear.max_iterations,
        )?;
        let emissivities = radiation.log_emissivities.iter().zip(enclosure.surfaces())
            .map(|(&value, surface)| checked(value / surface.emissivity().value()))
            .collect::<Result<Vec<_>, _>>()?;
        poll(cx, convection.iterations)?;
        Ok(Self {
            regions: regions.iter().map(|name| (*name).to_owned()).collect(),
            convection,
            radiation_regions: enclosure.surfaces().iter().map(|s| s.name().to_owned()).collect(),
            emissivities,
        })
    }
}
