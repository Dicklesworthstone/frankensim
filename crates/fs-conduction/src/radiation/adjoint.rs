//! Complete steady FEM transpose for reflected gray-diffuse enclosure exchange.
//!
//! Original Robin and Neumann patches retain their physical laws. Convection
//! uses h(T - T_c) + q = h(T - (T_c - q/h)); prescribed-flux patches add q
//! directly. The complete cross-patch derivative is composed with the same
//! material/contact Jacobian at the unchanged accepted field.

use std::collections::BTreeSet;

use fs_exec::Cx;

use super::{
    EnclosureRadiosityLinearization, GrayDiffuseEnclosure, radiation_error,
    radiation_mesh_identity_with_poll, require_temperature,
};
use crate::{
    ConductionError, ConductionProblem, LinearConfig, ScalarField, ThermalBc, ThermalBoundary,
    ThermalInterfaces,
    adjoint::{MeanFluxFeedback, RobinGradient, RobinResponse},
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
    /// Total shared radiosity-transpose iterations, including patch Jacobian
    /// construction and the final emissivity contraction.
    pub radiosity_iterations: usize,
    /// Largest true relative residual among those radiosity transposes.
    pub radiosity_relative_residual: f64,
}

fn invalid(what: impl Into<String>) -> ConductionError {
    radiation_error(
        "steady-enclosure-adjoint",
        what,
        "retain the exact enclosure, original uniform convection or prescribed-flux laws, unchanged field, and explicit derivative budgets",
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
/// uniform Robin or Neumann region. Radiating Robin regions must be included
/// in regions; Neumann patches are independent of those controls. Other selected
/// regions provide independent convection controls. Unselected convection is
/// fixed. The full smooth heterogeneous material and matching-contact tangent
/// is retained; emissivity and view factors do not depend on temperature.
///
/// Fresh radiosity must pass heat_tolerance_w and the unchanged field must
/// pass the original linear residual tolerance. No replacement primal solve
/// or derivative of the outer iteration is performed. Four nodal vectors per
/// selected Robin region must fit max_feedback_entries; a mixed-flux enclosure
/// additionally needs two nodal factors per radiation patch. Each order admits
/// at most 64 entries, including zero selected Robin controls.
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

/// Complete unchanged-field boundary representation of enclosure exchange.
///
/// This binding exposes the same physical Jacobian used by the full enclosure
/// dual, so prescribed-temperature lifts and other operator consumers cannot
/// silently freeze reflected radiation. Construction verifies the exact mesh,
/// surface coverage, fresh radiosity and coefficient budgets. Consumers must
/// still enforce their unchanged-field residual and complete transpose gates.
pub struct EnclosureRadiationLinearization {
    vertex_count: usize,
    boundary: ThermalBoundary,
    regions: Vec<String>,
    convective_mean_slopes: Vec<f64>,
    local_reference_slopes: Vec<f64>,
    feedback: Vec<f64>,
    physical: Vec<(f64, f64)>,
    ports: Vec<Option<usize>>,
    flux_feedback: Option<MeanFluxFeedback>,
    outward_flux: Vec<f64>,
    radiosity: EnclosureRadiosityLinearization,
    radiosity_iterations: usize,
    radiosity_relative_residual: f64,
    radiation_regions: Vec<String>,
    emissivities: Vec<f64>,
}

impl EnclosureRadiationLinearization {
    /// Fresh combined boundary: shifted Robin references or additive Neumann flux.
    pub fn boundary(&self) -> &ThermalBoundary { &self.boundary }

    /// Selected original Robin regions, in caller order.
    pub fn regions(&self) -> &[String] { &self.regions }

    /// Local derivative dh_i/d(mean_i T) of original convection.
    pub fn convective_mean_slopes(&self) -> &[f64] { &self.convective_mean_slopes }

    /// Local shifted-reference correction q_i (dh_i/dmean_i)/h_i^2.
    /// Add this diagonal to reference_feedback, exactly as RobinResponse does.
    pub fn local_reference_slopes(&self) -> &[f64] { &self.local_reference_slopes }

    /// Row-major dT_c_i/dmean_j - (dq_i/dmean_j)/h_i between Robin controls.
    /// Incoming air-reference diagonals remain included. Reflected terms
    /// involving a Neumann patch are supplied by flux_feedback; the local
    /// reference slopes supply the h-dependent correction.
    pub fn reference_feedback(&self) -> &[f64] { &self.feedback }

    /// Additional reflected-flux tangent involving prescribed-flux patches.
    /// These factors are independent of actual convective control ports.
    pub fn flux_feedback(&self) -> Option<&MeanFluxFeedback> { self.flux_feedback.as_ref() }

    /// Total Krylov iterations for all patch-flux transpose rows in this binding.
    pub fn radiosity_iterations(&self) -> usize { self.radiosity_iterations }

    /// Largest true relative residual of those patch-flux transpose rows.
    pub fn radiosity_relative_residual(&self) -> f64 { self.radiosity_relative_residual }

    /// Bind physical convection and reflected exchange at one unchanged field.
    ///
    /// The supplied original boundary must already contain h and T_c evaluated
    /// at temperature. convective_mean_slopes contains dh_i/dmean_i;
    /// reference_feedback contains all dT_c_i/dmean_j, including its diagonal.
    /// A radiation-shifted reference is not an original convective reference.
    #[allow(clippy::too_many_arguments)]
    pub fn bind(
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        linear: LinearConfig,
        temperature: &[f64],
        regions: &[&str],
        convective_mean_slopes_w_m2_k2: &[f64],
        reference_feedback: &[f64],
        enclosure: &GrayDiffuseEnclosure,
        heat_tolerance_w: f64,
        max_feedback_entries: usize,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        let n = problem.mesh.vertex_count();
        let m = regions.len();
        let p = enclosure.surfaces().len();
        if !(2..=64).contains(&p) || m > 64
            || n.checked_mul(4).and_then(|v| v.checked_mul(m))
                .is_none_or(|v| v > max_feedback_entries)
        {
            return Err(invalid("enclosure derivative exceeds its factor-entry or 64-region budget"));
        }
        if temperature.len() != n
            || convective_mean_slopes_w_m2_k2.len() != m
            || reference_feedback.len() != m * m
            || !(linear.tolerance.is_finite() && linear.tolerance > 0.0 && linear.tolerance < 1.0)
            || linear.max_iterations == 0 || linear.restart == 0
        {
            return Err(invalid("enclosure derivative needs complete fields, slopes and explicit linear limits"));
        }
        for (i, &v) in temperature.iter()
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
            let region = problem.boundary.region_names().iter()
                .position(|name| name == surface.name())
                .ok_or_else(|| invalid("enclosure patch has no same-named original boundary"))?;
            let port = match problem.boundary.conditions()[region] {
                ThermalBc::Robin {
                    htc: ScalarField::Uniform(_), t_ref: ScalarField::Uniform(_),
                } => Some(regions.iter().position(|&name| name == surface.name())
                    .ok_or_else(|| invalid("each radiating Robin patch needs its original selected control"))?),
                ThermalBc::Neumann { outward_flux: ScalarField::Uniform(flux) } => {
                    checked(flux)?;
                    None
                }
                _ => return Err(invalid("enclosure derivatives require uniform Robin or prescribed-flux patches")),
            };
            if counts[region] != surface.face_slots().len() {
                return Err(invalid("enclosure surface must exactly cover its original boundary region"));
            }
            let mut integral = 0.0;
            for (i, &slot) in surface.face_slots().iter().enumerate() {
                if i % 256 == 0 { poll(cx, i)?; }
                if problem.boundary.region_for(slot) != Some(region) {
                    return Err(invalid("enclosure surface must exactly cover its original boundary region"));
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
        let mut neumann_additions = Vec::with_capacity(p);
        let mut flux_derivatives = vec![0.0; p * p];
        let mut radiosity_iterations = 0usize;
        let mut radiosity_relative_residual = 0.0_f64;
        for (i, &port) in ports.iter().enumerate() {
            poll(cx, i)?;
            let flux = checked(fresh.net_outward_flux_w_m2[i])?;
            if let Some(port) = port {
                let (h, reference) = physical[port];
                replacements.push((selected[port], h, checked(reference - flux / h)?));
                local_reference_slopes[port] =
                    checked((flux / h) * (convective_mean_slopes_w_m2_k2[port] / h))?;
            } else {
                let region = problem.boundary.region_names().iter()
                    .position(|name| name == enclosure.surfaces()[i].name())
                    .ok_or_else(|| invalid("enclosure prescribed-flux region disappeared"))?;
                neumann_additions.push((region, flux));
            }

            // A flux-basis pullback supplies row i of the full nonlocal
            // patch Jacobian through the shared radiosity transpose owner.
            let mut weights = vec![0.0; p];
            weights[i] = 1.0;
            let partials = radiosity.pullback(cx, &weights, linear.tolerance, linear.max_iterations)?;
            radiosity_iterations = radiosity_iterations.checked_add(partials.iterations)
                .ok_or_else(|| invalid("enclosure radiosity transpose iteration count overflow"))?;
            radiosity_relative_residual =
                radiosity_relative_residual.max(checked(partials.relative_residual)?);
            for (j, &other) in ports.iter().enumerate() {
                if let (Some(port), Some(other)) = (port, other) {
                    feedback[port * m + other] = checked(
                        feedback[port * m + other] - partials.temperatures[j] / physical[port].0,
                    )?;
                } else {
                    // No physical convection coefficient exists on a flux
                    // patch. Compose its remaining G*dq/dmean*W directly.
                    flux_derivatives[i * p + j] = checked(partials.temperatures[j])?;
                }
            }
        }
        let boundary = problem.boundary.with_uniform_robin_replacements(&replacements)?
            .with_uniform_neumann_additions(&neumann_additions)?;
        let flux_feedback = if ports.iter().any(Option::is_none) {
            let robin_entries = n.checked_mul(m).and_then(|v| v.checked_mul(4))
                .ok_or_else(|| invalid("enclosure Robin factor allowance overflow"))?;
            let remaining = max_feedback_entries.checked_sub(robin_entries)
                .ok_or_else(|| invalid("enclosure Robin factors exhaust the feedback allowance"))?;
            let names: Vec<_> = enclosure.surfaces().iter().map(|surface| surface.name()).collect();
            Some(MeanFluxFeedback::new(
                cx, problem.mesh, problem.boundary, &names, &flux_derivatives, remaining,
            )?)
        } else { None };
        poll(cx, 0)?;
        Ok(Self {
            vertex_count: n,
            boundary,
            regions: regions.iter().map(|name| (*name).to_owned()).collect(),
            convective_mean_slopes: convective_mean_slopes_w_m2_k2.to_vec(),
            local_reference_slopes,
            feedback,
            physical,
            ports,
            flux_feedback,
            outward_flux: fresh.net_outward_flux_w_m2,
            radiosity,
            radiosity_iterations,
            radiosity_relative_residual,
            radiation_regions: enclosure.surfaces().iter().map(|s| s.name().to_owned()).collect(),
            emissivities: enclosure.surfaces().iter().map(|s| s.emissivity().value()).collect(),
        })
    }
    /// Restore original physical convection controls and contract emissivity.
    ///
    /// convection must come from the complete thermal Jacobian built with this
    /// binding's combined boundary and mean/reference feedback. A transient
    /// caller may include its storage tangent and backward dual right-hand side;
    /// no steady PDE solve or steady physical residual test occurs here. The
    /// caller owns primal validity and the binding of that dual to its state.
    /// Returned reference/log(h) bars are per unit assembled nodal load; a
    /// transient caller supplies its time-step weighting exactly once.
    pub fn remap_gradient(
        &self, cx: &Cx<'_>, mut convection: RobinGradient, linear: LinearConfig,
    ) -> Result<EnclosureRadiationGradient, ConductionError> {
        poll(cx, 0)?;
        if convection.references.len() != self.regions.len()
            || convection.log_htc.len() != self.regions.len()
            || convection.nodal_load.len() != self.vertex_count
            || !(linear.tolerance.is_finite() && linear.tolerance > 0.0 && linear.tolerance < 1.0)
            || linear.max_iterations == 0 || linear.restart == 0
            || !convection.relative_residual.is_finite()
            || convection.relative_residual < 0.0
            || convection.relative_residual >= linear.tolerance
            || self.radiosity_relative_residual >= linear.tolerance
        {
            return Err(invalid("enclosure control remapping requires the complete converged thermal dual and original trace order"));
        }
        for (i, &value) in convection.references.iter().chain(&convection.log_htc)
            .chain(&convection.nodal_load).enumerate()
        {
            if i % 256 == 0 { poll(cx, i)?; }
            checked(value)?;
        }
        // Restore physical convection controls before returning. The mean
        // feedback is already in the dual and must not be applied a second time.
        let direct_flux = self.flux_feedback.as_ref()
            .map(|feedback| feedback.flux_pullback(cx, &convection.nodal_load)).transpose()?;
        let mut flux_weights = vec![0.0; self.ports.len()];
        for (i, &port) in self.ports.iter().enumerate() {
            poll(cx, i)?;
            if let Some(port) = port {
                let h = self.physical[port].0;
                let flux = self.outward_flux[i];
                flux_weights[i] = checked(-convection.references[port] / h)?;
                convection.log_htc[port] =
                    checked(convection.log_htc[port] + convection.references[port] * (flux / h))?;
            } else {
                flux_weights[i] = checked(-direct_flux.as_ref()
                    .ok_or_else(|| invalid("prescribed-flux control lost its trace factors"))?[i])?;
            }
        }
        let radiation = self.radiosity.pullback(
            cx, &flux_weights, linear.tolerance, linear.max_iterations,
        )?;
        let emissivities = radiation.log_emissivities.iter().zip(&self.emissivities)
            .map(|(&value, &emissivity)| checked(value / emissivity))
            .collect::<Result<Vec<_>, _>>()?;
        poll(cx, convection.iterations)?;
        Ok(EnclosureRadiationGradient {
            regions: self.regions.clone(),
            convection,
            radiation_regions: self.radiation_regions.clone(),
            emissivities,
            radiosity_iterations: self.radiosity_iterations.checked_add(radiation.iterations)
                .ok_or_else(|| invalid("enclosure radiosity transpose iteration count overflow"))?,
            radiosity_relative_residual: self.radiosity_relative_residual
                .max(checked(radiation.relative_residual)?),
        })
    }

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
        let bound = EnclosureRadiationLinearization::bind(
            cx, problem, linear, temperature, regions, convective_mean_slopes_w_m2_k2,
            reference_feedback, enclosure, heat_tolerance_w, max_feedback_entries,
        )?;
        let combined_problem = ConductionProblem { boundary: bound.boundary(), ..problem };
        let convection = RobinResponse::pullback_mean_robin_with_reference_and_flux_feedback(
            cx, combined_problem, interfaces, linear, temperature, regions,
            bound.convective_mean_slopes(), bound.local_reference_slopes(),
            bound.reference_feedback(), bound.flux_feedback(), nodal_weights, max_feedback_entries,
        )?;

        bound.remap_gradient(cx, convection, linear)
    }
}
