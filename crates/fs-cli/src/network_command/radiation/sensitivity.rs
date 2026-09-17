//! Total steady derivative of the implemented mean-patch radiation closure.
//! Reconstruct the exact accepted combined rows, then use the same feedback
//! transpose as storage endpoints. Finite-tolerance roots are not certificates.
use super::*;
use fs_airflow::graph::thermal::coupled_transport::sensitivity::CoupledObjective;

// Keep the response-level implementation available to transient reconstruction
// without introducing a second radiation derivative or public numerical model.
#[path = "feedback.rs"]
pub(super) mod feedback;
#[path = "trajectory.rs"]
mod trajectory;

/// Exact inputs to the final steady solid evaluation, retained on gradient runs.
pub(super) struct Binding {
    pub(super) htc: BTreeMap<String, f64>,
    pub(super) references: Vec<f64>,
    pub(super) driving: Vec<f64>,
    pub(super) radiative_htc: Vec<f64>,
    pub(super) config: SolveConfig,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn pullback(
    policy: &Policy, request: &Request, cx: &Cx<'_>, network: &TransportNetwork<'_>,
    inner: &Inner, references: &[f64], htc: &BTreeMap<String, f64>,
    objective: &objective::ObjectiveState, derivations: &[convection::Derived],
) -> Result<(fan_gradient::CoolingGradient, String)> {
    poll(cx)?;
    let names = network.regions();
    let binding = inner.binding.as_ref().ok_or_else(|| bad("missing radiation derivative binding"))?;
    let boundary = request.boundary(&names, &binding.references, &binding.htc)?;
    let material = fs_conduction::ConductivityModel::isotropic_declared(request.conductivity).map_err(producer)?;
    let uniform = ScalarField::Uniform(request.source);
    let source = request.solid_data.nodal_source.as_ref().unwrap_or(&uniform);
    let problem = ConductionProblem { mesh: &request.mesh, boundary: &boundary,
        material: &material, element_materials: request.solid_data.element_materials.as_ref(), source };
    let solid = match &request.contacts {
        Some(contact) => RobinLinearization::new_with_interfaces(cx, problem, &contact.interfaces,
            binding.config.clone(), &names),
        None => RobinLinearization::new(cx, problem, binding.config.clone(), &names),
    }.map_err(producer)?;
    if solid.temperature().len() != inner.conduction.temperature.len()
        || solid.temperature().iter().zip(&inner.conduction.temperature)
            .any(|(a,b)| a.to_bits() != b.to_bits()) {
        return Err(producer("radiation adjoint reconstruction did not reproduce the accepted solid field"));
    }
    let walls = solid.wall_means(cx, solid.temperature()).map_err(producer)?;
    let air = network.linearize(cx, &walls).map_err(producer)?;
    let mut seed = CoupledObjective { nodal_temperatures: vec![0.0; solid.temperature().len()],
        wall_temperatures: vec![0.0; names.len()], solid_heat_rates: vec![0.0; names.len()],
        air: air.zero_objective() };
    objective.seed(&mut seed);
    let point = feedback::Point { htc: binding.htc.clone(), references: binding.references.clone(),
        driving: binding.driving.clone(), radiative_htc: binding.radiative_htc.clone() };
    let result = feedback::pullback(policy, request, cx, network, &solid, &point, references, htc,
        &seed.nodal_temperatures, &seed.wall_temperatures, derivations)?;
    let mut patches = Vec::new();
    for derivative in &result.patches {
        let patch = &policy.patches[&derivative.surface];
        patches.push(format!("{{\"surface\":{},\"dobjective_dlog_emissivity_k\":{},\"dobjective_demissivity_k\":{},\"dobjective_dambient_temperature\":{}}}",
            quote(&derivative.surface), num(derivative.log_emissivity)?,
            num(finite(derivative.log_emissivity/patch.emissivity.value(), "absolute emissivity derivative")?)?,
            num(derivative.ambient_temperature)?));
    }
    let gradient = result.thermal;
    let report = format!("{{\"method\":\"implicit-wall-feedback-adjoint\",\"iterations\":{},\"equation_residual\":{},\"equation_threshold\":{},\"reconstruction_solid_solves\":1,\"surfaces\":[{}],\"scope\":\"total steady derivative of the admitted mean-patch Robin closure, including k-prime, contacts, radiation-coefficient feedback and mixed-air feedback; full consistent face integrals, not products of mean jumps; one exact-input reconstruction must reproduce the accepted field bits; no iteration-history AD, pointwise T^4 substitution, transient derivative or error certificate; maximum objectives retain selected-active-vertex semantics; emissivity bounds restrict admissible perturbation directions\"}}",
        gradient.iterations, num(gradient.interface_residual)?, num(result.threshold)?, patches.join(","));
    poll(cx)?;
    Ok((gradient,report))
}
