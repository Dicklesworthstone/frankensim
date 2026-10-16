//! Total enclosure/air feedback transpose, not a frozen-radiosity derivative.
//! Only the accepted solid callback is reconstructed. Both coupled fixed-point
//! histories are excluded from differentiation; their physical equations remain.
use super::*;
use fs_conduction::adjoint::robin::RobinResponse;
use fs_airflow::graph::thermal::coupled_transport::sensitivity::CoupledObjective;
use fs_couple::iqn_ils::IqnIls;

pub(super) struct Gradient {
    pub thermal: fan_gradient::CoolingGradient,
    pub log_emissivities: Vec<f64>,
    pub threshold: f64,
    pub radiosity_iterations: usize,
    pub radiosity_residual: f64,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn pullback(owner:&Enclosure,request:&Request,cx:&Cx<'_>,
    network:&TransportNetwork<'_>,enclosure:&GrayDiffuseEnclosure,
    accepted:&Exchange<ConductionSolution>,references:&[f64],htc:&BTreeMap<String,f64>,
    objective:&objective::ObjectiveState,derivations:&[convection::Derived])
    ->Result<(fan_gradient::CoolingGradient,String)> {
    poll(cx)?;
    let names=network.regions();
    let boundary=request.boundary(&names,&accepted.shifted_references,htc)?;
    let material=fs_conduction::ConductivityModel::isotropic_declared(request.conductivity).map_err(producer)?;
    let uniform=ScalarField::Uniform(request.source);
    let source=request.solid_data.nodal_source.as_ref().unwrap_or(&uniform);
    let mut config=SolveConfig::default();
    // Exactly the initial guess/configuration used by the accepted callback.
    let initial=references.iter().try_fold(0.0,|s,&v|
        finite(s+v/references.len() as f64,"enclosure reconstruction initial"))?;
    config.initial=InitialGuess::Uniform(initial);
    config.linear.tolerance=request.limits.relative;
    config.linear.max_iterations=request.limits.linear;
    config.stop.residual_rtol=request.limits.relative;
    config.stop.step_atol=0.0;
    let problem=ConductionProblem {mesh:&request.mesh,boundary:&boundary,material:&material,
        element_materials:request.solid_data.element_materials.as_ref(),source};
    let solid=match &request.contacts {
        Some(contact)=>RobinLinearization::new_with_interfaces(cx,problem,&contact.interfaces,config,&names),
        None=>RobinLinearization::new(cx,problem,config,&names),
    }.map_err(producer)?;
    if solid.temperature().len()!=accepted.solid.temperature.len()
        || solid.temperature().iter().zip(&accepted.solid.temperature).any(|(a,b)|a.to_bits()!=b.to_bits()) {
        return Err(producer("enclosure adjoint reconstruction changed accepted temperature bits"));
    }
    let walls=solid.wall_means(cx,solid.temperature()).map_err(producer)?;
    let air=network.linearize(cx,&walls).map_err(producer)?;
    let mut seed=CoupledObjective {nodal_temperatures:vec![0.0;solid.temperature().len()],
        wall_temperatures:vec![0.0;names.len()],solid_heat_rates:vec![0.0;names.len()],
        air:air.zero_objective()};
    objective.seed(&mut seed);
    let result=response(request,cx,network,enclosure,&solid,&accepted.shifted_references,
        references,htc,&seed.nodal_temperatures,&seed.wall_temperatures,derivations)?;
    let rows=owner.surfaces.iter().zip(&result.log_emissivities).map(|(patch,&value)|Ok(format!(
        "{{\"surface\":{},\"dobjective_dlog_emissivity_k\":{},\"dobjective_demissivity_k\":{}}}",
        quote(&patch.name),num(value)?,num(finite(value/patch.epsilon,"absolute enclosure emissivity gradient")?)?)))
        .collect::<Result<Vec<_>>>()?.join(",");
    let report=format!("{{\"method\":\"implicit-enclosure-air-feedback-adjoint\",\"iterations\":{},\"equation_residual\":{},\"equation_threshold\":{},\"radiosity_krylov_iterations\":{},\"max_radiosity_relative_residual\":{},\"reconstruction_solid_solves\":1,\"surfaces\":[{}],\"scope\":\"steady temperature objective at the accepted discrete field; total material K-prime, contact, air mixing and reflected radiation feedback; fixed view-factor matrix and patch partition; exact accepted boundary rows must replay temperature bits; no differentiation of coupling iterations, ambient reservoir substitution, geometry derivative or error certificate; maximum objectives retain selected-active-vertex semantics and emissivity boundaries restrict perturbations\"}}",
        result.thermal.iterations,num(result.thermal.interface_residual)?,num(result.threshold)?,
        result.radiosity_iterations,num(result.radiosity_residual)?,rows);
    Ok((result.thermal,report))
}

/// Shared response-level equations. A storage response can use these equations
/// too, but trajectory admission/reconstruction is separately owned and gated.
#[allow(clippy::too_many_arguments)]
pub(super) fn response(request:&Request,cx:&Cx<'_>,network:&TransportNetwork<'_>,
    enclosure:&GrayDiffuseEnclosure,solid:&RobinResponse,shifted:&[f64],
    references:&[f64],htc:&BTreeMap<String,f64>,nodal_weights:&[f64],wall_weights:&[f64],
    derivations:&[convection::Derived])->Result<Gradient> {
    poll(cx)?;
    let names=network.regions();let n=names.len();
    if n==0 || references.len()!=n || shifted.len()!=n || wall_weights.len()!=n
        || solid.ports().len()!=n || nodal_weights.len()!=solid.temperature().len() {
        return Err(bad("enclosure derivative needs complete accepted ports and temperature objective"));
    }
    let walls=solid.wall_means(cx,solid.temperature()).map_err(producer)?;
    let air=network.linearize(cx,&walls).map_err(producer)?;
    for (i,port) in solid.ports().iter().enumerate() {
        let h=*htc.get(names[i]).ok_or_else(||bad("missing enclosure convection coefficient"))?;
        if port.name!=names[i] || port.reference_k.to_bits()!=shifted[i].to_bits()
            || port.htc_w_m2_k.to_bits()!=h.to_bits()
            || finite(air.primal().reference_temperatures_k[i]-references[i],"enclosure air binding")?.abs()>request.limits.temperature {
            return Err(producer("enclosure adjoint is not bound to the accepted solid/air rows"));
        }
    }
    let slots=enclosure.surfaces().iter().map(|s|names.iter().position(|n|*n==s.name())
        .ok_or_else(||bad("enclosure derivative patch has no air port"))).collect::<Result<Vec<_>>>()?;
    let temperatures=slots.iter().map(|&slot|walls[slot]).collect::<Vec<_>>();
    let radiation=radiosity_adjoint::Linearization::new(cx,enclosure,&temperatures,
        request.limits.heat/(slots.len() as f64+1.0)).map_err(producer)?;
    let mut current=vec![0.0;n];let zero_heat=vec![0.0;n];
    let mut accelerator=IqnIls::new(n,acceleration::POLICY).map_err(producer)?;
    let mut last=f64::INFINITY;
    let mut radiosity_iterations=0_usize;let mut radiosity_residual=0.0_f64;
    for iteration in 1..=request.limits.derivative {
        poll(cx)?;
        let weights=current.iter().zip(wall_weights).map(|(&a,&b)|finite(a+b,"enclosure objective seed"))
            .collect::<Result<Vec<_>>>()?;
        let gradient=solid.pullback(cx,nodal_weights,&weights,&zero_heat).map_err(producer)?;
        // r_shift = r_air - q_rad/h. Radiation flux is NOT air heat.
        let flux_weights=slots.iter().map(|&slot|finite(-gradient.references[slot]/htc[names[slot]],
            "enclosure flux pullback")).collect::<Result<Vec<_>>>()?;
        let reflected=radiation.pullback(cx,&flux_weights,request.limits.relative,request.limits.linear)
            .map_err(producer)?;
        radiosity_iterations=radiosity_iterations.checked_add(reflected.iterations)
            .ok_or_else(||producer("enclosure derivative work overflow"))?;
        radiosity_residual=radiosity_residual.max(reflected.relative_residual);
        let mut air_weights=air.zero_objective();
        air_weights.references=gradient.references.clone();
        let transported=air.pullback(cx,&air_weights).map_err(producer)?;
        let mut proposal=transported.walls.clone();
        for (&slot,&value) in slots.iter().zip(&reflected.temperatures) {
            proposal[slot]=finite(proposal[slot]+value,"enclosure wall feedback")?;
        }
        let mut residual=0.0_f64;let mut scale=0.0_f64;
        for (&old,&next) in current.iter().zip(&proposal) {
            residual=residual.max(finite(next-old,"enclosure transpose residual")?.abs());
            scale=scale.max(old.abs()).max(next.abs());
        }
        let threshold=finite(request.limits.relative+request.limits.relative*scale,"enclosure transpose threshold")?;
        last=residual;poll(cx)?;
        if residual<=threshold {
            let log_htc=(0..n).map(|i|finite(gradient.log_htc[i]
                +(references[i]-shifted[i])*gradient.references[i]+transported.log_conductances[i],
                "enclosure total convection derivative")).collect::<Result<Vec<_>>>()?;
            let flow_scale=transported.log_conductances.iter().try_fold(0.0,|sum,&v|
                finite(sum-v,"enclosure air capacity derivative"))?;
            let thermal=CoupledGradient {inlets:transported.inlets,log_htc,nodal_load:gradient.nodal_load,
                interface_adjoint:air_weights.references,interface_residual:residual,iterations:iteration};
            let thermal=fan_gradient::from_flow_response(request,cx,&names,derivations,thermal,flow_scale)?;
            return Ok(Gradient {thermal,log_emissivities:reflected.log_emissivities,threshold,
                radiosity_iterations,radiosity_residual});
        }
        if iteration<request.limits.derivative {
            current=accelerator.step(&current,&proposal,request.limits.relaxation).map_err(producer)?.values;
        }
    }
    Err(Failure {code:"cooling-network-radiation-budget",message:format!(
        "enclosure total adjoint exhausted {} sweeps; original wall-equation residual {last}; no partial gradient published",
        request.limits.derivative)})
}
