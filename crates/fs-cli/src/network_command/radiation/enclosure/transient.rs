//! Implicit closed-enclosure exchange at one physical time endpoint.
//! The callback owns storage, materials and contacts. Every invocation uses
//! the SAME old temperature field; radiosity iterations never advance time.
use super::*;
use fs_conduction::transient::backward_euler::{BackwardEuler,NonlinearStepConfig,
    StepConfig,StepLinearization,StepSolution};
use fs_airflow::graph::thermal::coupled_transport::sensitivity::CoupledObjective;
use crate::network_command::radiation::sensitivity::feedback;

impl Response for StepSolution {
    fn temperature(&self) -> &[f64] { &self.temperature }
    fn fluxes(&self) -> &[fs_conduction::RobinFlux] { &self.robin_fluxes }
    fn robin_out_w(&self) -> f64 { self.robin_out_w }
}

pub(in crate::network_command::radiation) struct EndpointExchangeHeat {
    pub outward_w: f64,
    pub applied_w: f64,
    pub max_mismatch_w: f64,
    pub report: String,
}

impl Enclosure {
    #[allow(clippy::too_many_arguments)]
    pub(in crate::network_command::radiation) fn advance_endpoint(
        &self, policy: &Policy, request: &Request, cx: &Cx<'_>,
        names: &[&str], references: &[f64], htc: &BTreeMap<String, f64>,
        old: &[f64], solve: impl FnMut(&ThermalBoundary) -> Result<StepSolution>,
    ) -> Result<(StepSolution, Vec<SolidRegionState>)> {
        let accepted=self.endpoint_with_rows(policy,request,cx,names,references,htc,old,solve)?;
        Ok((accepted.solid,accepted.states))
    }

    #[allow(clippy::too_many_arguments)]
    fn endpoint_with_rows(&self,policy:&Policy,request:&Request,cx:&Cx<'_>,
        names:&[&str],references:&[f64],htc:&BTreeMap<String,f64>,old:&[f64],
        solve:impl FnMut(&ThermalBoundary)->Result<StepSolution>) -> Result<Exchange<StepSolution>> {
        poll(cx)?;
        if old.len() != request.mesh.vertex_count()
            || old.iter().any(|t| !t.is_finite() || *t <= 0.0)
        {
            return Err(bad("enclosure endpoint requires a complete positive previous temperature field"));
        }
        let enclosure = self.bind(request, cx)?;
        let driving = enclosure.surfaces().iter().map(|surface|
            surface.mean_temperature(&request.mesh, old).map_err(producer))
            .collect::<Result<Vec<_>>>()?;
        let accepted = exchange(policy, request, cx, &enclosure, names,
            references, htc, driving, solve)?;
        poll(cx)?;
        Ok(accepted)
    }

    /// Re-evaluate the accepted boundary loop, then linearize its EXACT final
    /// shifted rows. Both fields must replay the tape before reverse work begins.
    #[allow(clippy::too_many_arguments,clippy::type_complexity)]
    pub(in crate::network_command::radiation) fn reconstruct_endpoint<'a>(
        &self,policy:&Policy,request:&Request,cx:&Cx<'_>,engine:&'a BackwardEuler<'_>,
        network:&TransportNetwork<'_>,references:&[f64],htc:&BTreeMap<String,f64>,
        old:&[f64],source:&ScalarField,dt:f64,config:StepConfig,
        nonlinear:Option<NonlinearStepConfig>,expected:&[f64],
    )->Result<(StepLinearization<'a>,Vec<SolidRegionState>,Vec<f64>,usize)> {
        poll(cx)?;
        let names=network.regions();
        let material=fs_conduction::ConductivityModel::isotropic_declared(request.conductivity).map_err(producer)?;
        let interfaces=request.contacts.as_ref().map(|contact|&contact.interfaces);
        let accepted=self.endpoint_with_rows(policy,request,cx,&names,references,htc,old,|boundary| {
            let problem=ConductionProblem {mesh:&request.mesh,boundary,material:&material,
                element_materials:request.solid_data.element_materials.as_ref(),source};
            match nonlinear {
                Some(policy)=>engine.advance_nonlinear(cx,problem,interfaces,old,dt,config,policy).map(|s|s.step),
                None=>engine.advance(cx,problem,interfaces,old,dt,config),
            }.map_err(|error|match error {
                fs_conduction::ConductionError::NotConverged {..}=>Failure {
                    code:"cooling-network-transient-budget",message:error.to_string()},
                other=>producer(other),
            })
        })?;
        check_field(&accepted.solid.temperature,expected)?;
        let boundary=request.boundary(&names,&accepted.shifted_references,htc)?;
        let step=engine.linearize_step(cx,ConductionProblem {mesh:&request.mesh,boundary:&boundary,
            material:&material,element_materials:request.solid_data.element_materials.as_ref(),source},
            interfaces,old,dt,config,nonlinear,&names).map_err(producer)?;
        check_field(&step.primal().temperature,expected)?;
        self.endpoint_heat(request,cx,&accepted.states,step.primal())?;
        let solves=accepted.iterations.checked_add(1).ok_or_else(||producer("enclosure reconstruction count overflow"))?;
        poll(cx)?;
        Ok((step,accepted.states,accepted.shifted_references,solves))
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::network_command::radiation) fn pullback_endpoint(
        &self,request:&Request,cx:&Cx<'_>,network:&TransportNetwork<'_>,
        step:&StepLinearization<'_>,shifted:&[f64],references:&[f64],htc:&BTreeMap<String,f64>,
        weights:&CoupledObjective,derived:&[convection::Derived],
    )->Result<feedback::Gradient> {
        let bound=self.bind(request,cx)?;
        let result=super::sensitivity::response(request,cx,network,&bound,step,shifted,references,htc,
            &weights.nodal_temperatures,&weights.wall_temperatures,derived)?;
        let patches=self.surfaces.iter().zip(result.log_emissivities).map(|(surface,value)|
            feedback::PatchGradient {surface:surface.name.clone(),log_emissivity:value,
                // Internal shared accumulator slot only; the enclosure report
                // has no ambient control and verifies this slot stays zero.
                ambient_temperature:0.0}).collect();
        Ok(feedback::Gradient {thermal:result.thermal,patches,threshold:result.threshold})
    }

    pub(in crate::network_command::radiation) fn zero_trajectory_gradient(&self)->BTreeMap<String,[f64;2]> {
        self.surfaces.iter().map(|s|(s.name.clone(),[0.0;2])).collect()
    }

    pub(in crate::network_command::radiation) fn trajectory_gradient_report(
        &self,gradients:&BTreeMap<String,[f64;2]>,
    )->Result<String> {
        if gradients.len()!=self.surfaces.len() {return Err(bad("enclosure trajectory gradient patch count changed"));}
        let rows=self.surfaces.iter().map(|s| {
            let values=gradients.get(&s.name).ok_or_else(||bad("missing enclosure trajectory finish gradient"))?;
            if values[1]!=0.0 {return Err(bad("closed enclosure has no surroundings-temperature control"));}
            Ok(format!("{{\"surface\":{},\"dtemperature_dlog_emissivity_k\":{},\"dtemperature_demissivity_k\":{}}}",
                quote(&s.name),num(values[0])?,num(finite(values[0]/s.epsilon,"absolute trajectory finish gradient")?)?))
        }).collect::<Result<Vec<_>>>()?.join(",");
        Ok(format!("{{\"method\":\"implicit-enclosure-air-feedback\",\"surfaces\":[{rows}],\"scope\":\"each constant emissivity changes throughout the complete fixed trajectory and all cycles; total new-temperature reflected radiation, solid material/storage and mixed-air feedback; fixed view-factor matrix and patch geometry; no ambient reservoir, timestep-selection derivative, unique max derivative at ties or continuous-time peak certificate\"}}"))
    }

    /// Independently re-evaluate heat from an accepted response. This also
    /// runs for adaptive trial endpoints, but only the enclosing time driver
    /// may accumulate the result into physical history.
    pub(in crate::network_command::radiation) fn endpoint_heat(
        &self, request: &Request, cx: &Cx<'_>, states: &[SolidRegionState],
        step: &StepSolution,
    ) -> Result<EndpointExchangeHeat> {
        poll(cx)?;
        if states.len() != request.surfaces.len() || states.len() != step.robin_fluxes.len()
            || step.temperature.len() != request.mesh.vertex_count()
            || step.temperature.iter().any(|t| !t.is_finite() || *t <= 0.0)
        {
            return Err(bad("enclosure endpoint has incomplete boundary rows or physical temperatures"));
        }
        let enclosure = self.bind(request, cx)?;
        let means = enclosure.surfaces().iter().map(|surface|
            surface.mean_temperature(&request.mesh, &step.temperature).map_err(producer))
            .collect::<Result<Vec<_>>>()?;
        let tolerance = request.limits.heat / (self.surfaces.len() as f64 + 1.0);
        let fresh = radiosity(cx, &enclosure, &means, tolerance)?;
        let mut applied = vec![0.0; self.surfaces.len()];
        let mut seen = BTreeSet::new();
        let mut assembled = 0.0;
        let mut applied_total = 0.0;
        let mut max_mismatch = 0.0_f64;
        for state in states {
            poll(cx)?;
            if !seen.insert(state.region.as_str()) {
                return Err(bad("duplicate enclosure endpoint boundary row"));
            }
            let flux = step.robin_fluxes.iter().find(|f| f.region == state.region)
                .ok_or_else(|| bad("enclosure endpoint boundary identity mismatch"))?;
            assembled = finite(assembled + flux.heat_rate_w, "enclosure endpoint boundary total")?;
            let actual = finite(flux.heat_rate_w - state.heat_rate_w,
                "enclosure endpoint applied radiation")?;
            if let Some(index) = self.surfaces.iter().position(|s| s.name == state.region) {
                applied[index] = actual;
                applied_total = finite(applied_total + actual, "enclosure endpoint applied heat sum")?;
                let mismatch = finite(actual - fresh.net_outward_heat_w[index],
                    "enclosure endpoint nonlinear patch mismatch")?.abs();
                max_mismatch = max_mismatch.max(mismatch);
            } else if actual.abs() > tolerance {
                return Err(producer("nonradiating endpoint patch has unexplained boundary heat"));
            }
        }
        if self.surfaces.iter().any(|s| !seen.contains(s.name.as_str())) {
            return Err(bad("enclosure endpoint omitted a radiating patch"));
        }
        if max_mismatch > tolerance
            || finite(assembled - step.robin_out_w, "enclosure endpoint boundary decomposition")?.abs()
                > request.limits.heat
            || finite(applied_total - fresh.enclosure_energy_closure_w,
                "enclosure endpoint internal radiation closure")?.abs() > request.limits.heat
        {
            return Err(producer("enclosure endpoint radiation or boundary energy does not close"));
        }
        let report = self.report(&fresh, &applied, max_mismatch, None, None)?;
        let prefix = report.strip_suffix('}').ok_or_else(|| bad("enclosure endpoint report framing"))?;
        let report = format!("{prefix},\"temporal_scope\":\"final-accepted-endpoint\",\"heat_semantics\":\"signed exchange between solved patches; radiative_out_w is the closed-enclosure residual, not an external heat sink\"}}");
        poll(cx)?;
        Ok(EndpointExchangeHeat {
            outward_w: fresh.enclosure_energy_closure_w,
            applied_w: applied_total,
            max_mismatch_w: max_mismatch,
            report,
        })
    }
}

fn check_field(actual:&[f64],expected:&[f64])->Result<()> {
    if actual.len()!=expected.len() || actual.iter().zip(expected)
        .any(|(a,b)|!a.is_finite() || a.to_bits()!=b.to_bits()) {
        return Err(producer("enclosure adjoint reconstruction changed accepted temperature bits"));
    }
    Ok(())
}
