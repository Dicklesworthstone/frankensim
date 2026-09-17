//! Radiation-aware reconstruction for the existing chronological adjoint tape.
//! No per-frame matrix or extra physical state is retained. Replay reconstructs
//! the inner boundary iteration at the saved air references, not the whole
//! coupled trajectory. Every replayed FEM solve is counted as reverse work.
use super::*;
use fs_conduction::transient::backward_euler::{
    BackwardEuler, NonlinearStepConfig, StepConfig, StepLinearization,
};

enum Point {
    Reservoir(feedback::Point),
    Enclosure { shifted_references: Vec<f64> },
}

pub(in crate::network_command) struct Reconstructed<'a> {
    pub step: StepLinearization<'a>,
    pub states: Vec<SolidRegionState>,
    pub solid_solves: usize,
    point: Point,
}

fn check_field(actual: &[f64], expected: &[f64]) -> Result<()> {
    if actual.len() != expected.len() || actual.iter().zip(expected)
        .any(|(a,b)| !a.is_finite() || a.to_bits() != b.to_bits()) {
        return Err(producer("radiative adjoint reconstruction changed accepted temperature bits"));
    }
    Ok(())
}

impl Policy {
    #[allow(clippy::too_many_arguments)]
    pub(in crate::network_command) fn reconstruct_endpoint<'a>(
        &self, request: &Request, cx: &Cx<'_>, engine: &'a BackwardEuler<'_>,
        network: &TransportNetwork<'_>, references: &[f64], htc: &BTreeMap<String,f64>,
        old: &[f64], source: &ScalarField, dt: f64, config: StepConfig,
        nonlinear: Option<NonlinearStepConfig>, expected: &[f64],
    ) -> Result<Reconstructed<'a>> {
        poll(cx)?;
        if let Some(enclosure)=&self.enclosure {
            let (step,states,shifted_references,solid_solves)=enclosure.reconstruct_endpoint(
                self,request,cx,engine,network,references,htc,old,source,dt,config,nonlinear,expected)?;
            return Ok(Reconstructed {step,states,solid_solves,point:Point::Enclosure {shifted_references}});
        }
        let names = network.regions();
        let material = fs_conduction::ConductivityModel::isotropic_declared(request.conductivity)
            .map_err(producer)?;
        let interfaces = request.contacts.as_ref().map(|contact| &contact.interfaces);
        let mut solid_solves = 0_usize;
        let (primal, states, rows) = self.advance_endpoint_bound(
            request,cx,&names,references,htc,old,|boundary| {
                let problem = ConductionProblem { mesh: &request.mesh, boundary, material: &material,
                    element_materials: request.solid_data.element_materials.as_ref(), source };
                let solved = match nonlinear {
                    Some(policy) => engine.advance_nonlinear(cx,problem,interfaces,old,dt,config,policy)
                        .map(|result| result.step),
                    None => engine.advance(cx,problem,interfaces,old,dt,config),
                }.map_err(|error| match error {
                    fs_conduction::ConductionError::NotConverged { .. } => Failure {
                        code: "cooling-network-transient-budget", message: error.to_string(),
                    },
                    other => producer(other),
                })?;
                solid_solves = solid_solves.checked_add(1)
                    .ok_or_else(|| producer("radiation reconstruction work overflow"))?;
                Ok(solved)
            },
        )?;
        check_field(&primal.temperature,expected)?;
        let boundary = request.boundary(&names,&rows.references,&rows.htc)?;
        let step = engine.linearize_step(cx,ConductionProblem { mesh: &request.mesh,
            boundary: &boundary, material: &material,
            element_materials: request.solid_data.element_materials.as_ref(), source },
            interfaces,old,dt,config,nonlinear,&names).map_err(producer)?;
        solid_solves = solid_solves.checked_add(1)
            .ok_or_else(|| producer("radiation reconstruction work overflow"))?;
        check_field(&step.primal().temperature,expected)?;
        self.endpoint_heat(request,cx,&states,step.primal())?;
        poll(cx)?;
        Ok(Reconstructed { step, states, solid_solves, point: Point::Reservoir(feedback::Point {
            htc: rows.htc, references: rows.references, driving: rows.driving,
            radiative_htc: rows.radiative_htc,
        }) })
    }

    /// Shared patch controls act throughout the complete fixed trajectory.
    pub(in crate::network_command) fn zero_trajectory_gradient(&self) -> BTreeMap<String,[f64;2]> {
        if let Some(enclosure)=&self.enclosure {return enclosure.zero_trajectory_gradient();}
        self.patches.keys().map(|name| (name.clone(),[0.0;2])).collect()
    }

    pub(in crate::network_command) fn trajectory_gradient_report(
        &self, gradients: &BTreeMap<String,[f64;2]>,
    ) -> Result<String> {
        if let Some(enclosure)=&self.enclosure {return enclosure.trajectory_gradient_report(gradients);}
        if gradients.len() != self.patches.len() {
            return Err(bad("trajectory radiation gradient has the wrong patch set"));
        }
        let mut rows = Vec::with_capacity(self.patches.len());
        for (name,patch) in &self.patches {
            let values = gradients.get(name).ok_or_else(|| bad("missing trajectory radiation derivative"))?;
            rows.push(format!("{{\"surface\":{},\"dtemperature_dlog_emissivity_k\":{},\"dtemperature_demissivity_k\":{},\"dtemperature_dambient_temperature\":{}}}",
                quote(name),num(values[0])?,num(finite(values[0]/patch.emissivity.value(),
                    "absolute trajectory emissivity derivative")?)?,num(values[1])?));
        }
        Ok(format!("{{\"method\":\"implicit-mean-patch-feedback\",\"surfaces\":[{}],\"scope\":\"each constant patch emissivity or surroundings temperature changes throughout the complete trajectory, including earlier cycles; new-temperature radiation feedback is included before carrying the storage adjoint to prior fields; complete consistent-face Robin derivative, not pointwise T(x)^4 or a frozen coefficient; emissivity-domain boundaries restrict perturbation directions; no derivative of timestep selection, solver iterations, physical-validation or continuous-time maximum certificate\"}}",rows.join(",")))
    }
}

impl Reconstructed<'_> {
    /// Use the existing objective seed vocabulary with correctly sized vectors.
    pub(in crate::network_command) fn zero_objective(
        &self, cx: &Cx<'_>, network: &TransportNetwork<'_>,
    ) -> Result<CoupledObjective> {
        let walls = self.step.wall_means(cx,self.step.temperature()).map_err(producer)?;
        let air = network.linearize(cx,&walls).map_err(producer)?;
        Ok(CoupledObjective { nodal_temperatures: vec![0.0; self.step.temperature().len()],
            wall_temperatures: vec![0.0; self.step.ports().len()],
            solid_heat_rates: vec![0.0; self.step.ports().len()], air: air.zero_objective() })
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::network_command) fn pullback(
        &self, policy: &Policy, request: &Request, cx: &Cx<'_>, network: &TransportNetwork<'_>,
        references: &[f64], htc: &BTreeMap<String,f64>,
        weights: &CoupledObjective, derived: &[convection::Derived],
    ) -> Result<feedback::Gradient> {
        if weights.solid_heat_rates.iter().any(|&w| w != 0.0) {
            return Err(bad("radiative trajectory adjoints admit temperature objectives only"));
        }
        match &self.point {
            Point::Reservoir(point)=>feedback::pullback(policy,request,cx,network,&self.step,point,references,htc,
                &weights.nodal_temperatures,&weights.wall_temperatures,derived),
            Point::Enclosure {shifted_references}=>policy.enclosure.as_ref()
                .ok_or_else(||bad("reconstructed enclosure no longer matches its radiation policy"))?
                .pullback_endpoint(request,cx,network,&self.step,shifted_references,references,htc,weights,derived),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn replay_refuses_changed_bits_and_nonfinite_fields() {
        assert!(check_field(&[300.0],&[300.0]).is_ok());
        assert!(check_field(&[f64::from_bits(300.0_f64.to_bits()+1)],&[300.0]).is_err());
        assert!(check_field(&[300.0,301.0],&[300.0]).is_err());
        assert!(check_field(&[f64::NAN],&[f64::NAN]).is_err());
    }
}
