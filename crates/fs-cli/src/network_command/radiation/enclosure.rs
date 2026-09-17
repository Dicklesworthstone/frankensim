//! Closed, gray-diffuse patch exchange through the existing radiosity solver.
//! View factors are supplied, not inferred from this mesh. Each patch emits at
//! its area-mean temperature and applies a UNIFORM radiative flux. Convection
//! remains pointwise P1 Robin transfer to the actual air reference.
use super::*;
use fs_conduction::radiation::{GrayDiffuseEnclosure, RadiationSurface, RadiosityReport,
    ViewFactorEvidence, ViewFactorMatrix, ViewFactorTolerance};

mod transient;
mod radiosity_adjoint;
mod sensitivity;

#[derive(Debug)]
struct SurfaceSpec {
    name: String,
    epsilon: f64,
    source: String,
}

#[derive(Debug)]
pub(super) struct Enclosure {
    surfaces: Vec<SurfaceSpec>,
    factors: ViewFactorMatrix,
    evidence_json: String,
}

impl Enclosure {
    pub(super) fn parse(value: &J, available: &[Surface]) -> Result<Self> {
        object(value, &["surfaces", "view_factors", "row_sum_tolerance",
            "reciprocity_relative_tolerance", "evidence"], "radiation.enclosure")?;
        let rows = array(get(value,"surfaces")?,"enclosure.surfaces",64)?;
        if rows.len() < 2 { return Err(bad("a closed enclosure needs at least two surfaces")); }
        let mut surfaces = Vec::with_capacity(rows.len());
        let mut names = BTreeSet::new();
        let mut areas = Vec::with_capacity(rows.len());
        for row in rows {
            object(row,&["surface","emissivity","source"],"enclosure surface")?;
            let name = string(get(row,"surface")?,"enclosure surface name")?;
            let surface = available.iter().find(|s|s.name==name)
                .ok_or_else(||bad(format!("unknown enclosure cooling surface {name}")))?;
            if !names.insert(name.clone()) { return Err(bad("duplicate enclosure surface")); }
            let epsilon = positive(get(row,"emissivity")?,"enclosure emissivity")?;
            if epsilon > 1.0 { return Err(bad("enclosure emissivity must lie in (0,1]")); }
            let source = string(get(row,"source")?,"enclosure emissivity source")?;
            surfaces.push(SurfaceSpec { name,epsilon,source });
            areas.push(surface.area);
        }
        let matrix = array(get(value,"view_factors")?,"enclosure.view_factors",rows.len())?;
        if matrix.len()!=rows.len() { return Err(bad("one view-factor row per enclosure surface required")); }
        let factors = matrix.iter().map(|row| {
            let row = array(row,"view-factor row",rows.len())?;
            if row.len()!=rows.len() { return Err(bad("view-factor matrix must be square")); }
            row.iter().map(|v|number(v,"view factor")).collect::<Result<Vec<_>>>()
        }).collect::<Result<Vec<_>>>()?;
        let tolerance = ViewFactorTolerance {
            row_sum_abs: tolerance(get(value,"row_sum_tolerance")?)?,
            reciprocity_rel: tolerance(get(value,"reciprocity_relative_tolerance")?)?,
        };
        let evidence = get(value,"evidence")?;
        let (evidence,evidence_json) = match get(evidence,"kind")?.as_str() {
            Some("analytic") => {
                object(evidence,&["kind","geometry"],"view-factor evidence")?;
                let geometry = string(get(evidence,"geometry")?,"analytic geometry/formula")?;
                let json = format!("{{\"kind\":\"analytic\",\"geometry\":{}}}",quote(&geometry));
                (ViewFactorEvidence::Analytic {geometry},json)
            }
            Some("external-qmc") => {
                object(evidence,&["kind","seed","samples","generator"],"view-factor evidence")?;
                let seed = string(get(evidence,"seed")?,"QMC seed")?.parse::<u64>()
                    .map_err(|_|bad("QMC seed must be a decimal u64 string"))?;
                let samples = string(get(evidence,"samples")?,"QMC samples")?.parse::<u64>()
                    .map_err(|_|bad("QMC samples must be a positive decimal u64 string"))?;
                let generator = string(get(evidence,"generator")?,"QMC generator")?;
                let json = format!("{{\"kind\":\"external-qmc\",\"seed\":{},\"samples\":{},\"generator\":{}}}",
                    quote(&seed.to_string()),quote(&samples.to_string()),quote(&generator));
                (ViewFactorEvidence::ExternalQmc {seed,samples,generator},json)
            }
            _ => return Err(bad("enclosure evidence must explicitly name analytic geometry or external-qmc generation")),
        };
        // Canonicalize BOTH axes with the surface list. Reordering named
        // patches and their corresponding matrix rows/columns is not new physics.
        let mut order: Vec<_> = (0..surfaces.len()).collect();
        order.sort_by(|&a,&b|surfaces[a].name.cmp(&surfaces[b].name));
        let areas = order.iter().map(|&i|areas[i]).collect();
        let factors = order.iter().map(|&i|order.iter().map(|&j|factors[i][j]).collect()).collect();
        surfaces.sort_by(|a,b|a.name.cmp(&b.name));
        let factors = ViewFactorMatrix::admit(areas,factors,evidence,tolerance).map_err(producer)?;
        Ok(Self {surfaces,factors,evidence_json})
    }

    pub(super) fn bind(&self, request: &Request, cx: &Cx<'_>) -> Result<GrayDiffuseEnclosure> {
        let mut surfaces = Vec::with_capacity(self.surfaces.len());
        for spec in &self.surfaces {
            poll(cx)?;
            let surface = request.surfaces.iter().find(|s|s.name==spec.name)
                .ok_or_else(||bad("enclosure surface disappeared"))?;
            // The caller declared a constant emissivity; the inline card has
            // no invented measured uncertainty or temperature dependence.
            let epsilon = declared_emissivity(&spec.name,spec.epsilon,300.0,&spec.source)?;
            surfaces.push(RadiationSurface::new(&request.mesh,spec.name.clone(),
                |face|surface.faces.contains(&face.vertices),epsilon).map_err(producer)?);
        }
        GrayDiffuseEnclosure::new(surfaces,self.factors.clone()).map_err(producer)
    }

    pub(super) fn evaluate(&self, policy: &Policy, request: &Request, cx: &Cx<'_>,
        flow: &GraphSolution, declared: &BTreeMap<String,f64>, want_gradient: bool)
        -> Result<fan_speed::ThermalEvaluation> {
        let enclosure = self.bind(request,cx)?;
        let (htc,convection) = convection::resolve(request,cx,flow,declared)?;
        let network = request.transport(cx,flow,&htc)?;
        let names = network.regions();
        let material = fs_conduction::ConductivityModel::isotropic_declared(request.conductivity).map_err(producer)?;
        let uniform = ScalarField::Uniform(request.source);
        let source = request.solid_data.nodal_source.as_ref().unwrap_or(&uniform);
        let gate = ConjugateConfig { max_iterations:request.limits.coupling,
            temperature_tolerance_k:request.limits.temperature,
            balance_tolerance_w:request.limits.heat/(names.len() as f64+1.0),
            balance_relative_tolerance:0.0,
            relaxation:Relaxation::Fixed {omega:request.limits.relaxation} };
        let mut last = None;
        let mut failure = None;
        let mut solid_solves = 0_usize;
        let coupled = solve_coupled_transport(cx,&network,&gate,|cx,references| {
            let initial = references.iter().try_fold(0.0,|s,&v|
                finite(s+v/references.len() as f64,"enclosure initial temperature"));
            let result = initial.and_then(|initial| {
                let driving = vec![initial;self.surfaces.len()];
                exchange(policy,request,cx,&enclosure,&names,references,&htc,driving,|boundary| {
                    let mut config = SolveConfig::default();
                    config.initial = InitialGuess::Uniform(initial);
                    config.linear.tolerance = request.limits.relative;
                    config.linear.max_iterations = request.limits.linear;
                    config.stop.residual_rtol = request.limits.relative;
                    config.stop.step_atol = 0.0;
                    let problem = ConductionProblem {mesh:&request.mesh,boundary,material:&material,
                        element_materials:request.solid_data.element_materials.as_ref(),source};
                    match &request.contacts {
                        Some(c)=>fs_conduction::solve::solve_with_interfaces(cx,problem,&c.interfaces,config),
                        None=>fs_conduction::solve::solve(cx,problem,config),
                    }.map_err(producer)
                })
            });
            match result {
                Ok(result) => {
                    solid_solves = match solid_solves.checked_add(result.iterations) {
                        Some(n)=>n,
                        None=>{failure=Some(producer("enclosure solve count overflow"));
                            return Err(AirflowError::Cancelled {iteration:0,references_k:references.to_vec()});}
                    };
                    let states = result.states.clone(); last=Some(result); Ok(states)
                }
                Err(error)=>{failure=Some(error);Err(AirflowError::Cancelled {iteration:0,references_k:references.to_vec()})}
            }
        });
        if let Some(error)=failure {return Err(error);}
        let coupled=coupled.map_err(producer)?;
        let result=last.ok_or_else(||bad("enclosure coupling returned without a solid"))?;
        for derived in &convection {derived.check_direction(&coupled.solid,request.limits.heat)?;}
        let source=finite(result.solid.report.energy.source_w,"enclosure source")?;
        let robin=finite(result.solid.report.energy.robin_out_w,"enclosure boundary heat")?;
        let convective=coupled.solid.iter().try_fold(0.0,|s,row|finite(s+row.heat_rate_w,"enclosure convective total"))?;
        let radiation=result.radiosity.enclosure_energy_closure_w;
        for residual in [source-robin,source-convective-radiation,
            convective-coupled.transport.external_heat_gain_w] {
            if finite(residual,"enclosure whole-system balance")?.abs()>request.limits.heat {
                return Err(producer("enclosure source, solid boundary and air energy do not close"));
            }
        }
        if let Some(power)=&request.solid_data.power {
            if (source-power.delivered_total_w()).abs()>request.limits.heat {
                return Err(producer("enclosure source differs from component power map"));
            }
        }
        let objective_state=request.objective.evaluate(cx,&result.solid.temperature,&coupled.solid)?;
        let (gradient,adjoint)=if want_gradient {
            let (gradient,report)=sensitivity::pullback(self,request,cx,&network,&enclosure,&result,
                &coupled.reference_temperatures_k,&htc,&objective_state,&convection)?;
            (Some(gradient),if request.gradient {report} else {"null".into()})
        } else {(None,"null".into())};
        let reconstruction_solves=usize::from(gradient.is_some());
        let total_solves=solid_solves.checked_add(reconstruction_solves)
            .ok_or_else(||producer("enclosure total solve count overflow"))?;
        let report=self.report(&result.radiosity,&result.applied_w,result.max_mismatch_w,
            Some(result.iterations),Some(total_solves))?;
        let prefix=report.strip_suffix('}').ok_or_else(||bad("enclosure report framing"))?;
        let report=format!("{prefix},\"forward_solid_solves\":{solid_solves},\"reconstruction_solid_solves\":{reconstruction_solves},\"adjoint\":{adjoint}}}");
        let contact_fluxes=request.contacts.as_ref().map(|c|c.interfaces.fluxes(&result.solid.temperature)
            .map_err(producer)).transpose()?.unwrap_or_default();
        let ordered_htc=names.iter().map(|name|htc[*name]).collect();
        poll(cx)?;
        Ok(fan_speed::ThermalEvaluation {solid_solves:total_solves,radiation:Some(report),value:Evaluation {
            coupled,temperatures:result.solid.temperature,gradient,objective:objective_state.value,
            objective_state,robin_total_w:robin,source_total_w:source,htc:ordered_htc,convection,contact_fluxes,
        }})
    }

    pub(super) fn report(&self, report: &RadiosityReport, applied: &[f64], mismatch: f64,
        iterations: Option<usize>, solves: Option<usize>) -> Result<String> {
        if applied.len()!=self.surfaces.len() {return Err(bad("enclosure report arity mismatch"));}
        let rows=self.surfaces.iter().enumerate().map(|(i,s)|Ok(format!(
            "{{\"surface\":{},\"emissivity\":{},\"area_m2\":{},\"mean_temperature_k\":{},\"radiosity_w_m2\":{},\"irradiation_w_m2\":{},\"outward_heat_w\":{},\"applied_heat_w\":{},\"source\":{}}}",
            quote(&s.name),num(s.epsilon)?,num(self.factors.areas_m2()[i])?,
            num(report.surface_temperatures_k[i])?,num(report.radiosity_w_m2[i])?,num(report.irradiation_w_m2[i])?,
            num(report.net_outward_heat_w[i])?,num(applied[i])?,quote(&s.source)))).collect::<Result<Vec<_>>>()?.join(",");
        let matrix=self.factors.factors().iter().map(|row|numbers(row)).collect::<Result<Vec<_>>>()?.join(",");
        Ok(format!("{{\"model\":\"closed-gray-diffuse-enclosure\",\"radiative_out_w\":{},\"max_nonlinear_mismatch_w\":{},\"radiosity_residual_w_m2\":{},\"iterations\":{},\"total_solid_solves\":{},\"view_factors\":[{}],\"view_factor_evidence\":{},\"max_reciprocity_residual\":{},\"surfaces\":[{}],\"scope\":\"closed fixed view-factor patch model; each area-mean temperature drives emission and each radiative flux is uniform on its patch; internal exchange, not heat lost to ambient or air; supplied analytic/QMC provenance is retained, not independently verified against mesh geometry; no occlusion computation, pointwise T(x)^4 law, radiation absorption in air, or physical validation\"}}",
            num(report.enclosure_energy_closure_w)?,num(mismatch)?,num(report.linear_residual_max_w_m2)?,
            iterations.map_or_else(||"null".into(),|n|n.to_string()),
            solves.map_or_else(||"null".into(),|n|n.to_string()),matrix,self.evidence_json,
            num(report.max_reciprocity_residual)?,rows))
    }
}

fn tolerance(value: &J) -> Result<f64> {
    let value=number(value,"view-factor tolerance")?;
    if !(0.0..=1e-8).contains(&value) {return Err(bad("view-factor tolerances must lie in [0,1e-8]"));}
    Ok(value)
}

pub(super) trait Response {
    fn temperature(&self)->&[f64];
    fn fluxes(&self)->&[fs_conduction::RobinFlux];
    fn robin_out_w(&self)->f64;
}
impl Response for ConductionSolution {
    fn temperature(&self)->&[f64] {&self.temperature}
    fn fluxes(&self)->&[fs_conduction::RobinFlux] {&self.report.robin_fluxes}
    fn robin_out_w(&self)->f64 {self.report.energy.robin_out_w}
}

pub(super) struct Exchange<T> {
    pub solid:T,
    pub states:Vec<SolidRegionState>,
    pub radiosity:RadiosityReport,
    pub applied_w:Vec<f64>,
    /// Exact boundary inputs of the accepted callback, not recovered from
    /// rounded patch watts or recomputed at a slightly different temperature.
    pub shifted_references:Vec<f64>,
    pub max_mismatch_w:f64,
    pub iterations:usize,
}

/// All callback invocations solve the same physical time/state, if transient.
#[allow(clippy::too_many_arguments)]
pub(super) fn exchange<T:Response>(policy:&Policy,request:&Request,cx:&Cx<'_>,
    enclosure:&GrayDiffuseEnclosure,names:&[&str],references:&[f64],htc:&BTreeMap<String,f64>,
    mut driving:Vec<f64>,mut solve:impl FnMut(&ThermalBoundary)->Result<T>)->Result<Exchange<T>> {
    if names.len()!=references.len() {return Err(bad("enclosure reference arity mismatch"));}
    let budget=request.limits.heat/(enclosure.surfaces().len() as f64+1.0);
    for iteration in 0..policy.max_iterations {
        poll(cx)?;
        let applied=radiosity(cx,enclosure,&driving,budget)?;
        let mut shifted=references.to_vec();
        for (i,surface) in enclosure.surfaces().iter().enumerate() {
            let slot=names.iter().position(|name|*name==surface.name()).ok_or_else(||bad("enclosure patch has no cooling row"))?;
            // h*(T-(T_air-q_rad/h)) = h*(T-T_air)+q_rad pointwise.
            // This algebraic reference need not be a physical temperature;
            // it is NEVER sent to air transport or reported as its inlet.
            shifted[slot]=finite(references[slot]-applied.net_outward_flux_w_m2[i]/htc[surface.name()],
                "enclosure shifted Robin reference")?;
        }
        let boundary=request.boundary(names,&shifted,htc)?;
        let solid=solve(&boundary)?;
        if solid.temperature().iter().any(|t|!t.is_finite()||*t<=0.0) {
            return Err(producer("enclosure solid produced a nonpositive physical temperature"));
        }
        let actual=enclosure.surfaces().iter().map(|surface|
            surface.mean_temperature(&request.mesh,solid.temperature()).map_err(producer)).collect::<Result<Vec<_>>>()?;
        let recomputed=radiosity(cx,enclosure,&actual,budget)?;
        let mut states=Vec::with_capacity(names.len());
        let mut assembled=0.0;
        let mut max_mismatch=0.0_f64;
        for (slot,&name) in names.iter().enumerate() {
            let flux=solid.fluxes().iter().find(|f|f.region==name).ok_or_else(||bad("enclosure callback omitted a Robin row"))?;
            let air=finite(htc[name]*flux.area_m2*(flux.mean_wall_temperature_k-references[slot]),"enclosure convective flux")?;
            let expected=enclosure.surfaces().iter().position(|s|s.name()==name)
                .map_or(0.0,|i|applied.net_outward_heat_w[i]);
            if finite(flux.heat_rate_w-air-expected,"enclosure boundary decomposition")?.abs()>budget {
                return Err(producer("enclosure boundary heat does not split into convection and the applied patch flux"));
            }
            assembled=finite(assembled+flux.heat_rate_w,"enclosure assembled boundary sum")?;
            states.push(SolidRegionState {region:name.into(),area_m2:flux.area_m2,
                mean_wall_temperature_k:flux.mean_wall_temperature_k,heat_rate_w:air,
                mean_reference_temperature_k:Some(references[slot])});
        }
        if finite(assembled-solid.robin_out_w(),"enclosure Robin total")?.abs()>request.limits.heat {
            return Err(producer("enclosure Robin regions do not partition the solid boundary"));
        }
        let mut change=0.0_f64;
        for i in 0..actual.len() {
            change=change.max(finite(actual[i]-driving[i],"enclosure temperature residual")?.abs());
            max_mismatch=max_mismatch.max(finite(recomputed.net_outward_heat_w[i]-applied.net_outward_heat_w[i],
                "enclosure nonlinear heat residual")?.abs());
            driving[i]=finite((1.0-policy.relaxation)*driving[i]+policy.relaxation*actual[i],"relaxed enclosure temperature")?;
        }
        if change<=policy.tolerance_k && max_mismatch<=budget {
            poll(cx)?;
            return Ok(Exchange {solid,states,radiosity:recomputed,applied_w:applied.net_outward_heat_w,
                shifted_references:shifted,max_mismatch_w:max_mismatch,iterations:iteration+1});
        }
    }
    Err(Failure {code:"cooling-network-radiation-budget",message:format!(
        "closed enclosure did not meet both temperature and patch-watt gates in {} solid solves; no partial result published",policy.max_iterations)})
}

pub(super) fn radiosity(cx:&Cx<'_>,enclosure:&GrayDiffuseEnclosure,temperatures:&[f64],
    tolerance_w:f64)->Result<RadiosityReport> {
    poll(cx)?;
    let report=enclosure.solve(cx,temperatures).map_err(producer)?;
    for &v in report.radiosity_w_m2.iter().chain(&report.irradiation_w_m2) {
        if !v.is_finite() || v<0.0 {return Err(producer("nonphysical enclosure radiosity/irradiation"));}
    }
    for &q in &report.net_outward_heat_w {finite(q,"enclosure patch heat")?;}
    let area=enclosure.surfaces().iter().map(|s|s.area_m2()).fold(0.0_f64,f64::max);
    if finite(report.linear_residual_max_w_m2*area,"enclosure radiosity residual")?>tolerance_w
        || finite(report.enclosure_energy_closure_w,"closed enclosure heat sum")?.abs()>tolerance_w {
        return Err(producer("radiosity solve failed its equation residual or closed-enclosure energy gate"));
    }
    Ok(report)
}
