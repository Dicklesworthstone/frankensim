//! Gray surface-to-surroundings radiation alongside convective cooling.
//!
//! Each selected surface is a lumped radiating patch: its AREA-MEAN P1
//! temperature drives eps*sigma*A*(T_mean^4 - T_surroundings^4). A positive
//! secant Robin coefficient is recomputed until this law agrees with the
//! solved patch heat. This is not the integral of pointwise T(x)^4 for a
//! nonuniform surface, nor a fixed small-temperature-departure linearization.
//! The surroundings are a large isothermal black reservoir with view factor
//! one; no inter-surface exchange, occlusion or participating medium is modeled.
//!
//! Alternatively, `enclosure` selects explicitly supplied closed view factors
//! and gray-diffuse reflection between solved surfaces. Its uniform patch flux
//! and internal heat accounting remain distinct from the reservoir model.
//!
//! Existing FEM, material, contact, air transport and cancellation producers
//! remain the numerical owners. Radiative heat never enters an air branch.
//! Both models support implicit transient endpoints and total steady adjoints.
//! Fixed-grid trajectory adjoints retain each model's complete feedback.
use super::*;
use fs_conduction::{ConductionSolution, SurfaceEmissivity, STEFAN_BOLTZMANN_W_M2_K4,
    SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, AmbientRadiationPatch,
    AmbientRadiationConfig, solve_with_ambient_radiation};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, SelectionPolicy, UncertaintyModel};

mod endpoint;
mod sensitivity;
mod enclosure;

#[derive(Debug)]
struct Patch {
    surface: String,
    ambient_k: f64,
    emissivity: SurfaceEmissivity,
    source: String,
}

impl Patch {
    fn model(&self) -> Result<AmbientRadiationPatch> {
        AmbientRadiationPatch::new(&self.surface, self.emissivity.clone(), self.ambient_k)
            .map_err(producer)
    }

    fn coefficient(&self, temperature: f64) -> Result<f64> {
        self.model()?.secant_coefficient_w_m2_k(temperature).map_err(producer)
    }
}

#[derive(Debug)]
pub(super) struct Policy {
    patches: BTreeMap<String, Patch>,
    enclosure: Option<enclosure::Enclosure>,
    max_iterations: usize,
    tolerance_k: f64,
    relaxation: f64,
}

#[derive(Debug)]
struct PatchHeat {
    surface: String,
    mean_k: f64,
    secant_h: f64,
    applied_w: f64,
    nonlinear_w: f64,
}

struct Inner {
    conduction: ConductionSolution,
    convective: Vec<SolidRegionState>,
    heats: Vec<PatchHeat>,
    iterations: usize,
    max_change_k: f64,
    max_mismatch_w: f64,
    binding: Option<sensitivity::Binding>,
}

fn finite(value: f64, stage: &str) -> Result<f64> {
    if value.is_finite() { Ok(value) } else { Err(producer(format!("nonfinite {stage}"))) }
}

impl Policy {
    pub(super) fn parse(value: &J, _root: &J, surfaces: &[Surface]) -> Result<Self> {
        object(value, &["max_iterations", "temperature_tolerance_k", "relaxation", "surfaces", "enclosure"], "radiation")?;
        let max_iterations = count(get(value, "max_iterations")?, "radiation.max_iterations", 1000)?;
        let tolerance_k = positive(get(value, "temperature_tolerance_k")?, "radiation.temperature_tolerance_k")?;
        let relaxation = positive(get(value, "relaxation")?, "radiation.relaxation")?;
        if relaxation > 1.0 { return Err(bad("radiation.relaxation must be in (0,1]")); }
        if let Some(enclosure) = value.get("enclosure") {
            if value.get("surfaces").is_some() {
                return Err(bad("choose radiation.surfaces or radiation.enclosure, not both"));
            }
            return Ok(Self {patches:BTreeMap::new(),enclosure:Some(enclosure::Enclosure::parse(enclosure,surfaces)?),
                max_iterations,tolerance_k,relaxation});
        }
        // Mesh studies use this complete producer, preserving the patch names
        // and partition while refining their faces. The transient parser still
        // owns fixed-grid/fixed-cycle derivative admission.
        let rows = array(get(value, "surfaces")?, "radiation.surfaces", surfaces.len())?;
        if rows.is_empty() { return Err(bad("radiation.surfaces must be nonempty")); }
        let mut patches = BTreeMap::new();
        for row in rows {
            object(row, &["surface", "emissivity", "ambient_temperature_k", "source"], "radiation surface")?;
            let surface = string(get(row, "surface")?, "radiation.surface")?;
            if !surfaces.iter().any(|s| s.name == surface) {
                return Err(bad(format!("unknown radiation cooling surface {surface}")));
            }
            let epsilon = positive(get(row, "emissivity")?, "radiation.emissivity")?;
            if epsilon > 1.0 { return Err(bad("radiation.emissivity must be in (0,1]")); }
            let ambient_k = positive(get(row, "ambient_temperature_k")?, "radiation.ambient_temperature_k")?;
            let source = string(get(row, "source")?, "radiation.source")?;
            let emissivity = declared_emissivity(&surface, epsilon, ambient_k, &source)?;
            let patch = Patch { surface: surface.clone(), ambient_k, emissivity, source };
            patch.coefficient(ambient_k)?;
            if patches.insert(surface, patch).is_some() {
                return Err(bad("a cooling surface has multiple radiation owners"));
            }
        }
        Ok(Self { patches, enclosure:None, max_iterations, tolerance_k, relaxation })
    }

    #[allow(clippy::too_many_arguments)]
    fn solid(&self, request: &Request, cx: &Cx<'_>, names: &[&str], references: &[f64],
        htc: &BTreeMap<String, f64>, want_gradient: bool) -> Result<Inner>
    {
        if names.is_empty() || names.len() != references.len() {
            return Err(bad("radiation requires nonempty matching air-reference rows"));
        }
        let boundary = request.boundary(names, references, htc)?;
        let material = fs_conduction::ConductivityModel::isotropic_declared(request.conductivity).map_err(producer)?;
        let uniform = ScalarField::Uniform(request.source);
        let source = request.solid_data.nodal_source.as_ref().unwrap_or(&uniform);
        let mut config = SolveConfig::default();
        config.initial = InitialGuess::Uniform(references.iter().sum::<f64>() / references.len() as f64);
        config.linear.tolerance = request.limits.relative;
        config.linear.max_iterations = self.limits_linear(request);
        config.stop.residual_rtol = request.limits.relative;
        config.stop.step_atol = 0.0;
        let patches = self.patches.values().map(Patch::model).collect::<Result<Vec<_>>>()?;
        let problem = ConductionProblem { mesh: &request.mesh, boundary: &boundary,
            material: &material, element_materials: request.solid_data.element_materials.as_ref(), source };
        // The canonical .fsim producer and this lab use exactly the same
        // radiation/contact/material solve, heat split and convergence gates.
        let result = solve_with_ambient_radiation(cx, problem,
            request.contacts.as_ref().map(|contact| &contact.interfaces), &patches, config,
            AmbientRadiationConfig { max_iterations: self.max_iterations,
                temperature_tolerance_k: self.tolerance_k, balance_tolerance_w: request.limits.heat,
                balance_relative_tolerance: 0.0, relaxation: self.relaxation },
        ).map_err(|error| match error {
            fs_conduction::ConductionError::AmbientRadiationNotConverged { .. } => Failure {
                code: "cooling-network-radiation-budget", message: error.to_string(),
            },
            other => producer(other),
        })?;
        poll(cx)?;
        let convective = names.iter().map(|name| {
            result.convective_robin_fluxes.iter().find(|flux| flux.region == *name)
                .map(SolidRegionState::from_robin_flux)
                .ok_or_else(|| bad(format!("radiation solid response lacks {name}")))
        }).collect::<Result<Vec<_>>>()?;
        let heats = result.radiation.patches.iter().map(|row| {
            Ok(PatchHeat { surface: row.patch.region().to_string(),
                mean_k: row.mean_surface_temperature_k,
                secant_h: row.patch.secant_coefficient_w_m2_k(row.mean_surface_temperature_k).map_err(producer)?,
                applied_w: row.applied_heat_w, nonlinear_w: row.nonlinear_heat_w })
        }).collect::<Result<Vec<_>>>()?;
        let binding = if want_gradient {
            let mut combined_h = htc.clone();
            let mut combined_ref = references.to_vec();
            let mut driving = references.to_vec();
            let mut applied_h = vec![0.0; names.len()];
            for (index, &name) in names.iter().enumerate() {
                let slot = result.combined_boundary.region_names().iter().position(|region| region == name)
                    .ok_or_else(|| bad("missing accepted radiation boundary"))?;
                let ThermalBc::Robin { htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(reference) }
                    = &result.combined_boundary.conditions()[slot]
                else { return Err(bad("accepted radiation boundary is not uniform Robin")); };
                combined_h.insert(name.to_string(), *h);
                combined_ref[index] = *reference;
                if let Some(row) = result.radiation.patches.iter().find(|row| row.patch.region() == name) {
                    driving[index] = row.driving_temperature_k;
                    applied_h[index] = row.applied_coefficient_w_m2_k;
                }
            }
            Some(sensitivity::Binding { htc: combined_h, references: combined_ref, driving,
                radiative_htc: applied_h, config: result.final_solve_config })
        } else { None };
        Ok(Inner { conduction: result.conduction, convective, heats,
            iterations: result.radiation.iterations,
            max_change_k: result.radiation.max_temperature_change_k,
            max_mismatch_w: result.radiation.max_heat_mismatch_w, binding })
    }

    fn limits_linear(&self, request: &Request) -> usize { request.limits.linear }

    pub(super) fn solve(&self, request: &Request, cx: &Cx<'_>) -> Result<String> {
        poll(cx)?;
        if let Some(schedule) = &request.transient {
            return super::transient::solve(request, cx, schedule);
        }
        if let Some(design) = &request.fan_speed_design {
            return fan_speed::solve_with(request, cx, design, |cx,flow,htc,want_gradient|
                self.evaluate(request,cx,flow,htc,want_gradient));
        }
        let flow = request.flow(cx)?;
        let result = if let Some(design) = &request.design {
            let designed = design::solve_with(request, cx, &flow, design, |cx, flow, htc, gradient|
                self.evaluate(request, cx, flow, htc, gradient))?;
            design::attach(designed.passing.render(request, &flow)?, &designed)?
        } else {
            let declared = request.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
            let evaluated = self.evaluate(request,cx,&flow,&declared,request.gradient)?;
            evaluated.render(request,&flow)?
        };
        poll(cx)?;
        match &request.fan { Some(fan) => fan.attach(result,&flow,fan.speed_ratio), None => Ok(result) }
    }

    /// Complete physical evaluation for fan candidates and mesh studies alike.
    /// The optional total derivative may be consumed internally by a marker.
    pub(super) fn evaluate(&self, request: &Request, cx: &Cx<'_>, flow: &GraphSolution,
        declared: &BTreeMap<String,f64>, want_gradient: bool) -> Result<fan_speed::ThermalEvaluation> {
        poll(cx)?;
        if let Some(enclosure)=&self.enclosure {
            return enclosure.evaluate(self,request,cx,flow,declared,want_gradient);
        }
        let (htc, convection) = convection::resolve(request, cx, flow, declared)?;
        let network = request.transport(cx, flow, &htc)?;
        let names = network.regions();
        let gate = ConjugateConfig { max_iterations: request.limits.coupling,
            temperature_tolerance_k: request.limits.temperature, balance_tolerance_w: request.limits.heat,
            balance_relative_tolerance: 0.0, relaxation: Relaxation::Fixed { omega: request.limits.relaxation } };
        let mut failure = None;
        let mut last = None;
        let mut solid_solves = 0_usize;
        let coupled = solve_coupled_transport(cx, &network, &gate, |cx, references| {
            match self.solid(request, cx, &names, references, &htc, want_gradient) {
                Ok(inner) => {
                    solid_solves = match solid_solves.checked_add(inner.iterations) {
                        Some(n) => n,
                        None => { failure = Some(producer("radiation solve count overflow"));
                            return Err(AirflowError::Cancelled { iteration: 0, references_k: references.to_vec() }); }
                    };
                    let response = inner.convective.clone();
                    last = Some(inner);
                    Ok(response)
                }
                Err(error) => { failure = Some(error);
                    Err(AirflowError::Cancelled { iteration: 0, references_k: references.to_vec() }) }
            }
        });
        if let Some(error) = failure { return Err(error); }
        let coupled = coupled.map_err(producer)?;
        let inner = last.ok_or_else(|| bad("radiation coupling produced no field"))?;
        for derived in &convection { derived.check_direction(&coupled.solid, request.limits.heat)?; }
        let source = finite(inner.conduction.report.energy.source_w, "radiative source")?;
        let robin = finite(inner.conduction.report.energy.robin_out_w, "radiative Robin heat")?;
        let convective = finite(coupled.solid.iter().map(|s| s.heat_rate_w).sum(), "total convective heat")?;
        let radiative = finite(inner.heats.iter().map(|h| h.nonlinear_w).sum(), "total radiative heat")?;
        let balance = finite(source - convective - radiative, "radiative whole-system balance")?;
        if balance.abs() > request.limits.heat || finite(source - robin, "assembled whole-solid balance")?.abs() > request.limits.heat {
            return Err(producer("source does not close against air plus recomputed nonlinear radiation"));
        }
        if let Some(power) = &request.solid_data.power {
            if finite(source - power.delivered_total_w(), "radiative source map balance")?.abs() > request.limits.heat {
                return Err(producer("radiative source disagrees with component power map"));
            }
        }
        let objective_state = request.objective.evaluate(cx, &inner.conduction.temperature, &coupled.solid)?;
        let (gradient, adjoint) = if want_gradient {
            let (gradient, report) = sensitivity::pullback(self, request, cx, &network, &inner,
                &coupled.reference_temperatures_k, &htc, &objective_state, &convection)?;
            (Some(gradient), if request.gradient { report } else { "null".into() })
        } else { (None,"null".into()) };
        let reconstruction_solves = usize::from(gradient.is_some());
        let total_solves = solid_solves.checked_add(reconstruction_solves)
            .ok_or_else(|| producer("radiation work count overflow"))?;
        let temperatures = inner.conduction.temperature;
        let contact_fluxes = request.contacts.as_ref().map(|contact|
            contact.interfaces.fluxes(&temperatures).map_err(producer)).transpose()?.unwrap_or_default();
        let evaluated = Evaluation { objective: objective_state.value, objective_state,
            coupled, temperatures, gradient, robin_total_w: robin, source_total_w: source,
            htc: names.iter().map(|name| htc[*name]).collect(), convection, contact_fluxes };
        let rows = inner.heats.iter().map(|heat| {
            let patch = &self.patches[&heat.surface];
            Ok(format!("{{\"surface\":{},\"emissivity\":{},\"ambient_temperature_k\":{},\"mean_temperature_k\":{},\"secant_htc_w_m2_k\":{},\"applied_heat_w\":{},\"nonlinear_heat_w\":{},\"source\":{}}}",
                quote(&patch.surface), num(patch.emissivity.value())?, num(patch.ambient_k)?, num(heat.mean_k)?,
                num(heat.secant_h)?, num(heat.applied_w)?, num(heat.nonlinear_w)?, quote(&patch.source)))
        }).collect::<Result<Vec<_>>>()?.join(",");
        let report = format!("{{\"model\":\"surface-mean-gray-to-isothermal-surroundings\",\"radiative_out_w\":{},\"convective_out_w\":{},\"energy_residual_w\":{},\"solid_solves\":{total_solves},\"forward_solid_solves\":{solid_solves},\"reconstruction_solid_solves\":{reconstruction_solves},\"final_inner_iterations\":{},\"final_temperature_change_k\":{},\"max_nonlinear_heat_mismatch_w\":{},\"surfaces\":[{}],\"adjoint\":{adjoint},\"scope\":\"caller-declared constant gray emissivity; fourth power of each surface's area-mean temperature, not pointwise T^4 integration; unit view factor to an isothermal black reservoir; convection and radiation share the declared faces but only convective heat enters the air; requested steady derivatives include radiative and air feedback; no enclosure reflection, occlusion, participating medium or physical validation\"}}",
            num(radiative)?, num(convective)?, num(balance)?, inner.iterations,
            num(inner.max_change_k)?, num(inner.max_mismatch_w)?, rows);
        poll(cx)?;
        Ok(fan_speed::ThermalEvaluation { value:evaluated, radiation:Some(report), solid_solves:total_solves })
    }
}

fn declared_emissivity(surface: &str, value: f64, temperature: f64, source: &str) -> Result<SurfaceEmissivity> {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim { key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Scalar { value, dims: EMISSIVITY_DIMS },
        validity: ValidityDomain::unconstrained(), uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity, observations: Vec::new(),
        provenance: Provenance { source: format!("caller-declared cooling radiation {surface}: {source}"),
            license: "unspecified".into(), artifact: None } }).map_err(producer)?;
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "caller-unspecified".into(),
        phase: "solid".into(), process: "caller-declared surface finish".into(), revision: 0 },
        claims, Vec::new()).map_err(producer)?;
    SurfaceEmissivity::from_card(surface, &card, temperature, SelectionPolicy::SingleClaimOnly).map_err(producer)
}
