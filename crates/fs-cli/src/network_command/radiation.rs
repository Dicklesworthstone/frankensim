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
//! Existing FEM, material, contact, air transport and cancellation producers
//! remain the numerical owners. Radiative heat never enters an air branch.
//! Transient endpoints use the same law implicitly with fixed old solid state.
use super::*;
use fs_conduction::{ConductionSolution, SurfaceEmissivity, STEFAN_BOLTZMANN_W_M2_K4,
    SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, SelectionPolicy, UncertaintyModel};

mod endpoint;

#[derive(Debug)]
struct Patch {
    surface: String,
    ambient_k: f64,
    emissivity: SurfaceEmissivity,
    source: String,
}

impl Patch {
    /// Factor the fourth-power difference rather than subtracting nearly equal
    /// fourth powers. At equal temperatures this is the finite tangent limit.
    fn coefficient(&self, temperature: f64) -> Result<f64> {
        if !(temperature.is_finite() && temperature > 0.0) {
            return Err(producer("radiation requires positive absolute surface temperature"));
        }
        let ambient = self.ambient_k;
        let value = self.emissivity.value() * STEFAN_BOLTZMANN_W_M2_K4
            * (temperature + ambient) * (temperature * temperature + ambient * ambient);
        if !(value.is_finite() && value > 0.0) {
            return Err(producer("radiation secant coefficient is not finite and positive"));
        }
        Ok(value)
    }
}

#[derive(Debug)]
pub(super) struct Policy {
    patches: BTreeMap<String, Patch>,
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
}

fn finite(value: f64, stage: &str) -> Result<f64> {
    if value.is_finite() { Ok(value) } else { Err(producer(format!("nonfinite {stage}"))) }
}

impl Policy {
    pub(super) fn parse(value: &J, root: &J, surfaces: &[Surface]) -> Result<Self> {
        object(value, &["max_iterations", "temperature_tolerance_k", "relaxation", "surfaces"], "radiation")?;
        // Forward transients share the implicit endpoint law and heat ledger.
        // Their derivatives do not: never supply a frozen-radiation adjoint.
        if ["mesh_convergence", "design", "fan_speed_design"].iter()
            .any(|key| root.get(key).is_some())
            || get(get(root, "objective")?, "gradient")? != &J::Bool(false)
        {
            return Err(bad("radiation requires a primal request without gradients, mesh studies or steady design searches"));
        }
        if root.get("transient").is_some_and(|schedule| schedule.get("adjoint").is_some()) {
            return Err(bad("radiative transient adjoints are not implemented; omit transient.adjoint rather than freezing radiation"));
        }
        let max_iterations = count(get(value, "max_iterations")?, "radiation.max_iterations", 1000)?;
        let tolerance_k = positive(get(value, "temperature_tolerance_k")?, "radiation.temperature_tolerance_k")?;
        let relaxation = positive(get(value, "relaxation")?, "radiation.relaxation")?;
        if relaxation > 1.0 { return Err(bad("radiation.relaxation must be in (0,1]")); }
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
        Ok(Self { patches, max_iterations, tolerance_k, relaxation })
    }

    /// One inner solve at fixed AIR references. The combined Robin reference is
    /// only an assembly device; the air callback receives convective heat alone.
    fn solid(&self, request: &Request, cx: &Cx<'_>, names: &[&str], references: &[f64],
        htc: &BTreeMap<String, f64>) -> Result<Inner>
    {
        if names.len() != references.len() { return Err(bad("radiation air-reference arity mismatch")); }
        let mut driving = references.to_vec();
        let material = fs_conduction::ConductivityModel::isotropic_declared(request.conductivity).map_err(producer)?;
        let uniform = ScalarField::Uniform(request.source);
        let source = request.solid_data.nodal_source.as_ref().unwrap_or(&uniform);
        let mut initial = InitialGuess::Uniform(references.iter().sum::<f64>() / references.len() as f64);
        for iteration in 0..self.max_iterations {
            poll(cx)?;
            let mut combined_h = htc.clone();
            let mut combined_ref = references.to_vec();
            let mut applied_h = vec![0.0; names.len()];
            for (i, &name) in names.iter().enumerate() {
                if let Some(patch) = self.patches.get(name) {
                    let h_rad = patch.coefficient(driving[i])?;
                    let h_air = htc[name];
                    let h = finite(h_air + h_rad, "combined radiation/convection coefficient")?;
                    // Convex weights avoid multiplying temperature by a large h.
                    combined_ref[i] = finite((h_air / h) * references[i]
                        + (h_rad / h) * patch.ambient_k, "combined Robin reference")?;
                    combined_h.insert(name.to_string(), h);
                    applied_h[i] = h_rad;
                }
            }
            let boundary = request.boundary(names, &combined_ref, &combined_h)?;
            let mut config = SolveConfig::default();
            config.initial = initial;
            config.linear.tolerance = request.limits.relative;
            config.linear.max_iterations = request.limits.linear;
            config.stop.residual_rtol = request.limits.relative;
            config.stop.step_atol = 0.0;
            let problem = ConductionProblem { mesh: &request.mesh, boundary: &boundary,
                material: &material, element_materials: request.solid_data.element_materials.as_ref(), source };
            let conduction = match &request.contacts {
                Some(contact) => fs_conduction::solve::solve_with_interfaces(cx, problem, &contact.interfaces, config),
                None => fs_conduction::solve::solve(cx, problem, config),
            }.map_err(producer)?;
            poll(cx)?;
            let mut convective = Vec::with_capacity(names.len());
            let mut heats = Vec::with_capacity(self.patches.len());
            let mut max_change = 0.0_f64;
            let mut max_mismatch = 0.0_f64;
            let mut combined_total = 0.0;
            for (i, &name) in names.iter().enumerate() {
                let flux = conduction.report.robin_fluxes.iter().find(|f| f.region == name)
                    .ok_or_else(|| bad(format!("radiation solid response lacks {name}")))?;
                let mean = flux.mean_wall_temperature_k;
                if !(mean.is_finite() && mean > 0.0) {
                    return Err(producer("radiation solid solve returned a nonpositive surface temperature"));
                }
                let q_air = finite(htc[name] * flux.area_m2 * (mean - references[i]), "convective heat")?;
                let mut q_applied = 0.0;
                if let Some(patch) = self.patches.get(name) {
                    let h = patch.coefficient(mean)?;
                    let q = finite(h * flux.area_m2 * (mean - patch.ambient_k), "nonlinear radiative heat")?;
                    q_applied = finite(applied_h[i] * flux.area_m2 * (mean - patch.ambient_k), "applied radiative heat")?;
                    max_change = max_change.max(finite(mean - driving[i], "radiation temperature residual")?.abs());
                    max_mismatch = max_mismatch.max(finite(q - q_applied, "radiation heat residual")?.abs());
                    heats.push(PatchHeat { surface: name.to_string(), mean_k: mean,
                        secant_h: h, applied_w: q_applied, nonlinear_w: q });
                    driving[i] = finite((1.0 - self.relaxation) * driving[i] + self.relaxation * mean,
                        "relaxed radiation temperature")?;
                }
                // Independently accumulated FEM Robin heat must split into the
                // exact two mechanisms, even before the nonlinear iteration closes.
                if finite(flux.heat_rate_w - q_air - q_applied, "Robin heat split")?.abs() > request.limits.heat {
                    return Err(producer("radiation/convection split disagrees with assembled Robin heat"));
                }
                combined_total = finite(combined_total + q_air + q_applied, "combined boundary heat")?;
                convective.push(SolidRegionState { region: name.to_string(), area_m2: flux.area_m2,
                    mean_wall_temperature_k: mean, heat_rate_w: q_air,
                    mean_reference_temperature_k: Some(references[i]) });
            }
            if finite(combined_total - conduction.report.energy.robin_out_w, "Robin decomposition")?.abs() > request.limits.heat {
                return Err(producer("radiative solid boundary decomposition failed"));
            }
            if max_change <= self.tolerance_k && max_mismatch <= request.limits.heat {
                return Ok(Inner { conduction, convective, heats, iterations: iteration + 1,
                    max_change_k: max_change, max_mismatch_w: max_mismatch });
            }
            // These models have only Robin/natural/contact rows; all nodal DOFs
            // are free. Reuse the field as a guess, never as accepted physics.
            initial = InitialGuess::Free(conduction.temperature);
        }
        Err(Failure { code: "cooling-network-radiation-budget", message: format!(
            "radiation did not satisfy both temperature and nonlinear watt gates within {} solid solves at the current air reference; no partial result published", self.max_iterations) })
    }

    pub(super) fn solve(&self, request: &Request, cx: &Cx<'_>) -> Result<String> {
        poll(cx)?;
        if let Some(schedule) = &request.transient {
            return super::transient::solve(request, cx, schedule);
        }
        let flow = request.flow(cx)?;
        let declared = request.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let (htc, convection) = convection::resolve(request, cx, &flow, &declared)?;
        let network = request.transport(cx, &flow, &htc)?;
        let names = network.regions();
        let gate = ConjugateConfig { max_iterations: request.limits.coupling,
            temperature_tolerance_k: request.limits.temperature, balance_tolerance_w: request.limits.heat,
            balance_relative_tolerance: 0.0, relaxation: Relaxation::Fixed { omega: request.limits.relaxation } };
        let mut failure = None;
        let mut last = None;
        let mut solid_solves = 0_usize;
        let coupled = solve_coupled_transport(cx, &network, &gate, |cx, references| {
            match self.solid(request, cx, &names, references, &htc) {
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
        let temperatures = inner.conduction.temperature;
        let objective_state = request.objective.evaluate(cx, &temperatures, &coupled.solid)?;
        let contact_fluxes = request.contacts.as_ref().map(|contact|
            contact.interfaces.fluxes(&temperatures).map_err(producer)).transpose()?.unwrap_or_default();
        let evaluated = Evaluation { objective: objective_state.value, objective_state,
            coupled, temperatures, gradient: None, robin_total_w: robin, source_total_w: source,
            htc: names.iter().map(|name| htc[*name]).collect(), convection, contact_fluxes };
        let result = render(request, &flow, &evaluated)?;
        let result = match &request.fan { Some(fan) => fan.attach(result, &flow, fan.speed_ratio)?, None => result };
        let rows = inner.heats.iter().map(|heat| {
            let patch = &self.patches[&heat.surface];
            Ok(format!("{{\"surface\":{},\"emissivity\":{},\"ambient_temperature_k\":{},\"mean_temperature_k\":{},\"secant_htc_w_m2_k\":{},\"applied_heat_w\":{},\"nonlinear_heat_w\":{},\"source\":{}}}",
                quote(&patch.surface), num(patch.emissivity.value())?, num(patch.ambient_k)?, num(heat.mean_k)?,
                num(heat.secant_h)?, num(heat.applied_w)?, num(heat.nonlinear_w)?, quote(&patch.source)))
        }).collect::<Result<Vec<_>>>()?.join(",");
        let prefix = result.strip_suffix("}\n").ok_or_else(|| bad("radiation result framing"))?;
        poll(cx)?;
        Ok(format!("{prefix},\"radiation\":{{\"model\":\"surface-mean-gray-to-isothermal-surroundings\",\"radiative_out_w\":{},\"convective_out_w\":{},\"energy_residual_w\":{},\"solid_solves\":{},\"final_inner_iterations\":{},\"final_temperature_change_k\":{},\"max_nonlinear_heat_mismatch_w\":{},\"surfaces\":[{}],\"scope\":\"caller-declared constant gray emissivity; fourth power of each surface's area-mean temperature, not pointwise T^4 integration; unit view factor to an isothermal black reservoir; convection and radiation share the declared faces but only convective heat enters the air; no enclosure reflection, occlusion, participating medium, radiation gradient or physical validation\"}}}}\n",
            num(radiative)?, num(convective)?, num(balance)?, solid_solves, inner.iterations,
            num(inner.max_change_k)?, num(inner.max_mismatch_w)?, rows))
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
