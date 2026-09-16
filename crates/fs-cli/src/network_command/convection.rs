//! Flow-derived duct coefficients using fs-convection's existing validity gates.
//! Geometry and transport properties are declarations. Every coefficient is
//! frozen during one solid/air solve and recomputed when hydraulic flow changes.

use super::*;
use fs_convection::{CorrelationId, CorrelationInputs, NusseltEvaluation, ThermalConductivity,
    ThermalDirection, EvaluationSourceRegion};
use fs_qty::Length;

#[derive(Debug, Clone)]
pub(super) struct Law {
    card: CorrelationId,
    diameter: f64,
    area: f64,
    length: f64,
    viscosity: f64,
    conductivity: f64,
    aspect: Option<f64>,
    direction: Option<ThermalDirection>,
    source: String,
}

impl Law {
    fn parse(value: &J) -> Result<Self> {
        object(value, &["card", "hydraulic_diameter_m", "flow_area_m2", "channel_length_m",
            "dynamic_viscosity_pa_s", "thermal_conductivity_w_m_k", "aspect_ratio",
            "thermal_direction", "source"], "surface.convection")?;
        let name = string(get(value, "card")?, "convection.card")?;
        let card = CorrelationId::ALL.into_iter().find(|id| id.name() == name)
            .ok_or_else(|| bad("unknown convection card"))?;
        if !matches!(card, CorrelationId::CircularDuctLaminarCwt | CorrelationId::CircularDuctHausen
            | CorrelationId::RectangularDuctLaminarCwt | CorrelationId::RectangularDuctLaminarCwtDevelopingPr072
            | CorrelationId::DittusBoelter | CorrelationId::Gnielinski) {
            return Err(bad("this adapter admits forced duct CWT, Hausen, Dittus-Boelter or Gnielinski cards, not external/natural/constant-flux laws"));
        }
        let rectangular = matches!(card, CorrelationId::RectangularDuctLaminarCwt
            | CorrelationId::RectangularDuctLaminarCwtDevelopingPr072);
        let aspect = value.get("aspect_ratio").map(|v| positive(v, "aspect_ratio")).transpose()?;
        if rectangular != aspect.is_some() || aspect.is_some_and(|a| a > 1.0) {
            return Err(bad("rectangular cards require aspect_ratio in (0,1]; other duct cards must omit it"));
        }
        let direction = value.get("thermal_direction").map(|v| match v.as_str() {
            Some("heating-fluid") => Ok(ThermalDirection::HeatingFluid),
            Some("cooling-fluid") => Ok(ThermalDirection::CoolingFluid),
            _ => Err(bad("thermal_direction must be heating-fluid or cooling-fluid")),
        }).transpose()?;
        if (card == CorrelationId::DittusBoelter) != direction.is_some() {
            return Err(bad("thermal_direction is required exactly for Dittus-Boelter"));
        }
        Ok(Self { card, aspect, direction,
            diameter: positive(get(value, "hydraulic_diameter_m")?, "hydraulic_diameter_m")?,
            area: positive(get(value, "flow_area_m2")?, "flow_area_m2")?,
            length: positive(get(value, "channel_length_m")?, "channel_length_m")?,
            viscosity: positive(get(value, "dynamic_viscosity_pa_s")?, "dynamic_viscosity_pa_s")?,
            conductivity: positive(get(value, "thermal_conductivity_w_m_k")?, "thermal_conductivity_w_m_k")?,
            source: string(get(value, "source")?, "convection.source")?,
        })
    }
}

/// The inactive declared-h slot is zero only for an explicitly derived law.
/// Zero is never a usable coefficient: resolve() must replace it before FEM.
pub(super) fn parse(surface: &J) -> Result<(f64, Option<Law>)> {
    match (surface.get("htc_w_m2_k"), surface.get("convection")) {
        (Some(h), None) => Ok((positive(h, "htc_w_m2_k")?, None)),
        (None, Some(law)) => Ok((0.0, Some(Law::parse(law)?))),
        _ => Err(bad("every surface requires exactly one of htc_w_m2_k or convection")),
    }
}

#[derive(Debug)]
pub(super) struct Derived {
    pub(super) surface: String,
    branch: String,
    flow: f64,
    velocity: f64,
    law: Law,
    pub(super) nu: NusseltEvaluation,
    h: f64,
}

pub(super) fn resolve(request: &Request, cx: &Cx<'_>, flow: &GraphSolution,
    declared: &BTreeMap<String, f64>) -> Result<(BTreeMap<String, f64>, Vec<Derived>)> {
    if declared.len() != request.surfaces.len() || flow.branches.len() != request.region_paths.len() {
        return Err(bad("coefficient or hydraulic branch arity mismatch"));
    }
    let mut coefficients = BTreeMap::new();
    let mut derivations = Vec::new();
    for surface in &request.surfaces {
        poll(cx)?;
        let supplied = *declared.get(&surface.name).ok_or_else(|| bad("missing surface coefficient slot"))?;
        let Some(law) = &surface.convection else {
            if !(supplied.is_finite() && supplied > 0.0) { return Err(bad("declared coefficient must be finite and positive")); }
            coefficients.insert(surface.name.clone(), supplied);
            continue;
        };
        if supplied != 0.0 { return Err(bad("a derived convection coefficient cannot be overridden by effective-h sizing")); }
        let index = request.region_paths.iter().position(|path| path.contains(&surface.name))
            .ok_or_else(|| bad("convection surface has no hydraulic owner"))?;
        let branch = &flow.branches[index];
        let velocity = branch.flow.value().abs() / law.area;
        let reynolds = request.air.density.value() * velocity * law.diameter / law.viscosity;
        // Pr uses the SAME cp that determines transported heat capacity.
        let prandtl = law.viscosity * request.air.specific_heat_j_kg_k / law.conductivity;
        let mut groups = CorrelationInputs::forced(reynolds, prandtl).with_length_ratio(law.length / law.diameter);
        if let Some(aspect) = law.aspect { groups = groups.with_aspect_ratio(aspect); }
        if let Some(direction) = law.direction { groups = groups.with_direction(direction); }
        let nu = fs_convection::evaluate(law.card, groups).map_err(|e| producer(format!("surface {} on branch {}: {e}", surface.name, branch.loss.name)))?;
        if !nu.evidence().model.in_domain { return Err(producer("convection card returned an out-of-domain result")); }
        let h = nu.heat_transfer_coefficient(ThermalConductivity::new(law.conductivity), Length::new(law.diameter))
            .map_err(producer)?.value.value();
        coefficients.insert(surface.name.clone(), h);
        derivations.push(Derived { surface: surface.name.clone(), branch: branch.loss.name.clone(),
            flow: branch.flow.value(), velocity, law: law.clone(), nu, h });
    }
    poll(cx)?;
    Ok((coefficients, derivations))
}

impl Derived {
    pub(super) fn check_direction(&self, states: &[SolidRegionState], tolerance: f64) -> Result<()> {
        if let Some(direction) = self.law.direction {
            let heat = states.iter().find(|s| s.region == self.surface)
                .ok_or_else(|| bad("derived coefficient has no solid result"))?.heat_rate_w;
            let disagrees = match direction { ThermalDirection::HeatingFluid => heat < -tolerance,
                ThermalDirection::CoolingFluid => heat > tolerance };
            if disagrees { return Err(producer(format!("surface {} solved opposite to its declared Dittus-Boelter thermal direction", self.surface))); }
        }
        Ok(())
    }

    pub(super) fn render(&self) -> Result<String> {
        let groups = self.nu.groups().iter().map(|(name, value)| Ok(format!("{}:{}", quote(name), num(*value)?)))
            .collect::<Result<Vec<_>>>()?.join(",");
        let region = match self.nu.source_region() {
            EvaluationSourceRegion::CitedFormula => "cited-formula",
            EvaluationSourceRegion::TableSliceInterpolation => "table-slice-interpolation",
            EvaluationSourceRegion::DeclaredEngineeringBridge => "declared-engineering-bridge",
        };
        Ok(format!("{{\"surface\":{},\"branch\":{},\"card\":{},\"flow_m3_s\":{},\"speed_m_s\":{},\"groups\":{{{groups}}},\"nusselt\":{},\"htc_w_m2_k\":{},\"hydraulic_diameter_m\":{},\"flow_area_m2\":{},\"channel_length_m\":{},\"dynamic_viscosity_pa_s\":{},\"thermal_conductivity_w_m_k\":{},\"input_source\":{},\"formula_source\":{},\"formula_source_id\":{},\"source_region\":{},\"discrepancy_basis\":{},\"thermal_direction\":{},\"in_domain\":true,\"uncertainty_propagated\":false}}",
            quote(&self.surface), quote(&self.branch), quote(self.law.card.name()), num(self.flow)?, num(self.velocity)?,
            num(self.nu.evidence().value)?, num(self.h)?, num(self.law.diameter)?, num(self.law.area)?, num(self.law.length)?,
            num(self.law.viscosity)?, num(self.law.conductivity)?, quote(&self.law.source),
            quote(self.nu.card().source.citation), quote(self.nu.card().source.identifier), quote(region),
            quote(&format!("{:?}", self.nu.card().discrepancy_basis)),
            self.law.direction.map_or_else(|| "null".into(), |d| quote(match d {
                ThermalDirection::HeatingFluid => "heating-fluid", ThermalDirection::CoolingFluid => "cooling-fluid" }))))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::tests::{with_cx, close};
    const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/fan-correlated-hotspot.json"));

    #[test]
    fn solved_branch_reynolds_drives_the_existing_hausen_card() {
        let r = Request::parse(FIXTURE).unwrap();
        with_cx(|cx| {
            let flow = r.flow(cx).unwrap();
            let slots = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
            let (_, derived) = resolve(&r, cx, &flow, &slots).unwrap();
            for d in &derived {
                let re = 1.2 * d.flow.abs() / d.law.area * d.law.diameter / d.law.viscosity;
                let pr = d.law.viscosity * 1007.0 / d.law.conductivity;
                let gz = re * pr / (d.law.length / d.law.diameter);
                let nu = 3.66 + 0.0668 * gz / (1.0 + 0.04 * gz.powf(2.0 / 3.0));
                close(d.h, nu * d.law.conductivity / d.law.diameter, 1e-8);
                J::parse(&d.render().unwrap()).unwrap();
            }
            let automatic = r.evaluate(cx, &flow, &slots, false).unwrap();
            let mut explicit = Request::parse(FIXTURE).unwrap();
            for surface in &mut explicit.surfaces {
                surface.h = derived.iter().find(|d| d.surface == surface.name).unwrap().h;
                surface.convection = None;
            }
            let h = explicit.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
            let manual = explicit.evaluate(cx, &flow, &h, false).unwrap();
            for (a, b) in automatic.temperatures.iter().zip(&manual.temperatures) { close(*a, *b, 1e-8); }
        });
    }

    #[test]
    fn convection_is_recomputed_for_speed_and_is_orientation_independent() {
        let mut r = Request::parse(FIXTURE).unwrap();
        with_cx(|cx| {
            let slots = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
            let fan = r.fan.as_ref().unwrap();
            let slow = fan.solve(cx, &r.graph, r.limits, 0.5).unwrap();
            let fast = fan.solve(cx, &r.graph, r.limits, 1.5).unwrap();
            let (a, _) = resolve(&r, cx, &slow, &slots).unwrap();
            let (b, _) = resolve(&r, cx, &fast, &slots).unwrap();
            assert!(b["first-face"] > a["first-face"]);
            let mut edges = r.graph.branches().to_vec();
            let edge = &mut edges[0]; std::mem::swap(&mut edge.from, &mut edge.to);
            r.region_paths[0].reverse(); r.graph = LossGraph::new(3, edges).unwrap();
            let reversed = r.fan.as_ref().unwrap().solve(cx, &r.graph, r.limits, 1.5).unwrap();
            let (c, _) = resolve(&r, cx, &reversed, &slots).unwrap();
            close(c["first-face"], b["first-face"], 1e-8);
        });
    }

    #[test]
    fn missing_conflicting_and_out_of_domain_coefficient_models_refuse() {
        assert!(Request::parse(&FIXTURE.replace("\"convection\": {", "\"htc_w_m2_k\": 50, \"convection\": {")).is_err());
        assert!(Request::parse(&FIXTURE.replace("convection.circular-duct-hausen-developing", "convection.churchill-chu-vertical-plate")).is_err());
        let mut r = Request::parse(FIXTURE).unwrap();
        with_cx(|cx| {
            let flow = r.flow(cx).unwrap();
            let mut slots = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect::<BTreeMap<_, _>>();
            slots.insert("first-face".into(), 99.0);
            assert!(resolve(&r, cx, &flow, &slots).is_err());
            slots.insert("first-face".into(), 0.0);
            r.surfaces[0].convection.as_mut().unwrap().area *= 1e-6;
            assert!(resolve(&r, cx, &flow, &slots).is_err());
        });
    }
}