//! Implicit closed-enclosure exchange at one physical time endpoint.
//! The callback owns storage, materials and contacts. Every invocation uses
//! the SAME old temperature field; radiosity iterations never advance time.
use super::*;
use fs_conduction::transient::backward_euler::StepSolution;

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
        // exchange recomputes radiosity at the NEW solved means and checks
        // both the unrelaxed temperature residual and every patch's watts.
        // The old field is only a starting guess for this boundary iteration.
        let accepted = exchange(policy, request, cx, &enclosure, names,
            references, htc, driving, solve)?;
        poll(cx)?;
        Ok((accepted.solid, accepted.states))
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
            // Convection was evaluated at the PHYSICAL air reference. The
            // difference is the uniform radiation flux actually assembled.
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
        // No inner-iteration count is reconstructed from the final field.
        // The time driver already counts ALL callback solves, including trials.
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
