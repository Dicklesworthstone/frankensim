//! Implicit endpoint exchange for the declared radiation model.
//! The callback owns backward Euler, heat capacity, k(T), contacts and source.
//! Every callback invocation must use the SAME previous accepted physical field.
use super::*;
use fs_conduction::transient::backward_euler::StepSolution;

/// Exact rows of the accepted callback, not coefficients recomputed later at a
/// slightly different, finite-tolerance wall. Retained only during adjoint replay.
pub(in crate::network_command::radiation) struct SolvedRows {
    pub htc: BTreeMap<String, f64>,
    pub references: Vec<f64>,
    pub driving: Vec<f64>,
    pub radiative_htc: Vec<f64>,
}

/// Recomputed endpoint heat, never accumulated until its timestep is accepted.
pub(in crate::network_command) struct EndpointHeat {
    pub outward_w: f64,
    applied_w: f64,
    max_mismatch_w: f64,
    rows: Vec<PatchHeat>,
    enclosure_report: Option<String>,
}

impl Policy {
    /// Iterate only radiation at fixed air references. No trial temperature is
    /// returned as the previous physical state for another radiation iteration.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::network_command) fn advance_endpoint(
        &self, request: &Request, cx: &Cx<'_>, names: &[&str], references: &[f64],
        htc: &BTreeMap<String, f64>, old: &[f64],
        solve: impl FnMut(&ThermalBoundary) -> Result<StepSolution>,
    ) -> Result<(StepSolution, Vec<SolidRegionState>)> {
        if let Some(enclosure) = &self.enclosure {
            return enclosure.advance_endpoint(self, request, cx, names, references, htc, old, solve);
        }
        self.endpoint_inner(request,cx,names,references,htc,old,false,solve)
            .map(|(step,states,_)| (step,states))
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::network_command::radiation) fn advance_endpoint_bound(
        &self, request: &Request, cx: &Cx<'_>, names: &[&str], references: &[f64],
        htc: &BTreeMap<String, f64>, old: &[f64],
        solve: impl FnMut(&ThermalBoundary) -> Result<StepSolution>,
    ) -> Result<(StepSolution, Vec<SolidRegionState>, SolvedRows)> {
        if self.enclosure.is_some() {
            return Err(bad("enclosure-radiation trajectory adjoints are not implemented; no frozen-radiosity derivative is returned"));
        }
        let (step,states,rows) = self.endpoint_inner(request,cx,names,references,htc,old,true,solve)?;
        Ok((step,states,rows.ok_or_else(|| bad("radiative replay did not retain its accepted rows"))?))
    }

    #[allow(clippy::too_many_arguments)]
    fn endpoint_inner(
        &self, request: &Request, cx: &Cx<'_>, names: &[&str], references: &[f64],
        htc: &BTreeMap<String, f64>, old: &[f64], retain_rows: bool,
        mut solve: impl FnMut(&ThermalBoundary) -> Result<StepSolution>,
    ) -> Result<(StepSolution, Vec<SolidRegionState>, Option<SolvedRows>)> {
        if names.len() != references.len() || old.len() != request.mesh.vertex_count()
            || old.iter().any(|t| !t.is_finite() || *t <= 0.0)
        {
            return Err(bad("radiative endpoint requires complete positive physical history and air references"));
        }
        let mut driving = references.to_vec();
        for (i, &name) in names.iter().enumerate() {
            if !self.patches.contains_key(name) { continue; }
            let surface = request.surfaces.iter().find(|s| s.name == name)
                .ok_or_else(|| bad("unknown endpoint radiation surface"))?;
            let mut mean = 0.0;
            for face in request.mesh.boundary() {
                poll(cx)?;
                if surface.faces.contains(&face.vertices) {
                    for vertex in face.vertices {
                        mean = finite(mean + (face.area / surface.area / 3.0) * old[vertex as usize],
                            "radiation initial surface mean")?;
                    }
                }
            }
            self.patches[name].coefficient(mean)?;
            driving[i] = mean;
        }
        for _ in 0..self.max_iterations {
            poll(cx)?;
            let mut combined_h = htc.clone();
            let mut combined_ref = references.to_vec();
            let mut applied_h = vec![0.0; names.len()];
            for (i, &name) in names.iter().enumerate() {
                if let Some(patch) = self.patches.get(name) {
                    let hr = patch.coefficient(driving[i])?;
                    let ha = *htc.get(name).ok_or_else(|| bad("missing endpoint air coefficient"))?;
                    let total = finite(ha + hr, "endpoint combined coefficient")?;
                    combined_h.insert(name.to_string(), total);
                    combined_ref[i] = finite((ha / total) * references[i]
                        + (hr / total) * patch.ambient_k, "endpoint combined reference")?;
                    applied_h[i] = hr;
                }
            }
            let boundary = request.boundary(names, &combined_ref, &combined_h)?;
            let rows = retain_rows.then(|| SolvedRows { htc: combined_h, references: combined_ref,
                driving: driving.clone(), radiative_htc: applied_h });
            let step = solve(&boundary)?;
            poll(cx)?;
            let mut states = Vec::with_capacity(names.len());
            let mut max_change = 0.0_f64;
            for (i, &name) in names.iter().enumerate() {
                let flux = step.robin_fluxes.iter().find(|f| f.region == name)
                    .ok_or_else(|| bad("radiative endpoint lacks an assembled Robin region"))?;
                let mean = flux.mean_wall_temperature_k;
                if !(mean.is_finite() && mean > 0.0) {
                    return Err(producer("radiative endpoint has a nonpositive wall temperature"));
                }
                if self.patches.contains_key(name) {
                    max_change = max_change.max(finite(mean - driving[i], "endpoint radiation update")?.abs());
                    driving[i] = finite((1.0 - self.relaxation) * driving[i] + self.relaxation * mean,
                        "endpoint relaxed radiation temperature")?;
                }
                let heat = finite(htc[name] * flux.area_m2 * (mean - references[i]), "endpoint convective heat")?;
                states.push(SolidRegionState { region: name.to_string(), area_m2: flux.area_m2,
                    mean_wall_temperature_k: mean, heat_rate_w: heat,
                    mean_reference_temperature_k: Some(references[i]) });
            }
            let heat = self.collect_endpoint_heat(request, cx, &states, &step)?;
            if max_change <= self.tolerance_k
                && heat.max_mismatch_w <= request.limits.heat / (self.patches.len() as f64 + 1.0)
            {
                poll(cx)?;
                return Ok((step, states, rows));
            }
        }
        Err(Failure { code: "cooling-network-radiation-budget", message: format!(
            "radiative time endpoint did not meet temperature and nonlinear watt gates within {} solid solves; previous physical history was not advanced",
            self.max_iterations) })
    }

    pub(in crate::network_command) fn endpoint_heat(
        &self, request: &Request, cx: &Cx<'_>, states: &[SolidRegionState], step: &StepSolution,
    ) -> Result<EndpointHeat> {
        let heat = self.collect_endpoint_heat(request, cx, states, step)?;
        if heat.max_mismatch_w > request.limits.heat / (self.patches.len() as f64 + 1.0) {
            return Err(producer("accepted radiative endpoint has an unclosed nonlinear patch flux"));
        }
        Ok(heat)
    }

    fn collect_endpoint_heat(
        &self, request: &Request, cx: &Cx<'_>, states: &[SolidRegionState], step: &StepSolution,
    ) -> Result<EndpointHeat> {
        if let Some(enclosure) = &self.enclosure {
            // The enclosure producer checks its stricter per-patch budget and
            // closed heat sum. It must not be treated as an empty reservoir list.
            let heat = enclosure.endpoint_heat(request, cx, states, step)?;
            return Ok(EndpointHeat { outward_w: heat.outward_w, applied_w: heat.applied_w,
                max_mismatch_w: heat.max_mismatch_w, rows: Vec::new(),
                enclosure_report: Some(heat.report) });
        }
        if states.len() != step.robin_fluxes.len() || states.len() != request.surfaces.len() {
            return Err(bad("endpoint radiation/Robin decomposition arity mismatch"));
        }
        let mut rows = Vec::with_capacity(self.patches.len());
        let mut outward_w = 0.0;
        let mut applied_w = 0.0;
        let mut max_mismatch_w = 0.0_f64;
        let mut assembled = 0.0;
        for state in states {
            poll(cx)?;
            let flux = step.robin_fluxes.iter().find(|f| f.region == state.region)
                .ok_or_else(|| bad("endpoint radiation/Robin identity mismatch"))?;
            assembled = finite(assembled + flux.heat_rate_w, "endpoint assembled Robin sum")?;
            let applied = finite(flux.heat_rate_w - state.heat_rate_w, "endpoint radiative heat split")?;
            if let Some(patch) = self.patches.get(&state.region) {
                let h = patch.coefficient(state.mean_wall_temperature_k)?;
                let q = finite(h * state.area_m2 * (state.mean_wall_temperature_k - patch.ambient_k),
                    "endpoint nonlinear radiative heat")?;
                max_mismatch_w = max_mismatch_w.max(finite(q - applied, "endpoint radiation mismatch")?.abs());
                outward_w = finite(outward_w + q, "endpoint nonlinear radiation sum")?;
                applied_w = finite(applied_w + applied, "endpoint applied radiation sum")?;
                rows.push(PatchHeat { surface: state.region.clone(), mean_k: state.mean_wall_temperature_k,
                    secant_h: h, applied_w: applied, nonlinear_w: q });
            } else if applied.abs() > request.limits.heat {
                return Err(producer("nonradiating endpoint region has unexplained boundary heat"));
            }
        }
        if rows.len() != self.patches.len()
            || finite(assembled - step.robin_out_w, "endpoint Robin decomposition")?.abs() > request.limits.heat
        {
            return Err(producer("radiative endpoint boundary partition does not close"));
        }
        Ok(EndpointHeat { outward_w, applied_w, max_mismatch_w, rows, enclosure_report: None })
    }

    pub(in crate::network_command) fn endpoint_report(&self, heat: &EndpointHeat) -> Result<String> {
        if let Some(report) = &heat.enclosure_report { return Ok(report.clone()); }
        let rows = heat.rows.iter().map(|heat| {
            let patch = &self.patches[&heat.surface];
            Ok(format!("{{\"surface\":{},\"emissivity\":{},\"ambient_temperature_k\":{},\"mean_temperature_k\":{},\"secant_htc_w_m2_k\":{},\"applied_heat_w\":{},\"nonlinear_heat_w\":{},\"source\":{}}}",
                quote(&patch.surface), num(patch.emissivity.value())?, num(patch.ambient_k)?, num(heat.mean_k)?,
                num(heat.secant_h)?, num(heat.applied_w)?, num(heat.nonlinear_w)?, quote(&patch.source)))
        }).collect::<Result<Vec<_>>>()?.join(",");
        Ok(format!("{{\"model\":\"mean-patch-gray-to-isothermal-surroundings\",\"temporal_scope\":\"final-accepted-endpoint\",\"radiative_out_w\":{},\"applied_radiative_out_w\":{},\"max_nonlinear_mismatch_w\":{},\"patches\":[{}],\"scope\":\"implicit new-temperature radiation alongside convection; mean-patch Robin closure, not pointwise T(x)^4 integration; no fluid radiation absorption, enclosure reflections, moving surroundings or continuous-time peak bound; trajectory derivatives, when requested, are reported separately\"}}",
            num(heat.outward_w)?, num(heat.applied_w)?, num(heat.max_mismatch_w)?, rows))
    }
}
