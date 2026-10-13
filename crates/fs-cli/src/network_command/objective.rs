//! Mean and exact discrete P1 temperature-maximum objectives. Ordinary maxima
//! retain their active-vertex derivative. A design-only smooth upper envelope
//! is labeled separately and retains the actual maximum for feasibility checks.

use fs_airflow::graph::thermal::coupled_transport::sensitivity::CoupledObjective;
use super::*;

mod envelope;
use envelope::Envelope;

#[derive(Debug)]
pub(super) enum Objective {
    MeanWall(String),
    Maximum { kind: &'static str, region: Option<String>, vertices: Vec<usize>, smoothing_k: Option<f64> },
}

#[derive(Debug)]
pub(super) struct ObjectiveState {
    pub value: f64,
    wall_index: Option<usize>,
    pub vertex: Option<usize>,
    pub ties: usize,
    pub separation_k: Option<f64>,
    envelope: Option<Envelope>,
}

impl Objective {
    pub fn parse(value: &J, surfaces: &[Surface], mesh: &ConductionMesh) -> Result<Self> {
        object(value, &["mean_wall_region", "max_wall_region", "max_solid_temperature", "max_vertices", "gradient"], "objective")?;
        let keys = ["mean_wall_region", "max_wall_region", "max_solid_temperature", "max_vertices"];
        let selected: Vec<_> = keys.iter().copied().filter(|key| value.get(key).is_some()).collect();
        if selected.len() != 1 { return Err(bad("objective requires exactly one mean or maximum selector")); }
        let key = selected[0];
        let mut region = None;
        let vertices = match key {
            "mean_wall_region" | "max_wall_region" => {
                let name = string(get(value, key)?, key)?;
                let surface = surfaces.iter().find(|s| s.name == name)
                    .ok_or_else(|| bad("objective names an unknown solid surface"))?;
                if key == "mean_wall_region" { return Ok(Self::MeanWall(name)); }
                region = Some(name);
                surface.faces.iter().flatten().map(|&v| v as usize).collect::<BTreeSet<_>>()
            }
            "max_solid_temperature" => {
                if !boolean(get(value, key)?, key)? { return Err(bad("max_solid_temperature must be true")); }
                (0..mesh.vertex_count()).collect()
            }
            _ => {
                let mut vertices = BTreeSet::new();
                for entry in array(get(value, key)?, key, mesh.vertex_count())? {
                    let vertex = integer(entry, "maximum vertex", mesh.vertex_count() - 1)?;
                    if !vertices.insert(vertex) { return Err(bad("max_vertices repeats a vertex")); }
                }
                vertices
            }
        };
        if vertices.is_empty() { return Err(bad("temperature maximum requires a nonempty vertex set")); }
        Ok(Self::Maximum { kind: key, region, vertices: vertices.into_iter().collect(), smoothing_k: None })
    }

    pub fn is_mean(&self) -> bool { matches!(self, Self::MeanWall(_)) }
    pub fn kind(&self) -> &'static str {
        match self { Self::MeanWall(_) => "mean_wall_region", Self::Maximum { kind, .. } => kind }
    }
    pub fn region(&self) -> Option<&str> {
        match self { Self::MeanWall(name) => Some(name), Self::Maximum { region, .. } => region.as_deref() }
    }

    /// Only a design policy can enable smoothing; ordinary model objectives
    /// continue to return the exact maximum and selected-active-vertex slope.
    #[cfg(any(feature = "cooling-design", test))]
    pub fn enable_upper_envelope(&mut self, smoothing_k: f64) -> Result<()> {
        if !smoothing_k.is_finite() || smoothing_k <= 0.0 {
            return Err(bad("maximum_smoothing_k must be finite and positive"));
        }
        match self {
            Self::Maximum { smoothing_k: slot, .. } => { *slot = Some(smoothing_k); Ok(()) }
            Self::MeanWall(_) => Err(bad("maximum_smoothing_k requires a maximum-temperature selector")),
        }
    }

    pub fn evaluate(&self, cx: &Cx<'_>, temperatures: &[f64], walls: &[SolidRegionState]) -> Result<ObjectiveState> {
        poll(cx)?;
        match self {
            Self::MeanWall(name) => {
                let index = walls.iter().position(|s| s.region == *name).ok_or_else(|| bad("missing objective surface"))?;
                let value = walls[index].mean_wall_temperature_k;
                if !value.is_finite() { return Err(producer("nonfinite mean-wall objective")); }
                Ok(ObjectiveState { value, wall_index: Some(index), vertex: None, ties: 0, separation_k: None, envelope: None })
            }
            Self::Maximum { vertices, smoothing_k, .. } => {
                let mut largest = f64::NEG_INFINITY;
                let mut second = f64::NEG_INFINITY;
                let mut winner = None;
                let mut ties = 0;
                for &vertex in vertices {
                    poll(cx)?;
                    let value = *temperatures.get(vertex).ok_or_else(|| bad("maximum vertex is outside the solved field"))?;
                    if !value.is_finite() { return Err(producer("nonfinite temperature in maximum objective")); }
                    if value > largest {
                        second = largest; largest = value; winner = Some(vertex); ties = 1;
                    } else if value == largest {
                        ties += 1;
                        winner = Some(winner.map_or(vertex, |prior: usize| prior.min(vertex)));
                    } else { second = second.max(value); }
                }
                if winner.is_none() { return Err(bad("empty maximum vertex selection")); }
                let separation = if ties > 1 { Some(0.0) }
                    else if second.is_finite() { Some(largest - second) } else { None };
                let envelope = smoothing_k.map(|tau| Envelope::new(cx, temperatures, vertices, largest, tau)).transpose()?;
                let value = envelope.as_ref().map_or(largest, |e| e.upper_k);
                Ok(ObjectiveState { value, wall_index: None, vertex: winner, ties,
                    separation_k: separation, envelope })
            }
        }
    }

    pub fn render(&self, state: &ObjectiveState, mesh: &ConductionMesh) -> Result<String> {
        let vertex = state.vertex.map_or_else(|| "null".into(), |v| v.to_string());
        let position = state.vertex.map(|v| numbers(&mesh.positions()[v])).transpose()?.unwrap_or_else(|| "null".into());
        let derivative = if state.envelope.is_some() {
            "smooth log-sum-exp upper-envelope derivative on the declared maximum selection; not the derivative of the raw maximum"
        } else if self.is_mean() { "coupled area-mean derivative" }
            else { "selected-active-vertex derivative; exact ties use the lowest vertex id, not a unique max gradient; near-tie stability is not certified" };
        let kind = if state.envelope.is_some() { "smooth_upper_envelope" } else { self.kind() };
        let extra = match &state.envelope {
            None => String::new(),
            Some(e) => format!(",\"maximum_selection\":{},\"actual_maximum_k\":{},\"maximum_smoothing_k\":{},\"analytic_bias_bound_k\":{},\"envelope_scope\":\"max <= upper <= max + smoothing*ln(selected vertex count) in exact arithmetic; floating-point evaluation, not an interval or continuum certificate\"",
                quote(self.kind()), num(e.maximum_k)?, num(e.smoothing_k)?, num(e.bias_bound_k)?),
        };
        Ok(format!("{{\"kind\":{},\"region\":{},\"value_k\":{},\"active_vertex\":{vertex},\"position_m\":{position},\"exact_tie_count\":{},\"separation_k\":{},\"derivative_semantics\":{}{extra}}}",
            quote(kind), self.region().map_or_else(|| "null".into(), quote), num(state.value)?,
            state.ties, optional(state.separation_k)?, quote(derivative)))
    }
}

impl ObjectiveState {
    /// The mean or actual selected maximum, never its smooth upper envelope.
    #[cfg(any(feature = "cooling-design", test))]
    pub fn physical_value(&self) -> f64 {
        self.envelope.as_ref().map_or(self.value, |e| e.maximum_k)
    }

    pub fn seed(&self, weights: &mut CoupledObjective) {
        if let Some(envelope) = &self.envelope {
            for &(vertex, weight) in &envelope.weights { weights.nodal_temperatures[vertex] = weight; }
        } else {
            if let Some(index) = self.wall_index { weights.wall_temperatures[index] = 1.0; }
            if let Some(vertex) = self.vertex { weights.nodal_temperatures[vertex] = 1.0; }
        }
    }
}

#[cfg(test)]
mod tests;
