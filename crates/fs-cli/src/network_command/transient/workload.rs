//! Independent per-component power histories. Footprints stay fixed while watts
//! change. Reuse the production PowerMap projection, including zero-base loads;
//! never rescale a previously combined nodal source to guess its contributors.

use super::*;
use fs_conduction::{ComponentPower, PowerMap, PowerUncertainty};

#[derive(Debug, Clone)]
pub(super) enum Workload {
    Scale(f64),
    Components(BTreeMap<String, f64>),
}

pub(super) struct PreparedLoad {
    pub source: ScalarField,
    pub expected_power_w: Option<f64>,
}

impl Workload {
    pub(super) fn parse(interval: &J) -> Result<Self> {
        match (interval.get("power_scale"), interval.get("component_powers_w")) {
            (Some(scale), None) => {
                let scale = number(scale, "power_scale")?;
                if scale < 0.0 { return Err(bad("power_scale must be nonnegative")); }
                Ok(Self::Scale(scale))
            }
            (None, Some(value)) => {
                let rows = value.as_object().ok_or_else(|| bad("component_powers_w must be an object mapping names to watts"))?;
                if rows.is_empty() || rows.len() > 4096 {
                    return Err(bad("component_powers_w requires 1 to 4096 named powers"));
                }
                let mut powers = BTreeMap::new();
                for (name, watts) in rows {
                    let name = string(&J::Str(name.clone()), "workload component name")?;
                    let watts = number(watts, "component power in watts")?;
                    if watts < 0.0 { return Err(bad("component powers must be nonnegative")); }
                    if powers.insert(name, watts).is_some() { return Err(bad("duplicate workload component")); }
                }
                Ok(Self::Components(powers))
            }
            _ => Err(bad("each interval requires exactly one of power_scale and component_powers_w")),
        }
    }

    pub(super) fn validate(&self, request: &Request, cx: &Cx<'_>) -> Result<()> {
        poll(cx)?;
        match self {
            Self::Scale(scale) => {
                if !scale.is_finite() || *scale < 0.0 { return Err(bad("invalid workload multiplier")); }
            }
            Self::Components(powers) => {
                let base = request.solid_data.component_map.as_ref()
                    .ok_or_else(|| bad("named workloads require solid.component_power with fixed footprints"))?;
                if powers.len() != base.components().len() {
                    return Err(bad("component_powers_w must name every declared component exactly once, including zero-power components"));
                }
                for component in base.components() {
                    poll(cx)?;
                    let value = powers.get(component.name()).ok_or_else(|| bad(format!(
                        "workload is missing component {}; unknown replacement names are not admitted", component.name())))?;
                    if !value.is_finite() || *value < 0.0 { return Err(bad("invalid component workload")); }
                }
            }
        }
        Ok(())
    }

    pub(super) fn prepare(&self, request: &Request, cx: &Cx<'_>) -> Result<PreparedLoad> {
        self.validate(request, cx)?;
        match self {
            Self::Scale(scale) => Ok(PreparedLoad {
                source: scaled_source(request, *scale)?,
                expected_power_w: request.solid_data.power.as_ref()
                    .map(|audit| finite(audit.delivered_total_w() * *scale)).transpose()?,
            }),
            Self::Components(powers) => {
                let base = request.solid_data.component_map.as_ref()
                    .ok_or_else(|| bad("missing component footprints"))?;
                let mut components = Vec::with_capacity(base.components().len());
                let mut total = 0.0;
                // Base components are in production PowerMap name order. The
                // same summation order sets the explicit total for admission.
                for component in base.components() {
                    poll(cx)?;
                    let watts = *powers.get(component.name()).ok_or_else(|| bad("missing named workload"))?;
                    total = finite(total + watts)?;
                    components.push(ComponentPower::new(component.name(), watts, PowerUncertainty::Unstated,
                        component.vertices().to_vec()).map_err(producer)?);
                }
                let map = PowerMap::new(components, total).map_err(producer)?;
                let (source, audit) = map.volumetric_source(&request.mesh, 0.0).map_err(producer)?;
                source.validate("scheduled component source", request.mesh.vertex_count()).map_err(producer)?;
                let delivered = finite(audit.delivered_total_w())?;
                if (delivered - total).abs() > request.limits.heat {
                    return Err(producer("component workload projection missed the declared watt budget"));
                }
                poll(cx)?;
                Ok(PreparedLoad { source, expected_power_w: Some(delivered) })
            }
        }
    }

    pub(super) fn render(&self) -> Result<String> {
        match self {
            Self::Scale(scale) => Ok(format!("\"power_scale\":{},\"component_powers_w\":null", num(*scale)?)),
            Self::Components(powers) => {
                let fields = powers.iter().map(|(name, &watts)| Ok(format!("{}:{}", quote(name), num(watts)?)))
                    .collect::<Result<Vec<_>>>()?.join(",");
                Ok(format!("\"power_scale\":null,\"component_powers_w\":{{{fields}}}"))
            }
        }
    }
}

#[cfg(test)]
mod tests;
