//! Uniform prescribed baths share the complete accepted physical transpose.
//! Each constrained node contributes once, including a selected fixed-node
//! objective. Declared baths sharing vertices have no independent coordinate.

use std::collections::BTreeMap;
use std::time::Instant;

use fs_conduction::{ConductionProblem, ScalarField, ThermalBc};
use fs_exec::Cx;
use fs_project::{ProjectSpec, ThermalBoundaryCondition};

use super::{SolveRefusal, bad, finite, json_string, number, poll};

struct Bath {
    target: String,
    temperature_k: f64,
    vertices: Vec<usize>,
    mapped: bool,
    shared: bool,
    derivative: f64,
}

pub(in crate::solve::transient::adjoint) struct Controls {
    baths: Vec<Bath>,
    vertices: usize,
    steps: usize,
}

impl Controls {
    pub(in crate::solve::transient::adjoint) fn new(
        cx: &Cx<'_>,
        spec: &ProjectSpec,
        problem: ConductionProblem<'_>,
        deadline: Option<(Instant, f64)>,
    ) -> Result<Self, SolveRefusal> {
        poll(cx, deadline)?;
        let setup = spec.cooling.as_ref().and_then(|cooling| cooling.conduction.as_ref())
            .ok_or_else(|| bad("prescribed bath controls lost their conduction declaration"))?;
        let mut baths = Vec::new();
        baths.try_reserve_exact(setup.boundaries.len())
            .map_err(|_| bad("prescribed bath-control metadata allocation refused"))?;
        let mut regions = BTreeMap::new();
        for declared in &setup.boundaries {
            poll(cx, deadline)?;
            let ThermalBoundaryCondition::FixedTemperature { temperature }
                = &declared.condition else { continue; };
            let region = problem.boundary.region_names().iter()
                .position(|name| name == &declared.target)
                .ok_or_else(|| bad("prescribed bath control has no physical boundary identity"))?;
            let ThermalBc::Dirichlet { temperature: ScalarField::Uniform(actual) }
                = &problem.boundary.conditions()[region]
                else { return Err(bad("prescribed bath control changed its uniform boundary law")); };
            if *actual != temperature.value || !actual.is_finite() || *actual <= 0.0 {
                return Err(bad("prescribed bath control differs from its physical temperature"));
            }
            if regions.insert(region, baths.len()).is_some() {
                return Err(bad("prescribed bath controls have duplicate declaration ownership"));
            }
            baths.push(Bath {
                target: declared.target.clone(),
                temperature_k: *actual,
                vertices: Vec::new(),
                mapped: false,
                shared: false,
                derivative: 0.0,
            });
        }
        let n = problem.mesh.vertex_count();
        if baths.is_empty() {
            if !problem.boundary.dirichlet().is_empty() {
                return Err(bad("a prescribed boundary has no declared control"));
            }
            return Ok(Self { baths, vertices: n, steps: 0 });
        }
        let mut owners = Vec::<Option<usize>>::new();
        owners.try_reserve_exact(n)
            .map_err(|_| bad("prescribed bath-control owner allocation refused"))?;
        owners.resize(n, None);
        for (slot, face) in problem.mesh.boundary().iter().enumerate() {
            if slot % 256 == 0 { poll(cx, deadline)?; }
            let Some(region) = problem.boundary.region_for(slot) else { continue; };
            if !matches!(problem.boundary.conditions()[region], ThermalBc::Dirichlet { .. }) {
                continue;
            }
            let control = *regions.get(&region)
                .ok_or_else(|| bad("a prescribed boundary has no declared control"))?;
            baths[control].mapped = true;
            for &vertex in &face.vertices {
                let owner = owners.get_mut(vertex as usize)
                    .ok_or_else(|| bad("prescribed bath control lost its vertex layout"))?;
                if let Some(previous) = *owner {
                    if previous != control {
                        baths[previous].shared = true;
                        baths[control].shared = true;
                    }
                } else {
                    *owner = Some(control);
                }
            }
        }
        // Walking the nodal owner map deduplicates all triangles of one bath.
        // Shared declarations keep identities/reasons but no false derivative.
        for (vertex, owner) in owners.into_iter().enumerate() {
            if vertex % 256 == 0 { poll(cx, deadline)?; }
            if let Some(control) = owner {
                if !baths[control].shared {
                    baths[control].vertices.try_reserve(1)
                        .map_err(|_| bad("prescribed bath trace allocation refused"))?;
                    baths[control].vertices.push(vertex);
                }
            }
        }
        if baths.iter().any(|bath| !bath.mapped || (!bath.shared && bath.vertices.is_empty())) {
            return Err(bad("a declared bath has no retained trace vertices"));
        }
        poll(cx, deadline)?;
        Ok(Self { baths, vertices: n, steps: 0 })
    }

    pub(in crate::solve::transient::adjoint) fn has_controls(&self) -> bool {
        !self.baths.is_empty()
    }

    pub(in crate::solve::transient::adjoint) fn add(
        &mut self,
        cx: &Cx<'_>,
        nodal_bars: &[f64],
        deadline: Option<(Instant, f64)>,
    ) -> Result<(), SolveRefusal> {
        poll(cx, deadline)?;
        if nodal_bars.len() != self.vertices {
            return Err(bad("prescribed bath derivative lost its nodal layout"));
        }
        for bath in &mut self.baths {
            for (index, &vertex) in bath.vertices.iter().enumerate() {
                if index % 256 == 0 { poll(cx, deadline)?; }
                bath.derivative = finite(bath.derivative + nodal_bars[vertex])?;
            }
        }
        self.steps = self.steps.checked_add(1)
            .ok_or_else(|| bad("prescribed bath-control step count overflows"))?;
        poll(cx, deadline)
    }

    pub(super) fn fragment(&self, expected_steps: usize) -> Result<String, SolveRefusal> {
        self.fragment_with_scope(expected_steps, "Each retained constrained vertex contributes once, including the direct selected-node objective and complete conductive, contact and reservoir feedback. Initial enthalpy is held fixed; an initial sampled peak has zero future bath sensitivity. Plateau ambiguity, chart endpoints and slope-changing branches refuse classical derivatives. Shared prescribed vertices do not admit independent boundary coordinates.")
    }

    pub(in crate::solve::transient::adjoint) fn fragment_for_capacity(
        &self, expected_steps: usize,
    ) -> Result<String, SolveRefusal> {
        self.fragment_with_scope(expected_steps, "Each retained constrained vertex contributes once, including the direct selected-node objective and the complete conductive, contact, natural-convection, transported-air and radiation lift. The same checked watt-coordinate dual is used without another timestep factor or solve. Initial temperatures and the row-sum-lumped capacity are held fixed; an initial sampled peak has zero future bath sensitivity. Shared prescribed vertices do not admit independent boundary coordinates.")
    }

    fn fragment_with_scope(&self, expected_steps: usize, scope: &str) -> Result<String, SolveRefusal> {
        if self.steps != expected_steps {
            return Err(bad("prescribed bath controls do not cover the selected reverse prefix"));
        }
        if self.baths.is_empty() { return Ok(String::new()); }
        let rows = self.baths.iter().map(|bath| {
            let state = if bath.shared {
                format!(
                    "\"status\":\"unsupported\",\"reason\":{}",
                    json_string("this trace shares prescribed vertices with another fixed-temperature declaration; independent variation would impose conflicting values"),
                )
            } else {
                format!(
                    "\"status\":\"available\",\"mapped_vertices\":{},\"dtemperature_dprescribed_temperature_k_per_k\":{}",
                    bath.vertices.len(), number(bath.derivative)?,
                )
            };
            Ok(format!(
                "{{\"target\":{},\"temperature_k\":{},{state}}}",
                json_string(&bath.target), number(bath.temperature_k)?,
            ))
        }).collect::<Result<Vec<_>, SolveRefusal>>()?;
        Ok(format!(
            ",\"prescribed_temperature\":{{\"control_steps\":{},\"parameter_scope\":\"independent-uniform-bath-temperatures-held-throughout-selected-trajectory\",\"scope\":{},\"boundaries\":[{}]}}",
            self.steps, json_string(scope), rows.join(","),
        ))
    }
}
