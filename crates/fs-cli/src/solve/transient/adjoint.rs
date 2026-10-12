//! Opt-in reverse differentiation of the accepted native fine time grid.
//! Reconstruction must reproduce the retained endpoint before any derivative
//! is returned. Forward temperatures, energy and decisions remain unchanged.

use std::time::Instant;
use fs_exec::Cx;
use fs_project::{ProjectSpec, ThermalBoundaryCondition};
use super::{SolveRefusal, bad, poll};

mod reverse;
pub(super) use reverse::finish;
mod enthalpy;
pub(super) use enthalpy::{
    EnthalpyTape, finish as finish_enthalpy, storage_bytes as enthalpy_storage_bytes,
    radiation_storage_bytes as enthalpy_radiation_storage_bytes,
    convection_storage_bytes as enthalpy_convection_storage_bytes,
};

#[cfg(test)]
mod tests;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Objective { Final, SampledPeak }

impl Objective {
    fn name(self) -> &'static str {
        match self {
            Self::Final => "temperature-final-adjoint",
            Self::SampledPeak => "temperature-sampled-peak-adjoint",
        }
    }
}

pub(super) fn requested(spec: &ProjectSpec) -> Result<Option<Objective>, SolveRefusal> {
    let mut selected = None;
    for output in spec.outputs.iter().flatten() {
        let objective = match output.name.as_str() {
            "temperature-final-adjoint" => Objective::Final,
            "temperature-sampled-peak-adjoint" => Objective::SampledPeak,
            _ => continue,
        };
        if selected.is_some() || output.kind != "report" || output.region.is_some() {
            return Err(bad("select one transient adjoint report without an output region; it uses the existing temperature-max requirement region"));
        }
        selected = Some(objective);
    }
    if selected.is_some() {
        let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
            .filter(|setup| setup.transient.is_some())
            .ok_or_else(|| bad("a transient adjoint report requires finite-time conduction"))?;
        let enthalpy = setup.transient.as_ref().is_some_and(|time| time.enthalpy.is_some());
        if !enthalpy && setup.boundaries.iter().any(|b| matches!(
            b.condition, ThermalBoundaryCondition::NaturalConvection { .. }
                | ThermalBoundaryCondition::AirflowConvection { .. })) {
            return Err(bad("native capacity-storage adjoints require prescribed convection laws; use total-enthalpy storage for complete natural and airflow feedback"));
        }
        if enthalpy && setup.boundaries.iter().any(|b| matches!(
            &b.condition, ThermalBoundaryCondition::NaturalConvection { correlation, .. }
                if correlation != crate::solve::natural::DIFFERENTIABLE_CARD)) {
            return Err(bad("native natural-convection enthalpy adjoints require the differentiated Churchill-Chu vertical-plate card"));
        }
        if setup.radiation.is_some() && !enthalpy {
            return Err(bad("native capacity-storage adjoints do not yet differentiate radiation; total-enthalpy storage has its own complete radiative endpoint owner"));
        }
        if setup.regions.len() > 64 || setup.boundaries.len() > 64 {
            return Err(bad("native transient adjoints admit at most 64 volume regions and boundary laws"));
        }
    }
    Ok(selected)
}

pub(super) fn is_output(name: &str) -> bool {
    matches!(name, "temperature-final-adjoint" | "temperature-sampled-peak-adjoint")
}

fn finite(value: f64) -> Result<f64, SolveRefusal> {
    if value.is_finite() { Ok(value) } else { Err(bad("nonfinite native transient adjoint arithmetic")) }
}

/// Additional retained-state, source-basis, reverse-scratch and output allowance.
/// The forward driver adds this to its existing estimate; neither is an RSS proof.
pub(super) fn storage_bytes(n: usize, steps: usize, regions: usize) -> Result<usize, SolveRefusal> {
    let states = steps.checked_add(1).ok_or_else(|| bad("adjoint state count overflows"))?;
    states.checked_mul(n).and_then(|v| v.checked_mul(8))
        .and_then(|b| states.checked_mul(std::mem::size_of::<Vec<f64>>())
            .and_then(|v| b.checked_add(v)))
        .and_then(|b| regions.checked_add(24).and_then(|r| r.checked_mul(n))
            .and_then(|v| v.checked_mul(8)).and_then(|v| b.checked_add(v)))
        .and_then(|b| regions.checked_mul(192).and_then(|r| r.checked_add(1024))
            .and_then(|r| r.checked_mul(steps)).and_then(|r| b.checked_add(r)))
        .ok_or_else(|| bad("native transient adjoint storage estimate overflows"))
}

pub(super) struct Tape {
    states: Vec<Vec<f64>>,
    n: usize,
    limit: usize,
}

impl Tape {
    pub(super) fn new(cx: &Cx<'_>, initial: &[f64], steps: usize,
        deadline: Option<(Instant, f64)>) -> Result<Self, SolveRefusal> {
        let limit = steps.checked_add(1).ok_or_else(|| bad("adjoint state count overflows"))?;
        let mut states = Vec::new();
        states.try_reserve_exact(limit).map_err(|_| bad("adjoint checkpoint allocation refused"))?;
        let mut tape = Self { states, n: initial.len(), limit };
        tape.record(cx, initial, deadline)?;
        Ok(tape)
    }

    pub(super) fn record(&mut self, cx: &Cx<'_>, values: &[f64],
        deadline: Option<(Instant, f64)>) -> Result<(), SolveRefusal> {
        poll(cx, deadline)?;
        if values.len() != self.n || self.states.len() == self.limit {
            return Err(bad("adjoint checkpoint shape or state-count allowance differs"));
        }
        let mut state = Vec::new();
        state.try_reserve_exact(self.n).map_err(|_| bad("adjoint field allocation refused"))?;
        for (i, &value) in values.iter().enumerate() {
            if i % 256 == 0 { poll(cx, deadline)?; }
            state.push(finite(value)?);
        }
        poll(cx, deadline)?;
        self.states.push(state);
        Ok(())
    }
}

#[derive(Debug)]
struct Selected { state: usize, vertex: usize, value: f64, ties: usize }

fn select(cx: &Cx<'_>, tape: &Tape, objective: Objective, vertices: &[usize],
    deadline: Option<(Instant, f64)>) -> Result<Selected, SolveRefusal> {
    if tape.states.len() != tape.limit || vertices.is_empty() {
        return Err(bad("transient adjoint needs a complete fine trajectory and a nonempty goal region"));
    }
    let first = if objective == Objective::Final { tape.states.len() - 1 } else { 0 };
    let mut selected: Option<Selected> = None;
    for index in first..tape.states.len() {
        for (slot, &vertex) in vertices.iter().enumerate() {
            if slot % 256 == 0 { poll(cx, deadline)?; }
            let value = *tape.states[index].get(vertex)
                .ok_or_else(|| bad("goal vertex is outside the retained transient field"))?;
            match &mut selected {
                Some(old) if value == old.value => old.ties += 1,
                Some(old) if value < old.value => {},
                _ => selected = Some(Selected { state: index, vertex, value, ties: 1 }),
            }
        }
    }
    selected.ok_or_else(|| bad("transient adjoint has no selected objective"))
}
