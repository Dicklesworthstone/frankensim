//! Conservative upwind enthalpy transport on a prescribed closed mass-flow graph.
//!
//! This is an observer of the original enthalpy timeline, not another thermal
//! integrator. Cell mass, volume and material identity stay fixed. Each flowing
//! component must use one immutable caloric law and balance mass at every cell.
//! Species mixing, density/volume changes, pressure work and phase-dependent
//! mobility need their own models; no velocity or melt motion is inferred here.

use std::collections::{BTreeMap, BTreeSet};
use super::{CellId, DepositionError, DepositionObserver, DepositionSimulation, Event, HeatTransfer};

/// One physical connection, independent of conductive-contact identity/area.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MassFlow {
    /// Stable connection identity. Parallel connections need distinct IDs.
    pub id: u64,
    /// First active thermal control volume.
    pub a: CellId,
    /// Second active thermal control volume.
    pub b: CellId,
    /// Signed mass rate [kg/s], positive from a to b.
    pub kg_s: f64,
}

/// Numerical observations for the last successful advance, not new heat input.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct AdvectionAudit {
    /// Sum of absolute mass transported across all connections [kg]. A mass
    /// parcel may cross several connections; this is NOT added/stored mass.
    pub transported_mass_kg: f64,
    /// Sum of absolute connection enthalpy transfers [J]. This depends on the
    /// common enthalpy datum; it is NOT irreversible dissipation or external heat.
    pub absolute_enthalpy_transfer_j: f64,
    /// Largest outgoing mass / stored cell mass on an accepted substep.
    pub maximum_courant: f64,
    /// Actual accepted thermal substeps carrying the prescribed flow graph.
    pub substeps: usize,
}

#[derive(Clone)]
struct Pending {
    time_s: f64,
    decay_per_s: f64,
    mass_kg_s: f64,
    absolute_power_w: f64,
}

/// Add prescribed flow to ANY existing thermal source/history observer.
///
/// The wrapped observer still owns its sources, conductances, events and work
/// accounts. Advection contributes pairwise INTERNAL watts and a combined CFL
/// restriction to that very same step. Material T(h), including sensible-heat
/// curvature and latent plateaux, is always evaluated by the original timeline.
/// Recreate this wrapper at an accepted endpoint to change the prescribed flow.
#[derive(Clone)]
pub struct Advected<O> {
    inner: O,
    flows: Vec<MassFlow>,
    pending: Option<Pending>,
    audit: AdvectionAudit,
}

fn finite(x: f64) -> Result<f64, DepositionError> {
    if x.is_finite() { Ok(x) } else { Err(DepositionError::Arithmetic) }
}
fn invalid(s: &'static str) -> DepositionError { DepositionError::InvalidInput(s) }

impl<O> Advected<O> {
    /// Admit bounded, uniquely named connections. Sort by stable ID so input
    /// order does not choose floating-point reduction order. Active-cell, common
    /// material and mass-balance admission happens at each actual thermal step.
    pub fn new(inner: O, mut flows: Vec<MassFlow>, max_flows: usize) -> Result<Self, DepositionError> {
        if flows.len() > max_flows { return Err(DepositionError::Budget("advective connections")); }
        let mut ids = BTreeSet::new();
        for flow in &flows {
            if flow.a == flow.b || !flow.kg_s.is_finite() || !ids.insert(flow.id) {
                return Err(invalid("invalid or duplicate advective connection"));
            }
        }
        flows.sort_by_key(|f| f.id);
        Ok(Self { inner, flows, pending: None, audit: AdvectionAudit::default() })
    }
    /// Original source/history owner, including any accepted material receipts.
    pub fn inner(&self) -> &O { &self.inner }
    /// Consume the wrapper without erasing the original source/history type.
    pub fn into_inner(self) -> O { self.inner }
    /// Immutable prescribed connection data in stable-ID order.
    pub fn flows(&self) -> &[MassFlow] { &self.flows }
    /// Observations reset by begin(), on the tentative observer only.
    pub fn audit(&self) -> AdvectionAudit { self.audit }
}

impl<O: DepositionObserver> DepositionObserver for Advected<O> {
    type Error = O::Error;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.inner.begin()?;
        self.pending = None;
        self.audit = AdvectionAudit::default();
        Ok(())
    }
    fn sources(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner.sources(thermal, check)
    }
    fn conductances(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.inner.conductances(thermal, check)
    }
    fn heat_transfers(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<Vec<HeatTransfer>, Self::Error> {
        check()?;
        let mut transfers = self.inner.heat_transfers(thermal, check)?;
        // Separate incoming/outgoing sums avoid cancellation in the admission
        // scale. Small roundoff is admitted, not an adjustable mass-leak budget.
        let mut rates: BTreeMap<CellId, [f64; 2]> = BTreeMap::new();
        let mut mass_kg_s = 0.0;
        let mut absolute_power_w = 0.0;
        for (index, flow) in self.flows.iter().enumerate() {
            if index % 256 == 0 { check()?; }
            let a = thermal.cells().get(&flow.a).ok_or(invalid("inactive advective cell"))?;
            let b = thermal.cells().get(&flow.b).ok_or(invalid("inactive advective cell"))?;
            if flow.kg_s == 0.0 { continue; }
            if a.material_id != b.material_id {
                return Err(invalid("advection across different material IDs requires a composition model").into());
            }
            let (up, down) = if flow.kg_s > 0.0 { (a, b) } else { (b, a) };
            let rate = flow.kg_s.abs();
            let power_w = finite(rate * up.specific_h_j_kg)?;
            let out = rates.entry(up.id).or_default();
            out[1] = finite(out[1] + rate)?;
            let incoming = rates.entry(down.id).or_default();
            incoming[0] = finite(incoming[0] + rate)?;
            mass_kg_s = finite(mass_kg_s + rate)?;
            absolute_power_w = finite(absolute_power_w + power_w.abs())?;
            transfers.push(HeatTransfer { a: up.id, b: down.id, power_w });
        }
        let mut decay_per_s: f64 = 0.0;
        for (index, (&id, &[incoming, outgoing])) in rates.iter().enumerate() {
            if index % 256 == 0 { check()?; }
            let scale = incoming.max(outgoing);
            if (incoming-outgoing).abs() > 64.0*f64::EPSILON*scale {
                return Err(invalid("fixed-mass advection requires balanced incoming and outgoing mass at every cell").into());
            }
            let decay = finite(outgoing / thermal.cells()[&id].mass_kg)?;
            if outgoing != 0.0 && decay == 0.0 {
                return Err(DepositionError::Arithmetic.into());
            }
            decay_per_s = decay_per_s.max(decay);
        }
        check()?;
        self.pending = Some(Pending { time_s: thermal.time_s(), decay_per_s,
            mass_kg_s, absolute_power_w });
        Ok(transfers)
    }
    fn limit_timestep_s(&mut self, thermal: &DepositionSimulation<'_>, proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<f64, Self::Error> {
        let dt = self.inner.limit_timestep_s(thermal, proposed_s, check)?;
        if !dt.is_finite() || dt <= 0.0 || dt > proposed_s {
            return Err(invalid("wrapped observer returned an invalid timestep bound").into());
        }
        let pending = self.pending.as_ref().ok_or(invalid("advection has no prepared source epoch"))?;
        if pending.time_s != thermal.time_s() {
            return Err(invalid("advection source epoch does not match thermal state").into());
        }
        if pending.decay_per_s == 0.0 { return Ok(dt); }
        // dt <= 0.8/a from the existing conductive kernel; b is the largest
        // advective outgoing mass/mass. The harmonic limit obeys dt <= 0.8/(a+b),
        // unlike independently imposing two CFL limits (which can overshoot).
        // Stable reciprocal evaluation avoids overflowing dt*b or 1/tiny_dt.
        let advective_dt = 0.8 / pending.decay_per_s;
        let result = if dt <= advective_dt {
            dt / (1.0 + dt/advective_dt)
        } else {
            advective_dt / (1.0 + advective_dt/dt)
        };
        if !result.is_finite() || result <= 0.0 {
            return Err(DepositionError::Arithmetic.into());
        }
        Ok(result)
    }
    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.inner.event(event, thermal)
    }
    fn step(&mut self, start_s: f64, dt: f64, old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>, check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<(), Self::Error> {
        check()?;
        let pending = self.pending.take().ok_or(invalid("advection step has no prepared flow"))?;
        if pending.time_s != start_s || !dt.is_finite() || dt <= 0.0 {
            return Err(invalid("advection step/source epoch mismatch").into());
        }
        self.audit.transported_mass_kg = finite(self.audit.transported_mass_kg + finite(dt*pending.mass_kg_s)?)?;
        self.audit.absolute_enthalpy_transfer_j = finite(self.audit.absolute_enthalpy_transfer_j
            + finite(dt*pending.absolute_power_w)?)?;
        self.audit.maximum_courant = self.audit.maximum_courant.max(finite(dt*pending.decay_per_s)?);
        self.audit.substeps = self.audit.substeps.checked_add(1).ok_or(DepositionError::Arithmetic)?;
        self.inner.step(start_s, dt, old, thermal, check)
    }
}

#[cfg(test)]
#[path = "advection/tests.rs"]
mod tests;
