//! Staggered mechanical separation/reclosure coupled to real thermal contact.
//!
//! A bonded original patch retains its supplied thermal conductance until the
//! cohesive law fails. A failed patch opens at a declared normal-gap threshold
//! and recontacts at a separate closure threshold. Closed failed patches use a
//! separately supplied conductance; recontact NEVER resets cohesive damage.
//! Partial damage is not interpreted as lost geometric area or scaled heat flow.
//!
//! Each accepted mechanical endpoint changes the NEXT thermal interval. This
//! does not localize a fracture time inside the preceding step, solve mechanical
//! equilibrium, move geometry, or model fracture heating. Fixed local frames,
//! small tangential travel and uniform original contact area are explicit limits.

use std::collections::BTreeMap;
use fs_exec::Cx;
use fs_material::cohesive::{CohesiveError, CohesivePoint};
use crate::transient::deposition::{CellId, ContactUpdate, DepositionError,
    DepositionObserver, DepositionSimulation, EnergyAudit, Event, StepConfig, ThermalLink};

/// Explicit reduced closure model, not a measured universal adhesion relation.
#[derive(Debug, Clone, Copy)]
pub struct ContactPolicy {
    /// A failed patch opens at or above this positive normal gap, m.
    pub open_gap_m: f64,
    /// An open failed patch recloses at or below this gap in [0, open_gap), m.
    pub reclose_gap_m: f64,
    /// Independent total W/K after debonding while still/re-again closed.
    pub reclosed_conductance_w_k: f64,
    /// Fixed-area/frame applicability limit for the local tangential jump, m.
    pub maximum_slip_m: f64,
}

/// One interface at the initial accepted mechanical/thermal stage. The point
/// may already carry damage; it is retained, never replaced by a virgin point.
#[derive(Debug, Clone, Copy)]
pub struct Interface<'a> {
    /// Original area/endpoints and the separately supplied bonded conductance.
    pub reference: ThermalLink,
    pub point: CohesivePoint<'a>,
    pub closed: bool,
    pub policy: ContactPolicy,
}
impl Interface<'_> {
    fn next_closed(&self, point: CohesivePoint<'_>) -> Result<bool, Error> {
        let p = self.policy;
        let j = point.jump_m();
        // Normalize before squaring to avoid under/overflow of metre squares.
        let a = j[1].abs() / p.maximum_slip_m;
        let b = j[2].abs() / p.maximum_slip_m;
        if !a.is_finite() || !b.is_finite() || a*a+b*b > 1.0 {
            return Err(Error::Invalid("tangential travel exceeds the fixed-contact-area model"));
        }
        Ok(if !point.failed() { true }
            else if j[0] >= p.open_gap_m { false }
            else if j[0] <= p.reclose_gap_m { true }
            else { self.closed })
    }

    fn thermal_link(&self) -> Option<ThermalLink> {
        self.closed.then(|| ThermalLink {
            conductance_w_k: if self.point.failed() { self.policy.reclosed_conductance_w_k }
                else { self.reference.conductance_w_k },
            ..self.reference
        })
    }
}

#[derive(Debug, Clone)]
pub enum Error {
    Thermal(DepositionError),
    Material { interface: u64, error: CohesiveError },
    Invalid(&'static str),
}
impl From<DepositionError> for Error {
    fn from(error: DepositionError) -> Self { Self::Thermal(error) }
}
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "fracture/thermal contact refused: {self:?}")
    }
}
impl std::error::Error for Error {}

/// Observations for this endpoint, not a closed thermodynamic balance.
#[derive(Debug, Clone, Copy, Default)]
pub struct ContactReport {
    pub opened: usize,
    pub reclosed: usize,
    pub changed_conductances: usize,
    /// Incremental fracture dissipation over the bound original areas, J.
    /// Reported only: not automatically added to thermal enthalpy.
    pub fracture_dissipation_j: f64,
}

#[derive(Clone)]
struct Contacts<'a>(BTreeMap<u64, Interface<'a>>);
impl DepositionObserver for Contacts<'_> {
    type Error = Error;
    fn event(&mut self, event: &Event, _: &DepositionSimulation<'_>) -> Result<(), Error> {
        let id = match event { Event::Contact(link) => Some(link.id),
            Event::Separation(id) => Some(*id), _ => None };
        if id.is_some_and(|id| self.0.contains_key(&id)) {
            return Err(Error::Invalid("scheduled contact change conflicts with mechanical ownership"));
        }
        Ok(())
    }
    fn step(&mut self, _: f64, _: f64, _: &BTreeMap<CellId, f64>, _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Error> {
        checkpoint()?;
        Ok(())
    }
}

/// Accepted thermal state and mechanical contact endpoints. Thermal transport
/// uses the previous accepted closure, then `synchronize` changes the following
/// interval. A caller stages a structural solve and commits that proposal only
/// after synchronization succeeds. Arbitrary material-law changes are refused.
#[derive(Clone)]
pub struct FractureCooling<'t, 'm> {
    thermal: DepositionSimulation<'t>,
    contacts: Contacts<'m>,
}
impl<'t, 'm> FractureCooling<'t, 'm> {
    /// Bind already activated thermal cells and exactly matching current links.
    /// Already-open/damaged interfaces may be imported with their real history.
    /// Construction and all sweeps are bounded by `max_interfaces`; cloning is
    /// bounded by the thermal owner's budget, but is not internally cancellable.
    pub fn new(thermal: DepositionSimulation<'t>, interfaces: Vec<Interface<'m>>,
        max_interfaces: usize, mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        checkpoint()?;
        if interfaces.is_empty() || interfaces.len() > max_interfaces {
            return Err(Error::Invalid("interface count is empty or exceeds the declared budget"));
        }
        let mut contacts = BTreeMap::new();
        let mut updates = Vec::with_capacity(interfaces.len());
        for interface in interfaces {
            checkpoint()?;
            let link = interface.reference;
            let p = interface.policy;
            if !positive(link.area_m2) || !nonnegative(link.conductance_w_k)
                || link.a == link.b || !thermal.cells().contains_key(&link.a)
                || !thermal.cells().contains_key(&link.b)
                || !positive(p.open_gap_m) || !nonnegative(p.reclose_gap_m)
                || p.reclose_gap_m >= p.open_gap_m
                || p.open_gap_m > interface.point.law().maximum_separation_m()
                || !nonnegative(p.reclosed_conductance_w_k) || !positive(p.maximum_slip_m)
            { return Err(Error::Invalid("invalid contact binding or explicit closure policy")); }
            if [link.a,link.b].iter().any(|id| link.area_m2 > thermal.cells()[id].geometry.surface_area_m2) {
                return Err(Error::Invalid("original contact area exceeds an incident cell surface"));
            }
            if interface.next_closed(interface.point)? != interface.closed {
                return Err(Error::Invalid("initial closure contradicts the supplied mechanical endpoint"));
            }
            let expected = interface.thermal_link();
            updates.push(ContactUpdate { id: link.id, expected, replacement: expected });
            if contacts.insert(link.id, interface).is_some() {
                return Err(Error::Invalid("duplicate mechanical interface identity"));
            }
        }
        // Also checks that due scheduled events have been processed first.
        let thermal = thermal.replace_contacts(&updates, &mut checkpoint)?;
        Ok(Self { thermal, contacts: Contacts(contacts) })
    }

    pub fn thermal(&self) -> &DepositionSimulation<'t> { &self.thermal }
    pub fn interfaces(&self) -> &BTreeMap<u64, Interface<'m>> { &self.contacts.0 }

    /// Advance heat with accepted closure held fixed. A scheduled event touching
    /// any mechanically owned link refuses, even if another equal-time event
    /// would restore its old value. Unrelated births/boundary events still work.
    pub fn advance_to(&self, target_s: f64, config: StepConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, EnergyAudit), Error> {
        let (thermal, contacts, audit) = self.thermal.advance_observed(
            &self.contacts, target_s, config, checkpoint)?;
        Ok((Self { thermal, contacts }, audit))
    }

    /// Accept exactly one proposed point for every bound ID at this thermal time.
    /// Points from `fs-solid::cohesive` can be passed directly. Re-evaluating their
    /// straight jump increments against the retained history verifies the SAME
    /// damage path and immutable law; it is not a second equilibrium solve.
    /// Split paths at intermediate reversals: a hidden larger excursion refuses.
    ///
    /// Missing/extra/duplicate IDs, stale time, foreign laws, lost damage,
    /// thermal-contact conflicts, excessive slip, or cancellation refuse the
    /// whole proposal. Reclosing a failed interface restores only heat transfer,
    /// never a tensile/shear bond. No thermal time passes in this operation.
    pub fn synchronize(&self, time_s: f64, proposed: &[(u64, CohesivePoint<'m>)],
        local_energy_tolerance_j_m2: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, ContactReport), Error> {
        checkpoint()?;
        if time_s != self.thermal.time_s() || proposed.len() != self.contacts.0.len()
            || !positive(local_energy_tolerance_j_m2)
        { return Err(Error::Invalid("require current thermal time, complete endpoints and positive energy gate")); }
        let mut contacts = BTreeMap::new();
        let mut updates = Vec::with_capacity(proposed.len());
        let mut report = ContactReport::default();
        // Sorted accumulation makes output independent of caller arrival order.
        let mut ordered = BTreeMap::new();
        for &(id, point) in proposed {
            checkpoint()?;
            if !self.contacts.0.contains_key(&id) || ordered.insert(id, point).is_some() {
                return Err(Error::Invalid("unknown or duplicate mechanical endpoint"));
            }
        }
        for (id, point) in ordered {
            checkpoint()?;
            let old = self.contacts.0[&id];
            if !std::ptr::eq(old.point.law(), point.law()) {
                return Err(Error::Invalid("mechanical endpoint belongs to a different immutable law"));
            }
            let trial = old.point.trial(point.jump_m(), local_energy_tolerance_j_m2)
                .map_err(|error| Error::Material { interface: id, error })?;
            if trial.point != point {
                return Err(Error::Invalid("mechanical endpoint lost history or contains an unresolved excursion"));
            }
            let closed = old.next_closed(point)?;
            let next = Interface { point, closed, ..old };
            let expected = old.thermal_link();
            let replacement = next.thermal_link();
            if old.closed && !closed { report.opened += 1; }
            if !old.closed && closed { report.reclosed += 1; }
            if old.closed && closed && expected != replacement { report.changed_conductances += 1; }
            report.fracture_dissipation_j += old.reference.area_m2 * trial.energy.dissipated;
            if !nonnegative(report.fracture_dissipation_j) {
                return Err(Error::Invalid("unrepresentable aggregate fracture dissipation"));
            }
            updates.push(ContactUpdate { id, expected, replacement });
            contacts.insert(id, next);
        }
        let thermal = self.thermal.replace_contacts(&updates, &mut checkpoint)?;
        checkpoint()?;
        Ok((Self { thermal, contacts: Contacts(contacts) }, report))
    }

    /// Executor-scoped endpoint synchronization with identical atomic semantics.
    pub fn synchronize_in(&self, cx: &Cx<'_>, time_s: f64,
        proposed: &[(u64, CohesivePoint<'m>)], local_energy_tolerance_j_m2: f64,
    ) -> Result<(Self, ContactReport), Error> {
        self.synchronize(time_s, proposed, local_energy_tolerance_j_m2,
            || cx.checkpoint().map_err(|_| DepositionError::Cancelled))
    }
}
fn positive(v: f64) -> bool { v.is_finite() && v > 0.0 }
fn nonnegative(v: f64) -> bool { v.is_finite() && v >= 0.0 }

#[cfg(test)]
mod tests;
