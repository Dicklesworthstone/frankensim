//! Compose independent physical owners on the SAME enthalpy substeps.
//!
//! Source rates add; they are not applied by running separate thermal solves.
//! Each link has at most one conductance override owner. Both owners see the
//! same start state, events, accepted duration and end state. Nest this pair to
//! combine more domains without adding a scheduler, history or caloric law.
//!
//! Atomicity comes from `DepositionSimulation::advance_observed`, which clones
//! the whole pair. As with every observer, Clone must not share mutable physical
//! history or irreversible external effects. Composition is explicit staggering,
//! not a strongly coupled nonlinear solve or an accuracy/stability certificate.

use std::collections::BTreeMap;
use crate::enthalpy::{CellId, DepositionError, DepositionObserver, DepositionSimulation, Event};

/// Two independent observers; their order is explicit and deterministic.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledObservers<A, B> {
    /// First physical owner, also first in event and substep notification order.
    pub first: A,
    /// Second physical owner. It sees the same thermal state, not the first's heat alone.
    pub second: B,
}

impl<A, B> CoupledObservers<A, B> {
    /// Combine owners without evaluating or cloning their state.
    pub const fn new(first: A, second: B) -> Self { Self { first, second } }

    /// Recover the two accepted histories after a successful common advance.
    pub fn into_parts(self) -> (A, B) { (self.first, self.second) }
}

/// Retain the failing physical owner instead of erasing its error into text.
#[derive(Clone, Debug, PartialEq)]
pub enum ObserverCouplingError<A, B> {
    /// Shared thermal, field-admission, arithmetic or cancellation failure.
    Thermal(DepositionError),
    /// Failure of the first physical owner.
    First(A),
    /// Failure of the second physical owner.
    Second(B),
    /// Two owners tried to replace the same link's total conductance.
    /// Even equal values refuse: ownership must not depend on coincident data.
    ConflictingConductance(u64),
}
impl<A, B> From<DepositionError> for ObserverCouplingError<A, B> {
    fn from(error: DepositionError) -> Self { Self::Thermal(error) }
}
impl<A: core::fmt::Display, B: core::fmt::Display> core::fmt::Display for ObserverCouplingError<A, B> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Thermal(e) => write!(f, "coupled thermal transport: {e}"),
            Self::First(e) => write!(f, "first thermal observer: {e}"),
            Self::Second(e) => write!(f, "second thermal observer: {e}"),
            Self::ConflictingConductance(id) => write!(f, "multiple conductance owners for thermal link {id}"),
        }
    }
}
impl<A: std::error::Error + 'static, B: std::error::Error + 'static> std::error::Error for ObserverCouplingError<A, B> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Thermal(e) => Some(e), Self::First(e) => Some(e), Self::Second(e) => Some(e),
            Self::ConflictingConductance(_) => None,
        }
    }
}

fn admit_sources(
    values: &BTreeMap<CellId, f64>, thermal: &DepositionSimulation<'_>,
    check: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), DepositionError> {
    if values.len() > thermal.cells().len() {
        return Err(DepositionError::InvalidInput("too many composed heat-source cells"));
    }
    for (index, (id, watts)) in values.iter().enumerate() {
        if index % 256 == 0 { check()?; }
        if !thermal.cells().contains_key(id) {
            return Err(DepositionError::InvalidInput("inactive composed heat-source cell"));
        }
        if !watts.is_finite() { return Err(DepositionError::Arithmetic); }
    }
    Ok(())
}
fn admit_conductances(
    values: &BTreeMap<u64, f64>, thermal: &DepositionSimulation<'_>,
    check: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), DepositionError> {
    if values.len() > thermal.links().len() {
        return Err(DepositionError::InvalidInput("too many composed thermal links"));
    }
    for (index, (id, value)) in values.iter().enumerate() {
        if index % 256 == 0 { check()?; }
        if !thermal.links().contains_key(id) || !value.is_finite() || *value < 0.0 {
            return Err(DepositionError::InvalidInput("invalid composed conductance or inactive link"));
        }
    }
    Ok(())
}

impl<A: DepositionObserver, B: DepositionObserver> DepositionObserver for CoupledObservers<A, B> {
    type Error = ObserverCouplingError<A::Error, B::Error>;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.first.begin().map_err(ObserverCouplingError::First)?;
        self.second.begin().map_err(ObserverCouplingError::Second)
    }

    fn sources(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check()?;
        let mut result = self.first.sources(thermal, check).map_err(ObserverCouplingError::First)?;
        admit_sources(&result, thermal, check)?;
        check()?;
        let other = self.second.sources(thermal, check).map_err(ObserverCouplingError::Second)?;
        admit_sources(&other, thermal, check)?;
        for (index, (id, watts)) in other.into_iter().enumerate() {
            if index % 256 == 0 { check()?; }
            match result.get_mut(&id) {
                Some(total) if watts != 0.0 => {
                    *total += watts;
                    if !total.is_finite() { return Err(DepositionError::Arithmetic.into()); }
                }
                Some(_) => {},
                None => { result.insert(id, watts); },
            }
        }
        check()?;
        Ok(result)
    }

    fn conductances(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<u64, f64>, Self::Error> {
        check()?;
        let mut result = self.first.conductances(thermal, check).map_err(ObserverCouplingError::First)?;
        admit_conductances(&result, thermal, check)?;
        check()?;
        let other = self.second.conductances(thermal, check).map_err(ObserverCouplingError::Second)?;
        admit_conductances(&other, thermal, check)?;
        for (index, (id, value)) in other.into_iter().enumerate() {
            if index % 256 == 0 { check()?; }
            if result.insert(id, value).is_some() {
                return Err(ObserverCouplingError::ConflictingConductance(id));
            }
        }
        check()?;
        Ok(result)
    }

    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.first.event(event, thermal).map_err(ObserverCouplingError::First)?;
        self.second.event(event, thermal).map_err(ObserverCouplingError::Second)
    }

    fn step(&mut self, start_s: f64, duration_s: f64,
        old: &BTreeMap<CellId, f64>, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Self::Error> {
        check()?;
        self.first.step(start_s, duration_s, old, thermal, check).map_err(ObserverCouplingError::First)?;
        check()?;
        self.second.step(start_s, duration_s, old, thermal, check).map_err(ObserverCouplingError::Second)?;
        check()?;
        Ok(())
    }
}

#[cfg(test)]
#[path = "observers/tests.rs"]
mod tests;
