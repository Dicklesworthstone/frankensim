//! Bond development on the actual explicit-deposition thermal timeline.
//!
//! Each bound interface has its own material history and actual contact area.
//! Temperature is an explicitly selected weighted mean of the two cell values,
//! linearly interpolated on each accepted thermal substep. This is a reduced
//! interface-temperature assumption, not a resolved through-contact field.
//! Exposure does not feed heat back, change conductance or imply bond strength.

use std::collections::BTreeMap;
use fs_exec::Cx;
use fs_material::healing::{ContactHistory, HealingError, HealingLaw};
use fs_material::healing::strength::{FrozenBond, StrengthMap};
use crate::transient::deposition::{CellId, DepositionBudget, DepositionError,
    DepositionObserver, DepositionSimulation, EnergyAudit, Event, MaterialId,
    StepConfig, ThermalLaw, ThermalLink, TimedEvent};

/// Select which real thermal interface is allowed to develop a mechanical bond.
/// Unbound thermal interfaces conduct heat but acquire no implicit bond.
#[derive(Debug, Clone, Copy)]
pub struct BondBinding<'a> {
    pub interface_id: u64,
    pub a: CellId,
    pub b: CellId,
    pub law: &'a dyn HealingLaw,
    /// T_interface = weight_a*T_a + (1-weight_a)*T_b, with weight in [0,1].
    pub temperature_weight_a: f64,
}

/// Retained memory and per-advance constitutive-work admission.
#[derive(Debug, Clone, Copy)]
pub struct BondBudget {
    /// Includes declarations and retained separated-contact histories.
    pub max_interfaces: usize,
    /// Number of history updates in one complete candidate advance.
    pub max_updates: usize,
}

#[derive(Debug, Clone)]
pub enum BondingError {
    Thermal(DepositionError),
    Healing { interface: u64, error: HealingError },
    Invalid(&'static str),
    Budget,
}
impl From<DepositionError> for BondingError {
    fn from(error: DepositionError) -> Self { Self::Thermal(error) }
}
impl core::fmt::Display for BondingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f,"deposition bond history refused: {self:?}")
    }
}
impl std::error::Error for BondingError {}

/// One original-area contact cohort. Closed histories are retained and never
/// automatically healed/rebonded. A new contact requires a NEW interface ID.
#[derive(Debug, Clone)]
pub struct BondRecord<'a> {
    contact: ThermalLink,
    first_contact_s: f64,
    separated_at_s: Option<f64>,
    history: ContactHistory<'a>,
}
impl<'a> BondRecord<'a> {
    pub fn contact(&self) -> &ThermalLink { &self.contact }
    pub fn first_contact_s(&self) -> f64 { self.first_contact_s }
    pub fn separated_at_s(&self) -> Option<f64> { self.separated_at_s }
    pub fn history(&self) -> &ContactHistory<'a> { &self.history }
}

/// A mechanical handoff carries actual original area and the explicit strength
/// map. It is for a new mechanical stage, never resetting prior cohesive damage.
#[derive(Debug, Clone)]
pub struct FrozenInterface {
    pub interface_id: u64,
    pub a: CellId,
    pub b: CellId,
    pub area_m2: f64,
    pub bond: FrozenBond,
}

#[derive(Debug, Clone)]
struct BondTracker<'a> {
    bindings: BTreeMap<u64, BondBinding<'a>>,
    records: BTreeMap<u64, BondRecord<'a>>,
    budget: BondBudget,
    updates: usize,
}
impl DepositionObserver for BondTracker<'_> {
    type Error = BondingError;
    fn begin(&mut self) -> Result<(),Self::Error> { self.updates = 0; Ok(()) }

    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>) -> Result<(),Self::Error> {
        match event {
            Event::Contact(link) => {
                let Some(binding) = self.bindings.get(&link.id) else { return Ok(()); };
                if binding.a != link.a || binding.b != link.b {
                    return Err(BondingError::Invalid("contact endpoints differ from healing binding"));
                }
                // Thermal IDs may be reused after separation; a mechanical
                // material cohort may not silently change area or lose history.
                if self.records.contains_key(&link.id) {
                    return Err(BondingError::Invalid("recontact requires a new bond interface identity"));
                }
                if self.records.len() >= self.budget.max_interfaces { return Err(BondingError::Budget); }
                let history = ContactHistory::new(binding.law)
                    .map_err(|error| BondingError::Healing { interface: link.id, error })?;
                self.records.insert(link.id, BondRecord { contact: *link,
                    first_contact_s: thermal.time_s(), separated_at_s: None, history });
            }
            Event::Separation(id) => {
                if let Some(record) = self.records.get_mut(id) {
                    record.separated_at_s = Some(thermal.time_s());
                }
            }
            _ => {}
        }
        Ok(())
    }

    fn step(&mut self, _start_s: f64, duration_s: f64,
        old_temperatures_k: &BTreeMap<CellId,f64>, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(),DepositionError>,
    ) -> Result<(),Self::Error> {
        for (index, (&id,record)) in self.records.iter_mut().enumerate() {
            if index % 256 == 0 { checkpoint()?; }
            if record.separated_at_s.is_some() { continue; }
            if thermal.links().get(&id) != Some(&record.contact) {
                return Err(BondingError::Invalid("active thermal/bond interface state differs"));
            }
            if self.updates >= self.budget.max_updates { return Err(BondingError::Budget); }
            let weight = self.bindings[&id].temperature_weight_a;
            let old_a = *old_temperatures_k.get(&record.contact.a)
                .ok_or(BondingError::Invalid("missing pre-step contact temperature"))?;
            let old_b = *old_temperatures_k.get(&record.contact.b)
                .ok_or(BondingError::Invalid("missing pre-step contact temperature"))?;
            let new_a = thermal.temperature_k(record.contact.a)?;
            let new_b = thermal.temperature_k(record.contact.b)?;
            let start_k = weight*old_a + (1.0-weight)*old_b;
            let end_k = weight*new_a + (1.0-weight)*new_b;
            record.history = record.history.advance(start_k,end_k,duration_s)
                .map_err(|error| BondingError::Healing { interface: id, error })?;
            self.updates += 1;
        }
        checkpoint()?;
        Ok(())
    }
}

/// Thermal and bond states are advanced and published together. Only immutable
/// access is exposed so detached, stale thermal state cannot age these bonds.
#[derive(Clone)]
pub struct BondedDeposition<'a> {
    thermal: DepositionSimulation<'a>,
    tracker: BondTracker<'a>,
}
impl<'a> BondedDeposition<'a> {
    /// Starts empty at t=0; use `advance_to(0)` to process initial births/contact.
    pub fn new(schedule: Vec<TimedEvent>, materials: &[(MaterialId,&'a dyn ThermalLaw)],
        thermal_budget: DepositionBudget, bindings: Vec<BondBinding<'a>>, bond_budget: BondBudget,
    ) -> Result<Self,BondingError> {
        if bindings.len() > bond_budget.max_interfaces { return Err(BondingError::Budget); }
        let mut bound = BTreeMap::new();
        for binding in bindings {
            if binding.a == binding.b || !binding.temperature_weight_a.is_finite()
                || !(0.0..=1.0).contains(&binding.temperature_weight_a) {
                return Err(BondingError::Invalid("healing requires distinct endpoints and weight in [0,1]"));
            }
            ContactHistory::new(binding.law).map_err(|error|
                BondingError::Healing { interface: binding.interface_id, error })?;
            if bound.insert(binding.interface_id,binding).is_some() {
                return Err(BondingError::Invalid("duplicate healing interface binding"));
            }
        }
        Ok(Self { thermal: DepositionSimulation::new(schedule,materials,thermal_budget)?,
            tracker: BondTracker { bindings: bound, records: BTreeMap::new(), budget: bond_budget, updates: 0 } })
    }
    pub fn thermal(&self) -> &DepositionSimulation<'a> { &self.thermal }
    pub fn records(&self) -> &BTreeMap<u64,BondRecord<'a>> { &self.tracker.records }

    /// Checkpoint-driven form; actual event/substep order is the thermal owner's.
    pub fn advance_to(&self, target_s: f64, config: StepConfig,
        checkpoint: impl FnMut() -> Result<(),DepositionError>,
    ) -> Result<(Self,EnergyAudit,usize),BondingError> {
        let (thermal,tracker,audit) = self.thermal.advance_observed(&self.tracker,target_s,config,checkpoint)?;
        let updates = tracker.updates;
        Ok((Self { thermal,tracker },audit,updates))
    }
    /// Executor integration; failed/cancelled steps commit neither heat nor dose.
    pub fn advance_in(&self, cx: &Cx<'_>, target_s: f64, config: StepConfig)
        -> Result<(Self,EnergyAudit,usize),BondingError> {
        self.advance_to(target_s,config,|| cx.checkpoint().map_err(|_| DepositionError::Cancelled))
    }
    pub fn freeze_bond(&self, id: u64, map: &StrengthMap) -> Result<FrozenInterface,BondingError> {
        let record = self.records().get(&id)
            .ok_or(BondingError::Invalid("an interface that never contacted has no bond"))?;
        if record.separated_at_s.is_some() {
            return Err(BondingError::Invalid("cannot create a bonded mechanical stage from a separated interface"));
        }
        let bond = map.freeze(&record.history).map_err(|error| BondingError::Healing { interface: id,error })?;
        Ok(FrozenInterface { interface_id: id, a: record.contact.a, b: record.contact.b,
            area_m2: record.contact.area_m2, bond })
    }
}

#[cfg(test)]
mod tests;
