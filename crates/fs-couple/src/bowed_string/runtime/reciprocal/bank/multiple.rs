//! Several material-derived bow contacts on the SAME retained string/body state.
//!
//! Contact shapes come from the selected material strings, not supplied modal
//! vectors. The original simultaneous-friction owner computes all reactions,
//! work, sliding loss and pressure together. This module adds no force law or
//! time integration. Adding contacts preserves existing vibration and clocks.
//!
//! Loads and surface speeds are prescribed independently per contact. A single
//! bow spanning several strings still needs an external normal-load model;
//! this adapter does not infer hair compliance or redistribute normal forces.
//! Every contact station is fixed. The underlying bounded nonlinear solve may
//! refuse a step; global uniqueness and measured instrument timbre are not claimed.

use super::*;
use crate::render::schedule::force::coupled::contact::friction::ModalFrictionPort;
use crate::render::schedule::force::coupled::contact::friction::multiple::{
    MultiFrictionConfig, MultiFrictionFrame, MultiFrictionModalSystem,
};

/// One physical contact on a selected material string and its original law.
/// The law retains caller authority; string material data do not certify rosin.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialBowContact {
    /// Original string index and fixed interior station along its full span.
    pub location: StringBowLocation,
    /// Unmodified dry-friction law for this prescribed moving surface.
    pub law: FrictionLaw,
    /// Original per-contact physical and nonlinear-solve limits.
    pub friction: ModalFrictionConfig,
}

/// Retained material, body, cavity and gas bindings plus all contact inputs.
#[derive(Clone, Debug)]
pub struct MultiBowStringBankSource {
    bank: MaterialStringBankSource,
    bows: Vec<MaterialBowContact>,
    joint: MultiFrictionConfig,
}
impl MultiBowStringBankSource {
    /// Original strings, layers, regional/uniform body, cavity, gas and supports.
    #[must_use]
    pub const fn bank(&self) -> &MaterialStringBankSource { &self.bank }
    /// Actual drive order: the original bow FIRST, then the added contacts.
    #[must_use]
    pub fn bows(&self) -> &[MaterialBowContact] { &self.bows }
    /// Original simultaneous-friction owner's declared setup and sweep limits.
    #[must_use]
    pub const fn joint(&self) -> &MultiFrictionConfig { &self.joint }
}

/// One mechanical state and pressure observer with jointly solved bow contacts.
/// A released bow removes only traction; no string or body vibration is muted.
pub struct MultiBowStringBank {
    source: MultiBowStringBankSource,
    cards: Vec<BowedStringCard>,
    system: MultiFrictionModalSystem,
}
impl ReciprocalStringBank {
    /// Add material-derived contacts without reducing or resetting any component.
    ///
    /// The original bow remains port zero. Additional locations use the original
    /// string declaration order and its actual material-derived displacement
    /// basis. They may contact different strings or distinct stations on the
    /// same string. All string, bridge, body and cavity reactions participate in
    /// the original owner's simultaneous solve, not sequential voice steps.
    ///
    /// Existing modal state, accepted clock and physical pressure weights survive
    /// the conversion. Empty additions use the original single-port step exactly.
    /// No source claims, state histories or observer normalizations are replaced.
    ///
    /// # Errors
    /// Unknown strings, endpoint/nonfinite stations, port/setup limits, invalid
    /// friction laws or per-contact budgets, numerical response or cancellation.
    pub fn with_additional_bows(
        self, additional: Vec<MaterialBowContact>, joint: MultiFrictionConfig,
        gate: &CancelGate,
    ) -> Result<MultiBowStringBank, ReciprocalStringError> {
        poll(gate)?;
        let count = additional.len().checked_add(1)
            .ok_or_else(|| invalid("material bow count overflow"))?;
        // Enforce the existing owner's port cap before allocating modal shapes.
        if !(1..=64).contains(&joint.max_ports) || count > joint.max_ports {
            return Err(invalid("material bow contacts exceed the joint port budget"));
        }
        for bow in &additional {
            if bow.location.string >= self.cards.len() || !interior(bow.location.station_fraction) {
                return Err(invalid("each added bow must name an actual material string and interior station"));
            }
        }
        let mut ports = Vec::with_capacity(additional.len());
        for bow in &additional {
            poll(gate)?;
            ports.push((ModalFrictionPort {
                left: ModalAttachment { component: bow.location.string,
                    shapes: shape(&self.cards[bow.location.string], bow.location.station_fraction) },
                right: None, law: bow.law.clone(),
            }, bow.friction));
        }
        let mut bows = Vec::with_capacity(count);
        bows.push(MaterialBowContact { location: self.source.options.bow,
            law: self.system.port().law.clone(), friction: self.source.options.friction });
        bows.extend(additional);
        let system = self.system.with_additional_ports(ports, joint, gate)
            .map_err(ReciprocalStringError::Friction)?;
        poll(gate)?;
        Ok(MultiBowStringBank { source: MultiBowStringBankSource {
            bank: self.source, bows, joint,
        }, cards: self.cards, system })
    }
}
impl MultiBowStringBank {
    /// Original material/geometry bindings and every contact in drive order.
    #[must_use]
    pub const fn source(&self) -> &MultiBowStringBankSource { &self.source }
    /// Original material-derived coefficients; source receipts remain attached.
    #[must_use]
    pub fn cards(&self) -> &[BowedStringCard] { &self.cards }
    /// Read-only state, physical pressure weights and per-contact work reports.
    #[must_use]
    pub const fn system(&self) -> &MultiFrictionModalSystem { &self.system }
    /// One physical sample for ALL strings, bodies and contacts. The drive slice
    /// must match `source().bows()` exactly; zero load releases only that bow.
    pub fn step(&mut self, external: &[f64], drives: &[FrictionDrive])
        -> Result<&MultiFrictionFrame, ModalFrictionError> { self.system.step(external, drives) }
    /// Same atomic physical sample with the original owner's cancellation polls.
    pub fn step_under_gate(&mut self, external: &[f64], drives: &[FrictionDrive], gate: &CancelGate)
        -> Result<&MultiFrictionFrame, ModalFrictionError> { self.system.step_under_gate(external, drives, gate) }
    /// Transfer the original sources together with the complete accepted state.
    #[must_use]
    pub fn into_parts(self) -> (MultiBowStringBankSource, MultiFrictionModalSystem) { (self.source, self.system) }
}

/// Common gesture-clock playback of all contacts and the shared pressure observer.
pub mod performance;
