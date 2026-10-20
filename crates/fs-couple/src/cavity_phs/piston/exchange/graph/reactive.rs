//! Shared-snapshot flow, declared chemistry and material wall motion.
//!
//! Each nozzle sees the same accepted state. The existing conservative graph
//! first transports species and enthalpy. Each selected volume then reacts ONCE
//! using its actual mixed state, at fixed mass/volume/total internal energy.
//! Finally the existing piston owner advances that reacted gas and its wall.
//! Source formation energy stays inside internal energy, never an extra heat
//! credit. The original all-owner balance spans all three stages.
//!
//! This is a first-order flow -> reaction -> motion operator split over ONE
//! duration. Incoming species can react now; products leave on the next flow
//! interval. Rates are not re-evaluated in mechanical Newton trials. Admission
//! and conservation do not certify time accuracy, positivity at arbitrary dt,
//! stiff kinetics, flames, shocks, or a measured combustion mechanism.
//! Control-volume convention: https://cantera.org/stable/reference/reactors/ideal-gas-reactor.html

use super::{BTreeMap, ClosedGasParcel, Error, GasEnergySolve, GasPiston,
    MovingGasCandidate, MovingGasDrive, MovingGasSolve, MovingGasStep,
    SnapshotNozzleLink, check};
use crate::cavity_phs::piston::reactive::GasReactionError;
// Reuse the original declared reaction binding. No second source/kinetic schema.
pub use crate::cavity_phs::piston::thermal::boundary::ports::network::material::reactive::ChamberReaction;
pub use crate::cavity_phs::piston::reactive::GasReactionStep;

/// Complete chemical and flow/motion proposal. Inputs and earlier port windows
/// remain untouched, including if a later volume or final balance refuses.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactiveMovingGasCandidate<P = GasPiston> {
    motion: MovingGasCandidate<P>,
    reactions: BTreeMap<u64, GasReactionStep>,
}
impl<P> ReactiveMovingGasCandidate<P> {
    /// All moving chambers retaining their original physical owner.
    #[must_use]
    pub fn pistons(&self) -> &BTreeMap<u64, P> { self.motion.pistons() }
    /// Finite fixed-volume gas stores, which may also have declared reactions.
    #[must_use]
    pub fn reservoirs(&self) -> &BTreeMap<u64, ClosedGasParcel> { self.motion.reservoirs() }
    /// Actual post-flow reaction evaluations, including zero-extent results.
    #[must_use]
    pub fn reactions(&self) -> &BTreeMap<u64, GasReactionStep> { &self.reactions }
    /// Flow and wall reports plus the fresh ORIGINAL-to-final whole balance.
    /// `flow` accounts for transport only; `pistons` start after chemistry.
    /// Gas storage and total residual span flow, chemistry AND wall motion.
    #[must_use]
    pub const fn step(&self) -> &MovingGasStep { self.motion.step() }
    /// Accept all states together. Reaction extent and enthalpy are already
    /// applied; consuming these reports must not deposit either again.
    #[must_use]
    pub fn into_parts(self) -> (BTreeMap<u64, P>, BTreeMap<u64, ClosedGasParcel>,
        MovingGasStep, BTreeMap<u64, GasReactionStep>) {
        let (pistons, reservoirs, step) = self.motion.into_parts();
        (pistons, reservoirs, step, self.reactions)
    }
    /// Retain a pre-existing material/inertia binding around each numeric wall.
    pub(crate) fn map_pistons<Q>(self, wrap: impl FnMut(u64, P) -> Q)
        -> ReactiveMovingGasCandidate<Q>
    {
        ReactiveMovingGasCandidate {
            motion: self.motion.map_pistons(wrap), reactions: self.reactions,
        }
    }
}

/// Couple an optional declared reaction to each moving OR fixed gas volume.
/// The existing `ChamberReaction::chamber` is the volume's graph id here.
/// Duplicate/foreign bindings refuse; unlisted volumes stay nonreacting.
/// Catalogs must cover every species that can enter the selected volume, with
/// exactly matching caloric, transport and source-pack definitions. Products
/// need not be initially present. No species or material is inferred by name.
///
/// Flow, species, chemistry consumption/entropy, wall-clock/impulse and global
/// energy gates all remain active. Every chemical step uses `duration_s`, not
/// an independently advancing clock. A source/rate discontinuity requires a
/// caller-selected split. Shorten and refine intervals after a refusal.
///
/// Gross mass exchange still makes a chamber's two-port wall/heat-only audit
/// unavailable; its enclosing graph balance remains independent. A closed
/// reacting chamber retains a complete one-tick account including its actual
/// chemical entropy change. Empty reactions are the original numeric path.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    pistons: &BTreeMap<u64, GasPiston>, reservoirs: &BTreeMap<u64, ClosedGasParcel>,
    links: &[SnapshotNozzleLink], duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>, drives: &BTreeMap<u64, MovingGasDrive>,
    reactions: &[ChamberReaction<'_>], controls: MovingGasSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<ReactiveMovingGasCandidate, Error> {
    let (motion, reactions) = super::staged_candidate(pistons, reservoirs, links,
        duration_s, endpoints, drives, controls, reactions, cancelled)?;
    Ok(ReactiveMovingGasCandidate { motion, reactions })
}

pub(super) fn bindings<'a, 's>(
    initial: &BTreeMap<u64, ClosedGasParcel>, reactions: &'a [ChamberReaction<'s>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<BTreeMap<u64, &'a ChamberReaction<'s>>, Error> {
    if reactions.len() > initial.len() {
        return Err(Error::Invalid("at most one reaction per existing gas volume required"));
    }
    let mut ordered = BTreeMap::new();
    for binding in reactions {
        check(cancelled)?;
        if !initial.contains_key(&binding.chamber)
            || ordered.insert(binding.chamber, binding).is_some() {
            return Err(Error::Invalid("unique reaction bindings to existing gas volumes required"));
        }
    }
    Ok(ordered)
}

pub(super) fn react(
    gases: &mut BTreeMap<u64, ClosedGasParcel>,
    reactions: &BTreeMap<u64, &ChamberReaction<'_>>, duration_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<BTreeMap<u64, GasReactionStep>, Error> {
    let mut reports = BTreeMap::new();
    for (&id, binding) in reactions {
        check(cancelled)?;
        let gas = gases.get(&id).ok_or(Error::Invalid("flow omitted a reacting volume"))?;
        let proposal = binding.reaction.parcel_candidate(gas, duration_s,
            binding.controls, &mut *cancelled).map_err(|error| match error {
                GasReactionError::Cancelled => Error::Cancelled,
                error => Error::Reaction { id, error },
            })?;
        if proposal.gas().mass_kg().to_bits() != gas.mass_kg().to_bits()
            || proposal.gas().volume_m3().to_bits() != gas.volume_m3().to_bits() {
            return Err(Error::Invalid("closed chemical stage changed transported mass or volume"));
        }
        let (gas, report) = proposal.into_parts();
        gases.insert(id, gas);
        reports.insert(id, report);
    }
    Ok(reports)
}
