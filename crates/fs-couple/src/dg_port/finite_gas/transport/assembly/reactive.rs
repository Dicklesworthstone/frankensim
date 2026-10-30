//! Original homogeneous chemistry in a continuing finite-gas acoustic assembly.
//!
//! Material transport samples one unchanged network. Each selected mechanism
//! then solves its actual post-transport inventory, external heat and proposed
//! gas-to-wave work at one endpoint. The existing acoustic work residual closes
//! that proposal; no reaction extent, inventory or acoustic history accumulates
//! across trial evaluations. Species formation energy stays in original total
//! gas internal energy, and mechanical work is never treated as heat entropy.
//!
//! This remains a first-order split from transport and the frozen-medium native
//! acoustic sample. Arriving species may react now; reaction products become
//! available to outgoing transport next interval. It does not resolve a flame,
//! chemical acoustic relaxation, or acoustic continuity across a restriction.

/// Reciprocal finite material heat exchange in the same reacting acoustic step.
pub mod thermal;

use super::{FiniteGasAcousticInput, FiniteGasAcousticNode, FiniteGasAcousticStep};
use crate::bernoulli_aperture::wall::material::{
    FiniteGasDgTransportError, FiniteGasDgTransportInput, FiniteGasDgTransportSolve,
    FiniteGasDgTransportStep, FiniteGasTransportEvolution, FiniteGasTransportNode,
    finite_gas_transport_candidate_with_evolution,
};
use crate::cavity_phs::piston::reactive::{
    GasReaction, GasReactionError, GasReactionSolve,
    boundary::{BoundaryGasReactionStep, GasReactionBoundary},
};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::composition::GasCompositionStep;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::transport::GasTransportPaths;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasHeatStep};
use fs_thermochem::SpeciesId;
use std::collections::BTreeMap;

/// One source-bound mechanism for the parcel identified by the binding map key.
/// Unlisted parcels retain their original frozen-composition caloric evolution.
/// The catalog must include every incoming species and every reaction product
/// with the same original caloric, transport and source-pack identities.
#[derive(Clone, Copy, Debug)]
pub struct FiniteGasAcousticReaction<'a> {
    pub reaction: &'a GasReaction,
    pub controls: GasReactionSolve,
}

/// Existing transport/acoustic controls and an independent stoichiometric gate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveFiniteGasAcousticTransportSolve {
    pub transport: FiniteGasDgTransportSolve,
    /// Maximum |sum(n_end - n_start) - sum(nu * extent)| for EACH species [mol].
    /// This spans transport and every reaction, rather than summing local bounds.
    pub species_tolerance_mol: f64,
}

/// The actual accepted gas endpoint behind an original native acoustic report.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteGasAcousticReactionSource {
    /// The unchanged original heat/work caloric owner for an unbound mechanism.
    Caloric(GasCompositionStep),
    /// Original endpoint chemistry, with heat and gas-to-wave work distinguished.
    Reaction(BoundaryGasReactionStep),
}

/// Original finite-reservoir evolution; no acoustic work is invented here.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteGasAcousticReactionReservoir {
    Caloric(GasHeatStep),
    Reaction(BoundaryGasReactionStep),
}

/// Native transport and physical-owner reports, including accepted chemistry.
#[derive(Debug)]
pub struct ReactiveFiniteGasAcousticTransportStep {
    /// Its original energy account now spans transport, chemistry, external heat
    /// and complete native-owner storage. Formation energy is not a heat source.
    pub transport: FiniteGasDgTransportStep<
        FiniteGasAcousticStep<FiniteGasAcousticReactionSource>,
        FiniteGasAcousticReactionReservoir,
    >,
    pub maximum_species_residual_mol: f64,
}

impl ReactiveFiniteGasAcousticTransportStep {
    /// Accepted original chemical evaluation for either a bound acoustic owner
    /// or a finite reservoir. No rejected work-trial report is exposed here.
    #[must_use]
    pub fn reaction(&self, parcel: u64) -> Option<&BoundaryGasReactionStep> {
        if let Some(step) = self.transport.acoustics.get(&parcel) {
            return match FiniteGasAcousticNode::source::<FiniteGasAcousticReactionSource>(step) {
                FiniteGasAcousticReactionSource::Caloric(_) => None,
                FiniteGasAcousticReactionSource::Reaction(step) => Some(step),
            };
        }
        match self.transport.reservoirs.get(&parcel) {
            Some(FiniteGasAcousticReactionReservoir::Reaction(step)) => Some(step),
            Some(FiniteGasAcousticReactionReservoir::Caloric(_)) | None => None,
        }
    }
}

/// Complete unpublished endpoint set. Gas and all native owners are accepted
/// together only after chemistry, source, clock, balance and cancellation gates.
pub struct ReactiveFiniteGasAcousticTransportCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    nodes: BTreeMap<u64, FiniteGasAcousticNode>,
    step: ReactiveFiniteGasAcousticTransportStep,
}

impl ReactiveFiniteGasAcousticTransportCandidate {
    #[must_use]
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }

    #[must_use]
    pub fn nodes(&self) -> &BTreeMap<u64, FiniteGasAcousticNode> {
        &self.nodes
    }

    #[must_use]
    pub const fn step(&self) -> &ReactiveFiniteGasAcousticTransportStep {
        &self.step
    }

    /// Consume the complete proposal without changing any external owner.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, ClosedGasParcel>,
        BTreeMap<u64, FiniteGasAcousticNode>,
        ReactiveFiniteGasAcousticTransportStep,
    ) {
        (self.parcels, self.nodes, self.step)
    }
}

/// Original physical refusals and the independent whole-network species gate.
#[derive(Debug)]
pub enum ReactiveFiniteGasAcousticTransportError {
    Invalid(&'static str),
    Native(FiniteGasDgTransportError<GasReactionError>),
    SpeciesBalance {
        maximum_residual_mol: f64,
        tolerance_mol: f64,
    },
    Cancelled,
}

impl core::fmt::Display for ReactiveFiniteGasAcousticTransportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reacting finite-gas acoustic transport: {self:?}")
    }
}
impl std::error::Error for ReactiveFiniteGasAcousticTransportError {}
impl From<FiniteGasDgTransportError<GasReactionError>> for ReactiveFiniteGasAcousticTransportError {
    fn from(error: FiniteGasDgTransportError<GasReactionError>) -> Self {
        match error {
            FiniteGasDgTransportError::Cancelled => Self::Cancelled,
            error => Self::Native(error),
        }
    }
}

use ReactiveFiniteGasAcousticTransportError as Error;

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

struct Evolution<'a, 'r> {
    reactions: &'a BTreeMap<u64, FiniteGasAcousticReaction<'r>>,
}

impl FiniteGasTransportEvolution for Evolution<'_, '_> {
    type Source = FiniteGasAcousticReactionSource;
    type Reservoir = FiniteGasAcousticReactionReservoir;
    type Error = GasReactionError;

    fn is_cancelled(error: &Self::Error) -> bool {
        matches!(
            error,
            GasReactionError::Cancelled | GasReactionError::Gas(GasEnergyError::Cancelled)
        )
    }

    fn admit<C: FnMut() -> bool>(
        &self,
        id: u64,
        gas: &ClosedGasParcel,
        input: &FiniteGasDgTransportInput,
        duration_s: f64,
        cancelled: &mut C,
    ) -> Result<(), Self::Error> {
        if cancelled() {
            return Err(GasReactionError::Cancelled);
        }
        if let Some(binding) = self.reactions.get(&id) {
            let controls = binding.controls;
            if controls.gas != input.gas
                || !positive(duration_s)
                || !positive(controls.maximum_consumed_fraction)
                || controls.maximum_consumed_fraction >= 1.0
                || !positive(controls.species_tolerance_mol)
                || !positive(controls.entropy_tolerance_j_k)
                || binding.reaction.implicit_mechanism().is_none()
                || binding.reaction.implicit().is_some()
                || binding.reaction.subcycling().is_some()
            {
                return Err(GasReactionError::Invalid(
                    "shared caloric controls and one bounded coupled implicit chemistry interval required",
                ));
            }
            // Refuse an incompatible original source before material transport.
            // The native endpoint repeats this check on the actual mixed gas.
            binding.reaction.inventory(gas)?;
        }
        if cancelled() {
            Err(GasReactionError::Cancelled)
        } else {
            Ok(())
        }
    }

    fn node_candidate<C: FnMut() -> bool>(
        &self,
        id: u64,
        gas: &ClosedGasParcel,
        input: &FiniteGasDgTransportInput,
        duration_s: f64,
        work_by_gas_j: f64,
        cancelled: &mut C,
    ) -> Result<(ClosedGasParcel, Self::Source), Self::Error> {
        if let Some(binding) = self.reactions.get(&id) {
            // Every probe starts from the SAME transported gas. The existing
            // native acoustic residual determines actual gas-to-wave work;
            // only its accepted endpoint receives the entropy inequality below.
            let candidate = binding.reaction.mechanical_implicit_trial(
                gas,
                duration_s,
                GasReactionBoundary {
                    final_volume_m3: gas.volume_m3(),
                    heat_j: input.supplied_heat_j,
                    work_by_gas_j,
                },
                binding.controls,
                &mut *cancelled,
            )?;
            let (gas, step) = candidate.into_parts();
            Ok((gas, FiniteGasAcousticReactionSource::Reaction(step)))
        } else {
            let mut next = gas.clone();
            let step = next.replace_model_with_heat_work_volume(
                gas.model().clone(),
                gas.volume_m3(),
                input.supplied_heat_j,
                work_by_gas_j,
                input.gas,
                cancelled,
            )?;
            Ok((next, FiniteGasAcousticReactionSource::Caloric(step)))
        }
    }

    fn reservoir_candidate<C: FnMut() -> bool>(
        &self,
        id: u64,
        gas: &ClosedGasParcel,
        input: &FiniteGasDgTransportInput,
        duration_s: f64,
        cancelled: &mut C,
    ) -> Result<(ClosedGasParcel, Self::Reservoir), Self::Error> {
        if let Some(binding) = self.reactions.get(&id) {
            let candidate = binding.reaction.boundary_implicit_candidate(
                gas,
                duration_s,
                GasReactionBoundary {
                    final_volume_m3: gas.volume_m3(),
                    heat_j: input.supplied_heat_j,
                    work_by_gas_j: 0.0,
                },
                binding.controls,
                &mut *cancelled,
            )?;
            let (gas, step) = candidate.into_parts();
            Ok((gas, FiniteGasAcousticReactionReservoir::Reaction(step)))
        } else {
            let mut next = gas.clone();
            let step = next.add_heat(input.supplied_heat_j, input.gas, cancelled)?;
            Ok((next, FiniteGasAcousticReactionReservoir::Caloric(step)))
        }
    }

    fn admit_node_endpoint(&self, id: u64, source: &Self::Source) -> Result<(), Self::Error> {
        match (self.reactions.get(&id), source) {
            (None, FiniteGasAcousticReactionSource::Caloric(_)) => Ok(()),
            (Some(binding), FiniteGasAcousticReactionSource::Reaction(step)) => {
                // Work transports no entropy. This is the original endpoint
                // Delta S - Q/T_end, including backward-Euler numerical entropy.
                let change_j_k = step.entropy_residual_j_k;
                let tolerance_j_k = binding.controls.entropy_tolerance_j_k;
                if !change_j_k.is_finite() {
                    return Err(GasReactionError::Invalid(
                        "accepted chemical heat/work entropy is unrepresentable",
                    ));
                }
                if change_j_k < -tolerance_j_k {
                    return Err(GasReactionError::EntropyDecrease {
                        change_j_k,
                        tolerance_j_k,
                    });
                }
                Ok(())
            }
            _ => Err(GasReactionError::Invalid(
                "accepted chemical source must match its original parcel binding",
            )),
        }
    }
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}

impl Sum {
    fn add(&mut self, term: f64) -> Result<(), Error> {
        let adjusted = term - self.correction;
        let next = self.value + adjusted;
        let correction = (next - self.value) - adjusted;
        if ![term, adjusted, next, correction]
            .iter()
            .all(|value| value.is_finite())
        {
            return Err(Error::Invalid(
                "whole-network chemical amount is unrepresentable",
            ));
        }
        self.value = next;
        self.correction = correction;
        Ok(())
    }
}

fn species_residual<'a>(
    initial: &BTreeMap<u64, ClosedGasParcel>,
    final_parcels: &BTreeMap<u64, ClosedGasParcel>,
    reactions: &BTreeMap<u64, FiniteGasAcousticReaction<'_>>,
    mut reaction_step: impl FnMut(u64) -> Option<&'a BoundaryGasReactionStep>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let mut residuals: BTreeMap<SpeciesId, Sum> = BTreeMap::new();
    for (parcels, sign) in [(initial, -1.0), (final_parcels, 1.0)] {
        for gas in parcels.values() {
            poll(cancelled)?;
            let total = gas.mass_kg()
                / gas
                    .state()
                    .thermodynamics()
                    .properties()
                    .molar_mass()
                    .value();
            if !positive(total) {
                return Err(Error::Invalid(
                    "whole-network total gas amount is unrepresentable",
                ));
            }
            let model = gas.model().thermodynamics();
            for (species, fraction) in model
                .components()
                .iter()
                .zip(model.mole_composition().fractions())
            {
                poll(cancelled)?;
                let amount = total * fraction;
                if !positive(amount) {
                    return Err(Error::Invalid(
                        "whole-network species amount is unrepresentable",
                    ));
                }
                residuals
                    .entry(species.species().clone())
                    .or_default()
                    .add(sign * amount)?;
            }
        }
    }
    for &id in reactions.keys() {
        poll(cancelled)?;
        let reaction = reaction_step(id).ok_or(Error::Invalid(
            "complete endpoint lost an original chemical binding",
        ))?;
        for species in &reaction.reaction.species {
            poll(cancelled)?;
            residuals
                .entry(species.species.clone())
                .or_default()
                .add(-species.stoichiometric_change_mol)?;
        }
    }
    Ok(residuals.values().fold(0.0_f64, |maximum, residual| {
        maximum.max(residual.value.abs())
    }))
}

/// Preview one original transport snapshot followed by same-endpoint chemistry,
/// external heat and native acoustic parameter-work feedback.
///
/// Sparse bindings may select wall ducts, complete loads, moving apertures and
/// finite reservoirs. Selected mechanisms require the original coupled implicit
/// policy without chemistry-only subcycling. Unselected parcels use exactly the
/// original caloric calls. Every final flow-driving-force, native source/history,
/// chemical entropy and complete energy/mass/species gate remains active.
///
/// The independent energy account includes formation energy once, signed exterior
/// exchange and complete native mechanical/acoustic storage. Native unpartitioned
/// loss remains exported; it is not guessed to heat a reacting chamber. The
/// accepted chemical entropy residual is stage-local and does not certify the
/// complete acoustic assembly's entropy production.
#[allow(clippy::too_many_arguments)]
pub fn finite_gas_reacting_acoustic_transport_candidate_with_paths<C: FnMut() -> bool>(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    nodes: &BTreeMap<u64, FiniteGasAcousticNode>,
    inputs: &BTreeMap<u64, FiniteGasAcousticInput>,
    reactions: &BTreeMap<u64, FiniteGasAcousticReaction<'_>>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: ReactiveFiniteGasAcousticTransportSolve,
    cancelled: &mut C,
) -> Result<ReactiveFiniteGasAcousticTransportCandidate, Error> {
    poll(cancelled)?;
    if !positive(controls.species_tolerance_mol)
        || reactions.len() > parcels.len()
        || reactions.keys().any(|id| !parcels.contains_key(id))
    {
        return Err(Error::Invalid(
            "positive whole-network species tolerance and reactions bound to existing parcels required",
        ));
    }
    let candidate = finite_gas_transport_candidate_with_evolution(
        parcels,
        nodes,
        inputs,
        paths,
        duration_s,
        controls.transport,
        &Evolution { reactions },
        cancelled,
    )?;
    let (next_parcels, next_nodes, transport) = candidate.into_parts();
    let mut step = ReactiveFiniteGasAcousticTransportStep {
        transport,
        maximum_species_residual_mol: 0.0,
    };
    step.maximum_species_residual_mol = species_residual(
        parcels,
        &next_parcels,
        reactions,
        |id| step.reaction(id),
        cancelled,
    )?;
    if step.maximum_species_residual_mol > controls.species_tolerance_mol {
        return Err(Error::SpeciesBalance {
            maximum_residual_mol: step.maximum_species_residual_mol,
            tolerance_mol: controls.species_tolerance_mol,
        });
    }
    poll(cancelled)?;
    Ok(ReactiveFiniteGasAcousticTransportCandidate {
        parcels: next_parcels,
        nodes: next_nodes,
        step,
    })
}

/// Atomically publish every gas and original physical owner after the complete
/// reactive proposal passes. A later-node refusal or final cancellation leaves
/// all accepted parcels, sources, histories and clocks untouched.
#[allow(clippy::too_many_arguments)]
pub fn advance_reacting_finite_gas_acoustic_transport_with_paths<C: FnMut() -> bool>(
    parcels: &mut BTreeMap<u64, ClosedGasParcel>,
    nodes: &mut BTreeMap<u64, FiniteGasAcousticNode>,
    inputs: &BTreeMap<u64, FiniteGasAcousticInput>,
    reactions: &BTreeMap<u64, FiniteGasAcousticReaction<'_>>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    controls: ReactiveFiniteGasAcousticTransportSolve,
    cancelled: &mut C,
) -> Result<ReactiveFiniteGasAcousticTransportStep, Error> {
    let candidate = finite_gas_reacting_acoustic_transport_candidate_with_paths(
        parcels, nodes, inputs, reactions, paths, duration_s, controls, cancelled,
    )?;
    poll(cancelled)?;
    let (next_parcels, next_nodes, step) = candidate.into_parts();
    *parcels = next_parcels;
    *nodes = next_nodes;
    Ok(step)
}
