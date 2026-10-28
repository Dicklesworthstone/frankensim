//! Source-bound reactions in finite chambers driven by prescribed volumes.
//!
//! Original owners advance one interval in flow -> chemistry -> motion order.
//! Every restriction sees the same accepted network. Incoming species then
//! react in the actual mixed parcel, preserving total formation internal
//! energy. The resulting composition, temperature and pressure drive the
//! original midpoint pressure/volume work. No reaction enthalpy is added as
//! external heat; products become available to outflow on the next interval.
//!
//! This is a first-order operator split even when the chemical owner uses an
//! implicit or adaptive policy. Refining chemistry alone does not refine the
//! flow/motion splitting. It is not a resolved flame, ignition-time certificate,
//! a simultaneous chemical/mechanical endpoint or a measured engine model.

use super::{
    BTreeMap, ClosedGasParcel, Error, GasEnergySolve, PrescribedGasCandidate, PrescribedGasDrive,
    PrescribedGasSolve, PrescribedGasStep, SnapshotNozzleLink, Sum, poll, positive,
};
use crate::cavity_phs::piston::reactive::{
    GasReaction, GasReactionError, GasReactionSolve, GasReactionStep,
};
use fs_thermochem::SpeciesId;

/// Bind an existing complete mechanism to one graph volume. Unlisted volumes
/// remain frozen. The immutable catalog must cover every species that can flow
/// into that volume, with the same caloric, transport and source-pack identities.
/// Explicit, implicit and chemical-subcycling policies remain with this owner.
#[derive(Clone, Copy, Debug)]
pub struct DrivenGasReaction<'a> {
    pub chamber: u64,
    pub reaction: &'a GasReaction,
    pub controls: GasReactionSolve,
}

/// Existing flow/motion controls plus an independent whole-network molar gate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactivePrescribedGasSolve {
    pub network: PrescribedGasSolve,
    /// Maximum |sum(n_end - n_start - nu*xi)| for EACH species [mol].
    /// This encloses flow, all chemical channels and volume work; individual
    /// chemical and transport owners retain their original admission bounds.
    pub species_tolerance_mol: f64,
}

/// The complete proposal, including source-resolved chemical evaluations.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactivePrescribedGasCandidate {
    motion: PrescribedGasCandidate,
    reactions: BTreeMap<u64, GasReactionStep>,
    maximum_species_residual_mol: f64,
}
impl ReactivePrescribedGasCandidate {
    /// All finite gas stores, with reaction and pressure work already applied.
    #[must_use]
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        self.motion.parcels()
    }
    /// Source evaluations, extents and original chemical admission reports.
    #[must_use]
    pub fn reactions(&self) -> &BTreeMap<u64, GasReactionStep> {
        &self.reactions
    }
    /// The original report now measures accepted-input to final storage across
    /// ALL three stages. Local heat/volume reports start after chemistry.
    #[must_use]
    pub const fn step(&self) -> &PrescribedGasStep {
        self.motion.step()
    }
    /// Measured maximum complete-interval stoichiometric residual [mol].
    #[must_use]
    pub const fn maximum_species_residual_mol(&self) -> f64 {
        self.maximum_species_residual_mol
    }
    /// Accept all endpoints together. Extents and heat/work are already applied.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, ClosedGasParcel>,
        PrescribedGasStep,
        BTreeMap<u64, GasReactionStep>,
        f64,
    ) {
        let (parcels, step) = self.motion.into_parts();
        (
            parcels,
            step,
            self.reactions,
            self.maximum_species_residual_mol,
        )
    }
}

/// Couple declared homogeneous chemistry to moving OR held finite volumes.
///
/// One physical duration drives original nozzle mass/enthalpy transfer, each
/// selected reaction, and the prescribed volume/heat endpoint. All reactions
/// consume their own post-flow inventory once; mechanical root trials never
/// consume another amount. Source domains, exact reaction conservation,
/// detailed balance, consumption/entropy checks, midpoint work and complete
/// energy/mass/species gates remain active. Empty bindings reproduce the
/// original nonreacting path, with an additional independent species audit.
///
/// A failed later chamber, inventory allowance, energy/species gate or
/// cancellation leaves EVERY input untouched. All proposed owners must be
/// published together with their mechanical drive's reciprocal work debit.
/// Neither formation energy nor transported enthalpy is a second heat supply.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    links: &[SnapshotNozzleLink],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    drives: &BTreeMap<u64, PrescribedGasDrive>,
    reactions: &[DrivenGasReaction<'_>],
    controls: ReactivePrescribedGasSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ReactivePrescribedGasCandidate, Error> {
    poll(&mut cancelled)?;
    if parcels.is_empty()
        || parcels.len() > super::MAX_SNAPSHOT_PARCELS
        || !positive(controls.species_tolerance_mol)
        || reactions.len() > parcels.len()
    {
        return Err(Error::Invalid(
            "positive species tolerance and at most one reaction per gas volume required",
        ));
    }
    let mut ordered = BTreeMap::new();
    for binding in reactions {
        poll(&mut cancelled)?;
        if !parcels.contains_key(&binding.chamber)
            || ordered.insert(binding.chamber, binding).is_some()
        {
            return Err(Error::Invalid(
                "unique declared reactions for existing gas volumes required",
            ));
        }
    }
    let (motion, reactions) = super::staged_candidate(
        parcels,
        links,
        duration_s,
        endpoints,
        drives,
        controls.network,
        |gases, cancelled| {
            let mut reports = BTreeMap::new();
            for (&id, binding) in &ordered {
                poll(cancelled)?;
                let before = gases
                    .get(&id)
                    .ok_or(Error::Invalid("original flow omitted a reacting volume"))?;
                let candidate = binding
                    .reaction
                    .parcel_candidate(before, duration_s, binding.controls, &mut *cancelled)
                    .map_err(|error| match error {
                        GasReactionError::Cancelled => Error::Cancelled,
                        error => Error::Reaction { id, error },
                    })?;
                if candidate.gas().mass_kg().to_bits() != before.mass_kg().to_bits()
                    || candidate.gas().volume_m3().to_bits() != before.volume_m3().to_bits()
                {
                    return Err(Error::Invalid(
                        "closed chemistry changed the transported mass or current volume",
                    ));
                }
                let (gas, report) = candidate.into_parts();
                gases.insert(id, gas);
                reports.insert(id, report);
            }
            Ok(reports)
        },
        &mut cancelled,
    )?;
    let mut residuals: BTreeMap<SpeciesId, Sum> = BTreeMap::new();
    for (initial, sign) in [(parcels, -1.0), (motion.parcels(), 1.0)] {
        for gas in initial.values() {
            poll(&mut cancelled)?;
            let model = gas.model().thermodynamics();
            let total = gas.mass_kg()
                / gas
                    .state()
                    .thermodynamics()
                    .properties()
                    .molar_mass()
                    .value();
            if !positive(total) {
                return Err(Error::Invalid(
                    "whole-network gas amount is unrepresentable",
                ));
            }
            for (species, fraction) in model
                .components()
                .iter()
                .zip(model.mole_composition().fractions())
            {
                poll(&mut cancelled)?;
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
    for reaction in reactions.values() {
        for species in &reaction.species {
            poll(&mut cancelled)?;
            residuals
                .entry(species.species.clone())
                .or_default()
                .add(-species.stoichiometric_change_mol)?;
        }
    }
    let maximum_species_residual_mol = residuals.values().fold(0.0_f64, |maximum, residual| {
        maximum.max(residual.value.abs())
    });
    if maximum_species_residual_mol > controls.species_tolerance_mol {
        return Err(Error::SpeciesBalance {
            maximum_residual_mol: maximum_species_residual_mol,
        });
    }
    poll(&mut cancelled)?;
    Ok(ReactivePrescribedGasCandidate {
        motion,
        reactions,
        maximum_species_residual_mol,
    })
}
