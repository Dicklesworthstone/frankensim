//! Concurrent reaction channels on one original species/energy inventory.
//!
//! Every channel evaluates its existing activity-based forward/reverse law at
//! the SAME accepted state. Sum nu_ir*dt*V*r_r before the original gas owner
//! reconstructs temperature/composition once. A channel cannot consume species
//! just produced by another channel in this explicit interval. Gross consumed
//! amounts across signed net channel extents are bounded by the original stock.
//! This is the default first-order explicit path. The optional
//! `implicit::coupled` path evaluates all channels at one solved endpoint.

use super::*;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_thermochem::ReactionId;
use material_rate::MaterialRateScale;

/// Hard channel bound per homogeneous mechanism, additional to species bounds.
pub const MAX_COUPLED_REACTIONS: usize = 32;

/// One admitted channel without another copy of the shared species catalog.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct Channel {
    kinetics: DeclaredStoichiometricMassActionV1,
    rate_material: Option<MaterialRateScale>,
}
impl Channel {
    fn id(&self) -> &ReactionId { self.kinetics.reverse_rate_closure().equilibrium().reaction() }
}

/// Original kinetic evaluation and extent for one additional reaction channel.
/// Extents from different stoichiometric columns must NOT be summed as if they
/// described one reaction. The enclosing species changes are already combined.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionChannelStep {
    pub reaction_id: ReactionId,
    pub rate: DeclaredStoichiometricMassActionEvaluationV1,
    pub rate_material: Option<ResolvedMaterialStatePoint>,
    pub extent_mol: f64,
}

impl GasReaction {
    /// Combine already-admitted reactions over ONE identical complete catalog.
    /// Channels must select columns of the SAME certified stoichiometric matrix.
    /// Each retains its own reverse-equilibrium closure
    /// and optional material rate. Bind rate sources BEFORE combining channels.
    /// Nested combinations are flattened; no inventories or time histories are
    /// copied into the mechanism. Bind a shared subcycling policy AFTER combining.
    /// Duplicate reaction IDs and source conflicts
    /// refuse. Canonical reaction-ID order removes caller-order dependence.
    ///
    /// A singleton preserves the original object and numerical path. For more
    /// than one channel, legacy `kinetics`, `rate_material` and the report's
    /// `rate`/`extent_mol` refer to the first canonical channel. Use `channels`
    /// and `reaction_extents` to observe the whole mechanism. All existing
    /// parcel, piston, network and thermal/electrical entry points accept it.
    pub fn combine(reactions: Vec<Self>) -> Result<Self, GasReactionError> {
        if reactions.is_empty() || reactions.len() > MAX_COUPLED_REACTIONS {
            return Err(invalid("a coupled mechanism needs a bounded nonempty reaction list"));
        }
        let count = reactions.iter().try_fold(0usize, |n, r|
            n.checked_add(1 + r.additional_channels.len()))
            .ok_or_else(|| invalid("reaction mechanism size overflow"))?;
        if count == 0 || count > MAX_COUPLED_REACTIONS {
            return Err(invalid("a coupled mechanism needs 1..=32 reaction channels"));
        }
        if reactions.len() == 1 {
            return Ok(reactions.into_iter().next().expect("nonempty admitted mechanism"));
        }
        if reactions.iter().any(|r| r.subcycling.is_some() || r.implicit.is_some() || r.implicit_mechanism.is_some()) {
            return Err(invalid("combine channels before binding one shared integration policy"));
        }
        let catalog = &reactions[0].catalog;
        if reactions.iter().any(|r| &r.catalog != catalog) {
            return Err(invalid("reaction channels require one identical complete source catalog"));
        }
        let authority = reactions[0].kinetics.reverse_rate_closure().equilibrium();
        if reactions.iter().any(|r| {
            let other = r.kinetics.reverse_rate_closure().equilibrium();
            other.stoichiometric_matrix().identity() != authority.stoichiometric_matrix().identity()
                || other.conservation_certificate().identity() != authority.conservation_certificate().identity()
        }) {
            return Err(invalid("coupled channels require one stoichiometric matrix and conservation certificate"));
        }
        let mut inputs = reactions.into_iter();
        let first = inputs.next().expect("nonempty admitted mechanism");
        if count == 1 { return Ok(first); }
        let GasReaction { catalog, kinetics, rate_material, additional_channels, subcycling: _, implicit: _, implicit_mechanism: _ } = first;
        let mut channels = Vec::with_capacity(count);
        channels.push(Channel { kinetics, rate_material });
        channels.extend(additional_channels);
        for reaction in inputs {
            channels.push(Channel { kinetics: reaction.kinetics, rate_material: reaction.rate_material });
            channels.extend(reaction.additional_channels);
        }
        channels.sort_by(|a, b| a.id().cmp(b.id()));
        if channels.windows(2).any(|pair| pair[0].id() == pair[1].id()) {
            return Err(invalid("reaction mechanism contains a duplicate reaction identity"));
        }
        let mut channels = channels.into_iter();
        let primary = channels.next().expect("admitted channel count");
        Ok(Self { catalog, kinetics: primary.kinetics, rate_material: primary.rate_material,
            additional_channels: channels.collect(), subcycling: None, implicit: None, implicit_mechanism: None })
    }

    /// Original immutable kinetic laws and source bindings, in canonical order.
    #[must_use]
    pub fn channels(&self) -> impl Iterator<Item = (&DeclaredStoichiometricMassActionV1, Option<&MaterialRateScale>)> {
        std::iter::once((&self.kinetics, self.rate_material.as_ref()))
            .chain(self.additional_channels.iter().map(|c| (&c.kinetics, c.rate_material.as_ref())))
    }
}

impl GasReactionStep {
    /// The original reaction identities and distinct extents, primary first.
    #[must_use]
    pub fn reaction_extents(&self) -> impl Iterator<Item = (&ReactionId, f64)> {
        std::iter::once((&self.primary_reaction, self.extent_mol))
            .chain(self.additional_reactions.iter().map(|r| (&r.reaction_id, r.extent_mol)))
    }
    /// True even when the primary channel is idle or species changes cancel.
    #[must_use]
    pub fn has_progress(&self) -> bool {
        self.reaction_extents().any(|(_, extent)| extent != 0.0)
            || self.subcycling.as_ref().is_some_and(|r|
                r.steps.iter().any(|s| s.reaction.has_progress()))
    }
}

pub(super) struct Progress {
    pub primary: ReactionChannelStep,
    pub additional: Vec<ReactionChannelStep>,
    pub changes: Vec<f64>,
    pub expected: Vec<f64>,
    pub consumed: f64,
}

#[derive(Clone, Copy, Default)]
pub(crate) struct Sum { pub(crate) value: f64, correction: f64 }
impl Sum {
    pub(crate) fn add(&mut self, term: f64) -> Result<(), GasReactionError> {
        let adjusted = term - self.correction;
        let next = self.value + adjusted;
        let correction = (next - self.value) - adjusted;
        if ![term, adjusted, next, correction].iter().all(|x| x.is_finite()) {
            return Err(invalid("reaction mechanism species accumulation is unrepresentable"));
        }
        self.value = next; self.correction = correction;
        Ok(())
    }
}

impl GasReaction {
    /// Shared canonical activities and material-rate query for explicit initial
    /// states and implicit endpoint trials. No caller supplies partial pressures.
    fn channel_activities(
        &self, kinetics: &DeclaredStoichiometricMassActionV1,
        gas: &ClosedGasParcel, amounts: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<IdealGasSpeciesActivityV1>, GasReactionError> {
        let components = self.catalog.thermodynamics().components();
        let total: f64 = amounts.iter().sum();
        if !positive(total) { return Err(invalid("reaction total amount is unrepresentable")); }
        let equilibrium = kinetics.reverse_rate_closure().equilibrium();
        let pressure_ratio = gas.state().state().pressure / equilibrium.convention().reference_pressure().value();
        let mut activities = Vec::with_capacity(equilibrium.terms().len());
        for term in equilibrium.terms() {
            poll(cancelled)?;
            let index = components.binary_search_by(|s| s.species().cmp(term.species()))
                .map_err(|_| invalid("bound reaction species disappeared"))?;
            let value = (amounts[index] / total) * pressure_ratio;
            if amounts[index] > 0.0 && !positive(value) {
                return Err(invalid("positive species partial pressure is unrepresentable"));
            }
            activities.push(IdealGasSpeciesActivityV1::new(term.species().clone(), Dimensionless::new(value))?);
        }
        poll(cancelled)?;
        Ok(activities)
    }

    pub(super) fn evaluate_channel(
        &self, kinetics: &DeclaredStoichiometricMassActionV1,
        source: Option<&MaterialRateScale>, gas: &ClosedGasParcel, amounts: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(DeclaredStoichiometricMassActionEvaluationV1, Option<ResolvedMaterialStatePoint>), GasReactionError> {
        let activities = self.channel_activities(kinetics, gas, amounts, cancelled)?;
        material_rate::evaluate_rate_for(kinetics, source, gas, activities, cancelled)
    }

    pub(super) fn evaluate_channel_trial(
        &self, kinetics: &DeclaredStoichiometricMassActionV1,
        source: Option<&MaterialRateScale>, gas: &ClosedGasParcel, amounts: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<fs_thermochem::equilibrium::NetReactionProgressRateV1>, GasReactionError> {
        let activities = self.channel_activities(kinetics, gas, amounts, cancelled)?;
        material_rate::evaluate_rate_for_trial(kinetics, source, gas, activities, cancelled)
    }

    /// Shared single/multiple-channel rate and inventory operation. No gas
    /// state is reconstructed or modified until every channel has succeeded.
    pub(super) fn progress(
        &self, gas: &ClosedGasParcel, initial: &[f64], duration_s: f64,
        controls: GasReactionSolve, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Progress, GasReactionError> {
        let components = self.catalog.thermodynamics().components();
        let mut changes = vec![Sum::default(); components.len()];
        let mut withdrawn = vec![Sum::default(); components.len()];
        let mut reports = Vec::with_capacity(1 + self.additional_channels.len());
        for (kinetics, source) in self.channels() {
            poll(cancelled)?;
            let equilibrium = kinetics.reverse_rate_closure().equilibrium();
            let terms = equilibrium.terms();
            let (rate, rate_material) = self.evaluate_channel(kinetics, source, gas, initial, cancelled)?;
            let progress = rate.net_progress_rate().ok_or_else(||
                invalid("reaction rate is log-only; no finite extent can be integrated"))?.value();
            let extent = (progress * gas.volume_m3()) * duration_s;
            if !extent.is_finite() || (progress != 0.0 && extent == 0.0) {
                return Err(invalid("reaction extent overflowed or lost a nonzero rate"));
            }
            for term in terms {
                poll(cancelled)?;
                let index = components.binary_search_by(|s| s.species().cmp(term.species()))
                    .map_err(|_| invalid("bound reaction species disappeared"))?;
                let change = term.coefficient() as f64 * extent;
                changes[index].add(change)?;
                if change < 0.0 { withdrawn[index].add(-change)?; }
            }
            reports.push(ReactionChannelStep { reaction_id: equilibrium.reaction().clone(),
                rate, rate_material, extent_mol: extent });
        }
        let mut consumed: f64 = 0.0;
        let mut expected = Vec::with_capacity(initial.len());
        for (i, &amount) in initial.iter().enumerate() {
            poll(cancelled)?;
            // Production in a different channel cannot replenish this allowance.
            if withdrawn[i].value > 0.0 {
                let fraction = withdrawn[i].value / amount;
                if !fraction.is_finite() || fraction > controls.maximum_consumed_fraction {
                    return Err(GasReactionError::ConsumptionLimit { species: components[i].species().clone(),
                        fraction, maximum: controls.maximum_consumed_fraction });
                }
                consumed = consumed.max(fraction);
            }
            let next = amount + changes[i].value;
            if !next.is_finite() || next < 0.0 || (changes[i].value != 0.0 && next == amount) {
                return Err(invalid("reaction lost a species increment or made an invalid inventory"));
            }
            expected.push(next);
        }
        let mut reports = reports.into_iter();
        Ok(Progress { primary: reports.next().expect("every reaction has a primary channel"),
            additional: reports.collect(), changes: changes.into_iter().map(|x| x.value).collect(), expected, consumed })
    }
}


