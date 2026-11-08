//! Original-source entropy over every finite molecular storage owner.

use super::{
    AdaptiveNetworkDissipationControl, AdaptiveNetworkDissipationReport, DissipationControl,
    DissipationError, MixtureNetwork1d, NetworkDissipationOperator, NetworkDissipationReport,
    Time, TransportLaw,
};
use crate::multicomponent::cell;
use crate::multicomponent::dissipation::entropy::EntropyTracker;
use crate::multicomponent::dissipation::{MolecularEntropyControl, MolecularEntropyReport};

mod stabilized;

impl<L: TransportLaw> NetworkDissipationOperator<L> {
    /// Advance a closed molecular network to declared local accuracy with
    /// source-conjugate diffusion and finite-interval entropy admission.
    ///
    /// The original adaptive RKL2 full/two-half estimator controls every species,
    /// temperature, velocity and signed port transfer. Every probe uses the same
    /// source-conjugate Maxwell--Stefan, divided enthalpy, Fourier, reciprocal
    /// cross-effect and viscous-work laws at branch interiors and finite ports.
    ///
    /// Each accepted paired-half endpoint and the direct original-to-final
    /// interval must satisfy the original entropy and species-continuity bounds.
    /// The discarded full-step estimate and internal half/stage states are not
    /// separately admitted physical entropy intervals. Entropy refusals halve
    /// the complete probe from its unchanged anchor; temporal refusals retain
    /// the native local-error controller. All probes share one actual RHS budget.
    /// Native endpoint temperatures supply entropy without another EOS inverse,
    /// constitutive RHS or per-probe allocation.
    ///
    /// Geometry is stationary, closed and constant-area. Every active face
    /// requires common species support. Baths, moving volumes, reactions and
    /// advection are outside this entry. Measured entropy and estimated molecular
    /// local accuracy do not certify a nonlinear invariant domain or global
    /// temporal error. Refusal preserves all accepted fields, clocks and reports.
    pub fn advance_adaptive_entropy_consistent_for(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: AdaptiveNetworkDissipationControl,
        entropy: MolecularEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(AdaptiveNetworkDissipationReport, MolecularEntropyReport), DissipationError>
    {
        controls.validate()?;
        self.admit_plenum_volumes(network, &mut cancelled)?;
        let mut account = initial_account(network, entropy, &mut cancelled)?;
        self.trial_report = None;
        let report = super::adaptive::advance::<true, L>(
            self,
            network,
            duration,
            controls,
            None,
            &mut 0,
            Some(&mut account),
            &mut cancelled,
        )?;
        self.trial_report = Some(report.stabilized.transport);
        self.publish_trial_report();
        Ok((report, account.report()))
    }

    /// Advance source-conjugate molecular transport through finite plenums.
    ///
    /// Interior branch faces and branch/plenum ports use the SAME original
    /// divided caloric enthalpies, logarithmic chemical efforts, Maxwell--Stefan
    /// law, Fourier heat, authored reciprocal Soret/Dufour terms and arithmetic
    /// velocity in viscous work as the axial entropy-consistent owner. Port
    /// gradients retain their authored plenum mixing lengths. The stationary
    /// plenum support receives axial impulse but does no mechanical work.
    /// No separate mixing or viscous heat source is added.
    ///
    /// Every actual synchronized SSPRK2 substep and the direct original-to-final
    /// interval must pass one source-entropy and individual-species continuity
    /// allowance. Internal ports are inside this accounting boundary: both
    /// participating inventories contribute and no internal entropy supply is
    /// subtracted. Only the source-reference effect of numerical species drift
    /// is removed from admission; the complete source entropy remains reported.
    ///
    /// Species must be positive on both sides of each active face or absent on
    /// both. One-sided absence refuses without a concentration floor or another
    /// flux. Exterior ends retain stationary zero molecular heat/species flux;
    /// periodic branches and isolated finite plenums are admitted. This entry
    /// does not compose baths, chemistry, advection, moving volumes or RKL2.
    ///
    /// The flow clock is unchanged. Any refusal or cancellation preserves all
    /// network fields, geometry and accepted molecular/advection receipts.
    /// This is measured finite-step admission of the original ideal-mixture
    /// model, not a fully discrete entropy or invariant-domain theorem.
    pub fn advance_entropy_consistent_for(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: DissipationControl,
        entropy: MolecularEntropyControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(NetworkDissipationReport, MolecularEntropyReport), DissipationError> {
        self.advance_entropy::<true>(network, duration, controls, entropy, cancelled)
    }

    /// Measure the original engineering molecular faces against the same
    /// complete finite-network entropy account. When admission succeeds, native
    /// states and transfer receipts are unchanged. A negative spatial entropy
    /// rate can exhaust step halving; reducing time does not repair an incorrect
    /// finite-temperature caloric/chemical pairing. Use
    /// [`Self::advance_entropy_consistent_for`] for source-conjugate faces.
    pub fn advance_for_with_entropy(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: DissipationControl,
        entropy: MolecularEntropyControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(NetworkDissipationReport, MolecularEntropyReport), DissipationError> {
        self.advance_entropy::<false>(network, duration, controls, entropy, cancelled)
    }

    fn advance_entropy<const ENTROPY_FACES: bool>(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: DissipationControl,
        entropy: MolecularEntropyControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(NetworkDissipationReport, MolecularEntropyReport), DissipationError> {
        let report = self.advance_entropy_for_trial::<ENTROPY_FACES>(
            network, duration, controls, entropy, cancelled,
        )?;
        self.publish_trial_report();
        Ok(report)
    }

    /// A complete composing owner publishes these molecular receipts only after
    /// its own final state, source and entropy checks succeed.
    pub(in crate::multicomponent::network) fn advance_entropy_for_trial<
        const ENTROPY_FACES: bool,
    >(
        &mut self,
        network: &mut MixtureNetwork1d,
        duration: Time,
        controls: DissipationControl,
        entropy: MolecularEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(NetworkDissipationReport, MolecularEntropyReport), DissipationError> {
        controls.validate()?;
        self.admit_plenum_volumes(network, &mut cancelled)?;
        let mut account = initial_account(network, entropy, &mut cancelled)?;
        self.trial_report = None;
        let report = super::evolve::advance_with_entropy::<ENTROPY_FACES, L>(
            self,
            network,
            duration,
            controls,
            &mut account,
            &mut cancelled,
        )?;
        self.trial_report = Some(report);
        Ok((report, account.report()))
    }
}

/// One original storage basis for either native molecular time integrator.
fn initial_account(
    network: &MixtureNetwork1d,
    controls: MolecularEntropyControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<EntropyTracker, DissipationError> {
    let gas = network.gas();
    let n = gas.species_count();
    let m = n + 2;
    let count = network
        .branches
        .iter()
        .map(|branch| branch.cell_count())
        .sum::<usize>()
        + network.plenums.len();
    let cells = network
        .branches
        .iter()
        .flat_map(|branch| {
            branch
                .states
                .chunks_exact(m)
                .zip(branch.mesh.volumes())
                .map(move |(q, &volume)| (cell(q, n), volume))
        })
        .chain(
            network
                .plenums
                .iter()
                .map(|plenum| (plenum.state(), plenum.volume_m3())),
        );
    Ok(EntropyTracker::from_cells(
        gas, count, controls, None, cells, cancelled,
    )?)
}

/// Reuse the actual candidate inversions, including kinetic-to-internal
/// conversion. The closed account has no exterior flux to contract.
pub(super) fn collect_candidate<L: TransportLaw>(
    owner: &NetworkDissipationOperator<L>,
    account: &mut EntropyTracker,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DissipationError> {
    let m = owner.gas.species_count() + 2;
    let volumes = owner
        .branches
        .iter()
        .flat_map(|shape| shape.mesh.volumes().iter())
        .chain(owner.plenum_volumes.iter());
    for (index, ((old, new), &volume)) in owner
        .work
        .current
        .chunks_exact(m)
        .zip(owner.work.candidate.chunks_exact(m))
        .zip(volumes)
        .enumerate()
    {
        super::poll(cancelled)?;
        account.candidate_cell(
            &owner.gas,
            index,
            old,
            new,
            owner.work.primitive[index].temperature,
            volume,
            cancelled,
        )?;
    }
    Ok(())
}

/// The caller's fields are still original. Reduce directly before the native
/// finish can publish any branch, plenum, port or molecular report.
pub(super) fn finish<L: TransportLaw>(
    owner: &NetworkDissipationOperator<L>,
    network: &MixtureNetwork1d,
    account: &mut EntropyTracker,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DissipationError> {
    let n = owner.gas.species_count();
    let m = n + 2;
    let original = network
        .branches
        .iter()
        .flat_map(|branch| {
            branch
                .states
                .chunks_exact(m)
                .zip(branch.mesh.volumes())
                .map(move |(q, &volume)| (cell(q, n), volume))
        })
        .chain(
            network
                .plenums
                .iter()
                .map(|plenum| (plenum.state(), plenum.volume_m3())),
        );
    let cells = original
        .zip(owner.work.current.chunks_exact(m))
        .map(|((old, volume), new)| (old, cell(new, n), volume));
    account.finish_cells(&owner.gas, cells, cancelled)?;
    Ok(())
}
