//! Source and reservoir entropy on original branch and compact-volume storage.

use super::MixtureNetwork1d;
use crate::multicomponent::dissipation::entropy::EntropyTracker;
use crate::multicomponent::{
    Boundaries, Boundary, Cell, Error, ExteriorMixtureEntropyReport, MixtureEntropyControl,
    MixtureEntropyReport, StepControl, cell, poll,
};
use fs_qty::Time;

impl MixtureNetwork1d {
    /// Advance the existing closed network with one extensive entropy gate.
    ///
    /// Every branch cell and finite plenum contributes its original-mixture
    /// thermal and composition entropy in one common source-reference basis.
    /// Internal characteristic-port transfers are not exterior entropy supplies:
    /// their two actual storage owners are already inside the accounting boundary.
    /// Each individual species continuity residual is bounded separately before
    /// removing only its numerical source-reference contribution.
    ///
    /// The native SSP RK2 transport, port closure, CFL, reconstruction and energy
    /// balance remain authoritative. A candidate whose entropy decreases beyond
    /// the declared allowance retries in the same bounded loop, with no state,
    /// temperature, energy or species correction. The returned duration may be
    /// shorter than requested. Failure or cancellation publishes no field,
    /// clock, plenum state or accepted exchange/report changes.
    ///
    /// Unconnected ends must be stationary reflecting walls; isolated periodic
    /// branches are allowed. Open reservoirs, molecular/chemical sources and
    /// moving branch walls require additional entropy supply and are excluded.
    /// This is measured finite-step admission of the engineering flux, not a
    /// newly proved entropy-stable flux or a global invariant-domain theorem.
    /// A negative semidiscrete entropy rate may still exhaust the retry budget.
    /// The entropy allowance applies to this accepted step; callers accumulating
    /// many steps must retain their own complete-interval budget.
    pub fn step_with_entropy(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        entropy: MixtureEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(super::NetworkStepReport, MixtureEntropyReport), Error> {
        let mut account = self.entropy_account(entropy, None, &mut cancelled)?;
        let report = self.step_impl(
            maximum_duration,
            controls,
            None,
            None,
            Some(&mut account),
            cancelled,
        )?;
        Ok((report, account.report()))
    }

    /// Advance prescribed chamber motion with the same closed-network entropy
    /// account. One signed dV/dt is required per plenum, as in
    /// [`Self::step_with_volume_rates`]. Actual endpoint volumes enter extensive
    /// entropy and species inventories; compression/expansion pressure work is
    /// mechanical energy exchange and supplies no entropy. The numerical
    /// irreversibility of the accepted time step remains visible in the report.
    /// Zero rates execute the rigid entropy path with identical arithmetic.
    pub fn step_with_volume_rates_and_entropy(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        volume_rates_m3_s: &[f64],
        entropy: MixtureEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(super::NetworkStepReport, MixtureEntropyReport), Error> {
        let moving = self.admit_volume_rates(volume_rates_m3_s, &mut cancelled)?;
        let mut account = self.entropy_account(entropy, None, &mut cancelled)?;
        let report = self.step_impl(
            maximum_duration,
            controls,
            moving.then_some(volume_rates_m3_s),
            None,
            Some(&mut account),
            cancelled,
        )?;
        Ok((report, account.report()))
    }

    /// Advance the native pressure-reservoir network with complete entropy
    /// admission over all branch cells and finite plenums.
    ///
    /// Reservoir placement, stationarity, source compatibility and characteristic
    /// fluxes are those of [`Self::step_with_reservoirs`]. Every unconnected end
    /// without a reservoir must be a stationary reflecting wall. Internal finite
    /// plenums and periodic faces stay inside the accounting boundary.
    ///
    /// The original exterior entropy variables contract each actual signed SSP
    /// species, momentum and referenced-energy transfer, including the advective
    /// entropy potential. For these stationary reservoirs this is the reservoir's
    /// heat/chemical effort on its total enthalpy exchange. Pressure-flow work is
    /// already in transported enthalpy; no additional p*dV or heat is added.
    /// Outflow uses the prescribed receiving reservoir's thermodynamics, so the
    /// measured generation includes irreversibility of that reservoir coupling.
    /// A nonzero species transfer at zero reservoir activity refuses at either
    /// stage, even if opposite stage transfers would cancel.
    ///
    /// The returned boundary entries sum all left and all right physical ends,
    /// respectively, in retained branch order. Every species residual is checked
    /// against original inventory plus absolute throughput at BOTH stages and
    /// all external ends. Only the source-reference effect of that numerical
    /// residual is removed from the single global entropy gate; actual entropy
    /// transported by species remains in the state and boundary reports.
    ///
    /// Rejected attempts restart the entire native network. No partial fields,
    /// clocks, geometry or accepted reports publish on refusal/cancellation.
    /// All-None sources execute the closed entropy path. This is measured
    /// finite-step admission of the existing engineering flux, not a fully
    /// discrete entropy or invariant-domain theorem. The allowance is per call.
    pub fn step_with_reservoirs_and_entropy(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        reservoirs: &[[Option<Cell<'_>>; 2]],
        entropy: MixtureEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(super::NetworkStepReport, ExteriorMixtureEntropyReport), Error> {
        let active = self.admit_reservoirs(reservoirs, &mut cancelled)?;
        let reservoirs = active.then_some(reservoirs);
        let mut account = self.entropy_account(entropy, reservoirs, &mut cancelled)?;
        let report = self.step_impl(
            maximum_duration,
            controls,
            None,
            reservoirs,
            Some(&mut account),
            cancelled,
        )?;
        Ok((report, account.exterior_report()))
    }

    /// Combine the same prescribed-reservoir entropy account with simultaneous
    /// finite-plenum motion. Native -p*dV pressure work changes extensive energy
    /// and the actual caloric state, but supplies no entropy. Accepted endpoint
    /// volumes enter the state account; transported exterior enthalpy remains
    /// separate from the network report's mechanical boundary work.
    ///
    /// All-None sources retain [`Self::step_with_volume_rates_and_entropy`];
    /// all-zero rates retain [`Self::step_with_reservoirs_and_entropy`].
    pub fn step_with_reservoirs_and_volume_rates_and_entropy(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        reservoirs: &[[Option<Cell<'_>>; 2]],
        volume_rates_m3_s: &[f64],
        entropy: MixtureEntropyControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(super::NetworkStepReport, ExteriorMixtureEntropyReport), Error> {
        let active = self.admit_reservoirs(reservoirs, &mut cancelled)?;
        let moving = self.admit_volume_rates(volume_rates_m3_s, &mut cancelled)?;
        let reservoirs = active.then_some(reservoirs);
        let mut account = self.entropy_account(entropy, reservoirs, &mut cancelled)?;
        let report = self.step_impl(
            maximum_duration,
            controls,
            moving.then_some(volume_rates_m3_s),
            reservoirs,
            Some(&mut account),
            cancelled,
        )?;
        Ok((report, account.exterior_report()))
    }

    /// Topology, not the presence of a reservoir, selects the physical boundary.
    /// Unowned exterior walls must still carry no species or energy at either
    /// SSP stage; internal ports are excluded even when their transfers are large.
    pub(super) fn entropy_exterior_ends(&self, branch: usize) -> [bool; 2] {
        if self.branches[branch].boundaries == Boundaries::Periodic {
            [false; 2]
        } else {
            self.endpoint_connections[branch].map(|connection| connection.is_none())
        }
    }

    pub(super) fn entropy_account(
        &self,
        controls: MixtureEntropyControl,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EntropyTracker, Error> {
        poll(cancelled)?;
        for (index, branch) in self.branches.iter().enumerate() {
            poll(cancelled)?;
            if let Boundaries::Ends { left, right } = branch.boundaries {
                for (end, boundary) in [left, right].into_iter().enumerate() {
                    if self.endpoint_connections[index][end].is_none()
                        && reservoirs.is_none_or(|sources| sources[index][end].is_none())
                        && boundary != Boundary::Reflecting
                    {
                        return Err(Error::EntropyDomain(
                            "network entropy requires a reservoir or stationary reflecting exterior end",
                        ));
                    }
                }
            }
        }
        let gas = self.gas();
        let n = gas.species_count();
        let m = n + 2;
        let count = self.plenums.len()
            + self
                .branches
                .iter()
                .map(|branch| branch.cell_count())
                .sum::<usize>();
        let cells = self
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
                self.plenums
                    .iter()
                    .map(|plenum| (plenum.state(), plenum.volume_m3())),
            );
        let mut account = EntropyTracker::from_cells(gas, count, controls, None, cells, cancelled)?;
        if let Some(reservoirs) = reservoirs {
            let count = (0..self.branches.len())
                .filter(|&index| self.entropy_exterior_ends(index).into_iter().any(|end| end))
                .count();
            let exterior = reservoirs.iter().enumerate().filter_map(|(index, pair)| {
                self.entropy_exterior_ends(index)
                    .into_iter()
                    .any(|end| end)
                    .then(|| {
                        core::array::from_fn(|end| {
                            pair[end]
                                .map(|state| (state, self.work.reservoir_primitive[index][end]))
                        })
                    })
            });
            account.prepare_exterior_pairs(gas, count, exterior, cancelled)?;
        }
        Ok(account)
    }
}
