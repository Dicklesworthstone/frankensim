//! Closed-network entropy admission on original branch and compact-volume storage.

use super::MixtureNetwork1d;
use crate::multicomponent::dissipation::entropy::EntropyTracker;
use crate::multicomponent::{
    Boundaries, Boundary, Error, MixtureEntropyControl, MixtureEntropyReport, StepControl, cell,
    poll,
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
        let mut account = self.closed_entropy_account(entropy, &mut cancelled)?;
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
        let mut account = self.closed_entropy_account(entropy, &mut cancelled)?;
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

    fn closed_entropy_account(
        &self,
        controls: MixtureEntropyControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EntropyTracker, Error> {
        poll(cancelled)?;
        for (index, branch) in self.branches.iter().enumerate() {
            poll(cancelled)?;
            if let Boundaries::Ends { left, right } = branch.boundaries {
                for (end, boundary) in [left, right].into_iter().enumerate() {
                    if self.endpoint_connections[index][end].is_none()
                        && boundary != Boundary::Reflecting
                    {
                        return Err(Error::EntropyDomain(
                            "closed network entropy requires stationary reflecting exterior ends",
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
        EntropyTracker::from_cells(gas, count, controls, None, cells, cancelled)
    }
}
