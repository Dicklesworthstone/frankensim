//! One potential convention for resolved branches and compact, stationary nodes.

use super::{
    Accumulator, BranchEnd, Conserved, Error, Euler1d, GasNetwork1d, NetworkTotals, Plenum,
    PlenumConnection, poll,
};

/// Specific gravitational potentials [J/kg] in the retained topology order.
///
/// A hydrostatic branch uses `phi(x) = offset - g * (x - x_ref)` with its
/// original reference. A branch without gravity uses the constant `offset`.
/// Each compact plenum is at its declared single potential, not a hydrostatic
/// column: its incident faces must agree with that value. A common additive
/// shift changes the energy datum and cannot change the evolving gas state.
#[derive(Clone, Debug, PartialEq)]
pub struct NetworkGravity {
    /// One offset for each branch's original local zero [J/kg].
    pub branch_offsets_j_kg: Vec<f64>,
    /// One specific potential for each acoustically compact plenum [J/kg].
    pub plenum_potentials_j_kg: Vec<f64>,
}

/// Accepted energy balance in the declared gravitational potential convention.
/// Gas storage and gas boundary exchanges remain in `NetworkStepReport`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NetworkGravityReport {
    /// Common zero subtracted before potential arithmetic [J/kg]. The first
    /// branch offset supplies this conditioning reference, without changing physics.
    pub reference_potential_j_kg: f64,
    /// Initial midpoint/compact storage of `mass * (potential - reference)` [J].
    pub before_centered_potential_energy: f64,
    /// Accepted potential storage in the same centered coordinate [J].
    pub after_centered_potential_energy: f64,
    /// External `(potential - reference) * mass` transport, positive inward [J].
    pub boundary_centered_potential_energy_change: f64,
    /// Actual gravitational work on the resolved branch gas [J].
    pub gravity_work: f64,
    /// Gas-plus-centered-potential change minus its external transport [J].
    /// This conditioned residual retains small gravitational work even under
    /// a large arbitrary common datum. Internal port transfer cancels.
    pub conservation_residual: f64,
    /// Residual in the caller's absolute potential coordinate [J], evaluated as
    /// `conservation_residual + reference_potential_j_kg * mass_residual`.
    /// Absolute before/after storage similarly adds `reference * total_mass`;
    /// differencing those large levels directly would discard resolved work.
    pub absolute_conservation_residual: f64,
}

impl NetworkGravity {
    pub(super) fn reference(&self) -> f64 {
        self.branch_offsets_j_kg[0]
    }

    pub(super) fn branch_potential(&self, index: usize, branch: &Euler1d, x: f64) -> f64 {
        self.branch_offsets_j_kg[index]
            + branch.hydrostatic.as_ref().map_or(0.0, |b| b.potential(x))
    }

    pub(super) fn centered_branch_potential(&self, index: usize, branch: &Euler1d, x: f64) -> f64 {
        (self.branch_offsets_j_kg[index] - self.reference())
            + branch.hydrostatic.as_ref().map_or(0.0, |b| b.potential(x))
    }

    pub(super) fn admit(
        &self,
        branches: &[Euler1d],
        plenums: &[Plenum],
        connections: &[PlenumConnection],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        if self.branch_offsets_j_kg.len() != branches.len()
            || self.plenum_potentials_j_kg.len() != plenums.len()
            || !self
                .branch_offsets_j_kg
                .iter()
                .chain(&self.plenum_potentials_j_kg)
                .all(|x| x.is_finite())
        {
            return Err(Error::Invalid(
                "network potential count or non-finite datum",
            ));
        }
        if !self
            .branch_offsets_j_kg
            .iter()
            .chain(&self.plenum_potentials_j_kg)
            .all(|x| (x - self.reference()).is_finite())
        {
            return Err(Error::Invalid(
                "unrepresentable relative network potentials",
            ));
        }
        for (i, branch) in branches.iter().enumerate() {
            for &x in branch.mesh.positions() {
                poll(cancelled)?;
                if !self.branch_potential(i, branch, x).is_finite() {
                    return Err(Error::Invalid("unrepresentable network face potential"));
                }
            }
        }
        for connection in connections {
            poll(cancelled)?;
            let branch = &branches[connection.branch];
            let face = match connection.end {
                BranchEnd::Left => 0,
                BranchEnd::Right => branch.states.len(),
            };
            let local = branch
                .hydrostatic
                .as_ref()
                .map_or(0.0, |b| b.potential(branch.mesh.positions()[face]));
            // Subtract the shared datum BEFORE adding the local head. In
            // particular, an enormous common additive zero must not enlarge
            // the tolerance enough to admit a real height mismatch.
            let translation = self.branch_offsets_j_kg[connection.branch]
                - self.plenum_potentials_j_kg[connection.plenum];
            let mismatch = local + translation;
            let scale = local.abs().max(translation.abs()).max(f64::MIN_POSITIVE);
            if !translation.is_finite()
                || !mismatch.is_finite()
                || mismatch.abs() > 64.0 * f64::EPSILON * scale
            {
                return Err(Error::Invalid(
                    "connected face and compact plenum potentials differ",
                ));
            }
        }
        // Admission also rejects a finite datum whose extensive storage cannot
        // be represented. Future accepted stages perform the same check.
        let storage = self.storage(branches, plenums, None, cancelled)?;
        if !(storage.energy + self.reference() * storage.mass).is_finite() {
            return Err(Error::Invalid(
                "unrepresentable absolute network potential energy",
            ));
        }
        Ok(())
    }

    pub(super) fn storage(
        &self,
        branches: &[Euler1d],
        plenums: &[Plenum],
        candidate_plenums: Option<&[Conserved]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NetworkTotals, Error> {
        let mut sum = Accumulator::default();
        for (b, branch) in branches.iter().enumerate() {
            let states = if candidate_plenums.is_some() {
                &branch.work.candidate
            } else {
                &branch.states
            };
            let offset = self.branch_offsets_j_kg[b] - self.reference();
            let mut local = Accumulator::default();
            for (i, state) in states.iter().enumerate() {
                poll(cancelled)?;
                let potential = offset
                    + branch
                        .hydrostatic
                        .as_ref()
                        .map_or(0.0, |balance| balance.potential(branch.mesh.center(i)));
                let mass = state.density * branch.mesh.volumes()[i];
                local.add(mass, mass * potential)?;
            }
            let local = local.totals();
            sum.add(local.mass, local.energy)?;
        }
        for (i, plenum) in plenums.iter().enumerate() {
            poll(cancelled)?;
            let state = candidate_plenums.map_or(plenum.state, |states| states[i]);
            let mass = state.density * plenum.volume_m3;
            sum.add(
                mass,
                mass * (self.plenum_potentials_j_kg[i] - self.reference()),
            )?;
        }
        Ok(sum.totals())
    }
}

impl GasNetwork1d {
    /// Current discrete gravitational energy [J], zero without a declaration.
    /// Physical total energy is `totals().absolute_energy(gas())` plus this value.
    pub fn gravitational_potential_energy(&self) -> f64 {
        let storage = self
            .centered_potential_storage(false, &mut || false)
            .expect("admitted potential storage and an infallible callback cannot refuse");
        storage.energy + self.gravity.as_ref().map_or(0.0, NetworkGravity::reference) * storage.mass
    }

    pub(super) fn centered_potential_storage(
        &self,
        candidates: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NetworkTotals, Error> {
        let Some(gravity) = &self.gravity else {
            return Ok(NetworkTotals::default());
        };
        gravity.storage(
            &self.branches,
            &self.plenums,
            candidates.then_some(self.work.plenum_candidate.as_slice()),
            cancelled,
        )
    }
}

#[cfg(test)]
mod tests;
