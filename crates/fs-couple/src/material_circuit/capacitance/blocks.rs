//! Explicit fixed-charge replacement of complete reciprocal electric storage.
//!
//! The native represented potential owns the work: delta H is one half of
//! q^T delta(C^-1) q. Subtract inverse coefficients BEFORE charge products;
//! differences of endpoint energies or voltages can erase real cross work.

use super::{add, finite, invalid};
use crate::material_circuit::storage::{MaterialStorageCandidate, StorageParameterWork};
use crate::material_circuit::{
    MAX_MATERIAL_CIRCUIT_BRANCHES, MaterialCircuit, MaterialCircuitError, poll,
};
use fs_phs::circuit::Branch;
use fs_phs::circuit::capacitance::CapacitanceBlock;

impl MaterialCircuit {
    /// Replace complete capacitance matrices at fixed original charge and flux.
    ///
    /// Each existing block retains the same SET of original capacitor branches.
    /// Block-list order and coordinate order inside a block may change. Adding,
    /// removing, splitting or regrouping blocks requires a separate topology /
    /// storage-authority owner and refuses here. Empty-to-empty is a no-op.
    ///
    /// An unchanged physical C, compared in original branch coordinates, keeps
    /// its ORIGINAL prepared inverse and order. Thus a reordered equivalent
    /// source cannot turn refactorization roundoff into physical work. Changed
    /// blocks use their admitted represented inverse; no matrix fit is inferred.
    /// Scalar material, nonlinear curve or dielectric source ownership cannot
    /// overlap the block. Unrelated storage, original source ports, contacts,
    /// resistor bindings and all charge/flux history remain retained.
    ///
    /// Per-branch values are signed q_i v_i/2 contributions to the complete
    /// coupled potential; individual terms need not be positive. Their work is
    /// accumulated from inverse-coefficient differences, not by subtracting
    /// those reported endpoint contributions. Positive parameter work must be
    /// supplied by the changing field/material owner, negative work returned.
    /// It is reversible work, never resistor heat. Native initialization refuses
    /// incompatible ideal clamps rather than changing charge impulsively.
    /// Refusal or cancellation leaves the accepted circuit unchanged.
    pub fn rebind_capacitance_blocks_candidate(
        &self,
        blocks: Vec<CapacitanceBlock>,
        inputs: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialStorageCandidate, MaterialCircuitError> {
        poll(&mut cancelled)?;
        let branch_count = self.graph.branches.len();
        if branch_count > MAX_MATERIAL_CIRCUIT_BRANCHES
            || blocks.len() != self.capacitance_blocks.len()
            || blocks.len() > branch_count
            || blocks
                .iter()
                .try_fold(0usize, |sum, block| sum.checked_add(block.branches().len()))
                .is_none_or(|count| count > branch_count)
        {
            return Err(invalid(
                "capacitance updates require the same bounded block memberships",
            ));
        }
        let mut owned = [false; MAX_MATERIAL_CIRCUIT_BRANCHES];
        for block in &blocks {
            poll(&mut cancelled)?;
            for &branch in block.branches() {
                if branch >= branch_count
                    || owned[branch]
                    || !matches!(self.graph.branches[branch].2, Branch::Capacitor { .. })
                {
                    return Err(invalid(
                        "capacitance updates require distinct original capacitor branches",
                    ));
                }
                if self.storage.iter().any(|source| source.branch == branch)
                    || self
                        .capacitor_curves
                        .iter()
                        .any(|curve| curve.branch == branch)
                    || self
                        .dielectric_capacitors
                        .iter()
                        .any(|source| source.branch() == branch)
                {
                    return Err(invalid(
                        "a coupled-capacitance update cannot replace retained scalar material authority",
                    ));
                }
                owned[branch] = true;
            }
        }
        let blocks = canonical_blocks(&self.capacitance_blocks, &blocks, &mut cancelled)?;
        poll(&mut cancelled)?;
        let mut candidate = Self::new_with_capacitance_blocks(
            self.graph.clone(),
            self.resistors.clone(),
            self.contacts.clone(),
            self.mutual.clone(),
            blocks,
            self.capacitor_curves.to_vec(),
            self.fluxes_wb.clone(),
            self.charges_c.clone(),
            inputs,
        )?;
        candidate.storage = self.storage.clone();
        candidate.magnetic_couplings = self.magnetic_couplings.clone();
        candidate.dielectric_capacitors = self.dielectric_capacitors.clone();
        poll(&mut cancelled)?;

        let mut slots = [usize::MAX; MAX_MATERIAL_CIRCUIT_BRANCHES];
        let mut slot = 0;
        for (branch, &(_, _, element)) in self.graph.branches.iter().enumerate() {
            if matches!(element, Branch::Capacitor { .. }) {
                slots[branch] = slot;
                slot += 1;
            }
        }
        let mut work = Vec::with_capacity(owned[..branch_count].iter().filter(|&&v| v).count());
        for (before, after) in self
            .capacitance_blocks
            .iter()
            .zip(candidate.capacitance_blocks.iter())
        {
            let positions = corresponding_positions(before, after)?;
            let n = before.branches().len();
            for (i, &branch) in before.branches().iter().enumerate() {
                poll(&mut cancelled)?;
                let q = self.charges_c[slots[branch]];
                let (mut old_voltage, mut old_correction) = (0.0, 0.0);
                let (mut new_voltage, mut new_correction) = (0.0, 0.0);
                let (mut change_voltage, mut change_correction) = (0.0, 0.0);
                for (j, &other_branch) in before.branches().iter().enumerate() {
                    let other_q = self.charges_c[slots[other_branch]];
                    let original = before.inverse_matrix_per_f()[i * n + j];
                    let proposed = after.inverse_matrix_per_f()[positions[i] * n + positions[j]];
                    add(
                        &mut old_voltage,
                        &mut old_correction,
                        finite(original * other_q)?,
                    )?;
                    add(
                        &mut new_voltage,
                        &mut new_correction,
                        finite(proposed * other_q)?,
                    )?;
                    let coefficient_change = finite(proposed - original)?;
                    add(
                        &mut change_voltage,
                        &mut change_correction,
                        finite(coefficient_change * other_q)?,
                    )?;
                }
                let contribution = |voltage: f64| -> Result<f64, MaterialCircuitError> {
                    finite(0.5 * finite(q * voltage)?)
                };
                work.push(StorageParameterWork {
                    branch,
                    stored_before_j: contribution(finite(old_voltage + old_correction)?)?,
                    stored_after_j: contribution(finite(new_voltage + new_correction)?)?,
                    supplied_work_j: contribution(finite(change_voltage + change_correction)?)?,
                });
            }
        }
        work.sort_by_key(|entry| entry.branch);
        let (mut total, mut correction) = (0.0, 0.0);
        for entry in &work {
            poll(&mut cancelled)?;
            add(&mut total, &mut correction, entry.supplied_work_j)?;
        }
        poll(&mut cancelled)?;
        Ok(MaterialStorageCandidate {
            circuit: candidate,
            work,
            parameter_work_j: finite(total + correction)?,
        })
    }
}

fn corresponding_positions(
    before: &CapacitanceBlock,
    after: &CapacitanceBlock,
) -> Result<Vec<usize>, MaterialCircuitError> {
    if before.branches().len() != after.branches().len() {
        return Err(invalid(
            "capacitance updates cannot change block membership",
        ));
    }
    before
        .branches()
        .iter()
        .map(|branch| {
            after
                .branches()
                .iter()
                .position(|other| other == branch)
                .ok_or_else(|| invalid("capacitance updates cannot change block membership"))
        })
        .collect()
}

fn canonical_blocks(
    original: &[CapacitanceBlock],
    proposed: &[CapacitanceBlock],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<CapacitanceBlock>, MaterialCircuitError> {
    let mut used = vec![false; proposed.len()];
    let mut canonical = Vec::with_capacity(original.len());
    for before in original {
        poll(cancelled)?;
        let (index, after) = proposed
            .iter()
            .enumerate()
            .find(|(index, after)| {
                !used[*index]
                    && before.branches().len() == after.branches().len()
                    && before
                        .branches()
                        .iter()
                        .all(|branch| after.branches().contains(branch))
            })
            .ok_or_else(|| {
                invalid("capacitance updates cannot add, remove or regroup block branches")
            })?;
        used[index] = true;
        let positions = corresponding_positions(before, after)?;
        let n = before.branches().len();
        let mut unchanged = true;
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..n {
                unchanged &= before.matrix_f()[i * n + j]
                    == after.matrix_f()[positions[i] * n + positions[j]];
            }
        }
        canonical.push(if unchanged {
            before.clone()
        } else {
            after.clone()
        });
    }
    Ok(canonical)
}
