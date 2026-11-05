//! Conservative donor-cell transport in the original fixed-volume frame.
//!
//! Each face retains BOTH directional volume rates. Opposing streams can mix
//! species even when their net volume flux is zero. For component s the signed
//! canonical flux is `c * (Q_forward*x_i,s - Q_reverse*x_j,s)` [mol/s]. A closed
//! divergence-free volume flow preserves each cell's original cV capacity.
//! Unequal molecular masses do not make this velocity barycentric: these rates
//! cannot, by themselves, close a variable-density momentum equation.
//!
//! Synchronous upwind Euler substeps satisfy the outgoing-volume CFL limit.
//! This is first order, with numerical mixing that must be checked by refinement;
//! it is not an interface reconstruction or a higher-order transport scheme.
//! The original molecular/gradient diffusion remains a separate operator.
//! Advection can increase bulk/interfacial free energy through imposed stirring
//! and numerical mixing. Its actual change is reported without imposing the
//! pure-diffusion dissipation inequality or claiming a mechanical work balance.

use super::solve::add_amount;
use super::{
    BinaryMiscibilityError as Error, BinaryMiscibilityNetwork, BinaryMiscibilityState, Result, Sum,
    finite, nonnegative, poll, positive,
};
use std::sync::Arc;

/// Nonnegative volume rates in the network's canonical face order [m³/s].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BinaryMiscibilityAdvectionFace {
    /// From the face's smaller cell index to its larger cell index.
    pub forward_m3_s: f64,
    /// From the face's larger cell index to its smaller cell index.
    pub reverse_m3_s: f64,
}

/// Immutable prescribed flow, prepared once for the original field/geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityAdvection {
    network: BinaryMiscibilityNetwork,
    face_rates: Vec<BinaryMiscibilityAdvectionFace>,
    // Outgoing donor fractions per second, avoiding divisions during substeps.
    donor_rates_s: Vec<[f64; 2]>,
    maximum_outgoing_rate_s: f64,
    maximum_relative_volume_divergence: f64,
}
impl BinaryMiscibilityAdvection {
    pub fn face_rates(&self) -> &[BinaryMiscibilityAdvectionFace] {
        &self.face_rates
    }
    /// Largest outgoing volume rate divided by the original cell volume [1/s].
    pub const fn maximum_outgoing_rate_s(&self) -> f64 {
        self.maximum_outgoing_rate_s
    }
    /// Actual largest relative divergence at preparation, never projected away.
    pub const fn maximum_relative_volume_divergence(&self) -> f64 {
        self.maximum_relative_volume_divergence
    }
}

/// Bounded work and conservation accuracy for a complete transport interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryMiscibilityAdvectionControls {
    /// Maximum outgoing volume fraction per synchronous substep, in (0,1].
    pub maximum_courant: f64,
    /// At most one million substeps can be admitted, independent of time scale.
    pub maximum_substeps: usize,
    /// Whole-interval per-species and per-cell capacity allowance [mol].
    /// This allowance is not multiplied by faces or substeps.
    pub amount_balance_tolerance_mol: f64,
}
impl Default for BinaryMiscibilityAdvectionControls {
    fn default() -> Self {
        Self {
            maximum_courant: 0.5,
            maximum_substeps: 65_536,
            amount_balance_tolerance_mol: 1e-12,
        }
    }
}

/// Actual transport and chemical/interfacial energy change over one interval.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityAdvectionStep {
    pub duration_s: f64,
    pub substeps: usize,
    pub maximum_courant: f64,
    /// Signed net component transfers, in canonical face order [mol]. Opposing
    /// streams are integrated separately before their net transfer is reported.
    pub transferred_component_moles: Vec<[f64; 2]>,
    pub component_residual_mol: [f64; 2],
    pub maximum_capacity_residual_mol: f64,
    pub free_energy_before_j: f64,
    pub free_energy_after_j: f64,
    /// Direct increments retain small changes beside large bulk energies.
    /// Either sign is admissible; this is not dissipated heat or stirring work.
    pub free_energy_change_j: f64,
}

/// Immutable all-or-nothing transport result. The original state is unchanged.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityAdvectionCandidate {
    state: BinaryMiscibilityState,
    step: BinaryMiscibilityAdvectionStep,
}
impl BinaryMiscibilityAdvectionCandidate {
    pub const fn state(&self) -> &BinaryMiscibilityState {
        &self.state
    }
    pub const fn step(&self) -> &BinaryMiscibilityAdvectionStep {
        &self.step
    }
    pub fn into_parts(self) -> (BinaryMiscibilityState, BinaryMiscibilityAdvectionStep) {
        (self.state, self.step)
    }
}

impl BinaryMiscibilityNetwork {
    /// Bind both directional rates to this original closed finite-volume graph.
    /// Relative divergence is `abs(in-out)/max(in,out)` per cell (zero at rest).
    /// Its allowance is numerical only, in [0,1e-12]. No flow projection or
    /// inventory normalization repairs a mismatch; later capacity admission
    /// still applies to every accepted state, including roundoff accumulation.
    pub fn prepare_advection(
        &self,
        face_rates: Vec<BinaryMiscibilityAdvectionFace>,
        maximum_relative_volume_divergence: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityAdvection> {
        poll(&mut cancelled)?;
        if face_rates.len() != self.owner.faces.len()
            || !maximum_relative_volume_divergence.is_finite()
            || !(0.0..=1e-12).contains(&maximum_relative_volume_divergence)
        {
            return Err(Error::Invalid(
                "one directional rate pair per canonical face and bounded divergence allowance required",
            ));
        }
        let mut inflow = vec![Sum::default(); self.owner.cells.len()];
        let mut outflow = vec![Sum::default(); self.owner.cells.len()];
        let mut donor_rates_s = Vec::with_capacity(face_rates.len());
        for (face, rate) in self.owner.faces.iter().zip(&face_rates) {
            poll(&mut cancelled)?;
            if !nonnegative(rate.forward_m3_s)
                || !nonnegative(rate.reverse_m3_s)
                || (face.physical.area_m2 == 0.0
                    && (rate.forward_m3_s != 0.0 || rate.reverse_m3_s != 0.0))
            {
                return Err(Error::Invalid(
                    "finite nonnegative directional rates and no transport through a zero-area face required",
                ));
            }
            let [i, j] = face.physical.cells;
            outflow[i].add(rate.forward_m3_s)?;
            inflow[j].add(rate.forward_m3_s)?;
            outflow[j].add(rate.reverse_m3_s)?;
            inflow[i].add(rate.reverse_m3_s)?;
            let mut donor = [0.0; 2];
            for (slot, (cell, q)) in [(i, rate.forward_m3_s), (j, rate.reverse_m3_s)]
                .into_iter()
                .enumerate()
            {
                donor[slot] = finite(q / self.owner.cells[cell].volume_m3, "donor volume rate")?;
                if q > 0.0 && donor[slot] == 0.0 {
                    return Err(Error::Arithmetic("positive donor volume rate underflows"));
                }
            }
            donor_rates_s.push(donor);
        }
        let mut maximum_outgoing_rate_s: f64 = 0.0;
        let mut maximum_divergence: f64 = 0.0;
        for (cell, (&incoming, &outgoing)) in inflow.iter().zip(&outflow).enumerate() {
            poll(&mut cancelled)?;
            let incoming = incoming.value()?;
            let outgoing = outgoing.value()?;
            let scale = incoming.max(outgoing);
            let relative = if scale == 0.0 {
                0.0
            } else {
                finite(
                    (incoming - outgoing).abs() / scale,
                    "relative volume divergence",
                )?
            };
            if relative > maximum_relative_volume_divergence {
                return Err(Error::VolumeDivergence {
                    cell,
                    inflow_m3_s: incoming,
                    outflow_m3_s: outgoing,
                    relative_residual: relative,
                });
            }
            maximum_divergence = maximum_divergence.max(relative);
            maximum_outgoing_rate_s = maximum_outgoing_rate_s.max(finite(
                outgoing / self.owner.cells[cell].volume_m3,
                "cell outgoing volume rate",
            )?);
        }
        poll(&mut cancelled)?;
        Ok(BinaryMiscibilityAdvection {
            network: self.clone(),
            face_rates,
            donor_rates_s,
            maximum_outgoing_rate_s,
            maximum_relative_volume_divergence: maximum_divergence,
        })
    }

    /// Transport both original species for one physical interval. All outgoing
    /// face transfers use the SAME preceding state, so face order cannot expose
    /// freshly received material to another donor within the same substep.
    /// Work is O(cells+faces) per substep; scratch arrays are reused throughout.
    /// The graph's original 512-byte/cell and /face workspace admission covers
    /// the prepared rates, candidate, transfer totals and transaction scratch.
    pub fn advect_candidate(
        &self,
        initial: &BinaryMiscibilityState,
        flow: &BinaryMiscibilityAdvection,
        duration_s: f64,
        controls: BinaryMiscibilityAdvectionControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityAdvectionCandidate> {
        poll(&mut cancelled)?;
        self.admit(initial)?;
        if !Arc::ptr_eq(&self.owner, &flow.network.owner) && self.owner != flow.network.owner {
            return Err(Error::Invalid(
                "advection requires its original field model and graph",
            ));
        }
        if !positive(duration_s)
            || !positive(controls.maximum_courant)
            || controls.maximum_courant > 1.0
            || !(1..=1_000_000).contains(&controls.maximum_substeps)
            || !nonnegative(controls.amount_balance_tolerance_mol)
        {
            return Err(Error::Invalid(
                "positive time, CFL in (0,1], bounded substeps and finite conservation allowance required",
            ));
        }
        let elapsed_s = finite(initial.elapsed_s + duration_s, "advection elapsed time")?;
        if elapsed_s <= initial.elapsed_s {
            return Err(Error::Arithmetic(
                "advection time increment is unrepresentable",
            ));
        }
        let required = (duration_s * flow.maximum_outgoing_rate_s / controls.maximum_courant)
            .ceil()
            .max(1.0);
        if !required.is_finite() || required > controls.maximum_substeps as f64 {
            return Err(Error::AdvectionWorkLimit {
                required_substeps: required as usize,
                maximum_substeps: controls.maximum_substeps,
            });
        }
        let substeps = required as usize;
        let before = self.free_energy_j(initial, &mut cancelled)?;
        let mut state = initial.clone();
        let mut changes = vec![[Sum::default(); 2]; self.owner.cells.len()];
        let mut transfers = vec![[Sum::default(); 2]; self.owner.faces.len()];
        let mut advanced_s = 0.0;
        let mut maximum_courant: f64 = 0.0;
        for substep in 0..substeps {
            poll(&mut cancelled)?;
            let dt = (duration_s - advanced_s) / (substeps - substep) as f64;
            if !positive(dt) || advanced_s + dt <= advanced_s {
                return Err(Error::Arithmetic("advection substep is unrepresentable"));
            }
            let courant = finite(
                dt * flow.maximum_outgoing_rate_s,
                "advection Courant number",
            )?;
            if courant > controls.maximum_courant * (1.0 + 16.0 * f64::EPSILON) {
                return Err(Error::Arithmetic(
                    "advection substep exceeds outgoing-volume CFL",
                ));
            }
            maximum_courant = maximum_courant.max(courant);
            changes.fill([Sum::default(); 2]);
            for (index, (face, donor)) in
                self.owner.faces.iter().zip(&flow.donor_rates_s).enumerate()
            {
                poll(&mut cancelled)?;
                let [i, j] = face.physical.cells;
                let fraction = [donor[0] * dt, donor[1] * dt];
                for direction in 0..2 {
                    if donor[direction] > 0.0 && !positive(fraction[direction]) {
                        return Err(Error::Arithmetic(
                            "positive transported volume fraction underflows",
                        ));
                    }
                }
                for species in 0..2 {
                    let amounts = [state.amounts[i][species], state.amounts[j][species]];
                    let moved = [amounts[0] * fraction[0], amounts[1] * fraction[1]];
                    for direction in 0..2 {
                        if !nonnegative(moved[direction])
                            || (amounts[direction] > 0.0
                                && fraction[direction] > 0.0
                                && moved[direction] == 0.0)
                        {
                            return Err(Error::Arithmetic(
                                "component advection is unrepresentable",
                            ));
                        }
                    }
                    // Retain cancellation residuals of opposing streams in the
                    // same two-word reduction used for cell increments/receipts.
                    changes[i][species].add(-moved[0])?;
                    changes[i][species].add(moved[1])?;
                    changes[j][species].add(moved[0])?;
                    changes[j][species].add(-moved[1])?;
                    transfers[index][species].add(moved[0])?;
                    transfers[index][species].add(-moved[1])?;
                }
            }
            for (cell, change) in changes.iter().enumerate() {
                poll(&mut cancelled)?;
                for (species, delta) in change.iter().enumerate() {
                    let (amount, remainder) = add_amount(
                        state.amounts[cell][species],
                        state.amount_remainders[cell][species],
                        delta.main,
                    )?;
                    let (amount, remainder) = add_amount(amount, remainder, delta.correction)?;
                    if !nonnegative(amount) {
                        return Err(Error::Arithmetic("advection exhausted a donor component"));
                    }
                    state.amounts[cell][species] = amount;
                    state.amount_remainders[cell][species] = remainder;
                }
                state.fractions[cell] =
                    self.fractions(state.amounts[cell], self.owner.capacities_mol[cell])?;
            }
            advanced_s = if substep + 1 == substeps {
                duration_s
            } else {
                advanced_s + dt
            };
        }
        let mut balances = [Sum::default(); 2];
        let mut capacity_residual: f64 = 0.0;
        for (cell, (old, new)) in initial.amounts.iter().zip(&state.amounts).enumerate() {
            poll(&mut cancelled)?;
            for species in 0..2 {
                balances[species].add(finite(
                    new[species] - old[species],
                    "advected species balance",
                )?)?;
            }
            capacity_residual =
                capacity_residual.max((new[0] + new[1] - self.owner.capacities_mol[cell]).abs());
        }
        let balances = [balances[0].value()?, balances[1].value()?];
        if balances
            .iter()
            .any(|x| x.abs() > controls.amount_balance_tolerance_mol)
            || capacity_residual > controls.amount_balance_tolerance_mol
        {
            return Err(Error::Conservation {
                component_residual_mol: balances,
                capacity_residual_mol: capacity_residual,
            });
        }
        let after = self.free_energy_j(&state, &mut cancelled)?;
        let change = self.free_energy_change(initial, &state, &mut cancelled)?;
        let transferred_component_moles = transfers
            .into_iter()
            .map(|pair| {
                poll(&mut cancelled)?;
                Ok([pair[0].value()?, pair[1].value()?])
            })
            .collect::<Result<Vec<_>>>()?;
        state.elapsed_s = elapsed_s;
        poll(&mut cancelled)?;
        Ok(BinaryMiscibilityAdvectionCandidate {
            state,
            step: BinaryMiscibilityAdvectionStep {
                duration_s,
                substeps,
                maximum_courant,
                transferred_component_moles,
                component_residual_mol: balances,
                maximum_capacity_residual_mol: capacity_residual,
                free_energy_before_j: before,
                free_energy_after_j: after,
                free_energy_change_j: change,
            },
        })
    }
}

#[cfg(test)]
mod tests;
