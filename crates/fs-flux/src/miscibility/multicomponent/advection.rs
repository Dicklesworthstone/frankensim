//! Prescribed stirring of every original component in the fixed-volume frame.
//!
//! Both directional face rates are retained: opposing streams mix liquids even
//! when net volume flow is zero. Each synchronous donor step transports every
//! component independently, including its compensated low word. Closed volume
//! continuity and the outgoing-volume CFL bound preserve cell capacities.
//!
//! The original homogeneous free energy is convex in this numerical lane, so
//! donor mixing must pass its actual whole-interval free-energy descent gate.
//! This is first-order numerical mixing, without interface reconstruction or a
//! mechanical-work/heat account. The volume velocity is not a mass-averaged
//! velocity when the original molecular masses differ.

use super::{
    MulticomponentMiscibilityError as Error, MulticomponentMiscibilityField,
    MulticomponentMiscibilityNetwork, MulticomponentMiscibilitySolve, Result, Sum, finite,
    nonnegative, poll, positive,
};
use std::sync::Arc;

/// Nonnegative forward/reverse volume rates in canonical physical face order.
pub use crate::miscibility::advection::BinaryMiscibilityAdvectionFace as MulticomponentMiscibilityAdvectionFace;

/// Immutable prescribed flow bound to the complete original numerical owner.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentMiscibilityAdvection {
    network: MulticomponentMiscibilityNetwork,
    face_rates: Arc<[MulticomponentMiscibilityAdvectionFace]>,
    donor_rates_s: Arc<[[f64; 2]]>,
    maximum_outgoing_rate_s: f64,
    maximum_relative_volume_divergence: f64,
}
impl MulticomponentMiscibilityAdvection {
    pub fn face_rates(&self) -> &[MulticomponentMiscibilityAdvectionFace] {
        &self.face_rates
    }
    /// Largest total outgoing rate divided by its original cell volume [1/s].
    pub const fn maximum_outgoing_rate_s(&self) -> f64 {
        self.maximum_outgoing_rate_s
    }
    /// Largest actual preparation residual; no continuity projection is applied.
    pub const fn maximum_relative_volume_divergence(&self) -> f64 {
        self.maximum_relative_volume_divergence
    }

    pub(super) fn admit(&self, network: &MulticomponentMiscibilityNetwork) -> Result<()> {
        if !Arc::ptr_eq(&self.network.owner, &network.owner) && self.network != *network {
            return Err(Error::Invalid(
                "advection requires its original complete field model and graph",
            ));
        }
        Ok(())
    }

    pub(super) fn substeps(
        &self,
        duration_s: f64,
        controls: MulticomponentMiscibilityAdvectionControls,
    ) -> Result<usize> {
        admit_controls(controls)?;
        if !positive(duration_s) {
            return Err(Error::Invalid("positive finite advection interval required"));
        }
        let required = (duration_s * self.maximum_outgoing_rate_s / controls.maximum_courant)
            .ceil()
            .max(1.0);
        if !required.is_finite() || required > controls.maximum_substeps as f64 {
            return Err(Error::AdvectionWorkLimit {
                required_substeps: required as usize,
                maximum_substeps: controls.maximum_substeps,
            });
        }
        Ok(required as usize)
    }
}

/// Bounded donor work and original-law accuracy for one complete interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentMiscibilityAdvectionControls {
    /// Largest total outgoing volume fraction per substep, in (0,1].
    pub maximum_courant: f64,
    pub maximum_substeps: usize,
    /// Applies once to every species balance and every original cell capacity.
    pub amount_balance_tolerance_mol: f64,
    /// Whole-interval actual convex-law free-energy allowance, in joules.
    pub free_energy_tolerance_j: f64,
}
impl Default for MulticomponentMiscibilityAdvectionControls {
    fn default() -> Self {
        Self {
            maximum_courant: 0.5,
            maximum_substeps: 65_536,
            amount_balance_tolerance_mol: 1e-12,
            free_energy_tolerance_j: 1e-9,
        }
    }
}

pub(super) fn admit_controls(controls: MulticomponentMiscibilityAdvectionControls) -> Result<()> {
    if !positive(controls.maximum_courant)
        || controls.maximum_courant > 1.0
        || !(1..=1_000_000).contains(&controls.maximum_substeps)
        || !nonnegative(controls.amount_balance_tolerance_mol)
        || !nonnegative(controls.free_energy_tolerance_j)
    {
        return Err(Error::Invalid(
            "CFL in (0,1], bounded substeps and finite conservation/free-energy allowances required",
        ));
    }
    Ok(())
}

/// Integrated physical face transfers and an independently audited endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentMiscibilityAdvectionStep {
    pub duration_s: f64,
    pub substeps: usize,
    pub maximum_courant: f64,
    /// Signed canonical-face then canonical-species transfers [mol]. Each
    /// opposing stream is accumulated before their net transfer is reported.
    pub transferred_mol: Vec<Vec<f64>>,
    pub component_residual_mol: Vec<f64>,
    pub maximum_capacity_residual_mol: f64,
    /// Direct change of the original homogeneous Gibbs law, not dissipated heat.
    pub free_energy_change_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentMiscibilityAdvectionCandidate {
    field: MulticomponentMiscibilityField,
    step: MulticomponentMiscibilityAdvectionStep,
}
impl MulticomponentMiscibilityAdvectionCandidate {
    pub const fn field(&self) -> &MulticomponentMiscibilityField {
        &self.field
    }
    pub const fn step(&self) -> &MulticomponentMiscibilityAdvectionStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        MulticomponentMiscibilityField,
        MulticomponentMiscibilityAdvectionStep,
    ) {
        (self.field, self.step)
    }
}

impl MulticomponentMiscibilityNetwork {
    /// Prepare a closed flow on the original canonical faces. The relative
    /// divergence allowance is numerical only, in [0,1e-12], and every accepted
    /// endpoint must independently satisfy its full molar-capacity allowance.
    pub fn prepare_advection(
        &self,
        face_rates: Vec<MulticomponentMiscibilityAdvectionFace>,
        maximum_relative_volume_divergence: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentMiscibilityAdvection> {
        poll(&mut cancelled)?;
        if face_rates.len() != self.faces().len()
            || !maximum_relative_volume_divergence.is_finite()
            || !(0.0..=1e-12).contains(&maximum_relative_volume_divergence)
        {
            return Err(Error::Invalid(
                "one directional rate pair per canonical face and bounded divergence allowance required",
            ));
        }
        let mut incoming = vec![Sum::default(); self.cells().len()];
        let mut outgoing = vec![Sum::default(); self.cells().len()];
        let mut donor_rates_s = Vec::with_capacity(face_rates.len());
        for (face, rate) in self.faces().iter().zip(&face_rates) {
            poll(&mut cancelled)?;
            if !nonnegative(rate.forward_m3_s)
                || !nonnegative(rate.reverse_m3_s)
                || (face.area_m2 == 0.0
                    && (rate.forward_m3_s != 0.0 || rate.reverse_m3_s != 0.0))
            {
                return Err(Error::Invalid(
                    "finite nonnegative rates and no flow through a zero-area face required",
                ));
            }
            let [i, j] = face.cells;
            outgoing[i].add(rate.forward_m3_s)?;
            incoming[j].add(rate.forward_m3_s)?;
            outgoing[j].add(rate.reverse_m3_s)?;
            incoming[i].add(rate.reverse_m3_s)?;
            let mut donor = [0.0; 2];
            for (direction, (cell, rate)) in [(i, rate.forward_m3_s), (j, rate.reverse_m3_s)]
                .into_iter()
                .enumerate()
            {
                donor[direction] = finite(rate / self.cells()[cell].volume_m3, "donor volume rate")?;
                if rate > 0.0 && donor[direction] == 0.0 {
                    return Err(Error::Arithmetic("positive donor volume rate underflows"));
                }
            }
            donor_rates_s.push(donor);
        }
        let mut maximum_outgoing_rate_s = 0.0_f64;
        let mut maximum_divergence = 0.0_f64;
        for (cell, (inflow, outflow)) in incoming.iter().zip(&outgoing).enumerate() {
            poll(&mut cancelled)?;
            let mut residual = Sum::default();
            for value in [inflow.high, inflow.low, -outflow.high, -outflow.low] {
                residual.add(value)?;
            }
            let inflow = inflow.value()?;
            let outflow = outflow.value()?;
            let scale = inflow.max(outflow);
            let relative = if scale == 0.0 {
                0.0
            } else {
                finite(residual.value()?.abs() / scale, "relative volume divergence")?
            };
            if relative > maximum_relative_volume_divergence {
                return Err(Error::VolumeDivergence {
                    cell,
                    inflow_m3_s: inflow,
                    outflow_m3_s: outflow,
                    relative_residual: relative,
                });
            }
            maximum_divergence = maximum_divergence.max(relative);
            maximum_outgoing_rate_s = maximum_outgoing_rate_s.max(finite(
                outflow / self.cells()[cell].volume_m3,
                "total outgoing volume rate",
            )?);
        }
        poll(&mut cancelled)?;
        Ok(MulticomponentMiscibilityAdvection {
            network: self.clone(),
            face_rates: face_rates.into(),
            donor_rates_s: donor_rates_s.into(),
            maximum_outgoing_rate_s,
            maximum_relative_volume_divergence: maximum_divergence,
        })
    }
}

impl MulticomponentMiscibilityField {
    /// Transport every component over one closed physical interval. All faces
    /// read the same preceding donor state before any cell update is applied.
    /// Scratch is reused across substeps within the original numerical payload
    /// allowance. Fractions are refreshed once at the complete endpoint because
    /// donor transport consumes independent inventories directly.
    pub fn advect_candidate(
        &self,
        flow: &MulticomponentMiscibilityAdvection,
        duration_s: f64,
        controls: MulticomponentMiscibilityAdvectionControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentMiscibilityAdvectionCandidate> {
        poll(&mut cancelled)?;
        flow.admit(self.network())?;
        let substeps = flow.substeps(duration_s, controls)?;
        let elapsed_s = finite(self.elapsed_s + duration_s, "advected elapsed time")?;
        if elapsed_s <= self.elapsed_s {
            return Err(Error::Arithmetic("advection interval is unrepresentable"));
        }
        let n = self.network.species_count();
        let mut field = self.clone();
        let mut changes = vec![Sum::default(); self.amounts.len()];
        let mut transfers = vec![Sum::default(); self.network.faces().len() * n];
        let mut advanced_s = 0.0;
        let mut maximum_courant = 0.0_f64;
        for substep in 0..substeps {
            poll(&mut cancelled)?;
            let dt = (duration_s - advanced_s) / (substeps - substep) as f64;
            if !positive(dt) || advanced_s + dt <= advanced_s {
                return Err(Error::Arithmetic("advection substep is unrepresentable"));
            }
            let courant = finite(dt * flow.maximum_outgoing_rate_s, "advection Courant number")?;
            if courant > controls.maximum_courant * (1.0 + 16.0 * f64::EPSILON) {
                return Err(Error::Arithmetic("advection substep exceeds outgoing-volume CFL"));
            }
            maximum_courant = maximum_courant.max(courant);
            changes.fill(Sum::default());
            for (face_index, (face, donor)) in self.network.faces().iter()
                .zip(flow.donor_rates_s.iter()).enumerate()
            {
                poll(&mut cancelled)?;
                let [i, j] = face.cells;
                for (direction, (from, to)) in [(i, j), (j, i)].into_iter().enumerate() {
                    let fraction = finite(donor[direction] * dt, "transported donor fraction")?;
                    if donor[direction] > 0.0 && !positive(fraction) {
                        return Err(Error::Arithmetic("positive donor fraction underflows"));
                    }
                    if fraction == 0.0 { continue; }
                    let sign = if direction == 0 { 1.0 } else { -1.0 };
                    for species in 0..n {
                        let from = from * n + species;
                        let to = to * n + species;
                        let high = finite(field.amounts[from] * fraction, "advected component amount")?;
                        let low = finite(field.amount_remainders[from] * fraction, "advected component remainder")?;
                        if high == 0.0 {
                            return Err(Error::Arithmetic("positive component transfer underflows"));
                        }
                        for moved in [high, low] {
                            changes[from].add(-moved)?;
                            changes[to].add(moved)?;
                            transfers[face_index * n + species].add(sign * moved)?;
                        }
                    }
                }
            }
            for cell in 0..self.network.cells().len() {
                poll(&mut cancelled)?;
                for species in 0..n {
                    let index = cell * n + species;
                    let mut value = (field.amounts[index], field.amount_remainders[index]);
                    for change in [changes[index].high, changes[index].low] {
                        value = super::super::solve::add_amount(value.0, value.1, change)
                            .map_err(|_| Error::Arithmetic("conservative advected amount update"))?;
                    }
                    if !positive(value.0) {
                        return Err(Error::Positivity { cell, species, trial_mol: value.0 });
                    }
                    field.amounts[index] = value.0;
                    field.amount_remainders[index] = value.1;
                }
            }
            advanced_s = if substep + 1 == substeps { duration_s } else { advanced_s + dt };
        }
        for cell in 0..self.network.cells().len() {
            poll(&mut cancelled)?;
            let fractions = self.network.fractions(&field.amounts[cell * n..(cell + 1) * n], &mut cancelled)?;
            field.fractions[cell * n..(cell + 1) * n].copy_from_slice(&fractions[..n]);
        }
        field.elapsed_s = elapsed_s;
        let (component_residual_mol, maximum_capacity_residual_mol, free_energy_change_j) = self.balances(
            &field,
            MulticomponentMiscibilitySolve {
                amount_balance_tolerance_mol: controls.amount_balance_tolerance_mol,
                free_energy_tolerance_j: controls.free_energy_tolerance_j,
                ..MulticomponentMiscibilitySolve::default()
            },
            &mut cancelled,
        )?;
        let mut transferred_mol = Vec::with_capacity(self.network.faces().len());
        for row in transfers.chunks_exact(n) {
            poll(&mut cancelled)?;
            transferred_mol.push(row.iter().map(Sum::value).collect::<Result<Vec<_>>>()?);
        }
        poll(&mut cancelled)?;
        Ok(MulticomponentMiscibilityAdvectionCandidate {
            field,
            step: MulticomponentMiscibilityAdvectionStep {
                duration_s,
                substeps,
                maximum_courant,
                transferred_mol,
                component_residual_mol,
                maximum_capacity_residual_mol,
                free_energy_change_j,
            },
        })
    }
}

#[cfg(test)]
mod tests;
