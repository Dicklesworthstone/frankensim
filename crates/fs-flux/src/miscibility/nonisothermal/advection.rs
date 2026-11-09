//! Conservative donor-volume transport of BOTH species and original enthalpy.
//!
//! A stream carries the donor's homogeneous parcel enthalpy, including its
//! nonideal mixing contribution. The same outgoing volume fraction multiplies
//! both component amounts and that enthalpy. Every debit has an opposite credit;
//! all faces sample the preceding substep before any receiving cell is updated.
//! Temperature is recovered from the original caloric model, never averaged.
//! Mixing heat is therefore already in the energy account, not a second source.
//!
//! This is a first-order upwind discretization of incompressible enthalpy
//! advection (see <https://mooseframework.inl.gov/modules/navier_stokes/insfv.html>).
//! Closed divergence-free prescribed flow, equal fixed partial molar volumes,
//! positive common species support, convex homogeneous thermodynamics and zero
//! square-gradient energy retain the parent lane's scope. Numerical mixing
//! produces entropy; spatial/time refinement is necessary. There is no inferred
//! stirring power, viscous heating or momentum feedback.

use super::*;
use crate::miscibility::advection::{
    BinaryMiscibilityAdvection, BinaryMiscibilityAdvectionControls,
};
use crate::miscibility::solve::add_amount;

/// Whole-interval conservative receipts in the original canonical face order.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalAdvectionStep {
    pub duration_s: f64,
    pub substeps: usize,
    pub maximum_courant: f64,
    /// Signed A and B amounts from the smaller to the larger cell index [mol].
    pub transferred_component_moles: Vec<[f64; 2]>,
    /// Signed total homogeneous donor enthalpy in that same orientation [J].
    pub transferred_energy_j: Vec<f64>,
    pub caloric_iterations: u64,
    pub component_residual_mol: [f64; 2],
    pub maximum_capacity_residual_mol: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_per_k: f64,
}

/// No incomplete advection sweep or failed caloric inverse is published.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalAdvectionCandidate {
    field: BinaryNonisothermalField,
    step: BinaryNonisothermalAdvectionStep,
}
impl BinaryNonisothermalAdvectionCandidate {
    pub const fn field(&self) -> &BinaryNonisothermalField {
        &self.field
    }
    pub const fn step(&self) -> &BinaryNonisothermalAdvectionStep {
        &self.step
    }
    pub fn into_parts(self) -> (BinaryNonisothermalField, BinaryNonisothermalAdvectionStep) {
        (self.field, self.step)
    }
}

impl BinaryNonisothermalField {
    /// Stir the original local-temperature field with the prepared volume flow.
    /// Original material diffusivity is qualified at each actual substep state.
    /// The smaller amount allowance and the original energy/entropy allowances
    /// apply ONCE to the complete interval. Temperature inverses share, rather
    /// than multiply, the global energy budget. Work is linear in cells/faces
    /// per CFL substep, with bounded original-caloric scalar inversions.
    pub fn advect_candidate(
        &self,
        flow: &BinaryMiscibilityAdvection,
        duration_s: f64,
        advection: BinaryMiscibilityAdvectionControls,
        mut controls: BinaryNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryNonisothermalAdvectionCandidate> {
        poll(&mut cancelled)?;
        admit_controls(controls)?;
        flow.admit(self.network())?;
        let substeps = flow.substeps(duration_s, advection)?;
        let target = finite(
            self.elapsed_s() + duration_s,
            "thermal advection elapsed time",
        )?;
        if target <= self.elapsed_s() {
            return Err(BinaryNonisothermalError::Invalid(
                "thermal advection time increment is unrepresentable",
            ));
        }
        controls.amount_balance_tolerance_mol = controls
            .amount_balance_tolerance_mol
            .min(advection.amount_balance_tolerance_mol);
        let cells = self.temperatures.len();
        let faces = self.owner.faces.len();
        let local = BinaryLiquidTransportSolveV1 {
            temperature_bracket: self.owner.temperature_bracket,
            temperature_tolerance_k: controls.temperature_tolerance_k,
            energy_tolerance_j: controls.energy_tolerance_j
                / (4.0 * cells as f64 * substeps as f64),
            maximum_iterations: controls.maximum_root_iterations,
        };
        // The parent's conservative workspace admission covers these reused
        // arrays and the retained field; no per-face or per-substep allocation.
        let mut field = self.clone();
        let mut states = Vec::with_capacity(cells);
        let mut changes = vec![[Sum::default(); 2]; cells];
        let mut heat = vec![Sum::default(); cells];
        let mut transfers = vec![[Sum::default(); 2]; faces];
        let mut energy = vec![Sum::default(); faces];
        let mut advanced = 0.0;
        let mut maximum_courant = 0.0_f64;
        let mut caloric_iterations = 0u64;
        for substep in 0..substeps {
            poll(&mut cancelled)?;
            let dt = (duration_s - advanced) / (substeps - substep) as f64;
            if !positive(dt) || advanced + dt <= advanced {
                return Err(BinaryNonisothermalError::Invalid(
                    "thermal advection substep is unrepresentable",
                ));
            }
            let courant = finite(dt * flow.maximum_outgoing_rate_s(), "thermal advection CFL")?;
            if courant > advection.maximum_courant * (1.0 + 16.0 * f64::EPSILON) {
                return Err(BinaryNonisothermalError::Invalid(
                    "thermal advection exceeds outgoing-volume CFL",
                ));
            }
            maximum_courant = maximum_courant.max(courant);
            states.clear();
            for cell in 0..cells {
                poll(&mut cancelled)?;
                states.push(field.caloric_state(cell)?);
            }
            changes.fill([Sum::default(); 2]);
            heat.fill(Sum::default());
            for (face, donor) in flow.donor_rates_s().iter().enumerate() {
                poll(&mut cancelled)?;
                let [i, j] = self.network().owner.faces[face].physical.cells;
                for (direction, (from, to, sign)) in
                    [(i, j, 1.0), (j, i, -1.0)].into_iter().enumerate()
                {
                    let fraction = finite(donor[direction] * dt, "transported donor fraction")?;
                    if donor[direction] > 0.0 && !positive(fraction) {
                        return Err(BinaryNonisothermalError::Invalid(
                            "positive thermal donor fraction underflows",
                        ));
                    }
                    if fraction == 0.0 {
                        continue;
                    }
                    for species in 0..2 {
                        let amount = field.state.amounts[from][species];
                        let moved = finite(amount * fraction, "advected component amount")?;
                        if !positive(moved) {
                            return Err(BinaryNonisothermalError::Invalid(
                                "positive advected component is unrepresentable",
                            ));
                        }
                        changes[from][species].add(-moved)?;
                        changes[to][species].add(moved)?;
                        transfers[face][species].add(sign * moved)?;
                    }
                    let enthalpy = states[from].enthalpy_j();
                    let moved = finite(enthalpy * fraction, "advected original enthalpy")?;
                    if enthalpy != 0.0 && moved == 0.0 {
                        return Err(BinaryNonisothermalError::Invalid(
                            "nonzero advected enthalpy underflows",
                        ));
                    }
                    heat[from].add(-moved)?;
                    heat[to].add(moved)?;
                    energy[face].add(sign * moved)?;
                }
            }
            // All donor samples precede every cell update. In particular a
            // circulation cannot carry newly received material a second face.
            for cell in 0..cells {
                poll(&mut cancelled)?;
                let old_amounts = field.state.amounts[cell];
                for species in 0..2 {
                    let delta = changes[cell][species];
                    let (amount, remainder) = add_amount(
                        field.state.amounts[cell][species],
                        field.state.amount_remainders[cell][species],
                        delta.main,
                    )?;
                    let (amount, remainder) = add_amount(amount, remainder, delta.correction)?;
                    if !positive(amount) {
                        return Err(BinaryNonisothermalError::Invalid(
                            "thermal advection lost positive common species support",
                        ));
                    }
                    field.state.amounts[cell][species] = amount;
                    field.state.amount_remainders[cell][species] = remainder;
                }
                field.state.fractions[cell] = self.network().fractions(
                    field.state.amounts[cell],
                    self.network().owner.capacities_mol[cell],
                )?;
                let enthalpy_change = heat[cell].value()?;
                // Exact rest/uniform states do not accumulate inverse-solve
                // roundoff or consume caloric iterations.
                if field.state.amounts[cell] != old_amounts || enthalpy_change != 0.0 {
                    let recovered = self.caloric_model().homogeneous_transport_candidate(
                        &states[cell],
                        field.state.amounts[cell],
                        enthalpy_change,
                        local,
                        &mut cancelled,
                    )?;
                    caloric_iterations += recovered.iterations() as u64;
                    field.temperatures[cell] = recovered.final_state().temperature();
                }
            }
            field.qualify(&mut cancelled)?;
            advanced = if substep + 1 == substeps {
                duration_s
            } else {
                advanced + dt
            };
        }
        field.state.elapsed_s = target;
        let (
            component_residual_mol,
            maximum_capacity_residual_mol,
            energy_residual_j,
            entropy_change_j_per_k,
        ) = self.balances(&field, controls, &mut cancelled)?;
        let mut transferred_component_moles = Vec::with_capacity(faces);
        let mut transferred_energy_j = Vec::with_capacity(faces);
        for face in 0..faces {
            poll(&mut cancelled)?;
            transferred_component_moles
                .push([transfers[face][0].value()?, transfers[face][1].value()?]);
            transferred_energy_j.push(energy[face].value()?);
        }
        poll(&mut cancelled)?;
        Ok(BinaryNonisothermalAdvectionCandidate {
            field,
            step: BinaryNonisothermalAdvectionStep {
                duration_s,
                substeps,
                maximum_courant,
                transferred_component_moles,
                transferred_energy_j,
                caloric_iterations,
                component_residual_mol,
                maximum_capacity_residual_mol,
                energy_residual_j,
                entropy_change_j_per_k,
            },
        })
    }
}

#[cfg(test)]
pub(super) mod tests;
