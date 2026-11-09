//! Conservative donor-volume transport of original species and parcel enthalpy.
//!
//! Opposing streams retain distinct donor compositions and temperatures. Every
//! substep samples all donors before changing any receiver, moves every amount
//! and its low word independently, and recovers temperature from the original
//! changed-inventory caloric inverse. Nonideal mixing heat is already included
//! in the transported enthalpy and is never added as another heat source.
//!
//! This is first-order upwind mixing on the original closed incompressible
//! volume flow. Its acceptance condition is complete-interval H/S and species
//! balance, not Gibbs energy at the network's reference temperature. Large
//! absolute caloric datums that cannot resolve transported heat within the
//! offered allowance explicitly refuse. The arithmetic resolution estimate is
//! an admission guard, not an interval certificate or source uncertainty bound.

use super::*;
use crate::miscibility::multicomponent::advection::{
    MulticomponentMiscibilityAdvection, MulticomponentMiscibilityAdvectionControls,
};

/// Outgoing-volume CFL, bounded work and whole-interval species allowance.
/// Thermal energy/entropy allowances belong to `MulticomponentNonisothermalSolve`.
pub use crate::miscibility::advection::BinaryMiscibilityAdvectionControls as MulticomponentNonisothermalAdvectionControls;

/// Complete physical receipts in canonical face/species order.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalAdvectionStep {
    pub duration_s: f64,
    pub substeps: usize,
    pub maximum_courant: f64,
    /// Positive from the smaller to the larger cell index, in mol.
    pub transferred_mol: Vec<Vec<f64>>,
    /// Original absolute homogeneous donor enthalpy, in that orientation [J].
    pub transferred_energy_j: Vec<f64>,
    pub caloric_iterations: u64,
    pub component_residual_mol: Vec<f64>,
    pub maximum_capacity_residual_mol: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_per_k: f64,
    /// Accumulated arithmetic resolution allowance for sampled donor H and
    /// its directed transfers. Accepted within one quarter of the energy
    /// budget; not measured dissipation or a source-data uncertainty estimate.
    pub parcel_energy_resolution_bound_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalAdvectionCandidate {
    field: MulticomponentNonisothermalField,
    step: MulticomponentNonisothermalAdvectionStep,
}
impl MulticomponentNonisothermalAdvectionCandidate {
    pub const fn field(&self) -> &MulticomponentNonisothermalField {
        &self.field
    }
    pub const fn step(&self) -> &MulticomponentNonisothermalAdvectionStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        MulticomponentNonisothermalField,
        MulticomponentNonisothermalAdvectionStep,
    ) {
        (self.field, self.step)
    }
}

impl MulticomponentNonisothermalField {
    /// Advect each original component and its donor enthalpy over one interval.
    /// The smaller species tolerance applies once to the complete interval;
    /// local caloric solves and precision checks share the energy allowance.
    /// No isothermal free-energy acceptance or temperature average is used.
    pub fn advect_candidate(
        &self,
        flow: &MulticomponentMiscibilityAdvection,
        duration_s: f64,
        advection: MulticomponentNonisothermalAdvectionControls,
        mut controls: MulticomponentNonisothermalSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentNonisothermalAdvectionCandidate> {
        poll(&mut cancelled)?;
        admit_controls(controls)?;
        flow.admit(self.network())?;
        // Only flow admission and CFL planning are shared with the isothermal
        // owner; none of its thermodynamic evolution or G gate is invoked.
        let substeps = flow.substeps(
            duration_s,
            MulticomponentMiscibilityAdvectionControls {
                maximum_courant: advection.maximum_courant,
                maximum_substeps: advection.maximum_substeps,
                amount_balance_tolerance_mol: advection.amount_balance_tolerance_mol,
                free_energy_tolerance_j: 0.0,
            },
        )?;
        let target = finite(self.elapsed_s() + duration_s, "thermal donor elapsed time")?;
        if target <= self.elapsed_s() {
            return Err(MulticomponentNonisothermalError::Invalid(
                "thermal advection interval is unrepresentable",
            ));
        }
        controls.amount_balance_tolerance_mol = controls
            .amount_balance_tolerance_mol
            .min(advection.amount_balance_tolerance_mol);
        let cells = self.network().cells().len();
        let faces = self.network().faces().len();
        let n = self.network().species_count();
        let local_energy_budget =
            controls.energy_tolerance_j / (4.0 * cells as f64 * substeps as f64);
        let local = MulticomponentLiquidTransportSolveV1 {
            temperature_bracket: self.temperature_bracket(),
            temperature_tolerance_k: controls.temperature_tolerance_k,
            energy_tolerance_j: local_energy_budget,
            maximum_iterations: controls.maximum_caloric_iterations,
        };
        let parcel_budget = 0.25 * controls.energy_tolerance_j;
        // Reuse two private fields and all numerical accumulation buffers.
        // Their payload is covered by the parent's thermal workspace admission.
        let mut field = self.clone();
        let mut next = self.clone();
        let mut states = Vec::with_capacity(cells);
        let mut sample_resolution = vec![0.0; cells];
        let mut changes = vec![Sum::default(); cells * n];
        let mut heat = vec![Sum::default(); cells];
        let mut transfers = vec![Sum::default(); faces * n];
        let mut energy = vec![Sum::default(); faces];
        let mut advanced = 0.0;
        let mut maximum_courant = 0.0_f64;
        let mut caloric_iterations = 0_u64;
        let mut parcel_energy_resolution_bound_j = 0.0;
        for substep in 0..substeps {
            poll(&mut cancelled)?;
            let dt = (duration_s - advanced) / (substeps - substep) as f64;
            if !positive(dt) || advanced + dt <= advanced {
                return Err(MulticomponentNonisothermalError::Invalid(
                    "thermal donor substep is unrepresentable",
                ));
            }
            let courant = finite(dt * flow.maximum_outgoing_rate_s(), "thermal donor CFL")?;
            if courant > advection.maximum_courant * (1.0 + 16.0 * f64::EPSILON) {
                return Err(MulticomponentNonisothermalError::Invalid(
                    "thermal donor substep exceeds outgoing-volume CFL",
                ));
            }
            maximum_courant = maximum_courant.max(courant);
            states.clear();
            for (cell, resolution) in sample_resolution.iter_mut().enumerate() {
                poll(&mut cancelled)?;
                let state = field.caloric_state(cell, &mut cancelled)?;
                *resolution = parcel_sample_resolution(
                    self.caloric_model(),
                    self.owner.maximum_inverse_temperature_magnitude_k,
                    &state,
                    &mut cancelled,
                )?;
                states.push(state);
            }
            changes.fill(Sum::default());
            heat.fill(Sum::default());
            for (face, donor) in flow.donor_rates_s().iter().enumerate() {
                poll(&mut cancelled)?;
                let [i, j] = self.network().faces()[face].cells;
                for (direction, (from, to, sign)) in
                    [(i, j, 1.0), (j, i, -1.0)].into_iter().enumerate()
                {
                    let fraction = finite(donor[direction] * dt, "transported donor fraction")?;
                    if donor[direction] > 0.0 && !positive(fraction) {
                        return Err(MulticomponentNonisothermalError::Invalid(
                            "positive thermal donor fraction underflows",
                        ));
                    }
                    if fraction == 0.0 {
                        continue;
                    }
                    for species in 0..n {
                        poll(&mut cancelled)?;
                        let from_index = from * n + species;
                        let to_index = to * n + species;
                        let high = finite(
                            field.state.amounts[from_index] * fraction,
                            "advected component amount",
                        )?;
                        let low = finite(
                            field.state.amount_remainders[from_index] * fraction,
                            "advected component remainder",
                        )?;
                        if !positive(high)
                            || (field.state.amount_remainders[from_index] != 0.0 && low == 0.0)
                        {
                            return Err(MulticomponentNonisothermalError::Invalid(
                                "nonzero thermal donor component transfer underflows",
                            ));
                        }
                        for moved in [high, low] {
                            changes[from_index].add(-moved)?;
                            changes[to_index].add(moved)?;
                            transfers[face * n + species].add(sign * moved)?;
                        }
                    }
                    let enthalpy = states[from].enthalpy_j();
                    let moved = finite(enthalpy * fraction, "advected original parcel enthalpy")?;
                    if enthalpy != 0.0 && moved == 0.0 {
                        return Err(MulticomponentNonisothermalError::Invalid(
                            "nonzero advected parcel enthalpy underflows",
                        ));
                    }
                    // Include donor sampling, fraction/product rounding and
                    // both local accumulations. Accumulation is compensated;
                    // this conservative magnitude allowance does not assume
                    // cancellation between large opposite enthalpy streams.
                    let stream_resolution = upper_sum(
                        upper_product(2.0 * fraction, sample_resolution[from]),
                        upper_product(32.0 * f64::EPSILON, moved.abs()),
                    );
                    parcel_energy_resolution_bound_j =
                        upper_sum(parcel_energy_resolution_bound_j, stream_resolution);
                    if parcel_energy_resolution_bound_j > parcel_budget {
                        return Err(MulticomponentNonisothermalError::UnresolvedParcelEnthalpy {
                            bound_j: parcel_energy_resolution_bound_j,
                            tolerance_j: parcel_budget,
                        });
                    }
                    heat[from].add(-moved)?;
                    heat[to].add(moved)?;
                    energy[face].add(sign * moved)?;
                }
            }
            // Every flux above used the preceding substep. Receivers never
            // become donors during the same sweep, including on cycles.
            for cell in 0..cells {
                poll(&mut cancelled)?;
                for species in 0..n {
                    let index = cell * n + species;
                    let mut value = (
                        field.state.amounts[index],
                        field.state.amount_remainders[index],
                    );
                    for change in [changes[index].high, changes[index].low] {
                        value = crate::miscibility::solve::add_amount(value.0, value.1, change)
                            .map_err(|_| {
                                MulticomponentMiscibilityError::Arithmetic(
                                    "conservative thermal donor amount update",
                                )
                            })?;
                    }
                    if !positive(value.0) {
                        return Err(MulticomponentMiscibilityError::Positivity {
                            cell,
                            species,
                            trial_mol: value.0,
                        }
                        .into());
                    }
                    next.state.amounts[index] = value.0;
                    next.state.amount_remainders[index] = value.1;
                }
                let range = cell * n..(cell + 1) * n;
                let enthalpy_change = heat[cell].value()?;
                let recovered;
                let new_state = if next.state.amounts[range.clone()]
                    == field.state.amounts[range.clone()]
                    && enthalpy_change == 0.0
                {
                    next.temperatures[cell] = field.temperatures[cell];
                    next.state.fractions[range.clone()]
                        .copy_from_slice(&field.state.fractions[range]);
                    &states[cell]
                } else {
                    recovered = self.caloric_model().homogeneous_transport_candidate(
                        &states[cell],
                        &next.state.amounts[range.clone()],
                        enthalpy_change,
                        local,
                        &mut cancelled,
                    )?;
                    caloric_iterations = caloric_iterations
                        .checked_add(recovered.iterations() as u64)
                        .ok_or(MulticomponentMiscibilityError::Arithmetic(
                            "thermal donor caloric-work count",
                        ))?;
                    next.temperatures[cell] = recovered.final_state().temperature();
                    next.state.fractions[range]
                        .copy_from_slice(recovered.final_state().mixing().mole_fractions());
                    recovered.final_state()
                };
                precision::admit_inventory_resolution(
                    &field,
                    &next,
                    cell,
                    &states[cell],
                    new_state,
                    local_energy_budget,
                    &mut cancelled,
                )?;
            }
            core::mem::swap(&mut field, &mut next);
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
        let mut transferred_mol = Vec::with_capacity(faces);
        let mut transferred_energy_j = Vec::with_capacity(faces);
        for face in 0..faces {
            poll(&mut cancelled)?;
            transferred_mol.push(
                transfers[face * n..(face + 1) * n]
                    .iter()
                    .map(Sum::value)
                    .collect::<core::result::Result<Vec<_>, MulticomponentMiscibilityError>>()?,
            );
            transferred_energy_j.push(energy[face].value()?);
        }
        poll(&mut cancelled)?;
        Ok(MulticomponentNonisothermalAdvectionCandidate {
            field,
            step: MulticomponentNonisothermalAdvectionStep {
                duration_s,
                substeps,
                maximum_courant,
                transferred_mol,
                transferred_energy_j,
                caloric_iterations,
                component_residual_mol,
                maximum_capacity_residual_mol,
                energy_residual_j,
                entropy_change_j_per_k,
                parcel_energy_resolution_bound_j,
            },
        })
    }
}

fn parcel_sample_resolution(
    model: &MulticomponentLiquidCaloricModelV1,
    maximum_inverse_temperature_magnitude_k: f64,
    state: &MulticomponentCaloricStateV1,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let mut magnitude = 0.0;
    for (&amount, source) in state.component_moles().iter().zip(model.heat_capacities()) {
        poll(cancelled)?;
        let h = source
            .enthalpy_change_j_per_mol(source.reference_temperature(), state.temperature())
            .map_err(MulticomponentCaloricErrorV1::from)?;
        magnitude = upper_sum(magnitude, upper_product(amount, h.abs()));
    }
    // Bound the original sum of absolute B-pair terms, not the possibly
    // cancelled net excess H. For normalized positive x, sum_i<j x_i*x_j<=1/2.
    // Balance factors before multiplication to avoid needless overflow in
    // max|B|*N when their dimensionless prefactor remains representable.
    let mut factors = [
        0.5 * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
        maximum_inverse_temperature_magnitude_k,
        state.total_moles(),
    ];
    factors.sort_by(f64::total_cmp);
    magnitude = upper_sum(
        magnitude,
        upper_product(upper_product(factors[0], factors[2]), factors[1]),
    );
    if magnitude == 0.0 {
        return Ok(0.0);
    }
    // Allow for the original per-component products and compensated sum,
    // final f64 publication, and the finite original source samples. This
    // guards arithmetic resolution; it is not a new pure-Cp uncertainty model.
    let allowance = 16.0 * (state.component_moles().len() + 1) as f64 * f64::EPSILON;
    let h = state.enthalpy_j().abs();
    Ok(upper_sum(
        upper_product(allowance, magnitude),
        h.next_up() - h,
    ))
}

fn upper_product(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 {
        0.0
    } else {
        (a * b).next_up()
    }
}
fn upper_sum(a: f64, b: f64) -> f64 {
    if a == 0.0 {
        b
    } else if b == 0.0 {
        a
    } else {
        (a + b).next_up()
    }
}

#[cfg(test)]
pub(super) mod tests;
