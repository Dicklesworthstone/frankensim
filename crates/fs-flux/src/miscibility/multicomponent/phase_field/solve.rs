//! Original full-Maxwell--Stefan transport of total chemical potentials.
//!
//! A canonical edge pass assembles the whole graph's gradient potential before
//! any molecular flux is evaluated. The simultaneous all-species update reuses
//! the ordinary inventory kernel, then applies the interface owner's total
//! free-energy gate instead of the zero-gradient homogeneous gate.

use super::*;
use super::super::solve as bulk_solve;

pub(super) fn step(
    initial: &MulticomponentPhaseField,
    duration_s: f64,
    controls: MulticomponentMiscibilitySolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MulticomponentPhaseFieldCandidate> {
    poll(&mut cancelled)?;
    admit_controls(controls)?;
    let elapsed = finite(initial.elapsed_s() + duration_s, "phase-field elapsed time")?;
    if !positive(duration_s) || elapsed <= initial.elapsed_s() {
        return Err(MulticomponentMiscibilityError::Invalid(
            "positive representable phase-field interval required",
        ));
    }
    let rates = rates(initial, controls, &mut cancelled)?;
    candidate_from_rates(initial, duration_s, &rates, controls, &mut cancelled)
}

/// Initial rates are reusable between full and first-half adaptive trials.
pub(super) fn rates(
    initial: &MulticomponentPhaseField,
    controls: MulticomponentMiscibilitySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<bulk_solve::Rates> {
    poll(cancelled)?;
    admit_controls(controls)?;
    let network = initial.network();
    let n = network.species_count();
    let mut gradients = vec![Sum::default(); initial.bulk.fractions.len()];
    for (face_index, face) in network.faces().iter().enumerate() {
        poll(cancelled)?;
        let coefficient = network.edge_coefficients_j[face_index];
        if coefficient == 0.0 {
            continue;
        }
        let [i, j] = face.cells;
        for species in 0..n {
            let left = i * n + species;
            let right = j * n + species;
            let jump = initial.bulk.fractions[right] - initial.bulk.fractions[left];
            let derivative = bulk_solve::product3(-jump, coefficient, 2.0)?;
            if jump != 0.0 && derivative == 0.0 {
                return Err(MulticomponentMiscibilityError::Arithmetic(
                    "active interfacial chemical force underflows",
                ));
            }
            gradients[left].add(derivative)?;
            gradients[right].add(-derivative)?;
        }
    }
    let rt = finite(R * network.temperature().value(), "phase-field RT")?;
    let mut reduced = Vec::with_capacity(gradients.len());
    for cell in 0..network.cells().len() {
        poll(cancelled)?;
        // Fractions describe the independently admitted, published amounts.
        // Use their actual total, not a nominal capacity that can differ by an
        // admitted roundoff-sized residual.
        let total = amount_sum(&initial.bulk.amounts[cell * n..(cell + 1) * n])?;
        for species in 0..n {
            let derivative = gradients[cell * n + species].value()?;
            // Divisions in the same direction cannot hide a representable
            // final result. Opposing scales must combine first so an
            // intermediate underflow/overflow does not erase a finite effort.
            let scaled = if (total >= 1.0) == (rt >= 1.0) {
                (derivative / total) / rt
            } else {
                derivative / (total * rt)
            };
            let value = finite(
                scaled,
                "reduced multicomponent gradient chemical potential",
            )?;
            if derivative != 0.0 && value == 0.0 {
                return Err(MulticomponentMiscibilityError::Arithmetic(
                    "active gradient chemical potential underflows",
                ));
            }
            reduced.push(value);
        }
    }
    poll(cancelled)?;
    bulk_solve::rates_with_gradient(&initial.bulk, Some(&reduced), controls, cancelled)
}

pub(super) fn candidate_from_rates(
    initial: &MulticomponentPhaseField,
    duration_s: f64,
    rates: &bulk_solve::Rates,
    controls: MulticomponentMiscibilitySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<MulticomponentPhaseFieldCandidate> {
    poll(cancelled)?;
    let candidate = bulk_solve::conservative_candidate_from_rates(
        &initial.bulk,
        duration_s,
        rates,
        controls,
        cancelled,
    )?;
    let (bulk, ordinary) = candidate.into_parts();
    let field = MulticomponentPhaseField {
        network: initial.network.clone(),
        bulk,
    };
    let interfacial_energy_change_j =
        interfacial_energy_change_j(initial, &field, cancelled)?;
    let free_energy_change_j = admit_total_energy(
        ordinary.free_energy_change_j,
        interfacial_energy_change_j,
        controls,
    )?;
    poll(cancelled)?;
    Ok(MulticomponentPhaseFieldCandidate {
        field,
        step: MulticomponentPhaseFieldStep {
            duration_s,
            transferred_mol: ordinary.transferred_mol,
            maxwell_stefan_evaluations: ordinary.maxwell_stefan_evaluations,
            maximum_gradient_residual_per_m: ordinary.maximum_gradient_residual_per_m,
            maximum_molar_frame_residual_mol_m2_s:
                ordinary.maximum_molar_frame_residual_mol_m2_s,
            component_residual_mol: ordinary.component_residual_mol,
            maximum_capacity_residual_mol: ordinary.maximum_capacity_residual_mol,
            bulk_free_energy_change_j: ordinary.free_energy_change_j,
            interfacial_energy_change_j,
            free_energy_change_j,
            dissipation_rate_w: ordinary.dissipation_rate_w,
        },
    })
}

/// Per-edge storage is also the non-cancelling interface error coordinate.
pub(super) fn edge_energy_j(field: &MulticomponentPhaseField, face_index: usize) -> Result<f64> {
    let network = field.network();
    let coefficient = network.edge_coefficients_j[face_index];
    if coefficient == 0.0 {
        return Ok(0.0);
    }
    let [i, j] = network.faces()[face_index].cells;
    let n = network.species_count();
    let mut energy = Sum::default();
    for species in 0..n {
        let jump = field.bulk.fractions[j * n + species]
            - field.bulk.fractions[i * n + species];
        let term = bulk_solve::product3(jump.abs(), jump.abs(), coefficient)?;
        if jump != 0.0 && term == 0.0 {
            return Err(MulticomponentMiscibilityError::Arithmetic(
                "active interface energy underflows",
            ));
        }
        energy.add(term)?;
    }
    energy.value()
}

/// Direct quadratic-energy increment, preserving small endpoint changes
/// without subtracting two large stored edge energies.
pub(super) fn edge_energy_change_j(
    initial: &MulticomponentPhaseField,
    next: &MulticomponentPhaseField,
    face_index: usize,
) -> Result<f64> {
    let network = initial.network();
    let coefficient = network.edge_coefficients_j[face_index];
    if coefficient == 0.0 {
        return Ok(0.0);
    }
    let [i, j] = network.faces()[face_index].cells;
    let n = network.species_count();
    let mut change = Sum::default();
    for species in 0..n {
        let left = i * n + species;
        let right = j * n + species;
        let old_jump = initial.bulk.fractions[right] - initial.bulk.fractions[left];
        let new_jump = next.bulk.fractions[right] - next.bulk.fractions[left];
        let mut difference = Sum::default();
        difference.add(next.bulk.fractions[right] - initial.bulk.fractions[right])?;
        difference.add(-(next.bulk.fractions[left] - initial.bulk.fractions[left]))?;
        let difference = difference.value()?;
        let sum = finite(old_jump + new_jump, "interface midpoint jump")?;
        let product = bulk_solve::product3(difference, sum.abs(), coefficient)?;
        let term = if sum < 0.0 { -product } else { product };
        if difference != 0.0 && sum != 0.0 && term == 0.0 {
            return Err(MulticomponentMiscibilityError::Arithmetic(
                "active interface energy increment underflows",
            ));
        }
        change.add(term)?;
    }
    change.value()
}

pub(super) fn interfacial_energy_change_j(
    initial: &MulticomponentPhaseField,
    next: &MulticomponentPhaseField,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let mut change = Sum::default();
    let mut next_storage = Sum::default();
    for face in 0..initial.network.faces().len() {
        poll(cancelled)?;
        change.add(edge_energy_change_j(initial, next, face)?)?;
        // The stable increment is the acceptance quantity; independently
        // retain representable nonnegative absolute interface storage.
        next_storage.add(edge_energy_j(next, face)?)?;
    }
    next_storage.value()?;
    poll(cancelled)?;
    change.value()
}

pub(super) fn admit_total_energy(
    bulk: f64,
    interfacial: f64,
    controls: MulticomponentMiscibilitySolve,
) -> Result<f64> {
    let mut total = Sum::default();
    total.add(bulk)?;
    total.add(interfacial)?;
    let change_j = total.value()?;
    if change_j > controls.free_energy_tolerance_j {
        return Err(MulticomponentMiscibilityError::FreeEnergy {
            change_j,
            tolerance_j: controls.free_energy_tolerance_j,
        });
    }
    Ok(change_j)
}
