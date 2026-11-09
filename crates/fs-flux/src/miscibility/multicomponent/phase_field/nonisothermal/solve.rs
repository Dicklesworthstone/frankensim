//! Entropy-conjugate full-Maxwell--Stefan rates and half-edge energy recovery.

use super::*;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

pub(super) fn step(
    initial: &MulticomponentNonisothermalPhaseField,
    duration_s: f64,
    controls: MulticomponentNonisothermalSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MulticomponentNonisothermalPhaseFieldCandidate> {
    poll(&mut cancelled)?;
    admit_controls(controls)?;
    let elapsed = finite(
        initial.elapsed_s() + duration_s,
        "thermal phase-field elapsed time",
    )?;
    if !positive(duration_s) || elapsed <= initial.elapsed_s() {
        return Err(MulticomponentNonisothermalError::Invalid(
            "positive representable thermal phase-field interval required",
        ));
    }
    let rates = rates(initial, controls, &mut cancelled)?;
    candidate_from_rates(initial, duration_s, &rates, controls, &mut cancelled)
}

/// Whole-graph entropy-conjugate gradients precede every original friction
/// solve. The immutable rates are reusable between full and first-half trials.
pub(super) fn rates(
    initial: &MulticomponentNonisothermalPhaseField,
    controls: MulticomponentNonisothermalSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<thermal_solve::Rates> {
    poll(cancelled)?;
    admit_controls(controls)?;
    let network = initial.network();
    let n = network.species_count();
    let fractions = initial.mole_fractions();
    let mut gradients = vec![Sum::default(); fractions.len()];
    for (face_index, face) in network.faces().iter().enumerate() {
        poll(cancelled)?;
        let coefficient = network.edge_coefficients_j[face_index];
        if coefficient == 0.0 {
            continue;
        }
        let [i, j] = face.cells;
        let left_temperature = initial.temperatures()[i].value();
        let right_temperature = initial.temperatures()[j].value();
        let low = left_temperature.min(right_temperature);
        let high = left_temperature.max(right_temperature);
        // beta_i+beta_j = (1+T_min/T_max)/T_min. This form avoids
        // constructing a reciprocal-temperature sum or T_i*T_j.
        let inverse_weight = 1.0 + low / high;
        for species in 0..n {
            let left = i * n + species;
            let right = j * n + species;
            let jump = fractions[right] - fractions[left];
            let derivative =
                spatial_solve::pair_force(-jump, coefficient, inverse_weight, low)?;
            gradients[left].add(derivative)?;
            gradients[right].add(-derivative)?;
        }
    }
    let amounts = initial.component_moles();
    let mut reduced = Vec::with_capacity(gradients.len());
    for cell in 0..network.cells().len() {
        poll(cancelled)?;
        let total = amount_sum(&amounts[cell * n..(cell + 1) * n])?;
        for species in 0..n {
            let derivative = gradients[cell * n + species].value()?;
            // R>1 in the declared SI units. Opposing total/R scales combine
            // before division to retain representable reduced potentials.
            let scaled = if total >= 1.0 {
                (derivative / total) / R
            } else {
                derivative / (total * R)
            };
            let value = finite(
                scaled,
                "reduced interface entropy-conjugate species potential",
            )?;
            if derivative != 0.0 && value == 0.0 {
                return Err(MulticomponentMiscibilityError::Arithmetic(
                    "active thermal interface potential underflows",
                )
                .into());
            }
            reduced.push(value);
        }
    }
    poll(cancelled)?;
    thermal_solve::rates_with_gradient(&initial.bulk, Some(&reduced), controls, cancelled)
}

pub(super) fn candidate_from_rates(
    initial: &MulticomponentNonisothermalPhaseField,
    duration_s: f64,
    rates: &thermal_solve::Rates,
    controls: MulticomponentNonisothermalSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<MulticomponentNonisothermalPhaseFieldCandidate> {
    poll(cancelled)?;
    let mut trial = thermal_solve::conservative_trial_from_rates(
        &initial.bulk,
        duration_s,
        rates,
        controls,
        cancelled,
    )?;
    // The complete simultaneous inventory update must determine composition
    // before any interface allocation or local temperature inverse is used.
    thermal_solve::prepare_trial_composition(&mut trial, cancelled)?;
    for (face_index, face) in initial.network.faces().iter().enumerate() {
        poll(cancelled)?;
        let change = interface_solve::edge_energy_change_from_fractions_j(
            initial.network(),
            initial.mole_fractions(),
            trial.field.mole_fractions(),
            face_index,
        )?;
        let half = finite(0.5 * change, "thermal half-edge energy increment")?;
        if change != 0.0 && half == 0.0 {
            return Err(MulticomponentMiscibilityError::Arithmetic(
                "active half-edge energy allocation underflows",
            )
            .into());
        }
        // Paired transport updates total local U. Removing each endpoint's
        // new interface allocation supplies the original bulk H inverse,
        // including for spectators whose own species amounts did not change.
        let [i, j] = face.cells;
        trial.energy_changes[i].add(-half)?;
        trial.energy_changes[j].add(-half)?;
    }
    let ordinary = thermal_solve::recover_candidate(
        &initial.bulk,
        trial,
        rates,
        controls,
        cancelled,
    )?;
    let (bulk, ordinary) = ordinary.into_parts();
    let field = MulticomponentNonisothermalPhaseField {
        network: initial.network.clone(),
        bulk,
        workspace_bytes: initial.workspace_bytes,
    };
    // Audit the actual recovered endpoint, not just the pre-inverse target.
    // The original recovery already checked inventories and caloric precision.
    let interfacial_energy_change_j =
        interfacial_energy_change_j(initial, &field, cancelled)?;
    let energy_residual_j = admit_total_energy_entropy(
        ordinary.energy_residual_j,
        interfacial_energy_change_j,
        ordinary.entropy_change_j_per_k,
        controls,
    )?;
    poll(cancelled)?;
    Ok(MulticomponentNonisothermalPhaseFieldCandidate {
        field,
        step: MulticomponentNonisothermalPhaseFieldStep {
            duration_s,
            transferred_mol: ordinary.transferred_mol,
            transferred_energy_j: ordinary.transferred_energy_j,
            maxwell_stefan_evaluations: ordinary.maxwell_stefan_evaluations,
            maximum_gradient_residual_per_m: ordinary.maximum_gradient_residual_per_m,
            maximum_molar_frame_residual_mol_m2_s:
                ordinary.maximum_molar_frame_residual_mol_m2_s,
            caloric_iterations: ordinary.caloric_iterations,
            component_residual_mol: ordinary.component_residual_mol,
            maximum_capacity_residual_mol: ordinary.maximum_capacity_residual_mol,
            energy_residual_j,
            entropy_change_j_per_k: ordinary.entropy_change_j_per_k,
            homogeneous_enthalpy_change_j: ordinary.energy_residual_j,
            interfacial_energy_change_j,
            entropy_production_rate_w_per_k: ordinary.entropy_production_rate_w_per_k,
        },
    })
}

/// Reuse the exact isothermal phase-field storage law without cloning fields.
pub(super) fn edge_energy_j(
    field: &MulticomponentNonisothermalPhaseField,
    face_index: usize,
) -> Result<f64> {
    Ok(interface_solve::edge_energy_from_fractions_j(
        field.network(),
        field.mole_fractions(),
        face_index,
    )?)
}

/// The same direct quadratic increment serves allocation, audits and the
/// non-cancelling per-edge adaptive error account.
pub(super) fn edge_energy_change_j(
    initial: &MulticomponentNonisothermalPhaseField,
    next: &MulticomponentNonisothermalPhaseField,
    face_index: usize,
) -> Result<f64> {
    Ok(interface_solve::edge_energy_change_from_fractions_j(
        initial.network(),
        initial.mole_fractions(),
        next.mole_fractions(),
        face_index,
    )?)
}

pub(super) fn interfacial_energy_change_j(
    initial: &MulticomponentNonisothermalPhaseField,
    next: &MulticomponentNonisothermalPhaseField,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let mut change = Sum::default();
    let mut storage = Sum::default();
    for face in 0..initial.network.faces().len() {
        poll(cancelled)?;
        change.add(edge_energy_change_j(initial, next, face)?)?;
        storage.add(edge_energy_j(next, face)?)?;
    }
    storage.value()?;
    poll(cancelled)?;
    Ok(change.value()?)
}

pub(super) fn admit_total_energy_entropy(
    homogeneous: f64,
    interfacial: f64,
    entropy: f64,
    controls: MulticomponentNonisothermalSolve,
) -> Result<f64> {
    let mut total = Sum::default();
    total.add(homogeneous)?;
    total.add(interfacial)?;
    let energy = total.value()?;
    thermal::admit_energy_entropy(energy, entropy, controls)?;
    Ok(energy)
}
