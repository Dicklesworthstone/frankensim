//! Synchronous original-force evaluations and conservative explicit updates.

use super::*;

/// Face fluxes use one original field. Keeping this stage separate allows the
/// adaptive owner to reuse identical initial rates for full and half trials.
pub(super) struct Rates {
    fluxes: Vec<f64>,
    evaluations: usize,
    maximum_gradient_residual_per_m: f64,
    maximum_molar_frame_residual_mol_m2_s: f64,
    dissipation_rate_w: f64,
}

pub(super) fn step(
    initial: &MulticomponentMiscibilityField,
    duration_s: f64,
    controls: MulticomponentMiscibilitySolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MulticomponentMiscibilityCandidate> {
    poll(&mut cancelled)?;
    admit_controls(controls)?;
    admit_duration(initial, duration_s)?;
    let rates = rates(initial, controls, &mut cancelled)?;
    candidate_from_rates(initial, duration_s, &rates, controls, &mut cancelled)
}

fn admit_duration(initial: &MulticomponentMiscibilityField, duration_s: f64) -> Result<f64> {
    let elapsed = finite(
        initial.elapsed_s + duration_s,
        "multicomponent elapsed time",
    )?;
    if !positive(duration_s) || elapsed <= initial.elapsed_s {
        return Err(MulticomponentMiscibilityError::Invalid(
            "positive representable physical interval required",
        ));
    }
    Ok(elapsed)
}

pub(super) fn rates(
    initial: &MulticomponentMiscibilityField,
    controls: MulticomponentMiscibilitySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Rates> {
    poll(cancelled)?;
    admit_controls(controls)?;
    let network = initial.network();
    let n = network.species_count();
    let mut fluxes = vec![0.0; network.faces().len() * n];
    let mut workspace = network.transport().workspace()?;
    let mut evaluations = 0usize;
    let mut maximum_gradient_residual_per_m = 0.0_f64;
    let mut maximum_molar_frame_residual_mol_m2_s = 0.0_f64;
    let mut dissipation = Sum::default();
    for (face_index, face) in network.faces().iter().enumerate() {
        poll(cancelled)?;
        if face.area_m2 == 0.0 {
            continue;
        }
        let [i, j] = face.cells;
        let left = &initial.fractions[i * n..(i + 1) * n];
        let right = &initial.fractions[j * n..(j + 1) * n];
        let mut mean = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        let mut jumps = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        let mut driving = [Sum::default(); MAX_REGULAR_SOLUTION_SPECIES_V1];
        for s in 0..n {
            poll(cancelled)?;
            // Add before halving so coincident subnormal traces survive.
            mean[s] = 0.5 * (left[s] + right[s]);
            if !positive(mean[s]) {
                return Err(MulticomponentMiscibilityError::Arithmetic(
                    "positive face composition",
                ));
            }
            let mut jump = Sum::default();
            jump.add(log_ratio(right[s], left[s]))?;
            for t in 0..n {
                jump.add(finite(
                    network.owner.chi[s * n + t] * (right[t] - left[t]),
                    "original excess chemical jump",
                )?)?;
            }
            // The common excess-Gibbs term in every chemical potential is a
            // gauge; remove it analytically before forming the balanced force.
            jumps[s] = jump.value()?;
        }
        for s in 0..n {
            poll(cancelled)?;
            for t in s + 1..n {
                let pair = finite(
                    (((jumps[s] - jumps[t]) * mean[s]) * mean[t]) / face.distance_m,
                    "balanced chemical-potential force pair",
                )?;
                driving[s].add(pair)?;
                driving[t].add(-pair)?;
            }
        }
        let mut forces = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        for s in 0..n {
            forces[s] = driving[s].value()?;
        }
        let result = network.transport().evaluate_molar_driving_forces(
            &mean[..n],
            &forces[..n],
            network.concentration(),
            &network.owner.diffusivity_matrices[face.diffusivity_matrix],
            controls.maxwell_stefan,
            &mut workspace,
            &mut *cancelled,
        )?;
        evaluations += 1;
        maximum_gradient_residual_per_m =
            maximum_gradient_residual_per_m.max(result.maximum_gradient_residual_per_m());
        maximum_molar_frame_residual_mol_m2_s =
            maximum_molar_frame_residual_mol_m2_s.max(result.molar_frame_residual_mol_m2_s().abs());
        let mut common = Sum::default();
        for s in 0..n {
            common.add(mean[s] * jumps[s])?;
        }
        let common = common.value()?;
        let mut power = Sum::default();
        for s in 0..n {
            let flux = result.molar_fluxes()[s];
            fluxes[face_index * n + s] = flux;
            power.add(finite(
                -flux * (jumps[s] - common),
                "chemical face dissipation",
            )?)?;
        }
        let power = finite(
            (power.value()? * face.area_m2) * (R * network.temperature().value()),
            "chemical dissipation rate",
        )?;
        if power < 0.0 {
            return Err(MulticomponentMiscibilityError::Arithmetic(
                "negative original chemical face dissipation",
            ));
        }
        dissipation.add(power)?;
    }
    poll(cancelled)?;
    Ok(Rates {
        fluxes,
        evaluations,
        maximum_gradient_residual_per_m,
        maximum_molar_frame_residual_mol_m2_s,
        dissipation_rate_w: dissipation.value()?,
    })
}

pub(super) fn candidate_from_rates(
    initial: &MulticomponentMiscibilityField,
    duration_s: f64,
    rates: &Rates,
    controls: MulticomponentMiscibilitySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<MulticomponentMiscibilityCandidate> {
    poll(cancelled)?;
    admit_controls(controls)?;
    let elapsed = admit_duration(initial, duration_s)?;
    let network = initial.network();
    let n = network.species_count();
    let mut field = initial.clone();
    let mut transferred_mol = Vec::with_capacity(network.faces().len());
    for (face_index, face) in network.faces().iter().enumerate() {
        poll(cancelled)?;
        let [i, j] = face.cells;
        let mut transfers = Vec::with_capacity(n);
        for s in 0..n {
            let flux = rates.fluxes[face_index * n + s];
            let transfer = product3(flux, face.area_m2, duration_s)?;
            if flux != 0.0 && face.area_m2 != 0.0 && transfer == 0.0 {
                return Err(MulticomponentMiscibilityError::Arithmetic(
                    "positive-duration species transfer underflows",
                ));
            }
            for (cell, change) in [(i, -transfer), (j, transfer)] {
                let index = cell * n + s;
                (field.amounts[index], field.amount_remainders[index]) =
                    super::super::solve::add_amount(
                        field.amounts[index],
                        field.amount_remainders[index],
                        change,
                    )
                    .map_err(|_| {
                        MulticomponentMiscibilityError::Arithmetic(
                            "conservative multicomponent amount update",
                        )
                    })?;
            }
            transfers.push(transfer);
        }
        transferred_mol.push(transfers);
    }
    // All source rates were simultaneous. Temporary partial accumulations are
    // never evaluated as physical cells; incoming and outgoing faces finish
    // before positivity is tested, with every low word retained.
    for cell in 0..network.cells().len() {
        poll(cancelled)?;
        for species in 0..n {
            let trial_mol = field.amounts[cell * n + species];
            if !positive(trial_mol) {
                return Err(MulticomponentMiscibilityError::Positivity {
                    cell,
                    species,
                    trial_mol,
                });
            }
        }
        let fractions = network.fractions(&field.amounts[cell * n..(cell + 1) * n], cancelled)?;
        field.fractions[cell * n..(cell + 1) * n].copy_from_slice(&fractions[..n]);
    }
    field.elapsed_s = elapsed;
    let (component_residual_mol, maximum_capacity_residual_mol, free_energy_change_j) =
        initial.balances(&field, controls, cancelled)?;
    poll(cancelled)?;
    Ok(MulticomponentMiscibilityCandidate {
        field,
        step: MulticomponentMiscibilityStep {
            duration_s,
            transferred_mol,
            maxwell_stefan_evaluations: rates.evaluations,
            maximum_gradient_residual_per_m: rates.maximum_gradient_residual_per_m,
            maximum_molar_frame_residual_mol_m2_s: rates.maximum_molar_frame_residual_mol_m2_s,
            component_residual_mol,
            maximum_capacity_residual_mol,
            free_energy_change_j,
            dissipation_rate_w: rates.dissipation_rate_w,
        },
    })
}

fn product3(a: f64, b: f64, c: f64) -> Result<f64> {
    let mut factors = [a.abs(), b, c];
    factors.sort_by(f64::total_cmp);
    finite(
        ((factors[0] * factors[2]) * factors[1]).copysign(a),
        "paired species transfer",
    )
}

/// a ln(a/b)-a+b, retaining the quadratic response near coincident amounts.
fn relative_entropy(a: f64, b: f64) -> Result<f64> {
    if a == b {
        return Ok(0.0);
    }
    let difference = a - b;
    let r = difference / b;
    let value = if r.abs() < 0.125 {
        let mut polynomial = 0.0;
        for k in (2..=22).rev() {
            let coefficient = if k % 2 == 0 { 1.0 } else { -1.0 } / (k * (k - 1)) as f64;
            polynomial = coefficient + r * polynomial;
        }
        (difference * r) * polynomial
    } else {
        a * log_ratio(a, b) - difference
    };
    finite(value, "original ideal relative entropy")
}

/// Increment the original degree-one Gibbs law algebraically. The ideal part
/// uses its relative-entropy remainder, and each excess pair uses a symmetric
/// product difference. No two large absolute Gibbs energies are subtracted.
pub(super) fn free_energy_change(
    initial: &MulticomponentMiscibilityField,
    next: &MulticomponentMiscibilityField,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let n = initial.network.species_count();
    let mut reduced = Sum::default();
    for cell in 0..initial.network.cells().len() {
        poll(cancelled)?;
        let old = &initial.amounts[cell * n..(cell + 1) * n];
        let new = &next.amounts[cell * n..(cell + 1) * n];
        let old_total = amount_sum(old)?;
        let new_total = amount_sum(new)?;
        for s in 0..n {
            poll(cancelled)?;
            reduced.add(finite(
                (new[s] - old[s]) * log_ratio(old[s], old_total),
                "ideal chemical increment",
            )?)?;
            reduced.add(relative_entropy(new[s], old[s])?)?;
            for t in s + 1..n {
                let old_x = old[t] / old_total;
                let new_x = new[t] / new_total;
                let change = (new[s] - old[s]) * (0.5 * old_x + 0.5 * new_x)
                    + (new_x - old_x) * (0.5 * old[s] + 0.5 * new[s]);
                reduced.add(finite(
                    initial.network.owner.chi[s * n + t] * change,
                    "original excess pair increment",
                )?)?;
            }
        }
        reduced.add(-relative_entropy(new_total, old_total)?)?;
    }
    finite(
        (reduced.value()? * R) * initial.network.temperature().value(),
        "original total free-energy change",
    )
}
