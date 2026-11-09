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
impl Rates {
    pub(super) const fn evaluations(&self) -> usize {
        self.evaluations
    }
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
    rates_with_gradient(initial, None, controls, cancelled)
}

/// The interface owner supplies cell-major extra chemical potentials divided by
/// RT. Bulk forces, the molar gauge, pair friction and receipts remain one kernel.
pub(super) fn rates_with_gradient(
    initial: &MulticomponentMiscibilityField,
    extra_reduced_potentials: Option<&[f64]>,
    controls: MulticomponentMiscibilitySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Rates> {
    poll(cancelled)?;
    admit_controls(controls)?;
    if let Some(potentials) = extra_reduced_potentials {
        if potentials.len() != initial.fractions.len() {
            return Err(MulticomponentMiscibilityError::Invalid(
                "one extra chemical potential per cell and original species required",
            ));
        }
        for row in potentials.chunks(initial.network.species_count()) {
            poll(cancelled)?;
            if row.iter().any(|value| !value.is_finite()) {
                return Err(MulticomponentMiscibilityError::Arithmetic(
                    "finite interfacial chemical potentials required",
                ));
            }
        }
    }
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
            jump.add(ideal_chemical_jump(left, right, s)?)?;
            for t in 0..n {
                jump.add(finite(
                    network.owner.chi[s * n + t] * (right[t] - left[t]),
                    "original excess chemical jump",
                )?)?;
            }
            if let Some(potentials) = extra_reduced_potentials {
                jump.add(finite(
                    potentials[j * n + s] - potentials[i * n + s],
                    "interfacial chemical jump",
                )?)?;
            }
            // The common excess-Gibbs term in every chemical potential is a
            // gauge; remove it analytically before forming the balanced force.
            jumps[s] = jump.value()?;
        }
        for s in 0..n {
            poll(cancelled)?;
            for t in s + 1..n {
                let pair = pair_force(jumps[s] - jumps[t], mean[s], mean[t], face.distance_m)?;
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
    let candidate =
        conservative_candidate_from_rates(initial, duration_s, rates, controls, cancelled)?;
    if candidate.step.free_energy_change_j > controls.free_energy_tolerance_j {
        return Err(MulticomponentMiscibilityError::FreeEnergy {
            change_j: candidate.step.free_energy_change_j,
            tolerance_j: controls.free_energy_tolerance_j,
        });
    }
    poll(cancelled)?;
    Ok(candidate)
}

/// Shared paired inventory update and homogeneous energy increment. The caller
/// must apply its own complete thermodynamic gate before publishing this private
/// candidate: bulk Gibbs for homogeneous liquids, bulk plus interface energy
/// for the phase-field owner. Applying a bulk-only gate here would incorrectly
/// refuse physical smoothing that spends interface energy to increase bulk G.
pub(super) fn conservative_candidate_from_rates(
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
    let (component_residual_mol, maximum_capacity_residual_mol) =
        initial.inventory_balances(&field, controls.amount_balance_tolerance_mol, cancelled)?;
    let free_energy_change_j = free_energy_change(initial, &field, cancelled)?;
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

pub(super) fn product3(a: f64, b: f64, c: f64) -> Result<f64> {
    let mut factors = [a.abs(), b, c];
    factors.sort_by(f64::total_cmp);
    finite(
        ((factors[0] * factors[2]) * factors[1]).copysign(a),
        "paired species transfer",
    )
}

/// Retain the original law's abundant-carrier complement convention. A carrier
/// may round to one while the other independent components remain real traces.
/// Differencing those other components also retains nearby carrier potentials
/// without subtracting two almost equal absolute logarithms.
pub(super) fn ideal_chemical_jump(left: &[f64], right: &[f64], species: usize) -> Result<f64> {
    if left[species] <= 0.5 && right[species] <= 0.5 {
        return finite(
            log_ratio(right[species], left[species]),
            "ideal chemical jump",
        );
    }
    let mut left_other = Sum::default();
    let mut right_other = Sum::default();
    let mut difference = Sum::default();
    for s in 0..left.len() {
        if s != species {
            left_other.add(left[s])?;
            right_other.add(right[s])?;
            difference.add(left[s] - right[s])?;
        }
    }
    let left_other = left_other.value()?;
    let right_other = right_other.value()?;
    let value = if left[species] > 0.5 && right[species] > 0.5 {
        let difference = difference.value()?;
        let relative = difference.abs() / (1.0 - left_other.max(right_other));
        let rounded = 1.0 + relative;
        let logarithm = if rounded == 1.0 {
            relative
        } else {
            ln(rounded) * (relative / (rounded - 1.0))
        };
        logarithm.copysign(difference)
    } else {
        let logarithm = |x: f64, others: f64| {
            if x <= 0.5 {
                ln(x)
            } else if others < 1e-4 {
                // Same bounded trace expansion as the original solution law.
                -others
                    * (1.0 + others * (0.5 + others * (1.0 / 3.0 + others * (0.25 + others * 0.2))))
            } else {
                ln(1.0 - others)
            }
        };
        logarithm(right[species], right_other) - logarithm(left[species], left_other)
    };
    finite(value, "original carrier chemical jump")
}

/// Preserve a representable force even when a trace product underflows before
/// division by distance. Ordinary values keep the original operation order;
/// otherwise exact binary decomposition postpones range rounding until the end.
/// A truly unrepresentable nonzero pair refuses rather than becoming zero force.
pub(super) fn pair_force(contrast: f64, left: f64, right: f64, distance: f64) -> Result<f64> {
    if !contrast.is_finite() || !positive(left) || !positive(right) || !positive(distance) {
        return Err(MulticomponentMiscibilityError::Arithmetic(
            "finite chemical force factors required",
        ));
    }
    if contrast == 0.0 {
        return Ok(contrast);
    }
    let mut direct = contrast;
    let mut ordinary = direct.is_normal();
    for factor in [left, right] {
        direct *= factor;
        ordinary &= direct.is_normal();
    }
    direct /= distance;
    ordinary &= direct.is_normal();
    if ordinary {
        return Ok(direct);
    }
    let mut mantissa = 1.0;
    let mut exponent = 0_i32;
    for value in [contrast.abs(), left, right] {
        let (part, power) = binary_parts(value);
        mantissa *= part;
        exponent += power;
        if mantissa >= 2.0 {
            mantissa *= 0.5;
            exponent += 1;
        }
    }
    let (part, power) = binary_parts(distance);
    mantissa /= part;
    exponent -= power;
    if mantissa < 1.0 {
        mantissa *= 2.0;
        exponent -= 1;
    }
    let magnitude = if exponent > 1023 {
        f64::INFINITY
    } else if exponent >= -1022 {
        mantissa * f64::from_bits(((exponent + 1023) as u64) << 52)
    } else if exponent >= -1075 {
        (mantissa * f64::from_bits(((exponent + 2045) as u64) << 52)) * f64::MIN_POSITIVE
    } else {
        0.0
    };
    if !positive(magnitude) {
        return Err(MulticomponentMiscibilityError::Arithmetic(
            "nonzero chemical force is unrepresentable",
        ));
    }
    Ok(magnitude.copysign(contrast))
}

/// Exact m*2^e with 1<=m<2 for each admitted normal or subnormal factor.
fn binary_parts(value: f64) -> (f64, i32) {
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & 0x000f_ffff_ffff_ffff;
    if exponent != 0 {
        (f64::from_bits((1023_u64 << 52) | fraction), exponent - 1023)
    } else {
        let highest = 63 - fraction.leading_zeros();
        (
            fraction as f64 / (1_u64 << highest) as f64,
            highest as i32 - 1074,
        )
    }
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
