//! Bound the energy omitted by representing original calorics with high words.
//!
//! Species transport retains independent compensated amounts, while the
//! original caloric API accepts one f64 per amount. A very large sensible
//! reference can turn an otherwise tiny omitted amount into resolvable energy.
//! Refuse that regime instead of interpreting transported reference energy as
//! heat. This is a range check, not a replacement two-word thermodynamic law.

use super::*;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

pub(super) fn admit_inventory_resolution(
    initial: &MulticomponentNonisothermalField,
    next: &MulticomponentNonisothermalField,
    cell: usize,
    old: &MulticomponentCaloricStateV1,
    new: &MulticomponentCaloricStateV1,
    tolerance_j: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    poll(cancelled)?;
    let before = endpoint_bound(initial, cell, old, cancelled)?;
    let after = endpoint_bound(next, cell, new, cancelled)?;
    let bound_j = upper_sum(before, after);
    if bound_j > tolerance_j {
        return Err(
            MulticomponentNonisothermalError::UnresolvedCaloricInventory {
                cell,
                bound_j,
                tolerance_j,
            },
        );
    }
    Ok(())
}

fn endpoint_bound(
    field: &MulticomponentNonisothermalField,
    cell: usize,
    state: &MulticomponentCaloricStateV1,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let n = field.network().species_count();
    let remainders = &field.state.amount_remainders[cell * n..(cell + 1) * n];
    if remainders.iter().all(|&r| r == 0.0) {
        return Ok(0.0);
    }
    let mut remainder_sum = 0.0;
    let mut pure_bound = 0.0;
    for ((&high, &low), source) in state
        .component_moles()
        .iter()
        .zip(remainders)
        .zip(field.caloric_model().heat_capacities())
    {
        poll(cancelled)?;
        // Positive endpoints imply positive amounts all along the straight
        // segment from the sampled inventory to its compensated inventory.
        if !low.is_finite() || !positive(high) || !positive(high + low) {
            return Err(MulticomponentNonisothermalError::Invalid(
                "positive finite compensated caloric component support required",
            ));
        }
        if low == 0.0 {
            continue;
        }
        let magnitude = low.abs();
        remainder_sum = upper_sum(remainder_sum, magnitude);
        let sensible = source
            .enthalpy_change_j_per_mol(source.reference_temperature(), state.temperature())
            .map_err(MulticomponentCaloricErrorV1::from)?;
        pure_bound = upper_sum(pure_bound, upper_product(magnitude, sensible.abs()));
    }
    let maximum_b = field.owner.maximum_inverse_temperature_magnitude_k;
    // For every positive composition, |partial H_mix_i| is at most
    // R*(|sum_j B_ij*x_j|+|sum_j<k B_jk*x_j*x_k|) <= 1.5*R*max|B|.
    // The mean-value theorem then bounds the full nonlinear mixing change;
    // there is no small-remainder linearization or assumed cancellation.
    let mut factors = [remainder_sum, maximum_b, (1.5 * R).next_up()];
    factors.sort_by(f64::total_cmp);
    let mixing_bound = upper_product(upper_product(factors[0], factors[2]), factors[1]);
    Ok(upper_sum(pure_bound, mixing_bound))
}

// Every term is nonnegative. Outward rounding keeps the bound conservative,
// including subnormal products, without depending on cancellation. Ordered
// three-factor products above avoid avoidable intermediate range loss.
fn upper_product(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 {
        0.0
    } else {
        let product = a * b;
        if product.is_finite() {
            product.next_up()
        } else {
            f64::INFINITY
        }
    }
}

fn upper_sum(a: f64, b: f64) -> f64 {
    if a == 0.0 {
        b
    } else if b == 0.0 {
        a
    } else {
        let sum = a + b;
        if sum.is_finite() {
            sum.next_up()
        } else {
            f64::INFINITY
        }
    }
}
