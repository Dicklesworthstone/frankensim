//! A bounded numerical starting point in the original transfer/charge box.
//!
//! This is only a feasibility construction. The original finite-rate equations
//! and both physical inventories still own acceptance in the caller.

use super::super::{CAPACITY, OsmoticNetworkError, Result};

fn unavailable<T>() -> Result<T> {
    Err(OsmoticNetworkError::Unbracketed.into())
}

fn pivot_transfer(values: &[f64], valences: &[i32], pivot: usize) -> Option<f64> {
    let z = f64::from(valences[pivot]);
    let (mut sum, mut correction) = (0.0, 0.0);
    for (i, &value) in values.iter().enumerate() {
        if i != pivot {
            super::super::compensated_add(
                &mut sum,
                &mut correction,
                -(f64::from(valences[i]) / z) * value,
            );
        }
    }
    sum.is_finite().then_some(sum)
}

/// Retain the preferred coordinates when possible, otherwise blend toward a
/// feasible charge extremum. Blocked and absent channels arrive as [0, 0].
/// Recover the pivot from the same charge-nullspace relation as the integrator;
/// no counterion, concentration floor, or equilibrium partition is supplied.
pub(super) fn feasible_seed(
    lower: &[f64],
    upper: &[f64],
    valences: &[i32],
    pivot: Option<usize>,
    transfers: &mut [f64],
) -> Result<()> {
    let n = transfers.len();
    if n > CAPACITY
        || lower.len() != n
        || upper.len() != n
        || valences.len() != n
        || (0..n).any(|i| {
            !lower[i].is_finite()
                || !upper[i].is_finite()
                || !transfers[i].is_finite()
                || lower[i] > upper[i]
        })
    {
        return unavailable();
    }
    for i in 0..n {
        transfers[i] = transfers[i].clamp(lower[i], upper[i]);
    }
    let Some(pivot) = pivot else {
        return if (0..n).all(|i| valences[i] == 0 || transfers[i] == 0.0) {
            Ok(())
        } else {
            unavailable()
        };
    };
    if pivot >= n || valences[pivot] == 0 {
        return unavailable();
    }
    let initial =
        pivot_transfer(transfers, valences, pivot).ok_or(OsmoticNetworkError::Unbracketed)?;
    if (lower[pivot]..=upper[pivot]).contains(&initial) {
        transfers[pivot] = initial;
        return Ok(());
    }

    // Only charged free coordinates contribute. Normalize by the pivot charge
    // before summation, avoiding an unnecessary large-valence charge scale.
    let increase = initial < lower[pivot];
    let mut extreme = [0.0; CAPACITY];
    extreme[..n].copy_from_slice(transfers);
    for i in 0..n {
        if i != pivot && valences[i] != 0 {
            let positive = valences[i].is_positive() != valences[pivot].is_positive();
            extreme[i] = if increase == positive {
                upper[i]
            } else {
                lower[i]
            };
        }
    }
    let end =
        pivot_transfer(&extreme[..n], valences, pivot).ok_or(OsmoticNetworkError::Unbracketed)?;
    let lo = lower[pivot].max(initial.min(end));
    let hi = upper[pivot].min(initial.max(end));
    if lo > hi || end == initial {
        return unavailable();
    }
    // An interior pivot avoids turning a rounded cap into a false refusal.
    // A singleton intersection retains its exact boundary semantics.
    let target = if lo.is_sign_positive() == hi.is_sign_positive() {
        lo + 0.5 * (hi - lo)
    } else {
        0.5 * lo + 0.5 * hi
    };
    let scale = initial.abs().max(end.abs()).max(target.abs());
    let weight =
        ((target / scale - initial / scale) / (end / scale - initial / scale)).clamp(0.0, 1.0);
    if !weight.is_finite() {
        return unavailable();
    }
    for i in 0..n {
        if i != pivot && valences[i] != 0 {
            transfers[i] =
                ((1.0 - weight) * transfers[i] + weight * extreme[i]).clamp(lower[i], upper[i]);
        }
    }
    let value =
        pivot_transfer(transfers, valences, pivot).ok_or(OsmoticNetworkError::Unbracketed)?;
    if !(lower[pivot]..=upper[pivot]).contains(&value) {
        return unavailable();
    }
    transfers[pivot] = value;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cap_driven_seed_retains_charge_and_original_bounds() {
        let lower = [0.2, -0.4, -0.1, 0.0, -0.3];
        let upper = [0.5, 0.2, 0.1, 0.0, 0.3];
        let z = [1, 1, -1, 2, 0];
        let mut seed = [0.0, 0.0, 0.0, 0.0, 0.17];
        feasible_seed(&lower, &upper, &z, Some(2), &mut seed).unwrap();
        assert!((0..seed.len()).all(|i| (lower[i]..=upper[i]).contains(&seed[i])));
        assert_eq!(seed[0] + seed[1] - seed[2], 0.0);
        assert_eq!(seed[3], 0.0);
        assert_eq!(seed[4], 0.17);
    }

    #[test]
    fn mixed_valence_feasibility_is_independent_of_equilibrium() {
        let mut seed = [0.0; 3];
        feasible_seed(
            &[0.1, -0.6, 0.0],
            &[0.4, -0.2, 0.0],
            &[2, 1, -1],
            Some(0),
            &mut seed,
        )
        .unwrap();
        assert_eq!(2.0 * seed[0] + seed[1] - seed[2], 0.0);
        assert_eq!(seed[2], 0.0);
        assert!(
            feasible_seed(
                &[0.2, 0.2, 0.0],
                &[0.3, 0.3, 0.0],
                &[1, 1, -1],
                Some(0),
                &mut seed,
            )
            .is_err()
        );
    }

    #[test]
    fn singleton_trace_boundary_is_not_rounded_to_zero() {
        let trace = f64::from_bits(1);
        let mut seed = [0.0; 2];
        feasible_seed(&[trace, 0.0], &[trace, trace], &[1, -1], Some(0), &mut seed).unwrap();
        assert_eq!(seed, [trace, trace]);
    }
}
