//! Positivity of the complete source pencil, including prescribed coordinates.
//! A positive free block alone cannot admit negative support inertia or storage.

use super::{SupportedSolidError, invalid, poll, zeros};
use fs_exec::Cx;

/// Numerical positivity admission; the original operators are never changed.
/// M must have strictly positive Cholesky pivots. K may have rigid modes:
/// complete diagonal pivoting retains each resolved positive direction, and
/// the entire unresolved Schur residual must fit the declared relative norm.
/// This is an ordinary floating-point check, not an interval certificate.
pub(super) fn admit_operators(
    stiffness: &[f64],
    mass: &[f64],
    n: usize,
    relative_tolerance: f64,
    cx: &Cx<'_>,
) -> Result<(), SupportedSolidError> {
    poll(cx)?;
    let square = n
        .checked_mul(n)
        .ok_or_else(|| invalid("operator shape overflow"))?;
    if n == 0 || stiffness.len() != square || mass.len() != square {
        return Err(invalid("complete source operator shape required"));
    }
    if !relative_tolerance.is_finite() || !(0.0..1.0).contains(&relative_tolerance) {
        return Err(invalid(
            "finite relative positivity allowance below one required",
        ));
    }
    let mut mass_scale = 0.0_f64;
    let mut stiffness_scale = 0.0_f64;
    for row in 0..n {
        poll(cx)?;
        for column in 0..n {
            let index = row * n + column;
            if !mass[index].is_finite() || !stiffness[index].is_finite() {
                return Err(invalid("finite complete operators required"));
            }
            mass_scale = mass_scale.max(mass[index].abs());
            stiffness_scale = stiffness_scale.max(stiffness[index].abs());
        }
    }
    if mass_scale == 0.0 {
        return Err(invalid("original complete mass must be positive definite"));
    }
    // Strict SPD admission on the original-dimensional mass, normalized only
    // for arithmetic. The original retained mass and physical energy are intact.
    // Each inner product is one vector-sized cancellation tile.
    {
        let mut factor = zeros(square)?;
        for row in 0..n {
            for column in 0..=row {
                poll(cx)?;
                let mut value = mass[row * n + column] / mass_scale;
                for k in 0..column {
                    value = (-factor[row * n + k]).mul_add(factor[column * n + k], value);
                }
                if row == column {
                    if !value.is_finite() || value <= 0.0 {
                        return Err(invalid("original complete mass must be positive definite"));
                    }
                    factor[row * n + column] = value.sqrt();
                } else {
                    let entry = value / factor[column * n + column];
                    if !entry.is_finite() {
                        return Err(invalid("unrepresentable complete mass factor"));
                    }
                    factor[row * n + column] = entry;
                }
            }
        }
    }
    if stiffness_scale == 0.0 {
        return poll(cx);
    }
    let mut residual = zeros(square)?;
    for row in 0..n {
        poll(cx)?;
        for column in 0..n {
            residual[row * n + column] = stiffness[row * n + column] / stiffness_scale;
        }
    }
    let mut column = zeros(n)?;
    for first in 0..n {
        let mut infinity_norm = 0.0_f64;
        let mut pivot = first;
        for row in first..n {
            poll(cx)?;
            let mut sum = 0.0;
            let mut correction = 0.0;
            for j in first..n {
                let add = residual[row * n + j].abs() - correction;
                let next = sum + add;
                correction = (next - sum) - add;
                sum = next;
            }
            if !sum.is_finite() {
                return Err(invalid("unrepresentable stiffness residual norm"));
            }
            infinity_norm = infinity_norm.max(sum);
            if residual[row * n + row] > residual[pivot * n + pivot] {
                pivot = row;
            }
        }
        // Test the whole residual, not its diagonal or individual entries:
        // a zero diagonal with nonzero coupling is indefinite and must refuse.
        if infinity_norm <= relative_tolerance {
            return poll(cx);
        }
        if residual[pivot * n + pivot] <= 0.0 {
            return Err(SupportedSolidError::Operator {
                defect: infinity_norm,
                allowance: relative_tolerance,
            });
        }
        if pivot != first {
            for j in first..n {
                residual.swap(first * n + j, pivot * n + j);
            }
            for i in first..n {
                residual.swap(i * n + first, i * n + pivot);
            }
        }
        let root = residual[first * n + first].sqrt();
        for row in first + 1..n {
            column[row] = residual[row * n + first] / root;
            if !column[row].is_finite() {
                return Err(invalid("unrepresentable complete stiffness factor"));
            }
        }
        for row in first + 1..n {
            poll(cx)?;
            for j in first + 1..=row {
                let value = (-column[row]).mul_add(column[j], residual[row * n + j]);
                if !value.is_finite() {
                    return Err(invalid("unrepresentable stiffness Schur residual"));
                }
                residual[row * n + j] = value;
                residual[j * n + row] = value;
            }
        }
    }
    poll(cx)
}
