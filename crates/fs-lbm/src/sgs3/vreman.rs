//! Vreman (2004), Physics of Fluids 16, 3670, doi:10.1063/1.1785131,
//! equations (5)-(8). Here G[i][j] = du_i/dx_j (the paper's alpha transpose).
//! beta = G diag(Delta_x^2, Delta_y^2, Delta_z^2) G^T.
//!
//! Cauchy-Binet writes beta's second invariant as the sum of the squares
//! of nine weighted 2x2 minors of G. This avoids cancellation between large
//! Gram products and needs no negative-invariant clipping. Coordinate filter
//! widths stay distinct; rotating only the flow on a fixed anisotropic mesh
//! is not a symmetry. The coefficient cv is LINEAR, unlike Cs^2 or Cw^2.

use super::SgsError;
use fs_math::det;

/// Allocation-free local Vreman kinematic eddy viscosity, in length^2/time.
///
/// `gradient[i][j] = du_i/dx_j`; `filter_widths[j]` belongs to derivative
/// coordinate j, NOT velocity component i. The full gradient is required.
/// A zero coefficient or zero/rank-one gradient returns zero. Lengths must
/// all be positive and finite, including for the disabled control.
///
/// Power-of-two normalization, a scaled minor norm and exponent-tracked
/// dimensional reconstruction avoid intermediate fourth powers. A nonzero
/// quantity lost during normalization or final dimensional under/overflow
/// refuses instead of silently becoming a laminar fallback. No averaging,
/// test filtering, wall-distance field, heap allocation, or viscosity cap
/// is introduced here; the owning solver retains its explicit ratio gate.
///
/// This is an algebraic LES closure, not a transition/separation prediction
/// or a replacement for unresolved-wall physics. Published Vreman validation
/// does not validate the caller's discretization or resolution.
pub fn vreman_eddy_viscosity(
    cv: f64,
    filter_widths: [f64; 3],
    gradient: [[f64; 3]; 3],
) -> Result<f64, SgsError> {
    if !cv.is_finite() || cv < 0.0 {
        return Err(SgsError::InvalidCoefficient);
    }
    if filter_widths.iter().any(|d| !d.is_finite() || *d <= 0.0) {
        return Err(SgsError::InvalidFilterWidth);
    }
    let mut maximum: f64 = 0.0;
    for row in gradient {
        for value in row {
            if !value.is_finite() {
                return Err(SgsError::NonFiniteGradient);
            }
            maximum = maximum.max(value.abs());
        }
    }
    if maximum == 0.0 || cv == 0.0 {
        return Ok(0.0);
    }
    let g_scale = power_of_two(parts(maximum).1);
    let mut g = [[0.0; 3]; 3];
    let mut norm2 = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            g[i][j] = gradient[i][j] / g_scale;
            if gradient[i][j] != 0.0 && g[i][j] == 0.0 {
                return Err(SgsError::NonFiniteViscosity);
            }
            norm2 += g[i][j] * g[i][j];
        }
    }
    let pairs = [(0, 1), (0, 2), (1, 2)];
    let mut minors = [(0.0, 0); 9];
    let mut largest_exponent = i32::MIN;
    for (row, &(i, j)) in pairs.iter().enumerate() {
        for (column, &(k, l)) in pairs.iter().enumerate() {
            let a = g[i][k];
            let b = g[j][l];
            let c = g[i][l];
            let e = g[j][k];
            let p = a * b;
            let q = c * e;
            if (p == 0.0 && a != 0.0 && b != 0.0)
                || (q == 0.0 && c != 0.0 && e != 0.0)
            {
                return Err(SgsError::NonFiniteViscosity);
            }
            // Recover both rounded-product residuals. Explicit mul_add has
            // fused semantics; no compiler contraction policy is required.
            let minor = ((p - q) + (a.mul_add(b, -p) - c.mul_add(e, -q))).abs();
            if minor != 0.0 {
                let weighted = product_parts([minor, filter_widths[k], filter_widths[l]]);
                largest_exponent = largest_exponent.max(weighted.1);
                minors[3 * row + column] = weighted;
            }
        }
    }
    if largest_exponent == i32::MIN {
        return Ok(0.0);
    }
    let mut squared = 0.0;
    for (mantissa, exponent) in minors {
        if mantissa != 0.0 && exponent - largest_exponent >= -1074 {
            let value = mantissa * power_of_two(exponent - largest_exponent);
            squared += value * value;
        }
    }
    let rate = det::sqrt(squared / norm2);
    positive_product([cv, g_scale, rate], largest_exponent)
}

// Exact decomposition of a finite positive float, including subnormals.
// x = mantissa * 2^exponent, with 1 <= mantissa < 2.
fn parts(x: f64) -> (f64, i32) {
    let bits = x.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    if exponent == 0 {
        let (mantissa, exponent) = parts(x * 4_503_599_627_370_496.0);
        (mantissa, exponent - 52)
    } else {
        (f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | (1023u64 << 52)), exponent - 1023)
    }
}

fn power_of_two(exponent: i32) -> f64 {
    if exponent >= -1022 {
        f64::from_bits(((exponent + 1023) as u64) << 52)
    } else {
        f64::from_bits(1u64 << (exponent + 1074))
    }
}

fn product_parts<const N: usize>(factors: [f64; N]) -> (f64, i32) {
    let mut mantissa = 1.0;
    let mut exponent = 0;
    for value in factors {
        let (m, e) = parts(value);
        mantissa *= m;
        exponent += e;
    }
    let (mantissa, carry) = parts(mantissa);
    (mantissa, exponent + carry)
}

fn positive_product(factors: [f64; 3], extra_exponent: i32) -> Result<f64, SgsError> {
    let (mantissa, exponent) = product_parts(factors);
    let exponent = exponent + extra_exponent;
    if !(-1075..=1023).contains(&exponent) {
        return Err(SgsError::NonFiniteViscosity);
    }
    let value = if exponent < -1022 {
        (mantissa * power_of_two(exponent + 1022)) * f64::MIN_POSITIVE
    } else {
        mantissa * power_of_two(exponent)
    };
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(SgsError::NonFiniteViscosity)
    }
}

#[cfg(test)]
mod tests;
