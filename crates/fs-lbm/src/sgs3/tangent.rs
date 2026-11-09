//! Analytic full-gradient derivatives of the original local SGS closures.
//!
//! The primal value always comes from the existing closure. The derivative
//! differentiates its mathematical invariants, not a finite-difference probe
//! or a frozen scalar-strain approximation. WALE follows Nicoud & Ducros
//! (1999), doi:10.1023/A:1009995426001, equations (10) and (13); Vreman follows
//! doi:10.1063/1.1785131, equations (5)-(8), in the parent's transposed-gradient
//! convention. Filter widths and coefficients are held fixed.

use super::vreman::{parts, positive_product, power_of_two, product_parts};
use super::{SgsError, SgsModel, eddy_viscosity, vreman_eddy_viscosity};
use fs_math::det;

/// Kinematic eddy viscosity and its response to all nine velocity derivatives.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SgsViscosityTangent {
    /// Exactly the original closure's value, in length squared per time.
    pub viscosity: f64,
    /// `derivative[i][j] = d(nu_t)/d(G[i][j])`, in length squared, where
    /// `G[i][j] = d(u_i)/d(x_j)`. This includes antisymmetric perturbations.
    ///
    /// At zero strain for Smagorinsky, zero/rank-one gradients for Vreman,
    /// and the all-zero WALE gradient, this is the ZERO generalized-derivative
    /// selection. The viscosity is generally not classically differentiable
    /// there. In particular, rank-one Vreman flow can have nonzero one-sided
    /// derivatives: zero does not mean that neighboring flow is laminar.
    /// Away from these sets this is the ordinary analytic derivative. WALE's
    /// derivative at a nonzero gradient with zero traceless square is zero.
    pub derivative: [[f64; 3]; 3],
}

/// Evaluate the original isotropic closure and its full-gradient tangent.
///
/// Input admission and the primal viscosity are unchanged. Tangents add an
/// explicit refusal if a nonzero derivative cannot be represented as a finite
/// nonzero `f64`; no coefficient, gradient, or derivative is clipped. The
/// singular-state selection is documented in [`SgsViscosityTangent`].
///
/// The model's implicit stress needs both terms of
/// `d(2 nu_t S) = 2 nu_t dS + 2 S (derivative:dG)`. The derivative alone is
/// not a stress Jacobian, nor does it supply a nonlinear convergence guarantee.
pub fn eddy_viscosity_tangent(
    model: SgsModel,
    filter_width: f64,
    gradient: [[f64; 3]; 3],
) -> Result<SgsViscosityTangent, SgsError> {
    if let SgsModel::Vreman { cv } = model {
        return vreman_eddy_viscosity_tangent(cv, [filter_width; 3], gradient);
    }
    let viscosity = eddy_viscosity(model, filter_width, gradient)?;
    let mut response = SgsViscosityTangent {
        viscosity,
        derivative: [[0.0; 3]; 3],
    };
    let scale = maximum(gradient);
    let coefficient = model.coefficient();
    if scale == 0.0 || coefficient == 0.0 {
        return Ok(response);
    }
    // Homogeneity cancels the positive normalization scale from d(nu)/dG.
    // It is held constant during differentiation, so max's branch is irrelevant.
    let g = gradient.map(|row| row.map(|value| value / scale));
    let strain = core::array::from_fn(|i| core::array::from_fn(|j| 0.5 * (g[i][j] + g[j][i])));
    if matches!(model, SgsModel::Smagorinsky { .. }) {
        let strain_scale = maximum(strain);
        if strain_scale == 0.0 {
            return Ok(response);
        }
        let normalized = strain.map(|row| row.map(|value| value / strain_scale));
        let rate = det::sqrt(2.0 * squared_norm(normalized));
        for i in 0..3 {
            for j in 0..3 {
                response.derivative[i][j] = signed_product(
                    2.0 * normalized[i][j] / rate,
                    [coefficient, coefficient, filter_width, filter_width],
                )?;
            }
        }
        return Ok(response);
    }
    let square: [[f64; 3]; 3] =
        core::array::from_fn(|i| core::array::from_fn(|j| (0..3).map(|k| g[i][k] * g[k][j]).sum()));
    let trace_third = (square[0][0] + square[1][1] + square[2][2]) / 3.0;
    let deviator = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            0.5 * (square[i][j] + square[j][i]) - if i == j { trace_third } else { 0.0 }
        })
    });
    let d_scale = maximum(deviator);
    if d_scale == 0.0 {
        return Ok(response);
    }
    let d_normalized = deviator.map(|row| row.map(|value| value / d_scale));
    let d_norm = det::sqrt(squared_norm(d_normalized));
    let d_unit = d_normalized.map(|row| row.map(|value| value / d_norm));
    let t = d_scale * d_norm; // t = sqrt(Sd:Sd), not its square.
    let a = squared_norm(strain);
    let a_three_halves = a * det::sqrt(a);
    let b_five_fourths = t * t * det::sqrt(t);
    let denominator = a * a_three_halves + b_five_fourths;
    if !(denominator.is_finite() && denominator > 0.0) {
        return Err(SgsError::NonFiniteViscosity);
    }
    let b_fraction = b_five_fourths / denominator;
    let unit_trace = d_unit[0][0] + d_unit[1][1] + d_unit[2][2];
    for i in 0..3 {
        for j in 0..3 {
            // d(t)/dG_ij = (Sd/t):d(dev(sym(G^2)))/dG_ij.
            let dt = (0..3)
                .map(|k| d_unit[i][k] * g[j][k] + g[k][i] * d_unit[k][j])
                .sum::<f64>()
                - (2.0 / 3.0) * unit_trace * g[j][i];
            // Factor t^2 out before dimensional reconstruction. A very small
            // near-wall derivative can remain representable with large widths.
            let direction = ((3.0 - 2.5 * b_fraction) * dt
                - 5.0 * (t / denominator) * a_three_halves * strain[i][j])
                / denominator;
            response.derivative[i][j] = signed_product(
                direction,
                [coefficient, coefficient, filter_width, filter_width, t, t],
            )?;
        }
    }
    Ok(response)
}

/// Evaluate anisotropic Vreman viscosity and its full-gradient tangent.
///
/// Each width belongs to the derivative coordinate (the gradient's COLUMN),
/// just as in [`vreman_eddy_viscosity`]. Weighted minors retain separate widths
/// and the original compensated products. Exponent tracking avoids squaring
/// dimensional widths or reconstructing an overflowing viscosity/gradient ratio.
/// See [`SgsViscosityTangent`] for the nondifferentiable rank-one selection.
pub fn vreman_eddy_viscosity_tangent(
    cv: f64,
    filter_widths: [f64; 3],
    gradient: [[f64; 3]; 3],
) -> Result<SgsViscosityTangent, SgsError> {
    let viscosity = vreman_eddy_viscosity(cv, filter_widths, gradient)?;
    let mut response = SgsViscosityTangent {
        viscosity,
        derivative: [[0.0; 3]; 3],
    };
    let scale = maximum(gradient);
    if scale == 0.0 || cv == 0.0 {
        return Ok(response);
    }
    let g_scale = power_of_two(parts(scale).1);
    let g = gradient.map(|row| row.map(|value| value / g_scale));
    let norm2 = squared_norm(g);
    let norm = det::sqrt(norm2);
    let pairs = [(0, 1), (0, 2), (1, 2)];
    let mut minors = [(0.0, 0); 9];
    let mut largest = i32::MIN;
    for (row, &(i, j)) in pairs.iter().enumerate() {
        for (column, &(k, l)) in pairs.iter().enumerate() {
            let p = g[i][k] * g[j][l];
            let q = g[i][l] * g[j][k];
            let minor = (p - q) + (g[i][k].mul_add(g[j][l], -p) - g[i][l].mul_add(g[j][k], -q));
            if minor != 0.0 {
                let (mantissa, exponent) =
                    product_parts([minor.abs(), filter_widths[k], filter_widths[l]]);
                minors[3 * row + column] = (mantissa.copysign(minor), exponent);
                largest = largest.max(exponent);
            }
        }
    }
    if largest == i32::MIN {
        return Ok(response);
    }
    let mut squared = 0.0;
    for (mantissa, exponent) in minors {
        if mantissa != 0.0 && exponent - largest >= -1074 {
            let value = mantissa * power_of_two(exponent - largest);
            squared += value * value;
        }
    }
    let minor_norm = det::sqrt(squared);
    for a in 0..3 {
        for b in 0..3 {
            // Four minors contain each entry. A fifth term differentiates
            // 1/||G||. Keep signed products scaled until AFTER their sum.
            let mut terms = [(0.0, 0); 5];
            let mut count = 0;
            for (row, &(i, j)) in pairs.iter().enumerate() {
                for (column, &(k, l)) in pairs.iter().enumerate() {
                    let dm = if a == i && b == k {
                        g[j][l]
                    } else if a == j && b == l {
                        g[i][k]
                    } else if a == i && b == l {
                        -g[j][k]
                    } else if a == j && b == k {
                        -g[i][l]
                    } else {
                        0.0
                    };
                    let (mantissa, exponent) = minors[3 * row + column];
                    if mantissa != 0.0 && dm != 0.0 {
                        let (m, e) = product_parts([
                            cv,
                            mantissa.abs(),
                            filter_widths[k],
                            filter_widths[l],
                            dm.abs(),
                            1.0 / (minor_norm * norm),
                        ]);
                        terms[count] = (m.copysign(mantissa * dm), e + exponent - largest);
                        count += 1;
                    }
                }
            }
            if g[a][b] != 0.0 {
                let (m, e) = product_parts([cv, minor_norm, g[a][b].abs(), 1.0 / (norm * norm2)]);
                terms[count] = (m.copysign(-g[a][b]), e + largest);
            }
            response.derivative[a][b] = signed_sum(terms)?;
        }
    }
    Ok(response)
}

fn maximum(tensor: [[f64; 3]; 3]) -> f64 {
    tensor.iter().flatten().fold(0.0_f64, |a, b| a.max(b.abs()))
}

fn squared_norm(tensor: [[f64; 3]; 3]) -> f64 {
    tensor.iter().flatten().map(|value| value * value).sum()
}

fn signed_product<const N: usize>(value: f64, factors: [f64; N]) -> Result<f64, SgsError> {
    if !value.is_finite() {
        return Err(SgsError::NonFiniteViscosity);
    }
    if value == 0.0 {
        return Ok(0.0);
    }
    let (m, e) = product_parts(factors);
    positive_product([value.abs(), m, 1.0], e).map(|result| result.copysign(value))
}

fn signed_sum(mut terms: [(f64, i32); 5]) -> Result<f64, SgsError> {
    loop {
        let Some(largest) = terms
            .iter()
            .filter(|(m, _)| *m != 0.0)
            .map(|(_, e)| *e)
            .max()
        else {
            return Ok(0.0);
        };
        let mut sum = 0.0_f64;
        let mut correction = 0.0_f64;
        for term in &mut terms {
            if term.0 != 0.0 && term.1 - largest >= -1074 {
                let value = term.0 * power_of_two(term.1 - largest);
                let next = sum + value;
                correction += if sum.abs() >= value.abs() {
                    (sum - next) + value
                } else {
                    (value - next) + sum
                };
                sum = next;
                *term = (0.0, 0);
            }
        }
        let value = sum + correction;
        if value != 0.0 {
            return positive_product([value.abs(), 1.0, 1.0], largest)
                .map(|result| result.copysign(value));
        }
        // Exact cancellation of the large terms must not erase smaller terms
        // whose exponent range exceeded a single normalized f64 window.
    }
}

#[cfg(test)]
mod tests;
