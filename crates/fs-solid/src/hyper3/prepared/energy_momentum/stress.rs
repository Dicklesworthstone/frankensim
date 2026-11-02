//! Objective strain-space discrete gradient using the original material law.
//!
//! Cbar=(C0+C1)/2 is evaluated through F*=chol(Cbar)^T. It is NOT C(Fbar).
//! Symmetric algorithmic second Piola stress supplies potential work within
//! arithmetic/continuous-limit tolerance and zero midpoint internal torque. The analytic
//! tangent differentiates this same map with respect to endpoint F1.

use fs_la::factor::{Cholesky, Lu, cholesky, lu};
use fs_material::MaterialError;
use fs_material::hyper::cards::HyperelasticMaterial;

pub(super) struct AlgorithmicStress {
    pub(super) piola: [f64; 9],
    pub(super) tangent_after: Option<[[f64; 9]; 9]>,
}

fn invalid(what: &str) -> MaterialError {
    MaterialError::State { what: what.into() }
}
fn dot(a: &[f64; 9], b: &[f64; 9]) -> f64 {
    a.iter().zip(b).fold(0.0, |s, (x, y)| x.mul_add(*y, s))
}
fn transpose(a: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| a[3 * (ij % 3) + ij / 3])
}
fn product(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3).fold(0.0, |s, k| {
            a[3 * (ij / 3) + k].mul_add(b[3 * k + ij % 3], s)
        })
    })
}
fn symmetric(a: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| 0.5 * a[ij] + 0.5 * a[3 * (ij % 3) + ij / 3])
}
fn left_solve(factor: &Lu, rhs: &[f64; 9]) -> [f64; 9] {
    let mut result = [0.0; 9];
    for column in 0..3 {
        let mut x = core::array::from_fn::<_, 3, _>(|row| rhs[3 * row + column]);
        factor.solve(&mut x);
        for row in 0..3 {
            result[3 * row + column] = x[row];
        }
    }
    result
}

// Differentiate the native factor C=L L^T; return d(F*)=d(L^T).
fn cholesky_direction(factor: &Cholesky, dc: &[f64; 9]) -> [f64; 9] {
    let mut dl = [0.0; 9];
    for column in 0..3 {
        let mut diagonal = dc[3 * column + column];
        for k in 0..column {
            diagonal -= 2.0 * factor.l(column, k) * dl[3 * column + k];
        }
        dl[3 * column + column] = diagonal / (2.0 * factor.l(column, column));
        for row in column + 1..3 {
            let mut value = dc[3 * row + column];
            for k in 0..column {
                value -=
                    dl[3 * row + k] * factor.l(column, k) + factor.l(row, k) * dl[3 * column + k];
            }
            value -= factor.l(row, column) * dl[3 * column + column];
            dl[3 * row + column] = value / factor.l(column, column);
        }
    }
    transpose(&dl)
}

pub(super) fn evaluate<L: HyperelasticMaterial>(
    material: &L,
    before: &[f64; 9],
    after: &[f64; 9],
    tangent: bool,
) -> Result<AlgorithmicStress, MaterialError> {
    // Both endpoints retain the original source and principal-stretch gates.
    material.first_piola(before)?;
    let p_after = material.first_piola(after)?;
    let energy_change = material.energy_density(after) - material.energy_density(before);
    let midpoint = core::array::from_fn(|i| 0.5 * before[i] + 0.5 * after[i]);
    let delta_f = core::array::from_fn(|i| after[i] - before[i]);
    let c0 = product(&transpose(before), before);
    let c1 = product(&transpose(after), after);
    let cbar = core::array::from_fn(|i| 0.5 * c0[i] + 0.5 * c1[i]);
    // Bilinear identity avoids subtracting two nearly equal absolute metrics.
    let half_delta = product(&transpose(&midpoint), &delta_f);
    let delta_c = core::array::from_fn(|ij| half_delta[ij] + half_delta[3 * (ij % 3) + ij / 3]);
    if !energy_change.is_finite() || cbar.iter().chain(&delta_c).any(|v| !v.is_finite()) {
        return Err(invalid(
            "unrepresentable discrete-gradient metric or energy",
        ));
    }
    let factor =
        cholesky(&cbar, 3).map_err(|_| invalid("positive mean material metric required"))?;
    let representative = core::array::from_fn(|ij| factor.l(ij % 3, ij / 3));
    let inverse =
        lu(&representative, 3).map_err(|_| invalid("singular material metric representative"))?;
    let p = material.first_piola(&representative)?;
    let raw_s = left_solve(&inverse, &p);
    let s = symmetric(&raw_s);
    let denominator = dot(&delta_c, &delta_c);
    let metric_norm = dot(&cbar, &cbar);
    if !denominator.is_finite() || !metric_norm.is_finite() || s.iter().any(|v| !v.is_finite()) {
        return Err(invalid("unrepresentable discrete-gradient stress"));
    }
    // Below this relative strain threshold the energy quotient loses digits.
    // Use its continuous limit; the caller still gates measured energy error.
    let resolved = denominator > 64.0 * f64::EPSILON * metric_norm;
    let correction = if resolved {
        (2.0 * energy_change - dot(&s, &delta_c)) / denominator
    } else {
        0.0
    };
    let algorithmic_s = core::array::from_fn(|i| correction.mul_add(delta_c[i], s[i]));
    let piola = product(&midpoint, &algorithmic_s);
    if piola.iter().any(|v| !v.is_finite()) {
        return Err(invalid("unrepresentable algorithmic Piola stress"));
    }
    let tangent_after = if tangent {
        let material_tangent = material.consistent_tangent(&representative)?;
        let mut result = [[0.0; 9]; 9];
        for column in 0..9 {
            let mut direction = [0.0; 9];
            direction[column] = 1.0;
            let term = product(&transpose(after), &direction);
            let d_delta_c = core::array::from_fn(|ij| term[ij] + term[3 * (ij % 3) + ij / 3]);
            let d_cbar = d_delta_c.map(|x| 0.5 * x);
            let d_representative = cholesky_direction(&factor, &d_cbar);
            let d_p: [f64; 9] =
                core::array::from_fn(|row| dot(&material_tangent[row], &d_representative));
            let frame_term = product(&d_representative, &raw_s);
            let d_raw_s = left_solve(&inverse, &core::array::from_fn(|i| d_p[i] - frame_term[i]));
            let d_s = symmetric(&d_raw_s);
            let d_correction = if resolved {
                let numerator = 2.0 * p_after[column] - dot(&d_s, &delta_c) - dot(&s, &d_delta_c);
                (numerator - correction * 2.0 * dot(&delta_c, &d_delta_c)) / denominator
            } else {
                0.0
            };
            let d_algorithmic_s = core::array::from_fn(|i| {
                d_correction.mul_add(delta_c[i], correction.mul_add(d_delta_c[i], d_s[i]))
            });
            let first = product(&direction, &algorithmic_s);
            let second = product(&midpoint, &d_algorithmic_s);
            for row in 0..9 {
                result[row][column] = 0.5 * first[row] + second[row];
            }
        }
        if result.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable algorithmic material tangent"));
        }
        Some(result)
    } else {
        None
    };
    Ok(AlgorithmicStress {
        piola,
        tangent_after,
    })
}

#[cfg(test)]
mod tests;
