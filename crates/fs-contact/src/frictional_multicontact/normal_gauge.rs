//! Global motion determinacy with a normal-load gauge.
//!
//! Let p*, v* be a checked root and c_i = mu_i*s_i*, with
//! s_i* in the subgradient of |v_ti*| and p_ti* = -p_ni* c_i.
//! For ANY other root, norm-subgradient monotonicity and p_ni >= 0 give
//!   delta p_ti dot delta v_ti <= -delta p_ni c_i dot delta v_ti.
//! Normal complementarity contributes a nonpositive term. Thus
//!   E + delta p^T N^T C T W delta p <= 0, E = delta p^T W delta p.
//! Write delta p = Q Lambda^-1/2 z + Q_null a on the admitted range.
//! If Q_null^T N^T C T Q = 0, the left side is
//!   z^T [I + sym(Lambda^-1/2 Q^T N^T C T Q Lambda^1/2)] z.
//! A positive energy margin forces z=0 and hence identical W*p. This is a
//! comparison with all roots, NOT uniqueness of a selected sliding branch or
//! of a local Jacobian. Normal loads can still differ without bound in sign-
//! unconstrained algebra; the original unilateral law remains authoritative.
//!
//! Coplanar supports with coherent traction are useful instances: null normal
//! loads have zero resultant and moments, which annihilate the affine point
//! velocity field. A numerical coherent representative is only a search aid.
//! It is used only after the original W image, Coulomb law and work are checked
//! again, and still receives no unique-pressure or local-heat authority.

use super::{
    FrictionalContactConfig, FrictionalContactError, PsdVelocityCriterion, PsdVelocityWitness,
    poll, range::RangeAdmission,
};

fn bounded_direction(mut c: [f64; 2], mu: f64) -> [f64; 2] {
    let norm = c[0].hypot(c[1]);
    if norm > mu {
        // Choose a valid norm-subgradient witness, not a clipped physical
        // impulse. Its mismatch against the original impulse is checked below.
        c[0] *= mu / norm;
        c[1] *= mu / norm;
    }
    c
}

pub(super) fn direction_witness(
    range: &RangeAdmission,
    mu: &[f64],
    point: &[f64],
    velocity: &[f64],
    proposed: Option<&[[f64; 2]]>,
    config: FrictionalContactConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<PsdVelocityWitness>, FrictionalContactError> {
    let n = mu.len();
    let d = range.dimension;
    let r = range.retained.len();
    let mut directions = Vec::with_capacity(n);
    let mut traction_residual = 0.0_f64;
    let mut subgradient_velocity_residual = 0.0_f64;
    for i in 0..n {
        poll(cancelled)?;
        let a = 3 * i;
        let speed = velocity[a + 1].hypot(velocity[a + 2]);
        let c = if let Some(proposed) = proposed {
            bounded_direction(proposed[i], mu[i])
        } else if speed > config.velocity_tolerance {
            [
                mu[i] * (velocity[a + 1] / speed),
                mu[i] * (velocity[a + 2] / speed),
            ]
        } else if point[a] > 0.0 {
            bounded_direction([-point[a + 1] / point[a], -point[a + 2] / point[a]], mu[i])
        } else {
            [0.0; 2]
        };
        if c.iter().any(|v| !v.is_finite()) || point[a] < 0.0 {
            return Ok(None);
        }
        traction_residual = traction_residual
            .max((point[a + 1] + point[a] * c[0]).hypot(point[a + 2] + point[a] * c[1]));
        let subgradient = (mu[i] * speed - c[0] * velocity[a + 1] - c[1] * velocity[a + 2]).abs();
        subgradient_velocity_residual = subgradient_velocity_residual.max(subgradient);
        if mu[i] > 0.0
            && subgradient / mu[i] > 2.0 * config.velocity_tolerance + range.zero_band * speed
        {
            return Ok(None);
        }
        if speed > config.velocity_tolerance {
            // Even an unloaded sliding contact must use its actual slip
            // direction. A free disk subgradient is allowed only at sticking.
            let direction_error = (c[0] - mu[i] * velocity[a + 1] / speed)
                .hypot(c[1] - mu[i] * velocity[a + 2] / speed);
            if direction_error > range.zero_band * mu[i] {
                return Ok(None);
            }
        }
        directions.push(c);
    }
    if traction_residual > config.impulse_tolerance
        || !traction_residual.is_finite()
        || !subgradient_velocity_residual.is_finite()
    {
        return Ok(None);
    }

    let c_scale = directions
        .iter()
        .map(|c| c[0].hypot(c[1]))
        .fold(0.0, f64::max);
    if !c_scale.is_finite() {
        return Ok(None);
    }
    let normalized: Vec<_> = directions
        .iter()
        .map(|c| {
            if c_scale > 0.0 {
                [c[0] / c_scale, c[1] / c_scale]
            } else {
                [0.0; 2]
            }
        })
        .collect();
    let q = &range.vectors;
    let mut coupling = vec![0.0; d * r];
    for i in 0..d {
        poll(cancelled)?;
        for (b, &k) in range.retained.iter().enumerate() {
            for j in 0..n {
                coupling[i * r + b] += q[3 * j * d + i]
                    * (normalized[j][0] * q[(3 * j + 1) * d + k]
                        + normalized[j][1] * q[(3 * j + 2) * d + k]);
            }
        }
    }
    // Test the null identity relative to c, before Lambda weighting. Tiny
    // friction or a weak retained eigenvalue cannot hide a nonzero coupling.
    let mut null_residual = 0.0_f64;
    if c_scale > 0.0 {
        for &i in &range.null {
            for b in 0..r {
                null_residual = null_residual.max(coupling[i * r + b].abs());
            }
        }
    }
    if !null_residual.is_finite() || null_residual > range.zero_band {
        return Ok(None);
    }
    let mut energy = vec![0.0; r * r];
    for (a, &i) in range.retained.iter().enumerate() {
        poll(cancelled)?;
        for (b, &j) in range.retained.iter().enumerate() {
            let forward = c_scale
                * coupling[i * r + b]
                * (range.eigenvalues[j] / range.eigenvalues[i]).sqrt();
            let backward = c_scale
                * coupling[j * r + a]
                * (range.eigenvalues[i] / range.eigenvalues[j]).sqrt();
            energy[a * r + b] = 0.5 * forward + 0.5 * backward + if a == b { 1.0 } else { 0.0 };
        }
    }
    if energy.iter().any(|v| !v.is_finite()) {
        return Ok(None);
    }
    poll(cancelled)?;
    let (values, basis) = fs_la::eigen::jacobi_eigh(&energy, r);
    poll(cancelled)?;
    if values.iter().chain(&basis).any(|v| !v.is_finite()) {
        return Ok(None);
    }
    let scale = energy
        .chunks_exact(r)
        .map(|row| row.iter().map(|v| v.abs()).sum::<f64>())
        .fold(0.0, f64::max);
    let mut spectral_residual = 0.0_f64;
    let mut orthogonality_residual = 0.0_f64;
    for i in 0..r {
        poll(cancelled)?;
        for j in 0..r {
            let mut image = 0.0;
            let mut gram = 0.0;
            for k in 0..r {
                image += energy[i * r + k] * basis[k * r + j];
                gram += basis[k * r + i] * basis[k * r + j];
            }
            spectral_residual = spectral_residual.max((image - values[j] * basis[i * r + j]).abs());
            orthogonality_residual =
                orthogonality_residual.max((gram - if i == j { 1.0 } else { 0.0 }).abs());
        }
    }
    if spectral_residual > range.zero_band * scale || orthogonality_residual > range.zero_band {
        return Ok(None);
    }
    let margin = values[0]
        - r as f64 * (spectral_residual + scale * orthogonality_residual)
        - 512.0 * f64::EPSILON * d as f64 * scale;
    if !margin.is_finite() || margin <= 1e-8 {
        return Ok(None);
    }
    Ok(Some(range.witness(
        mu,
        (1.0 - margin).max(0.0),
        PsdVelocityCriterion::CandidateFrictionDirection {
            weighted_directions: directions,
            nullspace_coupling_residual: null_residual,
            energy_coercivity_lower_bound: margin,
            energy_spectral_residual: spectral_residual,
            energy_orthogonality_residual: orthogonality_residual,
            traction_residual,
            subgradient_velocity_residual,
        },
    )))
}

/// Try an algebraic representative with common tangential / normal traction
/// magnitude along the free-slip directions. The proposal is covariant under
/// independent changes of each tangent basis. It is deliberately not a motion
/// model: nonuniform slip may fail the original-image, law or global-direction
/// checks. Minimum-norm pressure is not a material law.
pub(super) fn coherent_representative(
    range: &RangeAdmission,
    original: &[f64],
    point: &[f64],
    free: &[f64],
    config: FrictionalContactConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<(Vec<f64>, Vec<[f64; 2]>)>, FrictionalContactError> {
    let d = range.dimension;
    let n = d / 3;
    let mut total_normal = 0.0;
    let mut total_traction = 0.0;
    let mut directions = Vec::with_capacity(n);
    for i in 0..n {
        poll(cancelled)?;
        let a = 3 * i;
        let speed = free[a + 1].hypot(free[a + 2]);
        let shear = point[a + 1].hypot(point[a + 2]);
        let direction = if speed > 0.0 {
            [free[a + 1] / speed, free[a + 2] / speed]
        } else if shear > 0.0 {
            [-point[a + 1] / shear, -point[a + 2] / shear]
        } else {
            [0.0; 2]
        };
        total_normal += point[a];
        total_traction -= point[a + 1] * direction[0] + point[a + 2] * direction[1];
        directions.push(direction);
    }
    if !total_normal.is_finite() || !total_traction.is_finite() || total_normal <= 0.0 {
        return Ok(None);
    }
    let magnitude = total_traction / total_normal;
    if !magnitude.is_finite() {
        return Ok(None);
    }
    for direction in &mut directions {
        direction[0] *= magnitude;
        direction[1] *= magnitude;
    }
    let mut matrix = vec![0.0; d * n];
    let mut image = vec![0.0; d];
    for i in 0..d {
        poll(cancelled)?;
        image[i] = fs_solver::dot(&range.normalized[i * d..(i + 1) * d], point);
        for j in 0..n {
            matrix[i * n + j] = range.normalized[i * d + 3 * j]
                - directions[j][0] * range.normalized[i * d + 3 * j + 1]
                - directions[j][1] * range.normalized[i * d + 3 * j + 2];
        }
    }
    let scale = matrix.iter().map(|v| v.abs()).fold(0.0, f64::max);
    if !scale.is_finite() || scale <= 0.0 || matrix.iter().chain(&image).any(|v| !v.is_finite()) {
        return Ok(None);
    }
    for value in &mut matrix {
        *value /= scale;
    }
    poll(cancelled)?;
    let factor = fs_la::factor::svd_jacobi(&matrix, d, n);
    poll(cancelled)?;
    if factor
        .u
        .iter()
        .chain(&factor.v)
        .chain(&factor.sigma)
        .any(|v| !v.is_finite())
    {
        return Ok(None);
    }
    let cutoff = 128.0 * f64::EPSILON * d as f64 * factor.sigma[0];
    let mut normal = vec![0.0; n];
    for k in 0..n {
        poll(cancelled)?;
        if factor.sigma[k] <= cutoff {
            continue;
        }
        let component =
            (0..d).map(|i| factor.u[i * n + k] * image[i]).sum::<f64>() / factor.sigma[k] / scale;
        for j in 0..n {
            normal[j] += factor.v[j * n + k] * component;
        }
    }
    if normal.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Ok(None);
    }
    let candidate: Vec<_> = normal
        .iter()
        .zip(&directions)
        .flat_map(|(&p, c)| [p, -c[0] * p, -c[1] * p])
        .collect();
    if candidate.iter().any(|v| !v.is_finite()) {
        return Ok(None);
    }
    // Check the ORIGINAL mobility, not the SVD's truncated factor. A proposed
    // gauge change must preserve each contact velocity to roundoff and remain
    // inside the caller's physical velocity tolerance.
    for i in 0..d {
        poll(cancelled)?;
        let mut difference = 0.0;
        let mut magnitude = 0.0;
        for j in 0..d {
            difference += original[i * d + j] * (candidate[j] - point[j]);
            magnitude += original[i * d + j].abs() * (candidate[j].abs() + point[j].abs());
        }
        let allowed = (512.0 * f64::EPSILON * d as f64 * magnitude).min(config.velocity_tolerance);
        if !difference.is_finite() || !magnitude.is_finite() || difference.abs() > allowed {
            return Ok(None);
        }
    }
    Ok(Some((candidate, directions)))
}
