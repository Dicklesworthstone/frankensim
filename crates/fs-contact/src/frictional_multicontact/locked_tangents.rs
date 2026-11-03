//! Exact Coulomb identities for structurally locked tangent pairs.
//!
//! A zero tangent row and column of the ORIGINAL active mobility fixes that
//! velocity independently of every impulse. At nonzero fixed slip the disk
//! law therefore requires p_t = -mu*p_n*v_t/|v_t|. At zero slip it leaves the
//! whole static disk available; selecting a representative never identifies
//! its local traction. No small mobility is discarded or replaced here.

use super::{
    FrictionalContactConfig, FrictionalContactError, LockedTangentWorkAdmission,
    LockedTangentWorkSubstitution, PsdImpulseWorkWitness, poll, range::RangeAdmission,
};
use fs_math::dd::Dd;

/// These new traction/work identities require complete tangent-pair locks.
/// The original range owner independently admits single constrained axes when
/// its existing global motion and work proofs apply; that scope is unchanged.
pub(super) fn has_locked_pair(matrix: &[f64], dimension: usize) -> bool {
    dimension > 0
        && dimension % 3 == 0
        && (0..dimension / 3).any(|contact| is_locked(matrix, dimension, contact))
}

fn is_locked(matrix: &[f64], dimension: usize, contact: usize) -> bool {
    let tangent = 3 * contact + 1;
    if matrix[tangent * dimension + tangent] != 0.0
        || matrix[(tangent + 1) * dimension + tangent + 1] != 0.0
    {
        return false;
    }
    (tangent..=tangent + 1).all(|row| {
        (0..dimension).all(|column| {
            matrix[row * dimension + column] == 0.0 && matrix[column * dimension + row] == 0.0
        })
    })
}

/// Close the exact slip law before the unchanged original-matrix audit. Even
/// nonzero slip below the numerical velocity tolerance has its actual traction
/// direction. Changing these impulse coordinates cannot change ANY active
/// velocity, because their complete original columns are zero.
pub(super) fn complete_tractions(
    matrix: &[f64],
    velocity: &[f64],
    mu: &[f64],
    point: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), FrictionalContactError> {
    let dimension = velocity.len();
    for (contact, &coefficient) in mu.iter().enumerate() {
        poll(cancelled)?;
        if !is_locked(matrix, dimension, contact) {
            continue;
        }
        let normal = 3 * contact;
        let speed = velocity[normal + 1].hypot(velocity[normal + 2]);
        // The original natural map uses this same nonnegative disk radius;
        // the original normal impulse itself is never changed or clipped.
        let radius = coefficient * point[normal].max(0.0);
        if !speed.is_finite() || !radius.is_finite() {
            return Err(FrictionalContactError::NonFinite);
        }
        if speed > 0.0 {
            for axis in 1..=2 {
                point[normal + axis] = -radius * (velocity[normal + axis] / speed);
            }
        }
    }
    poll(cancelled)
}

/// Work differences between two law roots are v_before dot delta_p. For a
/// locked sliding pair the exact traction identity replaces its tangent terms
/// by -mu*|v_t|*delta_p_n. This changes only the work covector used for the
/// gauge proof; the physical mobility, impulses and velocities are untouched.
fn work_covector(
    matrix: &[f64],
    velocity: &[f64],
    mu: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<Dd>, Vec<LockedTangentWorkSubstitution>), FrictionalContactError> {
    let dimension = velocity.len();
    let mut covector: Vec<_> = velocity.iter().copied().map(Dd::from_f64).collect();
    let mut substitutions = Vec::new();
    for (contact, &coefficient) in mu.iter().enumerate() {
        poll(cancelled)?;
        if coefficient == 0.0 || !is_locked(matrix, dimension, contact) {
            continue;
        }
        let normal = 3 * contact;
        let speed = velocity[normal + 1].hypot(velocity[normal + 2]);
        if !speed.is_finite() {
            return Err(FrictionalContactError::NonFinite);
        }
        if speed == 0.0 {
            continue;
        }
        let loss_per_normal = Dd::from_f64(coefficient) * Dd::from_f64(speed);
        if !finite_dd(loss_per_normal) || loss_per_normal == Dd::ZERO {
            return Err(FrictionalContactError::Invalid(
                "representable locked-tangent work covector required",
            ));
        }
        let direction = [
            coefficient * (velocity[normal + 1] / speed),
            coefficient * (velocity[normal + 2] / speed),
        ];
        covector[normal] = covector[normal] - loss_per_normal;
        covector[normal + 1] = Dd::ZERO;
        covector[normal + 2] = Dd::ZERO;
        if !finite_dd(covector[normal]) || direction.iter().any(|v| !v.is_finite()) {
            return Err(FrictionalContactError::NonFinite);
        }
        substitutions.push(LockedTangentWorkSubstitution {
            active_contact: contact,
            weighted_slip_direction: direction,
            normal_velocity_shift: -loss_per_normal.hi,
            normal_velocity_shift_low: -loss_per_normal.lo,
        });
    }
    Ok((covector, substitutions))
}

fn finite_dd(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}

/// Work-specific admission uses the ORIGINAL W, never an absolute velocity
/// tolerance. A cached spectral range supplies only a candidate row-space
/// dual; its actual original-matrix residual is recomputed in DD.
pub(super) fn work_witness(
    range: &RangeAdmission,
    matrix: &[f64],
    velocity: &[f64],
    mu: &[f64],
    point: &[f64],
    config: FrictionalContactConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<PsdImpulseWorkWitness>, FrictionalContactError> {
    if !has_locked_pair(matrix, velocity.len()) {
        return Ok(None);
    }
    let (mut covector, substitutions) = work_covector(matrix, velocity, mu, cancelled)?;
    if substitutions.is_empty() {
        return Ok(None);
    }
    let d = range.dimension;
    let mut zero_friction_coordinates = Vec::new();
    for (contact, &coefficient) in mu.iter().enumerate() {
        if coefficient == 0.0 {
            for coordinate in [3 * contact + 1, 3 * contact + 2] {
                covector[coordinate] = Dd::ZERO;
                zero_friction_coordinates.push(coordinate);
            }
        }
    }
    let mut dual = vec![Dd::ZERO; d];
    for &mode in &range.retained {
        poll(cancelled)?;
        let mut component = Dd::ZERO;
        for row in 0..d {
            component = component + covector[row] * Dd::from_f64(range.vectors[row * d + mode]);
        }
        let coefficient =
            component / Dd::from_f64(range.eigenvalues[mode]) / Dd::from_f64(range.matrix_scale);
        for row in 0..d {
            dual[row] = dual[row] + coefficient * Dd::from_f64(range.vectors[row * d + mode]);
        }
    }
    let mut residual = 0.0_f64;
    let mut arithmetic_scale = 0.0_f64;
    for column in 0..d {
        poll(cancelled)?;
        let mut image = Dd::ZERO;
        let mut scale = covector[column].abs();
        for row in 0..d {
            let term = Dd::from_f64(matrix[row * d + column]) * dual[row];
            image = image + term;
            scale = scale + term.abs();
        }
        let difference = covector[column] - image;
        if !finite_dd(difference) || !finite_dd(scale) {
            return Err(FrictionalContactError::NonFinite);
        }
        residual = residual.hypot(difference.hi.abs() + difference.lo.abs());
        arithmetic_scale = arithmetic_scale.hypot(scale.to_f64());
    }
    // DD signed products protect small friction coefficients beside a large
    // incoming normal speed. Retain a conservative DD arithmetic allowance;
    // unlike a velocity tolerance it is converted to joules using EVERY root.
    residual = residual * (1.0 + 64.0 * f64::EPSILON * d as f64)
        + 256.0 * f64::EPSILON * f64::EPSILON * d as f64 * arithmetic_scale;

    let mut work_scale = 0.0;
    for row in 0..d {
        poll(cancelled)?;
        let mut after = velocity[row];
        for column in 0..d {
            after += matrix[row * d + column] * point[column];
        }
        work_scale += (point[row] * (0.5 * velocity[row] + 0.5 * after)).abs();
    }
    let allowed_work = config.energy_tolerance + 64.0 * f64::EPSILON * work_scale;
    if !residual.is_finite() || !allowed_work.is_finite() {
        return Err(FrictionalContactError::NonFinite);
    }
    let mut impulse_bound = cone_bound(matrix, point, point, mu, cancelled)?;
    let normal_probe: Vec<_> = (0..d)
        .map(|coordinate| if coordinate % 3 == 0 { 1.0 } else { 0.0 })
        .collect();
    if let Some(bound) = cone_bound(matrix, point, &normal_probe, mu, cancelled)? {
        impulse_bound = Some(impulse_bound.map_or(bound, |current| current.min(bound)));
    }
    let (maximum_work_difference, allowed_velocity) = match impulse_bound {
        Some(bound) if bound > 0.0 => (
            2.0 * (bound * residual),
            (0.5 * (allowed_work / bound)).min(f64::MAX),
        ),
        Some(_) => (0.0, residual),
        None if residual == 0.0 => (0.0, 0.0),
        None => {
            return Err(FrictionalContactError::UnresolvedWork {
                range_velocity_residual: residual,
                allowed: 0.0,
            });
        }
    };
    if maximum_work_difference > allowed_work {
        return Err(FrictionalContactError::UnresolvedWork {
            range_velocity_residual: residual,
            allowed: allowed_velocity,
        });
    }
    poll(cancelled)?;
    Ok(Some(PsdImpulseWorkWitness {
        impulse_coordinates: (0..d).collect(),
        rank: range.retained.len(),
        locked_tangent_substitutions: substitutions,
        locked_tangent_work_admission: Some(LockedTangentWorkAdmission {
            zero_friction_coordinates,
            impulse_norm_bound: impulse_bound,
            maximum_work_difference,
            allowed_work_difference: allowed_work,
        }),
        range_velocity_residual: residual,
        allowed_range_velocity_residual: allowed_velocity,
        spectral_residual: range.spectral_residual,
        orthogonality_residual: range.orthogonality_residual,
    }))
}

/// For EVERY same-image root p, (W^T*x) dot p = x dot (W*p0) = D.
/// With g_k=(W^T*x)_n-mu_k*|(W^T*x)_t|>0 and |p_t|<=mu*p_n,
/// sum p_n <= D/min(g), hence |p| <= hypot(1,max(mu))*D/min(g).
/// The candidate norm alone is never used as a bound on another root.
fn cone_bound(
    matrix: &[f64],
    point: &[f64],
    probe: &[f64],
    mu: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<f64>, FrictionalContactError> {
    let d = point.len();
    let mut dual = vec![Dd::ZERO; d];
    let mut dual_error = vec![0.0; d];
    for column in 0..d {
        poll(cancelled)?;
        let mut scale = Dd::ZERO;
        for row in 0..d {
            let term = Dd::from_f64(matrix[row * d + column]) * Dd::from_f64(probe[row]);
            dual[column] = dual[column] + term;
            scale = scale + term.abs();
        }
        // Cancellation can leave a tiny dual component after huge signed
        // products. Its arithmetic bound belongs to the ORIGINAL terms.
        dual_error[column] =
            (Dd::from_f64(256.0 * f64::EPSILON * f64::EPSILON * d as f64) * scale).to_f64();
        if !finite_dd(scale) || !dual_error[column].is_finite() {
            return Ok(None);
        }
    }
    if dual.iter().any(|value| !finite_dd(*value)) {
        return Ok(None);
    }
    let roundoff = 128.0 * f64::EPSILON * d as f64;
    let mut minimum_gap = f64::INFINITY;
    for (contact, &coefficient) in mu.iter().enumerate() {
        poll(cancelled)?;
        let normal = 3 * contact;
        let tangent = (dual[normal + 1].hi.abs()
            + dual[normal + 1].lo.abs()
            + dual_error[normal + 1])
            .hypot(dual[normal + 2].hi.abs() + dual[normal + 2].lo.abs() + dual_error[normal + 2]);
        let friction = Dd::from_f64(coefficient) * Dd::from_f64(tangent);
        let gap = (dual[normal] - friction).to_f64()
            - dual_error[normal]
            - roundoff * (dual[normal].abs() + friction.abs()).to_f64();
        if !gap.is_finite() || gap <= 0.0 {
            return Ok(None);
        }
        minimum_gap = minimum_gap.min(gap);
    }
    let mut pairing = Dd::ZERO;
    let mut scale = Dd::ZERO;
    let mut pairing_error = Dd::ZERO;
    for (index, (value, &impulse)) in dual.iter().zip(point).enumerate() {
        let term = *value * Dd::from_f64(impulse);
        pairing = pairing + term;
        scale = scale + term.abs();
        pairing_error =
            pairing_error + Dd::from_f64(dual_error[index]) * Dd::from_f64(impulse.abs());
    }
    let upper_pairing = pairing + Dd::from_f64(roundoff) * scale + pairing_error;
    if !finite_dd(upper_pairing) || upper_pairing.lt(Dd::ZERO) {
        return Ok(None);
    }
    let coefficient = mu.iter().copied().fold(0.0, f64::max);
    let bound = ((upper_pairing / Dd::from_f64(minimum_gap))
        * Dd::from_f64(1.0_f64.hypot(coefficient)))
    .to_f64()
        * (1.0 + roundoff);
    Ok(
        (bound.is_finite() && (bound > 0.0 || point.iter().all(|impulse| *impulse == 0.0)))
            .then_some(bound),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frictional_multicontact::range;

    #[test]
    fn exact_root_retains_sub_ulp_locked_friction_in_work_admission() {
        let mobility = 2.0_f64.powi(-50);
        let coefficient = 0.1;
        let slip = 1e-16;
        let mut matrix = [0.0; 36];
        for row in [0, 3] {
            for column in [0, 3] {
                matrix[6 * row + column] = mobility;
            }
        }
        let velocity = [-1.0, slip, 0.0, -1.0, 0.0, 0.0];
        let mu = [coefficient; 2];
        let pressure = 1.0 / mobility;
        let point = [pressure, -coefficient * pressure, 0.0, pressure, 0.0, 0.0];
        // This analytic root avoids an enormous-pressure natural-map rounding
        // plateau. Both normal speeds become exactly +1 (elastic targets),
        // both fixed tangents obey their original Coulomb disks, and its
        // strictly positive loss is coefficient*slip*pressure. Redistributing
        // normal pressure from one contact to the other leaves W*p unchanged
        // but changes that loss continuously from zero to the value below.
        assert_eq!(-1.0 + mobility * (point[0] + point[3]), 1.0);
        assert_eq!(point[1], -coefficient * point[0]);
        assert!(coefficient * slip * (2.0 * pressure) > 1e-4);
        assert_eq!(-1.0 - coefficient * slip, -1.0);

        let (covector, substitutions) =
            work_covector(&matrix, &velocity, &mu, &mut || false).unwrap();
        assert_eq!(substitutions.len(), 1);
        assert_eq!(covector[0].hi, -1.0);
        assert!(covector[0].lo < 0.0);
        assert_eq!(
            (covector[0] - covector[3]).to_f64(),
            -(Dd::from_f64(coefficient) * Dd::from_f64(slip)).to_f64()
        );
        let config = FrictionalContactConfig::default();
        let admitted = range::admit(&matrix, 6, config.matrix, &mut || false)
            .unwrap()
            .unwrap();
        let result = work_witness(
            &admitted,
            &matrix,
            &velocity,
            &mu,
            &point,
            config,
            &mut || false,
        );
        assert!(
            matches!(result, Err(FrictionalContactError::UnresolvedWork { .. })),
            "exact passive root lost its sub-ULP physical work gauge: {result:?}"
        );
    }
}
