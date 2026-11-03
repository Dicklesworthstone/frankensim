//! Gyrostat drift with carrier momentum retained separately from stored spin.
//!
//! k = I*omega is authoritative here. The old total-momentum API remains
//! available, but forming k+h and later subtracting h cannot retain a carrier
//! smaller than the spacing of that sum.

use super::{
    HeldTorque, LieStepError, So3, So3Tangent, Vec3, admit_rotor, finite, rigid_body_momentum_step,
};
use fs_ga::Quat;
use fs_math::{dd::Dd, det};

/// Free gyrostat drift retaining carrier momentum k = I*omega separately.
///
/// The rotor momentum is held in the body's principal frame. Each axis flow
/// uses the same exact Hamiltonian subflow as [`super::momentum_step`], with
/// k' = F^-1 k + (F^-1 - Id) h evaluated without subtracting large rotations.
/// Zero stored spin delegates to the original rigid-body map. Exactly
/// spherical/aligned axisymmetric carriers retain their exact commuting flows;
/// general carriers use the original symmetric five-axis composition.
///
/// This representation preserves small carrier motion beside a much larger
/// rotor. It does not recover information already lost in a supplied total L.
///
/// # Errors
/// Refuses the same inadmissible inertia, duration and group operations as
/// [`super::momentum_step`], or non-finite carrier momentum/transport.
pub fn momentum_step(
    mut rotation: So3,
    mut carrier: Vec3,
    inertia: Vec3,
    rotor: Vec3,
    duration: f64,
) -> Result<(So3, Vec3), LieStepError> {
    admit_rotor(rotor)?;
    if rotor == Vec3::default() {
        return rigid_body_momentum_step(rotation, carrier, inertia, duration);
    }
    if !duration.is_finite() || duration == 0.0 {
        return Err(invalid("gyrostat step size"));
    }
    if ![inertia.x, inertia.y, inertia.z]
        .into_iter()
        .all(|value| value.is_finite() && value > 0.0)
    {
        return Err(invalid("gyrostat reduced inertia"));
    }
    admit_carrier(carrier)?;
    if let Some(state) = symmetric_step(rotation, carrier, inertia, rotor, duration)? {
        return Ok(state);
    }
    let moments = [inertia.x, inertia.y, inertia.z];
    for (axis, step) in [
        (0, duration * 0.5),
        (1, duration * 0.5),
        (2, duration),
        (1, duration * 0.5),
        (0, duration * 0.5),
    ] {
        let component = [carrier.x, carrier.y, carrier.z][axis];
        let angle = step * (component / moments[axis]);
        if angle == 0.0 {
            continue;
        }
        let mut tangent = [0.0; 3];
        tangent[axis] = angle;
        let increment = So3::exp(So3Tangent::new(vec(tangent)))?;
        carrier = transport(increment, carrier, rotor)?;
        rotation = rotation.compose(increment)?;
    }
    Ok((rotation, carrier))
}

/// Change of a vector under a rotation, without rotate(value) - value.
/// q is allowed either quaternion sign; both products are sign invariant.
fn rotation_change(rotation: So3, value: Vec3) -> Vec3 {
    let q = rotation.as_quat();
    let v = Vec3::new(q.x, q.y, q.z);
    let norm_squared = q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z;
    let t = v.cross(value).scale(2.0 / norm_squared);
    t.scale(q.w) + v.cross(t)
}

fn transport(increment: So3, carrier: Vec3, rotor: Vec3) -> Result<Vec3, LieStepError> {
    let inverse = increment.inverse();
    let next = inverse.rotate(carrier)? + rotation_change(inverse, rotor);
    admit_carrier(next)?;
    Ok(next)
}

/// Half of the rotor-to-carrier transfer, retaining the endpoint subtraction
/// through the compensated kick. Rounding this difference to f64 first can
/// discard the entire carrier under a balanced external/internal impulse.
pub(super) fn half_transfer(before: Vec3, after: Vec3) -> [Dd; 3] {
    let before = [before.x, before.y, before.z];
    let after = [after.x, after.y, after.z];
    core::array::from_fn(|i| {
        Dd::from_f64(0.5) * Dd::from_f64(before[i]) - Dd::from_f64(0.5) * Dd::from_f64(after[i])
    })
}

/// Accumulate the actual external impulse and internal transfer before adding
/// them to the carrier. Opposing large terms must not erase a small k first.
pub(super) fn kick(
    carrier: Vec3,
    torque: Vec3,
    half_duration: f64,
    half_transfer: [Dd; 3],
) -> Result<Vec3, LieStepError> {
    if !finite(torque) {
        return Err(invalid("rigid-body torque"));
    }
    let k = [carrier.x, carrier.y, carrier.z];
    let t = [torque.x, torque.y, torque.z];
    let next = vec(core::array::from_fn(|i| {
        let impulse = Dd::from_f64(half_duration) * Dd::from_f64(t[i]);
        (Dd::from_f64(k[i]) + (impulse + half_transfer[i])).to_f64()
    }));
    admit_carrier(next)?;
    Ok(next)
}

pub(super) fn torque_at(torque: HeldTorque, rotation: So3) -> Result<Vec3, LieStepError> {
    match torque {
        HeldTorque::Body(value) => Ok(value),
        HeldTorque::Space(value) => Ok(rotation.inverse().rotate(value)?),
    }
}

fn symmetric_step(
    rotation: So3,
    carrier: Vec3,
    inertia: Vec3,
    rotor: Vec3,
    duration: f64,
) -> Result<Option<(So3, Vec3)>, LieStepError> {
    if carrier == Vec3::default() {
        return Ok(Some((rotation, carrier)));
    }
    let moments = [inertia.x, inertia.y, inertia.z];
    let stored = [rotor.x, rotor.y, rotor.z];
    let k = [carrier.x, carrier.y, carrier.z];
    for axis in 0..3 {
        let a = (axis + 1) % 3;
        let b = (axis + 2) % 3;
        if moments[a] == moments[b] && stored[a] == 0.0 && stored[b] == 0.0 {
            if k[a] == 0.0 && k[b] == 0.0 {
                // Exact axial motion depends only on the carrier rate, not on
                // the difference between two very large rotor phases.
                let mut tangent = [0.0; 3];
                tangent[axis] = duration * (k[axis] / moments[axis]);
                let increment = So3::exp(So3Tangent::new(vec(tangent)))?;
                return Ok(Some((rotation.compose(increment)?, carrier)));
            }
            let mut direction = [0.0; 3];
            direction[axis] = stored[axis].signum();
            return commuting_step(
                rotation,
                carrier,
                rotor,
                vec(direction),
                stored[axis].abs(),
                moments[a],
                moments[axis],
                duration,
            )
            .map(Some);
        }
    }
    if moments[0] == moments[1] && moments[1] == moments[2] {
        let magnitude = scaled_norm(rotor);
        if !(magnitude.is_finite() && magnitude > 0.0) {
            return Err(invalid("gyrostat stored momentum magnitude"));
        }
        let direction = Vec3::new(
            rotor.x / magnitude,
            rotor.y / magnitude,
            rotor.z / magnitude,
        );
        return commuting_step(
            rotation, carrier, rotor, direction, magnitude, moments[0], moments[0], duration,
        )
        .map(Some);
    }
    Ok(None)
}

// Exactly H = |L|^2/(2A) + H_residual, with rotor along n and axial
// inertia C. For the spherical case A=C and n may be any rotor direction.
#[allow(clippy::too_many_arguments)]
fn commuting_step(
    rotation: So3,
    carrier: Vec3,
    rotor: Vec3,
    direction: Vec3,
    stored_magnitude: f64,
    transverse_inertia: f64,
    axial_inertia: f64,
    duration: f64,
) -> Result<(So3, Vec3), LieStepError> {
    let a = transverse_inertia;
    let c = axial_inertia;
    let parallel = carrier.dot(direction);
    let perpendicular = carrier - direction.scale(parallel);
    let total_parallel = stored_magnitude + parallel;
    let transverse = scaled_norm(perpendicular);
    let total_magnitude = det::hypot(total_parallel, transverse);
    if !total_magnitude.is_finite() {
        return Err(invalid("gyrostat symmetric momentum magnitude"));
    }

    // The carrier precesses about n at h/A + k_parallel*(1/A-1/C).
    // Keep those phases separate: the second can matter even if adding it
    // to the first would round it away.
    let spin_phase = So3::exp(So3Tangent::new(Vec3::new(
        duration * (rotor.x / a),
        duration * (rotor.y / a),
        duration * (rotor.z / a),
    )))?;
    let differential_rate = if a <= c {
        (parallel / a) * ((c - a) / c)
    } else {
        (parallel / c) * ((c - a) / a)
    };
    let precession = if differential_rate == 0.0 {
        spin_phase
    } else {
        spin_phase.compose(So3::exp(So3Tangent::new(
            direction.scale(duration * differential_rate),
        ))?)?
    };
    // h is invariant under this exact residual rotation. Rotating k itself
    // avoids amplifying a tiny attitude-roundoff error by a huge h.
    let next_carrier = precession.rotate(carrier)?;
    admit_carrier(next_carrier)?;

    let increment = if total_parallel <= 0.0 {
        // No small-axis subtraction in the antiparallel/zero-total case.
        // Retain the ordinary exact product, including L=0.
        let total = carrier + rotor;
        let spatial = So3::exp(So3Tangent::new(Vec3::new(
            duration * (total.x / a),
            duration * (total.y / a),
            duration * (total.z / a),
        )))?;
        spatial.compose(precession.inverse())?
    } else {
        // |L|-L_parallel = |k_perp|^2/(|L|+L_parallel), evaluated without
        // squaring large values or subtracting nearly equal magnitudes.
        let transverse_correction =
            (transverse / total_magnitude) * transverse / (1.0 + total_parallel / total_magnitude);
        let d_parallel = -transverse_correction / total_magnitude;
        let axis_difference = Vec3::new(
            perpendicular.x / total_magnitude,
            perpendicular.y / total_magnitude,
            perpendicular.z / total_magnitude,
        ) + direction.scale(d_parallel);
        let cross_difference = Vec3::new(
            perpendicular.x / total_magnitude,
            perpendicular.y / total_magnitude,
            perpendicular.z / total_magnitude,
        )
        .cross(direction);
        // a_phase-b_phase, never obtained by subtracting large phases.
        let delta = 0.5 * duration * (transverse_correction / a + parallel / c);
        let sin_delta = det::sin(delta);
        let cos_delta = det::cos(delta);
        let qb = precession.as_quat();
        let sin_b = Vec3::new(qb.x, qb.y, qb.z).dot(direction);
        let cos_b = qb.w;
        let sin_a = sin_b * cos_delta + cos_b * sin_delta;
        let ss = sin_a * sin_b;
        let sc = sin_a * cos_b;
        let vector =
            direction.scale(sin_delta) + axis_difference.scale(sc) - cross_difference.scale(ss);
        So3::try_from_quat(Quat {
            w: cos_delta + ss * d_parallel,
            x: vector.x,
            y: vector.y,
            z: vector.z,
        })?
    };
    Ok((rotation.compose(increment)?, next_carrier))
}

fn scaled_norm(value: Vec3) -> f64 {
    let scale = value.x.abs().max(value.y.abs()).max(value.z.abs());
    if scale == 0.0 {
        return 0.0;
    }
    let x = value.x / scale;
    let y = value.y / scale;
    let z = value.z / scale;
    scale * det::sqrt(x * x + y * y + z * z)
}

fn vec(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}

fn admit_carrier(value: Vec3) -> Result<(), LieStepError> {
    if finite(value) {
        Ok(())
    } else {
        Err(invalid("gyrostat carrier momentum"))
    }
}

fn invalid(context: &'static str) -> LieStepError {
    LieStepError::InvalidParameter { context }
}
