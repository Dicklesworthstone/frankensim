//! Endpoint force evaluation for time- and configuration-dependent torques.

use super::{LieStepError, compensated, gyrostat::HeldTorque, held};
use fs_ga::{So3, Vec3};
use fs_math::dd::Dd;

/// Second-order rigid rotation under an orientation- and time-dependent torque.
///
/// `rotation` maps the body's principal axes to space; `omega_body` and the
/// positive principal moments `inertia` use that body frame. `torque_at` is
/// evaluated at `(time, rotation)` before the first momentum half-kick and
/// at `(time + h, next_rotation)` before the second. Each returned
/// [`HeldTorque`] specifies the frame of that *sample*, not a load frozen
/// for the step. Torques act about the center of mass.
///
/// This supports gravity-gradient, torsional-potential, magnetic-alignment,
/// and prescribed time-varying loads without freezing their dependence on
/// attitude. The symmetric kicks enclose the existing compensated free
/// rigid-body drift. For smooth configuration/time forces the map is second
/// order; conservative autonomous forces additionally retain the symmetric
/// Hamiltonian splitting. Spatial angular momentum changes by the
/// trapezoidal sum of the two **spatial** endpoint torques to roundoff.
/// Energy is not exactly conserved.
///
/// Exactly two callback evaluations are made on success. Impulses, frame
/// transforms and momentum products retain DD tails until the returned
/// angular velocity rounds to f64. There is no allocation, finite-difference
/// force estimate, or nonlinear solve. Signed steps support reversal when
/// the caller supplies the same deterministic force law and endpoint time.
///
/// The force law must depend only on the supplied attitude and time (and
/// immutable parameters). Velocity-dependent drag, moving mass, rotor
/// coupling, discontinuities and contact events need their own integration
/// treatment; capturing a stale velocity or force in the closure does not
/// produce a second-order solution of those models. This API does not
/// differentiate the callback or supply a discrete adjoint.
///
/// # Errors
/// Refuses nonfinite time, nonfinite/zero duration, an unrepresentable or
/// unchanged endpoint time, nonpositive/nonfinite inertia, and nonfinite
/// angular velocity or torque. Arithmetic and `So3` refusals propagate.
/// A callback can refuse either stage with [`LieStepError`]. Input state is
/// passed by value and remains unchanged on every refusal; callback side
/// effects, if any, are the caller's responsibility.
pub fn rigid_body_torque_function_step<F>(
    rotation: So3,
    omega_body: Vec3,
    inertia: Vec3,
    time: f64,
    h: f64,
    torque_at: F,
) -> Result<(So3, Vec3), LieStepError>
where
    F: FnMut(f64, So3) -> Result<HeldTorque, LieStepError>,
{
    let end_time = endpoint_time(time, h)?;
    let momentum = held::momentum_from_velocity(omega_body, inertia)?;
    let (next_rotation, next_momentum) = step(
        rotation,
        momentum,
        time,
        end_time,
        h,
        torque_at,
        |orientation, kicked| {
            compensated::rigid_body_momentum_step(orientation, kicked, inertia, h)
        },
    )?;
    Ok((
        next_rotation,
        held::velocity_from_momentum(next_momentum, inertia)?,
    ))
}

/// The same force evaluation in a cached, nonprincipal material frame.
/// The tensor owner has already admitted the positive principal moments.
pub(super) fn step_in_frame<F>(
    rotation: So3,
    omega_body: Vec3,
    inertia: Vec3,
    principal_to_body: So3,
    time: f64,
    h: f64,
    torque_at: F,
) -> Result<(So3, Vec3), LieStepError>
where
    F: FnMut(f64, So3) -> Result<HeldTorque, LieStepError>,
{
    let end_time = endpoint_time(time, h)?;
    let omega = held::finite_vector(omega_body, "rigid-body angular velocity")?;
    let omega_principal = compensated::spatial_vector_to_body(principal_to_body, omega)?;
    let moments = [inertia.x, inertia.y, inertia.z];
    let mut momentum_principal = [Dd::ZERO; 3];
    for d in 0..3 {
        momentum_principal[d] = compensated::checked(
            omega_principal[d] * Dd::from_f64(moments[d]),
            "rigid-body angular momentum",
        )?;
    }
    let momentum =
        compensated::spatial_vector_to_body(principal_to_body.inverse(), momentum_principal)?;
    let (next_rotation, next_momentum) = step(
        rotation,
        momentum,
        time,
        end_time,
        h,
        torque_at,
        |orientation, kicked| {
            compensated::rigid_body_momentum_step_in_frame(
                orientation,
                kicked,
                inertia,
                principal_to_body,
                h,
            )
        },
    )?;
    let next_principal = compensated::spatial_vector_to_body(principal_to_body, next_momentum)?;
    let mut next_omega_principal = [Dd::ZERO; 3];
    for d in 0..3 {
        next_omega_principal[d] = compensated::divide(next_principal[d], Dd::from_f64(moments[d]))?;
    }
    let next_omega =
        compensated::spatial_vector_to_body(principal_to_body.inverse(), next_omega_principal)?;
    let [x, y, z] = next_omega.map(Dd::to_f64);
    Ok((next_rotation, Vec3::new(x, y, z)))
}

fn endpoint_time(time: f64, h: f64) -> Result<f64, LieStepError> {
    if !time.is_finite() || !h.is_finite() || h == 0.0 {
        return Err(LieStepError::InvalidParameter {
            context: "configuration torque time/duration",
        });
    }
    let end_time = time + h;
    if !end_time.is_finite() || end_time == time {
        return Err(LieStepError::InvalidParameter {
            context: "configuration torque endpoint time",
        });
    }
    Ok(end_time)
}

fn torque_body(orientation: So3, torque: HeldTorque) -> Result<[Dd; 3], LieStepError> {
    match torque {
        HeldTorque::Body(value) => held::finite_vector(value, "rigid-body torque"),
        HeldTorque::Space(value) => compensated::spatial_vector_to_body(
            orientation,
            held::finite_vector(value, "rigid-body spatial torque")?,
        ),
    }
}

fn step<F, D>(
    rotation: So3,
    momentum: [Dd; 3],
    time: f64,
    end_time: f64,
    h: f64,
    mut torque_at: F,
    drift: D,
) -> Result<(So3, [Dd; 3]), LieStepError>
where
    F: FnMut(f64, So3) -> Result<HeldTorque, LieStepError>,
    D: FnOnce(So3, [Dd; 3]) -> Result<(So3, [Dd; 3]), LieStepError>,
{
    let start_torque = torque_body(rotation, torque_at(time, rotation)?)?;
    let kicked = held::half_kick(momentum, start_torque, h)?;
    let (next_rotation, drifted) = drift(rotation, kicked)?;
    let end_torque = torque_body(next_rotation, torque_at(end_time, next_rotation)?)?;
    let accepted = held::half_kick(drifted, end_torque, h)?;
    Ok((next_rotation, accepted))
}
