//! Instantaneous capture of a point parcel by an original native rigid body.
//!
//! The parcel has translational momentum and zero intrinsic inertia. It sticks
//! at the supplied material point. The new COM shifts, while body axes and all
//! material points retain their physical meaning. This is not a finite-volume
//! fluid, slosh, contact-discovery or finite-duration capture model.

use super::{
    AssemblyBalance, AssemblyControls, AssemblyError, ENTRIES, Sum, VectorSum, array, finite, poll,
    vector,
};
use crate::{MassProperties, Pose, RigidBodyState, Vec3};

/// Incoming mass at an explicitly identified material point, in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointMassCapture {
    /// Nonnegative incoming mass [kg]. Zero is an admitted identity event.
    pub mass_kg: f64,
    /// Capture point relative to the original COM, in original body axes [m].
    pub attachment_arm_body_m: Vec3,
    /// Incoming parcel velocity in the body's world inertial frame [m/s].
    pub velocity_world_m_s: Vec3,
}

/// One atomic capture; the original state and mass properties are unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointMassCaptureResult {
    /// Native endpoint at the new COM, with the original body orientation.
    pub state: RigidBodyState,
    /// Full central inertia and total mass, in the original body axes.
    pub properties: MassProperties,
    /// New COM relative to the old COM [m]. Subtract this from every original
    /// body point before continuing with the returned native state.
    pub center_shift_body_m: Vec3,
    /// Internal impulse received by the original body [N s]. The incoming
    /// parcel receives its negative; this is not the composite momentum gain.
    pub captured_impulse_world_ns: Vec3,
    /// Nonnegative mechanical energy removed by sticking [J]. Its thermal
    /// destination is caller-owned; no energy is silently credited to a fluid.
    pub removed_relative_energy_j: f64,
    /// Actual native momentum and polarized kinetic-energy balance residuals.
    pub balance: AssemblyBalance,
}

/// Capture a point parcel without inventing its radius or rotational inertia.
///
/// For old mass M, incoming mass m and arm r, the reduced mass is
/// mu=M*m/(M+m), the COM shift is m*r/(M+m), and the new central inertia is
/// I + mu*(|r|^2*1-r*r^T). Original body axes are preserved. Incoming momentum
/// is included once; no external impulse is applied.
///
/// The internal capture impulse solves the original point mobility, including
/// full rotational recoil. Its positive quadratic work gives the sticking
/// loss without subtracting large kinetic energies. The balance gate uses
/// actual native endpoints, a translating frame and polarized spin changes.
///
/// Positive capture consumes two part slots and one internal impulse slot in
/// the existing assembly controls. Zero mass consumes one part and no impulse.
/// All inputs remain borrowed; cancellation or refusal publishes no endpoint.
pub fn capture_point_mass(
    body: &RigidBodyState,
    properties: MassProperties,
    incoming: PointMassCapture,
    controls: AssemblyControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PointMassCaptureResult, AssemblyError> {
    poll(&mut cancelled)?;
    properties.validate()?;
    if !incoming.mass_kg.is_finite() || incoming.mass_kg < 0.0 {
        return Err(AssemblyError::Invalid(
            "finite nonnegative captured mass required",
        ));
    }
    vector(incoming.attachment_arm_body_m)?;
    vector(incoming.velocity_world_m_s)?;
    let active = usize::from(incoming.mass_kg > 0.0);
    controls.validate(1 + active, active)?;
    let velocity = body.center_of_mass_velocity_world(properties)?;
    let spin = body.angular_momentum_body();
    let omega = properties.angular_velocity_body_checked(spin)?;
    let zero_balance = AssemblyBalance {
        linear_momentum_residual_world_ns: Vec3::ZERO,
        angular_momentum_residual_world_nms: Vec3::ZERO,
        energy_residual_j: 0.0,
    };
    if active == 0 {
        poll(&mut cancelled)?;
        return Ok(PointMassCaptureResult {
            state: *body,
            properties,
            center_shift_body_m: Vec3::ZERO,
            captured_impulse_world_ns: Vec3::ZERO,
            removed_relative_energy_j: 0.0,
            balance: zero_balance,
        });
    }
    let mass = properties.mass();
    let added = incoming.mass_kg;
    let total = finite(mass + added)?;
    if total <= mass || total <= added {
        return Err(AssemblyError::Invalid(
            "both masses must remain represented in the total",
        ));
    }
    let alpha = added / total;
    let beta = mass / total;
    let reduced = finite(if mass <= added {
        mass * alpha
    } else {
        added * beta
    })?;
    if reduced <= 0.0 {
        return Err(AssemblyError::Invalid(
            "positive reduced capture mass required",
        ));
    }
    let root_mass = reduced.sqrt();
    let arm = incoming.attachment_arm_body_m;
    let scaled_arm = vector(arm.scale(root_mass))?;
    let orientation = body.pose().orientation();
    let relative_world = vector(incoming.velocity_world_m_s.sub(velocity))?;
    let relative = vector(orientation.rotate_world_to_body(relative_world))?;
    let slip = vector(relative.sub(omega.cross(arm)))?;

    // A = 1 + mu*C^T*C, with C^T*C = -[r]x I_old^-1 [r]x.
    // Reuse the native inverse-inertia factor and bounded 3-by-3 solver.
    let axes = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];
    let factors = axes.map(|axis| properties.angular_mobility_factor(scaled_arm.cross(axis)));
    let mut mobility = [[0.0; 3]; 3];
    for (row, col) in ENTRIES {
        let value = finite(factors[row].dot(factors[col]) + if row == col { 1.0 } else { 0.0 })?;
        mobility[row][col] = value;
        mobility[col][row] = value;
    }
    poll(&mut cancelled)?;
    let solution = crate::solve_rotation_3(mobility, array(slip))?;
    let impulse_over_root =
        vector(Vec3::new(solution[0], solution[1], solution[2]).scale(root_mass))?;
    let impulse = vector(impulse_over_root.scale(root_mass))?;
    let rotational_response =
        vector(properties.angular_mobility_factor(scaled_arm.cross(impulse_over_root)))?;
    let loss = finite(
        0.5 * (impulse_over_root.dot(impulse_over_root)
            + rotational_response.dot(rotational_response)),
    )?;
    let capture_work_defect = finite(0.5 * slip.dot(impulse) - loss)?;
    if loss < 0.0 || capture_work_defect.abs() > controls.energy_tolerance_j {
        return Err(AssemblyError::Invalid(
            "unresolved point-capture mobility work",
        ));
    }

    let center_shift = vector(arm.scale(alpha))?;
    let next_properties = if arm == Vec3::ZERO {
        MassProperties::from_principal_frame(
            total,
            Vec3::ZERO,
            properties.principal_inertia_body(),
            properties.principal_axes_to_body(),
        )?
    } else {
        let mut inertia = properties.inertia_tensor_body();
        let r = array(scaled_arm);
        for (row, col) in ENTRIES {
            let shift = if row == col {
                r[(row + 1) % 3].powi(2) + r[(row + 2) % 3].powi(2)
            } else {
                -r[row] * r[col]
            };
            let value = finite(inertia.m[3 * row + col] + shift)?;
            inertia.m[3 * row + col] = value;
            inertia.m[3 * col + row] = value;
        }
        MassProperties::from_inertia_tensor(total, Vec3::ZERO, inertia)?
    };
    let imported = vector(incoming.velocity_world_m_s.scale(added))?;
    let mut linear = VectorSum::default();
    linear.add(body.linear_momentum_world());
    linear.add(imported);
    let next_spin = vector(spin.add(scaled_arm.cross(relative.scale(root_mass))))?;
    let next = RigidBodyState::new(
        Pose::new(
            body.pose().point_world_from_body(center_shift)?,
            orientation,
        )?,
        vector(linear.value())?,
        next_spin,
    )?;
    let next_omega = next_properties.angular_velocity_body_checked(next_spin)?;
    next.center_of_mass_velocity_world(next_properties)?;

    // Recover actual native residuals before publication. Removing the common
    // translation before cross products keeps small spin and capture work.
    let before_relative = relative_momentum(body.linear_momentum_world(), mass, velocity)?;
    let after_relative = relative_momentum(next.linear_momentum_world(), total, velocity)?;
    let imported_relative = vector(relative_world.scale(added))?;
    let world_arm = vector(orientation.rotate_body_to_world(arm))?;
    let actual_shift = vector(
        next.pose()
            .position_world()
            .sub(body.pose().position_world()),
    )?;
    let first_moment_error = vector(actual_shift.scale(total).sub(world_arm.scale(added)))?;
    let mut linear_error = VectorSum::default();
    linear_error.add(next.linear_momentum_world());
    linear_error.add(body.linear_momentum_world().scale(-1.0));
    linear_error.add(imported.scale(-1.0));
    let mut angular_error = VectorSum::default();
    angular_error.add(orientation.rotate_body_to_world(next_spin.sub(spin)));
    angular_error.add(actual_shift.cross(after_relative));
    angular_error.add(world_arm.cross(imported_relative).scale(-1.0));
    angular_error.add(first_moment_error.cross(velocity));
    let mut energy_error = Sum::default();
    energy_error.add(0.5 * after_relative.dot(after_relative.scale(total.recip())));
    energy_error.add(-0.5 * before_relative.dot(before_relative.scale(mass.recip())));
    energy_error.add(-0.5 * relative_world.dot(imported_relative));
    energy_error.add(0.5 * next_spin.sub(spin).dot(next_omega));
    energy_error.add(0.5 * spin.dot(next_omega.sub(omega)));
    energy_error.add(loss);
    let balance = AssemblyBalance {
        linear_momentum_residual_world_ns: vector(linear_error.value())?,
        angular_momentum_residual_world_nms: vector(angular_error.value())?,
        energy_residual_j: finite(energy_error.value())?,
    };
    controls.admit(balance)?;
    poll(&mut cancelled)?;
    Ok(PointMassCaptureResult {
        state: next,
        properties: next_properties,
        center_shift_body_m: center_shift,
        captured_impulse_world_ns: vector(orientation.rotate_body_to_world(impulse))?,
        removed_relative_energy_j: loss,
        balance,
    })
}

fn relative_momentum(momentum: Vec3, mass: f64, boost: Vec3) -> Result<Vec3, AssemblyError> {
    vector(Vec3::new(
        (-mass).mul_add(boost.x, momentum.x),
        (-mass).mul_add(boost.y, momentum.y),
        (-mass).mul_add(boost.z, momentum.z),
    ))
}
