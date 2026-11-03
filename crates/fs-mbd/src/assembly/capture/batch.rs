//! Simultaneous rigid sticking from one original geometry and velocity snapshot.
//!
//! All incoming parcels join one final rigid velocity field. Their interactions
//! share the body's translation and full central inertia, so the batch is not
//! evaluated as independent single-parcel impulses or independent losses.

use super::{
    AssemblyBalance, AssemblyControls, AssemblyError, ENTRIES, MassProperties, PointMassCapture,
    PointMassCaptureResult, Pose, RigidBodyState, Sum, Vec3, VectorSum, array, finite, poll,
    relative_momentum, vector,
};

/// Capture all supplied point parcels in one instantaneous physical event.
///
/// Every arm is measured from the SAME original COM in its body axes. This
/// preserves each original material location without sequential rebasing.
/// The result reuses the native single-capture shape: its reported impulse is
/// the aggregate impulse on the original body, and its sticking loss covers
/// the whole event. Original inventory, event choice and heat ownership remain
/// with the caller.
///
/// The reductions require O(N) work and constant auxiliary storage. One
/// aggregate central inertia is admitted through the existing native owner;
/// there is no parcel-by-parcel factorization or N-by-N mobility matrix.
/// A single entry delegates to the original single-parcel path.
///
/// The part budget includes the original body and every supplied descriptor,
/// including zero-mass entries, bounding all scans before they begin. Each
/// positive parcel consumes one impulse slot. Empty or all-zero input is an
/// exact identity after validation. A positive aggregate mass addition that
/// rounds away, nonfinite quantities, failed balances, and cancellation return
/// no endpoint and leave all borrowed inputs unchanged.
pub fn capture_point_masses(
    body: &RigidBodyState,
    properties: MassProperties,
    incoming: &[PointMassCapture],
    controls: AssemblyControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PointMassCaptureResult, AssemblyError> {
    poll(&mut cancelled)?;
    let parts = incoming
        .len()
        .checked_add(1)
        .ok_or(AssemblyError::Invalid("capture part count overflow"))?;
    controls.validate(parts, 0)?;
    if let [single] = incoming {
        return super::capture_point_mass(body, properties, *single, controls, cancelled);
    }
    properties.validate()?;
    let mass = properties.mass();
    let orientation = body.pose().orientation();
    let velocity = body.center_of_mass_velocity_world(properties)?;
    let spin = body.angular_momentum_body();
    let omega = properties.angular_velocity_body_checked(spin)?;
    let mut added = Sum::default();
    let mut active = 0;
    let mut all_at_center = true;
    for parcel in incoming {
        poll(&mut cancelled)?;
        if !parcel.mass_kg.is_finite() || parcel.mass_kg < 0.0 {
            return Err(AssemblyError::Invalid(
                "finite nonnegative captured mass required",
            ));
        }
        vector(parcel.attachment_arm_body_m)?;
        vector(parcel.velocity_world_m_s)?;
        if parcel.mass_kg > 0.0 {
            added.add(parcel.mass_kg);
            active += 1;
            all_at_center &= parcel.attachment_arm_body_m == Vec3::ZERO;
        }
    }
    controls.validate(parts, active)?;
    if active == 0 {
        poll(&mut cancelled)?;
        return Ok(PointMassCaptureResult {
            state: *body,
            properties,
            center_shift_body_m: Vec3::ZERO,
            captured_impulse_world_ns: Vec3::ZERO,
            removed_relative_energy_j: 0.0,
            balance: AssemblyBalance {
                linear_momentum_residual_world_ns: Vec3::ZERO,
                angular_momentum_residual_world_nms: Vec3::ZERO,
                energy_residual_j: 0.0,
            },
        });
    }
    let added = finite(added.value())?;
    let total = finite(mass + added)?;
    if total <= mass || total <= added {
        return Err(AssemblyError::Invalid(
            "body and incoming mass must remain represented",
        ));
    }

    // Normalized first moments avoid forming a large m*r before dividing by T.
    let mut center = VectorSum::default();
    for parcel in incoming {
        poll(&mut cancelled)?;
        if parcel.mass_kg > 0.0 {
            center.add(vector(
                parcel.attachment_arm_body_m.scale(parcel.mass_kg / total),
            )?);
        }
    }
    let center_shift = vector(center.value())?;
    let before_relative = relative_momentum(body.linear_momentum_world(), mass, velocity)?;
    let before_body_relative = vector(orientation.rotate_world_to_body(before_relative))?;
    let mut inertia = properties.inertia_tensor_body();
    let mut entries = [Sum::default(); 6];
    for (entry, (row, col)) in entries.iter_mut().zip(ENTRIES) {
        entry.add(inertia.m[3 * row + col]);
    }
    add_parallel_axis(&mut entries, center_shift, mass)?;
    let mut angular = VectorSum::default();
    angular.add(spin);
    angular.add(vector(
        center_shift.cross(before_body_relative).scale(-1.0),
    )?);
    let mut linear = VectorSum::default();
    linear.add(body.linear_momentum_world());
    for parcel in incoming {
        poll(&mut cancelled)?;
        if parcel.mass_kg == 0.0 {
            continue;
        }
        let arm = vector(parcel.attachment_arm_body_m.sub(center_shift))?;
        add_parallel_axis(&mut entries, arm, parcel.mass_kg)?;
        let relative = vector(
            orientation.rotate_world_to_body(vector(parcel.velocity_world_m_s.sub(velocity))?),
        )?;
        angular.add(vector(arm.cross(vector(relative.scale(parcel.mass_kg))?))?);
        linear.add(vector(parcel.velocity_world_m_s.scale(parcel.mass_kg))?);
    }
    for (entry, (row, col)) in entries.into_iter().zip(ENTRIES) {
        let value = finite(entry.value())?;
        inertia.m[3 * row + col] = value;
        inertia.m[3 * col + row] = value;
    }
    let next_properties = if all_at_center {
        MassProperties::from_principal_frame(
            total,
            Vec3::ZERO,
            properties.principal_inertia_body(),
            properties.principal_axes_to_body(),
        )?
    } else {
        MassProperties::from_inertia_tensor(total, Vec3::ZERO, inertia)?
    };
    let next_spin = vector(angular.value())?;
    let next = RigidBodyState::new(
        Pose::new(
            body.pose().point_world_from_body(center_shift)?,
            orientation,
        )?,
        vector(linear.value())?,
        next_spin,
    )?;
    next.center_of_mass_velocity_world(next_properties)?;
    let next_omega = next_properties.angular_velocity_body_checked(next_spin)?;
    let after_relative = relative_momentum(next.linear_momentum_world(), total, velocity)?;
    let after_body_velocity =
        vector(orientation.rotate_world_to_body(vector(after_relative.scale(total.recip()))?))?;
    let original_velocity_change = vector(
        after_body_velocity
            .sub(next_omega.cross(center_shift))
            .sub(before_body_relative.scale(mass.recip())),
    )?;

    // Positive final-rigid-motion projection norm. Shared translation and spin
    // are solved once; independent per-hit losses would miss their coupling.
    let mut removed = Sum::default();
    removed.add(scaled_energy(original_velocity_change, mass)?);
    let delta_omega = vector(properties.to_principal(next_omega.sub(omega)))?;
    let moments = properties.principal_inertia_body();
    let spin_residual = vector(Vec3::new(
        delta_omega.x * moments.x.sqrt(),
        delta_omega.y * moments.y.sqrt(),
        delta_omega.z * moments.z.sqrt(),
    ))?;
    removed.add(finite(0.5 * spin_residual.dot(spin_residual))?);

    // Use actual returned momentum and COM displacement in the balance audit.
    let actual_shift = vector(
        next.pose()
            .position_world()
            .sub(body.pose().position_world()),
    )?;
    let mut linear_error = VectorSum::default();
    linear_error.add(next.linear_momentum_world());
    linear_error.add(body.linear_momentum_world().scale(-1.0));
    let mut angular_error = VectorSum::default();
    angular_error.add(orientation.rotate_body_to_world(next_spin.sub(spin)));
    angular_error.add(actual_shift.cross(after_relative));
    let mut first_moment_error = VectorSum::default();
    first_moment_error.add(vector(actual_shift.scale(total))?);
    let mut energy_error = Sum::default();
    energy_error.add(0.5 * after_relative.dot(after_relative.scale(total.recip())));
    energy_error.add(-0.5 * before_relative.dot(before_relative.scale(mass.recip())));
    energy_error.add(0.5 * next_spin.sub(spin).dot(next_omega));
    energy_error.add(0.5 * spin.dot(next_omega.sub(omega)));
    for parcel in incoming {
        poll(&mut cancelled)?;
        if parcel.mass_kg == 0.0 {
            continue;
        }
        let world_arm = vector(orientation.rotate_body_to_world(parcel.attachment_arm_body_m))?;
        let relative_world = vector(parcel.velocity_world_m_s.sub(velocity))?;
        let imported_relative = vector(relative_world.scale(parcel.mass_kg))?;
        let relative_body = vector(orientation.rotate_world_to_body(relative_world))?;
        let arm = vector(parcel.attachment_arm_body_m.sub(center_shift))?;
        let final_velocity = vector(after_body_velocity.add(next_omega.cross(arm)))?;
        removed.add(scaled_energy(
            vector(relative_body.sub(final_velocity))?,
            parcel.mass_kg,
        )?);
        linear_error.add(vector(parcel.velocity_world_m_s.scale(-parcel.mass_kg))?);
        angular_error.add(vector(world_arm.cross(imported_relative).scale(-1.0))?);
        first_moment_error.add(vector(world_arm.scale(-parcel.mass_kg))?);
        energy_error.add(-scaled_energy(relative_world, parcel.mass_kg)?);
    }
    let removed = finite(removed.value())?;
    if removed < 0.0 {
        return Err(AssemblyError::Invalid("negative captured relative energy"));
    }
    angular_error.add(vector(first_moment_error.value())?.cross(velocity));
    energy_error.add(removed);
    let balance = AssemblyBalance {
        linear_momentum_residual_world_ns: vector(linear_error.value())?,
        angular_momentum_residual_world_nms: vector(angular_error.value())?,
        energy_residual_j: finite(energy_error.value())?,
    };
    controls.admit(balance)?;
    let impulse =
        vector(orientation.rotate_body_to_world(vector(original_velocity_change.scale(mass))?))?;
    poll(&mut cancelled)?;
    Ok(PointMassCaptureResult {
        state: next,
        properties: next_properties,
        center_shift_body_m: center_shift,
        captured_impulse_world_ns: impulse,
        removed_relative_energy_j: removed,
        balance,
    })
}

fn scaled_energy(velocity: Vec3, mass: f64) -> Result<f64, AssemblyError> {
    let scaled = vector(velocity.scale(mass.sqrt()))?;
    finite(0.5 * scaled.dot(scaled))
}

// Positive shifts about the final COM avoid subtracting nearly equal global
// inertia tensors when incoming material is far from the original COM.
fn add_parallel_axis(entries: &mut [Sum; 6], arm: Vec3, mass: f64) -> Result<(), AssemblyError> {
    let r = array(vector(arm.scale(mass.sqrt()))?);
    for (entry, (row, col)) in entries.iter_mut().zip(ENTRIES) {
        entry.add(finite(if row == col {
            r[(row + 1) % 3].powi(2) + r[(row + 2) % 3].powi(2)
        } else {
            -r[row] * r[col]
        })?);
    }
    Ok(())
}
