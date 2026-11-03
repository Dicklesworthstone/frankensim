//! Simultaneous point-parcel detachment from one original rigid snapshot.

use super::{PointMassRelease, ReleasedPointMass, relative_momentum};
use crate::assembly::{
    AssemblyBalance, AssemblyControls, AssemblyError, ENTRIES, Sum, VectorSum, array, finite, poll,
    vector,
};
use crate::{MassProperties, Pose, RigidBodyState, Vec3};

/// One simultaneous release, retaining the original input order of positive masses.
#[derive(Clone, Debug, PartialEq)]
pub struct PointMassReleaseBatchResult {
    /// Retained body's native endpoint at its new COM, in the original body axes.
    pub state: RigidBodyState,
    /// Retained native mass and full central inertia.
    pub properties: MassProperties,
    /// New COM relative to the original COM [m]. Subtract this from every
    /// original body point before continuing with the returned native state.
    pub center_shift_body_m: Vec3,
    /// Canonical free parcels in positive-mass input order. Zero-mass entries
    /// require zero impulse and produce no parcel.
    pub parcels: Vec<ReleasedPointMass>,
    /// Work supplied by ALL separation impulses together [J]. Shared retained-
    /// body translational and rotational cross terms are included exactly once.
    pub supplied_work_j: f64,
    /// Actual returned linear, spin-plus-orbital and energy/work residuals.
    pub balance: AssemblyBalance,
}

/// Release several existing point parcels at one physical instant.
///
/// Every input arm and inherited velocity refers to the SAME original body,
/// before any mass removal or impulse. With retained mass M and COM shift a,
/// the central inertia is I - sum(m_i Q(r_i)) - M Q(a), where
/// Q(r)=|r|²*1-r*r^T. Each J_i acts on its parcel and -J_i on the retained body
/// at the original physical point. The retained angular impulse is
/// -sum((r_i-a) cross J_i).
///
/// Separation work is sum(|J_i|²/(2m_i)) + |sum(J_i)|²/(2M)
/// + torque^T I_retained^-1 torque/2. Summing independent single-release work
/// would omit the shared-body cross terms. Initial-spin terms cancel because
/// both sides of every impulse share the same original point velocity.
///
/// The caller identifies material already included in the original inertia
/// and owns the event decision and energy source. This is instantaneous point
/// detachment, without finite parcel extent, fracture or continuous mass flow.
///
/// Fixed-order compensated passes take O(inputs) work and one native retained-
/// inertia admission. Only returned parcels need O(inputs) storage. Each input
/// slot consumes one part-budget slot in addition to the retained body, even
/// when its mass is zero; this bounds scanning before allocation. Each nonzero
/// impulse consumes one impulse slot. Empty or all-zero input is an identity.
/// All inputs are borrowed; cancellation or refusal publishes no partial batch.
pub fn release_point_masses(
    body: &RigidBodyState,
    properties: MassProperties,
    outgoing: &[PointMassRelease],
    controls: AssemblyControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PointMassReleaseBatchResult, AssemblyError> {
    poll(&mut cancelled)?;
    let parts = outgoing
        .len()
        .checked_add(1)
        .ok_or(AssemblyError::Invalid("point-release input count overflow"))?;
    controls.validate(parts, 0)?;
    properties.validate()?;
    let velocity = body.center_of_mass_velocity_world(properties)?;
    let spin = body.angular_momentum_body();
    let omega = properties.angular_velocity_body_checked(spin)?;
    let mut removed = Sum::default();
    let mut active = 0;
    let mut impulses = 0;
    let mut changed_inertia = false;
    for input in outgoing {
        poll(&mut cancelled)?;
        if !input.mass_kg.is_finite() || input.mass_kg < 0.0 {
            return Err(AssemblyError::Invalid(
                "finite nonnegative released parcel masses required",
            ));
        }
        vector(input.attachment_arm_body_m)?;
        vector(input.impulse_on_parcel_body_ns)?;
        let kicked = input.impulse_on_parcel_body_ns != Vec3::ZERO;
        if input.mass_kg == 0.0 {
            if kicked {
                return Err(AssemblyError::Invalid(
                    "zero mass cannot receive a separation impulse",
                ));
            }
            continue;
        }
        active += 1;
        impulses += usize::from(kicked);
        changed_inertia |= input.attachment_arm_body_m != Vec3::ZERO;
        removed.add(input.mass_kg);
    }
    controls.validate(parts, impulses)?;
    if active == 0 {
        poll(&mut cancelled)?;
        return Ok(PointMassReleaseBatchResult {
            state: *body,
            properties,
            center_shift_body_m: Vec3::ZERO,
            parcels: Vec::new(),
            supplied_work_j: 0.0,
            balance: AssemblyBalance {
                linear_momentum_residual_world_ns: Vec3::ZERO,
                angular_momentum_residual_world_nms: Vec3::ZERO,
                energy_residual_j: 0.0,
            },
        });
    }
    let total = properties.mass();
    let removed = finite(removed.value())?;
    let retained = finite(total - removed)?;
    if retained <= 0.0 || retained >= total {
        return Err(AssemblyError::Invalid(
            "positive represented retained and released masses required",
        ));
    }

    // Weight by the retained mass before multiplying positions: a raw
    // sum(m_i*r_i) can overflow even when the COM shift is representable.
    let mut center_shift = VectorSum::default();
    let original_inertia = properties.inertia_tensor_body();
    let mut inertia = [Sum::default(); 6];
    for (entry, &(row, col)) in inertia.iter_mut().zip(&ENTRIES) {
        entry.add(original_inertia.m[3 * row + col]);
    }
    for input in outgoing {
        poll(&mut cancelled)?;
        if input.mass_kg == 0.0 {
            continue;
        }
        center_shift.add(vector(
            input.attachment_arm_body_m.scale(-input.mass_kg / retained),
        )?);
        subtract_parallel_axis(
            &mut inertia,
            vector(input.attachment_arm_body_m.scale(input.mass_kg.sqrt()))?,
        )?;
    }
    let center_shift = vector(center_shift.value())?;
    let next_properties = if changed_inertia {
        subtract_parallel_axis(&mut inertia, vector(center_shift.scale(retained.sqrt()))?)?;
        let mut tensor = original_inertia;
        for (entry, &(row, col)) in inertia.iter().zip(&ENTRIES) {
            let value = finite(entry.value())?;
            tensor.m[3 * row + col] = value;
            tensor.m[3 * col + row] = value;
        }
        MassProperties::from_inertia_tensor(retained, Vec3::ZERO, tensor)?
    } else {
        // Central point mass has no rotational inertia. Preserve the original
        // principal representation and spin instead of reconstructing either.
        MassProperties::from_principal_frame(
            retained,
            Vec3::ZERO,
            properties.principal_inertia_body(),
            properties.principal_axes_to_body(),
        )?
    };

    let orientation = body.pose().orientation();
    let mut parcels = Vec::with_capacity(active);
    let mut retained_linear = VectorSum::default();
    retained_linear.add(body.linear_momentum_world());
    let mut total_impulse = VectorSum::default();
    let mut angular_impulse = VectorSum::default();
    let mut work = Sum::default();
    for input in outgoing {
        poll(&mut cancelled)?;
        if input.mass_kg == 0.0 {
            continue;
        }
        let arm = input.attachment_arm_body_m;
        let impulse = input.impulse_on_parcel_body_ns;
        let relative = vector(omega.cross(arm).scale(input.mass_kg).add(impulse))?;
        let relative_world = vector(orientation.rotate_body_to_world(relative))?;
        let parcel = ReleasedPointMass {
            mass_kg: input.mass_kg,
            position_world_m: body.pose().point_world_from_body(arm)?,
            linear_momentum_world_ns: boosted_momentum(input.mass_kg, velocity, relative_world)?,
        };
        parcel.velocity_world_m_s()?;
        retained_linear.add(parcel.linear_momentum_world_ns.scale(-1.0));
        total_impulse.add(impulse);
        angular_impulse.add(vector(arm.sub(center_shift).cross(impulse).scale(-1.0))?);
        work.add(half_squared_norm(vector(
            impulse.scale(input.mass_kg.sqrt().recip()),
        )?)?);
        parcels.push(parcel);
    }
    let total_impulse = vector(total_impulse.value())?;
    let angular_impulse = vector(angular_impulse.value())?;
    let inherited_spin = if changed_inertia {
        next_properties.angular_momentum_body_checked(omega)?
    } else {
        spin
    };
    let next_spin = vector(inherited_spin.add(angular_impulse))?;
    let next = RigidBodyState::new(
        Pose::new(
            body.pose().point_world_from_body(center_shift)?,
            orientation,
        )?,
        vector(retained_linear.value())?,
        next_spin,
    )?;
    next.center_of_mass_velocity_world(next_properties)?;
    let next_omega = next_properties.angular_velocity_body_checked(next_spin)?;
    work.add(half_squared_norm(vector(
        total_impulse.scale(retained.sqrt().recip()),
    )?)?);
    work.add(half_squared_norm(vector(
        next_properties.angular_mobility_factor(angular_impulse),
    )?)?);
    let work = finite(work.value())?;

    // Audit canonical returned values, not the intended paired impulses.
    // Common-frame momenta retain small recoil under a large uniform boost.
    let before_relative = relative_momentum(body.linear_momentum_world(), total, velocity)?;
    let after_relative = relative_momentum(next.linear_momentum_world(), retained, velocity)?;
    let actual_shift = vector(
        next.pose()
            .position_world()
            .sub(body.pose().position_world()),
    )?;
    let mut first_moment = VectorSum::default();
    first_moment.add(vector(actual_shift.scale(retained))?);
    let mut linear_error = VectorSum::default();
    linear_error.add(next.linear_momentum_world());
    linear_error.add(body.linear_momentum_world().scale(-1.0));
    let mut angular_error = VectorSum::default();
    angular_error.add(orientation.rotate_body_to_world(next_spin.sub(spin)));
    angular_error.add(vector(actual_shift.cross(after_relative))?);
    let mut energy_error = Sum::default();
    energy_error.add(half_squared_norm(vector(
        after_relative.scale(retained.sqrt().recip()),
    )?)?);
    energy_error.add(-half_squared_norm(vector(
        before_relative.scale(total.sqrt().recip()),
    )?)?);
    for parcel in &parcels {
        poll(&mut cancelled)?;
        let relative =
            relative_momentum(parcel.linear_momentum_world_ns, parcel.mass_kg, velocity)?;
        let actual_arm = vector(parcel.position_world_m.sub(body.pose().position_world()))?;
        first_moment.add(vector(actual_arm.scale(parcel.mass_kg))?);
        linear_error.add(parcel.linear_momentum_world_ns);
        angular_error.add(vector(actual_arm.cross(relative))?);
        energy_error.add(half_squared_norm(vector(
            relative.scale(parcel.mass_kg.sqrt().recip()),
        )?)?);
    }
    angular_error.add(vector(first_moment.value())?.cross(velocity));
    energy_error.add(0.5 * next_spin.sub(spin).dot(next_omega));
    energy_error.add(0.5 * spin.dot(next_omega.sub(omega)));
    energy_error.add(-work);
    let balance = AssemblyBalance {
        linear_momentum_residual_world_ns: vector(linear_error.value())?,
        angular_momentum_residual_world_nms: vector(angular_error.value())?,
        energy_residual_j: finite(energy_error.value())?,
    };
    controls.admit(balance)?;
    poll(&mut cancelled)?;
    Ok(PointMassReleaseBatchResult {
        state: next,
        properties: next_properties,
        center_shift_body_m: center_shift,
        parcels,
        supplied_work_j: work,
        balance,
    })
}

// Six independent symmetric entries; compensate the complete subtraction
// before sending one retained tensor through the original native admission.
fn subtract_parallel_axis(entries: &mut [Sum; 6], scaled_arm: Vec3) -> Result<(), AssemblyError> {
    let r = array(scaled_arm);
    for (entry, &(row, col)) in entries.iter_mut().zip(&ENTRIES) {
        let value = if row == col {
            r[(row + 1) % 3].powi(2) + r[(row + 2) % 3].powi(2)
        } else {
            -r[row] * r[col]
        };
        entry.add(-finite(value)?);
    }
    Ok(())
}

fn half_squared_norm(value: Vec3) -> Result<f64, AssemblyError> {
    finite(value.scale(0.5).dot(value))
}

fn boosted_momentum(mass: f64, boost: Vec3, relative: Vec3) -> Result<Vec3, AssemblyError> {
    vector(Vec3::new(
        mass.mul_add(boost.x, relative.x),
        mass.mul_add(boost.y, relative.y),
        mass.mul_add(boost.z, relative.z),
    ))
}
