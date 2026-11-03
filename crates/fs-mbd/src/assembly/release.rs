//! Release a point parcel with full rigid recoil and changing central inertia.
//!
//! The parcel must already belong to the supplied mass distribution. It has
//! zero intrinsic inertia and initially shares the body's rigid velocity at
//! the selected material point. The caller chooses the release event and owns
//! the source of any separation work; no fracture or continuous-flow law is
//! inferred.

use super::{
    AssemblyBalance, AssemblyControls, AssemblyError, ENTRIES, Sum, VectorSum, array, finite, poll,
    vector,
};
use crate::{MassProperties, Pose, RigidBodyState, Vec3};

/// Material removed from the original body's existing mass distribution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointMassRelease {
    /// Nonnegative parcel mass [kg], strictly below the original body mass.
    pub mass_kg: f64,
    /// Parcel location relative to the original COM, in original body axes [m].
    pub attachment_arm_body_m: Vec3,
    /// Internal force impulse received by the parcel, in original body axes
    /// [N s]. The retained body receives its negative at this same point.
    /// A point parcel has no intrinsic spin and cannot receive a free couple.
    pub impulse_on_parcel_body_ns: Vec3,
}

/// A free zero-intrinsic-inertia parcel, retaining its canonical momentum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReleasedPointMass {
    /// Positive mass [kg].
    pub mass_kg: f64,
    /// Actual returned position in the original world frame [m].
    pub position_world_m: Vec3,
    /// Actual returned linear momentum in the original world frame [N s].
    pub linear_momentum_world_ns: Vec3,
}

impl ReleasedPointMass {
    /// Recover velocity from the stored momentum, refusing invalid or
    /// unrepresentable public values rather than inventing a rigid inertia.
    pub fn velocity_world_m_s(self) -> Result<Vec3, AssemblyError> {
        if !self.mass_kg.is_finite() || self.mass_kg <= 0.0 {
            return Err(AssemblyError::Invalid(
                "positive released parcel mass required",
            ));
        }
        vector(Vec3::new(
            self.linear_momentum_world_ns.x / self.mass_kg,
            self.linear_momentum_world_ns.y / self.mass_kg,
            self.linear_momentum_world_ns.z / self.mass_kg,
        ))
    }
}

/// One atomic release, with no mutation of the original native owners.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointMassReleaseResult {
    /// Retained body's native endpoint at its new COM.
    pub state: RigidBodyState,
    /// Retained full central inertia and mass in the original body axes.
    pub properties: MassProperties,
    /// New COM relative to the old COM [m]. Subtract this from original body
    /// points before continuing with the returned native state.
    pub center_shift_body_m: Vec3,
    /// The released parcel; absent only for an exact zero-mass identity.
    pub parcel: Option<ReleasedPointMass>,
    /// Nonnegative mechanical work supplied by the separation impulse [J].
    /// Its energy source must be debited by the caller.
    pub supplied_work_j: f64,
    /// Actual returned momentum and common-frame energy/work residuals.
    pub balance: AssemblyBalance,
}

/// Detach existing point mass and optionally impart equal/opposite recoil.
///
/// For original mass T, parcel mass m, remaining mass M=T-m and original arm
/// r, the retained COM shift is a=-m*r/M and its central inertia is
/// I - (m*T/M)*(|r|^2*1-r*r^T). The remaining inertia must pass the original
/// native positive-definiteness and rigid-inertia triangle admission. This
/// cannot establish that the caller's original geometry contained the parcel.
///
/// Both outputs inherit the original rigid velocity field before a prescribed
/// internal impulse J acts on the parcel and -J on the retained body at the
/// same point. The retained angular impulse is -(r-a) cross J. Its nonnegative
/// separation work uses translational and rotational impulse quadratics; no
/// subtraction of large laboratory kinetic energies determines that work.
///
/// Positive release consumes two part slots. A nonzero prescribed impulse
/// consumes one impulse slot; zero-impulse detachment consumes none. Zero mass
/// with zero impulse is an exact one-part identity. Refusal or cancellation
/// leaves the borrowed body and mass properties unchanged.
pub fn release_point_mass(
    body: &RigidBodyState,
    properties: MassProperties,
    outgoing: PointMassRelease,
    controls: AssemblyControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PointMassReleaseResult, AssemblyError> {
    poll(&mut cancelled)?;
    properties.validate()?;
    if !outgoing.mass_kg.is_finite() || outgoing.mass_kg < 0.0 {
        return Err(AssemblyError::Invalid(
            "finite nonnegative released mass required",
        ));
    }
    let arm = vector(outgoing.attachment_arm_body_m)?;
    let impulse = vector(outgoing.impulse_on_parcel_body_ns)?;
    let active = usize::from(outgoing.mass_kg > 0.0);
    let kicked = usize::from(impulse != Vec3::ZERO);
    controls.validate(1 + active, kicked)?;
    if active == 0 && kicked != 0 {
        return Err(AssemblyError::Invalid(
            "zero mass cannot receive a separation impulse",
        ));
    }
    let velocity = body.center_of_mass_velocity_world(properties)?;
    let spin = body.angular_momentum_body();
    let omega = properties.angular_velocity_body_checked(spin)?;
    if active == 0 {
        poll(&mut cancelled)?;
        return Ok(PointMassReleaseResult {
            state: *body,
            properties,
            center_shift_body_m: Vec3::ZERO,
            parcel: None,
            supplied_work_j: 0.0,
            balance: AssemblyBalance {
                linear_momentum_residual_world_ns: Vec3::ZERO,
                angular_momentum_residual_world_nms: Vec3::ZERO,
                energy_residual_j: 0.0,
            },
        });
    }
    let total = properties.mass();
    let removed = outgoing.mass_kg;
    let retained = finite(total - removed)?;
    if retained <= 0.0 || retained >= total {
        return Err(AssemblyError::Invalid(
            "positive represented retained and released masses required",
        ));
    }
    let center_shift = vector(arm.scale(-removed / retained))?;
    let separation = vector(arm.sub(center_shift))?;
    let next_properties = if arm == Vec3::ZERO {
        MassProperties::from_principal_frame(
            retained,
            Vec3::ZERO,
            properties.principal_inertia_body(),
            properties.principal_axes_to_body(),
        )?
    } else {
        // Keep sqrt(m*T/M)*r representable without forming m*T first.
        let scaled_arm = vector(arm.scale(removed.sqrt() * (total / retained).sqrt()))?;
        let r = array(scaled_arm);
        let mut inertia = properties.inertia_tensor_body();
        for (row, col) in ENTRIES {
            let shift = if row == col {
                r[(row + 1) % 3].powi(2) + r[(row + 2) % 3].powi(2)
            } else {
                -r[row] * r[col]
            };
            let value = finite(inertia.m[3 * row + col] - shift)?;
            inertia.m[3 * row + col] = value;
            inertia.m[3 * col + row] = value;
        }
        MassProperties::from_inertia_tensor(retained, Vec3::ZERO, inertia)?
    };
    poll(&mut cancelled)?;
    let angular_impulse = vector(separation.cross(impulse).scale(-1.0))?;
    let inherited_spin = if arm == Vec3::ZERO {
        spin
    } else {
        next_properties.angular_momentum_body_checked(omega)?
    };
    let next_spin = vector(inherited_spin.add(angular_impulse))?;
    let orientation = body.pose().orientation();
    // Use the same paired relative momentum for both outputs. This retains
    // orbital momentum from the original spin as well as the added impulse.
    let relative = vector(omega.cross(arm).scale(removed).add(impulse))?;
    let relative_world = vector(orientation.rotate_body_to_world(relative))?;
    let next = RigidBodyState::new(
        Pose::new(
            body.pose().point_world_from_body(center_shift)?,
            orientation,
        )?,
        vector(velocity.scale(retained).sub(relative_world))?,
        next_spin,
    )?;
    let parcel = ReleasedPointMass {
        mass_kg: removed,
        position_world_m: body.pose().point_world_from_body(arm)?,
        linear_momentum_world_ns: vector(velocity.scale(removed).add(relative_world))?,
    };
    let next_omega = next_properties.angular_velocity_body_checked(next_spin)?;
    next.center_of_mass_velocity_world(next_properties)?;
    parcel.velocity_world_m_s()?;
    let parcel_response = vector(impulse.scale(removed.sqrt().recip()))?;
    let retained_response = vector(impulse.scale(retained.sqrt().recip()))?;
    let spin_response = vector(next_properties.angular_mobility_factor(angular_impulse))?;
    let work = finite(
        0.5 * (parcel_response.dot(parcel_response)
            + retained_response.dot(retained_response)
            + spin_response.dot(spin_response)),
    )?;

    // Audit actual stored momenta and both returned positions. FMA removes the
    // common boost before products, retaining small recoil and spin balances.
    let before_relative = relative_momentum(body.linear_momentum_world(), total, velocity)?;
    let after_relative = relative_momentum(next.linear_momentum_world(), retained, velocity)?;
    let parcel_relative = relative_momentum(parcel.linear_momentum_world_ns, removed, velocity)?;
    let actual_shift = vector(
        next.pose()
            .position_world()
            .sub(body.pose().position_world()),
    )?;
    let actual_arm = vector(parcel.position_world_m.sub(body.pose().position_world()))?;
    let mut first_moment = VectorSum::default();
    first_moment.add(actual_shift.scale(retained));
    first_moment.add(actual_arm.scale(removed));
    let mut linear_error = VectorSum::default();
    linear_error.add(next.linear_momentum_world());
    linear_error.add(parcel.linear_momentum_world_ns);
    linear_error.add(body.linear_momentum_world().scale(-1.0));
    let mut angular_error = VectorSum::default();
    angular_error.add(orientation.rotate_body_to_world(next_spin.sub(spin)));
    angular_error.add(actual_shift.cross(after_relative));
    angular_error.add(actual_arm.cross(parcel_relative));
    angular_error.add(vector(first_moment.value())?.cross(velocity));
    let mut energy_error = Sum::default();
    energy_error.add(0.5 * after_relative.dot(after_relative.scale(retained.recip())));
    energy_error.add(0.5 * parcel_relative.dot(parcel_relative.scale(removed.recip())));
    energy_error.add(-0.5 * before_relative.dot(before_relative.scale(total.recip())));
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
    Ok(PointMassReleaseResult {
        state: next,
        properties: next_properties,
        center_shift_body_m: center_shift,
        parcel: Some(parcel),
        supplied_work_j: work,
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
