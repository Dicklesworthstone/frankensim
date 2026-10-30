//! Original two-body inertia under external drive and directional dry friction.
//!
//! External forces and body torques are explicit, held wrenches. They are not
//! folded into an accepted pre-kick velocity or friction history. The original
//! contact mobility supplies the driven law; opposing integrated traction and
//! both external impulses then update exactly the original two bodies.
//!
//! The relative/contact and complementary body motions separate orthogonally
//! in the mass metric. This gives independent full-body external work from
//! the free body work minus free relative work plus actual native drive work.
//! The complementary motion includes common translation, rotation and normal
//! recoil, so none of their energy is mistaken for frictional heat.
//!
//! Poses, arms and inertia are frozen during this impulse phase. There is no
//! pose transport, geometric acceleration, normal constraint or thermal store.

use super::*;
use crate::Wrench;
use fs_tribo::tangential::elliptic::driven::{
    EllipticDrivenControls, EllipticDrivenError, EllipticDrivenInput, EllipticDrivenMode,
    EllipticDrivenReport,
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidDirectionalDrivenAudit {
    pub body_linear_momentum_residual_ns: [Vec3; 2],
    pub body_angular_momentum_residual_nms: [Vec3; 2],
    /// Total change minus the explicit external linear impulse.
    pub linear_momentum_residual_ns: Vec3,
    /// Total change about the world origin minus explicit external angular impulse.
    pub angular_momentum_residual_nms: Vec3,
    pub endpoint_velocity_residual_m_s: f64,
    /// Duration times W*q-a, checked independently of the native force solve.
    pub drive_velocity_residual_m_s: f64,
    pub drive_work_residual_j: f64,
    /// Actual full kinetic change minus relative plus complementary changes.
    pub relative_energy_residual_j: f64,
    /// Momentum-update work versus actual body kinetic change.
    pub work_energy_residual_j: f64,
    /// Delta K + physical friction heat - external work.
    pub total_energy_residual_j: f64,
}

#[derive(Clone, Debug)]
pub struct RigidDirectionalDrivenReport {
    pub driven: EllipticDrivenReport,
    pub mobility_per_kg: [[f64; 2]; 2],
    pub driving_acceleration_m_s2: [f64; 2],
    pub slip_velocity_before_m_s: [f64; 2],
    pub slip_velocity_after_m_s: [f64; 2],
    pub normal_velocity_before_m_s: f64,
    pub normal_velocity_after_m_s: f64,
    /// Friction impulse only; the external wrenches retain their own work.
    pub force_impulse_on_a_world_ns: Vec3,
    pub kinetic_energy_change_j: f64,
    pub external_work_j: f64,
    /// Algebraic work of all momentum increments against endpoint mean rates.
    /// This is not the physical work of an averaged curved friction force.
    pub momentum_update_work_j: f64,
    pub audit: RigidDirectionalDrivenAudit,
}

#[derive(Debug)]
pub enum RigidDirectionalDrivenError<E = EllipticDrivenError> {
    Mechanics(RigidDirectionalError),
    Driven(E),
    Audit(RigidDirectionalDrivenAudit),
    Invalid(&'static str),
    Cancelled,
}
impl<E> From<DynamicsError> for RigidDirectionalDrivenError<E> {
    fn from(e: DynamicsError) -> Self {
        Self::Mechanics(RigidDirectionalError::Dynamics(e))
    }
}
impl<E> From<RigidDirectionalError> for RigidDirectionalDrivenError<E> {
    fn from(e: RigidDirectionalError) -> Self {
        Self::Mechanics(e)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidDirectionalDrivenError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "driven directional rigid contact: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for RigidDirectionalDrivenError<E> {}

/// Native law with original history and explicit accepted sliding/resting mode.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_directional_driven(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidDirectionalContact,
    external: [Wrench; 2],
    law: &EllipticCoulomb,
    history: &mut EllipticCoulombState,
    mode: &mut EllipticDrivenMode,
    dt_s: f64,
    controls: RigidDirectionalControls,
    driven_controls: EllipticDrivenControls,
    gate: &CancelGate,
) -> Result<RigidDirectionalDrivenReport, RigidDirectionalDrivenError> {
    apply_rigid_directional_driven_with(
        bodies,
        properties,
        contact,
        external,
        history,
        mode,
        dt_s,
        controls,
        |state, input| {
            let report =
                law.advance_driven(state, input, driven_controls, &mut || gate.is_requested())?;
            let next = report.state.clone();
            Ok((report, next))
        },
        gate,
    )
}

/// Apply a source owner's driven trajectory to the original finite body pair.
///
/// `external` follows contact A/B order, not global body-array indexing. Forces
/// are world-frame and torques are in each original body's coordinates. The
/// callback receives the actual relative velocity and full mobility together
/// with their acceleration under these wrenches; it owns source/constitutive
/// admission, physical branch selection and bounded native integration.
///
/// Body momenta, callback history and accepted mode publish jointly after
/// independent response, impulse and work checks. A callback refusal, invalid
/// candidate or cancellation publishes none. Moving contacts cannot acquire
/// static friction merely because their velocity is small: accepted mode is an
/// explicit original owner, and the native event result owns its transition.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn apply_rigid_directional_driven_with<S, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    contact: RigidDirectionalContact,
    external: [Wrench; 2],
    history: &mut S,
    mode: &mut EllipticDrivenMode,
    dt_s: f64,
    controls: RigidDirectionalControls,
    mut resolve: impl FnMut(&S, &EllipticDrivenInput) -> Result<(EllipticDrivenReport, S), E>,
    gate: &CancelGate,
) -> Result<RigidDirectionalDrivenReport, RigidDirectionalDrivenError<E>> {
    poll_driven(gate)?;
    let binding = contact.binding;
    if bodies.len() != properties.len()
        || bodies.len() > controls.maximum_bodies
        || binding.body_a >= bodies.len()
        || binding.body_b >= bodies.len()
        || binding.body_a == binding.body_b
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !contact.normal_load_n.is_finite()
        || contact.normal_load_n < 0.0
        || [
            controls.geometry_tolerance_m,
            controls.velocity_tolerance_m_s,
            controls.linear_momentum_tolerance_ns,
            controls.angular_momentum_tolerance_nms,
            controls.energy_tolerance_j,
            controls.contact.maximum_force_n,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(RigidDirectionalDrivenError::Invalid(
            "bounded original pair, physical duration and finite gates required",
        ));
    }
    for wrench in external {
        wrench.validate()?;
    }
    let before = [bodies[binding.body_a], bodies[binding.body_b]];
    let props = [properties[binding.body_a], properties[binding.body_b]];
    let arms = [binding.arm_a_body, binding.arm_b_body];
    let axes = [
        vector(contact.frame.longitudinal()),
        vector(contact.frame.lateral()),
        vector(contact.frame.normal()),
    ];
    let old = points(before, props, arms)?;
    if magnitude(old[0].point_world.sub(old[1].point_world)) > controls.geometry_tolerance_m {
        return Err(RigidDirectionalDrivenError::Invalid(
            "original contact points do not coincide",
        ));
    }
    let mobility = response::<TriboError>(before, props, arms, axes)?;
    let initial = old[0].point_velocity_world.sub(old[1].point_velocity_world);
    let velocity = [initial.dot(axes[0]), initial.dot(axes[1])];
    let mut free = before;
    let mut point_acceleration = [Vec3::ZERO; 2];
    let mut external_linear = [Vec3::ZERO; 2];
    let mut external_angular = [Vec3::ZERO; 2];
    for i in 0..2 {
        external_linear[i] = external[i].force_world.scale(dt_s);
        external_angular[i] = external[i].torque_body.scale(dt_s);
        free[i] = RigidBodyState::new(
            before[i].pose(),
            before[i].linear_momentum_world().add(external_linear[i]),
            before[i].angular_momentum_body().add(external_angular[i]),
        )?;
        // The original full inverse inertia acts on torque to give frozen-frame
        // angular acceleration. No baseline-velocity subtraction is needed.
        let angular_acceleration =
            props[i].angular_velocity_body_checked(external[i].torque_body)?;
        point_acceleration[i] = external[i].force_world.scale(props[i].mass().recip()).add(
            before[i]
                .pose()
                .orientation()
                .rotate_body_to_world(angular_acceleration)
                .cross(old[i].arm_world),
        );
        if !point_acceleration[i].is_finite()
            || !external_linear[i].is_finite()
            || !external_angular[i].is_finite()
        {
            return Err(RigidDirectionalDrivenError::Invalid(
                "unrepresentable external wrench response",
            ));
        }
    }
    let acceleration = point_acceleration[0].sub(point_acceleration[1]);
    let drive = [
        checked(acceleration.dot(axes[0]))?,
        checked(acceleration.dot(axes[1]))?,
    ];
    let (resolved, next_history) = resolve(
        history,
        &EllipticDrivenInput {
            velocity_m_s: velocity,
            mobility_per_kg: mobility,
            driving_acceleration_m_s2: drive,
            mode: *mode,
            normal_n: contact.normal_load_n,
            dt_s,
            maximum_force_n: controls.contact.maximum_force_n,
            contact_temperature_k: contact.contact_temperature_k,
            nominal_pressure_pa: contact.nominal_pressure_pa,
        },
    )
    .map_err(RigidDirectionalDrivenError::Driven)?;
    poll_driven(gate)?;
    if !resolved.distance_m.is_finite()
        || resolved.distance_m < 0.0
        || !resolved.dissipated_j.is_finite()
        || resolved.dissipated_j < 0.0
        || resolved.mode == EllipticDrivenMode::Resting && resolved.velocity_after_m_s != [0.0; 2]
        || resolved
            .slip_increment_m
            .iter()
            .chain(&resolved.driving_force_n)
            .chain(&resolved.velocity_after_m_s)
            .any(|x| !x.is_finite())
    {
        return Err(RigidDirectionalDrivenError::Invalid(
            "finite native travel, drive and nonnegative heat required",
        ));
    }
    let impulse = axes[0]
        .scale(-resolved.traction_impulse_ns[0])
        .add(axes[1].scale(-resolved.traction_impulse_ns[1]));
    let mut after = free;
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        after[i] = free[i]
            .apply_impulse_at_body_point(props[i], arms[i], impulse.scale(sign))?
            .state_after;
    }
    let free_points = points(free, props, arms)?;
    let new = points(after, props, arms)?;
    let final_rate = new[0].point_velocity_world.sub(new[1].point_velocity_world);
    let endpoint = [final_rate.dot(axes[0]), final_rate.dot(axes[1])];
    let mut kinetic = 0.0;
    let mut free_work = 0.0;
    let mut update_work = 0.0;
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    let mut body_linear = [Vec3::ZERO; 2];
    let mut body_angular = [Vec3::ZERO; 2];
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let dp = after[i]
            .linear_momentum_world()
            .sub(before[i].linear_momentum_world());
        let dl = after[i]
            .angular_momentum_body()
            .sub(before[i].angular_momentum_body());
        let mean_v = old[i]
            .center_of_mass_velocity_world
            .scale(0.5)
            .add(new[i].center_of_mass_velocity_world.scale(0.5));
        let mean_w = old[i]
            .angular_velocity_body
            .scale(0.5)
            .add(new[i].angular_velocity_body.scale(0.5));
        kinetic = checked(kinetic + dp.dot(mean_v) + dl.dot(mean_w))?;
        free_work = checked(
            free_work
                + external_linear[i].dot(
                    old[i]
                        .center_of_mass_velocity_world
                        .scale(0.5)
                        .add(free_points[i].center_of_mass_velocity_world.scale(0.5)),
                )
                + external_angular[i].dot(
                    old[i]
                        .angular_velocity_body
                        .scale(0.5)
                        .add(free_points[i].angular_velocity_body.scale(0.5)),
                ),
        )?;
        let expected_linear = external_linear[i].add(impulse.scale(sign));
        let expected_angular = external_angular[i].add(
            arms[i].cross(
                before[i]
                    .pose()
                    .orientation()
                    .rotate_world_to_body(impulse.scale(sign)),
            ),
        );
        update_work =
            checked(update_work + expected_linear.dot(mean_v) + expected_angular.dot(mean_w))?;
        body_linear[i] = dp.sub(expected_linear);
        body_angular[i] = dl.sub(expected_angular);
        let dp_internal = dp.sub(external_linear[i]);
        linear = linear.add(dp_internal);
        angular = angular
            .add(before[i].pose().position_world().cross(dp_internal))
            .add(
                before[i]
                    .pose()
                    .orientation()
                    .rotate_body_to_world(dl.sub(external_angular[i])),
            );
    }
    let q = resolved.driving_force_n;
    let free_relative_work = checked(
        dt_s * (q[0] * (velocity[0] + 0.5 * dt_s * drive[0])
            + q[1] * (velocity[1] + 0.5 * dt_s * drive[1])),
    )?;
    let complementary_work = checked(free_work - free_relative_work)?;
    let actual_drive_work =
        checked(q[0] * resolved.slip_increment_m[0] + q[1] * resolved.slip_increment_m[1])?;
    let external_work = checked(complementary_work + actual_drive_work)?;
    let drive_error =
        [0, 1].map(|i| dt_s * (mobility[i][0] * q[0] + mobility[i][1] * q[1] - drive[i]));
    let audit = RigidDirectionalDrivenAudit {
        body_linear_momentum_residual_ns: body_linear,
        body_angular_momentum_residual_nms: body_angular,
        linear_momentum_residual_ns: linear,
        angular_momentum_residual_nms: angular,
        endpoint_velocity_residual_m_s: det::hypot(
            endpoint[0] - resolved.velocity_after_m_s[0],
            endpoint[1] - resolved.velocity_after_m_s[1],
        ),
        drive_velocity_residual_m_s: det::hypot(drive_error[0], drive_error[1]),
        drive_work_residual_j: resolved.drive_work_j - actual_drive_work,
        relative_energy_residual_j: kinetic - complementary_work - resolved.kinetic_energy_change_j,
        work_energy_residual_j: kinetic - update_work,
        total_energy_residual_j: kinetic + resolved.dissipated_j - external_work,
    };
    if [
        audit.endpoint_velocity_residual_m_s,
        audit.drive_velocity_residual_m_s,
    ]
    .into_iter()
    .any(|x| !x.is_finite() || x > controls.velocity_tolerance_m_s)
        || [
            audit.drive_work_residual_j,
            audit.relative_energy_residual_j,
            audit.work_energy_residual_j,
            audit.total_energy_residual_j,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x.abs() > controls.energy_tolerance_j)
        || body_linear
            .into_iter()
            .chain([linear])
            .any(|x| !x.is_finite() || magnitude(x) > controls.linear_momentum_tolerance_ns)
        || body_angular
            .into_iter()
            .chain([angular])
            .any(|x| !x.is_finite() || magnitude(x) > controls.angular_momentum_tolerance_nms)
    {
        return Err(RigidDirectionalDrivenError::Audit(audit));
    }
    let next_mode = resolved.mode;
    let report = RigidDirectionalDrivenReport {
        driven: resolved,
        mobility_per_kg: mobility,
        driving_acceleration_m_s2: drive,
        slip_velocity_before_m_s: velocity,
        slip_velocity_after_m_s: endpoint,
        normal_velocity_before_m_s: initial.dot(axes[2]),
        normal_velocity_after_m_s: final_rate.dot(axes[2]),
        force_impulse_on_a_world_ns: impulse,
        kinetic_energy_change_j: kinetic,
        external_work_j: external_work,
        momentum_update_work_j: update_work,
        audit,
    };
    poll_driven(gate)?;
    bodies[binding.body_a] = after[0];
    bodies[binding.body_b] = after[1];
    *history = next_history;
    *mode = next_mode;
    Ok(report)
}

fn checked<E>(x: f64) -> Result<f64, RigidDirectionalDrivenError<E>> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(RigidDirectionalDrivenError::Invalid(
            "unrepresentable driven rigid work or response",
        ))
    }
}
fn poll_driven<E>(gate: &CancelGate) -> Result<(), RigidDirectionalDrivenError<E>> {
    if gate.is_requested() {
        Err(RigidDirectionalDrivenError::Cancelled)
    } else {
        Ok(())
    }
}
