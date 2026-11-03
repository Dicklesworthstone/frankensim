//! Isolated normal impacts between a native rod attachment and a finite body.
//!
//! The collision owner supplies an oriented pair of material witnesses. This
//! frozen-pose event uses the existing scalar contact law and both participants'
//! original impulse owners. It does not select simultaneous contacts, advance
//! time, repair penetration or turn restitution loss into thermal heat.

/// Full normal/tangential Coulomb impacts on the same continuing participants.
pub mod friction;

use super::impulse::{RodImpulseControls, RodImpulseReport, RodPointImpulse};
use super::{
    RodMotion, RodMotionError, RodMotionInertia, add, checkpoint, cross, finite, rotate, sub,
};
use fs_contact::impact::{
    NormalImpactError, NormalImpactInput, NormalImpactResult, solve_normal_impact,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_mbd::{DynamicsError, ImpulseReceipt, MassProperties, RigidBodyState, Vec3};

/// Original attachment points and the separating normal directed onto the rod.
#[derive(Clone, Copy, Debug)]
pub struct RodRigidNormalContact {
    pub node: usize,
    /// Offset from the rod's MATERIAL NODE in its current material axes [m].
    pub point_body_m: [f64; 3],
    /// Offset from the rigid body's CENTER OF MASS in its body axes [m].
    pub rigid_arm_body_m: [f64; 3],
    /// Nonzero world direction; normalized by the original rigid inertia owner.
    pub normal_world: [f64; 3],
    pub activation_tolerance_m: f64,
    pub restitution: f64,
}

/// Independent absolute bounds on the represented native response.
#[derive(Clone, Copy, Debug)]
pub struct RodRigidImpactControls {
    /// Original impulse bounds, also applied to the combined momentum/work.
    pub impulse: RodImpulseControls,
    /// Maximum tangential separation between the two original witnesses [m].
    /// Its actual moment still enters the angular momentum acceptance check.
    pub geometry_tolerance_m: f64,
    pub velocity_tolerance_m_s: f64,
}

/// Physical response and both original impulse receipts, with no new time step.
#[derive(Clone, Debug)]
pub struct RodRigidImpactReport {
    pub contact: NormalImpactResult,
    pub gap_m: f64,
    pub unit_normal_world: [f64; 3],
    pub rod_inverse_mass_kg_inverse: f64,
    pub rigid_inverse_mass_kg_inverse: f64,
    pub relative_normal_speed_before_m_s: f64,
    pub relative_normal_speed_after_m_s: f64,
    pub endpoint_velocity_residual_m_s: f64,
    /// Actual pair kinetic increment, accumulated from stored momentum changes.
    /// Common translational work cancels before rounding back to f64.
    pub kinetic_energy_change_j: f64,
    pub relative_impulse_work_j: f64,
    pub work_energy_residual_j: f64,
    /// Actual kinetic change plus original restitution dissipation.
    pub energy_residual_j: f64,
    /// Actual pair momentum change minus the rod's stationary support impulse.
    pub linear_momentum_residual_n_s: [f64; 3],
    /// Same balance about `controls.impulse.balance_origin_m` [N m s].
    pub angular_momentum_residual_nm_s: [f64; 3],
    /// No impulse owner is invoked for an inactive or separating witness.
    pub rod_impulse: Option<RodImpulseReport>,
    pub rigid_impulse: Option<ImpulseReceipt>,
}

#[derive(Debug)]
pub enum RodRigidImpactError {
    Invalid(&'static str),
    Rod(RodMotionError),
    Rigid(DynamicsError),
    Contact(NormalImpactError),
    /// A native update did not resolve the requested physical event. No
    /// participant is published, even if only the final pair balance failed.
    Balance {
        velocity_defect_m_s: f64,
        linear_defect_n_s: f64,
        angular_defect_nm_s: f64,
        work_defect_j: f64,
    },
}
impl From<RodMotionError> for RodRigidImpactError {
    fn from(error: RodMotionError) -> Self {
        Self::Rod(error)
    }
}
impl From<DynamicsError> for RodRigidImpactError {
    fn from(error: DynamicsError) -> Self {
        Self::Rigid(error)
    }
}
impl From<NormalImpactError> for RodRigidImpactError {
    fn from(error: NormalImpactError) -> Self {
        Self::Contact(error)
    }
}
impl core::fmt::Display for RodRigidImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rod/rigid normal impact: {self:?}")
    }
}
impl std::error::Error for RodRigidImpactError {}

impl<I: RodMotionInertia> RodMotion<I> {
    /// Resolve one frictionless impact and return BOTH continuing native states.
    ///
    /// The rod's actual nodal centroid, full central inertia and existing fixed
    /// DOFs give its directional mobility. The finite rigid body contributes
    /// its original full-inertia point mobility. The original `fs-contact`
    /// Newton law chooses the scalar impulse, then the original native owners
    /// apply its equal/opposite actions. Stationary rod supports supply their
    /// separately reported physical reactions, without a surrogate large mass.
    ///
    /// Poses, elastic/material history, clocks and accepted-step counts remain
    /// unchanged. The returned states can resume their original time integrators.
    /// An inactive/separating contact retains both states bit for bit. Refusal
    /// or cancellation publishes neither state. Cost is O(rod nodes) for the
    /// existing impulse owner; mobility/contact preparation is constant work
    /// with no source evaluation, matrix factorization or time subdivision.
    pub fn apply_rigid_normal_impact(
        &self,
        rigid: RigidBodyState,
        rigid_inertia: MassProperties,
        witness: RodRigidNormalContact,
        controls: RodRigidImpactControls,
        gate: &CancelGate,
    ) -> Result<(Self, RigidBodyState, RodRigidImpactReport), RodRigidImpactError> {
        checkpoint(gate)?;
        if !controls.geometry_tolerance_m.is_finite()
            || controls.geometry_tolerance_m < 0.0
            || [
                controls.velocity_tolerance_m_s,
                controls.impulse.linear_tolerance_n_s,
                controls.impulse.angular_tolerance_nm_s,
                controls.impulse.work_tolerance_j,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(RodRigidImpactError::Invalid(
                "invalid impact acceptance bounds",
            ));
        }
        finite(&controls.impulse.balance_origin_m)?;
        finite(&witness.point_body_m)?;
        let node = witness.node;
        let k = *self
            .kinetics
            .get(node)
            .ok_or(RodRigidImpactError::Invalid("absent rod node"))?;
        let frame = self.pose.frames[node];
        let body_arm = vector(witness.rigid_arm_body_m);
        let body_mobility = rigid.directional_effective_mass(
            rigid_inertia,
            body_arm,
            vector(witness.normal_world),
        )?;
        let normal = components(body_mobility.unit_direction_world);
        let normal_body = rotate(frame.inverse(), normal)?;
        let central_arm = sub(witness.point_body_m, k.center);
        let torque_axis = cross(central_arm, normal_body);
        let mut rod_mobility = 0.0;
        for c in 0..3 {
            if !self.model.fixed_positions()[node][c] {
                rod_mobility += normal[c] * normal[c] / k.mass;
            }
        }
        if !self.model.fixed_frames()[node] {
            rod_mobility += dot(torque_axis, k.omega(torque_axis)?);
        }
        if !rod_mobility.is_finite() || rod_mobility < 0.0 {
            return Err(RodRigidImpactError::Invalid(
                "unrepresentable rod point mobility",
            ));
        }
        let body_before = rigid.point_kinematics(rigid_inertia, body_arm)?;
        let point = add(
            self.pose.positions_m[node],
            rotate(frame, witness.point_body_m)?,
        );
        let separation = sub(point, components(body_before.point_world));
        finite(&separation)?;
        let gap = dot(normal, separation);
        let tangent = std::array::from_fn(|c| separation[c] - gap * normal[c]);
        if norm(tangent) > controls.geometry_tolerance_m {
            return Err(RodRigidImpactError::Invalid(
                "contact witnesses must share their normal line",
            ));
        }
        let speed = self.relative_normal_speed(
            node,
            central_arm,
            rigid,
            rigid_inertia,
            witness.rigid_arm_body_m,
            normal,
        )?;
        let contact = solve_normal_impact(NormalImpactInput {
            gap,
            activation_tolerance: witness.activation_tolerance_m,
            relative_normal_speed: speed,
            inverse_effective_mass: rod_mobility + body_mobility.inverse_mass_kg_inverse,
            restitution: witness.restitution,
        })?;
        let mut report = RodRigidImpactReport {
            contact,
            gap_m: gap,
            unit_normal_world: normal,
            rod_inverse_mass_kg_inverse: rod_mobility,
            rigid_inverse_mass_kg_inverse: body_mobility.inverse_mass_kg_inverse,
            relative_normal_speed_before_m_s: speed,
            relative_normal_speed_after_m_s: speed,
            endpoint_velocity_residual_m_s: 0.0,
            kinetic_energy_change_j: 0.0,
            relative_impulse_work_j: 0.0,
            work_energy_residual_j: 0.0,
            energy_residual_j: 0.0,
            linear_momentum_residual_n_s: [0.0; 3],
            angular_momentum_residual_nm_s: [0.0; 3],
            rod_impulse: None,
            rigid_impulse: None,
        };
        if !contact.impacted {
            checkpoint(gate)?;
            return Ok((self.clone(), rigid, report));
        }
        let impulse = normal.map(|n| n * contact.normal_impulse);
        let (next_rod, rod_report) = self.apply_point_impulses(
            &[RodPointImpulse {
                node,
                point_body_m: witness.point_body_m,
                linear_world_n_s: impulse,
                couple_world_nm_s: [0.0; 3],
            }],
            controls.impulse,
            gate,
        )?;
        let rigid_report = rigid.apply_impulse_at_body_point(
            rigid_inertia,
            body_arm,
            vector(impulse.map(|x| -x)),
        )?;
        let next_rigid = rigid_report.state_after;
        checkpoint(gate)?;
        let after_speed = next_rod.relative_normal_speed(
            node,
            central_arm,
            next_rigid,
            rigid_inertia,
            witness.rigid_arm_body_m,
            normal,
        )?;
        report.relative_normal_speed_after_m_s = after_speed;
        report.endpoint_velocity_residual_m_s = after_speed - contact.relative_normal_speed_after;

        let actual = audit_native_pair(
            self,
            &next_rod,
            rigid,
            rigid_inertia,
            &rigid_report,
            &rod_report,
            node,
            impulse,
            controls.impulse,
        )?;
        let kinetic = actual.kinetic_energy_change_j;
        let relative_work = dd(contact.normal_impulse) * (dd(speed) + dd(after_speed)) * dd(0.5);
        report.kinetic_energy_change_j = kinetic.to_f64();
        report.relative_impulse_work_j = relative_work.to_f64();
        report.work_energy_residual_j = (kinetic - relative_work).to_f64();
        report.energy_residual_j = (kinetic + dd(contact.dissipated_energy)).to_f64();

        report.linear_momentum_residual_n_s = actual.linear_momentum_residual_n_s;
        report.angular_momentum_residual_nm_s = actual.angular_momentum_residual_nm_s;
        finite(&[
            report.endpoint_velocity_residual_m_s,
            report.kinetic_energy_change_j,
            report.relative_impulse_work_j,
            report.work_energy_residual_j,
            report.energy_residual_j,
        ])?;
        let linear_defect = actual.linear_defect_n_s;
        let angular_defect = actual.angular_defect_nm_s;
        let velocity_defect = report.endpoint_velocity_residual_m_s.abs();
        let work_defect = report
            .energy_residual_j
            .abs()
            .max(report.work_energy_residual_j.abs());
        if [
            linear_defect,
            angular_defect,
            velocity_defect,
            work_defect,
            report.kinetic_energy_change_j,
            report.relative_impulse_work_j,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || linear_defect > controls.impulse.linear_tolerance_n_s
            || angular_defect > controls.impulse.angular_tolerance_nm_s
            || velocity_defect > controls.velocity_tolerance_m_s
            || work_defect > controls.impulse.work_tolerance_j
        {
            return Err(RodRigidImpactError::Balance {
                velocity_defect_m_s: velocity_defect,
                linear_defect_n_s: linear_defect,
                angular_defect_nm_s: angular_defect,
                work_defect_j: work_defect,
            });
        }
        report.rod_impulse = Some(rod_report);
        report.rigid_impulse = Some(rigid_report);
        checkpoint(gate)?;
        Ok((next_rod, next_rigid, report))
    }

    fn relative_normal_speed(
        &self,
        node: usize,
        central_arm: [f64; 3],
        rigid: RigidBodyState,
        inertia: MassProperties,
        rigid_arm: [f64; 3],
        normal: [f64; 3],
    ) -> Result<f64, RodRigidImpactError> {
        let k = self.kinetics[node];
        let rod_rotational = rotate(
            self.pose.frames[node],
            cross(k.omega(self.spin[node])?, central_arm),
        )?;
        let rigid_omega = inertia.angular_velocity_body_checked(rigid.angular_momentum_body())?;
        let rigid_rotational = components(
            rigid
                .pose()
                .orientation()
                .rotate_body_to_world(rigid_omega.cross(vector(rigid_arm))),
        );
        let rigid_p = components(rigid.linear_momentum_world());
        // Subtract common translation before adding rotation: constructing two
        // large absolute point velocities would erase small angular motion.
        let mut speed = dd(0.0);
        for c in 0..3 {
            speed = speed
                + dd(normal[c])
                    * (dd(self.linear[node][c]) / dd(k.mass) - dd(rigid_p[c]) / dd(inertia.mass())
                        + dd(rod_rotational[c])
                        - dd(rigid_rotational[c]));
        }
        let speed = speed.to_f64();
        finite(&[speed])?;
        Ok(speed)
    }
}

/// Law-neutral balance of the two ORIGINAL native impulse updates. Normal and
/// frictional contacts share actual stored increments, eccentric support moments
/// and compensated work; no commanded-impulse loss identity replaces this audit.
struct NativePairAudit {
    kinetic_energy_change_j: Dd,
    linear_momentum_residual_n_s: [f64; 3],
    angular_momentum_residual_nm_s: [f64; 3],
    linear_defect_n_s: f64,
    angular_defect_nm_s: f64,
}

#[allow(clippy::too_many_arguments)]
fn audit_native_pair<I: RodMotionInertia>(
    rod: &RodMotion<I>,
    next_rod: &RodMotion<I>,
    rigid: RigidBodyState,
    rigid_inertia: MassProperties,
    rigid_report: &ImpulseReceipt,
    rod_report: &RodImpulseReport,
    node: usize,
    impulse: [f64; 3],
    controls: RodImpulseControls,
) -> Result<NativePairAudit, RodRigidImpactError> {
    let k = rod.kinetics[node];
    let frame = rod.pose.frames[node];
    let next_rigid = rigid_report.state_after;
    let rod_dp = difference(next_rod.linear[node], rod.linear[node]);
    let rod_ds = sub(next_rod.spin[node], rod.spin[node]);
    let body_p = components(rigid.linear_momentum_world());
    let next_body_p = components(next_rigid.linear_momentum_world());
    let body_dp = difference(next_body_p, body_p);
    let body_spin = components(rigid.angular_momentum_body());
    let next_body_spin = components(next_rigid.angular_momentum_body());
    let body_ds = sub(next_body_spin, body_spin);
    let old_omega = k.omega(rod.spin[node])?;
    let next_omega = k.omega(next_rod.spin[node])?;
    let body_omega = components(rigid_report.point_kinematics_before.angular_velocity_body);
    let next_body_omega = components(rigid_report.point_kinematics_after.angular_velocity_body);
    let mut kinetic = dd(0.0);
    for c in 0..3 {
        kinetic = kinetic
            + rod_dp[c] * (dd(rod.linear[node][c]) + dd(next_rod.linear[node][c])) * dd(0.5)
                / dd(k.mass)
            + body_dp[c] * (dd(body_p[c]) + dd(next_body_p[c])) * dd(0.5)
                / dd(rigid_inertia.mass())
            + (dd(next_rod.spin[node][c]) - dd(rod.spin[node][c]))
                * (dd(old_omega[c]) + dd(next_omega[c]))
                * dd(0.5)
            + (dd(next_body_spin[c]) - dd(body_spin[c]))
                * (dd(body_omega[c]) + dd(next_body_omega[c]))
                * dd(0.5);
    }
    let support = &rod_report.support_impulses[node];
    let linear =
        std::array::from_fn(|c| (rod_dp[c] + body_dp[c] - dd(support.linear_n_s[c])).to_f64());
    let origin = controls.balance_origin_m;
    let node_arm = sub(rod.pose.positions_m[node], origin);
    let centroid_arm = add(node_arm, rotate(frame, k.center)?);
    let rigid_arm = sub(components(rigid.pose().position_world()), origin);
    let rod_spin_world = rotate(frame, rod_ds)?;
    let rigid_spin_world = components(
        rigid
            .pose()
            .orientation()
            .rotate_body_to_world(vector(body_ds)),
    );
    let rod_orbital = cross_dd(centroid_arm, rod_dp);
    let body_orbital = cross_dd(rigid_arm, body_dp);
    let support_orbital = cross_dd(node_arm, support.linear_n_s.map(dd));
    let angular = std::array::from_fn(|c| {
        (rod_orbital[c] + body_orbital[c] + dd(rod_spin_world[c]) + dd(rigid_spin_world[c])
            - support_orbital[c]
            - dd(support.angular_nm_s[c]))
        .to_f64()
    });
    let body_linear_error = std::array::from_fn(|c| (body_dp[c] + dd(impulse[c])).to_f64());
    let expected_spin = components(rigid_report.angular_impulse_body);
    let body_spin_error = std::array::from_fn(|c| {
        (dd(next_body_spin[c]) - dd(body_spin[c]) - dd(expected_spin[c])).to_f64()
    });
    finite(&linear)?;
    finite(&angular)?;
    finite(&body_linear_error)?;
    finite(&body_spin_error)?;
    finite(&[kinetic.to_f64()])?;
    Ok(NativePairAudit {
        kinetic_energy_change_j: kinetic,
        linear_momentum_residual_n_s: linear,
        angular_momentum_residual_nm_s: angular,
        linear_defect_n_s: norm(linear).max(norm(body_linear_error)),
        angular_defect_nm_s: norm(angular).max(norm(body_spin_error)),
    })
}

fn dd(x: f64) -> Dd {
    Dd::from_f64(x)
}
fn difference(a: [f64; 3], b: [f64; 3]) -> [Dd; 3] {
    std::array::from_fn(|c| dd(a[c]) - dd(b[c]))
}
fn cross_dd(a: [f64; 3], b: [Dd; 3]) -> [Dd; 3] {
    [
        dd(a[1]) * b[2] - dd(a[2]) * b[1],
        dd(a[2]) * b[0] - dd(a[0]) * b[2],
        dd(a[0]) * b[1] - dd(a[1]) * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (dd(a[0]) * dd(b[0]) + dd(a[1]) * dd(b[1]) + dd(a[2]) * dd(b[2])).to_f64()
}
fn norm(a: [f64; 3]) -> f64 {
    fs_math::det::hypot(fs_math::det::hypot(a[0], a[1]), a[2])
}
fn vector(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn components(a: Vec3) -> [f64; 3] {
    [a.x, a.y, a.z]
}
