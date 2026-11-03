//! Reciprocal impacts between native rod attachments and a finite body.
//!
//! The collision owner supplies an oriented pair of material witnesses. This
//! frozen-pose event uses the existing contact laws and both participants'
//! original impulse owners. It does not advance time, repair penetration or
//! turn restitution/friction loss into thermal heat.

/// Coupled normal and two-direction Coulomb impulses at one original attachment.
pub mod friction;

/// Simultaneous normal impacts sharing the original rod and finite rigid body.
pub mod multiple;

/// Simultaneous normal and tangential impacts on the same native participants.
pub mod friction_multiple;

/// Mixed contact networks sharing the native rod and multiple finite bodies.
pub mod graph;

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

/// Actual pair changes after the original impulse owners have staged a response.
/// The same audit applies to scalar, frictional and aggregate point impulses.
#[derive(Clone, Copy, Debug, Default)]
pub struct RodRigidImpulseBalance {
    pub kinetic_energy_change_j: f64,
    pub relative_impulse_work_j: f64,
    pub work_energy_residual_j: f64,
    pub energy_residual_j: f64,
    pub linear_momentum_residual_n_s: [f64; 3],
    pub angular_momentum_residual_nm_s: [f64; 3],
}

pub(super) struct PreparedRodRigidPoint {
    pub(super) node: usize,
    pub(super) normal_world: [f64; 3],
    pub(super) central_arm_body_m: [f64; 3],
    pub(super) gap_m: f64,
    pub(super) separation_world_m: [f64; 3],
    pub(super) rod_inverse_mass_kg_inverse: f64,
    pub(super) rigid_inverse_mass_kg_inverse: f64,
    pub(super) relative_velocity_world_m_s: [Dd; 3],
}

/// Constitutive targets and aggregate native body response, without a second
/// mechanics update. Relative work uses the original impulse/velocity pairing.
pub(super) struct RigidImpulseAudit {
    pub(super) before: RigidBodyState,
    pub(super) after: RigidBodyState,
    pub(super) inertia: MassProperties,
    pub(super) expected_linear_world_n_s: [f64; 3],
    pub(super) expected_angular_body_nm_s: [f64; 3],
    pub(super) relative_work_j: Dd,
    pub(super) dissipated_energy_j: f64,
    pub(super) velocity_defect_m_s: f64,
}

/// One original body's actual update and the resultant requested from its owner.
/// Keeping these separate per participant prevents opposing unresolved kicks
/// from disappearing in the assembled momentum balance.
pub(super) struct RigidMomentumTransfer {
    pub(super) before: RigidBodyState,
    pub(super) after: RigidBodyState,
    pub(super) inertia: MassProperties,
    pub(super) expected_linear_world_n_s: [Dd; 3],
    pub(super) expected_angular_body_nm_s: [Dd; 3],
}

/// Work belongs to the complete simultaneous event, while stored momentum
/// resolution must be checked independently for every original participant.
pub(super) struct RodRigidSystemAudit<'a> {
    pub(super) bodies: &'a [RigidMomentumTransfer],
    pub(super) relative_work_j: Dd,
    pub(super) dissipated_energy_j: f64,
    pub(super) velocity_defect_m_s: f64,
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
        write!(f, "rod/rigid impact: {self:?}")
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
        let prepared = self.prepare_rigid_point(rigid, rigid_inertia, witness, controls, gate)?;
        let node = prepared.node;
        let normal = prepared.normal_world;
        let gap = prepared.gap_m;
        let speed = project_velocity(normal, prepared.relative_velocity_world_m_s)?;
        let contact = solve_normal_impact(NormalImpactInput {
            gap,
            activation_tolerance: witness.activation_tolerance_m,
            relative_normal_speed: speed,
            inverse_effective_mass: prepared.rod_inverse_mass_kg_inverse
                + prepared.rigid_inverse_mass_kg_inverse,
            restitution: witness.restitution,
        })?;
        let mut report = RodRigidImpactReport {
            contact,
            gap_m: gap,
            unit_normal_world: normal,
            rod_inverse_mass_kg_inverse: prepared.rod_inverse_mass_kg_inverse,
            rigid_inverse_mass_kg_inverse: prepared.rigid_inverse_mass_kg_inverse,
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
            vector(witness.rigid_arm_body_m),
            vector(impulse.map(|x| -x)),
        )?;
        let next_rigid = rigid_report.state_after;
        checkpoint(gate)?;
        let after_speed = project_velocity(
            normal,
            next_rod.relative_point_velocity(
                node,
                prepared.central_arm_body_m,
                next_rigid,
                rigid_inertia,
                witness.rigid_arm_body_m,
            )?,
        )?;
        report.relative_normal_speed_after_m_s = after_speed;
        report.endpoint_velocity_residual_m_s = after_speed - contact.relative_normal_speed_after;
        let relative_work = dd(contact.normal_impulse) * (dd(speed) + dd(after_speed)) * dd(0.5);
        let balance = self.audit_rigid_impulse(
            &next_rod,
            &rod_report,
            RigidImpulseAudit {
                before: rigid,
                after: next_rigid,
                inertia: rigid_inertia,
                expected_linear_world_n_s: impulse.map(|x| -x),
                expected_angular_body_nm_s: components(rigid_report.angular_impulse_body),
                relative_work_j: relative_work,
                dissipated_energy_j: contact.dissipated_energy,
                velocity_defect_m_s: report.endpoint_velocity_residual_m_s.abs(),
            },
            controls,
            gate,
        )?;
        report.kinetic_energy_change_j = balance.kinetic_energy_change_j;
        report.relative_impulse_work_j = balance.relative_impulse_work_j;
        report.work_energy_residual_j = balance.work_energy_residual_j;
        report.energy_residual_j = balance.energy_residual_j;
        report.linear_momentum_residual_n_s = balance.linear_momentum_residual_n_s;
        report.angular_momentum_residual_nm_s = balance.angular_momentum_residual_nm_s;
        report.rod_impulse = Some(rod_report);
        report.rigid_impulse = Some(rigid_report);
        checkpoint(gate)?;
        Ok((next_rod, next_rigid, report))
    }

    pub(super) fn prepare_rigid_point(
        &self,
        rigid: RigidBodyState,
        inertia: MassProperties,
        witness: RodRigidNormalContact,
        controls: RodRigidImpactControls,
        gate: &CancelGate,
    ) -> Result<PreparedRodRigidPoint, RodRigidImpactError> {
        checkpoint(gate)?;
        validate_controls(controls)?;
        finite(&witness.point_body_m)?;
        let node = witness.node;
        let k = *self
            .kinetics
            .get(node)
            .ok_or(RodRigidImpactError::Invalid("absent rod node"))?;
        let frame = self.pose.frames[node];
        let body_arm = vector(witness.rigid_arm_body_m);
        let body_mobility =
            rigid.directional_effective_mass(inertia, body_arm, vector(witness.normal_world))?;
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
        let body_before = rigid.point_kinematics(inertia, body_arm)?;
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
        let relative_velocity = self.relative_point_velocity(
            node,
            central_arm,
            rigid,
            inertia,
            witness.rigid_arm_body_m,
        )?;
        Ok(PreparedRodRigidPoint {
            node,
            normal_world: normal,
            central_arm_body_m: central_arm,
            gap_m: gap,
            separation_world_m: separation,
            rod_inverse_mass_kg_inverse: rod_mobility,
            rigid_inverse_mass_kg_inverse: body_mobility.inverse_mass_kg_inverse,
            relative_velocity_world_m_s: relative_velocity,
        })
    }

    /// Audit the actual stored increments of every affected native node and
    /// the shared finite body. Aggregate events must pass their SUMMED body
    /// impulse, rather than sequentially publishing per-contact body states.
    pub(super) fn audit_rigid_impulse(
        &self,
        next_rod: &Self,
        rod_report: &RodImpulseReport,
        body: RigidImpulseAudit,
        controls: RodRigidImpactControls,
        gate: &CancelGate,
    ) -> Result<RodRigidImpulseBalance, RodRigidImpactError> {
        self.audit_rigid_system_impulse(
            next_rod,
            Some(rod_report),
            RodRigidSystemAudit {
                bodies: &[RigidMomentumTransfer {
                    before: body.before,
                    after: body.after,
                    inertia: body.inertia,
                    expected_linear_world_n_s: body.expected_linear_world_n_s.map(dd),
                    expected_angular_body_nm_s: body.expected_angular_body_nm_s.map(dd),
                }],
                relative_work_j: body.relative_work_j,
                dissipated_energy_j: body.dissipated_energy_j,
                velocity_defect_m_s: body.velocity_defect_m_s,
            },
            controls,
            gate,
        )
    }

    /// Audit one native rod update and any number of original finite bodies.
    /// A rigid-only contact graph can omit the rod receipt only when the rod's
    /// momenta are unchanged. There is one kinetic-work account for the event,
    /// independent transfer bounds per body, and no pairwise double counting.
    pub(super) fn audit_rigid_system_impulse(
        &self,
        next_rod: &Self,
        rod_report: Option<&RodImpulseReport>,
        audit: RodRigidSystemAudit<'_>,
        controls: RodRigidImpactControls,
        gate: &CancelGate,
    ) -> Result<RodRigidImpulseBalance, RodRigidImpactError> {
        checkpoint(gate)?;
        if rod_report.is_none() && (self.linear != next_rod.linear || self.spin != next_rod.spin) {
            return Err(RodRigidImpactError::Invalid(
                "changed rod momentum requires its native impulse receipt",
            ));
        }
        let mut kinetic = dd(0.0);
        let mut linear = [dd(0.0); 3];
        let mut angular = [dd(0.0); 3];
        let origin = controls.impulse.balance_origin_m;
        if let Some(rod_report) = rod_report {
            if rod_report.support_impulses.len() != self.linear.len() {
                return Err(RodRigidImpactError::Invalid(
                    "incomplete rod support receipt",
                ));
            }
            for node in 0..self.linear.len() {
                checkpoint(gate)?;
                let support = &rod_report.support_impulses[node];
                if self.linear[node] == next_rod.linear[node]
                    && self.spin[node] == next_rod.spin[node]
                    && support.linear_n_s == [0.0; 3]
                    && support.angular_nm_s == [0.0; 3]
                {
                    continue;
                }
                let k = self.kinetics[node];
                let frame = self.pose.frames[node];
                let dp = difference(next_rod.linear[node], self.linear[node]);
                let ds = difference(next_rod.spin[node], self.spin[node]);
                let old_omega = k.omega(self.spin[node])?;
                let next_omega = k.omega(next_rod.spin[node])?;
                for c in 0..3 {
                    kinetic = kinetic
                        + dp[c]
                            * (dd(self.linear[node][c]) * dd(0.5)
                                + dd(next_rod.linear[node][c]) * dd(0.5))
                            / dd(k.mass)
                        + ds[c] * (dd(old_omega[c]) * dd(0.5) + dd(next_omega[c]) * dd(0.5));
                    linear[c] = linear[c] + dp[c] - dd(support.linear_n_s[c]);
                }
                let node_arm = sub(self.pose.positions_m[node], origin);
                let centroid_arm = add(node_arm, rotate(frame, k.center)?);
                let orbital = cross_dd(centroid_arm, dp);
                let spin_world = rotate(frame, ds.map(Dd::to_f64))?;
                let support_orbital = cross_dd(node_arm, support.linear_n_s.map(dd));
                for c in 0..3 {
                    angular[c] = angular[c] + orbital[c] + dd(spin_world[c])
                        - support_orbital[c]
                        - dd(support.angular_nm_s[c]);
                }
            }
        }
        let mut body_linear_defect: f64 = 0.0;
        let mut body_angular_defect: f64 = 0.0;
        for body in audit.bodies {
            checkpoint(gate)?;
            let p = components(body.before.linear_momentum_world());
            let next_p = components(body.after.linear_momentum_world());
            let spin = components(body.before.angular_momentum_body());
            let next_spin = components(body.after.angular_momentum_body());
            let dp = difference(next_p, p);
            let ds = difference(next_spin, spin);
            let omega = components(
                body.inertia
                    .angular_velocity_body_checked(body.before.angular_momentum_body())?,
            );
            let next_omega = components(
                body.inertia
                    .angular_velocity_body_checked(body.after.angular_momentum_body())?,
            );
            let arm = sub(components(body.before.pose().position_world()), origin);
            let orbital = cross_dd(arm, dp);
            let spin_world = components(
                body.before
                    .pose()
                    .orientation()
                    .rotate_body_to_world(vector(ds.map(Dd::to_f64))),
            );
            for c in 0..3 {
                kinetic = kinetic
                    + dp[c] * (dd(p[c]) * dd(0.5) + dd(next_p[c]) * dd(0.5))
                        / dd(body.inertia.mass())
                    + ds[c] * (dd(omega[c]) * dd(0.5) + dd(next_omega[c]) * dd(0.5));
                linear[c] = linear[c] + dp[c];
                angular[c] = angular[c] + orbital[c] + dd(spin_world[c]);
            }
            let body_linear_error =
                std::array::from_fn(|c| (dp[c] - body.expected_linear_world_n_s[c]).to_f64());
            let body_spin_error =
                std::array::from_fn(|c| (ds[c] - body.expected_angular_body_nm_s[c]).to_f64());
            finite(&body_linear_error)?;
            finite(&body_spin_error)?;
            body_linear_defect = body_linear_defect.max(norm(body_linear_error));
            body_angular_defect = body_angular_defect.max(norm(body_spin_error));
        }
        let report = RodRigidImpulseBalance {
            kinetic_energy_change_j: kinetic.to_f64(),
            relative_impulse_work_j: audit.relative_work_j.to_f64(),
            work_energy_residual_j: (kinetic - audit.relative_work_j).to_f64(),
            energy_residual_j: (kinetic + dd(audit.dissipated_energy_j)).to_f64(),
            linear_momentum_residual_n_s: linear.map(Dd::to_f64),
            angular_momentum_residual_nm_s: angular.map(Dd::to_f64),
        };
        finite(&report.linear_momentum_residual_n_s)?;
        finite(&report.angular_momentum_residual_nm_s)?;
        finite(&[
            audit.velocity_defect_m_s,
            report.kinetic_energy_change_j,
            report.relative_impulse_work_j,
            report.work_energy_residual_j,
            report.energy_residual_j,
        ])?;
        let linear_defect = norm(report.linear_momentum_residual_n_s).max(body_linear_defect);
        let angular_defect = norm(report.angular_momentum_residual_nm_s).max(body_angular_defect);
        let velocity_defect = audit.velocity_defect_m_s.abs();
        let work_defect = report
            .energy_residual_j
            .abs()
            .max(report.work_energy_residual_j.abs());
        if [linear_defect, angular_defect, velocity_defect, work_defect]
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
        checkpoint(gate)?;
        Ok(report)
    }

    pub(super) fn relative_point_velocity(
        &self,
        node: usize,
        central_arm: [f64; 3],
        rigid: RigidBodyState,
        inertia: MassProperties,
        rigid_arm: [f64; 3],
    ) -> Result<[Dd; 3], RodRigidImpactError> {
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
        let velocity: [Dd; 3] = std::array::from_fn(|c| {
            dd(self.linear[node][c]) / dd(k.mass) - dd(rigid_p[c]) / dd(inertia.mass())
                + dd(rod_rotational[c])
                - dd(rigid_rotational[c])
        });
        finite(&velocity.map(Dd::to_f64))?;
        Ok(velocity)
    }
}

fn dd(x: f64) -> Dd {
    Dd::from_f64(x)
}
pub(super) fn validate_controls(
    controls: RodRigidImpactControls,
) -> Result<(), RodRigidImpactError> {
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
    Ok(())
}
pub(super) fn project_velocity(
    axis: [f64; 3],
    velocity: [Dd; 3],
) -> Result<f64, RodRigidImpactError> {
    let speed = (dd(axis[0]) * velocity[0] + dd(axis[1]) * velocity[1] + dd(axis[2]) * velocity[2])
        .to_f64();
    finite(&[speed])?;
    Ok(speed)
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
