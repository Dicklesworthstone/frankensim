//! Simultaneous normal impacts on one native rod and one finite rigid body.
//!
//! The original contact owner solves the complete coupled mobility. All rod
//! attachments and the shared rigid body receive one aggregate native update;
//! an impulse representative is never mistaken for uniquely determined local
//! contact loads when redundant witnesses determine only the physical motion.

use super::super::impulse::{RodImpulseReport, RodPointImpulse};
use super::super::{RodMotion, RodMotionError, RodMotionInertia, checkpoint};
use super::{
    RigidImpulseAudit, RodRigidImpactControls, RodRigidImpactError, RodRigidImpulseBalance,
    RodRigidNormalContact, components, cross, dd, dot, finite, project_velocity, rotate,
    validate_controls, vector,
};
use fs_contact::multicontact::{
    MultiImpactConfig, MultiImpactError, MultiImpactResult, NormalContactImpact,
    solve_normal_impacts,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_mbd::{DynamicsError, MassProperties, MomentumImpulseReceipt, RigidBodyState};

/// One native event and its original simultaneous-contact law receipt.
#[derive(Clone, Debug)]
pub struct RodRigidMultiImpactReport {
    /// Includes the original UniqueImpulses/UniqueVelocity classification.
    /// Redundant local impulses remain a checked representative, not a unique
    /// pressure or heat partition. The aggregate native momentum is checked.
    pub contact: MultiImpactResult,
    /// Complete row-major contact mobility J M^-1 J^T [1/kg].
    pub delassus_kg_inverse: Vec<f64>,
    /// Endpoint speeds reconstructed from both actually stored native states.
    pub relative_normal_speeds_after_m_s: Vec<f64>,
    pub rod_impulse: Option<RodImpulseReport>,
    pub rigid_impulse: Option<MomentumImpulseReceipt>,
    /// No native impulse owner runs when every accepted impulse is zero.
    pub balance: Option<RodRigidImpulseBalance>,
}

#[derive(Debug)]
pub enum RodRigidMultiImpactError {
    Invalid(&'static str),
    Native(RodRigidImpactError),
    Contact(MultiImpactError),
}
impl From<RodRigidImpactError> for RodRigidMultiImpactError {
    fn from(error: RodRigidImpactError) -> Self {
        Self::Native(error)
    }
}
impl From<RodMotionError> for RodRigidMultiImpactError {
    fn from(error: RodMotionError) -> Self {
        Self::Native(error.into())
    }
}
impl From<DynamicsError> for RodRigidMultiImpactError {
    fn from(error: DynamicsError) -> Self {
        Self::Native(error.into())
    }
}
impl From<MultiImpactError> for RodRigidMultiImpactError {
    fn from(error: MultiImpactError) -> Self {
        Self::Contact(error)
    }
}
impl core::fmt::Display for RodRigidMultiImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "simultaneous rod/rigid normal impact: {self:?}")
    }
}
impl std::error::Error for RodRigidMultiImpactError {}

struct AngularResponse {
    rod_torque: [f64; 3],
    rod_velocity: [f64; 3],
    body_torque: [f64; 3],
    body_velocity: [f64; 3],
}

impl<I: RodMotionInertia> RodMotion<I> {
    /// Resolve all supplied normal contacts at one unchanged configuration.
    ///
    /// Every row sees the SAME finite rigid body and the original rod nodal
    /// inertia. Attachments on a shared rod node include their translational
    /// and rotational cross response; attachments on distinct nodes remain
    /// coupled through the finite body's full tensor inertia. Geometry-active
    /// rows include initially separating contacts, which another impact may
    /// otherwise drive into closing motion.
    ///
    /// The original normal-contact law admits restitution, complementarity,
    /// global relative work and determinacy. Its complete impulse vector is
    /// applied once through each original native owner. Actual endpoint flow,
    /// support-inclusive linear/angular momentum and pair kinetic work must
    /// pass before either continuing state is returned. No sequential impact
    /// order, penetration repair or changed restitution is inferred.
    ///
    /// The existing contact solver caps the set at 64 rows; its variable and
    /// iteration caps apply before allocation/solve. Preparation uses O(n^2)
    /// work/storage and only O(n) prelowered inverse-inertia actions. Native
    /// publication costs O(rod nodes + n), with one rod clone, no per-row state
    /// copies, no source evaluation and no new eigen or mass factorization.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_rigid_normal_impacts(
        &self,
        rigid: RigidBodyState,
        properties: MassProperties,
        witnesses: &[RodRigidNormalContact],
        controls: RodRigidImpactControls,
        contact_config: MultiImpactConfig,
        gate: &CancelGate,
    ) -> Result<(Self, RigidBodyState, RodRigidMultiImpactReport), RodRigidMultiImpactError> {
        checkpoint(gate)?;
        let count = witnesses.len();
        if count > 64 || count > contact_config.solver.max_variables {
            return Err(RodRigidMultiImpactError::Invalid(
                "normal contact count budget",
            ));
        }
        // Empty sets still admit the same physical acceptance bounds.
        validate_controls(controls)?;

        let mut points = Vec::with_capacity(count);
        let mut responses = Vec::with_capacity(count);
        let mut contacts = Vec::with_capacity(count);
        for &witness in witnesses {
            checkpoint(gate)?;
            let point = self.prepare_rigid_point(rigid, properties, witness, controls, gate)?;
            let node = point.node;
            let direction_body = rotate(self.pose.frames[node].inverse(), point.normal_world)?;
            let rod_torque = cross(point.central_arm_body_m, direction_body);
            finite(&rod_torque)?;
            let rod_velocity = if self.model.fixed_frames()[node] {
                [0.0; 3]
            } else {
                self.kinetics[node].omega(rod_torque)?
            };
            let body_direction = components(
                rigid
                    .pose()
                    .orientation()
                    .rotate_world_to_body(vector(point.normal_world)),
            );
            let body_torque = cross(witness.rigid_arm_body_m, body_direction);
            let body_velocity =
                components(properties.angular_velocity_body_checked(vector(body_torque))?);
            contacts.push(NormalContactImpact {
                gap: point.gap_m,
                activation_tolerance: witness.activation_tolerance_m,
                relative_normal_speed: project_velocity(
                    point.normal_world,
                    point.relative_velocity_world_m_s,
                )?,
                restitution: witness.restitution,
            });
            responses.push(AngularResponse {
                rod_torque,
                rod_velocity,
                body_torque,
                body_velocity,
            });
            points.push(point);
        }

        let mut delassus = vec![0.0; count * count];
        for i in 0..count {
            checkpoint(gate)?;
            for j in i..count {
                let a = &points[i];
                let b = &points[j];
                let ra = &responses[i];
                let rb = &responses[j];
                // Symmetric bilinear evaluation of the SAME original inverse
                // inertia, without replacing it by a diagonal point mass.
                let mut value = dd(dot(a.normal_world, b.normal_world)) / dd(properties.mass())
                    + dd(0.5)
                        * (dd(dot(ra.body_torque, rb.body_velocity))
                            + dd(dot(rb.body_torque, ra.body_velocity)));
                if a.node == b.node {
                    let mass = self.kinetics[a.node].mass;
                    for c in 0..3 {
                        if !self.model.fixed_positions()[a.node][c] {
                            value =
                                value + dd(a.normal_world[c]) * dd(b.normal_world[c]) / dd(mass);
                        }
                    }
                    value = value
                        + dd(0.5)
                            * (dd(dot(ra.rod_torque, rb.rod_velocity))
                                + dd(dot(rb.rod_torque, ra.rod_velocity)));
                }
                let value = value.to_f64();
                finite(&[value])?;
                delassus[i * count + j] = value;
                delassus[j * count + i] = value;
            }
        }
        let contact =
            solve_normal_impacts(&contacts, &delassus, contact_config, || gate.is_requested())?;
        if contact.normal_impulses.iter().all(|&p| p == 0.0) {
            let speeds = contacts.iter().map(|p| p.relative_normal_speed).collect();
            checkpoint(gate)?;
            return Ok((
                self.clone(),
                rigid,
                RodRigidMultiImpactReport {
                    contact,
                    delassus_kg_inverse: delassus,
                    relative_normal_speeds_after_m_s: speeds,
                    rod_impulse: None,
                    rigid_impulse: None,
                    balance: None,
                },
            ));
        }

        let mut rod_impulses = Vec::with_capacity(count);
        let mut body_linear = [Dd::ZERO; 3];
        let mut body_angular = [Dd::ZERO; 3];
        for (i, &p) in contact.normal_impulses.iter().enumerate() {
            checkpoint(gate)?;
            if p == 0.0 {
                continue;
            }
            let impulse = points[i].normal_world.map(|x| x * p);
            finite(&impulse)?;
            // Form the body's requested couple from the SAME rounded physical
            // point impulse that enters the rod owner, then aggregate once.
            let local = components(
                rigid
                    .pose()
                    .orientation()
                    .rotate_world_to_body(vector(impulse)),
            );
            let torque = super::cross_dd(witnesses[i].rigid_arm_body_m, local.map(dd));
            for c in 0..3 {
                body_linear[c] = body_linear[c] - dd(impulse[c]);
                body_angular[c] = body_angular[c] - torque[c];
            }
            rod_impulses.push(RodPointImpulse {
                node: witnesses[i].node,
                point_body_m: witnesses[i].point_body_m,
                linear_world_n_s: impulse,
                couple_world_nm_s: [0.0; 3],
            });
        }
        let body_linear = body_linear.map(Dd::to_f64);
        let body_angular = body_angular.map(Dd::to_f64);
        finite(&body_linear)?;
        finite(&body_angular)?;
        let (next_rod, rod_report) =
            self.apply_point_impulses(&rod_impulses, controls.impulse, gate)?;
        let rigid_report =
            rigid.apply_momentum_impulse(properties, vector(body_linear), vector(body_angular))?;
        let next_rigid = rigid_report.state_after;
        let mut speeds = Vec::with_capacity(count);
        let mut velocity_defect: f64 = 0.0;
        let mut relative_work = Dd::ZERO;
        for i in 0..count {
            checkpoint(gate)?;
            let after = next_rod.relative_point_velocity(
                points[i].node,
                points[i].central_arm_body_m,
                next_rigid,
                properties,
                witnesses[i].rigid_arm_body_m,
            )?;
            let after = project_velocity(points[i].normal_world, after)?;
            velocity_defect =
                velocity_defect.max((after - contact.relative_normal_speeds_after[i]).abs());
            relative_work = relative_work
                + dd(contact.normal_impulses[i])
                    * (dd(contacts[i].relative_normal_speed) + dd(after))
                    * dd(0.5);
            speeds.push(after);
        }
        let balance = self.audit_rigid_impulse(
            &next_rod,
            &rod_report,
            RigidImpulseAudit {
                before: rigid,
                after: next_rigid,
                inertia: properties,
                expected_linear_world_n_s: body_linear,
                expected_angular_body_nm_s: body_angular,
                relative_work_j: relative_work,
                dissipated_energy_j: contact.dissipated_energy,
                velocity_defect_m_s: velocity_defect,
            },
            controls,
            gate,
        )?;
        checkpoint(gate)?;
        Ok((
            next_rod,
            next_rigid,
            RodRigidMultiImpactReport {
                contact,
                delassus_kg_inverse: delassus,
                relative_normal_speeds_after_m_s: speeds,
                rod_impulse: Some(rod_report),
                rigid_impulse: Some(rigid_report),
                balance: Some(balance),
            },
        ))
    }
}
