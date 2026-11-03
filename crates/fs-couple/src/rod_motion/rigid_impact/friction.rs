//! Coupled Coulomb impact between an original rod attachment and a finite body.
//!
//! The complete normal/tangent mobility includes both full inertias, their
//! original point arms and the rod's admitted stationary supports. The existing
//! Signorini--Coulomb owner solves all three impulse components together. This
//! is one isolated frozen-pose event, not persistent contact or collision search.

use super::super::impulse::{RodImpulseReport, RodPointImpulse};
use super::super::{
    RodMotion, RodMotionError, RodMotionInertia, checkpoint, cross, finite, rotate,
};
use super::{
    RigidImpulseAudit, RodRigidImpactControls, RodRigidImpactError, RodRigidNormalContact,
    components, dd, dot, norm, project_velocity, vector,
};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactError, FrictionalContactResult,
    NormalVelocityLaw, solve_frictional_contacts,
};
use fs_exec::CancelGate;
use fs_mbd::{DynamicsError, ImpulseReceipt, MassProperties, RigidBodyState};

/// Original material witnesses and an explicitly selected interface coefficient.
#[derive(Clone, Copy, Debug)]
pub struct RodRigidFrictionalContact {
    pub geometry: RodRigidNormalContact,
    /// Nonzero world tangent hint, projected into the original contact plane.
    /// Changing this coordinate basis does not change isotropic Coulomb friction.
    pub tangent_world: [f64; 3],
    /// Finite, nonnegative Coulomb coefficient; no material/default is inferred.
    pub friction_coefficient: f64,
}

/// Existing physical transfer limits and the existing coupled-law solve limits.
#[derive(Clone, Copy, Debug)]
pub struct RodRigidFrictionalControls {
    pub impact: RodRigidImpactControls,
    pub contact: FrictionalContactConfig,
}

/// Native endpoints, actual work and support-inclusive world momentum balances.
#[derive(Clone, Debug)]
pub struct RodRigidFrictionalReport {
    pub contact: FrictionalContactResult,
    pub gap_m: f64,
    /// Right-handed world axes, normal first, then the two tangent directions.
    pub frame_world: [[f64; 3]; 3],
    /// Full normal/tangent point mobility of each participant [kg^-1].
    pub rod_inverse_mass_kg_inverse: [[f64; 3]; 3],
    pub rigid_inverse_mass_kg_inverse: [[f64; 3]; 3],
    pub relative_velocity_before_m_s: [f64; 3],
    pub relative_velocity_after_m_s: [f64; 3],
    /// Actual minus predicted endpoint velocities, in the same contact frame.
    pub endpoint_velocity_residual_m_s: [f64; 3],
    /// One world impulse ON the rod; the finite body receives its opposite.
    pub impulse_world_n_s: [f64; 3],
    pub kinetic_energy_change_j: f64,
    pub relative_impulse_work_j: f64,
    pub work_energy_residual_j: f64,
    /// Actual kinetic change plus TOTAL impact loss, not a friction heat split.
    pub energy_residual_j: f64,
    /// Actual combined change minus the rod's stationary support impulse.
    pub linear_momentum_residual_n_s: [f64; 3],
    /// Same balance about `controls.impact.impulse.balance_origin_m`.
    pub angular_momentum_residual_nm_s: [f64; 3],
    pub rod_impulse: Option<RodImpulseReport>,
    pub rigid_impulse: Option<ImpulseReceipt>,
}

/// A refusal publishes neither participant, including after a successful solve.
#[derive(Debug)]
pub enum RodRigidFrictionalError {
    Impact(RodRigidImpactError),
    Contact(FrictionalContactError),
}
impl From<RodRigidImpactError> for RodRigidFrictionalError {
    fn from(value: RodRigidImpactError) -> Self {
        Self::Impact(value)
    }
}
impl From<RodMotionError> for RodRigidFrictionalError {
    fn from(value: RodMotionError) -> Self {
        Self::Impact(value.into())
    }
}
impl From<DynamicsError> for RodRigidFrictionalError {
    fn from(value: DynamicsError) -> Self {
        Self::Impact(value.into())
    }
}
impl From<FrictionalContactError> for RodRigidFrictionalError {
    fn from(value: FrictionalContactError) -> Self {
        Self::Contact(value)
    }
}
impl core::fmt::Display for RodRigidFrictionalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rod/rigid Coulomb impact: {self:?}")
    }
}
impl std::error::Error for RodRigidFrictionalError {}

impl<I: RodMotionInertia> RodMotion<I> {
    /// Resolve one simultaneous normal/tangent impulse on BOTH native states.
    ///
    /// Full 3x3 mobility retains the coupling between restitution, slip and spin.
    /// The original non-associated Signorini--Coulomb law supplies the accepted
    /// vector; normal impulses are not solved first and followed by a friction
    /// clip. Both original impulse owners then apply one equal/opposite vector.
    /// The shared native pair audit checks actual momentum increments and work.
    ///
    /// A nonzero tangential impulse requires collocated original witnesses up
    /// to the geometry tolerance. Their actual separation still enters the
    /// angular balance; no point is relocated to hide a couple. Stationary rod
    /// supports retain their original physical reactions and do zero work.
    ///
    /// Poses, elastic/material histories, time and step counters stay unchanged.
    /// Returned states continue on their original integrators. A zero impulse
    /// preserves both states bit for bit. Cancellation, unresolved Coulomb
    /// determinacy, infeasible geometry or failed endpoint admission returns no
    /// candidate. Total impact loss is not an automatic thermal heat allocation.
    pub fn apply_rigid_frictional_impact(
        &self,
        rigid: RigidBodyState,
        rigid_inertia: MassProperties,
        witness: RodRigidFrictionalContact,
        controls: RodRigidFrictionalControls,
        gate: &CancelGate,
    ) -> Result<(Self, RigidBodyState, RodRigidFrictionalReport), RodRigidFrictionalError> {
        let limits = controls.impact;
        let point = witness.geometry;
        let prepared = self.prepare_rigid_point(rigid, rigid_inertia, point, limits, gate)?;
        let node = prepared.node;
        let k = self.kinetics[node];
        let pose = self.pose.frames[node];
        let arm = prepared.central_arm_body_m;
        let rigid_arm = vector(point.rigid_arm_body_m);
        let axes = contact_frame(prepared.normal_world, witness.tangent_world)?;
        let gap = prepared.gap_m;

        // J M^-1 J^T from ORIGINAL centroid and central tensors. Compute one
        // triangle and mirror it; no diagonal effective-mass approximation.
        let mut rod_torque = [[0.0; 3]; 3];
        let mut rod_rotation = [[0.0; 3]; 3];
        let mut body_torque = [[0.0; 3]; 3];
        let mut body_rotation = [[0.0; 3]; 3];
        for j in 0..3 {
            rod_torque[j] = cross(arm, rotate(pose.inverse(), axes[j])?);
            if !self.model.fixed_frames()[node] {
                rod_rotation[j] = k.omega(rod_torque[j])?;
            }
            body_torque[j] = components(
                rigid_arm.cross(
                    rigid
                        .pose()
                        .orientation()
                        .rotate_world_to_body(vector(axes[j])),
                ),
            );
            body_rotation[j] =
                components(rigid_inertia.angular_velocity_body_checked(vector(body_torque[j]))?);
        }
        let mut rod_mobility = [[0.0; 3]; 3];
        let mut body_mobility = [[0.0; 3]; 3];
        let mut mobility = [0.0; 9];
        for i in 0..3 {
            for j in i..3 {
                let mut translation = dd(0.0);
                for c in 0..3 {
                    if !self.model.fixed_positions()[node][c] {
                        translation = translation + dd(axes[i][c]) * dd(axes[j][c]) / dd(k.mass);
                    }
                }
                let rod_value = (translation + dd(dot(rod_torque[i], rod_rotation[j]))).to_f64();
                let body_value = (dd(dot(axes[i], axes[j])) / dd(rigid_inertia.mass())
                    + dd(dot(body_torque[i], body_rotation[j])))
                .to_f64();
                finite(&[rod_value, body_value])?;
                rod_mobility[i][j] = rod_value;
                rod_mobility[j][i] = rod_value;
                body_mobility[i][j] = body_value;
                body_mobility[j][i] = body_value;
                mobility[3 * i + j] = rod_value + body_value;
                mobility[3 * j + i] = mobility[3 * i + j];
            }
        }
        let mut before = [0.0; 3];
        for j in 0..3 {
            before[j] = project_velocity(axes[j], prepared.relative_velocity_world_m_s)?;
        }
        let contact = solve_frictional_contacts(
            &[FrictionalContact {
                id: 0,
                gap,
                activation_tolerance: point.activation_tolerance_m,
                relative_velocity: before,
                friction_coefficient: witness.friction_coefficient,
                normal_law: NormalVelocityLaw::NewtonRestitution {
                    coefficient: point.restitution,
                },
            }],
            &mobility,
            controls.contact,
            || gate.is_requested(),
        )?;
        let p = contact.impulses[0];
        let impulse = std::array::from_fn(|c| {
            (dd(axes[0][c]) * dd(p[0]) + dd(axes[1][c]) * dd(p[1]) + dd(axes[2][c]) * dd(p[2]))
                .to_f64()
        });
        finite(&impulse)?;
        let mut report = RodRigidFrictionalReport {
            contact,
            gap_m: gap,
            frame_world: axes,
            rod_inverse_mass_kg_inverse: rod_mobility,
            rigid_inverse_mass_kg_inverse: body_mobility,
            relative_velocity_before_m_s: before,
            relative_velocity_after_m_s: before,
            endpoint_velocity_residual_m_s: [0.0; 3],
            impulse_world_n_s: impulse,
            kinetic_energy_change_j: 0.0,
            relative_impulse_work_j: 0.0,
            work_energy_residual_j: 0.0,
            energy_residual_j: 0.0,
            linear_momentum_residual_n_s: [0.0; 3],
            angular_momentum_residual_nm_s: [0.0; 3],
            rod_impulse: None,
            rigid_impulse: None,
        };
        if p.iter().all(|&x| x == 0.0) {
            checkpoint(gate)?;
            return Ok((self.clone(), rigid, report));
        }
        if (p[1] != 0.0 || p[2] != 0.0)
            && norm(prepared.separation_world_m) > limits.geometry_tolerance_m
        {
            return Err(invalid(
                "tangential impulses require collocated original witnesses",
            ));
        }
        let (next_rod, rod_report) = self.apply_point_impulses(
            &[RodPointImpulse {
                node,
                point_body_m: point.point_body_m,
                linear_world_n_s: impulse,
                couple_world_nm_s: [0.0; 3],
            }],
            limits.impulse,
            gate,
        )?;
        let rigid_report = rigid.apply_impulse_at_body_point(
            rigid_inertia,
            rigid_arm,
            vector(impulse.map(|x| -x)),
        )?;
        let next_rigid = rigid_report.state_after;
        checkpoint(gate)?;
        let after_world = next_rod.relative_point_velocity(
            node,
            arm,
            next_rigid,
            rigid_inertia,
            point.rigid_arm_body_m,
        )?;
        let mut relative_work = dd(0.0);
        for j in 0..3 {
            let after = project_velocity(axes[j], after_world)?;
            report.relative_velocity_after_m_s[j] = after;
            report.endpoint_velocity_residual_m_s[j] =
                after - report.contact.relative_velocities_after[0][j];
            relative_work = relative_work + dd(p[j]) * (dd(before[j]) + dd(after)) * dd(0.5);
        }
        let velocity_defect = norm(report.endpoint_velocity_residual_m_s);
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
                dissipated_energy_j: report.contact.dissipated_energy,
                velocity_defect_m_s: velocity_defect,
            },
            limits,
            gate,
        )?;
        report.kinetic_energy_change_j = balance.kinetic_energy_change_j;
        report.relative_impulse_work_j = balance.relative_impulse_work_j;
        report.work_energy_residual_j = balance.work_energy_residual_j;
        report.energy_residual_j = balance.energy_residual_j;
        report.linear_momentum_residual_n_s = balance.linear_momentum_residual_n_s;
        report.angular_momentum_residual_nm_s = balance.angular_momentum_residual_nm_s;
        // Retain the published independent check against the constitutive
        // owner's work, in addition to the shared actual native pair audit.
        let law_work_defect = (relative_work - dd(report.contact.relative_work))
            .to_f64()
            .abs();
        finite(&[law_work_defect])?;
        if law_work_defect > limits.impulse.work_tolerance_j {
            return Err(RodRigidImpactError::Balance {
                velocity_defect_m_s: velocity_defect,
                linear_defect_n_s: norm(balance.linear_momentum_residual_n_s),
                angular_defect_nm_s: norm(balance.angular_momentum_residual_nm_s),
                work_defect_j: law_work_defect,
            }
            .into());
        }
        report.rod_impulse = Some(rod_report);
        report.rigid_impulse = Some(rigid_report);
        checkpoint(gate)?;
        Ok((next_rod, next_rigid, report))
    }
}

fn invalid(reason: &'static str) -> RodRigidFrictionalError {
    RodRigidImpactError::Invalid(reason).into()
}

pub(super) fn contact_frame(
    normal: [f64; 3],
    hint: [f64; 3],
) -> Result<[[f64; 3]; 3], RodRigidFrictionalError> {
    finite(&hint)?;
    let length = norm(hint);
    if !length.is_finite() || length == 0.0 {
        return Err(invalid("finite nonzero tangent direction required"));
    }
    let cross_axis = cross(normal, hint.map(|x| x / length));
    let length = norm(cross_axis);
    if !length.is_finite() || length == 0.0 {
        return Err(invalid(
            "tangent direction must not be parallel to the normal",
        ));
    }
    let second = cross_axis.map(|x| x / length);
    let tangent = cross(second, normal);
    finite(&tangent)?;
    Ok([normal, tangent, second])
}
