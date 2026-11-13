//! Simultaneous Coulomb impacts sharing one native rod and one finite body.
//!
//! The original non-associated contact owner solves the full normal/tangent
//! mobility. One aggregate native update per participant preserves both cross
//! response and support reactions without introducing an impact ordering.

use super::super::impulse::{RodImpulseReport, RodPointImpulse};
use super::super::{RodMotion, RodMotionInertia, checkpoint};
use super::friction::{
    RodRigidFrictionalContact, RodRigidFrictionalControls, RodRigidFrictionalError, contact_frame,
};
use super::{
    RigidImpulseAudit, RodRigidImpactError, RodRigidImpulseBalance, components, cross, cross_dd,
    dd, finite, norm, project_velocity, rotate, validate_controls, vector,
};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactResult, NormalVelocityLaw, solve_frictional_contacts,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_mbd::{MassProperties, MomentumImpulseReceipt, RigidBodyState};

/// One admitted event, with the original full matrix and constitutive witness.
#[derive(Clone, Debug)]
pub struct RodRigidFrictionalMultiImpactReport {
    /// Preserves UniqueImpulses versus UniqueVelocity and its total-work
    /// witness. Redundant impulse coordinates are a checked representative,
    /// not authoritative local loads, pressures or individual heat allocations.
    pub contact: FrictionalContactResult,
    /// Complete row-major 3N-by-3N mobility J M^-1 J^T [kg^-1], in original
    /// contact order, with normal/tangent-1/tangent-2 coordinates per contact.
    pub delassus_kg_inverse: Vec<f64>,
    pub gaps_m: Vec<f64>,
    /// Orthonormal right-handed world frames, normal first.
    pub frames_world: Vec<[[f64; 3]; 3]>,
    pub relative_velocities_before_m_s: Vec<[f64; 3]>,
    /// Reconstructed from both actually stored native endpoint states.
    pub relative_velocities_after_m_s: Vec<[f64; 3]>,
    /// Actual minus the original law's predicted endpoint velocities.
    pub endpoint_velocity_residuals_m_s: Vec<[f64; 3]>,
    /// World impulses on the rod. The finite body receives their opposites.
    pub impulses_world_n_s: Vec<[f64; 3]>,
    /// Absent for an all-zero solution; no native impulse owner runs then.
    pub rod_impulse: Option<RodImpulseReport>,
    pub rigid_impulse: Option<MomentumImpulseReceipt>,
    pub balance: Option<RodRigidImpulseBalance>,
}

struct DirectionResponse {
    rod_torque: [f64; 3],
    rod_rotation: [f64; 3],
    body_torque: [f64; 3],
    body_rotation: [f64; 3],
}

impl<I: RodMotionInertia> RodMotion<I> {
    /// Resolve all supplied normal and tangent contacts at one configuration.
    ///
    /// The complete 3N-by-3N mobility includes the original central full rod
    /// tensors, eccentric material attachments, admitted stationary supports
    /// and the shared body's original full tensor. Contacts on the same rod
    /// node share its translation and spin; every contact shares the finite
    /// body. The original Signorini--Coulomb owner resolves all directions
    /// together, including geometrically active initially separating contacts.
    ///
    /// The contact owner's determinacy, total-work and passivity admission is
    /// preserved. Accepted impulses are accumulated with compensated sums and
    /// applied through one original rod point-impulse call and one original
    /// rigid momentum-impulse call. Every actual endpoint triple and the total
    /// support-inclusive momentum and kinetic work must pass before either
    /// continuing state is returned. Tangential impulses require collocated
    /// original witnesses; no geometry or local heat partition is inferred.
    ///
    /// The original cap is 21 contacts (63 scalar rows), further bounded by
    /// the caller's contact and matrix budgets before allocation. Assembly is
    /// O(N^2), with O(N) original inverse-inertia actions. Native publication
    /// costs O(rod nodes + N), with one rod clone, no per-contact state copies,
    /// no time step and no new mass factorization or constitutive solver.
    pub fn apply_rigid_frictional_impacts(
        &self,
        rigid: RigidBodyState,
        properties: MassProperties,
        witnesses: &[RodRigidFrictionalContact],
        controls: RodRigidFrictionalControls,
        gate: &CancelGate,
    ) -> Result<(Self, RigidBodyState, RodRigidFrictionalMultiImpactReport), RodRigidFrictionalError>
    {
        checkpoint(gate)?;
        let count = witnesses.len();
        if count > 21
            || count > controls.contact.maximum_contacts
            || 3 * count > controls.contact.matrix.max_variables
        {
            return Err(RodRigidImpactError::Invalid("frictional contact count budget").into());
        }
        let limits = controls.impact;
        validate_controls(limits)?;
        let dimension = 3 * count;
        let mut points = Vec::with_capacity(count);
        let mut frames = Vec::with_capacity(count);
        let mut responses = Vec::with_capacity(dimension);
        let mut contacts = Vec::with_capacity(count);
        for (i, &witness) in witnesses.iter().enumerate() {
            checkpoint(gate)?;
            let geometry = witness.geometry;
            let point = self.prepare_rigid_point(rigid, properties, geometry, limits, gate)?;
            let frame = contact_frame(point.normal_world, witness.tangent_world)?;
            let node = point.node;
            let mut before = [0.0; 3];
            for a in 0..3 {
                let rod_torque = cross(
                    point.central_arm_body_m,
                    rotate(self.pose.frames[node].inverse(), frame[a])?,
                );
                finite(&rod_torque)?;
                let rod_rotation = if self.model.fixed_frames()[node] {
                    [0.0; 3]
                } else {
                    self.kinetics[node].omega(rod_torque)?
                };
                let body_direction = components(
                    rigid
                        .pose()
                        .orientation()
                        .rotate_world_to_body(vector(frame[a])),
                );
                let body_torque = cross(geometry.rigid_arm_body_m, body_direction);
                let body_rotation =
                    components(properties.angular_velocity_body_checked(vector(body_torque))?);
                responses.push(DirectionResponse {
                    rod_torque,
                    rod_rotation,
                    body_torque,
                    body_rotation,
                });
                before[a] = project_velocity(frame[a], point.relative_velocity_world_m_s)?;
            }
            contacts.push(FrictionalContact {
                id: i as u64,
                gap: point.gap_m,
                activation_tolerance: geometry.activation_tolerance_m,
                relative_velocity: before,
                friction_coefficient: witness.friction_coefficient,
                normal_law: NormalVelocityLaw::NewtonRestitution {
                    coefficient: geometry.restitution,
                },
            });
            points.push(point);
            frames.push(frame);
        }

        let mut delassus = vec![0.0; dimension * dimension];
        for row in 0..dimension {
            checkpoint(gate)?;
            let a = &responses[row];
            let first = &points[row / 3];
            let axis_a = frames[row / 3][row % 3];
            for col in row..dimension {
                let b = &responses[col];
                let second = &points[col / 3];
                let axis_b = frames[col / 3][col % 3];
                // The same native inverse tensor defines both directions.
                // Average the bilinear evaluations only to preserve symmetry
                // under roundoff, retaining every physical cross term.
                let mut value = dot_dd(axis_a, axis_b) / dd(properties.mass())
                    + dd(0.5)
                        * (dot_dd(a.body_torque, b.body_rotation)
                            + dot_dd(b.body_torque, a.body_rotation));
                if first.node == second.node {
                    let node = first.node;
                    for c in 0..3 {
                        if !self.model.fixed_positions()[node][c] {
                            value = value
                                + dd(axis_a[c]) * dd(axis_b[c]) / dd(self.kinetics[node].mass);
                        }
                    }
                    value = value
                        + dd(0.5)
                            * (dot_dd(a.rod_torque, b.rod_rotation)
                                + dot_dd(b.rod_torque, a.rod_rotation));
                }
                let value = value.to_f64();
                finite(&[value])?;
                delassus[row * dimension + col] = value;
                delassus[col * dimension + row] = value;
            }
        }
        let contact = solve_frictional_contacts(&contacts, &delassus, controls.contact, || {
            gate.is_requested()
        })?;
        let mut report = RodRigidFrictionalMultiImpactReport {
            gaps_m: contacts.iter().map(|c| c.gap).collect(),
            relative_velocities_before_m_s: contacts.iter().map(|c| c.relative_velocity).collect(),
            relative_velocities_after_m_s: contacts.iter().map(|c| c.relative_velocity).collect(),
            endpoint_velocity_residuals_m_s: vec![[0.0; 3]; count],
            impulses_world_n_s: vec![[0.0; 3]; count],
            contact,
            delassus_kg_inverse: delassus,
            frames_world: frames,
            rod_impulse: None,
            rigid_impulse: None,
            balance: None,
        };
        if report.contact.impulses.iter().flatten().all(|&p| p == 0.0) {
            checkpoint(gate)?;
            return Ok((self.clone(), rigid, report));
        }

        let mut rod_impulses = Vec::with_capacity(count);
        let mut body_linear = [Dd::ZERO; 3];
        let mut body_angular = [Dd::ZERO; 3];
        for i in 0..count {
            checkpoint(gate)?;
            let p = report.contact.impulses[i];
            if p.iter().all(|&x| x == 0.0) {
                continue;
            }
            if (p[1] != 0.0 || p[2] != 0.0)
                && norm(points[i].separation_world_m) > limits.geometry_tolerance_m
            {
                return Err(RodRigidImpactError::Invalid(
                    "tangential impulses require collocated original witnesses",
                )
                .into());
            }
            let axes = report.frames_world[i];
            let impulse = std::array::from_fn(|c| {
                (dd(axes[0][c]) * dd(p[0]) + dd(axes[1][c]) * dd(p[1]) + dd(axes[2][c]) * dd(p[2]))
                    .to_f64()
            });
            finite(&impulse)?;
            report.impulses_world_n_s[i] = impulse;
            // Request the body's couple from the SAME represented point
            // impulse given to the rod, then round the aggregate only once.
            let local = components(
                rigid
                    .pose()
                    .orientation()
                    .rotate_world_to_body(vector(impulse)),
            );
            let torque = cross_dd(witnesses[i].geometry.rigid_arm_body_m, local.map(dd));
            for c in 0..3 {
                body_linear[c] = body_linear[c] - dd(impulse[c]);
                body_angular[c] = body_angular[c] - torque[c];
            }
            rod_impulses.push(RodPointImpulse {
                node: witnesses[i].geometry.node,
                point_body_m: witnesses[i].geometry.point_body_m,
                linear_world_n_s: impulse,
                couple_world_nm_s: [0.0; 3],
            });
        }
        let body_linear = body_linear.map(Dd::to_f64);
        let body_angular = body_angular.map(Dd::to_f64);
        finite(&body_linear)?;
        finite(&body_angular)?;
        let (next_rod, rod_report) =
            self.apply_point_impulses(&rod_impulses, limits.impulse, gate)?;
        let rigid_report =
            rigid.apply_momentum_impulse(properties, vector(body_linear), vector(body_angular))?;
        let next_rigid = rigid_report.state_after;
        let mut relative_work = Dd::ZERO;
        let mut velocity_defect: f64 = 0.0;
        for i in 0..count {
            checkpoint(gate)?;
            let after_world = next_rod.relative_point_velocity(
                points[i].node,
                points[i].central_arm_body_m,
                next_rigid,
                properties,
                witnesses[i].geometry.rigid_arm_body_m,
            )?;
            for a in 0..3 {
                let after = project_velocity(report.frames_world[i][a], after_world)?;
                report.relative_velocities_after_m_s[i][a] = after;
                report.endpoint_velocity_residuals_m_s[i][a] =
                    after - report.contact.relative_velocities_after[i][a];
                relative_work = relative_work
                    + dd(report.contact.impulses[i][a])
                        * (dd(contacts[i].relative_velocity[a]) + dd(after))
                        * dd(0.5);
            }
            velocity_defect = velocity_defect.max(norm(report.endpoint_velocity_residuals_m_s[i]));
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
                dissipated_energy_j: report.contact.dissipated_energy,
                velocity_defect_m_s: velocity_defect,
            },
            limits,
            gate,
        )?;
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
        report.balance = Some(balance);
        checkpoint(gate)?;
        Ok((next_rod, next_rigid, report))
    }
}

fn dot_dd(a: [f64; 3], b: [f64; 3]) -> Dd {
    dd(a[0]) * dd(b[0]) + dd(a[1]) * dd(b[1]) + dd(a[2]) * dd(b[2])
}
