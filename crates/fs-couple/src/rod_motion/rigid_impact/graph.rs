//! Coupled impacts on a graph of native rod attachments and finite bodies.
//!
//! Rows retain their two signed physical endpoints. Shared original mass blocks
//! supply the cross response, so a body-body collision, a rod-body collision and
//! rod self-contact can participate in one contact law without a pair ordering.

pub mod friction;
pub mod network;

use self::network::{
    RodRigidNetworkBudget, RodRigidNetworkContact, RodRigidNetworkEndpoint,
    RodRigidNetworkTransfer, apply_rigid_normal_contact_network, validate_network_budget,
};
use super::super::impulse::{RodImpulseReport, RodPointImpulse};
use super::super::{RodMotion, RodMotionInertia, add, checkpoint, sub};
use super::multiple::RodRigidMultiImpactError;
use super::{
    RigidMomentumTransfer, RodMomentumTransfer, RodRigidImpactControls, RodRigidImpactError,
    RodRigidImpulseBalance, RodRigidSystemAudit, components, cross, cross_dd, dd, dot, finite,
    norm, project_velocity, rotate, validate_controls, vector,
};
use fs_contact::multicontact::{MultiImpactConfig, MultiImpactResult};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_mbd::{MassProperties, MomentumImpulseReceipt, RigidBodyState};

/// One original finite body. Its index is stable for this graph transaction.
#[derive(Clone, Copy, Debug)]
pub struct RodRigidParticipant {
    pub state: RigidBodyState,
    pub inertia: MassProperties,
}

/// A material attachment on the existing rod or an original finite body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RodRigidEndpoint {
    Rod {
        node: usize,
        /// Offset from the rod MATERIAL NODE in its material frame [m].
        point_body_m: [f64; 3],
    },
    Rigid {
        body: usize,
        /// Offset from the body's CENTER OF MASS in its body frame [m].
        arm_body_m: [f64; 3],
    },
}

/// Oriented physical witnesses. The positive endpoint receives +normal impulse.
#[derive(Clone, Copy, Debug)]
pub struct RodRigidGraphContact {
    pub positive: RodRigidEndpoint,
    pub negative: RodRigidEndpoint,
    /// Nonzero world direction, normalized without an overflow-prone square sum.
    pub normal_world: [f64; 3],
    pub activation_tolerance_m: f64,
    pub restitution: f64,
}

/// Bounds input/output body storage before any graph-sized allocation.
/// At most 128 bodies can occur in the original 64-row normal law. Zero is
/// admitted for a rod-only self-contact graph; unused bodies still count.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RodRigidGraphBudget {
    pub maximum_bodies: usize,
}
impl Default for RodRigidGraphBudget {
    fn default() -> Self {
        Self {
            maximum_bodies: 128,
        }
    }
}

/// Actual native publication common to normal and three-direction contact laws.
#[derive(Clone, Debug)]
pub struct RodRigidGraphTransfer<const D: usize> {
    /// Measured endpoints in each contact's declared directional frame.
    pub relative_velocities_after_m_s: Vec<[f64; D]>,
    /// Actual minus the original contact owner's predicted endpoint velocity.
    pub endpoint_velocity_residuals_m_s: Vec<[f64; D]>,
    /// One world impulse ON each positive endpoint; negative receives its opposite.
    pub impulses_world_n_s: Vec<[f64; 3]>,
    /// None when no nonzero rod point impulse was requested.
    pub rod_impulse: Option<RodImpulseReport>,
    /// In original body order. Unchanged resultants invoke no native update.
    pub rigid_impulses: Vec<Option<MomentumImpulseReceipt>>,
    /// None only for an all-zero contact solution.
    pub balance: Option<RodRigidImpulseBalance>,
}

/// Globally admitted normal law and the actual graph's continuing native motion.
#[derive(Clone, Debug)]
pub struct RodRigidGraphImpactReport {
    /// Preserves unique-impulse versus unique-velocity determinacy. Redundant
    /// coordinates remain one representative, not unique local pressure/heat.
    pub contact: MultiImpactResult,
    /// Complete row-major signed contact mobility J M^-1 J^T [kg^-1].
    pub delassus_kg_inverse: Vec<f64>,
    pub gaps_m: Vec<f64>,
    pub unit_normals_world: Vec<[f64; 3]>,
    pub relative_normal_speeds_before_m_s: Vec<f64>,
    pub transfer: RodRigidGraphTransfer<1>,
}

pub(super) struct GraphGeometry {
    pub(super) contacts: Vec<PreparedGraphContact>,
}

pub(super) struct PreparedGraphContact {
    pub(super) normal_world: [f64; 3],
    pub(super) gap_m: f64,
    pub(super) separation_world_m: [f64; 3],
    pub(super) relative_velocity_world_m_s: [Dd; 3],
    endpoints: [PreparedEndpoint; 2],
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Participant {
    Rod { rod: usize, node: usize },
    Rigid(usize),
}

#[derive(Clone, Copy)]
struct PreparedEndpoint {
    original: RodRigidNetworkEndpoint,
    owner: Participant,
    central_arm_body_m: [f64; 3],
    point_world_m: [f64; 3],
    translation_m_s: [Dd; 3],
    rotation_m_s: [f64; 3],
}

struct DirectionResponse {
    owner: Participant,
    direction_world: [f64; 3],
    torque_body: [f64; 3],
    rotation_body: [f64; 3],
}

/// The same endpoint, mobility and publication owner serves one or many rods.
/// Borrowing the original motions retains their model, material and time owners.
pub(super) struct ContactNetwork<'a, I: RodMotionInertia> {
    pub(super) rods: &'a [RodMotion<I>],
}

impl<I: RodMotionInertia> RodMotion<I> {
    /// Resolve a normal contact network in one original globally passive law.
    ///
    /// A contact can join a rod attachment to a finite body, two finite bodies,
    /// or distinct rod nodes. All shared original inertia blocks contribute to
    /// the complete signed response. The graph admits no contact between two
    /// attachments on the same rigid participant (including one rod node).
    ///
    /// Each affected participant receives one aggregate native update. All
    /// measured endpoint speeds, support-inclusive momentum and actual work
    /// must pass together before any continuing state is returned. Unaffected
    /// bodies retain their bits and order. The geometry/material state, clocks
    /// and original time-integration histories do not advance at this event.
    ///
    /// The original normal law bounds the solve to 64 rows. Assembly uses
    /// O(contacts^2) work/storage and O(contacts) native tensor actions; staging
    /// uses O(rod nodes + bodies + contacts) work without per-contact clones.
    /// Collision discovery, contact time, moving supports and heat destinations
    /// remain with their owners. Penetration is never repaired by this event.
    #[allow(clippy::too_many_arguments)]
    pub fn apply_rigid_normal_contact_graph(
        &self,
        bodies: &[RodRigidParticipant],
        witnesses: &[RodRigidGraphContact],
        controls: RodRigidImpactControls,
        contact_config: MultiImpactConfig,
        budget: RodRigidGraphBudget,
        gate: &CancelGate,
    ) -> Result<(Self, Vec<RigidBodyState>, RodRigidGraphImpactReport), RodRigidMultiImpactError>
    {
        checkpoint(gate)?;
        if witnesses.len() > 64 || witnesses.len() > contact_config.solver.max_variables {
            return Err(RodRigidMultiImpactError::Invalid(
                "normal graph contact budget",
            ));
        }
        validate_graph_budget(bodies.len(), budget)?;
        let contacts: Vec<_> = witnesses.iter().copied().map(Into::into).collect();
        let (mut next_rods, next_bodies, report) = apply_rigid_normal_contact_network(
            core::slice::from_ref(self),
            bodies,
            &contacts,
            controls,
            contact_config,
            RodRigidNetworkBudget {
                maximum_rods: 1,
                maximum_bodies: budget.maximum_bodies,
                maximum_rod_nodes: self.kinetics.len(),
            },
            gate,
        )?;
        let next_rod = next_rods
            .pop()
            .ok_or_else(|| invalid("absent native rod candidate"))?;
        Ok((
            next_rod,
            next_bodies,
            RodRigidGraphImpactReport {
                contact: report.contact,
                delassus_kg_inverse: report.delassus_kg_inverse,
                gaps_m: report.gaps_m,
                unit_normals_world: report.unit_normals_world,
                relative_normal_speeds_before_m_s: report.relative_normal_speeds_before_m_s,
                transfer: report.transfer.into_single_rod(),
            },
        ))
    }
}

impl<I: RodMotionInertia> ContactNetwork<'_, I> {
    pub(super) fn prepare_contact_graph(
        &self,
        bodies: &[RodRigidParticipant],
        witnesses: &[RodRigidNetworkContact],
        controls: RodRigidImpactControls,
        budget: RodRigidNetworkBudget,
        gate: &CancelGate,
    ) -> Result<GraphGeometry, RodRigidImpactError> {
        checkpoint(gate)?;
        validate_network_budget(self.rods, bodies.len(), budget)?;
        if witnesses.len() > 64 {
            return Err(invalid("contact graph row budget"));
        }
        validate_controls(controls)?;
        let mut contacts = Vec::with_capacity(witnesses.len());
        let mut event_time = None;
        for witness in witnesses {
            checkpoint(gate)?;
            let positive = self.graph_endpoint(bodies, witness.positive)?;
            let negative = self.graph_endpoint(bodies, witness.negative)?;
            // A contact network is one physical event. Unreferenced rods are
            // retained but do not make unrelated clocks constrain this event.
            for endpoint in [positive, negative] {
                if let Participant::Rod { rod, .. } = endpoint.owner {
                    let time = self.rods[rod].time_s;
                    if event_time.is_some_and(|previous| previous != time) {
                        return Err(invalid("contact network rod clocks differ"));
                    }
                    event_time = Some(time);
                }
            }
            if positive.owner == negative.owner {
                return Err(invalid(
                    "contact endpoints belong to the same rigid participant",
                ));
            }
            let normal = unit_direction(witness.normal_world)?;
            let separation = sub(positive.point_world_m, negative.point_world_m);
            finite(&separation)?;
            let gap = dot(normal, separation);
            finite(&[gap])?;
            let tangent = std::array::from_fn(|c| separation[c] - gap * normal[c]);
            if norm(tangent) > controls.geometry_tolerance_m {
                return Err(invalid(
                    "contact graph witnesses must share their normal line",
                ));
            }
            let relative_velocity = relative_velocity(positive, negative)?;
            contacts.push(PreparedGraphContact {
                normal_world: normal,
                gap_m: gap,
                separation_world_m: separation,
                relative_velocity_world_m_s: relative_velocity,
                endpoints: [positive, negative],
            });
        }
        Ok(GraphGeometry { contacts })
    }

    pub(super) fn contact_graph_mobility<const D: usize>(
        &self,
        bodies: &[RodRigidParticipant],
        geometry: &GraphGeometry,
        axes: &[[[f64; 3]; D]],
        gate: &CancelGate,
    ) -> Result<Vec<f64>, RodRigidImpactError> {
        checkpoint(gate)?;
        let count = geometry.contacts.len();
        if !matches!(D, 1 | 3) || axes.len() != count || count > 64 / D {
            return Err(invalid("contact graph directional matrix budget"));
        }
        let dimension = D * count;
        let mut responses = Vec::with_capacity(dimension);
        for (point, frame) in geometry.contacts.iter().zip(axes) {
            checkpoint(gate)?;
            for &axis in frame {
                finite(&axis)?;
                responses.push([
                    self.graph_direction(bodies, point.endpoints[0], axis)?,
                    self.graph_direction(bodies, point.endpoints[1], axis.map(|v| -v))?,
                ]);
            }
        }
        let mut matrix = vec![0.0; dimension * dimension];
        for i in 0..dimension {
            checkpoint(gate)?;
            for j in i..dimension {
                let mut value = Dd::ZERO;
                for a in &responses[i] {
                    for b in &responses[j] {
                        if a.owner != b.owner {
                            continue;
                        }
                        // Both actions use the same original symmetric tensor.
                        // Symmetric bilinear evaluation avoids contact-order
                        // dependence from the two finite inverse actions.
                        value = value
                            + dd(0.5)
                                * (dot_dd(a.torque_body, b.rotation_body)
                                    + dot_dd(b.torque_body, a.rotation_body));
                        match a.owner {
                            Participant::Rod { rod, node } => {
                                let rod = &self.rods[rod];
                                for c in 0..3 {
                                    if !rod.model.fixed_positions()[node][c] {
                                        value = value
                                            + dd(a.direction_world[c]) * dd(b.direction_world[c])
                                                / dd(rod.kinetics[node].mass);
                                    }
                                }
                            }
                            Participant::Rigid(body) => {
                                value = value
                                    + dot_dd(a.direction_world, b.direction_world)
                                        / dd(bodies[body].inertia.mass());
                            }
                        }
                    }
                }
                let value = value.to_f64();
                finite(&[value])?;
                matrix[i * dimension + j] = value;
                matrix[j * dimension + i] = value;
            }
        }
        Ok(matrix)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn stage_contact_graph<const D: usize>(
        &self,
        bodies: &[RodRigidParticipant],
        geometry: &GraphGeometry,
        axes: &[[[f64; 3]; D]],
        impulses: &[[f64; D]],
        predicted_after: &[[f64; D]],
        dissipated_energy_j: f64,
        law_relative_work_j: f64,
        controls: RodRigidImpactControls,
        gate: &CancelGate,
    ) -> Result<
        (
            Vec<RodMotion<I>>,
            Vec<RigidBodyState>,
            RodRigidNetworkTransfer<D>,
        ),
        RodRigidImpactError,
    > {
        checkpoint(gate)?;
        let count = geometry.contacts.len();
        if !matches!(D, 1 | 3)
            || axes.len() != count
            || impulses.len() != count
            || predicted_after.len() != count
        {
            return Err(invalid("contact graph impulse shape"));
        }
        let mut before = Vec::with_capacity(count);
        for (point, frame) in geometry.contacts.iter().zip(axes) {
            let mut velocity = [0.0; D];
            for a in 0..D {
                velocity[a] = project_velocity(frame[a], point.relative_velocity_world_m_s)?;
            }
            before.push(velocity);
        }
        let mut next_bodies: Vec<_> = bodies.iter().map(|body| body.state).collect();
        let mut transfer = RodRigidNetworkTransfer {
            relative_velocities_after_m_s: before.clone(),
            endpoint_velocity_residuals_m_s: vec![[0.0; D]; count],
            impulses_world_n_s: vec![[0.0; 3]; count],
            rod_impulses: vec![None; self.rods.len()],
            rigid_impulses: vec![None; bodies.len()],
            balance: None,
        };
        if impulses.iter().flatten().all(|&p| p == 0.0) {
            let mut rods = Vec::with_capacity(self.rods.len());
            for rod in self.rods {
                checkpoint(gate)?;
                rods.push(rod.clone());
            }
            checkpoint(gate)?;
            return Ok((rods, next_bodies, transfer));
        }
        // A batch belongs to a complete native trajectory, while each kinetic
        // block in the mobility is identified by BOTH its rod and node index.
        let mut rod_impulses = vec![Vec::new(); self.rods.len()];
        let mut body_linear = vec![[Dd::ZERO; 3]; bodies.len()];
        let mut body_angular = vec![[Dd::ZERO; 3]; bodies.len()];
        for i in 0..count {
            checkpoint(gate)?;
            finite(&impulses[i])?;
            finite(&predicted_after[i])?;
            if impulses[i].iter().all(|&p| p == 0.0) {
                continue;
            }
            if impulses[i][1..].iter().any(|&p| p != 0.0)
                && norm(geometry.contacts[i].separation_world_m) > controls.geometry_tolerance_m
            {
                return Err(invalid(
                    "tangential graph impulses require collocated witnesses",
                ));
            }
            let impulse = std::array::from_fn(|c| {
                (0..D)
                    .fold(Dd::ZERO, |sum, a| {
                        sum + dd(axes[i][a][c]) * dd(impulses[i][a])
                    })
                    .to_f64()
            });
            finite(&impulse)?;
            transfer.impulses_world_n_s[i] = impulse;
            for (side, endpoint) in geometry.contacts[i].endpoints.iter().enumerate() {
                // Both owners receive the same represented world vector.
                let applied = if side == 0 {
                    impulse
                } else {
                    impulse.map(|v| -v)
                };
                match endpoint.original {
                    RodRigidNetworkEndpoint::Rod {
                        rod,
                        node,
                        point_body_m,
                    } => {
                        rod_impulses[rod].push(RodPointImpulse {
                            node,
                            point_body_m,
                            linear_world_n_s: applied,
                            couple_world_nm_s: [0.0; 3],
                        });
                    }
                    RodRigidNetworkEndpoint::Rigid { body, arm_body_m } => {
                        let local = components(
                            bodies[body]
                                .state
                                .pose()
                                .orientation()
                                .rotate_world_to_body(vector(applied)),
                        );
                        let moment = cross_dd(arm_body_m, local.map(dd));
                        for c in 0..3 {
                            body_linear[body][c] = body_linear[body][c] + dd(applied[c]);
                            body_angular[body][c] = body_angular[body][c] + moment[c];
                        }
                    }
                }
            }
        }
        let mut next_rods = Vec::with_capacity(self.rods.len());
        for (index, (rod, impulses)) in self.rods.iter().zip(&rod_impulses).enumerate() {
            checkpoint(gate)?;
            let next = if impulses.is_empty() {
                rod.clone()
            } else {
                let (next, report) = rod.apply_point_impulses(impulses, controls.impulse, gate)?;
                transfer.rod_impulses[index] = Some(report);
                next
            };
            next_rods.push(next);
        }
        let mut body_transfers = Vec::with_capacity(bodies.len());
        for (index, body) in bodies.iter().enumerate() {
            checkpoint(gate)?;
            let linear = body_linear[index].map(Dd::to_f64);
            let angular = body_angular[index].map(Dd::to_f64);
            finite(&linear)?;
            finite(&angular)?;
            if linear == [0.0; 3] && angular == [0.0; 3] {
                continue;
            }
            let report =
                body.state
                    .apply_momentum_impulse(body.inertia, vector(linear), vector(angular))?;
            next_bodies[index] = report.state_after;
            body_transfers.push(RigidMomentumTransfer {
                before: body.state,
                after: report.state_after,
                inertia: body.inertia,
                expected_linear_world_n_s: body_linear[index],
                expected_angular_body_nm_s: body_angular[index],
            });
            transfer.rigid_impulses[index] = Some(report);
        }
        let candidates: Vec<_> = bodies
            .iter()
            .zip(&next_bodies)
            .map(|(body, &state)| RodRigidParticipant {
                state,
                inertia: body.inertia,
            })
            .collect();
        let mut relative_work = Dd::ZERO;
        let mut velocity_defect: f64 = 0.0;
        let candidate = ContactNetwork { rods: &next_rods };
        for i in 0..count {
            checkpoint(gate)?;
            let point = &geometry.contacts[i];
            let positive = candidate.graph_endpoint(&candidates, point.endpoints[0].original)?;
            let negative = candidate.graph_endpoint(&candidates, point.endpoints[1].original)?;
            let after = relative_velocity(positive, negative)?;
            for a in 0..D {
                let after = project_velocity(axes[i][a], after)?;
                let error = after - predicted_after[i][a];
                transfer.relative_velocities_after_m_s[i][a] = after;
                transfer.endpoint_velocity_residuals_m_s[i][a] = error;
                relative_work =
                    relative_work + dd(impulses[i][a]) * (dd(0.5 * before[i][a]) + dd(0.5 * after));
            }
            // At most three components; stable hypot does not overflow merely
            // because a representable residual's square is unrepresentable.
            let error = transfer.endpoint_velocity_residuals_m_s[i]
                .iter()
                .fold(0.0, |value, &e| fs_math::det::hypot(value, e));
            finite(&[error])?;
            velocity_defect = velocity_defect.max(error);
        }
        let rod_transfers: Vec<_> = self
            .rods
            .iter()
            .zip(&next_rods)
            .zip(&transfer.rod_impulses)
            .map(|((before, after), receipt)| RodMomentumTransfer {
                before,
                after,
                receipt: receipt.as_ref(),
            })
            .collect();
        let balance = RodMotion::<I>::audit_rod_rigid_network(
            &rod_transfers,
            RodRigidSystemAudit {
                bodies: &body_transfers,
                relative_work_j: relative_work,
                dissipated_energy_j,
                velocity_defect_m_s: velocity_defect,
            },
            controls,
            gate,
        )?;
        let law_work_defect = (relative_work - dd(law_relative_work_j)).to_f64().abs();
        finite(&[law_work_defect])?;
        if law_work_defect > controls.impulse.work_tolerance_j {
            return Err(RodRigidImpactError::Balance {
                velocity_defect_m_s: velocity_defect,
                linear_defect_n_s: norm(balance.linear_momentum_residual_n_s),
                angular_defect_nm_s: norm(balance.angular_momentum_residual_nm_s),
                work_defect_j: law_work_defect,
            });
        }
        transfer.balance = Some(balance);
        checkpoint(gate)?;
        Ok((next_rods, next_bodies, transfer))
    }

    fn graph_endpoint(
        &self,
        bodies: &[RodRigidParticipant],
        original: RodRigidNetworkEndpoint,
    ) -> Result<PreparedEndpoint, RodRigidImpactError> {
        let (owner, arm, point, translation, rotation) = match original {
            RodRigidNetworkEndpoint::Rod {
                rod,
                node,
                point_body_m,
            } => {
                finite(&point_body_m)?;
                let motion = self
                    .rods
                    .get(rod)
                    .ok_or_else(|| invalid("absent graph rod"))?;
                let k = *motion
                    .kinetics
                    .get(node)
                    .ok_or_else(|| invalid("absent graph rod node"))?;
                let frame = motion.pose.frames[node];
                let arm = sub(point_body_m, k.center);
                let point = add(motion.pose.positions_m[node], rotate(frame, point_body_m)?);
                let rotation = rotate(frame, cross(k.omega(motion.spin[node])?, arm))?;
                let translation = motion.linear[node].map(|p| dd(p) / dd(k.mass));
                (
                    Participant::Rod { rod, node },
                    arm,
                    point,
                    translation,
                    rotation,
                )
            }
            RodRigidNetworkEndpoint::Rigid { body, arm_body_m } => {
                finite(&arm_body_m)?;
                let body_value = *bodies
                    .get(body)
                    .ok_or_else(|| invalid("absent graph rigid body"))?;
                let kinematics = body_value
                    .state
                    .point_kinematics(body_value.inertia, vector(arm_body_m))?;
                let translation = components(body_value.state.linear_momentum_world())
                    .map(|p| dd(p) / dd(body_value.inertia.mass()));
                let rotation =
                    components(body_value.state.pose().orientation().rotate_body_to_world(
                        kinematics.angular_velocity_body.cross(vector(arm_body_m)),
                    ));
                (
                    Participant::Rigid(body),
                    arm_body_m,
                    components(kinematics.point_world),
                    translation,
                    rotation,
                )
            }
        };
        finite(&arm)?;
        finite(&point)?;
        finite(&translation.map(Dd::to_f64))?;
        finite(&rotation)?;
        Ok(PreparedEndpoint {
            original,
            owner,
            central_arm_body_m: arm,
            point_world_m: point,
            translation_m_s: translation,
            rotation_m_s: rotation,
        })
    }

    fn graph_direction(
        &self,
        bodies: &[RodRigidParticipant],
        endpoint: PreparedEndpoint,
        direction: [f64; 3],
    ) -> Result<DirectionResponse, RodRigidImpactError> {
        let (torque, rotation) = match endpoint.owner {
            Participant::Rod { rod, node } => {
                let rod = &self.rods[rod];
                let local = rotate(rod.pose.frames[node].inverse(), direction)?;
                let torque = cross(endpoint.central_arm_body_m, local);
                finite(&torque)?;
                let rotation = if rod.model.fixed_frames()[node] {
                    [0.0; 3]
                } else {
                    rod.kinetics[node].omega(torque)?
                };
                (torque, rotation)
            }
            Participant::Rigid(body) => {
                let body = bodies[body];
                let local = components(
                    body.state
                        .pose()
                        .orientation()
                        .rotate_world_to_body(vector(direction)),
                );
                let torque = cross(endpoint.central_arm_body_m, local);
                let rotation =
                    components(body.inertia.angular_velocity_body_checked(vector(torque))?);
                (torque, rotation)
            }
        };
        finite(&torque)?;
        finite(&rotation)?;
        Ok(DirectionResponse {
            owner: endpoint.owner,
            direction_world: direction,
            torque_body: torque,
            rotation_body: rotation,
        })
    }
}

fn invalid(message: &'static str) -> RodRigidImpactError {
    RodRigidImpactError::Invalid(message)
}

pub(super) fn validate_graph_budget(
    body_count: usize,
    budget: RodRigidGraphBudget,
) -> Result<(), RodRigidImpactError> {
    if budget.maximum_bodies > 128 || body_count > budget.maximum_bodies {
        return Err(invalid("contact graph body budget"));
    }
    Ok(())
}

fn relative_velocity(
    positive: PreparedEndpoint,
    negative: PreparedEndpoint,
) -> Result<[Dd; 3], RodRigidImpactError> {
    // Common translation cancels before either rotational contribution enters.
    let value = std::array::from_fn(|c| {
        positive.translation_m_s[c] - negative.translation_m_s[c] + dd(positive.rotation_m_s[c])
            - dd(negative.rotation_m_s[c])
    });
    finite(&value.map(Dd::to_f64))?;
    Ok(value)
}

fn dot_dd(a: [f64; 3], b: [f64; 3]) -> Dd {
    dd(a[0]) * dd(b[0]) + dd(a[1]) * dd(b[1]) + dd(a[2]) * dd(b[2])
}

fn unit_direction(direction: [f64; 3]) -> Result<[f64; 3], RodRigidImpactError> {
    finite(&direction)?;
    let scale = direction.into_iter().map(f64::abs).fold(0.0, f64::max);
    if scale == 0.0 {
        return Err(invalid("nonzero contact graph normal required"));
    }
    let scaled = direction.map(|v| v / scale);
    let length = norm(scaled);
    let normalized = scaled.map(|v| v / length);
    finite(&normalized)?;
    Ok(normalized)
}
