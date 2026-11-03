//! One contact event across independently evolving native rods and finite bodies.
//!
//! Rod indices select original trajectories; node indices select their original
//! kinetic blocks. The graph core also serves the one-rod APIs, so introducing
//! another rod never substitutes a surrogate rigid body or another contact law.

pub mod friction;

use super::super::multiple::RodRigidMultiImpactError;
use super::super::{
    RodRigidImpactControls, RodRigidImpactError, RodRigidImpulseBalance, project_velocity,
};
use super::{
    ContactNetwork, RodRigidEndpoint, RodRigidGraphBudget, RodRigidGraphContact,
    RodRigidGraphTransfer, RodRigidParticipant, validate_graph_budget,
};
use crate::rod_motion::impulse::RodImpulseReport;
use crate::rod_motion::{RodMotion, RodMotionInertia, checkpoint};
use fs_contact::multicontact::{
    MultiImpactConfig, MultiImpactResult, NormalContactImpact, solve_normal_impacts,
};
use fs_exec::CancelGate;
use fs_mbd::{MomentumImpulseReceipt, RigidBodyState};

/// An attachment on one original rod or finite rigid body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RodRigidNetworkEndpoint {
    Rod {
        /// Index into the original rod slice, independent of the node index.
        rod: usize,
        node: usize,
        /// Offset from this rod's MATERIAL NODE in its material frame [m].
        point_body_m: [f64; 3],
    },
    Rigid {
        body: usize,
        /// Offset from this body's CENTER OF MASS in its body frame [m].
        arm_body_m: [f64; 3],
    },
}

/// Oriented physical witnesses, with +impulse on the positive endpoint.
#[derive(Clone, Copy, Debug)]
pub struct RodRigidNetworkContact {
    pub positive: RodRigidNetworkEndpoint,
    pub negative: RodRigidNetworkEndpoint,
    pub normal_world: [f64; 3],
    pub activation_tolerance_m: f64,
    pub restitution: f64,
}

impl From<RodRigidEndpoint> for RodRigidNetworkEndpoint {
    fn from(endpoint: RodRigidEndpoint) -> Self {
        match endpoint {
            RodRigidEndpoint::Rod { node, point_body_m } => Self::Rod {
                rod: 0,
                node,
                point_body_m,
            },
            RodRigidEndpoint::Rigid { body, arm_body_m } => Self::Rigid { body, arm_body_m },
        }
    }
}

/// The published one-rod witness names rod zero in an indexed network.
impl From<RodRigidGraphContact> for RodRigidNetworkContact {
    fn from(contact: RodRigidGraphContact) -> Self {
        Self {
            positive: contact.positive.into(),
            negative: contact.negative.into(),
            normal_world: contact.normal_world,
            activation_tolerance_m: contact.activation_tolerance_m,
            restitution: contact.restitution,
        }
    }
}

/// Bounds participant arrays and copied motion storage before graph allocation.
/// Counts include untouched participants because their continuing states are
/// returned. Rod and body counts have separate hard ceilings of 128. The node
/// bound counts the sum across all supplied rods, with checked arithmetic.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RodRigidNetworkBudget {
    pub maximum_rods: usize,
    pub maximum_bodies: usize,
    pub maximum_rod_nodes: usize,
}

impl Default for RodRigidNetworkBudget {
    fn default() -> Self {
        Self {
            maximum_rods: 128,
            maximum_bodies: 128,
            maximum_rod_nodes: 1_048_576,
        }
    }
}

/// Measured native publication of the complete normal or frictional network.
#[derive(Clone, Debug)]
pub struct RodRigidNetworkTransfer<const D: usize> {
    pub relative_velocities_after_m_s: Vec<[f64; D]>,
    pub endpoint_velocity_residuals_m_s: Vec<[f64; D]>,
    /// The same represented world vector acts oppositely at both endpoints.
    pub impulses_world_n_s: Vec<[f64; 3]>,
    /// One native batch receipt per supplied rod, in original input order.
    /// None when no nonzero point impulse was requested from that trajectory.
    pub rod_impulses: Vec<Option<RodImpulseReport>>,
    /// One native resultant receipt per body, in original input order.
    pub rigid_impulses: Vec<Option<MomentumImpulseReceipt>>,
    /// One compensated account across ALL original trajectories and bodies.
    /// None only when the original contact solution is entirely zero.
    pub balance: Option<RodRigidImpulseBalance>,
}

impl<const D: usize> RodRigidNetworkTransfer<D> {
    /// Called only by adapters that supplied exactly one original rod.
    pub(super) fn into_single_rod(self) -> RodRigidGraphTransfer<D> {
        debug_assert_eq!(self.rod_impulses.len(), 1);
        RodRigidGraphTransfer {
            relative_velocities_after_m_s: self.relative_velocities_after_m_s,
            endpoint_velocity_residuals_m_s: self.endpoint_velocity_residuals_m_s,
            impulses_world_n_s: self.impulses_world_n_s,
            rod_impulse: self.rod_impulses.into_iter().next().flatten(),
            rigid_impulses: self.rigid_impulses,
            balance: self.balance,
        }
    }
}

/// Original simultaneous normal law together with every native participant.
#[derive(Clone, Debug)]
pub struct RodRigidNetworkImpactReport {
    /// Redundant loads retain the original UniqueVelocity classification.
    pub contact: MultiImpactResult,
    pub delassus_kg_inverse: Vec<f64>,
    pub gaps_m: Vec<f64>,
    pub unit_normals_world: Vec<[f64; 3]>,
    pub relative_normal_speeds_before_m_s: Vec<f64>,
    pub transfer: RodRigidNetworkTransfer<1>,
}

/// Resolve one coupled normal event across independent native trajectories.
///
/// Edges can join nodes on different rods, distinct nodes on one rod, a rod and
/// a finite body, or two finite bodies. Both indices identify a rod kinetic
/// block: node zero on two rods is never the same participant. The rods share
/// an inertia representation type, but retain their own material models,
/// densities, tensors, supports, poses, integration histories and step counts.
/// Every rod referenced by a witness must name the same physical time. Unused
/// rods retain their state and do not constrain the event clock.
///
/// The original signed graph response and simultaneous contact law supply one
/// globally admitted solution. Each affected rod receives one original point
/// impulse batch; each affected finite body receives one original resultant.
/// Actual endpoint speeds and combined support-inclusive momentum and kinetic
/// work must pass before any candidate is returned. No pose/time/material
/// advancement occurs at this fixed-configuration event. Penetration, moving
/// supports and unresolved transfer retain their original refusals.
///
/// The original 64-row contact limit applies. Assembly uses O(contacts squared)
/// storage/work with O(contacts) original inverse-inertia actions. Staging uses
/// O(total rod nodes + bodies + contacts), with no clone per contact. Empty rod
/// and body slices are allowed when witnesses do not reference them; unchanged
/// participants preserve their bits and all-zero events invoke no native
/// impulse owner. Collision discovery, event timing, persistent friction and
/// thermal partition remain separate owners.
pub fn apply_rigid_normal_contact_network<I: RodMotionInertia>(
    rods: &[RodMotion<I>],
    bodies: &[RodRigidParticipant],
    witnesses: &[RodRigidNetworkContact],
    controls: RodRigidImpactControls,
    contact_config: MultiImpactConfig,
    budget: RodRigidNetworkBudget,
    gate: &CancelGate,
) -> Result<
    (
        Vec<RodMotion<I>>,
        Vec<RigidBodyState>,
        RodRigidNetworkImpactReport,
    ),
    RodRigidMultiImpactError,
> {
    checkpoint(gate)?;
    if witnesses.len() > 64 || witnesses.len() > contact_config.solver.max_variables {
        return Err(RodRigidMultiImpactError::Invalid(
            "normal network contact budget",
        ));
    }
    let network = ContactNetwork { rods };
    let geometry = network.prepare_contact_graph(bodies, witnesses, controls, budget, gate)?;
    let axes: Vec<_> = geometry
        .contacts
        .iter()
        .map(|point| [point.normal_world])
        .collect();
    let delassus = network.contact_graph_mobility(bodies, &geometry, &axes, gate)?;
    let contacts: Vec<_> = geometry
        .contacts
        .iter()
        .zip(witnesses)
        .map(|(point, witness)| {
            Ok(NormalContactImpact {
                gap: point.gap_m,
                activation_tolerance: witness.activation_tolerance_m,
                relative_normal_speed: project_velocity(
                    point.normal_world,
                    point.relative_velocity_world_m_s,
                )?,
                restitution: witness.restitution,
            })
        })
        .collect::<Result<_, RodRigidImpactError>>()?;
    let contact =
        solve_normal_impacts(&contacts, &delassus, contact_config, || gate.is_requested())?;
    let impulses: Vec<_> = contact.normal_impulses.iter().map(|&p| [p]).collect();
    let predicted: Vec<_> = contact
        .relative_normal_speeds_after
        .iter()
        .map(|&v| [v])
        .collect();
    let (next_rods, next_bodies, transfer) = network.stage_contact_graph(
        bodies,
        &geometry,
        &axes,
        &impulses,
        &predicted,
        contact.dissipated_energy,
        contact.relative_work,
        controls,
        gate,
    )?;
    let report = RodRigidNetworkImpactReport {
        contact,
        delassus_kg_inverse: delassus,
        gaps_m: geometry.contacts.iter().map(|point| point.gap_m).collect(),
        unit_normals_world: geometry
            .contacts
            .iter()
            .map(|point| point.normal_world)
            .collect(),
        relative_normal_speeds_before_m_s: contacts
            .iter()
            .map(|point| point.relative_normal_speed)
            .collect(),
        transfer,
    };
    checkpoint(gate)?;
    Ok((next_rods, next_bodies, report))
}

pub(super) fn validate_network_budget<I: RodMotionInertia>(
    rods: &[RodMotion<I>],
    body_count: usize,
    budget: RodRigidNetworkBudget,
) -> Result<(), RodRigidImpactError> {
    validate_graph_budget(
        body_count,
        RodRigidGraphBudget {
            maximum_bodies: budget.maximum_bodies,
        },
    )?;
    if budget.maximum_rods > 128 || rods.len() > budget.maximum_rods {
        return Err(RodRigidImpactError::Invalid("contact network rod budget"));
    }
    let mut nodes = 0_usize;
    for rod in rods {
        nodes = nodes
            .checked_add(rod.kinetics.len())
            .filter(|&count| count <= budget.maximum_rod_nodes)
            .ok_or(RodRigidImpactError::Invalid(
                "contact network rod-node budget",
            ))?;
    }
    Ok(())
}
