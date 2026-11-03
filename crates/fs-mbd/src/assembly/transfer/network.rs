//! Simultaneous closed material handoffs through a finite-body network.
//!
//! Every outgoing parcel is released from its donor's ORIGINAL snapshot.
//! Only after all releases do receivers jointly capture their incoming parcels.
//! A cycle therefore never re-exports newly received material during this event.

use super::{PointMassTransfer, PointMassTransferControls, attachment_mismatch, combined_balance};
use crate::assembly::capture::{PointMassMomentumCapture, capture_point_masses_with_momenta};
use crate::assembly::release::{ReleasedPointMass, release_point_masses};
use crate::assembly::{AssemblyBalance, AssemblyError, Sum, finite, maximum, poll, vector};
use crate::{MassProperties, RigidBodyState, Vec3};

/// One directed handoff between distinct original body indices.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointMassTransferEdge {
    /// Original donor index.
    pub donor: usize,
    /// Original receiving index; must differ from the donor.
    pub receiver: usize,
    /// Mass, donor-frame impulse and the two original material-point arms.
    pub transfer: PointMassTransfer,
}

/// One endpoint in the original body order, including untouched bodies.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointMassTransferEndpoint {
    /// Actual native state after all outgoing release and incoming capture.
    pub state: RigidBodyState,
    /// Final native mass and full central inertia, in original body axes.
    pub properties: MassProperties,
    /// Combined COM shift in original body axes [m]; rebase retained points once.
    pub center_shift_body_m: Vec3,
    /// Work supplied by this body's simultaneous separation impulses [J].
    pub supplied_work_j: f64,
    /// Nonnegative loss of this body's joint incoming sticking event [J].
    pub removed_relative_energy_j: f64,
    /// Aggregate capture impulse on the retained post-release material [N s].
    pub capture_impulse_world_ns: Vec3,
    /// Native release balance; exactly zero for an untouched outgoing row.
    pub release_balance: AssemblyBalance,
    /// Native capture balance; exactly zero for an untouched incoming row.
    pub capture_balance: AssemblyBalance,
}

/// One atomic closed network event, with no additional free parcel inventory.
#[derive(Clone, Debug, PartialEq)]
pub struct PointMassTransferNetworkResult {
    /// Final states and per-body source/loss accounts in original body order.
    pub bodies: Vec<PointMassTransferEndpoint>,
    /// Sum of all jointly computed separation work [J]; source debits are external.
    pub supplied_work_j: f64,
    /// Sum of all receiver sticking losses [J]; thermal destinations are external.
    pub removed_relative_energy_j: f64,
    /// Maximum absolute attachment mismatch component actually checked [m].
    pub maximum_attachment_mismatch_m: f64,
    /// Actual final minus original total mass [kg].
    pub mass_residual_kg: f64,
    /// Common world origin for the actual angular residual [m].
    pub balance_origin_world_m: Vec3,
    /// Actual global momentum and kinetic/work/loss residuals.
    pub balance: AssemblyBalance,
}

/// Transfer point material through a closed directed graph in one physical event.
///
/// Every edge uses ORIGINAL donor/receiver body frames. Each donor releases all
/// its outgoing parcels simultaneously, retaining positive mass and physical
/// inertia. Each receiver then captures all incoming canonical momenta together,
/// with arms rebased by its own release COM shift. This two-phase definition
/// includes cycles and bodies that both give and receive, without edge-order
/// feedback or repeated per-parcel factorization.
///
/// Physical attachment coincidence is checked in local geometry and again at
/// the actual post-release receiver and parcel positions. Instantaneous handoff
/// cannot cross a spatial gap. Original material membership, species/composition,
/// contact discovery, finite flight, work-source debit and heat remain caller-owned.
///
/// The shared part cap counts all original bodies plus EVERY edge, even zeros.
/// Every positive capture and nonzero separation impulse consumes an impulse
/// slot. Two prefix-indexed adjacency tables give O(bodies + edges) routing and
/// bounded O(bodies + edges) storage. Incoming/outgoing descriptor buffers are
/// reused; each participating body admits at most one release and one capture
/// inertia. Fixed body and input order give deterministic native reductions.
///
/// Refusal, an invalid intermediate retained body, failed local/global balances,
/// or cancellation publishes no endpoints and leaves every borrowed body intact.
pub fn transfer_point_mass_network(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    edges: &[PointMassTransferEdge],
    controls: PointMassTransferControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PointMassTransferNetworkResult, AssemblyError> {
    poll(&mut cancelled)?;
    if bodies.is_empty() || bodies.len() != properties.len() {
        return Err(AssemblyError::Invalid(
            "matching nonempty transfer body/property axes required",
        ));
    }
    let parts = bodies
        .len()
        .checked_add(edges.len())
        .ok_or(AssemblyError::Invalid(
            "network transfer part count overflow",
        ))?;
    controls.assembly.validate(parts, 0)?;
    if [controls.attachment_tolerance_m, controls.mass_tolerance_kg]
        .into_iter()
        .any(|v| !v.is_finite() || v < 0.0)
    {
        return Err(AssemblyError::Invalid(
            "finite nonnegative network transfer tolerances required",
        ));
    }
    let mut endpoints = Vec::with_capacity(bodies.len());
    for (&body, &property) in bodies.iter().zip(properties) {
        poll(&mut cancelled)?;
        property.validate()?;
        body.center_of_mass_velocity_world(property)?;
        property.angular_velocity_body_checked(body.angular_momentum_body())?;
        endpoints.push(PointMassTransferEndpoint {
            state: body,
            properties: property,
            center_shift_body_m: Vec3::ZERO,
            supplied_work_j: 0.0,
            removed_relative_energy_j: 0.0,
            capture_impulse_world_ns: Vec3::ZERO,
            release_balance: zero_balance(),
            capture_balance: zero_balance(),
        });
    }
    let mut impulses = 0_usize;
    let mut active = 0;
    let mut anchor = bodies.len();
    let mut mismatch = 0.0_f64;
    for edge in edges {
        poll(&mut cancelled)?;
        if edge.donor >= bodies.len()
            || edge.receiver >= bodies.len()
            || edge.donor == edge.receiver
        {
            return Err(AssemblyError::Invalid(
                "transfer requires distinct existing body indices",
            ));
        }
        let transfer = edge.transfer;
        let parcel = transfer.release;
        if !parcel.mass_kg.is_finite() || parcel.mass_kg < 0.0 {
            return Err(AssemblyError::Invalid(
                "finite nonnegative transferred mass required",
            ));
        }
        vector(parcel.attachment_arm_body_m)?;
        vector(parcel.impulse_on_parcel_body_ns)?;
        vector(transfer.receiver_attachment_arm_body_m)?;
        if parcel.mass_kg == 0.0 {
            if parcel.impulse_on_parcel_body_ns != Vec3::ZERO {
                return Err(AssemblyError::Invalid(
                    "zero transferred mass requires zero impulse",
                ));
            }
            continue;
        }
        active += 1;
        anchor = anchor.min(edge.donor).min(edge.receiver);
        impulses = impulses
            .checked_add(1 + usize::from(parcel.impulse_on_parcel_body_ns != Vec3::ZERO))
            .ok_or(AssemblyError::Invalid(
                "network transfer impulse count overflow",
            ))?;
        mismatch = mismatch.max(attachment_mismatch(
            &bodies[edge.donor],
            &bodies[edge.receiver],
            &transfer,
        )?);
    }
    controls.assembly.validate(parts, impulses)?;
    if mismatch > controls.attachment_tolerance_m {
        return Err(AssemblyError::Invalid(
            "noncoincident network transfer attachments",
        ));
    }
    if active == 0 {
        poll(&mut cancelled)?;
        return Ok(PointMassTransferNetworkResult {
            bodies: endpoints,
            supplied_work_j: 0.0,
            removed_relative_energy_j: 0.0,
            maximum_attachment_mismatch_m: 0.0,
            mass_residual_kg: 0.0,
            balance_origin_world_m: bodies[0].pose().position_world(),
            balance: zero_balance(),
        });
    }

    let outgoing = Incidence::new(bodies.len(), edges, |edge| edge.donor, &mut cancelled)?;
    let incoming = Incidence::new(bodies.len(), edges, |edge| edge.receiver, &mut cancelled)?;
    let mut parcels: Vec<Option<ReleasedPointMass>> = vec![None; edges.len()];
    let mut releases = Vec::with_capacity(outgoing.maximum_degree);
    for (index, endpoint) in endpoints.iter_mut().enumerate() {
        poll(&mut cancelled)?;
        let row = outgoing.row(index);
        if row.is_empty() {
            continue;
        }
        releases.clear();
        for &edge in row {
            poll(&mut cancelled)?;
            releases.push(edges[edge].transfer.release);
        }
        let released = release_point_masses(
            &bodies[index],
            properties[index],
            &releases,
            controls.assembly,
            &mut cancelled,
        )?;
        if released.parcels.len() != row.len() {
            return Err(AssemblyError::Invalid(
                "native network release axis mismatch",
            ));
        }
        for (&edge, parcel) in row.iter().zip(released.parcels) {
            poll(&mut cancelled)?;
            parcels[edge] = Some(parcel);
        }
        endpoint.state = released.state;
        endpoint.properties = released.properties;
        endpoint.center_shift_body_m = released.center_shift_body_m;
        endpoint.supplied_work_j = released.supplied_work_j;
        endpoint.release_balance = released.balance;
    }
    let mut captures = Vec::with_capacity(incoming.maximum_degree);
    for (index, endpoint) in endpoints.iter_mut().enumerate() {
        poll(&mut cancelled)?;
        let row = incoming.row(index);
        if row.is_empty() {
            continue;
        }
        captures.clear();
        for &edge in row {
            poll(&mut cancelled)?;
            let parcel = parcels[edge].take().ok_or(AssemblyError::Invalid(
                "missing or repeated network parcel handoff",
            ))?;
            let arm = vector(
                edges[edge]
                    .transfer
                    .receiver_attachment_arm_body_m
                    .sub(endpoint.center_shift_body_m),
            )?;
            let point = endpoint.state.pose().point_world_from_body(arm)?;
            mismatch = mismatch.max(maximum(vector(point.sub(parcel.position_world_m))?));
            if mismatch > controls.attachment_tolerance_m {
                return Err(AssemblyError::Invalid(
                    "unrepresented network attachment coincidence",
                ));
            }
            captures.push(PointMassMomentumCapture {
                mass_kg: parcel.mass_kg,
                attachment_arm_body_m: arm,
                linear_momentum_world_ns: parcel.linear_momentum_world_ns,
            });
        }
        let captured = capture_point_masses_with_momenta(
            &endpoint.state,
            endpoint.properties,
            &captures,
            controls.assembly,
            &mut cancelled,
        )?;
        endpoint.state = captured.state;
        endpoint.properties = captured.properties;
        endpoint.center_shift_body_m = vector(
            endpoint
                .center_shift_body_m
                .add(captured.center_shift_body_m),
        )?;
        endpoint.removed_relative_energy_j = captured.removed_relative_energy_j;
        endpoint.capture_impulse_world_ns = captured.captured_impulse_world_ns;
        endpoint.capture_balance = captured.balance;
    }
    let mut work = Sum::default();
    let mut loss = Sum::default();
    for endpoint in &endpoints {
        poll(&mut cancelled)?;
        work.add(endpoint.supplied_work_j);
        loss.add(endpoint.removed_relative_energy_j);
    }
    let work = finite(work.value())?;
    let loss = finite(loss.value())?;
    let origin = bodies[anchor].pose().position_world();
    let boost = bodies[anchor].center_of_mass_velocity_world(properties[anchor])?;
    let (mass_residual, balance) = combined_balance(
        bodies.iter().copied().zip(properties.iter().copied()),
        endpoints
            .iter()
            .map(|endpoint| (endpoint.state, endpoint.properties)),
        origin,
        boost,
        work,
        loss,
        &mut cancelled,
    )?;
    if mass_residual.abs() > controls.mass_tolerance_kg {
        return Err(AssemblyError::Invalid(
            "combined network mass balance exceeds tolerance",
        ));
    }
    controls.assembly.admit(balance)?;
    poll(&mut cancelled)?;
    Ok(PointMassTransferNetworkResult {
        bodies: endpoints,
        supplied_work_j: work,
        removed_relative_energy_j: loss,
        maximum_attachment_mismatch_m: mismatch,
        mass_residual_kg: mass_residual,
        balance_origin_world_m: origin,
        balance,
    })
}

struct Incidence {
    offsets: Vec<usize>,
    indices: Vec<usize>,
    maximum_degree: usize,
}
impl Incidence {
    fn new(
        bodies: usize,
        edges: &[PointMassTransferEdge],
        owner: impl Fn(&PointMassTransferEdge) -> usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AssemblyError> {
        let length = bodies
            .checked_add(1)
            .ok_or(AssemblyError::Invalid("transfer adjacency size overflow"))?;
        let mut offsets = vec![0_usize; length];
        for edge in edges {
            poll(cancelled)?;
            if edge.transfer.release.mass_kg > 0.0 {
                offsets[owner(edge) + 1] += 1;
            }
        }
        let mut maximum_degree = 0;
        for index in 0..bodies {
            poll(cancelled)?;
            maximum_degree = maximum_degree.max(offsets[index + 1]);
            offsets[index + 1] += offsets[index];
        }
        let mut indices = vec![0; offsets[bodies]];
        let mut cursor = offsets[..bodies].to_vec();
        for (index, edge) in edges.iter().enumerate() {
            poll(cancelled)?;
            if edge.transfer.release.mass_kg == 0.0 {
                continue;
            }
            let owner = owner(edge);
            indices[cursor[owner]] = index;
            cursor[owner] += 1;
        }
        Ok(Self {
            offsets,
            indices,
            maximum_degree,
        })
    }
    fn row(&self, body: usize) -> &[usize] {
        &self.indices[self.offsets[body]..self.offsets[body + 1]]
    }
}
const fn zero_balance() -> AssemblyBalance {
    AssemblyBalance {
        linear_momentum_residual_world_ns: Vec3::ZERO,
        angular_momentum_residual_world_nms: Vec3::ZERO,
        energy_residual_j: 0.0,
    }
}
