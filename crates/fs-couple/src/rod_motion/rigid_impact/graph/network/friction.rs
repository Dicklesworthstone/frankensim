//! One coupled Coulomb event across independent native rods and finite bodies.
//!
//! The network owner supplies geometry, signed original mass response and all
//! native momentum staging. This adapter selects the original three-axis law.

use super::super::super::friction::{
    RodRigidFrictionalControls, RodRigidFrictionalError, contact_frame,
};
use super::super::super::{RodRigidImpactError, project_velocity};
use super::super::{ContactNetwork, RodRigidParticipant};
use super::{
    RodRigidNetworkBudget, RodRigidNetworkContact, RodRigidNetworkTransfer, validate_network_budget,
};
use crate::rod_motion::{RodMotion, RodMotionInertia, checkpoint};
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactResult, NormalVelocityLaw, solve_frictional_contacts,
};
use fs_exec::CancelGate;
use fs_mbd::RigidBodyState;

/// An oriented native network edge and its resolved event-local friction law.
#[derive(Clone, Copy, Debug)]
pub struct RodRigidFrictionalNetworkContact {
    pub geometry: RodRigidNetworkContact,
    /// World direction projected into the admitted contact plane.
    pub tangent_world: [f64; 3],
    /// Finite, nonnegative Coulomb coefficient; no material default is inferred.
    pub friction_coefficient: f64,
}

/// The original coupled law and the actual multi-participant native response.
/// Redundant impulse representatives do not determine unique local loads or
/// thermal allocations; the original determinacy and work witnesses are kept.
#[derive(Clone, Debug)]
pub struct RodRigidFrictionalNetworkReport {
    pub contact: FrictionalContactResult,
    /// Complete contact-major 3N-by-3N mobility [kg^-1], normal first.
    pub delassus_kg_inverse: Vec<f64>,
    pub gaps_m: Vec<f64>,
    pub frames_world: Vec<[[f64; 3]; 3]>,
    pub relative_velocities_before_m_s: Vec<[f64; 3]>,
    /// Per-rod and per-body native receipts, measured velocities and total work.
    pub transfer: RodRigidNetworkTransfer<3>,
}

/// Resolve simultaneous normal and tangent impulses across independent rods.
///
/// Edges may join different rods, distinct nodes of one rod, a rod and a finite
/// body, or two finite bodies. Shared participants supply their original signed
/// cross response, including full tensors, eccentric material attachments and
/// supported coordinates. Equal node indices on different rods remain distinct.
/// Participating rods must represent the same physical event time; unused rods
/// retain their own histories and impose no event-clock constraint.
///
/// The original non-associated Signorini--Coulomb law admits the whole network.
/// Common staging applies one original impulse batch per affected rod and one
/// resultant update per affected body. Actual endpoint triples, each original
/// transfer, support reactions, aggregate momentum and kinetic work must pass
/// together before any continuing state is returned. Tangential impulses need
/// collocated witnesses. Poses, constitutive histories and clocks do not advance.
///
/// At most 21 contacts / 63 directional rows enter the original law. Participant
/// and total-node budgets are checked before network arrays are allocated.
/// Assembly is quadratic in contacts with a linear number of tensor actions;
/// native staging is linear in total rod nodes, bodies and contacts, without
/// per-edge trajectory clones. Empty participant slices are admitted when no
/// witness refers to them. Zero-impulse events invoke no native impulse owner.
/// Collision discovery, persistent friction and thermal partition stay with
/// their existing owners; failed or unresolved events return no partial state.
pub fn apply_rigid_frictional_contact_network<I: RodMotionInertia>(
    rods: &[RodMotion<I>],
    bodies: &[RodRigidParticipant],
    witnesses: &[RodRigidFrictionalNetworkContact],
    controls: RodRigidFrictionalControls,
    budget: RodRigidNetworkBudget,
    gate: &CancelGate,
) -> Result<
    (
        Vec<RodMotion<I>>,
        Vec<RigidBodyState>,
        RodRigidFrictionalNetworkReport,
    ),
    RodRigidFrictionalError,
> {
    checkpoint(gate)?;
    let count = witnesses.len();
    if count > 21
        || count > controls.contact.maximum_contacts
        || 3 * count > controls.contact.matrix.max_variables
    {
        return Err(RodRigidImpactError::Invalid("frictional network contact budget").into());
    }
    validate_network_budget(rods, bodies.len(), budget)?;
    let network = ContactNetwork { rods };
    let edges: Vec<_> = witnesses.iter().map(|witness| witness.geometry).collect();
    let geometry = network.prepare_contact_graph(bodies, &edges, controls.impact, budget, gate)?;
    let mut frames = Vec::with_capacity(count);
    let mut contacts = Vec::with_capacity(count);
    for (i, witness) in witnesses.iter().enumerate() {
        checkpoint(gate)?;
        let point = &geometry.contacts[i];
        let frame = contact_frame(point.normal_world, witness.tangent_world)?;
        let mut before = [0.0; 3];
        for a in 0..3 {
            before[a] = project_velocity(frame[a], point.relative_velocity_world_m_s)?;
        }
        frames.push(frame);
        contacts.push(FrictionalContact {
            id: i as u64,
            gap: point.gap_m,
            activation_tolerance: witness.geometry.activation_tolerance_m,
            relative_velocity: before,
            friction_coefficient: witness.friction_coefficient,
            normal_law: NormalVelocityLaw::NewtonRestitution {
                coefficient: witness.geometry.restitution,
            },
        });
    }
    let delassus = network.contact_graph_mobility(bodies, &geometry, &frames, gate)?;
    let contact = solve_frictional_contacts(&contacts, &delassus, controls.contact, || {
        gate.is_requested()
    })?;
    let (rods, states, transfer) = network.stage_contact_graph(
        bodies,
        &geometry,
        &frames,
        &contact.impulses,
        &contact.relative_velocities_after,
        contact.dissipated_energy,
        contact.relative_work,
        controls.impact,
        gate,
    )?;
    let report = RodRigidFrictionalNetworkReport {
        contact,
        delassus_kg_inverse: delassus,
        gaps_m: contacts.iter().map(|point| point.gap).collect(),
        frames_world: frames,
        relative_velocities_before_m_s: contacts
            .iter()
            .map(|point| point.relative_velocity)
            .collect(),
        transfer,
    };
    checkpoint(gate)?;
    Ok((rods, states, report))
}
