//! Coupled normal and two-tangent Coulomb impulses on articulated trees.
//!
//! This adapter reuses the normal adapter's link Jacobians, original ABA
//! inverse-mass actions, accumulated impulse application, and physical endpoint
//! audit. The existing contact owner supplies the complete non-associated
//! Signorini--Coulomb solve; no independent link masses or sequential tangent
//! clipping are substituted for articulated recoil.

/// Simultaneous frictional impacts with rigid bilateral distance reactions.
pub mod constrained;

/// Directional friction with the same native articulated contact mechanics.
pub mod elliptic;

use super::{
    ArticulatedContactError, ArticulatedContactSystem, ContactCounterpart, EventResponse,
    EventWorkspace, LinkContactPoint, Prepared, finite, poll, prepare_point, unit_normal,
};
use crate::articulated::FreeFloatingImpulseResponse;
use crate::articulated::contact_geometry::prepare_contact_frame;
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactDeterminacy,
    FrictionalContactResult, NormalVelocityLaw, solve_frictional_contacts,
};
use fs_ga::Vec3;

/// One supplied geometric witness with an explicit ideal Coulomb law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArticulatedFrictionalContact<F = f64> {
    /// Stable distinct identity used to order assembly and impulse accumulation.
    pub id: u64,
    /// Original link point receiving the positive impulse.
    pub point: LinkContactPoint,
    /// Opposite link point or prescribed translating plane.
    pub other: ContactCounterpart,
    /// Nonzero outward normal from `other` toward `point`, normalized here.
    pub normal_world: Vec3,
    /// Nonparallel tangent hint; its normal component is discarded.
    pub tangent_world: Vec3,
    /// Maximum admitted nonnegative normal gap, in metres.
    pub activation_tolerance: f64,
    /// Constitutive payload of the selected entry point. The default scalar is
    /// one nonnegative static/kinetic coefficient; the directional entry point
    /// carries its principal pair and axis. No velocity dependence is inferred.
    pub friction_coefficient: F,
    /// Explicit normal endpoint law owned by `fs-contact`.
    pub normal_law: NormalVelocityLaw,
}

/// All native system responses and independently reconstructed balances.
#[derive(Debug, Clone)]
pub struct ArticulatedFrictionalContactResult<R = FrictionalContactResult> {
    /// Canonical contact order for every contact-space vector and matrix below.
    pub contact_ids: Vec<u64>,
    /// Updated states and physical momenta, in original system order.
    pub responses: Vec<FreeFloatingImpulseResponse>,
    /// Original contact owner result, including cone, complementarity and
    /// passivity residuals. A velocity-only witness does not identify local loads.
    pub contact: R,
    /// Actual returned-state relative velocities, [normal, tangent, n cross t].
    pub relative_velocities_after: Vec<[f64; 3]>,
    /// Full row-major 3N-by-3N J M^-1 J^T, in inverse kilograms.
    pub delassus: Vec<f64>,
    /// Signed work supplied by prescribed translating supports, in joules.
    pub plane_work: f64,
    /// Actual kinetic change minus plane work plus admitted loss, in joules.
    pub energy_balance_residual: f64,
    /// Actual global momentum change minus external support impulse, in N s.
    pub linear_momentum_residual: Vec3,
    /// Actual global angular change minus external support angular impulse,
    /// about the world origin, in N m s.
    pub angular_momentum_residual: Vec3,
    /// Maximum actual versus predicted point-velocity discrepancy, in m/s.
    pub endpoint_velocity_residual: f64,
}

/// Actual original material-point kinematics before an articulated impact.
///
/// The frame is the exact normalized [normal, tangent, normal cross tangent]
/// basis used by the native contact assembly. A plane counterpart point is the
/// body's normal projection onto that plane, not its arbitrary defining anchor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArticulatedImpactKinematics {
    pub original_index: usize,
    pub contact_id: u64,
    pub axes_world: [Vec3; 3],
    pub gap_m: f64,
    pub relative_velocity_m_s: [f64; 3],
    pub point_world: Vec3,
    pub counterpart_point_world: Vec3,
}

/// Resolve supplied frictional contacts over freely moving articulated trees.
///
/// The full native mobility retains normal/tangent, joint/base and cross-contact
/// coupling. One simultaneous solve precedes one accumulated impulse update per
/// system. Material contact points are never projected. Dynamic witnesses must
/// share a normal line; any finite-gap tangent couple still has to pass the
/// independent global angular-momentum gate.
///
/// Work is bounded by the existing contact owner's 21-contact cap and at most
/// 128 systems. Each row probes only the one or two involved articulated trees
/// with linear-time ABA; no generalized dense mass matrix is assembled. Fixed
/// poses, borrowed state, and private scratch give atomic refusal and retry.
///
/// This is an instantaneous point-contact law, with a single Coulomb coefficient
/// and explicit restitution/target semantics. Detection, time stepping, static
/// versus kinetic coefficient transitions, finite-patch spin friction and the
/// disposition of impact loss into heat or deformation remain separate owners.
pub fn resolve_articulated_frictional_contacts(
    systems: &[ArticulatedContactSystem<'_>],
    contacts: &[ArticulatedFrictionalContact],
    config: FrictionalContactConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedFrictionalContactResult, ArticulatedContactError> {
    let PreparedContactEvent {
        mut workspace,
        prepared,
        laws,
        order,
    } = prepare_event(systems, contacts, config, &mut cancelled)?;
    let delassus = workspace.delassus(systems, &prepared, &mut cancelled)?;
    let contact = solve_frictional_contacts(&laws, &delassus, config, &mut cancelled)?;
    admit_support_work(&contact, contacts, &order)?;
    let impulses: Vec<_> = contact.impulses.iter().flatten().copied().collect();
    let predicted_speeds: Vec<_> = contact
        .relative_velocities_after
        .iter()
        .flatten()
        .copied()
        .collect();
    let mechanics = workspace.apply(
        systems,
        &prepared,
        EventResponse {
            impulses: &impulses,
            predicted_speeds: &predicted_speeds,
            stored_energy_change: 0.0,
            dissipated_energy: contact.dissipated_energy,
            velocity_tolerance: config.velocity_tolerance,
            energy_tolerance: config.energy_tolerance,
        },
        &mut cancelled,
    )?;
    Ok(ArticulatedFrictionalContactResult {
        contact_ids: order.iter().map(|&index| contacts[index].id).collect(),
        responses: mechanics.responses,
        contact,
        relative_velocities_after: mechanics
            .speeds_after
            .chunks_exact(3)
            .map(|row| [row[0], row[1], row[2]])
            .collect(),
        delassus,
        plane_work: mechanics.plane_work,
        energy_balance_residual: mechanics.energy_balance_residual,
        linear_momentum_residual: mechanics.linear_momentum_residual,
        angular_momentum_residual: mechanics.angular_momentum_residual,
        endpoint_velocity_residual: mechanics.endpoint_velocity_residual,
    })
}

struct PreparedContactEvent<L = FrictionalContact> {
    workspace: EventWorkspace,
    prepared: Vec<Prepared>,
    laws: Vec<L>,
    order: Vec<usize>,
}

fn prepare_event(
    systems: &[ArticulatedContactSystem<'_>],
    contacts: &[ArticulatedFrictionalContact],
    config: FrictionalContactConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedContactEvent, ArticulatedContactError> {
    prepare_event_with(systems, contacts, config, cancelled, |query| {
        let c = contacts[query.original_index];
        Ok(FrictionalContact {
            id: c.id,
            gap: query.gap_m,
            activation_tolerance: c.activation_tolerance,
            relative_velocity: query.relative_velocity_m_s,
            friction_coefficient: c.friction_coefficient,
            normal_law: c.normal_law,
        })
    })
}

/// Assemble geometry and native point rows once, independent of coefficient
/// source. The law callback only receives original incoming kinematics; the
/// complete native response is still admitted by EventWorkspace::apply.
fn prepare_event_with<F: Copy, L, E: From<ArticulatedContactError>>(
    systems: &[ArticulatedContactSystem<'_>],
    contacts: &[ArticulatedFrictionalContact<F>],
    config: FrictionalContactConfig,
    cancelled: &mut impl FnMut() -> bool,
    mut make_law: impl FnMut(&ArticulatedImpactKinematics) -> Result<L, E>,
) -> Result<PreparedContactEvent<L>, E> {
    poll(cancelled)?;
    if contacts.len() > 21 || contacts.len() > config.maximum_contacts || systems.len() > 128 {
        return Err(ArticulatedContactError::Invalid(
            "articulated frictional event exceeds its size cap",
        )
        .into());
    }
    let mut order: Vec<_> = (0..contacts.len()).collect();
    order.sort_unstable_by_key(|&index| contacts[index].id);
    if order
        .windows(2)
        .any(|pair| contacts[pair[0]].id == contacts[pair[1]].id)
    {
        return Err(ArticulatedContactError::Invalid("duplicate contact identity").into());
    }
    let workspace = EventWorkspace::new(systems, cancelled)?;
    let mut prepared = Vec::with_capacity(3 * contacts.len());
    let mut laws = Vec::with_capacity(contacts.len());
    for &index in &order {
        poll(cancelled)?;
        let c = contacts[index];
        let geometry = prepare_contact_frame(
            c.point,
            c.other,
            c.normal_world,
            c.tangent_world,
            index,
            &workspace.before,
        )
        .map_err(ArticulatedContactError::from)?;
        for (axis_index, row) in geometry.rows.iter().enumerate() {
            prepared.push(Prepared {
                first: row.first,
                second: row.second,
                axis_world: row.axis_world,
                world_point: row.world_point,
                plane_speed: row.plane_speed,
                minimum_speed: (axis_index == 0 && row.gap <= c.activation_tolerance)
                    .then_some(0.0),
            });
        }
        laws.push(make_law(&ArticulatedImpactKinematics {
            original_index: index,
            contact_id: c.id,
            axes_world: geometry.rows.map(|row| row.axis_world),
            gap_m: geometry.rows[0].gap,
            relative_velocity_m_s: geometry.rows.map(|row| row.relative_normal_speed),
            point_world: geometry.rows[0].world_point,
            counterpart_point_world: geometry.counterpart_point_world,
        })?);
        poll(cancelled)?;
    }
    Ok(PreparedContactEvent {
        workspace,
        prepared,
        laws,
        order,
    })
}

fn admit_support_work<F>(
    contact: &FrictionalContactResult,
    contacts: &[ArticulatedFrictionalContact<F>],
    order: &[usize],
) -> Result<(), ArticulatedContactError> {
    if matches!(
        contact.determinacy,
        FrictionalContactDeterminacy::UniqueVelocity(_)
    ) {
        // For a native mechanical Gram matrix a null contact impulse changes
        // no body momentum. Its total support impulse is therefore zero, as
        // is its support work at one common translation. Independently moving
        // supports need an additional work-determinacy witness, not a guess.
        let mut common_support_velocity = None;
        for &active in &contact.active_indices {
            if let ContactCounterpart::PrescribedPlane { velocity_world, .. } =
                contacts[order[active]].other
            {
                if common_support_velocity.is_some_and(|velocity| velocity != velocity_world) {
                    return Err(ArticulatedContactError::UnresolvedSupportWork);
                }
                common_support_velocity = Some(velocity_world);
            }
        }
    }
    Ok(())
}
