//! Coupled normal and two-tangent Coulomb impulses on articulated trees.
//!
//! This adapter reuses the normal adapter's link Jacobians, original ABA
//! inverse-mass actions, accumulated impulse application, and physical endpoint
//! audit. The existing contact owner supplies the complete non-associated
//! Signorini--Coulomb solve; no independent link masses or sequential tangent
//! clipping are substituted for articulated recoil.

use super::{
    ArticulatedContactError, ArticulatedContactSystem, ContactCounterpart, EventResponse,
    EventWorkspace, LinkContactPoint, Prepared, finite, magnitude, poll, prepare_point,
    unit_normal, vector,
};
use crate::articulated::FreeFloatingImpulseResponse;
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactDeterminacy,
    FrictionalContactResult, NormalVelocityLaw, solve_frictional_contacts,
};
use fs_ga::Vec3;

/// One supplied geometric witness with an explicit ideal Coulomb law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArticulatedFrictionalContact {
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
    /// Common nonnegative static/kinetic Coulomb coefficient. This event law
    /// does not infer velocity dependence or finite-patch microslip.
    pub friction_coefficient: f64,
    /// Explicit normal endpoint law owned by `fs-contact`.
    pub normal_law: NormalVelocityLaw,
}

/// All native system responses and independently reconstructed balances.
#[derive(Debug, Clone)]
pub struct ArticulatedFrictionalContactResult {
    /// Canonical contact order for every contact-space vector and matrix below.
    pub contact_ids: Vec<u64>,
    /// Updated states and physical momenta, in original system order.
    pub responses: Vec<FreeFloatingImpulseResponse>,
    /// Original contact owner result, including cone, complementarity and
    /// passivity residuals. A velocity-only witness does not identify local loads.
    pub contact: FrictionalContactResult,
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
    poll(&mut cancelled)?;
    if contacts.len() > 21 || contacts.len() > config.maximum_contacts || systems.len() > 128 {
        return Err(ArticulatedContactError::Invalid(
            "articulated frictional event exceeds its size cap",
        ));
    }
    let mut order: Vec<_> = (0..contacts.len()).collect();
    order.sort_unstable_by_key(|&index| contacts[index].id);
    if order
        .windows(2)
        .any(|pair| contacts[pair[0]].id == contacts[pair[1]].id)
    {
        return Err(ArticulatedContactError::Invalid(
            "duplicate contact identity",
        ));
    }
    let mut workspace = EventWorkspace::new(systems, &mut cancelled)?;
    let mut prepared = Vec::with_capacity(3 * contacts.len());
    let mut laws = Vec::with_capacity(contacts.len());
    for &index in &order {
        poll(&mut cancelled)?;
        let c = contacts[index];
        let normal = unit_normal(c.normal_world)?;
        // Normalize before crossing to avoid overflow in arbitrary finite hints.
        // Two cross products remove the parallel part without magnifying the
        // subtraction error of nearly parallel unit vectors.
        let hint = unit_normal(c.tangent_world)?;
        let second_tangent = unit_normal(normal.cross(hint))?;
        let first_tangent = unit_normal(second_tangent.cross(normal))?;
        let axes = [normal, first_tangent, second_tangent];
        let mut gap = 0.0;
        let mut relative_velocity = [0.0; 3];
        for (axis_index, axis) in axes.into_iter().enumerate() {
            let (first, world_point) = prepare_point(c.point, axis, &workspace.before)?;
            let (second, other_point, plane_speed) = match c.other {
                ContactCounterpart::Body(point) => {
                    if point.system_index == c.point.system_index
                        && point.link_index == c.point.link_index
                    {
                        return Err(ArticulatedContactError::Invalid(
                            "contact requires distinct rigid links",
                        ));
                    }
                    let (term, position) =
                        prepare_point(point, axis.scale(-1.0), &workspace.before)?;
                    (Some(term), position, 0.0)
                }
                ContactCounterpart::PrescribedPlane {
                    point_world,
                    velocity_world,
                } => {
                    vector(point_world)?;
                    vector(velocity_world)?;
                    (None, point_world, finite(axis.dot(velocity_world))?)
                }
            };
            if axis_index == 0 {
                let separation = world_point - other_point;
                gap = finite(normal.dot(separation))?;
                if second.is_some() {
                    let tangent = separation - normal.scale(gap);
                    let scale = magnitude(world_point)
                        .max(magnitude(other_point))
                        .max(magnitude(separation));
                    if magnitude(tangent) > 256.0 * f64::EPSILON * scale {
                        return Err(ArticulatedContactError::NonCollinearContact { index });
                    }
                }
            }
            let row = Prepared {
                first,
                second,
                axis_world: axis,
                world_point,
                plane_speed,
                minimum_speed: (axis_index == 0 && gap <= c.activation_tolerance).then_some(0.0),
            };
            relative_velocity[axis_index] = finite(row.velocity(&workspace.before))?;
            prepared.push(row);
        }
        laws.push(FrictionalContact {
            id: c.id,
            gap,
            activation_tolerance: c.activation_tolerance,
            relative_velocity,
            friction_coefficient: c.friction_coefficient,
            normal_law: c.normal_law,
        });
    }
    let delassus = workspace.delassus(systems, &prepared, &mut cancelled)?;
    let contact = solve_frictional_contacts(&laws, &delassus, config, &mut cancelled)?;
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
