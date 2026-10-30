//! Directional Coulomb impulses through the original articulated mechanics.
//!
//! The native contact owner solves every normal and tangent impulse together.
//! Original ABA inverse-mass actions retain base/joint recoil, full link inertia
//! and shared-contact coupling. Actual articulated endpoint momentum and work
//! must pass the same admission as the scalar-friction path.

use super::{
    ArticulatedContactError, ArticulatedContactSystem, ArticulatedFrictionalContact,
    ArticulatedFrictionalContactResult, ArticulatedImpactKinematics, EventResponse,
    FrictionalContactConfig, PreparedContactEvent, admit_support_work, prepare_event_with,
};
use fs_contact::frictional_multicontact::elliptic::{
    EllipticFriction, EllipticFrictionalContact, EllipticFrictionalContactResult,
    solve_elliptic_frictional_contacts,
};

/// Original articulated geometry with one oriented principal friction pair.
pub type ArticulatedEllipticContact = ArticulatedFrictionalContact<EllipticFriction>;
/// Original geometry and normal law, resolved against a source at this event.
pub type ArticulatedImpactGeometry = ArticulatedFrictionalContact<()>;
/// Original mechanical receipts with physical ellipse impulses and coordinates.
pub type ArticulatedEllipticContactResult =
    ArticulatedFrictionalContactResult<EllipticFrictionalContactResult>;

/// One private active-source candidate, retaining its original input mapping.
#[derive(Debug, Clone)]
pub struct ArticulatedSourceCandidate<S> {
    pub contact_id: u64,
    pub original_index: usize,
    pub state: S,
}

/// All returned source states and mechanical responses belong to one accepted
/// event. Borrowed inputs are unchanged; the consumer publishes them together.
#[derive(Debug, Clone)]
pub struct ArticulatedEllipticSourceResult<S> {
    pub mechanical: ArticulatedEllipticContactResult,
    /// Geometrically active contacts only, in canonical contact-ID order.
    pub candidates: Vec<ArticulatedSourceCandidate<S>>,
}

#[derive(Debug)]
pub enum ArticulatedEllipticSourceError<E> {
    Mechanics(ArticulatedContactError),
    Source {
        contact_id: u64,
        original_index: usize,
        source: E,
    },
}
impl<E> From<ArticulatedContactError> for ArticulatedEllipticSourceError<E> {
    fn from(error: ArticulatedContactError) -> Self {
        Self::Mechanics(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for ArticulatedEllipticSourceError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source directional articulated impact: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for ArticulatedEllipticSourceError<E> {}

/// Resolve one supplied simultaneous directional-friction event.
///
/// The coefficient pair and principal axis use the contact's normalized tangent
/// frame. The original contact owner admits positive principal coefficients or
/// the fully frictionless pair; a singly zero coefficient is unsupported.
/// Determinacy evidence remains explicitly in transformed coordinates. Native
/// unique-velocity admission does not make individual pressures or support work
/// unique: independently moving redundant supports retain their refusal.
///
/// No contact detection, finite displacement, timestep, inferred pressure or
/// thermal allocation is added. The caller owns the fixed event configuration
/// and normal endpoint law. Every returned state is an immutable candidate.
pub fn resolve_articulated_elliptic_contacts(
    systems: &[ArticulatedContactSystem<'_>],
    contacts: &[ArticulatedEllipticContact],
    config: FrictionalContactConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedEllipticContactResult, ArticulatedContactError> {
    let event = prepare_event_with(systems, contacts, config, &mut cancelled, |query| {
        let c = contacts[query.original_index];
        Ok::<_, ArticulatedContactError>(EllipticFrictionalContact {
            id: c.id,
            gap: query.gap_m,
            activation_tolerance: c.activation_tolerance,
            relative_velocity: query.relative_velocity_m_s,
            normal_law: c.normal_law,
            friction: c.friction_coefficient,
        })
    })?;
    resolve_prepared(systems, contacts, config, event, &mut cancelled)
}

/// Sample original source histories at actual incoming point velocities, then
/// resolve the complete event through the original articulated owner.
///
/// Geometry and histories share input order. Each geometrically active contact
/// is queried once in canonical ID order. Separating active contacts still
/// require a law because another contact can close them through shared recoil.
/// Inactive contacts neither query nor return changed history. The callback
/// receives exact assembled axes and actual world contact points; prescribed
/// plane coordinates are normal projections rather than arbitrary anchors.
///
/// Coefficients are explicitly held at this preimpact source query; this is not
/// an inner velocity-dependent impact iteration. Existing objective transport
/// must express retained texture in the supplied frame before sampling.
///
/// All source candidates remain private until native impulses, actual endpoint
/// kinematics, momentum and work are admitted. Failure or cancellation returns
/// no mechanical or history candidate and leaves every borrowed input intact.
/// External callback state must be staged by its owner. No clone bound on S or
/// wholesale clone of accepted history is needed to stage this transaction.
#[allow(clippy::too_many_arguments)]
pub fn resolve_articulated_elliptic_contacts_with<S, E>(
    systems: &[ArticulatedContactSystem<'_>],
    contacts: &[ArticulatedImpactGeometry],
    histories: &[S],
    config: FrictionalContactConfig,
    mut resolve: impl FnMut(&S, &ArticulatedImpactKinematics) -> Result<(EllipticFriction, S), E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedEllipticSourceResult<S>, ArticulatedEllipticSourceError<E>> {
    if histories.len() != contacts.len() {
        return Err(ArticulatedContactError::Invalid(
            "one accepted source history per original contact required",
        )
        .into());
    }
    let mut candidates = Vec::new();
    let event = prepare_event_with(
        systems,
        contacts,
        config,
        &mut cancelled,
        |query| -> Result<_, ArticulatedEllipticSourceError<E>> {
            let original_index = query.original_index;
            let c = contacts[original_index];
            if query.gap_m < 0.0
                || !c.activation_tolerance.is_finite()
                || c.activation_tolerance < 0.0
            {
                return Err(ArticulatedContactError::Invalid(
                    "nonnegative source contact gap and activation tolerance required",
                )
                .into());
            }
            let friction = if query.gap_m <= c.activation_tolerance {
                let (friction, state) =
                    resolve(&histories[original_index], query).map_err(|source| {
                        ArticulatedEllipticSourceError::Source {
                            contact_id: c.id,
                            original_index,
                            source,
                        }
                    })?;
                candidates.push(ArticulatedSourceCandidate {
                    contact_id: c.id,
                    original_index,
                    state,
                });
                friction
            } else {
                // This geometric row is excluded by the native contact owner.
                // Its actual endpoint still includes the other contacts' recoil.
                EllipticFriction {
                    coefficients_mu: [0.0; 2],
                    principal_axis: [1.0, 0.0],
                }
            };
            Ok(EllipticFrictionalContact {
                id: c.id,
                gap: query.gap_m,
                activation_tolerance: c.activation_tolerance,
                relative_velocity: query.relative_velocity_m_s,
                normal_law: c.normal_law,
                friction,
            })
        },
    )?;
    let mechanical = resolve_prepared(systems, contacts, config, event, &mut cancelled)?;
    Ok(ArticulatedEllipticSourceResult {
        mechanical,
        candidates,
    })
}

fn resolve_prepared<F>(
    systems: &[ArticulatedContactSystem<'_>],
    contacts: &[ArticulatedFrictionalContact<F>],
    config: FrictionalContactConfig,
    event: PreparedContactEvent<EllipticFrictionalContact>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ArticulatedEllipticContactResult, ArticulatedContactError> {
    let PreparedContactEvent {
        mut workspace,
        prepared,
        laws,
        order,
    } = event;
    let delassus = workspace.delassus(systems, &prepared, cancelled)?;
    let contact = solve_elliptic_frictional_contacts(&laws, &delassus, config, &mut *cancelled)?;
    admit_support_work(&contact.transformed, contacts, &order)?;
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
        cancelled,
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
