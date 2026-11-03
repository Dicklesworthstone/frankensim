//! Simultaneous Coulomb impacts with rigid internal distance reactions.
//!
//! The original native impulse response builds the augmented Gram matrix
//! `[D B; B^T A]` for contact and bilateral rows. Condensing the independent
//! bilateral mobility gives `D_c = D - B A^-1 B^T`. The existing contact law
//! solves on this constrained mobility, then bilateral impulses are recovered
//! and ALL impulses are applied together by the original mechanics owner.
//! This is the auxiliary-constraint construction of Baraff (1996), section 8:
//! <https://www.cs.cmu.edu/~baraff/papers/sig96.pdf>.
//!
//! Input distances and rates must already be admitted. An instantaneous impact
//! preserves their tangent velocities; it does not repair drift or locate an
//! event. Actual finite constraint work is retained, and the native global
//! passivity, momentum and endpoint gates remain authoritative.

/// Spherical point closures and weld reactions sharing this impact solve.
pub mod bilateral;

/// Reusable bilateral geometry and normalized contact-mobility condensation.
pub mod operators;

use super::{
    ArticulatedContactError, ArticulatedContactSystem, ArticulatedFrictionalContact,
    ArticulatedFrictionalContactResult, EventResponse, FrictionalContactConfig, LinkContactPoint,
    PreparedContactEvent, admit_support_work, finite, poll, prepare_event,
    resolve_articulated_frictional_contacts, solve_frictional_contacts,
};
use crate::articulated::motion::MAX_FREE_FLOATING_STEP_LINKS;
use operators::{ConstrainedMobility, native_reference_mobility, prepare_bilateral_rows};

/// Maximum independently reacting bilateral distance rows in one impact.
pub const MAX_IMPACT_DISTANCE_CONSTRAINTS: usize = 32;

/// A positive rigid distance between material points, within or across trees.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImpactDistanceConstraint {
    /// First original material point; receives the negative range impulse.
    pub a: LinkContactPoint,
    /// Second distinct link; receives the positive range impulse.
    pub b: LinkContactPoint,
    /// Required positive distance, in metres.
    pub length_m: f64,
}

/// Physical initial/endpoint and independent-mobility admission.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConstrainedImpactConfig {
    /// Absolute allowed supplied distance error, in metres.
    pub position_tolerance_m: f64,
    /// Allowed supplied and returned range rate, in m/s.
    pub velocity_tolerance_mps: f64,
    /// Minimum mobility/free-endpoint ratio and squared unit-diagonal
    /// Cholesky pivot, in [4096 eps,1). A floating regularity guard only.
    pub minimum_scaled_pivot: f64,
}
impl Default for ConstrainedImpactConfig {
    fn default() -> Self {
        Self {
            position_tolerance_m: 1e-9,
            velocity_tolerance_mps: 1e-9,
            minimum_scaled_pivot: 1e-10,
        }
    }
}

/// One globally admitted native contact/constraint impulse event.
#[derive(Debug, Clone)]
pub struct ConstrainedArticulatedContactResult {
    /// Original contact result. Its mobility is the condensed constrained
    /// mobility; native responses include both contact and bilateral impulses.
    /// Its energy residual includes any numerical bilateral work.
    pub impact: ArticulatedFrictionalContactResult,
    /// Bilateral impulses in input constraint order, N s; positive separates.
    /// If the contact law identifies only velocity, these reactions likewise
    /// need not be unique. The original determinacy verdict remains explicit.
    pub distance_impulses_ns: Vec<f64>,
    /// Actual length errors, unchanged by this instantaneous event, in metres.
    pub length_residual_m: Vec<f64>,
    /// Actual supplied material-point range rates, m/s.
    pub velocity_before_mps: Vec<f64>,
    /// Actual range rates reconstructed after the accumulated native event, m/s.
    pub velocity_after_mps: Vec<f64>,
    /// Actual bilateral impulse work using average before/after range rates, J.
    /// Zero for exact tangency in real arithmetic; never labeled physical heat.
    pub constraint_work_j: f64,
    /// Smallest squared Cholesky pivot of normalized bilateral mobility.
    /// One for the empty-constraint delegation.
    pub minimum_scaled_pivot: f64,
}

/// Resolve frictional contacts while preserving admitted rigid distances.
///
/// Contact and constraint reactions share the full original articulated mass
/// response, including cross-tree constraints, internal joints, lever arms,
/// complete inertia tensors and free-base recoil. No sequential contact kick
/// or velocity projection is performed. Existing restitution, friction-cone,
/// nonunique-load, prescribed-support work and global passivity rules apply.
///
/// At most 21 contacts, 32 independent distances, 128 trees and 4096 total links
/// are admitted when constraints are present. Native response work scales with
/// involved links per row; only the small constraint block is factored. Empty
/// constraints delegate to the original contact entry point exactly.
///
/// Input tangency is a precondition, not silently repaired. Roundoff and any
/// admitted initial rate defect remain in the actual constraint-work and global
/// energy residual; the original energy gate may refuse them. This is an event
/// solver, without collision detection, continuous constraint integration,
/// compliant history, spherical/weld closure or reaction uniqueness guarantees.
///
/// # Errors
/// Refuses inconsistent geometry/rates, redundant or ill-conditioned constraints,
/// native state/contact-law/endpoint/balance failures, budgets or cancellation.
pub fn resolve_constrained_articulated_frictional_contacts(
    systems: &[ArticulatedContactSystem<'_>],
    contacts: &[ArticulatedFrictionalContact],
    constraints: &[ImpactDistanceConstraint],
    config: FrictionalContactConfig,
    constraint_config: ConstrainedImpactConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedArticulatedContactResult, ArticulatedContactError> {
    poll(&mut cancelled)?;
    let cfg = constraint_config;
    operators::validate_constraint_config(cfg)?;
    if constraints.is_empty() {
        return Ok(ConstrainedArticulatedContactResult {
            impact: resolve_articulated_frictional_contacts(systems, contacts, config, cancelled)?,
            distance_impulses_ns: Vec::new(),
            length_residual_m: Vec::new(),
            velocity_before_mps: Vec::new(),
            velocity_after_mps: Vec::new(),
            constraint_work_j: 0.,
            minimum_scaled_pivot: 1.,
        });
    }
    if constraints.len() > MAX_IMPACT_DISTANCE_CONSTRAINTS
        || systems.iter().fold(0usize, |count, s| {
            count.saturating_add(s.model.link_count())
        }) > MAX_FREE_FLOATING_STEP_LINKS
    {
        return Err(ArticulatedContactError::Invalid(
            "constrained-impact work bound exceeded",
        ));
    }
    let PreparedContactEvent {
        mut workspace,
        mut prepared,
        laws,
        order,
    } = prepare_event(systems, contacts, config, &mut cancelled)?;
    let n = prepared.len();
    let k = constraints.len();
    let dimension = n + k;
    let constraints: Vec<_> = constraints
        .iter()
        .copied()
        .map(bilateral::ImpactBilateralConstraint::Distance)
        .collect();
    let geometry = prepare_bilateral_rows(&constraints, cfg, &workspace.before, &mut cancelled)?;
    let gaps: Vec<_> = geometry.iter().map(|r| r.position_residual_m).collect();
    let rates: Vec<_> = geometry.iter().map(|r| r.velocity_before).collect();
    prepared.extend(geometry.iter().map(|r| r.prepared()));
    let gram = workspace.delassus(systems, &prepared, &mut cancelled)?;
    let reference = native_reference_mobility(systems, &prepared[n..], &mut cancelled)?;
    let projection = ConstrainedMobility::new(
        &gram,
        n,
        &reference,
        cfg.minimum_scaled_pivot,
        &mut cancelled,
    )?;
    let contact =
        solve_frictional_contacts(&laws, &projection.contact_mobility, config, &mut cancelled)?;
    admit_support_work(&contact, contacts, &order)?;
    let mut impulses: Vec<_> = contact.impulses.iter().flatten().copied().collect();
    let distance_impulses = projection.recover(&impulses, &mut cancelled)?;
    impulses.extend_from_slice(&distance_impulses);
    let mut predicted: Vec<_> = contact
        .relative_velocities_after
        .iter()
        .flatten()
        .copied()
        .collect();
    for i in 0..k {
        let increment: f64 = gram[(n + i) * dimension..(n + i + 1) * dimension]
            .iter()
            .zip(&impulses)
            .map(|(a, b)| a * b)
            .sum();
        predicted.push(finite(rates[i] + increment)?);
    }
    let mechanics = workspace.apply(
        systems,
        &prepared,
        EventResponse {
            impulses: &impulses,
            predicted_speeds: &predicted,
            stored_energy_change: 0.,
            dissipated_energy: contact.dissipated_energy,
            velocity_tolerance: config.velocity_tolerance.min(cfg.velocity_tolerance_mps),
            energy_tolerance: config.energy_tolerance,
        },
        &mut cancelled,
    )?;
    let after = mechanics.speeds_after[n..].to_vec();
    if after.iter().any(|v| v.abs() > cfg.velocity_tolerance_mps) {
        return Err(ArticulatedContactError::Unresolved(
            "actual post-impact rigid distance tangency",
        ));
    }
    let constraint_work = finite(
        distance_impulses
            .iter()
            .zip(&rates)
            .zip(&after)
            .map(|((impulse, before), after)| impulse * 0.5 * (before + after))
            .sum(),
    )?;
    poll(&mut cancelled)?;
    Ok(ConstrainedArticulatedContactResult {
        impact: ArticulatedFrictionalContactResult {
            contact_ids: order.iter().map(|&i| contacts[i].id).collect(),
            responses: mechanics.responses,
            contact,
            relative_velocities_after: mechanics.speeds_after[..n]
                .chunks_exact(3)
                .map(|v| [v[0], v[1], v[2]])
                .collect(),
            delassus: projection.contact_mobility,
            plane_work: mechanics.plane_work,
            support_work: mechanics.support_work,
            boundary_work: mechanics.boundary_work,
            energy_balance_residual: mechanics.energy_balance_residual,
            linear_momentum_residual: mechanics.linear_momentum_residual,
            angular_momentum_residual: mechanics.angular_momentum_residual,
            endpoint_velocity_residual: mechanics.endpoint_velocity_residual,
        },
        distance_impulses_ns: distance_impulses,
        length_residual_m: gaps,
        velocity_before_mps: rates,
        velocity_after_mps: after,
        constraint_work_j: constraint_work,
        minimum_scaled_pivot: projection.minimum_scaled_pivot,
    })
}
