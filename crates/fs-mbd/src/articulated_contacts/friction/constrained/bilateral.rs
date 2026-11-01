//! Spherical and welded bilateral reactions during simultaneous articulated impacts.
//!
//! This extends the distance-only constrained impact owner with material-point
//! coincidence and weld constraints. All bilateral rows are condensed through
//! the same native articulated inverse-mass action before the existing Coulomb
//! solve, then contact and constraint impulses are accumulated and applied once.
//! No post-impact velocity projection or dense generalized mass matrix is used.

use super::super::{
    ArticulatedContactError, ArticulatedContactSystem, ArticulatedFrictionalContact,
    ArticulatedFrictionalContactResult, EventResponse, FrictionalContactConfig, LinkContactPoint,
    PreparedContactEvent, admit_support_work, finite, poll, prepare_event,
    solve_frictional_contacts,
};
use super::operators::{
    ConstrainedMobility, native_reference_mobility, prepare_bilateral_rows,
    validate_constraint_config,
};
use super::{
    ConstrainedImpactConfig, ImpactDistanceConstraint,
    resolve_constrained_articulated_frictional_contacts,
};
use crate::articulated::motion::MAX_FREE_FLOATING_STEP_LINKS;

/// Maximum independently reacting scalar bilateral rows in one impact.
pub const MAX_IMPACT_BILATERAL_ROWS: usize = 32;

/// Coincidence of two material points while relative rotation remains free.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImpactPointConstraint {
    /// First original material point.
    pub a: LinkContactPoint,
    /// Second original material point on a distinct rigid link.
    pub b: LinkContactPoint,
}

/// Coincidence of two material points plus zero relative angular velocity.
///
/// The instantaneous event preserves the supplied relative orientation. It does
/// not define a new rest orientation or repair an inconsistent pose.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImpactWeldConstraint {
    /// First original material point and reference material frame.
    pub a: LinkContactPoint,
    /// Second original material point on a distinct rigid link.
    pub b: LinkContactPoint,
}

/// One bilateral constraint participating in the simultaneous impact solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ImpactBilateralConstraint {
    /// One scalar positive-distance row.
    Distance(ImpactDistanceConstraint),
    /// Three material-frame point-coincidence rows.
    PointCoincidence(ImpactPointConstraint),
    /// Three point rows plus three relative-angular-velocity rows.
    Weld(ImpactWeldConstraint),
}

/// Interpretation of one flattened scalar bilateral row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImpactBilateralRowKind {
    /// Positive range impulse for a distance constraint.
    Distance,
    /// Point-coincidence component in A's material frame.
    Point { component: usize },
    /// Weld point component in A's material frame.
    WeldPoint { component: usize },
    /// Weld relative-angular-velocity component in A's material frame.
    WeldAngular { component: usize },
}

/// One flattened row and its original caller constraint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ImpactBilateralRow {
    /// Index in the caller's bilateral-constraint slice.
    pub constraint: usize,
    /// Physical interpretation of this scalar row.
    pub kind: ImpactBilateralRowKind,
}

/// One globally admitted contact plus general bilateral impulse event.
#[derive(Debug, Clone)]
pub struct BilateralArticulatedContactResult {
    /// Original constrained contact result over the condensed contact mobility.
    pub impact: ArticulatedFrictionalContactResult,
    /// Flattened bilateral row metadata.
    pub rows: Vec<ImpactBilateralRow>,
    /// Signed impulse per flattened row, in N s.
    pub bilateral_impulses_ns: Vec<f64>,
    /// Initial geometric residual per row. Angular weld rows are zero because
    /// the event preserves the current supplied relative orientation.
    pub position_residual_m: Vec<f64>,
    /// Initial objective bilateral rate per row, in m/s for point/distance rows
    /// and rad/s for weld angular rows.
    pub velocity_before: Vec<f64>,
    /// Actual returned objective bilateral rate per row, with matching units.
    pub velocity_after: Vec<f64>,
    /// Actual bilateral impulse work using average before/after row rates, J.
    pub constraint_work_j: f64,
    /// Smallest squared Cholesky pivot of the normalized bilateral mobility.
    pub minimum_scaled_pivot: f64,
}

/// Resolve Coulomb impacts while preserving distance, spherical, and weld rows.
///
/// Spherical rows use A-material-frame point coordinates. Away from exact
/// coincidence their A-side reaction includes the geometric couple required by
/// the differential of `R_A^T(x_B-x_A)`, so internal reactions retain zero net
/// spatial impulse. Weld angular rows use equal/opposite pure couples. The
/// existing contact law sees only the Schur-condensed contact mobility.
///
/// Existing distance-only callers should keep using
/// `resolve_constrained_articulated_frictional_contacts`; this generalized
/// owner preserves that path exactly when every row is a distance.
///
/// # Errors
/// Refuses inconsistent supplied geometry/rates, dependent or ill-conditioned
/// bilateral rows, native contact/dynamics failures, exhausted bounds, or
/// cancellation. No pose or velocity is projected.
pub fn resolve_bilateral_articulated_frictional_contacts(
    systems: &[ArticulatedContactSystem<'_>],
    contacts: &[ArticulatedFrictionalContact],
    constraints: &[ImpactBilateralConstraint],
    config: FrictionalContactConfig,
    constraint_config: ConstrainedImpactConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<BilateralArticulatedContactResult, ArticulatedContactError> {
    poll(&mut cancelled)?;
    let cfg = constraint_config;
    validate_constraint_config(cfg)?;

    if constraints
        .iter()
        .all(|c| matches!(c, ImpactBilateralConstraint::Distance(_)))
    {
        let distances: Vec<_> = constraints
            .iter()
            .map(|c| match *c {
                ImpactBilateralConstraint::Distance(d) => d,
                _ => unreachable!(),
            })
            .collect();
        let old = resolve_constrained_articulated_frictional_contacts(
            systems, contacts, &distances, config, cfg, cancelled,
        )?;
        return Ok(BilateralArticulatedContactResult {
            impact: old.impact,
            rows: (0..distances.len())
                .map(|constraint| ImpactBilateralRow {
                    constraint,
                    kind: ImpactBilateralRowKind::Distance,
                })
                .collect(),
            bilateral_impulses_ns: old.distance_impulses_ns,
            position_residual_m: old.length_residual_m,
            velocity_before: old.velocity_before_mps,
            velocity_after: old.velocity_after_mps,
            constraint_work_j: old.constraint_work_j,
            minimum_scaled_pivot: old.minimum_scaled_pivot,
        });
    }

    let scalar_rows: usize = constraints
        .iter()
        .map(|c| match c {
            ImpactBilateralConstraint::Distance(_) => 1,
            ImpactBilateralConstraint::PointCoincidence(_) => 3,
            ImpactBilateralConstraint::Weld(_) => 6,
        })
        .sum();
    if scalar_rows > MAX_IMPACT_BILATERAL_ROWS
        || systems.iter().fold(0usize, |count, s| {
            count.saturating_add(s.model.link_count())
        }) > MAX_FREE_FLOATING_STEP_LINKS
    {
        return Err(ArticulatedContactError::Invalid(
            "bilateral-impact work bound exceeded",
        ));
    }

    let PreparedContactEvent {
        mut workspace,
        mut prepared,
        laws,
        order,
    } = prepare_event(systems, contacts, config, &mut cancelled)?;
    let n = prepared.len();
    let geometry = prepare_bilateral_rows(constraints, cfg, &workspace.before, &mut cancelled)?;
    let rows: Vec<_> = geometry.iter().map(|r| r.row).collect();
    let residuals: Vec<_> = geometry.iter().map(|r| r.position_residual_m).collect();
    let rates: Vec<_> = geometry.iter().map(|r| r.velocity_before).collect();
    prepared.extend(geometry.iter().map(|r| r.prepared()));

    let k = rows.len();
    if k == 0 {
        return Err(ArticulatedContactError::Invalid(
            "mixed bilateral impact unexpectedly has no rows",
        ));
    }
    let dimension = n + k;
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
    let bilateral_impulses = projection.recover(&impulses, &mut cancelled)?;
    impulses.extend_from_slice(&bilateral_impulses);

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
            .map(|(x, impulse)| x * impulse)
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
    if after
        .iter()
        .zip(&rows)
        .any(|(v, _)| v.abs() > cfg.velocity_tolerance_mps)
    {
        return Err(ArticulatedContactError::Unresolved(
            "actual post-impact bilateral tangency",
        ));
    }
    let constraint_work = finite(
        bilateral_impulses
            .iter()
            .zip(&rates)
            .zip(&after)
            .map(|((impulse, before), after)| impulse * 0.5 * (before + after))
            .sum(),
    )?;
    poll(&mut cancelled)?;
    Ok(BilateralArticulatedContactResult {
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
            energy_balance_residual: mechanics.energy_balance_residual,
            linear_momentum_residual: mechanics.linear_momentum_residual,
            angular_momentum_residual: mechanics.angular_momentum_residual,
            endpoint_velocity_residual: mechanics.endpoint_velocity_residual,
        },
        rows,
        bilateral_impulses_ns: bilateral_impulses,
        position_residual_m: residuals,
        velocity_before: rates,
        velocity_after: after,
        constraint_work_j: constraint_work,
        minimum_scaled_pivot: projection.minimum_scaled_pivot,
    })
}
