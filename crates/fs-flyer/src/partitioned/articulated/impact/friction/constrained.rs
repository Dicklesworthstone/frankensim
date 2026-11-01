//! Bilateral linkages during impacts with complete dry and attached-fluid inertia.
//!
//! The original native geometry supplies distance, spherical and welded rows.
//! Its normalized Schur owner condenses those rows through the same prepared
//! physical mass used by the contact law and the final accumulated impulse.

use super::{
    ArticulatedFrictionalContact, ArticulatedFrictionalImpact, ArticulatedImpactSystem,
    ImpactWorkspace, PointImpactRow, admit_support_work, contact_order, law_refusal, poll,
    prepare_contacts, refusal,
};
use crate::Refusal;
use crate::addedmass::AeroGeneralizedLoads;
use fs_contact::frictional_multicontact::{FrictionalContactConfig, solve_frictional_contacts};
use fs_mbd::articulated::generalized::GeneralizedKinematics;
use fs_mbd::articulated::motion::MAX_FREE_FLOATING_STEP_LINKS;
use fs_mbd::articulated_contacts::ArticulatedContactError;
pub use fs_mbd::articulated_contacts::friction::constrained::bilateral::{
    ImpactBilateralConstraint, ImpactBilateralRow, ImpactBilateralRowKind, ImpactPointConstraint,
    ImpactWeldConstraint, MAX_IMPACT_BILATERAL_ROWS,
};
use fs_mbd::articulated_contacts::friction::constrained::operators::{
    ConstrainedMobility, prepare_bilateral_rows, validate_constraint_config,
};
pub use fs_mbd::articulated_contacts::friction::constrained::{
    ConstrainedImpactConfig, ImpactDistanceConstraint,
};

/// A simultaneous contact and bilateral event accepted by the full inertia.
#[derive(Debug, Clone)]
pub struct ConstrainedArticulatedImpact {
    /// Original law and complete physical endpoint audit. Its Delassus matrix
    /// is the bilateral Schur complement; its row speeds cover contact rows.
    pub impact: ArticulatedFrictionalImpact,
    /// Original constraint index and physical interpretation of each reaction.
    pub rows: Vec<ImpactBilateralRow>,
    /// Signed reactions in row order: N s for distance/point rows and N m s
    /// for weld angular rows. Velocity-only contact determinacy does not make
    /// these individual reactions unique.
    pub bilateral_impulses: Vec<f64>,
    /// Supplied distance/point errors, in metres. Angular weld rows are zero
    /// because the instantaneous event preserves the given relative pose.
    pub position_residual_m: Vec<f64>,
    /// Supplied objective rates: m/s for distance/point, rad/s for angular rows.
    pub velocity_before: Vec<f64>,
    /// Actual returned rates measured from native link twists, in matching units.
    pub velocity_after: Vec<f64>,
    /// Actual bilateral work from mean endpoint rates, in J. This numerical
    /// tangency defect remains in the global energy audit and is never heat.
    pub constraint_work_j: f64,
    /// Minimum squared Cholesky pivot of normalized bilateral mobility.
    pub minimum_scaled_pivot: f64,
}

fn constraint_refusal(error: ArticulatedContactError) -> Refusal {
    let code = if matches!(error, ArticulatedContactError::Cancelled) {
        "articulated-friction-cancelled"
    } else {
        "articulated-constraint-invalid"
    };
    refusal(code, error)
}

fn finite(value: f64) -> Result<f64, Refusal> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(refusal(
            "articulated-constraint-invalid",
            "nonfinite constrained impact arithmetic",
        ))
    }
}

/// Resolve Coulomb impacts while preserving admitted distance, point and weld
/// constraints through each articulated tree's complete dry-plus-fluid inertia.
///
/// Contact and constraint covectors form one augmented mobility. The original
/// native constraint owner factors its normalized bilateral block once,
/// condenses the contact mobility and recovers every bilateral reaction. The
/// common impact workspace then applies all impulses together. It retains
/// base/joint, normal/tangent, cross-contact and linkage recoil without a
/// sequential velocity projection, penalty stiffness or fictitious mass.
///
/// Up to 21 contacts and 32 independent scalar bilateral rows (95 rows total)
/// are admitted, with at most 128 trees, 4096 links in total and the existing
/// per-tree limits. Only the bilateral rows measure uncancelled endpoint
/// mobility, using the same complete inertia, so physical added mass cannot be
/// mistaken for a redundant constraint by comparison with a dry-only scale.
///
/// Initial geometry and rates must be admitted. The inherited native
/// `velocity_tolerance_mps` value is a scalar threshold in each row's units:
/// m/s for range/point rows and rad/s for weld angular rows. No pose or rate is
/// repaired. Numerical reaction work remains subject to the full energy gate.
/// Empty constraints delegate to the existing unconstrained entry exactly.
///
/// Configuration and memory stay fixed. This is a caller-located event, not
/// continuous loop integration or collision discovery. Independent per-tree
/// fluid masses do not infer cross-tree hydrodynamic added-mass blocks.
pub fn resolve(
    systems: &[ArticulatedImpactSystem<'_>],
    contacts: &[ArticulatedFrictionalContact],
    constraints: &[ImpactBilateralConstraint],
    config: FrictionalContactConfig,
    constraint_config: ConstrainedImpactConfig,
    evaluate_fluid: impl FnMut(usize, &GeneralizedKinematics) -> Result<AeroGeneralizedLoads, Refusal>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ConstrainedArticulatedImpact, Refusal> {
    poll(&mut cancelled)?;
    let cfg = constraint_config;
    validate_constraint_config(cfg).map_err(constraint_refusal)?;
    if constraints.is_empty() {
        return Ok(ConstrainedArticulatedImpact {
            impact: super::resolve(systems, contacts, config, evaluate_fluid, cancelled)?,
            rows: Vec::new(),
            bilateral_impulses: Vec::new(),
            position_residual_m: Vec::new(),
            velocity_before: Vec::new(),
            velocity_after: Vec::new(),
            constraint_work_j: 0.,
            minimum_scaled_pivot: 1.,
        });
    }
    if constraints.len() > MAX_IMPACT_BILATERAL_ROWS
        || systems.iter().fold(0usize, |count, s| {
            count.saturating_add(s.model.link_count())
        }) > MAX_FREE_FLOATING_STEP_LINKS
    {
        return Err(refusal(
            "articulated-constraint-invalid",
            "constrained impact exceeds its work bound",
        ));
    }
    let k: usize = constraints
        .iter()
        .map(|c| match c {
            ImpactBilateralConstraint::Distance(_) => 1,
            ImpactBilateralConstraint::PointCoincidence(_) => 3,
            ImpactBilateralConstraint::Weld(_) => 6,
        })
        .sum();
    if k > MAX_IMPACT_BILATERAL_ROWS {
        return Err(refusal(
            "articulated-constraint-invalid",
            "more than 32 scalar bilateral rows",
        ));
    }
    let order = contact_order(contacts, config)?;
    let workspace = ImpactWorkspace::new(systems, evaluate_fluid, &mut cancelled)?;
    let kinematics = workspace.kinematics();
    let (mut rows, laws) = prepare_contacts(contacts, &order, &kinematics, &mut cancelled)?;
    let n = rows.len();
    let geometry = prepare_bilateral_rows(constraints, cfg, &kinematics, &mut cancelled)
        .map_err(constraint_refusal)?;
    rows.extend(geometry.iter().map(|g| PointImpactRow {
        first: g.first,
        second: Some(g.second),
        axis_world: g.axis_world,
        world_point: g.world_point,
        plane_speed: 0.,
        relative_speed_before: g.velocity_before,
    }));
    let prepared = workspace.project_bilateral_rows(&rows, n, &mut cancelled)?;
    let projection = ConstrainedMobility::new(
        &prepared.delassus,
        n,
        &prepared.bilateral_reference_mobility,
        cfg.minimum_scaled_pivot,
        &mut cancelled,
    )
    .map_err(constraint_refusal)?;
    let contact =
        solve_frictional_contacts(&laws, &projection.contact_mobility, config, &mut cancelled)
            .map_err(law_refusal)?;
    admit_support_work(&contact, contacts, &order)?;
    let mut impulses: Vec<_> = contact.impulses.iter().flatten().copied().collect();
    let bilateral_impulses = projection
        .recover(&impulses, &mut cancelled)
        .map_err(constraint_refusal)?;
    impulses.extend_from_slice(&bilateral_impulses);
    let mut predicted: Vec<_> = contact
        .relative_velocities_after
        .iter()
        .flatten()
        .copied()
        .collect();
    let dimension = n + k;
    for (i, g) in geometry.iter().enumerate() {
        poll(&mut cancelled)?;
        let increment: f64 = prepared.delassus[(n + i) * dimension..(n + i + 1) * dimension]
            .iter()
            .zip(&impulses)
            .map(|(a, impulse)| a * impulse)
            .sum();
        predicted.push(finite(g.velocity_before + increment)?);
    }
    let mut mechanics = workspace.apply(
        &prepared,
        &impulses,
        &predicted,
        contact.dissipated_energy,
        config.velocity_tolerance.min(cfg.velocity_tolerance_mps),
        config.energy_tolerance,
        &mut cancelled,
    )?;
    let after = mechanics.relative_speeds_after[n..].to_vec();
    if after.iter().any(|v| v.abs() > cfg.velocity_tolerance_mps) {
        return Err(refusal(
            "articulated-constraint-unresolved",
            "actual post-impact bilateral rate is not tangent",
        ));
    }
    let constraint_work = finite(
        bilateral_impulses
            .iter()
            .zip(&geometry)
            .zip(&after)
            .map(|((impulse, g), after)| impulse * 0.5 * (g.velocity_before + after))
            .sum(),
    )?;
    mechanics.relative_speeds_after.truncate(n);
    poll(&mut cancelled)?;
    Ok(ConstrainedArticulatedImpact {
        impact: ArticulatedFrictionalImpact {
            contact_ids: order.iter().map(|&i| contacts[i].id).collect(),
            contact,
            delassus: projection.contact_mobility,
            mechanics,
        },
        rows: geometry.iter().map(|g| g.row).collect(),
        bilateral_impulses,
        position_residual_m: geometry.iter().map(|g| g.position_residual_m).collect(),
        velocity_before: geometry.iter().map(|g| g.velocity_before).collect(),
        velocity_after: after,
        constraint_work_j: constraint_work,
        minimum_scaled_pivot: projection.minimum_scaled_pivot,
    })
}
