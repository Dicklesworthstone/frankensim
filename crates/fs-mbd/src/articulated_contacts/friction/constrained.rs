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

use super::{
    ArticulatedContactError, ArticulatedContactSystem, ArticulatedFrictionalContact,
    ArticulatedFrictionalContactResult, EventResponse, FrictionalContactConfig, LinkContactPoint,
    Prepared, PreparedContactEvent, admit_support_work, finite, poll, prepare_event, prepare_point,
    resolve_articulated_frictional_contacts, solve_frictional_contacts, unit_normal,
};
use crate::articulated::{motion::MAX_FREE_FLOATING_STEP_LINKS, solve_floating_base_system};
use fs_ga::Vec3;
use fs_la::factor::cholesky;

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
    if [
        cfg.position_tolerance_m,
        cfg.velocity_tolerance_mps,
        cfg.minimum_scaled_pivot,
    ]
    .iter()
    .any(|x| !x.is_finite() || *x <= 0.)
        || cfg.minimum_scaled_pivot < 4096. * f64::EPSILON
        || cfg.minimum_scaled_pivot >= 1.
    {
        return Err(ArticulatedContactError::Invalid(
            "invalid constrained-impact controls",
        ));
    }
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
    let mut gaps = Vec::with_capacity(k);
    let mut rates = Vec::with_capacity(k);
    for c in constraints {
        poll(&mut cancelled)?;
        if !c.length_m.is_finite()
            || c.length_m <= 0.
            || (c.a.system_index == c.b.system_index && c.a.link_index == c.b.link_index)
        {
            return Err(ArticulatedContactError::Invalid(
                "invalid impact distance endpoints/length",
            ));
        }
        // The native point owner validates all indices and material coordinates.
        let (_, a) = prepare_point(c.a, Vec3::new(1., 0., 0.), &workspace.before)?;
        let (_, b) = prepare_point(c.b, Vec3::new(1., 0., 0.), &workspace.before)?;
        let delta = b - a;
        let length = finite(delta.norm())?;
        let gap = finite(length - c.length_m)?;
        if gap.abs() > cfg.position_tolerance_m {
            return Err(ArticulatedContactError::Invalid(
                "initial rigid distance is inconsistent",
            ));
        }
        let normal = unit_normal(delta)?;
        let (first, world_point) = prepare_point(c.b, normal, &workspace.before)?;
        let (second, _) = prepare_point(c.a, normal.scale(-1.), &workspace.before)?;
        let row = Prepared {
            first,
            second: Some(second),
            axis_world: normal,
            world_point,
            plane_speed: 0.,
            minimum_speed: None,
        };
        let rate = finite(row.velocity(&workspace.before))?;
        if rate.abs() > cfg.velocity_tolerance_mps {
            return Err(ArticulatedContactError::Invalid(
                "initial rigid distance velocity is not tangent",
            ));
        }
        gaps.push(gap);
        rates.push(rate);
        prepared.push(row);
    }
    let gram = workspace.delassus(systems, &prepared, &mut cancelled)?;
    let mut scale = Vec::with_capacity(k);
    for i in 0..k {
        poll(&mut cancelled)?;
        let d = gram[(n + i) * dimension + n + i];
        let mut free = 0.;
        for term in prepared[n + i].terms() {
            let system = &systems[term.point.system_index];
            let inertia = system.model.links()[term.point.link_index].inertia();
            let matrix = super::super::dynamics(term.point.system_index, inertia.matrix())?;
            let response = super::super::dynamics(
                term.point.system_index,
                solve_floating_base_system(&matrix, term.direction),
            )?;
            free += term.direction.pairing(response);
        }
        if !free.is_finite() || free <= 0. || d <= 0. || d / free < cfg.minimum_scaled_pivot {
            return Err(ArticulatedContactError::Invalid(
                "rigid distance has unresolved native mobility",
            ));
        }
        scale.push(d.sqrt());
    }
    let mut a = vec![0.; k * k];
    for i in 0..k {
        a[i * k + i] = 1.;
        for j in 0..i {
            let x = (gram[(n + i) * dimension + n + j] / scale[i]) / scale[j];
            let y = (gram[(n + j) * dimension + n + i] / scale[j]) / scale[i];
            if !x.is_finite() || !y.is_finite() || (x - y).abs() > 1e-10 {
                return Err(ArticulatedContactError::Unresolved(
                    "bilateral mobility reciprocity",
                ));
            }
            a[i * k + j] = 0.5 * (x + y);
            a[j * k + i] = a[i * k + j];
        }
    }
    poll(&mut cancelled)?;
    let factor = cholesky(&a, k)
        .map_err(|_| ArticulatedContactError::Invalid("dependent impact distance constraints"))?;
    let mut minimum_pivot = 1.0_f64;
    for i in 0..k {
        let pivot = factor.l(i, i).powi(2);
        if !pivot.is_finite() || pivot < cfg.minimum_scaled_pivot {
            return Err(ArticulatedContactError::Invalid(
                "ill-conditioned impact distance constraints",
            ));
        }
        minimum_pivot = minimum_pivot.min(pivot);
    }
    // Each column is A^-1 B^T for one contact covector. The factor is reused.
    let mut recovery = vec![0.; n * k];
    for column in 0..n {
        poll(&mut cancelled)?;
        let rhs = &mut recovery[column * k..(column + 1) * k];
        for i in 0..k {
            rhs[i] = gram[(n + i) * dimension + column] / scale[i];
        }
        factor.solve(rhs);
        for i in 0..k {
            rhs[i] = finite(rhs[i] / scale[i])?;
        }
    }
    let mut mobility = vec![0.; n * n];
    for i in 0..n {
        poll(&mut cancelled)?;
        for j in 0..n {
            let correction: f64 = (0..k)
                .map(|r| gram[i * dimension + n + r] * recovery[j * k + r])
                .sum();
            mobility[i * n + j] = finite(gram[i * dimension + j] - correction)?;
        }
    }
    let contact = solve_frictional_contacts(&laws, &mobility, config, &mut cancelled)?;
    admit_support_work(&contact, contacts, &order)?;
    let mut impulses: Vec<_> = contact.impulses.iter().flatten().copied().collect();
    let distance_impulses: Vec<_> = (0..k)
        .map(|r| {
            -(0..n)
                .map(|i| recovery[i * k + r] * impulses[i])
                .sum::<f64>()
        })
        .collect();
    for &value in &distance_impulses {
        finite(value)?;
    }
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
            delassus: mobility,
            plane_work: mechanics.plane_work,
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
        minimum_scaled_pivot: minimum_pivot,
    })
}
