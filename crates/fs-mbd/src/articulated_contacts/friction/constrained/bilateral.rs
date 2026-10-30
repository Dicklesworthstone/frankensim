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
    Prepared, PreparedContactEvent, Term, admit_support_work, finite, poll, prepare_event,
    prepare_point, solve_frictional_contacts, unit_normal,
};
use super::{
    ConstrainedImpactConfig, ImpactDistanceConstraint,
    resolve_constrained_articulated_frictional_contacts,
};
use crate::articulated::{motion::MAX_FREE_FLOATING_STEP_LINKS, solve_floating_base_system};
use fs_ga::{Vec3, Wrench};
use fs_la::factor::cholesky;

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

    if constraints.iter().all(|c| matches!(c, ImpactBilateralConstraint::Distance(_))) {
        let distances: Vec<_> = constraints
            .iter()
            .map(|c| match *c {
                ImpactBilateralConstraint::Distance(d) => d,
                _ => unreachable!(),
            })
            .collect();
        let old = resolve_constrained_articulated_frictional_contacts(
            systems,
            contacts,
            &distances,
            config,
            cfg,
            cancelled,
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
    let mut rows = Vec::with_capacity(scalar_rows);
    let mut residuals = Vec::with_capacity(scalar_rows);
    let mut rates = Vec::with_capacity(scalar_rows);

    for (constraint_index, constraint) in constraints.iter().copied().enumerate() {
        poll(&mut cancelled)?;
        match constraint {
            ImpactBilateralConstraint::Distance(c) => {
                if !c.length_m.is_finite()
                    || c.length_m <= 0.
                    || same_link(c.a, c.b)
                {
                    return Err(ArticulatedContactError::Invalid(
                        "invalid impact distance endpoints/length",
                    ));
                }
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
                prepared.push(row);
                rows.push(ImpactBilateralRow {
                    constraint: constraint_index,
                    kind: ImpactBilateralRowKind::Distance,
                });
                residuals.push(gap);
                rates.push(rate);
            }
            ImpactBilateralConstraint::PointCoincidence(c) => {
                append_point_rows(
                    c.a,
                    c.b,
                    constraint_index,
                    false,
                    cfg,
                    &workspace.before,
                    &mut prepared,
                    &mut rows,
                    &mut residuals,
                    &mut rates,
                )?;
            }
            ImpactBilateralConstraint::Weld(c) => {
                append_point_rows(
                    c.a,
                    c.b,
                    constraint_index,
                    true,
                    cfg,
                    &workspace.before,
                    &mut prepared,
                    &mut rows,
                    &mut residuals,
                    &mut rates,
                )?;
                append_weld_angular_rows(
                    c.a,
                    c.b,
                    constraint_index,
                    cfg,
                    &workspace.before,
                    &mut prepared,
                    &mut rows,
                    &mut residuals,
                    &mut rates,
                )?;
            }
        }
    }

    let k = rows.len();
    if k == 0 {
        return Err(ArticulatedContactError::Invalid(
            "mixed bilateral impact unexpectedly has no rows",
        ));
    }
    let dimension = n + k;
    let gram = workspace.delassus(systems, &prepared, &mut cancelled)?;
    let mut scale = Vec::with_capacity(k);
    for i in 0..k {
        poll(&mut cancelled)?;
        let row = &prepared[n + i];
        let d = gram[(n + i) * dimension + n + i];
        let mut free = 0.;
        for term in row.terms() {
            let system = &systems[term.point.system_index];
            let inertia = system.model.links()[term.point.link_index].inertia();
            let matrix = super::super::super::dynamics(term.point.system_index, inertia.matrix())?;
            let response = super::super::super::dynamics(
                term.point.system_index,
                solve_floating_base_system(&matrix, term.direction),
            )?;
            free += term.direction.pairing(response);
        }
        if !free.is_finite() || free <= 0. || d <= 0. || d / free < cfg.minimum_scaled_pivot {
            return Err(ArticulatedContactError::Invalid(
                "bilateral row has unresolved native mobility",
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
        .map_err(|_| ArticulatedContactError::Invalid("dependent impact bilateral constraints"))?;
    let mut minimum_pivot = 1.0_f64;
    for i in 0..k {
        let pivot = factor.l(i, i).powi(2);
        if !pivot.is_finite() || pivot < cfg.minimum_scaled_pivot {
            return Err(ArticulatedContactError::Invalid(
                "ill-conditioned impact bilateral constraints",
            ));
        }
        minimum_pivot = minimum_pivot.min(pivot);
    }

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
    let bilateral_impulses: Vec<_> = (0..k)
        .map(|r| {
            -(0..n)
                .map(|i| recovery[i * k + r] * impulses[i])
                .sum::<f64>()
        })
        .collect();
    for &value in &bilateral_impulses {
        finite(value)?;
    }
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
            delassus: mobility,
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
        minimum_scaled_pivot: minimum_pivot,
    })
}

#[allow(clippy::too_many_arguments)]
fn append_point_rows(
    a: LinkContactPoint,
    b: LinkContactPoint,
    constraint_index: usize,
    weld: bool,
    cfg: ConstrainedImpactConfig,
    kinematics: &[crate::articulated::Kinematics],
    prepared: &mut Vec<Prepared>,
    rows: &mut Vec<ImpactBilateralRow>,
    residuals: &mut Vec<f64>,
    rates: &mut Vec<f64>,
) -> Result<(), ArticulatedContactError> {
    if same_link(a, b) {
        return Err(ArticulatedContactError::Invalid(
            "point constraint requires distinct rigid links",
        ));
    }
    let (_, a_world) = prepare_point(a, Vec3::new(1., 0., 0.), kinematics)?;
    let (_, b_world) = prepare_point(b, Vec3::new(1., 0., 0.), kinematics)?;
    let pose_a = kinematics[a.system_index].world_from_link[a.link_index];
    let delta_world = b_world - a_world;
    let delta_a = pose_a.rotation().inverse().rotate(delta_world)?;
    if delta_a.norm() > cfg.position_tolerance_m {
        return Err(ArticulatedContactError::Invalid(
            "initial point closure is inconsistent",
        ));
    }
    let basis = [
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        Vec3::new(0., 0., 1.),
    ];
    for (component, axis_a) in basis.into_iter().enumerate() {
        let axis_world = pose_a.rotation().rotate(axis_a)?;
        let (first, world_point) = prepare_point(b, axis_world, kinematics)?;
        let (mut second, _) = prepare_point(a, axis_world.scale(-1.), kinematics)?;
        // g = R_A^T(x_B-x_A). Its A-side covector contains the geometric
        // couple -(x_B-x_A) x F, which makes the internal row objective and
        // preserves zero net spatial impulse away from exact closure.
        let couple_world = delta_world.cross(axis_world).scale(-1.);
        second.direction.torque =
            second.direction.torque + pose_a.rotation().inverse().rotate(couple_world)?;
        let row = Prepared {
            first,
            second: Some(second),
            axis_world,
            world_point,
            plane_speed: 0.,
            minimum_speed: None,
        };
        let rate = finite(row.velocity(kinematics))?;
        if rate.abs() > cfg.velocity_tolerance_mps {
            return Err(ArticulatedContactError::Invalid(
                "initial point-closure velocity is not tangent",
            ));
        }
        prepared.push(row);
        rows.push(ImpactBilateralRow {
            constraint: constraint_index,
            kind: if weld {
                ImpactBilateralRowKind::WeldPoint { component }
            } else {
                ImpactBilateralRowKind::Point { component }
            },
        });
        residuals.push(component_value(delta_a, component));
        rates.push(rate);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn append_weld_angular_rows(
    a: LinkContactPoint,
    b: LinkContactPoint,
    constraint_index: usize,
    cfg: ConstrainedImpactConfig,
    kinematics: &[crate::articulated::Kinematics],
    prepared: &mut Vec<Prepared>,
    rows: &mut Vec<ImpactBilateralRow>,
    residuals: &mut Vec<f64>,
    rates: &mut Vec<f64>,
) -> Result<(), ArticulatedContactError> {
    let pose_a = kinematics
        .get(a.system_index)
        .and_then(|k| k.world_from_link.get(a.link_index))
        .copied()
        .ok_or(ArticulatedContactError::Invalid(
            "weld A endpoint is out of range",
        ))?;
    let pose_b = kinematics
        .get(b.system_index)
        .and_then(|k| k.world_from_link.get(b.link_index))
        .copied()
        .ok_or(ArticulatedContactError::Invalid(
            "weld B endpoint is out of range",
        ))?;
    let basis = [
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        Vec3::new(0., 0., 1.),
    ];
    for (component, axis_a) in basis.into_iter().enumerate() {
        let axis_world = pose_a.rotation().rotate(axis_a)?;
        let torque_a = pose_a.rotation().inverse().rotate(axis_world.scale(-1.))?;
        let torque_b = pose_b.rotation().inverse().rotate(axis_world)?;
        let row = Prepared {
            first: Term {
                point: b,
                direction: Wrench::new(torque_b, Vec3::default()),
            },
            second: Some(Term {
                point: a,
                direction: Wrench::new(torque_a, Vec3::default()),
            }),
            axis_world,
            world_point: pose_b.transform_point(b.point_body)?,
            plane_speed: 0.,
            minimum_speed: None,
        };
        let rate = finite(row.velocity(kinematics))?;
        if rate.abs() > cfg.velocity_tolerance_mps {
            return Err(ArticulatedContactError::Invalid(
                "initial weld angular velocity is not tangent",
            ));
        }
        prepared.push(row);
        rows.push(ImpactBilateralRow {
            constraint: constraint_index,
            kind: ImpactBilateralRowKind::WeldAngular { component },
        });
        residuals.push(0.);
        rates.push(rate);
    }
    Ok(())
}

fn same_link(a: LinkContactPoint, b: LinkContactPoint) -> bool {
    a.system_index == b.system_index && a.link_index == b.link_index
}

fn component_value(v: Vec3, component: usize) -> f64 {
    match component {
        0 => v.x,
        1 => v.y,
        _ => v.z,
    }
}
