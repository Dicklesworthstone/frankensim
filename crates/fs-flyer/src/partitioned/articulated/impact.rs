//! Simultaneous normal impulses through the complete dry and fluid inertia.
//!
//! Native material-point geometry supplies the rows. The existing global
//! `fs-contact` law chooses impulses from one full Delassus matrix, and the
//! prepared total mass applies their sum once. Force/bias terms are not
//! integrated during an instantaneous event. The caller owns contact selection
//! and event time; this API never advances or repairs configuration.

use super::{ArticulatedState, MAX_ARTICULATED_LINKS};
use crate::Refusal;
use crate::addedmass::joint::native_index;
use crate::addedmass::{AeroGeneralizedLoads, MAX_CONTROL_COORDS, PreparedArticulatedMass};
use fs_contact::multicontact::{
    MultiImpactConfig, MultiImpactError, MultiImpactResult, NormalContactImpact,
    solve_normal_impacts,
};
use fs_ga::{Twist, Vec3, Wrench};
use fs_mbd::articulated::contact_geometry::prepare_normal_contact;
pub use fs_mbd::articulated::contact_geometry::{
    ArticulatedNormalContact, ContactCounterpart, LinkContactPoint,
};
use fs_mbd::articulated::generalized::{
    GeneralizedError, GeneralizedKinematics, free_floating_terms,
};
use fs_mbd::articulated::{ArticulatedModel, BaseState, forward_kinematics};

/// Positive separating row in Flyer `[linear, angular, joints]` coordinates.
/// Gap/rate use metres and m/s for points, radians and rad/s for revolute stops;
/// conjugate impulses are respectively N s and N m s. Their pairing is work.
pub(crate) struct GeneralizedNormalRow {
    pub jacobian: Vec<f64>,
    pub gap: f64,
    pub activation_tolerance: f64,
    pub support_speed: f64,
    pub restitution: f64,
}

/// One admitted generalized velocity reset at unchanged configuration.
#[derive(Clone, Debug, PartialEq)]
pub struct GeneralizedImpactResponse {
    /// Updated velocity in Flyer coordinate order.
    pub velocity_after: Vec<f64>,
    /// Original global contact-law result, in caller row order.
    pub law: MultiImpactResult,
    /// Full row-major contact mobility `J M_total^-1 J^T`.
    pub delassus: Vec<f64>,
    /// Complete dry plus fluid kinetic energy before impact [J].
    pub kinetic_before: f64,
    /// Complete dry plus fluid kinetic energy after impact [J].
    pub kinetic_after: f64,
    /// Endpoint energy change, evaluated by stable quadratic polarization [J].
    pub kinetic_change: f64,
    /// Signed prescribed-support work [J].
    pub support_work: f64,
    /// Kinetic change minus support work plus mechanical loss [J].
    pub energy_balance_residual: f64,
}

/// Accepted point-contact event on one native articulated tree.
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedImpact {
    /// Same configuration/memory with all base and joint rates updated together.
    pub state: ArticulatedState,
    /// Shared total-mass reset, original contact law and work accounting.
    pub generalized: GeneralizedImpactResponse,
    /// Relative normal speeds measured by original FK on the returned state.
    pub relative_normal_speeds_after: Vec<f64>,
    /// Combined dry/fluid linear momentum defect after external plane impulse.
    pub linear_momentum_residual: Vec3,
    /// Combined world-origin angular momentum defect after plane impulse.
    pub angular_momentum_residual: Vec3,
}

fn refusal(code: &'static str, message: impl ToString) -> Refusal {
    Refusal {
        code,
        message: message.to_string(),
        ranked_repairs: vec![
            "inspect the event geometry, total inertia and original contact-law budgets".into(),
        ],
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(refusal(
            "articulated-impact-cancelled",
            "impact cancelled before publication",
        ))
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64, Refusal> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(refusal(
            "articulated-impact-invalid",
            "nonfinite impulse arithmetic",
        ))
    }
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn vec(v: &[f64]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn native(error: GeneralizedError) -> Refusal {
    if error == GeneralizedError::Cancelled {
        refusal("articulated-impact-cancelled", error)
    } else {
        refusal("articulated-impact-invalid", error)
    }
}

/// Shared reset bridge for points and generalized travel stops. Numerical speed
/// budgets use row-coordinate units; energy budgets are joules. Geometry/event
/// localization and final native rate admission remain caller-owned.
pub(crate) fn resolve_normal_rows(
    mass: &PreparedArticulatedMass,
    velocity_before: &[f64],
    rows: &[GeneralizedNormalRow],
    config: MultiImpactConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GeneralizedImpactResponse, Refusal> {
    poll(&mut cancelled)?;
    let n = mass.coordinate_count();
    if rows.len() > 64
        || rows.len() > config.solver.max_variables
        || velocity_before.len() != n
        || !velocity_before.iter().all(|x| x.is_finite())
        || rows.iter().any(|row| {
            row.jacobian.len() != n
                || !row.jacobian.iter().all(|x| x.is_finite())
                || !row.support_speed.is_finite()
        })
    {
        return Err(refusal(
            "articulated-impact-invalid",
            "invalid generalized contact rows or event size",
        ));
    }
    let count = rows.len();
    let mut columns = Vec::with_capacity(count);
    let mut laws = Vec::with_capacity(count);
    for row in rows {
        poll(&mut cancelled)?;
        columns.push(mass.solve(&row.jacobian)?);
        laws.push(NormalContactImpact {
            gap: row.gap,
            activation_tolerance: row.activation_tolerance,
            relative_normal_speed: finite(dot(&row.jacobian, velocity_before) - row.support_speed)?,
            restitution: row.restitution,
        });
    }
    let mut delassus = vec![0.0; count * count];
    for i in 0..count {
        poll(&mut cancelled)?;
        for j in 0..=i {
            // The single admitted mass factor defines a symmetric bilinear
            // form; average only the two evaluation orders' roundoff.
            let value = finite(
                0.5 * dot(&rows[i].jacobian, &columns[j])
                    + 0.5 * dot(&rows[j].jacobian, &columns[i]),
            )?;
            delassus[i * count + j] = value;
            delassus[j * count + i] = value;
        }
    }
    let law = solve_normal_impacts(&laws, &delassus, config, &mut cancelled).map_err(|error| {
        let code = if error == MultiImpactError::Cancelled {
            "articulated-impact-cancelled"
        } else {
            "articulated-impact-law"
        };
        refusal(code, error)
    })?;
    let mut impulse = vec![0.0; n];
    let mut support_work = 0.0;
    for (row, lambda) in rows.iter().zip(&law.normal_impulses) {
        poll(&mut cancelled)?;
        for (entry, j) in impulse.iter_mut().zip(&row.jacobian) {
            *entry = finite(*entry + lambda * j)?;
        }
        support_work = finite(support_work + lambda * row.support_speed)?;
    }
    let change = mass.solve(&impulse)?;
    let velocity_after = velocity_before
        .iter()
        .zip(change)
        .map(|(u, du)| finite(u + du))
        .collect::<Result<Vec<_>, _>>()?;
    for (i, row) in rows.iter().enumerate() {
        let actual = finite(dot(&row.jacobian, &velocity_after) - row.support_speed)?;
        let scale: f64 = row
            .jacobian
            .iter()
            .zip(velocity_before.iter().zip(&velocity_after))
            .map(|(j, (a, b))| j.abs() * (a.abs() + b.abs()))
            .sum();
        let allowed = finite(
            config.solver.gradient_tolerance
                + 256.0 * f64::EPSILON * (scale + row.support_speed.abs()),
        )?;
        if (actual - law.relative_normal_speeds_after[i]).abs() > allowed {
            return Err(refusal(
                "articulated-impact-unresolved",
                "actual velocity reset does not resolve contact-law speeds",
            ));
        }
    }
    let kinetic_before = mass.kinetic_energy(velocity_before)?;
    let kinetic_after = mass.kinetic_energy(&velocity_after)?;
    let midpoint: Vec<_> = velocity_before
        .iter()
        .zip(&velocity_after)
        .map(|(a, b)| f64::midpoint(*a, *b))
        .collect();
    let momentum_midpoint = mass.momentum(&midpoint)?;
    let delta: Vec<_> = velocity_after
        .iter()
        .zip(velocity_before)
        .map(|(b, a)| b - a)
        .collect();
    let kinetic_change = finite(dot(&delta, &momentum_midpoint))?;
    let scale: f64 = delta
        .iter()
        .zip(&momentum_midpoint)
        .map(|(du, p)| (du * p).abs())
        .sum();
    let energy_balance_residual = finite(kinetic_change - support_work + law.dissipated_energy)?;
    let allowed = finite(
        config.energy_tolerance
            + 512.0 * f64::EPSILON * (scale + support_work.abs() + law.dissipated_energy),
    )?;
    if energy_balance_residual.abs() > allowed {
        return Err(refusal(
            "articulated-impact-unresolved",
            "complete dry/fluid kinetic work does not close",
        ));
    }
    poll(&mut cancelled)?;
    Ok(GeneralizedImpactResponse {
        velocity_after,
        law,
        delassus,
        kinetic_before,
        kinetic_after,
        kinetic_change,
        support_work,
        energy_balance_residual,
    })
}

/// Resolve caller-located simultaneous normal contacts using total fluid inertia.
///
/// Contacts use original native definitions; all system indices must be zero.
/// Self contacts retain equal/opposite impulses on one line of action. Moving
/// planes supply explicit external spatial impulse and signed work. Poses,
/// joint coordinates and aero-memory are unchanged by the instantaneous reset.
///
/// The fluid callback uses the same native event map as dry inertia and must
/// have zero applied-load fields. Bias is not integrated. All Delassus columns
/// reuse one factor; impulses are accumulated before one rate update. Original
/// global restitution/passivity, native rate limits and endpoint balances can
/// refuse atomically. Friction, CCD, persistent contact and multiple independent
/// trees are outside this entry point.
pub fn resolve(
    model: &ArticulatedModel,
    state: &ArticulatedState,
    contacts: &[ArticulatedNormalContact],
    config: MultiImpactConfig,
    evaluate_fluid: impl FnOnce(&GeneralizedKinematics) -> Result<AeroGeneralizedLoads, Refusal>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedImpact, Refusal> {
    poll(&mut cancelled)?;
    if model.link_count() > MAX_ARTICULATED_LINKS
        || model.dof_count() > MAX_CONTROL_COORDS
        || contacts.len() > 64
        || contacts.len() > config.solver.max_variables
    {
        return Err(refusal(
            "articulated-impact-invalid",
            "articulated contact event exceeds its size cap",
        ));
    }
    let motion = &state.motion;
    let dry = free_floating_terms(
        model,
        motion.base,
        &motion.joint_position,
        &motion.joint_velocity,
        &mut cancelled,
    )
    .map_err(native)?;
    let fluid = evaluate_fluid(&dry.kinematics)?;
    poll(&mut cancelled)?;
    if fluid
        .q_rigid_nonaccel
        .iter()
        .chain(&fluid.q_control_nonaccel)
        .any(|x| *x != 0.0)
    {
        return Err(refusal(
            "coupled-load-ownership",
            "instantaneous fluid inertia must not contain applied force loads",
        ));
    }
    let mass = PreparedArticulatedMass::new(&dry, &fluid)?;
    let native_u = motion.base.twist_body.to_array();
    let mut velocity: Vec<_> = (0..6).map(|i| native_u[native_index(i)]).collect();
    velocity.extend(&motion.joint_velocity);
    let mut geometries = Vec::with_capacity(contacts.len());
    let mut rows = Vec::with_capacity(contacts.len());
    for (i, contact) in contacts.iter().enumerate() {
        poll(&mut cancelled)?;
        let geometry = prepare_normal_contact(
            contact,
            i,
            std::slice::from_ref(dry.kinematics.native_kinematics()),
        )
        .map_err(|error| refusal("articulated-impact-invalid", error))?;
        let mut jacobian = vec![0.0; mass.coordinate_count()];
        for term in core::iter::once(geometry.first).chain(geometry.second) {
            let columns = dry.kinematics.body_jacobian(term.point.link_index);
            for (i, entry) in jacobian.iter_mut().enumerate() {
                *entry += term.direction.pairing(columns[native_index(i)]);
            }
        }
        rows.push(GeneralizedNormalRow {
            jacobian,
            gap: geometry.gap,
            activation_tolerance: contact.activation_tolerance,
            support_speed: geometry.plane_speed,
            restitution: contact.restitution,
        });
        geometries.push(geometry);
    }
    let generalized = resolve_normal_rows(&mass, &velocity, &rows, config, &mut cancelled)?;
    let mut accepted = state.clone();
    let after = &generalized.velocity_after;
    accepted.motion.base.twist_body = Twist::new(vec(&after[3..6]), vec(&after[..3]));
    accepted.motion.joint_velocity.copy_from_slice(&after[6..]);
    let after_kinematics = forward_kinematics(
        model,
        BaseState::prescribed(
            accepted.motion.base.world_from_base,
            accepted.motion.base.twist_body,
            Twist::zero(),
        ),
        &accepted.motion.joint_position,
        &accepted.motion.joint_velocity,
    )
    .map_err(|error| refusal("articulated-impact-invalid", error))?;
    let mut speeds = Vec::with_capacity(contacts.len());
    let mut external_linear = Vec3::default();
    let mut external_angular = Vec3::default();
    let (mut linear_scale, mut angular_scale) = (0.0, 0.0);
    for (i, (geometry, lambda)) in geometries
        .iter()
        .zip(&generalized.law.normal_impulses)
        .enumerate()
    {
        poll(&mut cancelled)?;
        let speed = finite(
            core::iter::once(geometry.first)
                .chain(geometry.second)
                .map(|term| {
                    term.direction
                        .pairing(after_kinematics.body_twist[term.point.link_index])
                })
                .sum::<f64>()
                - geometry.plane_speed,
        )?;
        let allowed = finite(
            config.solver.gradient_tolerance
                + 256.0 * f64::EPSILON * (speed.abs() + geometry.relative_normal_speed.abs()),
        )?;
        if (speed - generalized.law.relative_normal_speeds_after[i]).abs() > allowed {
            return Err(refusal(
                "articulated-impact-unresolved",
                "native endpoint point speed does not match global contact law",
            ));
        }
        speeds.push(speed);
        for term in core::iter::once(geometry.first).chain(geometry.second) {
            let world = term
                .direction
                .transform_by(&dry.kinematics.world_from_link(term.point.link_index));
            linear_scale = finite(linear_scale + lambda.abs() * world.force.norm())?;
            angular_scale = finite(angular_scale + lambda.abs() * world.torque.norm())?;
        }
        if geometry.second.is_none() {
            let impulse = geometry.axis_world.scale(*lambda);
            external_linear = external_linear + impulse;
            external_angular = external_angular + geometry.world_point.cross(impulse);
        }
    }
    let delta: Vec<_> = after.iter().zip(&velocity).map(|(b, a)| b - a).collect();
    let dp = mass.momentum(&delta)?;
    let world_delta =
        Wrench::new(vec(&dp[3..6]), vec(&dp[..3])).transform_by(&motion.base.world_from_base);
    let linear_momentum_residual = world_delta.force - external_linear;
    let angular_momentum_residual = world_delta.torque - external_angular;
    if finite(linear_momentum_residual.norm())?
        > 2048.0 * f64::EPSILON * (linear_scale + world_delta.force.norm())
        || finite(angular_momentum_residual.norm())?
            > 2048.0 * f64::EPSILON * (angular_scale + world_delta.torque.norm())
    {
        return Err(refusal(
            "articulated-impact-unresolved",
            "complete dry/fluid spatial impulse does not close",
        ));
    }
    poll(&mut cancelled)?;
    Ok(ArticulatedImpact {
        state: accepted,
        generalized,
        relative_normal_speeds_after: speeds,
        linear_momentum_residual,
        angular_momentum_residual,
    })
}
