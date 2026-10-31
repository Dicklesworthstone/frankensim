//! Simultaneous normal impulses through the complete dry and fluid inertia.
//!
//! Native material-point geometry supplies the rows. The existing global
//! `fs-contact` law chooses impulses from one full Delassus matrix, and the
//! prepared total mass applies their sum once. Force/bias terms are not
//! integrated during an instantaneous event. The caller owns contact selection
//! and event time; this API never advances or repairs configuration.

#[path = "impact/mechanics.rs"]
pub(crate) mod mechanics;
pub use mechanics::{ArticulatedImpactSystem, ImpactMechanics};

/// Coupled Coulomb contacts through the same physical dry/fluid mechanics.
#[path = "impact/friction.rs"]
pub mod friction;

use super::ArticulatedState;
use crate::Refusal;
use crate::addedmass::{AeroGeneralizedLoads, PreparedArticulatedMass};
use fs_contact::multicontact::{
    MultiImpactConfig, MultiImpactError, MultiImpactResult, NormalContactImpact,
    solve_normal_impacts,
};
use fs_ga::Vec3;
use fs_mbd::articulated::ArticulatedModel;
use fs_mbd::articulated::contact_geometry::prepare_normal_contact;
pub use fs_mbd::articulated::contact_geometry::{
    ArticulatedNormalContact, ContactCounterpart, LinkContactPoint,
};
use fs_mbd::articulated::generalized::{GeneralizedError, GeneralizedKinematics};
use mechanics::{ImpactWorkspace, PointImpactRow};

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

/// Simultaneous normal contact on independent fluid-bearing articulated trees.
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedImpacts {
    /// All accepted states and complete global physical endpoint balances.
    pub mechanics: ImpactMechanics,
    /// Existing globally passive normal law in original contact order.
    pub law: MultiImpactResult,
    /// Full row-major contact mobility, including shared-tree cross contacts.
    pub delassus: Vec<f64>,
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
/// refuse atomically. [`resolve_many`] exposes the same mechanics for multiple
/// independent trees; this entry point retains its one-tree public API.
pub fn resolve(
    model: &ArticulatedModel,
    state: &ArticulatedState,
    contacts: &[ArticulatedNormalContact],
    config: MultiImpactConfig,
    evaluate_fluid: impl FnOnce(&GeneralizedKinematics) -> Result<AeroGeneralizedLoads, Refusal>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedImpact, Refusal> {
    let mut callback = Some(evaluate_fluid);
    let mut event = resolve_many(
        &[ArticulatedImpactSystem { model, state }],
        contacts,
        config,
        |_, map| {
            callback.take().ok_or_else(|| {
                refusal(
                    "articulated-impact-invalid",
                    "single-tree fluid callback repeated",
                )
            })?(map)
        },
        &mut cancelled,
    )?;
    let state = event.mechanics.states.pop().ok_or_else(|| {
        refusal(
            "articulated-impact-invalid",
            "missing single-tree impact state",
        )
    })?;
    let velocity_after = event.mechanics.velocities_after.pop().ok_or_else(|| {
        refusal(
            "articulated-impact-invalid",
            "missing single-tree impact velocity",
        )
    })?;
    let generalized = GeneralizedImpactResponse {
        velocity_after,
        law: event.law,
        delassus: event.delassus,
        kinetic_before: event.mechanics.kinetic_before,
        kinetic_after: event.mechanics.kinetic_after,
        kinetic_change: event.mechanics.kinetic_change,
        support_work: event.mechanics.support_work,
        energy_balance_residual: event.mechanics.energy_balance_residual,
    };
    Ok(ArticulatedImpact {
        state,
        generalized,
        relative_normal_speeds_after: event.mechanics.relative_speeds_after,
        linear_momentum_residual: event.mechanics.linear_momentum_residual,
        angular_momentum_residual: event.mechanics.angular_momentum_residual,
    })
}

/// Resolve a common normal-impact event across up to 128 articulated trees.
///
/// Each tree has its own complete dry and attached-fluid inertia, prepared
/// once from its native map. Contacts can join different trees, act internally
/// within one tree, or meet prescribed moving planes. Only shared incident
/// trees contribute to a contact cross block; no global dense body mass is
/// assembled. Impulses are solved together and accumulated before each tree's
/// one velocity reset. States and aerodynamic memories remain atomic on refusal.
///
/// The fluid callback is invoked once per supplied tree in original order.
/// Cross-tree hydrodynamic added-mass terms are not represented: this is the
/// independent-fluid assembly used by the attached-strip model. Friction, CCD
/// and persistent contact belong to their separate law/step consumers.
pub fn resolve_many(
    systems: &[ArticulatedImpactSystem<'_>],
    contacts: &[ArticulatedNormalContact],
    config: MultiImpactConfig,
    evaluate_fluid: impl FnMut(usize, &GeneralizedKinematics) -> Result<AeroGeneralizedLoads, Refusal>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedImpacts, Refusal> {
    poll(&mut cancelled)?;
    if contacts.len() > 64 || contacts.len() > config.solver.max_variables {
        return Err(refusal(
            "articulated-impact-invalid",
            "articulated contact event exceeds its row cap",
        ));
    }
    let workspace = ImpactWorkspace::new(systems, evaluate_fluid, &mut cancelled)?;
    let maps = workspace.kinematics();
    let mut geometry = Vec::with_capacity(contacts.len());
    for (index, contact) in contacts.iter().enumerate() {
        poll(&mut cancelled)?;
        geometry.push(
            prepare_normal_contact(contact, index, &maps)
                .map_err(|e| refusal("articulated-impact-invalid", e))?,
        );
    }
    let directions: Vec<_> = geometry.iter().copied().map(PointImpactRow::from).collect();
    let prepared = workspace.project_rows(&directions, &mut cancelled)?;
    let laws: Vec<_> = contacts
        .iter()
        .zip(&geometry)
        .enumerate()
        .map(|(i, (contact, g))| NormalContactImpact {
            gap: g.gap,
            activation_tolerance: contact.activation_tolerance,
            relative_normal_speed: prepared.relative_speeds_before[i],
            restitution: contact.restitution,
        })
        .collect();
    let law =
        solve_normal_impacts(&laws, &prepared.delassus, config, &mut cancelled).map_err(|e| {
            refusal(
                if e == MultiImpactError::Cancelled {
                    "articulated-impact-cancelled"
                } else {
                    "articulated-impact-law"
                },
                e,
            )
        })?;
    let mechanics = workspace.apply(
        &prepared,
        &law.normal_impulses,
        &law.relative_normal_speeds_after,
        law.dissipated_energy,
        config.solver.gradient_tolerance,
        config.energy_tolerance,
        &mut cancelled,
    )?;
    Ok(ArticulatedImpacts {
        mechanics,
        law,
        delassus: prepared.delassus,
    })
}
