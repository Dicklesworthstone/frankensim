//! Simultaneous Coulomb point impacts with complete dry and attached-fluid inertia.
//!
//! Native geometry supplies the contact frames; the original non-associated
//! contact law solves normal and tangent impulses together. The common impact
//! mechanics workspace applies one accumulated impulse per articulated tree.

use super::mechanics::{ImpactWorkspace, PointImpactRow};
use super::{ArticulatedImpactSystem, ImpactMechanics};
use crate::Refusal;
use crate::addedmass::AeroGeneralizedLoads;
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactDeterminacy,
    FrictionalContactError, FrictionalContactResult, solve_frictional_contacts,
};
use fs_mbd::articulated::contact_geometry::{ContactCounterpart, prepare_contact_frame};
use fs_mbd::articulated::generalized::GeneralizedKinematics;
pub use fs_mbd::articulated_contacts::friction::ArticulatedFrictionalContact;

/// A fully audited event, retaining the original contact-law determinacy result.
#[derive(Debug, Clone)]
pub struct ArticulatedFrictionalImpact {
    /// Stable ascending identity order for all contact-space vectors and rows.
    pub contact_ids: Vec<u64>,
    /// Original Coulomb result. A velocity witness does not identify local loads.
    pub contact: FrictionalContactResult,
    /// Original complete `J M_total^-1 J^T`, normal-first within each contact.
    pub delassus: Vec<f64>,
    /// Accepted states, measured row speeds, total kinetic work and spatial audits.
    pub mechanics: ImpactMechanics,
}

fn refusal(code: &'static str, message: impl ToString) -> Refusal {
    Refusal {
        code,
        message: message.to_string(),
        ranked_repairs: vec![
            "inspect the original contact geometry, inertia and Coulomb solve budgets".into(),
        ],
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(refusal(
            "articulated-friction-cancelled",
            "impact cancelled before publication",
        ))
    } else {
        Ok(())
    }
}

/// Resolve one caller-located normal/tangent event on independent articulated trees.
///
/// The fluid callback supplies each tree's inertia at its original incoming
/// state. Cross-tree hydrodynamic added-mass blocks are not inferred. The same
/// prepared physical mass factor serves all normal/tangent response columns
/// and the accumulated update; no impedance diagonal, continuous force or
/// inertial bias acts during an instantaneous impact.
///
/// Contact identities canonicalize assembly. The original semismooth Coulomb
/// owner retains full cross-contact mobility and distinguishes unique impulses
/// from determined velocity with a redundant impulse gauge. The latter gives
/// no local load or heat allocation. Independently translating active supports
/// refuse when their work is not established as gauge-independent.
///
/// Native point speeds, rate limits, combined dry/fluid momentum and total
/// kinetic work must all admit the actual returned states. Pose, configuration
/// and memory are unchanged. Detection, time stepping, constitutive coefficient
/// selection, finite-patch spin and heat/deformation partition remain separate.
pub fn resolve(
    systems: &[ArticulatedImpactSystem<'_>],
    contacts: &[ArticulatedFrictionalContact],
    config: FrictionalContactConfig,
    evaluate_fluid: impl FnMut(usize, &GeneralizedKinematics) -> Result<AeroGeneralizedLoads, Refusal>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedFrictionalImpact, Refusal> {
    poll(&mut cancelled)?;
    if contacts.len() > 21 || contacts.len() > config.maximum_contacts {
        return Err(refusal(
            "articulated-friction-invalid",
            "contact event exceeds its size cap",
        ));
    }
    let mut order: Vec<_> = (0..contacts.len()).collect();
    order.sort_unstable_by_key(|&index| contacts[index].id);
    if order
        .windows(2)
        .any(|pair| contacts[pair[0]].id == contacts[pair[1]].id)
    {
        return Err(refusal(
            "articulated-friction-invalid",
            "duplicate contact identity",
        ));
    }
    let workspace = ImpactWorkspace::new(systems, evaluate_fluid, &mut cancelled)?;
    let kinematics = workspace.kinematics();
    let mut rows = Vec::with_capacity(3 * contacts.len());
    let mut laws = Vec::with_capacity(contacts.len());
    for &index in &order {
        poll(&mut cancelled)?;
        let c = contacts[index];
        let geometry = prepare_contact_frame(
            c.point,
            c.other,
            c.normal_world,
            c.tangent_world,
            index,
            &kinematics,
        )
        .map_err(|error| refusal("articulated-friction-invalid", error))?;
        laws.push(FrictionalContact {
            id: c.id,
            gap: geometry.rows[0].gap,
            activation_tolerance: c.activation_tolerance,
            relative_velocity: geometry.rows.map(|row| row.relative_normal_speed),
            friction_coefficient: c.friction_coefficient,
            normal_law: c.normal_law,
        });
        rows.extend(geometry.rows.map(PointImpactRow::from));
    }
    let prepared = workspace.project_rows(&rows, &mut cancelled)?;
    let contact = solve_frictional_contacts(&laws, &prepared.delassus, config, &mut cancelled)
        .map_err(|error| {
            let code = if matches!(error, FrictionalContactError::Cancelled) {
                "articulated-friction-cancelled"
            } else {
                "articulated-friction-law"
            };
            refusal(code, error)
        })?;
    if matches!(
        contact.determinacy,
        FrictionalContactDeterminacy::UniqueVelocity(_)
    ) {
        // A null contact impulse changes no body's velocity, hence no total
        // external impulse. One common support translation therefore has
        // fixed total work; unequal support velocities need a further witness.
        let mut common_velocity = None;
        for &active in &contact.active_indices {
            if let ContactCounterpart::PrescribedPlane { velocity_world, .. } =
                contacts[order[active]].other
            {
                if common_velocity.is_some_and(|v| v != velocity_world) {
                    return Err(refusal(
                        "articulated-friction-support-work",
                        "redundant impulses do not establish work for independently moving supports",
                    ));
                }
                common_velocity = Some(velocity_world);
            }
        }
    }
    let impulses: Vec<_> = contact.impulses.iter().flatten().copied().collect();
    let predicted: Vec<_> = contact
        .relative_velocities_after
        .iter()
        .flatten()
        .copied()
        .collect();
    let mechanics = workspace.apply(
        &prepared,
        &impulses,
        &predicted,
        contact.dissipated_energy,
        config.velocity_tolerance,
        config.energy_tolerance,
        &mut cancelled,
    )?;
    poll(&mut cancelled)?;
    Ok(ArticulatedFrictionalImpact {
        contact_ids: order.iter().map(|&index| contacts[index].id).collect(),
        contact,
        delassus: prepared.delassus,
        mechanics,
    })
}
