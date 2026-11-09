//! Continuous force/spin drive, original Hertz mode and finite-cell heat.
//!
//! Explicit external work is a supply, including work on common body motion
//! during static contact. Only actual relative friction loss heats the cells.
//! Every source/thermal trial starts from the same accepted bodies and mode.

use super::*;
use fs_mbd::Wrench;
use fs_mbd::tangential_friction::hertz::driven::{
    PreparedRigidHertzDriven, RigidHertzDrivenError, RigidHertzDrivenReport,
};
use fs_tribo::sliding_patch::contact::driven::{HertzDrivenControls, HertzDrivenMode};

pub type RigidHertzDrivenHeatError<E> = RigidHertzHeatError<E, RigidHertzDrivenError>;
pub type HeatedRigidHertzDrivenCandidate<'m, P> =
    HeatedRigidHertzCandidate<'m, P, RigidHertzDrivenReport>;

/// Original finite-cell/source gates and independent driven trajectory controls.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzDrivenHeatControls {
    pub heat: RigidHertzHeatControls,
    pub driven: HertzDrivenControls,
}

/// Publish the actual accepted mode, native pair and thermal state together.
/// A resting mode is native mechanical state, not a velocity cutoff.
pub struct HeatedRigidHertzDrivenMarch<'m, P> {
    pub march: march::HeatedRigidHertzMarchCandidate<'m, P>,
    pub mode_after: HertzDrivenMode,
}

/// Propose a complete driven Hertz/source/thermal phase without mutating input.
///
/// `external` follows contact A/B order: forces are world vectors and torques
/// are body vectors about the original centres of mass. They act continuously;
/// no external pre-kick chooses the friction branch. Each thermal trial resolves
/// the same original static/kinetic claims at its common-face temperature, then
/// advances the captured mode and body response under that held source snapshot.
/// The returned native report owns the accepted mode and actual drive work.
///
/// Poses, inertia, normal load, radius and thermal footprint remain fixed, with
/// held closed normal motion admitted by the native owner. Spatial cells,
/// source/thermal windows and the driven curved maximum step need independent
/// refinement. No extra heat store, temperature-driven mode guess, rate-weakening
/// branch selection, or dissipation inferred from supplied drive work is added.
#[allow(clippy::too_many_arguments)]
pub fn propose_rigid_hertz_driven_heat<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzHeatContact<'_>,
    external: [Wrench; 2],
    mode: HertzDrivenMode,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidHertzDrivenHeatControls,
    energy_guess_j: &[f64],
    temperature_guess_k: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzDrivenCandidate<'m, P>, RigidHertzDrivenHeatError<E>> {
    solve_rigid_hertz_heat(
        contact,
        thermal,
        route,
        accepted_time_s,
        dt_s,
        controls.heat,
        energy_guess_j,
        temperature_guess_k,
        |temperature| refresh_source(contact, temperature),
        |source, duration| {
            let prepared = PreparedRigidHertzDriven::new(
                bodies,
                properties,
                RigidHertzContact {
                    binding: contact.binding,
                    frame: contact.frame,
                    law: source.law(),
                },
                external,
                mode,
                duration,
                controls.heat.mechanical,
                controls.driven,
                gate,
            )?;
            let point = prepared.contact_points_world()[0];
            Ok((prepared, [point.x, point.y, point.z]))
        },
        |prepared, source| {
            let candidate = prepared.advance(source.law(), gate)?;
            Ok(HertzHeatMechanics {
                bodies_after: candidate.bodies_after,
                dissipated_j: candidate.report.dissipated_energy_j,
                total_energy_residual_j: candidate.report.audit.total_energy_residual_j,
                report: candidate.report,
            })
        },
        resolve,
        gate,
    )
}

/// Resolve evolving source temperature and drive/heat work over a whole span.
///
/// Coarse and fine comparisons restart the same actual accepted mode. Only
/// the two physical half phases advance it, together with the original momenta
/// and finite caloric state. A later failure publishes none of the earlier
/// private phases. Drive work has its own coarse/fine energy comparison and
/// remains separate from friction heat in the compensated whole-span balance.
///
/// Source temperature is sampled per phase; thermal threshold-crossing timing
/// requires outer refinement. Curved driven motion has independent native time
/// discretization: outer comparisons need not reveal a shared inner-step error.
/// The local estimates do not certify a global error or a unique discontinuous
/// branch. The normal/geometry/footprint limits of the phase API remain in force.
#[allow(clippy::too_many_arguments)]
pub fn propose_rigid_hertz_driven_heat_march<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzHeatContact<'_>,
    external: [Wrench; 2],
    mode: HertzDrivenMode,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    duration_s: f64,
    phase_controls: RigidHertzDrivenHeatControls,
    controls: march::HertzHeatMarchControls,
    temperature_guess_k: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzDrivenMarch<'m, P>, march::HertzHeatMarchError<E>> {
    let (march, mode_after) = march::propose_hertz_heat_phases(
        bodies,
        properties,
        contact,
        thermal,
        mode,
        duration_s,
        phase_controls.heat,
        controls,
        temperature_guess_k,
        |b, props, pair, h, state, dt, heat, energy, temperature| {
            let start = h.time_s();
            let candidate = propose_rigid_hertz_driven_heat(
                b,
                props,
                pair,
                external,
                *state,
                h,
                route,
                start,
                dt,
                RigidHertzDrivenHeatControls {
                    heat,
                    ..phase_controls
                },
                energy,
                temperature,
                &resolve,
                gate,
            )
            .map_err(march::HertzHeatMarchError::Driven)?;
            let next = candidate.mechanical.driven.mode;
            Ok(march::stage(
                candidate,
                start,
                march::HertzHeatMarchMotion::Driven,
                next,
            ))
        },
        gate,
    )?;
    Ok(HeatedRigidHertzDrivenMarch { march, mode_after })
}
