//! Original local-speed Hertz curves under drive and temperature feedback.
//!
//! Every trial refreshes the full selected parameter bundle, then advances
//! the captured original body pair and accepted mode. Physical pressure/slip
//! work heats the existing finite cells; motor work remains a separate supply.

use super::*;
use crate::friction_heat::hertz::driven::HeatedRigidHertzDrivenMarch;
use fs_mbd::Wrench;
use fs_mbd::tangential_friction::hertz::rate_dependent::driven::{
    PreparedRigidHertzRateDriven, RigidHertzRateDrivenError, RigidHertzRateDrivenReport,
};
use fs_tribo::sliding_patch::contact::driven::HertzDrivenMode;
use fs_tribo::sliding_patch::rate_dependent::driven::HertzRateDrivenControls;

pub type RigidHertzRateDrivenHeatError<E> =
    RigidHertzHeatError<E, RigidHertzRateDrivenError, SourceHertzFrictionError>;
pub type HeatedRigidHertzRateDrivenCandidate<'m, P> =
    HeatedRigidHertzCandidate<'m, P, RigidHertzRateDrivenReport, SourceHertzFriction>;
pub type HeatedRigidHertzRateDrivenMarch<'m, P> =
    HeatedRigidHertzDrivenMarch<'m, P, SourceHertzFriction>;

/// Independent finite-cell, original pressure-quadrature and driven-time gates.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzRateDrivenHeatControls {
    pub heat: RigidHertzHeatControls,
    pub quadrature: HertzFrictionQuadrature,
    pub driven: HertzRateDrivenControls,
}

/// Couple one driven original local-speed phase to its finite-cell heat.
///
/// `external` follows contact A/B order, with world forces and body-frame
/// torques about the original centres of mass. The same source claims resolve
/// at each common-face trial temperature; their original curve is evaluated
/// across the spatial patch and throughout physical motion. No single effective
/// slip speed or fitted coefficient replaces that curve. Retained mode belongs
/// to the original mechanical trajectory and restarts unchanged for each trial.
///
/// Poses, inertia, normal load, radius and thermal footprint remain held; the
/// native owner admits only its closed-normal domain. Spatial pressure, native
/// motion and thermal windows have independent refinement controls. Source
/// parameter support alone does not validate a chosen curve over unmeasured
/// speeds, and a converged coupled root does not establish its uniqueness.
/// Admitted invariant affine paths retain their original exact arrest,
/// subsequent static hold or reversal. Their native event clock/bracket stays
/// attached to the candidate. Other unresolved rate-law events still refuse;
/// genuine Coulomb retains its established driven event/time semantics.
#[allow(clippy::too_many_arguments)]
pub fn propose_rigid_hertz_rate_driven_heat<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzRateHeatContact<'_>,
    external: [Wrench; 2],
    mode: HertzDrivenMode,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidHertzRateDrivenHeatControls,
    energy_guess_j: &[f64],
    temperature_guess_k: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzRateDrivenCandidate<'m, P>, RigidHertzRateDrivenHeatError<E>> {
    solve_rigid_hertz_heat(
        contact,
        thermal,
        route,
        accepted_time_s,
        dt_s,
        controls.heat,
        energy_guess_j,
        temperature_guess_k,
        |temperature| {
            contact.source.hertz_at(
                temperature,
                contact.nominal_pressure_pa,
                contact.normal_force_n,
                contact.contact_radius_m,
            )
        },
        |source, duration| {
            let prepared = PreparedRigidHertzRateDriven::new(
                bodies,
                properties,
                RigidHertzRateContact {
                    binding: contact.binding,
                    frame: contact.frame,
                    law: source.law(),
                },
                external,
                mode,
                duration,
                controls.heat.mechanical,
                controls.quadrature,
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
                dissipated_j: candidate.report.mechanical.dissipated_energy_j,
                total_energy_residual_j: candidate.report.mechanical.audit.total_energy_residual_j,
                report: candidate.report,
            })
        },
        resolve,
        gate,
    )
}

/// Resolve rate/temperature feedback, accepted mode, heat and drive work in time.
///
/// The existing Hertz full/two-half controller retains only original physical
/// fine states. It independently compares each body and cell, travel, friction
/// heat, external work, endpoint mode and whether a driven arrest occurred.
/// A later source, mechanical, thermal,
/// numerical, budget or cancellation refusal returns no accepted mutation.
/// Local estimates do not certify a global bound; each inner discretization
/// remains an independent refinement axis. Available inner travel/heat estimates
/// also enter the comparison, including their contribution to native drive
/// work. Analytic invariant paths declare their exact formula separately from
/// adaptive estimates. Delegated curved Coulomb motion requires independent
/// held-traction step refinement. Source variation still needs thermal-window
/// refinement; a native exact event time is conditional on its held source.
/// The phase's physical scope applies.
#[allow(clippy::too_many_arguments)]
pub fn propose_rigid_hertz_rate_driven_heat_march<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzRateHeatContact<'_>,
    external: [Wrench; 2],
    mode: HertzDrivenMode,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    duration_s: f64,
    phase_controls: RigidHertzRateDrivenHeatControls,
    controls: march::HertzHeatMarchControls,
    temperature_guess_k: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzRateDrivenMarch<'m, P>, march::HertzHeatMarchError<E>> {
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
            let candidate = propose_rigid_hertz_rate_driven_heat(
                b,
                props,
                pair,
                external,
                *state,
                h,
                route,
                start,
                dt,
                RigidHertzRateDrivenHeatControls {
                    heat,
                    ..phase_controls
                },
                energy,
                temperature,
                &resolve,
                gate,
            )
            .map_err(march::HertzHeatMarchError::RateDriven)?;
            let mode = candidate.mechanical.mechanical.driven.mode;
            Ok(march::stage(
                candidate,
                start,
                march::HertzHeatMarchMotion::RateDriven,
                mode,
            ))
        },
        gate,
    )?;
    Ok(HeatedRigidHertzDrivenMarch { march, mode_after })
}
