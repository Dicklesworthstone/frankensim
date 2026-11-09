//! Full original local-speed friction curves with finite-cell temperature feedback.
//!
//! Every thermal trial refreshes all selected curve parameters from the same
//! interface claims. Native Hertz quadrature evaluates that curve at the local
//! translational/spin speed; a centre-speed effective coefficient is never used.
//! Physical path loss enters the unchanged finite thermal network exactly once.

use super::*;
use fs_contact::interface_binding::thermal_friction::rate_dependent::{
    SourceHertzFriction, SourceHertzFrictionError, SourceThermalFriction,
};
use fs_mbd::tangential_friction::hertz::rate_dependent::{
    PreparedRigidHertzRateTransient, RigidHertzRateContact, RigidHertzRateTransientError,
    RigidHertzRateTransientReport,
};
use fs_mbd::tangential_friction::hertz::transient::HertzTransientControls;
use fs_tribo::sliding_patch::rate_dependent::HertzFrictionQuadrature;

pub type RigidHertzRateHeatContact<'a> = RigidHertzHeatContact<'a, SourceThermalFriction>;
pub type RigidHertzRateHeatError<E> =
    RigidHertzHeatError<E, RigidHertzRateTransientError, SourceHertzFrictionError>;
pub type HeatedRigidHertzRateCandidate<'m, P> =
    HeatedRigidHertzCandidate<'m, P, RigidHertzRateTransientReport, SourceHertzFriction>;

/// Independent original thermal, spatial-pressure and trajectory controls.
#[derive(Clone, Copy, Debug)]
pub struct RigidHertzRateHeatControls {
    pub heat: RigidHertzHeatControls,
    pub quadrature: HertzFrictionQuadrature,
    pub transient: HertzTransientControls,
}

/// Couple original local-speed Hertz recoil to the solved source temperature.
///
/// All curve parameters refresh at the common-face endpoint temperature, then
/// remain fixed for that physical phase. Spatial pressure quadrature and native
/// physical-time accuracy are independent of thermal-window refinement. The
/// source retains every pinned claim, typed role, fixed condition and authority.
/// Local speed belongs to the original spatial law, not to a single source
/// query coordinate. Source/domain, geometry, closed-normal, work, numerical or
/// cancellation failure returns no body, heat or source update.
///
/// Poses, normal load, radius and footprint are held. This passive branch adds
/// no external drive, static reaction selection, partial-slip history, moving
/// heat geometry or extra thermal storage. Only actual physical loss is heat.
#[allow(clippy::too_many_arguments)]
pub fn propose_rigid_hertz_rate_heat<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzRateHeatContact<'_>,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidHertzRateHeatControls,
    energy_guess_j: &[f64],
    temperature_guess_k: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzRateCandidate<'m, P>, RigidHertzRateHeatError<E>> {
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
            let prepared = PreparedRigidHertzRateTransient::new(
                bodies,
                properties,
                RigidHertzRateContact {
                    binding: contact.binding,
                    frame: contact.frame,
                    law: source.law(),
                },
                duration,
                controls.heat.mechanical,
                controls.quadrature,
                controls.transient,
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

/// Refine temperature feedback and heat release with the existing Hertz march.
///
/// Only accepted physical half phases survive; all original source and native
/// rate-domain receipts remain attached. The same whole-span compensated energy,
/// individual-cell/body errors, budgets and atomic-candidate semantics apply.
#[allow(clippy::too_many_arguments)]
pub fn propose_rigid_hertz_rate_heat_march<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzRateHeatContact<'_>,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    duration_s: f64,
    phase_controls: RigidHertzRateHeatControls,
    controls: march::HertzHeatMarchControls,
    temperature_guess_k: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<
    march::HeatedRigidHertzMarchCandidate<'m, P, SourceHertzFriction>,
    march::HertzHeatMarchError<E>,
> {
    march::propose_hertz_heat_phases(
        bodies,
        properties,
        contact,
        thermal,
        duration_s,
        phase_controls.heat,
        controls,
        temperature_guess_k,
        |b, props, pair, h, dt, heat, energy, temperature| {
            let start = h.time_s();
            let candidate = propose_rigid_hertz_rate_heat(
                b,
                props,
                pair,
                h,
                route,
                start,
                dt,
                RigidHertzRateHeatControls {
                    heat,
                    ..phase_controls
                },
                energy,
                temperature,
                &resolve,
                gate,
            )
            .map_err(march::HertzHeatMarchError::RateDependent)?;
            Ok(march::stage(
                candidate,
                start,
                march::HertzHeatMarchMotion::RateDependent,
            ))
        },
        gate,
    )
}
