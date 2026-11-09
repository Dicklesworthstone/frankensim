//! Curved original Hertz trajectories inside the existing finite-heat solve.
//!
//! Geometry and the full native mobility are prepared once. Every thermal
//! source trial restarts the same accepted pair, then integrates the original
//! pressure-resolved shear/spin traction with bounded physical-power quadrature.
//! Only that physical loss enters the original finite-cell thermal network.

use super::*;
use fs_mbd::tangential_friction::hertz::transient::{
    HertzTransientControls, PreparedRigidHertzTransient, RigidHertzTransientError,
    RigidHertzTransientReport,
};

/// Propose curved finite-body Hertz motion and source-temperature feedback.
///
/// The native owner retains the complete mobility and each law's actual force,
/// moment, path and finite-extinction tail. Its numerical work residual is
/// independently bounded and never added to heat. Original endpoint bodies,
/// thermal cells, friction source receipt and material payload leave the shared
/// transaction only after all native, thermal and combined energy checks pass.
///
/// This closed, held-geometry contact must preserve normal motion over the
/// complete path, as checked by the native normal-response bound. Poses,
/// normal load, radius and inertia remain fixed. Source coefficients are held
/// at the trial's solved endpoint temperature. Thermal transport is backward
/// Euler over the whole window, including any stationary remainder, and the
/// original face is area averaged. Thermal-window/spatial refinement remains
/// necessary for evolving coefficients and heat-release timing. Prefix error
/// estimates and a bounded analytic tail are not a validated flow enclosure.
#[allow(clippy::too_many_arguments)]
pub fn propose_rigid_hertz_transient_heat<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzHeatContact<'_>,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidHertzHeatControls,
    transient: HertzTransientControls,
    energy_guess_j: &[f64],
    temperature_guess_k: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<
    HeatedRigidHertzCandidate<'m, P, RigidHertzTransientReport>,
    RigidHertzHeatError<E, RigidHertzTransientError>,
> {
    solve_rigid_hertz_heat(
        contact,
        thermal,
        route,
        accepted_time_s,
        dt_s,
        controls,
        energy_guess_j,
        temperature_guess_k,
        |temperature| refresh_source(contact, temperature),
        |source, duration| {
            let prepared = PreparedRigidHertzTransient::new(
                bodies,
                properties,
                RigidHertzContact {
                    binding: contact.binding,
                    frame: contact.frame,
                    law: source.law(),
                },
                duration,
                controls.mechanical,
                transient,
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
