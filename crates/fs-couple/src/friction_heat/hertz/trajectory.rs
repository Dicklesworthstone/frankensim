//! Original finite-pair Hertz stopping events with finite-cell thermal feedback.
//!
//! Each nonlinear trial starts from the same accepted bodies. The original
//! source claims at the solved common-face temperature supply one held Hertz
//! law for the whole interval. Only its moving phase releases frictional heat;
//! the existing thermal network conducts throughout the stationary remainder.
//! No static-work residual, endpoint clipping or extra heat store is introduced.

use super::*;
use fs_mbd::tangential_friction::hertz::trajectory::{
    PreparedRigidHertzTrajectory, RigidHertzTrajectoryError, RigidHertzTrajectoryReport,
};

/// Couple exact invariant-direction Hertz arrest to the original thermal cells.
///
/// The native full-mobility direction check is authoritative: pure sliding,
/// spin and compatible mixed paths retain their actual stopping time and travel.
/// Curved paths, external drive, pose evolution and changing normal load remain
/// outside this branch. No body velocity is projected into an admitted path.
/// Both native normal endpoint rates must satisfy the original velocity gate;
/// their linear moving phase and constant dwell then retain that bound throughout.
/// Original thermal face centre, normal, radius and A/B ordering are checked by
/// the same owner as the held-midpoint heat adapter.
///
/// The face is massless and area averaged; each original cell retains its own
/// heat capacity, temperature and conductivity. A source snapshot is held at
/// the trial's endpoint temperature, while thermal transport is backward Euler
/// over the complete interval. Exact held-coefficient arrest is not an exact
/// solution for continuously changing temperature or heat-release timing:
/// thermal-window and spatial refinement are required. Material advection,
/// moving footprints and radial flash temperature are not inferred.
#[allow(clippy::too_many_arguments)]
pub fn propose_rigid_hertz_trajectory_heat<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzHeatContact<'_>,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidHertzHeatControls,
    energy_guess_j: &[f64],
    temperature_guess_k: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<
    HeatedRigidHertzCandidate<'m, P, RigidHertzTrajectoryReport>,
    RigidHertzHeatError<E, RigidHertzTrajectoryError>,
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
        |source, duration| {
            let prepared = PreparedRigidHertzTrajectory::new(
                bodies,
                properties,
                RigidHertzContact {
                    binding: contact.binding,
                    frame: contact.frame,
                    law: source.law(),
                },
                duration,
                controls.mechanical,
                gate,
            )?;
            let point = prepared.contact_points_world()[0];
            Ok((prepared, [point.x, point.y, point.z]))
        },
        |prepared, source| {
            let candidate = prepared.advance(source.law(), gate)?;
            if [
                candidate.report.normal_velocity_before_m_s,
                candidate.report.normal_velocity_after_m_s,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v.abs() > controls.mechanical.velocity_tolerance_m_s)
            {
                return Err(RigidHertzError::Invalid(
                    "fixed thermal face requires held closed normal motion throughout arrest",
                )
                .into());
            }
            Ok(HertzHeatMechanics {
                bodies_after: candidate.bodies_after,
                dissipated_j: candidate.report.relaxation.dissipated_energy_j,
                total_energy_residual_j: candidate.report.total_energy_residual_j,
                report: candidate.report,
            })
        },
        resolve,
        gate,
    )
}
