//! Explicit directional impulse friction at the original sphere contact events.
//!
//! Normal restitution and the event clock stay with the sphere/plane engine.
//! The impulse ellipse has its own explicitly declared coefficients and texture
//! history. It is instantaneous: it does not advance accumulated slip or borrow
//! the sustained law's static/kinetic distinction as an implicit impact model.
//! All impact loss remains separate from sustained traction-slip heat.

use super::*;
use fs_tribo::tangential::elliptic::impulse::{EllipticImpulseControls, EllipticImpulseError};

/// Original mechanical state at a detected impact, before applying its impulse.
#[derive(Clone, Copy, Debug)]
pub struct SpherePlaneDirectionalImpactPhase {
    pub input: EllipticImpulseInput,
    pub time_s: f64,
    pub body: RigidBodyState,
    pub frame: TangentFrame,
    pub normal_velocity_before_m_s: f64,
    pub normal_velocity_after_m_s: f64,
}

#[derive(Debug)]
pub enum SpherePlaneDirectionalLawError {
    Support(EllipticDrivenSegmentError<SpherePlaneError>),
    Impact(EllipticImpulseError),
}

/// Advance with separate native sustained and instantaneous directional laws.
///
/// `impact_law` must have equal static/kinetic principal coefficients, defining
/// one impulse ellipse. Its borrowed history supplies the retained material
/// axes and is never advanced. `history` belongs to sustained contact alone.
/// The plane owns the frame throughout; changing coordinates between calls
/// requires the original objective-history transport for both histories.
#[allow(clippy::too_many_arguments)]
pub fn advance_sphere_plane_directional_impacts(
    scene: &SpherePlane,
    state: SpherePlaneState,
    support_law: &EllipticCoulomb,
    history: &mut EllipticCoulombState,
    impact_law: &EllipticCoulomb,
    impact_history: &EllipticCoulombState,
    frame: TangentFrame,
    duration_s: f64,
    controls: SpherePlaneDirectionalControls,
    impact_controls: EllipticImpulseControls,
    cancelled: impl FnMut() -> bool,
) -> Result<
    SpherePlaneDirectionalAdvance,
    SpherePlaneDirectionalError<SpherePlaneDirectionalLawError>,
> {
    advance_sphere_plane_directional_with_impacts(
        scene,
        state,
        history,
        frame,
        duration_s,
        controls,
        |history, phase, visit, cancel| {
            let report = support_law
                .advance_driven_segments(
                    history,
                    &phase.input,
                    controls.driven,
                    &mut || cancel(),
                    &mut |segment| visit(segment),
                )
                .map_err(SpherePlaneDirectionalLawError::Support)?;
            let next = report.state.clone();
            Ok((report, next))
        },
        |history, phase, cancel| {
            let report = impact_law
                .solve_impulse(impact_history, &phase.input, impact_controls, &mut || {
                    cancel()
                })
                .map_err(SpherePlaneDirectionalLawError::Impact)?;
            Ok((report, history.clone()))
        },
        cancelled,
    )
}

/// Source-compatible support and impact callbacks sharing one accepted history.
///
/// The source owner explicitly chooses the instantaneous impact sampling rule
/// and returns its candidate history. There is no synthetic impact duration,
/// contact pressure, slip distance or heat allocation. Mechanical inputs must
/// retain their original meaning. The original body response independently
/// checks the applied impulse, actual point velocity and physical impact work.
///
/// Both callbacks operate on a private history; even failure after a preceding
/// impact or support phase leaves the caller's accepted history unchanged.
#[allow(clippy::too_many_arguments)]
pub fn advance_sphere_plane_directional_with_impacts<S: Clone, E>(
    scene: &SpherePlane,
    state: SpherePlaneState,
    history: &mut S,
    frame: TangentFrame,
    duration_s: f64,
    controls: SpherePlaneDirectionalControls,
    resolve_support: impl FnMut(
        &S,
        &SpherePlaneDirectionalPhase,
        &mut dyn FnMut(&EllipticDrivenSegment) -> Result<(), SpherePlaneError>,
        &mut dyn FnMut() -> bool,
    ) -> Result<(EllipticDrivenReport, S), E>,
    mut resolve_impact: impl FnMut(
        &S,
        &SpherePlaneDirectionalImpactPhase,
        &mut dyn FnMut() -> bool,
    ) -> Result<(EllipticImpulseReport, S), E>,
    cancelled: impl FnMut() -> bool,
) -> Result<SpherePlaneDirectionalAdvance, SpherePlaneDirectionalError<E>> {
    advance_directional_contact(
        scene,
        state,
        history,
        frame,
        duration_s,
        controls,
        resolve_support,
        |history, phase, cancel| resolve_impact(history, phase, cancel).map(Some),
        cancelled,
    )
}

pub(super) fn admit_impact(
    phase: &SpherePlaneDirectionalImpactPhase,
    report: &EllipticImpulseReport,
    controls: SpherePlaneDirectionalControls,
) -> Result<TangentialImpactAdvance, SpherePlaneError> {
    let before = phase.input.velocity_before_m_s;
    let w = phase.input.mobility_per_kg[0][0];
    let after = [0, 1].map(|i| w.mul_add(report.impulse_ns[i], before[i]));
    velocity_check(
        after,
        report.velocity_after_m_s,
        controls.velocity_tolerance_m_s,
    )?;
    let midpoint = [0, 1].map(|i| 0.5 * before[i] + 0.5 * after[i]);
    let physical_loss = finite(-dot(report.impulse_ns, midpoint))?;
    // The actual endpoint difference avoids cancellation between two much
    // larger kinetic energies when the impulse removes only a small amount.
    let delta_velocity = [after[0] - before[0], after[1] - before[1]];
    let kinetic_change = finite(dot(midpoint, delta_velocity) / w)?;
    if [
        report.dissipated_j,
        report.kinetic_energy_change_j,
        report.energy_residual_j,
    ]
    .iter()
    .any(|x| !x.is_finite())
        || report.dissipated_j < 0.0
        || physical_loss < -controls.work_tolerance_j
        || (physical_loss - report.dissipated_j).abs() > controls.work_tolerance_j
        || (kinetic_change - report.kinetic_energy_change_j).abs() > controls.work_tolerance_j
        || report.energy_residual_j.abs() > controls.work_tolerance_j
    {
        return Err(SpherePlaneError::Unresolved(
            "directional physical impact work",
        ));
    }
    Ok(TangentialImpactAdvance {
        impulse: report.impulse_ns,
        velocity_after: report.velocity_after_m_s,
        dissipated: report.dissipated_j,
    })
}
