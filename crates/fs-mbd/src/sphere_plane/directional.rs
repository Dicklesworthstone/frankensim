//! Plane-owned directional friction through the original moving-sphere engine.
//!
//! Every accepted held-traction segment advances original position, world
//! momenta and attitude. The original engine retains normal support, ballistic
//! timing, instantaneous impact, lift-off and complete work/momentum admission.
//! Support heat is integrated traction-slip work; impact loss remains separate.
//! No thermal receiver, contact area or impact-to-heat allocation is inferred.
//!
//! The caller owns a fixed texture frame on the prescribed plane. Sphere spin
//! does not rotate that texture. Across calls retain the same frame or use the
//! original objective-history transport before changing its coordinates.

use super::*;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::elliptic::driven::{
    EllipticDrivenControls, EllipticDrivenInput, EllipticDrivenMode, EllipticDrivenReport,
    EllipticDrivenSegment, EllipticDrivenSegmentError,
};
use fs_tribo::tangential::elliptic::impulse::{EllipticImpulseInput, EllipticImpulseReport};
use fs_tribo::tangential::elliptic::{EllipticCoulomb, EllipticCoulombState};
use std::cell::RefCell;

/// Explicit instantaneous directional friction through the same contact events.
pub mod impact;
use impact::SpherePlaneDirectionalImpactPhase;

#[derive(Clone, Copy, Debug)]
pub struct SpherePlaneDirectionalControls {
    /// Original geometry/event/attitude controls. max_sliding_substeps also
    /// bounds aggregate native support trials across all source phases.
    pub motion: SpherePlaneConfig,
    /// Independent accuracy axis for changes in source conditions.
    pub maximum_support_phase_s: f64,
    pub maximum_support_phases: usize,
    /// Accuracy and bounded work inside each held-source support phase.
    pub driven: EllipticDrivenControls,
    /// Actual sphere point velocity versus native segment endpoints.
    pub velocity_tolerance_m_s: f64,
    /// Independent native segment/aggregate physical work acceptance.
    pub work_tolerance_j: f64,
}

/// Actual accepted mechanical origin for a source-owned support phase.
#[derive(Clone, Copy, Debug)]
pub struct SpherePlaneDirectionalPhase {
    pub input: EllipticDrivenInput,
    pub time_s: f64,
    pub body: RigidBodyState,
    pub frame: TangentFrame,
}

#[derive(Clone, Debug)]
pub struct SpherePlaneDirectionalAdvance {
    pub motion: SpherePlaneAdvance,
    /// Sustained traction-slip heat only; no normal or tangential impact loss.
    pub sustained_friction_heat_j: f64,
    /// Original instantaneous normal plus tangential impact loss, unallocated.
    pub impact_loss_j: f64,
    pub support_distance_m: f64,
    /// Integrated relative contact displacement in the plane-owned frame.
    pub support_slip_m: [f64; 2],
    pub support_phases: usize,
    pub segments: usize,
}

#[derive(Debug)]
pub enum SpherePlaneDirectionalError<E = EllipticDrivenSegmentError<SpherePlaneError>> {
    Mechanics(SpherePlaneError),
    Source(E),
    Invalid(&'static str),
}
impl<E> From<SpherePlaneError> for SpherePlaneDirectionalError<E> {
    fn from(error: SpherePlaneError) -> Self {
        Self::Mechanics(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SpherePlaneDirectionalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "directional sphere/plane trajectory: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SpherePlaneDirectionalError<E> {}

/// Advance real sphere motion using the original native directional law.
#[allow(clippy::too_many_arguments)]
pub fn advance_sphere_plane_directional(
    scene: &SpherePlane,
    state: SpherePlaneState,
    law: &EllipticCoulomb,
    history: &mut EllipticCoulombState,
    frame: TangentFrame,
    duration_s: f64,
    controls: SpherePlaneDirectionalControls,
    cancelled: impl FnMut() -> bool,
) -> Result<SpherePlaneDirectionalAdvance, SpherePlaneDirectionalError> {
    advance_sphere_plane_directional_with(
        scene,
        state,
        history,
        frame,
        duration_s,
        controls,
        |history, phase, visit, cancel| {
            let report = law.advance_driven_segments(
                history,
                &phase.input,
                controls.driven,
                &mut || cancel(),
                &mut |segment| visit(segment),
            )?;
            let next = report.state.clone();
            Ok((report, next))
        },
        cancelled,
    )
}

/// Replace sustained friction while retaining the original trajectory owner.
///
/// The callback resolves the original accepted history and emits only the
/// accepted constant-traction segments. Source-query trials must not visit the
/// sphere; replay the accepted held coefficients through the same segmented
/// native method. The final report is independently checked against that replay.
/// Both the sphere and S remain private until every phase and the original
/// whole-trajectory geometry/momentum/work gates succeed.
///
/// The supplied input owns actual point velocity, isotropic sphere mobility,
/// external tangent acceleration, supported normal reaction and interval. It
/// leaves temperature/pressure None: a source consumer may populate those from
/// its declared physical owners, without altering the mechanical fields.
/// Each source phase has bounded length independently of native curved sliding
/// and attitude steps. Refinement of one axis does not certify another.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_sphere_plane_directional_with<S: Clone, E>(
    scene: &SpherePlane,
    state: SpherePlaneState,
    history: &mut S,
    frame: TangentFrame,
    duration_s: f64,
    controls: SpherePlaneDirectionalControls,
    resolve: impl FnMut(
        &S,
        &SpherePlaneDirectionalPhase,
        &mut dyn FnMut(&EllipticDrivenSegment) -> Result<(), SpherePlaneError>,
        &mut dyn FnMut() -> bool,
    ) -> Result<(EllipticDrivenReport, S), E>,
    cancelled: impl FnMut() -> bool,
) -> Result<SpherePlaneDirectionalAdvance, SpherePlaneDirectionalError<E>> {
    advance_directional_contact(
        scene,
        state,
        history,
        frame,
        duration_s,
        controls,
        resolve,
        |_, _, _| Ok(None),
        cancelled,
    )
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn advance_directional_contact<S: Clone, E>(
    scene: &SpherePlane,
    state: SpherePlaneState,
    history: &mut S,
    frame: TangentFrame,
    duration_s: f64,
    controls: SpherePlaneDirectionalControls,
    mut resolve: impl FnMut(
        &S,
        &SpherePlaneDirectionalPhase,
        &mut dyn FnMut(&EllipticDrivenSegment) -> Result<(), SpherePlaneError>,
        &mut dyn FnMut() -> bool,
    ) -> Result<(EllipticDrivenReport, S), E>,
    mut resolve_impact: impl FnMut(
        &S,
        &SpherePlaneDirectionalImpactPhase,
        &mut dyn FnMut() -> bool,
    ) -> Result<Option<(EllipticImpulseReport, S)>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<SpherePlaneDirectionalAdvance, SpherePlaneDirectionalError<E>> {
    if !controls.maximum_support_phase_s.is_finite()
        || controls.maximum_support_phase_s <= 0.0
        || controls.maximum_support_phases == 0
        || controls.maximum_support_phases > 65_536
        || [controls.velocity_tolerance_m_s, controls.work_tolerance_j]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(SpherePlaneDirectionalError::Invalid(
            "positive finite support phase limits and physical tolerances required",
        ));
    }
    let accepted_history = RefCell::new(history.clone());
    let source_error = RefCell::new(None);
    let mut support_phases = 0usize;
    let mut segments = 0usize;
    let mut support_distance_m = 0.0;
    let mut support_slip_m = [0.0; 2];
    let mut sustained_friction_heat_j = 0.0;
    let result = advance_with_contact(
        scene,
        state,
        // Unused by this support hook. Instantaneous tangent response comes
        // from the impact callback or the original scene-owned fallback.
        &FrictionLaw::Coulomb {
            static_mu: 0.0,
            kinetic_mu: 0.0,
        },
        duration_s,
        controls.motion,
        Some(frame),
        cancelled,
        |engine, body, time, remaining, normal, sphere_mode, cancelled| {
            if support_phases >= controls.maximum_support_phases {
                return Err(SpherePlaneError::Unresolved(
                    "directional support phase budget",
                ));
            }
            let dt = remaining.min(controls.maximum_support_phase_s);
            progress(time, dt)?;
            let (_, before) = engine.kinematics(body)?;
            let mode = if sphere_mode == SpherePlaneMode::Rolling || before == [0.0; 2] {
                EllipticDrivenMode::Resting
            } else {
                EllipticDrivenMode::Sliding
            };
            let w = engine.inverse_tangent_mass;
            let phase = SpherePlaneDirectionalPhase {
                input: EllipticDrivenInput {
                    velocity_m_s: before,
                    mobility_per_kg: [[w, 0.0], [0.0, w]],
                    driving_acceleration_m_s2: engine.acceleration_tangent,
                    mode,
                    normal_n: normal,
                    dt_s: dt,
                    maximum_force_n: controls.driven.solve.maximum_force_n,
                    contact_temperature_k: None,
                    nominal_pressure_pa: None,
                },
                time_s: time,
                body,
                frame,
            };
            let mut next = body;
            let mut elapsed = 0.0;
            let mut current_mode = mode;
            let mut impulse = [0.0; 2];
            let mut slip = [0.0; 2];
            let mut distance = 0.0;
            let mut loss = 0.0;
            let mut count = 0usize;
            let mut arrest = None;
            let cancel = RefCell::new(cancelled);
            let callback_result = resolve(
                &accepted_history.borrow(),
                &phase,
                &mut |segment| {
                    poll(&mut || (*cancel.borrow_mut())())?;
                    if segment.elapsed_before_s != elapsed
                        || !segment.dt_s.is_finite()
                        || segment.dt_s <= 0.0
                        || segment.dt_s > dt - elapsed
                        || segment.mode_before != current_mode
                        || [segment.distance_m, segment.dissipated_j]
                            .iter()
                            .any(|x| !x.is_finite() || *x < 0.0)
                        || segment
                            .traction_n
                            .iter()
                            .chain(&segment.velocity_before_m_s)
                            .chain(&segment.velocity_after_m_s)
                            .chain(&segment.slip_increment_m)
                            .any(|x| !x.is_finite())
                    {
                        return Err(SpherePlaneError::Invalid(
                            "inconsistent directional support segment",
                        ));
                    }
                    if count >= controls.driven.maximum_substeps {
                        return Err(SpherePlaneError::SlidingBudget);
                    }
                    let start = finite(time + elapsed)?;
                    progress(start, segment.dt_s)?;
                    let (_, v0) = engine.kinematics(next)?;
                    velocity_check(
                        v0,
                        segment.velocity_before_m_s,
                        controls.velocity_tolerance_m_s,
                    )?;
                    if current_mode == EllipticDrivenMode::Resting
                        && length(v0) > controls.driven.rest_velocity_tolerance_m_s
                    {
                        return Err(SpherePlaneError::Unresolved("accepted rolling residual"));
                    }
                    let force = engine.tangent_vector(segment.traction_n).scale(-1.0);
                    vector(force)?;
                    let advanced =
                        engine.segment(next, segment.dt_s, normal, force, &mut || {
                            (*cancel.borrow_mut())()
                        })?;
                    let (normal_after, v1) = engine.kinematics(advanced)?;
                    velocity_check(
                        v1,
                        segment.velocity_after_m_s,
                        controls.velocity_tolerance_m_s,
                    )?;
                    if normal_after.abs() > controls.velocity_tolerance_m_s
                        || segment.mode_after == EllipticDrivenMode::Resting
                            && segment.velocity_after_m_s != [0.0; 2]
                    {
                        return Err(SpherePlaneError::Unresolved(
                            "directional supported endpoint",
                        ));
                    }
                    let displacement = [0, 1].map(|i| 0.5 * segment.dt_s * (v0[i] + v1[i]));
                    velocity_check(
                        displacement,
                        segment.slip_increment_m,
                        segment.dt_s * controls.velocity_tolerance_m_s,
                    )?;
                    let physical_heat = finite(dot(segment.traction_n, displacement))?;
                    if physical_heat < -controls.work_tolerance_j
                        || (physical_heat - segment.dissipated_j).abs() > controls.work_tolerance_j
                        || segment.distance_m + segment.dt_s * controls.velocity_tolerance_m_s
                            < length(segment.slip_increment_m)
                    {
                        return Err(SpherePlaneError::Unresolved(
                            "directional physical support work",
                        ));
                    }
                    engine.admit_supported_path(next, advanced, start, segment.dt_s, force)?;
                    if current_mode == EllipticDrivenMode::Resting
                        && segment.mode_after == EllipticDrivenMode::Sliding
                    {
                        engine.event(start, SpherePlaneEventKind::Sliding, 0.0, Vec3::ZERO)?;
                    }
                    elapsed = if segment.dt_s == dt - elapsed {
                        dt
                    } else {
                        finite(elapsed + segment.dt_s)?
                    };
                    if current_mode == EllipticDrivenMode::Sliding
                        && segment.mode_after == EllipticDrivenMode::Resting
                    {
                        engine.event(
                            finite(time + elapsed)?,
                            SpherePlaneEventKind::Rolling,
                            0.0,
                            Vec3::ZERO,
                        )?;
                        if arrest.is_none() {
                            arrest = Some(elapsed);
                        }
                    }
                    for i in 0..2 {
                        impulse[i] = finite(impulse[i] + segment.dt_s * segment.traction_n[i])?;
                        slip[i] = finite(slip[i] + segment.slip_increment_m[i])?;
                    }
                    distance = finite(distance + segment.distance_m)?;
                    loss = finite(loss + segment.dissipated_j)?;
                    engine.ledger.loss = finite(engine.ledger.loss + segment.dissipated_j)?;
                    next = advanced;
                    current_mode = segment.mode_after;
                    count += 1;
                    Ok(())
                },
                &mut || (*cancel.borrow_mut())(),
            );
            let (native, next_history) = match callback_result {
                Ok(result) => result,
                Err(error) => {
                    *source_error.borrow_mut() = Some(error);
                    return Err(SpherePlaneError::Unresolved("directional source refused"));
                }
            };
            poll(&mut || (*cancel.borrow_mut())())?;
            if elapsed != dt
                || count == 0
                || native.mode != current_mode
                || native.first_arrest_time_s != arrest
                || [
                    native.distance_m,
                    native.dissipated_j,
                    native.drive_work_j,
                    native.kinetic_energy_change_j,
                    native.energy_residual_j,
                ]
                .iter()
                .any(|x| !x.is_finite())
            {
                return Err(SpherePlaneError::Unresolved(
                    "incomplete directional support replay",
                ));
            }
            let (_, after) = engine.kinematics(next)?;
            velocity_check(
                after,
                native.velocity_after_m_s,
                controls.velocity_tolerance_m_s,
            )?;
            velocity_check(
                impulse,
                native.traction_impulse_ns,
                controls.velocity_tolerance_m_s / w,
            )?;
            velocity_check(
                slip,
                native.slip_increment_m,
                dt * controls.velocity_tolerance_m_s,
            )?;
            let drive = engine.acceleration_tangent.map(|a| a / w);
            velocity_check(
                drive,
                native.driving_force_n,
                controls.velocity_tolerance_m_s / (dt * w),
            )?;
            let initial = if mode == EllipticDrivenMode::Resting {
                [0.0; 2]
            } else {
                before
            };
            let relative_kinetic = 0.5 * (dot(after, after) - dot(initial, initial)) / w;
            if (native.distance_m - distance).abs() > dt * controls.velocity_tolerance_m_s
                || (native.dissipated_j - loss).abs() > controls.work_tolerance_j
                || (native.drive_work_j - dot(drive, slip)).abs() > controls.work_tolerance_j
                || (native.kinetic_energy_change_j - relative_kinetic).abs()
                    > controls.work_tolerance_j
                || native.energy_residual_j.abs() > controls.work_tolerance_j
            {
                return Err(SpherePlaneError::Unresolved(
                    "directional aggregate work/trajectory",
                ));
            }
            engine.sliding_steps = engine
                .sliding_steps
                .checked_add(native.substeps)
                .ok_or(SpherePlaneError::SlidingBudget)?;
            if engine.sliding_steps > engine.config.max_sliding_substeps {
                return Err(SpherePlaneError::SlidingBudget);
            }
            segments = segments
                .checked_add(count)
                .ok_or(SpherePlaneError::SlidingBudget)?;
            support_phases += 1;
            support_distance_m = finite(support_distance_m + distance)?;
            sustained_friction_heat_j = finite(sustained_friction_heat_j + loss)?;
            for i in 0..2 {
                support_slip_m[i] = finite(support_slip_m[i] + slip[i])?;
            }
            *accepted_history.borrow_mut() = next_history;
            Ok(Some(SupportAdvance {
                body: next,
                mode: match current_mode {
                    EllipticDrivenMode::Resting => SpherePlaneMode::Rolling,
                    EllipticDrivenMode::Sliding => SpherePlaneMode::Sliding,
                },
                duration_s: dt,
            }))
        },
        |engine, body, time, normal_before, normal, velocity, cancelled| {
            let w = engine.inverse_tangent_mass;
            let phase = SpherePlaneDirectionalImpactPhase {
                input: EllipticImpulseInput {
                    normal_impulse_ns: normal.normal_impulse,
                    velocity_before_m_s: velocity,
                    mobility_per_kg: [[w, 0.0], [0.0, w]],
                },
                time_s: time,
                body,
                frame,
                normal_velocity_before_m_s: normal_before,
                normal_velocity_after_m_s: normal.relative_normal_speed_after,
            };
            let response = resolve_impact(&accepted_history.borrow(), &phase, &mut || cancelled());
            let (native, next_history) = match response {
                Ok(Some(response)) => response,
                Ok(None) => return Ok(None),
                Err(error) => {
                    *source_error.borrow_mut() = Some(error);
                    return Err(SpherePlaneError::Unresolved(
                        "directional impact source refused",
                    ));
                }
            };
            poll(cancelled)?;
            let accepted = impact::admit_impact(&phase, &native, controls)?;
            *accepted_history.borrow_mut() = next_history;
            Ok(Some(accepted))
        },
    );
    if let Some(error) = source_error.into_inner() {
        return Err(SpherePlaneDirectionalError::Source(error));
    }
    let (motion, impact_loss_j) = result?;
    // No fallible work or callback follows the original global admission.
    *history = accepted_history.into_inner();
    Ok(SpherePlaneDirectionalAdvance {
        motion,
        sustained_friction_heat_j,
        impact_loss_j,
        support_distance_m,
        support_slip_m,
        support_phases,
        segments,
    })
}

fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn length(v: [f64; 2]) -> f64 {
    v[0].hypot(v[1])
}
fn velocity_check(a: [f64; 2], b: [f64; 2], tolerance: f64) -> Result<(), SpherePlaneError> {
    if !tolerance.is_finite()
        || tolerance < 0.0
        || a.into_iter().chain(b).any(|x| !x.is_finite())
        || length([a[0] - b[0], a[1] - b[1]]) > tolerance
    {
        Err(SpherePlaneError::Unresolved(
            "directional trajectory vector residual",
        ))
    } else {
        Ok(())
    }
}
