//! Sustained finite-patch shear/spin through the original sphere trajectory.
//!
//! The original distributed Hertz wrench set owns joint static feasibility and
//! kinetic response. Normal support supplies its actual load; the footprint is
//! separately declared. A free normal-axis couple acts together with traction
//! in the original momentum, attitude and world-work update. No point-contact
//! rolling torque, effective friction ellipse or independent saturated spin
//! resistance substitutes for the Hertz surface.
//!
//! Pure spinning and collinear pure sliding have exact constant-wrench stopping
//! events. Mixed motion uses bounded velocity-implicit steps, with first-order
//! continuum convergence and timestep-boundary capture. Every accepted held
//! wrench opposes motion at both ends of its actual linear velocity path. Its
//! loss is actual midpoint work, not the projection's endpoint support pairing.
//! Static capture from mixed motion requires interior kinetic-load feasibility;
//! unresolved boundary/reversal trials refine or refuse. No velocity is reset.

use super::*;
use core::cell::Cell;
use core::convert::Infallible;
use fs_tribo::sliding_patch::contact::{HertzWrenchConfig, HertzWrenchContact, HertzWrenchError};

#[path = "hertz_support/step.rs"]
mod step;

/// Accepted native motion and its independent full-patch sticking mode.
/// Zero centre slip alone does not establish zero normal-axis spin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpherePlaneHertzState {
    motion: SpherePlaneState,
    zero_motion: bool,
    patch_sticking: bool,
}

impl SpherePlaneHertzState {
    /// Admit a new incoming state without assuming static capture. Exact zero
    /// relative shear and spin are classified on the first actual support phase.
    pub const fn new(motion: SpherePlaneState) -> Self {
        Self {
            motion,
            zero_motion: false,
            patch_sticking: false,
        }
    }

    pub const fn motion(self) -> SpherePlaneState {
        self.motion
    }

    pub const fn patch_sticking(self) -> bool {
        self.patch_sticking
    }

    /// Resolved zero shear/spin, including an instantaneous overdriven zero
    /// that cannot remain statically supported. This retains event continuation
    /// without replacing the native rounded momenta with zeros.
    pub const fn at_zero_motion(self) -> bool {
        self.zero_motion
    }
}

/// Original motion/constitutive budgets and independently bounded source phases.
#[derive(Clone, Copy, Debug)]
pub struct SpherePlaneHertzSupportControls {
    /// Mixed slip uses max_sliding_step_s. Accepted and rejected support trials
    /// together consume max_sliding_substeps; attitude has its separate bounds.
    pub motion: SpherePlaneConfig,
    pub wrench: HertzWrenchConfig,
    pub maximum_support_phase_s: f64,
    pub maximum_support_phases: usize,
    /// Native endpoint shear velocity and retained sticking residual [m/s].
    pub velocity_tolerance_m_s: f64,
    /// Native endpoint normal spin and retained sticking residual [rad/s].
    pub spin_tolerance_rad_s: f64,
    /// Physical segment work residual [J], independently of native solve work.
    pub work_tolerance_j: f64,
    /// Sum of absolute static constraint-work defects over the full advance [J].
    /// This does not permit many individually small static defects to accumulate.
    pub maximum_absolute_static_work_residual_j: f64,
}

/// Actual incoming mechanical state for one held-source support interval.
#[derive(Clone, Copy, Debug)]
pub struct SpherePlaneHertzSupportPhase {
    pub time_s: f64,
    pub duration_s: f64,
    pub body: RigidBodyState,
    pub contact_point_world_m: Vec3,
    pub normal_world: Vec3,
    pub tangent_axes_world: [Vec3; 2],
    pub normal_force_n: f64,
    pub slip_velocity_m_s: [f64; 2],
    pub spin_rad_s: f64,
    pub tangent_mobility_per_kg: f64,
    pub spin_mobility_per_kg_m2: f64,
    pub driving_tangent_acceleration_m_s2: [f64; 2],
    pub driving_spin_acceleration_rad_s2: f64,
    pub at_zero_motion: bool,
    pub patch_sticking: bool,
}

/// Full-patch transitions. Point-contact rolling events keep their old meaning.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpherePlaneHertzSupportTransition {
    pub time_s: f64,
    pub at_zero_motion: bool,
    pub patch_sticking: bool,
}

#[derive(Clone, Debug)]
pub struct SpherePlaneHertzSupportAdvance {
    pub state_after: SpherePlaneHertzState,
    pub motion: SpherePlaneAdvance,
    /// Actual shear plus spin work absorbed during persistent support [J].
    /// Already included in motion.dissipated_energy; no heat split is inferred.
    pub sustained_friction_loss_j: f64,
    pub impact_loss_j: f64,
    /// Sum of absolute bounded static work defects; the signed sum is retained
    /// in motion.static_contact_work_residual and its global energy balance.
    pub absolute_static_work_residual_j: f64,
    /// Actual signed centre-slip displacement and normal-axis rotation.
    pub support_slip_m: [f64; 2],
    pub support_spin_angle_rad: f64,
    /// Free contact angular impulse; also included in the native world balance.
    pub support_spin_impulse_world_nms: Vec3,
    pub support_phases: usize,
    pub support_segments: usize,
    pub transitions: Vec<SpherePlaneHertzSupportTransition>,
}

#[derive(Debug)]
pub enum SpherePlaneHertzSupportError<E = Infallible> {
    Mechanics(SpherePlaneError),
    Law(HertzWrenchError),
    Source(E),
}
impl<E> From<SpherePlaneError> for SpherePlaneHertzSupportError<E> {
    fn from(error: SpherePlaneError) -> Self {
        Self::Mechanics(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SpherePlaneHertzSupportError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-patch sphere/plane support: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SpherePlaneHertzSupportError<E> {}

/// Advance the real sphere with one held original finite-patch support law.
/// Its normal force must equal the actual support reaction when contact occurs.
/// The scene's original point impulse law remains the separate impact model.
#[allow(clippy::too_many_arguments)]
pub fn advance_sphere_plane_hertz_support(
    scene: &SpherePlane,
    state: SpherePlaneHertzState,
    torque_world_nm: Vec3,
    law: &HertzWrenchContact,
    duration_s: f64,
    controls: SpherePlaneHertzSupportControls,
    cancelled: impl FnMut() -> bool,
) -> Result<SpherePlaneHertzSupportAdvance, SpherePlaneHertzSupportError> {
    advance_sphere_plane_hertz_support_with(
        scene,
        state,
        torque_world_nm,
        &mut (),
        duration_s,
        controls,
        |_, _, _| Ok((*law, ())),
        cancelled,
    )
}

/// Resolve original load/radius/static/kinetic parameters at accepted physical
/// phase origins. The source returns a constitutive owner, never an arbitrary
/// wrench. Coefficients remain held within the declared source interval; native
/// slip refinement does not resolve a changing source. Source history and all
/// motion publish only after complete trajectory admission, including any later
/// impact, cancellation or source failure. A translating plane has zero spin;
/// its real work is force contracted with prescribed translational velocity.
#[allow(clippy::too_many_arguments)]
pub fn advance_sphere_plane_hertz_support_with<S: Clone, E>(
    scene: &SpherePlane,
    state: SpherePlaneHertzState,
    torque_world_nm: Vec3,
    history: &mut S,
    duration_s: f64,
    controls: SpherePlaneHertzSupportControls,
    resolve: impl FnMut(
        &S,
        &SpherePlaneHertzSupportPhase,
        &mut dyn FnMut() -> bool,
    ) -> Result<(HertzWrenchContact, S), E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SpherePlaneHertzSupportAdvance, SpherePlaneHertzSupportError<E>> {
    poll(&mut cancelled)?;
    validate_controls(controls)?;
    let (result, candidate) = candidate_with_impact(
        scene,
        state,
        torque_world_nm,
        history.clone(),
        duration_s,
        controls,
        resolve,
        |_, _, _, _, _, _, _| Ok(None),
        cancelled,
    )?;
    *history = candidate;
    Ok(result)
}

pub(super) fn validate_controls<E>(
    controls: SpherePlaneHertzSupportControls,
) -> Result<(), SpherePlaneHertzSupportError<E>> {
    controls
        .wrench
        .validate()
        .map_err(SpherePlaneHertzSupportError::Law)?;
    if !controls.maximum_support_phase_s.is_finite()
        || controls.maximum_support_phase_s <= 0.0
        || !(1..=65_536).contains(&controls.maximum_support_phases)
        || [
            controls.velocity_tolerance_m_s,
            controls.spin_tolerance_rad_s,
            controls.work_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || !controls.maximum_absolute_static_work_residual_j.is_finite()
        || controls.maximum_absolute_static_work_residual_j < 0.0
    {
        return Err(SpherePlaneError::Invalid("Hertz support phase or residual limits").into());
    }
    Ok(())
}

// Both public trajectories use this same support driver and original complete
// motion audit. Returning the staged source lets a coupled impact owner accept
// its independent history in the same transaction without cloning it twice.
#[allow(clippy::too_many_arguments)]
pub(super) fn candidate_with_impact<S, E, C: FnMut() -> bool>(
    scene: &SpherePlane,
    state: SpherePlaneHertzState,
    torque_world_nm: Vec3,
    mut candidate_history: S,
    duration_s: f64,
    controls: SpherePlaneHertzSupportControls,
    mut resolve: impl FnMut(
        &S,
        &SpherePlaneHertzSupportPhase,
        &mut dyn FnMut() -> bool,
    ) -> Result<(HertzWrenchContact, S), E>,
    mut impact: impl FnMut(
        &Engine,
        RigidBodyState,
        f64,
        f64,
        &NormalImpactResult,
        [f64; 2],
        &mut C,
    ) -> Result<Option<TangentialImpactAdvance>, SpherePlaneError>,
    mut cancelled: C,
) -> Result<(SpherePlaneHertzSupportAdvance, S), SpherePlaneHertzSupportError<E>> {
    poll(&mut cancelled)?;
    validate_controls(controls)?;
    let mut source_error = None;
    let resting = Cell::new(state.zero_motion);
    let sticking = Cell::new(state.patch_sticking);
    let mut loss = 0.0;
    let mut absolute_static_work = 0.0;
    let mut slip = [0.0; 2];
    let mut spin_angle = 0.0;
    let mut spin_impulse = Vec3::ZERO;
    let mut phases = 0usize;
    let mut segments = 0usize;
    let mut transitions = Vec::new();
    let result = advance_with_contact(
        scene,
        state.motion,
        torque_world_nm,
        &FrictionLaw::Coulomb {
            static_mu: 0.0,
            kinetic_mu: 0.0,
        },
        duration_s,
        controls.motion,
        None,
        &mut cancelled,
        |engine, mut body, time, remaining, normal, _, cancel| {
            if phases >= controls.maximum_support_phases {
                return Err(SpherePlaneError::Unresolved("Hertz source phase budget"));
            }
            let dt = remaining.min(controls.maximum_support_phase_s);
            progress(time, dt)?;
            let before = contact_velocity(engine, body)?;
            if before == [0.0; 3] {
                resting.set(true);
            }
            if resting.get() {
                velocity_check(before, [0.0; 3], controls)?;
            }
            let phase = SpherePlaneHertzSupportPhase {
                time_s: time,
                duration_s: dt,
                body,
                contact_point_world_m: body.pose().position_world().add(engine.arm),
                normal_world: engine.n,
                tangent_axes_world: engine.tangent,
                normal_force_n: normal,
                slip_velocity_m_s: [before[0], before[1]],
                spin_rad_s: before[2],
                tangent_mobility_per_kg: engine.inverse_tangent_mass,
                spin_mobility_per_kg_m2: engine.inertia.recip(),
                driving_tangent_acceleration_m_s2: engine.acceleration_tangent,
                driving_spin_acceleration_rad_s2: finite(
                    engine.external_torque_world.dot(engine.n) / engine.inertia,
                )?,
                at_zero_motion: resting.get(),
                patch_sticking: sticking.get(),
            };
            let (law, next_history) = match resolve(&candidate_history, &phase, &mut || cancel()) {
                Ok(result) => result,
                Err(error) => {
                    source_error = Some(SpherePlaneHertzSupportError::Source(error));
                    return Err(SpherePlaneError::Unresolved("Hertz support source refused"));
                }
            };
            poll(cancel)?;
            if law.parameters().normal_force_n != normal {
                return Err(SpherePlaneError::Invalid(
                    "Hertz law must retain actual supported load",
                ));
            }
            let drive = [
                phase.driving_tangent_acceleration_m_s2[0],
                phase.driving_tangent_acceleration_m_s2[1],
                phase.driving_spin_acceleration_rad_s2,
            ];
            let mut elapsed = 0.0;
            while elapsed < dt {
                poll(cancel)?;
                let v0 = contact_velocity(engine, body)?;
                let was_resting = resting.get();
                let was_sticking = sticking.get();
                if was_resting {
                    velocity_check(v0, [0.0; 3], controls)?;
                }
                let start = finite(time + elapsed)?;
                let trial = match step::resolve(
                    engine,
                    &law,
                    if was_resting { [0.0; 3] } else { v0 },
                    drive,
                    was_resting,
                    start,
                    dt - elapsed,
                    controls,
                    cancel,
                ) {
                    Ok(trial) => trial,
                    Err(step::Error::Mechanics(error)) => return Err(error),
                    Err(step::Error::Law(error)) => {
                        source_error = Some(SpherePlaneHertzSupportError::Law(error));
                        return Err(SpherePlaneError::Unresolved("Hertz support law refused"));
                    }
                };
                let force = engine.tangent_vector(trial.force_n);
                let next = engine.segment_with_spin(
                    body,
                    trial.dt_s,
                    normal,
                    force,
                    trial.spin_moment_nm,
                    cancel,
                )?;
                let v1 = contact_velocity(engine, next)?;
                velocity_check(v1, trial.velocity_after, controls)?;
                if trial.resting {
                    velocity_check(v1, [0.0; 3], controls)?;
                }
                let mut increment = [0.0; 3];
                for i in 0..3 {
                    increment[i] = finite(trial.dt_s * (0.5 * v0[i] + 0.5 * v1[i]))?;
                }
                let effort_impulse = Vec3::new(
                    trial.force_n[0] * trial.dt_s,
                    trial.force_n[1] * trial.dt_s,
                    trial.spin_moment_nm * trial.dt_s,
                );
                vector(effort_impulse)?;
                let absorbed = finite(
                    -crate::midpoint_impulse_pairing(
                        Vec3::new(v0[0], v0[1], v0[2]),
                        Vec3::new(v1[0], v1[1], v1[2]),
                        effort_impulse,
                    )
                    .to_f64(),
                )?;
                let native_work = finite(
                    -trial.force_n[0] * increment[0]
                        - trial.force_n[1] * increment[1]
                        - trial.spin_moment_nm * increment[2],
                )?;
                if (absorbed - native_work).abs() > controls.work_tolerance_j
                    || absorbed < 0.0
                        && !(was_resting
                            && trial.resting
                            && absorbed.abs() <= controls.work_tolerance_j)
                {
                    return Err(SpherePlaneError::Unresolved(
                        "Hertz support physical friction work",
                    ));
                }
                // Retained static roundoff is a bounded work residual, never
                // physical friction heat or a velocity-correction impulse.
                let physical_loss = if was_resting && trial.resting {
                    if absorbed.abs() > controls.work_tolerance_j {
                        return Err(SpherePlaneError::Unresolved("static Hertz work residual"));
                    }
                    absolute_static_work = finite(absolute_static_work + absorbed.abs())?;
                    if absolute_static_work > controls.maximum_absolute_static_work_residual_j {
                        return Err(SpherePlaneError::Unresolved(
                            "cumulative static Hertz work budget",
                        ));
                    }
                    engine.ledger.static_work = finite(engine.ledger.static_work - absorbed)?;
                    0.0
                } else {
                    absorbed
                };
                engine.admit_supported_path(body, next, start, trial.dt_s, force)?;
                engine.ledger.loss = finite(engine.ledger.loss + physical_loss)?;
                loss = finite(loss + physical_loss)?;
                for i in 0..2 {
                    slip[i] = finite(slip[i] + increment[i])?;
                }
                spin_angle = finite(spin_angle + increment[2])?;
                spin_impulse = spin_impulse.add(engine.n.scale(effort_impulse.z));
                vector(spin_impulse)?;
                elapsed = if trial.dt_s == dt - elapsed {
                    dt
                } else {
                    finite(elapsed + trial.dt_s)?
                };
                if was_resting != trial.resting || was_sticking != trial.sticking {
                    if transitions.len() >= controls.motion.max_events {
                        return Err(SpherePlaneError::EventBudget);
                    }
                    transitions.push(SpherePlaneHertzSupportTransition {
                        time_s: if was_resting {
                            start
                        } else {
                            finite(time + elapsed)?
                        },
                        at_zero_motion: trial.resting,
                        patch_sticking: trial.sticking,
                    });
                }
                resting.set(trial.resting);
                sticking.set(trial.sticking);
                body = next;
                segments += 1;
            }
            phases += 1;
            candidate_history = next_history;
            let velocity = contact_velocity(engine, body)?;
            let mode = if velocity[0] == 0.0 && velocity[1] == 0.0 || resting.get() {
                SpherePlaneMode::Rolling
            } else {
                SpherePlaneMode::Sliding
            };
            Ok(Some(SupportAdvance {
                body,
                mode,
                duration_s: dt,
            }))
        },
        |engine, body, time, normal_before, normal, slip, cancel| {
            // An impact is a distinct constitutive event. Its actual endpoint
            // is freshly classified on support; it cannot inherit capture.
            resting.set(false);
            sticking.set(false);
            impact(engine, body, time, normal_before, normal, slip, cancel)
        },
    );
    if let Some(error) = source_error {
        return Err(error);
    }
    let (motion, impact_loss_j) = result?;
    poll(&mut cancelled)?;
    let state_after = SpherePlaneHertzState {
        motion: motion.state_after,
        zero_motion: resting.get() && motion.state_after.mode != SpherePlaneMode::Flight,
        patch_sticking: sticking.get() && motion.state_after.mode != SpherePlaneMode::Flight,
    };
    Ok((
        SpherePlaneHertzSupportAdvance {
            state_after,
            motion,
            sustained_friction_loss_j: loss,
            impact_loss_j,
            absolute_static_work_residual_j: absolute_static_work,
            support_slip_m: slip,
            support_spin_angle_rad: spin_angle,
            support_spin_impulse_world_nms: spin_impulse,
            support_phases: phases,
            support_segments: segments,
            transitions,
        },
        candidate_history,
    ))
}

fn contact_velocity(
    engine: &mut Engine,
    body: RigidBodyState,
) -> Result<[f64; 3], SpherePlaneError> {
    let (_, tangent) = engine.kinematics(body)?;
    Ok([
        tangent[0],
        tangent[1],
        finite(world_spin(body).dot(engine.n) / engine.inertia)?,
    ])
}

fn velocity_check(
    actual: [f64; 3],
    target: [f64; 3],
    controls: SpherePlaneHertzSupportControls,
) -> Result<(), SpherePlaneError> {
    if actual.into_iter().chain(target).any(|x| !x.is_finite())
        || (actual[0] - target[0]).hypot(actual[1] - target[1]) > controls.velocity_tolerance_m_s
        || (actual[2] - target[2]).abs() > controls.spin_tolerance_rad_s
    {
        return Err(SpherePlaneError::Unresolved(
            "native Hertz support endpoint velocity",
        ));
    }
    Ok(())
}
