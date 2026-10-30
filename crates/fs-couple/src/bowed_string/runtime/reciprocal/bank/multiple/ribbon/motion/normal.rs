//! Normal motion of the same material-fiber supports as the axial stroke.
//!
//! In the translating frame, `mu w_tt + T/EI/loss terms = contact - mu Y_tt`.
//! The source sine basis supplies `Gamma_n = integral(mu phi_n ds)`. Adding
//! `Y(t)` only to contact clearance would accelerate a massless support while
//! silently leaving every fiber inertial frame behind. We supply both terms,
//! and retain `Y_t sum(Gamma_n qdot_n) + M Y_t^2/2` in laboratory energy.
//! The normal contact owner still computes its exact discrete actuator work.
//! The two translations commute; there is no rotating-frame or flexible-stick
//! approximation hidden here. This is the rigid-support continuum rung.

use super::*;
use crate::bowed_string::runtime::reciprocal::bank::multiple::performance::fiber_motion::FiberSupportProgram;
use crate::render::schedule::force::coupled::contact::friction::multiple::obstacles::motion::ObstacleMotionLimits;

/// Accepted normal support motion in the original fiber normal coordinate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiberNormalSupportKinematics {
    pub velocity_m_s: f64,
    pub displacement_m: f64,
}

/// Normal contribution to the combined axial/normal laboratory frame report.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiberNormalSupportFrame {
    pub velocity_m_s: f64,
    pub translation_m: f64,
    pub acceleration_m_s2: f64,
    /// Already included in the native external-work account [J].
    pub inertial_coordinate_work_j: f64,
    pub laboratory_kinetic_correction_j: f64,
    /// `-sum R_i sign_i delta_Y`, already in native obstacle work [J].
    /// Subtract this once when substituting the complete support work into a
    /// laboratory balance. Migration work remains in normal_attachment_work_j.
    pub contact_translation_work_j: f64,
}

pub(super) struct NormalSupportState {
    program: FiberSupportProgram,
    accepted: FiberNormalSupportKinematics,
    cursor: usize,
    samples: u64,
    participation: Vec<f64>,
    component_start: usize,
    cells: usize,
    contacts: Vec<FiberNormalContact>,
}

pub(super) struct NormalSupportStep {
    accepted: FiberNormalSupportKinematics,
    next: usize,
    acceleration_m_s2: f64,
    previous_displacement_moment: f64,
    previous_kinetic_correction_j: f64,
    displacement_increment_m: f64,
    pub(super) offsets_m: Option<Vec<f64>>,
}

impl MultiBowStringBank {
    /// Add an explicit normal support trajectory to the existing axial frame.
    /// `samples` bounds its finite physical-clock horizon. Contact translation,
    /// normal inertia, source motion and observer state accept atomically. A
    /// caller-supplied independent offset vector must keep zero placeholders at
    /// these fiber contacts; other stops retain their own trajectory.
    pub fn with_fiber_normal_support_program(
        mut self,
        program: FiberSupportProgram,
        samples: u64,
        gate: &CancelGate,
    ) -> Result<Self, RenderError> {
        if gate.is_requested() {
            return Err(RenderError::Coupled(ModalCouplingError::Cancelled));
        }
        program.validate(samples, self.source.bank().options().sample_rate_hz, gate)?;
        let motion = self
            .bow_stations
            .as_ref()
            .and_then(|s| s.fiber_support.as_ref())
            .ok_or_else(|| {
                RenderError::Friction(error(
                    "normal fiber support needs the original axial support frame",
                ))
            })?;
        if self.system.samples_rendered() != 0 || motion.normal.is_some() {
            return Err(RenderError::Friction(error(
                "normal fiber support binds once before physical playback",
            )));
        }
        let fiber = admit_sources(&self.source).map_err(RenderError::Friction)?;
        let participation = fiber
            .model()
            .translation_participation(FiberDirection::Normal);
        let component_start = fiber.component_start();
        let cells = fiber.model().resolution().cells;
        let contacts = motion.normal_contacts.clone();
        // One rigid normal translation cannot move differently oriented strips.
        let direction = contacts.first().map(|c| self.source.stops()[c.stop].side);
        if contacts
            .iter()
            .any(|c| Some(self.source.stops()[c.stop].side) != direction)
        {
            return Err(RenderError::Friction(error(
                "normal fiber supports require one common physical normal direction",
            )));
        }
        if let Some(direction) = direction {
            let source =
                &self.source.bank().strings()[self.source.bows()[0].location.string].source;
            let MaterialStringSource::WithPolarizations(polarized) = source else {
                return Err(RenderError::Friction(error(
                    "normal fiber translation needs an explicit transverse frame",
                )));
            };
            let stops::StringStopSide::Along(normal) = direction else {
                return Err(RenderError::Friction(error(
                    "normal fiber translation needs an explicit physical normal",
                )));
            };
            let n = normal.components();
            let t = polarized.directions().bow();
            if (n[0] * t[0] + n[1] * t[1]).abs() > 1e-12 {
                return Err(RenderError::Friction(error(
                    "normal fiber translation must be perpendicular to its axial stroke",
                )));
            }
        }
        let mut displacement_bound = 0.0;
        let mut momentum_bound = 0.0;
        for component in &self.system.components()[component_start..component_start + cells] {
            for gamma in &participation {
                displacement_bound +=
                    gamma.abs() * component.budget().maximum_abs_displacement_m_sqrt_kg;
                momentum_bound +=
                    gamma.abs() * component.budget().maximum_abs_velocity_m_sqrt_kg_per_s;
            }
        }
        let speed = program.maximum_speed_m_s;
        let correction_bound = speed * momentum_bound + 0.5 * motion.mass_kg * speed * speed;
        let inertial_bound = 2.0 * program.maximum_acceleration_m_s2 * displacement_bound;
        let normal_bound = contacts
            .iter()
            .map(|c| {
                self.source.stops()[c.stop].limits.maximum_force_n
                    * speed
                    * self.system.sample_period_s()
            })
            .sum::<f64>();
        if participation.iter().any(|x| !x.is_finite())
            || !(2.0 * correction_bound + inertial_bound + normal_bound).is_finite()
            || 2.0 * correction_bound + inertial_bound + normal_bound > f64::MAX / 16.0
        {
            return Err(RenderError::Friction(error(
                "normal fiber frame exceeds representable inertia/work bounds",
            )));
        }
        if !contacts.is_empty() {
            let limits = ObstacleMotionLimits {
                maximum_abs_offset_m: program.maximum_abs_translation_m,
                maximum_speed_m_s: speed,
            };
            if let Some(existing) = self.system.obstacle_motion_limits() {
                if existing.maximum_abs_offset_m < limits.maximum_abs_offset_m
                    || existing.maximum_speed_m_s < limits.maximum_speed_m_s
                {
                    return Err(RenderError::Friction(error(
                        "existing contact motion limits must cover the normal fiber program",
                    )));
                }
            } else {
                self.system = self
                    .system
                    .with_obstacle_motion(limits, gate)
                    .map_err(RenderError::Friction)?;
            }
        }
        let state = self.bow_stations.as_mut().expect("admitted support frame");
        let max_update_terms = state.limits.max_update_terms.min(program.max_update_terms);
        self.system
            .validate_friction_geometry_config(FrictionGeometryConfig { max_update_terms })
            .map_err(RenderError::Friction)?;
        state.limits.max_update_terms = max_update_terms;
        let accepted = FiberNormalSupportKinematics {
            velocity_m_s: program.initial_velocity_m_s,
            displacement_m: 0.0,
        };
        self.bow_stations
            .as_mut()
            .expect("admitted support frame")
            .fiber_support
            .as_mut()
            .expect("admitted fiber support")
            .normal = Some(NormalSupportState {
            program,
            accepted,
            cursor: 0,
            samples,
            participation,
            component_start,
            cells,
            contacts,
        });
        Ok(self)
    }

    #[must_use]
    pub fn fiber_normal_support_program(&self) -> Option<&FiberSupportProgram> {
        Some(
            &self
                .bow_stations
                .as_ref()?
                .fiber_support
                .as_ref()?
                .normal
                .as_ref()?
                .program,
        )
    }

    #[must_use]
    pub fn fiber_normal_support_kinematics(&self) -> Option<FiberNormalSupportKinematics> {
        Some(
            self.bow_stations
                .as_ref()?
                .fiber_support
                .as_ref()?
                .normal
                .as_ref()?
                .accepted,
        )
    }
}

impl NormalSupportState {
    fn moments(&self, system: &MultiFrictionModalSystem) -> (f64, f64) {
        system.components()[self.component_start..self.component_start + self.cells]
            .iter()
            .flat_map(|component| component.states().iter().zip(&self.participation))
            .fold((0.0, 0.0), |(q, v), (state, gamma)| {
                (
                    q + gamma * state.displacement_m_sqrt_kg,
                    v + gamma * state.velocity_m_sqrt_kg_per_s,
                )
            })
    }

    pub(super) fn prepare(
        &self,
        system: &MultiFrictionModalSystem,
        force_start: usize,
        mass_kg: f64,
        forces: &mut [f64],
        independent_offsets: Option<&[f64]>,
        gate: &CancelGate,
    ) -> Result<NormalSupportStep, ModalFrictionError> {
        let endpoint = system
            .samples_rendered()
            .checked_add(1)
            .ok_or_else(|| error("normal fiber clock overflow"))?;
        if endpoint > self.samples {
            return Err(error(
                "normal fiber support exceeds its declared finite horizon",
            ));
        }
        let (velocity, next) = self.program.sample_from(self.cursor, endpoint);
        let previous = self.accepted.velocity_m_s;
        let dt = system.sample_period_s();
        let acceleration = (velocity - previous) / dt;
        let increment = dt * (0.5 * previous + 0.5 * velocity);
        let displacement = self.accepted.displacement_m + increment;
        // Interpolation arithmetic can differ by ulps from the admitted segment
        // derivative. This tolerance admits only that rounding, never clipping.
        let acceleration_roundoff = 16.0 * f64::EPSILON * self.program.maximum_speed_m_s / dt;
        if !velocity.is_finite()
            || velocity.abs() > self.program.maximum_speed_m_s
            || !acceleration.is_finite()
            || acceleration.abs() > self.program.maximum_acceleration_m_s2 + acceleration_roundoff
            || !displacement.is_finite()
            || displacement.abs() > self.program.maximum_abs_translation_m
        {
            return Err(error(
                "normal fiber motion exceeds its declared physical limits",
            ));
        }
        let offsets_m = if self.contacts.is_empty() {
            None
        } else {
            let accepted = system
                .obstacle_offsets_m()
                .ok_or_else(|| error("normal fiber contact translation was not admitted"))?;
            if independent_offsets.is_some_and(|s| s.len() != accepted.len()) {
                return Err(error(
                    "normal fiber motion requires complete independent stop offsets",
                ));
            }
            let mut offsets = independent_offsets.unwrap_or(accepted).to_vec();
            for contact in &self.contacts {
                if independent_offsets.is_some_and(|s| s[contact.stop] != 0.0) {
                    return Err(error(
                        "normal fiber motion cannot share its contact with another translation drive",
                    ));
                }
                offsets[contact.stop] = contact.sign * displacement;
            }
            Some(offsets)
        };
        for cell in 0..self.cells {
            if gate.is_requested() {
                return Err(ModalCouplingError::Cancelled.into());
            }
            let start = force_start + cell * self.participation.len();
            for (force, gamma) in forces[start..start + self.participation.len()]
                .iter_mut()
                .zip(&self.participation)
            {
                *force -= acceleration * gamma;
                if !force.is_finite() {
                    return Err(error("normal fiber inertia is unrepresentable"));
                }
            }
        }
        let (previous_displacement_moment, momentum) = self.moments(system);
        let previous_kinetic_correction_j =
            previous * momentum + 0.5 * mass_kg * previous * previous;
        Ok(NormalSupportStep {
            accepted: FiberNormalSupportKinematics {
                velocity_m_s: velocity,
                displacement_m: displacement,
            },
            next,
            acceleration_m_s2: acceleration,
            previous_displacement_moment,
            previous_kinetic_correction_j,
            displacement_increment_m: displacement - self.accepted.displacement_m,
            offsets_m,
        })
    }

    pub(super) fn accept(
        &mut self,
        step: NormalSupportStep,
        system: &MultiFrictionModalSystem,
        mass_kg: f64,
    ) -> (FiberNormalSupportFrame, f64) {
        let (displacement, momentum) = self.moments(system);
        let velocity = step.accepted.velocity_m_s;
        let correction = velocity * momentum + 0.5 * mass_kg * velocity * velocity;
        let inertial = -step.acceleration_m_s2 * (displacement - step.previous_displacement_moment);
        let contact_work = self
            .contacts
            .iter()
            .map(|contact| {
                let reaction = system
                    .last_obstacle_frames()
                    .expect("accepted normal contacts")[contact.stop]
                    .normal_force_n;
                -reaction * contact.sign * step.displacement_increment_m
            })
            .sum();
        let report = FiberNormalSupportFrame {
            velocity_m_s: velocity,
            translation_m: step.accepted.displacement_m,
            acceleration_m_s2: step.acceleration_m_s2,
            inertial_coordinate_work_j: inertial,
            laboratory_kinetic_correction_j: correction,
            contact_translation_work_j: contact_work,
        };
        self.accepted = step.accepted;
        self.cursor = step.next;
        (report, step.previous_kinetic_correction_j)
    }
}
