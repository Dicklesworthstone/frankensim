//! Continuous velocity programs for the original moving-fiber support owner.
//!
//! This scheduler supplies endpoint velocities only. The existing ribbon motion
//! owner supplies source-derived inertia, moving material contact, midpoint
//! traction and laboratory work accounting. No second force or motion law is
//! applied. Canonical normal contacts share that owner; shear/thermal patch
//! history transport remains unsupported.

pub use super::super::ribbon::motion::FiberSupportFrame;
use super::super::ribbon::motion::FiberSupportMotionLimits;
use super::{MultiBowPerformance, render_invalid};
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::RenderError;
use crate::render::schedule::force::coupled::ModalCouplingError;
use crate::render::schedule::force::coupled::contact::friction::FrictionDrive;
use fs_exec::CancelGate;

/// A support-velocity endpoint on the original mechanical clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiberSupportKnot {
    /// Positive elapsed samples. Sample 1 is the END of the first interval.
    pub sample: u64,
    /// Signed axial speed of both fiber supports [m/s].
    pub velocity_m_s: f64,
}

/// Continuous piecewise-linear support velocity, with explicit physical limits.
/// A final speed is held; it does not stop automatically at the final knot.
#[derive(Clone, Debug, PartialEq)]
pub struct FiberSupportProgram {
    /// Initial support-frame speed, independent of any preexisting relative modes.
    pub initial_velocity_m_s: f64,
    /// Strictly increasing physical endpoints, followed by a held final speed.
    pub knots: Vec<FiberSupportKnot>,
    /// Maximum magnitude of the prescribed axial speed [m/s].
    pub maximum_speed_m_s: f64,
    /// Maximum magnitude of each velocity segment's acceleration [m/s^2].
    pub maximum_acceleration_m_s2: f64,
    /// Bound absolute support translation from its initial position [m].
    /// The original finite hair span independently limits the material crossing.
    pub maximum_abs_translation_m: f64,
    /// Shared storage/work bound on authored velocity knots.
    pub max_knots: usize,
    /// Original complete contact-map response rebuild budget per physical step.
    pub max_update_terms: usize,
}

fn poll(gate: &CancelGate) -> Result<(), RenderError> {
    if gate.is_requested() {
        Err(RenderError::Coupled(ModalCouplingError::Cancelled))
    } else {
        Ok(())
    }
}
impl FiberSupportProgram {
    /// Admit every segment, including its interior turning point and the held
    /// final velocity. This can run before source-pack I/O.
    pub fn validate(&self, samples: u64, rate: u32, gate: &CancelGate) -> Result<(), RenderError> {
        poll(gate)?;
        if samples == 0
            || samples > (1_u64 << 53)
            || rate == 0
            || self.knots.is_empty()
            || self.knots.len() > self.max_knots
            || self.max_update_terms == 0
            || [
                self.maximum_speed_m_s,
                self.maximum_acceleration_m_s2,
                self.maximum_abs_translation_m,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
            || !self.initial_velocity_m_s.is_finite()
            || self.initial_velocity_m_s.abs() > self.maximum_speed_m_s
        {
            return Err(render_invalid(
                "fiber support needs bounded velocity/acceleration/translation, a finite clock and nonempty knots",
            ));
        }
        let (mut start, mut velocity, mut position) = (0, self.initial_velocity_m_s, 0.0);
        for knot in &self.knots {
            poll(gate)?;
            if knot.sample <= start
                || knot.sample > samples
                || !knot.velocity_m_s.is_finite()
                || knot.velocity_m_s.abs() > self.maximum_speed_m_s
            {
                return Err(render_invalid(
                    "fiber support knots must be ordered finite velocities inside the performance",
                ));
            }
            let dt = (knot.sample - start) as f64 / f64::from(rate);
            let acceleration = (knot.velocity_m_s - velocity) / dt;
            if !acceleration.is_finite() || acceleration.abs() > self.maximum_acceleration_m_s2 {
                return Err(render_invalid(
                    "fiber support acceleration exceeds its declared bound",
                ));
            }
            // A reversal can exceed the position bound BETWEEN the endpoints.
            if velocity * knot.velocity_m_s < 0.0 {
                let turning = position - 0.5 * velocity * velocity / acceleration;
                if !turning.is_finite() || turning.abs() > self.maximum_abs_translation_m {
                    return Err(render_invalid(
                        "fiber support turning point exceeds the local translation bound",
                    ));
                }
            }
            position += 0.5 * (velocity + knot.velocity_m_s) * dt;
            if !position.is_finite() || position.abs() > self.maximum_abs_translation_m {
                return Err(render_invalid(
                    "fiber support exceeds the local translation bound",
                ));
            }
            (start, velocity) = (knot.sample, knot.velocity_m_s);
        }
        position += velocity * (samples - start) as f64 / f64::from(rate);
        if !position.is_finite() || position.abs() > self.maximum_abs_translation_m {
            return Err(render_invalid(
                "held final fiber support speed exceeds the local translation bound",
            ));
        }
        Ok(())
    }
}

pub(super) struct BoundFiberSupport {
    source: FiberSupportProgram,
    next: usize,
    pending_next: usize,
    pub(super) drives: Vec<FrictionDrive>,
}
impl BoundFiberSupport {
    pub(super) fn prepare(
        &mut self,
        endpoint: u64,
        drives: &[FrictionDrive],
    ) -> Result<(), RenderError> {
        let mut next = self.next;
        while self
            .source
            .knots
            .get(next)
            .is_some_and(|k| k.sample < endpoint)
        {
            next += 1;
        }
        let velocity = if let Some(right) = self.source.knots.get(next) {
            if endpoint == right.sample {
                right.velocity_m_s
            } else {
                let left = next.checked_sub(1).map(|i| self.source.knots[i]);
                let start = left.map_or(0, |k| k.sample);
                let initial = left.map_or(self.source.initial_velocity_m_s, |k| k.velocity_m_s);
                initial
                    + (right.velocity_m_s - initial) * (endpoint - start) as f64
                        / (right.sample - start) as f64
            }
        } else {
            self.source.knots[self.source.knots.len() - 1].velocity_m_s
        };
        if drives.len() != self.drives.len()
            || drives.iter().any(|d| d.speed_m_s != 0.0)
            || !velocity.is_finite()
        {
            return Err(render_invalid(
                "fiber support program requires complete drives and no competing tangential speed",
            ));
        }
        self.drives.copy_from_slice(drives);
        for drive in &mut self.drives {
            drive.speed_m_s = velocity;
        }
        self.pending_next = next;
        Ok(())
    }
    pub(super) fn commit(&mut self) {
        self.next = self.pending_next;
    }
}

impl MultiBowPerformance {
    /// Supply continuous attacks and reversals to the original moving-fiber owner.
    ///
    /// This explicit velocity program owns tangential support motion. All bow
    /// gesture speeds must be zero; their normal loads still act. Endpoint speed
    /// enters the SAME source-derived support inertia and migrating hair contact
    /// used by `with_fiber_support_motion`. Source-projected normal contacts
    /// follow that crossing and supply solved loads. Shear/thermal material
    /// history transport remains unsupported. Bind after optional bow-station motion,
    /// before absolute sample zero; the two support controls are exclusive.
    pub fn with_fiber_support_program(
        mut self,
        program: FiberSupportProgram,
        gate: &CancelGate,
    ) -> Result<Self, RenderError> {
        poll(gate)?;
        self.validate_sample_count(0)?;
        if self.bank.system.samples_rendered() != 0
            || self.fiber_support.is_some()
            || self.bank.fiber_support_motion_limits().is_some()
            || self.events.iter().any(|e| e.drive.speed_m_s != 0.0)
        {
            return Err(render_invalid(
                "fiber support motion binds before playback with zero competing gesture speeds",
            ));
        }
        program.validate(
            self.config.samples,
            self.bank.source.bank().options().sample_rate_hz,
            gate,
        )?;
        let maximum_step_m = program.maximum_speed_m_s * self.bank.system.sample_period_s();
        self.bank = self.bank.with_fiber_support_motion(
            program.initial_velocity_m_s,
            FiberSupportMotionLimits {
                maximum_speed_m_s: program.maximum_speed_m_s,
                maximum_acceleration_m_s2: program.maximum_acceleration_m_s2,
                maximum_step_m,
                max_update_terms: program.max_update_terms,
            },
            gate,
        )?;
        let bound = BoundFiberSupport {
            source: program,
            next: 0,
            pending_next: 0,
            drives: self.candidate.clone(),
        };
        poll(gate)?;
        self.fiber_support = Some(bound);
        Ok(self)
    }
    /// Original physical support trajectory and its admission bounds.
    pub fn fiber_support_program(&self) -> Option<&FiberSupportProgram> {
        self.fiber_support.as_ref().map(|s| &s.source)
    }
    /// Last accepted support state and laboratory mechanical/work conversion.
    pub fn fiber_support_frame(&self) -> Option<&FiberSupportFrame> {
        self.bank.fiber_support_frame()
    }
}
