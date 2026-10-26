//! Normal translation of an existing compliant support, not a new contact law.
//!
//! The original attachment closure becomes x = N^T q - s. The same held
//! compressive reaction applies -N R to the mechanics and supplies -R delta_s
//! of actuator work. The original potential secant/loss is evaluated on delta_x,
//! not on string/body motion alone. No source gap is mutated and no startup
//! overlap, mechanical motion or interface history is discarded.

use super::*;

/// Physical limits for all driven supports, checked on each candidate endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObstacleMotionLimits {
    /// Maximum absolute translation from the authored reference geometry [m].
    pub maximum_abs_offset_m: f64,
    /// Maximum average translation speed over one mechanical step [m/s].
    /// Endpoint positions define a linear path within that step; this is not
    /// a bound on an unspecified, more oscillatory sub-sample trajectory.
    pub maximum_speed_m_s: f64,
}
impl ObstacleMotionLimits {
    /// Admit positive finite physical limits; no clipping or normalization.
    pub fn validate(self) -> Result<(), ModalCouplingError> {
        if !self.maximum_abs_offset_m.is_finite() || self.maximum_abs_offset_m <= 0.0
            || !self.maximum_speed_m_s.is_finite() || self.maximum_speed_m_s <= 0.0 {
            return Err(invalid("normal motion requires positive finite translation and speed limits"));
        }
        Ok(())
    }
}

pub(super) struct MotionState {
    pub(super) accepted: Vec<f64>,
    pub(super) candidate: Vec<f64>,
    limits: ObstacleMotionLimits,
}

impl MultiFrictionModalSystem {
    /// Admit translation along each existing obstacle's fixed normal coordinate.
    ///
    /// Bind once before any physical sample, after the normal obstacles exist.
    /// Offset zero is the original authored clearance, including any original
    /// overlap. Existing modal motion and contact preload are retained exactly.
    /// This does not translate the contact station along a string, rotate its
    /// normal, solve an actuator, or replace prescribed bow normal loads.
    ///
    /// A moving support on a receiver adds s to that receiver's existing normal
    /// displacement. Both original attachment bodies still receive reciprocal
    /// reactions. The extra actuator work belongs only to the prescribed s.
    pub fn with_obstacle_motion(mut self, limits: ObstacleMotionLimits,
        gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        poll(Some(gate))?;
        limits.validate()?;
        if self.samples_rendered() != 0 {
            return Err(invalid("normal motion must bind before the first physical sample").into());
        }
        let obstacles = self.obstacles.as_mut()
            .ok_or_else(|| invalid("normal motion requires existing compliant obstacles"))?;
        if obstacles.motion.is_some() {
            return Err(invalid("normal motion cannot replace an existing binding").into());
        }
        if obstacles.has_residual() {
            return Err(invalid("normal motion does not carry the obstacles' residual compliance").into());
        }
        let count = obstacles.points.len();
        obstacles.motion = Some(MotionState {
            accepted: vec![0.0; count], candidate: vec![0.0; count], limits,
        });
        Ok(self)
    }

    /// Last accepted support offsets [m], in the original obstacle order.
    /// An error or cancelled trial cannot advance these positions.
    #[must_use]
    pub fn obstacle_offsets_m(&self) -> Option<&[f64]> {
        Some(&self.obstacles.as_ref()?.motion.as_ref()?.accepted)
    }

    /// Physical translation/speed limits retained with this realization.
    #[must_use]
    pub fn obstacle_motion_limits(&self) -> Option<ObstacleMotionLimits> {
        Some(self.obstacles.as_ref()?.motion.as_ref()?.limits)
    }

    /// Advance one sample with one final normal offset per attached obstacle.
    ///
    /// Positive offset opens the gap; negative offset presses into the bodies.
    /// Offsets are absolute translations relative to the original source, NOT
    /// per-step increments or new source gaps. Existing `step` calls hold the
    /// accepted offsets, so omitting a new command never snaps supports to zero.
    ///
    /// All bow and obstacle reactions are solved jointly. Normal actuator work
    /// and contact storage/loss use the same endpoints before state, pressure,
    /// contact reports, offsets and clock are published together. Refusal is
    /// atomic. Linear paths are the declared per-step motion; contact resolution
    /// and time refinement still matter and no rigid-impact limit is claimed.
    pub fn step_with_obstacle_offsets(&mut self, external: &[f64], drives: &[FrictionDrive],
        end_offsets_m: &[f64]) -> Result<&MultiFrictionFrame, ModalFrictionError>
    {
        self.step_inner(external, drives, Some(end_offsets_m), None, None)
    }

    /// Same complete physical transaction with cancellation polling.
    pub fn step_with_obstacle_offsets_under_gate(&mut self, external: &[f64], drives: &[FrictionDrive],
        end_offsets_m: &[f64], gate: &CancelGate) -> Result<&MultiFrictionFrame, ModalFrictionError>
    {
        self.step_inner(external, drives, Some(end_offsets_m), None, Some(gate))
    }
}

impl ObstacleSet {
    pub(in super::super) fn check_motion(&self, offsets: &[f64], dt: f64)
        -> Result<(), ModalCouplingError>
    {
        let motion = self.motion.as_ref()
            .ok_or_else(|| invalid("normal translation has not been admitted"))?;
        if offsets.len() != self.points.len() {
            return Err(invalid("normal motion needs one endpoint for every attached obstacle"));
        }
        for (&end, &start) in offsets.iter().zip(&motion.accepted) {
            limit("normal support offset", finite(end)?.abs(), motion.limits.maximum_abs_offset_m)?;
            let speed = finite(finite(end-start)?/dt)?;
            limit("normal support speed", speed.abs(), motion.limits.maximum_speed_m_s)?;
        }
        Ok(())
    }
}

