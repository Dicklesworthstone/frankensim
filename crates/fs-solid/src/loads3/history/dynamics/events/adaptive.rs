//! Local step-doubling estimates over event-resolved original Newmark steps.
//! Publish the two-half-step trajectory, never Richardson-extrapolated motion.
//! Agreement of finite samples can alias unresolved smooth forcing: the caller
//! must still set a maximum step that resolves its bandwidth. No global error,
//! continuum accuracy, or certificate is inferred from this local estimate.
use super::*;

/// Componentwise SI absolute tolerances and a change-relative tolerance.
/// The minimum applies to retried PAIR spans, not their individual half steps.
/// A target/event may initially shorten a pair below that retry floor.
#[derive(Debug, Clone, Copy)]
pub struct TransientTimeControl3 {
    pub absolute_displacement_m: f64,
    pub absolute_velocity_m_s: f64,
    pub absolute_acceleration_m_s2: f64,
    pub relative: f64,
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    /// Includes the initial attempt; capped at 32 (at most 96 native intervals).
    pub maximum_attempts: usize,
}
impl Default for TransientTimeControl3 {
    fn default() -> Self {
        Self { absolute_displacement_m: 1e-8, absolute_velocity_m_s: 1e-6,
            absolute_acceleration_m_s2: 1e-4, relative: 1e-4,
            minimum_step_s: 1e-9, maximum_step_s: 0.01, maximum_attempts: 16 }
    }
}
impl TransientTimeControl3 {
    fn validate(self) -> Result<(), AdaptiveHistoryError3> {
        if [self.absolute_displacement_m, self.absolute_velocity_m_s,
            self.absolute_acceleration_m_s2, self.minimum_step_s, self.maximum_step_s]
            .iter().any(|v| !v.is_finite() || *v <= 0.0)
            || !self.relative.is_finite() || !(0.0..1.0).contains(&self.relative)
            || self.minimum_step_s > self.maximum_step_s
            || !(1..=32).contains(&self.maximum_attempts) {
            return Err(AdaptiveHistoryError3::Invalid("invalid local motion-error or retry controls"));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, PartialEq)]
pub enum AdaptiveHistoryError3 {
    History(TransientHistoryError3),
    Invalid(&'static str),
    Rejected { attempts: usize, attempted_step_s: f64, error_ratio: f64 },
}
impl From<TransientHistoryError3> for AdaptiveHistoryError3 {
    fn from(e: TransientHistoryError3) -> Self { Self::History(e) }
}
impl std::fmt::Display for AdaptiveHistoryError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "adaptive tetrahedral history refused: {self:?}")
    }
}
impl std::error::Error for AdaptiveHistoryError3 {}

/// Both accepted physical reports remain available for exact work accounting.
/// The discarded coarse solve contributes no trajectory work or accepted step.
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveScheduledStep3 {
    pub first: ScheduledStep3,
    pub second: ScheduledStep3,
    /// max(|fine-coarse| / (3 * component allowance)); asymptotic local estimate.
    pub error_ratio: f64,
    pub rejected_attempts: usize,
    pub suggested_step_s: f64,
}
impl EventTransientRun3<'_, '_, '_, '_, '_, '_> {
    /// Advance ONE accepted pair, respecting all original event and solve gates.
    /// Compare one full interval with two half intervals from the same immutable
    /// checkpoint. Displacement, velocity AND acceleration must pass. Relative
    /// scaling uses changes from the accepted state, not an arbitrary offset.
    ///
    /// The pair consumes two original lifetime step slots. Failed trials consume
    /// bounded computation but no accepted slots or source work. Budget, source,
    /// energy and linear-solve failures propagate without weakening any gate.
    /// Cancellation or an exhausted time-error allowance publishes neither half.
    /// The full endpoint is never beyond the first authored event. Its right
    /// load jump is applied equally to the coarse and fine comparison states.
    pub fn advance_adaptive_to(&mut self, target_s: f64, proposed_step_s: f64,
        control: TransientTimeControl3, cx: &Cx<'_>)
        -> Result<AdaptiveScheduledStep3, AdaptiveHistoryError3> {
        checkpoint(cx)?; control.validate()?;
        if !proposed_step_s.is_finite() || proposed_step_s <= 0.0 {
            return Err(AdaptiveHistoryError3::Invalid("proposed step must be finite and positive"));
        }
        let start = self.run.state().time_s;
        let mut cap = proposed_step_s.min(control.maximum_step_s).max(control.minimum_step_s);
        for attempt in 1..=control.maximum_attempts {
            checkpoint(cx)?;
            let end = self.schedule.endpoint(start, target_s, cap)?;
            let span = end - start;
            let middle = start + 0.5 * span;
            if middle <= start || middle >= end || !middle.is_finite() {
                return Err(AdaptiveHistoryError3::Invalid("adaptive interval has no distinct representable midpoint"));
            }
            let mut coarse = self.clone();
            let full = coarse.advance_to(end, span, cx)?;
            let mut fine = self.clone();
            let first = fine.advance_to(middle, middle-start, cx)?;
            let second = fine.advance_to(end, end-middle, cx)?;
            if full.end_time_s != end || first.end_time_s != middle || second.end_time_s != end {
                return Err(AdaptiveHistoryError3::Invalid("adaptive comparison failed to reach identical authored endpoints"));
            }
            let error_ratio = discrepancy(self, &coarse, &fine, control, cx)?;
            if error_ratio <= 1.0 {
                let factor = if error_ratio == 0.0 { 2.0 } else { (0.8 / error_ratio.cbrt()).clamp(0.5, 2.0) };
                let suggested_step_s = (span * factor).min(control.maximum_step_s).max(control.minimum_step_s);
                checkpoint(cx)?;
                *self = fine;
                return Ok(AdaptiveScheduledStep3 { first, second, error_ratio,
                    rejected_attempts: attempt-1, suggested_step_s });
            }
            let next = (span * (0.8 / error_ratio.cbrt()).clamp(0.2, 0.8)).max(control.minimum_step_s);
            if attempt == control.maximum_attempts || next >= span || start + next <= start {
                return Err(AdaptiveHistoryError3::Rejected { attempts: attempt,
                    attempted_step_s: span, error_ratio });
            }
            cap = next;
        }
        Err(AdaptiveHistoryError3::Invalid("empty adaptive attempt budget"))
    }
}
fn discrepancy(base: &EventTransientRun3<'_, '_, '_, '_, '_, '_>,
    coarse: &EventTransientRun3<'_, '_, '_, '_, '_, '_>,
    fine: &EventTransientRun3<'_, '_, '_, '_, '_, '_>, control: TransientTimeControl3,
    cx: &Cx<'_>) -> Result<f64, AdaptiveHistoryError3> {
    let (b, c, f) = (base.run.state(), coarse.run.state(), fine.run.state());
    let mut ratio = 0.0_f64;
    for (old, full, half, absolute) in [
        (&b.displacement, &c.displacement, &f.displacement, control.absolute_displacement_m),
        (&b.velocity, &c.velocity, &f.velocity, control.absolute_velocity_m_s),
        (&b.acceleration, &c.acceleration, &f.acceleration, control.absolute_acceleration_m_s2),
    ] {
        for (i, ((o, c), f)) in old.iter().zip(full).zip(half).enumerate() {
            if i % 256 == 0 { checkpoint(cx)?; }
            let scale = (f-o).abs().max((c-o).abs());
            let tolerance = absolute + control.relative * scale;
            let error = ((f-c) / 3.0).abs() / tolerance;
            if !scale.is_finite() || !tolerance.is_finite() || !error.is_finite() {
                return Err(AdaptiveHistoryError3::Invalid("nonfinite motion-error comparison"));
            }
            ratio = ratio.max(error);
        }
    }
    checkpoint(cx)?; Ok(ratio)
}
