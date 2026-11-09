//! A prescribed physical clock on the original proportional loading path.
//!
//! Lambda varies affinely between two caller-declared times and loads. The
//! original loading owner alone admits the segment and resolves its first
//! boundary. Its accepted fraction locates that boundary on the same clock.
//! This does not integrate inertia, supply an impulse, release a fold, or
//! establish a quasistatic timescale from material data. A downstream owner
//! must retain the clock/history and explicitly authorize any separation.

use super::{
    JkrCombinedLoading, JkrCombinedLoadingError, JkrCombinedLoadingState, JkrCombinedLoadingStep,
    dd, finite, poll,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// Source-authenticated accepted prefix of an explicitly prescribed affine
/// ramp. The untraversed tail is a prescription, never an admitted contact path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrCombinedLoadingRamp {
    step: JkrCombinedLoadingStep,
    initial_time_s: f64,
    requested_final_time_s: f64,
    accepted_final_time_s: f64,
    boundary_time_s: Option<f64>,
    clock_residual_s: f64,
    maximum_clock_error_s: f64,
}
impl JkrCombinedLoadingRamp {
    /// The exact original step, including its private accepted history.
    pub const fn step(&self) -> &JkrCombinedLoadingStep {
        &self.step
    }
    pub const fn previous_state(&self) -> JkrCombinedLoadingState {
        self.step.previous_state()
    }
    pub const fn state(&self) -> JkrCombinedLoadingState {
        self.step.state()
    }
    pub const fn initial_time_s(&self) -> f64 {
        self.initial_time_s
    }
    pub const fn requested_final_time_s(&self) -> f64 {
        self.requested_final_time_s
    }
    pub const fn accepted_final_time_s(&self) -> f64 {
        self.accepted_final_time_s
    }
    /// A new arrival, including an outward attempt from an existing boundary.
    /// Holding the same load at a boundary does not create another arrival.
    pub const fn boundary_time_s(&self) -> Option<f64> {
        self.boundary_time_s
    }
    pub fn untraversed_time_s(&self) -> f64 {
        (dd(self.requested_final_time_s) - dd(self.accepted_final_time_s)).to_f64()
    }
    /// Stored boundary time minus its compensated affine-clock evaluation [s].
    /// This is a floating-point discrepancy, not an interval time certificate.
    pub const fn clock_residual_s(&self) -> f64 {
        self.clock_residual_s
    }
    pub const fn maximum_clock_error_s(&self) -> f64 {
        self.maximum_clock_error_s
    }

    /// Sample the requested affine prescription at an actual represented time.
    ///
    /// Both requested endpoints and the accepted source boundary retain their
    /// exact loads. Interior times use compensated original affine arithmetic;
    /// the boundary endpoint's clock discrepancy is separately exposed/gated.
    /// Samples after the accepted prefix do not authorize attachment. Sampling
    /// is constant-cost and does not query or mutate the constitutive source.
    pub fn load_factor_at_time(&self, time_s: f64) -> Result<f64, JkrCombinedLoadingError> {
        if !time_s.is_finite()
            || time_s < self.initial_time_s
            || time_s > self.requested_final_time_s
        {
            return Err(JkrCombinedLoadingError::Invalid(
                "sample time outside the prescribed loading interval",
            ));
        }
        if time_s == self.initial_time_s {
            return Ok(self.previous_state().load_factor());
        }
        if time_s == self.accepted_final_time_s {
            return Ok(self.state().load_factor());
        }
        if time_s == self.requested_final_time_s {
            return Ok(self.step.report().requested_load_factor);
        }
        let fraction = checked_dd(
            (dd(time_s) - dd(self.initial_time_s))
                / (dd(self.requested_final_time_s) - dd(self.initial_time_s)),
        )?;
        let initial = dd(self.previous_state().load_factor());
        let requested = dd(self.step.report().requested_load_factor);
        finite(
            checked_dd(initial + fraction * (requested - initial))?.to_f64(),
            "finite affine prescribed load factor",
        )
    }
}

impl JkrCombinedLoading<'_> {
    /// Traverse the original path with an explicitly prescribed affine clock.
    ///
    /// The caller owns chronological continuation: these finite nonnegative
    /// absolute times must match its accepted physical clock. The requested
    /// interval is strictly positive. All original source bounds, work budgets,
    /// reversible stores and detached history remain authoritative.
    ///
    /// A clipped interior load ends at its actual affine boundary time. If that
    /// time collapses to either requested endpoint in floating point, refuse
    /// instead of shifting the event. Otherwise its compensated rounding
    /// discrepancy must fit the caller's explicit absolute seconds tolerance.
    /// A request directed outward from an existing boundary stops at the
    /// initial time; a held boundary consumes the full requested held interval.
    ///
    /// The returned original step never releases energy by itself. Explicit
    /// `detach_at_fold` and any detached continuation belong to the downstream
    /// physical owner. No native motion, heat or contact conductance is inferred.
    pub fn advance_linear_ramp(
        &self,
        previous: JkrCombinedLoadingState,
        initial_time_s: f64,
        requested_final_time_s: f64,
        requested_load_factor: f64,
        maximum_clock_error_s: f64,
        gate: &CancelGate,
    ) -> Result<JkrCombinedLoadingRamp, JkrCombinedLoadingError> {
        self.admit(previous, gate)?;
        if !initial_time_s.is_finite()
            || initial_time_s < 0.0
            || !requested_final_time_s.is_finite()
            || requested_final_time_s <= initial_time_s
            || !maximum_clock_error_s.is_finite()
            || maximum_clock_error_s < 0.0
        {
            return Err(JkrCombinedLoadingError::Invalid(
                "finite ordered loading times and nonnegative clock tolerance required",
            ));
        }
        let step = self.advance_to_load_factor(previous, requested_load_factor, gate)?;
        let initial = previous.load_factor();
        let accepted = step.state().load_factor();
        let moved = requested_load_factor != initial;
        let arrived = moved && step.report().boundary.is_some();
        let (accepted_final_time_s, clock_residual_s) = if accepted == requested_load_factor {
            (requested_final_time_s, 0.0)
        } else if arrived && accepted == initial {
            (initial_time_s, 0.0)
        } else {
            if !arrived {
                return Err(JkrCombinedLoadingError::Invalid(
                    "clipped loading needs the original first boundary",
                ));
            }
            let fraction = checked_dd(
                (dd(accepted) - dd(initial)) / (dd(requested_load_factor) - dd(initial)),
            )?;
            if !Dd::ZERO.lt(fraction) || !fraction.lt(Dd::ONE) {
                return Err(JkrCombinedLoadingError::Invalid(
                    "unrepresentable interior loading fraction",
                ));
            }
            let affine_time = checked_dd(
                dd(initial_time_s) + fraction * (dd(requested_final_time_s) - dd(initial_time_s)),
            )?;
            let time = affine_time.to_f64();
            if time <= initial_time_s || time >= requested_final_time_s {
                return Err(JkrCombinedLoadingError::Invalid(
                    "interior loading boundary collapses onto an interval endpoint",
                ));
            }
            let residual = checked_dd(dd(time) - affine_time)?;
            if dd(maximum_clock_error_s).lt(residual.abs()) {
                return Err(JkrCombinedLoadingError::ClockRoundoff {
                    residual_s: residual.to_f64(),
                    tolerance_s: maximum_clock_error_s,
                });
            }
            (time, residual.to_f64())
        };
        poll(gate)?;
        Ok(JkrCombinedLoadingRamp {
            step,
            initial_time_s,
            requested_final_time_s,
            accepted_final_time_s,
            boundary_time_s: arrived.then_some(accepted_final_time_s),
            clock_residual_s,
            maximum_clock_error_s,
        })
    }
}

fn checked_dd(value: Dd) -> Result<Dd, JkrCombinedLoadingError> {
    if value.hi.is_finite() && value.lo.is_finite() {
        Ok(value)
    } else {
        Err(JkrCombinedLoadingError::Invalid(
            "unrepresentable affine loading clock arithmetic",
        ))
    }
}
