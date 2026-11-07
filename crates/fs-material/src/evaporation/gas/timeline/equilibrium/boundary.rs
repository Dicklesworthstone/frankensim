//! Complete freezing on the original network and finite-gas clock.
//!
//! The parent timeline retries the complete native-heat/phase proposal when a
//! phase event cuts its interval. Discarded network heat never survives with a
//! shorter vapor clock. The first all-solid native substep also bounds the
//! search, including a covered interface or disabled vapor transfer. Source
//! identity, event order, exposure and the original work/balance gates remain
//! owned by the existing timeline. Event localization does not remove the
//! first-order split's need for physical timestep refinement.

use super::{
    ClosedGasParcel, DepositionError, DepositionSimulation, EquilibriumTimelineSurface,
    GasEvaporationSolve, GasEvaporationTimelineCandidate, GasEvaporationTimelineError, StepConfig,
    ThermochemicalGasModel, advance_deposition,
};

/// Time consumed by the whole coupled proposal and time left for a new owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingTimelineBoundary {
    pub elapsed_duration_s: f64,
    pub remaining_duration_s: f64,
    pub freezing_completed: bool,
}

/// Original network, gas and fusion account at one accepted physical time.
#[derive(Clone)]
pub struct EquilibriumFreezingTimelineCandidate<'thermal, 'catalog> {
    exchange: GasEvaporationTimelineCandidate<'thermal, 'catalog>,
    boundary: FreezingTimelineBoundary,
}
impl<'thermal, 'catalog> EquilibriumFreezingTimelineCandidate<'thermal, 'catalog> {
    #[must_use]
    pub const fn exchange(&self) -> &GasEvaporationTimelineCandidate<'thermal, 'catalog> {
        &self.exchange
    }
    #[must_use]
    pub const fn boundary(&self) -> &FreezingTimelineBoundary {
        &self.boundary
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        GasEvaporationTimelineCandidate<'thermal, 'catalog>,
        FreezingTimelineBoundary,
    ) {
        (self.exchange, self.boundary)
    }
}

/// Advance at most to `target_s`, returning early at the original solid knot.
///
/// Initial and scheduled bodies still require a liquid seed. A completed body
/// belongs to an explicitly selected next interface owner; this API does not
/// infer solid coverage or switch continuum-film evaporation to sublimation.
/// Every rejected event trial spends the original native thermal work budget.
/// Exact-knot selection in a heat-only event spends the original thermal joule
/// allowance and appears in its stored-energy change and residual. Active
/// vapor transfer must reach that knot at the SAME duration as network heat.
/// Unrepresentable or unbracketed endpoints refuse without publishing state.
#[allow(clippy::too_many_arguments)]
pub fn advance_deposition_to_freezing<'thermal, 'catalog>(
    thermal: &DepositionSimulation<'thermal>,
    surface: EquilibriumTimelineSurface<'_>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    target_s: f64,
    thermal_config: StepConfig,
    phase_config: GasEvaporationSolve,
    checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<EquilibriumFreezingTimelineCandidate<'thermal, 'catalog>, GasEvaporationTimelineError> {
    let exchange = advance_deposition(
        thermal,
        surface.into(),
        gas,
        catalog,
        target_s,
        thermal_config,
        phase_config,
        None,
        true,
        checkpoint,
    )?;
    let boundary = FreezingTimelineBoundary {
        elapsed_duration_s: exchange.thermal().time_s() - thermal.time_s(),
        remaining_duration_s: target_s - exchange.thermal().time_s(),
        freezing_completed: exchange
            .audit()
            .phase_inventory
            .as_ref()
            .is_some_and(|phase| phase.final_state.liquid_mass_kg == 0.0),
    };
    Ok(EquilibriumFreezingTimelineCandidate { exchange, boundary })
}

// This is only an event bracket for the existing physical owners, not another
// transport solver. Its residuals are source-relative extensive energy gaps.
pub(in crate::evaporation::gas::timeline) struct FreezingSearch {
    lower: (f64, f64),
    upper: Option<(f64, f64)>,
    trials: usize,
}
impl FreezingSearch {
    pub(in crate::evaporation::gas::timeline) fn new(start_s: f64, headroom_j: f64) -> Self {
        Self {
            lower: (start_s, headroom_j),
            upper: None,
            trials: 0,
        }
    }
    pub(in crate::evaporation::gas::timeline) const fn is_bracketed(&self) -> bool {
        self.upper.is_some()
    }
    /// None accepts an ordinary positive-liquid full interval. Once an upper
    /// bound exists, neither a positive trial nor a clipped event is published.
    pub(in crate::evaporation::gas::timeline) fn observe(
        &mut self,
        time_s: f64,
        residual_j: f64,
    ) -> Result<Option<f64>, DepositionError> {
        if !time_s.is_finite() || !residual_j.is_finite() || self.lower.1 <= 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        if residual_j <= 0.0 {
            self.upper = Some((time_s, residual_j));
        } else if self.upper.is_some() {
            self.lower = (time_s, residual_j);
        } else {
            return Ok(None);
        }
        self.trials += 1;
        if self.trials > 64 {
            return Err(DepositionError::Budget(
                "freezing event localization trials",
            ));
        }
        let (hi, r_hi) = self.upper.ok_or(DepositionError::Arithmetic)?;
        let (lo, r_lo) = self.lower;
        if hi <= lo {
            return Err(DepositionError::Arithmetic);
        }
        // Scale before summing to avoid overflowing a perfectly finite bracket.
        // Keep a useful contraction when a nearly flat response defeats secant.
        let scale = r_lo.max(r_hi.abs());
        let fraction = ((r_lo / scale) / (r_lo / scale + r_hi.abs() / scale)).clamp(0.05, 0.95);
        let mut next = lo + (hi - lo) * fraction;
        if next <= lo || next >= hi {
            next = lo + (hi - lo) * 0.5;
        }
        if !next.is_finite() || next <= lo || next >= hi {
            return Err(DepositionError::InvalidInput(
                "freezing event has no representable common heat/phase interval",
            ));
        }
        Ok(Some(next))
    }
}

#[cfg(test)]
mod tests;
