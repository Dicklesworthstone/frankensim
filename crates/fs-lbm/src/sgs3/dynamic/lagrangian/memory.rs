//! Bounded invariant timescales for the original signed momentum history.

use super::{DynamicSgsError, HistoryUpdate, det};

/// Opt-in local momentum-memory law, based on Meneveau, Lund & Cabot (1996),
/// doi:10.1017/S0022112096007379:
/// T = theta Delta (I_LM I_MM / rho^4)^(-1/8).
///
/// Here I is the ADVECTED, accepted history. The density factor removes the
/// squared-density units of this implementation's Favre stress contractions;
/// this variable-density extension is not an incompressible calibration claim.
/// Use the original grid-filter width, not the wider Germano test filter.
///
/// Clamp T between minimum_time and the step config's memory_time. Zero,
/// undetermined or nonpositive signed numerator history uses that upper bound.
/// The initial fit retains the original local spatial average. Both bounds
/// and theta are caller-authored; none is a hidden model default. This bounded
/// signed-history policy differs from the paper's positivity mechanism.
///
/// T is frozen over each current first-order pathline update. This does not
/// add a backscatter-energy model, adaptive scalar memory, or DNS qualification.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveLagrangianMemory {
    /// Positive dimensionless timescale coefficient theta.
    pub coefficient: f64,
    /// Positive lower time bound, in the same units as sample/memory time.
    pub minimum_time: f64,
}

impl AdaptiveLagrangianMemory {
    pub fn validate(self) -> Result<(), DynamicSgsError> {
        if !self.coefficient.is_finite()
            || self.coefficient <= 0.0
            || !self.minimum_time.is_finite()
            || self.minimum_time <= 0.0
        {
            return Err(DynamicSgsError::Invalid(
                "positive finite adaptive Lagrangian coefficient and minimum time required",
            ));
        }
        Ok(())
    }

    pub(super) fn validate_upper_bound(self, maximum_time: f64) -> Result<(), DynamicSgsError> {
        self.validate()?;
        if !maximum_time.is_finite() || maximum_time < self.minimum_time {
            return Err(DynamicSgsError::Invalid(
                "adaptive Lagrangian minimum time exceeds the authored upper time bound",
            ));
        }
        Ok(())
    }
}

/// Only fixed-size per-preparation metadata; no per-cell timescale field.
pub(super) struct AdaptiveUpdate {
    config: AdaptiveLagrangianMemory,
    maximum_time: f64,
    log_length_coefficient: f64,
    log_minimum_time: f64,
    log_maximum_time: f64,
    previous_log: Option<f64>,
    current_log: Option<f64>,
    common_log: f64,
}

impl AdaptiveUpdate {
    pub(super) fn new(
        config: AdaptiveLagrangianMemory,
        maximum_time: f64,
        filter_width: f64,
        previous_log: Option<f64>,
        current_log: Option<f64>,
        update: &HistoryUpdate,
    ) -> Self {
        // Every local T is at most maximum_time. Absorb at least that much
        // decay into the shared scale BEFORE exponentiation. Strong old
        // turbulence can then decay into a weak new state without either the
        // huge original physical contractions or exp(-dt/T) underflowing.
        let old_bound = previous_log
            .map(|value| value - update.elapsed / maximum_time)
            .filter(|value| value.is_finite());
        let common_log = match (old_bound, current_log) {
            (Some(a), Some(b)) => a.max(b),
            (Some(a), None) | (None, Some(a)) => a,
            (None, None) => 0.0,
        };
        Self {
            config,
            maximum_time,
            log_length_coefficient: det::ln(config.coefficient) + det::ln(filter_width),
            log_minimum_time: det::ln(config.minimum_time),
            log_maximum_time: det::ln(maximum_time),
            previous_log,
            current_log,
            common_log,
        }
    }

    pub(super) const fn common_log(&self) -> f64 {
        self.common_log
    }

    fn memory_time(&self, advected: [f64; 2], density: f64) -> f64 {
        let Some(previous_log) = self.previous_log else {
            return self.maximum_time;
        };
        if advected[0] <= 0.0 || advected[1] <= 0.0 {
            return self.maximum_time;
        }
        // No product, fourth power or eighth root is formed. The PREVIOUS
        // scale belongs to these advected values; the current/common scale
        // would give a different, nonphysical memory when amplitudes change.
        let log_time = self.log_length_coefficient
            - 0.125
                * (det::ln(advected[0]) + det::ln(advected[1]) + 2.0 * previous_log
                    - 4.0 * det::ln(density));
        if log_time <= self.log_minimum_time {
            self.config.minimum_time
        } else if log_time >= self.log_maximum_time {
            self.maximum_time
        } else {
            det::exp(log_time)
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn pair(
        &self,
        update: &HistoryUpdate,
        cell: usize,
        velocity: [f64; 3],
        density: f64,
        history: &[[f64; 2]],
        contraction: &[[f64; 2]],
    ) -> Result<([f64; 2], f64), DynamicSgsError> {
        if !update.initialized {
            return Ok((
                update.pair(cell, velocity, history, contraction)?,
                self.maximum_time,
            ));
        }
        let advected = update.advected(cell, velocity, history)?;
        let memory_time = self.memory_time(advected, density);
        let relaxation = update.elapsed / memory_time;
        let fresh = -det::expm1(-relaxation);
        if fresh == 0.0 {
            return Err(DynamicSgsError::Invalid(
                "unrepresentable adaptive Lagrangian relaxation interval",
            ));
        }
        let old_log = self
            .previous_log
            .map(|value| value - relaxation)
            .filter(|value| value.is_finite());
        let old = old_log.map_or(0.0, |value| det::exp(value - self.common_log));
        let new = self
            .current_log
            .map_or(0.0, |value| det::exp(value - self.common_log));
        let current = contraction[cell];
        if (old_log.is_some() && old == 0.0 && advected[1] != 0.0)
            || (self.current_log.is_some() && new == 0.0 && current[1] != 0.0)
        {
            return Err(DynamicSgsError::Cell {
                cell,
                reason: "unrepresentable adaptive Lagrangian history scale separation",
            });
        }
        let pair = core::array::from_fn(|k| old * advected[k] + fresh * (new * current[k]));
        if pair.iter().any(|value| !value.is_finite())
            || pair[1] < 0.0
            || (pair[1] == 0.0
                && (pair[0] != 0.0
                    || (old > 0.0 && advected[1] != 0.0)
                    || (fresh > 0.0 && new > 0.0 && current[1] != 0.0)))
        {
            return Err(DynamicSgsError::Cell {
                cell,
                reason: "unrepresentable adaptive Lagrangian contraction",
            });
        }
        Ok((pair, memory_time))
    }
}

#[cfg(test)]
#[path = "memory/tests.rs"]
mod tests;
