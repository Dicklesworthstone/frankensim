//! Typed bridge from an accepted trajectory adjoint to scalar design search.
//! A local derivative suggests a candidate, never a feasibility decision.
use super::*;

#[derive(Debug, Clone, Copy)]
pub(super) struct DesignSensitivity {
    peak_k: f64,
    // At the evaluated candidate, scale its WHOLE actual source by exp(eps).
    log_power_k: f64,
    log_speed_k: Option<f64>,
}

impl DesignSensitivity {
    /// The interval gradients already include every causal repeated occurrence.
    /// Sum those shared controls once, not once again for each cycle.
    pub(super) fn from_intervals(peak_k: f64, powers: &[f64], speeds: &[Option<f64>]) -> Result<Self> {
        if powers.len() != speeds.len() || powers.is_empty() {
            return Err(bad("trajectory design derivative requires one row per interval"));
        }
        let mut log_power_k = 0.0;
        let mut log_speed_k = Some(0.0);
        for (&power, &speed) in powers.iter().zip(speeds) {
            log_power_k = finite(log_power_k + finite(power)?)?;
            if let Some(value) = speed { finite(value)?; }
            log_speed_k = match (log_speed_k, speed) {
                (Some(sum), Some(value)) => Some(finite(sum + value)?),
                _ => None,
            };
        }
        Ok(Self { peak_k: finite(peak_k)?, log_power_k, log_speed_k })
    }

    /// The retained derivative must describe this exact sampled-peak value.
    /// Convert dT/dln(m) to dT/dm at the OUTER multiplier, not at one.
    /// A zero-power candidate has no invertible relative scaling: withhold
    /// the slope there rather than dividing zero by zero or inventing one.
    pub(super) fn derivative(self, peak_k: f64, multiplier: f64, power: bool) -> Result<Option<f64>> {
        if self.peak_k.to_bits() != peak_k.to_bits() {
            return Err(producer("trajectory design derivative does not match the accepted sampled peak"));
        }
        if !(multiplier.is_finite() && multiplier >= 0.0) {
            return Err(bad("invalid trajectory design multiplier"));
        }
        if multiplier == 0.0 { return Ok(None); }
        let relative = if power { Some(self.log_power_k) } else { self.log_speed_k };
        // Overflow of an optional proposal slope is not a failed primal solve.
        // The bounded search still has its derivative-free bisection fallback.
        Ok(relative.map(|value| value / multiplier).filter(|value| value.is_finite()))
    }
}

/// Staying strictly inside the central 80% guarantees a surviving bracket
/// contracts by at least 10%, even at a peak-branch change or with a bad slope.
/// Use the original multiplier coordinate so zero workload remains admissible.
/// The expected sign is a proposal filter, NOT a monotonicity assumption used
/// to certify any unevaluated design.
pub(super) fn newton_proposal(low: f64, high: f64, at: f64, temperature: f64,
    derivative: Option<f64>, limit: f64, power: bool) -> Option<f64> {
    let slope = derivative.filter(|value| value.is_finite()
        && if power { *value > 0.0 } else { *value < 0.0 })?;
    if ![low, high, at, temperature, limit].into_iter().all(f64::is_finite)
        || !(low < high && at >= low && at <= high) { return None; }
    let candidate = at - (temperature - limit) / slope;
    let guard = 0.1 * (high - low);
    (candidate.is_finite() && candidate > low + guard && candidate < high - guard)
        .then_some(candidate)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn candidate_chain_rule_is_not_a_derivative_at_multiplier_one() {
        let gradient = DesignSensitivity::from_intervals(310.0, &[6.0, 2.0], &[Some(-1.0), Some(-0.5)]).unwrap();
        assert_eq!(gradient.derivative(310.0, 2.0, true).unwrap(), Some(4.0));
        assert_eq!(gradient.derivative(310.0, 0.5, false).unwrap(), Some(-3.0));
        assert_eq!(gradient.derivative(310.0, 0.0, true).unwrap(), None);
        assert!(gradient.derivative(310.0_f64.next_up(), 1.0, true).is_err());
    }

    #[test]
    fn unavailable_fan_gradient_does_not_erase_the_power_gradient() {
        let gradient = DesignSensitivity::from_intervals(310.0, &[3.0, 0.0], &[Some(-1.0), None]).unwrap();
        assert_eq!(gradient.derivative(310.0, 1.0, false).unwrap(), None);
        assert_eq!(gradient.derivative(310.0, 1.0, true).unwrap(), Some(3.0));
        assert!(DesignSensitivity::from_intervals(310.0, &[f64::NAN], &[None]).is_err());
        assert!(DesignSensitivity::from_intervals(310.0, &[1.0], &[]).is_err());
    }

    #[test]
    fn unusable_slopes_and_outside_trials_require_bisection() {
        for slope in [None, Some(0.0), Some(1.0), Some(f64::NAN), Some(-1e-300)] {
            assert_eq!(newton_proposal(0.5, 2.0, 0.5, 310.0, slope, 305.0, false), None);
        }
        assert_eq!(newton_proposal(0.0, 2.0, 2.0, 310.0, Some(-2.0), 305.0, true), None);
        assert_eq!(newton_proposal(0.5, 2.0, 2.0, 304.99999, Some(-10.0), 305.0, false), None);
    }

    #[test]
    fn both_design_directions_keep_an_evaluated_bracket() {
        for (at, value, slope, power) in [(0.5, 310.0, -5.0, false), (2.0, 310.0, 5.0, true)] {
            let candidate = newton_proposal(0.5, 2.0, at, value, Some(slope), 305.0, power).unwrap();
            assert!(candidate > 0.65 && candidate < 1.85);
            assert!((candidate - 0.5).max(2.0 - candidate) < 0.9 * 1.5);
        }
    }
}
