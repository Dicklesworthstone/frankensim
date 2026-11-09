//! One scalar history reduction for native momentum and attached storage.

/// The derivative is `(new - history(recent, older))/effective_timestep_s`.
/// Equal-step arithmetic deliberately retains the original momentum bits.
#[derive(Clone, Copy, Debug)]
pub(crate) struct StokesBdfStep {
    pub(crate) timestep_s: f64,
    pub(crate) effective_timestep_s: f64,
    history_weight: f64,
    equal_steps: bool,
}

impl StokesBdfStep {
    pub(crate) fn new(
        timestep_s: f64,
        previous_timestep_s: Option<f64>,
    ) -> Result<Self, &'static str> {
        let valid_interval =
            |dt: f64| dt.is_finite() && dt > 0.0 && dt.recip().is_finite() && dt.recip() > 0.0;
        if !valid_interval(timestep_s) || previous_timestep_s.is_some_and(|dt| !valid_interval(dt))
        {
            return Err("finite history and representable positive timestep required");
        }
        let Some(previous_timestep_s) = previous_timestep_s else {
            return Ok(Self {
                timestep_s,
                effective_timestep_s: timestep_s,
                history_weight: 0.0,
                equal_steps: false,
            });
        };
        let ratio = timestep_s / previous_timestep_s;
        if !ratio.is_finite() || ratio <= 0.0 || ratio > 2.0 {
            return Err("BDF2 needs a representable positive step ratio no greater than two");
        }
        let history_weight = ratio * ratio / (1.0 + 2.0 * ratio);
        if !history_weight.is_finite() || history_weight <= 0.0 {
            return Err("BDF2 history coefficient is not representable");
        }
        let equal_steps = timestep_s == previous_timestep_s;
        let effective_timestep_s = if equal_steps {
            (2.0 / 3.0) * timestep_s
        } else {
            timestep_s * ((1.0 + ratio) / (1.0 + 2.0 * ratio))
        };
        if !valid_interval(effective_timestep_s) {
            return Err("finite history and representable positive timestep required");
        }
        Ok(Self {
            timestep_s,
            effective_timestep_s,
            history_weight,
            equal_steps,
        })
    }

    /// BDF1 ignores `older`. Callers validate the resulting state as well as
    /// the inputs: finite inputs alone do not guarantee representable history.
    pub(crate) fn history(self, recent: f64, older: f64) -> f64 {
        if self.history_weight == 0.0 {
            return recent;
        }
        let difference = recent - older;
        let correction = if self.equal_steps {
            if difference.is_finite() {
                difference / 3.0
            } else {
                recent / 3.0 - older / 3.0
            }
        } else if difference.is_finite() {
            self.history_weight * difference
        } else {
            self.history_weight * recent - self.history_weight * older
        };
        recent + correction
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// G0: unequal-step BDF2 differentiates quadratics exactly.
    #[test]
    fn quadratic_history_and_original_equal_step_bits() {
        for (dt, old_dt) in [(0.1, 0.1), (0.07, 0.1), (0.2, 0.1)] {
            let step = StokesBdfStep::new(dt, Some(old_dt)).unwrap();
            let t: f64 = 0.4;
            let history = step.history(t * t, (t - old_dt).powi(2));
            let derivative = ((t + dt).powi(2) - history) / step.effective_timestep_s;
            assert!((derivative - 2.0 * (t + dt)).abs() < 1e-14);
            assert_eq!(step.history(1e200, 1e200).to_bits(), 1e200_f64.to_bits());
        }
        let equal = StokesBdfStep::new(0.3, Some(0.3)).unwrap();
        assert_eq!(equal.history(1.1, 0.3), 1.1 + (1.1 - 0.3) / 3.0);
        assert_eq!(equal.effective_timestep_s, (2.0 / 3.0) * 0.3);
        let startup = StokesBdfStep::new(0.3, None).unwrap();
        assert_eq!(startup.history(1.1, f64::NAN), 1.1);
    }

    #[test]
    fn invalid_and_unrepresentable_intervals_are_refused() {
        for dt in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
            assert!(StokesBdfStep::new(dt, None).is_err());
        }
        assert!(StokesBdfStep::new(0.21, Some(0.1)).is_err());
        assert!(StokesBdfStep::new(0.1, Some(0.0)).is_err());
        assert!(StokesBdfStep::new(1e-200, Some(1e100)).is_err());
    }
}
