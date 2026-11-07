//! Shared bounded step-doubling schedule; each physical owner retains its maps,
//! candidate buffers, error norm, work ledger and final publication gates.

use super::det;

#[derive(Clone, Copy)]
pub(super) enum Order {
    First,
    Second,
}

#[derive(Clone, Copy)]
pub(super) struct Limits {
    pub initial_step_s: f64,
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    pub maximum_trials: usize,
}

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Statistics {
    pub trials: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub maximum_accepted_normalized_error: f64,
    pub suggested_next_step_s: f64,
}

#[derive(Clone, Copy)]
pub(super) struct Trial {
    end_s: f64,
    pub span_s: f64,
    pub halves_s: [f64; 2],
    pub fine_error_factor: f64,
}

pub(super) enum Error {
    Invalid,
    TrialBudget { completed_trials: usize },
    TemporalAccuracy { step_s: f64, normalized_error: f64 },
    TimeResolution { elapsed_s: f64 },
}

pub(super) struct Schedule {
    duration_s: f64,
    elapsed_s: f64,
    limits: Limits,
    order: Order,
    pub statistics: Statistics,
}

impl Schedule {
    pub fn new(duration_s: f64, limits: Limits, order: Order) -> Result<Self, Error> {
        if !duration_s.is_finite()
            || duration_s < 0.0
            || [
                limits.minimum_step_s,
                limits.initial_step_s,
                limits.maximum_step_s,
            ]
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
            || limits.minimum_step_s > limits.initial_step_s
            || limits.initial_step_s > limits.maximum_step_s
            || limits.maximum_trials == 0
        {
            return Err(Error::Invalid);
        }
        Ok(Self {
            duration_s,
            elapsed_s: 0.0,
            limits,
            order,
            statistics: Statistics {
                suggested_next_step_s: limits.initial_step_s,
                ..Statistics::default()
            },
        })
    }

    pub fn complete(&self) -> bool {
        self.elapsed_s >= self.duration_s
    }

    /// Called before an owner's whole-comparison work reservation. It does not
    /// consume a trial, so an exhausted work allowance has no hidden probe.
    pub fn check_trial_budget(&self) -> Result<(), Error> {
        if self.statistics.trials >= self.limits.maximum_trials {
            Err(Error::TrialBudget {
                completed_trials: self.statistics.trials,
            })
        } else {
            Ok(())
        }
    }

    pub fn begin(&mut self) -> Result<Trial, Error> {
        self.check_trial_budget()?;
        let remaining = self.duration_s - self.elapsed_s;
        let end_s = if self.statistics.suggested_next_step_s >= remaining {
            self.duration_s
        } else {
            self.elapsed_s + self.statistics.suggested_next_step_s
        };
        let span_s = end_s - self.elapsed_s;
        let first = 0.5 * span_s;
        let second = span_s - first;
        if end_s <= self.elapsed_s || end_s > self.duration_s || first <= 0.0 || second <= 0.0 {
            return Err(Error::TimeResolution {
                elapsed_s: self.elapsed_s,
            });
        }
        // Compare the actual representable halves. For order p, the fine
        // truncation estimate is (a^(p+1)+b^(p+1)) divided by
        // ((a+b)^(p+1)-a^(p+1)-b^(p+1)), times the full/fine difference.
        let ratio = first.min(second) / first.max(second);
        let fine_error_factor = match self.order {
            Order::First => (1.0 + ratio * ratio) / (2.0 * ratio),
            Order::Second => (1.0 - ratio + ratio * ratio) / (3.0 * ratio),
        };
        self.statistics.trials += 1;
        Ok(Trial {
            end_s,
            span_s,
            halves_s: [first, second],
            fine_error_factor,
        })
    }

    /// True accepts the unextrapolated fine candidate. No physical state is
    /// owned here; the caller alone commits a successful local candidate.
    pub fn consider(&mut self, trial: Trial, error: f64) -> Result<bool, Error> {
        if error.is_nan() || error < 0.0 {
            return Err(Error::Invalid);
        }
        if error > 1.0 {
            self.statistics.rejected_intervals += 1;
            let smaller =
                (trial.span_s * self.factor(error).min(0.8)).max(self.limits.minimum_step_s);
            if trial.span_s <= self.limits.minimum_step_s
                || smaller >= trial.span_s
                || smaller >= self.statistics.suggested_next_step_s
            {
                return Err(Error::TemporalAccuracy {
                    step_s: trial.span_s,
                    normalized_error: error,
                });
            }
            self.statistics.suggested_next_step_s = smaller;
            return Ok(false);
        }
        self.statistics.accepted_intervals += 1;
        self.statistics.maximum_accepted_normalized_error =
            self.statistics.maximum_accepted_normalized_error.max(error);
        self.elapsed_s = trial.end_s;
        let growth = self.factor(error);
        self.statistics.suggested_next_step_s =
            if trial.span_s >= self.limits.maximum_step_s / growth {
                self.limits.maximum_step_s
            } else {
                (trial.span_s * growth).max(self.limits.minimum_step_s)
            };
        Ok(true)
    }

    fn factor(&self, error: f64) -> f64 {
        if error == 0.0 {
            2.0
        } else if error.is_infinite() {
            0.2
        } else {
            let exponent = match self.order {
                Order::First => -0.5,
                Order::Second => -1.0 / 3.0,
            };
            (0.9 * det::pow(error, exponent)).clamp(0.2, 2.0)
        }
    }
}
