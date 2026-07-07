//! Time histories: typed signals with units and interpolation CONTRACTS
//! (a table is meaningless until its interpolation rule is declared).
//! Every signal knows its dimensions; evaluation is deterministic.

use crate::ScenarioError;
use crate::scenario::Violation;
use fs_cheb::Cheb1;
use fs_qty::{Dims, QtyAny};

/// Interpolation contract for tabulated signals.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Interp {
    /// Piecewise-linear between samples.
    Linear,
    /// Previous-sample hold (step schedules).
    Hold,
}

/// A chebfun-backed profile/history with declared dimensions.
#[derive(Debug, Clone)]
pub struct ChebProfile {
    /// The function object (domain in the signal's independent variable).
    pub cheb: Cheb1,
    /// Dimensions of the VALUES the function produces.
    pub dims: Dims,
}

impl PartialEq for ChebProfile {
    fn eq(&self, other: &Self) -> bool {
        self.dims == other.dims
            && self.cheb.domain() == other.cheb.domain()
            && self.cheb.coeffs() == other.cheb.coeffs()
    }
}

/// A typed time history. Times are seconds (SI); values carry `Dims`.
#[derive(Debug, Clone, PartialEq)]
pub enum TimeSignal {
    /// Constant in time.
    Constant(QtyAny),
    /// Linear ramp, clamped outside `[t_start, t_end]` — the vessel tilt
    /// schedule `(ramp 0deg 65deg 3s)` is exactly this.
    Ramp {
        /// Ramp start (s).
        t_start: f64,
        /// Ramp end (s).
        t_end: f64,
        /// Value at and before `t_start`.
        from: QtyAny,
        /// Value at and after `t_end`.
        to: QtyAny,
    },
    /// Recorded/tabulated trace with a declared interpolation contract,
    /// clamped at the ends.
    Table {
        /// Sample times (s), strictly increasing.
        times: Vec<f64>,
        /// Sample values in SI base units.
        values: Vec<f64>,
        /// Dimensions of the values.
        dims: Dims,
        /// The interpolation contract.
        interp: Interp,
    },
    /// Spectrally-represented smooth history (chebfun function object).
    Chebfun(ChebProfile),
}

impl TimeSignal {
    /// The dimensions this signal produces.
    #[must_use]
    pub fn dims(&self) -> Dims {
        match self {
            TimeSignal::Constant(q) => q.dims,
            TimeSignal::Ramp { from, .. } => from.dims,
            TimeSignal::Table { dims, .. } | TimeSignal::Chebfun(ChebProfile { dims, .. }) => *dims,
        }
    }

    /// Evaluate at time `t` (seconds).
    ///
    /// # Errors
    /// [`ScenarioError::Evaluate`] for structurally bad signals (empty
    /// table, non-finite time).
    pub fn eval(&self, t: f64) -> Result<QtyAny, ScenarioError> {
        if !t.is_finite() {
            return Err(ScenarioError::Evaluate {
                what: format!("non-finite evaluation time {t}"),
            });
        }
        match self {
            TimeSignal::Constant(q) => Ok(*q),
            TimeSignal::Ramp {
                t_start,
                t_end,
                from,
                to,
            } => {
                let clamped = t.clamp(*t_start, *t_end);
                let alpha = if t_end > t_start {
                    (clamped - t_start) / (t_end - t_start)
                } else {
                    1.0
                };
                Ok(QtyAny::new(
                    from.value + alpha * (to.value - from.value),
                    from.dims,
                ))
            }
            TimeSignal::Table {
                times,
                values,
                dims,
                interp,
            } => {
                if times.is_empty() || times.len() != values.len() {
                    return Err(ScenarioError::Evaluate {
                        what: "table signal empty or length-mismatched".to_string(),
                    });
                }
                let v = match times.binary_search_by(|probe| probe.total_cmp(&t)) {
                    Ok(i) => values[i],
                    Err(0) => values[0],
                    Err(i) if i >= times.len() => values[values.len() - 1],
                    Err(i) => match interp {
                        Interp::Hold => values[i - 1],
                        Interp::Linear => {
                            let alpha = (t - times[i - 1]) / (times[i] - times[i - 1]);
                            values[i - 1] + alpha * (values[i] - values[i - 1])
                        }
                    },
                };
                Ok(QtyAny::new(v, *dims))
            }
            TimeSignal::Chebfun(profile) => {
                let (a, b) = profile.cheb.domain();
                Ok(QtyAny::new(profile.cheb.eval(t.clamp(a, b)), profile.dims))
            }
        }
    }

    /// Structural validation, accumulated as [`Violation`]s.
    pub fn check(&self, context: &str, out: &mut Vec<Violation>) {
        match self {
            TimeSignal::Constant(_) => {}
            TimeSignal::Ramp {
                t_start,
                t_end,
                from,
                to,
            } => {
                if t_end < t_start || !t_start.is_finite() || !t_end.is_finite() {
                    out.push(Violation {
                        code: "signal-ramp-times",
                        what: format!("{context}: ramp interval [{t_start}, {t_end}] is invalid"),
                        fix: "order the ramp times so t_start <= t_end and both are finite"
                            .to_string(),
                    });
                }
                if from.dims != to.dims {
                    out.push(Violation {
                        code: "signal-ramp-dims",
                        what: format!(
                            "{context}: ramp endpoints have different dimensions ({:?} vs {:?})",
                            from.dims.0, to.dims.0
                        ),
                        fix: "give both ramp endpoints the same physical dimensions".to_string(),
                    });
                }
            }
            TimeSignal::Table { times, values, .. } => {
                if times.is_empty() || times.len() != values.len() {
                    out.push(Violation {
                        code: "signal-table-shape",
                        what: format!(
                            "{context}: table has {} times and {} values",
                            times.len(),
                            values.len()
                        ),
                        fix: "supply one value per sample time (at least one sample)".to_string(),
                    });
                }
                if times.windows(2).any(|w| w[1] <= w[0]) {
                    out.push(Violation {
                        code: "signal-table-order",
                        what: format!("{context}: table times are not strictly increasing"),
                        fix: "sort the samples by time and remove duplicates".to_string(),
                    });
                }
            }
            TimeSignal::Chebfun(profile) => {
                let (a, b) = profile.cheb.domain();
                let ordered = a.is_finite() && b.is_finite() && a < b;
                if !ordered {
                    out.push(Violation {
                        code: "signal-chebfun-domain",
                        what: format!("{context}: chebfun domain [{a}, {b}] is empty"),
                        fix: "build the function object on a nonempty time interval".to_string(),
                    });
                }
            }
        }
    }
}
