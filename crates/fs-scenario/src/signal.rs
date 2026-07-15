//! Time histories: typed signals with units and interpolation CONTRACTS
//! (a table is meaningless until its interpolation rule is declared).
//! Every signal knows its dimensions; evaluation is deterministic.

use crate::ScenarioError;
use crate::scenario::Violation;
use fs_cheb::Cheb1;
use fs_qty::{Dims, QtyAny};

/// Smooth histories have no finite breakpoint set. Net-flux admission uses
/// this many equal panels (endpoints included) on each function's declared
/// domain. This is a bounded deterministic screen, not a proof between points.
const SMOOTH_NET_FLUX_VALIDATION_PANELS: u32 = 32;

fn interpolation_fraction(value: f64, start: f64, end: f64) -> f64 {
    let width = end - start;
    if width.is_finite() {
        (value - start) / width
    } else {
        // Opposite-sign finite endpoints can have an infinite direct
        // difference. Scaling first preserves their finite affine ratio.
        let scale = start.abs().max(end.abs()).max(value.abs());
        ((value / scale) - (start / scale)) / ((end / scale) - (start / scale))
    }
}

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

impl ChebProfile {
    /// Validate the function object's finite domain and coefficients.
    ///
    /// `fs-cheb` enforces these properties in its constructors today, but
    /// `ChebProfile` is a public authority-boundary type. Keeping the checks
    /// here makes downstream `Result` APIs fail closed even if another
    /// constructor is added later.
    pub(crate) fn check_with_checkpoint<E>(
        &self,
        context: &str,
        out: &mut Vec<Violation>,
        checkpoint: &mut impl FnMut(&'static str) -> Result<(), E>,
    ) -> Result<(), E> {
        let (a, b) = self.cheb.domain();
        if !(a.is_finite() && b.is_finite() && a < b) {
            out.push(Violation {
                code: "signal-chebfun-domain",
                what: format!("{context}: chebfun domain [{a}, {b}] is invalid"),
                fix: "build the function object on a finite, nonempty interval".to_string(),
            });
        }
        let mut coefficients_invalid = self.cheb.coeffs().is_empty();
        for coefficient in self.cheb.coeffs() {
            checkpoint("signal Chebyshev coefficients")?;
            coefficients_invalid |= !coefficient.is_finite();
        }
        if coefficients_invalid {
            out.push(Violation {
                code: "signal-chebfun-coefficients",
                what: format!("{context}: chebfun coefficients are empty or non-finite"),
                fix: "supply at least one finite Chebyshev coefficient".to_string(),
            });
        }
        Ok(())
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

    /// Dynamically sized scalar payload visited by semantic validation.
    /// `None` denotes impossible `usize` addition rather than wrapping the
    /// admission plan.
    pub(crate) fn validation_dynamic_scalars(&self) -> Option<usize> {
        match self {
            Self::Constant(_) | Self::Ramp { .. } => Some(0),
            Self::Table { times, values, .. } => times.len().checked_add(values.len()),
            Self::Chebfun(profile) => Some(profile.cheb.coeffs().len()),
        }
    }

    /// Raw checkpoint contribution before set-local sort/deduplication.
    pub(crate) fn net_flux_validation_time_count(&self) -> usize {
        match self {
            Self::Constant(_) => 0,
            Self::Ramp { .. } => 2,
            Self::Table { times, .. } => times.len(),
            Self::Chebfun(_) => usize::try_from(SMOOTH_NET_FLUX_VALIDATION_PANELS + 1)
                .expect("the fixed smooth checkpoint count fits usize"),
        }
    }

    /// Append the deterministic time checkpoints needed by net-flux
    /// compatibility validation.
    ///
    /// Ramp and table histories contribute their exact breakpoints. A smooth
    /// Chebfun history contributes a bounded uniform grid over its actual
    /// declared domain because it has no finite breakpoint set. Constants add
    /// no checkpoint; callers always add `t = 0` for the shared baseline.
    pub(crate) fn append_net_flux_validation_times(&self, times: &mut Vec<f64>) {
        match self {
            TimeSignal::Constant(_) => {}
            TimeSignal::Ramp { t_start, t_end, .. } => {
                times.push(*t_start);
                times.push(*t_end);
            }
            TimeSignal::Table { times: samples, .. } => times.extend(samples.iter().copied()),
            TimeSignal::Chebfun(profile) => {
                let (start, end) = profile.cheb.domain();
                for panel in 0..=SMOOTH_NET_FLUX_VALIDATION_PANELS {
                    let time = if panel == 0 {
                        start
                    } else if panel == SMOOTH_NET_FLUX_VALIDATION_PANELS {
                        end
                    } else {
                        let alpha = f64::from(panel) / f64::from(SMOOTH_NET_FLUX_VALIDATION_PANELS);
                        (1.0 - alpha) * start + alpha * end
                    };
                    times.push(time);
                }
            }
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
        let mut violations = Vec::new();
        self.check("signal evaluation", &mut violations);
        if let Some(first) = violations.first() {
            return Err(ScenarioError::Evaluate {
                what: first.what.clone(),
            });
        }
        let value = match self {
            TimeSignal::Constant(q) => *q,
            TimeSignal::Ramp {
                t_start,
                t_end,
                from,
                to,
            } => {
                if t <= *t_start {
                    *from
                } else if t >= *t_end {
                    *to
                } else {
                    let alpha = interpolation_fraction(t, *t_start, *t_end);
                    QtyAny::new((1.0 - alpha) * from.value + alpha * to.value, from.dims)
                }
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
                            let alpha = interpolation_fraction(t, times[i - 1], times[i]);
                            (1.0 - alpha) * values[i - 1] + alpha * values[i]
                        }
                    },
                };
                QtyAny::new(v, *dims)
            }
            TimeSignal::Chebfun(profile) => {
                let (a, b) = profile.cheb.domain();
                QtyAny::new(profile.cheb.eval(t.clamp(a, b)), profile.dims)
            }
        };
        if !value.value.is_finite() {
            return Err(ScenarioError::Evaluate {
                what: "signal evaluation produced a non-finite value".to_string(),
            });
        }
        Ok(value)
    }

    /// Structural validation, accumulated as [`Violation`]s.
    pub fn check(&self, context: &str, out: &mut Vec<Violation>) {
        let mut checkpoint = |_: &'static str| Ok::<(), core::convert::Infallible>(());
        match self.check_with_checkpoint(context, out, &mut checkpoint) {
            Ok(()) => {}
            Err(never) => match never {},
        }
    }

    pub(crate) fn check_with_checkpoint<E>(
        &self,
        context: &str,
        out: &mut Vec<Violation>,
        checkpoint: &mut impl FnMut(&'static str) -> Result<(), E>,
    ) -> Result<(), E> {
        match self {
            TimeSignal::Constant(q) => {
                if !q.value.is_finite() {
                    out.push(Violation {
                        code: "signal-value-nonfinite",
                        what: format!("{context}: constant value {} is non-finite", q.value),
                        fix: "replace the constant with a finite value".to_string(),
                    });
                }
            }
            TimeSignal::Ramp {
                t_start,
                t_end,
                from,
                to,
            } => {
                if !t_start.is_finite() || !t_end.is_finite() || t_start >= t_end {
                    out.push(Violation {
                        code: "signal-ramp-times",
                        what: format!("{context}: ramp interval [{t_start}, {t_end}] is invalid"),
                        fix: "order the ramp times so t_start < t_end and both are finite"
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
                if !from.value.is_finite() || !to.value.is_finite() {
                    out.push(Violation {
                        code: "signal-value-nonfinite",
                        what: format!(
                            "{context}: ramp endpoint values [{}, {}] must be finite",
                            from.value, to.value
                        ),
                        fix: "replace both ramp endpoints with finite values".to_string(),
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
                let mut time_nonfinite = false;
                let mut order_invalid = false;
                let mut previous: Option<f64> = None;
                for &time in times {
                    checkpoint("signal table times")?;
                    time_nonfinite |= !time.is_finite();
                    if let Some(prior) = previous {
                        order_invalid |= !(prior.is_finite() && time > prior);
                    }
                    previous = Some(time);
                }
                let mut value_nonfinite = false;
                for value in values {
                    checkpoint("signal table values")?;
                    value_nonfinite |= !value.is_finite();
                }
                if time_nonfinite {
                    out.push(Violation {
                        code: "signal-table-time-nonfinite",
                        what: format!("{context}: table contains a non-finite sample time"),
                        fix: "replace every sample time with a finite value".to_string(),
                    });
                }
                if value_nonfinite {
                    out.push(Violation {
                        code: "signal-value-nonfinite",
                        what: format!("{context}: table contains a non-finite sample value"),
                        fix: "replace every sample value with a finite value".to_string(),
                    });
                }
                if order_invalid {
                    out.push(Violation {
                        code: "signal-table-order",
                        what: format!("{context}: table times are not strictly increasing"),
                        fix: "sort the samples by time and remove duplicates".to_string(),
                    });
                }
            }
            TimeSignal::Chebfun(profile) => {
                profile.check_with_checkpoint(context, out, checkpoint)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod validation_internal_tests {
    use super::{Interp, TimeSignal};
    use fs_qty::Dims;

    #[test]
    fn table_scan_is_single_pass_checkpointed_and_order_stable() {
        let signal = TimeSignal::Table {
            times: vec![0.0, f64::NAN, -1.0],
            values: vec![1.0, f64::INFINITY, 3.0],
            dims: Dims::NONE,
            interp: Interp::Linear,
        };
        let mut public_findings = Vec::new();
        signal.check("table", &mut public_findings);

        let mut checkpointed_findings = Vec::new();
        let mut phases = Vec::new();
        signal
            .check_with_checkpoint("table", &mut checkpointed_findings, &mut |phase| {
                phases.push(phase);
                Ok::<(), core::convert::Infallible>(())
            })
            .expect("infallible checkpoint");

        assert_eq!(checkpointed_findings, public_findings);
        assert_eq!(
            public_findings
                .iter()
                .map(|finding| finding.code)
                .collect::<Vec<_>>(),
            [
                "signal-table-time-nonfinite",
                "signal-value-nonfinite",
                "signal-table-order",
            ]
        );
        assert_eq!(
            phases,
            [
                "signal table times",
                "signal table times",
                "signal table times",
                "signal table values",
                "signal table values",
                "signal table values",
            ]
        );

        let mut cancelled_findings = Vec::new();
        let mut time_polls = 0usize;
        let result = signal.check_with_checkpoint("table", &mut cancelled_findings, &mut |phase| {
            if phase == "signal table times" {
                time_polls += 1;
            }
            if time_polls == 2 {
                Err("cancelled")
            } else {
                Ok(())
            }
        });
        assert_eq!(result, Err("cancelled"));
        assert!(cancelled_findings.is_empty());
    }
}
