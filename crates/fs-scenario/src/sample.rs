//! Shared time sampling for the original immutable scalar and vector carriers.

use crate::ScenarioError;
use crate::payload::{OutsideDomainPolicy, SampleSource, TableInterpolation};
use fs_qty::QtyAny;

/// No sample is returned after a policy, arithmetic, frame or caller refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum PayloadEvaluationError {
    Unsupported(&'static str),
    InvalidTime,
    OutsideDomain {
        time_s: f64,
        first_s: f64,
        last_s: f64,
    },
    Arithmetic(&'static str),
    Frame(ScenarioError),
    Cancelled,
}
impl std::fmt::Display for PayloadEvaluationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "scenario payload evaluation refused: {self:?}")
    }
}
impl std::error::Error for PayloadEvaluationError {}

pub(crate) fn value<T, const N: usize>(
    source: &SampleSource<T>,
    time_s: f64,
    components: impl Fn(&T) -> [f64; N],
) -> Result<[f64; N], PayloadEvaluationError> {
    match source {
        SampleSource::Fixed(value) => Ok(components(value)),
        SampleSource::Distribution { .. } => Err(PayloadEvaluationError::Unsupported(
            "a distribution is not an explicit time history",
        )),
        SampleSource::Table {
            times,
            values,
            interpolation,
            outside,
        } => table_value(times, values, *interpolation, *outside, time_s, components),
    }
}

fn table_value<T, const N: usize>(
    times: &[QtyAny],
    values: &[T],
    interpolation: TableInterpolation,
    outside: OutsideDomainPolicy,
    time_s: f64,
    components: impl Fn(&T) -> [f64; N],
) -> Result<[f64; N], PayloadEvaluationError> {
    use PayloadEvaluationError as E;
    // Nonempty, increasing, finite time and equal shapes are constructor invariants.
    let first = times[0].value;
    let last = times[times.len() - 1].value;
    let mut time = time_s;
    if time < first || time > last {
        time = match outside {
            OutsideDomainPolicy::Refuse => {
                return Err(E::OutsideDomain {
                    time_s,
                    first_s: first,
                    last_s: last,
                });
            }
            OutsideDomainPolicy::Clamp => time.clamp(first, last),
            OutsideDomainPolicy::Periodic => {
                let period = last - first;
                let delta = time - first;
                if !period.is_finite() || period <= 0.0 || !delta.is_finite() {
                    return Err(E::Arithmetic("unrepresentable periodic time reduction"));
                }
                let phase = delta.rem_euclid(period);
                let wrapped = first + phase;
                if !wrapped.is_finite() {
                    return Err(E::Arithmetic("nonfinite wrapped time"));
                }
                wrapped.clamp(first, last)
            }
        };
    }
    let left = times.partition_point(|t| t.value <= time).saturating_sub(1);
    if left + 1 == times.len()
        || time == times[left].value
        || interpolation == TableInterpolation::StepLeft
    {
        return Ok(components(&values[left]));
    }
    let width = times[left + 1].value - times[left].value;
    if !width.is_finite() || width <= 0.0 {
        return Err(E::Arithmetic("unrepresentable interpolation span"));
    }
    let weight = (time - times[left].value) / width;
    if !weight.is_finite() || !(0.0..=1.0).contains(&weight) {
        return Err(E::Arithmetic("unrepresentable interpolation weight"));
    }
    let a = components(&values[left]);
    let b = components(&values[left + 1]);
    // Preserve the original vector evaluator's convex arithmetic and exact knots.
    let out = std::array::from_fn(|i| (1.0 - weight) * a[i] + weight * b[i]);
    if out.iter().all(|v| v.is_finite()) {
        Ok(out)
    } else {
        Err(E::Arithmetic("nonfinite interpolated component"))
    }
}
