//! Shared bounded caloric inverse; constitutive states remain with their owners.

use super::BinaryCaloricErrorV1;

pub(crate) struct TemperatureControl {
    pub bracket: [f64; 2],
    pub temperature_tolerance_k: f64,
    pub energy_tolerance_j: f64,
    pub maximum_iterations: usize,
}

pub(crate) struct TemperatureSample<S> {
    pub state: S,
    pub change_j: f64,
    pub residual_j: f64,
}

pub(crate) struct TemperatureRoot<S> {
    pub sample: TemperatureSample<S>,
    pub bracket: [f64; 2],
    pub iterations: usize,
}

/// The caller admits controls, source domains and heat before entering. Every
/// endpoint change comes from a sampled original residual. Optional source
/// tangents accelerate the solve; bisection continues advancing the far end.
pub(crate) fn solve<S, E: From<BinaryCaloricErrorV1>>(
    controls: TemperatureControl,
    mut sample: impl FnMut(f64) -> Result<TemperatureSample<S>, E>,
    mut heat_capacity: impl FnMut(&S) -> Option<f64>,
) -> Result<TemperatureRoot<S>, E> {
    let [mut lo, mut hi] = controls.bracket;
    let mut lower = sample(lo)?;
    let mut upper = sample(hi)?;
    if lower.residual_j > 0.0 || upper.residual_j < 0.0 {
        return Err(BinaryCaloricErrorV1::HeatOutsideBracket {
            residuals_j: [lower.residual_j, upper.residual_j],
        }
        .into());
    }
    let mut consecutive_newton_steps = 2;
    for iterations in 0..=controls.maximum_iterations {
        let lower_is_best = lower.residual_j.abs() <= upper.residual_j.abs();
        let best = if lower_is_best { &lower } else { &upper };
        if best.residual_j == 0.0
            || ((hi - lo) <= controls.temperature_tolerance_k
                && best.residual_j.abs() <= controls.energy_tolerance_j)
        {
            return Ok(TemperatureRoot {
                sample: if lower_is_best { lower } else { upper },
                bracket: [lo, hi],
                iterations,
            });
        }
        if iterations == controls.maximum_iterations {
            break;
        }
        let midpoint = 0.5 * lo + 0.5 * hi;
        if midpoint == lo || midpoint == hi {
            break;
        }
        let old_width = hi - lo;
        let old_residual = best.residual_j.abs();
        let mut next = midpoint;
        if consecutive_newton_steps < 2
            && let Some(capacity) = heat_capacity(&best.state)
        {
            let temperature = if lower_is_best { lo } else { hi };
            let correction = -best.residual_j / capacity;
            let mut proposal = temperature + correction;
            if correction.is_finite()
                && correction.abs() <= 0.25 * controls.temperature_tolerance_k
                && old_residual <= controls.energy_tolerance_j
            {
                // Establish a real narrow sign bracket instead of accepting a
                // Newton estimate whose opposite endpoint is still far away.
                let distance = 0.5 * controls.temperature_tolerance_k;
                proposal = if best.residual_j < 0.0 {
                    (temperature + distance).max(temperature.next_up())
                } else {
                    (temperature - distance).min(temperature.next_down())
                };
            }
            if proposal.is_finite() && proposal > lo && proposal < hi {
                next = proposal;
            }
        }
        let next_sample = sample(next)?;
        if next_sample.residual_j <= 0.0 {
            lo = next;
            lower = next_sample;
        } else {
            hi = next;
            upper = next_sample;
        }
        consecutive_newton_steps = if next == midpoint {
            0
        } else if hi - lo <= 0.5 * old_width
            || lower.residual_j.abs().min(upper.residual_j.abs()) <= 0.5 * old_residual
        {
            consecutive_newton_steps + 1
        } else {
            2
        };
    }
    Err(BinaryCaloricErrorV1::IterationLimit.into())
}
