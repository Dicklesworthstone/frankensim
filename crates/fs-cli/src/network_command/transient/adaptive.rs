//! Step doubling on the actual coupled backward-Euler map. The accepted path
//! consists of two half-steps; the coarse trial is never committed or used for
//! heat accounting. The discrepancy estimates local endpoint error, not a
//! global error enclosure, midpoint error bound, or continuous-time peak.

use super::*;
use fs_math::det;

type Endpoint = (CoupledTransportSolution, StepSolution);

#[derive(Debug, Clone, Copy)]
pub(super) struct Config {
    pub absolute_tolerance_k: f64,
    pub relative_tolerance: f64,
    pub minimum_trial_step_s: f64,
    pub max_trials: usize,
}

impl Config {
    pub(super) fn parse(value: &J, max_step_s: f64) -> Result<Self> {
        object(value, &["absolute_tolerance_k", "relative_tolerance", "minimum_trial_step_s", "max_trials"], "transient.adaptive")?;
        let config = Self {
            absolute_tolerance_k: positive(get(value, "absolute_tolerance_k")?, "adaptive.absolute_tolerance_k")?,
            relative_tolerance: number(get(value, "relative_tolerance")?, "adaptive.relative_tolerance")?,
            minimum_trial_step_s: positive(get(value, "minimum_trial_step_s")?, "minimum_trial_step_s")?,
            max_trials: count(get(value, "max_trials")?, "adaptive.max_trials", 100_000)?,
        };
        if !(0.0..1.0).contains(&config.relative_tolerance) || config.minimum_trial_step_s > max_step_s {
            return Err(bad("adaptive relative_tolerance must be in [0,1), and minimum_trial_step_s <= max_step_s"));
        }
        Ok(config)
    }
}

#[derive(Debug, Default)]
pub(super) struct Stats {
    pub trials: usize,
    pub rejected: usize,
    pub largest_accepted_ratio: f64,
}

impl Stats {
    pub(super) fn render(&self, config: Option<Config>) -> Result<String> {
        let Some(config) = config else { return Ok("null".into()); };
        Ok(format!("{{\"method\":\"backward-euler-step-doubling\",\"trials\":{},\"rejected_trials\":{},\"accepted_trials\":{},\"largest_accepted_error_ratio\":{},\"absolute_tolerance_k\":{},\"relative_temperature_change_tolerance\":{},\"minimum_trial_step_s\":{},\"max_trials\":{},\"scope\":\"full-field local endpoint discrepancy estimate; accepts two half-steps without extrapolation; no global, midpoint, inter-step-peak or continuum error bound\"}}",
            self.trials, self.rejected, self.trials - self.rejected, num(self.largest_accepted_ratio)?,
            num(config.absolute_tolerance_k)?, num(config.relative_tolerance)?, num(config.minimum_trial_step_s)?, config.max_trials))
    }
}

pub(super) struct Accepted {
    pub samples: [(f64, Endpoint); 2],
    pub error_ratio: f64,
    pub next_trial_s: f64,
}

/// Retry from the same old state until a complete pair is admitted. The caller
/// counts ALL solid work in `advance`, but records only returned samples.
/// Producer errors propagate unchanged; only an excessive error estimate retries.
pub(super) fn step(
    cx: &Cx<'_>, old: &[f64], time: f64, interval_end: f64, suggested_s: f64,
    max_step_s: f64, config: Config, stats: &mut Stats,
    mut advance: impl FnMut(&[f64], f64) -> Result<Endpoint>,
) -> Result<Accepted> {
    let mut proposed = suggested_s.min(max_step_s).min(interval_end - time);
    loop {
        poll(cx)?;
        if stats.trials >= config.max_trials { return Err(budget("adaptive trial budget exhausted; no partial trajectory published")); }
        let end = if proposed >= interval_end - time { interval_end } else { time + proposed };
        let width = end - time;
        let middle = time + 0.5 * width;
        if !(width.is_finite() && middle > time && middle < end && end <= interval_end) {
            return Err(resolution("time representation cannot split this adaptive trial"));
        }
        stats.trials += 1;
        let coarse = advance(old, width)?;
        poll(cx)?;
        let first = advance(old, middle - time)?;
        poll(cx)?;
        let second = advance(&first.1.temperature, end - middle)?;
        let ratio = error_ratio(cx, old, &coarse.1.temperature, &second.1.temperature, config)?;
        poll(cx)?;
        if ratio <= 1.0 {
            stats.largest_accepted_ratio = stats.largest_accepted_ratio.max(ratio);
            let factor = factor(ratio, true);
            // Avoid overflowing a capped suggestion for very long intervals.
            let next_trial_s = if width >= max_step_s / factor { max_step_s } else { width * factor };
            return Ok(Accepted { samples: [(middle, first), (end, second)], error_ratio: ratio,
                next_trial_s: next_trial_s.max(config.minimum_trial_step_s).min(max_step_s) });
        }
        stats.rejected += 1;
        let reduced = (width * factor(ratio, false)).max(config.minimum_trial_step_s);
        if !(reduced < width && time + reduced < end) {
            return Err(resolution("local temperature discrepancy still exceeds tolerance at the minimum representable/admitted trial step"));
        }
        // No rejected midpoint or field becomes history; retry without recursion.
        proposed = reduced;
    }
}

fn error_ratio(cx: &Cx<'_>, old: &[f64], coarse: &[f64], fine: &[f64], config: Config) -> Result<f64> {
    if old.is_empty() || old.len() != coarse.len() || old.len() != fine.len() {
        return Err(bad("adaptive temperature vectors must have the same nonzero length"));
    }
    let mut largest = 0.0_f64;
    for (index, ((&old, &coarse), &fine)) in old.iter().zip(coarse).zip(fine).enumerate() {
        if index % 512 == 0 { poll(cx)?; }
        finite(old)?; finite(coarse)?; finite(fine)?;
        // Relative tolerance scales the physical change, not the 300 K offset.
        let change = finite(coarse - old)?.abs().max(finite(fine - old)?.abs());
        let tolerance = finite(config.absolute_tolerance_k + config.relative_tolerance * change)?;
        let ratio = finite(finite(fine - coarse)?.abs() / tolerance)?;
        largest = largest.max(ratio);
    }
    Ok(largest)
}

fn factor(ratio: f64, accepted: bool) -> f64 {
    if ratio == 0.0 { return 2.0; }
    (0.9 / det::sqrt(ratio)).clamp(0.2, if accepted { 2.0 } else { 0.8 })
}
fn resolution(message: &str) -> Failure { Failure { code: "cooling-network-adaptive-resolution", message: message.into() } }

#[cfg(test)]
mod tests;
