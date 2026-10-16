//! PERF-REGRESSION CI (plan §14.4, bead fz2.4): performance
//! regressions are TEST FAILURES. This module is the statistics and
//! diagnosis layer over the roofline harness: DISPERSION-AWARE
//! tolerance bands (thermal jitter must not cry wolf), CUSUM
//! change-point alarms with calibrated thresholds, flame-graph-
//! equivalent attribution from the phase-annotated event stream (a
//! regression arrives WITH its own diagnosis), and the dashboard
//! one-liners ("what got slower this month, and why").

use std::collections::{BTreeMap, BTreeSet};

/// Maximum nights admitted to one kernel history.
pub const MAX_REGRESSION_HISTORY_NIGHTS: usize = 4_096;
/// Maximum kernels admitted to one dashboard query.
pub const MAX_REGRESSION_KERNELS: usize = 4_096;
/// Maximum aggregate nights admitted across one dashboard query.
///
/// This admits every kernel at the 14-night dashboard minimum while preventing
/// independent per-kernel history caps from multiplying into millions of rows.
pub const MAX_REGRESSION_DASHBOARD_NIGHTS: usize = 65_536;
/// Maximum aggregate phase observations admitted across one dashboard query.
pub const MAX_REGRESSION_DASHBOARD_PHASE_OBSERVATIONS: usize = 262_144;
/// Maximum phases admitted on one night.
pub const MAX_REGRESSION_PHASES_PER_NIGHT: usize = 1_024;
/// Maximum distinct phase identities admitted across one history.
pub const MAX_REGRESSION_UNIQUE_PHASES: usize = 4_096;
/// Maximum aggregate phase observations admitted across one history.
pub const MAX_REGRESSION_PHASE_OBSERVATIONS: usize = 262_144;
/// Maximum bytes in a kernel or phase identity.
pub const MAX_REGRESSION_NAME_BYTES: usize = 128;
/// Maximum samples admitted to one standardized residual stream.
pub const MAX_REGRESSION_SERIES_SAMPLES: usize = 1_000_000;

/// Invalid public regression-analysis input.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrendInputError {
    reason: String,
}

impl TrendInputError {
    /// Stable diagnostic for the invalid input.
    #[must_use]
    pub fn reason(&self) -> &str {
        &self.reason
    }
}

impl std::fmt::Display for TrendInputError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(&self.reason)
    }
}

impl std::error::Error for TrendInputError {}

/// One nightly observation of a kernel: attainment plus its
/// phase-annotated event stream (phase name → seconds).
#[derive(Debug, Clone, PartialEq)]
pub struct Night {
    /// Nightly index (logical time). Histories must be strictly increasing;
    /// gaps are permitted, but duplicates and reversals are invalid evidence.
    pub night: u64,
    /// Roofline attainment in [0, 1]-ish.
    pub attainment: f64,
    /// Phase durations (the flame-graph-equivalent source).
    pub phases: BTreeMap<String, f64>,
}

/// The dispersion-aware gate: a run fails when attainment drops more
/// than `k_sigma` baseline standard deviations below the baseline
/// mean — the statistical band, not a naive threshold.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GateSpec {
    /// Band width in baseline sigmas.
    pub k_sigma: f64,
    /// Minimum baseline nights before the gate arms.
    pub min_baseline: usize,
}

impl Default for GateSpec {
    fn default() -> Self {
        GateSpec {
            k_sigma: 4.0,
            min_baseline: 8,
        }
    }
}

fn normalized_mean_std(xs: &[f64], scale: f64) -> Option<(f64, f64)> {
    if xs.is_empty() {
        return Some((0.0, 0.0));
    }
    if !scale.is_finite() {
        return None;
    }
    if scale == 0.0 {
        return Some((0.0, 0.0));
    }
    #[allow(clippy::cast_precision_loss)]
    let n = xs.len() as f64;
    let mean = xs.iter().map(|value| value / scale).sum::<f64>() / n;
    let variance = xs
        .iter()
        .map(|value| {
            let delta = value / scale - mean;
            delta * delta
        })
        .sum::<f64>()
        / (n - 1.0).max(1.0);
    let sigma = variance.sqrt();
    (mean.is_finite() && sigma.is_finite()).then_some((mean, sigma))
}

fn standardized_score(baseline: &[f64], newest: f64) -> Option<f64> {
    let scale = baseline
        .iter()
        .map(|value| value.abs())
        .chain(std::iter::once(newest.abs()))
        .fold(0.0_f64, f64::max);
    let (mean, sigma) = normalized_mean_std(baseline, scale)?;
    let newest = if scale == 0.0 { 0.0 } else { newest / scale };
    let residual = newest - mean;
    if residual == 0.0 {
        return Some(0.0);
    }
    if sigma == 0.0 {
        return Some(f64::INFINITY.copysign(residual));
    }
    let z = residual / sigma;
    (!z.is_nan()).then_some(z)
}

fn normalized_window_means(left: &[f64], right: &[f64]) -> Option<(f64, f64)> {
    let scale = left
        .iter()
        .chain(right)
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    let (left_mean, _) = normalized_mean_std(left, scale)?;
    let (right_mean, _) = normalized_mean_std(right, scale)?;
    Some((left_mean, right_mean))
}

/// The gate verdict for the newest night against its baseline.
#[derive(Debug, Clone, PartialEq)]
pub enum GateVerdict {
    /// Within the band (or the gate is not yet armed).
    Green {
        /// Standardized score of the newest night.
        z: f64,
    },
    /// RED: the regression, with its own diagnosis attached.
    Red {
        /// Standardized drop (negative).
        z: f64,
        /// The flame-graph-level attribution: phases ranked by their
        /// share growth vs the last-green baseline (top offender
        /// first), as (phase, baseline share, regressed share).
        attribution: Vec<(String, f64, f64)>,
    },
    /// MALFORMED EVIDENCE (bead fz2.4.1): non-finite or negative
    /// inputs, or an unusable spec. A proof-bearing gate never
    /// represents bad data as Green — it says so, with a diagnosis.
    Invalid {
        /// What was malformed (structured, human-readable).
        reason: String,
    },
}

/// First flaw in a night's fields, if any (the fail-closed screen).
fn night_flaw(idx: usize, night: &Night) -> Option<String> {
    if !night.attainment.is_finite() || night.attainment < 0.0 {
        return Some(format!(
            "night {idx} (index in history): attainment {} is not finite and non-negative",
            night.attainment
        ));
    }
    if night.phases.len() > MAX_REGRESSION_PHASES_PER_NIGHT {
        return Some(format!(
            "night {idx}: phase count {} exceeds limit {MAX_REGRESSION_PHASES_PER_NIGHT}",
            night.phases.len()
        ));
    }
    for (phase, &secs) in &night.phases {
        if phase.is_empty()
            || phase.len() > MAX_REGRESSION_NAME_BYTES
            || !phase.bytes().all(|byte| byte.is_ascii_graphic())
        {
            return Some(format!(
                "night {idx}: phase identity must contain 1..={MAX_REGRESSION_NAME_BYTES} ASCII graphic bytes"
            ));
        }
        if !secs.is_finite() || secs < 0.0 {
            return Some(format!(
                "night {idx}: phase '{phase}' duration {secs} is not finite and non-negative"
            ));
        }
    }
    None
}

/// First malformed field or non-monotone logical-time transition in a history.
fn history_flaw(history: &[Night]) -> Option<String> {
    if history.len() > MAX_REGRESSION_HISTORY_NIGHTS {
        return Some(format!(
            "history length {} exceeds limit {MAX_REGRESSION_HISTORY_NIGHTS}",
            history.len()
        ));
    }
    let mut phase_observations = 0usize;
    let mut unique_phases = BTreeSet::new();
    for (idx, night) in history.iter().enumerate() {
        if let Some(reason) = night_flaw(idx, night) {
            return Some(reason);
        }
        phase_observations = match phase_observations.checked_add(night.phases.len()) {
            Some(observations) if observations <= MAX_REGRESSION_PHASE_OBSERVATIONS => observations,
            Some(observations) => {
                return Some(format!(
                    "phase observation count {observations} exceeds limit {MAX_REGRESSION_PHASE_OBSERVATIONS}"
                ));
            }
            None => return Some("phase observation count overflowed usize".to_string()),
        };
        for phase in night.phases.keys() {
            unique_phases.insert(phase.as_str());
            if unique_phases.len() > MAX_REGRESSION_UNIQUE_PHASES {
                return Some(format!(
                    "unique phase count {} exceeds limit {MAX_REGRESSION_UNIQUE_PHASES}",
                    unique_phases.len()
                ));
            }
        }
    }
    history.windows(2).enumerate().find_map(|(idx, pair)| {
        let previous = pair[0].night;
        let next = pair[1].night;
        (next <= previous).then(|| {
            format!(
                "logical night must increase strictly: history indices {idx} and {} contain {previous} then {next}",
                idx + 1
            )
        })
    })
}

/// First flaw in a spec, if any.
fn spec_flaw(spec: GateSpec) -> Option<String> {
    if !(spec.k_sigma.is_finite() && spec.k_sigma > 0.0) {
        return Some(format!(
            "spec.k_sigma {} is not finite and positive",
            spec.k_sigma
        ));
    }
    if spec.min_baseline < 2 {
        return Some(format!(
            "spec.min_baseline {} cannot support a dispersion estimate (need >= 2)",
            spec.min_baseline
        ));
    }
    if spec.min_baseline >= MAX_REGRESSION_HISTORY_NIGHTS {
        return Some(format!(
            "spec.min_baseline {} cannot fit a baseline plus newest night inside the {MAX_REGRESSION_HISTORY_NIGHTS}-night history limit",
            spec.min_baseline
        ));
    }
    None
}

/// Median phase shares for one comparison window. Nights where a phase is
/// absent contribute an implicit zero without materializing the full product.
fn phase_share_medians(window: &[Night]) -> BTreeMap<String, f64> {
    let mut shares_by_phase: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for night in window {
        let scale = night.phases.values().copied().fold(0.0_f64, f64::max);
        let scaled_total = if scale == 0.0 {
            0.0
        } else {
            night.phases.values().map(|secs| secs / scale).sum()
        };
        for (phase, secs) in &night.phases {
            let share = if scaled_total == 0.0 {
                0.0
            } else {
                (secs / scale) / scaled_total
            };
            shares_by_phase
                .entry(phase.clone())
                .or_default()
                .push(share);
        }
    }
    let median_index = window.len() / 2;
    shares_by_phase
        .into_iter()
        .map(|(phase, mut shares)| {
            let missing_zeros = window.len().saturating_sub(shares.len());
            shares.sort_by(f64::total_cmp);
            let median = if median_index < missing_zeros {
                0.0
            } else {
                shares
                    .get(median_index - missing_zeros)
                    .copied()
                    .unwrap_or(0.0)
            };
            (phase, median)
        })
        .collect()
}

fn attribution_between_windows(
    baseline: &[Night],
    comparison: &[Night],
) -> Vec<(String, f64, f64)> {
    let base_shares = phase_share_medians(baseline);
    let comparison_shares = phase_share_medians(comparison);
    let mut phases = BTreeSet::new();
    phases.extend(base_shares.keys().cloned());
    phases.extend(comparison_shares.keys().cloned());
    let mut attribution: Vec<(String, f64, f64)> = phases
        .into_iter()
        .map(|phase| {
            let base = base_shares.get(&phase).copied().unwrap_or(0.0);
            let comparison = comparison_shares.get(&phase).copied().unwrap_or(0.0);
            (phase, base, comparison)
        })
        .collect();
    attribution.sort_by(|a, b| (b.2 - b.1).total_cmp(&(a.2 - a.1)).then(a.0.cmp(&b.0)));
    attribution
}

/// Gate the newest night against the preceding baseline, attributing
/// any red to the phases whose SHARE of the total grew most — the
/// event stream reconstructing the flame-graph diff post hoc.
///
/// FAIL-CLOSED (bead fz2.4.1): non-finite or negative attainment or
/// phase durations or non-increasing logical time anywhere in the history, and
/// unusable specs (non-finite/non-positive k_sigma, or a baseline outside the
/// admitted range), return [`GateVerdict::Invalid`] — never Green. NaN can
/// otherwise flip the red predicate false silently.
#[must_use]
pub fn gate(history: &[Night], spec: GateSpec) -> GateVerdict {
    if let Some(reason) = spec_flaw(spec) {
        return GateVerdict::Invalid { reason };
    }
    if let Some(reason) = history_flaw(history) {
        return GateVerdict::Invalid { reason };
    }
    let n = history.len();
    let Some(required_nights) = spec.min_baseline.checked_add(1) else {
        return GateVerdict::Invalid {
            reason: "spec.min_baseline overflowed the required history length".to_string(),
        };
    };
    if n < required_nights {
        return GateVerdict::Green { z: 0.0 };
    }
    let (baseline, newest) = history.split_at(n - 1);
    let newest = &newest[0];
    let xs: Vec<f64> = baseline.iter().map(|b| b.attainment).collect();
    let Some(z) = standardized_score(&xs, newest.attainment) else {
        return GateVerdict::Invalid {
            reason: "baseline mean/dispersion is not representable".to_string(),
        };
    };
    if z >= -spec.k_sigma {
        return GateVerdict::Green { z };
    }
    GateVerdict::Red {
        z,
        attribution: attribution_between_windows(baseline, std::slice::from_ref(newest)),
    }
}

/// A one-sided CUSUM change-point detector for slow drifts the
/// per-night gate misses: alarm when the cumulative standardized
/// shortfall crosses `h`. `k` is the slack (drift smaller than k·σ per
/// night is absorbed — the noise-robustness knob).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cusum {
    /// Slack per observation (in sigmas).
    pub k: f64,
    /// Alarm threshold (in cumulative sigmas).
    pub h: f64,
}

impl Default for Cusum {
    fn default() -> Self {
        Cusum { k: 0.5, h: 8.0 }
    }
}

impl Cusum {
    /// Run over standardized residuals (baseline-calibrated z-scores);
    /// returns the first alarm index, if any.
    ///
    /// FAIL-CLOSED (bead fz2.4.1): a detector with a non-finite or
    /// non-positive threshold cannot certify quiet — it alarms at
    /// index 0; a non-finite residual alarms at ITS index (NaN would
    /// otherwise silently reset the shortfall via `max`, suppressing
    /// detection). Malformed data can force an alarm; it can never
    /// suppress one.
    #[must_use]
    pub fn first_alarm(&self, z_scores: &[f64]) -> Option<usize> {
        if !(self.k.is_finite() && self.k >= 0.0 && self.h.is_finite() && self.h > 0.0) {
            return if z_scores.is_empty() { None } else { Some(0) };
        }
        if z_scores.len() > MAX_REGRESSION_SERIES_SAMPLES {
            return Some(0);
        }
        let mut s = 0.0f64;
        for (i, &z) in z_scores.iter().enumerate() {
            if !z.is_finite() {
                return Some(i);
            }
            s = (s - z - self.k).max(0.0); // accumulate SHORTFALL
            if s > self.h {
                return Some(i);
            }
        }
        None
    }
}

/// Standardize a history against its own expanding baseline (each
/// night scored against the nights before it; the first `warmup`
/// nights score 0).
///
/// FAIL-CLOSED (bead fz2.4.1): from the first non-finite input onward
/// every output is −∞ (the worst possible shortfall), so poisoned
/// history can never enter the expanding baseline as ordinary data or
/// read as good performance — downstream CUSUM alarms instead.
/// Against exact zero dispersion, equality scores zero and a change saturates
/// at the sign-appropriate finite f64 bound; scale does not alter the verdict.
///
/// # Errors
/// Returns [`TrendInputError`] when the residual stream exceeds its deterministic
/// sample budget.
pub fn standardize(history: &[f64], warmup: usize) -> Result<Vec<f64>, TrendInputError> {
    if history.len() > MAX_REGRESSION_SERIES_SAMPLES {
        return Err(TrendInputError {
            reason: format!(
                "standardization sample count {} exceeds limit {MAX_REGRESSION_SERIES_SAMPLES}",
                history.len()
            ),
        });
    }
    let mut out = Vec::with_capacity(history.len());
    let mut count = 0usize;
    let mut scale = 0.0_f64;
    let mut mean_scaled = 0.0_f64;
    let mut m2_scaled = 0.0_f64;
    let mut poisoned = false;
    for (index, &value) in history.iter().enumerate() {
        if poisoned || !value.is_finite() {
            poisoned = true;
            out.push(f64::NEG_INFINITY);
            continue;
        }
        let z = if index < warmup || count == 0 {
            0.0
        } else {
            let score_scale = scale.max(value.abs());
            let scale_ratio = if score_scale == 0.0 {
                0.0
            } else {
                scale / score_scale
            };
            let mean = mean_scaled * scale_ratio;
            let sigma = if count < 2 {
                0.0
            } else {
                #[allow(clippy::cast_precision_loss)]
                let denominator = (count - 1) as f64;
                (m2_scaled / denominator).sqrt() * scale_ratio
            };
            let value_scaled = if score_scale == 0.0 {
                0.0
            } else {
                value / score_scale
            };
            let residual = value_scaled - mean;
            if residual == 0.0 {
                0.0
            } else if sigma == 0.0 {
                f64::MAX.copysign(residual)
            } else {
                let standardized = residual / sigma;
                if standardized.is_infinite() {
                    f64::MAX.copysign(standardized)
                } else {
                    standardized
                }
            }
        };
        if z.is_finite() {
            out.push(z);
        } else {
            poisoned = true;
            out.push(f64::NEG_INFINITY);
            continue;
        }

        let new_scale = scale.max(value.abs());
        if new_scale > scale {
            let ratio = if new_scale == 0.0 {
                0.0
            } else {
                scale / new_scale
            };
            mean_scaled *= ratio;
            m2_scaled *= ratio * ratio;
            scale = new_scale;
        }
        let value_scaled = if scale == 0.0 { 0.0 } else { value / scale };
        count += 1;
        #[allow(clippy::cast_precision_loss)]
        let count_f64 = count as f64;
        let delta = value_scaled - mean_scaled;
        mean_scaled += delta / count_f64;
        let delta_after = value_scaled - mean_scaled;
        m2_scaled += delta * delta_after;
    }
    Ok(out)
}

/// THE DASHBOARD ONE-LINER: "what got slower this month, and why" —
/// kernels whose trailing-week mean attainment dropped more than
/// `pct_floor` percent below their opening-week mean, each with its
/// top-offender phase from the same opening-week/trailing-week comparison.
///
/// FAIL-CLOSED (bead fz2.4.1): a kernel whose history contains
/// non-finite or negative fields or non-increasing logical time is reported
/// FIRST with an infinite drop and the flaw as its "why" — malformed evidence
/// is flagged loudest, never silently skipped and never allowed to poison the
/// trend arithmetic of valid kernels.
///
/// # Errors
/// Returns [`TrendInputError`] when `pct_floor` is non-finite or negative, a
/// kernel identity is malformed, or the dashboard exceeds an aggregate input
/// budget. A malformed individual history remains an explicit `INVALID` row in
/// the successful report so one bad kernel cannot hide valid regressions in its
/// neighbors.
pub fn slower_this_month(
    kernels: &BTreeMap<String, Vec<Night>>,
    pct_floor: f64,
) -> Result<Vec<(String, f64, String)>, TrendInputError> {
    if !pct_floor.is_finite() || pct_floor < 0.0 {
        return Err(TrendInputError {
            reason: format!("pct_floor {pct_floor} is not finite and non-negative"),
        });
    }
    if kernels.len() > MAX_REGRESSION_KERNELS {
        return Err(TrendInputError {
            reason: format!(
                "kernel count {} exceeds limit {MAX_REGRESSION_KERNELS}",
                kernels.len()
            ),
        });
    }
    for kernel in kernels.keys() {
        if kernel.is_empty()
            || kernel.len() > MAX_REGRESSION_NAME_BYTES
            || !kernel.bytes().all(|byte| byte.is_ascii_graphic())
        {
            return Err(TrendInputError {
                reason: format!(
                    "kernel identity must contain 1..={MAX_REGRESSION_NAME_BYTES} ASCII graphic bytes"
                ),
            });
        }
    }
    let mut dashboard_nights = 0usize;
    let mut dashboard_phase_observations = 0usize;
    for history in kernels.values() {
        dashboard_nights =
            dashboard_nights
                .checked_add(history.len())
                .ok_or_else(|| TrendInputError {
                    reason: "dashboard night count overflowed usize".to_string(),
                })?;
        if dashboard_nights > MAX_REGRESSION_DASHBOARD_NIGHTS {
            return Err(TrendInputError {
                reason: format!(
                    "dashboard night count {dashboard_nights} exceeds limit {MAX_REGRESSION_DASHBOARD_NIGHTS}"
                ),
            });
        }
        for night in history {
            dashboard_phase_observations = dashboard_phase_observations
                .checked_add(night.phases.len())
                .ok_or_else(|| TrendInputError {
                    reason: "dashboard phase observation count overflowed usize".to_string(),
                })?;
            if dashboard_phase_observations > MAX_REGRESSION_DASHBOARD_PHASE_OBSERVATIONS {
                return Err(TrendInputError {
                    reason: format!(
                        "dashboard phase observation count {dashboard_phase_observations} exceeds limit {MAX_REGRESSION_DASHBOARD_PHASE_OBSERVATIONS}"
                    ),
                });
            }
        }
    }
    let mut out = Vec::new();
    for (kernel, history) in kernels {
        if let Some(flaw) = history_flaw(history) {
            out.push((kernel.clone(), f64::INFINITY, format!("INVALID: {flaw}")));
            continue;
        }
        if history.len() < 14 {
            continue;
        }
        let head: Vec<f64> = history[..7].iter().map(|n| n.attainment).collect();
        let tail: Vec<f64> = history[history.len() - 7..]
            .iter()
            .map(|n| n.attainment)
            .collect();
        let Some((mu_head, mu_tail)) = normalized_window_means(&head, &tail) else {
            out.push((
                kernel.clone(),
                f64::INFINITY,
                "INVALID: trend mean/dispersion is not representable".to_string(),
            ));
            continue;
        };
        if mu_tail >= mu_head {
            continue;
        }
        debug_assert!(
            mu_head > 0.0,
            "a non-negative mean cannot decline below zero"
        );
        let drop_pct = (mu_head - mu_tail) / mu_head * 100.0;
        if drop_pct > pct_floor {
            let why = attribution_between_windows(&history[..7], &history[history.len() - 7..])
                .first()
                .map_or_else(|| "unattributed".to_string(), |(p, _, _)| p.clone());
            out.push((kernel.clone(), drop_pct, why));
        }
    }
    out.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
    Ok(out)
}
