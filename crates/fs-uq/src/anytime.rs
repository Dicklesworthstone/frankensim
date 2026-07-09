//! ANYTIME-VALID STOPPING (bead o5kc, Bet 5): every stochastic
//! estimate under an e-process confidence sequence — sample until the
//! CS is tight enough FOR THE DECISION AT HAND, then stop, validly,
//! automatically. Optional stopping is safe BY CONSTRUCTION (the CS
//! is valid at every stopping time), which is what lets the fragility
//! study stop itself the moment the estimate is decision-grade.

use fs_eproc::GaussianMixtureCs;

/// The stopped estimate: the point value, the anytime-valid interval
/// it stopped inside, and the samples it took to get there.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AnytimeEstimate {
    /// Sample mean at the stopping time.
    pub mean: f64,
    /// The confidence-sequence interval at the stop.
    pub lo: f64,
    /// Upper end.
    pub hi: f64,
    /// Samples consumed.
    pub n: u64,
    /// True iff the target half-width was reached within the cap.
    pub converged: bool,
}

/// Estimate a BOUNDED-[0,1] probability/mean with a sub-Gaussian
/// (σ = 1/2) mixture confidence sequence, stopping as soon as the CS
/// half-width is at most `half_width` (or at `max_n`). Valid at the
/// stopping time by construction — no peeking penalty.
pub fn estimate_probability_anytime(
    mut sample: impl FnMut(u64) -> f64,
    alpha: f64,
    half_width: f64,
    max_n: u64,
) -> AnytimeEstimate {
    // Bounded [0,1] variables are sub-Gaussian with sigma = 1/2
    // (Hoeffding), so the Gaussian-mixture CS applies.
    let mut cs = GaussianMixtureCs::new(0.5, 1.0, alpha);
    let mut sum = 0.0f64;
    let mut n = 0u64;
    while n < max_n {
        let x = sample(n);
        debug_assert!((0.0..=1.0).contains(&x), "bounded outcomes only");
        cs.observe(x);
        sum += x;
        n += 1;
        // fs-eproc's interval() returns (CENTER, RADIUS).
        if let Some((center, radius)) = cs.interval()
            && radius <= half_width
        {
            return AnytimeEstimate {
                mean: sum / n as f64,
                lo: (center - radius).max(0.0),
                hi: (center + radius).min(1.0),
                n,
                converged: true,
            };
        }
    }
    let (center, radius) = cs.interval().unwrap_or((0.5, 0.5));
    AnytimeEstimate {
        mean: if n > 0 { sum / n as f64 } else { 0.5 },
        lo: (center - radius).max(0.0),
        hi: (center + radius).min(1.0),
        n,
        converged: false,
    }
}

/// CVaR (expected shortfall) of loss samples at level `beta` — the
/// risk functional ASCENT's robust formulations consume (kept here as
/// the UQ-side entry point; fs-robust hosts the ASCENT-side twin).
#[must_use]
pub fn cvar(samples: &[f64], beta: f64) -> f64 {
    let mut sorted = samples.to_vec();
    sorted.sort_by(f64::total_cmp);
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    let cut = ((samples.len() as f64) * beta).floor() as usize;
    let tail = &sorted[cut.min(samples.len().saturating_sub(1))..];
    #[allow(clippy::cast_precision_loss)]
    {
        tail.iter().sum::<f64>() / tail.len().max(1) as f64
    }
}
