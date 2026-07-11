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
    assert!(
        half_width.is_finite() && half_width >= 0.0,
        "target half-width must be finite and non-negative"
    );
    // Bounded [0,1] variables are sub-Gaussian with sigma = 1/2
    // (Hoeffding), so the Gaussian-mixture CS applies.
    let mut cs = GaussianMixtureCs::new(0.5, 1.0, alpha);
    let mut sum = 0.0f64;
    let mut n = 0u64;
    while n < max_n {
        let x = sample(n);
        assert!(
            (0.0..=1.0).contains(&x),
            "probability observations must lie in [0,1]; got {x}"
        );
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
///
/// # Panics
/// If `samples` is empty, any sample is non-finite, or `beta` is not finite
/// and strictly between 0 and 1.
#[must_use]
pub fn cvar(samples: &[f64], beta: f64) -> f64 {
    assert!(!samples.is_empty(), "cvar needs at least one sample");
    assert!(
        beta.is_finite() && 0.0 < beta && beta < 1.0,
        "cvar beta must be finite and in (0, 1)"
    );
    assert!(
        samples.iter().all(|s| s.is_finite()),
        "cvar samples must be finite"
    );
    let mut sorted = samples.to_vec();
    sorted.sort_by(f64::total_cmp);
    // Standard finite-sample empirical CVaR (Acerbi–Tasche): the mean over the
    // worst n·(1−β) fraction of losses, applying a FRACTIONAL weight to the
    // boundary order statistic when n·β is not an integer. Averaging the
    // ⌈n(1−β)⌉ top samples at EQUAL weight (dividing by the rounded-up count)
    // systematically UNDER-reports the shortfall — the anti-conservative
    // direction for a risk functional feeding ASCENT robust formulations
    // (bead zsvk). For integer n·β the two coincide.
    #[allow(clippy::cast_precision_loss)]
    let n = samples.len() as f64;
    // 1-based rank of the boundary order statistic, m = ⌈nβ⌉ ∈ [1, n].
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    let m = (n * beta).ceil() as usize;
    let n_alpha = n * (1.0 - beta); // effective tail count (> 0 since β < 1)
    #[allow(clippy::cast_precision_loss)]
    let boundary_weight = (m as f64) - n * beta; // ⌈nβ⌉ − nβ ∈ [0, 1)
    let boundary = sorted[m - 1];
    let above: f64 = sorted[m..].iter().sum();
    (boundary_weight * boundary + above) / n_alpha
}
