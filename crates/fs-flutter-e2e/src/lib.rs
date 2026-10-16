//! fs-flutter-e2e — FlutterCert: a PROVEN fluid-structure stability boundary.
//! Layer: L4 (ASCENT).
//!
//! # The campaign
//!
//! Flutter analysis traditionally means sweeping a parameter and plotting where
//! a damping curve crosses zero — a picture, not a proof. This locates the
//! added-mass instability boundary and hands back a machine-checked certificate,
//! composing crates never designed to meet:
//!
//! - **The model.** A minimal 2-DOF coupled operator `A(μ) = [[−1, 1], [μ−1,
//!   −1]]`, where `μ` is the added-mass ratio. `trace A = −2 < 0` always;
//!   `det A = 2 − μ`, so the system is asymptotically stable iff `μ < 2` — the
//!   flutter boundary is `μ* = 2`.
//! - **The proof** ([`fs_sos`]): `lyapunov_certifies_stability(A(μ), I)` checks
//!   `P ≻ 0` and `−(AᵀP + PA) ≻ 0`. With `P = I` this reduces to the eigenvalues
//!   of `[[2, −μ], [−μ, 2]]` being positive, i.e. `μ < 2` — the certificate
//!   recovers the EXACT boundary and is `Verified`.
//! - **The independent cross-check.** The Lyapunov `P=I` proof is only a
//!   SUFFICIENT condition (it equals `max eig((Aᵀ+A)/2) < 0`, i.e. `−1+μ/2 < 0`).
//!   The necessary-AND-sufficient criterion is `A`'s ACTUAL eigenvalues
//!   `−1 ± √(μ−1)`, whose largest real part crosses zero at `μ = 2` too — a
//!   DIFFERENT curve of `μ` reaching the same boundary, so the certificate is
//!   TIGHT here. Separately, [`fs_spectral`] recomputes the symmetric-part
//!   abscissa (the Lyapunov condition) so its agreement with `fs-sos` is an
//!   implementation cross-check of the two crates.
//! - **The solver** ([`fs_couple`]): to actually COMPUTE the coupled response by
//!   a partitioned scheme, naive staggering diverges early (`μ ≥ 1`), but Aitken
//!   relaxation converges across the whole physically-stable range up to `μ*`.
//! - **Honest colors** ([`fs_evidence`]): the only `Verified` payload this
//!   campaign mints is [`FlutterReport::witness_decay_rate_color`], and it names
//!   exactly one quantity — the LARGEST eigenvalue real part at `witness_mu` —
//!   enclosed by outward-rounded [`fs_ivl`] arithmetic
//!   ([`spectral_abscissa_interval`]). It is NOT an enclosure of the whole
//!   spectrum: for `μ > 1` the second eigenvalue's real part `−1 − √(μ−1)` lies
//!   strictly below the enclosure and is deliberately outside the claim.
//!
//! Deterministic; no dependencies beyond the composed crates.

use fs_couple::{iterate_aitken, iterate_fixed_relaxation};
use fs_evidence::Color;
use fs_ivl::Interval;
use fs_sos::lyapunov_certifies_stability;
use fs_spectral::symmetric_eigenvalues;

/// The 2-DOF coupled operator `A(μ)`.
#[must_use]
pub fn operator(mu: f64) -> [[f64; 2]; 2] {
    [[-1.0, 1.0], [mu - 1.0, -1.0]]
}

/// The NUMERICAL abscissa: largest eigenvalue of the symmetric part `(Aᵀ+A)/2`,
/// via [`fs_spectral`]. This is EXACTLY the Lyapunov-`P=I` condition
/// (`−(Aᵀ+A) ≻ 0 ⇔ this < 0`), so its agreement with the `fs-sos` certificate is
/// an IMPLEMENTATION cross-check between the two crates — not a new method.
#[must_use]
pub fn numerical_abscissa(mu: f64) -> f64 {
    let sym = vec![vec![-1.0, mu / 2.0], vec![mu / 2.0, -1.0]];
    symmetric_eigenvalues(&sym).map_or(f64::INFINITY, |e| {
        e.into_iter().fold(f64::NEG_INFINITY, f64::max)
    })
}

/// The SPECTRAL abscissa: largest real part of `A(μ)`'s ACTUAL eigenvalues
/// `−1 ± √(μ−1)` (trace `−2`, det `2−μ`). This is the necessary-AND-sufficient
/// asymptotic-stability criterion — genuinely independent of the (merely
/// sufficient) quadratic-Lyapunov certificate, so its boundary agreeing with the
/// Lyapunov boundary shows the `P=I` certificate is TIGHT for this operator.
///
/// This is a round-to-nearest `f64` DIAGNOSTIC: `√` is correctly rounded and the
/// following `−1 +` rounds again, so the returned value can sit up to about one
/// ulp on either side of the exact abscissa. It therefore carries no bound
/// authority — use [`spectral_abscissa_interval`] wherever a certified endpoint
/// is published.
#[must_use]
pub fn spectral_abscissa(mu: f64) -> f64 {
    // μ < 1 → complex pair with real part −1; μ ≥ 1 → real, max = −1 + √(μ−1).
    -1.0 + (mu - 1.0).max(0.0).sqrt()
}

/// A CERTIFIED enclosure of the LARGEST eigenvalue real part of `A(μ)` — the
/// asymptotic decay rate — computed as an outward-rounded [`fs_ivl`] evaluation
/// of `−1 + √(max(μ−1, 0))`.
///
/// Every step is outward-rounded, including the `μ − 1` subtraction, so the
/// returned interval encloses the abscissa both of the ideal operator `A(μ)` and
/// of the `f64` matrix [`operator`] actually builds (whose `μ−1` entry is the
/// nearest rounding of the same difference, hence inside the same enclosure).
///
/// The claim is deliberately narrow. It is the LARGEST real part only; `A(μ)`'s
/// other eigenvalue has real part `−1 − √(μ−1)`, which for `μ > 1` lies strictly
/// BELOW this interval and is not enclosed by it. A non-finite `μ` yields
/// [`Interval::WHOLE`] — the no-claim answer — rather than a fabricated bound.
#[must_use]
pub fn spectral_abscissa_interval(mu: f64) -> Interval {
    if !mu.is_finite() {
        return Interval::WHOLE;
    }
    // Outward-rounded μ − 1: encloses the exact difference and the rounded one.
    let shifted = Interval::point(mu) - Interval::point(1.0);
    // `max(·, 0)`: for μ < 1 the eigenvalues are a complex pair whose real part
    // is exactly −1, so the radical is exactly zero. `Interval::sqrt` already
    // clips a zero-straddling interval, but refuses an entirely negative one.
    let radicand = if shifted.hi() < 0.0 {
        Interval::point(0.0)
    } else {
        shifted
    };
    Interval::point(-1.0) + radicand.sqrt()
}

/// One sampled operating point.
#[derive(Debug, Clone, Copy)]
pub struct Sample {
    /// Added-mass ratio.
    pub mu: f64,
    /// Lyapunov certificate (fs-sos, `P=I`): is the system provably stable?
    pub lyapunov_stable: bool,
    /// Numerical abscissa (fs-spectral, symmetric part) — the Lyapunov condition,
    /// for the implementation cross-check.
    pub numerical_abscissa: f64,
    /// Spectral abscissa (A's actual eigenvalues) — the independent criterion
    /// (`< 0` ⇔ asymptotically stable).
    pub spectral_abscissa: f64,
    /// Did a naive staggered partitioned solve converge?
    pub naive_converged: bool,
    /// Did the Aitken-relaxed partitioned solve converge?
    pub aitken_converged: bool,
}

/// The campaign report.
#[derive(Debug, Clone)]
pub struct FlutterReport {
    /// All sampled operating points.
    pub samples: Vec<Sample>,
    /// The largest sampled `μ` the Lyapunov proof certifies stable. This is a
    /// diagnostic lower-side sample, not by itself a boundary location.
    pub lyapunov_boundary: f64,
    /// The largest sampled `μ` the independent eigenvalue criterion calls
    /// stable. This is likewise not a boundary location without a transition.
    pub eigen_boundary: f64,
    /// The shared ordered `[stable_sample, unstable_sample]` transition bracket
    /// when both criteria witness the same transition and agree at every
    /// sample. `None` means this sweep did not locate the boundary.
    pub boundary_bracket: Option<[f64; 2]>,
    /// Do the Lyapunov and independent eigenvalue stability classifications
    /// agree at every sample?
    pub stability_classifications_agree: bool,
    /// Did both criteria witness the same ordered stable-to-unstable bracket?
    /// Equality of two co-truncated stable maxima is never sufficient.
    pub boundaries_agree: bool,
    /// Do the fs-sos Lyapunov flag and fs-spectral numerical abscissa agree at
    /// every sample (an implementation cross-check of the SAME condition)?
    pub impl_consistent: bool,
    /// The largest `μ` a naive partitioned solve converged at.
    pub naive_boundary: f64,
    /// The largest `μ` the Aitken partitioned solve converged at.
    pub aitken_boundary: f64,
    /// Aitken converges strictly past the naive solver's reach.
    pub aitken_beats_naive: bool,
    /// A witness `μ` inside the certified-stable range where naive fails but
    /// Aitken succeeds — with the `Verified` Lyapunov color.
    pub witness_mu: Option<f64>,
    /// `Verified{lo, hi}` enclosing ONE named quantity: the LARGEST eigenvalue
    /// real part of `A(witness_mu)` — the asymptotic decay rate — as returned by
    /// [`spectral_abscissa_interval`], whose endpoints are outward-rounded.
    ///
    /// It is NOT an enclosure of the operator's spectrum: for `μ > 1` the second
    /// eigenvalue's real part `−1 − √(μ−1)` lies strictly below `lo`. Stability
    /// is read off `hi < 0`; the sample's separate `lyapunov_stable` flag is the
    /// `fs-sos` certificate and is not folded into these endpoints.
    pub witness_decay_rate_color: Option<Color>,
}

/// Run the FlutterCert sweep over `μ ∈ [lo, hi]` with `steps` points.
///
/// # Panics
/// If `steps < 2`.
#[must_use]
pub fn run_campaign(lo: f64, hi: f64, steps: usize) -> FlutterReport {
    assert!(steps >= 2, "need at least two samples");
    let mut samples = Vec::with_capacity(steps);
    for k in 0..steps {
        let mu = hi.mul_add(
            k as f64 / (steps - 1) as f64,
            lo * (1.0 - k as f64 / (steps - 1) as f64),
        );
        let lyapunov_stable = lyapunov_certifies_stability(operator(mu), [[1.0, 0.0], [0.0, 1.0]]);
        // Partitioned interface solves (fixed-point of H(x) = −μx + c). The
        // naive iteration gets a GENEROUS cap so its non-convergence reflects the
        // FUNDAMENTAL divergence for μ ≥ 1 (contraction factor μ), not a
        // budget-limited slow decay near μ = 1.
        let naive = iterate_fixed_relaxation(mu, 1.0, 0.0, 1.0, 20_000, 1e-9);
        let aitken = iterate_aitken(mu, 1.0, 0.0, 0.5, 2.0, 300, 1e-9);
        samples.push(Sample {
            mu,
            lyapunov_stable,
            numerical_abscissa: numerical_abscissa(mu),
            spectral_abscissa: spectral_abscissa(mu),
            naive_converged: naive.converged,
            aitken_converged: aitken.converged,
        });
    }

    let last_true = |pred: &dyn Fn(&Sample) -> bool| -> f64 {
        samples
            .iter()
            .filter(|s| pred(s))
            .map(|s| s.mu)
            .fold(f64::NEG_INFINITY, f64::max)
    };
    let lyapunov_boundary = last_true(&|s| s.lyapunov_stable);
    let eigen_boundary = last_true(&|s| s.spectral_abscissa < 0.0);
    let impl_consistent = samples
        .iter()
        .all(|s| s.lyapunov_stable == (s.numerical_abscissa < 0.0));
    let stability_classifications_agree = samples
        .iter()
        .all(|s| s.lyapunov_stable == (s.spectral_abscissa < 0.0));
    let transition_bracket = |pred: &dyn Fn(&Sample) -> bool| -> Option<[f64; 2]> {
        if !samples
            .windows(2)
            .all(|window| window[0].mu.is_finite() && window[0].mu < window[1].mu)
        {
            return None;
        }
        let first_unstable = samples.iter().position(|sample| !pred(sample))?;
        if first_unstable == 0
            || !samples[..first_unstable].iter().all(pred)
            || !samples[first_unstable..].iter().all(|sample| !pred(sample))
        {
            return None;
        }
        Some([samples[first_unstable - 1].mu, samples[first_unstable].mu])
    };
    let lyapunov_bracket = transition_bracket(&|s| s.lyapunov_stable);
    let eigen_bracket = transition_bracket(&|s| s.spectral_abscissa < 0.0);
    let boundary_bracket = match (lyapunov_bracket, eigen_bracket) {
        (Some(lyapunov), Some(eigen)) if stability_classifications_agree && lyapunov == eigen => {
            Some(lyapunov)
        }
        _ => None,
    };
    let naive_boundary = last_true(&|s| s.naive_converged);
    let aitken_boundary = last_true(&|s| s.aitken_converged);

    // A witness: certified stable, naive fails, Aitken succeeds. The Verified
    // band names exactly ONE quantity — the largest eigenvalue real part (the
    // asymptotic decay rate) at this μ — and both endpoints come from the
    // outward-rounded interval evaluation, never from the round-to-nearest
    // `spectral_abscissa` diagnostic. A non-finite enclosure mints no color.
    let witness_sample = samples
        .iter()
        .find(|s| s.lyapunov_stable && !s.naive_converged && s.aitken_converged);
    let witness = witness_sample.map(|s| s.mu);
    let witness_decay_rate_color = witness_sample.and_then(|s| {
        let enclosure = spectral_abscissa_interval(s.mu);
        (enclosure.lo().is_finite() && enclosure.hi().is_finite()).then(|| {
            // declared-color-ok: outward-rounded enclosure of the named decay rate; admitted only at a consumer's authority boundary (6pf9)
            Color::Verified {
                lo: enclosure.lo(),
                hi: enclosure.hi(),
            }
        })
    });

    FlutterReport {
        boundaries_agree: boundary_bracket.is_some(),
        boundary_bracket,
        stability_classifications_agree,
        impl_consistent,
        aitken_beats_naive: aitken_boundary > naive_boundary + 1e-9,
        lyapunov_boundary,
        eigen_boundary,
        naive_boundary,
        aitken_boundary,
        witness_mu: witness,
        witness_decay_rate_color,
        samples,
    }
}
