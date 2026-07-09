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
//! - **Honest colors** ([`fs_evidence`]): the Lyapunov certificate is `Verified`.
//!
//! Deterministic; no dependencies beyond the composed crates.

use fs_couple::{iterate_aitken, iterate_fixed_relaxation};
use fs_evidence::Color;
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
#[must_use]
pub fn spectral_abscissa(mu: f64) -> f64 {
    // μ < 1 → complex pair with real part −1; μ ≥ 1 → real, max = −1 + √(μ−1).
    -1.0 + (mu - 1.0).max(0.0).sqrt()
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
    /// The largest `μ` the Lyapunov proof certifies stable (the proven boundary).
    pub lyapunov_boundary: f64,
    /// The largest `μ` the INDEPENDENT eigenvalue criterion calls stable.
    pub eigen_boundary: f64,
    /// Do the Lyapunov (sufficient) and eigenvalue (necessary+sufficient)
    /// boundaries agree — the P=I certificate is tight (a genuine cross-check)?
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
    /// The witness's certificate color (`Verified`).
    pub witness_color: Option<Color>,
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
    let naive_boundary = last_true(&|s| s.naive_converged);
    let aitken_boundary = last_true(&|s| s.aitken_converged);

    // A witness: certified stable, naive fails, Aitken succeeds. Its Verified
    // band is the spectral abscissa — a certified negative upper bound on the
    // eigenvalues' real parts (the exponential decay rate).
    let witness_sample = samples
        .iter()
        .find(|s| s.lyapunov_stable && !s.naive_converged && s.aitken_converged);
    let witness = witness_sample.map(|s| s.mu);
    let witness_color = witness_sample.map(|s| Color::Verified {
        lo: s.spectral_abscissa,
        hi: 0.0,
    });

    FlutterReport {
        boundaries_agree: (lyapunov_boundary - eigen_boundary).abs() < 1e-9,
        impl_consistent,
        aitken_beats_naive: aitken_boundary > naive_boundary + 1e-9,
        lyapunov_boundary,
        eigen_boundary,
        naive_boundary,
        aitken_boundary,
        witness_mu: witness,
        witness_color,
        samples,
    }
}
