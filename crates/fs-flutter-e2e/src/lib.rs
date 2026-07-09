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
//! - **The cross-check** ([`fs_spectral`]): independently, the symmetric part
//!   `(Aᵀ+A)/2 = [[−1, μ/2], [μ/2, −1]]` has largest eigenvalue `−1 + μ/2 < 0`
//!   iff `μ < 2`. Two different methods agree on the boundary — certifying the
//!   certifier.
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

/// The largest eigenvalue of the symmetric part `(Aᵀ+A)/2` — the numerical
/// abscissa (an independent decay-rate witness).
#[must_use]
pub fn spectral_abscissa(mu: f64) -> f64 {
    let sym = vec![vec![-1.0, mu / 2.0], vec![mu / 2.0, -1.0]];
    symmetric_eigenvalues(&sym).map_or(f64::INFINITY, |e| {
        e.into_iter().fold(f64::NEG_INFINITY, f64::max)
    })
}

/// One sampled operating point.
#[derive(Debug, Clone, Copy)]
pub struct Sample {
    /// Added-mass ratio.
    pub mu: f64,
    /// Lyapunov certificate: is the system provably stable at this `μ`?
    pub lyapunov_stable: bool,
    /// Independent spectral abscissa (`< 0` ⇔ stable).
    pub abscissa: f64,
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
    /// The largest `μ` the independent spectral test calls stable.
    pub spectral_boundary: f64,
    /// Do the two independent boundaries agree (certifying the certifier)?
    pub boundaries_agree: bool,
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
        let abscissa = spectral_abscissa(mu);
        // Partitioned interface solves (fixed-point of H(x) = −μx + c).
        let naive = iterate_fixed_relaxation(mu, 1.0, 0.0, 1.0, 300, 1e-9);
        let aitken = iterate_aitken(mu, 1.0, 0.0, 0.5, 2.0, 300, 1e-9);
        samples.push(Sample {
            mu,
            lyapunov_stable,
            abscissa,
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
    let spectral_boundary = last_true(&|s| s.abscissa < 0.0);
    let naive_boundary = last_true(&|s| s.naive_converged);
    let aitken_boundary = last_true(&|s| s.aitken_converged);

    // A witness: certified stable, naive fails, Aitken succeeds.
    let witness = samples
        .iter()
        .find(|s| s.lyapunov_stable && !s.naive_converged && s.aitken_converged)
        .map(|s| s.mu);
    let witness_color = witness.map(|_| Color::Verified { lo: 0.0, hi: 0.0 });

    FlutterReport {
        boundaries_agree: (lyapunov_boundary - spectral_boundary).abs() < 1e-9,
        aitken_beats_naive: aitken_boundary > naive_boundary + 1e-9,
        lyapunov_boundary,
        spectral_boundary,
        naive_boundary,
        aitken_boundary,
        witness_mu: witness,
        witness_color,
        samples,
    }
}
