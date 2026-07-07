//! fs-cheb — compute with FUNCTIONS as values (plan §6.5): smooth 1D
//! functions as adaptively truncated Chebyshev expansions with automatic
//! near-machine-precision degree selection, plus spectral collocation
//! differentiation matrices.
//!
//! Representation: coefficients over FIRST-KIND Chebyshev points (the
//! roots grid) — chosen deliberately so values ↔ coefficients is exactly
//! fs-fft's DCT-II/III pair (cross-ISA bit-deterministic by construction).
//! The Lobatto/DCT-I flavor, 2D low-rank, Fourier-periodic, colleague-
//! matrix roots, and the Orr–Sommerfeld complex eigenproblem stack are
//! recorded follow-up scope (the last needs a nonsymmetric/complex
//! eigensolver that does not exist in-workspace yet).
//!
//! Determinism: sampling grids, plateau detection, Clenshaw evaluation,
//! and rootfinding subdivision are all fixed-order arithmetic on strict
//! kernels — NO platform libm in any path that feeds function state
//! (workspace contract rule).

pub mod orr_sommerfeld;

use fs_fft::{dct2, dct3};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// A smooth function on [a, b] as a truncated Chebyshev series
/// f(x) ≈ Σ' cₖ·Tₖ(t(x)) with t the affine pullback to [−1, 1]
/// (the k = 0 term is halved — the DCT-II convention).
#[derive(Debug, Clone)]
pub struct Cheb1 {
    a: f64,
    b: f64,
    /// Chebyshev coefficients, c[0] stored UN-halved (Clenshaw applies
    /// the ½ convention at evaluation).
    coeffs: Vec<f64>,
}

/// Relative plateau threshold for adaptive truncation. Sits ABOVE the
/// DCT rounding floor (~n·eps effects at large n): 10·2⁻⁵² — chasing the
/// floor itself inflates degrees ~20× for oscillatory functions (measured
/// during bring-up: sin(20x) resolved at 1090 instead of ~45).
const PLATEAU_REL: f64 = 2.2e-15;

impl Cheb1 {
    /// Build adaptively from a scalar function on [a, b]: sample at
    /// first-kind Chebyshev grids of doubling size until the trailing
    /// quarter of coefficients sits at the machine-precision plateau,
    /// then truncate. Panics (structured) if `max_degree` cannot resolve
    /// the function (non-smooth input is a modeling error here).
    #[must_use]
    pub fn build<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, max_degree: usize) -> Cheb1 {
        assert!(a < b, "domain must satisfy a < b (got [{a}, {b}])");
        let mut n = 16usize;
        loop {
            let coeffs = Self::coeffs_at(f, a, b, n);
            let maxc = coeffs
                .iter()
                .fold(0.0f64, |m, &c| m.max(c.abs()))
                .max(f64::MIN_POSITIVE);
            let tail = &coeffs[3 * n / 4..];
            if tail.iter().all(|&c| c.abs() <= PLATEAU_REL * maxc) {
                // Truncate at the last coefficient above the plateau.
                let keep = coeffs
                    .iter()
                    .rposition(|&c| c.abs() > PLATEAU_REL * maxc)
                    .map_or(1, |p| p + 1);
                return Cheb1 {
                    a,
                    b,
                    coeffs: coeffs[..keep].to_vec(),
                };
            }
            n *= 2;
            assert!(
                n <= max_degree,
                "function not resolved at degree {max_degree} on [{a}, {b}] \
                 (non-smooth or too oscillatory; raise max_degree or split the domain)"
            );
        }
    }

    /// Coefficients from n samples at first-kind points via DCT-II:
    /// cⱼ = (2/n)·Σₖ f(xₖ)·cos(πj(2k+1)/(2n)).
    fn coeffs_at<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, n: usize) -> Vec<f64> {
        let vals = sample_first_kind(f, a, b, n);
        let mut c = dct2(&vals);
        let scale = 2.0 / n as f64;
        for v in &mut c {
            *v *= scale;
        }
        c
    }

    /// Construct directly from coefficients (c[0] un-halved convention).
    #[must_use]
    pub fn from_coeffs(a: f64, b: f64, coeffs: Vec<f64>) -> Cheb1 {
        assert!(a < b, "domain must satisfy a < b");
        assert!(!coeffs.is_empty(), "need at least one coefficient");
        Cheb1 { a, b, coeffs }
    }

    /// Degree (number of retained coefficients − 1).
    #[must_use]
    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    /// The domain.
    #[must_use]
    pub fn domain(&self) -> (f64, f64) {
        (self.a, self.b)
    }

    /// Coefficient view (c[0] un-halved).
    #[must_use]
    pub fn coeffs(&self) -> &[f64] {
        &self.coeffs
    }

    /// Evaluate by Clenshaw recurrence (fixed order, fused).
    #[must_use]
    pub fn eval(&self, x: f64) -> f64 {
        let t = 2.0 * (x - self.a) / (self.b - self.a) - 1.0;
        let (mut b1, mut b2) = (0.0f64, 0.0f64);
        for &c in self.coeffs.iter().skip(1).rev() {
            let b0 = (2.0 * t).mul_add(b1, c - b2);
            b2 = b1;
            b1 = b0;
        }
        // Σ' convention: half the k = 0 coefficient.
        t.mul_add(b1, 0.5f64.mul_add(self.coeffs[0], -b2))
    }

    /// Derivative as a new Chebyshev object (coefficient recurrence with
    /// the domain chain rule).
    #[must_use]
    pub fn differentiate(&self) -> Cheb1 {
        let n = self.coeffs.len();
        if n == 1 {
            return Cheb1 {
                a: self.a,
                b: self.b,
                coeffs: vec![0.0],
            };
        }
        // Series with halved-c0 semantics: work on the "true" coefficients.
        let mut d = vec![0.0f64; n];
        // Standard recurrence: d[k-1] = d[k+1] + 2k·c[k] (true c series),
        // where the stored c[0] is un-halved but T0's coefficient never
        // enters the derivative sums.
        for k in (1..n).rev() {
            let above = if k + 2 < n { d[k + 1] } else { 0.0 };
            d[k - 1] = (2.0 * k as f64).mul_add(self.coeffs[k], above);
        }
        // Chain rule for [a,b] → factor 2/(b−a); d[0] doubles under the
        // Σ' storage convention (stored un-halved).
        let scale = 2.0 / (self.b - self.a);
        let mut out: Vec<f64> = d[..n - 1].iter().map(|&v| v * scale).collect();
        if out.is_empty() {
            out.push(0.0);
        }
        Cheb1 {
            a: self.a,
            b: self.b,
            coeffs: out,
        }
    }

    /// Definite integral over the whole domain: only even coefficients
    /// contribute (∫₋₁¹ Tₖ = 2/(1−k²) for even k, else 0).
    #[must_use]
    pub fn integral(&self) -> f64 {
        let mut acc = self.coeffs[0]; // (½·c0)·2 = c0 with stored-un-halved
        for (k, &c) in self.coeffs.iter().enumerate().skip(2).step_by(2) {
            acc += 2.0 * c / (1.0 - (k as f64) * (k as f64));
        }
        acc * (self.b - self.a) / 2.0
    }

    /// Sum of two functions on the same domain.
    #[must_use]
    pub fn add(&self, o: &Cheb1) -> Cheb1 {
        assert!(
            (self.a - o.a).abs() < 1e-14 && (self.b - o.b).abs() < 1e-14,
            "domain mismatch"
        );
        let n = self.coeffs.len().max(o.coeffs.len());
        let mut coeffs = vec![0.0f64; n];
        for (i, c) in coeffs.iter_mut().enumerate() {
            *c = self.coeffs.get(i).copied().unwrap_or(0.0)
                + o.coeffs.get(i).copied().unwrap_or(0.0);
        }
        Cheb1 {
            a: self.a,
            b: self.b,
            coeffs,
        }
    }

    /// Product via resampling at a grid resolving the sum of degrees.
    #[must_use]
    pub fn mul(&self, o: &Cheb1) -> Cheb1 {
        assert!(
            (self.a - o.a).abs() < 1e-14 && (self.b - o.b).abs() < 1e-14,
            "domain mismatch"
        );
        let n = (self.coeffs.len() + o.coeffs.len())
            .next_power_of_two()
            .max(16);
        let f = |x: f64| self.eval(x) * o.eval(x);
        let coeffs = Cheb1::coeffs_at(&f, self.a, self.b, n);
        // Truncate at plateau.
        let maxc = coeffs
            .iter()
            .fold(0.0f64, |m, &c| m.max(c.abs()))
            .max(f64::MIN_POSITIVE);
        let keep = coeffs
            .iter()
            .rposition(|&c| c.abs() > PLATEAU_REL * maxc)
            .map_or(1, |p| p + 1);
        Cheb1 {
            a: self.a,
            b: self.b,
            coeffs: coeffs[..keep].to_vec(),
        }
    }

    /// All real roots in [a, b] by recursive subdivision on sign changes
    /// of the interpolant with Newton polish. v1 limitation (documented):
    /// roots of even multiplicity (no sign change) are not found —
    /// colleague-matrix rootfinding joins the follow-up bead.
    #[must_use]
    pub fn roots(&self) -> Vec<f64> {
        let mut out = Vec::new();
        // Scan a fine deterministic grid for sign changes.
        let samples = (8 * self.coeffs.len()).max(64);
        let mut prev_x = self.a;
        let mut prev_v = self.eval(prev_x);
        for k in 1..=samples {
            let x = self.a + (self.b - self.a) * (k as f64) / (samples as f64);
            let v = self.eval(x);
            if prev_v == 0.0 {
                out.push(prev_x);
            } else if prev_v * v < 0.0 {
                out.push(self.bisect_newton(prev_x, x));
            }
            prev_x = x;
            prev_v = v;
        }
        if prev_v == 0.0 {
            out.push(prev_x);
        }
        out
    }

    fn bisect_newton(&self, mut lo: f64, mut hi: f64) -> f64 {
        let d = self.differentiate();
        for _ in 0..40 {
            let mid = f64::midpoint(lo, hi);
            let v = self.eval(mid);
            if v == 0.0 {
                return mid;
            }
            if self.eval(lo) * v < 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        // Newton polish from the bisection estimate.
        let mut x = f64::midpoint(lo, hi);
        for _ in 0..4 {
            let dv = d.eval(x);
            if dv == 0.0 {
                break;
            }
            x -= self.eval(x) / dv;
            x = x.clamp(self.a, self.b);
        }
        x
    }
}

/// Sample f at the n first-kind Chebyshev points mapped to [a, b],
/// ordered k = 0..n (xₖ = cos(π(k+½)/n) descending in t).
fn sample_first_kind<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, n: usize) -> Vec<f64> {
    (0..n)
        .map(|k| {
            let theta = std::f64::consts::PI * (k as f64 + 0.5) / (n as f64);
            let t = fs_math::det::cos(theta);
            let x = f64::midpoint(a, b) + t * (b - a) / 2.0;
            f(x)
        })
        .collect()
}

/// Synthesis: values at the n first-kind points from coefficients
/// (inverse of the analysis map; used by tests and resampling).
#[must_use]
pub fn values_from_coeffs(coeffs: &[f64], n: usize) -> Vec<f64> {
    // DCT-III with the k = 0 halving convention: dct3 already applies it.
    let mut padded = coeffs.to_vec();
    padded.resize(n, 0.0);
    dct3(&padded)
}

// ---------------------------------------------------------------------------
// Spectral collocation (Chebyshev–Lobatto differentiation matrices)
// ---------------------------------------------------------------------------

/// The n+1 Chebyshev–Lobatto points on [−1, 1], DESCENDING (x₀ = 1),
/// the classical collocation ordering.
#[must_use]
pub fn lobatto_points(n: usize) -> Vec<f64> {
    (0..=n)
        .map(|j| fs_math::det::cos(std::f64::consts::PI * (j as f64) / (n as f64)))
        .collect()
}

/// The (n+1)×(n+1) first-derivative collocation matrix on the Lobatto
/// grid (Trefethen's construction, with the NEGATIVE-SUM TRICK on the
/// diagonal — the classic accuracy fix: rows must sum to zero exactly
/// because differentiation annihilates constants).
#[must_use]
pub fn diff_matrix(n: usize) -> Vec<f64> {
    assert!(n >= 1, "need at least two points");
    let x = lobatto_points(n);
    let m = n + 1;
    let c = |i: usize| -> f64 {
        let ci = if i == 0 || i == n { 2.0 } else { 1.0 };
        if i.is_multiple_of(2) { ci } else { -ci }
    };
    let mut d = vec![0.0f64; m * m];
    for i in 0..m {
        for j in 0..m {
            if i != j {
                d[i * m + j] = (c(i) / c(j)) / (x[i] - x[j]);
            }
        }
    }
    // Negative-sum trick: D[i][i] = −Σ_{j≠i} D[i][j].
    for i in 0..m {
        let mut s = 0.0f64;
        for j in 0..m {
            if j != i {
                s += d[i * m + j];
            }
        }
        d[i * m + i] = -s;
    }
    d
}

/// Smallest `k` eigenvalues of the Dirichlet problem −u″ = λu on [−1, 1]
/// by collocation: interior block of −D², solved by SHIFT-INVERTED power
/// iteration. Shifts come from a coarse FINITE-DIFFERENCE surrogate (a
/// symmetric tridiagonal whose spectrum `fs_la::eigen::jacobi_eigh`
/// handles) — deterministic, independent of the analytic answer, and it
/// sidesteps the missing general nonsymmetric eigensolver (that solver
/// is the Orr–Sommerfeld follow-up's first deliverable).
#[must_use]
pub fn dirichlet_laplace_eigs(n: usize, k: usize) -> Vec<f64> {
    let m = n + 1;
    let d = diff_matrix(n);
    // D² then take the interior (n−1)×(n−1) block, negated.
    let mut d2 = vec![0.0f64; m * m];
    for i in 0..m {
        for j in 0..m {
            let mut acc = 0.0f64;
            for l in 0..m {
                acc = d[i * m + l].mul_add(d[l * m + j], acc);
            }
            d2[i * m + j] = acc;
        }
    }
    let ni = n - 1;
    let mut a = vec![0.0f64; ni * ni];
    for i in 0..ni {
        for j in 0..ni {
            a[i * ni + j] = -d2[(i + 1) * m + (j + 1)];
        }
    }
    // FD surrogate on a uniform interior grid: (-1, 2, -1)/h² tridiag —
    // symmetric, so the landed dense Jacobi handles it. Its k smallest
    // eigenvalues approximate the true ones well enough to be shifts.
    let nf = 64usize;
    let h = 2.0 / (nf as f64 + 1.0);
    let mut fd = vec![0.0f64; nf * nf];
    for i in 0..nf {
        fd[i * nf + i] = 2.0 / (h * h);
        if i + 1 < nf {
            fd[i * nf + i + 1] = -1.0 / (h * h);
            fd[(i + 1) * nf + i] = -1.0 / (h * h);
        }
    }
    let (fd_eigs, _) = fs_la::eigen::jacobi_eigh(&fd, nf);
    let mut eigs = Vec::with_capacity(k);
    for &fd_est in fd_eigs.iter().take(k) {
        // Shift slightly BELOW the surrogate estimate (FD underestimates
        // continuum eigenvalues; the offset keeps the shifted matrix
        // definite and the iteration locked to the intended eigenvalue).
        let mu = fd_est * 0.95;
        let mut shifted = a.clone();
        for i in 0..ni {
            shifted[i * ni + i] -= mu;
        }
        let lu =
            fs_la::factor::lu(&shifted, ni).expect("shifted collocation operator is nonsingular");
        let mut v: Vec<f64> = (0..ni)
            .map(|i| 1.0 + 0.25 * (((i * 7 + 3) % 11) as f64))
            .collect();
        for _ in 0..100 {
            let nrm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            for x in &mut v {
                *x /= nrm;
            }
            lu.solve(&mut v);
        }
        // Rayleigh quotient λ = vᵀAv / vᵀv on the UNSHIFTED operator.
        let nrm2: f64 = v.iter().map(|x| x * x).sum();
        let mut av = vec![0.0f64; ni];
        for i in 0..ni {
            let mut acc = 0.0f64;
            for j in 0..ni {
                acc = a[i * ni + j].mul_add(v[j], acc);
            }
            av[i] = acc;
        }
        eigs.push(v.iter().zip(&av).map(|(x, y)| x * y).sum::<f64>() / nrm2);
    }
    eigs
}
