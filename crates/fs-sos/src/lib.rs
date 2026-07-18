//! fs-sos — proof-carrying optimization (sum-of-squares certificates). Layer: L4.
//!
//! A lower bound proven by SAMPLING can be wrong (it can miss the true minimum).
//! A SUM-OF-SQUARES identity `p(x) − γ = Σ qᵢ(x)²` instead proves `p(x) ≥ γ`.
//! This crate separates a tolerance-based coefficient diagnostic
//! ([`SosCertificate::verify`]) from APIs that may make value claims:
//! [`SosCertificate::certified_bound_global`] admits only exact non-constant
//! residual cancellation, while [`SosCertificate::certified_bound_on`] encloses
//! residual error on one finite radius.
//!
//! Included: univariate [`Poly`] arithmetic; [`certify_quadratic`] (a
//! completed-square `f64` result + its SOS certificate); [`is_psd`] (the
//! SDP-feasibility core, by an in-house Jacobi eigensolver); and
//! [`lyapunov_certifies_stability`] (an SOS/quadratic Lyapunov stability
//! certificate for a linear system).
//!
//! Deterministic for a fixed build and ISA; exact residual accounting uses
//! `fs-ivl` expansions and `fs-math` error-free transforms.

/// A univariate polynomial, coefficients ascending (`coeffs[i]` multiplies `xⁱ`).
#[derive(Debug, Clone, PartialEq)]
pub struct Poly {
    coeffs: Vec<f64>,
}

impl Poly {
    /// A polynomial from ascending coefficients (trailing near-zeros trimmed).
    #[must_use]
    pub fn new(coeffs: Vec<f64>) -> Poly {
        let mut p = Poly { coeffs };
        p.trim();
        p
    }

    /// A constant polynomial.
    #[must_use]
    pub fn constant(c: f64) -> Poly {
        Poly::new(vec![c])
    }

    fn trim(&mut self) {
        while self.coeffs.len() > 1 && self.coeffs.last().copied().unwrap_or(0.0).abs() < 1e-15 {
            self.coeffs.pop();
        }
    }

    /// The coefficients (ascending).
    #[must_use]
    pub fn coeffs(&self) -> &[f64] {
        &self.coeffs
    }

    /// The degree.
    #[must_use]
    pub fn degree(&self) -> usize {
        self.coeffs.len().saturating_sub(1)
    }

    /// Evaluate at `x` (Horner).
    #[must_use]
    pub fn eval(&self, x: f64) -> f64 {
        self.coeffs.iter().rev().fold(0.0, |acc, &c| acc * x + c)
    }

    /// Sum.
    #[must_use]
    pub fn add(&self, other: &Poly) -> Poly {
        let n = self.coeffs.len().max(other.coeffs.len());
        let coeffs = (0..n)
            .map(|i| {
                self.coeffs.get(i).copied().unwrap_or(0.0)
                    + other.coeffs.get(i).copied().unwrap_or(0.0)
            })
            .collect();
        Poly::new(coeffs)
    }

    /// Difference.
    #[must_use]
    pub fn sub(&self, other: &Poly) -> Poly {
        let n = self.coeffs.len().max(other.coeffs.len());
        let coeffs = (0..n)
            .map(|i| {
                self.coeffs.get(i).copied().unwrap_or(0.0)
                    - other.coeffs.get(i).copied().unwrap_or(0.0)
            })
            .collect();
        Poly::new(coeffs)
    }

    /// Product.
    #[must_use]
    pub fn mul(&self, other: &Poly) -> Poly {
        if self.coeffs.is_empty() || other.coeffs.is_empty() {
            return Poly::constant(0.0);
        }
        let mut coeffs = vec![0.0; self.coeffs.len() + other.coeffs.len() - 1];
        for (i, &a) in self.coeffs.iter().enumerate() {
            for (j, &b) in other.coeffs.iter().enumerate() {
                coeffs[i + j] += a * b;
            }
        }
        Poly::new(coeffs)
    }

    /// The largest absolute coefficient (a supremum-style norm; `0` for the zero
    /// polynomial).
    #[must_use]
    pub fn max_abs_coeff(&self) -> f64 {
        self.coeffs.iter().map(|c| c.abs()).fold(0.0, f64::max)
    }
}

/// The square `q·q`.
#[must_use]
pub fn square(q: &Poly) -> Poly {
    q.mul(q)
}

/// Advance a non-negative residual magnitude four ULPs toward +∞ without
/// crossing `f64::MAX` into a NaN bit pattern.
fn outward_nudge_magnitude(mut magnitude: f64) -> f64 {
    if magnitude.is_nan() {
        return f64::INFINITY;
    }
    for _ in 0..4 {
        magnitude = fs_math::next_up(magnitude);
    }
    magnitude
}

/// A sum-of-squares certificate: the claim `p(x) − lower_bound = Σ squaresᵢ(x)²`,
/// which (if it holds as a polynomial identity) PROVES `p(x) ≥ lower_bound`.
#[derive(Debug, Clone, PartialEq)]
pub struct SosCertificate {
    /// The polynomials whose squares sum to `p − lower_bound`.
    pub squares: Vec<Poly>,
    /// The certified lower bound.
    pub lower_bound: f64,
}

impl SosCertificate {
    /// The certificate residual: the largest absolute coefficient of
    /// `p − lower_bound − Σ squaresᵢ²` (nominally `0`).
    #[must_use]
    pub fn residual(&self, p: &Poly) -> f64 {
        let mut acc = p.sub(&Poly::constant(self.lower_bound));
        for q in &self.squares {
            acc = acc.sub(&square(q));
        }
        acc.max_abs_coeff()
    }

    /// Does the DECOMPOSITION IDENTITY hold within `tol` in coefficient
    /// ∞-norm? This is a statement about the algebra of the certificate,
    /// NOT a value bound: a residual with tiny coefficients can be
    /// arbitrarily negative far from the origin (r(x) = −tol − tol·x is
    /// −tol·(1+10⁶) at x = 10⁶), so passing `verify` must never be
    /// converted into a global `p ≥ lower_bound` claim (bead wa8i S2 —
    /// the old `certified_bound(tol)` did exactly that and forged
    /// certificates). Use [`SosCertificate::certified_bound_global`] or
    /// [`SosCertificate::certified_bound_on`] for value claims.
    #[must_use]
    pub fn verify(&self, p: &Poly, tol: f64) -> bool {
        self.residual(p) <= tol
    }

    /// The EXACT residual coefficients of `p − lower_bound − Σ qᵢ²` as
    /// fs-ivl expansions: every product and sum of f64 coefficients is
    /// exactly representable, so these are the true residual
    /// coefficients with NO rounding (the substrate for rigorous value
    /// claims).
    fn residual_expansions(&self, p: &Poly) -> Vec<Vec<f64>> {
        use fs_ivl::expansion::expansion_diff;
        use fs_math::eft::two_prod;
        // Degree bound of the residual.
        let deg_p = p.coeffs().len();
        let deg_sq = self
            .squares
            .iter()
            .map(|q| {
                let d = q.coeffs().len();
                if d == 0 { 0 } else { 2 * d - 1 }
            })
            .max()
            .unwrap_or(0);
        let len = deg_p.max(deg_sq).max(1);
        let mut acc: Vec<Vec<f64>> = (0..len)
            .map(|j| {
                let mut c = p.coeffs().get(j).copied().unwrap_or(0.0);
                if j == 0 {
                    // p₀ − lower_bound must itself be exact: use a
                    // two-component difference expansion.
                    let (hi, lo) = fs_ivl::expansion::two_diff(c, self.lower_bound);
                    return vec![lo, hi];
                }
                if c == 0.0 {
                    c = 0.0; // normalize −0
                }
                vec![c]
            })
            .collect();
        for q in &self.squares {
            let qc = q.coeffs();
            for (i, &qi) in qc.iter().enumerate() {
                for (j, &qj) in qc.iter().enumerate() {
                    let (hi, lo) = two_prod(qi, qj);
                    let prod = vec![lo, hi];
                    let k = i + j;
                    acc[k] = expansion_diff(&acc[k], &prod);
                }
            }
        }
        acc
    }

    /// GLOBAL certified bound: `Some(b)` only when every non-constant
    /// residual coefficient is EXACTLY zero (expansion sign 0), in which
    /// case `p − lower_bound = Σ qᵢ² + r₀` identically and
    /// `p(x) ≥ lower_bound − |r₀|` for ALL real x — a theorem, not a
    /// tolerance. Returns the bound degraded by the exact |r₀| envelope.
    #[must_use]
    pub fn certified_bound_global(&self, p: &Poly) -> Option<f64> {
        use fs_ivl::expansion::{estimate, expansion_sign};
        let res = self.residual_expansions(p);
        for e in res.iter().skip(1) {
            if expansion_sign(e) != 0 {
                return None;
            }
        }
        // |r₀| upper bound: expansions are exact; Σ|components| rounded
        // up encloses the magnitude.
        let r0 = res.first().map_or(0.0, |e| {
            let approx = estimate(e).abs();
            // A few-ulp outward nudge covers the summation rounding.
            outward_nudge_magnitude(approx)
        });
        Some(self.lower_bound - r0)
    }

    /// RADIUS-SCOPED certified bound: for any finite `radius > 0`,
    /// `p(x) ≥ returned` for all |x| ≤ radius — true for EVERY
    /// certificate by the triangle inequality
    /// (p = lb + Σq² + r ≥ lb − Σⱼ|rⱼ|·Rʲ), with the envelope computed
    /// through outward-rounded interval arithmetic. A lying certificate
    /// does not break soundness here: its large residual simply degrades
    /// the returned bound into uselessness (the lie becomes visible as
    /// distance from `lower_bound`, never as a false claim).
    #[must_use]
    pub fn certified_bound_on(&self, p: &Poly, radius: f64) -> Option<f64> {
        use fs_ivl::Interval;
        use fs_ivl::expansion::estimate;
        if !(radius.is_finite() && radius > 0.0) {
            return None;
        }
        let res = self.residual_expansions(p);
        let r_iv = Interval::new(radius, radius);
        let mut envelope = Interval::new(0.0, 0.0);
        let mut pow = Interval::new(1.0, 1.0);
        for (j, e) in res.iter().enumerate() {
            if j > 0 {
                pow = pow * r_iv;
            }
            let mag = estimate(e).abs();
            // A few-ulp outward nudge on the exact-expansion magnitude.
            let mag_up = outward_nudge_magnitude(mag);
            envelope = envelope + Interval::new(mag_up, mag_up) * pow;
        }
        Some(self.lower_bound - envelope.hi())
    }
}

/// The completed-square `f64` result for `a·x² + b·x + c` (`a > 0`) with its
/// SOS certificate:
/// `p(x) − (c − b²/4a) = (√a·x + b/2√a)²` in exact arithmetic.
///
/// Returns `None` when the coefficients or derived certificate values are
/// non-finite, or when `a <= 0` (not bounded below by a square).
#[must_use]
pub fn certify_quadratic(a: f64, b: f64, c: f64) -> Option<SosCertificate> {
    if !(a.is_finite() && a > 0.0 && b.is_finite() && c.is_finite()) {
        return None;
    }
    let lower_bound = c - b * b / (4.0 * a);
    let root_a = a.sqrt();
    // q(x) = √a·x + b/(2√a).
    let q0 = b / (2.0 * root_a);
    if !(lower_bound.is_finite() && root_a.is_finite() && q0.is_finite()) {
        return None;
    }
    let q = Poly::new(vec![q0, root_a]);
    Some(SosCertificate {
        squares: vec![q],
        lower_bound,
    })
}

/// Is the symmetric matrix positive semidefinite (min eigenvalue `>= −tol`)?
/// The feasibility core of the SDP the full Lasserre hierarchy solves.
#[must_use]
pub fn is_psd(matrix: &[Vec<f64>], tol: f64) -> bool {
    min_eigenvalue(matrix) >= -tol
}

/// Does the quadratic Lyapunov function `V(x) = xᵀPx` certify asymptotic
/// stability of the 2-D linear system `ẋ = Ax`? True iff `P ≻ 0` and
/// `−(AᵀP + PA) ≻ 0` (Lyapunov's theorem) — a sound SOS/quadratic stability
/// certificate. Finding such a `P` is the SDP (staged); this VERIFIES a candidate.
#[must_use]
pub fn lyapunov_certifies_stability(a: [[f64; 2]; 2], p: [[f64; 2]; 2]) -> bool {
    let pm = vec![vec![p[0][0], p[0][1]], vec![p[1][0], p[1][1]]];
    // Aᵀ P + P A.
    let at = [[a[0][0], a[1][0]], [a[0][1], a[1][1]]];
    let atp = matmul2(at, p);
    let pa = matmul2(p, a);
    let q = [
        [-(atp[0][0] + pa[0][0]), -(atp[0][1] + pa[0][1])],
        [-(atp[1][0] + pa[1][0]), -(atp[1][1] + pa[1][1])],
    ];
    let qm = vec![vec![q[0][0], q[0][1]], vec![q[1][0], q[1][1]]];
    // strict definiteness: min eigenvalue > 0 (small positive threshold).
    min_eigenvalue(&pm) > 1e-9 && min_eigenvalue(&qm) > 1e-9
}

fn matmul2(a: [[f64; 2]; 2], b: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
    [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ],
    ]
}

/// The smallest eigenvalue of a symmetric matrix, by cyclic Jacobi rotations.
// A dense symmetric eigen-kernel: `m[i][j]` is inherently 2D-indexed, so the
// index loops are the correct, readable form.
#[allow(clippy::needless_range_loop)]
fn min_eigenvalue(a: &[Vec<f64>]) -> f64 {
    let n = a.len();
    if n == 0 {
        return 0.0;
    }
    // Operate on the SYMMETRIC part (M+Mᵀ)/2. A quadratic form xᵀMx equals
    // xᵀ((M+Mᵀ)/2)x, and the cyclic-Jacobi kernel below is only valid for a
    // symmetric matrix — it reads the upper triangle (`m[p][q]`) and never
    // consults `m[q][p]`. Without this projection a NON-symmetric input yields
    // meaningless "eigenvalues" (e.g. it returns the raw diagonal when the
    // upper triangle is zero) and could FORGE a PSD / Lyapunov-stability
    // certificate — see `non_symmetric_input_cannot_forge_a_certificate`.
    // `midpoint(x, x) == x` exactly, so a symmetric input is unchanged bitwise.
    let mut m: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| f64::midpoint(a[i][j], a[j][i])).collect())
        .collect();
    for _ in 0..100 {
        let mut off = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                off += m[i][j] * m[i][j];
            }
        }
        if off <= 1e-28 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if m[p][q].abs() <= 1e-20 {
                    continue;
                }
                let theta = (m[q][q] - m[p][p]) / (2.0 * m[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (mkp, mkq) = (m[k][p], m[k][q]);
                    m[k][p] = c * mkp - s * mkq;
                    m[k][q] = s * mkp + c * mkq;
                }
                for k in 0..n {
                    let (mpk, mqk) = (m[p][k], m[q][k]);
                    m[p][k] = c * mpk - s * mqk;
                    m[q][k] = s * mpk + c * mqk;
                }
            }
        }
    }
    (0..n).map(|i| m[i][i]).fold(f64::INFINITY, f64::min)
}
