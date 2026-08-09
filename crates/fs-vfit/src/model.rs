//! Pole-residue rational models with a conjugate-closed-by-construction
//! representation and a real block state-space realization.
//!
//! Convention: Laplace variable `s`; frequency responses are sampled on
//! the positive imaginary axis `s = i·omega` (engineering `j·omega`).
//! A model is REAL (impulse response real) because complex poles are
//! stored only as the `Im > 0` member of a conjugate pair — the
//! conjugate term is implied, never stored, so conjugacy cannot drift.

use fs_math::c64::C64;

/// One real-rational term of a pole-residue expansion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PoleTerm {
    /// `residue / (s - pole)` with both entries real.
    Real {
        /// Real pole location (stable iff `pole < 0`).
        pole: f64,
        /// Real residue.
        residue: f64,
    },
    /// `residue/(s - pole) + conj(residue)/(s - conj(pole))` — the
    /// stored member has `pole.im > 0`; the conjugate is implied.
    Pair {
        /// Pole with strictly positive imaginary part.
        pole: C64,
        /// Residue attached to the stored (`Im > 0`) pole.
        residue: C64,
    },
}

impl PoleTerm {
    /// State dimension this term contributes (1 for real, 2 for pair).
    #[must_use]
    pub fn state_dim(&self) -> usize {
        match self {
            PoleTerm::Real { .. } => 1,
            PoleTerm::Pair { .. } => 2,
        }
    }

    /// Term value at complex frequency `s` (pair terms include BOTH
    /// conjugate members).
    #[must_use]
    pub fn eval(&self, s: C64) -> C64 {
        match *self {
            PoleTerm::Real { pole, residue } => {
                C64::from_re(residue) * (s - C64::from_re(pole)).recip()
            }
            PoleTerm::Pair { pole, residue } => {
                residue * (s - pole).recip() + residue.conj() * (s - pole.conj()).recip()
            }
        }
    }

    /// True iff every pole of the term is strictly in the open left
    /// half-plane.
    #[must_use]
    pub fn is_stable(&self) -> bool {
        match self {
            PoleTerm::Real { pole, .. } => *pole < 0.0,
            PoleTerm::Pair { pole, .. } => pole.re < 0.0,
        }
    }
}

/// A real rational model `H(s) = sum(terms) + d + s*e`.
///
/// `d` is the direct feedthrough and `e` the improper linear term. For
/// IMPEDANCE-form passivity both must be non-negative: `Re(i*omega*e) =
/// 0`, so the `e` term is LOSSLESS on the imaginary axis (an ideal
/// inductance) and contributes nothing to `Re H` — see
/// `passivity` for the descriptor-form statement.
#[derive(Debug, Clone, PartialEq)]
pub struct RationalModel {
    /// Pole-residue terms; pair terms count twice toward the order.
    pub terms: Vec<PoleTerm>,
    /// Direct (constant) term.
    pub d: f64,
    /// Improper linear term (`s*e`).
    pub e: f64,
}

/// Real state-space realization `x' = A x + B u`, `y = C x + D u +
/// E u'` (row-major `A`, `n = A` side). `E` carries the improper term.
#[derive(Debug, Clone, PartialEq)]
pub struct StateSpace {
    /// State dimension.
    pub n: usize,
    /// Row-major `n x n` block-diagonal state matrix.
    pub a: Vec<f64>,
    /// Input map, length `n`.
    pub b: Vec<f64>,
    /// Output map, length `n`.
    pub c: Vec<f64>,
    /// Direct feedthrough.
    pub d: f64,
    /// Improper (derivative feedthrough) term.
    pub e: f64,
}

impl RationalModel {
    /// Model order: total pole count (pairs count 2).
    #[must_use]
    pub fn order(&self) -> usize {
        self.terms.iter().map(PoleTerm::state_dim).sum()
    }

    /// Evaluate `H(s)`.
    #[must_use]
    pub fn eval(&self, s: C64) -> C64 {
        let mut acc = C64::from_re(self.d) + s.scale(self.e);
        for t in &self.terms {
            acc = acc + t.eval(s);
        }
        acc
    }

    /// Evaluate on the imaginary axis at angular frequency `omega`.
    #[must_use]
    pub fn eval_iw(&self, omega: f64) -> C64 {
        self.eval(C64::new(0.0, omega))
    }

    /// All poles, expanded (pairs contribute both members), in the
    /// solver's canonical `(re, im)` order.
    #[must_use]
    pub fn poles_expanded(&self) -> Vec<C64> {
        let mut out = Vec::with_capacity(self.order());
        for t in &self.terms {
            match *t {
                PoleTerm::Real { pole, .. } => out.push(C64::from_re(pole)),
                PoleTerm::Pair { pole, .. } => {
                    out.push(pole.conj());
                    out.push(pole);
                }
            }
        }
        out.sort_by(|x, y| x.re.total_cmp(&y.re).then(x.im.total_cmp(&y.im)));
        out
    }

    /// True iff every pole is strictly in the open left half-plane.
    #[must_use]
    pub fn is_stable(&self) -> bool {
        self.terms.iter().all(PoleTerm::is_stable)
    }

    /// Real block-diagonal state-space realization.
    ///
    /// Real pole `p`, residue `r`: 1x1 block `[p]`, `b = 1`, `c = r`.
    /// Pair `p = a + i*b_im`, `r = rho + i*sigma`: 2x2 block
    /// `[[a, b_im], [-b_im, a]]`, input `[1, 0]`, output
    /// `[2*rho, 2*sigma]`, because
    /// `r/(s-p) + conj` = `(2*rho*(s-a) - 2*sigma*b_im) / ((s-a)^2 + b_im^2)`
    /// and `c (sI-A)^{-1} b = (c1*(s-a) - c2*b_im) / ((s-a)^2 + b_im^2)`,
    /// so `c2 = +2*sigma` (verified against the independent LU route by
    /// the realization-parity test).
    #[must_use]
    pub fn state_space(&self) -> StateSpace {
        let n = self.order();
        let mut a = vec![0.0f64; n * n];
        let mut b = vec![0.0f64; n];
        let mut c = vec![0.0f64; n];
        let mut at = 0usize;
        for t in &self.terms {
            match *t {
                PoleTerm::Real { pole, residue } => {
                    a[at * n + at] = pole;
                    b[at] = 1.0;
                    c[at] = residue;
                    at += 1;
                }
                PoleTerm::Pair { pole, residue } => {
                    let (re, im) = (pole.re, pole.im);
                    a[at * n + at] = re;
                    a[at * n + at + 1] = im;
                    a[(at + 1) * n + at] = -im;
                    a[(at + 1) * n + at + 1] = re;
                    b[at] = 1.0;
                    b[at + 1] = 0.0;
                    c[at] = 2.0 * residue.re;
                    c[at + 1] = 2.0 * residue.im;
                    at += 2;
                }
            }
        }
        StateSpace {
            n,
            a,
            b,
            c,
            d: self.d,
            e: self.e,
        }
    }
}

impl StateSpace {
    /// Evaluate `C (sI - A)^{-1} B + D + s E` by a dense complex solve —
    /// the INDEPENDENT route the realization-parity test uses against
    /// [`RationalModel::eval`].
    ///
    /// # Errors
    /// Propagates the fs-la singularity signal when `sI - A` is
    /// singular (i.e. `s` is an eigenvalue of `A`).
    pub fn eval(&self, s: C64) -> Result<C64, fs_la::eigen_complex::EigFailure> {
        let n = self.n;
        let mut m = vec![C64::ZERO; n * n];
        for i in 0..n {
            for j in 0..n {
                m[i * n + j] = C64::from_re(-self.a[i * n + j]);
            }
            m[i * n + i] = m[i * n + i] + s;
        }
        let lu = fs_la::eigen_complex::lu_complex(&m, n)?;
        let mut x: Vec<C64> = self.b.iter().map(|&v| C64::from_re(v)).collect();
        lu.solve(&mut x);
        let mut acc = C64::from_re(self.d) + s.scale(self.e);
        for (ci, xi) in self.c.iter().zip(&x) {
            acc = acc + xi.scale(*ci);
        }
        Ok(acc)
    }
}
