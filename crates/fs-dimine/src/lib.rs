//! fs-dimine — dimensional knowledge mining (plan addendum, Proposal 9's
//! knowledge apex). Layer: L4.
//!
//! The certified corpus is not just answers — it is DATA. With every quantity
//! reduced to dimensionless groups (Buckingham-π, done upstream by fs-regime),
//! automated mining fits closed-form SCALING LAWS `y = C · Π πⱼ^{aⱼ}` and lets
//! the system accumulate engineering knowledge with pedigrees. Humans never do
//! this systematically because it is tedious; a swarm does it overnight.
//!
//! This crate fits a power law by PURE-RUST log-linear least squares (no
//! external symbolic-regression library, no Python — Franken-only): taking
//! logs turns `y = C · Π πⱼ^{aⱼ}` into a linear model `ln y = ln C + Σ aⱼ ln πⱼ`,
//! solved by the normal equations. A fit is honest about itself:
//! - it carries an **estimated-color** certificate ([`fs_evidence::Color`]) —
//!   a mined law is a conjecture, never a certified bound;
//! - it exposes its **fit significance** (`r²`), so a corpus with no real
//!   dimensionless structure yields no significant law rather than a
//!   hallucinated one;
//! - it refuses to EXTRAPOLATE: predictions outside the π-space support
//!   (per-coordinate range — exactly the convex hull in one dimension) are
//!   refused, not silently served.
//!
//! Deterministic and side-effect-free (pure arithmetic, no RNG).

pub use fs_evidence::Color;

/// One corpus sample: its dimensionless-group coordinates and the observed QoI.
#[derive(Debug, Clone, PartialEq)]
pub struct Sample {
    /// The dimensionless-group (π) coordinates (all must be positive).
    pub pi: Vec<f64>,
    /// The observed quantity of interest (must be positive).
    pub qoi: f64,
}

impl Sample {
    /// A sample.
    #[must_use]
    pub fn new(pi: Vec<f64>, qoi: f64) -> Sample {
        Sample { pi, qoi }
    }
}

/// A structured mining failure (a refusal that teaches).
#[derive(Debug, Clone, PartialEq)]
pub enum MineError {
    /// Fewer samples than the fit needs (`need = groups + 2`).
    TooFewSamples {
        /// Samples supplied.
        have: usize,
        /// Samples required.
        need: usize,
    },
    /// Samples disagree on the number of π-groups.
    DimMismatch {
        /// Expected group count.
        expected: usize,
        /// A sample's group count.
        found: usize,
    },
    /// A π-coordinate or QoI is not strictly positive (log undefined).
    NonPositive {
        /// A short description of what was non-positive.
        what: &'static str,
        /// The offending value.
        value: f64,
    },
    /// The design is rank-deficient (collinear π-groups) — no unique fit.
    Singular,
    /// A prediction point falls OUTSIDE the trained π-space support.
    Extrapolation {
        /// The out-of-range coordinate index.
        coord: usize,
        /// The requested value.
        value: f64,
        /// The trained range `[min, max]`.
        range: (f64, f64),
    },
}

/// A mined scaling law `y = C · Π πⱼ^{aⱼ}` with its fit quality, validity
/// envelope, and estimated-color certificate.
#[derive(Debug, Clone, PartialEq)]
pub struct MinedLaw {
    /// The multiplicative coefficient `C`.
    pub coefficient: f64,
    /// The exponents `aⱼ`, one per π-group.
    pub exponents: Vec<f64>,
    /// Coefficient of determination `r²` in log space (fit significance).
    pub r_squared: f64,
    /// The trained π-space support: `(min, max)` per group (the validity
    /// envelope — in 1D this IS the convex hull).
    pub envelope: Vec<(f64, f64)>,
    /// How many samples the law was fit on.
    pub samples: usize,
    /// The epistemic color — always estimated (a mined law is a conjecture).
    pub color: Color,
}

impl MinedLaw {
    /// Is the fit significant at `r2_threshold`? A law mined from noise has a
    /// low `r²` and is not significant.
    #[must_use]
    pub fn is_significant(&self, r2_threshold: f64) -> bool {
        self.r_squared >= r2_threshold
    }

    /// Predict `y` at a π-point, refusing to extrapolate beyond the trained
    /// support.
    ///
    /// # Errors
    /// [`MineError::DimMismatch`] on the wrong group count;
    /// [`MineError::NonPositive`] on a non-positive coordinate;
    /// [`MineError::Extrapolation`] outside the validity envelope.
    pub fn predict(&self, pi: &[f64]) -> Result<f64, MineError> {
        if pi.len() != self.exponents.len() {
            return Err(MineError::DimMismatch {
                expected: self.exponents.len(),
                found: pi.len(),
            });
        }
        let mut y = self.coefficient;
        for (j, &p) in pi.iter().enumerate() {
            if !(p.is_finite() && p > 0.0) {
                return Err(MineError::NonPositive {
                    what: "pi coordinate",
                    value: p,
                });
            }
            let (lo, hi) = self.envelope[j];
            if p < lo || p > hi {
                return Err(MineError::Extrapolation {
                    coord: j,
                    value: p,
                    range: (lo, hi),
                });
            }
            y *= p.powf(self.exponents[j]);
        }
        Ok(y)
    }
}

/// Fit a power-law scaling law over a corpus by log-linear least squares.
///
/// # Errors
/// See [`MineError`] — too few samples, dimension mismatch, a non-positive
/// value, or a rank-deficient (collinear) design.
pub fn fit_power_law(corpus: &[Sample]) -> Result<MinedLaw, MineError> {
    let m = corpus.first().map_or(0, |s| s.pi.len());
    let need = m + 2;
    if corpus.len() < need {
        return Err(MineError::TooFewSamples {
            have: corpus.len(),
            need,
        });
    }
    // Build the log-space design matrix rows [1, ln π₁, …, ln π_m] and targets
    // ln(qoi); validate positivity and dimension.
    let mut rows: Vec<Vec<f64>> = Vec::with_capacity(corpus.len());
    let mut targets: Vec<f64> = Vec::with_capacity(corpus.len());
    let mut envelope = vec![(f64::INFINITY, f64::NEG_INFINITY); m];
    for s in corpus {
        if s.pi.len() != m {
            return Err(MineError::DimMismatch {
                expected: m,
                found: s.pi.len(),
            });
        }
        if !(s.qoi.is_finite() && s.qoi > 0.0) {
            return Err(MineError::NonPositive {
                what: "qoi",
                value: s.qoi,
            });
        }
        let mut row = Vec::with_capacity(m + 1);
        row.push(1.0);
        for (j, &p) in s.pi.iter().enumerate() {
            if !(p.is_finite() && p > 0.0) {
                return Err(MineError::NonPositive {
                    what: "pi coordinate",
                    value: p,
                });
            }
            row.push(p.ln());
            envelope[j].0 = envelope[j].0.min(p);
            envelope[j].1 = envelope[j].1.max(p);
        }
        rows.push(row);
        targets.push(s.qoi.ln());
    }
    // Normal equations: (XᵀX) β = Xᵀy, an (m+1)×(m+1) system.
    let k = m + 1;
    let mut ata = vec![vec![0.0_f64; k]; k];
    let mut atb = vec![0.0_f64; k];
    for (row, &t) in rows.iter().zip(&targets) {
        for a in 0..k {
            atb[a] += row[a] * t;
            for b in 0..k {
                ata[a][b] += row[a] * row[b];
            }
        }
    }
    let beta = solve(ata, atb).ok_or(MineError::Singular)?;
    let coefficient = beta[0].exp();
    let exponents = beta[1..].to_vec();

    // Fit quality (r²) in log space + residual dispersion.
    let mean_t = targets.iter().sum::<f64>() / targets.len() as f64;
    let mut ss_res = 0.0;
    let mut ss_tot = 0.0;
    for (row, &t) in rows.iter().zip(&targets) {
        let pred: f64 = row.iter().zip(&beta).map(|(x, b)| x * b).sum();
        ss_res += (t - pred).powi(2);
        ss_tot += (t - mean_t).powi(2);
    }
    let r_squared = if ss_tot <= f64::EPSILON {
        // constant target: a perfect constant fit, else undefined — report 1
        // only when the residual is essentially zero.
        f64::from(u8::from(ss_res <= f64::EPSILON))
    } else {
        1.0 - ss_res / ss_tot
    };
    let dispersion = (ss_res / corpus.len() as f64).sqrt();

    Ok(MinedLaw {
        coefficient,
        exponents,
        r_squared,
        envelope,
        samples: corpus.len(),
        color: Color::Estimated {
            estimator: "power-law-mining".to_string(),
            dispersion,
        },
    })
}

/// Solve `A x = b` for a small dense system by Gaussian elimination with
/// partial pivoting. Returns `None` if the matrix is (near-)singular.
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        // partial pivot: the largest-magnitude entry in this column.
        let piv = (col..n).max_by(|&r1, &r2| {
            a[r1][col]
                .abs()
                .partial_cmp(&a[r2][col].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        if a[piv][col].abs() <= 1e-12 {
            return None; // singular / rank-deficient
        }
        a.swap(col, piv);
        b.swap(col, piv);
        // eliminate below, using a copy of the pivot row to avoid aliasing.
        let pivot = a[col].clone();
        let bcol = b[col];
        for r in (col + 1)..n {
            let f = a[r][col] / pivot[col];
            for (arc, pc) in a[r].iter_mut().zip(&pivot).skip(col) {
                *arc -= f * pc;
            }
            b[r] -= f * bcol;
        }
    }
    // back-substitution.
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let s: f64 = a[i][(i + 1)..n]
            .iter()
            .zip(&x[(i + 1)..n])
            .map(|(aij, xj)| aij * xj)
            .sum();
        x[i] = (b[i] - s) / a[i][i];
    }
    Some(x)
}
