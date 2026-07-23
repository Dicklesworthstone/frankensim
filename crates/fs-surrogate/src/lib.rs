//! fs-surrogate — learned accelerators with GUARANTEES. Layer: L4.
//!
//! Surrogates are 100–10,000× cheap physics — but ML PROPOSES, certified
//! numerics DISPOSES: a surrogate is permitted only inside a certified validity
//! band, and the system never silently trusts a hallucinated field. This v0 is
//! the self-contained core of that discipline:
//!
//! - [`Pod`] — a POD (proper orthogonal decomposition) reduced-order model via
//!   the method of snapshots: an orthonormal reduced basis capturing a target
//!   energy fraction, with a reduced-vs-full reconstruction error;
//! - [`conformal_band`] — a DISTRIBUTION-FREE split-conformal prediction band
//!   from calibration residuals, with empirical `(1−α)` coverage;
//! - [`certify_or_escalate`] — the mechanical policy: inside the validity domain
//!   AND the band narrow enough for the decision → USE the surrogate; otherwise
//!   ESCALATE to a certified solve.
//!
//! The default core is deterministic and dependency-free (including an
//! in-house symmetric eigensolver). Optional features integrate the bounded
//! abstraction ladder and graph-aware escalation routing without changing the
//! default API.

#[cfg(feature = "abstraction-ladder")]
pub mod ladder;

#[cfg(feature = "graph-escalation")]
pub mod escalation;

/// A structured failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SurrogateError {
    /// No snapshots supplied.
    NoSnapshots,
    /// Snapshots disagree on dimension.
    DimMismatch {
        /// Expected dimension.
        expected: usize,
        /// A snapshot's dimension.
        found: usize,
    },
    /// The energy threshold is not in `(0, 1]`.
    BadThreshold,
}

/// A POD reduced-order model.
#[derive(Debug, Clone, PartialEq)]
pub struct Pod {
    mean: Vec<f64>,
    modes: Vec<Vec<f64>>,
    singular_values: Vec<f64>,
    total_energy: f64,
}

impl Pod {
    /// The reduced rank (number of retained modes).
    #[must_use]
    pub fn rank(&self) -> usize {
        self.modes.len()
    }

    /// The fraction of snapshot energy the retained modes capture.
    #[must_use]
    pub fn energy_captured(&self) -> f64 {
        if self.total_energy <= 0.0 {
            return 1.0;
        }
        self.singular_values.iter().map(|s| s * s).sum::<f64>() / self.total_energy
    }

    /// Project a full-space vector onto the reduced coordinates.
    #[must_use]
    pub fn project(&self, x: &[f64]) -> Vec<f64> {
        self.modes
            .iter()
            .map(|m| {
                m.iter()
                    .zip(x)
                    .zip(&self.mean)
                    .map(|((mi, xi), mu)| mi * (xi - mu))
                    .sum()
            })
            .collect()
    }

    /// Reconstruct a full-space vector from reduced coordinates.
    #[must_use]
    pub fn reconstruct(&self, coords: &[f64]) -> Vec<f64> {
        let mut out = self.mean.clone();
        for (c, mode) in coords.iter().zip(&self.modes) {
            for (o, mi) in out.iter_mut().zip(mode) {
                *o += c * mi;
            }
        }
        out
    }

    /// The reduced-vs-full reconstruction error `||x − reconstruct(project(x))||`.
    #[must_use]
    pub fn reconstruction_error(&self, x: &[f64]) -> f64 {
        let r = self.reconstruct(&self.project(x));
        x.iter()
            .zip(&r)
            .map(|(xi, ri)| (xi - ri) * (xi - ri))
            .sum::<f64>()
            .sqrt()
    }
}

/// Build a POD reduced-order model from `snapshots`, retaining the fewest modes
/// that capture at least `energy_threshold` of the (mean-centered) snapshot
/// energy (method of snapshots).
///
/// # Errors
/// [`SurrogateError`] on empty / ragged snapshots or a bad threshold.
pub fn pod(snapshots: &[Vec<f64>], energy_threshold: f64) -> Result<Pod, SurrogateError> {
    let m = snapshots.len();
    if m == 0 {
        return Err(SurrogateError::NoSnapshots);
    }
    if !(energy_threshold > 0.0 && energy_threshold <= 1.0) {
        return Err(SurrogateError::BadThreshold);
    }
    let n = snapshots[0].len();
    for s in snapshots {
        if s.len() != n {
            return Err(SurrogateError::DimMismatch {
                expected: n,
                found: s.len(),
            });
        }
    }
    let mut mean = vec![0.0; n];
    for s in snapshots {
        for (mu, si) in mean.iter_mut().zip(s) {
            *mu += si / m as f64;
        }
    }
    // correlation matrix C = SᵀS  (m×m), S mean-centered.
    let mut c = vec![vec![0.0; m]; m];
    for i in 0..m {
        for j in i..m {
            let dot: f64 = (0..n)
                .map(|k| (snapshots[i][k] - mean[k]) * (snapshots[j][k] - mean[k]))
                .sum();
            c[i][j] = dot;
            c[j][i] = dot;
        }
    }
    let (lambda, vecs) = jacobi_eig(c);
    let total: f64 = lambda.iter().map(|l| l.max(0.0)).sum();

    // choose the rank capturing the target energy.
    let mut cum = 0.0;
    let mut rank = 0;
    for &l in &lambda {
        cum += l.max(0.0);
        rank += 1;
        if total <= 0.0 || cum / total >= energy_threshold {
            break;
        }
    }
    // POD modes φₖ = (S vₖ)/σₖ.
    let mut modes = Vec::new();
    let mut singular_values = Vec::new();
    for k in 0..rank {
        let lk = lambda[k];
        if lk <= 1e-12 {
            break;
        }
        let sk = lk.sqrt();
        let vk = &vecs[k];
        let mode: Vec<f64> = (0..n)
            .map(|row| {
                (0..m)
                    .map(|i| (snapshots[i][row] - mean[row]) * vk[i])
                    .sum::<f64>()
                    / sk
            })
            .collect();
        modes.push(mode);
        singular_values.push(sk);
    }
    Ok(Pod {
        mean,
        modes,
        singular_values,
        total_energy: total,
    })
}

/// A distribution-free split-conformal prediction band.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConformalBand {
    /// The band half-width.
    pub half_width: f64,
    /// The miscoverage level `α` (target coverage `1 − α`).
    pub alpha: f64,
}

impl ConformalBand {
    /// Does the band cover the truth around a prediction?
    #[must_use]
    pub fn covers(&self, prediction: f64, truth: f64) -> bool {
        (prediction - truth).abs() <= self.half_width
    }
}

/// Build a split-conformal prediction band from calibration residuals at
/// miscoverage `alpha`: the `⌈(1−α)(n+1)⌉`-th smallest absolute residual, which
/// gives distribution-free `(1−α)` marginal coverage.
///
/// Distribution-free `(1−α)` coverage requires that order statistic to EXIST,
/// i.e. `⌈(1−α)(n+1)⌉ ≤ n`, equivalently `α ≥ 1/(n+1)`. With fewer calibration
/// points than that, the only honest `(1−α)` band is unbounded: the function
/// returns `half_width = +∞`, which forces [`certify_or_escalate`] to escalate.
/// (A previous `.clamp(1, n)` silently returned the MAX residual there, whose
/// true coverage is only `n/(n+1) < 1−α` — a band that under-covers while
/// claiming its nominal guarantee.)
///
/// # Panics
/// If `residuals` is empty or `alpha ∉ (0, 1)`.
#[must_use]
pub fn conformal_band(residuals: &[f64], alpha: f64) -> ConformalBand {
    assert!(
        !residuals.is_empty(),
        "need at least one calibration residual"
    );
    assert!(alpha > 0.0 && alpha < 1.0, "alpha must be in (0, 1)");
    let mut r: Vec<f64> = residuals.iter().map(|x| x.abs()).collect();
    r.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = r.len();
    let rank = ((1.0 - alpha) * (n as f64 + 1.0)).ceil() as usize;
    // `(1−α)(n+1) > 0` so `rank ≥ 1`; only the upper end can fall off the sample.
    let half_width = if rank > n { f64::INFINITY } else { r[rank - 1] };
    ConformalBand { half_width, alpha }
}

/// Empirical coverage of a band on held-out `(prediction, truth)` pairs.
#[must_use]
pub fn empirical_coverage(band: &ConformalBand, pairs: &[(f64, f64)]) -> f64 {
    if pairs.is_empty() {
        return 1.0;
    }
    let hit = pairs.iter().filter(|(p, t)| band.covers(*p, *t)).count();
    hit as f64 / pairs.len() as f64
}

/// The certify-or-escalate verdict.
#[derive(Debug, Clone, PartialEq)]
pub enum Decision {
    /// Use the surrogate — its band is trustworthy and decision-relevant.
    UseSurrogate {
        /// The band half-width backing the decision.
        band_half_width: f64,
    },
    /// Escalate to a certified solve.
    Escalate {
        /// Why.
        reason: String,
    },
}

/// The mechanical certify-or-escalate policy: use the surrogate iff every
/// numeric policy input is valid, the query is inside the validity domain, and
/// the conformal band is at least as tight as the decision tolerance; otherwise
/// escalate.
///
/// [`ConformalBand`] has public fields, and an honestly unbounded conformal band
/// uses `half_width = +∞`. Revalidate the complete policy at this decision
/// boundary so caller-constructed malformed or unbounded values cannot
/// authorize a surrogate through IEEE-754 comparisons such as `+∞ <= +∞`.
#[must_use]
pub fn certify_or_escalate(
    band: &ConformalBand,
    in_validity_domain: bool,
    decision_tolerance: f64,
) -> Decision {
    if band.half_width == f64::INFINITY {
        return Decision::Escalate {
            reason: "conformal band is unbounded".to_string(),
        };
    }
    if !band.half_width.is_finite() || band.half_width < 0.0 {
        return Decision::Escalate {
            reason: "conformal band half-width must be finite and non-negative".to_string(),
        };
    }
    if !(band.alpha.is_finite() && band.alpha > 0.0 && band.alpha < 1.0) {
        return Decision::Escalate {
            reason: "conformal band alpha must be finite and in (0, 1)".to_string(),
        };
    }
    if !decision_tolerance.is_finite() || decision_tolerance < 0.0 {
        return Decision::Escalate {
            reason: "decision tolerance must be finite and non-negative".to_string(),
        };
    }
    if !in_validity_domain {
        return Decision::Escalate {
            reason: "query outside the surrogate's validity domain".to_string(),
        };
    }
    if band.half_width <= decision_tolerance {
        Decision::UseSurrogate {
            band_half_width: band.half_width,
        }
    } else {
        Decision::Escalate {
            reason: format!(
                "band half-width {:.3e} exceeds decision tolerance {:.3e}",
                band.half_width, decision_tolerance
            ),
        }
    }
}

/// Jacobi symmetric eigensolver returning `(eigenvalues desc, eigenvectors)`
/// where `eigenvectors[k]` is the unit eigenvector for `eigenvalues[k]`.
// A dense symmetric eigen-kernel: `a[i][j]` / `v[i][j]` are inherently
// 2D-indexed by row/column, so index loops are the correct, readable form.
#[allow(clippy::needless_range_loop)]
fn jacobi_eig(mut a: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut v = vec![vec![0.0; n]; n];
    for i in 0..n {
        v[i][i] = 1.0;
    }
    for _ in 0..100 {
        let mut off = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                off += a[i][j] * a[i][j];
            }
        }
        if off <= 1e-24 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p][q].abs() <= 1e-18 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let (vkp, vkq) = (v[k][p], v[k][q]);
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&i, &j| {
        a[j][j]
            .partial_cmp(&a[i][i])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let eigenvalues: Vec<f64> = idx.iter().map(|&i| a[i][i]).collect();
    let eigenvectors: Vec<Vec<f64>> = idx
        .iter()
        .map(|&col| (0..n).map(|row| v[row][col]).collect())
        .collect();
    (eigenvalues, eigenvectors)
}
