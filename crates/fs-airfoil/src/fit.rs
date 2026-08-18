//! Regime-partitioned tensor-product cubic B-spline residual machinery
//! (bead wf-root-guzez.5.1.2, E4.0b). Plan §5.2.1: residual surfaces over
//! (α, log Re, δ) ON TOP of the analytic baselines, with COEFFICIENT-
//! DIFFERENCE shape constraints declared per regime — global monotonicity
//! is never imposed where the physics is not monotone. Constraints are
//! enforced fail-closed: a fit that violates its declared constraints is
//! a typed refusal, never a silently projected surface.
//!
//! The spline basis is clamped uniform cubic (Cox–de Boor); a degenerate
//! axis (`n_coef == 1`) means "no dependence on this axis" (e.g. the δ
//! axis of a wing table). Fitting is penalized least squares via normal
//! equations and a dense Cholesky solve — table sizes are bounded by
//! admission, so the dense path is exact and deterministic.

use crate::Refusal;

/// Maximum coefficients per axis (admission cap; refuse above).
pub const MAX_COEF_PER_AXIS: usize = 64;
/// Minimum coefficients for a non-degenerate cubic axis.
pub const MIN_CUBIC_COEF: usize = 4;
/// Ridge regularization added to the normal equations (deterministic).
pub const FIT_RIDGE: f64 = 1.0e-10;
/// Tolerance on coefficient-difference constraints (float slack only).
pub const CONSTRAINT_TOL: f64 = 1.0e-12;

/// One spline axis: clamped uniform cubic knots on [lo, hi].
#[derive(Clone, Debug, PartialEq)]
pub struct BsplineAxis {
    /// Axis name ("alpha_rad", "log10_re", "delta_rad").
    pub name: &'static str,
    /// Domain lower bound.
    pub lo: f64,
    /// Domain upper bound.
    pub hi: f64,
    /// Coefficient count: 1 (degenerate, constant) or ≥ 4 (cubic).
    pub n_coef: usize,
}

impl BsplineAxis {
    /// Validate the axis against the admission caps.
    ///
    /// # Errors
    /// `axis-domain-invalid`, `axis-coef-count-invalid`.
    pub fn admit(&self) -> Result<(), Refusal> {
        if !(self.lo.is_finite() && self.hi.is_finite()) || self.lo >= self.hi {
            return Err(Refusal {
                code: "axis-domain-invalid",
                message: format!(
                    "axis {}: [{}, {}] is not a finite interval",
                    self.name, self.lo, self.hi
                ),
                ranked_repairs: vec!["declare lo < hi, both finite".into()],
            });
        }
        let n = self.n_coef;
        if n != 1 && !(MIN_CUBIC_COEF..=MAX_COEF_PER_AXIS).contains(&n) {
            return Err(Refusal {
                code: "axis-coef-count-invalid",
                message: format!(
                    "axis {}: n_coef {n} must be 1 (degenerate) or in [{MIN_CUBIC_COEF}, {MAX_COEF_PER_AXIS}]",
                    self.name
                ),
                ranked_repairs: vec!["use n_coef = 1 for a constant axis, else at least 4".into()],
            });
        }
        Ok(())
    }

    /// Clamped uniform knot value t[i] for i in 0..n+4.
    fn knot(&self, i: usize) -> f64 {
        let n = self.n_coef;
        let interior = n - 3; // number of spans
        if i < 4 {
            self.lo
        } else if i >= n {
            self.hi
        } else {
            let step = (self.hi - self.lo) / interior as f64;
            self.lo + (i - 3) as f64 * step
        }
    }

    /// The 4 nonzero cubic basis values at `x` and the first coefficient
    /// index they multiply. Degenerate axes return ([1,0,0,0], 0).
    fn basis(&self, x: f64) -> (usize, [f64; 4]) {
        if self.n_coef == 1 {
            return (0, [1.0, 0.0, 0.0, 0.0]);
        }
        let n = self.n_coef;
        let x = x.clamp(self.lo, self.hi);
        // Span k: t[k] <= x < t[k+1], k in [3, n-1]; x == hi → k = n-1.
        let interior = n - 3;
        let step = (self.hi - self.lo) / interior as f64;
        let raw = ((x - self.lo) / step) as usize;
        let k = (raw.min(interior - 1)) + 3;
        // Cox–de Boor (NURBS-book basis-funs), degree 3.
        let mut nvals = [1.0f64, 0.0, 0.0, 0.0];
        let mut left = [0.0f64; 4];
        let mut right = [0.0f64; 4];
        for j in 1..=3 {
            left[j] = x - self.knot(k + 1 - j);
            right[j] = self.knot(k + j) - x;
            let mut saved = 0.0;
            for r in 0..j {
                let denom = right[r + 1] + left[j - r];
                let temp = nvals[r] / denom;
                nvals[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            nvals[j] = saved;
        }
        (k - 3, nvals)
    }

    /// Greville abscissa of coefficient `i` (the x where a linear function's
    /// exact-reproduction coefficient lives).
    #[must_use]
    pub fn greville(&self, i: usize) -> f64 {
        if self.n_coef == 1 {
            return f64::midpoint(self.lo, self.hi);
        }
        (self.knot(i + 1) + self.knot(i + 2) + self.knot(i + 3)) / 3.0
    }
}

/// Shape-constraint direction along one axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DiffDirection {
    /// Adjacent coefficient differences must be ≥ −tol.
    NonDecreasing,
    /// Adjacent coefficient differences must be ≤ +tol.
    NonIncreasing,
}

/// A coefficient-difference constraint along one axis of one regime patch.
/// Local by construction: it binds THIS patch only (plan law — no global
/// monotonicity where the physics is not monotone).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DiffConstraint {
    /// Axis index (0 = α, 1 = log Re, 2 = δ).
    pub axis: usize,
    /// Required direction of coefficient differences along that axis.
    pub direction: DiffDirection,
}

/// A fitted residual surface: tensor-product cubic B-spline over three
/// axes with declared shape constraints.
#[derive(Clone, Debug, PartialEq)]
pub struct ResidualSurface {
    /// The three axes (α, log Re, δ); degenerate axes have `n_coef = 1`.
    pub axes: [BsplineAxis; 3],
    /// Coefficients, index = (i·n1 + j)·n2 + k.
    pub coef: Vec<f64>,
    /// Declared shape constraints (verified, fail-closed).
    pub constraints: Vec<DiffConstraint>,
}

/// One training sample for a residual fit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FitSample {
    /// (α, log Re, δ) query point.
    pub x: [f64; 3],
    /// Residual value (measured − baseline) at the point.
    pub y: f64,
}

fn cholesky_solve(a: &mut [f64], b: &mut [f64], n: usize) -> Result<(), Refusal> {
    // In-place Cholesky A = L·Lᵀ (lower), then two triangular solves.
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i * n + j];
            for k in 0..j {
                sum -= a[i * n + k] * a[j * n + k];
            }
            if i == j {
                if sum <= 0.0 {
                    return Err(Refusal {
                        code: "fit-normal-equations-singular",
                        message: format!(
                            "Cholesky pivot {sum:e} at row {i} — samples do not determine the surface"
                        ),
                        ranked_repairs: vec![
                            "supply samples covering every spline span".into(),
                            "reduce n_coef on sparsely sampled axes".into(),
                        ],
                    });
                }
                a[i * n + j] = sum.sqrt();
            } else {
                a[i * n + j] = sum / a[j * n + j];
            }
        }
    }
    for i in 0..n {
        let mut sum = b[i];
        for k in 0..i {
            sum -= a[i * n + k] * b[k];
        }
        b[i] = sum / a[i * n + i];
    }
    for i in (0..n).rev() {
        let mut sum = b[i];
        for k in i + 1..n {
            sum -= a[k * n + i] * b[k];
        }
        b[i] = sum / a[i * n + i];
    }
    Ok(())
}

impl ResidualSurface {
    /// Total coefficient count.
    #[must_use]
    pub fn n_total(&self) -> usize {
        self.axes[0].n_coef * self.axes[1].n_coef * self.axes[2].n_coef
    }

    /// Evaluate the surface at (α, log Re, δ). Queries clamp to the patch
    /// domain — DOMAIN admission (refusing out-of-domain) is the caller's
    /// contract via `admit_query`; the patch itself is total on its box.
    #[must_use]
    pub fn eval(&self, x: [f64; 3]) -> f64 {
        let (i0, b0) = self.axes[0].basis(x[0]);
        let (i1, b1) = self.axes[1].basis(x[1]);
        let (i2, b2) = self.axes[2].basis(x[2]);
        let (n1, n2) = (self.axes[1].n_coef, self.axes[2].n_coef);
        let width = |axis: &BsplineAxis| if axis.n_coef == 1 { 1 } else { 4 };
        let mut acc = 0.0;
        for (a, &w_a) in b0.iter().enumerate().take(width(&self.axes[0])) {
            for (b, &w_b) in b1.iter().enumerate().take(width(&self.axes[1])) {
                for (c, &w_c) in b2.iter().enumerate().take(width(&self.axes[2])) {
                    let idx = ((i0 + a) * n1 + (i1 + b)) * n2 + (i2 + c);
                    acc += w_a * w_b * w_c * self.coef[idx];
                }
            }
        }
        acc
    }

    /// Verify every declared coefficient-difference constraint.
    ///
    /// # Errors
    /// `fit-constraint-violated` naming the axis, direction, and location.
    pub fn verify_constraints(&self) -> Result<(), Refusal> {
        let (n0, n1, n2) = (
            self.axes[0].n_coef,
            self.axes[1].n_coef,
            self.axes[2].n_coef,
        );
        let idx = |i: usize, j: usize, k: usize| (i * n1 + j) * n2 + k;
        for con in &self.constraints {
            let n_axis = self.axes[con.axis].n_coef;
            for i in 0..n0 {
                for j in 0..n1 {
                    for k in 0..n2 {
                        let along = [i, j, k][con.axis];
                        if along + 1 >= n_axis {
                            continue;
                        }
                        let mut next = [i, j, k];
                        next[con.axis] += 1;
                        let d = self.coef[idx(next[0], next[1], next[2])] - self.coef[idx(i, j, k)];
                        let bad = match con.direction {
                            DiffDirection::NonDecreasing => d < -CONSTRAINT_TOL,
                            DiffDirection::NonIncreasing => d > CONSTRAINT_TOL,
                        };
                        if bad {
                            return Err(Refusal {
                                code: "fit-constraint-violated",
                                message: format!(
                                    "coefficient difference {d:e} along axis {} ({:?}) at ({i},{j},{k}) violates {:?}",
                                    con.axis, self.axes[con.axis].name, con.direction
                                ),
                                ranked_repairs: vec![
                                    "the data contradicts the declared shape constraint — re-examine the regime partition".into(),
                                    "drop the constraint ONLY with a physics justification recorded in the table provenance".into(),
                                ],
                            });
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// Penalized least-squares fit of a residual surface, then fail-closed
    /// constraint verification.
    ///
    /// # Errors
    /// Axis admission refusals; `insufficient-samples` (need ≥ n_total,
    /// tested at cap and cap−1); `non-finite-input` on any sample;
    /// `fit-normal-equations-singular`; `fit-constraint-violated`.
    pub fn fit(
        axes: [BsplineAxis; 3],
        constraints: Vec<DiffConstraint>,
        samples: &[FitSample],
    ) -> Result<Self, Refusal> {
        for axis in &axes {
            axis.admit()?;
        }
        for con in &constraints {
            if con.axis >= 3 {
                return Err(Refusal {
                    code: "constraint-axis-invalid",
                    message: format!("constraint axis {} out of range", con.axis),
                    ranked_repairs: vec!["axes are 0 = alpha, 1 = log Re, 2 = delta".into()],
                });
            }
        }
        let mut surface = ResidualSurface {
            axes,
            coef: Vec::new(),
            constraints,
        };
        let n = surface.n_total();
        if samples.len() < n {
            return Err(Refusal {
                code: "insufficient-samples",
                message: format!(
                    "{} samples cannot determine {n} coefficients",
                    samples.len()
                ),
                ranked_repairs: vec!["supply at least one sample per coefficient".into()],
            });
        }
        let mut ata = vec![0.0f64; n * n];
        let mut atb = vec![0.0f64; n];
        let (n1, n2) = (surface.axes[1].n_coef, surface.axes[2].n_coef);
        let width = |axis: &BsplineAxis| if axis.n_coef == 1 { 1 } else { 4 };
        for s in samples {
            if !(s.x.iter().all(|v| v.is_finite()) && s.y.is_finite()) {
                return Err(Refusal {
                    code: "non-finite-input",
                    message: format!("fit sample {s:?} contains a non-finite value"),
                    ranked_repairs: vec!["filter or repair the ingest upstream".into()],
                });
            }
            let (i0, b0) = surface.axes[0].basis(s.x[0]);
            let (i1, b1) = surface.axes[1].basis(s.x[1]);
            let (i2, b2) = surface.axes[2].basis(s.x[2]);
            // Nonzero pattern of this row: ≤ 64 entries.
            let mut cols = [(0usize, 0.0f64); 64];
            let mut m = 0;
            for (a, &w_a) in b0.iter().enumerate().take(width(&surface.axes[0])) {
                for (b, &w_b) in b1.iter().enumerate().take(width(&surface.axes[1])) {
                    for (c, &w_c) in b2.iter().enumerate().take(width(&surface.axes[2])) {
                        let idx = ((i0 + a) * n1 + (i1 + b)) * n2 + (i2 + c);
                        cols[m] = (idx, w_a * w_b * w_c);
                        m += 1;
                    }
                }
            }
            for &(ci, vi) in cols.iter().take(m) {
                atb[ci] += vi * s.y;
                for &(cj, vj) in cols.iter().take(m) {
                    ata[ci * n + cj] += vi * vj;
                }
            }
        }
        for i in 0..n {
            ata[i * n + i] += FIT_RIDGE;
        }
        cholesky_solve(&mut ata, &mut atb, n)?;
        surface.coef = atb;
        surface.verify_constraints()?;
        Ok(surface)
    }
}

/// Verify C⁰ continuity of two regime patches across their shared α
/// boundary (patch `a` ends where patch `b` begins): the shared face is
/// sampled on a grid and the maximum mismatch must stay within `tol`.
///
/// # Errors
/// `regime-boundary-mismatch` (patches do not abut),
/// `regime-boundary-discontinuity` (face mismatch beyond tol, location
/// and magnitude stated).
pub fn verify_regime_continuity(
    a: &ResidualSurface,
    b: &ResidualSurface,
    tol: f64,
) -> Result<f64, Refusal> {
    let boundary = a.axes[0].hi;
    if (boundary - b.axes[0].lo).abs() > 1.0e-12 {
        return Err(Refusal {
            code: "regime-boundary-mismatch",
            message: format!(
                "patch A ends at α = {boundary} but patch B begins at {}",
                b.axes[0].lo
            ),
            ranked_repairs: vec!["regime partitions must tile α without gaps or overlaps".into()],
        });
    }
    let grid = 9;
    let mut worst = 0.0f64;
    for j in 0..grid {
        for k in 0..grid {
            let fj = f64::from(j) / f64::from(grid - 1);
            let fk = f64::from(k) / f64::from(grid - 1);
            // Sample the shared face over the INTERSECTION of the secondary
            // domains (patches may declare different Re/δ boxes).
            let re = a.axes[1].lo.max(b.axes[1].lo)
                + fj * (a.axes[1].hi.min(b.axes[1].hi) - a.axes[1].lo.max(b.axes[1].lo));
            let de = a.axes[2].lo.max(b.axes[2].lo)
                + fk * (a.axes[2].hi.min(b.axes[2].hi) - a.axes[2].lo.max(b.axes[2].lo));
            let diff = (a.eval([boundary, re, de]) - b.eval([boundary, re, de])).abs();
            worst = worst.max(diff);
        }
    }
    if worst > tol {
        return Err(Refusal {
            code: "regime-boundary-discontinuity",
            message: format!(
                "shared-face mismatch {worst:e} exceeds tol {tol:e} at α = {boundary}"
            ),
            ranked_repairs: vec![
                "refit the patches with shared boundary samples".into(),
                "a real discontinuity needs an explicit transition regime, not a wider tol".into(),
            ],
        });
    }
    Ok(worst)
}
