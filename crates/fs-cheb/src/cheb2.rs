//! 2D low-rank function objects (bead kw89): Chebfun2-style adaptive
//! cross approximation. f(x, y) ≈ Σ_k u_k(x)·v_k(y)/p_k where each
//! cross comes from the residual's maximal-|value| pivot on a
//! deterministic sample grid, and the slice functions are Cheb1
//! interpolants. Separable functions are captured EXACTLY at their
//! rank (gated); smooth non-separable ones converge spectrally in
//! rank. Integration is the separable product of 1D integrals.

use crate::{Cheb1, affine_from_reference, stable_finite_product3, stable_finite_sum};
use fs_math::det;

/// Multiply one low-rank component without committing to an intermediate that
/// overflows or underflows even though the final three-factor product is
/// representable. The established `(u*v)*inv_pivot` order remains first so
/// ordinary-domain bits do not move.
fn stable_component_product(u: f64, v: f64, inverse_pivot: f64) -> f64 {
    stable_finite_product3(u, v, inverse_pivot, "Cheb2 component product")
}

fn stable_component_sum<F>(len: usize, component: F, operation: &str) -> f64
where
    F: Fn(usize) -> f64,
{
    stable_finite_sum(len, component, operation)
}

#[cfg(test)]
mod stable_product_tests {
    use super::{fixed_cheb, stable_component_product, stable_component_sum};

    #[test]
    fn rescaling_does_not_magnify_a_rounded_subnormal_intermediate() {
        let u = 1e-200_f64;
        let v = 3.7e-124_f64;
        let inverse_pivot = 1e300_f64;
        let stable_order = (u * inverse_pivot) * v;
        let underflowed_order = (u * v) * inverse_pivot;
        assert_ne!(stable_order.to_bits(), underflowed_order.to_bits());
        assert_eq!(
            stable_component_product(u, v, inverse_pivot).to_bits(),
            stable_order.to_bits()
        );
    }

    #[test]
    fn representable_sum_survives_an_overflowing_prefix() {
        let half_max = f64::MAX / 2.0;
        let values = [half_max, half_max, half_max, -half_max];
        let result = stable_component_sum(values.len(), |index| values[index], "test sum");
        assert!(result.is_finite());
        assert!((f64::MAX - result) / f64::MAX <= f64::EPSILON);
    }

    #[test]
    fn representable_residual_survives_finite_prefix_cancellation() {
        let residual = fs_math::det::powi(2.0, -54);
        let values = [1.0, residual, -1.0];
        assert_eq!(1.0 + residual, 1.0, "fixture must exercise absorption");
        assert_eq!(
            stable_component_sum(values.len(), |index| values[index], "test sum").to_bits(),
            residual.to_bits()
        );
    }

    #[test]
    fn fixed_dct_applies_its_scale_before_a_prefix_can_overflow() {
        let value = f64::MAX / 4.0;
        let interpolant = fixed_cheb(&|_| value, -1.0, 1.0, 16);
        assert_eq!(interpolant.coeffs()[0].to_bits(), (2.0 * value).to_bits());
        assert!(
            interpolant
                .coeffs()
                .iter()
                .all(|coefficient| coefficient.is_finite())
        );
    }
}

/// Fixed-resolution Chebyshev interpolant of a slice function: sample
/// at `n` first-kind points, direct DCT-II, light trailing truncation.
/// ACA residual slices carry ABSOLUTE cancellation noise (~1e-16 of
/// the ORIGINAL function's scale), so the adaptive builder's
/// machine-precision plateau test can never pass on late-rank slices
/// — fixed resolution is the standard Chebfun2 practice (measured:
/// the adaptive path panicked "not resolved" on rank-2 residuals).
fn fixed_cheb<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, n: usize) -> Cheb1 {
    assert!(n > 0, "fixed Cheb2 slice degree must be positive");
    let vals: Vec<f64> = (0..n)
        .map(|k| {
            let theta = std::f64::consts::PI * (2.0 * k as f64 + 1.0) / (2.0 * n as f64);
            let t = det::cos(theta);
            let y = f(affine_from_reference(t, a, b));
            assert!(y.is_finite(), "Cheb2 slice samples must be finite");
            y
        })
        .collect();
    let mut coeffs = vec![0.0f64; n];
    let dct_scale = 2.0 / n as f64;
    for (j, cj) in coeffs.iter_mut().enumerate() {
        *cj = stable_finite_sum(
            vals.len(),
            |k| {
                let angle =
                    std::f64::consts::PI * j as f64 * (2.0 * k as f64 + 1.0) / (2.0 * n as f64);
                stable_finite_product3(vals[k], det::cos(angle), dct_scale, "Cheb2 fixed DCT term")
            },
            "Cheb2 fixed DCT accumulation",
        );
    }
    let cmax = coeffs.iter().fold(0.0f64, |m, &c| m.max(c.abs()));
    let keep = coeffs
        .iter()
        .rposition(|&c| c.abs() > 1e-15 * cmax.max(f64::MIN_POSITIVE))
        .map_or(1, |p| p + 1);
    coeffs.truncate(keep);
    Cheb1::from_coeffs(a, b, coeffs)
}

/// A rank-`k` bivariate approximant on [a, b] × [c, d].
pub struct Cheb2 {
    /// x-slices u_k (Cheb1 on [a, b]).
    pub cols: Vec<Cheb1>,
    /// y-slices v_k (Cheb1 on [c, d]).
    pub rows: Vec<Cheb1>,
    /// Pivot scalings 1/e_k(x*, y*).
    pub inv_pivots: Vec<f64>,
    /// Worst residual sample at the stop (the accuracy ledger).
    pub residual: f64,
}

impl Cheb2 {
    fn assert_valid(&self) {
        assert!(
            self.cols.len() == self.rows.len() && self.cols.len() == self.inv_pivots.len(),
            "Cheb2 rank component lengths must match"
        );
        assert!(
            self.inv_pivots.iter().all(|p| p.is_finite() && *p != 0.0),
            "Cheb2 inverse pivots must be finite and non-zero"
        );
        if let Some(domain) = self.cols.first().map(Cheb1::domain) {
            assert!(
                self.cols.iter().all(|column| column.domain() == domain),
                "all Cheb2 columns must share one x-domain"
            );
        }
        if let Some(domain) = self.rows.first().map(Cheb1::domain) {
            assert!(
                self.rows.iter().all(|row| row.domain() == domain),
                "all Cheb2 rows must share one y-domain"
            );
        }
        assert!(
            self.residual.is_finite() && self.residual >= 0.0,
            "Cheb2 residual must be finite and non-negative"
        );
    }

    /// Build by adaptive cross approximation: pivots from an
    /// `ns × ns` deterministic Chebyshev-point sample grid, stop when
    /// the residual pivot falls below `tol`·(first pivot) or rank
    /// `max_rank` is hit.
    ///
    /// # Panics
    /// If a domain endpoint or tolerance is non-finite, either domain is
    /// degenerate, `tol < 0`, `max_rank == 0`, the derived deterministic sample
    /// grid or transform size overflows, a sampled/pivot value is non-finite,
    /// a pivot is zero where inversion is required, or the resulting public
    /// low-rank factors violate their shared-domain/shape invariants. These are
    /// legacy assertion-based preconditions; the budgeted fallible constructor
    /// is tracked separately.
    #[must_use]
    pub fn build<F: Fn(f64, f64) -> f64>(
        f: &F,
        domain: (f64, f64, f64, f64),
        tol: f64,
        max_rank: usize,
        max_degree: usize,
    ) -> Cheb2 {
        let (a, b, c, d) = domain;
        assert!(
            a.is_finite() && b.is_finite() && c.is_finite() && d.is_finite() && b > a && d > c,
            "finite non-degenerate Cheb2 domain required"
        );
        assert!(
            tol.is_finite() && tol >= 0.0,
            "Cheb2 tolerance must be finite and non-negative"
        );
        assert!(max_rank > 0, "Cheb2 max_rank must be positive");
        let degree_cap = max_degree.max(16);
        let ns = degree_cap.saturating_mul(2).max(33);
        let xs: Vec<f64> = (0..=ns)
            .map(|k| {
                let t = det::cos(std::f64::consts::PI * k as f64 / ns as f64);
                affine_from_reference(t, a, b)
            })
            .collect();
        let ys: Vec<f64> = (0..=ns)
            .map(|k| {
                let t = det::cos(std::f64::consts::PI * k as f64 / ns as f64);
                affine_from_reference(t, c, d)
            })
            .collect();
        let mut cols: Vec<Cheb1> = Vec::new();
        let mut rows: Vec<Cheb1> = Vec::new();
        let mut inv_pivots: Vec<f64> = Vec::new();
        let mut first_pivot = 0.0f64;
        let mut residual = 0.0f64;
        for _ in 0..max_rank {
            // Residual on the sample grid; deterministic argmax
            // (first-in-scan-order tie break).
            let approx = |x: f64, y: f64, cols: &[Cheb1], rows: &[Cheb1], ip: &[f64]| -> f64 {
                stable_component_sum(
                    cols.len(),
                    |k| stable_component_product(cols[k].eval(x), rows[k].eval(y), ip[k]),
                    "Cheb2 approximation accumulation",
                )
            };
            let mut best = (0usize, 0usize, 0.0f64);
            for (i, &x) in xs.iter().enumerate() {
                for (j, &y) in ys.iter().enumerate() {
                    let e = f(x, y) - approx(x, y, &cols, &rows, &inv_pivots);
                    assert!(e.is_finite(), "Cheb2 residual sample must be finite");
                    if e.abs() > best.2 {
                        best = (i, j, e.abs());
                    }
                }
            }
            residual = best.2;
            if first_pivot == 0.0 {
                first_pivot = best.2;
            }
            if best.2 <= tol * first_pivot || best.2 == 0.0 {
                break;
            }
            let (xp, yp) = (xs[best.0], ys[best.1]);
            let pivot = f(xp, yp) - approx(xp, yp, &cols, &rows, &inv_pivots);
            assert!(
                pivot.is_finite() && pivot != 0.0,
                "Cheb2 pivot must be finite and non-zero"
            );
            // Residual slices as Cheb1 interpolants.
            let u = fixed_cheb(
                &|x| f(x, yp) - approx(x, yp, &cols, &rows, &inv_pivots),
                a,
                b,
                degree_cap,
            );
            let v = fixed_cheb(
                &|y| f(xp, y) - approx(xp, y, &cols, &rows, &inv_pivots),
                c,
                d,
                degree_cap,
            );
            cols.push(u);
            rows.push(v);
            let inverse_pivot = 1.0 / pivot;
            assert!(
                inverse_pivot.is_finite(),
                "Cheb2 inverse pivot is not representable as finite f64"
            );
            inv_pivots.push(inverse_pivot);
        }
        Cheb2 {
            cols,
            rows,
            inv_pivots,
            residual,
        }
    }

    /// Rank of the approximant.
    #[must_use]
    pub fn rank(&self) -> usize {
        self.assert_valid();
        self.cols.len()
    }

    /// Evaluate at (x, y).
    #[must_use]
    pub fn eval(&self, x: f64, y: f64) -> f64 {
        self.assert_valid();
        assert!(
            x.is_finite() && y.is_finite(),
            "Cheb2 evaluation point must be finite"
        );
        stable_component_sum(
            self.cols.len(),
            |k| {
                stable_component_product(
                    self.cols[k].eval(x),
                    self.rows[k].eval(y),
                    self.inv_pivots[k],
                )
            },
            "Cheb2 evaluation accumulation",
        )
    }

    /// ∫∫ f over the domain: Σ_k (∫u_k)(∫v_k)/p_k — the separable
    /// payoff of the low-rank form.
    #[must_use]
    pub fn integral(&self) -> f64 {
        self.assert_valid();
        stable_component_sum(
            self.cols.len(),
            |k| {
                stable_component_product(
                    self.cols[k].integral(),
                    self.rows[k].integral(),
                    self.inv_pivots[k],
                )
            },
            "Cheb2 integral accumulation",
        )
    }
}
