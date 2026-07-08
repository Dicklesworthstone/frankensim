//! Wavelet-style tile thresholding (mechanism 4 of 4):
//! compression-as-adaptivity. A 2D Haar transform of a cell-average
//! field, with detail coefficients dropped wherever they are smaller
//! than the LOCAL budget — and the budget is DWR-weighted: where the
//! adjoint says the goal cannot see, accuracy is spent for storage.
//! The battery measures both directions: goal impact stays under the
//! budget while unweighted thresholding at the same compression ratio
//! does measurably worse.

/// Threshold outcome: the reconstructed field plus bookkeeping.
#[derive(Debug, Clone)]
pub struct ThresholdOutcome {
    /// Reconstructed (compressed) cell averages, row-major n×n.
    pub field: Vec<f64>,
    /// Detail coefficients kept.
    pub kept: usize,
    /// Detail coefficients total.
    pub total: usize,
}

impl ThresholdOutcome {
    /// Compression ratio (total / kept, counting the scaling coeff).
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn ratio(&self) -> f64 {
        (self.total + 1) as f64 / (self.kept + 1) as f64
    }
}

/// Standard 2D Haar (nonstandard/square decomposition) threshold of a
/// row-major `n×n` field (`n` a power of two). A detail coefficient at
/// scale `s` covering block `(bi, bj)` survives iff `|d|` exceeds the
/// MINIMUM of `budget` over its block — conservative: a coefficient is
/// only dropped where every covered cell tolerates it.
///
/// # Panics
/// If `n` is not a power of two or `field.len() != n*n`.
#[must_use]
pub fn haar_threshold(
    field: &[f64],
    n: usize,
    budget: &dyn Fn(usize, usize) -> f64,
) -> ThresholdOutcome {
    assert!(n.is_power_of_two(), "n must be a power of two");
    assert_eq!(field.len(), n * n, "field must be n×n");
    let mut a = field.to_vec();
    let mut kept = 0usize;
    let mut total = 0usize;
    // Precompute per-cell budgets once.
    let cell_budget: Vec<f64> = (0..n * n).map(|k| budget(k % n, k / n)).collect();
    let mut size = n;
    // Forward transform with in-place thresholding per level.
    while size > 1 {
        let half = size / 2;
        let block = n / half; // cells covered per coefficient axis
        let mut next = a.clone();
        for bj in 0..half {
            for bi in 0..half {
                let v00 = a[2 * bi + 2 * bj * n];
                let v10 = a[2 * bi + 1 + 2 * bj * n];
                let v01 = a[2 * bi + (2 * bj + 1) * n];
                let v11 = a[2 * bi + 1 + (2 * bj + 1) * n];
                let avg = 0.25 * (v00 + v10 + v01 + v11);
                let dx = 0.25 * (v00 - v10 + v01 - v11);
                let dy = 0.25 * (v00 + v10 - v01 - v11);
                let dd = 0.25 * (v00 - v10 - v01 + v11);
                // Minimum budget over the covered block.
                let mut bmin = f64::INFINITY;
                for cj in (bj * block)..((bj + 1) * block).min(n) {
                    for ci in (bi * block)..((bi + 1) * block).min(n) {
                        bmin = bmin.min(cell_budget[ci + cj * n]);
                    }
                }
                let mut keep = |d: f64| -> f64 {
                    total += 1;
                    if d.abs() > bmin {
                        kept += 1;
                        d
                    } else {
                        0.0
                    }
                };
                let (dx, dy, dd) = (keep(dx), keep(dy), keep(dd));
                next[bi + bj * n] = avg;
                // Stash thresholded details in the free quadrants.
                next[bi + half + bj * n] = dx;
                next[bi + (bj + half) * n] = dy;
                next[bi + half + (bj + half) * n] = dd;
            }
        }
        a = next;
        size = half;
    }
    // Inverse transform.
    let mut size = 1;
    while size < n {
        let half = size;
        let mut next = a.clone();
        for bj in 0..half {
            for bi in 0..half {
                let avg = a[bi + bj * n];
                let dx = a[bi + half + bj * n];
                let dy = a[bi + (bj + half) * n];
                let dd = a[bi + half + (bj + half) * n];
                next[2 * bi + 2 * bj * n] = avg + dx + dy + dd;
                next[2 * bi + 1 + 2 * bj * n] = avg - dx + dy - dd;
                next[2 * bi + (2 * bj + 1) * n] = avg + dx - dy - dd;
                next[2 * bi + 1 + (2 * bj + 1) * n] = avg - dx - dy + dd;
            }
        }
        a = next;
        size *= 2;
    }
    ThresholdOutcome {
        field: a,
        kept,
        total,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_budget_roundtrips_exactly() {
        let n = 8;
        #[allow(clippy::cast_precision_loss)]
        let field: Vec<f64> = (0..n * n).map(|k| (k as f64).sin()).collect();
        let out = haar_threshold(&field, n, &|_, _| 0.0);
        for (a, b) in field.iter().zip(&out.field) {
            assert!((a - b).abs() < 1e-12, "lossless roundtrip");
        }
        assert_eq!(out.kept, out.total);
    }

    #[test]
    fn infinite_budget_keeps_only_the_mean() {
        let n = 4;
        let field: Vec<f64> = (0..16).map(f64::from).collect();
        let out = haar_threshold(&field, n, &|_, _| f64::INFINITY);
        assert_eq!(out.kept, 0);
        let mean: f64 = field.iter().sum::<f64>() / 16.0;
        for v in &out.field {
            assert!((v - mean).abs() < 1e-12);
        }
    }
}
