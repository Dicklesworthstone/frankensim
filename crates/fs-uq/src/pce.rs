//! Polynomial chaos by regression: probabilists' Hermite basis in
//! Gaussian germs, total-degree truncation, least squares via fs-la
//! Cholesky normal equations on QMC designs. Mean and variance drop
//! out of the coefficients (orthonormal basis) — verified against
//! closed forms in the battery.

/// Probabilists' Hermite Heₖ, ORTHONORMALIZED (divided by √k!):
/// E[hᵢ(ξ)hⱼ(ξ)] = δᵢⱼ under ξ ~ N(0,1).
#[must_use]
pub fn hermite_orthonormal(k: usize, x: f64) -> f64 {
    // Three-term recurrence on the monic Heₖ, then normalize.
    let mut h_prev = 1.0f64;
    if k == 0 {
        return 1.0;
    }
    let mut h = x;
    for j in 1..k {
        let next = x.mul_add(h, -(j as f64) * h_prev);
        h_prev = h;
        h = next;
    }
    let mut fact = 1.0f64;
    for j in 2..=k {
        fact *= j as f64;
    }
    h / fs_math::det::sqrt(fact)
}

/// Number of total-degree basis terms: binomial(dim + degree, degree).
///
/// Returns `None` if the count does not fit in `usize`. This is a cheap
/// preflight for the regression's sample and storage requirements; it does
/// not enumerate a basis. Zero germs or degree zero need only a constant.
#[must_use]
pub fn total_degree_basis_size(dim: usize, degree: usize) -> Option<usize> {
    if dim == 0 || degree == 0 {
        return Some(1);
    }
    let n = dim.checked_add(degree)?;
    let k = dim.min(degree);
    let mut count = 1usize;
    for j in 1..=k {
        let numerator = n - k + j;
        // Cancel before multiplying: an intermediate product can overflow
        // even when the binomial coefficient itself is representable.
        let mut a = numerator;
        let mut b = j;
        while b != 0 {
            (a, b) = (b, a % b);
        }
        count /= j / a;
        count = count.checked_mul(numerator / a)?;
    }
    Some(count)
}

/// Enumerate weak compositions directly, in the existing graded-lex order.
/// Work is proportional to the output size, not the enclosing (p + 1)^d
/// Cartesian grid. The iterative successor also avoids a dimension-deep
/// recursion stack and a final sort.
fn total_degree_indices(d: usize, p: usize) -> Vec<Vec<usize>> {
    let count = total_degree_basis_size(d, p).expect("PCE basis size overflow");
    let mut out = Vec::with_capacity(count);
    let mut idx = vec![0usize; d];
    if d == 0 {
        out.push(idx);
        return out;
    }
    for degree in 0..=p {
        idx[d - 1] = degree;
        loop {
            out.push(idx.clone());
            // Move one unit left from the rightmost nonzero exponent,
            // placing the remaining tail in the last coordinate.
            let Some(j) = (1..d).rev().find(|&j| idx[j] != 0) else {
                break;
            };
            let tail = idx[j] - 1;
            idx[j] = 0;
            idx[j - 1] += 1;
            idx[d - 1] = tail;
        }
        idx.fill(0);
    }
    debug_assert_eq!(out.len(), count);
    out
}

/// A fitted PCE surrogate.
pub struct PceModel {
    /// Basis multi-indices (graded-lex).
    pub indices: Vec<Vec<usize>>,
    /// Coefficients (aligned with `indices`).
    pub coefficients: Vec<f64>,
    /// Germ dimension.
    pub dim: usize,
}

impl PceModel {
    /// Evaluate the surrogate at a germ point.
    #[must_use]
    pub fn eval(&self, xi: &[f64]) -> f64 {
        self.indices
            .iter()
            .zip(&self.coefficients)
            .map(|(m, c)| {
                let mut phi = 1.0f64;
                for (k, &x) in m.iter().zip(xi) {
                    phi *= hermite_orthonormal(*k, x);
                }
                c * phi
            })
            .sum()
    }

    /// Mean = coefficient of the constant basis function.
    #[must_use]
    pub fn mean(&self) -> f64 {
        self.coefficients[0]
    }

    /// Variance = Σ non-constant coefficients² (orthonormal basis).
    #[must_use]
    pub fn variance(&self) -> f64 {
        self.coefficients[1..].iter().map(|c| c * c).sum()
    }
}

/// Fit a total-degree-`p` PCE to samples (ξᵢ, yᵢ) by least squares
/// (normal equations + ridge 1e−12 through fs-la Cholesky). The
/// design should oversample the basis (n ≥ 2·|basis| is the usual
/// rule; asserted).
#[must_use]
pub fn fit_pce(xi: &[Vec<f64>], y: &[f64], p: usize) -> PceModel {
    let n = xi.len();
    assert!(n != 0, "PCE regression requires samples");
    assert_eq!(n, y.len(), "PCE sample/response count mismatch");
    let d = xi[0].len();
    // Refuse under-sampled/overflowing requests BEFORE basis construction.
    let m = total_degree_basis_size(d, p).expect("PCE basis size overflow");
    assert!(
        m <= n / 2,
        "PCE regression wants n >= 2*basis ({n} vs {m} basis functions)"
    );
    for x in xi {
        assert_eq!(x.len(), d, "PCE sample dimension mismatch");
        assert!(x.iter().all(|v| v.is_finite()), "PCE non-finite germ");
    }
    assert!(y.iter().all(|v| v.is_finite()), "PCE non-finite response");
    let max_entries = (isize::MAX as usize) / core::mem::size_of::<f64>();
    let design_len = n.checked_mul(m).filter(|&len| len <= max_entries)
        .expect("PCE design matrix size overflow");
    let gram_len = m.checked_mul(m).filter(|&len| len <= max_entries)
        .expect("PCE normal matrix size overflow");
    let indices = total_degree_indices(d, p);
    // Design matrix A (n×m).
    let mut a = vec![0.0f64; design_len];
    for (i, x) in xi.iter().enumerate() {
        for (j, idx) in indices.iter().enumerate() {
            let mut phi = 1.0f64;
            for (k, &xv) in idx.iter().zip(x) {
                phi *= hermite_orthonormal(*k, xv);
            }
            assert!(phi.is_finite(), "PCE non-finite basis value");
            a[i * m + j] = phi;
        }
    }
    // Normal equations AᵀA c = Aᵀy with tiny ridge.
    let mut ata = vec![0.0f64; gram_len];
    let mut aty = vec![0.0f64; m];
    for i in 0..n {
        for j in 0..m {
            aty[j] = a[i * m + j].mul_add(y[i], aty[j]);
            for k in 0..=j {
                ata[j * m + k] = a[i * m + j].mul_add(a[i * m + k], ata[j * m + k]);
            }
        }
    }
    for j in 0..m {
        for k in 0..j {
            ata[k * m + j] = ata[j * m + k];
        }
        ata[j * m + j] += 1e-12;
    }
    assert!(
        ata.iter().chain(&aty).all(|v| v.is_finite()),
        "PCE normal equations overflow"
    );
    let chol = fs_la::factor::cholesky(&ata, m).expect("PCE normal equations SPD");
    let mut c = aty;
    chol.solve(&mut c);
    assert!(c.iter().all(|v| v.is_finite()), "PCE non-finite fit");
    PceModel {
        indices,
        coefficients: c,
        dim: d,
    }
}

#[cfg(test)]
mod basis_tests {
    use super::{total_degree_basis_size, total_degree_indices};

    // Independent, deliberately small exhaustive reference for the old
    // Cartesian enumeration and sort contract. Never used on large inputs.
    fn cartesian_reference(d: usize, p: usize) -> Vec<Vec<usize>> {
        let mut rows = vec![Vec::new()];
        for _ in 0..d {
            rows = rows.into_iter().flat_map(|prefix| {
                (0..=p).map(move |k| {
                    let mut row = prefix.clone();
                    row.push(k);
                    row
                })
            }).collect();
        }
        rows.retain(|row| row.iter().sum::<usize>() <= p);
        rows.sort_by_key(|row| (row.iter().sum::<usize>(), row.clone()));
        rows
    }

    #[test]
    fn direct_enumeration_preserves_every_small_basis_and_order() {
        for d in 0..=5 {
            for p in 0..=5 {
                let actual = total_degree_indices(d, p);
                assert_eq!(actual, cartesian_reference(d, p), "d={d}, p={p}");
                assert_eq!(Some(actual.len()), total_degree_basis_size(d, p));
            }
        }
    }

    #[test]
    fn sparse_high_dimensional_bases_do_not_scan_a_cartesian_grid() {
        // The old scan visits 3^20 = 3,486,784,401 candidates here.
        let quadratic = total_degree_indices(20, 2);
        assert_eq!(quadratic.len(), 231);
        assert!(quadratic.iter().all(|row| row.len() == 20));
        let linear = total_degree_indices(256, 1);
        assert_eq!(linear.len(), 257);
        assert_eq!(total_degree_indices(0, usize::MAX), vec![Vec::<usize>::new()]);
    }

    #[test]
    fn basis_count_checks_overflow_without_false_intermediate_overflow() {
        assert_eq!(total_degree_basis_size(0, usize::MAX), Some(1));
        assert_eq!(total_degree_basis_size(usize::MAX, 0), Some(1));
        assert_eq!(total_degree_basis_size(usize::MAX, 1), None);
        assert_eq!(total_degree_basis_size(usize::MAX / 2, 2), None);
        assert_eq!(total_degree_basis_size(10, 10), Some(184_756));
        #[cfg(target_pointer_width = "64")]
        {
            assert_eq!(total_degree_basis_size(34, 33), Some(14_226_520_737_620_288_370));
            assert_eq!(total_degree_basis_size(34, 34), None);
        }
    }
}
