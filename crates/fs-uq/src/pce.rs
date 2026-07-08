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

/// Multi-index set of total degree ≤ p in d germs (graded-lex order —
/// deterministic basis enumeration).
fn total_degree_indices(d: usize, p: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut idx = vec![0usize; d];
    loop {
        let deg: usize = idx.iter().sum();
        if deg <= p {
            out.push(idx.clone());
        }
        // Odometer increment with degree cap p per variable.
        let mut pos = 0;
        loop {
            if pos == d {
                out.sort_by_key(|m| (m.iter().sum::<usize>(), m.clone()));
                return out;
            }
            idx[pos] += 1;
            if idx[pos] <= p {
                break;
            }
            idx[pos] = 0;
            pos += 1;
        }
    }
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
    let d = xi[0].len();
    let indices = total_degree_indices(d, p);
    let m = indices.len();
    assert!(
        n >= 2 * m,
        "PCE regression wants n >= 2*basis ({n} vs {m} basis functions)"
    );
    // Design matrix A (n×m).
    let mut a = vec![0.0f64; n * m];
    for (i, x) in xi.iter().enumerate() {
        for (j, idx) in indices.iter().enumerate() {
            let mut phi = 1.0f64;
            for (k, &xv) in idx.iter().zip(x) {
                phi *= hermite_orthonormal(*k, xv);
            }
            a[i * m + j] = phi;
        }
    }
    // Normal equations AᵀA c = Aᵀy with tiny ridge.
    let mut ata = vec![0.0f64; m * m];
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
    let chol = fs_la::factor::cholesky(&ata, m).expect("PCE normal equations SPD");
    let mut c = aty;
    chol.solve(&mut c);
    PceModel {
        indices,
        coefficients: c,
        dim: d,
    }
}
