//! Karhunen–Loève expansions of random fields on point sets: the
//! covariance matrix's eigendecomposition (fs-la Jacobi), truncated
//! at a caller-chosen captured-variance fraction — and the fraction
//! actually captured is REPORTED, not assumed. Field samples come
//! from scrambled-Sobol Gaussian germs (fs-bo Φ⁻¹): deterministic
//! per seed.

use fs_la::eigen::jacobi_eigh;

/// Stationary covariance families for field construction.
#[derive(Debug, Clone, Copy)]
pub enum CovarianceKind {
    /// Exponential: σ²·exp(−r/ℓ) (rough fields).
    Exponential,
    /// Squared-exponential: σ²·exp(−r²/2ℓ²) (smooth fields).
    SquaredExponential,
}

/// A truncated KL expansion on a fixed point set.
pub struct KlExpansion {
    /// Retained eigenvalues (descending).
    pub eigenvalues: Vec<f64>,
    /// Retained modes (row-major: mode k at `modes[k*n..(k+1)*n]`).
    pub modes: Vec<f64>,
    /// Fraction of total variance captured by the truncation — the
    /// EVIDENCE the bead requires.
    pub captured_variance: f64,
    /// Point count.
    pub n: usize,
}

fn covariance(kind: CovarianceKind, sigma2: f64, ell: f64, a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let d2: f64 = a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum();
    let r = fs_math::det::sqrt(d2);
    match kind {
        CovarianceKind::Exponential => sigma2 * fs_math::det::exp(-r / ell),
        CovarianceKind::SquaredExponential => sigma2 * fs_math::det::exp(-d2 / (2.0 * ell * ell)),
    }
}

impl KlExpansion {
    /// Build by dense eigendecomposition of the covariance matrix at
    /// `points`, truncating at the smallest K whose eigenvalue mass
    /// reaches `variance_target` (of the trace).
    #[must_use]
    pub fn build(
        points: &[[f64; 3]],
        kind: CovarianceKind,
        sigma2: f64,
        ell: f64,
        variance_target: f64,
    ) -> KlExpansion {
        let n = points.len();
        let mut c = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..=i {
                let v = covariance(kind, sigma2, ell, &points[i], &points[j]);
                c[i * n + j] = v;
                c[j * n + i] = v;
            }
        }
        let trace: f64 = (0..n).map(|i| c[i * n + i]).sum();
        // jacobi_eigh returns ascending; we retain from the top.
        let (vals, vecs) = jacobi_eigh(&c, n);
        let mut eigenvalues = Vec::new();
        let mut modes = Vec::new();
        let mut captured = 0.0f64;
        for k in (0..n).rev() {
            if captured / trace >= variance_target {
                break;
            }
            let lam = vals[k].max(0.0);
            captured += lam;
            eigenvalues.push(lam);
            // Column k of the row-major eigenvector matrix.
            modes.extend((0..n).map(|r| vecs[r * n + k]));
        }
        KlExpansion {
            eigenvalues,
            modes,
            captured_variance: captured / trace,
            n,
        }
    }

    /// Number of retained modes.
    #[must_use]
    pub fn order(&self) -> usize {
        self.eigenvalues.len()
    }

    /// Realize the field from a Gaussian germ vector ξ (one entry per
    /// retained mode): f = Σₖ √λₖ·ξₖ·φₖ.
    #[must_use]
    pub fn realize(&self, xi: &[f64]) -> Vec<f64> {
        assert_eq!(xi.len(), self.order());
        let mut field = vec![0.0f64; self.n];
        for (k, (&lam, &x)) in self.eigenvalues.iter().zip(xi).enumerate() {
            let scale = fs_math::det::sqrt(lam) * x;
            for (f, m) in field
                .iter_mut()
                .zip(&self.modes[k * self.n..(k + 1) * self.n])
            {
                *f = scale.mul_add(*m, *f);
            }
        }
        field
    }

    /// Deterministic Gaussian germs: scrambled Sobol through Φ⁻¹ on
    /// the LEADING ≤ MAX_SOBOL_DIM modes (the variance-dominant germ
    /// directions, where QMC pays; the embedded Joe–Kuo table caps at
    /// 10 dims — larger tables are fs-rand's recorded follow-up) and
    /// Philox normals for the tail modes.
    #[must_use]
    pub fn qmc_germs(&self, count: usize, seed: u64) -> Vec<Vec<f64>> {
        let k = self.order();
        let kq = k.min(fs_rand::qmc::MAX_SOBOL_DIM);
        let sobol = fs_rand::qmc::Sobol::scrambled(kq, seed);
        let mut pt = vec![0.0f64; kq];
        (0..count)
            .map(|s| {
                sobol.point(u32::try_from(s + 1).expect("count fits u32"), &mut pt);
                let mut germ: Vec<f64> = pt
                    .iter()
                    .map(|u| fs_bo::phi_inv(u.clamp(1e-12, 1.0 - 1e-12)))
                    .collect();
                if k > kq {
                    let mut tail = fs_rand::StreamKey {
                        seed,
                        kernel: 0x0517,
                        tile: u32::try_from(s).expect("count fits u32"),
                    }
                    .stream();
                    germ.extend((kq..k).map(|_| tail.next_normal()));
                }
                germ
            })
            .collect()
    }

    /// Reconstruction error of the truncated covariance:
    /// ‖Σₖ λₖφₖφₖᵀ − C‖_F / ‖C‖_F (the truncation-quality audit; equals
    /// √(Σ dropped λ²)/‖C‖_F for exact arithmetic).
    #[must_use]
    pub fn covariance_reconstruction_error(
        &self,
        points: &[[f64; 3]],
        kind: CovarianceKind,
        sigma2: f64,
        ell: f64,
    ) -> f64 {
        let n = self.n;
        let mut err2 = 0.0f64;
        let mut norm2 = 0.0f64;
        for i in 0..n {
            for j in 0..n {
                let c = covariance(kind, sigma2, ell, &points[i], &points[j]);
                let mut rec = 0.0f64;
                for (k, &lam) in self.eigenvalues.iter().enumerate() {
                    rec = (lam * self.modes[k * n + i]).mul_add(self.modes[k * n + j], rec);
                }
                err2 += (c - rec) * (c - rec);
                norm2 += c * c;
            }
        }
        fs_math::det::sqrt(err2 / norm2)
    }
}
