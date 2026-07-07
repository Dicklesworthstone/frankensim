//! Stiff machinery: a second-order IMEX step (implicit stiff LINEAR
//! part via a prefactored LU, explicit nonlinearity) and an exponential
//! Euler for u′ = A·u + N(u) with SYMMETRIC A (eigenbasis via the landed
//! Jacobi; φ₁ evaluated as expm1(x)/x — the cancellation-free form
//! fs-math's expm1 exists for). Krylov φ-actions for large nonsymmetric
//! A are recorded follow-up (needs Arnoldi).

use fs_la::eigen::jacobi_eigh;
use fs_la::factor::{Lu, lu};
use fs_math::det;

/// Prefactored operators for the IMEX-θ two-stage (ARS(2,2,2)-style)
/// scheme on u′ = L·u + N(u).
pub struct Imex2 {
    n: usize,
    h: f64,
    l_mat: Vec<f64>,
    solve_gamma: Lu,
    gamma: f64,
}

impl Imex2 {
    /// Build from the stiff linear operator (row-major n×n) and step h.
    ///
    /// # Panics
    /// If (I − γhL) is singular (h out of the scheme's range).
    #[must_use]
    pub fn new(l_mat: &[f64], n: usize, h: f64) -> Imex2 {
        let gamma = 1.0 - std::f64::consts::FRAC_1_SQRT_2; // ARS(2,2,2)
        let mut m = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..n {
                let id = if i == j { 1.0 } else { 0.0 };
                m[i * n + j] = (-gamma * h).mul_add(l_mat[i * n + j], id);
            }
        }
        let solve_gamma = lu(&m, n).expect("(I - gamma h L) must be nonsingular");
        Imex2 { n, h, l_mat: l_mat.to_vec(), solve_gamma, gamma }
    }
}

fn matvec(a: &[f64], n: usize, x: &[f64], out: &mut [f64]) {
    for i in 0..n {
        let mut acc = 0.0f64;
        for j in 0..n {
            acc = a[i * n + j].mul_add(x[j], acc);
        }
        out[i] = acc;
    }
}

/// One ARS(2,2,2) IMEX step: L treated implicitly (γ-diagonal, one
/// prefactored LU reused for both stages), N explicitly; second order
/// in BOTH parts and R(∞) = 0 on the stiff part. Updates `u` in place.
pub fn imex2_step<N: Fn(&[f64], &mut [f64])>(im: &Imex2, u: &mut [f64], nonlin: &N) {
    let (n, h, g) = (im.n, im.h, im.gamma);
    let mut nu = vec![0.0f64; n];
    nonlin(u, &mut nu);
    // Stage 1 (backward-Euler-γ on L, explicit N):
    // (I − γhL)·u₁ = u + γh·N(u).
    let mut rhs = vec![0.0f64; n];
    for i in 0..n {
        rhs[i] = h.mul_add(g * nu[i], u[i]);
    }
    im.solve_gamma.solve(&mut rhs);
    let u1 = rhs;
    // Stage 2 (ARS(2,2,2), stiffly accurate — u⁺ IS the last stage):
    // (I − γhL)·u⁺ = u + h·[δ·N(u) + (1−δ)·N(u₁) + (1−γ)·L·u₁],
    // δ = 1 − 1/(2γ). The (δ, 1−δ) explicit weights — NOT trapezoidal
    // (½, ½), which drops the nonlinear part to first order — satisfy
    // (1−δ)γ = ½, the h²·N′N order condition.
    let delta = 1.0 - 1.0 / (2.0 * g);
    let mut nu1 = vec![0.0f64; n];
    nonlin(&u1, &mut nu1);
    let mut lu1 = vec![0.0f64; n];
    matvec(&im.l_mat, n, &u1, &mut lu1);
    let mut rhs2 = vec![0.0f64; n];
    for i in 0..n {
        let nbar = delta.mul_add(nu[i], (1.0 - delta) * nu1[i]);
        rhs2[i] = u[i] + h * ((1.0 - g) * lu1[i] + nbar);
    }
    im.solve_gamma.solve(&mut rhs2);
    u.copy_from_slice(&rhs2);
}

/// Exponential Euler for u′ = A·u + N(u), SYMMETRIC A:
/// u⁺ = e^{hA}·u + h·φ₁(hA)·N(u), computed in A's eigenbasis with
/// φ₁(x) = expm1(x)/x (exact limit 1 at x = 0). EXACT for N ≡ 0.
pub struct ExpEuler {
    n: usize,
    h: f64,
    /// Eigenvectors (columns) of A.
    vecs: Vec<f64>,
    /// e^{hλ} per eigenvalue.
    exp_h: Vec<f64>,
    /// h·φ₁(hλ) per eigenvalue.
    hphi1: Vec<f64>,
}

impl ExpEuler {
    /// Build from symmetric A (row-major n×n) and step h.
    #[must_use]
    pub fn new(a: &[f64], n: usize, h: f64) -> ExpEuler {
        let (vals, vecs) = jacobi_eigh(a, n);
        let exp_h: Vec<f64> = vals.iter().map(|&l| det::exp(h * l)).collect();
        let hphi1: Vec<f64> = vals
            .iter()
            .map(|&l| {
                let x = h * l;
                if x.abs() < 1e-300 { h } else { h * (det::expm1(x) / x) }
            })
            .collect();
        ExpEuler { n, h, vecs, exp_h, hphi1 }
    }

    /// The step size.
    #[must_use]
    pub fn h(&self) -> f64 {
        self.h
    }

    /// One exponential-Euler step (u updated in place).
    pub fn step<N: Fn(&[f64], &mut [f64])>(&self, u: &mut [f64], nonlin: &N) {
        let n = self.n;
        let mut nu = vec![0.0f64; n];
        nonlin(u, &mut nu);
        // Transform to eigenbasis: û = Vᵀu, n̂ = VᵀN.
        let mut uh = vec![0.0f64; n];
        let mut nh = vec![0.0f64; n];
        for i in 0..n {
            let (mut au, mut an) = (0.0f64, 0.0f64);
            for j in 0..n {
                au = self.vecs[j * n + i].mul_add(u[j], au);
                an = self.vecs[j * n + i].mul_add(nu[j], an);
            }
            uh[i] = au;
            nh[i] = an;
        }
        // Apply the scalar filters and transform back.
        for i in 0..n {
            uh[i] = self.exp_h[i].mul_add(uh[i], self.hphi1[i] * nh[i]);
        }
        for j in 0..n {
            let mut acc = 0.0f64;
            for i in 0..n {
                acc = self.vecs[j * n + i].mul_add(uh[i], acc);
            }
            u[j] = acc;
        }
    }
}
