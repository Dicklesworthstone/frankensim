//! Generalized-α (Chung–Hulbert) for structural dynamics
//! M·q̈ + C·q̇ + K·q = f(t): second-order accurate with CONTROLLABLE
//! high-frequency dissipation via ρ∞ ∈ [0, 1] (ρ∞ = 1: no dissipation;
//! ρ∞ = 0: annihilate the highest mode in one step). The spectral
//! behavior is TESTED against theory in the battery, not just cited.

use fs_la::factor::{Lu, lu};

/// Prefactored generalized-α stepper for fixed (M, C, K, h).
pub struct GeneralizedAlpha {
    n: usize,
    h: f64,
    alpha_m: f64,
    alpha_f: f64,
    beta: f64,
    gamma: f64,
    m_mat: Vec<f64>,
    c_mat: Vec<f64>,
    k_mat: Vec<f64>,
    eff: Lu,
}

impl GeneralizedAlpha {
    /// Build from mass/damping/stiffness (row-major n×n) and ρ∞.
    ///
    /// # Panics
    /// Structured panic if the effective matrix is singular (a modeling
    /// error: h and the system matrices are incompatible).
    #[must_use]
    pub fn new(
        m_mat: &[f64],
        c_mat: &[f64],
        k_mat: &[f64],
        n: usize,
        h: f64,
        rho_inf: f64,
    ) -> GeneralizedAlpha {
        assert!((0.0..=1.0).contains(&rho_inf), "rho_inf in [0,1]");
        // Chung–Hulbert parameterization.
        let alpha_m = (2.0 * rho_inf - 1.0) / (rho_inf + 1.0);
        let alpha_f = rho_inf / (rho_inf + 1.0);
        let gamma = 0.5 - alpha_m + alpha_f;
        let beta = 0.25 * (1.0 - alpha_m + alpha_f) * (1.0 - alpha_m + alpha_f);
        // Effective matrix: (1−αm)/(βh²)·M + (1−αf)γ/(βh)·C + (1−αf)·K.
        let cm = (1.0 - alpha_m) / (beta * h * h);
        let cc = (1.0 - alpha_f) * gamma / (beta * h);
        let ck = 1.0 - alpha_f;
        let mut eff = vec![0.0f64; n * n];
        for i in 0..n * n {
            eff[i] = cm.mul_add(m_mat[i], cc.mul_add(c_mat[i], ck * k_mat[i]));
        }
        let eff = lu(&eff, n).expect("generalized-alpha effective matrix must be nonsingular");
        GeneralizedAlpha {
            n,
            h,
            alpha_m,
            alpha_f,
            beta,
            gamma,
            m_mat: m_mat.to_vec(),
            c_mat: c_mat.to_vec(),
            k_mat: k_mat.to_vec(),
            eff,
        }
    }

    /// The (γ, β) Newmark parameters in use (diagnostics).
    #[must_use]
    pub fn newmark(&self) -> (f64, f64) {
        (self.gamma, self.beta)
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

/// One generalized-α step: (q, v, a) at t → t+h with load `f_next`
/// evaluated at t + (1−αf)·h by the caller (constant loads just pass
/// the value). Updates in place.
pub fn galpha_step(
    ga: &GeneralizedAlpha,
    q: &mut [f64],
    v: &mut [f64],
    a: &mut [f64],
    f_next: &[f64],
) {
    let n = ga.n;
    let (h, am, af, beta, gamma) =
        (ga.h, ga.alpha_m, ga.alpha_f, ga.beta, ga.gamma);
    // Predictors (Newmark form).
    let cm = (1.0 - am) / (beta * h * h);
    let cc = (1.0 - af) * gamma / (beta * h);
    // RHS = (1−αf)·f_next + αf·(K·q-term folds) … assembled explicitly:
    // r = (1−αf) f + M·[cm·q + (1−am)/(βh)·v + ((1−am)/(2β) − 1)·a]
    //   + C·[cc·q + ((1−αf)γ/β − 1)·v + (1−αf)h·(γ/(2β) − 1)·a]
    //   − αf·K·q.
    let mv_c1 = (1.0 - am) / (beta * h);
    let mv_c2 = (1.0 - am) / (2.0 * beta) - 1.0;
    let cv_c1 = (1.0 - af) * gamma / beta - 1.0;
    let cv_c2 = (1.0 - af) * h * (gamma / (2.0 * beta) - 1.0);
    let mut tm = vec![0.0f64; n];
    let mut tc = vec![0.0f64; n];
    let mut tk = vec![0.0f64; n];
    let mvec: Vec<f64> = (0..n)
        .map(|i| cm.mul_add(q[i], mv_c1.mul_add(v[i], mv_c2 * a[i])))
        .collect();
    let cvec: Vec<f64> = (0..n)
        .map(|i| cc.mul_add(q[i], cv_c1.mul_add(v[i], cv_c2 * a[i])))
        .collect();
    matvec(&ga.m_mat, n, &mvec, &mut tm);
    matvec(&ga.c_mat, n, &cvec, &mut tc);
    matvec(&ga.k_mat, n, q, &mut tk);
    let mut rhs = vec![0.0f64; n];
    for i in 0..n {
        rhs[i] = (1.0 - af).mul_add(f_next[i], tm[i] + tc[i] - af * tk[i]);
    }
    ga.eff.solve(&mut rhs); // rhs now holds q_{n+1}
    // Newmark corrector for a_{n+1}, v_{n+1}.
    for i in 0..n {
        let dq = rhs[i] - q[i];
        let a_new = (dq / (beta * h * h)) - v[i] / (beta * h) - (0.5 / beta - 1.0) * a[i];
        let v_new = v[i] + h * ((1.0 - gamma) * a[i] + gamma * a_new);
        q[i] = rhs[i];
        v[i] = v_new;
        a[i] = a_new;
    }
}
