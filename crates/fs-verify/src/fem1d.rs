//! The 1D elliptic testbed: P1 finite elements for `−u″ = f` on (0,1)
//! with homogeneous Dirichlet data, polynomial manufactured solutions
//! (so flux antiderivatives are exact and Gauss quadrature is
//! MATHEMATICALLY exact — the rigor of the verifier's bound then rests
//! only on outward-rounded interval evaluation), and a Newton solver
//! for the nonlinear warm-start class.

use crate::interval::Iv;

/// A polynomial in monomial coefficients (`c[0] + c[1]x + …`).
#[derive(Debug, Clone, PartialEq)]
pub struct Poly(pub Vec<f64>);

impl Poly {
    /// Derivative.
    #[must_use]
    pub fn derive(&self) -> Poly {
        if self.0.len() <= 1 {
            return Poly(vec![0.0]);
        }
        Poly(
            self.0[1..]
                .iter()
                .enumerate()
                .map(|(k, &c)| c * (k + 1) as f64)
                .collect(),
        )
    }

    /// Antiderivative with zero constant term.
    #[must_use]
    pub fn antiderive(&self) -> Poly {
        let mut out = vec![0.0];
        out.extend(self.0.iter().enumerate().map(|(k, &c)| c / (k + 1) as f64));
        Poly(out)
    }

    /// Negate.
    #[must_use]
    pub fn neg(&self) -> Poly {
        Poly(self.0.iter().map(|c| -c).collect())
    }

    /// Horner evaluation (f64).
    #[must_use]
    pub fn eval(&self, x: f64) -> f64 {
        self.0.iter().rev().fold(0.0, |acc, &c| acc * x + c)
    }

    /// Horner evaluation with outward-rounded intervals.
    #[must_use]
    pub fn eval_iv(&self, x: Iv) -> Iv {
        self.0
            .iter()
            .rev()
            .fold(Iv::zero(), |acc, &c| acc.mul(x).add(Iv::point(c)))
    }

    /// Polynomial degree (0 for constants).
    #[must_use]
    pub fn degree(&self) -> usize {
        self.0.len().saturating_sub(1)
    }
}

/// A manufactured problem: `−u″ = f`, `u(0) = u(1) = 0`.
#[derive(Debug, Clone)]
pub struct MmsProblem {
    /// Name (ledger rows).
    pub name: String,
    /// The exact solution (must vanish at 0 and 1).
    pub u: Poly,
    /// The load `f = −u″`.
    pub f: Poly,
    /// `F(x) = ∫₀ˣ f` (exact antiderivative; the equilibrated flux is
    /// `σ = c − F` for a free constant c).
    pub big_f: Poly,
    /// Mesh nodes (ascending, first 0, last 1).
    pub mesh: Vec<f64>,
}

impl MmsProblem {
    /// Build from an exact polynomial solution and a mesh.
    #[must_use]
    pub fn new(name: &str, u: Poly, mesh: Vec<f64>) -> Self {
        let f = u.derive().derive().neg();
        let big_f = f.antiderive();
        MmsProblem {
            name: name.to_string(),
            u,
            f,
            big_f,
            mesh,
        }
    }
}

/// Solve `−u″ = f` with P1 elements on the problem's mesh (Thomas
/// algorithm; deterministic). Returns interior+boundary nodal values.
#[must_use]
pub fn solve_p1(problem: &MmsProblem) -> Vec<f64> {
    let m = &problem.mesh;
    let n = m.len();
    if n < 3 {
        return vec![0.0; n];
    }
    let interior = n - 2;
    // Tridiagonal stiffness + exact load via 5-pt Gauss per element.
    let mut diag = vec![0.0; interior];
    let mut off = vec![0.0; interior.saturating_sub(1)];
    let mut rhs = vec![0.0; interior];
    for e in 0..n - 1 {
        let (x0, x1) = (m[e], m[e + 1]);
        let h = x1 - x0;
        let k = 1.0 / h;
        // Assemble stiffness.
        for (a, b) in [(e, e), (e + 1, e + 1)] {
            if (1..=interior).contains(&a) && (1..=interior).contains(&b) {
                diag[a - 1] += k;
            }
        }
        if e >= 1 && e < interior {
            off[e - 1] -= k;
        }
        // Load: ∫ f φ_a over the element, exact Gauss.
        for (gx, gw) in gauss5(x0, x1) {
            let fv = problem.f.eval(gx);
            let phi_left = (x1 - gx) / h;
            let phi_right = (gx - x0) / h;
            if e >= 1 {
                rhs[e - 1] += gw * fv * phi_left;
            }
            if e < interior {
                rhs[e] += gw * fv * phi_right;
            }
        }
    }
    // Thomas solve.
    let mut c = vec![0.0; interior];
    let mut d = vec![0.0; interior];
    c[0] = if interior > 1 { off[0] / diag[0] } else { 0.0 };
    d[0] = rhs[0] / diag[0];
    for i in 1..interior {
        let m_ = diag[i] - off[i - 1] * c[i - 1];
        if i < interior - 1 {
            c[i] = off[i] / m_;
        }
        d[i] = (rhs[i] - off[i - 1] * d[i - 1]) / m_;
    }
    let mut x = vec![0.0; interior];
    x[interior - 1] = d[interior - 1];
    for i in (0..interior - 1).rev() {
        x[i] = d[i] - c[i] * x[i + 1];
    }
    let mut full = vec![0.0; n];
    full[1..=interior].copy_from_slice(&x);
    full
}

/// 5-point Gauss–Legendre nodes/weights mapped to `[x0, x1]` (exact
/// for polynomial degree ≤ 9).
#[must_use]
pub fn gauss5(x0: f64, x1: f64) -> [(f64, f64); 5] {
    const N: [f64; 5] = [
        -0.906_179_845_938_664,
        -0.538_469_310_105_683,
        0.0,
        0.538_469_310_105_683,
        0.906_179_845_938_664,
    ];
    const W: [f64; 5] = [
        0.236_926_885_056_189,
        0.478_628_670_499_366,
        0.568_888_888_888_889,
        0.478_628_670_499_366,
        0.236_926_885_056_189,
    ];
    let mid = f64::midpoint(x0, x1);
    let half = 0.5 * (x1 - x0);
    core::array::from_fn(|i| (mid + half * N[i], half * W[i]))
}

/// True energy-norm error `‖u′ − u_h′‖` (the ORACLE; high-resolution
/// f64 quadrature — the oracle needs accuracy, not rigor).
#[must_use]
pub fn true_energy_error(problem: &MmsProblem, candidate: &[f64]) -> f64 {
    let m = &problem.mesh;
    let du = problem.u.derive();
    let mut acc = 0.0;
    for e in 0..m.len() - 1 {
        let (x0, x1) = (m[e], m[e + 1]);
        let h = x1 - x0;
        let slope = (candidate[e + 1] - candidate[e]) / h;
        // Subdivide each element for oracle accuracy.
        let sub = 32;
        for s in 0..sub {
            let a = x0 + h * f64::from(s) / f64::from(sub);
            let b = x0 + h * f64::from(s + 1) / f64::from(sub);
            for (gx, gw) in gauss5(a, b) {
                let d = du.eval(gx) - slope;
                acc += gw * d * d;
            }
        }
    }
    acc.sqrt()
}

/// Newton solve of the toy NONLINEAR class `−u″ + u³ = f` starting
/// from `start`; returns (solution, iterations to ‖R‖ < 1e-10).
#[must_use]
#[allow(clippy::too_many_lines)] // assembly + Newton loop are one method
pub fn solve_nonlinear(problem: &MmsProblem, start: &[f64], max_iter: u32) -> (Vec<f64>, u32) {
    let m = &problem.mesh;
    let n = m.len();
    let interior = n - 2;
    let mut u = start.to_vec();
    let mut iters = 0;
    for _ in 0..max_iter {
        // Residual R_i = (A u)_i + (M u³)_i − b_i with lumped mass.
        let mut resid = vec![0.0; interior];
        let mut jac_diag = vec![0.0; interior];
        let mut jac_off = vec![0.0; interior.saturating_sub(1)];
        for e in 0..n - 1 {
            let (x0, x1) = (m[e], m[e + 1]);
            let h = x1 - x0;
            let k = 1.0 / h;
            let slope_term = (u[e + 1] - u[e]) * k;
            if e >= 1 {
                resid[e - 1] -= slope_term;
                jac_diag[e - 1] += k;
            }
            if e < interior {
                resid[e] += slope_term;
                jac_diag[e] += k;
            }
            if e >= 1 && e < interior {
                jac_off[e - 1] -= k;
            }
            // Lumped nonlinear + load terms.
            let hl = 0.5 * h;
            if e >= 1 {
                resid[e - 1] += hl * u[e].powi(3);
                jac_diag[e - 1] += hl * 3.0 * u[e] * u[e];
            }
            if e < interior {
                resid[e] += hl * u[e + 1].powi(3);
                jac_diag[e] += hl * 3.0 * u[e + 1] * u[e + 1];
            }
            for (gx, gw) in gauss5(x0, x1) {
                let fv = problem.f.eval(gx) + problem.u.eval(gx).powi(3);
                let phi_l = (x1 - gx) / h;
                let phi_r = (gx - x0) / h;
                if e >= 1 {
                    resid[e - 1] -= gw * fv * phi_l;
                }
                if e < interior {
                    resid[e] -= gw * fv * phi_r;
                }
            }
        }
        let norm: f64 = resid.iter().map(|r| r * r).sum::<f64>().sqrt();
        if norm < 1e-10 {
            break;
        }
        iters += 1;
        // Thomas solve J δ = −R.
        let rhs: Vec<f64> = resid.iter().map(|r| -r).collect();
        let mut c = vec![0.0; interior];
        let mut d = vec![0.0; interior];
        c[0] = if interior > 1 {
            jac_off[0] / jac_diag[0]
        } else {
            0.0
        };
        d[0] = rhs[0] / jac_diag[0];
        for i in 1..interior {
            let mm = jac_diag[i] - jac_off[i - 1] * c[i - 1];
            if i < interior - 1 {
                c[i] = jac_off[i] / mm;
            }
            d[i] = (rhs[i] - jac_off[i - 1] * d[i - 1]) / mm;
        }
        let mut delta = vec![0.0; interior];
        delta[interior - 1] = d[interior - 1];
        for i in (0..interior - 1).rev() {
            delta[i] = d[i] - c[i] * delta[i + 1];
        }
        for i in 0..interior {
            u[i + 1] += delta[i];
        }
    }
    (u, iters)
}
