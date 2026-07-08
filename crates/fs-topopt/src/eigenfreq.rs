//! Eigenfrequency objectives: K(ρ̄)φ = λ·M(ρ̄)φ with exact eigenvalue
//! derivatives through the full design chain, the MASS-INTERPOLATION
//! trap handled (SIMP stiffness with linear mass gives λ ∝ ρ^{p−1} → 0
//! spurious modes in void regions — below ρ = 0.1 the mass follows
//! c·ρ⁶, continuously matched, so void cells cannot host artificial
//! low modes), and CLUSTERED eigenvalues aggregated by a smooth
//! minimum (single-eigenvalue objectives are nonsmooth at crossings —
//! the classic trap; the aggregate has an exact weighted-combination
//! gradient, FD-verified near a designed crossing in the battery).

use crate::elasticity::DensityElasticity;
use crate::pipeline::DesignPipeline;

/// Mass interpolation m(ρ̄): linear above the void threshold, ρ⁶-class
/// below (continuous at the threshold).
#[must_use]
pub fn mass_interp(rho_bar: f64) -> f64 {
    const T: f64 = 0.1;
    let r = rho_bar.clamp(0.0, 1.0);
    if r >= T {
        r
    } else {
        // c·r⁶ with c = T / T⁶ = T^{-5} for continuity at T.
        let c = fs_math::det::pow(T, -5.0);
        c * fs_math::det::pow(r, 6.0)
    }
}

/// dm/dρ̄.
#[must_use]
pub fn mass_interp_derivative(rho_bar: f64) -> f64 {
    const T: f64 = 0.1;
    let r = rho_bar.clamp(0.0, 1.0);
    if r >= T {
        1.0
    } else {
        let c = fs_math::det::pow(T, -5.0);
        6.0 * c * fs_math::det::pow(r, 5.0)
    }
}

/// The lowest `count` eigenpairs of K(ρ̄)φ = λM(ρ̄)φ on the free dofs
/// (dense, fixture scale): M = LLᵀ Cholesky reduction to standard
/// form, fs-la Jacobi, back-transformed vectors M-normalized
/// (φᵀMφ = 1).
#[must_use]
pub fn lowest_eigenpairs(
    elasticity: &DensityElasticity,
    mass_densities: &[f64],
    count: usize,
) -> (Vec<f64>, Vec<Vec<f64>>, Vec<usize>) {
    let (kd, md, free_idx) = elasticity.assemble_dense(mass_densities);
    let f = free_idx.len();
    let chol = fs_la::factor::cholesky(&md, f).expect("mass matrix SPD");
    // A = L⁻¹ K L⁻ᵀ: column-solve then row-solve.
    // L⁻¹K: forward-substitute each column of K.
    let mut a = kd;
    // Forward substitution on columns: solve L·X = K (column-wise).
    for col in 0..f {
        for i in 0..f {
            let mut acc = a[i * f + col];
            for j in 0..i {
                acc = (-chol.l(i, j)).mul_add(a[j * f + col], acc);
            }
            a[i * f + col] = acc / chol.l(i, i);
        }
    }
    // X·L⁻ᵀ = (L⁻¹Xᵀ)ᵀ: forward substitution on rows.
    for row in 0..f {
        for i in 0..f {
            let mut acc = a[row * f + i];
            for j in 0..i {
                acc = (-chol.l(i, j)).mul_add(a[row * f + j], acc);
            }
            a[row * f + i] = acc / chol.l(i, i);
        }
    }
    let (vals, vecs) = fs_la::eigen::jacobi_eigh(&a, f);
    // Lowest `count`: ascending already. Back-transform y → φ = L⁻ᵀy.
    let mut lambdas = Vec::with_capacity(count);
    let mut phis = Vec::with_capacity(count);
    for k in 0..count.min(f) {
        lambdas.push(vals[k]);
        let mut y: Vec<f64> = (0..f).map(|r| vecs[r * f + k]).collect();
        // Solve Lᵀφ = y (back substitution).
        for i in (0..f).rev() {
            let mut acc = y[i];
            for j in i + 1..f {
                acc = (-chol.l(j, i)).mul_add(y[j], acc);
            }
            y[i] = acc / chol.l(i, i);
        }
        phis.push(y);
    }
    (lambdas, phis, free_idx)
}

/// Exact eigenvalue design gradient for eigenpair (λ, φ) — w.r.t. the
/// RAW design through the full pipeline chain:
/// dλ/dρ̄_c = φᵀ(E′·K_c − λ·m′·M_c)φ (φ M-normalized), then the
/// projection/filter pullback. `phi_full` is the free-dof vector
/// scattered to full length.
#[must_use]
pub fn eigenvalue_gradient(
    pipeline: &DesignPipeline,
    elasticity: &DensityElasticity,
    rho: &[f64],
    lambda: f64,
    phi_full: &[f64],
) -> Vec<f64> {
    let p = &pipeline.params;
    let (rho_tilde, rho_bar, _) = pipeline.forward(rho);
    let strain = elasticity.cell_energies(phi_full);
    let kinetic = elasticity.cell_kinetic(phi_full);
    // dλ/dρ̄ per cell.
    let dlam_drhobar: Vec<f64> = rho_bar
        .iter()
        .zip(strain.iter().zip(&kinetic))
        .map(|(&rb, (&s, &m))| {
            let rc = rb.clamp(0.0, 1.0);
            let dsimp = (1.0 - p.e_min) * p.penal * fs_math::det::pow(rc.max(1e-12), p.penal - 1.0);
            let dmass = mass_interp_derivative(rb);
            dsimp.mul_add(s, -(lambda * dmass * m))
        })
        .collect();
    // Chain through projection then the transposed filter.
    let chained: Vec<f64> = rho_tilde
        .iter()
        .zip(&dlam_drhobar)
        .map(|(&rt, &d)| d * crate::filter::heaviside_derivative(rt, p.beta, p.eta))
        .collect();
    pipeline.filter.apply_transpose(&chained)
}

/// Smooth-min aggregation of the lowest eigenvalues:
/// λ_agg = −(1/β)·ln Σ exp(−β·λ_k) ≤ min_k λ_k, with weights
/// w_k = exp(−βλ_k)/Σ — the exact aggregate gradient is Σ w_k·∇λ_k.
/// Returns (λ_agg, weights).
#[must_use]
pub fn smooth_min(lambdas: &[f64], beta: f64) -> (f64, Vec<f64>) {
    let lmin = lambdas.iter().copied().fold(f64::INFINITY, f64::min);
    let exps: Vec<f64> = lambdas
        .iter()
        .map(|l| fs_math::det::exp(-beta * (l - lmin)))
        .collect();
    let sum: f64 = exps.iter().sum();
    let agg = lmin - fs_math::det::ln(sum) / beta;
    let weights: Vec<f64> = exps.iter().map(|e| e / sum).collect();
    (agg, weights)
}

/// Scatter a free-dof vector to full dof length (zeros on fixed).
#[must_use]
pub fn scatter_full(phi_free: &[f64], free_idx: &[usize], n: usize) -> Vec<f64> {
    let mut full = vec![0.0f64; n];
    for (v, &d) in phi_free.iter().zip(free_idx) {
        full[d] = *v;
    }
    full
}

/// Full smooth-min eigenfrequency objective and its exact design
/// gradient through the whole chain: returns (λ_agg, dλ_agg/dρ).
pub fn eigenfrequency_objective(
    pipeline: &DesignPipeline,
    elasticity: &mut DensityElasticity,
    rho: &[f64],
    cluster: usize,
    agg_beta: f64,
) -> (f64, Vec<f64>) {
    let (_, rho_bar, moduli) = pipeline.forward(rho);
    elasticity.moduli = moduli;
    let mass: Vec<f64> = rho_bar.iter().map(|&r| mass_interp(r)).collect();
    let (lambdas, phis, free_idx) = lowest_eigenpairs(elasticity, &mass, cluster);
    let (agg, weights) = smooth_min(&lambdas, agg_beta);
    let n = elasticity.n();
    let mut grad = vec![0.0f64; rho.len()];
    for ((lam, phi), w) in lambdas.iter().zip(&phis).zip(&weights) {
        let phi_full = scatter_full(phi, &free_idx, n);
        let g = eigenvalue_gradient(pipeline, elasticity, rho, *lam, &phi_full);
        for (gi, gk) in grad.iter_mut().zip(&g) {
            *gi = w.mul_add(*gk, *gi);
        }
    }
    (agg, grad)
}
