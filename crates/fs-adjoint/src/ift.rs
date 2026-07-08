//! Matrix-free IFT adjoints: for R(u; p) = 0 solved to tolerance,
//! dJ/dp = ∂J/∂p − (∂R/∂p)ᵀ·λ with (∂R/∂u)ᵀ·λ = ∂J/∂u — ONE
//! transposed solve through fs-solver, never a differentiated Krylov
//! iteration (fs-ad's IFT contract, upgraded from its dense v1 to the
//! matrix-free stack). Includes the density-parameterized Poisson
//! problem (the SIMP volumetric chain rule) as the canonical
//! parameterized-operator fixture: K(ρ) = Σ_t ρ_t·K_t with exact
//! per-cell derivative contributions.

use fs_feec::ElementGeometry;
use fs_rep_mesh::TetComplex;
use fs_solver::{GmresState, LinearOp};

/// Evidence attached to a matrix-free IFT gradient.
#[derive(Debug, Clone)]
pub struct AdjointReport {
    /// Adjoint-solve iterations.
    pub iters: usize,
    /// Adjoint relative residual achieved.
    pub adjoint_residual: f64,
    /// Whether the adjoint solve converged.
    pub converged: bool,
}

/// dJ/dp at a solution of R(u; p) = 0, matrix-free.
///
/// - `jacobian`: ∂R/∂u at (u*, p) as an adjoint-equipped operator;
/// - `djdu`: ∂J/∂u (the adjoint right-hand side);
/// - `djdp`: ∂J/∂p (the explicit part; empty slice for none);
/// - `drdp_t`: λ ↦ (∂R/∂p)ᵀ·λ (the parameter-space pullback).
///
/// Returns (gradient, report). The adjoint system is solved with
/// transposed GMRES sharing the operator's infrastructure.
pub fn ift_gradient_matfree<A: LinearOp>(
    jacobian: &A,
    djdu: &[f64],
    djdp: &[f64],
    drdp_t: &dyn Fn(&[f64]) -> Vec<f64>,
    tol: f64,
    max_cycles: usize,
) -> (Vec<f64>, AdjointReport) {
    let mut st = GmresState::new(djdu, 30);
    let rep = st.run(jacobian, djdu, tol, max_cycles, true);
    let pullback = drdp_t(&st.x);
    let mut g = if djdp.is_empty() {
        vec![0.0f64; pullback.len()]
    } else {
        djdp.to_vec()
    };
    for (gi, pi) in g.iter_mut().zip(&pullback) {
        *gi -= pi;
    }
    (
        g,
        AdjointReport {
            iters: rep.iters,
            adjoint_residual: rep.rel_residual,
            converged: rep.converged,
        },
    )
}

/// Density-parameterized Poisson on a tet complex (the SIMP chain-rule
/// fixture, and a real operator family): K(ρ)·u = b with
/// K(ρ) = Σ_t ρ_t·|V_t|·G_t (G_t the P1 gradient Gram — the exact
/// per-cell stiffness), homogeneous Dirichlet on the unit-cube
/// boundary. Every derivative below is EXACT (the assembly is linear
/// in ρ).
pub struct DensityPoisson<'c> {
    complex: &'c TetComplex,
    geo: ElementGeometry,
    interior: Vec<usize>,
    slot: Vec<usize>,
    /// Current densities (one per tet, > 0).
    pub rho: Vec<f64>,
}

impl<'c> DensityPoisson<'c> {
    /// Build for a complex with unit-cube positions.
    #[must_use]
    pub fn new(complex: &'c TetComplex, positions: &[[f64; 3]], rho: Vec<f64>) -> Self {
        assert_eq!(rho.len(), complex.tets.len(), "one density per tet");
        let geo = fs_feec::element_geometry(complex, positions);
        let interior: Vec<usize> = (0..positions.len())
            .filter(|&v| !fs_feec::on_unit_cube_boundary(positions[v]))
            .collect();
        let mut slot = vec![usize::MAX; positions.len()];
        for (i, &v) in interior.iter().enumerate() {
            slot[v] = i;
        }
        DensityPoisson {
            complex,
            geo,
            interior,
            slot,
            rho,
        }
    }

    /// Interior dof count.
    #[must_use]
    pub fn n(&self) -> usize {
        self.interior.len()
    }

    /// y = K(ρ)·x on interior dofs (matrix-free per-cell apply).
    fn apply_k(&self, rho: &[f64], x: &[f64], y: &mut [f64]) {
        y.fill(0.0);
        for (t, tet) in self.complex.tets.iter().enumerate() {
            let w = rho[t] * self.geo.vol_signed[t].abs();
            let gram = &self.geo.gram[t];
            for (a, &va) in tet.iter().enumerate() {
                let ia = self.slot[va as usize];
                if ia == usize::MAX {
                    continue;
                }
                let mut acc = 0.0f64;
                for (bb, &vb) in tet.iter().enumerate() {
                    let ib = self.slot[vb as usize];
                    if ib != usize::MAX {
                        acc = gram[a][bb].mul_add(x[ib], acc);
                    }
                }
                y[ia] = (w * acc).mul_add(1.0, y[ia]);
            }
        }
    }

    /// Apply K(density)·x for an ARBITRARY (possibly signed) density
    /// vector — the operator family is linear in ρ, so this is the
    /// directional operator K(v) the second-order adjoints need.
    #[must_use]
    pub fn apply_density(&self, density: &[f64], x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0f64; self.n()];
        self.apply_k(density, x, &mut y);
        y
    }

    /// Per-cell chain rule: out[t] = λᵀ·K_t·u (the exact
    /// ∂(K(ρ)u)/∂ρ_t pullback — the SIMP volumetric derivative).
    #[must_use]
    pub fn density_pullback(&self, lambda: &[f64], u: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0f64; self.complex.tets.len()];
        for (t, tet) in self.complex.tets.iter().enumerate() {
            let w = self.geo.vol_signed[t].abs();
            let gram = &self.geo.gram[t];
            let mut acc = 0.0f64;
            for (a, &va) in tet.iter().enumerate() {
                let ia = self.slot[va as usize];
                if ia == usize::MAX {
                    continue;
                }
                for (bb, &vb) in tet.iter().enumerate() {
                    let ib = self.slot[vb as usize];
                    if ib != usize::MAX {
                        acc += lambda[ia] * gram[a][bb] * u[ib];
                    }
                }
            }
            out[t] = w * acc;
        }
        out
    }
}

/// The K(ρ) operator view (symmetric — its transpose is itself).
pub struct DensityOp<'a, 'c> {
    problem: &'a DensityPoisson<'c>,
}

impl<'a, 'c> DensityOp<'a, 'c> {
    /// Borrow the operator at the problem's current densities.
    #[must_use]
    pub fn new(problem: &'a DensityPoisson<'c>) -> Self {
        DensityOp { problem }
    }
}

impl LinearOp for DensityOp<'_, '_> {
    fn n(&self) -> usize {
        self.problem.n()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.problem.apply_k(&self.problem.rho, x, y);
    }
}
