//! The density hygiene stack: Helmholtz PDE filtering (one solver,
//! two jobs — this IS the Sobolev/Riesz solve from fs-adjoint, applied
//! to densities) and Heaviside projection with β continuation. Both
//! stages carry EXACT derivatives: the filter is linear (its chain
//! rule is the TRANSPOSED map — scatterᵀ ∘ solveᵀ ∘ gatherᵀ, and the
//! Helmholtz operator is symmetric so solveᵀ = solve), the projection
//! is pointwise with a closed-form slope.

use fs_rep_mesh::TetComplex;
use fs_sparse::Csr;

/// Helmholtz density filter on a tet complex: cell densities are
/// volume-scattered to vertices, smoothed by (M + r²K)⁻¹M, and
/// volume-averaged back to cells. Linear, symmetric in the
/// volume-weighted inner product, mesh-independent length scale r.
pub struct DensityFilter {
    mass: Csr,
    helmholtz: Csr,
    /// Filter radius r (the length scale).
    pub radius: f64,
    /// Scatter matrix rows: for each cell, its 4 vertices with weight
    /// |V_c|/4 (volume-weighted vertex averaging).
    cells: Vec<[u32; 4]>,
    /// Cell volumes.
    vol: Vec<f64>,
    /// Vertex lumped volumes (Σ incident |V|/4).
    vertex_vol: Vec<f64>,
    nv: usize,
}

impl DensityFilter {
    /// Build for a complex (full vertex space — no Dirichlet here;
    /// natural boundary conditions are the CORRECT filter behavior:
    /// densities near the boundary are not pulled to zero).
    #[must_use]
    pub fn new(complex: &TetComplex, positions: &[[f64; 3]], radius: f64) -> DensityFilter {
        assert!(
            radius.is_finite() && radius >= 0.0,
            "density-filter radius must be finite and nonnegative"
        );
        let radius_squared = radius * radius;
        assert!(
            radius_squared.is_finite(),
            "density-filter squared radius must remain finite"
        );
        let geo = fs_feec::element_geometry(complex, positions);
        let mass = fs_feec::mass_matrix(complex, &geo, 0);
        let stiff = fs_feec::stiffness(
            &fs_feec::incidence_to_csr(&complex.d0()),
            &fs_feec::mass_matrix(complex, &geo, 1),
        );
        let nv = complex.vertex_count;
        let mut vertex_vol = vec![0.0f64; nv];
        let vol: Vec<f64> = geo.vol_signed.iter().map(|v| v.abs()).collect();
        assert!(
            vol.iter().all(|volume| volume.is_finite() && *volume > 0.0),
            "density filter requires finite positive cell volumes"
        );
        for (t, tet) in complex.tets.iter().enumerate() {
            for &v in tet {
                vertex_vol[v as usize] += vol[t] / 4.0;
            }
        }
        assert!(
            vertex_vol
                .iter()
                .all(|volume| volume.is_finite() && *volume > 0.0),
            "density filter requires every vertex to have finite positive incident volume"
        );
        // Assemble (M + r²K) ONCE.
        let mut coo = fs_sparse::Coo::new(nv, nv);
        for r in 0..nv {
            let (cols, vals) = mass.row(r);
            for (&c, &v) in cols.iter().zip(vals) {
                coo.push(r, c, v);
            }
            let (cols, vals) = stiff.row(r);
            for (&c, &v) in cols.iter().zip(vals) {
                coo.push(r, c, radius_squared * v);
            }
        }
        DensityFilter {
            mass,
            helmholtz: coo.assemble(),
            radius,
            cells: complex.tets.clone(),
            vol,
            vertex_vol,
            nv,
        }
    }

    /// Solve (M + r²K)·x = rhs on the FULL vertex space (natural BCs).
    fn helmholtz_solve(&self, rhs: &[f64]) -> Vec<f64> {
        let a = fs_solver::CsrOp::symmetric(self.helmholtz.clone());
        let mut st = fs_solver::CgState::new(&a, &fs_sparse::precond::IdentityPrecond, rhs);
        let rep = st.run(&a, &fs_sparse::precond::IdentityPrecond, 1e-12, 20_000);
        assert!(rep.converged, "Helmholtz filter solve failed: {rep:?}");
        st.x
    }

    /// FORWARD filter: cell densities → filtered cell densities.
    #[must_use]
    pub fn apply(&self, rho: &[f64]) -> Vec<f64> {
        assert_eq!(
            rho.len(),
            self.cells.len(),
            "density-filter input must contain one finite value per cell"
        );
        assert!(
            rho.iter().all(|value| value.is_finite()),
            "density-filter input must contain only finite values"
        );
        // Scatter: vertex value = Σ_c∋v (|V_c|/4)·ρ_c / vertex_vol.
        let mut vtx = vec![0.0f64; self.nv];
        for (c, tet) in self.cells.iter().enumerate() {
            for &v in tet {
                vtx[v as usize] += self.vol[c] / 4.0 * rho[c];
            }
        }
        for (x, w) in vtx.iter_mut().zip(&self.vertex_vol) {
            *x /= w;
        }
        let mut rhs = vec![0.0f64; self.nv];
        self.mass.spmv(&vtx, &mut rhs);
        let smooth = self.helmholtz_solve(&rhs);
        // Gather: cell value = vertex average.
        self.cells
            .iter()
            .map(|tet| tet.iter().map(|&v| smooth[v as usize]).sum::<f64>() / 4.0)
            .collect()
    }

    /// TRANSPOSED filter (the chain-rule pullback of `apply`): maps a
    /// sensitivity w.r.t. FILTERED densities back to raw densities.
    /// gatherᵀ (cell → vertices /4), solveᵀ = solve of the transposed
    /// composition M·(M + r²K)⁻¹ (both symmetric — but the ORDER
    /// matters: forward is solve(M·scatter(ρ)), so the pullback is
    /// scatterᵀ(M·solve(gatherᵀ(g)))).
    #[must_use]
    pub fn apply_transpose(&self, g_filtered: &[f64]) -> Vec<f64> {
        assert_eq!(
            g_filtered.len(),
            self.cells.len(),
            "density-filter pullback must contain one finite value per cell"
        );
        assert!(
            g_filtered.iter().all(|value| value.is_finite()),
            "density-filter pullback must contain only finite values"
        );
        // gatherᵀ: vertex accumulation of cell sensitivities /4.
        let mut vtx = vec![0.0f64; self.nv];
        for (c, tet) in self.cells.iter().enumerate() {
            for &v in tet {
                vtx[v as usize] += g_filtered[c] / 4.0;
            }
        }
        // Transpose of x ↦ (M + r²K)⁻¹·M·x is M·(M + r²K)⁻¹ (both
        // factors symmetric): solve first, THEN mass-multiply.
        let sol = self.helmholtz_solve(&vtx);
        let mut mx = vec![0.0f64; self.nv];
        self.mass.spmv(&sol, &mut mx);
        // scatterᵀ: back to cells with the forward scatter weights.
        self.cells
            .iter()
            .enumerate()
            .map(|(c, tet)| {
                tet.iter()
                    .map(|&v| self.vol[c] / 4.0 / self.vertex_vol[v as usize] * mx[v as usize])
                    .sum::<f64>()
            })
            .collect()
    }
}

/// Heaviside projection with threshold η and sharpness β:
/// ρ̄ = (tanh(βη) + tanh(β(ρ̃−η))) / (tanh(βη) + tanh(β(1−η))).
/// Monotone in ρ̃, exact 0→0 and 1→1, → step function as β → ∞.
#[must_use]
pub fn heaviside(rho_tilde: f64, beta: f64, eta: f64) -> f64 {
    assert_valid_projection_inputs(rho_tilde, beta, eta);
    // The analytic beta -> 0 limit is the identity map. At and below one
    // machine epsilon, evaluating the quotient directly can round both tanh
    // terms to zero even though the limit is well defined.
    if beta <= f64::EPSILON {
        return rho_tilde;
    }
    let denom = tanh(beta * eta) + tanh(beta * (1.0 - eta));
    (tanh(beta * eta) + tanh(beta * (rho_tilde - eta))) / denom
}

/// dρ̄/dρ̃ — the closed-form projection slope.
#[must_use]
pub fn heaviside_derivative(rho_tilde: f64, beta: f64, eta: f64) -> f64 {
    assert_valid_projection_inputs(rho_tilde, beta, eta);
    if beta <= f64::EPSILON {
        return 1.0;
    }
    let denom = tanh(beta * eta) + tanh(beta * (1.0 - eta));
    let t = tanh(beta * (rho_tilde - eta));
    beta * (1.0 - t * t) / denom
}

fn assert_valid_projection_inputs(rho_tilde: f64, beta: f64, eta: f64) {
    assert!(rho_tilde.is_finite(), "filtered density must be finite");
    assert!(
        beta.is_finite() && beta >= 0.0,
        "Heaviside sharpness beta must be finite and nonnegative"
    );
    assert!(
        eta.is_finite() && (0.0..=1.0).contains(&eta),
        "Heaviside threshold eta must be finite and lie in [0, 1]"
    );
}

/// tanh through the strict exp kernel (no platform libm in the
/// pipeline — golden-hash discipline).
fn tanh(x: f64) -> f64 {
    if x > 20.0 {
        return 1.0;
    }
    if x < -20.0 {
        return -1.0;
    }
    let e2 = fs_math::det::exp(2.0 * x);
    (e2 - 1.0) / (e2 + 1.0)
}
