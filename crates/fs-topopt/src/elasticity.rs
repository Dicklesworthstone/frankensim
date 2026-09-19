//! Density-parameterized vector elasticity: K(ρ̄) = Σ_c E(ρ̄_c)·K_c
//! with per-cell UNIT-modulus stiffness blocks K_c = |V_c|·B_aᵀCB_b
//! (fs-material tangent at zero strain, fs-feec barycentric
//! gradients) kept SEPARATE so the SIMP chain rule is exact:
//! ∂(Ku)/∂ρ̄_c = E′(ρ̄_c)·K_c·u — the fs-adjoint `DensityPoisson`
//! pattern lifted to 3 dofs per vertex. Compliance is self-adjoint
//! (λ = u), so sensitivities cost ZERO extra solves — stated, used,
//! and FD-verified in the battery.

use fs_material::{IsotropicElastic, SmallStrainLaw};
use fs_rep_mesh::TetComplex;
use fs_solver::LinearOp;

/// The density-elasticity problem on a tet complex: per-cell 12×12
/// unit-modulus stiffness blocks + Dirichlet mask on vector dofs.
pub struct DensityElasticity {
    /// Per-cell 12×12 unit-modulus element stiffness (row-major).
    ke: Vec<[f64; 144]>,
    /// Per-cell 12×12 unit-density consistent mass (row-major):
    /// the P1 tet pattern (|V|/20)·(1 + δ_ab) per displacement
    /// component — kept separate like `ke` for exact chain rules.
    me: Vec<[f64; 144]>,
    /// Per-cell C·B products (6×12 row-major): unit-modulus stress
    /// recovery σ = CB·u_local (constant per P1 tet).
    cb: Vec<[f64; 72]>,
    /// Cell → 4 vertex ids.
    tets: Vec<[u32; 4]>,
    /// Vector-dof count (3 per vertex).
    n: usize,
    /// Free-dof mask (false = Dirichlet-fixed).
    free: Vec<bool>,
    /// Current SIMP moduli E(ρ̄_c) (set per apply).
    pub moduli: Vec<f64>,
}

impl DensityElasticity {
    /// Build from a complex, positions, material, and a Dirichlet
    /// predicate over vertex positions (all 3 components fixed).
    ///
    /// # Panics
    /// On invalid material parameters.
    #[must_use]
    pub fn new(
        complex: &TetComplex,
        positions: &[[f64; 3]],
        youngs: f64,
        poisson: f64,
        fixed: &dyn Fn([f64; 3]) -> bool,
    ) -> DensityElasticity {
        let law = IsotropicElastic::new(youngs, poisson, 1.0).expect("valid material");
        let c = law.tangent(&[0.0; 6], &());
        let geo = fs_feec::element_geometry(complex, positions);
        let mut ke = Vec::with_capacity(complex.tets.len());
        let mut cb = Vec::with_capacity(complex.tets.len());
        for (t, _tet) in complex.tets.iter().enumerate() {
            let vol = geo.vol_signed[t].abs();
            // Bᵀ rows per node (3×6), from ∇λ_a.
            let bt = |a: usize| -> [[f64; 6]; 3] {
                let g = geo.grads[t][a];
                [
                    [g[0], 0.0, 0.0, g[1], 0.0, g[2]],
                    [0.0, g[1], 0.0, g[0], g[2], 0.0],
                    [0.0, 0.0, g[2], 0.0, g[1], g[0]],
                ]
            };
            let mut k = [0.0f64; 144];
            for a in 0..4 {
                let bta = bt(a);
                for b in 0..4 {
                    let btb = bt(b);
                    for (i, bai) in bta.iter().enumerate() {
                        for (j, bbj) in btb.iter().enumerate() {
                            let mut acc = 0.0f64;
                            for (p, baip) in bai.iter().enumerate() {
                                for (q, bbjq) in bbj.iter().enumerate() {
                                    acc = (baip * c[p][q]).mul_add(*bbjq, acc);
                                }
                            }
                            k[(3 * a + i) * 12 + (3 * b + j)] = vol * acc;
                        }
                    }
                }
            }
            ke.push(k);
            // CB (6×12): stress = C·(B·u); B rows are the bt blocks
            // transposed.
            let mut cbm = [0.0f64; 72];
            for bnode in 0..4 {
                let btb = bt(bnode);
                for si in 0..6 {
                    // σ_si column (3b+comp) = Σ_q C[si][q]·B[q][3b+comp],
                    // and B row q at that column is btb[comp][q].
                    for (comp, btbc) in btb.iter().enumerate() {
                        let mut v = 0.0f64;
                        for q in 0..6 {
                            v = c[si][q].mul_add(btbc[q], v);
                        }
                        cbm[si * 12 + 3 * bnode + comp] = v;
                    }
                }
            }
            cb.push(cbm);
        }
        // Consistent unit-density mass blocks.
        let mut me = Vec::with_capacity(complex.tets.len());
        for (t, _tet) in complex.tets.iter().enumerate() {
            let vol = geo.vol_signed[t].abs();
            let mut m = [0.0f64; 144];
            for a in 0..4 {
                for bb in 0..4 {
                    let w = vol / 20.0 * if a == bb { 2.0 } else { 1.0 };
                    for comp in 0..3 {
                        m[(3 * a + comp) * 12 + (3 * bb + comp)] = w;
                    }
                }
            }
            me.push(m);
        }
        let n = 3 * complex.vertex_count;
        let mut free = vec![true; n];
        for (v, &p) in positions.iter().enumerate() {
            if fixed(p) {
                for comp in 0..3 {
                    free[3 * v + comp] = false;
                }
            }
        }
        DensityElasticity {
            ke,
            me,
            cb,
            tets: complex.tets.clone(),
            n,
            free,
            moduli: vec![1.0; complex.tets.len()],
        }
    }

    /// Vector-dof count.
    #[must_use]
    pub fn n(&self) -> usize {
        self.n
    }

    /// Cell count.
    #[must_use]
    pub fn cells(&self) -> usize {
        self.ke.len()
    }

    /// The free-dof mask.
    #[must_use]
    pub fn free(&self) -> &[bool] {
        &self.free
    }

    /// Assemble DENSE reduced (K(moduli), M(densities)) on the free
    /// dofs (fixture-scale eigenproblems; row-major f×f each).
    #[must_use]
    pub fn assemble_dense(&self, densities: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<usize>) {
        let free_idx: Vec<usize> = (0..self.n).filter(|&d| self.free[d]).collect();
        let mut slot = vec![usize::MAX; self.n];
        for (i, &d) in free_idx.iter().enumerate() {
            slot[d] = i;
        }
        let f = free_idx.len();
        let mut kd = vec![0.0f64; f * f];
        let mut md = vec![0.0f64; f * f];
        for (c, tet) in self.tets.iter().enumerate() {
            let (e, rho) = (self.moduli[c], densities[c]);
            for a in 0..4 {
                for comp_a in 0..3 {
                    let da = 3 * tet[a] as usize + comp_a;
                    if slot[da] == usize::MAX {
                        continue;
                    }
                    for (bb, &vb) in tet.iter().enumerate() {
                        for comp_b in 0..3 {
                            let db = 3 * vb as usize + comp_b;
                            if slot[db] == usize::MAX {
                                continue;
                            }
                            let row = 3 * a + comp_a;
                            let col = 3 * bb + comp_b;
                            kd[slot[da] * f + slot[db]] += e * self.ke[c][row * 12 + col];
                            md[slot[da] * f + slot[db]] += rho * self.me[c][row * 12 + col];
                        }
                    }
                }
            }
        }
        (kd, md, free_idx)
    }

    /// Per-cell UNIT-MODULUS stress vectors σ = CB·u_local (6 Voigt
    /// components per cell; multiply by the stress interpolation
    /// outside) and the von Mises scalar per cell.
    #[must_use]
    pub fn cell_stress(&self, u: &[f64]) -> Vec<[f64; 6]> {
        let mut out = Vec::with_capacity(self.cb.len());
        for (cbm, tet) in self.cb.iter().zip(&self.tets) {
            let mut ul = [0.0f64; 12];
            for (a, &v) in tet.iter().enumerate() {
                for comp in 0..3 {
                    let d = 3 * v as usize + comp;
                    ul[3 * a + comp] = if self.free[d] { u[d] } else { 0.0 };
                }
            }
            let mut sig = [0.0f64; 6];
            for (si, sv) in sig.iter_mut().enumerate() {
                let mut acc = 0.0f64;
                for (j, uj) in ul.iter().enumerate() {
                    acc = cbm[si * 12 + j].mul_add(*uj, acc);
                }
                *sv = acc;
            }
            out.push(sig);
        }
        out
    }

    /// Pullback of per-cell stress sensitivities to the dof vector:
    /// given dJ/dσ per cell (6 Voigt each), returns Σ_c (CB_c)ᵀ·dJ/dσ_c
    /// scattered to free dofs — the RHS of the stress adjoint solve.
    #[must_use]
    pub fn stress_pullback(&self, dj_dsigma: &[[f64; 6]]) -> Vec<f64> {
        let mut rhs = vec![0.0f64; self.n];
        for ((cbm, tet), ds) in self.cb.iter().zip(&self.tets).zip(dj_dsigma) {
            for (a, &v) in tet.iter().enumerate() {
                for comp in 0..3 {
                    let d = 3 * v as usize + comp;
                    if !self.free[d] {
                        continue;
                    }
                    let col = 3 * a + comp;
                    let mut acc = 0.0f64;
                    for (si, dsi) in ds.iter().enumerate() {
                        acc = cbm[si * 12 + col].mul_add(*dsi, acc);
                    }
                    rhs[d] += acc;
                }
            }
        }
        rhs
    }

    /// Per-cell KINETIC quadratic form m_c = uᵀ·M_c·u (unit density —
    /// multiply by the mass-interpolation slope outside).
    #[must_use]
    pub fn cell_kinetic(&self, u: &[f64]) -> Vec<f64> {
        self.cell_quadratic(&self.me, u)
    }

    /// Per-cell strain energy density contribution: e_c = uᵀ·K_c·u
    /// (UNIT modulus — multiply by E′ outside for the chain rule).
    #[must_use]
    pub fn cell_energies(&self, u: &[f64]) -> Vec<f64> {
        self.cell_quadratic(&self.ke, u)
    }

    fn cell_quadratic(&self, blocks: &[[f64; 144]], u: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0f64; blocks.len()];
        for (c, (k, tet)) in blocks.iter().zip(&self.tets).enumerate() {
            let mut ul = [0.0f64; 12];
            for (a, &v) in tet.iter().enumerate() {
                for comp in 0..3 {
                    let d = 3 * v as usize + comp;
                    ul[3 * a + comp] = if self.free[d] { u[d] } else { 0.0 };
                }
            }
            let mut acc = 0.0f64;
            for i in 0..12 {
                for j in 0..12 {
                    acc += ul[i] * k[i * 12 + j] * ul[j];
                }
            }
            out[c] = acc;
        }
        out
    }
}

impl LinearOp for DensityElasticity {
    fn n(&self) -> usize {
        self.n
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y.fill(0.0);
        for (c, (k, tet)) in self.ke.iter().zip(&self.tets).enumerate() {
            let e = self.moduli[c];
            let mut xl = [0.0f64; 12];
            for (a, &v) in tet.iter().enumerate() {
                for comp in 0..3 {
                    let d = 3 * v as usize + comp;
                    xl[3 * a + comp] = if self.free[d] { x[d] } else { 0.0 };
                }
            }
            for (a, &v) in tet.iter().enumerate() {
                for comp in 0..3 {
                    let d = 3 * v as usize + comp;
                    if !self.free[d] {
                        continue;
                    }
                    let row = 3 * a + comp;
                    let mut acc = 0.0f64;
                    for (j, xlj) in xl.iter().enumerate() {
                        acc = k[row * 12 + j].mul_add(*xlj, acc);
                    }
                    y[d] = (e * acc).mul_add(1.0, y[d]);
                }
            }
        }
        // Identity on fixed dofs (keeps the operator SPD on the full
        // vector space).
        for (i, yi) in y.iter_mut().enumerate() {
            if !self.free[i] {
                *yi = x[i];
            }
        }
    }
}

impl DensityElasticity {
    /// Apply the consistent mass on the free displacement space without a
    /// global matrix. Fixed inputs are ignored and fixed outputs are ZERO:
    /// adding the stiffness operator's identity here would introduce false
    /// eigenmodes at constrained nodes. `densities` are physical mass weights,
    /// not stiffness moduli and not a lumped-mass approximation.
    ///
    /// # Panics
    /// On incompatible shapes or non-finite/negative mass weights.
    pub fn apply_mass(&self, densities: &[f64], x: &[f64], y: &mut [f64]) {
        self.assert_mass_weights(densities);
        assert_eq!(x.len(), self.n, "mass input shape mismatch");
        assert_eq!(y.len(), self.n, "mass output shape mismatch");
        y.fill(0.0);
        for ((mass, tet), &density) in self.me.iter().zip(&self.tets).zip(densities) {
            let mut local = [0.0; 12];
            for (a, &vertex) in tet.iter().enumerate() {
                for component in 0..3 {
                    let dof = 3 * vertex as usize + component;
                    if self.free[dof] {
                        local[3 * a + component] = x[dof];
                    }
                }
            }
            for (a, &vertex) in tet.iter().enumerate() {
                for component in 0..3 {
                    let dof = 3 * vertex as usize + component;
                    if self.free[dof] {
                        let row = 3 * a + component;
                        let mut value = 0.0;
                        for (column, &input) in local.iter().enumerate() {
                            value = mass[12 * row + column].mul_add(input, value);
                        }
                        y[dof] = density.mul_add(value, y[dof]);
                    }
                }
            }
        }
    }

    /// Diagonal of the CURRENT stiffness operator, for matrix-free Jacobi
    /// preconditioning. Fixed entries are one, matching `LinearOp::apply`.
    #[must_use]
    pub fn stiffness_diagonal(&self) -> Vec<f64> {
        assert_eq!(self.moduli.len(), self.cells(), "stiffness weight shape mismatch");
        assert!(self.moduli.iter().all(|e| e.is_finite() && *e > 0.0),
            "stiffness weights must be finite and positive");
        let mut diagonal = self.block_diagonal(&self.ke, &self.moduli);
        for (value, &free) in diagonal.iter_mut().zip(&self.free) {
            if !free {
                *value = 1.0;
            }
        }
        diagonal
    }

    /// Diagonal of the consistent mass. Fixed entries are zero. This is a
    /// diagnostic/preconditioning quantity, NOT a replacement mass operator.
    #[must_use]
    pub fn mass_diagonal(&self, densities: &[f64]) -> Vec<f64> {
        self.assert_mass_weights(densities);
        self.block_diagonal(&self.me, densities)
    }

    fn assert_mass_weights(&self, densities: &[f64]) {
        assert_eq!(densities.len(), self.cells(), "one mass weight per cell is required");
        assert!(densities.iter().all(|r| r.is_finite() && *r >= 0.0),
            "mass weights must be finite and nonnegative");
    }

    fn block_diagonal(&self, blocks: &[[f64; 144]], weights: &[f64]) -> Vec<f64> {
        let mut diagonal = vec![0.0; self.n];
        for ((block, tet), &weight) in blocks.iter().zip(&self.tets).zip(weights) {
            for (a, &vertex) in tet.iter().enumerate() {
                for component in 0..3 {
                    let dof = 3 * vertex as usize + component;
                    if self.free[dof] {
                        let local = 3 * a + component;
                        diagonal[dof] = weight.mul_add(block[13 * local], diagonal[dof]);
                    }
                }
            }
        }
        diagonal
    }
}

#[cfg(test)]
mod matrix_free_mass_tests {
    use super::*;

    #[test]
    fn consistent_mass_and_diagonals_match_dense_assembly() {
        let (mesh, positions) = fs_feec::kuhn_cube(2);
        let mut elasticity = DensityElasticity::new(&mesh, &positions, 1.0, 0.3,
            &|p| p[0] == 0.0);
        let weights: Vec<f64> = (0..elasticity.cells())
            .map(|i| 0.2 + (i % 7) as f64 / 10.0).collect();
        elasticity.moduli = weights.iter().map(|r| r * r * r).collect();
        let (stiffness, mass, free) = elasticity.assemble_dense(&weights);
        let n = free.len();
        let input: Vec<f64> = (0..elasticity.n()).map(|i| (i % 11) as f64 - 5.0).collect();
        let mut actual = vec![f64::NAN; elasticity.n()];
        elasticity.apply_mass(&weights, &input, &mut actual);
        let md = elasticity.mass_diagonal(&weights);
        let kd = elasticity.stiffness_diagonal();
        for (row, &dof) in free.iter().enumerate() {
            let expected: f64 = free.iter().enumerate()
                .map(|(column, &j)| mass[row * n + column] * input[j]).sum();
            assert!((actual[dof] - expected).abs() < 1e-13);
            assert!((md[dof] - mass[row * n + row]).abs() < 1e-13);
            assert!((kd[dof] - stiffness[row * n + row]).abs() < 1e-13);
        }
        for (dof, &is_free) in elasticity.free().iter().enumerate() {
            if !is_free {
                assert_eq!(actual[dof], 0.0);
                assert_eq!(md[dof], 0.0);
                assert_eq!(kd[dof], 1.0);
            }
        }
        let energy: f64 = input.iter().zip(&actual).map(|(x, y)| x * y).sum();
        let cell_energy: f64 = elasticity.cell_kinetic(&input).iter().zip(&weights)
            .map(|(kinetic, weight)| kinetic * weight).sum();
        assert!((energy - cell_energy).abs() < 1e-12);
        assert!(energy > 0.0);
    }

    #[test]
    fn zero_mass_and_fixed_only_inputs_cannot_create_modes() {
        let (mesh, positions) = fs_feec::kuhn_cube(1);
        let elasticity = DensityElasticity::new(&mesh, &positions, 1.0, 0.3,
            &|p| p[0] == 0.0);
        let mut output = vec![17.0; elasticity.n()];
        elasticity.apply_mass(&vec![0.0; elasticity.cells()],
            &vec![1.0; elasticity.n()], &mut output);
        assert!(output.iter().all(|y| *y == 0.0));
        let fixed_only: Vec<f64> = elasticity.free().iter()
            .map(|free| if *free { 0.0 } else { 3.0 }).collect();
        elasticity.apply_mass(&vec![1.0; elasticity.cells()], &fixed_only, &mut output);
        assert!(output.iter().all(|y| *y == 0.0));
    }

    #[test]
    #[should_panic(expected = "mass weights must be finite and nonnegative")]
    fn invalid_mass_is_rejected() {
        let (mesh, positions) = fs_feec::kuhn_cube(1);
        let elasticity = DensityElasticity::new(&mesh, &positions, 1.0, 0.3, &|_| false);
        let mut weights = vec![1.0; elasticity.cells()];
        weights[0] = f64::NAN;
        let _ = elasticity.mass_diagonal(&weights);
    }
}
