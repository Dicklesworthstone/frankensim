//! The incompressible solver: BDM1–P0 saddle systems with interior-
//! penalty viscosity (full jumps — the normal component's jump is zero
//! by conformity, so the penalty acts tangentially exactly as the
//! H(div)-DG theory prescribes), fully weak Dirichlet velocity BCs,
//! upwinded DG convection on the normal-continuous advecting field
//! (w·n is single-valued on faces — H(div)'s gift to upwinding),
//! Picard iteration for steady NS, IMEX BDF1 for transients, and the
//! discrete adjoint of the linearized system. Dirichlet velocity BCs
//! split by component: u·n STRONG (boundary-edge BDM dofs are exactly
//! the normal moments — identity rows), tangential weak via SIP;
//! pressure pinned by replacing cell 0's continuity row (implied by
//! the others plus the prescribed net flux).

use crate::bdm::{CellBasis, cell_basis, eval_basis, tri_quad};
use crate::trimesh::TriMesh;
use fs_la::factor::lu;
use fs_solver::krylov::GmresState;
use fs_solver::op::CsrOp;
use fs_sparse::{Coo, Csr};

fn positive_viscosity(nu: f64) -> f64 {
    assert!(
        nu.is_finite() && nu > 0.0,
        "viscosity must be positive and finite"
    );
    nu
}

/// Unit-viscosity traction, retaining the original Laplacian path's order.
fn viscous_normal(grad: [[f64; 2]; 2], n: [f64; 2], symmetric: bool) -> [f64; 2] {
    if symmetric {
        let shear = grad[0][1] + grad[1][0];
        [
            2.0 * grad[0][0] * n[0] + shear * n[1],
            shear * n[0] + 2.0 * grad[1][1] * n[1],
        ]
    } else {
        [
            grad[0][0] * n[0] + grad[0][1] * n[1],
            grad[1][0] * n[0] + grad[1][1] * n[1],
        ]
    }
}

/// Solver parameters.
#[derive(Debug, Clone, Copy)]
pub struct FluxParams {
    /// Kinematic viscosity.
    pub nu: f64,
    /// SIP penalty constant (applied as σν/h per face).
    pub sigma: f64,
    /// Krylov tolerance.
    pub tol: f64,
    /// Krylov iteration cap.
    pub max_iters: usize,
}

impl Default for FluxParams {
    fn default() -> Self {
        FluxParams {
            nu: 1.0,
            sigma: 10.0,
            tol: 1e-10,
            max_iters: 40_000,
        }
    }
}

/// A discrete solution: velocity dofs (2 per edge) then pressures.
pub struct FluxSolution {
    /// Full dof vector.
    pub x: Vec<f64>,
    /// Velocity dof count (2 × edges).
    pub n_u: usize,
    /// GMRES iterations of the last solve.
    pub iters: usize,
    /// Final relative residual.
    pub rel_residual: f64,
}

/// The assembled operator pieces for one mesh.
pub struct FluxSystem<'m> {
    /// The mesh.
    pub mesh: &'m TriMesh,
    /// Per-cell bases (built once).
    pub bases: Vec<CellBasis>,
    /// Velocity dof count.
    pub n_u: usize,
    /// Total dof count (velocity + pressure).
    pub n: usize,
}

impl<'m> FluxSystem<'m> {
    /// Precompute the per-cell bases.
    #[must_use]
    pub fn new(mesh: &'m TriMesh) -> FluxSystem<'m> {
        let bases = (0..mesh.tris.len()).map(|t| cell_basis(mesh, t)).collect();
        let n_u = 2 * mesh.edges.len();
        FluxSystem {
            mesh,
            bases,
            n_u,
            n: n_u + mesh.tris.len(),
        }
    }

    /// Global velocity dof of (triangle-local dof i) on triangle `t`.
    fn dof(&self, t: usize, i: usize) -> usize {
        let (e, _) = self.mesh.tri_edges[t][i / 2];
        2 * e + (i % 2)
    }

    /// Assemble the Stokes operator and rhs for body force `f` and
    /// boundary velocity `g` (fully weak). `advect` optionally adds
    /// upwinded convection with the given BDM1 field (Picard).
    ///
    /// Returns (matrix, rhs).
    #[must_use]
    #[allow(clippy::too_many_lines)] // cell + face + boundary assembly, one narrative
    pub fn assemble(
        &self,
        params: FluxParams,
        f: &dyn Fn([f64; 2]) -> [f64; 2],
        g: &dyn Fn([f64; 2]) -> [f64; 2],
        advect: Option<&[f64]>,
    ) -> (Csr, Vec<f64>) {
        self.assemble_transient(params, f, g, advect, 0.0, None)
    }

    /// Assemble with an optional BDF1 mass term: `inv_dt`*M added to
    /// the velocity block and `inv_dt`*M*`uprev` to the rhs.
    #[must_use]
    #[allow(clippy::too_many_lines)] // cell + face + boundary assembly, one narrative
    pub fn assemble_transient(
        &self,
        params: FluxParams,
        f: &dyn Fn([f64; 2]) -> [f64; 2],
        g: &dyn Fn([f64; 2]) -> [f64; 2],
        advect: Option<&[f64]>,
        inv_dt: f64,
        uprev: Option<&[f64]>,
    ) -> (Csr, Vec<f64>) {
        self.assemble_operator(params, None, f, g, advect, inv_dt, uprev)
    }

    /// Assemble `-div(2 nu D(u)) + grad(p) = f`, with
    /// `D(u) = (grad(u) + grad(u)^T)/2`. The viscosity callback receives
    /// the cell index and physical quadrature point, so jumps between
    /// materials have unambiguous one-sided traces. It returns positive
    /// finite kinematic viscosity [m²/s] for a constant reference density;
    /// `params.nu` is unused. Temperature-dependent dynamic viscosity must
    /// be divided by that same reference density, not a local density.
    ///
    /// Both volume stresses and face tractions use the symmetric gradient.
    /// Weighted SIP uses harmonic viscosity/normal-cell-height penalties
    /// across interfaces. Viscosity is evaluated once per quadrature point,
    /// outside the basis-pair loops. The existing constant-viscosity
    /// vector-Laplacian entry point remains available through `assemble`.
    ///
    /// # Panics
    /// If viscosity or the SIP penalty is not positive and finite.
    #[must_use]
    pub fn assemble_variable_viscosity(
        &self,
        params: FluxParams,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        f: &dyn Fn([f64; 2]) -> [f64; 2],
        g: &dyn Fn([f64; 2]) -> [f64; 2],
        advect: Option<&[f64]>,
    ) -> (Csr, Vec<f64>) {
        self.assemble_viscous_transient(params, viscosity, f, g, advect, 0.0, None)
    }

    /// Variable-viscosity symmetric-stress assembly with the same BDF1
    /// mass term and previous-velocity load as `assemble_transient`.
    /// The caller supplies the viscosity at the implicit time level (or a
    /// declared lagged field), including any nonlinear constitutive update.
    ///
    /// # Panics
    /// If viscosity or the SIP penalty is not positive and finite.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_viscous_transient(
        &self,
        params: FluxParams,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        f: &dyn Fn([f64; 2]) -> [f64; 2],
        g: &dyn Fn([f64; 2]) -> [f64; 2],
        advect: Option<&[f64]>,
        inv_dt: f64,
        uprev: Option<&[f64]>,
    ) -> (Csr, Vec<f64>) {
        assert!(
            params.sigma.is_finite() && params.sigma > 0.0,
            "SIP penalty must be positive and finite"
        );
        self.assemble_operator(params, Some(viscosity), f, g, advect, inv_dt, uprev)
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn assemble_operator(
        &self,
        params: FluxParams,
        viscosity: Option<&dyn Fn(usize, [f64; 2]) -> f64>,
        f: &dyn Fn([f64; 2]) -> [f64; 2],
        g: &dyn Fn([f64; 2]) -> [f64; 2],
        advect: Option<&[f64]>,
        inv_dt: f64,
        uprev: Option<&[f64]>,
    ) -> (Csr, Vec<f64>) {
        let mesh = self.mesh;
        let nu = params.nu;
        let symmetric_stress = viscosity.is_some();
        let mut coo = Coo::new(self.n, self.n);
        let mut rhs = vec![0.0f64; self.n];
        // Constrained rows: boundary-edge normal moments are EXACTLY
        // the BDM dofs, so u·n = g·n is enforced strongly as identity
        // rows (every remaining basis function then has v·n ≡ 0 on the
        // whole boundary by dof duality — the ∮p(v·n) consistency term
        // vanishes identically). The pressure constant is back in the
        // null space, so cell 0's continuity row is REPLACED by the
        // pin p_0 = 0: its div-freeness is implied by the other cells
        // plus the prescribed net flux, not corrupted by an additive
        // penalty.
        let mut con: Vec<Option<f64>> = vec![None; self.n];
        for (e, edge) in mesh.edges.iter().enumerate() {
            if edge.tris.1 != usize::MAX {
                continue;
            }
            let (va, vb) = (mesh.verts[edge.verts.0], mesh.verts[edge.verts.1]);
            let nrm = edge.normal;
            let mut m0 = 0.0;
            let mut m1 = 0.0;
            for (gx, w) in crate::bdm::edge_gauss_pub(va, vb) {
                let gv = g(gx);
                let gn = gv[0] * nrm[0] + gv[1] * nrm[1];
                let sl = ((gx[0] - va[0]) * (vb[0] - va[0]) + (gx[1] - va[1]) * (vb[1] - va[1]))
                    / (edge.len * edge.len)
                    - 0.5;
                m0 += w * gn / edge.len;
                m1 += w * gn * sl / edge.len;
            }
            con[2 * e] = Some(m0);
            con[2 * e + 1] = Some(m1);
        }
        con[self.n_u] = Some(0.0);
        // --- Cell terms: viscous stress, div coupling, load.
        for t in 0..mesh.tris.len() {
            let basis = &self.bases[t];
            let area = mesh.areas[t];
            let p: [[f64; 2]; 3] = core::array::from_fn(|k| mesh.verts[mesh.tris[t][k]]);
            // BDM1 gradients are constant: integrate the constitutive
            // coefficient only once, then reuse it for all 36 pairs.
            let integrated_nu = viscosity.map_or(nu * area, |field| {
                tri_quad(p, area)
                    .into_iter()
                    .map(|(q, w)| w * positive_viscosity(field(t, q)))
                    .sum::<f64>()
            });
            for i in 0..6 {
                let gi = self.dof(t, i);
                // Load.
                let mut fi = 0.0;
                for (q, w) in tri_quad(p, area) {
                    let v = eval_basis(basis, i, q);
                    let fv = f(q);
                    fi += w * (fv[0] * v[0] + fv[1] * v[1]);
                }
                rhs[gi] += fi;
                if con[gi].is_none() {
                    for j in 0..6 {
                        let gj = self.dof(t, j);
                        let mut visc = 0.0;
                        for r in 0..2 {
                            for c in 0..2 {
                                visc += basis.grad[i][r][c] * basis.grad[j][r][c];
                                if symmetric_stress {
                                    visc += basis.grad[i][r][c] * basis.grad[j][c][r];
                                }
                            }
                        }
                        coo.push(gi, gj, integrated_nu * visc);
                    }
                }
                // Div coupling: b(u, q) = −∫ q div u; symmetric block.
                let pdof = self.n_u + t;
                if con[pdof].is_none() {
                    coo.push(pdof, gi, -basis.div[i] * area);
                }
                if con[gi].is_none() {
                    coo.push(gi, pdof, -basis.div[i] * area);
                }
                // BDF1 mass: (1/dt)(u, v) and (1/dt)(u_prev, v).
                if inv_dt > 0.0 && con[gi].is_none() {
                    for j in 0..6 {
                        let gj = self.dof(t, j);
                        let mut mij = 0.0;
                        for (q, w) in tri_quad(p, area) {
                            let bi = eval_basis(basis, i, q);
                            let bj = eval_basis(basis, j, q);
                            mij += w * (bi[0] * bj[0] + bi[1] * bj[1]);
                        }
                        coo.push(gi, gj, inv_dt * mij);
                        if let Some(up) = uprev {
                            rhs[gi] += inv_dt * mij * up[gj];
                        }
                    }
                }
            }
        }
        // --- Face terms.
        for edge in &mesh.edges {
            let (t0, t1) = edge.tris;
            let h = edge.len;
            let n = edge.normal;
            let (va, vb) = (mesh.verts[edge.verts.0], mesh.verts[edge.verts.1]);
            let gauss = crate::bdm::edge_gauss_pub(va, vb);
            if t1 == usize::MAX {
                // Boundary: weak Dirichlet (exterior value = g).
                for (gx, w) in gauss {
                    let gval = g(gx);
                    let local_nu = viscosity.map_or(nu, |field| positive_viscosity(field(t0, gx)));
                    let normal_height = if symmetric_stress {
                        2.0 * mesh.areas[t0] / h
                    } else {
                        h
                    };
                    let penalty = params.sigma * local_nu / normal_height;
                    for i in 0..6 {
                        let gi = self.dof(t0, i);
                        let vi = eval_basis(&self.bases[t0], i, gx);
                        let gni = viscous_normal(self.bases[t0].grad[i], n, symmetric_stress);
                        for j in 0..6 {
                            let gj = self.dof(t0, j);
                            let vj = eval_basis(&self.bases[t0], j, gx);
                            let gnj = viscous_normal(self.bases[t0].grad[j], n, symmetric_stress);
                            let val = -local_nu * (gnj[0] * vi[0] + gnj[1] * vi[1])
                                - local_nu * (gni[0] * vj[0] + gni[1] * vj[1])
                                + penalty * (vi[0] * vj[0] + vi[1] * vj[1]);
                            if con[gi].is_none() {
                                coo.push(gi, gj, w * val);
                            }
                        }
                        // rhs consistency + penalty with g.
                        let val = -local_nu * (gni[0] * gval[0] + gni[1] * gval[1])
                            + penalty * (vi[0] * gval[0] + vi[1] * gval[1]);
                        rhs[gi] += w * val;
                    }
                }
            } else {
                // Interior SIP with full jumps (normal jump ≡ 0 by
                // conformity): sides (t0, +), (t1, −).
                let cells = [t0, t1];
                let side_sign = [1.0f64, -1.0f64];
                for (gx, w) in gauss {
                    let (traction_scale, penalty) =
                        viscosity.map_or(([0.5 * nu; 2], params.sigma * nu / h), |field| {
                            let nus = cells.map(|t| positive_viscosity(field(t, gx)));
                            let delta = core::array::from_fn::<_, 2, _>(|s| {
                                nus[s] * h / (2.0 * mesh.areas[cells[s]])
                            });
                            // Opposite-side weights avoid locking the
                            // less viscous material to its stiff neighbor.
                            // Scale before adding to avoid delta0*delta1.
                            let scale = delta[0].max(delta[1]);
                            let d = delta.map(|value| value / scale);
                            let weights = [d[1] / (d[0] + d[1]), d[0] / (d[0] + d[1])];
                            let harmonic = 2.0 * delta[0].min(delta[1]) / (1.0 + d[0].min(d[1]));
                            (
                                [weights[0] * nus[0], weights[1] * nus[1]],
                                params.sigma * harmonic,
                            )
                        });
                    // Weighted one-sided viscous tractions.
                    for (sa, &ta) in cells.iter().enumerate() {
                        for i in 0..6 {
                            let gi = self.dof(ta, i);
                            let va_ = eval_basis(&self.bases[ta], i, gx);
                            let ga_ = viscous_normal(self.bases[ta].grad[i], n, symmetric_stress);
                            for (sb, &tb) in cells.iter().enumerate() {
                                for j in 0..6 {
                                    let gj = self.dof(tb, j);
                                    let vb_ = eval_basis(&self.bases[tb], j, gx);
                                    let gb_ =
                                        viscous_normal(self.bases[tb].grad[j], n, symmetric_stress);
                                    // −{τ(u)n}·[[v]] − {τ(v)n}·[[u]] + penalty [[u]]·[[v]]
                                    let jump_v = side_sign[sa];
                                    let jump_u = side_sign[sb];
                                    let val = -traction_scale[sb]
                                        * (gb_[0] * va_[0] + gb_[1] * va_[1])
                                        * jump_v
                                        - traction_scale[sa]
                                            * (ga_[0] * vb_[0] + ga_[1] * vb_[1])
                                            * jump_u
                                        + penalty
                                            * (va_[0] * vb_[0] + va_[1] * vb_[1])
                                            * jump_v
                                            * jump_u;
                                    if con[gi].is_none() {
                                        coo.push(gi, gj, w * val);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        // --- Convection (Picard): −∫(u⊗w):∇v + face upwinding.
        if let Some(wfield) = advect {
            self.add_convection(&mut coo, &mut rhs, wfield, g, &con);
        }
        // Constrained rows last: identity diagonal, prescribed value.
        for (d, c) in con.iter().enumerate() {
            if let Some(v) = c {
                coo.push(d, d, 1.0);
                rhs[d] = *v;
            }
        }
        (coo.assemble(), rhs)
    }

    /// Upwinded DG convection linearized about `wfield`.
    fn add_convection(
        &self,
        coo: &mut Coo,
        rhs: &mut [f64],
        wfield: &[f64],
        g: &dyn Fn([f64; 2]) -> [f64; 2],
        con: &[Option<f64>],
    ) {
        let mesh = self.mesh;
        let w_at = |t: usize, x: [f64; 2]| -> [f64; 2] {
            let mut v = [0.0f64; 2];
            for i in 0..6 {
                let d = self.dof(t, i);
                let b = eval_basis(&self.bases[t], i, x);
                v[0] += wfield[d] * b[0];
                v[1] += wfield[d] * b[1];
            }
            v
        };
        // Cell term: −∫ (u ⊗ w) : ∇v.
        for t in 0..mesh.tris.len() {
            let p: [[f64; 2]; 3] = core::array::from_fn(|k| mesh.verts[mesh.tris[t][k]]);
            for (q, wq) in tri_quad(p, mesh.areas[t]) {
                let wv = w_at(t, q);
                for i in 0..6 {
                    let gi = self.dof(t, i);
                    let gv = &self.bases[t].grad[i];
                    for j in 0..6 {
                        let gj = self.dof(t, j);
                        let uj = eval_basis(&self.bases[t], j, q);
                        // (u ⊗ w) : ∇v = Σ_rc u_r w_c ∂v_r/∂x_c
                        let val = uj[0] * (wv[0] * gv[0][0] + wv[1] * gv[0][1])
                            + uj[1] * (wv[0] * gv[1][0] + wv[1] * gv[1][1]);
                        if con[gi].is_none() {
                            coo.push(gi, gj, -wq * val);
                        }
                    }
                }
            }
        }
        // Face upwinding: ∫ (w·n) u_up · [[v]].
        for edge in &mesh.edges {
            let (t0, t1) = edge.tris;
            let n = edge.normal;
            let (va, vb) = (mesh.verts[edge.verts.0], mesh.verts[edge.verts.1]);
            for (gx, wq) in crate::bdm::edge_gauss_pub(va, vb) {
                let wv = w_at(t0, gx); // normal-continuous: either side works
                let wn = wv[0] * n[0] + wv[1] * n[1];
                if t1 == usize::MAX {
                    // Boundary: inflow (wn < 0) brings g, outflow takes
                    // the interior trace.
                    for i in 0..6 {
                        let gi = self.dof(t0, i);
                        let vi = eval_basis(&self.bases[t0], i, gx);
                        if con[gi].is_some() {
                            continue;
                        }
                        if wn >= 0.0 {
                            for j in 0..6 {
                                let gj = self.dof(t0, j);
                                let uj = eval_basis(&self.bases[t0], j, gx);
                                coo.push(gi, gj, wq * wn * (uj[0] * vi[0] + uj[1] * vi[1]));
                            }
                        } else {
                            let gval = g(gx);
                            rhs[gi] -= wq * wn * (gval[0] * vi[0] + gval[1] * vi[1]);
                        }
                    }
                    continue;
                }
                // Interior: upwind side by sign(wn) (wn is outward from
                // t0 by the global normal convention).
                let up = if wn >= 0.0 { t0 } else { t1 };
                let cells = [t0, t1];
                let side_sign = [1.0f64, -1.0f64];
                for (sa, &ta) in cells.iter().enumerate() {
                    for i in 0..6 {
                        let gi = self.dof(ta, i);
                        let vi = eval_basis(&self.bases[ta], i, gx);
                        for j in 0..6 {
                            let gj = self.dof(up, j);
                            let uj = eval_basis(&self.bases[up], j, gx);
                            if con[gi].is_none() {
                                coo.push(
                                    gi,
                                    gj,
                                    wq * wn * (uj[0] * vi[0] + uj[1] * vi[1]) * side_sign[sa],
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    /// Solve an assembled system: dense LU below the fixture cutoff
    /// (exact, deterministic — the saddle is unpreconditioned poison
    /// for plain Krylov), GMRES restarts above it.
    ///
    /// # Panics
    /// If the dense factorization hits a zero pivot (singular system).
    #[must_use]
    pub fn solve_linear(&self, a: &Csr, rhs: &[f64], params: FluxParams) -> FluxSolution {
        if self.n <= 1500 {
            let flat = a.to_dense();
            let f = lu(&flat, self.n).expect("saddle system nonsingular");
            let mut x = rhs.to_vec();
            f.solve(&mut x);
            let mut r2 = 0.0;
            let mut b2 = 0.0;
            let mut ax = vec![0.0f64; self.n];
            // Viscosity contrasts badly scale the velocity block against
            // the exact continuity rows. Reuse the LU for two residual
            // corrections: O(n²) work after O(n³) factorization, no extra
            // buffer, and no loosening of the divergence tolerance.
            for _ in 0..2 {
                a.spmv(&x, &mut ax);
                for i in 0..self.n {
                    ax[i] = rhs[i] - ax[i];
                }
                f.solve(&mut ax);
                for i in 0..self.n {
                    x[i] += ax[i];
                }
            }
            a.spmv(&x, &mut ax);
            for i in 0..self.n {
                r2 += (ax[i] - rhs[i]).powi(2);
                b2 += rhs[i].powi(2);
            }
            return FluxSolution {
                x,
                n_u: self.n_u,
                iters: 0,
                rel_residual: (r2 / b2.max(1e-300)).sqrt(),
            };
        }
        let op = CsrOp::general(a.clone());
        let mut st = GmresState::new(rhs, 100);
        let cycles = params.max_iters / 100;
        let _ = st.run(&op, rhs, params.tol, cycles.max(1), false);
        FluxSolution {
            x: st.x.clone(),
            n_u: self.n_u,
            iters: st.iters,
            rel_residual: st.rel_residual(),
        }
    }

    /// Solve the TRANSPOSED system (discrete adjoint), dense path.
    ///
    /// # Panics
    /// On singular systems or above the dense cutoff.
    #[must_use]
    pub fn solve_adjoint(&self, a: &Csr, j: &[f64]) -> Vec<f64> {
        assert!(self.n <= 1500, "adjoint fixture path is dense-only");
        let flat = a.to_dense();
        let mut flat_t = vec![0.0f64; self.n * self.n];
        for r in 0..self.n {
            for c in 0..self.n {
                flat_t[c * self.n + r] = flat[r * self.n + c];
            }
        }
        let f = lu(&flat_t, self.n).expect("adjoint system nonsingular");
        let mut lam = j.to_vec();
        f.solve(&mut lam);
        lam
    }

    /// Picard iteration for steady NS: Stokes start, then relinearized
    /// convection until the velocity update stalls. Returns (solution,
    /// picard iterations, last relative update).
    #[must_use]
    pub fn picard(
        &self,
        params: FluxParams,
        f: &dyn Fn([f64; 2]) -> [f64; 2],
        g: &dyn Fn([f64; 2]) -> [f64; 2],
        max_picard: usize,
        picard_tol: f64,
    ) -> (FluxSolution, usize, f64) {
        let (a, rhs) = self.assemble(params, f, g, None);
        let mut sol = self.solve_linear(&a, &rhs, params);
        let mut update = f64::INFINITY;
        let mut its = 0;
        for k in 0..max_picard {
            let (a, rhs) = self.assemble(params, f, g, Some(&sol.x[..self.n_u]));
            let next = self.solve_linear(&a, &rhs, params);
            let mut d2 = 0.0;
            let mut n2 = 0.0;
            for i in 0..self.n_u {
                d2 += (next.x[i] - sol.x[i]).powi(2);
                n2 += next.x[i].powi(2);
            }
            update = (d2 / n2.max(1e-300)).sqrt();
            sol = next;
            its = k + 1;
            if update < picard_tol {
                break;
            }
        }
        (sol, its, update)
    }

    /// Objective strain-rate invariant `sqrt(2 D:D)` [1/s] in each cell.
    /// BDM1 velocity is affine, so its gradient and rate are cellwise
    /// constant. Rigid translation and rigid rotation contribute no rate.
    ///
    /// # Panics
    /// If `velocity` is shorter than the velocity dof count.
    #[must_use]
    pub fn shear_rates(&self, velocity: &[f64]) -> Vec<f64> {
        assert!(velocity.len() >= self.n_u, "velocity dof count mismatch");
        self.bases
            .iter()
            .enumerate()
            .map(|(t, basis)| {
                let mut grad = [[0.0; 2]; 2];
                for i in 0..6 {
                    let u = velocity[self.dof(t, i)];
                    for (r, row) in grad.iter_mut().enumerate() {
                        for (c, value) in row.iter_mut().enumerate() {
                            *value += u * basis.grad[i][r][c];
                        }
                    }
                }
                (core::f64::consts::SQRT_2 * grad[0][0])
                    .hypot(core::f64::consts::SQRT_2 * grad[1][1])
                    .hypot(grad[0][1] + grad[1][0])
            })
            .collect()
    }

    /// Steady generalized-Newtonian Navier–Stokes: update both apparent
    /// viscosity and convection from the previous velocity until the
    /// relative velocity update meets `picard_tol`.
    ///
    /// `law(cell, point, shear_rate)` returns kinematic viscosity. The
    /// callback can capture a temperature/material field and evaluate a
    /// Carreau or regularized power law without a solver-specific fit.
    /// Its zero-shear value initializes the Stokes solve. Every iteration
    /// uses the full symmetric-stress operator, including face tractions.
    /// Returns the solution, iteration count, and measured last update;
    /// exhausting the budget does not claim nonlinear convergence.
    ///
    /// # Panics
    /// If constitutive values, penalty, or the convergence controls are invalid.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn picard_rheology(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> f64,
        f: &dyn Fn([f64; 2]) -> [f64; 2],
        g: &dyn Fn([f64; 2]) -> [f64; 2],
        max_picard: usize,
        picard_tol: f64,
    ) -> (FluxSolution, usize, f64) {
        assert!(
            max_picard > 0 && picard_tol.is_finite() && picard_tol > 0.0,
            "positive Picard budget and finite positive tolerance required"
        );
        let (a, rhs) = self.assemble_variable_viscosity(
            params,
            &|cell, point| law(cell, point, 0.0),
            f,
            g,
            None,
        );
        let mut sol = self.solve_linear(&a, &rhs, params);
        let mut update = f64::INFINITY;
        let mut its = 0;
        for k in 0..max_picard {
            let rates = self.shear_rates(&sol.x);
            let (a, rhs) = self.assemble_variable_viscosity(
                params,
                &|cell, point| law(cell, point, rates[cell]),
                f,
                g,
                Some(&sol.x[..self.n_u]),
            );
            let next = self.solve_linear(&a, &rhs, params);
            let mut d2 = 0.0;
            let mut n2 = 0.0;
            for i in 0..self.n_u {
                d2 += (next.x[i] - sol.x[i]).powi(2);
                n2 += next.x[i].powi(2);
            }
            update = (d2 / n2.max(1e-300)).sqrt();
            sol = next;
            its = k + 1;
            if update < picard_tol {
                break;
            }
        }
        (sol, its, update)
    }

    /// One IMEX BDF1 step: implicit Stokes, convection lagged at the
    /// previous velocity (the H(div) flux w·n is single-valued, so the
    /// upwinding is well-posed without iteration).
    #[must_use]
    pub fn bdf1_step(
        &self,
        params: FluxParams,
        f: &dyn Fn([f64; 2]) -> [f64; 2],
        g: &dyn Fn([f64; 2]) -> [f64; 2],
        uprev: &[f64],
        dt: f64,
    ) -> FluxSolution {
        let (a, rhs) = self.assemble_transient(
            params,
            f,
            g,
            Some(&uprev[..self.n_u]),
            1.0 / dt,
            Some(uprev),
        );
        self.solve_linear(&a, &rhs, params)
    }

    /// Cell pressures shifted to zero area-weighted mean.
    #[must_use]
    pub fn pressure(&self, x: &[f64]) -> Vec<f64> {
        let mut total = 0.0;
        let mut mean = 0.0;
        for t in 0..self.mesh.tris.len() {
            total += self.mesh.areas[t];
            mean += self.mesh.areas[t] * x[self.n_u + t];
        }
        mean /= total;
        (0..self.mesh.tris.len())
            .map(|t| x[self.n_u + t] - mean)
            .collect()
    }

    /// Velocity at a point (containing-cell search by nearest
    /// centroid — fixture-scale linear scan).
    #[must_use]
    pub fn velocity_at(&self, x: &[f64], pt: [f64; 2]) -> [f64; 2] {
        let mut best = 0;
        let mut bd = f64::INFINITY;
        for (t, c) in self.mesh.centroids.iter().enumerate() {
            let d = (c[0] - pt[0]).powi(2) + (c[1] - pt[1]).powi(2);
            if d < bd {
                bd = d;
                best = t;
            }
        }
        let mut v = [0.0f64; 2];
        for i in 0..6 {
            let b = eval_basis(&self.bases[best], i, pt);
            v[0] += x[self.dof(best, i)] * b[0];
            v[1] += x[self.dof(best, i)] * b[1];
        }
        v
    }

    /// Kinetic energy of a discrete velocity.
    #[must_use]
    pub fn kinetic_energy(&self, x: &[f64]) -> f64 {
        let mut ke = 0.0;
        for t in 0..self.mesh.tris.len() {
            let p: [[f64; 2]; 3] = core::array::from_fn(|k| self.mesh.verts[self.mesh.tris[t][k]]);
            for (q, w) in tri_quad(p, self.mesh.areas[t]) {
                let mut uh = [0.0f64; 2];
                for i in 0..6 {
                    let b = eval_basis(&self.bases[t], i, q);
                    uh[0] += x[self.dof(t, i)] * b[0];
                    uh[1] += x[self.dof(t, i)] * b[1];
                }
                ke += 0.5 * w * (uh[0] * uh[0] + uh[1] * uh[1]);
            }
        }
        ke
    }

    /// L2 velocity error and worst per-cell divergence.
    #[must_use]
    pub fn velocity_error(&self, x: &[f64], exact: &dyn Fn([f64; 2]) -> [f64; 2]) -> (f64, f64) {
        let mesh = self.mesh;
        let mut l2 = 0.0;
        let mut worst_div = 0.0f64;
        for t in 0..mesh.tris.len() {
            let p: [[f64; 2]; 3] = core::array::from_fn(|k| mesh.verts[mesh.tris[t][k]]);
            let mut div = 0.0;
            for i in 0..6 {
                div += x[self.dof(t, i)] * self.bases[t].div[i];
            }
            worst_div = worst_div.max(div.abs());
            for (q, w) in tri_quad(p, mesh.areas[t]) {
                let mut uh = [0.0f64; 2];
                for i in 0..6 {
                    let b = eval_basis(&self.bases[t], i, q);
                    uh[0] += x[self.dof(t, i)] * b[0];
                    uh[1] += x[self.dof(t, i)] * b[1];
                }
                let ue = exact(q);
                l2 += w * ((ue[0] - uh[0]).powi(2) + (ue[1] - uh[1]).powi(2));
            }
        }
        (l2.sqrt(), worst_div)
    }
}
