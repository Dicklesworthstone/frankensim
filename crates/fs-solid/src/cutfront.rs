//! The CutFEM frontend: vector Q1 linear elasticity directly on an SDF
//! over an fs-cutfem background quadtree — the topology-optimization
//! path (topo-simp solves THIS physics on these grids with zero
//! meshing). Reuses fs-cutfem's certified classification and cut
//! quadrature verbatim; adds the VECTOR Nitsche boundary operator
//! (full traction consistency terms, penalty scaled by (λ+2μ)/h tied
//! to the certified cut cells) and the componentwise ghost penalty.
//!
//! v1 surface: UNIFORM background grids (every leaf at one level); the
//! graded/hanging machinery is proven in fs-cutfem's scalar battery
//! and its vector lift is a recorded no-claim.

use crate::SolidError;
use crate::linear::{Jacobi, PlaneKind, lame};
use fs_cutfem::quad::{CutRules, cut_cell_rules, tensor_gauss};
use fs_cutfem::{CutSdf, Quadtree};
use fs_solver::krylov::CgState;
use fs_solver::op::CsrOp;
use fs_sparse::Coo;
use std::collections::{BTreeMap, BTreeSet};

/// A cell key pair (ghost-face identity).
type FaceKey = ((u32, u32, u32), (u32, u32, u32));

/// A CutFEM linear-elasticity problem on `Ω = {φ < 0}`.
pub struct CutElasticity<'a> {
    /// The (uniform) background quadtree.
    pub grid: &'a Quadtree,
    /// The level set.
    pub sdf: &'a dyn CutSdf,
    /// Young's modulus.
    pub youngs: f64,
    /// Poisson ratio.
    pub poisson: f64,
    /// Nitsche constant β (applied as β(λ+2μ)/h).
    pub nitsche_beta: f64,
    /// Ghost constant γ (applied as γ(λ+2μ)h per face).
    pub ghost_gamma: f64,
    /// Cut-quadrature depth.
    pub quad_depth: u32,
    /// Strong zero-Dirichlet clamp on design-box boundary nodes
    /// satisfying the predicate (topology-optimization frontends clamp
    /// on the box edge, where Nitsche-on-Γ cannot reach).
    pub clamp: Option<&'a dyn Fn(f64, f64) -> bool>,
    /// Dead traction on design-box boundary edges of active cells.
    pub boundary_traction: Option<&'a dyn Fn(f64, f64) -> [f64; 2]>,
    /// Treat Γ as TRACTION-FREE (natural) instead of Nitsche
    /// Dirichlet — the topology-optimization void boundary. Ghost
    /// penalty still stabilizes the cut cells; support must come from
    /// clamps.
    pub traction_free_interface: bool,
}

/// The CutFEM solution: nodal displacements plus error measurement.
pub struct CutSolution {
    nodal: BTreeMap<(u32, u32), [f64; 2]>,
    active: Vec<(u32, u32, u32)>,
    rules: BTreeMap<(u32, u32, u32), CutRules>,
    /// CG iterations.
    pub iters: usize,
}

impl CutSolution {
    /// Nodal displacements (topology-optimization consumers sample
    /// energy densities and load work from these).
    #[must_use]
    pub fn nodal(&self) -> &BTreeMap<(u32, u32), [f64; 2]> {
        &self.nodal
    }

    /// The active cells of the solve.
    #[must_use]
    pub fn active_cells(&self) -> &[(u32, u32, u32)] {
        &self.active
    }
}

impl CutElasticity<'_> {
    /// Assemble and solve `−div σ(u) = f` in Ω, `u = g` on Γ (Nitsche).
    ///
    /// # Errors
    /// [`SolidError::SolveFailed`].
    #[allow(clippy::too_many_lines)] // classify + bulk + Nitsche + ghost, one deterministic sweep
    pub fn solve(
        &self,
        f: &dyn Fn(f64, f64) -> [f64; 2],
        g: &dyn Fn(f64, f64) -> [f64; 2],
    ) -> Result<CutSolution, SolidError> {
        let (lambda, mu) = lame(self.youngs, self.poisson, PlaneKind::Strain);
        let mut active: Vec<(u32, u32, u32)> = Vec::new();
        let mut cut: BTreeSet<(u32, u32, u32)> = BTreeSet::new();
        let mut rules: BTreeMap<(u32, u32, u32), CutRules> = BTreeMap::new();
        for c in self.grid.leaves() {
            let (lo, hi) = self.grid.rect(c);
            let iv = self.sdf.enclose(lo, hi);
            if iv.hi() < 0.0 {
                active.push(c);
            } else if iv.lo() <= 0.0 {
                let r = cut_cell_rules(self.sdf, lo, hi, self.quad_depth);
                let area = (hi[0] - lo[0]) * (hi[1] - lo[1]);
                let w: f64 = r.bulk.iter().map(|&(_, w)| w).sum();
                if w >= 1e-12 * area {
                    active.push(c);
                    cut.insert(c);
                    rules.insert(c, r);
                }
            }
        }
        let active_set: BTreeSet<_> = active.iter().copied().collect();
        let mut node_ids: BTreeMap<(u32, u32), usize> = BTreeMap::new();
        for &c in &active {
            for n in self.grid.corner_nodes(c) {
                let next = node_ids.len();
                node_ids.entry(n).or_insert(next);
            }
        }
        let ndof = 2 * node_ids.len();
        let ext = self.grid.node_extent();
        let mut clamped = vec![false; ndof];
        if let Some(pred) = self.clamp {
            for (&n, &id) in &node_ids {
                if (n.0 == 0 || n.0 == ext || n.1 == 0 || n.1 == ext) && {
                    let p = self.grid.node_pos(n);
                    pred(p[0], p[1])
                } {
                    clamped[2 * id] = true;
                    clamped[2 * id + 1] = true;
                }
            }
        }
        let mut coo = Coo::new(ndof, ndof);
        let mut rhs = vec![0.0f64; ndof];
        for &c in &active {
            let (lo, hi) = self.grid.rect(c);
            let corners = self.grid.corner_nodes(c);
            let ids: Vec<usize> = corners.iter().map(|n| node_ids[n]).collect();
            let h = self.grid.cell_h(c);
            let mut k = [[0.0f64; 8]; 8];
            let mut fl = [0.0f64; 8];
            let own_rule;
            let bulk: &[([f64; 2], f64)] = if cut.contains(&c) {
                &rules[&c].bulk
            } else {
                own_rule = {
                    let mut v = Vec::with_capacity(9);
                    tensor_gauss(lo, hi, &mut v);
                    v
                };
                &own_rule
            };
            for &(p, w) in bulk {
                let (nv, gr) = q1(lo, hi, p);
                let fv = f(p[0], p[1]);
                for a in 0..4 {
                    for ca in 0..2 {
                        let ba = b_row(gr[a], ca);
                        let db = d_mul(lambda, mu, ba);
                        for b in 0..4 {
                            for cb in 0..2 {
                                let bb = b_row(gr[b], cb);
                                k[2 * a + ca][2 * b + cb] +=
                                    w * (db[0] * bb[0] + db[1] * bb[1] + db[2] * bb[2]);
                            }
                        }
                        fl[2 * a + ca] += w * nv[a] * fv[ca];
                    }
                }
            }
            if cut.contains(&c) && !self.traction_free_interface {
                let pen = self.nitsche_beta * (lambda + 2.0 * mu) / h;
                for &(p, w, nrm) in &rules[&c].iface {
                    let (nv, gr) = q1(lo, hi, p);
                    let gval = g(p[0], p[1]);
                    // Traction of shape (a, i): σ(N_a e_i) · n.
                    let trac = |a: usize, i: usize| -> [f64; 2] {
                        let gvec = gr[a];
                        let gn = gvec[0] * nrm[0] + gvec[1] * nrm[1];
                        let mut t = [0.0f64; 2];
                        for kk in 0..2 {
                            t[kk] = lambda * gvec[i] * nrm[kk] + mu * nrm[i] * gvec[kk];
                        }
                        t[i] += mu * gn;
                        t
                    };
                    for a in 0..4 {
                        for ca in 0..2 {
                            let ta = trac(a, ca);
                            for b in 0..4 {
                                for cb in 0..2 {
                                    let tb = trac(b, cb);
                                    let mut v = pen * nv[a] * nv[b] * f64::from(ca == cb);
                                    v -= ta[cb] * nv[b];
                                    v -= nv[a] * tb[ca];
                                    k[2 * a + ca][2 * b + cb] += w * v;
                                }
                            }
                            let mut r = pen * nv[a] * gval[ca];
                            r -= ta[0] * gval[0] + ta[1] * gval[1];
                            fl[2 * a + ca] += w * r;
                        }
                    }
                }
            }
            for a in 0..8 {
                let ia = 2 * ids[a / 2] + a % 2;
                if clamped[ia] {
                    continue;
                }
                rhs[ia] += fl[a];
                for b in 0..8 {
                    let ib = 2 * ids[b / 2] + b % 2;
                    if k[a][b] != 0.0 && !clamped[ib] {
                        coo.push(ia, ib, k[a][b]);
                    }
                }
            }
        }
        // Design-box boundary tractions: integrate over box-edge faces
        // of ACTIVE cells (2-pt Gauss per edge), loading material only.
        if let Some(tr) = self.boundary_traction {
            for &c in &active {
                let (lv, i, j) = c;
                let nmax = 1u32 << lv;
                let corners = self.grid.corner_nodes(c);
                let edges: [(bool, [usize; 2]); 4] = [
                    (j == 0, [0, 1]),
                    (i + 1 == nmax, [1, 2]),
                    (j + 1 == nmax, [2, 3]),
                    (i == 0, [3, 0]),
                ];
                for (on_boundary, corner_ids) in edges {
                    if !on_boundary {
                        continue;
                    }
                    let pa = self.grid.node_pos(corners[corner_ids[0]]);
                    let pb = self.grid.node_pos(corners[corner_ids[1]]);
                    let len = (pb[0] - pa[0]).hypot(pb[1] - pa[1]);
                    let gpt = 0.5 / 3.0f64.sqrt();
                    for t in [0.5 - gpt, 0.5 + gpt] {
                        let p = [pa[0] + t * (pb[0] - pa[0]), pa[1] + t * (pb[1] - pa[1])];
                        if self.sdf.value(p) > 0.0 {
                            continue;
                        }
                        let tv = tr(p[0], p[1]);
                        let w = 0.5 * len;
                        for (ci, shape) in [(corner_ids[0], 1.0 - t), (corner_ids[1], t)] {
                            let id = node_ids[&corners[ci]];
                            for (comp, tc) in tv.iter().enumerate() {
                                let dof = 2 * id + comp;
                                if !clamped[dof] {
                                    rhs[dof] += w * shape * tc;
                                }
                            }
                        }
                    }
                }
            }
        }
        for (dof, is_clamped) in clamped.iter().enumerate() {
            if *is_clamped {
                coo.push(dof, dof, 1.0);
            }
        }
        // Componentwise ghost penalty on equal-level faces of cut cells.
        if self.ghost_gamma > 0.0 {
            let mut seen: BTreeSet<FaceKey> = BTreeSet::new();
            for &c in &cut {
                let (lv, i, j) = c;
                let nmax = 1u32 << lv;
                let neighbors = [
                    (i > 0).then(|| (lv, i - 1, j)),
                    (i + 1 < nmax).then_some((lv, i + 1, j)),
                    (j > 0).then(|| (lv, i, j - 1)),
                    (j + 1 < nmax).then_some((lv, i, j + 1)),
                ];
                for nb in neighbors.into_iter().flatten() {
                    if !active_set.contains(&nb) {
                        continue;
                    }
                    let key = if c < nb { (c, nb) } else { (nb, c) };
                    if !seen.insert(key) {
                        continue;
                    }
                    self.ghost_face(
                        key.0,
                        key.1,
                        lambda + 2.0 * mu,
                        &node_ids,
                        &clamped,
                        &mut coo,
                    );
                }
            }
        }
        let a = coo.assemble();
        let m = Jacobi::new(&a);
        let op = CsrOp::symmetric(a);
        let mut st = CgState::new(&op, &m, &rhs);
        let _ = st.run(&op, &m, 1e-12, 60_000);
        let rr = st.rel_residual();
        if !rr.is_finite() || rr > 1e-8 {
            return Err(SolidError::SolveFailed {
                iters: st.iters,
                rel_residual: rr,
            });
        }
        let nodal = node_ids
            .iter()
            .map(|(&n, &id)| (n, [st.x[2 * id], st.x[2 * id + 1]]))
            .collect();
        Ok(CutSolution {
            nodal,
            active,
            rules,
            iters: st.iters,
        })
    }

    #[allow(clippy::too_many_arguments)] // one face, one scatter
    fn ghost_face(
        &self,
        ca: (u32, u32, u32),
        cb: (u32, u32, u32),
        stiff: f64,
        node_ids: &BTreeMap<(u32, u32), usize>,
        clamped: &[bool],
        coo: &mut Coo,
    ) {
        let (lo_a, hi_a) = self.grid.rect(ca);
        let (lo_b, hi_b) = self.grid.rect(cb);
        let h = self.grid.cell_h(ca);
        let axis = usize::from(ca.1 == cb.1); // same i → horizontal face
        let (t0, t1) = if axis == 0 {
            (lo_a[1], hi_a[1])
        } else {
            (lo_a[0], hi_a[0])
        };
        let nrm = if axis == 0 { [1.0, 0.0] } else { [0.0, 1.0] };
        let xf = if axis == 0 { hi_a[0] } else { hi_a[1] };
        let corners_a = self.grid.corner_nodes(ca);
        let corners_b = self.grid.corner_nodes(cb);
        let gpt = 0.5 / 3.0f64.sqrt();
        let wq = 0.5 * (t1 - t0);
        let mut jump: BTreeMap<(u32, u32), [f64; 2]> = BTreeMap::new();
        for (qi, t) in [0.5 - gpt, 0.5 + gpt].into_iter().enumerate() {
            let tv = t0 + t * (t1 - t0);
            let p = if axis == 0 { [xf, tv] } else { [tv, xf] };
            let (_, gra) = q1(lo_a, hi_a, p);
            let (_, grb) = q1(lo_b, hi_b, p);
            for a in 0..4 {
                jump.entry(corners_a[a]).or_default()[qi] +=
                    gra[a][0] * nrm[0] + gra[a][1] * nrm[1];
                jump.entry(corners_b[a]).or_default()[qi] -=
                    grb[a][0] * nrm[0] + grb[a][1] * nrm[1];
            }
        }
        let scale = self.ghost_gamma * stiff * h * wq;
        let entries: Vec<((u32, u32), [f64; 2])> = jump.into_iter().collect();
        for (na, ja) in &entries {
            for (nb, jb) in &entries {
                let v = scale * (ja[0] * jb[0] + ja[1] * jb[1]);
                if v == 0.0 {
                    continue;
                }
                for c in 0..2 {
                    let (ia, ib) = (2 * node_ids[na] + c, 2 * node_ids[nb] + c);
                    if !clamped[ia] && !clamped[ib] {
                        coo.push(ia, ib, v);
                    }
                }
            }
        }
    }

    /// L2/H1 displacement errors with one-deeper cut quadrature.
    #[must_use]
    pub fn l2_h1_error(
        &self,
        sol: &CutSolution,
        exact: &dyn Fn(f64, f64) -> [f64; 2],
        grad_exact: &dyn Fn(f64, f64) -> [[f64; 2]; 2],
    ) -> (f64, f64) {
        let mut l2 = 0.0f64;
        let mut h1 = 0.0f64;
        for &c in &sol.active {
            let (lo, hi) = self.grid.rect(c);
            let corners = self.grid.corner_nodes(c);
            let vals: Vec<[f64; 2]> = corners.iter().map(|n| sol.nodal[n]).collect();
            let refined;
            let rule: &[([f64; 2], f64)] = if sol.rules.contains_key(&c) {
                refined = cut_cell_rules(self.sdf, lo, hi, self.quad_depth + 1).bulk;
                &refined
            } else {
                refined = {
                    let mut v = Vec::with_capacity(9);
                    tensor_gauss(lo, hi, &mut v);
                    v
                };
                &refined
            };
            for &(p, w) in rule {
                let (nv, gr) = q1(lo, hi, p);
                let mut uh = [0.0f64; 2];
                let mut guh = [[0.0f64; 2]; 2];
                for a in 0..4 {
                    for cc in 0..2 {
                        uh[cc] += nv[a] * vals[a][cc];
                        guh[cc][0] += gr[a][0] * vals[a][cc];
                        guh[cc][1] += gr[a][1] * vals[a][cc];
                    }
                }
                let ue = exact(p[0], p[1]);
                let ge = grad_exact(p[0], p[1]);
                for cc in 0..2 {
                    let e = ue[cc] - uh[cc];
                    l2 += w * e * e;
                    for r in 0..2 {
                        let d = ge[cc][r] - guh[cc][r];
                        h1 += w * d * d;
                    }
                }
            }
        }
        (l2.max(0.0).sqrt(), h1.max(0.0).sqrt())
    }
}

fn b_row(g: [f64; 2], c: usize) -> [f64; 3] {
    if c == 0 {
        [g[0], 0.0, g[1]]
    } else {
        [0.0, g[1], g[0]]
    }
}

fn d_mul(lambda: f64, mu: f64, b: [f64; 3]) -> [f64; 3] {
    [
        (lambda + 2.0 * mu) * b[0] + lambda * b[1],
        lambda * b[0] + (lambda + 2.0 * mu) * b[1],
        mu * b[2],
    ]
}

/// Q1 shapes on an axis-aligned cell (fs-cutfem corner order).
fn q1(lo: [f64; 2], hi: [f64; 2], p: [f64; 2]) -> ([f64; 4], [[f64; 2]; 4]) {
    let hx = hi[0] - lo[0];
    let hy = hi[1] - lo[1];
    let xi = (p[0] - lo[0]) / hx;
    let et = (p[1] - lo[1]) / hy;
    (
        [
            (1.0 - xi) * (1.0 - et),
            xi * (1.0 - et),
            xi * et,
            (1.0 - xi) * et,
        ],
        [
            [-(1.0 - et) / hx, -(1.0 - xi) / hy],
            [(1.0 - et) / hx, -xi / hy],
            [et / hx, xi / hy],
            [-et / hx, (1.0 - xi) / hy],
        ],
    )
}
