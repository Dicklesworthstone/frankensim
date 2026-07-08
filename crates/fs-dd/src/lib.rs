//! fs-dd — DOMAIN DECOMPOSITION (plan §8.9, bead tfz.11; [M] — behind
//! the `bddc` / `sheaf-coarse` features until Gauntlet evidence
//! promotes them). BDDC substructuring on the 2-D model problem:
//! corners-primal coarse spaces with weighted interface averaging, the
//! SHEAF-HARMONIC edge enrichment (the same cellular-sheaf machinery as
//! watertightness — Bet 11 earning its keep twice), coefficient-jump
//! robustness, and CCD-aligned partitioning metrics. Layer L3.
//!
//! The model problem is the variable-coefficient 5-point Laplacian on
//! the unit square (spectrally equivalent to P1 FEM — the setting BDDC
//! theory addresses): `n = s·m` cells per side, `s²` square subdomains
//! of `m²` cells. The solver is CG on the INTERFACE Schur system with
//! the BDDC preconditioner; condition numbers are estimated from the
//! CG Lanczos coefficients.
#![cfg(feature = "bddc")]

use fs_math::det;
use std::collections::BTreeMap;

/// The structured 2-D decomposition: `s × s` subdomains of `m × m`
/// cells, `n = s·m` cells per side, Dirichlet outer boundary.
#[derive(Debug, Clone)]
pub struct Decomposition {
    /// Subdomains per side.
    pub s: usize,
    /// Cells per subdomain side (H/h).
    pub m: usize,
    /// Per-cell coefficient field, row-major `n × n`.
    pub rho: Vec<f64>,
}

/// Node index helpers on the `(n+1)²` lattice (Dirichlet rim removed
/// later; we index the FULL lattice and mask).
impl Decomposition {
    /// Uniform-coefficient decomposition.
    #[must_use]
    pub fn uniform(s: usize, m: usize) -> Decomposition {
        Decomposition {
            s,
            m,
            rho: vec![1.0; s * m * s * m],
        }
    }

    /// Checkerboard coefficients: `rho_hi` on subdomains with odd
    /// (i + j), 1 elsewhere — the classic jump-robustness fixture.
    #[must_use]
    pub fn checkerboard(s: usize, m: usize, rho_hi: f64) -> Decomposition {
        let n = s * m;
        let mut rho = vec![1.0; n * n];
        for cy in 0..n {
            for cx in 0..n {
                let (si, sj) = (cx / m, cy / m);
                if (si + sj) % 2 == 1 {
                    rho[cy * n + cx] = rho_hi;
                }
            }
        }
        Decomposition { s, m, rho }
    }

    fn n(&self) -> usize {
        self.s * self.m
    }

    fn np(&self) -> usize {
        self.n() + 1
    }

    fn node(&self, x: usize, y: usize) -> usize {
        y * self.np() + x
    }

    /// True for nodes on the outer Dirichlet rim.
    fn on_rim(&self, x: usize, y: usize) -> bool {
        x == 0 || y == 0 || x == self.n() || y == self.n()
    }

    /// Harmonic mean of the four adjacent cell coefficients at an edge
    /// between nodes (the standard variable-coefficient stencil).
    fn edge_coeff(&self, ax: usize, ay: usize, bx: usize, by: usize) -> f64 {
        let n = self.n();
        let mut acc = 0.0f64;
        let mut cnt = 0.0f64;
        // The two cells flanking the edge (clamped at the boundary).
        let cells: [(Option<usize>, Option<usize>); 2] = if ay == by {
            let cx = Some(ax.min(bx));
            [(cx, ay.checked_sub(1)), (cx, Some(ay))]
        } else {
            let cy = Some(ay.min(by));
            [(ax.checked_sub(1), cy), (Some(ax), cy)]
        };
        for (cx, cy) in cells {
            if let (Some(cx), Some(cy)) = (cx, cy)
                && cx < n
                && cy < n
            {
                acc += self.rho[cy * n + cx];
                cnt += 1.0;
            }
        }
        if cnt > 0.0 { acc / cnt } else { 1.0 }
    }

    /// Apply the global variable-coefficient 5-point operator to a
    /// full-lattice vector (rim entries treated as zero/Dirichlet) —
    /// the whole-system oracle tests verify the Schur path against.
    #[must_use]
    pub fn apply_global(&self, x: &[f64]) -> Vec<f64> {
        let np = self.np();
        let mut out = vec![0.0f64; np * np];
        for y in 1..self.n() {
            for xx in 1..self.n() {
                let c = self.node(xx, y);
                let mut acc = 0.0f64;
                for (nx, ny) in [(xx - 1, y), (xx + 1, y), (xx, y - 1), (xx, y + 1)] {
                    let w = self.edge_coeff(xx, y, nx, ny);
                    let nb = if self.on_rim(nx, ny) {
                        0.0
                    } else {
                        x[self.node(nx, ny)]
                    };
                    acc += w * (x[c] - nb);
                }
                out[c] = acc;
            }
        }
        out
    }

    /// Node classification: which subdomain(s) a node touches.
    fn touching_subdomains(&self, x: usize, y: usize) -> Vec<(usize, usize)> {
        let mut out = Vec::with_capacity(4);
        let m = self.m;
        let s = self.s;
        for sj in 0..s {
            for si in 0..s {
                let (x0, x1) = (si * m, (si + 1) * m);
                let (y0, y1) = (sj * m, (sj + 1) * m);
                if x >= x0 && x <= x1 && y >= y0 && y <= y1 {
                    out.push((si, sj));
                }
            }
        }
        out
    }

    /// Interface nodes: non-rim nodes shared by ≥2 subdomains.
    /// Corners: shared by 4 (the subdomain cross points).
    fn classify(&self) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
        let mut interior = Vec::new();
        let mut interface = Vec::new();
        let mut corners = Vec::new();
        for y in 0..=self.n() {
            for x in 0..=self.n() {
                if self.on_rim(x, y) {
                    continue;
                }
                let t = self.touching_subdomains(x, y).len();
                let idx = self.node(x, y);
                match t {
                    1 => interior.push(idx),
                    2 => interface.push(idx),
                    _ => corners.push(idx),
                }
            }
        }
        (interior, interface, corners)
    }
}

/// Dense SPD Cholesky solve (small local systems).
fn cholesky_factor(a: &mut [Vec<f64>]) {
    let n = a.len();
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i][j];
            {
                let (ri, rj) = (&a[i], &a[j]);
                for (x, y) in ri[..j].iter().zip(&rj[..j]) {
                    sum -= x * y;
                }
            }
            if i == j {
                assert!(sum > 0.0, "local matrix must be SPD");
                a[i][i] = det::sqrt(sum);
            } else {
                a[i][j] = sum / a[j][j];
            }
        }
    }
}

fn cholesky_solve(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut y = b.to_vec();
    for i in 0..n {
        for k in 0..i {
            y[i] -= l[i][k] * y[k];
        }
        y[i] /= l[i][i];
    }
    for i in (0..n).rev() {
        for k in (i + 1)..n {
            y[i] -= l[k][i] * y[k];
        }
        y[i] /= l[i][i];
    }
    y
}

/// One subdomain's local data: node lists and factored local matrix
/// (interior block for Schur solves; full local for Neumann solves).
struct Subdomain {
    /// Global node ids owned (interior of this subdomain).
    interior: Vec<usize>,
    /// Global interface/corner node ids on this subdomain's boundary
    /// (excluding the outer Dirichlet rim).
    boundary: Vec<usize>,
    /// Factored K_II (interior-interior).
    l_ii: Vec<Vec<f64>>,
    /// K_IB (interior × boundary).
    k_ib: Vec<Vec<f64>>,
    /// K_BB (boundary × boundary).
    k_bb: Vec<Vec<f64>>,
}

/// The BDDC solver context for one decomposition.
pub struct Bddc {
    decomp: Decomposition,
    subs: Vec<Subdomain>,
    /// All interface+corner ("gamma") nodes, sorted.
    gamma: Vec<usize>,
    /// gamma node -> position.
    gpos: BTreeMap<usize, usize>,
    /// Corner nodes (subset of gamma), sorted.
    corners: Vec<usize>,
    /// Multiplicity weights per gamma node (counting functions).
    weight: Vec<f64>,
    /// Extra coarse edge-average vectors (sheaf enrichment), each a
    /// gamma-sized vector.
    edge_modes: Vec<Vec<f64>>,
    /// Per subdomain: free (non-corner) boundary positions + the
    /// FACTORED local Schur block on them (computed once).
    local_schur: Vec<(Vec<usize>, Vec<Vec<f64>>)>,
    /// Per subdomain: orthonormal deflation vectors (its edges' average
    /// modes on the free boundary) — local corrections are projected
    /// AWAY from the coarse edge space so the two do not fight.
    local_deflate: Vec<Vec<Vec<f64>>>,
    /// The coarse basis (corners + edge modes) and the FACTORED coarse
    /// Galerkin matrix (computed once).
    coarse_basis: Vec<Vec<f64>>,
    coarse_l: Vec<Vec<f64>>,
}

fn local_stiffness(d: &Decomposition, nodes: &[usize]) -> Vec<Vec<f64>> {
    let np = d.np();
    let pos: BTreeMap<usize, usize> = nodes.iter().copied().zip(0..).collect();
    let mut a = vec![vec![0.0f64; nodes.len()]; nodes.len()];
    for (&g, &i) in &pos {
        let (x, y) = (g % np, g / np);
        for (nx, ny) in [
            (x.wrapping_sub(1), y),
            (x + 1, y),
            (x, y.wrapping_sub(1)),
            (x, y + 1),
        ] {
            if nx > d.n() || ny > d.n() {
                continue;
            }
            let w = d.edge_coeff(x, y, nx, ny);
            a[i][i] += w;
            if !d.on_rim(nx, ny)
                && let Some(&j) = pos.get(&d.node(nx, ny))
            {
                a[i][j] -= w;
            }
        }
    }
    a
}

impl Bddc {
    /// Build the substructured context (factorizations per subdomain).
    #[must_use]
    #[allow(clippy::too_many_lines)]
    pub fn new(decomp: Decomposition, with_edge_modes: bool) -> Bddc {
        let (_, interface, corners) = decomp.classify();
        let mut gamma: Vec<usize> = interface.iter().chain(&corners).copied().collect();
        gamma.sort_unstable();
        let gpos: BTreeMap<usize, usize> = gamma.iter().copied().zip(0..).collect();
        let np = decomp.np();
        let mut subs = Vec::new();
        for sj in 0..decomp.s {
            for si in 0..decomp.s {
                let (x0, x1) = (si * decomp.m, (si + 1) * decomp.m);
                let (y0, y1) = (sj * decomp.m, (sj + 1) * decomp.m);
                let mut interior = Vec::new();
                let mut boundary = Vec::new();
                for y in y0..=y1 {
                    for x in x0..=x1 {
                        if decomp.on_rim(x, y) {
                            continue;
                        }
                        let g = decomp.node(x, y);
                        let strictly_inside = x > x0 && x < x1 && y > y0 && y < y1;
                        if strictly_inside {
                            interior.push(g);
                        } else {
                            boundary.push(g);
                        }
                    }
                }
                // Local blocks.
                let all: Vec<usize> = interior.iter().chain(&boundary).copied().collect();
                let a = local_stiffness(&decomp, &all);
                let ni = interior.len();
                let nb = boundary.len();
                let mut k_ii = vec![vec![0.0; ni]; ni];
                let mut k_ib = vec![vec![0.0; nb]; ni];
                let mut k_bb = vec![vec![0.0; nb]; nb];
                for i in 0..ni {
                    for j in 0..ni {
                        k_ii[i][j] = a[i][j];
                    }
                    for j in 0..nb {
                        k_ib[i][j] = a[i][ni + j];
                    }
                }
                for i in 0..nb {
                    for j in 0..nb {
                        k_bb[i][j] = a[ni + i][ni + j];
                    }
                }
                if ni > 0 {
                    cholesky_factor(&mut k_ii);
                }
                subs.push(Subdomain {
                    interior,
                    boundary,
                    l_ii: k_ii,
                    k_ib,
                    k_bb,
                });
            }
        }
        // Counting weights: 1 / (number of subdomains sharing the node).
        let mut weight = vec![0.0f64; gamma.len()];
        for sub in &subs {
            for &g in &sub.boundary {
                weight[gpos[&g]] += 1.0;
            }
        }
        for w in &mut weight {
            *w = 1.0 / *w;
        }
        // Sheaf-derived edge modes: one average constraint per open
        // interface edge (nodes shared by exactly 2 subdomains, grouped
        // by the subdomain pair).
        let mut edge_modes = Vec::new();
        if with_edge_modes {
            let mut groups: BTreeMap<(usize, usize, usize, usize), Vec<usize>> = BTreeMap::new();
            for &g in &interface {
                let (x, y) = (g % np, g / np);
                let t = decomp.touching_subdomains(x, y);
                if t.len() == 2 {
                    let key = (t[0].0, t[0].1, t[1].0, t[1].1);
                    groups.entry(key).or_default().push(g);
                }
            }
            for nodes in groups.values() {
                let mut v = vec![0.0f64; gamma.len()];
                #[allow(clippy::cast_precision_loss)]
                let inv = 1.0 / nodes.len() as f64;
                for &g in nodes {
                    v[gpos[&g]] = inv;
                }
                edge_modes.push(v);
            }
        }
        let mut ctx = Bddc {
            decomp,
            subs,
            gamma,
            gpos,
            corners,
            weight,
            edge_modes,
            local_schur: Vec::new(),
            local_deflate: Vec::new(),
            coarse_basis: Vec::new(),
            coarse_l: Vec::new(),
        };
        ctx.factor_local_schur();
        ctx.build_local_deflation();
        ctx.factor_coarse();
        ctx
    }

    /// Factor each subdomain's free-boundary Schur block once.
    fn factor_local_schur(&mut self) {
        let corner_set: std::collections::BTreeSet<usize> = self.corners.iter().copied().collect();
        let mut out = Vec::with_capacity(self.subs.len());
        for sub in &self.subs {
            let nb = sub.boundary.len();
            let ni = sub.interior.len();
            let free: Vec<usize> = (0..nb)
                .filter(|&i| !corner_set.contains(&sub.boundary[i]))
                .collect();
            let nf = free.len();
            let mut s_loc = vec![vec![0.0f64; nf]; nf];
            for (a, &fa) in free.iter().enumerate() {
                let mut col: Vec<f64> = sub.k_bb.iter().map(|row| row[fa]).collect();
                if ni > 0 {
                    let rhs: Vec<f64> = sub.k_ib.iter().map(|row| row[fa]).collect();
                    let w = cholesky_solve(&sub.l_ii, &rhs);
                    for (i, c) in col.iter_mut().enumerate() {
                        let acc: f64 = sub.k_ib.iter().zip(&w).map(|(row, wk)| row[i] * wk).sum();
                        *c -= acc;
                    }
                }
                for (b, &fb) in free.iter().enumerate() {
                    s_loc[b][a] = col[fb];
                }
            }
            for (i, row) in s_loc.iter_mut().enumerate() {
                row[i] += 1e-12;
            }
            if nf > 0 {
                cholesky_factor(&mut s_loc);
            }
            out.push((free, s_loc));
        }
        self.local_schur = out;
    }

    /// Per-subdomain deflation: each subdomain's slice of every edge
    /// mode, restricted to its free boundary and orthonormalized.
    fn build_local_deflation(&mut self) {
        let mut out = Vec::with_capacity(self.subs.len());
        for (sub, (free, _)) in self.subs.iter().zip(&self.local_schur) {
            let mut vecs: Vec<Vec<f64>> = Vec::new();
            for mode in &self.edge_modes {
                let mut v: Vec<f64> = free
                    .iter()
                    .map(|&i| mode[self.gpos[&sub.boundary[i]]])
                    .collect();
                // Gram–Schmidt against accepted vectors.
                for b in &vecs {
                    let proj: f64 = v.iter().zip(b).map(|(a, c)| a * c).sum();
                    for (vi, bi) in v.iter_mut().zip(b) {
                        *vi -= proj * bi;
                    }
                }
                let norm: f64 = det::sqrt(v.iter().map(|x| x * x).sum());
                if norm > 1e-12 {
                    for vi in &mut v {
                        *vi /= norm;
                    }
                    vecs.push(v);
                }
            }
            out.push(vecs);
        }
        self.local_deflate = out;
    }

    /// Assemble + factor the coarse Galerkin matrix once.
    fn factor_coarse(&mut self) {
        let n = self.gamma.len();
        let mut basis: Vec<Vec<f64>> =
            Vec::with_capacity(self.corners.len() + self.edge_modes.len());
        for &c in &self.corners {
            let mut v = vec![0.0f64; n];
            v[self.gpos[&c]] = 1.0;
            basis.push(v);
        }
        basis.extend(self.edge_modes.iter().cloned());
        if basis.is_empty() {
            return;
        }
        let nc = basis.len();
        let sv: Vec<Vec<f64>> = basis.iter().map(|v| self.schur_apply(v)).collect();
        let mut sc = vec![vec![0.0f64; nc]; nc];
        for i in 0..nc {
            for j in 0..nc {
                sc[i][j] = basis[i].iter().zip(&sv[j]).map(|(a, b)| a * b).sum();
            }
            sc[i][i] += 1e-12;
        }
        cholesky_factor(&mut sc);
        self.coarse_basis = basis;
        self.coarse_l = sc;
    }

    /// The interface Schur operator `S x_Γ` via subdomain interior
    /// solves (never assembled).
    #[must_use]
    pub fn schur_apply(&self, x_gamma: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0f64; x_gamma.len()];
        for sub in &self.subs {
            let nb = sub.boundary.len();
            let ni = sub.interior.len();
            let xb: Vec<f64> = sub.boundary.iter().map(|g| x_gamma[self.gpos[g]]).collect();
            // K_BB x_B
            let mut yb = vec![0.0f64; nb];
            for i in 0..nb {
                for j in 0..nb {
                    yb[i] += sub.k_bb[i][j] * xb[j];
                }
            }
            if ni > 0 {
                // w = K_II^{-1} K_IB x_B ; y_B -= K_IB^T w
                let mut rhs = vec![0.0f64; ni];
                for i in 0..ni {
                    for j in 0..nb {
                        rhs[i] += sub.k_ib[i][j] * xb[j];
                    }
                }
                let w = cholesky_solve(&sub.l_ii, &rhs);
                for j in 0..nb {
                    let mut acc = 0.0f64;
                    for i in 0..ni {
                        acc += sub.k_ib[i][j] * w[i];
                    }
                    yb[j] -= acc;
                }
            }
            for (i, &g) in sub.boundary.iter().enumerate() {
                out[self.gpos[&g]] += yb[i];
            }
        }
        out
    }

    /// The BDDC preconditioner application (corners + optional
    /// sheaf-edge coarse space; weighted local corrections). All
    /// factorizations were computed once at construction.
    #[must_use]
    pub fn precondition(&self, r: &[f64]) -> Vec<f64> {
        // Weighted residual.
        let rw: Vec<f64> = r.iter().zip(&self.weight).map(|(a, w)| a * w).collect();
        let mut z = vec![0.0f64; r.len()];
        // LOCAL corrections from the pre-factored free-boundary blocks,
        // DEFLATED against the coarse edge modes (P S_loc^{-1} P r).
        for ((sub, (free, l_loc)), deflate) in self
            .subs
            .iter()
            .zip(&self.local_schur)
            .zip(&self.local_deflate)
        {
            if free.is_empty() {
                continue;
            }
            let mut rhs: Vec<f64> = free
                .iter()
                .map(|&i| rw[self.gpos[&sub.boundary[i]]])
                .collect();
            for v in deflate {
                let proj: f64 = rhs.iter().zip(v).map(|(a, b)| a * b).sum();
                for (ri, vi) in rhs.iter_mut().zip(v) {
                    *ri -= proj * vi;
                }
            }
            let mut sol = cholesky_solve(l_loc, &rhs);
            for v in deflate {
                let proj: f64 = sol.iter().zip(v).map(|(a, b)| a * b).sum();
                for (si, vi) in sol.iter_mut().zip(v) {
                    *si -= proj * vi;
                }
            }
            for (a, &fa) in free.iter().enumerate() {
                z[self.gpos[&sub.boundary[fa]]] += sol[a];
            }
        }
        // COARSE correction from the pre-factored Galerkin matrix.
        if !self.coarse_basis.is_empty() {
            let rhs: Vec<f64> = self
                .coarse_basis
                .iter()
                .map(|v| v.iter().zip(&rw).map(|(a, b)| a * b).sum())
                .collect();
            let coef = cholesky_solve(&self.coarse_l, &rhs);
            for (k, v) in self.coarse_basis.iter().enumerate() {
                for (zi, vi) in z.iter_mut().zip(v) {
                    *zi += coef[k] * vi;
                }
            }
        }
        // Re-weight the output (the transpose of the input weighting).
        z.iter().zip(&self.weight).map(|(a, w)| a * w).collect()
    }

    /// Preconditioned CG on the Schur system with a Lanczos condition
    /// estimate. Returns (iterations, kappa_estimate).
    #[must_use]
    pub fn solve_cg(&self, b_gamma: &[f64], tol: f64, max_iter: usize) -> (usize, f64) {
        let n = b_gamma.len();
        let mut x = vec![0.0f64; n];
        let mut r = b_gamma.to_vec();
        let mut z = self.precondition(&r);
        let mut p = z.clone();
        let mut rz: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
        let b_norm: f64 = det::sqrt(b_gamma.iter().map(|v| v * v).sum());
        let mut alphas: Vec<f64> = Vec::new();
        let mut betas: Vec<f64> = Vec::new();
        let mut iters = 0usize;
        for _ in 0..max_iter {
            let sp = self.schur_apply(&p);
            let pap: f64 = p.iter().zip(&sp).map(|(a, b)| a * b).sum();
            let alpha = rz / pap;
            alphas.push(alpha);
            for i in 0..n {
                x[i] += alpha * p[i];
                r[i] -= alpha * sp[i];
            }
            iters += 1;
            let rn: f64 = det::sqrt(r.iter().map(|v| v * v).sum());
            if rn <= tol * b_norm {
                break;
            }
            z = self.precondition(&r);
            let rz_new: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
            let beta = rz_new / rz;
            betas.push(beta);
            rz = rz_new;
            for i in 0..n {
                p[i] = z[i] + beta * p[i];
            }
        }
        (iters, lanczos_kappa(&alphas, &betas))
    }

    /// Gamma dimension (for reports).
    #[must_use]
    pub fn gamma_len(&self) -> usize {
        self.gamma.len()
    }

    /// Coarse dimension (corners + edge modes).
    #[must_use]
    pub fn coarse_dim(&self) -> usize {
        self.corners.len() + self.edge_modes.len()
    }

    /// The CCD-locality metric for a partition mapping subdomains to
    /// `ccds` islands (row-major blocks): the fraction of interface
    /// NODES whose two subdomains land on the same island — the
    /// topological locality win of CCD-aligned partitioning.
    #[must_use]
    pub fn ccd_locality(&self, ccds: usize) -> f64 {
        let np = self.decomp.np();
        let island = |si: usize, sj: usize| -> usize {
            let per = self.decomp.s.div_ceil(ccds);
            (sj / per) * ccds + (si / per).min(ccds - 1)
        };
        let mut shared = 0usize;
        let mut total = 0usize;
        for &g in &self.gamma {
            let (x, y) = (g % np, g / np);
            let t = self.decomp.touching_subdomains(x, y);
            if t.len() == 2 {
                total += 1;
                if island(t[0].0, t[0].1) == island(t[1].0, t[1].1) {
                    shared += 1;
                }
            }
        }
        if total == 0 {
            1.0
        } else {
            #[allow(clippy::cast_precision_loss)]
            {
                shared as f64 / total as f64
            }
        }
    }
}

/// Condition estimate from CG's Lanczos tridiagonal (standard
/// coefficients-to-eigenvalues route, power/inverse-free).
fn lanczos_kappa(alphas: &[f64], betas: &[f64]) -> f64 {
    let k = alphas.len();
    if k == 0 {
        return 1.0;
    }
    // Tridiagonal entries (Golub–Van Loan): d_i, e_i from alpha/beta.
    let mut d = vec![0.0f64; k];
    let mut e = vec![0.0f64; k.saturating_sub(1)];
    for i in 0..k {
        d[i] = 1.0 / alphas[i];
        if i > 0 {
            d[i] += betas[i - 1] / alphas[i - 1];
        }
        if i + 1 < k {
            e[i] = det::sqrt(betas[i]) / alphas[i];
        }
    }
    // Eigenvalue extremes by bisection on the Sturm sequence.
    let radius: f64 = d
        .iter()
        .enumerate()
        .map(|(i, &di)| {
            let mut r = di.abs();
            if i > 0 {
                r += e[i - 1].abs();
            }
            if i < e.len() {
                r += e[i].abs();
            }
            r
        })
        .fold(0.0, f64::max);
    let count_below = |lam: f64| -> usize {
        let mut count = 0usize;
        let mut q = d[0] - lam;
        if q < 0.0 {
            count += 1;
        }
        for i in 1..k {
            let e2 = e[i - 1] * e[i - 1];
            q = d[i] - lam - e2 / if q.abs() < 1e-300 { 1e-300 } else { q };
            if q < 0.0 {
                count += 1;
            }
        }
        count
    };
    let bisect = |target: usize| -> f64 {
        let (mut lo, mut hi) = (-radius, radius);
        for _ in 0..80 {
            let mid = f64::midpoint(lo, hi);
            if count_below(mid) > target {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        f64::midpoint(lo, hi)
    };
    let lam_min = bisect(0).max(1e-300);
    let lam_max = bisect(k - 1);
    (lam_max / lam_min).max(1.0)
}

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
