//! The CutFEM discretization: Q1 elements on the ACTIVE cells (inside
//! + cut) of a quadtree background grid, embedded Dirichlet conditions
//! by symmetric Nitsche, and ghost-penalty stabilization on the faces
//! of cut cells.
//!
//! - Nitsche: a(u,v) += −∫_Γ ∂ₙu v − ∫_Γ u ∂ₙv + (β/h)∫_Γ u v, with
//!   the data terms mirrored into the load. The penalty scaling β/h is
//!   tied to the CERTIFIED cut geometry: h is the background cell size
//!   of the (certified-Cut) cell carrying the interface piece, and the
//!   ghost penalty is what makes a moderate constant β sound
//!   INDEPENDENT of how the interface cuts the cell (the conditioning
//!   battery measures exactly this).
//! - Ghost penalty: γ_g Σ_F h ∫_F [∂ₙu][∂ₙv] over faces between two
//!   active equal-level cells where at least one is cut — restoring
//!   λ_min against arbitrarily small cuts. The equal-level requirement
//!   is a build-time contract ([`crate::CutFemError::CutBandNotUniform`]);
//!   [`Quadtree::refine_toward_interface`] establishes it.
//! - Constraints (hanging, aggregation, strong outer Dirichlet) are
//!   eliminated at scatter time through affine node expansions, so the
//!   assembled system is SPD on the free DOFs and fs-solver CG applies
//!   unmodified.

use crate::CutFemError;
use crate::agg::{self, AggPolicy};
use crate::grid::{CellKey, NodeKey, Quadtree};
use crate::quad::{CutRules, cut_cell_rules, tensor_gauss};
use crate::sdf::CutSdf;
use fs_solver::krylov::CgState;
use fs_solver::op::CsrOp;
use fs_sparse::precond::IdentityPrecond;
use fs_sparse::{Coo, Csr};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

/// Certified cell classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CellClass {
    /// Certified φ < 0 over the whole cell.
    Inside,
    /// Certified φ > 0 over the whole cell (or a cut cell whose bulk
    /// rule vanished — logged as dropped).
    Outside,
    /// The certified enclosure straddles zero.
    Cut,
}

/// Discretization parameters.
#[derive(Debug, Clone, Copy)]
pub struct FemParams {
    /// Nitsche penalty constant β (the applied penalty is β/h).
    pub nitsche_beta: f64,
    /// Ghost-penalty constant γ_g; 0 disables ghost penalty.
    pub ghost_gamma: f64,
    /// Cut-quadrature subdivision depth.
    pub quad_depth: u32,
    /// Aggregated-element fallback (None disables).
    pub agg: Option<AggPolicy>,
    /// Constrain active nodes on ∂[0,1]² strongly to the Dirichlet
    /// data (for domains that reach the background boundary).
    pub strong_outer: bool,
    /// CG relative-residual target.
    pub solver_tol: f64,
    /// CG iteration cap.
    pub solver_max_iters: usize,
}

impl Default for FemParams {
    fn default() -> Self {
        FemParams {
            nitsche_beta: 10.0,
            ghost_gamma: 0.5,
            quad_depth: 3,
            agg: None,
            strong_outer: false,
            solver_tol: 1e-12,
            solver_max_iters: 60_000,
        }
    }
}

/// Build statistics (conformance logs).
#[derive(Debug, Clone, Copy, Default)]
pub struct BuildStats {
    /// Certified-inside cells.
    pub inside: usize,
    /// Certified-outside cells.
    pub outside: usize,
    /// Cut cells kept.
    pub cut: usize,
    /// Cut-classified cells dropped for vanishing bulk support.
    pub dropped: usize,
    /// Hanging-node constraints.
    pub hanging: usize,
    /// Aggregation constraints.
    pub aggregated: usize,
    /// Ghost faces assembled.
    pub ghost_faces: usize,
}

impl BuildStats {
    /// One JSON object (ledger-style log row).
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::new();
        let _ = write!(
            s,
            "{{\"inside\":{},\"outside\":{},\"cut\":{},\"dropped\":{},\
             \"hanging\":{},\"aggregated\":{},\"ghost_faces\":{}}}",
            self.inside,
            self.outside,
            self.cut,
            self.dropped,
            self.hanging,
            self.aggregated,
            self.ghost_faces
        );
        s
    }
}

enum Con {
    /// Affine combination of other nodes (hanging, aggregation).
    Affine(Vec<(NodeKey, f64)>),
    /// Strong Dirichlet: value = data at the node position.
    Strong,
}

/// A node's resolution into free DOFs and strong data samples.
#[derive(Debug, Clone, Default)]
struct Expansion {
    free: Vec<(usize, f64)>,
    strong: Vec<(NodeKey, f64)>,
}

/// The discrete solution.
#[derive(Debug, Clone)]
pub struct Solution {
    /// Free-DOF coefficients.
    pub free: Vec<f64>,
    /// Values at every mesh node (constraints applied).
    pub nodal: BTreeMap<NodeKey, f64>,
    /// CG iterations.
    pub iters: usize,
    /// Final relative residual.
    pub rel_residual: f64,
}

/// A built CutFEM space over a background quadtree.
pub struct Space<'g> {
    grid: &'g Quadtree,
    params: FemParams,
    class: BTreeMap<CellKey, CellClass>,
    rules: BTreeMap<CellKey, CutRules>,
    frac: BTreeMap<CellKey, f64>,
    active: BTreeSet<CellKey>,
    free: BTreeMap<NodeKey, usize>,
    expansions: BTreeMap<NodeKey, Expansion>,
    ghost_faces: Vec<(CellKey, CellKey, u8)>,
    agg_log: Vec<String>,
    stats: BuildStats,
}

impl<'g> Space<'g> {
    /// Classify, build quadrature, resolve constraints, enumerate
    /// ghost faces.
    ///
    /// # Errors
    /// Teaching errors: [`CutFemError::EmptyDomain`],
    /// [`CutFemError::CutBandNotUniform`],
    /// [`CutFemError::AggregationNoAnchor`],
    /// [`CutFemError::ConstraintCycle`].
    #[allow(clippy::too_many_lines)] // one linear build pipeline; splitting would smear the invariants
    pub fn build(
        grid: &'g Quadtree,
        sdf: &dyn CutSdf,
        params: FemParams,
    ) -> Result<Space<'g>, CutFemError> {
        let mut class = BTreeMap::new();
        let mut rules = BTreeMap::new();
        let mut frac = BTreeMap::new();
        let mut stats = BuildStats::default();
        for c in grid.leaves() {
            let (lo, hi) = grid.rect(c);
            let iv = sdf.enclose(lo, hi);
            let cl = if iv.hi() < 0.0 {
                CellClass::Inside
            } else if iv.lo() > 0.0 {
                CellClass::Outside
            } else {
                let r = cut_cell_rules(sdf, lo, hi, params.quad_depth);
                let area = (hi[0] - lo[0]) * (hi[1] - lo[1]);
                let w: f64 = r.bulk.iter().map(|&(_, w)| w).sum();
                if w < 1e-12 * area {
                    stats.dropped += 1;
                    CellClass::Outside
                } else {
                    frac.insert(c, w / area);
                    rules.insert(c, r);
                    CellClass::Cut
                }
            };
            match cl {
                CellClass::Inside => stats.inside += 1,
                CellClass::Outside => stats.outside += 1,
                CellClass::Cut => stats.cut += 1,
            }
            class.insert(c, cl);
        }
        let active: BTreeSet<CellKey> = class
            .iter()
            .filter(|&(_, &cl)| cl != CellClass::Outside)
            .map(|(&c, _)| c)
            .collect();
        if active.is_empty() {
            return Err(CutFemError::EmptyDomain);
        }
        let mut nodes: BTreeSet<NodeKey> = BTreeSet::new();
        for &c in &active {
            for n in grid.corner_nodes(c) {
                nodes.insert(n);
            }
        }
        // Constraints: hanging first (kinematic continuity), then
        // aggregation, then strong outer Dirichlet.
        let mut cons: BTreeMap<NodeKey, Con> = BTreeMap::new();
        for (m, ends) in grid.hanging_constraints(&active, &nodes) {
            cons.insert(m, Con::Affine(ends.to_vec()));
            stats.hanging += 1;
        }
        let mut agg_log = Vec::new();
        if let Some(policy) = params.agg {
            let constrained: BTreeSet<NodeKey> = cons.keys().copied().collect();
            let outcome = agg::aggregate(grid, &class, &frac, &active, &constrained, policy)?;
            for (n, terms) in outcome.constraints {
                cons.insert(n, Con::Affine(terms));
                stats.aggregated += 1;
            }
            agg_log = outcome.log;
        }
        if params.strong_outer {
            let ext = grid.node_extent();
            for &n in &nodes {
                if (n.0 == 0 || n.0 == ext || n.1 == 0 || n.1 == ext)
                    && !cons.contains_key(&n)
                {
                    cons.insert(n, Con::Strong);
                }
            }
        }
        let mut free: BTreeMap<NodeKey, usize> = BTreeMap::new();
        for &n in &nodes {
            if !cons.contains_key(&n) {
                let id = free.len();
                free.insert(n, id);
            }
        }
        let mut expansions: BTreeMap<NodeKey, Expansion> = BTreeMap::new();
        for &n in &nodes {
            let mut stack = BTreeSet::new();
            expand_node(n, &cons, &free, &mut expansions, &mut stack)?;
        }
        // Ghost faces: between equal-level active cells, at least one
        // cut. A differently-leveled active neighbor of a cut cell is
        // a build refusal — the interface band must be uniform.
        let mut ghost_faces = Vec::new();
        let mut seen: BTreeSet<(CellKey, CellKey)> = BTreeSet::new();
        for (&c, &cl) in &class {
            if cl != CellClass::Cut {
                continue;
            }
            for dir in 0..4u8 {
                let Some(nb) = grid.covering_neighbor(c, dir) else {
                    continue;
                };
                if !active.contains(&nb) {
                    continue;
                }
                if nb.0 != c.0 {
                    return Err(CutFemError::CutBandNotUniform {
                        cell: c,
                        neighbor: nb,
                    });
                }
                let key = if c < nb { (c, nb) } else { (nb, c) };
                let axis = u8::from(dir >= 2);
                if seen.insert(key) {
                    ghost_faces.push((key.0, key.1, axis));
                }
            }
        }
        stats.ghost_faces = ghost_faces.len();
        Ok(Space {
            grid,
            params,
            class,
            rules,
            frac,
            active,
            free,
            expansions,
            ghost_faces,
            agg_log,
            stats,
        })
    }

    /// Free-DOF count.
    #[must_use]
    pub fn dof_count(&self) -> usize {
        self.free.len()
    }

    /// Mesh-node count (free + constrained).
    #[must_use]
    pub fn node_count(&self) -> usize {
        self.expansions.len()
    }

    /// Build statistics.
    #[must_use]
    pub fn stats(&self) -> BuildStats {
        self.stats
    }

    /// Aggregation policy log rows (JSON objects).
    #[must_use]
    pub fn agg_log(&self) -> &[String] {
        &self.agg_log
    }

    /// A cell's classification (`None` for unknown keys).
    #[must_use]
    pub fn class_of(&self, c: CellKey) -> Option<CellClass> {
        self.class.get(&c).copied()
    }

    /// A cut cell's certified inside-area fraction.
    #[must_use]
    pub fn fraction_of(&self, c: CellKey) -> Option<f64> {
        self.frac.get(&c).copied()
    }

    /// Assemble the stiffness matrix and load vector for
    /// `−Δu = f` in Ω, `u = g` on Γ (Nitsche) and on the strong outer
    /// boundary.
    #[must_use]
    #[allow(clippy::too_many_lines)] // bulk + Nitsche + ghost in one deterministic sweep
    pub fn assemble(
        &self,
        f: &dyn Fn(f64, f64) -> f64,
        g: &dyn Fn(f64, f64) -> f64,
    ) -> (Csr, Vec<f64>) {
        let nf = self.free.len();
        let mut coo = Coo::new(nf, nf);
        let mut rhs = vec![0.0f64; nf];
        let consts: BTreeMap<NodeKey, f64> = self
            .expansions
            .iter()
            .map(|(n, e)| {
                let c: f64 = e
                    .strong
                    .iter()
                    .map(|&(sn, w)| {
                        let p = self.grid.node_pos(sn);
                        w * g(p[0], p[1])
                    })
                    .sum();
                (*n, c)
            })
            .collect();
        for &c in &self.active {
            let (lo, hi) = self.grid.rect(c);
            let corners = self.grid.corner_nodes(c);
            let h = self.grid.cell_h(c);
            let cl = self.class[&c];
            let mut k = [[0.0f64; 4]; 4];
            let mut fl = [0.0f64; 4];
            let inside_rule;
            let bulk: &[([f64; 2], f64)] = if cl == CellClass::Inside {
                inside_rule = {
                    let mut v = Vec::with_capacity(9);
                    tensor_gauss(lo, hi, &mut v);
                    v
                };
                &inside_rule
            } else {
                &self.rules[&c].bulk
            };
            for &(p, w) in bulk {
                let (nv, gr) = q1(lo, hi, p);
                let fv = f(p[0], p[1]);
                for a in 0..4 {
                    for b in 0..4 {
                        k[a][b] += w * (gr[a][0] * gr[b][0] + gr[a][1] * gr[b][1]);
                    }
                    fl[a] += w * fv * nv[a];
                }
            }
            if cl == CellClass::Cut {
                let pen = self.params.nitsche_beta / h;
                for &(p, w, nrm) in &self.rules[&c].iface {
                    let (nv, gr) = q1(lo, hi, p);
                    let gval = g(p[0], p[1]);
                    let dn = [
                        gr[0][0] * nrm[0] + gr[0][1] * nrm[1],
                        gr[1][0] * nrm[0] + gr[1][1] * nrm[1],
                        gr[2][0] * nrm[0] + gr[2][1] * nrm[1],
                        gr[3][0] * nrm[0] + gr[3][1] * nrm[1],
                    ];
                    for a in 0..4 {
                        for b in 0..4 {
                            k[a][b] +=
                                w * (-(dn[a] * nv[b]) - nv[a] * dn[b] + pen * nv[a] * nv[b]);
                        }
                        fl[a] += w * (-dn[a] * gval + pen * nv[a] * gval);
                    }
                }
            }
            self.scatter(&mut coo, &mut rhs, &corners, &k, &fl, &consts);
        }
        if self.params.ghost_gamma > 0.0 {
            for &(ca, cb, axis) in &self.ghost_faces {
                self.ghost_face(&mut coo, &mut rhs, ca, cb, axis, &consts);
            }
        }
        (coo.assemble(), rhs)
    }

    fn scatter(
        &self,
        coo: &mut Coo,
        rhs: &mut [f64],
        corners: &[NodeKey; 4],
        k: &[[f64; 4]; 4],
        fl: &[f64; 4],
        consts: &BTreeMap<NodeKey, f64>,
    ) {
        for a in 0..4 {
            let ea = &self.expansions[&corners[a]];
            for &(ia, wa) in &ea.free {
                rhs[ia] += wa * fl[a];
                for b in 0..4 {
                    let kab = k[a][b];
                    if kab == 0.0 {
                        continue;
                    }
                    for &(ib, wb) in &self.expansions[&corners[b]].free {
                        coo.push(ia, ib, wa * kab * wb);
                    }
                    let cb = consts[&corners[b]];
                    if cb != 0.0 {
                        rhs[ia] -= wa * kab * cb;
                    }
                }
            }
        }
    }

    /// One ghost face: γ_g · h · ∫_F [∂ₙu][∂ₙv] ds with 2-point Gauss
    /// (exact: the normal-derivative jump of Q1 is linear along the
    /// face).
    fn ghost_face(
        &self,
        coo: &mut Coo,
        rhs: &mut [f64],
        ca: CellKey,
        cb: CellKey,
        axis: u8,
        consts: &BTreeMap<NodeKey, f64>,
    ) {
        let (lo_a, hi_a) = self.grid.rect(ca);
        let (lo_b, hi_b) = self.grid.rect(cb);
        let h = self.grid.cell_h(ca);
        // Key order puts `ca` on the low side of the shared face.
        let (t0, t1) = if axis == 0 {
            (lo_a[1], hi_a[1])
        } else {
            (lo_a[0], hi_a[0])
        };
        let nrm = if axis == 0 { [1.0, 0.0] } else { [0.0, 1.0] };
        let xf = if axis == 0 { hi_a[0] } else { hi_a[1] };
        let ca_corners = self.grid.corner_nodes(ca);
        let cb_corners = self.grid.corner_nodes(cb);
        let gpt = 0.5 / 3.0f64.sqrt();
        let mut jump: BTreeMap<NodeKey, [f64; 2]> = BTreeMap::new();
        let wq = 0.5 * (t1 - t0);
        for (qi, t) in [0.5 - gpt, 0.5 + gpt].into_iter().enumerate() {
            let tv = t0 + t * (t1 - t0);
            let p = if axis == 0 { [xf, tv] } else { [tv, xf] };
            let (_, gra) = q1(lo_a, hi_a, p);
            let (_, grb) = q1(lo_b, hi_b, p);
            for a in 0..4 {
                jump.entry(ca_corners[a]).or_default()[qi] +=
                    gra[a][0] * nrm[0] + gra[a][1] * nrm[1];
            }
            for b in 0..4 {
                jump.entry(cb_corners[b]).or_default()[qi] -=
                    grb[b][0] * nrm[0] + grb[b][1] * nrm[1];
            }
        }
        let scale = self.params.ghost_gamma * h * wq;
        let entries: Vec<(NodeKey, [f64; 2])> = jump.into_iter().collect();
        for (na, ja) in &entries {
            let ea = &self.expansions[na];
            for (nb, jb) in &entries {
                let v = scale * (ja[0] * jb[0] + ja[1] * jb[1]);
                if v == 0.0 {
                    continue;
                }
                for &(ia, wa) in &ea.free {
                    for &(ib, wb) in &self.expansions[nb].free {
                        coo.push(ia, ib, wa * v * wb);
                    }
                    let cb = consts[nb];
                    if cb != 0.0 {
                        rhs[ia] -= wa * v * cb;
                    }
                }
            }
        }
    }

    /// Assemble and CG-solve; returns the solution with nodal values
    /// expanded through every constraint.
    ///
    /// # Errors
    /// [`CutFemError::SolveNotConverged`] if the residual gate (1e-8)
    /// is missed.
    pub fn solve(
        &self,
        f: &dyn Fn(f64, f64) -> f64,
        g: &dyn Fn(f64, f64) -> f64,
    ) -> Result<Solution, CutFemError> {
        let (a, b) = self.assemble(f, g);
        let op = CsrOp::symmetric(a);
        let m = IdentityPrecond;
        let mut st = CgState::new(&op, &m, &b);
        let _ = st.run(&op, &m, self.params.solver_tol, self.params.solver_max_iters);
        let rr = st.rel_residual();
        if !rr.is_finite() || rr > 1e-8 {
            return Err(CutFemError::SolveNotConverged {
                iters: st.iters,
                rel_residual: rr,
            });
        }
        let nodal = self.nodal_values(&st.x, g);
        Ok(Solution {
            free: st.x,
            nodal,
            iters: st.iters,
            rel_residual: rr,
        })
    }

    /// Expand free-DOF coefficients to values at every mesh node.
    #[must_use]
    pub fn nodal_values(
        &self,
        free_x: &[f64],
        g: &dyn Fn(f64, f64) -> f64,
    ) -> BTreeMap<NodeKey, f64> {
        self.expansions
            .iter()
            .map(|(n, e)| {
                let mut v: f64 = e.free.iter().map(|&(i, w)| w * free_x[i]).sum();
                v += e
                    .strong
                    .iter()
                    .map(|&(sn, w)| {
                        let p = self.grid.node_pos(sn);
                        w * g(p[0], p[1])
                    })
                    .sum::<f64>();
                (*n, v)
            })
            .collect()
    }

    /// L2 and H1-seminorm errors against an exact solution, integrated
    /// with one-deeper cut quadrature (so the measurement error is
    /// dominated by the discretization, not the meter).
    #[must_use]
    pub fn l2_h1_error(
        &self,
        sdf: &dyn CutSdf,
        exact: &dyn Fn(f64, f64) -> f64,
        grad_exact: &dyn Fn(f64, f64) -> [f64; 2],
        nodal: &BTreeMap<NodeKey, f64>,
    ) -> (f64, f64) {
        let mut l2 = 0.0f64;
        let mut h1 = 0.0f64;
        for &c in &self.active {
            let (lo, hi) = self.grid.rect(c);
            let corners = self.grid.corner_nodes(c);
            let vals = [
                nodal[&corners[0]],
                nodal[&corners[1]],
                nodal[&corners[2]],
                nodal[&corners[3]],
            ];
            let refined;
            let rule: &[([f64; 2], f64)] = if self.class[&c] == CellClass::Inside {
                refined = {
                    let mut v = Vec::with_capacity(9);
                    tensor_gauss(lo, hi, &mut v);
                    v
                };
                &refined
            } else {
                refined = cut_cell_rules(sdf, lo, hi, self.params.quad_depth + 1).bulk;
                &refined
            };
            for &(p, w) in rule {
                let (nv, gr) = q1(lo, hi, p);
                let mut uh = 0.0;
                let mut guh = [0.0, 0.0];
                for a in 0..4 {
                    uh += nv[a] * vals[a];
                    guh[0] += gr[a][0] * vals[a];
                    guh[1] += gr[a][1] * vals[a];
                }
                let e = exact(p[0], p[1]) - uh;
                let ge = grad_exact(p[0], p[1]);
                let gex = ge[0] - guh[0];
                let gey = ge[1] - guh[1];
                l2 += w * e * e;
                h1 += w * (gex * gex + gey * gey);
            }
        }
        (l2.max(0.0).sqrt(), h1.max(0.0).sqrt())
    }
}

fn expand_node(
    n: NodeKey,
    cons: &BTreeMap<NodeKey, Con>,
    free: &BTreeMap<NodeKey, usize>,
    memo: &mut BTreeMap<NodeKey, Expansion>,
    stack: &mut BTreeSet<NodeKey>,
) -> Result<(), CutFemError> {
    if memo.contains_key(&n) {
        return Ok(());
    }
    if !stack.insert(n) {
        return Err(CutFemError::ConstraintCycle { node: n });
    }
    let e = match cons.get(&n) {
        None => Expansion {
            free: vec![(free[&n], 1.0)],
            strong: Vec::new(),
        },
        Some(Con::Strong) => Expansion {
            free: Vec::new(),
            strong: vec![(n, 1.0)],
        },
        Some(Con::Affine(terms)) => {
            let terms = terms.clone();
            let mut fr: BTreeMap<usize, f64> = BTreeMap::new();
            let mut st: BTreeMap<NodeKey, f64> = BTreeMap::new();
            for (child, w) in terms {
                expand_node(child, cons, free, memo, stack)?;
                let ce = memo.get(&child).expect("just expanded").clone();
                for (id, cw) in ce.free {
                    *fr.entry(id).or_insert(0.0) += w * cw;
                }
                for (sn, sw) in ce.strong {
                    *st.entry(sn).or_insert(0.0) += w * sw;
                }
            }
            Expansion {
                free: fr.into_iter().collect(),
                strong: st.into_iter().collect(),
            }
        }
    };
    stack.remove(&n);
    memo.insert(n, e);
    Ok(())
}

/// Q1 shape values and gradients on an axis-aligned cell, corner order
/// (0,0), (1,0), (1,1), (0,1) — matching [`Quadtree::corner_nodes`].
fn q1(lo: [f64; 2], hi: [f64; 2], p: [f64; 2]) -> ([f64; 4], [[f64; 2]; 4]) {
    let hx = hi[0] - lo[0];
    let hy = hi[1] - lo[1];
    let xi = (p[0] - lo[0]) / hx;
    let et = (p[1] - lo[1]) / hy;
    let n = [
        (1.0 - xi) * (1.0 - et),
        xi * (1.0 - et),
        xi * et,
        (1.0 - xi) * et,
    ];
    let g = [
        [-(1.0 - et) / hx, -(1.0 - xi) / hy],
        [(1.0 - et) / hx, -xi / hy],
        [et / hx, xi / hy],
        [-et / hx, (1.0 - xi) / hy],
    ];
    (n, g)
}
