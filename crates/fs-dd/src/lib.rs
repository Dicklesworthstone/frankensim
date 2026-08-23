//! fs-dd — DOMAIN DECOMPOSITION (plan §8.9, bead tfz.11; [M] — behind
//! the `bddc` / `sheaf-coarse` features until Gauntlet evidence
//! promotes them). BDDC substructuring on the 2-D model problem:
//! corners-primal coarse spaces with coefficient-scaled interface averaging,
//! the SHEAF-HARMONIC edge enrichment (the same cellular-sheaf machinery as
//! watertightness — Bet 11 earning its keep twice), subdomain-aligned
//! coefficient-jump evidence, and CCD-aligned partitioning metrics. Layer L3.
//!
//! The model problem is the variable-coefficient 5-point Laplacian on
//! the unit square (spectrally equivalent to P1 FEM — the setting BDDC
//! theory addresses): `n = s·m` cells per side, `s²` square subdomains
//! of `m²` cells. The solver is CG on the INTERFACE Schur system with
//! the BDDC preconditioner. Condition diagnostics come only from validated
//! CG/Lanczos coefficients; conditioning claims require explicit convergence.
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

    /// Harmonic mean of the adjacent cell coefficients flanking an edge
    /// between nodes (the standard variable-coefficient transmissibility).
    fn edge_coeff(&self, ax: usize, ay: usize, bx: usize, by: usize) -> f64 {
        let n = self.n();
        let mut coefficients = [0.0f64; 2];
        let mut count = 0usize;
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
                let coefficient = self.rho[cy * n + cx];
                assert!(
                    coefficient.is_finite() && coefficient > 0.0,
                    "cell coefficients must be finite and positive"
                );
                coefficients[count] = coefficient;
                count += 1;
            }
        }
        // Identical coefficients make the harmonic mean bitwise-equal to either
        // input; an epsilon would silently corrupt exact uniform edges.
        #[allow(clippy::float_cmp)]
        match count {
            0 => 1.0,
            1 => coefficients[0],
            2 if coefficients[0] == coefficients[1] => coefficients[0],
            2 => {
                let lo = coefficients[0].min(coefficients[1]);
                let hi = coefficients[0].max(coefficients[1]);
                // Algebraically 2ab/(a+b), but scaled by the larger input so
                // positive subnormal coefficients do not vanish through
                // overflowing reciprocals and finite large inputs do not
                // overflow through the product.
                lo / (0.5 * (1.0 + lo / hi))
            }
            _ => unreachable!("a lattice edge has at most two flanking cells"),
        }
    }

    /// Number of closed subdomains containing both endpoints of a lattice
    /// edge. An edge along an internal interface has multiplicity two; all
    /// other edges have multiplicity one.
    fn edge_multiplicity(&self, ax: usize, ay: usize, bx: usize, by: usize) -> usize {
        debug_assert_eq!(ax.abs_diff(bx) + ay.abs_diff(by), 1);
        let interface_coordinate = if ay == by { ay } else { ax };
        if self.m > 0
            && interface_coordinate > 0
            && interface_coordinate < self.n()
            && interface_coordinate % self.m == 0
        {
            2
        } else {
            1
        }
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
    /// Standard nodal rho-scaling on `boundary`, normalized across every
    /// subdomain sharing the corresponding global node.
    boundary_weight: Vec<f64>,
    /// Factored K_II (interior-interior).
    l_ii: Vec<Vec<f64>>,
    /// K_IB (interior × boundary).
    k_ib: Vec<Vec<f64>>,
    /// K_BB (boundary × boundary).
    k_bb: Vec<Vec<f64>>,
}

/// A preconditioned-CG solve report.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CgReport {
    /// Why the iteration stopped.
    pub termination: CgTermination,
    /// Number of completed CG iterations.
    pub iterations: usize,
    /// Last finite residual norm divided by the normalized RHS norm.
    ///
    /// On numerical breakdown this is the last confirmed finite value before
    /// the rejected recurrence step, not a value synthesized from poisoned
    /// arithmetic.
    pub relative_residual: f64,
    /// Lanczos/Ritz condition estimate when at least one valid iteration
    /// supplied a finite positive-definite tridiagonal. A one-dimensional
    /// projection transparently reports `kappa = 1` with
    /// `krylov_dimension = 1`; zero steps have no spectral evidence.
    ///
    /// This diagnostic is not a certified enclosure and is not claim-worthy
    /// unless [`CgTermination::Converged`] also holds.
    pub condition_estimate: Option<LanczosConditionEstimate>,
}

/// Extreme-Ritz condition diagnostic from a valid CG/Lanczos prefix.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LanczosConditionEstimate {
    /// Ratio `ritz_max / ritz_min`.
    pub kappa: f64,
    /// Smallest Ritz value of the normalized projected operator.
    pub ritz_min: f64,
    /// Largest Ritz value of the normalized projected operator.
    pub ritz_max: f64,
    /// Dimension of the RHS-dependent Krylov projection.
    pub krylov_dimension: usize,
}

/// Termination disposition for [`Bddc::solve_cg`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CgTermination {
    /// The admitted relative-residual tolerance was met.
    Converged,
    /// The explicit iteration budget was exhausted with a finite recurrence.
    IterationLimit,
    /// The SPD or finite-arithmetic recurrence broke down.
    Breakdown(CgBreakdown),
}

/// A fail-closed preconditioned-CG recurrence failure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CgBreakdown {
    /// `rᵀ M⁻¹ r` was non-finite.
    NonFinitePreconditionedResidual {
        /// Completed iterations when detected.
        iteration: usize,
        /// Rejected inner product.
        value: f64,
    },
    /// `rᵀ M⁻¹ r` was not strictly positive for a non-converged residual.
    NonPositivePreconditionedResidual {
        /// Completed iterations when detected.
        iteration: usize,
        /// Rejected inner product.
        value: f64,
    },
    /// `pᵀ S p` was non-finite.
    NonFiniteCurvature {
        /// Completed iterations before the attempted step.
        iteration: usize,
        /// Rejected curvature.
        value: f64,
    },
    /// `pᵀ S p` was not strictly positive, violating the SPD contract.
    NonPositiveCurvature {
        /// Completed iterations before the attempted step.
        iteration: usize,
        /// Rejected curvature.
        value: f64,
    },
    /// The CG step length was non-finite or non-positive.
    InvalidStep {
        /// Completed iterations before the attempted step.
        iteration: usize,
        /// Rejected `alpha`.
        value: f64,
    },
    /// A normalized solution component became non-finite during the update.
    NonFiniteSolution {
        /// Completed iterations before the attempted step.
        iteration: usize,
        /// Rejected solution component index.
        index: usize,
        /// Rejected component.
        value: f64,
    },
    /// A residual component became non-finite during the update.
    NonFiniteResidual {
        /// Completed iterations before the attempted step.
        iteration: usize,
        /// Rejected residual component index.
        index: usize,
        /// Rejected component.
        value: f64,
    },
    /// The residual norm or relative residual became non-finite.
    NonFiniteResidualNorm {
        /// Completed iterations including the rejected residual update.
        iteration: usize,
        /// Rejected norm or ratio.
        value: f64,
    },
    /// The CG direction recurrence coefficient was non-finite or non-positive.
    InvalidRecurrence {
        /// Completed iterations when detected.
        iteration: usize,
        /// Rejected `beta`.
        value: f64,
    },
    /// A search-direction component became non-finite.
    NonFiniteDirection {
        /// Completed iterations when detected.
        iteration: usize,
        /// Rejected direction component index.
        index: usize,
        /// Rejected component.
        value: f64,
    },
}

/// Admission failures for [`Bddc::solve_cg`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CgError {
    /// The RHS length does not match the retained interface dimension.
    DimensionMismatch {
        /// Required interface dimension.
        expected: usize,
        /// Supplied RHS length.
        actual: usize,
    },
    /// A supplied RHS component is non-finite.
    NonFiniteRhs {
        /// Rejected component index.
        index: usize,
        /// Rejected component.
        value: f64,
    },
    /// Relative tolerance must be finite and non-negative.
    InvalidTolerance {
        /// Rejected tolerance.
        tolerance: f64,
    },
}

impl core::fmt::Display for CgError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::DimensionMismatch { expected, actual } => write!(
                f,
                "CG RHS length must match the {expected}-entry interface dimension (got {actual})"
            ),
            Self::NonFiniteRhs { index, value } => {
                write!(f, "CG RHS entry {index} must be finite (got {value})")
            }
            Self::InvalidTolerance { tolerance } => write!(
                f,
                "CG relative tolerance must be finite and non-negative (got {tolerance})"
            ),
        }
    }
}

impl std::error::Error for CgError {}

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

#[derive(Clone, Copy)]
struct SubdomainBounds {
    x0: usize,
    x1: usize,
    y0: usize,
    y1: usize,
}

impl SubdomainBounds {
    fn contains(self, x: usize, y: usize) -> bool {
        x >= self.x0 && x <= self.x1 && y >= self.y0 && y <= self.y1
    }
}

/// Standard nodal rho-scaling numerator for one subdomain trace.
///
/// For a boundary basis function, rho-scaling uses the largest coefficient
/// among the cells in this subdomain that touch the node. This is exactly the
/// subdomain coefficient on the checkerboard fixture while remaining defined
/// for the public arbitrary per-cell coefficient field.
fn local_nodal_coefficient(d: &Decomposition, bounds: SubdomainBounds, x: usize, y: usize) -> f64 {
    let n = d.n();
    let mut coefficient = 0.0_f64;
    for cy in [y.checked_sub(1), Some(y)].into_iter().flatten() {
        for cx in [x.checked_sub(1), Some(x)].into_iter().flatten() {
            if cx >= bounds.x0
                && cx < bounds.x1
                && cy >= bounds.y0
                && cy < bounds.y1
                && cx < n
                && cy < n
            {
                let local = d.rho[cy * n + cx];
                assert!(
                    local.is_finite() && local > 0.0,
                    "cell coefficients must be finite and positive"
                );
                coefficient = coefficient.max(local);
            }
        }
    }
    assert!(
        coefficient.is_finite() && coefficient > 0.0,
        "each retained subdomain-boundary node must touch a local cell"
    );
    coefficient
}

fn local_stiffness(d: &Decomposition, nodes: &[usize], bounds: SubdomainBounds) -> Vec<Vec<f64>> {
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
            if nx > d.n() || ny > d.n() || !bounds.contains(nx, ny) {
                continue;
            }
            let multiplicity = d.edge_multiplicity(x, y, nx, ny);
            assert!(multiplicity > 0, "local edge must have an owning subdomain");
            #[allow(clippy::cast_precision_loss)]
            let w = d.edge_coeff(x, y, nx, ny) / multiplicity as f64;
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
                let bounds = SubdomainBounds { x0, x1, y0, y1 };
                let a = local_stiffness(&decomp, &all, bounds);
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
                let boundary_weight = boundary
                    .iter()
                    .map(|&g| local_nodal_coefficient(&decomp, bounds, g % np, g / np))
                    .collect();
                subs.push(Subdomain {
                    interior,
                    boundary,
                    boundary_weight,
                    l_ii: k_ii,
                    k_ib,
                    k_bb,
                });
            }
        }
        // Standard rho-scaling:
        //
        //   delta_i(g) = rho_i(g) / sum_{j containing g} rho_j(g).
        //
        // Form the ratio through a per-node maximum so the denominator cannot
        // overflow when every admitted coefficient is finite but very large.
        let mut max_weight = vec![0.0f64; gamma.len()];
        for sub in &subs {
            for (&g, &raw) in sub.boundary.iter().zip(&sub.boundary_weight) {
                let slot = gpos[&g];
                max_weight[slot] = max_weight[slot].max(raw);
            }
        }
        assert!(
            max_weight
                .iter()
                .all(|weight| weight.is_finite() && *weight > 0.0),
            "every gamma node must have a finite positive rho-scaling numerator"
        );
        let mut scaled_sum = vec![0.0f64; gamma.len()];
        for sub in &subs {
            for (&g, &raw) in sub.boundary.iter().zip(&sub.boundary_weight) {
                let slot = gpos[&g];
                scaled_sum[slot] += raw / max_weight[slot];
            }
        }
        assert!(
            scaled_sum.iter().all(|sum| sum.is_finite() && *sum >= 1.0),
            "every gamma node must have a finite rho-scaling denominator"
        );
        for sub in &mut subs {
            for (&g, weight) in sub.boundary.iter().zip(&mut sub.boundary_weight) {
                let slot = gpos[&g];
                *weight = (*weight / max_weight[slot]) / scaled_sum[slot];
            }
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
            let mut yb: Vec<f64> = sub
                .k_bb
                .iter()
                .map(|row| row.iter().zip(&xb).map(|(a, b)| a * b).sum())
                .collect();
            if ni > 0 {
                // w = K_II^{-1} K_IB x_B ; y_B -= K_IB^T w
                let rhs: Vec<f64> = sub
                    .k_ib
                    .iter()
                    .map(|row| row.iter().zip(&xb).map(|(a, b)| a * b).sum())
                    .collect();
                let w = cholesky_solve(&sub.l_ii, &rhs);
                for (j, y) in yb.iter_mut().enumerate() {
                    let acc: f64 = sub.k_ib.iter().zip(&w).map(|(row, wi)| row[j] * wi).sum();
                    *y -= acc;
                }
            }
            let _ = nb;
            for (i, &g) in sub.boundary.iter().enumerate() {
                out[self.gpos[&g]] += yb[i];
            }
        }
        out
    }

    /// The BDDC preconditioner application: nodal rho-scaled local corrections
    /// plus the corners + optional sheaf-edge continuous coarse correction.
    /// All factorizations were computed once at construction.
    #[must_use]
    pub fn precondition(&self, r: &[f64]) -> Vec<f64> {
        let mut z = vec![0.0f64; r.len()];
        // LOCAL corrections from the pre-factored free-boundary blocks,
        // DEFLATED against the coarse edge modes. The same subdomain-local
        // diagonal scaling acts on gather and scatter:
        //
        //   sum_i R_i^T D_i P_i S_i^-1 P_i D_i R_i.
        //
        // Using the same D_i on both sides preserves symmetry.
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
                .map(|&i| sub.boundary_weight[i] * r[self.gpos[&sub.boundary[i]]])
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
                z[self.gpos[&sub.boundary[fa]]] += sub.boundary_weight[fa] * sol[a];
            }
        }
        // COARSE correction from the pre-factored global Galerkin matrix.
        // The basis is already continuous on gamma. Since the D_i form a
        // partition of unity, their gather/scatter sums cancel on this global
        // component; applying one fabricated global weight would damp it.
        if !self.coarse_basis.is_empty() {
            let rhs: Vec<f64> = self
                .coarse_basis
                .iter()
                .map(|v| v.iter().zip(r).map(|(a, b)| a * b).sum())
                .collect();
            let coef = cholesky_solve(&self.coarse_l, &rhs);
            for (k, v) in self.coarse_basis.iter().enumerate() {
                for (zi, vi) in z.iter_mut().zip(v) {
                    *zi += coef[k] * vi;
                }
            }
        }
        z
    }

    /// Preconditioned CG on the Schur system with explicit termination,
    /// residual, and RHS-dependent Lanczos/Ritz diagnostics.
    ///
    /// The finite RHS is normalized by its largest magnitude before any dot
    /// product. A zero iteration budget is a valid request: a zero RHS (or an
    /// already-admitted loose tolerance) is converged at iteration zero, while
    /// a non-converged RHS reports [`CgTermination::IterationLimit`].
    ///
    /// # Errors
    /// [`CgError`] for caller-invalid dimension, RHS, or tolerance. Numerical
    /// recurrence failures are successful reports with
    /// [`CgTermination::Breakdown`] so the last finite residual is retained.
    #[allow(clippy::too_many_lines)] // One CG driver keeps recurrence, breakdown attribution, and Lanczos bookkeeping in lockstep; splitting would let them drift apart.
    pub fn solve_cg(
        &self,
        b_gamma: &[f64],
        tol: f64,
        max_iter: usize,
    ) -> Result<CgReport, CgError> {
        if b_gamma.len() != self.gamma.len() {
            return Err(CgError::DimensionMismatch {
                expected: self.gamma.len(),
                actual: b_gamma.len(),
            });
        }
        if let Some((index, value)) = b_gamma
            .iter()
            .copied()
            .enumerate()
            .find(|(_, value)| !value.is_finite())
        {
            return Err(CgError::NonFiniteRhs { index, value });
        }
        if !tol.is_finite() || tol < 0.0 {
            return Err(CgError::InvalidTolerance { tolerance: tol });
        }
        let n = b_gamma.len();
        let b_scale = b_gamma
            .iter()
            .map(|value| value.abs())
            .fold(0.0f64, f64::max);
        if b_scale == 0.0 {
            return Ok(CgReport {
                termination: CgTermination::Converged,
                iterations: 0,
                relative_residual: 0.0,
                condition_estimate: None,
            });
        }
        // CG is homogeneous in the RHS, and this API returns diagnostics rather
        // than a solution vector. Normalize before dot products so a finite
        // nonzero RHS cannot disappear through squaring underflow or overflow.
        let scaled_b: Vec<f64> = b_gamma.iter().map(|value| value / b_scale).collect();
        let b_norm = scaled_l2_norm(&scaled_b);
        if !b_norm.is_finite() || b_norm <= 0.0 {
            return Ok(cg_breakdown_report(
                0,
                1.0,
                &[],
                &[],
                CgBreakdown::NonFiniteResidualNorm {
                    iteration: 0,
                    value: b_norm,
                },
            ));
        }
        let mut relative_residual = 1.0;
        if relative_residual <= tol {
            return Ok(CgReport {
                termination: CgTermination::Converged,
                iterations: 0,
                relative_residual,
                condition_estimate: None,
            });
        }
        if max_iter == 0 {
            return Ok(CgReport {
                termination: CgTermination::IterationLimit,
                iterations: 0,
                relative_residual,
                condition_estimate: None,
            });
        }
        let mut x = vec![0.0f64; n];
        let mut r = scaled_b.clone();
        let mut z = self.precondition(&r);
        let mut p = z.clone();
        let mut rz = dot(&r, &z);
        let mut alphas: Vec<f64> = Vec::new();
        let mut betas: Vec<f64> = Vec::new();
        let mut iters = 0usize;
        if !rz.is_finite() {
            return Ok(cg_breakdown_report(
                iters,
                relative_residual,
                &alphas,
                &betas,
                CgBreakdown::NonFinitePreconditionedResidual {
                    iteration: iters,
                    value: rz,
                },
            ));
        }
        if rz <= 0.0 {
            return Ok(cg_breakdown_report(
                iters,
                relative_residual,
                &alphas,
                &betas,
                CgBreakdown::NonPositivePreconditionedResidual {
                    iteration: iters,
                    value: rz,
                },
            ));
        }
        loop {
            let sp = self.schur_apply(&p);
            let pap = dot(&p, &sp);
            if !pap.is_finite() {
                return Ok(cg_breakdown_report(
                    iters,
                    relative_residual,
                    &alphas,
                    &betas,
                    CgBreakdown::NonFiniteCurvature {
                        iteration: iters,
                        value: pap,
                    },
                ));
            }
            if pap <= 0.0 {
                return Ok(cg_breakdown_report(
                    iters,
                    relative_residual,
                    &alphas,
                    &betas,
                    CgBreakdown::NonPositiveCurvature {
                        iteration: iters,
                        value: pap,
                    },
                ));
            }
            let alpha = rz / pap;
            if !alpha.is_finite() || alpha <= 0.0 {
                return Ok(cg_breakdown_report(
                    iters,
                    relative_residual,
                    &alphas,
                    &betas,
                    CgBreakdown::InvalidStep {
                        iteration: iters,
                        value: alpha,
                    },
                ));
            }
            for i in 0..n {
                let next_x = alpha.mul_add(p[i], x[i]);
                if !next_x.is_finite() {
                    return Ok(cg_breakdown_report(
                        iters,
                        relative_residual,
                        &alphas,
                        &betas,
                        CgBreakdown::NonFiniteSolution {
                            iteration: iters,
                            index: i,
                            value: next_x,
                        },
                    ));
                }
                x[i] = next_x;
                r[i] -= alpha * sp[i];
                if !r[i].is_finite() {
                    return Ok(cg_breakdown_report(
                        iters,
                        relative_residual,
                        &alphas,
                        &betas,
                        CgBreakdown::NonFiniteResidual {
                            iteration: iters,
                            index: i,
                            value: r[i],
                        },
                    ));
                }
            }
            iters += 1;
            let residual_norm = scaled_l2_norm(&r);
            let next_relative_residual = residual_norm / b_norm;
            if !residual_norm.is_finite() || !next_relative_residual.is_finite() {
                return Ok(cg_breakdown_report(
                    iters,
                    relative_residual,
                    &alphas,
                    &betas,
                    CgBreakdown::NonFiniteResidualNorm {
                        iteration: iters,
                        value: next_relative_residual,
                    },
                ));
            }
            relative_residual = next_relative_residual;
            alphas.push(alpha);
            let recursive_converged = if tol == 0.0 {
                residual_norm == 0.0
            } else {
                relative_residual <= tol
            };
            if recursive_converged || iters == max_iter {
                let applied = self.schur_apply(&x);
                let mut true_r = Vec::with_capacity(n);
                for (index, (rhs, image)) in scaled_b.iter().zip(&applied).enumerate() {
                    let value = rhs - image;
                    if !value.is_finite() {
                        return Ok(cg_breakdown_report(
                            iters,
                            relative_residual,
                            &alphas,
                            &betas,
                            CgBreakdown::NonFiniteResidual {
                                iteration: iters,
                                index,
                                value,
                            },
                        ));
                    }
                    true_r.push(value);
                }
                let true_norm = scaled_l2_norm(&true_r);
                let true_relative_residual = true_norm / b_norm;
                if !true_norm.is_finite() || !true_relative_residual.is_finite() {
                    return Ok(cg_breakdown_report(
                        iters,
                        relative_residual,
                        &alphas,
                        &betas,
                        CgBreakdown::NonFiniteResidualNorm {
                            iteration: iters,
                            value: true_relative_residual,
                        },
                    ));
                }
                relative_residual = true_relative_residual;
                let true_converged = if tol == 0.0 {
                    true_norm == 0.0
                } else {
                    relative_residual <= tol
                };
                if true_converged || iters == max_iter {
                    return Ok(CgReport {
                        termination: if true_converged {
                            CgTermination::Converged
                        } else {
                            CgTermination::IterationLimit
                        },
                        iterations: iters,
                        relative_residual,
                        condition_estimate: lanczos_condition_estimate(&alphas, &betas),
                    });
                }

                // The recursive residual met the tolerance but the independent
                // true residual did not. Restart from the verified residual;
                // the old Lanczos chain no longer describes the recurrence.
                r = true_r;
                alphas.clear();
                betas.clear();
                z = self.precondition(&r);
                p.clone_from(&z);
                rz = dot(&r, &z);
                if !rz.is_finite() {
                    return Ok(cg_breakdown_report(
                        iters,
                        relative_residual,
                        &alphas,
                        &betas,
                        CgBreakdown::NonFinitePreconditionedResidual {
                            iteration: iters,
                            value: rz,
                        },
                    ));
                }
                if rz <= 0.0 {
                    return Ok(cg_breakdown_report(
                        iters,
                        relative_residual,
                        &alphas,
                        &betas,
                        CgBreakdown::NonPositivePreconditionedResidual {
                            iteration: iters,
                            value: rz,
                        },
                    ));
                }
                continue;
            }
            z = self.precondition(&r);
            let rz_new = dot(&r, &z);
            if !rz_new.is_finite() {
                return Ok(cg_breakdown_report(
                    iters,
                    relative_residual,
                    &alphas,
                    &betas,
                    CgBreakdown::NonFinitePreconditionedResidual {
                        iteration: iters,
                        value: rz_new,
                    },
                ));
            }
            if rz_new <= 0.0 {
                return Ok(cg_breakdown_report(
                    iters,
                    relative_residual,
                    &alphas,
                    &betas,
                    CgBreakdown::NonPositivePreconditionedResidual {
                        iteration: iters,
                        value: rz_new,
                    },
                ));
            }
            let beta = rz_new / rz;
            if !beta.is_finite() || beta <= 0.0 {
                return Ok(cg_breakdown_report(
                    iters,
                    relative_residual,
                    &alphas,
                    &betas,
                    CgBreakdown::InvalidRecurrence {
                        iteration: iters,
                        value: beta,
                    },
                ));
            }
            for i in 0..n {
                let next = beta.mul_add(p[i], z[i]);
                if !next.is_finite() {
                    return Ok(cg_breakdown_report(
                        iters,
                        relative_residual,
                        &alphas,
                        &betas,
                        CgBreakdown::NonFiniteDirection {
                            iteration: iters,
                            index: i,
                            value: next,
                        },
                    ));
                }
                p[i] = next;
            }
            betas.push(beta);
            rz = rz_new;
        }
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
        // Zero islands is a degenerate request: `div_ceil(0)` divides by zero
        // and `ccds - 1` underflows. No island can be shared, so the locality
        // win is 0 — return it rather than panic (keep the method total).
        if ccds == 0 {
            return 0.0;
        }
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

fn dot(lhs: &[f64], rhs: &[f64]) -> f64 {
    lhs.iter().zip(rhs).map(|(a, b)| a * b).sum()
}

/// Scaled sum-of-squares norm (the BLAS `nrm2` recurrence).
fn scaled_l2_norm(values: &[f64]) -> f64 {
    let mut scale = 0.0f64;
    let mut sum_squares = 1.0f64;
    for value in values {
        let magnitude = value.abs();
        if !magnitude.is_finite() {
            return f64::NAN;
        }
        if magnitude == 0.0 {
            continue;
        }
        if scale < magnitude {
            let ratio = scale / magnitude;
            sum_squares = 1.0 + sum_squares * ratio * ratio;
            scale = magnitude;
        } else {
            let ratio = magnitude / scale;
            sum_squares += ratio * ratio;
        }
    }
    if scale == 0.0 {
        0.0
    } else {
        scale * det::sqrt(sum_squares)
    }
}

fn cg_breakdown_report(
    iterations: usize,
    relative_residual: f64,
    alphas: &[f64],
    betas: &[f64],
    breakdown: CgBreakdown,
) -> CgReport {
    CgReport {
        termination: CgTermination::Breakdown(breakdown),
        iterations,
        relative_residual,
        condition_estimate: lanczos_prefix_condition_estimate(alphas, betas),
    }
}

/// Retain only a completely determined tridiagonal when a later CG step
/// breaks down. Between steps, CG has already formed the beta leading into the
/// next (rejected) alpha; that trailing beta is not part of the preceding
/// `T_k` and must not erase the otherwise valid prefix.
fn lanczos_prefix_condition_estimate(
    alphas: &[f64],
    betas: &[f64],
) -> Option<LanczosConditionEstimate> {
    let complete_betas = alphas.len().checked_sub(1)?;
    match betas.len() {
        len if len == complete_betas => lanczos_condition_estimate(alphas, betas),
        len if len == alphas.len() => lanczos_condition_estimate(alphas, &betas[..complete_betas]),
        _ => None,
    }
}

/// Condition estimate from a validated CG/Lanczos tridiagonal (standard
/// coefficients-to-extreme-Ritz-values route, power/inverse-free).
fn lanczos_condition_estimate(alphas: &[f64], betas: &[f64]) -> Option<LanczosConditionEstimate> {
    let k = alphas.len();
    if k == 0 || betas.len() != k - 1 {
        return None;
    }
    if alphas
        .iter()
        .any(|alpha| !alpha.is_finite() || *alpha <= 0.0)
        || betas.iter().any(|beta| !beta.is_finite() || *beta <= 0.0)
    {
        return None;
    }
    // Tridiagonal entries (Golub–Van Loan): d_i, e_i from alpha/beta.
    let mut d = vec![0.0f64; k];
    let mut e = vec![0.0f64; k - 1];
    for i in 0..k {
        d[i] = 1.0 / alphas[i];
        if i > 0 {
            d[i] += betas[i - 1] / alphas[i - 1];
        }
        if i + 1 < k {
            e[i] = det::sqrt(betas[i]) / alphas[i];
        }
        if !d[i].is_finite() || d[i] <= 0.0 || (i < e.len() && (!e[i].is_finite() || e[i] <= 0.0)) {
            return None;
        }
    }
    let scale = d
        .iter()
        .chain(&e)
        .map(|value| value.abs())
        .fold(0.0f64, f64::max);
    if !scale.is_finite() || scale <= 0.0 {
        return None;
    }
    for value in d.iter_mut().chain(&mut e) {
        *value /= scale;
        if !value.is_finite() || *value <= 0.0 {
            return None;
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
    if !radius.is_finite() || radius <= 0.0 {
        return None;
    }
    let ritz_min = bisect_ritz(&d, &e, radius, 0)?;
    let ritz_max = bisect_ritz(&d, &e, radius, k - 1)?;
    if !ritz_min.is_finite() || !ritz_max.is_finite() || ritz_min <= 0.0 || ritz_max < ritz_min {
        return None;
    }
    let kappa = ritz_max / ritz_min;
    if !kappa.is_finite() || kappa < 1.0 {
        return None;
    }
    Some(LanczosConditionEstimate {
        kappa,
        ritz_min,
        ritz_max,
        krylov_dimension: k,
    })
}

fn sturm_count_below(d: &[f64], e: &[f64], lambda: f64) -> Option<usize> {
    let mut count = 0usize;
    let mut pivot = d[0] - lambda;
    if !pivot.is_finite() {
        return None;
    }
    if pivot < 0.0 {
        count += 1;
    }
    for i in 1..d.len() {
        let denominator = if pivot.abs() < f64::MIN_POSITIVE {
            if pivot.is_sign_negative() {
                -f64::MIN_POSITIVE
            } else {
                f64::MIN_POSITIVE
            }
        } else {
            pivot
        };
        pivot = d[i] - lambda - e[i - 1] * e[i - 1] / denominator;
        if !pivot.is_finite() {
            return None;
        }
        if pivot < 0.0 {
            count += 1;
        }
    }
    Some(count)
}

fn bisect_ritz(d: &[f64], e: &[f64], radius: f64, target: usize) -> Option<f64> {
    let (mut lo, mut hi) = (-radius, radius);
    // 1,100 halvings cover the complete binary64 exponent range; the adjacent
    // endpoint check normally stops much earlier.
    for _ in 0..1_100 {
        let mid = f64::midpoint(lo, hi);
        // Fixed-point termination of binary bisection IS bitwise identity;
        // an epsilon would stop on representable-but-unconverged midpoints.
        #[allow(clippy::float_cmp)]
        if mid == lo || mid == hi {
            break;
        }
        if sturm_count_below(d, e, mid)? > target {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let estimate = f64::midpoint(lo, hi);
    estimate.is_finite().then_some(estimate)
}

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
#[allow(clippy::float_cmp)] // Tests pin analytic kappa values (single-Ritz projections are exactly 1) bit-exactly.
#[cfg(test)]
mod tests {
    use super::{
        Bddc, Decomposition, lanczos_condition_estimate, lanczos_prefix_condition_estimate,
    };

    #[test]
    fn rho_scaling_is_a_coefficient_aware_partition_of_unity() {
        let uniform = Bddc::new(Decomposition::uniform(4, 4), false);
        for sub in &uniform.subs {
            for (&g, &weight) in sub.boundary.iter().zip(&sub.boundary_weight) {
                let multiplicity = uniform
                    .subs
                    .iter()
                    .filter(|candidate| candidate.boundary.contains(&g))
                    .count();
                let expected =
                    1.0 / f64::from(u32::try_from(multiplicity).expect("at most four subdomains"));
                assert_eq!(
                    weight.to_bits(),
                    expected.to_bits(),
                    "uniform rho-scaling must reduce exactly to counting at node {g}"
                );
            }
        }

        let jump = Bddc::new(Decomposition::checkerboard(4, 4, 1.0e6), false);
        let mut sums = vec![0.0_f64; jump.gamma.len()];
        let mut smallest = 1.0_f64;
        let mut largest = 0.0_f64;
        for sub in &jump.subs {
            for (&g, &weight) in sub.boundary.iter().zip(&sub.boundary_weight) {
                assert!(
                    weight.is_finite() && (0.0..=1.0).contains(&weight),
                    "node {g} has invalid rho weight {weight}"
                );
                sums[jump.gpos[&g]] += weight;
                smallest = smallest.min(weight);
                largest = largest.max(weight);
            }
        }
        for (slot, sum) in sums.into_iter().enumerate() {
            assert!(
                (sum - 1.0).abs() <= 8.0 * f64::EPSILON,
                "gamma slot {slot} weights sum to {sum}, not one"
            );
        }
        assert!(
            smallest < 2.0e-6 && largest > 1.0 - 2.0e-6,
            "the 1e6 jump must bias shared edge nodes: [{smallest}, {largest}]"
        );

        // Spatial falsifier for an invalid one-number-per-subdomain shortcut.
        // At interface node (2,1), the incident maxima are 8 on the left and
        // 5 on the right. The much larger coefficient at non-incident cell
        // (0,0) must not affect this node's weights.
        let mut spatial = Decomposition::uniform(2, 2);
        let n = spatial.s * spatial.m;
        spatial.rho[0] = 1.0e6;
        spatial.rho[1] = 8.0;
        spatial.rho[n + 1] = 2.0;
        spatial.rho[2] = 5.0;
        spatial.rho[n + 2] = 3.0;
        let spatial = Bddc::new(spatial, false);
        let node = spatial.decomp.node(2, 1);
        let left_slot = spatial.subs[0]
            .boundary
            .iter()
            .position(|candidate| *candidate == node)
            .expect("left trace owns target node");
        let right_slot = spatial.subs[1]
            .boundary
            .iter()
            .position(|candidate| *candidate == node)
            .expect("right trace owns target node");
        let left = spatial.subs[0].boundary_weight[left_slot];
        let right = spatial.subs[1].boundary_weight[right_slot];
        assert!((left - 8.0 / 13.0).abs() <= f64::EPSILON, "{left}");
        assert!((right - 5.0 / 13.0).abs() <= f64::EPSILON, "{right}");
    }

    #[test]
    fn lanczos_condition_estimate_rejects_malformed_evidence() {
        assert_eq!(lanczos_condition_estimate(&[], &[]), None);
        assert_eq!(lanczos_condition_estimate(&[1.0, 1.0], &[]), None);
        assert_eq!(lanczos_condition_estimate(&[0.0], &[]), None);
        assert_eq!(lanczos_condition_estimate(&[f64::NAN], &[]), None);
        assert_eq!(lanczos_condition_estimate(&[f64::from_bits(1)], &[]), None);
        assert_eq!(lanczos_condition_estimate(&[1.0, 1.0], &[-0.25]), None);
        assert_eq!(
            lanczos_condition_estimate(&[1.0, 1.0], &[f64::INFINITY]),
            None
        );
    }

    #[test]
    fn lanczos_condition_estimate_recovers_known_ritz_ratio() {
        let one = lanczos_condition_estimate(&[0.5], &[]).expect("one valid Ritz value");
        assert_eq!(one.kappa, 1.0);
        assert_eq!(one.krylov_dimension, 1);

        // alpha=[1/2, 2/3], beta=[1/4] reconstructs [[2, 1], [1, 2]],
        // whose eigenvalues are 1 and 3 (the internal normalization preserves
        // their ratio).
        let estimate = lanczos_condition_estimate(&[0.5, 2.0 / 3.0], &[0.25])
            .expect("finite positive-definite two-step tridiagonal");
        assert_eq!(estimate.krylov_dimension, 2);
        assert!((estimate.kappa - 3.0).abs() < 1e-12, "{estimate:?}");
        assert!(estimate.ritz_min > 0.0);
        assert!(estimate.ritz_max > estimate.ritz_min);
    }

    #[test]
    fn breakdown_estimate_ignores_only_one_pending_beta() {
        let prefix = lanczos_prefix_condition_estimate(&[0.5], &[0.25])
            .expect("one accepted alpha remains a valid one-dimensional prefix");
        assert_eq!(prefix.kappa, 1.0);
        assert_eq!(prefix.krylov_dimension, 1);
        assert_eq!(
            lanczos_prefix_condition_estimate(&[0.5], &[0.25, 0.5]),
            None,
            "more than one unmatched beta is malformed evidence"
        );
    }
}

/// Crate version, re-exported for provenance stamping.
