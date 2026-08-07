//! Sparse symmetric direct factorization: supernodal multifrontal LDLᵀ with
//! AMD fill-reducing ordering, restricted (within-supernode) 1×1/2×2
//! threshold pivoting, and Sylvester inertia (bead
//! frankensim-fsim-sparse-direct-4a38j, musical-acoustics program).
//!
//! WHY THIS EXISTS — shift-invert modal analysis factors (K − σM), which is
//! symmetric INDEFINITE for interior shifts, and the inertia of D is exactly
//! the eigenvalue count below σ that spectrum-slicing certification needs
//! (Sylvester's law). Cholesky alone cannot serve that consumer; iterative
//! solvers cannot produce inertia at all.
//!
//! DESIGN — three explicit stages with separate cost models:
//! 1. [`amd_order`]: approximate-minimum-degree ordering (Amestoy–Davis–Duff
//!    quotient-graph algorithm; v1 without supervariable merging — an
//!    ordering-QUALITY boundary, never a correctness one).
//! 2. [`SymbolicLdlt::analyze`]: elimination tree, postorder, per-column L
//!    structure, and fundamental supernodes (exact structure-match test, so
//!    within-supernode column exchanges provably preserve the symbolic
//!    factorization). Reusable across every matrix with the SAME pattern —
//!    each shift σ in (K − σM) refactors numerically but never re-analyzes.
//! 3. [`SymbolicLdlt::factor`]: multifrontal numeric factorization. Each
//!    supernode assembles a dense frontal matrix (original entries plus
//!    children's Schur updates), factors its fully-summed block with
//!    threshold pivoting restricted to that block (Duff–Reid growth test,
//!    1×1 and symmetric 2×2 pivots), and passes the Schur complement to its
//!    assembly-tree parent.
//!
//! Pivot breakdown REFUSES with a named error instead of perturbing: a
//! perturbed factorization would silently corrupt the inertia count that
//! downstream certification treats as authority (refusal-not-clamp,
//! workspace law). Delayed pivots (migrating a failed column to the parent
//! front) are the recorded follow-up; restricted pivoting is not a backward
//! stability guarantee for adversarial indefinite inputs.
//!
//! Determinism: every stage is sequential with index-ordered tie-breaking —
//! repeat factorization of the same matrix is bitwise identical (tested).

use crate::Csr;

/// Fill-reducing ordering strategy for [`SymbolicLdlt::analyze`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DirectOrdering {
    /// Identity permutation — the fill-comparison baseline and the debug aid.
    Natural,
    /// Approximate minimum degree ([`amd_order`]). The default.
    Amd,
}

/// Options for the numeric factorization stage.
#[derive(Debug, Clone, Copy)]
pub struct LdltOptions {
    /// Relative pivot threshold `u` in (0, 0.5] (Duff–Reid). A 1×1 pivot
    /// `a_kk` is accepted only when `|a_kk| ≥ u·γ_k` (γ_k = largest
    /// off-diagonal magnitude in the candidate column within the front); a
    /// 2×2 pivot must bound the multiplier growth by `1/u`. Larger `u` is
    /// more stable, smaller `u` accepts more pivots. Default `0.01`.
    pub pivot_threshold: f64,
    /// Permit symmetric 2×2 pivots (required for saddle-point/zero-diagonal
    /// indefinite matrices). Default `true`; disabling exists so tests can
    /// prove the 2×2 path is load-bearing.
    pub allow_2x2: bool,
    /// Relative zero-pivot floor: a pivot (or 2×2 block determinant, squared
    /// scale) below `zero_pivot_rel · max|A|` (resp. its square) refuses as
    /// numerically singular. Default `256·ε ≈ 5.7e-14`.
    pub zero_pivot_rel: f64,
}

impl Default for LdltOptions {
    fn default() -> LdltOptions {
        LdltOptions {
            pivot_threshold: 0.01,
            allow_2x2: true,
            zero_pivot_rel: 256.0 * f64::EPSILON,
        }
    }
}

/// Sylvester inertia of the factored matrix: eigenvalue sign counts of D
/// (equal to those of A by Sylvester's law of inertia). `zero` is always 0
/// here — a numerically singular pivot refuses instead of recording a zero
/// eigenvalue claim.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Inertia {
    /// Number of positive eigenvalues.
    pub positive: usize,
    /// Number of negative eigenvalues.
    pub negative: usize,
}

/// Work/memory/pivot accounting for one numeric factorization.
#[derive(Debug, Clone, Copy)]
pub struct FactorStats {
    /// Matrix dimension.
    pub n: usize,
    /// Stored nonzeros of the input matrix (both triangles as given).
    pub nnz_a: usize,
    /// Nonzeros of L including the unit diagonal.
    pub nnz_l: usize,
    /// `nnz_l / nnz_a` (1.0 when `nnz_a` is 0).
    pub fill_ratio: f64,
    /// Approximate floating-point operation count of the numeric phase.
    pub flops: u64,
    /// Peak bytes of live frontal + update-matrix storage.
    pub peak_front_bytes: usize,
    /// Number of supernodes.
    pub supernodes: usize,
    /// Count of 1×1 pivots taken.
    pub pivots_1x1: usize,
    /// Count of 2×2 pivots taken.
    pub pivots_2x2: usize,
    /// Largest frontal matrix dimension encountered.
    pub max_front_dim: usize,
}

/// Typed refusals of the direct-factorization stages. Display strings carry
/// stable `FS-SPARSE-DIRECT-*` codes for log matching.
#[derive(Debug, Clone, PartialEq)]
pub enum LdltError {
    /// The matrix is not square.
    NotSquare {
        /// Rows of the offending matrix.
        nrows: usize,
        /// Columns of the offending matrix.
        ncols: usize,
    },
    /// The sparsity pattern is not structurally symmetric; `(row, col)` is a
    /// witness entry whose mirror is absent.
    StructurallyAsymmetric {
        /// Witness row.
        row: usize,
        /// Witness column.
        col: usize,
    },
    /// A stored value is NaN or infinite. Non-finite input would contaminate
    /// every downstream claim, so it fails closed here.
    NonFiniteEntry {
        /// Row of the non-finite entry.
        row: usize,
        /// Column of the non-finite entry.
        col: usize,
    },
    /// `factor` was called with a matrix whose pattern differs from the one
    /// this symbolic analysis was computed for.
    PatternMismatch,
    /// Invalid [`LdltOptions`] field (threshold outside (0, 0.5], or a
    /// non-positive zero-pivot floor).
    InvalidOptions {
        /// Name of the offending option field.
        field: &'static str,
    },
    /// No acceptable 1×1 or 2×2 pivot exists inside a supernode's
    /// fully-summed block: the matrix is numerically singular there, or
    /// restricted pivoting is insufficient for it (delayed pivots are the
    /// recorded follow-up). The factorization REFUSES rather than perturbs.
    PivotBreakdown {
        /// Column in the ORIGINAL index space where elimination stopped.
        orig_col: usize,
        /// The same column in elimination (permuted) order.
        permuted_col: usize,
        /// Best candidate pivot magnitude seen in the block.
        best_abs: f64,
        /// The magnitude the threshold test demanded.
        required: f64,
        /// Whether 2×2 pivots were permitted during the search.
        two_by_two_enabled: bool,
    },
}

impl core::fmt::Display for LdltError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            LdltError::NotSquare { nrows, ncols } => {
                write!(f, "FS-SPARSE-DIRECT-NOT-SQUARE: {nrows}x{ncols} matrix")
            }
            LdltError::StructurallyAsymmetric { row, col } => write!(
                f,
                "FS-SPARSE-DIRECT-ASYMMETRIC-PATTERN: entry ({row},{col}) has no stored mirror"
            ),
            LdltError::NonFiniteEntry { row, col } => write!(
                f,
                "FS-SPARSE-DIRECT-NON-FINITE: entry ({row},{col}) is NaN or infinite"
            ),
            LdltError::PatternMismatch => write!(
                f,
                "FS-SPARSE-DIRECT-PATTERN-MISMATCH: matrix pattern differs from the analyzed one"
            ),
            LdltError::InvalidOptions { field } => {
                write!(f, "FS-SPARSE-DIRECT-INVALID-OPTIONS: {field}")
            }
            LdltError::PivotBreakdown {
                orig_col,
                permuted_col,
                best_abs,
                required,
                two_by_two_enabled,
            } => write!(
                f,
                "FS-SPARSE-DIRECT-PIVOT-BREAKDOWN: no acceptable pivot at original column \
                 {orig_col} (elimination position {permuted_col}); best |candidate| {best_abs:.3e} \
                 < required {required:.3e}; 2x2 pivots {}",
                if *two_by_two_enabled {
                    "enabled"
                } else {
                    "disabled"
                }
            ),
        }
    }
}

impl std::error::Error for LdltError {}

// ---------------------------------------------------------------------------
// AMD ordering
// ---------------------------------------------------------------------------

/// Approximate minimum degree ordering (Amestoy–Davis–Duff) on the symmetric
/// pattern of `a` (the pattern is symmetrized internally; values are
/// ignored). Returns `perm` with `perm[k]` = the original index eliminated
/// k-th. Deterministic: minimum approximate degree with smallest-index
/// tie-breaking.
///
/// v1 boundary: no supervariable (indistinguishable-node) merging and no
/// dense-row special case — both affect ordering quality and speed on some
/// graphs, never correctness. Panics on a non-square matrix (programmer
/// error; [`SymbolicLdlt::analyze`] refuses with a typed error first).
#[must_use]
#[allow(clippy::needless_range_loop)] // index-parallel arrays throughout
pub fn amd_order(a: &Csr) -> Vec<usize> {
    assert_eq!(a.nrows(), a.ncols(), "amd_order: matrix must be square");
    let n = a.nrows();
    if n == 0 {
        return Vec::new();
    }
    // Symmetrized adjacency without the diagonal.
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for r in 0..n {
        let (cols, _) = a.row(r);
        for &c in cols {
            if c != r {
                adj[r].push(c);
            }
        }
    }
    for r in 0..n {
        let (cols, _) = a.row(r);
        for &c in cols {
            if c != r && !a_contains(a, c, r) {
                adj[c].push(r);
            }
        }
    }
    for l in &mut adj {
        l.sort_unstable();
        l.dedup();
    }

    // Quotient graph state. Element ids reuse the pivot variable's id.
    let mut elems_of: Vec<Vec<usize>> = vec![Vec::new(); n]; // elements adjacent to variable
    let mut elem_vars: Vec<Vec<usize>> = vec![Vec::new(); n]; // element -> boundary variables
    let mut elem_alive = vec![false; n];
    let mut var_alive = vec![true; n];
    let mut degree: Vec<usize> = adj.iter().map(Vec::len).collect();

    // Lazy min-heap of (degree, index).
    let mut heap: std::collections::BinaryHeap<core::cmp::Reverse<(usize, usize)>> =
        (0..n).map(|i| core::cmp::Reverse((degree[i], i))).collect();

    let mut perm = Vec::with_capacity(n);
    let mut marked = vec![usize::MAX; n]; // stamp array for Lp membership
    let mut w_stamp = vec![usize::MAX; n]; // stamp array for |Le \ Lp| counters
    let mut w_val = vec![0usize; n];

    for step in 0..n {
        // Pop until a live entry whose recorded degree is current.
        let p = loop {
            let core::cmp::Reverse((d, i)) = heap.pop().expect("amd: heap exhausted early");
            if var_alive[i] && degree[i] == d {
                break i;
            }
        };
        var_alive[p] = false;
        perm.push(p);

        // Lp = (adj[p] ∪ ⋃_{e ∈ E_p} Le) \ {p}, live variables only.
        let mut lp: Vec<usize> = adj[p].iter().copied().filter(|&v| var_alive[v]).collect();
        for &e in &elems_of[p] {
            if elem_alive[e] {
                lp.extend(
                    elem_vars[e]
                        .iter()
                        .copied()
                        .filter(|&v| var_alive[v] && v != p),
                );
            }
        }
        lp.sort_unstable();
        lp.dedup();
        for &e in &elems_of[p] {
            elem_alive[e] = false; // absorbed into the new element p
        }
        for &v in &lp {
            marked[v] = step;
        }

        // |Le \ Lp| for every live element touching Lp, via the classic
        // decrement pass (i ∈ E_e ∩ Lp ⇔ e ∈ E_i for i ∈ Lp).
        for &i in &lp {
            for &e in &elems_of[i] {
                if elem_alive[e] {
                    if w_stamp[e] != step {
                        w_stamp[e] = step;
                        w_val[e] = elem_vars[e].iter().filter(|&&v| var_alive[v]).count();
                    }
                    w_val[e] -= 1;
                }
            }
        }
        // Absorb elements entirely contained in Lp.
        for &i in &lp {
            for &e in &elems_of[i] {
                if elem_alive[e] && w_stamp[e] == step && w_val[e] == 0 {
                    elem_alive[e] = false;
                }
            }
        }

        let remaining = n - step - 1;
        for &i in &lp {
            // Prune variable adjacency: members of Lp now reach i through
            // element p; dead variables drop.
            adj[i].retain(|&v| var_alive[v] && marked[v] != step);
            elems_of[i].retain(|&e| elem_alive[e]);
            elems_of[i].push(p);
            // Approximate external degree (ADD bound).
            let mut d = adj[i].len() + (lp.len() - 1);
            for &e in &elems_of[i] {
                if e != p && w_stamp[e] == step {
                    d += w_val[e];
                }
            }
            let d = d.min(remaining);
            degree[i] = d;
            heap.push(core::cmp::Reverse((d, i)));
        }
        elem_vars[p] = lp;
        elem_alive[p] = true;
        elems_of[p].clear();
    }
    perm
}

fn a_contains(a: &Csr, r: usize, c: usize) -> bool {
    let (cols, _) = a.row(r);
    cols.binary_search(&c).is_ok()
}

// ---------------------------------------------------------------------------
// Symbolic analysis
// ---------------------------------------------------------------------------

const NONE: usize = usize::MAX;

#[derive(Debug, Clone)]
struct SnodeSym {
    /// First column (elimination order) of this supernode.
    col_start: usize,
    /// Number of columns.
    width: usize,
    /// Below-block row structure (elimination-order ids, ascending) — shared
    /// by every column of the supernode (exact-match certified).
    panel: Vec<usize>,
    /// Assembly-tree parent supernode, or `NONE`.
    parent: usize,
    /// Children supernodes (ascending).
    children: Vec<usize>,
}

/// Reusable symbolic factorization of a symmetric sparsity pattern:
/// permutation, elimination tree, supernode partition, and L structure.
/// One `analyze` serves every numeric [`SymbolicLdlt::factor`] with the
/// same pattern — the shift ladder (K − σᵢM) analyzes once and factors per
/// shift.
#[derive(Debug, Clone)]
pub struct SymbolicLdlt {
    n: usize,
    /// `perm[k]` = original index eliminated k-th (AMD ∘ etree postorder).
    perm: Vec<usize>,
    /// Inverse of `perm`.
    iperm: Vec<usize>,
    snodes: Vec<SnodeSym>,
    /// nnz(L) including the unit diagonal.
    nnz_l: usize,
    /// Exact pattern the analysis was computed for (mismatch guard).
    pat_row_ptr: Vec<usize>,
    pat_col_idx: Vec<usize>,
}

impl SymbolicLdlt {
    /// Analyze the symmetric pattern of `a`: validate squareness and
    /// structural symmetry, order (AMD or natural), build the elimination
    /// tree and its postorder, compute per-column L structures, and certify
    /// fundamental supernodes by EXACT structure comparison.
    pub fn analyze(a: &Csr, ordering: DirectOrdering) -> Result<SymbolicLdlt, LdltError> {
        if a.nrows() != a.ncols() {
            return Err(LdltError::NotSquare {
                nrows: a.nrows(),
                ncols: a.ncols(),
            });
        }
        let n = a.nrows();
        // Structural symmetry: every (r,c) needs its mirror.
        for r in 0..n {
            let (cols, _) = a.row(r);
            for &c in cols {
                if !a_contains(a, c, r) {
                    return Err(LdltError::StructurallyAsymmetric { row: r, col: c });
                }
            }
        }
        let base_perm = match ordering {
            DirectOrdering::Natural => (0..n).collect::<Vec<usize>>(),
            DirectOrdering::Amd => amd_order(a),
        };
        let mut base_iperm = vec![0usize; n];
        for (k, &o) in base_perm.iter().enumerate() {
            base_iperm[o] = k;
        }
        // Lower-triangular adjacency in base elimination order.
        let ladj = build_lower_adj(a, &base_iperm);
        let parent0 = etree(&ladj);
        let post = postorder(&parent0);
        // Compose: final elimination order.
        let perm: Vec<usize> = post.iter().map(|&k| base_perm[k]).collect();
        let mut iperm = vec![0usize; n];
        for (k, &o) in perm.iter().enumerate() {
            iperm[o] = k;
        }
        let ladj = build_lower_adj(a, &iperm);
        let parent = etree(&ladj); // relabeled tree, recomputed for simplicity
        let colstruct = col_structs(&ladj, &parent);

        // Fundamental supernodes, certified by exact structure match:
        // col j may join col j+1's supernode iff struct(j) = {j+1} ∪ struct(j+1).
        let mut snodes: Vec<SnodeSym> = Vec::new();
        let mut snode_of = vec![NONE; n.max(1)];
        let mut j = 0usize;
        while j < n {
            let start = j;
            while j + 1 < n
                && parent[j] == j + 1
                && colstruct[j].len() == colstruct[j + 1].len() + 1
                && colstruct[j].first() == Some(&(j + 1))
                && colstruct[j][1..] == colstruct[j + 1][..]
            {
                j += 1;
            }
            let width = j - start + 1;
            let panel: Vec<usize> = colstruct[start][width - 1..].to_vec();
            let id = snodes.len();
            for c in start..=j {
                snode_of[c] = id;
            }
            snodes.push(SnodeSym {
                col_start: start,
                width,
                panel,
                parent: NONE,
                children: Vec::new(),
            });
            j += 1;
        }
        for id in 0..snodes.len() {
            let last = snodes[id].col_start + snodes[id].width - 1;
            let p = parent[last];
            if p != NONE {
                let ps = snode_of[p];
                snodes[id].parent = ps;
                snodes[ps].children.push(id);
            }
        }
        let nnz_l = n + colstruct.iter().map(Vec::len).sum::<usize>();
        let (pat_row_ptr, pat_col_idx) = pattern_of(a);
        Ok(SymbolicLdlt {
            n,
            perm,
            iperm,
            snodes,
            nnz_l,
            pat_row_ptr,
            pat_col_idx,
        })
    }

    /// Dimension of the analyzed pattern.
    #[must_use]
    pub fn n(&self) -> usize {
        self.n
    }

    /// nnz(L) including the unit diagonal (a pure function of the pattern
    /// and ordering — available before any numeric work, so callers can
    /// budget memory).
    #[must_use]
    pub fn nnz_l(&self) -> usize {
        self.nnz_l
    }

    /// Number of supernodes.
    #[must_use]
    pub fn supernode_count(&self) -> usize {
        self.snodes.len()
    }

    /// The elimination permutation (`perm[k]` = original index eliminated
    /// k-th).
    #[must_use]
    pub fn permutation(&self) -> &[usize] {
        &self.perm
    }

    /// Numeric multifrontal factorization of `a` (which must carry exactly
    /// the analyzed pattern) as P·A·Pᵀ = L·D·Lᵀ with unit-lower L and D of
    /// 1×1/2×2 blocks. Refuses (never perturbs) on pivot breakdown.
    #[allow(clippy::too_many_lines)] // one coherent numeric pipeline; splitting obscures the data flow
    pub fn factor(&self, a: &Csr, opts: &LdltOptions) -> Result<LdltFactor, LdltError> {
        if !(opts.pivot_threshold > 0.0 && opts.pivot_threshold <= 0.5) {
            return Err(LdltError::InvalidOptions {
                field: "pivot_threshold must be in (0, 0.5]",
            });
        }
        if !(opts.zero_pivot_rel > 0.0) {
            return Err(LdltError::InvalidOptions {
                field: "zero_pivot_rel must be positive",
            });
        }
        let (rp, ci) = pattern_of(a);
        if rp != self.pat_row_ptr || ci != self.pat_col_idx {
            return Err(LdltError::PatternMismatch);
        }
        let n = self.n;
        let mut max_abs = 0.0f64;
        for r in 0..n {
            let (cols, vals) = a.row(r);
            for (&c, &v) in cols.iter().zip(vals) {
                if !v.is_finite() {
                    return Err(LdltError::NonFiniteEntry { row: r, col: c });
                }
                max_abs = max_abs.max(v.abs());
            }
        }
        let zero_tol = opts.zero_pivot_rel * max_abs.max(f64::MIN_POSITIVE);

        let mut stats = FactorStats {
            n,
            nnz_a: a.nnz(),
            nnz_l: self.nnz_l,
            fill_ratio: if a.nnz() == 0 {
                1.0
            } else {
                self.nnz_l as f64 / a.nnz() as f64
            },
            flops: 0,
            peak_front_bytes: 0,
            supernodes: self.snodes.len(),
            pivots_1x1: 0,
            pivots_2x2: 0,
            max_front_dim: 0,
        };
        let mut inertia = Inertia {
            positive: 0,
            negative: 0,
        };
        let mut snodes_num: Vec<SnodeNum> = Vec::with_capacity(self.snodes.len());
        let mut updates: Vec<Option<Update>> = vec![None; self.snodes.len()];
        let mut live_update_bytes = 0usize;

        for (sid, sym) in self.snodes.iter().enumerate() {
            let w = sym.width;
            let p = sym.panel.len();
            let m = w + p;
            stats.max_front_dim = stats.max_front_dim.max(m);
            // Front variables (elimination-order ids), ascending: block ∪ panel.
            let mut front: Vec<usize> = Vec::with_capacity(m);
            front.extend(sym.col_start..sym.col_start + w);
            front.extend(sym.panel.iter().copied());
            let mut f = vec![0.0f64; m * m];
            let front_bytes = m * m * 8;
            stats.peak_front_bytes = stats.peak_front_bytes.max(front_bytes + live_update_bytes);

            // Scatter original entries of the block columns (values read
            // from the row of the eliminated column; the input contract is a
            // numerically symmetric matrix).
            for c in 0..w {
                let j = sym.col_start + c;
                let orig = self.perm[j];
                let (cols, vals) = a.row(orig);
                for (&oc, &v) in cols.iter().zip(vals) {
                    let pi = self.iperm[oc];
                    if pi >= j {
                        let li = front
                            .binary_search(&pi)
                            .expect("symbolic invariant: neighbor outside front");
                        f[li + c * m] += v;
                        if li != c {
                            f[c + li * m] += v;
                        }
                    }
                }
            }
            // Extend-add the children's Schur updates.
            for &ch in &sym.children {
                let up = updates[ch].take().expect("child update consumed twice");
                live_update_bytes -= up.data.len() * 8;
                let map: Vec<usize> = up
                    .vars
                    .iter()
                    .map(|v| {
                        front
                            .binary_search(v)
                            .expect("assembly invariant: child var outside parent front")
                    })
                    .collect();
                let q = up.vars.len();
                for cj in 0..q {
                    let dst_c = map[cj];
                    for ri in 0..q {
                        f[map[ri] + dst_c * m] += up.data[ri + cj * q];
                    }
                }
            }

            // Factor the fully-summed w×w block with restricted pivoting.
            let mut local2var = front.clone();
            let mut dblocks: Vec<DBlock> = Vec::with_capacity(w);
            let mut k = 0usize;
            while k < w {
                let choice = choose_pivot(&f, m, w, k, opts, zero_tol);
                match choice {
                    PivotChoice::One(q) => {
                        if q != k {
                            swap_sym(&mut f, m, k, q);
                            local2var.swap(k, q);
                        }
                        let d = f[k + k * m];
                        if d > 0.0 {
                            inertia.positive += 1;
                        } else {
                            inertia.negative += 1;
                        }
                        stats.pivots_1x1 += 1;
                        dblocks.push(DBlock::Single(d));
                        for i in k + 1..m {
                            f[i + k * m] /= d;
                        }
                        for jc in k + 1..m {
                            let cj = d * f[jc + k * m];
                            for i in jc..m {
                                let t = f[i + k * m] * cj;
                                f[i + jc * m] -= t;
                                if i != jc {
                                    f[jc + i * m] = f[i + jc * m];
                                }
                            }
                        }
                        stats.flops += ((m - k) as u64) * ((m - k) as u64);
                        k += 1;
                    }
                    PivotChoice::Two(pq, qq) => {
                        if pq != k {
                            swap_sym(&mut f, m, k, pq);
                            local2var.swap(k, pq);
                        }
                        // After the first swap, the partner may have moved.
                        let qq = if qq == k { pq } else { qq };
                        if qq != k + 1 {
                            swap_sym(&mut f, m, k + 1, qq);
                            local2var.swap(k + 1, qq);
                        }
                        let e11 = f[k + k * m];
                        let e21 = f[(k + 1) + k * m];
                        let e22 = f[(k + 1) + (k + 1) * m];
                        let det = e11 * e22 - e21 * e21;
                        if det < 0.0 {
                            inertia.positive += 1;
                            inertia.negative += 1;
                        } else if e11 + e22 > 0.0 {
                            inertia.positive += 2;
                        } else {
                            inertia.negative += 2;
                        }
                        stats.pivots_2x2 += 1;
                        dblocks.push(DBlock::Pair(e11, e21, e22));
                        let inv11 = e22 / det;
                        let inv21 = -e21 / det;
                        let inv22 = e11 / det;
                        let acol: Vec<f64> = (k + 2..m).map(|i| f[i + k * m]).collect();
                        let bcol: Vec<f64> = (k + 2..m).map(|i| f[i + (k + 1) * m]).collect();
                        for (idx, i) in (k + 2..m).enumerate() {
                            let l1 = acol[idx].mul_add(inv11, bcol[idx] * inv21);
                            let l2 = acol[idx].mul_add(inv21, bcol[idx] * inv22);
                            f[i + k * m] = l1;
                            f[i + (k + 1) * m] = l2;
                        }
                        for (jdx, jc) in (k + 2..m).enumerate() {
                            let aj = acol[jdx];
                            let bj = bcol[jdx];
                            for i in jc..m {
                                let t = f[i + k * m].mul_add(aj, f[i + (k + 1) * m] * bj);
                                f[i + jc * m] -= t;
                                if i != jc {
                                    f[jc + i * m] = f[i + jc * m];
                                }
                            }
                        }
                        stats.flops += 2 * ((m - k) as u64) * ((m - k) as u64);
                        k += 2;
                    }
                    PivotChoice::Fail { best_abs, required } => {
                        return Err(LdltError::PivotBreakdown {
                            orig_col: self.perm[local2var[k]],
                            permuted_col: local2var[k],
                            best_abs,
                            required,
                            two_by_two_enabled: opts.allow_2x2,
                        });
                    }
                }
            }

            // Harvest L (unit-diagonal dense trapezoid, col-major m×w).
            let mut l = vec![0.0f64; m * w];
            for c in 0..w {
                l[c + c * m] = 1.0;
                for i in c + 1..m {
                    l[i + c * m] = f[i + c * m];
                }
            }
            // Schur complement → parent update. A root supernode provably
            // has an empty panel (struct(last col) empty ⟺ no etree parent).
            debug_assert!(sym.parent != NONE || p == 0, "root supernode with panel");
            if sym.parent != NONE {
                let mut data = vec![0.0f64; p * p];
                for cj in 0..p {
                    for ri in 0..p {
                        data[ri + cj * p] = f[(w + ri) + (w + cj) * m];
                    }
                }
                live_update_bytes += data.len() * 8;
                stats.peak_front_bytes = stats.peak_front_bytes.max(live_update_bytes);
                updates[sid] = Some(Update {
                    vars: sym.panel.clone(),
                    data,
                });
            }
            snodes_num.push(SnodeNum {
                block_vars: local2var[..w].to_vec(),
                panel_vars: local2var[w..].to_vec(),
                l,
                d: dblocks,
            });
        }

        Ok(LdltFactor {
            n,
            perm: self.perm.clone(),
            snodes: snodes_num,
            inertia,
            stats,
        })
    }
}

fn pattern_of(a: &Csr) -> (Vec<usize>, Vec<usize>) {
    let n = a.nrows();
    let mut row_ptr = Vec::with_capacity(n + 1);
    let mut col_idx = Vec::with_capacity(a.nnz());
    row_ptr.push(0);
    for r in 0..n {
        let (cols, _) = a.row(r);
        col_idx.extend_from_slice(cols);
        row_ptr.push(col_idx.len());
    }
    (row_ptr, col_idx)
}

/// Lower-triangular adjacency in elimination order: for permuted row `i`,
/// the sorted permuted neighbors `j < i`.
fn build_lower_adj(a: &Csr, iperm: &[usize]) -> Vec<Vec<usize>> {
    let n = a.nrows();
    let mut ladj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for r in 0..n {
        let (cols, _) = a.row(r);
        let pi = iperm[r];
        for &c in cols {
            let pj = iperm[c];
            if pj < pi {
                ladj[pi].push(pj);
            }
        }
    }
    for l in &mut ladj {
        l.sort_unstable();
    }
    ladj
}

/// Elimination tree (Liu's algorithm with path compression).
fn etree(ladj: &[Vec<usize>]) -> Vec<usize> {
    let n = ladj.len();
    let mut parent = vec![NONE; n];
    let mut ancestor = vec![NONE; n];
    for i in 0..n {
        for &j in &ladj[i] {
            let mut t = j;
            while ancestor[t] != NONE && ancestor[t] != i {
                let next = ancestor[t];
                ancestor[t] = i;
                t = next;
            }
            if ancestor[t] == NONE {
                ancestor[t] = i;
                parent[t] = i;
            }
        }
    }
    parent
}

/// Iterative etree postorder with ascending-child determinism.
fn postorder(parent: &[usize]) -> Vec<usize> {
    let n = parent.len();
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut roots = Vec::new();
    for j in 0..n {
        if parent[j] == NONE {
            roots.push(j);
        } else {
            children[parent[j]].push(j);
        }
    }
    let mut post = Vec::with_capacity(n);
    let mut stack: Vec<(usize, usize)> = Vec::new();
    for &r in &roots {
        stack.push((r, 0));
        while let Some(&(node, next)) = stack.last() {
            if next < children[node].len() {
                stack.last_mut().expect("nonempty stack").1 += 1;
                stack.push((children[node][next], 0));
            } else {
                post.push(node);
                stack.pop();
            }
        }
    }
    post
}

/// Per-column below-diagonal L structure via the row-wise etree traversal
/// (each row i is added to the structure of every column on the path from
/// its lower neighbors up to i). Rows arrive in ascending order, so each
/// structure vector is sorted.
fn col_structs(ladj: &[Vec<usize>], parent: &[usize]) -> Vec<Vec<usize>> {
    let n = ladj.len();
    let mut cs: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut stamp = vec![NONE; n];
    for i in 0..n {
        stamp[i] = i;
        for &j in &ladj[i] {
            let mut t = j;
            while t != i && stamp[t] != i {
                cs[t].push(i);
                stamp[t] = i;
                let p = parent[t];
                if p == NONE {
                    break;
                }
                t = p;
            }
        }
    }
    cs
}

enum PivotChoice {
    One(usize),
    Two(usize, usize),
    Fail { best_abs: f64, required: f64 },
}

/// Largest off-diagonal magnitude of local column `c` among rows `k..m`.
fn gamma_of(f: &[f64], m: usize, k: usize, c: usize) -> f64 {
    let mut g = 0.0f64;
    for i in k..m {
        if i != c {
            g = g.max(f[i + c * m].abs());
        }
    }
    g
}

/// Deterministic restricted pivot search inside the fully-summed block
/// `k..w` of an m×m front: first acceptable 1×1 (current column, then
/// ascending scan), else the lexicographically first acceptable 2×2 under
/// the Duff–Reid growth bound, else a typed failure report.
fn choose_pivot(
    f: &[f64],
    m: usize,
    w: usize,
    k: usize,
    opts: &LdltOptions,
    zero_tol: f64,
) -> PivotChoice {
    let u = opts.pivot_threshold;
    let mut best_abs = 0.0f64;
    let mut required = zero_tol;
    for q in k..w {
        let dq = f[q + q * m].abs();
        let gq = gamma_of(f, m, k, q);
        best_abs = best_abs.max(dq);
        required = required.max(u * gq).max(zero_tol);
        if dq > zero_tol && dq >= u * gq {
            return PivotChoice::One(q);
        }
    }
    if opts.allow_2x2 && w - k >= 2 {
        let det_tol = zero_tol * zero_tol;
        for p in k..w {
            for q in p + 1..w {
                let e11 = f[p + p * m];
                let e21 = f[q + p * m];
                let e22 = f[q + q * m];
                let det = e11 * e22 - e21 * e21;
                if det.abs() <= det_tol {
                    continue;
                }
                let gp = gamma_pair(f, m, k, p, q);
                let gq = gamma_pair(f, m, k, q, p);
                let g1 = e22.abs().mul_add(gp, e21.abs() * gq) / det.abs();
                let g2 = e21.abs().mul_add(gp, e11.abs() * gq) / det.abs();
                if g1.max(g2) <= 1.0 / u {
                    return PivotChoice::Two(p, q);
                }
            }
        }
    }
    PivotChoice::Fail { best_abs, required }
}

/// Largest magnitude in local column `c` among rows `k..m` excluding both
/// pivot-pair rows.
fn gamma_pair(f: &[f64], m: usize, k: usize, c: usize, other: usize) -> f64 {
    let mut g = 0.0f64;
    for i in k..m {
        if i != c && i != other {
            g = g.max(f[i + c * m].abs());
        }
    }
    g
}

/// Symmetric row+column swap of a full col-major m×m matrix.
fn swap_sym(f: &mut [f64], m: usize, a: usize, b: usize) {
    if a == b {
        return;
    }
    for r in 0..m {
        f.swap(r + a * m, r + b * m);
    }
    for c in 0..m {
        f.swap(a + c * m, b + c * m);
    }
}

#[derive(Debug, Clone)]
struct Update {
    vars: Vec<usize>,
    data: Vec<f64>,
}

#[derive(Debug, Clone, Copy)]
enum DBlock {
    Single(f64),
    Pair(f64, f64, f64), // (d11, d21, d22)
}

#[derive(Debug, Clone)]
struct SnodeNum {
    /// Block columns in PIVOT order (elimination-order ids).
    block_vars: Vec<usize>,
    /// Panel rows (elimination-order ids, ascending).
    panel_vars: Vec<usize>,
    /// Dense unit-diagonal trapezoid, col-major (block+panel) × width.
    l: Vec<f64>,
    /// D blocks in pivot order.
    d: Vec<DBlock>,
}

/// Numeric LDLᵀ factorization: P·A·Pᵀ = L·D·Lᵀ. Produced by
/// [`SymbolicLdlt::factor`]; serves any number of right-hand sides via
/// [`LdltFactor::solve`] and exposes the Sylvester [`Inertia`] that
/// spectrum-slicing certification consumes.
#[derive(Debug, Clone)]
pub struct LdltFactor {
    n: usize,
    perm: Vec<usize>,
    snodes: Vec<SnodeNum>,
    inertia: Inertia,
    stats: FactorStats,
}

impl LdltFactor {
    /// Sylvester inertia of A (sign counts of D). `zero` is structurally
    /// absent: singular pivots refuse during factorization.
    #[must_use]
    pub fn inertia(&self) -> Inertia {
        self.inertia
    }

    /// Work/memory/pivot accounting of the factorization.
    #[must_use]
    pub fn stats(&self) -> &FactorStats {
        &self.stats
    }

    /// Solve A·x = b. Panics on length mismatch (programmer error).
    #[must_use]
    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        assert_eq!(b.len(), self.n, "solve: rhs length must equal n");
        // y indexed by elimination-order id.
        let mut y = vec![0.0f64; self.n];
        for (k, &o) in self.perm.iter().enumerate() {
            y[k] = b[o];
        }
        // Forward: L z = P b.
        for s in &self.snodes {
            let w = s.block_vars.len();
            let m = w + s.panel_vars.len();
            let mut v: Vec<f64> = Vec::with_capacity(m);
            v.extend(s.block_vars.iter().map(|&i| y[i]));
            v.extend(s.panel_vars.iter().map(|&i| y[i]));
            for k in 0..w {
                let vk = v[k];
                for i in k + 1..m {
                    v[i] = s.l[i + k * m].mul_add(-vk, v[i]);
                }
            }
            for (k, &i) in s.block_vars.iter().enumerate() {
                y[i] = v[k];
            }
            for (k, &i) in s.panel_vars.iter().enumerate() {
                y[i] = v[w + k];
            }
        }
        // Diagonal: D w = z.
        for s in &self.snodes {
            let mut k = 0usize;
            for blk in &s.d {
                match *blk {
                    DBlock::Single(d) => {
                        y[s.block_vars[k]] /= d;
                        k += 1;
                    }
                    DBlock::Pair(d11, d21, d22) => {
                        let det = d11 * d22 - d21 * d21;
                        let z1 = y[s.block_vars[k]];
                        let z2 = y[s.block_vars[k + 1]];
                        y[s.block_vars[k]] = (d22 * z1 - d21 * z2) / det;
                        y[s.block_vars[k + 1]] = (d11 * z2 - d21 * z1) / det;
                        k += 2;
                    }
                }
            }
        }
        // Backward: Lᵀ x = w.
        for s in self.snodes.iter().rev() {
            let w = s.block_vars.len();
            let m = w + s.panel_vars.len();
            let mut v: Vec<f64> = Vec::with_capacity(m);
            v.extend(s.block_vars.iter().map(|&i| y[i]));
            v.extend(s.panel_vars.iter().map(|&i| y[i]));
            for k in (0..w).rev() {
                let mut acc = v[k];
                for i in k + 1..m {
                    acc = s.l[i + k * m].mul_add(-v[i], acc);
                }
                v[k] = acc;
            }
            for (k, &i) in s.block_vars.iter().enumerate() {
                y[i] = v[k];
            }
        }
        // Un-permute.
        let mut x = vec![0.0f64; self.n];
        for (k, &o) in self.perm.iter().enumerate() {
            x[o] = y[k];
        }
        x
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Coo;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    }

    /// tridiag(-1, 2, -1): eigenvalues 2 − 2cos(kπ/(n+1)), k = 1..n.
    fn tridiag(n: usize) -> Csr {
        let mut coo = Coo::new(n, n);
        for i in 0..n {
            coo.push(i, i, 2.0);
            if i > 0 {
                coo.push(i, i - 1, -1.0);
            }
            if i + 1 < n {
                coo.push(i, i + 1, -1.0);
            }
        }
        coo.assemble()
    }

    fn shifted(a: &Csr, sigma: f64) -> Csr {
        let n = a.nrows();
        let mut coo = Coo::new(n, n);
        for r in 0..n {
            let (cols, vals) = a.row(r);
            for (&c, &v) in cols.iter().zip(vals) {
                coo.push(r, c, v);
            }
        }
        for i in 0..n {
            coo.push(i, i, -sigma);
        }
        coo.assemble()
    }

    /// Random sparse symmetric matrix with a stored diagonal.
    fn random_sym(n: usize, extra_per_row: usize, seed: u64) -> Csr {
        let mut s = seed;
        let mut coo = Coo::new(n, n);
        for i in 0..n {
            coo.push(i, i, 4.0 + lcg(&mut s));
        }
        for i in 0..n {
            for _ in 0..extra_per_row {
                let j = ((lcg(&mut s) + 0.5) * n as f64) as usize % n;
                if i != j {
                    let v = lcg(&mut s);
                    coo.push(i, j, v);
                    coo.push(j, i, v);
                }
            }
        }
        coo.assemble()
    }

    fn factor(a: &Csr, ordering: DirectOrdering, opts: &LdltOptions) -> LdltFactor {
        let sym = SymbolicLdlt::analyze(a, ordering).expect("analyze");
        sym.factor(a, opts).expect("factor")
    }

    fn max_abs(a: &Csr) -> f64 {
        let mut m = 0.0f64;
        for r in 0..a.nrows() {
            let (_, vals) = a.row(r);
            for &v in vals {
                m = m.max(v.abs());
            }
        }
        m
    }

    /// Residual gate ‖A·x − b‖∞ ≤ tol·‖A‖·(‖x‖∞ + ‖b‖∞-ish).
    fn check_solve(a: &Csr, f: &LdltFactor, seed: u64, tol_scale: f64) {
        let n = a.nrows();
        let mut s = seed;
        let x_true: Vec<f64> = (0..n).map(|_| lcg(&mut s)).collect();
        let mut b = vec![0.0; n];
        a.spmv(&x_true, &mut b);
        let x = f.solve(&b);
        let mut ax = vec![0.0; n];
        a.spmv(&x, &mut ax);
        let xin: f64 = x.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
        let tol = tol_scale * (n as f64) * f64::EPSILON * max_abs(a) * xin.max(1.0);
        for i in 0..n {
            assert!(
                (ax[i] - b[i]).abs() <= tol,
                "solve residual row {i}: {} vs tol {tol}",
                (ax[i] - b[i]).abs()
            );
        }
    }

    #[test]
    fn etree_matches_hand_example() {
        // Davis, Direct Methods, ch. 4 style: pattern chosen so the tree is
        // hand-checkable. Lower entries: (1,0),(2,1),(3,0),(3,2).
        let mut coo = Coo::new(4, 4);
        for i in 0..4 {
            coo.push(i, i, 4.0);
        }
        for (r, c) in [(1usize, 0usize), (2, 1), (3, 0), (3, 2)] {
            coo.push(r, c, -1.0);
            coo.push(c, r, -1.0);
        }
        let a = coo.assemble();
        let iperm: Vec<usize> = (0..4).collect();
        let ladj = build_lower_adj(&a, &iperm);
        let parent = etree(&ladj);
        // col 0: first lower neighbor is row 1 → parent 1; col 1 → 2; col 2 → 3.
        assert_eq!(parent, vec![1, 2, 3, NONE]);
    }

    #[test]
    fn amd_is_a_permutation_on_assorted_graphs() {
        for a in [
            tridiag(17),
            random_sym(40, 3, 7),
            crate::tests::laplacian_2d(8),
        ] {
            let p = amd_order(&a);
            let n = a.nrows();
            let mut seen = vec![false; n];
            for &v in &p {
                assert!(v < n && !seen[v], "amd emitted a non-permutation");
                seen[v] = true;
            }
        }
    }

    #[test]
    fn amd_beats_natural_ordering_on_the_arrow_matrix() {
        // Arrow pointing the WRONG way: dense first row/column. Natural
        // ordering eliminates column 0 first and fills completely; any
        // minimum-degree ordering eliminates the pendant nodes first and
        // keeps L linear in n.
        let n = 200;
        let mut coo = Coo::new(n, n);
        for i in 0..n {
            coo.push(i, i, 4.0);
        }
        for j in 1..n {
            coo.push(0, j, -1.0);
            coo.push(j, 0, -1.0);
        }
        let a = coo.assemble();
        let nat = SymbolicLdlt::analyze(&a, DirectOrdering::Natural).unwrap();
        let amd = SymbolicLdlt::analyze(&a, DirectOrdering::Amd).unwrap();
        assert!(
            nat.nnz_l() > n * (n - 1) / 4,
            "natural ordering should fill catastrophically (got {})",
            nat.nnz_l()
        );
        assert!(
            amd.nnz_l() < 4 * n,
            "amd must keep the arrow sparse (got {})",
            amd.nnz_l()
        );
        println!(
            "{{\"suite\":\"fs-sparse-direct\",\"case\":\"amd-vs-natural-arrow\",\"nnz_l_natural\":{},\"nnz_l_amd\":{},\"verdict\":\"pass\"}}",
            nat.nnz_l(),
            amd.nnz_l()
        );
    }

    #[test]
    fn spd_tridiagonal_solve_and_inertia() {
        let n = 50;
        let a = tridiag(n);
        let f = factor(&a, DirectOrdering::Amd, &LdltOptions::default());
        assert_eq!(
            f.inertia(),
            Inertia {
                positive: n,
                negative: 0
            }
        );
        assert_eq!(f.stats().pivots_2x2, 0, "SPD must factor with 1x1 pivots");
        check_solve(&a, &f, 11, 8.0);
    }

    #[test]
    fn zero_diagonal_indefinite_needs_the_2x2_path() {
        // [[0,1],[1,0]] — eigenvalues ±1; no 1×1 pivot exists.
        let mut coo = Coo::new(2, 2);
        coo.push(0, 1, 1.0);
        coo.push(1, 0, 1.0);
        let a = coo.assemble();
        let f = factor(&a, DirectOrdering::Natural, &LdltOptions::default());
        assert_eq!(
            f.inertia(),
            Inertia {
                positive: 1,
                negative: 1
            }
        );
        assert_eq!(f.stats().pivots_2x2, 1);
        check_solve(&a, &f, 3, 8.0);

        // MUTATION CONTROL: with 2×2 pivots disabled, the same matrix must
        // REFUSE with the named breakdown — proving the 2×2 path is
        // load-bearing, not decorative.
        let sym = SymbolicLdlt::analyze(&a, DirectOrdering::Natural).unwrap();
        let opts = LdltOptions {
            allow_2x2: false,
            ..LdltOptions::default()
        };
        let err = sym.factor(&a, &opts).unwrap_err();
        assert!(
            matches!(err, LdltError::PivotBreakdown { .. }),
            "expected pivot breakdown, got {err}"
        );
        assert!(err.to_string().contains("FS-SPARSE-DIRECT-PIVOT-BREAKDOWN"));
    }

    #[test]
    fn saddle_point_kkt_inertia() {
        // [[A, Bᵀ],[B, 0]] with A = 2·I₃ SPD, B = [[1,0,0],[0,1,0]] full
        // row rank → inertia (3, 2) by the standard KKT result.
        let mut coo = Coo::new(5, 5);
        for i in 0..3 {
            coo.push(i, i, 2.0);
        }
        for (bi, &c) in [0usize, 1usize].iter().enumerate() {
            coo.push(3 + bi, c, 1.0);
            coo.push(c, 3 + bi, 1.0);
        }
        let a = coo.assemble();
        let f = factor(&a, DirectOrdering::Amd, &LdltOptions::default());
        assert_eq!(
            f.inertia(),
            Inertia {
                positive: 3,
                negative: 2
            }
        );
        check_solve(&a, &f, 5, 8.0);
    }

    #[test]
    fn inertia_certifies_analytic_eigenvalue_counts_across_shifts() {
        // THE spectrum-slicing seed property: for K = tridiag(-1,2,-1) the
        // eigenvalues are λ_k = 2 − 2cos(kπ/(n+1)), so the inertia of
        // K − σI must equal the analytic counts on both sides of σ.
        let n = 60;
        let k_mat = tridiag(n);
        let sym = SymbolicLdlt::analyze(&k_mat, DirectOrdering::Amd).unwrap();
        for sigma in [0.1, 0.5, 1.0, 2.0, 3.7] {
            let a = shifted(&k_mat, sigma);
            let f = sym.factor(&a, &LdltOptions::default()).expect("factor");
            let below = (1..=n)
                .filter(|&k| {
                    2.0 - 2.0 * (k as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos() < sigma
                })
                .count();
            assert_eq!(
                f.inertia(),
                Inertia {
                    positive: n - below,
                    negative: below
                },
                "inertia mismatch at sigma={sigma}"
            );
        }
        println!(
            "{{\"suite\":\"fs-sparse-direct\",\"case\":\"inertia-vs-analytic\",\"shifts\":5,\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn reconstruction_and_solve_on_random_symmetric_matrices() {
        for (n, extra, seed) in [(12usize, 2usize, 1u64), (40, 3, 2), (80, 5, 3)] {
            let a = random_sym(n, extra, seed);
            let f = factor(&a, DirectOrdering::Amd, &LdltOptions::default());
            assert_eq!(f.inertia().positive + f.inertia().negative, n);
            check_solve(&a, &f, seed ^ 0xABCD, 64.0);
        }
    }

    #[test]
    fn indefinite_random_matrices_solve_within_gate() {
        // Indefinite: flip diagonal signs on alternating rows.
        for (n, seed) in [(30usize, 21u64), (70, 22)] {
            let base = random_sym(n, 3, seed);
            let mut coo = Coo::new(n, n);
            for r in 0..n {
                let (cols, vals) = base.row(r);
                for (&c, &v) in cols.iter().zip(vals) {
                    if r == c && r % 2 == 1 {
                        coo.push(r, c, -v);
                    } else {
                        coo.push(r, c, v);
                    }
                }
            }
            let a = coo.assemble();
            let f = factor(&a, DirectOrdering::Amd, &LdltOptions::default());
            assert!(f.inertia().negative >= n / 4, "should be indefinite");
            check_solve(&a, &f, seed ^ 0x77, 64.0);
        }
    }

    #[test]
    fn factorization_is_bitwise_deterministic() {
        let a = random_sym(60, 4, 0xDE7);
        let sym = SymbolicLdlt::analyze(&a, DirectOrdering::Amd).unwrap();
        let f1 = sym.factor(&a, &LdltOptions::default()).unwrap();
        let f2 = sym.factor(&a, &LdltOptions::default()).unwrap();
        assert_eq!(f1.perm, f2.perm);
        for (s1, s2) in f1.snodes.iter().zip(&f2.snodes) {
            assert_eq!(s1.block_vars, s2.block_vars);
            assert!(
                s1.l.iter()
                    .zip(&s2.l)
                    .all(|(x, y)| x.to_bits() == y.to_bits()),
                "L values must be bitwise identical"
            );
        }
        let b: Vec<f64> = (0..60).map(|i| (i as f64) - 30.0).collect();
        let x1 = f1.solve(&b);
        let x2 = f2.solve(&b);
        assert!(x1.iter().zip(&x2).all(|(a, b)| a.to_bits() == b.to_bits()));
        println!(
            "{{\"suite\":\"fs-sparse-direct\",\"case\":\"bitwise-determinism\",\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn symbolic_reuse_rejects_pattern_mismatch() {
        let a = tridiag(10);
        let sym = SymbolicLdlt::analyze(&a, DirectOrdering::Amd).unwrap();
        // Same pattern, different values (a shift) → accepted.
        let shifted_a = shifted(&a, 0.25);
        assert!(sym.factor(&shifted_a, &LdltOptions::default()).is_ok());
        // Different pattern → typed refusal.
        let b = tridiag(11);
        assert_eq!(
            sym.factor(&b, &LdltOptions::default()).unwrap_err(),
            LdltError::PatternMismatch
        );
    }

    #[test]
    fn named_refusals_fire() {
        // Not square.
        let mut coo = Coo::new(2, 3);
        coo.push(0, 0, 1.0);
        let rect = coo.assemble();
        assert!(matches!(
            SymbolicLdlt::analyze(&rect, DirectOrdering::Amd),
            Err(LdltError::NotSquare { nrows: 2, ncols: 3 })
        ));
        // Structurally asymmetric.
        let mut coo = Coo::new(2, 2);
        coo.push(0, 0, 1.0);
        coo.push(1, 1, 1.0);
        coo.push(0, 1, 1.0);
        let asym = coo.assemble();
        assert!(matches!(
            SymbolicLdlt::analyze(&asym, DirectOrdering::Amd),
            Err(LdltError::StructurallyAsymmetric { row: 0, col: 1 })
        ));
        // Non-finite entry.
        let mut coo = Coo::new(2, 2);
        coo.push(0, 0, 1.0);
        coo.push(1, 1, f64::NAN);
        let nan = coo.assemble();
        let sym = SymbolicLdlt::analyze(&nan, DirectOrdering::Amd).unwrap();
        assert!(matches!(
            sym.factor(&nan, &LdltOptions::default()),
            Err(LdltError::NonFiniteEntry { row: 1, col: 1 })
        ));
        // Numerically singular (an exactly zero matrix with stored pattern).
        let mut coo = Coo::new(3, 3);
        for i in 0..3 {
            coo.push(i, i, 0.0);
        }
        let zero = coo.assemble();
        let sym = SymbolicLdlt::analyze(&zero, DirectOrdering::Amd).unwrap();
        assert!(matches!(
            sym.factor(&zero, &LdltOptions::default()),
            Err(LdltError::PivotBreakdown { .. })
        ));
        // Invalid options.
        let a = tridiag(3);
        let sym = SymbolicLdlt::analyze(&a, DirectOrdering::Amd).unwrap();
        let bad = LdltOptions {
            pivot_threshold: 0.9,
            ..LdltOptions::default()
        };
        assert!(matches!(
            sym.factor(&a, &bad),
            Err(LdltError::InvalidOptions { .. })
        ));
    }

    #[test]
    fn empty_and_singleton_matrices() {
        let empty = Coo::new(0, 0).assemble();
        let f = factor(&empty, DirectOrdering::Amd, &LdltOptions::default());
        assert_eq!(f.solve(&[]), Vec::<f64>::new());
        assert_eq!(
            f.inertia(),
            Inertia {
                positive: 0,
                negative: 0
            }
        );

        let mut coo = Coo::new(1, 1);
        coo.push(0, 0, -3.0);
        let one = coo.assemble();
        let f = factor(&one, DirectOrdering::Natural, &LdltOptions::default());
        assert_eq!(
            f.inertia(),
            Inertia {
                positive: 0,
                negative: 1
            }
        );
        let x = f.solve(&[6.0]);
        assert!((x[0] + 2.0).abs() < 1e-15);
    }

    #[test]
    fn grid_laplacian_factors_under_both_orderings_with_fill_recorded() {
        let a = crate::tests::laplacian_2d(16); // 256 unknowns
        let nat = SymbolicLdlt::analyze(&a, DirectOrdering::Natural).unwrap();
        let amd = SymbolicLdlt::analyze(&a, DirectOrdering::Amd).unwrap();
        let f = amd.factor(&a, &LdltOptions::default()).unwrap();
        assert_eq!(f.inertia().negative, 0, "grid Laplacian is SPD");
        check_solve(&a, &f, 31, 64.0);
        println!(
            "{{\"suite\":\"fs-sparse-direct\",\"case\":\"grid-fill\",\"n\":256,\"nnz_l_natural\":{},\"nnz_l_amd\":{},\"supernodes\":{},\"max_front\":{},\"verdict\":\"pass\"}}",
            nat.nnz_l(),
            amd.nnz_l(),
            f.stats().supernodes,
            f.stats().max_front_dim
        );
    }
}
