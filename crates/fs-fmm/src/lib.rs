//! fs-fmm — kernel-independent black-box FMM (plan §8.3 [F], bead
//! tfz.20): the Fong–Darve Chebyshev scheme. No analytic multipole
//! expansions anywhere — every translation operator is polynomial
//! interpolation, so ANY kernel that is smooth away from the diagonal
//! works unchanged, and the accuracy knob is one number: the
//! interpolation order p.
//!
//! Layer: L2. The pipeline: an octree over the point cloud (leaf
//! capacity bounded); UPWARD, source strengths anterpolate to each
//! leaf's p³ Chebyshev grid (P2M) and coarsen grid-to-grid (M2M);
//! ACROSS, well-separated cells exchange through the kernel evaluated
//! between their grids (M2L, the only place the kernel is touched);
//! DOWNWARD, locals refine grid-to-grid (L2L) and evaluate at targets
//! (L2P); adjacent cells run the direct sum (P2P). The conformance
//! battery measures accuracy against the direct oracle as p sweeps,
//! translation invariance (G3), and the scaling trend.
//!
//! Determinism: BTree-keyed octree, fixed traversal orders,
//! straight-line IEEE arithmetic.

use std::collections::BTreeMap;

/// Lowest supported Chebyshev interpolation order.
pub const MIN_INTERPOLATION_ORDER: usize = 2;
/// Highest admitted interpolation order for this unfused v1 implementation.
pub const MAX_INTERPOLATION_ORDER: usize = 12;
/// Hard point-count ceiling for one in-memory FMM instance.
pub const MAX_POINTS: usize = 1_048_576;
/// Conservative ceiling on simultaneously materialized per-cell coefficients.
pub const MAX_COEFFICIENT_SLOTS: usize = 16_777_216;
/// Ceiling on same-level cell-pair admission tests across one pass.
pub const MAX_CELL_PAIR_SCANS: usize = 33_554_432;
/// Conservative ceiling on M2L kernel-grid multiply-add work.
pub const MAX_TRANSLATION_WORK: usize = 4_000_000_000;
/// Largest admitted distribution-aware direct near-field pair count.
pub const MAX_NEAR_FIELD_PAIR_WORK: usize = 268_435_456;
/// Largest cloud admitted by the quadratic direct oracle.
pub const MAX_DIRECT_POINTS: usize = 16_384;
const MAX_ABS_COORDINATE: f64 = 1.0e100;

/// Invalid input or an inadmissible v1 work request.
#[derive(Debug, Clone, PartialEq)]
pub enum FmmError {
    /// A point cloud must contain at least one source/target.
    EmptyPointCloud,
    /// The point-count ceiling was exceeded.
    TooManyPoints {
        /// Requested points.
        count: usize,
        /// Admitted maximum.
        max: usize,
    },
    /// A point coordinate was non-finite or outside the stable bounding range.
    InvalidPoint {
        /// Point index.
        index: usize,
        /// Coordinate axis.
        axis: usize,
        /// Rejected coordinate.
        value: f64,
    },
    /// The interpolation order was outside the supported range.
    InvalidOrder {
        /// Requested order.
        order: usize,
        /// Supported minimum.
        min: usize,
        /// Supported maximum.
        max: usize,
    },
    /// Leaf capacity must be nonzero and no larger than the point ceiling.
    InvalidLeafCapacity {
        /// Requested capacity.
        capacity: usize,
        /// Admitted maximum.
        max: usize,
    },
    /// The conservative cell-grid work envelope was exceeded.
    WorkEnvelopeExceeded {
        /// Bounded work category.
        operation: &'static str,
        /// Conservative work count.
        requested: usize,
        /// Admitted maximum.
        max: usize,
    },
    /// There must be exactly one finite charge per point.
    InvalidCharges {
        /// Point count.
        expected: usize,
        /// Charge count.
        actual: usize,
        /// First non-finite charge, when lengths matched.
        non_finite_index: Option<usize>,
    },
    /// Kernel evaluation produced a non-finite potential.
    NonFiniteOutput {
        /// First invalid output index.
        index: usize,
    },
}

impl core::fmt::Display for FmmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyPointCloud => f.write_str("FMM point cloud must not be empty"),
            Self::TooManyPoints { count, max } => {
                write!(
                    f,
                    "FMM point count {count} exceeds the admitted maximum {max}"
                )
            }
            Self::InvalidPoint { index, axis, value } => write!(
                f,
                "FMM point {index} coordinate {axis} is invalid ({value}); coordinates must be finite with absolute value <= {MAX_ABS_COORDINATE}"
            ),
            Self::InvalidOrder { order, min, max } => {
                write!(
                    f,
                    "FMM order {order} is outside the supported range {min}..={max}"
                )
            }
            Self::InvalidLeafCapacity { capacity, max } => write!(
                f,
                "FMM leaf capacity {capacity} is invalid; expected 1..={max}"
            ),
            Self::WorkEnvelopeExceeded {
                operation,
                requested,
                max,
            } => write!(
                f,
                "FMM {operation} work {requested} exceeds the admitted maximum {max}"
            ),
            Self::InvalidCharges {
                expected,
                actual,
                non_finite_index,
            } => {
                if let Some(index) = non_finite_index {
                    write!(f, "FMM charge {index} is non-finite")
                } else {
                    write!(f, "FMM expected {expected} charges, received {actual}")
                }
            }
            Self::NonFiniteOutput { index } => {
                write!(
                    f,
                    "FMM kernel produced a non-finite potential at index {index}"
                )
            }
        }
    }
}

impl std::error::Error for FmmError {}

/// A smooth-off-diagonal kernel K(x, y).
pub trait Kernel {
    /// Evaluate K(target, source).
    fn eval(&self, x: [f64; 3], y: [f64; 3]) -> f64;
}

/// The 3D Laplace single-layer kernel 1/(4π|x−y|) (zero at the
/// diagonal by convention — self-interaction is the caller's story).
#[derive(Debug, Clone, Copy)]
pub struct Laplace3d;

impl Kernel for Laplace3d {
    fn eval(&self, x: [f64; 3], y: [f64; 3]) -> f64 {
        let d = ((x[0] - y[0]).powi(2) + (x[1] - y[1]).powi(2) + (x[2] - y[2]).powi(2)).sqrt();
        if d < 1e-300 {
            0.0
        } else {
            1.0 / (4.0 * std::f64::consts::PI * d)
        }
    }
}

/// Chebyshev nodes of the first kind on [-1, 1], order p.
fn cheb_nodes(p: usize) -> Vec<f64> {
    (0..p)
        .map(|k| {
            #[allow(clippy::cast_precision_loss)]
            let t = (2.0 * k as f64 + 1.0) / (2.0 * p as f64) * std::f64::consts::PI;
            -t.cos()
        })
        .collect()
}

/// Lagrange basis values at `x` over the Chebyshev nodes (stable
/// barycentric form, first-kind weights (−1)^k·sin term).
fn lagrange_at(nodes: &[f64], x: f64) -> Vec<f64> {
    let p = nodes.len();
    // Exact hit → delta.
    for (k, &n) in nodes.iter().enumerate() {
        if (x - n).abs() < 1e-14 {
            let mut v = vec![0.0; p];
            v[k] = 1.0;
            return v;
        }
    }
    let mut w = vec![0.0f64; p];
    for (k, wk) in w.iter_mut().enumerate() {
        #[allow(clippy::cast_precision_loss)]
        let t = (2.0 * k as f64 + 1.0) / (2.0 * p as f64) * std::f64::consts::PI;
        *wk = if k % 2 == 0 { 1.0 } else { -1.0 } * t.sin();
    }
    let mut num = vec![0.0f64; p];
    let mut den = 0.0;
    for k in 0..p {
        let c = w[k] / (x - nodes[k]);
        num[k] = c;
        den += c;
    }
    num.iter().map(|v| v / den).collect()
}

type CellKey = (u32, u32, u32, u32); // (level, i, j, k)

struct Cell {
    /// Point indices for leaves (empty for internal cells).
    points: Vec<usize>,
    /// Is this a leaf?
    leaf: bool,
}

/// The FMM engine for one point cloud (sources = targets, the BEM
/// matvec shape).
pub struct Fmm<'k> {
    kernel: &'k dyn Kernel,
    points: Vec<[f64; 3]>,
    order: usize,
    max_level: u32,
    cells: BTreeMap<CellKey, Cell>,
    /// Bounding cube (lo corner, side).
    lo: [f64; 3],
    side: f64,
}

impl<'k> Fmm<'k> {
    /// Validate a construction request without allocating the octree.
    ///
    /// The coefficient bound is conservative: every point is charged for one
    /// occupied cell at every tree level. The real tree can only contain fewer
    /// cells, so admission cannot underestimate grid storage.
    pub fn validate_request(
        points: &[[f64; 3]],
        order: usize,
        leaf_cap: usize,
    ) -> Result<(), FmmError> {
        if points.is_empty() {
            return Err(FmmError::EmptyPointCloud);
        }
        if points.len() > MAX_POINTS {
            return Err(FmmError::TooManyPoints {
                count: points.len(),
                max: MAX_POINTS,
            });
        }
        if !(MIN_INTERPOLATION_ORDER..=MAX_INTERPOLATION_ORDER).contains(&order) {
            return Err(FmmError::InvalidOrder {
                order,
                min: MIN_INTERPOLATION_ORDER,
                max: MAX_INTERPOLATION_ORDER,
            });
        }
        if leaf_cap == 0 || leaf_cap > MAX_POINTS {
            return Err(FmmError::InvalidLeafCapacity {
                capacity: leaf_cap,
                max: MAX_POINTS,
            });
        }
        for (index, point) in points.iter().enumerate() {
            for (axis, &value) in point.iter().enumerate() {
                if !value.is_finite() || value.abs() > MAX_ABS_COORDINATE {
                    return Err(FmmError::InvalidPoint { index, axis, value });
                }
            }
        }

        #[allow(clippy::cast_precision_loss)]
        let mut max_level =
            ((points.len() as f64 / leaf_cap as f64).max(1.0).ln() / (8.0f64).ln()).ceil() as usize;
        max_level = max_level.clamp(1, 6);
        let mut cell_bound = 0usize;
        let mut nonroot_cell_bound = 0usize;
        let mut pair_scans = 0usize;
        for level in 0..=max_level {
            let cells_at_level = 8usize
                .checked_pow(level as u32)
                .unwrap_or(usize::MAX)
                .min(points.len());
            cell_bound =
                cell_bound
                    .checked_add(cells_at_level)
                    .ok_or(FmmError::WorkEnvelopeExceeded {
                        operation: "cell-count",
                        requested: usize::MAX,
                        max: MAX_COEFFICIENT_SLOTS,
                    })?;
            if level > 0 {
                nonroot_cell_bound = nonroot_cell_bound.checked_add(cells_at_level).ok_or(
                    FmmError::WorkEnvelopeExceeded {
                        operation: "cell-count",
                        requested: usize::MAX,
                        max: MAX_COEFFICIENT_SLOTS,
                    },
                )?;
                pair_scans = pair_scans
                    .checked_add(
                        cells_at_level
                            .checked_mul(cells_at_level)
                            .unwrap_or(usize::MAX),
                    )
                    .ok_or(FmmError::WorkEnvelopeExceeded {
                        operation: "cell-pair scan",
                        requested: usize::MAX,
                        max: MAX_CELL_PAIR_SCANS,
                    })?;
            }
        }
        if pair_scans > MAX_CELL_PAIR_SCANS {
            return Err(FmmError::WorkEnvelopeExceeded {
                operation: "cell-pair scan",
                requested: pair_scans,
                max: MAX_CELL_PAIR_SCANS,
            });
        }
        let grid_size = order.checked_pow(3).ok_or(FmmError::WorkEnvelopeExceeded {
            operation: "coefficient-slot",
            requested: usize::MAX,
            max: MAX_COEFFICIENT_SLOTS,
        })?;
        let requested = cell_bound
            .checked_mul(grid_size)
            .and_then(|slots| slots.checked_mul(2))
            .ok_or(FmmError::WorkEnvelopeExceeded {
                operation: "coefficient-slot",
                requested: usize::MAX,
                max: MAX_COEFFICIENT_SLOTS,
            })?;
        if requested > MAX_COEFFICIENT_SLOTS {
            return Err(FmmError::WorkEnvelopeExceeded {
                operation: "coefficient-slot",
                requested,
                max: MAX_COEFFICIENT_SLOTS,
            });
        }
        let translation_work = nonroot_cell_bound
            .checked_mul(189)
            .and_then(|work| order.checked_pow(6).and_then(|p6| work.checked_mul(p6)))
            .ok_or(FmmError::WorkEnvelopeExceeded {
                operation: "M2L translation",
                requested: usize::MAX,
                max: MAX_TRANSLATION_WORK,
            })?;
        if translation_work > MAX_TRANSLATION_WORK {
            return Err(FmmError::WorkEnvelopeExceeded {
                operation: "M2L translation",
                requested: translation_work,
                max: MAX_TRANSLATION_WORK,
            });
        }
        Ok(())
    }

    /// Build the octree: UNIFORM depth chosen from N/leaf_cap (empty
    /// cells omitted) — on a uniform tree, "adjacent leaves run P2P,
    /// first-separated ancestors run M2L" partitions every pair
    /// EXACTLY ONCE (the adaptive U/V/W/X machinery is a recorded
    /// successor).
    ///
    /// # Errors
    /// Returns [`FmmError`] when geometry, tuning parameters, or the
    /// conservative work envelope are invalid.
    pub fn new(
        kernel: &'k dyn Kernel,
        points: Vec<[f64; 3]>,
        order: usize,
        leaf_cap: usize,
    ) -> Result<Fmm<'k>, FmmError> {
        Self::validate_request(&points, order, leaf_cap)?;
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in &points {
            for c in 0..3 {
                lo[c] = lo[c].min(p[c]);
                hi[c] = hi[c].max(p[c]);
            }
        }
        let side = (hi[0] - lo[0])
            .max(hi[1] - lo[1])
            .max(hi[2] - lo[2])
            .max(1e-12)
            * (1.0 + 1e-9);
        // Uniform depth: ~leaf_cap points per leaf for a uniform cloud.
        #[allow(clippy::cast_precision_loss)]
        let mut max_level =
            ((points.len() as f64 / leaf_cap as f64).max(1.0).ln() / (8.0f64).ln()).ceil() as u32;
        max_level = max_level.clamp(1, 6);
        let n_side = 1u32 << max_level;
        let child_side = side / f64::from(n_side);
        let mut cells: BTreeMap<CellKey, Cell> = BTreeMap::new();
        for (idx, p) in points.iter().enumerate() {
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let f = |c: usize| -> u32 {
                (((p[c] - lo[c]) / child_side).floor().max(0.0) as u32).min(n_side - 1)
            };
            let key = (max_level, f(0), f(1), f(2));
            cells
                .entry(key)
                .or_insert_with(|| Cell {
                    points: Vec::new(),
                    leaf: true,
                })
                .points
                .push(idx);
        }
        // Register ancestors of every nonempty leaf.
        let leaf_keys: Vec<CellKey> = cells.keys().copied().collect();
        let mut near_field_pairs = 0usize;
        for key in &leaf_keys {
            let target_count = cells[key].points.len();
            for di in -1i64..=1 {
                for dj in -1i64..=1 {
                    for dk in -1i64..=1 {
                        let (ni, nj, nk) = (
                            i64::from(key.1) + di,
                            i64::from(key.2) + dj,
                            i64::from(key.3) + dk,
                        );
                        if ni < 0 || nj < 0 || nk < 0 {
                            continue;
                        }
                        #[allow(clippy::cast_sign_loss)]
                        let neighbor = (key.0, ni as u32, nj as u32, nk as u32);
                        let Some(source_cell) = cells.get(&neighbor) else {
                            continue;
                        };
                        let pairs = target_count
                            .checked_mul(source_cell.points.len())
                            .unwrap_or(usize::MAX);
                        near_field_pairs = near_field_pairs.checked_add(pairs).ok_or(
                            FmmError::WorkEnvelopeExceeded {
                                operation: "near-field pair",
                                requested: usize::MAX,
                                max: MAX_NEAR_FIELD_PAIR_WORK,
                            },
                        )?;
                    }
                }
            }
        }
        if near_field_pairs > MAX_NEAR_FIELD_PAIR_WORK {
            return Err(FmmError::WorkEnvelopeExceeded {
                operation: "near-field pair",
                requested: near_field_pairs,
                max: MAX_NEAR_FIELD_PAIR_WORK,
            });
        }
        for key in leaf_keys {
            let (mut lv, mut i, mut j, mut k) = key;
            while lv > 0 {
                lv -= 1;
                i >>= 1;
                j >>= 1;
                k >>= 1;
                cells.entry((lv, i, j, k)).or_insert_with(|| Cell {
                    points: Vec::new(),
                    leaf: false,
                });
            }
        }
        Ok(Fmm {
            kernel,
            points,
            order,
            max_level,
            cells,
            lo,
            side,
        })
    }

    fn cell_box(&self, key: CellKey) -> ([f64; 3], f64) {
        let (lv, i, j, k) = key;
        let s = self.side / f64::from(1u32 << lv);
        (
            [
                self.lo[0] + f64::from(i) * s,
                self.lo[1] + f64::from(j) * s,
                self.lo[2] + f64::from(k) * s,
            ],
            s,
        )
    }

    /// The p³ Chebyshev grid of a cell (global coordinates).
    fn cell_grid(&self, key: CellKey) -> Vec<[f64; 3]> {
        let (blo, s) = self.cell_box(key);
        let nodes = cheb_nodes(self.order);
        let map = |t: f64, c: usize| blo[c] + 0.5 * s * (t + 1.0);
        let mut out = Vec::with_capacity(self.order.pow(3));
        for &tz in &nodes {
            for &ty in &nodes {
                for &tx in &nodes {
                    out.push([map(tx, 0), map(ty, 1), map(tz, 2)]);
                }
            }
        }
        out
    }

    /// Tensor Lagrange weights of a point in a cell's grid.
    fn interp_weights(&self, key: CellKey, p: [f64; 3]) -> Vec<f64> {
        let (blo, s) = self.cell_box(key);
        let nodes = cheb_nodes(self.order);
        let local = |c: usize| ((p[c] - blo[c]) / s * 2.0 - 1.0).clamp(-1.0, 1.0);
        let lx = lagrange_at(&nodes, local(0));
        let ly = lagrange_at(&nodes, local(1));
        let lz = lagrange_at(&nodes, local(2));
        let mut out = Vec::with_capacity(self.order.pow(3));
        for wz in &lz {
            for wy in &ly {
                for wx in &lx {
                    out.push(wx * wy * wz);
                }
            }
        }
        out
    }

    /// Are two same-level cells adjacent (or identical)?
    fn adjacent(a: CellKey, b: CellKey) -> bool {
        a.0 == b.0 && a.1.abs_diff(b.1) <= 1 && a.2.abs_diff(b.2) <= 1 && a.3.abs_diff(b.3) <= 1
    }

    /// Evaluate potentials at every point: FMM (upward, M2L, downward,
    /// near-field direct).
    #[allow(clippy::too_many_lines)] // the four classic passes, one narrative
    pub fn potentials(&self, charges: &[f64]) -> Result<Vec<f64>, FmmError> {
        self.validate_charges(charges)?;
        let np = self.order.pow(3);
        // UPWARD: P2M at leaves, M2M to parents.
        let mut multipole: BTreeMap<CellKey, Vec<f64>> = BTreeMap::new();
        for (&key, cell) in self.cells.iter().rev() {
            let mut m = vec![0.0f64; np];
            if cell.leaf {
                for &idx in &cell.points {
                    let w = self.interp_weights(key, self.points[idx]);
                    for (mi, wi) in m.iter_mut().zip(&w) {
                        *mi += charges[idx] * wi;
                    }
                }
            } else {
                // M2M: children's grids anterpolate into this grid.
                for di in 0..2u32 {
                    for dj in 0..2u32 {
                        for dk in 0..2u32 {
                            let child = (key.0 + 1, 2 * key.1 + di, 2 * key.2 + dj, 2 * key.3 + dk);
                            let Some(cm) = multipole.get(&child) else {
                                continue;
                            };
                            let grid = self.cell_grid(child);
                            for (g, &cmv) in grid.iter().zip(cm) {
                                if cmv == 0.0 {
                                    continue;
                                }
                                let w = self.interp_weights(key, *g);
                                for (mi, wi) in m.iter_mut().zip(&w) {
                                    *mi += cmv * wi;
                                }
                            }
                        }
                    }
                }
            }
            multipole.insert(key, m);
        }
        // M2L: per cell, the interaction list = same-level cells whose
        // parents are adjacent but which are not themselves adjacent.
        let mut local: BTreeMap<CellKey, Vec<f64>> = BTreeMap::new();
        for &key in self.cells.keys() {
            local.insert(key, vec![0.0f64; np]);
        }
        let keys: Vec<CellKey> = self.cells.keys().copied().collect();
        // Group by level for interaction-list construction.
        let mut by_level: BTreeMap<u32, Vec<CellKey>> = BTreeMap::new();
        for &k in &keys {
            by_level.entry(k.0).or_default().push(k);
        }
        for (&lv, level_keys) in &by_level {
            if lv == 0 {
                continue;
            }
            for &a in level_keys {
                let pa = (a.0 - 1, a.1 >> 1, a.2 >> 1, a.3 >> 1);
                let ga = self.cell_grid(a);
                let la = local.get_mut(&a).expect("local").clone();
                let mut acc = la;
                for &b in level_keys {
                    if Self::adjacent(a, b) {
                        continue;
                    }
                    let pb = (b.0 - 1, b.1 >> 1, b.2 >> 1, b.3 >> 1);
                    if !Self::adjacent(pa, pb) {
                        continue;
                    }
                    let mb = &multipole[&b];
                    let gb = self.cell_grid(b);
                    for (ti, tp) in ga.iter().enumerate() {
                        let mut s = 0.0;
                        for (sj, sp) in gb.iter().enumerate() {
                            if mb[sj] != 0.0 {
                                s += self.kernel.eval(*tp, *sp) * mb[sj];
                            }
                        }
                        acc[ti] += s;
                    }
                }
                local.insert(a, acc);
            }
        }
        // DOWNWARD: L2L parent → child, then L2P + near-field P2P.
        let mut out = vec![0.0f64; self.points.len()];
        for &key in &keys {
            if key.0 > 0 {
                let parent = (key.0 - 1, key.1 >> 1, key.2 >> 1, key.3 >> 1);
                let lp = local[&parent].clone();
                let grid = self.cell_grid(key);
                let mut add = vec![0.0f64; np];
                let pgrid_weights: Vec<Vec<f64>> = grid
                    .iter()
                    .map(|g| self.interp_weights(parent, *g))
                    .collect();
                for (ti, w) in pgrid_weights.iter().enumerate() {
                    let mut s = 0.0;
                    for (wi, lv) in w.iter().zip(&lp) {
                        s += wi * lv;
                    }
                    add[ti] = s;
                }
                let l = local.get_mut(&key).expect("local");
                for (li, ai) in l.iter_mut().zip(&add) {
                    *li += ai;
                }
            }
            let cell = &self.cells[&key];
            if !cell.leaf {
                continue;
            }
            // L2P.
            let l = &local[&key];
            for &idx in &cell.points {
                let w = self.interp_weights(key, self.points[idx]);
                let mut s = 0.0;
                for (wi, li) in w.iter().zip(l) {
                    s += wi * li;
                }
                out[idx] += s;
            }
            // P2P: adjacent leaves (uniform depth ⇒ same level; the
            // M2L lists cover everything else exactly once).
            for di in -1i64..=1 {
                for dj in -1i64..=1 {
                    for dk in -1i64..=1 {
                        let (ni, nj, nk) = (
                            i64::from(key.1) + di,
                            i64::from(key.2) + dj,
                            i64::from(key.3) + dk,
                        );
                        if ni < 0 || nj < 0 || nk < 0 {
                            continue;
                        }
                        #[allow(clippy::cast_sign_loss)]
                        let okey = (key.0, ni as u32, nj as u32, nk as u32);
                        let Some(ocell) = self.cells.get(&okey) else {
                            continue;
                        };
                        if !ocell.leaf {
                            continue;
                        }
                        for &t in &cell.points {
                            let mut s = 0.0;
                            for &sidx in &ocell.points {
                                if sidx != t {
                                    s += self.kernel.eval(self.points[t], self.points[sidx])
                                        * charges[sidx];
                                }
                            }
                            out[t] += s;
                        }
                    }
                }
            }
        }
        if let Some(index) = out.iter().position(|value| !value.is_finite()) {
            return Err(FmmError::NonFiniteOutput { index });
        }
        Ok(out)
    }

    /// The direct O(N²) oracle.
    pub fn direct(&self, charges: &[f64]) -> Result<Vec<f64>, FmmError> {
        self.validate_charges(charges)?;
        let n = self.points.len();
        if n > MAX_DIRECT_POINTS {
            return Err(FmmError::WorkEnvelopeExceeded {
                operation: "direct pair",
                requested: n.checked_mul(n).unwrap_or(usize::MAX),
                max: MAX_DIRECT_POINTS * MAX_DIRECT_POINTS,
            });
        }
        let mut out = vec![0.0f64; n];
        for (t, slot) in out.iter_mut().enumerate() {
            let mut s = 0.0;
            for (sidx, &q) in charges.iter().enumerate() {
                if sidx != t {
                    s += self.kernel.eval(self.points[t], self.points[sidx]) * q;
                }
            }
            *slot = s;
        }
        if let Some(index) = out.iter().position(|value| !value.is_finite()) {
            return Err(FmmError::NonFiniteOutput { index });
        }
        Ok(out)
    }

    fn validate_charges(&self, charges: &[f64]) -> Result<(), FmmError> {
        if charges.len() != self.points.len() {
            return Err(FmmError::InvalidCharges {
                expected: self.points.len(),
                actual: charges.len(),
                non_finite_index: None,
            });
        }
        if let Some(index) = charges.iter().position(|charge| !charge.is_finite()) {
            return Err(FmmError::InvalidCharges {
                expected: self.points.len(),
                actual: charges.len(),
                non_finite_index: Some(index),
            });
        }
        Ok(())
    }

    /// Octree statistics (ledger row).
    #[must_use]
    pub fn stats(&self) -> String {
        let leaves = self.cells.values().filter(|c| c.leaf).count();
        format!(
            "{{\"cells\":{},\"leaves\":{leaves},\"max_level\":{},\"order\":{}}}",
            self.cells.len(),
            self.max_level,
            self.order
        )
    }
}
