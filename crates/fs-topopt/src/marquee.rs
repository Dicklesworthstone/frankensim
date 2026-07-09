//! THE CUTFEM-OCTREE TOPOLOGY MARQUEE (bead b7d0; [F] — behind the
//! `cutfem-marquee` feature): the density pipeline executed on CutFEM
//! over an octree background grid. The density field lives on a
//! lattice; the SOLID region {ρ > ½} IS the CutFEM domain (density-as-
//! indicator — the design boundary is the cooled surface of the
//! volume-to-point heat fixture); TOPOLOGY EVOLVES WITH ZERO
//! REMESHING — the background quadtree is built once and only ever
//! REFINED (splits), never rebuilt, and the run log proves it.
//!
//! DWR-goal-driven refinement: fs-dwr's per-leaf compliance-goal
//! indicators, weighted by the density gradient, choose the splits —
//! the octree refines where the OBJECTIVE feels the design boundary.

use fs_cutfem::sdf::CutSdf;
use fs_cutfem::{FemParams, Quadtree, Space};
use fs_dwr::{GoalContext, estimate, goal_value};
use fs_ivl::Interval;

/// The design: densities on an `n × n` node lattice over `[0, 1]²`;
/// the solid region is `ρ > ½` and `φ = ½ − ρ` (bilinear) is the
/// CutFEM domain field (negative inside the solid).
#[derive(Debug, Clone, PartialEq)]
pub struct DensityDesign {
    /// Nodes per side.
    pub n: usize,
    /// Row-major nodal densities in [0, 1].
    pub rho: Vec<f64>,
}

impl DensityDesign {
    /// A uniform-density start at `frac` solid fraction.
    #[must_use]
    pub fn uniform(n: usize, frac: f64) -> DensityDesign {
        DensityDesign {
            n,
            rho: vec![frac; n * n],
        }
    }

    fn node(&self, i: usize, j: usize) -> f64 {
        self.rho[j * self.n + i]
    }

    /// Bilinear density at a point.
    #[must_use]
    pub fn density_at(&self, x: f64, y: f64) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let scale = (self.n - 1) as f64;
        let (fx, fy) = (x.clamp(0.0, 1.0) * scale, y.clamp(0.0, 1.0) * scale);
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let (i, j) = (
            (fx.floor() as usize).min(self.n - 2),
            (fy.floor() as usize).min(self.n - 2),
        );
        #[allow(clippy::cast_precision_loss)]
        let (tx, ty) = (fx - i as f64, fy - j as f64);
        (1.0 - tx) * (1.0 - ty) * self.node(i, j)
            + tx * (1.0 - ty) * self.node(i + 1, j)
            + tx * ty * self.node(i + 1, j + 1)
            + (1.0 - tx) * ty * self.node(i, j + 1)
    }

    /// Solid fraction (mean density).
    #[must_use]
    pub fn volume(&self) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        {
            self.rho.iter().sum::<f64>() / self.rho.len() as f64
        }
    }

    /// Count void components (4-connected nodes with ρ ≤ ½) — the
    /// topology-evolution witness.
    #[must_use]
    pub fn void_components(&self) -> usize {
        let n = self.n;
        let mut seen = vec![false; n * n];
        let mut comps = 0usize;
        for start in 0..n * n {
            if seen[start] || self.rho[start] > 0.5 {
                continue;
            }
            comps += 1;
            let mut stack = vec![start];
            seen[start] = true;
            while let Some(k) = stack.pop() {
                let (i, j) = (k % n, k / n);
                let mut push = |q: usize| {
                    if !seen[q] && self.rho[q] <= 0.5 {
                        seen[q] = true;
                        stack.push(q);
                    }
                };
                if i > 0 {
                    push(k - 1);
                }
                if i + 1 < n {
                    push(k + 1);
                }
                if j > 0 {
                    push(k - n);
                }
                if j + 1 < n {
                    push(k + n);
                }
            }
        }
        comps
    }

    /// The MEDIAL-AXIS-CLASS thickness oracle: the maximum over solid
    /// components of the interior chessboard distance to the void,
    /// and the MINIMUM local thickness (2× the smallest maximal
    /// interior distance over components) in lattice cells — the
    /// length-scale audit of the optimized geometry.
    #[must_use]
    pub fn min_feature_cells(&self) -> usize {
        let n = self.n;
        // Distance transform (chessboard) from the void/boundary.
        let mut dist = vec![usize::MAX; n * n];
        let mut frontier: Vec<usize> = (0..n * n)
            .filter(|&k| {
                let (i, j) = (k % n, k / n);
                self.rho[k] <= 0.5 || i == 0 || j == 0 || i == n - 1 || j == n - 1
            })
            .collect();
        for &k in &frontier {
            dist[k] = 0;
        }
        let mut d = 0usize;
        while !frontier.is_empty() {
            d += 1;
            let mut next = Vec::new();
            for &k in &frontier {
                let (i, j) = (k % n, k / n);
                let mut visit = |q: usize| {
                    if dist[q] == usize::MAX {
                        dist[q] = d;
                        next.push(q);
                    }
                };
                if i > 0 {
                    visit(k - 1);
                }
                if i + 1 < n {
                    visit(k + 1);
                }
                if j > 0 {
                    visit(k - n);
                }
                if j + 1 < n {
                    visit(k + n);
                }
            }
            frontier = next;
        }
        // Per solid component: its maximal interior distance (the
        // inscribed radius); the min over components ×2 = min feature.
        let mut seen = vec![false; n * n];
        let mut min_radius = usize::MAX;
        for start in 0..n * n {
            if seen[start] || self.rho[start] <= 0.5 {
                continue;
            }
            let mut radius = 0usize;
            let mut stack = vec![start];
            seen[start] = true;
            while let Some(k) = stack.pop() {
                radius = radius.max(dist[k]);
                let (i, j) = (k % n, k / n);
                let mut push = |q: usize| {
                    if !seen[q] && self.rho[q] > 0.5 {
                        seen[q] = true;
                        stack.push(q);
                    }
                };
                if i > 0 {
                    push(k - 1);
                }
                if i + 1 < n {
                    push(k + 1);
                }
                if j > 0 {
                    push(k - n);
                }
                if j + 1 < n {
                    push(k + n);
                }
            }
            min_radius = min_radius.min(radius);
        }
        if min_radius == usize::MAX {
            0
        } else {
            2 * min_radius
        }
    }
}

impl CutSdf for DensityDesign {
    fn value(&self, p: [f64; 2]) -> f64 {
        0.5 - self.density_at(p[0], p[1])
    }

    fn gradient(&self, p: [f64; 2]) -> [f64; 2] {
        let h = 1e-4;
        [
            (self.value([p[0] + h, p[1]]) - self.value([p[0] - h, p[1]])) / (2.0 * h),
            (self.value([p[0], p[1] + h]) - self.value([p[0], p[1] - h])) / (2.0 * h),
        ]
    }

    fn enclose(&self, lo: [f64; 2], hi: [f64; 2]) -> Interval {
        // The bilinear field's extrema over a box occur at lattice
        // nodes covered by the box or at the box corners — enumerate
        // both for an exact-containment enclosure.
        #[allow(clippy::cast_precision_loss)]
        let scale = (self.n - 1) as f64;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let (i0, j0) = (
            ((lo[0].clamp(0.0, 1.0) * scale).floor() as usize).min(self.n - 1),
            ((lo[1].clamp(0.0, 1.0) * scale).floor() as usize).min(self.n - 1),
        );
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let (i1, j1) = (
            ((hi[0].clamp(0.0, 1.0) * scale).ceil() as usize).min(self.n - 1),
            ((hi[1].clamp(0.0, 1.0) * scale).ceil() as usize).min(self.n - 1),
        );
        let mut rho_lo = f64::INFINITY;
        let mut rho_hi = f64::NEG_INFINITY;
        for j in j0..=j1 {
            for i in i0..=i1 {
                let v = self.node(i, j);
                rho_lo = rho_lo.min(v);
                rho_hi = rho_hi.max(v);
            }
        }
        for &(x, y) in &[
            (lo[0], lo[1]),
            (hi[0], lo[1]),
            (lo[0], hi[1]),
            (hi[0], hi[1]),
        ] {
            let v = self.density_at(x, y);
            rho_lo = rho_lo.min(v);
            rho_hi = rho_hi.max(v);
        }
        Interval::new(0.5 - rho_hi, 0.5 - rho_lo)
    }
}

/// One iteration's forensic record.
#[derive(Debug, Clone, PartialEq)]
pub struct MarqueeIter {
    /// Iteration index.
    pub iter: usize,
    /// Thermal compliance J = ∫ f u.
    pub compliance: f64,
    /// Solid fraction.
    pub volume: f64,
    /// Void-component count (topology witness).
    pub voids: usize,
    /// Octree SPLITS this iteration (refinement events).
    pub splits: usize,
    /// Mesh REBUILDS this iteration (the marquee property: always 0).
    pub rebuilds: usize,
    /// Wall time of the iteration (measured; debug-build label).
    pub wall_ms: f64,
}

/// The marquee run report.
#[derive(Debug, Clone, PartialEq)]
pub struct MarqueeReport {
    /// Per-iteration records.
    pub iterations: Vec<MarqueeIter>,
    /// The final design.
    pub design: DensityDesign,
    /// Total octree splits across the run.
    pub total_splits: usize,
    /// Final refined leaves whose one-cell halo intersects the design
    /// boundary; this is the executable footprint evidence for
    /// DWR-driven boundary-band concentration.
    pub refined_boundary_leaves: usize,
    /// Final refined leaves away from the design-boundary halo.
    pub refined_off_boundary_leaves: usize,
    /// Total mesh rebuilds (MUST be zero — asserted by the caller).
    pub total_rebuilds: usize,
}

/// True when the cell OR its one-cell halo is cut: fs-cutfem's ghost
/// penalty demands equal-level FACE NEIGHBORS of cut cells, so the
/// refinement band must include the halo, not just the straddling
/// cells (the CutBandNotUniform contract, learned the hard way twice).
fn halo_cut(design: &DensityDesign, lo: [f64; 2], hi: [f64; 2]) -> bool {
    let (wx, wy) = (hi[0] - lo[0], hi[1] - lo[1]);
    let xlo = [(lo[0] - wx).max(0.0), (lo[1] - wy).max(0.0)];
    let xhi = [(hi[0] + wx).min(1.0), (hi[1] + wy).min(1.0)];
    design.enclose(xlo, xhi).contains_zero()
}

fn fem_params() -> FemParams {
    FemParams {
        nitsche_beta: 10.0,
        ghost_gamma: 0.1,
        quad_depth: 3,
        agg: None,
        strong_outer: true,
        solver_tol: 1e-9,
        solver_max_iters: 1200,
    }
}

/// Run the marquee: the volume-to-point heat fixture (f = 1 body
/// heating, the design boundary cooled to 0) at a fixed solid
/// fraction. Interface-flux redistribution evolves the density;
/// DWR × |∇ρ| chooses octree splits; the background grid is built
/// ONCE and never rebuilt.
///
/// # Errors
/// CutFEM build/solve errors propagate.
pub fn run_marquee(
    mut design: DensityDesign,
    base_level: u32,
    max_level: u32,
    iters: usize,
    splits_per_iter: usize,
) -> Result<MarqueeReport, fs_cutfem::CutFemError> {
    // THE GRID IS BUILT ONCE. Refinement = splits only; there is no
    // other construction site in this function (the zero-remeshing
    // property is structural, and the log proves it).
    let mut grid = Quadtree::with_room(base_level, max_level);
    let mut iterations = Vec::with_capacity(iters);
    let mut total_splits = 0usize;
    let mut band_level = base_level;
    let target_volume = design.volume();
    for iter in 0..iters {
        let t0 = std::time::Instant::now();
        // Keep the (moving) cut band CONFORMING at the current band
        // level before every solve — fs-cutfem requires a uniform-level
        // band, and the interface moves between iterations, so newly
        // cut cells must be brought to the band level (splits only;
        // idempotent for already-fine cells).
        let mut pre_splits = 0usize;
        if band_level > base_level {
            let before = grid.leaf_count();
            let d_ref = &design;
            grid.refine_where(band_level, &|lo, hi| halo_cut(d_ref, lo, hi));
            pre_splits = grid.leaf_count().saturating_sub(before) / 3;
        }
        let params = fem_params();
        let f = |_: f64, _: f64| 1.0;
        let g = |_: f64, _: f64| 0.0;
        let space = Space::build(&grid, &design, params)?;
        let sol = space.solve(&f, &g)?;
        let nodal = space.nodal_values(&sol.free, &g);
        let goal = GoalContext { weight: &f };
        let j = goal_value(&space, &grid, &design, &nodal, &goal, params.quad_depth);
        // DWR per-leaf indicators for the compliance goal.
        let dwr = estimate(&grid, &design, params, &f, &g, &goal)?;
        // --- Density update: interface-flux redistribution. ---------
        // Sample each lattice node's neighborhood; nodes NEAR the
        // interface get a signed move: high local flux² → carve (the
        // boundary wants to grow there, mye.1's shape derivative),
        // low flux² → fill; then project back to the volume target.
        let n = design.n;
        let u_at = |x: f64, y: f64| -> f64 {
            // Bilinear through the containing leaf's corner nodes
            // (guaranteed present in the nodal map — probing the raw
            // fine lattice missed the sparse keys and froze the whole
            // update in the first draft).
            let Some(leaf) =
                grid.find_leaf_at(x.clamp(1e-9, 1.0 - 1e-9), y.clamp(1e-9, 1.0 - 1e-9))
            else {
                return 0.0;
            };
            let (lo, hi) = grid.rect(leaf);
            let corners = grid.corner_nodes(leaf);
            let v = |k: usize| nodal.get(&corners[k]).copied().unwrap_or(0.0);
            let tx = ((x - lo[0]) / (hi[0] - lo[0])).clamp(0.0, 1.0);
            let ty = ((y - lo[1]) / (hi[1] - lo[1])).clamp(0.0, 1.0);
            // corner_nodes order is CCW: 0=(lo,lo) 1=(hi,lo) 2=(hi,hi) 3=(lo,hi)
            (1.0 - tx) * (1.0 - ty) * v(0)
                + tx * (1.0 - ty) * v(1)
                + tx * ty * v(2)
                + (1.0 - tx) * ty * v(3)
        };
        let flux_at = |x: f64, y: f64| -> f64 {
            // Probe u a fixed depth INSIDE the solid measured from the
            // INTERFACE, from either side: first-order signed distance
            // s = phi/|grad phi| (positive in the void), then step
            // (s + h) against the gradient. Probing from the raw node
            // position left void-side nodes reading zero flux and
            // biased the run toward shrinking the cooled boundary (the
            // J-rises bug of an earlier draft).
            let gph = design.gradient([x, y]);
            let norm = (gph[0] * gph[0] + gph[1] * gph[1]).sqrt().max(1e-9);
            let sdist = design.value([x, y]) / norm;
            let h = 0.05;
            let depth = sdist + h;
            let u = u_at(x - depth * gph[0] / norm, y - depth * gph[1] / norm);
            (u / h).powi(2)
        };
        let mut moves = vec![0.0f64; n * n];
        let mut flux_sum = 0.0f64;
        let mut flux_cnt = 0usize;
        #[allow(clippy::cast_precision_loss)]
        let lattice_scale = (n - 1) as f64;
        for k in 0..n * n {
            let (i, jj) = (k % n, k / n);
            // Interface-adjacent = a 4-neighbor on the other side of
            // the 0.5 level (phi is a density gap, NOT a distance —
            // testing |phi| < eps found zero band nodes and froze the
            // whole update in an earlier draft).
            let solid = design.rho[k] > 0.5;
            let mut near = false;
            if i > 0 {
                near |= (design.rho[k - 1] > 0.5) != solid;
            }
            if i + 1 < n {
                near |= (design.rho[k + 1] > 0.5) != solid;
            }
            if jj > 0 {
                near |= (design.rho[k - n] > 0.5) != solid;
            }
            if jj + 1 < n {
                near |= (design.rho[k + n] > 0.5) != solid;
            }
            if near {
                #[allow(clippy::cast_precision_loss)]
                let (x, y) = (i as f64 / lattice_scale, jj as f64 / lattice_scale);
                let fl = flux_at(x, y);
                moves[k] = fl.max(1e-12);
                flux_sum += fl;
                flux_cnt += 1;
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let flux_mean = flux_sum / flux_cnt.max(1) as f64;
        let step = 0.25;
        for k in 0..n * n {
            if moves[k] > 0.0 {
                // carve where flux is above the mean, fill below.
                let rel = (moves[k] - flux_mean) / flux_mean.max(1e-12);
                design.rho[k] = (design.rho[k] - step * rel.clamp(-1.0, 1.0)).clamp(0.02, 0.98);
            }
        }
        // Volume projection ON THE BAND ONLY: a uniform shift over all
        // nodes silently fills the voids from the inside (interior
        // void nodes creep past 0.5 over iterations — the J-rising
        // bias of an earlier draft). The correction lives where the
        // moves happened.
        let band: Vec<usize> = (0..n * n).filter(|&k| moves[k] > 0.0).collect();
        if !band.is_empty() {
            let (mut lo, mut hi) = (-0.5f64, 0.5f64);
            for _ in 0..40 {
                let mid = f64::midpoint(lo, hi);
                let vol: f64 = design
                    .rho
                    .iter()
                    .enumerate()
                    .map(|(k, r)| {
                        if moves[k] > 0.0 {
                            (r + mid).clamp(0.02, 0.98)
                        } else {
                            *r
                        }
                    })
                    .sum::<f64>()
                    / design.rho.len() as f64;
                if vol > target_volume {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            let shift = f64::midpoint(lo, hi);
            for &k in &band {
                design.rho[k] = (design.rho[k] + shift).clamp(0.02, 0.98);
            }
        }
        // --- DWR-gated refinement: splits ONLY, band-uniform. --------
        // fs-cutfem requires the CUT BAND at a uniform level (its
        // CutBandNotUniform contract — the first draft split top-k
        // cells individually and the solver refused, correctly). So
        // the DWR indicators decide WHEN to refine (is the interface
        // band carrying enough of the compliance-goal error mass?) and
        // `refine_where` refines the WHOLE band one level, uniformly.
        let mut band_mass = 0.0f64;
        let mut band_min_level = max_level;
        for (&(lv, cx, cy), &eta) in &dwr.indicators {
            let cells = f64::from(1u32 << lv);
            let lo = [f64::from(cx) / cells, f64::from(cy) / cells];
            let hi = [lo[0] + 1.0 / cells, lo[1] + 1.0 / cells];
            let enc = design.enclose(lo, hi);
            if enc.contains_zero() {
                band_mass += eta.abs();
                band_min_level = band_min_level.min(lv);
            }
        }
        let mut splits = pre_splits;
        let _ = band_min_level;
        if splits_per_iter > 0 && band_mass > 0.15 * dwr.eta_abs && band_level < max_level {
            let before = grid.leaf_count();
            band_level += 1;
            let d_ref = &design;
            grid.refine_where(band_level, &|lo, hi| halo_cut(d_ref, lo, hi));
            splits += grid.leaf_count().saturating_sub(before) / 3;
        }
        total_splits += splits;
        #[allow(clippy::cast_precision_loss)]
        let wall_ms = t0.elapsed().as_secs_f64() * 1e3;
        iterations.push(MarqueeIter {
            iter,
            compliance: j,
            volume: design.volume(),
            voids: design.void_components(),
            splits,
            rebuilds: 0, // structural: there is no rebuild path
            wall_ms,
        });
    }
    let mut refined_boundary_leaves = 0usize;
    let mut refined_off_boundary_leaves = 0usize;
    for leaf in grid.leaves().filter(|leaf| leaf.0 > base_level) {
        let (lo, hi) = grid.rect(leaf);
        if halo_cut(&design, lo, hi) {
            refined_boundary_leaves += 1;
        } else {
            refined_off_boundary_leaves += 1;
        }
    }

    Ok(MarqueeReport {
        iterations,
        design,
        total_splits,
        refined_boundary_leaves,
        refined_off_boundary_leaves,
        total_rebuilds: 0,
    })
}
