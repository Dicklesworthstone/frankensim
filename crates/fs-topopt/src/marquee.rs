//! THE CUTFEM-QUADTREE TOPOLOGY MARQUEE (the 2D analogue of the planned
//! octree lane; bead b7d0; [F] — behind the
//! `cutfem-marquee` feature): the density pipeline executed on CutFEM
//! over a quadtree background grid. The density field lives on a
//! lattice; the SOLID region {ρ > ½} IS the CutFEM domain (density-as-
//! indicator — the design boundary is the cooled surface of the
//! volume-to-point heat fixture); TOPOLOGY EVOLVES WITH ZERO
//! REMESHING — the background quadtree is built once and only ever
//! REFINED (splits), never rebuilt, and the run log proves it.
//!
//! DWR-goal-driven refinement: fs-dwr's per-leaf compliance-goal
//! indicators gate one-level refinement of the cut band and its ghost-
//! penalty halo — the quadtree refines when the OBJECTIVE's estimated
//! error mass is concentrated on exactly that refinement target.

use fs_cutfem::sdf::CutSdf;
use fs_cutfem::{CellKey, FemParams, Quadtree, ScalarSample, Space};
use fs_dwr::{GoalContext, estimate, goal_value};
use fs_ivl::Interval;
use std::cell::RefCell;
use std::collections::BTreeMap;

/// Version of the estimator-time mass gate / next-design band policy.
///
/// Bump this whenever the mass partition, strict comparison, level-advance
/// semantics, or spatial target changes; reports carry the version so old
/// evidence cannot be mistaken for evidence from a revised policy.
///
/// Version 2 measures mass on the region the policy actually refines (cut
/// cells plus their one-cell ghost halo, [`in_dwr_refinement_band`]) and adds
/// a concentration clause. Version 1 counted only zero-straddling cells, so
/// it measured a different region than it refined: on the seeded marquee at
/// level 4 the cut cells held 11% of the mass on 10% of the cells, while
/// cut + halo held 48% on 26%, and the largest indicators sat one cell off
/// the interface. The gate never fired and the marquee never refined.
pub const DWR_CUT_BAND_POLICY_VERSION: u16 = 2;

/// Strict fraction of absolute DWR mass that must lie on estimator-time
/// refinement-band leaves before the band receives one more level. The band
/// must also hold strictly more than its area share of the mass, so a band
/// that covers most of the domain cannot pass on area alone.
pub const DWR_CUT_BAND_MASS_GATE: f64 = 0.15;

/// Exhaustive reason emitted by the versioned DWR cut-band policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DwrBandDecision {
    /// The caller explicitly disabled DWR-driven band advancement.
    Disabled,
    /// Every supplied indicator had zero absolute mass.
    ZeroMass,
    /// Positive mass existed, but its estimator-time refinement-band share
    /// did not strictly exceed both [`DWR_CUT_BAND_MASS_GATE`] and the band's
    /// area share.
    GateNotMet,
    /// The mass gate passed, but the current band was already at the analysis
    /// grid's maximum level.
    LevelHeadroomExhausted,
    /// The gate passed with level headroom, authorizing exactly one level.
    Advanced,
}

impl DwrBandDecision {
    /// Whether this decision authorizes a one-level band advance.
    #[must_use]
    pub const fn is_advanced(self) -> bool {
        matches!(self, Self::Advanced)
    }
}

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

fn lattice_len(n: usize) -> usize {
    assert!(n >= 2, "density lattice needs at least 2 nodes per side");
    n.checked_mul(n)
        .expect("density lattice size overflows usize")
}

impl DensityDesign {
    /// A uniform-density start at `frac` solid fraction.
    ///
    /// # Panics
    /// If `n < 2`, `n * n` overflows, or `frac` is not finite and in
    /// `[0, 1]`.
    #[must_use]
    pub fn uniform(n: usize, frac: f64) -> DensityDesign {
        assert!(
            (0.0..=1.0).contains(&frac),
            "uniform density fraction must be finite and in [0, 1]"
        );
        let len = lattice_len(n);
        DensityDesign {
            n,
            rho: vec![frac; len],
        }
    }

    fn assert_shape(&self) {
        let expected = lattice_len(self.n);
        assert_eq!(
            self.rho.len(),
            expected,
            "density lattice length must equal n*n"
        );
    }

    fn node(&self, i: usize, j: usize) -> f64 {
        self.rho[j * self.n + i]
    }

    fn density_at_valid_shape(&self, x: f64, y: f64) -> f64 {
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

    /// Bilinear density at a point.
    #[must_use]
    pub fn density_at(&self, x: f64, y: f64) -> f64 {
        self.assert_shape();
        self.density_at_valid_shape(x, y)
    }

    /// Solid fraction (mean density).
    #[must_use]
    pub fn volume(&self) -> f64 {
        self.assert_shape();
        #[allow(clippy::cast_precision_loss)]
        {
            self.rho.iter().sum::<f64>() / self.rho.len() as f64
        }
    }

    /// Deterministic non-cryptographic witness over the exact lattice state.
    ///
    /// This is forensic correlation metadata, not a collision-resistant
    /// content address or a substitute for retaining the design artifact.
    #[must_use]
    pub fn state_witness(&self) -> u64 {
        self.assert_shape();
        let mut hash = 0xcbf2_9ce4_8422_2325u64;
        let mut mix_word = |word: u64| {
            for byte in word.to_le_bytes() {
                hash ^= u64::from(byte);
                hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
            }
        };
        mix_word(u64::try_from(self.n).expect("supported usize width fits u64"));
        for &density in &self.rho {
            mix_word(density.to_bits());
        }
        hash
    }

    /// Count void components (4-connected nodes with ρ ≤ ½) — the
    /// topology-evolution witness.
    #[must_use]
    pub fn void_components(&self) -> usize {
        self.assert_shape();
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
        self.assert_shape();
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
        0.5 - self.density_at_valid_shape(p[0], p[1])
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
            let v = self.density_at_valid_shape(x, y);
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
    /// Non-cryptographic witness for the exact design used by the solve and
    /// DWR estimator in this row.
    pub analysis_design_witness: u64,
    /// Thermal compliance `J = integral(f u)` on the analysis design.
    pub analysis_compliance: f64,
    /// Solid fraction of the analysis design.
    pub analysis_volume: f64,
    /// Void-component count of the analysis design.
    pub analysis_voids: usize,
    /// Non-cryptographic witness for the post-update design targeted by this
    /// row's refinement and consumed by the next solve.
    pub target_design_witness: u64,
    /// Solid fraction of the post-update target design.
    pub target_volume: f64,
    /// Void-component count of the post-update target design.
    pub target_voids: usize,
    /// Quadtree SPLITS this iteration: DWR-authorized next-design refinement
    /// plus post-update motion conformance.
    pub splits: usize,
    /// DWR mass-gate evidence from the exact design snapshot that produced the
    /// indicators, plus ONLY the policy-authorized structural splits applied
    /// to the post-update design. Motion-only splits are recorded separately.
    pub dwr_refinement: DwrBandRefinement,
    /// Splits needed after the density update to migrate an existing band to
    /// the new boundary when the DWR gate does not advance its level. This
    /// makes the final reported footprint conforming even when no next
    /// iteration exists to perform pre-solve repair.
    pub post_update_conformance_splits: usize,
    /// Leaf count after both policy-authorized and motion-conformance splits
    /// have targeted the post-update design.
    pub target_grid_leaf_count: usize,
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
    /// Witness matching the returned final design.
    pub final_design_witness: u64,
    /// Thermal compliance re-solved on the returned final design and final
    /// conformed grid. Unlike the last iteration's analysis compliance, this
    /// value describes `design` itself.
    pub final_compliance: f64,
    /// Leaf count of the final grid on which `final_compliance` was solved.
    pub final_grid_leaf_count: usize,
    /// Total quadtree splits across the run.
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

/// Evidence from one estimator-agnostic cut-band refinement decision.
///
/// The indicator source may be scalar heat or vector elasticity. This helper
/// only applies the marquee's shared planning policy; it does not claim that a
/// consumer can re-solve on the resulting graded grid.
#[derive(Debug, Clone, PartialEq)]
pub struct DwrBandRefinement {
    /// Versioned policy that produced this decision.
    pub policy_version: u16,
    /// Sum of absolute indicator mass on zero-straddling cells of the
    /// estimator's analysis snapshot (forensic; a subset of `band_mass`).
    pub cut_mass: f64,
    /// Sum of absolute indicator mass on refinement-band leaves (cut cells
    /// and their one-cell halo) of the analysis snapshot. The gate reads this.
    pub band_mass: f64,
    /// Summed area of the refinement-band leaves that carry indicators.
    pub band_area: f64,
    /// Summed area of every leaf that carries an indicator.
    pub total_area: f64,
    /// Sum of absolute indicator mass over every supplied cell.
    pub total_mass: f64,
    /// Band level before this decision.
    pub previous_level: u32,
    /// Band level after this decision.
    pub band_level: u32,
    /// Maximum level of the exact analysis-grid generation.
    pub analysis_max_level: u32,
    /// Exhaustive policy outcome; this distinguishes disabled, zero-mass,
    /// gate-failed, headroom-exhausted, and advancing decisions.
    pub decision: DwrBandDecision,
    /// Target-grid leaf count before applying the authorized refinement.
    pub leaves_before: usize,
    /// Target-grid leaf count after applying the authorized refinement.
    pub leaves_after: usize,
    /// Actual quadtree split count, including balance and halo splits, on the
    /// target design. The public one-shot helper uses one SDF for both the
    /// analysis snapshot and target; `run_marquee` deliberately targets the
    /// post-update design while retaining the pre-update mass evidence.
    pub splits: usize,
}

/// Validated, immutable authorization to advance the global cut-band level.
///
/// This private receipt is the time-level boundary between estimation and
/// design evolution: its masses are classified against the SAME geometry that
/// produced the DWR indicators. Applying the receipt may target the next
/// design, but cannot silently reclassify old indicators against that new SDF.
#[derive(Debug, Clone, PartialEq)]
struct DwrBandAdvanceReceipt {
    policy_version: u16,
    cut_mass: f64,
    band_mass: f64,
    band_area: f64,
    total_area: f64,
    total_mass: f64,
    previous_level: u32,
    decision: DwrBandDecision,
    analysis_max_level: u32,
    analysis_leaves: Vec<CellKey>,
}

/// True when the cell OR its one-cell halo is cut: fs-cutfem's ghost
/// penalty demands equal-level FACE NEIGHBORS of cut cells, so the
/// refinement band must include the halo, not just the straddling
/// cells (the CutBandNotUniform contract, learned the hard way twice).
fn halo_cut(sdf: &dyn CutSdf, lo: [f64; 2], hi: [f64; 2]) -> bool {
    sdf.enclose_halo(lo, hi).contains_zero()
}

/// Whether leaf `[lo, hi]` belongs to the DWR refinement band: the cell or
/// its one-cell halo straddles the zero level set. This is the exact
/// predicate the policy refines and, since policy version 2, the partition
/// its mass gate measures.
#[must_use]
pub fn in_dwr_refinement_band(sdf: &dyn CutSdf, lo: [f64; 2], hi: [f64; 2]) -> bool {
    halo_cut(sdf, lo, hi)
}

trait HaloEnclosure {
    fn enclose_halo(&self, lo: [f64; 2], hi: [f64; 2]) -> Interval;
}

impl<T: CutSdf + ?Sized> HaloEnclosure for T {
    fn enclose_halo(&self, lo: [f64; 2], hi: [f64; 2]) -> Interval {
        let (wx, wy) = (hi[0] - lo[0], hi[1] - lo[1]);
        let xlo = [(lo[0] - wx).max(0.0), (lo[1] - wy).max(0.0)];
        let xhi = [(hi[0] + wx).min(1.0), (hi[1] + wy).min(1.0)];
        self.enclose(xlo, xhi)
    }
}

/// Plan a fail-closed cut-band/halo refinement on a private grid clone.
fn plan_halo_refinement(
    grid: &Quadtree,
    sdf: &dyn CutSdf,
    target_level: u32,
) -> Result<(Quadtree, usize), fs_cutfem::CutFemError> {
    let max_level = grid.max_level();
    if target_level > max_level {
        return Err(fs_cutfem::CutFemError::InvalidFemInput {
            what: format!("cut-band target level {target_level} exceeds grid maximum {max_level}"),
        });
    }
    let mut planned = grid.clone();
    let leaves_before = planned.leaf_count();
    let enclosure_error = RefCell::new(None::<String>);
    planned.refine_where(target_level, &|lo, hi| {
        if enclosure_error.borrow().is_some() {
            return false;
        }
        let enclosure = sdf.enclose_halo(lo, hi);
        if enclosure.lo().is_finite() && enclosure.hi().is_finite() {
            enclosure.contains_zero()
        } else {
            *enclosure_error.borrow_mut() = Some(format!(
                "cut-band halo SDF enclosure for cell {lo:?}..{hi:?} is non-finite: [{}, {}]",
                enclosure.lo(),
                enclosure.hi()
            ));
            false
        }
    });
    if let Some(what) = enclosure_error.into_inner() {
        return Err(fs_cutfem::CutFemError::InvalidFemInput { what });
    }
    let leaves_after = planned.leaf_count();
    let added_leaves = leaves_after.checked_sub(leaves_before).ok_or_else(|| {
        fs_cutfem::CutFemError::InvalidFemInput {
            what: "cut-band refinement unexpectedly reduced the target leaf count".to_string(),
        }
    })?;
    if !added_leaves.is_multiple_of(3) {
        return Err(fs_cutfem::CutFemError::InvalidFemInput {
            what: format!(
                "cut-band refinement leaf delta {added_leaves} is not a quadtree split multiple"
            ),
        });
    }
    Ok((planned, added_leaves / 3))
}

fn conform_halo_to_level(
    grid: &mut Quadtree,
    sdf: &dyn CutSdf,
    target_level: u32,
) -> Result<usize, fs_cutfem::CutFemError> {
    let (planned, splits) = plan_halo_refinement(grid, sdf, target_level)?;
    *grid = planned;
    Ok(splits)
}

/// Validate and bind the mass decision to one exact analysis-grid generation.
fn classify_dwr_cut_band(
    grid: &Quadtree,
    sdf: &dyn CutSdf,
    indicators: &BTreeMap<CellKey, f64>,
    band_level: u32,
    enabled: bool,
) -> Result<DwrBandAdvanceReceipt, fs_cutfem::CutFemError> {
    let analysis_max_level = grid.max_level();
    if band_level > analysis_max_level {
        return Err(fs_cutfem::CutFemError::InvalidFemInput {
            what: format!("DWR band level {band_level} exceeds grid maximum {analysis_max_level}"),
        });
    }

    let mut total_mass = 0.0f64;
    let mut cut_mass = 0.0f64;
    let mut band_mass = 0.0f64;
    let mut band_area = 0.0f64;
    let mut total_area = 0.0f64;
    for (&cell, &eta) in indicators {
        if !grid.is_leaf(cell) {
            return Err(fs_cutfem::CutFemError::InvalidFemInput {
                what: format!("DWR indicator key {cell:?} is not a current grid leaf"),
            });
        }
        if !eta.is_finite() {
            return Err(fs_cutfem::CutFemError::InvalidFemInput {
                what: format!("DWR indicator for cell {cell:?} is non-finite: {eta}"),
            });
        }
        let (lo, hi) = grid.rect(cell);
        let enclosure = sdf.enclose(lo, hi);
        if !(enclosure.lo().is_finite() && enclosure.hi().is_finite()) {
            return Err(fs_cutfem::CutFemError::InvalidFemInput {
                what: format!(
                    "DWR SDF enclosure for cell {cell:?} is non-finite: [{}, {}]",
                    enclosure.lo(),
                    enclosure.hi()
                ),
            });
        }
        total_mass += eta.abs();
        if !total_mass.is_finite() {
            return Err(fs_cutfem::CutFemError::InvalidFemInput {
                what: "DWR total indicator mass is non-finite".to_string(),
            });
        }
        if enclosure.contains_zero() {
            cut_mass += eta.abs();
            if !cut_mass.is_finite() {
                return Err(fs_cutfem::CutFemError::InvalidFemInput {
                    what: "DWR cut-band indicator mass is non-finite".to_string(),
                });
            }
        }
        let halo = sdf.enclose_halo(lo, hi);
        if !(halo.lo().is_finite() && halo.hi().is_finite()) {
            return Err(fs_cutfem::CutFemError::InvalidFemInput {
                what: format!(
                    "DWR halo SDF enclosure for cell {cell:?} is non-finite: [{}, {}]",
                    halo.lo(),
                    halo.hi()
                ),
            });
        }
        let area = (hi[0] - lo[0]) * (hi[1] - lo[1]);
        total_area += area;
        if halo.contains_zero() {
            band_mass += eta.abs();
            band_area += area;
            if !band_mass.is_finite() {
                return Err(fs_cutfem::CutFemError::InvalidFemInput {
                    what: "DWR refinement-band indicator mass is non-finite".to_string(),
                });
            }
        }
    }
    // Share of mass versus share of area, cross-multiplied to stay exact.
    let concentrated = band_mass * total_area > band_area * total_mass;
    let decision = if !enabled {
        DwrBandDecision::Disabled
    } else if total_mass == 0.0 {
        DwrBandDecision::ZeroMass
    } else if band_mass <= DWR_CUT_BAND_MASS_GATE * total_mass || !concentrated {
        DwrBandDecision::GateNotMet
    } else if band_level >= analysis_max_level {
        DwrBandDecision::LevelHeadroomExhausted
    } else {
        DwrBandDecision::Advanced
    };

    Ok(DwrBandAdvanceReceipt {
        policy_version: DWR_CUT_BAND_POLICY_VERSION,
        cut_mass,
        band_mass,
        band_area,
        total_area,
        total_mass,
        previous_level: band_level,
        decision,
        analysis_max_level,
        analysis_leaves: grid.leaves().collect(),
    })
}

fn apply_dwr_band_receipt(
    grid: &mut Quadtree,
    target_sdf: &dyn CutSdf,
    receipt: DwrBandAdvanceReceipt,
    band_level: &mut u32,
) -> Result<DwrBandRefinement, fs_cutfem::CutFemError> {
    let DwrBandAdvanceReceipt {
        policy_version,
        cut_mass,
        band_mass,
        band_area,
        total_area,
        total_mass,
        previous_level,
        decision,
        analysis_max_level,
        analysis_leaves,
    } = receipt;
    if policy_version != DWR_CUT_BAND_POLICY_VERSION {
        return Err(fs_cutfem::CutFemError::InvalidFemInput {
            what: format!(
                "DWR band-advance receipt policy version {policy_version} does not match current version {DWR_CUT_BAND_POLICY_VERSION}"
            ),
        });
    }
    if grid.max_level() != analysis_max_level || !grid.leaves().eq(analysis_leaves) {
        return Err(fs_cutfem::CutFemError::InvalidFemInput {
            what: "DWR band-advance receipt does not match the current grid generation".to_string(),
        });
    }
    let current_level = *band_level;
    if current_level != previous_level {
        return Err(fs_cutfem::CutFemError::InvalidFemInput {
            what: format!(
                "DWR band-advance receipt was issued at level {previous_level} but current level is {current_level}"
            ),
        });
    }
    let leaves_before = grid.leaf_count();
    let (leaves_after, splits) = if decision.is_advanced() {
        let target_level = previous_level + 1;
        let (planned, splits) = plan_halo_refinement(grid, target_sdf, target_level)?;
        let leaves_after = planned.leaf_count();
        *grid = planned;
        *band_level = target_level;
        (leaves_after, splits)
    } else {
        (leaves_before, 0)
    };

    Ok(DwrBandRefinement {
        policy_version,
        cut_mass,
        band_mass,
        band_area,
        total_area,
        total_mass,
        previous_level,
        band_level: *band_level,
        analysis_max_level,
        decision,
        leaves_before,
        leaves_after,
        splits,
    })
}

/// Apply the marquee's shared DWR cut-band policy once.
///
/// The total marking mass is recomputed from `indicators`; callers cannot pass
/// an inconsistent denominator. A zero total, a disabled policy, or exhausted
/// level headroom is a deterministic no-op. Only
/// [`DwrBandDecision::Advanced`] advances the whole cut band and its one-cell
/// halo by exactly one level. The reported
/// `splits` is the actual structural count, not a requested split budget.
/// This one-shot API uses `sdf` for both mass classification and spatial
/// targeting. `run_marquee` uses the same private receipt machinery to bind
/// classification to the solved design while targeting the updated design.
///
/// # Errors
///
/// Returns [`fs_cutfem::CutFemError::InvalidFemInput`] without mutating the
/// grid or `band_level` when the level is out of range, an indicator is
/// non-finite, an indicator key is not a current leaf, an accumulated mass is
/// non-finite, or an SDF enclosure queried by the policy is non-finite.
pub fn refine_dwr_cut_band(
    grid: &mut Quadtree,
    sdf: &dyn CutSdf,
    indicators: &BTreeMap<CellKey, f64>,
    band_level: &mut u32,
    enabled: bool,
) -> Result<DwrBandRefinement, fs_cutfem::CutFemError> {
    let receipt = classify_dwr_cut_band(grid, sdf, indicators, *band_level, enabled)?;
    apply_dwr_band_receipt(grid, sdf, receipt, band_level)
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

fn evaluate_compliance(
    grid: &Quadtree,
    design: &DensityDesign,
) -> Result<f64, fs_cutfem::CutFemError> {
    let params = fem_params();
    let f = |_: f64, _: f64| 1.0;
    let g = |_: f64, _: f64| 0.0;
    let space = Space::build(grid, design, params)?;
    let solution = space.solve(&f, &g)?;
    let nodal = space.nodal_values(&solution.free, &g);
    goal_value(&space, &nodal, &GoalContext { weight: &f })
}

/// Run the marquee: the volume-to-point heat fixture (f = 1 body
/// heating, the design boundary cooled to 0) at a fixed solid
/// fraction. Interface-flux redistribution evolves the density;
/// the DWR cut-band mass gate enables at most one band-level advance
/// per iteration; the background grid is built ONCE and never rebuilt.
///
/// # Errors
/// CutFEM build/solve errors propagate.
#[allow(clippy::too_many_lines)] // one linear study loop: solve, grade, update, project, refine
pub fn run_marquee(
    mut design: DensityDesign,
    base_level: u32,
    max_level: u32,
    iters: usize,
    enable_band_refinement: bool,
) -> Result<MarqueeReport, fs_cutfem::CutFemError> {
    design.assert_shape();
    // THE GRID IS BUILT ONCE. Refinement = splits only; there is no
    // other construction site in this function (the zero-remeshing
    // property is structural, and the log proves it).
    let mut grid = Quadtree::with_room(base_level, max_level);
    let mut iterations = Vec::with_capacity(iters);
    let mut total_splits = 0usize;
    let mut band_level = base_level;
    let volume_constraint = design.volume();
    for iter in 0..iters {
        let t0 = std::time::Instant::now();
        let analysis_design_witness = design.state_witness();
        let analysis_volume = design.volume();
        let analysis_voids = design.void_components();
        // The previous iteration conformed this exact design after its update.
        // Do not silently repair it here: `Space::build` is the independent
        // fail-closed consumer that detects any broken cut-band invariant.
        let params = fem_params();
        let f = |_: f64, _: f64| 1.0;
        let g = |_: f64, _: f64| 0.0;
        let space = Space::build(&grid, &design, params)?;
        let sol = space.solve(&f, &g)?;
        let nodal = space.nodal_values(&sol.free, &g);
        let goal = GoalContext { weight: &f };
        let j = goal_value(&space, &nodal, &goal)?;
        // DWR per-leaf indicators for the compliance goal.
        let dwr = estimate(&grid, &design, params, &f, &g, &goal)?;
        // Authenticate the marking mass NOW, against the exact design
        // snapshot used by `estimate`. The density update below moves the
        // interface; classifying these old indicators against that new SDF
        // would mix time levels and can spuriously suppress the mass gate.
        let dwr_receipt = classify_dwr_cut_band(
            &grid,
            &design,
            &dwr.indicators,
            band_level,
            enable_band_refinement,
        )?;
        // --- Density update: interface-flux redistribution. ---------
        // Sample each lattice node's neighborhood; nodes NEAR the
        // interface get a signed move: high local flux² → carve (the
        // boundary wants to grow there, mye.1's shape derivative),
        // low flux² → fill; then project back to the volume target.
        let n = design.n;
        let u_at = |x: f64, y: f64| -> Result<f64, fs_cutfem::CutFemError> {
            // Bilinear through the containing leaf's corner nodes via
            // the canonical fail-closed sampler (ay40): missing or
            // non-finite active evidence refuses instead of reading as
            // a plausible zero (probing the raw fine lattice missed
            // the sparse keys and froze the whole update in the first
            // draft; zero-filling then hid exactly that class of bug).
            // A certified-Outside leaf reads the homogeneous Dirichlet
            // exterior u = 0 explicitly. The rim clamp keeps probes
            // inside the half-open background box.
            let p = [x.clamp(1e-9, 1.0 - 1e-9), y.clamp(1e-9, 1.0 - 1e-9)];
            match space.sample_scalar(&nodal, p)? {
                ScalarSample::Active(v) => Ok(v),
                ScalarSample::CertifiedOutside => Ok(0.0),
            }
        };
        let flux_at = |x: f64, y: f64| -> Result<f64, fs_cutfem::CutFemError> {
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
            let u = u_at(x - depth * gph[0] / norm, y - depth * gph[1] / norm)?;
            Ok((u / h).powi(2))
        };
        let mut moves = vec![0.0f64; n * n];
        let mut flux_sum = 0.0f64;
        let mut flux_cnt = 0usize;
        #[allow(clippy::cast_precision_loss)]
        let lattice_scale = (n - 1) as f64;
        for (k, slot) in moves.iter_mut().enumerate() {
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
                let fl = flux_at(x, y)?;
                *slot = fl.max(1e-12);
                flux_sum += fl;
                flux_cnt += 1;
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let flux_mean = flux_sum / flux_cnt.max(1) as f64;
        let step = 0.25;
        for (rho, &m_k) in design.rho.iter_mut().zip(&moves) {
            if m_k > 0.0 {
                // carve where flux is above the mean, fill below.
                let rel = (m_k - flux_mean) / flux_mean.max(1e-12);
                *rho = (*rho - step * rel.clamp(-1.0, 1.0)).clamp(0.02, 0.98);
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
                if vol > volume_constraint {
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
        // cells individually and the solver refused, correctly). The
        // estimator-agnostic helper is also the integration surface for
        // vector compliance indicators; it applies planning policy only.
        // The decision is the pre-update receipt above, while its structural
        // target is the post-update design that the NEXT solve will consume.
        // Thus both the DWR evidence and the final moving-band footprint have
        // honest, explicit time levels.
        let refinement = apply_dwr_band_receipt(&mut grid, &design, dwr_receipt, &mut band_level)?;
        // A non-advancing receipt still leaves an EXISTING fine band to migrate
        // after the interface moves. Do that now, including on the final
        // iteration, instead of relying on a next-iteration preflight that may
        // never run. An advancing receipt already performed this exact target
        // refinement, so only non-advancing outcomes need the motion pass.
        let post_splits = if band_level > base_level && !refinement.decision.is_advanced() {
            conform_halo_to_level(&mut grid, &design, band_level)?
        } else {
            0
        };
        let splits = refinement.splits + post_splits;
        total_splits += splits;
        let target_grid_leaf_count = grid.leaf_count();
        let target_design_witness = design.state_witness();
        let target_volume = design.volume();
        let target_voids = design.void_components();
        #[allow(clippy::cast_precision_loss)]
        let wall_ms = t0.elapsed().as_secs_f64() * 1e3;
        iterations.push(MarqueeIter {
            iter,
            analysis_design_witness,
            analysis_compliance: j,
            analysis_volume,
            analysis_voids,
            target_design_witness,
            target_volume,
            target_voids,
            splits,
            dwr_refinement: refinement,
            post_update_conformance_splits: post_splits,
            target_grid_leaf_count,
            rebuilds: 0, // structural: there is no rebuild path
            wall_ms,
        });
    }
    let final_design_witness = design.state_witness();
    let final_compliance = evaluate_compliance(&grid, &design)?;
    let final_grid_leaf_count = grid.leaf_count();
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
        final_design_witness,
        final_compliance,
        final_grid_leaf_count,
        total_splits,
        refined_boundary_leaves,
        refined_off_boundary_leaves,
        total_rebuilds: 0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_cutfem::sdf::HalfPlane;

    struct NonFiniteTarget;

    impl CutSdf for NonFiniteTarget {
        fn value(&self, _point: [f64; 2]) -> f64 {
            0.0
        }

        fn gradient(&self, _point: [f64; 2]) -> [f64; 2] {
            [1.0, 0.0]
        }

        fn enclose(&self, _lo: [f64; 2], _hi: [f64; 2]) -> Interval {
            Interval::WHOLE
        }
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "the regression keeps the analysis/target discriminator, topology oracle, and stale-receipt atomicity in one narrative"
    )]
    fn dwr_receipt_keeps_analysis_gate_when_target_boundary_moves() {
        let analysis_sdf = HalfPlane {
            normal: [1.0, 0.0],
            offset: 0.125,
        };
        let target_sdf = HalfPlane {
            normal: [1.0, 0.0],
            offset: 0.875,
        };
        let make_grid = || Quadtree::with_room(2, 3);
        let analysis_grid = make_grid();
        let indicators: BTreeMap<_, _> = analysis_grid
            .leaves()
            .map(|cell| {
                let (lo, hi) = analysis_grid.rect(cell);
                let eta = if analysis_sdf.enclose(lo, hi).contains_zero() {
                    1.0
                } else {
                    0.001
                };
                (cell, eta)
            })
            .collect();

        let receipt = classify_dwr_cut_band(&analysis_grid, &analysis_sdf, &indicators, 2, true)
            .expect("analysis snapshot is valid");
        assert_eq!(
            receipt.decision,
            DwrBandDecision::Advanced,
            "analysis boundary carries the DWR mass"
        );
        let analysis_cut_mass = receipt.cut_mass;
        let analysis_total_mass = receipt.total_mass;

        let mut transitioned_grid = make_grid();
        let mut transitioned_level = 2;
        let transitioned = apply_dwr_band_receipt(
            &mut transitioned_grid,
            &target_sdf,
            receipt,
            &mut transitioned_level,
        )
        .expect("exact-grid receipt may target the next design");
        assert_eq!(transitioned.decision, DwrBandDecision::Advanced);
        assert_eq!(transitioned.cut_mass.to_bits(), analysis_cut_mass.to_bits());
        assert_eq!(
            transitioned.total_mass.to_bits(),
            analysis_total_mass.to_bits()
        );
        assert!(transitioned.splits > 0);
        assert_eq!(transitioned_level, 3);

        let mut expected_target_grid = make_grid();
        expected_target_grid.refine_where(3, &|lo, hi| halo_cut(&target_sdf, lo, hi));
        assert_eq!(
            transitioned_grid.leaves().collect::<Vec<_>>(),
            expected_target_grid.leaves().collect::<Vec<_>>(),
            "the receipt authorizes the level while the target SDF selects the spatial band"
        );

        let mut incorrectly_reclassified_grid = make_grid();
        let mut incorrectly_reclassified_level = 2;
        let incorrectly_reclassified = refine_dwr_cut_band(
            &mut incorrectly_reclassified_grid,
            &target_sdf,
            &indicators,
            &mut incorrectly_reclassified_level,
            true,
        )
        .expect("the moved target is itself finite and valid");
        assert!(
            incorrectly_reclassified.cut_mass
                <= DWR_CUT_BAND_MASS_GATE * incorrectly_reclassified.total_mass,
            "the fixture must discriminate the two geometry time levels"
        );
        assert_eq!(
            incorrectly_reclassified.decision,
            DwrBandDecision::GateNotMet
        );
        assert_eq!(incorrectly_reclassified.splits, 0);
        assert_eq!(incorrectly_reclassified_level, 2);
        assert_eq!(
            incorrectly_reclassified_grid.leaves().collect::<Vec<_>>(),
            analysis_grid.leaves().collect::<Vec<_>>()
        );

        let stale_receipt =
            classify_dwr_cut_band(&analysis_grid, &analysis_sdf, &indicators, 2, true)
                .expect("second analysis receipt");
        let mut stale_grid = analysis_grid.clone();
        let first_leaf = stale_grid
            .leaves()
            .next()
            .expect("uniform grid is nonempty");
        stale_grid.split(first_leaf);
        let stale_before: Vec<_> = stale_grid.leaves().collect();
        let mut stale_level = 2;
        let error = apply_dwr_band_receipt(
            &mut stale_grid,
            &target_sdf,
            stale_receipt,
            &mut stale_level,
        )
        .expect_err("receipt cannot cross a grid-generation boundary");
        assert!(matches!(
            error,
            fs_cutfem::CutFemError::InvalidFemInput { .. }
        ));
        assert_eq!(stale_grid.leaves().collect::<Vec<_>>(), stale_before);
        assert_eq!(stale_level, 2);

        let stale_level_receipt =
            classify_dwr_cut_band(&analysis_grid, &analysis_sdf, &indicators, 2, true)
                .expect("third analysis receipt");
        let mut wrong_level_grid = analysis_grid.clone();
        let wrong_level_before: Vec<_> = wrong_level_grid.leaves().collect();
        let mut wrong_level = 3;
        apply_dwr_band_receipt(
            &mut wrong_level_grid,
            &target_sdf,
            stale_level_receipt,
            &mut wrong_level,
        )
        .expect_err("receipt cannot cross a band-level boundary");
        assert_eq!(
            wrong_level_grid.leaves().collect::<Vec<_>>(),
            wrong_level_before
        );
        assert_eq!(wrong_level, 3);

        let mut stale_policy_receipt =
            classify_dwr_cut_band(&analysis_grid, &analysis_sdf, &indicators, 2, true)
                .expect("fourth analysis receipt");
        stale_policy_receipt.policy_version += 1;
        let mut wrong_policy_grid = analysis_grid.clone();
        let wrong_policy_before: Vec<_> = wrong_policy_grid.leaves().collect();
        let mut wrong_policy_level = 2;
        apply_dwr_band_receipt(
            &mut wrong_policy_grid,
            &target_sdf,
            stale_policy_receipt,
            &mut wrong_policy_level,
        )
        .expect_err("receipt cannot cross a policy-version boundary");
        assert_eq!(
            wrong_policy_grid.leaves().collect::<Vec<_>>(),
            wrong_policy_before
        );
        assert_eq!(wrong_policy_level, 2);

        let nonfinite_receipt =
            classify_dwr_cut_band(&analysis_grid, &analysis_sdf, &indicators, 2, true)
                .expect("fifth analysis receipt");
        let mut nonfinite_grid = analysis_grid.clone();
        let nonfinite_before: Vec<_> = nonfinite_grid.leaves().collect();
        let mut nonfinite_level = 2;
        apply_dwr_band_receipt(
            &mut nonfinite_grid,
            &NonFiniteTarget,
            nonfinite_receipt,
            &mut nonfinite_level,
        )
        .expect_err("non-finite target refuses before mutation");
        assert_eq!(
            nonfinite_grid.leaves().collect::<Vec<_>>(),
            nonfinite_before
        );
        assert_eq!(nonfinite_level, 2);
    }

    #[test]
    fn existing_band_conforms_to_moved_target_without_an_advance() {
        let analysis_sdf = HalfPlane {
            normal: [1.0, 0.0],
            offset: 0.125,
        };
        let target_sdf = HalfPlane {
            normal: [1.0, 0.0],
            offset: 0.875,
        };
        let mut grid = Quadtree::with_room(2, 3);
        let initial_splits =
            conform_halo_to_level(&mut grid, &analysis_sdf, 3).expect("initial authorized band");
        assert!(initial_splits > 0);

        let indicators: BTreeMap<_, _> = grid
            .leaves()
            .map(|cell| {
                let (lo, hi) = grid.rect(cell);
                let eta = if analysis_sdf.enclose(lo, hi).contains_zero() {
                    1.0
                } else {
                    1.0 / 64.0
                };
                (cell, eta)
            })
            .collect();
        let receipt = classify_dwr_cut_band(&grid, &analysis_sdf, &indicators, 3, true)
            .expect("max-level analysis receipt");
        assert_eq!(receipt.decision, DwrBandDecision::LevelHeadroomExhausted);
        let mut band_level = 3;
        let refinement = apply_dwr_band_receipt(&mut grid, &target_sdf, receipt, &mut band_level)
            .expect("non-advancing receipt is a structural no-op");
        assert_eq!(refinement.splits, 0);

        let before_motion = grid.clone();
        let motion_splits = if refinement.decision.is_advanced() {
            0
        } else {
            conform_halo_to_level(&mut grid, &target_sdf, band_level)
                .expect("moved target conformance")
        };
        assert!(
            motion_splits > 0,
            "disjoint target band requires new leaves"
        );
        assert_eq!(
            grid.leaf_count() - before_motion.leaf_count(),
            3 * motion_splits
        );

        let mut expected = before_motion;
        expected.refine_where(3, &|lo, hi| halo_cut(&target_sdf, lo, hi));
        assert_eq!(
            grid.leaves().collect::<Vec<_>>(),
            expected.leaves().collect::<Vec<_>>(),
            "motion conformance equals the independently derived target halo"
        );
        Space::build(&grid, &target_sdf, fem_params())
            .expect("the next scalar solve accepts the moved, conformed band");
    }
}
