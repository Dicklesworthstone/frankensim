//! Nested plane-wave/Chebyshev bases for the outgoing Helmholtz channels.
//!
//! Messner, Schanz & Darve, JCP 231 (2012), 1175-1196,
//! doi:10.1016/j.jcp.2011.09.027. The centered basis and interlevel
//! re-expansion are also described in Borm & Melenk, arXiv:1510.07189,
//! section 2.2. No panel quadrature lives here.
use super::{ComplexBudget, FmmC64, Operator, finite, kernels, validate_output};
use crate::{
    Cell, CellKey, ComplexKernel, Fmm, FmmError, HelmholtzAdjointDouble, HelmholtzDouble,
    HelmholtzHypersingular, Laplace3d, MAX_BB_HELMHOLTZ_KD, MAX_CELL_PAIR_SCANS,
    MAX_COMPLEX_COEFFICIENT_SLOTS, MAX_INTERPOLATION_ORDER, MAX_NEAR_FIELD_PAIR_WORK, MAX_POINTS,
    MAX_TRANSLATION_WORK, MIN_INTERPOLATION_ORDER, point_leaves,
};
use fs_math::{c64::C64, det};
use std::collections::{BTreeMap, BTreeSet};

/// Only active cones allocate coefficients; this caps the angular grid itself.
pub const MAX_DIRECTIONAL_CONES: usize = 16_777_216;
/// Explicit bound on a sparse directional tree's geometric refinement.
pub const MAX_DIRECTIONAL_DEPTH: u32 = 10;
/// A coupling cannot have rank above its uncompressed tensor-grid dimension.
pub const MAX_DIRECTIONAL_RANK: usize = MAX_INTERPOLATION_ORDER.pow(3);
/// Largest declarable coefficient budget, in C64 units: its f64 count and
/// byte size stay representable.
pub const MAX_DECLARED_COEFFICIENT_SLOTS: usize = usize::MAX / (2 * core::mem::size_of::<C64>());
const ZERO_CONE: Cone = (6, 0, 0);
type Cone = (usize, usize, usize); // cube face, first bin, second bin
type Expansion = (CellKey, Cone);

/// Explicit directional discretization, accuracy design and resource choices.
#[derive(Clone, Debug, PartialEq)]
pub struct DirectionalOptions {
    /// Requested relative field accuracy on the independently qualified domain.
    pub tolerance: f64,
    /// Maximum size of an active level's angular grid, not a dense allocation.
    pub max_cones_per_level: usize,
    /// Chebyshev nodes per coordinate direction (2..=12).
    pub order: usize,
    /// Low-frequency bases have direction zero below this cell k*side.
    pub low_freq_threshold: f64,
    /// Angular resolution multiplier; it never silently changes.
    pub cone_resolution: usize,
    /// Occupancy target; frequency refinement can create smaller leaves.
    pub leaf_capacity: usize,
    /// Relative Frobenius residual target for each demodulated M2L block.
    pub compression_tolerance: f64,
    /// Coupling rank budget; only needed when a shared block is compressed.
    pub max_rank: usize,
    /// Simultaneous pass coefficients, retained coupling factors and setup
    /// scratch, in C64 units. Defaults to [`MAX_COMPLEX_COEFFICIENT_SLOTS`].
    /// A request whose root kD is at most `low_freq_threshold` delegates to
    /// the unchanged W5a owner under `w5a_budget` instead.
    pub max_coefficient_slots: usize,
    /// Counted coupling-setup products. Construction polls no cancellation,
    /// so this also bounds one uncancellable setup. Defaults to
    /// [`MAX_TRANSLATION_WORK`]; the same W5a delegation applies.
    pub max_setup_work: usize,
    /// W5a ceilings for a request whose root kD is at most
    /// `low_freq_threshold`, which delegates whole to W5a at this order and
    /// leaf capacity. Defaults to [`ComplexBudget::V1`].
    pub w5a_budget: ComplexBudget,
    /// Counted near-field point pairs. This route caches none, so every pass
    /// evaluates them directly. Defaults to [`MAX_NEAR_FIELD_PAIR_WORK`]; a
    /// plan with more pairs is refused before any coupling setup.
    pub max_near_pairs: usize,
}
impl DirectionalOptions {
    /// Order selected by the explicit 4^-p design indicator, not a universal
    /// relative-error certificate for arbitrary cancelling source fields.
    pub fn for_tolerance(tolerance: f64) -> Result<Self, FmmError> {
        let order = required_order(tolerance)?;
        Ok(Self {
            tolerance,
            order,
            max_cones_per_level: MAX_DIRECTIONAL_CONES,
            low_freq_threshold: 4.0,
            cone_resolution: 1,
            leaf_capacity: 64,
            compression_tolerance: tolerance * 0.01,
            max_rank: MAX_DIRECTIONAL_RANK,
            max_coefficient_slots: MAX_COMPLEX_COEFFICIENT_SLOTS,
            max_setup_work: MAX_TRANSLATION_WORK,
            w5a_budget: ComplexBudget::V1,
            max_near_pairs: MAX_NEAR_FIELD_PAIR_WORK,
        })
    }
    /// Validate parameters without building geometry or evaluating a kernel.
    pub fn validate(&self) -> Result<(), FmmError> {
        let required = required_order(self.tolerance)?;
        if !self.compression_tolerance.is_finite()
            || self.compression_tolerance <= 0.0
            || self.compression_tolerance >= 1.0
        {
            return Err(FmmError::ComplexRequest {
                reason: "M2L compression tolerance must be finite and in (0,1)",
            });
        }
        budget(
            "directional coupling rank",
            self.max_rank,
            MAX_DIRECTIONAL_RANK,
        )?;
        if !(MIN_INTERPOLATION_ORDER..=MAX_INTERPOLATION_ORDER).contains(&self.order) {
            return Err(FmmError::InvalidOrder {
                order: self.order,
                min: MIN_INTERPOLATION_ORDER,
                max: MAX_INTERPOLATION_ORDER,
            });
        }
        if self.order < required {
            return Err(FmmError::DirectionalAccuracy {
                tolerance: self.tolerance,
                order: self.order,
                required_order: required,
            });
        }
        if !self.low_freq_threshold.is_finite()
            || self.low_freq_threshold <= 0.0
            || self.low_freq_threshold > MAX_BB_HELMHOLTZ_KD
            || !(1..=4).contains(&self.cone_resolution)
        {
            return Err(FmmError::ComplexRequest {
                reason: "directional low-frequency threshold must be in (0,10], cone resolution in 1..=4",
            });
        }
        if self.leaf_capacity == 0 || self.leaf_capacity > MAX_POINTS {
            return Err(FmmError::InvalidLeafCapacity {
                capacity: self.leaf_capacity,
                max: MAX_POINTS,
            });
        }
        budget(
            "directional cone grid",
            self.max_cones_per_level,
            MAX_DIRECTIONAL_CONES,
        )?;
        budget(
            "directional declared coefficient-slot",
            self.max_coefficient_slots,
            MAX_DECLARED_COEFFICIENT_SLOTS,
        )?;
        if self.max_coefficient_slots == 0 || self.max_setup_work == 0 || self.max_near_pairs == 0 {
            return Err(FmmError::ComplexRequest {
                reason: "directional coefficient-slot, setup-work and near-pair budgets must be positive",
            });
        }
        Ok(())
    }
}
fn required_order(tol: f64) -> Result<usize, FmmError> {
    if !tol.is_finite() || tol <= 0.0 || tol >= 1.0 {
        return Err(FmmError::ComplexRequest {
            reason: "directional tolerance must be finite and in (0,1)",
        });
    }
    let required = ((-det::ln(tol) / det::ln(4.0)).ceil() as usize).max(MIN_INTERPOLATION_ORDER);
    if required > MAX_INTERPOLATION_ORDER {
        return Err(FmmError::DirectionalAccuracy {
            tolerance: tol,
            order: MAX_INTERPOLATION_ORDER,
            required_order: required,
        });
    }
    Ok(required)
}
fn budget(operation: &'static str, requested: usize, max: usize) -> Result<(), FmmError> {
    if requested > max {
        Err(FmmError::WorkEnvelopeExceeded {
            operation,
            requested,
            max,
        })
    } else {
        Ok(())
    }
}

/// Counted directional work/storage. Timings and RSS are measured by callers.
#[derive(Clone, Debug)]
pub struct DirectionalStats {
    /// Nominal angular grid sizes on used levels, zero on unused levels.
    pub cones_per_level: Vec<usize>,
    /// Active source cell/cone bases.
    pub source_expansions: usize,
    /// Active target cell/cone bases.
    pub target_expansions: usize,
    /// Accepted ordered far-cell pairs.
    pub m2l_pairs: usize,
    /// Accepted pairs with a nonzero plane-wave direction.
    pub high_frequency_pairs: usize,
    /// Exact leaf point pairs, including diagonal slots subsequently skipped.
    pub near_pairs: usize,
    /// Dual-tree admission tests.
    pub pair_scans: usize,
    /// Complex coupling products per matvec.
    pub m2l_work: usize,
    /// Interlevel complex weighted transfers per matvec.
    pub transfer_work: usize,
    /// Source/target leaf interpolation products per matvec.
    pub leaf_work: usize,
    /// Reserved pass coefficients and retained coupling factors in C64 units.
    pub coefficient_slots: usize,
    /// Retained directional numeric/index payload; excludes base/allocator overhead.
    pub memory_bytes: usize,
    /// Tensor order used in this plan.
    pub order: usize,
    /// Requested relative accuracy, separate from any measured oracle error.
    pub requested_tolerance: f64,
    /// Parameter-design indicator 4^-p. Actual error is measured by the
    /// independent oracle ladders; this is not an enclosure.
    pub interpolation_indicator: f64,
    /// Number of unique coupling keys (compressed or uncompressed).
    pub coupling_blocks: usize,
    /// Largest retained compressed rank; zero if no shared coupling was needed.
    pub maximum_rank: usize,
    /// Largest observed relative grid-block compression residual.
    pub compression_residual: f64,
    /// Counted coupling construction and residual-update products.
    pub setup_work: usize,
}
impl DirectionalStats {
    /// Counted matvec products, near pairs and tree tests; setup is separate.
    #[must_use]
    pub fn work(&self) -> usize {
        self.m2l_work + self.transfer_work + self.leaf_work + self.near_pairs + self.pair_scans
    }
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"directional\":true,\"cones_per_level\":{:?},\"source_expansions\":{},\"target_expansions\":{},\"m2l_pairs\":{},\"high_frequency_pairs\":{},\"near_pairs\":{},\"pair_scans\":{},\"m2l_work\":{},\"transfer_work\":{},\"leaf_work\":{},\"work\":{},\"complex_coefficient_slots\":{},\"memory_bytes\":{},\"order\":{},\"requested_tolerance\":{},\"interpolation_indicator\":{},\"coupling_blocks\":{},\"maximum_rank\":{},\"compression_residual\":{},\"setup_work\":{}}}",
            self.cones_per_level,
            self.source_expansions,
            self.target_expansions,
            self.m2l_pairs,
            self.high_frequency_pairs,
            self.near_pairs,
            self.pair_scans,
            self.m2l_work,
            self.transfer_work,
            self.leaf_work,
            self.work(),
            self.coefficient_slots,
            self.memory_bytes,
            self.order,
            self.requested_tolerance,
            self.interpolation_indicator,
            self.coupling_blocks,
            self.maximum_rank,
            self.compression_residual,
            self.setup_work
        )
    }
}

pub(super) struct DirectionalPlan {
    pub(super) stats: DirectionalStats,
    options: DirectionalOptions,
    bins: Vec<usize>,
    pairs: Vec<(CellKey, CellKey, Cone)>,
    near: Vec<(CellKey, CellKey)>,
    sources: BTreeSet<Expansion>,
    targets: BTreeSet<Expansion>,
    k: f64,
    coupling_keys: Vec<CouplingKey>,
    couplings: BTreeMap<CouplingKey, Coupling>,
}
enum Coupling {
    TensorGrid,
    Compressed(Vec<LowRank>),
}
type CouplingKey = (CellKey, CellKey, Cone);
struct LowRank {
    u: Vec<Vec<C64>>,
    v: Vec<Vec<C64>>,
    residual: f64,
}
impl LowRank {
    fn rank(&self) -> usize {
        self.u.len()
    }
    fn slots(&self) -> usize {
        self.u.iter().map(Vec::len).sum::<usize>() + self.v.iter().map(Vec::len).sum::<usize>()
    }
}

impl<'k> FmmC64<'k> {
    /// Construct a directional scalar/gradient operator on a source-target cloud.
    pub fn new_directional(
        kernel: &'k dyn ComplexKernel,
        points: &[[f64; 3]],
        options: DirectionalOptions,
    ) -> Result<Self, FmmError> {
        Self::build_directional(Operator::Scalar(kernel), points.to_vec(), None, options)
    }
    /// Source-normal derivative with directional interpolation of its channels.
    pub fn new_double_directional(
        kernel: HelmholtzDouble,
        points: Vec<[f64; 3]>,
        normals: Vec<[f64; 3]>,
        options: DirectionalOptions,
    ) -> Result<Self, FmmError> {
        Self::build_directional(Operator::Double(kernel), points, Some(normals), options)
    }
    /// Target-normal derivative with directional interpolation of its channels.
    pub fn new_adjoint_double_directional(
        kernel: HelmholtzAdjointDouble,
        points: Vec<[f64; 3]>,
        normals: Vec<[f64; 3]>,
        options: DirectionalOptions,
    ) -> Result<Self, FmmError> {
        Self::build_directional(Operator::Adjoint(kernel), points, Some(normals), options)
    }
    /// Off-diagonal mixed derivative; self/finite-part terms remain with BEM.
    pub fn new_hypersingular_directional(
        kernel: HelmholtzHypersingular,
        points: Vec<[f64; 3]>,
        normals: Vec<[f64; 3]>,
        options: DirectionalOptions,
    ) -> Result<Self, FmmError> {
        Self::build_directional(Operator::Mixed(kernel), points, Some(normals), options)
    }
    fn build_directional(
        operator: Operator<'k>,
        points: Vec<[f64; 3]>,
        normals: Option<Vec<[f64; 3]>>,
        options: DirectionalOptions,
    ) -> Result<Self, FmmError> {
        options.validate()?;
        // This validates coordinates/count without imposing the old all-cell-
        // pairs scan or non-directional translation bound on a different pass.
        // Actual directional scans/work/storage are admitted separately below.
        Fmm::validate_request(&points, 2, points.len())?;
        let k = operator.wavenumber().ok_or(FmmError::ComplexRequest {
            reason: "directional kernel must declare its Helmholtz wavenumber",
        })?;
        if !k.is_finite() || k < 0.0 {
            return Err(FmmError::ComplexRequest {
                reason: "directional wavenumber must be finite and nonnegative",
            });
        }
        if let Some(ns) = &normals {
            if ns.len() != points.len()
                || ns.iter().any(|n| {
                    n.iter().any(|v| !v.is_finite()) || (kernels::dot(*n, *n) - 1.0).abs() > 1e-12
                })
            {
                return Err(FmmError::ComplexRequest {
                    reason: "directional normal field requires one finite unit normal per point",
                });
            }
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in &points {
            for a in 0..3 {
                lo[a] = lo[a].min(p[a]);
                hi[a] = hi[a].max(p[a]);
            }
        }
        let side = (hi[0] - lo[0])
            .max(hi[1] - lo[1])
            .max(hi[2] - lo[2])
            .max(1e-12)
            * (1.0 + 1e-9);
        let kd = k * side;
        if !kd.is_finite() {
            return Err(FmmError::ComplexRequest {
                reason: "nonfinite directional root kD",
            });
        }
        // Entire low-frequency requests use the W5a code without modifying its
        // tree, pair partition, interpolation, reductions or arithmetic.
        if kd <= options.low_freq_threshold {
            return Self::build(
                operator,
                points,
                normals,
                options.order,
                options.leaf_capacity,
                options.w5a_budget,
            );
        }
        // d_axis G = (x_axis-y_axis) G'(r)/r. Interpolate only the radial
        // remainder and form the affine factor through source/target moments.
        // Directly interpolating d_axis G loses a degree and can give a large
        // relative error for transverse components near their angular zero.
        // This conversion follows the unchanged whole-request W5a branch.
        let operator = match operator {
            Operator::Scalar(kernel) => kernel
                .helmholtz_gradient()
                .map_or(Operator::Scalar(kernel), Operator::Gradient),
            other => other,
        };
        let mut depth = 1u32;
        while 8usize.saturating_pow(depth) < points.len().div_ceil(options.leaf_capacity)
            || kd / f64::from(1u32 << depth) > options.low_freq_threshold
        {
            depth += 1;
            budget(
                "directional depth",
                depth as usize,
                MAX_DIRECTIONAL_DEPTH as usize,
            )?;
        }
        let n_side = 1u32 << depth;
        let mut cells = point_leaves(&points, lo, side / f64::from(n_side), depth, n_side);
        let leaves: Vec<_> = cells.keys().copied().collect();
        for mut key in leaves {
            while key.0 > 0 {
                key = (key.0 - 1, key.1 >> 1, key.2 >> 1, key.3 >> 1);
                cells.entry(key).or_insert_with(|| Cell {
                    points: Vec::new(),
                    leaf: false,
                });
            }
        }
        let geometry = Fmm {
            kernel: &Laplace3d,
            points,
            order: options.order,
            max_level: depth,
            cells,
            lo,
            side,
        };
        let (mut nodes, mut barycentric) = (Vec::new(), Vec::new());
        for j in 0..options.order {
            let angle = (2 * j + 1) as f64 * core::f64::consts::PI / (2 * options.order) as f64;
            let (s, c) = det::sin_cos(angle);
            nodes.push(-c);
            barycentric.push(if j % 2 == 0 { s } else { -s });
        }
        let mut op = Self {
            geometry,
            operator,
            normals,
            nodes,
            barycentric,
            interactions: Vec::new(),
            translations: Vec::new(),
            near: None,
            coefficient_slots: 0,
            directional: None,
            workers: 1,
            m2l_rows: None,
        };
        let plan = DirectionalPlan::build(&op, options, k)?;
        op.coefficient_slots = plan.stats.coefficient_slots;
        op.directional = Some(plan);
        Ok(op)
    }
    #[must_use]
    /// None means the whole request delegated to the unchanged W5a owner.
    pub fn directional_stats(&self) -> Option<&DirectionalStats> {
        self.directional.as_ref().map(|p| &p.stats)
    }
    /// Independent sampled direct oracle, excluding self by index. At most
    /// 1000 deterministic caller-selected targets; no large N*N allocation.
    pub fn direct_at_targets(
        &self,
        charges: &[C64],
        targets: &[usize],
    ) -> Result<Vec<C64>, FmmError> {
        self.validate_charges(charges)?;
        budget("sampled-direct target", targets.len(), 1000)?;
        budget(
            "sampled-direct point pair",
            targets.len().saturating_mul(charges.len()),
            MAX_NEAR_FIELD_PAIR_WORK,
        )?;
        if targets.iter().any(|&i| i >= charges.len()) {
            return Err(FmmError::ComplexRequest {
                reason: "sampled direct target index outside point cloud",
            });
        }
        let out = targets
            .iter()
            .map(|&t| {
                charges
                    .iter()
                    .enumerate()
                    .filter(|(s, _)| *s != t)
                    .fold(C64::ZERO, |sum, (s, &q)| sum + self.pair(t, s) * q)
            })
            .collect();
        validate_output(out)
    }
}

fn children(op: &FmmC64<'_>, key: CellKey) -> Vec<CellKey> {
    let mut out = Vec::new();
    for di in 0..2 {
        for dj in 0..2 {
            for dk in 0..2 {
                let child = (key.0 + 1, 2 * key.1 + di, 2 * key.2 + dj, 2 * key.3 + dk);
                if op.geometry.cells.contains_key(&child) {
                    out.push(child);
                }
            }
        }
    }
    out
}
fn center(op: &FmmC64<'_>, key: CellKey) -> [f64; 3] {
    let (lo, h) = op.geometry.cell_box(key);
    lo.map(|v| v + 0.5 * h)
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn wave(k: f64, u: [f64; 3], d: [f64; 3]) -> C64 {
    let phase = k * kernels::dot(u, d);
    if phase == 0.0 {
        C64::ONE
    } else {
        let (s, c) = det::sin_cos(phase);
        C64::new(c, s)
    }
}
fn modulate(value: C64, phase: C64) -> C64 {
    // In zero-direction bases this is exactly W5a's original arithmetic.
    if phase == C64::ONE {
        value
    } else {
        value * phase
    }
}
fn cone(u: [f64; 3], q: usize) -> Cone {
    if q == 0 {
        return ZERO_CONE;
    }
    let mut axis = 0;
    for a in 1..3 {
        if u[a].abs() > u[axis].abs() {
            axis = a;
        }
    }
    let bin = |a: usize| {
        (((u[a] / u[axis].abs() + 1.0) * 0.5 * q as f64)
            .floor()
            .max(0.0) as usize)
            .min(q - 1)
    };
    let a = (axis + 1) % 3;
    let b = (axis + 2) % 3;
    (2 * axis + usize::from(u[axis] < 0.0), bin(a), bin(b))
}
fn direction(c: Cone, q: usize) -> [f64; 3] {
    if c == ZERO_CONE {
        return [0.0; 3];
    }
    let axis = c.0 / 2;
    let mut u = [0.0; 3];
    u[axis] = if c.0 % 2 == 0 { 1.0 } else { -1.0 };
    u[(axis + 1) % 3] = 2.0 * (c.1 as f64 + 0.5) / q as f64 - 1.0;
    u[(axis + 2) % 3] = 2.0 * (c.2 as f64 + 0.5) / q as f64 - 1.0;
    let r = det::sqrt(kernels::dot(u, u));
    u.map(|v| v / r)
}

impl DirectionalPlan {
    fn build(op: &FmmC64<'_>, options: DirectionalOptions, k: f64) -> Result<Self, FmmError> {
        let mut bins = Vec::new();
        let mut counts = Vec::new();
        for level in 0..=op.geometry.max_level {
            let kh = k * op.geometry.side / f64::from(1u32 << level);
            let q = if kh <= options.low_freq_threshold {
                0
            } else {
                (kh * options.cone_resolution as f64).ceil() as usize
            };
            bins.push(q);
            // Empty levels allocate no directional basis. Admit an angular
            // grid only when an interaction or nested basis actually uses it.
            counts.push(0);
        }
        let mut plan = Self {
            stats: DirectionalStats {
                cones_per_level: counts,
                source_expansions: 0,
                target_expansions: 0,
                m2l_pairs: 0,
                high_frequency_pairs: 0,
                near_pairs: 0,
                pair_scans: 0,
                m2l_work: 0,
                transfer_work: 0,
                leaf_work: 0,
                coefficient_slots: 0,
                memory_bytes: 0,
                order: options.order,
                requested_tolerance: options.tolerance,
                interpolation_indicator: det::powi(4.0, -(options.order as i32)),
                coupling_blocks: 0,
                maximum_rank: 0,
                compression_residual: 0.0,
                setup_work: 0,
            },
            options,
            bins,
            pairs: Vec::new(),
            near: Vec::new(),
            sources: BTreeSet::new(),
            targets: BTreeSet::new(),
            k,
            coupling_keys: Vec::new(),
            couplings: BTreeMap::new(),
        };
        let root = (0, 0, 0, 0);
        let mut stack = vec![(root, root)];
        let (ns, nt) = op.operator.dimensions();
        let p3 = op.nodes.len().pow(3);
        while let Some((a, b)) = stack.pop() {
            plan.stats.pair_scans += 1;
            budget(
                "directional cell-pair scan",
                plan.stats.pair_scans,
                MAX_CELL_PAIR_SCANS,
            )?;
            let h = op.geometry.cell_box(a).1;
            let d = sub(center(op, a), center(op, b));
            let distance = det::sqrt(kernels::dot(d, d));
            let gap = distance - det::sqrt(3.0) * h;
            // Geometric separation plus the parabolic condition k*diam^2/gap
            // controls curvature of r-u.r. Nonadmissible pairs descend, and
            // only unseparated leaf pairs use P2P.
            let separated =
                distance >= 4.0 * h && (plan.bins[a.0 as usize] == 0 || 3.0 * k * h * h <= gap);
            if separated {
                let c = cone(d, plan.bins[a.0 as usize]);
                plan.admit_cone_level(a.0)?;
                plan.pairs.push((a, b, c));
                plan.sources.insert((b, c));
                plan.targets.insert((a, c));
                if c != ZERO_CONE {
                    plan.stats.high_frequency_pairs += 1;
                }
            } else if op.geometry.cells[&a].leaf {
                let pairs = op.geometry.cells[&a]
                    .points
                    .len()
                    .saturating_mul(op.geometry.cells[&b].points.len());
                plan.stats.near_pairs = plan.stats.near_pairs.saturating_add(pairs);
                budget(
                    "directional near-field point pair",
                    plan.stats.near_pairs,
                    plan.options.max_near_pairs,
                )?;
                plan.near.push((a, b));
            } else {
                let ca = children(op, a);
                let cb = children(op, b);
                for &aa in ca.iter().rev() {
                    for &bb in cb.iter().rev() {
                        stack.push((aa, bb));
                    }
                }
            }
        }
        let (sw, sl) = plan.close_basis(op, true, ns)?;
        let (tw, tl) = plan.close_basis(op, false, nt)?;
        plan.stats.transfer_work = sw.saturating_add(tw);
        plan.stats.leaf_work = sl.saturating_add(tl);
        budget(
            "directional total translation",
            plan.stats
                .m2l_work
                .saturating_add(plan.stats.transfer_work)
                .saturating_add(plan.stats.leaf_work),
            MAX_TRANSLATION_WORK,
        )?;
        plan.stats.coefficient_slots = (plan.sources.len().saturating_mul(ns)
            + plan.targets.len().saturating_mul(nt)
            + 4 * (ns + nt))
            .saturating_mul(p3);
        FmmC64::validate_coefficient_slots_within(
            plan.stats.coefficient_slots,
            plan.options.max_coefficient_slots,
        )?;
        plan.prepare_couplings(op)?;
        budget(
            "directional total translation",
            plan.stats
                .m2l_work
                .saturating_add(plan.stats.transfer_work)
                .saturating_add(plan.stats.leaf_work),
            MAX_TRANSLATION_WORK,
        )?;
        plan.stats.source_expansions = plan.sources.len();
        plan.stats.target_expansions = plan.targets.len();
        plan.stats.m2l_pairs = plan.pairs.len();
        plan.stats.memory_bytes = plan.stats.coefficient_slots * core::mem::size_of::<C64>()
            + (plan.sources.len() + plan.targets.len()) * core::mem::size_of::<Expansion>()
            + plan.pairs.len() * core::mem::size_of::<(CellKey, CellKey, Cone)>()
            + plan.near.len() * core::mem::size_of::<(CellKey, CellKey)>()
            + plan.coupling_keys.len() * core::mem::size_of::<CouplingKey>();
        Ok(plan)
    }
    fn close_basis(
        &mut self,
        op: &FmmC64<'_>,
        source: bool,
        channels: usize,
    ) -> Result<(usize, usize), FmmError> {
        let records = if source {
            &mut self.sources
        } else {
            &mut self.targets
        };
        let mut queue: Vec<_> = records.iter().copied().collect();
        let mut pos = 0;
        let (mut work, mut leaf_work) = (0usize, 0usize);
        let p3 = op.nodes.len().pow(3);
        while pos < queue.len() {
            let (key, c) = queue[pos];
            pos += 1;
            let q = self.bins[key.0 as usize];
            let count = if q == 0 {
                1
            } else {
                6usize.saturating_mul(q.saturating_mul(q))
            };
            budget(
                "directional cones per level",
                count,
                self.options.max_cones_per_level,
            )?;
            self.stats.cones_per_level[key.0 as usize] = count;
            if op.geometry.cells[&key].leaf {
                leaf_work = leaf_work.saturating_add(
                    op.geometry.cells[&key]
                        .points
                        .len()
                        .saturating_mul(p3 * channels),
                );
            } else {
                let u = direction(c, self.bins[key.0 as usize]);
                for child in children(op, key) {
                    work = work.saturating_add(p3 * p3 * channels);
                    let v = cone(u, self.bins[child.0 as usize]);
                    if records.insert((child, v)) {
                        queue.push((child, v));
                    }
                }
            }
            budget(
                "directional basis coefficient-slot",
                records.len().saturating_mul(p3 * channels),
                self.options.max_coefficient_slots,
            )?;
            budget("directional basis transfer", work, MAX_TRANSLATION_WORK)?;
        }
        Ok((work, leaf_work))
    }
    fn admit_cone_level(&mut self, level: u32) -> Result<(), FmmError> {
        let q = self.bins[level as usize];
        let count = if q == 0 {
            1
        } else {
            6usize.saturating_mul(q.saturating_mul(q))
        };
        budget(
            "directional cones per level",
            count,
            self.options.max_cones_per_level,
        )?;
        self.stats.cones_per_level[level as usize] = count;
        Ok(())
    }
    fn coupling_key(&self, op: &FmmC64<'_>, a: CellKey, b: CellKey, c: Cone) -> CouplingKey {
        let invariant = match &op.operator {
            Operator::Scalar(k) => k.translation_invariant(),
            _ => true,
        };
        if !invariant {
            return (a, b, c);
        }
        let zig = |d: i64| {
            if d >= 0 {
                (2 * d) as u32
            } else {
                (-2 * d - 1) as u32
            }
        };
        (
            (a.0, 0, 0, 0),
            (
                b.0,
                zig(i64::from(a.1) - i64::from(b.1)),
                zig(i64::from(a.2) - i64::from(b.2)),
                zig(i64::from(a.3) - i64::from(b.3)),
            ),
            c,
        )
    }
    fn prepare_couplings(&mut self, op: &FmmC64<'_>) -> Result<(), FmmError> {
        let (ns, nt) = op.operator.dimensions();
        let n = op.nodes.len().pow(3);
        // Includes simultaneous pass coefficients, every retained factor,
        // one residual matrix, and rank-vector scratch before allocation.
        let mut retained = self.stats.coefficient_slots;
        let mut uses = BTreeMap::new();
        for &(a, b, c) in &self.pairs {
            *uses.entry(self.coupling_key(op, a, b, c)).or_insert(0usize) += 1;
        }
        for index in 0..self.pairs.len() {
            let (a, b, c) = self.pairs[index];
            let key = self.coupling_key(op, a, b, c);
            if !self.couplings.contains_key(&key) {
                if uses[&key] < 2 {
                    // A one-use block keeps W5a's matrix-free tensor M2L.
                    // No compression setup/storage is paid for no reuse.
                    self.couplings.insert(key, Coupling::TensorGrid);
                } else {
                    let u = direction(c, self.bins[a.0 as usize]);
                    let ca = center(op, a);
                    let cb = center(op, b);
                    let ga = op.grid(a);
                    let gb = op.grid(b);
                    let mut channels = Vec::new();
                    for channel in 0..ns * nt {
                        FmmC64::validate_coefficient_slots_within(
                            retained.saturating_add(n * n).saturating_add(2 * n),
                            self.options.max_coefficient_slots,
                        )?;
                        self.stats.setup_work = self.stats.setup_work.saturating_add(n * n);
                        budget(
                            "directional coupling setup",
                            self.stats.setup_work,
                            self.options.max_setup_work,
                        )?;
                        let mut matrix = vec![C64::ZERO; n * n];
                        for (i, &x) in ga.iter().enumerate() {
                            for (j, &y) in gb.iter().enumerate() {
                                let value = modulate(
                                    op.operator.block(x, y)[channel],
                                    wave(-self.k, u, sub(sub(x, ca), sub(y, cb))),
                                );
                                if !finite(value) {
                                    return Err(FmmError::NonFiniteOutput { index: i });
                                }
                                matrix[i * n + j] = value;
                            }
                        }
                        let factor = compress(
                            matrix,
                            n,
                            self.options.compression_tolerance,
                            self.options.max_rank,
                            Envelope {
                                retained,
                                max_slots: self.options.max_coefficient_slots,
                                max_setup: self.options.max_setup_work,
                            },
                            &mut self.stats.setup_work,
                        )?;
                        retained = retained.saturating_add(factor.slots());
                        self.stats.maximum_rank = self.stats.maximum_rank.max(factor.rank());
                        self.stats.compression_residual =
                            self.stats.compression_residual.max(factor.residual);
                        channels.push(factor);
                    }
                    self.couplings.insert(key, Coupling::Compressed(channels));
                }
            }
            let work = match &self.couplings[&key] {
                Coupling::TensorGrid => n * n * ns * nt,
                Coupling::Compressed(channels) => {
                    channels.iter().map(|f| 2 * n * f.rank()).sum::<usize>()
                }
            };
            self.stats.m2l_work = self.stats.m2l_work.saturating_add(work);
            budget(
                "directional compressed M2L",
                self.stats.m2l_work,
                MAX_TRANSLATION_WORK,
            )?;
            self.coupling_keys.push(key);
        }
        self.stats.coupling_blocks = self.couplings.len();
        self.stats.coefficient_slots = retained;
        Ok(())
    }
    pub(super) fn potentials(
        &self,
        op: &FmmC64<'_>,
        charges: &[C64],
    ) -> Result<Vec<C64>, FmmError> {
        let (ns, nt) = op.operator.dimensions();
        let p3 = op.nodes.len().pow(3);
        let mut multipoles: BTreeMap<Expansion, Vec<C64>> = BTreeMap::new();
        for &(key, c) in self.sources.iter().rev() {
            let u = direction(c, self.bins[key.0 as usize]);
            let origin = center(op, key);
            let mut m = vec![C64::ZERO; p3 * ns];
            if op.geometry.cells[&key].leaf {
                for &idx in &op.geometry.cells[&key].points {
                    let w = op.weights(key, op.geometry.points[idx]);
                    let phase = wave(-self.k, u, sub(op.geometry.points[idx], origin));
                    for a in 0..ns {
                        let q = modulate(charges[idx].scale(op.source_weight(idx, a)), phase);
                        for (j, &weight) in w.iter().enumerate() {
                            m[ns * j + a] = m[ns * j + a] + q.scale(weight);
                        }
                    }
                }
            } else {
                for child in children(op, key) {
                    let cc = cone(u, self.bins[child.0 as usize]);
                    let v = direction(cc, self.bins[child.0 as usize]);
                    let cm = &multipoles[&(child, cc)];
                    let child_origin = center(op, child);
                    for (j, point) in op.grid(child).into_iter().enumerate() {
                        let weights = op.weights(key, point);
                        let phase = wave(self.k, v, sub(point, child_origin))
                            * wave(-self.k, u, sub(point, origin));
                        for (i, &w) in weights.iter().enumerate() {
                            for a in 0..ns {
                                m[ns * i + a] =
                                    m[ns * i + a] + modulate(cm[ns * j + a], phase).scale(w);
                            }
                        }
                    }
                }
            }
            multipoles.insert((key, c), m);
        }
        let mut locals: BTreeMap<Expansion, Vec<C64>> = self
            .targets
            .iter()
            .map(|&r| (r, vec![C64::ZERO; p3 * nt]))
            .collect();
        for (&(a, b, c), key) in self.pairs.iter().zip(&self.coupling_keys) {
            let m = &multipoles[&(b, c)];
            let local = locals.get_mut(&(a, c)).expect("admitted target expansion");
            match &self.couplings[key] {
                Coupling::TensorGrid => {
                    let u = direction(c, self.bins[a.0 as usize]);
                    let ca = center(op, a);
                    let cb = center(op, b);
                    let ga = op.grid(a);
                    let gb = op.grid(b);
                    for (i, &x) in ga.iter().enumerate() {
                        for (j, &y) in gb.iter().enumerate() {
                            let block = op.operator.block(x, y);
                            if block[..ns * nt].iter().any(|&z| !finite(z)) {
                                return Err(FmmError::NonFiniteOutput { index: i });
                            }
                            let phase = wave(-self.k, u, sub(sub(x, ca), sub(y, cb)));
                            for ta in 0..nt {
                                for sa in 0..ns {
                                    local[nt * i + ta] = local[nt * i + ta]
                                        + modulate(block[ns * ta + sa], phase) * m[ns * j + sa];
                                }
                            }
                        }
                    }
                }
                Coupling::Compressed(channels) => {
                    for (channel, factor) in channels.iter().enumerate() {
                        let (ta, sa) = (channel / ns, channel % ns);
                        for (u, v) in factor.u.iter().zip(&factor.v) {
                            let mut value = C64::ZERO;
                            for (j, &vj) in v.iter().enumerate() {
                                value = value + vj * m[ns * j + sa];
                            }
                            for (i, &ui) in u.iter().enumerate() {
                                local[nt * i + ta] = local[nt * i + ta] + ui * value;
                            }
                        }
                    }
                }
            }
        }
        let mut out = vec![C64::ZERO; charges.len()];
        for &(key, c) in &self.targets {
            let u = direction(c, self.bins[key.0 as usize]);
            let origin = center(op, key);
            if !op.geometry.cells[&key].leaf {
                let lp = locals[&(key, c)].clone();
                for child in children(op, key) {
                    let cc = cone(u, self.bins[child.0 as usize]);
                    let v = direction(cc, self.bins[child.0 as usize]);
                    let child_origin = center(op, child);
                    let local = locals
                        .get_mut(&(child, cc))
                        .expect("admitted child expansion");
                    for (i, point) in op.grid(child).into_iter().enumerate() {
                        let weights = op.weights(key, point);
                        let phase = wave(self.k, u, sub(point, origin))
                            * wave(-self.k, v, sub(point, child_origin));
                        for (j, &w) in weights.iter().enumerate() {
                            for a in 0..nt {
                                local[nt * i + a] =
                                    local[nt * i + a] + modulate(lp[nt * j + a], phase).scale(w);
                            }
                        }
                    }
                }
            } else {
                let local = &locals[&(key, c)];
                for &idx in &op.geometry.cells[&key].points {
                    let weights = op.weights(key, op.geometry.points[idx]);
                    let phase = wave(self.k, u, sub(op.geometry.points[idx], origin));
                    for (j, &w) in weights.iter().enumerate() {
                        for a in 0..nt {
                            out[idx] = out[idx]
                                + modulate(local[nt * j + a], phase)
                                    .scale(w * op.target_weight(idx, a));
                        }
                    }
                }
            }
        }
        // The same exact normal-factored point kernel as W5a; no quadrature,
        // self term, charge modification or global direct fallback.
        for &(a, b) in &self.near {
            for &t in &op.geometry.cells[&a].points {
                for &s in &op.geometry.cells[&b].points {
                    if t != s {
                        out[t] = out[t] + op.pair(t, s) * charges[s];
                    }
                }
            }
        }
        validate_output(out)
    }
}

fn frobenius(values: &[C64]) -> f64 {
    let scale = values
        .iter()
        .map(|v| v.re.abs().max(v.im.abs()))
        .fold(0.0, f64::max);
    if scale == 0.0 {
        0.0
    } else {
        scale
            * det::sqrt(
                values
                    .iter()
                    .map(|v| C64::new(v.re / scale, v.im / scale).norm_sq())
                    .sum::<f64>(),
            )
    }
}
/// Coefficient slots already retained before a block is compressed, and the
/// declared slot and setup-work ceilings it must stay within.
#[derive(Clone, Copy)]
struct Envelope {
    retained: usize,
    max_slots: usize,
    max_setup: usize,
}
fn compress(
    mut residual: Vec<C64>,
    n: usize,
    tolerance: f64,
    max_rank: usize,
    envelope: Envelope,
    setup_work: &mut usize,
) -> Result<LowRank, FmmError> {
    let mut factor = LowRank {
        u: Vec::new(),
        v: Vec::new(),
        residual: 0.0,
    };
    let scale = residual
        .iter()
        .map(|v| v.re.abs().max(v.im.abs()))
        .fold(0.0, f64::max);
    if scale == 0.0 {
        return Ok(factor);
    }
    // Normalize before norms/pivot reciprocals, including subnormal or very
    // large finite kernels. An infinite norm must never accept a zero rank.
    for z in &mut residual {
        *z = C64::new(z.re / scale, z.im / scale);
    }
    let norm = frobenius(&residual);
    let mut current = norm;
    while current > tolerance * norm {
        let requested = factor.rank() + 1;
        budget("directional coupling rank", requested, max_rank.min(n))?;
        FmmC64::validate_coefficient_slots_within(
            envelope
                .retained
                .saturating_add(n * n)
                .saturating_add(2 * n * requested),
            envelope.max_slots,
        )?;
        *setup_work = setup_work.saturating_add(n * n);
        budget(
            "directional coupling setup",
            *setup_work,
            envelope.max_setup,
        )?;
        // Full residual pivoting, fixed lowest-index ties. Unlike heuristic
        // partial ACA stopping, the entire residual matrix is measured.
        let mut pivot_index = 0;
        let mut largest = 0.0;
        for (i, value) in residual.iter().enumerate() {
            let size = value.abs();
            if size > largest {
                largest = size;
                pivot_index = i;
            }
        }
        if largest == 0.0 {
            break;
        }
        let (row, col) = (pivot_index / n, pivot_index % n);
        let inverse = residual[pivot_index].recip();
        let u: Vec<_> = (0..n).map(|i| residual[i * n + col]).collect();
        let v: Vec<_> = (0..n).map(|j| residual[row * n + j] * inverse).collect();
        for i in 0..n {
            for j in 0..n {
                residual[i * n + j] = residual[i * n + j] - u[i] * v[j];
            }
        }
        current = frobenius(&residual);
        if !current.is_finite() {
            return Err(FmmError::NonFiniteOutput { index: pivot_index });
        }
        factor.u.push(u);
        factor.v.push(v);
    }
    factor.residual = current / norm;
    for (i, u) in factor.u.iter_mut().enumerate() {
        for z in u {
            *z = z.scale(scale);
            if !finite(*z) {
                return Err(FmmError::NonFiniteOutput { index: i });
            }
        }
    }
    Ok(factor)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn v1(retained: usize) -> Envelope {
        Envelope {
            retained,
            max_slots: MAX_COMPLEX_COEFFICIENT_SLOTS,
            max_setup: MAX_TRANSLATION_WORK,
        }
    }
    #[test]
    fn directional_compression_singular_value_ladder_and_caps() {
        let n = 8;
        let mut matrix = vec![C64::ZERO; n * n];
        for i in 0..n {
            matrix[i * n + i] = C64::from_re(det::powi(2.0, -(3 * i as i32)));
        }
        let total = (0..n).map(|i| det::powi(2.0, -(6 * i as i32))).sum::<f64>();
        for (rung, tol) in [0.2, 0.02, 0.002].into_iter().enumerate() {
            let rank = rung + 1;
            let mut work = 0;
            let f = compress(matrix.clone(), n, tol, 8, v1(0), &mut work).unwrap();
            let expected = det::sqrt(
                (rank..n)
                    .map(|i| det::powi(2.0, -(6 * i as i32)))
                    .sum::<f64>()
                    / total,
            );
            let error = det::sqrt(
                (0..n)
                    .flat_map(|i| (0..n).map(move |j| (i, j)))
                    .map(|(i, j)| {
                        let approximated =
                            f.u.iter()
                                .zip(&f.v)
                                .fold(C64::ZERO, |sum, (u, v)| sum + u[i] * v[j]);
                        (matrix[i * n + j] - approximated).norm_sq()
                    })
                    .sum::<f64>()
                    / total,
            );
            assert_eq!(f.rank(), rank);
            assert!((error - expected).abs() < 1e-15);
            assert!((f.residual - expected).abs() < 1e-15);
        }
        let mut work = 0;
        for scale in [1e-200, 1e200] {
            let scaled: Vec<_> = matrix.iter().map(|z| z.scale(scale)).collect();
            let f = compress(scaled, n, 0.02, 8, v1(0), &mut 0).unwrap();
            assert_eq!(f.rank(), 2);
            let first =
                f.u.iter()
                    .zip(&f.v)
                    .fold(C64::ZERO, |s, (u, v)| s + u[0] * v[0]);
            assert!((first.re / scale - 1.0).abs() < 1e-14 && first.im == 0.0);
        }
        assert!(matches!(
            compress(matrix.clone(), n, 0.002, 2, v1(0), &mut work),
            Err(FmmError::WorkEnvelopeExceeded {
                operation: "directional coupling rank",
                requested: 3,
                max: 2
            })
        ));
        let retained = MAX_COMPLEX_COEFFICIENT_SLOTS - n * n - 2 * n;
        assert!(compress(matrix.clone(), n, 0.2, 1, v1(retained), &mut 0).is_ok());
        assert!(matches!(
            compress(matrix.clone(), n, 0.2, 1, v1(retained + 1), &mut 0),
            Err(FmmError::WorkEnvelopeExceeded { .. })
        ));
        assert!(
            compress(
                matrix.clone(),
                n,
                0.2,
                1,
                v1(0),
                &mut (MAX_TRANSLATION_WORK - n * n)
            )
            .is_ok()
        );
        assert!(matches!(
            compress(
                matrix.clone(),
                n,
                0.2,
                1,
                v1(0),
                &mut (MAX_TRANSLATION_WORK - n * n + 1)
            ),
            Err(FmmError::WorkEnvelopeExceeded {
                operation: "directional coupling setup",
                ..
            })
        ));
        // A declared envelope moves both caps, below or above the v1 ceilings.
        let declared = |max_slots, max_setup| Envelope {
            retained: 0,
            max_slots,
            max_setup,
        };
        let slots = n * n + 2 * n;
        assert!(compress(matrix.clone(), n, 0.2, 1, declared(slots, n * n), &mut 0).is_ok());
        assert!(matches!(
            compress(matrix.clone(), n, 0.2, 1, declared(slots - 1, n * n), &mut 0),
            Err(FmmError::WorkEnvelopeExceeded {
                operation: "complex coefficient-slot (two f64 per C64)",
                requested,
                max,
            }) if requested == 2 * slots && max == 2 * (slots - 1)
        ));
        assert!(matches!(
            compress(matrix.clone(), n, 0.2, 1, declared(slots, n * n - 1), &mut 0),
            Err(FmmError::WorkEnvelopeExceeded {
                operation: "directional coupling setup",
                requested,
                max,
            }) if requested == n * n && max == n * n - 1
        ));
        let mut past_v1 = MAX_TRANSLATION_WORK;
        let above = MAX_TRANSLATION_WORK + n * n;
        assert!(compress(matrix, n, 0.2, 1, declared(slots, above), &mut past_v1).is_ok());
        assert_eq!(past_v1, above);
    }
    #[test]
    fn directional_pair_partition_and_low_frequency_identity() {
        let points: Vec<_> = (0..20)
            .map(|i| {
                [
                    i as f64 / 20.0,
                    ((i * 3) % 7) as f64 / 20.0,
                    ((i * 5) % 11) as f64 / 20.0,
                ]
            })
            .collect();
        let kernel = crate::HelmholtzSingle { k: 130.0 };
        let mut options = DirectionalOptions::for_tolerance(1e-3).unwrap();
        options.low_freq_threshold = 1.0;
        let op = FmmC64::new_directional(&kernel, &points, options).unwrap();
        let plan = op.directional.as_ref().unwrap();
        let descendants = |root: CellKey| {
            op.geometry
                .cells
                .iter()
                .filter(|(key, cell)| {
                    cell.leaf
                        && key.0 >= root.0
                        && (key.1 >> (key.0 - root.0)) == root.1
                        && (key.2 >> (key.0 - root.0)) == root.2
                        && (key.3 >> (key.0 - root.0)) == root.3
                })
                .flat_map(|(_, cell)| cell.points.iter().copied())
                .collect::<Vec<_>>()
        };
        let mut cover = vec![0; points.len() * points.len()];
        for (a, b) in plan
            .pairs
            .iter()
            .map(|&(a, b, _)| (a, b))
            .chain(plan.near.iter().copied())
        {
            for i in descendants(a) {
                for j in descendants(b) {
                    cover[i * points.len() + j] += 1;
                }
            }
        }
        assert!(cover.iter().all(|&v| v == 1));
        for q in [1, 2, 7, 32] {
            for u in [
                [1., 0., 0.],
                [-1., 0.2, 0.3],
                [0.1, -1., 0.2],
                [0.4, 0.5, 1.],
            ] {
                let c = cone(u, q);
                let v = direction(c, q);
                assert!((kernels::dot(v, v) - 1.0).abs() < 1e-14);
                assert_eq!(cone(v, q), c);
            }
        }
        let low = crate::HelmholtzSingle { k: 0.25 };
        let q = vec![C64::ONE; points.len()];
        let opts = DirectionalOptions::for_tolerance(1e-3).unwrap();
        let regular = FmmC64::new(&low, points.clone(), opts.order, opts.leaf_capacity).unwrap();
        let directional = FmmC64::new_directional(&low, &points, opts).unwrap();
        assert!(directional.directional_stats().is_none());
        let bits = |v: Vec<C64>| {
            v.into_iter()
                .map(|z| (z.re.to_bits(), z.im.to_bits()))
                .collect::<Vec<_>>()
        };
        assert_eq!(
            bits(regular.potentials(&q).unwrap()),
            bits(directional.potentials(&q).unwrap())
        );
    }
}
