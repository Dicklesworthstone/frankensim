//! Bounded 3-D bulk quadrature on a raw implicit domain `{phi < 0}`.
//!
//! Certified box signs select full/empty cells. Cut boxes subdivide, then use
//! a derivative-certified monotone height direction and positive tensor Gauss
//! weights. A cut height is enclosed, including crossings exactly on an edge;
//! an uncertain point sign never becomes a sampled root or an empty domain.
//!
//! The rule is NUMERICAL, not an enclosure of arbitrary integrals. Separately,
//! `volume_bounds` encloses domain volume by counting every unresolved cut box
//! in full. It does not confuse tiny root brackets with a base-quadrature error
//! bound. No surface rule, octree balancing, or high-order convergence claim.

use std::ops::ControlFlow;

use crate::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;

const AXES: [HeightAxis; 3] = [HeightAxis::X, HeightAxis::Y, HeightAxis::Z];
const GAUSS: [(f64, f64); 3] = [
    (-0.774_596_669_241_483_4, 0.555_555_555_555_555_6),
    (0.0, 0.888_888_888_888_889),
    (0.774_596_669_241_483_4, 0.555_555_555_555_555_6),
];

/// Limits for the entire family of calls sharing one control.
#[derive(Debug, Clone, Copy)]
pub struct QuadratureOptions3 {
    /// Octant subdivisions before the monotone-height rule (at most 12).
    pub depth: u32,
    /// Maximum visited boxes, including fully classified boxes.
    pub max_boxes: usize,
    /// Maximum emitted bulk points, across all calls and failed trials.
    pub max_points: usize,
    /// Maximum interval contraction steps in each height line (1..=96).
    pub root_iterations: usize,
    /// Root-bracket width relative to the cut box's height span.
    pub root_relative_tolerance: f64,
}

impl Default for QuadratureOptions3 {
    fn default() -> Self {
        Self {
            depth: 2, max_boxes: 1_000_000, max_points: 4_000_000,
            root_iterations: 64, root_relative_tolerance: 1e-10,
        }
    }
}

/// Work already spent, including work discarded by a returned error.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct QuadratureWork3 {
    /// Visited boxes.
    pub boxes: usize,
    /// Calls into the implicit-field interval producers.
    pub field_evaluations: usize,
    /// Emitted points; discarded candidates still consume the allowance.
    pub points: usize,
}

/// No partial rule is returned on any of these outcomes.
#[derive(Debug, Clone, PartialEq)]
pub enum QuadratureError3 {
    /// Invalid or unrepresentable arithmetic/configuration.
    Invalid(&'static str),
    /// Caller stopped at a box, field, root, or publication boundary.
    Cancelled,
    /// Shared box allowance exhausted.
    BoxBudget,
    /// Shared point allowance exhausted.
    PointBudget,
    /// Finite derivative enclosures do not prove any monotone direction.
    UnresolvedCell(HexCell),
    /// The retained height bracket did not reach the requested resolution.
    RootResolution(HexCell),
}

impl std::fmt::Display for QuadratureError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "3-D cut quadrature refused: {self:?}")
    }
}
impl std::error::Error for QuadratureError3 {}

/// Reusable cumulative quadrature budget and cooperative cancellation context.
pub struct QuadratureControl3<'a> {
    options: QuadratureOptions3,
    work: QuadratureWork3,
    callback: &'a mut dyn FnMut(QuadratureWork3) -> ControlFlow<()>,
}

impl<'a> QuadratureControl3<'a> {
    /// Validate settings before any field callback or quadrature allocation.
    pub fn new(
        options: QuadratureOptions3,
        callback: &'a mut dyn FnMut(QuadratureWork3) -> ControlFlow<()>,
    ) -> Result<Self, QuadratureError3> {
        if options.depth > 12 || !(1..=96).contains(&options.root_iterations)
            || !options.root_relative_tolerance.is_finite()
            || options.root_relative_tolerance <= 0.0
            || options.root_relative_tolerance >= 1.0 {
            return Err(QuadratureError3::Invalid("invalid depth/root budget or tolerance"));
        }
        Ok(Self { options, work: QuadratureWork3::default(), callback })
    }

    /// Consumed work is observable even after cancellation or refusal.
    #[must_use]
    pub const fn work(&self) -> QuadratureWork3 { self.work }

    pub(crate) fn poll(&mut self) -> Result<(), QuadratureError3> {
        match (self.callback)(self.work) {
            ControlFlow::Continue(()) => Ok(()),
            ControlFlow::Break(()) => Err(QuadratureError3::Cancelled),
        }
    }

    fn field(&mut self, evaluate: impl FnOnce() -> Interval) -> Result<Interval, QuadratureError3> {
        self.poll()?;
        self.work.field_evaluations = self.work.field_evaluations.checked_add(1)
            .ok_or(QuadratureError3::Invalid("field counter overflow"))?;
        let interval = evaluate();
        self.poll()?;
        if !interval.lo().is_finite() || !interval.hi().is_finite() || interval.lo() > interval.hi() {
            return Err(QuadratureError3::Invalid("nonfinite or unordered field enclosure"));
        }
        Ok(interval)
    }
}

/// Retained positive bulk rule and an independent conservative volume bracket.
#[derive(Debug, Clone)]
pub struct CutRules3 {
    bulk: Vec<([f64; 3], f64)>,
    volume_bounds: Interval,
    cut_boxes: usize,
}

impl CutRules3 {
    /// Global positions and positive volume weights. Weights are numerical.
    #[must_use]
    pub fn bulk(&self) -> &[([f64; 3], f64)] { &self.bulk }
    /// Domain-volume enclosure conditional on the supplied SDF enclosures.
    /// Cut boxes contribute [0, box volume], not their sampled Gauss sum.
    #[must_use]
    pub const fn volume_bounds(&self) -> Interval { self.volume_bounds }
    /// Numerically integrated material volume, distinct from the enclosure.
    #[must_use]
    pub fn volume(&self) -> f64 { self.bulk.iter().map(|(_, w)| w).sum() }
    /// Leaf boxes whose geometry needed a height rule.
    #[must_use]
    pub const fn cut_boxes(&self) -> usize { self.cut_boxes }
}

/// Construct a 3-D bulk rule; no partial points escape on a returned failure.
/// The same `control` may be used across an entire Cartesian background grid.
pub fn cut_cell_rules3(
    sdf: &dyn CutSdf3, cell: HexCell, control: &mut QuadratureControl3<'_>,
) -> Result<CutRules3, QuadratureError3> {
    control.poll()?;
    let mut result = CutRules3 { bulk: Vec::new(), volume_bounds: Interval::new(0.0, 0.0), cut_boxes: 0 };
    visit(sdf, cell, control.options.depth, control, &mut result)?;
    // Outward additions around exact zero can produce a negative subnormal.
    result.volume_bounds = Interval::new(result.volume_bounds.lo().max(0.0), result.volume_bounds.hi());
    if !result.volume().is_finite() || !result.volume_bounds.hi().is_finite() {
        return Err(QuadratureError3::Invalid("quadrature volume overflow"));
    }
    control.poll()?;
    Ok(result)
}

fn box_volume(cell: HexCell) -> Result<Interval, QuadratureError3> {
    let (lo, hi) = (cell.lo(), cell.hi());
    let mut result = Interval::new(1.0, 1.0);
    for axis in 0..3 {
        let span = Interval::new(hi[axis], hi[axis]) - Interval::new(lo[axis], lo[axis]);
        if !span.hi().is_finite() || span.lo() <= 0.0 {
            return Err(QuadratureError3::Invalid("unrepresentable box span"));
        }
        result = result * span;
    }
    if !result.hi().is_finite() || result.lo() <= 0.0 {
        return Err(QuadratureError3::Invalid("unrepresentable box volume"));
    }
    Ok(result)
}

fn point(out: &mut CutRules3, p: [f64; 3], weight: f64, control: &mut QuadratureControl3<'_>)
    -> Result<(), QuadratureError3> {
    if !p.iter().all(|v| v.is_finite()) || !weight.is_finite() || weight <= 0.0 {
        return Err(QuadratureError3::Invalid("unrepresentable Gauss point/weight"));
    }
    if control.work.points >= control.options.max_points { return Err(QuadratureError3::PointBudget); }
    control.work.points += 1;
    out.bulk.push((p, weight));
    Ok(())
}

fn full_box(cell: HexCell, control: &mut QuadratureControl3<'_>, out: &mut CutRules3)
    -> Result<(), QuadratureError3> {
    let (lo, hi) = (cell.lo(), cell.hi());
    let m = std::array::from_fn::<_, 3, _>(|i| f64::midpoint(lo[i], hi[i]));
    let s = std::array::from_fn::<_, 3, _>(|i| 0.5 * (hi[i] - lo[i]));
    for (x, wx) in GAUSS { for (y, wy) in GAUSS { for (z, wz) in GAUSS {
        point(out, [m[0] + s[0]*x, m[1] + s[1]*y, m[2] + s[2]*z],
            wx*wy*wz*s[0]*s[1]*s[2], control)?;
    } } }
    Ok(())
}

fn visit(sdf: &dyn CutSdf3, cell: HexCell, depth: u32, control: &mut QuadratureControl3<'_>, out: &mut CutRules3)
    -> Result<(), QuadratureError3> {
    control.poll()?;
    if control.work.boxes >= control.options.max_boxes { return Err(QuadratureError3::BoxBudget); }
    control.work.boxes += 1;
    let volume = box_volume(cell)?;
    let (lo, hi) = (cell.lo(), cell.hi());
    let sign = control.field(|| sdf.enclose(lo, hi))?;
    if sign.lo() >= 0.0 { return Ok(()); }
    if sign.hi() < 0.0 {
        out.volume_bounds = out.volume_bounds + volume;
        return full_box(cell, control, out);
    }
    if depth > 0 {
        let m = std::array::from_fn::<_, 3, _>(|a| f64::midpoint(lo[a], hi[a]));
        if (0..3).any(|a| m[a] <= lo[a] || m[a] >= hi[a]) {
            return Err(QuadratureError3::Invalid("subdivision cannot advance"));
        }
        for octant in 0..8 {
            let a = std::array::from_fn(|i| if octant & (1 << i) == 0 { lo[i] } else { m[i] });
            let b = std::array::from_fn(|i| if octant & (1 << i) == 0 { m[i] } else { hi[i] });
            let child = HexCell::try_new(a, b).map_err(|_| QuadratureError3::Invalid("invalid subdivision"))?;
            visit(sdf, child, depth - 1, control, out)?;
        }
        return Ok(());
    }
    let mut selected = None;
    let mut best = 0.0;
    for (a, axis) in AXES.into_iter().enumerate() {
        // Unusable derivatives disqualify only that direction. A field
        // enclosure itself must still be finite to classify the box.
        control.poll()?;
        control.work.field_evaluations = control.work.field_evaluations.checked_add(1)
            .ok_or(QuadratureError3::Invalid("field counter overflow"))?;
        let d = sdf.derivative_enclose(lo, hi, axis);
        control.poll()?;
        let bound = if d.lo().is_finite() && d.hi().is_finite() && d.lo() <= d.hi() {
            if d.lo() > 0.0 { d.lo() } else if d.hi() < 0.0 { -d.hi() } else { 0.0 }
        } else { 0.0 };
        if bound > best { best = bound; selected = Some((a, d.lo() > 0.0)); }
    }
    let (axis, increasing) = selected.ok_or(QuadratureError3::UnresolvedCell(cell))?;
    out.cut_boxes += 1;
    out.volume_bounds = out.volume_bounds + Interval::new(0.0, volume.hi());
    let bases: Vec<usize> = (0..3).filter(|&a| a != axis).collect();
    let (a, b) = (bases[0], bases[1]);
    let ma = f64::midpoint(lo[a], hi[a]); let mb = f64::midpoint(lo[b], hi[b]);
    let sa = 0.5 * (hi[a] - lo[a]); let sb = 0.5 * (hi[b] - lo[b]);
    for (ga, wa) in GAUSS { for (gb, wb) in GAUSS {
        control.poll()?;
        let mut p = lo; p[a] = ma + sa * ga; p[b] = mb + sb * gb;
        let height = height_cut(sdf, cell, p, axis, increasing, best, control)?;
        let (bottom, top) = if increasing { (lo[axis], height) } else { (height, hi[axis]) };
        if top <= bottom { continue; }
        let center = f64::midpoint(bottom, top); let half = 0.5 * (top - bottom);
        for (g, w) in GAUSS {
            p[axis] = center + half * g;
            point(out, p, wa * wb * w * sa * sb * half, control)?;
        }
    } }
    Ok(())
}

// Enclose the CLAMPED height threshold, not a root-count assertion. Monotonicity
// proves an empty/full line as well as a crossing. When a point interval straddles
// zero, the certified minimum slope bounds the distance to the threshold.
fn height_cut(sdf: &dyn CutSdf3, cell: HexCell, mut p: [f64; 3], axis: usize,
    increasing: bool, slope: f64, control: &mut QuadratureControl3<'_>) -> Result<f64, QuadratureError3> {
    let (mut lower, mut upper) = (cell.lo()[axis], cell.hi()[axis]);
    let tolerance = (upper - lower) * control.options.root_relative_tolerance;
    if !tolerance.is_finite() || tolerance <= 0.0 { return Err(QuadratureError3::Invalid("root tolerance under/overflow")); }
    p[axis] = lower;
    let left = control.field(|| sdf.enclose(p, p))?;
    p[axis] = upper;
    let right = control.field(|| sdf.enclose(p, p))?;
    if (increasing && left.lo() >= 0.0) || (!increasing && left.hi() <= 0.0) {
        return Ok(lower);
    }
    if (increasing && right.hi() <= 0.0) || (!increasing && right.lo() >= 0.0) {
        return Ok(upper);
    }
    for _ in 0..control.options.root_iterations {
        control.poll()?;
        if upper - lower <= tolerance { return Ok(f64::midpoint(lower, upper)); }
        let mid = f64::midpoint(lower, upper);
        if mid <= lower || mid >= upper { break; }
        p[axis] = mid;
        let iv = control.field(|| sdf.enclose(p, p))?;
        if iv.lo() == 0.0 && iv.hi() == 0.0 { return Ok(mid); }
        if iv.hi() < 0.0 {
            if increasing { lower = mid; } else { upper = mid; }
        } else if iv.lo() > 0.0 {
            if increasing { upper = mid; } else { lower = mid; }
        } else {
            let radius = (iv.lo().abs().max(iv.hi().abs()).next_up() / slope).next_up();
            if !radius.is_finite() { break; }
            let new_lower = lower.max((mid - radius).next_down());
            let new_upper = upper.min((mid + radius).next_up());
            if new_lower > new_upper { return Err(QuadratureError3::Invalid("inconsistent field/derivative enclosures")); }
            if new_lower <= lower && new_upper >= upper { break; }
            lower = new_lower; upper = new_upper;
        }
    }
    if upper - lower <= tolerance { Ok(f64::midpoint(lower, upper)) }
    else { Err(QuadratureError3::RootResolution(cell)) }
}
