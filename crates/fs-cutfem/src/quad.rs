//! Cut-cell quadrature: accurate integration on implicitly-defined
//! subdomains, with error control by subdivision depth.
//!
//! The scheme is tessellation-based with CERTIFIED routing: a cut cell
//! recursively subdivides; every sub-cell is re-classified through the
//! SDF's interval enclosure, so fully-inside sub-cells get exact tensor
//! Gauss rules and fully-outside sub-cells vanish — only genuinely cut
//! sub-cells reach the marching-squares finale, where edge crossings
//! are located by bisection (resolution ~2⁻⁵⁰ of the sub-cell) and the
//! inside region is polygonized and integrated with a degree-2-exact
//! midpoint triangle rule.
//!
//! Error control: for a linear level set the crossings, the polygon,
//! and hence ALL quadratic moments are exact (the G0 exactness
//! battery). For curved interfaces the geometric error is
//! O((h/2^depth)²) per cell — quadratic convergence in depth,
//! measured and asserted by the conformance suite. Features smaller
//! than the finest subdivision that never change a corner sign can be
//! missed by the QUADRATURE (bounded by one sub-cell area) — never by
//! the CLASSIFICATION, which is interval-certified at every level; the
//! saddle/blob ambiguity spends a bounded extra-subdivision budget
//! before falling back to corner-sign integration.

use crate::sdf::CutSdf;
use fs_ivl::Interval;

/// Bulk and interface rules for one cut cell, in global coordinates.
/// Bulk weights sum to the inside area (up to the documented error);
/// interface weights sum to the interface length, each point carrying
/// the OUTWARD unit normal (∇φ direction).
#[derive(Debug, Clone, Default)]
pub struct CutRules {
    /// Bulk points: (position, weight).
    pub bulk: Vec<([f64; 2], f64)>,
    /// Interface points: (position, weight, outward unit normal).
    pub iface: Vec<([f64; 2], f64, [f64; 2])>,
}

/// A certified implicit field over an axis-aligned 3-D cell.
///
/// `enclose` must contain every field value in its box, and
/// `derivative_enclose` must contain every partial derivative in its box.
/// The vertical-line isolator relies on continuity plus a derivative interval
/// excluding zero to prove that a line has at most one crossing; it refuses
/// instead of inferring a root from sampled values.
pub trait CutSdf3 {
    /// Field value at a point; negative values are inside the domain.
    fn value(&self, point: [f64; 3]) -> f64;
    /// Certified enclosure of the field over `[lo, hi]`.
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval;
    /// Certified enclosure of the partial derivative along `axis`.
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], axis: HeightAxis) -> Interval;
}

/// A finite, nondegenerate axis-aligned hexahedral cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HexCell {
    /// Lower coordinates.
    lo: [f64; 3],
    /// Upper coordinates.
    hi: [f64; 3],
}

impl HexCell {
    /// Construct a cell only when every axis is finite and has positive span.
    pub fn try_new(lo: [f64; 3], hi: [f64; 3]) -> Result<Self, CutQuadrature3dError> {
        for axis in 0..3 {
            if !(lo[axis].is_finite() && hi[axis].is_finite()) {
                return Err(CutQuadrature3dError::NonFiniteCell { axis });
            }
            if lo[axis] >= hi[axis] {
                return Err(CutQuadrature3dError::NonPositiveCellSpan { axis });
            }
        }
        Ok(Self { lo, hi })
    }

    /// Validated lower coordinates.
    #[must_use]
    pub const fn lo(self) -> [f64; 3] {
        self.lo
    }

    /// Validated upper coordinates.
    #[must_use]
    pub const fn hi(self) -> [f64; 3] {
        self.hi
    }
}

/// Candidate height directions for Saye-style dimension reduction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HeightAxis {
    /// Use x as the vertical coordinate.
    X,
    /// Use y as the vertical coordinate.
    Y,
    /// Use z as the vertical coordinate.
    Z,
}

impl HeightAxis {
    const ALL: [Self; 3] = [Self::X, Self::Y, Self::Z];

    const fn index(self) -> usize {
        match self {
            Self::X => 0,
            Self::Y => 1,
            Self::Z => 2,
        }
    }

    const fn base_axes(self) -> [usize; 2] {
        match self {
            Self::X => [1, 2],
            Self::Y => [0, 2],
            Self::Z => [0, 1],
        }
    }
}

/// Input refusals for the certified 3-D vertical-line primitive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CutQuadrature3dError {
    /// A cell coordinate was not finite.
    NonFiniteCell { axis: usize },
    /// A cell coordinate did not have strictly positive span.
    NonPositiveCellSpan { axis: usize },
    /// A base-line coordinate was not finite.
    NonFiniteBase { coordinate: usize },
    /// A base-line coordinate lies outside the selected cell.
    BaseOutsideCell { axis: usize },
}

impl core::fmt::Display for CutQuadrature3dError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonFiniteCell { axis } => write!(f, "3-D cut cell axis {axis} is non-finite"),
            Self::NonPositiveCellSpan { axis } => {
                write!(f, "3-D cut cell axis {axis} has non-positive span")
            }
            Self::NonFiniteBase { coordinate } => {
                write!(f, "3-D cut-line base coordinate {coordinate} is non-finite")
            }
            Self::BaseOutsideCell { axis } => {
                write!(f, "3-D cut-line base lies outside cell axis {axis}")
            }
        }
    }
}

impl std::error::Error for CutQuadrature3dError {}

/// Fail-closed outcome of one certified height line.
#[derive(Debug, Clone, PartialEq)]
pub enum VerticalLineRoot {
    /// The derivative certificate proves the line cannot cross zero.
    NoIntersection { height_axis: HeightAxis },
    /// Exactly one root is proven and retained as a closed interval.
    CertifiedRoot {
        /// Selected height coordinate.
        height_axis: HeightAxis,
        /// Certified enclosure of the unique crossing.
        enclosure: Interval,
    },
    /// No root count is claimed; the caller must subdivide or refuse.
    Ambiguous,
}

/// Isolate the unique crossing on one Saye-style height line when possible.
///
/// The height axis is selected deterministically from derivative enclosures:
/// the largest certified lower bound on `|∂φ|` wins, with X/Y/Z as the tie
/// order. `base` is ordered by the two remaining coordinates (Y,Z for X;
/// X,Z for Y; X,Y for Z). If no finite derivative enclosure excludes zero, or
/// any sampled line point cannot be sign-certified, this returns `Ambiguous`.
/// In particular, an endpoint enclosure containing zero is not proof that the
/// endpoint is a root: it could contain either strict sign.
/// It therefore supplies a certified root primitive to a future 3-D tensor
/// quadrature driver without claiming volume/surface quadrature by itself.
pub fn isolate_certified_height_root(
    sdf: &dyn CutSdf3,
    cell: HexCell,
    base: [f64; 2],
    max_bisections: u32,
) -> Result<VerticalLineRoot, CutQuadrature3dError> {
    let Some(height_axis) = select_height_axis(sdf, cell) else {
        return Ok(VerticalLineRoot::Ambiguous);
    };
    let base_axes = height_axis.base_axes();
    for coordinate in 0..2 {
        let axis = base_axes[coordinate];
        if !base[coordinate].is_finite() {
            return Err(CutQuadrature3dError::NonFiniteBase { coordinate });
        }
        if base[coordinate] < cell.lo[axis] || base[coordinate] > cell.hi[axis] {
            return Err(CutQuadrature3dError::BaseOutsideCell { axis });
        }
    }

    let axis = height_axis.index();
    let derivative = sdf.derivative_enclose(cell.lo, cell.hi, height_axis);
    if !is_finite_interval(derivative) {
        return Ok(VerticalLineRoot::Ambiguous);
    }
    let increasing = derivative.lo() > 0.0;
    let mut lower = cell.lo[axis];
    let mut upper = cell.hi[axis];
    let lower_sign = sign_of(sdf.enclose(
        line_point(height_axis, base, lower),
        line_point(height_axis, base, lower),
    ));
    let upper_sign = sign_of(sdf.enclose(
        line_point(height_axis, base, upper),
        line_point(height_axis, base, upper),
    ));

    let has_crossing = match (increasing, lower_sign, upper_sign) {
        (true, Some(Sign::Positive), _) | (true, _, Some(Sign::Negative)) => {
            return Ok(VerticalLineRoot::NoIntersection { height_axis });
        }
        (false, Some(Sign::Negative), _) | (false, _, Some(Sign::Positive)) => {
            return Ok(VerticalLineRoot::NoIntersection { height_axis });
        }
        (true, Some(Sign::Negative), Some(Sign::Positive))
        | (false, Some(Sign::Positive), Some(Sign::Negative)) => true,
        _ => false,
    };
    if !has_crossing {
        return Ok(VerticalLineRoot::Ambiguous);
    }

    let lower_is_negative = matches!(lower_sign, Some(Sign::Negative));
    for _ in 0..max_bisections {
        let midpoint = f64::midpoint(lower, upper);
        if midpoint == lower || midpoint == upper {
            break;
        }
        match sign_of(sdf.enclose(
            line_point(height_axis, base, midpoint),
            line_point(height_axis, base, midpoint),
        )) {
            Some(Sign::Negative) if lower_is_negative => lower = midpoint,
            Some(Sign::Positive) if lower_is_negative => upper = midpoint,
            Some(Sign::Positive) => lower = midpoint,
            Some(Sign::Negative) => upper = midpoint,
            None => return Ok(VerticalLineRoot::Ambiguous),
        }
    }
    Ok(VerticalLineRoot::CertifiedRoot {
        height_axis,
        enclosure: Interval::new(lower, upper),
    })
}

#[derive(Clone, Copy)]
enum Sign {
    Negative,
    Positive,
}

fn sign_of(interval: Interval) -> Option<Sign> {
    if !is_finite_interval(interval) {
        None
    } else if interval.hi() < 0.0 {
        Some(Sign::Negative)
    } else if interval.lo() > 0.0 {
        Some(Sign::Positive)
    } else {
        None
    }
}

fn is_finite_interval(interval: Interval) -> bool {
    interval.lo().is_finite() && interval.hi().is_finite()
}

fn select_height_axis(sdf: &dyn CutSdf3, cell: HexCell) -> Option<HeightAxis> {
    let mut selected = None;
    let mut lower_bound = 0.0;
    for axis in HeightAxis::ALL {
        let derivative = sdf.derivative_enclose(cell.lo, cell.hi, axis);
        if !is_finite_interval(derivative) {
            // An unusable derivative certificate rules out only this height
            // direction. Another axis can still prove strict monotonicity.
            continue;
        }
        let candidate = if derivative.lo() > 0.0 {
            derivative.lo()
        } else if derivative.hi() < 0.0 {
            -derivative.hi()
        } else {
            0.0
        };
        if candidate > lower_bound {
            selected = Some(axis);
            lower_bound = candidate;
        }
    }
    selected
}

fn line_point(height_axis: HeightAxis, base: [f64; 2], height: f64) -> [f64; 3] {
    match height_axis {
        HeightAxis::X => [height, base[0], base[1]],
        HeightAxis::Y => [base[0], height, base[1]],
        HeightAxis::Z => [base[0], base[1], height],
    }
}

/// 3-point Gauss–Legendre on [-1, 1].
const G3: [(f64, f64); 3] = [
    (-0.774_596_669_241_483_4, 0.555_555_555_555_555_6),
    (0.0, 0.888_888_888_888_889),
    (0.774_596_669_241_483_4, 0.555_555_555_555_555_6),
];

/// Extra subdivision budget for saddle/blob ambiguity at depth 0.
const EXTRA: u32 = 2;

/// Push a 3×3 tensor Gauss rule for the full box (degree-5 exact per
/// axis — exact for every integrand this crate assembles).
pub fn tensor_gauss(lo: [f64; 2], hi: [f64; 2], out: &mut Vec<([f64; 2], f64)>) {
    let mx = f64::midpoint(lo[0], hi[0]);
    let my = f64::midpoint(lo[1], hi[1]);
    let sx = 0.5 * (hi[0] - lo[0]);
    let sy = 0.5 * (hi[1] - lo[1]);
    for &(gx, wx) in &G3 {
        for &(gy, wy) in &G3 {
            out.push(([mx + sx * gx, my + sy * gy], wx * wy * sx * sy));
        }
    }
}

/// Build the bulk + interface rules for one (certified-Cut) cell.
#[must_use]
pub fn cut_cell_rules(sdf: &dyn CutSdf, lo: [f64; 2], hi: [f64; 2], depth: u32) -> CutRules {
    let mut rules = CutRules::default();
    worker(sdf, lo, hi, depth, EXTRA, &mut rules);
    rules
}

fn worker(
    sdf: &dyn CutSdf,
    lo: [f64; 2],
    hi: [f64; 2],
    depth: u32,
    extra: u32,
    out: &mut CutRules,
) {
    let iv = sdf.enclose(lo, hi);
    if iv.hi() < 0.0 {
        tensor_gauss(lo, hi, &mut out.bulk);
        return;
    }
    if iv.lo() > 0.0 {
        return;
    }
    if depth > 0 {
        recurse(sdf, lo, hi, depth - 1, extra, out);
        return;
    }
    // Finest level: marching squares on the corner signs.
    let corners = [
        [lo[0], lo[1]],
        [hi[0], lo[1]],
        [hi[0], hi[1]],
        [lo[0], hi[1]],
    ];
    let phi: Vec<f64> = corners.iter().map(|&p| sdf.value(p)).collect();
    let inside: Vec<bool> = phi.iter().map(|&v| v <= 0.0).collect();
    let mut crossings: [Option<[f64; 2]>; 4] = [None; 4];
    let mut ncross = 0;
    for e in 0..4 {
        if inside[e] != inside[(e + 1) % 4] {
            crossings[e] = Some(bisect_crossing(sdf, corners[e], corners[(e + 1) % 4]));
            ncross += 1;
        }
    }
    if ncross == 0 {
        // Interval says "maybe cut" but no corner sign change: a
        // sub-resolution feature (tangency or interior blob). Spend the
        // extra budget, then fall back to the corner-sign picture.
        if extra > 0 {
            recurse(sdf, lo, hi, 0, extra - 1, out);
        } else if inside[0] {
            tensor_gauss(lo, hi, &mut out.bulk);
        }
        return;
    }
    if ncross == 4 {
        // Saddle: prefer more resolution; if exhausted, resolve the
        // connectivity by the center sign.
        if extra > 0 {
            recurse(sdf, lo, hi, 0, extra - 1, out);
            return;
        }
        saddle_rules(sdf, &corners, &inside, &crossings, out);
        return;
    }
    // Regular case (2 crossings): walk the boundary counterclockwise,
    // emitting inside corners and crossings — a simple polygon.
    let mut poly: Vec<[f64; 2]> = Vec::with_capacity(6);
    for e in 0..4 {
        if inside[e] {
            poly.push(corners[e]);
        }
        if let Some(x) = crossings[e] {
            poly.push(x);
        }
    }
    polygon_rule(&poly, &mut out.bulk);
    let xs: Vec<[f64; 2]> = crossings.iter().flatten().copied().collect();
    chord_rule(sdf, xs[0], xs[1], &mut out.iface);
}

fn recurse(
    sdf: &dyn CutSdf,
    lo: [f64; 2],
    hi: [f64; 2],
    depth: u32,
    extra: u32,
    out: &mut CutRules,
) {
    let mx = f64::midpoint(lo[0], hi[0]);
    let my = f64::midpoint(lo[1], hi[1]);
    worker(sdf, lo, [mx, my], depth, extra, out);
    worker(sdf, [mx, lo[1]], [hi[0], my], depth, extra, out);
    worker(sdf, [lo[0], my], [mx, hi[1]], depth, extra, out);
    worker(sdf, [mx, my], hi, depth, extra, out);
}

/// Saddle finale: the center sign decides which crossing pairs bound
/// the inside region. Center inside → the hexagon walk (both chords
/// cut off outside corner lobes); center outside → two inside corner
/// triangles.
fn saddle_rules(
    sdf: &dyn CutSdf,
    corners: &[[f64; 2]; 4],
    inside: &[bool],
    crossings: &[Option<[f64; 2]>; 4],
    out: &mut CutRules,
) {
    let center = [
        f64::midpoint(corners[0][0], corners[2][0]),
        f64::midpoint(corners[0][1], corners[2][1]),
    ];
    let center_in = sdf.value(center) <= 0.0;
    if center_in {
        let mut poly: Vec<[f64; 2]> = Vec::with_capacity(6);
        for e in 0..4 {
            if inside[e] {
                poly.push(corners[e]);
            }
            if let Some(x) = crossings[e] {
                poly.push(x);
            }
        }
        polygon_rule(&poly, &mut out.bulk);
        for k in 0..4 {
            if !inside[k] {
                let a = crossings[(k + 3) % 4].expect("saddle has all crossings");
                let b = crossings[k].expect("saddle has all crossings");
                chord_rule(sdf, a, b, &mut out.iface);
            }
        }
    } else {
        for k in 0..4 {
            if inside[k] {
                let a = crossings[(k + 3) % 4].expect("saddle has all crossings");
                let b = crossings[k].expect("saddle has all crossings");
                polygon_rule(&[a, corners[k], b], &mut out.bulk);
                chord_rule(sdf, a, b, &mut out.iface);
            }
        }
    }
}

/// Bisection for the interface crossing on a segment with a corner
/// sign change (~2⁻⁵⁰ of the segment; exact for linear φ up to
/// roundoff).
fn bisect_crossing(sdf: &dyn CutSdf, a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    let sa = sdf.value(a) > 0.0;
    let (mut t0, mut t1) = (0.0f64, 1.0f64);
    for _ in 0..50 {
        let tm = f64::midpoint(t0, t1);
        let p = [a[0] + tm * (b[0] - a[0]), a[1] + tm * (b[1] - a[1])];
        if (sdf.value(p) > 0.0) == sa {
            t0 = tm;
        } else {
            t1 = tm;
        }
    }
    let t = f64::midpoint(t0, t1);
    [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]
}

/// Fan-triangulate a simple polygon (walk order) and push the
/// degree-2-exact midpoint rule per triangle. Signed areas make the
/// fan correct for non-convex walk polygons.
fn polygon_rule(poly: &[[f64; 2]], out: &mut Vec<([f64; 2], f64)>) {
    if poly.len() < 3 {
        return;
    }
    for k in 1..poly.len() - 1 {
        let (p, q, r) = (poly[0], poly[k], poly[k + 1]);
        let signed_area = 0.5 * ((q[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (q[1] - p[1]));
        if signed_area == 0.0 {
            continue;
        }
        let w = signed_area / 3.0;
        out.push(([f64::midpoint(p[0], q[0]), f64::midpoint(p[1], q[1])], w));
        out.push(([f64::midpoint(q[0], r[0]), f64::midpoint(q[1], r[1])], w));
        out.push(([f64::midpoint(r[0], p[0]), f64::midpoint(r[1], p[1])], w));
    }
}

/// 2-point Gauss along an interface chord; weights carry the chord
/// length, normals come from the (normalized) SDF gradient — outward
/// of Ω by the negative-inside convention.
fn chord_rule(
    sdf: &dyn CutSdf,
    a: [f64; 2],
    b: [f64; 2],
    out: &mut Vec<([f64; 2], f64, [f64; 2])>,
) {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let len = (dx * dx + dy * dy).sqrt();
    if len == 0.0 {
        return;
    }
    let g = 0.5 / 3.0f64.sqrt();
    for t in [0.5 - g, 0.5 + g] {
        let p = [a[0] + t * dx, a[1] + t * dy];
        let grad = sdf.gradient(p);
        let gn = (grad[0] * grad[0] + grad[1] * grad[1]).sqrt();
        let normal = if gn > 1e-300 {
            [grad[0] / gn, grad[1] / gn]
        } else {
            // Degenerate gradient: fall back to the chord perpendicular,
            // signed by a φ probe.
            let perp = [dy / len, -dx / len];
            let eps = 1e-9 * len.max(1e-12);
            let plus = sdf.value([p[0] + eps * perp[0], p[1] + eps * perp[1]]);
            let minus = sdf.value([p[0] - eps * perp[0], p[1] - eps * perp[1]]);
            if plus >= minus {
                perp
            } else {
                [-perp[0], -perp[1]]
            }
        };
        out.push((p, 0.5 * len, normal));
    }
}
