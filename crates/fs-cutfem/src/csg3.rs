//! Bounded constructive implicit domains for three-dimensional CutFEM.
//!
//! This is a geometry-to-physics adapter, not a distance oracle. Quadrics use
//! smooth, length-valued algebraic fields; boxes use a max-of-planes field.
//! Negative means material. Hard Booleans preserve the exact represented sets;
//! a positive blend width deliberately changes the boundary. Neither field
//! magnitudes nor the numerical quadrature are promoted to distance/error bounds.
//!
//! Field and partial-derivative ranges use outward-rounded interval arithmetic.
//! At hard creases the derivative range encloses all one-sided slopes. A strict
//! sign still proves coordinate-line monotonicity of the continuous, locally
//! Lipschitz field; an interval containing zero leaves the existing isolator to
//! subdivide/refuse. No smooth shape derivative is claimed at such a crease.

use crate::{CutSdf3, HeightAxis};
use fs_ivl::Interval;

/// Maximum work per field/derivative evaluation; shared children are evaluated once.
pub const MAX_CSG3_NODES: usize = 128;

/// A node in one [`CsgBuilder3`]. Handles must be used with their original builder.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CsgNode3(usize);

/// Constructive set operation. `Difference` means left minus right.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CsgOp3 {
    /// Union of the represented interiors.
    Union,
    /// Intersection of the represented interiors.
    Intersection,
    /// Remove the right interior from the left.
    Difference,
}

/// Invalid construction; no partially admitted domain is returned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CsgError3(pub &'static str);
impl core::fmt::Display for CsgError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(self.0)
    }
}
impl std::error::Error for CsgError3 {}

#[derive(Debug, Clone, Copy)]
enum Node {
    Plane { normal: [f64; 3], offset: f64 },
    Ellipsoid { center: [f64; 3], radii: [f64; 3] },
    Cylinder { axis: usize, center: [f64; 3], radius: f64 },
    Box { center: [f64; 3], half: [f64; 3] },
    Boolean { op: CsgOp3, a: usize, b: usize, blend: f64 },
}

/// A topologically ordered, allocation-bounded constructive-domain builder.
#[derive(Debug, Default)]
pub struct CsgBuilder3 {
    nodes: Vec<Node>,
}

fn finite(values: &[f64]) -> bool { values.iter().all(|v| v.is_finite()) }
fn lengths(values: &[f64]) -> bool {
    values.iter().all(|v| v.is_finite() && (1e-100..=1e100).contains(v))
}
fn axis_index(axis: HeightAxis) -> usize {
    match axis { HeightAxis::X => 0, HeightAxis::Y => 1, HeightAxis::Z => 2 }
}

impl CsgBuilder3 {
    /// Start an empty domain recipe.
    #[must_use]
    pub fn new() -> Self { Self::default() }

    fn push(&mut self, node: Node) -> Result<CsgNode3, CsgError3> {
        if self.nodes.len() == MAX_CSG3_NODES {
            return Err(CsgError3("constructive domain exceeds 128 nodes"));
        }
        let id = CsgNode3(self.nodes.len());
        self.nodes.push(node);
        Ok(id)
    }

    /// Half-space `normal dot point - offset < 0`. The finite nonzero normal is
    /// interpreted literally, NOT normalized or given exact-distance authority.
    pub fn half_space(&mut self, normal: [f64; 3], offset: f64) -> Result<CsgNode3, CsgError3> {
        if !finite(&normal) || !offset.is_finite() || normal.iter().all(|v| *v == 0.0) {
            return Err(CsgError3("half-space needs a finite nonzero normal and finite offset"));
        }
        self.push(Node::Plane { normal, offset })
    }

    /// Ellipsoid with finite center and positive semi-axes in `[1e-100,1e100]`.
    /// The field is `min(radii)/2 * (sum(((p-center)/radii)^2) - 1)`.
    pub fn ellipsoid(&mut self, center: [f64; 3], radii: [f64; 3]) -> Result<CsgNode3, CsgError3> {
        if !finite(&center) || !lengths(&radii) {
            return Err(CsgError3("ellipsoid needs a finite center and admitted positive semi-axes"));
        }
        self.push(Node::Ellipsoid { center, radii })
    }

    /// Sphere represented by the smooth algebraic field, not a distance sample.
    pub fn sphere(&mut self, center: [f64; 3], radius: f64) -> Result<CsgNode3, CsgError3> {
        self.ellipsoid(center, [radius; 3])
    }

    /// Infinite circular cylinder, clipped by the caller's background domain or
    /// another Boolean. The axis coordinate of `center` does not affect the set.
    pub fn cylinder(&mut self, axis: HeightAxis, center: [f64; 3], radius: f64) -> Result<CsgNode3, CsgError3> {
        if !finite(&center) || !lengths(&[radius]) {
            return Err(CsgError3("cylinder needs a finite center and admitted positive radius"));
        }
        self.push(Node::Cylinder { axis: axis_index(axis), center, radius })
    }

    /// Box represented by `max(abs(point-center)-half_extents)`.
    pub fn box_region(&mut self, center: [f64; 3], half_extents: [f64; 3]) -> Result<CsgNode3, CsgError3> {
        if !finite(&center) || !lengths(&half_extents) {
            return Err(CsgError3("box needs a finite center and admitted positive half extents"));
        }
        self.push(Node::Box { center, half: half_extents })
    }

    /// Combine two existing nodes. Zero width uses exact hard min/max; positive
    /// width uses a C1 polynomial blend in FIELD units, changing the geometry.
    /// A blend does not smooth pre-existing creases inside either child.
    pub fn combine(&mut self, op: CsgOp3, a: CsgNode3, b: CsgNode3, blend: f64) -> Result<CsgNode3, CsgError3> {
        if a.0 >= self.nodes.len() || b.0 >= self.nodes.len() {
            return Err(CsgError3("Boolean children must precede their parent in this builder"));
        }
        if blend != 0.0 && !lengths(&[blend]) {
            return Err(CsgError3("blend width must be zero or lie in [1e-100,1e100]"));
        }
        self.push(Node::Boolean { op, a: a.0, b: b.0, blend })
    }

    /// Seal the immutable recipe with an explicitly selected root.
    pub fn finish(self, root: CsgNode3) -> Result<CsgDomain3, CsgError3> {
        if root.0 >= self.nodes.len() { return Err(CsgError3("missing constructive root")); }
        Ok(CsgDomain3 { nodes: self.nodes, root: root.0 })
    }
}

/// An immutable constructive field accepted by all existing [`CutSdf3`] consumers.
#[derive(Debug, Clone)]
pub struct CsgDomain3 {
    nodes: Vec<Node>,
    root: usize,
}

fn point(x: f64) -> Interval { Interval::new(x, x) }
fn whole() -> Interval { Interval::new(f64::NEG_INFINITY, f64::INFINITY) }
fn square(x: Interval) -> Interval {
    let product = x * x;
    if x.lo() <= 0.0 && x.hi() >= 0.0 { Interval::new(0.0, product.hi()) }
    else { product }
}
fn absolute(x: Interval) -> Interval {
    if x.lo() >= 0.0 { x }
    else if x.hi() <= 0.0 { point(0.0) - x }
    else { Interval::new(0.0, (-x.lo()).max(x.hi())) }
}
fn hull(a: Interval, b: Interval) -> Interval {
    Interval::new(a.lo().min(b.lo()), a.hi().max(b.hi()))
}
fn clamp_unit(x: Interval) -> Interval {
    Interval::new(x.lo().clamp(0.0, 1.0), x.hi().clamp(0.0, 1.0))
}

#[derive(Clone, Copy)]
struct Range {
    value: Interval,
    derivative: Interval,
}
impl Range {
    fn negate(self) -> Self {
        Self { value: point(0.0) - self.value, derivative: point(0.0) - self.derivative }
    }
    fn checked(self) -> Self {
        if [self.value.lo(), self.value.hi(), self.derivative.lo(), self.derivative.hi()]
            .iter().all(|v| v.is_finite()) { self }
        else { Self { value: whole(), derivative: whole() } }
    }
}
fn minimum(a: Range, b: Range, blend: f64) -> Range {
    // Do not let f64::min/max hide an invalid/overflowed child.
    if ![a.value.lo(), a.value.hi(), b.value.lo(), b.value.hi(),
        a.derivative.lo(), a.derivative.hi(), b.derivative.lo(), b.derivative.hi()]
        .iter().all(|v| v.is_finite()) {
        return Range { value: whole(), derivative: whole() };
    }
    if blend == 0.0 {
        if a.value.hi() < b.value.lo() { return a; }
        if b.value.hi() < a.value.lo() { return b; }
        return Range {
            value: Interval::new(a.value.lo().min(b.value.lo()), a.value.hi().min(b.value.hi())),
            derivative: hull(a.derivative, b.derivative),
        };
    }
    let k = point(blend);
    let gap = b.value - a.value;
    if gap.lo() >= blend { return a; }
    if gap.hi() <= -blend { return b; }
    let h = clamp_unit((k - absolute(gap)) / k);
    let t = clamp_unit(point(0.5) + gap / (point(2.0) * k));
    Range {
        value: Interval::new(a.value.lo().min(b.value.lo()), a.value.hi().min(b.value.hi()))
            - point(0.25) * k * square(h),
        // The derivatives of t cancel algebraically in the polynomial blend.
        derivative: t * a.derivative + (point(1.0) - t) * b.derivative,
    }.checked()
}
fn minimum_value(a: f64, b: f64, blend: f64) -> f64 {
    if !a.is_finite() || !b.is_finite() { return f64::NAN; }
    if blend == 0.0 { return a.min(b); }
    let h = ((blend - (a - b).abs()) / blend).max(0.0);
    a.min(b) - 0.25 * blend * h * h
}

impl Node {
    fn range(self, lo: [f64; 3], hi: [f64; 3], axis: Option<usize>, values: &[Range]) -> Range {
        let x: [Interval; 3] = std::array::from_fn(|i| Interval::new(lo[i], hi[i]));
        match self {
            Self::Plane { normal, offset } => Range {
                value: (0..3).fold(point(0.0), |sum, i| sum + point(normal[i]) * x[i]) - point(offset),
                derivative: point(axis.map_or(0.0, |i| normal[i])),
            },
            Self::Ellipsoid { center, radii } => {
                let q: [Interval; 3] = std::array::from_fn(|i| (x[i] - point(center[i])) / point(radii[i]));
                let scale = radii.into_iter().fold(f64::INFINITY, f64::min);
                Range {
                    value: point(0.5 * scale) * (q.into_iter().fold(point(0.0), |sum, v| sum + square(v)) - point(1.0)),
                    derivative: axis.map_or(point(0.0), |i| point(scale) * q[i] / point(radii[i])),
                }
            }
            Self::Cylinder { axis: along, center, radius } => {
                let q: [Interval; 3] = std::array::from_fn(|i| {
                    if i == along { point(0.0) } else { (x[i] - point(center[i])) / point(radius) }
                });
                Range {
                    value: point(0.5 * radius) * (q.into_iter().fold(point(0.0), |sum, v| sum + square(v)) - point(1.0)),
                    derivative: axis.map_or(point(0.0), |i| q[i]),
                }
            }
            Self::Box { center, half } => {
                let face = |i: usize| {
                    let d = x[i] - point(center[i]);
                    let derivative = if axis != Some(i) { point(0.0) }
                        else if d.lo() > 0.0 { point(1.0) }
                        else if d.hi() < 0.0 { point(-1.0) }
                        else { Interval::new(-1.0, 1.0) };
                    Range { value: absolute(d) - point(half[i]), derivative }
                };
                let mut out = face(0);
                for i in 1..3 { out = minimum(out.negate(), face(i).negate(), 0.0).negate(); }
                out
            }
            Self::Boolean { op, a, b, blend } => match op {
                CsgOp3::Union => minimum(values[a], values[b], blend),
                CsgOp3::Intersection => minimum(values[a].negate(), values[b].negate(), blend).negate(),
                CsgOp3::Difference => minimum(values[a].negate(), values[b], blend).negate(),
            },
        }.checked()
    }

    fn value(self, p: [f64; 3], values: &[f64]) -> f64 {
        match self {
            Self::Plane { normal, offset } => (0..3).fold(0.0, |sum, i| sum + normal[i] * p[i]) - offset,
            Self::Ellipsoid { center, radii } => {
                let sum = (0..3).fold(0.0, |sum, i| {
                    let q = (p[i] - center[i]) / radii[i]; sum + q * q
                });
                0.5 * radii.into_iter().fold(f64::INFINITY, f64::min) * (sum - 1.0)
            }
            Self::Cylinder { axis, center, radius } => {
                let sum = (0..3).filter(|i| *i != axis).fold(0.0, |sum, i| {
                    let q = (p[i] - center[i]) / radius; sum + q * q
                });
                0.5 * radius * (sum - 1.0)
            }
            Self::Box { center, half } => (0..3).map(|i| (p[i] - center[i]).abs() - half[i])
                .fold(f64::NEG_INFINITY, f64::max),
            Self::Boolean { op, a, b, blend } => match op {
                CsgOp3::Union => minimum_value(values[a], values[b], blend),
                CsgOp3::Intersection => -minimum_value(-values[a], -values[b], blend),
                CsgOp3::Difference => -minimum_value(-values[a], values[b], blend),
            },
        }
    }
}

impl CsgDomain3 {
    /// Retained bounded recipe size, including explicitly shared nodes.
    #[must_use]
    pub fn node_count(&self) -> usize { self.nodes.len() }

    fn range(&self, lo: [f64; 3], hi: [f64; 3], axis: Option<usize>) -> Range {
        if !(0..3).all(|i| lo[i].is_finite() && hi[i].is_finite() && lo[i] <= hi[i]) {
            return Range { value: whole(), derivative: whole() };
        }
        let mut values = [Range { value: point(0.0), derivative: point(0.0) }; MAX_CSG3_NODES];
        for (i, node) in self.nodes.iter().enumerate().take(self.root + 1) {
            values[i] = node.range(lo, hi, axis, &values);
        }
        values[self.root]
    }
}
impl CutSdf3 for CsgDomain3 {
    fn value(&self, p: [f64; 3]) -> f64 {
        if !finite(&p) { return f64::NAN; }
        let mut values = [0.0; MAX_CSG3_NODES];
        for (i, node) in self.nodes.iter().enumerate().take(self.root + 1) {
            values[i] = node.value(p, &values);
        }
        values[self.root]
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval { self.range(lo, hi, None).value }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], axis: HeightAxis) -> Interval {
        self.range(lo, hi, Some(axis_index(axis))).derivative
    }
}

#[cfg(test)]
#[path = "csg3_tests.rs"]
mod tests;
