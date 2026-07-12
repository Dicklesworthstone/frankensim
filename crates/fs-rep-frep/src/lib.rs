//! fs-rep-frep — function-representation charts (plan §7.2). Layer: L2.
//!
//! A CSG **DAG** over implicit primitives, evaluated as one scalar field:
//! the function IS the shape (the region is `{ p : f(p) < 0 }`). Three
//! evaluators are AUTO-DERIVED from the same DAG — value+gradient (exact
//! chain rule), Lipschitz bound (per-node composition rules), and interval
//! range over a box (per-node inclusion rules) — the trio that feeds
//! certified sphere tracing (§10.2), certified inside/outside, and shape
//! optimization.
//!
//! Booleans come in two labeled flavors:
//! - [`BoolStyle::Hard`] — `min`/`max`. Exact set operations, but the
//!   derivative is DISCONTINUOUS across the `a = b` crease, which POISONS
//!   gradient-based shape optimization (the plan calls this out). Kept for
//!   non-optimization paths.
//! - [`BoolStyle::Blend`] — the quadratic-polynomial R-function blend
//!   (smooth min with radius `r`): C¹ everywhere, weights sum to one, and
//!   the radius is a DESIGN LEVER (`(fillet :r 3mm)` is a blend radius).
//!
//! Honesty of the field: valid sphere/half-space/cylinder/box primitives and
//! ring tori carry exact-distance semantics; horn/spindle tori, offsets, and
//! Boolean composites explicitly downgrade magnitude exactness. Every valid
//! node preserves its recorded Lipschitz composition rule and vanishes on its
//! own region's boundary — hence
//! `|f(p)| ≤ dist(p, ∂Ω)` with the EXACT sign. That one-sided conservative
//! bound is precisely the sphere-tracing safety contract; the magnitude
//! is NOT the exact distance once a Boolean or erosion is involved, so
//! composite samples carry an `Estimate` certificate (see CONTRACT.md
//! no-claims). Only primitive paths proven by `is_exact` stay `Exact`.

mod ival;

use fs_evidence::NumericalCertificate;
use fs_exec::Cx;
use fs_geom::{
    Aabb, BettiBounds, Chart, ChartSample, Differentiability, Point3, TraceStepClaim, Vec3,
};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Half-extent standing in for an unbounded support axis (half-spaces,
/// infinite cylinders). Intersections shrink it back to honest bounds.
pub const UNBOUNDED_HALF: f64 = 1.0e12;

/// A node handle inside one [`Frep`] DAG (indices are a topological
/// order: children always precede parents).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct NodeId(pub u32);

/// Which set operation a Boolean node performs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoolOp {
    /// `a ∪ b`.
    Union,
    /// `a ∩ b`.
    Intersect,
    /// `a \ b`.
    Difference,
}

/// How the Boolean is realized (see module docs for the trade).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoolStyle {
    /// Exact `min`/`max` — derivative discontinuity at `a = b`.
    Hard,
    /// C¹ quadratic R-function blend; `radius > 0` is a design lever.
    Blend {
        /// Blend radius (world units).
        radius: f64,
    },
}

/// One DAG node. Numeric fields are addressable as design parameters
/// (see [`Frep::params`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Node {
    /// Exact sphere SDF `|p − c| − r`.
    Sphere {
        /// Center.
        center: Point3,
        /// Radius (> 0).
        radius: f64,
    },
    /// Half-space `n·p − offset` (unit `n`; inside where negative).
    HalfSpace {
        /// Unit outward normal.
        normal: Vec3,
        /// Plane offset along `normal`.
        offset: f64,
    },
    /// Exact axis-aligned box SDF (corner-distance formula).
    BoxPrim {
        /// Center.
        center: Point3,
        /// Half extents (all > 0).
        half: Vec3,
    },
    /// Torus implicit field, axis +z through `center`. It is an exact signed
    /// distance for ring tori (`major > minor`); horn/spindle tori retain the
    /// same exact sign/zero set but only the Lipschitz-implicit trace claim.
    Torus {
        /// Center.
        center: Point3,
        /// Major (ring) radius (> 0).
        major: f64,
        /// Minor (tube) radius (> 0).
        minor: f64,
    },
    /// Infinite cylinder along +z through `center` (exact SDF).
    Cylinder {
        /// A point on the axis.
        center: Point3,
        /// Radius (> 0).
        radius: f64,
    },
    /// Rigid translation of a child region.
    Translate {
        /// Child node.
        child: NodeId,
        /// Translation vector.
        offset: Vec3,
    },
    /// Rigid rotation about `axis` (unit) through the origin by `angle`
    /// radians. (GA motors join with fs-ga; axis-angle is the v1 lever.)
    Rotate {
        /// Child node.
        child: NodeId,
        /// Unit rotation axis.
        axis: Vec3,
        /// Angle in radians.
        angle: f64,
    },
    /// Uniform scale by `factor > 0` (SDF-preserving: `s·f(p/s)`).
    Scale {
        /// Child node.
        child: NodeId,
        /// Scale factor (> 0).
        factor: f64,
    },
    /// Offset surface `f − distance`. Its sign/zero set is exact, but without
    /// a reach certificate neither dilation nor erosion is promoted to an
    /// exact signed distance.
    Offset {
        /// Child node.
        child: NodeId,
        /// Signed offset distance.
        distance: f64,
    },
    /// Boolean combination of two children.
    Bool {
        /// Operation.
        op: BoolOp,
        /// Hard or blended.
        style: BoolStyle,
        /// Left operand.
        a: NodeId,
        /// Right operand.
        b: NodeId,
    },
}

/// Teaching error for malformed DAG construction or parameter edits.
#[derive(Debug, Clone, PartialEq)]
pub enum FrepError {
    /// A length/radius/scale that must be strictly positive was not.
    NonPositive {
        /// Which field.
        field: &'static str,
        /// The offending value.
        value: f64,
    },
    /// A direction vector was too short to normalize.
    ZeroVector {
        /// Which field.
        field: &'static str,
    },
    /// A referenced child does not exist (or the DAG is empty).
    BadNode {
        /// The offending id.
        id: u32,
    },
    /// `set_param` addressed a slot the node does not have.
    BadParam {
        /// Node id.
        node: u32,
        /// Slot index.
        slot: u8,
    },
}

impl core::fmt::Display for FrepError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            FrepError::NonPositive { field, value } => write!(
                f,
                "`{field}` must be strictly positive, got {value}; zero-thickness \
                 primitives have no interior and break the SDF sign convention"
            ),
            FrepError::ZeroVector { field } => write!(
                f,
                "`{field}` is too short to normalize; pass a direction with \
                 norm well above 1e-12"
            ),
            FrepError::BadNode { id } => write!(
                f,
                "node id {id} does not exist in this DAG; combine only ids \
                 returned by this builder"
            ),
            FrepError::BadParam { node, slot } => write!(
                f,
                "node {node} has no parameter slot {slot}; enumerate levers \
                 with Frep::params()"
            ),
        }
    }
}

impl std::error::Error for FrepError {}

/// Builder for a [`Frep`] DAG. Methods validate inputs (teaching errors)
/// and return [`NodeId`]s that later nodes may reference — sharing a
/// subexpression is just reusing its id.
#[derive(Debug, Default)]
pub struct FrepBuilder {
    nodes: Vec<Node>,
}

impl FrepBuilder {
    /// Empty builder.
    #[must_use]
    pub fn new() -> Self {
        FrepBuilder::default()
    }

    fn push(&mut self, node: Node) -> NodeId {
        self.nodes.push(node);
        NodeId((self.nodes.len() - 1) as u32)
    }

    fn check_child(&self, id: NodeId) -> Result<(), FrepError> {
        if (id.0 as usize) < self.nodes.len() {
            Ok(())
        } else {
            Err(FrepError::BadNode { id: id.0 })
        }
    }

    fn positive(field: &'static str, value: f64) -> Result<f64, FrepError> {
        if value > 0.0 && value.is_finite() {
            Ok(value)
        } else {
            Err(FrepError::NonPositive { field, value })
        }
    }

    fn unit(field: &'static str, v: Vec3) -> Result<Vec3, FrepError> {
        let n = v.norm();
        if n < 1e-12 {
            return Err(FrepError::ZeroVector { field });
        }
        Ok(v.scale(1.0 / n))
    }

    /// Sphere primitive.
    ///
    /// # Errors
    /// [`FrepError::NonPositive`] for `radius ≤ 0`.
    pub fn sphere(&mut self, center: Point3, radius: f64) -> Result<NodeId, FrepError> {
        let radius = Self::positive("radius", radius)?;
        Ok(self.push(Node::Sphere { center, radius }))
    }

    /// Half-space primitive (`normal` is normalized here).
    ///
    /// # Errors
    /// [`FrepError::ZeroVector`] for a degenerate normal.
    pub fn half_space(&mut self, normal: Vec3, offset: f64) -> Result<NodeId, FrepError> {
        let normal = Self::unit("normal", normal)?;
        Ok(self.push(Node::HalfSpace { normal, offset }))
    }

    /// Axis-aligned box primitive.
    ///
    /// # Errors
    /// [`FrepError::NonPositive`] for a non-positive half extent.
    pub fn box_prim(&mut self, center: Point3, half: Vec3) -> Result<NodeId, FrepError> {
        for (f, v) in [("half.x", half.x), ("half.y", half.y), ("half.z", half.z)] {
            Self::positive(f, v)?;
        }
        Ok(self.push(Node::BoxPrim { center, half }))
    }

    /// Torus primitive (axis +z through `center`).
    ///
    /// # Errors
    /// [`FrepError::NonPositive`] for non-positive radii.
    pub fn torus(&mut self, center: Point3, major: f64, minor: f64) -> Result<NodeId, FrepError> {
        let major = Self::positive("major", major)?;
        let minor = Self::positive("minor", minor)?;
        Ok(self.push(Node::Torus {
            center,
            major,
            minor,
        }))
    }

    /// Infinite cylinder along +z.
    ///
    /// # Errors
    /// [`FrepError::NonPositive`] for `radius ≤ 0`.
    pub fn cylinder(&mut self, center: Point3, radius: f64) -> Result<NodeId, FrepError> {
        let radius = Self::positive("radius", radius)?;
        Ok(self.push(Node::Cylinder { center, radius }))
    }

    /// Rigid translation.
    ///
    /// # Errors
    /// [`FrepError::BadNode`] for an unknown child.
    pub fn translate(&mut self, child: NodeId, offset: Vec3) -> Result<NodeId, FrepError> {
        self.check_child(child)?;
        Ok(self.push(Node::Translate { child, offset }))
    }

    /// Rigid rotation about a unit axis through the origin.
    ///
    /// # Errors
    /// [`FrepError::BadNode`] / [`FrepError::ZeroVector`].
    pub fn rotate(&mut self, child: NodeId, axis: Vec3, angle: f64) -> Result<NodeId, FrepError> {
        self.check_child(child)?;
        let axis = Self::unit("axis", axis)?;
        Ok(self.push(Node::Rotate { child, axis, angle }))
    }

    /// Uniform scale.
    ///
    /// # Errors
    /// [`FrepError::BadNode`] / [`FrepError::NonPositive`].
    pub fn scale(&mut self, child: NodeId, factor: f64) -> Result<NodeId, FrepError> {
        self.check_child(child)?;
        let factor = Self::positive("factor", factor)?;
        Ok(self.push(Node::Scale { child, factor }))
    }

    /// Offset surface (dilate/erode).
    ///
    /// # Errors
    /// [`FrepError::BadNode`] for an unknown child.
    pub fn offset(&mut self, child: NodeId, distance: f64) -> Result<NodeId, FrepError> {
        self.check_child(child)?;
        Ok(self.push(Node::Offset { child, distance }))
    }

    /// Boolean node (any op, hard or blend).
    ///
    /// # Errors
    /// [`FrepError::BadNode`] / [`FrepError::NonPositive`] (blend radius).
    pub fn boolean(
        &mut self,
        op: BoolOp,
        style: BoolStyle,
        a: NodeId,
        b: NodeId,
    ) -> Result<NodeId, FrepError> {
        self.check_child(a)?;
        self.check_child(b)?;
        if let BoolStyle::Blend { radius } = style {
            Self::positive("blend radius", radius)?;
        }
        Ok(self.push(Node::Bool { op, style, a, b }))
    }

    /// Finish, declaring `root` the shape.
    ///
    /// # Errors
    /// [`FrepError::BadNode`] for an unknown root.
    pub fn finish(self, root: NodeId) -> Result<Frep, FrepError> {
        if (root.0 as usize) >= self.nodes.len() {
            return Err(FrepError::BadNode { id: root.0 });
        }
        let mut frep = Frep {
            nodes: self.nodes,
            root,
            support: Aabb::new(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0)),
        };
        frep.support = frep.support_of(frep.root);
        Ok(frep)
    }
}

/// A finished F-rep shape: the DAG plus its root and cached support box.
#[derive(Debug, Clone)]
pub struct Frep {
    nodes: Vec<Node>,
    root: NodeId,
    support: Aabb,
}

/// Addresses one numeric design lever inside the DAG.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParamId {
    /// Owning node.
    pub node: NodeId,
    /// Slot within that node (see [`Frep::params`] for names).
    pub slot: u8,
}

/// Rotate `p` by `angle` around unit `axis` (Rodrigues).
pub(crate) fn rotate_vec(v: Vec3, axis: Vec3, angle: f64) -> Vec3 {
    let (s, c) = angle.sin_cos();
    let kv = Vec3::new(
        axis.y * v.z - axis.z * v.y,
        axis.z * v.x - axis.x * v.z,
        axis.x * v.y - axis.y * v.x,
    );
    let kd = axis.dot(v);
    Vec3::new(
        v.x * c + kv.x * s + axis.x * kd * (1.0 - c),
        v.y * c + kv.y * s + axis.y * kd * (1.0 - c),
        v.z * c + kv.z * s + axis.z * kd * (1.0 - c),
    )
}

/// Quadratic smooth min (C¹): `min(a,b) − r·h²/4`, `h = max(r−|a−b|,0)/r`.
#[must_use]
pub fn smin(a: f64, b: f64, r: f64) -> f64 {
    let h = (r - (a - b).abs()).max(0.0) / r;
    a.min(b) - r * h * h * 0.25
}

/// Blend weights `(wa, wb)` of [`smin`]: nonnegative, `wa + wb = 1`,
/// continuous across `a = b` (the C¹ property frep-003 verifies).
#[must_use]
pub fn smin_weights(a: f64, b: f64, r: f64) -> (f64, f64) {
    let u = a - b;
    if u.abs() >= r {
        if a <= b { (1.0, 0.0) } else { (0.0, 1.0) }
    } else {
        let h = (r - u.abs()) / r;
        if a <= b {
            (1.0 - 0.5 * h, 0.5 * h)
        } else {
            (0.5 * h, 1.0 - 0.5 * h)
        }
    }
}

/// Per-op sign flips: every Boolean routes through ONE smooth/hard min via
/// `sr · min(sa·a, sb·b)` — union (+,+,+), intersect (−,−,−),
/// difference (−,+,−): `a \ b = −min(−a, b) = max(a, −b)`.
pub(crate) fn bool_signs(op: BoolOp) -> (f64, f64, f64) {
    match op {
        BoolOp::Union => (1.0, 1.0, 1.0),
        BoolOp::Intersect => (-1.0, -1.0, -1.0),
        BoolOp::Difference => (-1.0, 1.0, -1.0),
    }
}

impl Frep {
    /// The node list (topological order; read-only).
    #[must_use]
    pub fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    /// The root node.
    #[must_use]
    pub fn root(&self) -> NodeId {
        self.root
    }

    /// Field value at `p` (see module docs for what the value certifies).
    #[must_use]
    pub fn value(&self, p: Point3) -> f64 {
        self.value_at(self.root, p)
    }

    fn value_at(&self, id: NodeId, p: Point3) -> f64 {
        match self.nodes[id.0 as usize] {
            Node::Sphere { center, radius } => p.delta_from(center).norm() - radius,
            Node::HalfSpace { normal, offset } => {
                normal.dot(p.delta_from(Point3::new(0.0, 0.0, 0.0))) - offset
            }
            Node::BoxPrim { center, half } => {
                let d = p.delta_from(center);
                let q = Vec3::new(d.x.abs() - half.x, d.y.abs() - half.y, d.z.abs() - half.z);
                let outside = Vec3::new(q.x.max(0.0), q.y.max(0.0), q.z.max(0.0));
                outside.norm() + q.x.max(q.y).max(q.z).min(0.0)
            }
            Node::Torus {
                center,
                major,
                minor,
            } => {
                let d = p.delta_from(center);
                let ring = d.x.hypot(d.y) - major;
                ring.hypot(d.z) - minor
            }
            Node::Cylinder { center, radius } => {
                let d = p.delta_from(center);
                d.x.hypot(d.y) - radius
            }
            Node::Translate { child, offset } => self.value_at(
                child,
                Point3::new(p.x - offset.x, p.y - offset.y, p.z - offset.z),
            ),
            Node::Rotate { child, axis, angle } => {
                let v = rotate_vec(p.delta_from(Point3::new(0.0, 0.0, 0.0)), axis, -angle);
                self.value_at(child, Point3::new(v.x, v.y, v.z))
            }
            Node::Scale { child, factor } => {
                factor * self.value_at(child, Point3::new(p.x / factor, p.y / factor, p.z / factor))
            }
            Node::Offset { child, distance } => self.value_at(child, p) - distance,
            Node::Bool { op, style, a, b } => {
                let (sa, sb, sr) = bool_signs(op);
                let (fa, fb) = (sa * self.value_at(a, p), sb * self.value_at(b, p));
                match style {
                    BoolStyle::Hard => sr * fa.min(fb),
                    BoolStyle::Blend { radius } => sr * smin(fa, fb, radius),
                }
            }
        }
    }

    /// Field value and chain-rule gradient at `p`. The gradient is `None`
    /// where a primitive has no claim (medial points) — the honest gap
    /// propagates instead of being papered over. On a HARD crease the
    /// selected branch's gradient is returned (a subgradient choice, ties
    /// to the left operand — the discontinuity frep-003 exhibits).
    #[must_use]
    pub fn value_grad(&self, p: Point3) -> (f64, Option<Vec3>) {
        self.vg_at(self.root, p)
    }

    #[allow(clippy::too_many_lines)] // one match arm per node kind: splitting hides the chain rule
    #[allow(clippy::float_cmp)] // compares against exactly-constructed values (max-of, hard-branch 0.0 weights)
    fn vg_at(&self, id: NodeId, p: Point3) -> (f64, Option<Vec3>) {
        match self.nodes[id.0 as usize] {
            Node::Sphere { center, radius } => {
                let d = p.delta_from(center);
                let n = d.norm();
                let g = if n > 1e-12 {
                    Some(d.scale(1.0 / n))
                } else {
                    None
                };
                (n - radius, g)
            }
            Node::HalfSpace { normal, offset } => (
                normal.dot(p.delta_from(Point3::new(0.0, 0.0, 0.0))) - offset,
                Some(normal),
            ),
            Node::BoxPrim { center, half } => {
                let d = p.delta_from(center);
                let q = Vec3::new(d.x.abs() - half.x, d.y.abs() - half.y, d.z.abs() - half.z);
                let outside = Vec3::new(q.x.max(0.0), q.y.max(0.0), q.z.max(0.0));
                let on = outside.norm();
                if on > 1e-12 {
                    let g = Vec3::new(
                        d.x.signum() * outside.x / on,
                        d.y.signum() * outside.y / on,
                        d.z.signum() * outside.z / on,
                    );
                    (on + q.x.max(q.y).max(q.z).min(0.0), Some(g))
                } else {
                    // Inside: the max-coordinate face normal.
                    let m = q.x.max(q.y).max(q.z);
                    let g = if m == q.x {
                        Vec3::new(d.x.signum(), 0.0, 0.0)
                    } else if m == q.y {
                        Vec3::new(0.0, d.y.signum(), 0.0)
                    } else {
                        Vec3::new(0.0, 0.0, d.z.signum())
                    };
                    (m, Some(g))
                }
            }
            Node::Torus {
                center,
                major,
                minor,
            } => {
                let d = p.delta_from(center);
                let s = d.x.hypot(d.y);
                let ring = s - major;
                let m = ring.hypot(d.z);
                if s < 1e-12 || m < 1e-12 {
                    (m - minor, None) // axis or core circle: medial, no claim
                } else {
                    let g = Vec3::new((ring / m) * (d.x / s), (ring / m) * (d.y / s), d.z / m);
                    (m - minor, Some(g))
                }
            }
            Node::Cylinder { center, radius } => {
                let d = p.delta_from(center);
                let s = d.x.hypot(d.y);
                if s < 1e-12 {
                    (s - radius, None)
                } else {
                    (s - radius, Some(Vec3::new(d.x / s, d.y / s, 0.0)))
                }
            }
            Node::Translate { child, offset } => self.vg_at(
                child,
                Point3::new(p.x - offset.x, p.y - offset.y, p.z - offset.z),
            ),
            Node::Rotate { child, axis, angle } => {
                let v = rotate_vec(p.delta_from(Point3::new(0.0, 0.0, 0.0)), axis, -angle);
                let (f, g) = self.vg_at(child, Point3::new(v.x, v.y, v.z));
                // d/dp f(R⁻¹p) = R · ∇f (orthogonal: R⁻ᵀ = R).
                (f, g.map(|g| rotate_vec(g, axis, angle)))
            }
            Node::Scale { child, factor } => {
                let (f, g) =
                    self.vg_at(child, Point3::new(p.x / factor, p.y / factor, p.z / factor));
                (factor * f, g) // ∇[s·f(p/s)] = ∇f(p/s)
            }
            Node::Offset { child, distance } => {
                let (f, g) = self.vg_at(child, p);
                (f - distance, g)
            }
            Node::Bool { op, style, a, b } => {
                let (sa, sb, sr) = bool_signs(op);
                let (fa, ga) = self.vg_at(a, p);
                let (fb, gb) = self.vg_at(b, p);
                let (fa, fb) = (sa * fa, sb * fb);
                let (v, wa, wb) = match style {
                    BoolStyle::Hard => {
                        let v = fa.min(fb);
                        if fa <= fb {
                            (v, 1.0, 0.0)
                        } else {
                            (v, 0.0, 1.0)
                        }
                    }
                    BoolStyle::Blend { radius } => {
                        let (wa, wb) = smin_weights(fa, fb, radius);
                        (smin(fa, fb, radius), wa, wb)
                    }
                };
                let g = match (ga, gb) {
                    // A zero-weight operand's missing gradient is irrelevant.
                    (Some(ga), Some(gb)) => Some(Vec3::new(
                        sr * (wa * sa * ga.x + wb * sb * gb.x),
                        sr * (wa * sa * ga.y + wb * sb * gb.y),
                        sr * (wa * sa * ga.z + wb * sb * gb.z),
                    )),
                    (Some(ga), None) if wb == 0.0 => Some(ga.scale(sr * wa * sa)),
                    (None, Some(gb)) if wa == 0.0 => Some(gb.scale(sr * wb * sb)),
                    _ => None,
                };
                (sr * v, g)
            }
        }
    }

    /// Composed Lipschitz bound of the field (valid everywhere; frep-002
    /// hunts for violations). Primitive fields are unit-Lipschitz; rigid
    /// motions, uniform scale, offsets, and both Boolean styles all
    /// preserve `L ≤ max(L_children)` (the blend's weights are convex).
    #[must_use]
    pub fn lipschitz(&self) -> f64 {
        self.lipschitz_of(self.root)
    }

    fn lipschitz_of(&self, id: NodeId) -> f64 {
        match self.nodes[id.0 as usize] {
            Node::Sphere { .. }
            | Node::HalfSpace { .. }
            | Node::BoxPrim { .. }
            | Node::Torus { .. }
            | Node::Cylinder { .. } => 1.0,
            Node::Translate { child, .. }
            | Node::Rotate { child, .. }
            | Node::Scale { child, .. }
            | Node::Offset { child, .. } => self.lipschitz_of(child),
            Node::Bool { a, b, .. } => self.lipschitz_of(a).max(self.lipschitz_of(b)),
        }
    }

    /// True when any hard Boolean (or box edge) can break C¹ — the
    /// differentiability class the chart advertises.
    fn has_kink(&self) -> bool {
        self.nodes.iter().any(|n| {
            matches!(
                n,
                Node::Bool {
                    style: BoolStyle::Hard,
                    ..
                } | Node::BoxPrim { .. }
            )
        })
    }

    /// True when the field is the EXACT signed distance (a valid ring primitive
    /// under rigid motion and uniform scale, with no Boolean or offset that
    /// would require a reach certificate).
    fn is_exact(&self) -> bool {
        self.nodes.iter().all(|n| match n {
            Node::Bool { .. } => false,
            Node::Offset { .. } => false,
            Node::Torus { major, minor, .. } => major > minor,
            _ => true,
        })
    }

    fn support_of(&self, id: NodeId) -> Aabb {
        match self.nodes[id.0 as usize] {
            Node::Sphere { center, radius } => Aabb::new(
                center.offset(Vec3::new(-radius, -radius, -radius)),
                center.offset(Vec3::new(radius, radius, radius)),
            ),
            Node::HalfSpace { .. } => {
                let u = UNBOUNDED_HALF;
                Aabb::new(Point3::new(-u, -u, -u), Point3::new(u, u, u))
            }
            Node::BoxPrim { center, half } => Aabb::new(
                center.offset(Vec3::new(-half.x, -half.y, -half.z)),
                center.offset(half),
            ),
            Node::Torus {
                center,
                major,
                minor,
            } => {
                let r = major + minor;
                Aabb::new(
                    center.offset(Vec3::new(-r, -r, -minor)),
                    center.offset(Vec3::new(r, r, minor)),
                )
            }
            Node::Cylinder { center, radius } => Aabb::new(
                center.offset(Vec3::new(-radius, -radius, -UNBOUNDED_HALF)),
                center.offset(Vec3::new(radius, radius, UNBOUNDED_HALF)),
            ),
            Node::Translate { child, offset } => {
                let b = self.support_of(child);
                Aabb::new(b.min.offset(offset), b.max.offset(offset))
            }
            Node::Rotate { child, axis, angle } => {
                let b = self.support_of(child);
                let mut out: Option<Aabb> = None;
                for corner in corners(&b) {
                    let v = rotate_vec(corner.delta_from(Point3::new(0.0, 0.0, 0.0)), axis, angle);
                    let p = Point3::new(v.x, v.y, v.z);
                    let cell = Aabb::new(p, p);
                    out = Some(match out {
                        Some(acc) => acc.union(&cell),
                        None => cell,
                    });
                }
                out.expect("a box has corners")
            }
            Node::Scale { child, factor } => {
                let b = self.support_of(child);
                Aabb::new(
                    Point3::new(b.min.x * factor, b.min.y * factor, b.min.z * factor),
                    Point3::new(b.max.x * factor, b.max.y * factor, b.max.z * factor),
                )
            }
            Node::Offset { child, distance } => self.support_of(child).inflate(distance.max(0.0)),
            Node::Bool { op, style, a, b } => {
                let (ba, bb) = (self.support_of(a), self.support_of(b));
                let base = match op {
                    BoolOp::Union => ba.union(&bb),
                    BoolOp::Intersect => intersect_boxes(&ba, &bb),
                    BoolOp::Difference => ba,
                };
                match style {
                    BoolStyle::Hard => base,
                    // The blend adds at most r/4 of material past the hard
                    // surface (deficit of smin vs min is ≤ r/4).
                    BoolStyle::Blend { radius } => base.inflate(radius * 0.25),
                }
            }
        }
    }

    /// Enumerate every design lever: `(id, name, current value)`.
    #[must_use]
    pub fn params(&self) -> Vec<(ParamId, &'static str, f64)> {
        let mut out = Vec::new();
        for (i, node) in self.nodes.iter().enumerate() {
            let node_id = NodeId(i as u32);
            for (slot, (name, value)) in param_slots(node).into_iter().enumerate() {
                out.push((
                    ParamId {
                        node: node_id,
                        slot: slot as u8,
                    },
                    name,
                    value,
                ));
            }
        }
        out
    }

    /// Set one lever (validated like the builder; support box refreshed).
    ///
    /// # Errors
    /// [`FrepError::BadParam`] / [`FrepError::NonPositive`].
    pub fn set_param(&mut self, id: ParamId, value: f64) -> Result<(), FrepError> {
        let node = self
            .nodes
            .get_mut(id.node.0 as usize)
            .ok_or(FrepError::BadNode { id: id.node.0 })?;
        write_slot(node, id.node.0, id.slot, value)?;
        self.support = self.support_of(self.root);
        Ok(())
    }

    /// Jacobian action `∂f(p)/∂θ` for one lever, by symmetric finite
    /// difference on a cloned DAG. HONEST v1: FD with `h = 1e-6·max(1,|θ|)`
    /// (exact parameter adjoints join with fs-xform — contract no-claim).
    ///
    /// # Errors
    /// [`FrepError::BadParam`] / [`FrepError::NonPositive`] (a lever whose
    /// bump would leave the valid domain).
    pub fn d_value_d_param(&self, p: Point3, id: ParamId) -> Result<f64, FrepError> {
        let (_, _, theta) = self
            .params()
            .into_iter()
            .find(|(pid, _, _)| *pid == id)
            .ok_or(FrepError::BadParam {
                node: id.node.0,
                slot: id.slot,
            })?;
        let h = 1e-6 * theta.abs().max(1.0);
        let mut hi = self.clone();
        hi.set_param(id, theta + h)?;
        let mut lo = self.clone();
        lo.set_param(id, theta - h)?;
        Ok((hi.value(p) - lo.value(p)) / (2.0 * h))
    }
}

pub(crate) fn corners(b: &Aabb) -> [Point3; 8] {
    let (lo, hi) = (b.min, b.max);
    [
        Point3::new(lo.x, lo.y, lo.z),
        Point3::new(hi.x, lo.y, lo.z),
        Point3::new(lo.x, hi.y, lo.z),
        Point3::new(hi.x, hi.y, lo.z),
        Point3::new(lo.x, lo.y, hi.z),
        Point3::new(hi.x, lo.y, hi.z),
        Point3::new(lo.x, hi.y, hi.z),
        Point3::new(hi.x, hi.y, hi.z),
    ]
}

fn intersect_boxes(a: &Aabb, b: &Aabb) -> Aabb {
    let min = Point3::new(
        a.min.x.max(b.min.x),
        a.min.y.max(b.min.y),
        a.min.z.max(b.min.z),
    );
    let max = Point3::new(
        a.max.x.min(b.max.x).max(min.x),
        a.max.y.min(b.max.y).max(min.y),
        a.max.z.min(b.max.z).max(min.z),
    );
    Aabb::new(min, max)
}

/// The named parameter slots of a node, in slot order.
fn param_slots(node: &Node) -> Vec<(&'static str, f64)> {
    match *node {
        Node::Sphere { center, radius } | Node::Cylinder { center, radius } => vec![
            ("center.x", center.x),
            ("center.y", center.y),
            ("center.z", center.z),
            ("radius", radius),
        ],
        Node::HalfSpace { offset, .. } => vec![("offset", offset)],
        Node::BoxPrim { center, half } => vec![
            ("center.x", center.x),
            ("center.y", center.y),
            ("center.z", center.z),
            ("half.x", half.x),
            ("half.y", half.y),
            ("half.z", half.z),
        ],
        Node::Torus {
            center,
            major,
            minor,
        } => vec![
            ("center.x", center.x),
            ("center.y", center.y),
            ("center.z", center.z),
            ("major", major),
            ("minor", minor),
        ],
        Node::Translate { offset, .. } => vec![
            ("offset.x", offset.x),
            ("offset.y", offset.y),
            ("offset.z", offset.z),
        ],
        Node::Rotate { angle, .. } => vec![("angle", angle)],
        Node::Scale { factor, .. } => vec![("factor", factor)],
        Node::Offset { distance, .. } => vec![("distance", distance)],
        Node::Bool { style, .. } => match style {
            BoolStyle::Hard => vec![],
            BoolStyle::Blend { radius } => vec![("blend radius", radius)],
        },
    }
}

#[allow(clippy::too_many_lines)] // slot dispatch is one flat table; splitting obscures it
fn write_slot(node: &mut Node, node_ix: u32, slot: u8, value: f64) -> Result<(), FrepError> {
    let positive = |field: &'static str| -> Result<f64, FrepError> {
        if value > 0.0 && value.is_finite() {
            Ok(value)
        } else {
            Err(FrepError::NonPositive { field, value })
        }
    };

    match node {
        Node::Sphere { center, radius } | Node::Cylinder { center, radius } => match slot {
            0 => center.x = value,
            1 => center.y = value,
            2 => center.z = value,
            3 => *radius = positive("radius")?,
            _ => {
                return Err(FrepError::BadParam {
                    node: node_ix,
                    slot,
                });
            }
        },
        Node::HalfSpace { offset, .. } => match slot {
            0 => *offset = value,
            _ => {
                return Err(FrepError::BadParam {
                    node: node_ix,
                    slot,
                });
            }
        },
        Node::BoxPrim { center, half } => match slot {
            0 => center.x = value,
            1 => center.y = value,
            2 => center.z = value,
            3 => half.x = positive("half.x")?,
            4 => half.y = positive("half.y")?,
            5 => half.z = positive("half.z")?,
            _ => {
                return Err(FrepError::BadParam {
                    node: node_ix,
                    slot,
                });
            }
        },
        Node::Torus {
            center,
            major,
            minor,
        } => match slot {
            0 => center.x = value,
            1 => center.y = value,
            2 => center.z = value,
            3 => *major = positive("major")?,
            4 => *minor = positive("minor")?,
            _ => {
                return Err(FrepError::BadParam {
                    node: node_ix,
                    slot,
                });
            }
        },
        Node::Translate { offset, .. } => match slot {
            0 => offset.x = value,
            1 => offset.y = value,
            2 => offset.z = value,
            _ => {
                return Err(FrepError::BadParam {
                    node: node_ix,
                    slot,
                });
            }
        },
        Node::Rotate { angle, .. } => match slot {
            0 => *angle = value,
            _ => {
                return Err(FrepError::BadParam {
                    node: node_ix,
                    slot,
                });
            }
        },
        Node::Scale { factor, .. } => match slot {
            0 => *factor = positive("factor")?,
            _ => {
                return Err(FrepError::BadParam {
                    node: node_ix,
                    slot,
                });
            }
        },
        Node::Offset { distance, .. } => match slot {
            0 => *distance = value,
            _ => {
                return Err(FrepError::BadParam {
                    node: node_ix,
                    slot,
                });
            }
        },
        Node::Bool { style, .. } => match (style, slot) {
            (BoolStyle::Blend { radius }, 0) => *radius = positive("blend radius")?,
            _ => {
                return Err(FrepError::BadParam {
                    node: node_ix,
                    slot,
                });
            }
        },
    }
    Ok(())
}

impl Chart for Frep {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let (f, gradient) = self.value_grad(x);
        ChartSample {
            signed_distance: f,
            gradient,
            lipschitz: Some(self.lipschitz()),
            // Composite fields are a conservative bound, not the exact
            // distance (module docs): sign exact, |f| ≤ true distance.
            error: if self.is_exact() {
                NumericalCertificate::exact(f)
            } else {
                NumericalCertificate::estimate(f, f)
            },
        }
    }

    fn support(&self) -> Aabb {
        self.support
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        if self.is_exact() {
            TraceStepClaim::ExactDistance
        } else {
            TraceStepClaim::LipschitzImplicit
        }
    }

    fn topology_hint(&self) -> BettiBounds {
        BettiBounds::unknown()
    }

    fn name(&self) -> &'static str {
        "frep/csg"
    }

    fn differentiability(&self) -> Differentiability {
        if self.has_kink() {
            Differentiability::C0
        } else {
            Differentiability::C1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }

    #[test]
    fn builder_teaches_on_bad_input() {
        let mut b = FrepBuilder::new();
        let err = b.sphere(Point3::new(0.0, 0.0, 0.0), -1.0).unwrap_err();
        assert!(err.to_string().contains("strictly positive"), "{err}");
        let err = b.half_space(Vec3::new(0.0, 0.0, 0.0), 1.0).unwrap_err();
        assert!(err.to_string().contains("normalize"), "{err}");
        let err = b.finish(NodeId(7)).unwrap_err();
        assert!(err.to_string().contains("does not exist"), "{err}");
    }

    #[test]
    fn blend_self_union_is_quarter_radius_dilation() {
        // smin(a, a, r) = a − r/4 exactly (h = 1): a bitwise metamorphic law.
        let mut b = FrepBuilder::new();
        let s = b.sphere(Point3::new(0.1, -0.2, 0.3), 0.9).unwrap();
        let u = b
            .boolean(BoolOp::Union, BoolStyle::Blend { radius: 0.25 }, s, s)
            .unwrap();
        let f = b.finish(u).unwrap();
        let p = Point3::new(0.7, 0.4, -0.5);
        let plain = p.delta_from(Point3::new(0.1, -0.2, 0.3)).norm() - 0.9;
        assert_eq!(
            f.value(p).to_bits(),
            (plain - 0.25 * 0.25).to_bits(),
            "a ∪_blend a = offset(a, r/4) bitwise"
        );
    }
}
