//! Symmetry ENFORCED BY CONSTRUCTION: the quotient chart folds every
//! evaluation point into the group's fundamental domain before asking
//! the inner chart, so the presented shape is the symmetric orbit of
//! the inner design for ARBITRARY lever values — a symmetry constraint
//! that is structurally impossible to violate.
//!
//! Honesty note: `f(fold(p))` is the EXACT field of the symmetrized
//! shape when the inner design stays inside the fundamental domain;
//! designs that spill across the seam still produce an exactly
//! symmetric shape, but the field magnitude near the seam becomes a
//! conservative bound (sign remains exact) — the same contract every
//! CSG field in this workspace carries.

use fs_exec::Cx;
use fs_geom::{Aabb, Chart, ChartSample, Point3, Vec3};

/// The supported symmetry groups (v1: the plan's named trio).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SymmetryGroup {
    /// Reflection across the plane `x = 0`.
    ReflectX,
    /// Cyclic rotation `C_n` about the +z axis.
    Cyclic {
        /// Order (≥ 2).
        n: u32,
    },
    /// Translational repetition along x with the given period.
    Periodic {
        /// Period (> 0).
        period: f64,
    },
}

impl SymmetryGroup {
    /// Fold `p` into the fundamental domain.
    #[must_use]
    pub fn fold(&self, p: Point3) -> Point3 {
        match *self {
            SymmetryGroup::ReflectX => Point3::new(p.x.abs(), p.y, p.z),
            SymmetryGroup::Cyclic { n } => {
                let sector = core::f64::consts::TAU / f64::from(n);
                let r = p.x.hypot(p.y);
                if r < 1e-300 {
                    return p;
                }
                let theta = p.y.atan2(p.x).rem_euclid(sector);
                Point3::new(r * theta.cos(), r * theta.sin(), p.z)
            }
            SymmetryGroup::Periodic { period } => Point3::new(p.x.rem_euclid(period), p.y, p.z),
        }
    }

    /// Push a gradient at the folded point back to the query point
    /// (the chain rule through the fold, valid in the domain interior).
    #[must_use]
    pub fn unfold_gradient(&self, p: Point3, g: Vec3) -> Vec3 {
        match *self {
            SymmetryGroup::ReflectX => {
                if p.x < 0.0 {
                    Vec3::new(-g.x, g.y, g.z)
                } else {
                    g
                }
            }
            SymmetryGroup::Cyclic { n } => {
                // The fold is a rotation by −k·sector; rotate g back.
                let sector = core::f64::consts::TAU / f64::from(n);
                let theta = p.y.atan2(p.x);
                let k = (theta.rem_euclid(core::f64::consts::TAU) / sector).floor();
                let ang = k * sector;
                let (s, c) = ang.sin_cos();
                Vec3::new(c * g.x - s * g.y, s * g.x + c * g.y, g.z)
            }
            SymmetryGroup::Periodic { .. } => g,
        }
    }

    /// A representative set of nontrivial group elements applied to a
    /// point (the audit's probes).
    #[must_use]
    pub fn orbit(&self, p: Point3) -> Vec<Point3> {
        match *self {
            SymmetryGroup::ReflectX => vec![Point3::new(-p.x, p.y, p.z)],
            SymmetryGroup::Cyclic { n } => (1..n)
                .map(|k| {
                    let ang = core::f64::consts::TAU * f64::from(k) / f64::from(n);
                    let (s, c) = ang.sin_cos();
                    Point3::new(c * p.x - s * p.y, s * p.x + c * p.y, p.z)
                })
                .collect(),
            SymmetryGroup::Periodic { period } => vec![
                Point3::new(p.x + period, p.y, p.z),
                Point3::new(p.x - 2.0 * period, p.y, p.z),
            ],
        }
    }
}

/// The quotient chart: inner design levers act on the fundamental
/// domain; the presented shape is the symmetric orbit.
pub struct QuotientChart<'a> {
    /// The inner (fundamental-domain) chart.
    pub inner: &'a dyn Chart,
    /// The group.
    pub group: SymmetryGroup,
}

impl Chart for QuotientChart<'_> {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        let folded = self.group.fold(x);
        let mut s = self.inner.eval(folded, cx);
        s.gradient = s.gradient.map(|g| self.group.unfold_gradient(x, g));
        s
    }

    fn support(&self) -> Aabb {
        // The orbit of the inner support: conservative bounding box.
        let b = self.inner.support();
        match self.group {
            SymmetryGroup::ReflectX => Aabb::new(
                Point3::new(-b.max.x.abs().max(b.min.x.abs()), b.min.y, b.min.z),
                Point3::new(b.max.x.abs().max(b.min.x.abs()), b.max.y, b.max.z),
            ),
            SymmetryGroup::Cyclic { .. } => {
                let r = [b.min.x.abs(), b.max.x.abs(), b.min.y.abs(), b.max.y.abs()]
                    .into_iter()
                    .fold(0.0f64, f64::max)
                    * core::f64::consts::SQRT_2;
                Aabb::new(Point3::new(-r, -r, b.min.z), Point3::new(r, r, b.max.z))
            }
            SymmetryGroup::Periodic { .. } => {
                // Unbounded along x in principle; report a wide box.
                let u = 1.0e12;
                Aabb::new(
                    Point3::new(-u, b.min.y, b.min.z),
                    Point3::new(u, b.max.y, b.max.z),
                )
            }
        }
    }

    fn name(&self) -> &'static str {
        "geocon/quotient"
    }

    fn differentiability(&self) -> fs_geom::Differentiability {
        // The seam introduces C0 creases even for smooth inners.
        fs_geom::Differentiability::C0
    }
}
