//! fs-geom — the Region/Chart abstraction (plan §7.1): the geometry
//! kernel's founding move. Layer: L2.
//!
//! An abstract `Region` — semantically a measurable subset of ℝ³ with
//! piecewise-smooth boundary — is NEVER stored directly. It is PRESENTED
//! through [`Chart`]s: concrete representations answering signed-distance
//! queries with a value, a gradient where one exists, a certified local
//! Lipschitz bound, and a DECLARED error model relative to the abstract
//! region (an [`fs_evidence::NumericalCertificate`]).
//!
//! Three consequences no conventional kernel delivers (plan §7.1):
//! - "the same shape held three ways" is a normal, coherent state:
//!   [`Region`] holds multiple charts with provenance, and AGREEMENT
//!   BETWEEN CHARTS IS A CHECKABLE PROPOSITION
//!   ([`Region::check_agreement`]) with localized diagnostics, not an
//!   assumption;
//! - every conversion's error is explicit and composable: [`Convert`]
//!   returns a [`Certified`] receipt (fs-evidence) feeding the Error
//!   Ledger (Decalogue P4);
//! - no chart type is privileged, ever — the Rep Router (a later bead)
//!   picks per OPERATION from declared capabilities.
//!
//! Cancellation: every query takes `&fs_exec::Cx` — geometry is
//! interruptible at bounded strides like any other kernel (Decalogue P7).
//!
//! Object safety note: plan Appendix B sketches `Chart { type Param; ... }`.
//! `Region` must hold heterogeneous charts (`Arc<dyn Chart>`), so the
//! design-lever handle lives on the [`DesignChart`] subtrait instead —
//! same contract, object-safe core (fs-xform builds on `DesignChart`).

use fs_evidence::NumericalCertificate;
use fs_exec::Cx;

mod convert;
pub mod fixtures;
pub mod router;
mod region;

pub use convert::{Convert, ConvertDiag, ErrBudget, SampledSdf};
pub use router::{
    Binding, ChainOutcome, ConverterSpec, CostOracle, EdgeOutcome, EdgeRunner, ErrorModel,
    ExecuteError, MemoryCostOracle, RouteCandidate, RouteExplanation, RoutePlan, RouteRefusal,
    RouteRequest, Router, RouterError,
};
pub use region::{AgreementConfig, AgreementReport, Disagreement, Region, RegionChart};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// A point in ℝ³. Minimal geometry-local type (fs-la owns real linear
/// algebra; these exist so charts need no L1 dependency).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3 {
    /// x coordinate.
    pub x: f64,
    /// y coordinate.
    pub y: f64,
    /// z coordinate.
    pub z: f64,
}

impl Point3 {
    /// Construct a point.
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Point3 { x, y, z }
    }

    /// Difference vector `self - other`.
    #[must_use]
    pub fn delta_from(self, other: Point3) -> Vec3 {
        Vec3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    /// Translate by a vector.
    #[must_use]
    pub fn offset(self, v: Vec3) -> Point3 {
        Point3 {
            x: self.x + v.x,
            y: self.y + v.y,
            z: self.z + v.z,
        }
    }
}

/// A vector in ℝ³ (see [`Point3`]'s scope note).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    /// x component.
    pub x: f64,
    /// y component.
    pub y: f64,
    /// z component.
    pub z: f64,
}

impl Vec3 {
    /// Construct a vector.
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    /// Euclidean norm.
    #[must_use]
    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }

    /// Dot product.
    #[must_use]
    pub fn dot(self, o: Vec3) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Scale by a scalar.
    #[must_use]
    pub fn scale(self, s: f64) -> Vec3 {
        Vec3 {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }
}

/// An axis-aligned bounding box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb {
    /// Componentwise minimum corner.
    pub min: Point3,
    /// Componentwise maximum corner.
    pub max: Point3,
}

impl Aabb {
    /// Box from corners (normalized componentwise, total function).
    #[must_use]
    pub fn new(a: Point3, b: Point3) -> Self {
        Aabb {
            min: Point3::new(a.x.min(b.x), a.y.min(b.y), a.z.min(b.z)),
            max: Point3::new(a.x.max(b.x), a.y.max(b.y), a.z.max(b.z)),
        }
    }

    /// True when `p` lies inside or on the boundary.
    #[must_use]
    pub fn contains(&self, p: Point3) -> bool {
        p.x >= self.min.x
            && p.x <= self.max.x
            && p.y >= self.min.y
            && p.y <= self.max.y
            && p.z >= self.min.z
            && p.z <= self.max.z
    }

    /// Smallest box containing both.
    #[must_use]
    pub fn union(&self, other: &Aabb) -> Aabb {
        Aabb::new(
            Point3::new(
                self.min.x.min(other.min.x),
                self.min.y.min(other.min.y),
                self.min.z.min(other.min.z),
            ),
            Point3::new(
                self.max.x.max(other.max.x),
                self.max.y.max(other.max.y),
                self.max.z.max(other.max.z),
            ),
        )
    }

    /// Grow outward by `pad` on every side.
    #[must_use]
    pub fn inflate(&self, pad: f64) -> Aabb {
        Aabb::new(
            Point3::new(self.min.x - pad, self.min.y - pad, self.min.z - pad),
            Point3::new(self.max.x + pad, self.max.y + pad, self.max.z + pad),
        )
    }
}

/// Betti-number bounds `(lower, upper)` per dimension — the topology hint
/// charts may advertise without proving (certificates are a later bead;
/// `unknown()` is the honest default).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BettiBounds {
    /// Connected components.
    pub b0: (u32, u32),
    /// Tunnels/handles.
    pub b1: (u32, u32),
    /// Enclosed voids.
    pub b2: (u32, u32),
}

impl BettiBounds {
    /// No topology claim at all.
    #[must_use]
    pub const fn unknown() -> Self {
        BettiBounds {
            b0: (0, u32::MAX),
            b1: (0, u32::MAX),
            b2: (0, u32::MAX),
        }
    }

    /// Exact known Betti numbers.
    #[must_use]
    pub const fn exact(b0: u32, b1: u32, b2: u32) -> Self {
        BettiBounds {
            b0: (b0, b0),
            b1: (b1, b1),
            b2: (b2, b2),
        }
    }
}

/// Differentiability class a chart advertises for its signed distance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Differentiability {
    /// Continuous only (gradients may be `None` anywhere).
    C0,
    /// Continuously differentiable away from the medial axis.
    C1,
    /// Smooth away from the medial axis.
    Smooth,
}

/// One signed-distance query's answer (plan Appendix B: value + gradient +
/// certified Lipschitz data + the declared error model).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChartSample {
    /// Signed distance to the region boundary (negative inside — the SDF
    /// convention every chart maps onto, whatever its native form).
    pub signed_distance: f64,
    /// Gradient of the signed distance where it exists (`None` on medial
    /// axis/edges or for C0 charts).
    pub gradient: Option<Vec3>,
    /// Certified LOCAL Lipschitz bound for the signed distance near the
    /// query (sphere-tracing fuel; `None` = no claim).
    pub lipschitz: Option<f64>,
    /// Declared error of `signed_distance` relative to the ABSTRACT region
    /// (fs-evidence certificate: exact charts say Exact, sampled charts say
    /// Enclosure, heuristics say Estimate).
    pub error: NumericalCertificate,
}

/// The chart contract: a concrete presentation of an abstract region.
/// Object-safe (see module docs for the `Param` note).
pub trait Chart: Send + Sync {
    /// Answer a signed-distance query. Implementations poll
    /// `cx.checkpoint()` at bounded strides inside expensive evaluations.
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample;

    /// A box guaranteed to contain the region (queries outside are
    /// positive-distance by definition).
    fn support(&self) -> Aabb;

    /// Topology bounds this chart is willing to state.
    fn topology_hint(&self) -> BettiBounds {
        BettiBounds::unknown()
    }

    /// Stable chart-kind name (provenance, reports, router tables).
    fn name(&self) -> &'static str;

    /// Advertised differentiability class of the signed distance.
    fn differentiability(&self) -> Differentiability {
        Differentiability::C0
    }

    /// Convenience: strict inside test via the SDF convention.
    fn inside(&self, x: Point3, cx: &Cx<'_>) -> bool {
        self.eval(x, cx).signed_distance < 0.0
    }
}

/// A chart with design levers: the differentiable map θ → Region handle
/// (plan §7.6; fs-xform builds the parameterization zoo on this).
pub trait DesignChart: Chart {
    /// The design-lever handle.
    type Param;

    /// Report the current lever value.
    fn param(&self) -> &Self::Param;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }

    #[test]
    fn aabb_laws() {
        let a = Aabb::new(Point3::new(1.0, 2.0, 3.0), Point3::new(-1.0, 0.0, 5.0));
        assert_eq!(a.min, Point3::new(-1.0, 0.0, 3.0), "corners normalize");
        assert!(a.contains(Point3::new(0.0, 1.0, 4.0)));
        assert!(!a.contains(Point3::new(2.0, 1.0, 4.0)));
        let b = Aabb::new(Point3::new(5.0, 5.0, 5.0), Point3::new(6.0, 6.0, 6.0));
        let u = a.union(&b);
        assert!(u.contains(Point3::new(5.5, 5.5, 5.5)) && u.contains(Point3::new(0.0, 1.0, 4.0)));
        assert!(a.inflate(1.0).contains(Point3::new(1.5, 0.5, 4.0)));
    }

    #[test]
    fn vec_ops_are_the_usual_ones() {
        let v = Point3::new(1.0, 2.0, 2.0).delta_from(Point3::new(0.0, 0.0, 0.0));
        assert!((v.norm() - 3.0).abs() < 1e-12);
        assert!((v.dot(Vec3::new(1.0, 0.0, 0.0)) - 1.0).abs() < 1e-12);
        assert_eq!(v.scale(2.0), Vec3::new(2.0, 4.0, 4.0));
    }

    #[test]
    fn betti_bounds_default_to_no_claim() {
        let u = BettiBounds::unknown();
        assert_eq!(u.b0, (0, u32::MAX));
        assert_eq!(BettiBounds::exact(1, 0, 1).b2, (1, 1));
    }
}
