//! Symmetry ENFORCED BY CONSTRUCTION: the quotient chart folds every
//! evaluation point into the group's fundamental domain before asking
//! the inner chart, so the presented shape is the symmetric orbit of
//! the inner design for ARBITRARY lever values — a symmetry constraint
//! that is structurally impossible to violate.
//!
//! Honesty note: `f(fold(p))` is a symmetric PRESENTATION, but a generic
//! fold can change abstract-distance magnitude and introduce a seam. The raw
//! quotient therefore carries no certified Lipschitz theorem and demotes
//! finite inner evidence to `Estimate`; `NoClaim` remains absorbing.

use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::Cx;
use fs_geom::{Aabb, Chart, ChartSample, Point3, Vec3};

/// Invalid public symmetry-group parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SymmetryError {
    /// A cyclic group must contain at least two rotations.
    InvalidOrder {
        /// Rejected group order.
        n: u32,
    },
    /// A translation period must be finite and strictly positive.
    InvalidPeriod {
        /// Raw bits of the rejected period.
        period_bits: u64,
    },
}

impl core::fmt::Display for SymmetryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidOrder { n } => {
                write!(f, "cyclic symmetry order must be at least 2, got {n}")
            }
            Self::InvalidPeriod { period_bits } => write!(
                f,
                "periodic symmetry period must be finite and positive, got f64 bits {period_bits:#018x}"
            ),
        }
    }
}

impl core::error::Error for SymmetryError {}

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
    /// Validated cyclic group constructor.
    pub fn cyclic(n: u32) -> Result<Self, SymmetryError> {
        let group = Self::Cyclic { n };
        group.validate()?;
        Ok(group)
    }

    /// Validated periodic group constructor.
    pub fn periodic(period: f64) -> Result<Self, SymmetryError> {
        let group = Self::Periodic { period };
        group.validate()?;
        Ok(group)
    }

    /// Validate even a directly constructed public enum value.
    pub fn validate(self) -> Result<(), SymmetryError> {
        match self {
            Self::ReflectX => Ok(()),
            Self::Cyclic { n } if n >= 2 => Ok(()),
            Self::Cyclic { n } => Err(SymmetryError::InvalidOrder { n }),
            Self::Periodic { period } if period.is_finite() && period > 0.0 => Ok(()),
            Self::Periodic { period } => Err(SymmetryError::InvalidPeriod {
                period_bits: period.to_bits(),
            }),
        }
    }

    /// Fold `p` into the fundamental domain.
    #[must_use]
    pub fn fold(&self, p: Point3) -> Point3 {
        if self.validate().is_err() {
            return Point3::new(f64::NAN, f64::NAN, f64::NAN);
        }
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
        if self.validate().is_err() {
            return Vec3::new(f64::NAN, f64::NAN, f64::NAN);
        }
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
        if self.validate().is_err() {
            return Vec::new();
        }
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
        if self.group.validate().is_err() {
            return ChartSample {
                signed_distance: f64::NAN,
                gradient: None,
                lipschitz: None,
                error: NumericalCertificate::no_claim(),
            };
        }
        let folded = self.group.fold(x);
        let mut s = self.inner.eval(folded, cx);
        s.gradient = s.gradient.map(|g| self.group.unfold_gradient(x, g));
        s.lipschitz = None;
        s.error = if quotient_input_certificate_is_valid(&s) {
            NumericalCertificate::estimate(
                s.error.lo.min(s.signed_distance),
                s.error.hi.max(s.signed_distance),
            )
        } else {
            NumericalCertificate::no_claim()
        };
        s
    }

    fn support(&self) -> Aabb {
        // The orbit of the inner support: conservative bounding box.
        let b = self.inner.support();
        if !b.is_well_formed() {
            return b;
        }
        if self.group.validate().is_err() {
            return Aabb {
                min: Point3::new(f64::NAN, b.min.y, b.min.z),
                max: b.max,
            };
        }
        match self.group {
            SymmetryGroup::ReflectX => Aabb::new(
                Point3::new(-b.max.x.abs().max(b.min.x.abs()), b.min.y, b.min.z),
                Point3::new(b.max.x.abs().max(b.min.x.abs()), b.max.y, b.max.z),
            ),
            SymmetryGroup::Cyclic { .. } => {
                let max_abs = [b.min.x.abs(), b.max.x.abs(), b.min.y.abs(), b.max.y.abs()]
                    .into_iter()
                    .fold(0.0f64, f64::max);
                // The farthest possible rotated corner is no farther than
                // sqrt(2) * max_abs. Both the irrational factor and product
                // are rounded upward so this AABB preserves Chart::support's
                // containment promise at floating-point boundaries.
                let radius = max_abs * core::f64::consts::SQRT_2.next_up();
                let r = if radius == 0.0 || radius.is_infinite() {
                    radius
                } else {
                    radius.next_up()
                };
                Aabb::new(Point3::new(-r, -r, b.min.z), Point3::new(r, r, b.max.z))
            }
            SymmetryGroup::Periodic { .. } => {
                // Translational repetition is genuinely unbounded along x.
                Aabb::new(
                    Point3::new(f64::NEG_INFINITY, b.min.y, b.min.z),
                    Point3::new(f64::INFINITY, b.max.y, b.max.z),
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

/// Quotienting weakens abstract-distance authority but must not repair a
/// malformed inner certificate merely because the nominal field is finite.
fn quotient_input_certificate_is_valid(sample: &ChartSample) -> bool {
    if !sample.signed_distance.is_finite()
        || !sample.error.lo.is_finite()
        || !sample.error.hi.is_finite()
        || sample.error.lo > sample.error.hi
    {
        return false;
    }
    match sample.error.kind {
        NumericalKind::Exact => {
            sample.error.lo.to_bits() == sample.signed_distance.to_bits()
                && sample.error.hi.to_bits() == sample.signed_distance.to_bits()
        }
        NumericalKind::Enclosure | NumericalKind::Estimate => {
            sample.error.lo <= sample.signed_distance && sample.signed_distance <= sample.error.hi
        }
        NumericalKind::NoClaim => false,
    }
}
