//! Certified motor tubes: piecewise Taylor-model enclosures of a
//! motor path, with deterministic double-cover handling, validated
//! chart transitions, and interval actions on points and boxes.

use crate::MotionError;
use crate::algebra::{BLADES, SCALAR, TmMv, homogeneous_point, point_to_mv};
use fs_exec::Cx;
use fs_ga::{Motor, Pga};
use fs_geom::{Aabb, Point3};
use fs_ivl::Interval;

/// How strongly an output is enclosed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnclosureClass {
    /// Rigorous Taylor-model enclosure (remainder included).
    Certified,
    /// Only checked by finite sampling; sampling falsifies, never
    /// proves. No production consumer may upgrade this.
    FalsifiedOnly,
}

/// Midpoint scalar magnitude below which the double-cover sign choice
/// refuses instead of guessing.
const SIGN_ANCHOR_TOL: f64 = 1e-9;

/// One tube segment: sixteen component models over one domain, plus
/// the rigorously computed versor-defect bound.
#[derive(Debug, Clone)]
pub struct MotorTubeSegment {
    mv: TmMv,
    defect: f64,
}

impl MotorTubeSegment {
    /// Seal a segment: canonicalize the double-cover sign and record
    /// the versor defect. The sign anchor is the component midpoint at
    /// the domain's low endpoint, scanned in fixed blade order
    /// starting from the scalar; all-tiny anchors refuse as ambiguous.
    pub fn seal(mv: TmMv) -> Result<MotorTubeSegment, MotionError> {
        let t0 = Interval::point(mv.domain().lo());
        let at_start = mv.eval_all(t0)?;
        let mut flip: Option<bool> = None;
        let scalar_mid = at_start[SCALAR].midpoint();
        if scalar_mid.abs() > SIGN_ANCHOR_TOL {
            flip = Some(scalar_mid < 0.0);
        } else {
            for enc in at_start.iter().skip(1) {
                let mid = enc.midpoint();
                if mid.abs() > SIGN_ANCHOR_TOL {
                    flip = Some(mid < 0.0);
                    break;
                }
            }
        }
        let Some(flip) = flip else {
            return Err(MotionError::DoubleCoverAmbiguous {
                at: mv.domain().lo(),
            });
        };
        let mv = if flip { mv.negate()? } else { mv };
        let defect = mv.versor_defect()?;
        Ok(MotorTubeSegment { mv, defect })
    }

    /// The segment's time domain.
    #[must_use]
    pub fn domain(&self) -> Interval {
        self.mv.domain()
    }

    /// Upper bound of `‖M M̃ − 1‖∞` over the domain.
    #[must_use]
    pub fn defect(&self) -> f64 {
        self.defect
    }

    /// Component enclosures over a subinterval.
    pub fn components_over(&self, t: Interval) -> Result<[Interval; BLADES], MotionError> {
        self.mv.eval_all(t)
    }

    /// Crate-internal access to the component models.
    pub(crate) fn components(&self) -> &TmMv {
        &self.mv
    }

}

/// Honesty data returned alongside a pointwise motor evaluation.
#[derive(Debug, Clone, Copy)]
pub struct PathSample {
    /// Largest component-enclosure width at the queried time.
    pub max_enclosure_width: f64,
    /// The owning segment's versor-defect bound.
    pub defect: f64,
}

/// Interval enclosure of a moving point, with class and defect.
#[derive(Debug, Clone, Copy)]
pub struct PointActionEnclosure {
    /// Componentwise enclosure of `M(t)·x` over the queried span.
    pub coords: [Interval; 3],
    /// Enclosure strength.
    pub class: EnclosureClass,
    /// Versor-defect bound of the segments involved.
    pub defect: f64,
}

/// Interval enclosure of a moving box.
#[derive(Debug, Clone, Copy)]
pub struct BoxActionEnclosure {
    /// An AABB containing the image of the box for every queried time.
    pub bounds: Aabb,
    /// Enclosure strength.
    pub class: EnclosureClass,
    /// Versor-defect bound of the segments involved.
    pub defect: f64,
}

/// A piecewise certified enclosure of a motor path.
#[derive(Debug, Clone)]
pub struct CertifiedMotorTube {
    segments: Vec<MotorTubeSegment>,
}

impl CertifiedMotorTube {
    /// Assemble a tube from contiguous sealed segments, validating
    /// every interior chart transition BEFORE any consumer takes a
    /// logarithm: adjacent component enclosures at the shared boundary
    /// must intersect, and the representative vectors must have a
    /// positive dot product (a non-positive dot product is a
    /// double-cover flip across the boundary).
    pub fn from_segments(
        segments: Vec<MotorTubeSegment>,
    ) -> Result<CertifiedMotorTube, MotionError> {
        if segments.is_empty() {
            return Err(MotionError::EmptyTimeDomain);
        }
        for pair in segments.windows(2) {
            let t = pair[0].domain().hi();
            if pair[1].domain().lo() != t {
                return Err(MotionError::ChartTransition {
                    at: t,
                    dot: f64::NAN,
                });
            }
            let left = pair[0].components_over(Interval::point(t))?;
            let right = pair[1].components_over(Interval::point(t))?;
            let mut dot = 0.0f64;
            for (l, r) in left.iter().zip(right.iter()) {
                if l.lo() > r.hi() || r.lo() > l.hi() {
                    return Err(MotionError::ChartTransition { at: t, dot: f64::NAN });
                }
                dot += l.midpoint() * r.midpoint();
            }
            if !(dot > 0.0) {
                return Err(MotionError::ChartTransition { at: t, dot });
            }
        }
        Ok(CertifiedMotorTube { segments })
    }

    /// The tube's full time domain.
    #[must_use]
    pub fn domain(&self) -> Interval {
        let lo = self.segments[0].domain().lo();
        let hi = self.segments[self.segments.len() - 1].domain().hi();
        Interval::new(lo, hi)
    }

    /// The segments (read-only).
    #[must_use]
    pub fn segments(&self) -> &[MotorTubeSegment] {
        &self.segments
    }

    /// Worst versor-defect bound across all segments.
    #[must_use]
    pub fn defect(&self) -> f64 {
        self.segments
            .iter()
            .map(MotorTubeSegment::defect)
            .fold(0.0, f64::max)
    }

    fn segment_for(&self, t: f64) -> Result<&MotorTubeSegment, MotionError> {
        let domain = self.domain();
        if !(t >= domain.lo() && t <= domain.hi()) {
            return Err(MotionError::OutOfDomain {
                lo: t,
                hi: t,
                domain_lo: domain.lo(),
                domain_hi: domain.hi(),
            });
        }
        for seg in &self.segments {
            if t <= seg.domain().hi() {
                return Ok(seg);
            }
        }
        Ok(&self.segments[self.segments.len() - 1])
    }

    /// Segments overlapping a query span, with the span clamped to
    /// each segment's domain.
    fn segments_over(
        &self,
        span: Interval,
    ) -> Result<Vec<(&MotorTubeSegment, Interval)>, MotionError> {
        let domain = self.domain();
        if !domain.encloses(span) {
            return Err(MotionError::OutOfDomain {
                lo: span.lo(),
                hi: span.hi(),
                domain_lo: domain.lo(),
                domain_hi: domain.hi(),
            });
        }
        let mut out = Vec::new();
        for seg in &self.segments {
            let d = seg.domain();
            let lo = span.lo().max(d.lo());
            let hi = span.hi().min(d.hi());
            if lo <= hi {
                out.push((seg, Interval::new(lo, hi)));
            }
        }
        Ok(out)
    }

    /// Certified enclosure of `M(t)·x` for all `t` in `span`.
    pub fn point_action_over(
        &self,
        x: Point3,
        span: Interval,
        cx: &Cx<'_>,
    ) -> Result<PointActionEnclosure, MotionError> {
        if !(x.x.is_finite() && x.y.is_finite() && x.z.is_finite()) {
            return Err(MotionError::NonFiniteInput { what: "point" });
        }
        let parts = self.segments_over(span)?;
        let mut coords: Option<[Interval; 3]> = None;
        let mut defect = 0.0f64;
        for (seg, sub) in parts {
            cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
            let p = TmMv::constant(&point_to_mv(x.x, x.y, x.z), seg.domain(), seg.mv.order())?;
            let sandwich = seg.mv.gp(&p)?.gp(&seg.mv.reverse()?)?;
            let enc = sandwich.eval_all(sub)?;
            let pt = homogeneous_point(&enc)?;
            defect = defect.max(seg.defect());
            coords = Some(match coords {
                None => pt,
                Some(prev) => [prev[0].hull(pt[0]), prev[1].hull(pt[1]), prev[2].hull(pt[2])],
            });
        }
        let coords = coords.ok_or(MotionError::EmptyTimeDomain)?;
        Ok(PointActionEnclosure {
            coords,
            class: EnclosureClass::Certified,
            defect,
        })
    }

    /// Certified enclosure of the image of an AABB for all `t` in
    /// `span`. The action is affine in `x` at each fixed time, so the
    /// hull of the eight corner enclosures contains the image of every
    /// point of the box.
    pub fn box_action_over(
        &self,
        b: &Aabb,
        span: Interval,
        cx: &Cx<'_>,
    ) -> Result<BoxActionEnclosure, MotionError> {
        let corners = [
            Point3::new(b.min.x, b.min.y, b.min.z),
            Point3::new(b.min.x, b.min.y, b.max.z),
            Point3::new(b.min.x, b.max.y, b.min.z),
            Point3::new(b.min.x, b.max.y, b.max.z),
            Point3::new(b.max.x, b.min.y, b.min.z),
            Point3::new(b.max.x, b.min.y, b.max.z),
            Point3::new(b.max.x, b.max.y, b.min.z),
            Point3::new(b.max.x, b.max.y, b.max.z),
        ];
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        let mut defect = 0.0f64;
        for corner in corners {
            cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
            let enc = self.point_action_over(corner, span, cx)?;
            for a in 0..3 {
                lo[a] = lo[a].min(enc.coords[a].lo());
                hi[a] = hi[a].max(enc.coords[a].hi());
            }
            defect = defect.max(enc.defect);
        }
        Ok(BoxActionEnclosure {
            bounds: Aabb::new(
                Point3::new(lo[0], lo[1], lo[2]),
                Point3::new(hi[0], hi[1], hi[2]),
            ),
            class: EnclosureClass::Certified,
            defect,
        })
    }

    /// The point-evaluation view.
    #[must_use]
    pub fn path(&self) -> MotorPath<'_> {
        MotorPath { tube: self }
    }
}

/// Point-evaluation view of a tube: representative motors at single
/// times, with honesty data.
#[derive(Debug, Clone, Copy)]
pub struct MotorPath<'a> {
    tube: &'a CertifiedMotorTube,
}

impl MotorPath<'_> {
    /// The representative motor at `t` (componentwise enclosure
    /// midpoints) plus the enclosure width and segment defect. The
    /// motor is NOT renormalized; callers deciding to renormalize own
    /// that policy and its drift ledger entry.
    pub fn motor_at(&self, t: f64) -> Result<(Motor, PathSample), MotionError> {
        if !t.is_finite() {
            return Err(MotionError::NonFiniteInput { what: "time" });
        }
        let seg = self.tube.segment_for(t)?;
        let enc = seg.components_over(Interval::point(t))?;
        let mut arr = [0.0f64; BLADES];
        let mut width = 0.0f64;
        for (slot, e) in arr.iter_mut().zip(enc.iter()) {
            *slot = e.midpoint();
            width = width.max(e.width());
        }
        let mut pga = Pga::zero();
        pga.0 = arr;
        Ok((
            Motor(pga),
            PathSample {
                max_enclosure_width: width,
                defect: seg.defect(),
            },
        ))
    }

    /// The tube's domain.
    #[must_use]
    pub fn domain(&self) -> Interval {
        self.tube.domain()
    }
}

/// Builder contract: higher layers (frame trees, MBD trajectories)
/// lower their motion descriptions into certified tubes through this
/// trait, keeping the dependency direction downward.
pub trait LowerToMotorTube {
    /// Lower into a tube over `domain` using `segments` pieces of
    /// Taylor order `order`.
    fn lower_to_motor_tube(
        &self,
        domain: Interval,
        order: usize,
        segments: usize,
    ) -> Result<CertifiedMotorTube, MotionError>;
}
