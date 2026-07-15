//! Moving geometry: a base chart plus a body-to-world motor tube.
//!
//! A moving object does NOT implement the timeless `Chart` contract.
//! [`SpacetimeChart`] exposes `snapshot(t)` (an immutable frozen-time
//! view that DOES implement `Chart`, with deliberately weakened
//! claims) and `eval_over` (a certified field enclosure along a time
//! span, currently for exact-distance base charts only).

use crate::MotionError;
use crate::algebra::{TmMv, homogeneous_point, point_to_mv};
use crate::tube::{CertifiedMotorTube, EnclosureClass};
use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::Cx;
use fs_ga::{Motor, Point as GaPoint};
use fs_geom::{Aabb, BettiBounds, Chart, ChartSample, Point3, TraceStepClaim};
use fs_ivl::Interval;

/// A certified field enclosure over a time span.
#[derive(Debug, Clone, Copy)]
pub struct FieldEnclosure {
    /// Enclosure of the base field value along the pulled-back
    /// trajectory of the query point.
    pub value: Interval,
    /// Enclosure strength.
    pub class: EnclosureClass,
    /// Versor-defect bound of the tube segments involved.
    pub defect: f64,
}

/// A base chart moving through space along a body-to-world tube.
#[derive(Debug)]
pub struct SpacetimeChart<C> {
    base: C,
    tube: CertifiedMotorTube,
}

impl<C: Chart> SpacetimeChart<C> {
    /// Bind a base chart to a body-to-world motor tube.
    #[must_use]
    pub fn new(base: C, tube: CertifiedMotorTube) -> Self {
        SpacetimeChart { base, tube }
    }

    /// The base chart.
    #[must_use]
    pub fn base(&self) -> &C {
        &self.base
    }

    /// The body-to-world tube.
    #[must_use]
    pub fn tube(&self) -> &CertifiedMotorTube {
        &self.tube
    }

    /// Freeze the motion at `t`: an immutable snapshot implementing
    /// `Chart` with time and path provenance recorded and claims
    /// deliberately weakened (see the crate contract).
    pub fn snapshot(&self, t: f64, cx: &Cx<'_>) -> Result<MotionSnapshot<'_, C>, MotionError> {
        let (world_from_body, sample) = self.tube.path().motor_at(t)?;
        let body_from_world = world_from_body.reverse();
        // Transport the support box; an unbounded base support stays
        // unbounded rather than being silently clipped.
        let base_support = self.base.support();
        let finite = [
            base_support.min.x,
            base_support.min.y,
            base_support.min.z,
            base_support.max.x,
            base_support.max.y,
            base_support.max.z,
        ]
        .iter()
        .all(|v| v.is_finite());
        let support = if finite {
            self.tube
                .box_action_over(&base_support, Interval::point(t), cx)?
                .bounds
        } else {
            Aabb::WHOLE_SPACE
        };
        Ok(MotionSnapshot {
            base: &self.base,
            body_from_world,
            support,
            time: t,
            defect: sample.defect,
        })
    }

    /// Certified enclosure of the base field value at world point `x`
    /// over the whole time span. Requires the base chart to claim
    /// `TraceStepClaim::ExactDistance` (exact signed distance is
    /// globally 1-Lipschitz, which turns the pulled-back position
    /// enclosure radius into a field enclosure radius); every other
    /// claim refuses.
    pub fn eval_over(
        &self,
        x: Point3,
        span: Interval,
        cx: &Cx<'_>,
    ) -> Result<FieldEnclosure, MotionError> {
        if self.base.trace_step_claim() != TraceStepClaim::ExactDistance {
            return Err(MotionError::UnsupportedBaseClaim);
        }
        if !(x.x.is_finite() && x.y.is_finite() && x.z.is_finite()) {
            return Err(MotionError::NonFiniteInput { what: "point" });
        }
        // Pull-back enclosure Q ∋ M(t)⁻¹·x for all t in span. The
        // sandwich with the reversed components plus homogeneous
        // division IS the inverse action for versor directions
        // (uniform scaling cancels in the weight).
        let mut coords: Option<[Interval; 3]> = None;
        let mut defect = 0.0f64;
        let domain = self.tube.domain();
        if !domain.encloses(span) {
            return Err(MotionError::OutOfDomain {
                lo: span.lo(),
                hi: span.hi(),
                domain_lo: domain.lo(),
                domain_hi: domain.hi(),
            });
        }
        for seg in self.tube.segments() {
            let d = seg.domain();
            let lo = span.lo().max(d.lo());
            let hi = span.hi().min(d.hi());
            if lo > hi {
                continue;
            }
            cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
            let sub = Interval::new(lo, hi);
            let mv = seg.components();
            let p = TmMv::constant(&point_to_mv(x.x, x.y, x.z), mv.domain(), mv.order())?;
            let sandwich = mv.reverse()?.gp(&p)?.gp(mv)?;
            let enc = sandwich.eval_all(sub)?;
            let pt = homogeneous_point(&enc)?;
            defect = defect.max(seg.defect());
            coords = Some(match coords {
                None => pt,
                Some(prev) => [prev[0].hull(pt[0]), prev[1].hull(pt[1]), prev[2].hull(pt[2])],
            });
        }
        let coords = coords.ok_or(MotionError::EmptyTimeDomain)?;
        // Base sample at the box midpoint; its certificate must be a
        // rigorous enclosure of the signed distance.
        let mid = Point3::new(
            coords[0].midpoint(),
            coords[1].midpoint(),
            coords[2].midpoint(),
        );
        let sample = self.base.eval(mid, cx);
        let value_enclosure = match sample.error.kind {
            NumericalKind::Exact | NumericalKind::Enclosure => {
                Interval::new(sample.error.lo, sample.error.hi)
            }
            NumericalKind::Estimate | NumericalKind::NoClaim => {
                return Err(MotionError::UncertifiedBaseSample);
            }
        };
        // Outward radius of the pull-back box around its midpoint.
        let mut r2 = Interval::point(0.0);
        for c in coords {
            let dev = (c - Interval::point(c.midpoint())).abs_bound();
            r2 = r2 + dev * dev;
        }
        let r = r2.sqrt().hi();
        let value = value_enclosure + Interval::new(-r, r);
        Ok(FieldEnclosure {
            value,
            class: EnclosureClass::Certified,
            defect,
        })
    }
}

/// An immutable frozen-time view of a moving chart. Provenance (time,
/// defect) is recorded; ray-stepping claims, gradients, Lipschitz
/// data, and sample-error certificates are deliberately dropped (see
/// the crate contract's no-claim boundaries).
#[derive(Debug)]
pub struct MotionSnapshot<'a, C> {
    base: &'a C,
    body_from_world: Motor,
    support: Aabb,
    time: f64,
    defect: f64,
}

impl<C> MotionSnapshot<'_, C> {
    /// The frozen time.
    #[must_use]
    pub fn time(&self) -> f64 {
        self.time
    }

    /// The tube's versor-defect bound at the frozen segment.
    #[must_use]
    pub fn defect(&self) -> f64 {
        self.defect
    }
}

impl<C: Chart> Chart for MotionSnapshot<'_, C> {
    fn name(&self) -> &'static str {
        "motion-snapshot"
    }

    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        let pulled = self.body_from_world.transform_point(GaPoint {
            x: x.x,
            y: x.y,
            z: x.z,
        });
        match pulled {
            Ok(q) => {
                let base = self.base.eval(Point3::new(q.x, q.y, q.z), cx);
                ChartSample {
                    signed_distance: base.signed_distance,
                    gradient: None,
                    lipschitz: None,
                    error: NumericalCertificate::no_claim(),
                }
            }
            // Unreachable for finite points under near-unit motors;
            // still fail closed rather than fabricate a value.
            Err(_) => ChartSample {
                signed_distance: f64::NAN,
                gradient: None,
                lipschitz: None,
                error: NumericalCertificate::no_claim(),
            },
        }
    }

    fn support(&self) -> Aabb {
        self.support
    }

    fn topology_hint(&self) -> BettiBounds {
        // An invertible near-rigid pull-back is a homeomorphism of the
        // zero set, so the base bounds transport.
        self.base.topology_hint()
    }
}
