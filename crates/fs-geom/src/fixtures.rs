//! Analytic fixture charts: the shared test vocabulary for ALL of MORPH
//! (this bead's testing mandate). Exact signed distances with unit
//! Lipschitz constants — the ground truth every representation bead
//! (rep-sdf, rep-mesh, router) measures itself against. PUBLIC on purpose:
//! downstream conformance suites import these instead of inventing their
//! own spheres.

use crate::{
    Aabb, BettiBounds, Chart, ChartSample, Differentiability, Point3, TraceStepClaim, Vec3,
};
use fs_evidence::NumericalCertificate;
use fs_exec::Cx;

/// Exact sphere SDF: `|x - c| - r`.
#[derive(Debug, Clone, Copy)]
pub struct SphereChart {
    /// Center.
    pub center: Point3,
    /// Radius (> 0).
    pub radius: f64,
}

impl Chart for SphereChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let d = x.delta_from(self.center);
        let dist = d.norm();
        let gradient = if dist > 1e-12 {
            Some(d.scale(1.0 / dist))
        } else {
            None // the center is the medial axis: no gradient claim
        };
        ChartSample {
            signed_distance: dist - self.radius,
            gradient,
            lipschitz: Some(1.0),
            error: NumericalCertificate::exact(dist - self.radius),
        }
    }

    fn support(&self) -> Aabb {
        let r = self.radius;
        Aabb::new(
            self.center.offset(Vec3::new(-r, -r, -r)),
            self.center.offset(Vec3::new(r, r, r)),
        )
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        if self.radius.is_finite() && self.radius > 0.0 {
            TraceStepClaim::ExactDistance
        } else {
            TraceStepClaim::NoClaim
        }
    }

    fn topology_hint(&self) -> BettiBounds {
        BettiBounds::exact(1, 0, 1)
    }

    fn name(&self) -> &'static str {
        "fixture/sphere"
    }

    fn differentiability(&self) -> Differentiability {
        Differentiability::Smooth
    }
}

/// Exact axis-aligned box SDF (the standard corner-distance formula).
#[derive(Debug, Clone, Copy)]
pub struct BoxChart {
    /// The box.
    pub aabb: Aabb,
}

impl BoxChart {
    fn is_solid_box(&self) -> bool {
        let min = self.aabb.min;
        let max = self.aabb.max;
        min.x.is_finite()
            && min.y.is_finite()
            && min.z.is_finite()
            && max.x.is_finite()
            && max.y.is_finite()
            && max.z.is_finite()
            && max.x > min.x
            && max.y > min.y
            && max.z > min.z
    }
}

impl Chart for BoxChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let c = Point3::new(
            f64::midpoint(self.aabb.min.x, self.aabb.max.x),
            f64::midpoint(self.aabb.min.y, self.aabb.max.y),
            f64::midpoint(self.aabb.min.z, self.aabb.max.z),
        );
        let h = Vec3::new(
            0.5 * (self.aabb.max.x - self.aabb.min.x),
            0.5 * (self.aabb.max.y - self.aabb.min.y),
            0.5 * (self.aabb.max.z - self.aabb.min.z),
        );
        let p = x.delta_from(c);
        let q = Vec3::new(p.x.abs() - h.x, p.y.abs() - h.y, p.z.abs() - h.z);
        let outside = Vec3::new(q.x.max(0.0), q.y.max(0.0), q.z.max(0.0));
        let sd = outside.norm() + q.x.max(q.y.max(q.z)).min(0.0);
        ChartSample {
            signed_distance: sd,
            gradient: None, // piecewise form; gradients per-face arrive with rep-sdf
            lipschitz: Some(1.0),
            error: NumericalCertificate::exact(sd),
        }
    }

    fn support(&self) -> Aabb {
        self.aabb
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        if self.is_solid_box() {
            TraceStepClaim::ExactDistance
        } else {
            TraceStepClaim::NoClaim
        }
    }

    fn topology_hint(&self) -> BettiBounds {
        if self.is_solid_box() {
            BettiBounds::exact(1, 0, 1)
        } else {
            BettiBounds::unknown()
        }
    }

    fn name(&self) -> &'static str {
        "fixture/box"
    }
}

/// Torus implicit field (major radius `major`, tube radius `minor`, z-axis).
/// It is an exact signed distance only for the ring case `major > minor`.
#[derive(Debug, Clone, Copy)]
pub struct TorusChart {
    /// Center of the torus.
    pub center: Point3,
    /// Major (ring) radius.
    pub major: f64,
    /// Minor (tube) radius.
    pub minor: f64,
}

impl TorusChart {
    fn is_exact_distance(&self) -> bool {
        self.major.is_finite()
            && self.minor.is_finite()
            && self.major > 0.0
            && self.minor > 0.0
            && self.major > self.minor
    }
}

impl Chart for TorusChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let p = x.delta_from(self.center);
        let ring = (p.x * p.x + p.y * p.y).sqrt() - self.major;
        let sd = (ring * ring + p.z * p.z).sqrt() - self.minor;
        ChartSample {
            signed_distance: sd,
            gradient: None, // analytic gradient lands with rep-frep
            lipschitz: Some(1.0),
            error: if self.is_exact_distance() {
                NumericalCertificate::exact(sd)
            } else {
                NumericalCertificate::estimate(sd, sd)
            },
        }
    }

    fn support(&self) -> Aabb {
        let r = self.major + self.minor;
        Aabb::new(
            self.center.offset(Vec3::new(-r, -r, -self.minor)),
            self.center.offset(Vec3::new(r, r, self.minor)),
        )
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        if self.is_exact_distance() {
            TraceStepClaim::ExactDistance
        } else if self.major.is_finite()
            && self.minor.is_finite()
            && self.major > 0.0
            && self.minor > 0.0
        {
            TraceStepClaim::LipschitzImplicit
        } else {
            TraceStepClaim::NoClaim
        }
    }

    fn topology_hint(&self) -> BettiBounds {
        if self.is_exact_distance() {
            BettiBounds::exact(1, 1, 1)
        } else {
            BettiBounds::unknown()
        }
    }

    fn name(&self) -> &'static str {
        "fixture/torus"
    }
}

/// A DELIBERATELY WRONG chart for agreement-detection tests: a sphere
/// whose signed distance is biased by `bias` while its error model LIES
/// (declares Exact). The agreement checker must catch it — a chart that
/// mis-declares its error is exactly the failure mode plan §7.1's
/// "checkable proposition" exists for.
#[derive(Debug, Clone, Copy)]
pub struct LyingSphereChart {
    /// The honest geometry.
    pub sphere: SphereChart,
    /// The undeclared bias added to every signed distance.
    pub bias: f64,
}

impl Chart for LyingSphereChart {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        let honest = self.sphere.eval(x, cx);
        ChartSample {
            signed_distance: honest.signed_distance + self.bias,
            error: NumericalCertificate::exact(honest.signed_distance + self.bias),
            ..honest
        }
    }

    fn support(&self) -> Aabb {
        self.sphere.support().inflate(self.bias.abs())
    }

    fn name(&self) -> &'static str {
        "fixture/lying-sphere"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use asupersync::types::Budget;
    use fs_exec::{CancelGate, ExecMode, StreamKey};

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    #[test]
    fn fixture_sdfs_hit_known_values() {
        with_cx(|cx| {
            let s = SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 2.0,
            };
            assert!((s.eval(Point3::new(3.0, 0.0, 0.0), cx).signed_distance - 1.0).abs() < 1e-12);
            assert!((s.eval(Point3::new(0.0, 1.0, 0.0), cx).signed_distance + 1.0).abs() < 1e-12);
            assert!(s.inside(Point3::new(0.0, 1.0, 0.0), cx));

            let b = BoxChart {
                aabb: Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0)),
            };
            assert!((b.eval(Point3::new(2.0, 0.0, 0.0), cx).signed_distance - 1.0).abs() < 1e-12);
            assert!((b.eval(Point3::new(0.0, 0.0, 0.0), cx).signed_distance + 1.0).abs() < 1e-12);
            assert_eq!(b.trace_step_claim(), TraceStepClaim::ExactDistance);

            let degenerate = BoxChart {
                aabb: Aabb::new(Point3::new(-1.0, -1.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
            };
            assert_eq!(degenerate.trace_step_claim(), TraceStepClaim::NoClaim);
            assert_eq!(degenerate.topology_hint(), BettiBounds::unknown());

            let t = TorusChart {
                center: Point3::new(0.0, 0.0, 0.0),
                major: 3.0,
                minor: 1.0,
            };
            assert!((t.eval(Point3::new(3.0, 0.0, 0.0), cx).signed_distance + 1.0).abs() < 1e-12);
            assert!((t.eval(Point3::new(5.0, 0.0, 0.0), cx).signed_distance - 1.0).abs() < 1e-12);
        });
    }
}
