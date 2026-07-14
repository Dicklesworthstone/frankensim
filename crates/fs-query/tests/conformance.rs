//! fs-query conformance suite (CONTRACT.md: any reimplementation must
//! pass). Multi-chart AGREEMENT for closest point and raycasts,
//! tracer safety vs a dense oracle including tangent rays, offsets and
//! ball-Minkowski exactness, certified separation bounds, the
//! Estimate-authority thickness marcher on graded fixtures with the medial cross-check and
//! the design-lever subgradient, and curvature convergence at the
//! documented order per chart class. JSON-line verdicts; seeded cases
//! carry seeds.

use asupersync::types::Budget;
use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Aabb, Chart, ChartSample, Point3, SamplingDomainError, TraceStepClaim, Vec3};
use fs_query::{
    CurvatureClass, OffsetChart, QueryError, SEPARATION_MAX_CHART_SAMPLES, SeparationScope,
    closest_point, closest_point_clipped, curvature, medial_poles, min_thickness,
    min_thickness_clipped, minkowski_ball, raycast, separation, separation_clipped, thickness_at,
    thickness_at_clipped,
};
use fs_rep_frep::{BoolOp, BoolStyle, FrepBuilder};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-query/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\
         \"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "case {case}: {detail}");
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    fn unit(&mut self) -> f64 {
        ((self.next() >> 11) as f64) / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x9E4,
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

#[derive(Clone, Copy)]
enum MalformedPointKind {
    Value,
    Gradient,
}

struct MalformedPointChart {
    kind: MalformedPointKind,
}

struct OffsetEvidenceChart {
    certificate: NumericalCertificate,
}

impl Chart for MalformedPointChart {
    fn eval(&self, _point: Point3, _cx: &Cx<'_>) -> ChartSample {
        ChartSample {
            signed_distance: if matches!(self.kind, MalformedPointKind::Value) {
                f64::NAN
            } else {
                0.0
            },
            gradient: Some(if matches!(self.kind, MalformedPointKind::Gradient) {
                Vec3::new(f64::INFINITY, 0.0, 0.0)
            } else {
                Vec3::new(1.0, 0.0, 0.0)
            }),
            lipschitz: None,
            error: NumericalCertificate::no_claim(),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "test/malformed-point"
    }
}

impl Chart for OffsetEvidenceChart {
    fn eval(&self, _point: Point3, _cx: &Cx<'_>) -> ChartSample {
        ChartSample {
            signed_distance: 0.0,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: self.certificate,
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "test/offset-evidence"
    }
}

struct HugeNewtonChart;

impl Chart for HugeNewtonChart {
    fn eval(&self, _point: Point3, _cx: &Cx<'_>) -> ChartSample {
        ChartSample {
            signed_distance: f64::MAX,
            gradient: Some(Vec3::new(f64::MIN_POSITIVE, 0.0, 0.0)),
            lipschitz: None,
            error: NumericalCertificate::no_claim(),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "test/huge-newton"
    }
}

struct ExtremeFiniteDifferenceChart;

impl Chart for ExtremeFiniteDifferenceChart {
    fn eval(&self, _point: Point3, _cx: &Cx<'_>) -> ChartSample {
        ChartSample {
            signed_distance: 1.0,
            gradient: None,
            lipschitz: None,
            error: NumericalCertificate::no_claim(),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(
            Point3::new(1.7e308, 1.7e308, 1.7e308),
            Point3::new(f64::MAX, f64::MAX, f64::MAX),
        )
    }

    fn name(&self) -> &'static str {
        "test/extreme-finite-difference"
    }
}

struct ExplosiveCurvatureChart;

impl Chart for ExplosiveCurvatureChart {
    fn eval(&self, point: Point3, _cx: &Cx<'_>) -> ChartSample {
        let at_origin = point.x.to_bits() == 0 && point.y.to_bits() == 0 && point.z.to_bits() == 0;
        let value = if at_origin { 0.0 } else { f64::MAX };
        ChartSample {
            signed_distance: value,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: None,
            error: NumericalCertificate::estimate(value, value),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "test/explosive-curvature"
    }
}

fn frep_sphere(center: Point3, r: f64) -> fs_rep_frep::Frep {
    let mut b = FrepBuilder::new();
    let s = b.sphere(center, r).expect("sphere");
    b.finish(s).expect("frep")
}

/// gq-001 — closest point agrees across chart representations of the
/// SAME sphere within per-chart certificates, residuals are honest,
/// and the answer is translation-equivariant (G3).
#[test]
fn gq_001_closest_point_agreement() {
    with_cx(|cx| {
        let exact = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let frep = frep_sphere(Point3::new(0.0, 0.0, 0.0), 1.0);
        let tiled = fs_rep_sdf::TiledSdf::build(&exact, 0.05, cx).expect("tiled");
        let mesh = fs_rep_mesh::MeshChart::new(fs_rep_mesh::shapes::icosphere(
            Point3::new(0.0, 0.0, 0.0),
            1.0,
            4,
        ));
        // Per-chart tolerance = each chart's OWN certificate: exact and
        // F-rep are exact; the tiled grid declares its interpolation
        // bound; the mesh chart agrees at its faceting scale.
        let charts: Vec<(&dyn Chart, f64)> = vec![
            (&exact, 1e-9),
            (&frep, 1e-9),
            (&tiled, tiled.bound() * 1.05),
            (&mesh, 3e-2),
        ];
        let mut rng = Lcg(0x1001_2026_0707_0011);
        let mut worst = vec![0.0f64; charts.len()];
        for _ in 0..80 {
            // Shell radii keep queries inside every chart's valid
            // domain (the tiled grid ends at its support box).
            let z = rng.range(-1.0, 1.0);
            let sq = (1.0f64 - z * z).max(0.0).sqrt();
            let th = rng.range(0.0, core::f64::consts::TAU);
            let rad = rng.range(0.3, 1.35);
            let p = Point3::new(rad * sq * th.cos(), rad * sq * th.sin(), rad * z);
            let truth = {
                let d = p.delta_from(Point3::new(0.0, 0.0, 0.0));
                let n = d.norm();
                Point3::new(d.x / n, d.y / n, d.z / n)
            };
            for (k, (chart, _)) in charts.iter().enumerate() {
                let cp = closest_point(*chart, p, cx).expect("closest");
                worst[k] = worst[k]
                    .max(cp.point.delta_from(truth).norm())
                    .max(cp.residual);
            }
        }
        let agree = charts
            .iter()
            .enumerate()
            .all(|(k, (_, tol))| worst[k] < *tol);
        // G3: shifted sphere, shifted queries, identical geometry.
        let shifted = SphereChart {
            center: Point3::new(0.5, 0.25, -0.375),
            radius: 1.0,
        };
        let q = Point3::new(1.7, 0.3, -0.2);
        let base = closest_point(&exact, q, cx).expect("base");
        let moved =
            closest_point(&shifted, q.offset(Vec3::new(0.5, 0.25, -0.375)), cx).expect("moved");
        let g3 = moved
            .point
            .delta_from(base.point.offset(Vec3::new(0.5, 0.25, -0.375)))
            .norm()
            < 1e-9;
        verdict(
            "gq-001",
            agree && g3,
            &format!(
                "closest point agrees with the analytic answer across exact/F-rep/\
                 tiled-SDF/mesh charts (worst errors {:?} within per-chart \
                 certificates) with honest residuals, and is translation-equivariant; \
                 seed 0x1001_2026_0707_0011",
                worst.iter().map(|w| format!("{w:.1e}")).collect::<Vec<_>>()
            ),
        );
    });
}

/// gq-001a — closest-point queries refuse malformed producer samples,
/// overflowing Newton updates, overflowing finite-difference points, and
/// producer-requested cancellation instead of publishing NaN/Inf answers.
#[test]
#[allow(clippy::too_many_lines)] // One fail-closed campaign shares the malformed-producer matrix.
fn gq_001a_closest_point_fails_closed_on_nonfinite_paths() {
    let malformed_value = with_cx(|cx| {
        closest_point(
            &MalformedPointChart {
                kind: MalformedPointKind::Value,
            },
            Point3::new(0.0, 0.0, 0.0),
            cx,
        )
    });
    let malformed_gradient = with_cx(|cx| {
        closest_point(
            &MalformedPointChart {
                kind: MalformedPointKind::Gradient,
            },
            Point3::new(0.0, 0.0, 0.0),
            cx,
        )
    });
    let newton_overflow =
        with_cx(|cx| closest_point(&HugeNewtonChart, Point3::new(0.0, 0.0, 0.0), cx));
    let fd_overflow = with_cx(|cx| {
        closest_point(
            &ExtremeFiniteDifferenceChart,
            Point3::new(f64::MAX, f64::MAX, f64::MAX),
            cx,
        )
    });
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let cancelled = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x9E4,
                kernel_id: 11,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        closest_point(
            &CancellingSeparationSphere {
                gate: &gate,
                during: CancelDuring::Eval,
            },
            Point3::new(2.0, 0.0, 0.0),
            &cx,
        )
    });
    let medial_boundary = fs_rep_mesh::Soup {
        positions: vec![
            Point3::new(-1.0, -1.0, -1.0),
            Point3::new(1.0, -1.0, 1.0),
            Point3::new(-1.0, 1.0, 1.0),
            Point3::new(1.0, 1.0, -1.0),
        ],
        triangles: vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
    };
    let medial_gate = CancelGate::new();
    let medial_pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let cancelled_medial = medial_pool.scope(|arena| {
        let cx = Cx::new(
            &medial_gate,
            arena,
            StreamKey {
                seed: 0x9E4,
                kernel_id: 6,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        medial_poles(
            &CancellingSeparationSphere {
                gate: &medial_gate,
                during: CancelDuring::Eval,
            },
            &medial_boundary,
            0.0,
            &cx,
        )
    });
    let invalid_medial_boundary = with_cx(|cx| {
        medial_poles(
            &SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            },
            &fs_rep_mesh::Soup {
                positions: vec![Point3::new(0.0, 0.0, 0.0)],
                triangles: vec![[0, 1, 0]],
            },
            1.0,
            cx,
        )
    });
    let overflowing_medial_threshold = with_cx(|cx| {
        medial_poles(
            &SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            },
            &medial_boundary,
            f64::MAX,
            cx,
        )
    });
    verdict(
        "gq-001a",
        matches!(malformed_value, Err(QueryError::InvalidPointSample { .. }))
            && matches!(
                malformed_gradient,
                Err(QueryError::InvalidPointSample { .. })
            )
            && matches!(
                newton_overflow,
                Err(QueryError::InvalidPointArithmetic { .. })
            )
            && matches!(fd_overflow, Err(QueryError::InvalidPointSample { .. }))
            && matches!(cancelled, Err(QueryError::Cancelled))
            && matches!(cancelled_medial, Err(QueryError::Cancelled))
            && matches!(
                invalid_medial_boundary,
                Err(QueryError::InvalidBoundaryIndex {
                    triangle: 0,
                    corner: 1,
                    index: 1,
                    positions: 1,
                })
            )
            && matches!(
                overflowing_medial_threshold,
                Err(QueryError::InvalidPointArithmetic { .. })
            ),
        "closest point rejects nonfinite producer output, extreme FD coordinates, and \
         overflowing Newton/medial arithmetic; malformed public Soup indices refuse before \
         Delaunay, and cancellation requested inside closest-point or medial-pole eval wins \
         before publication",
    );
}

/// gq-002 — raycast: analytic agreement across chart types, and SAFETY
/// on CSG (no tunneling vs a dense oracle), including tangent rays.
#[test]
fn gq_002_raycast_safety() {
    with_cx(|cx| {
        let exact = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let frep = frep_sphere(Point3::new(0.0, 0.0, 0.0), 1.0);
        let charts: Vec<(&dyn Chart, f64)> = vec![(&exact, 1e-6), (&frep, 1e-6)];
        let mut rng = Lcg(0x1001_2026_0707_0012);
        let mut agree = true;
        for _ in 0..100 {
            let z = rng.range(-1.0, 1.0);
            let s = (1.0 - z * z).max(0.0).sqrt();
            let th = rng.range(0.0, core::f64::consts::TAU);
            let o = Point3::new(3.0 * s * th.cos(), 3.0 * s * th.sin(), 3.0 * z);
            let dir = Point3::new(0.0, 0.0, 0.0).delta_from(o);
            // Analytic first hit: |o| − 1 along the unit direction.
            let expect_t = o.delta_from(Point3::new(0.0, 0.0, 0.0)).norm() - 1.0;
            for (chart, tol) in &charts {
                let hit = raycast(*chart, o, dir, 10.0, cx)
                    .expect("cast")
                    .expect("hits");
                agree &= (hit.t - expect_t).abs() < 1e-3 + tol;
            }
        }
        // Tangent rays: aimed exactly at radius 1 (grazing) and 1.05
        // (clean miss): no tunneling, misses classified.
        let graze_o = Point3::new(-3.0, 1.0, 0.0);
        let graze = raycast(&frep, graze_o, Vec3::new(1.0, 0.0, 0.0), 10.0, cx);
        let graze_ok = match graze {
            Ok(Some(h)) => {
                (h.point.delta_from(Point3::new(0.0, 0.0, 0.0)).norm() - 1.0).abs() < 1e-3
            }
            Err(QueryError::UnresolvedTrace { .. }) => true,
            Ok(None) => false,
            Err(error) => panic!("unexpected grazing-ray failure: {error}"),
        };
        let miss = raycast(
            &frep,
            Point3::new(-3.0, 1.05, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            10.0,
            cx,
        )
        .expect("miss");
        // CSG safety vs dense oracle.
        let csg = {
            let mut b = FrepBuilder::new();
            let s1 = b.sphere(Point3::new(-0.4, 0.0, 0.0), 0.8).expect("s1");
            let s2 = b.sphere(Point3::new(0.4, 0.0, 0.0), 0.8).expect("s2");
            let u = b
                .boolean(BoolOp::Difference, BoolStyle::Blend { radius: 0.2 }, s1, s2)
                .expect("u");
            b.finish(u).expect("frep")
        };
        let mut safety = true;
        for _ in 0..200 {
            let z = rng.range(-1.0, 1.0);
            let s = (1.0 - z * z).max(0.0).sqrt();
            let th = rng.range(0.0, core::f64::consts::TAU);
            let o = Point3::new(2.5 * s * th.cos(), 2.5 * s * th.sin(), 2.5 * z);
            let target = Point3::new(
                rng.range(-0.5, 0.5),
                rng.range(-0.5, 0.5),
                rng.range(-0.5, 0.5),
            );
            let dir = target.delta_from(o);
            let dn = dir.norm();
            let d = dir.scale(1.0 / dn);
            let (hit, clean_miss) = match raycast(&csg, o, d, 6.0, cx) {
                Ok(Some(hit)) => (Some(hit), false),
                Ok(None) => (None, true),
                Err(QueryError::UnresolvedTrace { .. }) => (None, false),
                Err(error) => panic!("unexpected CSG raycast failure: {error}"),
            };
            // Dense oracle: first sign change.
            let mut oracle = None;
            let mut prev = csg.value(o);
            for i in 1..=1200 {
                let t = 6.0 * f64::from(i) / 1200.0;
                let v = csg.value(o.offset(d.scale(t)));
                if prev >= 0.0 && v < 0.0 {
                    oracle = Some(t);
                    break;
                }
                prev = v;
            }
            if let (Some(t_true), Some(h)) = (oracle, hit) {
                safety &= h.t <= t_true + 1e-3;
            }
            if oracle.is_some() && clean_miss {
                safety = false;
            }
            // Oracle-hit + explicit UnresolvedTrace is incomplete, not unsafe;
            // an oracle hit plus a clean certified miss is a failure.
        }
        verdict(
            "gq-002",
            agree && graze_ok && miss.is_none() && safety,
            "raycasts match the analytic sphere across chart types, tangent rays \
             never tunnel (grazes land on the surface or approach; the 1.05 offset \
             misses cleanly), and the CSG tracer never claims a hit past the dense \
             oracle over 200 rays; seed 0x1001_2026_0707_0012",
        );
    });
}

#[derive(Debug, Clone, Copy)]
enum EndpointFailure {
    MissingLipschitz,
    UncertifiedValue,
    ZeroStraddlingValue,
}

struct EndpointFailureChart {
    failure: EndpointFailure,
}

impl Chart for EndpointFailureChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let at_endpoint = x.x >= 0.5;
        ChartSample {
            signed_distance: if at_endpoint
                && matches!(self.failure, EndpointFailure::ZeroStraddlingValue)
            {
                1e-8
            } else {
                1.0
            },
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: if at_endpoint && matches!(self.failure, EndpointFailure::MissingLipschitz) {
                None
            } else {
                Some(1.0)
            },
            error: NumericalCertificate::estimate(1.0, 1.0),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-2.0, -2.0, -2.0), Point3::new(2.0, 2.0, 2.0))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::LipschitzImplicit
    }

    fn trace_value_enclosure(
        &self,
        x: Point3,
        _sample: &ChartSample,
        _cx: &Cx<'_>,
    ) -> NumericalCertificate {
        if x.x >= 0.5 && matches!(self.failure, EndpointFailure::UncertifiedValue) {
            NumericalCertificate::estimate(1.0, 1.0)
        } else if x.x >= 0.5 && matches!(self.failure, EndpointFailure::ZeroStraddlingValue) {
            NumericalCertificate::enclosure(-1e-8, 1e-8)
        } else {
            NumericalCertificate::exact(1.0)
        }
    }

    fn name(&self) -> &'static str {
        "test/endpoint-failure"
    }
}

/// gq-002c — the generic tracer validates every point, including the caller's
/// bounded endpoint. Missing local evidence and an Estimate cannot be silently
/// inherited from the origin or reported as a clean miss.
#[test]
fn gq_002c_raycast_validates_each_sample_and_tmax() {
    with_cx(|cx| {
        let cast = |failure| {
            raycast(
                &EndpointFailureChart { failure },
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                0.5,
                cx,
            )
        };
        let missing = cast(EndpointFailure::MissingLipschitz);
        let uncertified = cast(EndpointFailure::UncertifiedValue);
        let straddling = cast(EndpointFailure::ZeroStraddlingValue);
        let subnormal_direction = raycast(
            &EndpointFailureChart {
                failure: EndpointFailure::UncertifiedValue,
            },
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(f64::from_bits(1), 0.0, 0.0),
            0.5,
            cx,
        );
        verdict(
            "gq-002c",
            matches!(missing, Err(QueryError::NoLipschitz))
                && matches!(uncertified, Err(QueryError::InvalidTraceSample { .. }))
                && matches!(straddling, Err(QueryError::UnresolvedTrace { .. }))
                && matches!(
                    subnormal_direction,
                    Err(QueryError::InvalidTraceSample { .. })
                ),
            "raycast revalidates the local Lipschitz bound and Exact/Enclosure trace \
             evidence at tmax; a zero-straddling endpoint remains unresolved, and finite \
             nonzero subnormal directions normalize without overflow",
        );
    });
}

#[derive(Debug, Clone, Copy)]
enum CancelDuring {
    Eval,
    TraceEnclosure,
}

struct CancellingTraceChart<'a> {
    gate: &'a CancelGate,
    during: CancelDuring,
}

impl Chart for CancellingTraceChart<'_> {
    fn eval(&self, _x: Point3, _cx: &Cx<'_>) -> ChartSample {
        if matches!(self.during, CancelDuring::Eval) {
            self.gate.request();
        }
        ChartSample {
            signed_distance: 0.0,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: NumericalCertificate::exact(0.0),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn trace_value_enclosure(
        &self,
        _x: Point3,
        _sample: &ChartSample,
        _cx: &Cx<'_>,
    ) -> NumericalCertificate {
        if matches!(self.during, CancelDuring::TraceEnclosure) {
            self.gate.request();
        }
        NumericalCertificate::exact(0.0)
    }

    fn name(&self) -> &'static str {
        "test/cancelling-trace"
    }
}

/// gq-002d — cancellation requested inside a chart producer wins over that
/// same call's otherwise-valid hit evidence.
#[test]
fn gq_002d_raycast_rechecks_cancellation_after_chart_calls() {
    let run = |during| {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x9E4,
                    kernel_id: 2,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            raycast(
                &CancellingTraceChart {
                    gate: &gate,
                    during,
                },
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                1.0,
                &cx,
            )
        })
    };
    let after_eval = run(CancelDuring::Eval);
    let after_trace = run(CancelDuring::TraceEnclosure);
    verdict(
        "gq-002d",
        matches!(after_eval, Err(QueryError::Cancelled))
            && matches!(after_trace, Err(QueryError::Cancelled)),
        "cancellation requested during eval or trace_value_enclosure is observed before \
         either producer can authorize a hit",
    );
}

struct LooseLipschitzChart;

impl Chart for LooseLipschitzChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let value = 1e-12 * (x.x - 1.0);
        ChartSample {
            signed_distance: value,
            gradient: Some(Vec3::new(1e-12, 0.0, 0.0)),
            // Valid but deliberately loose: the true field Lipschitz constant
            // is 1e-12, so 1.0 remains an upper bound.
            lipschitz: Some(1.0),
            error: NumericalCertificate::estimate(value, value),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-2.0, -1.0, -1.0), Point3::new(2.0, 1.0, 1.0))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::LipschitzImplicit
    }

    fn trace_value_enclosure(
        &self,
        _x: Point3,
        sample: &ChartSample,
        _cx: &Cx<'_>,
    ) -> NumericalCertificate {
        NumericalCertificate::enclosure(
            sample.signed_distance.next_down(),
            sample.signed_distance.next_up(),
        )
    }

    fn name(&self) -> &'static str {
        "test/loose-lipschitz"
    }
}

/// gq-002e — `|f|/L` is a safe-step lower bound, not a proximity upper
/// bound. A loose valid L must not turn a point one world unit from the zero
/// set into an immediate geometric hit.
#[test]
fn gq_002e_lipschitz_implicit_residual_cannot_authorize_hit() {
    with_cx(|cx| {
        let result = raycast(
            &LooseLipschitzChart,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            1.0,
            cx,
        );
        verdict(
            "gq-002e",
            matches!(result, Err(QueryError::UnresolvedTrace { .. })),
            "a valid loose Lipschitz upper bound supports conservative marching but \
             cannot promote a 1e-12 normalized residual one unit from zero into RayHit",
        );
    });
}

/// gq-002b — raycast fails CLOSED on a chart that reports a Lipschitz value but
/// makes no tunneling-safe trace claim. Regression: raycast admitted any
/// `Some(lipschitz)` chart and stepped by φ/L, so an enclosure/heuristic chart
/// (dense SDF, mesh — `Some(lipschitz)` but `NoClaim`) whose reported distance
/// overshoots the true one would tunnel through the surface. Exact and
/// Lipschitz-implicit charts still march conservatively (gq-002/e).
#[test]
fn gq_002b_raycast_refuses_no_claim_charts() {
    with_cx(|cx| {
        // A degenerate SphereChart (radius 0) reports `lipschitz: Some(1.0)` but
        // `trace_step_claim() == NoClaim` — exactly the Some(lipschitz)+NoClaim
        // shape that TiledSdf / MeshChart present to a generic tracer.
        let no_claim = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 0.0,
        };
        assert_eq!(
            no_claim.trace_step_claim(),
            fs_geom::TraceStepClaim::NoClaim,
            "fixture precondition: radius-0 sphere makes no trace claim"
        );
        assert!(
            no_claim
                .eval(Point3::new(3.0, 0.0, 0.0), cx)
                .lipschitz
                .is_some(),
            "fixture precondition: it still reports a Lipschitz value"
        );
        let r = raycast(
            &no_claim,
            Point3::new(-3.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            10.0,
            cx,
        );
        verdict(
            "gq-002b",
            matches!(r, Err(QueryError::NoTraceClaim)),
            "raycast fails closed (NoTraceClaim) on a Some(lipschitz)+NoClaim chart \
             instead of tunneling; typed charts retain their certified march path",
        );
    });
}

/// gq-003 — offsets and the ball-Minkowski identity: offset spheres
/// are spheres (across chart types), erosion shrinks, and
/// minkowski_ball IS the offset (exact by construction).
#[test]
fn gq_003_offset_minkowski() {
    with_cx(|cx| {
        let exact = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let frep = frep_sphere(Point3::new(0.0, 0.0, 0.0), 1.0);
        let mut rng = Lcg(0x1001_2026_0707_0013);
        let mut ok = true;
        for (chart, tol) in [(&exact as &dyn Chart, 1e-12), (&frep as &dyn Chart, 1e-12)] {
            let grown = OffsetChart::new(chart, 0.3).expect("finite dilation");
            let eroded = OffsetChart::new(chart, -0.2).expect("finite erosion");
            let mink = minkowski_ball(chart, 0.3).expect("finite ball radius");
            for _ in 0..100 {
                let p = Point3::new(
                    rng.range(-2.0, 2.0),
                    rng.range(-2.0, 2.0),
                    rng.range(-2.0, 2.0),
                );
                let d = p.delta_from(Point3::new(0.0, 0.0, 0.0)).norm();
                ok &= (grown.eval(p, cx).signed_distance - (d - 1.3)).abs() < tol;
                ok &= (eroded.eval(p, cx).signed_distance - (d - 0.8)).abs() < tol;
                ok &= grown.eval(p, cx).signed_distance.to_bits()
                    == mink.eval(p, cx).signed_distance.to_bits();
            }
        }
        // Offset charts retain differential queries: closest point on the
        // grown sphere lands at radius 1.3. Generic raycast remains an honest
        // NoClaim without a reach/proximity theorem for the offset level set.
        let grown = OffsetChart::new(&exact, 0.3).expect("finite dilation");
        let transformed = grown.eval(Point3::new(2.0, 0.0, 0.0), cx);
        let transformed_authority = transformed.error.kind == NumericalKind::Estimate
            && transformed.error.lo <= transformed.signed_distance
            && transformed.signed_distance <= transformed.error.hi;
        let banded = OffsetChart::new(
            &OffsetEvidenceChart {
                certificate: NumericalCertificate::estimate(-100.0, 100.0),
            },
            0.3,
        )
        .expect("finite radius")
        .eval(Point3::new(0.0, 0.0, 0.0), cx);
        let preserved_band = banded.error.kind == NumericalKind::Estimate
            && banded.error.lo <= -100.3
            && banded.error.hi >= 99.7
            && banded.error.lo < banded.signed_distance
            && banded.signed_distance < banded.error.hi;
        let overflowing_band = OffsetChart::new(
            &OffsetEvidenceChart {
                certificate: NumericalCertificate::estimate(-f64::MAX, f64::MAX),
            },
            f64::MAX,
        )
        .expect("maximum finite radius")
        .eval(Point3::new(0.0, 0.0, 0.0), cx);
        let malformed = OffsetChart::new(
            &OffsetEvidenceChart {
                certificate: NumericalCertificate {
                    kind: NumericalKind::Exact,
                    lo: -1.0,
                    hi: 1.0,
                },
            },
            0.3,
        )
        .expect("finite radius")
        .eval(Point3::new(0.0, 0.0, 0.0), cx);
        let invalid_radii = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY]
            .into_iter()
            .all(|radius| {
                matches!(
                    OffsetChart::new(&exact, radius),
                    Err(QueryError::InvalidOffsetRadius { radius_bits })
                        if radius_bits == radius.to_bits()
                )
            });
        let cp = closest_point(&grown, Point3::new(2.0, 0.5, -0.3), cx).expect("cp");
        let on_grown = (cp.point.delta_from(Point3::new(0.0, 0.0, 0.0)).norm() - 1.3).abs() < 1e-9;
        verdict(
            "gq-003",
            ok && on_grown
                && transformed_authority
                && preserved_band
                && overflowing_band.error.kind == NumericalKind::NoClaim
                && malformed.error.kind == NumericalKind::NoClaim
                && invalid_radii,
            "offset spheres are exactly spheres of the summed radius across chart \
             types, erosion shrinks exactly, minkowski_ball is BITWISE the offset \
             chart, valid inner evidence is outward-translated into an Estimate preserving the \
             full band and containing the new nominal, malformed or overflowing evidence and \
             nonfinite radii fail closed, and offset charts \
             retain closest-point queries; \
             seed 0x1001_2026_0707_0013",
        );
    });
}

/// gq-004 — certified separation: the rigorous lower bound brackets
/// the analytic separation of two spheres, tracks as they approach,
/// and the clearance field dominates the separation everywhere (G0).
#[test]
fn gq_004_separation_certified() {
    with_cx(|cx| {
        let mut rng = Lcg(0x1001_2026_0707_0014);
        let mut ok = true;
        let mut worst_gap = 0.0f64;
        for gap in [1.0, 0.5, 0.2, 0.05] {
            let a = SphereChart {
                center: Point3::new(-(1.0 + gap / 2.0), 0.0, 0.0),
                radius: 1.0,
            };
            let b = SphereChart {
                center: Point3::new(1.0 + gap / 2.0, 0.0, 0.0),
                radius: 1.0,
            };
            let sep = separation(&a, &b, 24, cx).expect("separation");
            ok &= sep.lower_bound <= gap + 1e-9 && sep.observed >= gap - 1e-9;
            ok &= sep.observed - sep.lower_bound < 0.6; // slack is finite and stated
            worst_gap = worst_gap.max(sep.observed - gap);
            // Field law: c(p) ≥ separation at random points.
            let field = fs_query::ClearanceField { a: &a, b: &b };
            for _ in 0..50 {
                let p = Point3::new(
                    rng.range(-3.0, 3.0),
                    rng.range(-2.0, 2.0),
                    rng.range(-2.0, 2.0),
                );
                ok &= field.value(p, cx) >= gap - 1e-9;
            }
        }
        verdict(
            "gq-004",
            ok,
            &format!(
                "separation brackets hold across gaps 1.0 -> 0.05 (true separation in \
                 [lower_bound, observed], observed within {worst_gap:.1e} above \
                 truth), and the clearance field dominates the separation everywhere; \
                 seed 0x1001_2026_0707_0014"
            ),
        );
    });
}

/// gq-004a — neither plausible local Lipschitz/enclosure fields nor malformed
/// per-sample evidence can mint a rigorous separation bracket.
#[test]
fn gq_004a_separation_requires_global_exact_distance_authority() {
    with_cx(|cx| {
        let exact = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let forged = separation(&ForgedLocalSeparationChart, &exact, 2, cx);
        let estimate = separation(
            &MalformedSeparationChart {
                certificate: NumericalCertificate::estimate(-1.0, 1.0),
            },
            &exact,
            2,
            cx,
        );
        let nonfinite = separation(
            &MalformedSeparationChart {
                certificate: NumericalCertificate {
                    kind: NumericalKind::Enclosure,
                    lo: f64::NAN,
                    hi: 1.0,
                },
            },
            &exact,
            2,
            cx,
        );
        verdict(
            "gq-004a",
            matches!(
                forged,
                Err(QueryError::SeparationRequiresExactDistance {
                    input: "a",
                    claim: TraceStepClaim::NoClaim,
                })
            ) && matches!(estimate, Err(QueryError::InvalidTraceSample { .. }))
                && matches!(nonfinite, Err(QueryError::InvalidTraceSample { .. })),
            "local Lipschitz/enclosure fields do not upgrade NoClaim, and ExactDistance \
             inputs must still supply finite rigorous per-sample trace enclosures",
        );
    });
}

struct CancellingSeparationSphere<'a> {
    gate: &'a CancelGate,
    during: CancelDuring,
}

impl Chart for CancellingSeparationSphere<'_> {
    fn eval(&self, point: Point3, cx: &Cx<'_>) -> ChartSample {
        let sample = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        }
        .eval(point, cx);
        if matches!(self.during, CancelDuring::Eval) {
            self.gate.request();
        }
        sample
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn trace_value_enclosure(
        &self,
        _point: Point3,
        sample: &ChartSample,
        _cx: &Cx<'_>,
    ) -> NumericalCertificate {
        if matches!(self.during, CancelDuring::TraceEnclosure) {
            self.gate.request();
        }
        sample.error
    }

    fn name(&self) -> &'static str {
        "test/cancelling-separation-sphere"
    }
}

/// gq-004b — cancellation requested by either separation producer wins
/// immediately over that producer's otherwise-rigorous evidence.
#[test]
fn gq_004b_separation_rechecks_producer_cancellation() {
    let run = |during| {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x9E4,
                    kernel_id: 4,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            let chart = CancellingSeparationSphere {
                gate: &gate,
                during,
            };
            separation(&chart, &chart, 2, &cx)
        })
    };
    verdict(
        "gq-004b",
        matches!(run(CancelDuring::Eval), Err(QueryError::Cancelled))
            && matches!(
                run(CancelDuring::TraceEnclosure),
                Err(QueryError::Cancelled)
            ),
        "separation checkpoints directly after eval and trace_value_enclosure, so \
         producer-requested cancellation wins before bracket authority",
    );
}

/// gq-005 — the thickness estimator agrees with a graded slab, finds the
/// dumbbell neck, cross-checks against medial poles, and responds to a
/// DESIGN LEVER with the right subgradient (differentiable-friendly).
#[test]
#[allow(clippy::too_many_lines)] // slab, dumbbell, medial, and lever are one story
fn gq_005_thickness_oracle() {
    with_cx(|cx| {
        // Graded slab via F-rep: |z| ≤ (t0 + g·x)/2 within a box.
        // Analytic thickness at (x, y, 0-top): t0 + g·x.
        let slab = |t0: f64, g: f64| -> fs_rep_frep::Frep {
            let mut b = FrepBuilder::new();
            let bx = b
                .box_prim(Point3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 1.0, 1.5))
                .expect("box");
            // Half-space pair z ≤ (t0+g·x)/2 and −z ≤ (t0+g·x)/2 via
            // rotated half-spaces: n·p ≤ d with n = (−g/2, 0, 1)/|·|.
            let nz = 1.0f64;
            let nx = -g / 2.0;
            let nn = (nx * nx + nz * nz).sqrt();
            let top = b
                .half_space(Vec3::new(nx / nn, 0.0, nz / nn), t0 / (2.0 * nn))
                .expect("top");
            let bot = b
                .half_space(Vec3::new(nx / nn, 0.0, -nz / nn), t0 / (2.0 * nn))
                .expect("bot");
            let both = b
                .boolean(BoolOp::Intersect, BoolStyle::Hard, top, bot)
                .expect("b");
            let root = b
                .boolean(BoolOp::Intersect, BoolStyle::Hard, bx, both)
                .expect("r");
            b.finish(root).expect("frep")
        };
        let s = slab(0.4, 0.1);
        let mut ok = true;
        for x in [-1.0, 0.0, 1.0] {
            let expect = 0.4 + 0.1 * x;
            // Top surface point at height (t0+g·x)/2... solve: z = d(x).
            let z = expect / 2.0;
            // Project to be exactly on the boundary first.
            let cp = closest_point(&s, Point3::new(x, 0.0, z + 1e-3), cx).expect("cp");
            let t = thickness_at(&s, cp.point, cx).expect("thickness");
            // The inward-normal chord of a wedge is thickness/cos(tilt);
            // tilt is atan(g/2) — tiny; accept 1% relative.
            ok &=
                (t.value - expect).abs() / expect < 0.01 && t.authority == NumericalKind::Estimate;
        }
        // Dumbbell: two balls joined by a thin neck (hard union).
        let dumbbell = |neck_r: f64| -> fs_rep_frep::Frep {
            let mut b = FrepBuilder::new();
            let s1 = b.sphere(Point3::new(-1.2, 0.0, 0.0), 0.8).expect("s1");
            let s2 = b.sphere(Point3::new(1.2, 0.0, 0.0), 0.8).expect("s2");
            let neck = b
                .cylinder(Point3::new(0.0, 0.0, 0.0), neck_r)
                .expect("neck");
            // Cylinder is along z; rotate to x: rotate about y by 90°.
            let neck = b
                .rotate(neck, Vec3::new(0.0, 1.0, 0.0), core::f64::consts::FRAC_PI_2)
                .expect("rot");
            // Bound the infinite cylinder to the joint region.
            let span = b
                .box_prim(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.2, 0.5, 0.5))
                .expect("span");
            let neck = b
                .boolean(BoolOp::Intersect, BoolStyle::Hard, neck, span)
                .expect("n");
            let uni = b
                .boolean(BoolOp::Union, BoolStyle::Hard, s1, s2)
                .expect("u");
            let root = b
                .boolean(BoolOp::Union, BoolStyle::Hard, uni, neck)
                .expect("root");
            b.finish(root).expect("frep")
        };
        let d = dumbbell(0.15);
        // Boundary samples on the neck: points at radius 0.15 around x=0.
        let mut samples = Vec::new();
        for k in 0..16 {
            let th = core::f64::consts::TAU * f64::from(k) / 16.0;
            samples.push(Point3::new(0.0, 0.15 * th.cos(), 0.15 * th.sin()));
        }
        let minimum = min_thickness(&d, &samples, cx).expect("min thickness");
        let min_t = minimum.value;
        let skipped = minimum.skipped;
        let neck_ok = (min_t - 0.3).abs() < 0.01
            && skipped == 0
            && minimum.authority == NumericalKind::Estimate;
        // Medial cross-check on the slab: poles' 2r matches thickness.
        let (hull, _) =
            fs_rep_mesh::dual_contour(&s, fs_rep_mesh::DcOptions::sharp(0.1), cx).expect("dc");
        let poles = medial_poles(&s, &hull, 1.2, cx).expect("poles");
        let mid_pole = poles
            .iter()
            .filter(|(p, _)| p.x.abs() < 0.4 && p.y.abs() < 0.4)
            .map(|(_, r)| 2.0 * r)
            .fold(f64::INFINITY, f64::min);
        let medial_agrees = (mid_pole - 0.4).abs() < 0.08;
        // Design-lever subgradient: d(min neck thickness)/d(neck_r) ≈ 2.
        let h = 1e-4;
        let t_hi = min_thickness(
            &dumbbell(0.15 + h),
            &{
                let mut v = Vec::new();
                for k in 0..16 {
                    let th = core::f64::consts::TAU * f64::from(k) / 16.0;
                    v.push(Point3::new(
                        0.0,
                        (0.15 + h) * th.cos(),
                        (0.15 + h) * th.sin(),
                    ));
                }
                v
            },
            cx,
        )
        .expect("hi")
        .value;
        let t_lo = min_thickness(
            &dumbbell(0.15 - h),
            &{
                let mut v = Vec::new();
                for k in 0..16 {
                    let th = core::f64::consts::TAU * f64::from(k) / 16.0;
                    v.push(Point3::new(
                        0.0,
                        (0.15 - h) * th.cos(),
                        (0.15 - h) * th.sin(),
                    ));
                }
                v
            },
            cx,
        )
        .expect("lo")
        .value;
        let subgrad = (t_hi - t_lo) / (2.0 * h);
        let lever_ok = (subgrad - 2.0).abs() < 1e-3;
        let mut em = fs_obs::Emitter::new("fs-query/conformance", "gq-005/thickness");
        let line = em
            .emit(
                fs_obs::Severity::Info,
                fs_obs::EventKind::Custom {
                    name: "query-thickness-oracle".to_string(),
                    json: format!(
                        "{{\"slab_ok\":{ok},\"neck_min\":{min_t:.4},\"medial_2r\":{mid_pole:.4},\
                         \"lever_subgradient\":{subgrad:.4}}}"
                    ),
                },
                None,
            )
            .to_jsonl();
        fs_obs::validate_line(&line).expect("thickness event validates");
        println!("{line}");
        verdict(
            "gq-005",
            ok && neck_ok && medial_agrees && lever_ok,
            &format!(
                "the graded slab reads its analytic thickness at three stations \
                 (1% rel), the dumbbell neck minimum is 2x the neck radius \
                 ({min_t:.3} vs 0.300, 0 skipped), medial poles cross-check the slab \
                 core (2r = {mid_pole:.3} vs 0.4), and the design-lever subgradient \
                 is {subgrad:.4} (analytic 2) — differentiable-friendly, demonstrated"
            ),
        );
    });
}

struct NonfiniteThicknessChart;

struct UnrepresentableThicknessSlab {
    min_x: f64,
    max_x: f64,
}

impl Chart for NonfiniteThicknessChart {
    fn eval(&self, _point: Point3, _cx: &Cx<'_>) -> ChartSample {
        ChartSample {
            signed_distance: 0.0,
            gradient: Some(Vec3::new(f64::NAN, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: NumericalCertificate::no_claim(),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "test/nonfinite-thickness"
    }
}

impl Chart for UnrepresentableThicknessSlab {
    fn eval(&self, point: Point3, _cx: &Cx<'_>) -> ChartSample {
        let value = if point.x <= self.min_x || point.x >= self.max_x {
            0.0
        } else {
            -1.0
        };
        ChartSample {
            signed_distance: value,
            // The left wall's outward normal is -x, so the inward march is +x.
            gradient: Some(Vec3::new(-1.0, 0.0, 0.0)),
            lipschitz: None,
            error: NumericalCertificate::no_claim(),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(
            Point3::new(self.min_x, -1.0, -1.0),
            Point3::new(self.max_x, 1.0, 1.0),
        )
    }

    fn name(&self) -> &'static str {
        "test/unrepresentable-thickness-slab"
    }
}

/// gq-005a — malformed producer values fail closed, empty aggregates do not
/// return an infinite nominal, and producer cancellation wins after `eval`.
#[test]
fn gq_005a_thickness_authority_and_validation() {
    let malformed =
        with_cx(|cx| thickness_at(&NonfiniteThicknessChart, Point3::new(0.0, 0.0, 0.0), cx));
    let empty = with_cx(|cx| {
        min_thickness(
            &SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            },
            &[],
            cx,
        )
    });
    let translated = 4_503_599_627_370_496.0;
    let no_spatial_progress = with_cx(|cx| {
        thickness_at(
            &UnrepresentableThicknessSlab {
                min_x: translated,
                max_x: translated.next_up(),
            },
            Point3::new(translated, 0.0, 0.0),
            cx,
        )
    });
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let cancelled = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x9E4,
                kernel_id: 5,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        thickness_at(
            &CancellingSeparationSphere {
                gate: &gate,
                during: CancelDuring::Eval,
            },
            Point3::new(1.0, 0.0, 0.0),
            &cx,
        )
    });
    verdict(
        "gq-005a",
        matches!(malformed, Err(QueryError::InvalidThicknessSample { .. }))
            && matches!(empty, Err(QueryError::NoThicknessSamples { skipped: 0 }))
            && matches!(
                no_spatial_progress,
                Err(QueryError::InvalidThicknessArithmetic { .. })
            )
            && matches!(cancelled, Err(QueryError::Cancelled)),
        "thickness refuses nonfinite gradients, empty aggregates, and positive parametric \
         steps that make no representable geometric progress at a translated one-ulp wall; \
         cancellation requested inside thickness eval wins before publication",
    );
}

/// gq-006 — curvature: analytic values on sphere and torus, measured
/// O(h²) convergence for SecondOrder charts, the documented class per
/// chart family, and rotation invariance (G3).
#[test]
fn gq_006_curvature_convergence() {
    with_cx(|cx| {
        let exact = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 2.0,
        };
        let frep = frep_sphere(Point3::new(0.0, 0.0, 0.0), 2.0);
        let p = Point3::new(2.0, 0.0, 0.0);
        // Convergence order: errors at h and h/2 for the frep chart.
        let mut orders = Vec::new();
        for chart in [&exact as &dyn Chart, &frep as &dyn Chart] {
            let e = |h: f64| -> f64 {
                let c = curvature(chart, p, h, cx).expect("curv");
                (c.mean - 0.5).abs()
            };
            let (e1, e2) = (e(0.02), e(0.01));
            orders.push((e2 / e1.max(1e-300)).log2().abs());
        }
        let order_ok = orders.iter().all(|o| (*o - 2.0).abs() < 0.7);
        // Torus principal curvatures at the outer equator: 1/r and
        // 1/(R+r) — signs per outward convention.
        let torus = {
            let mut b = FrepBuilder::new();
            let t = b
                .torus(Point3::new(0.0, 0.0, 0.0), 1.0, 0.3)
                .expect("torus");
            b.finish(t).expect("frep")
        };
        let tp = Point3::new(1.3, 0.0, 0.0);
        let tc = curvature(&torus, tp, 0.01, cx).expect("torus curv");
        let (k1, k2) = (tc.principal[0], tc.principal[1]);
        let torus_ok = (k1 - 1.0 / 1.3).abs().min((k2 - 1.0 / 1.3).abs()) < 1e-2
            && (k1 - 1.0 / 0.3).abs().min((k2 - 1.0 / 0.3).abs()) < 1e-2
            && (tc.gaussian - (1.0 / 0.3) * (1.0 / 1.3)).abs() < 0.05;
        // Classes documented per chart family.
        let tiled = fs_rep_sdf::TiledSdf::build(&exact, 0.08, cx).expect("tiled");
        let mesh = fs_rep_mesh::MeshChart::new(fs_rep_mesh::shapes::icosphere(
            Point3::new(0.0, 0.0, 0.0),
            2.0,
            4,
        ));
        let classes_ok = curvature(&frep, p, 0.01, cx).expect("c").class
            == CurvatureClass::SecondOrder
            && fs_query::curvature_class(&tiled) == CurvatureClass::GridLimited
            && fs_query::curvature_class(&mesh) == CurvatureClass::Estimate;
        // Grid-limited chart still lands near truth at its own scale.
        let ct = curvature(&tiled, p, 0.08, cx).expect("tiled curv");
        let tiled_ok = (ct.mean - 0.5).abs() < 0.08;
        // G3: rotation invariance of curvature scalars (frep rotated).
        let rot = {
            let mut b = FrepBuilder::new();
            let s = b.sphere(Point3::new(0.0, 0.0, 0.0), 2.0).expect("s");
            let r = b.rotate(s, Vec3::new(0.3, -0.5, 0.8), 0.7).expect("rot");
            b.finish(r).expect("frep")
        };
        let cr = curvature(&rot, p, 0.01, cx).expect("rot curv");
        let g3 = (cr.mean - 0.5).abs() < 1e-4 && (cr.gaussian - 0.25).abs() < 1e-3;
        let mut em = fs_obs::Emitter::new("fs-query/conformance", "gq-006/curvature");
        let line = em
            .emit(
                fs_obs::Severity::Info,
                fs_obs::EventKind::Custom {
                    name: "query-curvature-convergence".to_string(),
                    json: format!(
                        "{{\"orders\":[{:.2},{:.2}],\"torus_k\":[{k1:.4},{k2:.4}],\
                         \"tiled_mean\":{:.4}}}",
                        orders[0], orders[1], ct.mean
                    ),
                },
                None,
            )
            .to_jsonl();
        fs_obs::validate_line(&line).expect("curvature event validates");
        println!("{line}");
        verdict(
            "gq-006",
            order_ok && torus_ok && classes_ok && tiled_ok && g3,
            &format!(
                "mean curvature converges at measured order ~2 on SecondOrder charts \
                 ({:.2}, {:.2}), torus principal curvatures hit 1/r and 1/(R+r), \
                 accuracy classes are documented per family (grid-limited lands \
                 within its own scale), and curvature scalars are rotation-invariant",
                orders[0], orders[1]
            ),
        );
    });
}

/// gq-006a — curvature stencils validate every producer result and every
/// arithmetic stage, and producer cancellation wins before publication.
#[test]
fn gq_006a_curvature_fails_closed_on_nonfinite_paths() {
    let malformed = with_cx(|cx| {
        curvature(
            &MalformedPointChart {
                kind: MalformedPointKind::Gradient,
            },
            Point3::new(0.0, 0.0, 0.0),
            0.01,
            cx,
        )
    });
    let overflow = with_cx(|cx| {
        curvature(
            &ExplosiveCurvatureChart,
            Point3::new(0.0, 0.0, 0.0),
            0.01,
            cx,
        )
    });
    let extreme = with_cx(|cx| {
        curvature(
            &ExtendedPlane {
                offset: 1.0e308,
                orientation: 1.0,
                analytic_gradient: false,
                support: Aabb::new(
                    Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
                    Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
                ),
            },
            Point3::new(1.0e308, 0.0, 0.0),
            0.01,
            cx,
        )
    });
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let cancelled = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x9E4,
                kernel_id: 12,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        curvature(
            &CancellingSeparationSphere {
                gate: &gate,
                during: CancelDuring::Eval,
            },
            Point3::new(1.0, 0.0, 0.0),
            0.01,
            &cx,
        )
    });
    verdict(
        "gq-006a",
        matches!(malformed, Err(QueryError::InvalidPointSample { .. }))
            && matches!(overflow, Err(QueryError::InvalidPointArithmetic { .. }))
            && matches!(extreme, Err(QueryError::InvalidPointArithmetic { .. }))
            && matches!(cancelled, Err(QueryError::Cancelled)),
        "curvature refuses malformed samples, overflowing stencil arithmetic, and \
         extreme coordinates where h makes no representable progress; cancellation \
         requested inside eval wins before finite scalars are published",
    );
}

#[derive(Clone, Copy)]
struct ExtendedPlane {
    offset: f64,
    orientation: f64,
    analytic_gradient: bool,
    support: Aabb,
}

impl Chart for ExtendedPlane {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let value = self.orientation * (x.x - self.offset);
        ChartSample {
            signed_distance: value,
            gradient: self
                .analytic_gradient
                .then_some(Vec3::new(self.orientation, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: NumericalCertificate::exact(value),
        }
    }

    fn support(&self) -> Aabb {
        self.support
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn name(&self) -> &'static str {
        "test/extended-plane"
    }
}

struct InfiniteCylinder {
    radius: f64,
}

impl Chart for InfiniteCylinder {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let rho = x.x.hypot(x.y);
        let value = rho - self.radius;
        ChartSample {
            signed_distance: value,
            gradient: (rho > 0.0).then_some(Vec3::new(x.x / rho, x.y / rho, 0.0)),
            lipschitz: Some(1.0),
            error: NumericalCertificate::exact(value),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(
            Point3::new(-self.radius, -self.radius, f64::NEG_INFINITY),
            Point3::new(self.radius, self.radius, f64::INFINITY),
        )
    }

    fn name(&self) -> &'static str {
        "test/infinite-cylinder"
    }
}

struct PanicEvalChart;

impl Chart for PanicEvalChart {
    fn eval(&self, _x: Point3, _cx: &Cx<'_>) -> ChartSample {
        panic!("separation must reject an overflowing work count before chart evaluation")
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn name(&self) -> &'static str {
        "test/panic-eval"
    }
}

struct ForgedLocalSeparationChart;

impl Chart for ForgedLocalSeparationChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        ChartSample {
            signed_distance: x.x,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            // These locally plausible fields deliberately do not upgrade the
            // trait's default global NoClaim theorem.
            lipschitz: Some(1.0),
            error: NumericalCertificate::enclosure(x.x.next_down(), x.x.next_up()),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn trace_value_enclosure(
        &self,
        _x: Point3,
        sample: &ChartSample,
        _cx: &Cx<'_>,
    ) -> NumericalCertificate {
        NumericalCertificate::enclosure(
            sample.signed_distance.next_down(),
            sample.signed_distance.next_up(),
        )
    }

    fn name(&self) -> &'static str {
        "test/forged-local-separation"
    }
}

struct MalformedSeparationChart {
    certificate: NumericalCertificate,
}

impl Chart for MalformedSeparationChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        ChartSample {
            signed_distance: x.x,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: NumericalCertificate::exact(x.x),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn trace_value_enclosure(
        &self,
        _x: Point3,
        _sample: &ChartSample,
        _cx: &Cx<'_>,
    ) -> NumericalCertificate {
        self.certificate
    }

    fn name(&self) -> &'static str {
        "test/malformed-separation"
    }
}

/// gq-007 — every support-derived sampler resolves a finite domain before
/// span/count arithmetic. Analytic point queries remain valid on honest
/// unbounded charts, finite-difference fallbacks use caller scale, clipped
/// separation carries local authority, and aggregate thickness cannot hide a
/// domain refusal as a skipped local sample.
#[test]
#[allow(clippy::too_many_lines)] // One admission campaign compares all support-derived samplers.
fn gq_007_unbounded_sampling_admission() {
    with_cx(|cx| {
        let unbounded = Aabb::new(
            Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
            Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
        );
        let analytic = ExtendedPlane {
            offset: 0.0,
            orientation: 1.0,
            analytic_gradient: true,
            support: unbounded,
        };
        let numerical = ExtendedPlane {
            analytic_gradient: false,
            ..analytic
        };
        let clip = Aabb::new(Point3::new(-2.0, -2.0, -2.0), Point3::new(3.0, 2.0, 2.0));

        let analytic_cp = closest_point(&analytic, Point3::new(2.0, 0.25, -0.5), cx)
            .expect("analytic closest point needs no finite sampling domain");
        let fd_default = closest_point(&numerical, Point3::new(2.0, 0.25, -0.5), cx);
        let fd_clipped = closest_point_clipped(&numerical, Point3::new(2.0, 0.25, -0.5), clip, cx)
            .expect("explicit clip supplies the FD scale");
        let closest_ok = analytic_cp.residual < 1e-12
            && matches!(
                fd_default,
                Err(QueryError::SamplingDomain(
                    SamplingDomainError::UnboundedSupport { .. }
                ))
            )
            && fd_clipped.residual < 1e-10;

        let flat = curvature(&numerical, Point3::new(0.0, 0.2, -0.4), 0.01, cx)
            .expect("curvature uses its explicit h for the missing gradient");
        let invalid_steps = [0.0, -0.01, f64::MIN_POSITIVE, f64::INFINITY, f64::NAN]
            .into_iter()
            .all(|h| {
                matches!(
                    curvature(&analytic, Point3::new(0.0, 0.0, 0.0), h, cx),
                    Err(QueryError::InvalidFiniteDifferenceStep { .. })
                )
            });
        let curvature_ok = flat.mean.abs() < 1e-12
            && flat.gaussian.abs() < 1e-12
            && flat.principal.iter().all(|value| value.abs() < 1e-12)
            && invalid_steps;

        let opposite = ExtendedPlane {
            offset: 1.0,
            orientation: -1.0,
            analytic_gradient: true,
            support: unbounded,
        };
        let separation_default = separation(&analytic, &opposite, 12, cx);
        let separation_local = separation_clipped(&analytic, &opposite, 12, clip, cx)
            .expect("clip resolves both extended supports");
        let malformed = ExtendedPlane {
            support: Aabb::new(
                Point3::new(f64::NAN, f64::NEG_INFINITY, f64::NEG_INFINITY),
                Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
            ),
            ..analytic
        };
        let malformed_pair = separation_clipped(&malformed, &opposite, 4, clip, cx);
        let too_large = separation(&PanicEvalChart, &PanicEvalChart, u32::MAX, cx);
        let over_cap = separation(&PanicEvalChart, &PanicEvalChart, 100, cx);
        let huge_clip = Aabb::new(
            Point3::new(-2.0e307, -2.0e307, -2.0e307),
            Point3::new(2.0e307, 2.0e307, 2.0e307),
        );
        let huge_local = separation_clipped(&analytic, &opposite, 2, huge_clip, cx)
            .expect("ratio-first coordinates keep an extreme finite clip representable");
        let separation_ok = matches!(
            separation_default,
            Err(QueryError::SamplingDomain(
                SamplingDomainError::UnboundedSupport { .. }
            ))
        ) && separation_local.scope == SeparationScope::ClippedLocal
            && separation_local.domain == clip
            && clip.contains(separation_local.witness)
            && separation_local.observed.is_finite()
            && separation_local.lower_bound.is_finite()
            && matches!(
                malformed_pair,
                Err(QueryError::SamplingDomain(
                    SamplingDomainError::InvalidSupport { .. }
                ))
            )
            && matches!(too_large, Err(QueryError::SamplingGridTooLarge { .. }))
            && matches!(
                over_cap,
                Err(QueryError::SamplingWorkLimitExceeded {
                    limit: SEPARATION_MAX_CHART_SAMPLES,
                    ..
                })
            )
            && huge_local.observed.is_finite()
            && huge_local.lower_bound.is_finite()
            && point_is_finite_for_test(huge_local.witness);

        let cylinder = InfiniteCylinder { radius: 1.0 };
        let boundary = Point3::new(1.0, 0.0, 0.0);
        let thickness_default = thickness_at(&cylinder, boundary, cx);
        let thickness_local = thickness_at_clipped(&cylinder, boundary, clip, cx)
            .expect("clip bounds the axial march");
        let anisotropic_clip = Aabb::new(
            Point3::new(-2.0, -1.0e6, -1.0e6),
            Point3::new(2.0, 1.0e6, 1.0e6),
        );
        let anisotropic_local = thickness_at_clipped(&cylinder, boundary, anisotropic_clip, cx)
            .expect("transverse clip scale cannot skip the opposite wall");
        let minimum_default = min_thickness(&cylinder, &[boundary], cx);
        let minimum_local = min_thickness_clipped(&cylinder, &[boundary], clip, cx)
            .expect("clip bounds the aggregate thickness query");
        let thickness_ok = matches!(
            thickness_default,
            Err(QueryError::SamplingDomain(
                SamplingDomainError::UnboundedSupport { .. }
            ))
        ) && (thickness_local.value - 2.0).abs() < 1e-9
            && (anisotropic_local.value - 2.0).abs() < 1e-9
            && matches!(
                minimum_default,
                Err(QueryError::SamplingDomain(
                    SamplingDomainError::UnboundedSupport { .. }
                ))
            )
            && thickness_local.authority == NumericalKind::Estimate
            && anisotropic_local.authority == NumericalKind::Estimate
            && (minimum_local.value - 2.0).abs() < 1e-9
            && minimum_local.skipped == 0
            && minimum_local.authority == NumericalKind::Estimate;

        verdict(
            "gq-007",
            closest_ok && curvature_ok && separation_ok && thickness_ok,
            "analytic point queries bypass finite-domain admission; FD fallbacks use \
             explicit caller scale; malformed, unbounded, and overflowing separation \
             domains refuse before evaluation; clipped separation is marked local; and \
             thickness domain failures propagate instead of becoming skipped samples",
        );
    });
}

fn point_is_finite_for_test(point: Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
}

/// Teaching-refusal spot checks.
#[test]
fn refusals_teach() {
    with_cx(|cx| {
        let exact = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let off = thickness_at(&exact, Point3::new(1.5, 0.0, 0.0), cx).expect_err("off boundary");
        assert!(matches!(off, QueryError::NotOnBoundary { .. }), "{off}");
        assert!(off.to_string().contains("project"), "{off}");
    });
}
