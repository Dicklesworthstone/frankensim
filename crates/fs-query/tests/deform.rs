//! Deformation-hook battery (bead rjnd, part 6).
//!
//! - gh-001 G0/G5: a scaled+translated sphere — the composed field
//!   keeps sign and zero set, the claimed `|f|/L` safe step never
//!   exceeds the true current-configuration distance, and replay is
//!   bit-identical.
//! - gh-002 G0: the trace enclosure passes through rigorously and
//!   collapses to no-claim on Estimate-class reference evidence.
//! - gh-003 G0/G4: weaker reference claims, malformed Lipschitz
//!   bounds, non-finite queries, and broken pull-backs all refuse —
//!   through typed errors or no-claim samples.
//! - gh-004 G0: the deformed chart cannot launder a distance claim
//!   into exact-distance consumers.
//!
//! Aggregate outcomes use canonical fs-obs events; evaluated cases carry
//! the shared execution seed and the constructor-only case uses zero.

use asupersync::types::Budget;
use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Aabb, Chart, ChartSample, Point3, TraceStepClaim, Vec3};
use fs_query::{DeformationMap, DeformedChart, ImplicitGapOracle, OffsetChart, QueryError};

const EXECUTION_SEED: u64 = 0xDEF0;

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new("fs-query/deform", case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: "fs-query/deform".to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("deformation verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("deformation verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: EXECUTION_SEED,
                kernel_id: 19,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// Uniform scale + translation: forward `F(p) = s·p + d`.
struct ScaleShift {
    scale: f64,
    shift: Vec3,
}

impl DeformationMap for ScaleShift {
    fn pull_back(&self, x: Point3) -> Point3 {
        Point3::new(
            (x.x - self.shift.x) / self.scale,
            (x.y - self.shift.y) / self.scale,
            (x.z - self.shift.z) / self.scale,
        )
    }

    fn pull_back_lipschitz(&self) -> f64 {
        (1.0 / self.scale).next_up()
    }

    fn name(&self) -> &'static str {
        "test/scale-shift"
    }
}

/// A pull-back that produces non-finite points.
struct BrokenMap;

impl DeformationMap for BrokenMap {
    fn pull_back(&self, _x: Point3) -> Point3 {
        Point3::new(f64::NAN, 0.0, 0.0)
    }

    fn pull_back_lipschitz(&self) -> f64 {
        1.0
    }

    fn name(&self) -> &'static str {
        "test/broken-map"
    }
}

/// A map whose declared Lipschitz bound is unusable.
struct LyingMap(f64);

impl DeformationMap for LyingMap {
    fn pull_back(&self, x: Point3) -> Point3 {
        x
    }

    fn pull_back_lipschitz(&self) -> f64 {
        self.0
    }

    fn name(&self) -> &'static str {
        "test/lying-map"
    }
}

/// Claims ExactDistance but serves Estimate-class evidence.
struct EstimatingSphere(SphereChart);

impl Chart for EstimatingSphere {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        let mut s = self.0.eval(x, cx);
        s.error = NumericalCertificate::estimate(s.signed_distance, s.signed_distance);
        s
    }

    fn support(&self) -> Aabb {
        self.0.support()
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn name(&self) -> &'static str {
        "test/estimating-sphere"
    }
}

fn reference_sphere() -> SphereChart {
    SphereChart {
        center: Point3::new(0.25, 0.0, 0.0),
        radius: 0.5,
    }
}

fn wide_box() -> Aabb {
    Aabb::new(Point3::new(-8.0, -8.0, -8.0), Point3::new(8.0, 8.0, 8.0))
}

#[test]
fn gh_001_sign_zero_set_and_safe_steps_survive_the_deformation() {
    let reference = reference_sphere();
    let map = ScaleShift {
        scale: 2.0,
        shift: Vec3::new(1.0, -0.5, 0.25),
    };
    let deformed = DeformedChart::new(&reference, &map, wide_box()).expect("hook");
    // Current-configuration surface: sphere at F(c) with radius s·r.
    let current_center = [1.0 + 2.0 * 0.25, -0.5, 0.25];
    let current_radius = 1.0;
    let probes = [
        [0.0, 0.0, 0.0],
        [1.5, -0.5, 0.25],
        [2.9, -0.5, 0.25],
        [1.5, 0.9, 0.25],
        [1.5, -0.5, -1.4],
        [-2.0, 3.0, 1.0],
        [1.6, -0.4, 0.3],
        [3.4, 1.2, -0.7],
    ];
    let mut checked = 0usize;
    with_cx(|cx| {
        for q in probes {
            let x = Point3::new(q[0], q[1], q[2]);
            let first = deformed.eval(x, cx);
            let second = deformed.eval(x, cx);
            assert_eq!(
                first.signed_distance.to_bits(),
                second.signed_distance.to_bits(),
                "replay is bit-identical"
            );
            let dx = [
                q[0] - current_center[0],
                q[1] - current_center[1],
                q[2] - current_center[2],
            ];
            let true_dist = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt() - current_radius;
            if true_dist.abs() > 1e-9 {
                assert_eq!(
                    first.signed_distance > 0.0,
                    true_dist > 0.0,
                    "sign must transfer at {q:?} (field {}, true {true_dist})",
                    first.signed_distance
                );
            }
            let l = first.lipschitz.expect("composed Lipschitz claim");
            let safe_step = first.signed_distance.abs() / l;
            assert!(
                safe_step <= true_dist.abs() + 1e-9,
                "safe step {safe_step} exceeds the true distance {true_dist} at {q:?}"
            );
            checked += 1;
        }
    });
    verdict(
        "gh-001",
        checked == probes.len(),
        &format!("{checked} probes: sign transfer + certified safe steps"),
        EXECUTION_SEED,
    );
}

#[test]
fn gh_002_enclosures_pass_through_or_collapse() {
    let reference = reference_sphere();
    let map = ScaleShift {
        scale: 2.0,
        shift: Vec3::new(1.0, -0.5, 0.25),
    };
    let deformed = DeformedChart::new(&reference, &map, wide_box()).expect("hook");
    let lying_reference = EstimatingSphere(reference_sphere());
    let lying =
        DeformedChart::new(&lying_reference, &map, wide_box()).expect("claims pass construction");
    let x = Point3::new(0.2, 0.1, 0.0);
    let (honest_kind, honest_contains, lying_kind) = with_cx(|cx| {
        let sample = deformed.eval(x, cx);
        let enclosure = deformed.trace_value_enclosure(x, &sample, cx);
        let bad_sample = lying.eval(x, cx);
        let bad = lying.trace_value_enclosure(x, &bad_sample, cx);
        (
            enclosure.kind,
            enclosure.lo <= sample.signed_distance && sample.signed_distance <= enclosure.hi,
            bad.kind,
        )
    });
    assert!(
        matches!(honest_kind, NumericalKind::Exact | NumericalKind::Enclosure),
        "the reference enclosure passes through rigorously"
    );
    assert!(honest_contains, "the enclosure contains the field value");
    assert_eq!(
        lying_kind,
        NumericalKind::NoClaim,
        "Estimate-class reference evidence collapses to no-claim"
    );
    verdict(
        "gh-002",
        true,
        "pass-through rigorous, collapse honest",
        EXECUTION_SEED,
    );
}

#[test]
fn gh_003_refusals_fail_closed() {
    let reference = reference_sphere();
    let map = ScaleShift {
        scale: 2.0,
        shift: Vec3::new(0.0, 0.0, 0.0),
    };
    // Weaker reference claim refuses at construction.
    let offset = OffsetChart::new(&reference, 0.1).expect("finite radius");
    assert!(matches!(
        DeformedChart::new(&offset, &map, wide_box()),
        Err(QueryError::DeformationRequiresExactDistance { .. })
    ));
    // Unusable Lipschitz bounds refuse at construction.
    for l in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            matches!(
                DeformedChart::new(&reference, &LyingMap(l), wide_box()),
                Err(QueryError::DeformationInvalidMap { .. })
            ),
            "pull-back Lipschitz {l} must refuse"
        );
    }
    // Non-finite queries and broken pull-backs yield no-claim samples.
    let deformed = DeformedChart::new(&reference, &map, wide_box()).expect("hook");
    let broken = DeformedChart::new(&reference, &BrokenMap, wide_box()).expect("hook");
    let (nan_query, broken_sample) = with_cx(|cx| {
        (
            deformed.eval(Point3::new(f64::NAN, 0.0, 0.0), cx),
            broken.eval(Point3::new(0.0, 0.0, 0.0), cx),
        )
    });
    for (label, sample) in [("nan-query", nan_query), ("broken-map", broken_sample)] {
        assert_eq!(
            sample.error.kind,
            NumericalKind::NoClaim,
            "{label} must produce a no-claim sample"
        );
        assert!(sample.lipschitz.is_none(), "{label} offers no step fuel");
        assert!(
            sample.signed_distance.is_nan(),
            "{label} publishes no plausible value"
        );
    }
    verdict(
        "gh-003",
        true,
        "construction and sample refusals all fail closed",
        EXECUTION_SEED,
    );
}

#[test]
fn gh_004_no_distance_claim_launders_through() {
    let reference = reference_sphere();
    let map = ScaleShift {
        scale: 2.0,
        shift: Vec3::new(1.0, -0.5, 0.25),
    };
    let deformed = DeformedChart::new(&reference, &map, wide_box()).expect("hook");
    assert_eq!(
        deformed.trace_step_claim(),
        TraceStepClaim::LipschitzImplicit
    );
    // Exact-distance consumers must refuse the deformed chart.
    let other = reference_sphere();
    let refused = ImplicitGapOracle::new(&deformed, &other);
    assert!(matches!(
        refused,
        Err(QueryError::SeparationRequiresExactDistance { input: "a", .. })
    ));
    verdict(
        "gh-004",
        true,
        "LipschitzImplicit claim holds; exact-distance consumers refuse the wrapper",
        0,
    );
}
