//! Implicit gap-oracle battery (bead rjnd, part 4 + the part-7
//! pointwise overlap-inradius witness).
//!
//! - gg-001 G0/G5: disjoint spheres — the sum enclosure contains the
//!   analytic value, the separation upper bound really bounds the true
//!   set distance, the contact axis points along the gap, and replay
//!   is bit-identical.
//! - gg-002 G0: overlapping spheres — a certified common-ball radius
//!   appears, never exceeding the analytic depth, and no separation
//!   bound is offered.
//! - gg-003 G0/G4: capability refusal at construction, malformed
//!   evidence refusal per sample, non-finite points, cancellation.
//! - gg-004 G0: honestly absent gradients yield no normal claim.

use asupersync::types::Budget;
use fs_evidence::NumericalCertificate;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Aabb, Chart, ChartSample, Point3, TraceStepClaim};
use fs_query::{ImplicitGapOracle, OffsetChart, QueryError};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-query/gap\",\"case\":\"{case}\",\"verdict\":\"{}\",\
         \"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
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
                seed: 0x6A9,
                kernel_id: 16,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn sphere(x: f64, r: f64) -> SphereChart {
    SphereChart {
        center: Point3::new(x, 0.0, 0.0),
        radius: r,
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

#[test]
fn gg_001_disjoint_spheres_bound_the_true_distance() {
    let a = sphere(-1.0, 0.5);
    let b = sphere(1.0, 0.5);
    let oracle = ImplicitGapOracle::new(&a, &b).expect("exact-distance pair");
    let mid = Point3::new(0.0, 0.0, 0.0);
    let (first, second) = with_cx(|cx| {
        (
            oracle.gap_at(mid, cx).expect("gap sample"),
            oracle.gap_at(mid, cx).expect("replayed sample"),
        )
    });
    // φ_A(0) + φ_B(0) = 0.5 + 0.5 = 1.0; the true set distance is 1.0.
    assert!(
        first.sum_lo <= 1.0 && 1.0 <= first.sum_hi,
        "sum enclosure [{}, {}] must contain 1.0",
        first.sum_lo,
        first.sum_hi
    );
    let upper = first.separation_upper.expect("outside both bodies");
    assert!(
        upper >= 1.0 - 1e-9,
        "the upper bound {upper} may not undercut the true distance 1.0"
    );
    assert!(first.overlap_inradius.is_none());
    let normal = first.normal.expect("both gradients exist at the midpoint");
    assert!(
        normal[0].abs() > 0.999 && normal[1].abs() < 1e-9 && normal[2].abs() < 1e-9,
        "the contact axis lies along x, got {normal:?}"
    );
    assert_eq!(first.sum_lo.to_bits(), second.sum_lo.to_bits());
    assert_eq!(first.sum_hi.to_bits(), second.sum_hi.to_bits());
    verdict(
        "gg-001",
        true,
        &format!(
            "sum [{:.9}, {:.9}] contains 1.0; upper {:.9}; axis {normal:?}",
            first.sum_lo, first.sum_hi, upper
        ),
    );
}

#[test]
fn gg_002_overlapping_spheres_yield_a_certified_witness() {
    let a = sphere(-0.25, 1.0);
    let b = sphere(0.25, 1.0);
    let oracle = ImplicitGapOracle::new(&a, &b).expect("exact-distance pair");
    let mid = Point3::new(0.0, 0.0, 0.0);
    let sample = with_cx(|cx| oracle.gap_at(mid, cx)).expect("overlap sample");
    // φ_A(0) = φ_B(0) = 0.25 - 1.0 = -0.75.
    let witness = sample.overlap_inradius.expect("certified common ball");
    assert!(
        witness <= 0.75 && witness > 0.75 - 1e-9,
        "witness radius {witness} must certify but never exceed the analytic 0.75"
    );
    assert!(
        sample.separation_upper.is_none(),
        "an interior point offers no separation bound"
    );
    verdict(
        "gg-002",
        true,
        &format!("common-ball witness {witness:.9} ≤ 0.75"),
    );
}

#[test]
fn gg_003_refusals_fail_closed() {
    let a = sphere(-1.0, 0.5);
    let b = sphere(1.0, 0.5);

    // Weaker trace claim refuses at construction, naming the input.
    let offset = OffsetChart::new(&a, 0.1).expect("finite radius");
    let weaker = ImplicitGapOracle::new(&offset, &b);
    assert!(matches!(
        weaker,
        Err(QueryError::SeparationRequiresExactDistance { input: "a", .. })
    ));
    let weaker_b = ImplicitGapOracle::new(&a, &offset);
    assert!(matches!(
        weaker_b,
        Err(QueryError::SeparationRequiresExactDistance { input: "b", .. })
    ));

    // Estimate-class evidence refuses per sample despite the claim.
    let lying = EstimatingSphere(sphere(1.0, 0.5));
    let oracle = ImplicitGapOracle::new(&a, &lying).expect("claims pass construction");
    let refused = with_cx(|cx| oracle.gap_at(Point3::new(0.0, 0.0, 0.0), cx));
    assert!(matches!(
        refused,
        Err(QueryError::InvalidTraceSample { .. })
    ));

    // Non-finite query point refuses before any chart call.
    let honest = ImplicitGapOracle::new(&a, &b).expect("pair");
    let nan_point = with_cx(|cx| honest.gap_at(Point3::new(f64::NAN, 0.0, 0.0), cx));
    assert!(matches!(
        nan_point,
        Err(QueryError::InvalidPointSample { .. })
    ));

    // Cancellation surfaces as the typed refusal.
    let gate = CancelGate::new();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let cancelled = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x6A9,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        honest.gap_at(Point3::new(0.0, 0.0, 0.0), &cx)
    });
    assert!(matches!(cancelled, Err(QueryError::Cancelled)));
    verdict(
        "gg-003",
        true,
        "claim/evidence/point/cancellation refusals all typed",
    );
}

#[test]
fn gg_004_absent_gradients_make_no_normal_claim() {
    let a = sphere(-1.0, 0.5);
    let b = sphere(1.0, 0.5);
    let oracle = ImplicitGapOracle::new(&a, &b).expect("pair");
    // At A's center the sphere fixture honestly declines a gradient.
    let sample = with_cx(|cx| oracle.gap_at(Point3::new(-1.0, 0.0, 0.0), cx)).expect("sample");
    assert!(
        sample.normal.is_none(),
        "a declined gradient must yield no contact-axis claim"
    );
    verdict("gg-004", true, "no gradient, no normal claim");
}
