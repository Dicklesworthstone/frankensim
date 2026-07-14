//! Certified convex-separation battery (bead rjnd, part 2).
//!
//! - gc-001 G0: sphere-sphere and sphere-box enclosures contain the
//!   analytic distances, with proven separation and tight widths
//!   (smooth pairs converge fast).
//! - gc-002 G0: box-box face-gap enclosure contains the analytic gap
//!   and proves separation (nonsmooth pair: width honest, not tight).
//! - gc-003 G0: touching and overlapping pairs keep 0 inside the
//!   enclosure and never claim separation.
//! - gc-004 G5: identical inputs replay bit-identical enclosures.
//! - gc-005 G0/G4: constructor refusals, zero budget, and cancellation
//!   all fail closed with the named typed error.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Aabb, Point3};
use fs_query::{ConvexBox, ConvexSphere, QueryError, convex_separation};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-query/convex\",\"case\":\"{case}\",\"verdict\":\"{}\",\
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
                seed: 0xC0F,
                kernel_id: 12,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn sphere(x: f64, y: f64, z: f64, r: f64) -> ConvexSphere {
    ConvexSphere::new(Point3::new(x, y, z), r).expect("valid sphere")
}

fn boxx(min: [f64; 3], max: [f64; 3]) -> ConvexBox {
    ConvexBox::new(Aabb::new(
        Point3::new(min[0], min[1], min[2]),
        Point3::new(max[0], max[1], max[2]),
    ))
    .expect("valid box")
}

#[test]
fn gc_001_smooth_pairs_contain_analytic_distances() {
    // Spheres on a skew axis: distance = |c1-c2| - r1 - r2.
    let a = sphere(-0.5, -0.25, 0.125, 0.25);
    let b = sphere(0.7, 0.35, -0.375, 0.375);
    let truth = {
        let d: [f64; 3] = [1.2, 0.6, -0.5];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() - 0.25 - 0.375
    };
    let sep = with_cx(|cx| convex_separation(&a, &b, 256, cx)).expect("sphere pair");
    assert!(
        sep.lo <= truth && truth <= sep.hi,
        "sphere-sphere: [{}, {}] must contain {truth}",
        sep.lo,
        sep.hi
    );
    assert!(sep.separation_proven, "clear gap must be proven");
    assert!(
        sep.hi - sep.lo < 1e-9,
        "smooth pair converges tight, got width {}",
        sep.hi - sep.lo
    );

    // Sphere vs box along +x: gap = (box min x) - (sphere max x).
    let s = sphere(-0.5, 0.0, 0.0, 0.25);
    let bx = boxx([0.5, -1.0, -1.0], [1.5, 1.0, 1.0]);
    let truth_sb = 0.5 - (-0.5 + 0.25);
    let sep_sb = with_cx(|cx| convex_separation(&s, &bx, 256, cx)).expect("sphere-box");
    assert!(
        sep_sb.lo <= truth_sb && truth_sb <= sep_sb.hi,
        "sphere-box: [{}, {}] must contain {truth_sb}",
        sep_sb.lo,
        sep_sb.hi
    );
    assert!(sep_sb.separation_proven);
    verdict(
        "gc-001",
        true,
        &format!(
            "sphere-sphere width {:.3e}, sphere-box [{:.6}, {:.6}] vs {truth_sb}",
            sep.hi - sep.lo,
            sep_sb.lo,
            sep_sb.hi
        ),
    );
}

#[test]
fn gc_002_box_box_face_gap() {
    let a = boxx([-1.0, -0.5, -0.5], [-0.25, 0.5, 0.5]);
    let b = boxx([0.35, -0.4, -0.4], [1.0, 0.6, 0.6]);
    let truth = 0.35 - (-0.25);
    let sep = with_cx(|cx| convex_separation(&a, &b, 4096, cx)).expect("box pair");
    assert!(
        sep.lo <= truth && truth <= sep.hi,
        "box-box: [{}, {}] must contain {truth}",
        sep.lo,
        sep.hi
    );
    assert!(sep.separation_proven, "face gap must be proven");
    assert!(
        sep.hi - sep.lo < 0.05,
        "nonsmooth width stays honest but bounded, got {}",
        sep.hi - sep.lo
    );
    verdict(
        "gc-002",
        true,
        &format!(
            "box-box [{:.6}, {:.6}] contains {truth}, width {:.3e}",
            sep.lo,
            sep.hi,
            sep.hi - sep.lo
        ),
    );
}

#[test]
fn gc_003_touching_and_overlap_never_claim_separation() {
    // Exactly touching spheres: distance is 0.
    let a = sphere(-0.5, 0.0, 0.0, 0.5);
    let b = sphere(0.5, 0.0, 0.0, 0.5);
    let touch = with_cx(|cx| convex_separation(&a, &b, 256, cx)).expect("touching");
    assert!(
        touch.lo <= 0.0 + 1e-12 && 0.0 <= touch.hi,
        "touching: [{}, {}] must contain 0",
        touch.lo,
        touch.hi
    );
    assert!(!touch.separation_proven, "touching must not prove a gap");

    // Overlapping boxes: distance is 0.
    let c = boxx([-1.0, -1.0, -1.0], [0.25, 1.0, 1.0]);
    let d = boxx([-0.25, -0.9, -0.9], [1.0, 0.9, 0.9]);
    let overlap = with_cx(|cx| convex_separation(&c, &d, 1024, cx)).expect("overlap");
    assert!(
        overlap.lo <= 0.0 && overlap.hi >= 0.0,
        "overlap: [{}, {}] must contain 0 with lo clamped at 0",
        overlap.lo,
        overlap.hi
    );
    assert!(!overlap.separation_proven);
    verdict(
        "gc-003",
        true,
        &format!(
            "touching hi {:.3e}; overlap hi {:.3e}; neither proves separation",
            touch.hi, overlap.hi
        ),
    );
}

#[test]
fn gc_004_replay_is_bit_identical() {
    let a = sphere(-0.5, -0.25, 0.125, 0.25);
    let b = boxx([0.35, -0.4, -0.4], [1.0, 0.6, 0.6]);
    let first = with_cx(|cx| convex_separation(&a, &b, 512, cx)).expect("first run");
    let second = with_cx(|cx| convex_separation(&a, &b, 512, cx)).expect("second run");
    assert_eq!(first.lo.to_bits(), second.lo.to_bits());
    assert_eq!(first.hi.to_bits(), second.hi.to_bits());
    assert_eq!(first.iterations, second.iterations);
    assert_eq!(
        first.witness_a.map(f64::to_bits),
        second.witness_a.map(f64::to_bits)
    );
    assert_eq!(
        first.witness_b.map(f64::to_bits),
        second.witness_b.map(f64::to_bits)
    );
    verdict("gc-004", true, "identical inputs replay bit-identically");
}

#[test]
fn gc_005_refusals_fail_closed() {
    let bad_sphere = ConvexSphere::new(Point3::new(f64::NAN, 0.0, 0.0), 0.5);
    let flat_sphere = ConvexSphere::new(Point3::new(0.0, 0.0, 0.0), 0.0);
    let bad_box = ConvexBox::new(Aabb {
        min: Point3::new(0.0, 0.0, 0.0),
        max: Point3::new(0.0, 1.0, 1.0),
    });
    assert!(matches!(
        bad_sphere,
        Err(QueryError::ConvexInvalidShape { .. })
    ));
    assert!(matches!(
        flat_sphere,
        Err(QueryError::ConvexInvalidShape { .. })
    ));
    assert!(matches!(
        bad_box,
        Err(QueryError::ConvexInvalidShape { .. })
    ));

    let a = sphere(-0.5, 0.0, 0.0, 0.25);
    let b = sphere(0.5, 0.0, 0.0, 0.25);
    let zero_budget = with_cx(|cx| convex_separation(&a, &b, 0, cx));
    assert!(matches!(
        zero_budget,
        Err(QueryError::ConvexInvalidShape { .. })
    ));

    let gate = CancelGate::new();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let cancelled = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xC0F,
                kernel_id: 13,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        convex_separation(&a, &b, 256, &cx)
    });
    assert!(matches!(cancelled, Err(QueryError::Cancelled)));
    verdict(
        "gc-005",
        true,
        "invalid shapes, zero budget, and cancellation refuse typed",
    );
}
