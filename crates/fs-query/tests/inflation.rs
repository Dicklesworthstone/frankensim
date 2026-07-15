//! Contact-bound inflation conformance (bead fugfk).
//!
//! The carrier admits only rigorous conversion/route receipts. Once admitted,
//! increasing its radius may only widen separation/gap/codimension brackets,
//! shrink overlap witnesses, add CCD candidates, and weaken moment-cell
//! classification. Exact zero is a bit-neutral compatibility path.

use std::collections::BTreeMap;

use asupersync::types::Budget;
use fs_evidence::{Certified, Evidence, ProvenanceHash};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::{BoxChart, SphereChart};
use fs_geom::{
    Aabb, Chart, ConverterSpec, EdgeOutcome, EdgeRunner, ErrorModel, MemoryCostOracle, Point3,
    RouteRequest, Router,
};
use fs_query::{
    CodimGap, CodimThickness, ContactInflation, ConvexSeparation, ConvexSphere, GapSample,
    GeometricMoments, ImplicitGapOracle, MomentEnclosure, ccd_candidates,
    ccd_candidates_with_inflation, codim_gap, codim_gap_with_inflation, convex_separation,
    convex_separation_with_inflation, geometric_moments, geometric_moments_with_inflation,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x1F1A_7100,
                kernel_id: 23,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn conversion_receipt(qoi: f64, hi: f64) -> Certified<f64> {
    Evidence::enclosed(
        qoi,
        0.0,
        hi,
        ProvenanceHash::of_bytes(b"fs-query/inflation/conversion"),
    )
    .certified()
    .expect("finite rigorous conversion receipt")
}

struct FixedCertifiedRunner {
    error: f64,
}

impl EdgeRunner for FixedCertifiedRunner {
    fn run(&self, _cx: &Cx<'_>) -> Result<EdgeOutcome, String> {
        EdgeOutcome::certified(conversion_receipt(self.error, self.error), 0.001)
            .map_err(|error| error.to_string())
    }
}

fn routed_inflation(radius: f64, cx: &Cx<'_>) -> ContactInflation {
    let mut router = Router::new();
    router
        .register(ConverterSpec {
            from: "native".to_string(),
            to: "converted".to_string(),
            name: "native->converted/inflation-test".to_string(),
            base_cost_s: 0.001,
            error: ErrorModel::AdditiveAbs(radius),
            certified: true,
        })
        .expect("valid converter spec");
    let mut oracle = MemoryCostOracle::new();
    let plan = router
        .plan(
            &RouteRequest {
                from: "native".to_string(),
                to: "converted".to_string(),
                scale: 1.0,
                max_abs_error: radius,
                max_cost_s: 1.0,
            },
            &oracle,
        )
        .expect("one-edge route");
    let runners: BTreeMap<String, Box<dyn EdgeRunner>> = BTreeMap::from([(
        "native->converted/inflation-test".to_string(),
        Box::new(FixedCertifiedRunner { error: radius }) as Box<dyn EdgeRunner>,
    )]);
    let outcome = router
        .execute(&plan, &runners, &mut oracle, cx)
        .expect("certified route executes");
    ContactInflation::from_route(&outcome)
        .expect("sealed rigorous route receipt admits its upper endpoint")
}

fn certified_inflation(radius: f64) -> ContactInflation {
    ContactInflation::from_conversion(&conversion_receipt(radius, radius))
        .expect("certified absolute-error receipt")
}

fn sphere_chart(x: f64, radius: f64) -> SphereChart {
    SphereChart {
        center: Point3::new(x, 0.0, 0.0),
        radius,
    }
}

fn convex_sphere(x: f64, radius: f64) -> ConvexSphere {
    ConvexSphere::new(Point3::new(x, 0.0, 0.0), radius).expect("valid convex sphere")
}

fn tetra(origin: [f64; 3], scale: f64) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let [x, y, z] = origin;
    (
        vec![
            [x, y, z],
            [x + scale, y, z],
            [x, y + scale, z],
            [x, y, z + scale],
        ],
        vec![[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]],
    )
}

fn assert_convex_bits_eq(actual: ConvexSeparation, expected: ConvexSeparation) {
    assert_eq!(actual.lo.to_bits(), expected.lo.to_bits());
    assert_eq!(actual.hi.to_bits(), expected.hi.to_bits());
    assert_eq!(actual.separation_proven, expected.separation_proven);
    assert_eq!(actual.iterations, expected.iterations);
    assert_eq!(
        actual.witness_a.map(f64::to_bits),
        expected.witness_a.map(f64::to_bits)
    );
    assert_eq!(
        actual.witness_b.map(f64::to_bits),
        expected.witness_b.map(f64::to_bits)
    );
}

fn assert_gap_bits_eq(actual: GapSample, expected: GapSample) {
    assert_eq!(actual.sum_lo.to_bits(), expected.sum_lo.to_bits());
    assert_eq!(actual.sum_hi.to_bits(), expected.sum_hi.to_bits());
    assert_eq!(
        actual.separation_upper.map(f64::to_bits),
        expected.separation_upper.map(f64::to_bits)
    );
    assert_eq!(
        actual.overlap_inradius.map(f64::to_bits),
        expected.overlap_inradius.map(f64::to_bits)
    );
    assert_eq!(
        actual.normal.map(|v| v.map(f64::to_bits)),
        expected.normal.map(|v| v.map(f64::to_bits))
    );
}

fn assert_codim_bits_eq(actual: CodimGap, expected: CodimGap) {
    assert_eq!(actual.lo.to_bits(), expected.lo.to_bits());
    assert_eq!(actual.hi.to_bits(), expected.hi.to_bits());
    assert_eq!(actual.verdict, expected.verdict);
}

fn enclosure_contains(outer: MomentEnclosure, inner: MomentEnclosure, label: &str) {
    assert!(
        outer.lo <= inner.lo && inner.hi <= outer.hi,
        "{label}: [{}, {}] must contain [{}, {}]",
        outer.lo,
        outer.hi,
        inner.lo,
        inner.hi
    );
}

fn moment_enclosures(m: &GeometricMoments) -> [MomentEnclosure; 10] {
    [
        m.volume,
        m.first[0],
        m.first[1],
        m.first[2],
        m.second.xx,
        m.second.yy,
        m.second.zz,
        m.second.xy,
        m.second.xz,
        m.second.yz,
    ]
}

fn assert_moment_bits_eq(actual: GeometricMoments, expected: GeometricMoments) {
    for (actual, expected) in moment_enclosures(&actual)
        .into_iter()
        .zip(moment_enclosures(&expected))
    {
        assert_eq!(actual.lo.to_bits(), expected.lo.to_bits());
        assert_eq!(actual.hi.to_bits(), expected.hi.to_bits());
    }
    assert_eq!(actual.h.to_bits(), expected.h.to_bits());
    assert_eq!(actual.sure_cells, expected.sure_cells);
    assert_eq!(actual.band_cells, expected.band_cells);
}

#[test]
fn fi_001_carriers_require_rigorous_receipts_and_compose_outward() {
    let zero = ContactInflation::exact_zero();
    assert_eq!(zero.radius().to_bits(), 0.0f64.to_bits());

    // The conservative upper endpoint, rather than the nominal QoI, is the
    // radius that must travel into contact bounds.
    let conversion = ContactInflation::from_conversion(&conversion_receipt(0.0625, 0.125))
        .expect("certified conversion receipt");
    let motion = ContactInflation::from_motion(&conversion_receipt(0.03125, 0.0625))
        .expect("certified motion-error receipt");
    let route = with_cx(|cx| routed_inflation(0.25, cx));
    assert_eq!(conversion.radius().to_bits(), 0.125f64.to_bits());
    assert_eq!(motion.radius().to_bits(), 0.0625f64.to_bits());
    assert_eq!(route.radius().to_bits(), 0.25f64.to_bits());
    assert_eq!(
        zero.compose(route)
            .expect("zero identity")
            .radius()
            .to_bits(),
        route.radius().to_bits()
    );
    assert_eq!(
        conversion
            .compose(zero)
            .expect("zero identity")
            .radius()
            .to_bits(),
        conversion.radius().to_bits()
    );

    let composed = conversion.compose(route).expect("finite radius sum");
    assert_eq!(
        composed.radius().to_bits(),
        (0.125f64 + 0.25).next_up().to_bits(),
        "positive carrier composition rounds outward exactly once"
    );
}

#[test]
fn fi_002_exact_zero_is_bit_neutral_for_every_query_family() {
    let zero = ContactInflation::exact_zero();
    let convex_a = convex_sphere(-1.0, 0.5);
    let convex_b = convex_sphere(1.0, 0.5);
    let chart_a = sphere_chart(-1.0, 0.5);
    let chart_b = sphere_chart(1.0, 0.5);
    let thickness_a = CodimThickness::new(0.125).expect("thickness");
    let thickness_b = CodimThickness::new(0.25).expect("thickness");
    let (a_pos, triangles) = tetra([0.0, 0.0, 0.0], 1.0);
    let (b_pos, _) = tetra([1.3, 0.0, 0.0], 1.0);
    let features_a = fs_query::FeatureComplex::from_triangles(&a_pos, &triangles).expect("a");
    let features_b = fs_query::FeatureComplex::from_triangles(&b_pos, &triangles).expect("b");
    let box_chart = BoxChart {
        aabb: Aabb::new(
            Point3::new(-0.25, -0.25, -0.25),
            Point3::new(0.25, 0.25, 0.25),
        ),
    };
    let domain = Aabb::new(Point3::new(-0.5, -0.5, -0.5), Point3::new(0.5, 0.5, 0.5));

    with_cx(|cx| {
        let convex = convex_separation(&convex_a, &convex_b, 256, cx).expect("convex");
        let convex_zero =
            convex_separation_with_inflation(&convex_a, &convex_b, 256, cx, zero, zero)
                .expect("zero-inflated convex");
        assert_convex_bits_eq(convex_zero, convex);

        let gap = ImplicitGapOracle::new(&chart_a, &chart_b)
            .expect("gap oracle")
            .gap_at(Point3::new(0.0, 0.0, 0.0), cx)
            .expect("gap");
        let gap_zero = ImplicitGapOracle::new_with_inflation(&chart_a, &chart_b, zero, zero)
            .expect("zero-inflated gap oracle")
            .gap_at(Point3::new(0.0, 0.0, 0.0), cx)
            .expect("zero-inflated gap");
        assert_gap_bits_eq(gap_zero, gap);

        let codim = codim_gap(1.0, 1.25, thickness_a, thickness_b).expect("codim");
        let codim_zero = codim_gap_with_inflation(1.0, 1.25, thickness_a, thickness_b, zero, zero)
            .expect("zero-inflated codim");
        assert_codim_bits_eq(codim_zero, codim);

        let ccd =
            ccd_candidates(&features_a, &features_b, 0.2, 0.2, 10_000, cx).expect("ccd candidates");
        let ccd_zero = ccd_candidates_with_inflation(
            &features_a,
            &features_b,
            0.2,
            0.2,
            10_000,
            cx,
            zero,
            zero,
        )
        .expect("zero-inflated ccd");
        assert_eq!(ccd_zero, ccd);

        let moments = geometric_moments(&box_chart, &domain, 0.125, cx).expect("moments");
        let moments_zero = geometric_moments_with_inflation(&box_chart, &domain, 0.125, cx, zero)
            .expect("zero-inflated moments");
        assert_moment_bits_eq(moments_zero, moments);
    });
}

#[test]
fn fi_003_convex_gap_and_codim_bounds_move_outward_monotonically() {
    let zero = ContactInflation::exact_zero();
    let small = certified_inflation(0.03125);
    let large = certified_inflation(0.25);
    let convex_a = convex_sphere(-1.0, 0.5);
    let convex_b = convex_sphere(1.0, 0.5);
    let outside_a = sphere_chart(-1.0, 0.5);
    let outside_b = sphere_chart(1.0, 0.5);
    let overlap_a = sphere_chart(-0.25, 1.0);
    let overlap_b = sphere_chart(0.25, 1.0);
    let midpoint = Point3::new(0.0, 0.0, 0.0);
    let thickness_a = CodimThickness::new(0.125).expect("thickness");
    let thickness_b = CodimThickness::new(0.25).expect("thickness");

    with_cx(|cx| {
        let nominal = convex_separation(&convex_a, &convex_b, 256, cx).expect("convex");
        let widened = convex_separation_with_inflation(&convex_a, &convex_b, 256, cx, large, zero)
            .expect("inflated convex");
        assert_eq!(
            widened.lo.to_bits(),
            (nominal.lo - large.radius()).next_down().max(0.0).to_bits()
        );
        assert_eq!(
            widened.hi.to_bits(),
            (nominal.hi + large.radius()).next_up().to_bits()
        );
        let convex_small =
            convex_separation_with_inflation(&convex_a, &convex_b, 256, cx, small, zero)
                .expect("small convex inflation");
        assert!(widened.lo <= convex_small.lo && widened.hi >= convex_small.hi);

        let nominal_gap = ImplicitGapOracle::new(&outside_a, &outside_b)
            .expect("gap oracle")
            .gap_at(midpoint, cx)
            .expect("gap");
        let large_gap = ImplicitGapOracle::new_with_inflation(&outside_a, &outside_b, large, zero)
            .expect("inflated gap oracle")
            .gap_at(midpoint, cx)
            .expect("inflated gap");
        assert_eq!(
            large_gap.sum_lo.to_bits(),
            (nominal_gap.sum_lo - large.radius()).next_down().to_bits()
        );
        assert_eq!(
            large_gap.sum_hi.to_bits(),
            (nominal_gap.sum_hi + large.radius()).next_up().to_bits()
        );
        assert_eq!(
            large_gap.separation_upper.map(f64::to_bits),
            Some(large_gap.sum_hi.to_bits())
        );
        let small_gap = ImplicitGapOracle::new_with_inflation(&outside_a, &outside_b, small, zero)
            .expect("small gap oracle")
            .gap_at(midpoint, cx)
            .expect("small gap");
        assert!(large_gap.sum_lo <= small_gap.sum_lo && large_gap.sum_hi >= small_gap.sum_hi);

        let nominal_overlap = ImplicitGapOracle::new(&overlap_a, &overlap_b)
            .expect("overlap oracle")
            .gap_at(midpoint, cx)
            .expect("overlap gap");
        let overlap_sample = overlap_a.eval(midpoint, cx);
        let overlap_enclosure = overlap_a.trace_value_enclosure(midpoint, &overlap_sample, cx);
        let expected_witness = (-(overlap_enclosure.hi + large.radius()).next_up()).next_down();
        let large_overlap =
            ImplicitGapOracle::new_with_inflation(&overlap_a, &overlap_b, large, zero)
                .expect("inflated overlap oracle")
                .gap_at(midpoint, cx)
                .expect("inflated overlap gap");
        assert_eq!(
            large_overlap.overlap_inradius.map(f64::to_bits),
            Some(expected_witness.to_bits()),
            "the active chart's overlap witness shrinks by its outward inflation"
        );
        assert!(
            large_overlap.overlap_inradius.expect("witness remains")
                < nominal_overlap.overlap_inradius.expect("nominal witness")
        );

        let codim_small =
            codim_gap_with_inflation(1.0, 1.25, thickness_a, thickness_b, small, zero)
                .expect("small codim inflation");
        let codim_large =
            codim_gap_with_inflation(1.0, 1.25, thickness_a, thickness_b, large, zero)
                .expect("large codim inflation");
        let codim_expected = codim_gap(
            (1.0 - large.radius()).next_down().max(0.0),
            (1.25 + large.radius()).next_up(),
            thickness_a,
            thickness_b,
        )
        .expect("manually widened distance enclosure");
        assert_codim_bits_eq(codim_large, codim_expected);
        assert!(codim_large.lo <= codim_small.lo && codim_large.hi >= codim_small.hi);
    });
}

#[test]
fn fi_004_ccd_and_moment_classification_are_monotone_in_inflation() {
    let zero = ContactInflation::exact_zero();
    let small = certified_inflation(0.03125);
    let large = certified_inflation(0.25);
    let (a_pos, triangles) = tetra([0.0, 0.0, 0.0], 1.0);
    let (b_pos, _) = tetra([1.3, 0.0, 0.0], 1.0);
    let features_a = fs_query::FeatureComplex::from_triangles(&a_pos, &triangles).expect("a");
    let features_b = fs_query::FeatureComplex::from_triangles(&b_pos, &triangles).expect("b");
    let chart = BoxChart {
        aabb: Aabb::new(
            Point3::new(-0.25, -0.25, -0.25),
            Point3::new(0.25, 0.25, 0.25),
        ),
    };
    let domain = Aabb::new(Point3::new(-0.5, -0.5, -0.5), Point3::new(0.5, 0.5, 0.5));

    with_cx(|cx| {
        let narrow = ccd_candidates_with_inflation(
            &features_a,
            &features_b,
            0.0,
            0.0,
            10_000,
            cx,
            small,
            zero,
        )
        .expect("narrow ccd window");
        let wide = ccd_candidates_with_inflation(
            &features_a,
            &features_b,
            0.0,
            0.0,
            10_000,
            cx,
            large,
            large,
        )
        .expect("wide ccd window");
        for pair in &narrow {
            assert!(wide.contains(pair), "wider CCD radius lost {pair:?}");
        }
        assert!(
            wide.contains(&(1, 0)),
            "0.5 total representation inflation must admit the 0.3-gap pair"
        );

        let nominal = geometric_moments_with_inflation(&chart, &domain, 0.125, cx, zero)
            .expect("nominal moments");
        let small_moments = geometric_moments_with_inflation(&chart, &domain, 0.125, cx, small)
            .expect("small moment inflation");
        let large_moments = geometric_moments_with_inflation(&chart, &domain, 0.125, cx, large)
            .expect("large moment inflation");

        for (index, (small_enclosure, large_enclosure)) in moment_enclosures(&small_moments)
            .into_iter()
            .zip(moment_enclosures(&large_moments))
            .enumerate()
        {
            enclosure_contains(
                large_enclosure,
                small_enclosure,
                &format!("moment component {index}"),
            );
        }
        assert!(large_moments.sure_cells <= small_moments.sure_cells);
        assert!(large_moments.band_cells >= small_moments.band_cells);
        assert!(small_moments.sure_cells <= nominal.sure_cells);
        assert!(small_moments.band_cells >= nominal.band_cells);
        assert!(
            large_moments.sure_cells < nominal.sure_cells,
            "the chosen bounded grid must observe stricter classification"
        );
    });
}
