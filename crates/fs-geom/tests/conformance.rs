//! fs-geom conformance suite (CONTRACT.md: any reimplementation must
//! pass). G0 trait laws on the fixture charts, the agreement checker's
//! detection acceptance, conversion receipts with empirical containment,
//! cancellation, and deterministic reports. JSON-line verdicts; seeded
//! cases carry seeds.

use asupersync::types::Budget;
use fs_evidence::ProvenanceHash;
use fs_exec::{CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::{BoxChart, LyingSphereChart, SphereChart, TorusChart};
use fs_geom::{
    Aabb, AgreementConfig, AgreementScope, AgreementStatus, AgreementUnknownReason, Axis, Chart,
    ClippedChart, Convert, ConvertDiag, ErrBudget, Point3, Region, SamplingDomain,
    SamplingDomainError, Vec3,
};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-geom/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\
         \"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "case {case}: {detail}");
}

struct Lcg(u64);

impl Lcg {
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / (1u64 << 53) as f64
    }
}

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x9E0,
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

fn charts() -> Vec<Box<dyn Chart>> {
    vec![
        Box::new(SphereChart {
            center: Point3::new(0.3, -0.2, 0.1),
            radius: 1.4,
        }),
        Box::new(BoxChart {
            aabb: Aabb::new(Point3::new(-1.0, -0.8, -1.2), Point3::new(0.9, 1.1, 0.7)),
        }),
        Box::new(TorusChart {
            center: Point3::new(0.0, 0.0, 0.0),
            major: 2.0,
            minor: 0.6,
        }),
    ]
}

#[test]
fn geo_001_g0_trait_laws_on_the_fixture_zoo() {
    const SEED: u64 = 0x0901_2026_0706_1AB5;
    let gate = CancelGate::new();
    let mut rng = Lcg(SEED);
    let mut checked = 0u64;
    with_cx(&gate, |cx| {
        for chart in charts() {
            let support = chart.support().inflate(2.0);
            for _ in 0..4_000 {
                let p = Point3::new(
                    support.min.x + (support.max.x - support.min.x) * rng.unit(),
                    support.min.y + (support.max.y - support.min.y) * rng.unit(),
                    support.min.z + (support.max.z - support.min.z) * rng.unit(),
                );
                let s = chart.eval(p, cx);
                // Law 1: inside ⇔ sd < 0.
                assert_eq!(chart.inside(p, cx), s.signed_distance < 0.0);
                // Law 2: support() actually bounds the region.
                if !chart.support().contains(p) {
                    assert!(
                        s.signed_distance > -1e-9,
                        "{}: negative sd outside support at {p:?}",
                        chart.name()
                    );
                }
                // Law 3: certified Lipschitz bound holds along a random ray.
                if let Some(l) = s.lipschitz {
                    let q = p.offset(Vec3::new(
                        (rng.unit() - 0.5) * 0.6,
                        (rng.unit() - 0.5) * 0.6,
                        (rng.unit() - 0.5) * 0.6,
                    ));
                    let sq = chart.eval(q, cx);
                    let step = q.delta_from(p).norm();
                    assert!(
                        (sq.signed_distance - s.signed_distance).abs() <= l * step + 1e-9,
                        "{}: Lipschitz {l} violated",
                        chart.name()
                    );
                }
                // Law 4: gradients (where claimed) match central FD and
                // have unit norm for exact SDFs.
                if let Some(g) = s.gradient {
                    assert!((g.norm() - 1.0).abs() < 1e-6, "{}", chart.name());
                    let h = 1e-6;
                    let fd = Vec3::new(
                        (chart
                            .eval(p.offset(Vec3::new(h, 0.0, 0.0)), cx)
                            .signed_distance
                            - chart
                                .eval(p.offset(Vec3::new(-h, 0.0, 0.0)), cx)
                                .signed_distance)
                            / (2.0 * h),
                        (chart
                            .eval(p.offset(Vec3::new(0.0, h, 0.0)), cx)
                            .signed_distance
                            - chart
                                .eval(p.offset(Vec3::new(0.0, -h, 0.0)), cx)
                                .signed_distance)
                            / (2.0 * h),
                        (chart
                            .eval(p.offset(Vec3::new(0.0, 0.0, h)), cx)
                            .signed_distance
                            - chart
                                .eval(p.offset(Vec3::new(0.0, 0.0, -h)), cx)
                                .signed_distance)
                            / (2.0 * h),
                    );
                    assert!(
                        g.sub_v(fd).norm() < 1e-4,
                        "{}: gradient vs FD",
                        chart.name()
                    );
                }
                checked += 1;
            }
        }
    });
    verdict(
        "geo-001",
        checked == 12_000,
        &format!(
            "trait laws hold over {checked} seeded queries on sphere/box/torus (seed {SEED:#x})"
        ),
    );
}

// Local extension: Vec3 difference (kept out of the tiny public surface).
trait VecExt {
    fn sub_v(self, o: Vec3) -> Vec3;
}

impl VecExt for Vec3 {
    fn sub_v(self, o: Vec3) -> Vec3 {
        Vec3::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

#[test]
fn geo_002_multi_chart_region_agrees_within_composed_bounds() {
    let sphere = SphereChart {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.5,
    };
    let gate = CancelGate::new();
    let (agreed, json_stable) = with_cx(&gate, |cx| {
        // Certified<T> is opaque (gp3.2.1): taking the value OUT is an
        // explicit downgrade to plain Evidence.
        let sampled = sphere
            .convert(ErrBudget { abs_sd_error: 0.08 }, cx)
            .expect("feasible budget")
            .into_evidence();
        let region = Region::from_chart(Arc::new(sphere), ProvenanceHash::of_bytes(b"exact"))
            .with_chart(Arc::new(sampled.value), sampled.provenance);
        let cfg = AgreementConfig::default();
        let r1 = region.check_agreement(&cfg, cx).expect("not cancelled");
        let r2 = region.check_agreement(&cfg, cx).expect("not cancelled");
        (
            r1.status == AgreementStatus::Agreed,
            r1.to_json() == r2.to_json(),
        )
    });
    verdict(
        "geo-002",
        agreed && json_stable,
        "exact sphere and its sampled conversion agree within composed declared bounds; \
         seeded reports replay identically (G5)",
    );
}

#[test]
fn geo_003_disagreement_is_detected_with_localized_diagnostics() {
    let sphere = SphereChart {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.5,
    };
    let gate = CancelGate::new();
    let report = with_cx(&gate, |cx| {
        let region = Region::from_chart(Arc::new(sphere), ProvenanceHash::of_bytes(b"honest"))
            .with_chart(
                Arc::new(LyingSphereChart { sphere, bias: 0.03 }),
                ProvenanceHash::of_bytes(b"liar"),
            );
        region
            .check_agreement(&AgreementConfig::default(), cx)
            .expect("not cancelled")
    });
    let localized = report.status == AgreementStatus::Disagreed
        && !report.disagreements.is_empty()
        && report.disagreements.iter().all(|d| {
            (d.gap - 0.03).abs() < 1e-9
                && (d.chart_a == "fixture/lying-sphere" || d.chart_b == "fixture/lying-sphere")
        });
    verdict(
        "geo-003",
        localized,
        &format!(
            "a 0.03 undeclared bias is caught and localized ({} diagnostics, worst excess \
             {:.4}); report: {}",
            report.disagreements.len(),
            report.worst_excess.expect("valid comparisons"),
            report.to_json()
        ),
    );
}

#[test]
fn geo_004_conversion_receipts_are_rigorous_and_refusals_teach() {
    const SEED: u64 = 0x0904_2026_0706_C0F0;
    let sphere = SphereChart {
        center: Point3::new(0.1, 0.2, -0.1),
        radius: 1.2,
    };
    let gate = CancelGate::new();
    let (contained, receipt_bound, refusal) = with_cx(&gate, |cx| {
        let certified = sphere
            .convert(ErrBudget { abs_sd_error: 0.05 }, cx)
            .expect("feasible");
        // Empirical containment: |sampled - exact| ≤ receipt bound over
        // seeded points inside the sampled box (G0 law of the receipt).
        let mut rng = Lcg(SEED);
        let box_ = certified.value.support();
        let mut worst = 0.0f64;
        for _ in 0..10_000 {
            let p = Point3::new(
                box_.min.x + (box_.max.x - box_.min.x) * rng.unit(),
                box_.min.y + (box_.max.y - box_.min.y) * rng.unit(),
                box_.min.z + (box_.max.z - box_.min.z) * rng.unit(),
            );
            let err = (certified.value.eval(p, cx).signed_distance
                - sphere.eval(p, cx).signed_distance)
                .abs();
            worst = worst.max(err);
        }
        let contained = worst <= certified.qoi;
        // Infeasible budget refuses BEFORE running, with ranked fixes.
        let refusal = sphere.convert(ErrBudget { abs_sd_error: 1e-6 }, cx);
        (contained, certified.qoi, refusal)
    });
    let teaches = matches!(&refusal, Err(ConvertDiag::BudgetInfeasible { .. }))
        && refusal
            .as_ref()
            .err()
            .is_some_and(|e| e.to_string().contains("Fixes (ranked)"));
    verdict(
        "geo-004",
        contained && teaches,
        &format!(
            "sampled-sdf receipt bound {receipt_bound:.4} contains the empirical error over \
             10k seeded points (seed {SEED:#x}); infeasible budgets refuse with ranked fixes"
        ),
    );
}

#[test]
fn geo_004c_zero_error_budget_refuses_instead_of_overflowing() {
    // Regression: abs_sd_error = 0.0 asks for infinite resolution, so
    // h_needed = 0 and edge / h_needed = +∞. The old
    // `(+∞).ceil() as u32 + 1` OVERFLOWED — a debug/test panic (release: wraps
    // to need_resolution = 0, a nonsensical diagnostic) — defeating the
    // BudgetInfeasible refusal and breaking the "no panic crosses the
    // boundary" contract. It must now refuse cleanly with a sane resolution.
    let sphere = SphereChart {
        center: Point3::new(0.1, 0.2, -0.1),
        radius: 1.2,
    };
    let gate = CancelGate::new();
    let refusal = with_cx(&gate, |cx| {
        sphere.convert(ErrBudget { abs_sd_error: 0.0 }, cx)
    });
    let ok = matches!(&refusal, Err(ConvertDiag::InvalidBudget { .. }));
    assert!(
        ok,
        "a zero error budget must refuse as invalid before support/evaluation/count math: \
         {refusal:?}"
    );
    verdict(
        "geo-004c",
        ok,
        "zero error budget refuses before evaluation or integer count arithmetic",
    );
}

#[test]
fn geo_004j_translated_tiny_exact_domain_keeps_outward_enclosures() {
    // Regression for proof arithmetic at a scale where adding an ideal grid
    // fraction to the translated minimum loses many low bits. The conversion
    // must interpolate against the ACTUAL retained nodes and every published
    // enclosure must contain an independently evaluated exact-distance chart,
    // including the outward extension beyond sampled support.
    let base = 1.0e12;
    let width = 0.031_25;
    let source = BoxChart {
        aabb: Aabb::new(
            Point3::new(base, base + width, base - width),
            Point3::new(base + width, base + 2.0 * width, base),
        ),
    };
    let budget = 0.05;
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let converted = source
            .convert(
                ErrBudget {
                    abs_sd_error: budget,
                },
                cx,
            )
            .expect("translated tiny exact-distance domain is representable");
        assert!(converted.qoi <= budget);
        for axis in [Axis::X, Axis::Y, Axis::Z] {
            let nodes = converted.value.axis_nodes(axis);
            assert_eq!(nodes.len(), converted.value.resolution() as usize);
            assert!(nodes.windows(2).all(|pair| pair[0] < pair[1]));
        }
        let support = converted.value.support();
        let point_at = |tx: f64, ty: f64, tz: f64| {
            Point3::new(
                support.min.x + (support.max.x - support.min.x) * tx,
                support.min.y + (support.max.y - support.min.y) * ty,
                support.min.z + (support.max.z - support.min.z) * tz,
            )
        };
        let mut points = Vec::new();
        for t in [0.0, 0.031_25, 0.2, 0.5, 0.812_5, 1.0] {
            points.push(point_at(t, 1.0 - t, (3.0 * t).fract()));
        }
        points.extend([
            support.min,
            support.max,
            Point3::new(
                support.max.x + 0.02,
                f64::midpoint(support.min.y, support.max.y),
                f64::midpoint(support.min.z, support.max.z),
            ),
            Point3::new(
                f64::midpoint(support.min.x, support.max.x),
                support.min.y - 0.02,
                support.max.z + 0.02,
            ),
        ]);

        for point in points {
            let independent = source.eval(point, cx);
            let sampled = converted.value.eval(point, cx);
            assert_eq!(sampled.error.kind, fs_evidence::NumericalKind::Enclosure);
            assert!(sampled.signed_distance.is_finite(), "at {point:?}");
            assert!(
                sampled.error.lo <= independent.error.lo
                    && independent.error.hi <= sampled.error.hi,
                "sampled interval [{}, {}] failed to contain independent source interval \
                 [{}, {}] at {point:?}",
                sampled.error.lo,
                sampled.error.hi,
                independent.error.lo,
                independent.error.hi
            );
        }
        let non_finite_query = converted.value.eval(Point3::new(f64::NAN, base, base), cx);
        assert!(non_finite_query.signed_distance.is_nan());
        assert_eq!(
            non_finite_query.error.kind,
            fs_evidence::NumericalKind::NoClaim
        );
    });
}

#[test]
fn geo_004k_conversion_refuses_a_nonrepresentable_dense_grid() {
    // Around 1e16 adjacent f64 values are two units apart. Padding this
    // one-ulp box for a unit budget leaves only three representable coordinates
    // per axis, while the full-cell-diagonal proof requires a denser grid.
    // Coincident ideal nodes must be a typed refusal, never silently sampled.
    let min = 1.0e16_f64;
    let max = min.next_up();
    let source = BoxChart {
        aabb: Aabb::new(Point3::new(min, min, min), Point3::new(max, max, max)),
    };
    let gate = CancelGate::new();
    let refusal = with_cx(&gate, |cx| {
        source.convert(ErrBudget { abs_sd_error: 1.0 }, cx)
    });
    assert!(
        matches!(
            refusal,
            Err(ConvertDiag::UnrepresentableGrid {
                resolution,
                min_bits,
                max_bits,
                ..
            }) if resolution > 3 && min_bits < max_bits
        ),
        "a dense grid with coincident representable nodes must refuse: {refusal:?}"
    );
}

struct UnboundedHalfSpace;

impl Chart for UnboundedHalfSpace {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> fs_geom::ChartSample {
        fs_geom::ChartSample {
            signed_distance: x.x,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: fs_evidence::NumericalCertificate::exact(x.x),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::WHOLE_SPACE
    }

    fn name(&self) -> &'static str {
        "test/unbounded-half-space"
    }
}

#[test]
fn geo_004d_region_agreement_requires_an_explicit_finite_scope() {
    let region = Region::from_chart(
        Arc::new(UnboundedHalfSpace),
        ProvenanceHash::of_bytes(b"unbounded-a"),
    )
    .with_chart(
        Arc::new(UnboundedHalfSpace),
        ProvenanceHash::of_bytes(b"unbounded-b"),
    );
    let clip = Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0));
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let unresolved = region
            .check_agreement(&AgreementConfig::default(), cx)
            .expect("not cancelled");
        assert_eq!(unresolved.status, AgreementStatus::Unknown);
        assert_eq!(unresolved.scope, AgreementScope::GlobalSupport);
        assert!(unresolved.sampling_domain.is_none());
        assert!(unresolved.unknowns.iter().any(|unknown| matches!(
            unknown.reason,
            AgreementUnknownReason::SamplingDomain(SamplingDomainError::UnboundedSupport { .. })
        )));

        let local = region
            .check_agreement(
                &AgreementConfig {
                    sampling_clip: Some(clip),
                    ..AgreementConfig::default()
                },
                cx,
            )
            .expect("not cancelled");
        assert_eq!(local.status, AgreementStatus::Agreed);
        assert_eq!(local.scope, AgreementScope::ExplicitClip);
        assert_eq!(local.sampling_domain, Some(clip));
    });
}

#[test]
fn geo_004e_unbounded_conversion_requires_a_geometric_clip() {
    let source = UnboundedHalfSpace;
    let clip = Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0));
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let global: Result<fs_evidence::Certified<fs_geom::SampledSdf>, _> =
            source.convert(ErrBudget { abs_sd_error: 0.25 }, cx);
        assert!(matches!(
            global,
            Err(ConvertDiag::SamplingDomain(
                SamplingDomainError::UnboundedSupport { .. }
            ))
        ));

        let clipped = ClippedChart::new(&source, clip).expect("finite clip is admissible");
        assert_eq!(clipped.support(), clip);
        assert!(
            clipped
                .eval(Point3::new(-0.5, 0.0, 0.0), cx)
                .signed_distance
                < 0.0
        );
        let outside = clipped.eval(Point3::new(-2.0, 0.0, 0.0), cx);
        assert!(
            outside.signed_distance > 0.0,
            "clip participates in the field"
        );
        assert_eq!(outside.error.kind, fs_evidence::NumericalKind::NoClaim);

        let local: fs_evidence::Evidence<fs_geom::SampledSdf> = source
            .convert_clipped(ErrBudget { abs_sd_error: 0.25 }, clip, cx)
            .expect("finite clipped composite is convertible");
        assert!(local.value.support().is_finite());
        assert!(local.value.nominal_field_bound().is_finite());
        assert_eq!(
            local.value.abstract_distance_kind(),
            fs_evidence::NumericalKind::NoClaim
        );
        assert_eq!(local.value.abstract_distance_bound(), None);
        assert_eq!(local.numerical.kind, fs_evidence::NumericalKind::NoClaim);
        assert_eq!(
            local.value.eval(Point3::new(-0.5, 0.0, 0.0), cx).error.kind,
            fs_evidence::NumericalKind::NoClaim
        );
        assert!(matches!(
            local.clone().certified(),
            Err(fs_evidence::CertifyError::NotRigorous {
                kind: fs_evidence::NumericalKind::NoClaim
            })
        ));
    });
}

struct BoundedNoClaim;

impl Chart for BoundedNoClaim {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> fs_geom::ChartSample {
        fs_geom::ChartSample {
            signed_distance: x.x,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: fs_evidence::NumericalCertificate::no_claim(),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "test/bounded-no-claim"
    }
}

struct BoundedLocalOnly;

impl Chart for BoundedLocalOnly {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> fs_geom::ChartSample {
        fs_geom::ChartSample {
            signed_distance: x.x,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: fs_evidence::NumericalCertificate::enclosure(x.x, x.x),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "test/bounded-local-only"
    }
}

#[test]
fn geo_004i_bounded_no_claim_source_cannot_be_laundered_by_sampling() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let source = BoundedNoClaim;
        let refusal = source.convert(ErrBudget { abs_sd_error: 0.5 }, cx);
        assert!(matches!(
            refusal,
            Err(ConvertDiag::NoAbstractDistanceClaim {
                kind: fs_evidence::NumericalKind::NoClaim
            })
        ));

        let evidence = source
            .convert_with_domain(ErrBudget { abs_sd_error: 0.5 }, None, cx)
            .expect("plain field evidence remains available");
        assert_eq!(evidence.numerical.kind, fs_evidence::NumericalKind::NoClaim);
        assert_eq!(
            evidence
                .value
                .eval(Point3::new(0.0, 0.0, 0.0), cx)
                .error
                .kind,
            fs_evidence::NumericalKind::NoClaim
        );

        let local_only = BoundedLocalOnly;
        assert!(matches!(
            local_only.convert(ErrBudget { abs_sd_error: 0.5 }, cx),
            Err(ConvertDiag::NoAbstractDistanceClaim {
                kind: fs_evidence::NumericalKind::Estimate
            })
        ));
        let local_evidence = local_only
            .convert_with_domain(ErrBudget { abs_sd_error: 0.5 }, None, cx)
            .expect("nominal estimate remains available");
        assert_eq!(
            local_evidence.numerical.kind,
            fs_evidence::NumericalKind::Estimate
        );
        assert_eq!(
            local_evidence.value.abstract_distance_kind(),
            fs_evidence::NumericalKind::Estimate
        );
    });
}

#[test]
fn geo_004f_sampling_domain_rejects_nan_and_span_overflow() {
    let malformed = Aabb::new(
        Point3::new(f64::NAN, -1.0, -1.0),
        Point3::new(1.0, 1.0, 1.0),
    );
    assert!(matches!(
        SamplingDomain::admit(malformed, None),
        Err(SamplingDomainError::InvalidSupport { .. })
    ));

    let overflowing = Aabb::new(
        Point3::new(-f64::MAX, -1.0, -1.0),
        Point3::new(f64::MAX, 1.0, 1.0),
    );
    assert!(matches!(
        SamplingDomain::admit(overflowing, None),
        Err(SamplingDomainError::NonFiniteSpan { .. })
    ));
}

/// An HONEST chart whose LOCAL Lipschitz varies: 1 near the origin, 50 in a
/// fast-oscillating shell beyond `x = 1`. It never lies about any local bound.
struct VaryingLipschitz;

impl Chart for VaryingLipschitz {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> fs_geom::ChartSample {
        let (sd, lip) = if x.x <= 1.0 {
            (x.x, 1.0)
        } else {
            (1.0 + (50.0 * (x.x - 1.0)).sin(), 50.0)
        };
        fs_geom::ChartSample {
            signed_distance: sd,
            gradient: None,
            lipschitz: Some(lip),
            error: fs_evidence::NumericalCertificate::exact(sd),
        }
    }
    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-2.0, -2.0, -2.0), Point3::new(2.0, 2.0, 2.0))
    }
    fn name(&self) -> &'static str {
        "varying-lipschitz-probe"
    }
}

#[test]
fn geo_004b_convert_bounds_global_lipschitz_not_the_center() {
    // The box-center probe sees local Lipschitz 1, but the sampling grid lands
    // in the 50-Lipschitz shell (|x.x| up to 2). Sampling only the center would
    // ship a sampled-SDF receipt understating the trilinear error ~18x; the
    // grid-max bound must instead REFUSE this budget rather than certify a false
    // enclosure (bead obnw F2). (Resolution stays well under the cap, so the
    // grid-Lipschitz check — not the resolution cap — is what refuses.)
    let gate = CancelGate::new();
    let refusal = with_cx(&gate, |cx| {
        VaryingLipschitz.convert(ErrBudget { abs_sd_error: 0.5 }, cx)
    });
    assert!(
        matches!(refusal, Err(ConvertDiag::BudgetInfeasible { .. })),
        "convert must refuse (grid Lipschitz 50 >> center 1), got {refusal:?}"
    );
}

struct NonFiniteGridChart;

impl Chart for NonFiniteGridChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> fs_geom::ChartSample {
        let signed_distance = if x == Point3::new(0.0, 0.0, 0.0) {
            0.0
        } else {
            f64::NAN
        };
        fs_geom::ChartSample {
            signed_distance,
            gradient: None,
            lipschitz: Some(1.0),
            error: fs_evidence::NumericalCertificate::no_claim(),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "test/non-finite-grid"
    }
}

#[test]
fn geo_004g_conversion_refuses_non_finite_grid_values() {
    let gate = CancelGate::new();
    let refusal = with_cx(&gate, |cx| {
        NonFiniteGridChart.convert(ErrBudget { abs_sd_error: 0.5 }, cx)
    });
    assert!(
        matches!(
            refusal,
            Err(ConvertDiag::NonFiniteSignedDistance { value_bits, .. })
                if value_bits == f64::NAN.to_bits()
        ),
        "a non-finite source field must never produce a certified conversion: {refusal:?}"
    );
}

struct CancellingGridChart<'a> {
    gate: &'a CancelGate,
    evaluations: AtomicUsize,
}

impl Chart for CancellingGridChart<'_> {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> fs_geom::ChartSample {
        if self.evaluations.fetch_add(1, Ordering::Relaxed) == 1 {
            self.gate.request();
        }
        fs_geom::ChartSample {
            signed_distance: x.x,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: fs_evidence::NumericalCertificate::exact(x.x),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "test/cancelling-grid"
    }
}

#[test]
fn geo_004h_conversion_grid_polls_cancellation_directly() {
    let gate = CancelGate::new();
    let chart = CancellingGridChart {
        gate: &gate,
        evaluations: AtomicUsize::new(0),
    };
    let refusal = with_cx(&gate, |cx| {
        chart.convert(ErrBudget { abs_sd_error: 0.5 }, cx)
    });
    assert!(
        matches!(
            refusal,
            Err(ConvertDiag::Cancelled {
                stage: "sampling-grid",
                completed_samples: 0,
            })
        ),
        "the converter itself must observe cancellation triggered by a non-polling source: \
         {refusal:?}"
    );
}

#[test]
fn geo_005_geometry_is_cancellable() {
    let sphere = SphereChart {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let gate = CancelGate::new();
    gate.request();
    let outcome = with_cx(&gate, |cx| {
        let region = Region::from_chart(Arc::new(sphere), ProvenanceHash::of_bytes(b"a"))
            .with_chart(Arc::new(sphere), ProvenanceHash::of_bytes(b"b"));
        region.check_agreement(&AgreementConfig::default(), cx)
    });
    verdict(
        "geo-005",
        outcome == Err(Cancelled),
        "agreement checking observes a pre-requested gate and returns Cancelled promptly (P7)",
    );
}
