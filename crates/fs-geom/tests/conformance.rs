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
    Aabb, AgreementConfig, Chart, Convert, ConvertDiag, ErrBudget, Point3, Region, Vec3,
};
use std::sync::Arc;

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
        let sampled = sphere
            .convert(ErrBudget { abs_sd_error: 0.08 }, cx)
            .expect("feasible budget");
        let region = Region::from_chart(Arc::new(sphere), ProvenanceHash::of_bytes(b"exact"))
            .with_chart(Arc::new(sampled.value), sampled.provenance);
        let cfg = AgreementConfig::default();
        let r1 = region.check_agreement(&cfg, cx).expect("not cancelled");
        let r2 = region.check_agreement(&cfg, cx).expect("not cancelled");
        (r1.agreed, r1.to_json() == r2.to_json())
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
    let localized = !report.agreed
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
            report.worst_gap,
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
