//! Semantic-diff conformance (the lmp4.10 bead; runs under the
//! `semantic-diff` feature). Acceptance: on branch pairs with known
//! causal edits the diff localizes field differences to the correct
//! region/quantity and attributes them to the correct ops — a TWO-edit
//! pair attributes BOTH edits ranked by measured contribution; entities
//! without stable IDs degrade to a FLAGGED geometric fallback and the R3
//! fraction is measured; ID stability survives topology-changing edits;
//! the diff is invariant under re-indexing and rigid motion (G3).
#![cfg(feature = "semantic-diff")]

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::diff::{IdentifiedPatch, SemanticDiffError, semantic_diff, semantic_diff_clipped};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Aabb, Chart, EntityId, IdTransform, IdentityMap, Point3, SamplingDomainError, Vec3};
use std::sync::atomic::{AtomicUsize, Ordering};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-geom/diff\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    with_gate_cx(&gate, f)
}

fn with_gate_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
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

fn sphere(cx: f64, r: f64) -> SphereChart {
    SphereChart {
        center: Point3::new(cx, 0.0, 0.0),
        radius: r,
    }
}

const HULL: EntityId = EntityId(7);
const KEEL: EntityId = EntityId(8);

struct UnboundedPlane;

impl Chart for UnboundedPlane {
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
        "test/unbounded-plane"
    }
}

#[test]
fn df_000_unbounded_comparison_requires_an_explicit_clip() {
    with_cx(|cx| {
        let plane = UnboundedPlane;
        let a = [IdentifiedPatch {
            id: Some(HULL),
            chart: &plane,
        }];
        let b = [IdentifiedPatch {
            id: Some(HULL),
            chart: &plane,
        }];
        let identity = IdentityMap::new();
        assert!(matches!(
            semantic_diff(&a, &b, &identity, &[], &[], 1e-6, cx),
            Err(SemanticDiffError::SamplingDomain {
                source: SamplingDomainError::UnboundedSupport { .. },
                ..
            })
        ));
        let clip = Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0));
        let local = semantic_diff_clipped(&a, &b, &identity, &[], &[], 1e-6, clip, cx)
            .expect("finite clip admits unbounded field comparison");
        assert!(local.objects.is_empty());
        assert_eq!(local.sampling_clip, Some(clip));
    });
}

struct CancellingDiffChart<'a> {
    gate: &'a CancelGate,
    evaluations: AtomicUsize,
}

impl Chart for CancellingDiffChart<'_> {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> fs_geom::ChartSample {
        if self.evaluations.fetch_add(1, Ordering::Relaxed) == 0 {
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
        "test/cancelling-diff"
    }
}

#[test]
fn df_000b_field_sampler_polls_cancellation_directly() {
    let gate = CancelGate::new();
    let chart = CancellingDiffChart {
        gate: &gate,
        evaluations: AtomicUsize::new(0),
    };
    with_gate_cx(&gate, |cx| {
        let a = [IdentifiedPatch {
            id: Some(HULL),
            chart: &chart,
        }];
        let b = [IdentifiedPatch {
            id: Some(HULL),
            chart: &chart,
        }];
        let refusal = semantic_diff(&a, &b, &IdentityMap::new(), &[], &[], 1e-6, cx);
        assert!(
            matches!(
                refusal,
                Err(SemanticDiffError::Cancelled {
                    stage: "field-sampling",
                    entity: Some(HULL),
                    completed_draws: 0,
                })
            ),
            "semantic diff must poll independently of chart implementations: {refusal:?}"
        );
    });
}

struct MalformedDiffChart {
    value: f64,
}

impl Chart for MalformedDiffChart {
    fn eval(&self, _x: Point3, _cx: &Cx<'_>) -> fs_geom::ChartSample {
        fs_geom::ChartSample {
            signed_distance: self.value,
            gradient: None,
            lipschitz: None,
            error: fs_evidence::NumericalCertificate::no_claim(),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "test/malformed-diff"
    }
}

#[test]
fn df_000c_invalid_tolerance_and_nonfinite_fields_refuse_instead_of_false_clean() {
    with_cx(|cx| {
        let finite = sphere(0.0, 1.0);
        let invalid = MalformedDiffChart { value: f64::NAN };
        let world_a = [IdentifiedPatch {
            id: Some(HULL),
            chart: &invalid,
        }];
        let world_b = [IdentifiedPatch {
            id: Some(HULL),
            chart: &finite,
        }];
        assert!(matches!(
            semantic_diff(&world_a, &world_b, &IdentityMap::new(), &[], &[], 1e-6, cx,),
            Err(SemanticDiffError::InvalidFieldSample {
                entity: Some(HULL),
                side: fs_geom::diff::SemanticDiffSide::A,
                ..
            })
        ));
        assert_eq!(
            semantic_diff(
                &world_b,
                &world_b,
                &IdentityMap::new(),
                &[],
                &[],
                f64::NAN,
                cx,
            ),
            Err(SemanticDiffError::InvalidTolerance {
                value_bits: f64::NAN.to_bits(),
            })
        );
    });
}

#[test]
fn df_001_two_edit_pair_attributes_both_ranked() {
    with_cx(|cx| {
        // World A: hull sphere r=1.0 at x=0; keel sphere r=0.5 at x=3
        // (untouched control).
        let hull_a = sphere(0.0, 1.0);
        let keel = sphere(3.0, 0.5);
        // Edit op 101 (small): radius 1.0 -> 1.004.
        let hull_g1 = sphere(0.0, 1.004);
        // Edit op 102 (large): radius 1.004 -> 1.05.
        let hull_b = sphere(0.0, 1.05);
        let world_a = vec![
            IdentifiedPatch {
                id: Some(HULL),
                chart: &hull_a,
            },
            IdentifiedPatch {
                id: Some(KEEL),
                chart: &keel,
            },
        ];
        let gen1 = vec![
            IdentifiedPatch {
                id: Some(HULL),
                chart: &hull_g1,
            },
            IdentifiedPatch {
                id: Some(KEEL),
                chart: &keel,
            },
        ];
        let gen2 = vec![
            IdentifiedPatch {
                id: Some(HULL),
                chart: &hull_b,
            },
            IdentifiedPatch {
                id: Some(KEEL),
                chart: &keel,
            },
        ];
        let world_b = vec![
            IdentifiedPatch {
                id: Some(HULL),
                chart: &hull_b,
            },
            IdentifiedPatch {
                id: Some(KEEL),
                chart: &keel,
            },
        ];
        let mut identity = IdentityMap::new();
        identity.record(101, vec![IdTransform::Preserved(HULL)]);
        identity.record(102, vec![IdTransform::Preserved(HULL)]);
        let report = semantic_diff(
            &world_a,
            &world_b,
            &identity,
            &[101, 102],
            &[gen1, gen2],
            1e-6,
            cx,
        )
        .expect("bounded semantic-diff domain");
        assert_eq!(report.sampling_clip, None);
        // Exactly one finding: the hull. The keel (untouched) is quiet.
        assert_eq!(report.objects.len(), 1, "only the edited entity differs");
        let obj = &report.objects[0];
        assert_eq!(obj.entity, Some(HULL), "localized to the right entity");
        assert_eq!(obj.quantity, "signed-distance");
        assert!(
            obj.magnitude > 0.04 && obj.magnitude < 0.06,
            "total magnitude ~ 0.05: {}",
            obj.magnitude
        );
        // BOTH edits attributed, RANKED by measured contribution:
        // op 102 (0.046) above op 101 (0.004).
        assert_eq!(obj.causes.len(), 2, "both causal edits present");
        assert_eq!(obj.causes[0].0, 102, "the larger edit ranks first");
        assert_eq!(obj.causes[1].0, 101);
        assert!(
            obj.causes[0].1 > 10.0 * obj.causes[1].1,
            "contributions measured: {:?}",
            obj.causes
        );
        assert!(obj.attributed);
        assert!(
            report.fallback_fraction.abs() < f64::EPSILON,
            "no fallbacks here"
        );
        verdict(
            "df-001",
            "two-edit pair: localized to the hull, both ops attributed, ranked 102 > 101 \
             by measured contribution",
        );
    });
}

#[test]
fn df_002_fallback_is_flagged_and_measured() {
    with_cx(|cx| {
        // One identified pair + one LEGACY (unidentified) pair that
        // genuinely differs.
        let a1 = sphere(0.0, 1.0);
        let b1 = sphere(0.0, 1.0);
        let a2 = sphere(3.0, 0.5);
        let b2 = sphere(3.0, 0.52);
        let world_a = vec![
            IdentifiedPatch {
                id: Some(HULL),
                chart: &a1,
            },
            IdentifiedPatch {
                id: None,
                chart: &a2,
            },
        ];
        let world_b = vec![
            IdentifiedPatch {
                id: Some(HULL),
                chart: &b1,
            },
            IdentifiedPatch {
                id: None,
                chart: &b2,
            },
        ];
        let identity = IdentityMap::new();
        let report = semantic_diff(&world_a, &world_b, &identity, &[], &[], 1e-6, cx)
            .expect("bounded semantic-diff domain");
        // The legacy difference is found but FLAGGED unattributed.
        let fallback: Vec<_> = report.objects.iter().filter(|o| !o.attributed).collect();
        assert_eq!(fallback.len(), 1, "the legacy pair's difference is found");
        assert!(fallback[0].entity.is_none());
        assert!(
            fallback[0].causes.is_empty(),
            "no causal claim without identity"
        );
        // The R3 early-warning metric: 1 fallback of 2 comparisons.
        assert!(
            (report.fallback_fraction - 0.5).abs() < 1e-12,
            "fallback fraction measured: {}",
            report.fallback_fraction
        );
        verdict(
            "df-002",
            "unidentified entities degrade to flagged geometric comparison; R3 fraction \
             = 0.5 measured",
        );
    });
}

#[test]
fn df_003_id_stability_through_topology_changing_edits() {
    // The R3 stress battery on the IdentityMap: replace, split, merge —
    // attribution walks ancestry across all of them.
    let mut identity = IdentityMap::new();
    let (a, b, c, d, e) = (
        EntityId(1),
        EntityId(2),
        EntityId(3),
        EntityId(4),
        EntityId(5),
    );
    identity.record(10, vec![IdTransform::Created(a)]);
    identity.record(20, vec![IdTransform::Replaced(a, b)]); // re-fit
    identity.record(30, vec![IdTransform::Split(b, vec![c, d])]); // boolean cut
    identity.record(40, vec![IdTransform::Merged(vec![c, d], e)]); // weld
    identity.record(50, vec![IdTransform::Preserved(e)]); // param tweak
    // Ops touching the FINAL entity include the whole ancestry chain.
    assert_eq!(
        identity.ops_touching(e),
        vec![10, 20, 30, 40, 50],
        "attribution walks replace/split/merge ancestry"
    );
    // An unrelated entity touches nothing.
    assert!(identity.ops_touching(EntityId(99)).is_empty());
    // A mid-chain entity sees its history but not later unrelated ops.
    let mid = identity.ops_touching(c);
    assert!(mid.contains(&30) && mid.contains(&10), "ancestry: {mid:?}");
    verdict(
        "df-003",
        "IDs survive replace/split/merge; ops_touching returns the full ancestry chain",
    );
}

#[test]
fn df_004_invariance_reindex_and_rigid_motion() {
    with_cx(|cx| {
        let hull_a = sphere(0.0, 1.0);
        let hull_b = sphere(0.0, 1.02);
        let keel = sphere(3.0, 0.5);
        let mut identity = IdentityMap::new();
        identity.record(7, vec![IdTransform::Preserved(HULL)]);
        // Baseline order.
        let wa1 = vec![
            IdentifiedPatch {
                id: Some(HULL),
                chart: &hull_a,
            },
            IdentifiedPatch {
                id: Some(KEEL),
                chart: &keel,
            },
        ];
        let wb1 = vec![
            IdentifiedPatch {
                id: Some(HULL),
                chart: &hull_b,
            },
            IdentifiedPatch {
                id: Some(KEEL),
                chart: &keel,
            },
        ];
        // Re-indexed order (patches swapped): ID keying makes it exact.
        let wa2 = vec![
            IdentifiedPatch {
                id: Some(KEEL),
                chart: &keel,
            },
            IdentifiedPatch {
                id: Some(HULL),
                chart: &hull_a,
            },
        ];
        let wb2 = vec![
            IdentifiedPatch {
                id: Some(KEEL),
                chart: &keel,
            },
            IdentifiedPatch {
                id: Some(HULL),
                chart: &hull_b,
            },
        ];
        let r1 = semantic_diff(&wa1, &wb1, &identity, &[7], &[], 1e-6, cx)
            .expect("bounded semantic-diff domain");
        let r2 = semantic_diff(&wa2, &wb2, &identity, &[7], &[], 1e-6, cx)
            .expect("bounded semantic-diff domain");
        assert_eq!(r1.objects.len(), 1);
        assert_eq!(
            r1.objects[0].magnitude.to_bits(),
            r2.objects[0].magnitude.to_bits(),
            "re-indexing invariance is exact (ID-keyed, geometry-seeded)"
        );
        // Rigid motion: translate BOTH worlds by the same offset —
        // magnitude agrees to tolerance (samples move with the boxes).
        let hull_a_m = sphere(10.0, 1.0);
        let hull_b_m = sphere(10.0, 1.02);
        let keel_m = sphere(13.0, 0.5);
        let wam = vec![
            IdentifiedPatch {
                id: Some(HULL),
                chart: &hull_a_m,
            },
            IdentifiedPatch {
                id: Some(KEEL),
                chart: &keel_m,
            },
        ];
        let wbm = vec![
            IdentifiedPatch {
                id: Some(HULL),
                chart: &hull_b_m,
            },
            IdentifiedPatch {
                id: Some(KEEL),
                chart: &keel_m,
            },
        ];
        let rm = semantic_diff(&wam, &wbm, &identity, &[7], &[], 1e-6, cx)
            .expect("bounded semantic-diff domain");
        assert_eq!(rm.objects.len(), 1);
        assert!(
            (rm.objects[0].magnitude - r1.objects[0].magnitude).abs() < 5e-3,
            "rigid motion preserves the physics diff: {} vs {}",
            rm.objects[0].magnitude,
            r1.objects[0].magnitude
        );
        verdict(
            "df-004",
            "re-index invariance bitwise; translation invariance to 5e-3 (G3 metamorphic)",
        );
    });
}

#[test]
fn df_005_created_deleted_and_filtering() {
    with_cx(|cx| {
        let hull_a = sphere(0.0, 1.0);
        let hull_b = sphere(0.0, 1.05);
        let old_fin = sphere(6.0, 0.3);
        let new_wing = sphere(9.0, 0.4);
        let world_a = vec![
            IdentifiedPatch {
                id: Some(HULL),
                chart: &hull_a,
            },
            IdentifiedPatch {
                id: Some(EntityId(20)),
                chart: &old_fin,
            },
        ];
        let world_b = vec![
            IdentifiedPatch {
                id: Some(HULL),
                chart: &hull_b,
            },
            IdentifiedPatch {
                id: Some(EntityId(21)),
                chart: &new_wing,
            },
        ];
        let mut identity = IdentityMap::new();
        identity.record(1, vec![IdTransform::Preserved(HULL)]);
        let report = semantic_diff(&world_a, &world_b, &identity, &[1], &[], 1e-6, cx)
            .expect("bounded semantic-diff domain");
        assert_eq!(report.only_a, vec![EntityId(20)], "deleted entity reported");
        assert_eq!(report.only_b, vec![EntityId(21)], "created entity reported");
        // Filtering: magnitude floor and region window.
        let all = report.filter(None, Some("signed-distance"), 0.0);
        assert_eq!(all.len(), 1);
        assert!(
            report.filter(None, None, 1.0).is_empty(),
            "magnitude floor filters"
        );
        let far_window = fs_geom::Aabb::new(
            Point3::new(100.0, 100.0, 100.0),
            Point3::new(101.0, 101.0, 101.0),
        );
        assert!(
            report.filter(Some(&far_window), None, 0.0).is_empty(),
            "region window filters"
        );
        verdict(
            "df-005",
            "created/deleted entities reported; region/quantity/magnitude filters work",
        );
    });
}
