//! Sheaf-certificate conformance (the wqd.13 bead). Acceptance: verdicts
//! correct with correct localization on seeded watertight/leaky
//! multi-chart fixtures; interval verification SOUND (no PASS on a truly
//! leaky retained sample); δδ = 0 bitwise; the retained two-chart mismatch
//! bound is invariant under swapping chart order (exact), and the fixture is
//! stable under rigid motion (tolerance-level); the adversarial seam zoo
//! (near-tangent, T-junction)
//! behaves honestly; the graph-gauge diagnostic remains explicitly separate
//! from the feature-gated cohomology classifier used by merge semantics.

use asupersync::types::Budget;
use fs_evidence::NumericalKind;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::{BoxChart, SphereChart};
use fs_geom::{
    Aabb, Chart, Interface, InterfaceSample, OUTSIDE_RAY_MAX_EVALUATIONS, OutsideRaySampleError,
    OutsideRaySampleReport, Point3, RayEndpoint, SamplingDomainError, SheafAlgebraError,
    SheafBuildError, SheafBuildProgressUnit, SheafComplex, SheafVerdict, TripleCell, Vec3,
    validate_outside_ray_samples,
};
use fs_ivl::Interval;
use std::sync::atomic::{AtomicUsize, Ordering};

fn verdict_line(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-geom/sheaf\",\"case\":\"{case}\",\"verdict\":\"pass\",\
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

/// A rigid-motion wrapper: evaluates the inner chart in rotated+shifted
/// coordinates (a rotation about z by a fixed angle plus a translation).
struct Moved<C: Chart> {
    inner: C,
    cos: f64,
    sin: f64,
    shift: [f64; 3],
}

impl<C: Chart> Moved<C> {
    fn new(inner: C, angle: f64, shift: [f64; 3]) -> Self {
        Moved {
            inner,
            cos: angle.cos(),
            sin: angle.sin(),
            shift,
        }
    }

    fn map(&self, p: Point3) -> Point3 {
        // Inverse motion: un-shift, un-rotate.
        let q = [
            p.x - self.shift[0],
            p.y - self.shift[1],
            p.z - self.shift[2],
        ];
        Point3::new(
            self.cos * q[0] + self.sin * q[1],
            -self.sin * q[0] + self.cos * q[1],
            q[2],
        )
    }
}

impl<C: Chart> Chart for Moved<C> {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> fs_geom::ChartSample {
        self.inner.eval(self.map(x), cx)
    }

    fn support(&self) -> fs_geom::Aabb {
        // Conservative: rotate the support's corners, box them.
        let s = self.inner.support();
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for &x in &[s.min.x, s.max.x] {
            for &y in &[s.min.y, s.max.y] {
                for &z in &[s.min.z, s.max.z] {
                    // Forward motion of the corner.
                    let fx = self.cos * x - self.sin * y + self.shift[0];
                    let fy = self.sin * x + self.cos * y + self.shift[1];
                    let fz = z + self.shift[2];
                    min[0] = min[0].min(fx);
                    min[1] = min[1].min(fy);
                    min[2] = min[2].min(fz);
                    max[0] = max[0].max(fx);
                    max[1] = max[1].max(fy);
                    max[2] = max[2].max(fz);
                }
            }
        }
        fs_geom::Aabb::new(
            Point3::new(min[0], min[1], min[2]),
            Point3::new(max[0], max[1], max[2]),
        )
    }

    fn name(&self) -> &'static str {
        "test/moved"
    }
}

/// Two identical spheres offset so their supports overlap in a lens
/// around the shared surface — the WATERTIGHT fixture (same abstract
/// region presented twice).
fn watertight_pair() -> (SphereChart, SphereChart) {
    let s = SphereChart {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    (s, s)
}

/// The LEAKY variant: the second chart's radius is off by delta.
fn leaky_pair(delta: f64) -> (SphereChart, SphereChart) {
    let (a, mut b) = watertight_pair();
    b.radius += delta;
    (a, b)
}

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
fn sh_000_unbounded_interfaces_require_a_preflighted_clip() {
    with_cx(|cx| {
        let plane = UnboundedPlane;
        let charts: Vec<&dyn Chart> = vec![&plane, &plane];
        assert!(matches!(
            SheafComplex::from_charts(&charts, cx),
            Err(SheafBuildError::SamplingDomain {
                source: SamplingDomainError::UnboundedSupport { .. },
                ..
            })
        ));

        let clip = Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0));
        let local = SheafComplex::from_charts_clipped(&charts, clip, cx)
            .expect("explicit finite clip admits the unbounded interface");
        assert_eq!(local.n_patches, 2);
        assert_eq!(local.sampling_clip, Some(clip));

        let degenerate = Aabb::new(Point3::new(0.0, -1.0, -1.0), Point3::new(0.0, 1.0, 1.0));
        assert!(matches!(
            SheafComplex::from_charts_clipped(&[], degenerate, cx),
            Err(SheafBuildError::SamplingClip {
                source: SamplingDomainError::DegenerateDomain { .. }
            })
        ));
    });
}

struct EstimatedSphere(SphereChart);

impl Chart for EstimatedSphere {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> fs_geom::ChartSample {
        let mut sample = self.0.eval(x, cx);
        sample.error = fs_evidence::NumericalCertificate::estimate(
            sample.signed_distance - 1e-12,
            sample.signed_distance + 1e-12,
        );
        sample
    }

    fn support(&self) -> Aabb {
        self.0.support()
    }

    fn name(&self) -> &'static str {
        "test/estimated-sphere"
    }
}

struct AlternatingAuthoritySphere {
    inner: SphereChart,
    evaluations: AtomicUsize,
}

impl Chart for AlternatingAuthoritySphere {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> fs_geom::ChartSample {
        let mut sample = self.inner.eval(x, cx);
        if self.evaluations.fetch_add(1, Ordering::Relaxed) % 2 == 0 {
            sample.error = fs_evidence::NumericalCertificate::no_claim();
        }
        sample
    }

    fn support(&self) -> Aabb {
        self.inner.support()
    }

    fn name(&self) -> &'static str {
        "test/alternating-authority-sphere"
    }
}

#[test]
fn sh_000d_one_proven_sample_leak_survives_same_interface_indeterminacy() {
    with_cx(|cx| {
        let exact = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let mixed = AlternatingAuthoritySphere {
            inner: SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.01,
            },
            evaluations: AtomicUsize::new(0),
        };
        let charts: Vec<&dyn Chart> = vec![&exact, &mixed];
        let complex = SheafComplex::from_charts(&charts, cx).expect("bounded mixed interface");
        assert!(complex.interfaces[0].samples.len() > 1);
        let evidence = complex.watertightness(1e-3);
        match evidence.value {
            SheafVerdict::Fail {
                interface_violations,
                gauge_fit_share,
            } => {
                assert_eq!(interface_violations.len(), 1);
                assert!(interface_violations[0].1 > 1e-3);
                assert_eq!(
                    gauge_fit_share, None,
                    "an unrelated indeterminate sample disables the aggregate gauge fit"
                );
            }
            other => {
                panic!("one rigorous above-tolerance sample proves a localized leak: {other:?}")
            }
        }
        assert_eq!(evidence.numerical.kind, NumericalKind::NoClaim);
    });
}

#[test]
fn sh_000d_estimates_cannot_certify_or_falsify_watertightness() {
    with_cx(|cx| {
        let estimated = EstimatedSphere(SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        });
        let charts: Vec<&dyn Chart> = vec![&estimated, &estimated];
        let complex = SheafComplex::from_charts(&charts, cx).expect("bounded sheaf domain");
        assert!(!complex.interfaces[0].samples.is_empty());
        let evidence = complex.watertightness(1e-9);
        assert!(matches!(&evidence.value, SheafVerdict::Unknown { .. }));
        assert_eq!(
            evidence.numerical.kind,
            NumericalKind::NoClaim,
            "indeterminate chart authority must not become an enclosure receipt"
        );
        match &evidence.value {
            SheafVerdict::Unknown { reported_bounds } => assert!(
                reported_bounds.iter().any(|(_, _, hi)| hi.is_infinite()),
                "the unknown interface must retain its unbounded mismatch report"
            ),
            _ => unreachable!(),
        }
    });
}

struct PartlyNonFiniteSphere {
    inner: SphereChart,
    evaluations: AtomicUsize,
}

impl Chart for PartlyNonFiniteSphere {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> fs_geom::ChartSample {
        let mut sample = self.inner.eval(x, cx);
        if self.evaluations.fetch_add(1, Ordering::Relaxed) == 1 {
            sample.signed_distance = f64::NAN;
            sample.error = fs_evidence::NumericalCertificate::no_claim();
        }
        sample
    }

    fn support(&self) -> Aabb {
        self.inner.support()
    }

    fn name(&self) -> &'static str {
        "test/partly-nonfinite-sphere"
    }
}

#[test]
fn sh_000e_partly_nonfinite_producer_cannot_be_skipped_into_a_pass() {
    with_cx(|cx| {
        let chart = PartlyNonFiniteSphere {
            inner: SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            },
            evaluations: AtomicUsize::new(0),
        };
        let charts: Vec<&dyn Chart> = vec![&chart, &chart];
        assert!(matches!(
            SheafComplex::from_charts(&charts, cx),
            Err(SheafBuildError::NonFiniteSample {
                patches: (0, 1),
                chart: 1,
                value_bits,
                completed_draws: 1,
                ..
            }) if f64::from_bits(value_bits).is_nan()
        ));
    });
}

#[test]
fn sh_000b_clip_scope_is_retained_and_bound_into_provenance() {
    with_cx(|cx| {
        let (a, b) = watertight_pair();
        let charts: Vec<&dyn Chart> = vec![&a, &b];
        let clip = a.support();
        let wider_clip = clip.inflate(1.0);
        let global = SheafComplex::from_charts(&charts, cx).expect("bounded global support");
        let local = SheafComplex::from_charts_clipped(&charts, clip, cx)
            .expect("support-sized clip is admissible");
        let wider_local = SheafComplex::from_charts_clipped(&charts, wider_clip, cx)
            .expect("wider clip is admissible");

        assert_eq!(global.sampling_clip, None);
        assert_eq!(local.sampling_clip, Some(clip));
        assert_eq!(wider_local.sampling_clip, Some(wider_clip));
        assert_eq!(
            global.interfaces[0].samples, local.interfaces[0].samples,
            "an equal effective domain should isolate provenance scope binding"
        );
        assert_eq!(
            global.interfaces[0].samples, wider_local.interfaces[0].samples,
            "an enclosing clip should leave the effective sampled domain unchanged"
        );

        let global_evidence = global.watertightness(1e-9);
        let local_evidence = local.watertightness(1e-9);
        let wider_evidence = wider_local.watertightness(1e-9);
        assert_ne!(
            global_evidence.provenance, local_evidence.provenance,
            "local evidence must not share global provenance"
        );
        assert_ne!(
            local_evidence.provenance, wider_evidence.provenance,
            "exact clip endpoint bits must participate in provenance"
        );
    });
}

struct CancellingInterfaceChart<'a> {
    gate: &'a CancelGate,
    evaluations: AtomicUsize,
}

impl Chart for CancellingInterfaceChart<'_> {
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
        "test/cancelling-interface"
    }
}

#[test]
fn sh_000c_interface_sampler_polls_cancellation_directly() {
    let gate = CancelGate::new();
    let chart = CancellingInterfaceChart {
        gate: &gate,
        evaluations: AtomicUsize::new(0),
    };
    with_gate_cx(&gate, |cx| {
        let charts: Vec<&dyn Chart> = vec![&chart, &chart];
        let refusal = SheafComplex::from_charts(&charts, cx);
        assert!(
            matches!(
                refusal,
                Err(SheafBuildError::Cancelled {
                    stage: "interface-sampling",
                    patches: Some((0, 1)),
                    completed_work: 0,
                    unit: SheafBuildProgressUnit::InterfaceDraws,
                })
            ),
            "the sheaf sampler itself must observe cancellation from a non-polling chart: \
             {refusal:?}"
        );
    });
}

#[test]
fn sh_001_verdicts_and_localization() {
    with_cx(|cx| {
        // Watertight: identical spheres agree exactly on the seam.
        let (a, b) = watertight_pair();
        let charts: Vec<&dyn Chart> = vec![&a, &b];
        let complex = SheafComplex::from_charts(&charts, cx).expect("bounded sheaf domain");
        assert_eq!(complex.interfaces.len(), 1, "one shared interface");
        assert!(!complex.interfaces[0].samples.is_empty());
        let ev = complex.watertightness(1e-9);
        let bounds = complex.mismatch_bounds().expect("bounded mismatch output");
        assert!(!bounds.is_empty());
        assert!(bounds.iter().all(|bound| bound.determinate()));
        assert!(bounds.iter().all(|bound| bound.all_within(1e-9)));
        assert!(bounds.iter().all(|bound| !bound.proven_leak(1e-9)));
        match &ev.value {
            SheafVerdict::Pass {
                worst_mismatch,
                margins,
            } => {
                assert!(*worst_mismatch <= 1e-9);
                assert_eq!(margins.len(), 1);
            }
            other => panic!("identical charts must certify: {other:?}"),
        }
        assert_eq!(
            ev.numerical.kind,
            NumericalKind::Enclosure,
            "interval-verified"
        );
        let raw_view: &SheafComplex = &complex;
        let downgraded = raw_view.watertightness(1e-9);
        assert!(matches!(downgraded.value, SheafVerdict::Unknown { .. }));
        assert_eq!(downgraded.numerical.kind, NumericalKind::NoClaim);
        assert_ne!(
            ev.provenance, downgraded.provenance,
            "admitted assessment and an explicit raw-view downgrade must bind distinct origins"
        );
        // Leaky: radius off by 1e-2 — FAIL, localized to the (0,1) seam,
        // with the mismatch magnitude ~ delta. This two-vertex, one-edge
        // adjacency graph is a tree, so its first graph cohomology is zero:
        // the failure is an interface violation, not an H¹ obstruction.
        let (a2, b2) = leaky_pair(1e-2);
        let charts2: Vec<&dyn Chart> = vec![&a2, &b2];
        let complex2 = SheafComplex::from_charts(&charts2, cx).expect("bounded sheaf domain");
        let ev2 = complex2.watertightness(1e-4);
        let leaky_bounds = complex2
            .mismatch_bounds()
            .expect("bounded leaky mismatch output");
        assert!(
            leaky_bounds.iter().any(|bound| bound.proven_leak(1e-4)),
            "strict tolerance derives a leak from the retained numeric enclosure"
        );
        assert!(
            leaky_bounds.iter().all(|bound| bound.all_within(1.0)),
            "the same context-free bounds can be re-evaluated under a loose tolerance without a detached stale boolean"
        );
        match &ev2.value {
            SheafVerdict::Fail {
                interface_violations,
                gauge_fit_share,
            } => {
                assert_eq!(
                    complex2.n_patches,
                    complex2.interfaces.len() + 1,
                    "the connected adjacency graph is a tree and has H1 = 0"
                );
                assert_eq!(
                    interface_violations[0].0,
                    (0, 1),
                    "offending interface named"
                );
                assert!(
                    interface_violations[0].1 > 5e-3 && interface_violations[0].1 < 2e-2,
                    "mismatch magnitude ~ delta: {}",
                    interface_violations[0].1
                );
                // A pure radius offset is a CONSTANT mismatch on the seam —
                // graph-gauge explained, without a topology claim.
                assert!(
                    gauge_fit_share.is_some_and(|share| share > 0.9),
                    "constant seam mismatch is graph-gauge explained: \
                     {gauge_fit_share:?}"
                );
            }
            other => panic!("leaky seam must fail: {other:?}"),
        }
        verdict_line(
            "sh-001",
            "identical charts PASS; radius leak FAILs at the named seam with ~delta \
             magnitude and >0.9 graph-gauge explained share, without an H1 claim",
        );
    });
}

#[test]
fn sh_001b_overflowing_section_diagnostic_is_explicitly_unavailable() {
    let sample = |mismatch: f64| InterfaceSample {
        point: Point3::new(0.0, 0.0, 0.0),
        values: [Interval::point(0.0), Interval::point(mismatch)],
    };
    let complex = SheafComplex {
        sampling_clip: None,
        n_patches: 3,
        interfaces: vec![
            Interface {
                patches: (0, 1),
                samples: vec![sample(1e200)],
            },
            Interface {
                patches: (0, 2),
                samples: vec![sample(1e200)],
            },
            Interface {
                patches: (1, 2),
                samples: vec![sample(1e200)],
            },
        ],
        triples: Vec::new(),
    };
    assert_eq!(
        complex.section_solve(),
        Err(SheafAlgebraError::NumericalOverflow {
            stage: "section-mean-square",
        }),
        "an unrepresentable diagnostic must refuse rather than return non-finite success"
    );
    let evidence = complex.watertightness(1.0);
    assert!(matches!(evidence.value, SheafVerdict::Unknown { .. }));
    assert_eq!(evidence.numerical.kind, NumericalKind::NoClaim);
}

#[test]
fn sh_001c_three_patch_pure_gauge_mismatch_is_fully_explained() {
    let sample = |left: f64, right: f64| InterfaceSample {
        point: Point3::new(0.0, 0.0, 0.0),
        values: [Interval::point(left), Interval::point(right)],
    };
    // Patch potentials [0, 1, 3] induce edge mismatches [1, 2, 3]. The
    // triangle has a 2-cell and the loop sum is exactly zero, so this is a
    // pure graph gauge rather than a structural/topological witness.
    let complex = SheafComplex {
        sampling_clip: None,
        n_patches: 3,
        interfaces: vec![
            Interface {
                patches: (0, 1),
                samples: vec![sample(0.0, 1.0)],
            },
            Interface {
                patches: (0, 2),
                samples: vec![sample(0.0, 3.0)],
            },
            Interface {
                patches: (1, 2),
                samples: vec![sample(1.0, 3.0)],
            },
        ],
        triples: vec![TripleCell {
            patches: (0, 1, 2),
            samples: 1,
        }],
    };

    let potentials = [0.0, 1.0, 3.0];
    let d0 = complex.delta0_edges().expect("valid edge incidence");
    let mismatch: Vec<f64> = (0..d0.nrows())
        .map(|row| {
            let (columns, values) = d0.row(row);
            columns
                .iter()
                .zip(values)
                .map(|(&column, &value)| value * potentials[column])
                .sum()
        })
        .collect();
    assert_eq!(mismatch, vec![1.0, 3.0, 2.0]);
    let d1 = complex.delta1().expect("valid triple incidence");
    let (columns, values) = d1.row(0);
    let loop_sum: f64 = columns
        .iter()
        .zip(values)
        .map(|(&column, &value)| value * mismatch[column])
        .sum();
    assert_eq!(
        loop_sum.to_bits(),
        0.0f64.to_bits(),
        "the retained pure-gauge cochain must be closed bitwise"
    );

    let (_, raw, residual) = complex.section_solve().expect("valid diagnostic complex");
    assert!(raw > 0.0);
    assert!(
        residual <= raw * 1e-24,
        "pure gauge mismatch must have no edge-mean fit residual: raw={raw}, residual={residual}"
    );
    let raw_evidence = complex.watertightness(1e-6);
    assert!(matches!(raw_evidence.value, SheafVerdict::Unknown { .. }));
    assert_eq!(raw_evidence.numerical.kind, NumericalKind::NoClaim);
}

#[test]
fn sh_001c2_gauge_share_measures_sample_energy_not_only_edge_means() {
    let complex = SheafComplex {
        sampling_clip: None,
        n_patches: 2,
        interfaces: vec![Interface {
            patches: (0, 1),
            samples: vec![
                InterfaceSample {
                    point: Point3::new(0.0, 0.0, 0.0),
                    values: [Interval::point(0.0), Interval::point(100.0)],
                },
                InterfaceSample {
                    point: Point3::new(1.0, 0.0, 0.0),
                    values: [Interval::point(0.0), Interval::point(-98.0)],
                },
            ],
        }],
        triples: Vec::new(),
    };
    let (_, raw, residual) = complex.section_solve().expect("valid diagnostic complex");
    assert!((raw - 9_802.0).abs() <= f64::EPSILON * raw);
    assert!((residual - 9_801.0).abs() <= f64::EPSILON * residual);
    let raw_evidence = complex.watertightness(1.0);
    assert!(matches!(raw_evidence.value, SheafVerdict::Unknown { .. }));
    assert_eq!(raw_evidence.numerical.kind, NumericalKind::NoClaim);
}

#[test]
fn sh_001d_raw_interval_parts_cannot_mint_verdict_authority() {
    let complex_with_mismatch = |mismatch: Interval| SheafComplex {
        sampling_clip: None,
        n_patches: 2,
        interfaces: vec![Interface {
            patches: (0, 1),
            samples: vec![InterfaceSample {
                point: Point3::new(0.0, 0.0, 0.0),
                values: [Interval::point(0.0), mismatch],
            }],
        }],
        triples: Vec::new(),
    };

    let pass = complex_with_mismatch(Interval::new(0.25, 0.75)).watertightness(1.0);
    assert!(matches!(pass.value, SheafVerdict::Unknown { .. }));
    assert_eq!(pass.numerical.kind, NumericalKind::NoClaim);

    let fail = complex_with_mismatch(Interval::new(1.25, 1.5)).watertightness(1.0);
    assert!(matches!(fail.value, SheafVerdict::Unknown { .. }));
    assert_eq!(fail.numerical.kind, NumericalKind::NoClaim);

    let unknown = complex_with_mismatch(Interval::new(0.75, 1.25)).watertightness(1.0);
    assert!(matches!(unknown.value, SheafVerdict::Unknown { .. }));
    assert_eq!(unknown.numerical.kind, NumericalKind::NoClaim);
}

#[test]
fn sh_001e_malformed_or_oversized_raw_algebra_refuses_without_panicking() {
    let malformed = SheafComplex {
        sampling_clip: None,
        n_patches: 2,
        interfaces: vec![Interface {
            patches: (0, 2),
            samples: vec![InterfaceSample {
                point: Point3::new(0.0, 0.0, 0.0),
                values: [Interval::point(0.0), Interval::point(1.0)],
            }],
        }],
        triples: Vec::new(),
    };
    assert_eq!(malformed.delta0(), Err(SheafAlgebraError::InvalidStructure));
    assert_eq!(
        malformed.delta0_edges(),
        Err(SheafAlgebraError::InvalidStructure)
    );
    assert_eq!(malformed.delta1(), Err(SheafAlgebraError::InvalidStructure));
    assert_eq!(
        malformed.section_solve(),
        Err(SheafAlgebraError::InvalidStructure)
    );

    let oversized = SheafComplex {
        sampling_clip: None,
        n_patches: usize::MAX,
        interfaces: Vec::new(),
        triples: Vec::new(),
    };
    assert!(matches!(
        oversized.section_solve(),
        Err(SheafAlgebraError::WorkLimit {
            stage: "patches",
            ..
        })
    ));

    let indeterminate = SheafComplex {
        sampling_clip: None,
        n_patches: 2,
        interfaces: vec![Interface {
            patches: (0, 1),
            samples: vec![InterfaceSample {
                point: Point3::new(0.0, 0.0, 0.0),
                values: [Interval::point(0.0), Interval::WHOLE],
            }],
        }],
        triples: Vec::new(),
    };
    assert_eq!(
        indeterminate.section_solve(),
        Err(SheafAlgebraError::IndeterminateSampleValue {
            interface: 0,
            sample: 0,
        })
    );

    let overflowing = SheafComplex {
        sampling_clip: None,
        n_patches: 2,
        interfaces: vec![Interface {
            patches: (0, 1),
            samples: vec![InterfaceSample {
                point: Point3::new(0.0, 0.0, 0.0),
                values: [Interval::point(-f64::MAX), Interval::point(f64::MAX)],
            }],
        }],
        triples: Vec::new(),
    };
    assert_eq!(
        overflowing.section_solve(),
        Err(SheafAlgebraError::NumericalOverflow {
            stage: "section-edge-mismatch",
        })
    );
}

#[test]
fn sh_001e_public_structure_cannot_vacuously_pass_or_alias_provenance() {
    let empty_interface = SheafComplex {
        sampling_clip: None,
        n_patches: 2,
        interfaces: vec![Interface {
            patches: (0, 1),
            samples: Vec::new(),
        }],
        triples: Vec::new(),
    };
    let empty = empty_interface.watertightness(1.0);
    assert!(matches!(empty.value, SheafVerdict::Unknown { .. }));
    assert_eq!(empty.numerical.kind, NumericalKind::NoClaim);
    assert!(
        empty.qoi.is_infinite(),
        "an empty sampled interface must not expose a false-clean zero QoI"
    );

    let malformed_scope = SheafComplex {
        sampling_clip: Some(Aabb {
            min: Point3::new(1.0, 0.0, 0.0),
            max: Point3::new(0.0, 1.0, 1.0),
        }),
        n_patches: 2,
        interfaces: vec![Interface {
            patches: (0, 1),
            samples: vec![InterfaceSample {
                point: Point3::new(0.5, 0.5, 0.5),
                values: [Interval::point(0.0), Interval::point(2.0)],
            }],
        }],
        triples: Vec::new(),
    }
    .watertightness(1.0);
    assert!(matches!(
        malformed_scope.value,
        SheafVerdict::Unknown { .. }
    ));
    assert_eq!(malformed_scope.numerical.kind, NumericalKind::NoClaim);

    let out_of_scope = SheafComplex {
        sampling_clip: Some(Aabb::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 1.0),
        )),
        n_patches: 2,
        interfaces: vec![Interface {
            patches: (0, 1),
            samples: vec![InterfaceSample {
                point: Point3::new(2.0, 0.5, 0.5),
                values: [Interval::point(0.0), Interval::point(2.0)],
            }],
        }],
        triples: Vec::new(),
    }
    .watertightness(1.0);
    assert!(matches!(out_of_scope.value, SheafVerdict::Unknown { .. }));
    assert_eq!(out_of_scope.numerical.kind, NumericalKind::NoClaim);

    let sample = |point: Point3, value: Interval| InterfaceSample {
        point,
        values: [Interval::point(0.0), value],
    };
    let mixed = SheafComplex {
        sampling_clip: None,
        n_patches: 3,
        interfaces: vec![
            Interface {
                patches: (0, 1),
                samples: vec![sample(Point3::new(0.0, 0.0, 0.0), Interval::new(2.0, 3.0))],
            },
            Interface {
                patches: (1, 2),
                samples: vec![sample(Point3::new(1.0, 0.0, 0.0), Interval::WHOLE)],
            },
        ],
        triples: Vec::new(),
    };
    let mixed_evidence = mixed.watertightness(1.0);
    assert!(matches!(mixed_evidence.value, SheafVerdict::Unknown { .. }));
    assert_eq!(mixed_evidence.numerical.kind, NumericalKind::NoClaim);

    let evidence_at = |x: f64| {
        SheafComplex {
            sampling_clip: None,
            n_patches: 2,
            interfaces: vec![Interface {
                patches: (0, 1),
                samples: vec![sample(Point3::new(x, 0.0, 0.0), Interval::new(0.25, 0.5))],
            }],
            triples: Vec::new(),
        }
        .watertightness(1.0)
    };
    assert_ne!(
        evidence_at(0.0).provenance,
        evidence_at(1.0).provenance,
        "sample coordinates remain provenance-bound diagnostic inputs"
    );
}

#[test]
fn sh_002_delta_delta_is_zero_bitwise() {
    with_cx(|cx| {
        // Three pairwise-overlapping boxes: the current discovery path retains
        // a candidate triangle from the pairwise-interface clique. It is
        // sufficient for the algebraic δδ fixture, but it is not verified
        // common-intersection or aligned-restriction evidence.
        let boxes: Vec<BoxChart> = (0..3)
            .map(|i| {
                let base = f64::from(i) * 0.4;
                BoxChart {
                    aabb: fs_geom::Aabb::new(
                        Point3::new(base - 0.5, -0.5, -0.5),
                        Point3::new(base + 0.5, 0.5, 0.5),
                    ),
                }
            })
            .collect();
        let charts: Vec<&dyn Chart> = boxes.iter().map(|b| b as &dyn Chart).collect();
        let complex = SheafComplex::from_charts(&charts, cx).expect("bounded sheaf domain");
        assert!(!complex.interfaces.is_empty());
        assert!(
            !complex.triples.is_empty(),
            "fixture must produce a candidate triangle cell; adjust geometry"
        );
        let d0 = complex.delta0_edges().expect("valid edge incidence");
        let d1 = complex.delta1().expect("valid triple incidence");
        let sampled_d0 = complex.delta0().expect("valid sampled incidence");
        // δ¹ · δ⁰ computed densely (test scale): every entry EXACTLY 0.0.
        let (rows, mid, cols) = (
            complex.triples.len(),
            complex.interfaces.len(),
            complex.n_patches,
        );
        assert_eq!(d0.nrows(), mid);
        assert_eq!(d0.ncols(), cols);
        assert_eq!(d1.nrows(), rows);
        assert_eq!(d1.ncols(), mid);
        assert_eq!(d1.ncols(), d0.nrows());
        assert_eq!(
            sampled_d0.nrows(),
            complex
                .interfaces
                .iter()
                .map(|interface| interface.samples.len())
                .sum::<usize>(),
            "sample-row restriction incidence is distinct from edge-level delta0"
        );
        assert_eq!(sampled_d0.ncols(), cols);
        for r in 0..rows {
            for c in 0..cols {
                let mut acc = 0.0f64;
                let (d1_cols, d1_vals) = d1.row(r);
                for (k, &e) in d1_cols.iter().enumerate() {
                    let (d0_cols, d0_vals) = d0.row(e);
                    for (j, &p) in d0_cols.iter().enumerate() {
                        if p == c {
                            acc += d1_vals[k] * d0_vals[j];
                        }
                    }
                }
                assert_eq!(acc.to_bits(), 0.0f64.to_bits(), "δδ must be bitwise zero");
            }
        }
        verdict_line(
            "sh-002",
            "δ¹·δ⁰ == 0 bitwise on a retained pairwise-clique 2-cell; no common-overlap claim",
        );
    });
}

#[test]
fn sh_003_pair_bound_swap_exact_rigid_tolerance() {
    with_cx(|cx| {
        let (a, b) = leaky_pair(2e-2);
        // Swapping this two-chart pair preserves the localized mismatch bound
        // exactly because geometry-derived sample seeds are index-free. Full
        // evidence/provenance still binds numeric patch labels.
        let fwd: Vec<&dyn Chart> = vec![&a, &b];
        let rev: Vec<&dyn Chart> = vec![&b, &a];
        let v1 = SheafComplex::from_charts(&fwd, cx)
            .expect("bounded sheaf domain")
            .watertightness(1e-4);
        let v2 = SheafComplex::from_charts(&rev, cx)
            .expect("bounded sheaf domain")
            .watertightness(1e-4);
        let key = |v: &SheafVerdict| match v {
            SheafVerdict::Fail {
                interface_violations,
                ..
            } => interface_violations[0].1,
            _ => f64::NAN,
        };
        assert_eq!(
            key(&v1.value).to_bits(),
            key(&v2.value).to_bits(),
            "the swapped two-chart mismatch bound is exact"
        );
        // Rigid motion: rotate+translate BOTH charts — verdict class and
        // magnitude agree to tolerance (samples differ, physics doesn't).
        let ma = Moved::new(a, 0.7, [3.0, -1.0, 0.5]);
        let mb = Moved::new(b, 0.7, [3.0, -1.0, 0.5]);
        let moved: Vec<&dyn Chart> = vec![&ma, &mb];
        let v3 = SheafComplex::from_charts(&moved, cx)
            .expect("bounded sheaf domain")
            .watertightness(1e-4);
        match (&v1.value, &v3.value) {
            (
                SheafVerdict::Fail {
                    interface_violations: o1,
                    ..
                },
                SheafVerdict::Fail {
                    interface_violations: o3,
                    ..
                },
            ) => {
                assert!(
                    (o1[0].1 - o3[0].1).abs() < 5e-3,
                    "rigid motion preserves the leak magnitude: {} vs {}",
                    o1[0].1,
                    o3[0].1
                );
            }
            other => panic!("verdict class must survive rigid motion: {other:?}"),
        }
        verdict_line(
            "sh-003",
            "two-chart swapped-pair bound bitwise; rigid-motion invariance to 5e-3",
        );
    });
}

#[test]
fn sh_004_adversarial_seams_and_soundness() {
    with_cx(|cx| {
        // T-junction stress case: three distinct boxes with nearby faces.
        let mk = |cx: f64, cy: f64| BoxChart {
            aabb: fs_geom::Aabb::new(
                Point3::new(cx - 0.5, cy - 0.5, -0.5),
                Point3::new(cx + 0.5, cy + 0.5, 0.5),
            ),
        };
        let b1 = mk(0.0, 0.0);
        let b2 = mk(0.9, 0.0);
        let b3 = mk(0.45, 0.9);
        let charts: Vec<&dyn Chart> = vec![&b1, &b2, &b3];
        let complex = SheafComplex::from_charts(&charts, cx).expect("bounded sheaf domain");
        // Adjacent distinct box SDFs disagree off the shared face, so this is
        // genuinely adversarial: the exercised tight-tolerance verdict must
        // not be a false PASS.
        let ev = complex.watertightness(1e-12);
        assert!(
            !matches!(ev.value, SheafVerdict::Pass { .. }),
            "distinct-box seams must not certify at 1e-12"
        );
        // Near-separated spheres (distinct surfaces with a small gap): NOT a
        // false PASS.
        let s1 = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let s2 = SphereChart {
            center: Point3::new(2.05, 0.0, 0.0),
            radius: 1.0,
        };
        let separated: Vec<&dyn Chart> = vec![&s1, &s2];
        let separated_complex =
            SheafComplex::from_charts(&separated, cx).expect("disjoint supports are admissible");
        // The 0.05 gap means NO overlap interface is discovered, so the
        // interface-agreement check has gathered zero evidence. It must NOT
        // report a positive PASS from an empty bound set (bead obnw: `all()`
        // on the empty set was vacuously true) — the honest verdict is Unknown.
        let kv = separated_complex.watertightness(1e-9);
        assert!(
            matches!(kv.value, SheafVerdict::Unknown { .. }),
            "near-separated distinct surfaces (no interface) must be Unknown, not a \
             false PASS, got {:?}",
            kv.value
        );
        // Sign-sequence replay: this validates the bounded sampler and
        // outside-endpoint contract only. Outside/outside boolean toggles are
        // necessarily even, so this is not independent topology evidence.
        let (wa, wb) = watertight_pair();
        let watertight: Vec<&dyn Chart> = vec![&wa, &wb];
        let rays = [
            (Point3::new(-3.0, 0.01, 0.02), Point3::new(3.0, 0.01, 0.02)),
            (Point3::new(0.02, -3.0, 0.01), Point3::new(0.02, 3.0, 0.01)),
            (Point3::new(-2.5, -2.5, 0.0), Point3::new(2.5, 2.5, 0.0)),
        ];
        let ray_report = validate_outside_ray_samples(&watertight, &rays, 4001, cx)
            .expect("finite outside rays fit the public work cap");
        assert_eq!(ray_report.rays, rays.len());
        assert_eq!(ray_report.sample_points, rays.len() * 4002);
        assert_eq!(ray_report.chart_evaluations, rays.len() * 4002 * 2);
        assert_eq!(
            ray_report.toggles, 6,
            "each sphere-crossing ray must retain two nonzero transitions"
        );
        verdict_line(
            "sh-004",
            "these distinct-surface fixtures do not falsely certify; the legacy sign scan is \
             retained only as a bounded replay diagnostic",
        );
    });
}

struct ConstantRayChart {
    value: f64,
}

struct CertifiedConstantRayChart {
    value: f64,
    error: fs_evidence::NumericalCertificate,
    support: Aabb,
}

struct CountingOverlapChart<'a> {
    evaluations: &'a AtomicUsize,
}

impl Chart for CountingOverlapChart<'_> {
    fn eval(&self, _x: Point3, _cx: &Cx<'_>) -> fs_geom::ChartSample {
        self.evaluations.fetch_add(1, Ordering::SeqCst);
        fs_geom::ChartSample {
            signed_distance: 1.0,
            gradient: None,
            lipschitz: None,
            error: fs_evidence::NumericalCertificate::no_claim(),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "test/counting-overlap"
    }
}

impl Chart for CertifiedConstantRayChart {
    fn eval(&self, _x: Point3, _cx: &Cx<'_>) -> fs_geom::ChartSample {
        fs_geom::ChartSample {
            signed_distance: self.value,
            gradient: None,
            lipschitz: None,
            error: self.error,
        }
    }

    fn support(&self) -> Aabb {
        self.support
    }

    fn name(&self) -> &'static str {
        "test/certified-constant-ray"
    }
}

impl Chart for ConstantRayChart {
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
        "test/constant-ray"
    }
}

#[test]
fn sh_004b_outside_ray_sampling_refuses_invalid_work_and_is_stable() {
    with_cx(|cx| {
        let outside = ConstantRayChart { value: 1.0 };
        let charts: Vec<&dyn Chart> = vec![&outside];
        let ray = [(Point3::new(-2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0))];

        assert_eq!(
            validate_outside_ray_samples(&[], &ray, 2, cx),
            Err(OutsideRaySampleError::EmptyCharts)
        );
        assert_eq!(
            validate_outside_ray_samples(&charts, &[], 2, cx),
            Err(OutsideRaySampleError::EmptyRays)
        );
        assert_eq!(
            validate_outside_ray_samples(&charts, &ray, 0, cx),
            Err(OutsideRaySampleError::InvalidSteps { steps: 0 })
        );
        assert!(matches!(
            validate_outside_ray_samples(&charts, &ray, OUTSIDE_RAY_MAX_EVALUATIONS, cx),
            Err(OutsideRaySampleError::WorkLimitExceeded { .. })
        ));
        let inside = ConstantRayChart { value: -1.0 };
        let inside_charts: Vec<&dyn Chart> = vec![&inside];
        assert!(matches!(
            validate_outside_ray_samples(&inside_charts, &ray, 2, cx),
            Err(OutsideRaySampleError::EndpointNotOutside {
                ray: 0,
                endpoint: RayEndpoint::Start,
                ..
            })
        ));
        let unproven_ray = [(Point3::new(-0.5, 0.0, 0.0), Point3::new(0.5, 0.0, 0.0))];
        assert!(matches!(
            validate_outside_ray_samples(&charts, &unproven_ray, 2, cx),
            Err(OutsideRaySampleError::EndpointOutsideUnproven {
                ray: 0,
                endpoint: RayEndpoint::Start,
                chart: 0,
                ..
            })
        ));
        let certified = CertifiedConstantRayChart {
            value: 1.0,
            error: fs_evidence::NumericalCertificate::enclosure(0.5, 1.5),
            support: Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0)),
        };
        let support_excluded = CertifiedConstantRayChart {
            value: 2.0,
            error: fs_evidence::NumericalCertificate::no_claim(),
            support: Aabb::new(Point3::new(-0.1, -0.1, -0.1), Point3::new(0.1, 0.1, 0.1)),
        };
        let mixed: Vec<&dyn Chart> = vec![&certified, &support_excluded];
        assert_eq!(
            validate_outside_ray_samples(&mixed, &unproven_ray, 2, cx),
            Ok(OutsideRaySampleReport {
                rays: 1,
                sample_points: 3,
                chart_evaluations: 6,
                toggles: 0,
            }),
            "one rigorous positive enclosure plus one support exclusion proves each endpoint"
        );
        let malformed = CertifiedConstantRayChart {
            value: 1.0,
            error: fs_evidence::NumericalCertificate {
                kind: NumericalKind::Exact,
                lo: 0.0,
                hi: 0.0,
            },
            support: certified.support,
        };
        assert!(matches!(
            validate_outside_ray_samples(&[&malformed as &dyn Chart], &unproven_ray, 2, cx,),
            Err(OutsideRaySampleError::EndpointOutsideUnproven { .. })
        ));
        assert!(matches!(
            validate_outside_ray_samples(
                &charts,
                &[(Point3::new(f64::NAN, 0.0, 0.0), ray[0].1)],
                2,
                cx,
            ),
            Err(OutsideRaySampleError::NonFiniteEndpoint {
                ray: 0,
                endpoint: RayEndpoint::Start,
                ..
            })
        ));

        let huge = [(
            Point3::new(-f64::MAX, -f64::MAX, -f64::MAX),
            Point3::new(f64::MAX, f64::MAX, f64::MAX),
        )];
        assert_eq!(
            validate_outside_ray_samples(&charts, &huge, 2, cx),
            Ok(OutsideRaySampleReport {
                rays: 1,
                sample_points: 3,
                chart_evaluations: 3,
                toggles: 0,
            })
        );

        let invalid = ConstantRayChart { value: f64::NAN };
        let invalid_charts: Vec<&dyn Chart> = vec![&invalid];
        assert!(matches!(
            validate_outside_ray_samples(&invalid_charts, &ray, 2, cx),
            Err(OutsideRaySampleError::NonFiniteSample {
                ray: 0,
                step: 0,
                chart: 0,
                ..
            })
        ));
    });
}

#[test]
fn sh_004d_dense_overlap_is_refused_before_chart_evaluation() {
    with_cx(|cx| {
        let evaluations = AtomicUsize::new(0);
        let chart = CountingOverlapChart {
            evaluations: &evaluations,
        };
        let charts = vec![&chart as &dyn Chart; 100];
        let refusal = SheafComplex::from_charts(&charts, cx);
        assert!(matches!(
            refusal,
            Err(SheafBuildError::BuildWorkLimit {
                stage: "interface-sampling-evaluations",
                requested: 16_781_312,
                cap: 16_777_216,
            })
        ));
        assert_eq!(
            evaluations.load(Ordering::SeqCst),
            0,
            "pair/evaluation admission must complete before any chart evaluation"
        );
    });
}

#[test]
fn sh_004c_outside_ray_sampling_observes_chart_cancellation() {
    let gate = CancelGate::new();
    let chart = CancellingInterfaceChart {
        gate: &gate,
        evaluations: AtomicUsize::new(0),
    };
    with_gate_cx(&gate, |cx| {
        let charts: Vec<&dyn Chart> = vec![&chart];
        let rays = [(Point3::new(2.0, 0.0, 0.0), Point3::new(3.0, 0.0, 0.0))];
        assert_eq!(
            validate_outside_ray_samples(&charts, &rays, 2, cx),
            Err(OutsideRaySampleError::Cancelled {
                completed_rays: 0,
                completed_points: 0,
                completed_chart_evaluations: 1,
            })
        );
    });
}

#[test]
fn sh_005_section_solve_reports_graph_gauge_fit_only() {
    with_cx(|cx| {
        // A chain of three spheres where the middle chart carries a
        // constant radial offset: per-patch graph offsets should fit nearly
        // all of the sampled edge-mean mismatch. This numerical diagnostic
        // does not itself classify a cochain or feed merge semantics.
        let a = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let b = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.015,
        };
        let charts: Vec<&dyn Chart> = vec![&a, &b];
        let complex = SheafComplex::from_charts(&charts, cx).expect("bounded sheaf domain");
        let (offsets, raw, residual) = complex
            .section_solve()
            .expect("valid disconnected diagnostic complex");
        assert!(raw > 0.0, "the leak is visible pre-gauge");
        assert!(
            residual < raw * 0.01,
            "a constant edge-mean offset is graph-gauge fitted: raw {raw}, residual {residual}"
        );
        assert!(
            (offsets[1] - offsets[0] - 0.015).abs() < 1e-3,
            "the recovered gauge cancels the radius delta: {offsets:?}"
        );
        verdict_line(
            "sh-005",
            "constant seam offsets are fitted by the graph-gauge section diagnostic without \
             a topology or merge-classification claim",
        );
    });
}
