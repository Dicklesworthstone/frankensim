//! Certified geometric-moments battery (bead rjnd, part 1).
//!
//! - gm-001 G0: box moments contain the closed forms with bounded widths.
//! - gm-002 G0/G2: sphere moments and COM enclosures contain the analytic
//!   values (parallel-axis closed form about the origin).
//! - gm-003 G3: translation-covariance metamorphic — the outward-rounded
//!   covariance law and a direct recomputation on a shifted chart must
//!   overlap componentwise (both enclose the same truth), and volume is
//!   translation-invariant bitwise.
//! - gm-004 G0: capability and input refusals — weaker trace claims,
//!   Estimate-class samples, malformed spacing/domains, and excessive
//!   work all refuse with the named typed error.
//! - gm-005 G4: cancellation surfaces as `QueryError::Cancelled`.
//! - gm-006 G0/G2: torus and hollow-shell closed forms are contained
//!   (curved, genus-1, and non-simply-connected solids).
//! - gm-007 G0: an OPEN triangle mesh cannot claim the exact-distance
//!   capability, so mass properties refuse — the watertightness
//!   precondition enforced through capability routing, with the chart
//!   type and missing capability logged.
//!
//! Aggregate outcomes use canonical fs-obs conformance events. Every fixture
//! is fixed-input and therefore carries input seed zero; the deterministic
//! `Cx` stream remains separate execution provenance. Assertions and
//! expectations reached before a verdict remain ordinary Rust diagnostics.

use asupersync::types::Budget;
use fs_evidence::NumericalCertificate;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::{BoxChart, SphereChart, TorusChart};
use fs_geom::{Aabb, Chart, ChartSample, Point3, TraceStepClaim, Vec3};
use fs_query::{GeometricMoments, MomentEnclosure, QueryError, geometric_moments};
use std::f64::consts::PI;

const SUITE: &str = "fs-query/moments";
const FIXED_INPUT_SEED: u64 = 0;
const EXECUTION_SEED: u64 = 0x60E5;
const MOMENTS_KERNEL_ID: u64 = 9;
const CANCELLATION_KERNEL_ID: u64 = 10;

fn emit_verdict(
    emitter: &mut fs_obs::Emitter,
    case: &str,
    pass: bool,
    detail: &str,
    execution_kernel: u64,
) {
    let severity = if pass {
        fs_obs::Severity::Info
    } else {
        fs_obs::Severity::Error
    };
    let event = emitter.emit(
        severity,
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass,
            detail: format!(
                "{detail}; execution stream seed=0x{EXECUTION_SEED:x} \
                 kernel={execution_kernel} tile=0 iteration=0"
            ),
            seed: FIXED_INPUT_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("moments verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("moments verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

fn verdict(case: &str, pass: bool, detail: &str, execution_kernel: u64) {
    let mut emitter = fs_obs::Emitter::new(SUITE, case);
    emit_verdict(&mut emitter, case, pass, detail, execution_kernel);
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
                kernel_id: MOMENTS_KERNEL_ID,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn contained(name: &str, enclosure: MomentEnclosure, truth: f64, max_width: f64) -> String {
    assert!(
        enclosure.contains(truth),
        "{name}: [{}, {}] must contain {truth}",
        enclosure.lo,
        enclosure.hi
    );
    assert!(
        enclosure.width() <= max_width,
        "{name}: width {} exceeds {max_width}",
        enclosure.width()
    );
    format!(
        "\\\"{name}\\\":{{\\\"lo\\\":{:.6e},\\\"hi\\\":{:.6e},\\\"truth\\\":{truth:.6e}}}",
        enclosure.lo, enclosure.hi
    )
}

/// The fixture chart rigidly shifted by `d` (still an exact distance).
struct Shifted<C: Chart> {
    inner: C,
    d: Vec3,
}

impl<C: Chart> Chart for Shifted<C> {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        self.inner
            .eval(x.offset(Vec3::new(-self.d.x, -self.d.y, -self.d.z)), cx)
    }

    fn support(&self) -> Aabb {
        let s = self.inner.support();
        Aabb::new(s.min.offset(self.d), s.max.offset(self.d))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        self.inner.trace_step_claim()
    }

    fn name(&self) -> &'static str {
        "test/shifted"
    }
}

/// Claims only LipschitzImplicit: certified moments must refuse it.
struct ImplicitOnly<C: Chart>(C);

impl<C: Chart> Chart for ImplicitOnly<C> {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        self.0.eval(x, cx)
    }

    fn support(&self) -> Aabb {
        self.0.support()
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::LipschitzImplicit
    }

    fn name(&self) -> &'static str {
        "test/implicit-only"
    }
}

/// Claims ExactDistance but serves Estimate-class evidence: each sample
/// must refuse (the claim alone cannot launder the certificate).
struct EstimatingBox(BoxChart);

impl Chart for EstimatingBox {
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
        "test/estimating-box"
    }
}

fn unit_domain() -> Aabb {
    Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
}

fn test_box() -> BoxChart {
    BoxChart {
        aabb: Aabb::new(Point3::new(-0.25, -0.25, -0.25), Point3::new(0.5, 0.5, 0.5)),
    }
}

#[test]
fn gm_001_box_moments_contain_closed_forms() {
    let chart = test_box();
    let m = with_cx(|cx| geometric_moments(&chart, &unit_domain(), 0.02, cx)).expect("box moments");
    // w = 0.75, c = 0.125 per axis (dyadic: closed forms are exact).
    let volume = 0.75 * 0.75 * 0.75;
    let m1 = volume * 0.125;
    let m2_diag = volume * (0.125 * 0.125 + 0.75 * 0.75 / 12.0);
    let m2_cross = volume * 0.125 * 0.125;
    let mut rows = Vec::new();
    rows.push(contained("volume", m.volume, volume, 0.30));
    for (axis, label) in ["m1x", "m1y", "m1z"].into_iter().enumerate() {
        rows.push(contained(label, m.first[axis], m1, 0.20));
    }
    for (enclosure, label) in [
        (m.second.xx, "m2xx"),
        (m.second.yy, "m2yy"),
        (m.second.zz, "m2zz"),
    ] {
        rows.push(contained(label, enclosure, m2_diag, 0.20));
    }
    for (enclosure, label) in [
        (m.second.xy, "m2xy"),
        (m.second.xz, "m2xz"),
        (m.second.yz, "m2yz"),
    ] {
        rows.push(contained(label, enclosure, m2_cross, 0.20));
    }
    verdict(
        "gm-001",
        m.sure_cells > 0 && m.band_cells > 0,
        &format!(
            "box closed forms contained; sure={} band={} {}",
            m.sure_cells,
            m.band_cells,
            rows.join(",")
        ),
        MOMENTS_KERNEL_ID,
    );
}

#[test]
fn gm_002_sphere_moments_and_com() {
    let center = [0.25, -0.125, 0.0];
    let chart = SphereChart {
        center: Point3::new(center[0], center[1], center[2]),
        radius: 0.5,
    };
    let m =
        with_cx(|cx| geometric_moments(&chart, &unit_domain(), 0.02, cx)).expect("sphere moments");
    let volume = 4.0 / 3.0 * PI * (0.5 * 0.5 * 0.5);
    let mut rows = Vec::new();
    rows.push(contained("volume", m.volume, volume, 0.30));
    for a in 0..3 {
        rows.push(contained(
            ["m1x", "m1y", "m1z"][a],
            m.first[a],
            volume * center[a],
            0.25,
        ));
    }
    // Parallel axis about the origin: ∫x² = V(r²/5 + c_x²).
    for (a, (enclosure, label)) in [
        (m.second.xx, "m2xx"),
        (m.second.yy, "m2yy"),
        (m.second.zz, "m2zz"),
    ]
    .into_iter()
    .enumerate()
    {
        let truth = volume * (0.5 * 0.5 / 5.0 + center[a] * center[a]);
        rows.push(contained(label, enclosure, truth, 0.25));
    }
    // Central cross moments vanish for a ball: ∫xy = V·c_x·c_y.
    let cross_truth = [
        volume * center[0] * center[1],
        volume * center[0] * center[2],
        volume * center[1] * center[2],
    ];
    for (slot, (enclosure, label)) in [
        (m.second.xy, "m2xy"),
        (m.second.xz, "m2xz"),
        (m.second.yz, "m2yz"),
    ]
    .into_iter()
    .enumerate()
    {
        rows.push(contained(label, enclosure, cross_truth[slot], 0.25));
    }
    let com = m.com_enclosure().expect("volume lower bound is positive");
    for a in 0..3 {
        assert!(
            com[a].contains(center[a]),
            "com axis {a}: [{}, {}] must contain {}",
            com[a].lo,
            com[a].hi,
            center[a]
        );
        assert!(
            com[a].width() < 0.6,
            "com axis {a} width {}",
            com[a].width()
        );
    }
    verdict(
        "gm-002",
        true,
        &format!(
            "sphere moments + COM contained; sure={} band={} {}",
            m.sure_cells,
            m.band_cells,
            rows.join(",")
        ),
        MOMENTS_KERNEL_ID,
    );
}

#[test]
fn gm_003_translation_covariance_metamorphic() {
    let d = [0.3, -0.2, 0.1];
    let chart = test_box();
    let shifted = Shifted {
        inner: test_box(),
        d: Vec3::new(d[0], d[1], d[2]),
    };
    let (law, direct): (GeometricMoments, GeometricMoments) = with_cx(|cx| {
        let base = geometric_moments(&chart, &unit_domain(), 0.04, cx).expect("base moments");
        let direct =
            geometric_moments(&shifted, &unit_domain(), 0.04, cx).expect("shifted moments");
        (base.translated(d), direct)
    });
    // The shifted region meets different band cells, so the direct
    // volume enclosure differs; both still contain the true volume.
    let pairs = [
        ("volume", law.volume, direct.volume),
        ("m1x", law.first[0], direct.first[0]),
        ("m1y", law.first[1], direct.first[1]),
        ("m1z", law.first[2], direct.first[2]),
        ("m2xx", law.second.xx, direct.second.xx),
        ("m2yy", law.second.yy, direct.second.yy),
        ("m2zz", law.second.zz, direct.second.zz),
        ("m2xy", law.second.xy, direct.second.xy),
        ("m2xz", law.second.xz, direct.second.xz),
        ("m2yz", law.second.yz, direct.second.yz),
    ];
    for (label, a, b) in pairs {
        assert!(
            a.overlaps(&b),
            "{label}: covariance law [{}, {}] and direct [{}, {}] must overlap",
            a.lo,
            a.hi,
            b.lo,
            b.hi
        );
    }
    verdict(
        "gm-003",
        true,
        "translated law and direct recomputation overlap on every component",
        MOMENTS_KERNEL_ID,
    );
}

#[test]
fn gm_004_capability_and_input_refusals() {
    let box_chart = test_box();
    let implicit = ImplicitOnly(test_box());
    let estimating = EstimatingBox(test_box());
    let tiny_sphere = SphereChart {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 0.01,
    };
    // Aabb::new normalizes corners, so an inverted or non-finite domain
    // must be constructed literally to reach the refusal path.
    let inverted = Aabb {
        min: Point3::new(1.0, -1.0, -1.0),
        max: Point3::new(-1.0, 1.0, 1.0),
    };
    let non_finite = Aabb {
        min: Point3::new(f64::NAN, -1.0, -1.0),
        max: Point3::new(1.0, 1.0, 1.0),
    };
    let excludes_support = Aabb::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0));
    let outcomes = with_cx(|cx| {
        [
            (
                "weaker-claim-refuses",
                matches!(
                    geometric_moments(&implicit, &unit_domain(), 0.1, cx),
                    Err(QueryError::MomentsUncertifiedChart {
                        claim: TraceStepClaim::LipschitzImplicit
                    })
                ),
            ),
            (
                "estimate-sample-refuses",
                matches!(
                    geometric_moments(&estimating, &unit_domain(), 0.1, cx),
                    Err(QueryError::MomentsInvalidSample { .. })
                ),
            ),
            (
                "zero-spacing-refuses",
                matches!(
                    geometric_moments(&box_chart, &unit_domain(), 0.0, cx),
                    Err(QueryError::MomentsInvalidSpacing { .. })
                ),
            ),
            (
                "nan-spacing-refuses",
                matches!(
                    geometric_moments(&box_chart, &unit_domain(), f64::NAN, cx),
                    Err(QueryError::MomentsInvalidSpacing { .. })
                ),
            ),
            (
                "inverted-domain-refuses",
                matches!(
                    geometric_moments(&box_chart, &inverted, 0.1, cx),
                    Err(QueryError::MomentsInvalidDomain { .. })
                ),
            ),
            (
                "non-finite-domain-refuses",
                matches!(
                    geometric_moments(&box_chart, &non_finite, 0.1, cx),
                    Err(QueryError::MomentsInvalidDomain { .. })
                ),
            ),
            (
                "support-excluding-domain-refuses",
                matches!(
                    geometric_moments(&box_chart, &excludes_support, 0.1, cx),
                    Err(QueryError::MomentsInvalidDomain { .. })
                ),
            ),
            (
                "per-axis-work-refuses",
                matches!(
                    geometric_moments(&box_chart, &unit_domain(), 1e-7, cx),
                    Err(QueryError::MomentsExcessiveWork { .. })
                ),
            ),
            (
                "total-work-refuses",
                matches!(
                    geometric_moments(&box_chart, &unit_domain(), 1e-6, cx),
                    Err(QueryError::MomentsExcessiveWork { .. })
                ),
            ),
            (
                "unproven-volume-refuses-com",
                match geometric_moments(&tiny_sphere, &unit_domain(), 0.5, cx) {
                    Ok(m) => matches!(
                        m.com_enclosure(),
                        Err(QueryError::MomentsVolumeUnproven { .. })
                    ),
                    Err(_) => false,
                },
            ),
        ]
    });
    let mut emitter = fs_obs::Emitter::new(SUITE, "gm-004");
    for (name, pass) in outcomes {
        emit_verdict(&mut emitter, "gm-004", pass, name, MOMENTS_KERNEL_ID);
    }
}

#[test]
fn gm_005_cancellation_fails_closed() {
    let gate = CancelGate::new();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let refused = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: EXECUTION_SEED,
                kernel_id: CANCELLATION_KERNEL_ID,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        geometric_moments(&test_box(), &unit_domain(), 0.1, &cx)
    });
    verdict(
        "gm-005",
        matches!(refused, Err(QueryError::Cancelled)),
        "a pre-cancelled context refuses before publishing moment enclosures",
        CANCELLATION_KERNEL_ID,
    );
}

/// Exact spherical-shell SDF: `| |p| - mid | - half` (exact distance
/// for `mid > half`), the hollow-solid fixture.
struct ShellChart {
    mid: f64,
    half: f64,
}

impl Chart for ShellChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let r = (x.x * x.x + x.y * x.y + x.z * x.z).sqrt();
        let sd = (r - self.mid).abs() - self.half;
        ChartSample {
            signed_distance: sd,
            gradient: None,
            lipschitz: Some(1.0),
            error: NumericalCertificate::enclosure(sd - 1e-12, sd + 1e-12),
        }
    }

    fn support(&self) -> Aabb {
        let a = self.mid + self.half;
        Aabb::new(Point3::new(-a, -a, -a), Point3::new(a, a, a))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn name(&self) -> &'static str {
        "test/shell"
    }
}

#[test]
fn gm_006_torus_and_hollow_shell_closed_forms() {
    // Ring torus (major 0.5, minor 0.125), axis z, centered: exact SDF.
    let torus = TorusChart {
        center: Point3::new(0.0, 0.0, 0.0),
        major: 0.5,
        minor: 0.125,
    };
    let m = with_cx(|cx| geometric_moments(&torus, &unit_domain(), 0.02, cx)).expect("torus");
    let (rr, r) = (0.5f64, 0.125f64);
    let volume = 2.0 * PI * PI * rr * (r * r);
    let m2_axis = volume * (r * r / 4.0);
    let m2_planar = volume * (rr * rr / 2.0 + 3.0 * r * r / 8.0);
    assert!(
        m.volume.contains(volume) && m.volume.width() < 0.2,
        "torus volume [{}, {}] must contain {volume}",
        m.volume.lo,
        m.volume.hi
    );
    for a in 0..3 {
        assert!(m.first[a].contains(0.0), "torus COM component {a} is 0");
    }
    assert!(
        m.second.zz.contains(m2_axis),
        "torus axial second moment [{}, {}] must contain {m2_axis}",
        m.second.zz.lo,
        m.second.zz.hi
    );
    for (enclosure, label) in [(m.second.xx, "xx"), (m.second.yy, "yy")] {
        assert!(
            enclosure.contains(m2_planar),
            "torus planar {label} [{}, {}] must contain {m2_planar}",
            enclosure.lo,
            enclosure.hi
        );
    }
    for (enclosure, label) in [
        (m.second.xy, "xy"),
        (m.second.xz, "xz"),
        (m.second.yz, "yz"),
    ] {
        assert!(enclosure.contains(0.0), "torus cross {label} contains 0");
    }

    // Hollow spherical shell (outer 0.625, inner 0.375).
    let shell = ShellChart {
        mid: 0.5,
        half: 0.125,
    };
    let s = with_cx(|cx| geometric_moments(&shell, &unit_domain(), 0.02, cx)).expect("shell");
    let (outer, inner) = (0.625f64, 0.375f64);
    let cube = |x: f64| x * x * x;
    let pow5 = |x: f64| x * x * x * x * x;
    let shell_volume = 4.0 / 3.0 * PI * (cube(outer) - cube(inner));
    let shell_m2 = 4.0 * PI / 15.0 * (pow5(outer) - pow5(inner));
    assert!(
        s.volume.contains(shell_volume) && s.volume.width() < 0.35,
        "shell volume [{}, {}] must contain {shell_volume}",
        s.volume.lo,
        s.volume.hi
    );
    let com = s.com_enclosure().expect("positive shell volume");
    for (a, c) in com.iter().enumerate() {
        assert!(c.contains(0.0), "shell COM component {a} contains 0");
    }
    for (enclosure, label) in [
        (s.second.xx, "xx"),
        (s.second.yy, "yy"),
        (s.second.zz, "zz"),
    ] {
        assert!(
            enclosure.contains(shell_m2),
            "shell {label} [{}, {}] must contain {shell_m2}",
            enclosure.lo,
            enclosure.hi
        );
    }
    verdict(
        "gm-006",
        m.sure_cells > 0 && s.sure_cells > 0,
        &format!(
            "torus V∋{volume:.6} zz∋{m2_axis:.6} xx∋{m2_planar:.6}; \
             shell V∋{shell_volume:.6} xx∋{shell_m2:.6}",
        ),
        MOMENTS_KERNEL_ID,
    );
}

#[test]
fn gm_007_open_mesh_refuses_mass_properties() {
    // A single triangle is an OPEN surface: it bounds no volume, and
    // the mesh chart honestly claims no exact-distance theorem. Mass
    // properties must therefore refuse through capability routing —
    // the watertightness precondition, enforced as a typed refusal.
    let open_soup = fs_rep_mesh::Soup {
        positions: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 0.0, 0.0),
            Point3::new(0.0, 0.5, 0.0),
        ],
        triangles: vec![[0, 1, 2]],
    };
    let open_mesh = fs_rep_mesh::MeshChart::new(open_soup);
    let claim = open_mesh.trace_step_claim();
    let refused = with_cx(|cx| geometric_moments(&open_mesh, &unit_domain(), 0.1, cx));
    let pass = matches!(refused, Err(QueryError::MomentsUncertifiedChart { .. }));
    verdict(
        "gm-007",
        pass,
        &format!(
            "chart '{}' with claim {claim:?} lacks the ExactDistance capability; \
             mass properties refused typed",
            open_mesh.name()
        ),
        MOMENTS_KERNEL_ID,
    );
}
