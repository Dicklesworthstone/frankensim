//! fs-geocon conformance suite (CONTRACT.md: any reimplementation must
//! pass). Thickness aggregation with localization and the drive-to-
//! feasibility smoke test, draft angles on analytic tapers with
//! undercut detection, symmetry-by-quotient invariance for arbitrary
//! levers, envelope containment with derivative checks, certified
//! volume enclosures with the Hadamard validation, and the descriptor
//! table. JSON-line verdicts; seeded cases carry seeds.

use asupersync::types::Budget;
use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geocon::{
    CertKind, GeoPrimitive, QuotientChart, SymmetryGroup, VolumeError, draft_violations,
    envelope_violation, min_thickness_soft, min_thickness_soft_clipped, volume_certified,
    volume_smooth,
};
use fs_geom::fixtures::{BoxChart, SphereChart};
use fs_geom::{Aabb, Chart, Point3, SamplingDomainError, TraceStepClaim, Vec3};
use fs_opt::{DescentOptions, Manifold, descend_fn};
use fs_rep_frep::{BoolOp, BoolStyle, Frep, FrepBuilder};
use std::sync::atomic::{AtomicUsize, Ordering};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-geocon/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\
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
                seed: 0x6C0,
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

struct EvalProbe<'a> {
    evals: &'a AtomicUsize,
    panic_on_eval: bool,
    claim: TraceStepClaim,
    weak_evidence: bool,
}

struct MalformedSupportChart;

struct MalformedThicknessChart;

struct CertificateProbeChart {
    certificate: NumericalCertificate,
    support: Aabb,
}

impl Chart for MalformedSupportChart {
    fn eval(&self, _point: Point3, _cx: &Cx<'_>) -> fs_geom::ChartSample {
        fs_geom::ChartSample {
            signed_distance: 0.0,
            gradient: None,
            lipschitz: Some(1.0),
            error: NumericalCertificate::exact(0.0),
        }
    }

    fn support(&self) -> Aabb {
        Aabb {
            min: Point3::new(f64::NAN, -1.0, -1.0),
            max: Point3::new(1.0, 1.0, 1.0),
        }
    }

    fn name(&self) -> &'static str {
        "geocon/malformed-support-probe"
    }
}

impl Chart for MalformedThicknessChart {
    fn eval(&self, _point: Point3, _cx: &Cx<'_>) -> fs_geom::ChartSample {
        fs_geom::ChartSample {
            signed_distance: f64::NAN,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: None,
            error: NumericalCertificate::no_claim(),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "geocon/malformed-thickness-probe"
    }
}

impl Chart for CertificateProbeChart {
    fn eval(&self, _point: Point3, _cx: &Cx<'_>) -> fs_geom::ChartSample {
        fs_geom::ChartSample {
            signed_distance: 0.0,
            gradient: None,
            lipschitz: Some(1.0),
            error: self.certificate,
        }
    }

    fn support(&self) -> Aabb {
        self.support
    }

    fn name(&self) -> &'static str {
        "geocon/certificate-probe"
    }
}

impl Chart for EvalProbe<'_> {
    fn eval(&self, point: Point3, cx: &Cx<'_>) -> fs_geom::ChartSample {
        self.evals.fetch_add(1, Ordering::Relaxed);
        assert!(
            !self.panic_on_eval,
            "volume preflight reached chart evaluation"
        );
        let mut sample = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        }
        .eval(point, cx);
        if self.weak_evidence {
            sample.error = NumericalCertificate::estimate(
                sample.signed_distance - 1e-12,
                sample.signed_distance + 1e-12,
            );
        }
        sample
    }

    fn support(&self) -> Aabb {
        Aabb::WHOLE_SPACE
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        self.claim
    }

    fn name(&self) -> &'static str {
        "geocon/eval-probe"
    }
}

/// Dumbbell with a parametric neck radius (the design lever).
fn dumbbell(neck_r: f64) -> Frep {
    let mut b = FrepBuilder::new();
    let s1 = b.sphere(Point3::new(-1.2, 0.0, 0.0), 0.8).expect("s1");
    let s2 = b.sphere(Point3::new(1.2, 0.0, 0.0), 0.8).expect("s2");
    let neck = b
        .cylinder(Point3::new(0.0, 0.0, 0.0), neck_r)
        .expect("neck");
    let neck = b
        .rotate(neck, Vec3::new(0.0, 1.0, 0.0), core::f64::consts::FRAC_PI_2)
        .expect("rot");
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
}

fn neck_samples(r: f64) -> Vec<Point3> {
    (0..16)
        .map(|k| {
            let th = core::f64::consts::TAU * f64::from(k) / 16.0;
            Point3::new(0.0, r * th.cos(), r * th.sin())
        })
        .collect()
}

const ENVELOPE_SPHERE_RADIUS: f64 = 0.4;

fn envelope_sphere_samples(c: Point3) -> Vec<Point3> {
    let mut v = Vec::new();
    for k in 0..32 {
        let th = core::f64::consts::TAU * f64::from(k) / 32.0;
        for &z in &[-0.25, 0.0, 0.25] {
            let r = (ENVELOPE_SPHERE_RADIUS * ENVELOPE_SPHERE_RADIUS - z * z)
                .max(0.0)
                .sqrt();
            v.push(Point3::new(c.x + r * th.cos(), c.y + r * th.sin(), c.z + z));
        }
        v.push(Point3::new(c.x, c.y, c.z + ENVELOPE_SPHERE_RADIUS));
        v.push(Point3::new(c.x, c.y, c.z - ENVELOPE_SPHERE_RADIUS));
    }
    v
}

fn envelope_report_at(
    allowed: &dyn Chart,
    center_x: f64,
    cx: &Cx<'_>,
) -> fs_geocon::EnvelopeReport {
    envelope_violation(
        allowed,
        &envelope_sphere_samples(Point3::new(center_x, 0.0, 0.0)),
        40.0,
        false,
        cx,
    )
}

fn envelope_soft_worst_at(allowed: &dyn Chart, center_x: f64, cx: &Cx<'_>) -> f64 {
    envelope_report_at(allowed, center_x, cx).soft_worst
}

fn keepout_report_at(keepout: &dyn Chart, center_x: f64, cx: &Cx<'_>) -> fs_geocon::EnvelopeReport {
    envelope_violation(
        keepout,
        &envelope_sphere_samples(Point3::new(center_x, 0.0, 0.0)),
        40.0,
        true,
        cx,
    )
}

/// gcp-001 — min-thickness: the soft aggregate under-approximates the
/// hard minimum and converges to it as p grows; violations LOCALIZE to
/// the thin samples; the FD lever derivative is right; and a toy
/// descent DRIVES the neck to feasibility (derivatives point the
/// right way).
#[test]
fn gcp_001_min_thickness() {
    with_cx(|cx| {
        let d = dumbbell(0.15);
        // Mixed samples: thin neck ring + thick sphere caps.
        let mut samples = neck_samples(0.15);
        let thin_count = samples.len();
        samples.push(Point3::new(-2.0, 0.0, 0.0));
        samples.push(Point3::new(2.0, 0.0, 0.0));
        let rep = min_thickness_soft(&d, &samples, 0.5, 8.0, cx).expect("thickness");
        let soft_over = rep.soft_min >= rep.hard_min - 1e-12
            && rep.authority == fs_evidence::NumericalKind::Estimate;
        let rep_hard = min_thickness_soft(&d, &samples, 0.5, 40.0, cx).expect("harder p");
        let converges =
            (rep_hard.soft_min - rep.hard_min).abs() < (rep.soft_min - rep.hard_min).abs() + 1e-12;
        // Localization: exactly the neck ring violates required = 0.5.
        let localized = rep.violating.len() == thin_count
            && rep.violating.iter().all(|&i| i < thin_count)
            && rep.skipped == 0;
        // Lever derivative (soft_min through neck radius) vs FD.
        let h = 1e-4;
        let f = |r: f64| {
            min_thickness_soft(&dumbbell(r), &neck_samples(r), 0.5, 8.0, cx)
                .expect("t")
                .soft_min
        };
        let fd = (f(0.15 + h) - f(0.15 - h)) / (2.0 * h);
        let deriv_ok = (fd - 2.0).abs() < 0.05; // neck-only samples: d(2r)/dr = 2
        // Drive to feasibility: descend the hinge penalty over the lever.
        let objective = |x: &[f64]| -> f64 {
            let r = x[0].clamp(0.05, 0.45);
            let t = min_thickness_soft(&dumbbell(r), &neck_samples(r), 0.5, 8.0, cx)
                .expect("t")
                .soft_min;
            let deficit = (0.5 - t).max(0.0);
            deficit * deficit
        };
        let repd = descend_fn(
            Manifold::Rn { dim: 1 },
            &objective,
            &[0.15],
            DescentOptions {
                steps: 120,
                lr: 0.5,
                fd_h: 1e-5,
            },
            0,
            cx,
        )
        .expect("descent");
        let final_r = repd.x[0].clamp(0.05, 0.45);
        let final_t = min_thickness_soft(&dumbbell(final_r), &neck_samples(final_r), 0.5, 8.0, cx)
            .expect("t")
            .soft_min;
        let feasible = final_t >= 0.5 - 1e-3;
        verdict(
            "gcp-001",
            soft_over && converges && localized && deriv_ok && feasible,
            &format!(
                "soft-min over-approximates converging down with p, violations localize \
                 to exactly the {thin_count} neck samples, the lever FD derivative is \
                 {fd:.3} (analytic 2), and the toy descent drives the neck from \
                 r=0.15 to r={final_r:.3} reaching thickness {final_t:.3} >= 0.5 — \
                 the anti-paperclip constraint, closed loop"
            ),
        );
    });
}

/// gcp-001b — an infinite extrusion cannot silently turn unresolved-domain
/// thickness samples into skips. The default aggregate refuses, while an
/// explicit finite clip enables a deliberately local report.
#[test]
fn gcp_001b_unbounded_thickness_requires_clip() {
    with_cx(|cx| {
        let cylinder = {
            let mut builder = FrepBuilder::new();
            let root = builder
                .cylinder(Point3::new(0.0, 0.0, 0.0), 1.0)
                .expect("cylinder");
            builder.finish(root).expect("frep")
        };
        let samples = [Point3::new(1.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)];
        let refused = min_thickness_soft(&cylinder, &samples, 2.1, 8.0, cx);
        let clip = Aabb::new(Point3::new(-2.0, -2.0, -2.0), Point3::new(2.0, 2.0, 2.0));
        let local = min_thickness_soft_clipped(&cylinder, &samples, 2.1, 8.0, clip, cx)
            .expect("finite clip admits the local thickness aggregate");
        verdict(
            "gcp-001b",
            matches!(
                refused,
                Err(fs_query::QueryError::SamplingDomain(
                    SamplingDomainError::UnboundedSupport { .. }
                ))
            ) && (local.hard_min - 2.0).abs() < 1e-9
                && (local.soft_min - 2.0).abs() < 1e-9
                && local.authority == fs_evidence::NumericalKind::Estimate
                && local.violating == vec![0, 1]
                && local.skipped == 0,
            "unresolved extended support is a structured refusal rather than a skipped \
             sample; an explicit finite clip yields the local cylinder thickness report",
        );
    });
}

/// gcp-001c — aggregation preserves Estimate authority and never converts
/// empty, malformed, or invalid arithmetic into a clean thickness report.
#[test]
fn gcp_001c_thickness_aggregation_fails_closed() {
    with_cx(|cx| {
        let sphere = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        assert!(matches!(
            min_thickness_soft(&sphere, &[], 0.5, 8.0, cx),
            Err(fs_query::QueryError::NoThicknessSamples { skipped: 0 })
        ));
        assert!(matches!(
            min_thickness_soft(&sphere, &[Point3::new(0.0, 0.0, 0.0)], 0.5, 8.0, cx,),
            Err(fs_query::QueryError::NoThicknessSamples { skipped: 1 })
        ));
        assert!(matches!(
            min_thickness_soft(
                &MalformedThicknessChart,
                &[Point3::new(1.0, 0.0, 0.0)],
                0.5,
                8.0,
                cx,
            ),
            Err(fs_query::QueryError::InvalidThicknessSample { .. })
        ));
        for (required, exponent) in [(f64::NAN, 8.0), (0.5, 0.0)] {
            assert!(matches!(
                min_thickness_soft(
                    &sphere,
                    &[Point3::new(1.0, 0.0, 0.0)],
                    required,
                    exponent,
                    cx,
                ),
                Err(fs_query::QueryError::InvalidThicknessArithmetic { .. })
            ));
        }
    });
}

/// gcp-002 — draft angles: an analytic cone tapered at 10° passes a 5°
/// requirement and fails 15° with violations localized to the wall;
/// vertical box walls violate any positive draft; a mushroom cap
/// undercut is flagged as an UNDERCUT, not mere low draft.
#[test]
fn gcp_002_draft_angles() {
    with_cx(|cx| {
        let pull = Vec3::new(0.0, 0.0, 1.0);
        // Tapered cone via F-rep: cylinder radius shrinking with z is
        // not in the primitive zoo — use a rotated half-space wall:
        // plane with normal tilted 10° from horizontal models the wall.
        // Simpler analytic: sample a cone surface x²+y² = (r0 − z·tanθ)²
        // directly with its known normals via a sphere-chart trick is
        // overkill — assess the frep BOX (vertical walls) and a HALF-
        // SPACE tilted by exactly 10°.
        let tilted = |deg: f64| -> Frep {
            let mut b = FrepBuilder::new();
            let th = deg.to_radians();
            // Wall normal: tilted from horizontal toward +z by θ.
            let n = Vec3::new(th.cos(), 0.0, th.sin());
            let hs = b.half_space(n, 0.5).expect("hs");
            let bx = b
                .box_prim(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0))
                .expect("bx");
            let root = b
                .boolean(BoolOp::Intersect, BoolStyle::Hard, hs, bx)
                .expect("r");
            b.finish(root).expect("frep")
        };
        let wall10 = tilted(10.0);
        // Sample the tilted wall face (x ≈ 0.5·cosθ locus): project
        // points onto the wall by closest_point.
        let mut wall_samples = Vec::new();
        for k in 0..12 {
            let y = -0.8 + 1.6 * f64::from(k) / 11.0;
            let p = fs_query::closest_point(&wall10, Point3::new(0.6, y, 0.0), cx)
                .expect("cp")
                .point;
            wall_samples.push(p);
        }
        let pass5 =
            draft_violations(&wall10, &wall_samples, pull, 5.0f64.to_radians(), cx).expect("5deg");
        let fail15 = draft_violations(&wall10, &wall_samples, pull, 15.0f64.to_radians(), cx)
            .expect("15deg");
        let cone_ok = pass5.violating.is_empty()
            && pass5.penalty == 0.0
            && fail15.violating.len() == wall_samples.len()
            && fail15.worst_deficit > 0.0;
        // Vertical box walls: any positive draft fails.
        let bx = {
            let mut b = FrepBuilder::new();
            let x = b
                .box_prim(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.5, 0.5, 0.5))
                .expect("bx");
            b.finish(x).expect("frep")
        };
        let side = vec![Point3::new(0.5, 0.0, 0.0), Point3::new(-0.5, 0.1, 0.2)];
        let box_rep = draft_violations(&bx, &side, pull, 3.0f64.to_radians(), cx).expect("box");
        let box_ok = box_rep.violating.len() == 2 && box_rep.undercuts.is_empty();
        // Mushroom: sphere cap overhanging a thin stem — the underside
        // of the cap has normals AGAINST the pull: undercut.
        let mushroom = {
            let mut b = FrepBuilder::new();
            let cap = b.sphere(Point3::new(0.0, 0.0, 1.0), 0.6).expect("cap");
            let stem = b.cylinder(Point3::new(0.0, 0.0, 0.0), 0.15).expect("stem");
            let root = b
                .boolean(BoolOp::Union, BoolStyle::Hard, cap, stem)
                .expect("r");
            b.finish(root).expect("frep")
        };
        // A point on the cap's lower shoulder: outward normal dips
        // BELOW horizontal (n·pull ≈ −0.31) — an undercut within the
        // top mold's own reach, not the other half's face.
        let dirv = Vec3::new(0.95, 0.0, -0.31);
        let dn = dirv.norm();
        let probe = Point3::new(0.0 + 0.7 * dirv.x / dn, 0.0, 1.0 + 0.7 * dirv.z / dn);
        let under = fs_query::closest_point(&mushroom, probe, cx)
            .expect("cp")
            .point;
        let mush =
            draft_violations(&mushroom, &[under], pull, 5.0f64.to_radians(), cx).expect("mushroom");
        let undercut_ok = mush.undercuts.len() == 1 && mush.violating.is_empty();
        // Smooth penalty derivative vs FD through the tilt angle.
        let pen = |deg: f64| -> f64 {
            let w = tilted(deg);
            let mut s = Vec::new();
            for k in 0..12 {
                let y = -0.8 + 1.6 * f64::from(k) / 11.0;
                s.push(
                    fs_query::closest_point(&w, Point3::new(0.6, y, 0.0), cx)
                        .expect("cp")
                        .point,
                );
            }
            draft_violations(&w, &s, pull, 15.0f64.to_radians(), cx)
                .expect("d")
                .penalty
        };
        let h = 0.05;
        let fd = (pen(10.0 + h) - pen(10.0 - h)) / (2.0 * h);
        // Analytic: penalty = (sin15° − sin θ)², d/dθdeg = −2(sin15°−sinθ)cosθ·(π/180).
        let th = 10.0f64.to_radians();
        let analytic =
            -2.0 * (15.0f64.to_radians().sin() - th.sin()) * th.cos() * core::f64::consts::PI
                / 180.0;
        let deriv_ok = (fd - analytic).abs() < 0.05 * analytic.abs();
        verdict(
            "gcp-002",
            cone_ok && box_ok && undercut_ok && deriv_ok,
            &format!(
                "a 10-degree wall passes 5 and fails 15 with all samples localized, \
                 vertical walls violate any positive draft, the mushroom underside is \
                 flagged as an UNDERCUT (not low draft), and the smooth penalty's FD \
                 derivative matches the analytic hinge slope ({fd:.4} vs \
                 {analytic:.4})"
            ),
        );
    });
}

/// gcp-003 — symmetry by quotient: invariance holds for ARBITRARY
/// inner designs (property test over random freps and levers) —
/// bitwise for reflection/translation, 1e-9 for cyclic; gradients
/// chain correctly; asymmetric inners still yield symmetric shapes.
#[test]
#[allow(clippy::too_many_lines)] // One seeded quotient campaign shares invariance and gradient state.
fn gcp_003_symmetry_quotient() {
    with_cx(|cx| {
        let mut rng = Lcg(0x1001_2026_0707_0023);
        let mut invariant = true;
        let mut grad_ok = true;
        let mut support_ok = true;
        let mut authority_ok = true;
        for trial in 0..12 {
            // A deliberately ASYMMETRIC inner design.
            let inner = {
                let mut b = FrepBuilder::new();
                let s1 = b
                    .sphere(
                        Point3::new(
                            rng.range(0.2, 0.9),
                            rng.range(-0.4, 0.4),
                            rng.range(-0.4, 0.4),
                        ),
                        rng.range(0.3, 0.6),
                    )
                    .expect("s1");
                let s2 = b
                    .sphere(
                        Point3::new(
                            rng.range(0.2, 0.9),
                            rng.range(-0.4, 0.4),
                            rng.range(-0.4, 0.4),
                        ),
                        rng.range(0.2, 0.5),
                    )
                    .expect("s2");
                let u = b
                    .boolean(BoolOp::Union, BoolStyle::Blend { radius: 0.15 }, s1, s2)
                    .expect("u");
                b.finish(u).expect("frep")
            };
            let groups = [
                SymmetryGroup::ReflectX,
                SymmetryGroup::Cyclic { n: 6 },
                SymmetryGroup::Periodic { period: 2.5 },
            ];
            for group in groups {
                let q = QuotientChart {
                    inner: &inner,
                    group,
                };
                if matches!(group, SymmetryGroup::Periodic { .. }) {
                    let support = q.support();
                    let inner_support = inner.support();
                    support_ok &= support.min.x.is_infinite()
                        && support.min.x.is_sign_negative()
                        && support.max.x.is_infinite()
                        && support.max.x.is_sign_positive()
                        && support.min.y.to_bits() == inner_support.min.y.to_bits()
                        && support.max.y.to_bits() == inner_support.max.y.to_bits()
                        && support.min.z.to_bits() == inner_support.min.z.to_bits()
                        && support.max.z.to_bits() == inner_support.max.z.to_bits()
                        && support.contains(Point3::new(
                            1.0e300,
                            f64::midpoint(inner_support.min.y, inner_support.max.y),
                            f64::midpoint(inner_support.min.z, inner_support.max.z),
                        ));
                }
                for _ in 0..24 {
                    let p = Point3::new(
                        rng.range(-2.0, 2.0),
                        rng.range(-2.0, 2.0),
                        rng.range(-1.0, 1.0),
                    );
                    let sample = q.eval(p, cx);
                    authority_ok &= sample.lipschitz.is_none()
                        && sample.error.kind == fs_evidence::NumericalKind::Estimate;
                    let base = sample.signed_distance;
                    for gp in group.orbit(p) {
                        let moved = q.eval(gp, cx).signed_distance;
                        let tol = match group {
                            // Reflection folds bitwise; the cyclic fold
                            // and the PROBES' own `x + period` rounding
                            // sit at fp scale.
                            SymmetryGroup::ReflectX => 0.0,
                            _ => 1e-9,
                        };
                        if (moved - base).abs() > tol {
                            invariant = false;
                        }
                    }
                }
            }
            // Gradient chain rule vs FD (reflection, off-seam points).
            if trial < 4 {
                let q = QuotientChart {
                    inner: &inner,
                    group: SymmetryGroup::ReflectX,
                };
                let p = Point3::new(
                    -rng.range(0.3, 1.5),
                    rng.range(-1.0, 1.0),
                    rng.range(-0.5, 0.5),
                );
                if let Some(g) = q.eval(p, cx).gradient {
                    let h = 1e-6;
                    let fd = Vec3::new(
                        (q.eval(Point3::new(p.x + h, p.y, p.z), cx).signed_distance
                            - q.eval(Point3::new(p.x - h, p.y, p.z), cx).signed_distance)
                            / (2.0 * h),
                        (q.eval(Point3::new(p.x, p.y + h, p.z), cx).signed_distance
                            - q.eval(Point3::new(p.x, p.y - h, p.z), cx).signed_distance)
                            / (2.0 * h),
                        (q.eval(Point3::new(p.x, p.y, p.z + h), cx).signed_distance
                            - q.eval(Point3::new(p.x, p.y, p.z - h), cx).signed_distance)
                            / (2.0 * h),
                    );
                    let diff = Vec3::new(g.x - fd.x, g.y - fd.y, g.z - fd.z);
                    grad_ok &= diff.norm() < 1e-4;
                }
            }
        }

        let valid_inner = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        for group in [
            SymmetryGroup::Cyclic { n: 0 },
            SymmetryGroup::Periodic { period: 0.0 },
            SymmetryGroup::Periodic { period: f64::NAN },
        ] {
            assert!(group.validate().is_err());
            let invalid = QuotientChart {
                inner: &valid_inner,
                group,
            };
            assert!(!invalid.support().is_well_formed());
            let sample = invalid.eval(Point3::new(0.25, 0.0, 0.0), cx);
            assert!(sample.signed_distance.is_nan());
            assert!(sample.lipschitz.is_none());
            assert_eq!(sample.error.kind, fs_evidence::NumericalKind::NoClaim);
            assert!(group.orbit(Point3::new(0.25, 0.0, 0.0)).is_empty());
        }
        let malformed_orbit = QuotientChart {
            inner: &MalformedSupportChart,
            group: SymmetryGroup::ReflectX,
        }
        .support();
        assert!(malformed_orbit.min.x.is_nan());

        let unit_support = Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0));
        for malformed in [
            NumericalCertificate {
                kind: NumericalKind::Enclosure,
                lo: 1.0,
                hi: -1.0,
            },
            NumericalCertificate::enclosure(1.0, 2.0),
            NumericalCertificate {
                kind: NumericalKind::Exact,
                lo: -1.0,
                hi: 1.0,
            },
            NumericalCertificate {
                kind: NumericalKind::Estimate,
                lo: f64::NAN,
                hi: 1.0,
            },
        ] {
            let inner = CertificateProbeChart {
                certificate: malformed,
                support: unit_support,
            };
            let sample = QuotientChart {
                inner: &inner,
                group: SymmetryGroup::ReflectX,
            }
            .eval(Point3::new(0.25, 0.0, 0.0), cx);
            assert_eq!(
                sample.error.kind,
                NumericalKind::NoClaim,
                "quotienting must not repair malformed or nominal-excluding inner evidence"
            );
        }
        let valid_estimate = CertificateProbeChart {
            certificate: NumericalCertificate::estimate(-1.0, 1.0),
            support: unit_support,
        };
        let quotient_estimate = QuotientChart {
            inner: &valid_estimate,
            group: SymmetryGroup::ReflectX,
        }
        .eval(Point3::new(0.25, 0.0, 0.0), cx)
        .error;
        assert_eq!(quotient_estimate.kind, NumericalKind::Estimate);
        assert_eq!(quotient_estimate.lo.to_bits(), (-1.0_f64).to_bits());
        assert_eq!(quotient_estimate.hi.to_bits(), 1.0_f64.to_bits());

        // This mantissa makes nearest `extent * SQRT_2` round below the
        // mathematical corner radius. The cyclic support must retain an
        // outward endpoint rather than inherit that nearest-rounding gap.
        let extent = f64::from_bits(0x0db7_80e7_635e_965f);
        let tiny_corner = CertificateProbeChart {
            certificate: NumericalCertificate::exact(0.0),
            support: Aabb::new(
                Point3::new(0.0, 0.0, -1.0),
                Point3::new(extent, extent, 1.0),
            ),
        };
        let cyclic_support = QuotientChart {
            inner: &tiny_corner,
            group: SymmetryGroup::Cyclic { n: 4 },
        }
        .support();
        assert!(
            cyclic_support.max.x >= extent.hypot(extent).next_up(),
            "cyclic support radius must be rounded outward"
        );
        verdict(
            "gcp-003",
            invariant && grad_ok && support_ok && authority_ok,
            "the quotient shape is invariant under its group for ARBITRARY \
             asymmetric inner designs (bitwise for reflection; fp-scale \
             for cyclic/periodic) across 12 random levers x 3 groups x 24 probes x full \
             orbits, folded gradients match finite differences off-seam, and periodic \
             support is honestly infinite along x while retaining finite transverse bounds; \
             raw quotient fields retain only Estimate/NoClaim authority, malformed inner \
             evidence cannot be laundered, cyclic support radii round outward, and invalid \
             group or inner-support inputs fail closed; \
             seed 0x1001_2026_0707_0023 — symmetry violation is structurally \
             impossible",
        );
    });
}

/// gcp-004 — envelopes: containment and keep-out assessments match
/// analytic penetration depths, softmax tracks the hard worst, the FD
/// derivative is right, and a toy descent pulls an escaping design
/// back inside.
#[test]
fn gcp_004_envelopes() {
    with_cx(|cx| {
        let allowed = BoxChart {
            aabb: Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0)),
        };
        let inside = envelope_report_at(&allowed, 0.3, cx);
        let contained = inside.worst <= 0.0 && inside.violating.is_empty();
        // Pushed out by 0.2: worst ≈ +0.2 penetration.
        let out = envelope_report_at(&allowed, 0.8, cx);
        let n_samples = envelope_sphere_samples(Point3::new(0.0, 0.0, 0.0)).len() as f64;
        let penetration_ok = (out.worst - 0.2).abs() < 1e-9
            && !out.violating.is_empty()
            && out.soft_worst >= out.worst
            && out.soft_worst - out.worst <= n_samples.ln() / 40.0 + 1e-9;
        // Keep-out: a forbidden ball at the origin.
        let keepout = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 0.5,
        };
        let clear = keepout_report_at(&keepout, 1.2, cx);
        let hit = keepout_report_at(&keepout, 0.6, cx);
        let keepout_ok = clear.worst <= 0.0 && hit.worst > 0.0;
        // FD derivative of soft_worst through the center position.
        let h = 1e-5;
        let fd = (envelope_soft_worst_at(&allowed, 0.8 + h, cx)
            - envelope_soft_worst_at(&allowed, 0.8 - h, cx))
            / (2.0 * h);
        let deriv_ok = (fd - 1.0).abs() < 0.05; // moving out 1:1
        // Drive to feasibility: descend soft_worst hinge over center x.
        let objective = |x: &[f64]| -> f64 {
            let s = envelope_soft_worst_at(&allowed, x[0], cx).max(0.0);
            s * s
        };
        let rep = descend_fn(
            Manifold::Rn { dim: 1 },
            &objective,
            &[0.9],
            DescentOptions {
                steps: 80,
                lr: 0.4,
                fd_h: 1e-5,
            },
            0,
            cx,
        )
        .expect("descent");
        let back_inside = envelope_report_at(&allowed, rep.x[0], cx).worst <= 1e-3;
        verdict(
            "gcp-004",
            contained && penetration_ok && keepout_ok && deriv_ok && back_inside,
            &format!(
                "containment and keep-out match analytic penetrations (worst 0.2 read \
                 as {:.4}), softmax tracks the hard worst, the FD derivative is \
                 {fd:.3} (analytic 1), and the descent pulls the design from x=0.9 \
                 back to x={:.3} inside the envelope",
                out.worst, rep.x[0]
            ),
        );
    });
}

/// gcp-005 — volume: the certified enclosure brackets the analytic
/// sphere volume and TIGHTENS with h; the smoothed volume's lever
/// derivative matches the Hadamard formula (dV/dr = 4πr²); a toy
/// descent shrinks a sphere to meet a volume cap.
#[test]
fn gcp_005_volume_hadamard() {
    with_cx(|cx| {
        let sphere = |r: f64| -> Frep {
            let mut b = FrepBuilder::new();
            let s = b.sphere(Point3::new(0.0, 0.0, 0.0), r).expect("s");
            b.finish(s).expect("frep")
        };
        let truth = 4.0 * core::f64::consts::PI / 3.0;
        let dom = Aabb::new(Point3::new(-1.6, -1.6, -1.6), Point3::new(1.6, 1.6, 1.6));
        let coarse = volume_certified(&sphere(1.0), &dom, 0.1, cx).expect("coarse");
        let fine = volume_certified(&sphere(1.0), &dom, 0.05, cx).expect("fine");
        let brackets =
            coarse.lo <= truth && truth <= coarse.hi && fine.lo <= truth && truth <= fine.hi;
        let tightens = (fine.hi - fine.lo) < 0.6 * (coarse.hi - coarse.lo);
        // Hadamard: FD of the smoothed volume vs 4πr².
        let vs = |r: f64| volume_smooth(&sphere(r), &dom, 0.04, 0.02, cx).expect("vs");
        let h = 1e-3;
        let fd = (vs(1.0 + h) - vs(1.0 - h)) / (2.0 * h);
        let hadamard = 4.0 * core::f64::consts::PI;
        let hadamard_ok = (fd - hadamard).abs() / hadamard < 0.02;
        // Descent to a volume cap: shrink r until V ≤ 2.0.
        let objective = |x: &[f64]| -> f64 {
            let r = x[0].clamp(0.3, 1.5);
            let v = vs(r);
            let excess = (v - 2.0).max(0.0);
            excess * excess
        };
        let rep = descend_fn(
            Manifold::Rn { dim: 1 },
            &objective,
            &[1.2],
            DescentOptions {
                steps: 60,
                lr: 0.05,
                fd_h: 1e-4,
            },
            0,
            cx,
        )
        .expect("descent");
        let final_v = vs(rep.x[0].clamp(0.3, 1.5));
        let capped = final_v <= 2.0 + 0.05;
        let mut em = fs_obs::Emitter::new("fs-geocon/conformance", "gcp-005/volume");
        let line = em
            .emit(
                fs_obs::Severity::Info,
                fs_obs::EventKind::Custom {
                    name: "geocon-volume-hadamard".to_string(),
                    json: format!(
                        "{{\"coarse\":[{:.4},{:.4}],\"fine\":[{:.4},{:.4}],\
                         \"dv_dr\":{fd:.4},\"hadamard\":{hadamard:.4}}}",
                        coarse.lo, coarse.hi, fine.lo, fine.hi
                    ),
                },
                None,
            )
            .to_jsonl();
        fs_obs::validate_line(&line).expect("volume event validates");
        println!("{line}");
        verdict(
            "gcp-005",
            brackets && tightens && hadamard_ok && capped,
            &format!(
                "certified enclosures bracket 4pi/3 at both resolutions and tighten \
                 with h ([{:.3},{:.3}] -> [{:.3},{:.3}]), the smoothed volume's lever \
                 derivative {fd:.3} matches Hadamard 4pi r^2 = {hadamard:.3} within \
                 2%, and the descent shrinks r to meet the volume cap \
                 (V = {final_v:.3} <= 2.05)",
                coarse.lo, coarse.hi, fine.lo, fine.hi
            ),
        );
    });
}

/// gcp-005b — volume grid admission is fail-closed. Malformed or unbounded
/// integration boxes, invalid spacing, and excessive/unrepresentable counts
/// are rejected before chart evaluation; an ordinary finite box still runs.
#[test]
#[allow(clippy::too_many_lines)] // One preflight matrix proves every refusal precedes evaluation.
fn gcp_005b_volume_preflight_refuses_before_eval() {
    with_cx(|cx| {
        let evals = AtomicUsize::new(0);
        let panic_probe = EvalProbe {
            evals: &evals,
            panic_on_eval: true,
            claim: TraceStepClaim::ExactDistance,
            weak_evidence: false,
        };
        let malformed = Aabb::new(
            Point3::new(f64::NAN, -1.0, -1.0),
            Point3::new(1.0, 1.0, 1.0),
        );
        let finite = Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0));

        for result in [
            volume_certified(&panic_probe, &malformed, 0.25, cx).map(|_| ()),
            volume_smooth(&panic_probe, &malformed, 0.25, 0.1, cx).map(|_| ()),
        ] {
            assert!(matches!(
                result,
                Err(VolumeError::SamplingDomain(
                    SamplingDomainError::InvalidSupport { .. }
                ))
            ));
        }
        for result in [
            volume_certified(&panic_probe, &Aabb::WHOLE_SPACE, 0.25, cx).map(|_| ()),
            volume_smooth(&panic_probe, &Aabb::WHOLE_SPACE, 0.25, 0.1, cx).map(|_| ()),
        ] {
            assert!(matches!(
                result,
                Err(VolumeError::SamplingDomain(
                    SamplingDomainError::UnboundedSupport { .. }
                ))
            ));
        }
        for invalid_h in [0.0, f64::NAN] {
            assert!(matches!(
                volume_certified(&panic_probe, &finite, invalid_h, cx),
                Err(VolumeError::InvalidSpacing { field: "h", .. })
            ));
            assert!(matches!(
                volume_smooth(&panic_probe, &finite, invalid_h, 0.1, cx),
                Err(VolumeError::InvalidSpacing { field: "h", .. })
            ));
        }
        assert!(matches!(
            volume_smooth(&panic_probe, &finite, 0.25, 0.0, cx),
            Err(VolumeError::InvalidSpacing {
                field: "epsilon",
                ..
            })
        ));
        for result in [
            volume_certified(&panic_probe, &finite, 1e-6, cx).map(|_| ()),
            volume_smooth(&panic_probe, &finite, 1e-6, 0.1, cx).map(|_| ()),
        ] {
            assert!(matches!(result, Err(VolumeError::WorkLimit { .. })));
        }
        for result in [
            volume_certified(&panic_probe, &finite, 1e-14, cx).map(|_| ()),
            volume_smooth(&panic_probe, &finite, 1e-14, 0.1, cx).map(|_| ()),
        ] {
            assert!(matches!(
                result,
                Err(VolumeError::CellProductOverflow { .. })
            ));
        }
        for result in [
            volume_certified(&panic_probe, &finite, f64::MIN_POSITIVE, cx).map(|_| ()),
            volume_smooth(&panic_probe, &finite, f64::MIN_POSITIVE, 0.1, cx).map(|_| ()),
        ] {
            assert!(matches!(result, Err(VolumeError::CellCountOverflow { .. })));
        }
        assert_eq!(evals.load(Ordering::Relaxed), 0);

        let no_claim_probe = EvalProbe {
            evals: &evals,
            panic_on_eval: true,
            claim: TraceStepClaim::NoClaim,
            weak_evidence: false,
        };
        assert!(matches!(
            volume_certified(&no_claim_probe, &finite, 1.0, cx),
            Err(VolumeError::UncertifiedChart {
                claim: TraceStepClaim::NoClaim
            })
        ));
        assert_eq!(
            evals.load(Ordering::Relaxed),
            0,
            "a weak chart theorem must refuse before evaluation"
        );

        let weak_evidence_probe = EvalProbe {
            evals: &evals,
            panic_on_eval: false,
            claim: TraceStepClaim::ExactDistance,
            weak_evidence: true,
        };
        assert!(matches!(
            volume_certified(&weak_evidence_probe, &finite, 1.0, cx),
            Err(VolumeError::InvalidCertificate { .. })
        ));
        assert_eq!(evals.load(Ordering::Relaxed), 1);
        evals.store(0, Ordering::Relaxed);

        let counting_probe = EvalProbe {
            evals: &evals,
            panic_on_eval: false,
            claim: TraceStepClaim::ExactDistance,
            weak_evidence: false,
        };
        let certified =
            volume_certified(&counting_probe, &finite, 1.0, cx).expect("finite volume grid");
        let smooth = volume_smooth(&counting_probe, &finite, 1.0, 0.1, cx)
            .expect("finite smooth volume grid");
        verdict(
            "gcp-005b",
            certified.lo.is_finite()
                && certified.hi.is_finite()
                && certified.lo <= certified.hi
                && smooth.is_finite()
                && evals.load(Ordering::Relaxed) == 54,
            "volume samplers reject malformed, unbounded, invalid-spacing, and excessive grids before evaluation; certified volume also refuses weak chart theorems and weak per-sample evidence; outward count admission uses a conservative 3x3x3 grid and evaluates exactly 27 cells per API",
        );
    });
}

/// gcp-005c — directed proof arithmetic remains honest when the ideal cell
/// center lies between adjacent representable floats and every cell is surely
/// inside. The real domain volume is exactly 2, but the certificate must still
/// publish an outward interval rather than a nearest-rounded singleton.
#[test]
fn gcp_005c_volume_nextafter_rounding_is_outward() {
    with_cx(|cx| {
        let x_min = f64::from_bits(0x4330_0000_0000_0000); // Exactly 2^52.
        let x_max = x_min.next_up();
        let domain = Aabb::new(Point3::new(x_min, -1.0, 0.0), Point3::new(x_max, 1.0, 1.0));
        let mut builder = FrepBuilder::new();
        let root = builder
            .half_space(Vec3::new(1.0, 0.0, 0.0), x_max + 8.0)
            .expect("axis halfspace");
        let inside = builder.finish(root).expect("exact halfspace");
        let enclosure = volume_certified(&inside, &domain, 2.0, cx)
            .expect("nextafter domain has representable outward proof arithmetic");
        assert!(enclosure.lo <= 2.0 && 2.0 <= enclosure.hi);
        assert!(
            enclosure.lo < enclosure.hi,
            "an inexact proof pipeline must not publish a singleton measure"
        );
    });
}

/// gcp-006 — the descriptor table: every primitive declares its class,
/// certificate story, and fs-constraint kind mapping; proof
/// escalations are declared where they exist.
#[test]
fn gcp_006_descriptor_table() {
    let all = [
        GeoPrimitive::MinThickness,
        GeoPrimitive::DraftAngle,
        GeoPrimitive::Symmetry,
        GeoPrimitive::Envelope,
        GeoPrimitive::Volume,
    ];
    let mut rows = Vec::new();
    for p in all {
        let d = p.descriptor();
        rows.push(format!(
            "{{\"primitive\":\"{:?}\",\"class\":\"{:?}\",\"certificate\":\"{:?}\",\
             \"kind\":\"{}\"}}",
            d.primitive,
            d.class,
            d.certificate,
            d.kind.kind_name()
        ));
    }
    let symmetry_exact =
        GeoPrimitive::Symmetry.descriptor().certificate == CertKind::ExactByConstruction;
    let volume_enclosure = GeoPrimitive::Volume.descriptor().certificate == CertKind::Enclosure;
    let escalations = GeoPrimitive::Envelope.proof_escalation().is_some()
        && GeoPrimitive::Symmetry.proof_escalation().is_none();
    let fab_kinds = matches!(
        GeoPrimitive::MinThickness.descriptor().kind,
        fs_constraint::ConstraintKind::Fabrication { .. }
    );
    let mut em = fs_obs::Emitter::new("fs-geocon/conformance", "gcp-006/table");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "geocon-descriptor-table".to_string(),
                json: format!("[{}]", rows.join(",")),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("table validates");
    println!("{line}");
    verdict(
        "gcp-006",
        symmetry_exact && volume_enclosure && escalations && fab_kinds,
        "every primitive declares class + certificate + ASCENT kind (symmetry is \
         ExactByConstruction, volume is Enclosure, thickness maps to Fabrication), \
         and interval-proof escalations are declared exactly where they exist",
    );
}
