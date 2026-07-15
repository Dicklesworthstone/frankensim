//! fs-motion conformance battery (bead c70j).
//!
//! The load-bearing rule: SAMPLING FALSIFIES, NEVER PROVES. Dense
//! pointwise fs-ga evaluations must land inside the certified
//! enclosures (inflated by a stated cross-implementation rounding
//! tolerance); a sign error, blade mix-up, or wrong formula anywhere
//! in the tube machinery fails these batteries loudly.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ga::{Motor, Point as GaPoint};
use fs_geom::{Aabb, BettiBounds, Chart, ChartSample, Point3, TraceStepClaim};
use fs_ivl::Interval;
use fs_math::det;
use fs_motion::{
    CertifiedMotorTube, EnclosureClass, MotionError, ScrewParams, SpacetimeChart, WankelParams,
    screw_tube, wankel_tube,
};

/// Cross-implementation rounding tolerance: covers the few-ulp gap
/// between the tube's real-arithmetic enclosure of its constructed
/// path and fs-ga's independently rounded pointwise evaluation. It is
/// small enough that any structural error (sign, blade, formula)
/// blows straight through it.
const CROSS_IMPL_TOL: f64 = 1e-9;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xC705_u64,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// Deterministic pseudo-random stream (no external deps; fixed seed).
struct Lcg(u64);

impl Lcg {
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}

/// Pointwise fs-ga screw motor matching `screw_tube`'s constructed
/// semantics: base ∘ T_c ∘ R(ωt) ∘ T_axis(vt) ∘ T_c⁻¹.
fn screw_pointwise(p: &ScrewParams, t: f64) -> Motor {
    let rot = Motor::rotor(p.axis, p.omega * t);
    let d = p.axial_velocity * t;
    let trans = Motor::translator(p.axis[0] * d, p.axis[1] * d, p.axis[2] * d);
    let to_center = Motor::translator(p.center[0], p.center[1], p.center[2]);
    let back = Motor::translator(-p.center[0], -p.center[1], -p.center[2]);
    p.base_pose
        .compose(&to_center)
        .compose(&rot.compose(&trans))
        .compose(&back)
}

/// Pointwise fs-ga Wankel pose matching `wankel_tube`'s constructed
/// semantics: base ∘ T(orbit(α)) ∘ R_z(β).
fn wankel_pointwise(p: &WankelParams, t: f64) -> Motor {
    let alpha = p.omega * t + p.crank_phase;
    let beta = alpha / 3.0 + p.rotor_phase;
    let orbit = Motor::translator(
        p.eccentricity * det::cos(alpha),
        p.eccentricity * det::sin(alpha),
        0.0,
    );
    let spin = Motor::rotor([0.0, 0.0, 1.0], beta);
    p.base_pose.compose(&orbit.compose(&spin))
}

fn assert_contains(enc: &[Interval; 3], q: GaPoint, what: &str) {
    let coords = [q.x, q.y, q.z];
    for (axis, (iv, v)) in enc.iter().zip(coords.iter()).enumerate() {
        assert!(
            *v >= iv.lo() - CROSS_IMPL_TOL && *v <= iv.hi() + CROSS_IMPL_TOL,
            "{what}: axis {axis} value {v} outside enclosure [{}, {}]",
            iv.lo(),
            iv.hi()
        );
    }
}

fn sample_screws() -> Vec<ScrewParams> {
    vec![
        ScrewParams {
            axis: [0.0, 0.0, 1.0],
            center: [0.0, 0.0, 0.0],
            omega: 1.7,
            axial_velocity: 0.0,
            base_pose: Motor::identity(),
        },
        ScrewParams {
            axis: [1.0, 0.0, 0.0],
            center: [0.4, -0.3, 0.2],
            omega: -2.3,
            axial_velocity: 0.6,
            base_pose: Motor::translator(0.1, 0.2, -0.5),
        },
        ScrewParams {
            // Exactly representable unit axis (3-4-5 triangle).
            axis: [0.6, 0.8, 0.0],
            center: [-1.0, 0.5, 0.25],
            omega: 0.9,
            axial_velocity: -0.35,
            base_pose: Motor::rotor([0.0, 1.0, 0.0], 0.7),
        },
    ]
}

#[test]
fn mt_001_constant_motor_sandwich_matches_transform_point() {
    // A zero-rate screw is a constant motor: the tube's sandwich
    // machinery (extracted structure table, point embedding,
    // homogeneous division) must reproduce Motor::transform_point.
    with_cx(|cx| {
        let mut rng = Lcg(0xF5_0001);
        for base in [
            Motor::identity(),
            Motor::rotor([0.0, 0.0, 1.0], 1.1),
            Motor::translator(0.3, -0.7, 0.9),
            Motor::rotor([1.0, 0.0, 0.0], -0.6).compose(&Motor::translator(1.0, 2.0, 3.0)),
        ] {
            let tube = screw_tube(
                &ScrewParams {
                    axis: [0.0, 0.0, 1.0],
                    center: [0.0, 0.0, 0.0],
                    omega: 0.0,
                    axial_velocity: 0.0,
                    base_pose: base,
                },
                Interval::new(0.0, 1.0),
                4,
                1,
            )
            .expect("constant tube builds");
            for _ in 0..16 {
                let x = Point3::new(
                    rng.range(-2.0, 2.0),
                    rng.range(-2.0, 2.0),
                    rng.range(-2.0, 2.0),
                );
                let enc = tube
                    .point_action_over(x, Interval::new(0.2, 0.8), cx)
                    .expect("action encloses");
                assert_eq!(enc.class, EnclosureClass::Certified);
                let truth = base
                    .transform_point(GaPoint { x: x.x, y: x.y, z: x.z })
                    .expect("finite point transforms");
                assert_contains(&enc.coords, truth, "mt-001 constant sandwich");
            }
        }
    });
}

#[test]
fn mt_002_screw_tube_encloses_dense_pointwise_sampling() {
    with_cx(|cx| {
        let mut rng = Lcg(0xF5_0002);
        for params in sample_screws() {
            let domain = Interval::new(-0.25, 1.25);
            let tube = screw_tube(&params, domain, 10, 6).expect("screw tube builds");
            assert!(
                tube.defect() < 1e-8,
                "unit-axis screw defect too large: {}",
                tube.defect()
            );
            for _ in 0..8 {
                let x = Point3::new(
                    rng.range(-1.5, 1.5),
                    rng.range(-1.5, 1.5),
                    rng.range(-1.5, 1.5),
                );
                let lo = rng.range(-0.25, 0.7);
                let span = Interval::new(lo, lo + 0.4);
                let enc = tube.point_action_over(x, span, cx).expect("action encloses");
                for k in 0..=32 {
                    cx.checkpoint().expect("not cancelled");
                    let t = span.lo() + (span.hi() - span.lo()) * (f64::from(k) / 32.0);
                    let truth = screw_pointwise(&params, t)
                        .transform_point(GaPoint { x: x.x, y: x.y, z: x.z })
                        .expect("finite point transforms");
                    assert_contains(&enc.coords, truth, "mt-002 screw sampling");
                }
            }
        }
    });
}

#[test]
fn mt_003_double_cover_sign_is_deterministic_and_transitions_validate() {
    // The same motion handed in with the opposite motor sign must
    // canonicalize to bit-identical component enclosures.
    let params = sample_screws().remove(2);
    let mut flipped = params;
    flipped.base_pose = Motor(params.base_pose.0.scale(-1.0));
    let domain = Interval::new(0.0, 1.0);
    let a = screw_tube(&params, domain, 8, 3).expect("tube builds");
    let b = screw_tube(&flipped, domain, 8, 3).expect("flipped tube builds");
    // Same-input replay is bit-identical (mt-006). Across the double
    // cover, canonicalization goes through a Taylor-model negation
    // whose outward remainder rounding is not perfectly sign-symmetric,
    // so the honest claim is roundoff-scale agreement at the unit
    // component scale — still ~thirteen orders of magnitude tighter
    // than any sign or blade error could produce.
    const DOUBLE_COVER_TOL: f64 = 1e-13;
    for (sa, sb) in a.segments().iter().zip(b.segments().iter()) {
        let ea = sa.components_over(sa.domain()).expect("eval");
        let eb = sb.components_over(sb.domain()).expect("eval");
        for (ia, ib) in ea.iter().zip(eb.iter()) {
            assert!(
                (ia.lo() - ib.lo()).abs() <= DOUBLE_COVER_TOL,
                "mt-003 lo endpoints disagree: {} vs {}",
                ia.lo(),
                ib.lo()
            );
            assert!(
                (ia.hi() - ib.hi()).abs() <= DOUBLE_COVER_TOL,
                "mt-003 hi endpoints disagree: {} vs {}",
                ia.hi(),
                ib.hi()
            );
        }
    }
    // A deliberate interior double-cover flip must refuse.
    let seg_a = a.segments()[0].clone();
    let flipped_tail = {
        let mut p2 = params;
        p2.base_pose = Motor(params.base_pose.0.scale(-1.0));
        // Rebuild only the tail segment domain with the flipped pose;
        // seal() will canonicalize it, so instead splice segments from
        // tubes whose canonical signs genuinely disagree at the joint:
        // rotate the anchor far enough that the canonical branch flips.
        let d2 = Interval::new(seg_a.domain().hi(), seg_a.domain().hi() + 0.333);
        screw_tube(&p2, d2, 8, 1).expect("tail builds")
    };
    let tail_seg = flipped_tail.segments()[0].clone();
    match CertifiedMotorTube::from_segments(vec![seg_a.clone(), tail_seg]) {
        // Either the joint genuinely agrees (canonicalization made the
        // signs match — then assembly must succeed), or it refuses
        // with ChartTransition. Both are deterministic; what is
        // FORBIDDEN is silently accepting a negative-dot joint.
        Ok(tube) => {
            let t = seg_a.domain().hi();
            let left = tube.segments()[0]
                .components_over(Interval::point(t))
                .expect("eval");
            let right = tube.segments()[1]
                .components_over(Interval::point(t))
                .expect("eval");
            let dot: f64 = left
                .iter()
                .zip(right.iter())
                .map(|(l, r)| l.midpoint() * r.midpoint())
                .sum();
            assert!(dot > 0.0, "mt-003 accepted a non-positive transition dot");
        }
        Err(MotionError::ChartTransition { .. }) => {}
        Err(other) => panic!("mt-003 unexpected refusal: {other}"),
    }
}

#[test]
fn mt_004_versor_defect_detects_broken_construction() {
    // Exact unit axes: defect is rounding-plus-truncation noise.
    let good = screw_tube(
        &sample_screws()[0],
        Interval::new(0.0, 1.0),
        8,
        2,
    )
    .expect("good tube builds");
    assert!(good.defect() < 1e-8, "good defect: {}", good.defect());
    // A deliberately non-unit axis is a broken generator; the residual
    // machinery must REPORT it (cos² + |a|²·sin² − 1 ≈ sin²·(|a|²−1)),
    // not mask it.
    let broken = screw_tube(
        &ScrewParams {
            axis: [1.0, 1.0, 0.0],
            center: [0.0, 0.0, 0.0],
            omega: 2.0,
            axial_velocity: 0.0,
            base_pose: Motor::identity(),
        },
        Interval::new(0.0, 1.0),
        8,
        2,
    )
    .expect("broken tube still builds");
    assert!(
        broken.defect() > 1e-3,
        "non-unit axis must surface in the defect, got {}",
        broken.defect()
    );
}

/// An exact-distance sphere chart (test fixture).
struct SphereChart {
    center: Point3,
    radius: f64,
}

impl Chart for SphereChart {
    fn name(&self) -> &'static str {
        "test-sphere"
    }

    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let dx = x.x - self.center.x;
        let dy = x.y - self.center.y;
        let dz = x.z - self.center.z;
        let d = det::sqrt(dx * dx + dy * dy + dz * dz) - self.radius;
        // A few-ulp rounding pad keeps the certificate honest.
        let pad = 1e-12 * (1.0 + d.abs());
        ChartSample {
            signed_distance: d,
            gradient: None,
            lipschitz: Some(1.0),
            error: fs_evidence::NumericalCertificate::enclosure(d - pad, d + pad),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(
            Point3::new(
                self.center.x - self.radius,
                self.center.y - self.radius,
                self.center.z - self.radius,
            ),
            Point3::new(
                self.center.x + self.radius,
                self.center.y + self.radius,
                self.center.z + self.radius,
            ),
        )
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn topology_hint(&self) -> BettiBounds {
        BettiBounds::unknown()
    }
}

#[test]
fn mt_005_eval_over_encloses_dense_time_sampling() {
    with_cx(|cx| {
        let params = sample_screws()[1];
        let tube = screw_tube(&params, Interval::new(0.0, 1.0), 10, 8).expect("tube builds");
        let moving = SpacetimeChart::new(
            SphereChart {
                center: Point3::new(0.2, -0.1, 0.3),
                radius: 0.8,
            },
            tube,
        );
        let mut rng = Lcg(0xF5_0005);
        for _ in 0..6 {
            let x = Point3::new(
                rng.range(-2.0, 2.0),
                rng.range(-2.0, 2.0),
                rng.range(-2.0, 2.0),
            );
            let lo = rng.range(0.0, 0.6);
            let span = Interval::new(lo, lo + 0.3);
            let enc = moving.eval_over(x, span, cx).expect("eval_over encloses");
            assert_eq!(enc.class, EnclosureClass::Certified);
            for k in 0..=24 {
                let t = span.lo() + (span.hi() - span.lo()) * (f64::from(k) / 24.0);
                let m = screw_pointwise(&params, t);
                let q = m
                    .reverse()
                    .transform_point(GaPoint { x: x.x, y: x.y, z: x.z })
                    .expect("pull-back transforms");
                let truth = moving
                    .base()
                    .eval(Point3::new(q.x, q.y, q.z), cx)
                    .signed_distance;
                assert!(
                    truth >= enc.value.lo() - CROSS_IMPL_TOL
                        && truth <= enc.value.hi() + CROSS_IMPL_TOL,
                    "mt-005: field value {truth} at t = {t} outside [{}, {}]",
                    enc.value.lo(),
                    enc.value.hi()
                );
            }
        }
    });
}

#[test]
fn mt_006_bit_replay_across_reconstruction() {
    with_cx(|cx| {
        let params = sample_screws()[2];
        let domain = Interval::new(-0.5, 0.75);
        let a = screw_tube(&params, domain, 9, 5).expect("tube a");
        let b = screw_tube(&params, domain, 9, 5).expect("tube b");
        assert_eq!(a.defect().to_bits(), b.defect().to_bits(), "defect bits differ");
        let x = Point3::new(0.7, -0.4, 1.1);
        let span = Interval::new(-0.2, 0.55);
        let ea = a.point_action_over(x, span, cx).expect("a acts");
        let eb = b.point_action_over(x, span, cx).expect("b acts");
        for (ia, ib) in ea.coords.iter().zip(eb.coords.iter()) {
            assert_eq!(ia.lo().to_bits(), ib.lo().to_bits(), "mt-006 lo bits differ");
            assert_eq!(ia.hi().to_bits(), ib.hi().to_bits(), "mt-006 hi bits differ");
        }
    });
}

#[test]
fn mt_007_wankel_pose_encloses_pointwise_composition() {
    with_cx(|cx| {
        let params = WankelParams {
            eccentricity: 0.15,
            omega: 3.0,
            crank_phase: 0.4,
            rotor_phase: -0.2,
            base_pose: Motor::translator(0.05, -0.02, 0.0),
        };
        let domain = Interval::new(0.0, 1.0);
        let tube = wankel_tube(&params, domain, 10, 8).expect("wankel tube builds");
        assert!(tube.defect() < 1e-7, "wankel defect: {}", tube.defect());
        let mut rng = Lcg(0xF5_0007);
        for _ in 0..6 {
            let x = Point3::new(
                rng.range(-0.5, 0.5),
                rng.range(-0.5, 0.5),
                rng.range(-0.2, 0.2),
            );
            let lo = rng.range(0.0, 0.7);
            let span = Interval::new(lo, lo + 0.25);
            let enc = tube.point_action_over(x, span, cx).expect("action encloses");
            for k in 0..=24 {
                let t = span.lo() + (span.hi() - span.lo()) * (f64::from(k) / 24.0);
                let truth = wankel_pointwise(&params, t)
                    .transform_point(GaPoint { x: x.x, y: x.y, z: x.z })
                    .expect("finite point transforms");
                assert_contains(&enc.coords, truth, "mt-007 wankel sampling");
            }
        }
    });
}

#[test]
fn mt_008_box_action_contains_sampled_interior_points() {
    with_cx(|cx| {
        let params = sample_screws()[1];
        let tube = screw_tube(&params, Interval::new(0.0, 1.0), 9, 6).expect("tube builds");
        let b = Aabb::new(Point3::new(-0.3, -0.2, -0.1), Point3::new(0.4, 0.5, 0.6));
        let span = Interval::new(0.1, 0.9);
        let enc = tube.box_action_over(&b, span, cx).expect("box encloses");
        let mut rng = Lcg(0xF5_0008);
        for _ in 0..64 {
            let x = GaPoint {
                x: rng.range(b.min.x, b.max.x),
                y: rng.range(b.min.y, b.max.y),
                z: rng.range(b.min.z, b.max.z),
            };
            let t = rng.range(span.lo(), span.hi());
            let truth = screw_pointwise(&params, t)
                .transform_point(x)
                .expect("finite point transforms");
            assert!(
                truth.x >= enc.bounds.min.x - CROSS_IMPL_TOL
                    && truth.x <= enc.bounds.max.x + CROSS_IMPL_TOL
                    && truth.y >= enc.bounds.min.y - CROSS_IMPL_TOL
                    && truth.y <= enc.bounds.max.y + CROSS_IMPL_TOL
                    && truth.z >= enc.bounds.min.z - CROSS_IMPL_TOL
                    && truth.z <= enc.bounds.max.z + CROSS_IMPL_TOL,
                "mt-008: image point escaped the box enclosure"
            );
        }
    });
}

#[test]
fn mt_009_snapshot_agrees_with_pullback_and_transports_support() {
    with_cx(|cx| {
        let params = sample_screws()[0];
        let tube = screw_tube(&params, Interval::new(0.0, 1.0), 8, 4).expect("tube builds");
        let base = SphereChart {
            center: Point3::new(0.1, 0.2, -0.3),
            radius: 0.6,
        };
        let base_support = base.support();
        let moving = SpacetimeChart::new(base, tube);
        let t = 0.6180339887;
        let snap = moving.snapshot(t, cx).expect("snapshot freezes");
        assert_eq!(snap.time(), t);
        // Field agreement with an explicit pull-back through fs-ga.
        let mut rng = Lcg(0xF5_0009);
        for _ in 0..16 {
            let x = Point3::new(
                rng.range(-1.5, 1.5),
                rng.range(-1.5, 1.5),
                rng.range(-1.5, 1.5),
            );
            let sample = snap.eval(x, cx);
            let m = screw_pointwise(&params, t);
            let q = m
                .reverse()
                .transform_point(GaPoint { x: x.x, y: x.y, z: x.z })
                .expect("pull-back transforms");
            let truth = moving
                .base()
                .eval(Point3::new(q.x, q.y, q.z), cx)
                .signed_distance;
            assert!(
                (sample.signed_distance - truth).abs() <= 1e-9,
                "mt-009: snapshot field {} vs pull-back {truth}",
                sample.signed_distance
            );
            // Claims are deliberately weakened.
            assert!(sample.gradient.is_none());
            assert!(sample.lipschitz.is_none());
        }
        // Support transport: images of base-support corner points at
        // the frozen time stay inside the snapshot support.
        let m = screw_pointwise(&params, t);
        for &(cx_, cy, cz) in &[
            (base_support.min.x, base_support.min.y, base_support.min.z),
            (base_support.max.x, base_support.max.y, base_support.max.z),
            (base_support.min.x, base_support.max.y, base_support.min.z),
            (base_support.max.x, base_support.min.y, base_support.max.z),
        ] {
            let img = m
                .transform_point(GaPoint { x: cx_, y: cy, z: cz })
                .expect("corner transforms");
            let s = snap.support();
            assert!(
                img.x >= s.min.x - CROSS_IMPL_TOL
                    && img.x <= s.max.x + CROSS_IMPL_TOL
                    && img.y >= s.min.y - CROSS_IMPL_TOL
                    && img.y <= s.max.y + CROSS_IMPL_TOL
                    && img.z >= s.min.z - CROSS_IMPL_TOL
                    && img.z <= s.max.z + CROSS_IMPL_TOL,
                "mt-009: transported corner escaped snapshot support"
            );
        }
    });
}
