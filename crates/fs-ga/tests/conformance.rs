//! fs-ga conformance suite (the wqd.22 bead). Acceptance: the generated
//! tables pass the G0 identity battery; motor sandwiches match a matrix
//! reference to ULP-level bounds; exp/log round-trip including
//! near-degenerate screws; interpolation is gimbal-free; façade
//! conversions are exact; CGA constructions satisfy their incidence and
//! tangency laws.

use fs_ga::{
    Cga, Mat34, Motor, Pga, Plane, Point, Quat, Vec3, cga, exp_bivector, motor_log,
    pga::{Line, axis_bivector, ideal_bivector},
};
use fs_propcheck::Shrink;

#[derive(Clone, Debug)]
struct GaCoeffs<const N: usize>([f64; N]);

impl<const N: usize> Shrink for GaCoeffs<N> {
    fn shrink_candidates(&self) -> Vec<Self> {
        let mut candidates = Vec::new();
        for (index, value) in self.0.iter().enumerate() {
            for candidate in value.shrink_candidates() {
                let mut coefficients = self.0;
                coefficients[index] = candidate;
                candidates.push(GaCoeffs(coefficients));
            }
        }
        candidates
    }
}

fn generate_coefficients<const N: usize>(stream: &mut fs_propcheck::Stream) -> GaCoeffs<N> {
    GaCoeffs(core::array::from_fn(|_| stream.f64_in(-1.0, 1.0)))
}

fn as_pga(coefficients: &GaCoeffs<16>) -> Pga {
    Pga(coefficients.0)
}

fn as_cga(coefficients: &GaCoeffs<32>) -> Cga {
    Cga(coefficients.0)
}

fn pga_residual_within(residual: &Pga, tolerance: f64) -> bool {
    // max_abs deliberately ignores unordered NaN comparisons, so finiteness
    // is a separate, load-bearing part of every generated numeric law.
    residual.0.iter().all(|value| value.is_finite()) && residual.max_abs() < tolerance
}

fn cga_residual_within(residual: &Cga, tolerance: f64) -> bool {
    residual.0.iter().all(|value| value.is_finite()) && residual.max_abs() < tolerance
}

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-ga/conformance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64
}

fn rand_unit_axis(seed: &mut u64) -> Vec3 {
    loop {
        let v = Vec3::new(
            lcg(seed) * 2.0 - 1.0,
            lcg(seed) * 2.0 - 1.0,
            lcg(seed) * 2.0 - 1.0,
        );
        let n = v.norm();
        if n > 0.1 {
            return v.scale(1.0 / n);
        }
    }
}

fn rand_pga(seed: &mut u64) -> Pga {
    let mut v = Pga::zero();
    for c in &mut v.0 {
        *c = lcg(seed) * 2.0 - 1.0;
    }
    v
}

fn rand_cga(seed: &mut u64) -> Cga {
    let mut v = Cga::zero();
    for c in &mut v.0 {
        *c = lcg(seed) * 2.0 - 1.0;
    }
    v
}

/// Rodrigues rotation — the independent matrix-land reference.
fn rodrigues(axis: Vec3, angle: f64, v: Vec3) -> Vec3 {
    let (s, c) = (fs_math::det::sin(angle), fs_math::det::cos(angle));
    v.scale(c) + axis.cross(v).scale(s) + axis.scale(axis.dot(v) * (1.0 - c))
}

#[test]
fn ga_001_identity_battery() {
    let mut seed = 0x6A_0001u64;
    for round in 0..200 {
        let (a, b, c) = (
            rand_pga(&mut seed),
            rand_pga(&mut seed),
            rand_pga(&mut seed),
        );
        // Associativity and distributivity of the geometric product.
        let assoc = a.gp(&b).gp(&c).sub(&a.gp(&b.gp(&c))).max_abs();
        assert!(assoc < 1e-12, "PGA associativity {assoc} at round {round}");
        let distrib = a.gp(&b.add(&c)).sub(&a.gp(&b).add(&a.gp(&c))).max_abs();
        assert!(distrib < 1e-12, "PGA distributivity {distrib}");
        // Grade projections partition the multivector.
        let mut sum = Pga::zero();
        for g in 0..=4 {
            sum = sum.add(&a.grade_part(g));
        }
        assert_eq!(sum, a, "grade projections must partition");
        // Reverse anti-automorphism.
        let anti = a
            .gp(&b)
            .reverse()
            .sub(&b.reverse().gp(&a.reverse()))
            .max_abs();
        assert!(anti < 1e-12, "reverse anti-automorphism {anti}");
        // Wedge antisymmetry on vectors.
        let v1 = rand_pga(&mut seed).grade_part(1);
        let v2 = rand_pga(&mut seed).grade_part(1);
        let sym = v1.wedge(&v2).add(&v2.wedge(&v1)).max_abs();
        assert!(sym < 1e-12, "vector wedge antisymmetry {sym}");
        // Same laws in CGA.
        let (x, y, z) = (
            rand_cga(&mut seed),
            rand_cga(&mut seed),
            rand_cga(&mut seed),
        );
        let assoc = x.gp(&y).gp(&z).sub(&x.gp(&y.gp(&z))).max_abs();
        assert!(assoc < 1e-11, "CGA associativity {assoc}");
        let anti = x
            .gp(&y)
            .reverse()
            .sub(&y.reverse().gp(&x.reverse()))
            .max_abs();
        assert!(anti < 1e-11, "CGA reverse anti-automorphism {anti}");
    }
    verdict(
        "ga-001",
        "200 fixed rounds: five PGA laws; associativity and reverse in CGA",
    );
}

/// G0 property adoption (bead frankensim-4nh8): the five dense PGA laws with
/// integrated shrinking. The fixed ga-001 loop above remains unchanged.
#[test]
fn generated_pga_identity_laws() {
    fs_propcheck::check(
        "pga-geometric-product-associates",
        0x6A_1001,
        400,
        |s| {
            (
                generate_coefficients::<16>(s),
                generate_coefficients::<16>(s),
                generate_coefficients::<16>(s),
            )
        },
        |(a, b, c)| {
            let (a, b, c) = (as_pga(a), as_pga(b), as_pga(c));
            pga_residual_within(&a.gp(&b).gp(&c).sub(&a.gp(&b.gp(&c))), 1e-12)
        },
    );
    fs_propcheck::check(
        "pga-geometric-product-left-distributes",
        0x6A_1002,
        400,
        |s| {
            (
                generate_coefficients::<16>(s),
                generate_coefficients::<16>(s),
                generate_coefficients::<16>(s),
            )
        },
        |(a, b, c)| {
            let (a, b, c) = (as_pga(a), as_pga(b), as_pga(c));
            pga_residual_within(&a.gp(&b.add(&c)).sub(&a.gp(&b).add(&a.gp(&c))), 1e-12)
        },
    );
    fs_propcheck::check(
        "pga-grade-projections-partition",
        0x6A_1003,
        400,
        generate_coefficients::<16>,
        |coefficients| {
            let value = as_pga(coefficients);
            let reconstructed =
                (0..=4).fold(Pga::zero(), |sum, grade| sum.add(&value.grade_part(grade)));
            value.0.iter().all(|coefficient| coefficient.is_finite()) && reconstructed == value
        },
    );
    fs_propcheck::check(
        "pga-reverse-is-antiautomorphism",
        0x6A_1004,
        400,
        |s| {
            (
                generate_coefficients::<16>(s),
                generate_coefficients::<16>(s),
            )
        },
        |(a, b)| {
            let (a, b) = (as_pga(a), as_pga(b));
            pga_residual_within(
                &a.gp(&b).reverse().sub(&b.reverse().gp(&a.reverse())),
                1e-12,
            )
        },
    );
    fs_propcheck::check(
        "pga-vector-wedge-is-antisymmetric",
        0x6A_1005,
        400,
        |s| {
            (
                generate_coefficients::<16>(s),
                generate_coefficients::<16>(s),
            )
        },
        |(a, b)| {
            let (a, b) = (as_pga(a).grade_part(1), as_pga(b).grade_part(1));
            pga_residual_within(&a.wedge(&b).add(&b.wedge(&a)), 1e-12)
        },
    );
    verdict(
        "ga-g0-propcheck-pga",
        "five PGA laws x 400 generated cases, shrink-armed",
    );
}

/// The same complete five-law battery for CGA. This supplies the three CGA
/// laws the historical fixed ga-001 loop did not actually exercise.
#[test]
fn generated_cga_identity_laws() {
    fs_propcheck::check(
        "cga-geometric-product-associates",
        0x6A_2001,
        400,
        |s| {
            (
                generate_coefficients::<32>(s),
                generate_coefficients::<32>(s),
                generate_coefficients::<32>(s),
            )
        },
        |(a, b, c)| {
            let (a, b, c) = (as_cga(a), as_cga(b), as_cga(c));
            cga_residual_within(&a.gp(&b).gp(&c).sub(&a.gp(&b.gp(&c))), 1e-11)
        },
    );
    fs_propcheck::check(
        "cga-geometric-product-left-distributes",
        0x6A_2002,
        400,
        |s| {
            (
                generate_coefficients::<32>(s),
                generate_coefficients::<32>(s),
                generate_coefficients::<32>(s),
            )
        },
        |(a, b, c)| {
            let (a, b, c) = (as_cga(a), as_cga(b), as_cga(c));
            cga_residual_within(&a.gp(&b.add(&c)).sub(&a.gp(&b).add(&a.gp(&c))), 1e-11)
        },
    );
    fs_propcheck::check(
        "cga-grade-projections-partition",
        0x6A_2003,
        400,
        generate_coefficients::<32>,
        |coefficients| {
            let value = as_cga(coefficients);
            let reconstructed =
                (0..=5).fold(Cga::zero(), |sum, grade| sum.add(&value.grade_part(grade)));
            value.0.iter().all(|coefficient| coefficient.is_finite()) && reconstructed == value
        },
    );
    fs_propcheck::check(
        "cga-reverse-is-antiautomorphism",
        0x6A_2004,
        400,
        |s| {
            (
                generate_coefficients::<32>(s),
                generate_coefficients::<32>(s),
            )
        },
        |(a, b)| {
            let (a, b) = (as_cga(a), as_cga(b));
            cga_residual_within(
                &a.gp(&b).reverse().sub(&b.reverse().gp(&a.reverse())),
                1e-11,
            )
        },
    );
    fs_propcheck::check(
        "cga-vector-wedge-is-antisymmetric",
        0x6A_2005,
        400,
        |s| {
            (
                generate_coefficients::<32>(s),
                generate_coefficients::<32>(s),
            )
        },
        |(a, b)| {
            let (a, b) = (as_cga(a).grade_part(1), as_cga(b).grade_part(1));
            cga_residual_within(&a.wedge(&b).add(&b.wedge(&a)), 1e-11)
        },
    );
    verdict(
        "ga-g0-propcheck-cga",
        "five CGA laws x 400 generated cases, shrink-armed",
    );
}

#[test]
fn ga_002_motor_sandwich_matches_matrix_reference() {
    let mut seed = 0x6A_0002u64;
    let mut worst = 0.0f64;
    for _ in 0..300 {
        let axis = rand_unit_axis(&mut seed);
        let angle = (lcg(&mut seed) * 2.0 - 1.0) * std::f64::consts::PI;
        let t = Vec3::new(
            lcg(&mut seed) * 20.0 - 10.0,
            lcg(&mut seed) * 20.0 - 10.0,
            lcg(&mut seed) * 20.0 - 10.0,
        );
        let q = Quat::from_axis_angle(axis, angle);
        let motor = Motor::from_parts(q, t);
        for _ in 0..4 {
            let p = Vec3::new(
                lcg(&mut seed) * 8.0 - 4.0,
                lcg(&mut seed) * 8.0 - 4.0,
                lcg(&mut seed) * 8.0 - 4.0,
            );
            let expected = rodrigues(axis, angle, p) + t;
            // Quaternion façade agrees with Rodrigues.
            let via_quat = q.rotate(p) + t;
            let qerr = (via_quat - expected).norm();
            assert!(qerr < 1e-12, "quat vs Rodrigues {qerr}");
            // Motor sandwich agrees with the matrix reference.
            let got = motor
                .transform_point(Point::new(p.x, p.y, p.z))
                .expect("Euclidean in, Euclidean out");
            let err = (Vec3::new(got.x, got.y, got.z) - expected).norm();
            worst = worst.max(err);
            assert!(err < 1e-11, "motor sandwich diverged from matrix: {err}");
            // The monomorphized kernel agrees with the dense 16-component
            // reference path.
            let dense = motor
                .transform_point_dense(Point::new(p.x, p.y, p.z))
                .expect("dense path");
            let kd = Vec3::new(got.x - dense.x, got.y - dense.y, got.z - dense.z).norm();
            assert!(kd < 1e-12, "kernel vs dense sandwich {kd}");
            // Lowered Mat34 agrees too.
            let mat = Mat34::from_motor(&motor).expect("motor lowers");
            let merr = (mat.apply(p) - expected).norm();
            assert!(merr < 1e-11, "Mat34 path diverged: {merr}");
        }
        // Tiny-rotation stability (near screw-axis singularity).
        let eps_motor = Motor::rotor([axis.x, axis.y, axis.z], 1e-13);
        let p0 = Point::new(1.0, 2.0, 3.0);
        let moved = eps_motor.transform_point(p0).expect("finite");
        let drift = Vec3::new(moved.x - 1.0, moved.y - 2.0, moved.z - 3.0).norm();
        assert!(
            drift < 1e-11,
            "tiny rotation must be near-identity: {drift}"
        );
    }
    println!(
        "{{\"suite\":\"fs-ga/conformance\",\"metric\":\"sandwich_vs_matrix_worst\",\
         \"value\":{worst:.3e}}}"
    );
    verdict(
        "ga-002",
        "1200 sandwiches vs Rodrigues+translation reference",
    );
}

#[test]
fn ga_003_exp_log_round_trip_and_gimbal_free_interpolation() {
    let mut seed = 0x6A_0003u64;
    for _ in 0..300 {
        // Random screw: bounded angle keeps the principal branch.
        let axis = rand_unit_axis(&mut seed);
        let angle = (lcg(&mut seed) * 2.0 - 1.0) * 3.0;
        let b = axis_bivector(axis.x, axis.y, axis.z)
            .scale(angle / 2.0)
            .add(&ideal_bivector(
                lcg(&mut seed) * 4.0 - 2.0,
                lcg(&mut seed) * 4.0 - 2.0,
                lcg(&mut seed) * 4.0 - 2.0,
            ));
        let m = exp_bivector(&b);
        assert!(m.unit_defect() < 1e-12, "exp must produce unit motors");
        let b_back = motor_log(&m);
        let mut err = b_back.sub(&b).max_abs();
        // log returns the principal branch; fold the sign ambiguity.
        err = err.min(b_back.add(&b).max_abs());
        assert!(err < 1e-10, "exp/log round trip {err}");
        // Pure translators round-trip exactly through the θ→0 branch.
        let bt = ideal_bivector(lcg(&mut seed) * 9.0, -3.0, lcg(&mut seed));
        let tr_back = motor_log(&exp_bivector(&bt));
        assert!(tr_back.sub(&bt).max_abs() < 1e-14, "translator round trip");
    }
    // Gimbal fixture: interpolate between attitudes that Euler-angle lerp
    // mangles (90° pitch — the classic lock pose). The motor path must
    // stay a unit rigid motion with constant angular speed throughout.
    let m0 = Motor::from_parts(
        Quat::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), std::f64::consts::FRAC_PI_2),
        Vec3::new(0.0, 0.0, 0.0),
    );
    let m1 = Motor::from_parts(
        Quat::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), std::f64::consts::PI * 0.9),
        Vec3::new(3.0, -1.0, 2.0),
    );
    let total = motor_log(&m0.reverse().compose(&m1));
    let mut prev = m0;
    let steps = 32;
    let mut speeds = Vec::new();
    for k in 1..=steps {
        let t = f64::from(k) / f64::from(steps);
        let mk = m0.slerp(&m1, t);
        assert!(mk.unit_defect() < 1e-12, "slerp left the motor manifold");
        // Angular increment between consecutive frames.
        let step_log = motor_log(&prev.reverse().compose(&mk));
        speeds.push(step_log.max_abs());
        prev = mk;
    }
    let (min_s, max_s) = speeds.iter().fold((f64::INFINITY, 0.0f64), |(lo, hi), &s| {
        (lo.min(s), hi.max(s))
    });
    assert!(
        (max_s - min_s) / max_s < 1e-9,
        "screw interpolation must be constant-speed: {min_s} vs {max_s}"
    );
    // Endpoints are met exactly (to ULP noise).
    let end = m0.slerp(&m1, 1.0);
    assert!(end.0.sub(&m1.0).max_abs() < 1e-12, "slerp endpoint");
    let _ = total;
    verdict(
        "ga-003",
        "300 exp/log round trips; 32-step slerp through the gimbal pose: unit, constant-speed",
    );
}

#[test]
fn ga_004_incidence_join_meet() {
    let mut seed = 0x6A_0004u64;
    for _ in 0..200 {
        let p = Point::new(
            lcg(&mut seed) * 10.0 - 5.0,
            lcg(&mut seed) * 10.0 - 5.0,
            lcg(&mut seed) * 10.0 - 5.0,
        );
        let q = Point::new(
            lcg(&mut seed) * 10.0 - 5.0,
            lcg(&mut seed) * 10.0 - 5.0,
            lcg(&mut seed) * 10.0 - 5.0,
        );
        // Join: both endpoints and every affine combination lie on the line.
        let line = Line::through(p, q);
        let lam = lcg(&mut seed) * 2.0 - 0.5;
        let mid = Point::new(
            p.x + lam * (q.x - p.x),
            p.y + lam * (q.y - p.y),
            p.z + lam * (q.z - p.z),
        );
        for pt in [p, q, mid] {
            assert!(line.incidence(pt) < 1e-9, "join incidence broke");
        }
        // Plane incidence measure equals the implicit equation.
        let plane = Plane::new(
            lcg(&mut seed) * 2.0 - 1.0,
            lcg(&mut seed) * 2.0 - 1.0,
            lcg(&mut seed) * 2.0 - 1.0,
            lcg(&mut seed) * 2.0 - 1.0,
        );
        let algebraic = plane.a * p.x + plane.b * p.y + plane.c * p.z + plane.d;
        assert!(
            (plane.incidence(p) - algebraic).abs() < 1e-10,
            "plane incidence must equal the implicit equation"
        );
        // Meet of two planes: every point on the meet line satisfies both.
        let p2 = Plane::new(
            lcg(&mut seed) * 2.0 - 1.0,
            lcg(&mut seed) * 2.0 - 1.0,
            lcg(&mut seed) * 2.0 - 1.0,
            lcg(&mut seed) * 2.0 - 1.0,
        );
        let meet = Line::meet(plane, p2);
        // Recover two points on the meet line by joining with coordinate
        // planes is overkill; instead verify the line is contained in both
        // planes: plane ∧ line = 0.
        for pl in [plane, p2] {
            let contained = pl.to_mv().wedge(&meet.0).max_abs();
            let scale = meet.0.max_abs().max(1e-300);
            assert!(
                contained / scale < 1e-10,
                "meet line must lie in both planes"
            );
        }
    }
    verdict(
        "ga-004",
        "200 rounds: join/meet incidence + implicit-equation agreement",
    );
}

#[test]
fn ga_005_cga_rounds_and_tangency() {
    let mut seed = 0x6A_0005u64;
    for _ in 0..100 {
        // up/down round trip.
        let p = Vec3::new(
            lcg(&mut seed) * 10.0 - 5.0,
            lcg(&mut seed) * 10.0 - 5.0,
            lcg(&mut seed) * 10.0 - 5.0,
        );
        let back = cga::down(&cga::up(p)).expect("finite point");
        assert!((back - p).norm() < 1e-12, "up/down round trip");
        // up(p) is null and up(a)·up(b) = −½‖a−b‖².
        let a = cga::up(p);
        assert!(a.gp(&a).scalar_part().abs() < 1e-10, "up(p) must be null");
        let q = Vec3::new(
            lcg(&mut seed) * 6.0,
            lcg(&mut seed) * 6.0,
            lcg(&mut seed) * 6.0,
        );
        let dot = cga::up(p).gp(&cga::up(q)).scalar_part();
        let want = -0.5 * (p - q).dot(p - q);
        assert!((dot - want).abs() < 1e-9, "conformal distance law");
        // Sphere through four of its own points: recover center + radius.
        let center = Vec3::new(
            lcg(&mut seed) * 4.0 - 2.0,
            lcg(&mut seed) * 4.0 - 2.0,
            lcg(&mut seed) * 4.0 - 2.0,
        );
        let r = 0.5 + lcg(&mut seed) * 3.0;
        let on_sphere = |seed: &mut u64| {
            let d = rand_unit_axis(seed);
            center + d.scale(r)
        };
        let (s1, s2, s3, s4) = (
            on_sphere(&mut seed),
            on_sphere(&mut seed),
            on_sphere(&mut seed),
            on_sphere(&mut seed),
        );
        let sphere = cga::sphere_through(s1, s2, s3, s4);
        if sphere.max_abs() < 1e-6 {
            continue; // nearly-degenerate sample; the blade says so
        }
        let (c_got, r_got) = cga::sphere_center_radius(&sphere).expect("round sphere");
        assert!((c_got - center).norm() < 1e-6 * r.max(1.0), "sphere center");
        assert!((r_got - r).abs() < 1e-6 * r.max(1.0), "sphere radius");
        // A fifth point on the sphere is incident with the direct blade.
        let s5 = on_sphere(&mut seed);
        let inc = cga::incidence(s5, &sphere) / sphere.max_abs();
        assert!(inc < 1e-8, "fifth-point sphere incidence {inc}");
        // Circle through three points contains them.
        let circle = cga::circle_through(s1, s2, s3);
        for pt in [s1, s2, s3] {
            assert!(
                cga::incidence(pt, &circle) / circle.max_abs().max(1e-300) < 1e-9,
                "circle incidence"
            );
        }
        // Plane through three points contains an affine combination.
        let plane = cga::plane_through(s1, s2, s3);
        let comb = s1 + (s2 - s1).scale(0.3) + (s3 - s1).scale(0.4);
        assert!(
            cga::incidence(comb, &plane) / plane.max_abs().max(1e-300) < 1e-9,
            "plane incidence"
        );
        // Tangency: externally tangent spheres have zero residual;
        // clearly separated ones do not.
        let dir = rand_unit_axis(&mut seed);
        let r2 = 0.5 + lcg(&mut seed) * 2.0;
        let tangent_center = center + dir.scale(r + r2);
        let sa = cga::dual_sphere(center, r);
        let sb = cga::dual_sphere(tangent_center, r2);
        assert!(
            cga::tangency_residual(&sa, &sb).abs() < 1e-9,
            "tangent spheres must have zero residual"
        );
        let far = cga::dual_sphere(center + dir.scale(r + r2 + 1.5), r2);
        assert!(
            cga::tangency_residual(&sa, &far).abs() > 1e-3,
            "separated spheres must not read as tangent"
        );
    }
    verdict(
        "ga-005",
        "CGA: null embedding, distance law, sphere/circle/plane construction + tangency",
    );
}

#[test]
fn ga_006_facades_are_exact() {
    let mut seed = 0x6A_0006u64;
    for _ in 0..300 {
        let axis = rand_unit_axis(&mut seed);
        let angle = (lcg(&mut seed) * 2.0 - 1.0) * std::f64::consts::PI;
        let q = Quat::from_axis_angle(axis, angle);
        // Quat → rotor → quat is a bitwise relabeling.
        let q2 = Quat::from_rotor(&q.to_rotor());
        assert_eq!(q.w.to_bits(), q2.w.to_bits());
        assert_eq!(q.x.to_bits(), q2.x.to_bits());
        assert_eq!(q.y.to_bits(), q2.y.to_bits());
        assert_eq!(q.z.to_bits(), q2.z.to_bits());
        // Motor round trip: quat exact, translation to tight ULP.
        let t = Vec3::new(
            lcg(&mut seed) * 20.0 - 10.0,
            lcg(&mut seed) * 20.0 - 10.0,
            lcg(&mut seed) * 20.0 - 10.0,
        );
        let (q3, t3) = Motor::from_parts(q, t).to_parts();
        assert_eq!(q.w.to_bits(), q3.w.to_bits(), "rotor part must be exact");
        assert_eq!(q.x.to_bits(), q3.x.to_bits());
        assert!(
            (t3 - t).norm() < 1e-13 * (1.0 + t.norm()),
            "translation ULP"
        );
    }
    verdict(
        "ga-006",
        "300 rounds: quat relabel bitwise, motor→(q,t) exact/ULP",
    );
}

#[test]
fn ga_007_versor_drift_and_renormalization_policy() {
    let mut seed = 0x6A_0007u64;
    let mut m = Motor::identity();
    let mut max_drift: f64 = 0.0;
    for step in 0..20_000 {
        let axis = rand_unit_axis(&mut seed);
        let inc = Motor::from_parts(
            Quat::from_axis_angle(axis, (lcg(&mut seed) - 0.5) * 0.02),
            Vec3::new(
                (lcg(&mut seed) - 0.5) * 0.01,
                (lcg(&mut seed) - 0.5) * 0.01,
                (lcg(&mut seed) - 0.5) * 0.01,
            ),
        );
        m = m.compose(&inc);
        if step % 64 == 63 {
            max_drift = max_drift.max(m.renormalize());
        }
    }
    let final_defect = m.unit_defect();
    println!(
        "{{\"suite\":\"fs-ga/conformance\",\"metric\":\"versor_drift\",\
         \"max_between_renorms\":{max_drift:.3e},\"final_defect\":{final_defect:.3e}}}"
    );
    assert!(
        max_drift < 1e-11,
        "drift between renormalizations exceeded policy bound: {max_drift}"
    );
    assert!(final_defect < 1e-12, "final motor not unit: {final_defect}");
    verdict(
        "ga-007",
        "20k-product chain, renorm every 64: drift within declared bounds",
    );
}
