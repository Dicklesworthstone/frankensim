//! fs-query conformance suite (CONTRACT.md: any reimplementation must
//! pass). Multi-chart AGREEMENT for closest point and raycasts,
//! tracer safety vs a dense oracle including tangent rays, offsets and
//! ball-Minkowski exactness, certified separation bounds, the
//! thickness oracle on graded fixtures with the medial cross-check and
//! the design-lever subgradient, and curvature convergence at the
//! documented order per chart class. JSON-line verdicts; seeded cases
//! carry seeds.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Chart, Point3, Vec3};
use fs_query::{
    CurvatureClass, OffsetChart, QueryError, closest_point, curvature, medial_poles, min_thickness,
    minkowski_ball, raycast, separation, thickness_at,
};
use fs_rep_frep::{BoolOp, BoolStyle, FrepBuilder};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-query/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\
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
                seed: 0x9E4,
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

fn frep_sphere(center: Point3, r: f64) -> fs_rep_frep::Frep {
    let mut b = FrepBuilder::new();
    let s = b.sphere(center, r).expect("sphere");
    b.finish(s).expect("frep")
}

/// gq-001 — closest point agrees across chart representations of the
/// SAME sphere within per-chart certificates, residuals are honest,
/// and the answer is translation-equivariant (G3).
#[test]
fn gq_001_closest_point_agreement() {
    with_cx(|cx| {
        let exact = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let frep = frep_sphere(Point3::new(0.0, 0.0, 0.0), 1.0);
        let tiled = fs_rep_sdf::TiledSdf::build(&exact, 0.05, cx).expect("tiled");
        let mesh = fs_rep_mesh::MeshChart::new(fs_rep_mesh::shapes::icosphere(
            Point3::new(0.0, 0.0, 0.0),
            1.0,
            4,
        ));
        // Per-chart tolerance = each chart's OWN certificate: exact and
        // F-rep are exact; the tiled grid declares its interpolation
        // bound; the mesh chart agrees at its faceting scale.
        let charts: Vec<(&dyn Chart, f64)> = vec![
            (&exact, 1e-9),
            (&frep, 1e-9),
            (&tiled, tiled.bound() * 1.05),
            (&mesh, 3e-2),
        ];
        let mut rng = Lcg(0x1001_2026_0707_0011);
        let mut worst = vec![0.0f64; charts.len()];
        for _ in 0..80 {
            // Shell radii keep queries inside every chart's valid
            // domain (the tiled grid ends at its support box).
            let z = rng.range(-1.0, 1.0);
            let sq = (1.0f64 - z * z).max(0.0).sqrt();
            let th = rng.range(0.0, core::f64::consts::TAU);
            let rad = rng.range(0.3, 1.35);
            let p = Point3::new(rad * sq * th.cos(), rad * sq * th.sin(), rad * z);
            let truth = {
                let d = p.delta_from(Point3::new(0.0, 0.0, 0.0));
                let n = d.norm();
                Point3::new(d.x / n, d.y / n, d.z / n)
            };
            for (k, (chart, _)) in charts.iter().enumerate() {
                let cp = closest_point(*chart, p, cx).expect("closest");
                worst[k] = worst[k]
                    .max(cp.point.delta_from(truth).norm())
                    .max(cp.residual);
            }
        }
        let agree = charts
            .iter()
            .enumerate()
            .all(|(k, (_, tol))| worst[k] < *tol);
        // G3: shifted sphere, shifted queries, identical geometry.
        let shifted = SphereChart {
            center: Point3::new(0.5, 0.25, -0.375),
            radius: 1.0,
        };
        let q = Point3::new(1.7, 0.3, -0.2);
        let base = closest_point(&exact, q, cx).expect("base");
        let moved =
            closest_point(&shifted, q.offset(Vec3::new(0.5, 0.25, -0.375)), cx).expect("moved");
        let g3 = moved
            .point
            .delta_from(base.point.offset(Vec3::new(0.5, 0.25, -0.375)))
            .norm()
            < 1e-9;
        verdict(
            "gq-001",
            agree && g3,
            &format!(
                "closest point agrees with the analytic answer across exact/F-rep/\
                 tiled-SDF/mesh charts (worst errors {:?} within per-chart \
                 certificates) with honest residuals, and is translation-equivariant; \
                 seed 0x1001_2026_0707_0011",
                worst.iter().map(|w| format!("{w:.1e}")).collect::<Vec<_>>()
            ),
        );
    });
}

/// gq-002 — raycast: analytic agreement across chart types, and SAFETY
/// on CSG (no tunneling vs a dense oracle), including tangent rays.
#[test]
fn gq_002_raycast_safety() {
    with_cx(|cx| {
        let exact = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let frep = frep_sphere(Point3::new(0.0, 0.0, 0.0), 1.0);
        let charts: Vec<(&dyn Chart, f64)> = vec![(&exact, 1e-6), (&frep, 1e-6)];
        let mut rng = Lcg(0x1001_2026_0707_0012);
        let mut agree = true;
        for _ in 0..100 {
            let z = rng.range(-1.0, 1.0);
            let s = (1.0 - z * z).max(0.0).sqrt();
            let th = rng.range(0.0, core::f64::consts::TAU);
            let o = Point3::new(3.0 * s * th.cos(), 3.0 * s * th.sin(), 3.0 * z);
            let dir = Point3::new(0.0, 0.0, 0.0).delta_from(o);
            // Analytic first hit: |o| − 1 along the unit direction.
            let expect_t = o.delta_from(Point3::new(0.0, 0.0, 0.0)).norm() - 1.0;
            for (chart, tol) in &charts {
                let hit = raycast(*chart, o, dir, 10.0, cx)
                    .expect("cast")
                    .expect("hits");
                agree &= (hit.t - expect_t).abs() < 1e-3 + tol;
            }
        }
        // Tangent rays: aimed exactly at radius 1 (grazing) and 1.05
        // (clean miss): no tunneling, misses classified.
        let graze_o = Point3::new(-3.0, 1.0, 0.0);
        let graze = raycast(&frep, graze_o, Vec3::new(1.0, 0.0, 0.0), 10.0, cx).expect("graze");
        let graze_ok = match graze {
            Some(h) => (h.point.delta_from(Point3::new(0.0, 0.0, 0.0)).norm() - 1.0).abs() < 1e-3,
            None => true, // grazing may legitimately exhaust its budget
        };
        let miss = raycast(
            &frep,
            Point3::new(-3.0, 1.05, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            10.0,
            cx,
        )
        .expect("miss");
        // CSG safety vs dense oracle.
        let csg = {
            let mut b = FrepBuilder::new();
            let s1 = b.sphere(Point3::new(-0.4, 0.0, 0.0), 0.8).expect("s1");
            let s2 = b.sphere(Point3::new(0.4, 0.0, 0.0), 0.8).expect("s2");
            let u = b
                .boolean(BoolOp::Difference, BoolStyle::Blend { radius: 0.2 }, s1, s2)
                .expect("u");
            b.finish(u).expect("frep")
        };
        let mut safety = true;
        for _ in 0..200 {
            let z = rng.range(-1.0, 1.0);
            let s = (1.0 - z * z).max(0.0).sqrt();
            let th = rng.range(0.0, core::f64::consts::TAU);
            let o = Point3::new(2.5 * s * th.cos(), 2.5 * s * th.sin(), 2.5 * z);
            let target = Point3::new(
                rng.range(-0.5, 0.5),
                rng.range(-0.5, 0.5),
                rng.range(-0.5, 0.5),
            );
            let dir = target.delta_from(o);
            let dn = dir.norm();
            let d = dir.scale(1.0 / dn);
            let hit = raycast(&csg, o, d, 6.0, cx).expect("cast");
            // Dense oracle: first sign change.
            let mut oracle = None;
            let mut prev = csg.value(o);
            for i in 1..=1200 {
                let t = 6.0 * f64::from(i) / 1200.0;
                let v = csg.value(o.offset(d.scale(t)));
                if prev >= 0.0 && v < 0.0 {
                    oracle = Some(t);
                    break;
                }
                prev = v;
            }
            if let (Some(t_true), Some(h)) = (oracle, hit) {
                safety &= h.t <= t_true + 1e-3;
            }
            // Oracle-hit + tracer-miss = grazing budget exhaustion:
            // incomplete, not unsafe. Oracle-miss cases carry no claim.
        }
        verdict(
            "gq-002",
            agree && graze_ok && miss.is_none() && safety,
            "raycasts match the analytic sphere across chart types, tangent rays \
             never tunnel (grazes land on the surface or approach; the 1.05 offset \
             misses cleanly), and the CSG tracer never claims a hit past the dense \
             oracle over 200 rays; seed 0x1001_2026_0707_0012",
        );
    });
}

/// gq-003 — offsets and the ball-Minkowski identity: offset spheres
/// are spheres (across chart types), erosion shrinks, and
/// minkowski_ball IS the offset (exact by construction).
#[test]
fn gq_003_offset_minkowski() {
    with_cx(|cx| {
        let exact = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let frep = frep_sphere(Point3::new(0.0, 0.0, 0.0), 1.0);
        let mut rng = Lcg(0x1001_2026_0707_0013);
        let mut ok = true;
        for (chart, tol) in [(&exact as &dyn Chart, 1e-12), (&frep as &dyn Chart, 1e-12)] {
            let grown = OffsetChart::new(chart, 0.3);
            let eroded = OffsetChart::new(chart, -0.2);
            let mink = minkowski_ball(chart, 0.3);
            for _ in 0..100 {
                let p = Point3::new(
                    rng.range(-2.0, 2.0),
                    rng.range(-2.0, 2.0),
                    rng.range(-2.0, 2.0),
                );
                let d = p.delta_from(Point3::new(0.0, 0.0, 0.0)).norm();
                ok &= (grown.eval(p, cx).signed_distance - (d - 1.3)).abs() < tol;
                ok &= (eroded.eval(p, cx).signed_distance - (d - 0.8)).abs() < tol;
                ok &= grown.eval(p, cx).signed_distance.to_bits()
                    == mink.eval(p, cx).signed_distance.to_bits();
            }
        }
        // Offset charts remain queryable: closest point on the grown
        // sphere lands at radius 1.3.
        let grown = OffsetChart::new(&exact, 0.3);
        let cp = closest_point(&grown, Point3::new(2.0, 0.5, -0.3), cx).expect("cp");
        let on_grown = (cp.point.delta_from(Point3::new(0.0, 0.0, 0.0)).norm() - 1.3).abs() < 1e-9;
        verdict(
            "gq-003",
            ok && on_grown,
            "offset spheres are exactly spheres of the summed radius across chart \
             types, erosion shrinks exactly, minkowski_ball is BITWISE the offset \
             chart, and offset charts stay fully queryable; \
             seed 0x1001_2026_0707_0013",
        );
    });
}

/// gq-004 — certified separation: the rigorous lower bound brackets
/// the analytic separation of two spheres, tracks as they approach,
/// and the clearance field dominates the separation everywhere (G0).
#[test]
fn gq_004_separation_certified() {
    with_cx(|cx| {
        let mut rng = Lcg(0x1001_2026_0707_0014);
        let mut ok = true;
        let mut worst_gap = 0.0f64;
        for gap in [1.0, 0.5, 0.2, 0.05] {
            let a = SphereChart {
                center: Point3::new(-(1.0 + gap / 2.0), 0.0, 0.0),
                radius: 1.0,
            };
            let b = SphereChart {
                center: Point3::new(1.0 + gap / 2.0, 0.0, 0.0),
                radius: 1.0,
            };
            let sep = separation(&a, &b, 24, cx).expect("separation");
            ok &= sep.lower_bound <= gap + 1e-9 && sep.observed >= gap - 1e-9;
            ok &= sep.observed - sep.lower_bound < 0.6; // slack is finite and stated
            worst_gap = worst_gap.max(sep.observed - gap);
            // Field law: c(p) ≥ separation at random points.
            let field = fs_query::ClearanceField { a: &a, b: &b };
            for _ in 0..50 {
                let p = Point3::new(
                    rng.range(-3.0, 3.0),
                    rng.range(-2.0, 2.0),
                    rng.range(-2.0, 2.0),
                );
                ok &= field.value(p, cx) >= gap - 1e-9;
            }
        }
        verdict(
            "gq-004",
            ok,
            &format!(
                "separation brackets hold across gaps 1.0 -> 0.05 (true separation in \
                 [lower_bound, observed], observed within {worst_gap:.1e} above \
                 truth), and the clearance field dominates the separation everywhere; \
                 seed 0x1001_2026_0707_0014"
            ),
        );
    });
}

/// gq-005 — the thickness oracle: exact on a graded slab, finds the
/// dumbbell neck, cross-checks against medial poles, and responds to a
/// DESIGN LEVER with the right subgradient (differentiable-friendly).
#[test]
#[allow(clippy::too_many_lines)] // slab, dumbbell, medial, and lever are one story
fn gq_005_thickness_oracle() {
    with_cx(|cx| {
        // Graded slab via F-rep: |z| ≤ (t0 + g·x)/2 within a box.
        // Thickness at (x, y, 0-top): t0 + g·x exactly.
        let slab = |t0: f64, g: f64| -> fs_rep_frep::Frep {
            let mut b = FrepBuilder::new();
            let bx = b
                .box_prim(Point3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 1.0, 1.5))
                .expect("box");
            // Half-space pair z ≤ (t0+g·x)/2 and −z ≤ (t0+g·x)/2 via
            // rotated half-spaces: n·p ≤ d with n = (−g/2, 0, 1)/|·|.
            let nz = 1.0f64;
            let nx = -g / 2.0;
            let nn = (nx * nx + nz * nz).sqrt();
            let top = b
                .half_space(Vec3::new(nx / nn, 0.0, nz / nn), t0 / (2.0 * nn))
                .expect("top");
            let bot = b
                .half_space(Vec3::new(nx / nn, 0.0, -nz / nn), t0 / (2.0 * nn))
                .expect("bot");
            let both = b
                .boolean(BoolOp::Intersect, BoolStyle::Hard, top, bot)
                .expect("b");
            let root = b
                .boolean(BoolOp::Intersect, BoolStyle::Hard, bx, both)
                .expect("r");
            b.finish(root).expect("frep")
        };
        let s = slab(0.4, 0.1);
        let mut ok = true;
        for x in [-1.0, 0.0, 1.0] {
            let expect = 0.4 + 0.1 * x;
            // Top surface point at height (t0+g·x)/2... solve: z = d(x).
            let z = expect / 2.0;
            // Project to be exactly on the boundary first.
            let cp = closest_point(&s, Point3::new(x, 0.0, z + 1e-3), cx).expect("cp");
            let t = thickness_at(&s, cp.point, cx).expect("thickness");
            // The inward-normal chord of a wedge is thickness/cos(tilt);
            // tilt is atan(g/2) — tiny; accept 1% relative.
            ok &= (t.value - expect).abs() / expect < 0.01;
        }
        // Dumbbell: two balls joined by a thin neck (hard union).
        let dumbbell = |neck_r: f64| -> fs_rep_frep::Frep {
            let mut b = FrepBuilder::new();
            let s1 = b.sphere(Point3::new(-1.2, 0.0, 0.0), 0.8).expect("s1");
            let s2 = b.sphere(Point3::new(1.2, 0.0, 0.0), 0.8).expect("s2");
            let neck = b
                .cylinder(Point3::new(0.0, 0.0, 0.0), neck_r)
                .expect("neck");
            // Cylinder is along z; rotate to x: rotate about y by 90°.
            let neck = b
                .rotate(neck, Vec3::new(0.0, 1.0, 0.0), core::f64::consts::FRAC_PI_2)
                .expect("rot");
            // Bound the infinite cylinder to the joint region.
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
        };
        let d = dumbbell(0.15);
        // Boundary samples on the neck: points at radius 0.15 around x=0.
        let mut samples = Vec::new();
        for k in 0..16 {
            let th = core::f64::consts::TAU * f64::from(k) / 16.0;
            samples.push(Point3::new(0.0, 0.15 * th.cos(), 0.15 * th.sin()));
        }
        let (min_t, skipped) = min_thickness(&d, &samples, cx).expect("min thickness");
        let neck_ok = (min_t - 0.3).abs() < 0.01 && skipped == 0;
        // Medial cross-check on the slab: poles' 2r matches thickness.
        let (hull, _) =
            fs_rep_mesh::dual_contour(&s, fs_rep_mesh::DcOptions::sharp(0.1), cx).expect("dc");
        let poles = medial_poles(&s, &hull, 1.2, cx).expect("poles");
        let mid_pole = poles
            .iter()
            .filter(|(p, _)| p.x.abs() < 0.4 && p.y.abs() < 0.4)
            .map(|(_, r)| 2.0 * r)
            .fold(f64::INFINITY, f64::min);
        let medial_agrees = (mid_pole - 0.4).abs() < 0.08;
        // Design-lever subgradient: d(min neck thickness)/d(neck_r) ≈ 2.
        let h = 1e-4;
        let (t_hi, _) = min_thickness(
            &dumbbell(0.15 + h),
            &{
                let mut v = Vec::new();
                for k in 0..16 {
                    let th = core::f64::consts::TAU * f64::from(k) / 16.0;
                    v.push(Point3::new(
                        0.0,
                        (0.15 + h) * th.cos(),
                        (0.15 + h) * th.sin(),
                    ));
                }
                v
            },
            cx,
        )
        .expect("hi");
        let (t_lo, _) = min_thickness(
            &dumbbell(0.15 - h),
            &{
                let mut v = Vec::new();
                for k in 0..16 {
                    let th = core::f64::consts::TAU * f64::from(k) / 16.0;
                    v.push(Point3::new(
                        0.0,
                        (0.15 - h) * th.cos(),
                        (0.15 - h) * th.sin(),
                    ));
                }
                v
            },
            cx,
        )
        .expect("lo");
        let subgrad = (t_hi - t_lo) / (2.0 * h);
        let lever_ok = (subgrad - 2.0).abs() < 1e-3;
        let mut em = fs_obs::Emitter::new("fs-query/conformance", "gq-005/thickness");
        let line = em
            .emit(
                fs_obs::Severity::Info,
                fs_obs::EventKind::Custom {
                    name: "query-thickness-oracle".to_string(),
                    json: format!(
                        "{{\"slab_ok\":{ok},\"neck_min\":{min_t:.4},\"medial_2r\":{mid_pole:.4},\
                         \"lever_subgradient\":{subgrad:.4}}}"
                    ),
                },
                None,
            )
            .to_jsonl();
        fs_obs::validate_line(&line).expect("thickness event validates");
        println!("{line}");
        verdict(
            "gq-005",
            ok && neck_ok && medial_agrees && lever_ok,
            &format!(
                "the graded slab reads its analytic thickness at three stations \
                 (1% rel), the dumbbell neck minimum is 2x the neck radius \
                 ({min_t:.3} vs 0.300, 0 skipped), medial poles cross-check the slab \
                 core (2r = {mid_pole:.3} vs 0.4), and the design-lever subgradient \
                 is {subgrad:.4} (analytic 2) — differentiable-friendly, demonstrated"
            ),
        );
    });
}

/// gq-006 — curvature: analytic values on sphere and torus, measured
/// O(h²) convergence for SecondOrder charts, the documented class per
/// chart family, and rotation invariance (G3).
#[test]
fn gq_006_curvature_convergence() {
    with_cx(|cx| {
        let exact = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 2.0,
        };
        let frep = frep_sphere(Point3::new(0.0, 0.0, 0.0), 2.0);
        let p = Point3::new(2.0, 0.0, 0.0);
        // Convergence order: errors at h and h/2 for the frep chart.
        let mut orders = Vec::new();
        for chart in [&exact as &dyn Chart, &frep as &dyn Chart] {
            let e = |h: f64| -> f64 {
                let c = curvature(chart, p, h, cx).expect("curv");
                (c.mean - 0.5).abs()
            };
            let (e1, e2) = (e(0.02), e(0.01));
            orders.push((e2 / e1.max(1e-300)).log2().abs());
        }
        let order_ok = orders.iter().all(|o| (*o - 2.0).abs() < 0.7);
        // Torus principal curvatures at the outer equator: 1/r and
        // 1/(R+r) — signs per outward convention.
        let torus = {
            let mut b = FrepBuilder::new();
            let t = b
                .torus(Point3::new(0.0, 0.0, 0.0), 1.0, 0.3)
                .expect("torus");
            b.finish(t).expect("frep")
        };
        let tp = Point3::new(1.3, 0.0, 0.0);
        let tc = curvature(&torus, tp, 0.01, cx).expect("torus curv");
        let (k1, k2) = (tc.principal[0], tc.principal[1]);
        let torus_ok = (k1 - 1.0 / 1.3).abs().min((k2 - 1.0 / 1.3).abs()) < 1e-2
            && (k1 - 1.0 / 0.3).abs().min((k2 - 1.0 / 0.3).abs()) < 1e-2
            && (tc.gaussian - (1.0 / 0.3) * (1.0 / 1.3)).abs() < 0.05;
        // Classes documented per chart family.
        let tiled = fs_rep_sdf::TiledSdf::build(&exact, 0.08, cx).expect("tiled");
        let mesh = fs_rep_mesh::MeshChart::new(fs_rep_mesh::shapes::icosphere(
            Point3::new(0.0, 0.0, 0.0),
            2.0,
            4,
        ));
        let classes_ok = curvature(&frep, p, 0.01, cx).expect("c").class
            == CurvatureClass::SecondOrder
            && fs_query::curvature_class(&tiled) == CurvatureClass::GridLimited
            && fs_query::curvature_class(&mesh) == CurvatureClass::Estimate;
        // Grid-limited chart still lands near truth at its own scale.
        let ct = curvature(&tiled, p, 0.08, cx).expect("tiled curv");
        let tiled_ok = (ct.mean - 0.5).abs() < 0.08;
        // G3: rotation invariance of curvature scalars (frep rotated).
        let rot = {
            let mut b = FrepBuilder::new();
            let s = b.sphere(Point3::new(0.0, 0.0, 0.0), 2.0).expect("s");
            let r = b.rotate(s, Vec3::new(0.3, -0.5, 0.8), 0.7).expect("rot");
            b.finish(r).expect("frep")
        };
        let cr = curvature(&rot, p, 0.01, cx).expect("rot curv");
        let g3 = (cr.mean - 0.5).abs() < 1e-4 && (cr.gaussian - 0.25).abs() < 1e-3;
        let mut em = fs_obs::Emitter::new("fs-query/conformance", "gq-006/curvature");
        let line = em
            .emit(
                fs_obs::Severity::Info,
                fs_obs::EventKind::Custom {
                    name: "query-curvature-convergence".to_string(),
                    json: format!(
                        "{{\"orders\":[{:.2},{:.2}],\"torus_k\":[{k1:.4},{k2:.4}],\
                         \"tiled_mean\":{:.4}}}",
                        orders[0], orders[1], ct.mean
                    ),
                },
                None,
            )
            .to_jsonl();
        fs_obs::validate_line(&line).expect("curvature event validates");
        println!("{line}");
        verdict(
            "gq-006",
            order_ok && torus_ok && classes_ok && tiled_ok && g3,
            &format!(
                "mean curvature converges at measured order ~2 on SecondOrder charts \
                 ({:.2}, {:.2}), torus principal curvatures hit 1/r and 1/(R+r), \
                 accuracy classes are documented per family (grid-limited lands \
                 within its own scale), and curvature scalars are rotation-invariant",
                orders[0], orders[1]
            ),
        );
    });
}

/// Teaching-refusal spot checks.
#[test]
fn refusals_teach() {
    with_cx(|cx| {
        let exact = SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let off = thickness_at(&exact, Point3::new(1.5, 0.0, 0.0), cx).expect_err("off boundary");
        assert!(matches!(off, QueryError::NotOnBoundary { .. }), "{off}");
        assert!(off.to_string().contains("project"), "{off}");
    });
}
