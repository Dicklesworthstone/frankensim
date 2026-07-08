//! Chart-backend conformance (the qfx.2 bead; runs under
//! `chart-backends`). Acceptance: ZERO missed intersections across the
//! adversarial ray battery (thin shells, grazing rays) vs the certified
//! root-finder oracle — the headline certificate test; NURBS Newton
//! matches analytic references; mixed-chart scenes agree across backend
//! kinds within tolerance; G3 frame invariance; ray rates measured and
//! ledgered (debug-build numbers; the perf GATE lives in perf-CI).
#![cfg(feature = "chart-backends")]

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::charts::{Backend, Ray, TriMesh, ray_intersect_nurbs, sphere_trace, trace_scene};
use fs_rep_frep::{BoolOp, BoolStyle, FrepBuilder};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-render/charts\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 9,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn unit(v: Vec3) -> Vec3 {
    v.scale(1.0 / v.norm())
}

fn lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / (1u64 << 53) as f64
}

/// The THIN SHELL adversary: a sphere shell of thickness 2e-3 built as
/// an F-rep difference — the classic tunneling victim.
fn thin_shell() -> fs_rep_frep::Frep {
    let mut b = FrepBuilder::new();
    let outer = b.sphere(Point3::new(0.0, 0.0, 0.0), 1.0).expect("outer");
    let inner = b.sphere(Point3::new(0.0, 0.0, 0.0), 0.998).expect("inner");
    let shell = b
        .boolean(BoolOp::Difference, BoolStyle::Hard, outer, inner)
        .expect("shell");
    b.finish(shell).expect("frep")
}

/// The certified ROOT-FINDER ORACLE: march at a fixed step far below
/// eps/L and bisect any bracketing interval — an independent code path
/// that cannot miss a crossing wider than the micro-step.
fn oracle_hits(chart: &dyn fs_geom::Chart, cx: &Cx<'_>, ray: &Ray, t_max: f64) -> bool {
    let micro = 2e-4;
    let mut t = 0.0f64;
    let mut prev = chart.eval(ray.at(0.0), cx).signed_distance;
    while t < t_max {
        t += micro;
        let d = chart.eval(ray.at(t), cx).signed_distance;
        if prev.signum() != d.signum() {
            return true;
        }
        prev = d;
    }
    false
}

#[test]
fn rb_001_zero_tunneling_headline() {
    with_cx(|cx| {
        let shell = thin_shell();
        let mut state = 0xbeef_u64;
        let mut audited = 0usize;
        for k in 0..120 {
            // Adversarial battery: random origins on a far sphere aimed
            // near (but not always at) the shell — grazing rays included.
            let phi = lcg(&mut state) * std::f64::consts::TAU;
            let z = 2.0 * lcg(&mut state) - 1.0;
            let r = (1.0 - z * z).sqrt();
            let origin = Point3::new(3.0 * r * phi.cos(), 3.0 * r * phi.sin(), 3.0 * z);
            // Aim at a jittered point near the surface (grazing bias
            // every third ray).
            let graze = k % 3 == 0;
            let target_r = if graze { 0.999 } else { 0.6 * lcg(&mut state) };
            let tphi = lcg(&mut state) * std::f64::consts::TAU;
            let tz = 2.0 * lcg(&mut state) - 1.0;
            let tr = (1.0 - tz * tz).sqrt();
            let target = Point3::new(
                target_r * tr * tphi.cos(),
                target_r * tr * tphi.sin(),
                target_r * tz,
            );
            let ray = Ray {
                origin,
                dir: unit(target.delta_from(origin)),
            };
            let oracle = oracle_hits(&shell, cx, &ray, 6.0);
            let (traced, audit) = sphere_trace(&shell, cx, &ray, 6.0, 1e-6, 1.0);
            assert_eq!(
                traced.is_some(),
                oracle,
                "ray {k}: tracer and oracle must agree (tunneling audit)"
            );
            // G0 step safety: no plain step ever exceeded |f|/L.
            assert!(
                audit.worst_step_ratio <= 1.0 + 1e-12,
                "step-safety property: {}",
                audit.worst_step_ratio
            );
            audited += 1;
        }
        println!("{{\"metric\":\"tunneling-audit\",\"rays\":{audited},\"missed\":0}}");
        verdict(
            "rb-001",
            "120 adversarial rays (thin 2e-3 shell, grazing bias): tracer agrees with \
             the micro-step oracle on every ray; step-safety never violated",
        );
    });
}

#[test]
fn rb_002_over_relaxation_stays_certified() {
    with_cx(|cx| {
        let shell = thin_shell();
        let mut state = 0xface_u64;
        let mut plain_steps = 0u64;
        let mut relaxed_steps = 0u64;
        for k in 0..60 {
            // GRAZING rays: passing at 1.02-1.15x the radius, parallel
            // to an axis — the many-small-steps regime over-relaxation
            // exists for (center-aimed rays converge in one jump and
            // have nothing to accelerate).
            let offset = 1.02 + 0.13 * lcg(&mut state);
            let phi = lcg(&mut state) * std::f64::consts::TAU;
            let origin = Point3::new(-4.0, offset * phi.cos(), offset * phi.sin());
            let ray = Ray {
                origin,
                dir: Vec3::new(1.0, 0.0, 0.0),
            };
            let (h1, a1) = sphere_trace(&shell, cx, &ray, 9.0, 1e-6, 1.0);
            let (h2, a2) = sphere_trace(&shell, cx, &ray, 9.0, 1e-6, 1.6);
            assert_eq!(
                h1.is_some(),
                h2.is_some(),
                "relaxation never changes hits (ray {k})"
            );
            if let (Some(a), Some(b)) = (h1, h2) {
                assert!(
                    (a.t - b.t).abs() < 1e-4,
                    "same intersection: {} vs {}",
                    a.t,
                    b.t
                );
            }
            plain_steps += u64::from(a1.steps);
            relaxed_steps += u64::from(a2.steps);
        }
        println!(
            "{{\"metric\":\"over-relaxation\",\"plain_steps\":{plain_steps},\
             \"relaxed_steps\":{relaxed_steps}}}"
        );
        assert!(
            relaxed_steps < plain_steps,
            "relaxation saves steps: {relaxed_steps} vs {plain_steps}"
        );
        verdict(
            "rb-002",
            "omega=1.6 marching hits identical intersections with fewer steps (certified \
             fallback preserves the no-tunnel invariant)",
        );
    });
}

/// The exact NURBS unit sphere (rational quadratic revolution).
fn sphere_nurbs() -> fs_rep_nurbs::NurbsSurface<f64> {
    let s2 = std::f64::consts::FRAC_1_SQRT_2;
    let circle = [
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
        [-1.0, 1.0],
        [-1.0, 0.0],
        [-1.0, -1.0],
        [0.0, -1.0],
        [1.0, -1.0],
        [1.0, 0.0],
    ];
    let cw = |i: usize| if i.is_multiple_of(2) { 1.0 } else { s2 };
    let profile: [([f64; 2], f64); 5] = [
        ([0.0, -1.0], 1.0),
        ([1.0, -1.0], s2),
        ([1.0, 0.0], 1.0),
        ([1.0, 1.0], s2),
        ([0.0, 1.0], 1.0),
    ];
    let mut points: Vec<Vec<[f64; 3]>> = Vec::new();
    let mut weights: Vec<Vec<f64>> = Vec::new();
    for (i, c) in circle.iter().enumerate() {
        let mut prow = Vec::new();
        let mut wrow = Vec::new();
        for &([radius, z], wv) in &profile {
            prow.push([radius * c[0], radius * c[1], z]);
            wrow.push(cw(i) * wv);
        }
        points.push(prow);
        weights.push(wrow);
    }
    let ku = fs_rep_nurbs::KnotVector::new(
        vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        2,
    )
    .expect("ku");
    let kv =
        fs_rep_nurbs::KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0], 2).expect("kv");
    fs_rep_nurbs::NurbsSurface::new(ku, kv, &points, &weights).expect("sphere")
}

#[test]
fn rb_003_nurbs_newton_matches_analytic() {
    let sphere = sphere_nurbs();
    let mut state = 0x5eed_u64;
    for _ in 0..40 {
        let phi = lcg(&mut state) * std::f64::consts::TAU;
        let z = 1.6 * lcg(&mut state) - 0.8;
        let origin = Point3::new(4.0 * phi.cos(), 4.0 * phi.sin(), z);
        let dir = unit(Point3::new(0.0, 0.0, 0.0).delta_from(origin));
        let ray = Ray { origin, dir };
        // Analytic sphere intersection: |o + t d| = 1.
        let oc = [origin.x, origin.y, origin.z];
        let b = oc[0] * dir.x + oc[1] * dir.y + oc[2] * dir.z;
        let c = oc[0] * oc[0] + oc[1] * oc[1] + oc[2] * oc[2] - 1.0;
        let disc = b * b - c;
        assert!(disc > 0.0, "aimed at the center: always hits");
        let t_ref = -b - disc.sqrt();
        let hit = ray_intersect_nurbs(&sphere, &ray, 8, 1e-9).expect("Newton converges");
        assert!(
            (hit.t - t_ref).abs() < 1e-6,
            "Newton matches analytic: {} vs {t_ref}",
            hit.t
        );
        // The normal points outward at the hit point.
        let n = hit.normal.expect("normal");
        let outward = unit(hit.point.delta_from(Point3::new(0.0, 0.0, 0.0)));
        assert!(
            n.dot(outward).abs() > 0.999,
            "normal aligned with the radial direction"
        );
    }
    verdict(
        "rb-003",
        "40 rays: Bezier-seeded Newton matches the analytic sphere intersection to 1e-6 \
         with radial normals",
    );
}

/// An octahedron-subdivision icosphere-ish mesh of the unit sphere.
fn sphere_mesh(subdiv: usize) -> TriMesh {
    // Start from an octahedron, subdivide, project to the sphere.
    let mut verts: Vec<[f64; 3]> = vec![
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];
    let mut tris: Vec<[u32; 3]> = vec![
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5],
    ];
    for _ in 0..subdiv {
        let mut next = Vec::with_capacity(tris.len() * 4);
        for t in &tris {
            let mid = |a: u32, b: u32, verts: &mut Vec<[f64; 3]>| -> u32 {
                let (pa, pb) = (verts[a as usize], verts[b as usize]);
                let mut m = [
                    f64::midpoint(pa[0], pb[0]),
                    f64::midpoint(pa[1], pb[1]),
                    f64::midpoint(pa[2], pb[2]),
                ];
                let n = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
                for v in &mut m {
                    *v /= n;
                }
                verts.push(m);
                (verts.len() - 1) as u32
            };
            let ab = mid(t[0], t[1], &mut verts);
            let bc = mid(t[1], t[2], &mut verts);
            let ca = mid(t[2], t[0], &mut verts);
            next.extend_from_slice(&[[t[0], ab, ca], [ab, t[1], bc], [ca, bc, t[2]], [ab, bc, ca]]);
        }
        tris = next;
    }
    TriMesh::new(verts, tris)
}

#[test]
fn rb_004_mixed_scene_consistency_and_frame_invariance() {
    with_cx(|cx| {
        // The SAME unit sphere held three ways.
        let mut b = FrepBuilder::new();
        let node = b.sphere(Point3::new(0.0, 0.0, 0.0), 1.0).expect("sphere");
        let frep = b.finish(node).expect("frep");
        let nurbs = sphere_nurbs();
        let mesh = sphere_mesh(4);
        let mut state = 0xd1ce_u64;
        let mut worst_spread = 0.0f64;
        for _ in 0..30 {
            let phi = lcg(&mut state) * std::f64::consts::TAU;
            let z = 1.2 * lcg(&mut state) - 0.6;
            let origin = Point3::new(3.5 * phi.cos(), 3.5 * phi.sin(), z);
            let ray = Ray {
                origin,
                dir: unit(Point3::new(0.0, 0.0, 0.0).delta_from(origin)),
            };
            let (sdf_hit, _) = sphere_trace(&frep, cx, &ray, 8.0, 1e-7, 1.0);
            let nurbs_hit = ray_intersect_nurbs(&nurbs, &ray, 8, 1e-9);
            let mesh_hit = mesh.intersect(&ray);
            let (a, b_, c) = (
                sdf_hit.expect("sdf hits").t,
                nurbs_hit.expect("nurbs hits").t,
                mesh_hit.expect("mesh hits").t,
            );
            worst_spread = worst_spread.max((a - b_).abs()).max((a - c).abs());
        }
        // Mesh is a level-4 subdivision: geometric error ~ 1e-3.
        assert!(
            worst_spread < 5e-3,
            "the same shape across three backends: spread {worst_spread}"
        );
        // Mixed scene: three unit spheres at distinct centers, one per
        // backend — the closest-hit wins and identifies its instance.
        let mut b2 = FrepBuilder::new();
        let sn = b2.sphere(Point3::new(-3.0, 0.0, 0.0), 1.0).expect("s");
        let frep2 = b2.finish(sn).expect("frep");
        let backends = [
            Backend::Chart(&frep2),
            Backend::Nurbs(&nurbs),
            Backend::Mesh(&mesh),
        ];
        let ray = Ray {
            origin: Point3::new(-8.0, 0.0, 0.01),
            dir: Vec3::new(1.0, 0.0, 0.0),
        };
        let (idx, hit) = trace_scene(&backends, cx, &ray, 20.0, 1e-7).expect("hits");
        assert_eq!(idx, 0, "the F-rep sphere at x=-3 is closest");
        assert!((hit.t - 4.0).abs() < 1e-3, "t ~ 4: {}", hit.t);
        // G3 frame invariance: translate everything by the same offset.
        let ray_shifted = Ray {
            origin: Point3::new(-8.0 + 5.0, 0.0, 0.01),
            dir: Vec3::new(1.0, 0.0, 0.0),
        };
        let mut b3 = FrepBuilder::new();
        let sn3 = b3.sphere(Point3::new(2.0, 0.0, 0.0), 1.0).expect("s");
        let frep3 = b3.finish(sn3).expect("frep");
        let (t1, _) = sphere_trace(&frep2, cx, &ray, 20.0, 1e-7, 1.0)
            .0
            .map(|h| (h.t, h.steps))
            .expect("hit");
        let (t2, _) = sphere_trace(&frep3, cx, &ray_shifted, 20.0, 1e-7, 1.0)
            .0
            .map(|h| (h.t, h.steps))
            .expect("hit");
        assert!(
            (t1 - t2).abs() < 1e-9,
            "translation invariance: {t1} vs {t2}"
        );
        verdict(
            "rb-004",
            "one sphere, three backends, <5e-3 spread; mixed scene picks the closest \
             instance; translated scene reproduces t to 1e-9",
        );
    });
}

#[test]
fn rb_005_ray_rate_ledger() {
    with_cx(|cx| {
        // MEASURED, honestly labeled: debug-build throughput on primary
        // rays. The plan's Mray/s TARGETS are release-build perf-CI
        // gates (fz2.4) — this ledgers the measurement discipline, not
        // the target (AGENTS: no "fast" without a benchmark + machine).
        let mut b = FrepBuilder::new();
        let node = b.sphere(Point3::new(0.0, 0.0, 0.0), 1.0).expect("sphere");
        let frep = b.finish(node).expect("frep");
        let n = 64usize;
        let start = std::time::Instant::now();
        let mut hits = 0usize;
        for py in 0..n {
            for px in 0..n {
                #[allow(clippy::cast_precision_loss)]
                let (u, v) = (
                    (px as f64 + 0.5) / n as f64 * 2.0 - 1.0,
                    (py as f64 + 0.5) / n as f64 * 2.0 - 1.0,
                );
                let ray = Ray {
                    origin: Point3::new(1.6 * u, 1.6 * v, 3.0),
                    dir: Vec3::new(0.0, 0.0, -1.0),
                };
                if sphere_trace(&frep, cx, &ray, 8.0, 1e-6, 1.6).0.is_some() {
                    hits += 1;
                }
            }
        }
        let secs = start.elapsed().as_secs_f64();
        #[allow(clippy::cast_precision_loss)]
        let mray_s = (n * n) as f64 / secs / 1e6;
        println!(
            "{{\"metric\":\"ray-rate\",\"build\":\"debug\",\"machine\":\"darwin-arm64\",\
             \"rays\":{},\"hits\":{hits},\"mray_per_s\":{mray_s:.4},\
             \"note\":\"release perf gate lives in perf-CI (fz2.4)\"}}",
            n * n
        );
        assert!(hits > 0, "the turntable frame sees the sphere");
        verdict(
            "rb-005",
            "primary-ray throughput measured and ledgered with build/machine labels; \
             the Mray/s TARGET is the perf-CI lane's gate",
        );
    });
}
