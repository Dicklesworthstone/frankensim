//! fs-rep-sdf conformance suite (CONTRACT.md: any reimplementation must
//! pass). Fixture reproduction within declared bounds, C¹ seamlessness,
//! eikonal evidence, the VDB-vs-oracle property battery, narrow-band
//! drift, and sphere tracing. JSON-line verdicts; seeded cases carry
//! seeds.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::{BoxChart, SphereChart, TorusChart};
use fs_geom::{Aabb, Chart, Point3, Vec3};
use fs_rep_sdf::{AdaptiveSdf, NarrowBand, TiledSdf, VdbGrid};
use std::collections::BTreeMap;

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-rep-sdf/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\
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

    fn below(&mut self, n: u64) -> u64 {
        (self.next() >> 32) % n
    }
}

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x5DF,
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

fn sphere() -> SphereChart {
    SphereChart {
        center: Point3::new(0.1, -0.2, 0.05),
        radius: 1.3,
    }
}

#[test]
fn rsdf_001_fixtures_reproduced_within_declared_bounds() {
    const SEED: u64 = 0x5DF1_2026_0706_0001;
    let gate = CancelGate::new();
    let mut rng = Lcg(SEED);
    let mut worst_ratio = 0.0f64;
    with_cx(&gate, |cx| {
        let sources: Vec<Box<dyn Chart>> = vec![
            Box::new(sphere()),
            Box::new(BoxChart {
                aabb: Aabb::new(Point3::new(-1.0, -0.7, -0.9), Point3::new(0.8, 1.0, 0.6)),
            }),
            Box::new(TorusChart {
                center: Point3::new(0.0, 0.0, 0.0),
                major: 1.6,
                minor: 0.5,
            }),
        ];
        for source in &sources {
            let grid = TiledSdf::build(source.as_ref(), 0.05, cx).expect("build");
            let support = grid.support();
            for _ in 0..4_000 {
                let p = Point3::new(
                    support.min.x + (support.max.x - support.min.x) * rng.unit(),
                    support.min.y + (support.max.y - support.min.y) * rng.unit(),
                    support.min.z + (support.max.z - support.min.z) * rng.unit(),
                );
                let err = (grid.eval(p, cx).signed_distance
                    - source.eval(p, cx).signed_distance)
                    .abs();
                worst_ratio = worst_ratio.max(err / grid.bound());
            }
        }
    });
    verdict(
        "rsdf-001",
        worst_ratio <= 1.0,
        &format!(
            "sphere/box/torus reproduced within the declared enclosure over 12k seeded points \
             (worst error / bound = {worst_ratio:.3}, seed {SEED:#x})"
        ),
    );
}

#[test]
fn rsdf_002_c1_reconstruction_is_seamless_and_gradients_match_fd() {
    let gate = CancelGate::new();
    let (seam_ok, grad_ok) = with_cx(&gate, |cx| {
        let src = sphere();
        let grid = TiledSdf::build(&src, 0.04, cx).expect("build");
        // Walk a line crossing many tile boundaries; value and gradient
        // must vary continuously (no seams at 8-sample boundaries).
        let mut prev: Option<(f64, Vec3)> = None;
        let mut seam_ok = true;
        let step = 0.003;
        for i in 0..800 {
            let p = Point3::new(-1.2 + f64::from(i) * step, 0.11, -0.07);
            let s = grid.eval(p, cx);
            let g = s.gradient.expect("in-box gradient");
            if let Some((pv, pg)) = prev {
                seam_ok &= (s.signed_distance - pv).abs() < 0.05;
                let dg = Vec3::new(g.x - pg.x, g.y - pg.y, g.z - pg.z).norm();
                seam_ok &= dg < 0.15; // C¹: gradient jumps stay small
            }
            prev = Some((s.signed_distance, g));
        }
        // Gradient vs central FD at scattered points.
        let mut grad_ok = true;
        let h = 1e-5;
        for i in 0..200 {
            let p = Point3::new(
                -0.9 + f64::from(i % 20) * 0.09,
                -0.8 + f64::from((i / 20) % 10) * 0.17,
                0.03,
            );
            let g = grid.eval(p, cx).gradient.expect("gradient");
            let f = |q: Point3| grid.eval(q, cx).signed_distance;
            let fd = Vec3::new(
                (f(p.offset(Vec3::new(h, 0.0, 0.0))) - f(p.offset(Vec3::new(-h, 0.0, 0.0))))
                    / (2.0 * h),
                (f(p.offset(Vec3::new(0.0, h, 0.0))) - f(p.offset(Vec3::new(0.0, -h, 0.0))))
                    / (2.0 * h),
                (f(p.offset(Vec3::new(0.0, 0.0, h))) - f(p.offset(Vec3::new(0.0, 0.0, -h))))
                    / (2.0 * h),
            );
            grad_ok &= Vec3::new(g.x - fd.x, g.y - fd.y, g.z - fd.z).norm() < 1e-3;
        }
        (seam_ok, grad_ok)
    });
    verdict(
        "rsdf-002",
        seam_ok && grad_ok,
        "triquadratic reconstruction is continuous across tile boundaries with C1 gradients \
         matching central differences (G0)",
    );
}

#[test]
fn rsdf_003_eikonal_evidence_is_measured_and_ledgered() {
    let gate = CancelGate::new();
    let stats = with_cx(&gate, |cx| {
        let grid = TiledSdf::build(&sphere(), 0.05, cx).expect("build");
        grid.measure_eikonal(0x5DF3, 4_000, cx)
    });
    let mut em = fs_obs::Emitter::new("fs-rep-sdf/conformance", "rsdf-003/eikonal");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "rep-sdf-eikonal-stats".to_string(),
                json: format!(
                    "{{\"mean_abs_dev\":{:.6},\"max_abs_dev\":{:.6},\"probes\":{}}}",
                    stats.mean_abs_dev, stats.max_abs_dev, stats.probes
                ),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("eikonal event validates");
    println!("{line}");
    verdict(
        "rsdf-003",
        stats.mean_abs_dev < 0.05 && stats.max_abs_dev < 0.6 && stats.probes == 4_000,
        &format!(
            "eikonal deviation measured and ledgered: mean {:.4}, max {:.4} (a sampled sphere \
             is nearly a distance field; the deviation is EVIDENCE, not a certificate)",
            stats.mean_abs_dev, stats.max_abs_dev
        ),
    );
}

#[test]
fn rsdf_004_vdb_matches_the_oracle_and_reports_footprint() {
    const SEED: u64 = 0x5DF4_2026_0706_0BDB;
    let mut rng = Lcg(SEED);
    let mut vdb: VdbGrid<f32> = VdbGrid::new(-1.0);
    let mut oracle: BTreeMap<[i32; 3], f32> = BTreeMap::new();
    // Clustered actives (three blobs) + scattered strays, interleaved
    // reads — the sparse-use shape.
    for blob in 0..3i32 {
        let base = [blob * 40 - 30, blob * 25 - 20, blob * 15 - 10];
        for _ in 0..4_000 {
            let c = [
                base[0] + (rng.below(16) as i32) - 8,
                base[1] + (rng.below(16) as i32) - 8,
                base[2] + (rng.below(16) as i32) - 8,
            ];
            let v = (rng.unit() * 4.0 - 2.0) as f32;
            vdb.set(c, v);
            oracle.insert(c, v);
        }
    }
    for _ in 0..500 {
        let c = [
            (rng.below(4_000) as i32) - 2_000,
            (rng.below(4_000) as i32) - 2_000,
            (rng.below(4_000) as i32) - 2_000,
        ];
        let v = rng.unit() as f32;
        vdb.set(c, v);
        oracle.insert(c, v);
    }
    // Exact agreement: actives, values, iteration set, misses.
    let mut agree = vdb.active_count() == oracle.len() as u64;
    for (&c, &v) in &oracle {
        agree &= vdb.is_active(c) && (vdb.get(c) - v).abs() == 0.0;
    }
    let iterated: BTreeMap<[i32; 3], f32> = vdb.iter_active().collect();
    agree &= iterated == oracle;
    for _ in 0..2_000 {
        let c = [
            (rng.below(8_000) as i32) - 4_000,
            (rng.below(8_000) as i32) - 4_000,
            (rng.below(8_000) as i32) - 4_000,
        ];
        if !oracle.contains_key(&c) {
            agree &= !vdb.is_active(c) && vdb.get(c) == -1.0;
        }
    }
    // Dilate grows, erode shrinks back below the dilated count.
    let before = vdb.active_count();
    vdb.dilate();
    let dilated = vdb.active_count();
    vdb.erode();
    let eroded = vdb.active_count();
    let topo_ok = dilated > before && eroded < dilated;
    let stats = vdb.memory_stats();
    let mut em = fs_obs::Emitter::new("fs-rep-sdf/conformance", "rsdf-004/vdb");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "rep-sdf-vdb-stats".to_string(),
                json: stats.to_json(),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("vdb stats validate");
    println!("{line}");
    verdict(
        "rsdf-004",
        agree && topo_ok,
        &format!(
            "VDB agrees exactly with the BTreeMap oracle over ~12.5k clustered+stray actives \
             (seed {SEED:#x}); dilate/erode behave; footprint ledgered: {}",
            stats.to_json()
        ),
    );
}

#[test]
fn rsdf_005_narrow_band_advects_with_bounded_drift_and_reinit_helps() {
    let gate = CancelGate::new();
    let (drift, dev_before, dev_after, stats_json) = with_cx(&gate, |cx| {
        let h = 0.06;
        let mut band =
            NarrowBand::from_chart(&sphere(), h, 5, cx).expect("band build not cancelled");
        // Translate at constant velocity; the zero crossing must follow.
        let v = Vec3::new(0.05, 0.0, 0.0);
        let (dt, steps) = (0.5, 6);
        for _ in 0..steps {
            band.advect(|_| v, dt);
        }
        let moved = v.scale(dt * f64::from(steps));
        let expected_center = sphere().center.offset(moved);
        // φ at points on the ANALYTIC translated boundary should be ~0.
        let mut worst = 0.0f64;
        for i in 0..24 {
            let a = f64::from(i) * core::f64::consts::TAU / 24.0;
            let p = Point3::new(
                expected_center.x + 1.3 * a.cos(),
                expected_center.y + 1.3 * a.sin(),
                expected_center.z,
            );
            if let Some(phi) = band.interpolate(p) {
                worst = worst.max(phi.abs());
            }
        }
        let dev_before = band.stats().mean_eikonal_dev;
        band.reinitialize(20);
        let s = band.stats();
        (worst, dev_before, s.mean_eikonal_dev, s.to_json())
    });
    let mut em = fs_obs::Emitter::new("fs-rep-sdf/conformance", "rsdf-005/band");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "rep-sdf-band-stats".to_string(),
                json: stats_json.clone(),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("band stats validate");
    println!("{line}");
    // First-order semi-Lagrangian on h=0.06: allow a few cells of drift,
    // measured and ledgered (tightening is the topo-levelset bead's).
    verdict(
        "rsdf-005",
        drift < 0.15 && dev_after <= dev_before + 1e-9,
        &format!(
            "sphere band translated 6 steps: zero-crossing drift {drift:.4} (< 0.15); \
             reinit does not worsen the eikonal residual ({dev_before:.4} -> {dev_after:.4}); \
             band stats ledgered"
        ),
    );
}

#[test]
fn rsdf_006_sphere_tracing_respects_its_own_bounds_and_adaptive_builds() {
    let gate = CancelGate::new();
    let (hits_ok, misses_ok, adaptive_json, adaptive_ok) = with_cx(&gate, |cx| {
        let src = sphere();
        let grid = TiledSdf::build(&src, 0.04, cx).expect("build");
        // Rays at the sphere from scattered directions: hit t must match
        // the analytic intersection within the chart's bound + one step.
        let mut hits_ok = true;
        for i in 0..16 {
            let a = f64::from(i) * core::f64::consts::TAU / 16.0;
            let origin = Point3::new(
                src.center.x + 4.0 * a.cos(),
                src.center.y + 4.0 * a.sin(),
                src.center.z,
            );
            let dir = src.center.delta_from(origin);
            let t = grid.raycast(origin, dir, 10.0, cx);
            let expected = dir.norm() - src.radius;
            match t {
                Some(t) => {
                    hits_ok &= (t - expected).abs() < grid.bound() * 4.0 + 0.02;
                }
                None => hits_ok = false,
            }
        }
        // Rays that miss must miss.
        let misses_ok = grid
            .raycast(
                Point3::new(4.0, 4.0, 4.0),
                Vec3::new(1.0, 0.0, 0.0),
                10.0,
                cx,
            )
            .is_none();
        // Adaptive octree: builds under tolerance, residual ledgered,
        // reproduces the source within its Estimate band on probe points.
        let ad = AdaptiveSdf::build(&src, 0.02, 6, cx).expect("adaptive build");
        let stats = ad.stats();
        let mut ok = stats.residual <= 0.02 + 1e-12;
        let mut rng = Lcg(0x5DF6);
        for _ in 0..2_000 {
            let s = ad.support();
            let p = Point3::new(
                s.min.x + (s.max.x - s.min.x) * rng.unit(),
                s.min.y + (s.max.y - s.min.y) * rng.unit(),
                s.min.z + (s.max.z - s.min.z) * rng.unit(),
            );
            let err = (ad.eval(p, cx).signed_distance - src.eval(p, cx).signed_distance).abs();
            // Probed residual is Estimate-grade: allow slack between probes.
            ok &= err < 0.08;
        }
        (hits_ok, misses_ok, stats.to_json(), ok)
    });
    let mut em = fs_obs::Emitter::new("fs-rep-sdf/conformance", "rsdf-006/adaptive");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "rep-sdf-adaptive-stats".to_string(),
                json: adaptive_json.clone(),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("adaptive stats validate");
    println!("{line}");
    verdict(
        "rsdf-006",
        hits_ok && misses_ok && adaptive_ok,
        &format!(
            "sphere tracing hits within bound and misses cleanly; adaptive octree residual \
             ledgered: {adaptive_json}"
        ),
    );
}
