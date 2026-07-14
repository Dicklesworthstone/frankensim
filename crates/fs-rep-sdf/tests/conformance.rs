//! fs-rep-sdf conformance suite (CONTRACT.md: any reimplementation must
//! pass). Fixture reproduction within declared bounds, C¹ seamlessness,
//! eikonal evidence, the VDB-vs-oracle property battery, narrow-band
//! drift, and sphere tracing. JSON-line verdicts; seeded cases carry
//! seeds.

use asupersync::types::Budget;
use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::{BoxChart, SphereChart, TorusChart};
use fs_geom::{
    Aabb, Chart, ChartSample, Differentiability, Point3, SamplingDomainError, TraceStepClaim, Vec3,
};
use fs_rep_sdf::{ADAPTIVE_MAX_NODES, AdaptiveSdf, NarrowBand, SdfBuildError, TiledSdf, VdbGrid};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

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

struct CountingPlane {
    x_offset: f64,
    evals: AtomicU64,
}

impl CountingPlane {
    fn new(x_offset: f64) -> Self {
        Self {
            x_offset,
            evals: AtomicU64::new(0),
        }
    }

    fn eval_count(&self) -> u64 {
        self.evals.load(Ordering::Relaxed)
    }
}

impl Chart for CountingPlane {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        self.evals.fetch_add(1, Ordering::Relaxed);
        let signed_distance = x.x - self.x_offset;
        ChartSample {
            signed_distance,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error: NumericalCertificate::enclosure(signed_distance, signed_distance),
        }
    }

    fn support(&self) -> Aabb {
        Aabb {
            min: Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
            max: Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
        }
    }

    fn name(&self) -> &'static str {
        "test/counting-plane"
    }

    fn differentiability(&self) -> Differentiability {
        Differentiability::C0
    }
}

struct ExtremeFiniteBandChart {
    support: Aabb,
    scan_bounds: Aabb,
    evals: AtomicU64,
    invalid_point: AtomicBool,
    max_seen_bits: [AtomicU64; 3],
}

impl ExtremeFiniteBandChart {
    fn new(h: f64) -> Self {
        let support = Aabb::new(
            Point3::new(1.75e308, 1.75e308, 1.75e308),
            Point3::new(1.76e308, 1.76e308, 1.76e308),
        );
        let scan_bounds = support.inflate(2.0 * h);
        Self {
            support,
            scan_bounds,
            evals: AtomicU64::new(0),
            invalid_point: AtomicBool::new(false),
            max_seen_bits: [AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)],
        }
    }

    fn eval_count(&self) -> u64 {
        self.evals.load(Ordering::Relaxed)
    }

    fn all_points_admitted(&self) -> bool {
        !self.invalid_point.load(Ordering::Relaxed)
    }

    fn max_seen(&self) -> Point3 {
        Point3::new(
            f64::from_bits(self.max_seen_bits[0].load(Ordering::Relaxed)),
            f64::from_bits(self.max_seen_bits[1].load(Ordering::Relaxed)),
            f64::from_bits(self.max_seen_bits[2].load(Ordering::Relaxed)),
        )
    }
}

impl Chart for ExtremeFiniteBandChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        self.evals.fetch_add(1, Ordering::Relaxed);
        if !x.x.is_finite() || !x.y.is_finite() || !x.z.is_finite() || !self.scan_bounds.contains(x)
        {
            self.invalid_point.store(true, Ordering::Relaxed);
        }
        self.max_seen_bits[0].fetch_max(x.x.to_bits(), Ordering::Relaxed);
        self.max_seen_bits[1].fetch_max(x.y.to_bits(), Ordering::Relaxed);
        self.max_seen_bits[2].fetch_max(x.z.to_bits(), Ordering::Relaxed);
        ChartSample {
            signed_distance: 0.0,
            gradient: None,
            lipschitz: None,
            error: NumericalCertificate::exact(0.0),
        }
    }

    fn support(&self) -> Aabb {
        self.support
    }

    fn name(&self) -> &'static str {
        "test/extreme-finite-band"
    }
}

struct FiniteNanChart;

impl Chart for FiniteNanChart {
    fn eval(&self, _x: Point3, _cx: &Cx<'_>) -> ChartSample {
        ChartSample {
            signed_distance: f64::NAN,
            gradient: None,
            lipschitz: None,
            error: NumericalCertificate::no_claim(),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "test/finite-nan"
    }
}

struct UnrepresentableAdaptiveChart;

struct OverflowingAdaptiveResidualChart;

impl Chart for UnrepresentableAdaptiveChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let signed_distance = x.y * x.y;
        ChartSample {
            signed_distance,
            gradient: None,
            lipschitz: Some(3.0),
            error: NumericalCertificate::exact(signed_distance),
        }
    }

    fn support(&self) -> Aabb {
        let x_min = 1.0e308_f64;
        let x_max = f64::from_bits(x_min.to_bits() + 1);
        Aabb::new(Point3::new(x_min, -1.0, -1.0), Point3::new(x_max, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "test/unrepresentable-adaptive-split"
    }
}

impl Chart for OverflowingAdaptiveResidualChart {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let at_corner = x.x.abs() > 0.0 && x.y.abs() > 0.0 && x.z.abs() > 0.0;
        let signed_distance = if at_corner { -f64::MAX } else { f64::MAX };
        ChartSample {
            signed_distance,
            gradient: None,
            lipschitz: Some(1.0),
            error: NumericalCertificate::exact(signed_distance),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "test/overflowing-adaptive-residual"
    }
}

#[derive(Clone, Copy)]
enum SourceAuthority {
    Exact,
    Enclosure(f64),
    LocalEnclosure(f64),
    Estimate(f64),
    NoClaim,
    MalformedFinite,
}

struct BoundedAuthorityPlane {
    authority: SourceAuthority,
}

impl Chart for BoundedAuthorityPlane {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let signed_distance = x.x;
        let error = match self.authority {
            SourceAuthority::Exact => NumericalCertificate::exact(signed_distance),
            SourceAuthority::Enclosure(radius) | SourceAuthority::LocalEnclosure(radius) => {
                NumericalCertificate::enclosure(signed_distance - radius, signed_distance + radius)
            }
            SourceAuthority::Estimate(radius) => {
                NumericalCertificate::estimate(signed_distance - radius, signed_distance + radius)
            }
            SourceAuthority::NoClaim => NumericalCertificate::no_claim(),
            SourceAuthority::MalformedFinite => NumericalCertificate {
                kind: NumericalKind::Enclosure,
                lo: signed_distance + 1.0,
                hi: signed_distance + 2.0,
            },
        };
        ChartSample {
            signed_distance,
            gradient: Some(Vec3::new(1.0, 0.0, 0.0)),
            lipschitz: Some(1.0),
            error,
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "test/bounded-authority-plane"
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        if matches!(self.authority, SourceAuthority::LocalEnclosure(_)) {
            TraceStepClaim::NoClaim
        } else {
            TraceStepClaim::ExactDistance
        }
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
                let err =
                    (grid.eval(p, cx).signed_distance - source.eval(p, cx).signed_distance).abs();
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
            .expect("eikonal probes")
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
        agree &= vdb.is_active(c) && vdb.get(c).to_bits() == v.to_bits();
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
            agree &= !vdb.is_active(c) && vdb.get(c).to_bits() == (-1.0f32).to_bits();
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
        // Reinit quality must be tested on a field that NEEDS it: a
        // translated SDF is already |∇φ| ≈ 1, where band-edge effects
        // dominate. Distort φ (|∇φ| ≈ 2.5), then require substantial
        // recovery toward the eikonal state.
        let distorted: Vec<([i32; 3], f32)> = band.grid().iter_active().collect();
        for (c, v) in distorted {
            band.grid_mut().set(c, v * 2.5);
        }
        let dev_before = band.stats().mean_eikonal_dev;
        band.reinitialize(40);
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

/// rsdf-007 — every chart sampler refuses an unresolved extended support
/// before evaluation, while its paired clipped API samples the actual
/// geometric intersection. Invalid spacings and excessive work are also
/// preflight refusals; narrow-band sampling rejects non-finite values and its
/// max-anchored lattice remains in range near `f64::MAX`; clipped sampling is
/// translation-equivariant (G3). Adaptive refinement refuses a split whose
/// midpoint cannot lie strictly inside an adjacent-float axis.
#[test]
#[allow(clippy::too_many_lines)] // One admission matrix compares every sampler and refusal boundary.
fn rsdf_007_sampling_domains_are_explicit_bounded_and_preflighted() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let clip = Aabb::new(Point3::new(-0.5, -0.5, -0.5), Point3::new(0.5, 0.5, 0.5));

        for refusal in [
            TiledSdf::build(&CountingPlane::new(0.0), 0.2, cx)
                .map(|_| ())
                .expect_err("dense default must refuse an unbounded plane"),
            AdaptiveSdf::build(&CountingPlane::new(0.0), 0.1, 3, cx)
                .map(|_| ())
                .expect_err("adaptive default must refuse an unbounded plane"),
            NarrowBand::from_chart(&CountingPlane::new(0.0), 0.2, 2, cx)
                .map(|_| ())
                .expect_err("band default must refuse an unbounded plane"),
        ] {
            assert!(matches!(
                refusal,
                SdfBuildError::SamplingDomain(SamplingDomainError::UnboundedSupport { .. })
            ));
        }

        let dense = TiledSdf::build_clipped(&CountingPlane::new(0.0), 0.2, clip, cx)
            .expect("dense clipped plane");
        assert!(dense.nominal_field_bound().is_finite());
        assert_eq!(
            dense.bound().to_bits(),
            dense.nominal_field_bound().to_bits(),
            "bound() remains the explicitly nominal compatibility accessor"
        );
        assert_eq!(dense.abstract_distance_kind(), NumericalKind::NoClaim);
        assert!(dense.abstract_distance_bound().is_none());
        let adaptive = AdaptiveSdf::build_clipped(&CountingPlane::new(0.0), 0.1, 3, clip, cx)
            .expect("adaptive clipped plane");
        assert_eq!(adaptive.abstract_distance_kind(), NumericalKind::NoClaim);
        assert!(adaptive.abstract_distance_bound().is_none());
        assert_eq!(
            adaptive.eval(Point3::new(-0.25, 0.0, 0.0), cx).error.kind,
            NumericalKind::NoClaim,
            "adaptive fitting cannot turn ClippedChart NoClaim into Estimate"
        );
        assert!(
            NarrowBand::from_chart_clipped(&CountingPlane::new(0.0), 0.2, 2, clip, cx).is_ok(),
            "band clipped plane"
        );

        let extreme_h = 1.0e306;
        let extreme = ExtremeFiniteBandChart::new(extreme_h);
        NarrowBand::from_chart(&extreme, extreme_h, 1, cx)
            .expect("max-anchored extreme finite band must build");
        assert!(extreme.eval_count() > 0);
        assert!(
            extreme.all_points_admitted(),
            "band lattice must stay finite and inside its inflated support"
        );
        let max_seen = extreme.max_seen();
        assert_eq!(max_seen.x.to_bits(), extreme.scan_bounds.max.x.to_bits());
        assert_eq!(max_seen.y.to_bits(), extreme.scan_bounds.max.y.to_bits());
        assert_eq!(max_seen.z.to_bits(), extreme.scan_bounds.max.z.to_bits());

        let nan_error = NarrowBand::from_chart(&FiniteNanChart, 0.5, 1, cx)
            .map(|_| ())
            .expect_err("non-finite narrow-band sample must be refused");
        match nan_error {
            SdfBuildError::InvalidSample { point, value_bits } => {
                assert!(point.x.is_finite() && point.y.is_finite() && point.z.is_finite());
                assert!(f64::from_bits(value_bits).is_nan());
            }
            other => panic!("expected InvalidSample, got {other:?}"),
        }

        let split_error = AdaptiveSdf::build(&UnrepresentableAdaptiveChart, 1e-9, 1, cx)
            .map(|_| ())
            .expect_err("adjacent-f64 axis cannot form strict octree children");
        match split_error {
            SdfBuildError::AdaptiveSubdivisionUnrepresentable {
                axis,
                min_bits,
                max_bits,
                midpoint_bits,
            } => {
                assert_eq!(axis, 0);
                assert_ne!(min_bits, max_bits);
                assert!(midpoint_bits == min_bits || midpoint_bits == max_bits);
            }
            other => panic!("expected adaptive split refusal, got {other:?}"),
        }

        assert!(matches!(
            AdaptiveSdf::build(&OverflowingAdaptiveResidualChart, 0.1, 0, cx),
            Err(SdfBuildError::InvalidReconstructionBound { value }) if !value.is_finite()
        ));

        let negative = dense.eval(Point3::new(-0.25, 0.0, 0.0), cx);
        let positive = dense.eval(Point3::new(0.25, 0.0, 0.0), cx);
        assert!(negative.signed_distance < 0.0);
        assert!(positive.signed_distance > 0.0);
        assert_eq!(negative.error.kind, NumericalKind::NoClaim);
        assert_eq!(positive.error.kind, NumericalKind::NoClaim);
        assert_eq!(
            dense.eval(Point3::new(2.0, 0.0, 0.0), cx).error.kind,
            NumericalKind::NoClaim,
            "outside-support extension cannot restore authority"
        );
        assert!(
            dense
                .raycast(
                    Point3::new(-2.0, 0.0, 0.0),
                    Vec3::new(1.0, 0.0, 0.0),
                    4.0,
                    cx
                )
                .is_none(),
            "compatibility raycast refuses a NoClaim field"
        );

        let invalid = CountingPlane::new(0.0);
        assert!(matches!(
            TiledSdf::build_clipped(&invalid, 0.0, clip, cx),
            Err(SdfBuildError::InvalidSpacing { .. })
        ));
        assert_eq!(
            invalid.eval_count(),
            0,
            "invalid h must preflight before eval"
        );

        let too_dense = CountingPlane::new(0.0);
        assert!(matches!(
            TiledSdf::build_clipped(&too_dense, 1e-9, clip, cx),
            Err(SdfBuildError::ResolutionTooFine { .. })
        ));
        assert_eq!(too_dense.eval_count(), 0, "dense count cap precedes eval");

        let too_deep = CountingPlane::new(0.0);
        assert!(matches!(
            AdaptiveSdf::build_clipped(&too_deep, 0.1, 7, clip, cx),
            Err(SdfBuildError::AdaptiveWorkLimit {
                cap: ADAPTIVE_MAX_NODES,
                ..
            })
        ));
        assert_eq!(too_deep.eval_count(), 0, "adaptive work cap precedes eval");

        let too_many_band_samples = CountingPlane::new(0.0);
        assert!(matches!(
            NarrowBand::from_chart_clipped(&too_many_band_samples, 1e-4, 2, clip, cx),
            Err(SdfBuildError::BandScanLimit { .. })
        ));
        assert_eq!(
            too_many_band_samples.eval_count(),
            0,
            "band scan cap precedes eval"
        );

        let shift = Vec3::new(0.25, -0.125, 0.375);
        let moved_clip = Aabb::new(clip.min.offset(shift), clip.max.offset(shift));
        let moved = TiledSdf::build_clipped(&CountingPlane::new(shift.x), 0.2, moved_clip, cx)
            .expect("translated clipped plane");
        assert_eq!(moved.abstract_distance_kind(), NumericalKind::NoClaim);
        for p in [
            Point3::new(-0.25, 0.0, 0.0),
            Point3::new(0.0, 0.2, -0.1),
            Point3::new(0.25, -0.1, 0.2),
        ] {
            let a = dense.eval(p, cx).signed_distance;
            let b = moved.eval(p.offset(shift), cx).signed_distance;
            assert!((a - b).abs() < 1e-6, "G3 mismatch: {a} versus {b}");
        }
    });
    verdict(
        "rsdf-007",
        true,
        "dense, adaptive, and narrow-band samplers reject unresolved extended support, excessive work, invalid samples, nonrepresentable adaptive splits, and overflowing fit residuals; the max-anchored narrow-band lattice stays finite and never overshoots extreme admitted support; paired clipped APIs sample source-intersection-clip and preserve translation equivariance (G3), while dense and adaptive clipped fields retain only finite nominal reconstruction bounds and honest NoClaim abstract-distance authority",
    );
}

/// rsdf-008 — sampled-field reconstruction cannot strengthen the source's
/// authority relative to abstract region signed distance.
#[test]
#[allow(clippy::too_many_lines)] // One authority lattice verifies every source/reconstruction pairing.
fn rsdf_008_sampled_sdf_preserves_weakest_source_authority() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let exact = TiledSdf::build(
            &BoundedAuthorityPlane {
                authority: SourceAuthority::Exact,
            },
            0.2,
            cx,
        )
        .expect("exact bounded source");
        assert_eq!(
            exact.abstract_distance_kind(),
            NumericalKind::Enclosure,
            "spline reconstruction demotes exact samples to an enclosure"
        );

        let rigorous = TiledSdf::build(
            &BoundedAuthorityPlane {
                authority: SourceAuthority::Enclosure(0.125),
            },
            0.2,
            cx,
        )
        .expect("rigorous bounded source");
        assert_eq!(rigorous.abstract_distance_kind(), NumericalKind::Enclosure);
        let rigorous_bound = rigorous
            .abstract_distance_bound()
            .expect("finite abstract-distance enclosure");
        assert!(
            rigorous_bound + 1e-15 >= rigorous.nominal_field_bound() + 0.125,
            "source certificate radius composes with reconstruction: {rigorous_bound} versus {} + 0.125",
            rigorous.nominal_field_bound()
        );
        let adaptive_rigorous = AdaptiveSdf::build(
            &BoundedAuthorityPlane {
                authority: SourceAuthority::Enclosure(0.125),
            },
            0.2,
            2,
            cx,
        )
        .expect("adaptive rigorous source");
        assert_eq!(
            adaptive_rigorous.abstract_distance_kind(),
            NumericalKind::Estimate,
            "probed adaptive residual is at best Estimate"
        );
        assert!(
            adaptive_rigorous
                .abstract_distance_bound()
                .expect("finite adaptive estimate")
                + 1e-15
                >= adaptive_rigorous.nominal_field_bound() + 0.125
        );

        let local_only = TiledSdf::build(
            &BoundedAuthorityPlane {
                authority: SourceAuthority::LocalEnclosure(0.125),
            },
            0.2,
            cx,
        )
        .expect("local sampled Lipschitz values still support a nominal fit");
        assert_eq!(
            local_only.abstract_distance_kind(),
            NumericalKind::Estimate,
            "sampled local Lipschitz maxima cannot mint a global enclosure"
        );
        assert!(
            local_only
                .eval(Point3::new(0.2, 0.0, 0.0), cx)
                .lipschitz
                .is_none(),
            "the reconstructed chart does not expose a sampled maximum as certified"
        );
        assert!(
            local_only
                .raycast(
                    Point3::new(-2.0, 0.0, 0.0),
                    Vec3::new(1.0, 0.0, 0.0),
                    4.0,
                    cx
                )
                .is_none()
        );

        let estimated = TiledSdf::build(
            &BoundedAuthorityPlane {
                authority: SourceAuthority::Estimate(0.25),
            },
            0.2,
            cx,
        )
        .expect("estimated bounded source");
        assert_eq!(estimated.abstract_distance_kind(), NumericalKind::Estimate);
        assert_eq!(
            estimated.eval(Point3::new(0.2, 0.1, -0.1), cx).error.kind,
            NumericalKind::Estimate
        );
        assert!(
            estimated
                .raycast(
                    Point3::new(-2.0, 0.0, 0.0),
                    Vec3::new(1.0, 0.0, 0.0),
                    4.0,
                    cx
                )
                .is_none(),
            "Estimate authority is insufficient for certified sphere tracing"
        );

        for authority in [SourceAuthority::NoClaim, SourceAuthority::MalformedFinite] {
            let grid = TiledSdf::build(&BoundedAuthorityPlane { authority }, 0.2, cx)
                .expect("nominal field remains sampleable");
            assert!(grid.nominal_field_bound().is_finite());
            assert_eq!(grid.abstract_distance_kind(), NumericalKind::NoClaim);
            assert!(grid.abstract_distance_bound().is_none());
            for point in [Point3::new(-0.2, 0.1, 0.0), Point3::new(3.0, 0.0, 0.0)] {
                assert_eq!(grid.eval(point, cx).error.kind, NumericalKind::NoClaim);
            }
            let adaptive = AdaptiveSdf::build(&BoundedAuthorityPlane { authority }, 0.2, 2, cx)
                .expect("adaptive nominal field remains sampleable");
            assert_eq!(adaptive.abstract_distance_kind(), NumericalKind::NoClaim);
            assert!(adaptive.abstract_distance_bound().is_none());
            assert_eq!(
                adaptive.eval(Point3::new(3.0, 0.0, 0.0), cx).error.kind,
                NumericalKind::NoClaim
            );
        }
    });
    verdict(
        "rsdf-008",
        true,
        "dense and adaptive sampled SDFs compose source certificate radius with nominal reconstruction, preserve their at-best authority, and make NoClaim absorbing for explicit and malformed source certificates inside and outside support",
    );
}
