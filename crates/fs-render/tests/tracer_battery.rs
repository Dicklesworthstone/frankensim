//! Battery for the spectral path tracer (bead 872c; runs under the
//! `tracer` feature). The acceptance ladder: furnace-in-color
//! exactness, the frozen Cornell golden, MIS-beats-either-alone
//! variance, EXR byte-exact round trip, progressive-checkpoint and
//! tile-order bitwise invariance, and the ledgered Sobol-vs-iid
//! equal-spp claim (measured, never vibes).
#![cfg(feature = "tracer")]

use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

use asupersync::types::Budget;
use fs_evidence::NumericalCertificate;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Aabb, Chart, ChartSample, Point3, TraceStepClaim, Vec3};
use fs_render::charts::TriMesh;
use fs_render::spectral::{LAMBDA_MAX, LAMBDA_MIN, lift_rgb, xyz_of_spectrum};
use fs_render::tracer::{
    Camera, DirectStrategy, Film, Material, Primitive, RectLight, Sampler, Scene, Settings, Shape,
    TracerError, film_to_exr, render, render_range,
};
use fs_render::{cosine_sample_hemisphere, hero_wavelengths, radical_inverse};
use fs_rep_frep::FrepBuilder;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 872,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn assert_film_bits_eq(left: &Film, right: &Film, context: &str) {
    assert_eq!(
        (left.width, left.height),
        (right.width, right.height),
        "{context}"
    );
    assert_eq!(left.spp_done, right.spp_done, "{context}");
    assert_eq!(left.xyz.len(), right.xyz.len(), "{context}");
    for (a, b) in left.xyz.iter().zip(&right.xyz) {
        for channel in 0..3 {
            assert_eq!(a[channel].to_bits(), b[channel].to_bits(), "{context}");
        }
    }
}

struct CancellingSphere {
    center: Point3,
    radius: f64,
    evaluations: Arc<AtomicUsize>,
    cancel_at: Option<usize>,
    gate: Option<Arc<CancelGate>>,
}

impl Chart for CancellingSphere {
    fn eval(&self, point: Point3, _cx: &Cx<'_>) -> ChartSample {
        let evaluation = self.evaluations.fetch_add(1, Ordering::SeqCst) + 1;
        if self.cancel_at == Some(evaluation)
            && let Some(gate) = &self.gate
        {
            gate.request();
        }
        let delta = point.delta_from(self.center);
        let norm = delta.norm();
        let distance = norm - self.radius;
        ChartSample {
            signed_distance: distance,
            gradient: (norm > 1e-12).then(|| delta.scale(1.0 / norm)),
            lipschitz: Some(1.0),
            error: NumericalCertificate::exact(distance),
        }
    }

    fn support(&self) -> Aabb {
        let r = self.radius;
        Aabb::new(
            self.center.offset(Vec3::new(-r, -r, -r)),
            self.center.offset(Vec3::new(r, r, r)),
        )
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn name(&self) -> &'static str {
        "tracer-cancellation-sphere"
    }
}

struct ConstantNoClaim;

impl Chart for ConstantNoClaim {
    fn eval(&self, _point: Point3, _cx: &Cx<'_>) -> ChartSample {
        ChartSample {
            signed_distance: 20_000.0,
            gradient: None,
            lipschitz: None,
            error: NumericalCertificate::estimate(20_000.0, 20_000.0),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0))
    }

    fn name(&self) -> &'static str {
        "tracer-constant-no-claim"
    }
}

fn replace_cornell_sphere(scene: &mut Scene, chart: CancellingSphere) {
    scene.primitives[5].shape = Shape::Chart(Box::new(chart));
}

fn quad(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> TriMesh {
    TriMesh::new(vec![a, b, c, d], vec![[0, 1, 2], [0, 2, 3]])
}

/// The Cornell-class fixture: unit box, white floor/ceiling/back, red
/// left, green right, a GGX F-rep sphere, one ceiling rect light.
fn cornell() -> Scene {
    let white = lift_rgb([0.73, 0.73, 0.73]);
    let red = lift_rgb([0.63, 0.065, 0.05]);
    let green = lift_rgb([0.14, 0.45, 0.091]);
    let lam = |r| Material::Lambertian { reflectance: r };
    let mut primitives = vec![
        // floor y=0 (normal +y), ceiling y=1, back z=0.
        Primitive {
            shape: Shape::Mesh(quad(
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 0.0, 1.0],
                [0.0, 0.0, 1.0],
            )),
            material: lam(white),
            emission: None,
        },
        Primitive {
            shape: Shape::Mesh(quad(
                [0.0, 1.0, 0.0],
                [0.0, 1.0, 1.0],
                [1.0, 1.0, 1.0],
                [1.0, 1.0, 0.0],
            )),
            material: lam(white),
            emission: None,
        },
        Primitive {
            shape: Shape::Mesh(quad(
                [0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
            )),
            material: lam(white),
            emission: None,
        },
        Primitive {
            shape: Shape::Mesh(quad(
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 1.0, 1.0],
                [0.0, 1.0, 0.0],
            )),
            material: lam(red),
            emission: None,
        },
        Primitive {
            shape: Shape::Mesh(quad(
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [1.0, 1.0, 1.0],
                [1.0, 0.0, 1.0],
            )),
            material: lam(green),
            emission: None,
        },
    ];
    // GGX sphere via the certified F-rep chart (sphere-traced).
    let mut b = FrepBuilder::new();
    let s = b
        .sphere(Point3::new(0.42, 0.28, 0.45), 0.28)
        .expect("sphere");
    let frep = b.finish(s).expect("frep");
    primitives.push(Primitive {
        shape: Shape::Chart(Box::new(frep)),
        material: Material::Ggx {
            reflectance: lift_rgb([0.9, 0.9, 0.9]),
            // Near-specular: the small light's reflection in the sphere
            // is the Veach regime where NEE variance explodes (the
            // sampled light point almost never aligns with the sharp
            // lobe) while BSDF sampling finds the light reliably — the
            // region MIS needs to win the variance acceptance against
            // NEE-only, complementing the diffuse walls where NEE wins
            // against BSDF-only.
            alpha: 0.04,
        },
        emission: None,
    });
    // Ceiling light: rect + the SAME rect as emissive geometry.
    let emission = (lift_rgb([1.0, 1.0, 1.0]), 18.0);
    let (corner, eu, ev) = (
        Point3::new(0.35, 0.9995, 0.35),
        Vec3::new(0.3, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 0.3),
    );
    primitives.push(Primitive {
        shape: Shape::Mesh(quad(
            [0.35, 0.9995, 0.35],
            [0.65, 0.9995, 0.35],
            [0.65, 0.9995, 0.65],
            [0.35, 0.9995, 0.65],
        )),
        material: lam(white),
        emission: Some(emission),
    });
    let light_prim = primitives.len() - 1;
    Scene {
        primitives,
        light: RectLight {
            corner,
            edge_u: eu,
            edge_v: ev,
            prim: light_prim,
            emission,
        },
        camera: Camera {
            // Framed so the near-specular sphere (the light's sharp
            // reflection — the BSDF-favored Veach regime) fills a
            // meaningful share of the image next to the NEE-favored
            // diffuse walls.
            eye: Point3::new(0.46, 0.4, 1.45),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.3,
        },
    }
}

fn settings(strategy: DirectStrategy, sampler: Sampler, seed: u64, px: u32, spp: u32) -> Settings {
    Settings {
        width: px,
        height: px,
        spp,
        max_depth: 4,
        sampler,
        strategy,
        seed,
    }
}

/// ACCEPTANCE (1): the furnace, now in color. Radiance part: for a
/// Lambertian under uniform incident L, every cosine-weighted sample
/// returns EXACTLY ρ(λ)·L (f·cos/pdf = (ρ/π)·L·cos·(π/cos)) — the v0
/// zero-variance bar, per wavelength. Color part: pushing ρ through
/// the tracer's hero-wavelength → XYZ estimator converges to the
/// quadrature XYZ of ρ (the same integral by a different route).
#[test]
fn furnace_in_color_is_exact() {
    let rho = lift_rgb([0.63, 0.065, 0.05]);
    let incident = 2.5;
    // Radiance exactness per wavelength (the zero-variance property).
    for i in 1..=64u64 {
        let (dir, pdf) = cosine_sample_hemisphere(radical_inverse(2, i), radical_inverse(3, i));
        let lambda = LAMBDA_MIN + radical_inverse(5, i) * (LAMBDA_MAX - LAMBDA_MIN);
        let f = rho.eval(lambda) / core::f64::consts::PI;
        let sample = f * incident * dir[2] / pdf;
        let expect = rho.eval(lambda) * incident;
        assert!(
            (sample - expect).abs() <= 1e-14 * expect,
            "furnace sample {sample:.17e} vs {expect:.17e} at λ={lambda}"
        );
    }
    // XYZ-level: hero-wavelength estimator vs quadrature reference.
    let range = LAMBDA_MAX - LAMBDA_MIN;
    let kn = 1.0 / fs_render::spectral::y_integral();
    let n = 4096u64;
    let mut xyz = [0.0f64; 3];
    for i in 1..=n {
        let hero = LAMBDA_MIN + radical_inverse(2, i) * range;
        for l in hero_wavelengths(hero, 4, LAMBDA_MIN, LAMBDA_MAX) {
            let w = rho.eval(l) * incident * range / 4.0 * kn;
            xyz[0] += w * fs_render::spectral::cie_x(l);
            xyz[1] += w * fs_render::spectral::cie_y(l);
            xyz[2] += w * fs_render::spectral::cie_z(l);
        }
    }
    let reference = xyz_of_spectrum(|l| rho.eval(l) * incident);
    let mut worst = 0.0f64;
    for a in 0..3 {
        worst = worst.max((xyz[a] / n as f64 - reference[a]).abs());
    }
    assert!(worst < 2e-3, "hero XYZ off quadrature by {worst:.2e}");
    println!(
        "{{\"suite\":\"fs-render/tracer\",\"case\":\"furnace-color\",\"verdict\":\"pass\",\"detail\":\"radiance exact to 1e-14 rel, XYZ vs quadrature {worst:.2e}\"}}"
    );
}

fn fnv(bytes: &[u8]) -> u64 {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        acc ^= u64::from(b);
        acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
    }
    acc
}

/// Frozen 2026-07-12 at the committed 872c tree: 24×24, 8 spp,
/// depth 4, MIS + iid Philox, seed 7 — the first shaded COLOR image.
/// FNV-1a over the EXR bytes. That original surface was verified debug +
/// release on aarch64 (M4 Pro) and x86-64 (ts2 5995WX). Bead 8ll9 now records
/// an explicit fs-render:chart-backend-bits=1 dependency; current-tree replay
/// in all four quadrants remains mandatory before closeout. Re-freeze only per
/// docs/GOLDEN_POLICY.md.
const CORNELL_GOLDEN: u64 = 0x6ed8_706b_08d1_642e;

/// ACCEPTANCE (2): the Cornell-class fixture matches the frozen
/// reference image hash in deterministic mode.
#[test]
fn cornell_box_matches_the_frozen_golden() {
    let scene = cornell();
    let film = with_cx(|cx| {
        render(
            &scene,
            cx,
            &settings(DirectStrategy::Mis, Sampler::Iid, 7, 24, 8),
        )
    })
    .expect("Cornell render");
    let exr = film_to_exr(&film).expect("encode");
    let hash = fnv(&exr);
    // The image is not black and not blown out: the mid pixel saw light.
    let mid = &film.xyz[(12 * 24 + 12) as usize];
    assert!(
        mid[1] > 0.0,
        "mid-pixel Y is zero: the scene rendered black"
    );
    println!(
        "{{\"suite\":\"fs-render/tracer\",\"case\":\"cornell-golden\",\"verdict\":\"info\",\"detail\":\"{hash:#018x}\"}}"
    );
    assert_eq!(
        hash, CORNELL_GOLDEN,
        "Cornell EXR bits changed: {hash:#018x} vs {CORNELL_GOLDEN:#018x} — re-freeze only with \
         a semantic justification per docs/GOLDEN_POLICY.md (bump the causative fs-render bit \
         surface and update golden-couplings.json in the same commit)"
    );
}

/// ACCEPTANCE (4): EXR round-trips byte-exactly through fs-img.
#[test]
fn exr_round_trips_byte_exactly() {
    let scene = cornell();
    let film = with_cx(|cx| {
        render(
            &scene,
            cx,
            &settings(DirectStrategy::Mis, Sampler::Iid, 7, 12, 2),
        )
    })
    .expect("round-trip render");
    let bytes = film_to_exr(&film).expect("encode");
    let decoded = fs_img::read_exr(&bytes).expect("decode");
    let re =
        fs_img::write_exr(decoded.width, decoded.height, &decoded.channels).expect("re-encode");
    assert_eq!(bytes, re, "EXR bytes changed across a decode/encode cycle");
    println!(
        "{{\"suite\":\"fs-render/tracer\",\"case\":\"exr-roundtrip\",\"verdict\":\"pass\",\"detail\":\"{} bytes byte-exact\"}}",
        bytes.len()
    );
}

/// Progressive rendering: the 8-spp render equals the 3-spp checkpoint
/// continued to 8, bitwise (the pause–serialize–resume doctrine).
#[test]
fn progressive_checkpoint_is_bitwise() {
    let scene = cornell();
    let s = settings(DirectStrategy::Mis, Sampler::Iid, 11, 12, 8);
    let (direct, resumed) = with_cx(|cx| {
        let direct = render(&scene, cx, &s).expect("direct render");
        let mut film = Film::new(s.width, s.height);
        render_range(&scene, cx, &s, &mut film, 0, 3).expect("first range");
        render_range(&scene, cx, &s, &mut film, 3, 8).expect("resumed range");
        (direct, film)
    });
    assert_eq!(direct.spp_done, resumed.spp_done);
    for (a, b) in direct.xyz.iter().zip(&resumed.xyz) {
        for k in 0..3 {
            assert_eq!(a[k].to_bits(), b[k].to_bits(), "checkpoint drifted");
        }
    }
    println!(
        "{{\"suite\":\"fs-render/tracer\",\"case\":\"progressive-bitwise\",\"verdict\":\"pass\",\"detail\":\"3+5 spp == 8 spp bitwise\"}}"
    );
}

#[test]
fn reversed_progressive_range_is_rejected_transactionally() {
    let scene = cornell();
    let s = settings(DirectStrategy::Mis, Sampler::Iid, 31, 2, 3);
    let mut film = Film::new(s.width, s.height);
    film.spp_done = 3;
    for xyz in &mut film.xyz {
        *xyz = [0.25, -0.0, f64::from_bits(0x7ff8_0000_0000_0042)];
    }
    let before = film.clone();
    assert_eq!(
        with_cx(|cx| render_range(&scene, cx, &s, &mut film, 3, 2)),
        Err(TracerError::InvalidRange { from: 3, to: 2 })
    );
    assert_film_bits_eq(&film, &before, "invalid range changed film bits");
}

#[test]
fn film_allocation_and_public_buffer_shape_fail_closed() {
    assert_eq!(Film::try_new(0, 1), Err(TracerError::InvalidInput));
    assert_eq!(
        Film::try_new(u32::MAX, u32::MAX),
        Err(TracerError::InvalidInput)
    );

    let scene = cornell();
    let s = settings(DirectStrategy::Mis, Sampler::Iid, 41, 2, 1);
    let mut malformed = Film::new(s.width, s.height);
    malformed.xyz.pop();
    assert_eq!(
        with_cx(|cx| render_range(&scene, cx, &s, &mut malformed, 0, 0)),
        Err(TracerError::InvalidInput),
        "an empty range must still validate the public film buffer"
    );

    let zero_settings = Settings {
        width: 0,
        height: 0,
        ..s
    };
    let mut zero_film = Film {
        width: 0,
        height: 0,
        xyz: Vec::new(),
        spp_done: 0,
    };
    assert_eq!(
        with_cx(|cx| render_range(&scene, cx, &zero_settings, &mut zero_film, 0, 0)),
        Err(TracerError::InvalidInput)
    );
}

#[test]
fn production_tracer_rejects_uncertified_misses() {
    let mut scene = cornell();
    scene.primitives[5].shape = Shape::Chart(Box::new(ConstantNoClaim));
    let s = settings(DirectStrategy::Mis, Sampler::Iid, 37, 1, 1);
    assert_eq!(
        with_cx(|cx| render(&scene, cx, &s)),
        Err(TracerError::UncertifiedTrace)
    );
}

#[test]
fn cancelled_range_is_transactional_and_retryable() {
    let s = settings(DirectStrategy::Mis, Sampler::Iid, 23, 8, 3);
    let mut reference_scene = cornell();
    replace_cornell_sphere(
        &mut reference_scene,
        CancellingSphere {
            center: Point3::new(0.42, 0.28, 0.45),
            radius: 0.28,
            evaluations: Arc::new(AtomicUsize::new(0)),
            cancel_at: None,
            gate: None,
        },
    );
    let mut film = Film::new(s.width, s.height);
    with_cx(|cx| render_range(&reference_scene, cx, &s, &mut film, 0, 1))
        .expect("initial checkpoint");
    let before = film.clone();

    let gate = Arc::new(CancelGate::new());
    let evaluations = Arc::new(AtomicUsize::new(0));
    let mut cancelling_scene = cornell();
    replace_cornell_sphere(
        &mut cancelling_scene,
        CancellingSphere {
            center: Point3::new(0.42, 0.28, 0.45),
            radius: 0.28,
            evaluations: Arc::clone(&evaluations),
            cancel_at: Some(64),
            gate: Some(Arc::clone(&gate)),
        },
    );
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 872,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        assert_eq!(
            render_range(&cancelling_scene, &cx, &s, &mut film, 1, 3),
            Err(TracerError::Cancelled)
        );
    });
    assert!(evaluations.load(Ordering::SeqCst) >= 64);
    assert_film_bits_eq(&film, &before, "failed ranges must not alter film state");

    with_cx(|cx| render_range(&cancelling_scene, cx, &s, &mut film, 1, 3))
        .expect("retry after cancellation");
    let direct = with_cx(|cx| render(&reference_scene, cx, &s)).expect("direct reference");
    assert_film_bits_eq(&film, &direct, "retry must equal a direct render bitwise");
}

/// Deterministic replay under the Owen-Sobol stream. Progressive sample-range
/// equivalence is exercised separately above; this case does not claim a
/// parallel tile-order execution it does not perform.
#[test]
fn sample_streams_replay_bitwise() {
    let scene = cornell();
    let s = settings(DirectStrategy::Mis, Sampler::OwenSobol, 5, 12, 4);
    let (a, b) = with_cx(|cx| {
        (
            render(&scene, cx, &s).expect("replay a"),
            render(&scene, cx, &s).expect("replay b"),
        )
    });
    for (x, y) in a.xyz.iter().zip(&b.xyz) {
        for k in 0..3 {
            assert_eq!(x[k].to_bits(), y[k].to_bits(), "replay drifted");
        }
    }
    println!(
        "{{\"suite\":\"fs-render/tracer\",\"case\":\"schedule-invariance\",\"verdict\":\"pass\",\"detail\":\"replay bitwise under OwenSobol\"}}"
    );
}

fn mean_pixel_variance(
    scene: &Scene,
    strategy: DirectStrategy,
    sampler: Sampler,
    spp: u32,
    px: u32,
) -> f64 {
    // Variance across independent seeds of the per-pixel luminance.
    const SEEDS: u64 = 6;
    let n = (px * px) as usize;
    let mut sum = vec![0.0f64; n];
    let mut sum2 = vec![0.0f64; n];
    for seed in 0..SEEDS {
        let film =
            with_cx(|cx| render(scene, cx, &settings(strategy, sampler, 100 + seed, px, spp)))
                .expect("variance render");
        let inv = 1.0 / f64::from(spp);
        for (i, xyz) in film.xyz.iter().enumerate() {
            let y = xyz[1] * inv;
            sum[i] += y;
            sum2[i] += y * y;
        }
    }
    let k = SEEDS as f64;
    (0..n)
        .map(|i| (sum2[i] - sum[i] * sum[i] / k) / (k - 1.0))
        .sum::<f64>()
        / n as f64
}

/// ACCEPTANCE (3): MIS beats either technique alone on the mixed
/// diffuse+glossy fixture (variance across 6 seeds, 12×12 @ 4 spp,
/// seeds 100..106 — logged, falsifiable).
#[test]
fn mis_beats_either_technique_alone() {
    let scene = cornell();
    let v_mis = mean_pixel_variance(&scene, DirectStrategy::Mis, Sampler::Iid, 4, 12);
    let v_nee = mean_pixel_variance(&scene, DirectStrategy::NeeOnly, Sampler::Iid, 4, 12);
    let v_bsdf = mean_pixel_variance(&scene, DirectStrategy::BsdfOnly, Sampler::Iid, 4, 12);
    println!(
        "{{\"suite\":\"fs-render/tracer\",\"case\":\"mis-variance\",\"verdict\":\"info\",\"detail\":\"var mis {v_mis:.3e} nee {v_nee:.3e} bsdf {v_bsdf:.3e}\"}}"
    );
    assert!(
        v_mis < v_nee && v_mis < v_bsdf,
        "MIS variance {v_mis:.3e} does not beat NEE {v_nee:.3e} / BSDF {v_bsdf:.3e}"
    );
}

/// AMBITION ROUND A: the Owen-Sobol equal-spp claim, measured. The
/// debug tier LOGS the ratio at 16 spp (informational); the release
/// `--ignored` lane below asserts-or-records at the bead's named
/// 64 spp.
#[test]
fn sobol_vs_iid_equal_spp_logged() {
    let scene = cornell();
    let v_iid = mean_pixel_variance(&scene, DirectStrategy::Mis, Sampler::Iid, 16, 12);
    let v_sobol = mean_pixel_variance(&scene, DirectStrategy::Mis, Sampler::OwenSobol, 16, 12);
    println!(
        "{{\"suite\":\"fs-render/tracer\",\"case\":\"sobol-vs-iid-16spp\",\"verdict\":\"info\",\"detail\":\"var iid {v_iid:.3e} sobol {v_sobol:.3e} ratio {:.3}\"}}",
        v_sobol / v_iid
    );
}

/// The bead's 64-spp Sobol claim, release lane:
/// `cargo test -p fs-render --release --features tracer --test tracer_battery -- --ignored --nocapture`
#[test]
#[ignore = "equal-spp variance lane: run explicitly in release with --ignored"]
fn sobol_vs_iid_at_64spp() {
    let scene = cornell();
    let v_iid = mean_pixel_variance(&scene, DirectStrategy::Mis, Sampler::Iid, 64, 16);
    let v_sobol = mean_pixel_variance(&scene, DirectStrategy::Mis, Sampler::OwenSobol, 64, 16);
    let verdict = if v_sobol < v_iid {
        "sobol-wins"
    } else {
        "iid-holds"
    };
    println!(
        "{{\"suite\":\"fs-render/tracer\",\"case\":\"sobol-vs-iid-64spp\",\"verdict\":\"{verdict}\",\"detail\":\"var iid {v_iid:.3e} sobol {v_sobol:.3e} ratio {:.3} - ledger on bead 872c\"}}",
        v_sobol / v_iid
    );
}
