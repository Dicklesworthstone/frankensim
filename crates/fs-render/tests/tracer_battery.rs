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
use fs_blake3::hash_domain;
use fs_evidence::NumericalCertificate;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Aabb, Chart, ChartSample, Point3, TraceStepClaim, Vec3};
use fs_render::animated_instances::{
    AnimatedGeometryInstance, RigidTransformTrajectory, TransformKeyframe,
};
use fs_render::camera::{
    AnimatedCamera, Aperture, CameraError, CameraKeyframe, CameraShot, CutSide, PhysicalCamera,
};
use fs_render::charts::TriMesh;
use fs_render::instances::{GeometryInstance, RigidTransform, SharedGeometry};
use fs_render::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution, ShutterInterval};
use fs_render::motion_vectors::StableFeatureIdentity;
use fs_render::spectral::{LAMBDA_MAX, LAMBDA_MIN, lift_rgb, xyz_of_spectrum};
use fs_render::tracer::{
    Camera, DirectStrategy, Film, FilmTimeMode, Material, Primitive, RectLight, Sampler, Scene,
    Settings, Shape, TracerError, film_to_exr, render, render_cinematic, render_cinematic_range,
    render_motion, render_motion_range, render_range, trace_cinematic_pixel_sample,
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

fn assert_film_state_bits_eq(left: &Film, right: &Film, context: &str) {
    assert_film_bits_eq(left, right, context);
    assert_eq!(left.time_mode, right.time_mode, "{context}");
}

struct CancellingSphere {
    center: Point3,
    radius: f64,
    evaluations: Arc<AtomicUsize>,
    cancel_at: Option<usize>,
    gate: Option<Arc<CancelGate>>,
}

impl Chart for CancellingSphere {
    fn eval(&self, point: Point3, cx: &Cx<'_>) -> ChartSample {
        let evaluation = self.evaluations.fetch_add(1, Ordering::SeqCst) + 1;
        if self.cancel_at == Some(evaluation)
            && let Some(gate) = &self.gate
        {
            gate.request();
        }
        SphereChart {
            center: self.center,
            radius: self.radius,
        }
        .eval(point, cx)
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
        lights: vec![RectLight {
            corner,
            edge_u: eu,
            edge_v: ev,
            prim: light_prim,
            emission,
        }],
        environment: None,
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

fn motion_shutter(frame_time_s: f64, duration_s: f64) -> ShutterInterval {
    ShutterInterval::resolve(
        frame_time_s,
        duration_s,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::StratifiedCounterV1 { strata: 4_096 },
        ShotTimeBounds::try_new(0.0, 1.0).expect("motion-fixture shot bounds"),
    )
    .expect("motion-fixture shutter")
}

fn emissive_motion_scene(animated: bool) -> Scene {
    let local_mesh = quad(
        [-0.5, -2.0, 0.0],
        [0.5, -2.0, 0.0],
        [0.5, 2.0, 0.0],
        [-0.5, 2.0, 0.0],
    );
    let geometry = SharedGeometry::mesh(local_mesh);
    let geometry_identity = hash_domain(
        "org.frankensim.test.motion-blur-emissive-quad.v1",
        b"unit-width-quad",
    );
    let shape = if animated {
        let start = TransformKeyframe::try_new(
            0.0,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [-1.0, 0.0, 0.0])
                .expect("start transform"),
            [2.0, 0.0, 0.0],
        )
        .expect("start keyframe");
        let end = TransformKeyframe::try_new(
            1.0,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).expect("end transform"),
            [2.0, 0.0, 0.0],
        )
        .expect("end keyframe");
        let trajectory =
            RigidTransformTrajectory::try_new(vec![start, end]).expect("linear trajectory");
        Shape::AnimatedInstance(
            AnimatedGeometryInstance::try_new(101, geometry_identity, geometry, trajectory)
                .expect("animated instance"),
        )
    } else {
        Shape::Instance(
            GeometryInstance::try_new(101, geometry_identity, geometry, RigidTransform::identity())
                .expect("static instance"),
        )
    };
    emissive_instance_scene(shape, Point3::new(0.0, 0.0, 2.0))
}

fn absolute_time_translation_scene(animated: bool, start_s: f64, end_s: f64) -> Scene {
    let half_width = 0.125;
    let local_mesh = quad(
        [-half_width, -2.0, 0.0],
        [half_width, -2.0, 0.0],
        [half_width, 2.0, 0.0],
        [-half_width, 2.0, 0.0],
    );
    let geometry = SharedGeometry::mesh(local_mesh);
    let geometry_identity = hash_domain(
        "org.frankensim.test.absolute-time-translation.v1",
        &half_width.to_bits().to_le_bytes(),
    );
    let shape = if animated {
        let velocity_x = 2.0 / (end_s - start_s);
        let start = TransformKeyframe::try_new(
            start_s,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [-1.0, 0.0, 0.0])
                .expect("absolute-time start transform"),
            [velocity_x, 0.0, 0.0],
        )
        .expect("absolute-time start keyframe");
        let end = TransformKeyframe::try_new(
            end_s,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [1.0, 0.0, 0.0])
                .expect("absolute-time end transform"),
            [velocity_x, 0.0, 0.0],
        )
        .expect("absolute-time end keyframe");
        let trajectory =
            RigidTransformTrajectory::try_new(vec![start, end]).expect("absolute-time trajectory");
        Shape::AnimatedInstance(
            AnimatedGeometryInstance::try_new(103, geometry_identity, geometry, trajectory)
                .expect("absolute-time animated instance"),
        )
    } else {
        Shape::Instance(
            GeometryInstance::try_new(103, geometry_identity, geometry, RigidTransform::identity())
                .expect("absolute-time static instance"),
        )
    };
    emissive_instance_scene(shape, Point3::new(0.0, 0.0, 2.0))
}

fn emissive_spin_scene(animated: bool) -> Scene {
    let local_mesh = quad(
        [-1.0, -0.1, 0.0],
        [1.0, -0.1, 0.0],
        [1.0, 0.1, 0.0],
        [-1.0, 0.1, 0.0],
    );
    let geometry = SharedGeometry::mesh(local_mesh);
    let geometry_identity = hash_domain(
        "org.frankensim.test.motion-blur-emissive-spinner.v1",
        b"thin-spinning-quad",
    );
    let shape = if animated {
        let start = TransformKeyframe::try_new(0.0, RigidTransform::identity(), [0.0; 3])
            .expect("spin start keyframe");
        let half_angle = core::f64::consts::FRAC_PI_4;
        let end = TransformKeyframe::try_new(
            1.0,
            RigidTransform::try_new([0.0, 0.0, half_angle.sin(), half_angle.cos()], [0.0; 3])
                .expect("spin end transform"),
            [0.0; 3],
        )
        .expect("spin end keyframe");
        let trajectory =
            RigidTransformTrajectory::try_new(vec![start, end]).expect("constant spin trajectory");
        Shape::AnimatedInstance(
            AnimatedGeometryInstance::try_new(102, geometry_identity, geometry, trajectory)
                .expect("spinning instance"),
        )
    } else {
        Shape::Instance(
            GeometryInstance::try_new(102, geometry_identity, geometry, RigidTransform::identity())
                .expect("static spinner"),
        )
    };
    emissive_instance_scene(shape, Point3::new(0.75, 0.0, 2.0))
}

fn emissive_instance_scene(shape: Shape, eye: Point3) -> Scene {
    let white = lift_rgb([0.8, 0.8, 0.8]);
    Scene {
        primitives: vec![
            Primitive {
                shape,
                material: Material::Lambertian { reflectance: white },
                emission: Some((white, 1.0)),
            },
            Primitive {
                shape: Shape::Mesh(quad(
                    [10.0, -0.5, 0.0],
                    [11.0, -0.5, 0.0],
                    [11.0, 0.5, 0.0],
                    [10.0, 0.5, 0.0],
                )),
                material: Material::Lambertian { reflectance: white },
                emission: Some((white, 1.0)),
            },
        ],
        lights: vec![RectLight {
            corner: Point3::new(10.0, -0.5, 0.0),
            edge_u: Vec3::new(1.0, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 1.0, 0.0),
            prim: 1,
            emission: (white, 1.0),
        }],
        environment: None,
        camera: Camera {
            eye,
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.0,
        },
    }
}

fn motion_settings(spp: u32) -> Settings {
    Settings {
        width: 1,
        height: 1,
        spp,
        max_depth: 1,
        sampler: Sampler::Iid,
        strategy: DirectStrategy::Mis,
        seed: 0x6d6f_7469_6f6e,
    }
}

fn physical_from_legacy(
    camera: &Camera,
    focus_distance_m: f64,
    aperture: Aperture,
) -> PhysicalCamera {
    PhysicalCamera::try_legacy_compatible(
        camera.eye,
        camera.forward,
        camera.up,
        camera.half_tan,
        focus_distance_m,
        aperture,
    )
    .expect("valid legacy-equivalent physical camera")
}

fn static_cinematic_camera(camera: PhysicalCamera) -> AnimatedCamera {
    AnimatedCamera::try_static(701, 0.0, 1.0, camera).expect("static cinematic camera")
}

/// Directly visible emissive cards at two depths make lens-origin changes
/// observable without conflating the camera test with indirect-light noise.
fn depth_varying_emissive_scene() -> Scene {
    let near = lift_rgb([1.0, 0.12, 0.04]);
    let far = lift_rgb([0.04, 0.3, 1.0]);
    let white = lift_rgb([1.0, 1.0, 1.0]);
    Scene {
        primitives: vec![
            Primitive {
                shape: Shape::Mesh(quad(
                    [-0.32, -0.32, 1.0],
                    [0.32, -0.32, 1.0],
                    [0.32, 0.32, 1.0],
                    [-0.32, 0.32, 1.0],
                )),
                material: Material::Lambertian { reflectance: near },
                emission: Some((near, 2.0)),
            },
            Primitive {
                shape: Shape::Mesh(quad(
                    [-0.9, -0.9, 0.0],
                    [0.9, -0.9, 0.0],
                    [0.9, 0.9, 0.0],
                    [-0.9, 0.9, 0.0],
                )),
                material: Material::Lambertian { reflectance: far },
                emission: Some((far, 1.0)),
            },
            Primitive {
                shape: Shape::Mesh(quad(
                    [10.0, -0.5, 0.0],
                    [11.0, -0.5, 0.0],
                    [11.0, 0.5, 0.0],
                    [10.0, 0.5, 0.0],
                )),
                material: Material::Lambertian { reflectance: white },
                emission: Some((white, 1.0)),
            },
        ],
        lights: vec![RectLight {
            corner: Point3::new(10.0, -0.5, 0.0),
            edge_u: Vec3::new(1.0, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 1.0, 0.0),
            prim: 2,
            emission: (white, 1.0),
        }],
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 2.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.5,
        },
    }
}

/// A centered emissive square used to measure ideal thin-lens occupancy.
/// `half_extent` is in world metres and `target_z` selects focused (`0`) or
/// halfway-to-camera (`1`) placement for the canonical two-metre focus.
fn focus_probe_scene(target_z: f64, half_extent: f64) -> Scene {
    let white = lift_rgb([1.0, 1.0, 1.0]);
    Scene {
        primitives: vec![
            Primitive {
                shape: Shape::Mesh(quad(
                    [-half_extent, -half_extent, target_z],
                    [half_extent, -half_extent, target_z],
                    [half_extent, half_extent, target_z],
                    [-half_extent, half_extent, target_z],
                )),
                material: Material::Lambertian { reflectance: white },
                emission: Some((white, 1.0)),
            },
            Primitive {
                shape: Shape::Mesh(quad(
                    [10.0, -0.5, 0.0],
                    [11.0, -0.5, 0.0],
                    [11.0, 0.5, 0.0],
                    [10.0, 0.5, 0.0],
                )),
                material: Material::Lambertian { reflectance: white },
                emission: Some((white, 1.0)),
            },
        ],
        lights: vec![RectLight {
            corner: Point3::new(10.0, -0.5, 0.0),
            edge_u: Vec3::new(1.0, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 1.0, 0.0),
            prim: 1,
            emission: (white, 1.0),
        }],
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 2.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.0,
        },
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

/// Re-frozen 2026-08-06 for tracer bit-semantics v2: reflective GGX moved
/// from NDF sampling to a matched isotropic visible-normal sampler. The BSDF
/// integral is unchanged in expectation, but deterministic sample directions,
/// PDFs, MIS weights, image bits, and variance deliberately changed. Fixture:
/// 24×24, 8 spp, depth 4, MIS + iid Philox, seed 7; FNV-1a over the EXR bytes.
/// This freeze also records the current chart-backend-bits=11 dependency.
/// Re-freeze only per docs/GOLDEN_POLICY.md.
const CORNELL_GOLDEN: u64 = 0xe89b_e51c_b59b_21cc;

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
    assert_film_state_bits_eq(&film, &before, "invalid range changed film bits");
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
        time_mode: FilmTimeMode::Uninitialized,
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
    assert_film_state_bits_eq(&film, &before, "failed ranges must not alter film state");

    with_cx(|cx| render_range(&cancelling_scene, cx, &s, &mut film, 1, 3))
        .expect("retry after cancellation");
    let direct = with_cx(|cx| render(&reference_scene, cx, &s)).expect("direct reference");
    assert_film_bits_eq(&film, &direct, "retry must equal a direct render bitwise");
}

#[test]
fn motion_path_zero_shutter_is_static_exact_and_missing_time_refuses() {
    let settings = motion_settings(256);
    let static_scene = emissive_motion_scene(false);
    let animated_scene = emissive_motion_scene(true);
    let static_film = with_cx(|cx| render(&static_scene, cx, &settings)).expect("static reference");
    let static_with_time =
        with_cx(|cx| render_motion(&static_scene, cx, &settings, motion_shutter(0.0, 1.0)))
            .expect("timed static render");
    assert_film_bits_eq(
        &static_with_time,
        &static_film,
        "drawing shutter time perturbed an existing static sample dimension",
    );
    let zero_shutter = motion_shutter(0.5, 0.0);
    let motion_film = with_cx(|cx| render_motion(&animated_scene, cx, &settings, zero_shutter))
        .expect("zero-width motion render");
    assert_film_bits_eq(
        &motion_film,
        &static_film,
        "zero-width shutter must equal the matching static pose",
    );

    let mut refused_film = Film::new(1, 1);
    let before = refused_film.clone();
    assert_eq!(
        with_cx(|cx| render_range(&animated_scene, cx, &settings, &mut refused_film, 0, 1)),
        Err(TracerError::MissingRayTime)
    );
    assert_film_state_bits_eq(
        &refused_film,
        &before,
        "missing ray time published a partial film",
    );
}

#[test]
fn constant_velocity_motion_blur_matches_analytic_occupancy_and_replays_progressively() {
    let settings = motion_settings(4_096);
    let static_scene = emissive_motion_scene(false);
    let animated_scene = emissive_motion_scene(true);
    let full_shutter = motion_shutter(0.0, 1.0);
    let static_film = with_cx(|cx| render(&static_scene, cx, &settings)).expect("static reference");
    let blurred = with_cx(|cx| render_motion(&animated_scene, cx, &settings, full_shutter))
        .expect("high-rate motion render");

    // The unit-width quad center moves linearly from x=-1 to x=1. The
    // center camera ray is covered exactly for t in [0.25, 0.75], so the
    // analytic temporal occupancy is one half. Spectral wavelength samples
    // are independently keyed, hence finite-sample channel ratios approach
    // (rather than being forced to) that physical occupancy.
    for channel in 0..3 {
        let ratio = blurred.xyz[0][channel] / static_film.xyz[0][channel];
        assert!(
            (ratio - 0.5).abs() < 0.025,
            "channel {channel} blur ratio {ratio} missed analytic occupancy 0.5"
        );
    }

    let frozen_start =
        with_cx(|cx| render_motion(&animated_scene, cx, &settings, motion_shutter(0.0, 0.0)))
            .expect("frozen start render");
    assert!(
        frozen_start.xyz[0].iter().all(|value| *value == 0.0),
        "the off-axis start pose unexpectedly covered the center ray"
    );
    assert!(
        blurred.xyz[0].iter().all(|value| *value > 0.0),
        "temporal integration collapsed to a frozen endpoint"
    );

    let mut progressive = Film::new(1, 1);
    with_cx(|cx| {
        render_motion_range(
            &animated_scene,
            cx,
            &settings,
            &mut progressive,
            0,
            997,
            full_shutter,
        )?;
        render_motion_range(
            &animated_scene,
            cx,
            &settings,
            &mut progressive,
            997,
            4_096,
            full_shutter,
        )
    })
    .expect("progressive motion render");
    assert_film_bits_eq(
        &progressive,
        &blurred,
        "motion sample partitions changed replay bits",
    );
}

#[test]
fn constant_spin_motion_blur_matches_analytic_angular_envelope() {
    let settings = motion_settings(4_096);
    let static_scene = emissive_spin_scene(false);
    let animated_scene = emissive_spin_scene(true);
    let static_film = with_cx(|cx| render(&static_scene, cx, &settings)).expect("static spinner");
    let blurred =
        with_cx(|cx| render_motion(&animated_scene, cx, &settings, motion_shutter(0.0, 1.0)))
            .expect("spinning motion render");

    // At the fixed camera ray (x=0.75,y=0), the thin local rectangle remains
    // covered while |0.75 sin(theta)| <= 0.1. The admitted SLERP is constant
    // angular speed from 0 to pi/2, so its exact occupancy fraction is
    // asin(0.1/0.75)/(pi/2).
    let analytic_occupancy = (0.1_f64 / 0.75).asin() / core::f64::consts::FRAC_PI_2;
    for channel in 0..3 {
        let ratio = blurred.xyz[0][channel] / static_film.xyz[0][channel];
        assert!(
            (ratio - analytic_occupancy).abs() < 0.02,
            "channel {channel} spin ratio {ratio} missed analytic occupancy {analytic_occupancy}"
        );
    }
}

#[test]
fn motion_shutter_outside_trajectory_refuses_before_film_publication() {
    let scene = emissive_motion_scene(true);
    let settings = motion_settings(4);
    let outside = ShutterInterval::resolve(
        1.0,
        0.5,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::UniformCounterV1,
        ShotTimeBounds::try_new(0.0, 2.0).expect("extended shot"),
    )
    .expect("otherwise valid shutter");
    let mut film = Film::new(1, 1);
    let before = film.clone();
    assert_eq!(
        with_cx(|cx| render_motion_range(&scene, cx, &settings, &mut film, 0, 4, outside)),
        Err(TracerError::MotionOutsideTrajectory)
    );
    assert_film_state_bits_eq(&film, &before, "out-of-trajectory shutter changed film");
}

#[test]
fn progressive_motion_checkpoint_rejects_shutter_changes_transactionally() {
    let scene = emissive_motion_scene(true);
    let settings = motion_settings(2);
    let mut film = Film::new(1, 1);
    let full = motion_shutter(0.0, 1.0);
    with_cx(|cx| render_motion_range(&scene, cx, &settings, &mut film, 0, 1, full))
        .expect("first motion partition");
    assert_eq!(
        film.time_mode,
        FilmTimeMode::Motion {
            shutter: full,
            stream_identity: settings.seed,
        }
    );
    let before = film.clone();
    assert_eq!(
        with_cx(|cx| {
            render_motion_range(
                &scene,
                cx,
                &settings,
                &mut film,
                1,
                2,
                motion_shutter(0.5, 0.0),
            )
        }),
        Err(TracerError::ProgressiveTimeModeMismatch)
    );
    assert_eq!(film, before);
}

#[test]
fn progressive_motion_checkpoint_rejects_shutter_stream_changes_transactionally() {
    let scene = emissive_motion_scene(true);
    let settings = motion_settings(2);
    let shutter = motion_shutter(0.0, 1.0);
    let mut film = Film::new(1, 1);
    with_cx(|cx| render_motion_range(&scene, cx, &settings, &mut film, 0, 1, shutter))
        .expect("first seeded motion partition");
    let before = film.clone();
    let mut reseeded = settings;
    reseeded.seed ^= 0x5eed;

    assert_eq!(
        with_cx(|cx| { render_motion_range(&scene, cx, &reseeded, &mut film, 1, 2, shutter) }),
        Err(TracerError::ProgressiveTimeModeMismatch)
    );
    assert_eq!(film, before);
}

#[test]
fn animated_nee_light_refuses_until_time_dependent_light_sampling_exists() {
    let mut scene = emissive_motion_scene(true);
    scene.lights[0].prim = 0;
    let settings = motion_settings(1);
    assert_eq!(
        with_cx(|cx| render_motion(&scene, cx, &settings, motion_shutter(0.5, 0.0))),
        Err(TracerError::AnimatedLightUnsupported)
    );
}

#[test]
fn cinematic_pinhole_matches_legacy_bits_for_iid_and_owen_sobol() {
    for sampler in [Sampler::Iid, Sampler::OwenSobol] {
        let scene = cornell();
        let camera = static_cinematic_camera(physical_from_legacy(
            &scene.camera,
            1.0,
            Aperture::try_circular(0.0).expect("pinhole aperture"),
        ));
        let settings = settings(DirectStrategy::Mis, sampler, 0x0070_696e_686f_6c65, 5, 3);
        let shutter = motion_shutter(0.5, 0.0);
        let (legacy, cinematic) = with_cx(|cx| {
            (
                render(&scene, cx, &settings).expect("legacy pinhole render"),
                render_cinematic(&scene, &camera, CutSide::After, cx, &settings, shutter)
                    .expect("cinematic pinhole render"),
            )
        });
        assert_film_bits_eq(
            &cinematic,
            &legacy,
            "zero-aperture cinematic path changed legacy XYZ bits",
        );
    }
}

#[test]
fn cinematic_primary_record_is_the_exact_beauty_sample_hit() {
    let scene = emissive_motion_scene(true);
    let camera = static_cinematic_camera(physical_from_legacy(
        &scene.camera,
        2.0,
        Aperture::try_circular(0.0).expect("pinhole aperture"),
    ));
    let settings = motion_settings(1);
    let shutter = motion_shutter(0.375, 0.0);
    let (film, sample) = with_cx(|cx| {
        let film = render_cinematic(&scene, &camera, CutSide::After, cx, &settings, shutter)?;
        let sample = trace_cinematic_pixel_sample(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            shutter,
            0,
            0,
        )?;
        Ok::<_, TracerError>((film, sample))
    })
    .expect("aligned cinematic sample");

    assert_eq!(sample.xyz.map(f64::to_bits), film.xyz[0].map(f64::to_bits));
    assert_eq!(sample.absolute_time_s.to_bits(), 0.375_f64.to_bits());
    let primary = sample
        .primary
        .expect("center ray accepts the instanced quad");
    assert_eq!(primary.primitive_index, 0);
    assert_eq!(
        primary.material_identity,
        scene.primitives[0].material.content_identity()
    );
    let surface = primary
        .surface
        .expect("instance retains local correspondence");
    assert_eq!(surface.identity().object_id(), 101);
    assert_eq!(
        surface.identity().material_identity(),
        primary.material_identity
    );
    assert!(matches!(
        surface.identity().feature(),
        StableFeatureIdentity::MeshTriangle(_)
    ));
    assert_eq!(primary.hit.point, Point3::new(0.0, 0.0, 0.0));
    assert_eq!(
        with_cx(|cx| trace_cinematic_pixel_sample(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            shutter,
            1,
            0,
        )),
        Err(TracerError::InvalidInput)
    );
}

#[test]
fn cinematic_depth_of_field_changes_depth_fixture_and_replays_progressively() {
    let half_extent = 0.05;
    let aperture_radius = 0.3;
    let scene = focus_probe_scene(1.0, half_extent);
    let pinhole = static_cinematic_camera(physical_from_legacy(
        &scene.camera,
        2.0,
        Aperture::try_circular(0.0).expect("pinhole aperture"),
    ));
    let finite_aperture = static_cinematic_camera(physical_from_legacy(
        &scene.camera,
        2.0,
        Aperture::try_circular(aperture_radius).expect("finite circular aperture"),
    ));
    let settings = Settings {
        width: 1,
        height: 1,
        spp: 4_096,
        max_depth: 1,
        sampler: Sampler::Iid,
        strategy: DirectStrategy::Mis,
        seed: 0x6465_7074_682d_6f66,
    };
    let shutter = motion_shutter(0.5, 0.0);
    let (pinhole_film, finite_aperture_film) = with_cx(|cx| {
        (
            render_cinematic(&scene, &pinhole, CutSide::After, cx, &settings, shutter)
                .expect("pinhole depth-fixture render"),
            render_cinematic(
                &scene,
                &finite_aperture,
                CutSide::After,
                cx,
                &settings,
                shutter,
            )
            .expect("finite-aperture depth-fixture render"),
        )
    });
    assert!(
        pinhole_film
            .xyz
            .iter()
            .flatten()
            .any(|channel| *channel > 0.0),
        "depth fixture produced no pinhole radiance"
    );
    // At z=1 the ray is halfway between a lens point and the z=0 focus
    // point. The target therefore admits lens coordinates in a centered
    // square of half-extent 2h. That square lies wholly inside the radius-r
    // disk, so its exact area occupancy is (4h)^2/(pi*r^2).
    let expected_occupancy = 16.0 * half_extent * half_extent
        / (core::f64::consts::PI * aperture_radius * aperture_radius);
    for channel in 0..3 {
        let ratio = finite_aperture_film.xyz[0][channel] / pinhole_film.xyz[0][channel];
        assert!(
            (ratio - expected_occupancy).abs() < 0.03,
            "channel {channel} defocus occupancy {ratio} missed analytic {expected_occupancy}"
        );
    }

    let focused_scene = focus_probe_scene(0.0, half_extent);
    let focused_pinhole = static_cinematic_camera(physical_from_legacy(
        &focused_scene.camera,
        2.0,
        Aperture::try_circular(0.0).expect("focused pinhole aperture"),
    ));
    let focused_aperture = static_cinematic_camera(physical_from_legacy(
        &focused_scene.camera,
        2.0,
        Aperture::try_circular(aperture_radius).expect("focused finite aperture"),
    ));
    let (focused_pinhole_film, focused_aperture_film) = with_cx(|cx| {
        (
            render_cinematic(
                &focused_scene,
                &focused_pinhole,
                CutSide::After,
                cx,
                &settings,
                shutter,
            )
            .expect("focused pinhole probe"),
            render_cinematic(
                &focused_scene,
                &focused_aperture,
                CutSide::After,
                cx,
                &settings,
                shutter,
            )
            .expect("focused finite-aperture probe"),
        )
    });
    assert_film_bits_eq(
        &focused_aperture_film,
        &focused_pinhole_film,
        "an on-focus emitter changed under finite aperture",
    );

    let mut progressive = Film::new(settings.width, settings.height);
    with_cx(|cx| {
        render_cinematic_range(
            &scene,
            &finite_aperture,
            CutSide::After,
            cx,
            &settings,
            &mut progressive,
            0,
            997,
            shutter,
        )?;
        render_cinematic_range(
            &scene,
            &finite_aperture,
            CutSide::After,
            cx,
            &settings,
            &mut progressive,
            997,
            settings.spp,
            shutter,
        )
    })
    .expect("progressive finite-aperture render");
    assert_film_state_bits_eq(
        &progressive,
        &finite_aperture_film,
        "lens samples changed across progressive partitions",
    );
}

#[test]
fn moving_keyframed_camera_replays_progressively_bitwise() {
    let scene = depth_varying_emissive_scene();
    let first_legacy = Camera {
        eye: Point3::new(-0.25, 0.0, 2.0),
        forward: scene.camera.forward,
        up: scene.camera.up,
        half_tan: scene.camera.half_tan,
    };
    let last_legacy = Camera {
        eye: Point3::new(0.25, 0.0, 2.0),
        forward: scene.camera.forward,
        up: scene.camera.up,
        half_tan: scene.camera.half_tan,
    };
    let shot = CameraShot::try_new(
        702,
        0.0,
        1.0,
        vec![
            CameraKeyframe::try_new(
                0.0,
                physical_from_legacy(
                    &first_legacy,
                    2.0,
                    Aperture::try_circular(0.0).expect("pinhole aperture"),
                ),
            )
            .expect("first camera keyframe"),
            CameraKeyframe::try_new(
                1.0,
                physical_from_legacy(
                    &last_legacy,
                    2.0,
                    Aperture::try_circular(0.0).expect("pinhole aperture"),
                ),
            )
            .expect("last camera keyframe"),
        ],
    )
    .expect("moving camera shot");
    let camera = AnimatedCamera::try_new(vec![shot]).expect("moving animated camera");
    let settings = Settings {
        width: 6,
        height: 6,
        spp: 11,
        max_depth: 1,
        sampler: Sampler::OwenSobol,
        strategy: DirectStrategy::Mis,
        seed: 0x6d6f_7669_6e67_6361,
    };
    let shutter = motion_shutter(0.0, 1.0);
    let full =
        with_cx(|cx| render_cinematic(&scene, &camera, CutSide::After, cx, &settings, shutter))
            .expect("one-shot moving-camera render");
    let mut progressive = Film::new(settings.width, settings.height);
    with_cx(|cx| {
        render_cinematic_range(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            &mut progressive,
            0,
            4,
            shutter,
        )?;
        render_cinematic_range(
            &scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            &mut progressive,
            4,
            settings.spp,
            shutter,
        )
    })
    .expect("progressive moving-camera render");
    assert_film_state_bits_eq(
        &progressive,
        &full,
        "moving-camera evaluation changed across progressive partitions",
    );
}

#[test]
fn moving_camera_and_geometry_share_one_absolute_path_time() {
    let start_s = 7.0;
    let end_s = 9.0;
    let animated_scene = absolute_time_translation_scene(true, start_s, end_s);
    let first = Camera {
        eye: Point3::new(-1.0, 0.0, 2.0),
        forward: animated_scene.camera.forward,
        up: animated_scene.camera.up,
        half_tan: animated_scene.camera.half_tan,
    };
    let last = Camera {
        eye: Point3::new(1.0, 0.0, 2.0),
        forward: animated_scene.camera.forward,
        up: animated_scene.camera.up,
        half_tan: animated_scene.camera.half_tan,
    };
    let camera = AnimatedCamera::try_new(vec![
        CameraShot::try_new(
            703,
            start_s,
            end_s,
            vec![
                CameraKeyframe::try_new(
                    start_s,
                    physical_from_legacy(
                        &first,
                        2.0,
                        Aperture::try_circular(0.0).expect("first pinhole"),
                    ),
                )
                .expect("first tracking keyframe"),
                CameraKeyframe::try_new(
                    end_s,
                    physical_from_legacy(
                        &last,
                        2.0,
                        Aperture::try_circular(0.0).expect("last pinhole"),
                    ),
                )
                .expect("last tracking keyframe"),
            ],
        )
        .expect("tracking camera shot"),
    ])
    .expect("tracking camera timeline");
    let settings = motion_settings(4_096);
    let shutter = ShutterInterval::resolve(
        start_s,
        end_s - start_s,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::StratifiedCounterV1 { strata: 4_096 },
        ShotTimeBounds::try_new(start_s, end_s).expect("absolute-time shot bounds"),
    )
    .expect("absolute-time shutter");
    let joint = with_cx(|cx| {
        render_cinematic(
            &animated_scene,
            &camera,
            CutSide::After,
            cx,
            &settings,
            shutter,
        )
    })
    .expect("joint camera-and-geometry motion render");

    // Object and camera both translate from x=-1 to x=+1 at the same path
    // time, so their relative pose is the static centered reference for every
    // sample. This is exact emissive occupancy, not a visual-difference proxy.
    let static_scene = absolute_time_translation_scene(false, start_s, end_s);
    let static_camera = AnimatedCamera::try_static(
        704,
        start_s,
        end_s,
        physical_from_legacy(
            &static_scene.camera,
            2.0,
            Aperture::try_circular(0.0).expect("static reference pinhole"),
        ),
    )
    .expect("static reference camera");
    let static_reference = with_cx(|cx| {
        render_cinematic(
            &static_scene,
            &static_camera,
            CutSide::After,
            cx,
            &settings,
            shutter,
        )
    })
    .expect("static relative-pose reference");
    assert_film_bits_eq(
        &joint,
        &static_reference,
        "camera and geometry did not evaluate at the same path time",
    );

    // The control keeps the camera fixed while the narrow object sweeps past.
    // A nonzero but incomplete occupancy proves the fixture would detect an
    // unshared, normalized-as-absolute, or frozen camera time.
    let fixed_camera_control = with_cx(|cx| {
        render_cinematic(
            &animated_scene,
            &static_camera,
            CutSide::After,
            cx,
            &settings,
            shutter,
        )
    })
    .expect("fixed-camera sensitivity control");
    for channel in 0..3 {
        assert!(
            fixed_camera_control.xyz[0][channel] > 0.0
                && fixed_camera_control.xyz[0][channel] < static_reference.xyz[0][channel],
            "channel {channel} fixed-camera control did not expose relative motion"
        );
    }
}

#[test]
fn cinematic_cut_and_out_of_shot_refusals_are_transactional() {
    let scene = depth_varying_emissive_scene();
    let physical = physical_from_legacy(
        &scene.camera,
        2.0,
        Aperture::try_circular(0.0).expect("pinhole aperture"),
    );
    let camera = AnimatedCamera::try_new(vec![
        CameraShot::try_new(
            801,
            0.0,
            0.5,
            vec![CameraKeyframe::try_new(0.0, physical.clone()).expect("first shot keyframe")],
        )
        .expect("first cut shot"),
        CameraShot::try_new(
            802,
            0.5,
            1.0,
            vec![CameraKeyframe::try_new(0.5, physical).expect("second shot keyframe")],
        )
        .expect("second cut shot"),
    ])
    .expect("hard-cut camera timeline");
    let settings = Settings {
        width: 2,
        height: 2,
        spp: 2,
        max_depth: 1,
        sampler: Sampler::Iid,
        strategy: DirectStrategy::Mis,
        seed: 0x6375_742d_7265_6675,
    };
    let crossing_cut = ShutterInterval::resolve(
        0.4,
        0.2,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::UniformCounterV1,
        ShotTimeBounds::try_new(0.0, 2.0).expect("admission test bounds"),
    )
    .expect("crossing-cut shutter syntax");
    let outside_timeline = ShutterInterval::resolve(
        1.25,
        0.0,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::UniformCounterV1,
        ShotTimeBounds::try_new(0.0, 2.0).expect("extended admission test bounds"),
    )
    .expect("out-of-shot shutter syntax");
    let mut film = Film::new(settings.width, settings.height);
    for xyz in &mut film.xyz {
        *xyz = [0.25, -0.0, f64::from_bits(0x7ff8_0000_0000_0701)];
    }
    let before = film.clone();

    for (shutter, context) in [
        (crossing_cut, "crossing-cut exposure"),
        (outside_timeline, "out-of-shot exposure"),
    ] {
        assert_eq!(
            with_cx(|cx| {
                render_cinematic_range(
                    &scene,
                    &camera,
                    CutSide::After,
                    cx,
                    &settings,
                    &mut film,
                    0,
                    settings.spp,
                    shutter,
                )
            }),
            Err(TracerError::Camera(CameraError::ShutterCrossesCut)),
            "{context} was not refused during admission"
        );
        assert_film_state_bits_eq(&film, &before, "{context} changed film state");
    }
}

#[test]
fn progressive_cinematic_checkpoint_binds_cut_side_and_camera_path() {
    let scene = depth_varying_emissive_scene();
    let physical = physical_from_legacy(
        &scene.camera,
        2.0,
        Aperture::try_circular(0.0).expect("pinhole aperture"),
    );
    let camera = AnimatedCamera::try_new(vec![
        CameraShot::try_new(
            811,
            0.0,
            0.5,
            vec![CameraKeyframe::try_new(0.0, physical.clone()).expect("outgoing keyframe")],
        )
        .expect("outgoing shot"),
        CameraShot::try_new(
            812,
            0.5,
            1.0,
            vec![CameraKeyframe::try_new(0.5, physical).expect("entering keyframe")],
        )
        .expect("entering shot"),
    ])
    .expect("hard-cut camera timeline");
    let settings = Settings {
        width: 2,
        height: 2,
        spp: 2,
        max_depth: 1,
        sampler: Sampler::Iid,
        strategy: DirectStrategy::Mis,
        seed: 0x6375_742d_6269_6e64,
    };
    let cut_instant = motion_shutter(0.5, 0.0);

    let mut cut_film = Film::new(settings.width, settings.height);
    with_cx(|cx| {
        render_cinematic_range(
            &scene,
            &camera,
            CutSide::Before,
            cx,
            &settings,
            &mut cut_film,
            0,
            1,
            cut_instant,
        )
    })
    .expect("outgoing cut-side partition");
    assert_eq!(
        cut_film.time_mode,
        FilmTimeMode::Cinematic {
            shutter: cut_instant,
            stream_identity: settings.seed,
            shot_id: 811,
        }
    );
    let before_side_switch = cut_film.clone();
    assert_eq!(
        with_cx(|cx| {
            render_cinematic_range(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings,
                &mut cut_film,
                1,
                2,
                cut_instant,
            )
        }),
        Err(TracerError::ProgressiveTimeModeMismatch)
    );
    assert_film_state_bits_eq(
        &cut_film,
        &before_side_switch,
        "cut-side mismatch changed film state",
    );

    let mut path_film = Film::new(settings.width, settings.height);
    with_cx(|cx| render_motion_range(&scene, cx, &settings, &mut path_film, 0, 1, cut_instant))
        .expect("legacy-motion partition");
    let before_path_switch = path_film.clone();
    assert_eq!(
        with_cx(|cx| {
            render_cinematic_range(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings,
                &mut path_film,
                1,
                2,
                cut_instant,
            )
        }),
        Err(TracerError::ProgressiveTimeModeMismatch)
    );
    assert_film_state_bits_eq(
        &path_film,
        &before_path_switch,
        "camera-path mismatch changed film state",
    );
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
    let v_power_mis = mean_pixel_variance(&scene, DirectStrategy::PowerMis, Sampler::Iid, 4, 12);
    let v_nee = mean_pixel_variance(&scene, DirectStrategy::NeeOnly, Sampler::Iid, 4, 12);
    let v_bsdf = mean_pixel_variance(&scene, DirectStrategy::BsdfOnly, Sampler::Iid, 4, 12);
    println!(
        "{{\"suite\":\"fs-render/tracer\",\"case\":\"mis-variance\",\"verdict\":\"info\",\"detail\":\"var balance_mis {v_mis:.3e} power_mis {v_power_mis:.3e} nee {v_nee:.3e} bsdf {v_bsdf:.3e}\"}}"
    );
    assert!(
        v_mis < v_nee && v_mis < v_bsdf,
        "MIS variance {v_mis:.3e} does not beat NEE {v_nee:.3e} / BSDF {v_bsdf:.3e}"
    );
    assert!(
        v_power_mis < v_nee && v_power_mis < v_bsdf,
        "power-MIS variance {v_power_mis:.3e} does not beat NEE {v_nee:.3e} / BSDF {v_bsdf:.3e}"
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

/// Controlled equal-SPP comparison between the established pixel-only Owen
/// stream and the independently scrambled per-bounce path blocks. This is a
/// measured fixture result, not a universal variance ordering claim.
#[test]
fn full_path_owen_vs_pixel_owen_equal_spp_logged() {
    let scene = cornell();
    let v_pixel = mean_pixel_variance(&scene, DirectStrategy::Mis, Sampler::OwenSobol, 16, 12);
    let v_full_path = mean_pixel_variance(
        &scene,
        DirectStrategy::Mis,
        Sampler::OwenSobolFullPath,
        16,
        12,
    );
    assert!(v_pixel.is_finite() && v_pixel > 0.0);
    assert!(v_full_path.is_finite() && v_full_path > 0.0);
    println!(
        "{{\"suite\":\"fs-render/tracer\",\"case\":\"full-path-owen-vs-pixel-owen-16spp\",\"verdict\":\"info\",\"detail\":\"var pixel {v_pixel:.3e} full_path {v_full_path:.3e} ratio {:.3}\"}}",
        v_full_path / v_pixel
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
