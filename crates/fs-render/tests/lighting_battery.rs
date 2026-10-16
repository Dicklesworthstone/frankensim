//! Direct transport battery for deterministic multi-emitter and environment
//! lighting (bead `frankensim-h7xu5.4.4`). Failures print the seed, estimates,
//! and relative residual so a sampling defect can be replayed directly.

#![cfg(feature = "tracer")]

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::charts::TriMesh;
use fs_render::lighting::EnvironmentMap;
use fs_render::spectral::lift_rgb;
use fs_render::tracer::{
    Camera, DirectStrategy, Film, Material, Primitive, RectLight, Sampler, Scene, Settings, Shape,
    TracerError, render, render_range,
};

const SEED: u64 = 0x6c69_6768_7469_6e67;

fn with_cx<R>(cancelled: bool, operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: SEED,
                kernel_id: 0x4c49_4748,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn quad(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> TriMesh {
    TriMesh::new(vec![a, b, c, d], vec![[0, 1, 2], [0, 2, 3]])
}

fn settings(spp: u32, strategy: DirectStrategy) -> Settings {
    Settings {
        width: 1,
        height: 1,
        spp,
        max_depth: 1,
        sampler: Sampler::Iid,
        strategy,
        seed: SEED,
    }
}

fn camera(direction: Vec3, up: Vec3) -> Camera {
    Camera {
        eye: Point3::new(0.0, 0.5, 0.0),
        forward: direction,
        up,
        half_tan: 0.0,
    }
}

fn environment_scene(environment: EnvironmentMap) -> Scene {
    Scene {
        primitives: Vec::new(),
        lights: Vec::new(),
        environment: Some(environment),
        camera: camera(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)),
    }
}

fn studio_scene(include_left: bool, include_right: bool) -> Scene {
    let diffuse = lift_rgb([0.72, 0.72, 0.72]);
    let emission = (lift_rgb([1.0, 0.82, 0.62]), 12.0);
    let mut primitives = vec![Primitive {
        shape: Shape::Mesh(quad(
            [-0.5, 0.0, -0.5],
            [-0.5, 0.0, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.0, -0.5],
        )),
        material: Material::Lambertian {
            reflectance: diffuse,
        },
        emission: None,
    }];
    let mut lights = Vec::new();
    for (enabled, x0) in [(include_left, -1.1), (include_right, 0.7)] {
        if !enabled {
            continue;
        }
        let primitive = primitives.len();
        primitives.push(Primitive {
            shape: Shape::Mesh(quad(
                [x0, 1.0, -0.2],
                [x0 + 0.4, 1.0, -0.2],
                [x0 + 0.4, 1.0, 0.2],
                [x0, 1.0, 0.2],
            )),
            material: Material::Lambertian {
                reflectance: diffuse,
            },
            emission: Some(emission),
        });
        lights.push(RectLight {
            corner: Point3::new(x0, 1.0, -0.2),
            edge_u: Vec3::new(0.4, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 0.0, 0.4),
            prim: primitive,
            emission,
        });
    }
    Scene {
        primitives,
        lights,
        environment: None,
        camera: camera(Vec3::new(0.0, -1.0, 0.0), Vec3::new(0.0, 0.0, 1.0)),
    }
}

fn mean_y(film: &Film) -> f64 {
    film.xyz[0][1] / f64::from(film.spp_done)
}

fn assert_film_bits_eq(left: &Film, right: &Film, context: &str) {
    assert_eq!(left.spp_done, right.spp_done, "{context}: spp");
    assert_eq!(left.time_mode, right.time_mode, "{context}: time mode");
    assert_eq!(left.xyz.len(), right.xyz.len(), "{context}: film size");
    for (pixel, (left, right)) in left.xyz.iter().zip(&right.xyz).enumerate() {
        assert_eq!(
            left.map(f64::to_bits),
            right.map(f64::to_bits),
            "{context}: pixel={pixel}"
        );
    }
}

#[test]
fn g0_no_finite_emitter_refuses_before_sampling() {
    let scene = Scene {
        primitives: Vec::new(),
        lights: Vec::new(),
        environment: None,
        camera: camera(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)),
    };
    let result = with_cx(false, |cx| {
        render(&scene, cx, &settings(1, DirectStrategy::Mis))
    });
    assert!(
        matches!(result, Err(TracerError::Lighting(_))),
        "seed={SEED:#x}: an emitter-free scene must return a lighting refusal; got {result:?}"
    );
}

#[test]
fn g0_rectangle_metadata_must_match_the_named_emissive_geometry() {
    let mut scene = studio_scene(true, false);
    let primitive = scene.lights[0].prim;
    scene.lights[0].corner.x += 0.01;
    let result = with_cx(false, |cx| {
        render(&scene, cx, &settings(1, DirectStrategy::Mis))
    });
    assert_eq!(
        result,
        Err(TracerError::LightPrimitiveMismatch {
            light_primitive: primitive
        }),
        "seed={SEED:#x}: shifted rectangle metadata must not borrow the PDF/visibility identity of same-emission geometry"
    );
}

#[test]
fn g3_multi_light_order_is_bit_stable_and_matches_summed_single_light_estimates() {
    let spp = 16_384;
    let nee_settings = settings(spp, DirectStrategy::NeeOnly);
    let mut ordered = studio_scene(true, true);
    let ordered_film = with_cx(false, |cx| render(&ordered, cx, &nee_settings)).expect("ordered");
    let final_depth_mis = with_cx(false, |cx| {
        render(&ordered, cx, &settings(spp, DirectStrategy::Mis))
    })
    .expect("final-depth MIS");
    assert_film_bits_eq(
        &ordered_film,
        &final_depth_mis,
        "new lighting-v1 final-depth MIS has no untraced BSDF competitor",
    );
    ordered.lights.swap(0, 1);
    let reversed_film = with_cx(false, |cx| render(&ordered, cx, &nee_settings)).expect("reversed");
    assert_film_bits_eq(
        &ordered_film,
        &reversed_film,
        "construction-order-independent light identity",
    );

    let left = with_cx(false, |cx| {
        render(&studio_scene(true, false), cx, &nee_settings)
    })
    .expect("left reference");
    let right = with_cx(false, |cx| {
        render(&studio_scene(false, true), cx, &nee_settings)
    })
    .expect("right reference");
    let observed = mean_y(&ordered_film);
    let expected = mean_y(&left) + mean_y(&right);
    let relative = (observed - expected).abs() / expected.abs().max(1.0e-15);
    assert!(
        relative < 0.025,
        "seed={SEED:#x} spp={spp}: multi-light estimate={observed:.17e}, summed single-light reference={expected:.17e}, relative residual={relative:.6e}"
    );
}

#[test]
fn g3_constant_environment_is_rotation_invariant_and_hdr_transport_is_finite() {
    let pixels = vec![[10_000.0_f32, 250.0, 2.0]; 8];
    let unrotated =
        EnvironmentMap::try_from_linear_srgb(4, 2, pixels.clone(), 0.0).expect("constant HDR map");
    let rotated = EnvironmentMap::try_from_linear_srgb(4, 2, pixels, core::f64::consts::FRAC_PI_2)
        .expect("rotated constant HDR map");
    let settings = settings(32, DirectStrategy::Mis);
    let first = with_cx(false, |cx| {
        render(&environment_scene(unrotated), cx, &settings)
    })
    .expect("unrotated environment");
    let second = with_cx(false, |cx| {
        render(&environment_scene(rotated), cx, &settings)
    })
    .expect("rotated environment");
    assert_film_bits_eq(&first, &second, "constant environment rotation");
    assert!(
        first.xyz.iter().flatten().all(|value| value.is_finite()),
        "seed={SEED:#x}: high-dynamic-range environment produced NaN/Inf: {:?}",
        first.xyz
    );
    assert!(mean_y(&first) > 0.0, "HDR environment produced no radiance");
}

#[test]
fn g3_environment_rotation_changes_directional_radiance() {
    let pixels = vec![
        [8.0, 0.1, 0.1],
        [0.1, 8.0, 0.1],
        [0.1, 0.1, 8.0],
        [1.0, 1.0, 1.0],
        [8.0, 0.1, 0.1],
        [0.1, 8.0, 0.1],
        [0.1, 0.1, 8.0],
        [1.0, 1.0, 1.0],
    ];
    let first = EnvironmentMap::try_from_linear_srgb(4, 2, pixels.clone(), 0.0)
        .expect("directional environment");
    let second = EnvironmentMap::try_from_linear_srgb(4, 2, pixels, core::f64::consts::FRAC_PI_2)
        .expect("rotated directional environment");
    let settings = settings(64, DirectStrategy::Mis);
    let first = with_cx(false, |cx| render(&environment_scene(first), cx, &settings))
        .expect("first orientation");
    let second = with_cx(false, |cx| {
        render(&environment_scene(second), cx, &settings)
    })
    .expect("second orientation");
    assert_ne!(
        first.xyz[0].map(f64::to_bits),
        second.xyz[0].map(f64::to_bits),
        "seed={SEED:#x}: environment rotation was not consumed by miss radiance"
    );
}

#[test]
fn g4_mixed_lighting_progressive_replay_and_cancellation_are_transactional() {
    let environment = EnvironmentMap::try_from_linear_srgb(2, 1, vec![[0.25, 0.3, 0.4]; 2], 0.0)
        .expect("constant environment");
    let mut scene = studio_scene(true, true);
    scene.environment = Some(environment);
    let settings = settings(32, DirectStrategy::Mis);
    let straight = with_cx(false, |cx| render(&scene, cx, &settings)).expect("straight render");
    let split = with_cx(false, |cx| {
        let mut film = Film::new(1, 1);
        render_range(&scene, cx, &settings, &mut film, 0, 13).expect("first range");
        render_range(&scene, cx, &settings, &mut film, 13, 32).expect("second range");
        film
    });
    assert_film_bits_eq(&straight, &split, "mixed-light progressive replay");

    let mut film = Film::new(1, 1);
    let before = film.clone();
    let cancelled = with_cx(true, |cx| {
        render_range(&scene, cx, &settings, &mut film, 0, 32)
    });
    assert_eq!(cancelled, Err(TracerError::Cancelled));
    assert_film_bits_eq(&film, &before, "cancelled mixed-light transaction");
}
