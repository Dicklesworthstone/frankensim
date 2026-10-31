//! Focused G0/G3/E2E battery for spectral dielectric path transport.
#![cfg(feature = "tracer")]

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::charts::TriMesh;
use fs_render::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use fs_render::spectral::lift_rgb;
use fs_render::tracer::{
    Camera, DirectStrategy, Film, Material, Primitive, RectLight, Sampler, Scene, Settings, Shape,
    TracerError, render, render_range,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    with_gate(&gate, f)
}

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x0000_0047_4c41_5353,
                kernel_id: 42,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn material(index: f64, extinction_per_m: f64, surface: DielectricSurface) -> Material {
    Material::Dielectric {
        glass: DielectricGlass::new(
            CauchyIor::try_constant(index).expect("fixture IOR"),
            BeerLambertAbsorption::try_constant(extinction_per_m).expect("fixture extinction"),
            GlassProvenance::Custom,
        ),
        surface,
    }
}

fn box_mesh_at(center_x: f64, z_front: f64, z_back: f64, reverse_winding: bool) -> TriMesh {
    let half = 4.0;
    let vertices = vec![
        [center_x - half, -half, z_back],
        [center_x + half, -half, z_back],
        [center_x + half, half, z_back],
        [center_x - half, half, z_back],
        [center_x - half, -half, z_front],
        [center_x + half, -half, z_front],
        [center_x + half, half, z_front],
        [center_x - half, half, z_front],
    ];
    let mut triangles = vec![
        [0, 2, 1],
        [0, 3, 2],
        [4, 5, 6],
        [4, 6, 7],
        [0, 1, 5],
        [0, 5, 4],
        [1, 2, 6],
        [1, 6, 5],
        [2, 3, 7],
        [2, 7, 6],
        [3, 0, 4],
        [3, 4, 7],
    ];
    if reverse_winding {
        for triangle in &mut triangles {
            triangle.swap(1, 2);
        }
    }
    TriMesh::new(vertices, triangles)
}

fn emitter_mesh_at(center_x: f64, z: f64, half: f64) -> TriMesh {
    TriMesh::new(
        vec![
            [center_x - half, -half, z],
            [center_x + half, -half, z],
            [center_x + half, half, z],
            [center_x - half, half, z],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
}

fn scene_with_boundaries(
    boundaries: &[(f64, f64, Material, bool)],
    emitter_z: f64,
    emitter_half: f64,
) -> Scene {
    scene_with_boundaries_at(0.0, boundaries, emitter_z, emitter_half)
}

fn scene_with_boundaries_at(
    center_x: f64,
    boundaries: &[(f64, f64, Material, bool)],
    emitter_z: f64,
    emitter_half: f64,
) -> Scene {
    let white = lift_rgb([1.0, 1.0, 1.0]);
    let mut primitives = boundaries
        .iter()
        .map(|&(front, back, material, reverse)| Primitive {
            shape: Shape::Mesh(box_mesh_at(center_x, front, back, reverse)),
            material,
            emission: None,
        })
        .collect::<Vec<_>>();
    let emission = (white, 3.0);
    let light_primitive = primitives.len();
    primitives.push(Primitive {
        shape: Shape::Mesh(emitter_mesh_at(center_x, emitter_z, emitter_half)),
        material: Material::Lambertian { reflectance: white },
        emission: Some(emission),
    });
    Scene {
        primitives,
        lights: vec![RectLight {
            corner: Point3::new(center_x - emitter_half, -emitter_half, emitter_z),
            edge_u: Vec3::new(2.0 * emitter_half, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 2.0 * emitter_half, 0.0),
            prim: light_primitive,
            emission,
        }],
        environment: None,
        camera: Camera {
            eye: Point3::new(center_x, 0.0, 1.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.0,
        },
    }
}

fn settings(strategy: DirectStrategy, spp: u32, max_depth: u32) -> Settings {
    Settings {
        width: 1,
        height: 1,
        spp,
        max_depth,
        sampler: Sampler::Iid,
        strategy,
        seed: 0x4449_454c_4543_5452,
    }
}

fn assert_film_bits_eq(left: &Film, right: &Film, context: &str) {
    assert_eq!(left.spp_done, right.spp_done, "{context}");
    assert_eq!(left.time_mode, right.time_mode, "{context}");
    assert_eq!(left.xyz.len(), right.xyz.len(), "{context}");
    for (left_pixel, right_pixel) in left.xyz.iter().zip(&right.xyz) {
        for channel in 0..3 {
            assert_eq!(
                left_pixel[channel].to_bits(),
                right_pixel[channel].to_bits(),
                "{context}: channel {channel}"
            );
        }
    }
}

fn assert_channel_ratio(
    numerator: &Film,
    denominator: &Film,
    expected: f64,
    tolerance: f64,
    context: &str,
) {
    for channel in 0..3 {
        let ratio = numerator.xyz[0][channel] / denominator.xyz[0][channel];
        assert!(
            (ratio - expected).abs() <= tolerance,
            "{context}: channel {channel} ratio {ratio:.17e}, expected {expected:.17e}"
        );
    }
}

#[test]
fn equal_ior_lossless_slab_is_a_null_boundary_under_every_strategy() {
    let clear = material(1.0, 0.0, DielectricSurface::SMOOTH);
    let scene = scene_with_boundaries(&[(0.0, -0.2, clear, false)], -1.0, 2.0);
    let films = with_cx(|cx| {
        [
            DirectStrategy::Mis,
            DirectStrategy::BsdfOnly,
            DirectStrategy::NeeOnly,
        ]
        .map(|strategy| render(&scene, cx, &settings(strategy, 16, 3)).expect("null slab render"))
    });
    assert_film_bits_eq(&films[0], &films[1], "MIS versus BSDF-only delta path");
    assert_film_bits_eq(&films[0], &films[2], "MIS versus NEE-only delta path");

    let bare = scene_with_boundaries(&[], -1.0, 2.0);
    let bare_film = with_cx(|cx| {
        render(&bare, cx, &settings(DirectStrategy::Mis, 16, 3)).expect("bare emitter")
    });
    assert_film_bits_eq(&films[0], &bare_film, "equal-IOR slab versus ambient");
}

#[test]
fn translated_thin_slab_preserves_the_null_boundary_result() {
    let clear = material(1.0, 0.0, DielectricSurface::SMOOTH);
    let boundary = [(0.0, -1.0e-4, clear, false)];
    let local = scene_with_boundaries_at(0.0, &boundary, -1.0, 2.0);
    let translated = scene_with_boundaries_at(1.0e9, &boundary, -1.0, 2.0);
    let settings = settings(DirectStrategy::Mis, 8, 3);
    let (local_film, translated_film) = with_cx(|cx| {
        (
            render(&local, cx, &settings).expect("local thin slab"),
            render(&translated, cx, &settings).expect("translated thin slab"),
        )
    });

    assert_film_bits_eq(
        &translated_film,
        &local_film,
        "rigid translation changed a thin null-boundary slab",
    );
}

#[test]
fn thin_and_thick_slabs_follow_exact_beer_lambert_scaling() {
    let sigma = 2.0;
    let clear = material(1.0, 0.0, DielectricSurface::SMOOTH);
    let absorbing = material(1.0, sigma, DielectricSurface::SMOOTH);
    let (baseline, thin, thick) = with_cx(|cx| {
        let render_one = |scene: &Scene| {
            render(scene, cx, &settings(DirectStrategy::Mis, 32, 3)).expect("Beer slab")
        };
        (
            render_one(&scene_with_boundaries(
                &[(0.0, -0.5, clear, false)],
                -1.0,
                2.0,
            )),
            render_one(&scene_with_boundaries(
                &[(0.0, -0.5, absorbing, false)],
                -1.0,
                2.0,
            )),
            render_one(&scene_with_boundaries(
                &[(0.0, -1.0, absorbing, false)],
                -1.5,
                2.0,
            )),
        )
    });
    assert_channel_ratio(
        &thin,
        &baseline,
        fs_math::det::exp(-1.0),
        2.0e-14,
        "0.5 m slab",
    );
    assert_channel_ratio(
        &thick,
        &baseline,
        fs_math::det::exp(-2.0),
        2.0e-14,
        "1.0 m slab",
    );
}

#[test]
fn spectral_absorption_reaches_native_slab_transport_without_rgb_reconstruction() {
    let spectral_material = |knots: &[(f64, f64)]| Material::Dielectric {
        glass: DielectricGlass::new(
            CauchyIor::try_constant(1.0).expect("null-index fixture"),
            BeerLambertAbsorption::try_from_spectral_extinction(knots).expect("spectral source"),
            GlassProvenance::Custom,
        ),
        surface: DielectricSurface::SMOOTH,
    };
    let clear = material(1.0, 0.0, DielectricSurface::SMOOTH);
    let constant = material(1.0, 2.0, DielectricSurface::SMOOTH);
    let flat = spectral_material(&[(380.0, 2.0), (550.0, 2.0), (780.0, 2.0)]);
    let selective = spectral_material(&[(380.0, 0.0), (500.0, 0.0), (501.0, 20.0), (780.0, 20.0)]);
    let [clear_film, constant_film, flat_film, selective_film] = with_cx(|cx| {
        [clear, constant, flat, selective].map(|material| {
            render(
                &scene_with_boundaries(&[(0.0, -0.5, material, false)], -1.0, 2.0),
                cx,
                &settings(DirectStrategy::Mis, 128, 3),
            )
            .expect("spectral slab render")
        })
    });
    assert_channel_ratio(
        &flat_film,
        &constant_film,
        1.0,
        2.0e-14,
        "flat spectral equivalence",
    );
    assert_channel_ratio(
        &flat_film,
        &clear_film,
        (-1.0_f64).exp(),
        2.0e-14,
        "spectral path length",
    );
    let blue_ratio = selective_film.xyz[0][2] / clear_film.xyz[0][2];
    let green_ratio = selective_film.xyz[0][1] / clear_film.xyz[0][1];
    assert!(
        blue_ratio > 3.0 * green_ratio,
        "wavelength-selective absorption was flattened: Z={blue_ratio}, Y={green_ratio}"
    );
    let other_declaration = spectral_material(&[(380.0, 2.0), (780.0, 2.0)]);
    assert_ne!(
        flat.content_identity(),
        other_declaration.content_identity(),
        "material identity must retain the supplied spectral knots"
    );
    assert_ne!(flat.content_identity(), selective.content_identity());
}

#[test]
fn nested_equal_ior_media_apply_only_the_active_medium_per_segment() {
    let clear_outer = material(1.0, 0.0, DielectricSurface::SMOOTH);
    let outer = material(1.0, 2.0, DielectricSurface::SMOOTH);
    let inner = material(1.0, 5.0, DielectricSurface::SMOOTH);
    let (baseline, nested) = with_cx(|cx| {
        let settings = settings(DirectStrategy::Mis, 24, 5);
        (
            render(
                &scene_with_boundaries(
                    &[
                        (0.0, -0.6, clear_outer, false),
                        (-0.2, -0.4, clear_outer, false),
                    ],
                    -1.0,
                    2.0,
                ),
                cx,
                &settings,
            )
            .expect("clear nested baseline"),
            render(
                &scene_with_boundaries(
                    &[(0.0, -0.6, outer, false), (-0.2, -0.4, inner, false)],
                    -1.0,
                    2.0,
                ),
                cx,
                &settings,
            )
            .expect("absorbing nested media"),
        )
    });
    let expected = fs_math::det::exp(-(2.0 * 0.4 + 5.0 * 0.2));
    assert_channel_ratio(&nested, &baseline, expected, 3.0e-14, "nested media");
}

#[test]
fn reversed_and_non_lifo_boundaries_refuse_transactionally() {
    let glass = material(1.0, 0.0, DielectricSurface::SMOOTH);
    let cases = [
        scene_with_boundaries(&[(0.0, -0.5, glass, true)], -1.0, 2.0),
        scene_with_boundaries(
            &[(0.0, -0.6, glass, false), (-0.2, -0.8, glass, false)],
            -1.0,
            2.0,
        ),
    ];
    with_cx(|cx| {
        for scene in &cases {
            let settings = settings(DirectStrategy::Mis, 1, 6);
            let mut film = Film::new(1, 1);
            let before = film.clone();
            assert!(matches!(
                render_range(scene, cx, &settings, &mut film, 0, 1),
                Err(TracerError::MediumStackMismatch { .. })
            ));
            assert_eq!(film, before, "refusal partially published a film");
        }
    });
}

#[test]
fn g1_terminal_rough_refraction_and_reflection_keep_the_complete_direct_light() {
    // At the last vertex only the light-sampling technique is evaluated.
    // MIS must therefore reduce exactly to NEE, including for the legacy
    // one-rectangle lighting configuration. An unevaluated BSDF competitor
    // incorrectly darkens both sides of the interface.
    with_cx(|cx| {
        let surface = DielectricSurface::try_rough(0.2).unwrap();
        let glass = material(1.5, 0.7, surface);
        for transmission in [true, false] {
            let emitter_z = if transmission { -0.1 } else { 0.5 };
            let emitter_x = if transmission { 0.0 } else { 0.2 };
            let mut scene = scene_with_boundaries(&[(0.0, -0.2, glass, false)], emitter_z, 0.04);
            // The reflected light is off the primary camera ray, ensuring
            // this measures the rough interface instead of the emitter.
            scene.primitives[1].shape = Shape::Mesh(emitter_mesh_at(emitter_x, emitter_z, 0.04));
            scene.lights[0].corner.x += emitter_x;
            let reference = render(&scene, cx, &settings(DirectStrategy::NeeOnly, 64, 1)).unwrap();
            assert!(
                reference.xyz[0][1] > 0.0,
                "fixture must receive direct light"
            );
            for strategy in [DirectStrategy::Mis, DirectStrategy::PowerMis] {
                let film = render(&scene, cx, &settings(strategy, 64, 1)).unwrap();
                assert_film_bits_eq(&film, &reference, "terminal dielectric direct light");
            }
        }
    });
}

#[test]
fn g1_terminal_diffuse_receiver_after_glass_keeps_the_complete_direct_light() {
    with_cx(|cx| {
        let glass = material(1.5, 0.7, DielectricSurface::SMOOTH);
        let mut scene = scene_with_boundaries(&[(0.0, -0.2, glass, false)], -0.3, 0.04);
        // This emitter illuminates the receiver from inside the air gap. It
        // stays off the primary transmitted ray through the glass slab.
        scene.primitives[1].shape = Shape::Mesh(emitter_mesh_at(0.3, -0.3, 0.04));
        scene.lights[0].corner.x += 0.3;
        scene.primitives.push(Primitive {
            shape: Shape::Mesh(emitter_mesh_at(0.0, -0.5, 4.0)),
            material: Material::Lambertian {
                reflectance: lift_rgb([0.8; 3]),
            },
            emission: None,
        });
        let reference = render(&scene, cx, &settings(DirectStrategy::NeeOnly, 64, 3)).unwrap();
        assert!(
            reference.xyz[0][1] > 0.0,
            "glass must transmit light to the receiver"
        );
        for strategy in [DirectStrategy::Mis, DirectStrategy::PowerMis] {
            let film = render(&scene, cx, &settings(strategy, 64, 3)).unwrap();
            assert_film_bits_eq(&film, &reference, "terminal receiver after glass");
        }
    });
}

#[test]
fn rough_transmission_nee_is_finite_and_absorbed_inside_glass() {
    let surface = DielectricSurface::try_rough(0.2).expect("rough glass");
    let clear = material(1.5, 0.0, surface);
    let absorbed = material(1.5, 8.0, surface);
    let (clear_film, absorbed_film) = with_cx(|cx| {
        let settings = settings(DirectStrategy::Mis, 128, 1);
        (
            render(
                &scene_with_boundaries(&[(0.0, -0.2, clear, false)], -0.1, 0.04),
                cx,
                &settings,
            )
            .expect("clear rough interface"),
            render(
                &scene_with_boundaries(&[(0.0, -0.2, absorbed, false)], -0.1, 0.04),
                cx,
                &settings,
            )
            .expect("absorbing rough interface"),
        )
    });
    for film in [&clear_film, &absorbed_film] {
        assert!(film.xyz[0].iter().all(|value| value.is_finite()));
        assert!(film.xyz[0][1] > 0.0, "rough transmitted NEE was black");
    }
    assert!(
        absorbed_film.xyz[0][1] < clear_film.xyz[0][1],
        "target-medium Beer attenuation did not reduce rough NEE"
    );
}

#[test]
fn dielectric_progressive_replay_is_bitwise_and_cancellation_is_transactional() {
    let glass = material(1.0, 0.7, DielectricSurface::SMOOTH);
    let scene = scene_with_boundaries(&[(0.0, -0.3, glass, false)], -1.0, 2.0);
    let settings = settings(DirectStrategy::Mis, 12, 3);
    let (direct, resumed) = with_cx(|cx| {
        let direct = render(&scene, cx, &settings).expect("direct glass render");
        let mut resumed = Film::new(1, 1);
        render_range(&scene, cx, &settings, &mut resumed, 0, 5).expect("first partition");
        render_range(&scene, cx, &settings, &mut resumed, 5, 12).expect("second partition");
        (direct, resumed)
    });
    assert_film_bits_eq(&direct, &resumed, "glass progressive replay");

    let gate = CancelGate::new();
    gate.request();
    with_gate(&gate, |cx| {
        let mut film = Film::new(1, 1);
        let before = film.clone();
        assert_eq!(
            render_range(&scene, cx, &settings, &mut film, 0, 1),
            Err(TracerError::Cancelled)
        );
        assert_eq!(
            film, before,
            "cancelled glass render published partial data"
        );
    });
}
