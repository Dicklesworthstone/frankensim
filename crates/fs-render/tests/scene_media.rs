//! G1/G3/G4/G5: physical cell transport belongs to the active glass interior.
#![cfg(all(feature = "tracer", feature = "volumes"))]

#[path = "scene_media/exposure.rs"]
mod exposure;

use fs_exec::{Budget, CancelGate, Cx, ExecMode, RunId, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::charts::TriMesh;
use fs_render::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use fs_render::spectral::{LiftedSpectrum, lift_rgb, xyz_of_spectrum};
use fs_render::thermal_surface::GrayDiffuseThermalSurface;
use fs_render::tracer::{
    Camera, DirectStrategy, Film, InteriorMedium, Material, Primitive, RectLight,
    RenderExecutionConfig, Sampler, Scene, SceneMedia, SceneMediaBudget, SceneMediaFilm, Settings,
    Shape, TracerError, render, render_scene_media, render_scene_media_range,
    render_scene_media_range_with_execution, render_scene_media_with_execution,
};
use fs_render::volumes::VolumeGrid;
use fs_render::volumes::transport::scattering::{
    MassScatteringSpectrum, PhaseFunction, ScatteringGrid,
};
use fs_render::volumes::transport::timeline::PhysicalMedium;
use fs_render::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        operation(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 173,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}

fn media_budget() -> SceneMediaBudget {
    SceneMediaBudget {
        max_interior_bindings: 8,
        max_primitives: 16,
        max_light_triangles: 32,
    }
}

fn settings(strategy: DirectStrategy, spp: u32, max_depth: u32) -> Settings {
    Settings {
        width: 1,
        height: 1,
        spp,
        max_depth,
        sampler: Sampler::OwenSobolFullPath,
        strategy,
        seed: 173,
    }
}

fn glass(index: CauchyIor, extinction: f64, surface: DielectricSurface) -> Material {
    Material::Dielectric {
        glass: DielectricGlass::new(
            index,
            BeerLambertAbsorption::try_constant(extinction).unwrap(),
            GlassProvenance::Custom,
        ),
        surface,
    }
}

fn slab(front: f64, back: f64, half: f64, material: Material) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [-half, -half, back],
                [half, -half, back],
                [half, half, back],
                [-half, half, back],
                [-half, -half, front],
                [half, -half, front],
                [half, half, front],
                [-half, half, front],
            ],
            vec![
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
            ],
        )),
        material,
        emission: None,
    }
}

fn rectangle(corner: Point3, u: Vec3, v: Vec3) -> Primitive {
    let point = |a: f64, b: f64| {
        [
            corner.x + a * u.x + b * v.x,
            corner.y + a * u.y + b * v.y,
            corner.z + a * u.z + b * v.z,
        ]
    };
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                point(0.0, 0.0),
                point(1.0, 0.0),
                point(1.0, 1.0),
                point(0.0, 1.0),
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )),
        material: Material::Lambertian {
            reflectance: lift_rgb([0.0; 3]),
        },
        emission: None,
    }
}

fn empty_scene(primitives: Vec<Primitive>) -> Scene {
    Scene {
        primitives,
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 1.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.0,
        },
    }
}

fn add_light(scene: &mut Scene, corner: Point3, u: Vec3, v: Vec3) -> (LiftedSpectrum, f64) {
    let emission = (lift_rgb([1.0; 3]), 3.0);
    let prim = scene.primitives.len();
    let mut source = rectangle(corner, u, v);
    source.emission = Some(emission);
    scene.primitives.push(source);
    scene.lights.push(RectLight {
        corner,
        edge_u: u,
        edge_v: v,
        prim,
        emission,
    });
    emission
}

fn transmission_scene(
    boundaries: Vec<Primitive>,
    emitter_z: f64,
) -> (Scene, (LiftedSpectrum, f64)) {
    let mut scene = empty_scene(boundaries);
    let emission = add_light(
        &mut scene,
        Point3::new(-8.0, -8.0, emitter_z),
        Vec3::new(16.0, 0.0, 0.0),
        Vec3::new(0.0, 16.0, 0.0),
    );
    (scene, emission)
}

fn field<'a>(
    cx: &Cx<'_>,
    density: &'a [f64],
    origin: [f64; 3],
    extent: [f64; 3],
    absorption: MassAbsorptionSpectrum,
    temperature: Option<&'a [f64]>,
) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new(
            [1, 1, density.len()],
            density,
            origin,
            [extent[0], extent[1], extent[2] / density.len() as f64],
        ),
        absorption,
        temperature,
        AbsorbingGridBudget {
            max_grid_cells: 64,
            max_steps_per_segment: 64,
        },
    )
    .unwrap()
}

fn assert_mean(film: &Film, expected: [f64; 3], tolerance: f64) {
    for (channel, expected) in expected.into_iter().enumerate() {
        let actual = film.xyz[0][channel] / f64::from(film.spp_done);
        assert!(
            (actual / expected - 1.0).abs() <= tolerance,
            "channel {channel}: {actual:.12e}, expected {expected:.12e}, tolerance {tolerance}"
        );
    }
}

fn assert_ratio(actual: &Film, reference: &Film, expected: f64) {
    for channel in 0..3 {
        let ratio = actual.xyz[0][channel] / reference.xyz[0][channel];
        assert!(
            (ratio / expected - 1.0).abs() < 3.0e-13,
            "channel {channel}: ratio {ratio:.17e}, expected {expected:.17e}"
        );
    }
}

// Independent unpolarized Fresnel and Snell reference, without renderer APIs.
fn interface(index: f64, cosine: f64) -> (f64, f64) {
    let transmitted = (1.0 - (1.0 - cosine * cosine) / (index * index)).sqrt();
    let rs = (cosine - index * transmitted) / (cosine + index * transmitted);
    let rp = (index * cosine - transmitted) / (index * cosine + transmitted);
    (0.5 * (rs * rs + rp * rp), transmitted)
}

fn planck(temperature: f64, wavelength_nm: f64) -> f64 {
    let wavelength = wavelength_nm * 1.0e-9;
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0;
    let k = 1.380_649e-23;
    2.0 * h * c * c / (wavelength.powi(5) * (h * c / (wavelength * k * temperature)).exp_m1())
        * 1.0e-9
}

#[test]
fn g1_bound_absorption_and_vacuum_replace_homogeneous_glass_extinction() {
    with_cx(|cx, _| {
        let ior = CauchyIor::try_constant(1.5).unwrap();
        let make = |sigma| {
            transmission_scene(
                vec![slab(
                    0.0,
                    -0.4,
                    4.0,
                    glass(ior, sigma, DielectricSurface::SMOOTH),
                )],
                -1.0,
            )
            .0
        };
        let clear = make(0.0);
        let opaque_host = make(80.0);
        let settings = settings(DirectStrategy::BsdfOnly, 64, 3);
        let reference = render(&clear, cx, &settings).unwrap();
        let density = [1.0];
        for (coefficient, expected) in [(0.0, 1.0), (2.0, (-0.8_f64).exp())] {
            let grid = field(
                cx,
                &density,
                [-4.0, -4.0, -0.4],
                [8.0, 8.0, 0.4],
                MassAbsorptionSpectrum::constant(coefficient).unwrap(),
                None,
            );
            let bindings = [InteriorMedium {
                boundary_primitive: 0,
                medium: PhysicalMedium::Absorbing(&grid),
            }];
            let media =
                SceneMedia::try_new(cx, &opaque_host, None, &bindings, media_budget()).unwrap();
            assert_ratio(
                &render_scene_media(&media, cx, &settings).unwrap(),
                &reference,
                expected,
            );
        }
    });
}

#[test]
fn g1_refracted_heterogeneous_slab_and_ambient_match_independent_optical_depth() {
    with_cx(|cx, _| {
        let absorption =
            MassAbsorptionSpectrum::table(&[(380.0, 0.1), (580.0, 0.6), (780.0, 1.0)]).unwrap();
        let density = [0.25, 1.5, 3.0];
        let grid = field(
            cx,
            &density,
            [-4.0, -4.0, -0.6],
            [8.0, 8.0, 0.6],
            absorption,
            None,
        );
        let ambient_density = [0.8];
        let ambient = field(
            cx,
            &ambient_density,
            [-8.0, -8.0, -1.0],
            [16.0, 16.0, 2.0],
            MassAbsorptionSpectrum::constant(0.2).unwrap(),
            None,
        );
        let (mut scene, (source, scale)) = transmission_scene(
            vec![slab(
                0.0,
                -0.6,
                4.0,
                glass(
                    CauchyIor::try_new(1.4, 0.012, 0.0).unwrap(),
                    100.0,
                    DielectricSurface::SMOOTH,
                ),
            )],
            -1.0,
        );
        scene.camera.forward = Vec3::new(0.6, 0.0, -0.8);
        let bindings = [InteriorMedium {
            boundary_primitive: 0,
            medium: PhysicalMedium::Absorbing(&grid),
        }];
        let media = SceneMedia::try_new(
            cx,
            &scene,
            Some(PhysicalMedium::Absorbing(&ambient)),
            &bindings,
            media_budget(),
        )
        .unwrap();
        let expected = xyz_of_spectrum(|lambda| {
            let um = lambda * 1.0e-3;
            let (r, cos_t) = interface(1.4 + 0.012 / (um * um), 0.8);
            let interior_depth = (0.25 + 1.5 + 3.0) * 0.2 * absorption.at(lambda).unwrap() / cos_t;
            let ambient_depth = 0.8 * 0.2 * (2.0 - 0.6) / 0.8;
            source.eval(lambda)
                * scale
                * (1.0 - r).powi(2)
                * (-interior_depth - ambient_depth).exp()
        });
        assert_mean(
            &render_scene_media(&media, cx, &settings(DirectStrategy::Mis, 4096, 3)).unwrap(),
            expected,
            0.006,
        );
    });
}

#[test]
fn g1_nested_regions_restore_the_outer_field_after_inner_exit() {
    with_cx(|cx, _| {
        let density = [1.0];
        let outer = field(
            cx,
            &density,
            [-4.0, -4.0, -1.0],
            [8.0, 8.0, 1.0],
            MassAbsorptionSpectrum::constant(0.4).unwrap(),
            None,
        );
        let ambient = field(
            cx,
            &density,
            [-8.0, -8.0, -1.5],
            [16.0, 16.0, 2.5],
            MassAbsorptionSpectrum::constant(0.2).unwrap(),
            None,
        );
        let index = CauchyIor::try_constant(1.0).unwrap();
        let (scene, _) = transmission_scene(
            vec![
                slab(
                    0.0,
                    -1.0,
                    4.0,
                    glass(index, 100.0, DielectricSurface::SMOOTH),
                ),
                slab(
                    -0.3,
                    -0.7,
                    3.0,
                    glass(index, 200.0, DielectricSurface::SMOOTH),
                ),
            ],
            -1.5,
        );
        let bare = transmission_scene(Vec::new(), -1.5).0;
        let settings = settings(DirectStrategy::BsdfOnly, 64, 5);
        let reference = render(&bare, cx, &settings).unwrap();
        for inner_sigma in [0.0, 1.4] {
            let inner = field(
                cx,
                &density,
                [-3.0, -3.0, -0.7],
                [6.0, 6.0, 0.4],
                MassAbsorptionSpectrum::constant(inner_sigma).unwrap(),
                None,
            );
            let bindings = [
                InteriorMedium {
                    boundary_primitive: 0,
                    medium: PhysicalMedium::Absorbing(&outer),
                },
                InteriorMedium {
                    boundary_primitive: 1,
                    medium: PhysicalMedium::Absorbing(&inner),
                },
            ];
            let media = SceneMedia::try_new(
                cx,
                &scene,
                Some(PhysicalMedium::Absorbing(&ambient)),
                &bindings,
                media_budget(),
            )
            .unwrap();
            let expected = (-(0.4 * 0.6 + inner_sigma * 0.4 + 0.2 * 1.5)).exp();
            assert_ratio(
                &render_scene_media(&media, cx, &settings).unwrap(),
                &reference,
                expected,
            );
        }
    });
}

#[test]
fn g1_hot_dispersive_interior_obeys_external_planck_radiance() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [1800.0];
        let grid = field(
            cx,
            &density,
            [-4.0, -4.0, -0.6],
            [8.0, 8.0, 0.6],
            MassAbsorptionSpectrum::constant(30.0).unwrap(),
            Some(&temperature),
        );
        let scene = empty_scene(vec![slab(
            0.0,
            -0.6,
            4.0,
            glass(
                CauchyIor::try_new(1.4, 0.012, 0.0).unwrap(),
                90.0,
                DielectricSurface::SMOOTH,
            ),
        )]);
        let bindings = [InteriorMedium {
            boundary_primitive: 0,
            medium: PhysicalMedium::Absorbing(&grid),
        }];
        let media = SceneMedia::try_new(cx, &scene, None, &bindings, media_budget()).unwrap();
        // This depth sees the front transmission and complete first interior
        // source integral, but no internally reflected return segment. The
        // host n² Planck source cancels the camera's radiance eta factor.
        let expected = xyz_of_spectrum(|lambda| {
            let um = lambda * 1.0e-3;
            let (r, _) = interface(1.4 + 0.012 / (um * um), 1.0);
            planck(1800.0, lambda) * (1.0 - r) * -(-18.0_f64).exp_m1()
        });
        assert_mean(
            &render_scene_media(&media, cx, &settings(DirectStrategy::BsdfOnly, 4096, 2)).unwrap(),
            expected,
            0.006,
        );
    });
}

fn scattering_scene() -> (Scene, (LiftedSpectrum, f64)) {
    let mut scene = empty_scene(vec![
        slab(
            0.0,
            -1.0,
            4.0,
            glass(
                CauchyIor::try_constant(1.5).unwrap(),
                40.0,
                DielectricSurface::SMOOTH,
            ),
        ),
        rectangle(
            Point3::new(-4.0, -4.0, -0.8),
            Vec3::new(8.0, 0.0, 0.0),
            Vec3::new(0.0, 8.0, 0.0),
        ),
    ]);
    let emission = add_light(
        &mut scene,
        Point3::new(1.0, -0.5, -0.75),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 0.5),
    );
    (scene, emission)
}

// Independent single-scatter volume/area integral. t runs from the front
// boundary to a black stop. The source plane's projected solid angle is
// dA/r³, and both flight and shadow attenuation use sigma_a + sigma_s.
fn single_scatter_integral(bins: u32) -> f64 {
    let step_t = 0.8 / f64::from(bins);
    let step_y = 1.0 / f64::from(bins);
    let step_z = 0.5 / f64::from(bins);
    let mut integral = 0.0;
    for it in 0..bins {
        let t = (f64::from(it) + 0.5) * step_t;
        for iy in 0..bins {
            let y = -0.5 + (f64::from(iy) + 0.5) * step_y;
            for iz in 0..bins {
                let z = -0.75 + (f64::from(iz) + 0.5) * step_z;
                let distance = (1.0 + y * y + (z + t).powi(2)).sqrt();
                integral += (-0.8 * (t + distance)).exp() / distance.powi(3);
            }
        }
    }
    integral * step_t * step_y * step_z * 0.6 / (4.0 * core::f64::consts::PI)
}

#[test]
fn g1_scattering_inside_refractive_body_matches_independent_side_light_integral() {
    with_cx(|cx, _| {
        let density = [1.0];
        let absorption = field(
            cx,
            &density,
            [-4.0, -4.0, -1.0],
            [8.0, 8.0, 1.0],
            MassAbsorptionSpectrum::constant(0.2).unwrap(),
            None,
        );
        let scattering = ScatteringGrid::try_new(
            cx,
            absorption,
            MassScatteringSpectrum::constant(0.6).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let (scene, (source, scale)) = scattering_scene();
        let bindings = [InteriorMedium {
            boundary_primitive: 0,
            medium: PhysicalMedium::Scattering(&scattering),
        }];
        let media = SceneMedia::try_new(cx, &scene, None, &bindings, media_budget()).unwrap();
        let integral = single_scatter_integral(48);
        assert!((single_scatter_integral(24) / integral - 1.0).abs() < 0.001);
        let (r, _) = interface(1.5, 1.0);
        let expected = xyz_of_spectrum(|lambda| {
            source.eval(lambda) * scale * (1.0 - r) / 1.5_f64.powi(2) * integral
        });
        assert_mean(
            &render_scene_media(&media, cx, &settings(DirectStrategy::NeeOnly, 16_384, 2)).unwrap(),
            expected,
            0.035,
        );
    });
}

#[test]
fn g0_phase_connections_refuse_only_visible_thermal_surfaces_in_glass() {
    with_cx(|cx, _| {
        let density = [1.0];
        let absorption = field(
            cx,
            &density,
            [-4.0, -4.0, -1.0],
            [8.0, 8.0, 1.0],
            MassAbsorptionSpectrum::constant(0.2).unwrap(),
            None,
        );
        let scattering = ScatteringGrid::try_new(
            cx,
            absorption,
            MassScatteringSpectrum::constant(0.6).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let thermal = Material::GrayDiffuseThermal {
            surface: GrayDiffuseThermalSurface::try_new(
                1800.0,
                1.0,
                fs_blake3::ContentHash([7; 32]),
            )
            .unwrap(),
        };
        for inside in [true, false] {
            let (mut scene, _) = scattering_scene();
            scene.lights.clear();
            let corner_z = if inside { -0.75 } else { 0.2 };
            let mut source = rectangle(
                Point3::new(1.0, -0.5, corner_z),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(0.0, 0.0, 0.5),
            );
            source.material = thermal;
            scene.primitives[2] = source;
            let bindings = [InteriorMedium {
                boundary_primitive: 0,
                medium: PhysicalMedium::Scattering(&scattering),
            }];
            let media = SceneMedia::try_new(cx, &scene, None, &bindings, media_budget()).unwrap();
            let result = render_scene_media(&media, cx, &settings(DirectStrategy::NeeOnly, 64, 2));
            if inside {
                assert!(matches!(
                    result,
                    Err(TracerError::ThermalSurfaceMediumUnsupported)
                ));
            } else {
                assert_eq!(result.unwrap().xyz, vec![[0.0; 3]]);
            }
        }
    });
}

#[test]
fn g3_rough_boundary_sources_on_both_sides_keep_complementary_mis() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [2200.0];
        let ambient = field(
            cx,
            &density,
            [-1.0, -1.0, 1.0],
            [2.0, 2.0, 1.0],
            MassAbsorptionSpectrum::constant(0.7).unwrap(),
            Some(&temperature),
        );
        let interior = field(
            cx,
            &density,
            // Compact lateral support gives both source and refracted BSDF
            // proposals useful probability near the narrow transmission lobe.
            [-0.05, -0.05, -0.2],
            [0.1, 0.1, 0.2],
            MassAbsorptionSpectrum::constant(0.7).unwrap(),
            Some(&temperature),
        );
        let scene = {
            let mut scene = empty_scene(vec![
                slab(
                    0.0,
                    -1.0,
                    8.0,
                    glass(
                        CauchyIor::try_constant(1.5).unwrap(),
                        0.0,
                        DielectricSurface::try_rough(0.25).unwrap(),
                    ),
                ),
                rectangle(
                    Point3::new(-8.0, -8.0, -0.2),
                    Vec3::new(16.0, 0.0, 0.0),
                    Vec3::new(0.0, 16.0, 0.0),
                ),
            ]);
            scene.camera.eye.z = 0.5;
            scene
        };
        let bindings = [InteriorMedium {
            boundary_primitive: 0,
            medium: PhysicalMedium::Absorbing(&interior),
        }];
        for with_ambient in [false, true] {
            let media = SceneMedia::try_new(
                cx,
                &scene,
                with_ambient.then_some(PhysicalMedium::Absorbing(&ambient)),
                &bindings,
                media_budget(),
            )
            .unwrap();
            let mut policy = settings(DirectStrategy::BsdfOnly, 32_768, 2);
            let reference = render_scene_media(&media, cx, &policy).unwrap();
            let expected = reference.xyz[0].map(|value| value / f64::from(policy.spp));
            assert!(expected.iter().all(|value| *value > 0.0));
            for strategy in [
                DirectStrategy::NeeOnly,
                DirectStrategy::Mis,
                DirectStrategy::PowerMis,
            ] {
                policy.strategy = strategy;
                let film = render_scene_media(&media, cx, &policy).unwrap();
                for (channel, expected) in expected.into_iter().enumerate() {
                    let actual = film.xyz[0][channel] / f64::from(film.spp_done);
                    assert!(
                        (actual / expected - 1.0).abs() < 0.07,
                        "ambient={with_ambient}, {strategy:?}, channel {channel}: {actual:.12e}, expected {expected:.12e}"
                    );
                }
            }
        }
    });
}

#[test]
fn g3_g4_g5_interior_paths_replay_across_tiles_ranges_and_refusals() {
    with_cx(|cx, gate| {
        let density = [1.0];
        let temperature = [1600.0];
        let absorption = field(
            cx,
            &density,
            [-4.0, -4.0, -1.0],
            [8.0, 8.0, 1.0],
            MassAbsorptionSpectrum::constant(0.2).unwrap(),
            Some(&temperature),
        );
        let scattering = ScatteringGrid::try_new(
            cx,
            absorption,
            MassScatteringSpectrum::constant(0.7).unwrap(),
            PhaseFunction::HenyeyGreenstein { g: 0.3 },
        )
        .unwrap();
        let (mut scene, _) = scattering_scene();
        scene.camera.half_tan = 0.12;
        let bindings = [InteriorMedium {
            boundary_primitive: 0,
            medium: PhysicalMedium::Scattering(&scattering),
        }];
        let media = SceneMedia::try_new(cx, &scene, None, &bindings, media_budget()).unwrap();
        let mut policy = settings(DirectStrategy::PowerMis, 24, 6);
        policy.width = 3;
        policy.height = 2;
        let execution = RenderExecutionConfig::try_new(2, 1, 3, 8 << 20, RunId(713)).unwrap();
        for sampler in [Sampler::Iid, Sampler::OwenSobolFullPath] {
            policy.sampler = sampler;
            let serial = render_scene_media(&media, cx, &policy).unwrap();
            assert_eq!(
                render_scene_media_with_execution(&media, cx, &policy, &execution)
                    .unwrap()
                    .film,
                serial
            );
            let mut progressive = SceneMediaFilm::try_new(&media, &policy).unwrap();
            render_scene_media_range(&media, cx, &policy, &mut progressive, 0, 7).unwrap();
            render_scene_media_range_with_execution(
                &media,
                cx,
                &policy,
                &mut progressive,
                7,
                19,
                &execution,
            )
            .unwrap();
            render_scene_media_range(&media, cx, &policy, &mut progressive, 19, 24).unwrap();
            assert_eq!(progressive.film(), &serial);
        }
        let mut progressive = SceneMediaFilm::try_new(&media, &policy).unwrap();
        render_scene_media_range(&media, cx, &policy, &mut progressive, 0, 7).unwrap();
        let before = progressive.clone();
        let changed_density = [2.0];
        let changed = field(
            cx,
            &changed_density,
            [-4.0, -4.0, -1.0],
            [8.0, 8.0, 1.0],
            MassAbsorptionSpectrum::constant(0.2).unwrap(),
            Some(&temperature),
        );
        let changed_bindings = [InteriorMedium {
            boundary_primitive: 0,
            medium: PhysicalMedium::Absorbing(&changed),
        }];
        let changed_media =
            SceneMedia::try_new(cx, &scene, None, &changed_bindings, media_budget()).unwrap();
        assert!(matches!(
            render_scene_media_range(&changed_media, cx, &policy, &mut progressive, 7, 24),
            Err(TracerError::InvalidInput)
        ));
        assert_eq!(progressive, before);
        gate.request();
        assert!(matches!(
            render_scene_media_range(&media, cx, &policy, &mut progressive, 7, 24),
            Err(TracerError::Cancelled)
        ));
        assert_eq!(progressive, before);
    });
}
