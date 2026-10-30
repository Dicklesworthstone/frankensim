//! G1/G3/G4/G5: physical cell transport reaches camera and light paths.

use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::charts::TriMesh;
use fs_render::dielectric::{DielectricGlass, DielectricSurface};
use fs_render::lighting::{EnvironmentMap, RectLight};
use fs_render::spectral::{LiftedSpectrum, lift_rgb, xyz_of_spectrum};
use fs_render::thermal_surface::GrayDiffuseThermalSurface;
use fs_render::tracer::{
    AbsorbingFilm, Camera, DirectStrategy, Film, Material, Primitive, Sampler, Scene, Settings,
    Shape, TracerError, render, render_absorbing, render_absorbing_range,
};
use fs_render::volumes::VolumeGrid;
use fs_render::volumes::transport::{
    AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum, VolumeTransportError,
};

#[path = "absorbing_transport/scattering.rs"]
mod scattering;

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

fn budget() -> AbsorbingGridBudget {
    AbsorbingGridBudget {
        max_grid_cells: 64,
        max_steps_per_segment: 64,
    }
}

fn settings(strategy: DirectStrategy) -> Settings {
    Settings {
        width: 1,
        height: 1,
        spp: 8192,
        max_depth: 2,
        sampler: Sampler::OwenSobolFullPath,
        strategy,
        seed: 173,
    }
}

fn empty_scene() -> Scene {
    Scene {
        primitives: Vec::new(),
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 2.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.0,
        },
    }
}

fn environment_scene() -> (Scene, (LiftedSpectrum, f64)) {
    let mut scene = empty_scene();
    let environment = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
    let emission = environment.evaluate(scene.camera.forward).unwrap().emission;
    scene.environment = Some(environment);
    (scene, emission)
}

fn quad(z: f64, half_width: f64, material: Material) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [-half_width, -half_width, z],
                [half_width, -half_width, z],
                [half_width, half_width, z],
                [-half_width, half_width, z],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )),
        material,
        emission: None,
    }
}

fn lit_receiver() -> (Scene, LiftedSpectrum, (LiftedSpectrum, f64)) {
    let reflectance = lift_rgb([0.4; 3]);
    let emission = (lift_rgb([1.0; 3]), 3.0);
    let material = Material::Lambertian { reflectance };
    let mut scene = empty_scene();
    scene.camera.eye.z = 0.5;
    scene.primitives.push(quad(0.0, 10.0, material));
    let mut emitter = quad(2.0, 0.5, material);
    emitter.emission = Some(emission);
    scene.primitives.push(emitter);
    scene.lights.push(RectLight {
        corner: Point3::new(-0.5, -0.5, 2.0),
        edge_u: Vec3::new(1.0, 0.0, 0.0),
        edge_v: Vec3::new(0.0, 1.0, 0.0),
        prim: 1,
        emission,
    });
    (scene, reflectance, emission)
}

fn assert_mean(film: &Film, expected: [f64; 3], tolerance: f64) {
    for (channel, want) in expected.into_iter().enumerate() {
        let got = film.xyz[0][channel] / f64::from(film.spp_done);
        assert!(
            (got / want - 1.0).abs() <= tolerance,
            "channel {channel}: {got:.12e} vs {want:.12e}, tolerance {tolerance}"
        );
    }
}

// Direct SI Planck law, evaluated independently of the production log-domain
// implementation. These ordinary visible wavelengths/temperatures are far
// from its overflow/underflow branches. Per-metre radiance converts to per-nm.
fn planck(temperature_k: f64, wavelength_nm: f64) -> f64 {
    let wavelength_m = wavelength_nm * 1.0e-9;
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0;
    let k = 1.380_649e-23;
    2.0 * h * c * c / (wavelength_m.powi(5) * (h * c / (wavelength_m * k * temperature_k)).exp_m1())
        * 1.0e-9
}

fn shadow_square_integral(sigma: f64) -> f64 {
    let mut integral = 0.0;
    let n = 128;
    let dx = 1.0 / f64::from(n);
    for i in 0..n {
        for j in 0..n {
            let x = -0.5 + (f64::from(i) + 0.5) * dx;
            let y = -0.5 + (f64::from(j) + 0.5) * dx;
            let r2 = x * x + y * y + 4.0;
            let slab_length = 0.25 * r2.sqrt();
            integral += 4.0 / (r2 * r2) * (-sigma * slab_length).exp() * dx * dx;
        }
    }
    integral
}

#[test]
fn g1_camera_spectral_beer_lambert_matches_physical_slab() {
    with_cx(|cx, _| {
        let density = [2.0];
        let absorption =
            MassAbsorptionSpectrum::table(&[(380.0, 0.1), (580.0, 0.5), (780.0, 1.2)]).unwrap();
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 1.0]),
            absorption,
            None,
            budget(),
        )
        .unwrap();
        let (scene, (spectrum, scale)) = environment_scene();
        let film = render_absorbing(&scene, &medium, cx, &settings(DirectStrategy::Mis)).unwrap();
        let expected = xyz_of_spectrum(|lambda| {
            spectrum.eval(lambda) * scale * (-2.0 * absorption.at(lambda).unwrap()).exp()
        });
        assert_mean(&film, expected, 0.005);
    });
}

#[test]
fn g1_absorption_is_applied_to_direct_light_shadow_segments() {
    with_cx(|cx, _| {
        // The camera stops below the medium. Only the receiver-to-light rays
        // cross this 0.5 m slab; a post-render camera filter cannot pass.
        let density = [2.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-2.0, -2.0, 1.0], [4.0, 4.0, 0.5]),
            MassAbsorptionSpectrum::constant(0.7).unwrap(),
            None,
            budget(),
        )
        .unwrap();
        let (scene, reflectance, (emission, scale)) = lit_receiver();
        let mut settings = settings(DirectStrategy::NeeOnly);
        settings.max_depth = 1;
        let direct_reference = render_absorbing(&scene, &medium, cx, &settings).unwrap();
        let unattenuated = render(&scene, cx, &settings).unwrap();
        let integral = shadow_square_integral(1.4);
        let expected = xyz_of_spectrum(|lambda| {
            reflectance.eval(lambda) * emission.eval(lambda) * scale * integral
                / core::f64::consts::PI
        });
        assert_mean(&direct_reference, expected, 0.006);
        assert!(direct_reference.xyz[0][1] < 0.55 * unattenuated.xyz[0][1]);
        // Every strategy has the same complete NEE estimator at the final
        // admitted vertex. There is no BSDF continuation with which to share
        // its weight, even for a single rectangle's legacy surface scene.
        for strategy in [DirectStrategy::Mis, DirectStrategy::PowerMis] {
            settings.strategy = strategy;
            let film = render_absorbing(&scene, &medium, cx, &settings).unwrap();
            assert_eq!(film, direct_reference);
        }
    });
}

#[test]
fn g1_physical_thermal_light_shadow_uses_volume_absorption() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-2.0, -2.0, 1.0], [4.0, 4.0, 0.5]),
            MassAbsorptionSpectrum::constant(0.8).unwrap(),
            None,
            budget(),
        )
        .unwrap();
        let (mut scene, reflectance, _) = lit_receiver();
        scene.primitives[1].material = Material::GrayDiffuseThermal {
            surface: GrayDiffuseThermalSurface::try_new(2000.0, 0.7, ContentHash([85; 32]))
                .unwrap(),
        };
        scene.primitives[1].emission = None;
        scene.lights = Vec::new();
        let mut settings = settings(DirectStrategy::NeeOnly);
        settings.max_depth = 1;
        let film = render_absorbing(&scene, &medium, cx, &settings).unwrap();
        let integral = shadow_square_integral(0.8);
        assert_mean(
            &film,
            xyz_of_spectrum(|lambda| {
                reflectance.eval(lambda) * 0.7 * planck(2000.0, lambda) * integral
                    / core::f64::consts::PI
            }),
            0.006,
        );
    });
}

#[test]
fn g1_reflected_environment_paths_cross_the_same_physical_slab() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-1.0e6, -1.0e6, 1.0], [2.0e6, 2.0e6, 1.0]),
            MassAbsorptionSpectrum::constant(0.4).unwrap(),
            None,
            budget(),
        )
        .unwrap();
        let (mut scene, (emission, scale)) = environment_scene();
        scene.camera.eye.z = 0.25;
        let reflectance = lift_rgb([0.6; 3]);
        scene
            .primitives
            .push(quad(0.0, 10.0, Material::Lambertian { reflectance }));
        let mut settings = settings(DirectStrategy::BsdfOnly);
        settings.spp = 16_384;
        let film = render_absorbing(&scene, &medium, cx, &settings).unwrap();
        // Integrate over the cosine-weighted hemisphere independently of the
        // path sampler: E[T] = integral_0^1 2 mu exp(-tau/mu) dmu.
        let integral: f64 = (0..8192)
            .map(|i| {
                let mu = (f64::from(i) + 0.5) / 8192.0;
                2.0 * mu * (-0.4 / mu).exp() / 8192.0
            })
            .sum();
        assert_mean(
            &film,
            xyz_of_spectrum(|lambda| {
                reflectance.eval(lambda) * emission.eval(lambda) * scale * integral
            }),
            0.012,
        );
    });
}

#[test]
fn g1_hot_cells_emit_in_si_units_and_stop_at_opaque_geometry() {
    with_cx(|cx, _| {
        let density = [2.0];
        let temperatures = [1800.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 1.0]),
            MassAbsorptionSpectrum::constant(0.7).unwrap(),
            Some(&temperatures),
            budget(),
        )
        .unwrap();
        let mut settings = settings(DirectStrategy::Mis);
        settings.max_depth = 1;
        let expected = |length: f64| {
            xyz_of_spectrum(|lambda| planck(1800.0, lambda) * -(-1.4_f64 * length).exp_m1())
        };
        let scene = empty_scene();
        let film = render_absorbing(&scene, &medium, cx, &settings).unwrap();
        assert_mean(&film, expected(1.0), 0.005);
        for (surface_z, observed_length) in [(0.5, 0.5), (1.5, 0.0)] {
            let mut blocked = empty_scene();
            blocked.primitives.push(quad(
                surface_z,
                10.0,
                Material::Lambertian {
                    reflectance: lift_rgb([0.0; 3]),
                },
            ));
            let film = render_absorbing(&blocked, &medium, cx, &settings).unwrap();
            if observed_length == 0.0 {
                assert_eq!(film.xyz, vec![[0.0; 3]]);
            } else {
                assert_mean(&film, expected(observed_length), 0.005);
            }
        }
    });
}

#[test]
fn g5_vacuum_medium_preserves_original_surface_tracer_bits() {
    with_cx(|cx, _| {
        let density = [0.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-2.0, -2.0, 0.0], [4.0; 3]),
            MassAbsorptionSpectrum::constant(0.7).unwrap(),
            None,
            budget(),
        )
        .unwrap();
        let (scene, _, _) = lit_receiver();
        for strategy in [
            DirectStrategy::NeeOnly,
            DirectStrategy::BsdfOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
        ] {
            let mut settings = settings(strategy);
            settings.spp = 256;
            assert_eq!(
                render_absorbing(&scene, &medium, cx, &settings).unwrap(),
                render(&scene, cx, &settings).unwrap(),
            );
        }
    });
}

#[test]
fn g4_g5_progressive_volume_binding_replay_and_refusal_are_transactional() {
    with_cx(|cx, gate| {
        let density = [2.0, 2.0];
        let make = |rho, max_steps_per_segment| {
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1, 1, 2], rho, [-1.0, -1.0, 0.0], [2.0, 2.0, 0.5]),
                MassAbsorptionSpectrum::constant(0.7).unwrap(),
                None,
                AbsorbingGridBudget {
                    max_steps_per_segment,
                    ..budget()
                },
            )
            .unwrap()
        };
        let medium = make(&density, 64);
        let (scene, _) = environment_scene();
        let mut settings = settings(DirectStrategy::Mis);
        settings.spp = 256;
        let complete = render_absorbing(&scene, &medium, cx, &settings).unwrap();
        let mut film = AbsorbingFilm::try_new(&medium, &settings).unwrap();
        render_absorbing_range(&scene, &medium, cx, &settings, &mut film, 0, 73).unwrap();
        let checkpoint = film.clone();
        render_absorbing_range(&scene, &medium, cx, &settings, &mut film, 73, 256).unwrap();
        assert_eq!(film.film(), &complete);

        let mut film = checkpoint.clone();
        let changed_density = [3.0, 3.0];
        let changed = make(&changed_density, 64);
        assert_eq!(
            render_absorbing_range(&scene, &changed, cx, &settings, &mut film, 73, 256),
            Err(TracerError::InvalidInput),
        );
        let changed_settings = Settings {
            seed: settings.seed + 1,
            ..settings
        };
        assert_eq!(
            render_absorbing_range(&scene, &medium, cx, &changed_settings, &mut film, 73, 256),
            Err(TracerError::InvalidInput),
        );
        let insufficient_budget = make(&density, 1);
        assert_eq!(
            render_absorbing_range(
                &scene,
                &insufficient_budget,
                cx,
                &settings,
                &mut film,
                73,
                256
            ),
            Err(TracerError::Volume(VolumeTransportError::Budget)),
        );
        assert_eq!(film, checkpoint);
        gate.request();
        assert_eq!(
            render_absorbing_range(&scene, &medium, cx, &settings, &mut film, 73, 256),
            Err(TracerError::Cancelled),
        );
        assert_eq!(film, checkpoint);
    });
}

#[test]
fn g0_nonzero_volume_refuses_unintegrated_dielectric_connectors() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 0.0], [2.0; 3]),
            MassAbsorptionSpectrum::constant(0.7).unwrap(),
            None,
            budget(),
        )
        .unwrap();
        let (mut scene, _) = environment_scene();
        scene.primitives.push(quad(
            0.0,
            1.0,
            Material::Dielectric {
                glass: DielectricGlass::representative_borosilicate(),
                surface: DielectricSurface::SMOOTH,
            },
        ));
        assert_eq!(
            render_absorbing(&scene, &medium, cx, &settings(DirectStrategy::Mis)).unwrap_err(),
            TracerError::VolumeDielectricUnsupported,
        );
    });
}
