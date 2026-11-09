use super::*;
use crate::camera::{
    Aperture, CameraProjection, LensSample, LensSurface, PhysicalCamera, SequentialLens,
};
use crate::charts::TriMesh;
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::gradient_index::event::GradientEventOptions;
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::Stokes;
use crate::polarization::dielectric::fresnel_dielectric_polarized;
use crate::polarization::dielectric::scene::curved::camera::{
    EmissionSides, RadianceSpectrum, SourcePolarization, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::polarization::dielectric::transport::PolarizationFrame;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 13,
                kernel_id: 8,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn glass(n: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(n).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn camera() -> PhysicalCamera {
    let lens = SequentialLens::try_new(
        550.0,
        &[
            LensSurface::try_new(0.08, 20.0, 0.012, Some(glass(1.5, 2.0))).unwrap(),
            LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
        ],
    )
    .unwrap();
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
        1.0,
        Aperture::try_circular(0.001).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(lens)
    .unwrap()
}

fn plane(id: u64, z: f64) -> GeometryInstance {
    // The outward normal points toward the sensor, so +z enters inside.
    GeometryInstance::try_new(
        id,
        fs_blake3::hash_domain("sequential-camera-source", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![[-10.0, -10.0, z], [0.0, 10.0, z], [10.0, -10.0, z]],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn opaque(id: u64) -> SurfaceBinding {
    SurfaceBinding {
        object_id: id,
        interaction: Interaction::Opaque,
    }
}

fn emitter() -> SurfaceEmitter {
    SurfaceEmitter::new(
        7,
        RadianceSpectrum::try_constant(2.0).unwrap(),
        EmissionSides::Both,
    )
}

fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 8,
        max_distance_m: 4.0,
        max_total_steps: 16384,
        max_query_work: 1000000,
        event: GradientEventOptions {
            surface_tolerance_m: 1e-8,
            ..GradientEventOptions::default()
        },
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn lens_seed_retains_source_polarization_and_actual_scene_beer_distance() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(7, 1.0)]).unwrap();
        let outside = glass(1.0, 0.3);
        let media = [SceneMedium::Homogeneous(&outside)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(7)], 100000).unwrap();
        let physical = camera();
        for row in [
            Stokes::unpolarized(1.0),
            Stokes {
                i: 0.5,
                q: 0.3,
                u: 0.2,
                v: 0.1,
            },
            Stokes {
                i: 0.5,
                q: 0.0,
                u: 0.0,
                v: 0.5,
            },
        ] {
            let optical = physical
                .generate_sequential_spectral_beams_analyzed_from_tangent_offsets(
                    cx,
                    0.15,
                    0.05,
                    LensSample::CENTER,
                    &[550.0],
                    Vec3::new(0.0, 1.0, 0.0),
                    row,
                )
                .unwrap()[0]
                .unwrap();
            let source = Stokes {
                i: 1.0,
                q: 0.2,
                u: -0.3,
                v: 0.4,
            };
            let emitters = [emitter()
                .with_polarization(
                    SourcePolarization::try_constant(optical.frame().horizontal(), source).unwrap(),
                )
                .unwrap()];
            let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
            let result = scene
                .sample_sequential_ray_with_samples(cx, &optical, 0, settings(), |_| {
                    panic!("opaque source needs no sample")
                })
                .unwrap();
            let p = optical.power();
            let distance = (1.0 - optical.ray().origin.z) / optical.ray().dir.z;
            // Opposite travel frames share the horizontal axis. Reciprocal
            // contraction flips U alone; V retains the native physical sign.
            let expected = 2.0
                * (p.i * source.i + p.q * source.q - p.u * source.u + p.v * source.v)
                * (-0.3 * distance).exp();
            near(result.radiance_w_m2_sr_nm(), expected, 2e-7);
            near(result.path().traveled_m(), distance, 2e-7);
            assert_eq!(result.launch_phase_index(), 1.0);
            assert_ne!(optical.geometry_weight(), 1.0);
            // The scene sample excludes the rear-cap geometry multiplier.
            assert!(
                (result.radiance_w_m2_sr_nm() - expected * optical.geometry_weight()).abs() > 1e-5
            );
        }
    });
}

#[test]
fn lens_and_scene_fresnel_are_composed_before_intensity_projection() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(3, 0.5), plane(7, 1.0)]).unwrap();
        let outside = glass(1.0, 0.0);
        let inside = glass(1.5, 0.0);
        let media = [
            SceneMedium::Homogeneous(&outside),
            SceneMedium::Homogeneous(&inside),
        ];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 3,
                    interaction: Interaction::Dielectric {
                        outside: 0,
                        inside: 1,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                opaque(7),
            ],
            100000,
        )
        .unwrap();
        let emitters = [emitter()];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let optical = camera()
            .generate_sequential_spectral_beams_from_tangent_offsets(
                cx,
                0.25,
                0.0,
                LensSample::CENTER,
                &[550.0],
            )
            .unwrap()[0]
            .unwrap();
        let interface = fresnel_dielectric_polarized(optical.ray().dir.z, 1.0, 1.5).unwrap();
        let incidence_frame =
            PolarizationFrame::try_new(optical.ray().dir, Vec3::new(0.0, 1.0, 0.0)).unwrap();
        let power = optical.frame().express_in(optical.power(), incidence_frame);
        assert!(power.q.abs() > 1e-4);
        let survival = 1.0 - interface.scalar.reflectance;
        let transmitted = interface.transmission.apply_diagonal_jones(power).i;
        let result = scene
            .sample_sequential_ray_with_samples(cx, &optical, 0, settings(), |_| [0.5, 0.5, 0.99])
            .unwrap();
        let expected = 2.0 * transmitted / survival / (1.5 * 1.5);
        near(result.radiance_w_m2_sr_nm(), expected, 2e-7);
        let incorrectly_reset = 2.0
            * interface
                .transmission
                .apply_diagonal_jones(Stokes::unpolarized(power.i))
                .i
            / survival
            / (1.5 * 1.5);
        assert!((expected - incorrectly_reset).abs() > 1e-6);
    });
}

#[test]
fn original_exterior_admission_and_cancellation_precede_scene_transport() {
    with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![plane(7, 1.0)]).unwrap();
        let outside = glass(1.1, 0.0);
        let media = [SceneMedium::Homogeneous(&outside)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(7)], 100000).unwrap();
        let emitters = [emitter()];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let optical = camera()
            .generate_sequential_spectral_beams_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[550.0],
            )
            .unwrap()[0]
            .unwrap();
        assert_eq!(
            scene.sample_sequential_ray_with_samples(cx, &optical, 0, settings(), |_| panic!(
                "mismatched medium must refuse first"
            )),
            Err(CameraRadianceError::IncompatibleLensExterior)
        );
        gate.request();
        assert!(matches!(
            scene.sample_sequential_ray(cx, &optical, 0, settings(), 13, 1),
            Err(CameraRadianceError::Scene(_))
        ));
    });
}

#[test]
fn immersion_lens_uses_original_scene_medium_and_sensor_radiance_measure() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(7, 1.0)]).unwrap();
        let outside = glass(1.33, 0.3);
        let media = [SceneMedium::Homogeneous(&outside)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(7)], 100000).unwrap();
        let emitters = [emitter()];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let interface_power = |a: f64, b: f64| 4.0 * a * b / ((a + b) * (a + b));
        for sensor_index in [1.0, 1.2, 1.6] {
            let lens = SequentialLens::try_new_with_sensor_medium(
                550.0,
                Some(glass(sensor_index, 0.2)),
                &[
                    LensSurface::try_new(0.08, 20.0, 0.012, Some(glass(2.0, 0.7))).unwrap(),
                    LensSurface::try_new(0.09, -20.0, 0.012, Some(outside)).unwrap(),
                ],
            )
            .unwrap();
            let physical = camera().with_sequential_lens(lens).unwrap();
            let optical = physical
                .generate_sequential_spectral_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[550.0],
                )
                .unwrap()[0]
                .unwrap();
            let result = scene
                .sample_sequential_ray_with_samples(cx, &optical, 0, settings(), |_| {
                    panic!("opaque source needs no sampling")
                })
                .unwrap();
            // Normal-incidence Fresnel power, Beer loss in three physical
            // segments, and the sensor/source radiance endpoint ratio. The
            // scene exterior is already in the original native medium owner.
            let ratio = sensor_index / 1.33;
            let expected = 2.0
                * interface_power(sensor_index, 2.0)
                * interface_power(2.0, 1.33)
                * f64::exp(-0.2 * 0.08 - 0.7 * 0.01 - 0.3 * 0.91)
                * ratio
                * ratio;
            near(result.radiance_w_m2_sr_nm(), expected, 2e-7);
            near(result.path().traveled_m(), 0.91, 2e-7);
            assert_eq!(result.launch_phase_index(), 1.33);
            assert_eq!(optical.scene_phase_index(), 1.33);
            if sensor_index == 1.6 {
                // A reciprocal radiance weight may exceed one even though
                // every physical power-transmission coefficient is passive.
                assert!(optical.power().i > 1.0);
            }
            let wrong_exterior = glass(1.31, 0.3);
            let wrong_media = [SceneMedium::Homogeneous(&wrong_exterior)];
            let wrong_native =
                CurvedPolarizedScene::try_new(cx, &geometry, &wrong_media, vec![opaque(7)], 100000)
                    .unwrap();
            let wrong = CurvedCameraScene::try_new(cx, &wrong_native, &emitters, None).unwrap();
            assert_eq!(
                wrong.sample_sequential_ray_with_samples(cx, &optical, 0, settings(), |_| {
                    panic!("mismatched authored exterior must refuse before transport")
                }),
                Err(CameraRadianceError::IncompatibleLensExterior)
            );
        }
    });
}

#[cfg(feature = "volumes")]
#[test]
fn zero_scattering_preserves_complete_lens_seed_and_surface_sample_bits() {
    use crate::polarization::dielectric::scene::curved::camera::CameraAnalyzer;
    use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
    use crate::polarization::dielectric::scene::curved::participating::{
        CurvedRayleighMedium, CurvedRayleighScene, CurvedRayleighTraceSettings,
    };
    use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(7, 1.0)]).unwrap();
        let outside = glass(1.0, 0.2);
        let media = [SceneMedium::Homogeneous(&outside)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(7)], 100000).unwrap();
        let emitters = [emitter()];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let optical = camera()
            .generate_sequential_spectral_beams_analyzed_from_tangent_offsets(
                cx,
                0.15,
                0.05,
                LensSample::CENTER,
                &[550.0],
                Vec3::new(0.0, 1.0, 0.0),
                CameraAnalyzer::try_circular(Vec3::new(0.0, 1.0, 0.0), true)
                    .unwrap()
                    .row(),
            )
            .unwrap()[0]
            .unwrap();
        let expected = scene
            .sample_sequential_ray(cx, &optical, 0, settings(), 17, 9)
            .unwrap();
        let coefficient = MassScatteringSpectrum::constant(0.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&coefficient, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let coupled = CurvedRayleighCameraScene::try_new(cx, &scene, &scattering).unwrap();
        let settings = CurvedRayleighTraceSettings {
            surfaces: settings(),
            max_volume_events: 4,
        };
        let result = coupled
            .sample_sequential_ray(cx, &optical, 0, settings, 17, 9)
            .unwrap();
        assert_eq!(*result.radiance_sample(), expected);
        assert_eq!(result.volume_events(), 0);
        assert_eq!(result.free_flight_draws(), 0);
        let explicit = coupled
            .sample_sequential_ray_with_tracking_samples(
                cx,
                &optical,
                0,
                settings,
                |_| panic!("opaque source"),
                |_| panic!("zero scattering"),
            )
            .unwrap();
        assert_eq!(*explicit.radiance_sample(), expected);
    });
}

#[test]
fn first_lte_leg_keeps_analyzer_and_resolved_lens_intensity() {
    use crate::polarization::dielectric::scene::curved::camera::{
        CameraAnalyzer, IsothermalVolumeEmitter,
    };
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(7, 1.0)]).unwrap();
        let host = glass(1.0, 0.3);
        let media = [SceneMedium::Homogeneous(&host)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(7)], 100_000).unwrap();
        let source =
            IsothermalVolumeEmitter::try_new(0, 2300.0, fs_blake3::ContentHash([127; 32])).unwrap();
        let scene = CurvedCameraScene::try_new(cx, &native, &[], None)
            .unwrap()
            .with_isothermal_volume_emitters(cx, &[source])
            .unwrap();
        let mut limits = settings();
        // Resolve the exact-plane analytic source integral independently of
        // the default event-position error allowance.
        limits.event.surface_tolerance_m = 1e-10;
        limits.event.integration.position_tolerance_m = 1e-12;
        let wavelength = 550e-9_f64;
        let (h, c, k) = (6.626_070_15e-34, 299_792_458.0, 1.380_649e-23);
        let blackbody = 2.0 * h * c * c
            / (wavelength.powi(5) * (h * c / (wavelength * k * 2300.0)).exp_m1())
            * 1e-9;
        let ray = crate::charts::Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let ordinary = scene
            .sample_ray_with_samples(cx, ray, 550.0, 0, limits, |_| {
                panic!("opaque endpoint needs no sample")
            })
            .unwrap();
        let identity =
            CameraAnalyzer::try_new(Vec3::new(1.0, 0.0, 0.0), Stokes::unpolarized(1.0)).unwrap();
        assert_eq!(
            ordinary,
            scene
                .sample_ray_analyzed_with_samples(cx, ray, 550.0, 0, identity, limits, |_| panic!(
                    "opaque endpoint needs no sample"
                ),)
                .unwrap()
        );
        near(
            ordinary.radiance_w_m2_sr_nm() / blackbody,
            -(-0.3_f64).exp_m1(),
            1e-9,
        );
        let row = Stokes {
            i: 0.125,
            q: 0.075,
            u: 0.0,
            v: 0.0,
        };
        let analyzer = CameraAnalyzer::try_new(Vec3::new(1.0, 0.0, 0.0), row).unwrap();
        let analyzed = scene
            .sample_ray_analyzed_with_samples(cx, ray, 550.0, 0, analyzer, limits, |_| {
                panic!("opaque endpoint needs no sample")
            })
            .unwrap();
        near(
            analyzed.radiance_w_m2_sr_nm() / blackbody,
            row.i * -(-0.3_f64).exp_m1(),
            1e-9,
        );
        let physical = camera();
        for row in [Stokes::unpolarized(1.0), row] {
            let optical = physical
                .generate_sequential_spectral_beams_analyzed_from_tangent_offsets(
                    cx,
                    0.15,
                    0.05,
                    LensSample::CENTER,
                    &[550.0],
                    Vec3::new(0.0, 1.0, 0.0),
                    row,
                )
                .unwrap()[0]
                .unwrap();
            assert!(optical.power().i > 0.0 && optical.power().i < row.i);
            let actual = scene
                .sample_sequential_ray_with_samples(cx, &optical, 0, limits, |_| {
                    panic!("opaque endpoint needs no sample")
                })
                .unwrap();
            let distance = (1.0 - optical.ray().origin.z) / optical.ray().dir.z;
            let expected = optical.power().i * -(-0.3 * distance).exp_m1();
            near(actual.radiance_w_m2_sr_nm() / blackbody, expected, 1e-9);
            assert_ne!(optical.geometry_weight(), 1.0);
        }
    });
}
