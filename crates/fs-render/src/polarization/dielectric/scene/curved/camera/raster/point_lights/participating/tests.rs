use super::*;
use crate::charts::TriMesh;
use crate::polarization::dielectric::scene::curved::bundle::connection::images::{
    BundleImageOptions, BundleImageSeed,
};
use crate::polarization::dielectric::scene::curved::camera::point_lights::participating::images::{
    CurvedRayleighImageLightCameraScene, PointLightImageGroup, PointLightImageReceiver,
};
use crate::polarization::dielectric::scene::curved::participating::transmittance::images::CurvedTransmittanceImageWorkspace;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnectionOptions,
};
use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
use crate::polarization::dielectric::scene::curved::camera::point_lights::participating::MolecularPointLightBranch;
use crate::polarization::dielectric::scene::curved::camera::point_lights::{
    CurvedPointLightCameraScene, PointLightBranch, PointLightCameraSettings,
};
use crate::polarization::dielectric::scene::curved::camera::raster::tests::{
    bindings, camera, geometry, settings, with_cx,
};
use crate::polarization::dielectric::scene::curved::camera::{
    CurvedCameraScene, EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighMedium, CurvedRayleighScene, CurvedRayleighTraceSettings,
};
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedTraceSettings, SceneMedium,
};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::spectral::{SpectralReflectance, cie_x, cie_y, cie_z, y_integral};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_blake3::hash_domain;
use fs_geom::fixtures::SphereChart;
use fs_geom::{Point3, Vec3};

fn material() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(1.25).unwrap(),
        BeerLambertAbsorption::try_constant(0.1).unwrap(),
        GlassProvenance::Custom,
    )
}

fn trace_settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_distance_m: 10.0,
        ..settings()
    }
}

fn point_settings() -> RayleighPointLightCameraSettings {
    let mut connection = BundleConnectionOptions::default();
    connection.scene.trace = trace_settings();
    RayleighPointLightCameraSettings {
        point_lights: PointLightCameraSettings {
            surfaces: trace_settings(),
            connection,
            max_connections: 4,
            max_trace_evaluations: 96,
        },
        max_volume_events: 8,
        max_source_free_flight_draws: 256,
        max_source_candidates: 128,
    }
}

#[test]
fn g3_molecular_fixture_pixels_match_analytic_single_scattering_and_keep_camera_paths() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            GeometryInstance::try_new(
                9,
                hash_domain("molecular-point-pixel", b"black-enclosure"),
                SharedGeometry::chart(SphereChart {
                    center: Point3::new(0.0, 0.0, 0.0),
                    radius: 1.6,
                }),
                RigidTransform::identity(),
            )
            .unwrap(),
        ])
        .unwrap();
        let material = material();
        let media = [SceneMedium::Homogeneous(&material)];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 9,
                interaction: Interaction::Opaque,
            }],
            100_000,
        )
        .unwrap();
        let surface = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        let fixtures = CurvedPointLightCameraScene::try_new(cx, &surface, &[]).unwrap();
        let spectrum = MassScatteringSpectrum::constant(0.2).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &laws).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
        let source = PhotometricSource::try_new(
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.4, 0.0, 0.0]).unwrap(),
            1000.0,
            1.0,
            AngularDistribution::isotropic(),
            SourceSpectrum::try_from_bins(core::array::from_fn(|i| 1.0 + i as f64 / 80.0)).unwrap(),
        )
        .unwrap();
        let branches = [MolecularPointLightBranch {
            source: &source,
            receiver_medium: 0,
            initial_medium: 0,
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            policy: BundleBranchPolicy::Sequence(&[]),
            position_tolerance_m: 1e-6,
        }];
        let lit = CurvedRayleighPointLightCameraScene::try_new_with_volume_branches(
            cx, &molecular, &fixtures, &branches,
        )
        .unwrap();
        let unbound =
            CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &fixtures).unwrap();
        let image_seeds = [BundleImageSeed {
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            policy: BundleBranchPolicy::Sequence(&[]),
        }; 2];
        let groups = [PointLightImageGroup {
            source: &source,
            receiver: PointLightImageReceiver::MolecularMedium(0),
            initial_medium: 0,
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &image_seeds,
            image_options: BundleImageOptions {
                max_images: 2,
                ..BundleImageOptions::default()
            },
            position_tolerance_m: 1e-6,
        }];
        let images = CurvedRayleighImageLightCameraScene::try_new(cx, &unbound, &groups).unwrap();
        let mut workspace = CurvedTransmittanceImageWorkspace::try_new(2, 2).unwrap();
        let physical = camera(0.0);
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        // This seed gives one real collision in each lane inside the sphere;
        // every second flight exceeds its diameter, whatever phase direction
        // was sampled. The test reconstructs these distances independently.
        let seed = 120;
        let ordinary = raster
            .sample_pixel_rayleigh(
                cx,
                &molecular,
                0,
                CurvedRayleighTraceSettings {
                    surfaces: trace_settings(),
                    max_volume_events: 8,
                },
                seed,
                4,
                9,
            )
            .unwrap();
        let sample = raster
            .sample_pixel_rayleigh_point_lights(cx, &lit, 0, point_settings(), seed, 4, 9)
            .unwrap();
        let discovered = raster
            .sample_pixel_rayleigh_image_lights(
                cx, &images, 0, point_settings(), &mut workspace, seed, 4, 9,
            )
            .unwrap();
        let repeated = raster
            .sample_pixel_rayleigh_image_lights(
                cx, &images, 0, point_settings(), &mut workspace, seed, 4, 9,
            )
            .unwrap();
        assert_eq!(discovered, repeated);
        assert_eq!(discovered.volume_events(), sample.volume_events());
        assert_eq!(discovered.free_flight_draws(), sample.free_flight_draws());
        assert_eq!(discovered.source_free_flight_draws(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(discovered.source_candidates(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(sample.volume_events(), &[1; CAMERA_WAVELENGTHS]);
        assert_eq!(sample.free_flight_draws(), &[2; CAMERA_WAVELENGTHS]);
        assert_eq!(sample.source_free_flight_draws(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(sample.source_candidates(), &[0; CAMERA_WAVELENGTHS]);
        let pixel = sample.pixel_sample();
        assert_eq!(pixel.camera_ray(), ordinary.pixel_sample().camera_ray());
        let mut expected_xyz = [0.0; 3];
        for (lane, &wavelength) in pixel.wavelengths_nm().iter().enumerate() {
            let flight = philox4x32_10([4, 9, 0, 0x6372_7350 + lane as u32], [seed as u32, 0]);
            let distance = -(1.0 - unit(flight[0])).ln() / 0.2;
            assert!((0.05..1.4).contains(&distance));
            let vertex = pixel
                .camera_ray()
                .origin
                .offset(pixel.camera_ray().dir.scale(distance));
            let source_delta = vertex.delta_from(source.position_m());
            let source_distance = source_delta.dot(source_delta).sqrt();
            let source_direction = source_delta.scale(1.0 / source_distance);
            let cosine = -source_direction.dot(pixel.camera_ray().dir);
            let phase_per_sr = 3.0 * (1.0 + cosine * cosine) / (16.0 * core::f64::consts::PI);
            let intensity = source
                .spectral_intensity_w_sr_nm(source_direction, wavelength)
                .unwrap();
            // The selected camera collision already supplies beta_s and its
            // survival probability. Only source survival appears explicitly;
            // this direct term also has no continuation-phase PDF divisor.
            let expected = phase_per_sr * intensity / (source_distance * source_distance)
                * (-0.1 * distance - 0.3 * source_distance).exp();
            let actual = pixel.spectral_samples()[lane].radiance_w_m2_sr_nm();
            assert!((actual - expected).abs() <= 8e-6 * expected + 1e-14);
            assert_eq!(
                actual,
                sample.point_light_sample().fixture_radiance_w_m2_sr_nm()[lane]
            );
            assert_eq!(
                pixel.spectral_samples()[lane].path(),
                ordinary.pixel_sample().spectral_samples()[lane].path()
            );
            assert_eq!(sample.point_light_sample().connections()[lane], 1);
            assert!(sample.point_light_sample().trace_evaluations()[lane] >= 2);
            // Two seeds converge to one image. Its source power and survival
            // occur once, while both geometric searches remain in the work.
            let image_sample = &discovered.pixel_sample().spectral_samples()[lane];
            assert!((image_sample.radiance_w_m2_sr_nm() - expected).abs()
                <= 8e-6 * expected + 1e-14);
            assert_eq!(image_sample.path(), pixel.spectral_samples()[lane].path());
            assert_eq!(discovered.point_light_sample().connections()[lane], 1);
            assert!(discovered.point_light_sample().trace_evaluations()[lane]
                > sample.point_light_sample().trace_evaluations()[lane]);
            for (sum, observer) in expected_xyz.iter_mut().zip([
                cie_x(wavelength),
                cie_y(wavelength),
                cie_z(wavelength),
            ]) {
                *sum += expected * observer * 100.0 / y_integral();
            }
        }
        for (actual, expected) in pixel.xyz().into_iter().zip(expected_xyz) {
            assert!((actual - expected).abs() <= 8e-6 * expected + 1e-14);
        }
        for (actual, expected) in discovered.pixel_sample().xyz().into_iter().zip(expected_xyz) {
            assert!((actual - expected).abs() <= 8e-6 * expected + 1e-14);
        }
        let mut undersized = CurvedTransmittanceImageWorkspace::try_new(1, 1).unwrap();
        assert!(matches!(
            raster.sample_pixel_rayleigh_image_lights(
                cx, &images, 0, point_settings(), &mut undersized, seed, 4, 9,
            ),
            Err(PointLightRasterError::Transport { lane: 0, .. })
        ));
    });
}

#[test]
fn g3_empty_fixtures_preserve_original_molecular_pixel_and_camera_flight_streams() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let material = material();
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [SurfaceEmitter::new(
            7,
            RadianceSpectrum::try_constant(2.0).unwrap(),
            EmissionSides::Both,
        )];
        let surface = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let fixtures = CurvedPointLightCameraScene::try_new(cx, &surface, &[]).unwrap();
        let spectrum = MassScatteringSpectrum::constant(1e-12).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &laws).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
        let lit = CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &fixtures).unwrap();
        let physical = camera(0.01);
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        let ordinary = raster
            .sample_pixel_rayleigh(
                cx,
                &molecular,
                0,
                CurvedRayleighTraceSettings {
                    surfaces: trace_settings(),
                    max_volume_events: 8,
                },
                73,
                4,
                9,
            )
            .unwrap();
        let sample = raster
            .sample_pixel_rayleigh_point_lights(cx, &lit, 0, point_settings(), 73, 4, 9)
            .unwrap();
        assert_eq!(sample.pixel_sample(), ordinary.pixel_sample());
        assert_eq!(sample.volume_events(), ordinary.volume_events());
        assert_eq!(sample.free_flight_draws(), ordinary.free_flight_draws());
        assert_eq!(sample.source_free_flight_draws(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(sample.source_candidates(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(
            sample.point_light_sample().connections(),
            &[0; CAMERA_WAVELENGTHS]
        );
    });
}

#[test]
fn g1_spectral_pixels_apply_source_survival_once_and_preserve_conditional_camera_weight() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            GeometryInstance::try_new(
                7,
                hash_domain("participating-point-pixel", b"diffuse-plane"),
                SharedGeometry::mesh(TriMesh::new(
                    vec![[-3.0, -3.0, 2.0], [3.0, -3.0, 2.0], [0.0, 3.0, 2.0]],
                    vec![[0, 1, 2]],
                )),
                RigidTransform::identity(),
            )
            .unwrap(),
            GeometryInstance::try_new(
                9,
                hash_domain("participating-point-pixel", b"opaque-enclosure"),
                SharedGeometry::chart(SphereChart {
                    center: Point3::new(0.0, 0.0, 0.0),
                    radius: 4.0,
                }),
                RigidTransform::identity(),
            )
            .unwrap(),
        ])
        .unwrap();
        let material = material();
        let media = [SceneMedium::Homogeneous(&material)];
        let rho = 0.37;
        let reflectances = [SpectralReflectance::constant(rho).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![
                SurfaceBinding {
                    object_id: 7,
                    interaction: Interaction::DepolarizingDiffuse {
                        medium: 0,
                        reflectance: 0,
                    },
                },
                SurfaceBinding {
                    object_id: 9,
                    interaction: Interaction::Opaque,
                },
            ],
            100_000,
        )
        .unwrap();
        let surface = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        let source = PhotometricSource::try_new(
            RigidTransform::identity(),
            1000.0,
            1.0,
            AngularDistribution::isotropic(),
            SourceSpectrum::try_from_bins(core::array::from_fn(|i| 1.0 + i as f64 / 80.0)).unwrap(),
        )
        .unwrap();
        let branches = [PointLightBranch {
            source: &source,
            receiver_object_id: 7,
            initial_medium: 0,
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            policy: BundleBranchPolicy::TransmitOrReflectTir,
            position_tolerance_m: 1e-6,
        }];
        let fixtures = CurvedPointLightCameraScene::try_new(cx, &surface, &branches).unwrap();
        let spectrum = MassScatteringSpectrum::table(&[(380.0, 0.01), (780.0, 0.04)]).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &laws).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
        let lit = CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &fixtures).unwrap();
        let physical = camera(0.0);
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        // Select a deterministic fixture whose two camera flights in every
        // lane exceed the complete 10 m budget even at the largest rate.
        // This isolates conditional camera weight, not a Monte Carlo mean.
        let seed = (0..256_u64)
            .find(|&seed| {
                (0..CAMERA_WAVELENGTHS).all(|lane| {
                    (0..2).all(|flight| {
                        let block = philox4x32_10(
                            [4, 9, flight, 0x6372_7350 + lane as u32],
                            [seed as u32, (seed >> 32) as u32],
                        );
                        -(1.0 - unit(block[0])).ln() / 0.04 > 10.0
                    })
                })
            })
            .unwrap();
        let ordinary = raster
            .sample_pixel_rayleigh(
                cx,
                &molecular,
                0,
                CurvedRayleighTraceSettings {
                    surfaces: trace_settings(),
                    max_volume_events: 8,
                },
                seed,
                4,
                9,
            )
            .unwrap();
        let sample = raster
            .sample_pixel_rayleigh_point_lights(cx, &lit, 0, point_settings(), seed, 4, 9)
            .unwrap();
        let repeat = raster
            .sample_pixel_rayleigh_point_lights(cx, &lit, 0, point_settings(), seed, 4, 9)
            .unwrap();
        assert_eq!(sample, repeat);
        assert_eq!(sample.volume_events(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(sample.free_flight_draws(), &[2; CAMERA_WAVELENGTHS]);
        assert_eq!(sample.source_candidates(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(sample.source_free_flight_draws(), &[0; CAMERA_WAVELENGTHS]);
        let pixel = sample.pixel_sample();
        let direction = pixel.camera_ray().dir;
        let distance = 2.0 / direction.z;
        let mut expected_xyz = [0.0; 3];
        for (lane, &wavelength) in pixel.wavelengths_nm().iter().enumerate() {
            let intensity = source
                .spectral_intensity_w_sr_nm(direction, wavelength)
                .unwrap();
            let beta = 0.01 + 0.03 * (wavelength - 380.0) / 400.0;
            // Both legs have true absorption. Only the deterministic source
            // leg receives explicit scattering survival; camera survival is
            // already in the probability of the two sampled free flights.
            let expected = rho * intensity * direction.z
                / (core::f64::consts::PI * distance * distance)
                * (-(0.2 + beta) * distance).exp();
            let actual = pixel.spectral_samples()[lane].radiance_w_m2_sr_nm();
            assert!((actual - expected).abs() <= 3e-6 * expected + 1e-14);
            assert_eq!(
                actual,
                sample.point_light_sample().fixture_radiance_w_m2_sr_nm()[lane]
            );
            assert_eq!(
                pixel.spectral_samples()[lane].path(),
                ordinary.pixel_sample().spectral_samples()[lane].path()
            );
            assert_eq!(sample.point_light_sample().connections()[lane], 1);
            assert!(sample.point_light_sample().trace_evaluations()[lane] >= 2);
            for (sum, observer) in expected_xyz.iter_mut().zip([
                cie_x(wavelength),
                cie_y(wavelength),
                cie_z(wavelength),
            ]) {
                *sum += expected * observer * 100.0 / y_integral();
            }
        }
        for (actual, expected) in pixel.xyz().into_iter().zip(expected_xyz) {
            assert!((actual - expected).abs() <= 3e-6 * expected + 1e-14);
        }
        let mut exhausted = point_settings();
        exhausted.point_lights.max_trace_evaluations = 1;
        assert!(matches!(
            raster.sample_pixel_rayleigh_point_lights(cx, &lit, 0, exhausted, seed, 4, 9),
            Err(PointLightRasterError::Transport { lane: 0, .. })
        ));
    });
}
