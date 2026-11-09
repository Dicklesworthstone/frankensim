//! Full rough-glass detector composition on the original resolved native scene.

use super::*;
use crate::camera::{Aperture, CameraProjection, PhysicalCamera};
use crate::charts::TriMesh;
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnectionOptions,
};
use crate::polarization::dielectric::scene::curved::camera::CurvedCameraScene;
use crate::polarization::dielectric::scene::curved::camera::point_lights::{
    PointLightBranch, PointLightSourcePolarization,
};
use crate::polarization::dielectric::scene::curved::camera::raster::tests::{
    settings as trace_settings, with_cx,
};
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedTraceSettings, SceneMedium,
};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;
use crate::polarization::photometry::{PolarizedSource, SpectralOptics};
use crate::spectral::detector::QuantumEfficiencySensor;
use crate::spectral::{cie_x, cie_y, cie_z, y_integral};
use fs_blake3::hash_domain;
use fs_geom::{Point3, Vec3};

mod sequential;

fn material(index: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}
fn geometry() -> InstanceScene {
    InstanceScene::try_new(
        [(1_u64, 0.0), (8, -2.0), (9, 2.0)]
            .map(|(id, z)| {
                GeometryInstance::try_new(
                    id,
                    hash_domain("rough-glass-detector", &id.to_le_bytes()),
                    SharedGeometry::mesh(TriMesh::new(
                        vec![[-20.0, -20.0, z], [20.0, -20.0, z], [0.0, 20.0, z]],
                        vec![[0, 1, 2]],
                    )),
                    RigidTransform::identity(),
                )
                .unwrap()
            })
            .to_vec(),
    )
    .unwrap()
}
fn bindings() -> Vec<SurfaceBinding> {
    vec![
        SurfaceBinding {
            object_id: 1,
            interaction: Interaction::Dielectric {
                outside: 0,
                inside: 1,
                surface: DielectricSurface::try_anisotropic_rough(
                    0.35,
                    0.65,
                    Vec3::new(1.0, 0.0, 0.0),
                )
                .unwrap(),
            },
        },
        SurfaceBinding {
            object_id: 8,
            interaction: Interaction::Opaque,
        },
        SurfaceBinding {
            object_id: 9,
            interaction: Interaction::Opaque,
        },
    ]
}
fn settings() -> PointLightCameraSettings {
    let surfaces = CurvedTraceSettings {
        max_distance_m: 6.0,
        ..trace_settings()
    };
    let mut connection = BundleConnectionOptions::default();
    connection.scene.trace = surfaces;
    PointLightCameraSettings {
        surfaces,
        connection,
        max_connections: 8,
        max_trace_evaluations: 256,
    }
}
fn lamp(z: f64, angle: f64) -> PolarizedSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.2, 0.0, z]).unwrap(),
        1000.0,
        1.0,
        AngularDistribution::isotropic(),
        SourceSpectrum::try_from_bins(core::array::from_fn(|i| 1.0 + i as f64 / 80.0)).unwrap(),
    )
    .unwrap()
    .with_unpolarized_optics()
    .through(&SpectralOptics::linear_polarizer(angle).unwrap())
    .unwrap()
}
fn branch(source: &PhotometricSource, medium: usize) -> PointLightBranch<'_> {
    PointLightBranch {
        source,
        receiver_object_id: 1,
        initial_medium: medium,
        initial_direction: Vec3::new(0.0, 0.0, if medium == 0 { -1.0 } else { 1.0 }),
        source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
        policy: BundleBranchPolicy::Sequence(&[]),
        position_tolerance_m: 1e-6,
    }
}
fn analyzer() -> SpectralCameraAnalyzer {
    SpectralCameraAnalyzer::try_new(
        Vec3::new(1.0, 0.0, 0.0),
        vec![
            (
                LAMBDA_MIN,
                Stokes {
                    i: 0.45,
                    q: 0.1,
                    u: -0.06,
                    v: 0.03,
                },
            ),
            (
                LAMBDA_MAX,
                Stokes {
                    i: 0.65,
                    q: -0.2,
                    u: 0.08,
                    v: -0.04,
                },
            ),
        ],
    )
    .unwrap()
}
fn sensor() -> QuantumEfficiencySensor<3> {
    // Constant sum gives the same four wavelengths as one uniform channel;
    // the two live channels still have different spectral response curves.
    QuantumEfficiencySensor::try_new(vec![
        (LAMBDA_MIN, [0.1, 0.9, 0.0]),
        (LAMBDA_MAX, [0.9, 0.1, 0.0]),
    ])
    .unwrap()
}
fn qe(wavelength: f64) -> [f64; 3] {
    let first = 0.1 + 0.8 * (wavelength - LAMBDA_MIN) / (LAMBDA_MAX - LAMBDA_MIN);
    [first, 1.0 - first, 0.0]
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1e-12),
        "{actual:e} != {expected:e}"
    );
}
fn check_work<P>(pixel: &PointLightPixelSample<P>, lane: usize, sample: &PointLightCameraSample) {
    assert_eq!(
        pixel.fixture_radiance_w_m2_sr_nm()[lane],
        sample.fixture_radiance_w_m2_sr_nm()
    );
    assert_eq!(pixel.connections()[lane], u64::from(sample.connections()));
    assert_eq!(
        pixel.trace_evaluations()[lane],
        u64::from(sample.trace_evaluations())
    );
    assert_eq!(pixel.field_evaluations()[lane], sample.field_evaluations());
    assert_eq!(pixel.query_work()[lane], sample.query_work());
}

#[test]
fn g1_full_glass_detector_keeps_both_polarized_source_media_native_lanes_and_photon_units() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let materials = [material(1.0, 0.12), material(1.5, 0.24)];
        let media = [
            SceneMedium::Homogeneous(&materials[0]),
            SceneMedium::Homogeneous(&materials[1]),
        ];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let camera = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        let polarized = [lamp(1.0, 0.17), lamp(-1.0, 0.83)];
        let branches = [
            branch(polarized[0].intensity_source(), 0),
            branch(polarized[1].intensity_source(), 1),
        ];
        let polarization = polarized.each_ref().map(|source| {
            PointLightSourcePolarization::try_new(source.intensity_source(), source).unwrap()
        });
        let points = CurvedPointLightCameraScene::try_new(cx, &camera, &branches)
            .unwrap()
            .with_source_polarizations(cx, &polarization)
            .unwrap();
        let smith = SmithDielectricPointLightCameraScene::new(
            &points,
            SmithWalkSettings::try_new(256).unwrap(),
        );
        let analyzer = analyzer();
        let sensor = sensor();
        let one = QuantumEfficiencySensor::try_new(vec![(LAMBDA_MIN, [1.0]), (LAMBDA_MAX, [1.0])])
            .unwrap();
        let (seed, pixel, sample) = (0x1234_5678_0000_0049_u64, 0, 9);
        let render_cx = cx.with_stream_seed(seed);
        for medium in [0, 1] {
            let sign = if medium == 0 { 1.0 } else { -1.0 };
            let physical = PhysicalCamera::try_look_at(
                Point3::new(0.0, 0.0, 0.5 * sign),
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                CameraProjection::try_half_tangent(0.08).unwrap(),
                1.0,
                Aperture::try_circular(0.0).unwrap(),
            )
            .unwrap();
            let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
            let result = raster
                .sample_pixel_sensor_smith_dielectric_point_lights_spectrally_analyzed(
                    cx,
                    &smith,
                    &sensor,
                    medium,
                    &analyzer,
                    settings(),
                    seed,
                    pixel,
                    sample,
                )
                .unwrap();
            let repeat = raster
                .sample_pixel_sensor_smith_dielectric_point_lights_spectrally_analyzed(
                    cx,
                    &smith,
                    &sensor,
                    medium,
                    &analyzer,
                    settings(),
                    seed,
                    pixel,
                    sample,
                )
                .unwrap();
            assert_eq!(result, repeat);
            let single = raster
                .sample_pixel_sensor_smith_dielectric_point_lights_spectrally_analyzed(
                    cx,
                    &smith,
                    &one,
                    medium,
                    &analyzer,
                    settings(),
                    seed,
                    pixel,
                    sample,
                )
                .unwrap();
            assert_eq!(
                result.pixel_sample().spectral_samples(),
                single.pixel_sample().spectral_samples()
            );
            assert_eq!(result.trace_evaluations(), single.trace_evaluations());
            assert_eq!(result.field_evaluations(), single.field_evaluations());
            assert_eq!(result.query_work(), single.query_work());
            let mut energy = [0.0; 3];
            let mut electrons = [0.0; 3];
            for lane in 0..4 {
                let wavelength = result.pixel_sample().detector_wavelengths()[lane].wavelength_nm();
                let ray = result.pixel_sample().camera_rays()[lane];
                // Explicit original lane/path addresses, retaining both source
                // bindings and every complete microscopic continuation stream.
                let lane_seed = seed ^ (0x6c67_7300_u64 << 32) ^ ((lane as u64) << 30);
                let path = (u64::from(pixel) << 32) | u64::from(sample);
                let expected = smith
                    .sample_ray_analyzed(
                        &render_cx,
                        ray,
                        wavelength,
                        medium,
                        analyzer.at(wavelength).unwrap(),
                        settings(),
                        lane_seed,
                        path,
                    )
                    .unwrap();
                assert_eq!(
                    result.pixel_sample().spectral_samples()[lane],
                    *expected.radiance_sample()
                );
                check_work(&result, lane, &expected);
                assert_eq!(result.connections()[lane], 2);
                let radiance = expected.radiance_sample().radiance_w_m2_sr_nm();
                assert!(radiance > 0.0);
                for channel in 0..3 {
                    let value =
                        radiance * qe(wavelength)[channel] * (LAMBDA_MAX - LAMBDA_MIN) / 4.0;
                    energy[channel] += value;
                    electrons[channel] +=
                        value * wavelength * 1e-9 / (6.62607015e-34 * 299_792_458.0);
                }
            }
            for channel in 0..3 {
                near(
                    result.pixel_sample().qe_weighted_radiance_w_m2_sr()[channel],
                    energy[channel],
                );
                near(
                    result.pixel_sample().photoelectron_rate_per_m2_sr()[channel],
                    electrons[channel],
                );
            }

            let visible = raster
                .sample_pixel_smith_dielectric_point_lights_spectrally_analyzed(
                    cx,
                    &smith,
                    medium,
                    &analyzer,
                    settings(),
                    seed,
                    pixel,
                    sample,
                )
                .unwrap();
            let mut xyz = [0.0; 3];
            for lane in 0..4 {
                let wavelength = visible.pixel_sample().wavelengths_nm()[lane];
                let expected = smith
                    .sample_ray_analyzed(
                        &render_cx,
                        visible.pixel_sample().camera_rays()[lane],
                        wavelength,
                        medium,
                        analyzer.at(wavelength).unwrap(),
                        settings(),
                        seed ^ (0x6c67_7300_u64 << 32) ^ ((lane as u64) << 30),
                        (u64::from(pixel) << 32) | u64::from(sample),
                    )
                    .unwrap();
                assert_eq!(
                    visible.pixel_sample().spectral_samples()[lane],
                    *expected.radiance_sample()
                );
                check_work(&visible, lane, &expected);
                for (sum, response) in
                    xyz.iter_mut()
                        .zip([cie_x(wavelength), cie_y(wavelength), cie_z(wavelength)])
                {
                    *sum += expected.radiance_sample().radiance_w_m2_sr_nm()
                        * response
                        * (LAMBDA_MAX - LAMBDA_MIN)
                        / 4.0
                        / y_integral();
                }
            }
            for channel in 0..3 {
                near(visible.pixel_sample().xyz()[channel], xyz[channel]);
            }
            let constant = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
            for selected in [None, Some(constant)] {
                let detected = raster
                    .sample_pixel_sensor_smith_dielectric_point_lights(
                        cx,
                        &smith,
                        &one,
                        medium,
                        selected,
                        settings(),
                        seed,
                        pixel,
                        sample,
                    )
                    .unwrap();
                for lane in 0..4 {
                    let wavelength =
                        detected.pixel_sample().detector_wavelengths()[lane].wavelength_nm();
                    let ray = detected.pixel_sample().camera_rays()[lane];
                    let lane_seed = seed ^ (0x6c67_7300_u64 << 32) ^ ((lane as u64) << 30);
                    let path = (u64::from(pixel) << 32) | u64::from(sample);
                    let manual = match selected {
                        Some(a) => smith.sample_ray_analyzed(
                            &render_cx,
                            ray,
                            wavelength,
                            medium,
                            a,
                            settings(),
                            lane_seed,
                            path,
                        ),
                        None => smith.sample_ray(
                            &render_cx,
                            ray,
                            wavelength,
                            medium,
                            settings(),
                            lane_seed,
                            path,
                        ),
                    }
                    .unwrap();
                    assert_eq!(
                        detected.pixel_sample().spectral_samples()[lane],
                        *manual.radiance_sample()
                    );
                    check_work(&detected, lane, &manual);
                }
            }
            let mut exhausted = settings();
            exhausted.max_trace_evaluations = 0;
            assert!(matches!(
                raster.sample_pixel_sensor_smith_dielectric_point_lights(
                    cx, &smith, &sensor, medium, None, exhausted, seed, pixel, sample
                ),
                Err(PointLightRasterError::Transport { lane: 0, .. })
            ));
        }
        let physical = PhysicalCamera::try_look_at(
            Point3::new(0.0, 0.0, 0.5),
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_half_tangent(0.08).unwrap(),
            1.0,
            Aperture::try_circular(0.0).unwrap(),
        )
        .unwrap();
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        gate.request();
        assert!(
            raster
                .sample_pixel_sensor_smith_dielectric_point_lights_spectrally_analyzed(
                    cx,
                    &smith,
                    &sensor,
                    0,
                    &analyzer,
                    settings(),
                    seed,
                    pixel,
                    sample
                )
                .is_err()
        );
    });
}
