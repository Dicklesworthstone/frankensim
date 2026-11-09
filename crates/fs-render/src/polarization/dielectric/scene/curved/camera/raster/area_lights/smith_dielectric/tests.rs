//! Original finite-source partitions and full-glass detector work.

use super::*;
use crate::camera::{
    Aperture, CameraProjection, LensSample, LensSurface, PhysicalCamera, SequentialLens,
};
use crate::charts::TriMesh;
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::gradient_index::event::BoundarySide;
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::bundle::connection::BundleBranchPolicy;
use crate::polarization::dielectric::scene::curved::bundle::connection::images::{
    BundleImageOptions, BundleImageSeed,
};
use crate::polarization::dielectric::scene::curved::camera::area_lights::CurvedAreaEmitter;
use crate::polarization::dielectric::scene::curved::camera::area_lights::camera::AreaLightBinding;
use crate::polarization::dielectric::scene::curved::camera::raster::tests::{
    settings as trace_settings, with_cx,
};
use crate::polarization::dielectric::scene::curved::camera::{
    CurvedCameraScene, EmissionSides, RadianceSpectrum, SourcePolarization, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;
use crate::spectral::detector::QuantumEfficiencySensor;
use crate::spectral::{cie_x, cie_y, cie_z, y_integral};
use fs_blake3::hash_domain;
use fs_geom::{Point3, Vec3};

const SEED: u64 = 79;
const SAMPLE: u32 = 6;
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-11 * b.abs().max(1e-12),
        "{a:e} != {b:e}"
    );
}
fn material(index: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}
fn panel(id: u64, z: f64, radius: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("smith-area-detector", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![
                [-radius, -radius, z],
                [radius, -radius, z],
                [0.0, 2.0 * radius, z],
            ],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}
fn settings() -> AreaLightCameraSettings {
    let mut result = AreaLightCameraSettings::default();
    result.surfaces = trace_settings();
    result.connection.scene.trace = result.surfaces;
    result.max_area_samples = 8;
    result.max_trace_evaluations = 4096;
    result
}
fn physical_camera(stop: f64) -> PhysicalCamera {
    let lens_glass = DielectricGlass::new(
        CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(2.0).unwrap(),
        GlassProvenance::Custom,
    );
    let lens = SequentialLens::try_new(
        550.0,
        &[
            LensSurface::try_new(0.08, 20.0, 0.01, Some(lens_glass)).unwrap(),
            LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
            LensSurface::try_new(0.12, 0.0, stop, None).unwrap(),
        ],
    )
    .unwrap();
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
        1.0,
        Aperture::try_circular(0.01).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(lens)
    .unwrap()
}
fn packet(
    camera: &PhysicalCamera,
    cx: &Cx<'_>,
    wavelengths: &[f64; 4],
    response: Option<AnalyzerResponse<'_>>,
) -> [Option<SequentialCameraBeam>; 4] {
    let key = [SEED as u32, 0];
    let raster = philox4x32_10([0, SAMPLE, 0, 0x6372_7331], key);
    let pupil = philox4x32_10([0, SAMPLE, 0, 0x6372_7332], key);
    let half = camera.raster_vertical_half_tan();
    let x = (2.0 * unit(raster[0]) - 1.0) * half;
    let y = (1.0 - 2.0 * unit(raster[1])) * half;
    let lens = LensSample::try_new(unit(pupil[0]), unit(pupil[1])).unwrap();
    match response {
        Some(response) => camera
            .generate_sequential_spectral_beams_with_detector_rows_from_tangent_offsets(
                cx,
                x,
                y,
                lens,
                wavelengths,
                response.at(wavelengths[0]).unwrap().reference_axis(),
                &wavelengths.map(|w| response.at(w).unwrap().row()),
            )
            .unwrap(),
        None => camera
            .generate_sequential_spectral_beams_from_tangent_offsets(cx, x, y, lens, wavelengths)
            .unwrap(),
    }
}
fn check_work<P>(pixel: &AreaLightPixelSample<P>, lane: usize, sample: &AreaLightCameraSample) {
    assert_eq!(
        pixel.area_radiance_w_m2_sr_nm()[lane],
        sample.area_radiance_w_m2_sr_nm()
    );
    assert_eq!(pixel.area_samples()[lane], u64::from(sample.area_samples()));
    assert_eq!(pixel.images()[lane], u64::from(sample.images()));
    assert_eq!(
        pixel.trace_evaluations()[lane],
        u64::from(sample.trace_evaluations())
    );
    assert_eq!(pixel.field_evaluations()[lane], sample.field_evaluations());
    assert_eq!(pixel.query_work()[lane], sample.query_work());
}

#[test]
fn g1_full_glass_area_detectors_keep_both_source_sides_native_partition_and_photon_weights() {
    with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![
            panel(1, 1.0, 3.0),
            panel(9, -1.0, 0.3),
            panel(10, 3.0, 0.3),
        ])
        .unwrap();
        let materials = [material(1.0, 0.0), material(1.5, 0.2)];
        let media = [
            SceneMedium::Homogeneous(&materials[0]),
            SceneMedium::Homogeneous(&materials[1]),
        ];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        outside: 1,
                        inside: 0,
                        surface: DielectricSurface::try_anisotropic_rough(
                            0.35,
                            0.65,
                            Vec3::new(1.0, 0.0, 0.0),
                        )
                        .unwrap(),
                    },
                },
                SurfaceBinding {
                    object_id: 9,
                    interaction: Interaction::Opaque,
                },
                SurfaceBinding {
                    object_id: 10,
                    interaction: Interaction::Opaque,
                },
            ],
            100_000,
        )
        .unwrap();
        let emitters = [
            SurfaceEmitter::new(
                1,
                RadianceSpectrum::try_constant(0.5).unwrap(),
                EmissionSides::Both,
            ),
            SurfaceEmitter::new(
                9,
                RadianceSpectrum::try_new(vec![(LAMBDA_MIN, 4.0), (LAMBDA_MAX, 6.0)]).unwrap(),
                EmissionSides::Both,
            )
            .with_polarization(
                SourcePolarization::try_constant(
                    Vec3::new(0.0, 1.0, 0.0),
                    Stokes {
                        i: 1.0,
                        q: 0.3,
                        u: 0.2,
                        v: -0.1,
                    },
                )
                .unwrap(),
            ).unwrap(),
            SurfaceEmitter::new(
                10,
                RadianceSpectrum::try_new(vec![(LAMBDA_MIN, 7.0), (LAMBDA_MAX, 3.0)]).unwrap(),
                EmissionSides::Both,
            ),
        ];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let sources = [
            CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap(),
            CurvedAreaEmitter::from_camera(cx, &camera, 10, 1).unwrap(),
        ];
        let seeds = [BundleImageSeed {
            initial_direction: Vec3::new(1.0, 0.0, 0.0),
            policy: BundleBranchPolicy::TransmitOrReflectTir,
        }];
        let selected: [_; 2] = core::array::from_fn(|i| AreaLightBinding {
            source: &sources[i],
            receiver_object_id: 1,
            source_side: if i == 0 {
                BoundarySide::Positive
            } else {
                BoundarySide::Negative
            },
            initial_medium: i,
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds,
            image_options: BundleImageOptions {
                max_images: 1,
                ..BundleImageOptions::default()
            },
            position_tolerance_m: 2e-7,
        });
        let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &selected).unwrap();
        let smith = SmithDielectricAreaLightCameraScene::try_new(
            cx,
            &owner,
            SmithWalkSettings::try_new(256).unwrap(),
        )
        .unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(2, 1, 8).unwrap();
        let analyzer = SpectralCameraAnalyzer::try_new(
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
        .unwrap();
        let constant = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (LAMBDA_MIN, [0.25, 0.75, 0.0]),
            (LAMBDA_MAX, [0.25, 0.75, 0.0]),
        ])
        .unwrap();
        let one = QuantumEfficiencySensor::try_new(vec![(LAMBDA_MIN, [1.0]), (LAMBDA_MAX, [1.0])])
            .unwrap();
        let ideal = PhysicalCamera::try_look_at(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_half_tangent(0.08).unwrap(),
            1.0,
            Aperture::try_circular(0.0).unwrap(),
        )
        .unwrap();
        let raster = CurvedCameraRaster::try_new(&ideal, 1, 1).unwrap();
        let render_cx = cx.with_stream_seed(SEED);
        for response in [
            None,
            Some(AnalyzerResponse::Constant(constant)),
            Some(AnalyzerResponse::Spectral(&analyzer)),
        ] {
            let result = match response {
                Some(AnalyzerResponse::Spectral(a)) => raster
                    .sample_pixel_sensor_smith_dielectric_area_lights_spectrally_analyzed(
                        cx,
                        &smith,
                        &sensor,
                        0,
                        a,
                        settings(),
                        &mut workspace,
                        SEED,
                        0,
                        SAMPLE,
                    ),
                other => raster.sample_pixel_sensor_smith_dielectric_area_lights(
                    cx,
                    &smith,
                    &sensor,
                    0,
                    other.map(|r| r.at(550.0).unwrap()),
                    settings(),
                    &mut workspace,
                    SEED,
                    0,
                    SAMPLE,
                ),
            }
            .unwrap();
            let mut energy = [0.0; 3];
            let mut photons = [0.0; 3];
            for lane in 0..4 {
                let wavelength = result.pixel_sample().detector_wavelengths()[lane].wavelength_nm();
                let ray = result.pixel_sample().camera_rays()[lane];
                let seed = SEED ^ (0x6c67_7300_u64 << 32) ^ ((lane as u64) << 30);
                let expected = match response {
                    Some(a) => smith.sample_ray_analyzed(
                        &render_cx,
                        ray,
                        wavelength,
                        0,
                        a.at(wavelength).unwrap(),
                        settings(),
                        &mut workspace,
                        seed,
                        u64::from(SAMPLE),
                    ),
                    None => smith.sample_ray(
                        &render_cx,
                        ray,
                        wavelength,
                        0,
                        settings(),
                        &mut workspace,
                        seed,
                        u64::from(SAMPLE),
                    ),
                }
                .unwrap();
                assert_eq!(
                    result.pixel_sample().spectral_samples()[lane],
                    *expected.radiance_sample()
                );
                check_work(&result, lane, &expected);
                assert_eq!(result.area_samples()[lane], 2);
                assert_eq!(result.images()[lane], 2);
                assert!(
                    expected.radiance_sample().radiance_w_m2_sr_nm()
                        > expected.area_radiance_w_m2_sr_nm(),
                    "native receiver emission must survive the direct endpoint partition"
                );
                for channel in 0..3 {
                    let value = expected.radiance_sample().radiance_w_m2_sr_nm()
                        * [0.25, 0.75, 0.0][channel]
                        * (LAMBDA_MAX - LAMBDA_MIN)
                        / 4.0;
                    energy[channel] += value;
                    photons[channel] +=
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
                    photons[channel],
                );
            }
        }
        let multi = raster
            .sample_pixel_sensor_smith_dielectric_area_lights_spectrally_analyzed(
                cx,
                &smith,
                &sensor,
                0,
                &analyzer,
                settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let single = raster
            .sample_pixel_sensor_smith_dielectric_area_lights_spectrally_analyzed(
                cx,
                &smith,
                &one,
                0,
                &analyzer,
                settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert_eq!(
            multi.pixel_sample().spectral_samples(),
            single.pixel_sample().spectral_samples()
        );
        assert_eq!(multi.trace_evaluations(), single.trace_evaluations());
        assert_eq!(multi.query_work(), single.query_work());
        let visible = raster
            .sample_pixel_smith_dielectric_area_lights_spectrally_analyzed(
                cx,
                &smith,
                0,
                &analyzer,
                settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
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
                    0,
                    analyzer.at(wavelength).unwrap(),
                    settings(),
                    &mut workspace,
                    SEED ^ (0x6c67_7300_u64 << 32) ^ ((lane as u64) << 30),
                    u64::from(SAMPLE),
                )
                .unwrap();
            assert_eq!(
                visible.pixel_sample().spectral_samples()[lane],
                *expected.radiance_sample()
            );
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

        let physical = physical_camera(0.02);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        for response in [
            None,
            Some(AnalyzerResponse::Constant(constant)),
            Some(AnalyzerResponse::Spectral(&analyzer)),
        ] {
            let result = match response {
                Some(AnalyzerResponse::Spectral(a)) => raster.sample_pixel_sensor_sequential_smith_dielectric_area_lights_spectrally_analyzed(
                    cx, &smith, &sensor, 0, a, settings(), &mut workspace, SEED, 0, SAMPLE),
                other => raster.sample_pixel_sensor_sequential_smith_dielectric_area_lights(
                    cx, &smith, &sensor, 0, other.map(|r| r.at(550.0).unwrap()), settings(), &mut workspace, SEED, 0, SAMPLE),
            }.unwrap();
            let wavelengths = result
                .pixel_sample()
                .detector_wavelengths()
                .map(|draw| draw.wavelength_nm());
            let beams = packet(&physical, cx, &wavelengths, response);
            let mut energy = [0.0; 3];
            let mut photons = [0.0; 3];
            for lane in 0..4 {
                let beam = beams[lane].unwrap();
                let expected = smith
                    .sample_sequential_ray(
                        &render_cx,
                        &beam,
                        0,
                        settings(),
                        &mut workspace,
                        SEED ^ (0x6c67_7300_u64 << 32) ^ ((lane as u64) << 30),
                        u64::from(SAMPLE),
                    )
                    .unwrap();
                assert_eq!(
                    result.pixel_sample().spectral_samples()[lane],
                    Some(*expected.radiance_sample())
                );
                check_work(&result, lane, &expected);
                for channel in 0..3 {
                    let value = expected.radiance_sample().radiance_w_m2_sr_nm()
                        * beam.geometry_weight()
                        * [0.25, 0.75, 0.0][channel]
                        * (LAMBDA_MAX - LAMBDA_MIN)
                        / 4.0;
                    energy[channel] += value;
                    photons[channel] +=
                        value * wavelengths[lane] * 1e-9 / (6.62607015e-34 * 299_792_458.0);
                }
            }
            for channel in 0..3 {
                near(
                    result.pixel_sample().qe_weighted_radiance_w_m2_sr()[channel],
                    energy[channel],
                );
                near(
                    result.pixel_sample().photoelectron_rate_per_m2_sr()[channel],
                    photons[channel],
                );
            }
        }
        let visible = raster
            .sample_pixel_sequential_smith_dielectric_area_lights_spectrally_analyzed(
                cx,
                &smith,
                0,
                &analyzer,
                settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let beams = packet(
            &physical,
            cx,
            visible.pixel_sample().wavelengths_nm(),
            Some(AnalyzerResponse::Spectral(&analyzer)),
        );
        let mut xyz = [0.0; 3];
        for lane in 0..4 {
            let beam = beams[lane].unwrap();
            let expected = smith
                .sample_sequential_ray(
                    &render_cx,
                    &beam,
                    0,
                    settings(),
                    &mut workspace,
                    SEED ^ (0x6c67_7300_u64 << 32) ^ ((lane as u64) << 30),
                    u64::from(SAMPLE),
                )
                .unwrap();
            assert_eq!(
                visible.pixel_sample().spectral_samples()[lane],
                Some(*expected.radiance_sample())
            );
            check_work(&visible, lane, &expected);
            for (sum, response) in xyz.iter_mut().zip([
                cie_x(beam.wavelength_nm()),
                cie_y(beam.wavelength_nm()),
                cie_z(beam.wavelength_nm()),
            ]) {
                *sum += expected.radiance_sample().radiance_w_m2_sr_nm()
                    * beam.geometry_weight()
                    * response
                    * (LAMBDA_MAX - LAMBDA_MIN)
                    / 4.0
                    / y_integral();
            }
        }
        for channel in 0..3 {
            near(visible.pixel_sample().xyz()[channel], xyz[channel]);
        }
        let mut exhausted = settings();
        exhausted.max_trace_evaluations = 0;
        assert!(matches!(
            raster.sample_pixel_sensor_sequential_smith_dielectric_area_lights(
                cx,
                &smith,
                &sensor,
                0,
                None,
                exhausted,
                &mut workspace,
                SEED,
                0,
                SAMPLE
            ),
            Err(AreaLightRasterError::Transport { lane: 0, .. })
        ));
        let closed_camera = physical_camera(1e-12);
        let closed = CurvedCameraRaster::try_new(&closed_camera, 1, 1).unwrap();
        let dark = closed
            .sample_pixel_sensor_sequential_smith_dielectric_area_lights_spectrally_analyzed(
                cx,
                &smith,
                &sensor,
                0,
                &analyzer,
                settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert_eq!(dark.pixel_sample().camera_rays(), &[None; 4]);
        assert_eq!(
            dark.pixel_sample().qe_weighted_radiance_w_m2_sr(),
            &[0.0; 3]
        );
        assert_eq!(dark.area_samples(), &[0; 4]);
        assert_eq!(dark.images(), &[0; 4]);
        assert_eq!(dark.trace_evaluations(), &[0; 4]);
        let incomplete = SpectralCameraAnalyzer::try_new(
            Vec3::new(1.0, 0.0, 0.0),
            vec![
                (500.0, Stokes::unpolarized(0.5)),
                (600.0, Stokes::unpolarized(0.5)),
            ],
        )
        .unwrap();
        assert!(matches!(
            closed.sample_pixel_sensor_sequential_smith_dielectric_area_lights_spectrally_analyzed(
                cx,
                &smith,
                &sensor,
                0,
                &incomplete,
                settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE
            ),
            Err(AreaLightRasterError::Raster(CurvedRasterError::Transport(
                _
            )))
        ));
        gate.request();
        assert!(
            closed
                .sample_pixel_sensor_sequential_smith_dielectric_area_lights(
                    cx,
                    &smith,
                    &sensor,
                    0,
                    None,
                    settings(),
                    &mut workspace,
                    SEED,
                    0,
                    SAMPLE
                )
                .is_err()
        );
    });
}
