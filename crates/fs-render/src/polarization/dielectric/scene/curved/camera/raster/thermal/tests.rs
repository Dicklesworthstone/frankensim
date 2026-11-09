use super::super::tests::{camera, geometry, settings, with_cx};
use super::*;
use crate::camera::{Aperture, CameraProjection, LensSurface, SequentialLens};
use crate::conductor::{
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource, ConductorSurface,
};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::gradient_index::event::BoundarySide;
use crate::polarization::dielectric::scene::curved::camera::thermal::smith::SmithThermalEmitter;
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;
use crate::thermal_surface::{
    CoatingThermalEmission,
    smith::{SmithCoatingThermalEmission, SmithThermalError},
};
use fs_blake3::ContentHash;
use fs_geom::{Point3, Vec3};

const SEED: u64 = 79;
const SAMPLE: u32 = 6;

fn with_scene(cx: &Cx<'_>, index: f64, f: impl FnOnce(&SmithThermalCameraScene<'_, '_, '_>)) {
    let geometry = geometry();
    let host = DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(0.15).unwrap(),
        GlassProvenance::Custom,
    );
    let media = [SceneMedium::Homogeneous(&host)];
    let metal = ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.3, 3.0).unwrap()
        }),
        ConductorSource::try_new(ContentHash([111; 32]), ConductorDataStatus::Representative)
            .unwrap(),
    )
    .unwrap();
    let optics = [metal];
    let native = CurvedPolarizedScene::try_new_with_conductors(
        cx,
        &geometry,
        &media,
        &optics,
        vec![SurfaceBinding {
            object_id: 7,
            interaction: Interaction::RoughConductor {
                optics: 0,
                medium: 0,
                surface: ConductorSurface::try_rough(0.7).unwrap(),
            },
        }],
        100_000,
    )
    .unwrap();
    let camera = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
    let thermal =
        CoatingThermalEmission::try_opaque(&[], &[], metal, 2400.0, ContentHash([112; 32]))
            .unwrap();
    let sources = [SmithThermalEmitter::try_new(
        7,
        BoundarySide::Negative,
        SmithCoatingThermalEmission::try_new(thermal, metal).unwrap(),
    )
    .unwrap()];
    let smith = SmithThermalCameraScene::try_new(
        cx,
        &camera,
        &sources,
        SmithWalkSettings::try_new(1024).unwrap(),
    )
    .unwrap();
    f(&smith);
}

fn sensor() -> QuantumEfficiencySensor<2> {
    QuantumEfficiencySensor::try_new(vec![(LAMBDA_MIN, [0.25, 0.75]), (LAMBDA_MAX, [0.25, 0.75])])
        .unwrap()
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-12 * b.abs().max(1e-12),
        "{a:e} != {b:e}"
    );
}
fn lens(stop: f64) -> PhysicalCamera {
    let glass = DielectricGlass::new(
        CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(2.0).unwrap(),
        GlassProvenance::Custom,
    );
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
        1.0,
        Aperture::try_circular(0.01).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(
        SequentialLens::try_new(
            550.0,
            &[
                LensSurface::try_new(0.08, 20.0, 0.01, Some(glass)).unwrap(),
                LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
                LensSurface::try_new(0.12, 0.0, stop, None).unwrap(),
            ],
        )
        .unwrap(),
    )
    .unwrap()
}
fn packet(
    camera: &PhysicalCamera,
    cx: &Cx<'_>,
    wavelengths: &[f64; 4],
    analyzer: CameraAnalyzer,
) -> [Option<SequentialCameraBeam>; 4] {
    let key = [SEED as u32, 0];
    let raster = philox4x32_10([0, SAMPLE, 0, 0x6372_7331], key);
    let pupil = philox4x32_10([0, SAMPLE, 0, 0x6372_7332], key);
    let half = camera.raster_vertical_half_tan();
    camera
        .generate_sequential_spectral_beams_with_detector_rows_from_tangent_offsets(
            cx,
            (2.0 * unit(raster[0]) - 1.0) * half,
            (1.0 - 2.0 * unit(raster[1])) * half,
            LensSample::try_new(unit(pupil[0]), unit(pupil[1])).unwrap(),
            wavelengths,
            analyzer.reference_axis(),
            &[analyzer.row(); 4],
        )
        .unwrap()
}

#[test]
fn g1_smith_thermal_visible_and_qe_pixels_keep_original_lanes_and_photon_units() {
    with_cx(|cx, _| {
        with_scene(cx, 1.0, |scene| {
            let physical = camera(0.001);
            let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
            let sensor = sensor();
            let one =
                QuantumEfficiencySensor::try_new(vec![(LAMBDA_MIN, [1.0]), (LAMBDA_MAX, [1.0])])
                    .unwrap();
            let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.31).unwrap();
            let render_cx = cx.with_stream_seed(SEED);
            for selected in [None, Some(analyzer)] {
                let visible = raster
                    .sample_pixel_smith_thermal(cx, scene, 0, selected, settings(), SEED, 0, SAMPLE)
                    .unwrap();
                let detected = raster
                    .sample_pixel_sensor_smith_thermal(
                        cx,
                        scene,
                        &sensor,
                        0,
                        selected,
                        settings(),
                        SEED,
                        0,
                        SAMPLE,
                    )
                    .unwrap();
                assert_eq!(
                    detected,
                    raster
                        .sample_pixel_sensor_smith_thermal(
                            cx,
                            scene,
                            &sensor,
                            0,
                            selected,
                            settings(),
                            SEED,
                            0,
                            SAMPLE
                        )
                        .unwrap()
                );
                let scalar = raster
                    .sample_pixel_sensor_smith_thermal(
                        cx,
                        scene,
                        &one,
                        0,
                        selected,
                        settings(),
                        SEED,
                        0,
                        SAMPLE,
                    )
                    .unwrap();
                assert_eq!(detected.spectral_samples(), scalar.spectral_samples());
                let mut xyz = [0.0; 3];
                let mut energy = [0.0; 2];
                let mut electrons = [0.0; 2];
                for lane in 0..4 {
                    let manual = |ray, wavelength| {
                        match selected {
                            Some(a) => scene.sample_ray_analyzed(
                                &render_cx,
                                ray,
                                wavelength,
                                0,
                                a,
                                settings(),
                                SEED ^ (0x6c67_7300_u64 << 32) ^ ((lane as u64) << 30),
                                u64::from(SAMPLE),
                            ),
                            None => scene.sample_ray(
                                &render_cx,
                                ray,
                                wavelength,
                                0,
                                settings(),
                                SEED ^ (0x6c67_7300_u64 << 32) ^ ((lane as u64) << 30),
                                u64::from(SAMPLE),
                            ),
                        }
                        .unwrap()
                    };
                    let wavelength = visible.wavelengths_nm()[lane];
                    let expected = manual(visible.camera_rays()[lane], wavelength);
                    assert_eq!(visible.spectral_samples()[lane], expected);
                    assert!(expected.radiance_w_m2_sr_nm() > 0.0);
                    for (sum, observer) in xyz.iter_mut().zip([
                        cie_x(wavelength),
                        cie_y(wavelength),
                        cie_z(wavelength),
                    ]) {
                        *sum += expected.radiance_w_m2_sr_nm() * observer * 100.0 / y_integral();
                    }
                    let wavelength = detected.detector_wavelengths()[lane].wavelength_nm();
                    let expected = manual(detected.camera_rays()[lane], wavelength);
                    assert_eq!(detected.spectral_samples()[lane], expected);
                    for (channel, qe) in [0.25, 0.75].into_iter().enumerate() {
                        let value = expected.radiance_w_m2_sr_nm() * qe * 100.0;
                        energy[channel] += value;
                        electrons[channel] +=
                            value * wavelength * 1e-9 / (6.62607015e-34 * 299_792_458.0);
                    }
                }
                for i in 0..3 {
                    near(visible.xyz()[i], xyz[i]);
                }
                for i in 0..2 {
                    near(detected.qe_weighted_radiance_w_m2_sr()[i], energy[i]);
                    near(detected.photoelectron_rate_per_m2_sr()[i], electrons[i]);
                }
            }
        })
    });
}

#[test]
fn g1_g4_smith_thermal_lens_keeps_polarization_clipping_and_typed_failure() {
    with_cx(|cx, gate| {
        with_scene(cx, 1.0, |scene| {
            let physical = lens(0.02);
            let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
            let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.31).unwrap();
            let sensor = sensor();
            let visible = raster
                .sample_pixel_sequential_smith_thermal(
                    cx,
                    scene,
                    0,
                    Some(analyzer),
                    settings(),
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap();
            let detected = raster
                .sample_pixel_sensor_sequential_smith_thermal(
                    cx,
                    scene,
                    &sensor,
                    0,
                    Some(analyzer),
                    settings(),
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap();
            let render_cx = cx.with_stream_seed(SEED);
            let wavelengths = detected
                .detector_wavelengths()
                .map(|draw| draw.wavelength_nm());
            let beams = packet(&physical, &render_cx, &wavelengths, analyzer);
            for lane in 0..4 {
                let beam = beams[lane].unwrap();
                let manual = scene
                    .sample_sequential_ray(
                        &render_cx,
                        &beam,
                        0,
                        settings(),
                        SEED ^ (0x6c67_7300_u64 << 32) ^ ((lane as u64) << 30),
                        u64::from(SAMPLE),
                    )
                    .unwrap();
                assert_eq!(detected.spectral_samples()[lane], Some(manual));
                assert_eq!(detected.geometry_weights()[lane], beam.geometry_weight());
                assert_eq!(beam.group_delay_s(), None);
            }
            let visible_beams = packet(&physical, &render_cx, visible.wavelengths_nm(), analyzer);
            for lane in 0..4 {
                let manual = scene
                    .sample_sequential_ray(
                        &render_cx,
                        &visible_beams[lane].unwrap(),
                        0,
                        settings(),
                        lighting_seed(SEED, lane),
                        lighting_path(0, SAMPLE),
                    )
                    .unwrap();
                assert_eq!(visible.spectral_samples()[lane], Some(manual));
            }
            let radii = beams.map(|beam| {
                let ray = beam.unwrap().ray();
                ray.origin.x.hypot(ray.origin.y)
            });
            let lo = radii.into_iter().fold(f64::INFINITY, f64::min);
            let hi = radii.into_iter().fold(0.0_f64, f64::max);
            assert!(hi - lo > 1e-10);
            let stop = 0.5 * (lo + hi);
            let clipped_camera = lens(stop);
            let clipped = CurvedCameraRaster::try_new(&clipped_camera, 1, 1)
                .unwrap()
                .sample_pixel_sensor_sequential_smith_thermal(
                    cx,
                    scene,
                    &sensor,
                    0,
                    Some(analyzer),
                    settings(),
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap();
            let mut energy = [0.0; 2];
            let mut live = 0;
            for lane in 0..4 {
                let present = radii[lane] < stop;
                assert_eq!(clipped.spectral_samples()[lane].is_some(), present);
                if present {
                    live += 1;
                    assert_eq!(
                        clipped.spectral_samples()[lane],
                        detected.spectral_samples()[lane]
                    );
                    for (channel, qe) in [0.25, 0.75].into_iter().enumerate() {
                        energy[channel] += clipped.spectral_samples()[lane]
                            .unwrap()
                            .radiance_w_m2_sr_nm()
                            * clipped.geometry_weights()[lane]
                            * qe
                            * 100.0;
                    }
                }
            }
            assert!((1..4).contains(&live));
            for channel in 0..2 {
                near(
                    clipped.qe_weighted_radiance_w_m2_sr()[channel],
                    energy[channel],
                );
            }
        });
        with_scene(cx, 1.3, |scene| {
            let physical = camera(0.0);
            let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
            assert!(matches!(
                raster.sample_pixel_smith_thermal(cx, scene, 0, None, settings(), SEED, 0, SAMPLE),
                Err(SmithThermalRasterError::Transport {
                    lane: 0,
                    error: SmithThermalCameraError::Thermal {
                        object_id: 7,
                        error: SmithThermalError::VacuumRequired,
                    }
                })
            ));
            gate.request();
            assert!(matches!(
                raster.sample_pixel_sensor_smith_thermal(
                    cx,
                    scene,
                    &sensor(),
                    0,
                    None,
                    settings(),
                    SEED,
                    0,
                    SAMPLE
                ),
                Err(SmithThermalRasterError::Raster(_))
            ));
        });
    });
}
