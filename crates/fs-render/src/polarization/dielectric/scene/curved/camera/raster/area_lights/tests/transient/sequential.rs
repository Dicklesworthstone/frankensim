//! Original area-source estimator observed through a dispersive real lens.

use super::*;
use crate::camera::{
    Aperture, CameraProjection, LensSample, LensSurface, PhysicalCamera, SequentialCameraBeam,
    SequentialLens,
};
use crate::polarization::dielectric::scene::curved::camera::CameraAnalyzer;
use crate::spectral::{cie_x, cie_y, cie_z, y_integral};

const SEED: u64 = 79;
const SAMPLE: u32 = 6;
const LIGHT_SPEED: f64 = 299_792_458.0;

fn physical_camera(stop: Option<f64>) -> PhysicalCamera {
    let material = DielectricGlass::new(
        CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(2.0).unwrap(),
        GlassProvenance::Custom,
    );
    let mut surfaces = vec![
        LensSurface::try_new(0.08, 20.0, 0.01, Some(material)).unwrap(),
        LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
        LensSurface::try_new(0.1, 0.0, 0.02, None).unwrap(),
    ];
    if let Some(radius) = stop {
        surfaces.push(LensSurface::try_new(0.12, 0.0, radius, None).unwrap());
    }
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
        1.0,
        Aperture::try_circular(0.01).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(SequentialLens::try_new(550.0, &surfaces).unwrap())
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
        .generate_sequential_spectral_timed_beams_analyzed_from_tangent_offsets(
            cx,
            (2.0 * unit(raster[0]) - 1.0) * half,
            (1.0 - 2.0 * unit(raster[1])) * half,
            LensSample::try_new(unit(pupil[0]), unit(pupil[1])).unwrap(),
            wavelengths,
            analyzer.reference_axis(),
            analyzer.row(),
        )
        .unwrap()
}

fn vacuum() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(1.0).unwrap(),
        BeerLambertAbsorption::try_constant(0.2).unwrap(),
        GlassProvenance::Custom,
    )
}

fn bindings<'a, 's, 'm>(
    source: &'a CurvedAreaEmitter<'s, 'm, 's>,
    seeds: &'a [BundleImageSeed<'a>],
) -> [AreaLightBinding<'a, 's, 'm>; 1] {
    [AreaLightBinding {
        source,
        receiver_object_id: 1,
        source_side: BoundarySide::Positive,
        initial_medium: 0,
        source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
        seeds,
        image_options: BundleImageOptions {
            max_images: 1,
            ..BundleImageOptions::default()
        },
        position_tolerance_m: 2e-7,
    }]
}

#[test]
fn real_lens_area_pixels_match_retarded_plane_radiometry_and_native_lane_replay() {
    fixtures::with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            panel(1, 1.0, 10.0),
            panel(9, -1.0, 0.3),
            panel(10, -1.5, 50.0),
        ])
        .unwrap();
        let vacuum = vacuum();
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let reflectances = [SpectralReflectance::constant(0.6).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::DepolarizingDiffuse {
                        medium: 0,
                        reflectance: 0,
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
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_constant(4.0).unwrap(),
            EmissionSides::Front,
        )];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [BundleImageSeed {
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            policy: BundleBranchPolicy::TransmitOrReflectTir,
        }];
        let bindings = bindings(&source, &seeds);
        let lights = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(30e-9);
        let history = LinearHistory {
            epoch,
            offset: 1.0,
            slope: 1e7,
            calls: AtomicUsize::new(0),
            fail_at: AtomicUsize::new(0),
        };
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let timed = TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&history)])
            .unwrap()
            .with_analyzer(analyzer);
        let raw =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&history)]).unwrap();
        let observer = TransientCurvedAreaLightCameraScene::try_new(cx, &timed, &lights).unwrap();
        let raw_observer = TransientCurvedAreaLightCameraScene::try_new(cx, &raw, &lights).unwrap();
        let physical = physical_camera(None);
        let settings = || scene_settings(&geometry, &physical);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let result = raster
            .sample_pixel_sequential_area_lights_at_time(
                cx,
                &observer,
                0,
                clock,
                settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert_eq!(result.observation_time_s(), clock);
        let pixel = result.pixel_sample();
        let beams = packet(&physical, cx, pixel.wavelengths_nm(), analyzer);
        let key = [SEED as u32, 0];
        let render_cx = cx.with_stream_seed(SEED);
        let mut xyz = [0.0; 3];
        for lane in 0..4 {
            let beam = beams[lane].unwrap();
            assert_eq!(pixel.camera_rays()[lane], Some(beam.ray()));
            assert_eq!(pixel.geometry_weights()[lane], beam.geometry_weight());
            let ray = beam.ray();
            let camera_length = (1.0 - ray.origin.z) / ray.dir.z;
            let receiver = Point3::new(
                ray.origin.x + camera_length * ray.dir.x,
                ray.origin.y + camera_length * ray.dir.y,
                1.0,
            );
            let draw = philox4x32_10([0, SAMPLE, 0, 0x6372_6100 + lane as u32 * 64], key);
            let root = unit(draw[0]).sqrt();
            let v = unit(draw[1]);
            // Independent triangle coordinates, homogeneous Beer transport,
            // projected solid angle, cosine BSDF PDF and power-heuristic MIS.
            let point = Point3::new(
                -0.3 + 0.6 * root * (1.0 - v) + 0.3 * root * v,
                -0.3 + 0.9 * root * v,
                -1.0,
            );
            let displacement = receiver.delta_from(point);
            let distance_squared = displacement.dot(displacement);
            let distance = distance_squared.sqrt();
            let cosine = 2.0 / distance;
            let delay = beam.group_delay_s().unwrap()
                + Dd::from_f64(camera_length + distance) / Dd::from_f64(LIGHT_SPEED);
            let scale = 1.0 + 1e7 * ((clock - epoch) - delay).to_f64();
            let irradiance = 4.0 * 0.27 * cosine * cosine / distance_squared;
            let p_light = distance_squared / (0.27 * cosine);
            let p_bsdf = cosine / core::f64::consts::PI;
            let mis = p_light * p_light / (p_light * p_light + p_bsdf * p_bsdf);
            let expected = beam.power().i * (-0.2 * (camera_length + distance)).exp() * 0.6
                / core::f64::consts::PI
                * irradiance
                * mis
                * scale;
            near(
                result.area_sample().area_radiance_w_m2_sr_nm()[lane],
                expected,
            );
            assert_eq!(result.area_sample().area_samples()[lane], 1);
            assert_eq!(result.area_sample().images()[lane], 1);

            let manual = raw_observer
                .sample_sequential_ray_with_samples_at_time(
                    &render_cx,
                    &beam,
                    0,
                    clock,
                    settings(),
                    &mut workspace,
                    |event| {
                        let block =
                            philox4x32_10([0, SAMPLE, event, 0x6372_7340 + lane as u32], key);
                        core::array::from_fn(|i| unit(block[i]))
                    },
                    |event, binding| {
                        let block = philox4x32_10(
                            [
                                0,
                                SAMPLE,
                                event,
                                0x6372_6100 + lane as u32 * 64 + binding as u32,
                            ],
                            key,
                        );
                        [unit(block[0]), unit(block[1])]
                    },
                )
                .unwrap();
            assert_eq!(
                pixel.spectral_samples()[lane],
                Some(*manual.radiance_sample())
            );
            assert_eq!(
                result.area_sample().trace_evaluations()[lane],
                u64::from(manual.area_sample().trace_evaluations())
            );
            let weight = manual.radiance_sample().radiance_w_m2_sr_nm()
                * beam.geometry_weight()
                * (LAMBDA_MAX - LAMBDA_MIN)
                / 4.0
                / y_integral();
            for (sum, observer) in xyz.iter_mut().zip([
                cie_x(beam.wavelength_nm()),
                cie_y(beam.wavelength_nm()),
                cie_z(beam.wavelength_nm()),
            ]) {
                *sum += weight * observer;
            }
        }
        for channel in 0..3 {
            near(pixel.xyz()[channel], xyz[channel]);
        }

        let window = ExposureWindow::try_new(clock, clock + Dd::from_f64(4e-9)).unwrap();
        let exposed = raster
            .sample_pixel_sequential_area_lights_during_exposure(
                cx,
                &observer,
                0,
                window,
                settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let middle = raster
            .sample_pixel_sequential_area_lights_at_time(
                cx,
                &observer,
                0,
                clock + Dd::from_f64(2e-9),
                settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        for channel in 0..3 {
            near(
                exposed.pixel_sample().mean_pixel_sample().xyz()[channel],
                middle.pixel_sample().xyz()[channel],
            );
        }
        assert_eq!(
            exposed.trace_evaluations(),
            middle.area_sample().trace_evaluations()
        );
    });
}

#[test]
fn qe_real_lens_area_pixels_keep_clipped_probability_mass_and_sensor_admission() {
    fixtures::with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![
            panel(1, 1.0, 10.0),
            panel(9, -1.0, 0.3),
            panel(10, -1.5, 50.0),
        ])
        .unwrap();
        let vacuum = vacuum();
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let reflectances = [SpectralReflectance::constant(0.6).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::DepolarizingDiffuse {
                        medium: 0,
                        reflectance: 0,
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
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_constant(4.0).unwrap(),
            EmissionSides::Front,
        )];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [BundleImageSeed {
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            policy: BundleBranchPolicy::TransmitOrReflectTir,
        }];
        let bindings = bindings(&source, &seeds);
        let lights = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let clock = Dd::from_f64(2.0_f64.powi(30)) + Dd::from_f64(30e-9);
        let history = LinearHistory {
            epoch: clock,
            offset: 1.0,
            slope: 0.0,
            calls: AtomicUsize::new(0),
            fail_at: AtomicUsize::new(0),
        };
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let timed = TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&history)])
            .unwrap()
            .with_analyzer(analyzer);
        let observer = TransientCurvedAreaLightCameraScene::try_new(cx, &timed, &lights).unwrap();
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (LAMBDA_MIN, [0.25, 0.75, 0.0]),
            (LAMBDA_MAX, [0.25, 0.75, 0.0]),
        ])
        .unwrap();
        // Keep the same front plane and its vacuum path in both prescriptions.
        // Only the clear radius changes; absorbing scene length must not change.
        let physical = physical_camera(Some(0.02));
        let settings = || scene_settings(&geometry, &physical);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let full = raster
            .sample_pixel_sensor_sequential_area_lights_at_time(
                cx,
                &observer,
                &sensor,
                0,
                clock,
                settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let wavelengths = full
            .pixel_sample()
            .detector_wavelengths()
            .map(|draw| draw.wavelength_nm());
        let beams = packet(&physical, cx, &wavelengths, analyzer);
        let radii = beams.map(|beam| {
            let ray = beam.unwrap().ray();
            let distance = (0.12 - ray.origin.z) / ray.dir.z;
            (ray.origin.x + distance * ray.dir.x).hypot(ray.origin.y + distance * ray.dir.y)
        });
        let lo = radii.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = radii.iter().copied().fold(0.0, f64::max);
        assert!(hi - lo > 1e-10, "fixture needs chromatic pupil clipping");
        let stop = (lo + hi) * 0.5;
        let clipped_camera = physical_camera(Some(stop));
        let clipped = CurvedCameraRaster::try_new(&clipped_camera, 1, 1)
            .unwrap()
            .sample_pixel_sensor_sequential_area_lights_at_time(
                cx,
                &observer,
                &sensor,
                0,
                clock,
                settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let mut energy = [0.0; 3];
        let mut electrons = [0.0; 3];
        let mut surviving = 0;
        for lane in 0..4 {
            let present = radii[lane] < stop;
            assert_eq!(
                clipped.pixel_sample().spectral_samples()[lane].is_some(),
                present
            );
            if !present {
                assert_eq!(clipped.area_sample().area_samples()[lane], 0);
                assert_eq!(clipped.area_sample().trace_evaluations()[lane], 0);
                continue;
            }
            surviving += 1;
            let original = full.pixel_sample().spectral_samples()[lane].unwrap();
            near(
                clipped.pixel_sample().spectral_samples()[lane]
                    .unwrap()
                    .radiance_w_m2_sr_nm(),
                original.radiance_w_m2_sr_nm(),
            );
            for channel in 0..3 {
                let contribution = original.radiance_w_m2_sr_nm()
                    * beams[lane].unwrap().geometry_weight()
                    * [0.25, 0.75, 0.0][channel]
                    * (LAMBDA_MAX - LAMBDA_MIN)
                    / 4.0;
                energy[channel] += contribution;
                electrons[channel] +=
                    contribution * wavelengths[lane] * 1e-9 / (6.62607015e-34 * LIGHT_SPEED);
            }
        }
        assert!(surviving > 0 && surviving < 4);
        for channel in 0..3 {
            near(
                clipped.pixel_sample().qe_weighted_radiance_w_m2_sr()[channel],
                energy[channel],
            );
            near(
                clipped.pixel_sample().photoelectron_rate_per_m2_sr()[channel],
                electrons[channel],
            );
        }
        let window = ExposureWindow::try_new(clock, clock + Dd::from_f64(4e-9)).unwrap();
        let exposed = raster
            .sample_pixel_sensor_sequential_area_lights_during_exposure(
                cx,
                &observer,
                &sensor,
                0,
                window,
                settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert_eq!(
            exposed.trace_evaluations(),
            full.area_sample().trace_evaluations()
        );
        for channel in 0..3 {
            near(
                exposed.pixel_sample().qe_weighted_energy_j_m2_sr()[channel],
                full.pixel_sample().qe_weighted_radiance_w_m2_sr()[channel] * 4e-9,
            );
            near(
                exposed.pixel_sample().photoelectrons_per_m2_sr()[channel],
                full.pixel_sample().photoelectron_rate_per_m2_sr()[channel] * 4e-9,
            );
        }

        let closed_camera = physical_camera(Some(1e-12));
        let closed = CurvedCameraRaster::try_new(&closed_camera, 1, 1).unwrap();
        let before = history.calls.load(Ordering::Relaxed);
        let dark = closed
            .sample_pixel_sensor_sequential_area_lights_at_time(
                cx,
                &observer,
                &sensor,
                0,
                clock,
                settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert_eq!(dark.observation_time_s(), clock);
        assert_eq!(dark.pixel_sample().camera_rays(), &[None; 4]);
        assert_eq!(
            dark.pixel_sample().qe_weighted_radiance_w_m2_sr(),
            &[0.0; 3]
        );
        assert_eq!(dark.area_sample().area_samples(), &[0; 4]);
        assert_eq!(dark.area_sample().trace_evaluations(), &[0; 4]);
        assert_eq!(history.calls.load(Ordering::Relaxed), before);
        assert!(matches!(
            closed.sample_pixel_sequential_area_lights_at_time(
                cx,
                &observer,
                0,
                Dd::from_f64(f64::NAN),
                settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            ),
            Err(TransientAreaLightRasterError::Source(
                TransientCameraError::InvalidObservationTime
            ))
        ));
        history.fail_at.store(before + 1, Ordering::Relaxed);
        assert!(matches!(
            raster.sample_pixel_sensor_sequential_area_lights_at_time(
                cx,
                &observer,
                &sensor,
                0,
                clock,
                settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            ),
            Err(TransientAreaLightRasterError::Transport {
                lane: 0,
                error: TransientAreaLightCameraError::Source(TransientCameraError::History {
                    object_id: 9,
                    error: SurfaceHistoryError::OutsideTimeDomain,
                }),
            })
        ));
        gate.request();
        assert!(
            closed
                .sample_pixel_sensor_sequential_area_lights_at_time(
                    cx,
                    &observer,
                    &sensor,
                    0,
                    clock,
                    settings(),
                    &mut workspace,
                    SEED,
                    0,
                    SAMPLE,
                )
                .is_err()
        );
    });
}
