//! Actual full-Smith source paths retain clocks, shutter energy and admission.

use super::*;
use crate::camera::LensSample;
use crate::polarization::dielectric::scene::curved::camera::CameraAnalyzer;
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;

fn receiver() -> SurfaceBinding {
    SurfaceBinding {
        object_id: 1,
        interaction: Interaction::Dielectric {
            outside: 1,
            inside: 0,
            surface: DielectricSurface::try_anisotropic_rough(0.2, 0.4, Vec3::new(1.0, 0.0, 0.0))
                .unwrap(),
        },
    }
}

fn binding<'a, 's, 'm>(
    source: &'a CurvedAreaEmitter<'s, 'm, 's>,
    seeds: &'a [BundleImageSeed<'a>],
    transmitted: bool,
) -> AreaLightBinding<'a, 's, 'm> {
    AreaLightBinding {
        source,
        receiver_object_id: 1,
        source_side: if transmitted {
            BoundarySide::Negative
        } else {
            BoundarySide::Positive
        },
        initial_medium: usize::from(transmitted),
        source_transverse_axis: Vec3::new(1.0, 0.0, 0.0),
        seeds,
        image_options: image_options(1),
        position_tolerance_m: 2e-7,
    }
}

#[test]
fn g3_smith_source_histories_preserve_native_paths_and_use_both_bulk_and_lens_group_clocks() {
    with_cx(|cx| {
        let air = glass(1.0, 0.1);
        let material = DielectricGlass::new(
            CauchyIor::try_new(1.4, 0.01, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&material),
        ];
        let physical = sequential::physical_camera();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(40e-9);
        assert_eq!(clock.to_f64(), epoch.to_f64());
        for transmitted in [false, true] {
            let source_z = if transmitted { 3.0 } else { -1.0 };
            let geometry =
                InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, source_z, 0.3)])
                    .unwrap();
            let scene = CurvedPolarizedScene::try_new(
                cx,
                &geometry,
                &media,
                vec![receiver(), opaque(9)],
                100_000,
            )
            .unwrap();
            let emitters = [SurfaceEmitter::new(
                9,
                RadianceSpectrum::try_constant(LE).unwrap(),
                EmissionSides::Both,
            )];
            let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
            let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
            let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
            let bindings = [binding(&light, &seeds, transmitted)];
            let lights = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
            let history = History::new(|_| Ok(1.0));
            let timed =
                TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&history)])
                    .unwrap();
            let source = TransientCurvedAreaLightCameraScene::try_new(cx, &timed, &lights).unwrap();
            let observer = source
                .with_smith_dielectric(cx, SmithWalkSettings::try_new(64).unwrap())
                .unwrap();
            let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
            for wavelength in [450.0, 650.0] {
                let optical = sequential::optical(&physical, cx, wavelength);
                for lens in [false, true] {
                    history.queries.lock().unwrap().clear();
                    let baseline = if lens {
                        observer.source_area_camera().sample_sequential_ray(
                            cx,
                            &optical,
                            0,
                            settings(),
                            &mut workspace,
                            71,
                            5,
                        )
                    } else {
                        observer.source_area_camera().sample_ray(
                            cx,
                            camera_ray(),
                            wavelength,
                            0,
                            settings(),
                            &mut workspace,
                            71,
                            5,
                        )
                    }
                    .unwrap();
                    let actual = if lens {
                        observer.sample_sequential_ray_at_time(
                            cx,
                            &optical,
                            0,
                            clock,
                            settings(),
                            &mut workspace,
                            71,
                            5,
                        )
                    } else {
                        observer.sample_ray_at_time(
                            cx,
                            camera_ray(),
                            wavelength,
                            0,
                            clock,
                            settings(),
                            &mut workspace,
                            71,
                            5,
                        )
                    }
                    .unwrap();
                    assert_eq!(actual.area_sample(), &baseline);
                    assert_eq!(actual.observation_time_s(), clock);
                    assert_eq!((baseline.area_samples(), baseline.images()), (1, 1));
                    assert!(baseline.area_radiance_w_m2_sr_nm() > 0.0);
                    let queries = history.queries.lock().unwrap();
                    assert!(!queries.is_empty());
                    let group_index = if transmitted {
                        1.4 + 0.03 / (wavelength * 1e-3).powi(2)
                    } else {
                        1.0
                    };
                    // Independent straight macroscopic legs; the Smith height
                    // variable contributes no physical distance. The resolved
                    // lens adds 9 cm air and 1 cm of its original Cauchy glass.
                    let camera_group_m = if lens {
                        0.99 + 0.01 * (1.48 + 0.018 / (wavelength * 1e-3).powi(2))
                    } else {
                        1.0
                    };
                    for query in queries.iter() {
                        let delta = query.position_m.delta_from(Point3::new(0.0, 0.0, 1.0));
                        let distance = delta.norm();
                        let expected = (camera_group_m + group_index * distance) / C;
                        assert_eq!(query.observation_time_s, clock);
                        assert_eq!((query.object_id, query.wavelength_nm), (9, wavelength));
                        near(query.group_delay_s.to_f64(), expected, 3e-7);
                        near(
                            (query.emission_time_s - epoch).to_f64(),
                            40e-9 - expected,
                            3e-7,
                        );
                    }
                }
            }
        }
    });
}

#[test]
fn g1_smith_finite_shutter_preserves_pulse_energy_and_refuses_unqualified_observations() {
    with_cx(|cx| {
        let air = glass(1.0, 0.1);
        let material = glass(1.5, 0.2);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&material),
        ];
        let geometry =
            InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, 3.0, 0.3)]).unwrap();
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![receiver(), opaque(9)],
            100_000,
        )
        .unwrap();
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_constant(LE).unwrap(),
            EmissionSides::Both,
        )];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [seed(Vec3::new(0.0, 0.0, -1.0))];
        let bindings = [binding(&light, &seeds, true)];
        let lights = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let physical = sequential::physical_camera();
        let optical = sequential::optical(&physical, cx, 550.0);
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let pulse = PiecewiseLinearHistory::try_new(
            cx,
            epoch,
            &[
                (-1e-6, 0.0),
                (0.0, 0.0),
                (1e-9, 2.0),
                (2e-9, 0.0),
                (1e-6, 0.0),
            ],
        )
        .unwrap();
        let timed =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&pulse)]).unwrap();
        let source = TransientCurvedAreaLightCameraScene::try_new(cx, &timed, &lights).unwrap();
        let limits = SmithWalkSettings::try_new(64).unwrap();
        let observer = source.with_smith_dielectric(cx, limits).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let window = ExposureWindow::try_new(epoch, epoch + Dd::from_f64(100e-9)).unwrap();
        for lens in [false, true] {
            let baseline = if lens {
                observer.source_area_camera().sample_sequential_ray(
                    cx,
                    &optical,
                    0,
                    settings(),
                    &mut workspace,
                    71,
                    5,
                )
            } else {
                observer.source_area_camera().sample_ray(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    settings(),
                    &mut workspace,
                    71,
                    5,
                )
            }
            .unwrap();
            let actual = if lens {
                observer.sample_sequential_ray_during_exposure(
                    cx,
                    &optical,
                    0,
                    window,
                    settings(),
                    &mut workspace,
                    71,
                    5,
                )
            } else {
                observer.sample_ray_during_exposure(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    window,
                    settings(),
                    &mut workspace,
                    71,
                    5,
                )
            }
            .unwrap();
            let steady = baseline.radiance_sample().radiance_w_m2_sr_nm();
            assert!(steady > 0.0);
            // The entire triangular pulse is inside every retarded interval:
            // area 2 ns, shutter 100 ns. A midpoint sample would return zero.
            near(
                actual.mean_radiance_sample().radiance_w_m2_sr_nm(),
                steady * 0.02,
                2e-12,
            );
            near(
                actual.integrated_radiance_j_m2_sr_nm(),
                steady * 2e-9,
                2e-12,
            );
            assert_eq!(
                actual.mean_radiance_sample().path(),
                baseline.radiance_sample().path()
            );
            assert_eq!(actual.mean_area_sample().images(), baseline.images());
            assert_eq!(
                actual.mean_area_sample().trace_evaluations(),
                baseline.trace_evaluations()
            );
            assert_eq!(
                actual.mean_area_sample().field_evaluations(),
                baseline.field_evaluations()
            );
            assert_eq!(
                actual.mean_area_sample().query_work(),
                baseline.query_work()
            );
        }
        let untimed = physical
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
            observer
                .sample_sequential_ray_at_time(
                    cx,
                    &untimed,
                    0,
                    epoch,
                    settings(),
                    &mut workspace,
                    71,
                    5
                )
                .unwrap_err(),
            TransientAreaLightCameraError::Source(TransientCameraError::MissingLensGroupDelay)
        );
        let analyzed = TransientCurvedCameraScene::try_new(cx, &camera, &[])
            .unwrap()
            .with_analyzer(CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap());
        let analyzed_source =
            TransientCurvedAreaLightCameraScene::try_new(cx, &analyzed, &lights).unwrap();
        let analyzed_observer = analyzed_source.with_smith_dielectric(cx, limits).unwrap();
        assert_eq!(
            analyzed_observer
                .sample_sequential_ray_at_time(
                    cx,
                    &optical,
                    0,
                    epoch,
                    settings(),
                    &mut workspace,
                    71,
                    5
                )
                .unwrap_err(),
            TransientAreaLightCameraError::Source(TransientCameraError::AnalyzerMustPrecedeLens)
        );
        let failure = History::new(|_| Err(SurfaceHistoryError::OutsideTimeDomain));
        let failed =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&failure)]).unwrap();
        let failed_source =
            TransientCurvedAreaLightCameraScene::try_new(cx, &failed, &lights).unwrap();
        let failed_observer = failed_source.with_smith_dielectric(cx, limits).unwrap();
        assert_eq!(
            failed_observer
                .sample_ray_at_time(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    epoch,
                    settings(),
                    &mut workspace,
                    71,
                    5
                )
                .unwrap_err(),
            TransientAreaLightCameraError::Source(TransientCameraError::History {
                object_id: 9,
                error: SurfaceHistoryError::OutsideTimeDomain,
            })
        );
        assert_eq!(failure.queries.lock().unwrap().len(), 1);
        assert_eq!(
            failed_observer
                .sample_ray_during_exposure(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    window,
                    settings(),
                    &mut workspace,
                    71,
                    5
                )
                .unwrap_err(),
            TransientAreaLightCameraError::Source(TransientCameraError::History {
                object_id: 9,
                error: SurfaceHistoryError::ExposureUnsupported,
            })
        );
    });
}
