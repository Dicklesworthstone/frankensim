//! Original finite-area connections, per-image clocks and finite shutters.

use super::*;

mod sequential;
mod smith_dielectric;
use crate::gradient_index::bundle::{IndexHessian, IndexJet, SecondOrderIndexField};
use crate::gradient_index::{GradientIndexError, GradientIndexField, GradientIndexSample};
use crate::polarization::dielectric::scene::curved::camera::transient::area_lights::{
    TransientAreaLightCameraError, TransientCurvedAreaLightCameraScene,
};
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::ExposureWindow;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::piecewise_linear::PiecewiseLinearHistory;
use crate::polarization::dielectric::scene::curved::camera::transient::{
    SurfaceEmissionHistory, SurfaceEmissionHistoryBinding, SurfaceEmissionQuery,
    SurfaceHistoryError, TransientCameraError, TransientCurvedCameraScene,
};
use fs_blake3::ContentHash;
use fs_math::dd::Dd;
use std::sync::Mutex;

const C: f64 = 299_792_458.0;
const AREA_DRAW: [f64; 2] = [4.0 / 9.0, 0.5];

struct History<F> {
    evaluate: F,
    queries: Mutex<Vec<SurfaceEmissionQuery>>,
}

impl<F> History<F> {
    fn new(evaluate: F) -> Self {
        Self {
            evaluate,
            queries: Mutex::new(Vec::new()),
        }
    }
}

impl<F: Fn(SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> + Sync> SurfaceEmissionHistory
    for History<F>
{
    fn sample(&self, cx: &Cx<'_>, query: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        self.queries.lock().unwrap().push(query);
        (self.evaluate)(query)
    }
}

fn history_binding(history: &dyn SurfaceEmissionHistory) -> SurfaceEmissionHistoryBinding<'_> {
    SurfaceEmissionHistoryBinding {
        object_id: 9,
        source_identity: ContentHash([73; 32]),
        history,
    }
}

#[test]
fn g3_area_histories_preserve_native_mis_work_and_use_the_complete_group_clock() {
    with_cx(|cx| {
        let geometry =
            InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, -1.0, 0.3)])
                .unwrap();
        let material = DielectricGlass::new(
            CauchyIor::try_new(1.2, 0.01, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(1, 0), opaque(9)],
            100_000,
        )
        .unwrap();
        let emitters = [emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
        let bindings = [AreaLightBinding {
            source: &light,
            receiver_object_id: 1,
            source_side: BoundarySide::Positive,
            initial_medium: 0,
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds,
            image_options: image_options(1),
            position_tolerance_m: 2e-7,
        }];
        let lights = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let history = History::new(|_| Ok(1.0));
        let timed =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&history)]).unwrap();
        let observer = TransientCurvedAreaLightCameraScene::try_new(cx, &timed, &lights).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let baseline = lights
            .sample_ray_with_samples(
                cx,
                camera_ray(),
                550.0,
                0,
                settings(),
                &mut workspace,
                |_| [0.0; 3],
                |_, _| AREA_DRAW,
            )
            .unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(20e-9);
        assert_eq!(clock.to_f64(), epoch.to_f64());
        let mut surface_draws = Vec::new();
        let mut area_draws = Vec::new();
        let result = observer
            .sample_ray_with_samples_at_time(
                cx,
                camera_ray(),
                550.0,
                0,
                clock,
                settings(),
                &mut workspace,
                |event| {
                    surface_draws.push(event);
                    [0.0; 3]
                },
                |event, binding| {
                    area_draws.push((event, binding));
                    AREA_DRAW
                },
            )
            .unwrap();
        assert_eq!(result.area_sample(), &baseline);
        assert_eq!(surface_draws, [0]);
        assert_eq!(area_draws, [(0, 0)]);
        let queries = history.queries.lock().unwrap();
        assert_eq!(
            queries.len(),
            2,
            "one NEE image and one complementary emitter hit"
        );
        let group_m = 3.0 * (1.2 + 3.0 * 0.01 / 0.55_f64.powi(2));
        for query in queries.iter() {
            assert_eq!((query.object_id, query.wavelength_nm), (9, 550.0));
            assert_eq!(query.observation_time_s, clock);
            near(query.group_delay_s.to_f64() * C, group_m, 3e-7);
            near(
                (query.emission_time_s - epoch).to_f64(),
                20e-9 - group_m / C,
                3e-6,
            );
            assert!(query.position_m.x.abs() < 2e-7 && query.position_m.y.abs() < 2e-7);
            assert!((query.position_m.z + 1.0).abs() < 2e-7);
        }
        assert_eq!(
            queries[0].position_m,
            light.sample_area(AREA_DRAW).unwrap().position_m()
        );
        let other_camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let foreign = TransientCurvedCameraScene::try_new(cx, &other_camera, &[]).unwrap();
        assert!(matches!(
            TransientCurvedAreaLightCameraScene::try_new(cx, &foreign, &lights),
            Err(TransientAreaLightCameraError::CameraIdentityMismatch)
        ));
    });
}

#[test]
fn g1_area_shutter_integrates_the_retarded_pulse_on_both_mis_strategies() {
    with_cx(|cx| {
        let geometry =
            InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, -1.0, 0.3)])
                .unwrap();
        let material = glass(1.3, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(1, 0), opaque(9)],
            100_000,
        )
        .unwrap();
        let emitters = [emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
        let bindings = [AreaLightBinding {
            source: &light,
            receiver_object_id: 1,
            source_side: BoundarySide::Positive,
            initial_medium: 0,
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds,
            image_options: image_options(1),
            position_tolerance_m: 2e-7,
        }];
        let lights = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
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
        let observer = TransientCurvedAreaLightCameraScene::try_new(cx, &timed, &lights).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let baseline = lights
            .sample_ray_with_samples(
                cx,
                camera_ray(),
                550.0,
                0,
                settings(),
                &mut workspace,
                |_| [0.0; 3],
                |_, _| AREA_DRAW,
            )
            .unwrap();
        let arrival = epoch + Dd::from_f64(3.9) / Dd::from_f64(C);
        let window = ExposureWindow::try_new(
            arrival - Dd::from_f64(0.5e-9),
            arrival + Dd::from_f64(2.5e-9),
        )
        .unwrap();
        let result = observer
            .sample_ray_with_samples_during_exposure(
                cx,
                camera_ray(),
                550.0,
                0,
                window,
                settings(),
                &mut workspace,
                |_| [0.0; 3],
                |_, _| AREA_DRAW,
            )
            .unwrap();
        near(
            result.mean_radiance_sample().radiance_w_m2_sr_nm(),
            baseline.radiance_sample().radiance_w_m2_sr_nm() * (2.0 / 3.0),
            2e-12,
        );
        near(
            result.mean_area_sample().area_radiance_w_m2_sr_nm(),
            baseline.area_radiance_w_m2_sr_nm() * (2.0 / 3.0),
            2e-12,
        );
        near(
            result.integrated_radiance_j_m2_sr_nm(),
            baseline.radiance_sample().radiance_w_m2_sr_nm() * 2e-9,
            2e-12,
        );
        assert_eq!(
            result.mean_area_sample().trace_evaluations(),
            baseline.trace_evaluations()
        );
        assert_eq!(
            result.mean_radiance_sample().path(),
            baseline.radiance_sample().path()
        );
    });
}

#[test]
fn g1_each_refracted_area_image_retains_its_own_source_clock() {
    with_cx(|cx| {
        let receiver = Point3::new(0.03, 0.0, 2.5);
        let geometry = InstanceScene::try_new(vec![
            panel(1, receiver.x, receiver.z, 1.0),
            GeometryInstance::try_new(
                2,
                hash_domain("transient-area-ball", b"unit"),
                SharedGeometry::chart(SphereChart {
                    center: Point3::new(0.0, 0.0, 0.0),
                    radius: 1.0,
                }),
                RigidTransform::identity(),
            )
            .unwrap(),
            panel(9, 0.0, -3.0, 0.3),
        ])
        .unwrap();
        let lens = glass(1.5, 0.2);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![
                diffuse(1, 1),
                SurfaceBinding {
                    object_id: 2,
                    interaction: Interaction::Dielectric {
                        outside: 1,
                        inside: 0,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let emitters = [emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let source = light.sample_area(AREA_DRAW).unwrap();
        let seeds = [
            seed(Vec3::new((-0.18_f64).sin(), 0.0, (-0.18_f64).cos())),
            seed(Vec3::new(0.14_f64.sin(), 0.0, 0.14_f64.cos())),
        ];
        let target = BundleConnectionTarget::try_diffuse_receiver(
            1,
            receiver,
            Vec3::new(0.0, 0.0, -1.0),
            2e-7,
        )
        .unwrap();
        let mut delays = [0.0; 2];
        let mut outgoing = Vec3::new(0.0, 0.0, -1.0);
        for index in 0..2 {
            let connection = scene
                .connect_bundle_from_surface_with_policy(
                    cx,
                    BundleConnectionLaunch {
                        ray: Ray {
                            origin: source.position_m(),
                            dir: seeds[index].initial_direction,
                        },
                        transverse_axis: Vec3::new(0.0, 1.0, 0.0),
                        wavelength_nm: 550.0,
                        normalized_source_stokes: Stokes::unpolarized(1.0),
                        initial_medium: 1,
                    },
                    BundleSurfaceSource::try_new(9, BoundarySide::Positive).unwrap(),
                    target,
                    settings().connection,
                    seeds[index].policy,
                )
                .unwrap();
            let propagation = connection.path().bundle().propagation().propagation();
            delays[index] = (1.0 + propagation.group_optical_path_length_m.unwrap()) / C;
            outgoing = propagation.ray.dir.scale(-1.0);
        }
        let separation = (delays[0] - delays[1]).abs();
        assert!(
            separation > 1e-13,
            "fixture must resolve two distinct optical clocks"
        );
        let draws = |event| {
            if event == 0 {
                [
                    outgoing.x * outgoing.x + outgoing.y * outgoing.y,
                    if outgoing.x >= 0.0 { 0.25 } else { 0.75 },
                    0.5,
                ]
            } else {
                [0.0, 0.0, 0.5]
            }
        };
        let ray = Ray {
            origin: Point3::new(receiver.x, 0.0, 1.5),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let binding = AreaLightBinding {
            source: &light,
            receiver_object_id: 1,
            source_side: BoundarySide::Positive,
            initial_medium: 1,
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds,
            image_options: image_options(2),
            position_tolerance_m: 2e-7,
        };
        let full_bindings = [binding];
        let first_bindings = [AreaLightBinding {
            seeds: &seeds[..1],
            image_options: image_options(1),
            ..binding
        }];
        let lights = CurvedAreaLightCameraScene::try_new(cx, &camera, &full_bindings).unwrap();
        let first = CurvedAreaLightCameraScene::try_new(cx, &camera, &first_bindings).unwrap();
        let history = History::new(|query: SurfaceEmissionQuery| {
            Ok(
                if (query.emission_time_s.to_f64() + delays[0]).abs() < separation * 0.25 {
                    1.0
                } else {
                    0.0
                },
            )
        });
        let timed =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&history)]).unwrap();
        let observer = TransientCurvedAreaLightCameraScene::try_new(cx, &timed, &lights).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(2, 2, 8).unwrap();
        let baseline = first
            .sample_ray_with_samples(
                cx,
                ray,
                550.0,
                1,
                settings(),
                &mut workspace,
                draws,
                |_, _| AREA_DRAW,
            )
            .unwrap();
        let result = observer
            .sample_ray_with_samples_at_time(
                cx,
                ray,
                550.0,
                1,
                Dd::ZERO,
                settings(),
                &mut workspace,
                draws,
                |_, _| AREA_DRAW,
            )
            .unwrap();
        assert_eq!(result.area_sample().images(), 2);
        near(
            result.area_sample().area_radiance_w_m2_sr_nm(),
            baseline.area_radiance_w_m2_sr_nm(),
            2e-6,
        );
        let queries = history.queries.lock().unwrap();
        assert_eq!(
            queries.len(),
            3,
            "each NEE image plus the original emitter hit"
        );
        for expected in delays {
            assert!(
                queries[..2]
                    .iter()
                    .any(|q| (q.group_delay_s.to_f64() - expected).abs() < 2e-15)
            );
        }
        assert!(
            queries[..2]
                .iter()
                .all(|q| q.position_m == source.position_m())
        );
    });
}

#[test]
fn g0_missing_connection_group_delay_or_source_failure_refuses_the_whole_observation() {
    struct Field {
        missing_group: bool,
    }
    impl GradientIndexField for Field {
        fn sample(
            &self,
            _: Point3,
            _: f64,
        ) -> Result<GradientIndexSample, GradientIndexError> {
            // A trace that starts with a group index and then loses it is refused
            // (GradientIndexError::GroupIndexUnavailable, gradient_index.rs rhs). A field
            // without one only below z = 0 therefore refused even the steady baseline: the
            // camera's first leg (z = 0 to the receiver at z = 1) keeps one, and its diffuse
            // continuation from the receiver runs back down into z < 0. The missing field has no
            // group index anywhere, so neither the camera leg nor the source connection carries
            // a group delay.
            GradientIndexSample::try_new(
                1.3,
                Vec3::new(0.0, 0.0, 0.0),
                0.0,
                if self.missing_group { None } else { Some(1.4) },
            )
        }
    }
    impl SecondOrderIndexField for Field {
        fn sample_second_order(
            &self,
            position: Point3,
            wavelength: f64,
        ) -> Result<IndexJet, GradientIndexError> {
            Ok(IndexJet::new(
                self.sample(position, wavelength)?,
                IndexHessian::ZERO,
            ))
        }
    }
    with_cx(|cx| {
        for missing_group in [false, true] {
            let field = Field { missing_group };
            let geometry =
                InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, -1.0, 0.3)])
                    .unwrap();
            let media = [SceneMedium::Differentiable {
                field: &field,
                identity: ContentHash([81; 32]),
            }];
            let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
            let scene = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &[],
                &[],
                &reflectances,
                vec![diffuse(1, 0), opaque(9)],
                100_000,
            )
            .unwrap();
            let emitters = [emitter(9)];
            let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
            let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
            let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
            let bindings = [AreaLightBinding {
                source: &light,
                receiver_object_id: 1,
                source_side: BoundarySide::Positive,
                initial_medium: 0,
                source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
                seeds: &seeds,
                image_options: image_options(1),
                position_tolerance_m: 2e-7,
            }];
            let lights = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
            let history = History::new(|_| Err(SurfaceHistoryError::OutsideTimeDomain));
            let timed =
                TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&history)])
                    .unwrap();
            let observer =
                TransientCurvedAreaLightCameraScene::try_new(cx, &timed, &lights).unwrap();
            let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
            let baseline = lights
                .sample_ray_with_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    settings(),
                    &mut workspace,
                    |_| [0.0; 3],
                    |_, _| AREA_DRAW,
                )
                .unwrap();
            assert!(baseline.area_radiance_w_m2_sr_nm() > 0.0);
            let result = observer.sample_ray_with_samples_at_time(
                cx,
                camera_ray(),
                550.0,
                0,
                Dd::ZERO,
                settings(),
                &mut workspace,
                |_| [0.0; 3],
                |_, _| AREA_DRAW,
            );
            let expected = if missing_group {
                TransientCameraError::MissingGroupDelay { object_id: 9 }
            } else {
                TransientCameraError::History {
                    object_id: 9,
                    error: SurfaceHistoryError::OutsideTimeDomain,
                }
            };
            assert_eq!(result, Err(TransientAreaLightCameraError::Source(expected)));
            assert_eq!(
                history.queries.lock().unwrap().len(),
                usize::from(!missing_group)
            );
        }
    });
}
