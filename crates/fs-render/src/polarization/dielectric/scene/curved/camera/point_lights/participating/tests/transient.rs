//! Retarded point histories preserve native collisions, survival and source images.

use super::*;
use crate::polarization::dielectric::scene::curved::camera::transient::{TransientCurvedCameraScene, SurfaceHistoryError};
use crate::polarization::dielectric::scene::curved::camera::transient::point_lights::*;
use crate::polarization::dielectric::scene::curved::camera::transient::point_lights::participating::*;
use crate::polarization::dielectric::scene::curved::camera::transient::point_lights::participating::images::TransientCurvedRayleighImageLightCameraScene;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::ExposureWindow;
use crate::polarization::dielectric::scene::curved::camera::point_lights::participating::images::{CurvedRayleighImageLightCameraScene, PointLightImageGroup, PointLightImageReceiver};
use crate::polarization::dielectric::scene::curved::camera::point_lights::PointLightSourcePolarization;
use crate::polarization::dielectric::scene::curved::bundle::connection::images::{BundleImageOptions, BundleImageSeed};
use crate::polarization::dielectric::scene::curved::participating::transmittance::images::CurvedTransmittanceImageWorkspace;
use crate::dielectric::DielectricSurface;
use crate::polarization::photometry::SpectralOptics;
use fs_geom::fixtures::SphereChart;
use fs_math::dd::Dd;
use std::sync::Mutex;

const C: f64 = 299_792_458.0;
const ID: ContentHash = ContentHash([126; 32]);
struct History<'a> {
    source: &'a PhotometricSource,
    epoch: Dd,
    offset: f64,
    slope: f64,
    calls: Mutex<Vec<(Vec3, Dd, Dd, Dd)>>,
}
impl PointEmissionHistory for History<'_> {
    fn sample(
        &self,
        cx: &Cx<'_>,
        query: PointEmissionQuery<'_>,
    ) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        assert!(core::ptr::eq(query.source, self.source));
        assert_eq!(query.source_identity, ID);
        self.calls.lock().unwrap().push((
            query.outgoing_direction,
            query.observation_time_s,
            query.emission_time_s,
            query.group_delay_s,
        ));
        Ok(self.offset + self.slope * (query.emission_time_s - self.epoch).to_f64())
    }
    fn average(
        &self,
        cx: &Cx<'_>,
        query: PointExposureQuery<'_>,
    ) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        assert!(core::ptr::eq(query.opening.source, self.source));
        assert_eq!(query.opening.source_identity, ID);
        Ok(self.offset
            + self.slope
                * (query.opening.emission_time_s + query.duration_s * Dd::from_f64(0.5)
                    - self.epoch)
                    .to_f64())
    }
}
fn history_binding<'a>(
    source: &'a PhotometricSource,
    history: &'a dyn PointEmissionHistory,
) -> PointEmissionHistoryBinding<'a> {
    PointEmissionHistoryBinding {
        source,
        source_identity: ID,
        history,
    }
}

#[test]
fn g1_molecular_point_histories_follow_real_collisions_and_independent_source_survival() {
    with_cx(|cx| {
        let geometry = geometry();
        let material = material();
        let field = AffineIndexField::try_new(
            &material,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            bounds(),
        )
        .unwrap();
        let media = [SceneMedium::Differentiable {
            field: &field,
            identity: FIELD_ID,
        }];
        let reflectance = [SpectralReflectance::constant(REFLECTANCE).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectance,
            surfaces(),
            100_000,
        )
        .unwrap();
        let emission = emitters();
        let camera = CurvedCameraScene::try_new(cx, &scene, &emission, None).unwrap();
        let fixtures = CurvedPointLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let lamp = lamp();
        let branches = [volume_branch(&lamp)];
        let law = LocalLaw { rate: 0.25 };
        let volumes = [Some(CurvedRayleighBinding::Field(&law))];
        let transport = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let lit = CurvedRayleighPointLightCameraScene::try_new_with_volume_branches(
            cx, &molecular, &fixtures, &branches,
        )
        .unwrap();
        let flight = |ordinal| match ordinal {
            0 => [1.0 - (-0.1_f64).exp(), 0.0, 0.2, 0.9],
            1 => [1.0 - (-0.15_f64).exp(), 0.0, 0.2, 0.1],
            _ => [0.99, 0.0, 0.2, 0.1],
        };
        let ordinary = molecular
            .sample_ray_with_tracking_samples(
                cx,
                ray(),
                550.0,
                0,
                camera_settings(),
                |_| [0.0; 3],
                flight,
            )
            .unwrap();
        let mut camera_draws = Vec::new();
        let mut source_draws = Vec::new();
        let sample = lit
            .sample_ray_with_tracking_samples(
                cx,
                ray(),
                550.0,
                0,
                settings(),
                |_| [0.0; 3],
                |ordinal| {
                    camera_draws.push(ordinal);
                    flight(ordinal)
                },
                |connection, ordinal| {
                    source_draws.push((connection, ordinal));
                    1.0 - (-(if ordinal == 0 { 0.2_f64 } else { 3.0_f64 })).exp()
                },
            )
            .unwrap();
        assert_eq!(
            sample.radiance_sample().path(),
            ordinary.radiance_sample().path()
        );
        assert_eq!(camera_draws, vec![0, 1, 2]);
        assert_eq!(source_draws, vec![(0, 0), (0, 1)]);
        assert_eq!(sample.volume_events(), 1);
        assert_eq!(
            sample.connections(),
            1,
            "the rejected first candidate contributes no source term"
        );
        assert_eq!(sample.source_candidates(), 1);
        let distance = 1.0625_f64.sqrt();
        let cosine_squared = 0.25_f64.powi(2) / distance.powi(2);
        let delta = 0.03;
        let phase_per_sr = 3.0 * ((1.0 + delta) + (1.0 - delta) * cosine_squared)
            / (8.0 * core::f64::consts::PI * (2.0 + delta));
        let intensity = lamp.flux_lm() * lamp.spectrum().bins_w_per_lumen_nm()[0]
            / (4.0 * core::f64::consts::PI);
        let expected = phase_per_sr * intensity / distance.powi(2)
            * (-ABSORPTION * (0.25 + distance)).exp()
            * 0.75;
        near(sample.fixture_radiance_w_m2_sr_nm(), expected);

        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(20e-9);
        let history = History {
            source: &lamp,
            epoch,
            offset: 2.0,
            slope: 0.0,
            calls: Mutex::new(Vec::new()),
        };
        let histories = [history_binding(&lamp, &history)];
        let timed = TransientCurvedCameraScene::try_new(cx, &camera, &[]).unwrap();
        let observer =
            TransientCurvedRayleighPointLightCameraScene::try_new(cx, &timed, &lit, &histories)
                .unwrap();
        let mut timed_camera_draws = Vec::new();
        let mut timed_source_draws = Vec::new();
        let result = observer
            .sample_ray_with_tracking_samples_at_time(
                cx,
                ray(),
                550.0,
                0,
                clock,
                settings(),
                |_| [0.0; 3],
                |ordinal| {
                    timed_camera_draws.push(ordinal);
                    flight(ordinal)
                },
                |connection, ordinal| {
                    timed_source_draws.push((connection, ordinal));
                    1.0 - (-(if ordinal == 0 { 0.2_f64 } else { 3.0_f64 })).exp()
                },
            )
            .unwrap();
        let actual = result.point_light_sample();
        near(actual.fixture_radiance_w_m2_sr_nm(), 2.0 * expected);
        near(
            actual.radiance_sample().radiance_w_m2_sr_nm(),
            sample.radiance_sample().radiance_w_m2_sr_nm() + sample.fixture_radiance_w_m2_sr_nm(),
        );
        assert_eq!(
            actual.radiance_sample().path(),
            sample.radiance_sample().path()
        );
        assert_eq!(actual.trace_evaluations(), sample.trace_evaluations());
        assert_eq!(actual.field_evaluations(), sample.field_evaluations());
        assert_eq!(actual.query_work(), sample.query_work());
        assert_eq!(actual.volume_events(), sample.volume_events());
        assert_eq!(actual.free_flight_draws(), sample.free_flight_draws());
        assert_eq!(
            actual.source_free_flight_draws(),
            sample.source_free_flight_draws()
        );
        assert_eq!(actual.source_candidates(), sample.source_candidates());
        assert_eq!(timed_camera_draws, camera_draws);
        assert_eq!(timed_source_draws, source_draws);
        let queries = history.calls.lock().unwrap();
        assert_eq!(
            queries.len(),
            1,
            "null candidates are not emitting vertices"
        );
        assert_eq!(queries[0].1, clock);
        assert_eq!(queries[0].2, clock - queries[0].3);
        near(queries[0].3.to_f64(), 1.2 * (0.25 + distance) / C);
        drop(queries);
        // A source ratio-tracking zero remains zero before history evaluation.
        let zero_law = LocalLaw { rate: 1.0 };
        let zero_volumes = [Some(CurvedRayleighBinding::Field(&zero_law))];
        let zero_transport =
            CurvedRayleighScene::try_new_with_fields(cx, &scene, &zero_volumes).unwrap();
        let zero_camera = CurvedRayleighCameraScene::try_new(cx, &camera, &zero_transport).unwrap();
        let zero_lights = CurvedRayleighPointLightCameraScene::try_new_with_volume_branches(
            cx,
            &zero_camera,
            &fixtures,
            &branches,
        )
        .unwrap();
        let zero_observer = TransientCurvedRayleighPointLightCameraScene::try_new(
            cx,
            &timed,
            &zero_lights,
            &histories,
        )
        .unwrap();
        let zero = zero_observer
            .sample_ray_with_tracking_samples_at_time(
                cx,
                ray(),
                550.0,
                0,
                clock,
                settings(),
                |_| [0.0; 3],
                |ordinal| {
                    [
                        if ordinal == 0 {
                            1.0 - (-0.25_f64).exp()
                        } else {
                            0.99
                        },
                        0.0,
                        0.2,
                        0.1,
                    ]
                },
                |_, ordinal| 1.0 - (-(if ordinal == 0 { 0.2_f64 } else { 3.0_f64 })).exp(),
            )
            .unwrap();
        assert_eq!(zero.point_light_sample().fixture_radiance_w_m2_sr_nm(), 0.0);
        assert!(zero.radiance_sample().radiance_w_m2_sr_nm() > 0.0);
        assert_eq!(zero.point_light_sample().source_candidates(), 1);
        assert_eq!(history.calls.lock().unwrap().len(), 1);
    });
}

#[test]
fn g1_polarized_ball_lens_images_retain_individual_source_clocks_and_finite_shutters() {
    with_cx(|cx| {
        // These three finite seeds resolve the three regular ball-lens images
        // independently checked by the native image-discovery regressions.
        let target = Point3::new(0.03, 0.0, 2.5);
        let geometry = InstanceScene::try_new(vec![
            GeometryInstance::try_new(
                1,
                hash_domain("image-point-camera-test", b"unit-sphere"),
                SharedGeometry::chart(SphereChart {
                    center: Point3::new(0.0, 0.0, 0.0),
                    radius: 1.0,
                }),
                RigidTransform::identity(),
            )
            .unwrap(),
            GeometryInstance::try_new(
                9,
                hash_domain("image-point-camera-test", b"diffuse-receiver"),
                SharedGeometry::mesh(TriMesh::new(
                    vec![
                        [target.x - 0.5, -0.5, target.z],
                        [target.x + 0.5, -0.5, target.z],
                        [target.x, 0.5, target.z],
                    ],
                    vec![[0, 1, 2]],
                )),
                RigidTransform::identity(),
            )
            .unwrap(),
        ])
        .unwrap();
        let glass = |index, absorption| {
            DielectricGlass::new(
                CauchyIor::try_constant(index).unwrap(),
                BeerLambertAbsorption::try_constant(absorption).unwrap(),
                GlassProvenance::Custom,
            )
        };
        let lens = glass(1.5, 0.2);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let reflectance = [SpectralReflectance::constant(REFLECTANCE).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectance,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        outside: 1,
                        inside: 0,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                SurfaceBinding {
                    object_id: 9,
                    interaction: Interaction::DepolarizingDiffuse {
                        medium: 1,
                        reflectance: 0,
                    },
                },
            ],
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let volumes = [
            Some(CurvedRayleighMedium::try_new(&spectrum, 0.3, &phase).unwrap()),
            None,
        ];
        let transport = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let source = PhotometricSource::try_new(
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -3.0]).unwrap(),
            100.0,
            1.0,
            AngularDistribution::isotropic(),
            SourceSpectrum::try_from_bins([1.0; SPECTRAL_BINS]).unwrap(),
        )
        .unwrap();
        let polarized = source
            .with_unpolarized_optics()
            .through(&SpectralOptics::linear_polarizer(0.3).unwrap())
            .unwrap();
        let source = polarized.intensity_source();
        let polarizations = [PointLightSourcePolarization::try_new(source, &polarized).unwrap()];
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(30e-9);
        let history = History {
            source,
            epoch,
            offset: 1.0,
            slope: 1e8,
            calls: Mutex::new(Vec::new()),
        };
        let histories = [history_binding(source, &history)];
        let timed = TransientCurvedCameraScene::try_new(cx, &camera, &[]).unwrap();
        let mut steady_sum = 0.0;
        let seed = |angle: f64| BundleImageSeed {
            initial_direction: Vec3::new(angle.sin(), 0.0, angle.cos()),
            policy: BundleBranchPolicy::TransmitOrReflectTir,
        };
        let seeds = [seed(-0.18), seed(0.065), seed(0.14), seed(-0.18)];
        let ray = Ray {
            origin: Point3::new(target.x, 0.0, 2.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let mut limits = settings();
        limits.point_lights.surfaces.max_distance_m = 8.0;
        limits.point_lights.surfaces.event.surface_tolerance_m = 1e-8;
        limits.point_lights.connection.scene.trace = limits.point_lights.surfaces;
        limits.point_lights.max_connections = 3;
        limits.point_lights.max_trace_evaluations = 512;
        limits.max_source_free_flight_draws = 0;
        limits.max_source_candidates = 0;

        // The diffuse continuation is oblique enough to miss the lens. Source
        // replay must not consume or alter these same camera variates.
        let singles = [seeds[0], seeds[1], seeds[2]].map(|seed| {
            let branches = [PointLightBranch {
                source,
                receiver_object_id: 9,
                initial_medium: 1,
                initial_direction: seed.initial_direction,
                source_transverse_axis: Vec3::new(1.0, 0.0, 0.0),
                policy: seed.policy,
                position_tolerance_m: 2e-7,
            }];
            let fixtures = CurvedPointLightCameraScene::try_new(cx, &camera, &branches)
                .unwrap()
                .with_source_polarizations(cx, &polarizations)
                .unwrap();
            let single =
                CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &fixtures).unwrap();
            let steady = single
                .sample_ray_with_tracking_samples(
                    cx,
                    ray,
                    550.0,
                    1,
                    limits,
                    |_| [0.8; 3],
                    |_| [0.99, 0.5, 0.2, 0.9],
                    |_, _| panic!("uniform source scattering must consume no tracking variates"),
                )
                .unwrap();

            steady_sum += steady.fixture_radiance_w_m2_sr_nm();
            let observer = TransientCurvedRayleighPointLightCameraScene::try_new(
                cx, &timed, &single, &histories,
            )
            .unwrap();
            let observed = observer
                .sample_ray_with_tracking_samples_at_time(
                    cx,
                    ray,
                    550.0,
                    1,
                    clock,
                    limits,
                    |_| [0.8; 3],
                    |_| [0.99, 0.5, 0.2, 0.9],
                    |_, _| panic!("uniform source survival consumes no random draws"),
                )
                .unwrap();
            let query = *history.calls.lock().unwrap().last().unwrap();
            near(
                observed.point_light_sample().fixture_radiance_w_m2_sr_nm(),
                steady.fixture_radiance_w_m2_sr_nm() * (1.0 + 1e8 * (query.2 - epoch).to_f64()),
            );
            assert_eq!(
                observed.radiance_sample().path(),
                steady.radiance_sample().path()
            );
            assert_eq!(
                observed.point_light_sample().trace_evaluations(),
                steady.trace_evaluations()
            );
            *observed.point_light_sample()
        });
        let single_queries = history.calls.lock().unwrap().clone();
        assert_eq!(single_queries.len(), 3);
        for i in 0..3 {
            assert_eq!(single_queries[i].1, clock);
            assert_eq!(single_queries[i].2, clock - single_queries[i].3);
            assert!(single_queries[i].3.to_f64() > 5.5 / C);
            for j in 0..i {
                assert!(
                    (single_queries[i].3 - single_queries[j].3).abs().to_f64() > 0.0,
                    "fixture needs distinct optical-image group delays"
                );
            }
        }
        history.calls.lock().unwrap().clear();
        for single in &singles {
            assert!(single.fixture_radiance_w_m2_sr_nm() > 0.0);
            assert_eq!(single.connections(), 1);
            assert_eq!(single.volume_events(), 0);
            assert_eq!(single.free_flight_draws(), 0);
            assert_eq!(
                single.radiance_sample().path(),
                singles[0].radiance_sample().path()
            );
        }
        let expected = singles
            .iter()
            .map(|sample| sample.fixture_radiance_w_m2_sr_nm())
            .sum::<f64>();
        let fixtures = CurvedPointLightCameraScene::try_new(cx, &camera, &[])
            .unwrap()
            .with_source_polarizations(cx, &polarizations)
            .unwrap();
        let owner =
            CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &fixtures).unwrap();
        let groups = [PointLightImageGroup {
            source,
            receiver: PointLightImageReceiver::DiffuseObject(9),
            initial_medium: 1,
            source_transverse_axis: Vec3::new(1.0, 0.0, 0.0),
            seeds: &seeds,
            image_options: BundleImageOptions {
                max_images: 3,
                min_separation_rad: 0.01,
                duplicate_tolerance_rad: 1e-5,
                max_trace_evaluations: 512,
            },
            position_tolerance_m: 2e-7,
        }];
        let images = CurvedRayleighImageLightCameraScene::try_new(cx, &owner, &groups).unwrap();
        let observer =
            TransientCurvedRayleighImageLightCameraScene::try_new(cx, &timed, &images, &histories)
                .unwrap();
        let mut workspace = CurvedTransmittanceImageWorkspace::try_new(3, 4).unwrap();
        let sample = |limits, workspace: &mut CurvedTransmittanceImageWorkspace| {
            observer
                .sample_ray_with_tracking_samples_at_time(
                    cx,
                    ray,
                    550.0,
                    1,
                    clock,
                    limits,
                    workspace,
                    |_| [0.8; 3],
                    |_| [0.99, 0.5, 0.2, 0.9],
                    |_, _| panic!("uniform source scattering must consume no tracking variates"),
                )
                .map(|observed| *observed.point_light_sample())
                .unwrap()
        };
        let combined = sample(limits, &mut workspace);
        let image_queries = history.calls.lock().unwrap().clone();
        assert_eq!(
            image_queries.len(),
            3,
            "duplicate seeds must not duplicate source histories"
        );
        for image in &image_queries {
            assert_eq!(image.1, clock);
            assert_eq!(image.2, clock - image.3);
            let matching = single_queries
                .iter()
                .min_by(|a, b| {
                    (a.3 - image.3)
                        .abs()
                        .to_f64()
                        .total_cmp(&(b.3 - image.3).abs().to_f64())
                })
                .unwrap();
            near(image.3.to_f64(), matching.3.to_f64());
            near(image.0.x, matching.0.x);
            near(image.0.z, matching.0.z);
        }
        assert_eq!(combined.connections(), 3);
        assert_eq!(combined.source_free_flight_draws(), 0);
        assert_eq!(combined.source_candidates(), 0);
        assert_eq!(
            combined.radiance_sample().path(),
            singles[0].radiance_sample().path()
        );
        near(combined.fixture_radiance_w_m2_sr_nm(), expected);
        near(combined.radiance_sample().radiance_w_m2_sr_nm(), expected);
        // Each reference includes one uniform survival replay. The duplicate
        // contributes its geometric solve, but neither a replay nor power.
        let expected_traces = singles.iter().map(|s| s.trace_evaluations()).sum::<u32>()
            + singles[0].trace_evaluations()
            - 1;
        assert_eq!(combined.trace_evaluations(), expected_traces);

        // Cached contractions require no further traces after the native
        // search and its three replays exactly exhaust the shared budget.
        limits.point_lights.max_trace_evaluations = combined.trace_evaluations();
        let exact = sample(limits, &mut workspace);
        assert_eq!(exact, combined);

        let window = ExposureWindow::try_new(clock, clock + Dd::from_f64(4e-9)).unwrap();
        let exposed = observer
            .sample_ray_with_tracking_samples_during_exposure(
                cx,
                ray,
                550.0,
                1,
                window,
                limits,
                &mut workspace,
                |_| [0.8; 3],
                |_| [0.99, 0.5, 0.2, 0.9],
                |_, _| panic!("source shutter adds no tracking draws"),
            )
            .unwrap();
        let mean = expected + steady_sum * 1e8 * 2e-9;
        near(exposed.mean_radiance_sample().radiance_w_m2_sr_nm(), mean);
        near(exposed.integrated_radiance_j_m2_sr_nm(), mean * 4e-9);
        assert_eq!(
            exposed.mean_point_light_sample().trace_evaluations(),
            combined.trace_evaluations()
        );
        assert_eq!(
            exposed.mean_radiance_sample().path(),
            combined.radiance_sample().path()
        );
    });
}

#[test]
fn g1_timed_full_stokes_lens_prefix_reaches_the_original_molecular_point_vertex() {
    use crate::camera::{
        Aperture, CameraProjection, LensSample, LensSurface, PhysicalCamera, SequentialLens,
    };
    use crate::polarization::dielectric::scene::curved::camera::CameraAnalyzer;
    with_cx(|cx| {
        let geometry = geometry();
        let material = DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(ABSORPTION).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectance = [SpectralReflectance::constant(REFLECTANCE).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectance,
            surfaces(),
            100_000,
        )
        .unwrap();
        let emission = emitters();
        let camera = CurvedCameraScene::try_new(cx, &native, &emission, None).unwrap();
        let fixtures = CurvedPointLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let lamp = lamp();
        let distance = 1.0225_f64.sqrt();
        let branches = [MolecularPointLightBranch {
            initial_direction: Vec3::new(-1.0, 0.0, -0.15).scale(1.0 / distance),
            ..volume_branch(&lamp)
        }];
        let spectrum = MassScatteringSpectrum::constant(0.4).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let transport = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let owner = CurvedRayleighPointLightCameraScene::try_new_with_volume_branches(
            cx, &molecular, &fixtures, &branches,
        )
        .unwrap();
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
                LensSurface::try_new(0.1, 0.0, 0.02, None).unwrap(),
            ],
        )
        .unwrap();
        let physical = PhysicalCamera::try_look_at(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
            1.0,
            Aperture::try_circular(0.01).unwrap(),
        )
        .unwrap()
        .with_sequential_lens(lens)
        .unwrap();
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let beams = physical
            .generate_sequential_spectral_timed_beams_analyzed_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[450.0, 650.0],
                analyzer.reference_axis(),
                analyzer.row(),
            )
            .unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(30e-9);
        let history = History {
            source: &lamp,
            epoch,
            offset: 2.0,
            slope: 0.0,
            calls: Mutex::new(Vec::new()),
        };
        let histories = [history_binding(&lamp, &history)];
        let timed = TransientCurvedCameraScene::try_new(cx, &camera, &[]).unwrap();
        let observer =
            TransientCurvedRayleighPointLightCameraScene::try_new(cx, &timed, &owner, &histories)
                .unwrap();
        let flight = |ordinal| {
            [
                if ordinal == 0 {
                    1.0 - (-0.1_f64).exp()
                } else {
                    0.99
                },
                0.0,
                0.2,
                0.0,
            ]
        };
        for optical in beams.into_iter().flatten() {
            assert!(optical.power().q != 0.0 || optical.power().u != 0.0);
            let baseline = owner
                .sample_sequential_ray_with_tracking_samples(
                    cx,
                    &optical,
                    0,
                    settings(),
                    |_| panic!("camera ends at opaque wall"),
                    flight,
                    |_, _| panic!("uniform source survival consumes no draws"),
                )
                .unwrap();
            let result = observer
                .sample_sequential_ray_with_tracking_samples_at_time(
                    cx,
                    &optical,
                    0,
                    clock,
                    settings(),
                    |_| panic!("camera ends at opaque wall"),
                    flight,
                    |_, _| panic!("uniform source survival consumes no draws"),
                )
                .unwrap();
            assert_eq!(baseline.volume_events(), 1);
            near(
                result.point_light_sample().fixture_radiance_w_m2_sr_nm(),
                2.0 * baseline.fixture_radiance_w_m2_sr_nm(),
            );
            near(
                result.radiance_sample().radiance_w_m2_sr_nm(),
                baseline.radiance_sample().radiance_w_m2_sr_nm()
                    + baseline.fixture_radiance_w_m2_sr_nm(),
            );
            assert_eq!(
                result.radiance_sample().path(),
                baseline.radiance_sample().path()
            );
            assert_eq!(
                result.point_light_sample().free_flight_draws(),
                baseline.free_flight_draws()
            );
            assert_eq!(
                result.point_light_sample().trace_evaluations(),
                baseline.trace_evaluations()
            );
            let query = *history.calls.lock().unwrap().last().unwrap();
            let expected_delay =
                optical.group_delay_s().unwrap() + Dd::from_f64(0.25 + distance) / Dd::from_f64(C);
            near(query.3.to_f64(), expected_delay.to_f64());
            assert_eq!(query.1, clock);
            assert_eq!(query.2, clock - query.3);
        }
    });
}
