use super::super::{
    CurvedPointLightCameraScene, PointLightBranch, PointLightCameraError, PointLightCameraSettings,
};
use super::{
    CurvedRayleighPointLightCameraScene, MolecularPointLightBranch,
    RayleighPointLightCameraSettings,
};
use crate::charts::{Ray, TriMesh};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::{AffineIndexField, GradientIndexOptions};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::photometry::source::{
    AngularDistribution, PhotometricSource, SPECTRAL_BINS, SourceSpectrum,
};
use crate::polarization::dielectric::scene::curved::bundle::CurvedBundleSettings;
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnectionOptions,
};
use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
use crate::polarization::dielectric::scene::curved::camera::{
    CurvedCameraScene, EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighBinding, CurvedRayleighField, CurvedRayleighMedium, CurvedRayleighSample,
    CurvedRayleighScene, CurvedRayleighTraceSettings,
};
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedTraceSettings, SceneMedium,
};
use crate::polarization::dielectric::scene::{Interaction, PathTermination, SurfaceBinding};
use crate::spectral::SpectralReflectance;
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Aabb, Point3, Vec3};

const FIELD_ID: ContentHash = ContentHash([179; 32]);
const REFLECTANCE: f64 = 0.4;
const ABSORPTION: f64 = 0.2;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 43,
                kernel_id: 8,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn material() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(1.2).unwrap(),
        BeerLambertAbsorption::try_constant(ABSORPTION).unwrap(),
        GlassProvenance::Custom,
    )
}

fn geometry() -> InstanceScene {
    InstanceScene::try_new(
        [(1, 1.0), (9, 0.0)]
            .map(|(id, z)| {
                GeometryInstance::try_new(
                    id,
                    hash_domain("participating-point-camera-test", &id.to_le_bytes()),
                    SharedGeometry::mesh(TriMesh::new(
                        vec![[-2.0, -2.0, z], [2.0, -2.0, z], [0.0, 2.0, z]],
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

fn surfaces() -> Vec<SurfaceBinding> {
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
    ]
}

fn emitters() -> [SurfaceEmitter; 2] {
    [(1, 2.0), (9, 3.0)].map(|(id, value)| {
        SurfaceEmitter::new(
            id,
            RadianceSpectrum::try_constant(value).unwrap(),
            EmissionSides::Both,
        )
    })
}

fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.5),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}

fn lamp() -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [1.0, 0.0, 0.5]).unwrap(),
        100.0,
        1.0,
        AngularDistribution::isotropic(),
        SourceSpectrum::try_from_bins([1.0; SPECTRAL_BINS]).unwrap(),
    )
    .unwrap()
}

fn branch(source: &PhotometricSource) -> PointLightBranch<'_> {
    PointLightBranch {
        source,
        receiver_object_id: 1,
        initial_medium: 0,
        initial_direction: Vec3::new(-1.0, 0.0, 0.5).scale(1.0 / 1.25_f64.sqrt()),
        source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
        policy: BundleBranchPolicy::Sequence(&[]),
        position_tolerance_m: 2e-6,
    }
}

fn volume_branch(source: &PhotometricSource) -> MolecularPointLightBranch<'_> {
    MolecularPointLightBranch {
        source,
        receiver_medium: 0,
        initial_medium: 0,
        initial_direction: Vec3::new(-1.0, 0.0, 0.25).scale(1.0 / 1.0625_f64.sqrt()),
        source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
        policy: BundleBranchPolicy::Sequence(&[]),
        position_tolerance_m: 2e-6,
    }
}

fn settings() -> RayleighPointLightCameraSettings {
    let surfaces = CurvedTraceSettings {
        max_surface_events: 8,
        max_distance_m: 6.0,
        max_total_steps: 16_384,
        max_query_work: 1_000_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-11,
                direction_tolerance: 1e-12,
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-7,
            ..GradientEventOptions::default()
        },
    };
    RayleighPointLightCameraSettings {
        point_lights: PointLightCameraSettings {
            surfaces,
            connection: BundleConnectionOptions {
                scene: CurvedBundleSettings {
                    trace: surfaces,
                    position_derivative_tolerance_m: 1e-11,
                    direction_derivative_tolerance: 1e-12,
                    ..CurvedBundleSettings::default()
                },
                ..BundleConnectionOptions::default()
            },
            max_connections: 8,
            max_trace_evaluations: 256,
        },
        max_volume_events: 8,
        max_source_free_flight_draws: 64,
        max_source_candidates: 64,
    }
}

fn camera_settings() -> CurvedRayleighTraceSettings {
    CurvedRayleighTraceSettings {
        surfaces: settings().point_lights.surfaces,
        max_volume_events: 8,
    }
}

fn bounds() -> Aabb {
    Aabb {
        min: Point3::new(-2.0, -2.0, -1.0),
        max: Point3::new(3.0, 2.0, 8.0),
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-6 * expected.abs().max(1e-10),
        "{actual:.17e} != {expected:.17e}"
    );
}

// Independent isotropic calibration, receiver cosine and both physical Beer
// legs. The camera prefix is conditional on its already sampled free flight.
fn clear_fixture_radiance(source: &PhotometricSource) -> f64 {
    let distance = 1.25_f64.sqrt();
    let intensity = source.flux_lm() * source.spectrum().bins_w_per_lumen_nm()[0]
        / (4.0 * core::f64::consts::PI);
    REFLECTANCE / core::f64::consts::PI * intensity * 0.5 / distance.powi(3)
        * (-ABSORPTION * (distance + 0.5)).exp()
}

struct LocalLaw {
    rate: f64,
}
impl CurvedRayleighField for LocalLaw {
    fn source_field_identity(&self) -> ContentHash {
        FIELD_ID
    }
    fn scattering_majorant_per_m(&self, _: f64) -> Result<f64, VolumeTransportError> {
        Ok(1.0)
    }
    fn sample_scattering(
        &self,
        _: Point3,
        _: f64,
    ) -> Result<CurvedRayleighSample, VolumeTransportError> {
        Ok(CurvedRayleighSample {
            scattering_per_m: self.rate,
            depolarization_ratio: 0.03,
        })
    }
}

#[test]
fn g1_fixture_survival_does_not_repeat_the_sampled_camera_prefix_survival() {
    with_cx(|cx| {
        let geometry = geometry();
        let material = material();
        let media = [SceneMedium::Homogeneous(&material)];
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
        let lamp = lamp();
        let branches = [branch(&lamp)];
        let fixtures = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        let clear = fixtures
            .sample_ray_with_samples(cx, ray(), 550.0, 0, settings().point_lights, |_| [0.0; 3])
            .unwrap();
        let phase = PhaseFunction::Rayleigh;
        for beta in [0.0, 0.7] {
            let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
            let volumes = [Some(
                CurvedRayleighMedium::try_new(&spectrum, beta, &phase).unwrap(),
            )];
            let transport = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
            let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
            let lit =
                CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &fixtures).unwrap();
            let ordinary = molecular
                .sample_ray_with_tracking_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    camera_settings(),
                    |_| [0.0; 3],
                    |_| [0.99, 0.5, 0.2, 0.9],
                )
                .unwrap();
            let mut limits = settings();
            limits.max_source_free_flight_draws = 0;
            limits.max_source_candidates = 0;
            let sample = lit
                .sample_ray_with_tracking_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    limits,
                    |_| [0.0; 3],
                    |_| [0.99, 0.5, 0.2, 0.9],
                    |_, _| panic!("uniform source scattering consumes no shadow RNG"),
                )
                .unwrap();
            assert_eq!(
                sample.radiance_sample().path(),
                ordinary.radiance_sample().path()
            );
            assert_eq!(sample.volume_events(), 0);
            assert_eq!(sample.free_flight_draws(), if beta == 0.0 { 0 } else { 2 });
            assert_eq!(sample.source_free_flight_draws(), 0);
            assert_eq!(sample.source_candidates(), 0);
            assert_eq!(sample.connections(), 1);
            let expected = clear_fixture_radiance(&lamp) * (-beta * 1.25_f64.sqrt()).exp();
            near(sample.fixture_radiance_w_m2_sr_nm(), expected);
            near(
                ordinary.radiance_sample().radiance_w_m2_sr_nm(),
                2.0 * (-0.1_f64).exp() + 3.0 * REFLECTANCE * (-0.3_f64).exp(),
            );
            near(
                sample.radiance_sample().radiance_w_m2_sr_nm(),
                ordinary.radiance_sample().radiance_w_m2_sr_nm() + expected,
            );
            if beta == 0.0 {
                assert_eq!(sample.radiance_sample(), clear.radiance_sample());
                assert_eq!(
                    sample.fixture_radiance_w_m2_sr_nm(),
                    clear.fixture_radiance_w_m2_sr_nm()
                );
                assert_eq!(sample.trace_evaluations(), clear.trace_evaluations());
            }
        }
    });
}

#[test]
fn g1_stochastic_shadow_zero_keeps_positive_emission_and_the_original_camera_draws() {
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
        let lamp = lamp();
        let branches = [branch(&lamp)];
        let fixtures = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        let law = LocalLaw { rate: 1.0 };
        let volumes = [Some(CurvedRayleighBinding::Field(&law))];
        let transport = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let lit = CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &fixtures).unwrap();
        let ordinary = molecular
            .sample_ray_with_tracking_samples(
                cx,
                ray(),
                550.0,
                0,
                camera_settings(),
                |_| [0.0; 3],
                |_| [0.99, 0.5, 0.2, 0.9],
            )
            .unwrap();
        for first_flight in [0.3_f64, 3.0] {
            let mut surface_draws = Vec::new();
            let mut camera_draws = Vec::new();
            let sample = lit
                .sample_ray_with_tracking_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    settings(),
                    |ordinal| {
                        surface_draws.push(ordinal);
                        [0.0; 3]
                    },
                    |ordinal| {
                        camera_draws.push(ordinal);
                        [0.99, 0.5, 0.2, 0.9]
                    },
                    |connection, ordinal| {
                        assert_eq!(connection, 0);
                        1.0 - (-(if ordinal == 0 { first_flight } else { 3.0 })).exp()
                    },
                )
                .unwrap();
            assert_eq!(surface_draws, vec![0]);
            assert_eq!(camera_draws, vec![0, 1]);
            assert_eq!(
                sample.radiance_sample().path(),
                ordinary.radiance_sample().path()
            );
            assert_eq!(
                sample.radiance_sample().path().termination(),
                PathTermination::Opaque
            );
            assert_eq!(sample.connections(), 1);
            assert_eq!(sample.volume_events(), 0);
            if first_flight < 1.25_f64.sqrt() {
                assert!(sample.source_candidates() >= 1);
                assert_eq!(sample.fixture_radiance_w_m2_sr_nm(), 0.0);
                assert_eq!(sample.radiance_sample(), ordinary.radiance_sample());
                assert!(sample.radiance_sample().radiance_w_m2_sr_nm() > 0.0);
            } else {
                assert_eq!(sample.source_candidates(), 0);
                near(
                    sample.fixture_radiance_w_m2_sr_nm(),
                    clear_fixture_radiance(&lamp),
                );
            }
        }
    });
}

#[test]
fn g0_source_stream_ordinals_and_tracking_limits_are_cumulative_across_fixtures() {
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
        let lamps = [lamp(), lamp()];
        let branches = [branch(&lamps[0]), branch(&lamps[1])];
        let fixtures = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        let law = LocalLaw { rate: 0.5 };
        let volumes = [Some(CurvedRayleighBinding::Field(&law))];
        let transport = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let lit = CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &fixtures).unwrap();
        let source_sample = |ordinal| 1.0 - (-(if ordinal == 0 { 0.3_f64 } else { 3.0 })).exp();
        let mut source_draws = Vec::new();
        let complete = lit
            .sample_ray_with_tracking_samples(
                cx,
                ray(),
                550.0,
                0,
                settings(),
                |_| [0.0; 3],
                |_| [0.99, 0.5, 0.2, 0.9],
                |connection, ordinal| {
                    source_draws.push((connection, ordinal));
                    source_sample(ordinal)
                },
            )
            .unwrap();
        assert_eq!(source_draws, vec![(0, 0), (0, 1), (1, 0), (1, 1)]);
        assert_eq!(complete.connections(), 2);
        assert_eq!(complete.source_free_flight_draws(), 4);
        assert_eq!(complete.source_candidates(), 2);
        // Each authored source has one factor 1-beta/M=.5. The independent
        // calibrated sum is therefore exactly one clear fixture contribution.
        near(
            complete.fixture_radiance_w_m2_sr_nm(),
            clear_fixture_radiance(&lamps[0]),
        );
        for (draw_limit, candidate_limit, expected) in [
            (3, 64, PointLightCameraError::SourceDrawLimit),
            (64, 1, PointLightCameraError::SourceCandidateLimit),
        ] {
            let mut limited = settings();
            limited.max_source_free_flight_draws = draw_limit;
            limited.max_source_candidates = candidate_limit;
            let mut draws = Vec::new();
            assert_eq!(
                lit.sample_ray_with_tracking_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    limited,
                    |_| panic!("failure during fixture observation precedes diffuse continuation"),
                    |_| [0.99, 0.5, 0.2, 0.9],
                    |connection, ordinal| {
                        draws.push((connection, ordinal));
                        source_sample(ordinal)
                    },
                )
                .unwrap_err(),
                expected
            );
            assert!(
                draws.contains(&(1, 0)),
                "the first complete fixture alone must fit the shared limit"
            );
        }
    });
}

#[test]
fn g1_point_fixture_at_molecular_collision_uses_transverse_area_and_one_phase_factor() {
    with_cx(|cx| {
        let geometry = geometry();
        let material = material();
        let media = [SceneMedium::Homogeneous(&material)];
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
        let phase = PhaseFunction::Rayleigh;
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        for beta in [0.2_f64, 0.9] {
            let volumes = [Some(
                CurvedRayleighMedium::try_new(&spectrum, beta, &phase).unwrap(),
            )];
            let transport = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
            let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
            let lit = CurvedRayleighPointLightCameraScene::try_new_with_volume_branches(
                cx, &molecular, &fixtures, &branches,
            )
            .unwrap();
            let flight = |ordinal| {
                [
                    if ordinal == 0 {
                        1.0 - (-beta * 0.25).exp()
                    } else {
                        0.99
                    },
                    0.0,
                    0.2,
                ]
            };
            let ordinary = molecular
                .sample_ray_with_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    camera_settings(),
                    |_| panic!("backscattered camera path hits the original opaque wall"),
                    flight,
                )
                .unwrap();
            let mut limits = settings();
            limits.max_source_free_flight_draws = 0;
            limits.max_source_candidates = 0;
            let sample = lit
                .sample_ray_with_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    limits,
                    |_| panic!("volume NEE does not create a surface interaction"),
                    flight,
                    |_, _| panic!("uniform source survival is analytic"),
                )
                .unwrap();
            assert_eq!(
                sample.radiance_sample().path(),
                ordinary.radiance_sample().path()
            );
            assert_eq!(sample.volume_events(), 1);
            assert_eq!(sample.free_flight_draws(), 2);
            assert_eq!(sample.connections(), 1);
            assert_eq!(sample.source_free_flight_draws(), 0);
            assert_eq!(sample.source_candidates(), 0);
            let distance = 1.0625_f64.sqrt();
            let cosine_squared = 0.25_f64.powi(2) / distance.powi(2);
            let phase_per_sr = 3.0 * (1.0 + cosine_squared) / (16.0 * core::f64::consts::PI);
            let intensity = lamp.flux_lm() * lamp.spectrum().bins_w_per_lumen_nm()[0]
                / (4.0 * core::f64::consts::PI);
            // Camera survival and beta are in the sampled event probability.
            // Only true camera absorption and source (absorption+scattering)
            // enter this CONDITIONAL radiance. No Lambertian pi or surface
            // cosine belongs at the molecular receiving point.
            let expected = phase_per_sr * intensity / distance.powi(2)
                * (-ABSORPTION * 0.25 - (ABSORPTION + beta) * distance).exp();
            near(sample.fixture_radiance_w_m2_sr_nm(), expected);
            near(
                sample.radiance_sample().radiance_w_m2_sr_nm(),
                ordinary.radiance_sample().radiance_w_m2_sr_nm() + expected,
            );
        }
    });
}

#[test]
fn g1_null_candidates_do_not_emit_and_local_phase_is_shared_with_camera_continuation() {
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
        let duplicate = [branches[0], branches[0]];
        assert!(matches!(
            CurvedRayleighPointLightCameraScene::try_new_with_volume_branches(
                cx, &molecular, &fixtures, &duplicate,
            ),
            Err(PointLightCameraError::DuplicateBinding {
                first: 0,
                second: 1
            })
        ));
    });
}
