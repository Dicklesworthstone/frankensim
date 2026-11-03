use super::super::{
    CurvedRayleighBinding, CurvedRayleighField, CurvedRayleighMedium, CurvedRayleighSample,
    CurvedRayleighScene,
};
use super::{CurvedTransmittanceError, CurvedTransmittanceOptions};
use crate::charts::{Ray, TriMesh};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricEvent, DielectricGlass, DielectricSurface,
    GlassProvenance,
};
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::{AffineIndexField, GradientIndexOptions};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::bundle::CurvedBundleSettings;
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnectionError, BundleConnectionLaunch, BundleConnectionOptions,
    BundleConnectionTarget,
};
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedTraceSettings, SceneMedium,
};
use crate::polarization::dielectric::scene::{Interaction, PathTermination, SurfaceBinding};
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Aabb, Point3, Vec3};
use std::sync::Mutex;

const FIELD_ID: ContentHash = ContentHash([178; 32]);

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
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
        );
        f(&cx, &gate);
    });
}

fn glass(index: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn plane(id: u64, center: Point3) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("curved-source-survival-test", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![
                [center.x - 3.0, -3.0, center.z],
                [center.x + 3.0, -3.0, center.z],
                [center.x, 3.0, center.z],
            ],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn detector() -> SurfaceBinding {
    SurfaceBinding {
        object_id: 9,
        interaction: Interaction::Detector {
            reference_axis: Vec3::new(0.0, 1.0, 0.0),
        },
    }
}

fn launch(direction: Vec3) -> BundleConnectionLaunch {
    BundleConnectionLaunch {
        ray: Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: direction,
        },
        transverse_axis: Vec3::new(0.0, 1.0, 0.0),
        wavelength_nm: 550.0,
        normalized_source_stokes: Stokes::unpolarized(1.0),
        initial_medium: 0,
    }
}

fn target(point: Point3) -> BundleConnectionTarget {
    BundleConnectionTarget::try_new(9, point, Vec3::new(0.0, 0.0, -1.0), 2e-7).unwrap()
}

fn options() -> BundleConnectionOptions {
    BundleConnectionOptions {
        scene: CurvedBundleSettings {
            trace: CurvedTraceSettings {
                max_surface_events: 8,
                max_distance_m: 8.0,
                max_total_steps: 16_384,
                max_query_work: 1_000_000,
                event: GradientEventOptions {
                    integration: GradientIndexOptions {
                        position_tolerance_m: 1e-11,
                        direction_tolerance: 1e-12,
                        optical_path_tolerance_m: 1e-11,
                        optical_depth_tolerance: 1e-12,
                        relative_tolerance: 1e-11,
                        ..GradientIndexOptions::default()
                    },
                    surface_tolerance_m: 1e-8,
                    ..GradientEventOptions::default()
                },
            },
            position_derivative_tolerance_m: 1e-10,
            direction_derivative_tolerance: 1e-11,
            ..CurvedBundleSettings::default()
        },
        ..BundleConnectionOptions::default()
    }
}

fn bounds() -> Aabb {
    Aabb {
        min: Point3::new(-0.4, -3.0, -1.0),
        max: Point3::new(4.0, 3.0, 10.0),
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn flight_sample(majorant: f64, distance: f64) -> f64 {
    1.0 - (-majorant * distance).exp()
}

#[test]
fn g1_volume_receiver_survival_ends_at_the_native_point_without_a_surface_event() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, Point3::new(0.0, 0.0, 4.0))]).unwrap();
        let material = glass(1.2, 0.17);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        let spectrum = MassScatteringSpectrum::constant(0.4).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
        let receiver = BundleConnectionTarget::try_volume_receiver(
            0,
            Point3::new(0.0, 0.0, 1.3),
            Vec3::new(0.0, 0.0, -1.0),
            2e-7,
        )
        .unwrap();
        let received = scattering
            .connect_bundle_with_transmittance_samples(
                cx,
                launch(Vec3::new(0.0, 0.0, 1.0)),
                receiver,
                options(),
                BundleBranchPolicy::Sequence(&[]),
                CurvedTransmittanceOptions {
                    max_free_flight_draws: 0,
                    max_tracking_candidates: 0,
                },
                |_| panic!("uniform source survival must not draw random numbers"),
            )
            .unwrap();
        let path = received.connection().path();
        assert!(path.is_volume_receiver());
        assert!(!path.is_diffuse_receiver());
        assert_eq!(path.medium_index(), 0);
        assert_eq!(path.surface_events(), 0);
        assert!(path.last_contact().is_none());
        assert_eq!(received.free_flight_draws(), 0);
        assert_eq!(received.tracking_candidates(), 0);
        assert_eq!(
            received.trace_evaluations(),
            received.connection().trace_evaluations() + 1
        );
        let propagation = path.bundle().propagation().propagation();
        near(propagation.distance_m, 1.3, 2e-8);
        near(propagation.optical_depth, 0.17 * 1.3, 2e-8);
        near(received.log_scattering_survival(), -0.4 * 1.3, 2e-8);
    });
}

#[test]
fn g1_zero_ratio_sample_still_validates_the_curved_volume_receiver() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, Point3::new(0.0, 0.0, 4.0))]).unwrap();
        let material = glass(1.2, 0.15);
        let field = AffineIndexField::try_new(
            &material,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.4, 0.0, 0.0),
            bounds(),
        )
        .unwrap();
        let media = [SceneMedium::Differentiable {
            field: &field,
            identity: FIELD_ID,
        }];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        let law = LocalLaw::constant(1.0, 1.0);
        let volumes = [Some(CurvedRayleighBinding::Field(&law))];
        let scattering = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
        let receiver = BundleConnectionTarget::try_volume_receiver(
            0,
            affine_point(1.2),
            Vec3::new(0.0, 0.0, -1.0),
            2e-7,
        )
        .unwrap();
        let received = scattering
            .connect_bundle_with_transmittance_samples(
                cx,
                launch(Vec3::new(0.0, 0.0, 1.0)),
                receiver,
                options(),
                BundleBranchPolicy::Sequence(&[]),
                CurvedTransmittanceOptions::default(),
                |ordinal| {
                    assert_eq!(ordinal, 0, "exact zero stops stochastic work, not geometry");
                    flight_sample(1.0, 0.25)
                },
            )
            .unwrap();
        assert!(received.ratio_tracking_zero());
        assert_eq!(received.log_scattering_survival(), f64::NEG_INFINITY);
        assert_eq!(received.tracking_candidates(), 1);
        assert_eq!(received.free_flight_draws(), 1);
        assert!(received.connection().path().is_volume_receiver());
        assert_eq!(received.connection().path().surface_events(), 0);
        assert!(received.connection().residual_m() <= receiver.position_tolerance_m());
        near(
            received
                .connection()
                .path()
                .bundle()
                .propagation()
                .propagation()
                .distance_m,
            1.2,
            2e-7,
        );
    });
}

// Closed-form ray for n(x)=n0+g*x launched along +z. The parameter is physical
// arc length; neither the chord nor optical/group path length gives this point.
fn affine_point(arc: f64) -> Point3 {
    let (n0, gradient) = (1.2_f64, 0.4_f64);
    Point3::new(
        ((n0 * n0 + gradient * gradient * arc * arc).sqrt() - n0) / gradient,
        0.0,
        n0 / gradient * (gradient * arc / n0).asinh(),
    )
}

struct LocalLaw {
    majorant: f64,
    low: f64,
    high: f64,
    queries: Mutex<Vec<Point3>>,
}

impl LocalLaw {
    fn constant(majorant: f64, local: f64) -> Self {
        Self {
            majorant,
            low: local,
            high: local,
            queries: Mutex::new(Vec::new()),
        }
    }
}

impl CurvedRayleighField for LocalLaw {
    fn source_field_identity(&self) -> ContentHash {
        FIELD_ID
    }

    fn scattering_majorant_per_m(&self, wavelength_nm: f64) -> Result<f64, VolumeTransportError> {
        assert_eq!(wavelength_nm, 550.0);
        Ok(self.majorant)
    }

    fn sample_scattering(
        &self,
        point: Point3,
        wavelength_nm: f64,
    ) -> Result<CurvedRayleighSample, VolumeTransportError> {
        assert_eq!(wavelength_nm, 550.0);
        self.queries.lock().unwrap().push(point);
        Ok(CurvedRayleighSample {
            scattering_per_m: if point.x < 0.08 { self.low } else { self.high },
            depolarization_ratio: 0.03,
        })
    }
}

#[test]
fn g1_uniform_scattering_survival_is_separate_from_native_beer_and_clear_connection() {
    with_cx(|cx, _| {
        let point = Point3::new(0.0, 0.0, 2.0);
        let geometry = InstanceScene::try_new(vec![plane(9, point)]).unwrap();
        let material = glass(1.3, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        let initial = launch(Vec3::new(0.0, 0.0, 1.0));
        let receiver = target(point);
        let policy = BundleBranchPolicy::Sequence(&[]);
        let original = scene
            .connect_bundle_with_policy(cx, initial, receiver, options(), policy)
            .unwrap();
        let spectrum = MassScatteringSpectrum::constant(0.3).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let no_work = CurvedTransmittanceOptions {
            max_free_flight_draws: 0,
            max_tracking_candidates: 0,
        };
        for density in [0.0, 2.0] {
            let volumes = [Some(
                CurvedRayleighMedium::try_new(&spectrum, density, &phase).unwrap(),
            )];
            let scattering = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
            let received = scattering
                .connect_bundle_with_transmittance_samples(
                    cx,
                    initial,
                    receiver,
                    options(),
                    policy,
                    no_work,
                    |_| panic!("uniform source survival is analytic and consumes no RNG"),
                )
                .unwrap();
            assert_eq!(*received.connection(), original);
            assert_eq!(received.free_flight_draws(), 0);
            assert_eq!(received.tracking_candidates(), 0);
            assert!(!received.ratio_tracking_zero());
            near(
                received.log_scattering_survival(),
                -0.3 * density * 2.0,
                1e-8,
            );
            let bundle = received.connection().path().bundle();
            near(bundle.propagation().propagation().optical_depth, 0.4, 1e-8);
            let incident = bundle
                .point_source_irradiance(receiver.normal(), options().projection, 2.0)
                .unwrap();
            near(incident, 0.5 * (-0.4_f64).exp(), 1e-8);
            near(
                incident * received.log_scattering_survival().exp(),
                0.5 * (-0.4 - 0.6 * density).exp(),
                1e-8,
            );
            assert_eq!(
                received.trace_evaluations(),
                original.trace_evaluations() + u32::from(density > 0.0)
            );
        }
    });
}

#[test]
fn g1_ratio_candidates_follow_the_solved_curved_arc_and_keep_the_original_bundle() {
    with_cx(|cx, _| {
        let point = affine_point(2.0);
        let geometry = InstanceScene::try_new(vec![plane(9, point)]).unwrap();
        let material = glass(1.2, 0.3);
        let field = AffineIndexField::try_new(
            &material,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.4, 0.0, 0.0),
            bounds(),
        )
        .unwrap();
        let media = [SceneMedium::Differentiable {
            field: &field,
            identity: FIELD_ID,
        }];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        let law = LocalLaw {
            majorant: 2.0,
            low: 0.5,
            high: 1.0,
            queries: Mutex::new(Vec::new()),
        };
        let volumes = [Some(CurvedRayleighBinding::Field(&law))];
        let scattering = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
        let direction = Vec3::new(0.03, 0.0, 1.0);
        let initial = launch(direction.scale(1.0 / direction.dot(direction).sqrt()));
        let receiver = target(point);
        let policy = BundleBranchPolicy::Sequence(&[]);
        let original = scene
            .connect_bundle_with_policy(cx, initial, receiver, options(), policy)
            .unwrap();
        assert!(original.iterations() > 0);
        let received = scattering
            .connect_bundle_with_transmittance_samples(
                cx,
                initial,
                receiver,
                options(),
                policy,
                CurvedTransmittanceOptions::default(),
                |ordinal| flight_sample(2.0, [0.4, 0.9, 3.0][ordinal as usize]),
            )
            .unwrap();
        assert_eq!(*received.connection(), original);
        assert_eq!(received.free_flight_draws(), 3);
        assert_eq!(received.tracking_candidates(), 2);
        assert_eq!(
            received.trace_evaluations(),
            original.trace_evaluations() + 1
        );
        assert!(!received.ratio_tracking_zero());
        let queries = law.queries.lock().unwrap();
        assert_eq!(queries.len(), 2, "Newton trials must not query scattering");
        for (actual, arc) in queries.iter().zip([0.4, 1.3]) {
            let expected = affine_point(arc);
            near(actual.x, expected.x, 5e-7);
            near(actual.y, expected.y, 5e-7);
            near(actual.z, expected.z, 5e-7);
        }
        // Ratio tracking with zero control component (Novak et al., 2014,
        // https://doi.org/10.1145/2661229.2661292): factors at the two actual
        // Poisson candidates are .75 and .5, with NO extra majorant survival.
        // The first curved point is below x=.08 while the endpoint chord is not.
        near(received.log_scattering_survival().exp(), 0.375, 1e-12);
        near(
            received
                .connection()
                .path()
                .bundle()
                .propagation()
                .propagation()
                .optical_depth,
            0.6,
            5e-7,
        );
    });
}

#[test]
fn g1_a_full_majorant_candidate_is_a_legitimate_zero_sample_not_a_failed_connection() {
    with_cx(|cx, _| {
        let point = Point3::new(0.0, 0.0, 2.0);
        let geometry = InstanceScene::try_new(vec![plane(9, point)]).unwrap();
        let material = glass(1.2, 0.15);
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
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        let law = LocalLaw::constant(1.0, 1.0);
        let volumes = [Some(CurvedRayleighBinding::Field(&law))];
        let scattering = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
        for first_flight in [0.3, 3.0] {
            let received = scattering
                .connect_bundle_with_transmittance_samples(
                    cx,
                    launch(Vec3::new(0.0, 0.0, 1.0)),
                    target(point),
                    options(),
                    BundleBranchPolicy::Sequence(&[]),
                    CurvedTransmittanceOptions::default(),
                    |ordinal| flight_sample(1.0, if ordinal == 0 { first_flight } else { 3.0 }),
                )
                .unwrap();
            let zero = first_flight < 2.0;
            assert_eq!(received.ratio_tracking_zero(), zero);
            assert_eq!(
                received.log_scattering_survival(),
                if zero { f64::NEG_INFINITY } else { 0.0 }
            );
            assert_eq!(
                received.connection().path().termination(),
                PathTermination::Detected
            );
            assert_eq!(received.connection().path().surface_events(), 1);
            near(
                received
                    .connection()
                    .path()
                    .bundle()
                    .propagation()
                    .propagation()
                    .distance_m,
                2.0,
                2e-8,
            );
            assert!(
                received
                    .connection()
                    .path()
                    .bundle()
                    .point_source_irradiance(target(point).normal(), options().projection, 1.0)
                    .unwrap()
                    > 0.0
            );
        }
    });
}

#[test]
fn g1_transmitted_source_branch_integrates_each_original_medium_without_changing_fresnel_or_area() {
    with_cx(|cx, _| {
        let point = Point3::new(0.0, 0.0, 3.0);
        let geometry =
            InstanceScene::try_new(vec![plane(1, Point3::new(0.0, 0.0, 1.0)), plane(9, point)])
                .unwrap();
        let first = glass(1.2, 0.1);
        let second = glass(1.5, 0.3);
        let media = [
            SceneMedium::Homogeneous(&first),
            SceneMedium::Homogeneous(&second),
        ];
        let surfaces = vec![
            SurfaceBinding {
                object_id: 1,
                interaction: Interaction::Dielectric {
                    outside: 1,
                    inside: 0,
                    surface: DielectricSurface::SMOOTH,
                },
            },
            detector(),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, surfaces, 100_000).unwrap();
        let a = MassScatteringSpectrum::constant(0.2).unwrap();
        let b = MassScatteringSpectrum::constant(0.7).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [
            Some(CurvedRayleighMedium::try_new(&a, 2.0, &phase).unwrap()),
            Some(CurvedRayleighMedium::try_new(&b, 1.0, &phase).unwrap()),
        ];
        let scattering = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
        let initial = launch(Vec3::new(0.0, 0.0, 1.0));
        let receiver = target(point);
        let plan = [DielectricEvent::Transmission];
        let policy = BundleBranchPolicy::Sequence(&plan);
        let original = scene
            .connect_bundle_with_policy(cx, initial, receiver, options(), policy)
            .unwrap();
        let received = scattering
            .connect_bundle_with_transmittance_samples(
                cx,
                initial,
                receiver,
                options(),
                policy,
                CurvedTransmittanceOptions::default(),
                |_| panic!("uniform physical segments need no Poisson sample"),
            )
            .unwrap();
        assert_eq!(*received.connection(), original);
        assert_eq!(received.connection().path().surface_events(), 2);
        let bundle = received.connection().path().bundle();
        near(received.log_scattering_survival(), -1.8, 1e-7);
        near(bundle.propagation().propagation().optical_depth, 0.7, 1e-7);
        near(bundle.interface_log_power().exp(), 80.0 / 81.0, 1e-12);
        let area = (1.0_f64 + 1.2 / 1.5 * 2.0).powi(2);
        near(
            bundle
                .propagation()
                .project_at_endpoint(receiver.normal(), options().projection)
                .unwrap()
                .area_per_solid_angle_m2(),
            area,
            5e-7,
        );
        near(
            bundle
                .point_source_irradiance(receiver.normal(), options().projection, 1.0)
                .unwrap()
                * received.log_scattering_survival().exp(),
            (80.0 / 81.0) * (-2.5_f64).exp() / area,
            1e-8,
        );
        assert!(matches!(
            scattering.connect_bundle_with_transmittance_samples(
                cx,
                initial,
                receiver,
                options(),
                BundleBranchPolicy::Sequence(&[DielectricEvent::Reflection]),
                CurvedTransmittanceOptions::default(),
                |_| panic!("failed deterministic branch must not consume shadow RNG"),
            ),
            Err(CurvedTransmittanceError::Connection(_))
        ));
    });
}

#[test]
fn g0_source_tracking_refuses_exceeded_majorants_and_invalid_variates() {
    with_cx(|cx, _| {
        let point = Point3::new(0.0, 0.0, 2.0);
        let geometry = InstanceScene::try_new(vec![plane(9, point)]).unwrap();
        let material = glass(1.2, 0.0);
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
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        for (majorant, local) in [(1.0, 1.01), (1.0, -0.1), (1.0, f64::NAN), (f64::NAN, 0.5)] {
            let law = LocalLaw::constant(majorant, local);
            let volumes = [Some(CurvedRayleighBinding::Field(&law))];
            let scattering =
                CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
            assert!(matches!(
                scattering.connect_bundle_with_transmittance_samples(
                    cx,
                    launch(Vec3::new(0.0, 0.0, 1.0)),
                    target(point),
                    options(),
                    BundleBranchPolicy::Sequence(&[]),
                    CurvedTransmittanceOptions::default(),
                    |_| flight_sample(1.0, 0.3),
                ),
                Err(CurvedTransmittanceError::Volume(
                    VolumeTransportError::ScatteringSpectrum
                ))
            ));
        }
        let law = LocalLaw::constant(1.0, 0.5);
        let volumes = [Some(CurvedRayleighBinding::Field(&law))];
        let scattering = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
        for sample in [f64::NAN, -0.1, 1.0] {
            assert!(matches!(
                scattering.connect_bundle_with_transmittance_samples(
                    cx,
                    launch(Vec3::new(0.0, 0.0, 1.0)),
                    target(point),
                    options(),
                    BundleBranchPolicy::Sequence(&[]),
                    CurvedTransmittanceOptions::default(),
                    |_| sample,
                ),
                Err(CurvedTransmittanceError::Volume(
                    VolumeTransportError::Sample
                ))
            ));
        }
        // Vacuum density cannot extend the original molecular source's
        // wavelength support, even though the optical field admits 650 nm.
        let spectrum = MassScatteringSpectrum::source_table(
            &[(500.0, 0.5), (600.0, 0.5)],
            ContentHash([180; 32]),
        )
        .unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 0.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
        let mut initial = launch(Vec3::new(0.0, 0.0, 1.0));
        initial.wavelength_nm = 650.0;
        assert_eq!(
            scattering
                .connect_bundle_with_transmittance_samples(
                    cx,
                    initial,
                    target(point),
                    options(),
                    BundleBranchPolicy::Sequence(&[]),
                    CurvedTransmittanceOptions::default(),
                    |_| panic!("source-domain refusal precedes any shadow draw"),
                )
                .unwrap_err(),
            CurvedTransmittanceError::Volume(VolumeTransportError::Wavelength)
        );
    });
}

#[test]
fn g0_source_tracking_limits_and_cancellation_refuse_partial_estimates() {
    with_cx(|cx, gate| {
        let point = Point3::new(0.0, 0.0, 2.0);
        let geometry = InstanceScene::try_new(vec![plane(9, point)]).unwrap();
        let material = glass(1.2, 0.0);
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
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        let law = LocalLaw::constant(2.0, 0.5);
        let volumes = [Some(CurvedRayleighBinding::Field(&law))];
        let scattering = CurvedRayleighScene::try_new_with_fields(cx, &scene, &volumes).unwrap();
        let initial = launch(Vec3::new(0.0, 0.0, 1.0));
        let receiver = target(point);
        let policy = BundleBranchPolicy::Sequence(&[]);
        for (tracking, expected) in [
            (
                CurvedTransmittanceOptions {
                    max_free_flight_draws: 1,
                    max_tracking_candidates: 10,
                },
                CurvedTransmittanceError::DrawLimit,
            ),
            (
                CurvedTransmittanceOptions {
                    max_free_flight_draws: 10,
                    max_tracking_candidates: 1,
                },
                CurvedTransmittanceError::CandidateLimit,
            ),
        ] {
            assert_eq!(
                scattering
                    .connect_bundle_with_transmittance_samples(
                        cx,
                        initial,
                        receiver,
                        options(),
                        policy,
                        tracking,
                        |_| flight_sample(2.0, 0.3),
                    )
                    .unwrap_err(),
                expected
            );
        }
        let mut limited = options();
        limited.max_trace_evaluations = 1;
        assert!(matches!(
            scattering.connect_bundle_with_transmittance_samples(
                cx,
                initial,
                receiver,
                limited,
                policy,
                CurvedTransmittanceOptions::default(),
                |_| panic!("no replay budget must not consume a source sample"),
            ),
            Err(CurvedTransmittanceError::Connection(
                BundleConnectionError::TraceBudgetExceeded
            ))
        ));
        gate.request();
        assert!(matches!(
            scattering.connect_bundle_with_transmittance_samples(
                cx,
                initial,
                receiver,
                options(),
                policy,
                CurvedTransmittanceOptions::default(),
                |_| panic!("cancelled work must not consume a source sample"),
            ),
            Err(CurvedTransmittanceError::Cancelled)
                | Err(CurvedTransmittanceError::Connection(
                    BundleConnectionError::Cancelled
                ))
        ));
    });
}

mod images;
