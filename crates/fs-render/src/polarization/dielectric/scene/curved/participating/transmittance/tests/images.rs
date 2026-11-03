//! Native multi-image source survival: one solve per seed, one replay per image.

use super::*;
use crate::polarization::dielectric::scene::curved::bundle::connection::images::{
    BundleImageOptions, BundleImageSeed, BundleImageSeedOutcome,
};
use crate::polarization::dielectric::scene::curved::participating::CurvedDensityField;
use crate::polarization::dielectric::scene::curved::participating::transmittance::images::{
    CurvedTransmittanceImageError, CurvedTransmittanceImageWorkspace,
};
use fs_geom::fixtures::SphereChart;

fn lens_launch() -> BundleConnectionLaunch {
    BundleConnectionLaunch {
        ray: Ray {
            origin: Point3::new(0.0, 0.0, -3.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        },
        initial_medium: 1,
        ..launch(Vec3::new(0.0, 0.0, 1.0))
    }
}

fn lens_target() -> BundleConnectionTarget {
    target(Point3::new(0.03, 0.0, 2.5))
}

fn image_seed(angle: f64) -> BundleImageSeed<'static> {
    BundleImageSeed {
        initial_direction: Vec3::new(angle.sin(), 0.0, angle.cos()),
        policy: BundleBranchPolicy::TransmitOrReflectTir,
    }
}

fn lens_seeds() -> [BundleImageSeed<'static>; 6] {
    let mut seeds = [-0.18, 0.065, 0.14, 0.075, -0.185, 0.125].map(image_seed);
    // Equivalent physical history expressed by a different policy must still
    // share the image, including when it supplies the better representative.
    seeds[3].policy = BundleBranchPolicy::Sequence(&[
        DielectricEvent::Transmission,
        DielectricEvent::Transmission,
    ]);
    seeds
}

fn image_options() -> BundleImageOptions {
    BundleImageOptions {
        max_images: 3,
        max_trace_evaluations: 512,
        min_separation_rad: 0.01,
        duplicate_tolerance_rad: 1e-5,
        ..BundleImageOptions::default()
    }
}

fn with_lens(f: impl FnOnce(&Cx<'_>, &CurvedPolarizedScene<'_>)) {
    with_cx(|cx, _| {
        let sphere = GeometryInstance::try_new(
            1,
            hash_domain("molecular-image-survival-test", b"native-unit-sphere"),
            SharedGeometry::chart(SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            }),
            RigidTransform::identity(),
        )
        .unwrap();
        let geometry =
            InstanceScene::try_new(vec![sphere, plane(9, lens_target().position_m())]).unwrap();
        let lens = glass(1.5, 0.2);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        outside: 1,
                        inside: 0,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                detector(),
            ],
            100_000,
        )
        .unwrap();
        f(cx, &scene);
    });
}

#[test]
fn g1_native_lens_images_replay_once_and_reuse_prepared_survival_storage() {
    with_lens(|cx, scene| {
        let spectrum = MassScatteringSpectrum::constant(0.4).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let seeds = lens_seeds();
        let mut workspace = CurvedTransmittanceImageWorkspace::try_new(3, seeds.len()).unwrap();
        let no_random_work = CurvedTransmittanceOptions {
            max_free_flight_draws: 0,
            max_tracking_candidates: 0,
        };
        let mut stored_address = None;
        // Clear keeps the native solve exactly, while uniform scattering uses
        // the same actual geometric arc without consuming any random samples.
        for density in [1.0, 0.0, 1.0] {
            let law = CurvedRayleighMedium::try_new(&spectrum, density, &phase).unwrap();
            let media = [Some(law), Some(law)];
            let scattering = CurvedRayleighScene::try_new(cx, scene, &media).unwrap();
            let catalog = scattering
                .discover_bundle_images_with_transmittance_samples(
                    cx,
                    lens_launch(),
                    lens_target(),
                    options(),
                    &seeds,
                    image_options(),
                    no_random_work,
                    &mut workspace,
                    |_, _| panic!("uniform or clear images must not draw source samples"),
                )
                .unwrap();
            assert_eq!(catalog.images().len(), 3);
            let geometric = catalog.geometric_catalog();
            assert_eq!(geometric.seed_results().len(), seeds.len());
            for (seed, image) in [(3, 1), (4, 0), (5, 2)] {
                assert!(matches!(
                    geometric.seed_results()[seed].outcome(),
                    BundleImageSeedOutcome::Duplicate { image: found } if *found == image
                ));
            }
            assert_eq!(catalog.free_flight_draws(), 0);
            assert_eq!(catalog.tracking_candidates(), 0);
            assert_eq!(
                catalog.trace_evaluations(),
                geometric.trace_evaluations() + if density > 0.0 { 3 } else { 0 }
            );
            let mut fields = geometric.field_evaluations();
            let mut queries = geometric.query_work();
            for (survival, image) in catalog.images().iter().zip(geometric.images()) {
                assert_eq!(survival.connection(), image.connection());
                let path = survival.connection().path();
                assert_eq!(path.termination(), PathTermination::Detected);
                assert_eq!(path.surface_events(), 3);
                assert!(!survival.ratio_tracking_zero());
                near(
                    survival.log_scattering_survival(),
                    -0.4 * density * path.bundle().propagation().propagation().distance_m,
                    5e-8,
                );
                fields += survival.field_evaluations() - image.connection().field_evaluations();
                queries += survival.query_work() - image.connection().query_work();
            }
            assert_eq!(catalog.field_evaluations(), fields);
            assert_eq!(catalog.query_work(), queries);
            if let Some(address) = stored_address {
                assert_eq!(catalog.images().as_ptr(), address);
            }
            stored_address = Some(catalog.images().as_ptr());
        }
        assert_eq!(workspace.image_capacity(), 3);
        assert_eq!(workspace.seed_capacity(), seeds.len());
    });
}

struct HalfDensity;

impl CurvedDensityField for HalfDensity {
    fn density_kg_m3(&self, _: Point3) -> Result<f64, VolumeTransportError> {
        Ok(0.5)
    }
}

#[test]
fn g4_source_work_limits_are_cumulative_across_distinct_images_and_fail_closed() {
    with_lens(|cx, scene| {
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let law = CurvedRayleighMedium::try_new_heterogeneous(
            &spectrum,
            &HalfDensity,
            1.0,
            ContentHash([179; 32]),
            &phase,
        )
        .unwrap();
        let media = [Some(law), Some(law)];
        let scattering = CurvedRayleighScene::try_new(cx, scene, &media).unwrap();
        let seeds = lens_seeds();
        let mut workspace = CurvedTransmittanceImageWorkspace::try_new(3, seeds.len()).unwrap();
        let sample = |flight| flight_sample(1.0, if flight == 0 { 0.25 } else { 16.0 });
        // Native discovery must finish before even the first survival draw.
        let too_few_images = BundleImageOptions {
            max_images: 2,
            ..image_options()
        };
        assert!(matches!(
            scattering.discover_bundle_images_with_transmittance_samples(
                cx,
                lens_launch(),
                lens_target(),
                options(),
                &seeds,
                too_few_images,
                CurvedTransmittanceOptions::default(),
                &mut workspace,
                |_, _| panic!("failed discovery must not draw source samples"),
            ),
            Err(CurvedTransmittanceImageError::Discovery(_))
        ));
        let mut keys = Vec::new();
        let (geometric_traces, total_draws) = {
            let catalog = scattering
                .discover_bundle_images_with_transmittance_samples(
                    cx,
                    lens_launch(),
                    lens_target(),
                    options(),
                    &seeds,
                    image_options(),
                    CurvedTransmittanceOptions::default(),
                    &mut workspace,
                    |image, flight| {
                        keys.push((image, flight));
                        sample(flight)
                    },
                )
                .unwrap();
            assert_eq!(catalog.images().len(), 3);
            assert_eq!(catalog.tracking_candidates(), 3);
            let mut offset = 0_usize;
            for (image, survival) in catalog.images().iter().enumerate() {
                assert_eq!(survival.tracking_candidates(), 1);
                near(survival.log_scattering_survival().exp(), 0.5, 1e-12);
                for flight in 0..survival.free_flight_draws() {
                    assert_eq!(keys[offset], (image as u32, flight));
                    offset += 1;
                }
            }
            assert_eq!(keys.len(), offset);
            assert_eq!(catalog.free_flight_draws() as usize, offset);
            assert_eq!(
                catalog.trace_evaluations(),
                catalog.geometric_catalog().trace_evaluations() + 3
            );
            (
                catalog.geometric_catalog().trace_evaluations(),
                catalog.free_flight_draws(),
            )
        };
        // Exactly the complete discovery's work leaves no replay budget,
        // even though each individual connection's own trace cap is sufficient.
        let no_replay = BundleImageOptions {
            max_trace_evaluations: geometric_traces,
            ..image_options()
        };
        assert_eq!(
            scattering
                .discover_bundle_images_with_transmittance_samples(
                    cx,
                    lens_launch(),
                    lens_target(),
                    options(),
                    &seeds,
                    no_replay,
                    CurvedTransmittanceOptions::default(),
                    &mut workspace,
                    |_, _| panic!("an exhausted aggregate trace cap must precede RNG"),
                )
                .unwrap_err(),
            CurvedTransmittanceImageError::TraceBudgetExceeded { image: 0 }
        );
        for (tracking, expected) in [
            (
                CurvedTransmittanceOptions {
                    max_free_flight_draws: total_draws - 1,
                    max_tracking_candidates: 3,
                },
                CurvedTransmittanceError::DrawLimit,
            ),
            (
                CurvedTransmittanceOptions {
                    max_free_flight_draws: total_draws,
                    max_tracking_candidates: 2,
                },
                CurvedTransmittanceError::CandidateLimit,
            ),
        ] {
            assert_eq!(
                scattering
                    .discover_bundle_images_with_transmittance_samples(
                        cx,
                        lens_launch(),
                        lens_target(),
                        options(),
                        &seeds,
                        image_options(),
                        tracking,
                        &mut workspace,
                        |_, flight| sample(flight),
                    )
                    .unwrap_err(),
                CurvedTransmittanceImageError::Replay {
                    image: 2,
                    error: expected,
                }
            );
        }
        // Failure never publishes partial images and does not poison reuse.
        let complete = scattering
            .discover_bundle_images_with_transmittance_samples(
                cx,
                lens_launch(),
                lens_target(),
                options(),
                &seeds,
                image_options(),
                CurvedTransmittanceOptions {
                    max_free_flight_draws: total_draws,
                    max_tracking_candidates: 3,
                },
                &mut workspace,
                |_, flight| sample(flight),
            )
            .unwrap();
        assert_eq!(complete.images().len(), 3);
        assert_eq!(complete.free_flight_draws(), total_draws);
        assert_eq!(complete.tracking_candidates(), 3);
    });
}

#[test]
fn g1_duplicate_curved_seeds_share_one_independent_native_survival_replay() {
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
        let seeds = [image_seed(0.03), image_seed(0.0)];
        let mut workspace = CurvedTransmittanceImageWorkspace::try_new(1, 2).unwrap();
        let mut draws = Vec::new();
        let catalog = scattering
            .discover_bundle_images_with_transmittance_samples(
                cx,
                launch(Vec3::new(0.0, 0.0, 1.0)),
                target(point),
                options(),
                &seeds,
                BundleImageOptions {
                    max_images: 1,
                    ..image_options()
                },
                CurvedTransmittanceOptions {
                    max_free_flight_draws: 3,
                    max_tracking_candidates: 2,
                },
                &mut workspace,
                |image, flight| {
                    draws.push((image, flight));
                    flight_sample(2.0, [0.4, 0.9, 3.0][flight as usize])
                },
            )
            .unwrap();
        assert_eq!(catalog.images().len(), 1);
        assert!(matches!(
            catalog.geometric_catalog().seed_results()[1].outcome(),
            BundleImageSeedOutcome::Duplicate { image: 0 }
        ));
        assert_eq!(draws, [(0, 0), (0, 1), (0, 2)]);
        assert_eq!(catalog.free_flight_draws(), 3);
        assert_eq!(catalog.tracking_candidates(), 2);
        assert_eq!(
            catalog.trace_evaluations(),
            catalog.geometric_catalog().trace_evaluations() + 1
        );
        let survival = &catalog.images()[0];
        assert_eq!(
            survival.connection(),
            catalog.geometric_catalog().images()[0].connection()
        );
        near(survival.log_scattering_survival().exp(), 0.375, 1e-12);
        let queries = law.queries.lock().unwrap();
        assert_eq!(queries.len(), 2, "duplicate shooting seeds must not retrack");
        for (actual, arc) in queries.iter().zip([0.4, 1.3]) {
            let expected = affine_point(arc);
            near(actual.x, expected.x, 5e-7);
            near(actual.y, expected.y, 5e-7);
            near(actual.z, expected.z, 5e-7);
        }
        near(
            survival
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
