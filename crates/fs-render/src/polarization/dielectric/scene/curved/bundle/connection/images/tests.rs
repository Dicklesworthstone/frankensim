use super::super::super::tests::{
    ball_reference, bindings, close, close_vec, glass, plane, reference_columns, settings, sphere,
    with_cx,
};
use super::super::{
    BundleBranchPolicy, BundleConnectionError, BundleConnectionLaunch, BundleConnectionOptions,
    BundleConnectionTarget,
};
use super::*;
use crate::charts::Ray;
use crate::dielectric::DielectricEvent;
use crate::gradient_index::bundle::ProjectionOptions;
use crate::instances::InstanceScene;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::PathTermination;
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use fs_exec::{CancelGate, Cx};
use fs_geom::{Point3, Vec3};

fn ray(angle: f64) -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, -3.0),
        dir: Vec3::new(angle.sin(), 0.0, angle.cos()),
    }
}

fn launch() -> BundleConnectionLaunch {
    BundleConnectionLaunch {
        ray: ray(0.0),
        transverse_axis: Vec3::new(1.0, 0.0, 0.0),
        wavelength_nm: 550.0,
        normalized_source_stokes: Stokes::unpolarized(1.0),
        initial_medium: 1,
    }
}

fn target() -> BundleConnectionTarget {
    BundleConnectionTarget::try_new(
        9,
        Point3::new(0.03, 0.0, 2.5),
        Vec3::new(0.0, 0.0, -1.0),
        2e-7,
    )
    .unwrap()
}

fn connection_options() -> BundleConnectionOptions {
    BundleConnectionOptions {
        scene: settings(),
        ..BundleConnectionOptions::default()
    }
}

fn image_options() -> BundleImageOptions {
    BundleImageOptions {
        min_separation_rad: 0.01,
        duplicate_tolerance_rad: 1e-5,
        ..BundleImageOptions::default()
    }
}

fn seed(angle: f64) -> BundleImageSeed<'static> {
    BundleImageSeed {
        initial_direction: ray(angle).dir,
        policy: BundleBranchPolicy::TransmitOrReflectTir,
    }
}

fn with_ball(f: impl FnOnce(&Cx<'_>, &CancelGate, &CurvedPolarizedScene<'_>)) {
    with_cx(|cx, gate| {
        let geometry =
            InstanceScene::try_new(vec![sphere(), plane(9, target().position_m(), 0.5)]).unwrap();
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
            bindings(Vec3::new(0.0, 1.0, 0.0)),
            100_000,
        )
        .unwrap();
        f(cx, gate, &scene);
    });
}

/// Scalar bisection of independent quadratic intersections and vector Snell.
/// None of the native connection solver, event walk or analytic derivatives
/// participates in establishing these three different roots.
fn reference_root(mut lo: f64, mut hi: f64) -> Ray {
    let residual = |angle| ball_reference(ray(angle), 1.5).0.origin.x - 0.03;
    let mut at_lo = residual(lo);
    assert!(at_lo * residual(hi) < 0.0);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        let value = residual(mid);
        if value == 0.0 {
            return ray(mid);
        }
        if value * at_lo > 0.0 {
            lo = mid;
            at_lo = value;
        } else {
            hi = mid;
        }
    }
    ray(0.5 * (lo + hi))
}

#[test]
fn g1_multistart_ball_lens_retains_three_regular_same_history_images() {
    with_ball(|cx, _, scene| {
        let seeds = [seed(-0.18), seed(0.065), seed(0.14)];
        let catalog = scene
            .discover_bundle_images(
                cx,
                launch(),
                target(),
                connection_options(),
                &seeds,
                image_options(),
            )
            .unwrap();
        assert_eq!(catalog.images().len(), 3);
        let expected = [
            reference_root(-0.20, -0.18),
            reference_root(0.06, 0.08),
            reference_root(0.12, 0.14),
        ];
        let signature = catalog.images()[0].connection().path().branch_signature();
        let mut signed_areas = [0.0; 3];
        for (index, (image, expected_ray)) in catalog.images().iter().zip(expected).enumerate() {
            let connection = image.connection();
            let path = connection.path();
            assert_eq!(image.first_seed(), index);
            assert_eq!(path.branch_signature(), signature);
            assert_eq!(path.surface_events(), 3);
            assert_eq!(path.termination(), PathTermination::Detected);
            assert!(connection.iterations() > 0);
            assert!(connection.residual_m() <= target().position_tolerance_m());
            assert!(image.angular_resolution_rad().is_finite());
            assert!(image.angular_resolution_rad() < image_options().duplicate_tolerance_rad);
            let geometric = path.bundle().propagation();
            let actual_ray = geometric.initial_ray();
            close_vec(actual_ray.dir, expected_ray.dir, 2e-6);
            let (endpoint, outside, inside, transfer) = ball_reference(actual_ray, 1.5);
            close_vec(geometric.propagation().ray.dir, endpoint.dir, 2e-7);
            close_vec(
                geometric
                    .propagation()
                    .ray
                    .origin
                    .delta_from(endpoint.origin),
                Vec3::new(0.0, 0.0, 0.0),
                2e-7,
            );
            close(geometric.propagation().distance_m, outside + inside, 2e-7);
            close(geometric.propagation().optical_depth, 0.2 * inside, 5e-8);
            close(path.bundle().interface_log_power().exp(), transfer, 2e-7);
            let columns = reference_columns(actual_ray, |r| ball_reference(r, 1.5).0.origin);
            let signed_area = -(columns[0].x * columns[1].y - columns[0].y * columns[1].x);
            let projection = geometric
                .project_at_endpoint(target().normal(), ProjectionOptions::default())
                .unwrap();
            for (actual, expected) in projection.angular_derivatives().iter().zip(columns) {
                close_vec(actual.position_m_per_rad, expected, 5e-6);
            }
            close(
                projection.signed_area_per_solid_angle_m2(),
                signed_area,
                5e-6,
            );
            signed_areas[index] = signed_area;
            let expected_irradiance = transfer * (-0.2 * inside).exp() / signed_area.abs();
            close(
                path.bundle()
                    .point_source_irradiance(target().normal(), ProjectionOptions::default(), 1.0)
                    .unwrap(),
                expected_irradiance,
                5e-5 * expected_irradiance,
            );
        }
        // One radial image changes parity. Taking abs(det J) retains
        // all positive powers; a signed sum would incorrectly cancel images.
        assert!(signed_areas[0] * signed_areas[1] > 0.0);
        assert!(signed_areas[1] * signed_areas[2] < 0.0);
    });
}

#[test]
fn g3_multistart_deduplicates_equivalent_seeds_and_reports_real_misses_and_work() {
    with_ball(|cx, _, scene| {
        let mut seeds = [
            seed(-0.18),
            seed(0.065),
            seed(0.14),
            seed(0.075),
            seed(-0.185),
            seed(0.125),
            seed(0.0),
        ];
        seeds[3].policy = BundleBranchPolicy::Sequence(&[
            DielectricEvent::Transmission,
            DielectricEvent::Transmission,
        ]);
        seeds[6].initial_direction = Vec3::new(0.0, 1.0, 0.0);
        let discover = || {
            scene
                .discover_bundle_images(
                    cx,
                    launch(),
                    target(),
                    connection_options(),
                    &seeds,
                    image_options(),
                )
                .unwrap()
        };
        let catalog = discover();
        assert_eq!(catalog.images().len(), 3);
        let outcomes = catalog.seed_results();
        assert_eq!(outcomes.len(), seeds.len());
        for (seed, image) in [(0, 0), (1, 1), (2, 2)] {
            assert!(
                matches!(outcomes[seed].outcome(), BundleImageSeedOutcome::NewImage { image: i } if *i == image)
            );
        }
        for (seed, image) in [(3, 1), (4, 0), (5, 2)] {
            assert!(
                matches!(outcomes[seed].outcome(), BundleImageSeedOutcome::Duplicate { image: i } if *i == image)
            );
        }
        assert!(matches!(
            outcomes[6].outcome(),
            BundleImageSeedOutcome::Miss(BundleConnectionError::ReceiverNotReached(
                PathTermination::Escaped
            ))
        ));
        assert_eq!(
            catalog.trace_evaluations(),
            outcomes.iter().map(|o| o.trace_evaluations()).sum::<u32>()
        );
        assert_eq!(
            catalog.field_evaluations(),
            outcomes.iter().map(|o| o.field_evaluations()).sum::<u64>()
        );
        assert_eq!(
            catalog.query_work(),
            outcomes.iter().map(|o| o.query_work()).sum::<u64>()
        );
        assert!(outcomes[6].trace_evaluations() > 0);
        assert!(
            catalog.trace_evaluations()
                > catalog
                    .images()
                    .iter()
                    .map(|i| i.connection().trace_evaluations())
                    .sum::<u32>()
        );
        let replay = discover();
        assert_eq!(catalog.trace_evaluations(), replay.trace_evaluations());
        assert_eq!(catalog.field_evaluations(), replay.field_evaluations());
        assert_eq!(catalog.query_work(), replay.query_work());
        for (actual, repeated) in catalog.images().iter().zip(replay.images()) {
            assert_eq!(actual.connection(), repeated.connection());
        }
    });
}

#[test]
fn g0_multistart_refuses_unresolved_image_separation_and_invalid_admission() {
    with_ball(|cx, _, scene| {
        let discover = |seeds: &[BundleImageSeed<'_>], options| {
            scene.discover_bundle_images(
                cx,
                launch(),
                target(),
                connection_options(),
                seeds,
                options,
            )
        };
        assert!(matches!(
            discover(&[], image_options()),
            Err(BundleImageError::EmptySeeds)
        ));
        assert!(matches!(
            discover(&[seed(0.065); MAX_BUNDLE_IMAGE_SEEDS + 1], image_options()),
            Err(BundleImageError::TooManySeeds)
        ));
        assert!(matches!(
            discover(
                &[seed(0.065)],
                BundleImageOptions {
                    duplicate_tolerance_rad: f64::NAN,
                    ..image_options()
                }
            ),
            Err(BundleImageError::InvalidOptions)
        ));
        // These two exact Snell roots are 0.058 rad apart: neither a duplicate
        // nor a resolved pair under the requested 0.1 rad separation allowance.
        assert!(matches!(
            discover(&[seed(0.065), seed(0.14)], BundleImageOptions {
                min_separation_rad: 0.1,
                ..image_options()
            }),
            Err(BundleImageError::AmbiguousImage { seed: 1, image: 0, separation_rad })
                if separation_rad > 0.05 && separation_rad < 0.07
        ));
        // Opposite native footprint parities prove these are different regular
        // images even when an overly broad duplicate tolerance overlaps them.
        // They must remain unresolved, never silently lose one image's power.
        assert!(matches!(
            discover(
                &[seed(0.065), seed(0.14)],
                BundleImageOptions {
                    duplicate_tolerance_rad: 0.08,
                    min_separation_rad: 0.1,
                    ..image_options()
                }
            ),
            Err(BundleImageError::AmbiguousImage {
                seed: 1,
                image: 0,
                ..
            })
        ));
        let no_direction = BundleImageSeed {
            initial_direction: Vec3::new(0.0, 0.0, 0.0),
            ..seed(0.0)
        };
        assert!(matches!(
            discover(&[no_direction], image_options()),
            Err(BundleImageError::Seed { seed: 0, .. })
        ));
    });
}

#[test]
fn g4_multistart_enforces_total_work_capacity_and_cancellation_without_partial_catalogs() {
    with_ball(|cx, gate, scene| {
        let seeds = [seed(-0.18), seed(0.065), seed(0.14)];
        let discover = |options| {
            scene.discover_bundle_images(
                cx,
                launch(),
                target(),
                connection_options(),
                &seeds,
                options,
            )
        };
        assert!(matches!(
            discover(BundleImageOptions {
                max_trace_evaluations: 1,
                ..image_options()
            }),
            Err(BundleImageError::TraceBudgetExceeded { seed: 0 })
        ));
        assert!(matches!(
            discover(BundleImageOptions {
                max_images: 2,
                ..image_options()
            }),
            Err(BundleImageError::ImageCapacity {
                seed: 2,
                capacity: 2
            })
        ));
        // Per-seed nonconvergence is reported explicitly with its spent work.
        // It is not a complete catalog and does not consume another seed's cap.
        let limited = scene
            .discover_bundle_images(
                cx,
                launch(),
                target(),
                BundleConnectionOptions {
                    max_trace_evaluations: 1,
                    ..connection_options()
                },
                &seeds[..1],
                image_options(),
            )
            .unwrap();
        assert!(limited.images().is_empty());
        assert_eq!(limited.trace_evaluations(), 1);
        assert!(matches!(
            limited.seed_results()[0].outcome(),
            BundleImageSeedOutcome::Miss(BundleConnectionError::TraceBudgetExceeded)
        ));
        gate.request();
        assert!(matches!(
            discover(image_options()),
            Err(BundleImageError::Seed {
                seed: 0,
                error: BundleConnectionError::Cancelled
            })
        ));
    });
}
