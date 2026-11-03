//! Finite light proposals may omit resolved geometry, never unresolved work.

use super::*;
use super::super::super::tests::{glass, plane, settings, with_cx};
use crate::dielectric::DielectricSurface;
use crate::gradient_index::event::BoundarySide;
use crate::instances::InstanceScene;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::{Interaction, SceneError, SurfaceBinding};
use crate::polarization::dielectric::scene::curved::SceneMedium;
use fs_geom::Point3;

fn source_binding() -> SurfaceBinding {
    SurfaceBinding { object_id: 2, interaction: Interaction::Dielectric {
        outside: 0, inside: 1, surface: DielectricSurface::SMOOTH,
    } }
}

fn detector(id: u64) -> SurfaceBinding {
    SurfaceBinding { object_id: id, interaction: Interaction::Detector {
        reference_axis: Vec3::new(0.0, 1.0, 0.0),
    } }
}

fn launch() -> BundleConnectionLaunch {
    BundleConnectionLaunch {
        ray: Ray { origin: Point3::new(0.0, 0.0, 0.0), dir: Vec3::new(0.0, 0.0, 1.0) },
        transverse_axis: Vec3::new(0.0, 1.0, 0.0),
        wavelength_nm: 550.0,
        normalized_source_stokes: Stokes::unpolarized(1.0),
        initial_medium: 0,
    }
}

fn seed(direction: Vec3) -> BundleImageSeed<'static> {
    BundleImageSeed { initial_direction: direction, policy: BundleBranchPolicy::Sequence(&[]) }
}

fn target(normal: Vec3) -> BundleConnectionTarget {
    BundleConnectionTarget::try_new(9, Point3::new(0.0, 0.0, 2.0), normal, 2e-7).unwrap()
}

fn options() -> BundleConnectionOptions {
    BundleConnectionOptions { scene: settings(), ..BundleConnectionOptions::default() }
}

#[test]
fn g1_source_side_support_preserves_native_owner_medium_and_actual_work() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            plane(2, Point3::new(0.0, 0.0, 0.0), 1.0),
            plane(9, Point3::new(0.0, 0.0, 2.0), 0.25),
        ]).unwrap();
        let air = glass(1.0, 0.0);
        // Equal optical values do not merge original scene-medium identities.
        let media = [SceneMedium::Homogeneous(&air), SceneMedium::Homogeneous(&air)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![source_binding(), detector(9)], 100_000).unwrap();
        let source = BundleSurfaceSource::try_new(2, BoundarySide::Positive).unwrap();
        let seeds = [seed(Vec3::new(0.0, 0.0, -1.0)), seed(Vec3::new(0.0, 0.0, 1.0))];
        let receiver = target(Vec3::new(0.0, 0.0, -1.0));
        let mut workspace = BundleImageWorkspace::try_new(1, 2).unwrap();
        let catalog = scene.discover_bundle_images_from_surface_with_workspace(cx, launch(),
            source, receiver, options(), &seeds, BundleImageOptions { max_images: 1,
                ..BundleImageOptions::default() }, &mut workspace).unwrap();
        catalog.validate_proposal_support().unwrap();
        assert_eq!(catalog.images().len(), 1);
        assert_eq!(catalog.images()[0].first_seed(), 1);
        let excluded = &catalog.seed_results()[0];
        let BundleImageSeedOutcome::Miss(BundleConnectionError::Scene(
            CurvedBundleSceneError::SourceSideExcluded {
                object_id, field_evaluations, query_work,
            }
        )) = excluded.outcome() else { panic!("native wrong-side seed must be explicit"); };
        assert_eq!(*object_id, 2);
        assert!(*field_evaluations > 0 && *query_work > 0);
        assert_eq!(excluded.trace_evaluations(), 1);
        assert_eq!(excluded.field_evaluations(), *field_evaluations);
        assert_eq!(excluded.query_work(), *query_work);
        assert_eq!(catalog.trace_evaluations(), 2);
        assert_eq!(catalog.field_evaluations(), catalog.seed_results().iter()
            .map(BundleImageSeedResult::field_evaluations).sum::<u64>());
        assert_eq!(catalog.query_work(), catalog.seed_results().iter()
            .map(BundleImageSeedResult::query_work).sum::<u64>());

        let empty = scene.discover_bundle_images_from_surface_with_workspace(cx, launch(),
            source, receiver, options(), &seeds[..1], BundleImageOptions { max_images: 1,
                ..BundleImageOptions::default() }, &mut workspace).unwrap();
        assert!(empty.images().is_empty());
        empty.validate_proposal_support().unwrap();
        // Source-side exclusion cannot hide an invalid selected-side medium.
        let bad = scene.discover_bundle_images_from_surface_with_workspace(cx,
            BundleConnectionLaunch { initial_medium: 1, ..launch() }, source, receiver,
            options(), &seeds[..1], BundleImageOptions { max_images: 1,
                ..BundleImageOptions::default() }, &mut workspace).unwrap_err();
        assert_eq!(bad, BundleImageError::Seed { seed: 0,
            error: BundleConnectionError::Scene(CurvedBundleSceneError::Scene(
                SceneError::MediumMismatch { object_id: 2, expected: 0, actual: 1 }.into(),
            )),
        });
        let wrong_source = BundleSurfaceSource::try_new(9, BoundarySide::Positive).unwrap();
        assert!(matches!(scene.discover_bundle_images_from_surface_with_workspace(cx,
            launch(), wrong_source, receiver, options(), &seeds[..1], BundleImageOptions {
                max_images: 1, ..BundleImageOptions::default()
            }, &mut workspace), Err(BundleImageError::Seed { error:
                BundleConnectionError::Scene(CurvedBundleSceneError::InvalidSurfaceSource {
                    object_id: 9,
                }), .. })));
    });
}

#[test]
fn g1_completed_occlusion_escape_and_wrong_receiver_define_only_proposal_support() {
    with_cx(|cx, _| {
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air), SceneMedium::Homogeneous(&air)];
        for interaction in [Interaction::Opaque, detector(5).interaction] {
            let geometry = InstanceScene::try_new(vec![
                plane(2, Point3::new(0.0, 0.0, 0.0), 1.0),
                plane(5, Point3::new(0.0, 0.0, 1.0), 0.25),
                plane(9, Point3::new(0.0, 0.0, 2.0), 0.25),
            ]).unwrap();
            let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
                vec![source_binding(), SurfaceBinding { object_id: 5, interaction }, detector(9)],
                100_000).unwrap();
            let mut workspace = BundleImageWorkspace::try_new(1, 1).unwrap();
            let catalog = scene.discover_bundle_images_from_surface_with_workspace(cx, launch(),
                BundleSurfaceSource::try_new(2, BoundarySide::Positive).unwrap(),
                target(Vec3::new(0.0, 0.0, -1.0)), options(),
                &[seed(Vec3::new(0.0, 0.0, 1.0))], BundleImageOptions {
                    max_images: 1, ..BundleImageOptions::default()
                }, &mut workspace).unwrap();
            assert!(catalog.images().is_empty());
            let expected = if matches!(interaction, Interaction::Opaque) {
                BundleConnectionError::ReceiverNotReached(PathTermination::Opaque)
            } else { BundleConnectionError::WrongDetector };
            assert_eq!(catalog.seed_results()[0].outcome(), &BundleImageSeedOutcome::Miss(expected));
            catalog.validate_proposal_support().unwrap();
            assert!(catalog.field_evaluations() > 0 && catalog.query_work() > 0);
        }
        let geometry = InstanceScene::try_new(vec![
            plane(2, Point3::new(0.0, 0.0, 0.0), 1.0),
            plane(9, Point3::new(0.0, 0.0, 2.0), 0.25),
        ]).unwrap();
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![source_binding(), detector(9)], 100_000).unwrap();
        let mut workspace = BundleImageWorkspace::try_new(1, 1).unwrap();
        for (direction, normal, expected) in [
            (Vec3::new(1.0, 0.0, 1.0), Vec3::new(0.0, 0.0, -1.0),
                BundleConnectionError::ReceiverNotReached(PathTermination::Escaped)),
            (Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 0.0, 1.0),
                BundleConnectionError::BackFacingReceiver),
        ] {
            let catalog = scene.discover_bundle_images_from_surface_with_workspace(cx, launch(),
                BundleSurfaceSource::try_new(2, BoundarySide::Positive).unwrap(), target(normal),
                options(), &[seed(direction)], BundleImageOptions {
                    max_images: 1, ..BundleImageOptions::default()
                }, &mut workspace).unwrap();
            assert!(catalog.images().is_empty());
            assert_eq!(catalog.seed_results()[0].outcome(), &BundleImageSeedOutcome::Miss(expected));
            catalog.validate_proposal_support().unwrap();
        }
    });
}

#[test]
fn g0_unresolved_native_work_never_becomes_empty_light_support() {
    for error in [
        BundleConnectionError::ReceiverNotReached(PathTermination::DistanceLimit),
        BundleConnectionError::ReceiverNotReached(PathTermination::SurfaceLimit),
        BundleConnectionError::ReceiverNotReached(PathTermination::Extinguished),
        BundleConnectionError::BranchChanged,
        BundleConnectionError::BranchPlanNotConsumed { planned: 2, consumed: 1 },
        BundleConnectionError::ContactResolution { residual_m: 1e-4 },
        BundleConnectionError::IterationLimit { residual_m: 1e-4 },
        BundleConnectionError::TraceBudgetExceeded,
        BundleConnectionError::NoDescent { residual_m: 1e-4 },
        BundleConnectionError::Projection(BundleError::SingularProjection),
        BundleConnectionError::Projection(BundleError::GrazingProjection),
    ] {
        let catalog = BundleImageCatalog {
            images: Vec::new(),
            seeds: vec![BundleImageSeedResult { outcome: BundleImageSeedOutcome::Miss(error.clone()),
                work: Work::default() }],
            work: Work::default(),
        };
        assert_eq!(catalog.validate_proposal_support(),
            Err(BundleImageError::Seed { seed: 0, error }));
    }
}
