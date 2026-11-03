//! Surface-launched molecular source survival on the original native walk.

use super::*;
use crate::gradient_index::event::BoundarySide;
use crate::polarization::dielectric::scene::curved::bundle::connection::BundleSurfaceSource;
use crate::polarization::dielectric::scene::curved::bundle::connection::images::{
    BundleImageOptions, BundleImageSeed,
};
use crate::polarization::dielectric::scene::curved::participating::CurvedDensityField;
use crate::polarization::dielectric::scene::curved::participating::transmittance::images::CurvedTransmittanceImageWorkspace;

fn with_surface_source(f: impl FnOnce(&Cx<'_>, &CurvedPolarizedScene<'_>)) {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            plane(1, Point3::new(0.0, 0.0, 0.0)),
            plane(9, Point3::new(0.0, 0.0, 2.0)),
        ])
        .unwrap();
        let material = glass(1.2, 0.17);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Opaque,
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
fn g1_surface_source_images_keep_native_origin_and_one_uniform_survival_replay() {
    with_surface_source(|cx, scene| {
        let spectrum = MassScatteringSpectrum::constant(0.4).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, scene, &volumes).unwrap();
        let source = BundleSurfaceSource::try_new(1, BoundarySide::Positive).unwrap();
        let launch = launch(Vec3::new(0.0, 0.0, 1.0));
        let target = target(Point3::new(0.0, 0.0, 2.0));
        let tracking = CurvedTransmittanceOptions {
            max_free_flight_draws: 0,
            max_tracking_candidates: 0,
        };
        let single = scattering
            .connect_bundle_from_surface_with_transmittance_samples(
                cx,
                launch,
                source,
                target,
                options(),
                BundleBranchPolicy::Sequence(&[]),
                tracking,
                |_| panic!("uniform survival must consume no samples"),
            )
            .unwrap();
        let seeds = [BundleImageSeed {
            initial_direction: launch.ray.dir,
            policy: BundleBranchPolicy::Sequence(&[]),
        }; 2];
        let mut workspace = CurvedTransmittanceImageWorkspace::try_new(1, 2).unwrap();
        let catalog = scattering
            .discover_bundle_images_from_surface_with_transmittance_samples(
                cx,
                launch,
                source,
                target,
                options(),
                &seeds,
                BundleImageOptions {
                    max_images: 1,
                    ..BundleImageOptions::default()
                },
                tracking,
                &mut workspace,
                |_, _| panic!("uniform images must consume no samples"),
            )
            .unwrap();
        assert_eq!(catalog.images().len(), 1);
        let survival = &catalog.images()[0];
        assert_eq!(survival.connection(), single.connection());
        assert_eq!(catalog.free_flight_draws(), 0);
        assert_eq!(catalog.tracking_candidates(), 0);
        assert_eq!(
            catalog.trace_evaluations(),
            catalog.geometric_catalog().trace_evaluations() + 1
        );
        let path = survival.connection().path();
        assert_eq!(path.termination(), PathTermination::Detected);
        assert_eq!(path.surface_events(), 1, "source departure is not a surface event");
        assert_eq!(path.bundle().propagation().initial_ray().origin, launch.ray.origin);
        let propagation = path.bundle().propagation().propagation();
        near(propagation.distance_m, 2.0, 3e-7);
        near(propagation.optical_depth, 0.34, 6e-8);
        near(survival.log_scattering_survival(), -0.8, 2e-7);
        assert!(!survival.ratio_tracking_zero());
    });
}

struct FullDensity;

impl CurvedDensityField for FullDensity {
    fn density_kg_m3(&self, _: Point3) -> Result<f64, VolumeTransportError> {
        Ok(1.0)
    }
}

#[test]
fn g3_zero_ratio_surface_replay_still_resolves_the_original_receiver() {
    with_surface_source(|cx, scene| {
        let spectrum = MassScatteringSpectrum::constant(2.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let density = FullDensity;
        let volumes = [Some(
            CurvedRayleighMedium::try_new_heterogeneous(
                &spectrum,
                &density,
                1.0,
                ContentHash([181; 32]),
                &phase,
            )
            .unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, scene, &volumes).unwrap();
        let source = BundleSurfaceSource::try_new(1, BoundarySide::Positive).unwrap();
        let mut draws = 0;
        let result = scattering
            .connect_bundle_from_surface_with_transmittance_samples(
                cx,
                launch(Vec3::new(0.0, 0.0, 1.0)),
                source,
                target(Point3::new(0.0, 0.0, 2.0)),
                options(),
                BundleBranchPolicy::Sequence(&[]),
                CurvedTransmittanceOptions {
                    max_free_flight_draws: 1,
                    max_tracking_candidates: 1,
                },
                |draw| {
                    assert_eq!(draw, 0);
                    draws += 1;
                    flight_sample(2.0, 0.4)
                },
            )
            .unwrap();
        assert_eq!(draws, 1);
        assert_eq!(result.free_flight_draws(), 1);
        assert_eq!(result.tracking_candidates(), 1);
        assert!(result.ratio_tracking_zero());
        assert_eq!(result.log_scattering_survival(), f64::NEG_INFINITY);
        assert_eq!(result.connection().path().termination(), PathTermination::Detected);
        assert_eq!(result.connection().path().surface_events(), 1);
        assert_eq!(
            result.connection().path().bundle().propagation().initial_ray().origin,
            Point3::new(0.0, 0.0, 0.0)
        );
        near(
            result.connection().path().bundle().propagation().propagation().distance_m,
            2.0,
            3e-7,
        );
        assert_eq!(result.trace_evaluations(), result.connection().trace_evaluations() + 1);
    });
}
