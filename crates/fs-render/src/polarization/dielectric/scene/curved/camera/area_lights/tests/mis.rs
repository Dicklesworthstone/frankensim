//! Analytic area-to-direction reciprocity, including an oblique receiver.

use super::*;
use crate::spectral::SpectralReflectance;

#[test]
fn g1_reciprocal_area_pdf_cancels_receiver_obliquity_and_common_host_index() {
    with_cx(|cx| {
        let point = Point3::new(0.0, 0.0, 2.0);
        let receiver_normal = Vec3::new(0.6, 0.0, -0.8);
        let geometry = InstanceScene::try_new(vec![
            source(0.0),
            mesh(9, vec![[-0.8, -1.0, 1.4], [0.8, -1.0, 2.6], [0.0, 1.0, 2.0]],
                vec![[0, 1, 2]], RigidTransform::identity()),
        ]).unwrap();
        let reflectance = [SpectralReflectance::constant(0.6).unwrap()];
        let emitter = SurfaceEmitter::new(2,
            super::super::super::RadianceSpectrum::try_constant(4.0).unwrap(),
            EmissionSides::Front);
        for (index, absorption) in [(1.0, 0.0), (1.7, 0.5)] {
            let material = glass(index, absorption);
            let media = [SceneMedium::Homogeneous(&material)];
            let scene = CurvedPolarizedScene::try_new_with_materials(
                cx, &geometry, &media, &[], &[], &reflectance,
                vec![
                    SurfaceBinding { object_id: 2, interaction: Interaction::Opaque },
                    SurfaceBinding { object_id: 9, interaction: Interaction::DepolarizingDiffuse {
                        medium: 0, reflectance: 0,
                    } },
                ],
                100_000,
            ).unwrap();
            let panel = CurvedAreaEmitter::try_new(cx, &scene, &emitter, 1).unwrap();
            let sample = panel.sample_area([4.0 / 9.0, 0.5]).unwrap();
            let target = BundleConnectionTarget::try_diffuse_receiver(
                9, point, receiver_normal, 2e-7,
            ).unwrap();
            let connected = scene.connect_bundle_from_surface_with_policy(
                cx,
                BundleConnectionLaunch {
                    ray: Ray { origin: sample.position_m(), dir: Vec3::new(0.0, 0.0, 1.0) },
                    transverse_axis: Vec3::new(0.0, 1.0, 0.0),
                    wavelength_nm: 550.0,
                    normalized_source_stokes: Stokes::unpolarized(1.0),
                    initial_medium: 0,
                },
                BundleSurfaceSource::try_new(2, BoundarySide::Positive).unwrap(),
                target, options(), BundleBranchPolicy::Sequence(&[]),
            ).unwrap();
            let weights = panel.mis_weights(sample, &connected, target, options().projection).unwrap();
            // The oblique receiving footprint is r²/cos(receiver); reciprocal
            // etendue cancels that cosine, leaving r²/cos(source) = 4.
            let light_pdf = 4.0 / 0.27;
            let continuation_pdf = 0.8 / core::f64::consts::PI;
            let light_squared = light_pdf * light_pdf;
            let continuation_squared = continuation_pdf * continuation_pdf;
            near(weights.light_log_pdf_solid_angle().exp(), light_pdf, 3e-6);
            near(weights.continuation_log_pdf_solid_angle().exp(), continuation_pdf, 3e-6);
            near(weights.light_weight(),
                light_squared / (light_squared + continuation_squared), 3e-6);
            near(weights.continuation_weight(),
                continuation_squared / (light_squared + continuation_squared), 3e-6);
            near(weights.light_weight() + weights.continuation_weight(), 1.0, 1e-15);
            // Reverse lookup keeps the same represented area density and
            // actual integrated launch point, independently of patch numbering.
            let recovered = panel.sample_at_contact(
                connected.path().source_contact().unwrap(), sample.position_m(),
            ).unwrap();
            assert_eq!(recovered.position_m(), sample.position_m());
            assert_eq!(recovered.triangle_index(), sample.triangle_index());
            assert_eq!(recovered.pdf_area_per_m2(), sample.pdf_area_per_m2());
        }
    });
}

#[test]
fn g1_volume_area_pdf_uses_transverse_etendue_not_the_newton_section() {
    with_cx(|cx| {
        let point = Point3::new(0.3, 0.0, 2.0);
        let distance = (0.3_f64 * 0.3 + 4.0).sqrt();
        let direction = Vec3::new(0.3 / distance, 0.0, 2.0 / distance);
        let geometry = InstanceScene::try_new(vec![source(0.0)]).unwrap();
        let emitter = SurfaceEmitter::new(2,
            super::super::super::RadianceSpectrum::try_constant(4.0).unwrap(),
            EmissionSides::Front);
        // An incoming camera ray along +X scatters toward the original source.
        // Rayleigh's scalar PDF is independent of the transported Stokes power.
        let cosine_phase = -direction.x;
        let phase_pdf = 3.0 * (1.0 + cosine_phase * cosine_phase)
            / (16.0 * core::f64::consts::PI);
        let expected_light = distance * distance / (0.27 * direction.z);
        let expected_continuation = phase_pdf * phase_pdf
            / (expected_light * expected_light + phase_pdf * phase_pdf);
        let mut first = None;
        for (index, absorption) in [(1.0, 0.0), (1.7, 0.5)] {
            let material = glass(index, absorption);
            let media = [SceneMedium::Homogeneous(&material)];
            let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media,
                vec![SurfaceBinding { object_id: 2, interaction: Interaction::Opaque }],
                100_000).unwrap();
            let panel = CurvedAreaEmitter::try_new(cx, &scene, &emitter, 1).unwrap();
            let sample = panel.sample_area([4.0 / 9.0, 0.5]).unwrap();
            for normal in [direction.scale(-1.0), Vec3::new(0.6, 0.0, -0.8)] {
                let target = BundleConnectionTarget::try_volume_receiver(0, point, normal, 2e-7)
                    .unwrap();
                let connected = scene.connect_bundle_from_surface_with_policy(
                    cx,
                    BundleConnectionLaunch {
                        ray: Ray { origin: sample.position_m(), dir: direction },
                        transverse_axis: Vec3::new(0.0, 1.0, 0.0),
                        wavelength_nm: 550.0,
                        normalized_source_stokes: Stokes::unpolarized(1.0),
                        initial_medium: 0,
                    },
                    BundleSurfaceSource::try_new(2, BoundarySide::Positive).unwrap(),
                    target, options(), BundleBranchPolicy::Sequence(&[]),
                ).unwrap();
                let weights = panel.volume_mis_weights(sample, &connected, target,
                    options().projection, phase_pdf).unwrap();
                near(weights.light_log_pdf_solid_angle().exp(), expected_light, 3e-6);
                near(weights.continuation_log_pdf_solid_angle().exp(), phase_pdf, 3e-6);
                near(weights.continuation_weight(), expected_continuation, 3e-6);
                near(weights.light_weight() + weights.continuation_weight(), 1.0, 1e-15);
                if let Some(previous) = first {
                    near(weights.continuation_weight(), previous, 3e-6);
                } else {
                    first = Some(weights.continuation_weight());
                }
                let zero = panel.volume_mis_weights(sample, &connected, target,
                    options().projection, 0.0).unwrap();
                assert_eq!(zero.light_weight(), 1.0);
                assert_eq!(zero.continuation_weight(), 0.0);
                assert_eq!(zero.continuation_log_pdf_solid_angle(), f64::NEG_INFINITY);
                for invalid in [-1.0, f64::NAN, f64::INFINITY] {
                    assert!(matches!(panel.volume_mis_weights(sample, &connected, target,
                        options().projection, invalid), Err(AreaEmitterError::NumericalRange)));
                }
                let wrong_medium = BundleConnectionTarget::try_volume_receiver(1, point, normal, 2e-7)
                    .unwrap();
                assert!(matches!(panel.volume_mis_weights(sample, &connected, wrong_medium,
                    options().projection, phase_pdf), Err(AreaEmitterError::SourceContactMismatch)));
                let wrong_point = BundleConnectionTarget::try_volume_receiver(0,
                    Point3::new(point.x + 0.1, point.y, point.z), normal, 2e-7).unwrap();
                assert!(matches!(panel.volume_mis_weights(sample, &connected, wrong_point,
                    options().projection, phase_pdf), Err(AreaEmitterError::SourceContactMismatch)));
            }
        }
    });
}

#[test]
fn g1_rough_curved_receiver_uses_native_area_on_both_material_sides() {
    use crate::dielectric::DielectricSurface;
    use fs_geom::fixtures::SphereChart;

    with_cx(|cx| {
        let outside = glass(1.3, 0.2);
        let inside = glass(1.7, 0.4);
        let media = [SceneMedium::Homogeneous(&outside), SceneMedium::Homogeneous(&inside)];
        let height = (1.0_f64 - 0.3 * 0.3).sqrt();
        // The first source is outside the sphere; the second is at its centre.
        // Neither connection crosses the rough interface: it is the receiver.
        for (source_z, medium, target_z, native_cosine) in [
            (0.0, 0, 3.0 - height, None),
            (3.0, 1, 3.0 + height, Some(1.0)),
        ] {
            let geometry = InstanceScene::try_new(vec![source(source_z),
                GeometryInstance::try_new(9,
                    hash_domain("rough-receiver-mis-sphere", b"native"),
                    SharedGeometry::chart(SphereChart {
                        center: Point3::new(0.0, 0.0, 3.0), radius: 1.0,
                    }), RigidTransform::identity()).unwrap(),
            ]).unwrap();
            let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![
                SurfaceBinding { object_id: 2, interaction: Interaction::Opaque },
                SurfaceBinding { object_id: 9, interaction: Interaction::Dielectric {
                    outside: 0, inside: 1, surface: DielectricSurface::try_rough(0.3).unwrap(),
                } },
            ], 100_000).unwrap();
            let emitter = SurfaceEmitter::new(2,
                super::super::super::RadianceSpectrum::try_constant(4.0).unwrap(),
                EmissionSides::Front);
            let panel = CurvedAreaEmitter::try_new(cx, &scene, &emitter, 1).unwrap();
            let sample = panel.sample_area([4.0 / 9.0, 0.5]).unwrap();
            let point = Point3::new(0.3, 0.0, target_z);
            let displacement = point.delta_from(sample.position_m());
            let distance = displacement.norm();
            let direction = displacement.scale(1.0 / distance);
            let area_density = distance * distance / (AREA_FOR_TEST * direction.z);
            let continuation_pdf = 0.17;
            let expected_hit = continuation_pdf * continuation_pdf
                / (area_density * area_density + continuation_pdf * continuation_pdf);
            for section in [direction.scale(-1.0), Vec3::new(0.6, 0.0, -0.8)] {
                let target = BundleConnectionTarget::try_surface_receiver(9, point, section, 2e-7)
                    .unwrap();
                let connected = scene.connect_bundle_from_surface_with_policy(cx,
                    BundleConnectionLaunch {
                        ray: Ray { origin: sample.position_m(), dir: direction },
                        transverse_axis: Vec3::new(0.0, 1.0, 0.0), wavelength_nm: 550.0,
                        normalized_source_stokes: Stokes::unpolarized(1.0), initial_medium: medium,
                    }, BundleSurfaceSource::try_new(2, BoundarySide::Positive).unwrap(),
                    target, options(), BundleBranchPolicy::Sequence(&[])).unwrap();
                assert!(connected.path().is_surface_receiver());
                assert_eq!(connected.path().medium_index(), medium);
                assert_eq!(connected.path().bundle().reverse_branch_log_probability(), 0.0);
                let contact = connected.path().last_contact().unwrap();
                let cosine = contact.witness.normal.dot(direction).abs();
                if let Some(expected) = native_cosine { near(cosine, expected, 2e-6); }
                let projection = connected.receiver_projection(target, options().projection).unwrap();
                near(projection.area_per_solid_angle_m2(), distance * distance / cosine, 3e-6);
                let weights = panel.surface_mis_weights(sample, &connected, target,
                    options().projection, continuation_pdf).unwrap();
                near(weights.light_log_pdf_solid_angle().exp(), area_density, 3e-6);
                near(weights.continuation_weight(), expected_hit, 3e-6);
                near(weights.light_weight() + weights.continuation_weight(), 1.0, 1e-15);
                let zero = panel.surface_mis_weights(sample, &connected, target,
                    options().projection, 0.0).unwrap();
                assert_eq!((zero.light_weight(), zero.continuation_weight()), (1.0, 0.0));
                let diffuse_target = BundleConnectionTarget::try_diffuse_receiver(
                    9, point, section, 2e-7).unwrap();
                assert!(matches!(panel.surface_mis_weights(sample, &connected, diffuse_target,
                    options().projection, continuation_pdf), Err(AreaEmitterError::SourceContactMismatch)));
            }
        }
    });
}

const AREA_FOR_TEST: f64 = 0.27;

