//! Updating mounted snapshots preserves native scene ownership and SI outputs.

use super::*;
use fs_render::animated_instances::{
    AnimatedGeometryInstance, RigidTransformTrajectory, TransformKeyframe,
};

#[derive(Debug, PartialEq)]
enum PlacementState {
    Static(RigidTransform),
    Animated(Vec<TransformKeyframe>),
}

struct MountedState {
    id: u64,
    geometry_id: fs_blake3::ContentHash,
    geometry: SharedGeometry,
    placement: PlacementState,
    material: Material,
    emission: Option<(LiftedSpectrum, f64)>,
}

fn mounted_state(primitives: &[Primitive]) -> Vec<MountedState> {
    primitives
        .iter()
        .map(|primitive| {
            let (id, geometry_id, geometry, placement) = match &primitive.shape {
                Shape::Instance(body) => (
                    body.object_id(),
                    body.geometry_identity(),
                    body.geometry().clone(),
                    PlacementState::Static(body.transform()),
                ),
                Shape::AnimatedInstance(body) => (
                    body.object_id(),
                    body.geometry_identity(),
                    body.geometry().clone(),
                    PlacementState::Animated(body.trajectory().keyframes().to_vec()),
                ),
                _ => panic!("this fixture uses native instances"),
            };
            MountedState {
                id,
                geometry_id,
                geometry,
                placement,
                material: primitive.material,
                emission: primitive.emission,
            }
        })
        .collect()
}

fn assert_unchanged(primitives: &[Primitive], before: &[MountedState]) {
    let after = mounted_state(primitives);
    assert_eq!(after.len(), before.len());
    for (actual, expected) in after.iter().zip(before) {
        assert_eq!(actual.id, expected.id);
        assert_eq!(actual.geometry_id, expected.geometry_id);
        assert!(actual.geometry.ptr_eq(&expected.geometry));
        assert_eq!(actual.placement, expected.placement);
        assert_eq!(actual.material, expected.material);
        assert_eq!(actual.emission, expected.emission);
    }
}

fn primitive(body: &FixtureGeometry) -> Primitive {
    Primitive {
        shape: Shape::Instance(body.instance.clone()),
        material: body.material,
        emission: None,
    }
}

#[test]
fn reordered_updates_reach_camera_and_probe_transport_with_supplemental_emission() {
    with_cx(|cx, _| {
        let mut authored = vec![
            finite_fixture(41, placement([-4.0, 0.0, 1.0])),
            finite_fixture(42, placement([4.0, 0.0, 1.0])),
        ];
        for item in &mut authored {
            item.emitters.push(FixtureEmitter {
                id: EmitterId(1),
                source: isotropic([0.0, 0.0, 0.4], 20.0, 0.5),
            });
        }
        let previous = assemble(&authored, cx);
        let unrelated = FixtureGeometry {
            instance: GeometryInstance::try_new(
                90001,
                fs_blake3::hash_domain("org.frankensim.test.fixture-update", b"supplemental"),
                SharedGeometry::mesh(square(0.25)),
                placement([0.0, 0.0, 1.0]),
            )
            .unwrap(),
            material: gray(),
        };
        let mut updated = native_scene(vec![primitive(&unrelated)]);
        previous
            .append_geometry_to(&mut updated.primitives, cx)
            .unwrap();
        updated.primitives.swap(1, 3);
        let extra_before = mounted_state(&updated.primitives[..1]);

        authored[0].transform =
            RigidTransform::try_new([0.0, FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2], [4.0, 0.0, 3.0])
                .unwrap();
        authored[0].optics = FixtureOpticalState::try_new(0.6, [0.5; SPECTRAL_BINS]).unwrap();
        authored[0].emitters[0].source = authored[0].emitters[0].source.dimmed(0.5).unwrap();
        authored[0].surfaces[1].source = authored[0].surfaces[1].source.dimmed(0.25).unwrap();
        let body = &mut authored[0].geometry[1];
        body.instance = GeometryInstance::try_new(
            body.instance.object_id(),
            fs_blake3::hash_domain("org.frankensim.test.fixture-update", b"larger-aperture"),
            SharedGeometry::mesh(square(0.3)),
            body.instance.transform(),
        )
        .unwrap();
        body.material = Material::Lambertian {
            reflectance: LiftedSpectrum { c: [0.25; 3] },
        };
        authored[1].transform =
            RigidTransform::try_new([0.0, 0.0, FRAC_1_SQRT_2, FRAC_1_SQRT_2], [-4.0, 0.0, 3.0])
                .unwrap();
        authored[1].optics = FixtureOpticalState::try_new(0.4, [1.0; SPECTRAL_BINS]).unwrap();
        authored.reverse();
        let replacement = assemble(&authored, cx);
        assert_eq!(
            replacement
                .replace_geometry_in(&previous, &mut updated.primitives, cx)
                .unwrap(),
            [1, 4, 3, 2]
        );
        assert_unchanged(&updated.primitives[..1], &extra_before);
        close(replacement.sources()[1].flux_lm(), 6.0 * PI);
        close(replacement.surface_sources()[3].flux_lm(), 0.75 * PI);
        close(previous.surface_sources()[1].flux_lm(), 2.0 * PI);
        assert!(matches!(
            PathLuminanceJob::new_with_rig(
                &updated,
                &previous,
                &[],
                PathMedium::Vacuum,
                TransportOptions::default(),
                settings(1),
                cx,
            ),
            Err(PhotometryError::Fixture(
                FixtureError::SurfaceGeometryMismatch { .. }
            ))
        ));

        let mut fresh = native_scene(vec![primitive(&unrelated)]);
        replacement
            .append_geometry_to(&mut fresh.primitives, cx)
            .unwrap();
        // Match the caller's retained scene order, independently of rig order.
        fresh.primitives.swap(2, 4);
        assert_unchanged(&updated.primitives, &mounted_state(&fresh.primitives));
        let supplemental = [PhotometricAreaEmitter::try_new(
            0,
            9.0 * PI,
            0.5,
            SourceSpectrum::equal_energy(),
            EmissionSidedness::Front,
        )
        .unwrap()];
        let density = [1.0];
        let haze = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-10.0, -3.0, 0.0], [20.0, 6.0, 8.0]),
                MassAbsorptionSpectrum::constant(0.04).unwrap(),
                None,
                AbsorbingGridBudget {
                    max_grid_cells: 16,
                    max_steps_per_segment: 64,
                },
            )
            .unwrap(),
            MassScatteringSpectrum::constant(0.06).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        for index in [0, 2] {
            let frame = replacement.geometry()[index].instance.transform();
            let camera = PhysicalCamera::try_look_at(
                frame.transform_point(Point3::new(0.0, 0.0, 2.0)),
                frame.transform_point(Point3::new(0.0, 0.0, 0.0)),
                frame.transform_vector(Vec3::new(0.0, 1.0, 0.0)),
                CameraProjection::try_vertical_fov(0.03).unwrap(),
                2.0,
                Aperture::try_circular(0.01).unwrap(),
            )
            .unwrap();
            let mut actual = PathLuminanceJob::new_with_rig_camera(
                &updated,
                &camera,
                &replacement,
                &supplemental,
                PathMedium::Scattering(&haze),
                TransportOptions::default(),
                settings(128),
                cx,
            )
            .unwrap();
            let mut expected = PathLuminanceJob::new_with_rig_camera(
                &fresh,
                &camera,
                &replacement,
                &supplemental,
                PathMedium::Scattering(&haze),
                TransportOptions::default(),
                settings(128),
                cx,
            )
            .unwrap();
            actual.advance(cx, usize::MAX).unwrap();
            expected.advance(cx, usize::MAX).unwrap();
            let value = actual.pixel_estimate(0).unwrap();
            assert_eq!(
                path_bits(value),
                path_bits(expected.pixel_estimate(0).unwrap())
            );
            assert!(value.luminance_cd_m2 > 0.0);
            assert!(value.total_error_bound_cd_m2.is_none());
            assert_eq!(
                actual.emitter_identities(),
                Some(replacement.emitter_identities())
            );
            assert_eq!(
                actual.surface_emitter_identities(),
                Some(replacement.surface_emitter_identities())
            );
        }
        let probes = [probe([0.0, 0.0, 6.0], [0.0, 0.0, -1.0])];
        let sampling = IlluminanceSampling {
            samples_per_probe: 128,
            max_depth: 3,
            sampler: Sampler::Iid,
            seed: 173,
        };
        let mut actual = PathIlluminanceJob::new_with_rig(
            &updated,
            &replacement,
            &supplemental,
            PathMedium::Scattering(&haze),
            &probes,
            TransportOptions::default(),
            sampling,
            cx,
        )
        .unwrap();
        let mut expected = PathIlluminanceJob::new_with_rig(
            &fresh,
            &replacement,
            &supplemental,
            PathMedium::Scattering(&haze),
            &probes,
            TransportOptions::default(),
            sampling,
            cx,
        )
        .unwrap();
        actual.advance(cx, usize::MAX).unwrap();
        expected.advance(cx, usize::MAX).unwrap();
        let value = actual.probe_estimate(0).unwrap();
        assert_eq!(Some(value), expected.probe_estimate(0));
        assert!(value.direct_lux > 0.0 && value.continuous_lux > 0.0);
        assert!(value.total_error_bound_lux.is_none());
        assert_eq!(supplemental[0].primitive_index(), 0);
    });
}

#[test]
fn late_mount_refusals_and_cancellation_leave_every_primitive_unchanged() {
    with_cx(|cx, gate| {
        let mut authored = finite_fixture(71, placement([0.0, 0.0, 1.0]));
        // The late target is non-emitting housing: every mounted body must
        // match the prior snapshot, not just those with finite-source bindings.
        authored.surfaces.truncate(1);
        let previous = assemble(&[authored.clone()], cx);
        let mut changed = authored;
        changed.transform = placement([3.0, 0.0, 2.0]);
        let replacement = assemble(&[changed.clone()], cx);
        let last = &previous.geometry()[1];
        for refusal in 0..8 {
            let mut primitives = Vec::new();
            previous.append_geometry_to(&mut primitives, cx).unwrap();
            match refusal {
                0 | 1 | 4 | 6 => {
                    primitives[1].shape = Shape::Instance(
                        GeometryInstance::try_new(
                            if refusal == 6 {
                                9999
                            } else {
                                last.instance.object_id()
                            },
                            if refusal == 4 {
                                fs_blake3::hash_domain(
                                    "org.frankensim.test.fixture-update",
                                    b"foreign",
                                )
                            } else {
                                last.instance.geometry_identity()
                            },
                            if refusal == 1 {
                                SharedGeometry::mesh(square(0.25))
                            } else {
                                last.instance.geometry().clone()
                            },
                            if refusal == 0 {
                                placement([9.0, 0.0, 1.0])
                            } else {
                                last.instance.transform()
                            },
                        )
                        .unwrap(),
                    );
                }
                2 => {
                    primitives[1].material = Material::Lambertian {
                        reflectance: LiftedSpectrum { c: [0.2; 3] },
                    }
                }
                3 => primitives[1].emission = Some((LiftedSpectrum { c: [1.0; 3] }, 2.0)),
                5 => {
                    let trajectory = RigidTransformTrajectory::try_new(vec![
                        TransformKeyframe::try_new(0.0, last.instance.transform(), [0.0; 3])
                            .unwrap(),
                    ])
                    .unwrap();
                    primitives[1].shape = Shape::AnimatedInstance(
                        AnimatedGeometryInstance::try_new(
                            last.instance.object_id(),
                            last.instance.geometry_identity(),
                            last.instance.geometry().clone(),
                            trajectory,
                        )
                        .unwrap(),
                    );
                }
                7 => primitives.push(primitive(last)),
                _ => unreachable!(),
            }
            let before = mounted_state(&primitives);
            let error = replacement
                .replace_geometry_in(&previous, &mut primitives, cx)
                .unwrap_err();
            match refusal {
                6 => assert!(matches!(
                    error,
                    FixtureError::GeometryNotMounted {
                        fixture: FixtureId(71),
                        object_id: 711
                    }
                )),
                7 => assert!(matches!(error, FixtureError::DuplicateGeometryId(711))),
                _ => assert!(matches!(
                    error,
                    FixtureError::GeometrySnapshotMismatch {
                        fixture: FixtureId(71),
                        object_id: 711
                    }
                )),
            }
            assert_unchanged(&primitives, &before);
        }
        let body = &mut changed.geometry[1];
        body.instance = GeometryInstance::try_new(
            712,
            body.instance.geometry_identity(),
            body.instance.geometry().clone(),
            body.instance.transform(),
        )
        .unwrap();
        let different_set = assemble(&[changed], cx);
        let mut primitives = Vec::new();
        previous.append_geometry_to(&mut primitives, cx).unwrap();
        let before = mounted_state(&primitives);
        assert!(matches!(
            different_set.replace_geometry_in(&previous, &mut primitives, cx),
            Err(FixtureError::GeometrySetMismatch {
                previous: Some(711),
                replacement: Some(712)
            })
        ));
        assert_unchanged(&primitives, &before);
        gate.request();
        assert!(matches!(
            replacement.replace_geometry_in(&previous, &mut primitives, cx),
            Err(FixtureError::Cancelled)
        ));
        assert_unchanged(&primitives, &before);
        with_cx(|resumed, _| {
            assert_eq!(
                replacement
                    .replace_geometry_in(&previous, &mut primitives, resumed)
                    .unwrap(),
                [0, 1]
            );
            let committed = mounted_state(&primitives);
            assert!(matches!(
                replacement.replace_geometry_in(&previous, &mut primitives, resumed),
                Err(FixtureError::GeometrySnapshotMismatch {
                    fixture: FixtureId(71),
                    object_id: 710
                })
            ));
            assert_unchanged(&primitives, &committed);
        });
    });
}

#[test]
fn three_hundred_fixture_updates_roundtrip_without_reindexing_or_flux_accumulation() {
    with_cx(|cx, _| {
        let mut fixtures: Vec<_> = (0..300)
            .map(|i| {
                finite_fixture(
                    1000 + i,
                    placement([(i % 20) as f64 * 2.0, (i / 20) as f64 * 2.0, 1.0]),
                )
            })
            .collect();
        let previous = assemble(&fixtures, cx);
        let mut primitives = Vec::new();
        previous.append_geometry_to(&mut primitives, cx).unwrap();
        primitives.reverse();
        let before = mounted_state(&primitives);
        for (index, fixture) in fixtures.iter_mut().enumerate() {
            fixture.transform = RigidTransform::try_new(
                [0.0, 0.0, FRAC_1_SQRT_2, FRAC_1_SQRT_2],
                [(index % 20) as f64 * 2.0, (index / 20) as f64 * 2.0, 3.0],
            )
            .unwrap();
            fixture.optics = FixtureOpticalState::try_new(
                [0.0, 0.25, 0.5, 1.0][index % 4],
                [0.5; SPECTRAL_BINS],
            )
            .unwrap();
        }
        fixtures.reverse();
        let replacement = assemble(&fixtures, cx);
        let slots = replacement
            .replace_geometry_in(&previous, &mut primitives, cx)
            .unwrap();
        assert_eq!(slots, (0..600).map(|index| index ^ 1).collect::<Vec<_>>());
        assert_eq!(primitives.len(), 600);
        assert_eq!(replacement.surface_emitter_identities().len(), 600);
        for (body, &slot) in replacement.geometry().iter().zip(&slots) {
            let Shape::Instance(mounted) = &primitives[slot].shape else {
                panic!("native instance")
            };
            assert_eq!(mounted.frame_identity(), body.instance.frame_identity());
            assert!(mounted.geometry().ptr_eq(body.instance.geometry()));
        }
        close(
            replacement
                .fixtures()
                .iter()
                .map(|fixture| fixture.flux_lm)
                .sum(),
            3937.5 * PI,
        );
        close(
            previous
                .fixtures()
                .iter()
                .map(|fixture| fixture.flux_lm)
                .sum(),
            3600.0 * PI,
        );
        assert_eq!(
            previous
                .replace_geometry_in(&replacement, &mut primitives, cx)
                .unwrap(),
            (0..600).rev().collect::<Vec<_>>()
        );
        assert_unchanged(&primitives, &before);
    });
}

#[test]
fn empty_rigs_accept_an_unchanged_native_scene() {
    with_cx(|cx, _| {
        let previous = assemble(&[], cx);
        let replacement = assemble(
            &[fixture(
                17,
                placement([0.0, 0.0, 4.0]),
                vec![isotropic([0.0; 3], 5.0, 0.5)],
            )],
            cx,
        );
        let unrelated = assemble(&[finite_fixture(18, placement([4.0, 0.0, 1.0]))], cx);
        let mut primitives = Vec::new();
        unrelated.append_geometry_to(&mut primitives, cx).unwrap();
        let before = mounted_state(&primitives);
        assert!(
            replacement
                .replace_geometry_in(&previous, &mut primitives, cx)
                .unwrap()
                .is_empty()
        );
        assert_unchanged(&primitives, &before);
    });
}
