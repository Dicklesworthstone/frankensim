//! Fixture iris/shutter preparation through the existing native light field.
//! The cell-average mask and large rig do not certify physical error or speed.
use super::*;
use fs_render::photometry::area::{
    FramingShutter, GoboMask, GoboMaskError, ParaxialBeamSplitter, ParaxialBoxDiffuser,
    PolygonalIris, ProjectorBeam,
};

fn projector(gobo: GoboTransmission, s: f64, roll: f64) -> ThinLensProjector {
    let frame =
        RigidTransform::try_new([0.0, 0.0, (roll * 0.5).sin(), (roll * 0.5).cos()], [0.0; 3])
            .unwrap();
    ThinLensProjector::try_new(frame, s, s * 3.0 / (3.0 + s), gobo)
        .unwrap()
        .with_beam_splitter(
            ParaxialBeamSplitter::try_new(vec![
                ProjectorBeam::try_new([0.0; 2], 0.5).unwrap(),
                ProjectorBeam::try_new([0.02, -0.01], 0.25).unwrap(),
            ])
            .unwrap(),
        )
        .unwrap()
        .with_diffuser(ParaxialBoxDiffuser::try_new([0.01, 0.02], 0.5, 0.5).unwrap())
        .unwrap()
}

fn authored(id: u64, transform: RigidTransform) -> PhotometricFixture {
    let mut input = finite_fixture(id, transform);
    for (index, emitter) in input.surfaces.iter_mut().enumerate() {
        let values = if index == 0 {
            vec![1.0, 0.25, 0.75, 0.5]
        } else {
            vec![0.5, 1.0, 0.125, 0.875]
        };
        let gobo =
            GoboTransmission::try_new(0.12 - index as f64 * 0.02, 0.08, 2, 2, values).unwrap();
        emitter.source = emitter
            .source
            .clone()
            .with_projector(projector(
                gobo,
                0.1 + index as f64 * 0.02,
                index as f64 * 0.2,
            ))
            .unwrap()
            .filtered([0.75; SPECTRAL_BINS])
            .unwrap();
    }
    input
}

#[test]
fn g3_transformed_iris_and_shutters_retain_individual_optics_in_camera_and_probe_paths() {
    with_cx(|cx, _| {
        let turn =
            RigidTransform::try_new([0.0, FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2], [2.0, -1.0, 3.0])
                .unwrap();
        let mut input = authored(811, turn);
        let mask = GoboMask::try_new(
            Some(PolygonalIris::try_new(0.035, 8, 0.2).unwrap()),
            vec![
                FramingShutter::try_new(0.4, 0.034).unwrap(),
                FramingShutter::try_new(PI + 0.25, 0.032).unwrap(),
            ],
        )
        .unwrap();
        input.optics = input
            .optics
            .with_lens_scale(1.5, 1.5)
            .unwrap()
            .with_gobo_mask(mask.clone(), [3, 2])
            .unwrap();
        let rig = assemble(&[input.clone()], cx);
        let mut native = receiver_scene();
        rig.append_geometry_to(&mut native.primitives, cx).unwrap();
        let mut explicit = Vec::new();
        for (index, emitter) in input.surfaces.iter().enumerate() {
            let original = emitter.source.projector().unwrap();
            // Prepare each original map independently of fixture cache/selection.
            let gobo = original.gobo().masked(&mask, [3, 2], cx).unwrap();
            assert_eq!((gobo.columns(), gobo.rows()), (6, 4));
            assert_eq!(gobo.width_m(), original.gobo().width_m());
            assert_eq!(gobo.height_m(), original.gobo().height_m());
            assert!(gobo.mean_transmission() > 0.0);
            assert!(gobo.mean_transmission() < original.gobo().mean_transmission());
            let optics = ThinLensProjector::try_new(
                original.frame(),
                original.gobo_distance_m() * 1.5,
                original.focal_length_m() * 1.5,
                gobo.clone(),
            )
            .unwrap()
            .with_beam_splitter(original.beam_splitter().unwrap().clone())
            .unwrap()
            .with_diffuser(*original.diffuser().unwrap())
            .unwrap();
            let source = PhotometricAreaEmitter::try_new(
                index + 1,
                emitter.source.clear_flux_lm(),
                1.0,
                emitter.source.spectrum().clone(),
                emitter.source.sidedness(),
            )
            .unwrap()
            .with_projector(optics)
            .unwrap()
            .dimmed(0.5)
            .unwrap()
            .filtered([0.4; SPECTRAL_BINS])
            .unwrap();
            let actual = &rig.surface_sources()[index];
            close(
                actual.flux_lm(),
                emitter.source.clear_flux_lm() * 0.5 * 0.4 * gobo.mean_transmission() * 0.75 * 0.5,
            );
            assert_eq!(actual.projector().unwrap().gobo(), &gobo);
            assert_eq!(actual.spectrum(), emitter.source.spectrum());
            assert_eq!(actual.projector().unwrap().frame(), original.frame());
            explicit.push(source);
        }
        assert_ne!(
            rig.surface_sources()[0].projector().unwrap().gobo(),
            rig.surface_sources()[1].projector().unwrap().gobo()
        );
        let frame = rig.geometry()[0].instance.transform();
        native.camera.eye = frame.transform_point(Point3::new(0.1, 0.0, 4.5));
        native.camera.forward = frame.transform_vector(Vec3::new(0.0, 0.0, -1.0));
        native.camera.up = frame.transform_vector(Vec3::new(0.0, 1.0, 0.0));
        let bindings = PhotometricEmitters {
            points: &[],
            surfaces: &explicit,
        };
        let mut grouped = paths(&rig, &native, PathMedium::Vacuum, 128, cx);
        let mut flat = PathLuminanceJob::new_with_emitters(
            &native,
            bindings,
            PathMedium::Vacuum,
            TransportOptions::default(),
            settings(128),
            cx,
        )
        .unwrap();
        grouped.advance(cx, usize::MAX).unwrap();
        flat.advance(cx, usize::MAX).unwrap();
        let actual = grouped.pixel_estimate(0).unwrap();
        assert!(actual.luminance_cd_m2 > 0.0);
        assert!(actual.total_error_bound_cd_m2.is_none());
        for (actual, expected) in actual
            .xyz
            .into_iter()
            .zip(flat.pixel_estimate(0).unwrap().xyz)
        {
            close(actual, expected);
        }
        let probes = [IlluminanceProbe::try_new(
            input.transform.transform_point(Point3::new(0.0, 0.0, 4.5)),
            native.camera.forward,
        )
        .unwrap()];
        let sampling = IlluminanceSampling {
            samples_per_probe: 128,
            max_depth: 1,
            sampler: Sampler::OwenSobol,
            seed: 173,
        };
        let mut grouped = PathIlluminanceJob::new_with_rig(
            &native,
            &rig,
            &[],
            PathMedium::Vacuum,
            &probes,
            TransportOptions::default(),
            sampling,
            cx,
        )
        .unwrap();
        let mut flat = PathIlluminanceJob::new_with_emitters(
            &native,
            bindings,
            PathMedium::Vacuum,
            &probes,
            TransportOptions::default(),
            sampling,
            cx,
        )
        .unwrap();
        grouped.advance(cx, usize::MAX).unwrap();
        flat.advance(cx, usize::MAX).unwrap();
        let actual = grouped.probe_estimate(0).unwrap();
        assert!(actual.continuous_lux > 0.0);
        close(
            actual.illuminance_lux,
            flat.probe_estimate(0).unwrap().illuminance_lux,
        );
        assert!(actual.total_error_bound_lux.is_none());
        assert_eq!(
            grouped.surface_emitter_identities(),
            Some(rig.surface_emitter_identities())
        );
    });
}

#[test]
fn g4_closed_masks_retain_identities_refuse_unsupported_sources_and_cancel_atomically() {
    with_cx(|cx, gate| {
        let closed = GoboMask::try_new(
            Some(PolygonalIris::try_new(0.0, 8, 0.0).unwrap()),
            Vec::new(),
        )
        .unwrap();
        let mut input = authored(812, placement([0.0, 0.0, 1.0]));
        input.optics = input.optics.with_gobo_mask(closed.clone(), [2, 2]).unwrap();
        let rig = assemble(&[input.clone()], cx);
        assert_eq!(rig.fixtures()[0].flux_lm, 0.0);
        assert_eq!(
            rig.surface_emitter_identities(),
            &[
                FixtureEmitterIdentity {
                    fixture: FixtureId(812),
                    emitter: EmitterId(11)
                },
                FixtureEmitterIdentity {
                    fixture: FixtureId(812),
                    emitter: EmitterId(12)
                },
            ]
        );
        for source in rig.surface_sources() {
            assert_eq!(source.flux_lm(), 0.0);
            assert!(source.clear_flux_lm() > 0.0);
        }
        let mut native = receiver_scene();
        rig.append_geometry_to(&mut native.primitives, cx).unwrap();
        native.camera.eye = Point3::new(-0.75, 0.0, 4.0);
        let mut image = paths(&rig, &native, PathMedium::Vacuum, 4, cx);
        image.advance(cx, usize::MAX).unwrap();
        assert_eq!(image.pixel_estimate(0).unwrap().luminance_cd_m2, 0.0);
        for mut unsupported in [
            fixture(
                813,
                RigidTransform::identity(),
                vec![isotropic([0.0; 3], 1.0, 0.0)],
            ),
            finite_fixture(814, RigidTransform::identity()),
        ] {
            let emitter = if unsupported.emitters.is_empty() {
                EmitterId(11)
            } else {
                EmitterId(1)
            };
            let expected = FixtureEmitterIdentity {
                fixture: unsupported.id,
                emitter,
            };
            unsupported.optics = FixtureOpticalState::try_new(0.0, [1.0; SPECTRAL_BINS])
                .unwrap()
                .with_gobo_mask(closed.clone(), [1, 1])
                .unwrap();
            assert!(
                matches!(PhotometricRig::try_new(&[unsupported], FixtureAssemblyLimits::default(), cx),
                Err(FixtureError::UnsupportedProjectorControl { identity }) if identity == expected)
            );
        }
        gate.request();
        let original = input.surfaces[0].source.projector().unwrap().gobo();
        assert!(matches!(
            original.masked(&closed, [2, 2], cx),
            Err(GoboMaskError::Cancelled)
        ));
        assert!(matches!(
            PhotometricRig::try_new(&[input.clone()], FixtureAssemblyLimits::default(), cx),
            Err(FixtureError::Cancelled)
        ));
        assert!(original.mean_transmission() > 0.0);
        assert!(input.surfaces[0].source.flux_lm() > 0.0);
        with_cx(|resumed, _| {
            let replay = assemble(&[input], resumed);
            assert_eq!(
                replay.surface_emitter_identities(),
                rig.surface_emitter_identities()
            );
            for (actual, expected) in replay.surface_sources().iter().zip(rig.surface_sources()) {
                assert_eq!(
                    actual.clear_flux_lm().to_bits(),
                    expected.clear_flux_lm().to_bits()
                );
                assert_eq!(actual.flux_lm(), 0.0);
            }
        });
    });
}

#[test]
fn g1_three_hundred_pixel_fixtures_preserve_analytic_mask_flux_and_native_transport() {
    with_cx(|cx, _| {
        let mesh = SharedGeometry::mesh(square(0.0005));
        let identity =
            fs_blake3::hash_domain("org.frankensim.test.fixture-aperture", b"millimetre-square");
        let original = projector(
            GoboTransmission::try_new(0.04, 0.04, 1, 1, vec![0.25]).unwrap(),
            0.1,
            0.0,
        );
        let gobo = GoboTransmission::try_new(0.08, 0.08, 32, 32, vec![0.5; 1024]).unwrap();
        // The square iris passes one quarter of the original rectangle. Its
        // left half passes the shutter, giving 1/8 open area and 1/16 mean.
        let mask = GoboMask::try_new(
            Some(PolygonalIris::try_new(0.02, 4, 0.0).unwrap()),
            vec![FramingShutter::try_new(0.0, 0.0).unwrap()],
        )
        .unwrap();
        let mut fixtures = Vec::new();
        for fixture_index in 0..300_u64 {
            let dimmer = if fixture_index % 2 == 0 { 0.5 } else { 1.0 };
            let mut item = fixture(
                fixture_index + 1,
                placement([fixture_index as f64 * 0.1 - 15.0, 0.0, 0.0]),
                Vec::new(),
            );
            item.optics = FixtureOpticalState::try_new(dimmer, [0.5; SPECTRAL_BINS])
                .unwrap()
                .with_gobo(gobo.clone())
                .with_gobo_mask(mask.clone(), [2, 2])
                .unwrap();
            for emitter_index in 0..24_usize {
                item.geometry.push(FixtureGeometry {
                    instance: GeometryInstance::try_new(
                        fixture_index * 24 + emitter_index as u64 + 1,
                        identity,
                        mesh.clone(),
                        placement([emitter_index as f64 * 0.002, 0.0, 0.0]),
                    )
                    .unwrap(),
                    material: gray(),
                });
                item.surfaces.push(FixtureSurfaceEmitter {
                    id: EmitterId(emitter_index as u64 + 1),
                    source: PhotometricAreaEmitter::try_new(
                        emitter_index,
                        32.0,
                        (emitter_index % 4) as f64 * 0.25,
                        SourceSpectrum::equal_energy(),
                        EmissionSidedness::Front,
                    )
                    .unwrap()
                    .with_projector(original.clone())
                    .unwrap(),
                });
            }
            fixtures.push(item);
        }
        let rig = assemble(&fixtures, cx);
        assert_eq!(rig.surface_sources().len(), 7200);
        for (index, fixture) in rig.fixtures().iter().enumerate() {
            assert_eq!(fixture.surfaces, index * 24..(index + 1) * 24);
            close(fixture.flux_lm, if index % 2 == 0 { 1.6875 } else { 3.375 });
            assert_eq!(rig.surface_sources()[index * 24].flux_lm(), 0.0);
            close(
                rig.surface_sources()[index * 24 + 1]
                    .projector()
                    .unwrap()
                    .gobo()
                    .mean_transmission(),
                0.0625,
            );
        }
        close(
            rig.fixtures().iter().map(|fixture| fixture.flux_lm).sum(),
            759.375,
        );
        let mut native = native_scene(Vec::new());
        rig.append_geometry_to(&mut native.primitives, cx).unwrap();
        let probes = [probe([15.1, 0.0, 200.0], [0.0, 0.0, -1.0])];
        let mut job = PathIlluminanceJob::new_with_rig(
            &native,
            &rig,
            &[],
            PathMedium::Vacuum,
            &probes,
            TransportOptions::default(),
            IlluminanceSampling {
                samples_per_probe: 16,
                max_depth: 1,
                sampler: Sampler::OwenSobol,
                seed: 173,
            },
            cx,
        )
        .unwrap();
        job.advance(cx, usize::MAX).unwrap();
        let value = job.probe_estimate(0).unwrap();
        assert!(value.continuous_lux.is_finite() && value.continuous_lux > 0.0);
        assert_eq!(value.direct_lux, 0.0);
        assert!(value.total_error_bound_lux.is_none());
        assert_eq!(
            job.surface_emitter_identities(),
            Some(rig.surface_emitter_identities())
        );
    });
}
