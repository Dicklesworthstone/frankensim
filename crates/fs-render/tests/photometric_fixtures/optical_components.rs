//! Shared native optics retain per-emitter calibration and physical transport.
//! Large-rig coverage asserts functionality, not throughput or total error bounds.
use super::*;
use fs_render::photometry::area::{ParaxialBeamSplitter, ParaxialBoxDiffuser, ProjectorBeam};
use fs_render::photometry::fixture::FixtureOpticOverride::{Remove, Replace, Retain};

fn splitter(transmission: f64) -> ParaxialBeamSplitter {
    ParaxialBeamSplitter::try_new(vec![
        ProjectorBeam::try_new([0.0, 0.0], transmission * 0.5).unwrap(),
        ProjectorBeam::try_new([0.02, -0.01], transmission * 0.5).unwrap(),
    ])
    .unwrap()
}

fn authored(id: u64, transform: RigidTransform) -> PhotometricFixture {
    let mut input = finite_fixture(id, transform);
    for (index, emitter) in input.surfaces.iter_mut().enumerate() {
        let s = 0.1 + index as f64 * 0.02;
        let angle = index as f64 * 0.1;
        let optics = ThinLensProjector::try_new(
            RigidTransform::try_new([0.0, 0.0, angle.sin(), angle.cos()], [0.0; 3]).unwrap(),
            s,
            s * 3.0 / (3.0 + s),
            GoboTransmission::try_new(0.04, 0.03, 2, 1, vec![1.0 - index as f64 * 0.25, 0.25])
                .unwrap(),
        )
        .unwrap()
        .with_beam_splitter(splitter(0.75))
        .unwrap()
        .with_diffuser(ParaxialBoxDiffuser::try_new([0.01, 0.02], 0.5, 0.8).unwrap())
        .unwrap();
        emitter.source = emitter
            .source
            .clone()
            .with_projector(optics)
            .unwrap()
            .filtered([0.5 / (index + 1) as f64; SPECTRAL_BINS])
            .unwrap();
    }
    input
}

#[test]
fn g3_shared_components_and_lens_scales_match_native_camera_and_probe_paths() {
    with_cx(|cx, _| {
        let turn =
            RigidTransform::try_new([0.0, FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2], [2.0, -1.0, 3.0])
                .unwrap();
        let mut input = authored(404, turn);
        let gobo = GoboTransmission::try_new(0.12, 0.08, 2, 2, vec![1.0, 0.0, 0.5, 0.5]).unwrap();
        let split = splitter(0.5);
        let frost = ParaxialBoxDiffuser::try_new([0.015, 0.025], 0.6, 0.75).unwrap();
        input.optics = input
            .optics
            .with_lens_scale(1.5, 1.5)
            .unwrap()
            .with_gobo(gobo.clone())
            .with_beam_splitter(Replace(split.clone()))
            .with_diffuser(Replace(frost));
        let rig = assemble(&[input.clone()], cx);
        let mut native = receiver_scene();
        rig.append_geometry_to(&mut native.primitives, cx).unwrap();
        let mut explicit = Vec::new();
        for (index, emitter) in input.surfaces.iter().enumerate() {
            let original = emitter.source.projector().unwrap();
            // Author the final native state explicitly from clear calibration,
            // without the fixture's shared-control selection or assembly path.
            let optics = ThinLensProjector::try_new(
                original.frame(),
                original.gobo_distance_m() * 1.5,
                original.focal_length_m() * 1.5,
                gobo.clone(),
            )
            .unwrap()
            .with_beam_splitter(split.clone())
            .unwrap()
            .with_diffuser(frost)
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
            let assembled = &rig.surface_sources()[index];
            close(
                assembled.clear_flux_lm(),
                emitter.source.clear_flux_lm() * 0.5 * 0.4,
            );
            close(
                assembled.flux_lm(),
                assembled.clear_flux_lm() * 0.5 * 0.5 * 0.75,
            );
            assert_eq!(assembled.spectrum(), emitter.source.spectrum());
            assert_eq!(assembled.projector().unwrap().frame(), original.frame());
            assert_eq!(assembled.projector().unwrap().gobo(), &gobo);
            close(
                assembled.projector().unwrap().focus_distance_m().unwrap(),
                original.focus_distance_m().unwrap() * 1.5,
            );
            explicit.push(source);
        }
        assert_ne!(
            rig.surface_sources()[0].spectrum(),
            rig.surface_sources()[1].spectrum()
        );
        close(
            rig.fixtures()[0].flux_lm,
            explicit.iter().map(PhotometricAreaEmitter::flux_lm).sum(),
        );
        let frame = rig.geometry()[0].instance.transform();
        let position = frame.transform_point(Point3::new(0.1, 0.0, 4.5));
        let camera = PhysicalCamera::try_look_at(
            position,
            frame.transform_point(Point3::new(0.0, 0.0, 0.0)),
            frame.transform_vector(Vec3::new(0.0, 1.0, 0.0)),
            CameraProjection::try_vertical_fov(1.0e-5).unwrap(),
            4.5,
            Aperture::try_circular(0.0).unwrap(),
        )
        .unwrap();
        let policy = Settings {
            max_depth: 1,
            sampler: Sampler::OwenSobol,
            ..settings(256)
        };
        let mut grouped = PathLuminanceJob::new_with_rig_camera(
            &native,
            &camera,
            &rig,
            &[],
            PathMedium::Vacuum,
            TransportOptions::default(),
            policy,
            cx,
        )
        .unwrap();
        let mut flat = PathLuminanceJob::new_with_camera(
            &native,
            &camera,
            PhotometricEmitters {
                points: &[],
                surfaces: &explicit,
            },
            PathMedium::Vacuum,
            TransportOptions::default(),
            policy,
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
        assert_eq!(
            grouped.surface_emitter_identities(),
            Some(rig.surface_emitter_identities())
        );
        let probes = [IlluminanceProbe::try_new(
            input.transform.transform_point(Point3::new(0.0, 0.0, 4.5)),
            frame.transform_vector(Vec3::new(0.0, 0.0, -1.0)),
        )
        .unwrap()];
        let sampling = IlluminanceSampling {
            samples_per_probe: 256,
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
            PhotometricEmitters {
                points: &[],
                surfaces: &explicit,
            },
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
    });
}

#[test]
fn g0_retain_remove_and_replace_use_clear_calibration_without_restoring_filter_losses() {
    with_cx(|cx, _| {
        let input = authored(405, RigidTransform::identity());
        let original = assemble(&[input.clone()], cx);
        let mut retained = input.clone();
        retained.optics = retained
            .optics
            .with_beam_splitter(Retain)
            .with_diffuser(Retain);
        let retained = assemble(&[retained], cx);
        let mut removed = input.clone();
        removed.optics = removed
            .optics
            .with_beam_splitter(Remove)
            .with_diffuser(Remove);
        let removed = assemble(&[removed], cx);
        let mut replaced = input.clone();
        replaced.optics = replaced
            .optics
            .with_gobo(GoboTransmission::try_new(0.04, 0.03, 1, 1, vec![1.0]).unwrap())
            .with_beam_splitter(Remove)
            .with_diffuser(Remove);
        let replaced = assemble(&[replaced], cx);
        for index in 0..input.surfaces.len() {
            let original = &original.surface_sources()[index];
            let kept = &retained.surface_sources()[index];
            let bare = &removed.surface_sources()[index];
            let open = &replaced.surface_sources()[index];
            assert_eq!(original.flux_lm().to_bits(), kept.flux_lm().to_bits());
            assert_eq!(
                original.projector().unwrap().gobo(),
                kept.projector().unwrap().gobo()
            );
            assert!(bare.projector().unwrap().beam_splitter().is_none());
            assert!(bare.projector().unwrap().diffuser().is_none());
            let before = &input.surfaces[index].source;
            close(
                bare.flux_lm(),
                before.clear_flux_lm()
                    * 0.5
                    * 0.4
                    * before.projector().unwrap().gobo().mean_transmission(),
            );
            close(open.flux_lm(), before.clear_flux_lm() * 0.5 * 0.4);
            for source in [kept, bare, open] {
                assert_eq!(
                    source.clear_flux_lm().to_bits(),
                    original.clear_flux_lm().to_bits()
                );
                assert_eq!(source.spectrum(), before.spectrum());
                assert_eq!(
                    source.projector().unwrap().frame(),
                    before.projector().unwrap().frame()
                );
            }
            assert!(before.projector().unwrap().beam_splitter().is_some());
            assert!(before.projector().unwrap().diffuser().is_some());
            assert!(open.flux_lm() > original.flux_lm());
        }
    });
}

#[test]
fn g4_shared_optics_refuse_dark_unsupported_sources_and_preserve_cancelled_paths() {
    with_cx(|cx, gate| {
        for optics in [
            FixtureOpticalState::default()
                .with_gobo(GoboTransmission::try_new(0.02, 0.02, 1, 1, vec![0.0]).unwrap()),
            FixtureOpticalState::default().with_beam_splitter(Remove),
            FixtureOpticalState::default().with_diffuser(Remove),
        ] {
            let mut point = fixture(
                505,
                RigidTransform::identity(),
                vec![isotropic([0.0; 3], 1.0, 0.0)],
            );
            point.optics = optics.clone();
            assert!(matches!(
                PhotometricRig::try_new(&[point], FixtureAssemblyLimits::default(), cx),
                Err(FixtureError::UnsupportedProjectorControl {
                    identity: FixtureEmitterIdentity {
                        fixture: FixtureId(505),
                        emitter: EmitterId(1)
                    }
                })
            ));
            let mut area = finite_fixture(606, RigidTransform::identity());
            area.optics = optics;
            area.surfaces[0].source = area.surfaces[0].source.dimmed(0.0).unwrap();
            assert!(matches!(
                PhotometricRig::try_new(&[area], FixtureAssemblyLimits::default(), cx),
                Err(FixtureError::UnsupportedProjectorControl {
                    identity: FixtureEmitterIdentity {
                        fixture: FixtureId(606),
                        emitter: EmitterId(11)
                    }
                })
            ));
        }
        let mut input = authored(707, placement([0.0, 0.0, 1.0]));
        input.optics = input
            .optics
            .with_diffuser(Remove)
            .with_beam_splitter(Replace(splitter(0.5)));
        let rig = assemble(&[input.clone()], cx);
        let mut native = receiver_scene();
        rig.append_geometry_to(&mut native.primitives, cx).unwrap();
        native.camera.eye = Point3::new(-0.75, 0.0, 4.0);
        let mut full = paths(&rig, &native, PathMedium::Vacuum, 32, cx);
        full.advance(cx, usize::MAX).unwrap();
        assert!(full.pixel_estimate(0).unwrap().luminance_cd_m2 > 0.0);
        let mut split = paths(&rig, &native, PathMedium::Vacuum, 32, cx);
        split.advance(cx, 5).unwrap();
        let prefix = path_bits(split.pixel_estimate(0).unwrap());
        gate.request();
        assert!(matches!(
            split.advance(cx, 4),
            Err(PhotometryError::Cancelled)
        ));
        assert_eq!(path_bits(split.pixel_estimate(0).unwrap()), prefix);
        assert!(matches!(
            PhotometricRig::try_new(&[input.clone()], FixtureAssemblyLimits::default(), cx),
            Err(FixtureError::Cancelled)
        ));
        with_cx(|resumed, _| {
            while !split.is_complete() {
                split.advance(resumed, 7).unwrap();
            }
            let replay = assemble(&[input], resumed);
            for (actual, expected) in replay.surface_sources().iter().zip(rig.surface_sources()) {
                assert_eq!(actual.flux_lm().to_bits(), expected.flux_lm().to_bits());
            }
        });
        assert_eq!(
            path_bits(split.pixel_estimate(0).unwrap()),
            path_bits(full.pixel_estimate(0).unwrap())
        );
        assert_eq!(
            split.surface_emitter_identities(),
            Some(rig.surface_emitter_identities())
        );
    });
}

#[test]
fn g3_three_hundred_pixel_fixtures_share_optics_and_reach_native_illuminance() {
    with_cx(|cx, _| {
        let mesh = SharedGeometry::mesh(square(0.0005));
        let identity =
            fs_blake3::hash_domain("org.frankensim.test.fixture-aperture", b"millimetre-square");
        let original = ThinLensProjector::try_new(
            RigidTransform::identity(),
            0.1,
            0.05,
            GoboTransmission::try_new(0.04, 0.04, 1, 1, vec![0.25]).unwrap(),
        )
        .unwrap();
        let gobo = GoboTransmission::try_new(0.08, 0.08, 32, 32, vec![0.5; 1024]).unwrap();
        let frost = ParaxialBoxDiffuser::try_new([0.01, 0.02], 0.5, 0.5).unwrap();
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
                .with_beam_splitter(Replace(splitter(0.75)))
                .with_diffuser(Replace(frost));
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
        assert_eq!(rig.surface_emitter_identities().len(), 7200);
        assert_eq!(rig.geometry().len(), 7200);
        for (index, fixture) in rig.fixtures().iter().enumerate() {
            assert_eq!(fixture.surfaces, index * 24..(index + 1) * 24);
            close(fixture.flux_lm, if index % 2 == 0 { 13.5 } else { 27.0 });
            assert_eq!(rig.surface_sources()[index * 24].flux_lm(), 0.0);
            assert_eq!(
                rig.surface_emitter_identities()[index * 24],
                FixtureEmitterIdentity {
                    fixture: FixtureId(index as u64 + 1),
                    emitter: EmitterId(1)
                }
            );
        }
        close(
            rig.fixtures().iter().map(|fixture| fixture.flux_lm).sum(),
            6075.0,
        );
        let mut native = native_scene(Vec::new());
        rig.append_geometry_to(&mut native.primitives, cx).unwrap();
        let probes = [probe([0.0, 0.0, 100.0], [0.0, 0.0, -1.0])];
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
        assert!(value.illuminance_lux.is_finite() && value.illuminance_lux > 0.0);
        assert_eq!(value.direct_lux, 0.0);
        assert!(value.total_error_bound_lux.is_none());
        assert_eq!(
            job.surface_emitter_identities(),
            Some(rig.surface_emitter_identities())
        );
    });
}
