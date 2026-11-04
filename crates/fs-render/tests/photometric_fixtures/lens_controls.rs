//! Shared physical lens controls through the original native projector paths.
use super::*;
use fs_render::photometry::area::{ParaxialBeamSplitter, ParaxialBoxDiffuser, ProjectorBeam};
use fs_render::photometry::source::SourceError;

fn projector(index: usize) -> ThinLensProjector {
    let s = 0.1 + index as f64 * 0.02;
    ThinLensProjector::try_new(
        RigidTransform::try_new(
            [
                0.0,
                0.0,
                (index as f64 * 0.1).sin(),
                (index as f64 * 0.1).cos(),
            ],
            [0.0; 3],
        )
        .unwrap(),
        s,
        s * 3.0 / (3.0 + s),
        GoboTransmission::try_new(0.04, 0.03, 2, 1, vec![1.0 - index as f64 * 0.2, 0.25]).unwrap(),
    )
    .unwrap()
    .with_beam_splitter(
        ParaxialBeamSplitter::try_new(vec![
            ProjectorBeam::try_new([0.0, 0.0], 0.5).unwrap(),
            ProjectorBeam::try_new([0.03, -0.01], 0.25).unwrap(),
        ])
        .unwrap(),
    )
    .unwrap()
    .with_diffuser(ParaxialBoxDiffuser::try_new([0.01, 0.02], 0.5, 0.8).unwrap())
    .unwrap()
}

#[test]
fn shared_lens_scale_preserves_individual_optics_and_matches_native_camera_and_probe() {
    with_cx(|cx, _| {
        let turn =
            RigidTransform::try_new([0.0, FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2], [2.0, -1.0, 3.0])
                .unwrap();
        let mut authored = finite_fixture(404, turn);
        for (index, emitter) in authored.surfaces.iter_mut().enumerate() {
            emitter.source = emitter
                .source
                .clone()
                .with_projector(projector(index))
                .unwrap();
        }
        let original = assemble(&[authored.clone()], cx);
        let identity = assemble(
            &[PhotometricFixture {
                optics: authored.optics.clone().with_lens_scale(1.0, 1.0).unwrap(),
                ..authored.clone()
            }],
            cx,
        );
        authored.optics = authored.optics.with_lens_scale(2.0, 2.0).unwrap();
        assert_eq!(authored.optics.lens_scale(), [2.0; 2]);
        let rig = assemble(&[authored.clone()], cx);
        let mut native = receiver_scene();
        rig.append_geometry_to(&mut native.primitives, cx).unwrap();
        let mut explicit = Vec::new();
        for (index, emitter) in authored.surfaces.iter().enumerate() {
            let before = original.surface_sources()[index].projector().unwrap();
            let after = rig.surface_sources()[index].projector().unwrap();
            assert_eq!(after.frame(), before.frame());
            assert_eq!(after.diffuser(), before.diffuser());
            assert_eq!(
                after.beam_splitter().unwrap().beams(),
                before.beam_splitter().unwrap().beams()
            );
            assert_eq!(
                after.gobo().mean_transmission(),
                before.gobo().mean_transmission()
            );
            close(
                after.focus_distance_m().unwrap(),
                2.0 * before.focus_distance_m().unwrap(),
            );
            for source in [
                &identity.surface_sources()[index],
                &rig.surface_sources()[index],
            ] {
                assert_eq!(
                    source.clear_flux_lm().to_bits(),
                    original.surface_sources()[index].clear_flux_lm().to_bits()
                );
                assert_eq!(
                    source.flux_lm().to_bits(),
                    original.surface_sources()[index].flux_lm().to_bits()
                );
                assert_eq!(
                    source.spectrum(),
                    original.surface_sources()[index].spectrum()
                );
            }
            let p = emitter.source.projector().unwrap();
            // Independently construct the new physical distances and prepared
            // frost, instead of using the fixture or distance-refresh helpers.
            let changed = ThinLensProjector::try_new(
                p.frame(),
                2.0 * p.gobo_distance_m(),
                2.0 * p.focal_length_m(),
                p.gobo().clone(),
            )
            .unwrap()
            .with_beam_splitter(p.beam_splitter().unwrap().clone())
            .unwrap()
            .with_diffuser(*p.diffuser().unwrap())
            .unwrap();
            explicit.push(
                PhotometricAreaEmitter::try_new(
                    index + 1,
                    emitter.source.clear_flux_lm(),
                    0.5,
                    emitter.source.spectrum().clone(),
                    emitter.source.sidedness(),
                )
                .unwrap()
                .filtered([0.4; SPECTRAL_BINS])
                .unwrap()
                .with_projector(changed)
                .unwrap(),
            );
        }
        assert_ne!(
            rig.surface_sources()[0]
                .projector()
                .unwrap()
                .gobo_distance_m(),
            rig.surface_sources()[1]
                .projector()
                .unwrap()
                .gobo_distance_m()
        );
        let source_frame = rig.geometry()[0].instance.transform();
        let camera = PhysicalCamera::try_look_at(
            source_frame.transform_point(Point3::new(0.2, 0.0, 6.0)),
            source_frame.transform_point(Point3::new(0.0, 0.0, 0.0)),
            source_frame.transform_vector(Vec3::new(0.0, 1.0, 0.0)),
            CameraProjection::try_vertical_fov(1.0e-5).unwrap(),
            6.0,
            Aperture::try_circular(0.0).unwrap(),
        )
        .unwrap();
        let policy = Settings {
            max_depth: 1,
            sampler: Sampler::OwenSobol,
            ..settings(1024)
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
        assert!(grouped.pixel_estimate(0).unwrap().luminance_cd_m2 > 0.0);
        for (actual, expected) in grouped
            .pixel_estimate(0)
            .unwrap()
            .xyz
            .into_iter()
            .zip(flat.pixel_estimate(0).unwrap().xyz)
        {
            close(actual, expected);
        }
        let probes = [IlluminanceProbe::try_new(
            source_frame.transform_point(Point3::new(0.2, 0.0, 6.0)),
            source_frame.transform_vector(Vec3::new(0.0, 0.0, -1.0)),
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
        let value = grouped.probe_estimate(0).unwrap();
        assert!(value.continuous_lux > 0.0);
        close(
            value.illuminance_lux,
            flat.probe_estimate(0).unwrap().illuminance_lux,
        );
        assert!(value.total_error_bound_lux.is_none());
    });
}

#[test]
fn unsupported_dark_sources_and_unrepresentable_lens_scales_refuse_with_identity() {
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            FixtureOpticalState::default()
                .with_lens_scale(invalid, 1.0)
                .is_err()
        );
        assert!(
            FixtureOpticalState::default()
                .with_lens_scale(1.0, invalid)
                .is_err()
        );
    }
    with_cx(|cx, _| {
        let mut point = fixture(
            505,
            RigidTransform::identity(),
            vec![isotropic([0.0; 3], 1.0, 0.0)],
        );
        point.optics = point.optics.with_lens_scale(2.0, 2.0).unwrap();
        assert!(matches!(
            PhotometricRig::try_new(&[point], FixtureAssemblyLimits::default(), cx),
            Err(FixtureError::UnsupportedLensControl {
                identity: FixtureEmitterIdentity {
                    fixture: FixtureId(505),
                    emitter: EmitterId(1),
                }
            })
        ));
        let mut lambertian = finite_fixture(606, RigidTransform::identity());
        lambertian.optics = FixtureOpticalState::try_new(0.0, [1.0; SPECTRAL_BINS])
            .unwrap()
            .with_lens_scale(2.0, 2.0)
            .unwrap();
        assert!(matches!(
            PhotometricRig::try_new(&[lambertian], FixtureAssemblyLimits::default(), cx),
            Err(FixtureError::UnsupportedLensControl {
                identity: FixtureEmitterIdentity {
                    fixture: FixtureId(606),
                    emitter: EmitterId(11),
                }
            })
        ));
        let mut valid = finite_fixture(707, RigidTransform::identity());
        for emitter in &mut valid.surfaces {
            let optic = ThinLensProjector::try_new(
                RigidTransform::identity(),
                2.0,
                0.5,
                GoboTransmission::try_new(0.1, 0.1, 1, 1, vec![1.0]).unwrap(),
            )
            .unwrap();
            emitter.source = emitter.source.clone().with_projector(optic).unwrap();
        }
        for (s, f) in [(f64::MAX, 1.0), (1.0, f64::from_bits(1))] {
            let mut input = valid.clone();
            input.optics = input.optics.with_lens_scale(s, f).unwrap();
            assert!(matches!(
                PhotometricRig::try_new(&[input], FixtureAssemblyLimits::default(), cx),
                Err(FixtureError::Source {
                    identity: FixtureEmitterIdentity {
                        fixture: FixtureId(707),
                        emitter: EmitterId(11),
                    },
                    error: SourceError::NumericalRange
                })
            ));
        }
        let mut empty = fixture(808, RigidTransform::identity(), Vec::new());
        empty.optics = empty.optics.with_lens_scale(2.0, 0.5).unwrap();
        assert!(assemble(&[empty], cx).surface_sources().is_empty());
    });
}
