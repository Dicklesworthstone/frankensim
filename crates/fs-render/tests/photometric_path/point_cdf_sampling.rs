//! Streaming point proposals through grouped native spectral transport.
//! These are functional regressions, not timing or total-error qualifications.

use super::*;
use core::f64::consts::FRAC_1_SQRT_2;
use fs_exec::RunId;
use fs_render::camera::{Aperture, CameraProjection, PhysicalCamera};
use fs_render::photometry::IlluminanceProbe;
use fs_render::photometry::fixture::{
    EmitterId, FixtureAssemblyLimits, FixtureEmitter, FixtureEmitterIdentity, FixtureId,
    FixtureOpticalState, PhotometricFixture, PhotometricRig,
};
use fs_render::photometry::path::{
    IlluminanceSampling, MAX_POINT_CONNECTIONS_PER_VERTEX, PathIlluminanceJob, PointFixtureSampling,
};
use fs_render::tracer::RenderExecutionConfig;

fn cdf(connections_per_vertex: u32) -> PointFixtureSampling {
    PointFixtureSampling::ImportanceCdf {
        connections_per_vertex,
    }
}

fn fixture(
    id: u64,
    transform: RigidTransform,
    optics: FixtureOpticalState,
    sources: Vec<PhotometricSource>,
) -> PhotometricFixture {
    PhotometricFixture {
        id: FixtureId(id),
        transform,
        optics,
        emitters: sources
            .into_iter()
            .enumerate()
            .map(|(index, source)| FixtureEmitter {
                id: EmitterId(index as u64 + 1),
                source,
            })
            .collect(),
        geometry: Vec::new(),
        surfaces: Vec::new(),
    }
}

fn assemble(fixtures: &[PhotometricFixture], cx: &Cx<'_>) -> PhotometricRig {
    PhotometricRig::try_new(fixtures, FixtureAssemblyLimits::default(), cx).unwrap()
}

fn image_job<'scene, 'field>(
    scene: &'scene Scene,
    rig: &'scene PhotometricRig,
    medium: PathMedium<'scene, 'field>,
    settings: Settings,
    policy: PointFixtureSampling,
    cx: &Cx<'_>,
) -> PathLuminanceJob<'scene, 'field> {
    PathLuminanceJob::new_with_rig(
        scene,
        rig,
        &[],
        medium,
        TransportOptions::default(),
        settings,
        cx,
    )
    .unwrap()
    .with_point_sampling(policy)
    .unwrap()
}

fn probe_job<'scene, 'field>(
    scene: &'scene Scene,
    rig: &'scene PhotometricRig,
    medium: PathMedium<'scene, 'field>,
    probes: &'scene [IlluminanceProbe],
    settings: Settings,
    policy: PointFixtureSampling,
    cx: &Cx<'_>,
) -> PathIlluminanceJob<'scene, 'field> {
    PathIlluminanceJob::new_with_rig(
        scene,
        rig,
        &[],
        medium,
        probes,
        TransportOptions::default(),
        IlluminanceSampling {
            samples_per_probe: settings.spp,
            max_depth: settings.max_depth,
            sampler: settings.sampler,
            seed: settings.seed,
        },
        cx,
    )
    .unwrap()
    .with_point_sampling(policy)
    .unwrap()
}

type SampleBits = ([u64; 3], u32, Option<u64>);

fn image_bits(job: &PathLuminanceJob<'_, '_>) -> Vec<Option<SampleBits>> {
    job.estimates()
        .map(|value| value.map(estimate_bits))
        .collect()
}

fn probe_bits(job: &PathIlluminanceJob<'_, '_>) -> Vec<Option<SampleBits>> {
    job.estimates()
        .map(|value| {
            value.map(|v| {
                (
                    [v.direct_lux, v.continuous_lux, v.illuminance_lux].map(f64::to_bits),
                    v.samples,
                    v.estimated_standard_error_lux.map(f64::to_bits),
                )
            })
        })
        .collect()
}

fn probes() -> [IlluminanceProbe; 2] {
    [0.0, 0.5].map(|x| {
        IlluminanceProbe::try_new(Point3::new(x, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0)).unwrap()
    })
}

fn haze<'a>(cx: &Cx<'_>, density: &'a [f64; 8]) -> ScatteringGrid<'a> {
    ScatteringGrid::try_new(
        cx,
        AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([2; 3], density, [-6.0, -6.0, 0.0], [6.0, 6.0, 1.0]),
            MassAbsorptionSpectrum::table(&[(380.0, 0.02), (780.0, 0.09)]).unwrap(),
            None,
            medium_budget(),
        )
        .unwrap(),
        MassScatteringSpectrum::constant(0.08).unwrap(),
        PhaseFunction::henyey_greenstein(0.35).unwrap(),
    )
    .unwrap()
}

#[test]
fn cdf_single_candidate_preserves_calibrated_camera_and_probe_paths_through_haze() {
    with_cx(|cx, _| {
        let scene = receiver_scene();
        let spectrum = SourceSpectrum::try_from_bins(std::array::from_fn(|bin| {
            0.2 + 2.0 * bin as f64 / (SPECTRAL_BINS - 1) as f64
        }))
        .unwrap();
        let active = PhotometricSource::try_new(
            RigidTransform::try_new([1.0, 0.0, 0.0, 0.0], [0.25, 0.0, 0.0]).unwrap(),
            100.0,
            0.6,
            AngularDistribution::cosine_power(1.0, FRAC_PI_2).unwrap(),
            spectrum,
        )
        .unwrap();
        let away = PhotometricSource::try_new(
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [-0.25, 0.0, 0.0]).unwrap(),
            100_000.0,
            1.0,
            AngularDistribution::uniform_cone(0.2).unwrap(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap();
        let mut authored = [fixture(
            41,
            RigidTransform::try_new([0.0, 0.0, FRAC_1_SQRT_2, FRAC_1_SQRT_2], [0.0, 0.0, 3.0])
                .unwrap(),
            FixtureOpticalState::try_new(0.5, [0.4; SPECTRAL_BINS]).unwrap(),
            vec![
                active,
                away,
                isotropic([0.0, 0.0, -3.0], 1.0, SourceSpectrum::equal_energy())
                    .dimmed(0.0)
                    .unwrap(),
            ],
        )];
        let rig = assemble(&authored, cx);
        close(
            rig.sources()[0].flux_lm(),
            12.0,
            2.0e-13,
            "retained calibration",
        );
        assert_eq!(rig.emitter_identities().len(), 3);
        let density = [0.2, 1.0, 0.6, 0.4, 1.4, 0.3, 0.7, 1.1];
        let medium = haze(cx, &density);
        let camera = PhysicalCamera::try_look_at(
            Point3::new(0.0, 0.0, 2.0),
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_vertical_fov(0.15).unwrap(),
            2.0,
            Aperture::try_circular(0.01).unwrap(),
        )
        .unwrap();
        let settings = Settings {
            sampler: Sampler::Iid,
            ..settings(2, 64)
        };
        let image = |policy| {
            PathLuminanceJob::new_with_rig_camera(
                &scene,
                &camera,
                &rig,
                &[],
                PathMedium::Scattering(&medium),
                TransportOptions::default(),
                settings,
                cx,
            )
            .unwrap()
            .with_point_sampling(policy)
            .unwrap()
        };
        let probes = probes();
        let probe = |policy| {
            probe_job(
                &scene,
                &rig,
                PathMedium::Scattering(&medium),
                &probes,
                settings,
                policy,
                cx,
            )
        };
        let mut reference_image = image(PointFixtureSampling::All);
        let mut reference_probe = probe(PointFixtureSampling::All);
        reference_image.advance(cx, usize::MAX).unwrap();
        reference_probe.advance(cx, usize::MAX).unwrap();
        for connections in [1, MAX_POINT_CONNECTIONS_PER_VERTEX] {
            let mut sampled_image = image(cdf(connections));
            let mut sampled_probe = probe(cdf(connections));
            sampled_image.advance(cx, usize::MAX).unwrap();
            sampled_probe.advance(cx, usize::MAX).unwrap();
            assert_eq!(image_bits(&sampled_image), image_bits(&reference_image));
            assert_eq!(probe_bits(&sampled_probe), probe_bits(&reference_probe));
            let value = sampled_image.pixel_estimate(0).unwrap();
            assert!(value.luminance_cd_m2 > 0.0 && value.total_error_bound_cd_m2.is_none());
            assert_eq!(
                sampled_image.emitter_identities(),
                Some(rig.emitter_identities())
            );
            assert!(sampled_probe.probe_estimate(0).unwrap().direct_lux > 0.0);
        }
        authored[0].optics = FixtureOpticalState::try_new(0.0, [1.0; SPECTRAL_BINS]).unwrap();
        let dark_rig = assemble(&authored, cx);
        let mut dark = image_job(&scene, &dark_rig, PathMedium::Vacuum, settings, cdf(32), cx);
        dark.advance(cx, usize::MAX).unwrap();
        assert_eq!(dark.pixel_estimate(0).unwrap().luminance_cd_m2, 0.0);
    });
}

#[test]
fn cdf_mixed_rig_replays_flat_parallel_and_cancelled_paths_with_policy_reports() {
    for sampler in [Sampler::Iid, Sampler::OwenSobol] {
        with_cx(|cx, gate| {
            let mut scene = receiver_scene();
            scene.camera.half_tan = 0.2;
            scene.primitives[0].material = Material::Lambertian {
                reflectance: LiftedSpectrum {
                    c: [0.0, 0.5, -0.2],
                },
            };
            scene
                .primitives
                .push(floor(5.0, 8.0, LiftedSpectrum { c: [0.0; 3] }));
            let inputs: Vec<_> = (0..2)
                .map(|index| {
                    let spectrum = SourceSpectrum::try_from_bins(std::array::from_fn(|bin| {
                        let x = bin as f64 / (SPECTRAL_BINS - 1) as f64;
                        0.2 + if index == 0 { 2.0 * x } else { 2.0 * (1.0 - x) }
                    }))
                    .unwrap();
                    fixture(
                        51 + index,
                        RigidTransform::try_new(
                            [0.0, 0.0, 0.0, 1.0],
                            [index as f64 * 2.0 - 1.0, 0.0, 3.0],
                        )
                        .unwrap(),
                        FixtureOpticalState::try_new(
                            0.4 + 0.2 * index as f64,
                            [0.8; SPECTRAL_BINS],
                        )
                        .unwrap(),
                        vec![
                            isotropic([-0.2, 0.0, 0.0], 80.0, spectrum)
                                .dimmed(0.5)
                                .unwrap(),
                            PhotometricSource::try_new(
                                RigidTransform::try_new([1.0, 0.0, 0.0, 0.0], [0.2, 0.1, 0.2])
                                    .unwrap(),
                                70.0,
                                0.25,
                                AngularDistribution::cosine_power(1.0, FRAC_PI_2).unwrap(),
                                SourceSpectrum::equal_energy(),
                            )
                            .unwrap(),
                        ],
                    )
                })
                .collect();
            let rig = assemble(&inputs, cx);
            let density = [0.2, 1.0, 0.6, 0.4, 1.4, 0.3, 0.7, 1.1];
            let medium = haze(cx, &density);
            let settings = Settings {
                width: 3,
                height: 2,
                sampler,
                ..settings(3, 16)
            };
            let probes = probes();
            let image = || {
                image_job(
                    &scene,
                    &rig,
                    PathMedium::Scattering(&medium),
                    settings,
                    cdf(8),
                    cx,
                )
            };
            let probe = |policy| {
                probe_job(
                    &scene,
                    &rig,
                    PathMedium::Scattering(&medium),
                    &probes,
                    settings,
                    policy,
                    cx,
                )
            };
            let mut reference_image = image();
            let mut reference_probe = probe(cdf(8));
            reference_image.advance(cx, usize::MAX).unwrap();
            reference_probe.advance(cx, usize::MAX).unwrap();
            let mut flat = PathLuminanceJob::new(
                &scene,
                rig.sources(),
                PathMedium::Scattering(&medium),
                TransportOptions::default(),
                settings,
                cx,
            )
            .unwrap()
            .with_point_sampling(cdf(8))
            .unwrap();
            flat.advance(cx, usize::MAX).unwrap();
            assert_eq!(image_bits(&flat), image_bits(&reference_image));
            let mut full_direct = probe(PointFixtureSampling::All);
            full_direct.advance(cx, probes.len()).unwrap();
            for index in 0..probes.len() {
                assert_eq!(
                    full_direct
                        .probe_estimate(index)
                        .unwrap()
                        .direct_lux
                        .to_bits(),
                    reference_probe
                        .probe_estimate(index)
                        .unwrap()
                        .direct_lux
                        .to_bits()
                );
            }
            let execution = RenderExecutionConfig::try_new(2, 1, 4, 64 << 20, RunId(2301)).unwrap();
            let mut image = image();
            let mut probe = probe(cdf(8));
            image.advance(cx, 3).unwrap();
            probe.advance(cx, 1).unwrap();
            image.advance_with_execution(cx, 7, &execution).unwrap();
            probe.advance_with_execution(cx, 3, &execution).unwrap();
            let image_prefix = image_bits(&image);
            let probe_prefix = probe_bits(&probe);
            assert!(image.set_point_sampling(PointFixtureSampling::All).is_err());
            assert!(probe.set_point_sampling(PointFixtureSampling::All).is_err());
            gate.request();
            assert!(matches!(
                image.advance_with_execution(cx, 20, &execution),
                Err(PhotometryError::Cancelled)
            ));
            assert!(matches!(
                probe.advance(cx, 20),
                Err(PhotometryError::Cancelled)
            ));
            assert_eq!(image_bits(&image), image_prefix);
            assert_eq!(probe_bits(&probe), probe_prefix);
            with_cx(|resumed, _| {
                while !image.is_complete() {
                    image
                        .advance_with_execution(resumed, 7, &execution)
                        .unwrap();
                    image.advance(resumed, 2).unwrap();
                }
                while !probe.is_complete() {
                    probe.advance(resumed, 3).unwrap();
                    probe
                        .advance_with_execution(resumed, 5, &execution)
                        .unwrap();
                }
            });
            assert_eq!(image_bits(&image), image_bits(&reference_image));
            assert_eq!(probe_bits(&probe), probe_bits(&reference_probe));
            assert_eq!(image.emitter_identities(), Some(rig.emitter_identities()));
            assert_eq!(probe.emitter_identities(), Some(rig.emitter_identities()));
            assert_eq!(image.error_budget_cd_m2(0).unwrap().point_sampling, cdf(8));
            assert_eq!(probe.error_budget_lux(0).unwrap().point_sampling, cdf(8));
            for value in image.estimates().flatten() {
                assert_eq!(
                    value.estimated_standard_error_cd_m2.is_some(),
                    sampler == Sampler::Iid
                );
                assert!(value.total_error_bound_cd_m2.is_none());
            }
        });
    }
}

#[test]
fn cdf_invalid_policy_and_unselected_sources_refuse_before_path_publication() {
    with_cx(|cx, _| {
        let scene = receiver_scene();
        let build = |sources| {
            assemble(
                &[fixture(
                    71,
                    RigidTransform::identity(),
                    FixtureOpticalState::default(),
                    sources,
                )],
                cx,
            )
        };
        let rig = build(vec![isotropic(
            [0.0, 0.0, 3.0],
            90.0,
            SourceSpectrum::equal_energy(),
        )]);
        let settings = settings(1, 8);
        let probes = probes();
        let mut image = image_job(&scene, &rig, PathMedium::Vacuum, settings, cdf(1), cx);
        let mut probe = probe_job(
            &scene,
            &rig,
            PathMedium::Vacuum,
            &probes,
            settings,
            cdf(1),
            cx,
        );
        for invalid in [0, MAX_POINT_CONNECTIONS_PER_VERTEX + 1] {
            assert!(matches!(
                image.set_point_sampling(cdf(invalid)),
                Err(PhotometryError::InvalidOptions)
            ));
            assert!(matches!(
                probe.set_point_sampling(cdf(invalid)),
                Err(PhotometryError::InvalidOptions)
            ));
            assert_eq!(image.point_sampling(), cdf(1));
            assert_eq!(probe.point_sampling(), cdf(1));
        }
        for refusal in 0..3 {
            let unselected = if refusal == 2 {
                isotropic([0.0, 0.0, 3.0], 1.0e-30, SourceSpectrum::equal_energy())
            } else {
                PhotometricSource::try_new(
                    RigidTransform::try_new(
                        [0.0, 0.0, 0.0, 1.0],
                        if refusal == 0 {
                            [20.0, 0.0, 1.0]
                        } else {
                            [0.0, 0.0, 1.0e-8]
                        },
                    )
                    .unwrap(),
                    1.0e-20,
                    1.0,
                    AngularDistribution::uniform_cone(0.2).unwrap(),
                    SourceSpectrum::equal_energy(),
                )
                .unwrap()
            };
            // A zero-intensity beam still owes full geometric admission;
            // positive importance below represented support also refuses.
            let rig = build(vec![
                isotropic([0.0, 0.0, 3.0], 90.0, SourceSpectrum::equal_energy()),
                unselected,
            ]);
            let mut job = PathLuminanceJob::new_with_rig(
                &scene,
                &rig,
                &[],
                PathMedium::Vacuum,
                TransportOptions {
                    max_distance_m: 10.0,
                    ..TransportOptions::default()
                },
                settings,
                cx,
            )
            .unwrap()
            .with_point_sampling(cdf(1))
            .unwrap();
            assert!(matches!(
                job.advance(cx, 1),
                Err(PhotometryError::Geometry(TracerError::FixtureConnection {
                    source: 1,
                    ..
                }))
            ));
            assert_eq!(job.completed_path_samples(), 0);
            assert!(job.pixel_estimate(0).is_none());
            assert_eq!(
                job.emitter_identities().unwrap()[1],
                FixtureEmitterIdentity {
                    fixture: FixtureId(71),
                    emitter: EmitterId(2)
                }
            );
        }
    });
}

#[test]
fn cdf_three_hundred_multi_emitter_fixtures_preserve_units_and_native_source_sum() {
    with_cx(|cx, _| {
        let mut expected_lux = 0.0;
        let fixtures: Vec<_> = (0..300)
            .map(|index| {
                let x = (index % 20) as f64 - 9.5;
                let y = (index / 20) as f64 - 7.0;
                let shared_dimmer = [0.0, 0.25, 0.5, 1.0][index % 4];
                let sources = (0..24)
                    .map(|emitter| {
                        let local_x = (emitter as f64 - 11.5) * 0.02;
                        let dimmer = [0.25, 0.5, 1.0][emitter % 3];
                        let r2 = (x + local_x) * (x + local_x) + y * y + 64.0;
                        expected_lux +=
                            10.0 * dimmer * shared_dimmer * 0.4 * 8.0 / (r2 * r2.sqrt());
                        isotropic([local_x, 0.0, 0.0], 10.0, SourceSpectrum::equal_energy())
                            .dimmed(dimmer)
                            .unwrap()
                    })
                    .collect();
                fixture(
                    index as u64 + 1,
                    RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [x, y, 8.0]).unwrap(),
                    FixtureOpticalState::try_new(shared_dimmer, [0.4; SPECTRAL_BINS]).unwrap(),
                    sources,
                )
            })
            .collect();
        let rig = assemble(&fixtures, cx);
        assert_eq!(rig.sources().len(), 7200);
        assert!(
            rig.fixtures()
                .iter()
                .all(|fixture| fixture.sources.len() == 24)
        );
        assert_eq!(
            rig.emitter_identities()[7199],
            FixtureEmitterIdentity {
                fixture: FixtureId(300),
                emitter: EmitterId(24)
            }
        );
        let scene = receiver_scene();
        let settings = settings(1, 8);
        let mut all = image_job(
            &scene,
            &rig,
            PathMedium::Vacuum,
            settings,
            PointFixtureSampling::All,
            cx,
        );
        all.advance(cx, usize::MAX).unwrap();
        let reference = all.pixel_estimate(0).unwrap();
        for connections in [1, MAX_POINT_CONNECTIONS_PER_VERTEX] {
            let mut sampled = image_job(
                &scene,
                &rig,
                PathMedium::Vacuum,
                settings,
                cdf(connections),
                cx,
            );
            sampled.advance(cx, usize::MAX).unwrap();
            let value = sampled.pixel_estimate(0).unwrap();
            // In this clear gray scene C_j is proportional to its proposal
            // weight. The same wavelengths isolate selection arithmetic;
            // this tolerance is not a physical-error or convergence bound.
            for (actual, expected) in value.xyz.into_iter().zip(reference.xyz) {
                close(actual, expected, 2.0e-10, "large-rig source sum");
            }
            assert!(value.total_error_bound_cd_m2.is_none());
        }
        let probe_points = [probes()[0]];
        let mut probe = probe_job(
            &scene,
            &rig,
            PathMedium::Vacuum,
            &probe_points,
            settings,
            cdf(32),
            cx,
        );
        probe.advance(cx, 1).unwrap();
        close(
            probe.probe_estimate(0).unwrap().direct_lux,
            expected_lux,
            2.0e-12,
            "independent 7200-emitter lux sum",
        );
    });
}
