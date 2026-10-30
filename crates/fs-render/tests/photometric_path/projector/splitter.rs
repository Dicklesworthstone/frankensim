//! Independent references for passive paraxial beam copies and spectral filters.
//! These tests describe the declared exit light field, not Snell refraction
//! through a physical prism or wavelength-dependent diffraction.

use super::*;
use crate::spectral_reference::cie_reference_xyz;
use fs_render::photometry::area::{ParaxialBeamSplitter, ProjectorBeam};

type Beam = ([f64; 2], f64);

fn split_emitter(case: &ProjectorCase, beams: &[Beam]) -> PhotometricAreaEmitter {
    let source = case.emitter(1.0, RigidTransform::identity());
    let splitter = ParaxialBeamSplitter::try_new(
        beams
            .iter()
            .map(|&(offset, throughput)| ProjectorBeam::try_new(offset, throughput).unwrap())
            .collect(),
    )
    .unwrap();
    let projector = source
        .projector()
        .unwrap()
        .clone()
        .with_beam_splitter(splitter)
        .unwrap();
    source.with_projector(projector).unwrap()
}

fn split_lux(case: &ProjectorCase, beams: &[Beam], point: [f64; 2], plane: f64) -> f64 {
    // In physical screen coordinates, each paraxial slope offset translates
    // the original finite-pupil image by D*d, at focus and out of focus.
    beams
        .iter()
        .map(|&(offset, throughput)| {
            throughput
                * case.lux(
                    [point[0] - plane * offset[0], point[1] - plane * offset[1]],
                    plane,
                )
        })
        .sum()
}

#[test]
fn g1_splitter_focused_images_retain_shift_sign_weights_and_total_flux() {
    with_cx(|cx, _| {
        let focus = 3.0;
        let case = ProjectorCase::new(focus, [2, 2], vec![1.0, 0.25, 0.5, 0.0]);
        let beams = [([-0.5, 0.0], 0.2), ([0.0, 0.0], 0.3), ([0.5, 0.0], 0.25)];
        let scene = case.scene();
        let sources = [split_emitter(&case, &beams)];
        let cell_area =
            case.gobo_size[0] * case.gobo_size[1] * (focus / case.distance).powi(2) / 4.0;
        let expected_flux = case.clear_flux * (1.0 + 0.25 + 0.5) / 4.0 * 0.75;
        close(
            sources[0].flux_lm(),
            expected_flux,
            1.0e-15,
            "passive split flux",
        );
        close(
            sources[0].clear_flux_lm(),
            case.clear_flux,
            1.0e-15,
            "splitter preserves the clear calibration",
        );
        let mut integrated_flux = 0.0;
        // The copies do not overlap. Their twelve constant screen cells
        // partition every transmitted ray, including the three black cells.
        for (offset, _) in beams {
            for gobo in [
                [-0.01, -0.005],
                [0.01, -0.005],
                [-0.01, 0.005],
                [0.01, 0.005],
            ] {
                let point = [
                    -focus / case.distance * gobo[0] + focus * offset[0],
                    -focus / case.distance * gobo[1] + focus * offset[1],
                ];
                let result = probe(
                    &scene,
                    &sources,
                    Point3::new(point[0], point[1], focus),
                    Vec3::new(0.0, 0.0, -1.0),
                    8192,
                    cx,
                );
                close(
                    result.illuminance_lux,
                    split_lux(&case, &beams, point, focus),
                    0.003,
                    "weighted translated focused cell",
                );
                integrated_flux += result.illuminance_lux * cell_area;
            }
        }
        close(
            integrated_flux,
            expected_flux,
            0.003,
            "integrated split screen flux",
        );
    });
}

#[test]
fn g1_splitter_defocus_overlap_and_absorbed_power_match_pupil_integrals() {
    with_cx(|cx, _| {
        let case = ProjectorCase::new(3.0, [2, 1], vec![1.0, 0.0]);
        let beams = [([-0.06, 0.0], 0.2), ([0.0, 0.0], 0.35), ([0.06, 0.0], 0.25)];
        let scene = case.scene();
        let sources = [split_emitter(&case, &beams)];
        for plane in [1.5, 6.0] {
            let c = 1.0 - case.distance / case.focal_length + case.distance / plane;
            let blur = plane / case.distance * c.abs() * case.half_pupil;
            for fraction in [-0.5, 0.0, 0.5] {
                let point = [fraction * blur, 0.0];
                let expected = split_lux(&case, &beams, point, plane);
                assert!(expected > 0.0);
                let result = probe(
                    &scene,
                    &sources,
                    Point3::new(point[0], point[1], plane),
                    Vec3::new(0.0, 0.0, -1.0),
                    32_768,
                    cx,
                );
                close(
                    result.illuminance_lux,
                    expected,
                    0.035,
                    "defocused beam-copy overlap",
                );
            }
        }
        for point in [[0.2, 0.0], [-0.4, 0.0]] {
            let result = probe(
                &scene,
                &sources,
                Point3::new(point[0], point[1], 3.0),
                Vec3::new(0.0, 0.0, -1.0),
                8192,
                cx,
            );
            close(
                result.illuminance_lux,
                split_lux(&case, &beams, point, 3.0),
                0.003,
                "overlapping clear copies and shared black support",
            );
        }
        let opaque = [split_emitter(
            &case,
            &[([-0.06, 0.0], 0.0), ([0.06, 0.0], 0.0)],
        )];
        assert_eq!(opaque[0].flux_lm(), 0.0);
        assert_eq!(opaque[0].clear_flux_lm(), case.clear_flux);
        assert_eq!(
            probe(
                &scene,
                &opaque,
                Point3::new(0.2, 0.0, 3.0),
                Vec3::new(0.0, 0.0, -1.0),
                32,
                cx,
            )
            .illuminance_lux,
            0.0
        );
        assert_eq!(
            pixel(&scene, &opaque, PathMedium::Vacuum, settings(1, 32), cx).luminance_cd_m2,
            0.0
        );
    });
}

#[test]
fn g1_splitter_nee_bsdf_and_mis_share_the_absolute_surface_reference() {
    with_cx(|cx, _| {
        let mut case = ProjectorCase::new(1.0, [1, 1], vec![1.0]);
        case.half_pupil = 0.5;
        case.gobo_size = [0.2, 0.2];
        case.distance = 0.4;
        case.focal_length = 0.4 / 1.4;
        let beams = [
            ([-0.12, 0.0], 0.2),
            ([0.04, 0.05], 0.35),
            ([0.14, -0.07], 0.25),
        ];
        let mut scene = case.scene();
        scene.camera.eye = Point3::new(0.0, 0.0, 0.75);
        scene.camera.forward = Vec3::new(0.0, 0.0, 1.0);
        scene
            .primitives
            .push(floor(1.0, 2.0, LiftedSpectrum { c: [0.0; 3] }));
        let sources = [split_emitter(&case, &beams)];
        let expected = 0.5 * split_lux(&case, &beams, [0.0, 0.0], 1.0) / PI;
        for strategy in [
            DirectStrategy::NeeOnly,
            DirectStrategy::BsdfOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
        ] {
            let value = pixel(
                &scene,
                &sources,
                PathMedium::Vacuum,
                Settings {
                    strategy,
                    sampler: Sampler::Iid,
                    ..settings(2, 65_536)
                },
                cx,
            );
            let error = value.estimated_standard_error_cd_m2.unwrap();
            assert!(
                (value.luminance_cd_m2 - expected).abs() <= 6.0 * error + 0.003 * expected,
                "{strategy:?}: observed={} reference={expected} estimated_se={error}",
                value.luminance_cd_m2
            );
            close(
                value.luminance_cd_m2,
                expected,
                0.04,
                "splitter surface estimators",
            );
        }
    });
}

#[test]
fn g1_finite_spectral_filters_preserve_passband_power_optics_and_absorption() {
    with_cx(|cx, _| {
        let case = ProjectorCase::new(3.0, [1, 1], vec![1.0]);
        let beams = [([0.04, 0.0], 0.4), ([-0.03, 0.02], 0.35)];
        let projector = split_emitter(&case, &beams).projector().unwrap().clone();
        let mut scene = case.scene();
        let outgoing = Vec3::new(0.04, 0.0, 1.0);
        let outgoing = outgoing.scale(1.0 / outgoing.norm());
        scene.camera.eye = Point3::new(0.0, 0.0, 0.0).offset(outgoing.scale(2.0));
        scene.camera.forward = outgoing.scale(-1.0);
        let density = [1.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-2.0, -2.0, 0.0], [4.0, 4.0, 3.0]),
            MassAbsorptionSpectrum::table(&[(380.0, 0.1), (780.0, 1.1)]).unwrap(),
            None,
            medium_budget(),
        )
        .unwrap();
        let mut bins = [0.0; SPECTRAL_BINS];
        bins[12] = 1.0;
        bins[60] = 1.0;
        let spectrum = SourceSpectrum::try_from_bins(bins).unwrap();
        let photopic = |bin| cie_reference_xyz(wavelength_nm(bin))[1];
        let observer_sum = photopic(12) + photopic(60);
        for projected in [false, true] {
            let original = PhotometricAreaEmitter::try_new(
                0,
                case.clear_flux,
                1.0,
                spectrum.clone(),
                EmissionSidedness::Front,
            )
            .unwrap();
            let original = if projected {
                original.with_projector(projector.clone()).unwrap()
            } else {
                original
            };
            let clear_luminance = case.clear_flux / (4.0 * case.half_pupil.powi(2))
                * if projected {
                    case.distance.powi(2)
                        / (case.gobo_size[0] * case.gobo_size[1] * outgoing.z.powi(4))
                        * 0.75
                } else {
                    1.0 / PI
                };
            for (bin, transmission) in [(12, 0.6), (60, 0.8)] {
                let mut filter = [0.0; SPECTRAL_BINS];
                filter[bin] = transmission;
                let filtered = original.filtered(filter).unwrap();
                let fraction = transmission * photopic(bin) / observer_sum;
                close(
                    filtered.flux_lm(),
                    original.flux_lm() * fraction,
                    2.0e-14,
                    "filtered output flux",
                );
                close(
                    filtered.clear_flux_lm(),
                    original.clear_flux_lm() * fraction,
                    2.0e-14,
                    "filtered clear calibration",
                );
                let point_source = PhotometricSource::try_new(
                    RigidTransform::identity(),
                    original.flux_lm(),
                    1.0,
                    AngularDistribution::isotropic(),
                    spectrum.clone(),
                )
                .unwrap()
                .filtered(filter)
                .unwrap();
                close(
                    point_source.flux_lm(),
                    filtered.flux_lm(),
                    2.0e-14,
                    "point and finite sources retain the same filtered lumens",
                );
                let before_optics = PhotometricAreaEmitter::try_new(
                    0,
                    case.clear_flux,
                    1.0,
                    spectrum.clone(),
                    EmissionSidedness::Front,
                )
                .unwrap()
                .filtered(filter)
                .unwrap();
                let before_optics = if projected {
                    before_optics.with_projector(projector.clone()).unwrap()
                } else {
                    before_optics
                };
                close(
                    before_optics.flux_lm(),
                    filtered.flux_lm(),
                    2.0e-14,
                    "filtering commutes with the passive optical losses",
                );
                for (sample, transmission) in filter.iter().enumerate() {
                    let expected_power = original.flux_lm()
                        * original.spectrum().bins_w_per_lumen_nm()[sample]
                        * transmission;
                    let actual_power =
                        filtered.flux_lm() * filtered.spectrum().bins_w_per_lumen_nm()[sample];
                    close(
                        actual_power,
                        expected_power,
                        2.0e-14,
                        "per-bin passive filter product",
                    );
                    close(
                        point_source.flux_lm()
                            * point_source.spectrum().bins_w_per_lumen_nm()[sample],
                        actual_power,
                        2.0e-14,
                        "point and finite source spectral power",
                    );
                    close(
                        before_optics.flux_lm()
                            * before_optics.spectrum().bins_w_per_lumen_nm()[sample],
                        actual_power,
                        2.0e-14,
                        "filter and optics order preserves each spectral bin",
                    );
                }
                let replacement = case
                    .emitter(1.0, RigidTransform::identity())
                    .projector()
                    .unwrap()
                    .clone();
                let replaced = filtered.clone().with_projector(replacement).unwrap();
                assert_eq!(replaced.clear_flux_lm(), filtered.clear_flux_lm());
                close(
                    replaced.flux_lm(),
                    filtered.clear_flux_lm(),
                    2.0e-14,
                    "replacement clear optics retain the filtered clear calibration",
                );
                assert_eq!(filtered.projector().is_some(), projected);
                if let Some(retained) = filtered.projector() {
                    assert_eq!(retained.focus_distance_m(), projector.focus_distance_m());
                    let retained = retained.beam_splitter().unwrap().beams();
                    assert_eq!(retained.len(), beams.len());
                    for (actual, &(offset, throughput)) in retained.iter().zip(&beams) {
                        assert_eq!(actual.slope_offset(), offset);
                        assert_eq!(actual.throughput(), throughput);
                    }
                }
                let value = pixel(
                    &scene,
                    &[filtered],
                    PathMedium::Absorbing(&medium),
                    settings(1, 65_536),
                    cx,
                );
                let kappa = 0.1 + (wavelength_nm(bin) - 380.0) / 400.0;
                close(
                    value.luminance_cd_m2,
                    clear_luminance * fraction * (-2.0 * kappa).exp(),
                    0.002,
                    "filtered source-hit spectrum through physical absorption",
                );
            }
            let black = original.filtered([0.0; SPECTRAL_BINS]).unwrap();
            assert_eq!(black.flux_lm(), 0.0);
            assert_eq!(black.clear_flux_lm(), 0.0);
            assert_eq!(black.projector().is_some(), projected);
            assert_eq!(black.spectrum(), original.spectrum());
        }
    });
}

#[test]
fn g1_filtered_splitter_haze_matches_independent_single_scattering() {
    with_cx(|cx, _| {
        let mut case = ProjectorCase::new(3.0, [1, 1], vec![1.0]);
        case.gobo_size = [0.08, 0.08];
        let beams = [
            ([-0.02, 0.0], 0.2),
            ([0.0, 0.01], 0.3),
            ([0.02, -0.01], 0.25),
        ];
        let mut scene = case.scene();
        scene.camera.eye = Point3::new(0.1, 0.0, 2.0);
        let filter = std::array::from_fn(|bin| if wavelength_nm(bin) < 550.0 { 0.2 } else { 0.7 });
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for (bin, transmission) in filter.iter().enumerate() {
            let weight = cie_reference_xyz(wavelength_nm(bin))[1];
            numerator += transmission * weight;
            denominator += weight;
        }
        let fraction = numerator / denominator;
        let sources = [split_emitter(&case, &beams).filtered(filter).unwrap()];
        let density = [1.0];
        let medium = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-2.0, -2.0, 0.5], [4.0, 4.0, 1.0]),
                MassAbsorptionSpectrum::constant(0.2).unwrap(),
                None,
                medium_budget(),
            )
            .unwrap(),
            MassScatteringSpectrum::constant(0.4).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        // Every pupil/event pair stays inside all three clear gobo copies:
        // |g_x| < 0.034 m and |g_y| < 0.013 m, both below the 0.04 m edge.
        // Therefore the independent unsplit volume integral scales by the
        // total transmitted ray power, including the non-gray filter.
        let reference = single_scatter_reference(&case, 64);
        close(
            single_scatter_reference(&case, 32),
            reference,
            0.0003,
            "haze quadrature refinement",
        );
        let value = pixel(
            &scene,
            &sources,
            PathMedium::Scattering(&medium),
            settings(1, 65_536),
            cx,
        );
        close(
            value.luminance_cd_m2,
            reference * 0.75 * fraction,
            0.035,
            "spectrally filtered multi-beam haze with both Beer legs",
        );
    });
}

#[test]
fn g4_g5_splitter_filter_cancellation_and_mixed_execution_preserve_bits() {
    with_cx(|cx, gate| {
        let case = ProjectorCase::new(3.0, [2, 2], vec![1.0, 0.25, 0.5, 0.0]);
        let scene = screen_scene(&case);
        let sources = [split_emitter(
            &case,
            &[
                ([-0.05, 0.0], 0.2),
                ([0.0, 0.03], 0.3),
                ([0.06, -0.02], 0.25),
            ],
        )
        .filtered(std::array::from_fn(|bin| if bin < 40 { 0.3 } else { 0.8 }))
        .unwrap()];
        let new_job = || {
            PathLuminanceJob::new_with_emitters(
                &scene,
                PhotometricEmitters {
                    points: &[],
                    surfaces: &sources,
                },
                PathMedium::Vacuum,
                TransportOptions::default(),
                Settings {
                    width: 4,
                    height: 3,
                    sampler: Sampler::Iid,
                    ..settings(2, 64)
                },
                cx,
            )
            .unwrap()
        };
        let mut reference = new_job();
        reference.advance(cx, usize::MAX).unwrap();
        assert!(
            reference
                .estimates()
                .flatten()
                .any(|value| value.luminance_cd_m2 > 0.0)
        );
        let mut resumed = new_job();
        assert_eq!(resumed.advance(cx, 7).unwrap(), 7);
        let before = image_bits(&resumed);
        gate.request();
        let execution =
            RenderExecutionConfig::try_new(2, 2, 2, 64 << 20, RunId(0x73706c31)).unwrap();
        assert!(matches!(
            resumed.advance_with_execution(cx, 19, &execution),
            Err(PhotometryError::Cancelled)
        ));
        assert_eq!(image_bits(&resumed), before);
        with_cx(|next, _| {
            let execution =
                RenderExecutionConfig::try_new(3, 1, 4, 64 << 20, RunId(0x73706c32)).unwrap();
            assert_eq!(
                resumed
                    .advance_with_execution(next, 23, &execution)
                    .unwrap(),
                23
            );
            resumed.advance(next, usize::MAX).unwrap();
        });
        assert!(resumed.is_complete());
        assert_eq!(image_bits(&resumed), image_bits(&reference));
    });
}

#[test]
#[ignore = "run explicitly to compare unsplit, three-beam and eight-beam projector throughput"]
fn photometric_splitter_throughput() {
    with_cx(|cx, _| {
        let case = ProjectorCase::new(3.0, [2, 2], vec![1.0, 0.25, 0.5, 0.0]);
        let scene = screen_scene(&case);
        let execution =
            RenderExecutionConfig::try_new(4, 4, 4, 64 << 20, RunId(0x73706c33)).unwrap();
        for count in [0_usize, 3, 8] {
            let sources = if count == 0 {
                [case.emitter(1.0, RigidTransform::identity())]
            } else {
                let beams: Vec<_> = (0..count)
                    .map(|index| {
                        (
                            [0.16 * index as f64 / (count - 1) as f64 - 0.08, 0.0],
                            1.0 / count as f64,
                        )
                    })
                    .collect();
                [split_emitter(&case, &beams)]
            };
            let mut replay = None;
            let mut times = [0.0; 3];
            for (repeat, elapsed) in times.iter_mut().enumerate() {
                let mut job = PathLuminanceJob::new_with_emitters(
                    &scene,
                    PhotometricEmitters {
                        points: &[],
                        surfaces: &sources,
                    },
                    PathMedium::Vacuum,
                    TransportOptions::default(),
                    Settings {
                        width: 16,
                        height: 16,
                        sampler: Sampler::Iid,
                        ..settings(2, 64)
                    },
                    cx,
                )
                .unwrap();
                let start = std::time::Instant::now();
                job.advance_with_execution(cx, usize::MAX, &execution)
                    .unwrap();
                *elapsed = start.elapsed().as_secs_f64();
                assert!(job.is_complete());
                assert_eq!(job.completed_path_samples(), 16_384);
                let bits = image_bits(&job);
                if let Some(expected) = &replay {
                    assert_eq!(&bits, expected, "beam_count={count} repeat={repeat}");
                } else {
                    assert!(
                        job.estimates()
                            .flatten()
                            .any(|value| value.luminance_cd_m2 > 0.0)
                    );
                    replay = Some(bits);
                }
            }
            times.sort_by(f64::total_cmp);
            eprintln!(
                "photometric_splitter_throughput beam_count={count} workers=4 paths=16384 repeats=3 seconds_median={:.6} paths_per_second={:.1} seconds_sorted={times:?}",
                times[1],
                16_384.0 / times[1]
            );
        }
    });
}
