//! Independent power-law references for a passive paraxial box diffuser.
//! Pupil and slope deflections are independent uniform variables in this
//! declared model. Their exact CDF supplies profiles without the production
//! gobo rectangle integrator, projector evaluator, source sampler or BSDF.

use super::*;
use crate::spectral_reference::cie_reference_xyz;
use fs_render::photometry::area::{ParaxialBeamSplitter, ParaxialBoxDiffuser, ProjectorBeam};

type Beam = ([f64; 2], f64);

#[derive(Clone, Copy)]
struct Frost {
    slope: [f64; 2],
    scattered: f64,
    transmission: f64,
}

fn source(
    case: &ProjectorCase,
    frost: Option<Frost>,
    beams: &[Beam],
    frame: RigidTransform,
) -> PhotometricAreaEmitter {
    let source = case.emitter(1.0, frame);
    let mut optics = source.projector().unwrap().clone();
    if !beams.is_empty() {
        optics = optics
            .with_beam_splitter(
                ParaxialBeamSplitter::try_new(
                    beams
                        .iter()
                        .map(|&(offset, power)| ProjectorBeam::try_new(offset, power).unwrap())
                        .collect(),
                )
                .unwrap(),
            )
            .unwrap();
    }
    if let Some(frost) = frost {
        optics = optics
            .with_diffuser(
                ParaxialBoxDiffuser::try_new(
                    frost.slope.map(f64::atan),
                    frost.scattered,
                    frost.transmission,
                )
                .unwrap(),
            )
            .unwrap();
    }
    source.with_projector(optics).unwrap()
}

// CDF of the sum of two independent centred uniforms. The triangular tails
// and optional flat centre follow by integrating interval overlap. The
// degenerate limit is the original uniform (or a point mass), not a fit.
fn uniform_sum_cdf(x: f64, a: f64, b: f64) -> f64 {
    let wide = a.max(b);
    let narrow = a.min(b);
    if wide == 0.0 {
        return if x <= 0.0 { 0.0 } else { 1.0 };
    }
    if narrow == 0.0 {
        return (0.5 + x / (2.0 * wide)).clamp(0.0, 1.0);
    }
    if x <= -wide - narrow {
        0.0
    } else if x < -wide + narrow {
        (x + wide + narrow).powi(2) / (8.0 * wide * narrow)
    } else if x <= wide - narrow {
        0.5 + x / (2.0 * wide)
    } else if x < wide + narrow {
        1.0 - (wide + narrow - x).powi(2) / (8.0 * wide * narrow)
    } else {
        1.0
    }
}

fn mask_probability(
    case: &ProjectorCase,
    centre: [f64; 2],
    pupil_half: f64,
    diffuser_half: [f64; 2],
) -> f64 {
    let mut result = 0.0;
    for row in 0..case.gobo_shape[1] {
        for column in 0..case.gobo_shape[0] {
            let mut probability = 1.0;
            for (axis, cell) in [column, row].into_iter().enumerate() {
                let step = case.gobo_size[axis] / case.gobo_shape[axis] as f64;
                let low = -0.5 * case.gobo_size[axis] + cell as f64 * step;
                let high = low + step;
                probability *=
                    uniform_sum_cdf(high - centre[axis], pupil_half, diffuser_half[axis])
                        - uniform_sum_cdf(low - centre[axis], pupil_half, diffuser_half[axis]);
            }
            result += probability * case.transmission[row * case.gobo_shape[0] + column];
        }
    }
    result
}

fn reference_lux(
    case: &ProjectorCase,
    frost: Frost,
    beams: &[Beam],
    point: [f64; 2],
    plane: f64,
) -> f64 {
    // Screen position p = D*(a*(1/s-1/f) - g/s + delta) + a.
    // Therefore g = c*a - s*p/D + s*delta. The two X/Y probabilities
    // factor, and dA_screen/dA_gobo=(D/s)^2 supplies the SI lux scale.
    let c = 1.0 - case.distance / case.focal_length + case.distance / plane;
    let pupil_half = if c.abs() < 1.0e-14 {
        0.0 // This fixture's exact focus, with only rounded lens distances.
    } else {
        c.abs() * case.half_pupil
    };
    let default = [([0.0; 2], 1.0)];
    let beams = if beams.is_empty() { &default } else { beams };
    let mut coverage = 0.0;
    for &(offset, power) in beams {
        let centre =
            core::array::from_fn(|axis| -case.distance * (point[axis] / plane - offset[axis]));
        let clear = mask_probability(case, centre, pupil_half, [0.0; 2]);
        let diffused = mask_probability(
            case,
            centre,
            pupil_half,
            frost.slope.map(|spread| case.distance * spread),
        );
        coverage += power * ((1.0 - frost.scattered) * clear + frost.scattered * diffused);
    }
    case.clear_flux * case.distance.powi(2)
        / (case.gobo_size[0] * case.gobo_size[1] * plane.powi(2))
        * frost.transmission
        * coverage
}

#[test]
fn g1_frost_focus_softens_declared_edges_and_preserves_passive_screen_flux() {
    with_cx(|cx, _| {
        let focus = 3.0;
        let case = ProjectorCase::new(focus, [1, 1], vec![1.0]);
        let scene = case.scene();
        let frost = Frost {
            slope: [0.07, 0.04],
            scattered: 1.0,
            transmission: 0.73,
        };
        let half_image = case
            .gobo_size
            .map(|size| 0.5 * size * focus / case.distance);
        let blur = frost.slope.map(|spread| focus * spread);
        let full = case.clear_flux * case.distance.powi(2)
            / (case.gobo_size[0] * case.gobo_size[1] * focus.powi(2));
        let sources = [source(&case, Some(frost), &[], RigidTransform::identity())];
        // The full-frost outer edge is a linear ramp of half-width q*tan(h).
        // Choose both axes and non-boundary outside points to detect swapped
        // anisotropic spreads, shifts, lost power, or a hard-edge substitute.
        for (point, fraction) in [
            ([half_image[0] - 0.5 * blur[0], 0.0], 0.75),
            ([half_image[0], 0.0], 0.5),
            ([half_image[0] + 0.5 * blur[0], 0.0], 0.25),
            ([0.0, -half_image[1] - 0.5 * blur[1]], 0.25),
            ([half_image[0] + 1.1 * blur[0], 0.0], 0.0),
        ] {
            let expected = full * frost.transmission * fraction;
            close(
                reference_lux(&case, frost, &[], point, focus),
                expected,
                2.0e-13,
                "independent focused box profile",
            );
            close(
                probe(
                    &scene,
                    &sources,
                    Point3::new(point[0], point[1], focus),
                    Vec3::new(0.0, 0.0, -1.0),
                    8192,
                    cx,
                )
                .illuminance_lux,
                expected,
                0.006,
                "native softened focus edge",
            );
        }
        // At focus the two-dimensional profile is bilinear between these
        // knots, including the retained sharp core. One midpoint per cell
        // integrates each polynomial exactly in the stated optical model.
        for scattered in [0.35, 1.0] {
            let frost = Frost { scattered, ..frost };
            let sources = [source(&case, Some(frost), &[], RigidTransform::identity())];
            let knots: [[f64; 6]; 2] = core::array::from_fn(|axis| {
                let a = half_image[axis];
                let b = blur[axis];
                [-a - b, -a, -a + b, a - b, a, a + b]
            });
            let mut reference_flux = 0.0;
            let mut observed_flux = 0.0;
            for x in knots[0].windows(2) {
                for y in knots[1].windows(2) {
                    let point = [0.5 * (x[0] + x[1]), 0.5 * (y[0] + y[1])];
                    let area = (x[1] - x[0]) * (y[1] - y[0]);
                    reference_flux += reference_lux(&case, frost, &[], point, focus) * area;
                    observed_flux += probe(
                        &scene,
                        &sources,
                        Point3::new(point[0], point[1], focus),
                        Vec3::new(0.0, 0.0, -1.0),
                        2048,
                        cx,
                    )
                    .illuminance_lux
                        * area;
                }
            }
            let expected = case.clear_flux * frost.transmission;
            close(
                reference_flux,
                expected,
                2.0e-13,
                "exact integrated frost power",
            );
            close(
                sources[0].flux_lm(),
                expected,
                2.0e-15,
                "reported passive frost flux",
            );
            close(
                observed_flux,
                expected,
                0.008,
                "native integrated frost power",
            );
        }
    });
}

#[test]
fn g1_frost_defocus_matches_independent_sum_of_pupil_and_slope_uniforms() {
    with_cx(|cx, _| {
        let case = ProjectorCase::new(3.0, [2, 1], vec![1.0, 0.0]);
        let scene = case.scene();
        let frost = Frost {
            slope: [0.035, 0.02],
            scattered: 0.65,
            transmission: 0.8,
        };
        let sources = [source(&case, Some(frost), &[], RigidTransform::identity())];
        for plane in [1.5, 6.0] {
            let c = 1.0 - case.distance / case.focal_length + case.distance / plane;
            let half_blur = plane / case.distance
                * (c.abs() * case.half_pupil + case.distance * frost.slope[0]);
            for fraction in [-0.8, -0.2, 0.2, 0.8] {
                let point = [fraction * half_blur, 0.0];
                let expected = reference_lux(&case, frost, &[], point, plane);
                assert!(expected > 0.0);
                close(
                    probe(
                        &scene,
                        &sources,
                        Point3::new(point[0], point[1], plane),
                        Vec3::new(0.0, 0.0, -1.0),
                        32_768,
                        cx,
                    )
                    .illuminance_lux,
                    expected,
                    0.03,
                    "independent defocused pupil-plus-diffuser distribution",
                );
            }
        }
    });
}

#[test]
fn g1_frost_splitter_filter_and_rigid_frame_preserve_absolute_power() {
    with_cx(|cx, _| {
        let case = ProjectorCase::new(3.0, [2, 2], vec![1.0, 0.25, 0.5, 0.0]);
        let frost = Frost {
            slope: [0.08, 0.05],
            scattered: 0.7,
            transmission: 0.6,
        };
        let beams = [
            ([-0.06, 0.01], 0.2),
            ([0.02, -0.03], 0.35),
            ([0.08, 0.04], 0.25),
        ];
        let filter = core::array::from_fn(|bin| if wavelength_nm(bin) < 550.0 { 0.2 } else { 0.7 });
        let y_sum = (0..SPECTRAL_BINS)
            .map(|bin| cie_reference_xyz(wavelength_nm(bin))[1])
            .sum::<f64>();
        let filter_fraction = filter
            .iter()
            .enumerate()
            .map(|(bin, transmission)| transmission * cie_reference_xyz(wavelength_nm(bin))[1])
            .sum::<f64>()
            / y_sum;
        let frame = RigidTransform::try_new(
            [
                0.0,
                0.0,
                core::f64::consts::FRAC_1_SQRT_2,
                core::f64::consts::FRAC_1_SQRT_2,
            ],
            [0.3, -0.2, 0.4],
        )
        .unwrap();
        let placement = RigidTransform::try_new(
            [
                0.0,
                core::f64::consts::FRAC_1_SQRT_2,
                0.0,
                core::f64::consts::FRAC_1_SQRT_2,
            ],
            [7.0, -3.0, 4.0],
        )
        .unwrap();
        let h = case.half_pupil;
        let vertices = [[-h, -h], [h, -h], [h, h], [-h, h]]
            .map(|p| {
                let p = frame.transform_point(Point3::new(p[0], p[1], 0.0));
                [p.x, p.y, p.z]
            })
            .to_vec();
        let instance = GeometryInstance::try_new(
            93,
            fs_blake3::hash_domain("org.frankensim.test.projector", b"frost-rigid-aperture"),
            SharedGeometry::mesh(TriMesh::new(vertices, vec![[0, 1, 2], [0, 2, 3]])),
            placement,
        )
        .unwrap();
        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        scene.primitives.push(Primitive {
            shape: Shape::Instance(instance),
            material: Material::Lambertian {
                reflectance: LiftedSpectrum { c: [0.0; 3] },
            },
            emission: None,
        });
        let unfiltered = source(&case, Some(frost), &beams, frame);
        let sources = [unfiltered.filtered(filter).unwrap()];
        let expected_flux = case.clear_flux * 0.4375 * 0.8 * frost.transmission * filter_fraction;
        close(
            sources[0].flux_lm(),
            expected_flux,
            3.0e-14,
            "combined absolute optical losses",
        );
        close(
            sources[0].clear_flux_lm(),
            case.clear_flux * filter_fraction,
            3.0e-14,
            "diffuser preserves filtered clear calibration",
        );
        let filtered_first = case
            .emitter(1.0, frame)
            .filtered(filter)
            .unwrap()
            .with_projector(unfiltered.projector().unwrap().clone())
            .unwrap();
        close(
            filtered_first.flux_lm(),
            expected_flux,
            3.0e-14,
            "filter/optics construction order",
        );
        // Keep both probes away from every translated sharp gobo edge.
        // At p_y=q*d_y the retained sharp core is discontinuous, and rounded
        // focus/rigid coordinates can select either adjacent transmission.
        for point in [[0.35, 0.11], [-0.45, -0.17]] {
            let world = placement
                .transform_point(frame.transform_point(Point3::new(point[0], point[1], 3.0)));
            let normal =
                placement.transform_vector(frame.transform_vector(Vec3::new(0.0, 0.0, -1.0)));
            close(
                probe(&scene, &sources, world, normal, 32_768, cx).illuminance_lux,
                reference_lux(&case, frost, &beams, point, 3.0) * filter_fraction,
                0.006,
                "native filtered frost and splitter in original rigid optical frame",
            );
        }
    });
}

fn haze_reference(case: &ProjectorCase, frost: Frost, cells: usize) -> f64 {
    let da = 2.0 * case.half_pupil / cells as f64;
    let dz = 1.0 / (2 * cells) as f64;
    let mut sum = 0.0;
    for iz in 0..2 * cells {
        let z = 0.5 + (iz as f64 + 0.5) * dz;
        for ix in 0..cells {
            let a = -case.half_pupil + (ix as f64 + 0.5) * da;
            for iy in 0..cells {
                let b = -case.half_pupil + (iy as f64 + 0.5) * da;
                let centre = [
                    a * (1.0 - case.distance / case.focal_length) - case.distance * (0.1 - a) / z,
                    b * (1.0 - case.distance / case.focal_length) + case.distance * b / z,
                ];
                let mask =
                    mask_probability(case, centre, 0.0, frost.slope.map(|h| case.distance * h));
                let r = ((0.1 - a).powi(2) + b * b + z * z).sqrt();
                let medium_distance = (1.5 - z) + (z - 0.5) * r / z;
                // L*cos(source)/r² = Phi*s²/(A*G) * mask * r/z³.
                // Camera and source paths each retain their own Beer length.
                sum += mask * r / z.powi(3) * (-0.6 * medium_distance).exp() * da * da * dz;
            }
        }
    }
    case.clear_flux * case.distance.powi(2)
        / (4.0 * case.half_pupil.powi(2) * case.gobo_size[0] * case.gobo_size[1])
        * frost.transmission
        * 0.4
        / (4.0 * PI)
        * sum
}

#[test]
fn g1_frost_reuses_surface_mis_and_homogeneous_haze_transport() {
    with_cx(|cx, _| {
        let mut case = ProjectorCase::new(1.0, [1, 1], vec![1.0]);
        case.half_pupil = 0.5;
        case.gobo_size = [0.2, 0.2];
        case.distance = 0.4;
        case.focal_length = 0.4 / 1.4;
        let frost = Frost {
            slope: [0.17, 0.12],
            scattered: 0.75,
            transmission: 0.82,
        };
        let mut scene = case.scene();
        scene.camera.eye = Point3::new(0.23, 0.0, 0.75);
        scene.camera.forward = Vec3::new(0.0, 0.0, 1.0);
        scene
            .primitives
            .push(floor(1.0, 2.0, LiftedSpectrum { c: [0.0; 3] }));
        let sources = [source(&case, Some(frost), &[], RigidTransform::identity())];
        let expected = 0.5 * reference_lux(&case, frost, &[], [0.23, 0.0], 1.0) / PI;
        for strategy in [
            DirectStrategy::NeeOnly,
            DirectStrategy::BsdfOnly,
            DirectStrategy::Mis,
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
            let se = value.estimated_standard_error_cd_m2.unwrap();
            // A fixed regression rule using estimated SE; no confidence or
            // finite-depth/model-error certificate follows from this check.
            assert!(
                (value.luminance_cd_m2 - expected).abs() <= 6.0 * se + 0.003 * expected,
                "{strategy:?}: observed={} expected={expected} estimated_se={se}",
                value.luminance_cd_m2
            );
            close(
                value.luminance_cd_m2,
                expected,
                0.04,
                "frosted source rho*E/pi",
            );
        }

        let case = ProjectorCase::new(3.0, [1, 1], vec![1.0]);
        let frost = Frost {
            slope: [0.07, 0.05],
            scattered: 1.0,
            transmission: 0.8,
        };
        let mut scene = case.scene();
        scene.camera.eye = Point3::new(0.1, 0.0, 2.0);
        let sources = [source(&case, Some(frost), &[], RigidTransform::identity())];
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
        let expected = haze_reference(&case, frost, 64);
        close(
            haze_reference(&case, frost, 32),
            expected,
            0.001,
            "independent frosted haze quadrature refinement",
        );
        let value = pixel(
            &scene,
            &sources,
            PathMedium::Scattering(&medium),
            Settings {
                sampler: Sampler::Iid,
                ..settings(1, 65_536)
            },
            cx,
        );
        let se = value.estimated_standard_error_cd_m2.unwrap();
        assert!(
            (value.luminance_cd_m2 - expected).abs() <= 6.0 * se + 0.001 * expected,
            "frosted haze: observed={} expected={expected} estimated_se={se}",
            value.luminance_cd_m2
        );
        close(
            value.luminance_cd_m2,
            expected,
            0.04,
            "frosted source with both Beer legs",
        );
    });
}

#[test]
fn g4_g5_frost_images_and_probes_preserve_bits_across_cancel_and_mixed_execution() {
    with_cx(|cx, gate| {
        let case = ProjectorCase::new(3.0, [2, 2], vec![1.0, 0.25, 0.5, 0.0]);
        let scene = screen_scene(&case);
        let frost = Frost {
            slope: [0.08, 0.05],
            scattered: 0.7,
            transmission: 0.9,
        };
        let sources = [source(
            &case,
            Some(frost),
            &[
                ([-0.05, 0.0], 0.2),
                ([0.0, 0.03], 0.3),
                ([0.06, -0.02], 0.25),
            ],
            RigidTransform::identity(),
        )
        .filtered(core::array::from_fn(|bin| if bin < 40 { 0.3 } else { 0.8 }))
        .unwrap()];
        let density = [0.1, 1.0, 0.5, 1.4];
        let medium = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([2, 2, 1], &density, [-3.0, -3.0, 0.2], [3.0, 3.0, 2.6]),
                MassAbsorptionSpectrum::constant(0.03).unwrap(),
                None,
                medium_budget(),
            )
            .unwrap(),
            MassScatteringSpectrum::constant(0.1).unwrap(),
            PhaseFunction::henyey_greenstein(0.35).unwrap(),
        )
        .unwrap();
        let emitters = PhotometricEmitters {
            points: &[],
            surfaces: &sources,
        };
        let probes = [-0.3, 0.0, 0.3].map(|x| {
            IlluminanceProbe::try_new(Point3::new(x, 0.1, 3.0), Vec3::new(0.0, 0.0, -1.0)).unwrap()
        });
        let new_image = || {
            PathLuminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Scattering(&medium),
                TransportOptions::default(),
                Settings {
                    width: 4,
                    height: 3,
                    sampler: Sampler::Iid,
                    ..settings(3, 32)
                },
                cx,
            )
            .unwrap()
        };
        let new_probes = || {
            PathIlluminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Scattering(&medium),
                &probes,
                TransportOptions::default(),
                IlluminanceSampling {
                    samples_per_probe: 32,
                    max_depth: 3,
                    sampler: Sampler::Iid,
                    seed: 173,
                },
                cx,
            )
            .unwrap()
        };
        let mut image_reference = new_image();
        let mut probe_reference = new_probes();
        image_reference.advance(cx, usize::MAX).unwrap();
        probe_reference.advance(cx, usize::MAX).unwrap();
        assert!(
            image_reference
                .estimates()
                .flatten()
                .any(|v| v.luminance_cd_m2 > 0.0)
        );
        assert!(
            probe_reference
                .estimates()
                .flatten()
                .all(|v| v.continuous_lux > 0.0)
        );
        let mut image = new_image();
        let mut probe_job = new_probes();
        let execution =
            RenderExecutionConfig::try_new(2, 2, 2, 64 << 20, RunId(0x66727331)).unwrap();
        image.advance(cx, 7).unwrap();
        image.advance_with_execution(cx, 13, &execution).unwrap();
        probe_job.advance_with_execution(cx, 5, &execution).unwrap();
        let image_prefix = image_bits(&image);
        let probe_prefix = probe_bits(&probe_job);
        gate.request();
        assert!(matches!(
            image.advance_with_execution(cx, 19, &execution),
            Err(PhotometryError::Cancelled)
        ));
        assert!(matches!(
            probe_job.advance(cx, 17),
            Err(PhotometryError::Cancelled)
        ));
        assert_eq!(image_bits(&image), image_prefix);
        assert_eq!(probe_bits(&probe_job), probe_prefix);
        assert_eq!(image.completed_path_samples(), 20);
        assert_eq!(probe_job.completed_path_samples(), 5);
        with_cx(|next, _| {
            let execution =
                RenderExecutionConfig::try_new(3, 1, 4, 64 << 20, RunId(0x66727332)).unwrap();
            image.advance_with_execution(next, 23, &execution).unwrap();
            image.advance(next, usize::MAX).unwrap();
            probe_job.advance(next, 7).unwrap();
            probe_job
                .advance_with_execution(next, usize::MAX, &execution)
                .unwrap();
        });
        assert!(image.is_complete() && probe_job.is_complete());
        assert_eq!(image_bits(&image), image_bits(&image_reference));
        assert_eq!(probe_bits(&probe_job), probe_bits(&probe_reference));
        assert!(
            image
                .estimates()
                .flatten()
                .all(|v| v.total_error_bound_cd_m2.is_none())
        );
        assert!(
            probe_job
                .estimates()
                .flatten()
                .all(|v| v.total_error_bound_lux.is_none())
        );
    });
}

#[test]
#[ignore = "host-specific clear and diffused 64x64 gobo, three-beam projector throughput"]
fn photometric_frost_throughput() {
    with_cx(|cx, _| {
        let mask = (0..64 * 64)
            .map(|index| match ((index % 64) / 4 + (index / 64) / 4) % 3 {
                0 => 0.0,
                1 => 0.25,
                _ => 1.0,
            })
            .collect();
        let case = ProjectorCase::new(3.0, [64, 64], mask);
        let scene = screen_scene(&case);
        let beams = [
            ([-0.06, 0.0], 0.2),
            ([0.0, 0.02], 0.4),
            ([0.06, -0.01], 0.2),
        ];
        let frost = Frost {
            slope: [0.04, 0.02],
            scattered: 0.85,
            transmission: 0.9,
        };
        let execution =
            RenderExecutionConfig::try_new(4, 4, 4, 64 << 20, RunId(0x66727333)).unwrap();
        for (label, diffusion) in [("clear", None), ("box_frost", Some(frost))] {
            let sources = [source(&case, diffusion, &beams, RigidTransform::identity())];
            let mut replay = None;
            let mut times = [0.0; 3];
            let mut mean = 0.0;
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
                    assert_eq!(&bits, expected, "{label}, repeat {repeat}");
                } else {
                    mean = job
                        .estimates()
                        .flatten()
                        .map(|v| v.luminance_cd_m2)
                        .sum::<f64>()
                        / 256.0;
                    assert!(mean > 0.0);
                    replay = Some(bits);
                }
            }
            times.sort_by(f64::total_cmp);
            eprintln!(
                "photometric_frost_throughput model={label} gobo=64x64 beams=3 workers=4 paths=16384 repeats=3 seconds_median={:.6} paths_per_second={:.1} mean_cd_m2={mean:.9} seconds_sorted={times:?}",
                times[1],
                16_384.0 / times[1]
            );
        }
    });
}
