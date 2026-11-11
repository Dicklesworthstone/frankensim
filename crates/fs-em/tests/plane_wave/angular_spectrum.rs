use super::*;
use fs_em::plane_wave::{
    angular_spectrum::{
        AngularSpectrumAperture, AngularSpectrumApertureLimits, AngularSpectrumError,
        AngularSpectrumGrid, AngularSpectrumLimits, AngularSpectrumSide as Side,
        PlanarAngularSpectrum,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_math::det;

#[path = "angular_spectrum/carrier.rs"]
mod carrier;
#[path = "angular_spectrum/aplanatic.rs"]
mod aplanatic;
#[path = "angular_spectrum/carrier_correlation.rs"]
mod carrier_correlation;
#[path = "angular_spectrum/oversampled.rs"]
mod oversampled;
#[path = "angular_spectrum/probes.rs"]
mod probes;
#[path = "angular_spectrum/transverse.rs"]
mod transverse;

fn limits() -> AngularSpectrumLimits {
    AngularSpectrumLimits {
        maximum_grid_points: 4096,
        maximum_workspace_bytes: 4 << 20,
        maximum_fft_work: 1 << 20,
    }
}
fn phase(angle: f64) -> C64 {
    C64::new(det::cos(angle), det::sin(angle))
}
fn relative(a: f64, b: f64) {
    near(a, b, 2e-11 * a.abs().max(b.abs()).max(1e-20));
}

#[test]
fn uniform_vector_beam_reduces_to_original_films_and_signed_poynting_work() {
    let incident = material(1.0, 1.0);
    let film = material(4.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&incident),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&film),
            thickness_m: 0.37,
        }],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 8,
        height: 4,
        spacing_m: [8.0, 8.0],
    };
    let input = vec![[C64::new(0.7, 0.2), C64::new(-0.1, 0.9)]; grid.width * grid.height];
    let beam = PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
    assert!(std::ptr::eq(beam.source(), &stack));
    let response = stack.response(1.0, Polarization::Te).unwrap();
    let area = grid.width as f64 * grid.height as f64 * 64.0;
    let p = 0.5 * area * (input[0][0].norm_sq() + input[0][1].norm_sq());
    relative(beam.power().incident_w, p);
    relative(
        beam.power().reflected_w,
        p * response.reflected_power_fraction,
    );
    relative(
        beam.power().transmitted_w,
        p * response.transmitted_power_fraction,
    );
    assert_eq!(beam.power().film_absorbed_w, 0.0);
    assert!(beam.power().balance_defect_w().abs() < 3e-12 * p);
    for distance in [0.0, 0.4, 7.0] {
        for side in [
            Side::Incident,
            Side::Reflected,
            Side::TotalIncidentSide,
            Side::Transmitted,
        ] {
            let field = beam.field(side, distance, || false).unwrap();
            let forward = phase(distance);
            let backward = phase(-distance) * response.reflected_tangential_e;
            let gain = match side {
                Side::Incident => forward,
                Side::Reflected => backward,
                Side::TotalIncidentSide => forward + backward,
                Side::Transmitted => phase(-distance) * response.transmitted_tangential_e,
            };
            for value in field.electric_v_per_m() {
                complex_near(value[0], input[0][0] * gain, 2e-13);
                complex_near(value[1], input[0][1] * gain, 2e-13);
                complex_near(value[2], C64::ZERO, 2e-13);
            }
            let expected = match side {
                Side::Incident => p,
                Side::Reflected => -beam.power().reflected_w,
                Side::TotalIncidentSide => p - beam.power().reflected_w,
                Side::Transmitted => beam.power().transmitted_w,
            };
            relative(field.normal_power_w(), expected);
        }
    }
}

#[test]
fn spatial_diffraction_matches_independent_direct_fourier_series() {
    // G1/G3: odd, prime and composite axes must preserve their physical
    // periods. Padding a source to the next power of two changes these modes.
    for [width, height] in [[8, 4], [5, 3], [6, 5], [7, 11]] {
        check_spatial_diffraction_against_direct_series(width, height);
    }
}

fn check_spatial_diffraction_against_direct_series(width: usize, height: usize) {
    let material = material(1.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &material,
        DielectricResponse::Lossless(&material),
        1.0,
        &[],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width,
        height,
        spacing_m: [8.0, 9.0],
    };
    let n = grid.width * grid.height;
    // Two localized, phase-shifted openings. They define the band-limited
    // periodic field, not a claim of unrepresented hard-aperture high modes.
    let mut input = vec![[C64::ZERO; 2]; n];
    input[grid.width + 2] = [C64::ONE, C64::new(0.2, 0.4)];
    input[3 * grid.width - 2] = [C64::new(0.1, 0.8), C64::new(-0.3, 0.1)];
    let beam = PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
    for distance in [0.0, 3.0, 41.0] {
        let field = beam.field(Side::Transmitted, distance, || false).unwrap();
        let mut reference_e = vec![[C64::ZERO; 3]; n];
        let mut reference_h = reference_e.clone();
        for fy in 0..grid.height {
            for fx in 0..grid.width {
                let mx = if fx <= (grid.width - 1) / 2 {
                    fx as f64
                } else {
                    fx as f64 - grid.width as f64
                };
                let my = if fy <= (grid.height - 1) / 2 {
                    fy as f64
                } else {
                    fy as f64 - grid.height as f64
                };
                let kx = core::f64::consts::TAU * mx / (grid.width as f64 * grid.spacing_m[0]);
                let ky = core::f64::consts::TAU * my / (grid.height as f64 * grid.spacing_m[1]);
                let kz = det::sqrt(1.0 - kx * kx - ky * ky);
                let mut a = [C64::ZERO; 2];
                for y in 0..grid.height {
                    for x in 0..grid.width {
                        let weight = phase(
                            kx * x as f64 * grid.spacing_m[0] + ky * y as f64 * grid.spacing_m[1],
                        )
                        .scale(1.0 / n as f64);
                        for c in 0..2 {
                            a[c] = a[c] + input[y * grid.width + x][c] * weight;
                        }
                    }
                }
                let e = [
                    a[0],
                    a[1],
                    -(a[0].scale(kx) + a[1].scale(ky)).scale(1.0 / kz),
                ];
                let h = [
                    e[2].scale(ky) - e[1].scale(kz),
                    e[0].scale(kz) - e[2].scale(kx),
                    e[1].scale(kx) - e[0].scale(ky),
                ];
                for y in 0..grid.height {
                    for x in 0..grid.width {
                        let weight = phase(
                            -(kx * x as f64 * grid.spacing_m[0]
                                + ky * y as f64 * grid.spacing_m[1]
                                + kz * distance),
                        );
                        for c in 0..3 {
                            reference_e[y * grid.width + x][c] =
                                reference_e[y * grid.width + x][c] + e[c] * weight;
                            reference_h[y * grid.width + x][c] =
                                reference_h[y * grid.width + x][c] + h[c] * weight;
                        }
                    }
                }
            }
        }
        for i in 0..n {
            for c in 0..3 {
                complex_near(field.electric_v_per_m()[i][c], reference_e[i][c], 8e-14);
                complex_near(field.magnetic_a_per_m()[i][c], reference_h[i][c], 8e-14);
            }
        }
        relative(field.normal_power_w(), beam.power().incident_w);
        if distance == 0.0 {
            for i in 0..n {
                for c in 0..2 {
                    complex_near(field.electric_v_per_m()[i][c], input[i][c], 4e-15);
                }
            }
        } else {
            assert!(field.electric_v_per_m()[0][0].abs() > 1e-6);
        }
    }
}

#[test]
fn transmitted_evanescent_beam_decays_and_carries_no_spurious_normal_power() {
    let dense = material(4.0, 1.0);
    let sparse = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&dense, DielectricResponse::Lossless(&sparse), 1.0, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 4,
        height: 1,
        spacing_m: [2.5, 1.0],
    };
    // Exact alternating samples isolate the negative Nyquist TE mode.
    let input: Vec<_> = (0..4)
        .map(|i| [C64::ZERO, C64::from_re(if i % 2 == 0 { 1.0 } else { -1.0 })])
        .collect();
    let beam = PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
    relative(beam.power().reflected_w, beam.power().incident_w);
    assert!(beam.power().transmitted_w.abs() < 1e-13);
    let near_field = beam.field(Side::Transmitted, 0.0, || false).unwrap();
    let farther = beam.field(Side::Transmitted, 2.0, || false).unwrap();
    let alpha = det::sqrt((core::f64::consts::PI / 2.5).powi(2) - 1.0);
    for i in 0..4 {
        complex_near(
            farther.electric_v_per_m()[i][1],
            near_field.electric_v_per_m()[i][1].scale(det::exp(-2.0 * alpha)),
            1e-13,
        );
    }
    assert!(near_field.normal_power_w().abs() < 1e-13);
    assert!(farther.normal_power_w().abs() < 1e-13);
}

#[test]
fn transmitted_evanescent_tail_preserves_all_maxwell_components_and_probes() {
    use fs_em::plane_wave::angular_spectrum::{
        AngularSpectrumObservation, AngularSpectrumProbeLimits,
    };

    // G1/G3: physically large incident fields remain finite at the face, and
    // all six Maxwell components have representable tails after exp(-800)
    // itself has underflowed. The source is lossless, so decay creates no heat.
    let dense = material(4.0, 1.0);
    let sparse = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&dense, DielectricResponse::Lossless(&sparse), 1.0, &[])
            .unwrap();
    let beam = PlanarAngularSpectrum::prepare_modes_with_carrier(
        &stack,
        AngularSpectrumGrid {
            width: 1,
            height: 1,
            spacing_m: [1.0, 1.0],
        },
        [1.5, 0.0],
        &[[C64::new(1e100, -2e100), C64::new(-0.5e100, 3e100)]],
        limits(),
        || false,
    )
    .unwrap();
    let face = beam.field(Side::Transmitted, 0.0, || false).unwrap();
    let cosine = (1.0_f64 - 0.75 * 0.75).sqrt();
    let alpha = -stack
        .response(cosine, Polarization::Te)
        .unwrap()
        .exit_normal_wavenumber_rad_m
        .im;
    for (exposure, surviving) in [
        (720.0, 2.032_230_802_424_293_2e-213),
        (800.0, 3.667_874_584_177_687e-248),
    ] {
        let depth = exposure / alpha;
        let field = beam.field(Side::Transmitted, depth, || false).unwrap();
        let probe = beam
            .probe_fields(
                &[AngularSpectrumObservation {
                    side: Side::Transmitted,
                    distance_from_face_m: depth,
                    position_m: [0.0, 0.0],
                }],
                AngularSpectrumProbeLimits {
                    maximum_points: 1,
                    maximum_mode_evaluations: 1,
                    maximum_workspace_bytes: 4 << 20,
                },
                || false,
            )
            .unwrap();
        for (actual, original) in field.electric_v_per_m()[0]
            .iter()
            .chain(&field.magnetic_a_per_m()[0])
            .zip(
                face.electric_v_per_m()[0]
                    .iter()
                    .chain(&face.magnetic_a_per_m()[0]),
            )
        {
            for (actual, original) in [actual.re, actual.im]
                .into_iter()
                .zip([original.re, original.im])
            {
                let expected = (original / 1e100) * surviving;
                if expected == 0.0 {
                    assert_eq!(actual, 0.0);
                } else {
                    assert!((actual / expected - 1.0).abs() < 5e-13);
                }
            }
        }
        assert_eq!(probe[0].electric_v_per_m, field.electric_v_per_m()[0]);
        assert_eq!(probe[0].magnetic_a_per_m, field.magnetic_a_per_m()[0]);
    }
}

#[test]
fn critical_refraction_preserves_longitudinal_tm_and_finite_te_fields_through_films() {
    // G1: this negative Nyquist mode's admitted cosine rounds to exactly 3/4.
    // Both it and epsilon=7/16 are dyadic, so the outgoing discriminant is
    // exactly zero without relying on a rounded product to manufacture it.
    let cosine = 0.75;
    let epsilon_exit = 0.4375;
    let grid = AngularSpectrumGrid {
        width: 4,
        height: 1,
        spacing_m: [core::f64::consts::PI / det::sqrt(epsilon_exit), 1.0],
    };
    let sine = core::f64::consts::PI / grid.spacing_m[0];
    let incident = material(1.0, 1.0);
    let exit = material(epsilon_exit, 1.0);
    let film = material(1.7, 1.0);
    let amplitude = [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)];
    let input: Vec<_> = (0..4)
        .map(|i| amplitude.map(|a| a.scale(if i % 2 == 0 { 1.0 } else { -1.0 })))
        .collect();
    for thickness in [0.0, 0.43] {
        let stack = PlanarDielectricLayers::try_new(
            &incident,
            DielectricResponse::Lossless(&exit),
            1.0,
            &[DielectricLayer {
                response: DielectricResponse::Lossless(&film),
                thickness_m: thickness,
            }],
        )
        .unwrap();
        let te = stack.response(cosine, Polarization::Te).unwrap();
        let tm = stack.response(cosine, Polarization::Tm).unwrap();
        assert_eq!(te.exit_normal_wavenumber_rad_m, C64::ZERO);
        assert_eq!(tm.exit_normal_wavenumber_rad_m, C64::ZERO);
        assert_eq!(te.transmitted_conjugate_h, C64::ZERO);
        assert_eq!(tm.transmitted_tangential_e, C64::ZERO);
        // Independent one-film transfer matrix with the critical load:
        // TE has H_exit=0; TM has E_exit=0. These retain film phase exactly.
        let film_kz = (1.7 - sine * sine).sqrt();
        let phase = film_kz * thickness;
        let te_gain = C64::from_re(2.0) / C64::new(phase.cos(), phase.sin() * film_kz / cosine);
        let tm_h_gain =
            C64::from_re(2.0) / C64::new(phase.cos(), phase.sin() * film_kz / (1.7 * cosine));
        complex_near(te.transmitted_tangential_e, te_gain, 3e-14);
        complex_near(tm.transmitted_conjugate_h, tm_h_gain, 3e-14);
        let beam =
            PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
        relative(beam.power().reflected_w, beam.power().incident_w);
        assert_eq!(beam.power().transmitted_w, 0.0);
        let face = beam.field(Side::Transmitted, 0.0, || false).unwrap();
        for distance in [0.0, 2.0, 1e6] {
            let field = beam.field(Side::Transmitted, distance, || false).unwrap();
            assert!(field.normal_power_w().abs() < 1e-13);
            for (i, &sample) in input.iter().enumerate() {
                let e = field.electric_v_per_m()[i];
                let h = field.magnetic_a_per_m()[i];
                complex_near(e[0], C64::ZERO, 1e-14);
                complex_near(e[1], sample[1] * te_gain, 2e-13);
                complex_near(
                    e[2],
                    (sample[0] * tm_h_gain).scale(sine / (cosine * epsilon_exit)),
                    2e-13,
                );
                complex_near(h[0], C64::ZERO, 1e-14);
                complex_near(h[1], (sample[0] * tm_h_gain).scale(1.0 / cosine), 2e-13);
                complex_near(h[2], (sample[1] * te_gain).scale(-sine), 2e-13);
                // Normal flux is zero, but the critical wave still transports
                // finite energy parallel to the interface along negative x.
                let tangential_flux = 0.5 * (e[1] * h[2].conj() - e[2] * h[1].conj()).re;
                assert!(tangential_flux < -0.1);
                assert_eq!(e, face.electric_v_per_m()[i]);
                assert_eq!(h, face.magnetic_a_per_m()[i]);
            }
        }
        if thickness == 0.0 {
            let entrance = beam.field(Side::TotalIncidentSide, 0.0, || false).unwrap();
            for i in 0..4 {
                for c in 0..2 {
                    complex_near(
                        entrance.electric_v_per_m()[i][c],
                        face.electric_v_per_m()[i][c],
                        2e-13,
                    );
                    complex_near(
                        entrance.magnetic_a_per_m()[i][c],
                        face.magnetic_a_per_m()[i][c],
                        2e-13,
                    );
                }
                complex_near(
                    entrance.electric_v_per_m()[i][2],
                    face.electric_v_per_m()[i][2].scale(epsilon_exit),
                    2e-13,
                );
                complex_near(
                    entrance.magnetic_a_per_m()[i][2],
                    face.magnetic_a_per_m()[i][2],
                    2e-13,
                );
            }
        }
    }
}

#[test]
fn propagating_and_evanescent_transmitted_fields_converge_to_the_critical_wave() {
    // G3: both physical one-sided limits approach the same finite vector field.
    // The admitted cosine is exactly 3/4 and its critical contrast is 7/16.
    let critical = 0.4375;
    let grid = AngularSpectrumGrid {
        width: 2,
        height: 1,
        spacing_m: [core::f64::consts::PI / det::sqrt(critical), 1.0],
    };
    let incident = material(1.0, 1.0);
    let input = [
        [C64::ONE, C64::new(0.3, 0.2)],
        [C64::from_re(-1.0), C64::new(-0.3, -0.2)],
    ];
    let field = |contrast: f64| {
        let exit = material(contrast, 1.0);
        let stack = PlanarDielectricLayers::try_new(
            &incident,
            DielectricResponse::Lossless(&exit),
            1.0,
            &[],
        )
        .unwrap();
        PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false)
            .unwrap()
            .field(Side::Transmitted, 0.7, || false)
            .unwrap()
    };
    let reference = field(critical);
    for side in [-1.0, 1.0] {
        let mut previous = f64::INFINITY;
        for delta in [1e-6, 1e-8, 1e-10, 1e-12] {
            let nearby = field(critical + side * delta);
            let error = nearby
                .electric_v_per_m()
                .iter()
                .flatten()
                .chain(nearby.magnetic_a_per_m().iter().flatten())
                .zip(
                    reference
                        .electric_v_per_m()
                        .iter()
                        .flatten()
                        .chain(reference.magnetic_a_per_m().iter().flatten()),
                )
                .map(|(&actual, &expected)| (actual - expected).abs())
                .fold(0.0, f64::max);
            assert!(
                error < previous * 0.2,
                "critical-limit error {error}, previous {previous}"
            );
            previous = error;
        }
        assert!(previous < 1e-4);
    }
}

#[test]
fn original_debye_films_and_exit_keep_absorption_separate_from_propagation_loss() {
    let incident = material(E, M);
    let state = source(
        &[
            (EPS_NAME, EPS, 2.0 * E),
            (MU_NAME, MU, M),
            ("delta", EPS, E),
            ("tau", Time::DIMS, 1e-9),
        ],
        ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        QueryPoint::new().with("T", 300.0).unwrap(),
    );
    let dielectric = ResolvedDebyeCapacitor::try_from_material(
        state,
        EPS_NAME,
        &[DebyePoleProperties {
            permittivity_increment: "delta".into(),
            relaxation_time: "tau".into(),
        }],
        0.01,
        0.001,
    )
    .unwrap();
    let stack = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Debye(&dielectric),
        1e9,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&dielectric),
            thickness_m: 0.2,
        }],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 4,
        height: 4,
        spacing_m: [2.0, 2.0],
    };
    let input = vec![[C64::ONE, C64::new(0.2, 0.7)]; 16];
    let beam = PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
    assert!(beam.power().film_absorbed_w > 0.0);
    assert!(beam.power().balance_defect_w().abs() < 1e-12 * beam.power().incident_w);
    let exit = beam.field(Side::Transmitted, 0.0, || false).unwrap();
    let later = beam.field(Side::Transmitted, 0.3, || false).unwrap();
    relative(exit.normal_power_w(), beam.power().transmitted_w);
    let kz = stack
        .response(1.0, Polarization::Te)
        .unwrap()
        .exit_normal_wavenumber_rad_m;
    relative(
        later.normal_power_w(),
        exit.normal_power_w() * det::exp(2.0 * kz.im * 0.3),
    );
    assert!(later.normal_power_w() < exit.normal_power_w());
}

#[test]
fn budgets_light_cone_and_cancellation_refuse_without_mutating_the_beam() {
    let material = material(1.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &material,
        DielectricResponse::Lossless(&material),
        1.0,
        &[],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 4,
        height: 4,
        spacing_m: [8.0, 8.0],
    };
    let input = vec![[C64::ONE, C64::ZERO]; 16];
    let beam = PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
    let before = beam.field(Side::Transmitted, 0.4, || false).unwrap();
    assert!(matches!(
        beam.field(Side::Transmitted, 0.4, || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    let after = beam.field(Side::Transmitted, 0.4, || false).unwrap();
    assert_eq!(before.electric_v_per_m(), after.electric_v_per_m());
    for constrained in [
        AngularSpectrumLimits {
            maximum_grid_points: 15,
            ..limits()
        },
        AngularSpectrumLimits {
            maximum_workspace_bytes: beam.workspace_bytes() - 1,
            ..limits()
        },
        AngularSpectrumLimits {
            maximum_fft_work: beam.field_fft_work() - 1,
            ..limits()
        },
    ] {
        assert!(
            PlanarAngularSpectrum::prepare(&stack, grid, &input, constrained, || false).is_err()
        );
    }
    let unsupported: Vec<_> = (0..16)
        .map(|i| [C64::from_re(if i % 2 == 0 { 1.0 } else { -1.0 }), C64::ZERO])
        .collect();
    assert!(matches!(
        PlanarAngularSpectrum::prepare(
            &stack,
            AngularSpectrumGrid {
                spacing_m: [1.0, 1.0],
                ..grid
            },
            &unsupported,
            limits(),
            || false
        ),
        Err(AngularSpectrumError::IncidentLightCone { .. })
    ));
    assert!(
        PlanarAngularSpectrum::prepare(&stack, grid, &input[..15], limits(), || false).is_err()
    );
    assert!(beam.field(Side::Transmitted, -0.1, || false).is_err());
    let mut polls = 0;
    assert!(matches!(
        beam.field(Side::Transmitted, 0.4, || {
            polls += 1;
            polls > 5
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
}

#[test]
fn finite_detector_resolves_original_beam_fringe_width_position_and_weak_dark_flux() {
    let material = material(1.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &material,
        DielectricResponse::Lossless(&material),
        1.0,
        &[],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 4,
        height: 1,
        spacing_m: [8.0; 2],
    };
    let q = core::f64::consts::TAU / 32.0;
    let kz = det::sqrt(1.0 - q * q);
    let bounds = AngularSpectrumApertureLimits {
        maximum_mode_pairs: 20,
        maximum_workspace_bytes: 1 << 20,
    };
    for sign in [1.0, -1.0] {
        let samples: Vec<_> = (0..4)
            .map(|i| {
                let x = i as f64 * 8.0;
                [C64::ZERO, phase(-q * x) + phase(q * x).scale(sign)]
            })
            .collect();
        let beam =
            PlanarAngularSpectrum::prepare(&stack, grid, &samples, limits(), || false).unwrap();
        for (center, width) in [(0.0, 0.013), (1.3, 5.0), (0.0, 16.0)] {
            let detector =
                AngularSpectrumAperture::try_new([center - width / 2.0, -0.4], [width, 2.0])
                    .unwrap();
            let actual = beam
                .measure_aperture(Side::Transmitted, 0.37, detector, bounds, || false)
                .unwrap();
            let expected =
                kz * (1.0 + sign * (2.0 * q * center).cos() * (q * width).sin() / (q * width));
            near(actual.mean_poynting_w_m2[2], expected, 3e-14);
            near(actual.normal_power_w, expected * detector.area_m2(), 1e-12);
        }
        if sign < 0.0 {
            let width = 1e-7;
            let detector =
                AngularSpectrumAperture::try_new([-width / 2.0, 0.0], [width, 2.0]).unwrap();
            let actual = beam
                .measure_aperture(Side::Incident, 0.0, detector, bounds, || false)
                .unwrap();
            let expected = kz * (q * width).powi(2) / 6.0;
            near(actual.mean_poynting_w_m2[2], expected, expected * 3e-7);
            assert!(actual.normal_power_w > 0.0);
        }
    }
}

#[test]
fn finite_detector_keeps_original_reflection_tir_and_whole_cell_power() {
    let incident = material(1.0, 1.0);
    let exit = material(0.04, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&incident, DielectricResponse::Lossless(&exit), 1.0, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 4,
        height: 2,
        spacing_m: [8.0; 2],
    };
    let samples: Vec<_> = (0..8)
        .map(|i| {
            [
                C64::new(0.3 + i as f64 * 0.1, -0.2),
                C64::new(-0.7, 0.1 * i as f64),
            ]
        })
        .collect();
    let beam = PlanarAngularSpectrum::prepare(&stack, grid, &samples, limits(), || false).unwrap();
    let bounds = AngularSpectrumApertureLimits {
        maximum_mode_pairs: 48,
        maximum_workspace_bytes: 1 << 20,
    };
    let detector = AngularSpectrumAperture::try_new([-3.0, 0.7], [32.0, 16.0]).unwrap();
    for side in [
        Side::Incident,
        Side::Reflected,
        Side::TotalIncidentSide,
        Side::Transmitted,
    ] {
        let sampled = beam.field(side, 0.63, || false).unwrap();
        let observed = beam
            .measure_aperture(side, 0.63, detector, bounds, || false)
            .unwrap();
        relative(observed.normal_power_w, sampled.normal_power_w());
        if side == Side::Reflected {
            assert!(observed.normal_power_w < 0.0);
        }
    }
    assert!(matches!(
        beam.measure_aperture(Side::Transmitted, 0.0, detector, bounds, || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(
        beam.measure_aperture(Side::Transmitted, -1.0, detector, bounds, || false)
            .is_err()
    );
    assert!(
        beam.measure_aperture(
            Side::Transmitted,
            0.0,
            detector,
            AngularSpectrumApertureLimits {
                maximum_mode_pairs: 47,
                ..bounds
            },
            || false
        )
        .is_err()
    );
    let before = beam
        .measure_aperture(Side::Transmitted, 0.0, detector, bounds, || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        beam.measure_aperture(Side::Transmitted, 0.0, detector, bounds, || {
            polls += 1;
            polls > 20
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert_eq!(
        beam.measure_aperture(Side::Transmitted, 0.0, detector, bounds, || false)
            .unwrap(),
        before
    );
}
