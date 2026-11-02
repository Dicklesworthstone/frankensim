//! Joint-face, finite-window interference through the native two-port owner.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumIlluminationLimits, AngularSpectrumIlluminationSide as Face,
    PreparedAngularSpectrumIllumination,
};

fn illuminated<'b, 's, 'm>(
    left: &'b PlanarAngularSpectrum<'s, 'm>,
    right: [C64; 2],
) -> PreparedAngularSpectrumIllumination<'b, 's, 'm> {
    left.prepare_illumination(
        &[right; 4],
        AngularSpectrumIlluminationLimits {
            maximum_layer_mode_pairs: 4096,
            maximum_workspace_bytes: 8 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap()
}

fn phase(angle: f64) -> C64 {
    let (s, c) = angle.sin_cos();
    C64::new(c, s)
}

fn local(field: &AngularSpectrumLayerField) -> [C64; 6] {
    let e = field.field().electric_v_per_m()[0];
    let h = field.field().magnetic_a_per_m()[0];
    [e[0], e[1], e[2], h[0], h[1], h[2]]
}

fn moments(point: &fs_em::plane_wave::angular_spectrum::AngularSpectrumExposurePoint) -> [f64; 9] {
    let e = point.electric_mean_squared_v2_m2;
    let h = point.magnetic_mean_squared_a2_m2;
    let p = point.poynting_w_m2;
    [e[0], e[1], e[2], h[0], h[1], h[2], p[0], p[1], p[2]]
}

// Independent real-signal Simpson integration, with an explicitly enumerated
// proper phase amplitude. The last three channels use the original fixture
// ODE tau P_dot+P=Delta E, eps_inf=2, mu=1, Delta=.7, tau=.4, sigma=.3.
// They never use a production temporal kernel or normalized loss transfer.
fn time_reference(
    lines: &[([C64; 6], [C64; 6], f64)],
    window: HarmonicWindow,
    phases: &[C64],
) -> [f64; 12] {
    let steps = 2048;
    let mut sums = [0.0; 12];
    for &z in phases {
        for step in 0..=steps {
            let t = window.start_s() + window.duration_s() * step as f64 / steps as f64;
            let mut f = [0.0; 6];
            let mut df = [0.0; 6];
            let mut p = [0.0; 3];
            let mut dp = [0.0; 3];
            for (mean, factor, omega) in lines {
                let carrier = phase(omega * t);
                for c in 0..6 {
                    let amplitude = mean[c] + z * factor[c];
                    f[c] += (amplitude * carrier).re;
                    df[c] += (C64::new(-amplitude.im, amplitude.re).scale(*omega) * carrier).re;
                    if c < 3 {
                        let polarization = amplitude.scale(0.7) / C64::new(1.0, omega * 0.4);
                        p[c] += (polarization * carrier).re;
                        dp[c] += (C64::new(-polarization.im, polarization.re).scale(*omega)
                            * carrier)
                            .re;
                    }
                }
            }
            let weight = if step == 0 || step == steps {
                1.0
            } else if step % 2 == 0 {
                2.0
            } else {
                4.0
            };
            for c in 0..3 {
                sums[c] += weight * f[c] * f[c];
                sums[c + 3] += weight * f[c + 3] * f[c + 3];
                let (a, b) = ((c + 1) % 3, (c + 2) % 3);
                sums[c + 6] += weight * (f[a] * f[b + 3] - f[b] * f[a + 3]);
                sums[9] += weight * 0.3 * f[c] * f[c];
                sums[10] += weight * (0.4 / 0.7) * dp[c] * dp[c];
                sums[11] +=
                    weight * (2.0 * f[c] * df[c] + f[c + 3] * df[c + 3] + p[c] * dp[c] / 0.7);
            }
        }
    }
    sums.map(|v| v / (3.0 * steps as f64 * phases.len() as f64))
}

// Independent homogeneous epsilon=mu=1 Maxwell solution. The right input is
// phased at z=d; both electric inputs use the same Cartesian x/y axes.
fn transparent_field(
    face: Face,
    distance: f64,
    thickness: f64,
    omega: f64,
    left: [C64; 2],
    right: [C64; 2],
) -> [C64; 6] {
    let z = if matches!(
        face,
        Face::IncomingLeft | Face::OutgoingLeft | Face::TotalLeft
    ) {
        -distance
    } else {
        thickness + distance
    };
    let forward = matches!(
        face,
        Face::IncomingLeft | Face::OutgoingRight | Face::TotalLeft | Face::TotalRight
    );
    let backward = matches!(
        face,
        Face::OutgoingLeft | Face::IncomingRight | Face::TotalLeft | Face::TotalRight
    );
    let l = if forward {
        left.map(|v| v * phase(-omega * z))
    } else {
        [C64::ZERO; 2]
    };
    let r = if backward {
        right.map(|v| v * phase(omega * (z - thickness)))
    } else {
        [C64::ZERO; 2]
    };
    [
        l[0] + r[0],
        l[1] + r[1],
        C64::ZERO,
        -l[1] + r[1],
        l[0] - r[0],
        C64::ZERO,
    ]
}

#[test]
fn g1_counterpropagating_two_frequency_fields_match_real_maxwell_waves_on_both_faces() {
    let medium = material(1.0, 1.0);
    let frequencies = [2.0, 2.3];
    let thickness = 0.31;
    let stacks = frequencies.map(|omega| {
        PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&medium),
            omega,
            &[DielectricLayer {
                response: DielectricResponse::Lossless(&medium),
                thickness_m: thickness,
            }],
        )
        .unwrap()
    });
    let l = [
        [C64::new(0.8, 0.2), C64::new(-0.1, 0.3)],
        [C64::new(-0.3, 0.6), C64::new(0.4, 0.1)],
    ];
    let r = [
        [C64::new(-0.2, 0.4), C64::new(0.3, -0.2)],
        [C64::new(0.7, -0.1), C64::new(-0.1, 0.2)],
    ];
    let left = [constant(&stacks[0], l[0]), constant(&stacks[1], l[1])];
    let lines = [illuminated(&left[0], r[0]), illuminated(&left[1], r[1])];
    let bank = AngularSpectrumExposure::try_new_illumination(&left[0], &lines, limits(), || false)
        .unwrap();
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    for face in [
        Face::IncomingLeft,
        Face::OutgoingLeft,
        Face::TotalLeft,
        Face::IncomingRight,
        Face::OutgoingRight,
        Face::TotalRight,
    ] {
        let reference = (0..2)
            .map(|i| {
                (
                    transparent_field(face, 0.23, thickness, frequencies[i], l[i], r[i]),
                    [C64::ZERO; 6],
                    frequencies[i],
                )
            })
            .collect::<Vec<_>>();
        let expected = time_reference(&reference, window, &[C64::ONE]);
        let observed = bank
            .illumination_field(face, 0.23, window, || false)
            .unwrap();
        for point in observed.points() {
            for (actual, expected) in moments(point).into_iter().zip(expected) {
                near(actual, expected, 3e-12 * expected.abs().max(1.0));
            }
        }
        if face == Face::IncomingRight {
            assert!(observed.normal_power_w() < 0.0);
        }
    }
    for (side, face) in [
        (Side::Incident, Face::IncomingLeft),
        (Side::Reflected, Face::OutgoingLeft),
        (Side::TotalIncidentSide, Face::TotalLeft),
        (Side::Transmitted, Face::OutgoingRight),
    ] {
        assert_eq!(
            bank.field(side, 0.23, window, || false).unwrap().points(),
            bank.illumination_field(face, 0.23, window, || false)
                .unwrap()
                .points()
        );
    }
}

#[test]
fn g1_two_sided_dispersive_film_heat_and_storage_close_the_finite_window_balance() {
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let response = DielectricResponse::DebyeWithOhmic {
        dielectric: &dielectric,
        conductor: &conductor,
    };
    let frequencies = [2.0, 2.3];
    let stacks = frequencies.map(|omega| {
        PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&medium),
            omega,
            &[DielectricLayer {
                response,
                thickness_m: 0.31,
            }],
        )
        .unwrap()
    });
    let left = [
        constant(&stacks[0], [C64::new(0.8, 0.2), C64::ZERO]),
        constant(&stacks[1], [C64::new(-0.3, 0.6), C64::ZERO]),
    ];
    let lines = [
        illuminated(&left[0], [C64::new(-0.2, 0.4), C64::ZERO]),
        illuminated(&left[1], [C64::new(0.7, -0.1), C64::ZERO]),
    ];
    let bank = AngularSpectrumExposure::try_new_illumination(&left[0], &lines, limits(), || false)
        .unwrap();
    let prepared = bank.prepare_layer(0, heat_limits(), || false).unwrap();
    assert_eq!(prepared.frequency_count(), 2);
    assert_eq!(
        prepared.preparation_layer_mode_pairs(),
        0,
        "native two-port caches must be reused"
    );
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    let phasors = lines
        .iter()
        .zip(frequencies)
        .map(|(line, omega)| {
            (
                local(&line.layer_field(0, 0.155, || false).unwrap()),
                [C64::ZERO; 6],
                omega,
            )
        })
        .collect::<Vec<_>>();
    let expected = time_reference(&phasors, window, &[C64::ONE]);
    let observed = prepared.field(0.155, window, || false).unwrap();
    for (point, heat) in observed.field().points().iter().zip(observed.dissipation()) {
        for (actual, expected) in moments(point).into_iter().zip(expected) {
            near(actual, expected, 3e-12 * expected.abs().max(1.0));
        }
        for (actual, expected) in [
            heat.ohmic_heat_w_m3,
            heat.polarization_heat_w_m3,
            heat.storage_rate_w_m3,
        ]
        .into_iter()
        .zip(&expected[9..])
        {
            near(actual, *expected, 3e-11 * expected.abs().max(1.0));
        }
        assert!(heat.total_heat_w_m3 > 0.0);
        // Independent two-interface solution plus 2048-step real-time Simpson
        // quadrature; these anchors do not reuse native two-port fields.
        near(heat.total_heat_w_m3, 0.25445485203009016, 2e-11);
        near(heat.storage_rate_w_m3, -2.600697727169648, 2e-10);
    }
    let mut cancelled_polls = 0;
    assert!(matches!(
        prepared.field(0.155, window, || {
            cancelled_polls += 1;
            cancelled_polls > 6
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert_eq!(
        observed.dissipation(),
        prepared
            .field(0.155, window, || false)
            .unwrap()
            .dissipation()
    );
    let intervals = 256;
    let dz = 0.31 / intervals as f64;
    let mut rate = 0.0;
    for i in 0..=intervals {
        let z = if i == intervals { 0.31 } else { i as f64 * dz };
        let field = prepared.field(z, window, || false).unwrap();
        let weight = if i == 0 || i == intervals {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        rate +=
            weight * (field.absorbed_power_per_depth_w_m() + field.storage_rate_per_depth_w_m());
    }
    rate *= dz / 3.0;
    let boundary = bank
        .illumination_field(Face::TotalLeft, 0.0, window, || false)
        .unwrap()
        .normal_power_w()
        - bank
            .illumination_field(Face::TotalRight, 0.0, window, || false)
            .unwrap()
            .normal_power_w();
    near(rate, boundary, 2e-9 * boundary.abs().max(1.0));
}

#[test]
fn g3_joint_face_factors_add_to_one_sided_mean_without_losing_cross_face_covariance() {
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let response = DielectricResponse::DebyeWithOhmic {
        dielectric: &dielectric,
        conductor: &conductor,
    };
    let frequencies = [2.0, 2.3];
    let stacks = frequencies.map(|omega| {
        PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&medium),
            omega,
            &[DielectricLayer {
                response,
                thickness_m: 0.31,
            }],
        )
        .unwrap()
    });
    let mean_lines = [
        constant(&stacks[0], [C64::new(0.1, -0.2), C64::ZERO]),
        constant(&stacks[1], [C64::new(0.2, 0.1), C64::ZERO]),
    ];
    let left = [
        constant(&stacks[0], [C64::new(0.8, 0.2), C64::ZERO]),
        constant(&stacks[1], [C64::new(-0.3, 0.6), C64::ZERO]),
    ];
    let lines = [
        illuminated(&left[0], [C64::new(-0.2, 0.4), C64::ZERO]),
        illuminated(&left[1], [C64::new(0.7, -0.1), C64::ZERO]),
    ];
    let mean =
        AngularSpectrumExposure::try_new(&mean_lines[0], &mean_lines, limits(), || false).unwrap();
    let factors =
        [
            AngularSpectrumExposure::try_new_illumination(&left[0], &lines, limits(), || false)
                .unwrap(),
        ];
    let ensemble = AngularSpectrumExposureEnsemble::try_new(
        &left[0],
        Some(&mean),
        &factors,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let prepared = ensemble.prepare_layer(0, heat_limits(), || false).unwrap();
    assert_eq!(
        prepared.preparation_layer_mode_pairs(),
        16,
        "only the two mean frequencies need a cache"
    );
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    let phasors = (0..2)
        .map(|i| {
            let mean_field = mean_lines[i]
                .prepare_layer(
                    0,
                    AngularSpectrumLayerLimits {
                        maximum_layer_mode_pairs: 4096,
                        maximum_workspace_bytes: 8 << 20,
                        maximum_fft_work: 1 << 20,
                    },
                    || false,
                )
                .unwrap()
                .field(0.155, || false)
                .unwrap();
            (
                local(&mean_field),
                local(&lines[i].layer_field(0, 0.155, || false).unwrap()),
                frequencies[i],
            )
        })
        .collect::<Vec<_>>();
    let phases = [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)];
    let expected = time_reference(&phasors, window, &phases);
    let observed = prepared.field(0.155, window, || false).unwrap();
    for (point, heat) in observed.field().points().iter().zip(observed.dissipation()) {
        for (actual, expected) in moments(point).into_iter().zip(expected) {
            near(actual, expected, 4e-12 * expected.abs().max(1.0));
        }
        for (actual, expected) in [
            heat.ohmic_heat_w_m3,
            heat.polarization_heat_w_m3,
            heat.storage_rate_w_m3,
        ]
        .into_iter()
        .zip(&expected[9..])
        {
            near(actual, *expected, 4e-11 * expected.abs().max(1.0));
        }
    }
    let exterior = ensemble
        .illumination_field(Face::TotalRight, 0.23, window, || false)
        .unwrap();
    let exterior_lines = (0..2)
        .map(|i| {
            let m = mean_lines[i]
                .field(Side::Transmitted, 0.23, || false)
                .unwrap();
            let f = lines[i].field(Face::TotalRight, 0.23, || false).unwrap();
            let vector = |e: &[C64; 3], h: &[C64; 3]| [e[0], e[1], e[2], h[0], h[1], h[2]];
            (
                vector(&m.electric_v_per_m()[0], &m.magnetic_a_per_m()[0]),
                vector(&f.electric_v_per_m()[0], &f.magnetic_a_per_m()[0]),
                frequencies[i],
            )
        })
        .collect::<Vec<_>>();
    let expected = time_reference(&exterior_lines, window, &phases);
    for point in exterior.points() {
        for (actual, expected) in moments(point).into_iter().zip(expected) {
            near(actual, expected, 4e-12 * expected.abs().max(1.0));
        }
    }
}

#[test]
fn g3_equal_frequency_joint_sources_cancel_before_fft_but_independent_factors_do_not() {
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let response = DielectricResponse::DebyeWithOhmic {
        dielectric: &dielectric,
        conductor: &conductor,
    };
    let stack = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&medium),
        2.0,
        &[DielectricLayer {
            response,
            thickness_m: 0.31,
        }],
    )
    .unwrap();
    let amplitude = [C64::new(0.8, 0.2), C64::new(-0.1, 0.3)];
    let right = [C64::new(-0.2, 0.4), C64::new(0.2, 0.1)];
    let left = [
        constant(&stack, amplitude),
        constant(&stack, amplitude.map(|v| -v)),
    ];
    let lines = [
        illuminated(&left[0], right),
        illuminated(&left[1], right.map(|v| -v)),
    ];
    let bank = AngularSpectrumExposure::try_new_illumination(&left[0], &lines, limits(), || false)
        .unwrap();
    assert_eq!(bank.frequency_count(), 1);
    assert_eq!(bank.field_fft_work(), left[0].field_fft_work());
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    for face in [Face::TotalLeft, Face::TotalRight] {
        let field = bank
            .illumination_field(face, 0.0, window, || false)
            .unwrap();
        assert!(field.points().iter().all(|p| moments(p) == [0.0; 9]));
    }
    let heat = bank
        .prepare_layer(0, heat_limits(), || false)
        .unwrap()
        .field(0.155, window, || false)
        .unwrap();
    assert!(heat.dissipation().iter().all(|p| p.total_heat_w_m3 == 0.0));
    let factors = lines
        .iter()
        .map(|line| {
            AngularSpectrumExposure::try_new_illumination(
                &left[0],
                std::slice::from_ref(line),
                limits(),
                || false,
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    let ensemble = AngularSpectrumExposureEnsemble::try_new(
        &left[0],
        None,
        &factors,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    assert!(
        ensemble
            .prepare_layer(0, heat_limits(), || false)
            .unwrap()
            .field(0.155, window, || false)
            .unwrap()
            .dissipation()[0]
            .total_heat_w_m3
            > 0.0
    );
    let dark =
        AngularSpectrumExposure::try_new_illumination(&left[0], &[], limits(), || false).unwrap();
    assert_eq!(
        dark.illumination_field(Face::TotalRight, 0.0, window, || false)
            .unwrap()
            .normal_power_w(),
        0.0
    );
}

#[test]
fn g4_two_sided_exposure_retains_source_identity_budgets_and_atomic_cancellation() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 2.0, &[])
            .unwrap();
    let left = constant(&stack, [C64::ONE, C64::ZERO]);
    let lines = [illuminated(&left, [C64::new(0.2, 0.3), C64::ZERO])];
    let bank =
        AngularSpectrumExposure::try_new_illumination(&left, &lines, limits(), || false).unwrap();
    let exact = AngularSpectrumExposureLimits {
        maximum_lines: 1,
        maximum_point_line_pairs: 4,
        maximum_point_frequency_pairs: bank.point_frequency_pairs(),
        maximum_workspace_bytes: bank.workspace_bytes(),
        maximum_fft_work: bank.field_fft_work(),
    };
    assert!(AngularSpectrumExposure::try_new_illumination(&left, &lines, exact, || false).is_ok());
    for cap in [
        AngularSpectrumExposureLimits {
            maximum_lines: 0,
            ..exact
        },
        AngularSpectrumExposureLimits {
            maximum_point_line_pairs: 3,
            ..exact
        },
        AngularSpectrumExposureLimits {
            maximum_point_frequency_pairs: exact.maximum_point_frequency_pairs - 1,
            ..exact
        },
        AngularSpectrumExposureLimits {
            maximum_workspace_bytes: exact.maximum_workspace_bytes - 1,
            ..exact
        },
        AngularSpectrumExposureLimits {
            maximum_fft_work: exact.maximum_fft_work - 1,
            ..exact
        },
    ] {
        assert!(
            AngularSpectrumExposure::try_new_illumination(&left, &lines, cap, || false).is_err()
        );
    }
    assert!(matches!(
        AngularSpectrumExposure::try_new_illumination(&left, &lines, limits(), || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    let replacement = material(1.0, 1.0);
    let other_stack = PlanarDielectricLayers::try_new(
        &replacement,
        DielectricResponse::Lossless(&replacement),
        2.0,
        &[],
    )
    .unwrap();
    let other_left = constant(&other_stack, [C64::ONE, C64::ZERO]);
    let wrong = [illuminated(&other_left, [C64::ZERO; 2])];
    assert!(
        AngularSpectrumExposure::try_new_illumination(&left, &wrong, limits(), || false).is_err()
    );
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    for distance in [-1.0, f64::NAN] {
        assert!(
            bank.illumination_field(Face::TotalRight, distance, window, || false)
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        bank.illumination_field(Face::TotalRight, 0.0, window, || {
            polls += 1;
            polls > 4
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    let retry = bank
        .illumination_field(Face::TotalRight, 0.0, window, || false)
        .unwrap();
    assert_eq!(
        retry.points(),
        bank.illumination_field(Face::TotalRight, 0.0, window, || false)
            .unwrap()
            .points()
    );
    let one = [constant(&stack, [C64::ONE, C64::ZERO])];
    let original = AngularSpectrumExposure::try_new(&one[0], &one, limits(), || false).unwrap();
    assert_eq!(
        original
            .illumination_field(Face::IncomingRight, 0.0, window, || false)
            .unwrap()
            .normal_power_w(),
        0.0
    );
    assert_eq!(
        original
            .field(Side::Transmitted, 0.0, window, || false)
            .unwrap()
            .points(),
        original
            .illumination_field(Face::TotalRight, 0.0, window, || false)
            .unwrap()
            .points()
    );
}
