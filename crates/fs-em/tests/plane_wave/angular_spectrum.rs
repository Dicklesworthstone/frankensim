use super::*;
use fs_em::plane_wave::{
    angular_spectrum::{
        AngularSpectrumError, AngularSpectrumGrid, AngularSpectrumLimits,
        AngularSpectrumSide as Side, PlanarAngularSpectrum,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_math::det;

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
    let material = material(1.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &material,
        DielectricResponse::Lossless(&material),
        1.0,
        &[],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 8,
        height: 4,
        spacing_m: [8.0, 9.0],
    };
    let n = grid.width * grid.height;
    // Two localized, phase-shifted openings. They define the band-limited
    // periodic field, not a claim of unrepresented hard-aperture high modes.
    let mut input = vec![[C64::ZERO; 2]; n];
    input[grid.width + 2] = [C64::ONE, C64::new(0.2, 0.4)];
    input[2 * grid.width + 5] = [C64::new(0.1, 0.8), C64::new(-0.3, 0.1)];
    let beam = PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
    for distance in [0.0, 3.0, 41.0] {
        let field = beam.field(Side::Transmitted, distance, || false).unwrap();
        let mut reference_e = vec![[C64::ZERO; 3]; n];
        let mut reference_h = reference_e.clone();
        for fy in 0..grid.height {
            for fx in 0..grid.width {
                let mx = if fx < grid.width / 2 {
                    fx as f64
                } else {
                    fx as f64 - grid.width as f64
                };
                let my = if fy < grid.height / 2 {
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
    assert!(matches!(
        PlanarAngularSpectrum::prepare(
            &stack,
            AngularSpectrumGrid {
                spacing_m: [1.0, 1.0],
                ..grid
            },
            &input,
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
