//! Authored real/nonproper source factors through the native Maxwell owners.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumExposureStressLimits, AngularSpectrumFactorStatistics as Statistics,
    AngularSpectrumIlluminationLimits, AngularSpectrumIlluminationSide as Face,
    PreparedAngularSpectrumIllumination,
};
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn stress_limits() -> AngularSpectrumExposureStressLimits {
    AngularSpectrumExposureStressLimits {
        maximum_workspace_bytes: 16 << 20,
    }
}

fn illuminate<'b, 's, 'm>(
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

fn values(point: &fs_em::plane_wave::angular_spectrum::AngularSpectrumExposurePoint) -> [f64; 9] {
    let e = point.electric_mean_squared_v2_m2;
    let h = point.magnetic_mean_squared_a2_m2;
    let s = point.poynting_w_m2;
    [e[0], e[1], e[2], h[0], h[1], h[2], s[0], s[1], s[2]]
}

fn stress_values(stress: [[f64; 3]; 3]) -> [f64; 6] {
    [
        stress[0][0],
        stress[1][1],
        stress[2][2],
        stress[0][1],
        stress[0][2],
        stress[1][2],
    ]
}

fn assert_close(a: &[f64], b: &[f64], tolerance: f64) {
    for (&x, &y) in a.iter().zip(b) {
        assert!(
            (x - y).abs() <= tolerance * y.abs().max(1.0),
            "{x:e} != {y:e}"
        );
    }
}

/// Independent positive-time real-field quadrature. One explicitly enumerated
/// coefficient spans every frequency and both illumination faces of the bank.
fn real_reference(
    fields: &[[C64; 6]],
    omega: &[f64],
    coefficients: &[C64],
    window: HarmonicWindow,
) -> [f64; 15] {
    let steps = 2048;
    let mut result = [0.0; 15];
    for &coefficient in coefficients {
        for k in 0..=steps {
            let u = k as f64 / steps as f64;
            let t = window.start_s() + window.duration_s() * u;
            let mut field = [0.0; 6];
            for (phasor, &frequency) in fields.iter().zip(omega) {
                let (sin, cos) = (frequency * t).sin_cos();
                for c in 0..6 {
                    let a = phasor[c] * coefficient;
                    field[c] += a.re * cos - a.im * sin;
                }
            }
            let simpson = if k == 0 || k == steps {
                1.0
            } else if k % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let shutter = match window.profile() {
                HarmonicWindowProfile::Uniform => 1.0,
                HarmonicWindowProfile::Triangular => 2.0 * (1.0 - (2.0 * u - 1.0).abs()),
            };
            let weight = simpson * shutter;
            let mut raw = [0.0; 15];
            for c in 0..6 {
                raw[c] = field[c] * field[c];
            }
            for c in 0..3 {
                let (a, b) = ((c + 1) % 3, (c + 2) % 3);
                raw[c + 6] = field[a] * field[b + 3] - field[b] * field[a + 3];
            }
            let trace = raw[..6].iter().sum::<f64>();
            for (c, (i, j)) in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
                .into_iter()
                .enumerate()
            {
                raw[c + 9] = field[i] * field[j] + field[i + 3] * field[j + 3]
                    - if i == j { 0.5 * trace } else { 0.0 };
            }
            for (sum, value) in result.iter_mut().zip(raw) {
                *sum += weight * value;
            }
        }
    }
    result.map(|sum| sum / (3.0 * steps as f64 * coefficients.len() as f64))
}

#[test]
fn g1_phase_locked_real_factor_retains_subcycle_field_and_stress() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 2.0, &[])
            .unwrap();
    let lines = [constant(&stack, [C64::new(0.0, 1.0), C64::ZERO])];
    let factors =
        [AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap()];
    let real = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &lines[0],
        None,
        &factors,
        &[Statistics::Real],
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let proper = AngularSpectrumExposureEnsemble::try_new(
        &lines[0],
        None,
        &factors,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    assert_eq!(real.factor_statistics(0), Some(Statistics::Real));
    assert_eq!(proper.factor_statistics(0), Some(Statistics::ProperComplex));
    assert_eq!(real.factor_statistics(1), None);
    let window =
        HarmonicWindow::with_profile(0.0, 1e-10, HarmonicWindowProfile::Triangular).unwrap();
    let deterministic = factors[0]
        .field(Side::Incident, 0.0, window, || false)
        .unwrap();
    let phase_locked = real.field(Side::Incident, 0.0, window, || false).unwrap();
    let circular = proper.field(Side::Incident, 0.0, window, || false).unwrap();
    assert_eq!(phase_locked.points(), deterministic.points());
    // E[t²] for this triangular clock is 7*T²/24, so E[sin²(2t)]~7*T²/6.
    let expected = 7e-20 / 6.0;
    let actual = phase_locked.points()[0].electric_mean_squared_v2_m2[0];
    assert!((actual / expected - 1.0).abs() < 3e-13);
    near(
        circular.points()[0].electric_mean_squared_v2_m2[0],
        0.5,
        2e-15,
    );
    let a = real
        .stress_field(Side::Incident, 0.0, window, stress_limits(), || false)
        .unwrap();
    let b = factors[0]
        .stress_field(Side::Incident, 0.0, window, stress_limits(), || false)
        .unwrap();
    assert_eq!(a.points(), b.points());
    assert!((a.points()[0].maxwell_stress_pa[2][2] / -expected - 1.0).abs() < 3e-13);
    let proper_stress = proper
        .stress_field(Side::Incident, 0.0, window, stress_limits(), || false)
        .unwrap();
    near(
        proper_stress.points()[0].maxwell_stress_pa[2][2],
        -0.5,
        2e-15,
    );
}

#[test]
fn g1_mixed_factors_preserve_joint_face_field_heat_storage_and_stress() {
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
        constant(&stacks[0], [C64::new(0.8, 0.2), C64::new(-0.1, 0.3)]),
        constant(&stacks[1], [C64::new(-0.3, 0.6), C64::new(0.4, 0.1)]),
    ];
    let lines = [
        illuminate(&left[0], [C64::new(-0.2, 0.4), C64::new(0.3, -0.2)]),
        illuminate(&left[1], [C64::new(0.7, -0.1), C64::new(-0.1, 0.2)]),
    ];
    let mean = AngularSpectrumExposure::try_new_illumination(&left[0], &lines, limits(), || false)
        .unwrap();
    let factors = [
        AngularSpectrumExposure::try_new_illumination(&left[0], &lines, limits(), || false)
            .unwrap(),
        AngularSpectrumExposure::try_new_illumination(&left[0], &lines, limits(), || false)
            .unwrap(),
    ];
    let statistics = [Statistics::Real, Statistics::ProperComplex];
    let ensemble = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &left[0],
        Some(&mean),
        &factors,
        &statistics,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let layer = ensemble.prepare_layer(0, heat_limits(), || false).unwrap();
    // Original two-port illumination already owns both source columns: no new
    // TE/TM film solve is introduced by the real/proper statistical choice.
    assert_eq!(layer.preparation_layer_mode_pairs(), 0);
    assert_eq!(layer.frequency_count(), 2);
    let phases = [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)];
    let coefficients: Vec<_> = [-1.0, 1.0]
        .into_iter()
        .flat_map(|real| {
            phases
                .into_iter()
                .map(move |z| C64::from_re(1.0 + real) + z)
        })
        .collect();
    let external = lines
        .iter()
        .map(|line| {
            let f = line.field(Face::TotalLeft, 0.13, || false).unwrap();
            let e = f.electric_v_per_m()[0];
            let h = f.magnetic_a_per_m()[0];
            [e[0], e[1], e[2], h[0], h[1], h[2]]
        })
        .collect::<Vec<_>>();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.12, profile).unwrap();
        let expected = real_reference(&external, &frequencies, &coefficients, window);
        let observed = ensemble
            .illumination_field(Face::TotalLeft, 0.13, window, || false)
            .unwrap();
        let stress = ensemble
            .illumination_stress_field(Face::TotalLeft, 0.13, window, stress_limits(), || false)
            .unwrap();
        for point in observed.points() {
            assert_close(&values(point), &expected[..9], 2e-11);
        }
        for point in stress.points() {
            assert_close(
                &stress_values(point.maxwell_stress_pa),
                &expected[9..],
                2e-11,
            );
        }
    }
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    let native = lines
        .iter()
        .map(|line| line.layer_field(0, 0.155, || false).unwrap())
        .collect::<Vec<_>>();
    let expected = heat_time_reference(&native, &frequencies, 0.17, 0.12, &coefficients);
    let observed = layer.field(0.155, window, || false).unwrap();
    for q in observed.dissipation() {
        assert_close(
            &[
                q.ohmic_heat_w_m3,
                q.polarization_heat_w_m3,
                q.storage_rate_w_m3,
            ],
            &expected,
            3e-11,
        );
    }
    let real_only = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &left[0],
        None,
        &factors[..1],
        &statistics[..1],
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let real_heat = real_only
        .prepare_layer(0, heat_limits(), || false)
        .unwrap()
        .field(0.155, window, || false)
        .unwrap();
    let deterministic_heat = mean
        .prepare_layer(0, heat_limits(), || false)
        .unwrap()
        .field(0.155, window, || false)
        .unwrap();
    assert_eq!(
        real_heat.field().points(),
        deterministic_heat.field().points()
    );
    assert_eq!(real_heat.dissipation(), deterministic_heat.dissipation());
    let real_resultant = real_only
        .boundary_stress_resultant_n([0.0, 0.0], window, stress_limits(), || false)
        .unwrap();
    let deterministic_resultant = mean
        .boundary_stress_resultant_n([0.0, 0.0], window, stress_limits(), || false)
        .unwrap();
    assert_eq!(real_resultant, deterministic_resultant);
}

#[test]
fn g3_proper_factor_equals_two_real_quadratures_across_all_observers() {
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let response = DielectricResponse::DebyeWithOhmic {
        dielectric: &dielectric,
        conductor: &conductor,
    };
    let stacks = [2.0, 2.3].map(|omega| {
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
    let amplitudes = [
        [C64::new(0.8, 0.2), C64::new(-0.1, 0.3)],
        [C64::new(-0.3, 0.6), C64::new(0.4, 0.1)],
    ];
    let lines = [
        constant(&stacks[0], amplitudes[0]),
        constant(&stacks[1], amplitudes[1]),
    ];
    let rotated: [Vec<PlanarAngularSpectrum<'_, '_>>; 2] = core::array::from_fn(|r| {
        stacks
            .iter()
            .zip(amplitudes)
            .map(|(stack, amplitude)| {
                let z = if r == 0 {
                    C64::from_re(core::f64::consts::FRAC_1_SQRT_2)
                } else {
                    C64::new(0.0, core::f64::consts::FRAC_1_SQRT_2)
                };
                constant(stack, amplitude.map(|a| a * z))
            })
            .collect()
    });
    let proper_factors =
        [AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap()];
    let real_factors = [
        AngularSpectrumExposure::try_new(&lines[0], &rotated[0], limits(), || false).unwrap(),
        AngularSpectrumExposure::try_new(&lines[0], &rotated[1], limits(), || false).unwrap(),
    ];
    let proper = AngularSpectrumExposureEnsemble::try_new(
        &lines[0],
        None,
        &proper_factors,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let real = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &lines[0],
        None,
        &real_factors,
        &[Statistics::Real; 2],
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let window =
        HarmonicWindow::with_profile(0.17, 0.12, HarmonicWindowProfile::Triangular).unwrap();
    let a = real
        .field(Side::Transmitted, 0.1, window, || false)
        .unwrap();
    let b = proper
        .field(Side::Transmitted, 0.1, window, || false)
        .unwrap();
    for (a, b) in a.points().iter().zip(b.points()) {
        assert_close(&values(a), &values(b), 3e-12);
    }
    let a = real
        .stress_field(
            Side::TotalIncidentSide,
            0.1,
            window,
            stress_limits(),
            || false,
        )
        .unwrap();
    let b = proper
        .stress_field(
            Side::TotalIncidentSide,
            0.1,
            window,
            stress_limits(),
            || false,
        )
        .unwrap();
    for (a, b) in a.points().iter().zip(b.points()) {
        assert_close(
            &stress_values(a.maxwell_stress_pa),
            &stress_values(b.maxwell_stress_pa),
            3e-12,
        );
    }
    let a = real
        .prepare_layer(0, heat_limits(), || false)
        .unwrap()
        .field(0.155, window, || false)
        .unwrap();
    let b = proper
        .prepare_layer(0, heat_limits(), || false)
        .unwrap()
        .field(0.155, window, || false)
        .unwrap();
    for (a, b) in a.dissipation().iter().zip(b.dissipation()) {
        assert_close(
            &[a.total_heat_w_m3, a.storage_rate_w_m3],
            &[b.total_heat_w_m3, b.storage_rate_w_m3],
            3e-12,
        );
    }
}

#[test]
fn g4_statistics_shape_identity_budgets_darkness_and_cancellation_are_explicit() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 2.0, &[])
            .unwrap();
    let lines = [constant(&stack, [C64::ONE, C64::ZERO])];
    let factors =
        [AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap()];
    for statistics in [&[][..], &[Statistics::Real; 2][..]] {
        assert!(
            AngularSpectrumExposureEnsemble::try_new_with_statistics(
                &lines[0],
                None,
                &factors,
                statistics,
                ensemble_limits(),
                || false,
            )
            .is_err()
        );
    }
    let denied = AngularSpectrumExposureEnsembleLimits {
        maximum_factors: 0,
        ..ensemble_limits()
    };
    assert!(
        AngularSpectrumExposureEnsemble::try_new_with_statistics(
            &lines[0],
            None,
            &factors,
            &[Statistics::Real],
            denied,
            || false,
        )
        .is_err()
    );
    assert!(matches!(
        AngularSpectrumExposureEnsemble::try_new_with_statistics(
            &lines[0],
            None,
            &factors,
            &[Statistics::Real],
            ensemble_limits(),
            || true,
        ),
        Err(AngularSpectrumError::Cancelled)
    ));
    let other_medium = material(1.1, 1.0);
    let other_stack = PlanarDielectricLayers::try_new(
        &other_medium,
        DielectricResponse::Lossless(&other_medium),
        2.0,
        &[],
    )
    .unwrap();
    let other_lines = [constant(&other_stack, [C64::ONE, C64::ZERO])];
    let other = [
        AngularSpectrumExposure::try_new(&other_lines[0], &other_lines, limits(), || false)
            .unwrap(),
    ];
    assert!(
        AngularSpectrumExposureEnsemble::try_new_with_statistics(
            &lines[0],
            None,
            &other,
            &[Statistics::Real],
            ensemble_limits(),
            || false,
        )
        .is_err()
    );
    let real = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &lines[0],
        None,
        &factors,
        &[Statistics::Real],
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let window = HarmonicWindow::new(0.1, 0.2).unwrap();
    let mut calls = 0;
    assert!(matches!(
        real.field(Side::Incident, 0.0, window, || {
            calls += 1;
            calls > 5
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    let a = real.field(Side::Incident, 0.0, window, || false).unwrap();
    let b = real.field(Side::Incident, 0.0, window, || false).unwrap();
    assert_eq!(a.points(), b.points());
    let dark = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &lines[0],
        None,
        &[],
        &[],
        ensemble_limits(),
        || false,
    )
    .unwrap();
    assert_eq!(dark.factor_count(), 0);
    assert_eq!(dark.factor_statistics(0), None);
    let output = dark.field(Side::Incident, 0.0, window, || false).unwrap();
    assert!(
        output
            .points()
            .iter()
            .all(|point| values(point).iter().all(|&v| v == 0.0))
    );
}
