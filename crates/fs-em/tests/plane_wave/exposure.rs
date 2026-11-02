//! G1/G3 finite-window vector interference through the actual film operator.
mod two_sided;
mod stress;
use super::*;
use fs_em::plane_wave::{
    angular_spectrum::{
        AngularSpectrumCoherence, AngularSpectrumCoherenceLimits, AngularSpectrumError,
        AngularSpectrumExposure, AngularSpectrumExposureEnsemble,
        AngularSpectrumExposureEnsembleLimits, AngularSpectrumExposureLimits, AngularSpectrumGrid,
        AngularSpectrumExposureLayerLimits, AngularSpectrumLayerField, AngularSpectrumLayerLimits,
        AngularSpectrumLimits, AngularSpectrumSide as Side, PlanarAngularSpectrum,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductorState};

fn grid() -> AngularSpectrumGrid {
    AngularSpectrumGrid {
        width: 4,
        height: 1,
        spacing_m: [8.0, 3.0],
    }
}
fn limits() -> AngularSpectrumExposureLimits {
    AngularSpectrumExposureLimits {
        maximum_lines: 8,
        maximum_point_line_pairs: 1024,
        maximum_point_frequency_pairs: 1024,
        maximum_workspace_bytes: 2 << 20,
        maximum_fft_work: 1 << 20,
    }
}
fn ensemble_limits() -> AngularSpectrumExposureEnsembleLimits {
    AngularSpectrumExposureEnsembleLimits {
        maximum_factors: 8,
        maximum_point_line_pairs: 1024,
        maximum_point_frequency_pairs: 1024,
        maximum_workspace_bytes: 4 << 20,
        maximum_fft_work: 1 << 20,
    }
}
fn prepare<'s, 'm>(
    source: &'s PlanarDielectricLayers<'m>,
    samples: &[[C64; 2]],
) -> PlanarAngularSpectrum<'s, 'm> {
    PlanarAngularSpectrum::prepare(
        source,
        grid(),
        samples,
        AngularSpectrumLimits {
            maximum_grid_points: 128,
            maximum_workspace_bytes: 1 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap()
}
fn constant<'s, 'm>(
    source: &'s PlanarDielectricLayers<'m>,
    e: [C64; 2],
) -> PlanarAngularSpectrum<'s, 'm> {
    prepare(source, &[e; 4])
}

#[test]
fn g1_proper_spectral_covariance_and_mean_match_authored_ensemble_through_a_film() {
    let medium = material(E, M);
    let film = finite_layers::debye(2.0 * E, E, 1.0, 1.0);
    let films = [DielectricLayer {
        response: DielectricResponse::Debye(&film),
        thickness_m: 0.4,
    }];
    let stacks = [1.0e9, 1.17e9].map(|omega| {
        PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&medium),
            omega,
            &films,
        )
        .unwrap()
    });
    let samples: Vec<_> = (0..4)
        .map(|x| {
            let phase = -(x as f64) * core::f64::consts::FRAC_PI_2;
            let p = C64::new(phase.cos(), phase.sin());
            [C64::new(0.3, -0.2) * p, C64::new(0.7, 0.1) * p]
        })
        .collect();
    let mean_lines = [constant(&stacks[0], [C64::new(0.4, -0.1), C64::ZERO])];
    let first_lines = [
        constant(&stacks[0], [C64::new(0.8, 0.2), C64::new(0.1, -0.6)]),
        prepare(&stacks[1], &samples),
    ];
    let second_lines = [constant(
        &stacks[1],
        [C64::new(-0.2, 0.3), C64::new(0.4, 0.1)],
    )];
    let mean =
        AngularSpectrumExposure::try_new(&mean_lines[0], &mean_lines, limits(), || false).unwrap();
    let factors = [
        AngularSpectrumExposure::try_new(&first_lines[0], &first_lines, limits(), || false)
            .unwrap(),
        AngularSpectrumExposure::try_new(&second_lines[0], &second_lines, limits(), || false)
            .unwrap(),
    ];
    let bank = AngularSpectrumExposureEnsemble::try_new(
        &mean_lines[0],
        Some(&mean),
        &factors,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    assert_eq!(bank.factor_count(), 2);
    assert_eq!(bank.line_count(), 4);
    let (start, duration, distance) = (2e-9, 4e-9, 0.23);
    let result = bank
        .field(
            Side::Transmitted,
            distance,
            HarmonicWindow::new(start, duration).unwrap(),
            || false,
        )
        .unwrap();
    let lines = [
        &mean_lines[0],
        &first_lines[0],
        &first_lines[1],
        &second_lines[0],
    ];
    let fields = lines
        .iter()
        .map(|line| line.field(Side::Transmitted, distance, || false).unwrap())
        .collect::<Vec<_>>();
    let phases = [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)];
    let steps = 2048;
    for (point, observed) in result.points().iter().enumerate() {
        let mut expected = [0.0; 9];
        // Enumerate two independent, explicitly authored proper four-state
        // amplitudes and integrate their REAL fields plus the coherent mean.
        // No production temporal kernel or covariance contraction is reused.
        for z in phases {
            for w in phases {
                for step in 0..=steps {
                    let time = start + duration * step as f64 / steps as f64;
                    let mut real = [0.0; 6];
                    for ((line, field), coefficient) in
                        lines.iter().zip(&fields).zip([C64::ONE, z, z, w])
                    {
                        let (sin, cos) = (line.source().angular_frequency_rad_s() * time).sin_cos();
                        for (component, amplitude) in field.electric_v_per_m()[point]
                            .iter()
                            .chain(&field.magnetic_a_per_m()[point])
                            .enumerate()
                        {
                            let amplitude = *amplitude * coefficient;
                            real[component] += amplitude.re * cos - amplitude.im * sin;
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
                        expected[c] += weight * real[c] * real[c];
                        expected[c + 3] += weight * real[c + 3] * real[c + 3];
                        let (a, b) = ((c + 1) % 3, (c + 2) % 3);
                        expected[c + 6] += weight * (real[a] * real[b + 3] - real[b] * real[a + 3]);
                    }
                }
            }
        }
        for (actual, expected) in observed
            .electric_mean_squared_v2_m2
            .iter()
            .chain(&observed.magnetic_mean_squared_a2_m2)
            .chain(&observed.poynting_w_m2)
            .zip(expected.map(|x| x / (48.0 * steps as f64)))
        {
            near(*actual, expected, 3e-10 * expected.abs().max(1e-14));
        }
    }
}

#[test]
fn g3_proper_monochromatic_exposure_matches_existing_spatial_coherence() {
    let medium = material(E, M);
    let film = finite_layers::debye(2.0 * E, E, 1.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&medium),
        1e9,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&film),
            thickness_m: 0.3,
        }],
    )
    .unwrap();
    let lines = [
        constant(&stack, [C64::new(0.7, -0.2), C64::new(0.1, 0.5)]),
        constant(&stack, [C64::new(-0.3, 0.4), C64::new(0.6, -0.1)]),
    ];
    let spatial = AngularSpectrumCoherence::try_new(
        &lines[0],
        &lines,
        AngularSpectrumCoherenceLimits {
            maximum_modes: 8,
            maximum_point_mode_pairs: 1024,
            maximum_workspace_bytes: 2 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap();
    let factors = lines
        .iter()
        .map(|line| {
            AngularSpectrumExposure::try_new(line, std::slice::from_ref(line), limits(), || false)
                .unwrap()
        })
        .collect::<Vec<_>>();
    let ensemble = AngularSpectrumExposureEnsemble::try_new(
        &lines[0],
        None,
        &factors,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    for side in [Side::Incident, Side::Reflected, Side::Transmitted] {
        let reference = spatial.field(side, 0.23, || false).unwrap();
        for window in [
            HarmonicWindow::new(0.0, 1e-19).unwrap(),
            HarmonicWindow::new(2.0e-9, 1.0e-9).unwrap(),
        ] {
            let observed = ensemble.field(side, 0.23, window, || false).unwrap();
            for (actual, expected) in observed.points().iter().zip(reference.points()) {
                for c in 0..3 {
                    for (a, b) in [
                        (
                            actual.electric_mean_squared_v2_m2[c],
                            expected.electric_coherency_v2_m2[c][c].re,
                        ),
                        (
                            actual.magnetic_mean_squared_a2_m2[c],
                            expected.magnetic_rms_squared_a2_m2[c],
                        ),
                        (
                            actual.poynting_w_m2[c],
                            expected.complex_poynting_w_m2[c].re,
                        ),
                    ] {
                        near(a, b, 2e-13 * b.abs().max(1e-18));
                    }
                }
            }
            near(observed.normal_power_w(), reference.normal_power_w(), 1e-13);
        }
    }
}

#[test]
fn g1_dispersive_vector_interference_matches_independent_real_time_quadrature() {
    let medium = material(E, M);
    let film = finite_layers::debye(2.0 * E, E, 1.0, 1.0);
    let films = [DielectricLayer {
        response: DielectricResponse::Debye(&film),
        thickness_m: 0.4,
    }];
    let stacks = [1.0e9, 1.17e9].map(|omega| {
        PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&medium),
            omega,
            &films,
        )
        .unwrap()
    });
    let samples: Vec<_> = (0..4)
        .map(|x| {
            let phase = -(x as f64) * core::f64::consts::FRAC_PI_2;
            let p = C64::new(phase.cos(), phase.sin());
            [C64::new(0.3, -0.2) * p, C64::new(0.7, 0.1) * p]
        })
        .collect();
    let lines = [
        constant(&stacks[0], [C64::new(0.8, 0.2), C64::new(0.1, -0.6)]),
        prepare(&stacks[1], &samples),
    ];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let (start, duration, distance) = (2e-9, 4e-9, 0.23);
    let result = bank
        .field(
            Side::Transmitted,
            distance,
            HarmonicWindow::new(start, duration).unwrap(),
            || false,
        )
        .unwrap();
    let fields = lines
        .iter()
        .map(|line| line.field(Side::Transmitted, distance, || false).unwrap())
        .collect::<Vec<_>>();
    // Simpson integration of actual REAL E/H, independent of any harmonic
    // window/sinc calculation and without using a scalar impedance shortcut.
    let steps = 4096;
    for (point, observed) in result.points().iter().enumerate() {
        let mut expected = [0.0; 9];
        for step in 0..=steps {
            let time = start + duration * step as f64 / steps as f64;
            let mut real = [0.0; 6];
            for (line, field) in lines.iter().zip(&fields) {
                let phase = line.source().angular_frequency_rad_s() * time;
                let (sin, cos) = phase.sin_cos();
                for (component, z) in field.electric_v_per_m()[point]
                    .iter()
                    .chain(&field.magnetic_a_per_m()[point])
                    .enumerate()
                {
                    real[component] += z.re * cos - z.im * sin;
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
                expected[c] += weight * real[c] * real[c];
                expected[c + 3] += weight * real[c + 3] * real[c + 3];
                let (a, b) = ((c + 1) % 3, (c + 2) % 3);
                expected[c + 6] += weight * (real[a] * real[b + 3] - real[b] * real[a + 3]);
            }
        }
        for (actual, expected) in observed
            .electric_mean_squared_v2_m2
            .iter()
            .chain(&observed.magnetic_mean_squared_a2_m2)
            .chain(&observed.poynting_w_m2)
            .zip(expected.map(|x| x / (3.0 * steps as f64)))
        {
            near(*actual, expected, 2e-10 * expected.abs().max(1e-14));
        }
    }
    let expected_power = result
        .points()
        .iter()
        .map(|p| p.poynting_w_m2[2])
        .sum::<f64>()
        * grid().spacing_m[0]
        * grid().spacing_m[1];
    near(result.normal_power_w(), expected_power, 1e-14);
}

#[test]
fn g1_partial_spectral_coherence_controls_beats_without_inventing_carrier_phase() {
    let medium = material(1.0, 1.0);
    let stacks = [10.0, 11.0].map(|omega| {
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap()
    });
    let mean_lines = [constant(&stacks[0], [C64::ONE, C64::ZERO])];
    let mean =
        AngularSpectrumExposure::try_new(&mean_lines[0], &mean_lines, limits(), || false).unwrap();
    for rho in [0.0, 0.6, 1.0] {
        let first = [
            constant(&stacks[0], [C64::ONE, C64::ZERO]),
            constant(&stacks[1], [C64::new(0.6 * rho, 0.8 * rho), C64::ZERO]),
        ];
        let second = [constant(
            &stacks[1],
            [C64::from_re(fs_math::det::sqrt(1.0 - rho * rho)), C64::ZERO],
        )];
        let factors = [
            AngularSpectrumExposure::try_new(&first[0], &first, limits(), || false).unwrap(),
            AngularSpectrumExposure::try_new(&second[0], &second, limits(), || false).unwrap(),
        ];
        let ensemble = AngularSpectrumExposureEnsemble::try_new(
            &mean_lines[0],
            None,
            &factors,
            ensemble_limits(),
            || false,
        )
        .unwrap();
        let illuminated = AngularSpectrumExposureEnsemble::try_new(
            &mean_lines[0],
            Some(&mean),
            &factors,
            ensemble_limits(),
            || false,
        )
        .unwrap();
        let duration = 1e-10;
        for start in [0.0, core::f64::consts::FRAC_PI_2] {
            let window = HarmonicWindow::new(start, duration).unwrap();
            let short = ensemble
                .field(Side::Transmitted, 0.0, window, || false)
                .unwrap();
            // Beat angular frequency is 1 rad/s. The second factor keeps the
            // diagonal source powers fixed while rho controls their coherence.
            let midpoint = start + duration * 0.5;
            let sinc = (0.5 * duration).sin() / (0.5 * duration);
            let expected = 1.0 + rho * (0.6 * midpoint.cos() - 0.8 * midpoint.sin()) * sinc;
            near(short.points()[0].poynting_w_m2[2], expected, 2e-14);
            let total = illuminated
                .field(Side::Transmitted, 0.0, window, || false)
                .unwrap();
            let coherent = mean
                .field(Side::Transmitted, 0.0, window, || false)
                .unwrap();
            near(
                total.points()[0].poynting_w_m2[2],
                expected + coherent.points()[0].poynting_w_m2[2],
                3e-14,
            );
        }
        let long = ensemble
            .field(
                Side::Transmitted,
                0.0,
                HarmonicWindow::new(0.0, core::f64::consts::TAU).unwrap(),
                || false,
            )
            .unwrap();
        near(long.points()[0].poynting_w_m2[2], 1.0, 2e-14);
    }
}

#[test]
fn g4_ensemble_budgets_identity_cancellation_and_empty_factors_are_explicit() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 2.0, &[])
            .unwrap();
    let lines = [constant(&stack, [C64::ONE, C64::ZERO])];
    let factors =
        [AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap()];
    let admitted = AngularSpectrumExposureEnsemble::try_new(
        &lines[0],
        Some(&factors[0]),
        &factors,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    assert_eq!(admitted.point_line_pairs(), 8);
    assert_eq!(admitted.point_frequency_pairs(), 8);
    assert_eq!(admitted.field_fft_work(), 2 * factors[0].field_fft_work());
    for cap in [
        AngularSpectrumExposureEnsembleLimits {
            maximum_factors: 0,
            ..ensemble_limits()
        },
        AngularSpectrumExposureEnsembleLimits {
            maximum_point_line_pairs: admitted.point_line_pairs() - 1,
            ..ensemble_limits()
        },
        AngularSpectrumExposureEnsembleLimits {
            maximum_point_frequency_pairs: admitted.point_frequency_pairs() - 1,
            ..ensemble_limits()
        },
        AngularSpectrumExposureEnsembleLimits {
            maximum_workspace_bytes: admitted.workspace_bytes() - 1,
            ..ensemble_limits()
        },
        AngularSpectrumExposureEnsembleLimits {
            maximum_fft_work: admitted.field_fft_work() - 1,
            ..ensemble_limits()
        },
    ] {
        assert!(
            AngularSpectrumExposureEnsemble::try_new(
                &lines[0],
                Some(&factors[0]),
                &factors,
                cap,
                || false,
            )
            .is_err()
        );
    }
    let exact = AngularSpectrumExposureEnsembleLimits {
        maximum_factors: 1,
        maximum_point_line_pairs: admitted.point_line_pairs(),
        maximum_point_frequency_pairs: admitted.point_frequency_pairs(),
        maximum_workspace_bytes: admitted.workspace_bytes(),
        maximum_fft_work: admitted.field_fft_work(),
    };
    assert!(
        AngularSpectrumExposureEnsemble::try_new(
            &lines[0],
            Some(&factors[0]),
            &factors,
            exact,
            || false,
        )
        .is_ok()
    );
    let replacement = material(1.0, 1.0);
    let other = PlanarDielectricLayers::try_new(
        &replacement,
        DielectricResponse::Lossless(&replacement),
        2.0,
        &[],
    )
    .unwrap();
    let changed = [constant(&other, [C64::ONE, C64::ZERO])];
    let changed_banks =
        [AngularSpectrumExposure::try_new(&changed[0], &changed, limits(), || false).unwrap()];
    assert!(
        AngularSpectrumExposureEnsemble::try_new(
            &lines[0],
            None,
            &changed_banks,
            ensemble_limits(),
            || false,
        )
        .is_err()
    );
    assert!(
        AngularSpectrumExposureEnsemble::try_new(
            &lines[0],
            Some(&changed_banks[0]),
            &factors,
            ensemble_limits(),
            || false,
        )
        .is_err()
    );
    let window = HarmonicWindow::new(0.13, 0.71).unwrap();
    let mut polls = 0;
    assert!(matches!(
        admitted.field(Side::Transmitted, 0.0, window, || {
            polls += 1;
            polls > 5
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    let retry = admitted
        .field(Side::Transmitted, 0.0, window, || false)
        .unwrap();
    let replay = admitted
        .field(Side::Transmitted, 0.0, window, || false)
        .unwrap();
    assert_eq!(retry.points(), replay.points());
    let mean_only = AngularSpectrumExposureEnsemble::try_new(
        &lines[0],
        Some(&factors[0]),
        &[],
        ensemble_limits(),
        || false,
    )
    .unwrap()
    .field(Side::Transmitted, 0.0, window, || false)
    .unwrap();
    let original = factors[0]
        .field(Side::Transmitted, 0.0, window, || false)
        .unwrap();
    assert_eq!(mean_only.points(), original.points());
    let dark =
        AngularSpectrumExposureEnsemble::try_new(&lines[0], None, &[], ensemble_limits(), || false)
            .unwrap()
            .field(Side::Transmitted, 0.0, window, || false)
            .unwrap();
    assert_eq!(dark.normal_power_w(), 0.0);
    assert!(
        dark.points()
            .iter()
            .all(|p| p.electric_mean_squared_v2_m2 == [0.0; 3])
    );
}

#[test]
fn g3_equal_frequency_sources_cancel_before_powers_and_preserve_small_residuals() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 2.0, &[])
            .unwrap();
    let lines =
        [1e16, 1.0, -1e16].map(|amplitude| constant(&stack, [C64::from_re(amplitude), C64::ZERO]));
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    assert_eq!(bank.line_count(), 3);
    assert_eq!(bank.frequency_count(), 1);
    assert_eq!(bank.point_frequency_pairs(), 4);
    let window = HarmonicWindow::new(0.27, core::f64::consts::PI).unwrap();
    let result = bank
        .field(Side::Transmitted, 0.0, window, || false)
        .unwrap();
    for point in result.points() {
        near(point.electric_mean_squared_v2_m2[0], 0.5, 2e-14);
        near(point.poynting_w_m2[2], 0.5, 2e-14);
    }
    let null = [1.0, -1.0].map(|a| constant(&stack, [C64::from_re(a), C64::ZERO]));
    let dark = AngularSpectrumExposure::try_new(&null[0], &null, limits(), || false)
        .unwrap()
        .field(Side::Transmitted, 0.0, window, || false)
        .unwrap();
    assert_eq!(dark.normal_power_w(), 0.0);
    assert!(
        dark.points()
            .iter()
            .all(|p| p.electric_mean_squared_v2_m2 == [0.0; 3])
    );
    // Noncontiguous equal-frequency sources remain in one coherent group;
    // the independently ordered frequency is not folded into that field.
    let other =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 3.0, &[])
            .unwrap();
    let interleaved = [
        constant(&other, [C64::ONE, C64::ZERO]),
        constant(&stack, [C64::from_re(1e16), C64::ZERO]),
        constant(&other, [C64::ONE, C64::ZERO]),
        constant(&stack, [C64::from_re(-1e16), C64::ZERO]),
    ];
    let grouped =
        AngularSpectrumExposure::try_new(&interleaved[0], &interleaved, limits(), || false)
            .unwrap();
    assert_eq!(grouped.frequency_count(), 2);
    assert_eq!(
        grouped.field_fft_work(),
        2 * interleaved[0].field_fft_work()
    );
    let exposure = grouped
        .field(Side::Transmitted, 0.0, window, || false)
        .unwrap();
    for point in exposure.points() {
        near(point.electric_mean_squared_v2_m2[0], 2.0, 5e-14);
        near(point.poynting_w_m2[2], 2.0, 5e-14);
    }
}

#[test]
fn g1_short_shutters_retain_carrier_phase_and_long_shutters_wash_out_beats() {
    let medium = material(1.0, 1.0);
    let stacks = [10.0, 11.0].map(|omega| {
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap()
    });
    let lines = stacks
        .iter()
        .map(|stack| constant(stack, [C64::ONE, C64::ZERO]))
        .collect::<Vec<_>>();
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let short = bank
        .field(
            Side::Transmitted,
            0.0,
            HarmonicWindow::new(0.0, 1e-8).unwrap(),
            || false,
        )
        .unwrap();
    near(short.points()[0].poynting_w_m2[2], 4.0, 1e-13);
    let long = bank
        .field(
            Side::Transmitted,
            0.0,
            HarmonicWindow::new(0.0, core::f64::consts::TAU).unwrap(),
            || false,
        )
        .unwrap();
    near(long.points()[0].poynting_w_m2[2], 1.0, 1e-13);
    let sine = [constant(&stacks[0], [C64::new(0.0, 1.0), C64::ZERO])];
    let tiny = AngularSpectrumExposure::try_new(&sine[0], &sine, limits(), || false)
        .unwrap()
        .field(
            Side::Transmitted,
            0.0,
            HarmonicWindow::new(0.0, 1e-10).unwrap(),
            || false,
        )
        .unwrap();
    near(
        tiny.points()[0].electric_mean_squared_v2_m2[0] / (1e-18 / 3.0),
        1.0,
        1e-12,
    );
}

#[test]
fn g1_vector_polarization_and_signed_reflection_remain_physical() {
    let medium = material(1.0, 1.0);
    let dielectric = material(4.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&dielectric),
        2.0,
        &[],
    )
    .unwrap();
    let lines = [
        constant(&stack, [C64::ONE, C64::ZERO]),
        constant(&stack, [C64::ZERO, C64::ONE]),
    ];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let window = HarmonicWindow::new(0.2, core::f64::consts::PI).unwrap();
    let incident = bank.field(Side::Incident, 0.0, window, || false).unwrap();
    let reflected = bank.field(Side::Reflected, 0.0, window, || false).unwrap();
    near(incident.points()[0].poynting_w_m2[2], 1.0, 1e-14);
    near(reflected.points()[0].poynting_w_m2[2], -1.0 / 9.0, 1e-14);
}

#[test]
fn g4_source_budget_and_cancellation_refusals_do_not_publish_partial_exposures() {
    let medium = material(1.0, 1.0);
    let replacement = material(1.0, 1.0);
    let a =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 2.0, &[])
            .unwrap();
    let b = PlanarDielectricLayers::try_new(
        &replacement,
        DielectricResponse::Lossless(&replacement),
        2.0,
        &[],
    )
    .unwrap();
    let lines = [constant(&a, [C64::ONE, C64::ZERO])];
    let changed = [constant(&b, [C64::ONE, C64::ZERO])];
    assert!(AngularSpectrumExposure::try_new(&lines[0], &changed, limits(), || false).is_err());
    let admitted = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    for cap in [
        AngularSpectrumExposureLimits {
            maximum_lines: 0,
            ..limits()
        },
        AngularSpectrumExposureLimits {
            maximum_point_line_pairs: 3,
            ..limits()
        },
        AngularSpectrumExposureLimits {
            maximum_point_frequency_pairs: 3,
            ..limits()
        },
        AngularSpectrumExposureLimits {
            maximum_workspace_bytes: admitted.workspace_bytes() - 1,
            ..limits()
        },
        AngularSpectrumExposureLimits {
            maximum_fft_work: admitted.field_fft_work() - 1,
            ..limits()
        },
    ] {
        assert!(AngularSpectrumExposure::try_new(&lines[0], &lines, cap, || false).is_err());
    }
    let exact = AngularSpectrumExposureLimits {
        maximum_lines: 1,
        maximum_point_line_pairs: 4,
        maximum_point_frequency_pairs: admitted.point_frequency_pairs(),
        maximum_workspace_bytes: admitted.workspace_bytes(),
        maximum_fft_work: admitted.field_fft_work(),
    };
    assert!(AngularSpectrumExposure::try_new(&lines[0], &lines, exact, || false).is_ok());
    let window = HarmonicWindow::new(0.0, 1.0).unwrap();
    let mut polls = 0;
    assert!(matches!(
        admitted.field(Side::Transmitted, 0.0, window, || {
            polls += 1;
            polls > 4
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    let retry = admitted
        .field(Side::Transmitted, 0.0, window, || false)
        .unwrap();
    let replay = admitted
        .field(Side::Transmitted, 0.0, window, || false)
        .unwrap();
    assert_eq!(retry.points(), replay.points());
    let dark = AngularSpectrumExposure::try_new(&lines[0], &[], limits(), || false)
        .unwrap()
        .field(Side::Transmitted, 0.0, window, || false)
        .unwrap();
    assert_eq!(dark.normal_power_w(), 0.0);
}

#[test]
fn g3_lossy_spectral_lines_retain_original_ohmic_and_debye_owners() {
    let incident = material(E, M);
    let basis = ElectricalPropertyBasis::Conductivity;
    // Use the same complete source bundle for Debye polarization and Ohmic
    // conduction, as required by the original combined constitutive owner.
    let shared = source(
        &[
            (EPS_NAME, EPS, 2.0 * E),
            (MU_NAME, MU, M),
            ("delta", EPS, E),
            ("tau", Time::DIMS, 2e-9),
            (basis.property_name(), basis.dimensions(), 0.02),
        ],
        ValidityDomain::unconstrained(),
        QueryPoint::new(),
    );
    let conductor = ResolvedConductorState::try_from_material(shared.clone()).unwrap();
    let replacement_conductor = ResolvedConductorState::try_from_material(shared.clone()).unwrap();
    let make_dielectric = || {
        ResolvedDebyeCapacitor::try_from_material(
            shared.clone(),
            EPS_NAME,
            &[DebyePoleProperties {
                permittivity_increment: "delta".into(),
                relaxation_time: "tau".into(),
            }],
            0.03,
            0.002,
        )
        .unwrap()
    };
    let dielectric = make_dielectric();
    let replacement_dielectric = make_dielectric();
    let ohmic = DielectricResponse::Ohmic(&conductor);
    let combined = DielectricResponse::DebyeWithOhmic {
        dielectric: &dielectric,
        conductor: &conductor,
    };
    for (response, replacements) in [
        (
            ohmic,
            [
                DielectricResponse::Ohmic(&replacement_conductor),
                combined,
                DielectricResponse::Lossless(conductor.material()),
            ],
        ),
        (
            combined,
            [
                DielectricResponse::DebyeWithOhmic {
                    dielectric: &dielectric,
                    conductor: &replacement_conductor,
                },
                DielectricResponse::DebyeWithOhmic {
                    dielectric: &replacement_dielectric,
                    conductor: &conductor,
                },
                ohmic,
            ],
        ),
    ] {
        let stacks = [1e9, 1.2e9].map(|omega| {
            PlanarDielectricLayers::try_new(
                &incident,
                response,
                omega,
                &[DielectricLayer {
                    response,
                    thickness_m: 0.07,
                }],
            )
            .unwrap()
        });
        let lines = stacks
            .iter()
            .map(|stack| constant(stack, [C64::ONE, C64::new(0.0, 0.5)]))
            .collect::<Vec<_>>();
        let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
        let observed = bank
            .field(
                Side::Transmitted,
                0.03,
                HarmonicWindow::new(0.2e-9, 4e-9).unwrap(),
                || false,
            )
            .unwrap();
        assert!(observed.normal_power_w().is_finite());
        // Identical source values do not authorize a substituted model owner,
        // and the same owners do not authorize changing the constitutive law.
        for replacement in replacements {
            let changed = PlanarDielectricLayers::try_new(
                &incident,
                replacement,
                1.2e9,
                &[DielectricLayer {
                    response: replacement,
                    thickness_m: 0.07,
                }],
            )
            .unwrap();
            let changed_lines = [
                constant(&stacks[0], [C64::ONE, C64::ZERO]),
                constant(&changed, [C64::ONE, C64::ZERO]),
            ];
            assert!(
                AngularSpectrumExposure::try_new(
                    &changed_lines[0],
                    &changed_lines,
                    limits(),
                    || false,
                )
                .is_err()
            );
        }
    }
}

fn heat_limits() -> AngularSpectrumExposureLayerLimits {
    AngularSpectrumExposureLayerLimits {
        maximum_frequency_groups: 16,
        maximum_frequency_layer_mode_pairs: 4096,
        maximum_point_line_pairs: 4096,
        maximum_point_reservoir_frequency_pairs: 16384,
        maximum_workspace_bytes: 16 << 20,
        maximum_fft_work: 1 << 20,
    }
}

fn heat_material() -> (ResolvedDebyeCapacitor, ResolvedConductorState) {
    let basis = ElectricalPropertyBasis::Conductivity;
    let shared = source(
        &[
            (EPS_NAME, EPS, 2.0),
            (MU_NAME, MU, 1.0),
            ("delta", EPS, 0.7),
            ("tau", Time::DIMS, 0.4),
            (basis.property_name(), basis.dimensions(), 0.3),
        ],
        ValidityDomain::unconstrained(),
        QueryPoint::new(),
    );
    let conductor = ResolvedConductorState::try_from_material(shared.clone()).unwrap();
    let dielectric = ResolvedDebyeCapacitor::try_from_material(
        shared,
        EPS_NAME,
        &[DebyePoleProperties {
            permittivity_increment: "delta".into(),
            relaxation_time: "tau".into(),
        }],
        0.03,
        0.002,
    ).unwrap();
    assert_eq!(dielectric.high_frequency_permittivity_f_m(), 2.0);
    (dielectric, conductor)
}

// Independent real-time quadrature for the complete fixture below. Coefficients
// enumerate an explicitly authored common amplitude, not sampled random phases.
// Uses P=Delta E/(1+i omega tau) and actual E/H from the original interior owner,
// never the production normalized-reservoir or temporal-moment helpers.
fn heat_time_reference(
    fields: &[AngularSpectrumLayerField],
    frequencies: &[f64],
    start: f64,
    duration: f64,
    coefficients: &[C64],
) -> [f64; 3] {
    let steps = 2048;
    let mut sum = [0.0; 3];
    for &coefficient in coefficients {
        for step in 0..=steps {
            let time = start + duration * step as f64 / steps as f64;
            let mut e = [0.0; 3];
            let mut h = [0.0; 3];
            let mut ed = [0.0; 3];
            let mut hd = [0.0; 3];
            let mut p = [0.0; 3];
            let mut pd = [0.0; 3];
            for (field, &omega) in fields.iter().zip(frequencies) {
                let (sin, cos) = (omega * time).sin_cos();
                let real = |a: C64| a.re * cos - a.im * sin;
                let derivative = |a: C64| C64::new(-omega * a.im, omega * a.re);
                for c in 0..3 {
                    let a = field.field().electric_v_per_m()[0][c] * coefficient;
                    let b = field.field().magnetic_a_per_m()[0][c] * coefficient;
                    let polarization = a.scale(0.7) / C64::new(1.0, omega * 0.4);
                    e[c] += real(a);
                    h[c] += real(b);
                    ed[c] += real(derivative(a));
                    hd[c] += real(derivative(b));
                    p[c] += real(polarization);
                    pd[c] += real(derivative(polarization));
                }
            }
            let weight = if step == 0 || step == steps { 1.0 }
                else if step % 2 == 0 { 2.0 } else { 4.0 };
            for c in 0..3 {
                sum[0] += weight * 0.3 * e[c] * e[c];
                sum[1] += weight * (0.4 / 0.7) * pd[c] * pd[c];
                sum[2] += weight * (2.0 * e[c] * ed[c] + h[c] * hd[c] + p[c] * pd[c] / 0.7);
            }
        }
    }
    sum.map(|v| v / (3.0 * steps as f64 * coefficients.len() as f64))
}

#[test]
fn g1_short_window_film_heat_is_positive_while_stored_energy_flows_out() {
    // Complete SI reference fixture: host epsilon=mu=1, one 0.31 m film with
    // epsilon_inf=2, mu=1, Delta=0.7 F/m, tau=0.4 s, sigma=0.3 S/m.
    // Normal-incidence Ex peak phasors [0.8+0.2i, -0.3+0.6i] V/m at
    // omega=[2,2.3] rad/s. Window [0.17,0.29] s, sampled at z=0.155 m.
    //
    // Independent JavaScript reference, NOT execution of this Rust test:
    // 2048-step real-time Simpson quadrature and 512-step depth quadrature
    // of the separate two-interface multiple-reflection solution gave:
    // q=0.08872419055599473 W/m³, storage_rate=-2.130029947074984 W/m³;
    // integral(q dz)=0.03412361321730128 W/m²,
    // integral(storage_rate dz)=-0.5562140726132628 W/m²,
    // S_z(0)-S_z(d)=-0.522090459395133 W/m² (relative residual 1.59e-12).
    // The incorrect sum of cycle-mean line losses is 0.35214455718911153 W/m³.
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let response = DielectricResponse::DebyeWithOhmic {
        dielectric: &dielectric, conductor: &conductor,
    };
    let frequencies = [2.0, 2.3];
    let stacks = frequencies.map(|omega| PlanarDielectricLayers::try_new(
        &medium, DielectricResponse::Lossless(&medium), omega,
        &[DielectricLayer { response, thickness_m: 0.31 }],
    ).unwrap());
    let lines = [
        constant(&stacks[0], [C64::new(0.8, 0.2), C64::ZERO]),
        constant(&stacks[1], [C64::new(-0.3, 0.6), C64::ZERO]),
    ];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let layer = bank.prepare_layer(0, heat_limits(), || false).unwrap();
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    let observed = layer.field(0.155, window, || false).unwrap();
    let reference_fields = lines.iter().map(|line| line.prepare_layer(
        0, AngularSpectrumLayerLimits {
            maximum_layer_mode_pairs: 4096,
            maximum_workspace_bytes: 16 << 20,
            maximum_fft_work: 1 << 20,
        }, || false,
    ).unwrap().field(0.155, || false).unwrap()).collect::<Vec<_>>();
    let reference = heat_time_reference(&reference_fields, &frequencies, 0.17, 0.12, &[C64::ONE]);
    for value in observed.dissipation() {
        near(value.ohmic_heat_w_m3, reference[0], 2e-12);
        near(value.polarization_heat_w_m3, reference[1], 2e-12);
        near(value.storage_rate_w_m3, reference[2], 2e-11);
        near(value.total_heat_w_m3, 0.08872419055599473, 2e-11);
        assert!(value.total_heat_w_m3 > 0.0);
        assert!(value.storage_rate_w_m3 < -1.0);
    }
    let intervals = 512;
    let dz = 0.31 / intervals as f64;
    let mut heat = 0.0;
    let mut storage = 0.0;
    for i in 0..=intervals {
        let z = if i == intervals { 0.31 } else { i as f64 * dz };
        let field = layer.field(z, window, || false).unwrap();
        let weight = if i == 0 || i == intervals { 1.0 }
            else if i % 2 == 0 { 2.0 } else { 4.0 };
        heat += weight * field.absorbed_power_per_depth_w_m();
        storage += weight * field.storage_rate_per_depth_w_m();
    }
    heat *= dz / 3.0;
    storage *= dz / 3.0;
    let front = layer.field(0.0, window, || false).unwrap();
    let back = layer.field(0.31, window, || false).unwrap();
    let boundary = front.field().normal_power_w() - back.field().normal_power_w();
    assert!(heat > 0.0 && boundary < 0.0);
    near(heat + storage, boundary, 3e-10 * boundary.abs());
}

#[test]
fn g3_layer_proper_covariance_reuses_frequency_caches_and_adds_mean_heat() {
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let response = DielectricResponse::DebyeWithOhmic {
        dielectric: &dielectric, conductor: &conductor,
    };
    let frequencies = [2.0, 2.3];
    let stacks = frequencies.map(|omega| PlanarDielectricLayers::try_new(
        &medium, DielectricResponse::Lossless(&medium), omega,
        &[DielectricLayer { response, thickness_m: 0.31 }],
    ).unwrap());
    let lines = [
        constant(&stacks[0], [C64::new(0.8, 0.2), C64::ZERO]),
        constant(&stacks[1], [C64::new(-0.3, 0.6), C64::ZERO]),
    ];
    let mean = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let factors = [AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap()];
    let ensemble = AngularSpectrumExposureEnsemble::try_new(
        &lines[0], Some(&mean), &factors, ensemble_limits(), || false,
    ).unwrap();
    let layer = ensemble.prepare_layer(0, heat_limits(), || false).unwrap();
    assert_eq!(layer.frequency_count(), 2);
    // Two frequencies, four grid modes, two polarizations, one source layer:
    // the mean and covariance factor do not duplicate the constitutive solve.
    assert_eq!(layer.preparation_layer_mode_pairs(), 16);
    assert_eq!(layer.field_fft_work(), 4 * lines[0].field_fft_work());
    let fields = lines.iter().map(|line| line.prepare_layer(
        0, AngularSpectrumLayerLimits {
            maximum_layer_mode_pairs: 4096,
            maximum_workspace_bytes: 16 << 20,
            maximum_fft_work: 1 << 20,
        }, || false,
    ).unwrap().field(0.155, || false).unwrap()).collect::<Vec<_>>();
    // Mean plus one proper factor: 1+z with z equally likely 1,-1,i,-i.
    let reference = heat_time_reference(&fields, &frequencies, 0.17, 0.12, &[
        C64::from_re(2.0), C64::ZERO, C64::new(1.0, 1.0), C64::new(1.0, -1.0),
    ]);
    let observed = layer.field(0.155, HarmonicWindow::new(0.17, 0.12).unwrap(), || false).unwrap();
    for value in observed.dissipation() {
        near(value.ohmic_heat_w_m3, reference[0], 3e-12);
        near(value.polarization_heat_w_m3, reference[1], 3e-12);
        near(value.storage_rate_w_m3, reference[2], 3e-11);
    }
}

#[test]
fn g3_single_tone_finite_layer_heat_agrees_with_original_cycle_mean_absorption() {
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let response = DielectricResponse::DebyeWithOhmic {
        dielectric: &dielectric, conductor: &conductor,
    };
    let stack = PlanarDielectricLayers::try_new(
        &medium, DielectricResponse::Lossless(&medium), 2.0,
        &[DielectricLayer { response, thickness_m: 0.31 }],
    ).unwrap();
    let lines = [constant(&stack, [C64::new(0.8, 0.2), C64::new(-0.1, 0.3)])];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let layer = bank.prepare_layer(0, heat_limits(), || false).unwrap();
    let original = lines[0].prepare_layer(0, AngularSpectrumLayerLimits {
        maximum_layer_mode_pairs: 4096,
        maximum_workspace_bytes: 16 << 20,
        maximum_fft_work: 1 << 20,
    }, || false).unwrap();
    let window = HarmonicWindow::new(0.17, core::f64::consts::PI).unwrap();
    let intervals = 256;
    let dz = 0.31 / intervals as f64;
    let mut integrated_heat = 0.0;
    for i in 0..=intervals {
        let z = if i == intervals { 0.31 } else { i as f64 * dz };
        let observed = layer.field(z, window, || false).unwrap();
        let reference = original.field(z, || false).unwrap();
        for (point, expected) in observed.dissipation().iter()
            .zip(reference.absorbed_power_density_w_m3())
        {
            near(point.total_heat_w_m3, *expected, 3e-12 * expected.abs().max(1e-15));
            near(point.storage_rate_w_m3, 0.0, 3e-12);
        }
        let weight = if i == 0 || i == intervals { 1.0 }
            else if i % 2 == 0 { 2.0 } else { 4.0 };
        integrated_heat += weight * observed.absorbed_power_per_depth_w_m();
    }
    integrated_heat *= dz / 3.0;
    near(integrated_heat, original.absorbed_power_w(), 2e-9 * original.absorbed_power_w());
}

#[test]
fn g4_layer_exposure_budgets_depth_and_cancellation_remain_explicit() {
    let medium = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    let response = DielectricResponse::DebyeWithOhmic {
        dielectric: &dielectric, conductor: &conductor,
    };
    let stacks = [2.0, 2.3].map(|omega| PlanarDielectricLayers::try_new(
        &medium, DielectricResponse::Lossless(&medium), omega,
        &[DielectricLayer { response, thickness_m: 0.31 }],
    ).unwrap());
    let lines = stacks.iter().map(|stack| constant(stack, [C64::ONE, C64::ZERO]))
        .collect::<Vec<_>>();
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let admitted = bank.prepare_layer(0, heat_limits(), || false).unwrap();
    let exact = AngularSpectrumExposureLayerLimits {
        maximum_frequency_groups: 2,
        maximum_frequency_layer_mode_pairs: admitted.preparation_layer_mode_pairs(),
        maximum_point_line_pairs: admitted.point_line_pairs(),
        maximum_point_reservoir_frequency_pairs: admitted.point_reservoir_frequency_pairs(),
        maximum_workspace_bytes: admitted.workspace_bytes(),
        maximum_fft_work: admitted.field_fft_work(),
    };
    assert!(bank.prepare_layer(0, exact, || false).is_ok());
    for cap in [
        AngularSpectrumExposureLayerLimits { maximum_frequency_groups: 1, ..exact },
        AngularSpectrumExposureLayerLimits {
            maximum_frequency_layer_mode_pairs: exact.maximum_frequency_layer_mode_pairs - 1, ..exact
        },
        AngularSpectrumExposureLayerLimits {
            maximum_point_line_pairs: exact.maximum_point_line_pairs - 1, ..exact
        },
        AngularSpectrumExposureLayerLimits {
            maximum_point_reservoir_frequency_pairs: exact.maximum_point_reservoir_frequency_pairs - 1,
            ..exact
        },
        AngularSpectrumExposureLayerLimits {
            maximum_workspace_bytes: exact.maximum_workspace_bytes - 1, ..exact
        },
        AngularSpectrumExposureLayerLimits { maximum_fft_work: exact.maximum_fft_work - 1, ..exact },
    ] {
        assert!(bank.prepare_layer(0, cap, || false).is_err());
    }
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    for depth in [-1.0, 0.32, f64::NAN] {
        assert!(admitted.field(depth, window, || false).is_err());
    }
    let mut polls = 0;
    assert!(matches!(admitted.field(0.15, window, || {
        polls += 1;
        polls > 5
    }), Err(AngularSpectrumError::Cancelled)));
    let retry = admitted.field(0.15, window, || false).unwrap();
    let replay = admitted.field(0.15, window, || false).unwrap();
    assert_eq!(retry.field().points(), replay.field().points());
    assert_eq!(retry.dissipation(), replay.dissipation());
}
