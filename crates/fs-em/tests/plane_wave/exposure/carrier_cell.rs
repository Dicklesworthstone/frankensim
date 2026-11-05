//! Continuous cell and detector interference with an authored transverse carrier.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumAperture, AngularSpectrumExposureApertureLimits, AngularSpectrumExposureField,
    AngularSpectrumExposurePoint, AngularSpectrumExposureStressMomentumField,
    AngularSpectrumExposureStressMomentumLimits, AngularSpectrumFactorStatistics,
};
use fs_fft::harmonic_window::HarmonicWindowProfile;

const FREQUENCIES: [f64; 2] = [2.3, 3.7];
const CARRIER: [f64; 2] = [0.23, 0.17];

fn amplitudes() -> [C64; 2] {
    [C64::new(1.1, 0.2), C64::new(-0.3, 0.8)]
}

fn carrier_line<'s, 'm>(
    stack: &'s PlanarDielectricLayers<'m>,
    carrier: [f64; 2],
    amplitude: C64,
) -> PlanarAngularSpectrum<'s, 'm> {
    let kt = (carrier[0] * carrier[0] + carrier[1] * carrier[1]).sqrt();
    let electric = if kt == 0.0 {
        [C64::ZERO, amplitude]
    } else {
        [
            amplitude.scale(-carrier[1] / kt),
            amplitude.scale(carrier[0] / kt),
        ]
    };
    PlanarAngularSpectrum::prepare_with_carrier(
        stack,
        grid(),
        carrier,
        &[electric; 4],
        AngularSpectrumLimits {
            maximum_grid_points: 128,
            maximum_workspace_bytes: 1 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap()
}

fn sinc(x: f64) -> f64 {
    if x == 0.0 { 1.0 } else { x.sin() / x }
}

fn phase(angle: f64) -> C64 {
    let (s, c) = angle.sin_cos();
    C64::new(c, s)
}

fn positive_time_mean(omega: f64, window: HarmonicWindow) -> C64 {
    let response = match window.profile() {
        HarmonicWindowProfile::Uniform => sinc(0.5 * omega * window.duration_s()),
        HarmonicWindowProfile::Triangular => sinc(0.25 * omega * window.duration_s()).powi(2),
    };
    phase(omega * (window.start_s() + 0.5 * window.duration_s())).scale(response)
}

struct Reference {
    point: AngularSpectrumExposurePoint,
    stress: [[f64; 3]; 3],
    momentum_rate: [f64; 3],
}

// Independent positive-time product identity, with no native harmonic kernel,
// reconstructed field, Gram matrix or spatial quadrature. Both waves share
// the transverse carrier, so their difference-wavevector factor is exactly 1.
// The sum term averages exp(-2i*k.x) over the actual physical rectangle.
// weights = [deterministic/real contributions, proper contributions].
fn reference(
    window: HarmonicWindow,
    origin: [f64; 2],
    extent: [f64; 2],
    weights: [f64; 2],
) -> Reference {
    let kt = (CARRIER[0] * CARRIER[0] + CARRIER[1] * CARRIER[1]).sqrt();
    let te = [-CARRIER[1] / kt, CARRIER[0] / kt];
    let fields: [[C64; 6]; 2] = core::array::from_fn(|f| {
        let amplitude = amplitudes()[f];
        let omega = FREQUENCIES[f];
        let kz = (omega * omega - kt * kt).sqrt();
        [
            amplitude.scale(te[0]),
            amplitude.scale(te[1]),
            C64::ZERO,
            amplitude.scale(-kz * te[1] / omega),
            amplitude.scale(kz * te[0] / omega),
            amplitude.scale(kt / omega),
        ]
    });
    let spatial_sum = phase(
        -2.0 * (CARRIER[0] * (origin[0] + 0.5 * extent[0])
            + CARRIER[1] * (origin[1] + 0.5 * extent[1])),
    )
    .scale(sinc(CARRIER[0] * extent[0]) * sinc(CARRIER[1] * extent[1]));
    let mut second = [[0.0; 6]; 6];
    let mut rate = [[0.0; 6]; 6];
    for a in 0..2 {
        for b in 0..2 {
            let difference_frequency = FREQUENCIES[a] - FREQUENCIES[b];
            let sum_frequency = FREQUENCIES[a] + FREQUENCIES[b];
            let difference =
                positive_time_mean(difference_frequency, window).scale(weights[0] + weights[1]);
            let sum = positive_time_mean(sum_frequency, window) * spatial_sum;
            for i in 0..6 {
                for j in 0..6 {
                    let d = fields[a][i] * fields[b][j].conj() * difference;
                    let s = (fields[a][i] * fields[b][j] * sum).scale(weights[0]);
                    second[i][j] += 0.5 * (d.re + s.re);
                    // The window weights the instantaneous derivative. A
                    // triangular shutter is not an endpoint-difference rate.
                    rate[i][j] += 0.5
                        * ((C64::new(0.0, difference_frequency) * d).re
                            + (C64::new(0.0, sum_frequency) * s).re);
                }
            }
        }
    }
    let electric = core::array::from_fn(|i| second[i][i]);
    let magnetic = core::array::from_fn(|i| second[i + 3][i + 3]);
    let poynting = core::array::from_fn(|i| {
        let (a, b) = ((i + 1) % 3, (i + 2) % 3);
        second[a][b + 3] - second[b][a + 3]
    });
    let momentum_rate = core::array::from_fn(|i| {
        let (a, b) = ((i + 1) % 3, (i + 2) % 3);
        rate[a][b + 3] - rate[b][a + 3]
    });
    // This original fixture has epsilon=mu=1: g=D cross B=E cross H.
    let trace: f64 = electric.iter().chain(&magnetic).sum();
    let stress = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            second[i][j] + second[i + 3][j + 3] - if i == j { 0.5 * trace } else { 0.0 }
        })
    });
    Reference {
        point: AngularSpectrumExposurePoint {
            electric_mean_squared_v2_m2: electric,
            magnetic_mean_squared_a2_m2: magnetic,
            poynting_w_m2: poynting,
        },
        stress,
        momentum_rate,
    }
}

fn close(actual: f64, expected: f64) {
    near(actual, expected, 5e-12 * expected.abs().max(1.0));
}

fn compare_point(actual: AngularSpectrumExposurePoint, expected: AngularSpectrumExposurePoint) {
    for (actual, expected) in actual
        .electric_mean_squared_v2_m2
        .into_iter()
        .chain(actual.magnetic_mean_squared_a2_m2)
        .chain(actual.poynting_w_m2)
        .zip(
            expected
                .electric_mean_squared_v2_m2
                .into_iter()
                .chain(expected.magnetic_mean_squared_a2_m2)
                .chain(expected.poynting_w_m2),
        )
    {
        close(actual, expected);
    }
}

fn check_cell(
    field: &AngularSpectrumExposureField,
    joint: &AngularSpectrumExposureStressMomentumField,
    window: HarmonicWindow,
    weights: [f64; 2],
) {
    let cell = grid();
    let extent = [
        cell.width as f64 * cell.spacing_m[0],
        cell.height as f64 * cell.spacing_m[1],
    ];
    let area = extent[0] * extent[1];
    let expected = reference(window, [0.0; 2], extent, weights);
    close(
        field.normal_power_w(),
        area * expected.point.poynting_w_m2[2],
    );
    for (index, point) in field.points().iter().enumerate() {
        let origin = [
            (index % cell.width) as f64 * cell.spacing_m[0],
            (index / cell.width) as f64 * cell.spacing_m[1],
        ];
        let at = reference(window, origin, [0.0; 2], weights);
        compare_point(*point, at.point);
        for i in 0..3 {
            for j in 0..3 {
                close(
                    joint.stress().points()[index].maxwell_stress_pa[i][j],
                    at.stress[i][j],
                );
            }
            close(
                joint.momentum().points()[index].minkowski_density_n_s_m3[i],
                at.point.poynting_w_m2[i],
            );
            close(
                joint.momentum().points()[index].minkowski_rate_n_m3[i],
                at.momentum_rate[i],
            );
        }
    }
    for i in 0..3 {
        for j in 0..3 {
            close(
                joint.stress().integrated_stress_n()[i][j],
                area * expected.stress[i][j],
            );
        }
        close(
            joint.momentum().momentum_per_depth_n_s_m()[i],
            area * expected.point.poynting_w_m2[i],
        );
        close(
            joint.momentum().momentum_rate_per_depth_n_m()[i],
            area * expected.momentum_rate[i],
        );
    }
}

fn aperture_limits() -> AngularSpectrumExposureApertureLimits {
    AngularSpectrumExposureApertureLimits {
        maximum_mode_frequency_pairs: 4096,
        maximum_workspace_bytes: 8 << 20,
    }
}

#[test]
fn g1_common_carrier_real_fields_match_continuous_cell_and_rectangle_integrals() {
    let medium = material(1.0, 1.0);
    let stacks = FREQUENCIES.map(|omega| {
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap()
    });
    let lines: Vec<_> = stacks
        .iter()
        .zip(amplitudes())
        .map(|(stack, a)| carrier_line(stack, CARRIER, a))
        .collect();
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(-0.17, 0.19, profile).unwrap();
        let observed = bank.field(Side::Incident, 0.0, window, || false).unwrap();
        let joint = bank
            .stress_momentum_field(
                Side::Incident,
                0.0,
                window,
                AngularSpectrumExposureStressMomentumLimits {
                    maximum_workspace_bytes: bank.stress_momentum_workspace_bytes().unwrap(),
                },
                || false,
            )
            .unwrap();
        check_cell(&observed, &joint, window, [1.0, 0.0]);
        let sampled = observed
            .points()
            .iter()
            .map(|p| p.poynting_w_m2[2])
            .sum::<f64>()
            * grid().spacing_m[0]
            * grid().spacing_m[1];
        // Negative control: the carrier's real quadratic field is not periodic
        // on the envelope cell. Sampling then summing gives a ~20% wrong power.
        assert!(
            (sampled - observed.normal_power_w()).abs() > 0.1 * observed.normal_power_w().abs()
        );
        for (origin, extent) in [([0.0; 2], [32.0, 3.0]), ([-3.2, 1.7], [11.3, 2.6])] {
            let aperture = AngularSpectrumAperture::try_new(origin, extent).unwrap();
            let measured = bank
                .measure_aperture(
                    Side::Incident,
                    0.0,
                    aperture,
                    window,
                    aperture_limits(),
                    || false,
                )
                .unwrap();
            let expected = reference(window, origin, extent, [1.0, 0.0]);
            compare_point(measured.mean, expected.point);
            close(
                measured.normal_power_w,
                aperture.area_m2() * expected.point.poynting_w_m2[2],
            );
        }
    }
}

#[test]
fn g1_common_carrier_real_and_proper_factors_keep_their_distinct_sum_terms() {
    let medium = material(1.0, 1.0);
    let stacks = FREQUENCIES.map(|omega| {
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap()
    });
    let lines: Vec<_> = stacks
        .iter()
        .zip(amplitudes())
        .map(|(stack, a)| carrier_line(stack, CARRIER, a))
        .collect();
    let mean = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let factors = [
        AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap(),
        AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap(),
    ];
    let statistics = [
        AngularSpectrumFactorStatistics::Real,
        AngularSpectrumFactorStatistics::ProperComplex,
    ];
    let ensemble = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &lines[0],
        Some(&mean),
        &factors,
        &statistics,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let proper = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &lines[0],
        None,
        &factors[1..],
        &statistics[1..],
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let window =
        HarmonicWindow::with_profile(0.13, 0.37, HarmonicWindowProfile::Triangular).unwrap();
    let observed = ensemble
        .field(Side::Incident, 0.0, window, || false)
        .unwrap();
    let joint = ensemble
        .stress_momentum_field(
            Side::Incident,
            0.0,
            window,
            AngularSpectrumExposureStressMomentumLimits {
                maximum_workspace_bytes: ensemble.stress_momentum_workspace_bytes().unwrap(),
            },
            || false,
        )
        .unwrap();
    check_cell(&observed, &joint, window, [2.0, 1.0]);
    let origin = [-3.2, 1.7];
    let extent = [11.3, 2.6];
    let aperture = AngularSpectrumAperture::try_new(origin, extent).unwrap();
    let measured = ensemble
        .measure_aperture(
            Side::Incident,
            0.0,
            aperture,
            window,
            aperture_limits(),
            || false,
        )
        .unwrap();
    let expected = reference(window, origin, extent, [2.0, 1.0]);
    compare_point(measured.mean, expected.point);
    close(
        measured.normal_power_w,
        aperture.area_m2() * expected.point.poynting_w_m2[2],
    );
    let proper_field = proper.field(Side::Incident, 0.0, window, || false).unwrap();
    let proper_reference = reference(window, [0.0; 2], [32.0, 3.0], [0.0, 1.0]);
    for &point in proper_field.points() {
        // Proper statistics remove the common spatial carrier from every
        // second moment, while the two frequencies still beat coherently.
        compare_point(point, proper_reference.point);
    }
    close(
        proper_field.normal_power_w(),
        96.0 * proper_reference.point.poynting_w_m2[2],
    );
}

#[test]
fn g4_carrier_work_identity_and_cancellation_are_explicit_without_changing_zero_carrier_caps() {
    let medium = material(1.0, 1.0);
    let stacks = FREQUENCIES.map(|omega| {
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap()
    });
    let zero: Vec<_> = stacks
        .iter()
        .zip(amplitudes())
        .map(|(stack, a)| carrier_line(stack, [0.0; 2], a))
        .collect();
    let lines: Vec<_> = stacks
        .iter()
        .zip(amplitudes())
        .map(|(stack, a)| carrier_line(stack, CARRIER, a))
        .collect();
    let n = grid().width * grid().height;
    let frequencies = FREQUENCIES.len();
    let old_pairs = n * frequencies * frequencies;
    let carrier_pairs = old_pairs + frequencies * n * (grid().width + grid().height);
    let mut cap = limits();
    cap.maximum_point_frequency_pairs = old_pairs;
    let zero_bank = AngularSpectrumExposure::try_new(&zero[0], &zero, cap, || false).unwrap();
    assert_eq!(zero_bank.point_frequency_pairs(), old_pairs);
    assert!(AngularSpectrumExposure::try_new(&lines[0], &lines, cap, || false).is_err());
    cap.maximum_point_frequency_pairs = carrier_pairs;
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, cap, || false).unwrap();
    assert_eq!(bank.point_frequency_pairs(), carrier_pairs);
    cap.maximum_workspace_bytes = bank.workspace_bytes() - 1;
    assert!(AngularSpectrumExposure::try_new(&lines[0], &lines, cap, || false).is_err());
    cap.maximum_workspace_bytes += 1;
    assert!(AngularSpectrumExposure::try_new(&lines[0], &lines, cap, || false).is_ok());
    let other_carrier = [CARRIER[0] + 0.01, CARRIER[1]];
    let mixed = [
        carrier_line(&stacks[0], CARRIER, amplitudes()[0]),
        carrier_line(&stacks[1], other_carrier, amplitudes()[1]),
    ];
    assert!(AngularSpectrumExposure::try_new(&lines[0], &mixed, limits(), || false).is_err());
    let other_lines: Vec<_> = stacks
        .iter()
        .zip(amplitudes())
        .map(|(stack, a)| carrier_line(stack, other_carrier, a))
        .collect();
    let wrong_factors =
        [
            AngularSpectrumExposure::try_new(&other_lines[0], &other_lines, limits(), || false)
                .unwrap(),
        ];
    assert!(
        AngularSpectrumExposureEnsemble::try_new(
            &lines[0],
            Some(&bank),
            &wrong_factors,
            ensemble_limits(),
            || false,
        )
        .is_err()
    );
    assert!(matches!(
        AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    let window = HarmonicWindow::new(-0.17, 0.19).unwrap();
    let proper_law = [AngularSpectrumFactorStatistics::ProperComplex];
    let mut proper_cap = ensemble_limits();
    proper_cap.maximum_point_frequency_pairs = old_pairs;
    let proper = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &lines[0],
        None,
        core::slice::from_ref(&bank),
        &proper_law,
        proper_cap,
        || false,
    )
    .unwrap();
    assert_eq!(proper.point_frequency_pairs(), old_pairs);
    assert!(
        AngularSpectrumExposureEnsemble::try_new_with_statistics(
            &lines[0],
            None,
            core::slice::from_ref(&bank),
            &[AngularSpectrumFactorStatistics::Real],
            proper_cap,
            || false,
        )
        .is_err()
    );
    assert!(
        AngularSpectrumExposureEnsemble::try_new_with_statistics(
            &lines[0],
            Some(&bank),
            &[],
            &[],
            proper_cap,
            || false,
        )
        .is_err()
    );
    assert!(proper.field(Side::Incident, 0.0, window, || false).is_ok());
    let before = bank.field(Side::Incident, 0.0, window, || false).unwrap();
    let mut polls = 0;
    assert!(matches!(
        bank.field(Side::Incident, 0.0, window, || {
            polls += 1;
            polls == 9
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    let after = bank.field(Side::Incident, 0.0, window, || false).unwrap();
    assert_eq!(before.points(), after.points());
    assert_eq!(before.normal_power_w(), after.normal_power_w());
    assert!(
        zero_bank
            .field(Side::Incident, 0.0, window, || false)
            .is_ok()
    );
}

#[test]
fn g1_carrier_source_control_matches_actual_energy_gradient_and_hessian_and_rejects_foreign_carriers()
 {
    use fs_em::plane_wave::angular_spectrum::{
        AngularSpectrumApertureControlLimits, PreparedAngularSpectrumApertureEnergyControl,
    };

    let medium = material(1.0, 1.0);
    let stacks = FREQUENCIES.map(|omega| {
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap()
    });
    let source_amplitudes = [C64::new(0.4, 0.6), C64::new(-0.7, 0.2)];
    let fixed_lines: Vec<_> = stacks
        .iter()
        .zip(amplitudes())
        .map(|(stack, a)| carrier_line(stack, CARRIER, a))
        .collect();
    let source_lines: Vec<_> = stacks
        .iter()
        .zip(source_amplitudes)
        .map(|(stack, a)| carrier_line(stack, CARRIER, a))
        .collect();
    let fixed = AngularSpectrumExposure::try_new(&fixed_lines[0], &fixed_lines, limits(), || false)
        .unwrap();
    let controls =
        [
            AngularSpectrumExposure::try_new(&fixed_lines[0], &source_lines, limits(), || false)
                .unwrap(),
        ];
    let aperture = AngularSpectrumAperture::try_new([-3.2, 1.7], [11.3, 2.6]).unwrap();
    let window =
        HarmonicWindow::with_profile(0.13, 0.37, HarmonicWindowProfile::Triangular).unwrap();
    let control_limits = AngularSpectrumApertureControlLimits {
        maximum_controls: 4,
        maximum_point_line_pairs: 4096,
        maximum_source_field_entries: 4096,
        maximum_kernel_pairs: 4096,
        maximum_workspace_bytes: 16 << 20,
    };
    let actual = |coefficient: C64| {
        let lines: Vec<_> = stacks
            .iter()
            .enumerate()
            .map(|(i, stack)| {
                carrier_line(
                    stack,
                    CARRIER,
                    amplitudes()[i] + coefficient * source_amplitudes[i],
                )
            })
            .collect();
        let bank =
            AngularSpectrumExposure::try_new(&fixed_lines[0], &lines, limits(), || false).unwrap();
        bank.measure_aperture_stress(
            Side::Incident,
            0.0,
            aperture,
            window,
            aperture_limits(),
            || false,
        )
        .unwrap()
        .mean_field_energy_density_j_m3()
    };
    // Exact centered identities for a real quadratic in complex Cartesian
    // source coordinates. Finite unit steps avoid small-difference roundoff.
    let actual_gradient = |coefficient: C64| {
        C64::new(
            0.5 * (actual(coefficient + C64::ONE) - actual(coefficient - C64::ONE)),
            0.5 * (actual(coefficient + C64::new(0.0, 1.0))
                - actual(coefficient - C64::new(0.0, 1.0))),
        )
    };
    let mut prepared = PreparedAngularSpectrumApertureEnergyControl::new(
        &fixed,
        &controls,
        Side::Incident,
        0.0,
        aperture,
        window,
        control_limits,
        || false,
    )
    .unwrap();
    let coefficient = C64::new(0.35, -0.4);
    let result = prepared.value_gradient(&[coefficient], || false).unwrap();
    close(result.mean_field_energy_density_j_m3, actual(coefficient));
    let expected_gradient = actual_gradient(coefficient);
    close(result.source_gradient[0].re, expected_gradient.re);
    close(result.source_gradient[0].im, expected_gradient.im);

    let direction = C64::new(0.3, 0.8);
    let hessian = prepared.hessian_vector(&[direction], || false).unwrap()[0];
    let expected_hessian = (actual_gradient(coefficient + direction)
        - actual_gradient(coefficient - direction))
    .scale(0.5);
    close(hessian.re, expected_hessian.re);
    close(hessian.im, expected_hessian.im);

    // This is a valid source bank in its own original carrier frame, with
    // identical materials, layers, grid and frequency support. It must not
    // become a control in the fixed bank's different carrier frame.
    let different_carrier = [CARRIER[0] + 0.01, CARRIER[1]];
    let other_lines: Vec<_> = stacks
        .iter()
        .zip(source_amplitudes)
        .map(|(stack, a)| carrier_line(stack, different_carrier, a))
        .collect();
    let foreign =
        [
            AngularSpectrumExposure::try_new(&other_lines[0], &other_lines, limits(), || false)
                .unwrap(),
        ];
    assert!(
        PreparedAngularSpectrumApertureEnergyControl::new(
            &fixed,
            &foreign,
            Side::Incident,
            0.0,
            aperture,
            window,
            control_limits,
            || false,
        )
        .is_err()
    );
}

#[test]
fn g1_two_axis_signed_carrier_modes_match_dense_cell_power_stress_and_momentum() {
    use fs_em::plane_wave::angular_spectrum::AngularSpectrumExposureMomentumLimits;

    let medium = material(1.0, 1.0);
    let stacks = FREQUENCIES.map(|omega| {
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap()
    });
    let cell = AngularSpectrumGrid {
        width: 4,
        height: 4,
        spacing_m: [2.5, 2.1],
    };
    let extent = [
        cell.width as f64 * cell.spacing_m[0],
        cell.height as f64 * cell.spacing_m[1],
    ];
    // Authored Fourier modes occupy both signs in BOTH axes; -2 is the
    // negative Nyquist bin, not an inferred positive partner. All bins of this
    // small grid remain inside the lowest-frequency incident light cone.
    let modes = [
        ([0_i32, 0_i32], C64::new(0.7, 0.15)),
        ([1, -1], C64::new(-0.35, 0.4)),
        ([-1, 1], C64::new(0.25, -0.3)),
        ([-2, -2], C64::new(0.45, 0.2)),
    ];
    let lines: Vec<_> = stacks
        .iter()
        .enumerate()
        .map(|(frequency, stack)| {
            let samples: Vec<_> = (0..cell.width * cell.height)
                .map(|point| {
                    let (column, row) = (point % cell.width, point / cell.width);
                    let mut electric = [C64::ZERO; 2];
                    for (mode, &(bins, amplitude)) in modes.iter().enumerate() {
                        let kx = CARRIER[0] + core::f64::consts::TAU * bins[0] as f64 / extent[0];
                        let ky = CARRIER[1] + core::f64::consts::TAU * bins[1] as f64 / extent[1];
                        let kt = (kx * kx + ky * ky).sqrt();
                        let envelope = phase(
                            -core::f64::consts::TAU
                                * (bins[0] as f64 * column as f64 / cell.width as f64
                                    + bins[1] as f64 * row as f64 / cell.height as f64),
                        );
                        let signal = (amplitude * amplitudes()[frequency] * envelope)
                            .scale(1.0 + 0.17 * frequency as f64 * mode as f64);
                        electric[0] = electric[0] + signal.scale(-ky / kt);
                        electric[1] = electric[1] + signal.scale(kx / kt);
                    }
                    electric
                })
                .collect();
            PlanarAngularSpectrum::prepare_with_carrier(
                stack,
                cell,
                CARRIER,
                &samples,
                AngularSpectrumLimits {
                    maximum_grid_points: 128,
                    maximum_workspace_bytes: 1 << 20,
                    maximum_fft_work: 1 << 20,
                },
                || false,
            )
            .unwrap()
        })
        .collect();
    let factors =
        [AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap()];
    let aperture = AngularSpectrumAperture::try_new([0.0; 2], extent).unwrap();
    let area = aperture.area_m2();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.13, 0.37, profile).unwrap();
        let mut powers = [0.0; 2];
        for (law_index, law) in [
            AngularSpectrumFactorStatistics::Real,
            AngularSpectrumFactorStatistics::ProperComplex,
        ]
        .into_iter()
        .enumerate()
        {
            let laws = [law];
            let mut cap = ensemble_limits();
            cap.maximum_point_frequency_pairs = if law == AngularSpectrumFactorStatistics::Real {
                factors[0].point_frequency_pairs()
            } else {
                cell.width * cell.height * FREQUENCIES.len() * FREQUENCIES.len()
            };
            let ensemble = AngularSpectrumExposureEnsemble::try_new_with_statistics(
                &lines[0],
                None,
                &factors,
                &laws,
                cap,
                || false,
            )
            .unwrap();
            let field = ensemble
                .field(Side::Incident, 0.0, window, || false)
                .unwrap();
            let joint = ensemble
                .stress_momentum_field(
                    Side::Incident,
                    0.0,
                    window,
                    AngularSpectrumExposureStressMomentumLimits {
                        maximum_workspace_bytes: ensemble
                            .stress_momentum_workspace_bytes()
                            .unwrap(),
                    },
                    || false,
                )
                .unwrap();
            let momentum = ensemble
                .momentum_field(
                    Side::Incident,
                    0.0,
                    window,
                    AngularSpectrumExposureMomentumLimits {
                        maximum_workspace_bytes: ensemble.momentum_workspace_bytes().unwrap(),
                    },
                    || false,
                )
                .unwrap();
            // The preserved dense N²F² rectangle integrator is an independent
            // spatial algorithm; no tensor-product cell code builds this oracle.
            let dense = ensemble
                .measure_aperture_stress(
                    Side::Incident,
                    0.0,
                    aperture,
                    window,
                    aperture_limits(),
                    || false,
                )
                .unwrap();
            powers[law_index] = dense.normal_power_w().unwrap();
            close(field.normal_power_w(), powers[law_index]);
            let expected_stress = dense.integrated_stress_n().unwrap();
            let expected_flux = dense.mean_poynting_w_m2();
            if law == AngularSpectrumFactorStatistics::Real {
                let deterministic = factors[0]
                    .field(Side::Incident, 0.0, window, || false)
                    .unwrap();
                close(deterministic.normal_power_w(), powers[law_index]);
            }
            // Translating a fixed-profile shutter differentiates its average
            // into the weighted instantaneous rate, for BOTH window profiles.
            // This checks the rate through independent dense flux evaluations.
            let h = 2e-6;
            let before = ensemble
                .measure_aperture_stress(
                    Side::Incident,
                    0.0,
                    aperture,
                    HarmonicWindow::with_profile(
                        window.start_s() - h,
                        window.duration_s(),
                        profile,
                    )
                    .unwrap(),
                    aperture_limits(),
                    || false,
                )
                .unwrap()
                .mean_poynting_w_m2();
            let after = ensemble
                .measure_aperture_stress(
                    Side::Incident,
                    0.0,
                    aperture,
                    HarmonicWindow::with_profile(
                        window.start_s() + h,
                        window.duration_s(),
                        profile,
                    )
                    .unwrap(),
                    aperture_limits(),
                    || false,
                )
                .unwrap()
                .mean_poynting_w_m2();
            for i in 0..3 {
                for j in 0..3 {
                    close(
                        joint.stress().integrated_stress_n()[i][j],
                        expected_stress[i][j],
                    );
                }
                // epsilon=mu=1 makes the actual Minkowski density equal to S.
                close(
                    joint.momentum().momentum_per_depth_n_s_m()[i],
                    area * expected_flux[i],
                );
                close(
                    momentum.momentum_per_depth_n_s_m()[i],
                    area * expected_flux[i],
                );
                let expected_rate = area * (after[i] - before[i]) / (2.0 * h);
                near(
                    joint.momentum().momentum_rate_per_depth_n_m()[i],
                    expected_rate,
                    2e-6 * expected_rate.abs().max(1.0),
                );
                close(
                    momentum.momentum_rate_per_depth_n_m()[i],
                    joint.momentum().momentum_rate_per_depth_n_m()[i],
                );
            }
        }
        assert!((powers[0] - powers[1]).abs() > 1e-5);
    }
}
