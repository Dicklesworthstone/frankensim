use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

type Law = HarmonicExposureFactorStatistics;
type Coefficients = [[C64; 2]; 2];

fn limits() -> HarmonicExposureLimits {
    HarmonicExposureLimits {
        maximum_lines: 32,
        maximum_factors: 8,
        maximum_transfer_products: 8192,
        maximum_receiver_frequency_pairs: 8192,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn array_limits() -> HarmonicPressureArrayLimits {
    HarmonicPressureArrayLimits {
        maximum_channels: 8,
        maximum_taps: 32,
        maximum_prepare_products: 8192,
        maximum_query_products: 8192,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn readout_limits() -> HarmonicReadoutLimits {
    HarmonicReadoutLimits {
        maximum_receiver_pairs: 4,
        maximum_products: 8192,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn owner(omega: f64, frequency: usize) -> PreparedHarmonicObservation {
    PreparedHarmonicObservation {
        omega,
        medium: crate::helmholtz::Medium {
            density: 1.2,
            sound_speed: 343.0,
        },
        sources: vec![
            SourceId::StructuralMode("structural-mode".into()),
            SourceId::PointVolumeFlow("original-volume-flow".into()),
        ],
        receivers: (0..2)
            .map(|r| ReceiverSpec {
                id: format!("physical-{r}"),
                position_m: [0.2 * r as f64, 0.3, -0.1],
            })
            .collect(),
        producer_hash: ContentHash([41; 32]),
        // Synthetic original transfer rows exercise the full complex p/u
        // algebra. The time oracle below never calls a harmonic kernel.
        transfers: (0..4)
            .map(|row| {
                core::array::from_fn(|c| {
                    let q = (1 + row + 2 * c + frequency) as f64;
                    C64::new(0.07 * q - 0.23, 0.04 * q - 0.17 * c as f64)
                })
            })
            .collect(),
    }
}

struct Fixture {
    owners: [PreparedHarmonicObservation; 2],
    coefficients: [Coefficients; 3],
}

impl Fixture {
    fn new() -> Self {
        Self {
            owners: [owner(5.0, 0), owner(9.0, 1)],
            coefficients: [
                [
                    [C64::new(0.7, 0.2), C64::new(-0.1, 0.3)],
                    [C64::new(-0.4, 0.5), C64::new(0.2, -0.6)],
                ],
                [
                    [C64::new(0.3, -0.1), C64::new(0.8, 0.2)],
                    [C64::new(0.6, 0.8), C64::new(-0.3, 0.1)],
                ],
                [
                    [C64::new(-0.2, 0.9), C64::new(0.4, 0.6)],
                    [C64::new(0.4, -0.3), C64::new(-0.5, -0.2)],
                ],
            ],
        }
    }

    fn signal(&self, state: usize, receiver: usize, component: usize, time: f64) -> f64 {
        let real = if state % 2 == 0 { 1.0 } else { -1.0 };
        let phase = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)][state / 2];
        let mut value = 0.0;
        for f in 0..2 {
            let (sin, cos) = (self.owners[f].omega * time).sin_cos();
            for s in 0..2 {
                let coefficient = self.coefficients[0][f][s]
                    + self.coefficients[1][f][s] * C64::new(real, 0.0)
                    + self.coefficients[2][f][s] * phase;
                let a = self.owners[f].transfers[2 * s + receiver][component] * coefficient;
                value += a.re * cos + a.im * sin;
            }
        }
        value
    }
}

fn prepare<'a>(
    owners: &'a [PreparedHarmonicObservation; 2],
    coefficients: &[Coefficients],
    statistics: &[Law],
) -> PreparedHarmonicExposure<'a> {
    let banks: Vec<[HarmonicExposureLine<'_>; 2]> = coefficients
        .iter()
        .map(|c| {
            core::array::from_fn(|f| HarmonicExposureLine {
                observation: &owners[f],
                coefficients: &c[f],
            })
        })
        .collect();
    let factors: Vec<&[HarmonicExposureLine<'_>]> = banks[1..].iter().map(|b| &b[..]).collect();
    PreparedHarmonicExposure::new_with_statistics(
        &owners[0],
        "physical-clock",
        &banks[0],
        &factors,
        statistics,
        limits(),
    )
    .unwrap()
}

fn channels() -> [HarmonicPressureChannel; 3] {
    let tap = |receiver, weight, delay_s| HarmonicPressureTap {
        receiver,
        weight,
        delay_s,
    };
    [
        HarmonicPressureChannel {
            id: "sum".into(),
            taps: vec![tap(0, 0.6, 0.03), tap(1, 0.4, 0.0)],
        },
        HarmonicPressureChannel {
            id: "difference".into(),
            taps: vec![tap(0, 1.0, 0.0), tap(1, -0.7, 0.07)],
        },
        HarmonicPressureChannel {
            id: "repeated-tap".into(),
            taps: vec![tap(0, 0.8, 0.04), tap(0, -0.2, 0.09), tap(1, 0.5, 0.02)],
        },
    ]
}

fn pairs() -> [HarmonicReceiverPair; 4] {
    core::array::from_fn(|i| HarmonicReceiverPair {
        first: i / 2,
        second: i % 2,
    })
}

fn moment_matrix(p: HarmonicExposureMoments) -> [[f64; 4]; 4] {
    core::array::from_fn(|r| {
        core::array::from_fn(|c| match (r, c) {
            (0, 0) => p.mean_square_pressure_pa2,
            (0, c) => p.mean_intensity_w_m2[c - 1],
            (r, 0) => p.mean_intensity_w_m2[r - 1],
            (r, c) => p.velocity_second_m2_s2[r - 1][c - 1],
        })
    })
}

fn covariance_matrix(p: HarmonicReadoutCovariance) -> [[f64; 4]; 4] {
    core::array::from_fn(|r| {
        core::array::from_fn(|c| match (r, c) {
            (0, 0) => p.pressure_covariance_pa2,
            (0, c) => p.first_pressure_second_velocity_covariance_w_m2[c - 1],
            (r, 0) => p.first_velocity_second_pressure_covariance_w_m2[r - 1],
            (r, c) => p.velocity_covariance_m2_s2[r - 1][c - 1],
        })
    })
}

fn observed(exposure: &PreparedHarmonicExposure<'_>, window: HarmonicWindow) -> Vec<f64> {
    let mut values = Vec::new();
    let instantaneous = exposure.observe_window(window).unwrap();
    for &point in instantaneous.points() {
        values.push(point.mean_pressure_pa);
        values.extend(point.mean_particle_velocity_m_s);
        values.extend(moment_matrix(point).into_iter().flatten());
    }
    let readout = exposure
        .observe_readout(window, &pairs(), readout_limits())
        .unwrap();
    for mean in readout.means() {
        values.push(mean.pressure_pa);
        values.extend(mean.particle_velocity_m_s);
    }
    for &covariance in readout.covariances() {
        values.extend(covariance_matrix(covariance).into_iter().flatten());
    }
    let array = PreparedHarmonicPressureArray::new(exposure, &channels(), array_limits()).unwrap();
    assert_eq!(
        array.factor_statistics().collect::<Vec<_>>(),
        exposure.factor_statistics().collect::<Vec<_>>()
    );
    let output = array.observe_window(window).unwrap();
    assert_eq!(output.physical_receivers(), &exposure.reference.receivers);
    assert_eq!(output.sources(), &exposure.reference.sources);
    assert_eq!(output.producer_hash(), exposure.reference.producer_hash);
    assert_eq!(output.source_clock(), "physical-clock");
    for point in output.moments() {
        values.extend([point.mean_pressure_pa, point.mean_square_pressure_pa2]);
    }
    values.extend(output.averaged_covariance_pa2());
    values
}

fn close(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (i, (&a, &b)) in actual.iter().zip(expected).enumerate() {
        assert!(
            (a - b).abs() <= tolerance,
            "entry {i}: {a:e} != {b:e}, tolerance {tolerance:e}"
        );
    }
}

// Eight exact finite states give the declared real/proper joint moments.
// Independently integrate actual delayed real signals, then center the
// per-state averages. There are no ensemble or time samples in production.
fn direct_time(fixture: &Fixture, window: HarmonicWindow) -> Vec<f64> {
    let channels = channels();
    let steps = 2048;
    let mut averages = [[0.0; 11]; 8];
    let mut second = [[[0.0; 4]; 4]; 2];
    let mut array_second = [0.0; 3];
    for (state, average) in averages.iter_mut().enumerate() {
        for sample in 0..=steps {
            let u = sample as f64 / steps as f64;
            let t = window.start_s() + window.duration_s() * u;
            let simpson = if sample == 0 || sample == steps {
                1.0
            } else if sample % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let sensitivity = match window.profile() {
                HarmonicWindowProfile::Uniform => 1.0,
                HarmonicWindowProfile::Triangular => 2.0 * (1.0 - (2.0 * u - 1.0).abs()),
            };
            let weight = simpson * sensitivity / (3.0 * steps as f64);
            let fields: [f64; 8] = core::array::from_fn(|i| fixture.signal(state, i / 4, i % 4, t));
            for i in 0..8 {
                average[i] += weight * fields[i];
            }
            for receiver in 0..2 {
                for row in 0..4 {
                    for col in 0..4 {
                        second[receiver][row][col] +=
                            weight * fields[4 * receiver + row] * fields[4 * receiver + col] / 8.0;
                    }
                }
            }
            for (c, channel) in channels.iter().enumerate() {
                let y: f64 = channel
                    .taps
                    .iter()
                    .map(|tap| tap.weight * fixture.signal(state, tap.receiver, 0, t - tap.delay_s))
                    .sum();
                average[8 + c] += weight * y;
                array_second[c] += weight * y * y / 8.0;
            }
        }
    }
    let mean: [f64; 11] = core::array::from_fn(|c| averages.iter().map(|a| a[c] / 8.0).sum());
    let covariance = |a: usize, b: usize| {
        averages
            .iter()
            .map(|x| (x[a] - mean[a]) * (x[b] - mean[b]) / 8.0)
            .sum::<f64>()
    };
    let mut values = Vec::new();
    for receiver in 0..2 {
        values.extend_from_slice(&mean[4 * receiver..4 * receiver + 4]);
        values.extend(second[receiver].into_iter().flatten());
    }
    values.extend_from_slice(&mean[..8]);
    for pair in pairs() {
        for row in 0..4 {
            for col in 0..4 {
                values.push(covariance(4 * pair.first + row, 4 * pair.second + col));
            }
        }
    }
    for c in 0..3 {
        values.extend([mean[8 + c], array_second[c]]);
    }
    for row in 0..3 {
        for col in 0..3 {
            values.push(covariance(8 + row, 8 + col));
        }
    }
    values
}

#[test]
fn real_factor_retains_subperiod_pressure_velocity_and_pseudocovariance() {
    let mut reference = owner(2.0, 0);
    reference.transfers.fill([C64::ZERO; 4]);
    reference.transfers[0] = [C64::ONE, C64::new(0.0, 1.0), C64::ZERO, C64::ZERO];
    let coefficients = [C64::ONE, C64::ZERO];
    let line = [HarmonicExposureLine {
        observation: &reference,
        coefficients: &coefficients,
    }];
    let real = PreparedHarmonicExposure::new_with_statistics(
        &reference,
        "physical-clock",
        &[],
        &[&line],
        &[Law::Real],
        limits(),
    )
    .unwrap();
    let proper =
        PreparedHarmonicExposure::new(&reference, "physical-clock", &[], &[&line], limits())
            .unwrap();
    let window = HarmonicWindow::new(0.0, 0.5).unwrap();
    let point = real.observe_window(window).unwrap().points()[0];
    close(
        &[point.mean_pressure_pa, point.mean_particle_velocity_m_s[0]],
        &[0.0, 0.0],
        0.0,
    );
    close(
        &[
            point.mean_square_pressure_pa2,
            point.velocity_second_m2_s2[0][0],
            point.mean_intensity_w_m2[0],
        ],
        &[
            0.5 + 2.0_f64.sin() / 4.0,
            0.5 - 2.0_f64.sin() / 4.0,
            (1.0 - 2.0_f64.cos()) / 4.0,
        ],
        3e-15,
    );
    let circular = proper.observe_window(window).unwrap().points()[0];
    close(
        &[
            circular.mean_square_pressure_pa2,
            circular.velocity_second_m2_s2[0][0],
            circular.mean_intensity_w_m2[0],
        ],
        &[0.5, 0.5, 0.0],
        0.0,
    );
    let selected = [HarmonicReceiverPair {
        first: 0,
        second: 0,
    }];
    let covariance = real
        .observe_readout(window, &selected, readout_limits())
        .unwrap()
        .covariances()[0];
    close(
        &[
            covariance.pressure_covariance_pa2,
            covariance.first_pressure_second_velocity_covariance_w_m2[0],
        ],
        &[1.0_f64.sin().powi(2), 1.0_f64.sin() * (1.0 - 1.0_f64.cos())],
        3e-15,
    );
    let covariance = proper
        .observe_readout(window, &selected, readout_limits())
        .unwrap()
        .covariances()[0];
    close(
        &[
            covariance.pressure_covariance_pa2,
            covariance.first_pressure_second_velocity_covariance_w_m2[0],
        ],
        &[1.0 - 1.0_f64.cos(), 0.0],
        3e-15,
    );
}

#[test]
fn mixed_factors_match_finite_states_and_direct_delayed_time_signals() {
    let fixture = Fixture::new();
    let exposure = prepare(
        &fixture.owners,
        &fixture.coefficients,
        &[Law::Real, Law::ProperComplex],
    );
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.13, 0.37, profile).unwrap();
        close(
            &observed(&exposure, window),
            &direct_time(&fixture, window),
            2e-10,
        );
    }
    // The adapter borrows the physical owner, not temporary exposure/statistics.
    let array = {
        let temporary = prepare(
            &fixture.owners,
            &fixture.coefficients,
            &[Law::Real, Law::ProperComplex],
        );
        PreparedHarmonicPressureArray::new(&temporary, &channels(), array_limits()).unwrap()
    };
    assert_eq!(
        array.factor_statistics().collect::<Vec<_>>(),
        [Law::Real, Law::ProperComplex]
    );
    assert!(
        array
            .observe_window(HarmonicWindow::new(0.1, 0.2).unwrap())
            .is_ok()
    );
}

#[test]
fn proper_factor_equals_two_real_quadratures_for_all_observers() {
    let fixture = Fixture::new();
    let original = prepare(
        &fixture.owners,
        &fixture.coefficients,
        &[Law::Real, Law::ProperComplex],
    );
    let scale = core::f64::consts::FRAC_1_SQRT_2;
    let mut expanded = [
        fixture.coefficients[0],
        fixture.coefficients[1],
        fixture.coefficients[2],
        fixture.coefficients[2],
    ];
    for f in 0..2 {
        for s in 0..2 {
            expanded[2][f][s] = fixture.coefficients[2][f][s] * C64::new(scale, 0.0);
            expanded[3][f][s] = fixture.coefficients[2][f][s] * C64::new(0.0, scale);
        }
    }
    let quadratures = prepare(&fixture.owners, &expanded, &[Law::Real; 3]);
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(-0.17, 0.29, profile).unwrap();
        close(
            &observed(&original, window),
            &observed(&quadratures, window),
            8e-14,
        );
    }
}

#[test]
fn real_basis_rotation_and_mixed_factor_permutation_preserve_observables() {
    let fixture = Fixture::new();
    let original = prepare(&fixture.owners, &fixture.coefficients, &[Law::Real; 2]);
    let mut rotated = fixture.coefficients;
    for f in 0..2 {
        for s in 0..2 {
            rotated[1][f][s] = fixture.coefficients[1][f][s] * C64::new(0.6, 0.0)
                + fixture.coefficients[2][f][s] * C64::new(0.8, 0.0);
            rotated[2][f][s] = fixture.coefficients[1][f][s] * C64::new(-0.8, 0.0)
                + fixture.coefficients[2][f][s] * C64::new(0.6, 0.0);
        }
    }
    let mixed = prepare(
        &fixture.owners,
        &fixture.coefficients,
        &[Law::Real, Law::ProperComplex],
    );
    let permuted = [
        fixture.coefficients[0],
        fixture.coefficients[2],
        fixture.coefficients[1],
    ];
    let swapped = prepare(&fixture.owners, &permuted, &[Law::ProperComplex, Law::Real]);
    let rotated = prepare(&fixture.owners, &rotated, &[Law::Real; 2]);
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.23, 0.41, profile).unwrap();
        close(
            &observed(&original, window),
            &observed(&rotated, window),
            1e-13,
        );
        close(
            &observed(&mixed, window),
            &observed(&swapped, window),
            1e-14,
        );
    }
}

#[test]
fn clock_epoch_shift_preserves_nonproper_terms_and_delayed_readouts() {
    let fixture = Fixture::new();
    let original = prepare(
        &fixture.owners,
        &fixture.coefficients,
        &[Law::Real, Law::ProperComplex],
    );
    let shift = 1.125;
    let mut translated = fixture.coefficients;
    for bank in &mut translated {
        for f in 0..2 {
            let (sin, cos) = (fixture.owners[f].omega * shift).sin_cos();
            for coefficient in &mut bank[f] {
                *coefficient = *coefficient * C64::new(cos, sin);
            }
        }
    }
    let translated = prepare(
        &fixture.owners,
        &translated,
        &[Law::Real, Law::ProperComplex],
    );
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let original_window = HarmonicWindow::with_profile(0.125, 0.375, profile).unwrap();
        let shifted_window = HarmonicWindow::with_profile(0.125 + shift, 0.375, profile).unwrap();
        close(
            &observed(&original, original_window),
            &observed(&translated, shifted_window),
            1e-13,
        );
    }
}

#[test]
fn statistics_admission_cancellation_and_legacy_proper_meaning_are_preserved() {
    let fixture = Fixture::new();
    let banks: [[HarmonicExposureLine<'_>; 2]; 3] = core::array::from_fn(|b| {
        core::array::from_fn(|f| HarmonicExposureLine {
            observation: &fixture.owners[f],
            coefficients: &fixture.coefficients[b][f],
        })
    });
    let factors: [&[HarmonicExposureLine<'_>]; 2] = [&banks[1], &banks[2]];
    let explicit = prepare(
        &fixture.owners,
        &fixture.coefficients,
        &[Law::ProperComplex; 2],
    );
    let legacy = PreparedHarmonicExposure::new(
        &fixture.owners[0],
        "physical-clock",
        &banks[0],
        &factors,
        limits(),
    )
    .unwrap();
    let window = HarmonicWindow::with_profile(0.1, 0.4, HarmonicWindowProfile::Triangular).unwrap();
    close(
        &observed(&legacy, window),
        &observed(&explicit, window),
        0.0,
    );
    for laws in [&[][..], &[Law::Real][..], &[Law::Real; 3][..]] {
        assert!(matches!(
            PreparedHarmonicExposure::new_with_statistics(
                &fixture.owners[0],
                "physical-clock",
                &banks[0],
                &factors,
                laws,
                limits(),
            ),
            Err(TransferError::Invalid(_))
        ));
    }
    let mut tight = limits();
    tight.maximum_workspace_bytes = explicit.workspace_bytes() - 1;
    assert!(matches!(
        PreparedHarmonicExposure::new_with_statistics(
            &fixture.owners[0],
            "physical-clock",
            &banks[0],
            &factors,
            &[Law::Real; 2],
            tight,
        ),
        Err(TransferError::Cap { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        PreparedHarmonicExposure::new_with_statistics_cancellable(
            &fixture.owners[0],
            "physical-clock",
            &banks[0],
            &factors,
            &[Law::Real; 2],
            limits(),
            || {
                polls += 1;
                polls >= 12
            },
        ),
        Err(TransferError::Cancelled)
    ));
    let mixed = prepare(
        &fixture.owners,
        &fixture.coefficients,
        &[Law::Real, Law::ProperComplex],
    );
    let mut polls = 0;
    assert!(matches!(
        mixed.observe_window_cancellable(window, || {
            polls += 1;
            polls >= 5
        }),
        Err(TransferError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        mixed.observe_readout_cancellable(window, &pairs(), readout_limits(), || {
            polls += 1;
            polls >= 5
        }),
        Err(TransferError::Cancelled)
    ));
    let array = PreparedHarmonicPressureArray::new(&mixed, &channels(), array_limits()).unwrap();
    let mut polls = 0;
    assert!(matches!(
        array.observe_window_cancellable(window, || {
            polls += 1;
            polls >= 5
        }),
        Err(TransferError::Cancelled)
    ));
    close(
        &observed(&mixed, window),
        &observed(
            &prepare(
                &fixture.owners,
                &fixture.coefficients,
                &[Law::Real, Law::ProperComplex],
            ),
            window,
        ),
        0.0,
    );
    let dark = PreparedHarmonicExposure::new_with_statistics(
        &fixture.owners[0],
        "physical-clock",
        &[],
        &[&[], &[]],
        &[Law::Real, Law::ProperComplex],
        limits(),
    )
    .unwrap();
    assert_eq!(dark.factor_count(), 2);
    assert_eq!(
        dark.factor_statistics().collect::<Vec<_>>(),
        [Law::Real, Law::ProperComplex]
    );
    assert!(observed(&dark, window).iter().all(|&v| v == 0.0));
}
