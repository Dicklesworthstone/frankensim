use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn owner(omega: f64, pressure: &[C64]) -> PreparedHarmonicObservation {
    PreparedHarmonicObservation {
        omega,
        medium: crate::helmholtz::Medium {
            density: 1.2,
            sound_speed: 343.0,
        },
        sources: vec![SourceId::StructuralMode("original-source".into())],
        receivers: pressure
            .iter()
            .enumerate()
            .map(|(r, _)| ReceiverSpec {
                id: format!("physical-{r}"),
                position_m: [0.4 * r as f64, 0.0, 0.0],
            })
            .collect(),
        producer_hash: ContentHash([27; 32]),
        transfers: pressure
            .iter()
            .map(|&p| [p, C64::ZERO, C64::ZERO, C64::ZERO])
            .collect(),
    }
}

fn preparation() -> HarmonicExposureLimits {
    HarmonicExposureLimits {
        maximum_lines: 32,
        maximum_factors: 8,
        maximum_transfer_products: 8192,
        maximum_receiver_frequency_pairs: 8192,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn limits() -> HarmonicPressureArrayLimits {
    HarmonicPressureArrayLimits {
        maximum_channels: 16,
        maximum_taps: 64,
        maximum_prepare_products: 8192,
        maximum_query_products: 8192,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn tap(receiver: usize, weight: f64, delay_s: f64) -> HarmonicPressureTap {
    HarmonicPressureTap {
        receiver,
        weight,
        delay_s,
    }
}

fn channel(id: &str, taps: &[HarmonicPressureTap]) -> HarmonicPressureChannel {
    HarmonicPressureChannel {
        id: id.into(),
        taps: taps.to_vec(),
    }
}

fn line<'a>(
    reference: &'a PreparedHarmonicObservation,
    coefficients: &'a [C64],
) -> HarmonicExposureLine<'a> {
    HarmonicExposureLine {
        observation: reference,
        coefficients,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

struct Fixture {
    owners: [PreparedHarmonicObservation; 2],
    amplitudes: [[C64; 2]; 3],
}

impl Fixture {
    fn new() -> Self {
        Self {
            owners: [
                owner(7.0, &[C64::new(1.1, 0.3), C64::new(-0.2, 0.9)]),
                owner(11.0, &[C64::new(-0.3, 0.2), C64::new(0.8, -0.4)]),
            ],
            amplitudes: [
                [C64::new(0.7, 0.2), C64::new(-0.4, 0.5)],
                [C64::new(0.3, -0.1), C64::new(0.6, 0.8)],
                [C64::new(-0.2, 0.9), C64::new(0.4, -0.3)],
            ],
        }
    }
    fn exposure(&self) -> PreparedHarmonicExposure<'_> {
        let banks: [[HarmonicExposureLine<'_>; 2]; 3] = core::array::from_fn(|b| {
            [
                line(&self.owners[0], &self.amplitudes[b][..1]),
                line(&self.owners[1], &self.amplitudes[b][1..]),
            ]
        });
        PreparedHarmonicExposure::new(
            &self.owners[0],
            "emission-clock",
            &banks[0],
            &[&banks[1], &banks[2]],
            preparation(),
        )
        .unwrap()
    }
}

fn channels() -> Vec<HarmonicPressureChannel> {
    vec![
        channel("forward", &[tap(0, 0.6, 0.03), tap(1, 0.4, 0.0)]),
        channel("difference", &[tap(0, 1.0, 0.0), tap(1, -0.7, 0.07)]),
        channel(
            "multi-tap",
            &[tap(0, 0.8, 0.04), tap(0, -0.2, 0.09), tap(1, 0.5, 0.02)],
        ),
    ]
}

// Independent delayed real-signal integration for every exact proper-phase
// realization, with the output averaged BEFORE its readout covariance.
fn oracle(
    fixture: &Fixture,
    channels: &[HarmonicPressureChannel],
    window: HarmonicWindow,
) -> (Vec<HarmonicPressureChannelMoments>, Vec<f64>) {
    let phases = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)];
    let n = channels.len();
    let steps = 4096;
    let mut realizations = vec![vec![0.0; n]; 16];
    let mut moments = vec![HarmonicPressureChannelMoments::default(); n];
    for (state, averages) in realizations.iter_mut().enumerate() {
        for sample in 0..=steps {
            let u = sample as f64 / steps as f64;
            let time = window.start_s() + window.duration_s() * u;
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
            for (c, channel) in channels.iter().enumerate() {
                let mut value = 0.0;
                for tap in &channel.taps {
                    for f in 0..2 {
                        let amplitude = fixture.owners[f].transfers[tap.receiver][0]
                            * (fixture.amplitudes[0][f]
                                + phases[state % 4] * fixture.amplitudes[1][f]
                                + phases[state / 4] * fixture.amplitudes[2][f]);
                        let angle = fixture.owners[f].omega * (time - tap.delay_s);
                        let (sin, cos) = angle.sin_cos();
                        value += tap.weight * (amplitude.re * cos + amplitude.im * sin);
                    }
                }
                averages[c] += weight * value;
                moments[c].mean_square_pressure_pa2 += weight * value * value / 16.0;
            }
        }
    }
    for c in 0..n {
        moments[c].mean_pressure_pa = realizations.iter().map(|row| row[c]).sum::<f64>() / 16.0;
    }
    let covariance = (0..n * n)
        .map(|index| {
            let row = index / n;
            let col = index % n;
            realizations
                .iter()
                .map(|v| {
                    (v[row] - moments[row].mean_pressure_pa)
                        * (v[col] - moments[col].mean_pressure_pa)
                })
                .sum::<f64>()
                / 16.0
        })
        .collect();
    (moments, covariance)
}

#[test]
fn delayed_multifrequency_channels_match_direct_time_and_proper_ensemble() {
    let fixture = Fixture::new();
    let exposure = fixture.exposure();
    let channels = channels();
    let array = PreparedHarmonicPressureArray::new(&exposure, &channels, limits()).unwrap();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.43, profile).unwrap();
        let observed = array.observe_window(window).unwrap();
        let (expected, covariance) = oracle(&fixture, &channels, window);
        for (actual, expected) in observed.moments().iter().zip(expected) {
            close(actual.mean_pressure_pa, expected.mean_pressure_pa, 3e-11);
            close(
                actual.mean_square_pressure_pa2,
                expected.mean_square_pressure_pa2,
                3e-11,
            );
        }
        for (&actual, expected) in observed.averaged_covariance_pa2().iter().zip(covariance) {
            close(actual, expected, 3e-11);
        }
        for row in 0..channels.len() {
            for col in 0..channels.len() {
                assert_eq!(
                    observed.averaged_covariance_pa2()[row * channels.len() + col],
                    observed.averaged_covariance_pa2()[col * channels.len() + row]
                );
            }
        }
        assert_eq!(observed.physical_receivers(), exposure.reference.receivers);
        assert_eq!(observed.channels(), channels);
        assert_eq!(observed.producer_hash(), exposure.reference.producer_hash);
        assert_eq!(observed.source_clock(), "emission-clock");
    }
}

#[test]
fn causal_steering_aligns_a_propagating_wave_at_both_frequencies_and_places_a_null() {
    let travel_s = 0.17;
    let frequencies = [3.0_f64, 11.0];
    let owners = frequencies.map(|omega| {
        let angle = omega * travel_s;
        owner(omega, &[C64::ONE, C64::new(angle.cos(), angle.sin())])
    });
    let amplitudes = [[C64::new(1.0, 0.3)], [C64::new(0.6, -0.2)]];
    let mean = [
        line(&owners[0], &amplitudes[0]),
        line(&owners[1], &amplitudes[1]),
    ];
    let exposure =
        PreparedHarmonicExposure::new(&owners[0], "plane-wave-clock", &mean, &[], preparation())
            .unwrap();
    let channels = [
        channel("aligned", &[tap(0, 0.5, travel_s), tap(1, 0.5, 0.0)]),
        channel("null", &[tap(0, 0.5, travel_s), tap(1, -0.5, 0.0)]),
        channel("unsteered", &[tap(0, 0.5, 0.0), tap(1, 0.5, 0.0)]),
    ];
    let prepared = PreparedHarmonicPressureArray::new(&exposure, &channels, limits()).unwrap();
    let observed = prepared
        .observe_window(HarmonicWindow::new(0.0, core::f64::consts::TAU).unwrap())
        .unwrap();
    let total = 0.5 * (amplitudes[0][0].norm_sq() + amplitudes[1][0].norm_sq());
    close(observed.moments()[0].mean_square_pressure_pa2, total, 2e-13);
    assert!(observed.moments()[1].mean_square_pressure_pa2 < 1e-28);
    assert!(observed.moments()[2].mean_square_pressure_pa2 < total - 0.1);
    assert!(observed.averaged_covariance_pa2().iter().all(|&v| v == 0.0));
}

#[test]
fn identity_channel_matches_original_observers_and_full_period_readout_cancels() {
    let reference = owner(core::f64::consts::TAU, &[C64::new(0.7, -0.3)]);
    let mean_amplitude = [C64::new(0.4, 0.8)];
    let noise = [C64::new(-0.2, 0.7)];
    let mean = [line(&reference, &mean_amplitude)];
    let factors = [line(&reference, &noise)];
    let exposure =
        PreparedHarmonicExposure::new(&reference, "source", &mean, &[&factors], preparation())
            .unwrap();
    let channels = [channel("identity", &[tap(0, 1.0, 0.0)])];
    let array = PreparedHarmonicPressureArray::new(&exposure, &channels, limits()).unwrap();
    for window in [
        HarmonicWindow::new(0.0, 1.0).unwrap(),
        HarmonicWindow::with_profile(0.03, 0.21, HarmonicWindowProfile::Triangular).unwrap(),
    ] {
        let observed = array.observe_window(window).unwrap();
        let point = exposure.observe_window(window).unwrap();
        let readout = exposure
            .observe_readout(
                window,
                &[HarmonicReceiverPair {
                    first: 0,
                    second: 0,
                }],
                HarmonicReadoutLimits {
                    maximum_receiver_pairs: 4,
                    maximum_products: 4096,
                    maximum_workspace_bytes: 1_000_000,
                },
            )
            .unwrap();
        close(
            observed.moments()[0].mean_pressure_pa,
            point.points()[0].mean_pressure_pa,
            1e-15,
        );
        close(
            observed.moments()[0].mean_square_pressure_pa2,
            point.points()[0].mean_square_pressure_pa2,
            1e-15,
        );
        close(
            observed.averaged_covariance_pa2()[0],
            readout.covariances()[0].pressure_covariance_pa2,
            1e-15,
        );
    }
    let long = array
        .observe_window(HarmonicWindow::new(0.0, 1.0).unwrap())
        .unwrap();
    assert!(long.averaged_covariance_pa2()[0] < 1e-28);
    assert!(long.moments()[0].mean_square_pressure_pa2 > 0.2);
}

#[test]
fn delayed_taps_keep_weak_residuals_and_sub_ulp_cosine_decrements() {
    let reference = owner(1.0, &[C64::ONE]);
    let unit = [C64::ONE];
    let line = [line(&reference, &unit)];
    let exposure =
        PreparedHarmonicExposure::new(&reference, "source", &line, &[], preparation()).unwrap();
    let delay = 2.0_f64.powi(-40);
    let channels = [
        channel("tiny difference", &[tap(0, 1.0, 0.0), tap(0, -1.0, delay)]),
        channel(
            "weak survivor",
            &[
                tap(0, 1e16, 0.031),
                tap(0, 1.0, 0.031),
                tap(0, -1e16, 0.031),
            ],
        ),
        channel("unit", &[tap(0, 1.0, 0.031)]),
    ];
    let array = PreparedHarmonicPressureArray::new(&exposure, &channels, limits()).unwrap();
    let result = array
        .observe_window(HarmonicWindow::new(-0.5, 1.0).unwrap())
        .unwrap();
    let expected = 0.5 * delay * delay * (0.5_f64.sin() / 0.5);
    assert!(result.moments()[0].mean_pressure_pa > 0.0);
    close(
        result.moments()[0].mean_pressure_pa,
        expected,
        expected * 3e-14,
    );
    close(
        result.moments()[1].mean_pressure_pa,
        result.moments()[2].mean_pressure_pa,
        3e-15,
    );
    close(
        result.moments()[1].mean_square_pressure_pa2,
        result.moments()[2].mean_square_pressure_pa2,
        3e-15,
    );
}

#[test]
fn common_latency_and_channel_permutation_preserve_readout_physics() {
    let fixture = Fixture::new();
    let exposure = fixture.exposure();
    let channels = channels();
    let array = PreparedHarmonicPressureArray::new(&exposure, &channels, limits()).unwrap();
    let window =
        HarmonicWindow::with_profile(-0.13, 0.71, HarmonicWindowProfile::Triangular).unwrap();
    let result = array.observe_window(window).unwrap();
    let latency = 0.37;
    let order = [2, 0, 1];
    let translated: Vec<_> = order
        .iter()
        .map(|&c| {
            let mut channel = channels[c].clone();
            for tap in &mut channel.taps {
                tap.delay_s += latency;
            }
            channel
        })
        .collect();
    let shifted = PreparedHarmonicPressureArray::new(&exposure, &translated, limits()).unwrap();
    let shifted = shifted
        .observe_window(
            HarmonicWindow::with_profile(
                window.start_s() + latency,
                window.duration_s(),
                window.profile(),
            )
            .unwrap(),
        )
        .unwrap();
    for c in 0..3 {
        close(
            shifted.moments()[c].mean_pressure_pa,
            result.moments()[order[c]].mean_pressure_pa,
            3e-14,
        );
        close(
            shifted.moments()[c].mean_square_pressure_pa2,
            result.moments()[order[c]].mean_square_pressure_pa2,
            3e-14,
        );
        for d in 0..3 {
            close(
                shifted.averaged_covariance_pa2()[c * 3 + d],
                result.averaged_covariance_pa2()[order[c] * 3 + order[d]],
                3e-14,
            );
        }
    }
}

#[test]
fn bounds_identity_cancellation_and_explicit_dark_channels_are_atomic() {
    let fixture = Fixture::new();
    let exposure = fixture.exposure();
    let channels = channels();
    for tap in [
        tap(2, 1.0, 0.0),
        tap(0, f64::NAN, 0.0),
        tap(0, 1.0, -0.1),
        tap(0, 1.0, f64::INFINITY),
        tap(0, 1.0, f64::MAX),
    ] {
        assert!(
            PreparedHarmonicPressureArray::new(&exposure, &[channel("bad", &[tap])], limits())
                .is_err()
        );
    }
    assert!(
        PreparedHarmonicPressureArray::new(
            &exposure,
            &[channel("duplicate", &[]), channel("duplicate", &[]),],
            limits()
        )
        .is_err()
    );
    assert!(
        PreparedHarmonicPressureArray::new(&exposure, &[channel(&"x".repeat(129), &[])], limits())
            .is_err()
    );
    for choice in 0..5 {
        let mut budget = limits();
        match choice {
            0 => budget.maximum_channels = 0,
            1 => budget.maximum_taps = 0,
            2 => budget.maximum_prepare_products = 0,
            3 => budget.maximum_query_products = 0,
            _ => budget.maximum_workspace_bytes = 0,
        }
        assert!(matches!(
            PreparedHarmonicPressureArray::new(&exposure, &channels, budget),
            Err(TransferError::Cap { .. })
        ));
    }
    let mut checks = 0;
    assert!(matches!(
        PreparedHarmonicPressureArray::new_cancellable(&exposure, &channels, limits(), || {
            checks += 1;
            checks == 18
        }),
        Err(TransferError::Cancelled)
    ));
    let array = PreparedHarmonicPressureArray::new(&exposure, &channels, limits()).unwrap();
    let window = HarmonicWindow::new(0.0, 0.1).unwrap();
    checks = 0;
    assert!(matches!(
        array.observe_window_cancellable(window, || {
            checks += 1;
            checks == 7
        }),
        Err(TransferError::Cancelled)
    ));
    let a = array.observe_window(window).unwrap();
    let b = array.observe_window(window).unwrap();
    assert_eq!(a.moments(), b.moments());
    assert_eq!(a.averaged_covariance_pa2(), b.averaged_covariance_pa2());
    let dark =
        PreparedHarmonicPressureArray::new(&exposure, &[channel("dark", &[])], limits()).unwrap();
    let dark = dark.observe_window(window).unwrap();
    assert_eq!(dark.moments(), &[HarmonicPressureChannelMoments::default()]);
    assert_eq!(dark.averaged_covariance_pa2(), &[0.0]);
    let empty = PreparedHarmonicPressureArray::new(&exposure, &[], limits()).unwrap();
    assert!(empty.observe_window(window).unwrap().moments().is_empty());
}
