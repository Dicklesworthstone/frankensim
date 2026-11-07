use super::*;
use crate::harmonic::HarmonicResponse;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_math::c64::C64;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 97,
                kernel_id: 13,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn controls() -> HarmonicReceiverControls {
    HarmonicReceiverControls {
        maximum_receivers: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_projection_terms: 10_000,
        maximum_product_terms: 10_000,
    }
}

fn row(weights: &[f64], quantity: HarmonicReceiverQuantity) -> HarmonicReceiver<'_> {
    HarmonicReceiver { weights, quantity }
}

// Analytically specified particular motion. These tests exercise the receiver
// of an already solved response, not the dynamic-stiffness solver. The fixture
// is realizable with M=I, C=0, K=2I and F=(2-omega^2)q away from its poles.
fn tone(omega: f64, displacement: &[C64]) -> Tone {
    let velocity: Vec<_> = displacement
        .iter()
        .map(|q| C64::new(omega * q.im, -omega * q.re))
        .collect();
    Tone {
        omega,
        force: displacement
            .iter()
            .map(|q| q.scale(2.0 - omega * omega))
            .collect(),
        mass_velocity: velocity.iter().map(|&v| parts(v)).collect(),
        stiffness_displacement: displacement.iter().map(|q| parts(q.scale(2.0))).collect(),
        viscous_velocity: vec![[Dd::ZERO; 2]; displacement.len()],
        response: HarmonicResponse {
            displacement: displacement.to_vec(),
            velocity,
            relative_residual: 0.0,
            average_input_power_w: 0.0,
            average_viscous_power_w: 0.0,
            port_load: None,
            power_balance_residual_w: 0.0,
        },
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * (1.0 + b.abs()), "{a} != {b}");
}

#[test]
fn physical_rows_cancel_before_mean_square_and_acceleration_differentiation() {
    with_cx(&CancelGate::new(), |cx| {
        let source = PreparedHarmonicInterference {
            n: 3,
            tones: vec![tone(
                0.0,
                &[C64::from_re(1e16), C64::from_re(1.0), C64::from_re(-1e16)],
            )],
            maximum_relative_response_residual: 3e-15,
        };
        let rows = [
            row(&[1.0, 1.0, 1.0], HarmonicReceiverQuantity::Displacement),
            row(&[-1.0, -1.0, -1.0], HarmonicReceiverQuantity::Displacement),
            row(&[1.0, 1.0, 1.0], HarmonicReceiverQuantity::Acceleration),
        ];
        let receiver = source.prepare_receivers(&rows, controls(), cx).unwrap();
        let window = HarmonicWindow::new(0.2, 0.01).unwrap();
        let own = receiver.observe_pair(0, 0, window, cx).unwrap();
        assert_eq!(own.mean_left, 1.0);
        assert_eq!(own.mean_product, 1.0);
        assert_eq!(own.mean_stochastic_product, 0.0);
        assert_eq!(own.averaged_covariance, 0.0);
        let cross = receiver.observe_pair(0, 1, window, cx).unwrap();
        assert_eq!(cross.mean_right, -1.0);
        assert_eq!(cross.mean_product, -1.0);
        assert_eq!(
            receiver
                .observe_pair(2, 2, window, cx)
                .unwrap()
                .mean_product,
            0.0
        );
        assert_eq!(receiver.receiver_count(), 3);
        assert_eq!(receiver.maximum_relative_response_residual(), 3e-15);
    });
}

#[test]
fn quiet_variance_survives_large_coherent_bias_and_proper_dc_quadrature() {
    with_cx(&CancelGate::new(), |cx| {
        let mean = [tone(0.0, &[C64::from_re(1e14)])];
        let real = [tone(0.0, &[C64::new(0.25, 99.0)])];
        let proper = [tone(0.0, &[C64::new(0.0, 2.0)])];
        let banks = [
            (HarmonicFactorStatistics::Real, mean.as_slice()),
            (HarmonicFactorStatistics::Real, real.as_slice()),
            (HarmonicFactorStatistics::ProperComplex, proper.as_slice()),
        ];
        let receiver = PreparedHarmonicReceivers::prepare(
            1,
            &[row(&[1.0], HarmonicReceiverQuantity::Displacement)],
            banks.into_iter(),
            0.0,
            controls(),
            cx,
        )
        .unwrap();
        let value = receiver
            .observe_pair(0, 0, HarmonicWindow::new(3.0, 2.0).unwrap(), cx)
            .unwrap();
        assert_eq!(value.mean_left, 1e14);
        assert_eq!(value.mean_stochastic_product, 2.0625);
        assert_eq!(value.averaged_covariance, 2.0625);
        assert_eq!(value.mean_coherent_product, 1e28);
    });
}

#[test]
fn averaging_a_carrier_does_not_erase_instantaneous_motion_variance() {
    with_cx(&CancelGate::new(), |cx| {
        let proper = [tone(8.0, &[C64::new(1.0, 2.0)])];
        let banks = [
            (HarmonicFactorStatistics::Real, &[][..]),
            (HarmonicFactorStatistics::ProperComplex, proper.as_slice()),
        ];
        let receiver = PreparedHarmonicReceivers::prepare(
            1,
            &[row(&[1.0], HarmonicReceiverQuantity::Displacement)],
            banks.into_iter(),
            0.0,
            controls(),
            cx,
        )
        .unwrap();
        let window = HarmonicWindow::new(0.17, core::f64::consts::TAU / 8.0).unwrap();
        let value = receiver.observe_pair(0, 0, window, cx).unwrap();
        close(value.mean_stochastic_product, 2.5, 2e-15);
        assert!(value.averaged_covariance.abs() < 1e-28);
        assert_eq!(value.mean_left, 0.0);
    });
}

fn instantaneous(
    tones: &[Tone],
    receiver: HarmonicReceiver<'_>,
    multiplier: C64,
    time: f64,
) -> f64 {
    let mut result = 0.0;
    for tone in tones {
        let phase = C64::new((tone.omega * time).cos(), -(tone.omega * time).sin());
        for (i, &weight) in receiver.weights.iter().enumerate() {
            let amplitude = match receiver.quantity {
                HarmonicReceiverQuantity::Displacement => tone.response.displacement[i],
                HarmonicReceiverQuantity::Velocity => tone.response.velocity[i],
                HarmonicReceiverQuantity::Acceleration => {
                    let v = tone.response.velocity[i];
                    C64::new(tone.omega * v.im, -tone.omega * v.re)
                }
            };
            result += weight * (amplitude * multiplier * phase).re;
        }
    }
    result
}

#[test]
fn coherent_beats_and_mixed_factor_cross_moments_match_real_trajectory_ensemble() {
    use fs_fft::harmonic_window::HarmonicWindowProfile;
    with_cx(&CancelGate::new(), |cx| {
        let mean = [
            tone(4.0, &[C64::new(0.7, 0.2), C64::new(-0.3, 0.5)]),
            tone(4.7, &[C64::new(-0.4, 0.1), C64::new(0.6, -0.3)]),
        ];
        let real = [
            tone(4.2, &[C64::new(0.3, -0.5), C64::new(-0.2, 0.6)]),
            tone(5.1, &[C64::new(-0.2, 0.4), C64::new(0.8, 0.1)]),
        ];
        let proper = [
            tone(3.9, &[C64::new(0.8, 0.3), C64::new(0.4, 0.2)]),
            tone(4.5, &[C64::new(-0.3, 0.2), C64::new(0.1, -0.7)]),
        ];
        let rows = [
            row(&[0.7, -0.3], HarmonicReceiverQuantity::Displacement),
            row(&[0.4, 0.2], HarmonicReceiverQuantity::Acceleration),
            row(&[0.3, -0.8], HarmonicReceiverQuantity::Velocity),
        ];
        let banks = [
            (HarmonicFactorStatistics::Real, mean.as_slice()),
            (HarmonicFactorStatistics::Real, real.as_slice()),
            (HarmonicFactorStatistics::ProperComplex, proper.as_slice()),
        ];
        let receiver =
            PreparedHarmonicReceivers::prepare(2, &rows, banks.into_iter(), 0.0, controls(), cx)
                .unwrap();
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let window = HarmonicWindow::with_profile(0.37, 1.1, profile).unwrap();
            for (left, right) in [(0, 1), (1, 0), (2, 2)] {
                let actual = receiver.observe_pair(left, right, window, cx).unwrap();
                // Independent midpoint quadrature of eight equally likely real
                // trajectories: xi=+/-1, z=+/-1,+/-i. They have exactly the
                // declared real/proper joint second moments.
                let mut samples = Vec::new();
                let mut coherent_product = 0.0;
                for real_weight in [-1.0, 1.0] {
                    for proper_weight in [
                        C64::ONE,
                        C64::from_re(-1.0),
                        C64::new(0.0, 1.0),
                        C64::new(0.0, -1.0),
                    ] {
                        let mut average = [0.0; 3];
                        for step in 0..8_000 {
                            let u = (step as f64 + 0.5) / 8_000.0;
                            let time = window.start_s() + window.duration_s() * u;
                            let weight = match profile {
                                HarmonicWindowProfile::Uniform => 1.0,
                                HarmonicWindowProfile::Triangular => 4.0 * u.min(1.0 - u),
                            } / 8_000.0;
                            let pair = [left, right].map(|i| {
                                instantaneous(&mean, rows[i], C64::ONE, time)
                                    + instantaneous(&real, rows[i], C64::from_re(real_weight), time)
                                    + instantaneous(&proper, rows[i], proper_weight, time)
                            });
                            average[0] += weight * pair[0];
                            average[1] += weight * pair[1];
                            average[2] += weight * pair[0] * pair[1];
                            coherent_product += weight / 8.0
                                * instantaneous(&mean, rows[left], C64::ONE, time)
                                * instantaneous(&mean, rows[right], C64::ONE, time);
                        }
                        samples.push(average);
                    }
                }
                let expected: [f64; 3] =
                    core::array::from_fn(|i| samples.iter().map(|s| s[i]).sum::<f64>() / 8.0);
                let covariance = samples
                    .iter()
                    .map(|s| (s[0] - expected[0]) * (s[1] - expected[1]))
                    .sum::<f64>()
                    / 8.0;
                close(actual.mean_left, expected[0], 3e-7);
                close(actual.mean_right, expected[1], 3e-7);
                close(actual.mean_product, expected[2], 3e-7);
                close(actual.mean_coherent_product, coherent_product, 3e-7);
                close(
                    actual.mean_stochastic_product,
                    expected[2] - coherent_product,
                    3e-7,
                );
                close(actual.averaged_covariance, covariance, 3e-7);
            }
        }
    });
}

#[test]
fn explicit_limits_invalid_rows_indices_and_cancellation_refuse() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let source = PreparedHarmonicInterference {
            n: 1,
            tones: vec![tone(2.0, &[C64::ONE])],
            maximum_relative_response_residual: 0.0,
        };
        let rows = [row(&[1.0], HarmonicReceiverQuantity::Velocity)];
        for controls in [
            HarmonicReceiverControls {
                maximum_receivers: 0,
                ..controls()
            },
            HarmonicReceiverControls {
                maximum_payload_bytes: 0,
                ..controls()
            },
            HarmonicReceiverControls {
                maximum_projection_terms: 0,
                ..controls()
            },
            HarmonicReceiverControls {
                maximum_product_terms: 0,
                ..controls()
            },
        ] {
            assert!(source.prepare_receivers(&rows, controls, cx).is_err());
        }
        for invalid in [&[][..], &[f64::NAN][..], &[1.0, 2.0][..]] {
            assert!(
                source
                    .prepare_receivers(
                        &[row(invalid, HarmonicReceiverQuantity::Displacement)],
                        controls(),
                        cx,
                    )
                    .is_err()
            );
        }
        assert!(source.prepare_receivers(&[], controls(), cx).is_err());
        let receiver = source.prepare_receivers(&rows, controls(), cx).unwrap();
        let window = HarmonicWindow::new(0.0, 0.5).unwrap();
        assert!(receiver.observe_pair(1, 0, window, cx).is_err());
        gate.request();
        assert!(matches!(
            receiver.observe_pair(0, 0, window, cx),
            Err(HarmonicResponseError::Cancelled)
        ));
        assert!(matches!(
            source.prepare_receivers(&rows, controls(), cx),
            Err(HarmonicResponseError::Cancelled)
        ));
    });
}

#[test]
fn original_coupled_force_owners_feed_receivers_and_reconstruct_physical_energy() {
    use crate::harmonic::{
        HarmonicEnsembleControls, HarmonicForceFactor, HarmonicForceTone,
        HarmonicInterferenceControls, HarmonicResponseControls,
        PreparedHarmonicInterferenceEnsemble,
    };
    with_cx(&CancelGate::new(), |cx| {
        let mass = [2.0, 0.3, 0.3, 1.0];
        let stiffness = [5.0, -1.0, -1.0, 3.0];
        // The skew part changes the phase/coupling, but contributes no heat.
        let damping = [0.2, 0.1, -0.1, 0.2];
        let force_a = [C64::new(1.0, 0.3), C64::new(-0.4, 0.1)];
        let force_b = [C64::new(-0.3, 0.2), C64::new(0.6, -0.2)];
        let forces = [
            HarmonicForceTone {
                angular_frequency_rad_s: 0.8,
                force: &force_a,
            },
            HarmonicForceTone {
                angular_frequency_rad_s: 1.3,
                force: &force_b,
            },
        ];
        let response_controls = HarmonicResponseControls {
            maximum_coordinates: 2,
            maximum_payload_bytes: 1 << 20,
            maximum_relative_residual: 1e-12,
        };
        let coherent = PreparedHarmonicInterference::try_new(
            &mass,
            &damping,
            &stiffness,
            2,
            &forces,
            response_controls,
            HarmonicInterferenceControls {
                maximum_tones: 4,
                maximum_payload_bytes: 1 << 20,
                maximum_product_terms: 10_000,
            },
            cx,
        )
        .unwrap();
        let factors = [
            HarmonicForceFactor {
                tones: &forces[..1],
                statistics: HarmonicFactorStatistics::Real,
            },
            HarmonicForceFactor {
                tones: &forces,
                statistics: HarmonicFactorStatistics::ProperComplex,
            },
        ];
        let ensemble = PreparedHarmonicInterferenceEnsemble::try_new(
            &mass,
            &damping,
            &stiffness,
            2,
            &forces,
            &factors,
            response_controls,
            HarmonicEnsembleControls {
                maximum_factors: 2,
                maximum_tones: 8,
                maximum_payload_bytes: 1 << 20,
                maximum_product_terms: 10_000,
            },
            cx,
        )
        .unwrap();
        let rows = [
            row(&[1.0, 0.0], HarmonicReceiverQuantity::Velocity),
            row(&[0.0, 1.0], HarmonicReceiverQuantity::Velocity),
            row(&[1.0, 0.0], HarmonicReceiverQuantity::Displacement),
            row(&[0.0, 1.0], HarmonicReceiverQuantity::Displacement),
        ];
        let window = HarmonicWindow::new(0.27, 0.43).unwrap();
        for (receiver, energy) in [
            (
                coherent.prepare_receivers(&rows, controls(), cx).unwrap(),
                coherent.observe(window, cx).unwrap(),
            ),
            (
                ensemble.prepare_receivers(&rows, controls(), cx).unwrap(),
                ensemble.observe(window, cx).unwrap(),
            ),
        ] {
            for (offset, matrix, expected) in [
                (0, mass, energy.mean_kinetic_energy_j),
                (2, stiffness, energy.mean_elastic_energy_j),
            ] {
                let mut recovered = 0.0;
                for i in 0..2 {
                    for j in 0..2 {
                        recovered += 0.5
                            * matrix[2 * i + j]
                            * receiver
                                .observe_pair(offset + i, offset + j, window, cx)
                                .unwrap()
                                .mean_product;
                    }
                }
                close(recovered, expected, 3e-14);
            }
        }
    });
}
