use super::duffing_kinematic_spectrum_exact;
use crate::{KcStringParams, SosModalStorage, StressChannel, kirchhoff_carrier_string};
use fs_phs::Storage;

#[test]
fn g1_kinematic_seed_satisfies_the_original_modal_force_and_channel_admission() {
    let string = kirchhoff_carrier_string(
        &KcStringParams {
            length: 0.65,
            tension: 70.0,
            lin_density: 0.005,
            ea: 100_000.0,
        },
        3,
    )
    .unwrap();
    let generic = SosModalStorage {
        omegas: vec![2.0, 3.0],
        channels: vec![
            StressChannel {
                coefficient: 2.0,
                coupling: vec![1.2, 0.0, 0.0, -2.0],
            },
            StressChannel {
                coefficient: 0.4,
                coupling: vec![-1.0, 0.0, 0.0, 3.0],
            },
        ],
    };
    for (storage, amplitude) in [(&string, 0.001), (&generic, 0.7)] {
        for mode in 0..storage.omegas.len() {
            let spectrum = storage
                .invariant_axis_kinematic_spectrum(mode, amplitude, 24)
                .unwrap();
            let beta: f64 = storage
                .channels
                .iter()
                .map(|channel| {
                    let diagonal = channel.coupling[mode * storage.omegas.len() + mode];
                    channel.coefficient * diagonal * diagonal
                })
                .sum();
            assert_eq!(
                spectrum,
                duffing_kinematic_spectrum_exact(storage.omegas[mode], beta, amplitude, 24)
                    .unwrap()
            );
            assert_eq!(
                spectrum.displacement,
                storage
                    .invariant_axis_spectrum(mode, amplitude, 24)
                    .unwrap()
            );
            for phase in [0.0, 0.3, 1.1, 2.0] {
                let mut state = vec![0.0; 2 * storage.omegas.len()];
                let mut acceleration = 0.0;
                for j in 0..24 {
                    let (sine, cosine) = ((2 * j + 1) as f64 * phase).sin_cos();
                    state[2 * mode] += spectrum.displacement.cosine_amplitudes[j] * cosine;
                    state[2 * mode + 1] += spectrum.velocity_sine_amplitudes[j] * sine;
                    acceleration += spectrum.acceleration_cosine_amplitudes[j] * cosine;
                }
                let mut effort = vec![0.0; state.len()];
                storage.gradient(&state, &mut effort);
                let scale = amplitude * spectrum.displacement.angular_frequency.powi(2);
                assert!((effort[2 * mode] + acceleration).abs() / scale < 3e-12);
                assert_eq!(effort[2 * mode + 1], state[2 * mode + 1]);
                for other in 0..storage.omegas.len() {
                    if other != mode {
                        assert_eq!(effort[2 * other], 0.0);
                    }
                }
            }
        }
    }
    let coupled = SosModalStorage {
        omegas: vec![2.0, 3.0],
        channels: vec![StressChannel {
            coefficient: 1.0,
            coupling: vec![1.0, 1e-30, 1e-30, 1.0],
        }],
    };
    assert!(
        coupled
            .invariant_axis_kinematic_spectrum(0, 0.7, 8)
            .is_err()
    );
    assert!(
        generic
            .invariant_axis_kinematic_spectrum(2, 0.7, 8)
            .is_err()
    );
}
