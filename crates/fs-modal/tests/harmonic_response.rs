//! G1/G3 forced resonance: independent closed forms, reciprocal power and units.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicResponse,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 91,
                kernel_id: 21,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 16,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 1e-12,
    }
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-280),
        "{actual:?} != {expected:?}"
    );
}
fn scalar_close(actual: f64, expected: f64, tolerance: f64) {
    close(C64::from_re(actual), C64::from_re(expected), tolerance);
}

#[test]
fn g1_single_coordinate_resonance_phase_static_limit_and_power() {
    with_cx(&CancelGate::new(), |cx| {
        // m=5 kg, k=20 N/m -> omega_n=2 rad/s. At resonance the
        // velocity is F/c, in phase with the force, including its complex phase.
        let force = C64::new(7.0, -3.0);
        for frequency in [0.0, 0.1, 1.8, 2.0, 2.1, 10.0, 100.0] {
            let prepared = PreparedHarmonicResponse::try_new(
                &[5.0],
                &[2.0],
                &[20.0],
                1,
                frequency,
                controls(),
                cx,
            )
            .unwrap();
            let response = prepared.solve(&[force], cx).unwrap();
            let expected = force / C64::new(20.0 - 5.0 * frequency * frequency, -2.0 * frequency);
            close(response.displacement[0], expected, 3e-14);
            close(
                response.velocity[0],
                C64::new(expected.im, -expected.re).scale(frequency),
                3e-14,
            );
            let heat = response.velocity[0].norm_sq();
            scalar_close(response.average_viscous_power_w, heat, 3e-14);
            assert!(response.power_balance_residual_w.abs() < 2e-13 * heat.max(1e-30));
            assert!(response.relative_residual < 2e-15);
            if frequency == 2.0 {
                close(response.velocity[0], force.scale(0.5), 2e-15);
            }
        }
    });
}

#[test]
fn g1_high_q_modal_resonance_keeps_declared_frequency_and_linewidth() {
    with_cx(&CancelGate::new(), |cx| {
        let natural = 0.1_f64; // squaring then subtracting loses this exact resonance
        let damping = 2e-13;
        for offset in [-2e-13, 0.0, 2e-13] {
            let frequency = natural + offset;
            let response = PreparedHarmonicResponse::try_modal(
                &[natural],
                &[damping],
                frequency,
                controls(),
                cx,
            )
            .unwrap()
            .solve(&[C64::ONE], cx)
            .unwrap();
            let expected = C64::ONE
                / C64::new(
                    (natural - frequency) * (natural + frequency),
                    -frequency * damping,
                );
            close(response.displacement[0], expected, 4e-15);
            if offset == 0.0 {
                assert_eq!(response.displacement[0].re, 0.0);
                scalar_close(response.velocity[0].re, 1.0 / damping, 3e-15);
            }
            assert!(
                response.power_balance_residual_w.abs() <= 1e-14 * response.average_input_power_w
            );
        }
    });
}

#[test]
fn g1_clustered_modes_need_full_nonproportional_damping() {
    with_cx(&CancelGate::new(), |cx| {
        let frequencies = [100.0_f64, 100.1];
        let damping = [4.0, 3.8, 3.8, 4.0];
        let drive = 100.05;
        let force = [C64::ONE, C64::from_re(-1.0)];
        let full =
            PreparedHarmonicResponse::try_modal(&frequencies, &damping, drive, controls(), cx)
                .unwrap()
                .solve(&force, cx)
                .unwrap();
        // Independent inverse of the symmetric 2x2 dynamic stiffness.
        let a = C64::new(
            (frequencies[0] - drive) * (frequencies[0] + drive),
            -4.0 * drive,
        );
        let d = C64::new(
            (frequencies[1] - drive) * (frequencies[1] + drive),
            -4.0 * drive,
        );
        let b = C64::new(0.0, -3.8 * drive);
        let determinant = a * d - b * b;
        close(full.displacement[0], (d + b) / determinant, 4e-14);
        close(full.displacement[1], (-a - b) / determinant, 4e-14);
        let diagonal = PreparedHarmonicResponse::try_modal(
            &frequencies,
            &[4.0, 0.0, 0.0, 4.0],
            drive,
            controls(),
            cx,
        )
        .unwrap()
        .solve(&force, cx)
        .unwrap();
        assert!(full.displacement[0].abs() > 15.0 * diagonal.displacement[0].abs());
        assert!(full.average_input_power_w > 0.0);
        assert!(full.power_balance_residual_w.abs() < 5e-14 * full.average_input_power_w);
    });
}

#[test]
fn g1_antiresonance_is_port_cancellation_not_absent_modal_motion() {
    with_cx(&CancelGate::new(), |cx| {
        // Equal drive/observation projections. At omega=5, the two
        // receptances are -1/24 and +1/24, despite nonzero mode amplitudes.
        let response =
            PreparedHarmonicResponse::try_modal(&[1.0, 7.0], &[0.0; 4], 5.0, controls(), cx)
                .unwrap()
                .solve(&[C64::ONE; 2], cx)
                .unwrap();
        close(response.displacement[0], C64::from_re(-1.0 / 24.0), 1e-15);
        close(response.displacement[1], C64::from_re(1.0 / 24.0), 1e-15);
        assert_eq!(
            response.displacement[0] + response.displacement[1],
            C64::ZERO
        );
        assert_eq!(response.velocity[0] + response.velocity[1], C64::ZERO);
        assert_eq!(response.average_input_power_w, 0.0);
        assert_eq!(response.average_viscous_power_w, 0.0);
    });
}

#[test]
fn g3_full_mass_stiffness_response_is_reciprocal_and_coordinate_covariant() {
    with_cx(&CancelGate::new(), |cx| {
        let m = [2.0, 0.3, 0.3, 1.0];
        let c = [0.7, -0.2, -0.2, 0.3];
        let k = [13.0, -2.0, -2.0, 7.0];
        let frequency = 2.4;
        let prepared =
            PreparedHarmonicResponse::try_new(&m, &c, &k, 2, frequency, controls(), cx).unwrap();
        let first = prepared.solve(&[C64::ONE, C64::ZERO], cx).unwrap();
        let second = prepared.solve(&[C64::ZERO, C64::ONE], cx).unwrap();
        close(first.displacement[1], second.displacement[0], 3e-15);
        let force = [C64::new(1.0, 2.0), C64::new(-0.3, 0.4)];
        let reference = prepared.solve(&force, cx).unwrap();
        let a = C64::new(k[0] - frequency * frequency * m[0], -frequency * c[0]);
        let b = C64::new(k[1] - frequency * frequency * m[1], -frequency * c[1]);
        let d = C64::new(k[3] - frequency * frequency * m[3], -frequency * c[3]);
        let determinant = a * d - b * b;
        close(
            reference.displacement[0],
            (d * force[0] - b * force[1]) / determinant,
            4e-14,
        );
        close(
            reference.displacement[1],
            (a * force[1] - b * force[0]) / determinant,
            4e-14,
        );
        for scales in [[1e-50, 1e50], [1e4, 1e-3]] {
            let transform = |matrix: &[f64]| {
                // Form symmetric pairs once, as physical assemblers do;
                // reversing two rounded scale multiplications can differ by an ulp.
                let mut transformed = vec![0.0; 4];
                for i in 0..2 {
                    for j in 0..=i {
                        let value = matrix[i * 2 + j] * scales[i] * scales[j];
                        transformed[i * 2 + j] = value;
                        transformed[j * 2 + i] = value;
                    }
                }
                transformed
            };
            let scaled = PreparedHarmonicResponse::try_new(
                &transform(&m),
                &transform(&c),
                &transform(&k),
                2,
                frequency,
                controls(),
                cx,
            )
            .unwrap();
            let result = scaled
                .solve(&[force[0].scale(scales[0]), force[1].scale(scales[1])], cx)
                .unwrap();
            for i in 0..2 {
                close(
                    result.displacement[i].scale(scales[i]),
                    reference.displacement[i],
                    5e-14,
                );
                close(
                    result.velocity[i].scale(scales[i]),
                    reference.velocity[i],
                    5e-14,
                );
            }
            scalar_close(
                result.average_input_power_w,
                reference.average_input_power_w,
                5e-14,
            );
            scalar_close(
                result.average_viscous_power_w,
                reference.average_viscous_power_w,
                5e-14,
            );
        }
    });
}

#[test]
fn g1_free_coordinate_keeps_inertial_response_and_singular_limits_refuse() {
    with_cx(&CancelGate::new(), |cx| {
        let response =
            PreparedHarmonicResponse::try_new(&[3.0], &[0.0], &[0.0], 1, 2.0, controls(), cx)
                .unwrap()
                .solve(&[C64::from_re(6.0)], cx)
                .unwrap();
        assert_eq!(response.displacement[0], C64::from_re(-0.5));
        assert_eq!(response.velocity[0], C64::new(0.0, 1.0));
        for (natural, frequency) in [(0.0, 0.0), (0.1, 0.1), (4.0, 4.0)] {
            assert!(matches!(
                PreparedHarmonicResponse::try_modal(&[natural], &[0.0], frequency, controls(), cx),
                Err(HarmonicResponseError::Singular { .. })
            ));
        }
        // A damped collective mode cannot conceal the orthogonal lossless
        // pole in a degenerate cluster. Rank-one damping is PSD, but singular.
        assert!(matches!(
            PreparedHarmonicResponse::try_modal(&[0.1; 2], &[1.0; 4], 0.1, controls(), cx),
            Err(HarmonicResponseError::Singular { .. })
        ));
    });
}

#[test]
fn g0_indefinite_damping_cannot_hide_behind_positive_diagonal_or_large_modes() {
    with_cx(&CancelGate::new(), |cx| {
        assert!(matches!(
            PreparedHarmonicResponse::try_modal(
                &[1.0, 2.0],
                &[1.0, 2.0, 2.0, 1.0],
                0.5,
                controls(),
                cx
            ),
            Err(HarmonicResponseError::NotPositiveSemidefinite {
                matrix: "damping",
                ..
            })
        ));
        let damping = [1.0, 0.0, 0.0, 0.0, 1e-30, 2e-30, 0.0, 2e-30, 1e-30];
        assert!(matches!(
            PreparedHarmonicResponse::try_modal(&[1.0, 2.0, 3.0], &damping, 0.5, controls(), cx),
            Err(HarmonicResponseError::NotPositiveSemidefinite {
                matrix: "damping",
                ..
            })
        ));
        assert!(matches!(
            PreparedHarmonicResponse::try_new(
                &[1.0, 0.0, 0.0, 1.0],
                &[0.0; 4],
                &[1.0, 2.0, 2.0, 1.0],
                2,
                0.5,
                controls(),
                cx
            ),
            Err(HarmonicResponseError::NotPositiveSemidefinite {
                matrix: "stiffness",
                ..
            })
        ));
        assert!(matches!(
            PreparedHarmonicResponse::try_new(
                &[1.0, 2.0, 2.0, 1.0],
                &[0.0; 4],
                &[1.0, 0.0, 0.0, 1.0],
                2,
                0.5,
                controls(),
                cx
            ),
            Err(HarmonicResponseError::MassNotPositiveDefinite)
        ));
    });
}

#[test]
fn g0_budgets_cancellation_nonfinite_forces_and_replay_are_transactional() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let prepared =
            PreparedHarmonicResponse::try_modal(&[2.0], &[0.3], 1.9, controls(), cx).unwrap();
        let first = prepared.solve(&[C64::ONE], cx).unwrap();
        assert!(matches!(
            prepared.solve(&[C64::new(f64::NAN, 0.0)], cx),
            Err(HarmonicResponseError::Invalid(_))
        ));
        assert!(matches!(
            prepared.solve(&[], cx),
            Err(HarmonicResponseError::Invalid(_))
        ));
        let replay = prepared.solve(&[C64::ONE], cx).unwrap();
        assert_eq!(replay.displacement, first.displacement);
        assert_eq!(replay.velocity, first.velocity);
        assert_eq!(
            replay.average_viscous_power_w.to_bits(),
            first.average_viscous_power_w.to_bits()
        );
        let mut limited = controls();
        limited.maximum_payload_bytes = 1;
        assert!(matches!(
            PreparedHarmonicResponse::try_modal(&[2.0], &[0.3], 1.9, limited, cx),
            Err(HarmonicResponseError::Budget { .. })
        ));
        gate.request();
        assert!(matches!(
            prepared.solve(&[C64::ONE], cx),
            Err(HarmonicResponseError::Cancelled)
        ));
        assert!(matches!(
            PreparedHarmonicResponse::try_modal(&[2.0], &[0.3], 1.9, controls(), cx),
            Err(HarmonicResponseError::Cancelled)
        ));
    });
}
