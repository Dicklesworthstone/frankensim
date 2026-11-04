//! Compact loaded modal banks: bare poles, full port physics and original audits.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicFrequencyControls, HarmonicPortControls, HarmonicResponseControls,
    HarmonicResponseError, PreparedHarmonicResponse,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 91,
                kernel_id: 24,
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
        maximum_coordinates: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 1e-12,
    }
}
fn ports() -> HarmonicPortControls {
    HarmonicPortControls {
        maximum_ports: 4,
        maximum_projection_terms: 1_000_000,
        maximum_payload_bytes: 1 << 20,
    }
}
fn tangent_controls() -> HarmonicFrequencyControls {
    HarmonicFrequencyControls {
        maximum_additional_payload_bytes: 1 << 20,
        maximum_work: 1_000_000,
    }
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-250),
        "{actual:?} != {expected:?}"
    );
}
fn power_close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 5e-10 * expected.abs().max(1.0));
}

#[test]
fn large_loaded_bank_lifts_a_bare_pole_and_keeps_the_dc_load_tangent() {
    with_cx(&CancelGate::new(), |cx| {
        let n = 1025;
        let natural: Vec<_> = (0..n)
            .map(|i| if i == 0 { 2.0 } else { 3.0 + i as f64 / 8.0 })
            .collect();
        let rates: Vec<_> = (0..n).map(|i| if i == 0 { 0.0 } else { 0.03 }).collect();
        let map: Vec<_> = (0..n).map(|i| 1.0 / (i + 1) as f64).collect();
        let mut force = vec![C64::ZERO; n];
        force[0] = C64::new(1.0, -0.3);
        let bounded = HarmonicResponseControls {
            maximum_coordinates: n,
            maximum_payload_bytes: 272 * n,
            maximum_relative_residual: 1e-11,
        };
        let load_limits = HarmonicPortControls {
            maximum_ports: 1,
            maximum_projection_terms: 2_000_000,
            maximum_payload_bytes: 128 * n + 8192,
        };
        // Causal dashpot plus inertance: Z = c_load - i*w*m_load.
        let w = 2.0;
        let z = C64::new(0.7, -0.1 * w);
        let prepared = PreparedHarmonicResponse::try_diagonal_modal_with_port_impedance(
            &natural,
            &rates,
            &map,
            1,
            &[z],
            w,
            bounded,
            load_limits,
            cx,
        )
        .unwrap();
        let response = prepared.solve(&force, cx).unwrap();
        // At the undamped first bare pole, -i*w*B_0*Z*(B^T q)=F_0.
        // All remaining modal equations then give their own scalar motion.
        let port_displacement = force[0] / (C64::new(0.0, -w) * z).scale(map[0]);
        let mut expected = vec![C64::ZERO; n];
        let mut other_port_displacement = C64::ZERO;
        for i in 1..n {
            let d = C64::new((natural[i] - w) * (natural[i] + w), -w * rates[i]);
            expected[i] = -force[0].scale(map[i] / map[0]) / d;
            other_port_displacement = other_port_displacement + expected[i].scale(map[i]);
        }
        expected[0] = (port_displacement - other_port_displacement).scale(1.0 / map[0]);
        for (actual, expected) in response.displacement.iter().zip(&expected) {
            close(*actual, *expected, 3e-11);
        }
        let port = response.port_load.as_ref().unwrap();
        close(port.flow[0], force[0] / z.scale(map[0]), 3e-12);
        power_close(
            port.average_absorbed_power_w,
            0.5 * z.re * port.flow[0].norm_sq(),
        );
        assert!(response.relative_residual < bounded.maximum_relative_residual);
        assert!(
            response.power_balance_residual_w.abs()
                < 3e-10 * response.average_input_power_w.abs().max(1.0)
        );
        assert_eq!(
            prepared.solve(&force, cx).unwrap().displacement,
            response.displacement
        );

        // A finite velocity impedance leaves static motion unchanged, while
        // its real drag still enters the first frequency derivative at DC.
        let dc = PreparedHarmonicResponse::try_diagonal_modal_with_port_impedance(
            &natural,
            &rates,
            &map,
            1,
            &[C64::from_re(0.7)],
            0.0,
            bounded,
            load_limits,
            cx,
        )
        .unwrap()
        .solve_with_frequency_derivative(
            &force,
            None,
            Some(&[C64::new(0.0, -0.1)]),
            tangent_controls(),
            cx,
        )
        .unwrap();
        let q0 = force[0].scale(0.25);
        close(dc.response.displacement[0], q0, 3e-14);
        for i in 0..n {
            let expected = C64::new(0.0, 0.7 * map[i] * map[0] / natural[i].powi(2)) * q0;
            close(dc.displacement_derivative[i], expected, 5e-13);
        }
        assert_eq!(
            dc.response
                .port_load
                .as_ref()
                .unwrap()
                .average_absorbed_power_w,
            0.0
        );
    });
}

#[test]
fn full_mutual_rank_deficient_and_lossless_loads_match_dense_motion_and_tangents() {
    with_cx(&CancelGate::new(), |cx| {
        let map = [1.0, 0.3, -0.2, 1.1, 0.4, -0.5];
        let rates = [0.0, 0.0, 0.04];
        let mut full_damping = [0.0; 9];
        full_damping[8] = rates[2];
        let force = [C64::new(0.8, -0.1), C64::new(-0.3, 0.7), C64::new(0.2, 0.4)];
        let slope = [C64::new(0.1, -0.03); 3];
        for variant in 0..3 {
            let natural = [0.1, if variant == 1 { 0.3 } else { 0.1 }, 2.0];
            let (resistance, inertance) = match variant {
                // Passive nonreciprocal drag and full mutual inertance.
                0 => ([0.7, 0.15, -0.05, 0.8], [0.2, 0.03, 0.03, 0.1]),
                // Rank-one Z: the impedance must never be inverted.
                1 => ([0.5, 1.0, 1.0, 2.0], [0.1, 0.2, 0.2, 0.4]),
                _ => ([0.0; 4], [0.2, 0.03, 0.03, 0.1]),
            };
            for w in [0.0, 0.1 - 1e-12, 0.1, 0.1 + 1e-12, 0.7] {
                let z: [C64; 4] =
                    std::array::from_fn(|i| C64::new(resistance[i], -w * inertance[i]));
                let dz: [C64; 4] = std::array::from_fn(|i| C64::new(0.0, -inertance[i]));
                let compact = PreparedHarmonicResponse::try_diagonal_modal_with_port_impedance(
                    &natural,
                    &rates,
                    &map,
                    2,
                    &z,
                    w,
                    controls(),
                    ports(),
                    cx,
                )
                .unwrap()
                .solve_with_frequency_derivative(
                    &force,
                    Some(&slope),
                    Some(&dz),
                    tangent_controls(),
                    cx,
                )
                .unwrap();
                let dense = PreparedHarmonicResponse::try_modal_with_port_impedance(
                    &natural,
                    &full_damping,
                    &map,
                    2,
                    &z,
                    w,
                    controls(),
                    ports(),
                    cx,
                )
                .unwrap()
                .solve_with_frequency_derivative(
                    &force,
                    Some(&slope),
                    Some(&dz),
                    tangent_controls(),
                    cx,
                )
                .unwrap();
                for i in 0..3 {
                    close(
                        compact.response.displacement[i],
                        dense.response.displacement[i],
                        3e-10,
                    );
                    close(
                        compact.response.velocity[i],
                        dense.response.velocity[i],
                        3e-10,
                    );
                    close(
                        compact.displacement_derivative[i],
                        dense.displacement_derivative[i],
                        5e-10,
                    );
                    close(
                        compact.velocity_derivative[i],
                        dense.velocity_derivative[i],
                        5e-10,
                    );
                }
                let cp = compact.response.port_load.as_ref().unwrap();
                let dp = dense.response.port_load.as_ref().unwrap();
                let ct = compact.port_load.as_ref().unwrap();
                let dt = dense.port_load.as_ref().unwrap();
                for a in 0..2 {
                    close(cp.flow[a], dp.flow[a], 3e-10);
                    close(cp.effort[a], dp.effort[a], 3e-10);
                    close(ct.flow_derivative[a], dt.flow_derivative[a], 5e-10);
                    close(ct.effort_derivative[a], dt.effort_derivative[a], 5e-10);
                }
                power_close(cp.average_absorbed_power_w, dp.average_absorbed_power_w);
                power_close(
                    compact.response.average_viscous_power_w,
                    dense.response.average_viscous_power_w,
                );
                power_close(
                    compact.average_viscous_power_derivative,
                    dense.average_viscous_power_derivative,
                );
                power_close(
                    ct.average_absorbed_power_derivative,
                    dt.average_absorbed_power_derivative,
                );
                assert!(compact.response.relative_residual < controls().maximum_relative_residual);
                assert!(
                    compact.relative_derivative_residual < controls().maximum_relative_residual
                );
                power_close(compact.response.power_balance_residual_w, 0.0);
                power_close(compact.power_balance_residual_derivative, 0.0);
            }
        }
    });
}

#[test]
fn an_uncoupled_huge_response_does_not_erase_a_tiny_retained_force() {
    with_cx(&CancelGate::new(), |cx| {
        let natural = [1.0, 1.0];
        let rates = [1e-200, 0.0];
        let map = [0.0, 1e-100];
        let force = [C64::ONE, C64::from_re(1e-200)];
        let compact = PreparedHarmonicResponse::try_diagonal_modal_with_port_impedance(
            &natural,
            &rates,
            &map,
            1,
            &[C64::ONE],
            1.0,
            controls(),
            ports(),
            cx,
        )
        .unwrap()
        .solve(&force, cx)
        .unwrap();
        let dense = PreparedHarmonicResponse::try_modal_with_port_impedance(
            &natural,
            &[1e-200, 0.0, 0.0, 0.0],
            &map,
            1,
            &[C64::ONE],
            1.0,
            controls(),
            ports(),
            cx,
        )
        .unwrap()
        .solve(&force, cx)
        .unwrap();
        close(compact.displacement[0], C64::new(0.0, 1e200), 3e-14);
        close(compact.displacement[1], C64::new(0.0, 1.0), 3e-14);
        for i in 0..2 {
            close(compact.displacement[i], dense.displacement[i], 3e-14);
            close(compact.velocity[i], dense.velocity[i], 3e-14);
        }
        let cp = compact.port_load.as_ref().unwrap();
        let dp = dense.port_load.as_ref().unwrap();
        close(cp.flow[0], C64::from_re(1e-100), 3e-14);
        close(
            C64::from_re(cp.average_absorbed_power_w),
            C64::from_re(dp.average_absorbed_power_w),
            5e-14,
        );
        assert!(compact.relative_residual < 1e-12);
    });
}

#[test]
fn unsupported_poles_core_resources_and_cancellation_refuse() {
    with_cx(&CancelGate::new(), |cx| {
        assert!(matches!(
            PreparedHarmonicResponse::try_diagonal_modal_with_port_impedance(
                &[1.0, 1.0],
                &[0.0, 0.0],
                &[1.0, 1.0],
                1,
                &[C64::ONE],
                1.0,
                controls(),
                ports(),
                cx,
            ),
            Err(HarmonicResponseError::Singular { .. })
        ));
        assert!(matches!(
            PreparedHarmonicResponse::try_diagonal_modal_with_port_impedance(
                &[1.0],
                &[0.0],
                &[0.0],
                1,
                &[C64::ONE],
                1.0,
                controls(),
                ports(),
                cx,
            ),
            Err(HarmonicResponseError::Singular { .. })
        ));
        let small_payload = HarmonicPortControls {
            maximum_payload_bytes: 512, // admits source ports, not the retained core
            ..ports()
        };
        assert!(matches!(
            PreparedHarmonicResponse::try_diagonal_modal_with_port_impedance(
                &[1.0, 1.0000001],
                &[0.0, 0.0],
                &[1.0, 1.0],
                1,
                &[C64::ONE],
                1.0,
                controls(),
                small_payload,
                cx,
            ),
            Err(HarmonicResponseError::Budget {
                quantity: "compact load numerical payload bytes",
                ..
            })
        ));
        let small_work = HarmonicPortControls {
            maximum_projection_terms: 6, // original source projection fits exactly
            ..ports()
        };
        assert!(matches!(
            PreparedHarmonicResponse::try_diagonal_modal_with_port_impedance(
                &[1.0, 1.0000001],
                &[0.0, 0.0],
                &[1.0, 1.0],
                1,
                &[C64::ONE],
                1.0,
                controls(),
                small_work,
                cx,
            ),
            Err(HarmonicResponseError::WorkBudget { cap: 6, .. })
        ));
    });
    let gate = CancelGate::new();
    gate.request();
    with_cx(&gate, |cx| {
        assert!(matches!(
            PreparedHarmonicResponse::try_diagonal_modal_with_port_impedance(
                &[1.0],
                &[0.1],
                &[1.0],
                1,
                &[C64::ONE],
                1.0,
                controls(),
                ports(),
                cx,
            ),
            Err(HarmonicResponseError::Cancelled)
        ));
    });
}
