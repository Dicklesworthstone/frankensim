//! Complete changing-attachment/load jets against independently assembled physics.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicFrequencyParameterControls, HarmonicFrequencyParameterResponse, HarmonicPortControls,
    HarmonicPortFrequencyParameterDirection, HarmonicPortParameterDirection,
    HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicResponse,
};
use fs_modal::quadratic_modes::QuadraticParameterDirection;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 83,
                kernel_id: 257,
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
        maximum_relative_residual: 2e-11,
    }
}
fn port_controls() -> HarmonicPortControls {
    HarmonicPortControls {
        maximum_ports: 4,
        maximum_projection_terms: 100_000,
        maximum_payload_bytes: 1 << 20,
    }
}
fn limits() -> HarmonicFrequencyParameterControls {
    HarmonicFrequencyParameterControls {
        maximum_additional_payload_bytes: 1 << 20,
        maximum_work: 1_000_000,
    }
}
fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() < 8e-10 * expected.abs().max(1e-10),
        "{actual:?} != {expected:?}",
    );
}
fn close_power(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-10 * expected.abs().max(1e-8),
        "{actual:e} != {expected:e}",
    );
}

fn assemble(base: [f64; 4], b: [f64; 4], load: [f64; 4]) -> [f64; 4] {
    let mut result = base;
    for i in 0..2 {
        for j in 0..2 {
            for a in 0..2 {
                for c in 0..2 {
                    result[2 * i + j] += b[2 * i + a] * load[2 * a + c] * b[2 * j + c];
                }
            }
        }
    }
    result
}
// The physical mass/stiffness tensors are symmetric. Mirror their shared
// off-diagonal coefficient after independent assembly, so the reference owner
// does not receive roundoff asymmetry from different summation orders.
fn exactly_symmetric(mut matrix: [f64; 4]) -> [f64; 4] {
    let cross = 0.5 * (matrix[1] + matrix[2]);
    matrix[1] = cross;
    matrix[2] = cross;
    matrix
}

fn assemble_direction(
    base: [f64; 4],
    b: [f64; 4],
    db: [f64; 4],
    load: [f64; 4],
    dload: [f64; 4],
) -> [f64; 4] {
    let mut result = base;
    for i in 0..2 {
        for j in 0..2 {
            for a in 0..2 {
                for c in 0..2 {
                    result[2 * i + j] += db[2 * i + a] * load[2 * a + c] * b[2 * j + c]
                        + b[2 * i + a] * dload[2 * a + c] * b[2 * j + c]
                        + b[2 * i + a] * load[2 * a + c] * db[2 * j + c];
                }
            }
        }
    }
    result
}
fn compare_motion(
    actual: &HarmonicFrequencyParameterResponse,
    expected: &HarmonicFrequencyParameterResponse,
) {
    for i in 0..2 {
        close(
            actual
                .frequency_response
                .frequency_response
                .response
                .displacement[i],
            expected
                .frequency_response
                .frequency_response
                .response
                .displacement[i],
        );
        close(
            actual
                .frequency_response
                .frequency_response
                .displacement_derivative[i],
            expected
                .frequency_response
                .frequency_response
                .displacement_derivative[i],
        );
        close(
            actual.frequency_response.displacement_second_derivative[i],
            expected.frequency_response.displacement_second_derivative[i],
        );
        close(
            actual.parameter_response.displacement_derivative[i],
            expected.parameter_response.displacement_derivative[i],
        );
        close(
            actual.displacement_mixed_derivative[i],
            expected.displacement_mixed_derivative[i],
        );
        close(
            actual.velocity_mixed_derivative[i],
            expected.velocity_mixed_derivative[i],
        );
    }
    assert!(actual.relative_mixed_derivative_residual < 3e-12);
}

#[test]
fn dense_changing_load_matches_original_added_mass_damping_and_stiffness_in_extreme_port_units() {
    with_cx(&CancelGate::new(), |cx| {
        let m = [1.4, 0.1, 0.1, 0.9];
        // Positive symmetric damping plus a retained gyroscopic component.
        let c = [0.2, 0.13, -0.07, 0.25];
        let k = [4.0, -0.6, -0.6, 6.0];
        let dm = [0.12, 0.03, 0.03, -0.06];
        let dc = [0.04, -0.09, 0.07, -0.02];
        let dk = [0.3, -0.11, -0.11, 0.2];
        let b = [0.8, -0.2, 0.3, 1.1];
        let db = [0.1, 0.07, -0.08, 0.04];
        let lm = [0.4, 0.08, 0.08, 0.7];
        // The symmetric part is positive; the off-diagonal skew part must survive.
        let lc = [0.6, 0.4, -0.1, 0.7];
        let lk = [0.5, -0.12, -0.12, 0.9];
        let dlm = [0.05, -0.02, -0.02, -0.03];
        let dlc = [-0.03, 0.08, -0.04, 0.07];
        let dlk = [0.13, 0.06, 0.06, -0.08];
        let me = exactly_symmetric(assemble(m, b, lm));
        let ce = assemble(c, b, lc);
        let ke = exactly_symmetric(assemble(k, b, lk));
        let dme = exactly_symmetric(assemble_direction(dm, b, db, lm, dlm));
        let dce = assemble_direction(dc, b, db, lc, dlc);
        let dke = exactly_symmetric(assemble_direction(dk, b, db, lk, dlk));
        let f = [C64::new(1.7, -0.6), C64::new(-0.2, 0.4)];
        let fw = [C64::new(-0.2, 0.8), C64::new(0.6, -0.2)];
        let fww = [C64::new(0.3, -0.4), C64::new(-0.1, 0.3)];
        let ft = [C64::new(0.1, 0.7), C64::new(0.4, 0.2)];
        let fwt = [C64::new(-0.6, 0.2), C64::new(0.3, -0.1)];

        for w in [0.04, 0.4, 1.4, 2.0, 4.0] {
            // Z=R-i*w*L+i*S/w contributes exactly B*L*B^T to mass,
            // B*R*B^T to damping and B*S*B^T to stiffness.
            let z = core::array::from_fn::<_, 4, _>(|i| C64::new(lc[i], lk[i] / w - w * lm[i]));
            let zw = core::array::from_fn::<_, 4, _>(|i| C64::new(0.0, -lm[i] - lk[i] / (w * w)));
            let zww = core::array::from_fn::<_, 4, _>(|i| C64::new(0.0, 2.0 * lk[i] / (w * w * w)));
            let zt = core::array::from_fn::<_, 4, _>(|i| C64::new(dlc[i], dlk[i] / w - w * dlm[i]));
            let zwt =
                core::array::from_fn::<_, 4, _>(|i| C64::new(0.0, -dlm[i] - dlk[i] / (w * w)));
            let effective =
                PreparedHarmonicResponse::try_new(&me, &ce, &ke, 2, w, controls(), cx).unwrap();
            let expected = effective
                .solve_with_frequency_parameter_derivative(
                    &f,
                    QuadraticParameterDirection {
                        mass: &dme,
                        damping: &dce,
                        stiffness: &dke,
                    },
                    Some(&fw),
                    Some(&fww),
                    Some(&ft),
                    Some(&fwt),
                    limits(),
                    cx,
                )
                .unwrap();

            for exponents in [[0, 0], [400, -400]] {
                let units = exponents.map(|power| 2.0_f64.powi(power));
                let scaled_b = core::array::from_fn::<_, 4, _>(|i| b[i] * units[i % 2]);
                let scaled_db = core::array::from_fn::<_, 4, _>(|i| db[i] * units[i % 2]);
                let transform = |values: [C64; 4]| {
                    core::array::from_fn::<_, 4, _>(|i| {
                        values[i].scale(2.0_f64.powi(-exponents[i / 2] - exponents[i % 2]))
                    })
                };
                let (zs, zws, zwws, zts, zwts) = (
                    transform(z),
                    transform(zw),
                    transform(zww),
                    transform(zt),
                    transform(zwt),
                );
                let owner = PreparedHarmonicResponse::try_with_port_impedance(
                    &m,
                    &c,
                    &k,
                    2,
                    &scaled_b,
                    2,
                    &zs,
                    w,
                    controls(),
                    port_controls(),
                    cx,
                )
                .unwrap();
                let actual = owner
                    .solve_with_frequency_parameter_and_port_derivative(
                        &f,
                        QuadraticParameterDirection {
                            mass: &dm,
                            damping: &dc,
                            stiffness: &dk,
                        },
                        Some(&fw),
                        Some(&fww),
                        Some(&ft),
                        Some(&fwt),
                        &zws,
                        &zwws,
                        HarmonicPortFrequencyParameterDirection {
                            parameter: HarmonicPortParameterDirection {
                                force_map: &scaled_db,
                                impedance: &zts,
                            },
                            impedance_angular_parameter_derivative: &zwts,
                        },
                        limits(),
                        cx,
                    )
                    .unwrap();
                compare_motion(&actual, &expected);
                close_power(
                    actual.average_input_power_mixed_derivative,
                    expected.average_input_power_mixed_derivative,
                );
                close_power(
                    actual.average_viscous_power_mixed_derivative
                        + actual
                            .port_load
                            .as_ref()
                            .unwrap()
                            .average_absorbed_power_mixed_derivative,
                    expected.average_viscous_power_mixed_derivative,
                );
                let work_scale = actual.average_input_power_mixed_derivative.abs()
                    + actual.average_viscous_power_mixed_derivative.abs()
                    + actual
                        .port_load
                        .as_ref()
                        .unwrap()
                        .average_absorbed_power_mixed_derivative
                        .abs();
                assert!(
                    actual.power_balance_residual_mixed_derivative.abs()
                        < 3e-11 * work_scale.max(1e-10)
                );
            }
        }
    });
}

#[test]
fn changing_load_admission_is_complete_bounded_cancellable_and_reusable() {
    with_cx(&CancelGate::new(), |cx| {
        let owner = PreparedHarmonicResponse::try_with_port_impedance(
            &[1.0],
            &[0.1],
            &[1.0],
            1,
            &[1.0],
            1,
            &[C64::from_re(0.2)],
            0.8,
            controls(),
            port_controls(),
            cx,
        )
        .unwrap();
        let zero = [0.0];
        let complex_zero = [C64::ZERO];
        let direction = QuadraticParameterDirection {
            mass: &zero,
            damping: &zero,
            stiffness: &zero,
        };
        let port_direction = HarmonicPortFrequencyParameterDirection {
            parameter: HarmonicPortParameterDirection {
                force_map: &zero,
                impedance: &complex_zero,
            },
            impedance_angular_parameter_derivative: &complex_zero,
        };
        let solve = |port, controls, context: &Cx<'_>| {
            owner.solve_with_frequency_parameter_and_port_derivative(
                &[C64::ONE],
                direction,
                None,
                None,
                None,
                None,
                &complex_zero,
                &complex_zero,
                port,
                controls,
                context,
            )
        };
        let exact = HarmonicFrequencyParameterControls {
            maximum_additional_payload_bytes: 5000,
            maximum_work: 7680,
        };
        let result = solve(port_direction, exact, cx).unwrap();
        assert_eq!(result.displacement_mixed_derivative, vec![C64::ZERO]);
        assert_eq!(
            result.port_load.unwrap().flow_mixed_derivative,
            vec![C64::ZERO]
        );
        assert!(matches!(
            solve(
                port_direction,
                HarmonicFrequencyParameterControls {
                    maximum_additional_payload_bytes: 4999,
                    ..exact
                },
                cx
            ),
            Err(HarmonicResponseError::Budget {
                requested: 5000,
                ..
            })
        ));
        assert!(matches!(
            solve(
                port_direction,
                HarmonicFrequencyParameterControls {
                    maximum_work: 7679,
                    ..exact
                },
                cx
            ),
            Err(HarmonicResponseError::WorkBudget {
                requested: 7680,
                ..
            })
        ));
        let bad_real = [f64::INFINITY];
        let bad_complex = [C64::new(f64::NAN, 0.0)];
        for bad in [
            HarmonicPortFrequencyParameterDirection {
                parameter: HarmonicPortParameterDirection {
                    force_map: &[],
                    ..port_direction.parameter
                },
                ..port_direction
            },
            HarmonicPortFrequencyParameterDirection {
                parameter: HarmonicPortParameterDirection {
                    force_map: &bad_real,
                    ..port_direction.parameter
                },
                ..port_direction
            },
            HarmonicPortFrequencyParameterDirection {
                parameter: HarmonicPortParameterDirection {
                    impedance: &[],
                    ..port_direction.parameter
                },
                ..port_direction
            },
            HarmonicPortFrequencyParameterDirection {
                parameter: HarmonicPortParameterDirection {
                    impedance: &bad_complex,
                    ..port_direction.parameter
                },
                ..port_direction
            },
            HarmonicPortFrequencyParameterDirection {
                impedance_angular_parameter_derivative: &[],
                ..port_direction
            },
            HarmonicPortFrequencyParameterDirection {
                impedance_angular_parameter_derivative: &bad_complex,
                ..port_direction
            },
        ] {
            assert!(matches!(
                solve(bad, exact, cx),
                Err(HarmonicResponseError::Invalid(_))
            ));
        }
        let cancelled = CancelGate::new();
        cancelled.request();
        with_cx(&cancelled, |context| {
            assert!(matches!(
                solve(port_direction, exact, context),
                Err(HarmonicResponseError::Cancelled)
            ));
        });
        let baseline = owner.solve(&[C64::ONE], cx).unwrap();
        let after = solve(port_direction, exact, cx).unwrap();
        assert_eq!(
            after
                .frequency_response
                .frequency_response
                .response
                .displacement,
            baseline.displacement
        );

        let unloaded =
            PreparedHarmonicResponse::try_new(&[1.0], &[0.1], &[1.0], 1, 0.8, controls(), cx)
                .unwrap();
        assert!(matches!(
            unloaded.solve_with_frequency_parameter_and_port_derivative(
                &[C64::ONE],
                direction,
                None,
                None,
                None,
                None,
                &complex_zero,
                &complex_zero,
                port_direction,
                limits(),
                cx,
            ),
            Err(HarmonicResponseError::Invalid(_))
        ));
    });
}
