//! Physical coefficient sensitivities of forced motion, including dark ports.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicParameterControls, HarmonicPortControls, HarmonicResponseControls,
    HarmonicResponseError, PreparedHarmonicResponse, harmonic_mean_power_derivative_with_cancel,
};
use fs_modal::quadratic_modes::QuadraticParameterDirection;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 79,
                kernel_id: 25,
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
        maximum_relative_residual: 2e-12,
    }
}
fn tangent() -> HarmonicParameterControls {
    HarmonicParameterControls {
        maximum_additional_payload_bytes: 1 << 20,
        maximum_work: 1_000_000,
    }
}
fn ports() -> HarmonicPortControls {
    HarmonicPortControls {
        maximum_ports: 4,
        maximum_projection_terms: 100_000,
        maximum_payload_bytes: 1 << 20,
    }
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-280),
        "{a:?} != {b:?}"
    );
}
fn scalar(a: f64, b: f64, tolerance: f64) {
    close(C64::from_re(a), C64::from_re(b), tolerance);
}
fn pair(a: C64, b: C64) -> f64 {
    (a.conj() * b).re
}

#[test]
fn shared_mean_work_tangent_retains_reactive_and_mutual_cancellation_before_rounding() {
    let large = 1e200;
    let z = C64::new(large, large);
    let phase_slope = C64::new(-large, large);
    // A common phase rotation changes no real work, although each real
    // component product overflows if formed before its signed cancellation.
    assert_eq!(
        harmonic_mean_power_derivative_with_cancel(
            &[z], &[z], &[phase_slope], &[phase_slope], || false,
        ).unwrap(),
        0.0,
    );
    // Retain the complete mutual exchange before rounding individual ports.
    // The first two per-port slopes exceed f64, but their exact cancellation
    // leaves the third original port's finite contribution.
    assert_eq!(
        harmonic_mean_power_derivative_with_cancel(
            &[C64::from_re(large), C64::from_re(large), C64::from_re(2.0)],
            &[C64::ZERO; 3],
            &[C64::ZERO; 3],
            &[C64::from_re(large), C64::from_re(-large), C64::from_re(3.0)],
            || false,
        ).unwrap(),
        3.0,
    );
    for exponent in [-900, 900] {
        let flow = 2_f64.powi(exponent);
        let effort = 3.0 * 2_f64.powi(-exponent);
        let work = harmonic_mean_power_derivative_with_cancel(
            &[C64::from_re(flow)], &[C64::from_re(effort)],
            &[C64::from_re(2.0 * flow)], &[C64::from_re(-effort)], || false,
        ).unwrap();
        assert_eq!(work, 1.5);
    }
    let smallest = f64::from_bits(1);
    assert_eq!(
        harmonic_mean_power_derivative_with_cancel(
            &[C64::ZERO], &[C64::from_re(2.0)],
            &[C64::from_re(smallest)], &[C64::ZERO], || false,
        ).unwrap(),
        smallest,
    );
    assert!(harmonic_mean_power_derivative_with_cancel(
        &[C64::from_re(large)], &[C64::ZERO], &[C64::ZERO], &[C64::from_re(large)], || false,
    ).is_err());
    assert!(harmonic_mean_power_derivative_with_cancel(
        &[C64::ZERO], &[C64::from_re(0.5)], &[C64::from_re(smallest)], &[C64::ZERO], || false,
    ).is_err());
}

#[test]
fn shared_mean_work_tangent_checks_shape_finite_inputs_and_late_cancellation() {
    let z = [C64::ONE];
    assert_eq!(
        harmonic_mean_power_derivative_with_cancel(&[], &[], &[], &[], || false).unwrap(),
        0.0,
    );
    assert!(harmonic_mean_power_derivative_with_cancel(&z, &[], &z, &z, || false).is_err());
    assert!(harmonic_mean_power_derivative_with_cancel(
        &z, &z, &[C64::new(f64::NAN, 0.0)], &z, || false,
    ).is_err());
    let mut polls = 0;
    assert!(matches!(
        harmonic_mean_power_derivative_with_cancel(&z, &z, &z, &z, || {
            polls += 1;
            polls == 3
        }),
        Err(HarmonicResponseError::Cancelled),
    ));
}

#[test]
fn scalar_mass_stiffness_damping_and_source_slopes_match_closed_form() {
    with_cx(&CancelGate::new(), |cx| {
        let force = C64::new(1.4, -0.3);
        let df = C64::new(-0.2, 0.07);
        for w in [0.0, 1.8, 2.0, 2.1, 10.0] {
            let p =
                PreparedHarmonicResponse::try_new(&[3.0], &[0.7], &[12.0], 1, w, controls(), cx)
                    .unwrap();
            let result = p
                .solve_with_parameter_derivative(
                    &[force],
                    QuadraticParameterDirection {
                        mass: &[0.2],
                        damping: &[-0.1],
                        stiffness: &[0.8],
                    },
                    Some(&[df]),
                    tangent(),
                    cx,
                )
                .unwrap();
            let a = C64::new(12.0 - 3.0 * w * w, -0.7 * w);
            let da = C64::new(0.8 - 0.2 * w * w, 0.1 * w);
            let q = force / a;
            let dq = (df - da * q) / a;
            let v = C64::new(0.0, -w) * q;
            let dv = C64::new(0.0, -w) * dq;
            close(result.displacement_derivative[0], dq, 4e-13);
            close(result.velocity_derivative[0], dv, 4e-13);
            let input = 0.5 * (pair(df, v) + pair(force, dv));
            let heat = 0.7 * pair(v, dv) - 0.05 * pair(v, v);
            scalar(result.average_input_power_derivative, input, 5e-12);
            scalar(result.average_viscous_power_derivative, heat, 5e-12);
            assert!(result.power_balance_residual_derivative.abs() < 2e-12 * (1.0 + input.abs()));
            assert!(result.relative_derivative_residual < 1e-12);
            assert_eq!(
                result.response.displacement,
                p.solve(&[force], cx).unwrap().displacement
            );
        }
    });
}

#[test]
fn an_exact_antiresonance_has_a_nonzero_physical_sensitivity() {
    with_cx(&CancelGate::new(), |cx| {
        // A=[[0,-1],[-1,0]] at w=1: first coordinate is exactly dark,
        // while q_2=-1. Changing K_22 gives dq_1=-1, not zero or undefined.
        let p = PreparedHarmonicResponse::try_new(
            &[1.0, 0.0, 0.0, 2.0],
            &[0.0; 4],
            &[1.0, -1.0, -1.0, 2.0],
            2,
            1.0,
            controls(),
            cx,
        )
        .unwrap();
        let r = p
            .solve_with_parameter_derivative(
                &[C64::ONE, C64::ZERO],
                QuadraticParameterDirection {
                    mass: &[0.0; 4],
                    damping: &[0.0; 4],
                    stiffness: &[0.0, 0.0, 0.0, 1.0],
                },
                None,
                tangent(),
                cx,
            )
            .unwrap();
        assert_eq!(r.response.displacement[0], C64::ZERO);
        close(r.displacement_derivative[0], C64::from_re(-1.0), 1e-14);
        assert_eq!(r.displacement_derivative[1], C64::ZERO);
        assert_eq!(r.average_viscous_power_derivative, 0.0);
    });
}

#[test]
fn a_parameter_direction_can_couple_a_diagonal_base_without_creating_gyroscopic_heat() {
    with_cx(&CancelGate::new(), |cx| {
        let p = PreparedHarmonicResponse::try_diagonal_modal(
            &[2.0, 3.0],
            &[0.0, 0.0],
            1.0,
            controls(),
            cx,
        )
        .unwrap();
        let direction = QuadraticParameterDirection {
            mass: &[0.0; 4],
            damping: &[0.0, -2.0, 2.0, 0.0],
            stiffness: &[0.0; 4],
        };
        let r = p
            .solve_with_parameter_derivative(&[C64::ONE, C64::ZERO], direction, None, tangent(), cx)
            .unwrap();
        close(
            r.displacement_derivative[1],
            C64::new(0.0, 1.0 / 12.0),
            1e-14,
        );
        assert_eq!(r.displacement_derivative[0], C64::ZERO);
        assert_eq!(r.average_viscous_power_derivative, 0.0);
        assert_eq!(r.power_balance_residual_derivative, 0.0);
    });
}

#[test]
fn coupled_and_loaded_parameter_tangents_match_independent_perturbations_and_power() {
    with_cx(&CancelGate::new(), |cx| {
        let (m, c, k) = (
            [2.0, 0.3, 0.3, 1.0],
            [0.7, -0.2, 0.4, 0.5],
            [8.0, -1.0, -1.0, 5.0],
        );
        let (dm, dc, dk) = (
            [0.2, -0.03, -0.03, -0.1],
            [0.1, -0.2, 0.3, -0.05],
            [0.6, 0.1, 0.1, -0.2],
        );
        let force = [C64::new(1.2, -0.3), C64::new(-0.4, 0.2)];
        let df = [C64::new(-0.1, 0.04), C64::new(0.03, -0.07)];
        let prepare = |theta: f64| {
            let m = std::array::from_fn::<_, 4, _>(|i| m[i] + theta * dm[i]);
            let c = std::array::from_fn::<_, 4, _>(|i| c[i] + theta * dc[i]);
            let k = std::array::from_fn::<_, 4, _>(|i| k[i] + theta * dk[i]);
            PreparedHarmonicResponse::try_with_port_impedance(
                &m,
                &c,
                &k,
                2,
                &[1.0, 0.4],
                1,
                &[C64::new(0.6, -0.3)],
                1.7,
                controls(),
                ports(),
                cx,
            )
            .unwrap()
        };
        let p = prepare(0.0);
        let r = p
            .solve_with_parameter_derivative(
                &force,
                QuadraticParameterDirection {
                    mass: &dm,
                    damping: &dc,
                    stiffness: &dk,
                },
                Some(&df),
                tangent(),
                cx,
            )
            .unwrap();
        let h = 2e-5;
        let plus = prepare(h)
            .solve(
                &std::array::from_fn::<_, 2, _>(|i| force[i] + df[i].scale(h)),
                cx,
            )
            .unwrap();
        let minus = prepare(-h)
            .solve(
                &std::array::from_fn::<_, 2, _>(|i| force[i] - df[i].scale(h)),
                cx,
            )
            .unwrap();
        for i in 0..2 {
            close(
                r.displacement_derivative[i],
                (plus.displacement[i] - minus.displacement[i]).scale(0.5 / h),
                2e-8,
            );
            close(
                r.velocity_derivative[i],
                (plus.velocity[i] - minus.velocity[i]).scale(0.5 / h),
                2e-8,
            );
        }
        scalar(
            r.average_input_power_derivative,
            (plus.average_input_power_w - minus.average_input_power_w) / (2.0 * h),
            2e-8,
        );
        scalar(
            r.average_viscous_power_derivative,
            (plus.average_viscous_power_w - minus.average_viscous_power_w) / (2.0 * h),
            2e-8,
        );
        let slope = r.port_load.as_ref().unwrap();
        let pa = plus.port_load.unwrap();
        let pb = minus.port_load.unwrap();
        close(
            slope.flow_derivative[0],
            (pa.flow[0] - pb.flow[0]).scale(0.5 / h),
            2e-8,
        );
        scalar(
            slope.average_absorbed_power_derivative,
            (pa.average_absorbed_power_w - pb.average_absorbed_power_w) / (2.0 * h),
            2e-8,
        );
        assert!(r.power_balance_residual_derivative.abs() < 2e-12);
        assert!(r.relative_derivative_residual < 1e-12);
    });
}

#[test]
fn large_lossless_source_phase_rotations_keep_zero_input_and_boundary_work() {
    with_cx(&CancelGate::new(), |cx| {
        // A pure phase rotation leaves every physical cycle average unchanged.
        // Individual differentiated products are +/- O(1e400), while all
        // published motion, load variables and their summed powers are finite.
        let amplitude = 1e200;
        for loaded in [false, true] {
            let p = if loaded {
                PreparedHarmonicResponse::try_with_port_impedance(
                    &[1.0],
                    &[0.0],
                    &[2.0],
                    1,
                    &[1.0],
                    1,
                    &[C64::new(0.0, 1.0)],
                    1.0,
                    controls(),
                    ports(),
                    cx,
                )
            } else {
                PreparedHarmonicResponse::try_new(&[1.0], &[0.0], &[2.0], 1, 1.0, controls(), cx)
            }
            .unwrap();
            let result = p
                .solve_with_parameter_derivative(
                    &[C64::from_re(amplitude)],
                    QuadraticParameterDirection {
                        mass: &[0.0],
                        damping: &[0.0],
                        stiffness: &[0.0],
                    },
                    Some(&[C64::new(0.0, amplitude)]),
                    tangent(),
                    cx,
                )
                .unwrap();
            let displacement = if loaded { 0.5 * amplitude } else { amplitude };
            close(
                result.displacement_derivative[0],
                C64::new(0.0, displacement),
                1e-14,
            );
            close(
                result.velocity_derivative[0],
                C64::from_re(displacement),
                1e-14,
            );
            assert_eq!(result.response.average_input_power_w, 0.0);
            assert_eq!(result.average_input_power_derivative, 0.0);
            assert_eq!(result.average_viscous_power_derivative, 0.0);
            assert_eq!(result.power_balance_residual_derivative, 0.0);
            if loaded {
                let port = result.port_load.unwrap();
                close(port.flow_derivative[0], C64::from_re(displacement), 1e-14);
                close(
                    port.effort_derivative[0],
                    C64::new(0.0, displacement),
                    1e-14,
                );
                assert_eq!(port.average_absorbed_power_derivative, 0.0);
            }
        }
    });
}

#[test]
fn large_motion_in_a_damping_direction_nullspace_produces_no_heat_slope() {
    with_cx(&CancelGate::new(), |cx| {
        let p = PreparedHarmonicResponse::try_new(
            &[1.0, 0.0, 0.0, 1.0],
            &[0.0; 4],
            &[2.0, 0.0, 0.0, 2.0],
            2,
            1.0,
            controls(),
            cx,
        )
        .unwrap();
        // The signed motion lies in the exact nullspace of this PSD direction.
        // Entrywise squared powers overflow, but dC*v and v^H*dC*v are zero.
        let result = p
            .solve_with_parameter_derivative(
                &[C64::from_re(1e200), C64::from_re(-1e200)],
                QuadraticParameterDirection {
                    mass: &[0.0; 4],
                    damping: &[1.0; 4],
                    stiffness: &[0.0; 4],
                },
                None,
                tangent(),
                cx,
            )
            .unwrap();
        assert_eq!(result.displacement_derivative, [C64::ZERO; 2]);
        assert_eq!(result.velocity_derivative, [C64::ZERO; 2]);
        assert_eq!(result.average_input_power_derivative, 0.0);
        assert_eq!(result.average_viscous_power_derivative, 0.0);
        assert_eq!(result.power_balance_residual_derivative, 0.0);
    });
}

#[test]
fn physical_coordinate_units_preserve_the_full_parameter_direction() {
    with_cx(&CancelGate::new(), |cx| {
        let run = |s: [f64; 2]| {
            // Both transposed entries use identical arithmetic so the exact
            // symmetry contract survives this change of physical units.
            let scaled =
                |a: [f64; 4]| std::array::from_fn::<_, 4, _>(|i| a[i] / (s[i / 2] * s[i % 2]));
            let p = PreparedHarmonicResponse::try_new(
                &scaled([2.0, 0.2, 0.2, 1.0]),
                &scaled([0.6, 0.1, 0.1, 0.4]),
                &scaled([8.0, -1.0, -1.0, 5.0]),
                2,
                1.7,
                controls(),
                cx,
            )
            .unwrap();
            let dm = scaled([0.1, 0.02, 0.02, -0.1]);
            let dc = scaled([-0.03, 0.01, 0.01, 0.02]);
            let dk = scaled([0.2, -0.1, -0.1, 0.4]);
            p.solve_with_parameter_derivative(
                &[
                    C64::new(1.0, 0.3).scale(1.0 / s[0]),
                    C64::new(-0.2, 0.1).scale(1.0 / s[1]),
                ],
                QuadraticParameterDirection {
                    mass: &dm,
                    damping: &dc,
                    stiffness: &dk,
                },
                None,
                tangent(),
                cx,
            )
            .unwrap()
        };
        let reference = run([1.0, 1.0]);
        for s in [[1e-80, 1e80], [1e80, 1e-80]] {
            let r = run(s);
            for i in 0..2 {
                close(
                    r.displacement_derivative[i].scale(1.0 / s[i]),
                    reference.displacement_derivative[i],
                    2e-12,
                );
            }
            scalar(
                r.average_input_power_derivative,
                reference.average_input_power_derivative,
                2e-11,
            );
            scalar(
                r.average_viscous_power_derivative,
                reference.average_viscous_power_derivative,
                2e-11,
            );
        }
    });
}

#[test]
fn malformed_directions_resource_exhaustion_and_cancellation_preserve_the_factor() {
    with_cx(&CancelGate::new(), |cx| {
        let p = PreparedHarmonicResponse::try_diagonal_modal(&[2.0], &[0.2], 1.0, controls(), cx)
            .unwrap();
        let direction = QuadraticParameterDirection {
            mass: &[0.1],
            damping: &[-0.1],
            stiffness: &[0.4],
        };
        let original = p.solve(&[C64::ONE], cx).unwrap();
        for control in [
            HarmonicParameterControls {
                maximum_work: 0,
                ..tangent()
            },
            HarmonicParameterControls {
                maximum_additional_payload_bytes: 0,
                ..tangent()
            },
        ] {
            assert!(
                p.solve_with_parameter_derivative(&[C64::ONE], direction, None, control, cx)
                    .is_err()
            );
        }
        for bad in [&[][..], &[f64::NAN][..], &[f64::INFINITY][..]] {
            assert!(matches!(
                p.solve_with_parameter_derivative(
                    &[C64::ONE],
                    QuadraticParameterDirection {
                        mass: bad,
                        ..direction
                    },
                    None,
                    tangent(),
                    cx
                ),
                Err(HarmonicResponseError::Invalid(_))
            ));
        }
        assert!(
            p.solve_with_parameter_derivative(&[C64::ONE], direction, Some(&[]), tangent(), cx)
                .is_err()
        );
        let gate = CancelGate::new();
        gate.request();
        with_cx(&gate, |cancelled| {
            assert!(matches!(
                p.solve_with_parameter_derivative(
                    &[C64::ONE],
                    direction,
                    None,
                    tangent(),
                    cancelled
                ),
                Err(HarmonicResponseError::Cancelled)
            ))
        });
        assert_eq!(
            p.solve(&[C64::ONE], cx).unwrap().displacement,
            original.displacement
        );
    });
}
