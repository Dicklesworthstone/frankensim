//! Original-coordinate frequency tangents, causal loads and independent work.

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
                kernel_id: 23,
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
fn tangent_controls() -> HarmonicFrequencyControls {
    HarmonicFrequencyControls {
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
fn scalar_frequency_source_and_dc_kinematics_match_closed_form() {
    with_cx(&CancelGate::new(), |cx| {
        let f = C64::new(2.0, -0.7);
        let df = C64::new(-0.3, 0.8);
        for w in [0.0, 0.1, 1.8, 2.0, 2.3, 10.0, 100.0] {
            let p =
                PreparedHarmonicResponse::try_new(&[5.0], &[2.0], &[20.0], 1, w, controls(), cx)
                    .unwrap();
            let r = p
                .solve_with_frequency_derivative(&[f], Some(&[df]), None, tangent_controls(), cx)
                .unwrap();
            let a = C64::new(20.0 - 5.0 * w * w, -2.0 * w);
            let q = f / a;
            let dq = (df + C64::new(10.0 * w, 2.0) * q) / a;
            let v = C64::new(0.0, -w) * q;
            let dv = C64::new(0.0, -1.0) * (q + dq.scale(w));
            close(r.response.displacement[0], q, 5e-14);
            close(r.displacement_derivative[0], dq, 7e-14);
            close(r.velocity_derivative[0], dv, 7e-14);
            scalar(
                r.average_input_power_derivative,
                0.5 * (pair(df, v) + pair(f, dv)),
                2e-12,
            );
            scalar(r.average_viscous_power_derivative, 2.0 * pair(v, dv), 2e-12);
            assert!(
                r.power_balance_residual_derivative.abs()
                    < 2e-12 * (1.0 + r.average_viscous_power_derivative.abs())
            );
            assert!(r.relative_derivative_residual < 2e-14);
            if w == 0.0 {
                close(r.velocity_derivative[0], C64::new(0.0, -1.0) * q, 2e-15);
            }
        }
    });
}

#[test]
fn high_q_modal_derivative_preserves_factored_nominal_resonance() {
    with_cx(&CancelGate::new(), |cx| {
        let natural = 0.1;
        let damping = 2e-13;
        for delta in [-2e-13, 0.0, 2e-13] {
            let w = natural + delta;
            let p = PreparedHarmonicResponse::try_modal(&[natural], &[damping], w, controls(), cx)
                .unwrap();
            let r = p
                .solve_with_frequency_derivative(&[C64::ONE], None, None, tangent_controls(), cx)
                .unwrap();
            let a = C64::new((natural - w) * (natural + w), -w * damping);
            let q = C64::ONE / a;
            close(
                r.displacement_derivative[0],
                C64::new(2.0 * w, damping) * q / a,
                8e-15,
            );
            close(
                r.velocity_derivative[0],
                C64::new(0.0, -1.0) * (q + r.displacement_derivative[0].scale(w)),
                8e-15,
            );
            if delta == 0.0 {
                assert_eq!(r.response.displacement[0].re, 0.0);
                scalar((r.displacement_derivative[0] / q).im, 2.0 / damping, 8e-15);
            }
        }
    });
}

#[test]
fn causal_added_mass_and_spring_differentiate_the_loaded_frequency_not_frozen_z() {
    with_cx(&CancelGate::new(), |cx| {
        let f = C64::new(3.0, -0.7);
        let (load_c, load_m, load_k) = (0.4, 0.7, 0.2);
        let m: f64 = 2.0 + 9.0 * load_m;
        let k = 8.0 + 9.0 * load_k;
        let c = 0.1 + 9.0 * load_c;
        for w in [0.01, 0.7, (k / m).sqrt(), 2.0, 4.0, 100.0] {
            let z = C64::new(load_c, load_k / w - w * load_m);
            let dz = C64::new(0.0, -load_k / (w * w) - load_m);
            let p = PreparedHarmonicResponse::try_with_port_impedance(
                &[2.0],
                &[0.1],
                &[8.0],
                1,
                &[3.0],
                1,
                &[z],
                w,
                controls(),
                ports(),
                cx,
            )
            .unwrap();
            let r = p
                .solve_with_frequency_derivative(&[f], None, Some(&[dz]), tangent_controls(), cx)
                .unwrap();
            let a = C64::new(k - m * w * w, -c * w);
            let q = f / a;
            let dq = C64::new(2.0 * m * w, c) * q / a;
            let dv = C64::new(0.0, -1.0) * (q + dq.scale(w));
            close(r.displacement_derivative[0], dq, 2e-13);
            close(r.velocity_derivative[0], dv, 2e-13);
            let primal = r.response.port_load.as_ref().unwrap();
            let derivative = r.port_load.as_ref().unwrap();
            close(derivative.flow_derivative[0], dv.scale(3.0), 2e-13);
            close(
                derivative.effort_derivative[0],
                z * derivative.flow_derivative[0] + dz * primal.flow[0],
                2e-13,
            );
            // At the velocity maximum the signed heat slope vanishes. Bound
            // rounding by the original bilinear terms, not the cancelled slope.
            let heat_slope = load_c * pair(primal.flow[0], derivative.flow_derivative[0]);
            let heat_scale = load_c * primal.flow[0].abs() * derivative.flow_derivative[0].abs();
            assert!(
                (derivative.average_absorbed_power_derivative - heat_slope).abs()
                    <= 4e-13 * heat_scale.max(1e-280)
            );
            assert!(
                r.power_balance_residual_derivative.abs()
                    < 2e-11 * (1.0 + r.average_input_power_derivative.abs())
            );
            let frozen = p
                .solve_with_frequency_derivative(
                    &[f],
                    None,
                    Some(&[C64::ZERO]),
                    tangent_controls(),
                    cx,
                )
                .unwrap();
            assert!(
                (r.displacement_derivative[0] - frozen.displacement_derivative[0]).abs()
                    > 1e-5 * dq.abs()
            );
        }
        // Finite velocity loading contributes to q_w at DC even though it did
        // not influence the static solve or its numerical time scale.
        let p = PreparedHarmonicResponse::try_with_port_impedance(
            &[2.0],
            &[0.0],
            &[8.0],
            1,
            &[3.0],
            1,
            &[C64::from_re(1e250)],
            0.0,
            controls(),
            ports(),
            cx,
        )
        .unwrap();
        let r = p
            .solve_with_frequency_derivative(
                &[C64::ONE],
                None,
                Some(&[C64::ZERO]),
                tangent_controls(),
                cx,
            )
            .unwrap();
        close(
            r.displacement_derivative[0],
            C64::new(0.0, 9e250 / 64.0),
            3e-14,
        );
        close(r.velocity_derivative[0], C64::new(0.0, -0.125), 2e-15);
    });
}

#[test]
fn retained_maxwell_load_matches_an_independent_rational_derivative_and_work() {
    with_cx(&CancelGate::new(), |cx| {
        let (m, c, k, branch, tau) = (2.0, 0.3, 8.0, 5.0, 0.4);
        let f = C64::new(0.7, -0.2);
        let df = C64::new(-0.03, 0.04);
        for w in [0.0, 0.2, 1.0, 2.0, 3.5, 20.0] {
            let d = C64::new(1.0, -w * tau);
            let z = C64::from_re(branch * tau) / d;
            let dz = C64::new(0.0, branch * tau * tau) / (d * d);
            let p = PreparedHarmonicResponse::try_with_port_impedance(
                &[m],
                &[c],
                &[k],
                1,
                &[1.0],
                1,
                &[z],
                w,
                controls(),
                ports(),
                cx,
            )
            .unwrap();
            let r = p
                .solve_with_frequency_derivative(
                    &[f],
                    Some(&[df]),
                    Some(&[dz]),
                    tangent_controls(),
                    cx,
                )
                .unwrap();
            let bare = C64::new(k - m * w * w, -c * w);
            let denominator = bare * d + C64::new(0.0, -w * branch * tau);
            let denominator_derivative = C64::new(-2.0 * m * w, -c) * d
                + bare * C64::new(0.0, -tau)
                + C64::new(0.0, -branch * tau);
            let dq = df * d / denominator
                + f * (C64::new(0.0, -tau) * denominator - d * denominator_derivative)
                    / (denominator * denominator);
            close(r.displacement_derivative[0], dq, 7e-14);
            let primal = r.response.port_load.as_ref().unwrap();
            let derivative = r.port_load.as_ref().unwrap();
            let expected_heat = 0.5 * dz.re * primal.flow[0].norm_sq()
                + z.re * pair(primal.flow[0], derivative.flow_derivative[0]);
            scalar(
                derivative.average_absorbed_power_derivative,
                expected_heat,
                8e-13,
            );
            assert!(
                r.power_balance_residual_derivative.abs() < 2e-13 * (1.0 + expected_heat.abs())
            );
        }
    });
}

#[test]
fn full_mutual_nonproportional_and_nonreciprocal_tangents_match_cramers_rule() {
    with_cx(&CancelGate::new(), |cx| {
        let m = [2.0, 0.15, 0.15, 1.1];
        let c = [0.3, 0.04, 0.04, 0.2];
        let k = [9.0, -0.8, -0.8, 5.0];
        let b = [1.0, 0.2, -0.3, 1.5];
        let z = [
            C64::new(2.0, 0.4),
            C64::new(0.3, 0.2),
            C64::new(-0.1, 0.1),
            C64::new(1.0, -0.2),
        ];
        let dz = [
            C64::new(-0.1, -0.7),
            C64::new(0.2, 0.1),
            C64::new(0.03, -0.2),
            C64::new(0.1, 0.4),
        ];
        let f = [C64::new(0.8, -0.1), C64::new(-0.3, 0.7)];
        let df = [C64::new(0.1, 0.2), C64::new(0.05, -0.3)];
        for w in [0.0, 0.7, 1.9, 3.0, 10.0] {
            let p = PreparedHarmonicResponse::try_with_port_impedance(
                &m,
                &c,
                &k,
                2,
                &b,
                2,
                &z,
                w,
                controls(),
                ports(),
                cx,
            )
            .unwrap();
            let r = p
                .solve_with_frequency_derivative(&f, Some(&df), Some(&dz), tangent_controls(), cx)
                .unwrap();
            let mut a = [C64::ZERO; 4];
            let mut da = [C64::ZERO; 4];
            for i in 0..2 {
                for j in 0..2 {
                    let at = i * 2 + j;
                    a[at] = C64::new(k[at] - w * w * m[at], -w * c[at]);
                    da[at] = C64::new(-2.0 * w * m[at], -c[at]);
                    for p in 0..2 {
                        for q in 0..2 {
                            let map = b[i * 2 + p] * b[j * 2 + q];
                            a[at] = a[at] + C64::new(0.0, -w) * z[p * 2 + q].scale(map);
                            da[at] = da[at]
                                + C64::new(0.0, -1.0)
                                    * (z[p * 2 + q] + dz[p * 2 + q].scale(w)).scale(map);
                        }
                    }
                }
            }
            let determinant = a[0] * a[3] - a[1] * a[2];
            let ddet = da[0] * a[3] + a[0] * da[3] - da[1] * a[2] - a[1] * da[2];
            let numerator = [a[3] * f[0] - a[1] * f[1], a[0] * f[1] - a[2] * f[0]];
            let dn = [
                da[3] * f[0] + a[3] * df[0] - da[1] * f[1] - a[1] * df[1],
                da[0] * f[1] + a[0] * df[1] - da[2] * f[0] - a[2] * df[0],
            ];
            for i in 0..2 {
                close(
                    r.displacement_derivative[i],
                    (dn[i] * determinant - numerator[i] * ddet) / (determinant * determinant),
                    3e-13,
                );
            }
            assert!(r.relative_derivative_residual < 5e-14);
            assert!(
                r.power_balance_residual_derivative.abs()
                    < 4e-13 * (1.0 + r.average_input_power_derivative.abs())
            );
            if w == 1.9 {
                // q = T*q_new, B_new = T*B*S, Z_new = S^-1*Z*S^-1.
                // Frequencies stay fixed; flow_new=S*flow and effort_new=S^-1*effort.
                let coordinate = [2.0_f64.powi(200), 2.0_f64.powi(-200)];
                let port_unit = [2.0_f64.powi(-100), 2.0_f64.powi(100)];
                let congruence = |source: &[f64; 4]| -> [f64; 4] {
                    core::array::from_fn(|at| source[at] * coordinate[at / 2] * coordinate[at % 2])
                };
                let mapped: [f64; 4] =
                    core::array::from_fn(|at| b[at] * coordinate[at / 2] * port_unit[at % 2]);
                let load = |source: &[C64; 4]| -> [C64; 4] {
                    core::array::from_fn(|at| {
                        source[at]
                            .scale(1.0 / port_unit[at / 2])
                            .scale(1.0 / port_unit[at % 2])
                    })
                };
                let transformed = PreparedHarmonicResponse::try_with_port_impedance(
                    &congruence(&m),
                    &congruence(&c),
                    &congruence(&k),
                    2,
                    &mapped,
                    2,
                    &load(&z),
                    w,
                    controls(),
                    ports(),
                    cx,
                )
                .unwrap();
                let force: [C64; 2] = core::array::from_fn(|i| f[i].scale(coordinate[i]));
                let slope: [C64; 2] = core::array::from_fn(|i| df[i].scale(coordinate[i]));
                let changed = transformed
                    .solve_with_frequency_derivative(
                        &force,
                        Some(&slope),
                        Some(&load(&dz)),
                        tangent_controls(),
                        cx,
                    )
                    .unwrap();
                for i in 0..2 {
                    close(
                        changed.displacement_derivative[i].scale(coordinate[i]),
                        r.displacement_derivative[i],
                        3e-13,
                    );
                    close(
                        changed.port_load.as_ref().unwrap().flow_derivative[i]
                            .scale(1.0 / port_unit[i]),
                        r.port_load.as_ref().unwrap().flow_derivative[i],
                        3e-13,
                    );
                    close(
                        changed.port_load.as_ref().unwrap().effort_derivative[i]
                            .scale(port_unit[i]),
                        r.port_load.as_ref().unwrap().effort_derivative[i],
                        3e-13,
                    );
                }
                scalar(
                    changed.average_input_power_derivative,
                    r.average_input_power_derivative,
                    3e-13,
                );
                scalar(
                    changed.port_load.unwrap().average_absorbed_power_derivative,
                    r.port_load
                        .as_ref()
                        .unwrap()
                        .average_absorbed_power_derivative,
                    3e-13,
                );
            }
        }
    });
}

#[test]
fn exact_antiresonance_and_zero_source_have_regular_frequency_tangents() {
    with_cx(&CancelGate::new(), |cx| {
        let p = PreparedHarmonicResponse::try_new(
            &[1.0, 0.0, 0.0, 1.0],
            &[0.3, 0.0, 0.0, 0.0],
            &[7.0, -1.0, -1.0, 4.0],
            2,
            2.0,
            controls(),
            cx,
        )
        .unwrap();
        let r = p
            .solve_with_frequency_derivative(
                &[C64::ONE, C64::ZERO],
                None,
                None,
                tangent_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(r.response.displacement[0], C64::ZERO);
        close(r.displacement_derivative[0], C64::from_re(4.0), 3e-15);
        close(r.velocity_derivative[0], C64::new(0.0, -8.0), 3e-15);
        let df = [C64::new(0.2, 0.4), C64::new(-0.5, 0.1)];
        let zero = p
            .solve_with_frequency_derivative(
                &[C64::ZERO; 2],
                Some(&df),
                None,
                tangent_controls(),
                cx,
            )
            .unwrap();
        let forced = p.solve(&df, cx).unwrap();
        assert_eq!(zero.response.displacement, vec![C64::ZERO; 2]);
        assert_eq!(zero.displacement_derivative, forced.displacement);
        assert_eq!(zero.velocity_derivative, forced.velocity);
        assert_eq!(zero.average_input_power_derivative, 0.0);
    });
}

#[test]
fn physical_coordinate_and_time_units_preserve_frequency_derivatives() {
    with_cx(&CancelGate::new(), |cx| {
        let (m, c, k, w) = (2.0, 0.4, 18.0, 1.2);
        let f = C64::new(0.8, -0.3);
        let p = PreparedHarmonicResponse::try_new(&[m], &[c], &[k], 1, w, controls(), cx).unwrap();
        let reference = p
            .solve_with_frequency_derivative(&[f], None, None, tangent_controls(), cx)
            .unwrap();
        for exponent in [-300, 300] {
            let unit = 2.0_f64.powi(exponent);
            let p = PreparedHarmonicResponse::try_new(
                &[m * unit * unit],
                &[c * unit * unit],
                &[k * unit * unit],
                1,
                w,
                controls(),
                cx,
            )
            .unwrap();
            let r = p
                .solve_with_frequency_derivative(
                    &[f.scale(unit)],
                    None,
                    None,
                    tangent_controls(),
                    cx,
                )
                .unwrap();
            close(
                r.displacement_derivative[0].scale(unit),
                reference.displacement_derivative[0],
                4e-14,
            );
            close(
                r.velocity_derivative[0].scale(unit),
                reference.velocity_derivative[0],
                4e-14,
            );
            scalar(
                r.average_viscous_power_derivative,
                reference.average_viscous_power_derivative,
                4e-14,
            );
        }
        for exponent in [-400, 400] {
            let rate = 2.0_f64.powi(exponent);
            let p = PreparedHarmonicResponse::try_new(
                &[m / rate],
                &[c],
                &[k * rate],
                1,
                w * rate,
                controls(),
                cx,
            )
            .unwrap();
            let r = p
                .solve_with_frequency_derivative(&[f], None, None, tangent_controls(), cx)
                .unwrap();
            close(
                r.displacement_derivative[0].scale(rate).scale(rate),
                reference.displacement_derivative[0],
                4e-14,
            );
            close(
                r.velocity_derivative[0].scale(rate),
                reference.velocity_derivative[0],
                4e-14,
            );
            scalar(
                r.average_viscous_power_derivative * rate,
                reference.average_viscous_power_derivative,
                4e-14,
            );
        }
    });
}

#[test]
fn missing_load_slope_budgets_cancellation_and_bad_inputs_refuse_without_mutation() {
    with_cx(&CancelGate::new(), |cx| {
        let p = PreparedHarmonicResponse::try_with_port_impedance(
            &[2.0],
            &[0.3],
            &[8.0],
            1,
            &[1.0],
            1,
            &[C64::new(1.0, 0.2)],
            2.0,
            controls(),
            ports(),
            cx,
        )
        .unwrap();
        let force = [C64::new(0.8, -0.3)];
        let slope = [C64::new(0.1, -0.4)];
        let before = p.solve(&force, cx).unwrap();
        assert!(
            p.solve_with_frequency_derivative(&force, None, None, tangent_controls(), cx)
                .is_err()
        );
        assert!(
            p.solve_with_frequency_derivative(&force, None, Some(&[]), tangent_controls(), cx)
                .is_err()
        );
        assert!(
            p.solve_with_frequency_derivative(
                &force,
                None,
                Some(&[C64::new(f64::NAN, 0.0)]),
                tangent_controls(),
                cx
            )
            .is_err()
        );
        assert!(
            p.solve_with_frequency_derivative(
                &force,
                Some(&[]),
                Some(&slope),
                tangent_controls(),
                cx
            )
            .is_err()
        );
        let mut limited = tangent_controls();
        limited.maximum_additional_payload_bytes = 0;
        assert!(matches!(
            p.solve_with_frequency_derivative(&force, None, Some(&slope), limited, cx),
            Err(HarmonicResponseError::Budget { .. })
        ));
        limited = tangent_controls();
        limited.maximum_work = 0;
        assert!(matches!(
            p.solve_with_frequency_derivative(&force, None, Some(&slope), limited, cx),
            Err(HarmonicResponseError::WorkBudget { .. })
        ));
        let gate = CancelGate::new();
        gate.request();
        with_cx(&gate, |cancelled| {
            assert!(matches!(
                p.solve_with_frequency_derivative(
                    &force,
                    None,
                    Some(&slope),
                    tangent_controls(),
                    cancelled
                ),
                Err(HarmonicResponseError::Cancelled)
            ))
        });
        let first = p
            .solve_with_frequency_derivative(&force, None, Some(&slope), tangent_controls(), cx)
            .unwrap();
        let second = p
            .solve_with_frequency_derivative(&force, None, Some(&slope), tangent_controls(), cx)
            .unwrap();
        assert_eq!(
            first.displacement_derivative,
            second.displacement_derivative
        );
        assert_eq!(first.velocity_derivative, second.velocity_derivative);
        assert_eq!(first.response.displacement, before.displacement);
        assert_eq!(p.solve(&force, cx).unwrap().velocity, before.velocity);
        let unloaded =
            PreparedHarmonicResponse::try_new(&[1.0], &[0.2], &[4.0], 1, 1.0, controls(), cx)
                .unwrap();
        assert!(
            unloaded
                .solve_with_frequency_derivative(&force, None, Some(&slope), tangent_controls(), cx)
                .is_err()
        );
    });
}

#[test]
fn large_diagonal_modal_frequency_tangent_uses_linear_work_and_preserves_narrow_peaks() {
    with_cx(&CancelGate::new(), |cx| {
        let n = 129;
        let natural: Vec<_> = (0..n)
            .map(|i| if i == 0 { 0.1 } else { i as f64 + 1.0 })
            .collect();
        let mut damping = vec![0.0; n * n];
        for i in 0..n {
            damping[i * n + i] = if i == 0 {
                2e-13
            } else {
                0.05 + i as f64 / 100.0
            };
        }
        let force: Vec<_> = (0..n)
            .map(|i| C64::new(1.0 + (i % 3) as f64, -0.25))
            .collect();
        let slope: Vec<_> = (0..n)
            .map(|i| C64::new(0.2, (i % 4) as f64 / 8.0))
            .collect();
        let compact = HarmonicResponseControls {
            maximum_coordinates: n,
            maximum_payload_bytes: 272 * n,
            ..controls()
        };
        let linear = HarmonicFrequencyControls {
            maximum_additional_payload_bytes: 512 * n,
            maximum_work: 320 * n as u64,
        };
        for w in [0.0, 0.1 - 2e-13, 0.1, 0.1 + 2e-13, 9.25] {
            let prepared =
                PreparedHarmonicResponse::try_modal(&natural, &damping, w, compact, cx).unwrap();
            let result = prepared
                .solve_with_frequency_derivative(&force, Some(&slope), None, linear, cx)
                .unwrap();
            let mut expected_input = 0.0;
            let mut expected_viscous = 0.0;
            let mut power_envelope = 0.0;
            for i in 0..n {
                let c = damping[i * n + i];
                let a = C64::new((natural[i] - w) * (natural[i] + w), -w * c);
                let q = force[i] / a;
                let dq = (slope[i] + C64::new(2.0 * w, c) * q) / a;
                let v = C64::new(0.0, -w) * q;
                let dv = C64::new(0.0, -1.0) * (q + dq.scale(w));
                close(result.response.displacement[i], q, 6e-14);
                close(result.displacement_derivative[i], dq, 8e-14);
                close(result.velocity_derivative[i], dv, 8e-14);
                expected_input += 0.5 * (pair(slope[i], v) + pair(force[i], dv));
                expected_viscous += c * pair(v, dv);
                power_envelope += 0.5 * (slope[i].abs() * v.abs() + force[i].abs() * dv.abs())
                    + c * v.abs() * dv.abs();
            }
            // At an extremely narrow peak the signed power slope cancels
            // large terms; use their uncancelled physical envelope as scale.
            let power_tolerance = 2e-12 * power_envelope.max(1.0);
            assert!(
                (result.average_input_power_derivative - expected_input).abs() < power_tolerance
            );
            assert!(
                (result.average_viscous_power_derivative - expected_viscous).abs()
                    < power_tolerance
            );
            assert!(result.power_balance_residual_derivative.abs() < power_tolerance);
            assert!(result.relative_derivative_residual < 3e-14);
            let insufficient_work = HarmonicFrequencyControls {
                maximum_work: linear.maximum_work - 1,
                ..linear
            };
            assert!(matches!(
                prepared.solve_with_frequency_derivative(
                    &force, Some(&slope), None, insufficient_work, cx
                ),
                Err(HarmonicResponseError::WorkBudget { requested, .. })
                    if requested == linear.maximum_work
            ));
        }
    });
}

#[test]
fn compact_modal_and_dense_frequency_responses_agree_at_dc_and_resonance() {
    with_cx(&CancelGate::new(), |cx| {
        let natural = [1.0, 2.0, 3.0, 4.0];
        let mut mass = [0.0; 16];
        let mut damping = [0.0; 16];
        let mut stiffness = [0.0; 16];
        for i in 0..4 {
            mass[i * 4 + i] = 1.0;
            damping[i * 4 + i] = [0.125, 0.0, 0.25, 1e-9][i];
            stiffness[i * 4 + i] = natural[i] * natural[i];
        }
        let force = [C64::new(1.0, -0.2), C64::ONE, C64::new(-0.3, 0.7), C64::ONE];
        let slope = [C64::new(0.1, -0.3); 4];
        for w in [0.0, 0.3, 1.0, 3.1, 4.0] {
            let compact =
                PreparedHarmonicResponse::try_modal(&natural, &damping, w, controls(), cx)
                    .unwrap()
                    .solve_with_frequency_derivative(
                        &force,
                        Some(&slope),
                        None,
                        tangent_controls(),
                        cx,
                    )
                    .unwrap();
            let dense = PreparedHarmonicResponse::try_new(
                &mass,
                &damping,
                &stiffness,
                4,
                w,
                controls(),
                cx,
            )
            .unwrap()
            .solve_with_frequency_derivative(&force, Some(&slope), None, tangent_controls(), cx)
            .unwrap();
            for i in 0..4 {
                close(
                    compact.response.displacement[i],
                    dense.response.displacement[i],
                    6e-14,
                );
                close(
                    compact.displacement_derivative[i],
                    dense.displacement_derivative[i],
                    8e-14,
                );
                close(
                    compact.velocity_derivative[i],
                    dense.velocity_derivative[i],
                    8e-14,
                );
            }
            assert!(
                (compact.average_viscous_power_derivative - dense.average_viscous_power_derivative)
                    .abs()
                    <= 3e-12 * dense.average_viscous_power_derivative.abs().max(1.0)
            );
        }
    });
}

#[test]
fn tiny_off_diagonal_damping_remains_in_the_dc_frequency_derivative() {
    with_cx(&CancelGate::new(), |cx| {
        let coupling = 1e-18;
        let prepared = PreparedHarmonicResponse::try_modal(
            &[1.0, 2.0],
            &[0.25, coupling, coupling, 0.5],
            0.0,
            controls(),
            cx,
        )
        .unwrap();
        let result = prepared
            .solve_with_frequency_derivative(
                &[C64::ONE, C64::ZERO],
                None,
                None,
                tangent_controls(),
                cx,
            )
            .unwrap();
        // The DC dynamic operator is diagonal, but its derivative is not:
        // K q_w = i C q. Compactness must be decided from the physical pencil.
        assert_eq!(result.response.displacement[1], C64::ZERO);
        close(
            result.displacement_derivative[0],
            C64::new(0.0, 0.25),
            5e-14,
        );
        close(
            result.displacement_derivative[1],
            C64::new(0.0, coupling / 4.0),
            5e-14,
        );
        assert!(result.displacement_derivative[1].abs() > 0.0);
    });
}

#[test]
fn explicit_diagonal_modal_matches_matrix_input_and_frequency_tangents() {
    with_cx(&CancelGate::new(), |cx| {
        let natural = [0.1, 1.0, 2.0, 3.0];
        let rates = [2e-13, 0.125, 0.0, 0.5];
        let mut full_damping = [0.0; 16];
        for i in 0..4 {
            full_damping[i * 4 + i] = rates[i];
        }
        let force = [C64::ONE, C64::new(0.3, -0.7), C64::ONE, C64::new(-0.2, 0.5)];
        let slope = [C64::new(0.1, -0.05); 4];
        for w in [0.0, 0.1 - 2e-13, 0.1, 0.1 + 2e-13, 1.5] {
            let diagonal = PreparedHarmonicResponse::try_diagonal_modal(
                &natural, &rates, w, controls(), cx,
            )
            .unwrap()
            .solve_with_frequency_derivative(
                &force, Some(&slope), None, tangent_controls(), cx,
            )
            .unwrap();
            let matrix = PreparedHarmonicResponse::try_modal(
                &natural, &full_damping, w, controls(), cx,
            )
            .unwrap()
            .solve_with_frequency_derivative(
                &force, Some(&slope), None, tangent_controls(), cx,
            )
            .unwrap();
            assert_eq!(diagonal.response.displacement, matrix.response.displacement);
            assert_eq!(diagonal.response.velocity, matrix.response.velocity);
            assert_eq!(diagonal.displacement_derivative, matrix.displacement_derivative);
            assert_eq!(diagonal.velocity_derivative, matrix.velocity_derivative);
            assert_eq!(
                diagonal.response.average_viscous_power_w.to_bits(),
                matrix.response.average_viscous_power_w.to_bits()
            );
            assert_eq!(
                diagonal.average_viscous_power_derivative.to_bits(),
                matrix.average_viscous_power_derivative.to_bits()
            );
            assert_eq!(
                diagonal.relative_derivative_residual.to_bits(),
                matrix.relative_derivative_residual.to_bits()
            );
        }
    });
}
