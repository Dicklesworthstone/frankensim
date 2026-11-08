//! Original-coordinate second frequency derivatives: resonance, dispersion and work.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicFrequencyCurvatureControls, HarmonicPortControls, HarmonicResponseControls,
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
        maximum_coordinates: 256,
        maximum_payload_bytes: 1 << 24,
        maximum_relative_residual: 2e-11,
    }
}
fn curve() -> HarmonicFrequencyCurvatureControls {
    HarmonicFrequencyCurvatureControls {
        maximum_additional_payload_bytes: 1 << 24,
        maximum_work: 10_000_000,
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
fn pair(a: C64, b: C64) -> f64 {
    (a.conj() * b).re
}

#[test]
fn scalar_source_jet_and_dc_kinematics_match_closed_rational_derivative() {
    with_cx(&CancelGate::new(), |cx| {
        let (m, c, k) = (2.0_f64, 0.3, 8.0);
        let (f, df, ddf) = (
            C64::new(1.7, -0.6),
            C64::new(-0.2, 0.8),
            C64::new(0.3, -0.4),
        );
        for w in [0.0, 0.1, 1.999, 2.0, 2.001, 3.0, 25.0] {
            let prepared =
                PreparedHarmonicResponse::try_new(&[m], &[c], &[k], 1, w, controls(), cx).unwrap();
            let result = prepared
                .solve_with_frequency_curvature(
                    &[f],
                    Some(&[df]),
                    Some(&[ddf]),
                    None,
                    None,
                    curve(),
                    cx,
                )
                .unwrap();
            let d = C64::new(k - m * w * w, -c * w);
            let a = C64::new(2.0 * m * w, c);
            let q = f / d;
            let dq = df / d + f * a / (d * d);
            // Independent quotient-rule expression, not the implicit RHS.
            let ddq = ddf / d
                + (df * a).scale(2.0) / (d * d)
                + f.scale(2.0 * m) / (d * d)
                + (f * a * a).scale(2.0) / (d * d * d);
            let v = C64::new(0.0, -w) * q;
            let dv = C64::new(0.0, -1.0) * (q + dq.scale(w));
            let ddv = C64::new(0.0, -1.0) * (dq.scale(2.0) + ddq.scale(w));
            close(result.displacement_second_derivative[0], ddq, 3e-12);
            close(result.velocity_second_derivative[0], ddv, 3e-12);
            let source = 0.5 * (pair(ddf, v) + 2.0 * pair(df, dv) + pair(f, ddv));
            let viscous = c * (pair(dv, dv) + pair(v, ddv));
            let envelope = 1.0 + source.abs() + c * (dv.abs() * dv.abs() + v.abs() * ddv.abs());
            assert!(
                (result.average_input_power_second_derivative - source).abs() < 3e-12 * envelope
            );
            assert!(
                (result.average_viscous_power_second_derivative - viscous).abs() < 3e-12 * envelope
            );
            assert!(result.power_balance_residual_second_derivative.abs() < 5e-12 * envelope);
            assert!(result.relative_second_derivative_residual < 2e-13);
        }
    });
}

#[test]
fn dispersive_spring_and_inertance_curvature_matches_actual_added_coefficients() {
    with_cx(&CancelGate::new(), |cx| {
        let (lc, lm, lk) = (0.4, 0.7, 0.2);
        let (m, c, k) = (2.0 + 9.0 * lm, 0.1 + 9.0 * lc, 8.0 + 9.0 * lk);
        let f = C64::new(3.0, -0.7);
        for w in [0.01, 0.7, 1.1, 2.0, 4.0, 100.0] {
            let z = C64::new(lc, lk / w - w * lm);
            let dz = C64::new(0.0, -lk / (w * w) - lm);
            let ddz = C64::new(0.0, 2.0 * lk / (w * w * w));
            let loaded = PreparedHarmonicResponse::try_with_port_impedance(
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
            let actual = loaded
                .solve_with_frequency_curvature(
                    &[f],
                    None,
                    None,
                    Some(&[dz]),
                    Some(&[ddz]),
                    curve(),
                    cx,
                )
                .unwrap();
            let equivalent =
                PreparedHarmonicResponse::try_new(&[m], &[c], &[k], 1, w, controls(), cx)
                    .unwrap()
                    .solve_with_frequency_curvature(&[f], None, None, None, None, curve(), cx)
                    .unwrap();
            close(
                actual.displacement_second_derivative[0],
                equivalent.displacement_second_derivative[0],
                2e-11,
            );
            close(
                actual.velocity_second_derivative[0],
                equivalent.velocity_second_derivative[0],
                2e-11,
            );
            let primal = actual
                .frequency_response
                .response
                .port_load
                .as_ref()
                .unwrap();
            let slope = actual.frequency_response.port_load.as_ref().unwrap();
            let curvature = actual.port_load.as_ref().unwrap();
            let heat = lc
                * (pair(slope.flow_derivative[0], slope.flow_derivative[0])
                    + pair(primal.flow[0], curvature.flow_second_derivative[0]));
            let bound = 1.0
                + lc * (slope.flow_derivative[0].abs().powi(2)
                    + primal.flow[0].abs() * curvature.flow_second_derivative[0].abs());
            assert!(
                (curvature.average_absorbed_power_second_derivative - heat).abs() < 2e-11 * bound
            );
            assert!(actual.power_balance_residual_second_derivative.abs() < 3e-11 * bound);
        }
    });
}

#[test]
fn maxwell_memory_curvature_matches_independent_cubic_quotient() {
    with_cx(&CancelGate::new(), |cx| {
        let (m, c, k, kb, tau) = (1.7, 0.06, 4.3, 2.1, 0.8);
        let f = C64::new(0.7, 0.2);
        for w in [0.0, 0.25, 1.7, 2.0, 10.0] {
            let n = C64::new(1.0, -w * tau);
            let dn = C64::new(0.0, -tau);
            let z = C64::from_re(kb * tau) / n;
            let dz = C64::new(0.0, kb * tau * tau) / (n * n);
            let ddz = C64::from_re(-2.0 * kb * tau * tau * tau) / (n * n * n);
            let result = PreparedHarmonicResponse::try_with_port_impedance(
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
            .unwrap()
            .solve_with_frequency_curvature(
                &[f],
                None,
                None,
                Some(&[dz]),
                Some(&[ddz]),
                curve(),
                cx,
            )
            .unwrap();
            let s = C64::new(k - m * w * w, -c * w);
            let ds = C64::new(-2.0 * m * w, -c);
            let p = s * n - C64::new(0.0, w * kb * tau);
            let dp = ds * n + s * dn - C64::new(0.0, kb * tau);
            let ddp = n.scale(-2.0 * m) + (ds * dn).scale(2.0);
            let ddq = f
                * ((dn * dp).scale(-2.0) / (p * p) - n * ddp / (p * p)
                    + (n * dp * dp).scale(2.0) / (p * p * p));
            close(result.displacement_second_derivative[0], ddq, 3e-12);
            let bound = 1.0 + result.average_input_power_second_derivative.abs();
            assert!(result.power_balance_residual_second_derivative.abs() < 5e-12 * bound);
        }
    });
}

#[test]
fn gyroscopic_dc_curvature_keeps_signed_coupling_without_heat() {
    with_cx(&CancelGate::new(), |cx| {
        let f = [C64::new(0.7, -0.2), C64::new(-0.4, 0.8)];
        let g = 0.6;
        let result = PreparedHarmonicResponse::try_new(
            &[1.0, 0.0, 0.0, 1.0],
            &[0.0, g, -g, 0.0],
            &[4.0, 0.0, 0.0, 9.0],
            2,
            0.0,
            controls(),
            cx,
        )
        .unwrap()
        .solve_with_frequency_curvature(&f, None, None, None, None, curve(), cx)
        .unwrap();
        let q = [f[0].scale(0.25), f[1].scale(1.0 / 9.0)];
        close(
            result.displacement_second_derivative[0],
            q[0].scale(0.5 * (1.0 + g * g / 9.0)),
            2e-13,
        );
        close(
            result.displacement_second_derivative[1],
            q[1].scale((2.0 / 9.0) * (1.0 + g * g / 4.0)),
            2e-13,
        );
        close(
            result.velocity_second_derivative[0],
            q[1].scale(g / 2.0),
            2e-13,
        );
        close(
            result.velocity_second_derivative[1],
            q[0].scale(-2.0 * g / 9.0),
            2e-13,
        );
        assert_eq!(result.average_viscous_power_second_derivative, 0.0);
        assert!(result.average_input_power_second_derivative.abs() < 3e-17);
    });
}

#[test]
fn high_q_velocity_peak_retains_finite_extreme_power_curvature() {
    with_cx(&CancelGate::new(), |cx| {
        let c = 2.0_f64.powi(-100);
        let f = C64::from_re(2.0_f64.powi(350));
        let result =
            PreparedHarmonicResponse::try_new(&[1.0], &[c], &[1.0], 1, 1.0, controls(), cx)
                .unwrap()
                .solve_with_frequency_curvature(&[f], None, None, None, None, curve(), cx)
                .unwrap();
        // P(w)=|F|²*c*w²/[2((1-w²)²+c²*w²)]; P''(1)=-4|F|²/c³.
        let expected = -2.0_f64.powi(1002);
        assert!((result.average_input_power_second_derivative / expected - 1.0).abs() < 3e-13);
        assert!((result.average_viscous_power_second_derivative / expected - 1.0).abs() < 3e-13);
        assert!(result.power_balance_residual_second_derivative.abs() <= expected.abs() * 3e-13);
        assert!(result.relative_second_derivative_residual < 2e-13);
    });
}

#[test]
fn diagonal_bank_preserves_narrow_peaks_under_linear_curvature_work_budget() {
    with_cx(&CancelGate::new(), |cx| {
        let n = 129;
        let frequencies: Vec<_> = (0..n).map(|i| i as f64 + 1.0).collect();
        let damping = vec![2e-10; n];
        let f = vec![C64::ONE; n];
        let limits = HarmonicFrequencyCurvatureControls {
            maximum_additional_payload_bytes: 1536 * n,
            maximum_work: 1280 * n as u64,
        };
        for w in [0.0, 1.0, 1.0 + 2e-10, 7.25] {
            let prepared = PreparedHarmonicResponse::try_diagonal_modal(
                &frequencies,
                &damping,
                w,
                controls(),
                cx,
            )
            .unwrap();
            let result = prepared
                .solve_with_frequency_curvature(&f, None, None, None, None, limits, cx)
                .unwrap();
            for i in 0..n {
                let d = C64::new((frequencies[i] - w) * (frequencies[i] + w), -w * damping[i]);
                let a = C64::new(2.0 * w, damping[i]);
                let expected = C64::from_re(2.0) / (d * d) + (a * a).scale(2.0) / (d * d * d);
                close(result.displacement_second_derivative[i], expected, 3e-12);
            }
            assert!(
                matches!(prepared.solve_with_frequency_curvature(&f,None,None,None,None,
                HarmonicFrequencyCurvatureControls { maximum_work: limits.maximum_work-1, ..limits },cx),
                Err(HarmonicResponseError::WorkBudget{requested,..}) if requested==limits.maximum_work)
            );
        }
    });
}

#[test]
fn compact_mutual_loaded_curvature_matches_dense_original_pencil() {
    with_cx(&CancelGate::new(), |cx| {
        let natural = [1.0, 2.0, 3.0];
        let damping = [0.0, 0.04, 0.08];
        let b = [1.0, 0.2, -0.4, 0.6, 0.7, -0.3];
        let z = [
            C64::new(0.4, 0.1),
            C64::new(0.06, 0.02),
            C64::new(0.06, 0.02),
            C64::new(0.3, -0.2),
        ];
        let dz = [
            C64::new(0.03, -0.05),
            C64::new(-0.01, 0.02),
            C64::new(0.02, 0.01),
            C64::new(0.01, 0.07),
        ];
        let ddz = [
            C64::new(-0.04, 0.02),
            C64::new(0.01, -0.03),
            C64::new(-0.02, 0.04),
            C64::new(0.03, 0.01),
        ];
        let f = [C64::new(1.0, 0.2), C64::new(-0.3, 0.8), C64::new(0.7, -0.5)];
        for w in [0.0, 1.0, 1.000001, 2.2] {
            let compact = PreparedHarmonicResponse::try_diagonal_modal_with_port_impedance(
                &natural,
                &damping,
                &b,
                2,
                &z,
                w,
                controls(),
                ports(),
                cx,
            )
            .unwrap();
            let dense = PreparedHarmonicResponse::try_with_port_impedance(
                &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
                &[0.0, 0.0, 0.0, 0.0, 0.04, 0.0, 0.0, 0.0, 0.08],
                &[1.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 9.0],
                3,
                &b,
                2,
                &z,
                w,
                controls(),
                ports(),
                cx,
            )
            .unwrap();
            let a = compact
                .solve_with_frequency_curvature(&f, None, None, Some(&dz), Some(&ddz), curve(), cx)
                .unwrap();
            let b = dense
                .solve_with_frequency_curvature(&f, None, None, Some(&dz), Some(&ddz), curve(), cx)
                .unwrap();
            for i in 0..3 {
                close(
                    a.displacement_second_derivative[i],
                    b.displacement_second_derivative[i],
                    5e-11,
                );
                close(
                    a.velocity_second_derivative[i],
                    b.velocity_second_derivative[i],
                    5e-11,
                );
            }
            let scale = 1.0 + b.average_input_power_second_derivative.abs();
            assert!(
                (a.average_input_power_second_derivative - b.average_input_power_second_derivative)
                    .abs()
                    < 5e-11 * scale
            );
            assert!(a.power_balance_residual_second_derivative.abs() < 5e-11 * scale);
        }
    });
}

#[test]
fn curvature_refusals_leave_original_factor_and_response_reusable() {
    with_cx(&CancelGate::new(), |cx| {
        let prepared = PreparedHarmonicResponse::try_with_port_impedance(
            &[1.0],
            &[0.1],
            &[4.0],
            1,
            &[1.0],
            1,
            &[C64::ONE],
            1.0,
            controls(),
            ports(),
            cx,
        )
        .unwrap();
        let before = prepared.solve(&[C64::ONE], cx).unwrap();
        let zero = [C64::ZERO];
        for (dz, ddz) in [(None, Some(zero.as_slice())), (Some(zero.as_slice()), None)] {
            assert!(
                prepared
                    .solve_with_frequency_curvature(&[C64::ONE], None, None, dz, ddz, curve(), cx)
                    .is_err()
            );
        }
        assert!(
            prepared
                .solve_with_frequency_curvature(
                    &[C64::ONE],
                    None,
                    Some(&[]),
                    Some(&zero),
                    Some(&zero),
                    curve(),
                    cx
                )
                .is_err()
        );
        assert!(matches!(
            prepared.solve_with_frequency_curvature(
                &[C64::ONE],
                None,
                None,
                Some(&zero),
                Some(&zero),
                HarmonicFrequencyCurvatureControls {
                    maximum_additional_payload_bytes: 0,
                    ..curve()
                },
                cx
            ),
            Err(HarmonicResponseError::Budget { .. })
        ));
        let gate = CancelGate::new();
        gate.request();
        with_cx(&gate, |cancelled| {
            assert!(matches!(
                prepared.solve_with_frequency_curvature(
                    &[C64::ONE],
                    None,
                    None,
                    Some(&zero),
                    Some(&zero),
                    curve(),
                    cancelled
                ),
                Err(HarmonicResponseError::Cancelled)
            ))
        });
        let first = prepared
            .solve_with_frequency_curvature(
                &[C64::ONE],
                None,
                None,
                Some(&zero),
                Some(&zero),
                curve(),
                cx,
            )
            .unwrap();
        let second = prepared
            .solve_with_frequency_curvature(
                &[C64::ONE],
                None,
                None,
                Some(&zero),
                Some(&zero),
                curve(),
                cx,
            )
            .unwrap();
        assert_eq!(
            first.displacement_second_derivative,
            second.displacement_second_derivative
        );
        assert_eq!(
            first.frequency_response.response.displacement,
            before.displacement
        );
        assert_eq!(
            prepared.solve(&[C64::ONE], cx).unwrap().velocity,
            before.velocity
        );
        assert!(matches!(
            PreparedHarmonicResponse::try_new(&[1.0], &[0.0], &[1.0], 1, 1.0, controls(), cx),
            Err(HarmonicResponseError::Singular { .. })
        ));
    });
}
