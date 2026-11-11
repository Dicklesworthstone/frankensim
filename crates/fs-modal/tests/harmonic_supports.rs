//! G1/G3 prescribed harmonic supports: cross inertia, phase, reactions and work.

#[path = "harmonic_supports/impedance.rs"]
mod impedance;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicResponse,
    PreparedSupportedHarmonicResponse,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 913,
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
fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 4e-12 * expected.abs().max(1.0),
        "{actual:?} != {expected:?}"
    );
}
fn scalar(actual: f64, expected: f64) {
    close(C64::from_re(actual), C64::from_re(expected));
}
fn inner(a: C64, b: C64) -> f64 {
    a.re * b.re + a.im * b.im
}
fn entry(k: f64, m: f64, c: f64, omega: f64) -> C64 {
    C64::new(k - omega * omega * m, -omega * c)
}

const M: [f64; 4] = [2.0, 0.4, 0.4, 1.0];
const K: [f64; 4] = [9.0, -9.0, -9.0, 9.0];
const C: [f64; 4] = [0.7, -0.2, -0.2, 0.4];

#[test]
fn g1_moving_base_retains_cross_mass_phase_reaction_and_original_work() {
    with_cx(&CancelGate::new(), |cx| {
        let force = [C64::new(0.7, -0.9), C64::new(0.4, 0.2)];
        let base = C64::new(0.3, -0.15);
        for omega in [0.0, 0.4, 2.9, 3.0, 4.0] {
            let p = PreparedSupportedHarmonicResponse::try_new(
                &M,
                &C,
                &K,
                2,
                &[0],
                omega,
                controls(),
                cx,
            )
            .unwrap();
            assert_eq!(p.free_dofs(), &[1]);
            assert_eq!(p.support_dofs(), &[0]);
            assert_eq!(p.angular_frequency_rad_s(), omega);
            let a = |i| entry(K[i], M[i], C[i], omega);
            let q = (force[1] - a(2) * base) / a(3);
            let reaction = a(0) * base + a(1) * q - force[0];
            let expected_v = [
                C64::new(base.im, -base.re).scale(omega),
                C64::new(q.im, -q.re).scale(omega),
            ];
            let result = p.solve(&force, &[base], cx).unwrap();
            close(result.displacement[0], base);
            close(result.displacement[1], q);
            close(result.velocity[0], expected_v[0]);
            close(result.velocity[1], expected_v[1]);
            close(result.support_reaction[0], reaction);
            scalar(
                result.average_input_power_w,
                0.5 * (inner(force[0], expected_v[0]) + inner(force[1], expected_v[1])),
            );
            scalar(
                result.average_support_power_w,
                0.5 * inner(reaction, expected_v[0]),
            );
            let heat = 0.5
                * (C[0] * expected_v[0].norm_sq()
                    + C[3] * expected_v[1].norm_sq()
                    + 2.0 * C[1] * inner(expected_v[0], expected_v[1]));
            scalar(result.average_viscous_power_w, heat);
            assert!(result.relative_residual < 1e-13);
            assert!(result.power_balance_residual_w.abs() < 1e-12 * heat.max(1.0));
            if omega == 0.0 {
                close(result.support_reaction[0], -force[0] - force[1]);
                assert_eq!(result.average_input_power_w, 0.0);
                assert_eq!(result.average_support_power_w, 0.0);
                assert_eq!(result.average_viscous_power_w, 0.0);
            }
        }
    });
}

#[test]
fn g3_support_applied_force_changes_reaction_and_work_without_changing_motion() {
    with_cx(&CancelGate::new(), |cx| {
        let p =
            PreparedSupportedHarmonicResponse::try_new(&M, &C, &K, 2, &[0], 1.7, controls(), cx)
                .unwrap();
        let base = [C64::new(-0.2, 0.5)];
        let a = p.solve(&[C64::ZERO, C64::ONE], &base, cx).unwrap();
        let extra = C64::new(2.0, -3.0);
        let b = p.solve(&[extra, C64::ONE], &base, cx).unwrap();
        assert_eq!(a.displacement, b.displacement);
        close(b.support_reaction[0], a.support_reaction[0] - extra);
        let added_power = 0.5 * inner(extra, a.velocity[0]);
        scalar(
            b.average_input_power_w - a.average_input_power_w,
            added_power,
        );
        scalar(
            b.average_support_power_w - a.average_support_power_w,
            -added_power,
        );
        scalar(b.average_viscous_power_w, a.average_viscous_power_w);
    });
}

#[test]
fn g1_relative_support_phases_and_gyroscopic_couplings_keep_full_work() {
    with_cx(&CancelGate::new(), |cx| {
        let m = [2.0, 0.2, 0.0, 0.2, 1.0, 0.1, 0.0, 0.1, 3.0];
        let k = [4.0, -4.0, 0.0, -4.0, 9.0, -5.0, 0.0, -5.0, 5.0];
        let c = [0.3, -0.1, 0.6, -0.1, 0.4, -0.15, -0.6, -0.15, 0.5];
        let omega = 2.0;
        let p = PreparedSupportedHarmonicResponse::try_new(
            &m,
            &c,
            &k,
            3,
            &[0, 2],
            omega,
            controls(),
            cx,
        )
        .unwrap();
        let force = [C64::ZERO, C64::new(0.4, 0.3), C64::ZERO];
        for phase in [C64::ONE, -C64::ONE, C64::new(0.0, 1.0)] {
            let base = [C64::new(0.3, -0.1), phase.scale(0.25)];
            let a = |i| entry(k[i], m[i], c[i], omega);
            let q = (force[1] - a(3) * base[0] - a(5) * base[1]) / a(4);
            let r = p.solve(&force, &base, cx).unwrap();
            close(r.displacement[1], q);
            for (support, &i) in [0, 2].iter().enumerate() {
                let reaction =
                    (0..3).fold(-force[i], |sum, j| sum + a(i * 3 + j) * r.displacement[j]);
                close(r.support_reaction[support], reaction);
            }
            let mut heat = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    heat +=
                        0.25 * (c[i * 3 + j] + c[j * 3 + i]) * inner(r.velocity[i], r.velocity[j]);
                }
            }
            scalar(r.average_viscous_power_w, heat);
            scalar(r.average_input_power_w + r.average_support_power_w, heat);
            assert!(r.power_balance_residual_w.abs() < 2e-13);
        }
    });
}

#[test]
fn g3_stationary_and_empty_support_maps_match_the_existing_owner() {
    with_cx(&CancelGate::new(), |cx| {
        let omega = 1.3;
        let p =
            PreparedSupportedHarmonicResponse::try_new(&M, &C, &K, 2, &[0], omega, controls(), cx)
                .unwrap();
        let free =
            PreparedHarmonicResponse::try_new(&[M[3]], &[C[3]], &[K[3]], 1, omega, controls(), cx)
                .unwrap();
        let f = C64::new(0.5, -0.2);
        let a = p.solve(&[C64::ZERO, f], &[C64::ZERO], cx).unwrap();
        let b = free.solve(&[f], cx).unwrap();
        close(a.displacement[1], b.displacement[0]);
        scalar(a.average_viscous_power_w, b.average_viscous_power_w);
        assert_eq!(a.average_support_power_w, 0.0);
        let p =
            PreparedSupportedHarmonicResponse::try_new(&M, &C, &K, 2, &[], omega, controls(), cx)
                .unwrap();
        let original =
            PreparedHarmonicResponse::try_new(&M, &C, &K, 2, omega, controls(), cx).unwrap();
        let force = [C64::new(0.1, 0.4), f];
        let a = p.solve(&force, &[], cx).unwrap();
        let b = original.solve(&force, cx).unwrap();
        for i in 0..2 {
            close(a.displacement[i], b.displacement[i]);
        }
        scalar(a.average_input_power_w, b.average_input_power_w);
        scalar(a.average_viscous_power_w, b.average_viscous_power_w);
    });
}

#[test]
fn g1_all_prescribed_coordinates_need_no_invertible_full_operator() {
    with_cx(&CancelGate::new(), |cx| {
        for omega in [0.0, 2.1] {
            let p = PreparedSupportedHarmonicResponse::try_new(
                &M,
                &C,
                &K,
                2,
                &[0, 1],
                omega,
                controls(),
                cx,
            )
            .unwrap();
            let q = [C64::new(0.2, 0.1), C64::new(-0.3, 0.4)];
            let r = p.solve(&[C64::ZERO; 2], &q, cx).unwrap();
            assert_eq!(r.displacement, q);
            assert_eq!(r.relative_residual, 0.0);
            for i in 0..2 {
                let expected = (0..2).fold(C64::ZERO, |sum, j| {
                    sum + entry(K[i * 2 + j], M[i * 2 + j], C[i * 2 + j], omega) * q[j]
                });
                close(r.support_reaction[i], expected);
            }
            scalar(r.average_support_power_w, r.average_viscous_power_w);
        }
    });
}

#[test]
fn g3_support_motion_reactions_and_work_are_invariant_to_coordinate_units() {
    with_cx(&CancelGate::new(), |cx| {
        let units = [2.0_f64.powi(120), 2.0_f64.powi(-100)];
        let congruence = |matrix: [f64; 4]| {
            std::array::from_fn::<_, 4, _>(|i| matrix[i] * units[i / 2] * units[i % 2])
        };
        let base = C64::new(0.1, -0.3);
        let force = [C64::new(0.3, -0.1), C64::new(-0.2, 0.6)];
        let p =
            PreparedSupportedHarmonicResponse::try_new(&M, &C, &K, 2, &[0], 2.7, controls(), cx)
                .unwrap();
        let expected = p.solve(&force, &[base], cx).unwrap();
        let p = PreparedSupportedHarmonicResponse::try_new(
            &congruence(M),
            &congruence(C),
            &congruence(K),
            2,
            &[0],
            2.7,
            controls(),
            cx,
        )
        .unwrap();
        let result = p
            .solve(
                &[force[0].scale(units[0]), force[1].scale(units[1])],
                &[base.scale(1.0 / units[0])],
                cx,
            )
            .unwrap();
        for i in 0..2 {
            close(
                result.displacement[i].scale(units[i]),
                expected.displacement[i],
            );
            close(result.velocity[i].scale(units[i]), expected.velocity[i]);
        }
        close(
            result.support_reaction[0].scale(1.0 / units[0]),
            expected.support_reaction[0],
        );
        scalar(result.average_input_power_w, expected.average_input_power_w);
        scalar(
            result.average_support_power_w,
            expected.average_support_power_w,
        );
        scalar(
            result.average_viscous_power_w,
            expected.average_viscous_power_w,
        );
    });
}

#[test]
fn g0_full_source_admission_free_poles_budgets_and_cancellation_refuse() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        for supports in [&[0, 0][..], &[1, 0][..], &[2][..]] {
            assert!(
                PreparedSupportedHarmonicResponse::try_new(
                    &M,
                    &C,
                    &K,
                    2,
                    supports,
                    1.0,
                    controls(),
                    cx
                )
                .is_err()
            );
        }
        assert!(matches!(
            PreparedSupportedHarmonicResponse::try_new(
                &M,
                &[0.0; 4],
                &K,
                2,
                &[0],
                3.0,
                controls(),
                cx
            ),
            Err(HarmonicResponseError::Singular { .. })
        ));
        assert!(
            PreparedSupportedHarmonicResponse::try_new(
                &[-1.0, 0.0, 0.0, 1.0],
                &C,
                &K,
                2,
                &[0],
                1.0,
                controls(),
                cx
            )
            .is_err()
        );
        assert!(
            PreparedSupportedHarmonicResponse::try_new(
                &M,
                &[-1.0, 0.0, 0.0, 0.4],
                &K,
                2,
                &[0],
                1.0,
                controls(),
                cx
            )
            .is_err()
        );
        let tight = HarmonicResponseControls {
            maximum_payload_bytes: 320 * 4 + 512 * 2 - 1,
            ..controls()
        };
        assert!(matches!(
            PreparedSupportedHarmonicResponse::try_new(&M, &C, &K, 2, &[0], 1.0, tight, cx),
            Err(HarmonicResponseError::Budget { .. })
        ));
        let p =
            PreparedSupportedHarmonicResponse::try_new(&M, &C, &K, 2, &[0], 1.0, controls(), cx)
                .unwrap();
        assert!(p.solve(&[C64::ZERO], &[C64::ZERO], cx).is_err());
        assert!(p.solve(&[C64::ZERO; 2], &[], cx).is_err());
        assert!(
            p.solve(&[C64::ZERO; 2], &[C64::from_re(f64::NAN)], cx)
                .is_err()
        );
        let expected = p.solve(&[C64::ONE; 2], &[C64::ONE], cx).unwrap();
        let repeated = p.solve(&[C64::ONE; 2], &[C64::ONE], cx).unwrap();
        assert_eq!(expected.displacement, repeated.displacement);
        gate.request();
        assert!(matches!(
            p.solve(&[C64::ONE; 2], &[C64::ONE], cx),
            Err(HarmonicResponseError::Cancelled)
        ));
        assert!(matches!(
            PreparedSupportedHarmonicResponse::try_new(&M, &C, &K, 2, &[0], 1.0, controls(), cx),
            Err(HarmonicResponseError::Cancelled)
        ));
    });
}
