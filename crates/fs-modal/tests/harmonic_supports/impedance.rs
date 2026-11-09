//! Complete original-coordinate loads with prescribed support velocities.

use super::*;
use fs_modal::harmonic::HarmonicPortControls;

const B: [f64; 4] = [1.0, 0.2, -0.4, 1.3];

fn port_controls() -> HarmonicPortControls {
    HarmonicPortControls {
        maximum_ports: 4,
        maximum_projection_terms: 10_000,
        maximum_payload_bytes: 1 << 20,
    }
}

fn impedance() -> [C64; 4] {
    // Positive Hermitian part, with deliberately nonreciprocal mutual entries.
    [
        C64::new(2.0, 0.8),
        C64::new(0.3, 0.2),
        C64::new(-0.1, 0.2),
        C64::new(1.5, -0.4),
    ]
}

fn loaded_entry(i: usize, j: usize, omega: f64, z: &[C64; 4]) -> C64 {
    let mut load = C64::ZERO;
    for a in 0..2 {
        for b in 0..2 {
            load = load + z[2 * a + b].scale(B[2 * i + a] * B[2 * j + b]);
        }
    }
    entry(K[2 * i + j], M[2 * i + j], C[2 * i + j], omega)
        + C64::new(load.im, -load.re).scale(omega)
}

#[test]
fn mutual_load_keeps_moving_support_flow_reactions_and_separate_boundary_work() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let z = impedance();
        let support = C64::new(0.03, -0.04);
        let force = [C64::new(0.07, -0.02), C64::new(0.04, 0.09)];
        for omega in [0.3, 2.8, 4.7] {
            let prepared = PreparedSupportedHarmonicResponse::try_with_port_impedance(
                &M,
                &C,
                &K,
                2,
                &[0],
                &B,
                2,
                &z,
                omega,
                controls(),
                port_controls(),
                cx,
            )
            .unwrap();
            let q1 = (force[1] - loaded_entry(1, 0, omega, &z) * support)
                / loaded_entry(1, 1, omega, &z);
            let reaction = loaded_entry(0, 0, omega, &z) * support
                + loaded_entry(0, 1, omega, &z) * q1
                - force[0];
            let velocity = [support, q1].map(|q| C64::new(q.im, -q.re).scale(omega));
            let flow: [C64; 2] =
                core::array::from_fn(|a| velocity[0].scale(B[a]) + velocity[1].scale(B[2 + a]));
            let effort: [C64; 2] =
                core::array::from_fn(|a| z[2 * a] * flow[0] + z[2 * a + 1] * flow[1]);
            let absorbed = 0.5 * (inner(flow[0], effort[0]) + inner(flow[1], effort[1]));
            let heat = 0.5
                * (C[0] * velocity[0].norm_sq()
                    + C[3] * velocity[1].norm_sq()
                    + 2.0 * C[1] * inner(velocity[0], velocity[1]));
            let result = prepared.solve(&force, &[support], cx).unwrap();
            close(result.displacement[1], q1);
            close(result.support_reaction[0], reaction);
            let load = result.port_load.as_ref().unwrap();
            for a in 0..2 {
                close(load.flow[a], flow[a]);
                close(load.effort[a], effort[a]);
            }
            scalar(load.average_absorbed_power_w, absorbed);
            scalar(result.average_viscous_power_w, heat);
            scalar(
                result.average_support_power_w,
                0.5 * inner(reaction, velocity[0]),
            );
            scalar(
                result.average_input_power_w,
                0.5 * (inner(force[0], velocity[0]) + inner(force[1], velocity[1])),
            );
            assert!(absorbed > 0.0);
            scalar(result.power_balance_residual_w, 0.0);
            let repeated = prepared.solve(&force, &[support], cx).unwrap();
            assert_eq!(result.displacement, repeated.displacement);
            assert_eq!(result.support_reaction, repeated.support_reaction);
            assert_eq!(load.flow, repeated.port_load.unwrap().flow);
        }
    });
}

#[test]
fn empty_all_and_static_support_limits_keep_complete_loaded_power() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let z = impedance();
        let omega = 1.3;
        let force = [C64::new(0.07, -0.09), C64::new(0.04, 0.02)];
        let ordinary = PreparedHarmonicResponse::try_with_port_impedance(
            &M,
            &C,
            &K,
            2,
            &B,
            2,
            &z,
            omega,
            controls(),
            port_controls(),
            cx,
        )
        .unwrap()
        .solve(&force, cx)
        .unwrap();
        let empty = PreparedSupportedHarmonicResponse::try_with_port_impedance(
            &M,
            &C,
            &K,
            2,
            &[],
            &B,
            2,
            &z,
            omega,
            controls(),
            port_controls(),
            cx,
        )
        .unwrap()
        .solve(&force, &[], cx)
        .unwrap();
        for i in 0..2 {
            close(empty.displacement[i], ordinary.displacement[i]);
        }
        scalar(
            empty.port_load.unwrap().average_absorbed_power_w,
            ordinary.port_load.unwrap().average_absorbed_power_w,
        );
        let q = [C64::new(0.02, 0.03), C64::new(-0.01, 0.04)];
        let all = PreparedSupportedHarmonicResponse::try_with_port_impedance(
            &M,
            &C,
            &K,
            2,
            &[0, 1],
            &B,
            2,
            &z,
            omega,
            controls(),
            port_controls(),
            cx,
        )
        .unwrap()
        .solve(&force, &q, cx)
        .unwrap();
        for i in 0..2 {
            close(all.displacement[i], q[i]);
            close(
                all.support_reaction[i],
                loaded_entry(i, 0, omega, &z) * q[0] + loaded_entry(i, 1, omega, &z) * q[1]
                    - force[i],
            );
        }
        assert_eq!(all.relative_residual, 0.0);
        scalar(all.power_balance_residual_w, 0.0);
        let huge = [
            C64::from_re(1e290),
            C64::ZERO,
            C64::ZERO,
            C64::from_re(2e290),
        ];
        let static_result = PreparedSupportedHarmonicResponse::try_with_port_impedance(
            &M,
            &C,
            &K,
            2,
            &[0],
            &B,
            2,
            &huge,
            0.0,
            controls(),
            port_controls(),
            cx,
        )
        .unwrap()
        .solve(&force, &[q[0]], cx)
        .unwrap();
        let static_bare =
            PreparedSupportedHarmonicResponse::try_new(&M, &C, &K, 2, &[0], 0.0, controls(), cx)
                .unwrap()
                .solve(&force, &[q[0]], cx)
                .unwrap();
        assert_eq!(static_result.displacement, static_bare.displacement);
        assert_eq!(
            static_result.port_load.unwrap().average_absorbed_power_w,
            0.0
        );
        assert_eq!(static_result.power_balance_residual_w, 0.0);
    });
}

#[test]
fn a_physical_load_resolves_bare_free_resonance_before_support_elimination() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let identity = [1.0, 0.0, 0.0, 1.0];
        let zero = [0.0; 4];
        assert!(
            PreparedSupportedHarmonicResponse::try_new(
                &identity,
                &zero,
                &identity,
                2,
                &[0],
                1.0,
                controls(),
                cx,
            )
            .is_err()
        );
        let z = C64::new(2.0, 0.4);
        let prepared = PreparedSupportedHarmonicResponse::try_with_port_impedance(
            &identity,
            &zero,
            &identity,
            2,
            &[0],
            &[0.5, 1.0],
            1,
            &[z],
            1.0,
            controls(),
            port_controls(),
            cx,
        )
        .unwrap();
        let support = C64::new(0.03, -0.02);
        let force = [C64::new(0.1, 0.0), C64::new(0.0, 0.2)];
        let load = C64::new(z.im, -z.re);
        let expected = force[1] / load - support.scale(0.5);
        let result = prepared.solve(&force, &[support], cx).unwrap();
        close(result.displacement[1], expected);
        assert_eq!(result.average_viscous_power_w, 0.0);
        assert!(result.port_load.unwrap().average_absorbed_power_w > 0.0);
        scalar(result.power_balance_residual_w, 0.0);
    });
}

#[test]
fn loaded_support_admission_and_cancellation_leave_the_original_factor_usable() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let z = impedance();
        let passive = PreparedSupportedHarmonicResponse::try_with_port_impedance(
            &M,
            &C,
            &K,
            2,
            &[0],
            &B,
            2,
            &z,
            1.0,
            controls(),
            port_controls(),
            cx,
        )
        .unwrap();
        let negative = [C64::from_re(-1.0), C64::ZERO, C64::ZERO, C64::ONE];
        assert!(
            PreparedSupportedHarmonicResponse::try_with_port_impedance(
                &M,
                &C,
                &K,
                2,
                &[0],
                &B,
                2,
                &negative,
                1.0,
                controls(),
                port_controls(),
                cx,
            )
            .is_err()
        );
        for (map, ports, limits) in [
            (&B[..3], 2, port_controls()),
            (
                &B[..],
                2,
                HarmonicPortControls {
                    maximum_payload_bytes: 0,
                    ..port_controls()
                },
            ),
            (
                &B[..],
                2,
                HarmonicPortControls {
                    maximum_projection_terms: 0,
                    ..port_controls()
                },
            ),
        ] {
            assert!(
                PreparedSupportedHarmonicResponse::try_with_port_impedance(
                    &M,
                    &C,
                    &K,
                    2,
                    &[0],
                    map,
                    ports,
                    &z,
                    1.0,
                    controls(),
                    limits,
                    cx,
                )
                .is_err()
            );
        }
        let force = [C64::new(0.1, 0.0); 2];
        let q = [C64::new(0.02, 0.01)];
        let expected = passive.solve(&force, &q, cx).unwrap();
        assert!(passive.solve(&force, &[], cx).is_err());
        let retry = passive.solve(&force, &q, cx).unwrap();
        assert_eq!(expected.displacement, retry.displacement);
        gate.request();
        assert!(matches!(
            passive.solve(&force, &q, cx),
            Err(HarmonicResponseError::Cancelled)
        ));
    });
}
