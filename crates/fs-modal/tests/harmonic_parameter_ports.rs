//! Complete physical parameter directions through harmonic attachments and loads.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicFrequencyControls, HarmonicParameterControls, HarmonicPortControls,
    HarmonicPortParameterDirection, HarmonicResponseControls, HarmonicResponseError,
    PreparedHarmonicResponse,
};
use fs_modal::quadratic_modes::QuadraticParameterDirection;

fn with_cx(gate: &CancelGate, run: impl FnOnce(&Cx<'_>)) {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 83,
                kernel_id: 31,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 2e-12,
    }
}
fn ports() -> HarmonicPortControls {
    HarmonicPortControls {
        maximum_ports: 4,
        maximum_projection_terms: 1_000_000,
        maximum_payload_bytes: 1 << 20,
    }
}
fn tangent() -> HarmonicParameterControls {
    HarmonicParameterControls {
        maximum_additional_payload_bytes: 1 << 20,
        maximum_work: 1_000_000,
    }
}
fn relative(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-260),
        "{actual:?} != {expected:?}"
    );
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:?} != {expected:?}"
    );
}
fn power(actual: f64, expected: f64, tolerance: f64) {
    near(C64::from_re(actual), C64::from_re(expected), tolerance);
}
fn inner(a: C64, b: C64) -> f64 {
    (a.conj() * b).re
}

const M: [f64; 4] = [2.0, 0.3, 0.3, 1.0];
const C: [f64; 4] = [0.7, -0.2, 0.4, 0.5];
const K: [f64; 4] = [8.0, -1.0, -1.0, 5.0];
const B: [f64; 4] = [1.0, 0.2, -0.3, 0.8];
const Z: [C64; 4] = [
    C64::new(0.8, 0.4),
    C64::new(0.12, 0.2),
    C64::new(-0.08, 0.1),
    C64::new(0.6, -0.3),
];
const DM: [f64; 4] = [0.2, -0.03, -0.03, -0.1];
const DC: [f64; 4] = [0.1, -0.2, 0.3, -0.05];
const DK: [f64; 4] = [0.6, 0.1, 0.1, -0.2];
const DB: [f64; 4] = [0.1, -0.03, 0.05, 0.2];
const DZ: [C64; 4] = [
    C64::new(0.1, -0.2),
    C64::new(0.04, 0.05),
    C64::new(-0.07, 0.03),
    C64::new(-0.1, 0.08),
];
const F: [C64; 2] = [C64::new(1.2, -0.3), C64::new(-0.4, 0.2)];
const DF: [C64; 2] = [C64::new(-0.1, 0.04), C64::new(0.03, -0.07)];

fn mechanical() -> QuadraticParameterDirection<'static> {
    QuadraticParameterDirection {
        mass: &DM,
        damping: &DC,
        stiffness: &DK,
    }
}
fn force(theta: f64) -> [C64; 2] {
    std::array::from_fn(|i| F[i] + DF[i].scale(theta))
}
fn fixture(theta: f64, w: f64, scales: [f64; 2], cx: &Cx<'_>) -> PreparedHarmonicResponse {
    let m: [f64; 4] = std::array::from_fn(|i| M[i] + theta * DM[i]);
    let c: [f64; 4] = std::array::from_fn(|i| C[i] + theta * DC[i]);
    let k: [f64; 4] = std::array::from_fn(|i| K[i] + theta * DK[i]);
    let b: [f64; 4] = std::array::from_fn(|i| (B[i] + theta * DB[i]) * scales[i % 2]);
    let z: [C64; 4] = std::array::from_fn(|i| {
        (Z[i] + DZ[i].scale(theta))
            .scale(1.0 / scales[i / 2])
            .scale(1.0 / scales[i % 2])
    });
    PreparedHarmonicResponse::try_with_port_impedance(
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
    .unwrap()
}

#[test]
fn scalar_total_direction_matches_closed_form_including_dc_kinematics() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let (m, c, k, b) = (2.0, 0.4, 7.0, 1.3);
        let (dm, dc, dk, db) = (0.12, -0.08, 0.3, -0.15);
        let (z, dz) = (C64::new(0.6, -0.2), C64::new(0.07, 0.11));
        let (f, df) = (C64::new(0.8, -0.5), C64::new(-0.03, 0.04));
        for w in [0.0, 0.8, 1.87, 2.3] {
            let prepared = PreparedHarmonicResponse::try_with_port_impedance(
                &[m],
                &[c],
                &[k],
                1,
                &[b],
                1,
                &[z],
                w,
                controls(),
                ports(),
                cx,
            )
            .unwrap();
            for dw in [0.0, 0.7, -0.4] {
                let result = prepared
                    .solve_with_total_parameter_derivative(
                        &[f],
                        QuadraticParameterDirection {
                            mass: &[dm],
                            damping: &[dc],
                            stiffness: &[dk],
                        },
                        Some(&[df]),
                        Some(HarmonicPortParameterDirection {
                            force_map: &[db],
                            impedance: &[dz],
                        }),
                        dw,
                        tangent(),
                        cx,
                    )
                    .unwrap();
                let iw = C64::new(0.0, -w);
                let a = C64::from_re(k - w * w * m) + iw * (C64::from_re(c) + z.scale(b * b));
                let da = C64::from_re(dk - 2.0 * w * dw * m - w * w * dm)
                    + C64::new(0.0, -dw) * (C64::from_re(c) + z.scale(b * b))
                    + iw * (C64::from_re(dc) + z.scale(2.0 * b * db) + dz.scale(b * b));
                let q = f / a;
                let dq = (df - da * q) / a;
                let v = iw * q;
                let dv = iw * dq + C64::new(0.0, -dw) * q;
                let y = v.scale(b);
                let dy = dv.scale(b) + v.scale(db);
                let e = z * y;
                let de = z * dy + dz * y;
                relative(result.displacement_derivative[0], dq, 3e-12);
                relative(result.velocity_derivative[0], dv, 3e-12);
                let port = result.port_load.unwrap();
                relative(port.flow_derivative[0], dy, 3e-12);
                relative(port.effort_derivative[0], de, 3e-12);
                power(
                    result.average_input_power_derivative,
                    0.5 * (inner(df, v) + inner(f, dv)),
                    3e-12,
                );
                power(
                    result.average_viscous_power_derivative,
                    c * inner(v, dv) + 0.5 * dc * inner(v, v),
                    3e-12,
                );
                power(
                    port.average_absorbed_power_derivative,
                    0.5 * (inner(dy, e) + inner(y, de)),
                    3e-12,
                );
                power(result.power_balance_residual_derivative, 0.0, 3e-12);
                assert!(result.relative_derivative_residual < 2e-12);
            }
        }
    });
}

#[test]
fn coupled_nonreciprocal_total_direction_matches_independent_physical_perturbations() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let prepared = fixture(0.0, 1.7, [1.0; 2], cx);
        let result = prepared
            .solve_with_total_parameter_derivative(
                &F,
                mechanical(),
                Some(&DF),
                Some(HarmonicPortParameterDirection {
                    force_map: &DB,
                    impedance: &DZ,
                }),
                0.13,
                tangent(),
                cx,
            )
            .unwrap();
        let h = 1e-5;
        let plus = fixture(h, 1.7 + 0.13 * h, [1.0; 2], cx)
            .solve(&force(h), cx)
            .unwrap();
        let minus = fixture(-h, 1.7 - 0.13 * h, [1.0; 2], cx)
            .solve(&force(-h), cx)
            .unwrap();
        for i in 0..2 {
            near(
                result.displacement_derivative[i],
                (plus.displacement[i] - minus.displacement[i]).scale(0.5 / h),
                3e-8,
            );
            near(
                result.velocity_derivative[i],
                (plus.velocity[i] - minus.velocity[i]).scale(0.5 / h),
                3e-8,
            );
        }
        power(
            result.average_input_power_derivative,
            (plus.average_input_power_w - minus.average_input_power_w) / (2.0 * h),
            3e-8,
        );
        power(
            result.average_viscous_power_derivative,
            (plus.average_viscous_power_w - minus.average_viscous_power_w) / (2.0 * h),
            3e-8,
        );
        let (actual, pa, pb) = (
            result.port_load.unwrap(),
            plus.port_load.unwrap(),
            minus.port_load.unwrap(),
        );
        for a in 0..2 {
            near(
                actual.flow_derivative[a],
                (pa.flow[a] - pb.flow[a]).scale(0.5 / h),
                3e-8,
            );
            near(
                actual.effort_derivative[a],
                (pa.effort[a] - pb.effort[a]).scale(0.5 / h),
                3e-8,
            );
        }
        power(
            actual.average_absorbed_power_derivative,
            (pa.average_absorbed_power_w - pb.average_absorbed_power_w) / (2.0 * h),
            3e-8,
        );
        power(result.power_balance_residual_derivative, 0.0, 3e-12);
    });
}

#[test]
fn compact_loaded_bare_poles_keep_the_complete_attachment_and_frequency_direction() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let natural = [1.0, 1.0, 2.0];
        let rates = [0.0, 0.0, 0.04];
        let damping = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04];
        let map = [1.0, 0.3, -0.2, 1.1, 0.4, -0.5];
        let dmap = [0.1, -0.06, 0.03, 0.12, -0.05, 0.07];
        let dz = [
            C64::new(0.05, -0.03),
            C64::new(-0.01, 0.02),
            C64::new(0.04, -0.01),
            C64::new(-0.02, 0.06),
        ];
        let force = [C64::new(0.8, -0.1), C64::new(-0.3, 0.7), C64::new(0.2, 0.4)];
        let direction = QuadraticParameterDirection {
            mass: &[0.01, 0.0, 0.0, 0.0, -0.02, 0.0, 0.0, 0.0, 0.0],
            damping: &[0.0, -0.03, 0.0, 0.03, 0.0, 0.0, 0.0, 0.0, 0.02],
            stiffness: &[0.1, -0.02, 0.0, -0.02, -0.05, 0.0, 0.0, 0.0, 0.03],
        };
        let load_direction = HarmonicPortParameterDirection {
            force_map: &dmap,
            impedance: &dz,
        };
        let compact = PreparedHarmonicResponse::try_diagonal_modal_with_port_impedance(
            &natural,
            &rates,
            &map,
            2,
            &Z,
            1.0,
            controls(),
            ports(),
            cx,
        )
        .unwrap()
        .solve_with_total_parameter_derivative(
            &force,
            direction,
            None,
            Some(load_direction),
            0.2,
            tangent(),
            cx,
        )
        .unwrap();
        let dense = PreparedHarmonicResponse::try_modal_with_port_impedance(
            &natural,
            &damping,
            &map,
            2,
            &Z,
            1.0,
            controls(),
            ports(),
            cx,
        )
        .unwrap()
        .solve_with_total_parameter_derivative(
            &force,
            direction,
            None,
            Some(load_direction),
            0.2,
            tangent(),
            cx,
        )
        .unwrap();
        for i in 0..3 {
            relative(
                compact.displacement_derivative[i],
                dense.displacement_derivative[i],
                3e-10,
            );
            relative(
                compact.velocity_derivative[i],
                dense.velocity_derivative[i],
                3e-10,
            );
        }
        let (a, b) = (compact.port_load.unwrap(), dense.port_load.unwrap());
        for p in 0..2 {
            relative(a.flow_derivative[p], b.flow_derivative[p], 3e-10);
            relative(a.effort_derivative[p], b.effort_derivative[p], 3e-10);
        }
        power(
            a.average_absorbed_power_derivative,
            b.average_absorbed_power_derivative,
            3e-10,
        );
        power(compact.power_balance_residual_derivative, 0.0, 3e-10);
        assert!(compact.relative_derivative_residual < 2e-12);
    });
}

#[test]
fn independent_port_unit_scales_preserve_total_motion_and_work() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let run = |scales: [f64; 2]| {
            let db: [f64; 4] = std::array::from_fn(|i| DB[i] * scales[i % 2]);
            let dz: [C64; 4] = std::array::from_fn(|i| {
                DZ[i].scale(1.0 / scales[i / 2]).scale(1.0 / scales[i % 2])
            });
            fixture(0.0, 1.7, scales, cx)
                .solve_with_total_parameter_derivative(
                    &F,
                    mechanical(),
                    Some(&DF),
                    Some(HarmonicPortParameterDirection {
                        force_map: &db,
                        impedance: &dz,
                    }),
                    0.13,
                    tangent(),
                    cx,
                )
                .unwrap()
        };
        let reference = run([1.0; 2]);
        for scales in [[1e100, 1e-100], [1e-100, 1e100]] {
            let result = run(scales);
            for i in 0..2 {
                relative(
                    result.displacement_derivative[i],
                    reference.displacement_derivative[i],
                    2e-11,
                );
                relative(
                    result.velocity_derivative[i],
                    reference.velocity_derivative[i],
                    2e-11,
                );
            }
            let (actual, expected) = (
                result.port_load.unwrap(),
                reference.port_load.as_ref().unwrap(),
            );
            for a in 0..2 {
                relative(
                    actual.flow_derivative[a].scale(1.0 / scales[a]),
                    expected.flow_derivative[a],
                    2e-11,
                );
                relative(
                    actual.effort_derivative[a].scale(scales[a]),
                    expected.effort_derivative[a],
                    2e-11,
                );
            }
            power(
                actual.average_absorbed_power_derivative,
                expected.average_absorbed_power_derivative,
                2e-11,
            );
            power(
                result.average_input_power_derivative,
                reference.average_input_power_derivative,
                2e-11,
            );
            power(result.power_balance_residual_derivative, 0.0, 2e-11);
        }
    });
}

#[test]
fn changing_only_port_units_changes_port_variables_without_changing_the_physics() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let z = C64::new(0.6, -0.2);
        let prepared = PreparedHarmonicResponse::try_with_port_impedance(
            &[2.0],
            &[0.4],
            &[7.0],
            1,
            &[2.0],
            1,
            &[z],
            1.1,
            controls(),
            ports(),
            cx,
        )
        .unwrap();
        // B(theta)=B*(1+theta/8), Z(theta)=Z/(1+theta/8)^2
        // leaves B Z B^T unchanged to every order.
        let result = prepared
            .solve_with_parameter_and_port_derivative(
                &[C64::new(0.7, -0.3)],
                QuadraticParameterDirection {
                    mass: &[0.0],
                    damping: &[0.0],
                    stiffness: &[0.0],
                },
                None,
                HarmonicPortParameterDirection {
                    force_map: &[0.25],
                    impedance: &[z.scale(-0.25)],
                },
                tangent(),
                cx,
            )
            .unwrap();
        near(result.displacement_derivative[0], C64::ZERO, 2e-12);
        near(result.velocity_derivative[0], C64::ZERO, 2e-12);
        let primal = result.response.port_load.unwrap();
        let slope = result.port_load.unwrap();
        relative(slope.flow_derivative[0], primal.flow[0].scale(0.125), 2e-12);
        relative(
            slope.effort_derivative[0],
            primal.effort[0].scale(-0.125),
            2e-12,
        );
        power(slope.average_absorbed_power_derivative, 0.0, 2e-12);
        power(result.average_input_power_derivative, 0.0, 2e-12);
    });
}

#[test]
fn huge_attachment_products_cancel_before_their_physical_magnitude_is_restored() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let prepared = PreparedHarmonicResponse::try_with_port_impedance(
            &[1.0, 0.0, 0.0, 1.0],
            &[0.0; 4],
            &[2.0, 0.0, 0.0, 2.0],
            2,
            &[1.0, -1.0],
            1,
            &[C64::new(0.0, 1e-200)],
            1.0,
            controls(),
            ports(),
            cx,
        )
        .unwrap();
        let result = prepared
            .solve_with_parameter_and_port_derivative(
                &[C64::from_re(1e200), C64::from_re(-1e200)],
                QuadraticParameterDirection {
                    mass: &[0.0; 4],
                    damping: &[0.0; 4],
                    stiffness: &[0.0; 4],
                },
                None,
                HarmonicPortParameterDirection {
                    force_map: &[1e200, 1e200],
                    impedance: &[C64::ZERO],
                },
                tangent(),
                cx,
            )
            .unwrap();
        // dB^T v consists of +/-i*1e400, but is exactly zero. The remaining
        // force direction is -dB*effort=-2e200*(1,1), in B^T's nullspace.
        for &dq in &result.displacement_derivative {
            relative(dq, C64::from_re(-2e200), 2e-13);
        }
        let port = result.port_load.unwrap();
        assert_eq!(port.flow_derivative, [C64::ZERO]);
        assert_eq!(port.effort_derivative, [C64::ZERO]);
        assert_eq!(port.average_absorbed_power_derivative, 0.0);
        assert_eq!(result.average_input_power_derivative, 0.0);
        assert_eq!(result.average_viscous_power_derivative, 0.0);
    });
}

#[test]
fn total_frequency_direction_matches_the_existing_frequency_owner() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for kind in 0..3 {
            for w in [0.0, 0.4, 2.3] {
                let dz = if kind == 2 { DZ } else { [C64::ZERO; 4] };
                let z: [C64; 4] = std::array::from_fn(|i| Z[i] + dz[i].scale(w));
                let prepared = if kind == 0 {
                    PreparedHarmonicResponse::try_new(&M, &C, &K, 2, w, controls(), cx)
                } else {
                    PreparedHarmonicResponse::try_with_port_impedance(
                        &M,
                        &C,
                        &K,
                        2,
                        &B,
                        2,
                        &z,
                        w,
                        controls(),
                        ports(),
                        cx,
                    )
                }
                .unwrap();
                let total = prepared
                    .solve_with_total_parameter_derivative(
                        &F,
                        QuadraticParameterDirection {
                            mass: &[0.0; 4],
                            damping: &[0.0; 4],
                            stiffness: &[0.0; 4],
                        },
                        Some(&DF),
                        (kind == 2).then_some(HarmonicPortParameterDirection {
                            force_map: &[0.0; 4],
                            impedance: &dz,
                        }),
                        1.0,
                        tangent(),
                        cx,
                    )
                    .unwrap();
                let frequency = prepared
                    .solve_with_frequency_derivative(
                        &F,
                        Some(&DF),
                        (kind != 0).then_some(&dz[..]),
                        HarmonicFrequencyControls {
                            maximum_additional_payload_bytes: 1 << 20,
                            maximum_work: 1_000_000,
                        },
                        cx,
                    )
                    .unwrap();
                for i in 0..2 {
                    relative(
                        total.displacement_derivative[i],
                        frequency.displacement_derivative[i],
                        5e-12,
                    );
                    relative(
                        total.velocity_derivative[i],
                        frequency.velocity_derivative[i],
                        5e-12,
                    );
                }
                power(
                    total.average_input_power_derivative,
                    frequency.average_input_power_derivative,
                    5e-12,
                );
                power(
                    total.average_viscous_power_derivative,
                    frequency.average_viscous_power_derivative,
                    5e-12,
                );
                if let (Some(a), Some(b)) = (total.port_load, frequency.port_load) {
                    for p in 0..2 {
                        relative(a.flow_derivative[p], b.flow_derivative[p], 5e-12);
                        relative(a.effort_derivative[p], b.effort_derivative[p], 5e-12);
                    }
                    power(
                        a.average_absorbed_power_derivative,
                        b.average_absorbed_power_derivative,
                        5e-12,
                    );
                }
            }
        }
    });
}

#[test]
fn load_and_total_direction_refusals_preserve_the_prepared_factor() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let prepared = fixture(0.0, 1.7, [1.0; 2], cx);
        let original = prepared.solve(&F, cx).unwrap();
        for db in [
            &[][..],
            &[0.0; 3][..],
            &[f64::NAN; 4][..],
            &[f64::INFINITY; 4][..],
        ] {
            assert!(
                prepared
                    .solve_with_parameter_and_port_derivative(
                        &F,
                        mechanical(),
                        Some(&DF),
                        HarmonicPortParameterDirection {
                            force_map: db,
                            impedance: &DZ,
                        },
                        tangent(),
                        cx,
                    )
                    .is_err()
            );
        }
        for dz in [
            &[][..],
            &[C64::ZERO; 3][..],
            &[C64::new(f64::NAN, 0.0); 4][..],
        ] {
            assert!(
                prepared
                    .solve_with_parameter_and_port_derivative(
                        &F,
                        mechanical(),
                        Some(&DF),
                        HarmonicPortParameterDirection {
                            force_map: &DB,
                            impedance: dz,
                        },
                        tangent(),
                        cx,
                    )
                    .is_err()
            );
        }
        for dw in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(
                prepared
                    .solve_with_total_parameter_derivative(
                        &F,
                        mechanical(),
                        Some(&DF),
                        None,
                        dw,
                        tangent(),
                        cx,
                    )
                    .is_err()
            );
        }
        for limits in [
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
                prepared
                    .solve_with_parameter_and_port_derivative(
                        &F,
                        mechanical(),
                        Some(&DF),
                        HarmonicPortParameterDirection {
                            force_map: &DB,
                            impedance: &DZ,
                        },
                        limits,
                        cx,
                    )
                    .is_err()
            );
        }
        let unloaded =
            PreparedHarmonicResponse::try_new(&M, &C, &K, 2, 1.7, controls(), cx).unwrap();
        assert!(
            unloaded
                .solve_with_parameter_and_port_derivative(
                    &F,
                    mechanical(),
                    Some(&DF),
                    HarmonicPortParameterDirection {
                        force_map: &DB,
                        impedance: &DZ,
                    },
                    tangent(),
                    cx,
                )
                .is_err()
        );
        let fixed = prepared
            .solve_with_parameter_derivative(&F, mechanical(), Some(&DF), tangent(), cx)
            .unwrap();
        let explicit = prepared
            .solve_with_parameter_and_port_derivative(
                &F,
                mechanical(),
                Some(&DF),
                HarmonicPortParameterDirection {
                    force_map: &[0.0; 4],
                    impedance: &[C64::ZERO; 4],
                },
                tangent(),
                cx,
            )
            .unwrap();
        assert_eq!(
            fixed.displacement_derivative,
            explicit.displacement_derivative
        );
        assert_eq!(fixed.velocity_derivative, explicit.velocity_derivative);
        assert_eq!(
            fixed.average_input_power_derivative,
            explicit.average_input_power_derivative
        );
        assert_eq!(
            prepared.solve(&F, cx).unwrap().displacement,
            original.displacement
        );
        gate.request();
        assert!(matches!(
            prepared.solve_with_total_parameter_derivative(
                &F,
                mechanical(),
                None,
                None,
                1.0,
                tangent(),
                cx,
            ),
            Err(HarmonicResponseError::Cancelled)
        ));
    });
}
