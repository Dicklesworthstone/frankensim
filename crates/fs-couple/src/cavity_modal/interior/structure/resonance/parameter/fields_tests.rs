//! G1/G3: original coupled rows, pressure kinematics and physical work jets.
//! This exercises recovery from the actual canonical harmonic owner; FE geometry
//! and structural/acoustic quadrature assembly are separate integration coverage.

use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_modal::harmonic::{
    HarmonicFrequencyParameterControls, HarmonicResponseControls, PreparedHarmonicResponse,
};
use fs_modal::quadratic_modes::QuadraticParameterDirection;

type Jet = [Vec<C64>; 5];

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 61,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    });
}

struct Pencil {
    m: Vec<f64>,
    c: Vec<f64>,
    k: Vec<f64>,
    dm: Vec<f64>,
    dc: Vec<f64>,
    dk: Vec<f64>,
    b: Vec<f64>,
    db: Vec<f64>,
}

fn pencil() -> Pencil {
    Pencil {
        m: vec![
            1.25, 0., 0., 0., 0., 0.875, 0., 0., 0., 0., 0.75, 0., 0., 0., 0., 1.,
        ],
        c: vec![
            0.375, 0.125, -0.5, -0.125, -0.125, 0.25, 0.25, -0.75, 0.5, -0.25, 0.25, 0.0625, 0.125,
            0.75, 0.0625, 0.375,
        ],
        k: vec![
            5., -0.5, 0., 0., -0.5, 6., 0., 0., 0., 0., 3.25, -0.25, 0., 0., -0.25, 4.5,
        ],
        dm: vec![
            0.125, 0., 0., 0., 0., -0.0625, 0., 0., 0., 0., 0.03125, 0., 0., 0., 0., 0.0625,
        ],
        dc: vec![
            0.0625, 0.09375, -0.0625, 0.015625, -0.03125, -0.03125, -0.03125, -0.046875, 0.0625,
            0.03125, 0.03125, -0.015625, -0.015625, 0.046875, -0.015625, 0.0625,
        ],
        dk: vec![
            0.25, -0.0625, 0., 0., -0.0625, -0.125, 0., 0., 0., 0., 0.125, 0.03125, 0., 0.,
            0.03125, 0.0625,
        ],
        b: vec![0.5, -0.25, 0.125, 0.75],
        db: vec![0.0625, 0.03125, -0.015625, 0.046875],
    }
}

fn sources() -> Jet {
    [
        vec![
            C64::new(0.8, -0.3),
            C64::new(-0.2, 0.6),
            C64::new(0.15, 0.07),
            C64::new(-0.08, 0.11),
        ],
        vec![
            C64::new(-0.1, 0.05),
            C64::new(0.03, -0.08),
            C64::new(0.01, 0.02),
            C64::new(-0.02, 0.03),
        ],
        vec![
            C64::new(0.02, -0.04),
            C64::new(-0.03, 0.02),
            C64::new(0.005, -0.01),
            C64::new(0.01, 0.005),
        ],
        vec![
            C64::new(0.07, 0.09),
            C64::new(-0.04, 0.06),
            C64::new(0.02, -0.03),
            C64::new(-0.01, 0.015),
        ],
        vec![
            C64::new(-0.06, 0.02),
            C64::new(0.08, -0.01),
            C64::new(0.015, 0.01),
            C64::new(0.005, -0.02),
        ],
    ]
}

fn native(p: &Pencil, w: f64, f: &Jet, cx: &Cx<'_>) -> HarmonicFrequencyParameterResponse {
    PreparedHarmonicResponse::try_new(
        &p.m,
        &p.c,
        &p.k,
        4,
        w,
        HarmonicResponseControls {
            maximum_coordinates: 4,
            maximum_payload_bytes: 1 << 20,
            maximum_relative_residual: 2e-11,
        },
        cx,
    )
    .unwrap()
    .solve_with_frequency_parameter_derivative(
        &f[0],
        QuadraticParameterDirection {
            mass: &p.dm,
            damping: &p.dc,
            stiffness: &p.dk,
        },
        Some(&f[1]),
        Some(&f[2]),
        Some(&f[3]),
        Some(&f[4]),
        HarmonicFrequencyParameterControls {
            maximum_additional_payload_bytes: 1 << 20,
            maximum_work: 1_000_000,
        },
        cx,
    )
    .unwrap()
}

fn recovered(p: &Pencil, w: f64, f: &Jet, cx: &Cx<'_>) -> Recovered {
    recover(
        2,
        &p.b,
        &p.c,
        &p.dc,
        &p.db,
        field(f),
        native(p, w, f, cx),
        &mut || false,
    )
    .unwrap()
}

fn close(a: C64, b: C64) {
    assert!(
        (a - b).abs() <= 4e-10 * a.abs().max(b.abs()).max(1.0),
        "{a:?} != {b:?}"
    );
}

fn scalars(p: HarmonicPowerParameterJet) -> [f64; 5] {
    [
        p.sample.average_power_w,
        p.sample.angular_derivative,
        p.angular_second_derivative,
        p.parameter_derivative,
        p.angular_parameter_derivative,
    ]
}

fn dot(a: &[C64], b: &[C64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| (a.conj() * *b).re).sum()
}

// Independent ordinary-range product rules, deliberately outside the production
// normalized work helper. Both phasors and both mixed cross terms are retained.
fn power(a: &Jet, b: &Jet) -> [f64; 5] {
    [
        0.5 * dot(&a[0], &b[0]),
        0.5 * (dot(&a[1], &b[0]) + dot(&a[0], &b[1])),
        0.5 * (dot(&a[2], &b[0]) + 2. * dot(&a[1], &b[1]) + dot(&a[0], &b[2])),
        0.5 * (dot(&a[3], &b[0]) + dot(&a[0], &b[3])),
        0.5 * (dot(&a[4], &b[0]) + dot(&a[1], &b[3]) + dot(&a[3], &b[1]) + dot(&a[0], &b[4])),
    ]
}

fn compare_power(actual: HarmonicPowerParameterJet, expected: [f64; 5]) {
    for (a, b) in scalars(actual).into_iter().zip(expected) {
        close(C64::from_re(a), C64::from_re(b));
    }
}

fn states(r: &Recovered) -> (Jet, Jet) {
    let derivatives = [
        &r.angular_derivative,
        &r.angular_second_derivative,
        &r.parameter_derivative,
        &r.angular_parameter_derivative,
    ];
    let q = core::array::from_fn(|order| {
        if order == 0 {
            return r.response.displacement.clone();
        }
        let d = derivatives[order - 1];
        d.structural_displacement
            .iter()
            .chain(&d.pressure_potential_pa_s)
            .copied()
            .collect()
    });
    let v = core::array::from_fn(|order| {
        if order == 0 {
            return r.response.velocity.clone();
        }
        let d = derivatives[order - 1];
        d.structural_velocity
            .iter()
            .chain(&d.pressure_pa)
            .copied()
            .collect()
    });
    (q, v)
}

fn row(a: &[f64], x: &[C64], i: usize) -> C64 {
    (0..4).fold(C64::ZERO, |s, j| s + x[j].scale(a[4 * i + j]))
}

fn part(x: &Jet, range: core::ops::Range<usize>) -> Jet {
    core::array::from_fn(|order| x[order][range.clone()].to_vec())
}

fn heat_reference(p: &Pencil, velocity: &Jet, offset: usize) -> [f64; 5] {
    let mut effort: Jet = core::array::from_fn(|_| vec![C64::ZERO; 2]);
    for i in 0..2 {
        for j in 0..2 {
            let c = 0.5 * (p.c[4 * (offset + i) + offset + j] + p.c[4 * (offset + j) + offset + i]);
            let dc =
                0.5 * (p.dc[4 * (offset + i) + offset + j] + p.dc[4 * (offset + j) + offset + i]);
            for order in 0..5 {
                effort[order][i] = effort[order][i] + velocity[order][j].scale(c);
            }
            effort[3][i] = effort[3][i] + velocity[0][j].scale(dc);
            effort[4][i] = effort[4][i] + velocity[1][j].scale(dc);
        }
    }
    power(velocity, &effort)
}

#[test]
fn complete_original_rows_pressure_kinematics_interface_changes_and_work_match() {
    with_cx(|cx| {
        let p = pencil();
        let f = sources();
        for w in [0.0, 1.4, 2.1, 2.7] {
            let r = recovered(&p, w, &f, cx);
            let (q, v) = states(&r);
            for i in 0..4 {
                // Pressure is the acoustic velocity of the potential, including
                // the nonzero DC slope and every higher kinematic product rule.
                for (order, expected) in [
                    (0, C64::new(0., -w) * q[0][i]),
                    (1, C64::new(0., -1.) * (q[0][i] + q[1][i].scale(w))),
                    (
                        2,
                        C64::new(0., -1.) * (q[1][i].scale(2.) + q[2][i].scale(w)),
                    ),
                    (3, C64::new(0., -w) * q[3][i]),
                    (4, C64::new(0., -1.) * (q[3][i] + q[4][i].scale(w))),
                ] {
                    close(v[order][i], expected);
                }
                let dynamic = |order: usize| {
                    row(&p.k, &q[order], i) + row(&p.c, &v[order], i)
                        - row(&p.m, &q[order], i).scale(w * w)
                };
                let direction = |order: usize| {
                    row(&p.dk, &q[order], i) + row(&p.dc, &v[order], i)
                        - row(&p.dm, &q[order], i).scale(w * w)
                };
                close(dynamic(0), f[0][i]);
                close(dynamic(1) - row(&p.m, &q[0], i).scale(2. * w), f[1][i]);
                close(
                    dynamic(2) - row(&p.m, &q[1], i).scale(4. * w) - row(&p.m, &q[0], i).scale(2.),
                    f[2][i],
                );
                close(dynamic(3) + direction(0), f[3][i]);
                close(
                    dynamic(4) - row(&p.m, &q[3], i).scale(2. * w) + direction(1)
                        - row(&p.dm, &q[0], i).scale(2. * w),
                    f[4][i],
                );
            }
            let sv = part(&v, 0..2);
            let pressure = part(&v, 2..4);
            let derivatives = [
                &r.angular_derivative,
                &r.angular_second_derivative,
                &r.parameter_derivative,
                &r.angular_parameter_derivative,
            ];
            let force: Jet = core::array::from_fn(|o| {
                if o == 0 {
                    r.fluid_force.clone()
                } else {
                    derivatives[o - 1].fluid_force.clone()
                }
            });
            let flow: Jet = core::array::from_fn(|o| {
                if o == 0 {
                    r.displaced_flow.clone()
                } else {
                    derivatives[o - 1].displaced_flow_m3_s.clone()
                }
            });
            for order in 0..5 {
                for i in 0..2 {
                    let base_order = if order == 3 {
                        Some(0)
                    } else if order == 4 {
                        Some(1)
                    } else {
                        None
                    };
                    let mut expected_force = C64::ZERO;
                    let mut expected_flow = C64::ZERO;
                    for j in 0..2 {
                        expected_force = expected_force + pressure[order][j].scale(p.b[2 * j + i]);
                        expected_flow = expected_flow + sv[order][j].scale(p.b[2 * i + j]);
                        if let Some(o) = base_order {
                            expected_force = expected_force + pressure[o][j].scale(p.db[2 * j + i]);
                            expected_flow = expected_flow + sv[o][j].scale(p.db[2 * i + j]);
                        }
                    }
                    close(force[order][i], expected_force);
                    close(flow[order][i], expected_flow);
                }
            }
            compare_power(r.power.structural_input, power(&sv, &part(&f, 0..2)));
            compare_power(r.power.acoustic_input, power(&part(&f, 2..4), &pressure));
            compare_power(r.power.fluid_to_structure, power(&sv, &force));
            compare_power(r.power.fluid_to_structure, power(&flow, &pressure));
            compare_power(r.power.structural_dissipation, heat_reference(&p, &sv, 0));
            compare_power(r.power.wall_dissipation, heat_reference(&p, &pressure, 2));
            for balance in [
                r.power.structural_balance().unwrap(),
                r.power.acoustic_balance().unwrap(),
                r.power.original_balance,
            ] {
                compare_power(balance, [0.; 5]);
            }
            for residual in [
                r.relative_derivative_residual,
                r.relative_second_derivative_residual,
                r.relative_parameter_derivative_residual,
                r.relative_mixed_derivative_residual,
            ] {
                assert!(residual < 2e-11);
            }
        }
    });
}

#[test]
fn physical_frequency_dilation_moves_all_fields_and_work_with_their_actual_units() {
    with_cx(|cx| {
        let mut p = pencil();
        p.dm.fill(0.);
        p.dc = p.c.clone();
        p.dk = p.k.iter().map(|x| 2. * x).collect();
        p.db = p.b.clone();
        for w in [0., 1.6, 2.2] {
            let mut f = sources();
            for i in 0..4 {
                f[3][i] = f[0][i].scale(2.) - f[1][i].scale(w);
                f[4][i] = f[1][i] - f[2][i].scale(w);
            }
            let r = recovered(&p, w, &f, cx);
            let (q, v) = states(&r);
            for i in 0..4 {
                close(q[3][i], q[1][i].scale(-w));
                close(q[4][i], -q[1][i] - q[2][i].scale(w));
                close(v[3][i], v[0][i] - v[1][i].scale(w));
                close(v[4][i], v[2][i].scale(-w));
            }
            for i in 0..2 {
                close(
                    r.parameter_derivative.fluid_force[i],
                    r.fluid_force[i].scale(2.) - r.angular_derivative.fluid_force[i].scale(w),
                );
                close(
                    r.angular_parameter_derivative.fluid_force[i],
                    r.angular_derivative.fluid_force[i]
                        - r.angular_second_derivative.fluid_force[i].scale(w),
                );
                close(
                    r.parameter_derivative.displaced_flow_m3_s[i],
                    r.displaced_flow[i].scale(2.)
                        - r.angular_derivative.displaced_flow_m3_s[i].scale(w),
                );
            }
            for power in [
                r.power.structural_input,
                r.power.acoustic_input,
                r.power.fluid_to_structure,
                r.power.structural_dissipation,
                r.power.wall_dissipation,
            ] {
                let p = scalars(power);
                close(C64::from_re(p[3]), C64::from_re(3. * p[0] - w * p[1]));
                close(C64::from_re(p[4]), C64::from_re(2. * p[1] - w * p[2]));
            }
        }
    });
}

#[test]
fn physical_recovery_refuses_incomplete_jets_and_cancellation_without_partial_output() {
    with_cx(|cx| {
        let p = pencil();
        let f = sources();
        let jet = native(&p, 1.6, &f, cx);
        assert!(matches!(
            recover(
                2,
                &p.b,
                &p.c,
                &p.dc,
                &p.db,
                field(&f),
                jet.clone(),
                &mut || true
            ),
            Err(HarmonicResponseError::Cancelled)
        ));
        assert!(matches!(
            recover(
                2,
                &p.b[..3],
                &p.c,
                &p.dc,
                &p.db,
                field(&f),
                jet.clone(),
                &mut || false
            ),
            Err(HarmonicResponseError::Invalid(_))
        ));
        let mut malformed = jet.clone();
        malformed.velocity_mixed_derivative.pop();
        assert!(matches!(
            recover(
                2,
                &p.b,
                &p.c,
                &p.dc,
                &p.db,
                field(&f),
                malformed,
                &mut || false
            ),
            Err(HarmonicResponseError::Invalid(_))
        ));
        let mut incomplete = f.clone();
        incomplete[4].pop();
        assert!(matches!(
            recover(
                2,
                &p.b,
                &p.c,
                &p.dc,
                &p.db,
                field(&incomplete),
                jet,
                &mut || false
            ),
            Err(HarmonicResponseError::Invalid(_))
        ));
        let mut polls = 0;
        let result = recover(
            2,
            &p.b,
            &p.c,
            &p.dc,
            &p.db,
            field(&f),
            native(&p, 1.6, &f, cx),
            &mut || {
                polls += 1;
                polls == 7
            },
        );
        assert!(matches!(result, Err(HarmonicResponseError::Cancelled)));
        assert_eq!(polls, 7);
    });
}
