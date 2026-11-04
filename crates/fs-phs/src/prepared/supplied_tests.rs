use super::*;
use crate::{PortHamiltonian, QuadraticStorage};

fn system() -> PortHamiltonian {
    PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![0.0],
        vec![1.0],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
    )
    .unwrap()
}
fn hessian(_: &[f64], d: &[f64], out: &mut [f64]) -> bool {
    out.copy_from_slice(d);
    true
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-12, "{a} != {b}");
}

#[test]
fn driven_conductance_recovers_independent_midpoint_solution_and_both_powers() {
    let sys = system();
    for (x0, reservoir, drive) in [(0.2, 1.4, 0.03), (1.7, 0.6, -0.04), (0.5, -0.8, 0.1)] {
        let k = 0.7;
        let port = |_: &[f64], e: &[f64], out: &mut [f64]| {
            out[0] = k * (e[0] - reservoir);
            Some(NonlinearPortPower {
                supplied_power_w: -reservoir * out[0],
                dissipated_power_w: k * (e[0] - reservoir).powi(2),
            })
        };
        let tangent = |_: &[f64], _: &[f64], _: &[f64], de: &[f64], out: &mut [f64]| {
            out[0] = k * de[0];
            true
        };
        for analytic in [false, true] {
            let dt = 0.03;
            let mut work = StepWorkspace::new(&sys).unwrap();
            let (mut next, mut y) = ([0.0], [0.0]);
            let record = work
                .step_into_supplied_controlled(
                    &sys,
                    &[x0],
                    &[drive],
                    dt,
                    &mut next,
                    &mut y,
                    &port,
                    analytic.then_some(&hessian),
                    analytic.then_some(&tangent),
                    || false,
                )
                .unwrap();
            let expected =
                ((1.0 - 0.5 * dt * k) * x0 + dt * (k * reservoir + drive)) / (1.0 + 0.5 * dt * k);
            close(next[0], expected);
            let effort = 0.5 * (x0 + expected);
            close(record.dissipated, dt * k * (effort - reservoir).powi(2));
            close(
                record.supplied,
                dt * (drive * effort - reservoir * k * (effort - reservoir)),
            );
            close(record.delta_h + record.dissipated, record.supplied);
            if x0 > reservoir && reservoir > 0.0 {
                assert!(record.supplied < 0.0);
            }
        }
    }
}

#[test]
fn state_dependent_supply_reuses_analytic_newton_and_zero_supply_is_exactly_passive() {
    let sys = system();
    let mut fd = StepWorkspace::new(&sys).unwrap();
    let mut analytic = StepWorkspace::new(&sys).unwrap();
    analytic.set_backtracking_limit(4).unwrap();
    let reservoir = 1.2;
    let port = |x: &[f64], e: &[f64], out: &mut [f64]| {
        out[0] = (0.4 + x[0] * x[0]) * (e[0] - reservoir);
        Some(NonlinearPortPower {
            supplied_power_w: -reservoir * out[0],
            dissipated_power_w: (0.4 + x[0] * x[0]) * (e[0] - reservoir).powi(2),
        })
    };
    let tangent = |x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]| {
        out[0] = 2.0 * x[0] * dx[0] * (e[0] - reservoir) + (0.4 + x[0] * x[0]) * de[0];
        true
    };
    let (mut a, mut b, mut next, mut expected, mut y, mut ey) =
        ([0.2], [0.2], [0.0], [0.0], [0.0], [0.0]);
    let mut work_net = 0.0;
    for _ in 0..40 {
        let record = analytic
            .step_into_supplied_controlled(
                &sys,
                &a,
                &[0.02],
                0.01,
                &mut next,
                &mut y,
                &port,
                Some(&hessian),
                Some(&tangent),
                || false,
            )
            .unwrap();
        fd.step_into_supplied_controlled(
            &sys,
            &b,
            &[0.02],
            0.01,
            &mut expected,
            &mut ey,
            &port,
            None,
            None,
            || false,
        )
        .unwrap();
        close(next[0], expected[0]);
        close(record.balance_residual(), 0.0);
        work_net += record.supplied - record.dissipated;
        a = next;
        b = expected;
    }
    close(0.5 * (a[0] * a[0] - 0.2 * 0.2), work_net);

    let passive = |_: &[f64], e: &[f64], out: &mut [f64]| {
        out[0] = 0.7 * e[0];
        true
    };
    let zero_supply = |x: &[f64], e: &[f64], out: &mut [f64]| {
        if !passive(x, e, out) {
            return None;
        }
        power(e, out).ok().map(|loss| NonlinearPortPower {
            supplied_power_w: 0.0,
            dissipated_power_w: loss,
        })
    };
    let p = analytic
        .step_into_dissipative_controlled(
            &sys,
            &[0.4],
            &[0.03],
            0.02,
            &mut next,
            &mut y,
            &passive,
            None,
            None,
            || false,
        )
        .unwrap();
    let s = analytic
        .step_into_supplied_controlled(
            &sys,
            &[0.4],
            &[0.03],
            0.02,
            &mut expected,
            &mut ey,
            &zero_supply,
            None,
            None,
            || false,
        )
        .unwrap();
    assert_eq!(next, expected);
    assert_eq!(y, ey);
    assert_eq!(p, s);
    // A later ordinary call cannot reuse a prior supplied power.
    let ordinary = analytic
        .step_into(&sys, &[0.4], &[0.03], 0.02, &mut next, &mut y)
        .unwrap();
    close(ordinary.supplied, 0.02 * 0.03 * y[0]);
    assert_eq!(ordinary.dissipated, 0.0);
}

#[test]
fn invalid_source_loss_and_final_cancellation_leave_both_outputs_unchanged() {
    let sys = system();
    let mut work = StepWorkspace::new(&sys).unwrap();
    let (mut next, mut y) = ([31.0], [47.0]);
    for kind in 0..5 {
        let port = |_: &[f64], e: &[f64], out: &mut [f64]| {
            out[0] = -e[0];
            if kind == 0 {
                None
            } else {
                Some(NonlinearPortPower {
                    supplied_power_w: match kind {
                        1 => f64::NAN,
                        2 => 0.0,
                        _ => -1.0,
                    },
                    dissipated_power_w: if kind == 4 { -1.0 } else { 0.0 },
                })
            }
        };
        assert!(
            work.step_into_supplied_controlled(
                &sys,
                &[0.3],
                &[0.0],
                0.01,
                &mut next,
                &mut y,
                &port,
                None,
                None,
                || false
            )
            .is_err()
        );
        assert_eq!(next, [31.0]);
        assert_eq!(y, [47.0]);
        assert!(work.last_step_discrete_effort().is_none());
    }
    let port = |_: &[f64], e: &[f64], out: &mut [f64]| {
        out[0] = 0.7 * (e[0] - 1.2);
        Some(NonlinearPortPower {
            supplied_power_w: -1.2 * out[0],
            dissipated_power_w: 0.7 * (e[0] - 1.2).powi(2),
        })
    };
    let mut clean = StepWorkspace::new(&sys).unwrap();
    let (mut expected, mut ey) = ([0.0], [0.0]);
    let mut polls = 0;
    let receipt = clean
        .step_into_supplied_controlled(
            &sys,
            &[0.3],
            &[0.0],
            0.01,
            &mut expected,
            &mut ey,
            &port,
            None,
            None,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    assert_eq!(
        work.step_into_supplied_controlled(
            &sys,
            &[0.3],
            &[0.0],
            0.01,
            &mut next,
            &mut y,
            &port,
            None,
            None,
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(PreparedStepError::Cancelled)
    );
    assert_eq!(next, [31.0]);
    assert_eq!(y, [47.0]);
    let retry = work
        .step_into_supplied_controlled(
            &sys,
            &[0.3],
            &[0.0],
            0.01,
            &mut next,
            &mut y,
            &port,
            None,
            None,
            || false,
        )
        .unwrap();
    assert_eq!(next, expected);
    assert_eq!(y, ey);
    assert_eq!(retry, receipt);
}

#[test]
fn independent_positive_heat_survives_large_opposing_supply_and_storage_powers() {
    let sys = system();
    let mut work = StepWorkspace::new(&sys).unwrap();
    let port = |_: &[f64], e: &[f64], out: &mut [f64]| {
        out[0] = -1e16;
        Some(NonlinearPortPower {
            supplied_power_w: -e[0] * out[0] + 1.0,
            dissipated_power_w: 1.0,
        })
    };
    let (mut next, mut y) = ([0.0], [0.0]);
    let receipt = work
        .step_into_supplied_controlled(
            &sys,
            &[1.0],
            &[0.0],
            1e-18,
            &mut next,
            &mut y,
            &port,
            None,
            None,
            || false,
        )
        .unwrap();
    close(next[0], 1.01);
    assert_eq!(receipt.dissipated, 1e-18);
    assert!((receipt.supplied - 1e-18 * y[0] * 1e16).abs()
        <= 4.0*f64::EPSILON*receipt.supplied.abs());
    let bound = y[0].abs()*receipt.solver_residual + 32.0*f64::EPSILON
        *(receipt.delta_h.abs()+receipt.supplied.abs()+receipt.dissipated);
    assert!(receipt.balance_residual().abs() <= bound, "{receipt:?}, bound={bound}");
}
