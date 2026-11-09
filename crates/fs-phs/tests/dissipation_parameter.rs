//! G0/G1: original nonlinear loss parameters, full effort chain and weak work.
use fs_phs::harmonic::equilibrium::{LinearizationError, NonlinearDissipation};
use fs_phs::{PhsStructureParameterDirection, PortHamiltonian, QuadraticStorage};
use std::cell::Cell;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:e} != {b:e}"
    );
}
fn original(x: &[f64], e: &[f64], parameter: f64) -> [f64; 3] {
    let beta = 0.4 - 0.12 * parameter;
    let sigma = 0.3 + 0.2 * parameter;
    let eta = 0.2 + 0.3 * parameter;
    let w = e[1] + eta * e[2];
    let f = beta * (1.0 + sigma * x[0] * x[0]) * w.powi(3);
    [0.0, f, eta * f]
}
fn flow(x: &[f64], e: &[f64], out: &mut [f64]) -> bool {
    out.copy_from_slice(&original(x, e, 0.0));
    true
}
fn tangent(x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]) -> bool {
    let w = e[1] + 0.2 * e[2];
    let d = 0.4
        * (0.6 * x[0] * dx[0] * w.powi(3)
            + 3.0 * (1.0 + 0.3 * x[0] * x[0]) * w * w * (de[1] + 0.2 * de[2]));
    out.copy_from_slice(&[0.0, d, 0.2 * d]);
    true
}
fn direct(x: &[f64], e: &[f64], out: &mut [f64]) -> bool {
    let w = e[1] + 0.2 * e[2];
    let factor = 1.0 + 0.3 * x[0] * x[0];
    let d = (-0.12 * factor + 0.08 * x[0] * x[0]) * w.powi(3) + 1.2 * factor * w * w * 0.3 * e[2];
    out.copy_from_slice(&[0.0, d, 0.2 * d + 0.3 * 0.4 * factor * w.powi(3)]);
    true
}
fn law() -> NonlinearDissipation<'static> {
    NonlinearDissipation {
        flow: &flow,
        tangent: &tangent,
        max_tangent_actions: 3,
    }
}

#[test]
fn complete_original_loss_parameter_and_physical_field_match_neighboring_sources() {
    let calls = Cell::new(0);
    let action = |x: &[f64], e: &[f64], dx: &[f64], de: &[f64], out: &mut [f64]| {
        assert!(dx.iter().all(|v| *v == 0.0));
        calls.set(calls.get() + 1);
        tangent(x, e, dx, de, out)
    };
    let source = NonlinearDissipation {
        flow: &flow,
        tangent: &action,
        max_tangent_actions: 1,
    };
    let mut loss = source
        .prepare_parameter_direction(3, Some(&direct))
        .unwrap();
    let x = [0.8, -0.2, 0.4];
    let e = [2.1, -0.7, 0.2];
    let ep = [-0.4, 0.5, -0.8];
    let sample = loss.evaluate(&x, &e, &ep, || false).unwrap();
    assert_eq!(calls.get(), 1);
    let j = [0.0, 1.0, -0.2, -1.0, 0.0, 0.1, 0.2, -0.1, 0.0];
    let r = [0.1, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.3];
    let g = [0.0, 1.0, 0.4];
    let jp = [0.0, 0.2, 0.0, -0.2, 0.0, 0.15, 0.0, -0.15, 0.0];
    let rp = [0.03, 0.0, 0.0, 0.0, -0.02, 0.0, 0.0, 0.0, 0.04];
    let gp = [0.1, -0.05, 0.3];
    let system = PortHamiltonian::new(
        3,
        1,
        j.to_vec(),
        r.to_vec(),
        g.to_vec(),
        Box::new(
            QuadraticStorage::new(vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], 3).unwrap(),
        ),
    )
    .unwrap();
    let direction = system
        .prepare_parameter_direction(PhsStructureParameterDirection {
            interconnection: Some(&jp),
            dissipation: Some(&rp),
            ports: Some(&gp),
        })
        .unwrap();
    let mut result = [0.0; 3];
    direction
        .vector_field_with_dissipation_from_effort_into(
            &e,
            &ep,
            &[0.4],
            &[-0.2],
            sample.effort_chain,
            sample.direct_flow,
            &mut result,
        )
        .unwrap();
    let nearby = |p: f64| {
        let e = core::array::from_fn::<_, 3, _>(|i| e[i] + p * ep[i]);
        let d = original(&x, &e, p);
        let field = core::array::from_fn::<_, 3, _>(|i| {
            (0..3)
                .map(|k| {
                    (j[3 * i + k] + p * jp[3 * i + k] - r[3 * i + k] - p * rp[3 * i + k]) * e[k]
                })
                .sum::<f64>()
                + (g[i] + p * gp[i]) * (0.4 - 0.2 * p)
                - d[i]
        });
        let power = (0..3).map(|i| e[i] * d[i]).sum::<f64>();
        (field, power, d)
    };
    for h in [1e-4, 5e-5] {
        let (plus, qp, dp) = nearby(h);
        let (minus, qm, dm) = nearby(-h);
        for i in 0..3 {
            close(result[i], (plus[i] - minus[i]) / (2.0 * h), 2e-8);
            close(
                sample.effort_chain[i] + sample.direct_flow[i],
                (dp[i] - dm[i]) / (2.0 * h),
                2e-8,
            );
        }
        close(sample.power_derivative, (qp - qm) / (2.0 * h), 2e-8);
    }
    close(sample.power, nearby(0.0).1, 1e-15);
}

#[test]
fn fixed_constitutive_law_keeps_effort_chain_and_signed_parameter_heat() {
    let flow = |_: &[f64], e: &[f64], out: &mut [f64]| {
        out[0] = 0.7 * e[0].powi(3);
        true
    };
    let action = |_: &[f64], e: &[f64], _: &[f64], de: &[f64], out: &mut [f64]| {
        out[0] = 2.1 * e[0] * e[0] * de[0];
        true
    };
    let direct = |_: &[f64], e: &[f64], out: &mut [f64]| {
        out[0] = -2.0 * e[0].powi(3);
        true
    };
    let law = NonlinearDissipation {
        flow: &flow,
        tangent: &action,
        max_tangent_actions: 1,
    };
    let mut fixed = law.prepare_parameter_direction(1, None).unwrap();
    let result = fixed.evaluate(&[0.0], &[2.0], &[0.3], || false).unwrap();
    assert_eq!(result.direct_flow, &[0.0]);
    close(result.effort_chain[0], 2.1 * 4.0 * 0.3, 1e-15);
    close(result.power_derivative, 4.0 * 0.7 * 8.0 * 0.3, 1e-15);
    let mut changed = law.prepare_parameter_direction(1, Some(&direct)).unwrap();
    let result = changed.evaluate(&[0.0], &[2.0], &[0.3], || false).unwrap();
    assert!(result.power > 0.0 && result.power_derivative < 0.0);
    close(result.power_derivative, 4.0 * 0.7 * 8.0 * 0.3 - 32.0, 1e-15);
}

#[test]
fn weak_parameter_work_and_flow_survive_opposing_large_source_contributions() {
    const BIG: f64 = 18_014_398_509_481_984.0;
    let flow = |_: &[f64], e: &[f64], out: &mut [f64]| {
        out[0] = BIG * (2.0 - e[0]);
        true
    };
    let action = |_: &[f64], _: &[f64], _: &[f64], de: &[f64], out: &mut [f64]| {
        out[0] = -BIG * de[0];
        true
    };
    let direct = |_: &[f64], _: &[f64], out: &mut [f64]| {
        out[0] = 1.0;
        true
    };
    let law = NonlinearDissipation {
        flow: &flow,
        tangent: &action,
        max_tangent_actions: 1,
    };
    let mut prepared = law.prepare_parameter_direction(1, Some(&direct)).unwrap();
    let sample = prepared.evaluate(&[0.0], &[1.0], &[1.0], || false).unwrap();
    assert_eq!(sample.power, BIG);
    assert_eq!(sample.power_derivative, 1.0);
    let system = PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![0.0],
        vec![1.0],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
    )
    .unwrap();
    let direction = system
        .prepare_parameter_direction(PhsStructureParameterDirection::default())
        .unwrap();
    let mut result = [0.0];
    direction
        .vector_field_with_dissipation_from_effort_into(
            &[1.0],
            &[1.0],
            &[0.0],
            &[-BIG],
            sample.effort_chain,
            sample.direct_flow,
            &mut result,
        )
        .unwrap();
    assert_eq!(result, [-1.0]);
}

#[test]
fn incomplete_original_partials_and_cancellation_refuse_without_poisoning_reuse() {
    let complete = Cell::new(false);
    let partial = |_: &[f64], _: &[f64], out: &mut [f64]| {
        if complete.get() {
            out.fill(0.1);
        } else {
            out[0] = 0.1;
        }
        true
    };
    let mut prepared = law()
        .prepare_parameter_direction(3, Some(&partial))
        .unwrap();
    let (x, e, ep) = ([0.2, 0.1, 0.3], [0.3, -0.4, 0.2], [0.1, 0.3, -0.2]);
    assert!(matches!(
        prepared.evaluate(&x, &e, &ep, || false),
        Err(LinearizationError::DissipationParameterRefused { .. })
    ));
    complete.set(true);
    let mut polls = 0;
    assert!(matches!(
        prepared.evaluate(&x, &e, &ep, || {
            polls += 1;
            polls == 3
        }),
        Err(LinearizationError::Cancelled)
    ));
    let result = prepared.evaluate(&x, &e, &ep, || false).unwrap();
    assert_eq!(result.direct_flow, &[0.1; 3]);
    assert!(prepared.evaluate(&x, &e, &[0.1], || false).is_err());
    assert!(law().prepare_parameter_direction(usize::MAX, None).is_err());
    assert!(
        NonlinearDissipation {
            max_tangent_actions: 0,
            ..law()
        }
        .prepare_parameter_direction(3, None)
        .is_err()
    );
    let active = |_: &[f64], e: &[f64], out: &mut [f64]| {
        for (i, o) in out.iter_mut().enumerate() {
            *o = -e[i];
        }
        true
    };
    let mut prepared = NonlinearDissipation {
        flow: &active,
        ..law()
    }
    .prepare_parameter_direction(3, None)
    .unwrap();
    assert!(matches!(
        prepared.evaluate(&x, &e, &ep, || false),
        Err(LinearizationError::ActiveDissipation { .. })
    ));
}
