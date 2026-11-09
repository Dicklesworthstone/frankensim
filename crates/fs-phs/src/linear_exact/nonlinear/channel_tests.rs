use super::*;

fn split(n: usize) -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: n,
        max_omega_dt: 10.0,
    }
}
fn dense() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
        max_dense_bytes: usize::MAX,
        max_multiply_terms: u64::MAX,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn witness(n: usize, row: &[f64]) -> ResistanceWitness {
    let matrix = (0..n)
        .flat_map(|i| (0..n).map(move |j| row[i] * row[j]))
        .collect();
    ResistanceWitness::try_from_source_rows(n, matrix, row.to_vec()).unwrap()
}
fn rk4<const N: usize>(
    initial: [f64; N],
    dt: f64,
    count: usize,
    rhs: impl Fn(&[f64; N]) -> [f64; N],
) -> [f64; N] {
    let mut state = initial;
    let h = dt / count as f64;
    for _ in 0..count {
        let a = rhs(&state);
        let b = rhs(&std::array::from_fn(|i| state[i] + 0.5 * h * a[i]));
        let c = rhs(&std::array::from_fn(|i| state[i] + 0.5 * h * b[i]));
        let d = rhs(&std::array::from_fn(|i| state[i] + h * c[i]));
        for i in 0..N {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}

#[test]
fn full_cross_resistance_channels_and_tiny_relaxation_use_positive_integrals() {
    let rows = [[0.5, 0.25], [-0.2, 0.4]];
    let r: Vec<f64> = (0..2)
        .flat_map(|i| {
            let rows = rows;
            (0..2).map(move |j| rows.iter().map(|v| v[i] * v[j]).sum())
        })
        .collect();
    let channels: Vec<_> = rows.iter().map(|row| witness(2, row)).collect();
    let make = || {
        PortHamiltonian::new(
            2,
            1,
            vec![0.0, 1.0, -1.0, 0.0],
            r.clone(),
            vec![1.0, 0.2],
            Box::new(QuadraticStorage::new(vec![2.0, 0.3, 0.3, 1.2], 2).unwrap()),
        )
        .unwrap()
        .with_resistance_witness(
            ResistanceWitness::try_from_source_rows(2, r.clone(), rows.concat()).unwrap(),
        )
        .unwrap()
    };
    let dt = 0.37;
    let mut observed = PreparedNonlinearZoh::new_owned_with_resistance_channels(
        make(),
        dt,
        split(2),
        dense(),
        &channels,
    )
    .unwrap();
    let mut plain = PreparedNonlinearZoh::new_owned(make(), dt, split(2)).unwrap();
    let (mut next, mut output, mut heat) = ([0.0; 2], [0.0; 1], [0.0; 2]);
    let record = observed
        .step_observed_into(&[0.7, -0.3], &[0.17], &mut next, &mut output, &mut heat)
        .unwrap();
    let original = plain.step(&[0.7, -0.3], &[0.17]).unwrap();
    assert_eq!(next.as_slice(), original.x);
    assert_eq!(record.dissipated, original.dissipated);
    let reference = rk4([0.7, -0.3, 0.0, 0.0], dt, 2048, |x| {
        let e = [2.0 * x[0] + 0.3 * x[1], 0.3 * x[0] + 1.2 * x[1]];
        [
            e[1] - r[0] * e[0] - r[1] * e[1] + 0.17,
            -e[0] - r[2] * e[0] - r[3] * e[1] + 0.2 * 0.17,
            (rows[0][0] * e[0] + rows[0][1] * e[1]).powi(2),
            (rows[1][0] * e[0] + rows[1][1] * e[1]).powi(2),
        ]
    });
    for i in 0..2 {
        near(next[i], reference[i], 2e-12);
        near(heat[i], reference[2 + i], 2e-13);
    }
    near(heat.iter().sum(), record.dissipated, 2e-13);
    assert!(heat.iter().all(|q| *q > 0.0));

    let make_tiny = || {
        PortHamiltonian::new(
            1,
            0,
            vec![0.0],
            vec![5.0],
            vec![],
            Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
        )
        .unwrap()
    };
    let channel = witness(1, &[5.0_f64.sqrt()]);
    let mut tiny = PreparedNonlinearZoh::new_owned_with_resistance_channels(
        make_tiny(),
        4.0,
        split(1),
        dense(),
        &[channel],
    )
    .unwrap();
    let (mut x, mut y, mut q) = ([0.0], [], [0.0]);
    tiny.step_observed_into(&[1e-12], &[], &mut x, &mut y, &mut q)
        .unwrap();
    // Channel observation must preserve the original exponential endpoint
    // bit for bit. The analytical comparison allows 1e-12 relative error
    // after the original matrix exponential's scaling and squaring.
    let mut plain_tiny = PreparedNonlinearZoh::new_owned(make_tiny(), 4.0, split(1)).unwrap();
    let plain_record = plain_tiny.step(&[1e-12], &[]).unwrap();
    assert_eq!(x.as_slice(), plain_record.x);
    let analytical_endpoint = 1e-12 * (-20.0_f64).exp();
    near(x[0], analytical_endpoint, 1e-12 * analytical_endpoint.abs());
    near(q[0], 0.5e-24 * (1.0 - (-40.0_f64).exp()), 2e-36);
    assert!(q[0] > 0.0);
}

const Q: [f64; 9] = [6.25, 0.0, -1.5, 0.0, 1.0, 0.0, -1.5, 0.0, 1.0];
struct QuarticMemory;
impl Storage for QuarticMemory {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        2.0 * x[0] * x[0]
            + 0.5 * x[1] * x[1]
            + 0.5 * (1.5 * x[0] - x[2]).powi(2)
            + 0.75 * x[0].powi(4)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = 6.25 * x[0] - 1.5 * x[2] + 3.0 * x[0].powi(3);
        out[1] = x[1];
        out[2] = x[2] - 1.5 * x[0];
    }
    fn quadratic_part(&self) -> Option<super::super::QuadraticPart> {
        Some(super::super::QuadraticPart {
            matrix: Q.to_vec(),
            exact: false,
        })
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        let delta = b[0] - a[0];
        Some(quadratic_midpoint_remainder(&Q, a, b, mid) + 0.75 * mid[0] * delta.powi(3))
    }
}
fn nonlinear() -> PortHamiltonian {
    PortHamiltonian::new(
        3,
        1,
        vec![0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        vec![0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 5.0],
        vec![0.0, 1.0, 0.0],
        Box::new(QuarticMemory),
    )
    .unwrap()
}
fn nonlinear_channels() -> [ResistanceWitness; 2] {
    [
        witness(3, &[0.0, 0.4_f64.sqrt(), 0.0]),
        witness(3, &[0.0, 0.0, 5.0_f64.sqrt()]),
    ]
}
fn nonlinear_rhs(x: &[f64; 5], held: Option<[f64; 3]>) -> [f64; 5] {
    let f = held.unwrap_or([3.0 * x[0].powi(3), 0.0, 0.0]);
    let e = [
        6.25 * x[0] - 1.5 * x[2] + f[0],
        x[1] + f[1],
        x[2] - 1.5 * x[0] + f[2],
    ];
    [
        e[1],
        -e[0] - 0.4 * e[1] + 0.2,
        -5.0 * e[2],
        0.4 * e[1] * e[1],
        5.0 * e[2] * e[2],
    ]
}

#[test]
fn nonlinear_channels_follow_the_same_implicit_trajectory_and_preserve_motion() {
    let initial = [0.7, -0.3, 0.1];
    let dt = 0.18;
    let mut observed = PreparedNonlinearZoh::new_owned_with_resistance_channels(
        nonlinear(),
        dt,
        split(3),
        dense(),
        &nonlinear_channels(),
    )
    .unwrap();
    let mut plain = PreparedNonlinearZoh::new_owned(nonlinear(), dt, split(3)).unwrap();
    let (mut next, mut output, mut heat) = ([0.0; 3], [0.0; 1], [0.0; 2]);
    let record = observed
        .step_observed_into(&initial, &[0.2], &mut next, &mut output, &mut heat)
        .unwrap();
    let original = plain.step(&initial, &[0.2]).unwrap();
    assert_eq!(next.as_slice(), original.x);
    assert_eq!(record.dissipated, original.dissipated);
    let mid: [f64; 3] = std::array::from_fn(|i| 0.5 * (initial[i] + next[i]));
    let delta: [f64; 3] = std::array::from_fn(|i| next[i] - initial[i]);
    let correction = 0.75 * mid[0] * delta[0].powi(3) / delta.iter().map(|d| d * d).sum::<f64>();
    let held = [
        3.0 * mid[0].powi(3) + correction * delta[0],
        correction * delta[1],
        correction * delta[2],
    ];
    let reference = rk4(
        [initial[0], initial[1], initial[2], 0.0, 0.0],
        dt,
        2048,
        |x| nonlinear_rhs(x, Some(held)),
    );
    for i in 0..3 {
        near(next[i], reference[i], 4e-12);
    }
    for i in 0..2 {
        near(heat[i], reference[i + 3], 3e-13);
    }
    near(heat.iter().sum(), record.dissipated, 3e-13);
    let midpoint_heat = dt * 5.0 * (mid[2] - 1.5 * mid[0]).powi(2);
    assert!(
        (heat[1] - midpoint_heat).abs() > 1e-3,
        "endpoint/midpoint power must be distinguishable"
    );
}

#[test]
fn nonlinear_motion_and_channel_heat_converge_to_continuous_material_dynamics() {
    let initial = [0.7, -0.3, 0.1, 0.0, 0.0];
    let reference = rk4(initial, 0.4, 32768, |x| nonlinear_rhs(x, None));
    let mut errors = Vec::new();
    for count in [20, 40, 80] {
        let mut owner = PreparedNonlinearZoh::new_owned_with_resistance_channels(
            nonlinear(),
            0.4 / count as f64,
            split(3),
            dense(),
            &nonlinear_channels(),
        )
        .unwrap();
        let mut state = [initial[0], initial[1], initial[2]];
        let (mut next, mut output, mut heat) = ([0.0; 3], [0.0; 1], [0.0; 2]);
        let mut integrated = [0.0; 2];
        for _ in 0..count {
            owner
                .step_observed_into(&state, &[0.2], &mut next, &mut output, &mut heat)
                .unwrap();
            state = next;
            for i in 0..2 {
                integrated[i] += heat[i];
            }
        }
        let error = (0..3)
            .map(|i| (state[i] - reference[i]).abs())
            .chain((0..2).map(|i| (integrated[i] - reference[i + 3]).abs()))
            .fold(0.0_f64, f64::max);
        errors.push(error);
    }
    assert!(
        errors[1] < errors[0] / 3.7 && errors[2] < errors[1] / 3.7,
        "{errors:?}"
    );
    assert!(errors[2] < 2e-5, "{errors:?}");
}

#[test]
fn channel_admission_and_failed_steps_leave_all_outputs_untouched() {
    let channels = nonlinear_channels();
    let mut owner = PreparedNonlinearZoh::new_owned_with_resistance_channels(
        nonlinear(),
        0.18,
        split(3),
        dense(),
        &channels,
    )
    .unwrap();
    assert_eq!(owner.resistance_channel_count(), 2);
    owner.set_max_iterations(0).unwrap();
    assert_eq!(owner.iteration_limit(), 0);
    let (mut next, mut output, mut heat) = ([9.0; 3], [8.0; 1], [7.0; 2]);
    assert!(
        owner
            .step_observed_into(&[0.7, -0.3, 0.1], &[0.2], &mut next, &mut output, &mut heat)
            .is_err()
    );
    assert_eq!((next, output, heat), ([9.0; 3], [8.0], [7.0; 2]));
    owner.set_max_iterations(MAX_NEWTON).unwrap();
    assert!(matches!(
        owner.step_observed_into_controlled(
            &[0.7, -0.3, 0.1],
            &[0.2],
            &mut next,
            &mut output,
            &mut heat,
            || true,
        ),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!((next, output, heat), ([9.0; 3], [8.0], [7.0; 2]));
    assert!(
        owner
            .step_observed_into(
                &[0.7, -0.3, 0.1],
                &[0.2],
                &mut next,
                &mut output,
                &mut heat[..1],
            )
            .is_err()
    );
    assert_eq!((next, output, heat), ([9.0; 3], [8.0], [7.0; 2]));
    let mut byte_budget = dense();
    byte_budget.max_dense_bytes = 1;
    assert!(
        PreparedNonlinearZoh::new_owned_with_resistance_channels(
            nonlinear(),
            0.18,
            split(3),
            byte_budget,
            &channels,
        )
        .is_err()
    );
    let mut work_budget = dense();
    work_budget.max_multiply_terms = 1;
    assert!(
        PreparedNonlinearZoh::new_owned_with_resistance_channels(
            nonlinear(),
            0.18,
            split(3),
            work_budget,
            &channels,
        )
        .is_err()
    );
    assert!(
        PreparedNonlinearZoh::new_owned_with_resistance_channels(
            nonlinear(),
            0.18,
            split(3),
            dense(),
            &[witness(1, &[1.0])],
        )
        .is_err()
    );
    let mut fresh = PreparedNonlinearZoh::new_owned_with_resistance_channels(
        nonlinear(),
        0.18,
        split(3),
        dense(),
        &channels,
    )
    .unwrap();
    owner
        .step_observed_into(&[0.7, -0.3, 0.1], &[0.2], &mut next, &mut output, &mut heat)
        .unwrap();
    let (mut expected, mut expected_output, mut expected_heat) = ([0.0; 3], [0.0; 1], [0.0; 2]);
    fresh
        .step_observed_into(
            &[0.7, -0.3, 0.1],
            &[0.2],
            &mut expected,
            &mut expected_output,
            &mut expected_heat,
        )
        .unwrap();
    assert_eq!(
        (next, output, heat),
        (expected, expected_output, expected_heat)
    );
}
