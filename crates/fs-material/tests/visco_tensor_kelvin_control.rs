//! G1/G3/G4: original Kelvin/Maxwell stress loading, heat and history retention.
use fs_material::visco_tensor::kelvin_control::{
    PreparedTensorKelvinStressControl as Controller, TensorKelvinStressControlBudget as Budget,
};
use fs_material::visco_tensor::{TensorProny, TensorPronyBranch, TensorPronyError};
use fs_material::{Tangent6, Voigt};

#[path = "visco_tensor_kelvin_control/semidefinite.rs"]
mod semidefinite;
#[path = "visco_tensor_kelvin_control/temperature.rs"]
mod temperature;

const W: Voigt = [
    1.0,
    1.0,
    1.0,
    std::f64::consts::SQRT_2,
    std::f64::consts::SQRT_2,
    std::f64::consts::SQRT_2,
];
fn diagonal(values: Voigt) -> Tangent6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
}
fn dot(a: Voigt, b: Voigt) -> f64 {
    (0..6).map(|i| a[i] * b[i]).sum()
}
fn apply(a: &Tangent6, x: Voigt) -> Voigt {
    core::array::from_fn(|i| dot(a[i], x))
}
fn tensor(x: Voigt) -> Voigt {
    core::array::from_fn(|i| x[i] / W[i])
}
fn mandel(x: Voigt) -> Voigt {
    core::array::from_fn(|i| x[i] * W[i])
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn relative(actual: f64, expected: f64, tolerance: f64) {
    close(actual / expected, 1.0, tolerance);
}
fn kelvin(e: f64, viscosity: f64) -> TensorProny {
    TensorProny::new_mandel(diagonal([e; 6]), vec![], 0.1)
        .unwrap()
        .with_kelvin_viscosity_mandel(diagonal([viscosity; 6]))
        .unwrap()
}

#[test]
fn g1_kelvin_stress_jump_creep_ramp_and_recovery_match_closed_forms() {
    let (e, viscosity, stress, duration) = (400.0_f64, 75.0, 2.0, 0.7);
    let law = kelvin(e, viscosity);
    let mut history = law.history();
    let control = Controller::new(&history, Budget::default(), || false).unwrap();
    assert_eq!(control.state_coordinates(), 6);
    let target = [stress, 0.0, 0.0, 0.0, 0.0, 0.0];
    let jump = control
        .advance(&mut history, target, 0.0, || false)
        .unwrap();
    assert_eq!(jump.strain, [0.0; 6]);
    assert_eq!(jump.work, 0.0);
    assert_eq!(jump.dissipated, 0.0);
    relative(jump.strain_rate[0], stress / viscosity, 2e-13);
    let hold = control
        .advance(&mut history, target, duration, || false)
        .unwrap();
    let r = duration * e / viscosity;
    let decay = (-r).exp();
    let expected = stress / e * -(-r).exp_m1();
    relative(hold.strain[0], expected, 3e-13);
    relative(hold.strain_rate[0], stress / viscosity * decay, 3e-12);
    relative(hold.work, stress * expected, 3e-13);
    relative(
        hold.dissipated,
        stress * stress / (2.0 * e) * -(-2.0 * r).exp_m1(),
        3e-13,
    );
    close(hold.balance_residual, 0.0, 2e-17);
    let release = control
        .advance(&mut history, [0.0; 6], 0.0, || false)
        .unwrap();
    assert_eq!(release.strain, hold.strain);
    assert_eq!(release.work, 0.0);
    let recovery = control
        .advance(&mut history, [0.0; 6], duration, || false)
        .unwrap();
    relative(recovery.strain[0], expected * decay, 3e-12);
    assert_eq!(recovery.work, 0.0);
    relative(
        recovery.dissipated,
        0.5 * e * expected * expected * -(-2.0 * r).exp_m1(),
        3e-13,
    );

    let mut ramp = law.history();
    let ramp_control = Controller::new(&ramp, Budget::default(), || false).unwrap();
    let step = ramp_control
        .advance(&mut ramp, target, duration, || false)
        .unwrap();
    let b = -(-r).exp_m1() / r;
    let first = (1.0 - (1.0 + r) * decay) / (r * r);
    relative(step.strain[0], stress / e * (1.0 - b), 3e-13);
    relative(step.work, stress * stress / e * (0.5 - first), 3e-13);
    let heat = stress * stress / (e * r) * (1.0 - 2.0 * b + -(-2.0 * r).exp_m1() / (2.0 * r));
    relative(step.dissipated, heat, 3e-13);
    close(step.balance_residual, 0.0, 2e-17);
}

#[test]
fn g1_maxwell_and_kelvin_parallel_creep_matches_two_physical_rates_and_heat() {
    let (e, branch, viscosity, tau, stress, duration) =
        (300.0_f64, 1200.0_f64, 140.0, 0.65, 1.8, 2.4);
    let mut factor = [[0.0; 6]; 6];
    factor[0][0] = branch.sqrt();
    let law = TensorProny::new_mandel(
        diagonal([e; 6]),
        vec![TensorPronyBranch::from_gram_rows(factor, tau).unwrap()],
        0.1,
    )
    .unwrap()
    .with_kelvin_viscosity_mandel(diagonal([viscosity; 6]))
    .unwrap();
    let mut history = law.history();
    let control = Controller::new(&history, Budget::default(), || false).unwrap();
    let target = [stress, 0.0, 0.0, 0.0, 0.0, 0.0];
    control
        .advance(&mut history, target, 0.0, || false)
        .unwrap();
    let step = control
        .advance(&mut history, target, duration, || false)
        .unwrap();
    let a = (viscosity + tau * (e + branch)) / (viscosity * tau);
    let c = e / (viscosity * tau);
    let fast = 0.5 * (a + (a * a - 4.0 * c).sqrt());
    let slow = c / fast;
    let x = (stress / viscosity - fast * stress / e) / (fast - slow);
    let y = -stress / e - x;
    let strain = stress / e + x * (-slow * duration).exp() + y * (-fast * duration).exp();
    let velocities = [-slow * x, -fast * y];
    let stresses = [
        -e * x - viscosity * velocities[0],
        -e * y - viscosity * velocities[1],
    ];
    let integral = |rate: f64| -(-rate * duration).exp_m1() / rate;
    let mut heat = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            heat += (viscosity * velocities[i] * velocities[j]
                + stresses[i] * stresses[j] / (branch * tau))
                * integral([slow, fast][i] + [slow, fast][j]);
        }
    }
    relative(step.strain[0], strain, 3e-12);
    relative(
        step.strain_rate[0],
        velocities[0] * (-slow * duration).exp() + velocities[1] * (-fast * duration).exp(),
        3e-12,
    );
    relative(step.work, stress * strain, 3e-12);
    relative(step.dissipated, heat, 3e-12);
    close(step.balance_residual, 0.0, 2e-16);
}

fn gram(b: &Tangent6) -> Tangent6 {
    core::array::from_fn(|i| core::array::from_fn(|j| (0..6).map(|k| b[k][i] * b[k][j]).sum()))
}
fn coupled() -> (TensorProny, Tangent6, Tangent6, [Tangent6; 2]) {
    let mut eq = diagonal([1200.0, 1700.0, 1000.0, 900.0, 700.0, 1300.0]);
    eq[0][3] = -130.0;
    eq[3][0] = -130.0;
    eq[1][5] = 70.0;
    eq[5][1] = 70.0;
    let mut viscous = diagonal([12.0, 15.0, 10.0, 9.0, 8.0, 13.0]);
    viscous[0][1] = 3.0;
    viscous[1][4] = -2.0;
    viscous[2][5] = 4.0;
    viscous[3][0] = -2.0;
    viscous[4][2] = 1.5;
    let viscosity = gram(&viscous);
    let mut factors = [[[0.0; 6]; 6]; 2];
    factors[0][0] = [17.0, -3.0, 2.0, 7.0, 1.0, 0.0];
    factors[1][0] = [1.0, 12.0, -2.0, 6.0, 3.0, -4.0];
    factors[1][1] = [-3.0, 2.0, 8.0, -4.0, 0.0, 4.0];
    let stiffness = factors.map(|b| gram(&b));
    let branches = factors
        .into_iter()
        .zip([0.23, 1.1])
        .map(|(b, tau)| TensorPronyBranch::from_gram_rows(b, tau).unwrap())
        .collect();
    (
        TensorProny::new_mandel(eq, branches, 0.1)
            .unwrap()
            .with_kelvin_viscosity_mandel(viscosity)
            .unwrap(),
        eq,
        viscosity,
        stiffness,
    )
}

// Independent elimination on the original viscosity tensor; no prepared modal
// factor or Prony ramp is used in this physical-coordinate ODE oracle.
fn solve(mut a: Tangent6, mut rhs: Voigt) -> Voigt {
    for j in 0..6 {
        let pivot = (j..6)
            .max_by(|&i, &k| a[i][j].abs().total_cmp(&a[k][j].abs()))
            .unwrap();
        a.swap(j, pivot);
        rhs.swap(j, pivot);
        let row = a[j];
        for i in j + 1..6 {
            let multiplier = a[i][j] / row[j];
            for (value, pivot) in a[i][j..].iter_mut().zip(&row[j..]) {
                *value -= multiplier * pivot;
            }
            rhs[i] -= multiplier * rhs[j];
        }
    }
    let mut x = [0.0; 6];
    for i in (0..6).rev() {
        x[i] = (rhs[i] - (i + 1..6).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
    }
    x
}
// Original strain, two original elastic branch strains, accumulated work/heat.
fn derivative(
    y: &[f64; 20],
    stress: Voigt,
    eq: &Tangent6,
    viscosity: &Tangent6,
    branches: &[Tangent6; 2],
) -> [f64; 20] {
    let strain = core::array::from_fn(|i| y[i]);
    let equilibrium = apply(eq, strain);
    let mut rhs = core::array::from_fn(|i| stress[i] - equilibrium[i]);
    let mut heat = 0.0;
    for (j, tau) in [0.23, 1.1].into_iter().enumerate() {
        let q = core::array::from_fn(|i| y[6 + 6 * j + i]);
        let sigma = apply(&branches[j], q);
        for i in 0..6 {
            rhs[i] -= sigma[i];
        }
        heat += dot(q, sigma) / tau;
    }
    let rate = solve(*viscosity, rhs);
    let mut result = [0.0; 20];
    result[..6].copy_from_slice(&rate);
    for (j, tau) in [0.23, 1.1].into_iter().enumerate() {
        for i in 0..6 {
            result[6 + 6 * j + i] = rate[i] - y[6 + 6 * j + i] / tau;
        }
    }
    result[18] = dot(stress, rate);
    result[19] = heat + dot(rate, apply(viscosity, rate));
    result
}
#[allow(clippy::too_many_arguments)] // Independent original-coordinate ODE and its complete tensors.
fn integrate(
    y: &mut [f64; 20],
    start: Voigt,
    target: Voigt,
    dt: f64,
    count: usize,
    eq: &Tangent6,
    viscosity: &Tangent6,
    branches: &[Tangent6; 2],
) -> [Voigt; 2] {
    let h = dt / count as f64;
    let mut bounds = [core::array::from_fn(|i| y[i]); 2];
    for step in 0..count {
        let stress = |offset: f64| {
            core::array::from_fn(|i| {
                start[i] + (step as f64 + offset) / count as f64 * (target[i] - start[i])
            })
        };
        let k1 = derivative(y, stress(0.0), eq, viscosity, branches);
        let k2 = derivative(
            &core::array::from_fn(|i| y[i] + h * 0.5 * k1[i]),
            stress(0.5),
            eq,
            viscosity,
            branches,
        );
        let k3 = derivative(
            &core::array::from_fn(|i| y[i] + h * 0.5 * k2[i]),
            stress(0.5),
            eq,
            viscosity,
            branches,
        );
        let k4 = derivative(
            &core::array::from_fn(|i| y[i] + h * k3[i]),
            stress(1.0),
            eq,
            viscosity,
            branches,
        );
        for i in 0..20 {
            y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        for i in 0..6 {
            bounds[0][i] = bounds[0][i].min(y[i]);
            bounds[1][i] = bounds[1][i].max(y[i]);
        }
    }
    bounds
}

#[test]
fn g1_noncommuting_full_tensors_match_independent_physical_ode_refinement() {
    let (law, eq, viscosity, branches) = coupled();
    let mut history = law.history();
    let control = Controller::new(&history, Budget::default(), || false).unwrap();
    assert_eq!(control.state_coordinates(), 18);
    let targets = [
        ([1.5, -0.8, 0.6, 0.7, -0.5, 0.4], 0.42),
        ([-0.7, 1.1, -0.3, -0.4, 0.5, 0.2], 0.61),
        ([0.0; 6], 0.35),
    ];
    let mut reports = Vec::new();
    for (target, dt) in targets {
        reports.push(
            control
                .advance(&mut history, tensor(target), dt, || false)
                .unwrap(),
        );
    }
    let mut errors = [0.0_f64; 3];
    for (rung, count) in [32, 64, 128].into_iter().enumerate() {
        let mut y = [0.0; 20];
        let mut start = [0.0; 6];
        for ((target, dt), report) in targets.into_iter().zip(&reports) {
            let bounds = integrate(&mut y, start, target, dt, count, &eq, &viscosity, &branches);
            for i in 0..6 {
                errors[rung] = errors[rung].max((y[i] - mandel(report.strain)[i]).abs());
                if rung == 2 {
                    assert!(bounds[0][i] >= mandel(report.strain_bounds[0])[i] - 3e-10);
                    assert!(bounds[1][i] <= mandel(report.strain_bounds[1])[i] + 3e-10);
                }
            }
            start = target;
        }
        if rung == 2 {
            close(y[18], history.work(), 2e-10);
            close(y[19], history.dissipated(), 2e-10);
            let strain = core::array::from_fn(|i| y[i]);
            let branch_energy = (0..2)
                .map(|j| {
                    let q = core::array::from_fn(|i| y[6 + 6 * j + i]);
                    0.5 * dot(q, apply(&branches[j], q))
                })
                .sum::<f64>();
            close(branch_energy, history.branch_stored_energy(), 2e-12);
            close(
                0.5 * dot(strain, apply(&eq, strain)) + branch_energy,
                history.stored_energy(),
                2e-12,
            );
            let rate = derivative(&y, [0.0; 6], &eq, &viscosity, &branches);
            for (i, &rate) in rate.iter().take(6).enumerate() {
                close(rate, mandel(reports.last().unwrap().strain_rate)[i], 2e-9);
            }
        }
    }
    assert!(
        errors[1] < errors[0] / 8.0 && errors[2] < errors[1] / 8.0,
        "{errors:?}"
    );
    for report in reports {
        assert!(report.dissipated >= 0.0);
        close(report.balance_residual, 0.0, 3e-16);
    }
}

#[test]
fn g3_subdivision_rotation_and_strain_control_keep_complete_original_memory() {
    let (law, _, _, _) = coupled();
    let rotation = [[0.6, -0.8, 0.0], [0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];
    let rotate = |v| fs_material::tensor::rotate(&v, &rotation);
    let mut original = law.history();
    let mut rotated = law.rotated(rotation).unwrap().history();
    let preload = [0.001, -0.0004, 0.0007, 0.0003, -0.0002, 0.0005];
    original.advance(preload, 0.4, || false).unwrap();
    rotated.advance(rotate(preload), 0.4, || false).unwrap();
    let control = Controller::new(&original, Budget::default(), || false).unwrap();
    let rotated_control = Controller::new(&rotated, Budget::default(), || false).unwrap();
    let mut split = original.clone();
    let target = [-0.5, 0.7, -0.2, 0.3, 0.4, -0.6];
    let middle = core::array::from_fn(|i| 0.5 * original.stress()[i] + 0.5 * target[i]);
    control
        .advance(&mut original, target, 0.73, || false)
        .unwrap();
    control
        .advance(&mut split, middle, 0.365, || false)
        .unwrap();
    control
        .advance(&mut split, target, 0.365, || false)
        .unwrap();
    rotated_control
        .advance(&mut rotated, rotate(target), 0.73, || false)
        .unwrap();
    let expected = rotate(original.strain());
    for (i, &expected) in expected.iter().enumerate() {
        close(original.strain()[i], split.strain()[i], 3e-13);
        close(rotated.strain()[i], expected, 3e-13);
    }
    close(original.work(), split.work(), 3e-13);
    close(original.dissipated(), split.dissipated(), 3e-13);
    close(original.work(), rotated.work(), 3e-13);
    close(original.dissipated(), rotated.dissipated(), 3e-13);
    let strain = [0.0007, 0.0001, -0.0003, 0.0004, -0.0005, 0.0002];
    let a = original.advance(strain, 0.31, || false).unwrap();
    let b = split.advance(strain, 0.31, || false).unwrap();
    let c = rotated.advance(rotate(strain), 0.31, || false).unwrap();
    for i in 0..6 {
        close(a.stress[i], b.stress[i], 2e-10);
        close(c.stress[i], rotate(a.stress)[i], 2e-10);
    }
    let a = control
        .advance(&mut original, [0.0; 6], 0.21, || false)
        .unwrap();
    let b = control
        .advance(&mut split, [0.0; 6], 0.21, || false)
        .unwrap();
    for i in 0..6 {
        close(a.strain[i], b.strain[i], 3e-13);
    }
}

#[test]
fn g1_weak_and_saturated_damping_retain_physical_heat_and_unrounded_increments() {
    let law = kelvin(1.0, 1e100);
    let mut history = law.history();
    let control = Controller::new(&history, Budget::default(), || false).unwrap();
    let target = [1e100, 0.0, 0.0, 0.0, 0.0, 0.0];
    control
        .advance(&mut history, target, 0.0, || false)
        .unwrap();
    // lambda*dt underflows, but the physical strain and heat do not.
    let tiny = control
        .advance(&mut history, target, 1e-300, || false)
        .unwrap();
    relative(tiny.strain_increment[0], 1e-300, 4e-13);
    relative(tiny.dissipated, 1e-200, 4e-13);
    relative(tiny.work, 1e-200, 4e-13);

    let law = kelvin(400.0, 1e300);
    let mut history = law.history();
    history
        .advance([0.01, 0.0, 0.0, 0.0, 0.0, 0.0], 1e300, || false)
        .unwrap();
    let control = Controller::new(&history, Budget::default(), || false).unwrap();
    let target = [4.1, 0.0, 0.0, 0.0, 0.0, 0.0];
    control
        .advance(&mut history, target, 0.0, || false)
        .unwrap();
    let old = history.strain();
    let weak = control
        .advance(&mut history, target, 1.0, || false)
        .unwrap();
    assert_eq!(history.strain(), old);
    relative(weak.strain_increment[0], 1e-301, 5e-12);
    relative(weak.stored_energy_change, 4e-301, 5e-12);
    relative(weak.dissipated, 1e-302, 5e-12);

    let law = kelvin(1.0, 1e-300);
    let mut history = law.history();
    let control = Controller::new(&history, Budget::default(), || false).unwrap();
    let target = [1e-3, 0.0, 0.0, 0.0, 0.0, 0.0];
    control
        .advance(&mut history, target, 0.0, || false)
        .unwrap();
    let stiff = control
        .advance(&mut history, target, 1e300, || false)
        .unwrap();
    relative(stiff.strain[0], 1e-3, 4e-13);
    relative(stiff.work, 1e-6, 4e-13);
    relative(stiff.dissipated, 5e-7, 4e-13);
    let mut small = law.history();
    let control = Controller::new(&small, Budget::default(), || false).unwrap();
    let target = [1e-200, 0.0, 0.0, 0.0, 0.0, 0.0];
    control.advance(&mut small, target, 0.0, || false).unwrap();
    // The loss-coordinate integral underflows, while the physical strain does
    // not. Amplitude-first integration must retain that source response.
    let step = control
        .advance(&mut small, target, 1e-300, || false)
        .unwrap();
    relative(step.strain[0], 1e-200 * -(-1.0_f64).exp_m1(), 4e-13);
}

#[test]
fn g4_zero_viscosity_budgets_binding_domain_and_cancellation_refuse_atomically() {
    let elastic = TensorProny::new_mandel(diagonal([100.0; 6]), vec![], 0.1).unwrap();
    assert!(Controller::new(&elastic.history(), Budget::default(), || false).is_err());
    let zero = elastic
        .with_kelvin_viscosity_mandel(diagonal([0.0; 6]))
        .unwrap();
    assert!(Controller::new(&zero.history(), Budget::default(), || false).is_err());
    let rotation = [[0.6, -0.8, 0.0], [0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];
    assert!(
        Controller::new(
            &zero.rotated(rotation).unwrap().history(),
            Budget::default(),
            || false
        )
        .is_err()
    );
    let (law, _, _, _) = coupled();
    let history = law.history();
    for budget in [
        Budget {
            max_state_coordinates: 5,
            ..Budget::default()
        },
        Budget {
            max_dense_scalars: 1,
            ..Budget::default()
        },
        Budget {
            max_preparation_work: 1,
            ..Budget::default()
        },
    ] {
        assert!(Controller::new(&history, budget, || false).is_err());
    }
    assert_eq!(
        Controller::new(&history, Budget::default(), || true).unwrap_err(),
        TensorPronyError::Cancelled
    );
    let control = Controller::new(&history, Budget::default(), || false).unwrap();
    let mut independent = law.history();
    assert!(
        control
            .advance(&mut independent, [0.0; 6], 0.1, || false)
            .is_err()
    );
    assert_eq!(independent, history);
    for (stress, dt) in [
        ([f64::NAN; 6], 1.0),
        ([1.0; 6], -1.0),
        ([1.0; 6], f64::INFINITY),
        ([1e6; 6], 100.0),
    ] {
        let mut trial = history.clone();
        assert!(control.advance(&mut trial, stress, dt, || false).is_err());
        assert_eq!(trial, history);
    }
    for stop in [1, 3, 12, 20] {
        let mut trial = history.clone();
        let mut polls = 0;
        let error = control
            .advance(&mut trial, [0.1; 6], 0.2, || {
                polls += 1;
                polls == stop
            })
            .unwrap_err();
        assert_eq!(error, TensorPronyError::Cancelled);
        assert_eq!(trial, history);
    }
}
