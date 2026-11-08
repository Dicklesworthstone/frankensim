//! G1/G3/G4: physical tensor creep against closed forms and an independent ODE.
use fs_material::visco_tensor::stress_control::{
    PreparedTensorStressControl as Controller, TensorStressControlBudget as Budget,
};
use fs_material::visco_tensor::{TensorProny, TensorPronyBranch, TensorPronyError};
use fs_material::{Tangent6, Voigt};

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
fn apply(a: &Tangent6, x: Voigt) -> Voigt {
    core::array::from_fn(|i| (0..6).map(|j| a[i][j] * x[j]).sum())
}
fn dot(a: Voigt, b: Voigt) -> f64 {
    (0..6).map(|i| a[i] * b[i]).sum()
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

fn sls(e: f64, branch: f64, tau: f64, limit: f64) -> TensorProny {
    let mut factor = [[0.0; 6]; 6];
    factor[0][0] = branch.sqrt();
    TensorProny::new_mandel(
        diagonal([e; 6]),
        vec![TensorPronyBranch::from_gram_rows(factor, tau).unwrap()],
        limit,
    )
    .unwrap()
}

#[test]
fn g1_standard_linear_solid_creep_and_zero_load_recovery_close_work_and_heat() {
    let (e, b, tau, stress, hold) = (400.0, 1600.0, 0.6, 2.0, 2.3);
    let law = sls(e, b, tau, 0.1);
    let mut history = law.history();
    let control = Controller::new(&history, Budget::default(), || false).unwrap();
    assert_eq!(control.memory_coordinates(), 6);
    let target = [stress, 0.0, 0.0, 0.0, 0.0, 0.0];
    let jump = control
        .advance(&mut history, target, 0.0, || false)
        .unwrap();
    let initial = stress / (e + b);
    relative(jump.strain[0], initial, 3e-13);
    relative(jump.work, 0.5 * stress * initial, 3e-13);
    assert_eq!(jump.dissipated.to_bits(), 0.0_f64.to_bits());
    let rate = e / ((e + b) * tau);
    let fraction = -(-rate * hold).exp_m1();
    let increase = (stress / e - initial) * fraction;
    let creep = control
        .advance(&mut history, target, hold, || false)
        .unwrap();
    relative(creep.strain[0], initial + increase, 3e-12);
    relative(creep.work, stress * increase, 3e-12);
    let heat = b * stress * stress / (2.0 * e * (e + b)) * -(-2.0 * rate * hold).exp_m1();
    relative(creep.dissipated, heat, 3e-12);
    close(creep.balance_residual, 0.0, 2e-14);
    let release = control
        .advance(&mut history, [0.0; 6], 0.0, || false)
        .unwrap();
    relative(release.strain[0], increase, 3e-12);
    let recovery = control
        .advance(&mut history, [0.0; 6], 1.4, || false)
        .unwrap();
    relative(recovery.strain[0], increase * (-rate * 1.4).exp(), 3e-12);
    close(recovery.work, 0.0, 1e-18);
    relative(recovery.dissipated, -recovery.stored_energy_change, 3e-12);
    relative(
        history.work(),
        history.stored_energy() + history.dissipated(),
        3e-12,
    );
}

fn coupled() -> (TensorProny, Tangent6, [Tangent6; 2]) {
    let eq = diagonal([1200.0, 1700.0, 1000.0, 900.0, 700.0, 1300.0]);
    let mut factors = [[[0.0; 6]; 6]; 2];
    factors[0][0] = [17.0, -3.0, 2.0, 7.0, 1.0, 0.0];
    factors[1][0] = [1.0, 12.0, -2.0, 6.0, 3.0, -4.0];
    factors[1][1] = [-3.0, 2.0, 8.0, -4.0, 0.0, 4.0];
    let stiffness = factors.map(|b| {
        core::array::from_fn(|i| core::array::from_fn(|j| (0..6).map(|k| b[k][i] * b[k][j]).sum()))
    });
    let branches = factors
        .into_iter()
        .zip([0.23, 1.1])
        .map(|(b, tau)| TensorPronyBranch::from_gram_rows(b, tau).unwrap())
        .collect();
    (
        TensorProny::new_mandel(eq, branches, 0.1).unwrap(),
        eq,
        stiffness,
    )
}

// Independent partial-pivot elimination for the six original physical strain
// rates. No spectral/compliance preparation or Prony update enters this oracle.
fn solve(mut a: Tangent6, mut rhs: Voigt) -> Voigt {
    for j in 0..6 {
        let pivot = (j..6)
            .max_by(|&i, &k| a[i][j].abs().total_cmp(&a[k][j].abs()))
            .unwrap();
        a.swap(j, pivot);
        rhs.swap(j, pivot);
        let pivot_row = a[j];
        for i in j + 1..6 {
            let factor = a[i][j] / pivot_row[j];
            for (value, pivot) in a[i][j..].iter_mut().zip(&pivot_row[j..]) {
                *value -= factor * pivot;
            }
            rhs[i] -= factor * rhs[j];
        }
    }
    let mut x = [0.0; 6];
    for i in (0..6).rev() {
        x[i] = (rhs[i] - (i + 1..6).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
    }
    x
}

// State = original strain, two original elastic branch strains, work, heat.
fn derivative(
    y: &[f64; 20],
    stress: Voigt,
    slope: Voigt,
    c0: &Tangent6,
    branches: &[Tangent6; 2],
) -> [f64; 20] {
    let mut rhs = slope;
    let mut heat = 0.0;
    for (j, tau) in [0.23, 1.1].into_iter().enumerate() {
        let q = core::array::from_fn(|k| y[6 + 6 * j + k]);
        let sigma = apply(&branches[j], q);
        heat += dot(q, sigma) / tau;
        for k in 0..6 {
            rhs[k] += sigma[k] / tau;
        }
    }
    let rate = solve(*c0, rhs);
    let mut result = [0.0; 20];
    result[..6].copy_from_slice(&rate);
    for (j, tau) in [0.23, 1.1].into_iter().enumerate() {
        for k in 0..6 {
            result[6 + 6 * j + k] = rate[k] - y[6 + 6 * j + k] / tau;
        }
    }
    result[18] = dot(stress, rate);
    result[19] = heat;
    result
}

fn integrate(
    y: &mut [f64; 20],
    start: Voigt,
    target: Voigt,
    dt: f64,
    count: usize,
    c0: &Tangent6,
    branches: &[Tangent6; 2],
) -> [Voigt; 2] {
    let h = dt / count as f64;
    let slope = core::array::from_fn(|k| (target[k] - start[k]) / dt);
    let mut bounds = [core::array::from_fn(|k| y[k]); 2];
    for i in 0..count {
        let stress =
            |offset: f64| core::array::from_fn(|k| start[k] + (i as f64 + offset) * h * slope[k]);
        let k1 = derivative(y, stress(0.0), slope, c0, branches);
        let k2 = derivative(
            &core::array::from_fn(|k| y[k] + 0.5 * h * k1[k]),
            stress(0.5),
            slope,
            c0,
            branches,
        );
        let k3 = derivative(
            &core::array::from_fn(|k| y[k] + 0.5 * h * k2[k]),
            stress(0.5),
            slope,
            c0,
            branches,
        );
        let k4 = derivative(
            &core::array::from_fn(|k| y[k] + h * k3[k]),
            stress(1.0),
            slope,
            c0,
            branches,
        );
        for k in 0..20 {
            y[k] += h * (k1[k] + 2.0 * k2[k] + 2.0 * k3[k] + k4[k]) / 6.0;
        }
        for k in 0..6 {
            bounds[0][k] = bounds[0][k].min(y[k]);
            bounds[1][k] = bounds[1][k].max(y[k]);
        }
    }
    bounds
}

#[test]
fn g1_noncommuting_tensor_creep_matches_original_ode_work_heat_and_path_bounds() {
    let (law, eq, branches) = coupled();
    let c0 = core::array::from_fn(|i| {
        core::array::from_fn(|j| eq[i][j] + branches[0][i][j] + branches[1][i][j])
    });
    // These branch tensors genuinely do not share principal directions.
    let ab = apply(
        &branches[0],
        apply(&branches[1], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    );
    let ba = apply(
        &branches[1],
        apply(&branches[0], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    );
    assert!((ab[0] - ba[0]).abs() + (ab[1] - ba[1]).abs() > 1.0);
    let ramps = [
        (0.14, [1.0, -0.7, 0.2, 0.8, -0.4, 0.3]),
        (0.57, [1.0, -0.7, 0.2, 0.8, -0.4, 0.3]),
        (0.32, [-0.6, 0.5, -0.3, 0.2, 0.4, -0.5]),
        (0.83, [0.0; 6]),
    ];
    let mut history = law.history();
    let control = Controller::new(&history, Budget::default(), || false).unwrap();
    let reports: Vec<_> = ramps
        .iter()
        .map(|&(dt, stress)| {
            control
                .advance(&mut history, tensor(stress), dt, || false)
                .unwrap()
        })
        .collect();
    let mut errors = Vec::new();
    for count in [32, 64, 128] {
        let mut y = [0.0; 20];
        let mut start = [0.0; 6];
        for ((dt, target), report) in ramps.iter().zip(&reports) {
            let bounds = integrate(&mut y, start, *target, *dt, count, &c0, &branches);
            for k in 0..6 {
                assert!(bounds[0][k] / W[k] >= report.strain_bounds[0][k] - 2e-9);
                assert!(bounds[1][k] / W[k] <= report.strain_bounds[1][k] + 2e-9);
            }
            start = *target;
        }
        let strain = mandel(history.strain());
        let error = (0..6)
            .map(|k| (y[k] - strain[k]).abs())
            .fold(0.0_f64, f64::max)
            + (y[18] - history.work()).abs()
            + (y[19] - history.dissipated()).abs();
        errors.push(error);
        if count == 128 {
            assert!(error < 2e-9, "independent ODE discrepancy {error:e}");
            let e = core::array::from_fn(|k| y[k]);
            let stored = 0.5 * dot(e, apply(&eq, e))
                + (0..2)
                    .map(|j| {
                        let q = core::array::from_fn(|k| y[6 + 6 * j + k]);
                        0.5 * dot(q, apply(&branches[j], q))
                    })
                    .sum::<f64>();
            close(stored, history.stored_energy(), 2e-9);
        }
    }
    assert!(
        errors[1] < errors[0] / 8.0 && errors[2] < errors[1] / 8.0,
        "{errors:?}"
    );
    for report in reports {
        assert!(report.dissipated >= 0.0);
        close(report.balance_residual, 0.0, 2e-14);
    }
}

#[test]
fn g3_stress_ramp_subdivision_and_switching_preserve_the_original_memory() {
    let (law, _, _) = coupled();
    let rotation = [[0.6, -0.8, 0.0], [0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];
    let rotate = |v| fs_material::tensor::rotate(&v, &rotation);
    // The actual rotation API retains a NONTRIANGULAR equilibrium Gram factor.
    // Replacing this with a fresh new_mandel constructor would miss that case.
    let mut rotated = law.rotated(rotation).unwrap().history();
    let initial = [0.001, -0.0004, 0.0007, 0.0003, -0.0002, 0.0005];
    let mut whole = law.history();
    whole.advance(initial, 0.4, || false).unwrap();
    rotated.advance(rotate(initial), 0.4, || false).unwrap();
    let control = Controller::new(&whole, Budget::default(), || false).unwrap();
    let rotated_control = Controller::new(&rotated, Budget::default(), || false).unwrap();
    let mut split = whole.clone();
    let target = [-0.5, 0.7, -0.2, 0.3, 0.4, -0.6];
    let middle = core::array::from_fn(|k| 0.5 * whole.stress()[k] + 0.5 * target[k]);
    control.advance(&mut whole, target, 0.73, || false).unwrap();
    control
        .advance(&mut split, middle, 0.365, || false)
        .unwrap();
    control
        .advance(&mut split, target, 0.365, || false)
        .unwrap();
    rotated_control
        .advance(&mut rotated, rotate(target), 0.73, || false)
        .unwrap();
    let expected = rotate(whole.strain());
    for (k, &expected) in expected.iter().enumerate() {
        close(whole.strain()[k], split.strain()[k], 3e-13);
        close(rotated.strain()[k], expected, 3e-13);
    }
    close(whole.work(), split.work(), 3e-13);
    close(whole.dissipated(), split.dissipated(), 3e-13);
    close(whole.work(), rotated.work(), 3e-13);
    close(whole.dissipated(), rotated.dissipated(), 3e-13);
    let strain = [0.0007, 0.0001, -0.0003, 0.0004, -0.0005, 0.0002];
    let a = whole.advance(strain, 0.31, || false).unwrap();
    let b = split.advance(strain, 0.31, || false).unwrap();
    for k in 0..6 {
        close(a.stress[k], b.stress[k], 3e-10);
    }
    let a = control
        .advance(&mut whole, [0.0; 6], 0.21, || false)
        .unwrap();
    let b = control
        .advance(&mut split, [0.0; 6], 0.21, || false)
        .unwrap();
    for k in 0..6 {
        close(a.strain[k], b.strain[k], 3e-13);
    }
}

#[test]
fn g1_weak_equilibrium_saturated_creep_and_sub_ulp_heat_keep_physical_limits() {
    for (e, b, tau, time, stress) in [
        (1e-20, 1.0, 1.0, 1e20, 1e-23),
        (400.0, 1600.0, 1e-12, 2e3, 1.0),
        (400.0, 1600.0, 1e18, 1.0, 1.0),
        (400.0, 1600.0, 1e300, 1.0, 1.0),
    ] {
        let law = sls(e, b, tau, 0.02);
        let mut history = law.history();
        let control = Controller::new(&history, Budget::default(), || false).unwrap();
        let target = [stress, 0.0, 0.0, 0.0, 0.0, 0.0];
        control
            .advance(&mut history, target, 0.0, || false)
            .unwrap();
        let report = control
            .advance(&mut history, target, time, || false)
            .unwrap();
        let r = time * (e / (e + b)) / tau;
        let increase = stress * (b / (e + b)) / e * -(-r).exp_m1();
        let heat = (0.5 * stress * (stress / e)) * (b / (e + b)) * -(-2.0 * r).exp_m1();
        relative(report.strain_increment[0], increase, 2e-10);
        relative(report.dissipated, heat, 2e-10);
        relative(report.work, stress * increase, 2e-10);
    }
}

#[test]
fn g4_budgets_binding_kelvin_strain_and_cancellation_refuse_atomically() {
    let (law, _, _) = coupled();
    let initial = law.history();
    for budget in [
        Budget {
            max_memory_coordinates: 11,
            ..Budget::default()
        },
        Budget {
            max_dense_scalars: 0,
            ..Budget::default()
        },
        Budget {
            max_preparation_work: 0,
            ..Budget::default()
        },
    ] {
        assert!(Controller::new(&initial, budget, || false).is_err());
    }
    assert!(matches!(
        Controller::new(&initial, Budget::default(), || true),
        Err(TensorPronyError::Cancelled)
    ));
    let control = Controller::new(&initial, Budget::default(), || false).unwrap();
    let mut foreign = law.history();
    let before = foreign.clone();
    assert!(
        control
            .advance(&mut foreign, [0.0; 6], 0.1, || false)
            .is_err()
    );
    assert_eq!(foreign, before);
    let viscous = law
        .with_kelvin_viscosity_mandel(diagonal([1.0; 6]))
        .unwrap();
    assert!(Controller::new(&viscous.history(), Budget::default(), || false).is_err());
    let target = [1.0, -0.5, 0.2, 0.4, -0.3, 0.7];
    let mut accepted = initial.clone();
    let mut polls = 0;
    control
        .advance(&mut accepted, target, 0.7, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut history = initial.clone();
        let mut count = 0;
        assert!(matches!(
            control.advance(&mut history, target, 0.7, || {
                count += 1;
                count >= stop
            }),
            Err(TensorPronyError::Cancelled)
        ));
        assert_eq!(history, initial);
    }
    let mut history = initial.clone();
    assert!(
        control
            .advance(&mut history, [f64::NAN; 6], 0.1, || false)
            .is_err()
    );
    assert!(
        control
            .advance(&mut history, target, -0.1, || false)
            .is_err()
    );
    assert_eq!(history, initial);
    let limited = sls(100.0, 900.0, 0.2, 0.005);
    let mut history = limited.history();
    let control = Controller::new(&history, Budget::default(), || false).unwrap();
    let target = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    control
        .advance(&mut history, target, 0.0, || false)
        .unwrap();
    let before = history.clone();
    assert!(
        control
            .advance(&mut history, target, 10.0, || false)
            .is_err()
    );
    assert_eq!(history, before);
}
