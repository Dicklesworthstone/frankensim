//! G1/G3/G4: mixed hereditary loading against physical source equations.
use fs_material::visco_tensor::stress_control::{
    PreparedTensorMixedControl as Mixed, PreparedTensorStressControl as Stress,
    TensorStressControlBudget as Budget,
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
const PLANE: [bool; 6] = [false, false, true, false, true, true];
const TAUS: [f64; 2] = [0.23, 1.1];

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
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "actual={a:e}, expected={b:e}, tolerance={tolerance:e}"
    );
}

fn coupled() -> (TensorProny, Tangent6, [Tangent6; 2]) {
    let mut eq = diagonal([1200.0, 1700.0, 1000.0, 900.0, 700.0, 1300.0]);
    eq[0][2] = -130.0;
    eq[2][0] = -130.0;
    eq[1][4] = 90.0;
    eq[4][1] = 90.0;
    let mut factors = [[[0.0; 6]; 6]; 2];
    factors[0][0] = [17.0, -3.0, 14.0, 7.0, 1.0, -2.0];
    factors[1][0] = [1.0, 12.0, -2.0, 6.0, 9.0, -4.0];
    factors[1][1] = [-3.0, 2.0, 8.0, -4.0, 0.0, 11.0];
    let branches = factors.map(|b| {
        core::array::from_fn(|i| core::array::from_fn(|j| (0..6).map(|k| b[k][i] * b[k][j]).sum()))
    });
    let arms = factors
        .into_iter()
        .zip(TAUS)
        .map(|(b, tau)| TensorPronyBranch::from_gram_rows(b, tau).unwrap())
        .collect();
    (
        TensorProny::new_mandel(eq, arms, 0.1).unwrap(),
        eq,
        branches,
    )
}

#[test]
fn g1_plane_stress_creep_keeps_three_dimensional_memory_and_reaction_work() {
    let (ex, ez, a, b, tau, strain) = (1200.0, 600.0, 30.0, 24.0, 0.4, 0.0015);
    let mut factor = [[0.0; 6]; 6];
    factor[0][0] = a;
    factor[0][2] = b;
    let law = TensorProny::new_mandel(
        diagonal([ex, 900.0, ez, 800.0, 700.0, 1100.0]),
        vec![TensorPronyBranch::from_gram_rows(factor, tau).unwrap()],
        0.1,
    )
    .unwrap();
    let mut history = law.history();
    let control = Mixed::new(&history, PLANE, Budget::default(), || false).unwrap();
    assert_eq!(control.stress_components(), PLANE);
    assert_eq!(control.memory_coordinates(), 6);
    let target = [strain, 0.0, 0.0, 0.0, 0.0, 0.0];
    let jump = control
        .advance(&mut history, target, 0.0, || false)
        .unwrap();
    let gain = a * ez / (ez + b * b);
    let rate = ez / ((ez + b * b) * tau);
    let q0 = gain * strain;
    close(jump.strain[2], -b * q0 / ez, 2e-15);
    close(jump.stress[0], ex * strain + a * q0, 2e-12);
    close(jump.work, 0.5 * jump.stress[0] * strain, 2e-14);
    assert_eq!(jump.dissipated, 0.0);

    let dt = 0.73;
    let hold = control.advance(&mut history, target, dt, || false).unwrap();
    let q = q0 * (-rate * dt).exp();
    close(hold.strain[2], -b * q / ez, 2e-15);
    close(hold.stress[0], ex * strain + a * q, 2e-12);
    close(hold.work, 0.0, 2e-15);
    close(
        hold.dissipated,
        q0 * q0 / (2.0 * rate * tau) * -(-2.0 * rate * dt).exp_m1(),
        2e-14,
    );
    assert!((hold.strain[2] - jump.strain[2]).abs() > 1e-4);

    // The prescribed xx strain does work against the evolving reaction. The
    // free zz surface does zero work while its physical strain keeps changing.
    let next_strain = -0.0004;
    let duration = 0.31;
    let delta = next_strain - strain;
    let r = rate * duration;
    let phi = -(-r).exp_m1() / r;
    let q1 = (-r).exp() * q + gain * delta * phi;
    let mean_q = q * phi + gain * delta * (1.0 - phi) / r;
    let work = 0.5 * ex * (next_strain * next_strain - strain * strain) + a * delta * mean_q;
    let ramp = control
        .advance(
            &mut history,
            [next_strain, 0.0, 0.0, 0.0, 0.0, 0.0],
            duration,
            || false,
        )
        .unwrap();
    close(ramp.strain[2], -b * q1 / ez, 2e-15);
    close(ramp.stress[0], ex * next_strain + a * q1, 2e-12);
    close(ramp.work, work, 2e-14);
    for k in [2, 4, 5] {
        assert_eq!(ramp.stress[k], 0.0);
    }
    close(ramp.balance_residual, 0.0, 2e-14);
    let release = control
        .advance(&mut history, [0.0; 6], 0.0, || false)
        .unwrap();
    let remaining = q1 - gain * next_strain;
    close(release.strain[2], -b * remaining / ez, 2e-15);
    assert!(
        release.strain[2].abs() > 1e-5,
        "original memory must survive unloading"
    );
    let recover = control
        .advance(&mut history, [0.0; 6], 1.3, || false)
        .unwrap();
    close(
        recover.strain[2],
        release.strain[2] * (-rate * 1.3).exp(),
        2e-15,
    );
    close(recover.work, 0.0, 2e-15);
    close(recover.dissipated, -recover.stored_energy_change, 2e-14);
    close(
        history.work(),
        history.stored_energy() + history.dissipated(),
        2e-14,
    );
}

// Independent partial-pivot elimination. Stress-controlled rows use C0;
// strain-controlled rows are identity constraints on the original rates.
fn solve(mut a: Tangent6, mut rhs: Voigt) -> Voigt {
    for j in 0..6 {
        let pivot = (j..6)
            .max_by(|&i, &k| a[i][j].abs().total_cmp(&a[k][j].abs()))
            .unwrap();
        a.swap(j, pivot);
        rhs.swap(j, pivot);
        let row = a[j];
        for i in j + 1..6 {
            let factor = a[i][j] / row[j];
            for k in j..6 {
                a[i][k] -= factor * row[k];
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

// State: six strains, two ORIGINAL branch elastic strains, work and heat.
// Neither native spectral modes nor Prony scalar responses enter this oracle.
struct OriginalOde {
    eq: Tangent6,
    branches: [Tangent6; 2],
    c0: Tangent6,
    mask: [bool; 6],
}
impl OriginalOde {
    fn stress(&self, y: &[f64; 20]) -> Voigt {
        let mut stress = apply(&self.eq, core::array::from_fn(|k| y[k]));
        for (j, branch) in self.branches.iter().enumerate() {
            let arm = apply(branch, core::array::from_fn(|k| y[6 + 6 * j + k]));
            for k in 0..6 {
                stress[k] += arm[k];
            }
        }
        stress
    }
    fn derivative(&self, y: &[f64; 20], slope: Voigt) -> [f64; 20] {
        let mut rhs = slope;
        let mut heat = 0.0;
        for (j, tau) in TAUS.into_iter().enumerate() {
            let q = core::array::from_fn(|k| y[6 + 6 * j + k]);
            let stress = apply(&self.branches[j], q);
            heat += dot(q, stress) / tau;
            for k in 0..6 {
                if self.mask[k] {
                    rhs[k] += stress[k] / tau;
                }
            }
        }
        let mut matrix = self.c0;
        for k in 0..6 {
            if !self.mask[k] {
                matrix[k] = [0.0; 6];
                matrix[k][k] = 1.0;
            }
        }
        let rate = solve(matrix, rhs);
        let mut result = [0.0; 20];
        result[..6].copy_from_slice(&rate);
        for (j, tau) in TAUS.into_iter().enumerate() {
            for k in 0..6 {
                result[6 + 6 * j + k] = rate[k] - y[6 + 6 * j + k] / tau;
            }
        }
        result[18] = dot(self.stress(y), rate);
        result[19] = heat;
        result
    }
    fn integrate(
        &self,
        y: &mut [f64; 20],
        start: Voigt,
        target: Voigt,
        dt: f64,
        count: usize,
    ) -> [Voigt; 2] {
        let h = dt / count as f64;
        let slope = core::array::from_fn(|k| (target[k] - start[k]) / dt);
        let mut bounds = [core::array::from_fn(|k| y[k]); 2];
        for _ in 0..count {
            let k1 = self.derivative(y, slope);
            let k2 = self.derivative(&core::array::from_fn(|k| y[k] + 0.5 * h * k1[k]), slope);
            let k3 = self.derivative(&core::array::from_fn(|k| y[k] + 0.5 * h * k2[k]), slope);
            let k4 = self.derivative(&core::array::from_fn(|k| y[k] + h * k3[k]), slope);
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
}

#[test]
fn g1_noncommuting_mixed_history_matches_original_ode_strain_stress_work_and_heat() {
    let (law, eq, branches) = coupled();
    let mask = [true, false, true, false, false, true];
    let c0 = core::array::from_fn(|i| {
        core::array::from_fn(|j| eq[i][j] + branches[0][i][j] + branches[1][i][j])
    });
    let ode = OriginalOde {
        eq,
        branches,
        c0,
        mask,
    };
    let seed = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let ab = apply(&branches[0], apply(&branches[1], seed));
    let ba = apply(&branches[1], apply(&branches[0], seed));
    assert!((0..6).map(|k| (ab[k] - ba[k]).abs()).sum::<f64>() > 1.0);
    let ramps = [
        (0.13, [0.8, 0.0006, -0.25, 0.0003, -0.0002, 0.4]),
        (0.61, [0.8, 0.0006, -0.25, 0.0003, -0.0002, 0.4]),
        (0.29, [-0.3, -0.0004, 0.6, -0.0001, 0.0005, -0.2]),
        (0.77, [0.0; 6]),
    ];
    let mut history = law.history();
    let control = Mixed::new(&history, mask, Budget::default(), || false).unwrap();
    let reports: Vec<_> = ramps
        .iter()
        .map(|&(dt, target)| {
            control
                .advance(&mut history, tensor(target), dt, || false)
                .unwrap()
        })
        .collect();
    let mut errors = Vec::new();
    for count in [16, 32, 64] {
        let mut y = [0.0; 20];
        let mut start = [0.0; 6];
        for ((dt, target), report) in ramps.iter().zip(&reports) {
            let bounds = ode.integrate(&mut y, start, *target, *dt, count);
            for k in 0..6 {
                assert!(bounds[0][k] / W[k] >= report.strain_bounds[0][k] - 2e-8);
                assert!(bounds[1][k] / W[k] <= report.strain_bounds[1][k] + 2e-8);
                if mask[k] {
                    close(mandel(report.stress)[k], target[k], 2e-14);
                } else {
                    close(mandel(report.strain)[k], target[k], 2e-16);
                }
            }
            start = *target;
        }
        let strain = mandel(history.strain());
        let stress = mandel(history.stress());
        let source_stress = ode.stress(&y);
        let error = (0..6)
            .map(|k| (y[k] - strain[k]).abs() + (source_stress[k] - stress[k]).abs() / 1000.0)
            .sum::<f64>()
            + (y[18] - history.work()).abs()
            + (y[19] - history.dissipated()).abs();
        errors.push(error);
        if count == 64 {
            assert!(error < 2e-9, "independent ODE discrepancy {error:e}");
        }
    }
    assert!(
        errors[1] < errors[0] / 8.0 && errors[2] < errors[1] / 8.0,
        "{errors:?}"
    );
    for report in reports {
        assert!(report.dissipated >= 0.0);
        close(report.balance_residual, 0.0, 2e-13);
    }
    close(
        history.work(),
        history.stored_energy() + history.dissipated(),
        2e-13,
    );
}

#[test]
fn g3_all_stress_and_all_strain_partitions_match_the_existing_history_owners() {
    let (law, _, _) = coupled();
    let mut initial = law.history();
    initial
        .advance(
            [0.001, -0.0004, 0.0007, 0.0003, -0.0002, 0.0005],
            0.4,
            || false,
        )
        .unwrap();
    let stress = Stress::new(&initial, Budget::default(), || false).unwrap();
    let all_stress = Mixed::new(&initial, [true; 6], Budget::default(), || false).unwrap();
    let mut native = initial.clone();
    let mut mixed = initial.clone();
    let target = [-0.5, 0.7, -0.2, 0.3, 0.4, -0.6];
    let a = stress.advance(&mut native, target, 0.73, || false).unwrap();
    let b = all_stress
        .advance(&mut mixed, target, 0.73, || false)
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(native, mixed);
    let all_strain = Mixed::new(&mixed, [false; 6], Budget::default(), || false).unwrap();
    for (dt, strain) in [
        (0.23, [0.0007, 0.0001, -0.0003, 0.0004, -0.0005, 0.0002]),
        (0.0, [-0.0004, 0.0008, 0.0002, -0.0003, 0.0001, 0.0006]),
        (0.38, [-0.0004, 0.0008, 0.0002, -0.0003, 0.0001, 0.0006]),
    ] {
        let a = native.advance(strain, dt, || false).unwrap();
        let b = all_strain
            .advance(&mut mixed, strain, dt, || false)
            .unwrap();
        for k in 0..6 {
            close(a.stress[k], b.stress[k], 3e-11);
            close(native.strain()[k], mixed.strain()[k], 3e-15);
        }
        close(a.stored_energy, b.stored_energy, 3e-14);
        close(a.work, b.work, 3e-14);
        close(a.dissipated, b.dissipated, 3e-14);
        close(
            native.branch_stored_energy(),
            mixed.branch_stored_energy(),
            3e-14,
        );
    }
}

#[test]
fn g3_partition_switching_and_split_ramps_preserve_the_original_branch_memory() {
    let (law, _, _) = coupled();
    let mut whole = law.history();
    whole
        .advance(
            [0.001, -0.0004, 0.0007, 0.0003, -0.0002, 0.0005],
            0.4,
            || false,
        )
        .unwrap();
    let mut split = whole.clone();
    for (mask, dt, target) in [
        (PLANE, 0.73, [-0.0005, 0.0007, 0.0, 0.0003, 0.0, 0.0]),
        (
            [true, false, true, false, true, false],
            0.22,
            [0.6, -0.0002, -0.3, 0.0001, 0.2, -0.0004],
        ),
    ] {
        let control = Mixed::new(&whole, mask, Budget::default(), || false).unwrap();
        let start = core::array::from_fn::<_, 6, _>(|k| {
            if mask[k] {
                whole.stress()[k]
            } else {
                whole.strain()[k]
            }
        });
        let middle = core::array::from_fn(|k| 0.5 * start[k] + 0.5 * target[k]);
        control.advance(&mut whole, target, dt, || false).unwrap();
        control
            .advance(&mut split, middle, dt / 2.0, || false)
            .unwrap();
        control
            .advance(&mut split, target, dt / 2.0, || false)
            .unwrap();
        for k in 0..6 {
            close(whole.strain()[k], split.strain()[k], 3e-14);
            close(whole.stress()[k], split.stress()[k], 3e-11);
        }
        close(whole.work(), split.work(), 3e-14);
        close(whole.dissipated(), split.dissipated(), 3e-14);
    }
    let strain = [0.0003, 0.0005, -0.0004, 0.0001, -0.0002, 0.0007];
    let a = whole.advance(strain, 0.31, || false).unwrap();
    let b = split.advance(strain, 0.31, || false).unwrap();
    for k in 0..6 {
        close(a.stress[k], b.stress[k], 3e-11);
    }
    close(a.dissipated, b.dissipated, 3e-14);
}

#[test]
fn g4_mixed_control_budgets_binding_strain_and_cancellation_refuse_atomically() {
    let (law, _, _) = coupled();
    let initial = law.history();
    for budget in [
        Budget {
            max_memory_coordinates: 11,
            ..Budget::default()
        },
        Budget {
            max_preparation_work: 0,
            ..Budget::default()
        },
    ] {
        assert!(Mixed::new(&initial, PLANE, budget, || false).is_err());
    }
    assert!(matches!(
        Mixed::new(&initial, PLANE, Budget::default(), || true),
        Err(TensorPronyError::Cancelled)
    ));
    let control = Mixed::new(&initial, PLANE, Budget::default(), || false).unwrap();
    let mut foreign = law.history();
    let before = foreign.clone();
    assert!(
        control
            .advance(&mut foreign, [0.0; 6], 0.1, || false)
            .is_err()
    );
    assert_eq!(foreign, before);
    let target = [0.001, -0.0003, 0.0, 0.0002, 0.0, 0.0];
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
    for (value, dt) in [
        ([f64::NAN; 6], 0.1),
        (target, -0.1),
        ([0.2, 0.0, 0.0, 0.0, 0.0, 0.0], 0.1),
    ] {
        assert!(control.advance(&mut history, value, dt, || false).is_err());
        assert_eq!(history, initial);
    }
    let viscous = law
        .with_kelvin_viscosity_mandel(diagonal([1.0; 6]))
        .unwrap();
    assert!(Mixed::new(&viscous.history(), PLANE, Budget::default(), || false).is_err());
}
