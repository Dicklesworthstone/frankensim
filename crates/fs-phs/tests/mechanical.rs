//! Independent physical-variable oracles for the original M/C/K adapter.
use fs_phs::RelaxationBranch;
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_phs::mechanical::LinearMechanicalSystem;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance={tolerance:e}"
    );
}
fn budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 32,
        max_omega_dt: 100.0,
    }
}

#[test]
fn dimensional_damped_oscillator_matches_closed_motion_and_integrated_heat() {
    let model = LinearMechanicalSystem::new(&[2.0], &[3.0], &[8.0], &[1.0], 1, 1, vec![]).unwrap();
    let time = 0.7;
    let force = 1.0;
    let initial = model.state(&[0.4], &[-0.3], &[]).unwrap();
    let flow = PreparedExactZoh::new(model.system(), time, budget()).unwrap();
    let step = flow.step(&initial, &[force]).unwrap();
    let closed = |t: f64| {
        let decay = 0.75_f64;
        let frequency = (4.0 - decay * decay).sqrt();
        let a = 0.4 - force / 8.0;
        let b = (-0.3 + decay * a) / frequency;
        let (sin, cos) = (frequency * t).sin_cos();
        let scale = (-decay * t).exp();
        (
            force / 8.0 + scale * (a * cos + b * sin),
            scale * ((-decay * a + frequency * b) * cos + (-decay * b - frequency * a) * sin),
        )
    };
    let (q, v) = closed(time);
    close(step.x[0], q, 3e-13);
    close(model.velocities(&step.x).unwrap()[0], v, 3e-13);
    let mut heat = 0.0;
    for i in 0..=4096 {
        let weight = if i == 0 || i == 4096 {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        heat += weight * 3.0 * closed(time * i as f64 / 4096.0).1.powi(2);
    }
    heat *= time / (3.0 * 4096.0);
    close(step.dissipated, heat, 2e-12);
    close(step.supplied, force * (q - 0.4), 3e-13);
    close(step.balance_residual(), 0.0, 2e-13);
}

fn arms() -> Vec<RelaxationBranch> {
    vec![
        RelaxationBranch {
            projection: vec![1.0, -0.6],
            stiffness: 3.0,
            relaxation_time_s: 0.08,
        },
        RelaxationBranch {
            projection: vec![0.3, 1.0],
            stiffness: 2.0,
            relaxation_time_s: 0.6,
        },
    ]
}
fn coupled() -> LinearMechanicalSystem {
    LinearMechanicalSystem::new(
        &[2.0, 0.4, 0.4, 1.5],
        &[0.8, 1.55, -1.05, 0.6],
        &[7.0, -1.0, -1.0, 4.0],
        &[1.0, 0.2, -0.3, 0.8],
        2,
        2,
        arms(),
    )
    .unwrap()
}

// Original [q0,q1,v0,v1,s0,s1,Qvisc,Qmemory,W] equations. This reference
// never builds pHS matrices, modal coordinates, exponentials or Gramian heat.
fn original_derivative(x: [f64; 9]) -> [f64; 9] {
    let force = [0.62, -0.53]; // P * [0.7,-0.4]
    let r0 = force[0] - (7.0 * x[0] - x[1]) - (0.8 * x[2] + 1.55 * x[3]) - x[4] - 0.3 * x[5];
    let r1 = force[1] - (-x[0] + 4.0 * x[1]) - (-1.05 * x[2] + 0.6 * x[3]) + 0.6 * x[4] - x[5];
    [
        x[2],
        x[3],
        (1.5 * r0 - 0.4 * r1) / 2.84,
        (-0.4 * r0 + 2.0 * r1) / 2.84,
        3.0 * (x[2] - 0.6 * x[3]) - x[4] / 0.08,
        2.0 * (0.3 * x[2] + x[3]) - x[5] / 0.6,
        0.8 * x[2] * x[2] + 0.5 * x[2] * x[3] + 0.6 * x[3] * x[3],
        x[4] * x[4] / (3.0 * 0.08) + x[5] * x[5] / (2.0 * 0.6),
        force[0] * x[2] + force[1] * x[3],
    ]
}
fn reference(time: f64) -> [f64; 9] {
    let steps = 16000;
    let h = time / steps as f64;
    let mut x = [0.2, -0.15, 0.1, 0.25, 0.72, 0.02, 0.0, 0.0, 0.0];
    for _ in 0..steps {
        let a = original_derivative(x);
        let b = original_derivative(core::array::from_fn(|i| x[i] + 0.5 * h * a[i]));
        let c = original_derivative(core::array::from_fn(|i| x[i] + 0.5 * h * b[i]));
        let d = original_derivative(core::array::from_fn(|i| x[i] + h * c[i]));
        for i in 0..9 {
            x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

#[test]
fn consistent_mass_gyroscopic_and_maxwell_forced_motion_keep_original_power() {
    let model = coupled();
    let initial = model
        .state(&[0.2, -0.15], &[0.1, 0.25], &[0.05, -0.1])
        .unwrap();
    assert_eq!(model.coordinates(), 2);
    assert_eq!(model.system().state_dim(), 6);
    assert_eq!(model.branches().len(), 2);
    let stress = model.branch_stresses(&initial).unwrap();
    close(stress[0], 0.72, 2e-15);
    close(stress[1], 0.02, 2e-15);
    let y = model.system().output(&initial);
    close(y[0], 0.1 - 0.3 * 0.25, 2e-15);
    close(y[1], 0.2 * 0.1 + 0.8 * 0.25, 2e-15);
    let time = 0.8;
    let truth = reference(time);
    let whole = PreparedExactZoh::new(model.system(), time, budget())
        .unwrap()
        .step(&initial, &[0.7, -0.4])
        .unwrap();
    let velocity = model.velocities(&whole.x).unwrap();
    let stress = model.branch_stresses(&whole.x).unwrap();
    let actual = [
        whole.x[0],
        whole.x[1],
        velocity[0],
        velocity[1],
        stress[0],
        stress[1],
    ];
    for i in 0..6 {
        close(actual[i], truth[i], 2e-11);
    }
    close(whole.dissipated, truth[6] + truth[7], 2e-11);
    close(whole.supplied, truth[8], 2e-11);
    close(whole.balance_residual(), 0.0, 2e-12);
    assert!(truth[6] > 0.01 && truth[7] > 0.01);

    // The same continuing material state and force through the warm cache.
    let mut part = PreparedExactZoh::new(model.system(), time / 8.0, budget()).unwrap();
    let mut x = initial;
    let mut next = vec![0.0; 6];
    let mut output = [0.0; 2];
    let mut heat = 0.0;
    let mut work = 0.0;
    for _ in 0..8 {
        let step = part
            .step_into(&x, &[0.7, -0.4], &mut next, &mut output)
            .unwrap();
        heat += step.dissipated;
        work += step.supplied;
        core::mem::swap(&mut x, &mut next);
    }
    for i in 0..6 {
        close(x[i], whole.x[i], 2e-12);
    }
    close(heat, whole.dissipated, 2e-12);
    close(work, whole.supplied, 2e-12);
    let relaxed = model.relaxed_state(&[0.2, -0.15], &[0.1, 0.25]).unwrap();
    assert!(
        model
            .branch_stresses(&relaxed)
            .unwrap()
            .iter()
            .all(|s| *s == 0.0)
    );
}

#[test]
fn free_consistent_mass_and_gyroscopic_motion_have_no_fictitious_heat() {
    let model = LinearMechanicalSystem::new(
        &[2.0, 0.4, 0.4, 1.5],
        &[0.0; 4],
        &[0.0; 4],
        &[1.0, 2.0],
        2,
        1,
        vec![],
    )
    .unwrap();
    let initial = model.state(&[0.4, -0.2], &[0.3, 0.1], &[]).unwrap();
    let time = 0.7;
    let step = PreparedExactZoh::new(model.system(), time, budget())
        .unwrap()
        .step(&initial, &[0.6])
        .unwrap();
    let acceleration = [
        (1.5 * 0.6 - 0.4 * 1.2) / 2.84,
        (-0.4 * 0.6 + 2.0 * 1.2) / 2.84,
    ];
    let velocity = model.velocities(&step.x).unwrap();
    for i in 0..2 {
        close(
            step.x[i],
            [0.4, -0.2][i] + time * [0.3, 0.1][i] + 0.5 * time * time * acceleration[i],
            2e-13,
        );
        close(velocity[i], [0.3, 0.1][i] + time * acceleration[i], 2e-13);
    }
    assert_eq!(step.dissipated, 0.0);

    let gyro = LinearMechanicalSystem::new(
        &[1.0, 0.0, 0.0, 1.0],
        &[0.0, 3.0, -3.0, 0.0],
        &[0.0; 4],
        &[],
        2,
        0,
        vec![],
    )
    .unwrap();
    let initial = gyro.state(&[0.0; 2], &[1.0, 0.0], &[]).unwrap();
    let time = 0.4;
    let step = PreparedExactZoh::new(gyro.system(), time, budget())
        .unwrap()
        .step(&initial, &[])
        .unwrap();
    let (sin, cos) = (3.0_f64 * time).sin_cos();
    let velocity = gyro.velocities(&step.x).unwrap();
    close(velocity[0], cos, 2e-13);
    close(velocity[1], sin, 2e-13);
    close(step.x[0], sin / 3.0, 2e-13);
    close(step.x[1], (1.0 - cos) / 3.0, 2e-13);
    assert_eq!(step.dissipated, 0.0);
    close(step.delta_h, 0.0, 2e-13);
}

#[test]
fn hidden_memory_retains_its_independent_analytic_relaxation_heat() {
    let arm = RelaxationBranch {
        projection: vec![0.0],
        stiffness: 2.0,
        relaxation_time_s: 0.3,
    };
    let model = LinearMechanicalSystem::new(&[2.0], &[0.0], &[0.0], &[], 1, 0, vec![arm]).unwrap();
    let initial = model.state(&[0.4], &[0.3], &[1.0]).unwrap();
    let time = 0.7;
    let step = PreparedExactZoh::new(model.system(), time, budget())
        .unwrap()
        .step(&initial, &[])
        .unwrap();
    close(step.x[0], 0.4 + 0.3 * time, 2e-13);
    close(model.velocities(&step.x).unwrap()[0], 0.3, 2e-13);
    close(step.x[2], 2.0_f64.sqrt() * (-time / 0.3).exp(), 2e-13);
    close(step.dissipated, -(-2.0 * time / 0.3).exp_m1(), 2e-13);
    close(step.balance_residual(), 0.0, 2e-13);
}

#[test]
fn nonphysical_matrices_and_incomplete_original_strains_refuse() {
    for (mass, damping, stiffness) in [(-1.0, 1.0, 1.0), (1.0, -1.0, 1.0), (1.0, 1.0, -1.0)] {
        assert!(
            LinearMechanicalSystem::new(&[mass], &[damping], &[stiffness], &[1.0], 1, 1, vec![])
                .is_err()
        );
    }
    let wrong = RelaxationBranch {
        projection: vec![1.0],
        stiffness: 2.0,
        relaxation_time_s: 0.3,
    };
    assert!(
        LinearMechanicalSystem::new(
            &[1.0, 0.0, 0.0, 1.0],
            &[0.0; 4],
            &[0.0; 4],
            &[],
            2,
            0,
            vec![wrong]
        )
        .is_err()
    );
    assert!(
        LinearMechanicalSystem::new(
            &[1.0, 0.2, 0.0, 1.0],
            &[0.0; 4],
            &[0.0; 4],
            &[],
            2,
            0,
            vec![]
        )
        .is_err()
    );
    let model = coupled();
    assert!(model.state(&[0.0; 2], &[0.0; 2], &[0.0]).is_err());
    assert!(model.velocities(&[0.0; 4]).is_err());
}
