//! Actual transient rotor motion: physical impulses, whirl and bearing work.
use fs_ga::{Mat3, So3, Vec3};
use fs_mbd::articulated::SpatialInertia;
use fs_mbd::rotordynamics::rigid::RigidRotor;
use fs_mbd::rotordynamics::transient::RotorTransientState;
use fs_mbd::rotordynamics::{RotorBearingModel, RotorError};
use std::cell::Cell;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn dot_action(matrix: &[f64], left: &[f64], right: &[f64]) -> f64 {
    let n = left.len();
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| left[i] * matrix[i * n + j] * right[j])
                .sum::<f64>()
        })
        .sum()
}

#[test]
fn full_mass_impulse_then_forced_motion_matches_independent_solution() {
    let model = RotorBearingModel::new(
        2,
        &[2.0, 0.5, 0.5, 1.0],
        &[0.0; 4],
        &[0.0; 4],
        &[0.0; 4],
        || false,
    )
    .unwrap();
    let integrator = model
        .at_speed(9.0, || false)
        .unwrap()
        .midpoint_stepper(0.125, 1e-12, || false)
        .unwrap();
    let original = RotorTransientState {
        displacement: vec![1.0, -2.0],
        velocity: vec![0.3, -0.7],
    };
    let kick = integrator
        .apply_impulse(&original, &[1.0, 0.0], || false)
        .unwrap();
    assert_eq!(kick.state_after.displacement, original.displacement);
    near(kick.state_after.velocity[0], 0.3 + 4.0 / 7.0, 2e-14);
    near(kick.state_after.velocity[1], -0.7 - 2.0 / 7.0, 2e-14);
    near(kick.impulse_work_j, 0.3 + 2.0 / 7.0, 2e-14);
    near(kick.kinetic_energy_change_j, kick.impulse_work_j, 2e-14);

    let result = integrator
        .step(&kick.state_after, &[0.2, -0.4], || false)
        .unwrap();
    let h = integrator.duration_seconds();
    // M^-1 [0.2,-0.4] = [8/35,-18/35], without mass lumping.
    for (i, acceleration) in [8.0 / 35.0, -18.0 / 35.0].into_iter().enumerate() {
        near(
            result.state_after.velocity[i],
            kick.state_after.velocity[i] + h * acceleration,
            2e-14,
        );
        near(
            result.state_after.displacement[i],
            original.displacement[i]
                + h * kick.state_after.velocity[i]
                + 0.5 * h * h * acceleration,
            2e-14,
        );
    }
    near(result.kinetic_energy_change_j, result.input_work_j, 2e-14);
    near(result.work_energy_residual_j, 0.0, 2e-14);
    assert_eq!(
        integrator
            .step(&kick.state_after, &[0.2, -0.4], || false)
            .unwrap(),
        result
    );
}

fn free_rotor() -> RigidRotor {
    let inertia = SpatialInertia::new(
        4.0,
        Vec3::default(),
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0],
        },
    )
    .unwrap();
    RigidRotor::new(inertia, So3::identity(), &[], 1e-12, || false).unwrap()
}

fn whirl(spin: f64, steps: usize) -> RotorTransientState {
    let rotor = free_rotor();
    let h = 1.0 / steps as f64;
    let integrator = rotor
        .model()
        .at_speed(spin, || false)
        .unwrap()
        .midpoint_stepper(h, 1e-11, || false)
        .unwrap();
    let mut state = RotorTransientState {
        displacement: vec![0.0; 4],
        velocity: vec![0.0, 0.0, 0.1, 0.05],
    };
    let initial_energy = 0.5 * dot_action(rotor.model().mass(), &state.velocity, &state.velocity);
    for _ in 0..steps {
        let result = integrator.step(&state, &[0.0; 4], || false).unwrap();
        near(result.gyroscopic_work_j, 0.0, 2e-15);
        near(result.work_energy_residual_j, 0.0, 2e-15);
        state = result.state_after;
        // Physical transverse H = Jd theta_dot + Jp Omega [theta_y,-theta_x].
        near(
            2.0 * state.velocity[2] + 3.0 * spin * state.displacement[3],
            0.2,
            5e-14,
        );
        near(
            2.0 * state.velocity[3] - 3.0 * spin * state.displacement[2],
            0.1,
            5e-14,
        );
        near(
            0.5 * dot_action(rotor.model().mass(), &state.velocity, &state.velocity),
            initial_energy,
            5e-14,
        );
    }
    state
}

#[test]
fn native_rotor_whirl_preserves_angular_momentum_and_refines_at_second_order() {
    for spin in [-2.0, 2.0] {
        let frequency = 1.5_f64 * spin;
        let (sin, cos) = frequency.sin_cos();
        let exact = [
            (0.1 * sin + 0.05 * (cos - 1.0)) / frequency,
            (0.1 * (1.0 - cos) + 0.05 * sin) / frequency,
            0.1 * cos - 0.05 * sin,
            0.1 * sin + 0.05 * cos,
        ];
        let error = |steps| {
            let state = whirl(spin, steps);
            [
                state.displacement[2],
                state.displacement[3],
                state.velocity[2],
                state.velocity[3],
            ]
            .iter()
            .zip(exact)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
        };
        let coarse = error(20);
        let medium = error(40);
        let fine = error(80);
        assert!(coarse / medium > 3.95, "{coarse:e} / {medium:e}");
        assert!(medium / fine > 3.98, "{medium:e} / {fine:e}");
    }
}

fn damped_motion(steps: usize) -> RotorTransientState {
    let model = RotorBearingModel::new(1, &[2.0], &[0.4], &[8.0], &[0.0], || false).unwrap();
    let integrator = model
        .at_speed(12.0, || false)
        .unwrap()
        .midpoint_stepper(3.0 / steps as f64, 1e-11, || false)
        .unwrap();
    let mut state = RotorTransientState {
        displacement: vec![1.0],
        velocity: vec![0.0],
    };
    for _ in 0..steps {
        let before = state.velocity[0].powi(2) + 4.0 * state.displacement[0].powi(2);
        let result = integrator.step(&state, &[0.0], || false).unwrap();
        let next = &result.state_after;
        let after = next.velocity[0].powi(2) + 4.0 * next.displacement[0].powi(2);
        assert!(result.velocity_force_work_j >= 0.0);
        assert!(after <= before + 1e-13);
        near(after - before, -result.velocity_force_work_j, 2e-13);
        near(result.work_energy_residual_j, 0.0, 2e-13);
        state = result.state_after;
    }
    state
}

#[test]
fn passive_ring_down_has_physical_bearing_loss_and_analytic_refinement() {
    let decay = 0.1_f64;
    let frequency = (4.0 - decay * decay).sqrt();
    let (sin, cos) = (3.0 * frequency).sin_cos();
    let envelope = (-3.0 * decay).exp();
    let exact = [
        envelope * (cos + decay * sin / frequency),
        -4.0 * envelope * sin / frequency,
    ];
    let error = |steps| {
        let state = damped_motion(steps);
        (state.displacement[0] - exact[0]).hypot(state.velocity[0] - exact[1])
    };
    let coarse = error(150);
    let fine = error(300);
    assert!(coarse / fine > 3.95);
    assert!(fine < 0.0004, "{fine:e}");
}

#[test]
fn nonsymmetric_active_bearings_keep_their_original_force_and_signed_work() {
    let mass = [2.0, 0.2, 0.2, 1.0];
    let damping = [-1.0, 0.4, -0.2, -0.3];
    let stiffness = [2.0, 1.5, -0.5, 3.0];
    let gyro = [0.0, 0.4, -0.4, 0.0];
    let model = RotorBearingModel::new(2, &mass, &damping, &stiffness, &gyro, || false).unwrap();
    let spin = 7.0;
    let h = 0.02;
    let integrator = model
        .at_speed(spin, || false)
        .unwrap()
        .midpoint_stepper(h, 1e-12, || false)
        .unwrap();
    let state = RotorTransientState {
        displacement: vec![0.4, -0.2],
        velocity: vec![1.0, 0.3],
    };
    let force = [0.2, -0.1];
    let result = integrator.step(&state, &force, || false).unwrap();
    let q: Vec<_> = state
        .displacement
        .iter()
        .zip(&result.state_after.displacement)
        .map(|(a, b)| 0.5 * (a + b))
        .collect();
    let v: Vec<_> = state
        .velocity
        .iter()
        .zip(&result.state_after.velocity)
        .map(|(a, b)| 0.5 * (a + b))
        .collect();
    for i in 0..2 {
        let impulse: f64 = (0..2)
            .map(|j| {
                mass[2 * i + j] * (result.state_after.velocity[j] - state.velocity[j])
                    + h * (damping[2 * i + j] + spin * gyro[2 * i + j]) * v[j]
                    + h * stiffness[2 * i + j] * q[j]
            })
            .sum();
        near(impulse, h * force[i], 2e-14);
    }
    near(
        result.velocity_force_work_j,
        h * dot_action(&damping, &v, &v),
        2e-14,
    );
    near(
        result.stiffness_force_work_j,
        h * dot_action(&stiffness, &v, &q),
        2e-14,
    );
    assert!(
        result.velocity_force_work_j < 0.0,
        "active bearing work must remain signed"
    );
    near(result.gyroscopic_work_j, 0.0, 2e-14);
    near(result.work_energy_residual_j, 0.0, 2e-14);
}

#[test]
fn singular_unresolved_and_cancelled_motion_publish_no_state() {
    let unstable = RotorBearingModel::new(1, &[1.0], &[-2.0], &[0.0], &[0.0], || false).unwrap();
    assert!(matches!(
        unstable
            .at_speed(0.0, || false)
            .unwrap()
            .midpoint_stepper(1.0, 1e-12, || false),
        Err(RotorError::TransientLinear(_))
    ));
    let model = RotorBearingModel::new(1, &[1.0], &[0.0], &[1.0], &[0.0], || false).unwrap();
    let point = model.at_speed(0.0, || false).unwrap();
    assert!(point.midpoint_stepper(0.0, 1e-12, || false).is_err());
    let integrator = point.midpoint_stepper(0.1, 1e-12, || false).unwrap();
    let original = RotorTransientState {
        displacement: vec![0.2],
        velocity: vec![0.3],
    };
    let polls = Cell::new(0);
    integrator
        .step(&original, &[0.1], || {
            polls.set(polls.get() + 1);
            false
        })
        .unwrap();
    let count = polls.get();
    polls.set(0);
    assert!(matches!(
        integrator.step(&original, &[0.1], || {
            polls.set(polls.get() + 1);
            polls.get() == count
        }),
        Err(RotorError::Cancelled)
    ));
    assert!(matches!(
        integrator.apply_impulse(&original, &[0.1], || true),
        Err(RotorError::Cancelled)
    ));
    assert!(integrator.step(&original, &[], || false).is_err());
    let too_large = RotorTransientState {
        displacement: vec![0.0],
        velocity: vec![2.0_f64.powi(53)],
    };
    assert!(matches!(
        integrator.apply_impulse(&too_large, &[0.5], || false),
        Err(RotorError::Residual { .. })
    ));
    assert_eq!(original.displacement, vec![0.2]);
    assert_eq!(original.velocity, vec![0.3]);
}
