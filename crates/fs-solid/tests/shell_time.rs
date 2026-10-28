//! Shell transient integration, passive dissipation, and checkpoint battery (bead `frankensim-b8bxd.9.2`).

use fs_solid::shell_time::{
    DynamicState, SHELL_TIME_INTEGRATOR_SCHEMA_V1, ShellTimeCheckpoint, ShellTimeConfig,
    ShellTimeError, step_newmark,
};

fn sample_system() -> (usize, Vec<f64>, Vec<f64>, Vec<f64>) {
    // 2-DOF spring-mass system:
    // M = [[2, 0], [0, 1]]
    // K = [[6, -2], [-2, 4]]
    // C = 0.05 * M + 0.02 * K (Rayleigh damping)
    let n = 2;
    let mass = vec![2.0, 0.0, 0.0, 1.0];
    let stiffness = vec![6.0, -2.0, -2.0, 4.0];
    let mut damping = vec![0.0; 4];
    for i in 0..4 {
        damping[i] = 0.05 * mass[i] + 0.02 * stiffness[i];
    }
    (n, mass, damping, stiffness)
}

#[test]
fn shell_time_undamped_free_oscillation_conserves_energy() {
    let (n, mass, _, stiffness) = sample_system();
    let config = ShellTimeConfig {
        dt_s: 1.0e-3,
        duration_s: 0.1,
        max_steps: 1_000,
        beta: 0.25,
        gamma: 0.50,
    };

    let mut state = DynamicState {
        time_s: 0.0,
        step: 0,
        displacement: vec![1.0, 0.5],
        velocity: vec![0.0, 0.0],
        acceleration: vec![0.0, 0.0],
        kinetic_energy_j: 0.0,
        strain_energy_j: 3.5, // 0.5 * (1*6*1 + 2*1*(-2)*0.5 + 0.5*4*0.5) = 0.5 * (6 - 2 + 1) = 2.5 + ...
        total_energy_j: 3.5,
    };

    let f_zero = vec![0.0; n];
    let initial_energy = 0.5
        * (state.displacement[0]
            * (stiffness[0] * state.displacement[0] + stiffness[1] * state.displacement[1])
            + state.displacement[1]
                * (stiffness[2] * state.displacement[0] + stiffness[3] * state.displacement[1]));
    state.strain_energy_j = initial_energy;
    state.total_energy_j = initial_energy;

    let num_steps = 100;
    for _ in 0..num_steps {
        state = step_newmark(n, &mass, None, &stiffness, &state, &f_zero, &config)
            .expect("undamped step succeeds");
    }

    // Newmark average acceleration is unconditionally stable and symplectic-like; energy oscillates slightly around initial
    let energy_rel_diff = (state.total_energy_j - initial_energy).abs() / initial_energy;
    assert!(
        energy_rel_diff < 1.0e-2,
        "undamped energy relative drift {energy_rel_diff} must be < 1%"
    );
}

#[test]
fn shell_time_passive_damping_dissipates_energy_monotonically() {
    let (n, mass, damping, stiffness) = sample_system();
    let config = ShellTimeConfig {
        dt_s: 1.0e-3,
        duration_s: 0.1,
        max_steps: 1_000,
        beta: 0.25,
        gamma: 0.50,
    };

    let mut state = DynamicState {
        time_s: 0.0,
        step: 0,
        displacement: vec![1.0, 0.5],
        velocity: vec![0.5, -0.5],
        acceleration: vec![0.0, 0.0],
        kinetic_energy_j: 0.0,
        strain_energy_j: 0.0,
        total_energy_j: 0.0,
    };
    let initial_energy = 0.5
        * (mass[0] * state.velocity[0] * state.velocity[0]
            + mass[3] * state.velocity[1] * state.velocity[1])
        + 0.5
            * (state.displacement[0]
                * (stiffness[0] * state.displacement[0] + stiffness[1] * state.displacement[1])
                + state.displacement[1]
                    * (stiffness[2] * state.displacement[0]
                        + stiffness[3] * state.displacement[1]));
    state.total_energy_j = initial_energy;

    let f_zero = vec![0.0; n];
    let num_steps = 100;
    let mut prev_energy = f64::INFINITY;

    for _ in 0..num_steps {
        state = step_newmark(
            n,
            &mass,
            Some(&damping),
            &stiffness,
            &state,
            &f_zero,
            &config,
        )
        .expect("damped step succeeds");
        if prev_energy.is_finite() {
            assert!(
                state.total_energy_j <= prev_energy + 1.0e-12,
                "energy must be non-increasing under passive damping: {} > {}",
                state.total_energy_j,
                prev_energy
            );
        }
        prev_energy = state.total_energy_j;
    }

    assert!(
        state.total_energy_j < initial_energy,
        "energy must decay over time: {} < {}",
        state.total_energy_j,
        initial_energy
    );
}

#[test]
fn shell_time_checkpoint_and_resume_are_exact() {
    let (n, mass, damping, stiffness) = sample_system();
    let config = ShellTimeConfig {
        dt_s: 1.0e-3,
        duration_s: 0.05,
        max_steps: 1_000,
        beta: 0.25,
        gamma: 0.50,
    };

    let initial = DynamicState {
        time_s: 0.0,
        step: 0,
        displacement: vec![1.0, 0.5],
        velocity: vec![0.5, -0.5],
        acceleration: vec![0.0, 0.0],
        kinetic_energy_j: 0.0,
        strain_energy_j: 0.0,
        total_energy_j: 0.0,
    };

    let f_zero = vec![0.0; n];

    // Run 50 steps continuously
    let mut state_continuous = initial.clone();
    for _ in 0..50 {
        state_continuous = step_newmark(
            n,
            &mass,
            Some(&damping),
            &stiffness,
            &state_continuous,
            &f_zero,
            &config,
        )
        .expect("continuous step succeeds");
    }

    // Run 25 steps, take checkpoint, resume and run 25 steps
    let mut state_resumed = initial;
    for _ in 0..25 {
        state_resumed = step_newmark(
            n,
            &mass,
            Some(&damping),
            &stiffness,
            &state_resumed,
            &f_zero,
            &config,
        )
        .expect("first half step succeeds");
    }

    let checkpoint = ShellTimeCheckpoint::new(&state_resumed);
    assert_eq!(checkpoint.schema_version, SHELL_TIME_INTEGRATOR_SCHEMA_V1);
    checkpoint.verify().expect("checkpoint verification passes");

    for _ in 0..25 {
        state_resumed = step_newmark(
            n,
            &mass,
            Some(&damping),
            &stiffness,
            &state_resumed,
            &f_zero,
            &config,
        )
        .expect("second half step succeeds");
    }

    assert_eq!(state_continuous.step, state_resumed.step);
    assert_eq!(state_continuous.time_s, state_resumed.time_s);
    assert_eq!(state_continuous.displacement, state_resumed.displacement);
    assert_eq!(state_continuous.velocity, state_resumed.velocity);
    assert_eq!(state_continuous.acceleration, state_resumed.acceleration);
}

#[test]
fn shell_time_refuses_corrupt_checkpoint() {
    let state = DynamicState {
        time_s: 0.01,
        step: 10,
        displacement: vec![0.1, 0.2],
        velocity: vec![0.01, -0.02],
        acceleration: vec![0.0, 0.0],
        kinetic_energy_j: 0.01,
        strain_energy_j: 0.02,
        total_energy_j: 0.03,
    };

    let mut checkpoint = ShellTimeCheckpoint::new(&state);
    checkpoint.verify().expect("initial checkpoint is valid");

    // Tamper with displacement
    checkpoint.u[0] = 0.999;
    let err = checkpoint
        .verify()
        .expect_err("tampered checkpoint must fail verification");
    assert!(matches!(err, ShellTimeError::CorruptCheckpoint { .. }));
}

#[test]
fn shell_time_refuses_invalid_dt() {
    let (n, mass, _, stiffness) = sample_system();
    let config = ShellTimeConfig {
        dt_s: 0.0, // Invalid
        ..Default::default()
    };
    let state = DynamicState {
        time_s: 0.0,
        step: 0,
        displacement: vec![1.0, 0.5],
        velocity: vec![0.0, 0.0],
        acceleration: vec![0.0, 0.0],
        kinetic_energy_j: 0.0,
        strain_energy_j: 0.0,
        total_energy_j: 0.0,
    };

    let err = step_newmark(n, &mass, None, &stiffness, &state, &[0.0, 0.0], &config)
        .expect_err("must refuse dt <= 0");
    assert!(matches!(err, ShellTimeError::InvalidTimeStep { .. }));
}

/// G0/G3: the coordinate datum cannot create a force or change free momentum.
#[test]
fn free_drift_keeps_momentum_when_positions_are_large_or_steps_are_short() {
    for (position, dt) in [(1e9, 1e-6), (1e6, 1e-5), (1.0, 1e-12)] {
        let config = ShellTimeConfig {
            dt_s: dt,
            ..Default::default()
        };
        let mut state = DynamicState {
            time_s: 0.0,
            step: 0,
            displacement: vec![position],
            velocity: vec![1.0],
            acceleration: vec![0.0],
            kinetic_energy_j: 1.0,
            strain_energy_j: 0.0,
            total_energy_j: 1.0,
        };
        for _ in 0..64 {
            state = step_newmark(1, &[2.0], None, &[0.0], &state, &[0.0], &config).unwrap();
            assert_eq!(state.velocity, [1.0]);
            assert_eq!(state.acceleration, [0.0]);
            assert_eq!(state.kinetic_energy_j, 1.0);
        }
        // A ramp from zero force to six newtons contributes h*(0+6)/2.
        let next = step_newmark(1, &[2.0], None, &[0.0], &state, &[6.0], &config).unwrap();
        assert_eq!(next.acceleration, [3.0]);
        assert!(
            (2.0 * (next.velocity[0] - state.velocity[0]) - 3.0 * dt).abs() <= 4.0 * f64::EPSILON
        );
    }
}

/// G0: internal spring/damper forces cancel; the full consistent mass maps
/// velocities to physical momentum, and external work closes with dissipation.
#[test]
fn coupled_consistent_mass_obeys_impulse_and_force_work_balances() {
    let mass = [2.0, 0.5, 0.5, 1.25];
    let stiffness = [4.0, -4.0, -4.0, 4.0];
    let damping = [0.3, -0.3, -0.3, 0.3];
    let config = ShellTimeConfig {
        dt_s: 0.002,
        ..Default::default()
    };
    let energy = |q: &[f64], v: &[f64]| {
        0.5 * (v[0] * (2.0 * v[0] + 0.5 * v[1]) + v[1] * (0.5 * v[0] + 1.25 * v[1]))
            + 2.0 * (q[0] - q[1]).powi(2)
    };
    let mut force = [0.2, -0.4];
    let q = [0.1, -0.2];
    let v = [0.3, -0.5];
    let internal = 4.0 * (q[0] - q[1]) + 0.3 * (v[0] - v[1]);
    let rhs = [force[0] - internal, force[1] + internal];
    // Independent inverse of M; initial acceleration satisfies equilibrium.
    let a = [
        (1.25 * rhs[0] - 0.5 * rhs[1]) / 2.25,
        (2.0 * rhs[1] - 0.5 * rhs[0]) / 2.25,
    ];
    let mut state = DynamicState {
        time_s: 0.0,
        step: 0,
        displacement: q.to_vec(),
        velocity: v.to_vec(),
        acceleration: a.to_vec(),
        kinetic_energy_j: 0.0,
        strain_energy_j: 0.0,
        total_energy_j: energy(&q, &v),
    };
    let mut expected_momentum = 2.5 * v[0] + 1.75 * v[1];
    for step in 1..=500 {
        let next_force = [
            0.2 + 0.002 * f64::from(step),
            -0.4 + 0.001 * f64::from(step),
        ];
        let next = step_newmark(
            2,
            &mass,
            Some(&damping),
            &stiffness,
            &state,
            &next_force,
            &config,
        )
        .unwrap();
        expected_momentum +=
            0.5 * config.dt_s * (force[0] + force[1] + next_force[0] + next_force[1]);
        let momentum = 2.5 * next.velocity[0] + 1.75 * next.velocity[1];
        assert!((momentum - expected_momentum).abs() < 2e-12);
        let mean_v = [
            0.5 * (state.velocity[0] + next.velocity[0]),
            0.5 * (state.velocity[1] + next.velocity[1]),
        ];
        let work = config.dt_s
            * 0.5
            * ((force[0] + next_force[0]) * mean_v[0] + (force[1] + next_force[1]) * mean_v[1]);
        let loss = config.dt_s * 0.3 * (mean_v[0] - mean_v[1]).powi(2);
        let change = next.total_energy_j - energy(&state.displacement, &state.velocity);
        assert!((change - work + loss).abs() < 1e-12);
        state = next;
        force = next_force;
    }
}

/// G3: scaling all masses, stiffnesses, damping and forces leaves motion intact.
/// A dimensional absolute pivot threshold would incorrectly reject small bodies.
#[test]
fn inertial_dynamics_are_invariant_to_uniform_physical_scale() {
    let (n, mass, damping, stiffness) = sample_system();
    let config = ShellTimeConfig {
        dt_s: 0.1,
        ..Default::default()
    };
    let initial = DynamicState {
        time_s: 0.0,
        step: 0,
        displacement: vec![0.25, -0.125],
        velocity: vec![0.125, 0.25],
        acceleration: vec![-0.5, 0.75],
        kinetic_energy_j: 0.0,
        strain_energy_j: 0.0,
        total_energy_j: 0.0,
    };
    let force = [0.5, -0.25];
    let reference = step_newmark(
        n,
        &mass,
        Some(&damping),
        &stiffness,
        &initial,
        &force,
        &config,
    )
    .unwrap();
    for scale in [1e-24, 1e24] {
        let scaled = |values: &[f64]| values.iter().map(|x| scale * x).collect::<Vec<_>>();
        let actual = step_newmark(
            n,
            &scaled(&mass),
            Some(&scaled(&damping)),
            &scaled(&stiffness),
            &initial,
            &scaled(&force),
            &config,
        )
        .unwrap();
        for (a, b) in actual
            .displacement
            .iter()
            .chain(&actual.velocity)
            .chain(&actual.acceleration)
            .zip(
                reference
                    .displacement
                    .iter()
                    .chain(&reference.velocity)
                    .chain(&reference.acceleration),
            )
        {
            assert!((a - b).abs() < 1e-13);
        }
        assert!((actual.total_energy_j / scale - reference.total_energy_j).abs() < 1e-13);
    }
}

#[test]
fn newmark_refuses_invalid_dynamics_without_panicking_or_publishing() {
    let state = DynamicState {
        time_s: 0.0,
        step: 0,
        displacement: vec![0.0],
        velocity: vec![1.0],
        acceleration: vec![0.0],
        kinetic_energy_j: 0.5,
        strain_energy_j: 0.0,
        total_energy_j: 0.5,
    };
    let config = ShellTimeConfig::default();
    assert!(matches!(
        step_newmark(1, &[], None, &[0.0], &state, &[0.0], &config),
        Err(ShellTimeError::DimensionMismatch { .. })
    ));
    assert!(matches!(
        step_newmark(1, &[1.0], None, &[0.0], &state, &[], &config),
        Err(ShellTimeError::DimensionMismatch { .. })
    ));
    for (beta, gamma) in [(0.0, 0.5), (-0.25, 0.5), (0.25, f64::NAN)] {
        let invalid = ShellTimeConfig {
            beta,
            gamma,
            ..config.clone()
        };
        assert!(matches!(
            step_newmark(1, &[1.0], None, &[0.0], &state, &[0.0], &invalid),
            Err(ShellTimeError::InvalidParameter { .. })
        ));
    }
    assert!(step_newmark(1, &[1.0], None, &[0.0], &state, &[f64::INFINITY], &config).is_err());
    assert!(step_newmark(1, &[0.0], None, &[0.0], &state, &[0.0], &config).is_err());
    assert_eq!(
        step_newmark(1, &[1.0], None, &[0.0], &state, &[0.0], &config)
            .unwrap()
            .velocity,
        state.velocity
    );
}
