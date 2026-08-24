//! Shell transient integration, passive dissipation, and checkpoint battery (bead `frankensim-b8bxd.9.2`).

use fs_solid::shell_time::{
    step_newmark, DynamicState, ShellTimeCheckpoint, ShellTimeConfig, ShellTimeError,
    SHELL_TIME_INTEGRATOR_SCHEMA_V1,
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
    let initial_energy = 0.5 * (state.displacement[0] * (stiffness[0] * state.displacement[0] + stiffness[1] * state.displacement[1])
        + state.displacement[1] * (stiffness[2] * state.displacement[0] + stiffness[3] * state.displacement[1]));
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
                    * (stiffness[2] * state.displacement[0] + stiffness[3] * state.displacement[1]));
    state.total_energy_j = initial_energy;

    let f_zero = vec![0.0; n];
    let num_steps = 100;
    let mut prev_energy = f64::INFINITY;

    for _ in 0..num_steps {
        state = step_newmark(n, &mass, Some(&damping), &stiffness, &state, &f_zero, &config)
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
        state_continuous =
            step_newmark(n, &mass, Some(&damping), &stiffness, &state_continuous, &f_zero, &config)
                .expect("continuous step succeeds");
    }

    // Run 25 steps, take checkpoint, resume and run 25 steps
    let mut state_resumed = initial;
    for _ in 0..25 {
        state_resumed =
            step_newmark(n, &mass, Some(&damping), &stiffness, &state_resumed, &f_zero, &config)
                .expect("first half step succeeds");
    }

    let checkpoint = ShellTimeCheckpoint::new(&state_resumed);
    assert_eq!(checkpoint.schema_version, SHELL_TIME_INTEGRATOR_SCHEMA_V1);
    checkpoint.verify().expect("checkpoint verification passes");

    for _ in 0..25 {
        state_resumed =
            step_newmark(n, &mass, Some(&damping), &stiffness, &state_resumed, &f_zero, &config)
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
    let err = checkpoint.verify().expect_err("tampered checkpoint must fail verification");
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
