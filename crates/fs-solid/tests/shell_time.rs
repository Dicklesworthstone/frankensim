//! Shell transient integration, passive dissipation, and checkpoint battery (bead `frankensim-b8bxd.9.2`).

use fs_solid::shell_time::{
    DynamicState, PreparedNewmark, SHELL_TIME_INTEGRATOR_SCHEMA_V1, ShellTimeCheckpoint,
    ShellTimeConfig, ShellTimeError, step_newmark,
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

    let f_zero = vec![0.0; n];
    let mut prepared = PreparedNewmark::new(n, &mass, None, &stiffness, &config).unwrap();
    let mut state = prepared
        .initialize(0.0, &[1.0, 0.5], &[0.0, 0.0], &f_zero)
        .unwrap();
    let initial_energy = state.total_energy_j;
    assert_eq!(state.acceleration, [-2.5, 0.0]);

    let num_steps = 100;
    for _ in 0..num_steps {
        state = step_newmark(n, &mass, None, &stiffness, &state, &f_zero, &config)
            .expect("undamped step succeeds");
    }

    // G0: average acceleration conserves the quadratic energy from the very
    // first step when the initial acceleration satisfies the actual equation.
    let energy_rel_diff = (state.total_energy_j - initial_energy).abs() / initial_energy;
    assert!(
        energy_rel_diff < 1.0e-13,
        "undamped energy relative drift {energy_rel_diff}"
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

    let f_zero = vec![0.0; n];
    let mut prepared = PreparedNewmark::new(n, &mass, Some(&damping), &stiffness, &config).unwrap();
    let mut state = prepared
        .initialize(0.0, &[1.0, 0.5], &[0.5, -0.5], &f_zero)
        .unwrap();
    let initial_energy = state.total_energy_j;
    let num_steps = 100;
    let mut prev_energy = initial_energy;
    let mut heat = 0.0;

    for _ in 0..num_steps {
        state = prepared
            .step(&state, &f_zero)
            .expect("damped step succeeds");
        let work = prepared.last_step_work().unwrap();
        assert!(work.damping_work_j >= 0.0);
        assert!(work.algorithmic_energy_defect_j.abs() < 2.0e-15);
        heat += work.damping_work_j;
        assert!(
            state.total_energy_j <= prev_energy + 1.0e-12,
            "energy must be non-increasing under passive damping: {} > {}",
            state.total_energy_j,
            prev_energy
        );
        prev_energy = state.total_energy_j;
    }

    assert!((state.total_energy_j + heat - initial_energy).abs() < 1.0e-13);

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

fn scalar_state(u: f64, v: f64, a: f64) -> DynamicState {
    DynamicState {
        time_s: 0.0,
        step: 0,
        displacement: vec![u],
        velocity: vec![v],
        acceleration: vec![a],
        kinetic_energy_j: 0.0,
        strain_energy_j: 0.0,
        total_energy_j: 0.0,
    }
}

fn action_2d(matrix: &[f64; 4], x: &[f64]) -> [f64; 2] {
    [
        matrix[0] * x[0] + matrix[1] * x[1],
        matrix[2] * x[0] + matrix[3] * x[1],
    ]
}

fn dot_2d(x: &[f64], y: &[f64]) -> f64 {
    x[0] * y[0] + x[1] * y[1]
}

#[test]
fn shell_time_static_preload_preserves_damped_increment_forces() {
    // Linear superposition gives the same rate update about a force-balanced
    // static displacement. The old absolute predictor rounded away K*du, so
    // changing only that preload changed the damping/vibration dynamics.
    let mass = [2.0, 0.25, 0.25, 1.0];
    let stiffness = [6.0, -2.0, -2.0, 4.0];
    let config = ShellTimeConfig {
        dt_s: 0.1,
        ..Default::default()
    };
    for damping_scale in [0.0, 1.0, 100.0] {
        let damping = [
            0.4 * damping_scale,
            0.1 * damping_scale,
            0.1 * damping_scale,
            0.3 * damping_scale,
        ];
        let velocity = [0.75, -0.25];
        let cv = action_2d(&damping, &velocity);
        let det = mass[0] * mass[3] - mass[1] * mass[2];
        let acceleration = [
            (-mass[3] * cv[0] + mass[1] * cv[1]) / det,
            (-mass[0] * cv[1] + mass[2] * cv[0]) / det,
        ];
        let state = DynamicState {
            time_s: 0.0,
            step: 0,
            displacement: vec![0.0; 2],
            velocity: velocity.to_vec(),
            acceleration: acceleration.to_vec(),
            kinetic_energy_j: 0.0,
            strain_energy_j: 0.0,
            total_energy_j: 0.0,
        };
        let mut prepared =
            PreparedNewmark::new(2, &mass, Some(&damping), &stiffness, &config).unwrap();
        let reference = prepared.step(&state, &[0.0; 2]).unwrap();
        for offset in [1.0e4, 1.0e16] {
            let mut preloaded = state.clone();
            preloaded.displacement = vec![offset, -2.0 * offset];
            let force = action_2d(&stiffness, &preloaded.displacement);
            let shifted = prepared.step(&preloaded, &force).unwrap();
            for (actual, expected) in shifted
                .acceleration
                .iter()
                .chain(&shifted.velocity)
                .zip(reference.acceleration.iter().chain(&reference.velocity))
            {
                let tolerance = 8.0 * f64::EPSILON * expected.abs().max(1.0);
                assert!(
                    (actual - expected).abs() <= tolerance,
                    "preload {offset:e}, damping scale {damping_scale:e}: rate {actual:e} differs from {expected:e}"
                );
            }
            assert_eq!(
                shifted,
                step_newmark(
                    2,
                    &mass,
                    Some(&damping),
                    &stiffness,
                    &preloaded,
                    &force,
                    &config,
                )
                .unwrap()
            );
        }
    }
}

#[test]
fn shell_time_full_matrices_balance_forcing_work_and_damping_loss() {
    let mass = [2.0, 0.25, 0.25, 1.0];
    let stiffness = [6.0, -2.0, -2.0, 4.0];
    // Positive definite and not a linear combination of M and K.
    let damping = [0.4, 0.1, 0.1, 0.3];
    let config = ShellTimeConfig {
        dt_s: 0.01,
        ..Default::default()
    };
    let mut prepared = PreparedNewmark::new(2, &mass, Some(&damping), &stiffness, &config).unwrap();
    let force = |t: f64| [0.75 * (1.7 * t).cos(), -0.4 * (0.8 * t).sin()];
    let u = vec![0.3, -0.2];
    let v = vec![-0.1, 0.4];
    let ku = action_2d(&stiffness, &u);
    let cv = action_2d(&damping, &v);
    let f0 = force(0.0);
    let rhs = [f0[0] - ku[0] - cv[0], f0[1] - ku[1] - cv[1]];
    let det = mass[0] * mass[3] - mass[1] * mass[2];
    let a = vec![
        (mass[3] * rhs[0] - mass[1] * rhs[1]) / det,
        (mass[0] * rhs[1] - mass[2] * rhs[0]) / det,
    ];
    let kinetic = 0.5 * dot_2d(&v, &action_2d(&mass, &v));
    let strain = 0.5 * dot_2d(&u, &ku);
    let mut state = prepared.initialize(0.0, &u, &v, &f0).unwrap();
    for (actual, expected) in state.acceleration.iter().zip(a) {
        assert!((actual - expected).abs() < 3.0e-16);
    }
    assert_eq!(state.kinetic_energy_j, kinetic);
    assert_eq!(state.strain_energy_j, strain);
    assert!(prepared.last_step_work().is_none());
    for _ in 0..500 {
        let f_prev = force(state.time_s);
        let f_next = force(state.time_s + config.dt_s);
        let result = prepared.step_with_work(&state, &f_prev, &f_next).unwrap();
        let next = result.state;
        assert_eq!(prepared.last_step_work(), Some(&result.work));
        assert_eq!(
            next,
            step_newmark(
                2,
                &mass,
                Some(&damping),
                &stiffness,
                &state,
                &f_next,
                &config
            )
            .unwrap()
        );
        let mean_v = [
            0.5 * (state.velocity[0] + next.velocity[0]),
            0.5 * (state.velocity[1] + next.velocity[1]),
        ];
        let mean_f = [0.5 * (f_prev[0] + f_next[0]), 0.5 * (f_prev[1] + f_next[1])];
        let work = config.dt_s * dot_2d(&mean_f, &mean_v);
        let loss = config.dt_s * dot_2d(&mean_v, &action_2d(&damping, &mean_v));
        let balance = next.total_energy_j - state.total_energy_j - work + loss;
        assert!(balance.abs() < 2.0e-14, "step {}: {balance:e}", next.step);
        assert!((result.work.applied_work_j - work).abs() < 5.0e-17);
        assert!((result.work.damping_work_j - loss).abs() < 5.0e-17);
        assert!(result.work.damping_work_j >= 0.0);
        assert!(result.work.algorithmic_energy_defect_j.abs() < 2.0e-14);
        assert!(
            (result.work.mechanical_energy_change_j - (next.total_energy_j - state.total_energy_j))
                .abs()
                < 2.0e-14
        );
        let ma = action_2d(&mass, &next.acceleration);
        let cv = action_2d(&damping, &next.velocity);
        let ku = action_2d(&stiffness, &next.displacement);
        for i in 0..2 {
            assert!((ma[i] + cv[i] + ku[i] - f_next[i]).abs() < 3.0e-14);
        }
        if next.step == 250 {
            let checkpoint = ShellTimeCheckpoint::new(&next);
            checkpoint.verify().unwrap();
            let restored = DynamicState {
                time_s: checkpoint.time_s,
                step: checkpoint.step,
                displacement: checkpoint.u,
                velocity: checkpoint.v,
                acceleration: checkpoint.a,
                kinetic_energy_j: 0.0,
                strain_energy_j: 0.0,
                total_energy_j: 0.0,
            };
            let mut rebuilt =
                PreparedNewmark::new(2, &mass, Some(&damping), &stiffness, &config).unwrap();
            assert!(rebuilt.last_step_work().is_none());
            let following_force = force(next.time_s + config.dt_s);
            assert_eq!(
                prepared.step(&next, &following_force).unwrap(),
                rebuilt.step(&restored, &following_force).unwrap()
            );
            assert_eq!(prepared.last_step_work(), rebuilt.last_step_work());
        }
        state = next;
    }
}

#[test]
fn shell_time_forced_response_converges_at_second_order() {
    // Manufactured exact motion u=sin(2t), v=2cos(2t), a=-4sin(2t).
    // M=2, C=0.7, K=13 imply f=5sin(2t)+1.4cos(2t).
    let error = |steps: usize| {
        let config = ShellTimeConfig {
            dt_s: 1.0 / steps as f64,
            ..Default::default()
        };
        let mut prepared = PreparedNewmark::new(1, &[2.0], Some(&[0.7]), &[13.0], &config).unwrap();
        let mut state = scalar_state(0.0, 2.0, 0.0);
        for _ in 0..steps {
            let phase = 2.0 * (state.time_s + config.dt_s);
            state = prepared
                .step(&state, &[5.0 * phase.sin() + 1.4 * phase.cos()])
                .unwrap();
        }
        (state.displacement[0] - 2.0_f64.sin()).abs()
            + (state.velocity[0] - 2.0 * 2.0_f64.cos()).abs()
    };
    let coarse = error(50);
    let fine = error(100);
    assert!(
        (3.8..4.2).contains(&(coarse / fine)),
        "errors {coarse:e}, {fine:e}"
    );
    assert!(fine < 5.0e-4, "fine-step error {fine:e}");
}

#[test]
fn shell_time_viscous_heat_and_applied_work_converge_to_physical_integrals() {
    // G1: exact forced motion u=sin(2t). The discrete midpoint-power sums
    // converge to the independent continuous integrals; a zero discrete
    // energy defect alone would not establish this physical heat accuracy.
    let exact_heat = 1.4 + 0.35 * 4.0_f64.sin();
    let exact_work = exact_heat + 1.25 * (1.0 - 4.0_f64.cos());
    let errors = |steps: usize| {
        let config = ShellTimeConfig {
            dt_s: 1.0 / steps as f64,
            ..Default::default()
        };
        let mut prepared = PreparedNewmark::new(1, &[2.0], Some(&[0.7]), &[13.0], &config).unwrap();
        let mut state = prepared.initialize(0.0, &[0.0], &[2.0], &[1.4]).unwrap();
        let initial_energy = state.total_energy_j;
        let mut heat = 0.0;
        let mut work = 0.0;
        for _ in 0..steps {
            let phase = 2.0 * (state.time_s + config.dt_s);
            state = prepared
                .step(&state, &[5.0 * phase.sin() + 1.4 * phase.cos()])
                .unwrap();
            let interval = prepared.last_step_work().unwrap();
            assert!(interval.damping_work_j >= 0.0);
            heat += interval.damping_work_j;
            work += interval.applied_work_j;
        }
        assert!((state.total_energy_j - initial_energy + heat - work).abs() < 2.0e-13);
        [(heat - exact_heat).abs(), (work - exact_work).abs()]
    };
    for (coarse, fine) in errors(40).into_iter().zip(errors(80)) {
        assert!(
            (3.8..4.2).contains(&(coarse / fine)),
            "errors {coarse:e}, {fine:e}"
        );
        assert!(fine < 1.0e-3);
    }
}

#[test]
fn shell_time_numerical_energy_loss_is_not_viscous_heat() {
    // G0: a dissipative Newmark parameter choice removes energy even with no
    // physical damping. That signed defect must remain separate from heat.
    let config = ShellTimeConfig {
        dt_s: 0.1,
        beta: 0.3025,
        gamma: 0.6,
        ..Default::default()
    };
    let mut prepared = PreparedNewmark::new(1, &[1.0], None, &[4.0], &config).unwrap();
    let state = prepared.initialize(0.0, &[1.0], &[0.0], &[0.0]).unwrap();
    let result = prepared.step_with_work(&state, &[0.0], &[0.0]).unwrap();
    assert_eq!(result.work.applied_work_j, 0.0);
    assert_eq!(result.work.damping_work_j, 0.0);
    assert!(result.work.algorithmic_energy_defect_j < -1.0e-5);
    assert_eq!(
        result.work.algorithmic_energy_defect_j,
        result.work.mechanical_energy_change_j
    );

    // Operator admission deliberately does not certify passivity. A supplied
    // active C retains its signed work instead of silently clipping it to heat.
    let config = ShellTimeConfig {
        dt_s: 0.1,
        ..Default::default()
    };
    let mut active = PreparedNewmark::new(1, &[1.0], Some(&[-0.5]), &[0.0], &config).unwrap();
    let state = active.initialize(0.0, &[0.0], &[1.0], &[0.0]).unwrap();
    let result = active.step_with_work(&state, &[0.0], &[0.0]).unwrap();
    assert!(result.work.damping_work_j < 0.0);
    assert!(result.work.mechanical_energy_change_j > 0.0);
    assert!(result.work.algorithmic_energy_defect_j.abs() < 1.0e-15);
}

#[test]
fn shell_time_actual_old_force_exposes_inconsistent_start_and_preserves_accepted_work() {
    let config = ShellTimeConfig {
        dt_s: 0.1,
        ..Default::default()
    };
    let mut prepared = PreparedNewmark::new(1, &[1.0], None, &[4.0], &config).unwrap();
    let mut state = prepared.initialize(0.0, &[1.0], &[0.0], &[0.0]).unwrap();
    // Legacy callers may supply guessed acceleration and stale energy caches.
    // Explicit physical endpoint forces reveal that defect; work is never
    // manufactured by subtracting those caches or balancing the endpoint loss.
    state.acceleration[0] = 0.0;
    state.kinetic_energy_j = f64::NAN;
    state.strain_energy_j = f64::NAN;
    state.total_energy_j = f64::NAN;
    let result = prepared.step_with_work(&state, &[0.0], &[0.0]).unwrap();
    assert_eq!(result.work.applied_work_j, 0.0);
    assert_eq!(result.work.damping_work_j, 0.0);
    assert!(result.work.algorithmic_energy_defect_j < -0.01);
    assert_eq!(prepared.last_step_work(), Some(&result.work));

    assert!(
        prepared
            .step_with_work(&result.state, &[f64::NAN], &[0.0])
            .is_err()
    );
    assert_eq!(prepared.last_step_work(), Some(&result.work));
    assert!(prepared.initialize(0.0, &[], &[0.0], &[0.0]).is_err());
    assert_eq!(prepared.last_step_work(), Some(&result.work));
    prepared.initialize(0.0, &[1.0], &[0.0], &[0.0]).unwrap();
    assert!(prepared.last_step_work().is_none());
}

#[test]
fn shell_time_work_retains_small_quadratic_energy_change_on_large_preload() {
    // G3: endpoint energy subtraction loses part of this representable
    // displacement's energy. The interval calculation uses vector increments
    // and still exposes the energy defect caused by rounded absolute motion.
    let config = ShellTimeConfig {
        dt_s: 1.0e-8,
        ..Default::default()
    };
    let mut prepared = PreparedNewmark::new(1, &[1.0], None, &[1.0], &config).unwrap();
    let state = prepared
        .initialize(0.0, &[1.0e8], &[1.0], &[1.0e8])
        .unwrap();
    let result = prepared.step_with_work(&state, &[1.0e8], &[1.0e8]).unwrap();
    let du = result.state.displacement[0] - state.displacement[0];
    let expected_change = du * state.displacement[0] + 0.5 * du * du;
    assert!((result.work.mechanical_energy_change_j - expected_change).abs() < 1.0e-14);
    assert!(
        (result.work.mechanical_energy_change_j
            - (result.state.total_energy_j - state.total_energy_j))
            .abs()
            > 0.1
    );
    assert!((result.work.applied_work_j - 1.0).abs() < 1.0e-14);
}
