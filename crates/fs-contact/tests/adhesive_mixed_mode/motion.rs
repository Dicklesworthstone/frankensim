//! G1/G3/G4: coupled full-stick dynamics, work, independent ODE and refusal.
use super::*;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::adhesive::mixed_mode::{motion::*, path::JkrMixedModePathError};
use fs_exec::CancelGate;

const MASS: f64 = 1e-5;
fn scales(source: &SourceJkrMixedMode) -> (f64, f64, f64, f64) {
    let h = source.normal_source().indentation_scale_m();
    let force = PI * 0.05 * 0.01;
    let speed = (force * h / MASS).sqrt();
    (h, force, speed, h / speed)
}
fn budget(source: &SourceJkrMixedMode) -> JkrMixedModeMotionBudget {
    let (h, force, speed, time) = scales(source);
    JkrMixedModeMotionBudget {
        maximum_step_s: 0.5 * time,
        maximum_coordinate_increment_m: 2.0 * h,
        maximum_speed_m_s: 100.0 * speed,
        maximum_force_n: 1e4 * force,
        maximum_iterations: 48,
        maximum_evaluations: 512,
        maximum_backtracks: 24,
        coordinate_absolute_tolerance_m: 1e-14 * h,
        coordinate_relative_tolerance: 1e-13,
        velocity_absolute_tolerance_m_s: 1e-14 * speed,
        velocity_relative_tolerance: 1e-13,
        energy_absolute_tolerance_j: 1e-11 * force * h,
        energy_relative_tolerance: 1e-10,
    }
}
fn mobility() -> [[f64; 3]; 3] {
    [[1.0, 0.2, -0.15], [0.2, 0.8, 0.12], [-0.15, 0.12, 1.1]].map(|row| row.map(|v| v / MASS))
}
fn mul(a: [[f64; 3]; 3], x: [f64; 3]) -> [f64; 3] {
    a.map(|row| (0..3).map(|i| row[i] * x[i]).sum())
}
fn norm(a: [f64; 3]) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}
fn force(source: &SourceJkrMixedMode, q: [f64; 3]) -> [f64; 3] {
    let r = source
        .evaluate_displacement_controlled(q[0], [q[1], q[2]])
        .unwrap();
    [
        r.normal_force_n,
        r.tangential_force_n[0],
        r.tangential_force_n[1],
    ]
}

#[test]
fn g1_original_three_force_equilibrium_is_stationary_with_full_cross_mobility() {
    let (source, shear) = super::path::sourced_law(1e9);
    let (h, _, _, time) = scales(&source);
    let q = [2.0 * h, 0.3 * h, -0.2 * h];
    let applied = mul(mobility(), force(&source, q));
    let mut motion = JkrMixedModeMotion::new(&source, &shear, mobility(), budget(&source)).unwrap();
    let gate = CancelGate::new_clock_free();
    let state = motion.initialize(q, [0.0; 3], &gate).unwrap();
    let step = motion.advance(state, applied, 0.1 * time, &gate).unwrap();
    assert_eq!(step.state().position_m(), q);
    assert_eq!(step.initial_state(), state);
    assert_eq!(step.state().steps(), 1);
    assert_eq!(step.iterations, 0);
    close(
        norm(step.state().velocity_m_s()),
        0.0,
        scales(&source).2,
        2e-15,
    );
    assert_eq!(step.contact_work_j, 0.0);
    assert!(step.evaluations <= budget(&source).maximum_evaluations);
    assert_eq!(step.contact.model_identity(), source.identity());
}

// Independent explicit RK4 uses only original instantaneous forces, not the
// production discrete-gradient, analytic tangent, nonlinear solve or path work.
fn rk4(
    source: &SourceJkrMixedMode,
    w: [[f64; 3]; 3],
    initial: [f64; 6],
    external: [f64; 3],
    duration: f64,
    steps: usize,
) -> [f64; 6] {
    let rhs = |y: [f64; 6]| {
        let reaction = mul(w, force(source, [y[0], y[1], y[2]]));
        [
            y[3],
            y[4],
            y[5],
            external[0] - reaction[0],
            external[1] - reaction[1],
            external[2] - reaction[2],
        ]
    };
    let plus =
        |y: [f64; 6], k: [f64; 6], fraction: f64| core::array::from_fn(|i| y[i] + fraction * k[i]);
    let dt = duration / steps as f64;
    let mut y = initial;
    for _ in 0..steps {
        let k1 = rhs(y);
        let k2 = rhs(plus(y, k1, 0.5 * dt));
        let k3 = rhs(plus(y, k2, 0.5 * dt));
        let k4 = rhs(plus(y, k3, dt));
        y = core::array::from_fn(|i| {
            y[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
        });
    }
    y
}
fn run(
    source: &SourceJkrMixedMode,
    shear: &SourceAdhesiveShear,
    initial: [f64; 6],
    external: [f64; 3],
    duration: f64,
    steps: usize,
) -> (JkrMixedModeMotionState, f64, f64) {
    let mut motion = JkrMixedModeMotion::new(source, shear, mobility(), budget(source)).unwrap();
    let gate = CancelGate::new_clock_free();
    let mut state = motion
        .initialize(
            [initial[0], initial[1], initial[2]],
            [initial[3], initial[4], initial[5]],
            &gate,
        )
        .unwrap();
    let mut work = 0.0;
    let mut largest_residual = 0.0_f64;
    for _ in 0..steps {
        let next = motion
            .advance(state, external, duration / steps as f64, &gate)
            .unwrap();
        assert_eq!(next.initial_state(), state);
        work += next.applied_work_j;
        largest_residual = largest_residual.max(next.energy_residual_j.abs());
        state = next.state();
    }
    (state, work, largest_residual)
}

#[test]
fn g1_coupled_area_and_normal_recoil_converge_to_independent_original_force_ode() {
    let (source, shear) = super::path::sourced_law(1e9);
    let (h, force_scale, speed, time) = scales(&source);
    let initial = [
        2.0 * h,
        0.3 * h,
        -0.2 * h,
        0.1 * speed,
        0.4 * speed,
        -0.3 * speed,
    ];
    let q0 = [initial[0], initial[1], initial[2]];
    let original = source
        .evaluate_displacement_controlled(q0[0], [q0[1], q0[2]])
        .unwrap();
    let external = mul(mobility(), [original.normal_force_n, 0.0, 0.0]);
    let duration = 0.25 * time;
    let oracle = rk4(&source, mobility(), initial, external, duration, 512);
    let coarse = run(&source, &shear, initial, external, duration, 4);
    let fine = run(&source, &shear, initial, external, duration, 8);
    let error = |state: JkrMixedModeMotionState| {
        norm(core::array::from_fn(|i| {
            (state.position_m()[i] - oracle[i]) / h
        })) + norm(core::array::from_fn(|i| {
            (state.velocity_m_s()[i] - oracle[i + 3]) / speed
        }))
    };
    assert!(
        error(coarse.0) > 3.0 * error(fine.0),
        "coarse={}, fine={}",
        error(coarse.0),
        error(fine.0)
    );
    assert!(error(fine.0) < 2e-3);
    assert_eq!(fine.0.steps(), 8);
    close(fine.0.time_s(), duration, duration, 2e-15);
    let end = source
        .evaluate_displacement_controlled(
            fine.0.position_m()[0],
            [fine.0.position_m()[1], fine.0.position_m()[2]],
        )
        .unwrap();
    assert!(
        (end.contact_radius_m - original.contact_radius_m).abs() > 1e-5 * original.contact_radius_m
    );
    assert!((end.normal_force_n - original.normal_force_n).abs() > 1e-5 * force_scale);
    let mass = JkrMixedModeMotion::new(&source, &shear, mobility(), budget(&source))
        .unwrap()
        .effective_mass_matrix_kg();
    let kinetic = |v: [f64; 3]| 0.5 * (0..3).map(|i| v[i] * mul(mass, v)[i]).sum::<f64>();
    let balance = kinetic(fine.0.velocity_m_s()) - kinetic([initial[3], initial[4], initial[5]])
        + end.reversible_energy_j
        - original.reversible_energy_j
        - fine.1;
    close(balance, 0.0, force_scale * h, 2e-9);
    assert!(fine.2 < 1e-9 * force_scale * h);

    let diagonal = core::array::from_fn(|i| {
        core::array::from_fn(|j| if i == j { mobility()[i][j] } else { 0.0 })
    });
    let lost_coupling = rk4(&source, diagonal, initial, external, duration, 512);
    assert!((oracle[0] - lost_coupling[0]).abs() > 1e-4 * h);
}

#[test]
fn g3_elastic_step_reverses_without_hysteresis_and_rotates_with_its_tangent_frame() {
    let (source, shear) = super::path::sourced_law(1e9);
    let (h, _, speed, time) = scales(&source);
    let q = [2.0 * h, 0.3 * h, -0.2 * h];
    let v = [0.1 * speed, 0.4 * speed, -0.3 * speed];
    let external = mul(mobility(), [force(&source, q)[0], 0.0, 0.0]);
    let gate = CancelGate::new_clock_free();
    let mut motion = JkrMixedModeMotion::new(&source, &shear, mobility(), budget(&source)).unwrap();
    let start = motion.initialize(q, v, &gate).unwrap();
    let next = motion.advance(start, external, 0.1 * time, &gate).unwrap();
    let reverse = motion
        .initialize(
            next.state().position_m(),
            next.state().velocity_m_s().map(|x| -x),
            &gate,
        )
        .unwrap();
    let recovered = motion
        .advance(reverse, external, 0.1 * time, &gate)
        .unwrap();
    for i in 0..3 {
        close(recovered.state().position_m()[i], q[i], h, 2e-11);
        close(recovered.state().velocity_m_s()[i], -v[i], speed, 2e-11);
    }
    close(
        next.contact_work_j + recovered.contact_work_j,
        0.0,
        scales(&source).1 * h,
        2e-11,
    );

    let (s, c) = 0.43_f64.sin_cos();
    let rotation = [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]];
    let rotated_w = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            (0..3)
                .flat_map(|k| {
                    (0..3).map(move |l| rotation[i][k] * mobility()[k][l] * rotation[j][l])
                })
                .sum()
        })
    });
    let mut rotated = JkrMixedModeMotion::new(&source, &shear, rotated_w, budget(&source)).unwrap();
    let start = rotated
        .initialize(mul(rotation, q), mul(rotation, v), &gate)
        .unwrap();
    let result = rotated
        .advance(start, mul(rotation, external), 0.1 * time, &gate)
        .unwrap();
    for i in 0..3 {
        close(
            result.state().position_m()[i],
            mul(rotation, next.state().position_m())[i],
            h,
            3e-11,
        );
        close(
            result.state().velocity_m_s()[i],
            mul(rotation, next.state().velocity_m_s())[i],
            speed,
            3e-11,
        );
    }
    close(
        result.contact_work_j,
        next.contact_work_j,
        scales(&source).1 * h,
        3e-11,
    );
}

#[test]
fn g4_source_capacity_work_budget_policy_and_cancellation_refuse_atomically() {
    let (source, shear) = super::path::sourced_law(1e9);
    let (h, _, speed, time) = scales(&source);
    let q = [2.0 * h, 0.3 * h, -0.2 * h];
    let v = [0.1 * speed, 0.4 * speed, -0.3 * speed];
    let gate = CancelGate::new_clock_free();
    let mut motion = JkrMixedModeMotion::new(&source, &shear, mobility(), budget(&source)).unwrap();
    let initial = motion.initialize(q, v, &gate).unwrap();
    let reference = motion
        .advance(initial, [0.0; 3], 0.05 * time, &gate)
        .unwrap();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        motion.advance(initial, [0.0; 3], 0.05 * time, &cancelled),
        Err(JkrMixedModeMotionError::Cancelled)
    ));
    assert_eq!(
        motion
            .advance(initial, [0.0; 3], 0.05 * time, &gate)
            .unwrap(),
        reference
    );
    let mut changed = budget(&source);
    changed.maximum_step_s *= 0.5;
    let mut different = JkrMixedModeMotion::new(&source, &shear, mobility(), changed).unwrap();
    assert!(matches!(
        different.advance(initial, [0.0; 3], 0.05 * time, &gate),
        Err(JkrMixedModeMotionError::StateMismatch)
    ));
    let mut short = budget(&source);
    short.maximum_evaluations = 6;
    let mut limited = JkrMixedModeMotion::new(&source, &shear, mobility(), short).unwrap();
    let seed = limited.initialize(q, v, &gate).unwrap();
    assert!(matches!(
        limited.advance(seed, [0.0; 3], 0.05 * time, &gate),
        Err(JkrMixedModeMotionError::Budget
            | JkrMixedModeMotionError::Path(JkrMixedModePathError::Budget))
    ));
    let (weak_source, weak_shear) = super::path::sourced_law(1.0);
    let weak = JkrMixedModeMotion::new(&weak_source, &weak_shear, mobility(), budget(&weak_source))
        .unwrap();
    assert!(matches!(
        weak.initialize(q, v, &gate),
        Err(JkrMixedModeMotionError::Path(
            JkrMixedModePathError::StaticCapacity { .. }
        ))
    ));
    assert!(JkrMixedModeMotion::new(&source, &weak_shear, mobility(), budget(&source)).is_err());
    let mut indefinite = mobility();
    indefinite[0][0] = -indefinite[0][0];
    assert!(JkrMixedModeMotion::new(&source, &shear, indefinite, budget(&source)).is_err());
}
