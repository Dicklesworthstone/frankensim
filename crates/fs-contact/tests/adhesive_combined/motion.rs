//! Coupled four-coordinate dynamics: native mass metric, source work and hidden-turn refusal.
#[path = "motion/temperature.rs"]
mod temperature;

use super::source;
use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::mixed_mode::{
    SourceJkrMixedMode,
    combined::{motion::*, path::JkrCombinedPathError},
    torsion::motion::{JkrTorsionalMotion, JkrTorsionalMotionBudget},
};
use fs_exec::CancelGate;

const MASS: f64 = 1e-5;
const LOWER: [[f64; 4]; 4] = [
    [1.0, 0.0, 0.0, 0.0],
    [0.2, 1.0, 0.0, 0.0],
    [-0.15, 0.25, 1.0, 0.0],
    [0.1, -0.2, 0.3, 1.0],
];

fn scales(source: &SourceJkrMixedMode) -> ([f64; 4], f64, f64) {
    let normal = source.normal_source();
    let length = normal.indentation_scale_m();
    let force = PI * normal.normal_state().adhesion_energy_j_per_m2() * normal.effective_radius_m();
    let response = source
        .evaluate_combined_displacement_controlled([0.0; 4])
        .unwrap();
    let stiffness = response.coupled_stiffness.unwrap();
    let energy = force * length;
    let shear = (energy / stiffness[1][1]).sqrt();
    let angle = (energy / stiffness[3][3]).sqrt();
    (
        [length, shear, shear, angle],
        force,
        (MASS * length / force).sqrt(),
    )
}

fn metric_scale(source: &SourceJkrMixedMode) -> [f64; 4] {
    let (q, _, _) = scales(source);
    q.map(|value| (value / q[0]) / MASS.sqrt())
}

fn mobility(source: &SourceJkrMixedMode) -> [[f64; 4]; 4] {
    let scale = metric_scale(source);
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            scale[i] * scale[j] * (0..4).map(|k| LOWER[i][k] * LOWER[j][k]).sum::<f64>()
        })
    })
}

fn budget(source: &SourceJkrMixedMode) -> JkrCombinedMotionBudget {
    let (q, force, time) = scales(source);
    JkrCombinedMotionBudget {
        maximum_step_s: 0.5 * time,
        maximum_coordinate_increment: q.map(|x| 2.0 * x),
        maximum_rate: q.map(|x| 100.0 * x / time),
        maximum_effort: q.map(|x| 1e4 * force * q[0] / x),
        maximum_iterations: 48,
        maximum_evaluations: 512,
        maximum_backtracks: 24,
        coordinate_absolute_tolerance: q.map(|x| 1e-14 * x),
        coordinate_relative_tolerance: 1e-13,
        rate_absolute_tolerance: q.map(|x| 1e-14 * x / time),
        rate_relative_tolerance: 1e-13,
        energy_absolute_tolerance_j: 1e-11 * force * q[0],
        energy_relative_tolerance: 1e-10,
    }
}

fn start(source: &SourceJkrMixedMode) -> ([f64; 4], [f64; 4]) {
    let (scale, _, time) = scales(source);
    (
        [
            0.6 * scale[0],
            0.12 * scale[1],
            -0.08 * scale[2],
            0.1 * scale[3],
        ],
        core::array::from_fn(|i| [0.1, -0.05, 0.07, 0.09][i] * scale[i] / time),
    )
}
fn mul(w: [[f64; 4]; 4], value: [f64; 4]) -> [f64; 4] {
    w.map(|row| row.into_iter().zip(value).map(|(a, b)| a * b).sum())
}
fn force(source: &SourceJkrMixedMode, q: [f64; 4]) -> [f64; 4] {
    source
        .evaluate_combined_displacement_controlled(q)
        .unwrap()
        .force
}
fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * scale,
        "actual={a:e}, expected={b:e}, scale={scale:e}, tolerance={tolerance:e}"
    );
}

// The authored physical factor defines this independent kinetic oracle. It
// does not reuse the production Cholesky, inverse, normalization or linear solve.
fn kinetic(source: &SourceJkrMixedMode, v: [f64; 4]) -> f64 {
    let scale = metric_scale(source);
    let a = v[0] / scale[0];
    let b = v[1] / scale[1] - 0.2 * a;
    let c = v[2] / scale[2] + 0.15 * a - 0.25 * b;
    let d = v[3] / scale[3] - 0.1 * a + 0.2 * b - 0.3 * c;
    0.5 * (a * a + b * b + c * c + d * d)
}

#[test]
fn full_four_coordinate_mobility_preserves_original_work_and_actual_recoil() {
    let source = source();
    let gate = CancelGate::new();
    let (scale, f, time) = scales(&source);
    let w = mobility(&source);
    assert!(w[3][3] / w[0][0] > 1e6);
    for i in 0..4 {
        for j in 0..4 {
            assert_ne!(w[i][j], 0.0);
        }
    }
    let motion = JkrCombinedMotion::new(&source, w, budget(&source)).unwrap();
    let (q, v) = start(&source);
    let initial = motion.initial_state(q, v, &gate).unwrap();
    let applied = core::array::from_fn(|i| [0.1, -0.08, 0.03, 0.04][i] * f * scale[0] / scale[i]);
    let external = mul(w, applied);
    let dt = 0.06 * time;
    let step = motion.advance(initial, dt, external, &gate).unwrap();
    assert_eq!(motion.advance(initial, dt, external, &gate).unwrap(), step);
    assert_eq!(step.initial_state(), initial);
    assert_eq!(step.initial_contact, initial.contact());
    assert_eq!(step.contact, step.state().contact());
    assert_eq!(step.state().steps(), 1);
    assert_eq!(step.state().time_s(), dt);
    assert_eq!(step.contact.model_identity(), source.identity());
    assert!(step.iterations > 0);
    assert!(step.evaluations <= budget(&source).maximum_evaluations);
    let reaction = mul(w, step.mean_effort);
    let dq: [f64; 4] = core::array::from_fn(|i| step.state().coordinates()[i] - q[i]);
    for i in 0..4 {
        close(
            dq[i],
            0.5 * dt * (v[i] + step.state().rates()[i]),
            scale[i],
            2e-12,
        );
        close(
            step.state().rates()[i] - v[i],
            dt * (external[i] - reaction[i]),
            scale[i] / time,
            2e-12,
        );
        close(
            step.effort_impulse[i],
            step.mean_effort[i] * dt,
            f * scale[0] * dt / scale[i],
            2e-14,
        );
        close(
            step.effort_time_moment[i],
            0.5 * dt * step.effort_impulse[i],
            f * scale[0] * dt * dt / scale[i],
            2e-14,
        );
    }
    let energy = f * scale[0];
    close(
        step.contact_work_j,
        step.contact.reversible_energy_j - initial.contact().reversible_energy_j,
        energy,
        3e-12,
    );
    close(
        step.kinetic_energy_change_j,
        kinetic(&source, step.state().rates()) - kinetic(&source, v),
        energy,
        3e-12,
    );
    let expected_work: f64 = applied.into_iter().zip(dq).map(|(a, b)| a * b).sum();
    close(step.applied_work_j, expected_work, energy, 2e-12);
    close(
        step.kinetic_energy_change_j + step.contact_work_j,
        expected_work,
        energy,
        2e-11,
    );
    assert!(step.source_work_residual_j.abs() < 3e-12 * energy);
    for k in 0..=32 {
        let t = dt * k as f64 / 32.0;
        let point =
            core::array::from_fn(|i| q[i] + t * v[i] + 0.5 * t * t * (external[i] - reaction[i]));
        for (i, effort) in force(&source, point).into_iter().enumerate() {
            assert!(effort.abs() <= step.maximum_source_effort[i]);
        }
    }
}

// Independent explicit RK4 queries only original point forces. It shares no
// discrete gradient, analytic path tangent, quadrature or production integrator.
fn rk4(
    source: &SourceJkrMixedMode,
    w: [[f64; 4]; 4],
    mut state: [f64; 8],
    external: [f64; 4],
    duration: f64,
) -> [f64; 8] {
    let rhs = |y: [f64; 8]| {
        let reaction = mul(w, force(source, core::array::from_fn(|i| y[i])));
        core::array::from_fn(|i| {
            if i < 4 {
                y[i + 4]
            } else {
                external[i - 4] - reaction[i - 4]
            }
        })
    };
    let plus =
        |y: [f64; 8], slope: [f64; 8], dt: f64| core::array::from_fn(|i| y[i] + dt * slope[i]);
    let dt = duration / 512.0;
    for _ in 0..512 {
        let a = rhs(state);
        let b = rhs(plus(state, a, 0.5 * dt));
        let c = rhs(plus(state, b, 0.5 * dt));
        let d = rhs(plus(state, c, dt));
        state =
            core::array::from_fn(|i| state[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    state
}

#[test]
fn coupled_inertial_motion_converges_to_pointwise_ode_and_retraces_reversibly() {
    let source = source();
    let gate = CancelGate::new();
    let (scale, _, time) = scales(&source);
    let w = mobility(&source);
    let (q, v) = start(&source);
    let external = mul(w, [force(&source, q)[0], 0.0, 0.0, 0.0]);
    let duration = 0.3 * time;
    let initial_ode = core::array::from_fn(|i| if i < 4 { q[i] } else { v[i - 4] });
    let reference = rk4(&source, w, initial_ode, external, duration);
    let motion = JkrCombinedMotion::new(&source, w, budget(&source)).unwrap();
    let initial = motion.initial_state(q, v, &gate).unwrap();
    let run = |count: usize| {
        let mut state = initial;
        for _ in 0..count {
            state = motion
                .advance(state, duration / count as f64, external, &gate)
                .unwrap()
                .state();
        }
        state
    };
    let coarse = run(4);
    let fine = run(8);
    let error = |state: JkrCombinedMotionState| -> f64 {
        (0..4)
            .map(|i| {
                (state.coordinates()[i] - reference[i]).abs() / scale[i]
                    + time * (state.rates()[i] - reference[i + 4]).abs() / scale[i]
            })
            .sum()
    };
    assert!(
        error(coarse) > 3.0 * error(fine),
        "coarse={}, fine={}",
        error(coarse),
        error(fine)
    );
    assert!(error(fine) < 2e-3);
    assert_eq!(fine.steps(), 8);
    close(fine.time_s(), duration, duration, 4e-15);
    assert!(
        (fine.contact().contact_radius_m - initial.contact().contact_radius_m).abs()
            > 1e-6 * initial.contact().contact_radius_m
    );
    let mut returned = motion
        .initial_state(fine.coordinates(), fine.rates().map(|x| -x), &gate)
        .unwrap();
    for _ in 0..8 {
        returned = motion
            .advance(returned, duration / 8.0, external, &gate)
            .unwrap()
            .state();
    }
    for i in 0..4 {
        close(returned.coordinates()[i], q[i], scale[i], 3e-10);
        close(returned.rates()[i], -v[i], scale[i] / time, 3e-10);
    }
}

fn rotate(value: [f64; 4]) -> [f64; 4] {
    let (s, c) = 0.63_f64.sin_cos();
    [
        value[0],
        c * value[1] - s * value[2],
        s * value[1] + c * value[2],
        value[3],
    ]
}

#[test]
fn tangent_gauge_covariance_and_original_two_coordinate_limit_are_preserved() {
    let source = source();
    let gate = CancelGate::new();
    let (scale, f, time) = scales(&source);
    let w = mobility(&source);
    let (q, v) = start(&source);
    let dt = 0.04 * time;
    let external = mul(w, [0.1 * f, -0.05 * f, 0.03 * f, 0.0]);
    let original = JkrCombinedMotion::new(&source, w, budget(&source)).unwrap();
    let before = original.initial_state(q, v, &gate).unwrap();
    let result = original.advance(before, dt, external, &gate).unwrap();
    // Build Q W Q^T by transforming its columns and then its rows.
    let columns: [[f64; 4]; 4] =
        core::array::from_fn(|j| rotate(core::array::from_fn(|i| w[i][j])));
    let rotated_w = core::array::from_fn(|i| rotate(core::array::from_fn(|j| columns[j][i])));
    let rotated = JkrCombinedMotion::new(&source, rotated_w, budget(&source)).unwrap();
    let before_rotated = rotated.initial_state(rotate(q), rotate(v), &gate).unwrap();
    let result_rotated = rotated
        .advance(before_rotated, dt, rotate(external), &gate)
        .unwrap();
    for i in 0..4 {
        close(
            result_rotated.state().coordinates()[i],
            rotate(result.state().coordinates())[i],
            scale[i],
            3e-11,
        );
        close(
            result_rotated.state().rates()[i],
            rotate(result.state().rates())[i],
            scale[i] / time,
            3e-11,
        );
        close(
            result_rotated.mean_effort[i],
            rotate(result.mean_effort)[i],
            f * scale[0] / scale[i],
            3e-11,
        );
    }
    close(
        result_rotated.contact.contact_radius_m,
        result.contact.contact_radius_m,
        result.contact.contact_radius_m,
        3e-12,
    );
    close(
        result_rotated.contact_work_j,
        result.contact_work_j,
        f * scale[0],
        3e-12,
    );

    let mut block = w;
    for i in 0..4 {
        for j in 0..4 {
            if i != j && (i == 1 || i == 2 || j == 1 || j == 2) {
                block[i][j] = 0.0;
            }
        }
    }
    let b = budget(&source);
    let indices = [0, 3];
    let pair_budget = JkrTorsionalMotionBudget {
        maximum_step_s: b.maximum_step_s,
        maximum_coordinate_increment: indices.map(|i| b.maximum_coordinate_increment[i]),
        maximum_rate: indices.map(|i| b.maximum_rate[i]),
        maximum_effort: indices.map(|i| b.maximum_effort[i]),
        maximum_iterations: b.maximum_iterations,
        maximum_evaluations: b.maximum_evaluations,
        maximum_backtracks: b.maximum_backtracks,
        coordinate_absolute_tolerance: indices.map(|i| b.coordinate_absolute_tolerance[i]),
        coordinate_relative_tolerance: b.coordinate_relative_tolerance,
        rate_absolute_tolerance: indices.map(|i| b.rate_absolute_tolerance[i]),
        rate_relative_tolerance: b.rate_relative_tolerance,
        energy_absolute_tolerance_j: b.energy_absolute_tolerance_j,
        energy_relative_tolerance: b.energy_relative_tolerance,
    };
    let two_w = indices.map(|i| indices.map(|j| block[i][j]));
    let pair = JkrTorsionalMotion::new(&source, two_w, pair_budget).unwrap();
    let all = JkrCombinedMotion::new(&source, block, b).unwrap();
    let initial_pair = pair
        .initial_state([q[0], q[3]], [v[0], v[3]], &gate)
        .unwrap();
    let initial_all = all
        .initial_state([q[0], 0.0, 0.0, q[3]], [v[0], 0.0, 0.0, v[3]], &gate)
        .unwrap();
    let pair_step = pair.advance(initial_pair, dt, [0.0; 2], &gate).unwrap();
    let all_step = all.advance(initial_all, dt, [0.0; 4], &gate).unwrap();
    for (i, j) in indices.into_iter().enumerate() {
        close(
            all_step.state().coordinates()[j],
            pair_step.state().coordinates()[i],
            scale[j],
            3e-11,
        );
        close(
            all_step.state().rates()[j],
            pair_step.state().rates()[i],
            scale[j] / time,
            3e-11,
        );
    }
    assert_eq!(all_step.state().coordinates()[1..3], [0.0; 2]);
    assert_eq!(all_step.state().rates()[1..3], [0.0; 2]);
    close(
        all_step.contact_work_j,
        pair_step.contact_work_j,
        f * scale[0],
        3e-12,
    );
}

#[test]
fn hidden_joint_shear_twist_excursion_refuses_although_endpoints_and_each_mode_fit() {
    let source = source();
    let gate = CancelGate::new();
    let (scale, _, time) = scales(&source);
    let delta = -scale[0];
    let shear = 0.8 * source.critical_tangential_displacement_m(delta).unwrap();
    let twist = 0.8
        * source
            .torsional_displacement_limit(delta)
            .unwrap()
            .unwrap()
            .twist_angle_rad;
    let peak = [delta, 0.6 * shear, 0.8 * shear, twist];
    source
        .evaluate_combined_displacement_controlled([delta, peak[1], peak[2], 0.0])
        .unwrap();
    source
        .evaluate_combined_displacement_controlled([delta, 0.0, 0.0, twist])
        .unwrap();
    assert!(
        source
            .evaluate_combined_displacement_controlled(peak)
            .is_err()
    );
    let q = [delta, 0.0, 0.0, 0.0];
    source.combined_discrete_path(q, q, &gate).unwrap();
    let dt = 0.2 * time;
    let v = [
        0.0,
        4.0 * peak[1] / dt,
        4.0 * peak[2] / dt,
        4.0 * peak[3] / dt,
    ];
    let w = mobility(&source);
    let mut b = budget(&source);
    b.maximum_coordinate_increment = [2.0 * scale[0], 2.0 * shear, 2.0 * shear, 2.0 * twist];
    b.maximum_rate = core::array::from_fn(|i| b.maximum_rate[i].max(2.0 * v[i].abs()));
    b.maximum_effort = b.maximum_effort.map(|x| x * 100.0);
    let motion = JkrCombinedMotion::new(&source, w, b).unwrap();
    let initial = motion.initial_state(q, v, &gate).unwrap();
    let reaction = mul(w, force(&source, q));
    // q1=q0 and v1=-v0 satisfy the discrete equation. The actual quadratic
    // trajectory reaches the inadmissible combined peak at half the duration.
    let external = core::array::from_fn(|i| reaction[i] - 2.0 * v[i] / dt);
    assert!(matches!(
        motion.advance(initial, dt, external, &gate),
        Err(JkrCombinedMotionError::Path(
            JkrCombinedPathError::UnresolvedEnvelope
        ))
    ));
    assert_eq!(initial.coordinates(), q);
    assert_eq!(initial.rates(), v);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.steps(), 0);
}

#[test]
fn unresolved_inertia_budget_cancellation_and_wrong_owner_preserve_accepted_state() {
    let source = source();
    let gate = CancelGate::new();
    let (scale, _, time) = scales(&source);
    let w = mobility(&source);
    let mut nonsymmetric = w;
    nonsymmetric[0][3] += 0.2 * (w[0][0] * w[3][3]).sqrt();
    let singular = core::array::from_fn(|i| core::array::from_fn(|j| (w[i][i] * w[j][j]).sqrt()));
    let mut indefinite = w;
    indefinite[0][1] = 2.0 * (w[0][0] * w[1][1]).sqrt();
    indefinite[1][0] = indefinite[0][1];
    let mut zero = w;
    zero[2][2] = 0.0;
    for matrix in [nonsymmetric, singular, indefinite, zero] {
        assert!(matches!(
            JkrCombinedMotion::new(&source, matrix, budget(&source)),
            Err(JkrCombinedMotionError::Invalid(_))
        ));
    }
    let mut extreme = [[0.0; 4]; 4];
    for i in 0..4 {
        extreme[i][i] = [1e-120, 1e-40, 1e40, 1e120][i];
    }
    assert_eq!(
        JkrCombinedMotion::new(&source, extreme, budget(&source))
            .unwrap()
            .mobility(),
        extreme
    );

    let fold_q = [
        0.99 * source.normal_source().minimum_attached_indentation_m(),
        0.0,
        0.0,
        0.0,
    ];
    let stiffness = source
        .evaluate_combined_displacement_controlled(fold_q)
        .unwrap()
        .coupled_stiffness
        .unwrap()[0][0];
    assert!(stiffness < 0.0);
    let mut diagonal = [[0.0; 4]; 4];
    for i in 0..4 {
        diagonal[i][i] = w[i][i];
    }
    let dt = 4.0 / (-diagonal[0][0] * stiffness).sqrt();
    let mut coarse_budget = budget(&source);
    coarse_budget.maximum_step_s = dt;
    let coarse = JkrCombinedMotion::new(&source, diagonal, coarse_budget).unwrap();
    let folded = coarse.initial_state(fold_q, [0.0; 4], &gate).unwrap();
    assert_eq!(
        coarse.advance(folded, dt, mul(diagonal, force(&source, fold_q)), &gate),
        Err(JkrCombinedMotionError::InertialTangent)
    );

    let (q, v) = start(&source);
    let external = mul(w, force(&source, q));
    let motion = JkrCombinedMotion::new(&source, w, budget(&source)).unwrap();
    let state = motion.initial_state(q, v, &gate).unwrap();
    let successful = motion.advance(state, 0.05 * time, external, &gate).unwrap();
    let mut one_query = budget(&source);
    one_query.maximum_evaluations = 1;
    let limited = JkrCombinedMotion::new(&source, w, one_query).unwrap();
    let limited_state = limited.initial_state(q, v, &gate).unwrap();
    assert_eq!(
        limited.advance(limited_state, 0.05 * time, external, &gate),
        Err(JkrCombinedMotionError::Budget)
    );
    assert_eq!(
        limited.advance(state, 0.05 * time, external, &gate),
        Err(JkrCombinedMotionError::StateMismatch)
    );
    let other = SourceJkrMixedMode::new(
        source.normal_source().clone(),
        fs_contact::normal_patch::adhesive::mixed_mode::JkrMixedModeEnvelope {
            maximum_tangential_displacement_to_radius: 0.09,
            ..source.envelope()
        },
    )
    .unwrap();
    let other_motion = JkrCombinedMotion::new(&other, w, budget(&source)).unwrap();
    assert_eq!(
        other_motion.advance(state, 0.05 * time, external, &gate),
        Err(JkrCombinedMotionError::StateMismatch)
    );
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_eq!(
        motion.advance(state, 0.05 * time, external, &cancelled),
        Err(JkrCombinedMotionError::Cancelled)
    );
    assert_eq!(
        motion.advance(state, 0.05 * time, external, &gate).unwrap(),
        successful
    );
    assert_eq!(state.coordinates(), q);
    assert_eq!(state.rates(), v);
    assert_eq!(state.steps(), 0);
    assert!(
        successful
            .maximum_source_effort
            .iter()
            .all(|f| f.is_finite())
    );
    close(successful.state().time_s(), 0.05 * time, time, 2e-15);
}
