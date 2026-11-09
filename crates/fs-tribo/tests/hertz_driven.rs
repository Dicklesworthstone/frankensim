//! Powered full-mobility Hertz shear/spin: original static-set separation,
//! analytic mixed events, independent continuum refinement and mode honesty.
use core::f64::consts::PI;
use fs_tribo::sliding_patch::contact::driven::{
    HertzDrivenControls, HertzDrivenError, HertzDrivenInput, HertzDrivenMode, HertzDrivenReport,
};
use fs_tribo::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchError, HertzWrenchParameters,
};
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};

const R: f64 = 0.1;
type V = [f64; 3];
type M = [[f64; 3]; 3];

fn law(static_mu: f64, kinetic_mu: f64) -> HertzWrenchContact {
    HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: 2.0,
        contact_radius_m: R,
        static_mu,
        kinetic_mu,
    })
    .unwrap()
}
fn config() -> HertzWrenchConfig {
    HertzWrenchConfig {
        maximum_static_iterations: 96,
        maximum_iterations: 256,
        maximum_evaluations: 1024,
        force_tolerance_n: 2e-11,
        moment_tolerance_nm: 2e-12,
        displacement_tolerance_m: 1e-12,
        energy_tolerance_j: 2e-9,
        relative_matrix_tolerance: 1e-13,
        maximum_force_n: 20.0,
        maximum_spin_moment_nm: 2.0,
    }
}
fn controls() -> HertzDrivenControls {
    HertzDrivenControls {
        maximum_step_s: 0.01,
        minimum_step_s: 1e-14,
        maximum_attempts: 100_000,
        maximum_evaluations: 400_001,
        scaled_rate_tolerance_m_s: 2e-9,
        rest_scaled_rate_tolerance_m_s: 2e-10,
    }
}
fn full() -> M {
    [[2.0, 0.3, 0.4], [0.3, 1.7, -0.2], [0.4, -0.2, 3.5]]
}
fn make_input(
    q: V,
    matrix: M,
    acceleration: V,
    dt: f64,
    mode: HertzDrivenMode,
) -> HertzDrivenInput {
    HertzDrivenInput {
        velocity_m_s: [q[0], q[1]],
        spin_rate_rad_s: q[2] / R,
        scaled_mobility_per_kg: matrix,
        driving_acceleration_m_s2: [acceleration[0], acceleration[1]],
        driving_spin_acceleration_rad_s2: acceleration[2] / R,
        dt_s: dt,
        mode,
    }
}
fn run(law: &HertzWrenchContact, input: &HertzDrivenInput) -> HertzDrivenReport {
    law.advance_driven(input, config(), controls(), || false)
        .unwrap()
}
fn rate(report: &HertzDrivenReport) -> V {
    [
        report.velocity_after_m_s[0],
        report.velocity_after_m_s[1],
        R * report.spin_rate_after_rad_s,
    ]
}
fn impulse(report: &HertzDrivenReport) -> V {
    [
        report.force_impulse_ns[0],
        report.force_impulse_ns[1],
        report.spin_angular_impulse_nms / R,
    ]
}
fn travel(report: &HertzDrivenReport) -> V {
    [
        report.slip_increment_m[0],
        report.slip_increment_m[1],
        R * report.spin_increment_rad,
    ]
}
fn add(a: V, b: V) -> V {
    core::array::from_fn(|i| a[i] + b[i])
}
fn scale(a: V, s: f64) -> V {
    a.map(|v| v * s)
}
fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn multiply(a: M, v: V) -> V {
    a.map(|row| dot(row, v))
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance {tolerance:.3e}"
    );
}
fn close_vector(a: V, b: V, tolerance: f64) {
    for i in 0..3 {
        close(a[i], b[i], tolerance);
    }
}

#[test]
fn static_hold_and_joint_breakaway_use_the_original_hertz_support_body() {
    let law = law(0.8, 0.5);
    let matrix = full();
    let capacity = 1.6;
    // This demand lies in the inscribed ball of radius
    // C*(3*pi/16)/sqrt(1+(3*pi/16)^2), so membership is independent of
    // any implementation's projection or static classifier.
    let demand = scale([0.2, 0.1, 0.08], capacity);
    let held = run(
        &law,
        &make_input(
            [0.0; 3],
            matrix,
            multiply(matrix, demand),
            0.4,
            HertzDrivenMode::Resting,
        ),
    );
    assert_eq!(held.mode, HertzDrivenMode::Resting);
    assert_eq!(rate(&held), [0.0; 3]);
    assert_eq!(travel(&held), [0.0; 3]);
    close_vector(impulse(&held), scale(demand, -0.4), 2e-13);
    close_vector(held.driving_scaled_wrench_n, demand, 2e-13);
    close(held.dissipated_energy_j, 0.0, 0.0);
    close(held.drive_work_j, 0.0, 0.0);
    assert_eq!(held.first_arrest_time_s, None);

    // Exact pressure moments at scaled slip/spin (1,0,1) give support
    // (9*pi/32,0,9*pi/128)*C and plane Fx+M/a<=45*pi*C/128.
    // Each separate force/spin cap admits this drive, but the joint plane
    // rejects it. Independent force/torque saturation would falsely hold.
    let demand = scale([0.8, 0.0, 0.35], capacity);
    assert!(demand[0] < capacity);
    assert!(demand[2] < (3.0 * PI / 16.0) * capacity);
    assert!(demand[0] + demand[2] > (45.0 * PI / 128.0) * capacity);
    let released_input = make_input(
        [0.0; 3],
        matrix,
        multiply(matrix, demand),
        0.25,
        HertzDrivenMode::Resting,
    );
    let released = run(&law, &released_input);
    assert_eq!(released.mode, HertzDrivenMode::Sliding);
    assert!(rate(&released).iter().any(|v| v.abs() > 1e-3));
    assert!(released.dissipated_energy_j > 0.0);
    assert!(released.drive_work_j > released.dissipated_energy_j);
    assert_eq!(released.first_arrest_time_s, None);
    assert_eq!(released.accepted_steps, 1);
    // Homogeneity makes release from exact rest a constant-acceleration ray,
    // even with the original off-diagonal mobility. Phase splitting must
    // preserve that exact solution and retained mode.
    let first = run(
        &law,
        &HertzDrivenInput {
            dt_s: 0.1,
            ..released_input
        },
    );
    let second = run(
        &law,
        &make_input(
            rate(&first),
            matrix,
            multiply(matrix, demand),
            0.15,
            first.mode,
        ),
    );
    close_vector(rate(&second), rate(&released), 2e-9);
    close_vector(
        add(travel(&first), travel(&second)),
        travel(&released),
        2e-10,
    );
    close(
        first.dissipated_energy_j + second.dissipated_energy_j,
        released.dissipated_energy_j,
        2e-10,
    );
}

#[test]
fn mixed_shear_spin_ray_uses_exact_pressure_moments_and_then_holds() {
    let law = law(0.8, 0.5);
    let matrix = [[1.0, 0.0, 0.0], [0.0, 1.7, 0.0], [0.0, 0.0, 4.0]];
    let initial = 0.8;
    let drive = 0.2;
    let stop = initial / (9.0 * PI / 32.0 - drive);
    let dwell = 0.4;
    let report = run(
        &law,
        &make_input(
            [initial, 0.0, initial],
            matrix,
            [drive, 0.0, drive],
            stop + dwell,
            HertzDrivenMode::Sliding,
        ),
    );
    let wrench = [-9.0 * PI / 32.0, 0.0, -9.0 * PI / 128.0];
    let demand = [drive, 0.0, drive / 4.0];
    let distance = [0.5 * initial * stop, 0.0, 0.5 * initial * stop];
    assert_eq!(report.mode, HertzDrivenMode::Resting);
    assert_eq!(rate(&report), [0.0; 3]);
    assert!(report.arrest_time_is_exact);
    close(report.first_arrest_time_s.unwrap(), stop, 2e-14);
    close_vector(travel(&report), distance, 2e-14);
    close_vector(
        impulse(&report),
        add(scale(wrench, stop), scale(demand, -dwell)),
        2e-13,
    );
    close(report.dissipated_energy_j, -dot(wrench, distance), 2e-14);
    close(report.drive_work_j, dot(demand, distance), 2e-14);
    close(
        report.kinetic_energy_change_j,
        -0.625 * initial * initial,
        2e-14,
    );
    close(report.energy_residual_j, 0.0, 2e-14);
    assert_eq!(report.accepted_steps, 2);
}

#[test]
fn opposing_drive_stops_and_restarts_without_negative_heat_or_static_capture() {
    // N=2, mu_k=.5: |friction|=1. Demand -3 exceeds static capacity2.
    // Wxx=2 gives acceleration -8 before zero and -4 after reversal.
    let law = law(1.0, 0.5);
    let matrix = [[2.0, 0.0, 0.0], [0.0, 1.2, 0.0], [0.0, 0.0, 3.0]];
    let report = run(
        &law,
        &make_input(
            [0.6, 0.0, 0.0],
            matrix,
            [-6.0, 0.0, 0.0],
            0.2,
            HertzDrivenMode::Sliding,
        ),
    );
    assert_eq!(report.mode, HertzDrivenMode::Sliding);
    assert!(report.arrest_time_is_exact);
    close(report.first_arrest_time_s.unwrap(), 0.075, 1e-15);
    close_vector(rate(&report), [-0.5, 0.0, 0.0], 2e-12);
    close_vector(travel(&report), [-0.00875, 0.0, 0.0], 2e-13);
    close_vector(impulse(&report), [0.05, 0.0, 0.0], 2e-12);
    close(report.dissipated_energy_j, 0.05375, 2e-12);
    close(report.drive_work_j, 0.02625, 2e-12);
    close(report.kinetic_energy_change_j, -0.0275, 2e-12);
    close(report.energy_residual_j, 0.0, 2e-12);
}

#[test]
fn retained_static_mode_and_tiny_sliding_mode_keep_distinct_hysteresis() {
    let law = law(0.5, 0.25);
    let matrix = [[1.7, 0.0, 0.0], [0.0, 1.2, 0.0], [0.0, 0.0, 2.0]];
    let drive = [1.7 * 0.75, 0.0, 0.0];
    let rest = make_input([0.0; 3], matrix, drive, 0.3, HertzDrivenMode::Resting);
    let held = run(&law, &rest);
    assert_eq!(held.mode, HertzDrivenMode::Resting);
    close(held.dissipated_energy_j, 0.0, 0.0);
    let tiny = 1e-14;
    assert!(tiny < controls().rest_scaled_rate_tolerance_m_s);
    let sliding = make_input(
        [tiny, 0.0, 0.0],
        matrix,
        drive,
        0.3,
        HertzDrivenMode::Sliding,
    );
    let moving = run(&law, &sliding);
    assert_eq!(moving.mode, HertzDrivenMode::Sliding);
    assert_eq!(moving.first_arrest_time_s, None);
    close(rate(&moving)[0], tiny + 1.7 * (0.75 - 0.5) * 0.3, 2e-14);
    let zero_sliding = HertzDrivenInput {
        mode: HertzDrivenMode::Sliding,
        ..rest
    };
    assert!(matches!(
        law.advance_driven(&zero_sliding, config(), controls(), || false),
        Err(HertzDrivenError::Invalid(_))
    ));
    let false_rest = HertzDrivenInput {
        mode: HertzDrivenMode::Resting,
        velocity_m_s: [1e-4, 0.0],
        ..rest
    };
    assert!(matches!(
        law.advance_driven(&false_rest, config(), controls(), || false),
        Err(HertzDrivenError::Invalid(_))
    ));
    assert_eq!(moving, run(&law, &sliding));
}

fn reference(input: &HertzDrivenInput, steps: usize) -> [f64; 10] {
    // Independent fixed-grid RK4 of the original instantaneous pressure law.
    // No driven projection, event logic, midpoint heat or production controller
    // enters this oracle. Constitutive spatial closed forms are shared exactly.
    let patch = HertzSlidingPatch::new(HertzSlidingPatchParameters {
        normal_force_n: 2.0,
        contact_radius_m: R,
        friction_coefficient: 0.5,
    })
    .unwrap();
    let drive = [
        input.driving_acceleration_m_s2[0],
        input.driving_acceleration_m_s2[1],
        R * input.driving_spin_acceleration_rad_s2,
    ];
    let rhs = |state: [f64; 10]| {
        let w = patch.evaluate([state[0], state[1]], state[2] / R).unwrap();
        let wrench = [w.force_n[0], w.force_n[1], w.spin_moment_nm / R];
        let acceleration = add(drive, multiply(input.scaled_mobility_per_kg, wrench));
        [
            acceleration[0],
            acceleration[1],
            acceleration[2],
            wrench[0],
            wrench[1],
            wrench[2],
            state[0],
            state[1],
            state[2],
            w.dissipated_power_w,
        ]
    };
    let mut state = [0.0; 10];
    state[..3].copy_from_slice(&[
        input.velocity_m_s[0],
        input.velocity_m_s[1],
        R * input.spin_rate_rad_s,
    ]);
    let h = input.dt_s / steps as f64;
    for _ in 0..steps {
        let a = rhs(state);
        let b = rhs(core::array::from_fn(|i| state[i] + 0.5 * h * a[i]));
        let c = rhs(core::array::from_fn(|i| state[i] + 0.5 * h * b[i]));
        let d = rhs(core::array::from_fn(|i| state[i] + h * c[i]));
        state =
            core::array::from_fn(|i| state[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    state
}

#[test]
fn full_mobility_curved_drive_refines_to_independent_pressure_ode() {
    let law = law(0.8, 0.5);
    let demand = [0.1, -0.05, 0.08];
    let input = make_input(
        [0.8, -0.3, 0.5],
        full(),
        multiply(full(), demand),
        0.18,
        HertzDrivenMode::Sliding,
    );
    let coarse_reference = reference(&input, 2048);
    let expected = reference(&input, 4096);
    for i in 0..10 {
        close(coarse_reference[i], expected[i], 2e-12);
    }
    let mut errors = Vec::new();
    for steps in [8, 16, 32, 64] {
        let controls = HertzDrivenControls {
            maximum_step_s: input.dt_s / steps as f64,
            ..controls()
        };
        let report = law
            .advance_driven(&input, config(), controls, || false)
            .unwrap();
        assert_eq!(report.mode, HertzDrivenMode::Sliding);
        assert_eq!(report.first_arrest_time_s, None);
        assert!(report.dissipated_energy_j > 0.0);
        close(report.drive_work_j, dot(demand, travel(&report)), 2e-13);
        close(report.energy_residual_j, 0.0, 2e-9);
        let actual = [
            rate(&report)[0],
            rate(&report)[1],
            rate(&report)[2],
            impulse(&report)[0],
            impulse(&report)[1],
            impulse(&report)[2],
            travel(&report)[0],
            travel(&report)[1],
            travel(&report)[2],
            report.dissipated_energy_j,
        ];
        errors.push(
            (0..10)
                .map(|i| (actual[i] - expected[i]).abs())
                .fold(0.0_f64, f64::max),
        );
    }
    for pair in errors.windows(2) {
        assert!(
            pair[1] < 0.56 * pair[0],
            "first-order error ladder {errors:?}"
        );
    }
    assert!(errors[3] < 4e-4, "{errors:?}");
    eprintln!("driven original full-mobility continuum errors: {errors:?}");
}

#[test]
fn curved_subcritical_arrest_is_timestep_resolved_and_static_dwell_has_no_heat() {
    let law = law(0.8, 0.5);
    let drive = multiply(full(), [0.03, -0.02, 0.01]);
    let report = run(
        &law,
        &make_input(
            [0.16, -0.05, 0.1],
            full(),
            drive,
            0.5,
            HertzDrivenMode::Sliding,
        ),
    );
    assert_eq!(report.mode, HertzDrivenMode::Resting);
    assert!(!report.arrest_time_is_exact);
    assert!(
        report
            .first_arrest_time_s
            .is_some_and(|t| t > 0.0 && t < 0.5)
    );
    assert_eq!(rate(&report), [0.0; 3]);
    let dwell = run(&law, &make_input([0.0; 3], full(), drive, 1.5, report.mode));
    assert_eq!(dwell.mode, HertzDrivenMode::Resting);
    close(dwell.dissipated_energy_j, 0.0, 0.0);
    close(dwell.drive_work_j, 0.0, 0.0);
    assert_eq!(travel(&dwell), [0.0; 3]);
}

#[test]
fn budgets_effort_and_cancellation_refuse_without_mutating_contact_or_input() {
    let law = law(0.8, 0.5);
    let input = make_input(
        [0.8, -0.3, 0.5],
        full(),
        [0.1, 0.02, -0.01],
        0.08,
        HertzDrivenMode::Sliding,
    );
    let before = input;
    assert_eq!(
        law.advance_driven(&input, config(), controls(), || true),
        Err(HertzDrivenError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        law.advance_driven(&input, config(), controls(), || {
            polls += 1;
            polls > 12
        }),
        Err(HertzDrivenError::Cancelled)
    );
    assert_eq!(
        law.advance_driven(
            &input,
            config(),
            HertzDrivenControls {
                maximum_attempts: 1,
                ..controls()
            },
            || false
        ),
        Err(HertzDrivenError::AttemptBudget)
    );
    assert_eq!(
        law.advance_driven(
            &input,
            config(),
            HertzDrivenControls {
                maximum_evaluations: 1,
                ..controls()
            },
            || false
        ),
        Err(HertzDrivenError::EvaluationBudget)
    );
    let rest = make_input(
        [0.0; 3],
        full(),
        multiply(full(), [0.3, 0.1, 0.05]),
        0.4,
        HertzDrivenMode::Resting,
    );
    let error = law
        .advance_driven(
            &rest,
            HertzWrenchConfig {
                maximum_force_n: 0.01,
                ..config()
            },
            controls(),
            || false,
        )
        .unwrap_err();
    assert!(
        matches!(
            error,
            HertzDrivenError::Contact(HertzWrenchError::EffortBudget { .. })
        ),
        "{error:?}"
    );
    assert_eq!(input, before);
    assert_eq!(law.parameters().kinetic_mu, 0.5);
}

#[test]
fn rounded_zero_before_exact_stop_cannot_publish_sliding_zero_mode() {
    let law = law(0.5, 0.35);
    let stop = 0.3_f64 / 0.7;
    let immediately_before = f64::from_bits(stop.to_bits() - 1);
    assert!(immediately_before < stop);
    assert_eq!(0.3 - 0.7 * immediately_before, 0.0);
    let input = make_input(
        [0.3, 0.0, 0.0],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [0.0; 3],
        immediately_before,
        HertzDrivenMode::Sliding,
    );
    assert_eq!(
        law.advance_driven(&input, config(), controls(), || false),
        Err(HertzDrivenError::UnresolvedDirection)
    );
}
