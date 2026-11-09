//! Original driven rate laws: analytic terminal modes, unstable weakening,
//! clipping re-entry, independent pressure/time quadrature and honest events.
use core::f64::consts::PI;
use fs_tribo::FrictionLaw;
use fs_tribo::sliding_patch::contact::driven::{
    HertzDrivenControls, HertzDrivenInput, HertzDrivenMode,
};
use fs_tribo::sliding_patch::contact::transient::{HertzTransientControls, HertzTransientError};
use fs_tribo::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchParameters,
};
use fs_tribo::sliding_patch::rate_dependent::driven::{
    HertzRateDrivenControls, HertzRateDrivenError, HertzRateDrivenReport,
};
use fs_tribo::sliding_patch::rate_dependent::{HertzFrictionPatch, HertzFrictionQuadrature};

const A: f64 = 0.1;
const N: f64 = 2.0;
const K: f64 = 3.0 * PI / 16.0;
type V = [f64; 3];
type M = [[f64; 3]; 3];
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
        maximum_force_n: 100.0,
        maximum_spin_moment_nm: 10.0,
    }
}
fn controls() -> HertzRateDrivenControls {
    HertzRateDrivenControls {
        integration: HertzTransientControls {
            maximum_step_s: 0.08,
            minimum_step_s: 1e-14,
            maximum_attempts: 50_000,
            maximum_evaluations: 400_001,
            scaled_rate_tolerance_m_s: 2e-9,
            scaled_impulse_tolerance_ns: 2e-9,
            scaled_travel_tolerance_m: 8e-9,
            heat_tolerance_j: 2e-10,
            stopping_time_tolerance_s: 2e-7,
        },
        rest_scaled_rate_tolerance_m_s: 2e-10,
    }
}
fn grid(panels: usize) -> HertzFrictionQuadrature {
    HertzFrictionQuadrature {
        radial_panels: panels,
        angular_panels: panels,
    }
}
fn diagonal() -> M {
    [[1.7, 0.0, 0.0], [0.0, 1.3, 0.0], [0.0, 0.0, 2.4]]
}
fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn mul(m: M, v: V) -> V {
    m.map(|row| dot(row, v))
}
fn make_input(q: V, m: M, drive: V, dt: f64, mode: HertzDrivenMode) -> HertzDrivenInput {
    HertzDrivenInput {
        velocity_m_s: [q[0], q[1]],
        spin_rate_rad_s: q[2] / A,
        scaled_mobility_per_kg: m,
        driving_acceleration_m_s2: [drive[0], drive[1]],
        driving_spin_acceleration_rad_s2: drive[2] / A,
        dt_s: dt,
        mode,
    }
}
fn scalar(axis: usize, q: f64, demand: f64, dt: f64) -> HertzDrivenInput {
    let mut rate = [0.0; 3];
    rate[axis] = q;
    let mut drive = [0.0; 3];
    drive[axis] = diagonal()[axis][axis] * demand;
    make_input(
        rate,
        diagonal(),
        drive,
        dt,
        if q == 0.0 {
            HertzDrivenMode::Resting
        } else {
            HertzDrivenMode::Sliding
        },
    )
}
fn q(r: &HertzRateDrivenReport) -> V {
    [
        r.trajectory.velocity_after_m_s[0],
        r.trajectory.velocity_after_m_s[1],
        A * r.trajectory.spin_rate_after_rad_s,
    ]
}
fn impulse(r: &HertzRateDrivenReport) -> V {
    [
        r.trajectory.force_impulse_ns[0],
        r.trajectory.force_impulse_ns[1],
        r.trajectory.spin_angular_impulse_nms / A,
    ]
}
fn travel(r: &HertzRateDrivenReport) -> V {
    [
        r.trajectory.slip_increment_m[0],
        r.trajectory.slip_increment_m[1],
        A * r.trajectory.spin_increment_rad,
    ]
}
fn close(a: f64, b: f64, t: f64) {
    assert!(
        (a - b).abs() <= t,
        "{a:.17e} != {b:.17e}; tolerance {t:.3e}"
    );
}
fn run(law: FrictionLaw, input: HertzDrivenInput) -> HertzRateDrivenReport {
    HertzFrictionPatch::new(N, A, law)
        .unwrap()
        .advance_driven(&input, config(), grid(8), controls(), || false)
        .unwrap()
}
fn affine(mu0: f64, slope: f64) -> FrictionLaw {
    FrictionLaw::VelocityDependent {
        static_mu: 0.8,
        mu_zero: mu0,
        slope_per_speed: slope,
    }
}

#[test]
fn affine_terminal_translation_and_spin_include_right_limit_release_and_true_heat() {
    let (mu0, slope) = (0.35, 0.4);
    for axis in [0, 2] {
        let (dry, viscous) = if axis == 0 {
            (N * mu0, N * slope)
        } else {
            (N * mu0 * K, 0.4 * N * slope)
        };
        let demand = if axis == 0 { 2.0 } else { 1.2 };
        let terminal = (demand - dry) / viscous;
        let beta = diagonal()[axis][axis] * viscous;
        for initial in [0.0, 0.8] {
            for dt in [0.1, 3.0] {
                let delta = initial - terminal;
                let e = (-beta * dt).exp();
                let end = terminal + delta * e;
                let distance = terminal * dt + delta * (1.0 - e) / beta;
                let square = terminal * terminal * dt
                    + 2.0 * terminal * delta * (1.0 - e) / beta
                    + delta * delta * (1.0 - e * e) / (2.0 * beta);
                let r = run(affine(mu0, slope), scalar(axis, initial, demand, dt));
                close(q(&r)[axis], end, 8e-9);
                close(travel(&r)[axis], distance, 2e-8);
                close(
                    impulse(&r)[axis],
                    (end - initial) / diagonal()[axis][axis] - demand * dt,
                    8e-9,
                );
                close(
                    r.trajectory.dissipated_energy_j,
                    dry * distance + viscous * square,
                    1e-9,
                );
                close(r.trajectory.drive_work_j, demand * distance, 4e-8);
                close(
                    r.trajectory.kinetic_energy_change_j,
                    (end * end - initial * initial) / (2.0 * diagonal()[axis][axis]),
                    2e-8,
                );
                assert_eq!(r.trajectory.mode, HertzDrivenMode::Sliding);
                assert_eq!(r.trajectory.first_arrest_time_s, None);
                assert!(r.exact_invariant);
                assert!(!r.adaptive_error_estimates_available);
            }
        }
    }
}

#[test]
fn full_spd_viscosity_has_exact_shifted_modal_motion_work_and_no_false_arrest() {
    // C=Nb*diag(1,1,2/5), W=C^-1/2 R Lambda R^T C^-1/2.
    // z=C^1/2 q approaches z* by three independently known exponentials.
    let slope = 0.35;
    let damping = [N * slope, N * slope, 0.4 * N * slope];
    let roots = damping.map(f64::sqrt);
    let modes = [
        [2.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0],
        [2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0],
        [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
    ];
    let lambda = [0.4, 1.3, 3.0];
    let amplitudes = [0.7, -0.2, 0.4];
    let matrix: M = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            (0..3)
                .map(|k| modes[i][k] * lambda[k] * modes[j][k])
                .sum::<f64>()
                / (roots[i] * roots[j])
        })
    });
    let terminal = [0.25, -0.1, 0.4];
    let zstar: V = core::array::from_fn(|i| roots[i] * terminal[i]);
    let modal_star: V = core::array::from_fn(|k| (0..3).map(|i| modes[i][k] * zstar[i]).sum());
    let demand: V = core::array::from_fn(|i| damping[i] * terminal[i]);
    let initial: V = core::array::from_fn(|i| {
        terminal[i] + (0..3).map(|k| modes[i][k] * amplitudes[k]).sum::<f64>() / roots[i]
    });
    for dt in [0.1, 5.0, 20.0] {
        let e = lambda.map(|x| (-x * dt).exp());
        let expected: V = core::array::from_fn(|i| {
            terminal[i]
                + (0..3)
                    .map(|k| modes[i][k] * amplitudes[k] * e[k])
                    .sum::<f64>()
                    / roots[i]
        });
        let distance: V = core::array::from_fn(|i| {
            terminal[i] * dt
                + (0..3)
                    .map(|k| modes[i][k] * amplitudes[k] * (1.0 - e[k]) / lambda[k])
                    .sum::<f64>()
                    / roots[i]
        });
        let heat = dot(zstar, zstar) * dt
            + (0..3)
                .map(|k| {
                    2.0 * modal_star[k] * amplitudes[k] * (1.0 - e[k]) / lambda[k]
                        + amplitudes[k] * amplitudes[k] * (1.0 - e[k] * e[k]) / (2.0 * lambda[k])
                })
                .sum::<f64>();
        let kinetic = 0.5
            * (0..3)
                .map(|k| {
                    ((modal_star[k] + amplitudes[k] * e[k]).powi(2)
                        - (modal_star[k] + amplitudes[k]).powi(2))
                        / lambda[k]
                })
                .sum::<f64>();
        let r = run(
            affine(0.0, slope),
            make_input(
                initial,
                matrix,
                mul(matrix, demand),
                dt,
                HertzDrivenMode::Sliding,
            ),
        );
        for i in 0..3 {
            close(q(&r)[i], expected[i], 8e-9);
            close(travel(&r)[i], distance[i], 2e-8);
            close(impulse(&r)[i], -damping[i] * distance[i], 8e-9);
        }
        close(r.trajectory.dissipated_energy_j, heat, 2e-9);
        close(r.trajectory.drive_work_j, dot(demand, distance), 1e-8);
        close(r.trajectory.kinetic_energy_change_j, kinetic, 2e-9);
        assert_eq!(r.trajectory.first_arrest_time_s, None);
        assert_eq!(r.trajectory.mode, HertzDrivenMode::Sliding);
    }
    // The zero-drive limit remains asymptotic even with a nonzero static set.
    let exposure = 12.0_f64;
    let dt = exposure / (N * slope * diagonal()[0][0]);
    let r = run(affine(0.0, slope), scalar(0, 0.7, 0.0, dt));
    close(q(&r)[0], 0.7 * (-exposure).exp(), 8e-9);
    assert!(q(&r)[0] > 0.0);
    assert_eq!(r.trajectory.mode, HertzDrivenMode::Sliding);
    assert_eq!(r.trajectory.first_arrest_time_s, None);
}

#[test]
fn weakening_equilibrium_is_unstable_and_does_not_erase_static_mode_memory() {
    let mu0 = 0.6;
    let slope = -1.0;
    let demand = 0.7;
    let unstable = (N * mu0 - demand) / (N * (-slope));
    let exponent = diagonal()[0][0] * N * (-slope);
    for perturbation in [-0.01, 0.01] {
        let initial = unstable + perturbation;
        let dt = 0.4;
        let expected = unstable + perturbation * (exponent * dt).exp();
        let r = run(affine(mu0, slope), scalar(0, initial, demand, dt));
        close(q(&r)[0], expected, 8e-9);
        assert!((q(&r)[0] - unstable).abs() > 3.0 * perturbation.abs());
        assert_eq!(r.trajectory.mode, HertzDrivenMode::Sliding);
    }
    let held = run(affine(mu0, slope), scalar(0, 0.0, demand, 0.4));
    assert_eq!(held.trajectory.mode, HertzDrivenMode::Resting);
    assert_eq!(q(&held), [0.0; 3]);
    assert_eq!(held.trajectory.dissipated_energy_j, 0.0);
    // A static viscous interface may transmit force above its zero kinetic
    // capacity; the native whole-path bound must include that actual reaction.
    let held_viscous = run(affine(0.0, 0.3), scalar(0, 0.0, 1.1, 0.4));
    assert_eq!(held_viscous.trajectory.mode, HertzDrivenMode::Resting);
    assert!(held_viscous.maximum_scaled_wrench_norm_n >= 1.1);
    close(impulse(&held_viscous)[0], -1.1 * 0.4, 2e-14);
}

#[test]
fn clipped_zero_motion_reenters_actual_affine_traction_under_opposing_drive() {
    let mu0 = 0.4;
    let slope = -1.0;
    let initial = 0.5;
    let demand = -0.05;
    let mobility = diagonal()[0][0];
    let cutoff = mu0 / (-slope);
    let a = mobility * demand;
    let crossing = (cutoff - initial) / a;
    let active = 0.2;
    let k = mobility * N * (-slope);
    let unstable = (N * mu0 - demand) / (N * (-slope));
    let delta = cutoff - unstable;
    let e = (k * active).exp();
    let expected = unstable + delta * e;
    let coast_distance = initial * crossing + 0.5 * a * crossing * crossing;
    let active_distance = unstable * active + delta * (e - 1.0) / k;
    let square = unstable * unstable * active
        + 2.0 * unstable * delta * (e - 1.0) / k
        + delta * delta * (e * e - 1.0) / (2.0 * k);
    let heat = N * mu0 * active_distance + N * slope * square;
    let r = run(
        affine(mu0, slope),
        scalar(0, initial, demand, crossing + active),
    );
    close(q(&r)[0], expected, 1e-8);
    close(travel(&r)[0], coast_distance + active_distance, 2e-8);
    close(r.trajectory.dissipated_energy_j, heat, 2e-9);
    assert!(r.trajectory.dissipated_energy_j > 0.0);
    assert!(q(&r)[0] < cutoff);
    let coast = run(affine(mu0, slope), scalar(0, initial, 0.08, 0.5));
    close(q(&coast)[0], initial + mobility * 0.08 * 0.5, 2e-14);
    close(
        travel(&coast)[0],
        initial * 0.5 + 0.5 * mobility * 0.08 * 0.5 * 0.5,
        2e-14,
    );
    assert_eq!(coast.trajectory.dissipated_energy_j, 0.0);
}

// Independent Simpson integration of Hertz radial pressure and direct
// Stribeck coefficients; neither production Gauss nodes nor baseline/correction
// decomposition enter. RK4 advances rate, travel and actual pressure power.
fn spin_reference(initial: f64, demand: f64, dt: f64, radial: usize, steps: usize) -> [f64; 3] {
    let moment = |q: f64| {
        let dtheta = (PI / 2.0) / radial as f64;
        let mut sum = 0.0;
        for i in 0..=radial {
            let (s, c) = (i as f64 * dtheta).sin_cos();
            let speed = q * s;
            let mu = 0.1 + 0.8 * (-(speed / 0.08).powi(2)).exp() + 0.05 * speed;
            let weight = if i == 0 || i == radial {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            sum += weight * mu * s * s * c * c;
        }
        -N * dtheta * sum
    };
    let rhs = |state: [f64; 3]| {
        let w = moment(state[0]);
        [diagonal()[2][2] * (demand + w), state[0], -w * state[0]]
    };
    let h = dt / steps as f64;
    let mut state = [initial, 0.0, 0.0];
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
fn driven_stribeck_spin_refines_pressure_quadrature_independently_of_time_accuracy() {
    let (initial, demand, dt) = (0.9, 0.25, 0.3);
    let patch = HertzFrictionPatch::new(
        N,
        A,
        FrictionLaw::Stribeck {
            static_mu: 0.9,
            kinetic_mu: 0.1,
            characteristic_speed: 0.08,
            viscous_per_speed: 0.05,
        },
    )
    .unwrap();
    let expected = spin_reference(initial, demand, dt, 512, 512);
    let coarse = spin_reference(initial, demand, dt, 256, 256);
    for i in 0..3 {
        close(expected[i], coarse[i], 3e-10);
    }
    let mut errors = Vec::new();
    for panels in [4, 8, 16, 32] {
        let r = patch
            .advance_driven(
                &scalar(2, initial, demand, dt),
                config(),
                grid(panels),
                controls(),
                || false,
            )
            .unwrap();
        errors.push((q(&r)[2] - expected[0]).abs());
        if panels == 32 {
            close(q(&r)[2], expected[0], 5e-9);
            close(travel(&r)[2], expected[1], 1e-8);
            close(r.trajectory.dissipated_energy_j, expected[2], 1e-9);
            close(r.trajectory.drive_work_j, demand * expected[1], 3e-9);
        }
    }
    for pair in errors.windows(2) {
        assert!(pair[1] < 0.15 * pair[0], "{errors:?}");
    }
    eprintln!("driven Stribeck pressure errors4/8/16/32: {errors:?}");
}

#[test]
fn coulomb_delegation_and_unloaded_original_curve_bounds_remain_honest() {
    let law = FrictionLaw::Coulomb {
        static_mu: 0.8,
        kinetic_mu: 0.5,
    };
    let input = scalar(0, 0.3, 0.0, 1.0);
    let r = run(law, input);
    let c = controls();
    let i = c.integration;
    let exact = HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: N,
        contact_radius_m: A,
        static_mu: 0.8,
        kinetic_mu: 0.5,
    })
    .unwrap()
    .advance_driven(
        &input,
        config(),
        HertzDrivenControls {
            maximum_step_s: i.maximum_step_s,
            minimum_step_s: i.minimum_step_s,
            maximum_attempts: i.maximum_attempts,
            maximum_evaluations: i.maximum_evaluations,
            scaled_rate_tolerance_m_s: i.scaled_rate_tolerance_m_s,
            rest_scaled_rate_tolerance_m_s: c.rest_scaled_rate_tolerance_m_s,
        },
        || false,
    )
    .unwrap();
    assert_eq!(r.trajectory, exact);
    assert!(!r.adaptive_error_estimates_available);
    assert!(r.trajectory.arrest_time_is_exact);
    for law in [
        affine(0.4, 0.3),
        FrictionLaw::Stribeck {
            static_mu: 0.8,
            kinetic_mu: 0.2,
            characteristic_speed: 0.1,
            viscous_per_speed: 0.4,
        },
    ] {
        let patch = HertzFrictionPatch::new(0.0, A, law.clone()).unwrap();
        let input = make_input(
            [0.5, -0.2, 0.3],
            diagonal(),
            [-0.1, 0.3, 0.2],
            0.8,
            HertzDrivenMode::Sliding,
        );
        let r = patch
            .advance_driven(&input, config(), grid(4), controls(), || false)
            .unwrap();
        assert!(r.local_speed_bound_m_s.is_finite());
        assert_eq!(r.maximum_scaled_wrench_norm_n, 0.0);
        assert_eq!(r.trajectory.dissipated_energy_j, 0.0);
        let before = [0.5, -0.2, 0.3];
        let acceleration = [-0.1, 0.3, 0.2];
        for i in 0..3 {
            close(q(&r)[i], before[i] + input.dt_s * acceleration[i], 2e-14);
        }
        for i in 0..=100 {
            let speed = r.local_speed_bound_m_s * i as f64 / 100.0;
            let mu = if speed == 0.0 {
                law.zero_speed_kinetic_coefficient().unwrap()
            } else {
                law.kinetic_coefficient(speed).unwrap()
            };
            assert!(r.kinetic_coefficient_bounds[0] <= mu && mu <= r.kinetic_coefficient_bounds[1]);
        }
    }
}

#[test]
fn phase_splitting_and_signed_frame_change_preserve_rate_curve_work() {
    let patch = HertzFrictionPatch::new(N, A, affine(0.3, 0.4)).unwrap();
    let matrix = [[2.0, 0.3, 0.4], [0.3, 1.7, -0.2], [0.4, -0.2, 3.5]];
    let initial = [0.8, -0.3, 0.5];
    let drive = [0.2, -0.05, 0.1];
    let input = make_input(initial, matrix, drive, 0.18, HertzDrivenMode::Sliding);
    let whole = patch
        .advance_driven(&input, config(), grid(1), controls(), || false)
        .unwrap();
    let first = patch
        .advance_driven(
            &HertzDrivenInput {
                dt_s: 0.07,
                ..input
            },
            config(),
            grid(1),
            controls(),
            || false,
        )
        .unwrap();
    let second = patch
        .advance_driven(
            &make_input(q(&first), matrix, drive, 0.11, first.trajectory.mode),
            config(),
            grid(1),
            controls(),
            || false,
        )
        .unwrap();
    for i in 0..3 {
        close(q(&second)[i], q(&whole)[i], 8e-9);
        close(
            travel(&first)[i] + travel(&second)[i],
            travel(&whole)[i],
            2e-8,
        );
    }
    close(
        first.trajectory.dissipated_energy_j + second.trajectory.dissipated_energy_j,
        whole.trajectory.dissipated_energy_j,
        2e-9,
    );
    close(
        first.trajectory.drive_work_j + second.trajectory.drive_work_j,
        whole.trajectory.drive_work_j,
        2e-9,
    );
    let transform = |q: V| [-q[1], q[0], q[2]];
    let permutation = [1, 0, 2];
    let sign = [-1.0, 1.0, 1.0];
    let rotated: M = core::array::from_fn(|i| {
        core::array::from_fn(|j| sign[i] * sign[j] * matrix[permutation[i]][permutation[j]])
    });
    let r = patch
        .advance_driven(
            &make_input(
                transform(initial),
                rotated,
                transform(drive),
                0.18,
                HertzDrivenMode::Sliding,
            ),
            config(),
            grid(1),
            controls(),
            || false,
        )
        .unwrap();
    let expected = transform(q(&whole));
    for i in 0..3 {
        close(q(&r)[i], expected[i], 8e-9);
    }
    close(
        r.trajectory.dissipated_energy_j,
        whole.trajectory.dissipated_energy_j,
        2e-9,
    );
    close(
        r.trajectory.drive_work_j,
        whole.trajectory.drive_work_j,
        2e-9,
    );
    assert_eq!(
        whole,
        patch
            .advance_driven(&input, config(), grid(1), controls(), || false)
            .unwrap()
    );
}

#[test]
fn unresolved_arrest_cancellation_and_budgets_never_publish_candidate_mode() {
    let patch = HertzFrictionPatch::new(N, A, affine(0.4, 0.3)).unwrap();
    let snapshot = patch.clone();
    // Non-invariant curved arrest remains explicitly outside the analytic
    // event owner. Scalar affine braking is now a resolved physical event.
    let matrix = [[2.0, 0.3, 0.4], [0.3, 1.7, -0.2], [0.4, -0.2, 3.5]];
    let input = make_input(
        [0.1, -0.04, 0.06],
        matrix,
        [0.0; 3],
        1.0,
        HertzDrivenMode::Sliding,
    );
    let before = input;
    assert!(matches!(
        patch.advance_driven(&input, config(), grid(1), controls(), || false),
        Err(HertzRateDrivenError::UnresolvedDrivenEvent { .. })
    ));
    let moving = make_input(
        [0.8, -0.3, 0.5],
        matrix,
        [0.2, -0.05, 0.1],
        0.1,
        HertzDrivenMode::Sliding,
    );
    assert_eq!(
        patch.advance_driven(&moving, config(), grid(1), controls(), || true),
        Err(HertzRateDrivenError::Transient(
            HertzTransientError::Cancelled
        ))
    );
    let mut polls = 0;
    assert_eq!(
        patch.advance_driven(&moving, config(), grid(1), controls(), || {
            polls += 1;
            polls > 10
        }),
        Err(HertzRateDrivenError::Transient(
            HertzTransientError::Cancelled
        ))
    );
    let c = HertzRateDrivenControls {
        integration: HertzTransientControls {
            maximum_evaluations: 1,
            ..controls().integration
        },
        ..controls()
    };
    assert_eq!(
        patch.advance_driven(&moving, config(), grid(1), c, || false),
        Err(HertzRateDrivenError::Transient(
            HertzTransientError::EvaluationBudget
        ))
    );
    let c = HertzRateDrivenControls {
        integration: HertzTransientControls {
            maximum_attempts: 1,
            ..controls().integration
        },
        ..controls()
    };
    assert_eq!(
        patch.advance_driven(&moving, config(), grid(1), c, || false),
        Err(HertzRateDrivenError::Transient(
            HertzTransientError::AttemptBudget
        ))
    );
    assert_eq!(patch, snapshot);
    assert_eq!(input, before);
}

#[path = "hertz_rate_driven/invariant.rs"]
mod invariant;
