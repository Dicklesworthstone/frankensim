//! Original speed laws in finite-body Hertz recoil: analytic rate/pressure
//! moments, independent time integration, spatial refinement and honest arrest.
use core::f64::consts::PI;
use fs_tribo::FrictionLaw;
use fs_tribo::sliding_patch::contact::relaxation::HertzRelaxationInput;
use fs_tribo::sliding_patch::contact::transient::{HertzTransientControls, HertzTransientError};
use fs_tribo::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchParameters,
};
use fs_tribo::sliding_patch::rate_dependent::transient::{
    HertzRateTransientError, HertzRateTransientReport,
};
use fs_tribo::sliding_patch::rate_dependent::{HertzFrictionPatch, HertzFrictionQuadrature};

const A: f64 = 0.1;
const N: f64 = 2.0;
const K: f64 = 3.0 * PI / 16.0;
type V = [f64; 3];
type M = [[f64; 3]; 3];
fn config() -> HertzWrenchConfig {
    HertzWrenchConfig {
        maximum_static_iterations: 1,
        maximum_iterations: 1,
        maximum_evaluations: 1,
        force_tolerance_n: 1e-11,
        moment_tolerance_nm: 1e-12,
        displacement_tolerance_m: 1e-12,
        energy_tolerance_j: 2e-9,
        relative_matrix_tolerance: 1e-13,
        maximum_force_n: 100.0,
        maximum_spin_moment_nm: 100.0,
    }
}
fn controls() -> HertzTransientControls {
    HertzTransientControls {
        maximum_step_s: 0.08,
        minimum_step_s: 1e-16,
        maximum_attempts: 50_000,
        maximum_evaluations: 400_001,
        scaled_rate_tolerance_m_s: 2e-9,
        scaled_impulse_tolerance_ns: 2e-9,
        scaled_travel_tolerance_m: 8e-9,
        heat_tolerance_j: 2e-10,
        stopping_time_tolerance_s: 2e-7,
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
fn scalar_input(axis: usize, q: f64, dt: f64) -> HertzRelaxationInput {
    let mut rate = [0.0; 3];
    rate[axis] = q;
    HertzRelaxationInput {
        velocity_m_s: [rate[0], rate[1]],
        spin_rate_rad_s: rate[2] / A,
        scaled_mobility_per_kg: diagonal(),
        dt_s: dt,
    }
}
fn q(r: &HertzRateTransientReport) -> V {
    [
        r.trajectory.velocity_after_m_s[0],
        r.trajectory.velocity_after_m_s[1],
        A * r.trajectory.spin_rate_after_rad_s,
    ]
}
fn impulse(r: &HertzRateTransientReport) -> V {
    [
        r.trajectory.force_impulse_ns[0],
        r.trajectory.force_impulse_ns[1],
        r.trajectory.spin_angular_impulse_nms / A,
    ]
}
fn travel(r: &HertzRateTransientReport) -> V {
    [
        r.trajectory.slip_increment_m[0],
        r.trajectory.slip_increment_m[1],
        A * r.trajectory.spin_increment_rad,
    ]
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance {tolerance:e}"
    );
}
fn run(law: FrictionLaw, input: HertzRelaxationInput) -> HertzRateTransientReport {
    HertzFrictionPatch::new(N, A, law)
        .unwrap()
        .advance_transient(&input, config(), grid(8), controls(), || false)
        .unwrap()
}
fn scalar_check(r: &HertzRateTransientReport, axis: usize, initial: f64, end: f64, distance: f64) {
    let mobility = diagonal()[axis][axis];
    close(q(r)[axis], end, 5e-9);
    close(impulse(r)[axis], (end - initial) / mobility, 5e-9);
    close(travel(r)[axis], distance, 2e-8);
    close(
        r.trajectory.dissipated_energy_j,
        (initial * initial - end * end) / (2.0 * mobility),
        8e-10,
    );
    assert!(r.trajectory.work_residual_j.abs() <= config().energy_tolerance_j);
    assert!(r.trajectory.dissipated_energy_j >= 0.0);
}

#[test]
fn affine_translation_and_spin_follow_exact_exponential_motion_and_pressure_moments() {
    // For pure translation q'=-W*N*(mu0+b*q). For pure spin, the actual
    // Hertz first/second pressure moments give q'=-W*N*(k*mu0+(2/5)*b*q).
    // Evaluating an effective coefficient only at the centre misses this drag.
    for axis in [0, 2] {
        for slope in [-0.12, 0.35] {
            let (initial, mu0) = (0.7, 0.4);
            let alpha = diagonal()[axis][axis] * N * mu0 * if axis == 0 { 1.0 } else { K };
            let beta = diagonal()[axis][axis] * N * slope * if axis == 0 { 1.0 } else { 0.4 };
            let stop = (1.0 + beta * initial / alpha).ln() / beta;
            for dt in [0.2 * stop, 2.0 * stop] {
                let t = dt.min(stop);
                let decay = (-beta * t).exp();
                let end = if dt >= stop {
                    0.0
                } else {
                    (initial + alpha / beta) * decay - alpha / beta
                };
                let distance = (initial + alpha / beta) * (1.0 - decay) / beta - alpha * t / beta;
                let r = run(
                    FrictionLaw::VelocityDependent {
                        static_mu: 0.8,
                        mu_zero: mu0,
                        slope_per_speed: slope,
                    },
                    scalar_input(axis, initial, dt),
                );
                scalar_check(&r, axis, initial, end, distance);
                assert!(
                    !r.trajectory.exact_ray,
                    "a rate-dependent ray is not linear Coulomb decay"
                );
                if dt >= stop {
                    let bracket = r.trajectory.stopping_time_bracket_s.unwrap();
                    assert!(
                        bracket[0] - 2e-8 <= stop && stop <= bracket[1] + 2e-8,
                        "{bracket:?}, exact {stop}"
                    );
                    assert_eq!(q(&r), [0.0; 3]);
                } else {
                    assert!(r.trajectory.stopping_time_bracket_s.is_none());
                }
            }
        }
    }
}

#[test]
fn zero_asymptote_hyperbolic_and_shifted_power_still_stop_in_finite_time() {
    let (initial, mu, c, exponent) = (0.8_f64, 0.7, 0.3, 0.4);
    let acceleration = diagonal()[0][0] * N * mu;
    for n in [1.0, exponent] {
        let law = if n == 1.0 {
            FrictionLaw::Hyperbolic {
                static_mu: mu,
                kinetic_mu: 0.0,
                characteristic_speed: c,
            }
        } else {
            FrictionLaw::ShiftedPower {
                static_mu: mu,
                kinetic_mu: 0.0,
                characteristic_speed: c,
                exponent: n,
            }
        };
        let z0 = 1.0 + initial / c;
        let stop = c * (z0.powf(n + 1.0) - 1.0) / (acceleration * (n + 1.0));
        for dt in [0.35 * stop, 2.0 * stop] {
            let t = dt.min(stop);
            let z = if dt >= stop {
                1.0
            } else {
                (z0.powf(n + 1.0) - acceleration * (n + 1.0) * t / c).powf(1.0 / (n + 1.0))
            };
            let end = c * (z - 1.0);
            let distance = c * c / acceleration
                * ((z0.powf(n + 2.0) - z.powf(n + 2.0)) / (n + 2.0)
                    - (z0.powf(n + 1.0) - z.powf(n + 1.0)) / (n + 1.0));
            let r = run(law.clone(), scalar_input(0, initial, dt));
            scalar_check(&r, 0, initial, end, distance);
            assert!(r.kinetic_coefficient_bounds[0] > 0.0);
            if dt >= stop {
                let bracket = r.trajectory.stopping_time_bracket_s.unwrap();
                assert!(bracket[0] - 2e-8 <= stop && stop <= bracket[1] + 2e-8);
            }
        }
    }
}

fn add(a: V, b: V) -> V {
    core::array::from_fn(|i| a[i] + b[i])
}
fn scale(a: V, s: f64) -> V {
    a.map(|x| s * x)
}
fn mul(m: M, v: V) -> V {
    m.map(|row| row[0] * v[0] + row[1] * v[1] + row[2] * v[2])
}
fn reference(
    patch: &HertzFrictionPatch,
    input: HertzRelaxationInput,
    steps: usize,
) -> (V, V, V, f64) {
    // Independent fixed-grid RK4; no production error estimator, mass metric,
    // step controller, tail or integration helper is used.
    let h = input.dt_s / steps as f64;
    let mut q = [
        input.velocity_m_s[0],
        input.velocity_m_s[1],
        A * input.spin_rate_rad_s,
    ];
    let (mut j, mut x, mut heat) = ([0.0; 3], [0.0; 3], 0.0);
    let sample = |q: V| {
        let r = patch.evaluate([q[0], q[1]], q[2] / A, grid(8)).unwrap();
        (
            [r.force_n[0], r.force_n[1], r.spin_moment_nm / A],
            r.dissipated_power_w,
        )
    };
    for _ in 0..steps {
        let (a, pa) = sample(q);
        let qb = add(q, scale(mul(input.scaled_mobility_per_kg, a), 0.5 * h));
        let (b, pb) = sample(qb);
        let qc = add(q, scale(mul(input.scaled_mobility_per_kg, b), 0.5 * h));
        let (c, pc) = sample(qc);
        let qd = add(q, scale(mul(input.scaled_mobility_per_kg, c), h));
        let (d, pd) = sample(qd);
        let dj = core::array::from_fn(|i| h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
        x = add(
            x,
            core::array::from_fn(|i| h * (q[i] + 2.0 * qb[i] + 2.0 * qc[i] + qd[i]) / 6.0),
        );
        q = add(q, mul(input.scaled_mobility_per_kg, dj));
        j = add(j, dj);
        heat += h * (pa + 2.0 * pb + 2.0 * pc + pd) / 6.0;
    }
    (q, j, x, heat)
}
fn mixed_input() -> HertzRelaxationInput {
    HertzRelaxationInput {
        velocity_m_s: [0.8, -0.3],
        spin_rate_rad_s: 5.0,
        scaled_mobility_per_kg: [[2.0, 0.3, 0.4], [0.3, 1.7, -0.2], [0.4, -0.2, 3.5]],
        dt_s: 0.15,
    }
}

#[test]
fn full_mobility_rate_trajectory_matches_independent_refinement_and_split_history() {
    let patch = HertzFrictionPatch::new(
        N,
        A,
        FrictionLaw::VelocityDependent {
            static_mu: 0.8,
            mu_zero: 0.4,
            slope_per_speed: 0.3,
        },
    )
    .unwrap();
    let input = mixed_input();
    let coarse = reference(&patch, input, 512);
    let fine = reference(&patch, input, 1024);
    let r = patch
        .advance_transient(&input, config(), grid(1), controls(), || false)
        .unwrap();
    for i in 0..3 {
        close(coarse.0[i], fine.0[i], 2e-13);
        close(q(&r)[i], fine.0[i], 4e-9);
        close(impulse(&r)[i], fine.1[i], 3e-9);
        close(travel(&r)[i], fine.2[i], 8e-9);
    }
    close(r.trajectory.dissipated_energy_j, fine.3, 8e-10);
    let first = patch
        .advance_transient(
            &HertzRelaxationInput {
                dt_s: 0.06,
                ..input
            },
            config(),
            grid(1),
            controls(),
            || false,
        )
        .unwrap();
    let second = patch
        .advance_transient(
            &HertzRelaxationInput {
                velocity_m_s: first.trajectory.velocity_after_m_s,
                spin_rate_rad_s: first.trajectory.spin_rate_after_rad_s,
                dt_s: 0.09,
                ..input
            },
            config(),
            grid(1),
            controls(),
            || false,
        )
        .unwrap();
    for i in 0..3 {
        close(q(&second)[i], q(&r)[i], 5e-9);
        close(
            impulse(&first)[i] + impulse(&second)[i],
            impulse(&r)[i],
            5e-9,
        );
        close(travel(&first)[i] + travel(&second)[i], travel(&r)[i], 1e-8);
    }
    close(
        first.trajectory.dissipated_energy_j + second.trajectory.dissipated_energy_j,
        r.trajectory.dissipated_energy_j,
        1e-9,
    );
}

// Independent Simpson integration of the ORIGINAL Hertz radial pressure, with
// the Stribeck curve evaluated directly. Production uses a Coulomb baseline
// and Gauss correction; this oracle shares neither that decomposition nor nodes.
fn spin_reference(initial: f64, dt: f64, radial: usize, steps: usize) -> (f64, f64) {
    let h = dt / steps as f64;
    let moment = |q: f64| {
        let dtheta = (PI / 2.0) / radial as f64;
        let mut sum = 0.0;
        for i in 0..=radial {
            let theta = i as f64 * dtheta;
            let (s, c) = theta.sin_cos();
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
    let mut q = initial;
    let mut distance = 0.0;
    for _ in 0..steps {
        let a = diagonal()[2][2] * moment(q);
        let qb = q + 0.5 * h * a;
        let b = diagonal()[2][2] * moment(qb);
        let qc = q + 0.5 * h * b;
        let c = diagonal()[2][2] * moment(qc);
        let qd = q + h * c;
        let d = diagonal()[2][2] * moment(qd);
        distance += h * (q + 2.0 * qb + 2.0 * qc + qd) / 6.0;
        q += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    (q, distance)
}

#[test]
fn nonlinear_spin_converges_to_independent_pressure_and_time_quadrature() {
    let (initial, dt) = (0.9, 0.3);
    let law = FrictionLaw::Stribeck {
        static_mu: 0.9,
        kinetic_mu: 0.1,
        characteristic_speed: 0.08,
        viscous_per_speed: 0.05,
    };
    let patch = HertzFrictionPatch::new(N, A, law.clone()).unwrap();
    let expected = spin_reference(initial, dt, 512, 512);
    let coarse_reference = spin_reference(initial, dt, 256, 256);
    close(expected.0, coarse_reference.0, 2e-10);
    close(expected.1, coarse_reference.1, 2e-10);
    let mut errors = Vec::new();
    // This sharply localized low-speed Stribeck core needs resolved radial
    // pressure quadrature as well as a converged physical-time integration.
    for panels in [4, 8, 16, 32] {
        let r = patch
            .advance_transient(
                &scalar_input(2, initial, dt),
                config(),
                grid(panels),
                controls(),
                || false,
            )
            .unwrap();
        errors.push((q(&r)[2] - expected.0).abs());
        if panels == 32 {
            scalar_check(&r, 2, initial, expected.0, expected.1);
            for i in 0..=100 {
                let speed = r.local_speed_bound_m_s * i as f64 / 100.0;
                let mu = if i == 0 {
                    law.zero_speed_kinetic_coefficient().unwrap()
                } else {
                    law.kinetic_coefficient(speed).unwrap()
                };
                assert!(
                    mu >= r.kinetic_coefficient_bounds[0] && mu <= r.kinetic_coefficient_bounds[1]
                );
            }
        }
    }
    assert!(
        errors[0] > errors[1] && errors[1] > errors[2],
        "spatial errors {errors:?}"
    );
    assert!(errors[3] < errors[0] / 1000.0, "spatial errors {errors:?}");
    eprintln!("Stribeck Hertz trajectory radial-panel errors {errors:?}");
}

#[test]
fn zero_traction_coasts_but_noncoercive_drag_never_claims_finite_arrest() {
    let viscous = FrictionLaw::VelocityDependent {
        static_mu: 0.8,
        mu_zero: 0.0,
        slope_per_speed: 0.3,
    };
    let clipped = FrictionLaw::VelocityDependent {
        static_mu: 0.8,
        mu_zero: 0.4,
        slope_per_speed: -1.0,
    };
    for law in [viscous.clone(), clipped.clone()] {
        let patch = HertzFrictionPatch::new(N, A, law).unwrap();
        let input = scalar_input(2, 0.7, 0.1);
        let r = patch
            .advance_transient(&input, config(), grid(8), controls(), || false)
            .unwrap();
        assert_eq!(r.kinetic_coefficient_bounds[0], 0.0);
        assert!(r.trajectory.stopping_time_bracket_s.is_none());
        assert!(q(&r)[2] > 0.0 && q(&r)[2] < 0.7);
        assert!(r.trajectory.dissipated_energy_j > 0.0);
    }
    // All local speeds exceed the clipping cutoff even with spin: exact
    // autonomous coasting has no moving-to-rest event and generates no heat.
    let patch = HertzFrictionPatch::new(N, A, clipped).unwrap();
    let input = HertzRelaxationInput {
        velocity_m_s: [1.2, 0.0],
        spin_rate_rad_s: 2.0,
        dt_s: 5.0,
        scaled_mobility_per_kg: diagonal(),
    };
    let coast = patch
        .advance_transient(&input, config(), grid(8), controls(), || false)
        .unwrap();
    assert_eq!(coast.trajectory.velocity_after_m_s, input.velocity_m_s);
    assert_eq!(
        coast.trajectory.spin_rate_after_rad_s,
        input.spin_rate_rad_s
    );
    assert_eq!(coast.trajectory.slip_increment_m, [6.0, 0.0]);
    assert_eq!(coast.trajectory.spin_increment_rad, 10.0);
    assert!(coast.trajectory.stopping_time_bracket_s.is_none());
    assert_eq!(coast.maximum_scaled_wrench_norm_n, 0.0);
    assert_eq!(coast.trajectory.dissipated_energy_j, 0.0);
    for (normal, rate) in [(N, 0.0), (0.0, 0.7)] {
        let patch = HertzFrictionPatch::new(normal, A, viscous.clone()).unwrap();
        let r = patch
            .advance_transient(
                &scalar_input(0, rate, 1.0),
                config(),
                grid(2),
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(r.trajectory.velocity_after_m_s, [rate, 0.0]);
        assert_eq!(r.trajectory.evaluations, 1);
        assert_eq!(r.trajectory.stopping_time_bracket_s.is_some(), rate == 0.0);
    }
}

#[test]
fn original_coulomb_fast_path_and_frictionless_coast_are_preserved() {
    for mu in [0.0, 0.4] {
        let patch = HertzFrictionPatch::new(
            N,
            A,
            FrictionLaw::Coulomb {
                static_mu: 0.8,
                kinetic_mu: mu,
            },
        )
        .unwrap();
        let original = HertzWrenchContact::new(HertzWrenchParameters {
            normal_force_n: N,
            contact_radius_m: A,
            static_mu: 0.8,
            kinetic_mu: mu,
        })
        .unwrap();
        for input in [scalar_input(0, 0.7, 2.0), mixed_input()] {
            let actual = patch
                .advance_transient(&input, config(), grid(1), controls(), || false)
                .unwrap();
            let expected = original
                .advance_transient(&input, config(), controls(), || false)
                .unwrap();
            assert_eq!(actual.trajectory, expected);
            assert_eq!(actual.kinetic_coefficient_bounds, [mu; 2]);
        }
    }
}

#[test]
fn cancellation_effort_and_numerical_budgets_preserve_the_original_curve() {
    for law in [
        FrictionLaw::ShiftedPower {
            static_mu: 0.8,
            kinetic_mu: 0.2,
            characteristic_speed: 0.3,
            exponent: 0.4,
        },
        FrictionLaw::VelocityDependent {
            static_mu: 0.8,
            mu_zero: 0.0,
            slope_per_speed: 0.3,
        },
        FrictionLaw::VelocityDependent {
            static_mu: 0.8,
            mu_zero: 0.4,
            slope_per_speed: -1.0,
        },
    ] {
        let patch = HertzFrictionPatch::new(N, A, law).unwrap();
        let snapshot = patch.clone();
        let input = mixed_input();
        let baseline = patch
            .advance_transient(&input, config(), grid(8), controls(), || false)
            .unwrap();
        let mut polls = 0;
        let cancelled = patch.advance_transient(&input, config(), grid(8), controls(), || {
            polls += 1;
            polls >= 17
        });
        assert_eq!(
            cancelled,
            Err(HertzRateTransientError::Transient(
                HertzTransientError::Cancelled
            ))
        );
        let exhausted = patch.advance_transient(
            &input,
            config(),
            grid(8),
            HertzTransientControls {
                maximum_evaluations: 1,
                ..controls()
            },
            || false,
        );
        assert_eq!(
            exhausted,
            Err(HertzRateTransientError::Transient(
                HertzTransientError::EvaluationBudget
            ))
        );
        let effort = patch.advance_transient(
            &input,
            HertzWrenchConfig {
                maximum_force_n: 0.0,
                ..config()
            },
            grid(8),
            controls(),
            || false,
        );
        assert_eq!(
            effort,
            Err(HertzRateTransientError::Transient(
                HertzTransientError::EffortBudget
            ))
        );
        assert_eq!(patch, snapshot);
        assert_eq!(
            patch
                .advance_transient(&input, config(), grid(8), controls(), || false)
                .unwrap(),
            baseline
        );
    }
}

#[test]
fn pure_viscosity_retains_exponential_translation_and_spin_without_finite_arrest() {
    let slope = 0.3;
    let law = FrictionLaw::VelocityDependent {
        static_mu: 0.8,
        mu_zero: 0.0,
        slope_per_speed: slope,
    };
    let patch = HertzFrictionPatch::new(N, A, law).unwrap();
    let bounded = HertzTransientControls {
        maximum_attempts: 5_000,
        maximum_evaluations: 60_001,
        ..controls()
    };
    for axis in [0, 2] {
        let initial = 0.7;
        let decay_rate = N * slope * diagonal()[axis][axis] * if axis == 0 { 1.0 } else { 0.4 };
        for exposure in [0.01, 2.0, 25.0] {
            let dt = exposure / decay_rate;
            let end = initial * (-exposure).exp();
            let distance = initial * (-(-exposure).exp_m1()) / decay_rate;
            let r = patch
                .advance_transient(
                    &scalar_input(axis, initial, dt),
                    config(),
                    grid(1),
                    bounded,
                    || false,
                )
                .unwrap();
            scalar_check(&r, axis, initial, end, distance);
            assert!(
                q(&r)[axis] > 0.0,
                "viscosity must not manufacture finite arrest"
            );
            assert!(r.trajectory.stopping_time_bracket_s.is_none());
            assert_eq!(r.trajectory.stopping_tail_heat_j, 0.0);
            assert_eq!(r.trajectory.integrated_until_s, dt);
            assert!(
                r.trajectory.accepted_steps + r.trajectory.rejected_steps
                    < bounded.maximum_attempts
            );
        }
    }
}

#[test]
fn full_mobility_viscosity_matches_independent_modal_exponentials_and_heat() {
    // Construct W=D^-1/2 R Lambda R^T D^-1/2 using known orthonormal
    // columns. Thus z=D^1/2 q has independently known modal exponential
    // motion. No numerical eigenproblem or production integration is used
    // to obtain endpoints, impulse, travel or integrated physical power.
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
    let initial: V = core::array::from_fn(|i| {
        (0..3).map(|k| modes[i][k] * amplitudes[k]).sum::<f64>() / roots[i]
    });
    let patch = HertzFrictionPatch::new(
        N,
        A,
        FrictionLaw::VelocityDependent {
            static_mu: 0.8,
            mu_zero: 0.0,
            slope_per_speed: slope,
        },
    )
    .unwrap();
    let bounded = HertzTransientControls {
        maximum_attempts: 5_000,
        maximum_evaluations: 60_001,
        ..controls()
    };
    for dt in [0.1, 40.0] {
        let expected: V = core::array::from_fn(|i| {
            (0..3)
                .map(|k| modes[i][k] * amplitudes[k] * (-lambda[k] * dt).exp())
                .sum::<f64>()
                / roots[i]
        });
        let distance: V = core::array::from_fn(|i| {
            (0..3)
                .map(|k| modes[i][k] * amplitudes[k] * (-(-lambda[k] * dt).exp_m1()) / lambda[k])
                .sum::<f64>()
                / roots[i]
        });
        let heat = (0..3)
            .map(|k| {
                amplitudes[k] * amplitudes[k] * (-(-2.0 * lambda[k] * dt).exp_m1())
                    / (2.0 * lambda[k])
            })
            .sum::<f64>();
        let input = HertzRelaxationInput {
            velocity_m_s: [initial[0], initial[1]],
            spin_rate_rad_s: initial[2] / A,
            scaled_mobility_per_kg: matrix,
            dt_s: dt,
        };
        let r = patch
            .advance_transient(&input, config(), grid(1), bounded, || false)
            .unwrap();
        for i in 0..3 {
            close(q(&r)[i], expected[i], 4e-9);
            close(travel(&r)[i], distance[i], 1e-8);
            close(impulse(&r)[i], -damping[i] * distance[i], 4e-9);
        }
        close(r.trajectory.dissipated_energy_j, heat, 8e-10);
        assert!(q(&r).into_iter().all(|v| v > 0.0));
        assert!(r.trajectory.stopping_time_bracket_s.is_none());
        assert_eq!(r.trajectory.stopping_tail_heat_j, 0.0);
        assert!(r.trajectory.work_residual_j.abs() <= config().energy_tolerance_j);
        eprintln!(
            "viscous full-mobility dt={dt}: {} accepted, {} rejected, {} evaluations",
            r.trajectory.accepted_steps, r.trajectory.rejected_steps, r.trajectory.evaluations
        );
    }
}

#[test]
fn clipped_initial_envelope_enters_original_coercive_domain_and_stops_at_the_analytic_time() {
    // Actual initial translation lies below the clipping cutoff, while the
    // conservative whole-mobility admission ball extends above it. The former
    // global-coercivity admission refused this real finite-arrest trajectory.
    let (mu0, slope, initial) = (0.4_f64, -1.0, 0.35);
    let alpha = diagonal()[0][0] * N * mu0;
    let beta = diagonal()[0][0] * N * slope;
    let stop = (1.0 + beta * initial / alpha).ln() / beta;
    let patch = HertzFrictionPatch::new(
        N,
        A,
        FrictionLaw::VelocityDependent {
            static_mu: 0.8,
            mu_zero: mu0,
            slope_per_speed: slope,
        },
    )
    .unwrap();
    for dt in [0.2 * stop, 2.0 * stop] {
        let t = dt.min(stop);
        let decay = (-beta * t).exp();
        let end = if dt >= stop {
            0.0
        } else {
            (initial + alpha / beta) * decay - alpha / beta
        };
        let distance = (initial + alpha / beta) * (1.0 - decay) / beta - alpha * t / beta;
        let r = patch
            .advance_transient(
                &scalar_input(0, initial, dt),
                config(),
                grid(1),
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(r.kinetic_coefficient_bounds[0], 0.0);
        assert!(r.local_speed_bound_m_s > mu0 / -slope);
        scalar_check(&r, 0, initial, end, distance);
        if dt >= stop {
            let bracket = r.trajectory.stopping_time_bracket_s.unwrap();
            assert!(
                bracket[0] - 2e-8 <= stop && stop <= bracket[1] + 2e-8,
                "{bracket:?}, exact {stop}"
            );
            assert_eq!(q(&r), [0.0; 3]);
            assert!(r.trajectory.stopping_tail_heat_j > 0.0);
        } else {
            assert!(r.trajectory.stopping_time_bracket_s.is_none());
        }
    }
    let coast = patch
        .advance_transient(
            &scalar_input(0, 0.5, 10.0),
            config(),
            grid(1),
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(q(&coast), [0.5, 0.0, 0.0]);
    assert_eq!(coast.trajectory.slip_increment_m, [5.0, 0.0]);
    assert!(coast.trajectory.stopping_time_bracket_s.is_none());
    assert_eq!(coast.trajectory.dissipated_energy_j, 0.0);
}

fn clipped_spin_reference(initial: f64, dt: f64, radial: usize, steps: usize) -> (f64, f64) {
    // Independently integrate only the actual positive-traction disk, using
    // theta coordinates and Simpson pressure quadrature, then fixed-grid RK4.
    let h = dt / steps as f64;
    let acceleration = |q: f64| {
        let upper = (0.4 / q).min(1.0).asin();
        let dtheta = upper / radial as f64;
        let mut sum = 0.0;
        for i in 0..=radial {
            let (s, c) = (i as f64 * dtheta).sin_cos();
            let weight = if i == 0 || i == radial {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            sum += weight * (0.4 - q * s).max(0.0) * s * s * c * c;
        }
        -diagonal()[2][2] * N * dtheta * sum
    };
    let mut q = initial;
    let mut distance = 0.0;
    for _ in 0..steps {
        let a = acceleration(q);
        let qb = q + 0.5 * h * a;
        let b = acceleration(qb);
        let qc = q + 0.5 * h * b;
        let c = acceleration(qc);
        let qd = q + h * c;
        let d = acceleration(qd);
        distance += h * (q + 2.0 * qb + 2.0 * qc + qd) / 6.0;
        q += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    (q, distance)
}

#[test]
fn partially_clipped_spin_uses_original_active_pressure_region_through_time() {
    let (initial, dt) = (0.7, 0.3);
    let coarse = clipped_spin_reference(initial, dt, 256, 256);
    let fine = clipped_spin_reference(initial, dt, 512, 512);
    close(coarse.0, fine.0, 2e-10);
    close(coarse.1, fine.1, 2e-10);
    let patch = HertzFrictionPatch::new(
        N,
        A,
        FrictionLaw::VelocityDependent {
            static_mu: 0.8,
            mu_zero: 0.4,
            slope_per_speed: -1.0,
        },
    )
    .unwrap();
    let r = patch
        .advance_transient(
            &scalar_input(2, initial, dt),
            config(),
            grid(16),
            controls(),
            || false,
        )
        .unwrap();
    scalar_check(&r, 2, initial, fine.0, fine.1);
    assert_eq!(r.kinetic_coefficient_bounds[0], 0.0);
    assert!(r.trajectory.stopping_time_bracket_s.is_none());
    assert!(
        q(&r)[2] > 0.4,
        "fixture must retain its clipped outer region"
    );
}
