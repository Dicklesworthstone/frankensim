//! G1/G3: joint Hertz stiction/sliding, full body recoil and original local work.

use core::f64::consts::PI;

use fs_tribo::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchError, HertzWrenchInput, HertzWrenchMode,
    HertzWrenchParameters,
};
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};

type Vector = [f64; 3];
type Matrix = [[f64; 3]; 3];

fn law() -> HertzWrenchContact {
    HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: 10.0,
        contact_radius_m: 0.25,
        static_mu: 0.6,
        kinetic_mu: 0.4,
    })
    .unwrap()
}

fn config() -> HertzWrenchConfig {
    HertzWrenchConfig {
        maximum_static_iterations: 96,
        maximum_iterations: 256,
        maximum_evaluations: 512,
        force_tolerance_n: 2e-9,
        moment_tolerance_nm: 5e-10,
        displacement_tolerance_m: 2e-13,
        energy_tolerance_j: 2e-10,
        relative_matrix_tolerance: 1e-13,
        maximum_force_n: 20.0,
        maximum_spin_moment_nm: 5.0,
    }
}

fn matrix() -> Matrix {
    [
        [0.02, 0.005, -0.003],
        [0.005, 0.03, 0.002],
        [-0.003, 0.002, 0.04],
    ]
}

fn dot(a: Vector, b: Vector) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}

fn mul(a: Matrix, b: Vector) -> Vector {
    a.map(|row| dot(row, b))
}

fn add(a: Vector, b: Vector) -> Vector {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn norm(a: Vector) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}

fn input(free: Vector, compliance: Matrix, radius: f64) -> HertzWrenchInput {
    HertzWrenchInput {
        free_displacement_m: [free[0], free[1]],
        free_spin_angle_rad: free[2] / radius,
        scaled_compliance_m_per_n: compliance,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance {tolerance:e}"
    );
}

fn vector_close(a: Vector, b: Vector, tolerance: f64) {
    for (a, b) in a.into_iter().zip(b) {
        close(a, b, tolerance);
    }
}

#[test]
fn full_static_reaction_and_exact_pure_translation_spin_limits() {
    let law = law();
    let radius = law.parameters().contact_radius_m;
    let requested = [1.2, -0.7, 0.3];
    let held = law
        .advance(
            &input(mul(matrix(), requested), matrix(), radius),
            config(),
            || false,
        )
        .unwrap();
    assert_eq!(held.mode, HertzWrenchMode::Sticking);
    vector_close(
        [
            held.force_n[0],
            held.force_n[1],
            held.spin_moment_nm / radius,
        ],
        requested.map(|v| -v),
        1e-15,
    );
    assert_eq!(held.slip_displacement_m, [0.0; 2]);
    assert_eq!(held.spin_angle_rad, 0.0);
    assert_eq!(held.contact_work_j, 0.0);
    assert_eq!(held.dissipated_energy_j, 0.0);

    let diagonal = [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.3]];
    let translate = law
        .advance(&input([2.0, 0.0, 0.0], diagonal, radius), config(), || {
            false
        })
        .unwrap();
    assert_eq!(translate.mode, HertzWrenchMode::Sliding);
    assert_eq!(translate.force_n, [-4.0, 0.0]);
    assert_eq!(translate.spin_moment_nm, 0.0);
    close(translate.slip_displacement_m[0], 1.6, 2e-16);
    close(translate.dissipated_energy_j, 6.4, 1e-15);
    assert_eq!(translate.iterations, 0);

    let spin = law
        .advance(&input([0.0, 0.0, 3.0], diagonal, radius), config(), || {
            false
        })
        .unwrap();
    let spin_capacity = 3.0 * PI * 4.0 * radius / 16.0;
    assert_eq!(spin.force_n, [0.0; 2]);
    close(spin.spin_moment_nm, -spin_capacity, 2e-16);
    close(
        spin.spin_angle_rad,
        (3.0 - 0.3 * spin_capacity / radius) / radius,
        2e-15,
    );
    close(
        spin.dissipated_energy_j,
        spin_capacity * spin.spin_angle_rad,
        2e-15,
    );
    assert_eq!(spin.iterations, 0);
}

#[test]
fn static_priority_and_joint_surface_do_not_use_independent_force_torque_caps() {
    let law = law();
    let radius = law.parameters().contact_radius_m;
    let c = [[0.1, 0.0, 0.0], [0.0, 0.15, 0.0], [0.0, 0.0, 0.2]];
    // This reaction exceeds the kinetic limit but is a valid static one.
    let static_input = input(mul(c, [4.8, 0.0, 0.0]), c, radius);
    let held = law.advance(&static_input, config(), || false).unwrap();
    assert_eq!(held.mode, HertzWrenchMode::Sticking);
    close(held.force_n[0], -4.8, 1e-15);
    let mut limited = config();
    limited.maximum_force_n = 4.5;
    assert!(matches!(
        law.advance(&static_input, limited, || false),
        Err(HertzWrenchError::EffortBudget { .. })
    ));

    // Each component passes its pure static cap, but their pair is infeasible.
    let requested = [0.9 * 6.0, 0.0, 0.8 * (3.0 * PI / 16.0) * 6.0];
    let step = law
        .advance(&input(mul(c, requested), c, radius), config(), || false)
        .unwrap();
    assert_eq!(step.mode, HertzWrenchMode::Sliding);
    assert!(step.slip_displacement_m[0] > 0.0 && step.spin_angle_rad > 0.0);
    assert!(step.force_n[0].abs() < 4.0);
    assert!(step.spin_moment_nm.abs() < (3.0 * PI / 16.0) * 4.0 * radius);
    assert!(step.dissipated_energy_j > 0.0);
    close(step.work_residual_j, 0.0, config().energy_tolerance_j);
}

#[test]
fn anisotropic_translation_matches_independent_coulomb_kkt() {
    let law = law();
    let c = [[0.015, 0.006, 0.0], [0.006, 0.03, 0.0], [0.0, 0.0, 0.02]];
    let free = [0.3, -0.25, 0.0];
    let reaction = |lambda: f64| {
        let a = c[0][0] + lambda;
        let b = c[1][1] + lambda;
        let determinant = a * b - c[0][1] * c[0][1];
        [
            (b * free[0] - c[0][1] * free[1]) / determinant,
            (a * free[1] - c[0][1] * free[0]) / determinant,
            0.0,
        ]
    };
    let mut low = 0.0;
    let mut high = 1.0;
    for _ in 0..100 {
        let middle = (low + high) * 0.5;
        if norm(reaction(middle)) > 4.0 {
            low = middle;
        } else {
            high = middle;
        }
    }
    let expected = reaction((low + high) * 0.5);
    let step = law
        .advance(&input(free, c, 0.25), config(), || false)
        .unwrap();
    vector_close(
        [step.force_n[0], step.force_n[1], step.spin_moment_nm],
        expected.map(|v| -v),
        3e-9,
    );
    assert_eq!(step.spin_angle_rad, 0.0);
    assert!(step.force_residual_n <= config().force_tolerance_n);
    assert!(step.iterations > 0);
}

#[test]
fn exact_static_extrema_exclude_every_nonzero_cross_component() {
    let law = law();
    let c = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for required in [[6.0, 0.0, 1e-9], [1e-9, 0.0, (3.0 * PI / 16.0) * 6.0]] {
        let step = law
            .advance(&input(required, c, 0.25), config(), || false)
            .unwrap();
        assert_eq!(step.mode, HertzWrenchMode::Sliding);
        assert!(step.dissipated_energy_j > 0.0);
    }
}

// Independent spatial midpoint integration of local traction and local work.
// The production solver uses the analytic original kernel and no pressure grid.
fn spatial(law: HertzWrenchContact, motion: Vector, radial: usize, angular: usize) -> [f64; 4] {
    let p = law.parameters();
    let mut result = [0.0; 4];
    let dr = (PI / 2.0) / radial as f64;
    for i in 0..radial {
        let theta = (i as f64 + 0.5) * dr;
        let (sine, cosine) = theta.sin_cos();
        let r = p.contact_radius_m * sine;
        let weight =
            3.0 * p.normal_force_n * p.kinetic_mu * sine * cosine * cosine * dr / angular as f64;
        for j in 0..angular {
            let angle = 2.0 * PI * (j as f64 + 0.5) / angular as f64;
            let (s, c) = angle.sin_cos();
            let x = r * c;
            let y = r * s;
            let spin = motion[2] / p.contact_radius_m;
            let velocity = [motion[0] - spin * y, motion[1] + spin * x];
            let speed = velocity[0].hypot(velocity[1]);
            if speed > 0.0 {
                let traction = velocity.map(|v| weight * v / speed);
                result[0] += traction[0];
                result[1] += traction[1];
                result[2] += (x * traction[1] - y * traction[0]) / p.contact_radius_m;
                result[3] += weight * speed;
            }
        }
    }
    result
}

#[test]
fn full_coupled_recoil_and_local_heat_converge_to_independent_surface_quadrature() {
    let law = law();
    let actual = [0.2, -0.15, 0.07];
    let mut errors = Vec::new();
    for radial in [32, 64, 128] {
        let reference = spatial(law, actual, radial, 256);
        let wrench = [reference[0], reference[1], reference[2]];
        let free = add(actual, mul(matrix(), wrench));
        let step = law
            .advance(&input(free, matrix(), 0.25), config(), || false)
            .unwrap();
        let observed = [
            -step.force_n[0],
            -step.force_n[1],
            -step.spin_moment_nm / 0.25,
        ];
        let error = norm([
            observed[0] - wrench[0],
            observed[1] - wrench[1],
            observed[2] - wrench[2],
        ]);
        errors.push(error);
        assert!(step.iterations > 0);
        assert_eq!(step.compatibility_residual_m, [0.0; 3]);
        close(
            step.contact_work_j,
            -(step.force_n[0] * step.slip_displacement_m[0]
                + step.force_n[1] * step.slip_displacement_m[1]
                + step.spin_moment_nm * step.spin_angle_rad),
            2e-15,
        );
        close(step.work_residual_j, 0.0, config().energy_tolerance_j);
        if radial == 128 {
            vector_close(
                [
                    step.slip_displacement_m[0],
                    step.slip_displacement_m[1],
                    0.25 * step.spin_angle_rad,
                ],
                actual,
                3e-6,
            );
            close(step.dissipated_energy_j, reference[3], 4e-5);
        }
    }
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            order > 1.9 && order < 2.1,
            "errors {errors:?}, order {order}"
        );
    }
    assert!(errors[2] < 1e-4, "errors {errors:?}");
}

fn rotate_matrix(q: Matrix, matrix: Matrix) -> Matrix {
    let mut result = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            result[i][j] = dot(q[i], mul(matrix, q[j]));
        }
    }
    result
}

#[test]
fn coupled_solution_is_covariant_under_tangent_rotations_and_reflections() {
    let law = law();
    let free = [0.34, -0.2, 0.15];
    let base = law
        .advance(&input(free, matrix(), 0.25), config(), || false)
        .unwrap();
    let (s, c) = 0.713_f64.sin_cos();
    for q in [
        [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]],
        [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
    ] {
        let step = law
            .advance(
                &input(mul(q, free), rotate_matrix(q, matrix()), 0.25),
                config(),
                || false,
            )
            .unwrap();
        vector_close(
            [step.force_n[0], step.force_n[1], step.spin_moment_nm / 0.25],
            mul(
                q,
                [base.force_n[0], base.force_n[1], base.spin_moment_nm / 0.25],
            ),
            8e-9,
        );
        vector_close(
            [
                step.slip_displacement_m[0],
                step.slip_displacement_m[1],
                0.25 * step.spin_angle_rad,
            ],
            mul(
                q,
                [
                    base.slip_displacement_m[0],
                    base.slip_displacement_m[1],
                    0.25 * base.spin_angle_rad,
                ],
            ),
            3e-10,
        );
        close(step.dissipated_energy_j, base.dissipated_energy_j, 5e-10);
    }
}

#[test]
fn dimensional_similarity_preserves_force_torque_and_work_scalings() {
    let law = law();
    let free = [0.34, -0.2, 0.15];
    let base = law
        .advance(&input(free, matrix(), 0.25), config(), || false)
        .unwrap();
    for (length, force) in [(1000.0, 0.01), (0.001, 100.0)] {
        let mut parameters = law.parameters();
        parameters.contact_radius_m *= length;
        parameters.normal_force_n *= force;
        let scaled_law = HertzWrenchContact::new(parameters).unwrap();
        let scaled_c = matrix().map(|row| row.map(|v| v * length / force));
        let mut cfg = config();
        cfg.force_tolerance_n *= force;
        cfg.moment_tolerance_nm *= force * length;
        cfg.displacement_tolerance_m *= length;
        cfg.energy_tolerance_j *= force * length;
        cfg.maximum_force_n *= force;
        cfg.maximum_spin_moment_nm *= force * length;
        let scaled = scaled_law
            .advance(
                &input(
                    free.map(|v| v * length),
                    scaled_c,
                    parameters.contact_radius_m,
                ),
                cfg,
                || false,
            )
            .unwrap();
        close(scaled.force_n[0] / force, base.force_n[0], 8e-9);
        close(scaled.force_n[1] / force, base.force_n[1], 8e-9);
        close(
            scaled.spin_moment_nm / (force * length),
            base.spin_moment_nm,
            2e-9,
        );
        close(
            scaled.dissipated_energy_j / (force * length),
            base.dissipated_energy_j,
            5e-10,
        );
    }
}

#[test]
fn zero_body_response_rest_and_zero_capacities_keep_native_selection() {
    let law = law();
    let free = [0.34, -0.2, 0.15];
    let trial = input(free, [[0.0; 3]; 3], 0.25);
    let native = HertzSlidingPatch::new(HertzSlidingPatchParameters {
        normal_force_n: 10.0,
        contact_radius_m: 0.25,
        friction_coefficient: 0.4,
    })
    .unwrap()
    .evaluate(trial.free_displacement_m, trial.free_spin_angle_rad)
    .unwrap();
    let step = law.advance(&trial, config(), || false).unwrap();
    assert_eq!(step.slip_displacement_m, trial.free_displacement_m);
    assert_eq!(step.spin_angle_rad, trial.free_spin_angle_rad);
    close(step.force_n[0], native.force_n[0], 1e-15);
    close(step.force_n[1], native.force_n[1], 1e-15);
    close(step.spin_moment_nm, native.spin_moment_nm, 1e-16);
    close(step.dissipated_energy_j, native.dissipated_power_w, 3e-16);
    assert_eq!(step.evaluations, 1);
    let rest = law
        .advance(&input([0.0; 3], matrix(), 0.25), config(), || false)
        .unwrap();
    assert_eq!(rest.mode, HertzWrenchMode::Sticking);
    assert_eq!(rest.force_n, [0.0; 2]);
    for zero_load in [true, false] {
        let mut parameters = law.parameters();
        if zero_load {
            parameters.normal_force_n = 0.0;
        } else {
            parameters.static_mu = 0.0;
            parameters.kinetic_mu = 0.0;
        }
        let result = HertzWrenchContact::new(parameters)
            .unwrap()
            .advance(&input(free, matrix(), 0.25), config(), || false)
            .unwrap();
        assert_eq!(result.mode, HertzWrenchMode::Frictionless);
        assert_eq!(result.force_n, [0.0; 2]);
        assert_eq!(result.spin_moment_nm, 0.0);
        assert_eq!(result.dissipated_energy_j, 0.0);
    }
}

#[test]
fn refusal_cancellation_and_retry_are_atomic_without_surface_regularization() {
    let law = law();
    let trial = input([0.34, -0.2, 0.15], matrix(), 0.25);
    let saved_law = law;
    let saved_trial = trial;
    let clean = law.advance(&trial, config(), || false).unwrap();
    let mut calls = 0;
    assert_eq!(
        law.advance(&trial, config(), || {
            calls += 1;
            calls == 4
        }),
        Err(HertzWrenchError::Cancelled)
    );
    let mut short = config();
    short.maximum_iterations = 1;
    assert!(matches!(
        law.advance(&trial, short, || false),
        Err(HertzWrenchError::Unresolved { .. })
    ));
    short = config();
    short.maximum_evaluations = 1;
    assert!(law.advance(&trial, short, || false).is_err());
    // Static classification ambiguity is an explicit refusal, not kinetic onset.
    let required = [0.9 * 6.0, 0.0, 0.35 * (3.0 * PI / 16.0) * 6.0];
    short = config();
    short.maximum_static_iterations = 1;
    assert!(matches!(
        law.advance(
            &input(mul(matrix(), required), matrix(), 0.25),
            short,
            || false
        ),
        Err(HertzWrenchError::StaticBoundaryUnresolved { .. })
    ));
    for bad in [
        [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.0]],
        [[0.1, 0.2, 0.0], [0.2, 0.1, 0.0], [0.0, 0.0, 0.1]],
        [[0.1, 0.01, 0.0], [0.02, 0.1, 0.0], [0.0, 0.0, 0.1]],
        [[f64::NAN, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]],
    ] {
        assert!(matches!(
            law.advance(&input([0.34, -0.2, 0.15], bad, 0.25), config(), || false),
            Err(HertzWrenchError::Invalid(_))
        ));
    }
    let mut limited = config();
    limited.maximum_spin_moment_nm = 0.2;
    assert!(matches!(
        law.advance(
            &input(mul(matrix(), [0.0, 0.0, 1.0]), matrix(), 0.25),
            limited,
            || false
        ),
        Err(HertzWrenchError::EffortBudget { .. })
    ));
    assert_eq!(law, saved_law);
    assert_eq!(trial, saved_trial);
    assert_eq!(law.advance(&trial, config(), || false).unwrap(), clean);
    let tiny = HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: 1e-30,
        contact_radius_m: 1e-300,
        static_mu: 1.0,
        kinetic_mu: 1.0,
    })
    .unwrap();
    assert!(matches!(
        tiny.advance(
            &HertzWrenchInput {
                free_displacement_m: [0.0; 2],
                free_spin_angle_rad: 1e300,
                scaled_compliance_m_per_n: [[0.0; 3]; 3],
            },
            config(),
            || false
        ),
        Err(HertzWrenchError::Invalid("unrepresentable_spin_port"))
    ));
}
