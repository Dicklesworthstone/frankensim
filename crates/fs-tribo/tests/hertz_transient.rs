//! G1/G3: original Hertz traction, independently integrated curved motion,
//! analytical anisotropic-mobility arrest and bounded uncomputed travel.
use core::f64::consts::PI;
use fs_tribo::sliding_patch::contact::relaxation::HertzRelaxationInput;
use fs_tribo::sliding_patch::contact::transient::{
    HertzTransientControls, HertzTransientError, HertzTransientReport,
};
use fs_tribo::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchContact, HertzWrenchParameters,
};
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};

const RADIUS: f64 = 0.1;
type V = [f64; 3];
type M = [[f64; 3]; 3];
fn law() -> HertzWrenchContact {
    HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: 2.0,
        contact_radius_m: RADIUS,
        static_mu: 4.0,
        kinetic_mu: 0.5,
    })
    .unwrap()
}
fn patch() -> HertzSlidingPatch {
    HertzSlidingPatch::new(HertzSlidingPatchParameters {
        normal_force_n: 2.0,
        contact_radius_m: RADIUS,
        friction_coefficient: 0.5,
    })
    .unwrap()
}
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
        maximum_force_n: 2.0,
        maximum_spin_moment_nm: 1.0,
    }
}
fn controls() -> HertzTransientControls {
    HertzTransientControls {
        maximum_step_s: 0.1,
        minimum_step_s: 1e-16,
        maximum_attempts: 100_000,
        maximum_evaluations: 400_000,
        scaled_rate_tolerance_m_s: 2e-8,
        scaled_impulse_tolerance_ns: 2e-8,
        scaled_travel_tolerance_m: 2e-8,
        heat_tolerance_j: 2e-9,
        stopping_time_tolerance_s: 2e-7,
    }
}
fn matrix() -> M {
    [[2.0, 0.3, 0.4], [0.3, 1.7, -0.2], [0.4, -0.2, 3.5]]
}
fn input(dt_s: f64) -> HertzRelaxationInput {
    HertzRelaxationInput {
        velocity_m_s: [0.8, -0.3],
        spin_rate_rad_s: 5.0,
        scaled_mobility_per_kg: matrix(),
        dt_s,
    }
}
fn add(a: V, b: V) -> V {
    core::array::from_fn(|i| a[i] + b[i])
}
fn scale(a: V, s: f64) -> V {
    a.map(|v| s * v)
}
fn mul(m: M, v: V) -> V {
    m.map(|r| r[0] * v[0] + r[1] * v[1] + r[2] * v[2])
}
fn norm(v: V) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tol {tolerance}"
    );
}
fn rates(r: &HertzTransientReport) -> V {
    [
        r.velocity_after_m_s[0],
        r.velocity_after_m_s[1],
        RADIUS * r.spin_rate_after_rad_s,
    ]
}
fn impulses(r: &HertzTransientReport) -> V {
    [
        r.force_impulse_ns[0],
        r.force_impulse_ns[1],
        r.spin_angular_impulse_nms / RADIUS,
    ]
}
fn travel(r: &HertzTransientReport) -> V {
    [
        r.slip_increment_m[0],
        r.slip_increment_m[1],
        RADIUS * r.spin_increment_rad,
    ]
}

// Classical RK4 at an independently fixed physical grid. No production
// adaptive stages, matrix metric, extinction bounds or step decisions are used.
fn reference(input: HertzRelaxationInput, count: usize) -> (V, V, V, f64) {
    reference_for_patch(input, count, patch())
}
fn reference_for_patch(
    input: HertzRelaxationInput,
    count: usize,
    law: HertzSlidingPatch,
) -> (V, V, V, f64) {
    let radius = law.parameters().contact_radius_m;
    let h = input.dt_s / count as f64;
    let mut q = [
        input.velocity_m_s[0],
        input.velocity_m_s[1],
        radius * input.spin_rate_rad_s,
    ];
    let (mut impulse, mut travel, mut heat) = ([0.0; 3], [0.0; 3], 0.0);
    let sample = |q: V| {
        let value = law.evaluate([q[0], q[1]], q[2] / radius).unwrap();
        (
            [
                value.force_n[0],
                value.force_n[1],
                value.spin_moment_nm / radius,
            ],
            value.dissipated_power_w,
        )
    };
    for _ in 0..count {
        let (a, pa) = sample(q);
        let qb = add(q, scale(mul(input.scaled_mobility_per_kg, a), 0.5 * h));
        let (b, pb) = sample(qb);
        let qc = add(q, scale(mul(input.scaled_mobility_per_kg, b), 0.5 * h));
        let (c, pc) = sample(qc);
        let qd = add(q, scale(mul(input.scaled_mobility_per_kg, c), h));
        let (d, pd) = sample(qd);
        let j = core::array::from_fn(|i| h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
        travel = add(
            travel,
            core::array::from_fn(|i| h * (q[i] + 2.0 * qb[i] + 2.0 * qc[i] + qd[i]) / 6.0),
        );
        q = add(q, mul(input.scaled_mobility_per_kg, j));
        impulse = add(impulse, j);
        heat += h * (pa + 2.0 * pb + 2.0 * pc + pd) / 6.0;
    }
    (q, impulse, travel, heat)
}

#[test]
fn small_patch_strict_thermal_accuracy_finishes_within_the_original_work_budget() {
    // The original finite-body thermal fixture has fast shear and much slower
    // torsional recoil. Its strict controls previously consumed 100,000 BS3(2)
    // attempts before reaching even one eighth of this physical interval.
    let radius = (0.001 / PI).sqrt();
    let original = HertzSlidingPatch::new(HertzSlidingPatchParameters {
        normal_force_n: 10.0,
        contact_radius_m: radius,
        friction_coefficient: 0.4,
    })
    .unwrap();
    let law = HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: 10.0,
        contact_radius_m: radius,
        static_mu: 0.8,
        kinetic_mu: 0.4,
    })
    .unwrap();
    let input = HertzRelaxationInput {
        velocity_m_s: [0.1, 0.0],
        spin_rate_rad_s: 3.0,
        scaled_mobility_per_kg: [
            [
                1.0 / 2.0 + 1.0 / 3.0 + 0.2 * 0.2 / 0.10 + 0.15 * 0.15 / 0.20,
                0.0,
                0.0,
            ],
            [
                0.0,
                1.0 / 2.0 + 1.0 / 3.0 + 0.2 * 0.2 / 0.08 + 0.15 * 0.15 / 0.15,
                0.0,
            ],
            [0.0, 0.0, radius * radius * (1.0 / 0.12 + 1.0 / 0.25)],
        ],
        dt_s: 0.2,
    };
    let controls = HertzTransientControls {
        maximum_step_s: 0.02,
        minimum_step_s: 1e-15,
        maximum_attempts: 100_000,
        maximum_evaluations: 400_001,
        scaled_rate_tolerance_m_s: 1e-11,
        scaled_impulse_tolerance_ns: 1e-10,
        scaled_travel_tolerance_m: 1e-10,
        heat_tolerance_j: 1e-12,
        stopping_time_tolerance_s: 1e-7,
    };
    let config = HertzWrenchConfig {
        maximum_force_n: 20.0,
        maximum_spin_moment_nm: 2.0,
        energy_tolerance_j: 2e-11,
        ..config()
    };
    let reference = reference_for_patch(input, 16_384, original);
    let coarse_reference = reference_for_patch(input, 8_192, original);
    let result = law
        .advance_transient(&input, config, controls, || false)
        .unwrap();
    let rate = [
        result.velocity_after_m_s[0],
        result.velocity_after_m_s[1],
        radius * result.spin_rate_after_rad_s,
    ];
    let impulse = [
        result.force_impulse_ns[0],
        result.force_impulse_ns[1],
        result.spin_angular_impulse_nms / radius,
    ];
    let travel = [
        result.slip_increment_m[0],
        result.slip_increment_m[1],
        radius * result.spin_increment_rad,
    ];
    for i in 0..3 {
        close(reference.0[i], coarse_reference.0[i], 2e-13);
        close(rate[i], reference.0[i], controls.scaled_rate_tolerance_m_s);
        close(
            impulse[i],
            reference.1[i],
            controls.scaled_impulse_tolerance_ns,
        );
        close(
            travel[i],
            reference.2[i],
            controls.scaled_travel_tolerance_m,
        );
    }
    close(
        result.dissipated_energy_j,
        reference.3,
        controls.heat_tolerance_j,
    );
    assert_eq!(result.integrated_until_s, input.dt_s);
    assert!(result.stopping_time_bracket_s.is_none());
    assert!(result.accepted_steps + result.rejected_steps < controls.maximum_attempts);
    assert!(result.evaluations <= controls.maximum_evaluations);
    assert!(result.work_residual_j.abs() <= config.energy_tolerance_j);
}
fn inverse(m: M, b: V) -> V {
    let [a, d, g] = m[0];
    let [_, e, h] = m[1];
    let [_, _, i] = m[2];
    let determinant = a * (e * i - h * h) - d * (d * i - h * g) + g * (d * h - e * g);
    mul(
        [
            [e * i - h * h, g * h - d * i, d * h - e * g],
            [g * h - d * i, a * i - g * g, d * g - a * h],
            [d * h - e * g, d * g - a * h, a * e - d * d],
        ],
        b,
    )
    .map(|v| v / determinant)
}

#[test]
fn full_shear_spin_cross_response_matches_independent_physical_time_integration() {
    let input = input(0.2);
    let reference = reference(input, 16_384);
    let result = law()
        .advance_transient(&input, config(), controls(), || false)
        .unwrap();
    assert!(!result.exact_ray && result.accepted_steps > 1);
    assert!(result.stopping_time_bracket_s.is_none());
    for i in 0..3 {
        close(rates(&result)[i], reference.0[i], 2e-8);
        close(impulses(&result)[i], reference.1[i], 2e-8);
        close(travel(&result)[i], reference.2[i], 2e-8);
    }
    close(result.dissipated_energy_j, reference.3, 2e-9);
    assert!(result.work_residual_j.abs() <= config().energy_tolerance_j);
    assert!(result.prefix_work_residual_j.abs() <= config().energy_tolerance_j);
    assert!((rates(&result)[0] * 0.5 - rates(&result)[2] * 0.8).abs() > 1e-3);
}

#[test]
fn anisotropic_recoil_has_analytical_curved_stopping_time_and_travel() {
    // With C=1 and W=diag(1,2,3), q=(z,z²,0) and
    // dt=-sqrt(1+z²) dz. Its entire changing-direction path is integrable.
    let input = HertzRelaxationInput {
        velocity_m_s: [1.0, 1.0],
        spin_rate_rad_s: 0.0,
        scaled_mobility_per_kg: [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]],
        dt_s: 5.0,
    };
    let primitive = |z: f64| 0.5 * (z * (1.0 + z * z).sqrt() + z.asinh());
    let x_primitive = |z: f64| (1.0 + z * z).powf(1.5) / 3.0;
    let y_primitive = |z: f64| (z * (2.0 * z * z + 1.0) * (1.0 + z * z).sqrt() - z.asinh()) / 8.0;
    let expected_stop = primitive(1.0);
    let complete = law()
        .advance_transient(&input, config(), controls(), || false)
        .unwrap();
    let [lower, upper] = complete.stopping_time_bracket_s.unwrap();
    assert!(!complete.exact_ray && lower < upper);
    assert!(upper - lower <= controls().stopping_time_tolerance_s);
    // The prefix has estimated error; this reference tests the whole numerical
    // trajectory without claiming an interval-certified prefix.
    assert!((0.5 * (lower + upper) - expected_stop).abs() < controls().stopping_time_tolerance_s);
    close(
        complete.slip_increment_m[0],
        x_primitive(1.0) - x_primitive(0.0),
        3e-8,
    );
    close(
        complete.slip_increment_m[1],
        y_primitive(1.0) - y_primitive(0.0),
        3e-8,
    );
    close(
        complete.dissipated_energy_j,
        0.75,
        config().energy_tolerance_j,
    );
    assert_eq!(rates(&complete), [0.0; 3]);
    close(complete.force_impulse_ns[0], -1.0, 2e-14);
    close(complete.force_impulse_ns[1], -0.5, 2e-14);
    let z = 0.5;
    let partial = law()
        .advance_transient(
            &HertzRelaxationInput {
                dt_s: primitive(1.0) - primitive(z),
                ..input
            },
            config(),
            controls(),
            || false,
        )
        .unwrap();
    assert!(partial.stopping_time_bracket_s.is_none());
    close(partial.velocity_after_m_s[0], z, 2e-8);
    close(partial.velocity_after_m_s[1], z * z, 2e-8);
    close(
        partial.slip_increment_m[0],
        x_primitive(1.0) - x_primitive(z),
        2e-8,
    );
    close(
        partial.slip_increment_m[1],
        y_primitive(1.0) - y_primitive(z),
        2e-8,
    );
}

#[test]
fn full_matrix_extinction_preserves_physical_work_and_additional_dwell_is_free() {
    let input = input(20.0);
    let result = law()
        .advance_transient(&input, config(), controls(), || false)
        .unwrap();
    let longer = law()
        .advance_transient(
            &HertzRelaxationInput {
                dt_s: 200.0,
                ..input
            },
            config(),
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(result, longer);
    let [lower, upper] = result.stopping_time_bracket_s.unwrap();
    assert!(result.integrated_until_s < lower && lower < upper && upper < input.dt_s);
    assert!(upper - lower <= controls().stopping_time_tolerance_s);
    assert_eq!(rates(&result), [0.0; 3]);
    let q = [0.8, -0.3, 0.5];
    let expected = inverse(matrix(), q).map(|v| -v);
    for i in 0..3 {
        close(impulses(&result)[i], expected[i], 2e-14);
    }
    let energy = -0.5 * (q[0] * expected[0] + q[1] * expected[1] + q[2] * expected[2]);
    close(
        result.dissipated_energy_j,
        energy,
        config().energy_tolerance_j,
    );
    assert!(result.stopping_tail_heat_j > 0.0);
    assert_eq!(
        result.dissipated_energy_j,
        result.prefix_dissipated_energy_j + result.stopping_tail_heat_j
    );
    assert!(result.prefix_work_residual_j.abs() <= config().energy_tolerance_j);
    assert!(
        norm([
            result.tail_slip_bound_m[0],
            result.tail_slip_bound_m[1],
            RADIUS * result.tail_spin_bound_rad
        ]) + result.estimated_scaled_travel_error_m
            <= controls().scaled_travel_tolerance_m
    );
}

#[test]
fn refinement_and_split_continuation_converge_to_the_same_original_motion() {
    let input = input(0.2);
    let truth = reference(input, 16_384);
    let coarse_controls = HertzTransientControls {
        maximum_step_s: 0.3,
        scaled_rate_tolerance_m_s: 2e-5,
        scaled_impulse_tolerance_ns: 2e-5,
        scaled_travel_tolerance_m: 2e-5,
        heat_tolerance_j: 2e-6,
        ..controls()
    };
    let coarse = law()
        .advance_transient(
            &input,
            HertzWrenchConfig {
                energy_tolerance_j: 2e-6,
                ..config()
            },
            coarse_controls,
            || false,
        )
        .unwrap();
    let fine = law()
        .advance_transient(&input, config(), controls(), || false)
        .unwrap();
    let error = |r: &HertzTransientReport| norm(add(rates(r), scale(truth.0, -1.0)));
    assert!(error(&fine) < error(&coarse));
    assert!(fine.accepted_steps > coarse.accepted_steps);
    let first = law()
        .advance_transient(
            &HertzRelaxationInput {
                dt_s: 0.08,
                ..input
            },
            config(),
            controls(),
            || false,
        )
        .unwrap();
    let second = law()
        .advance_transient(
            &HertzRelaxationInput {
                velocity_m_s: first.velocity_after_m_s,
                spin_rate_rad_s: first.spin_rate_after_rad_s,
                dt_s: 0.12,
                ..input
            },
            config(),
            controls(),
            || false,
        )
        .unwrap();
    for i in 0..3 {
        close(rates(&second)[i], rates(&fine)[i], 4e-8);
        close(
            impulses(&first)[i] + impulses(&second)[i],
            impulses(&fine)[i],
            4e-8,
        );
        close(
            travel(&first)[i] + travel(&second)[i],
            travel(&fine)[i],
            4e-8,
        );
    }
    close(
        first.dissipated_energy_j + second.dissipated_energy_j,
        fine.dissipated_energy_j,
        4e-9,
    );
}

#[test]
fn tangent_rotation_and_velocity_reversal_preserve_the_physical_trajectory() {
    let input = input(0.2);
    let original = law()
        .advance_transient(&input, config(), controls(), || false)
        .unwrap();
    let rotation = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    let transformed = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            (0..3)
                .map(|a| {
                    (0..3)
                        .map(|b| {
                            rotation[i][a] * input.scaled_mobility_per_kg[a][b] * rotation[j][b]
                        })
                        .sum::<f64>()
                })
                .sum::<f64>()
        })
    });
    let rotated = law()
        .advance_transient(
            &HertzRelaxationInput {
                velocity_m_s: [0.3, 0.8],
                scaled_mobility_per_kg: transformed,
                ..input
            },
            config(),
            controls(),
            || false,
        )
        .unwrap();
    let reversed = law()
        .advance_transient(
            &HertzRelaxationInput {
                velocity_m_s: [-0.8, 0.3],
                spin_rate_rad_s: -5.0,
                ..input
            },
            config(),
            controls(),
            || false,
        )
        .unwrap();
    for extractor in [rates, impulses, travel] {
        let expected = mul(rotation, extractor(&original));
        for i in 0..3 {
            close(extractor(&rotated)[i], expected[i], 3e-8);
            close(extractor(&reversed)[i], -extractor(&original)[i], 2e-14);
        }
    }
    close(
        rotated.dissipated_energy_j,
        original.dissipated_energy_j,
        3e-9,
    );
    assert_eq!(reversed.dissipated_energy_j, original.dissipated_energy_j);
}

#[test]
fn hertz_pressure_moment_supplies_the_actual_extinction_norm_bounds() {
    let k = 3.0 * PI / 16.0;
    let upper = (1.0 + k * k).sqrt();
    let lower = k / upper;
    for speed in [0.0, 1e-10, 0.1, 1.0, 100.0] {
        for spin in [-1000.0, -10.0, -1e-10, 0.0, 1e-10, 10.0, 1000.0] {
            let q = [0.6 * speed, 0.8 * speed, RADIUS * spin];
            let result = patch().evaluate([q[0], q[1]], spin).unwrap();
            assert!(result.dissipated_power_w >= lower * norm(q) * (1.0 - 8.0 * f64::EPSILON));
            assert!(result.dissipated_power_w <= upper * norm(q) * (1.0 + 8.0 * f64::EPSILON));
        }
    }
}

#[test]
fn pure_modes_rest_and_frictionless_motion_keep_the_exact_fast_path() {
    for (velocity, spin, normal) in [
        ([0.4, 0.0], 0.0, 2.0),
        ([0.0, 0.0], 4.0, 2.0),
        ([0.0, 0.0], 0.0, 2.0),
        ([0.4, 0.2], 4.0, 0.0),
    ] {
        let input = HertzRelaxationInput {
            velocity_m_s: velocity,
            spin_rate_rad_s: spin,
            scaled_mobility_per_kg: [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]],
            dt_s: 100.0,
        };
        let law = HertzWrenchContact::new(HertzWrenchParameters {
            normal_force_n: normal,
            ..law().parameters()
        })
        .unwrap();
        let exact = law.relax_velocity(&input, config(), || false).unwrap();
        let transient = law
            .advance_transient(&input, config(), controls(), || false)
            .unwrap();
        assert!(transient.exact_ray);
        assert_eq!(transient.evaluations, 1);
        assert_eq!(transient.force_impulse_ns, exact.force_impulse_ns);
        assert_eq!(
            transient.spin_angular_impulse_nms,
            exact.spin_angular_impulse_nms
        );
        assert_eq!(transient.dissipated_energy_j, exact.dissipated_energy_j);
        assert_eq!(
            transient.stopping_time_bracket_s,
            exact.kinetic_phase_end_s.map(|t| [t, t])
        );
        assert_eq!(transient.tail_slip_bound_m, [0.0; 2]);
    }
}

#[test]
fn budget_and_cancellation_refusals_leave_the_original_law_and_input_replayable() {
    let input = input(0.2);
    let law = law();
    let before = (input, law);
    assert!(matches!(
        law.advance_transient(
            &input,
            config(),
            HertzTransientControls {
                maximum_evaluations: 1,
                ..controls()
            },
            || false
        ),
        Err(HertzTransientError::EvaluationBudget)
    ));
    assert!(matches!(
        law.advance_transient(
            &input,
            config(),
            HertzTransientControls {
                maximum_attempts: 1,
                ..controls()
            },
            || false
        ),
        Err(HertzTransientError::AttemptBudget)
    ));
    assert!(matches!(
        law.advance_transient(
            &input,
            config(),
            HertzTransientControls {
                maximum_step_s: 1.0,
                minimum_step_s: 0.5,
                ..controls()
            },
            || false
        ),
        Err(HertzTransientError::StepUnderflow)
    ));
    let mut polls = 0;
    assert!(matches!(
        law.advance_transient(&input, config(), controls(), || {
            polls += 1;
            polls == 20
        }),
        Err(HertzTransientError::Cancelled)
    ));
    assert_eq!(polls, 20);
    assert_eq!((input, law), before);
    let retry = law
        .advance_transient(&input, config(), controls(), || false)
        .unwrap();
    assert_eq!(
        retry,
        law.advance_transient(&input, config(), controls(), || false)
            .unwrap()
    );
    let different_static = HertzWrenchContact::new(HertzWrenchParameters {
        static_mu: 100.0,
        ..law.parameters()
    })
    .unwrap();
    assert_eq!(
        retry,
        different_static
            .advance_transient(&input, config(), controls(), || false)
            .unwrap()
    );
}
