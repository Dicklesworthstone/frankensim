//! Curved trajectories use actual uniform-disk traction, with separate numerical error.
use fs_tribo::adhesive::spinning::contact::{
    AdhesiveWrenchConfig, AdhesiveWrenchContact, AdhesiveWrenchParameters,
};
use fs_tribo::adhesive::spinning::relaxation::AdhesiveRelaxationInput;
use fs_tribo::adhesive::spinning::transient::{AdhesiveTransientControls, AdhesiveTransientError};

const RADIUS: f64 = 0.1;
type V = [f64; 3];
type M = [[f64; 3]; 3];
fn law() -> AdhesiveWrenchContact {
    AdhesiveWrenchContact::new(AdhesiveWrenchParameters {
        contact_radius_m: RADIUS,
        static_shear_strength_pa: 30.0,
        kinetic_shear_strength_pa: 4.0,
    })
    .unwrap()
}
fn config() -> AdhesiveWrenchConfig {
    AdhesiveWrenchConfig {
        maximum_static_iterations: 64,
        maximum_iterations: 128,
        maximum_evaluations: 256,
        force_tolerance_n: 1e-11,
        moment_tolerance_nm: 1e-12,
        displacement_tolerance_m: 1e-12,
        energy_tolerance_j: 2e-9,
        relative_matrix_tolerance: 1e-13,
        maximum_force_n: 1.0,
        maximum_spin_moment_nm: 1.0,
    }
}
fn controls() -> AdhesiveTransientControls {
    AdhesiveTransientControls {
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
fn input(dt_s: f64) -> AdhesiveRelaxationInput {
    AdhesiveRelaxationInput {
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

// Independent fixed-step classical RK4 reference. It uses the already-verified
// original disk integral, but none of the production adaptive/time/tail code.
fn reference(input: AdhesiveRelaxationInput, count: usize) -> (V, V, V, f64) {
    let law = law();
    let h = input.dt_s / count as f64;
    let mut q = [
        input.velocity_m_s[0],
        input.velocity_m_s[1],
        RADIUS * input.spin_rate_rad_s,
    ];
    let mut impulse = [0.0; 3];
    let mut travel = [0.0; 3];
    let mut heat = 0.0;
    let sample = |v: V| {
        let native = law
            .kinetic_patch()
            .evaluate([v[0], v[1]], v[2] / RADIUS)
            .unwrap();
        (
            [
                native.force_n[0],
                native.force_n[1],
                native.spin_moment_nm / RADIUS,
            ],
            native.dissipated_power_w,
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
fn full_offdiagonal_curved_motion_matches_an_independent_rk4_reference() {
    let input = input(0.8);
    let reference = reference(input, 16_384);
    let result = law()
        .advance_transient(&input, config(), controls(), || false)
        .unwrap();
    assert!(!result.exact_ray && result.accepted_steps > 1);
    assert!(result.stopping_time_bracket_s.is_none());
    let q = [
        result.velocity_after_m_s[0],
        result.velocity_after_m_s[1],
        RADIUS * result.spin_rate_after_rad_s,
    ];
    let j = [
        result.force_impulse_ns[0],
        result.force_impulse_ns[1],
        result.spin_angular_impulse_nms / RADIUS,
    ];
    let x = [
        result.slip_increment_m[0],
        result.slip_increment_m[1],
        RADIUS * result.spin_increment_rad,
    ];
    for i in 0..3 {
        close(q[i], reference.0[i], 2e-8);
        close(j[i], reference.1[i], 2e-8);
        close(x[i], reference.2[i], 2e-8);
    }
    close(result.dissipated_energy_j, reference.3, 2e-9);
    assert!(result.work_residual_j.abs() <= config().energy_tolerance_j);
    assert!(result.prefix_work_residual_j.abs() <= config().energy_tolerance_j);
    // Motion actually changes direction; a fixed-ray substitution cannot pass.
    assert!((q[0] * 0.5 - q[2] * 0.8).abs() > 1e-4);
}

#[test]
fn conservative_extinction_tail_closes_actual_impulse_and_energy_without_invented_travel() {
    let input = input(100.0);
    let result = law()
        .advance_transient(&input, config(), controls(), || false)
        .unwrap();
    let [lower, upper] = result.stopping_time_bracket_s.unwrap();
    assert!(!result.exact_ray && lower < upper && upper < input.dt_s);
    assert!(result.integrated_until_s < lower);
    assert!(upper - lower <= controls().stopping_time_tolerance_s);
    assert_eq!(result.velocity_after_m_s, [0.0; 2]);
    assert_eq!(result.spin_rate_after_rad_s, 0.0);
    let q = [0.8, -0.3, 0.5];
    let expected = inverse(matrix(), q).map(|v| -v);
    close(result.force_impulse_ns[0], expected[0], 2e-14);
    close(result.force_impulse_ns[1], expected[1], 2e-14);
    close(result.spin_angular_impulse_nms, RADIUS * expected[2], 2e-15);
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
    assert!(
        norm([
            result.tail_slip_bound_m[0],
            result.tail_slip_bound_m[1],
            RADIUS * result.tail_spin_bound_rad
        ]) + result.estimated_scaled_travel_error_m
            <= controls().scaled_travel_tolerance_m
    );
    // An earlier requested horizon remains moving; a small rate is not a stop.
    let before = law()
        .advance_transient(
            &AdhesiveRelaxationInput {
                dt_s: 0.8 * lower,
                ..input
            },
            config(),
            controls(),
            || false,
        )
        .unwrap();
    assert!(before.stopping_time_bracket_s.is_none());
    assert!(
        norm([
            before.velocity_after_m_s[0],
            before.velocity_after_m_s[1],
            RADIUS * before.spin_rate_after_rad_s
        ]) > 0.0
    );
}

#[test]
fn physical_time_and_accuracy_refinement_converge_without_changing_the_law() {
    let input = input(0.8);
    let truth = reference(input, 16_384);
    let coarse_control = AdhesiveTransientControls {
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
            AdhesiveWrenchConfig {
                energy_tolerance_j: 2e-6,
                ..config()
            },
            coarse_control,
            || false,
        )
        .unwrap();
    let fine = law()
        .advance_transient(&input, config(), controls(), || false)
        .unwrap();
    let error = |r: &fs_tribo::adhesive::spinning::transient::AdhesiveTransientReport| {
        norm([
            r.velocity_after_m_s[0] - truth.0[0],
            r.velocity_after_m_s[1] - truth.0[1],
            RADIUS * r.spin_rate_after_rad_s - truth.0[2],
        ])
    };
    assert!(error(&fine) < error(&coarse));
    assert!(fine.accepted_steps > coarse.accepted_steps);
    assert!(fine.estimated_scaled_rate_error_m_s < coarse.estimated_scaled_rate_error_m_s);
}

#[test]
fn pure_shear_and_spin_keep_the_exact_event_fast_path() {
    let law = law();
    for spin in [false, true] {
        let input = AdhesiveRelaxationInput {
            velocity_m_s: if spin { [0.0; 2] } else { [0.4, 0.0] },
            spin_rate_rad_s: if spin { 4.0 } else { 0.0 },
            scaled_mobility_per_kg: [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]],
            dt_s: 100.0,
        };
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
    // Original machine-roundoff ray admission must not bypass a stricter
    // transient endpoint-rate gate when the full response is only nearly a ray.
    let almost_ray = AdhesiveRelaxationInput {
        velocity_m_s: [1.0, 1.0],
        spin_rate_rad_s: 0.0,
        scaled_mobility_per_kg: [[1.0, 0.0, 0.0], [0.0, 1.0 + 1e-14, 0.0], [0.0, 0.0, 1.0]],
        dt_s: 100.0,
    };
    assert!(law.relax_velocity(&almost_ray, config(), || false).is_ok());
    assert_eq!(
        law.advance_transient(
            &almost_ray,
            config(),
            AdhesiveTransientControls {
                scaled_rate_tolerance_m_s: 1e-18,
                ..controls()
            },
            || false
        ),
        Err(AdhesiveTransientError::UnresolvedWork)
    );
}

#[test]
fn tangent_gauge_rotation_preserves_the_curved_physics() {
    let input = input(0.8);
    let original = law()
        .advance_transient(&input, config(), controls(), || false)
        .unwrap();
    let indices = [1, 0, 2];
    let signs = [-1.0, 1.0, 1.0];
    let rotated_input = AdhesiveRelaxationInput {
        velocity_m_s: [-input.velocity_m_s[1], input.velocity_m_s[0]],
        scaled_mobility_per_kg: core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                signs[i] * signs[j] * input.scaled_mobility_per_kg[indices[i]][indices[j]]
            })
        }),
        ..input
    };
    let rotated = law()
        .advance_transient(&rotated_input, config(), controls(), || false)
        .unwrap();
    close(
        rotated.force_impulse_ns[0],
        -original.force_impulse_ns[1],
        2e-9,
    );
    close(
        rotated.force_impulse_ns[1],
        original.force_impulse_ns[0],
        2e-9,
    );
    close(
        rotated.spin_angular_impulse_nms,
        original.spin_angular_impulse_nms,
        2e-10,
    );
    close(
        rotated.dissipated_energy_j,
        original.dissipated_energy_j,
        2e-9,
    );
}

#[test]
fn bounded_work_and_cancellation_refuse_atomically_and_retry_deterministically() {
    let law = law();
    let input = input(0.8);
    let expected = law
        .advance_transient(&input, config(), controls(), || false)
        .unwrap();
    for at in [0, 6, 20] {
        let mut calls = 0;
        assert_eq!(
            law.advance_transient(&input, config(), controls(), || {
                let stop = calls == at;
                calls += 1;
                stop
            }),
            Err(AdhesiveTransientError::Cancelled)
        );
        assert_eq!(
            law.advance_transient(&input, config(), controls(), || false)
                .unwrap(),
            expected
        );
    }
    assert_eq!(
        law.advance_transient(
            &input,
            config(),
            AdhesiveTransientControls {
                maximum_evaluations: 2,
                ..controls()
            },
            || false
        ),
        Err(AdhesiveTransientError::EvaluationBudget)
    );
    assert_eq!(
        law.advance_transient(
            &input,
            config(),
            AdhesiveTransientControls {
                maximum_attempts: 1,
                ..controls()
            },
            || false
        ),
        Err(AdhesiveTransientError::AttemptBudget)
    );
    let mut asymmetric = input;
    asymmetric.scaled_mobility_per_kg[0][1] += 1e-14;
    assert!(
        law.advance_transient(&asymmetric, config(), controls(), || false)
            .is_err()
    );
}
