use fs_contact::tangential::impact::multiple::{
    CoupledFrictionImpulseControls, CoupledFrictionImpulseError, CoupledFrictionImpulseInput,
    CoupledFrictionImpulseReport, MAX_FRICTION_IMPULSE_PORTS, solve_coupled_friction_impulses,
};
use fs_contact::tangential::impact::{
    CoulombImpulseInput, FrictionLaw, solve_isotropic_coulomb_impulse,
};

fn controls() -> CoupledFrictionImpulseControls {
    CoupledFrictionImpulseControls {
        maximum_iterations: 256,
        impulse_tolerance_ns: 2e-11,
        velocity_tolerance_m_s: 2e-11,
        energy_tolerance_j: 2e-10,
        relative_matrix_tolerance: 1e-12,
        maximum_impulse_ns: 100.0,
    }
}

fn solve(velocity: &[f64], matrix: &[f64], capacities: &[f64]) -> CoupledFrictionImpulseReport {
    solve_coupled_friction_impulses(
        &CoupledFrictionImpulseInput {
            velocity_before_m_s: velocity,
            mobility_per_kg: matrix,
            capacity_ns: capacities,
        },
        controls(),
        || false,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_scalar_limits_match_original_coulomb_without_a_small_velocity_dead_band() {
    for velocity in [-7.0, -0.02, 0.0, 1e-12, 0.5, 4.0] {
        for capacity in [0.0, 1e-14, 0.2, 20.0] {
            for mobility in [0.1, 1.0, 3.0] {
                let original = solve_isotropic_coulomb_impulse(
                    CoulombImpulseInput {
                        normal_impulse: 1.0,
                        relative_tangent_velocity: [velocity, 0.0],
                        inverse_tangent_mass: mobility,
                    },
                    &FrictionLaw::Coulomb {
                        static_mu: capacity,
                        kinetic_mu: capacity,
                    },
                )
                .unwrap();
                let result = solve(&[velocity], &[mobility], &[capacity]);
                near(
                    result.impulse_ns[0],
                    original.tangent_impulse[0],
                    2e-13 * original.tangent_impulse[0].abs().max(1e-300),
                );
                near(
                    result.velocity_after_m_s[0],
                    original.relative_tangent_velocity_after[0],
                    2e-13 * velocity.abs().max(1e-300),
                );
                near(
                    result.dissipated_j,
                    original.dissipated_energy,
                    2e-13 * original.dissipated_energy.max(1e-300),
                );
                if capacity > 0.0 && velocity != 0.0 {
                    assert_ne!(result.impulse_ns[0], 0.0);
                }
            }
        }
    }
}

#[test]
fn g1_three_finite_rotors_require_the_coupled_mixed_active_set() {
    // Two adjacent relative-spin ports on rotors with inertias 2,3,5.
    // Initial spins [4,0,3] give relative rates [-4,3]. The exact solution
    // saturates only the first joint; its impulse changes the second's stop.
    let w = [5.0 / 6.0, -1.0 / 3.0, -1.0 / 3.0, 8.0 / 15.0];
    let result = solve(&[-4.0, 3.0], &w, &[1.0, 10.0]);
    near(result.impulse_ns[0], 1.0, 1e-14);
    near(result.impulse_ns[1], -5.0, 3e-14);
    near(result.velocity_after_m_s[0], -1.5, 3e-14);
    near(result.velocity_after_m_s[1], 0.0, 3e-14);
    near(result.dissipated_j, 10.25, 5e-14);
    let p = result.impulse_ns;
    let after = [4.0 - p[0] / 2.0, (p[0] - p[1]) / 3.0, 3.0 + p[1] / 5.0];
    near(
        2.0 * after[0] + 3.0 * after[1] + 5.0 * after[2],
        23.0,
        3e-14,
    );
    let kinetic = after
        .iter()
        .zip([2.0, 3.0, 5.0])
        .map(|(w, i)| 0.5 * i * w * w)
        .sum::<f64>();
    near(38.5 - kinetic, result.dissipated_j, 5e-14);
    let reversed = solve(&[3.0, -4.0], &[w[3], w[2], w[1], w[0]], &[10.0, 1.0]);
    near(reversed.impulse_ns[1], p[0], 3e-14);
    near(reversed.impulse_ns[0], p[1], 3e-14);
    near(reversed.dissipated_j, result.dissipated_j, 3e-14);
}

#[test]
fn g1_bound_release_and_signed_port_work_preserve_the_global_event() {
    // The initial free solve first hits the WRONG lower face for port one.
    // Once the second port saturates, the first must be released into the box.
    let released = solve(&[-1.0, -2.0], &[1.0, 0.8, 0.8, 1.0], &[0.5, 1.0]);
    near(released.impulse_ns[0], 0.2, 3e-14);
    near(released.impulse_ns[1], 1.0, 3e-14);
    near(released.velocity_after_m_s[0], 0.0, 3e-14);
    near(released.velocity_after_m_s[1], -0.84, 3e-14);

    let result = solve(&[-1.0, -10.0], &[2.0, 1.0, 1.0, 2.0], &[10.0, 2.0]);
    near(result.impulse_ns[0], -0.5, 3e-14);
    near(result.impulse_ns[1], 2.0, 3e-14);
    near(result.velocity_after_m_s[0], 0.0, 3e-14);
    near(result.velocity_after_m_s[1], -6.5, 3e-14);
    let first_work = 0.5 * result.impulse_ns[0] * (-1.0 + result.velocity_after_m_s[0]);
    let second_work = 0.5 * result.impulse_ns[1] * (-10.0 + result.velocity_after_m_s[1]);
    near(first_work, 0.25, 2e-14);
    near(second_work, -16.5, 3e-14);
    near(result.dissipated_j, 16.25, 5e-14);
    near(result.kinetic_energy_change_j, -16.25, 5e-14);
}

#[test]
fn g1_manufactured_full_tensor_ports_cover_every_box_face() {
    let w = [2.0, 0.5, -0.25, 0.5, 1.5, 0.2, -0.25, 0.2, 1.25];
    let caps = [1.0, 2.0, 3.0];
    for pattern in 0..27 {
        let mut code = pattern;
        let mut expected = [0.0; 3];
        let mut outgoing = [0.0; 3];
        for i in 0..3 {
            let face = code % 3;
            code /= 3;
            (expected[i], outgoing[i]) = match face {
                0 => (-caps[i], 0.7 + i as f64),
                1 => (0.25 * (i + 1) as f64, 0.0),
                _ => (caps[i], -0.7 - i as f64),
            };
        }
        // The supplied target satisfies the original box KKT exactly by
        // construction; strict diagonal dominance makes the minimizer unique.
        let before: [f64; 3] = core::array::from_fn(|i| {
            outgoing[i] - (0..3).map(|j| w[3 * i + j] * expected[j]).sum::<f64>()
        });
        let result = solve(&before, &w, &caps);
        for i in 0..3 {
            near(result.impulse_ns[i], expected[i], 3e-13);
            near(result.velocity_after_m_s[i], outgoing[i], 5e-13);
        }
        let objective = |p: [f64; 3]| {
            (0..3)
                .map(|i| {
                    p[i] * (before[i] + 0.5 * (0..3).map(|j| w[3 * i + j] * p[j]).sum::<f64>())
                })
                .sum::<f64>()
        };
        for sample in 0..27 {
            let mut index = sample;
            let candidate = core::array::from_fn(|i| {
                let coordinate = (index % 3) as f64 - 1.0;
                index /= 3;
                caps[i] * coordinate
            });
            assert!(objective(candidate) >= objective(expected) - 3e-13);
        }
        near(result.dissipated_j, -objective(expected), 3e-13);
    }
}

#[test]
fn g3_conjugate_coordinate_changes_reversal_and_units_preserve_physics() {
    let w = [5.0 / 6.0, -1.0 / 3.0, -1.0 / 3.0, 8.0 / 15.0];
    let before = [-4.0, 3.0];
    let caps = [1.0, 10.0];
    let baseline = solve(&before, &w, &caps);
    for scales in [[1e-4, 1e4], [3.0, 0.25], [-2.0, 0.5]] {
        let matrix: [f64; 4] = core::array::from_fn(|i| scales[i / 2] * w[i] * scales[i % 2]);
        let velocity: [f64; 2] = core::array::from_fn(|i| scales[i] * before[i]);
        let capacities: [f64; 2] = core::array::from_fn(|i| caps[i] / scales[i].abs());
        let mut c = controls();
        c.maximum_impulse_ns = 2e5;
        c.impulse_tolerance_ns = 1e-8;
        c.velocity_tolerance_m_s = 1e-8;
        let result = solve_coupled_friction_impulses(
            &CoupledFrictionImpulseInput {
                velocity_before_m_s: &velocity,
                mobility_per_kg: &matrix,
                capacity_ns: &capacities,
            },
            c,
            || false,
        )
        .unwrap();
        for i in 0..2 {
            near(
                scales[i] * result.impulse_ns[i],
                baseline.impulse_ns[i],
                2e-13,
            );
            near(
                result.velocity_after_m_s[i] / scales[i],
                baseline.velocity_after_m_s[i],
                3e-13,
            );
        }
        near(result.dissipated_j, baseline.dissipated_j, 1e-12);
    }
    for (velocity_unit, impulse_unit) in [(1e-20, 1e20), (1e20, 1e-20), (3.0, 7.0)] {
        let matrix = w.map(|x| x * velocity_unit / impulse_unit);
        let velocity = before.map(|x| x * velocity_unit);
        let capacities = caps.map(|x| x * impulse_unit);
        let mut c = controls();
        c.maximum_impulse_ns *= impulse_unit;
        c.impulse_tolerance_ns *= impulse_unit;
        c.velocity_tolerance_m_s *= velocity_unit;
        c.energy_tolerance_j *= velocity_unit * impulse_unit;
        let result = solve_coupled_friction_impulses(
            &CoupledFrictionImpulseInput {
                velocity_before_m_s: &velocity,
                mobility_per_kg: &matrix,
                capacity_ns: &capacities,
            },
            c,
            || false,
        )
        .unwrap();
        for i in 0..2 {
            near(
                result.impulse_ns[i] / impulse_unit,
                baseline.impulse_ns[i],
                3e-13,
            );
            near(
                result.velocity_after_m_s[i] / velocity_unit,
                baseline.velocity_after_m_s[i],
                3e-13,
            );
        }
        near(
            result.dissipated_j / (velocity_unit * impulse_unit),
            baseline.dissipated_j,
            5e-13,
        );
    }
}

#[test]
fn g0_empty_zero_capacity_and_maximum_dimension_retain_exact_limits() {
    let empty = solve(&[], &[], &[]);
    assert_eq!(empty.count, 0);
    assert_eq!(empty.impulse_ns, [0.0; MAX_FRICTION_IMPULSE_PORTS]);
    assert_eq!(empty.dissipated_j.to_bits(), 0.0_f64.to_bits());
    assert_eq!(empty.iterations, 0);
    assert_eq!(empty.factorizations, 0);
    let unchanged = solve(&[1.25, -3.5], &[2.0, 0.25, 0.25, 1.0], &[0.0, 0.0]);
    assert_eq!(&unchanged.velocity_after_m_s[..2], &[1.25, -3.5]);
    assert_eq!(unchanged.impulse_ns, [0.0; MAX_FRICTION_IMPULSE_PORTS]);
    assert_eq!(unchanged.dissipated_j, 0.0);

    let n = MAX_FRICTION_IMPULSE_PORTS;
    let mut w = [0.0; MAX_FRICTION_IMPULSE_PORTS * MAX_FRICTION_IMPULSE_PORTS];
    let mut expected = [0.0_f64; MAX_FRICTION_IMPULSE_PORTS];
    let mut caps = [2.0; MAX_FRICTION_IMPULSE_PORTS];
    let mut outgoing = [0.0; MAX_FRICTION_IMPULSE_PORTS];
    for i in 0..n {
        w[n * i + i] = 2.0;
        if i > 0 {
            w[n * i + i - 1] = -0.25;
            w[n * (i - 1) + i] = -0.25;
        }
        expected[i] = if i % 2 == 0 { 0.25 } else { -0.5 };
        if i % 7 == 0 {
            caps[i] = 0.0;
            expected[i] = 0.0;
            outgoing[i] = 1.75;
        } else if i % 3 == 0 {
            caps[i] = expected[i].abs();
            outgoing[i] = -expected[i].signum();
        }
    }
    let before: [f64; MAX_FRICTION_IMPULSE_PORTS] = core::array::from_fn(|i| {
        outgoing[i] - (0..n).map(|j| w[n * i + j] * expected[j]).sum::<f64>()
    });
    let result = solve(&before, &w, &caps);
    assert_eq!(result.count, n);
    for i in 0..n {
        near(result.impulse_ns[i], expected[i], 1e-13);
        near(result.velocity_after_m_s[i], outgoing[i], 1e-13);
    }
    assert!(
        solve_coupled_friction_impulses(
            &CoupledFrictionImpulseInput {
                velocity_before_m_s: &[0.0; 33],
                mobility_per_kg: &[0.0; 33 * 33],
                capacity_ns: &[0.0; 33],
            },
            controls(),
            || false
        )
        .is_err()
    );
}

#[test]
fn g4_physical_capacity_is_not_a_numerical_budget_and_original_mobility_is_audited() {
    let velocity = [-4.0, 3.0];
    let w = [5.0 / 6.0, -1.0 / 3.0, -1.0 / 3.0, 8.0 / 15.0];
    let caps = [1.0, 10.0];
    let input = CoupledFrictionImpulseInput {
        velocity_before_m_s: &velocity,
        mobility_per_kg: &w,
        capacity_ns: &caps,
    };
    let mut effort = controls();
    effort.maximum_impulse_ns = 4.0;
    assert!(matches!(
        solve_coupled_friction_impulses(&input, effort, || false),
        Err(CoupledFrictionImpulseError::ImpulseBudget { index: 1, .. })
    ));
    let mut short = controls();
    short.maximum_iterations = 1;
    assert!(matches!(
        solve_coupled_friction_impulses(&input, short, || false),
        Err(CoupledFrictionImpulseError::IterationLimit { iterations: 1, .. })
    ));
    for invalid in [
        [1.0, 0.1, 0.0, 1.0],
        [1.0, 1.0, 1.0, 1.0],
        [1.0, 2.0, 2.0, 1.0],
        [1.0, f64::NAN, 0.0, 1.0],
    ] {
        assert!(
            solve_coupled_friction_impulses(
                &CoupledFrictionImpulseInput {
                    mobility_per_kg: &invalid,
                    ..input
                },
                controls(),
                || false
            )
            .is_err()
        );
    }
    let mut tight = controls();
    tight.impulse_tolerance_ns = 1e-16;
    assert!(matches!(
        solve_coupled_friction_impulses(
            &CoupledFrictionImpulseInput {
                velocity_before_m_s: &[-1.0, 0.0],
                mobility_per_kg: &[1.0, 1e-13, 0.0, 1.0],
                capacity_ns: &[2.0, 2.0],
            },
            tight,
            || false
        ),
        Err(CoupledFrictionImpulseError::Unresolved { .. })
    ));
    assert_eq!(velocity, [-4.0, 3.0]);
    assert_eq!(caps, [1.0, 10.0]);
}

#[test]
fn g4_every_observed_cancellation_checkpoint_returns_no_candidate() {
    let velocity = [-1.0, -2.0];
    let w = [1.0, 0.8, 0.8, 1.0];
    let caps = [0.5, 1.0];
    let input = CoupledFrictionImpulseInput {
        velocity_before_m_s: &velocity,
        mobility_per_kg: &w,
        capacity_ns: &caps,
    };
    let mut calls = 0;
    let baseline = solve_coupled_friction_impulses(&input, controls(), || {
        calls += 1;
        false
    })
    .unwrap();
    assert!(calls > 10);
    for stop in 1..=calls {
        let mut seen = 0;
        let result = solve_coupled_friction_impulses(&input, controls(), || {
            seen += 1;
            seen == stop
        });
        assert!(
            matches!(result, Err(CoupledFrictionImpulseError::Cancelled)),
            "checkpoint {stop}: {result:?}"
        );
        assert_eq!(velocity, [-1.0, -2.0]);
        assert_eq!(w, [1.0, 0.8, 0.8, 1.0]);
        assert_eq!(caps, [0.5, 1.0]);
    }
    assert_eq!(solve(&velocity, &w, &caps), baseline);
}
