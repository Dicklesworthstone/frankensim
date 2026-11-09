//! Original physical damping/Maxwell equations and transient calibration.
use fs_phs::linear_exact::{DenseExponentialBudget, LinearSplitBudget, PreparedExactZoh};
use fs_phs::mechanical::{LinearMechanicalSystem, MechanicalDampingDirection};
use fs_phs::{PreparedStepError, RelaxationBranch};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}
fn dense() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 64 * 1024 * 1024,
        max_multiply_terms: 1_000_000_000,
    }
}
fn split() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 32,
        max_omega_dt: 100.0,
    }
}
const MASS: [f64; 4] = [2.0, 0.4, 0.4, 1.5];
const DAMPING: [f64; 4] = [0.8, 1.55, -1.05, 0.6];
const DIRECTION: [f64; 4] = [0.3, -0.2, 0.4, 0.5];
const STIFFNESS: [f64; 4] = [7.0, -1.0, -1.0, 4.0];
const PORTS: [f64; 4] = [1.0, 0.2, -0.3, 0.8];
const PROJECTIONS: [[f64; 2]; 3] = [[1.0, -0.6], [0.3, 1.0], [0.0, 0.0]];
const BRANCH_K: [f64; 3] = [3.0, 2.0, 4.0];
const RATES: [f64; 3] = [12.5, 1.0 / 0.6, 1.0 / 0.6];
const RATE_DIRECTION: [f64; 3] = [1.2, -0.3, 0.8];

fn model(theta: f64) -> LinearMechanicalSystem {
    let damping = core::array::from_fn::<_, 4, _>(|i| DAMPING[i] + theta * DIRECTION[i]);
    let branches = (0..3)
        .map(|i| RelaxationBranch {
            projection: PROJECTIONS[i].to_vec(),
            stiffness: BRANCH_K[i],
            relaxation_time_s: 1.0 / (RATES[i] + theta * RATE_DIRECTION[i]),
        })
        .collect();
    LinearMechanicalSystem::new(&MASS, &damping, &STIFFNESS, &PORTS, 2, 2, branches).unwrap()
}
fn direction() -> MechanicalDampingDirection<'static> {
    MechanicalDampingDirection {
        damping: &DIRECTION,
        relaxation_rates: &RATE_DIRECTION,
    }
}
fn initial(source: &LinearMechanicalSystem) -> Vec<f64> {
    source
        .state(&[0.2, -0.15], &[0.1, 0.25], &[0.05, -0.1, 0.35])
        .unwrap()
}
fn native(source: &LinearMechanicalSystem, time: f64, state: &[f64], force: &[f64]) -> Vec<f64> {
    PreparedExactZoh::new(source.system(), time, split())
        .unwrap()
        .step(state, force)
        .unwrap()
        .x
}

#[test]
fn dimensional_free_and_critical_limits_have_analytic_damping_derivatives() {
    let time: f64 = 0.7;
    let mass: f64 = 2.0;
    let source =
        LinearMechanicalSystem::new(&[mass], &[0.0], &[0.0], &[1.0], 1, 1, vec![]).unwrap();
    let cache = source
        .prepare_damping_sensitivity(
            MechanicalDampingDirection {
                damping: &[1.0],
                relaxation_rates: &[],
            },
            time,
            dense(),
            || false,
        )
        .unwrap();
    let x = source.state(&[0.4], &[-0.3], &[]).unwrap();
    let mut tangent = [0.0; 2];
    cache
        .advance_into(&x, &[0.0; 2], &[3.0], &[0.0], &mut tangent)
        .unwrap();
    close(
        tangent[0],
        0.3 * time.powi(2) / (2.0 * mass) - 3.0 * time.powi(3) / (6.0 * mass * mass),
        3e-15,
    );
    close(
        tangent[1] / mass,
        0.3 * time / mass - 3.0 * time.powi(2) / (2.0 * mass * mass),
        3e-15,
    );
    let omega: f64 = 2.0;
    for offset in [-1e-9, 0.0, 1e-9] {
        let source = LinearMechanicalSystem::new(
            &[mass],
            &[2.0 * mass * omega + offset],
            &[mass * omega * omega],
            &[],
            1,
            0,
            vec![],
        )
        .unwrap();
        let cache = source
            .prepare_damping_sensitivity(
                MechanicalDampingDirection {
                    damping: &[1.0],
                    relaxation_rates: &[],
                },
                time,
                dense(),
                || false,
            )
            .unwrap();
        cache
            .advance_into(
                &source.state(&[1.0], &[0.0], &[]).unwrap(),
                &[0.0; 2],
                &[],
                &[],
                &mut tangent,
            )
            .unwrap();
        let scale = omega * omega * (-omega * time).exp() / (6.0 * mass);
        close(tangent[0], scale * time.powi(3), 2e-11);
        close(
            tangent[1] / mass,
            scale * (3.0 * time.powi(2) - omega * time.powi(3)),
            2e-11,
        );
    }
}

// Independent physical [q0,q1,v0,v1,s0,s1,s2] plus its directional derivative.
// This never builds a pHS generator, factor derivative or matrix exponential.
fn physical_rhs(x: &[f64; 14], force: [f64; 2], force_direction: [f64; 2]) -> [f64; 14] {
    let mut result = [0.0; 14];
    for tangent in 0..2 {
        let offset = 7 * tangent;
        result[offset] = x[offset + 2];
        result[offset + 1] = x[offset + 3];
        let u = if tangent == 0 { force } else { force_direction };
        let mut rhs = [
            PORTS[0] * u[0] + PORTS[1] * u[1],
            PORTS[2] * u[0] + PORTS[3] * u[1],
        ];
        for row in 0..2 {
            for column in 0..2 {
                rhs[row] -= STIFFNESS[2 * row + column] * x[offset + column]
                    + DAMPING[2 * row + column] * x[offset + 2 + column];
                if tangent == 1 {
                    rhs[row] -= DIRECTION[2 * row + column] * x[2 + column];
                }
            }
            for branch in 0..3 {
                rhs[row] -= PROJECTIONS[branch][row] * x[offset + 4 + branch];
            }
        }
        result[offset + 2] = (1.5 * rhs[0] - 0.4 * rhs[1]) / 2.84;
        result[offset + 3] = (-0.4 * rhs[0] + 2.0 * rhs[1]) / 2.84;
        for branch in 0..3 {
            result[offset + 4 + branch] = BRANCH_K[branch]
                * (PROJECTIONS[branch][0] * x[offset + 2] + PROJECTIONS[branch][1] * x[offset + 3])
                - RATES[branch] * x[offset + 4 + branch];
            if tangent == 1 {
                result[offset + 4 + branch] -= RATE_DIRECTION[branch] * x[4 + branch];
            }
        }
    }
    result
}

#[test]
fn full_consistent_mass_skew_damping_and_original_memory_match_variational_ode() {
    let source = model(0.0);
    let x = initial(&source);
    let x_theta = source
        .state(&[0.03, -0.02], &[-0.04, 0.06], &[0.02, -0.01, 0.07])
        .unwrap();
    let force = [0.7, -0.4];
    let force_theta = [-0.2, 0.3];
    let time = 0.8;
    let cache = source
        .prepare_damping_sensitivity(direction(), time, dense(), || false)
        .unwrap();
    assert_eq!(
        (cache.state_dim(), cache.port_dim(), cache.duration()),
        (7, 2, time)
    );
    let mut actual = [0.0; 7];
    cache
        .advance_into(&x, &x_theta, &force, &force_theta, &mut actual)
        .unwrap();
    let stress = source.branch_stresses(&x).unwrap();
    let stress_theta = source.branch_stresses(&x_theta).unwrap();
    let mut physical = [
        0.2,
        -0.15,
        0.1,
        0.25,
        stress[0],
        stress[1],
        stress[2],
        0.03,
        -0.02,
        -0.04,
        0.06,
        stress_theta[0],
        stress_theta[1],
        stress_theta[2],
    ];
    let steps = 8000;
    let h = time / steps as f64;
    for _ in 0..steps {
        let a = physical_rhs(&physical, force, force_theta);
        let b = physical_rhs(
            &core::array::from_fn(|i| physical[i] + 0.5 * h * a[i]),
            force,
            force_theta,
        );
        let c = physical_rhs(
            &core::array::from_fn(|i| physical[i] + 0.5 * h * b[i]),
            force,
            force_theta,
        );
        let d = physical_rhs(
            &core::array::from_fn(|i| physical[i] + h * c[i]),
            force,
            force_theta,
        );
        for i in 0..14 {
            physical[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    let velocity = source.velocities(&actual).unwrap();
    let stresses = source.branch_stresses(&actual).unwrap();
    let observed = [
        actual[0],
        actual[1],
        velocity[0],
        velocity[1],
        stresses[0],
        stresses[1],
        stresses[2],
    ];
    for i in 0..7 {
        close(observed[i], physical[7 + i], 3e-12);
    }
    // Third arm has zero projection but an independently decaying native history.
    let rho = RATES[2];
    close(
        actual[6],
        (-rho * time).exp() * (x_theta[6] - time * RATE_DIRECTION[2] * x[6]),
        2e-14,
    );
}

#[test]
fn exact_direction_composes_through_native_substeps_and_zero_direction() {
    let source = model(0.0);
    let state = initial(&source);
    let input = [0.7, -0.4];
    let whole = source
        .prepare_damping_sensitivity(direction(), 0.8, dense(), || false)
        .unwrap();
    let small = source
        .prepare_damping_sensitivity(direction(), 0.05, dense(), || false)
        .unwrap();
    let flow = PreparedExactZoh::new(source.system(), 0.05, split()).unwrap();
    let mut expected = [0.0; 7];
    whole
        .advance_into(&state, &[0.0; 7], &input, &[0.0; 2], &mut expected)
        .unwrap();
    let mut x = state.clone();
    let mut tangent = [0.0; 7];
    for _ in 0..16 {
        let mut next = [0.0; 7];
        small
            .advance_into(&x, &tangent, &input, &[0.0; 2], &mut next)
            .unwrap();
        x = flow.step(&x, &input).unwrap().x;
        tangent = next;
    }
    for i in 0..7 {
        close(tangent[i], expected[i], 8e-13);
    }
    // A local direction is signed, not itself a physical resistance. Parameter
    // unit changes alter the derivative units without changing source physics.
    for parameter_units in [-1.0, 1e-6, 1e6] {
        let damping = DIRECTION.map(|value| parameter_units * value);
        let rates = RATE_DIRECTION.map(|value| parameter_units * value);
        let rescaled = source
            .prepare_damping_sensitivity(
                MechanicalDampingDirection {
                    damping: &damping,
                    relaxation_rates: &rates,
                },
                0.8,
                dense(),
                || false,
            )
            .unwrap();
        rescaled
            .advance_into(&state, &[0.0; 7], &input, &[0.0; 2], &mut tangent)
            .unwrap();
        for i in 0..7 {
            close(tangent[i] / parameter_units, expected[i], 3e-13);
        }
    }
    let zero = source
        .prepare_damping_sensitivity(
            MechanicalDampingDirection {
                damping: &[0.0; 4],
                relaxation_rates: &[0.0; 3],
            },
            0.8,
            dense(),
            || false,
        )
        .unwrap();
    zero.advance_into(&state, &[0.0; 7], &input, &[0.0; 2], &mut tangent)
        .unwrap();
    assert_eq!(tangent, [0.0; 7]);
    let initial_theta = [0.03, -0.02, 0.04, 0.01, -0.02, 0.03, 0.01];
    let input_theta = [0.2, -0.3];
    zero.advance_into(&state, &initial_theta, &input, &input_theta, &mut tangent)
        .unwrap();
    let propagated = native(&source, 0.8, &initial_theta, &input_theta);
    for i in 0..7 {
        close(tangent[i], propagated[i], 4e-13);
    }
}

#[test]
fn reverse_state_force_and_parameter_actions_match_independent_five_point_audit() {
    let source = model(0.0);
    let x = initial(&source);
    let u = [0.7, -0.4];
    let lambda = [0.4, -0.7, 0.2, 0.9, -0.3, 0.8, 0.5];
    let cache = source
        .prepare_damping_sensitivity(direction(), 0.4, dense(), || false)
        .unwrap();
    let mut x_bar = [0.0; 7];
    let mut u_bar = [0.0; 2];
    let parameter = cache
        .pullback_into(&x, &u, &lambda, &mut x_bar, &mut u_bar)
        .unwrap();
    let objective = |theta: f64, state: &[f64], input: &[f64]| -> f64 {
        native(&model(theta), 0.4, state, input)
            .iter()
            .zip(lambda)
            .map(|(a, b)| a * b)
            .sum()
    };
    let fd = |f: &dyn Fn(f64) -> f64| {
        let h = 2e-4;
        (f(-2.0 * h) - 8.0 * f(-h) + 8.0 * f(h) - f(2.0 * h)) / (12.0 * h)
    };
    close(parameter, fd(&|theta| objective(theta, &x, &u)), 2e-10);
    for index in 0..7 {
        close(
            x_bar[index],
            fd(&|step| {
                let mut state = x.clone();
                state[index] += step;
                objective(0.0, &state, &u)
            }),
            2e-10,
        );
    }
    for index in 0..2 {
        close(
            u_bar[index],
            fd(&|step| {
                let mut input = u;
                input[index] += step;
                objective(0.0, &x, &input)
            }),
            2e-10,
        );
    }
    let x_theta = [0.03, -0.02, 0.04, 0.01, -0.02, 0.03, 0.01];
    let u_theta = [0.2, -0.3];
    let mut tangent = [0.0; 7];
    cache
        .advance_into(&x, &x_theta, &u, &u_theta, &mut tangent)
        .unwrap();
    let forward: f64 = tangent.iter().zip(lambda).map(|(a, b)| a * b).sum();
    let reverse = parameter
        + x_bar.iter().zip(x_theta).map(|(a, b)| a * b).sum::<f64>()
        + u_bar.iter().zip(u_theta).map(|(a, b)| a * b).sum::<f64>();
    close(forward, reverse, 3e-15);
}

fn forces() -> Vec<[f64; 2]> {
    (0..24)
        .map(|i| [0.7 * (0.37 * i as f64).cos(), -0.4 + 0.02 * i as f64])
        .collect()
}
fn observations(theta: f64, inputs: &[[f64; 2]]) -> Vec<Vec<f64>> {
    let source = model(theta);
    let flow = PreparedExactZoh::new(source.system(), 0.04, split()).unwrap();
    let mut x = initial(&source);
    inputs
        .iter()
        .map(|u| {
            x = flow.step(&x, u).unwrap().x;
            x.clone()
        })
        .collect()
}
fn fit_objective(theta: f64, inputs: &[[f64; 2]], target: &[Vec<f64>]) -> (f64, f64) {
    let source = model(theta);
    let cache = source
        .prepare_damping_sensitivity(direction(), 0.04, dense(), || false)
        .unwrap();
    let states = observations(theta, inputs);
    let x0 = initial(&source);
    let mut loss = 0.0;
    let mut gradient = 0.0;
    let mut adjoint = [0.0; 7];
    for step in (0..inputs.len()).rev() {
        for row in 0..7 {
            let residual = states[step][row] - target[step][row];
            loss += 0.5 * residual * residual;
            adjoint[row] += residual;
        }
        let mut previous = [0.0; 7];
        gradient += cache
            .pullback_into(
                if step == 0 { &x0 } else { &states[step - 1] },
                &inputs[step],
                &adjoint,
                &mut previous,
                &mut [0.0; 2],
            )
            .unwrap();
        adjoint = previous;
    }
    (loss, gradient)
}

#[test]
fn actual_transient_recovers_shared_correlated_damping_and_relaxation_parameter() {
    let input = forces();
    let truth = 0.73;
    let target = observations(truth, &input);
    let initial_loss = fit_objective(0.1, &input, &target).0;
    let (loss, gradient) = fit_objective(0.25, &input, &target);
    let h = 1e-4;
    let finite_difference = (fit_objective(0.25 - 2.0 * h, &input, &target).0
        - 8.0 * fit_objective(0.25 - h, &input, &target).0
        + 8.0 * fit_objective(0.25 + h, &input, &target).0
        - fit_objective(0.25 + 2.0 * h, &input, &target).0)
        / (12.0 * h);
    assert!(loss > 0.0 && gradient < 0.0);
    close(gradient, finite_difference, 2e-9);
    let (mut low, mut high) = (0.0, 1.4);
    assert!(fit_objective(low, &input, &target).1 < 0.0);
    assert!(fit_objective(high, &input, &target).1 > 0.0);
    // Test-only scalar bracket: production exposes physical actions, no optimizer.
    for _ in 0..38 {
        let mid = (low + high) * 0.5;
        if fit_objective(mid, &input, &target).1 > 0.0 {
            high = mid;
        } else {
            low = mid;
        }
    }
    let fitted = 0.5 * (low + high);
    close(fitted, truth, 2e-9);
    assert!(fit_objective(fitted, &input, &target).0 < initial_loss * 1e-18);
}

#[test]
fn original_physical_coordinate_units_preserve_the_damping_direction() {
    let source = model(0.0);
    let units = [1e-3, 1e3];
    let transform = |values: [f64; 4]| {
        core::array::from_fn::<_, 4, _>(|i| values[i] / units[i / 2] / units[i % 2])
    };
    let ports = core::array::from_fn::<_, 4, _>(|i| PORTS[i] / units[i / 2]);
    let arms = (0..3)
        .map(|i| RelaxationBranch {
            projection: vec![PROJECTIONS[i][0] / units[0], PROJECTIONS[i][1] / units[1]],
            stiffness: BRANCH_K[i],
            relaxation_time_s: 1.0 / RATES[i],
        })
        .collect();
    let rescaled = LinearMechanicalSystem::new(
        &transform(MASS),
        &transform(DAMPING),
        &transform(STIFFNESS),
        &ports,
        2,
        2,
        arms,
    )
    .unwrap();
    let x = initial(&source);
    let mut rescaled_x = x.clone();
    for i in 0..2 {
        rescaled_x[i] *= units[i];
        rescaled_x[2 + i] /= units[i];
    }
    let reference = source
        .prepare_damping_sensitivity(direction(), 0.2, dense(), || false)
        .unwrap();
    let scaled = rescaled
        .prepare_damping_sensitivity(
            MechanicalDampingDirection {
                damping: &transform(DIRECTION),
                relaxation_rates: &RATE_DIRECTION,
            },
            0.2,
            dense(),
            || false,
        )
        .unwrap();
    let mut expected = [0.0; 7];
    let mut actual = [0.0; 7];
    reference
        .advance_into(&x, &[0.0; 7], &[0.7, -0.4], &[0.0; 2], &mut expected)
        .unwrap();
    scaled
        .advance_into(&rescaled_x, &[0.0; 7], &[0.7, -0.4], &[0.0; 2], &mut actual)
        .unwrap();
    for i in 0..2 {
        actual[i] /= units[i];
        actual[2 + i] *= units[i];
    }
    for i in 0..7 {
        close(actual[i], expected[i], 3e-12);
    }
}

#[test]
fn complete_direction_resource_bounds_and_cancellation_refuse_without_state_mutation() {
    let source = model(0.0);
    let state = initial(&source);
    let original = state.clone();
    for invalid in [
        MechanicalDampingDirection {
            damping: &[0.0; 3],
            relaxation_rates: &RATE_DIRECTION,
        },
        MechanicalDampingDirection {
            damping: &DIRECTION,
            relaxation_rates: &[0.0; 2],
        },
        MechanicalDampingDirection {
            damping: &[f64::NAN; 4],
            relaxation_rates: &RATE_DIRECTION,
        },
    ] {
        assert!(
            source
                .prepare_damping_sensitivity(invalid, 0.2, dense(), || false)
                .is_err()
        );
    }
    for duration in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            source
                .prepare_damping_sensitivity(direction(), duration, dense(), || false)
                .is_err()
        );
    }
    for budget in [
        DenseExponentialBudget {
            max_dimension: 17,
            ..dense()
        },
        DenseExponentialBudget {
            max_dense_bytes: 1,
            ..dense()
        },
        DenseExponentialBudget {
            max_multiply_terms: 1,
            ..dense()
        },
    ] {
        assert!(
            source
                .prepare_damping_sensitivity(direction(), 0.2, budget, || false)
                .is_err()
        );
    }
    assert!(matches!(
        source.prepare_damping_sensitivity(direction(), 0.2, dense(), || true),
        Err(PreparedStepError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        source.prepare_damping_sensitivity(direction(), 0.2, dense(), || {
            polls += 1;
            polls == 8
        }),
        Err(PreparedStepError::Cancelled)
    ));
    let cache = source
        .prepare_damping_sensitivity(direction(), 0.2, dense(), || false)
        .unwrap();
    let mut out = [99.0; 7];
    assert!(
        cache
            .advance_into(&state, &[0.0; 7], &[0.0], &[0.0; 2], &mut out)
            .is_err()
    );
    assert_eq!(out, [99.0; 7]);
    assert!(
        cache
            .pullback_into(&state, &[0.0; 2], &[f64::NAN; 7], &mut out, &mut [0.0; 2])
            .is_err()
    );
    assert_eq!(out, [99.0; 7]);
    cache
        .advance_into(&state, &[0.0; 7], &[0.7, -0.4], &[0.0; 2], &mut out)
        .unwrap();
    let mut replay = [0.0; 7];
    cache
        .advance_into(&state, &[0.0; 7], &[0.7, -0.4], &[0.0; 2], &mut replay)
        .unwrap();
    assert_eq!(out, replay);
    assert_eq!(state, original);
}
