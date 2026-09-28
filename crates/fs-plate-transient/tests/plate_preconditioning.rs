//! G0/G3/G4/G5: physical forward scaling, bounded Krylov solves and replayed AD.
#[path = "../examples/preconditioned_plate.rs"]
#[allow(dead_code)]
mod comparison;

use comparison::moving;
use fs_plate::{PlateMesh, PlateModel};
use fs_plate_transient::{PlateDynamics, PlateLoad};
use fs_solver::{NewtonStallDiagnosis, StallDiagnosis};
use fs_time::galpha::{
    GeneralizedAlpha, SecondOrderOperatorWeights, SecondOrderProblem, SecondOrderState,
    TimeSolveError, galpha_step,
    initialization::{initial_workspace_components, second_order_acceleration_vjp},
    second_order_adjoint::{
        SecondOrderAdjointConfig,
        trajectory::{
            RecordedStructural, StructuralRecordingConfig, StructuralRecordingStatus,
            StructuralReplayBudget, StructuralTrajectoryError, StructuralTrajectoryModel,
            samples::StructuralSampleObjective,
        },
    },
};

fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn close_vector(actual: &[f64], expected: &[f64]) {
    let scale = expected.iter().fold(1e-12_f64, |m, v| m.max(v.abs()));
    let error = actual
        .iter()
        .zip(expected)
        .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
    assert!(error <= 2e-8 * scale, "error {error:e}, scale {scale:e}");
}

// Independent dense assembly of the complete physical M/C/K pencil and LU
// stepping. This does not call PlateDynamics or any preconditioner/adjoint.
fn dense_run(
    mesh: &PlateMesh,
    plate: &PlateModel,
    p: [f64; 5],
    step: f64,
    steps: usize,
    samples: &[usize],
) -> (SecondOrderState, f64) {
    let n = plate.free;
    let mass: Vec<f64> = plate.m.to_dense().iter().map(|v| p[1] * v).collect();
    let stiffness: Vec<f64> = plate.k.to_dense().iter().map(|v| p[0] * v).collect();
    let damping: Vec<f64> = mass
        .iter()
        .zip(&stiffness)
        .map(|(m, k)| p[2] * m + p[3] * k)
        .collect();
    let method = GeneralizedAlpha::new(&mass, &damping, &stiffness, n, step, comparison::RHO);
    let load = moving::MovingLoad {
        mesh,
        plate,
        amplitude: p[4],
    };
    let mut force = vec![0.0; n];
    load.forcing(0.0, &mut force).unwrap();
    let mut state = SecondOrderState::new(
        0.0,
        &vec![0.0; n],
        &vec![0.0; n],
        &(0..n)
            .map(|i| force[i] / mass[i * n + i])
            .collect::<Vec<_>>(),
    );
    let weights = moving::sensor_weights(mesh, plate).unwrap();
    let alpha_f = comparison::RHO / (1.0 + comparison::RHO);
    let mut values = Vec::new();
    for endpoint in 0..=steps {
        if endpoint != 0 {
            load.forcing(step.mul_add(1.0 - alpha_f, state.t), &mut force)
                .unwrap();
            galpha_step(
                &method,
                &mut state.q,
                &mut state.v,
                &mut state.a,
                &force,
            );
            state.t += step;
            state.steps += 1;
        }
        for &sample in samples {
            if sample == endpoint {
                values.push(0.5 * (dot(&weights, &state.q) / moving::DEFLECTION_SCALE_M).powi(2));
            }
        }
    }
    (state, values.into_iter().rev().sum())
}

#[test]
fn forward_hook_inverts_the_physical_mass_damping_and_tangent_diagonal() {
    let (mesh, plate) = comparison::build_plate(4).unwrap();
    let p = comparison::PARAMETERS;
    let load = moving::MovingLoad {
        mesh: &mesh,
        plate: &plate,
        amplitude: p[4],
    };
    let dynamics = PlateDynamics::new(
        &plate,
        moving::dynamics_parameters(p),
        &load,
        moving::dynamics_budget(&plate),
        &mut || false,
    )
    .unwrap();
    let q: Vec<f64> = (0..plate.free).map(|i| 1e-5 * (i as f64).sin()).collect();
    let rhs: Vec<f64> = (0..plate.free).map(|i| 0.2 + (i as f64).cos()).collect();
    let step = comparison::STEP;
    let rho = comparison::RHO;
    let alpha_m = (2.0 * rho - 1.0) / (1.0 + rho);
    let alpha_f = rho / (1.0 + rho);
    let gamma = 0.5 - alpha_m + alpha_f;
    let beta = 0.25 * (1.0 - alpha_m + alpha_f).powi(2);
    let policies = [
        SecondOrderOperatorWeights {
            mass: 1.0,
            damping: 0.0,
            tangent: 0.0,
        },
        SecondOrderOperatorWeights {
            mass: 0.0,
            damping: 1.0,
            tangent: 0.0,
        },
        SecondOrderOperatorWeights {
            mass: 0.0,
            damping: 0.0,
            tangent: 1.0,
        },
        SecondOrderOperatorWeights {
            mass: (1.0 - alpha_m) / (beta * step * step),
            damping: (1.0 - alpha_f) * gamma / (beta * step),
            tangent: 1.0 - alpha_f,
        },
    ];
    for weights in policies {
        let mut actual = vec![f64::NAN; plate.free];
        dynamics.preconditioner_apply(&q, weights, 3, 7, &rhs, &mut actual);
        for i in 0..plate.free {
            let mass = p[1] * plate.m.get(i, i);
            let stiffness = p[0] * plate.k.get(i, i);
            let damping = p[2] * mass + p[3] * stiffness;
            let expected = rhs[i]
                / (weights.mass * mass + weights.damping * damping + weights.tangent * stiffness);
            assert!((actual[i] - expected).abs() <= 8e-16 * expected.abs());
        }
        comparison::IdentityProblem(&dynamics)
            .preconditioner_apply(&q, weights, 3, 7, &rhs, &mut actual);
        assert_eq!(actual, rhs, "baseline must retain the trait's identity hook");
    }
    let mut refused = vec![0.0; plate.free];
    dynamics.preconditioner_apply(
        &q,
        SecondOrderOperatorWeights {
            mass: -1.0,
            damping: 0.0,
            tangent: 0.0,
        },
        0,
        0,
        &rhs,
        &mut refused,
    );
    assert!(refused.iter().all(|v| v.is_nan()));
}

#[test]
fn larger_plates_converge_with_short_restart_while_identity_exhausts_the_same_cap() {
    for elements in [8, 12] {
        let (mesh, plate) = comparison::build_plate(elements).unwrap();
        assert_eq!(plate.free, 3 * (elements - 1).pow(2));
        let p = comparison::PARAMETERS;
        let load = moving::MovingLoad {
            mesh: &mesh,
            plate: &plate,
            amplitude: p[4],
        };
        let dynamics = PlateDynamics::new(
            &plate,
            moving::dynamics_parameters(p),
            &load,
            moving::dynamics_budget(&plate),
            &mut || false,
        )
        .unwrap();
        let initial = moving::resting_state(&dynamics).unwrap();
        let method = comparison::method(
            plate.free,
            comparison::STEP,
            comparison::RESTART,
            comparison::CYCLES,
        );
        let mut forcing = vec![0.0; plate.free];
        dynamics
            .forcing(method.forcing_time(initial.t).unwrap(), &mut forcing)
            .unwrap();
        let mut actual = initial.clone();
        let report = method.step(&mut actual, &dynamics, &forcing).unwrap();
        let inner: usize = report
            .newton
            .history
            .iter()
            .map(|i| i.linear_iterations)
            .sum();
        eprintln!(
            "{elements}x{elements}: dofs={}, inner={inner}, residual={:e}",
            plate.free, report.newton.residual_norm
        );
        assert!(report.newton.converged && report.newton.residual_norm < 1e-11);
        assert!(inner > 0 && inner <= comparison::RESTART * comparison::CYCLES);
        assert!(
            report
                .newton
                .history
                .iter()
                .all(|i| i.linear_iterations <= comparison::RESTART * comparison::CYCLES)
        );
        let mut baseline = initial.clone();
        let failure = method.step(
            &mut baseline,
            &comparison::IdentityProblem(&dynamics),
            &forcing,
        );
        assert!(matches!(failure, Err(TimeSolveError::NotConverged(ref r))
            if r.diagnosis == Some(NewtonStallDiagnosis::LinearSolveFailed(StallDiagnosis::BudgetExhausted))));
        assert_eq!(baseline, initial, "a refused step must not publish any state");
        // Dense LU is a separate algorithm and retains all physical DOFs. The
        // smaller of these two larger meshes suffices for the endpoint oracle.
        if elements == 8 {
            let (expected, _) = dense_run(&mesh, &plate, p, comparison::STEP, 1, &[]);
            close_vector(&actual.q, &expected.q);
            close_vector(&actual.v, &expected.v);
            close_vector(&actual.a, &expected.a);
        }
        let mut cancelled = actual.clone();
        assert!(matches!(
            method.step_controlled(&mut cancelled, &dynamics, &forcing, &mut || true),
            Err(TimeSolveError::Cancelled)
        ));
        assert_eq!(cancelled, actual);
    }
}

struct Objective<'a>(&'a [f64]);
impl StructuralSampleObjective for Objective<'_> {
    fn evaluate(
        &self,
        _sample: usize,
        state: &SecondOrderState,
        bars: (&mut [f64], &mut [f64], &mut [f64]),
        parameters: &mut [f64],
    ) -> Result<f64, String> {
        let value = dot(self.0, &state.q) / moving::DEFLECTION_SCALE_M;
        for (bar, weight) in bars.0.iter_mut().zip(self.0) {
            *bar = value * weight / moving::DEFLECTION_SCALE_M;
        }
        bars.1.fill(0.0);
        bars.2.fill(0.0);
        parameters.fill(0.0);
        Ok(0.5 * value * value)
    }
}

#[test]
#[allow(clippy::too_many_lines)] // One complete sparse forward/replay/initial-chain experiment.
fn preconditioned_checkpointed_trajectory_matches_dense_finite_difference() {
    let (mesh, plate) = comparison::build_plate(8).unwrap();
    let p = comparison::PARAMETERS;
    let load = moving::MovingLoad {
        mesh: &mesh,
        plate: &plate,
        amplitude: p[4],
    };
    let dynamics = PlateDynamics::new(
        &plate,
        moving::dynamics_parameters(p),
        &load,
        moving::dynamics_budget(&plate),
        &mut || false,
    )
    .unwrap();
    let initial = moving::resting_state(&dynamics).unwrap();
    let method = comparison::method(
        plate.free,
        comparison::STEP,
        comparison::RESTART,
        comparison::CYCLES,
    );
    let adjoint = SecondOrderAdjointConfig {
        restart: comparison::RESTART,
        max_cycles: comparison::CYCLES,
        tolerance: 1e-12,
    };
    let mut tape = RecordedStructural::new(
        method,
        &dynamics,
        &initial,
        StructuralRecordingConfig {
            steps: moving::STEPS,
            adjoint,
            max_workspace_components: method.adjoint_workspace_components(5, adjoint).unwrap(),
        },
    )
    .unwrap();
    assert_eq!(
        tape.advance(5, moving::STEPS, &mut || false)
            .unwrap()
            .status,
        StructuralRecordingStatus::StepLimit
    );
    assert_eq!(
        tape.advance(moving::STEPS, moving::STEPS, &mut || false)
            .unwrap()
            .status,
        StructuralRecordingStatus::ReachedEnd
    );
    let weights = moving::sensor_weights(&mesh, &plate).unwrap();
    let jacobi = dynamics
        .effective_preconditioner(comparison::STEP, comparison::RHO, &mut || false)
        .unwrap();
    let budget = StructuralReplayBudget {
        checkpoints: tape.required_checkpoints(),
        forward_steps: 36,
    };
    let endpoint = tape.state().clone();
    let result = tape
        .pullback_samples(
            &moving::SAMPLES,
            moving::SAMPLES.len(),
            &Objective(&weights),
            &jacobi,
            budget,
            &mut || false,
        )
        .unwrap();
    assert_eq!(result.gradient.replayed_steps, 36);
    assert_eq!(result.gradient.peak_checkpoints, 4);
    assert_eq!(tape.state(), &endpoint);
    // A short replay allowance refuses atomically; a complete retry remains
    // bit-identical, including the initial-acceleration cotangent.
    assert!(matches!(
        tape.pullback_samples(
            &moving::SAMPLES,
            moving::SAMPLES.len(),
            &Objective(&weights),
            &jacobi,
            StructuralReplayBudget {
                forward_steps: 1,
                ..budget
            },
            &mut || false
        ),
        Err(StructuralTrajectoryError::ReplayLimit)
    ));
    assert_eq!(tape.state(), &endpoint);
    let retry = tape
        .pullback_samples(
            &moving::SAMPLES,
            moving::SAMPLES.len(),
            &Objective(&weights),
            &jacobi,
            budget,
            &mut || false,
        )
        .unwrap();
    assert_eq!(result.value.to_bits(), retry.value.to_bits());
    assert_eq!(result.gradient.parameters, retry.gradient.parameters);
    assert_eq!(result.gradient.initial_a, retry.gradient.initial_a);
    let mass = dynamics.mass_preconditioner();
    let mut forcing = vec![0.0; plate.free];
    dynamics.forcing(0.0, &mut forcing).unwrap();
    let initial_gradient = second_order_acceleration_vjp(
        &dynamics,
        &initial.q,
        &initial.v,
        &forcing,
        &result.gradient.initial_a,
        &mass,
        &mass,
        moving::initial_config(),
        initial_workspace_components(plate.free, 5, moving::initial_config()).unwrap(),
        &mut || false,
    )
    .unwrap();
    let mut force_gradient = vec![0.0; 5];
    dynamics
        .forcing_vjp(0.0, &initial_gradient.forcing, &mut force_gradient)
        .unwrap();
    let total: Vec<f64> = (0..5)
        .map(|i| result.gradient.parameters[i] + initial_gradient.parameters[i] + force_gradient[i])
        .collect();
    assert!(
        force_gradient[4].abs() > 1e-12,
        "the initial force chain must actually contribute"
    );
    let (expected, value) = dense_run(
        &mesh,
        &plate,
        p,
        comparison::STEP,
        moving::STEPS,
        &moving::SAMPLES,
    );
    close_vector(&endpoint.q, &expected.q);
    close_vector(&endpoint.v, &expected.v);
    close_vector(&endpoint.a, &expected.a);
    assert!((result.value - value).abs() < 1e-8 * value);
    for (i, delta) in [2e-5, 2e-5, 5e-3, 2e-8, 2e-5].into_iter().enumerate() {
        let (mut plus, mut minus) = (p, p);
        plus[i] += delta;
        minus[i] -= delta;
        let finite_difference = (dense_run(
            &mesh,
            &plate,
            plus,
            comparison::STEP,
            moving::STEPS,
            &moving::SAMPLES,
        )
        .1 - dense_run(
            &mesh,
            &plate,
            minus,
            comparison::STEP,
            moving::STEPS,
            &moving::SAMPLES,
        )
        .1) / (2.0 * delta);
        eprintln!(
            "parameter{i}: adjoint={:e}, dense_fd={finite_difference:e}",
            total[i]
        );
        assert!(finite_difference.abs() > 1e-12);
        assert!(
            (total[i] - finite_difference).abs() <= 5e-5 * finite_difference.abs(),
            "parameter {i}: adjoint {}, dense FD {finite_difference}",
            total[i]
        );
    }
}
