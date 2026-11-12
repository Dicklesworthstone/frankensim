//! Original exact heat derivatives, independent physical references and atomic control.

#![allow(clippy::float_cmp)]

use core::cell::Cell;
use std::rc::Rc;

use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh, QuadraticPart};
use fs_phs::{PortHamiltonian, PreparedStepError, QuadraticStorage, Storage};

fn budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 4,
        max_omega_dt: 0.04,
    }
}

fn scalar(resistance: f64) -> PortHamiltonian {
    PortHamiltonian::new(
        1,
        2,
        vec![0.0],
        vec![resistance],
        vec![2.0, -3.0],
        Box::new(QuadraticStorage::new(vec![4.0], 1).unwrap()),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} versus {expected:.17e}, tolerance {tolerance:.3e}"
    );
}

#[test]
fn full_two_port_derivative_matches_independent_scalar_exponential_and_heat() {
    let system = scalar(0.25);
    let base = PreparedExactZoh::new(&system, 0.5, budget()).unwrap();
    let mut workspace = base.workspace();
    for (duration, initial, input) in [
        (1.25, [0.12], [0.17, -0.23]),
        (0.25, [-0.4], [-0.3, 0.7]),
        (0.03125, [0.8], [0.7, 0.2]),
        (1.25, [0.12], [0.17, -0.23]),
    ] {
        let owner =
            PreparedExactZoh::new_controlled(&system, duration, budget(), || false).unwrap();
        let (mut next, mut output, mut gradient) = ([0.0], [0.0; 2], [f64::NAN; 2]);
        let record = owner
            .step_with_workspace_into_controlled(
                &mut workspace,
                &initial,
                &input,
                &mut next,
                &mut output,
                || false,
            )
            .unwrap();
        let magnitude = owner.port_work_magnitude(&workspace, &input).unwrap();
        assert!(magnitude >= record.supplied.abs() * (1.0 - 8.0 * f64::EPSILON));
        let heat = owner
            .heat_input_gradient_into_controlled(
                &mut workspace,
                &initial,
                &input,
                &mut gradient,
                || false,
            )
            .unwrap();
        assert_eq!(heat.to_bits(), record.dissipated.to_bits());

        // Independent physical equation x'=-x+2*u0-3*u1, H=2*x².
        let equilibrium = 2.0 * input[0] - 3.0 * input[1];
        let c = initial[0] - equilibrium;
        let a = -(-duration).exp_m1();
        let b = -(-2.0 * duration).exp_m1() / 2.0;
        let mean = equilibrium + c * a / duration;
        close(
            magnitude,
            duration * ((input[0] * 8.0 * mean).abs() + (input[1] * -12.0 * mean).abs()),
            3.0e-12,
        );
        let expected_heat =
            4.0 * (equilibrium * equilibrium * duration + 2.0 * equilibrium * c * a + c * c * b);
        close(next[0], equilibrium + c * (-duration).exp(), 3.0e-13);
        close(
            heat,
            expected_heat,
            3.0e-12 * expected_heat.abs().max(1.0e-3),
        );
        for (port, g) in [2.0, -3.0].into_iter().enumerate() {
            close(output[port], 4.0 * g * mean, 3.0e-12);
            let expected = 8.0 * g * (equilibrium * duration + (c - equilibrium) * a - c * b);
            close(
                gradient[port],
                expected,
                4.0e-12 * expected.abs().max(1.0e-3),
            );
            let (_, response) = owner.average_output_map(port).unwrap();
            for (column, drive) in [2.0, -3.0].into_iter().enumerate() {
                close(
                    response[column],
                    4.0 * g * drive * (1.0 - a / duration),
                    3.0e-12,
                );
            }
        }
        close(
            record.delta_h,
            2.0 * (next[0] * next[0] - initial[0] * initial[0]),
            3.0e-13,
        );
        close(
            record.supplied,
            duration * (input[0] * output[0] + input[1] * output[1]),
            3.0e-12,
        );
    }
}

fn integrate(f: impl Fn(f64) -> f64, duration: f64) -> f64 {
    let count = 4096;
    let h = duration / f64::from(count);
    let mut sum = f(0.0) + f(duration);
    for i in 1..count {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(f64::from(i) * h);
    }
    sum * h / 3.0
}

#[test]
fn positive_factor_retains_weak_heat_and_derivative_near_force_supported_silence() {
    let system = PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, 0.2],
        vec![0.0, 0.0, 1.0, -1.0],
        Box::new(QuadraticStorage::new(vec![1.0, 0.0, 0.0, 1.0], 2).unwrap()),
    )
    .unwrap();
    for duration in [0.1, 0.4] {
        let owner = PreparedExactZoh::new(&system, duration, budget()).unwrap();
        let mut workspace = owner.workspace();
        let mut gradient = [31.0, 37.0];
        let silence = owner
            .heat_input_gradient_into_controlled(
                &mut workspace,
                &[1.0, 0.0],
                &[2.0, 1.0],
                &mut gradient,
                || false,
            )
            .unwrap();
        assert_eq!(silence, 0.0);
        assert_eq!(gradient, [0.0, 0.0]);

        let epsilon = 1.0e-10;
        let heat = owner
            .heat_input_gradient_into_controlled(
                &mut workspace,
                &[1.0, epsilon],
                &[2.0, 1.0],
                &mut gradient,
                || false,
            )
            .unwrap();
        // Original oscillator: q''+0.2*q'+q=1, initial q=1, q'=epsilon.
        // The input-step sensitivity of velocity is exp(-0.1*t)*sin(w*t)/w.
        let omega = 0.99_f64.sqrt();
        let velocity = |t: f64| {
            epsilon * (-0.1 * t).exp() * ((omega * t).cos() - 0.1 / omega * (omega * t).sin())
        };
        let response = |t: f64| (-0.1 * t).exp() * (omega * t).sin() / omega;
        let expected_heat = integrate(|t| 0.2 * velocity(t).powi(2), duration);
        let expected_gradient = integrate(|t| 0.4 * velocity(t) * response(t), duration);
        assert!(heat > 0.0 && expected_gradient > 0.0);
        close(heat / expected_heat, 1.0, 3.0e-11);
        close(gradient[0] / expected_gradient, 1.0, 3.0e-11);
        close(gradient[1] / expected_gradient, -1.0, 3.0e-11);
    }
}

#[test]
fn coupled_storage_and_off_diagonal_resistance_keep_all_input_columns() {
    let system = PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.4, 0.1, 0.1, 0.5],
        vec![1.0, 0.2, -0.3, 0.7],
        Box::new(QuadraticStorage::new(vec![2.0, 0.3, 0.3, 3.0], 2).unwrap()),
    )
    .unwrap();
    let owner = PreparedExactZoh::new(&system, 0.2, budget()).unwrap();
    let mut workspace = owner.workspace();
    let initial = [0.4, -0.7];
    let input = [0.3, -0.2];
    let mut gradient = [0.0; 2];
    owner
        .heat_input_gradient_into_controlled(
            &mut workspace,
            &initial,
            &input,
            &mut gradient,
            || false,
        )
        .unwrap();
    for port in 0..2 {
        let (mut plus, mut minus) = (input, input);
        let delta = 1.0e-5;
        plus[port] += delta;
        minus[port] -= delta;
        let difference = (owner.step(&initial, &plus).unwrap().dissipated
            - owner.step(&initial, &minus).unwrap().dissipated)
            / (2.0 * delta);
        close(gradient[port], difference, 3.0e-10);
    }
}

#[test]
fn zero_duration_and_zero_resistance_give_zero_heat_derivatives() {
    for (resistance, duration) in [(0.25, 0.0), (0.0, 0.25)] {
        let owner =
            PreparedExactZoh::new_owned_controlled(scalar(resistance), duration, budget(), || {
                false
            })
            .unwrap();
        let mut workspace = owner.workspace();
        let mut gradient = [f64::NAN; 2];
        let heat = owner
            .heat_input_gradient_into_controlled(
                &mut workspace,
                &[0.3],
                &[0.7, -0.2],
                &mut gradient,
                || false,
            )
            .unwrap();
        assert_eq!(heat, 0.0);
        assert_eq!(gradient, [0.0, 0.0]);
    }
}

#[test]
fn cancellation_at_every_exposed_boundary_is_atomic_and_reusable() {
    let system = scalar(0.25);
    let calls = Cell::new(0);
    PreparedExactZoh::new_controlled(&system, 0.25, budget(), || {
        calls.set(calls.get() + 1);
        false
    })
    .unwrap();
    for stop in 1..=calls.get() {
        let count = Cell::new(0);
        assert!(matches!(
            PreparedExactZoh::new_controlled(&system, 0.25, budget(), || {
                count.set(count.get() + 1);
                count.get() == stop
            }),
            Err(PreparedStepError::Cancelled)
        ));
    }
    let owner = PreparedExactZoh::new(&system, 0.25, budget()).unwrap();
    let mut workspace = owner.workspace();
    let initial = [0.3];
    let input = [0.7, -0.2];
    let reference = owner.step(&initial, &input).unwrap();
    let mut reference_gradient = [0.0; 2];
    owner
        .heat_input_gradient_into_controlled(
            &mut workspace,
            &initial,
            &input,
            &mut reference_gradient,
            || false,
        )
        .unwrap();
    for tangent in [false, true] {
        calls.set(0);
        let (mut next, mut output) = ([0.0], [0.0; 2]);
        let mut count = || {
            calls.set(calls.get() + 1);
            false
        };
        if tangent {
            owner
                .heat_input_gradient_into_controlled(
                    &mut workspace,
                    &initial,
                    &input,
                    &mut output,
                    &mut count,
                )
                .unwrap();
        } else {
            owner
                .step_with_workspace_into_controlled(
                    &mut workspace,
                    &initial,
                    &input,
                    &mut next,
                    &mut output,
                    &mut count,
                )
                .unwrap();
        }
        for stop in 1..=calls.get() {
            let count = Cell::new(0);
            let cancelled = || {
                count.set(count.get() + 1);
                count.get() == stop
            };
            let (mut next, mut output) = ([17.0], [23.0, 29.0]);
            let result = if tangent {
                owner
                    .heat_input_gradient_into_controlled(
                        &mut workspace,
                        &initial,
                        &input,
                        &mut output,
                        cancelled,
                    )
                    .map(|_| ())
            } else {
                owner
                    .step_with_workspace_into_controlled(
                        &mut workspace,
                        &initial,
                        &input,
                        &mut next,
                        &mut output,
                        cancelled,
                    )
                    .map(|_| ())
            };
            assert!(matches!(result, Err(PreparedStepError::Cancelled)));
            assert_eq!(next, [17.0]);
            assert_eq!(output, [23.0, 29.0]);
            owner
                .step_with_workspace_into_controlled(
                    &mut workspace,
                    &initial,
                    &input,
                    &mut next,
                    &mut output,
                    || false,
                )
                .unwrap();
            assert_eq!(next.as_slice(), reference.x.as_slice());
            assert_eq!(output.as_slice(), reference.y.as_slice());
            owner
                .heat_input_gradient_into_controlled(
                    &mut workspace,
                    &initial,
                    &input,
                    &mut output,
                    || false,
                )
                .unwrap();
            assert_eq!(output, reference_gradient);
        }
    }
    let mut output = [23.0, 29.0];
    for (state, input) in [([f64::NAN], input), (initial, [f64::INFINITY, 0.0])] {
        assert!(matches!(
            owner.heat_input_gradient_into_controlled(
                &mut workspace,
                &state,
                &input,
                &mut output,
                || false
            ),
            Err(PreparedStepError::Solver(_))
        ));
        assert_eq!(output, [23.0, 29.0]);
    }
    assert!(
        owner
            .heat_input_gradient_into_controlled(
                &mut workspace,
                &initial,
                &input,
                &mut output[..1],
                || false
            )
            .is_err()
    );
    assert_eq!(output, [23.0, 29.0]);
    assert!(owner.port_work_magnitude(&workspace, &[0.0]).is_err());
    assert!(matches!(
        PreparedExactZoh::new_owned_controlled(
            scalar(0.25),
            0.1,
            LinearSplitBudget {
                max_modes: 0,
                max_omega_dt: 0.04
            },
            || false
        ),
        Err(PreparedStepError::Solver(_))
    ));
}

struct CancellingStorage {
    armed: Rc<Cell<bool>>,
    requested: Rc<Cell<bool>>,
}

impl CancellingStorage {
    fn request(&self) {
        if self.armed.replace(false) {
            self.requested.set(true);
        }
    }
}

impl Storage for CancellingStorage {
    fn hamiltonian(&self, state: &[f64]) -> f64 {
        0.5 * state[0] * state[0]
    }
    fn gradient(&self, state: &[f64], output: &mut [f64]) {
        output[0] = state[0];
        self.request();
    }
    fn quadratic_part(&self) -> Option<QuadraticPart> {
        self.request();
        Some(QuadraticPart {
            matrix: vec![1.0],
            exact: true,
        })
    }
}

#[test]
fn cancellation_requested_inside_original_storage_never_publishes_outputs() {
    let armed = Rc::new(Cell::new(false));
    let requested = Rc::new(Cell::new(false));
    let system = PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![0.25],
        vec![1.0],
        Box::new(CancellingStorage {
            armed: Rc::clone(&armed),
            requested: Rc::clone(&requested),
        }),
    )
    .unwrap();
    armed.set(true);
    assert!(matches!(
        PreparedExactZoh::new_controlled(&system, 0.2, budget(), || requested.get()),
        Err(PreparedStepError::Cancelled)
    ));
    assert!(requested.get() && !armed.get());
    requested.set(false);
    let owner =
        PreparedExactZoh::new_controlled(&system, 0.2, budget(), || requested.get()).unwrap();
    let mut workspace = owner.workspace();
    for tangent in [false, true] {
        armed.set(true);
        let (mut next, mut output) = ([17.0], [23.0]);
        let result = if tangent {
            owner
                .heat_input_gradient_into_controlled(
                    &mut workspace,
                    &[0.3],
                    &[0.7],
                    &mut output,
                    || requested.get(),
                )
                .map(|_| ())
        } else {
            owner
                .step_with_workspace_into_controlled(
                    &mut workspace,
                    &[0.3],
                    &[0.7],
                    &mut next,
                    &mut output,
                    || requested.get(),
                )
                .map(|_| ())
        };
        assert!(matches!(result, Err(PreparedStepError::Cancelled)));
        assert!(requested.get() && !armed.get());
        assert_eq!(next, [17.0]);
        assert_eq!(output, [23.0]);
        requested.set(false);
        owner
            .step_with_workspace_into_controlled(
                &mut workspace,
                &[0.3],
                &[0.7],
                &mut next,
                &mut output,
                || requested.get(),
            )
            .unwrap();
        owner
            .heat_input_gradient_into_controlled(
                &mut workspace,
                &[0.3],
                &[0.7],
                &mut output,
                || requested.get(),
            )
            .unwrap();
        assert!(output[0].is_finite());
    }
}
