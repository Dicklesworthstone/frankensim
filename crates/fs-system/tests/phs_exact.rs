//! Exact native physics, heat derivatives, duration changes and atomic trials.

#![allow(clippy::float_cmp)]

use core::cell::Cell;
use core::num::NonZeroUsize;
use std::rc::Rc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh, QuadraticPart};
use fs_phs::{PortHamiltonian, PreparedStepError, QuadraticStorage, Storage};
use fs_ports::component::conformance::{Case, Observed, check_component, check_port_tangent};
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal};
use fs_ports::{PortKind, StableId};
use fs_system::elements::HeatCapacity;
use fs_system::phs::PhsComponentError;
use fs_system::phs_exact::ExactPortHamiltonianComponent;
use fs_system::{AdaptiveConfig, AdaptiveError, NodeId, System, SystemBuilder};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 79,
                kernel_id: 23,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 64,
        max_omega_dt: 1.0,
    }
}

fn power(id: &str, role: PortRole) -> PortSpec {
    PortSpec::power(
        StableId::new(id).unwrap(),
        PortKind::MechanicalForceVelocity,
        role,
        NonZeroUsize::MIN,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} versus {expected}, tolerance {tolerance}"
    );
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn coupled_owner() -> PortHamiltonian {
    PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.4, 0.1, 0.1, 0.5],
        vec![1.0, 0.2, -0.3, 0.7],
        Box::new(QuadraticStorage::new(vec![2.0, 0.3, 0.3, 3.0], 2).unwrap()),
    )
    .unwrap()
}

fn scalar_owner() -> PortHamiltonian {
    PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![2.0],
        vec![1.0],
        Box::new(QuadraticStorage::new(vec![3.0], 1).unwrap()),
    )
    .unwrap()
}

#[test]
fn exact_multiport_proposals_match_original_owner_and_complete_tangents_conform() {
    let component = ExactPortHamiltonianComponent::new(
        coupled_owner(),
        vec![
            power("force", PortRole::EffortIn),
            power("motion", PortRole::FlowIn),
        ],
        budget(),
    )
    .unwrap();
    assert_eq!(component.class(), Class::Smooth);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for (old, dt, input) in [
            (vec![0.4, -0.7], 0.1, [0.3, -0.2, 310.0]),
            (vec![-0.3, 0.2], 0.025, [-0.7, 0.6, 290.0]),
            (vec![0.4, -0.7], 0.1, [0.3, -0.2, 310.0]),
        ] {
            let stale = vec![31.0, -11.0];
            let case = Case {
                state: &old,
                stale: &stale,
                dt_s: dt,
                inputs: &input,
            };
            let Observed::Stepped {
                receipt, outputs, ..
            } = check_component(&component, &case, cx).unwrap()
            else {
                panic!("exact proposal unexpectedly shortened");
            };
            check_port_tangent(&component, &case, 2.0e-8, 2.0e-10, cx).unwrap();
            let native = PreparedExactZoh::new(component.system(), dt, budget()).unwrap();
            let expected = native.step(&old, &input[..2]).unwrap();
            let (mut next, mut actual) = (Vec::new(), [f64::NAN; 3]);
            component
                .propose(&old, dt, &input, &mut next, &mut actual, cx)
                .unwrap();
            assert_eq!(bits(&next), bits(&expected.x));
            assert_eq!(bits(&outputs[..2]), bits(&expected.y));
            assert_eq!(bits(&actual), bits(&outputs));
            assert_eq!(
                receipt.dissipated_j.to_bits(),
                expected.dissipated.to_bits()
            );
            assert_eq!(
                receipt.stored_change_j.to_bits(),
                expected.delta_h.to_bits()
            );
            assert_eq!(
                receipt.numerical_loss_j.to_bits(),
                (expected.supplied - expected.dissipated - expected.delta_h).to_bits()
            );
            assert_eq!(outputs[2].to_bits(), (-expected.dissipated / dt).to_bits());
            assert_eq!(
                receipt.entropy_generated_j_per_k.to_bits(),
                (expected.dissipated / input[2]).to_bits()
            );
        }
    });
}

/// Independent closed integrals for x'=-6x+u, y=3x, D=integral 18x^2 dt.
fn scalar_oracle(initial: f64, input: f64, duration: f64) -> (f64, f64, f64) {
    let equilibrium = input / 6.0;
    let amplitude = initial - equilibrium;
    let decay = (-6.0 * duration).exp();
    let integral_decay = (1.0 - decay) / 6.0;
    let integral_decay_squared = (1.0 - decay * decay) / 12.0;
    let mean = equilibrium + amplitude * integral_decay / duration;
    let heat = 18.0
        * (equilibrium * equilibrium * duration
            + 2.0 * equilibrium * amplitude * integral_decay
            + amplitude * amplitude * integral_decay_squared);
    (equilibrium + amplitude * decay, mean, heat)
}

fn heated_device() -> (System, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let component = ExactPortHamiltonianComponent::new(
        scalar_owner(),
        vec![power("p", PortRole::EffortIn)],
        budget(),
    )
    .unwrap();
    let owner = builder
        .add_with_tangent("owner", component, vec![0.8])
        .unwrap();
    let body = builder
        .add_with_tangent("body", HeatCapacity::new(10.0).unwrap(), 300.0)
        .unwrap();
    builder
        .boundary(
            builder.port(owner, "p").unwrap(),
            Box::new(|_, _, values| values[0] = 0.9),
        )
        .unwrap();
    let heat = builder.port(owner, "heat").unwrap();
    builder
        .bond(builder.port(body, "heat").unwrap(), heat)
        .unwrap();
    let mut system = builder.seal().unwrap();
    assert!(system.guess(heat, &[300.0]));
    (system, owner, body)
}

#[test]
fn exact_driven_storage_heats_a_body_under_variable_durations() {
    let (mut system, owner, body) = heated_device();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for dt in [0.02, 0.05, 0.005, 0.025, 0.1] {
            system.step(dt, cx).unwrap();
        }
    });
    let (state, mean, heat) = scalar_oracle(0.8, 0.9, 0.2);
    close(system.time_s(), 0.2, 1.0e-15);
    close(system.state::<Vec<f64>>(owner).unwrap()[0], state, 2.0e-13);
    close(
        system.ledger().nodes[owner.index()].dissipated_j.value(),
        heat,
        2.0e-12,
    );
    close(
        system.ledger().nodes[body.index()].stored_change_j.value(),
        heat,
        2.0e-11,
    );
    close(
        10.0 * (system.state::<f64>(body).unwrap() - 300.0),
        heat,
        2.0e-10,
    );
    close(
        system.ledger().boundary_energy_in_j(),
        3.0 * 0.9 * mean * 0.2,
        2.0e-12,
    );
    close(system.ledger().energy_residual_j(), 0.0, 2.0e-11);
}

#[test]
fn heat_derivative_matches_driven_integral_and_temperature_column_is_zero() {
    let component = ExactPortHamiltonianComponent::new(
        scalar_owner(),
        vec![power("p", PortRole::EffortIn)],
        budget(),
    )
    .unwrap();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for dt in [0.2, 0.05, 0.1] {
            let (initial, input) = (0.8, 0.9);
            let decay = (-6.0_f64 * dt).exp();
            let a = (1.0 - decay) / 6.0;
            let b = (1.0 - decay * decay) / 12.0;
            let mean_derivative = 0.5 * (1.0 - a / dt);
            let heat_derivative =
                6.0 * ((input / 6.0) * (dt - a) + (initial - input / 6.0) * (a - b));
            let mut jvp = [31.0, 37.0];
            component
                .jvp(
                    &vec![initial],
                    dt,
                    &[input, 300.0],
                    &[1.0, 11.0],
                    &mut jvp,
                    cx,
                )
                .unwrap();
            close(jvp[0], mean_derivative, 2.0e-12);
            close(jvp[1], -heat_derivative / dt, 2.0e-11);
            let mut vjp = [41.0, 43.0];
            component
                .vjp(
                    &vec![initial],
                    dt,
                    &[input, 300.0],
                    &[0.4, -0.7],
                    &mut vjp,
                    cx,
                )
                .unwrap();
            close(
                vjp[0],
                0.4 * mean_derivative + 0.7 * heat_derivative / dt,
                2.0e-11,
            );
            assert_eq!(vjp[1].to_bits(), 0.0_f64.to_bits());
            component
                .jvp(
                    &vec![initial],
                    dt,
                    &[input, 300.0],
                    &[0.0, 1.0],
                    &mut jvp,
                    cx,
                )
                .unwrap();
            assert_eq!(jvp, [0.0, 0.0]);
        }
    });
}

struct ObservedStorage {
    inner: QuadraticStorage,
    preparations: Rc<Cell<usize>>,
    gate: Rc<CancelGate>,
    cancel_gradient: Rc<Cell<bool>>,
    support: Rc<Cell<f64>>,
}

impl Storage for ObservedStorage {
    fn hamiltonian(&self, state: &[f64]) -> f64 {
        if state[0].abs() > self.support.get() {
            f64::NAN
        } else {
            self.inner.hamiltonian(state)
        }
    }

    fn gradient(&self, state: &[f64], output: &mut [f64]) {
        self.inner.gradient(state, output);
        if self.cancel_gradient.replace(false) {
            self.gate.request();
        }
    }

    fn quadratic_part(&self) -> Option<QuadraticPart> {
        self.preparations.set(self.preparations.get() + 1);
        self.inner.quadratic_part()
    }

    fn exact_quadratic_matrix(&self) -> Option<&[f64]> {
        self.inner.exact_quadratic_matrix()
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], midpoint: &[f64]) -> Option<f64> {
        self.inner.midpoint_energy_remainder(a, b, midpoint)
    }
}

fn observed_component(
    preparations: Rc<Cell<usize>>,
    gate: Rc<CancelGate>,
    cancel_gradient: Rc<Cell<bool>>,
    support: Rc<Cell<f64>>,
) -> ExactPortHamiltonianComponent {
    let owner = PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![2.0],
        vec![1.0],
        Box::new(ObservedStorage {
            inner: QuadraticStorage::new(vec![3.0], 1).unwrap(),
            preparations,
            gate,
            cancel_gradient,
            support,
        }),
    )
    .unwrap();
    ExactPortHamiltonianComponent::new(owner, vec![power("p", PortRole::EffortIn)], budget())
        .unwrap()
}

#[test]
fn duration_cache_is_bounded_and_uses_exact_duration_bits() {
    let preparations = Rc::new(Cell::new(0));
    let component = observed_component(
        Rc::clone(&preparations),
        Rc::new(CancelGate::new_clock_free()),
        Rc::new(Cell::new(false)),
        Rc::new(Cell::new(f64::INFINITY)),
    );
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let first = 0.1_f64;
        let adjacent = f64::from_bits(first.to_bits() + 1);
        let propose = |dt| {
            let (mut next, mut output) = (vec![17.0], [23.0, 29.0]);
            component
                .propose(&vec![0.8], dt, &[0.9, 300.0], &mut next, &mut output, cx)
                .unwrap();
            (bits(&next), bits(&output))
        };
        let expected = propose(first);
        let warm = preparations.get();
        assert_eq!(propose(first), expected);
        assert_eq!(preparations.get(), warm);
        propose(adjacent);
        assert!(
            preparations.get() > warm,
            "adjacent IEEE durations have distinct caches"
        );
        let two = preparations.get();
        assert_eq!(propose(first), expected);
        assert_eq!(preparations.get(), two);
        propose(0.2);
        let three = preparations.get();
        assert!(three > two);
        propose(adjacent);
        assert!(
            preparations.get() > three,
            "the least-recent of two durations was evicted"
        );
    });
}

#[test]
fn refused_or_cancelled_trials_and_tangents_publish_nothing_and_retry_cleanly() {
    let gate = Rc::new(CancelGate::new_clock_free());
    let armed = Rc::new(Cell::new(false));
    let support = Rc::new(Cell::new(1.0));
    let component = observed_component(
        Rc::new(Cell::new(0)),
        Rc::clone(&gate),
        Rc::clone(&armed),
        Rc::clone(&support),
    );
    let (old, input) = (vec![0.8], [0.9, 300.0]);
    let (mut expected, mut expected_output) = (Vec::new(), [0.0; 2]);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        component
            .propose(&old, 0.1, &input, &mut expected, &mut expected_output, cx)
            .unwrap();
        for bad_input in [[50.0, 300.0], [0.9, 0.0], [f64::NAN, 300.0]] {
            let (mut next, mut output) = (vec![17.0], [23.0, 29.0]);
            assert!(
                component
                    .propose(&old, 0.1, &bad_input, &mut next, &mut output, cx)
                    .is_err()
            );
            assert_eq!(next, vec![17.0]);
            assert_eq!(output, [23.0, 29.0]);
            assert!(
                component
                    .jvp(&old, 0.1, &bad_input, &[1.0, 0.0], &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [23.0, 29.0]);
            assert!(
                component
                    .vjp(&old, 0.1, &bad_input, &[0.4, 0.7], &mut output, cx)
                    .is_err()
            );
            assert_eq!(output, [23.0, 29.0]);
        }
    });
    armed.set(true);
    let (mut next, mut output) = (vec![17.0], [23.0, 29.0]);
    with_cx(&gate, |cx| {
        assert!(matches!(
            component.propose(&old, 0.1, &input, &mut next, &mut output, cx),
            Err(PhsComponentError::Owner(PreparedStepError::Cancelled))
        ));
        assert!(
            !armed.get(),
            "cancel request came from the original storage callback"
        );
        assert_eq!(next, vec![17.0]);
        assert_eq!(output, [23.0, 29.0]);
        assert!(matches!(
            component.jvp(&old, 0.1, &input, &[1.0, 0.0], &mut output, cx),
            Err(PhsComponentError::Owner(PreparedStepError::Cancelled))
        ));
        assert_eq!(output, [23.0, 29.0]);
        assert!(matches!(
            component.vjp(&old, 0.1, &input, &[1.0, 1.0], &mut output, cx),
            Err(PhsComponentError::Owner(PreparedStepError::Cancelled))
        ));
        assert_eq!(output, [23.0, 29.0]);
    });
    with_cx(&CancelGate::new_clock_free(), |cx| {
        component
            .propose(&old, 0.1, &input, &mut next, &mut output, cx)
            .unwrap();
        assert_eq!(bits(&next), bits(&expected));
        assert_eq!(bits(&output), bits(&expected_output));
    });
}

#[test]
fn adaptive_rejection_and_restore_keep_no_exact_trial_history() {
    let (mut system, owner, body) = heated_device();
    let (mut clean, clean_owner, clean_body) = heated_device();
    let initial = system.checkpoint().unwrap();
    let config = AdaptiveConfig {
        order: 2,
        min_step_s: 1.0e-9,
        max_attempts: 8,
    };
    with_cx(&CancelGate::new_clock_free(), |cx| {
        assert!(matches!(
            system.step_adaptive(0.125, config, |_, _, _| f64::NAN, cx),
            Err(AdaptiveError::InvalidEstimate { .. })
        ));
        assert_eq!(system.time_s().to_bits(), initial.time_s().to_bits());
        assert_eq!(system.ledger(), initial.ledger());
        assert_eq!(
            system.state::<Vec<f64>>(owner),
            initial.state::<Vec<f64>>(owner)
        );
        let mut estimates = 0;
        // Force one rejected comparison to exercise discarded duration caches.
        // An isolated exact held-input model has no time-discretization error.
        let report = system
            .step_adaptive(
                0.125,
                config,
                |_, _, _| {
                    estimates += 1;
                    if estimates == 1 { 24.0 } else { 0.0 }
                },
                cx,
            )
            .unwrap();
        assert_eq!(report.attempts, 2);
        assert!(report.dt_s < 0.125);
        assert_eq!(system.steps(), 2);
        clean.step(report.dt_s * 0.5, cx).unwrap();
        clean.step(report.dt_s * 0.5, cx).unwrap();
        assert_eq!(
            bits(system.state::<Vec<f64>>(owner).unwrap()),
            bits(clean.state::<Vec<f64>>(clean_owner).unwrap())
        );
        assert_eq!(
            system.state::<f64>(body).unwrap().to_bits(),
            clean.state::<f64>(clean_body).unwrap().to_bits()
        );
        assert_eq!(system.ledger(), clean.ledger());
        system.restore(&initial).unwrap();
        system.step(0.2, cx).unwrap();
        let (expected, _, heat) = scalar_oracle(0.8, 0.9, 0.2);
        close(
            system.state::<Vec<f64>>(owner).unwrap()[0],
            expected,
            2.0e-13,
        );
        close(
            system.ledger().nodes[owner.index()].dissipated_j.value(),
            heat,
            2.0e-12,
        );
    });
}

#[test]
fn metadata_and_nonquadratic_source_refusals_remain_native() {
    assert!(ExactPortHamiltonianComponent::new(scalar_owner(), vec![], budget()).is_err());
    assert!(
        ExactPortHamiltonianComponent::new(
            scalar_owner(),
            vec![power("heat", PortRole::EffortIn)],
            budget()
        )
        .is_err()
    );
    assert!(
        ExactPortHamiltonianComponent::new(
            coupled_owner(),
            vec![power("p", PortRole::EffortIn), power("p", PortRole::FlowIn)],
            budget()
        )
        .is_err()
    );
    assert!(
        ExactPortHamiltonianComponent::new(
            scalar_owner(),
            vec![PortSpec::signal_in(
                StableId::new("p").unwrap(),
                NonZeroUsize::MIN
            )],
            budget()
        )
        .is_err()
    );
    assert!(matches!(
        ExactPortHamiltonianComponent::new(
            fs_phs::duffing_oscillator(1.0, 2.0, 5.0, 0.3).unwrap(),
            vec![power("p", PortRole::EffortIn)],
            budget(),
        ),
        Err(PhsComponentError::Owner(PreparedStepError::Solver(_)))
    ));
    assert!(matches!(
        ExactPortHamiltonianComponent::new(
            coupled_owner(),
            vec![power("p", PortRole::EffortIn), power("q", PortRole::FlowIn)],
            LinearSplitBudget {
                max_modes: 1,
                max_omega_dt: 1.0
            },
        ),
        Err(PhsComponentError::Owner(PreparedStepError::Solver(_)))
    ));
}

struct ShiftedQuadratic {
    inner: QuadraticStorage,
    allowance_factor: Option<f64>,
}

impl Storage for ShiftedQuadratic {
    fn hamiltonian(&self, state: &[f64]) -> f64 {
        self.inner.hamiltonian(state) - 5.0e15
    }

    fn gradient(&self, state: &[f64], output: &mut [f64]) {
        self.inner.gradient(state, output);
    }

    fn quadratic_part(&self) -> Option<QuadraticPart> {
        self.inner.quadratic_part()
    }

    fn exact_quadratic_matrix(&self) -> Option<&[f64]> {
        self.inner.exact_quadratic_matrix()
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], midpoint: &[f64]) -> Option<f64> {
        self.inner.midpoint_energy_remainder(a, b, midpoint)
    }

    fn energy_evaluation_error_bound_j(&self, state: &[f64]) -> Option<f64> {
        self.allowance_factor
            .map(|factor| factor * (self.inner.hamiltonian(state).abs() + 5.0e15))
    }
}

fn shifted_component(allowance_factor: Option<f64>) -> ExactPortHamiltonianComponent {
    let source = ShiftedQuadratic {
        inner: QuadraticStorage::new(vec![1.0], 1).unwrap(),
        allowance_factor,
    };
    let owner =
        PortHamiltonian::new(1, 1, vec![0.0], vec![0.0], vec![1.0], Box::new(source)).unwrap();
    ExactPortHamiltonianComponent::new(owner, vec![power("p", PortRole::EffortIn)], budget())
        .unwrap()
}

#[test]
fn source_endpoint_evaluation_allowance_does_not_replace_native_physical_balance() {
    let component = shifted_component(Some(8.0 * f64::EPSILON));
    let (old, stale, input) = (vec![1.0e8], vec![17.0], [1.0, 300.0]);
    let case = Case {
        state: &old,
        stale: &stale,
        dt_s: 1.0,
        inputs: &input,
    };
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let Observed::Stepped { receipt, .. } = check_component(&component, &case, cx).unwrap()
        else {
            panic!("step");
        };
        close(receipt.stored_change_j, 100_000_000.5, 1.0e-6);
        assert_eq!(receipt.dissipated_j, 0.0);
        close(receipt.numerical_loss_j, 0.0, 1.0e-6);
        assert!(receipt.energy_tolerance_j > 0.5);
        // The reported stable native change differs by 0.5 J from subtracting
        // rounded absolute Hamiltonians near the source's large energy datum.
        assert!(check_component(&shifted_component(None), &case, cx).is_err());
        for invalid_bound in [-f64::EPSILON, f64::NAN, f64::INFINITY] {
            let bad = shifted_component(Some(invalid_bound));
            let (mut next, mut output) = (vec![17.0], [23.0, 29.0]);
            assert!(
                bad.propose(&old, 1.0, &input, &mut next, &mut output, cx)
                    .is_err()
            );
            assert_eq!(next, vec![17.0]);
            assert_eq!(output, [23.0, 29.0]);
        }
    });
}
