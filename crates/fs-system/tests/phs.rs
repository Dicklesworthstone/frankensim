//! Original multiport physics, physical heat, variable steps and atomic reuse.

#![allow(clippy::float_cmp)]

use core::cell::Cell;
use core::num::NonZeroUsize;
use std::rc::Rc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_phs::{
    PortHamiltonian, PreparedStepError, QuadraticStorage, SolveBudget, StepWorkspace, Storage,
};
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{Class, Component, PortRole, PortSpec, Proposal};
use fs_ports::{PortKind, StableId};
use fs_system::elements::HeatCapacity;
use fs_system::phs::{PhsComponentError, PortHamiltonianComponent};
use fs_system::{
    AdaptiveConfig, AdaptiveError, ErrorNorm, ErrorTolerance, NodeId, System, SystemBuilder,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 43,
                kernel_id: 19,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
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

#[test]
fn multiport_owner_conforms_and_keeps_native_work_heat_and_state_at_each_duration() {
    let component = PortHamiltonianComponent::new(
        coupled_owner(),
        vec![
            power("force", PortRole::EffortIn),
            power("motion", PortRole::FlowIn),
        ],
    )
    .unwrap();
    assert_eq!(component.class(), Class::Smooth);
    let mut native = StepWorkspace::new(component.system()).unwrap();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for (old, dt, input) in [
            (vec![0.4, -0.7], 0.1, [0.3, -0.2, 310.0]),
            (vec![-0.3, 0.2], 0.025, [-0.7, 0.6, 290.0]),
            (vec![0.4, -0.7], 0.1, [0.3, -0.2, 310.0]),
        ] {
            let Observed::Stepped {
                receipt, outputs, ..
            } = check_component(
                &component,
                &Case {
                    state: &old,
                    stale: &vec![31.0, -11.0],
                    dt_s: dt,
                    inputs: &input,
                },
                cx,
            )
            .unwrap()
            else {
                panic!("unexpected shortening")
            };
            let (mut expected, mut output) = (vec![0.0; 2], vec![0.0; 2]);
            let record = native
                .step_into(
                    component.system(),
                    &old,
                    &input[..2],
                    dt,
                    &mut expected,
                    &mut output,
                )
                .unwrap();
            let (mut next, mut actual) = (Vec::new(), [f64::NAN; 3]);
            component
                .propose(&old, dt, &input, &mut next, &mut actual, cx)
                .unwrap();
            assert_eq!(bits(&next), bits(&expected));
            assert_eq!(bits(&outputs[..2]), bits(&output));
            assert_eq!(bits(&actual), bits(&outputs));
            assert!(record.dissipated > 0.0);
            assert_eq!(receipt.dissipated_j.to_bits(), record.dissipated.to_bits());
            assert_eq!(receipt.stored_change_j.to_bits(), record.delta_h.to_bits());
            assert_eq!(
                receipt.numerical_loss_j.to_bits(),
                (record.supplied - record.dissipated - record.delta_h).to_bits()
            );
            assert_eq!(outputs[2].to_bits(), (-record.dissipated / dt).to_bits());
            assert_eq!(
                receipt.entropy_generated_j_per_k.to_bits(),
                (record.dissipated / input[2]).to_bits()
            );
        }
    });
}

#[test]
fn nonlinear_original_storage_remains_the_source_of_energy_and_motion() {
    let component = PortHamiltonianComponent::new(
        fs_phs::duffing_oscillator(1.0, 2.0, 5.0, 0.3).unwrap(),
        vec![power("p", PortRole::EffortIn)],
    )
    .unwrap();
    let old = vec![0.7, -0.2];
    close(
        component.stored_energy_j(&old).unwrap(),
        old[0].powi(2) + 1.25 * old[0].powi(4) + 0.5 * old[1].powi(2),
        1.0e-15,
    );
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for dt in [0.1, 0.05, 0.025] {
            check_component(
                &component,
                &Case {
                    state: &old,
                    stale: &vec![-1.0, 0.4],
                    dt_s: dt,
                    inputs: &[0.8, 300.0],
                },
                cx,
            )
            .unwrap();
            let original = fs_phs::step(component.system(), &old, &[0.8], dt).unwrap();
            let (mut next, mut output) = (old.clone(), [0.0; 2]);
            let Proposal::Step(receipt) = component
                .propose(&old, dt, &[0.8, 300.0], &mut next, &mut output, cx)
                .unwrap()
            else {
                panic!("step")
            };
            for (actual, expected) in next.iter().zip(&original.x) {
                close(*actual, *expected, 2.0e-13);
            }
            close(output[0], original.y[0], 2.0e-13);
            close(receipt.dissipated_j, original.dissipated, 2.0e-13);
        }
    });
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

fn heated_device() -> (System, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let component =
        PortHamiltonianComponent::new(scalar_owner(), vec![power("p", PortRole::EffortIn)])
            .unwrap();
    let owner = builder.add("owner", component, vec![0.8]).unwrap();
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

fn run_heated(dt: f64) -> f64 {
    let (mut system, owner, body) = heated_device();
    let (mut coordinate, mut heat, mut work) = (0.8, 0.0, 0.0);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for _ in 0..(0.2 / dt).round() as usize {
            // Independently solved midpoint equation: x' = -6 x + 0.9.
            let next = ((1.0 - 3.0 * dt) * coordinate + 0.9 * dt) / (1.0 + 3.0 * dt);
            let mean = 0.5 * (coordinate + next);
            heat += 18.0 * dt * mean * mean;
            work += 2.7 * dt * mean;
            coordinate = next;
            system.step(dt, cx).unwrap();
        }
    });
    let actual = system.state::<Vec<f64>>(owner).unwrap()[0];
    close(actual, coordinate, 2.0e-13);
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
    close(system.ledger().boundary_energy_in_j(), work, 2.0e-12);
    close(system.ledger().energy_residual_j(), 0.0, 2.0e-11);
    // Continuous physical oracle, independent of the discrete/native owner.
    (actual - (0.15 + 0.65 * (-1.2_f64).exp())).abs()
}

#[test]
fn composed_source_storage_and_thermal_body_match_physics_at_second_order() {
    let coarse = run_heated(0.01);
    let fine = run_heated(0.005);
    assert!(
        fine < 0.26 * coarse && fine > 0.24 * coarse,
        "{coarse} -> {fine}"
    );
}

#[test]
fn adaptive_half_steps_and_rollback_do_not_retain_native_trial_history() {
    let (mut system, owner, body) = heated_device();
    let (mut clean, clean_owner, clean_body) = heated_device();
    let initial = system.checkpoint().unwrap();
    let config = AdaptiveConfig {
        order: 2,
        min_step_s: 1.0e-9,
        max_attempts: 16,
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
        let report = system
            .step_adaptive(
                0.125,
                config,
                |a, b, c| {
                    let mut norm = ErrorNorm::new();
                    norm.include(
                        a.state::<Vec<f64>>(owner).unwrap()[0],
                        b.state::<Vec<f64>>(owner).unwrap()[0],
                        c.state::<Vec<f64>>(owner).unwrap()[0],
                        ErrorTolerance {
                            absolute: 1.0e-6,
                            relative: 0.0,
                        },
                    );
                    norm.value()
                },
                cx,
            )
            .unwrap();
        assert!(report.attempts > 1);
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
    });
}

struct CancellingStorage {
    gate: Rc<CancelGate>,
    armed: Rc<Cell<bool>>,
}

impl Storage for CancellingStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * x[0] * x[0]
    }
    fn gradient(&self, x: &[f64], output: &mut [f64]) {
        output[0] = x[0];
        if self.armed.replace(false) {
            self.gate.request();
        }
    }
}

#[test]
fn cancellation_inside_the_original_solver_preserves_outputs_and_allows_clean_retry() {
    let gate = Rc::new(CancelGate::new_clock_free());
    let armed = Rc::new(Cell::new(false));
    let owner = PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![0.2],
        vec![1.0],
        Box::new(CancellingStorage {
            gate: Rc::clone(&gate),
            armed: Rc::clone(&armed),
        }),
    )
    .unwrap();
    let component =
        PortHamiltonianComponent::new(owner, vec![power("p", PortRole::EffortIn)]).unwrap();
    let (old, input) = (vec![0.4], [0.7, 300.0]);
    let (mut expected, mut expected_output) = (old.clone(), [0.0; 2]);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        component
            .propose(&old, 0.1, &input, &mut expected, &mut expected_output, cx)
            .unwrap();
    });
    armed.set(true);
    let (mut next, mut output) = (vec![17.0], [23.0, 29.0]);
    with_cx(&gate, |cx| {
        assert!(matches!(
            component.propose(&old, 0.1, &input, &mut next, &mut output, cx),
            Err(PhsComponentError::Owner(PreparedStepError::Cancelled))
        ));
    });
    assert!(!armed.get());
    assert_eq!(next, vec![17.0]);
    assert_eq!(output, [23.0, 29.0]);
    with_cx(&CancelGate::new_clock_free(), |cx| {
        component
            .propose(&old, 0.1, &input, &mut next, &mut output, cx)
            .unwrap();
    });
    assert_eq!(bits(&next), bits(&expected));
    assert_eq!(bits(&output), bits(&expected_output));
}

#[test]
fn metadata_native_refusals_and_unphysical_energy_are_not_hidden() {
    assert!(PortHamiltonianComponent::new(scalar_owner(), vec![]).is_err());
    assert!(
        PortHamiltonianComponent::new(scalar_owner(), vec![power("heat", PortRole::EffortIn)])
            .is_err()
    );
    assert!(
        PortHamiltonianComponent::new(
            scalar_owner(),
            vec![PortSpec::signal_in(
                StableId::new("p").unwrap(),
                NonZeroUsize::MIN
            )]
        )
        .is_err()
    );
    assert!(
        PortHamiltonianComponent::new(
            coupled_owner(),
            vec![power("p", PortRole::EffortIn), power("p", PortRole::FlowIn)]
        )
        .is_err()
    );
    let mut component =
        PortHamiltonianComponent::new(scalar_owner(), vec![power("p", PortRole::EffortIn)])
            .unwrap();
    component
        .set_solve_budget(SolveBudget {
            relative: 1.0e-13,
            absolute: 0.0,
            max_iterations: 0,
        })
        .unwrap();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let (mut next, mut output) = (vec![17.0], [23.0, 29.0]);
        assert!(matches!(
            component.propose(&vec![0.8], 0.1, &[0.9, 300.0], &mut next, &mut output, cx),
            Err(PhsComponentError::Owner(PreparedStepError::Solver(_)))
        ));
        assert_eq!(next, vec![17.0]);
        assert_eq!(output, [23.0, 29.0]);
        assert!(
            component
                .propose(&vec![0.8], 0.1, &[0.9, 0.0], &mut next, &mut output, cx)
                .is_err()
        );
    });
    // Deliberately bypass the native structure gate: positive J creates
    // energy while the state equation itself solves perfectly. The adapter
    // must not launder that energy creation as signed numerical loss.
    let mutant = PortHamiltonian::from_raw_parts(
        1,
        1,
        vec![1.0],
        vec![0.0],
        vec![1.0],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
    );
    let component =
        PortHamiltonianComponent::new(mutant, vec![power("p", PortRole::EffortIn)]).unwrap();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let (mut next, mut output) = (vec![17.0], [23.0, 29.0]);
        assert!(matches!(
            component.propose(&vec![1.0], 0.1, &[0.0, 300.0], &mut next, &mut output, cx),
            Err(PhsComponentError::Balance { .. })
        ));
        assert_eq!(next, vec![17.0]);
        assert_eq!(output, [23.0, 29.0]);
    });
}
