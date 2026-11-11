//! The sliding interface under a prescribed slip rate, on the composition
//! contract (`fs_tribo::tangential::component`, feature `component`).
//!
//! Run with `cargo test -p fs-tribo --features component --test
//! sliding_interface`. The dry law is the crate's own; the spring is a
//! synthetic law for the storage terms, not a material.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{Case, Failure, Observed, check_component};
use fs_ports::component::{Class, Component, Proposal, Receipt};
use fs_tribo::tangential::component::{PrescribedSlip, SlidingInterface, SlidingState};
use fs_tribo::{
    DryTangentialState, FrictionLaw, FrictionRegime, TangentialLaw, TangentialSolveConfig,
    TangentialStep, TangentialStepInput, TriboError,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 1,
            kernel_id: 2,
            tile: 3,
            iteration: 4,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

const LOAD_N: f64 = 10.0;

/// A normal load of 10 N at `slip_rate_m_s`, with the budgets of the
/// friction-heat coupling's own tests.
fn drive(slip_rate_m_s: [f64; 2]) -> PrescribedSlip {
    PrescribedSlip {
        normal_n: LOAD_N,
        slip_rate_m_s,
        nominal_pressure_pa: Some(10.0),
        solve: TangentialSolveConfig {
            max_iterations: 64,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-12,
            force_relative_tolerance: 1e-12,
        },
        work_tolerance_j: 1e-12,
    }
}

/// A coefficient that falls with temperature: `0.8 - 0.01 (T - 300 K)`.
fn mu(kelvin: f64) -> f64 {
    0.8 - 0.01 * (kelvin - 300.0)
}

fn weakening(kelvin: f64) -> Result<FrictionLaw, TriboError> {
    Ok(FrictionLaw::Coulomb {
        static_mu: mu(kelvin),
        kinetic_mu: mu(kelvin),
    })
}

/// One step outside a kernel: the candidate state, the receipt and the
/// outputs (heat rate, then the two tractions).
fn propose<C: Component>(
    component: &C,
    state: &C::State,
    dt_s: f64,
    inputs: &[f64],
) -> Result<(C::State, Receipt, Vec<f64>), C::Error> {
    let mut next = state.clone();
    let mut outputs = vec![f64::NAN; 3];
    let proposal =
        with_cx(|cx| component.propose(state, dt_s, inputs, &mut next, &mut outputs, cx))?;
    let Proposal::Step(receipt) = proposal else {
        panic!("the interface takes the step or refuses");
    };
    Ok((next, receipt, outputs))
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} vs {b:.17e}");
}

fn field(result: Result<impl Sized, TriboError>) -> &'static str {
    match result {
        Err(TriboError::InvalidInput { field }) => field,
        Err(other) => panic!("{other:?}"),
        Ok(_) => panic!("admitted"),
    }
}

#[test]
fn the_interface_conforms_and_its_heat_is_the_laws_dissipation() {
    for rate in [[1.0, 0.0], [0.6, -0.8], [-0.25, 0.0]] {
        let speed = f64::hypot(rate[0], rate[1]);
        let interface = SlidingInterface::new(weakening, drive(rate)).expect("interface");
        let names: Vec<&str> = interface.ports().iter().map(|port| port.id().as_str()).collect();
        assert_eq!(names, ["heat", "traction"]);
        assert_eq!(interface.class(), Class::Smooth);
        let start = interface.start(300.0).expect("start");
        assert_eq!(interface.stored_energy_j(&start), Ok(0.0));
        // A state with history, to fill the candidate slot with.
        let (stale, _, _) = propose(&interface, &start, 0.5, &[310.0]).expect("stale");
        for kelvin in [300.0, 320.0, 350.0] {
            for dt_s in [0.2, 0.25, 1.0e-3] {
                let case = Case {
                    state: &start,
                    stale: &stale,
                    dt_s,
                    inputs: &[kelvin],
                };
                let observed = with_cx(|cx| check_component(&interface, &case, cx));
                let Ok(Observed::Stepped {
                    receipt, outputs, ..
                }) = observed
                else {
                    panic!("{observed:?}");
                };
                // Coulomb sliding: traction mu N along the slip, all of
                // its work exported as heat.
                let force_n = mu(kelvin) * LOAD_N;
                let heat_j = force_n * speed * dt_s;
                close(outputs[0], -force_n * speed, 4.0 * f64::EPSILON * force_n * speed);
                for axis in 0..2 {
                    close(outputs[1 + axis], force_n * rate[axis] / speed, 4.0 * f64::EPSILON * force_n);
                }
                close(receipt.dissipated_j, heat_j, 4.0 * f64::EPSILON * heat_j);
                close(receipt.stored_change_j, -heat_j, 4.0 * f64::EPSILON * heat_j);
                assert_eq!(receipt.numerical_loss_j, 0.0);
                assert_eq!(receipt.entropy_generated_j_per_k, receipt.dissipated_j / kelvin);
                assert_eq!((receipt.residual, receipt.tolerance), (0.0, 1e-12));

                let (next, same, _) = propose(&interface, &start, dt_s, &[kelvin]).expect("step");
                assert_eq!(same, receipt);
                // The law's own checkpoint: the slip it solved and the
                // temperature it was given.
                let law_state: &DryTangentialState = next.tangential();
                assert_eq!(law_state.slip_m(), [rate[0] * dt_s, rate[1] * dt_s]);
                assert_eq!(law_state.contact_temperature_k(), Some(kelvin));
                // The work the motion did is the account, and nothing else
                // is stored.
                assert_eq!(next.work_j(), receipt.dissipated_j);
                assert_eq!(interface.stored_energy_j(&next), Ok(-receipt.dissipated_j));
            }
        }
    }
}

#[test]
fn the_account_runs_on_and_an_interface_at_rest_makes_no_heat() {
    let interface = SlidingInterface::new(weakening, drive([1.0, 0.0])).expect("interface");
    let mut state = interface.start(300.0).expect("start");
    let mut worked_j = 0.0;
    for step in 0..50_u32 {
        let kelvin = 300.0 + f64::from(step);
        let (next, receipt, outputs) = propose(&interface, &state, 0.01, &[kelvin]).expect("step");
        worked_j += receipt.dissipated_j;
        close(outputs[0] * 0.01, -receipt.dissipated_j, 1e-15);
        state = next;
    }
    // 50 steps of 0.01 s at 1 m/s under mu N from 8 N down to 3.1 N.
    close(worked_j, 0.01 * (50.0 * 8.0 - 0.1 * (49.0 * 50.0 / 2.0)), 1e-12);
    close(state.work_j(), worked_j, 1e-13);
    close(state.tangential().slip_m()[0], 0.5, 1e-15);

    let rest = SlidingInterface::new(weakening, drive([0.0, 0.0])).expect("interface");
    let start = rest.start(300.0).expect("start");
    let (next, receipt, outputs) = propose(&rest, &start, 0.2, &[330.0]).expect("step");
    assert_eq!(outputs, [0.0, 0.0, 0.0]);
    assert_eq!((receipt.dissipated_j, receipt.stored_change_j), (0.0, 0.0));
    assert_eq!(receipt.entropy_generated_j_per_k, 0.0);
    assert_eq!(next.tangential().slip_m(), [0.0; 2]);
    assert_eq!(next.tangential().contact_temperature_k(), Some(330.0));
}

/// A linear tangential spring along the first axis: it stores
/// `k x^2 / 2` and its traction over a step is the mean of the two ends.
/// Two defects can be switched on: it dissipates `unaccounted_j` that no
/// work pays for, and it reports a stored change that is `misreported_j`
/// more than its stored energies changed.
#[derive(Debug, Clone, PartialEq)]
struct Spring {
    stiffness_n_m: f64,
    unaccounted_j: f64,
    misreported_j: f64,
}

/// The sound spring of `stiffness_n_m`.
fn spring(stiffness_n_m: f64) -> Spring {
    Spring {
        stiffness_n_m,
        unaccounted_j: 0.0,
        misreported_j: 0.0,
    }
}

impl TangentialLaw for Spring {
    type State = f64;

    fn model_id(&self) -> &str {
        "test-linear-spring"
    }

    fn initial_state(&self) -> f64 {
        0.0
    }

    fn stored_energy_j(&self, stretch_m: &f64) -> f64 {
        0.5 * self.stiffness_n_m * stretch_m * stretch_m
    }

    fn advance(
        &self,
        stretch_m: &f64,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<f64>, TriboError> {
        let moved_m = input.slip_increment_m[0];
        let next = stretch_m + moved_m;
        let traction_n = self.stiffness_n_m * (stretch_m + 0.5 * moved_m);
        let stored_change_j = self.stored_energy_j(&next) - self.stored_energy_j(stretch_m);
        Ok(TangentialStep {
            traction_n: [traction_n, 0.0],
            slip_increment_m: [moved_m, 0.0],
            state: next,
            regime: FrictionRegime::Sticking,
            stored_change_j: stored_change_j + self.misreported_j,
            dissipated_j: self.unaccounted_j,
            closure_residual_j: stored_change_j + self.unaccounted_j - traction_n * moved_m,
            compatibility_residual_m: [0.0; 2],
            iterations: 0,
        })
    }
}

#[test]
fn storage_is_in_the_account_and_a_law_that_revalues_the_held_state_is_refused() {
    let sound = |_: f64| Ok(spring(400.0));
    let interface = SlidingInterface::new(sound, drive([0.5, 0.0])).expect("interface");
    let start = interface.start(300.0).expect("start");
    let (stale, _, _) = propose(&interface, &start, 0.3, &[300.0]).expect("stale");
    let mut state = start.clone();
    for step in 0..8_u32 {
        let case = Case {
            state: &state,
            stale: &stale,
            dt_s: 0.02,
            inputs: &[300.0 + f64::from(step)],
        };
        let observed = with_cx(|cx| check_component(&interface, &case, cx));
        let Ok(Observed::Stepped {
            receipt, outputs, ..
        }) = observed
        else {
            panic!("{observed:?}");
        };
        // Nothing dissipated: no heat, and the stored change is what the
        // spring gained less the work that paid for it.
        assert_eq!(outputs[0], 0.0);
        assert_eq!(receipt.dissipated_j, 0.0);
        assert!(receipt.stored_change_j.abs() <= 1e-13, "{}", receipt.stored_change_j);
        assert_eq!(receipt.numerical_loss_j, -receipt.stored_change_j);
        let (next, _, _) = propose(&interface, &state, 0.02, case.inputs).expect("step");
        state = next;
    }
    // 8 steps of 0.01 m: 0.08 m of stretch, 1.28 J in the spring, paid for
    // by 1.28 J of work, so the component holds nothing net.
    close(*state.tangential(), 0.08, 1e-16);
    close(state.work_j(), 1.28, 1e-14);
    close(interface.stored_energy_j(&state).expect("stored"), 0.0, 1e-14);

    // The same spring with a stiffness that follows the temperature: at
    // another temperature the held stretch is worth another energy.
    let softening = |kelvin: f64| Ok(spring(400.0 - (kelvin - 300.0)));
    let revaluing = SlidingInterface::new(softening, drive([0.5, 0.0])).expect("interface");
    let cold = revaluing.start(300.0).expect("start");
    let (stretched, _, _) = propose(&revaluing, &cold, 0.02, &[300.0]).expect("step");
    assert!(propose(&revaluing, &stretched, 0.02, &[300.0]).is_ok());
    assert_eq!(
        field(propose(&revaluing, &stretched, 0.02, &[310.0])),
        "sliding_interface_law_changed_stored_energy"
    );
    // From the unloaded state every stiffness agrees, and the step is taken.
    assert!(propose(&revaluing, &cold, 0.02, &[310.0]).is_ok());
    // However little the law revalues the held state, it is the owner's
    // refusal. 1e-8 K changes the spring's 0.02 J by 5e-13 J, inside the
    // work tolerance of 1e-12 J and far outside the rounding the receipt
    // claims. The component as first landed (27fb5fd106) compared the two
    // energies to the work tolerance, took this step, and the contract's
    // check answered `EndpointEnergy { measured_change_j:
    // -5.000028169277471e-13, stored_change_j: 0.0, tolerance_j:
    // 3.694822225874806e-15 }` (found by review).
    let barely = [300.000_000_01];
    let case = Case {
        state: &stretched,
        stale: &stretched,
        dt_s: 0.02,
        inputs: &barely,
    };
    let checked = with_cx(|cx| check_component(&revaluing, &case, cx));
    assert!(
        matches!(
            checked,
            Err(Failure::Owner(TriboError::InvalidInput {
                field: "sliding_interface_law_changed_stored_energy"
            }))
        ),
        "{checked:?}"
    );
    assert_eq!(
        field(propose(&revaluing, &stretched, 0.02, &barely)),
        "sliding_interface_law_changed_stored_energy"
    );

}

#[test]
fn a_law_whose_accounts_do_not_balance_is_refused() {
    // A law whose heat no work pays for breaks the work balance.
    let leaking = |_: f64| {
        Ok(Spring {
            unaccounted_j: 1.0e-6,
            ..spring(400.0)
        })
    };
    let leaky = SlidingInterface::new(leaking, drive([0.5, 0.0])).expect("interface");
    let unloaded = leaky.start(300.0).expect("start");
    assert_eq!(
        field(propose(&leaky, &unloaded, 0.02, &[300.0])),
        "sliding_interface_work_balance"
    );
    // So does one whose energies balance while the stored change it
    // reports is not the change of its stored energies: nothing else is
    // wrong with this step.
    let misreporting = |_: f64| {
        Ok(Spring {
            misreported_j: 1.0e-6,
            ..spring(400.0)
        })
    };
    let misreports = SlidingInterface::new(misreporting, drive([0.5, 0.0])).expect("interface");
    let at_rest = misreports.start(300.0).expect("start");
    assert_eq!(
        field(propose(&misreports, &at_rest, 0.02, &[300.0])),
        "sliding_interface_work_balance"
    );
    // And one that reports heat taken out of nothing, or no number at all.
    for unaccounted_j in [-1.0e-6, f64::NAN] {
        let taking = move |_: f64| {
            Ok(Spring {
                unaccounted_j,
                ..spring(400.0)
            })
        };
        let takes = SlidingInterface::new(taking, drive([0.5, 0.0])).expect("interface");
        let rest = takes.start(300.0).expect("start");
        assert_eq!(
            field(propose(&takes, &rest, 0.02, &[300.0])),
            "sliding_interface_work_balance"
        );
    }
    // A law whose stored energy is not a number has no state to start from.
    let unbounded = SlidingInterface::new(|_: f64| Ok(spring(f64::INFINITY)), drive([0.5, 0.0]))
        .expect("interface");
    assert_eq!(field(unbounded.start(300.0)), "stored_energy_j");
    // The conformance check carries the owner's refusal as it is.
    let case = Case {
        state: &unloaded,
        stale: &unloaded,
        dt_s: 0.02,
        inputs: &[300.0],
    };
    assert!(matches!(
        with_cx(|cx| check_component(&leaky, &case, cx)),
        Err(Failure::Owner(TriboError::InvalidInput {
            field: "sliding_interface_work_balance"
        }))
    ));
}

#[test]
fn it_refuses_what_it_cannot_take() {
    let interface = SlidingInterface::new(weakening, drive([1.0, 0.0])).expect("interface");
    let start = interface.start(300.0).expect("start");
    for kelvin in [0.0, -5.0, f64::NAN, f64::INFINITY] {
        assert_eq!(field(propose(&interface, &start, 0.2, &[kelvin])), "interface_temperature_k");
        assert_eq!(field(interface.start(kelvin)), "interface_temperature_k");
    }
    for dt_s in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(field(propose(&interface, &start, dt_s, &[300.0])), "dt_s");
    }
    assert_eq!(field(propose(&interface, &start, 0.2, &[])), "sliding_interface_inputs");
    assert_eq!(
        field(propose(&interface, &start, 0.2, &[300.0, 300.0])),
        "sliding_interface_inputs"
    );

    // The source's refusal and the law's own are carried unchanged, and
    // the state they were asked from is as it was.
    let refusing = |kelvin: f64| {
        if kelvin > 400.0 {
            Err(TriboError::OutsideApplicability {
                field: "temperature_kelvin",
                value: kelvin,
                minimum: 250.0,
                maximum: 400.0,
            })
        } else {
            weakening(kelvin)
        }
    };
    let bounded = SlidingInterface::new(refusing, drive([1.0, 0.0])).expect("interface");
    let held = bounded.start(300.0).expect("start");
    let before = held.clone();
    assert_eq!(
        propose(&bounded, &held, 0.2, &[450.0]).err(),
        Some(TriboError::OutsideApplicability {
            field: "temperature_kelvin",
            value: 450.0,
            minimum: 250.0,
            maximum: 400.0,
        })
    );
    // At 390 K the coefficient is negative, which the dry law refuses.
    assert_eq!(field(propose(&bounded, &held, 0.2, &[390.0])), "static_mu");
    assert_eq!(held, before);

    // An output slice that is not the heat rate and the two tractions is
    // refused with nothing written to it or to the candidate.
    for length in [0_usize, 1, 2, 4] {
        let mut candidate = start.clone();
        let mut outputs = vec![f64::NAN; length];
        let answered = with_cx(|cx| {
            interface.propose(&start, 0.2, &[300.0], &mut candidate, &mut outputs, cx)
        });
        assert_eq!(field(answered), "sliding_interface_outputs");
        assert_eq!(candidate, start);
        assert!(outputs.iter().all(|value| value.is_nan()), "{outputs:?}");
    }

    // What an interface cannot be held to.
    let built = |change: fn(&mut PrescribedSlip)| {
        let mut held_to = drive([1.0, 0.0]);
        change(&mut held_to);
        SlidingInterface::new(weakening, held_to).map(|_| ())
    };
    assert_eq!(field(built(|d| d.normal_n = -1.0)), "normal_force");
    assert_eq!(field(built(|d| d.normal_n = f64::NAN)), "normal_force");
    assert_eq!(field(built(|d| d.nominal_pressure_pa = Some(-1.0))), "contact_pressure_pa");
    assert_eq!(field(built(|d| d.slip_rate_m_s = [0.0, f64::INFINITY])), "slip_rate_m_s");
    assert_eq!(field(built(|d| d.work_tolerance_j = -1e-12)), "work_tolerance_j");
    assert_eq!(field(built(|d| d.work_tolerance_j = f64::NAN)), "work_tolerance_j");
    assert_eq!(field(built(|d| d.solve.max_iterations = 0)), "max_iterations");
    assert!(built(|d| d.work_tolerance_j = 0.0).is_ok());
    assert!(built(|d| d.nominal_pressure_pa = None).is_ok());
}

/// The state type is nameable by a caller that stores it.
#[test]
fn a_state_is_what_it_holds() {
    let interface = SlidingInterface::new(weakening, drive([1.0, 0.0])).expect("interface");
    let start: SlidingState<DryTangentialState> = interface.start(300.0).expect("start");
    let (a, _, _) = propose(&interface, &start, 0.2, &[320.0]).expect("step");
    let (b, _, _) = propose(&interface, &start, 0.2, &[320.0]).expect("step");
    let (c, _, _) = propose(&interface, &start, 0.2, &[321.0]).expect("step");
    assert_eq!(a, b);
    assert_ne!(a, c);
    assert_ne!(a, start);
    assert_eq!(interface.drive(), &drive([1.0, 0.0]));
}
