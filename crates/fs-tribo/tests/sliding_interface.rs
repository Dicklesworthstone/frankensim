//! The sliding interface under a prescribed slip rate and the friction
//! interface between two bodies, on the composition contract
//! (`fs_tribo::tangential::component`, feature `component`).
//!
//! Run with `cargo test -p fs-tribo --features component --test
//! sliding_interface`. The dry law is the crate's own; the spring is a
//! synthetic law for the storage terms, not a material.

// Exact float comparisons here assert what the components fix exactly: a
// force that is a product of two given numbers, a pair at rest with no slip
// and no heat, a state that is another's to the bit.
#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{
    Case, Failure, Observed, check_component, check_load_line,
};
use fs_ports::component::{Class, Component, LoadLine, Proposal, Receipt};
use fs_tribo::stribeck_friction::HystereticFriction;
use fs_tribo::tangential::component::{
    FrictionInterface, FrictionState, InterfaceLoad, PrescribedSlip, SlidingInterface, SlidingState,
};
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

// The friction interface between two bodies.

/// Static and kinetic coefficients of [`dry`]: 5 N and 3 N under the load.
const HOLDS: f64 = 0.5;
const SLIDES: f64 = 0.3;

// A law source: its type is the interface's.
#[allow(clippy::unnecessary_wraps)]
fn dry(_: f64) -> Result<FrictionLaw, TriboError> {
    Ok(FrictionLaw::Coulomb {
        static_mu: HOLDS,
        kinetic_mu: SLIDES,
    })
}

/// The load of 10 N with the budgets of [`drive`] and an impulse tolerance
/// of 1e-10 N s.
fn pressed() -> InterfaceLoad {
    let held_to = drive([0.0; 2]);
    InterfaceLoad {
        normal_n: LOAD_N,
        nominal_pressure_pa: held_to.nominal_pressure_pa,
        solve: held_to.solve,
        impulse_tolerance_n_s: 1e-10,
        work_tolerance_j: 1e-12,
    }
}

type LawAt = fn(f64) -> Result<FrictionLaw, TriboError>;
type Rubbing = FrictionInterface<FrictionLaw, LawAt>;
type Rubbed = FrictionState<DryTangentialState>;

fn rubbing(law_at: LawAt) -> Rubbing {
    FrictionInterface::new(law_at, pressed()).expect("interface")
}

/// What the interface answered to one proposal on a line.
#[derive(Debug)]
enum Answer {
    /// It took the step: the candidate, the receipt, the force at `a` [N]
    /// and the heat rate into the interface [W].
    Took {
        state: Rubbed,
        receipt: Receipt,
        force_n: f64,
        heat_w: f64,
    },
    /// It asked for a step of this length [s].
    Cut(f64),
}

fn on_line(
    interface: &Rubbing,
    state: &Rubbed,
    dt_s: f64,
    inputs: &[f64],
    compliance: &[f64],
) -> Result<Answer, TriboError> {
    let mut next = state.clone();
    let mut outputs = vec![f64::NAN; 3];
    let proposal = with_cx(|cx| {
        interface.propose_on_line(state, dt_s, inputs, compliance, &mut next, &mut outputs, cx)
    })?;
    Ok(match proposal {
        Proposal::Step(receipt) => {
            assert_eq!(outputs[1], -outputs[0]);
            Answer::Took {
                state: next,
                receipt,
                force_n: outputs[0],
                heat_w: outputs[2],
            }
        }
        Proposal::Shorten { dt_s } => Answer::Cut(dt_s),
    })
}

/// One proposal with port `a` on a body of `mass_kg` that a force of
/// `push_n` drives along the slip, and `b` on nothing that answers: on the
/// midpoint rule the line is `c = dt / (2 m)` with the released mean slip
/// `start + c push`.
fn on_a_body(interface: &Rubbing, state: &Rubbed, dt_s: f64, mass_kg: f64, push_n: f64) -> Answer {
    let c = 0.5 * dt_s / mass_kg;
    let released = state.slip_m_s() + c * push_n;
    on_line(interface, state, dt_s, &[released, 0.0, 300.0], &[c, 0.0, 0.0, 0.0]).expect("answer")
}

#[test]
fn at_rest_the_laws_own_map_holds_up_to_the_static_limit_and_slides_beyond_it() {
    let interface = rubbing(dry);
    let names: Vec<&str> = interface.ports().iter().map(|port| port.id().as_str()).collect();
    assert_eq!(names, ["a", "b", "heat"]);
    assert_eq!(interface.class(), Class::NonSmooth);
    assert_eq!(interface.line_ports(), [0, 1]);
    assert_eq!(interface.load(), &pressed());
    let rest = interface.start(300.0, 0.0).expect("start");
    assert!(rest.at_rest());
    assert_eq!(interface.stored_energy_j(&rest), Ok(0.0));
    let (dt_s, mass_kg) = (0.01, 2.0);
    // Up to 5 N the pair holds with the force the push needs, no slip and
    // no heat.
    for push_n in [0.0, 1.0, -4.0, 4.999, -4.999] {
        let Answer::Took {
            state,
            receipt,
            force_n,
            heat_w,
        } = on_a_body(&interface, &rest, dt_s, mass_kg, push_n)
        else {
            panic!("{push_n} N cut the step");
        };
        close(force_n, push_n, 8.0 * f64::EPSILON * 5.0);
        assert!(state.at_rest(), "{push_n} N");
        assert_eq!((heat_w, receipt.dissipated_j, receipt.stored_change_j), (0.0, 0.0, 0.0));
        assert_eq!(state.tangential().slip_m(), [0.0; 2]);
    }
    // Beyond it the pair slides against 3 N, from rest on the midpoint
    // rule: the mean slip is `c (push - 3)` and the step ends with twice
    // that.
    let c = 0.5 * dt_s / mass_kg;
    for push_n in [5.5_f64, -8.0, 50.0] {
        let Answer::Took {
            state,
            receipt,
            force_n,
            heat_w,
        } = on_a_body(&interface, &rest, dt_s, mass_kg, push_n)
        else {
            panic!("{push_n} N cut the step");
        };
        let way = push_n.signum();
        assert_eq!(force_n, SLIDES * LOAD_N * way);
        let mean = c * (push_n - 3.0 * way);
        close(state.slip_m_s(), 2.0 * mean, 8.0 * f64::EPSILON * mean.abs());
        close(heat_w, -3.0 * mean.abs(), 8.0 * f64::EPSILON * 3.0 * mean.abs());
        close(receipt.dissipated_j, dt_s * 3.0 * mean.abs(), 8.0 * f64::EPSILON * mean.abs());
        assert_eq!(receipt.entropy_generated_j_per_k, receipt.dissipated_j / 300.0);
        close(state.tangential().slip_m()[0], dt_s * mean, 8.0 * f64::EPSILON * mean.abs());
        assert_eq!(state.tangential().contact_temperature_k(), Some(300.0));
    }
}

#[test]
fn a_sliding_pair_slides_on_or_is_cut_at_the_stop_and_the_cut_step_ends_at_rest() {
    let interface = rubbing(dry);
    let mass_kg = 2.0;
    for start in [1.0_f64, -0.4] {
        let way = start.signum();
        let sliding = interface.start(300.0, start).expect("start");
        assert!(!sliding.at_rest());
        // Pushed along the slip, coasting, or braked too little to stop
        // inside the step: it slides on against 3 N, and the body's slip
        // changes by `dt (push - 3) / m`.
        for (dt_s, push_n) in [(0.01, 5.0 * way), (0.01, 0.0), (0.05, -2.0 * way)] {
            let Answer::Took {
                state,
                receipt,
                force_n,
                heat_w,
            } = on_a_body(&interface, &sliding, dt_s, mass_kg, push_n)
            else {
                panic!("{start} m/s, {push_n} N cut the step");
            };
            assert_eq!(force_n, SLIDES * LOAD_N * way);
            let ends = start + dt_s * (push_n - 3.0 * way) / mass_kg;
            close(state.slip_m_s(), ends, 8.0 * f64::EPSILON);
            let mean = 0.5 * (start + ends);
            close(heat_w, -3.0 * mean.abs(), 16.0 * f64::EPSILON * 3.0);
            close(receipt.dissipated_j, dt_s * 3.0 * mean.abs(), 16.0 * f64::EPSILON * 3.0 * dt_s);
            assert_eq!(receipt.numerical_loss_j, 0.0);
        }
        // Against 3 N and a push the slip stops after
        // `m abs(start) / (3 - way push)`. Asked for any step that ends
        // beyond that, the interface asks for the step that ends there, the
        // same from every step length. That step ends at rest, exactly,
        // with the heat of the stopping distance. A step that ends short
        // of it slides on.
        for push_n in [0.0, 2.0 * way, -7.0 * way] {
            let stop_s = mass_kg * start.abs() / (3.0 - way * push_n);
            let mut found_s = f64::NAN;
            for beyond in [1.0001, 1.5, 4.0, 1000.0] {
                let Answer::Cut(asked_s) =
                    on_a_body(&interface, &sliding, beyond * stop_s, mass_kg, push_n)
                else {
                    panic!("{start} m/s, {push_n} N, {beyond}: not cut");
                };
                close(asked_s, stop_s, 64.0 * f64::EPSILON * stop_s);
                found_s = asked_s;
            }
            let Answer::Took {
                state,
                receipt,
                force_n,
                heat_w,
            } = on_a_body(&interface, &sliding, found_s, mass_kg, push_n)
            else {
                panic!("{start} m/s, {push_n} N: the located step is cut again");
            };
            // At rest, under the law's own kinetic traction to the last
            // digit; what the stop left of the slip is rounding here.
            assert!(state.at_rest());
            assert!(state.slip_m_s().abs() <= 4.0 * f64::EPSILON, "{}", state.slip_m_s());
            assert_eq!(force_n, SLIDES * LOAD_N * way);
            close(heat_w, -1.5 * start.abs(), 1e-9);
            close(receipt.dissipated_j, stop_s * 1.5 * start.abs(), 1e-9);
            close(state.tangential().slip_m()[0], 0.5 * stop_s * start, 1e-15);
            let Answer::Took { state, .. } =
                on_a_body(&interface, &sliding, 0.999 * stop_s, mass_kg, push_n)
            else {
                panic!("{start} m/s, {push_n} N: a step short of the stop is cut");
            };
            assert!(!state.at_rest());
            assert!(state.slip_m_s() * way > 0.0, "{}", state.slip_m_s());
            close(state.slip_m_s(), 1.0e-3 * start, 1.0e-12);
        }
    }
}

#[test]
fn a_slip_below_the_impulse_tolerance_is_rest() {
    // On a body of 2 kg a slip of 5e-11 m/s is a momentum of 1e-10 N s,
    // the tolerance. Below it the pair is at rest and the static map
    // decides: a push of 1 N is held. Above it the pair slides, 3 N stop
    // it almost at once, and the step is cut there; the cut step ends at
    // rest.
    let interface = rubbing(dry);
    let Answer::Took { state, force_n, .. } =
        on_a_body(&interface, &interface.start(300.0, 4.0e-11).expect("start"), 0.01, 2.0, 1.0)
    else {
        panic!("a slip below the tolerance cut the step");
    };
    assert!(state.at_rest());
    close(force_n, 1.0, 1e-7);
    let creeping = interface.start(300.0, 6.0e-11).expect("start");
    let Answer::Cut(asked_s) = on_a_body(&interface, &creeping, 0.01, 2.0, 1.0) else {
        panic!("a slip above the tolerance was not cut");
    };
    close(asked_s, 2.0 * 6.0e-11 / (3.0 - 1.0), 1e-20);
    let Answer::Took { state, .. } = on_a_body(&interface, &creeping, asked_s, 2.0, 1.0) else {
        panic!("the located step is cut again");
    };
    assert!(state.at_rest());

    // The same tolerance is the band of a stop. A slip of 1 m/s on the
    // body of 2 kg stops after 2/3 s against 3 N. A step longer than that
    // by 1e-11 of itself ends beyond the stop by 2e-11 N s of momentum:
    // inside the tolerance, outside the law's force tolerance of 4e-12 N
    // alone. It is taken to end at rest under the law's own 3 N, to the
    // last digit, and what is left of the slip, 1e-11 m/s the other way,
    // is kept.
    let sliding = interface.start(300.0, 1.0).expect("start");
    let dt_s = 2.0 / 3.0 * (1.0 + 1.0e-11);
    let Answer::Took {
        state,
        receipt,
        force_n,
        heat_w,
    } = on_a_body(&interface, &sliding, dt_s, 2.0, 0.0)
    else {
        panic!("a step that ends 2e-11 N s beyond the stop was cut");
    };
    assert!(state.at_rest());
    assert_eq!(force_n, SLIDES * LOAD_N);
    close(state.slip_m_s(), -1.0e-11, 1.0e-15);
    let mean = 0.5 * (1.0 + state.slip_m_s());
    close(heat_w, -force_n * mean, 2.0 * f64::EPSILON);
    close(receipt.dissipated_j, dt_s * force_n * mean, 2.0 * f64::EPSILON);
    // The law's own accounts close for the step it was asked for.
    assert_eq!((receipt.residual, receipt.numerical_loss_j), (0.0, 0.0));
    // The first step at rest takes the rest of the slip out: the pair is
    // held with the force that ends that step at none, the push of 1 N and
    // half the leftover over `c`, and the record is then exactly zero.
    let (later_s, c) = (0.01_f64, 0.25 * 0.01);
    let Answer::Took {
        state: held,
        force_n,
        ..
    } = on_a_body(&interface, &state, later_s, 2.0, 1.0)
    else {
        panic!("the step after the stop was cut");
    };
    assert!(held.at_rest());
    assert_eq!(held.slip_m_s(), 0.0);
    close(force_n - 1.0, 0.5 * state.slip_m_s() / c, 1.0e-12);
    // The band of a stop is the impulse tolerance and the law's force
    // tolerance over the step, so a stop can leave a hair more than the
    // tolerance alone: here 1.015e-10 N s. The pair is at rest by the
    // state's mark, not by measuring that slip again, and the next step
    // holds it.
    let dt_s = 2.0 / 3.0 * (1.0 + 5.075e-11);
    let Answer::Took { state, .. } = on_a_body(&interface, &sliding, dt_s, 2.0, 0.0) else {
        panic!("a step that ends 1.015e-10 N s beyond the stop was cut");
    };
    assert!(state.at_rest());
    assert!(2.0 * state.slip_m_s().abs() > 1.0e-10, "{}", state.slip_m_s());
    let Answer::Took { state: held, .. } = on_a_body(&interface, &state, later_s, 2.0, 1.0) else {
        panic!("the step after that stop was cut");
    };
    assert!(held.at_rest());
    assert_eq!(held.slip_m_s(), 0.0);
    // And the step is the kernel's to audit: the work of the two ends is
    // the heat, at the flows on the line.
    let c = 0.25 * dt_s;
    let case = Case {
        state: &sliding,
        stale: &state,
        dt_s,
        inputs: &[1.0, 0.0, 300.0],
    };
    let audited = with_cx(|cx| check_load_line(&interface, &case, &[c, 0.0, 0.0, 0.0], cx));
    assert!(matches!(&audited, Ok(on_line) if matches!(on_line.observed, Observed::Stepped { .. })), "{audited:?}");
    // A step twice as far beyond the stop as the tolerance allows is cut.
    let Answer::Cut(asked_s) = on_a_body(&interface, &sliding, 2.0 / 3.0 * (1.0 + 2.0e-10), 2.0, 0.0)
    else {
        panic!("a step that ends 4e-10 N s beyond the stop was taken");
    };
    close(asked_s, 2.0 / 3.0, 1.0e-15);
    assert_eq!(field(FrictionInterface::new(dry as LawAt, {
        let mut load = pressed();
        load.impulse_tolerance_n_s = 0.0;
        load
    })), "impulse_tolerance_n_s");
}

#[test]
fn a_load_that_reaches_the_static_limit_inside_a_step_is_cut_there() {
    // At rest the law's static map decides for a whole step from the
    // step's mean load. Where the load reaches what the law holds inside
    // the step, the interface asks for the step that ends there first.
    // The load at a step's two ends is read off the straight line through
    // the means of the last step and of this one.
    let interface = rubbing(dry);
    let (dt_s, mass_kg) = (0.01_f64, 2.0);
    let rest = interface.start(300.0, 0.0).expect("start");
    let held_under = |state: &Rubbed, step_s: f64, push_n: f64| -> Rubbed {
        let Answer::Took { state, force_n, .. } = on_a_body(&interface, state, step_s, mass_kg, push_n)
        else {
            panic!("{push_n} N over {step_s} s was cut");
        };
        assert!(state.at_rest(), "{push_n} N over {step_s} s");
        close(force_n, push_n, 1e-12);
        state
    };
    // A run's first step has no step before it: by its mean, as before.
    let first = held_under(&rest, dt_s, 4.0);

    // From 4 N to a mean of 4.9 N over two steps of 0.01 s the load rises
    // at 90 N/s: 4.45 N at this step's start, 5.35 N at its end, and the
    // 5 N held are reached after 0.55 / 0.9 of it.
    let fraction = (5.0 - 4.45) / 0.9;
    let Answer::Cut(asked_s) = on_a_body(&interface, &first, dt_s, mass_kg, 4.9) else {
        panic!("a load that passes the limit inside the step was not cut");
    };
    close(asked_s, fraction * dt_s, 1e-15);
    // The step of that length has the mean of its own stretch of the
    // line: it is held, and ends at the limit.
    let at_the_limit = held_under(&first, asked_s, 4.45 + 0.5 * fraction * 0.9);
    // The load rises on: the next step starts at the limit, and the law's
    // map lets go from its start, against 3 N.
    let Answer::Took { state, force_n, .. } =
        on_a_body(&interface, &at_the_limit, dt_s, mass_kg, 5.0 + 90.0 * 0.5 * dt_s)
    else {
        panic!("the step at the limit was cut again");
    };
    assert!(!state.at_rest());
    assert_eq!(force_n, SLIDES * LOAD_N);

    // The last step's length is in the line: after 0.004 s at 4 N, a mean
    // of 4.9 N over 0.01 s is a rise of 0.9 N in 0.007 s.
    let short = held_under(&rest, 0.004, 4.0);
    let rate = 0.9 / 0.007;
    let (from, to) = (4.9 - 0.005 * rate, 4.9 + 0.005 * rate);
    let Answer::Cut(asked_s) = on_a_body(&interface, &short, dt_s, mass_kg, 4.9) else {
        panic!("not cut after a shorter step");
    };
    close(asked_s, (5.0 - from) / (to - from) * dt_s, 1e-15);

    // Inside the limit at both ends: held whole. Beyond it from the
    // step's start: the law's map on the mean, which lets go.
    held_under(&first, dt_s, 4.4);
    let Answer::Took { state, .. } = on_a_body(&interface, &first, dt_s, mass_kg, 7.0) else {
        panic!("a load beyond the limit from the start was cut");
    };
    assert!(!state.at_rest());
    // Falling through the limit the other way is the same event with the
    // other sign.
    let pulled = held_under(&rest, dt_s, -4.0);
    let Answer::Cut(asked_s) = on_a_body(&interface, &pulled, dt_s, mass_kg, -4.9) else {
        panic!("a load that passes minus the limit was not cut");
    };
    close(asked_s, fraction * dt_s, 1e-15);

    // An instant within the impulse tolerance of a step's end is that
    // end. The load passes the limit at a constant rate, so letting go a
    // time t late is a momentum of half the rate times t squared: with
    // 1e-10 N s and 90 N/s, 1.5e-6 s. A line that ends 5e-5 N beyond the
    // limit reaches it 5.6e-7 s before the end and is held whole; one
    // that ends 1e-3 N beyond, 1.1e-5 s before the end, is cut.
    let ending = |beyond_n: f64| on_a_body(&interface, &first, dt_s, mass_kg, (7.0 + beyond_n) / 1.5);
    assert!(matches!(ending(5.0e-5), Answer::Took { ref state, .. } if state.at_rest()));
    assert!(matches!(ending(1.0e-3), Answer::Cut(_)));
}

#[test]
fn a_given_end_that_jumps_starts_the_step_with_that_much_more_slip() {
    // `a` on a body of 2 kg, `b` prescribed: its row and column of the
    // line are zero. Both move at 0.2 m/s and the pair is at rest. The
    // first time an end is given a flow nothing is inferred.
    let interface = rubbing(dry);
    let (dt_s, mass_kg) = (0.01_f64, 2.0);
    let line = [0.5 * dt_s / mass_kg, 0.0, 0.0, 0.0];
    let rest = interface.start(300.0, 0.0).expect("start");
    let Ok(Answer::Took { state, force_n, .. }) =
        on_line(&interface, &rest, dt_s, &[0.2, 0.2, 300.0], &line)
    else {
        panic!("not taken");
    };
    assert!(state.at_rest());
    assert_eq!(force_n, 0.0);
    // `b` jumps to 0.5 m/s: the step starts with a slip of -0.3 m/s, and
    // the body is dragged after the side under 3 N, 1.5 m/s2.
    let Ok(Answer::Took {
        state: dragged,
        force_n,
        ..
    }) = on_line(&interface, &state, dt_s, &[0.2, 0.5, 300.0], &line)
    else {
        panic!("not taken");
    };
    assert_eq!(force_n, -SLIDES * LOAD_N);
    close(dragged.slip_m_s(), -0.3 + 0.015, 4.0 * f64::EPSILON);
    // Held at 0.5 m/s there is no further jump: the slip goes on from
    // where the last step left it. The body's released flow is its own
    // velocity, 0.2 m/s and what the first step added.
    let Ok(Answer::Took { state: next, .. }) =
        on_line(&interface, &dragged, dt_s, &[0.215, 0.5, 300.0], &line)
    else {
        panic!("not taken");
    };
    close(next.slip_m_s(), -0.3 + 0.030, 8.0 * f64::EPSILON);
    // Had the side not jumped, the pair would have stayed at rest.
    let Ok(Answer::Took { state: kept, .. }) =
        on_line(&interface, &state, dt_s, &[0.2, 0.2, 300.0], &line)
    else {
        panic!("not taken");
    };
    assert!(kept.at_rest());
    // A jump ends the rest, whatever its size: after one of 0.01 m/s,
    // which 2 N of the 5 N held could make up within the step, the pair
    // slides against 3 N, and the step is cut where the side is caught,
    // after two thirds of it.
    let Ok(Answer::Cut(asked_s)) = on_line(&interface, &state, dt_s, &[0.2, 0.21, 300.0], &line)
    else {
        panic!("a small jump was held");
    };
    close(asked_s, dt_s * 2.0 / 3.0, 1.0e-15);
}

#[test]
fn with_both_ends_given_it_is_the_prescribed_interface() {
    // No compliance: the mandatory tier. The slip is what the two flows
    // make it, the law is asked for its traction at that slip, and the
    // force and the law's checkpoint are those of `SlidingInterface` held
    // to the same slip rate; the heat is the same to rounding.
    let interface = rubbing(weakening);
    let rest = interface.start(300.0, 0.0).expect("start");
    for (a, b) in [(0.25, 0.0), (0.1, 1.1), (-0.5, -0.5)] {
        for (kelvin, dt_s) in [(300.0, 0.2), (320.0, 1.0e-3)] {
            let (next, receipt, outputs) =
                propose(&interface, &rest, dt_s, &[a, b, kelvin]).expect("step");
            let prescribed = SlidingInterface::new(weakening, drive([a - b, 0.0])).expect("held");
            let (held, booked, expected) =
                propose(&prescribed, &prescribed.start(300.0).expect("start"), dt_s, &[kelvin])
                    .expect("step");
            assert_eq!(outputs[0], expected[1], "{a} {b} {kelvin}");
            assert_eq!(outputs[1], -expected[1]);
            close(outputs[2], expected[0], 4.0 * f64::EPSILON * expected[0].abs());
            close(receipt.dissipated_j, booked.dissipated_j, 4.0 * f64::EPSILON * booked.dissipated_j);
            assert_eq!(next.tangential(), held.tangential());
            assert_eq!(next.slip_m_s(), a - b);
            assert_eq!(receipt.stored_change_j, 0.0);
        }
    }
}

#[test]
fn the_friction_interface_conforms_on_its_line_and_its_heat_is_the_work_of_its_ends() {
    let interface = rubbing(weakening);
    let stale = {
        let rest = interface.start(300.0, 0.0).expect("start");
        let Ok(Answer::Took { state, .. }) =
            on_line(&interface, &rest, 0.01, &[0.4, 0.1, 310.0], &[0.0; 4])
        else {
            panic!("stale");
        };
        state
    };
    // One end on a body, two bodies, and two bodies that answer each
    // other's force.
    let lines = [
        [2.5e-3, 0.0, 0.0, 0.0],
        [2.5e-3, 0.0, 0.0, 1.0e-2],
        [2.5e-3, -1.0e-3, -1.0e-3, 1.0e-2],
        [0.0; 4],
    ];
    let (mut stepped, mut cut) = (0, 0);
    for start in [0.0, 0.3, -0.3, 1.0e-3] {
        let state = interface.start(300.0, start).expect("start");
        for line in lines {
            for (a, b) in [(0.0, 0.0), (0.01, 0.0), (0.3, 0.0), (0.0, 0.3), (-0.2, 0.2), (0.31, 0.02)] {
                for kelvin in [300.0, 330.0] {
                    let case = Case {
                        state: &state,
                        stale: &stale,
                        dt_s: 0.01,
                        inputs: &[a, b, kelvin],
                    };
                    let on_line = with_cx(|cx| check_load_line(&interface, &case, &line, cx))
                        .unwrap_or_else(|failure| {
                            panic!("{start} {line:?} {a} {b} {kelvin}: {failure:?}")
                        });
                    match on_line.observed {
                        Observed::Shortened { .. } => cut += 1,
                        Observed::Stepped {
                            receipt, outputs, ..
                        } => {
                            stepped += 1;
                            // The heat is the force times the slip that
                            // crossed, and all of it is booked.
                            let slip = on_line.flows[0] - on_line.flows[1];
                            let work_w = outputs[0] * slip;
                            close(outputs[2], -work_w.max(0.0), 1e-14);
                            close(receipt.dissipated_j, 0.01 * work_w.max(0.0), 1e-16);
                            assert!(outputs[0].abs() <= mu(kelvin) * LOAD_N * (1.0 + 1e-9) + 1e-7);
                            assert!(work_w >= -1e-14, "{work_w}");
                        }
                    }
                }
            }
        }
    }
    assert!(stepped > 100 && cut > 10, "{stepped} stepped, {cut} cut");
    println!("the friction interface on its line: {stepped} points stepped and audited, {cut} cut");
}

/// A law with every verdict on record: the free slip rate it was asked at
/// [m/s] and what it says the answer cost it (its own count of iterations),
/// or `None` where it refused.
struct Counted<L>(L);

thread_local! {
    static VERDICTS: core::cell::RefCell<Vec<(f64, Option<usize>)>> =
        const { core::cell::RefCell::new(Vec::new()) };
}

impl<L: TangentialLaw> TangentialLaw for Counted<L> {
    type State = L::State;

    fn model_id(&self) -> &str {
        self.0.model_id()
    }

    fn initial_state(&self) -> Self::State {
        self.0.initial_state()
    }

    fn stored_energy_j(&self, state: &Self::State) -> f64 {
        self.0.stored_energy_j(state)
    }

    fn requires_physical_scalar_compliance(&self) -> bool {
        self.0.requires_physical_scalar_compliance()
    }

    fn advance(
        &self,
        state: &Self::State,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        let result = self.0.advance(state, input);
        let cost = result.as_ref().ok().map(|step| step.iterations);
        VERDICTS.with(|verdicts| {
            verdicts.borrow_mut().push((input.slip_increment_m[0] / input.dt_s, cost));
        });
        result
    }
}

type CountedAt<L> = fn(f64) -> Result<Counted<L>, TriboError>;

fn counted_dry(kelvin: f64) -> Result<Counted<FrictionLaw>, TriboError> {
    dry(kelvin).map(Counted)
}

/// The friction of the bowed-string fixture, with the branch memory of the
/// module: static 1.2, kinetic 0.35 with a Gaussian drop over 0.05 m/s. Its
/// verdict on a line is a search: it prepares the line, then finds a root.
fn searching(_: f64) -> Result<Counted<HystereticFriction>, TriboError> {
    HystereticFriction::new(FrictionLaw::Stribeck {
        static_mu: 1.2,
        kinetic_mu: 0.35,
        characteristic_speed: 0.05,
        viscous_per_speed: 0.0,
    })
    .map(Counted)
}

/// What an interface answered to one proposal: the step it asks for
/// instead, if it cuts; its candidate; the force it takes the step with.
type Asked<S> = (Option<f64>, FrictionState<S>, f64);

/// One proposal of a counted law's interface from `state`, with port `a` on
/// a line of compliance `c` under a load of `load_n` and `b` on nothing
/// that answers.
fn asked_under<L: TangentialLaw>(
    interface: &FrictionInterface<Counted<L>, CountedAt<L>>,
    state: &FrictionState<L::State>,
    dt_s: f64,
    c: f64,
    load_n: f64,
) -> Result<Asked<L::State>, TriboError> {
    let released = state.slip_m_s() + c * load_n;
    let mut next = state.clone();
    let mut outputs = vec![f64::NAN; 3];
    let proposal = with_cx(|cx| {
        interface.propose_on_line(
            state,
            dt_s,
            &[released, 0.0, 300.0],
            &[c, 0.0, 0.0, 0.0],
            &mut next,
            &mut outputs,
            cx,
        )
    })?;
    Ok(match proposal {
        Proposal::Step(_) => (None, next, outputs[0]),
        Proposal::Shorten { dt_s } => (Some(dt_s), next, f64::NAN),
    })
}

#[test]
fn the_instant_of_a_break_away_is_narrowed_as_far_as_the_tolerance_needs_of_a_law_that_searches() {
    // The step is cut where the law's verdict changes, and the verdict is
    // asked of the law at loads between the step's two ends. How near the
    // instant it is asked is the question here.
    //
    // A law that searches for its verdict cannot be asked arbitrarily near
    // what it holds: just beyond its limit its sliding root is next to
    // zero, and its own search for it grows with every halving. The step
    // below is the one in which this was found, on a bowed string (the
    // string's mobility over a step cut to 1.0063e-5 s, the bow's force,
    // the fixture's friction): narrowed to rounding, the interface asked
    // this law 51 times, the last answers cost it 123 and 125 of the 128
    // iterations it may spend, and the 51st was refused. So such a law is
    // asked no nearer than the impulse tolerance needs: letting go a time
    // `t` early is a momentum of half the load's rate times `t` squared.
    let normal_n = 0.455_417_220_019_903_76;
    // Under two tolerances, 1e-10 and 2e-10 N s. Eight halvings bring the
    // bracket inside either: letting go where it then stands leaves
    // 8.0e-11 N s at most, and 3.2e-10 one halving before. So a rule that
    // stopped at half the tolerance, or at twice it, is a call off under
    // one of the two.
    for tolerance_n_s in [1e-10, 2e-10] {
        let load = InterfaceLoad {
            normal_n,
            nominal_pressure_pa: None,
            solve: TangentialSolveConfig {
                max_iterations: 128,
                maximum_force_n: 1e6,
                force_absolute_tolerance_n: 1e-12,
                force_relative_tolerance: 1e-12,
            },
            impulse_tolerance_n_s: tolerance_n_s,
            work_tolerance_j: 1e-12,
        };
        let interface = FrictionInterface::new(searching as CountedAt<HystereticFriction>, load)
            .expect("interface");
        let (dt_s, c) = (1.006_301_376e-5, 0.170_996_036_286_401_67);
        // The step before, 8.7808e-6 s long on a mobility of 0.149559, is held
        // under -0.32449 N, and this one has a mean of 0.64796 N: on the line
        // through the two the load goes from 0.12865 N to 1.16727 N in this
        // step and passes the 0.5465 N the law holds. The step before was on
        // another mobility, so the law's state carries no line prepared for
        // this one, and every verdict prepares it again (27 iterations here).
        let (last_s, last_c, last_n, load_n) = (8.7808e-6, 0.149_559, -0.324_485_3, 0.647_96);
        let rest = interface.start(300.0, 0.0).expect("start");
        let (cut, first, force_n) =
            asked_under(&interface, &rest, last_s, last_c, last_n).expect("held");
        assert!(cut.is_none() && first.at_rest());
        close(force_n, last_n, 1e-12);

        VERDICTS.with(|verdicts| verdicts.borrow_mut().clear());
        let answer = asked_under(&interface, &first, dt_s, c, load_n);
        let verdicts = VERDICTS.with(|verdicts| verdicts.borrow().clone());
        let costs: Vec<usize> = verdicts.iter().filter_map(|(_, cost)| *cost).collect();
        let refused = verdicts.len() - costs.len();
        let rate_n_s = (load_n - last_n) / (0.5 * (dt_s + last_s));
        let (from_n, to_n) = (load_n - 0.5 * dt_s * rate_n_s, load_n + 0.5 * dt_s * rate_n_s);
        let swing_n_s = dt_s * (to_n - from_n);
        // The halvings after which the bracket is inside the tolerance's time.
        let (mut halvings, mut width) = (0_usize, 1.0_f64);
        while 0.5 * width * width * swing_n_s > tolerance_n_s {
            width *= 0.5;
            halvings += 1;
        }
        println!(
            "a break-away under a law that searches: asked {} times ({halvings} halvings bring the \
             bracket inside the tolerance's time), {refused} refused, its iterations from {:?} to \
             {:?} of the 128 it may spend",
            verdicts.len(),
            costs.iter().min(),
            costs.iter().max()
        );
        let (cut, ..) = answer.expect("the instant is located, not refused");
        let asked_s = cut.expect("the step is cut at the instant");
        assert_eq!(verdicts.len(), 2 + halvings);
        // Every verdict was a search, and each left a quarter of the law's cap
        // unspent (61 of 128 at most, measured).
        assert!(costs.iter().all(|cost| (1..=96).contains(cost)), "{costs:?}");
        // The instant, by the law's static capacity mu_s N: the cut is short of
        // it by less than the bracket, and by less than the tolerance allows.
        let instant = (1.2 * normal_n - from_n) / (to_n - from_n);
        let short = instant - asked_s / dt_s;
        assert!(short >= 0.0 && short <= width, "{short} of the step, bracket {width}");
        assert!(0.5 * short * short * swing_n_s <= tolerance_n_s);
        assert_eq!(halvings, 8);
    }

    // A law whose verdict costs it nothing is narrowed to rounding, as it
    // was: Coulomb, 5 N held, the load rising from 4 N to a mean of 4.9 N
    // over two steps of 0.01 s, cut 0.55 / 0.9 of the way through.
    let interface =
        FrictionInterface::new(counted_dry as CountedAt<FrictionLaw>, pressed()).expect("interface");
    let (dt_s, c) = (0.01_f64, 0.5 * 0.01 / 2.0);
    let rest = interface.start(300.0, 0.0).expect("start");
    let (cut, first, _) = asked_under(&interface, &rest, dt_s, c, 4.0).expect("held");
    assert!(cut.is_none() && first.at_rest());
    VERDICTS.with(|verdicts| verdicts.borrow_mut().clear());
    let (cut, ..) = asked_under(&interface, &first, dt_s, c, 4.9).expect("answer");
    let verdicts = VERDICTS.with(|verdicts| verdicts.borrow().clone());
    println!(
        "the same under Coulomb: asked {} times, every verdict free: {}",
        verdicts.len(),
        verdicts.iter().all(|(_, cost)| *cost == Some(0))
    );
    close(cut.expect("cut"), (5.0 - 4.45) / 0.9 * dt_s, 1e-15);
    assert!(verdicts.iter().all(|(_, cost)| *cost == Some(0)), "{verdicts:?}");
    assert!(verdicts.len() > 40, "{}", verdicts.len());
}

/// `mu(v) = 0.2 + 0.3 / (1 + v / 0.1 m/s)`: 5 N at rest, falling to 2 N.
#[allow(clippy::unnecessary_wraps)]
fn falling_with_speed(_: f64) -> Result<FrictionLaw, TriboError> {
    Ok(FrictionLaw::Hyperbolic {
        static_mu: 0.5,
        kinetic_mu: 0.2,
        characteristic_speed: 0.1,
    })
}

#[test]
fn a_traction_that_falls_with_speed_is_solved_on_the_line_while_sliding() {
    // While sliding, the mean speed `v` solves `v = released - c N mu(v)`.
    // The interface brackets it between half the starting slip and the
    // released speed; the same root by bisection here.
    let interface = rubbing(falling_with_speed);
    let traction = |speed: f64| LOAD_N * (0.2 + 0.3 / (1.0 + speed / 0.1));
    let (dt_s, mass_kg) = (0.01_f64, 2.0);
    let c = 0.5 * dt_s / mass_kg;
    for (start, push_n) in [(0.5, 0.0), (0.05, 6.0), (-0.2, -1.0), (-0.01, -30.0)] {
        let way = f64::signum(start);
        let released = start + c * push_n;
        let (mut low, mut high) = (0.5 * f64::abs(start), way * released);
        for _ in 0..200 {
            let middle = 0.5 * (low + high);
            if middle - (way * released - c * traction(middle)) < 0.0 {
                low = middle;
            } else {
                high = middle;
            }
        }
        let sliding = interface.start(300.0, start).expect("start");
        let Answer::Took {
            state,
            receipt,
            force_n,
            heat_w,
        } = on_a_body(&interface, &sliding, dt_s, mass_kg, push_n)
        else {
            panic!("{start} m/s, {push_n} N cut the step");
        };
        let mean = released - c * force_n;
        close(way * mean, low, 1e-13);
        close(way * force_n, traction(low), 1e-11);
        close(state.slip_m_s(), 2.0 * mean - start, 0.0);
        close(heat_w, -force_n * mean, 4.0 * f64::EPSILON * traction(low));
        // The law's own dissipation, at the speed it was asked at, agrees
        // with the work of the ends within the declared tolerance.
        assert!(receipt.residual <= receipt.tolerance, "{receipt:?}");
        assert_eq!(receipt.tolerance, 1e-12);
    }
    // A root budget of two evaluations is not enough for this law, and the
    // refusal is the law owner's.
    let mut tight = pressed();
    tight.solve.max_iterations = 2;
    let short: Rubbing = FrictionInterface::new(falling_with_speed as LawAt, tight).expect("interface");
    let sliding = short.start(300.0, 0.5).expect("start");
    let refusal = on_line(&short, &sliding, dt_s, &[0.5, 0.0, 300.0], &[c, 0.0, 0.0, 0.0]);
    assert!(
        matches!(refusal, Err(TriboError::TangentialSolve { iterations: 2, .. })),
        "{refusal:?}"
    );
}

/// The dry law with two defects that can be switched on: it reports
/// `unaccounted_j` of dissipation that no work pays for, and it values
/// every state at `valued_j` of stored energy.
#[derive(Debug, Clone, PartialEq)]
struct Doctored {
    unaccounted_j: f64,
    valued_j: f64,
}

impl TangentialLaw for Doctored {
    type State = DryTangentialState;

    fn model_id(&self) -> &'static str {
        "test-doctored-dry-law"
    }

    fn initial_state(&self) -> DryTangentialState {
        dry(300.0).expect("law").initial_state()
    }

    fn stored_energy_j(&self, _: &DryTangentialState) -> f64 {
        self.valued_j
    }

    fn advance(
        &self,
        state: &DryTangentialState,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<DryTangentialState>, TriboError> {
        let mut step = dry(300.0)?.advance(state, input)?;
        step.dissipated_j += self.unaccounted_j;
        Ok(step)
    }
}

#[test]
fn the_friction_interface_refuses_what_it_cannot_take() {
    let interface = rubbing(dry);
    let rest = interface.start(300.0, 0.0).expect("start");
    let line = [2.5e-3, 0.0, 0.0, 0.0];
    let asked = |dt_s: f64, inputs: &[f64], compliance: &[f64]| {
        field(on_line(&interface, &rest, dt_s, inputs, compliance))
    };
    assert_eq!(asked(0.01, &[0.0, 0.0, 0.0], &line), "interface_temperature_k");
    assert_eq!(asked(0.01, &[0.0, 0.0, f64::NAN], &line), "interface_temperature_k");
    assert_eq!(asked(0.0, &[0.0, 0.0, 300.0], &line), "dt_s");
    assert_eq!(asked(0.01, &[f64::NAN, 0.0, 300.0], &line), "friction_interface_inputs");
    assert_eq!(asked(0.01, &[0.0, 0.0, 300.0], &[f64::INFINITY, 0.0, 0.0, 0.0]), "friction_interface_inputs");
    assert_eq!(asked(0.01, &[0.0, 300.0], &line), "friction_interface_inputs");
    assert_eq!(asked(0.01, &[0.0, 0.0, 300.0], &[2.5e-3]), "friction_interface_inputs");
    assert_eq!(field(interface.start(0.0, 0.0)), "interface_temperature_k");
    assert_eq!(field(interface.start(300.0, f64::NAN)), "slip_m_s");
    // The source's refusal passes through unchanged, wherever the law is
    // asked for.
    let refusing: Rubbing = FrictionInterface::new(
        (|kelvin| if kelvin > 305.0 { Err(TriboError::InvalidInput { field: "too hot" }) } else { dry(kelvin) })
            as LawAt,
        pressed(),
    )
    .expect("interface");
    let start = refusing.start(300.0, 0.2).expect("start");
    assert!(on_line(&refusing, &start, 0.01, &[0.2, 0.0, 305.0], &line).is_ok());
    assert_eq!(field(on_line(&refusing, &start, 0.01, &[0.2, 0.0, 306.0], &line)), "too hot");
    assert_eq!(field(refusing.start(306.0, 0.0)), "too hot");
    // A law that stores energy is taken if it holds no slip at rest and
    // slides at a prescribed slip. The spring of the storage tests does
    // neither: at rest it gives way.
    let sprung = FrictionInterface::new(|_| Ok(spring(4.0)), pressed()).expect("interface");
    let start = sprung.start(300.0, 0.0).expect("start");
    let mut next = start.clone();
    let mut outputs = [f64::NAN; 3];
    let refusal = with_cx(|cx| {
        sprung.propose_on_line(&start, 0.01, &[0.2, 0.0, 300.0], &line, &mut next, &mut outputs, cx)
    });
    assert_eq!(field(refusal), "friction_interface_static_map_kept_slip");
    // A law whose own accounts for the step do not close within the
    // declared 1e-12 J is refused, at rest and sliding; inside it the
    // defect is the receipt's residual.
    for (unaccounted_j, flow) in [(1.0e-9, 0.0), (1.0e-9, 0.2), (1.0e-13, 0.0), (1.0e-13, 0.2)] {
        let leaky = FrictionInterface::new(
            move |_| {
                Ok(Doctored {
                    unaccounted_j,
                    valued_j: 0.0,
                })
            },
            pressed(),
        )
        .expect("interface");
        let start = leaky.start(300.0, flow).expect("start");
        let mut next = start.clone();
        let answer = with_cx(|cx| {
            leaky.propose_on_line(&start, 0.01, &[flow, 0.0, 300.0], &line, &mut next, &mut outputs, cx)
        });
        if unaccounted_j > 1.0e-12 {
            assert_eq!(field(answer), "friction_interface_work_balance");
        } else {
            let Ok(Proposal::Step(receipt)) = answer else {
                panic!("{answer:?}");
            };
            close(receipt.residual, unaccounted_j, 1.0e-16);
            assert_eq!(receipt.tolerance, 1.0e-12);
        }
    }
    // A source whose law gives the held state another stored energy than
    // it was left with is refused, by any amount.
    let revaluing = FrictionInterface::new(
        |kelvin: f64| {
            Ok(Doctored {
                unaccounted_j: 0.0,
                valued_j: 1.0e-9 * (kelvin - 300.0),
            })
        },
        pressed(),
    )
    .expect("interface");
    let start = revaluing.start(300.0, 0.0).expect("start");
    let mut next = start.clone();
    let mut asked = |kelvin: f64| {
        with_cx(|cx| {
            revaluing.propose_on_line(&start, 0.01, &[0.0, 0.0, kelvin], &line, &mut next, &mut outputs, cx)
        })
    };
    assert!(matches!(asked(300.0), Ok(Proposal::Step(_))));
    assert_eq!(field(asked(300.5)), "friction_interface_law_changed_stored_energy");
    // What an interface cannot be pressed with.
    let built = |change: fn(&mut InterfaceLoad)| {
        let mut load = pressed();
        change(&mut load);
        FrictionInterface::new(dry as LawAt, load).map(|_| ())
    };
    assert_eq!(field(built(|load| load.normal_n = -1.0)), "normal_force");
    assert_eq!(field(built(|load| load.nominal_pressure_pa = Some(f64::NAN))), "contact_pressure_pa");
    assert_eq!(field(built(|load| load.work_tolerance_j = -1.0)), "work_tolerance_j");
    assert_eq!(field(built(|load| load.impulse_tolerance_n_s = f64::INFINITY)), "impulse_tolerance_n_s");
    assert_eq!(field(built(|load| load.solve.force_relative_tolerance = 1.0)), "force_relative_tolerance");
    assert!(built(|load| load.normal_n = 0.0).is_ok());
}

/// What the end of a step moves per newton beside the told neighbour of the
/// tests below, as a multiple of what its mean moves: for a body under a
/// held force on the midpoint rule it is 2, for a modal string neither.
const TOLD: f64 = 1.917;

/// One proposal of an interface that reads `a`'s end, with `b` given at no
/// flow: the mean line `(released, c)` and the end line
/// `(released_end, c_end)` of the slip.
fn on_the_end_line(
    interface: &Rubbing,
    state: &Rubbed,
    dt_s: f64,
    (released, c): (f64, f64),
    (released_end, c_end): (f64, f64),
) -> Result<Answer, TriboError> {
    on_line(
        interface,
        state,
        dt_s,
        &[released, 0.0, 300.0, released_end],
        &[c, 0.0, 0.0, 0.0, c_end, 0.0],
    )
}

/// One proposal beside a neighbour at `a` that a force of `push_n` drives
/// along the slip and whose step is exact under a held force. From the slip
/// the last step ended with, a newton moves the step's mean slip by
/// `c = dt / (2 m)` and its end by `ratio c`.
fn beside(
    interface: &Rubbing,
    state: &Rubbed,
    dt_s: f64,
    mass_kg: f64,
    push_n: f64,
    ratio: f64,
) -> Answer {
    let c = 0.5 * dt_s / mass_kg;
    let start = state.slip_m_s();
    on_the_end_line(
        interface,
        state,
        dt_s,
        (start + c * push_n, c),
        (start + ratio * c * push_n, ratio * c),
    )
    .expect("answer")
}

#[test]
fn an_interface_that_reads_an_end_holds_a_pair_at_rest_at_no_end_slip_and_lets_go_on_its_line() {
    // `a` reads the end its neighbour tells; `b` is given. The end input
    // comes after `heat`, and the line names it.
    let interface = rubbing(dry).reading_ends([true, false]).expect("ends");
    let names: Vec<&str> = interface
        .ports()
        .iter()
        .map(|port| port.id().as_str())
        .collect();
    assert_eq!(names, ["a", "b", "heat", "a.end"]);
    assert_eq!(interface.line_ports(), [0, 1]);
    assert_eq!(interface.line_ends(), [Some(3), None]);
    assert_eq!(rubbing(dry).line_ends(), [None, None]);
    let rest = interface.start(300.0, 0.0).expect("start");
    let (dt_s, mass_kg) = (0.01, 2.0);
    let c = 0.5 * dt_s / mass_kg;
    // Up to 5 N the pair is held with the force the push needs, and the
    // step ends with no slip at all.
    for push_n in [0.0, 1.0, -4.0, 4.999, -4.999] {
        let Answer::Took {
            state,
            receipt,
            force_n,
            heat_w,
        } = beside(&interface, &rest, dt_s, mass_kg, push_n, TOLD)
        else {
            panic!("{push_n} N cut the step");
        };
        close(force_n, push_n, 8.0 * f64::EPSILON * 5.0);
        assert!(state.at_rest(), "{push_n} N");
        assert_eq!(state.slip_m_s(), 0.0);
        assert_eq!(state.tangential().slip_m(), [0.0; 2]);
        assert_eq!(
            (heat_w, receipt.dissipated_j, receipt.stored_change_j),
            (0.0, 0.0, 0.0)
        );
        close(receipt.numerical_loss_j, 0.0, 1e-18);
    }
    // A neighbour whose mean does not pass through the end line's rest, as
    // on a string, where a step's mean is not where its end is. The pair is
    // held with the end line's force all the same, the mean slip `r - c F`
    // is what the neighbour's mean does, and the work of the two ends,
    // either way, is the step's numerical loss: no heat, no entropy.
    for off_m_s in [1.0e-3, -1.0e-3] {
        let released = 4.0 * c + off_m_s;
        let end_line = (4.0 * TOLD * c, TOLD * c);
        let Ok(Answer::Took {
            state,
            receipt,
            force_n,
            heat_w,
        }) = on_the_end_line(&interface, &rest, dt_s, (released, c), end_line)
        else {
            panic!("{off_m_s} m/s: not taken");
        };
        close(force_n, 4.0, 8.0 * f64::EPSILON * 4.0);
        assert!(state.at_rest());
        assert_eq!(state.slip_m_s(), 0.0);
        let work_j = dt_s * force_n * (released - c * force_n);
        assert_eq!(receipt.numerical_loss_j, work_j);
        assert!(work_j * off_m_s > 0.0, "{work_j}");
        assert_eq!(
            (
                heat_w,
                receipt.dissipated_j,
                receipt.entropy_generated_j_per_k
            ),
            (0.0, 0.0, 0.0)
        );
    }
    // Beyond 5 N it lets go against 3 N, the traction the law solves on the
    // end line. The step ends where the neighbour says, `r_e - c_e F` (on
    // the midpoint rule, twice the mean less the start), and its heat is
    // the work of the two ends at the mean slip.
    for push_n in [6.0_f64, -8.0, 50.0] {
        let Answer::Took {
            state,
            receipt,
            force_n,
            heat_w,
        } = beside(&interface, &rest, dt_s, mass_kg, push_n, TOLD)
        else {
            panic!("{push_n} N cut the step");
        };
        let way = push_n.signum();
        assert!(!state.at_rest());
        assert_eq!(force_n, SLIDES * LOAD_N * way);
        let (mean, ends) = (c * (push_n - 3.0 * way), TOLD * c * (push_n - 3.0 * way));
        close(state.slip_m_s(), ends, 8.0 * f64::EPSILON * ends.abs());
        close(
            state.tangential().slip_m()[0],
            dt_s * ends,
            16.0 * f64::EPSILON * dt_s * ends.abs(),
        );
        close(
            heat_w,
            -3.0 * mean.abs(),
            16.0 * f64::EPSILON * 3.0 * mean.abs(),
        );
        let heat_j = dt_s * 3.0 * mean.abs();
        close(receipt.dissipated_j, heat_j, 16.0 * f64::EPSILON * heat_j);
        assert_eq!(receipt.numerical_loss_j, 0.0);
    }
}

#[test]
fn beside_a_told_end_a_sliding_pair_ends_its_steps_on_the_end_line_and_is_cut_where_it_stops() {
    let interface = rubbing(dry).reading_ends([true, false]).expect("ends");
    let mass_kg = 2.0;
    for start in [1.0_f64, -0.4] {
        let way = start.signum();
        let sliding = interface.start(300.0, start).expect("start");
        // Sliding on against 3 N: the law is asked at the step's mean slip,
        // `start + c (push - 3)`, and the step ends at
        // `start + c_e (push - 3)`.
        for (dt_s, push_n) in [(0.01, 5.0 * way), (0.01, 0.0), (0.05, -2.0 * way)] {
            let Answer::Took {
                state,
                receipt,
                force_n,
                heat_w,
            } = beside(&interface, &sliding, dt_s, mass_kg, push_n, TOLD)
            else {
                panic!("{start} m/s, {push_n} N cut the step");
            };
            let c = 0.5 * dt_s / mass_kg;
            assert_eq!(force_n, SLIDES * LOAD_N * way);
            let mean = start + c * (push_n - 3.0 * way);
            close(
                state.slip_m_s(),
                start + TOLD * c * (push_n - 3.0 * way),
                8.0 * f64::EPSILON,
            );
            close(
                state.tangential().slip_m()[0],
                dt_s * mean,
                8.0 * f64::EPSILON * dt_s,
            );
            close(heat_w, -3.0 * mean.abs(), 16.0 * f64::EPSILON * 3.0);
            close(
                receipt.dissipated_j,
                dt_s * 3.0 * mean.abs(),
                16.0 * f64::EPSILON * 3.0 * dt_s,
            );
            assert_eq!(receipt.numerical_loss_j, 0.0);
        }
        // The end slip reaches zero after `2 m abs(start) / (TOLD (3 - way
        // push))`, by the end line, which for this neighbour grows with the
        // step: the same instant from every step length that ends beyond
        // it. That step ends at rest under the law's own 3 N, asked at the
        // step's mean, `start (1 - 1 / TOLD)` (half the start on the
        // midpoint rule). A step that ends short of it slides on.
        for push_n in [0.0, 2.0 * way, -7.0 * way] {
            let stop_s = 2.0 * mass_kg * start.abs() / (TOLD * (3.0 - way * push_n));
            let mut found_s = f64::NAN;
            for beyond in [1.0001, 1.5, 4.0, 1000.0] {
                let Answer::Cut(asked_s) =
                    beside(&interface, &sliding, beyond * stop_s, mass_kg, push_n, TOLD)
                else {
                    panic!("{start} m/s, {push_n} N, {beyond}: not cut");
                };
                close(asked_s, stop_s, 64.0 * f64::EPSILON * stop_s);
                found_s = asked_s;
            }
            let Answer::Took {
                state,
                receipt,
                force_n,
                heat_w,
            } = beside(&interface, &sliding, found_s, mass_kg, push_n, TOLD)
            else {
                panic!("{start} m/s, {push_n} N: the located step is cut again");
            };
            assert!(state.at_rest());
            assert!(
                state.slip_m_s().abs() <= 8.0 * f64::EPSILON,
                "{}",
                state.slip_m_s()
            );
            assert_eq!(force_n, SLIDES * LOAD_N * way);
            let mean = start * (1.0 - 1.0 / TOLD);
            close(heat_w, -3.0 * mean.abs(), 1e-9);
            close(receipt.dissipated_j, found_s * 3.0 * mean.abs(), 1e-9);
            close(state.tangential().slip_m()[0], found_s * mean, 1e-15);
            // The next step at rest holds the pair at no slip at all: what
            // the stop left is in the neighbour's released end.
            let Answer::Took {
                state: held,
                force_n,
                ..
            } = beside(&interface, &state, 0.01, mass_kg, way, TOLD)
            else {
                panic!("{start} m/s, {push_n} N: the step after the stop was cut");
            };
            assert!(held.at_rest());
            assert_eq!(held.slip_m_s(), 0.0);
            close(force_n, way, 1e-12);
            let Answer::Took { state, .. } =
                beside(&interface, &sliding, 0.999 * stop_s, mass_kg, push_n, TOLD)
            else {
                panic!("{start} m/s, {push_n} N: a step short of the stop is cut");
            };
            assert!(!state.at_rest());
            close(state.slip_m_s(), 1.0e-3 * start, 1.0e-12);
        }
    }
}

#[test]
fn beside_a_told_end_a_break_away_is_located_from_the_load_on_the_end_line() {
    // As on the midpoint rule: from 4 N to a mean of 4.9 N over two steps
    // of 0.01 s the load rises at 90 N/s, and the 5 N held are reached
    // after 0.55 / 0.9 of the step. The load is read off the end line,
    // `(r_e - s0) / c_e`. The neighbour's mean line passes 0.3 N beyond it
    // (on the midpoint rows a load of 4.3 N and then 5.2 N, which lets go
    // from the step's start) and moves nothing.
    let interface = rubbing(dry).reading_ends([true, false]).expect("ends");
    let mass_kg = 2.0;
    let pushed = |state: &Rubbed, dt_s: f64, push_n: f64| {
        let c = 0.5 * dt_s / mass_kg;
        let end_line = (TOLD * c * push_n, TOLD * c);
        on_the_end_line(&interface, state, dt_s, (c * (push_n + 0.3), c), end_line).expect("answer")
    };
    let dt_s = 0.01;
    let rest = interface.start(300.0, 0.0).expect("start");
    let Answer::Took {
        state: first,
        force_n,
        ..
    } = pushed(&rest, dt_s, 4.0)
    else {
        panic!("4 N cut the step");
    };
    assert!(first.at_rest());
    close(force_n, 4.0, 1e-12);
    let fraction = (5.0 - 4.45) / 0.9;
    let Answer::Cut(asked_s) = pushed(&first, dt_s, 4.9) else {
        panic!("a load that passes the limit inside the step was not cut");
    };
    close(asked_s, fraction * dt_s, 1e-15);
    // The step of that length is held and ends at the limit; the next one
    // lets go from its start, against 3 N.
    let Answer::Took {
        state: at_the_limit,
        ..
    } = pushed(&first, asked_s, 4.45 + 0.5 * fraction * 0.9)
    else {
        panic!("the located step was cut again");
    };
    assert!(at_the_limit.at_rest());
    let Answer::Took { state, force_n, .. } = pushed(&at_the_limit, dt_s, 5.0 + 90.0 * 0.5 * dt_s)
    else {
        panic!("the step at the limit was cut");
    };
    assert!(!state.at_rest());
    assert_eq!(force_n, SLIDES * LOAD_N);
}

#[test]
fn a_told_end_that_ends_its_steps_by_the_midpoint_rule_is_the_midpoint_interface() {
    // A body under a held force, told: its end moves twice what its mean
    // does. Every verdict, force, slip record and located instant is then
    // the midpoint table's, to rounding: from rest, sliding, and after a
    // step whose load the next one is located from.
    let told = rubbing(dry).reading_ends([true, false]).expect("ends");
    let midpoint = rubbing(dry);
    let mass_kg = 2.0;
    let (mut taken, mut cut) = (0, 0);
    let mut compare = |from: (&Rubbed, &Rubbed), dt_s: f64, push_n: f64| match (
        beside(&told, from.0, dt_s, mass_kg, push_n, 2.0),
        on_a_body(&midpoint, from.1, dt_s, mass_kg, push_n),
    ) {
        (Answer::Cut(asked_s), Answer::Cut(expected_s)) => {
            close(asked_s, expected_s, 64.0 * f64::EPSILON * expected_s);
            cut += 1;
            None
        }
        (
            Answer::Took { state, force_n, .. },
            Answer::Took {
                state: expected,
                force_n: expected_n,
                ..
            },
        ) => {
            close(force_n, expected_n, 1e-12);
            close(state.slip_m_s(), expected.slip_m_s(), 1e-13);
            assert_eq!(state.at_rest(), expected.at_rest(), "{dt_s} s, {push_n} N");
            taken += 1;
            Some((state, expected))
        }
        (told, expected) => panic!("{dt_s} s, {push_n} N: told {told:?}, midpoint {expected:?}"),
    };
    for start in [0.0, 1.0, -0.4, 1.0e-3] {
        let from = (
            told.start(300.0, start).expect("start"),
            midpoint.start(300.0, start).expect("start"),
        );
        for dt_s in [0.01, 0.05, 0.5] {
            for push_n in [0.0, 1.0, -4.9, 6.0, -8.0, 2.0] {
                compare((&from.0, &from.1), dt_s, push_n);
            }
        }
    }
    let rest = (
        told.start(300.0, 0.0).expect("start"),
        midpoint.start(300.0, 0.0).expect("start"),
    );
    let (held, expected) = compare((&rest.0, &rest.1), 0.01, 4.0).expect("held");
    for push_n in [4.4, 4.9, 7.0, -4.9] {
        compare((&held, &expected), 0.01, push_n);
    }
    assert!(taken > 40 && cut > 10, "{taken} taken, {cut} cut");
    println!("a told end on the midpoint rule: {taken} steps taken and {cut} cut as on that rule");
}

#[test]
fn with_both_ends_told_the_end_line_is_the_difference_of_their_rows() {
    // `a` and `b` each read their end, and neither is given. The line hands
    // a row per end; taken along the efforts (F, -F) they are
    // `E_a = 4.8e-3 + 1.5e-3` and `E_b = -2.0e-3 - 1.9e-2`, so the end slip
    // moves by `c_e = E_a - E_b = 2.73e-2` per newton, and its released
    // value is that of `a`'s end less `b`'s.
    let interface = rubbing(dry).reading_ends([true, true]).expect("ends");
    let names: Vec<&str> = interface
        .ports()
        .iter()
        .map(|port| port.id().as_str())
        .collect();
    assert_eq!(names, ["a", "b", "heat", "a.end", "b.end"]);
    assert_eq!(interface.line_ends(), [Some(3), Some(4)]);
    let line = [
        2.5e-3, -1.0e-3, -1.0e-3, 1.0e-2, 4.8e-3, -1.5e-3, -2.0e-3, 1.9e-2,
    ];
    let (c, c_end, dt_s) = (1.45e-2, 2.73e-2, 0.01);
    // Both ends move at 0.1 m/s, and `a` by the slip besides.
    let asked = |state: &Rubbed, released: f64, released_end: f64| {
        let inputs = [0.1 + released, 0.1, 300.0, 0.1 + released_end, 0.1];
        on_line(&interface, state, dt_s, &inputs, &line).expect("answer")
    };
    let rest = interface.start(300.0, 0.0).expect("start");
    let Answer::Took { state, force_n, .. } = asked(&rest, 4.0 * c, 4.0 * c_end) else {
        panic!("4 N cut the step");
    };
    assert!(state.at_rest());
    assert_eq!(state.slip_m_s(), 0.0);
    close(force_n, 4.0, 1e-13);
    let Answer::Took {
        state,
        force_n,
        heat_w,
        ..
    } = asked(&rest, 6.0 * c, 6.0 * c_end)
    else {
        panic!("6 N cut the step");
    };
    assert!(!state.at_rest());
    assert_eq!(force_n, SLIDES * LOAD_N);
    close(state.slip_m_s(), 3.0 * c_end, 1e-15);
    close(heat_w, -3.0 * 3.0 * c, 1e-14);
    // Sliding at 0.5 m/s under a push of 1 N it slides on, and ends at
    // `0.5 + c_e (1 - 3)`.
    let sliding = interface.start(300.0, 0.5).expect("start");
    let Answer::Took {
        state,
        force_n,
        heat_w,
        ..
    } = asked(&sliding, 0.5 + c, 0.5 + c_end)
    else {
        panic!("sliding on was cut");
    };
    assert_eq!(force_n, SLIDES * LOAD_N);
    close(state.slip_m_s(), 0.5 - 2.0 * c_end, 1e-15);
    close(heat_w, -3.0 * (0.5 - 2.0 * c), 1e-14);
}

#[test]
fn the_told_interface_refuses_what_it_cannot_take() {
    assert_eq!(
        field(rubbing(dry).reading_ends([false; 2])),
        "friction_interface_end_inputs"
    );
    let twice = rubbing(dry)
        .reading_ends([true, false])
        .and_then(|told| told.reading_ends([false, true]));
    assert_eq!(field(twice), "friction_interface_end_inputs");
    let interface = rubbing(dry).reading_ends([true, false]).expect("ends");
    let rest = interface.start(300.0, 0.0).expect("start");
    let (dt_s, c) = (0.01, 2.5e-3);
    let line = [c, 0.0, 0.0, 0.0, TOLD * c, 0.0];
    let asked = |state: &Rubbed, inputs: &[f64], compliance: &[f64]| {
        field(on_line(&interface, state, dt_s, inputs, compliance))
    };
    // The end input and its row are counted, and finite.
    let miscounted: [(&[f64], &[f64]); 6] = [
        (&[0.01, 0.0, 300.0], &line),
        (&[0.01, 0.0, 300.0, 0.01, 0.0], &line),
        (&[0.01, 0.0, 300.0, 0.01], &line[..4]),
        (
            &[0.01, 0.0, 300.0, 0.01],
            &[c, 0.0, 0.0, 0.0, TOLD * c, 0.0, 0.0, 0.0],
        ),
        (&[0.01, 0.0, 300.0, f64::NAN], &line),
        (
            &[0.01, 0.0, 300.0, 0.01],
            &[c, 0.0, 0.0, 0.0, f64::INFINITY, 0.0],
        ),
    ];
    for (inputs, compliance) in miscounted {
        let refused = asked(&rest, inputs, compliance);
        assert_eq!(
            refused, "friction_interface_inputs",
            "{inputs:?} {compliance:?}"
        );
    }
    // `b` neither reads an end nor is given: its row of the line is not zero.
    let untold = [c, 0.0, 0.0, 1.0e-3, TOLD * c, 0.0];
    assert_eq!(
        asked(&rest, &[0.01, 0.0, 300.0, 0.01], &untold),
        "friction_interface_untold_end"
    );
    // An end line that takes nothing from the slip, or adds to it, where
    // the mean line takes from it.
    for c_end in [0.0, -TOLD * c] {
        let end_line = [c, 0.0, 0.0, 0.0, c_end, 0.0];
        let refused = asked(&rest, &[0.01, 0.0, 300.0, 0.01], &end_line);
        assert_eq!(refused, "friction_interface_end_compliance");
    }
    // A sliding pair that the end line has 400 N along and the mean line
    // barely moving: the stop step's mean would run against the slip. The
    // law holds 5 N, so the end cannot stop, and the step slides on under
    // its 3 N.
    let sliding = interface.start(300.0, 1.0).expect("start");
    let Ok(Answer::Took {
        state,
        force_n,
        heat_w,
        ..
    }) = on_line(
        &interface,
        &sliding,
        dt_s,
        &[0.01, 0.0, 300.0, 1.0],
        &[c, 0.0, 0.0, 0.0, c, 0.0],
    )
    else {
        panic!("400 N along was refused or cut");
    };
    assert_eq!(force_n, SLIDES * LOAD_N);
    close(state.slip_m_s(), 1.0 - 3.0 * c, 1e-15);
    close(heat_w, -3.0 * (0.01 - 3.0 * c), 1e-15);
    // At 4 N along, which the law holds, the stop step's mean would still
    // run against the slip: nothing says what that step is.
    let refused = asked(
        &sliding,
        &[0.004, 0.0, 300.0, 4.0 * c],
        &[c, 0.0, 0.0, 0.0, c, 0.0],
    );
    assert_eq!(refused, "friction_interface_stop_against_the_slip");
    // A pair that lets go from rest on the end line while its mean slip
    // runs against the force: the ends would take heat out of it.
    let refused = asked(&rest, &[3.0 * c - 0.5, 0.0, 300.0, 6.0 * TOLD * c], &line);
    assert_eq!(refused, "friction_interface_let_go_heat");
}

#[test]
fn the_told_interface_conforms_on_its_line_and_books_the_work_of_its_ends() {
    // On its line it is pure, its step is audited at the flows on the line
    // (the work of the two ends is what the receipt books, as heat or as a
    // defect), and with a zero line it is `propose`. Neighbours told at the
    // ratio of these tests and at the midpoint rule's, and at rest mean
    // lines that pass off the end line's rest.
    let interface = rubbing(dry).reading_ends([true, false]).expect("ends");
    let rest = interface.start(300.0, 0.0).expect("start");
    let (dt_s, c) = (0.01, 2.5e-3);
    let stale = {
        let inputs = [0.4, 0.1, 310.0, 0.4];
        let Ok(Answer::Took { state, .. }) = on_line(&interface, &rest, 0.01, &inputs, &[0.0; 6])
        else {
            panic!("stale");
        };
        state
    };
    let (mut stepped, mut cut) = (0, 0);
    for start in [0.0, 0.3, -0.3, 1.0e-3] {
        let state = interface.start(300.0, start).expect("start");
        let offs: &[f64] = if start == 0.0 {
            &[-1.0e-3, 0.0, 1.0e-3]
        } else {
            &[0.0]
        };
        let points = offs.iter().flat_map(|&off_m_s| {
            [TOLD, 2.0].into_iter().flat_map(move |ratio| {
                [0.0, 1.0, 4.9, -6.0, 12.0].map(|push_n| (ratio, push_n, off_m_s))
            })
        });
        for (ratio, push_n, off_m_s) in points {
            let inputs = [
                start + c * push_n + off_m_s,
                0.0,
                300.0,
                start + ratio * c * push_n,
            ];
            let case = Case {
                state: &state,
                stale: &stale,
                dt_s,
                inputs: &inputs,
            };
            let compliance = [c, 0.0, 0.0, 0.0, ratio * c, 0.0];
            let audited = with_cx(|cx| check_load_line(&interface, &case, &compliance, cx))
                .unwrap_or_else(|failure| {
                    panic!("{start} {ratio} {push_n} {off_m_s}: {failure:?}")
                });
            match audited.observed {
                Observed::Shortened { .. } => cut += 1,
                Observed::Stepped {
                    receipt, outputs, ..
                } => {
                    stepped += 1;
                    let work_j = dt_s * outputs[0] * (audited.flows[0] - audited.flows[1]);
                    close(
                        receipt.dissipated_j + receipt.numerical_loss_j,
                        work_j,
                        1e-16,
                    );
                    assert!(
                        outputs[0].abs() <= HOLDS * LOAD_N * (1.0 + 1e-12),
                        "{outputs:?}"
                    );
                }
            }
        }
    }
    assert!(stepped > 40 && cut >= 6, "{stepped} stepped, {cut} cut");
    println!("the told interface on its line: {stepped} points stepped and audited, {cut} cut");
}
