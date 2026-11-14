//! G0: the Component contract v0 on a reference implementation.
//!
//! The reference is a linear spring-mass-damper with a mechanical port and a
//! heat port, stepped by the implicit midpoint rule. For a linear system that
//! rule balances energy exactly: the work of the held force on the mean
//! velocity equals the stored change plus the damper loss, and the loss
//! leaves through the heat port. Five deliberately broken variants each fail
//! the law they break. A sixth keeps a cache in a `Cell`, keyed by the bits of
//! what it is proposed from, that `propose` never shows; it passes, as a cache
//! may. Without its damper the reference is reversible in time: it passes
//! `check_time_reversal` within bounds derived from its own arithmetic, and
//! is refused with its damper, with a reversal that keeps the velocity, and
//! with one that is not an involution. Through its drive port, the audit
//! admits a balance or an endpoint one subnormal off at a zero tolerance (the
//! f64 underflow allowance) and still refuses a defect of 1e-300 J.

use core::cell::Cell;
use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_ports::component::conformance::{
    Case, Failure, Observed, ReversalAllowance, check_component, check_port_tangent,
    check_time_reversal,
};
use fs_ports::component::{
    Carrier, Class, Component, ContractViolation, PortLayout, PortRole, PortSpec, PortTangent,
    Proposal, Receipt, StepAudit, audit_step,
};
use fs_ports::{PortKind, StableId};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 1,
            kernel_id: 2,
            tile: 3,
            iteration: 4,
        };
        let cx = Cx::new(gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Oscillator {
    position_m: f64,
    velocity_m_s: f64,
}

/// Which law, if any, the component breaks on purpose.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Fault {
    None,
    /// Reports only 90 % of the damper loss as heat.
    HeatLeak,
    /// Adds the stale candidate's position to the result.
    ReadsCandidate,
    /// Returns the VJP with the wrong sign.
    VjpSign,
    /// Counts its proposals in a `Cell` through `&self`, and moves the mean
    /// velocity it returns by that many units in its last place. It stands
    /// for anything a component carries between calls that reaches what a
    /// call returns, such as a warm start.
    CountsIntoOutput,
    /// Never polls `cx`, so it takes a step asked under a cancelled context.
    IgnoresCancellation,
    /// Not a fault: keeps, in a `Cell` through `&self`, the mean velocity of
    /// its last proposal under the bits of the accepted state, the step and
    /// the force, and returns it when all of them repeat. What it returns is
    /// what it would have computed, so `propose` never shows the cache.
    CachesItsLastMean,
}

/// The bits a cached mean velocity is keyed by: the accepted state's
/// position and velocity, the step and the force.
type MeanKey = [u64; 4];

#[derive(Debug)]
struct SpringMassDamper {
    mass_kg: f64,
    stiffness_n_m: f64,
    damping_n_s_m: f64,
    ports: [PortSpec; 2],
    fault: Fault,
    /// Proposals made, kept only by `CountsIntoOutput`.
    proposals: Cell<u32>,
    /// The last mean velocity and its key, kept only by `CachesItsLastMean`.
    cache: Cell<Option<(MeanKey, f64)>>,
    /// Proposals that took their mean velocity from `cache`.
    cache_hits: Cell<u32>,
}

impl SpringMassDamper {
    fn new(fault: Fault) -> Self {
        let drive = PortSpec::power(
            StableId::new("drive").expect("id"),
            PortKind::MechanicalForceVelocity,
            PortRole::EffortIn,
            NonZeroUsize::MIN,
        )
        .expect("power port");
        let heat =
            PortSpec::heat(StableId::new("heat").expect("id"), PortRole::EffortIn).expect("heat");
        Self {
            mass_kg: 0.25,
            stiffness_n_m: 400.0,
            damping_n_s_m: 1.5,
            ports: [drive, heat],
            fault,
            proposals: Cell::new(0),
            cache: Cell::new(None),
            cache_hits: Cell::new(0),
        }
    }

    /// `mean_velocity`'s mean, taken from the cache when the accepted state,
    /// the step and the force are those it was kept for, bit for bit.
    fn cached_mean(&self, old: &Oscillator, dt_s: f64, force_n: f64) -> f64 {
        let key = [
            old.position_m.to_bits(),
            old.velocity_m_s.to_bits(),
            dt_s.to_bits(),
            force_n.to_bits(),
        ];
        if let Some((kept, mean)) = self.cache.get()
            && kept == key
        {
            self.cache_hits.set(self.cache_hits.get() + 1);
            return mean;
        }
        let (mean, _) = self.mean_velocity(old, dt_s, force_n);
        self.cache.set(Some((key, mean)));
        mean
    }

    /// Mean velocity over the step and its derivative with respect to force.
    fn mean_velocity(&self, old: &Oscillator, dt_s: f64, force_n: f64) -> (f64, f64) {
        let half = 0.5 * dt_s / self.mass_kg;
        let denominator = 1.0 + half * self.damping_n_s_m + 0.5 * dt_s * half * self.stiffness_n_m;
        let mean = (old.velocity_m_s + half * (force_n - self.stiffness_n_m * old.position_m))
            / denominator;
        (mean, half / denominator)
    }
}

impl Component for SpringMassDamper {
    type State = Oscillator;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, state: &Oscillator) -> Result<f64, Cancelled> {
        Ok(0.5 * self.mass_kg * state.velocity_m_s * state.velocity_m_s
            + 0.5 * self.stiffness_n_m * state.position_m * state.position_m)
    }

    fn propose(
        &self,
        old: &Oscillator,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Oscillator,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        if self.fault != Fault::IgnoresCancellation {
            cx.checkpoint()?;
        }
        let (force_n, temperature_k) = (inputs[0], inputs[1]);
        let stale_position = new.position_m;
        let mean = if self.fault == Fault::CachesItsLastMean {
            self.cached_mean(old, dt_s, force_n)
        } else {
            self.mean_velocity(old, dt_s, force_n).0
        };
        let mid_position = old.position_m + 0.5 * dt_s * mean;
        new.position_m = old.position_m + dt_s * mean;
        new.velocity_m_s = 2.0 * mean - old.velocity_m_s;
        if self.fault == Fault::ReadsCandidate {
            new.position_m += 1.0e-3 * stale_position;
        }
        let loss_w = self.damping_n_s_m * mean * mean;
        let heat_out_w = if self.fault == Fault::HeatLeak {
            0.9 * loss_w
        } else {
            loss_w
        };
        outputs[0] = mean;
        if self.fault == Fault::CountsIntoOutput {
            let count = self.proposals.get() + 1;
            self.proposals.set(count);
            // `count` units in the last place of a finite, nonzero mean.
            outputs[0] = f64::from_bits(mean.to_bits() + u64::from(count));
        }
        outputs[1] = -heat_out_w;
        // Stored change without cancellation: m v_mean dv + k x_mid dx.
        let stored_change_j = self.mass_kg * mean * (new.velocity_m_s - old.velocity_m_s)
            + self.stiffness_n_m * mid_position * (dt_s * mean);
        let scale_j = (dt_s * force_n * mean).abs()
            + dt_s * loss_w
            + self.stored_energy_j(old)?
            + self.stored_energy_j(new)?;
        Ok(Proposal::Step(Receipt {
            stored_change_j,
            dissipated_j: dt_s * loss_w,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dt_s * loss_w / temperature_k,
            energy_tolerance_j: 16.0 * f64::EPSILON * scale_j,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

impl PortTangent for SpringMassDamper {
    fn jvp(
        &self,
        old: &Oscillator,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        cx.checkpoint()?;
        let (mean, slope) = self.mean_velocity(old, dt_s, inputs[0]);
        d_outputs[0] = slope * d_inputs[0];
        d_outputs[1] = -2.0 * self.damping_n_s_m * mean * slope * d_inputs[0];
        Ok(())
    }

    fn vjp(
        &self,
        old: &Oscillator,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        cx.checkpoint()?;
        let (mean, slope) = self.mean_velocity(old, dt_s, inputs[0]);
        let sign = if self.fault == Fault::VjpSign {
            -1.0
        } else {
            1.0
        };
        bar_inputs[0] = sign
            * (slope * bar_outputs[0] - 2.0 * self.damping_n_s_m * mean * slope * bar_outputs[1]);
        bar_inputs[1] = 0.0;
        Ok(())
    }
}

const START: Oscillator = Oscillator {
    position_m: 2.0e-3,
    velocity_m_s: -0.3,
};
const STALE: Oscillator = Oscillator {
    position_m: 7.0,
    velocity_m_s: 11.0,
};
const INPUTS: [f64; 2] = [1.25, 293.15];

fn case() -> Case<'static, Oscillator> {
    Case {
        state: &START,
        stale: &STALE,
        dt_s: 1.0e-3,
        inputs: &INPUTS,
    }
}

#[test]
fn reference_component_passes_conformance() {
    let gate = CancelGate::new_clock_free();
    let component = SpringMassDamper::new(Fault::None);
    let observed = with_cx(&gate, |cx| check_component(&component, &case(), cx)).expect("conforms");
    let Observed::Stepped {
        receipt,
        audit,
        outputs,
    } = observed
    else {
        panic!("the reference takes the step");
    };
    assert!(receipt.dissipated_j > 0.0);
    // The damper loss leaves as heat: heat-port energy in is minus the loss.
    let heat_in_j = component.ports()[1].energy_in_j(1.0e-3, &INPUTS[1..], &outputs[1..]);
    assert!((heat_in_j + receipt.dissipated_j).abs() <= receipt.energy_tolerance_j);
    // With a bound that does not come from the receipt: the net port energy
    // is the stored change to 1e-13 of the energy that crossed the ports.
    let work_in_j = component.ports()[0].energy_in_j(1.0e-3, &INPUTS[..1], &outputs[..1]);
    let crossed_j = work_in_j.abs() + heat_in_j.abs();
    assert!((audit.port_energy_in_j - receipt.stored_change_j).abs() <= 1.0e-13 * crossed_j);
    // Entropy leaves through the heat port at the port temperature.
    let entropy_out = -audit.port_entropy_in_j_per_k;
    assert!((entropy_out - receipt.entropy_generated_j_per_k).abs() <= 1.0e-18);
    with_cx(&gate, |cx| {
        check_port_tangent(&component, &case(), 1.0e-7, 1.0e-12, cx)
    })
    .expect("tangent conforms");
}

#[test]
fn energy_closes_over_many_steps_without_drift() {
    let gate = CancelGate::new_clock_free();
    let component = SpringMassDamper::new(Fault::None);
    let layout = PortLayout::new(component.ports()).expect("layout");
    let dt_s = 2.0e-4;
    let (mut state, mut next) = (START, START);
    let mut outputs = vec![0.0; layout.output_len()];
    let start_j = component.stored_energy_j(&state).expect("energy");
    let (mut work_j, mut heat_j) = (0.0, 0.0);
    with_cx(&gate, |cx| {
        for step in 0..20_000_u32 {
            let force_n = 0.8 * (f64::from(step) * dt_s * 180.0).sin();
            let inputs = [force_n, 300.0];
            let Proposal::Step(_) = component
                .propose(&state, dt_s, &inputs, &mut next, &mut outputs, cx)
                .expect("step")
            else {
                panic!("the reference never shortens");
            };
            work_j += component.ports()[0].energy_in_j(dt_s, &inputs[..1], &outputs[..1]);
            heat_j -= component.ports()[1].energy_in_j(dt_s, &inputs[1..], &outputs[1..]);
            core::mem::swap(&mut state, &mut next);
        }
    });
    let end_j = component.stored_energy_j(&state).expect("energy");
    let defect_j = work_j - heat_j - (end_j - start_j);
    let scale_j = work_j.abs() + heat_j.abs() + start_j + end_j;
    assert!(heat_j > 0.0);
    assert!(
        defect_j.abs() <= 1.0e-11 * scale_j,
        "defect {defect_j} J of {scale_j} J"
    );
}

#[test]
fn a_heat_leak_fails_the_first_law() {
    let gate = CancelGate::new_clock_free();
    let component = SpringMassDamper::new(Fault::HeatLeak);
    let failure = with_cx(&gate, |cx| check_component(&component, &case(), cx)).expect_err("leak");
    assert!(
        matches!(
            failure,
            Failure::Contract(ContractViolation::EnergyBalance { .. })
        ),
        "{failure}"
    );
}

#[test]
fn reading_the_candidate_slot_fails_purity() {
    let gate = CancelGate::new_clock_free();
    let component = SpringMassDamper::new(Fault::ReadsCandidate);
    let failure =
        with_cx(&gate, |cx| check_component(&component, &case(), cx)).expect_err("impure");
    let Failure::Law { law, .. } = failure else {
        panic!("expected a law failure, got {failure}");
    };
    assert_eq!(law, "propose reads nothing from the candidate slot");
}

#[test]
fn a_component_that_mutates_itself_through_a_cell_fails_purity() {
    let gate = CancelGate::new_clock_free();
    let component = SpringMassDamper::new(Fault::CountsIntoOutput);
    let failure = with_cx(&gate, |cx| check_component(&component, &case(), cx))
        .expect_err("a count that reaches the output");
    let Failure::Law { law, detail } = failure else {
        panic!("expected a law failure, got {failure}");
    };
    let proposals = component.proposals.get();
    eprintln!("refused after {proposals} proposals: {law} ({detail})");
    assert_eq!(
        law,
        "propose is pure: outputs are bitwise equal on repetition"
    );
    // The kit refuses at its repetition: the first two proposals.
    assert_eq!(proposals, 2);
}

#[test]
fn an_interior_cache_that_propose_never_shows_passes_the_kit() {
    let gate = CancelGate::new_clock_free();
    let component = SpringMassDamper::new(Fault::CachesItsLastMean);
    let observed = with_cx(&gate, |cx| check_component(&component, &case(), cx))
        .expect("a cache that propose never shows conforms");
    let hits = component.cache_hits.get();
    let reference = SpringMassDamper::new(Fault::None);
    let expected = with_cx(&gate, |cx| check_component(&reference, &case(), cx)).expect("conforms");
    // Debug prints each float as its shortest round trip, so equal text is
    // equal bits for every value here (none is NaN).
    let (observed, expected) = (format!("{observed:?}"), format!("{expected:?}"));
    eprintln!(
        "cache hits {hits}; observed equals the reference's text: {}",
        observed == expected
    );
    // The kit's repeated proposal read the cache, and saw nothing of it.
    assert!(hits >= 1, "the repeated proposal took the kept mean");
    assert_eq!(observed, expected);
}

#[test]
fn a_wrong_vjp_sign_fails_the_transpose_test() {
    let gate = CancelGate::new_clock_free();
    let component = SpringMassDamper::new(Fault::VjpSign);
    let failure = with_cx(&gate, |cx| {
        check_port_tangent(&component, &case(), 1.0e-7, 1.0e-12, cx)
    })
    .expect_err("wrong sign");
    let Failure::Law { law, .. } = failure else {
        panic!("expected a law failure, got {failure}");
    };
    assert_eq!(law, "vjp is the transpose of jvp");
}

#[test]
fn cancellation_reaches_the_caller_as_the_owner_error() {
    let gate = CancelGate::new_clock_free();
    gate.request();
    let component = SpringMassDamper::new(Fault::None);
    let failure =
        with_cx(&gate, |cx| check_component(&component, &case(), cx)).expect_err("cancelled");
    assert!(matches!(failure, Failure::Owner(Cancelled)), "{failure}");
}

#[test]
fn a_component_that_never_polls_cx_is_refused_at_the_cancellation_law() {
    let gate = CancelGate::new_clock_free();
    let component = SpringMassDamper::new(Fault::IgnoresCancellation);
    let failure = with_cx(&gate, |cx| check_component(&component, &case(), cx))
        .expect_err("a step taken under a cancelled context");
    let Failure::Law { law, detail } = failure else {
        panic!("expected a law failure, got {failure}");
    };
    eprintln!("{law} ({detail})");
    assert_eq!(
        law,
        "propose polls cx: a step asked under a cancelled context is refused"
    );
    assert!(
        detail.starts_with("a context cancelled before the call gave Step("),
        "{detail}"
    );
}

#[test]
fn an_audit_with_a_port_list_that_does_not_fit_the_layout_is_refused() {
    // Too many ports used to index past the layout and panic; too few were
    // silently left out of the energy sum.
    let component = SpringMassDamper::new(Fault::None);
    let ports = component.ports();
    let receipt = Receipt::lossless(0.0, 1.0);
    let short = PortLayout::new(&ports[..1]).expect("layout");
    let long = PortLayout::new(ports).expect("layout");
    for (given, layout, values) in [(ports, &short, 1_usize), (&ports[..1], &long, 2)] {
        let vector = vec![0.0; values];
        let refusal = audit_step(given, layout, 1.0e-3, &vector, &vector, &receipt, 0.0, 0.0);
        assert_eq!(
            refusal,
            Err(ContractViolation::VectorLength {
                vector: "ports",
                expected: layout.ports(),
                actual: given.len(),
            })
        );
    }
}

/// The audit of one step of `dt = 1 s` through the reference's drive port
/// alone, held at 1 N while it moves at `velocity`, under a receipt that
/// claims `stored_change` J and a zero tolerance, with stored energies 0 J
/// before and `stored_new` J after.
fn audit_at_a_zero_tolerance(
    velocity: f64,
    stored_change: f64,
    stored_new: f64,
) -> Result<StepAudit, ContractViolation> {
    let component = SpringMassDamper::new(Fault::None);
    let drive = &component.ports()[..1];
    let layout = PortLayout::new(drive).expect("layout");
    let receipt = Receipt::lossless(stored_change, 0.0);
    audit_step(
        drive,
        &layout,
        1.0,
        &[1.0],
        &[velocity],
        &receipt,
        0.0,
        stored_new,
    )
}

#[test]
fn a_balance_or_an_endpoint_one_subnormal_off_is_admitted_at_a_zero_tolerance() {
    // Below the normal range an f64 rounds by an absolute step, so the audit
    // allows the f64 underflow allowance on top of the receipt's tolerance.
    let subnormal = f64::from_bits(1);
    // 2 subnormals in, 1 stored: the balance is one subnormal off.
    let balance = audit_at_a_zero_tolerance(2.0 * subnormal, subnormal, subnormal)
        .expect("one subnormal off the balance");
    assert_eq!(balance.balance_defect_j.to_bits(), subnormal.to_bits());
    assert_eq!(balance.endpoint_defect_j.to_bits(), 0.0_f64.to_bits());
    // 2 in and 2 stored, but the endpoints measure 1: the endpoint is off.
    let endpoint = audit_at_a_zero_tolerance(2.0 * subnormal, 2.0 * subnormal, subnormal)
        .expect("one subnormal off the endpoints");
    assert_eq!(endpoint.balance_defect_j.to_bits(), 0.0_f64.to_bits());
    assert_eq!(endpoint.endpoint_defect_j.to_bits(), (-subnormal).to_bits());
}

#[test]
fn a_defect_of_1e_300_j_at_a_zero_tolerance_is_still_refused() {
    // The allowance is the smallest normal energy, 2.2e-308 J; a defect of
    // 1e-300 J is far above it and is refused by the law it breaks.
    assert_eq!(
        audit_at_a_zero_tolerance(1.0e-300, 0.0, 0.0),
        Err(ContractViolation::EnergyBalance {
            port_energy_in_j: 1.0e-300,
            stored_change_j: 0.0,
            numerical_loss_j: 0.0,
            tolerance_j: 0.0,
        })
    );
    assert_eq!(
        audit_at_a_zero_tolerance(1.0e-300, 1.0e-300, 0.0),
        Err(ContractViolation::EndpointEnergy {
            measured_change_j: 0.0,
            stored_change_j: 1.0e-300,
            tolerance_j: f64::MIN_POSITIVE,
        })
    );
}

#[test]
fn layout_orders_ports_and_refuses_duplicates() {
    let component = SpringMassDamper::new(Fault::None);
    let layout = PortLayout::new(component.ports()).expect("layout");
    assert_eq!(
        (layout.ports(), layout.input_len(), layout.output_len()),
        (2, 2, 2)
    );
    assert_eq!(layout.input_range(1), 1..2);
    assert_eq!(component.ports()[1].carrier(), Carrier::Heat);
    let twice = [component.ports()[0].clone(), component.ports()[0].clone()];
    assert!(matches!(
        PortLayout::new(&twice),
        Err(ContractViolation::DuplicatePort { .. })
    ));
    let signal = PortSpec::signal_out(StableId::new("probe").expect("id"), NonZeroUsize::MIN);
    let mixed = PortLayout::new(&[component.ports()[0].clone(), signal]).expect("layout");
    assert_eq!((mixed.input_len(), mixed.output_len()), (1, 2));
    assert!(mixed.input_range(1).is_empty());
}

#[test]
fn a_stream_port_books_the_enthalpy_and_mass_that_cross_it_on_either_side() {
    let id = || StableId::new("s").expect("id");
    let storage = PortSpec::stream(id(), PortRole::FlowIn, 2).expect("storage side");
    let transport = PortSpec::stream(id(), PortRole::EffortIn, 2).expect("transport side");
    assert_eq!((storage.width(), storage.carrier()), (4, Carrier::Stream));
    assert_eq!(
        PortSpec::stream(id(), PortRole::FlowIn, 0).expect("bare stream").width(),
        2
    );
    assert!(matches!(
        PortSpec::stream(id(), PortRole::SignalIn, 0),
        Err(ContractViolation::PortRole { .. })
    ));
    let layout = PortLayout::new(core::slice::from_ref(&storage)).expect("layout");
    assert_eq!((layout.input_len(), layout.output_len()), (4, 4));

    // What the storage is like, and what crosses into the transport element
    // (3 g/s out of it here): mass, enthalpy and two carried rates.
    let contents = [2.0e5, 3.0e5, 0.25, 0.75];
    let crossing = [-3.0e-3, -900.0, -7.5e-4, -2.25e-3];
    // Half a second, so the products below are exact.
    let dt = 0.5;
    let same = |got: f64, want: f64| assert!((got - want).abs() <= 0.0, "{got} vs {want}");
    // The storage receives the crossing; the transport element returns it.
    // The energy is the enthalpy rate, not a product of the two sides.
    same(storage.energy_in_j(dt, &crossing, &contents), -450.0);
    same(storage.mass_in_kg(dt, &crossing, &contents), -1.5e-3);
    same(transport.energy_in_j(dt, &contents, &crossing), -450.0);
    same(transport.mass_in_kg(dt, &contents, &crossing), -1.5e-3);
    // A stream's entropy is not on the wire, and only streams carry mass.
    same(storage.entropy_in_j_per_k(dt, &crossing, &contents), 0.0);
    let heat = PortSpec::heat(id(), PortRole::FlowIn).expect("heat");
    same(heat.mass_in_kg(dt, &[5.0], &[300.0]), 0.0);
    assert!(storage.energy_in_j(dt, &crossing[..3], &contents).is_nan());
    assert!(transport.mass_in_kg(dt, &contents, &crossing[..3]).is_nan());
}

/// Unit roundoff of f64.
const U: f64 = f64::EPSILON / 2.0;

/// `n u / (1 - n u)`: a product of `n` factors `(1 + d)`, `|d| <= u`, lies
/// within it of 1.
fn gamma(n: u32) -> f64 {
    let nu = f64::from(n) * U;
    nu / (1.0 - nu)
}

/// The reference without its damper: lossless, and reversible in time.
fn undamped() -> SpringMassDamper {
    let mut component = SpringMassDamper::new(Fault::None);
    component.damping_n_s_m = 0.0;
    component
}

/// Time reversal of the oscillator: the position is even, the velocity odd.
fn reversed(state: &Oscillator) -> Oscillator {
    Oscillator {
        position_m: state.position_m,
        velocity_m_s: -state.velocity_m_s,
    }
}

fn observed(state: &Oscillator) -> Vec<f64> {
    vec![state.position_m, state.velocity_m_s]
}

/// Bounds for the undamped reference's round trip at `case()`, derived from
/// the rounding of `mean_velocity` and `propose`, where every operation is
/// exact times `(1 + d)`, `|d| <= u`. The step is `h = 0.5 dt / m`,
/// `D = 1 + 0.5 dt h k` (computed within `gamma(3) D` of that), and the mean
/// is `fl(N / D)` with `N = fl(v + fl(h fl(F - fl(k x))))`, which is within
/// `gamma(4) (|v| + h |F| + h k |x|)` of its exact value.
/// - Output: in exact arithmetic, the step from `(x1, -v1)` returns minus the
///   forward mean, for any `h`. The rounded `x1` and `v1` are within
///   `gamma(1)` of `x0 + dt mean` and `2 mean - v0`, relative to `dt |mean| +
///   |x1|` and to `|v1|`. With the two numerators' and the two divisions'
///   rounding, the backward mean `b` is within
///   `gamma(8) (|mean| + |b| + (both numerators' sizes + h k dt |mean|) / D)`
///   of minus the forward one.
/// - The velocity of `R(x2)`, `-fl(2 b + v1)`, is within twice that of `v0`,
///   plus `gamma(1) (|v1| + |v2|)`.
/// - The position `fl(x1 + fl(dt b))` is within `dt` times it of `x0`, plus
///   `gamma(1) (dt |mean| + |x|)` for each of the two position updates.
/// - The heat output is `-((0 mean) mean)` both ways: zero, exactly.
/// - The energy `((m / 2) v) v + ((k / 2) x) x` moves by each coordinate's
///   gap times the sum of its two end values, plus `gamma(4)` of both ends'
///   terms.
fn reference_allowance(
    component: &SpringMassDamper,
) -> impl Fn(&Oscillator, &Oscillator, &Oscillator) -> ReversalAllowance + '_ {
    move |start, forward_end, back_end| {
        let (dt_s, force_n) = (1.0e-3, INPUTS[0]);
        let (mass, stiffness) = (component.mass_kg, component.stiffness_n_m);
        let half = 0.5 * dt_s / mass;
        let denominator = 1.0 + 0.5 * dt_s * half * stiffness;
        let mean = component.mean_velocity(start, dt_s, force_n).0.abs();
        let mean_back = component
            .mean_velocity(&reversed(forward_end), dt_s, force_n)
            .0
            .abs();
        let [x0, x1, x2] = [start, forward_end, back_end].map(|s| s.position_m.abs());
        let [v0, v1, v2] = [start, forward_end, back_end].map(|s| s.velocity_m_s.abs());
        let sizes =
            v0 + v1 + 2.0 * half * force_n.abs() + half * stiffness * (x0 + x1 + dt_s * mean);
        let output = gamma(8) * (mean + mean_back + sizes / denominator);
        let velocity = 2.0 * output + gamma(1) * (v1 + v2);
        let position = dt_s * output + gamma(1) * (dt_s * (mean + mean_back) + x1 + x2);
        let terms = 0.5 * mass * (v0 * v0 + v2 * v2) + 0.5 * stiffness * (x0 * x0 + x2 * x2);
        ReversalAllowance {
            outputs: vec![output, 0.0],
            observed: vec![position, velocity],
            energy_j: 0.5 * mass * velocity * (v0 + v2)
                + 0.5 * stiffness * position * (x0 + x2)
                + gamma(4) * terms,
        }
    }
}

#[test]
fn the_undamped_reference_reverses_in_time() {
    let gate = CancelGate::new_clock_free();
    let component = undamped();
    let allowance = reference_allowance(&component);
    let reversal = with_cx(&gate, |cx| {
        check_time_reversal(&component, &case(), &reversed, &observed, &allowance, cx)
    })
    .expect("the undamped reference reverses");
    eprintln!("{reversal:?}");
    // The heat port returns zero heat both ways.
    assert_eq!(reversal.output_gaps[1].to_bits(), 0);
}

/// The law `check_time_reversal` refuses `component` at under `reverse`.
fn refused_at(
    component: &SpringMassDamper,
    reverse: &dyn Fn(&Oscillator) -> Oscillator,
) -> &'static str {
    let gate = CancelGate::new_clock_free();
    let allowance = reference_allowance(component);
    let failure = with_cx(&gate, |cx| {
        check_time_reversal(component, &case(), reverse, &observed, &allowance, cx)
    })
    .expect_err("refused");
    let Failure::Law { law, detail } = failure else {
        panic!("expected a law failure, got {failure}");
    };
    eprintln!("refused: {law} ({detail})");
    law
}

#[test]
fn the_damped_reference_is_refused_reversal_because_it_dissipates() {
    let law = refused_at(&SpringMassDamper::new(Fault::None), &reversed);
    assert_eq!(
        law,
        "a reversible step dissipates nothing and generates no entropy"
    );
}

#[test]
fn a_reversal_that_keeps_the_velocity_fails_the_outputs_law() {
    let law = refused_at(&undamped(), &|state: &Oscillator| *state);
    assert_eq!(
        law,
        "the backward outputs are the forward outputs times their parity"
    );
}

#[test]
fn a_reversal_that_is_not_an_involution_is_refused() {
    let doubled = |state: &Oscillator| Oscillator {
        position_m: state.position_m,
        velocity_m_s: -2.0 * state.velocity_m_s,
    };
    let law = refused_at(&undamped(), &doubled);
    assert_eq!(law, "time reversal is an involution: R(R(x)) == x");
}
