//! Owner slice: the linear modal structure of fs-nlmodal on the composition kernel.
//!
//! `fs_nlmodal::component::ModalStructure` puts fs-modal's exact held-force step behind the
//! Component contract: attachments through which efforts act on the modes, and a heat port for
//! the modal damping. Here it passes the conformance check, is shown to add nothing to the
//! owner's step, and is run by the kernel through a strike and a free decay whose heat goes to a
//! sink, with the ledger closing on the owner's work and heat. Beside a load line, a 64-mode string
//! added with its tangent hands the line its own mobilities, and joins a line only where the line
//! is told how its step ends.
//! Run with `cargo test -p fs-system --features modal-slices --test modal_structure_slice`.

// Exact float comparisons here assert bitwise agreement with the owner.
#![allow(clippy::float_cmp)]

use core::num::NonZeroUsize;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_nlmodal::component::{
    ModalAcousticState, ModalAcousticTimeBudget, ModalAttachment, ModalStructure,
    ModalStructureRefusal,
};
use fs_ports::component::conformance::{Case, Observed, check_component};
use fs_ports::component::{Component, Proposal};
use fs_ports::{PortKind, StableId};
use fs_system::{RouteTemperature, SystemBuilder, heat_sink};

const TAU: f64 = core::f64::consts::TAU;
const RATE_HZ: u32 = 44_100;
/// The first three vertical modes of a bar [Hz] and one loss factor for all of them; any
/// positive values serve, these are of the size the vibraphone cards carry.
const MODES_HZ: [f64; 3] = [331.125_927_8, 1_324.505_956, 2_972.946_707];
const LOSS_FACTOR: f64 = 8.34e-5;
/// Mode shapes at the struck point [1/sqrt(kg)].
const STRIKE: [f64; 3] = [2.1, -1.4, 0.9];
const AIR_K: f64 = 293.15;

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

fn modes() -> Vec<(f64, f64)> {
    // A loss factor eta is a viscous damping ratio of eta / 2 at each mode's own frequency.
    MODES_HZ
        .iter()
        .map(|&f| (TAU * f, 0.5 * LOSS_FACTOR))
        .collect()
}

fn attachments() -> Vec<ModalAttachment> {
    let mut identity = vec![0.0; 9];
    for k in 0..3 {
        identity[4 * k] = 1.0;
    }
    vec![
        ModalAttachment {
            id: StableId::new("strike").expect("id"),
            kind: PortKind::MechanicalForceVelocity,
            width: NonZeroUsize::MIN,
            weights: STRIKE.to_vec(),
        },
        ModalAttachment {
            id: StableId::new("modal").expect("id"),
            kind: PortKind::MechanicalForceVelocity,
            width: NonZeroUsize::new(3).expect("three"),
            weights: identity,
        },
    ]
}

fn structure() -> ModalStructure {
    ModalStructure::from_modes(
        RATE_HZ,
        &modes(),
        ModalAcousticTimeBudget::audible_reference(),
        attachments(),
    )
    .expect("structure")
}

fn moving() -> Vec<ModalAcousticState> {
    [(1.0e-6, 2.0e-3), (-3.0e-7, 1.0e-3), (2.0e-8, -4.0e-4)]
        .iter()
        .map(|&(q, v)| ModalAcousticState {
            displacement_m_sqrt_kg: q,
            velocity_m_sqrt_kg_per_s: v,
        })
        .collect()
}

/// The generalized forces the structure forms from its rows, in its own order of operations.
fn generalized(inputs: &[f64]) -> Vec<f64> {
    let attachments = attachments();
    let weights: Vec<f64> = attachments.iter().flat_map(|a| a.weights.clone()).collect();
    let mut force = vec![0.0; 3];
    for (row, &effort) in inputs[..4].iter().enumerate() {
        for k in 0..3 {
            force[k] += weights[row * 3 + k] * effort;
        }
    }
    force
}

/// The owner's own step from `start` under the same generalized forces.
fn owner_step(
    component: &ModalStructure,
    start: &[ModalAcousticState],
    force: &[f64],
    dt_s: f64,
) -> (
    Vec<ModalAcousticState>,
    fs_nlmodal::component::ModalAcousticTimeModel,
    f64,
    f64,
    f64,
) {
    let mut model = component.model().clone();
    model.restore_states(start).expect("restore");
    let frame = model.step_duration(force, dt_s).expect("owner step");
    (
        model.states().to_vec(),
        model,
        frame.input_work_j,
        frame.viscous_dissipation_j,
        frame.total_modal_energy_j,
    )
}

#[test]
fn the_modal_structure_conforms_and_books_the_owners_work_and_heat() {
    let gate = CancelGate::new_clock_free();
    let component = structure();
    assert_eq!((component.modes(), component.ports().len()), (3, 3));
    let start = moving();
    let stale = component.at_rest();
    let sample = 1.0 / f64::from(RATE_HZ);
    for (dt_s, scale) in [(sample, 1.0), (0.5 * sample, -2.0), (3.7e-4, 0.5)] {
        let inputs = [
            0.8 * scale,
            0.01 * scale,
            -0.02 * scale,
            0.005 * scale,
            AIR_K,
        ];
        let case = Case {
            state: &start,
            stale: &stale,
            dt_s,
            inputs: &inputs,
        };
        let observed =
            with_cx(&gate, |cx| check_component(&component, &case, cx)).expect("conforms");
        let Observed::Stepped {
            receipt,
            audit,
            outputs,
        } = observed
        else {
            panic!("the structure takes the step");
        };
        let force = generalized(&inputs);
        let (_, _, work_j, heat_j, _) = owner_step(&component, &start, &force, dt_s);
        assert!(heat_j > 0.0, "the modes dissipate");
        assert!(work_j.abs() > 0.0, "the efforts do work");
        let paired_j: f64 = (0..4)
            .map(|row| dt_s * (inputs[row] * outputs[row]).abs())
            .sum();
        let exchanged_j = paired_j + heat_j;
        assert!(
            (audit.port_energy_in_j - (work_j - heat_j)).abs() <= 1.0e-12 * exchanged_j,
            "{} vs {} - {}",
            audit.port_energy_in_j,
            work_j,
            heat_j
        );
        assert_eq!(outputs[4], -heat_j / dt_s);
        assert_eq!(receipt.dissipated_j, heat_j);
        assert_eq!(receipt.entropy_generated_j_per_k, heat_j / AIR_K);
        assert!(
            receipt.numerical_loss_j.abs()
                <= 1.0e-12 * (exchanged_j + receipt.stored_change_j.abs()),
            "{}",
            receipt.numerical_loss_j
        );
    }
}

#[test]
fn the_wrapper_adds_nothing_to_the_owners_step() {
    let gate = CancelGate::new_clock_free();
    let component = structure();
    let start = moving();
    let inputs = [-0.3, 0.0, 0.04, -0.01, 300.0];
    let dt_s = 1.0 / f64::from(RATE_HZ);
    let mut next = Vec::new();
    let mut outputs = [f64::NAN; 5];
    let proposal = with_cx(&gate, |cx| {
        component.propose(&start, dt_s, &inputs, &mut next, &mut outputs, cx)
    })
    .expect("step");
    let Proposal::Step(receipt) = proposal else {
        panic!("the structure takes the step");
    };
    let force = generalized(&inputs);
    let (states, model, _, heat_j, total_j) = owner_step(&component, &start, &force, dt_s);
    assert_eq!(next, states);
    let moved: Vec<f64> = (0..3)
        .map(|k| next[k].displacement_m_sqrt_kg - start[k].displacement_m_sqrt_kg)
        .collect();
    let strike: f64 = (0..3).map(|k| STRIKE[k] * moved[k]).sum();
    assert_eq!(outputs[0], strike / dt_s);
    for k in 0..3 {
        // An identity row is a sum whose other terms are exact zeros.
        assert_eq!(outputs[1 + k], moved[k] / dt_s);
    }
    assert_eq!(outputs[4], -heat_j / dt_s);
    // The stored energy is the owner's own total, bit for bit.
    assert_eq!(model.energy_j(&next).expect("energy"), total_j);
    assert_eq!(component.stored_energy_j(&next).expect("energy"), total_j);
    assert_eq!(
        receipt.stored_change_j,
        total_j - component.stored_energy_j(&start).expect("energy")
    );
}

fn account<'a>(system: &'a fs_system::System, name: &str) -> &'a fs_system::ledger::NodeAccount {
    system
        .ledger()
        .nodes
        .iter()
        .find(|node| node.name == name)
        .expect("account")
}

/// The structure in a system: the strike port prescribed by `strike(t0)`, the modal port held at
/// zero, the heat sent to a sink at the air's temperature.
fn graph(
    start: Vec<ModalAcousticState>,
    strike: fn(f64) -> f64,
) -> (fs_system::System, fs_system::NodeId) {
    let mut builder = SystemBuilder::new();
    let node = builder.add("bar", structure(), start).expect("bar");
    let port = |builder: &SystemBuilder, name: &str| builder.port(node, name).expect("port");
    builder
        .boundary(
            port(&builder, "strike"),
            Box::new(move |t0, _dt, out: &mut [f64]| out[0] = strike(t0)),
        )
        .expect("strike boundary");
    builder
        .boundary(
            port(&builder, "modal"),
            Box::new(|_t0, _dt, out: &mut [f64]| out.fill(0.0)),
        )
        .expect("modal boundary");
    builder
        .heat_route(
            port(&builder, "heat"),
            &[],
            Some(heat_sink("air", AIR_K).expect("sink")),
            RouteTemperature::Fixed(AIR_K),
        )
        .expect("heat route");
    (builder.seal().expect("sealed"), node)
}

#[test]
fn a_free_decay_sends_its_heat_to_the_sink_at_the_loss_rows_rate() {
    let gate = CancelGate::new_clock_free();
    let start = moving();
    let (mut system, node) = graph(start.clone(), |_| 0.0);
    let component = structure();
    let dt_s = 1.0 / f64::from(RATE_HZ);
    let steps = 4_410_usize;
    let mut twin = start.clone();
    let mut heat_j = 0.0_f64;
    for _ in 0..steps {
        let (next, _, work_j, heat, _) = owner_step(&component, &twin, &[0.0; 3], dt_s);
        assert_eq!(work_j, 0.0);
        heat_j += heat;
        twin = next;
        with_cx(&gate, |cx| system.step(dt_s, cx)).expect("kernel step");
    }
    let reached = system
        .state::<Vec<ModalAcousticState>>(node)
        .expect("bar state");
    assert_eq!(reached, &twin);

    let before_j = component.stored_energy_j(&start).expect("energy");
    let after_j = component.stored_energy_j(&twin).expect("energy");
    assert!(
        heat_j > 1.0e-3 * before_j,
        "enough heat to measure: {heat_j} of {before_j} J"
    );
    assert!((before_j - after_j - heat_j).abs() <= 1.0e-10 * before_j);
    let ledger = system.ledger();
    assert_eq!(ledger.steps, steps as u64);
    assert!((ledger.sink_heat_j() - heat_j).abs() <= 1.0e-11 * heat_j);
    let bar = account(&system, "bar");
    assert!((bar.dissipated_j.value() - heat_j).abs() <= 1.0e-11 * heat_j);
    assert!((bar.stored_change_j.value() + heat_j).abs() <= 1.0e-10 * before_j);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-11 * before_j);

    // Each mode's energy decays at twice its damping rate, zeta * omega, up to the ripple of a
    // damped oscillator's energy within a cycle, at most zeta of it at each end (so 2 zeta / (1 - zeta)).
    let t = dt_s * steps as f64;
    for (k, &(omega, zeta)) in modes().iter().enumerate() {
        let energy = |s: &ModalAcousticState| {
            0.5 * (s.velocity_m_sqrt_kg_per_s * s.velocity_m_sqrt_kg_per_s
                + omega * omega * s.displacement_m_sqrt_kg * s.displacement_m_sqrt_kg)
        };
        let ratio = energy(&twin[k]) / energy(&start[k]);
        let law = (-2.0 * zeta * omega * t).exp();
        assert!(
            (ratio / law - 1.0).abs() <= 2.0 * zeta / (1.0 - zeta) + 1.0e-12,
            "mode {k}: {ratio} vs {law}"
        );
    }
}

fn pulse(t0: f64) -> f64 {
    // A half-sine force of 1 ms and 5 N peak, held at its value at each step's start.
    let width_s = 1.0e-3;
    if t0 < width_s {
        5.0 * (core::f64::consts::PI * t0 / width_s).sin()
    } else {
        0.0
    }
}

#[test]
fn a_strike_through_the_kernel_enters_as_boundary_work() {
    let gate = CancelGate::new_clock_free();
    let start = vec![ModalAcousticState::default(); 3];
    let (mut system, node) = graph(start.clone(), pulse);
    let component = structure();
    let dt_s = 1.0 / f64::from(RATE_HZ);
    let steps = 2_205_usize;
    let mut twin = start;
    let (mut work_j, mut heat_j) = (0.0_f64, 0.0_f64);
    let mut t0 = 0.0_f64;
    for step in 0..steps {
        let held = pulse(t0);
        let force: Vec<f64> = STRIKE.iter().map(|&s| s * held + 0.0).collect();
        let (next, _, work, heat, _) = owner_step(&component, &twin, &force, dt_s);
        work_j += work;
        heat_j += heat;
        twin = next;
        with_cx(&gate, |cx| system.step(dt_s, cx)).expect("kernel step");
        t0 = system.time_s();
        assert_eq!(system.steps(), step as u64 + 1);
    }
    let reached = system
        .state::<Vec<ModalAcousticState>>(node)
        .expect("bar state");
    assert_eq!(reached, &twin);
    let stored_j = component.stored_energy_j(&twin).expect("energy");
    assert!(work_j > 0.0 && stored_j > 0.0);
    let ledger = system.ledger();
    let exchanged_j = work_j.abs() + heat_j;
    assert!(
        (ledger.boundary_energy_in_j() - work_j).abs() <= 1.0e-11 * exchanged_j,
        "{} vs {}",
        ledger.boundary_energy_in_j(),
        work_j
    );
    assert!((ledger.sink_heat_j() - heat_j).abs() <= 1.0e-11 * exchanged_j);
    assert!((work_j - heat_j - stored_j).abs() <= 1.0e-10 * exchanged_j);
    assert!(ledger.energy_residual_j().abs() <= 1.0e-11 * exchanged_j);
}

#[test]
fn the_structure_refuses_what_it_cannot_represent() {
    let gate = CancelGate::new_clock_free();
    let budget = ModalAcousticTimeBudget::audible_reference();
    let mut short = attachments();
    short[0].weights.pop();
    assert_eq!(
        ModalStructure::from_modes(RATE_HZ, &modes(), budget, short).unwrap_err(),
        ModalStructureRefusal::Weights {
            attachment: 0,
            expected: 3,
            found: 2
        }
    );
    let mut bad = attachments();
    bad[1].weights[4] = f64::NAN;
    assert_eq!(
        ModalStructure::from_modes(RATE_HZ, &modes(), budget, bad).unwrap_err(),
        ModalStructureRefusal::NonFiniteWeight { attachment: 1 }
    );
    let mut clash = attachments();
    clash[0].id = StableId::new("heat").expect("id");
    assert!(matches!(
        ModalStructure::from_modes(RATE_HZ, &modes(), budget, clash).unwrap_err(),
        ModalStructureRefusal::Port(_)
    ));
    // A model with an observer transfer is refused in fs-nlmodal's own unit test, which can build
    // the complex transfer.

    let component = structure();
    let start = moving();
    let mut next = start.clone();
    let mut outputs = [0.0; 5];
    let mut inputs = [0.1, 0.0, 0.0, 0.0, 0.0];
    let refused = with_cx(&gate, |cx| {
        component.propose(&start, 1.0e-4, &inputs, &mut next, &mut outputs, cx)
    })
    .unwrap_err();
    assert_eq!(refused, ModalStructureRefusal::Temperature { kelvin: 0.0 });
    inputs[4] = AIR_K;
    let two = vec![start[0], start[1]];
    let refused = with_cx(&gate, |cx| {
        component.propose(&two, 1.0e-4, &inputs, &mut next, &mut outputs, cx)
    })
    .unwrap_err();
    assert_eq!(
        refused,
        ModalStructureRefusal::States {
            expected: 3,
            found: 2
        }
    );
    assert!(component.stored_energy_j(&two).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    let refused = with_cx(&cancelled, |cx| {
        component.propose(&start, 1.0e-4, &inputs, &mut next, &mut outputs, cx)
    })
    .unwrap_err();
    assert_eq!(refused, ModalStructureRefusal::Cancelled);
}

// A string beside a load line (bead rc-1109 .12.5, commit 2): the line it is handed, and how
// its step ends.

/// A load-line component of one port, `a`, that takes the end of the step on the signal input
/// `a.end` or takes none: it returns no effort and keeps what it is handed.
struct LineProbe {
    ports: Vec<fs_ports::component::PortSpec>,
    takes_end: bool,
}

impl LineProbe {
    fn new(takes_end: bool) -> Self {
        use fs_ports::component::{PortRole, PortSpec};
        let mut ports = vec![
            PortSpec::power(
                StableId::new("a").expect("id"),
                PortKind::MechanicalForceVelocity,
                PortRole::FlowIn,
                NonZeroUsize::MIN,
            )
            .expect("port"),
        ];
        if takes_end {
            ports.push(PortSpec::signal_in(
                StableId::new("a.end").expect("id"),
                NonZeroUsize::MIN,
            ));
        }
        Self { ports, takes_end }
    }
}

impl Component for LineProbe {
    type State = Vec<f64>;
    // It refuses nothing (`fs_system::elements` is behind the `constitutive` feature).
    type Error = core::convert::Infallible;

    fn ports(&self) -> &[fs_ports::component::PortSpec] {
        &self.ports
    }

    fn class(&self) -> fs_ports::component::Class {
        fs_ports::component::Class::NonSmooth
    }

    fn stored_energy_j(&self, _handed: &Vec<f64>) -> Result<f64, Self::Error> {
        Ok(0.0)
    }

    fn propose(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Vec<f64>,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        let zero = vec![0.0; if self.takes_end { 2 } else { 1 }];
        fs_ports::component::LoadLine::propose_on_line(
            self, old, dt_s, inputs, &zero, new, outputs, cx,
        )
    }
}

impl fs_ports::component::LoadLine for LineProbe {
    fn line_ports(&self) -> &[usize] {
        &[0]
    }

    fn line_ends(&self) -> &[Option<usize>] {
        if self.takes_end { &[Some(1)] } else { &[] }
    }

    fn propose_on_line(
        &self,
        _old: &Vec<f64>,
        _dt_s: f64,
        _inputs: &[f64],
        compliance: &[f64],
        new: &mut Vec<f64>,
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, Self::Error> {
        new.clear();
        new.extend_from_slice(compliance);
        outputs.fill(0.0);
        Ok(Proposal::Step(fs_ports::component::Receipt::lossless(
            0.0,
            f64::MIN_POSITIVE,
        )))
    }
}

const STRING_RATE_HZ: u32 = 48_000;
const STRING_MODES: usize = 64;
const STRING_HZ: f64 = 196.0;
const STRING_RATIO: f64 = 2.0e-4;
/// A string of 1 g bowed at a tenth of its length: mode shapes `sqrt(2 / m) sin(n pi x / L)`.
const STRING_KG: f64 = 1.0e-3;
const BOW_AT: f64 = 0.1;

fn string_modes() -> Vec<(f64, f64)> {
    (1..=STRING_MODES)
        .map(|n| (TAU * STRING_HZ * n as f64, STRING_RATIO))
        .collect()
}

fn bow_shapes() -> Vec<f64> {
    (1..=STRING_MODES)
        .map(|n| (2.0 / STRING_KG).sqrt() * (core::f64::consts::PI * n as f64 * BOW_AT).sin())
        .collect()
}

fn string(telling: bool) -> ModalStructure {
    let structure = ModalStructure::from_modes(
        STRING_RATE_HZ,
        &string_modes(),
        ModalAcousticTimeBudget::audible_reference(),
        vec![ModalAttachment {
            id: StableId::new("bow").expect("id"),
            kind: PortKind::MechanicalForceVelocity,
            width: NonZeroUsize::MIN,
            weights: bow_shapes(),
        }],
    )
    .expect("string");
    if telling {
        structure.publishing_ends().expect("ends")
    } else {
        structure
    }
}

fn ringing() -> Vec<ModalAcousticState> {
    (1..=STRING_MODES)
        .map(|n| ModalAcousticState {
            displacement_m_sqrt_kg: 1.0e-6 / n as f64,
            velocity_m_sqrt_kg_per_s: 1.0e-3 / (n * n) as f64,
        })
        .collect()
}

/// How the string joins the line.
#[derive(Clone, Copy, PartialEq)]
enum Joined {
    /// Added with its tangent.
    Declared,
    /// Added without it: the kernel differences it.
    Differenced,
}

/// A system of the string beside a probe on a line, sealed; its heat goes to a sink.
fn beside(
    telling: bool,
    joined: Joined,
    takes_end: bool,
) -> Result<(fs_system::System, fs_system::NodeId), fs_system::AdmissionError> {
    let mut builder = SystemBuilder::new();
    let structure = string(telling);
    let node = match joined {
        Joined::Declared => builder.add_with_tangent("string", structure, ringing()),
        Joined::Differenced => builder.add("string", structure, ringing()),
    }
    .expect("string");
    let probe = builder
        .add_on_line("probe", LineProbe::new(takes_end), Vec::new())
        .expect("probe");
    let bow = builder.port(node, "bow").expect("bow");
    let a = builder.port(probe, "a").expect("a");
    builder.common_flow(bow, &[a]).expect("bow on the line");
    builder
        .heat_route(
            builder.port(node, "heat").expect("heat"),
            &[],
            Some(heat_sink("air", AIR_K).expect("sink")),
            RouteTemperature::Fixed(AIR_K),
        )
        .expect("heat route");
    builder.seal().map(|system| (system, probe))
}

#[test]
fn a_line_beside_a_string_with_its_tangent_is_handed_the_strings_own_mobilities() {
    let dt_s = 1.0 / f64::from(STRING_RATE_HZ);
    let structure = string(true);
    let mut gains = vec![ModalAcousticState::default(); STRING_MODES];
    structure
        .model()
        .unit_force_response_into(dt_s, &mut gains)
        .expect("gains");
    // The line's two entries from the owner's gains: the mean flow's and the end's mobility.
    let (mut mean, mut end) = (0.0_f64, 0.0_f64);
    for (shape, gain) in bow_shapes().iter().zip(&gains) {
        mean += shape * shape * gain.displacement_m_sqrt_kg;
        end += shape * shape * gain.velocity_m_sqrt_kg_per_s;
    }
    mean /= dt_s;
    // The gains against the closed form of a damped oscillator's response from rest to a unit
    // force: an independent check that they are the modes' own.
    let mut closed_worst = 0.0_f64;
    for (&(omega, zeta), gain) in string_modes().iter().zip(&gains) {
        let rate = zeta * omega;
        let damped = omega * (1.0 - zeta * zeta).sqrt();
        let decay = (-rate * dt_s).exp();
        let (sine, cosine) = (damped * dt_s).sin_cos();
        let velocity = decay * sine / damped;
        let displacement = (1.0 - decay * (cosine + rate / damped * sine)) / (omega * omega);
        closed_worst = closed_worst
            .max(((gain.velocity_m_sqrt_kg_per_s - velocity) / velocity).abs())
            .max(((gain.displacement_m_sqrt_kg - displacement) / displacement).abs());
    }
    eprintln!("the owner's gains against the closed form: worst relative {closed_worst:e}");
    assert!(closed_worst <= 1.0e-9, "{closed_worst:e}");

    for joined in [Joined::Declared, Joined::Differenced] {
        let (mut system, probe) = beside(true, joined, true).expect("sealed");
        let report =
            with_cx(&CancelGate::new_clock_free(), |cx| system.step(dt_s, cx)).expect("step");
        let handed = system.state::<Vec<f64>>(probe).expect("probe").clone();
        let (off_mean, off_end) = ((handed[0] - mean) / mean, (handed[1] - end) / end);
        eprintln!(
            "{}: handed {:e} and {:e} for {mean:e} and {end:e}: relative {off_mean:e} and {off_end:e}, {} unresolved",
            if joined == Joined::Declared {
                "declared"
            } else {
                "differenced"
            },
            handed[0],
            handed[1],
            report.unresolved
        );
        if joined == Joined::Declared {
            assert!(off_mean.abs() <= 1.0e-14 && off_end.abs() <= 1.0e-14);
            assert_eq!(report.unresolved, 0);
        }
    }
}

#[test]
fn a_string_joins_a_line_only_where_the_line_is_told_how_its_step_ends() {
    let refused =
        |result: Result<(fs_system::System, fs_system::NodeId), fs_system::AdmissionError>,
         says: &str| match result {
            Err(fs_system::AdmissionError::LineEnd { port, reason }) => {
                assert!(reason.contains(says), "{port}: {reason}");
                assert!(port.contains("probe"), "{port}");
            }
            Err(other) => panic!("{other}"),
            Ok(_) => panic!("sealed: expected a refusal that says `{says}`"),
        };
    assert!(
        beside(true, Joined::Declared, true).is_ok(),
        "told, beside a port that takes it"
    );
    refused(
        beside(true, Joined::Declared, false),
        "the port takes no end input",
    );
    refused(
        beside(false, Joined::Declared, true),
        "does not tell its end",
    );
    refused(
        beside(false, Joined::Declared, false),
        "does not tell its end",
    );
}
