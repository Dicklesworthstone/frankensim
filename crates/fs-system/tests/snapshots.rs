//! Durable snapshots (bead rc-1109 .5.4, item 4): a system sealed in an
//! fs-exec v2 envelope resumes in another process on the trajectory it would
//! have followed, forks into systems that go on alone, and is refused, with
//! nothing changed, by a system of another problem.

use std::fmt::Write as _;
use std::io::{Read as _, Write as _};
use std::marker::PhantomData;
use std::process::{Command, Stdio};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::identity::{
    Admitted, AuthorityAdmitter, AuthorityRef, AuthorityVerifier, CanonicalLimits, CanonicalSchema,
    ContentId, ExternalAnchorRef, FieldSpec, IdentityReceipt, KeyPolicyId, Presented,
    StrongIdentity, Verified, VerifierId, WireType,
};
use fs_exec::solver::SolverStateV2;
use fs_exec::solver::snapshot_v2::{
    self, ExpectedResumeContextV2, PausedSnapshotBoundaryV2, SnapshotAlgorithmIdV2,
    SnapshotAuthoritySubjectIdV2, SnapshotBudgetStateIdV2, SnapshotDeterminismV2,
    SnapshotExecutionFingerprintIdV2, SnapshotLimitsV2, SnapshotPauseRequestIdV2,
    SnapshotProvenanceIdV2, SnapshotRngCounterIdV2, StateIdentityCharterV2,
};
use fs_exec::{Budget, CancelGate, Cx, DrainTracker, ExecMode, RunId, StreamKey};
use fs_ports::component::{Class, Component, EndOfStep, PortSpec, PortTangent, Proposal, Receipt};
use fs_system::elements::{
    Damper, Force, GasContents, GasVolume, HeatCapacity, Mass, Motion, Restriction, SourceAccount,
    Spring, Waveform,
};
use fs_system::ledger::Sum;
use fs_system::{NodeId, PortRef, RouteTemperature, System, SystemBuilder, SystemImage, heat_sink};

/// The variable that makes a run of this binary the resumed child.
const CHILD: &str = "FS_SYSTEM_SNAPSHOT_CHILD";

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let gate = CancelGate::new_clock_free();
        let key = StreamKey {
            seed: 13,
            kernel_id: 37,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

fn limits() -> SnapshotLimitsV2 {
    SnapshotLimitsV2::new(
        1 << 22,
        4_096,
        CanonicalLimits::new(16_384, 4_096, 32, 32, 64),
        1 << 16,
        1 << 22,
        4_096,
    )
}

/// The pause boundary of a run whose one worker drained, as the executor
/// mints it. The same run and counts give the same boundary in any process.
fn paused() -> PausedSnapshotBoundaryV2 {
    let gate = CancelGate::new_clock_free();
    let tracker = DrainTracker::new(RunId(17), &gate);
    let worker = tracker.register_worker().expect("worker registers");
    gate.request();
    drop(worker);
    let report = tracker.finalize().expect("run drained");
    PausedSnapshotBoundaryV2::from_drain_report(
        report,
        SnapshotPauseRequestIdV2::from_bytes([0x5A; 32]),
        1,
    )
}

/// The context a snapshot of `system` is sealed and opened under.
fn context(system: &System) -> ExpectedResumeContextV2 {
    ExpectedResumeContextV2::for_paused_state::<SystemImage>(
        SnapshotAlgorithmIdV2::from_bytes([0xA1; 32]),
        1,
        system.problem_id(),
        SnapshotRngCounterIdV2::from_bytes([0; 32]),
        SnapshotDeterminismV2::Deterministic,
        SnapshotExecutionFingerprintIdV2::from_bytes([0xE1; 32]),
        SnapshotBudgetStateIdV2::from_bytes([0xB1; 32]),
        SnapshotProvenanceIdV2::from_bytes([0xC1; 32]),
        paused(),
    )
}

/// Open by the root the sealer kept: the route inside one process.
fn opened_here(
    sealed: &snapshot_v2::SealedSnapshotV2,
) -> snapshot_v2::OpenedSnapshotV2<SystemImage> {
    SystemImage::unseal_v2_expected(sealed.bytes(), &sealed.expectation(), limits(), || false)
        .expect("the kept root opens the snapshot")
}

// ---------------------------------------------------------------------------
// Codecs the test declares for kernel elements and its own members
// ---------------------------------------------------------------------------

/// A state the test writes as plain words.
trait Words: Sized {
    const CHARTER: StateIdentityCharterV2;
    fn put(&self, words: &mut Vec<u64>);
    fn get(words: &[u64]) -> Option<Self>;
}

impl Words for Motion {
    const CHARTER: StateIdentityCharterV2 = StateIdentityCharterV2 {
        owner: "fs-system::tests::snapshots",
        state_family: "motion",
        schema_grammar: "position_m f64; velocity_m_s f64",
        codec_grammar: "two u64 words, the bits of position_m then velocity_m_s",
        codec_version: 1,
    };

    fn put(&self, words: &mut Vec<u64>) {
        words.extend([self.position_m.to_bits(), self.velocity_m_s.to_bits()]);
    }

    fn get(words: &[u64]) -> Option<Self> {
        match words {
            &[position, velocity] => Some(Self {
                position_m: f64::from_bits(position),
                velocity_m_s: f64::from_bits(velocity),
            }),
            _ => None,
        }
    }
}

impl Words for f64 {
    const CHARTER: StateIdentityCharterV2 = StateIdentityCharterV2 {
        owner: "fs-system::tests::snapshots",
        state_family: "one-f64",
        schema_grammar: "one f64",
        codec_grammar: "one u64 word, its bits",
        codec_version: 1,
    };

    fn put(&self, words: &mut Vec<u64>) {
        words.push(self.to_bits());
    }

    fn get(words: &[u64]) -> Option<Self> {
        match words {
            &[bits] => Some(f64::from_bits(bits)),
            _ => None,
        }
    }
}

impl Words for () {
    const CHARTER: StateIdentityCharterV2 = StateIdentityCharterV2 {
        owner: "fs-system::tests::snapshots",
        state_family: "no-state",
        schema_grammar: "()",
        codec_grammar: "no words",
        codec_version: 1,
    };

    fn put(&self, _words: &mut Vec<u64>) {}

    fn get(words: &[u64]) -> Option<Self> {
        words.is_empty().then_some(())
    }
}

impl Words for SourceAccount {
    const CHARTER: StateIdentityCharterV2 = StateIdentityCharterV2 {
        owner: "fs-system::tests::snapshots",
        state_family: "source-account",
        schema_grammar: "clock_s: a compensated sum, its running value and compensation; \
            energy_in_j f64",
        codec_grammar: "three u64 words: Sum::to_words of clock_s, then the bits of energy_in_j",
        codec_version: 1,
    };

    fn put(&self, words: &mut Vec<u64>) {
        words.extend(self.clock_s.to_words());
        words.push(self.energy_in_j.to_bits());
    }

    fn get(words: &[u64]) -> Option<Self> {
        match words {
            &[sum, compensation, energy] => Some(Self {
                clock_s: Sum::from_words([sum, compensation]),
                energy_in_j: f64::from_bits(energy),
            }),
            _ => None,
        }
    }
}

/// A kernel element with a codec for its state and the bits of its
/// parameters as its identity; every call is the element's own.
struct Coded<C> {
    element: C,
    identity: Vec<u64>,
}

impl<C> Coded<C> {
    fn new(element: C, parameters: &[f64]) -> Self {
        Self {
            element,
            identity: parameters.iter().map(|value| value.to_bits()).collect(),
        }
    }
}

impl<C> Component for Coded<C>
where
    C: Component,
    C::State: Words,
{
    type State = C::State;
    type Error = C::Error;

    fn ports(&self) -> &[PortSpec] {
        self.element.ports()
    }

    fn class(&self) -> Class {
        self.element.class()
    }

    fn stored_energy_j(&self, state: &C::State) -> Result<f64, C::Error> {
        self.element.stored_energy_j(state)
    }

    fn stored_mass_kg(&self, state: &C::State) -> Result<Option<f64>, C::Error> {
        self.element.stored_mass_kg(state)
    }

    fn stored_carried_quantity(
        &self,
        state: &C::State,
        carried: usize,
    ) -> Result<Option<f64>, C::Error> {
        self.element.stored_carried_quantity(state, carried)
    }

    fn limit_dt_s(&self, state: &C::State, proposed_s: f64) -> Result<f64, C::Error> {
        self.element.limit_dt_s(state, proposed_s)
    }

    fn next_event_s(&self, state: &C::State) -> Option<f64> {
        self.element.next_event_s(state)
    }

    fn end_of_step(&self, port: usize) -> EndOfStep {
        self.element.end_of_step(port)
    }

    fn propose(
        &self,
        old: &C::State,
        dt_s: f64,
        inputs: &[f64],
        new: &mut C::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, C::Error> {
        self.element.propose(old, dt_s, inputs, new, outputs, cx)
    }

    fn admit(&self, old: &C::State, new: &C::State, receipt: &Receipt) -> Result<(), C::Error> {
        self.element.admit(old, new, receipt)
    }

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(<C::State as Words>::CHARTER);

    fn encode_state(&self, state: &C::State, words: &mut Vec<u64>) -> bool {
        state.put(words);
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<C::State> {
        <C::State as Words>::get(words)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend(&self.identity);
    }
}

impl<C> PortTangent for Coded<C>
where
    C: PortTangent,
    C::State: Words,
{
    fn jvp(
        &self,
        old: &C::State,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), C::Error> {
        self.element.jvp(old, dt_s, inputs, d_inputs, d_outputs, cx)
    }

    fn vjp(
        &self,
        old: &C::State,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), C::Error> {
        self.element
            .vjp(old, dt_s, inputs, bar_outputs, bar_inputs, cx)
    }
}

#[derive(Debug, PartialEq, Eq)]
struct Refused;

/// A member without ports that records the bits of every accepted step
/// length: a state whose number of words grows with the run.
struct Recorder;

impl Component for Recorder {
    type State = Vec<u64>;
    type Error = Refused;

    fn ports(&self) -> &[PortSpec] {
        &[]
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, _state: &Vec<u64>) -> Result<f64, Refused> {
        Ok(0.0)
    }

    fn propose(
        &self,
        old: &Vec<u64>,
        dt_s: f64,
        _inputs: &[f64],
        new: &mut Vec<u64>,
        _outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, Refused> {
        new.clone_from(old);
        new.push(dt_s.to_bits());
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::tests::snapshots",
        state_family: "recorder",
        schema_grammar: "the bits of each accepted step length, in order",
        codec_grammar: "one u64 word per step",
        codec_version: 1,
    });

    fn encode_state(&self, state: &Vec<u64>, words: &mut Vec<u64>) -> bool {
        words.extend(state);
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<Vec<u64>> {
        Some(words.to_vec())
    }
}

/// The charter of a [`Window`].
trait Tag {
    const CHARTER: StateIdentityCharterV2;
}

enum TagA {}
enum TagB {}

impl Tag for TagA {
    const CHARTER: StateIdentityCharterV2 = StateIdentityCharterV2 {
        owner: "fs-system::tests::snapshots",
        state_family: "window-a",
        schema_grammar: "width u64 values",
        codec_grammar: "width u64 words",
        codec_version: 1,
    };
}

impl Tag for TagB {
    const CHARTER: StateIdentityCharterV2 = StateIdentityCharterV2 {
        owner: "fs-system::tests::snapshots",
        state_family: "window-b",
        schema_grammar: "width u64 values",
        codec_grammar: "width u64 words",
        codec_version: 1,
    };
}

/// A member holding `width` words that leaves its width out of its
/// identity on purpose, so that a system of another width is of the same
/// problem and only the member's own decode can refuse.
struct Window<T> {
    width: usize,
    tag: PhantomData<T>,
}

impl<T> Window<T> {
    fn new(width: usize) -> Self {
        Self {
            width,
            tag: PhantomData,
        }
    }
}

impl<T: Tag> Component for Window<T> {
    type State = Vec<u64>;
    type Error = Refused;

    fn ports(&self) -> &[PortSpec] {
        &[]
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, _state: &Vec<u64>) -> Result<f64, Refused> {
        Ok(0.0)
    }

    fn propose(
        &self,
        old: &Vec<u64>,
        _dt_s: f64,
        _inputs: &[f64],
        new: &mut Vec<u64>,
        _outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, Refused> {
        new.clone_from(old);
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(T::CHARTER);

    fn encode_state(&self, state: &Vec<u64>, words: &mut Vec<u64>) -> bool {
        words.extend(state);
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<Vec<u64>> {
        (words.len() == self.width).then(|| words.to_vec())
    }
}

// ---------------------------------------------------------------------------
// The two graphs
// ---------------------------------------------------------------------------

/// A graph that runs, shows everything it holds, and is rebuilt the same in
/// any process.
trait Graph {
    /// The step after which it is sealed, and the step it runs to.
    const SPLIT: usize;
    const END: usize;
    fn build() -> Self;
    fn system(&mut self) -> &mut System;
    fn dt(index: usize) -> f64;
    /// Everything the system shows, exactly: float bits and the ledger's
    /// `Debug`, whose sums print their running value and compensation.
    fn shown(&self) -> String;
}

/// Which members the driven, dissipative device is built with.
#[derive(Clone, Copy)]
struct Variant {
    stiffness: f64,
    window: usize,
    tag_b: bool,
    coded_body: bool,
}

impl Default for Variant {
    fn default() -> Self {
        Self {
            stiffness: 60.0,
            window: 2,
            tag_b: false,
            coded_body: true,
        }
    }
}

/// The device of tests/checkpoints.rs: a mass on a spring whose far end is
/// a moving support, driven by a sine force whose account holds a clock,
/// and a damper whose heat a route gives 0.8 to a body and the rest to the
/// ambient sink; with a recorder and a window.
struct Device {
    system: System,
    mass: NodeId,
    spring: NodeId,
    source: NodeId,
    body: NodeId,
    recorder: NodeId,
    window: NodeId,
    ports: Vec<PortRef>,
}

impl Device {
    fn with(variant: Variant) -> Self {
        let mut builder = SystemBuilder::new();
        let mass = builder
            .add_with_tangent(
                "mass",
                Coded::new(Mass::new(0.8).unwrap(), &[0.8]),
                Motion {
                    position_m: 0.07,
                    velocity_m_s: 0.34,
                },
            )
            .unwrap();
        let spring = builder
            .add_with_tangent(
                "spring",
                Coded::new(
                    Spring::new(variant.stiffness, 0.0).unwrap(),
                    &[variant.stiffness, 0.0],
                ),
                -0.07,
            )
            .unwrap();
        let damper = builder
            .add_with_tangent("damper", Coded::new(Damper::new(1.1).unwrap(), &[1.1]), ())
            .unwrap();
        let source = builder
            .add_with_tangent(
                "source",
                Coded::new(
                    Force::new(Waveform::sine(1.2, 1.7, 0.3, 0.2).unwrap()),
                    &[1.2, 1.7, 0.3, 0.2],
                ),
                SourceAccount::default(),
            )
            .unwrap();
        let body = if variant.coded_body {
            builder.add_with_tangent(
                "body",
                Coded::new(HeatCapacity::new(4.0).unwrap(), &[4.0]),
                300.0,
            )
        } else {
            builder.add_with_tangent("body", HeatCapacity::new(4.0).unwrap(), 300.0)
        }
        .unwrap();
        let recorder = builder.add("recorder", Recorder, Vec::new()).unwrap();
        let held = vec![7; variant.window];
        let window = if variant.tag_b {
            builder.add("window", Window::<TagB>::new(variant.window), held)
        } else {
            builder.add("window", Window::<TagA>::new(variant.window), held)
        }
        .unwrap();
        let ports = [
            (mass, "p"),
            (spring, "a"),
            (spring, "b"),
            (damper, "a"),
            (damper, "b"),
            (damper, "heat"),
            (body, "heat"),
            (source, "p"),
        ]
        .map(|(node, name)| builder.port(node, name).unwrap());
        builder
            .common_flow(ports[0], &[ports[1], ports[3], ports[7]])
            .unwrap();
        builder
            .boundary(
                ports[2],
                Box::new(|t0, dt, value| value[0] = 0.05 + 0.02 * (t0 + 0.5 * dt)),
            )
            .unwrap();
        builder
            .boundary(ports[4], Box::new(|_, _, value| value[0] = 0.0))
            .unwrap();
        builder
            .heat_route(
                ports[5],
                &[(ports[6], 0.8)],
                Some(heat_sink("ambient", 300.0).unwrap()),
                RouteTemperature::Fixed(300.0),
            )
            .unwrap();
        Self {
            system: builder.seal().unwrap(),
            mass,
            spring,
            source,
            body,
            recorder,
            window,
            ports: ports.to_vec(),
        }
    }
}

/// The bits of what every port in `ends` last received and returned.
fn wires(system: &System, ends: &[PortRef]) -> Vec<(Vec<u64>, Vec<u64>)> {
    let bits = |values: &[f64]| values.iter().map(|value| value.to_bits()).collect();
    ends.iter()
        .map(|&port| {
            (
                bits(system.last_input(port)),
                bits(system.last_output(port)),
            )
        })
        .collect()
}

impl Graph for Device {
    const SPLIT: usize = 73;
    const END: usize = 257;

    fn build() -> Self {
        Self::with(Variant::default())
    }

    fn system(&mut self) -> &mut System {
        &mut self.system
    }

    fn dt(index: usize) -> f64 {
        [0.003, 0.0017, 0.0023, 0.0009][index % 4]
    }

    fn shown(&self) -> String {
        let system = &self.system;
        let motion = system.state::<Motion>(self.mass).unwrap();
        let mut out = String::new();
        let _ = write!(
            out,
            "t {:016x} n {} x {:016x} v {:016x} spring {:016x} source {:?} body {:016x} \
             recorder {:?} window {:?} wires {:?} ledger {:?}",
            system.time_s().to_bits(),
            system.steps(),
            motion.position_m.to_bits(),
            motion.velocity_m_s.to_bits(),
            system.state::<f64>(self.spring).unwrap().to_bits(),
            system.state::<SourceAccount>(self.source).unwrap(),
            system.state::<f64>(self.body).unwrap().to_bits(),
            system.state::<Vec<u64>>(self.recorder).unwrap(),
            system.state::<Vec<u64>>(self.window).unwrap(),
            wires(system, &self.ports),
            system.ledger(),
        );
        out
    }
}

const AIR_J_KG_K: f64 = 287.0;
const RATIO: f64 = 1.4;

fn alternate_heat(t: f64) -> [f64; 2] {
    if t < 1.0 {
        [800.0, 0.0]
    } else if t < 2.0 {
        [0.0, 4_000.0]
    } else {
        [0.0, 0.0]
    }
}

/// The flow that turns round of tests/elements.rs: two vessels of air
/// carrying two quantities, joined by a restriction that reports the
/// entropy of its pressure drop; `a` is heated first and `b` after.
struct Vessels {
    system: System,
    nodes: [NodeId; 2],
    ends: Vec<PortRef>,
}

impl Vessels {
    fn with(gas_constant: f64, conductance: f64) -> Self {
        let air = || GasVolume::new(gas_constant, RATIO, 2).unwrap();
        let start = [
            air().contents(2.0e5, 300.0, 0.01, &[0.3, 0.0]).unwrap(),
            air().contents(2.0e5, 300.0, 0.03, &[0.0, 0.6]).unwrap(),
        ];
        let mut builder = SystemBuilder::new();
        let pipe = builder
            .add_with_tangent(
                "pipe",
                Restriction::for_gas(conductance, &air()).unwrap(),
                (),
            )
            .unwrap();
        let mut nodes = Vec::new();
        let mut ends = vec![
            builder.port(pipe, "a").unwrap(),
            builder.port(pipe, "b").unwrap(),
        ];
        for (side, (name, contents)) in ["a", "b"].into_iter().zip(start).enumerate() {
            let node = builder.add_with_tangent(name, air(), contents).unwrap();
            builder
                .bond(
                    builder.port(node, "stream").unwrap(),
                    builder.port(pipe, name).unwrap(),
                )
                .unwrap();
            let port = builder.port(node, "heat").unwrap();
            builder
                .boundary(
                    port,
                    Box::new(move |t0, dt, values| {
                        values[0] = alternate_heat(t0 + 0.5 * dt)[side];
                    }),
                )
                .unwrap();
            ends.extend([builder.port(node, "stream").unwrap(), port]);
            nodes.push(node);
        }
        Self {
            system: builder.seal().unwrap(),
            nodes: [nodes[0], nodes[1]],
            ends,
        }
    }
}

impl Graph for Vessels {
    /// Sealed at t = 1.25 s, while `b` is heated and the flow turns.
    const SPLIT: usize = 250;
    const END: usize = 600;

    fn build() -> Self {
        Self::with(AIR_J_KG_K, 2.0e-8)
    }

    fn system(&mut self) -> &mut System {
        &mut self.system
    }

    fn dt(_index: usize) -> f64 {
        3.0 / 600.0
    }

    fn shown(&self) -> String {
        let system = &self.system;
        let held: Vec<Vec<u64>> = self
            .nodes
            .iter()
            .map(|&node| {
                let contents = system.state::<GasContents>(node).unwrap();
                [contents.mass_kg, contents.energy_j, contents.volume_m3]
                    .iter()
                    .chain(&contents.carried_kg)
                    .map(|value| value.to_bits())
                    .collect()
            })
            .collect();
        let mut out = String::new();
        let _ = write!(
            out,
            "t {:016x} n {} held {:?} wires {:?} ledger {:?}",
            system.time_s().to_bits(),
            system.steps(),
            held,
            wires(system, &self.ends),
            system.ledger(),
        );
        out
    }
}

/// What the graph shows after each step from `start` to its end.
fn run<G: Graph>(graph: &mut G, start: usize, end: usize) -> Vec<String> {
    with_cx(|cx| {
        (start..end)
            .map(|index| {
                graph
                    .system()
                    .step(G::dt(index), cx)
                    .unwrap_or_else(|error| panic!("step {index}: {error:?}"));
                graph.shown()
            })
            .collect()
    })
}

/// What an uninterrupted run shows after its split and after every later
/// step, in that order.
fn uninterrupted<G: Graph>() -> Vec<String> {
    let mut graph = G::build();
    let mut shown = run(&mut graph, 0, G::SPLIT).split_off(G::SPLIT - 1);
    shown.extend(run(&mut graph, G::SPLIT, G::END));
    shown
}

// ---------------------------------------------------------------------------
// Kill and resume in another process
// ---------------------------------------------------------------------------

/// The authorized route's verifier: it accepts only the subject the parent
/// kept and passed on.
struct Kept([u8; 32]);

enum KeptSchema {}
impl CanonicalSchema for KeptSchema {
    const DOMAIN: &'static str = "org.frankensim.fs-system.tests.snapshots.verifier.v1";
    const NAME: &'static str = "kept-subject-verifier";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "the subject id the sealing process kept";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("version", WireType::U64)];
}

enum OpenPolicySchema {}
impl CanonicalSchema for OpenPolicySchema {
    const DOMAIN: &'static str = "org.frankensim.fs-system.tests.snapshots.policy.v1";
    const NAME: &'static str = "resume-in-a-child-policy";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "a child process resumes what its parent sealed";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("version", WireType::U64)];
}

type Subject<S> = AuthorityRef<SnapshotAuthoritySubjectIdV2, KeptSchema, OpenPolicySchema, S>;

impl AuthorityVerifier<SnapshotAuthoritySubjectIdV2, KeptSchema, OpenPolicySchema> for Kept {
    type Error = &'static str;

    fn verify(&self, presented: &Subject<Presented>) -> Result<(), &'static str> {
        if presented.receipt().id().as_bytes() == &self.0 {
            Ok(())
        } else {
            Err("the subject is not the one the sealing process kept")
        }
    }
}

struct OpenPolicy;

impl AuthorityAdmitter<SnapshotAuthoritySubjectIdV2, KeptSchema, OpenPolicySchema> for OpenPolicy {
    type Error = &'static str;

    fn admit(&self, _verified: &Subject<Verified>) -> Result<(), &'static str> {
        Ok(())
    }
}

/// Present the subject of the bytes received, verified against `kept`.
fn admitted(
    receipt: IdentityReceipt<SnapshotAuthoritySubjectIdV2>,
    kept: [u8; 32],
) -> Result<Subject<Admitted>, &'static str> {
    let verifier = VerifierId::<KeptSchema>::parse_slice(&[0x22; 32]).unwrap();
    let policy = KeyPolicyId::<OpenPolicySchema>::parse_slice(&[0x33; 32]).unwrap();
    let anchor = ExternalAnchorRef::presented(ContentId::of_bytes(&kept));
    Subject::<Presented>::present(receipt, anchor, verifier, policy)
        .verify(&Kept(kept))?
        .admit(&OpenPolicy)
}

/// The child: build the graph, refuse a wrong subject and a flipped byte
/// with nothing changed, resume, and print what it shows at every step.
fn resume_here<G: Graph>(kept: [u8; 32], bytes: &[u8]) {
    let mut graph = G::build();
    let fresh = graph.shown();
    let context = context(graph.system());
    let receipt = snapshot_v2::inspect(bytes, limits(), || false)
        .expect("the envelope is whole")
        .authority_subject_receipt();
    // Everything printed goes out at the end in one write, from a line of
    // its own: libtest has printed the test's name without a newline.
    let mut out = String::from("\n");
    let mut wrong = kept;
    wrong[0] ^= 1;
    let refused = admitted(receipt, wrong).expect_err("a subject the parent did not keep");
    assert_eq!(graph.shown(), fresh);
    let _ = writeln!(out, "REFUSED wrong subject: {refused}");
    let authority = admitted(receipt, kept).expect("the subject the parent kept");
    let mut flipped = bytes.to_vec();
    let last = flipped.len() - 1;
    flipped[last] ^= 1;
    let refused =
        SystemImage::unseal_v2_authorized(&flipped, &authority, &context, limits(), || false)
            .expect_err("a flipped byte is refused");
    assert_eq!(graph.shown(), fresh);
    let _ = writeln!(out, "REFUSED flipped byte: {refused:?}");
    let opened = SystemImage::unseal_v2_authorized(bytes, &authority, &context, limits(), || false)
        .expect("opened by the authorized route");
    graph.system().resume(&opened).expect("resumed");
    let _ = writeln!(out, "STEP {} {}", G::SPLIT, graph.shown());
    for (offset, shown) in run(&mut graph, G::SPLIT, G::END).into_iter().enumerate() {
        let _ = writeln!(out, "STEP {} {shown}", G::SPLIT + 1 + offset);
    }
    let mut stdout = std::io::stdout().lock();
    stdout.write_all(out.as_bytes()).unwrap();
    stdout.flush().unwrap();
}

#[test]
#[ignore = "the child process of the kill-and-resume tests, which run it with --ignored"]
fn child_resumes_from_stdin() {
    let graph = std::env::var(CHILD).expect("run only as the child of a kill-and-resume test");
    let mut input = Vec::new();
    std::io::stdin().read_to_end(&mut input).unwrap();
    let (kept, bytes) = input.split_at(32);
    let kept: [u8; 32] = kept.try_into().unwrap();
    match graph.as_str() {
        "device" => resume_here::<Device>(kept, bytes),
        "vessels" => resume_here::<Vessels>(kept, bytes),
        other => panic!("no graph {other}"),
    }
}

/// Seal at the split, drop the sealing system, and resume in a child that
/// receives only the envelope bytes and the subject id on its stdin.
fn across_processes<G: Graph>(name: &str) {
    let want = uninterrupted::<G>();
    let mut first = G::build();
    run(&mut first, 0, G::SPLIT);
    let system = first.system();
    let sealed = system
        .snapshot(&context(system), limits(), || false)
        .expect("sealed");
    drop(first);
    let kept = *sealed.authority_subject_receipt().id().as_bytes();
    let mut child = Command::new(std::env::current_exe().unwrap())
        .args([
            "child_resumes_from_stdin",
            "--exact",
            "--ignored",
            "--nocapture",
            "--test-threads",
            "1",
        ])
        .env(CHILD, name)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    {
        let mut stdin = child.stdin.take().unwrap();
        stdin.write_all(&kept).unwrap();
        stdin.write_all(sealed.bytes()).unwrap();
    }
    let output = child.wait_with_output().unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(
        output.status.success(),
        "child failed: {}\n{}",
        stdout,
        String::from_utf8_lossy(&output.stderr)
    );
    let refused: Vec<&str> = stdout
        .lines()
        .filter(|line| line.starts_with("REFUSED "))
        .collect();
    assert_eq!(refused.len(), 2, "{stdout}");
    let steps: Vec<&str> = stdout
        .lines()
        .filter_map(|line| line.strip_prefix("STEP "))
        .collect();
    assert_eq!(steps.len(), want.len(), "steps the child printed");
    for (offset, (got, shown)) in steps.iter().zip(&want).enumerate() {
        let index = G::SPLIT + offset;
        assert_eq!(*got, format!("{index} {shown}"), "step {index}");
    }
    eprintln!(
        "{name}: {} envelope bytes; the child refused [{}]; then {} steps from {} to {} the \
         same to the bit",
        sealed.bytes().len(),
        refused.join("; "),
        want.len() - 1,
        G::SPLIT,
        G::END
    );
}

#[test]
fn a_device_killed_and_resumed_in_another_process_follows_the_uninterrupted_run() {
    across_processes::<Device>("device");
}

#[test]
fn a_turning_gas_flow_killed_and_resumed_in_another_process_follows_the_uninterrupted_run() {
    across_processes::<Vessels>("vessels");
}

// ---------------------------------------------------------------------------
// Fork, rewind and refusals, in this process
// ---------------------------------------------------------------------------

#[test]
fn one_snapshot_forks_into_systems_that_go_on_alone_and_rewinds_its_writer() {
    let want = uninterrupted::<Device>();
    let mut writer = Device::build();
    run(&mut writer, 0, Device::SPLIT);
    let sealed = writer
        .system
        .snapshot(&context(&writer.system), limits(), || false)
        .expect("sealed");
    let opened = opened_here(&sealed);
    let (mut left, mut right) = (Device::build(), Device::build());
    left.system().resume(&opened).unwrap();
    right.system().resume(&opened).unwrap();
    assert_eq!(left.shown(), want[0]);
    assert_eq!(right.shown(), want[0]);
    // Stepping one leaves the other to the bit.
    let ahead = run(&mut left, Device::SPLIT, Device::SPLIT + 20);
    assert_eq!(right.shown(), want[0]);
    assert_eq!(ahead[..], want[1..=20]);
    // Both go on to the uninterrupted trajectory.
    assert_eq!(run(&mut right, Device::SPLIT, Device::END)[..], want[1..]);
    assert_eq!(
        run(&mut left, Device::SPLIT + 20, Device::END)[..],
        want[21..]
    );
    // The writer, run on to the end, is rewound by its own snapshot and
    // replays the same steps.
    assert_eq!(run(&mut writer, Device::SPLIT, Device::END)[..], want[1..]);
    writer.system().resume(&opened).unwrap();
    assert_eq!(writer.shown(), want[0]);
    assert_eq!(run(&mut writer, Device::SPLIT, Device::END)[..], want[1..]);
}

/// Resume `opened` into `graph` after it ran five steps: the refusal must
/// be `expected`, and the graph must show what it showed before.
fn refused<G: Graph>(
    mut graph: G,
    opened: &snapshot_v2::OpenedSnapshotV2<SystemImage>,
    expected: &str,
) {
    run(&mut graph, 0, 5);
    let before = graph.shown();
    let error = graph.system().resume(opened).expect_err(expected);
    assert_eq!(format!("{error:?}"), expected);
    assert_eq!(graph.shown(), before, "{expected}");
}

#[test]
fn every_refusal_leaves_the_system_as_it_was() {
    let mut writer = Device::build();
    run(&mut writer, 0, Device::SPLIT);
    let sealed = writer
        .system
        .snapshot(&context(&writer.system), limits(), || false)
        .expect("sealed");
    let opened = opened_here(&sealed);

    // Another graph.
    refused(Vessels::build(), &opened, "ForeignProblem");
    // The same graph with the spring stiffer by one ulp.
    let stiffer = Variant {
        stiffness: 60.0_f64.next_up(),
        ..Variant::default()
    };
    assert_ne!(
        Device::with(stiffer).system.problem_id(),
        writer.system.problem_id()
    );
    refused(Device::with(stiffer), &opened, "ForeignProblem");
    // A member whose codec has another charter.
    let other_tag = Variant {
        tag_b: true,
        ..Variant::default()
    };
    refused(Device::with(other_tag), &opened, "ForeignProblem");
    // A member without a codec cannot be resumed into or sealed.
    let uncoded = Variant {
        coded_body: false,
        ..Variant::default()
    };
    refused(Device::with(uncoded), &opened, "ForeignProblem");
    let plain = Device::with(uncoded);
    let error = plain
        .system
        .snapshot(&context(&plain.system), limits(), || false)
        .expect_err("a member without a codec");
    assert_eq!(format!("{error:?}"), r#"NoCodec { node: "body" }"#);
    // The last member refuses words of another width; the members before
    // it decoded, and none of them is published.
    let wider = Variant {
        window: 3,
        ..Variant::default()
    };
    assert_eq!(
        Device::with(wider).system.problem_id(),
        writer.system.problem_id()
    );
    refused(Device::with(wider), &opened, r#"Decode { node: "window" }"#);
    // A context of another problem seals nothing.
    let other = Vessels::build();
    let error = writer
        .system
        .snapshot(&context(&other.system), limits(), || false)
        .expect_err("another problem's context");
    assert_eq!(format!("{error:?}"), "ForeignProblem");

    // Truncated or flipped bytes never open.
    let bytes = sealed.bytes();
    let expectation = sealed.expectation();
    let open = |bytes: &[u8]| {
        SystemImage::unseal_v2_expected(bytes, &expectation, limits(), || false).err()
    };
    assert!(open(&bytes[..bytes.len() - 1]).is_some(), "truncated");
    for at in [0, bytes.len() / 2, bytes.len() - 1] {
        let mut flipped = bytes.to_vec();
        flipped[at] ^= 0x10;
        assert!(open(&flipped).is_some(), "byte {at} flipped");
    }
    // The writer itself is where it was.
    let mut again = Device::build();
    run(&mut again, 0, Device::SPLIT);
    assert_eq!(writer.shown(), again.shown());
}

#[test]
fn the_gas_elements_are_of_their_parameters_to_the_ulp() {
    let opened = {
        let mut writer = Vessels::build();
        run(&mut writer, 0, 40);
        let sealed = writer
            .system
            .snapshot(&context(&writer.system), limits(), || false)
            .expect("sealed");
        opened_here(&sealed)
    };
    refused(
        Vessels::with(AIR_J_KG_K, 2.0e-8_f64.next_up()),
        &opened,
        "ForeignProblem",
    );
    refused(
        Vessels::with(AIR_J_KG_K.next_down(), 2.0e-8),
        &opened,
        "ForeignProblem",
    );
    refused(Device::build(), &opened, "ForeignProblem");
    // The contents' codec: every bit back, and words no volume holds are
    // refused.
    let gas = GasVolume::new(AIR_J_KG_K, RATIO, 2).unwrap();
    let contents = gas.contents(2.0e5, 300.0, 0.01, &[0.3, 0.0]).unwrap();
    let mut words = Vec::new();
    assert!(gas.encode_state(&contents, &mut words));
    assert_eq!(words.len(), 5);
    let back = gas.decode_state(&words).unwrap();
    assert_eq!(format!("{back:?}"), format!("{contents:?}"));
    assert!(gas.decode_state(&words[..4]).is_none(), "too few words");
    let mut empty = words.clone();
    empty[0] = 0.0_f64.to_bits();
    assert!(gas.decode_state(&empty).is_none(), "no mass");
    let restriction = Restriction::for_gas(2.0e-8, &gas).unwrap();
    let mut none = Vec::new();
    assert!(restriction.encode_state(&(), &mut none));
    assert!(none.is_empty());
    assert!(
        restriction.decode_state(&[1]).is_none(),
        "a word for no state"
    );
}

#[test]
fn the_image_charter_derives_its_pinned_identities() {
    let charter = SystemImage::charter().expect("a chartered family");
    let hex = |bytes: &[u8; 32]| {
        bytes.iter().fold(String::new(), |mut out, byte| {
            let _ = write!(out, "{byte:02x}");
            out
        })
    };
    assert_eq!(
        SystemImage::STATE_TYPE_ID_V2,
        charter.state_type_id(),
        "derived {}",
        hex(charter.state_type_id().as_bytes())
    );
    assert_eq!(
        SystemImage::STATE_SCHEMA_ID_V2,
        charter.state_schema_id(),
        "derived {}",
        hex(charter.state_schema_id().as_bytes())
    );
    assert_eq!(
        SystemImage::STATE_CODEC_ID_V2,
        charter.state_codec_id(),
        "derived {}",
        hex(charter.state_codec_id().as_bytes())
    );
    assert!(snapshot_v2::verify_state_charter::<SystemImage>().is_ok());
}
