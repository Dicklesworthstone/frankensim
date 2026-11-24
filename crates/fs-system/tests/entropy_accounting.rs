//! G0: the kernel audits the stored entropy of every component that
//! declares it, at every step and before anything is published. The ledger
//! keeps one entropy account per declaring component; a component that
//! declares nothing has none. The accounts rewind with a checkpoint and
//! come back from a snapshot. The shipped `elements::Conductance` declares.

#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::identity::CanonicalLimits;
use fs_exec::solver::SolverStateV2;
use fs_exec::solver::snapshot_v2::{
    ExpectedResumeContextV2, PausedSnapshotBoundaryV2, SnapshotAlgorithmIdV2,
    SnapshotBudgetStateIdV2, SnapshotDeterminismV2, SnapshotExecutionFingerprintIdV2,
    SnapshotLimitsV2, SnapshotPauseRequestIdV2, SnapshotProvenanceIdV2, SnapshotRngCounterIdV2,
    StateIdentityCharterV2,
};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, DrainTracker, ExecMode, RunId, StreamKey};
use fs_ports::StableId;
use fs_ports::component::{
    Class, Component, ContractViolation, PortRole, PortSpec, Proposal, Receipt,
};
use fs_system::elements::Conductance;
use fs_system::ledger::{EntropyAccount, Ledger};
use fs_system::{NodeId, Phase, System, SystemBuilder, SystemError, SystemImage};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 41,
            kernel_id: 5,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

fn heat(name: &str, role: PortRole) -> PortSpec {
    PortSpec::heat(StableId::new(name).expect("id"), role).expect("heat port")
}

/// A body of heat capacity `capacity_j_k` at one temperature [K], its
/// state. Its heat port takes the heat rate `q` in and returns the end
/// temperature of a backward-Euler step, `T1 = T0 + dt q / C`. Declared, its
/// stored entropy is `C ln(T / reference_k)` times `skew`, and it books
/// `C ln(T1 / T0) - dt q / T1` as generated, which closes the account to
/// rounding when `skew` is 1.
///
/// Both logarithms are taken as `ln_1p` of a difference that is exact in
/// this range, so each rounds relative to the entropy it gives. The
/// kernel's tolerance is relative to the reported inventories and needs
/// that. `textbook` takes `ln` of the quotient instead, which near the
/// reference rounds at `C EPSILON`, far above the inventory's own size.
/// `bias_j_per_k` is added to every generation it books, and a body that
/// `refuses` answers every measurement of its entropy with `Cancelled`.
struct Body {
    ports: [PortSpec; 1],
    capacity_j_k: f64,
    reference_k: f64,
    skew: f64,
    bias_j_per_k: f64,
    declared: bool,
    textbook: bool,
    refuses: bool,
}

impl Body {
    fn new(capacity_j_k: f64) -> Self {
        Self {
            ports: [heat("heat", PortRole::FlowIn)],
            capacity_j_k,
            reference_k: 300.0,
            skew: 1.0,
            bias_j_per_k: 0.0,
            declared: true,
            textbook: false,
            refuses: false,
        }
    }

    /// `ln(to / from)`.
    fn log_ratio(&self, to: f64, from: f64) -> f64 {
        if self.textbook {
            (to / from).ln()
        } else {
            ((to - from) / from).ln_1p()
        }
    }
}

impl Component for Body {
    type State = f64;
    type Error = Cancelled;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::tests::entropy_accounting::Body",
        state_family: "temperature",
        schema_grammar: "temperature_k: f64",
        codec_grammar: "one u64 word: the bits of temperature_k",
        codec_version: 1,
    });

    fn encode_state(&self, state: &f64, words: &mut Vec<u64>) -> bool {
        words.push(state.to_bits());
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<f64> {
        let &[word] = words else { return None };
        Some(f64::from_bits(word))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([
            self.capacity_j_k.to_bits(),
            self.reference_k.to_bits(),
            self.skew.to_bits(),
            self.bias_j_per_k.to_bits(),
            u64::from(self.declared),
            u64::from(self.textbook),
            u64::from(self.refuses),
        ]);
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, kelvin: &f64) -> Result<f64, Cancelled> {
        Ok(self.capacity_j_k * kelvin)
    }

    fn stored_entropy_j_per_k(&self, kelvin: &f64) -> Result<Option<f64>, Cancelled> {
        if self.refuses {
            return Err(Cancelled);
        }
        Ok(self
            .declared
            .then(|| self.skew * self.capacity_j_k * self.log_ratio(*kelvin, self.reference_k)))
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        let heat_j = dt_s * inputs[0];
        *new = old + heat_j / self.capacity_j_k;
        outputs[0] = *new;
        Ok(Proposal::Step(Receipt {
            entropy_generated_j_per_k: self.capacity_j_k * self.log_ratio(*new, *old)
                - heat_j / *new
                + self.bias_j_per_k,
            ..Receipt::lossless(heat_j, 0.0)
        }))
    }
}

/// A wall of conductance `conductance_w_k` between heat ports `a` and `b`,
/// with no state: it takes `G (Ta - Tb)` in at `a` and gives it out at `b`,
/// and books `dt G (Ta - Tb)^2 / (Ta Tb)` as generated. Declared, it stores
/// no entropy, `Some(0.0)`.
struct Wall {
    ports: [PortSpec; 2],
    conductance_w_k: f64,
    declared: bool,
}

impl Wall {
    fn new(conductance_w_k: f64) -> Self {
        Self {
            ports: [heat("a", PortRole::EffortIn), heat("b", PortRole::EffortIn)],
            conductance_w_k,
            declared: true,
        }
    }
}

impl Component for Wall {
    type State = ();
    type Error = Cancelled;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::tests::entropy_accounting::Wall",
        state_family: "stateless",
        schema_grammar: "()",
        codec_grammar: "no words",
        codec_version: 1,
    });

    fn encode_state(&self, (): &(), _words: &mut Vec<u64>) -> bool {
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<()> {
        words.is_empty().then_some(())
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([self.conductance_w_k.to_bits(), u64::from(self.declared)]);
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, (): &()) -> Result<f64, Cancelled> {
        Ok(0.0)
    }

    fn stored_entropy_j_per_k(&self, (): &()) -> Result<Option<f64>, Cancelled> {
        Ok(self.declared.then_some(0.0))
    }

    fn propose(
        &self,
        (): &(),
        dt_s: f64,
        inputs: &[f64],
        _new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        let (hot, cold) = (inputs[0], inputs[1]);
        let rate_w = self.conductance_w_k * (hot - cold);
        outputs[0] = rate_w;
        outputs[1] = -rate_w;
        Ok(Proposal::Step(Receipt {
            entropy_generated_j_per_k: dt_s * self.conductance_w_k * (hot - cold) * (hot - cold)
                / (hot * cold),
            ..Receipt::lossless(0.0, 0.0)
        }))
    }
}

const HOT: (f64, f64) = (50.0, 400.0);
const COLD: (f64, f64) = (120.0, 300.0);
const DT_S: f64 = 0.05;

/// A hot and a cold body joined through a wall, all three declaring unless
/// `tune` changes them.
struct Pair {
    system: System,
    hot: NodeId,
    cold: NodeId,
    wall: NodeId,
}

fn pair(tune: impl FnOnce(&mut Body, &mut Body, &mut Wall)) -> Pair {
    let (mut hot_body, mut cold_body) = (Body::new(HOT.0), Body::new(COLD.0));
    let mut wall_member = Wall::new(3.0);
    tune(&mut hot_body, &mut cold_body, &mut wall_member);
    joined((hot_body, HOT.1), (cold_body, COLD.1), wall_member)
}

/// Two bodies, each at its starting temperature, joined through `wall`, a
/// stateless member with heat ports `a` and `b`.
fn joined<W: Component<State = ()> + 'static>(
    (hot_body, hot_k): (Body, f64),
    (cold_body, cold_k): (Body, f64),
    wall_member: W,
) -> Pair {
    let mut builder = SystemBuilder::new();
    let hot = builder.add("hot", hot_body, hot_k).expect("hot");
    let cold = builder.add("cold", cold_body, cold_k).expect("cold");
    let wall = builder.add("wall", wall_member, ()).expect("wall");
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).expect("port");
    let (a, b) = (port(&builder, wall, "a"), port(&builder, wall, "b"));
    let (h, c) = (port(&builder, hot, "heat"), port(&builder, cold, "heat"));
    builder.common_effort(h, &[a]).expect("hot side");
    builder.common_effort(c, &[b]).expect("cold side");
    let mut system = builder.seal().expect("sealed");
    assert!(system.guess(a, &[hot_k]));
    assert!(system.guess(b, &[cold_k]));
    Pair {
        system,
        hot,
        cold,
        wall,
    }
}

/// `steps` steps of `DT_S`, none of them shortened.
fn run(system: &mut System, steps: u32) {
    with_cx(|cx| {
        for index in 0..steps {
            let report = system.step(DT_S, cx).expect("step");
            assert_eq!(report.dt_s, DT_S, "step {index} was shortened");
            assert_eq!(report.shortenings, 0, "step {index}");
        }
    });
}

fn entropy<'a>(ledger: &'a Ledger, name: &str) -> Option<&'a EntropyAccount> {
    ledger
        .nodes
        .iter()
        .find(|node| node.name == name)
        .expect("node")
        .entropy
        .as_ref()
}

#[test]
fn every_declared_entropy_is_audited_and_the_pair_closes_its_entropy() {
    let Pair {
        mut system,
        hot,
        cold,
        ..
    } = pair(|_, _, _| {});
    run(&mut system, 200);
    assert_eq!(system.steps(), 200);
    let ledger = system.ledger();
    let mut generated = 0.0;
    let mut stored = 0.0;
    for name in ["hot", "cold", "wall"] {
        let account = entropy(ledger, name).expect("declared");
        assert_eq!(account.audited_steps, 200, "{name}");
        let defect = account.absolute_defect_j_per_k.value();
        assert!(defect <= 1.0e-10, "{name}: {defect:e} J/K");
        let port = account.port_entropy_in_j_per_k.value();
        let made = account.entropy_generated_j_per_k.value();
        let change = account.stored_change_j_per_k.value();
        let signed = account.balance_defect_j_per_k.value();
        let scale = port.abs() + made.abs() + change.abs();
        assert!(
            (port + made - change - signed).abs() <= 1.0e-12 * scale,
            "{name}: the four accounts disagree"
        );
        generated += made;
        stored += change;
    }
    // The accounts agree with the states and with the ledger's own sums.
    for (name, node, (capacity, start)) in [("hot", hot, HOT), ("cold", cold, COLD)] {
        let end = *system.state::<f64>(node).expect("state");
        let measured = capacity * (end / start).ln();
        let booked = entropy(ledger, name)
            .expect("declared")
            .stored_change_j_per_k
            .value();
        assert!(
            (booked - measured).abs() <= 1.0e-12 * measured.abs(),
            "{name}"
        );
    }
    let total = ledger.entropy_generated_j_per_k();
    assert!(
        (generated - total).abs() <= 1.0e-12 * total,
        "{generated} vs {total}"
    );
    // Heat crosses each junction at one temperature, so the closed pair's
    // stored change is what its members generated.
    assert!(generated > 0.0);
    assert!(
        (stored - generated).abs() <= 1.0e-9 * generated,
        "{stored} vs {generated}"
    );
}

#[test]
fn the_shipped_conductance_is_audited_and_its_generation_cancels_its_port_entropy() {
    let conductance = || Conductance::new(3.0).expect("conductance");
    assert_eq!(
        conductance().stored_entropy_j_per_k(&()).expect("measured"),
        Some(0.0)
    );
    // At the pair's temperatures, and with the hot body ten times hotter
    // than the cold one.
    for hot_k in [HOT.1, 3000.0] {
        let hot = (Body::new(HOT.0), hot_k);
        let mut shipped = joined(hot, (Body::new(COLD.0), COLD.1), conductance());
        let hot = (Body::new(HOT.0), hot_k);
        let mut tested = joined(hot, (Body::new(COLD.0), COLD.1), Wall::new(3.0));
        run(&mut shipped.system, 200);
        run(&mut tested.system, 200);
        // The test wall has the same law: the states and every entropy
        // account agree to the bit.
        let state = |pair: &Pair, node| pair.system.state::<f64>(node).expect("state").to_bits();
        for (a, b) in [(shipped.hot, tested.hot), (shipped.cold, tested.cold)] {
            assert_eq!(state(&shipped, a), state(&tested, b), "{hot_k} K");
        }
        let (ledger, reference) = (shipped.system.ledger(), tested.system.ledger());
        for name in ["hot", "cold", "wall"] {
            assert_eq!(entropy(ledger, name), entropy(reference, name), "{name}");
        }
        // It stores nothing, and it generates what its heat carries out.
        let account = entropy(ledger, "wall").expect("declared");
        assert_eq!(account.audited_steps, 200);
        assert_eq!(account.stored_change_j_per_k.value(), 0.0);
        let port = account.port_entropy_in_j_per_k.value();
        let made = account.entropy_generated_j_per_k.value();
        assert!(made > 0.0, "{hot_k} K");
        assert!(
            (port + made).abs() <= 1.0e-12 * made,
            "{hot_k} K: {port:e} in, {made:e} generated"
        );
    }
}

#[test]
fn an_undeclared_pair_has_no_entropy_account_and_the_same_trajectory() {
    let mut declared = pair(|_, _, _| {});
    let mut undeclared = pair(|hot, cold, wall| {
        (hot.declared, cold.declared, wall.declared) = (false, false, false);
    });
    run(&mut declared.system, 50);
    run(&mut undeclared.system, 50);
    for name in ["hot", "cold", "wall"] {
        assert!(
            entropy(undeclared.system.ledger(), name).is_none(),
            "{name}"
        );
    }
    for node in [declared.hot, declared.cold] {
        let (a, b) = (
            declared.system.state::<f64>(node).expect("state"),
            undeclared.system.state::<f64>(node).expect("state"),
        );
        assert_eq!(a.to_bits(), b.to_bits(), "the audit moves no state");
    }
    assert_eq!(
        declared.system.ledger().entropy_generated_j_per_k(),
        undeclared.system.ledger().entropy_generated_j_per_k()
    );
}

#[test]
fn a_skewed_entropy_is_refused_by_name_and_publishes_nothing() {
    let Pair {
        mut system,
        hot,
        wall,
        ..
    } = pair(|hot, _, _| hot.skew = 1.0 + 1.0e-6);
    let before = system.ledger().clone();
    let refused = with_cx(|cx| system.step(DT_S, cx)).expect_err("refused");
    assert!(
        matches!(
            &refused,
            SystemError::Contract {
                node,
                violation: ContractViolation::EntropyBalance { .. },
            } if node == "hot"
        ),
        "{refused:?}"
    );
    assert_eq!(system.steps(), 0);
    assert_eq!(*system.state::<f64>(hot).expect("state"), HOT.1);
    assert!(system.state::<()>(wall).is_some());
    assert_eq!(*system.ledger(), before, "nothing is published");
}

#[test]
fn the_textbook_logarithm_is_refused_near_its_reference() {
    // The cold body starts at its 300 K reference. Its first step is exact
    // (the entropy before it is ln 1 = 0, and the generation's ln is the
    // stored one). At the second, ln(T / 300) of a quotient within 1e-3 of
    // one rounds at about 120 J/K times EPSILON, 2.7e-14 J/K, while the
    // tolerance of an inventory near 0.1 J/K is 2.8e-15 J/K. (Measured
    // with both bodies textbook, the defect was 1.9e-14 J/K.)
    let Pair { mut system, .. } = pair(|_, cold, _| cold.textbook = true);
    let refused = with_cx(|cx| {
        system.step(DT_S, cx).expect("the first step is exact");
        system.step(DT_S, cx).expect_err("the second is refused")
    });
    assert!(
        matches!(
            &refused,
            SystemError::Contract {
                node,
                violation: ContractViolation::EntropyBalance { .. },
            } if node == "cold"
        ),
        "{refused:?}"
    );
    assert_eq!(system.steps(), 1);
}

#[test]
fn a_steady_defect_inside_the_allowance_is_booked_signed_and_absolute() {
    // The hot body's entropy stays between 4.6 and 14.4 J/K, so its smallest
    // allowance is about 128 EPSILON 4.6 J/K = 1.3e-13 J/K. A generation
    // short by 4e-14 J/K at every step passes, and leaves a defect of that
    // size, always of one sign, which the ledger books in both accounts.
    const BIAS: f64 = -4.0e-14;
    let mut pair = pair(|hot, _, _| hot.bias_j_per_k = BIAS);
    run(&mut pair.system, 200);
    let account = entropy(pair.system.ledger(), "hot").expect("declared");
    assert_eq!(account.audited_steps, 200);
    let (signed, absolute) = (
        account.balance_defect_j_per_k.value(),
        account.absolute_defect_j_per_k.value(),
    );
    let booked = 200.0 * BIAS;
    assert!(signed <= 0.5 * booked, "{signed:e}");
    assert!(absolute >= -0.5 * booked, "{absolute:e}");
    assert!(absolute >= signed.abs(), "{absolute:e} {signed:e}");
}

#[test]
fn a_refused_measurement_is_named_in_the_measure_phase() {
    let Pair {
        mut system, cold, ..
    } = pair(|_, cold, _| cold.refuses = true);
    let before = system.ledger().clone();
    let refused = with_cx(|cx| system.step(DT_S, cx)).expect_err("refused");
    assert!(
        matches!(
            &refused,
            SystemError::Component {
                node,
                phase: Phase::Measure,
                ..
            } if node == "cold"
        ),
        "{refused:?}"
    );
    assert_eq!(system.steps(), 0);
    assert_eq!(*system.state::<f64>(cold).expect("state"), COLD.1);
    assert_eq!(*system.ledger(), before, "nothing is published");
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

/// The pause boundary of a run whose one worker drained.
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

#[test]
fn the_entropy_accounts_rewind_with_a_checkpoint_and_come_back_from_a_snapshot() {
    let mut first = pair(|_, _, _| {});
    run(&mut first.system, 50);
    let checkpoint = first.system.checkpoint().expect("checkpoint");
    let sealed = first
        .system
        .snapshot(&context(&first.system), limits(), || false)
        .expect("sealed");
    run(&mut first.system, 50);
    let after = first.system.ledger().clone();
    assert_eq!(
        entropy(&after, "wall").expect("declared").audited_steps,
        100
    );

    first.system.restore(&checkpoint).expect("restored");
    assert_eq!(
        entropy(first.system.ledger(), "wall")
            .expect("declared")
            .audited_steps,
        50
    );
    run(&mut first.system, 50);
    assert_eq!(*first.system.ledger(), after, "a checkpoint rewinds them");

    let opened =
        SystemImage::unseal_v2_expected(sealed.bytes(), &sealed.expectation(), limits(), || false)
            .expect("opened");
    let mut second = pair(|_, _, _| {});
    second.system.resume(&opened).expect("resumed");
    assert_eq!(
        entropy(second.system.ledger(), "hot")
            .expect("declared")
            .audited_steps,
        50
    );
    run(&mut second.system, 50);
    assert_eq!(*second.system.ledger(), after, "a snapshot carries them");
}
