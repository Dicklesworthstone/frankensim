//! I01.2 battery (bead i94v.1.1.2): ConstitutiveGraph nodes adapted
//! into the compiler — provenance survives lowering exactly, tangent
//! claims are verified evidence, unsupported lanes route typed, and
//! state schemas refuse drift. Runs only with `constitutive-graph`
//! (ambition [F], default-off).
#![cfg(feature = "constitutive-graph")]

use fs_evidence::ValidityDomain;
use fs_matdb::LawId;
use fs_material::graph::{
    Differentiability, EnergyBehavior, GraphError, LawNode, NodeDeclaration, NodeOutput, Port,
    TimeParity,
};
use fs_opdsl::constitutive::{
    AdaptError, BoundConstitutiveNode, DifferentiabilityTag, PotentialChart, TangentLane,
};
use fs_qty::Dims;

fn port(name: &str, dims: Dims) -> Port {
    Port {
        name: name.to_string(),
        dims,
        parity: TimeParity::Even,
    }
}

fn declaration(
    law: &str,
    differentiability: Differentiability,
    energy: EnergyBehavior,
    tangent_claimed: bool,
    state_slots: Vec<String>,
) -> NodeDeclaration {
    NodeDeclaration {
        law: LawId(law.to_string()),
        law_version: 2,
        role: fs_material::graph::NodeRole::BulkTransport,
        inputs: vec![port("gradient", Dims([-1, 0, 0, 1, 0, 0]))],
        outputs: vec![port("flux", Dims([0, 1, -3, 0, 0, 0]))],
        state_slots,
        state_schema_version: 3,
        calibration: ValidityDomain::unconstrained(),
        differentiability,
        energy,
        tangent_claimed,
    }
}

/// Honest linear conduction: flux = -k * gradient, exact tangent -k.
struct HonestFourier {
    decl: NodeDeclaration,
    k: f64,
}

impl HonestFourier {
    fn new() -> HonestFourier {
        HonestFourier {
            decl: declaration(
                "fourier-iso",
                Differentiability::Smooth,
                EnergyBehavior::NonNegativeDissipation,
                true,
                Vec::new(),
            ),
            k: 2.5,
        }
    }
}

impl LawNode for HonestFourier {
    fn declaration(&self) -> &NodeDeclaration {
        &self.decl
    }
    fn evaluate(&self, _state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        Ok(NodeOutput {
            outputs: vec![-self.k * inputs[0]],
            next_state: Vec::new(),
            dissipation_rate: Some(self.k * inputs[0] * inputs[0]),
        })
    }
    fn tangent(&self, _state: &[f64], _inputs: &[f64]) -> Option<Vec<f64>> {
        Some(vec![-self.k])
    }
}

/// A LIAR: claims a consistent tangent, supplies the wrong sign.
struct LyingFourier(HonestFourier);

impl LawNode for LyingFourier {
    fn declaration(&self) -> &NodeDeclaration {
        &self.0.decl
    }
    fn evaluate(&self, state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        self.0.evaluate(state, inputs)
    }
    fn tangent(&self, _state: &[f64], _inputs: &[f64]) -> Option<Vec<f64>> {
        Some(vec![self.0.k]) // sign-flipped lie
    }
}

/// A nonlinear memory node whose state update is order-sensitive:
/// s' = s + input * (1 + s²) — the quadratic term breaks the
/// commutativity that a multiplicative accumulator would have.
struct NonlinearMemory {
    decl: NodeDeclaration,
}

impl NonlinearMemory {
    fn new() -> NonlinearMemory {
        NonlinearMemory {
            decl: declaration(
                "nonlinear-memory",
                Differentiability::PiecewiseSmooth,
                EnergyBehavior::Empirical,
                true,
                vec!["accumulated".to_string()],
            ),
        }
    }
}

impl LawNode for NonlinearMemory {
    fn declaration(&self) -> &NodeDeclaration {
        &self.decl
    }
    fn evaluate(&self, state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        let next = state[0] + inputs[0] * (1.0 + state[0] * state[0]);
        Ok(NodeOutput {
            outputs: vec![next],
            next_state: vec![next],
            dissipation_rate: None,
        })
    }
    fn tangent(&self, state: &[f64], _inputs: &[f64]) -> Option<Vec<f64>> {
        Some(vec![1.0 + state[0]])
    }
}

/// A smooth storage law with a free energy: psi = 0.5 * x^2.
struct QuadraticStorage {
    decl: NodeDeclaration,
}

impl QuadraticStorage {
    fn new() -> QuadraticStorage {
        QuadraticStorage {
            decl: declaration(
                "quadratic-storage",
                Differentiability::Smooth,
                EnergyBehavior::FreeEnergyStorage,
                true,
                Vec::new(),
            ),
        }
    }
}

impl LawNode for QuadraticStorage {
    fn declaration(&self) -> &NodeDeclaration {
        &self.decl
    }
    fn evaluate(&self, _state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        Ok(NodeOutput {
            outputs: vec![inputs[0]],
            next_state: Vec::new(),
            dissipation_rate: None,
        })
    }
    fn tangent(&self, _state: &[f64], _inputs: &[f64]) -> Option<Vec<f64>> {
        Some(vec![1.0])
    }
    fn free_energy(&self, _state: &[f64], inputs: &[f64]) -> Option<f64> {
        Some(0.5 * inputs[0] * inputs[0])
    }
}

/// Provenance survives lowering EXACTLY: every declared identity field
/// reappears verbatim in the compiler-owned receipt, and mutating any
/// one of them is detectable by receipt inequality.
#[test]
fn material_provenance_survives_binding_exactly() {
    let node = HonestFourier::new();
    let bound = BoundConstitutiveNode::bind(&node, None).expect("honest node binds");
    let receipt = bound.provenance().clone();
    assert_eq!(receipt.law, "fourier-iso");
    assert_eq!(receipt.law_version, 2);
    assert_eq!(receipt.state_schema_version, 3);
    assert_eq!(receipt.state_slots, 0);
    assert_eq!(receipt.input_dims, vec![Dims([-1, 0, 0, 1, 0, 0])]);
    assert_eq!(receipt.output_dims, vec![Dims([0, 1, -3, 0, 0, 0])]);
    assert_eq!(receipt.differentiability, DifferentiabilityTag::Smooth);
    assert_eq!(receipt.potential_chart, PotentialChart::Dissipation);

    // Receipt mutations are detectable (the missing-receipt mutation
    // lane): any single-field drift breaks equality.
    for mutate in 0..3 {
        let mut mutated = receipt.clone();
        match mutate {
            0 => mutated.law_version += 1,
            1 => mutated.state_schema_version += 1,
            _ => mutated.law.push('x'),
        }
        assert_ne!(mutated, receipt, "mutation {mutate} must be visible");
    }
}

/// A supplied tangent is evidence: the honest node earns the
/// Consistent lane; the liar is refused at binding with a typed error.
#[test]
fn tangent_claims_are_verified_evidence_not_authority() {
    let honest = HonestFourier::new();
    let bound = BoundConstitutiveNode::bind(&honest, None).expect("honest tangent verifies");
    assert_eq!(bound.lane(), TangentLane::Consistent);

    let liar = LyingFourier(HonestFourier::new());
    let refused = BoundConstitutiveNode::bind(&liar, None);
    assert!(
        matches!(refused, Err(AdaptError::TangentEvidenceRejected { ref law, .. }) if law == "fourier-iso"),
        "the lying tangent must be rejected at binding: {refused:?}"
    );
}

/// Routing is typed: a NonSmooth declaration claiming a tangent
/// refuses; an unclaimed tangent routes DerivativeFree, and asking
/// that lane for a tangent refuses instead of differentiating.
#[test]
fn unsupported_differentiability_routes_typed() {
    let mut nonsmooth = HonestFourier::new();
    nonsmooth.decl.differentiability = Differentiability::NonSmooth;
    let refused = BoundConstitutiveNode::bind(&nonsmooth, None);
    assert!(matches!(
        refused,
        Err(AdaptError::UnsupportedDifferentiability {
            requested: "consistent-tangent",
            ..
        })
    ));

    let mut unclaimed = HonestFourier::new();
    unclaimed.decl.tangent_claimed = false;
    let bound = BoundConstitutiveNode::bind(&unclaimed, None).expect("derivative-free binds");
    assert_eq!(bound.lane(), TangentLane::DerivativeFree);
    assert!(matches!(
        bound.tangent(&[1.0]),
        Err(AdaptError::UnsupportedDifferentiability {
            requested: "tangent",
            ..
        })
    ));
}

/// State-owning nodes demand explicit initialization; the state codec
/// refuses schema drift and round-trips under the bound version.
#[test]
fn state_initialization_and_codec_migration_refuse_drift() {
    let memory = NonlinearMemory::new();
    let refused = BoundConstitutiveNode::bind(&memory, None);
    assert!(matches!(
        refused,
        Err(AdaptError::MissingStateInitialization { state_slots: 1, .. })
    ));

    let mut bound =
        BoundConstitutiveNode::bind(&memory, Some(&[0.25])).expect("declared state binds");
    assert_eq!(bound.state(), &[0.25]);

    // Round trip under the bound schema version.
    bound
        .restore_state(3, &[0.5])
        .expect("same-version restore admits");
    assert_eq!(bound.state(), &[0.5]);

    // Drifted schema version refuses with both versions named.
    let drift = bound.restore_state(4, &[0.5]);
    assert!(matches!(
        drift,
        Err(AdaptError::StateSchemaDrift {
            bound: 3,
            offered: 4,
            ..
        })
    ));
}

/// History order is faithfully sequenced: the nonlinear memory node
/// reaches genuinely different states under reordered input histories
/// (the adapter must not commute what physics does not).
#[test]
fn history_reorder_reaches_different_states() {
    let memory = NonlinearMemory::new();
    let mut forward = BoundConstitutiveNode::bind(&memory, Some(&[0.0])).expect("binds");
    forward.evaluate(&[0.1]).expect("step 1");
    forward.evaluate(&[0.4]).expect("step 2");

    let mut reordered = BoundConstitutiveNode::bind(&memory, Some(&[0.0])).expect("binds");
    reordered.evaluate(&[0.4]).expect("step 1");
    reordered.evaluate(&[0.1]).expect("step 2");

    assert_ne!(
        forward.state()[0].to_bits(),
        reordered.state()[0].to_bits(),
        "nonlinear history must be order-sensitive through the adapter"
    );

    // And replay is deterministic: the same history is bitwise stable.
    let mut replay = BoundConstitutiveNode::bind(&memory, Some(&[0.0])).expect("binds");
    replay.evaluate(&[0.1]).expect("step 1");
    replay.evaluate(&[0.4]).expect("step 2");
    assert_eq!(forward.state()[0].to_bits(), replay.state()[0].to_bits());
}

/// The VJP is the exact transpose contraction of the verified tangent,
/// and the potential chart gates energy access: storage laws expose
/// free energy; dissipative and empirical charts do not.
#[test]
fn vjp_and_potential_chart_route_exactly() {
    let node = HonestFourier::new();
    let bound = BoundConstitutiveNode::bind(&node, None).expect("binds");
    let vjp = bound.vjp(&[0.7], &[2.0]).expect("vjp on verified lane");
    assert_eq!(vjp.len(), 1);
    assert!((vjp[0] - (2.0 * -2.5)).abs() < 1e-12);
    assert!(matches!(
        bound.vjp(&[0.7], &[1.0, 1.0]),
        Err(AdaptError::Evaluation { .. })
    ));
    // Dissipative chart: no free energy exposed.
    assert!(bound.free_energy(&[0.7]).is_none());

    let storage = QuadraticStorage::new();
    let bound_storage = BoundConstitutiveNode::bind(&storage, None).expect("binds");
    let psi = bound_storage.free_energy(&[3.0]).expect("storage chart");
    assert!((psi - 4.5).abs() < 1e-12);
    // Dissipation rate reported through evaluation survives.
    let mut flux = BoundConstitutiveNode::bind(&node, None).expect("binds");
    let (outputs, dissipation) = flux.evaluate(&[2.0]).expect("evaluates");
    assert!((outputs[0] + 5.0).abs() < 1e-12);
    assert!(dissipation.expect("reported") >= 0.0, "second-law fixture");
}

/// The hand-written escape hatch binds under the same gates and
/// retains its no-generated-consistency marker.
#[test]
fn hand_written_escape_hatch_retains_its_marker() {
    let node = HonestFourier::new();
    let bound = BoundConstitutiveNode::bind_hand_written(&node, None).expect("escape hatch binds");
    let escape = bound.escape().expect("marker retained");
    assert!(escape.no_generated_consistency_claim);
    // The same evidence gate applies on the escape path too.
    let liar = LyingFourier(HonestFourier::new());
    assert!(matches!(
        BoundConstitutiveNode::bind_hand_written(&liar, None),
        Err(AdaptError::TangentEvidenceRejected { .. })
    ));
}

// ---------------------------------------------------------------------------
// Slice 2: batched evaluation under Cx (request-drain-finalize).
// ---------------------------------------------------------------------------

use std::collections::BTreeMap;

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_opdsl::constitutive::{BatchPoint, BatchRun, evaluate_batch};

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x1012,
                kernel_id: 5,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn conduction_graph() -> fs_material::graph::ConstitutiveGraph {
    let mut graph = fs_material::graph::ConstitutiveGraph::new();
    graph
        .add_node("fourier", Box::new(HonestFourier::new()))
        .expect("admits");
    graph
}

fn batch_of(n: usize) -> Vec<BatchPoint> {
    let (version, empty_state) = {
        let graph = conduction_graph();
        let schema = graph.state_schema();
        (schema.version(), vec![0.0f64; schema.total_slots()])
    };
    (0..n)
        .map(|i| {
            let mut external = BTreeMap::new();
            external.insert(
                ("fourier".to_string(), "gradient".to_string()),
                0.1 * (i as f64 + 1.0),
            );
            BatchPoint {
                state_version: version,
                state: empty_state.clone(),
                external,
            }
        })
        .collect()
}

/// A healthy batch completes deterministically: every point evaluated
/// in order, values exact, replay bitwise.
#[test]
fn batched_evaluation_completes_deterministically() {
    let graph = conduction_graph();
    let batch = batch_of(40);
    let run = with_cx(false, |cx| evaluate_batch(&graph, &batch, cx)).expect("healthy batch");
    let BatchRun::Complete {
        outputs,
        points_evaluated,
    } = run
    else {
        panic!("uncancelled batch must complete");
    };
    assert_eq!(points_evaluated, 40);
    for (i, point) in outputs.iter().enumerate() {
        let expected = -2.5 * 0.1 * (i as f64 + 1.0);
        let got = point.outputs[&("fourier".to_string(), "flux".to_string())];
        assert!(
            (got - expected).abs() < 1e-12,
            "point {i}: {got} vs {expected}"
        );
        assert!(point.total_dissipation >= 0.0);
    }
    let replay = with_cx(false, |cx| evaluate_batch(&graph, &batch, cx)).expect("replay");
    let BatchRun::Complete {
        outputs: replayed, ..
    } = replay
    else {
        panic!("replay completes");
    };
    for (a, b) in outputs.iter().zip(&replayed) {
        for (key, value) in &a.outputs {
            assert_eq!(value.to_bits(), b.outputs[key].to_bits(), "bitwise replay");
        }
    }
}

/// A pre-cancelled context drains at the FIRST point boundary with an
/// empty completed prefix and resume cursor 0; resuming the remainder
/// under a healthy context is bitwise-equivalent to the uncancelled
/// run (request-drain-finalize + deterministic resume).
#[test]
fn cancellation_drains_at_point_boundaries_and_resumes_bitwise() {
    let graph = conduction_graph();
    let batch = batch_of(40);

    let cancelled = with_cx(true, |cx| evaluate_batch(&graph, &batch, cx)).expect("drains");
    let BatchRun::Cancelled {
        completed,
        resume_from,
    } = cancelled
    else {
        panic!("a pre-cancelled context must drain, not complete");
    };
    assert_eq!(resume_from, 0, "pre-cancel drains before the first point");
    assert!(completed.is_empty());

    // Resume: evaluate the remainder healthily and splice.
    let resumed = with_cx(false, |cx| {
        evaluate_batch(&graph, &batch[resume_from..], cx)
    })
    .expect("resume evaluates");
    let BatchRun::Complete { outputs: tail, .. } = resumed else {
        panic!("resume completes");
    };
    let full = with_cx(false, |cx| evaluate_batch(&graph, &batch, cx)).expect("oracle");
    let BatchRun::Complete {
        outputs: oracle, ..
    } = full
    else {
        panic!("oracle completes");
    };
    assert_eq!(completed.len() + tail.len(), oracle.len());
    for (spliced, reference) in completed.iter().chain(&tail).zip(&oracle) {
        for (key, value) in &reference.outputs {
            assert_eq!(
                spliced.outputs[key].to_bits(),
                value.to_bits(),
                "resumed run must be bitwise-equivalent to the uncancelled run"
            );
        }
    }
}

/// A defective point refuses the WHOLE batch typed (its index named);
/// no partial result set is published past a refusal.
#[test]
fn defective_batch_points_refuse_typed_with_index() {
    let graph = conduction_graph();
    let mut batch = batch_of(3);
    batch[1].external.clear(); // unfed port at point 1
    let refused = with_cx(false, |cx| evaluate_batch(&graph, &batch, cx));
    assert!(
        matches!(&refused, Err(AdaptError::Evaluation { law, .. }) if law == "batch point 1"),
        "the refusal must name the offending point: {refused:?}"
    );
}

/// A law node that requests cancellation during its N-th evaluation —
/// the deterministic mid-batch trip wire.
struct TrippingFourier {
    inner: HonestFourier,
    evals: std::cell::Cell<usize>,
    trip_at: usize,
    gate: &'static CancelGate,
}

impl LawNode for TrippingFourier {
    fn declaration(&self) -> &NodeDeclaration {
        &self.inner.decl
    }
    fn evaluate(&self, state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        let n = self.evals.get() + 1;
        self.evals.set(n);
        if n == self.trip_at {
            self.gate.request();
        }
        self.inner.evaluate(state, inputs)
    }
    fn tangent(&self, state: &[f64], inputs: &[f64]) -> Option<Vec<f64>> {
        self.inner.tangent(state, inputs)
    }
}

/// Deterministic MID-BATCH drain: cancellation requested during point
/// 20 is observed at the NEXT stride boundary (32) — the in-flight
/// stride drains whole, exactly 32 points complete, the resume cursor
/// lands on 32, and the completed prefix is bitwise-identical to the
/// oracle's prefix. No partial point is ever published.
#[test]
fn mid_batch_cancellation_drains_at_the_next_stride_boundary() {
    let gate: &'static CancelGate = Box::leak(Box::new(CancelGate::new()));
    let mut graph = fs_material::graph::ConstitutiveGraph::new();
    graph
        .add_node(
            "fourier",
            Box::new(TrippingFourier {
                inner: HonestFourier::new(),
                evals: std::cell::Cell::new(0),
                trip_at: 21, // fires while evaluating batch point index 20
                gate,
            }),
        )
        .expect("admits");
    let batch = batch_of(40);
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let run = pool
        .scope(|arena| {
            let cx = Cx::new(
                gate,
                arena,
                StreamKey {
                    seed: 0x1012,
                    kernel_id: 5,
                    tile: 0,
                    iteration: 1,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            evaluate_batch(&graph, &batch, &cx)
        })
        .expect("drains, not errors");
    let BatchRun::Cancelled {
        completed,
        resume_from,
    } = run
    else {
        panic!("mid-batch cancellation must drain at a stride boundary");
    };
    assert_eq!(resume_from, 32, "requested during point 20, observed at 32");
    assert_eq!(completed.len(), 32);
    let oracle_graph = conduction_graph();
    let oracle = with_cx(false, |cx| evaluate_batch(&oracle_graph, &batch, cx)).expect("oracle");
    let BatchRun::Complete { outputs, .. } = oracle else {
        panic!("oracle completes");
    };
    for (done, reference) in completed.iter().zip(&outputs) {
        for (key, value) in &reference.outputs {
            assert_eq!(done.outputs[key].to_bits(), value.to_bits());
        }
    }
}
