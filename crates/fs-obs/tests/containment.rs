//! Local execution-containment battery (i94v.7.3.2): deterministic
//! parentage under reorder/duplicate/interleave, aliasing refusal, lost
//! parents as explicit gaps, wrong-root embedding refusal, and wire-valid
//! JSONL projection.

use fs_obs::containment::{
    AttemptId, AttemptTree, CampaignRunId, CaseId, ContainmentContext, ContainmentError,
    ContainmentRecord, DsrRunId, ExecutionOpId, ExecutionScopeId, GlobalDagEmbedder, IdError,
    Ingest, JourneyId, LocalNodeId, LocalParent, SealedAttemptTree, ShardId, TileId,
    check_embedding_root,
};
use fs_obs::{lint_failure_record, validate_line};

fn ctx() -> ContainmentContext {
    ContainmentContext {
        dsr_run: Some(DsrRunId::new("dsr-7").unwrap()),
        campaign_run: Some(CampaignRunId::new("camp-1").unwrap()),
        shard: Some(ShardId::new("shard-3").unwrap()),
        journey: Some(JourneyId::new("thermal-fatigue").unwrap()),
        case: Some(CaseId::new("case-12").unwrap()),
    }
}

fn op(name: &str, parent: LocalParent, seq: u64) -> ContainmentRecord {
    ContainmentRecord {
        node: LocalNodeId::Op(ExecutionOpId::new(name).unwrap()),
        parent,
        seq,
        context: ctx(),
    }
}

fn scope(name: &str, parent: LocalParent, seq: u64) -> ContainmentRecord {
    ContainmentRecord {
        node: LocalNodeId::Scope(ExecutionScopeId::new(name).unwrap()),
        parent,
        seq,
        context: ctx(),
    }
}

fn tile(name: &str, parent: LocalParent, seq: u64) -> ContainmentRecord {
    ContainmentRecord {
        node: LocalNodeId::Tile(TileId::new(name).unwrap()),
        parent,
        seq,
        context: ctx(),
    }
}

fn node(record: &ContainmentRecord) -> LocalParent {
    LocalParent::Node(record.node.clone())
}

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-obs/containment\",\"case\":\"{case}\",\
         \"verdict\":\"pass\",\"detail\":\"{detail}\"}}"
    );
}

#[test]
fn typed_ids_validate_and_roles_never_confuse() {
    assert!(matches!(
        ExecutionOpId::new(""),
        Err(IdError {
            role: "ExecutionOpId",
            ..
        })
    ));
    assert!(matches!(TileId::new("has space"), Err(IdError { .. })));
    assert!(matches!(
        AttemptId::new("a".repeat(257)),
        Err(IdError { .. })
    ));

    // Same raw text under different roles is a DIFFERENT node identity.
    let as_op = LocalNodeId::Op(ExecutionOpId::new("x-1").unwrap());
    let as_tile = LocalNodeId::Tile(TileId::new("x-1").unwrap());
    assert_ne!(as_op, as_tile, "role tag is part of node identity");

    let root = AttemptId::new("attempt-9").unwrap();
    let mut tree = AttemptTree::new(root);
    tree.ingest(op("x-1", LocalParent::AttemptRoot, 0)).unwrap();
    tree.ingest(tile("x-1", LocalParent::AttemptRoot, 1))
        .expect("same raw text under a different role is not a redelivery");
    verdict(
        "typed-ids",
        "construction validates and role tags keep identical raw text non-confusable",
    );
}

#[test]
fn reordered_and_duplicate_delivery_converge_to_one_deterministic_tree() {
    let root = AttemptId::new("attempt-9").unwrap();
    let solve = op("solve", LocalParent::AttemptRoot, 0);
    let kernel = scope("kernel-lbm", node(&solve), 0);
    let t42 = tile("tile-42", node(&kernel), 0);
    let t43 = tile("tile-43", node(&kernel), 1);

    // In-order delivery.
    let mut a = AttemptTree::new(root.clone());
    for r in [&solve, &kernel, &t42, &t43] {
        assert_eq!(a.ingest(r.clone()).unwrap(), Ingest::Admitted);
    }

    // Fully reversed delivery buffers, then drains as parents arrive; plus
    // byte-identical redelivery of every record is idempotent.
    let mut b = AttemptTree::new(root);
    assert_eq!(b.ingest(t43.clone()).unwrap(), Ingest::Buffered);
    assert_eq!(b.ingest(t42.clone()).unwrap(), Ingest::Buffered);
    assert_eq!(b.ingest(kernel.clone()).unwrap(), Ingest::Buffered);
    assert_eq!(b.ingest(solve.clone()).unwrap(), Ingest::Admitted);
    for r in [&solve, &kernel, &t42, &t43] {
        assert_eq!(b.ingest(r.clone()).unwrap(), Ingest::Duplicate);
    }

    let a = a.seal();
    let b = b.seal();
    assert_eq!(a, b, "delivery order cannot move the sealed tree");
    assert!(a.gaps().is_empty(), "complete lineage seals gap-free");
    assert_eq!(a.nodes().len(), 4);
    verdict(
        "deterministic-parentage",
        "reversed + duplicated delivery seals bit-identically to in-order delivery",
    );
}

#[test]
fn aliasing_self_parentage_and_cycles_refuse() {
    let root = AttemptId::new("attempt-9").unwrap();
    let mut tree = AttemptTree::new(root);
    let solve = op("solve", LocalParent::AttemptRoot, 0);
    tree.ingest(solve.clone()).unwrap();

    // Same node id, different content: aliasing, refused.
    let mut aliased = solve.clone();
    aliased.seq = 5;
    assert!(matches!(
        tree.ingest(aliased),
        Err(ContainmentError::ConflictingRedelivery { .. })
    ));

    // Self-parentage refused.
    let selfish = scope(
        "s",
        LocalParent::Node(LocalNodeId::Scope(ExecutionScopeId::new("s").unwrap())),
        0,
    );
    assert!(matches!(
        tree.ingest(selfish),
        Err(ContainmentError::SelfParent { .. })
    ));

    // Two buffered records forming a cycle: the second refuses.
    let a = scope(
        "cycle-a",
        LocalParent::Node(LocalNodeId::Scope(
            ExecutionScopeId::new("cycle-b").unwrap(),
        )),
        0,
    );
    let b = scope(
        "cycle-b",
        LocalParent::Node(LocalNodeId::Scope(
            ExecutionScopeId::new("cycle-a").unwrap(),
        )),
        0,
    );
    assert_eq!(tree.ingest(a).unwrap(), Ingest::Buffered);
    assert!(matches!(
        tree.ingest(b),
        Err(ContainmentError::Cycle { .. })
    ));
    verdict(
        "refusals",
        "aliased redelivery, self-parentage, and containment cycles all refuse",
    );
}

#[test]
fn context_edges_never_stand_in_for_parentage() {
    // A record whose only relation is contextual (campaign/shard/journey)
    // still needs its primary parent; context cannot flatten into lineage.
    let root = AttemptId::new("attempt-9").unwrap();
    let mut tree = AttemptTree::new(root);
    let orphan = tile(
        "tile-1",
        LocalParent::Node(LocalNodeId::Scope(
            ExecutionScopeId::new("never-arrives").unwrap(),
        )),
        0,
    );
    assert_eq!(tree.ingest(orphan).unwrap(), Ingest::Buffered);
    let sealed = tree.seal();
    assert!(
        sealed.nodes().is_empty(),
        "rich context must not admit a node whose containment parent is unknown"
    );
    assert_eq!(sealed.gaps().len(), 1);
    let gap = &sealed.gaps()[0];
    assert_eq!(gap.node, "tile-1");
    assert_eq!(gap.missing_parent, "never-arrives");
    verdict(
        "context-vs-parentage",
        "contextual edges cannot alias containment; lost parent seals as an explicit gap",
    );
}

#[test]
fn worker_handoff_and_cancellation_races_stay_deterministic_and_honest() {
    // Worker handoff: worker-1 builds the solve subtree, hands off, and
    // worker-2 continues under the SAME propagated attempt root. Whatever
    // interleaving the two record streams arrive in, the sealed tree is
    // identical.
    let root = AttemptId::new("attempt-9").unwrap();
    let solve = op("solve", LocalParent::AttemptRoot, 0);
    let k1 = scope("kernel-lbm", node(&solve), 0);
    let w1_records = vec![solve.clone(), k1.clone(), tile("tile-1", node(&k1), 0)];
    let k2 = scope("kernel-adjoint", node(&solve), 1);
    let w2_records = vec![k2.clone(), tile("tile-2", node(&k2), 0)];

    let interleavings: Vec<Vec<usize>> = vec![
        vec![0, 0, 0, 1, 1], // worker-1 fully first
        vec![1, 1, 0, 0, 0], // worker-2's stream arrives first (reordered)
        vec![0, 1, 0, 1, 0], // alternating
    ];
    let mut sealed = Vec::new();
    for order in &interleavings {
        let mut tree = AttemptTree::new(root.clone());
        let (mut i1, mut i2) = (0usize, 0usize);
        for which in order {
            let record = if *which == 0 {
                let r = w1_records[i1].clone();
                i1 += 1;
                r
            } else {
                let r = w2_records[i2].clone();
                i2 += 1;
                r
            };
            tree.ingest(record).expect("handoff records admit");
        }
        sealed.push(tree.seal());
    }
    assert!(
        sealed.windows(2).all(|w| w[0] == w[1]),
        "worker handoff seals identically under every interleaving"
    );
    assert!(sealed[0].gaps().is_empty());

    // Cancellation race: the scope is cancelled while its tiles are in
    // flight, so some records never arrive. The seal never invents
    // lineage — the missing parent is an explicit gap, and a tile whose
    // parent DID arrive still admits even if delivered after the cancel.
    let mut tree = AttemptTree::new(root);
    tree.ingest(solve.clone()).unwrap();
    tree.ingest(k1.clone()).unwrap();
    // k2 was cancelled before its record flushed; its tile raced ahead.
    tree.ingest(tile("tile-2", node(&k2), 0)).unwrap();
    // tile-1's parent k1 arrived, so late delivery still admits.
    tree.ingest(tile("tile-1", node(&k1), 0)).unwrap();
    let sealed = tree.seal();
    assert_eq!(sealed.nodes().len(), 3, "solve, kernel-lbm, tile-1");
    assert_eq!(sealed.gaps().len(), 1, "tile-2's cancelled parent is a gap");
    assert_eq!(sealed.gaps()[0].node, "tile-2");
    assert_eq!(sealed.gaps()[0].missing_parent, "kernel-adjoint");
    verdict(
        "handoff-and-cancel-race",
        "handoff interleavings seal bit-identically; cancellation races surface as explicit gaps",
    );
}

#[test]
fn embedding_refuses_the_wrong_attempt_root() {
    struct Recorder {
        embedded: usize,
        under: AttemptId,
    }
    impl GlobalDagEmbedder for Recorder {
        type Error = ContainmentError;
        fn embed(&mut self, tree: &SealedAttemptTree) -> Result<(), Self::Error> {
            check_embedding_root(tree, &self.under)?;
            self.embedded += 1;
            Ok(())
        }
    }

    let mut tree = AttemptTree::new(AttemptId::new("attempt-9").unwrap());
    tree.ingest(op("solve", LocalParent::AttemptRoot, 0))
        .unwrap();
    let sealed = tree.seal();

    let mut wrong = Recorder {
        embedded: 0,
        under: AttemptId::new("attempt-10").unwrap(),
    };
    assert!(matches!(
        wrong.embed(&sealed),
        Err(ContainmentError::WrongAttemptRoot { .. })
    ));
    assert_eq!(wrong.embedded, 0);

    let mut right = Recorder {
        embedded: 0,
        under: AttemptId::new("attempt-9").unwrap(),
    };
    right.embed(&sealed).unwrap();
    assert_eq!(right.embedded, 1);
    verdict(
        "embedding-root",
        "V.3.8 embedding seam refuses a mismatched attempt root and admits the true one",
    );
}

#[test]
fn sealed_trees_project_to_wire_valid_lossless_jsonl() {
    let root = AttemptId::new("attempt-9").unwrap();
    let mut tree = AttemptTree::new(root.clone());
    let solve = op("solve", LocalParent::AttemptRoot, 0);
    let kernel = scope("kernel-lbm", node(&solve), 0);
    tree.ingest(solve).unwrap();
    tree.ingest(kernel.clone()).unwrap();
    // One lost parent so the gap projection is exercised too.
    tree.ingest(tile(
        "tile-9",
        LocalParent::Node(LocalNodeId::Scope(ExecutionScopeId::new("gone").unwrap())),
        0,
    ))
    .unwrap();
    let sealed = tree.seal();
    let events = sealed.to_events("study-x", "attempt-9");
    assert_eq!(events.len(), 3, "two nodes + one gap");
    for event in &events {
        let line = event.to_jsonl();
        validate_line(&line).expect("containment projection stays wire-valid");
        lint_failure_record(event).expect("containment records pass the failure lint");
    }
    let joined: Vec<String> = events.iter().map(fs_obs::Event::to_jsonl).collect();
    assert!(
        joined
            .iter()
            .any(|l| l.contains("\"kind\":\"containment_node\"")
                && l.contains("\"shard\":\"shard-3\"")),
        "node lines carry the full typed context"
    );
    assert!(
        joined
            .iter()
            .any(|l| l.contains("\"kind\":\"containment_gap\"")
                && l.contains("\"missing_parent\":\"gone\"")),
        "gap lines name the missing parent"
    );
    verdict(
        "jsonl-projection",
        "sealed nodes and gaps project to wire-valid typed containment lines with full context",
    );
}
