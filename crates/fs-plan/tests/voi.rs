//! VoI-query conformance (the knh1.6 bead; runs under `voi-queries`).
//! Acceptance: decision sensitivity from cached sweeps (near-free —
//! call counts prove it); rankings by sampled flip-fraction reduction per
//! dollar; the
//! probe menu unifies computational and physical experiments; the
//! ranking surfaces as the query hint and the probe scheduler; myopic
//! one-step only; the prospective-audit kill criterion demotes VoI
//! when recommendations stop outperforming. The wk4m cases add G0 authority
//! construction and G3 tamper/freshness refusals for ledger-backed receipts.
#![cfg(feature = "voi-queries")]

use std::cell::Cell;
use std::cell::RefCell;

use asupersync::types::Budget;
use fs_blake3::{ContentHash, hash_bytes};
use fs_ledger::vcs::Vcs;
use fs_ledger::{EdgeRole, ExecMode, FiveExplicits, Ledger, MAIN_BRANCH, OpOutcome};
use fs_plan::voi::{
    AuditReport, AuditVerdict, AuthenticatedRankedMenu, Cx, DecisionBudget,
    DecisionEvaluationPermit, DecisionOracle, LiveDecision, MAX_VOI_AUDIT_RECORDS,
    MAX_VOI_EVALUATIONS, MAX_VOI_GRID, MAX_VOI_NAME_BYTES, MAX_VOI_NODES, MAX_VOI_PROBES,
    MAX_VOI_WORK_UNITS, MatchedAuditRecord, Probe, ProbeKind, RankedMenu, UncertaintyNode,
    VOI_AUDIT_CONTEXT_IDENTITY_VERSION, VOI_AUDIT_OUTCOME_ARTIFACT_KIND,
    VOI_AUDIT_OUTCOME_ARTIFACT_METADATA, VOI_DECISION_ARTIFACT_KIND, VOI_MODEL_ARTIFACT_KIND,
    VOI_PROBE_CATALOG_ARTIFACT_KIND, VOI_PROBE_CATALOG_ARTIFACT_METADATA,
    VOI_RANKED_MENU_IDENTITY_VERSION, VOI_RANKED_SOURCE_IDENTITY_VERSION, VoiError,
    VoiLedgerContext, VoiLedgerReceipt, VoiScheduler, audit_scheduling as audit_scheduling_scoped,
    hint_for_query, rank_purchases as rank_purchases_impl, record_prospective_audit_receipt,
    record_ranked_menu_receipt, verify_prospective_audit_receipt, verify_ranked_menu_receipt,
    voi_audit_outcome_bytes, voi_decision_artifact_metadata, voi_model_artifact_metadata,
    voi_probe_catalog_bytes,
};

const POLICY_SCOPE: &str = "fixture-policy-v1";
const SNAPSHOT_ID: &str = "fixture-snapshot-v1";
const ORACLE_TEST_EVALUATIONS: usize = 5;
const ORACLE_TEST_WORK_UNITS: u64 = 5;
const AUTHENTICATED_AUDIT_WINS: usize = 16;

const LEDGER_FIXTURE_EXPLICITS: FiveExplicits<'static> = FiveExplicits {
    seed: b"fs-plan-voi-test-fixture-v1",
    versions: r#"{"fs-plan-voi-test-fixture":1}"#,
    budget: r#"{"artifacts":64}"#,
    capability: r#"{"ops":["voi.fixture"]}"#,
};

fn decision_budget() -> DecisionBudget {
    DecisionBudget::new(MAX_VOI_EVALUATIONS, MAX_VOI_WORK_UNITS)
        .expect("valid fixture decision budget")
}

fn rank_purchases_scoped(
    decision: &LiveDecision<'_>,
    nodes: &[UncertaintyNode],
    menu: &[Probe],
    grid: usize,
    policy_scope: &str,
    snapshot_id: &str,
) -> Result<RankedMenu, VoiError> {
    let cx = Cx::for_testing();
    rank_purchases_impl(
        &cx,
        decision,
        nodes,
        menu,
        grid,
        decision_budget(),
        policy_scope,
        snapshot_id,
    )
}

fn rank_purchases(
    decision: &LiveDecision<'_>,
    nodes: &[UncertaintyNode],
    menu: &[Probe],
    grid: usize,
) -> Result<RankedMenu, VoiError> {
    rank_purchases_scoped(decision, nodes, menu, grid, POLICY_SCOPE, SNAPSHOT_ID)
}

fn nominal_verdict(
    decision: &LiveDecision<'_>,
    nodes: &[UncertaintyNode],
) -> Result<bool, VoiError> {
    decision.nominal_verdict(&Cx::for_testing(), nodes, decision_budget())
}

fn flip_probability(
    decision: &LiveDecision<'_>,
    nodes: &[UncertaintyNode],
    node_idx: usize,
    grid: usize,
) -> Result<f64, VoiError> {
    decision.flip_probability(&Cx::for_testing(), nodes, node_idx, grid, decision_budget())
}

fn audit_scheduling(records: &[MatchedAuditRecord]) -> Result<AuditReport, VoiError> {
    audit_scheduling_scoped(POLICY_SCOPE, records)
}

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-plan/voi\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// The live decision: ship when margin = 2·drag_gap − mass_penalty +
/// 0.3·hazard > 0 (a cheap cached surrogate).
fn margin(v: &[f64]) -> f64 {
    2.0 * v[0] - v[1] + 0.3 * v[2]
}

fn nodes() -> Vec<UncertaintyNode> {
    vec![
        // Straddles the boundary: pivotal.
        UncertaintyNode {
            name: "drag-gap".to_string(),
            lo: -0.4,
            hi: 0.6,
            nominal: 0.1,
        },
        // Entirely on one side: irrelevant to the verdict.
        UncertaintyNode {
            name: "mass-penalty".to_string(),
            lo: 0.05,
            hi: 0.15,
            nominal: 0.1,
        },
        // Mildly pivotal.
        UncertaintyNode {
            name: "hazard".to_string(),
            lo: -1.0,
            hi: 1.0,
            nominal: 0.2,
        },
    ]
}

fn probe(name: &str, target: &str, cost: f64, shrink: f64) -> Probe {
    Probe {
        name: name.to_string(),
        target: target.to_string(),
        cost,
        shrink,
        kind: ProbeKind::Computational,
    }
}

fn ranked_menu(probes: &[Probe]) -> RankedMenu {
    ranked_menu_for_snapshot(probes, SNAPSHOT_ID)
}

fn ranked_menu_for_snapshot(probes: &[Probe], snapshot_id: &str) -> RankedMenu {
    let decision = LiveDecision {
        margin: &margin,
        arity: 3,
    };
    rank_purchases_scoped(&decision, &nodes(), probes, 64, POLICY_SCOPE, snapshot_id)
        .expect("valid ranked-menu fixture")
}

fn audit_record(index: usize, recommended_wins: bool) -> MatchedAuditRecord {
    MatchedAuditRecord::new(
        format!("audit-{index:04}"),
        "voi-choice",
        "agent-choice",
        format!("fixture/run-{index:04}"),
        10.0,
        10.0,
        recommended_wins,
        !recommended_wins,
    )
    .expect("valid matched-cost audit fixture")
}

fn winning_scheduler(budget: f64) -> VoiScheduler {
    let records: Vec<_> = (0..128).map(|index| audit_record(index, true)).collect();
    let mut scheduler = VoiScheduler::new(POLICY_SCOPE, budget).expect("valid scheduler fixture");
    for record in records {
        scheduler
            .observe_audit(record)
            .expect("valid prospective audit observation");
    }
    assert_eq!(
        scheduler.audit_report().expect("audit report").verdict(),
        AuditVerdict::KeepScheduling
    );
    scheduler
}

fn scheduler_menu() -> RankedMenu {
    ranked_menu(&[
        probe("a", "drag-gap", 10.0, 0.01),
        probe("b", "drag-gap", 5.0, 0.01),
    ])
}

struct LedgerAuthorityFixture {
    ledger: Ledger,
    vcs: Vcs,
    context: VoiLedgerContext,
    decision_artifact: ContentHash,
    model_artifact: ContentHash,
    catalog_artifact: ContentHash,
    outcomes: Vec<(ContentHash, ContentHash)>,
    menu: Vec<Probe>,
    records: Vec<MatchedAuditRecord>,
}

fn fixture_output(
    ledger: &Ledger,
    operation: i64,
    kind: &str,
    bytes: &[u8],
    meta: &str,
) -> ContentHash {
    let receipt = ledger
        .put_artifact(kind, bytes, Some(meta))
        .expect("store VoI authority fixture artifact");
    ledger
        .link(operation, &receipt.hash, EdgeRole::Out)
        .expect("link VoI authority fixture output");
    receipt.hash
}

fn ledger_authority_fixture() -> LedgerAuthorityFixture {
    let ledger = Ledger::open(":memory:").expect("open VoI authority fixture ledger");
    let menu = vec![
        probe("authenticated-a", "drag-gap", 10.0, 0.01),
        probe("authenticated-b", "drag-gap", 5.0, 0.01),
    ];
    let records = (0..AUTHENTICATED_AUDIT_WINS)
        .map(|index| audit_record(index, true))
        .collect::<Vec<_>>();
    let operation = ledger
        .begin_op_on(
            MAIN_BRANCH,
            ExecMode::Deterministic,
            Some(b"voi-authority-fixture".as_slice()),
            r#"{"op":"voi.fixture"}"#,
            &LEDGER_FIXTURE_EXPLICITS,
            1,
        )
        .expect("begin VoI authority fixture operation");
    let decision_meta =
        voi_decision_artifact_metadata("decision-v1").expect("decision artifact metadata");
    let model_meta = voi_model_artifact_metadata("model-v1").expect("model artifact metadata");
    let decision_artifact = fixture_output(
        &ledger,
        operation,
        VOI_DECISION_ARTIFACT_KIND,
        b"decision-v1-canonical-bytes",
        &decision_meta,
    );
    let model_artifact = fixture_output(
        &ledger,
        operation,
        VOI_MODEL_ARTIFACT_KIND,
        b"model-v1-canonical-bytes",
        &model_meta,
    );
    let catalog = voi_probe_catalog_bytes(&menu).expect("canonical complete probe catalog");
    let catalog_artifact = fixture_output(
        &ledger,
        operation,
        VOI_PROBE_CATALOG_ARTIFACT_KIND,
        &catalog,
        VOI_PROBE_CATALOG_ARTIFACT_METADATA,
    );
    let outcomes = records
        .iter()
        .map(|record| {
            let recommended =
                voi_audit_outcome_bytes(record, true).expect("recommended outcome bytes");
            let alternative =
                voi_audit_outcome_bytes(record, false).expect("alternative outcome bytes");
            (
                fixture_output(
                    &ledger,
                    operation,
                    VOI_AUDIT_OUTCOME_ARTIFACT_KIND,
                    &recommended,
                    VOI_AUDIT_OUTCOME_ARTIFACT_METADATA,
                ),
                fixture_output(
                    &ledger,
                    operation,
                    VOI_AUDIT_OUTCOME_ARTIFACT_KIND,
                    &alternative,
                    VOI_AUDIT_OUTCOME_ARTIFACT_METADATA,
                ),
            )
        })
        .collect::<Vec<_>>();
    ledger
        .finish_op(operation, OpOutcome::Ok, None, 2)
        .expect("finish VoI authority fixture operation");
    let mut vcs = Vcs::new();
    let snapshot = vcs
        .commit(&ledger, MAIN_BRANCH)
        .expect("commit VoI authority fixture snapshot");
    let context = VoiLedgerContext::new(
        "session-v1",
        "decision-v1",
        "model-v1",
        POLICY_SCOPE,
        MAIN_BRANCH,
        snapshot.root,
        10,
        20,
    )
    .expect("valid VoI ledger context");
    LedgerAuthorityFixture {
        ledger,
        vcs,
        context,
        decision_artifact,
        model_artifact,
        catalog_artifact,
        outcomes,
        menu,
        records,
    }
}

fn authenticate_ranked_menu(
    fixture: &LedgerAuthorityFixture,
) -> (VoiLedgerReceipt, AuthenticatedRankedMenu) {
    let ranked = ranked_menu_for_snapshot(&fixture.menu, &fixture.context.snapshot_root().to_hex());
    let receipt = record_ranked_menu_receipt(
        &fixture.ledger,
        &fixture.vcs,
        &ranked,
        &fixture.context,
        fixture.decision_artifact,
        fixture.model_artifact,
        fixture.catalog_artifact,
    )
    .expect("record ranked-menu authority receipt");
    let replayed = record_ranked_menu_receipt(
        &fixture.ledger,
        &fixture.vcs,
        &ranked,
        &fixture.context,
        fixture.decision_artifact,
        fixture.model_artifact,
        fixture.catalog_artifact,
    )
    .expect("exact ranked receipt retry is idempotent");
    assert_eq!(replayed, receipt);
    let authenticated = verify_ranked_menu_receipt(
        &fixture.ledger,
        &fixture.vcs,
        ranked,
        fixture.context.clone(),
        fixture.decision_artifact,
        fixture.model_artifact,
        fixture.catalog_artifact,
        receipt,
        10,
    )
    .expect("independently verify ranked-menu authority receipt");
    (receipt, authenticated)
}

fn authenticate_audits(fixture: &LedgerAuthorityFixture, scheduler: &mut VoiScheduler) {
    for (sequence, (record, (recommended, alternative))) in
        fixture.records.iter().zip(&fixture.outcomes).enumerate()
    {
        let receipt = record_prospective_audit_receipt(
            &fixture.ledger,
            &fixture.vcs,
            record,
            sequence,
            &fixture.context,
            fixture.decision_artifact,
            fixture.model_artifact,
            *recommended,
            *alternative,
        )
        .expect("record prospective audit receipt");
        let authenticated = verify_prospective_audit_receipt(
            &fixture.ledger,
            &fixture.vcs,
            record.clone(),
            sequence,
            fixture.context.clone(),
            fixture.decision_artifact,
            fixture.model_artifact,
            *recommended,
            *alternative,
            receipt,
            10,
        )
        .expect("independently verify prospective audit receipt");
        scheduler
            .observe_authenticated_audit(authenticated)
            .expect("append authenticated prospective audit");
    }
}

#[test]
fn voi_001_sensitivity_from_cached_sweeps() {
    let calls = Cell::new(0usize);
    let counting = |v: &[f64]| -> f64 {
        calls.set(calls.get() + 1);
        margin(v)
    };
    let decision = LiveDecision {
        margin: &counting,
        arity: 3,
    };
    let ns = nodes();
    let grid = 32;
    let p_pivotal = flip_probability(&decision, &ns, 0, grid).expect("valid pivotal sweep");
    let p_irrelevant = flip_probability(&decision, &ns, 1, grid).expect("valid irrelevant sweep");
    let p_mild = flip_probability(&decision, &ns, 2, grid).expect("valid mild sweep");
    println!(
        "{{\"metric\":\"sensitivity\",\"drag_gap\":{p_pivotal:.3},\
         \"mass_penalty\":{p_irrelevant:.3},\"hazard\":{p_mild:.3},\"surrogate_calls\":{}}}",
        calls.get()
    );
    assert!(
        p_pivotal > 0.3,
        "the straddling node is pivotal: {p_pivotal}"
    );
    assert_eq!(
        p_irrelevant.to_bits(),
        0.0f64.to_bits(),
        "the one-sided node cannot flip"
    );
    assert!(
        p_mild > 0.0 && p_mild < p_pivotal,
        "mildly pivotal: {p_mild}"
    );
    // NEAR-FREE: three sweeps cost 3·(grid + 1) surrogate calls — no
    // solver in the loop, no combinatorial sweep.
    assert!(
        calls.get() <= 3 * (grid + 1),
        "cached-sweep budget held: {} calls",
        calls.get()
    );
    verdict(
        "voi-001",
        "sampled flip fractions from cached interval sweeps: pivotal 0.4-class, one-sided \
         exactly 0, mild in between — at <= 3(grid+1) surrogate calls (near-free)",
    );
}

#[test]
fn voi_002_ranking_is_flip_prob_per_dollar() {
    let decision = LiveDecision {
        margin: &margin,
        arity: 3,
    };
    let ns = nodes();
    let menu = vec![
        // Cheap AND decisive: must rank first.
        Probe {
            name: "climb-rung-drag".to_string(),
            target: "drag-gap".to_string(),
            cost: 10.0,
            shrink: 0.2,
            kind: ProbeKind::Computational,
        },
        // Decisive but expensive: positive score, ranks below.
        Probe {
            name: "wind-tunnel-drag".to_string(),
            target: "drag-gap".to_string(),
            cost: 500.0,
            shrink: 0.05,
            kind: ProbeKind::Physical,
        },
        // Cheap but irrelevant: zero score, ranks last.
        Probe {
            name: "refine-mass-model".to_string(),
            target: "mass-penalty".to_string(),
            cost: 5.0,
            shrink: 0.1,
            kind: ProbeKind::Computational,
        },
    ];
    let ranked = rank_purchases(&decision, &ns, &menu, 64).expect("valid probe menu");
    let names: Vec<&str> = ranked.iter().map(|r| r.probe().name.as_str()).collect();
    println!(
        "{{\"metric\":\"ranking\",\"order\":{names:?},\"scores\":{:?}}}",
        ranked
            .iter()
            .map(|r| (r.score() * 1e4).round() / 1e4)
            .collect::<Vec<_>>()
    );
    assert_eq!(names[0], "climb-rung-drag", "cheap+decisive wins");
    assert_eq!(names[1], "wind-tunnel-drag", "decisive-but-pricey second");
    assert_eq!(names[2], "refine-mass-model", "irrelevant last");
    let first = ranked.get(0).expect("first ranked purchase");
    let second = ranked.get(1).expect("second ranked purchase");
    let third = ranked.get(2).expect("third ranked purchase");
    assert!(first.score() > second.score() && second.score() > 0.0);
    assert_eq!(
        third.score().to_bits(),
        0.0f64.to_bits(),
        "an irrelevant probe buys nothing"
    );
    verdict(
        "voi-002",
        "the ranking is sampled flip-fraction reduction per dollar: cheap+decisive > decisive+pricey > \
         irrelevant (score exactly 0)",
    );
}

#[test]
fn voi_003_menu_unifies_compute_and_physical() {
    // A physical anchor that shrinks the pivotal node far harder can
    // BEAT a cheap but weak computational probe — the menu prices
    // evidence, not its substrate.
    let decision = LiveDecision {
        margin: &margin,
        arity: 3,
    };
    let ns = nodes();
    let menu = vec![
        Probe {
            name: "weak-rung-climb".to_string(),
            target: "drag-gap".to_string(),
            cost: 8.0,
            shrink: 0.95, // barely tightens
            kind: ProbeKind::Computational,
        },
        Probe {
            name: "wind-tunnel-anchor".to_string(),
            target: "drag-gap".to_string(),
            cost: 60.0,
            shrink: 0.05, // near-eliminates the interval
            kind: ProbeKind::Physical,
        },
    ];
    let ranked = rank_purchases(&decision, &ns, &menu, 64).expect("valid mixed menu");
    let top = ranked.top().expect("nonempty sealed menu");
    let second = ranked.get(1).expect("second ranked purchase");
    assert_eq!(top.probe().kind, ProbeKind::Physical);
    assert!(
        top.score() > second.score(),
        "the physical anchor wins on flip-prob-per-dollar: {:.4} vs {:.4}",
        top.score(),
        second.score()
    );
    verdict(
        "voi-003",
        "a decisive wind-tunnel anchor outranks a weak rung climb despite 7.5x the \
         price — computational and physical evidence priced on one menu",
    );
}

#[test]
fn voi_004_surfacing_hint_and_scheduler() {
    let decision = LiveDecision {
        margin: &margin,
        arity: 3,
    };
    let ns = nodes();
    let menu = vec![
        Probe {
            name: "climb-rung-drag".to_string(),
            target: "drag-gap".to_string(),
            cost: 10.0,
            shrink: 0.2,
            kind: ProbeKind::Computational,
        },
        Probe {
            name: "hazard-samples".to_string(),
            target: "hazard".to_string(),
            cost: 25.0,
            shrink: 0.3,
            kind: ProbeKind::Computational,
        },
        Probe {
            name: "wind-tunnel-drag".to_string(),
            target: "drag-gap".to_string(),
            cost: 500.0,
            shrink: 0.05,
            kind: ProbeKind::Physical,
        },
    ];
    let ranked = rank_purchases(&decision, &ns, &menu, 64).expect("valid scheduling menu");
    // (i) The query-result hint (the Proposal-8 anytime shape, now
    // decision-priced).
    let hint = hint_for_query(&ranked);
    println!("{{\"metric\":\"hint\",\"hint\":{}}}", hint.to_json());
    let text = hint.render_text();
    assert!(text.contains("climb-rung-drag") && text.contains("$10"));
    assert!(text.contains("sampled flip fraction") && text.contains("64-point"));
    // (ii) Scheduling requires independently won audit authority and executes
    // one myopic purchase before the caller must update and rerank.
    // Source menu order carries no authority: sealing canonicalizes both rows
    // and the context identity.
    let mut reordered_menu = menu.clone();
    reordered_menu.reverse();
    let reordered =
        rank_purchases(&decision, &ns, &reordered_menu, 64).expect("reordered source menu");
    assert_eq!(reordered.context_id(), ranked.context_id());
    assert_eq!(reordered, ranked);
    let mut scheduler = winning_scheduler(40.0);
    let audit = scheduler.audit_report().expect("live audit snapshot");
    let scheduled = scheduler
        .schedule(ranked)
        .expect("valid finite schedule")
        .expect("one affordable purchase");
    assert_eq!(scheduled.purchase().probe().name, "climb-rung-drag");
    assert_eq!(scheduled.ranked_context_id(), reordered.context_id());
    assert_eq!(scheduled.ranked_grid(), 64);
    assert_eq!(scheduled.policy_scope(), POLICY_SCOPE);
    assert_eq!(scheduled.snapshot_id(), SNAPSHOT_ID);
    assert_eq!(scheduled.audit_context_id(), audit.audit_context_id());
    assert_eq!(scheduled.audit_observations(), audit.observations());
    assert_eq!(
        scheduled.audit_log_e_value().to_bits(),
        audit.log_e_value().to_bits()
    );
    assert_eq!(scheduled.budget_dollars().to_bits(), 40.0f64.to_bits());
    assert_eq!(
        scheduled.remaining_budget_dollars().to_bits(),
        30.0f64.to_bits()
    );
    verdict(
        "voi-004",
        "the structured hint is grid-qualified; a sealed context ignores source-menu order; \
         authority schedules exactly one purchase before mandatory reranking",
    );
}

#[test]
fn voi_005_prospective_audit_kill_criterion() {
    let empty = audit_scheduling(&[]).expect("empty audit reports safely");
    assert_eq!(empty.verdict(), AuditVerdict::DemoteToReporting);

    let insufficient = vec![audit_record(0, true)];
    let insufficient = audit_scheduling(&insufficient).expect("bounded prefix");
    assert_eq!(insufficient.verdict(), AuditVerdict::DemoteToReporting);

    let losing: Vec<_> = (0..128).map(|index| audit_record(index, false)).collect();
    let losing = audit_scheduling(&losing).expect("valid losing audit");
    assert_eq!(losing.verdict(), AuditVerdict::DemoteToReporting);

    let mut winning: Vec<_> = (0..128).map(|index| audit_record(index, true)).collect();
    let won = audit_scheduling(&winning).expect("valid winning audit");
    assert_eq!(won.verdict(), AuditVerdict::KeepScheduling);
    winning.reverse();
    let replay = audit_scheduling(&winning).expect("reversed chronology remains well formed");
    assert_ne!(replay.audit_context_id(), won.audit_context_id());
    assert_eq!(replay.log_e_value().to_bits(), won.log_e_value().to_bits());

    let mut wins_then_loss: Vec<_> = (0..11).map(|index| audit_record(index, true)).collect();
    wins_then_loss.push(audit_record(11, false));
    let mut loss_second = wins_then_loss.clone();
    let loss = loss_second.pop().expect("loss record");
    loss_second.insert(1, loss);
    let wins_first = audit_scheduling(&wins_then_loss).expect("prospective winning prefix");
    let loss_second = audit_scheduling(&loss_second).expect("different prospective order");
    assert_eq!(wins_first.verdict(), AuditVerdict::KeepScheduling);
    assert_eq!(loss_second.verdict(), AuditVerdict::DemoteToReporting);
    assert_ne!(
        wins_first.audit_context_id(),
        loss_second.audit_context_id()
    );
    verdict(
        "voi-005",
        "the fixed-alpha pairwise e-process mints authority only after a sufficient winning \
         matched-cost prefix; empty, short, and losing audits remain reporting-only",
    );
}

#[test]
fn voi_006_node_name_and_grid_boundaries() {
    let boundary_nodes: Vec<UncertaintyNode> = (0..MAX_VOI_NODES)
        .map(|index| UncertaintyNode {
            name: format!("node-{index}"),
            lo: 0.0,
            hi: 1.0,
            nominal: 0.5,
        })
        .collect();
    let wide_margin = |values: &[f64]| values[0] - 0.5;
    let boundary_decision = LiveDecision {
        margin: &wide_margin,
        arity: MAX_VOI_NODES,
    };
    assert!(
        nominal_verdict(&boundary_decision, &boundary_nodes).is_ok(),
        "the exact node-count boundary is admitted"
    );
    let mut too_many = boundary_nodes.clone();
    too_many.push(UncertaintyNode {
        name: "one-too-many".to_string(),
        lo: 0.0,
        hi: 1.0,
        nominal: 0.5,
    });
    let oversized_decision = LiveDecision {
        margin: &wide_margin,
        arity: too_many.len(),
    };
    assert!(matches!(
        nominal_verdict(&oversized_decision, &too_many),
        Err(VoiError::SizeLimit { .. })
    ));

    let exact_name = vec![UncertaintyNode {
        name: "n".repeat(MAX_VOI_NAME_BYTES),
        lo: 0.0,
        hi: 1.0,
        nominal: 0.5,
    }];
    let one = LiveDecision {
        margin: &wide_margin,
        arity: 1,
    };
    assert!(nominal_verdict(&one, &exact_name).is_ok());
    let mut long_name = exact_name;
    long_name[0].name.push('n');
    assert!(matches!(
        nominal_verdict(&one, &long_name),
        Err(VoiError::InvalidName { .. })
    ));

    let single = vec![UncertaintyNode {
        name: "x".to_string(),
        lo: -1.0,
        hi: 1.0,
        nominal: 0.0,
    }];
    assert!(
        flip_probability(&one, &single, 0, MAX_VOI_GRID).is_ok(),
        "the exact grid boundary is admitted"
    );
    for grid in [0, MAX_VOI_GRID + 1] {
        assert!(matches!(
            flip_probability(&one, &single, 0, grid),
            Err(VoiError::InvalidGrid { .. })
        ));
    }
    verdict(
        "voi-006",
        "node count, visible-ASCII name bytes, and sweep grid admit exact boundaries and refuse limit+1",
    );
}

#[test]
fn voi_007_menu_and_aggregate_evaluation_boundaries() {
    let calls = Cell::new(0usize);
    let counting = |values: &[f64]| {
        calls.set(calls.get() + 1);
        values[0]
    };
    let decision = LiveDecision {
        margin: &counting,
        arity: 1,
    };
    let ns = vec![UncertaintyNode {
        name: "x".to_string(),
        lo: -1.0,
        hi: 1.0,
        nominal: 0.0,
    }];
    let exact_evaluation_menu: Vec<Probe> = (0..818)
        .map(|index| probe(&format!("probe-{index}"), "x", index as f64 + 1.0, 0.5))
        .collect();
    assert_eq!(
        1 + 5 * (exact_evaluation_menu.len() + 1),
        MAX_VOI_EVALUATIONS
    );
    let ranked = rank_purchases(&decision, &ns, &exact_evaluation_menu, 5)
        .expect("exact evaluation boundary is admitted");
    assert_eq!(ranked.len(), exact_evaluation_menu.len());
    assert_eq!(calls.get(), MAX_VOI_EVALUATIONS);

    calls.set(0);
    let menu: Vec<Probe> = (0..MAX_VOI_PROBES)
        .map(|index| probe(&format!("max-probe-{index}"), "x", index as f64 + 1.0, 0.5))
        .collect();
    let ranked =
        rank_purchases(&decision, &ns, &menu, 3).expect("exact menu-count boundary is admitted");
    assert_eq!(ranked.len(), MAX_VOI_PROBES);
    assert_eq!(calls.get(), 1 + 3 * (MAX_VOI_PROBES + 1));

    calls.set(0);
    assert!(matches!(
        rank_purchases(&decision, &ns, &menu, 4),
        Err(VoiError::EvaluationLimitExceeded { .. })
    ));
    assert_eq!(calls.get(), 0, "aggregate refusal precedes callbacks");

    let mut oversized = menu;
    oversized.push(probe("one-too-many", "x", 1.0, 0.5));
    assert!(matches!(
        rank_purchases(&decision, &ns, &oversized, 1),
        Err(VoiError::SizeLimit { .. })
    ));
    assert_eq!(calls.get(), 0, "menu-size refusal precedes callbacks");
    verdict(
        "voi-007",
        "menu and aggregate evaluation caps admit exact limits and refuse before callback at limit+1",
    );
}

#[test]
fn voi_008_decision_interval_and_callback_refusals() {
    let calls = Cell::new(0usize);
    let counting = |values: &[f64]| {
        calls.set(calls.get() + 1);
        values.first().copied().unwrap_or(0.0)
    };
    let decision = LiveDecision {
        margin: &counting,
        arity: 1,
    };
    assert!(matches!(
        nominal_verdict(&decision, &[]),
        Err(VoiError::SizeLimit { .. })
    ));
    assert_eq!(calls.get(), 0);
    let ns = vec![UncertaintyNode {
        name: "x".to_string(),
        lo: 0.0,
        hi: 1.0,
        nominal: 0.5,
    }];
    let wrong_arity = LiveDecision {
        margin: &counting,
        arity: 2,
    };
    assert!(matches!(
        nominal_verdict(&wrong_arity, &ns),
        Err(VoiError::ArityMismatch { .. })
    ));
    assert!(matches!(
        flip_probability(&decision, &ns, 1, 4),
        Err(VoiError::NodeIndexOutOfRange { .. })
    ));

    let mut duplicate = ns.clone();
    duplicate.push(ns[0].clone());
    let duplicate_decision = LiveDecision {
        margin: &counting,
        arity: 2,
    };
    assert!(matches!(
        nominal_verdict(&duplicate_decision, &duplicate),
        Err(VoiError::DuplicateName { .. })
    ));
    for (lo, nominal, hi) in [
        (f64::NAN, 0.5, 1.0),
        (1.0, 0.5, 0.0),
        (0.0, 2.0, 1.0),
        (-f64::MAX, 0.0, f64::MAX),
    ] {
        let invalid = vec![UncertaintyNode {
            name: "x".to_string(),
            lo,
            hi,
            nominal,
        }];
        assert!(matches!(
            nominal_verdict(&decision, &invalid),
            Err(VoiError::InvalidInterval { .. })
        ));
    }

    let nan_margin = |_: &[f64]| f64::NAN;
    let malformed = LiveDecision {
        margin: &nan_margin,
        arity: 1,
    };
    assert!(matches!(
        nominal_verdict(&malformed, &ns),
        Err(VoiError::NonFiniteMargin { .. })
    ));
    let interior_nan = |values: &[f64]| {
        if values[0] > 0.5 { f64::INFINITY } else { 0.0 }
    };
    let malformed = LiveDecision {
        margin: &interior_nan,
        arity: 1,
    };
    assert!(matches!(
        flip_probability(&malformed, &ns, 0, 2),
        Err(VoiError::NonFiniteMargin { .. })
    ));
    verdict(
        "voi-008",
        "arity, identity, intervals, indices, duplicate nodes, and nonfinite callbacks refuse structurally",
    );
}

#[test]
fn voi_009_menu_targets_and_probe_economics_refuse_before_sweep() {
    let calls = Cell::new(0usize);
    let counting = |values: &[f64]| {
        calls.set(calls.get() + 1);
        values[0]
    };
    let decision = LiveDecision {
        margin: &counting,
        arity: 1,
    };
    let ns = vec![UncertaintyNode {
        name: "x".to_string(),
        lo: -1.0,
        hi: 1.0,
        nominal: 0.0,
    }];
    assert!(matches!(
        rank_purchases(&decision, &ns, &[], 4),
        Err(VoiError::SizeLimit { .. })
    ));
    let duplicate_nodes = vec![ns[0].clone(), ns[0].clone()];
    let duplicate_target_decision = LiveDecision {
        margin: &counting,
        arity: 2,
    };
    assert!(matches!(
        rank_purchases(
            &duplicate_target_decision,
            &duplicate_nodes,
            &[probe("ambiguous", "x", 1.0, 0.5)],
            4,
        ),
        Err(VoiError::DuplicateName { .. })
    ));
    let unknown = vec![probe("unknown", "missing", 1.0, 0.5)];
    let unknown_before = unknown.clone();
    assert!(matches!(
        rank_purchases(&decision, &ns, &unknown, 4),
        Err(VoiError::TargetResolution { matches: 0, .. })
    ));
    assert_eq!(unknown, unknown_before, "refusal does not mutate the menu");
    assert_eq!(calls.get(), 0, "unknown target refuses before sweep");

    let duplicated = vec![probe("same", "x", 1.0, 0.5), probe("same", "x", 2.0, 0.5)];
    assert!(matches!(
        rank_purchases(&decision, &ns, &duplicated, 4),
        Err(VoiError::DuplicateName { .. })
    ));
    for cost in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            rank_purchases(&decision, &ns, &[probe("bad-cost", "x", cost, 0.5)], 4),
            Err(VoiError::InvalidProbeValue { field: "cost", .. })
        ));
    }
    for shrink in [0.0, 1.0, -0.1, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            rank_purchases(&decision, &ns, &[probe("bad-shrink", "x", 1.0, shrink)], 4),
            Err(VoiError::InvalidProbeValue {
                field: "shrink",
                ..
            })
        ));
    }
    for bad_name in ["", " padded", "trailing "] {
        assert!(matches!(
            rank_purchases(&decision, &ns, &[probe(bad_name, "x", 1.0, 0.5)], 4),
            Err(VoiError::InvalidName { .. })
        ));
    }
    assert_eq!(calls.get(), 0, "all malformed menus refuse before sweep");
    verdict(
        "voi-009",
        "unknown/ambiguous targets, duplicate identities, and invalid probe economics refuse before evaluation",
    );
}

#[test]
fn voi_010_probe_name_and_score_arithmetic_boundaries() {
    let calls = Cell::new(0usize);
    let counting = |values: &[f64]| {
        calls.set(calls.get() + 1);
        values[0]
    };
    let decision = LiveDecision {
        margin: &counting,
        arity: 1,
    };
    let ns = vec![UncertaintyNode {
        name: "x".to_string(),
        lo: -1.0,
        hi: 1.0,
        nominal: 0.0,
    }];
    let exact = "p".repeat(MAX_VOI_NAME_BYTES);
    assert!(
        rank_purchases(&decision, &ns, &[probe(&exact, "x", 1.0, 0.5)], 4).is_ok(),
        "the exact probe-name byte boundary is admitted"
    );
    let before_long_name = calls.get();
    let too_long = "p".repeat(MAX_VOI_NAME_BYTES + 1);
    assert!(matches!(
        rank_purchases(&decision, &ns, &[probe(&too_long, "x", 1.0, 0.5)], 4),
        Err(VoiError::InvalidName { .. })
    ));
    assert_eq!(calls.get(), before_long_name);
    let too_long_target = "t".repeat(MAX_VOI_NAME_BYTES + 1);
    assert!(matches!(
        rank_purchases(
            &decision,
            &ns,
            &[probe("long-target", &too_long_target, 1.0, 0.5)],
            4,
        ),
        Err(VoiError::InvalidName { .. })
    ));

    let decision = LiveDecision {
        margin: &margin,
        arity: 3,
    };
    assert!(matches!(
        rank_purchases(
            &decision,
            &nodes(),
            &[probe("underflow-cost", "drag-gap", f64::from_bits(1), 0.01,)],
            64,
        ),
        Err(VoiError::ArithmeticRefusal { .. })
    ));
    verdict(
        "voi-010",
        "probe names admit the exact byte bound, limit+1 refuses before callback, and derived score overflow refuses",
    );
}

#[test]
fn voi_011_scheduler_is_transactional_and_budget_monotone() {
    let mut zero_budget = winning_scheduler(0.0);
    assert!(
        zero_budget
            .schedule(scheduler_menu())
            .expect("zero budget is valid")
            .is_none()
    );
    assert_eq!(zero_budget.consumed_snapshots(), 1);
    for budget in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            VoiScheduler::new(POLICY_SCOPE, budget),
            Err(VoiError::InvalidBudget { .. })
        ));
    }
    let mut unaudited = VoiScheduler::new(POLICY_SCOPE, 15.0).expect("valid scheduler");
    assert!(matches!(
        unaudited.schedule(scheduler_menu()),
        Err(VoiError::MissingSchedulingAuthority)
    ));
    let no_progress = ranked_menu(&[probe("tiny", "drag-gap", 1.0, 0.01)]);
    let mut huge_budget = winning_scheduler(f64::MAX);
    assert!(matches!(
        huge_budget.schedule(no_progress),
        Err(VoiError::ArithmeticRefusal { .. })
    ));
    assert_eq!(huge_budget.consumed_snapshots(), 0);

    let mut scheduler = winning_scheduler(15.0);
    let scheduled = scheduler
        .schedule(scheduler_menu())
        .expect("exact finite budget")
        .expect("one purchase");
    assert_eq!(scheduled.purchase().probe().name, "b");
    assert_eq!(
        scheduler.remaining_budget_dollars().to_bits(),
        10.0f64.to_bits()
    );
    assert!(matches!(
        scheduler.schedule(scheduler_menu()),
        Err(VoiError::RankingSnapshotAlreadyConsumed { .. })
    ));
    let second =
        ranked_menu_for_snapshot(&[probe("a", "drag-gap", 10.0, 0.01)], "fixture-snapshot-v2");
    let second = scheduler
        .schedule(second)
        .expect("fresh snapshot")
        .expect("remaining budget admits exact-cost purchase");
    assert_eq!(
        second.remaining_budget_dollars().to_bits(),
        0.0f64.to_bits()
    );
    assert_eq!(
        scheduler.remaining_budget_dollars().to_bits(),
        0.0f64.to_bits()
    );
    verdict(
        "voi-011",
        "one live scheduler gates audit authority, consumes each snapshot once, and owns a monotone cumulative budget",
    );
}

#[test]
fn voi_012_asymmetric_contraction_is_a_subset_and_preflights() {
    let samples = RefCell::new(Vec::new());
    let recording = |values: &[f64]| {
        samples.borrow_mut().push(values[0]);
        values[0] - 0.5
    };
    let decision = LiveDecision {
        margin: &recording,
        arity: 1,
    };
    let asymmetric = vec![UncertaintyNode {
        name: "x".to_string(),
        lo: 0.0,
        hi: 1.0,
        nominal: 0.9,
    }];
    rank_purchases(
        &decision,
        &asymmetric,
        &[probe("contract", "x", 1.0, 0.8)],
        2,
    )
    .expect("asymmetric interval contracts inside its support");
    let samples = samples.borrow();
    assert_eq!(samples.len(), 5, "one shared nominal plus two grid sweeps");
    let after = &samples[3..5];
    assert!(after.iter().all(|sample| (0.0..=1.0).contains(sample)));
    assert!((after[0] - 0.38).abs() < 1e-12, "left midpoint");
    assert!((after[1] - 0.78).abs() < 1e-12, "right midpoint");
    assert!((after[1] - after[0] - 0.4).abs() < 1e-12);

    let calls = Cell::new(0usize);
    let counting = |values: &[f64]| {
        calls.set(calls.get() + 1);
        values[0]
    };
    let decision = LiveDecision {
        margin: &counting,
        arity: 1,
    };
    let subnormal = vec![UncertaintyNode {
        name: "tiny".to_string(),
        lo: 0.0,
        hi: f64::from_bits(1),
        nominal: 0.0,
    }];
    assert!(matches!(
        rank_purchases(
            &decision,
            &subnormal,
            &[probe("underflow", "tiny", 1.0, 0.5)],
            2,
        ),
        Err(VoiError::ArithmeticRefusal { .. })
    ));
    assert_eq!(
        calls.get(),
        0,
        "derived intervals preflight before callbacks"
    );
}

#[test]
fn voi_013_context_and_sampled_zero_are_explicit() {
    let threshold = |values: &[f64]| values[0];
    let decision = LiveDecision {
        margin: &threshold,
        arity: 1,
    };
    let ns = vec![UncertaintyNode {
        name: "x".to_string(),
        lo: -1.0,
        hi: 1.0,
        nominal: 0.0,
    }];
    let menu = vec![probe("midpoint-alias", "x", 2.0, 0.5)];
    let grid_one = rank_purchases(&decision, &ns, &menu, 1).expect("grid one is estimated");
    let hint = hint_for_query(&grid_one);
    assert!(hint.purchase().is_none());
    assert!(hint.render_text().contains("no sampled purchase"));
    assert!(hint.render_text().contains("does not prove"));
    assert!(hint.to_json().contains("\"authoritative_zero\":false"));

    let grid_two = rank_purchases(&decision, &ns, &menu, 2).expect("second grid");
    assert_ne!(grid_one.context_id(), grid_two.context_id());
    let mut expanded_menu = menu.clone();
    expanded_menu.push(probe("second-probe", "x", 3.0, 0.25));
    let expanded =
        rank_purchases(&decision, &ns, &expanded_menu, 1).expect("expanded supplied menu");
    assert_ne!(grid_one.context_id(), expanded.context_id());
    let mut changed = ns.clone();
    changed[0].hi = 1.5;
    let changed = rank_purchases(&decision, &changed, &menu, 1).expect("changed snapshot");
    assert_ne!(grid_one.context_id(), changed.context_id());

    let shifted_margin = |values: &[f64]| values[0] + 0.75;
    let shifted = LiveDecision {
        margin: &shifted_margin,
        arity: 1,
    };
    let source_a = rank_purchases_scoped(&decision, &ns, &menu, 4, POLICY_SCOPE, SNAPSHOT_ID)
        .expect("first decision model");
    let source_b = rank_purchases_scoped(&shifted, &ns, &menu, 4, POLICY_SCOPE, SNAPSHOT_ID)
        .expect("second decision model");
    assert_eq!(source_a.source_context_id(), source_b.source_context_id());
    assert_ne!(
        source_a.context_id(),
        source_b.context_id(),
        "canonical ranking outputs participate in the final context identity"
    );

    let other_snapshot = rank_purchases_scoped(
        &decision,
        &ns,
        &menu,
        4,
        POLICY_SCOPE,
        "fixture-snapshot-other",
    )
    .expect("other snapshot");
    assert_ne!(
        source_a.source_context_id(),
        other_snapshot.source_context_id()
    );
    let other_policy =
        rank_purchases_scoped(&decision, &ns, &menu, 4, "fixture-policy-v2", SNAPSHOT_ID)
            .expect("other policy");
    assert_ne!(
        source_a.source_context_id(),
        other_policy.source_context_id()
    );
}

#[test]
fn voi_014_audit_validation_duplicates_and_work_bound() {
    assert!(matches!(
        MatchedAuditRecord::new("obs", "same", "same", "fixture/run", 1.0, 1.0, true, false,),
        Err(VoiError::InvalidAuditPair { .. })
    ));
    for (recommended, alternative) in [(1.0, 2.0), (f64::NAN, f64::NAN), (0.0, 0.0)] {
        assert!(matches!(
            MatchedAuditRecord::new(
                "obs",
                "recommended",
                "alternative",
                "fixture/run",
                recommended,
                alternative,
                true,
                false,
            ),
            Err(VoiError::InvalidAuditCost { .. })
        ));
    }
    assert!(matches!(
        MatchedAuditRecord::new(
            "bad\nobs",
            "recommended",
            "alternative",
            "fixture/run",
            1.0,
            1.0,
            true,
            false,
        ),
        Err(VoiError::InvalidName { .. })
    ));

    let duplicated = vec![audit_record(0, true), audit_record(0, true)];
    assert!(matches!(
        audit_scheduling(&duplicated),
        Err(VoiError::DuplicateAuditObservation { .. })
    ));
    let oversized: Vec<_> = (0..=MAX_VOI_AUDIT_RECORDS)
        .map(|index| audit_record(index, true))
        .collect();
    assert!(matches!(
        audit_scheduling(&oversized),
        Err(VoiError::SizeLimit { .. })
    ));
}

#[test]
fn voi_015_structured_hint_escapes_and_preserves_price() {
    let decision = LiveDecision {
        margin: &margin,
        arity: 3,
    };
    let quoted = probe("quoted\"probe", "drag-gap", 0.49, 0.01);
    let ranked = rank_purchases(&decision, &nodes(), &[quoted], 64).expect("visible ASCII quote");
    let hint = hint_for_query(&ranked);
    let text = hint.render_text();
    assert!(
        text.contains("$0.49"),
        "price is not rounded to zero: {text}"
    );
    assert!(
        text.contains("quoted\\\"probe"),
        "text escapes identity: {text}"
    );
    let json = hint.to_json();
    assert!(json.contains("quoted\\\"probe"));
    assert!(json.contains("\"cost_dollars\":0.49"));

    let calls = Cell::new(0usize);
    let counting = |values: &[f64]| {
        calls.set(calls.get() + 1);
        margin(values)
    };
    let decision = LiveDecision {
        margin: &counting,
        arity: 3,
    };
    assert!(matches!(
        rank_purchases(
            &decision,
            &nodes(),
            &[probe("line\nbreak", "drag-gap", 1.0, 0.5)],
            64,
        ),
        Err(VoiError::InvalidName { .. })
    ));
    assert_eq!(calls.get(), 0, "control characters refuse before callback");
}

#[test]
fn voi_016_live_scheduler_scopes_serializes_and_revokes_authority() {
    let decision = LiveDecision {
        margin: &margin,
        arity: 3,
    };
    let foreign = rank_purchases_scoped(
        &decision,
        &nodes(),
        &[probe("foreign", "drag-gap", 5.0, 0.01)],
        64,
        "foreign-policy-v1",
        "foreign-snapshot-v1",
    )
    .expect("valid foreign ranking");
    let mut scheduler = winning_scheduler(100.0);
    assert!(matches!(
        scheduler.schedule(foreign),
        Err(VoiError::PolicyScopeMismatch { .. })
    ));
    assert_eq!(scheduler.consumed_snapshots(), 0);
    assert_eq!(
        scheduler.remaining_budget_dollars().to_bits(),
        100.0f64.to_bits()
    );

    let mut demoted = false;
    for index in 128..MAX_VOI_AUDIT_RECORDS {
        scheduler
            .observe_audit(audit_record(index, false))
            .expect("append-only losing observation");
        if index % 16 == 15
            && scheduler.audit_report().expect("audit snapshot").verdict()
                == AuditVerdict::DemoteToReporting
        {
            demoted = true;
            break;
        }
    }
    assert!(
        demoted,
        "later prospective losses revoke live scheduling eligibility"
    );
    let after_demotion = ranked_menu_for_snapshot(
        &[probe("after-demotion", "drag-gap", 5.0, 0.01)],
        "fixture-snapshot-after-demotion",
    );
    assert!(matches!(
        scheduler.schedule(after_demotion),
        Err(VoiError::MissingSchedulingAuthority)
    ));

    let scheduler = std::sync::Arc::new(std::sync::Mutex::new(winning_scheduler(15.0)));
    let barrier = std::sync::Arc::new(std::sync::Barrier::new(3));
    let mut joins = Vec::new();
    for _ in 0..2 {
        let scheduler = std::sync::Arc::clone(&scheduler);
        let barrier = std::sync::Arc::clone(&barrier);
        let ranked = scheduler_menu();
        joins.push(std::thread::spawn(move || {
            barrier.wait();
            scheduler.lock().expect("scheduler lock").schedule(ranked)
        }));
    }
    barrier.wait();
    let outcomes = joins
        .into_iter()
        .map(|join| join.join().expect("scheduler thread"))
        .collect::<Vec<_>>();
    assert_eq!(
        outcomes.iter().filter(|outcome| outcome.is_ok()).count(),
        1,
        "exclusive mutation admits exactly one identical snapshot"
    );
    assert_eq!(
        outcomes
            .iter()
            .filter(|outcome| matches!(
                outcome,
                Err(VoiError::RankingSnapshotAlreadyConsumed { .. })
            ))
            .count(),
        1
    );
}

struct ControlledOracle<'a> {
    calls: &'a Cell<usize>,
    cancel_at: Option<usize>,
    refuse_at: Option<usize>,
    work_units_per_evaluation: u64,
}

impl DecisionOracle for ControlledOracle<'_> {
    fn arity(&self) -> usize {
        1
    }

    fn work_units_per_evaluation(&self) -> u64 {
        self.work_units_per_evaluation
    }

    fn evaluate(
        &self,
        cx: &Cx,
        permit: DecisionEvaluationPermit,
        values: &[f64],
    ) -> Result<f64, VoiError> {
        let call = self.calls.get() + 1;
        self.calls.set(call);
        assert_eq!(permit.ordinal() + 1, call);
        assert_eq!(permit.total_evaluations(), ORACLE_TEST_EVALUATIONS);
        assert_eq!(permit.charged_work_units(), self.work_units_per_evaluation);
        assert_eq!(
            permit.remaining_evaluations(),
            ORACLE_TEST_EVALUATIONS - call
        );
        assert!(
            permit.envelope().max_evaluations() >= permit.total_evaluations(),
            "the permit remains inside the caller envelope"
        );
        let completed_work = u64::try_from(call)
            .expect("fixture call count fits u64")
            .checked_mul(self.work_units_per_evaluation)
            .expect("fixture work accounting stays bounded");
        assert_eq!(
            permit.remaining_work_units(),
            permit.envelope().max_work_units() - completed_work
        );
        if self.cancel_at == Some(call) {
            cx.set_cancel_requested(true);
        }
        if self.refuse_at == Some(call) {
            return Err(VoiError::ArithmeticRefusal {
                operation: "fixture oracle refusal",
                subject: "x".to_string(),
            });
        }
        Ok(values[0] - 0.5)
    }
}

struct UnderdeclaringOracle<'a> {
    calls: &'a Cell<usize>,
}

impl DecisionOracle for UnderdeclaringOracle<'_> {
    fn arity(&self) -> usize {
        1
    }

    fn work_units_per_evaluation(&self) -> u64 {
        1
    }

    fn evaluate(
        &self,
        _cx: &Cx,
        permit: DecisionEvaluationPermit,
        _values: &[f64],
    ) -> Result<f64, VoiError> {
        self.calls.set(self.calls.get() + 1);
        Err(VoiError::WorkLimitExceeded {
            requested: 100,
            max: permit.charged_work_units(),
        })
    }
}

#[test]
#[allow(clippy::too_many_lines)] // one adversarial boundary/refusal matrix
fn voi_017_oracle_cancellation_and_work_budget_fail_closed() {
    assert_eq!(
        DecisionBudget::new(0, 1),
        Err(VoiError::InvalidEvaluationBudget {
            supplied: 0,
            max: MAX_VOI_EVALUATIONS,
        })
    );
    assert_eq!(
        DecisionBudget::new(1, 0),
        Err(VoiError::InvalidWorkBudget {
            supplied: 0,
            max: MAX_VOI_WORK_UNITS,
        })
    );
    assert_eq!(
        DecisionBudget::new(MAX_VOI_EVALUATIONS + 1, 1),
        Err(VoiError::InvalidEvaluationBudget {
            supplied: MAX_VOI_EVALUATIONS + 1,
            max: MAX_VOI_EVALUATIONS,
        })
    );
    assert_eq!(
        DecisionBudget::new(1, MAX_VOI_WORK_UNITS + 1),
        Err(VoiError::InvalidWorkBudget {
            supplied: MAX_VOI_WORK_UNITS + 1,
            max: MAX_VOI_WORK_UNITS,
        })
    );
    let nodes = [UncertaintyNode {
        name: "x".to_string(),
        lo: 0.0,
        hi: 1.0,
        nominal: 0.6,
    }];
    let menu = [probe("measure-x", "x", 1.0, 0.5)];
    let required_evaluations = ORACLE_TEST_EVALUATIONS;

    let pre_cancel_calls = Cell::new(0);
    let pre_cancel_oracle = ControlledOracle {
        calls: &pre_cancel_calls,
        cancel_at: None,
        refuse_at: None,
        work_units_per_evaluation: 1,
    };
    let cancelled = Cx::for_testing();
    cancelled.set_cancel_requested(true);
    let result = rank_purchases_impl(
        &cancelled,
        &pre_cancel_oracle,
        &nodes,
        &menu,
        2,
        DecisionBudget::new(required_evaluations, ORACLE_TEST_WORK_UNITS)
            .expect("exact computation budget"),
        POLICY_SCOPE,
        SNAPSHOT_ID,
    );
    assert_eq!(result, Err(VoiError::DecisionEvaluationCancelled));
    assert_eq!(
        pre_cancel_calls.get(),
        0,
        "pre-cancellation precedes every callback"
    );

    let exhausted_calls = Cell::new(0);
    let exhausted_oracle = ControlledOracle {
        calls: &exhausted_calls,
        cancel_at: None,
        refuse_at: None,
        work_units_per_evaluation: 1,
    };
    let exhausted = Cx::for_testing_with_budget(Budget::new().with_poll_quota(0));
    let result = rank_purchases_impl(
        &exhausted,
        &exhausted_oracle,
        &nodes,
        &menu,
        2,
        DecisionBudget::new(required_evaluations, ORACLE_TEST_WORK_UNITS)
            .expect("exact computation budget"),
        POLICY_SCOPE,
        SNAPSHOT_ID,
    );
    assert_eq!(result, Err(VoiError::DecisionEvaluationCancelled));
    assert_eq!(
        exhausted_calls.get(),
        0,
        "ambient asupersync budget refusal precedes every callback"
    );

    for cancel_at in 1..=required_evaluations {
        let calls = Cell::new(0);
        let oracle = ControlledOracle {
            calls: &calls,
            cancel_at: Some(cancel_at),
            refuse_at: None,
            work_units_per_evaluation: 1,
        };
        let result = rank_purchases_impl(
            &Cx::for_testing(),
            &oracle,
            &nodes,
            &menu,
            2,
            DecisionBudget::new(required_evaluations, ORACLE_TEST_WORK_UNITS)
                .expect("exact computation budget"),
            POLICY_SCOPE,
            SNAPSHOT_ID,
        );
        assert_eq!(result, Err(VoiError::DecisionEvaluationCancelled));
        assert_eq!(
            calls.get(),
            cancel_at,
            "cancel after oracle call {cancel_at}"
        );
    }

    let refusal_calls = Cell::new(0);
    let refusal_oracle = ControlledOracle {
        calls: &refusal_calls,
        cancel_at: None,
        refuse_at: Some(3),
        work_units_per_evaluation: 1,
    };
    let result = rank_purchases_impl(
        &Cx::for_testing(),
        &refusal_oracle,
        &nodes,
        &menu,
        2,
        DecisionBudget::new(required_evaluations, ORACLE_TEST_WORK_UNITS)
            .expect("exact computation budget"),
        POLICY_SCOPE,
        SNAPSHOT_ID,
    );
    assert!(matches!(
        result,
        Err(VoiError::ArithmeticRefusal {
            operation: "fixture oracle refusal",
            ..
        })
    ));
    assert_eq!(
        refusal_calls.get(),
        3,
        "oracle refusal returns no partial ranking"
    );

    let expensive_calls = Cell::new(0);
    let expensive_oracle = ControlledOracle {
        calls: &expensive_calls,
        cancel_at: None,
        refuse_at: None,
        work_units_per_evaluation: 100,
    };
    let result = rank_purchases_impl(
        &Cx::for_testing(),
        &expensive_oracle,
        &nodes,
        &menu,
        2,
        DecisionBudget::new(required_evaluations, 499).expect("bounded work budget"),
        POLICY_SCOPE,
        SNAPSHOT_ID,
    );
    assert_eq!(
        result,
        Err(VoiError::WorkLimitExceeded {
            requested: 500,
            max: 499,
        })
    );
    assert_eq!(
        expensive_calls.get(),
        0,
        "work refusal must precede every callback"
    );
    let expensive = rank_purchases_impl(
        &Cx::for_testing(),
        &expensive_oracle,
        &nodes,
        &menu,
        2,
        DecisionBudget::new(required_evaluations, 500).expect("exact expensive-work budget"),
        POLICY_SCOPE,
        SNAPSHOT_ID,
    )
    .expect("the exact declared-work envelope admits the expensive oracle");
    assert_eq!(expensive_calls.get(), required_evaluations);
    assert_eq!(expensive.computation().work_units(), 500);
    expensive_calls.set(0);

    let underdeclaring_calls = Cell::new(0);
    let underdeclaring = UnderdeclaringOracle {
        calls: &underdeclaring_calls,
    };
    let result = rank_purchases_impl(
        &Cx::for_testing(),
        &underdeclaring,
        &nodes,
        &menu,
        2,
        DecisionBudget::new(required_evaluations, ORACLE_TEST_WORK_UNITS)
            .expect("declared-work envelope"),
        POLICY_SCOPE,
        SNAPSHOT_ID,
    );
    assert_eq!(
        result,
        Err(VoiError::WorkLimitExceeded {
            requested: 100,
            max: 1,
        })
    );
    assert_eq!(
        underdeclaring_calls.get(),
        1,
        "a cooperative under-declaring oracle refuses before authority exists"
    );

    let invalid_calls = Cell::new(0);
    let invalid_oracle = ControlledOracle {
        calls: &invalid_calls,
        cancel_at: None,
        refuse_at: None,
        work_units_per_evaluation: 0,
    };
    let result = rank_purchases_impl(
        &Cx::for_testing(),
        &invalid_oracle,
        &nodes,
        &menu,
        2,
        DecisionBudget::new(required_evaluations, 500).expect("bounded work budget"),
        POLICY_SCOPE,
        SNAPSHOT_ID,
    );
    assert_eq!(
        result,
        Err(VoiError::InvalidOracleWorkUnits {
            work_units_per_evaluation: 0,
            max: MAX_VOI_WORK_UNITS,
        })
    );
    assert_eq!(
        invalid_calls.get(),
        0,
        "invalid oracle metadata must precede every callback"
    );

    let result = rank_purchases_impl(
        &Cx::for_testing(),
        &expensive_oracle,
        &nodes,
        &menu,
        2,
        DecisionBudget::new(required_evaluations - 1, 500).expect("evaluation budget"),
        POLICY_SCOPE,
        SNAPSHOT_ID,
    );
    assert_eq!(
        result,
        Err(VoiError::EvaluationLimitExceeded {
            requested: required_evaluations,
            max: required_evaluations - 1,
        })
    );
    assert_eq!(
        expensive_calls.get(),
        0,
        "evaluation refusal must precede every callback"
    );

    let calls = Cell::new(0);
    let oracle = ControlledOracle {
        calls: &calls,
        cancel_at: None,
        refuse_at: None,
        work_units_per_evaluation: 1,
    };
    let first = rank_purchases_impl(
        &Cx::for_testing(),
        &oracle,
        &nodes,
        &menu,
        2,
        DecisionBudget::new(5, 5).expect("exact budget"),
        POLICY_SCOPE,
        SNAPSHOT_ID,
    )
    .expect("exact budget admits the ranking");
    calls.set(0);
    let second = rank_purchases_impl(
        &Cx::for_testing(),
        &oracle,
        &nodes,
        &menu,
        2,
        DecisionBudget::new(6, 6).expect("larger budget"),
        POLICY_SCOPE,
        SNAPSHOT_ID,
    )
    .expect("larger budget admits the same numeric ranking");
    assert!(first.iter().eq(second.iter()));
    assert_ne!(first.source_context_id(), second.source_context_id());
    assert_ne!(first.context_id(), second.context_id());
    assert_eq!(first.computation().evaluations(), required_evaluations);
    assert_eq!(first.computation().work_units(), ORACLE_TEST_WORK_UNITS);
}

#[test]
fn voi_018_ledger_receipts_mint_authenticated_scheduling_authority() {
    let fixture = ledger_authority_fixture();
    let (ranked_receipt, ranked) = authenticate_ranked_menu(&fixture);
    let mut scheduler =
        VoiScheduler::new(POLICY_SCOPE, 15.0).expect("authenticated scheduler fixture");
    authenticate_audits(&fixture, &mut scheduler);
    assert_eq!(
        scheduler
            .audit_report()
            .expect("authenticated audit report")
            .verdict(),
        AuditVerdict::KeepScheduling
    );

    let authenticated = scheduler
        .schedule_authenticated(&fixture.ledger, &fixture.vcs, ranked, 10)
        .expect("all ledger authority independently reverifies")
        .expect("one authenticated affordable purchase");
    assert_eq!(authenticated.ranked_receipt(), ranked_receipt);
    assert_eq!(
        authenticated.audit_receipts().len(),
        AUTHENTICATED_AUDIT_WINS
    );
    assert_eq!(authenticated.valid_through_day(), 20);
    assert_eq!(
        authenticated.scheduled().purchase().probe().name,
        "authenticated-b"
    );
    assert_eq!(
        authenticated.scheduled().snapshot_id(),
        fixture.context.snapshot_root().to_hex()
    );
    assert_eq!(
        scheduler.remaining_budget_dollars().to_bits(),
        10.0f64.to_bits()
    );
    verdict(
        "voi-018",
        "sealed ranked-menu and prospective matched-outcome receipts bind exact ledger inputs, current snapshot, policy/version scope, expiry, and the authenticated budget transition",
    );
}

#[test]
fn voi_019_ledger_authority_refuses_forgery_expiry_outcome_swap_and_stale_head() {
    let mut fixture = ledger_authority_fixture();
    let ranked = ranked_menu_for_snapshot(&fixture.menu, &fixture.context.snapshot_root().to_hex());
    let ranked_receipt = record_ranked_menu_receipt(
        &fixture.ledger,
        &fixture.vcs,
        &ranked,
        &fixture.context,
        fixture.decision_artifact,
        fixture.model_artifact,
        fixture.catalog_artifact,
    )
    .expect("record ranked receipt");

    let forged = VoiLedgerReceipt::from_parts(
        ranked_receipt.operation(),
        hash_bytes(b"not-the-ranked-receipt"),
    );
    let forged_result = verify_ranked_menu_receipt(
        &fixture.ledger,
        &fixture.vcs,
        ranked_menu_for_snapshot(&fixture.menu, &fixture.context.snapshot_root().to_hex()),
        fixture.context.clone(),
        fixture.decision_artifact,
        fixture.model_artifact,
        fixture.catalog_artifact,
        forged,
        10,
    );
    assert!(matches!(
        forged_result,
        Err(VoiError::ReceiptMismatch {
            field: "receipt locator"
        })
    ));

    let wrong_model_context = VoiLedgerContext::new(
        "session-v1",
        "decision-v1",
        "model-v2",
        POLICY_SCOPE,
        MAIN_BRANCH,
        fixture.context.snapshot_root(),
        10,
        20,
    )
    .expect("well-formed but wrong model scope");
    let wrong_model_result = verify_ranked_menu_receipt(
        &fixture.ledger,
        &fixture.vcs,
        ranked_menu_for_snapshot(&fixture.menu, &fixture.context.snapshot_root().to_hex()),
        wrong_model_context,
        fixture.decision_artifact,
        fixture.model_artifact,
        fixture.catalog_artifact,
        ranked_receipt,
        10,
    );
    assert!(matches!(
        wrong_model_result,
        Err(VoiError::ReceiptMismatch {
            field: "artifact envelope"
        })
    ));

    let wrong_policy_context = VoiLedgerContext::new(
        "session-v1",
        "decision-v1",
        "model-v1",
        "fixture-policy-v2",
        MAIN_BRANCH,
        fixture.context.snapshot_root(),
        10,
        20,
    )
    .expect("well-formed but wrong policy scope");
    let wrong_policy_result = verify_ranked_menu_receipt(
        &fixture.ledger,
        &fixture.vcs,
        ranked_menu_for_snapshot(&fixture.menu, &fixture.context.snapshot_root().to_hex()),
        wrong_policy_context,
        fixture.decision_artifact,
        fixture.model_artifact,
        fixture.catalog_artifact,
        ranked_receipt,
        10,
    );
    assert!(matches!(
        wrong_policy_result,
        Err(VoiError::ReceiptMismatch {
            field: "ranked policy version"
        })
    ));

    let expired_result = verify_ranked_menu_receipt(
        &fixture.ledger,
        &fixture.vcs,
        ranked_menu_for_snapshot(&fixture.menu, &fixture.context.snapshot_root().to_hex()),
        fixture.context.clone(),
        fixture.decision_artifact,
        fixture.model_artifact,
        fixture.catalog_artifact,
        ranked_receipt,
        21,
    );
    assert!(matches!(
        expired_result,
        Err(VoiError::ReceiptExpired {
            expires_day: 20,
            current_day: 21
        })
    ));

    let record = fixture.records[0].clone();
    let (recommended, alternative) = fixture.outcomes[0];
    let audit_receipt = record_prospective_audit_receipt(
        &fixture.ledger,
        &fixture.vcs,
        &record,
        0,
        &fixture.context,
        fixture.decision_artifact,
        fixture.model_artifact,
        recommended,
        alternative,
    )
    .expect("record audit receipt");
    let swapped_result = verify_prospective_audit_receipt(
        &fixture.ledger,
        &fixture.vcs,
        record,
        0,
        fixture.context.clone(),
        fixture.decision_artifact,
        fixture.model_artifact,
        alternative,
        recommended,
        audit_receipt,
        10,
    );
    assert!(matches!(
        swapped_result,
        Err(VoiError::ReceiptMismatch {
            field: "artifact bytes"
        })
    ));

    let advanced = fixture
        .vcs
        .commit(&fixture.ledger, MAIN_BRANCH)
        .expect("advance ledger head after receipt issuance");
    assert_ne!(advanced.root, fixture.context.snapshot_root());
    let stale_result = verify_ranked_menu_receipt(
        &fixture.ledger,
        &fixture.vcs,
        ranked_menu_for_snapshot(&fixture.menu, &fixture.context.snapshot_root().to_hex()),
        fixture.context.clone(),
        fixture.decision_artifact,
        fixture.model_artifact,
        fixture.catalog_artifact,
        ranked_receipt,
        10,
    );
    assert!(matches!(
        stale_result,
        Err(VoiError::StaleLedgerSnapshot { .. })
    ));
    verdict(
        "voi-019",
        "forged locators, model/policy drift, expired contexts, swapped matched outcomes, and superseded ranked roots all fail closed before scheduler mutation",
    );
}

#[test]
fn voi_identity_versions_fail_closed() {
    let ranked = scheduler_menu();
    assert_eq!(
        ranked.source_identity_version(),
        VOI_RANKED_SOURCE_IDENTITY_VERSION
    );
    assert_eq!(ranked.identity_version(), VOI_RANKED_MENU_IDENTITY_VERSION);
    assert_eq!(
        ranked.admit_retained_identity_versions(
            VOI_RANKED_SOURCE_IDENTITY_VERSION,
            VOI_RANKED_MENU_IDENTITY_VERSION,
        ),
        Ok(())
    );

    for declared in [
        VOI_RANKED_SOURCE_IDENTITY_VERSION - 1,
        VOI_RANKED_SOURCE_IDENTITY_VERSION + 1,
    ] {
        assert_eq!(
            ranked.admit_retained_identity_versions(declared, VOI_RANKED_MENU_IDENTITY_VERSION,),
            Err(VoiError::UnsupportedIdentityVersion {
                identity: "VoI ranked source",
                declared,
                supported: VOI_RANKED_SOURCE_IDENTITY_VERSION,
            })
        );
    }
    for declared in [
        VOI_RANKED_MENU_IDENTITY_VERSION - 1,
        VOI_RANKED_MENU_IDENTITY_VERSION + 1,
    ] {
        assert_eq!(
            ranked.admit_retained_identity_versions(VOI_RANKED_SOURCE_IDENTITY_VERSION, declared,),
            Err(VoiError::UnsupportedIdentityVersion {
                identity: "VoI ranked menu",
                declared,
                supported: VOI_RANKED_MENU_IDENTITY_VERSION,
            })
        );
    }

    let report = audit_scheduling(&[audit_record(0, true)]).expect("valid retained audit root");
    assert_eq!(
        report.identity_version(),
        VOI_AUDIT_CONTEXT_IDENTITY_VERSION
    );
    assert_eq!(
        report.admit_retained_identity_version(VOI_AUDIT_CONTEXT_IDENTITY_VERSION),
        Ok(())
    );
    for declared in [
        VOI_AUDIT_CONTEXT_IDENTITY_VERSION - 1,
        VOI_AUDIT_CONTEXT_IDENTITY_VERSION + 1,
    ] {
        assert_eq!(
            report.admit_retained_identity_version(declared),
            Err(VoiError::UnsupportedIdentityVersion {
                identity: "VoI audit context",
                declared,
                supported: VOI_AUDIT_CONTEXT_IDENTITY_VERSION,
            })
        );
    }
}
