//! One-bet-per-lane cross-crate E2E (bead frankensim-ext-epic-gov-rjoq.6,
//! slice 2): the REAL no-mock composition the bead demands — two
//! independent theorem lanes plus competing preregistered mechanisms in
//! one comparison lane, driven through fs-govern admission, persisted
//! into a real FrankenSQLite-backed fs-ledger (events + content-addressed
//! artifacts), packaged as color-typed fs-package claims, re-checked by
//! the solver-free fs-checker, and REPLAYED: a fresh ledger fed the same
//! request sequence reproduces the retained decision log byte-for-byte,
//! and byte-identical retries stay idempotent. Terminal release happens
//! only against a finalization receipt whose artifact ACTUALLY exists in
//! the design ledger.

use fs_checker::check_against_root;
use fs_govern::{
    FinalizationReceipt, HeadToHeadCharter, IdempotencyKey, LaneCharter, LaneError, MechanismId,
    PortfolioLedger, PortfolioPolicy, ResourceEnvelope, TerminalKind,
};
use fs_ledger::{EventRow, Ledger};
use fs_package::{Claim, EvidencePackage, Provenance};

fn charter(statement: &str, falsifier: &str, class: &str) -> LaneCharter {
    LaneCharter::new(
        statement,
        "small-strain elasticity and quasi-static EM, polyhedral domains",
        &["isotropic material", "linear constitutive law"],
        "verified",
        "hand-checked FEEC reference",
        falsifier,
        class,
    )
    .expect("valid charter")
}

fn envelope(work: u64) -> ResourceEnvelope {
    ResourceEnvelope {
        work_units: work,
        memory_bytes: work * 1024,
        reviewer_slots: 1,
        falsification_capacity: 1,
    }
}

fn policy() -> PortfolioPolicy {
    PortfolioPolicy {
        global: ResourceEnvelope {
            work_units: 10_000,
            memory_bytes: 10_000 * 1024,
            reviewer_slots: 16,
            falsification_capacity: 16,
        },
        max_active_mechanisms: 8,
    }
}

/// One admission-side request, retained so the SAME sequence (same
/// idempotency keys, same payloads) can replay against a fresh
/// portfolio ledger.
enum Step {
    Preregister(HeadToHeadCharter, IdempotencyKey),
    Admit(LaneCharter, MechanismId, ResourceEnvelope, IdempotencyKey),
    Finalize(FinalizationReceipt, IdempotencyKey),
}

fn apply(ledger: &mut PortfolioLedger, step: &Step) -> Result<(), LaneError> {
    match step {
        Step::Preregister(charter, key) => ledger.preregister_comparison(charter.clone(), *key),
        Step::Admit(charter, mechanism, reservation, key) => {
            ledger.admit(charter, *mechanism, *reservation, *key)
        }
        Step::Finalize(receipt, key) => ledger.finalize(receipt, *key),
    }
}

#[test]
#[allow(clippy::too_many_lines)] // one auditable script IS the composition evidence
fn lanes_e2e_ledger_package_checker_replay() {
    let dir = std::env::temp_dir().join(format!("lanes-e2e-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("mkdir");
    let design_ledger =
        Ledger::open(dir.join("lanes.led").to_str().expect("utf8")).expect("design ledger");

    // ---- Lanes: two independent theorem lanes + one comparison lane.
    let elastic = charter(
        "graded CutFEM elasticity adapts at order p+1",
        "manufactured-solution refutation family",
        "elasticity-adaptivity",
    );
    let em = charter(
        "cohomology-preserving EM topology optimization keeps Betti witnesses",
        "adversarial-seam refutation family",
        "em-topology",
    );
    let gating = charter(
        "conformal surrogate gating preserves campaign coverage",
        "held-out coverage refutation family",
        "surrogate-gating",
    );
    let m_elastic = elastic
        .mechanism_id("equilibrated-flux adapt", 1)
        .expect("id");
    let m_em = em.mechanism_id("betti-witness opt", 1).expect("id");
    let c_split = gating.mechanism_id("split conformal", 1).expect("id");
    let c_jackknife = gating.mechanism_id("jackknife conformal", 1).expect("id");

    // The preregistration protocol is a REAL ledgered artifact.
    let prereg_bytes = b"h2h protocol: split vs jackknife conformal gating, shared budget 200";
    let prereg = design_ledger
        .put_artifact("proof-lane/preregistration", prereg_bytes, None)
        .expect("preregistration artifact");
    // The shared comparison budget must cover BOTH candidates on every
    // axis (each reserves one reviewer slot and one falsification run).
    let shared_budget = ResourceEnvelope {
        work_units: 200,
        memory_bytes: 200 * 1024,
        reviewer_slots: 2,
        falsification_capacity: 2,
    };
    let h2h = HeadToHeadCharter::new(&gating, &[c_split, c_jackknife], shared_budget, prereg.hash)
        .expect("comparison charter");

    // ---- The scripted request sequence (adversaries included).
    let rival = elastic.mechanism_id("rival adapt", 1).expect("id");
    let cosmetic_split = charter(
        "graded CutFEM elasticity adapts at order p+1, rephrased",
        "manufactured-solution refutation family",
        "elasticity-adaptivity",
    );
    let m_cosmetic = cosmetic_split.mechanism_id("cosmetic mech", 1).expect("id");
    let intruder = gating.mechanism_id("undeclared intruder", 1).expect("id");

    let mut script: Vec<Step> = vec![
        Step::Admit(
            elastic.clone(),
            m_elastic,
            envelope(100),
            IdempotencyKey::derive("admit-elastic"),
        ),
        Step::Admit(
            em.clone(),
            m_em,
            envelope(100),
            IdempotencyKey::derive("admit-em"),
        ),
        // Same-lane adversary: refused, recorded.
        Step::Admit(
            elastic.clone(),
            rival,
            envelope(10),
            IdempotencyKey::derive("admit-rival"),
        ),
        // Independence-class adversary on a different lane: refused.
        Step::Admit(
            cosmetic_split.clone(),
            m_cosmetic,
            envelope(10),
            IdempotencyKey::derive("admit-cosmetic"),
        ),
        Step::Preregister(h2h.clone(), IdempotencyKey::derive("prereg-gating")),
        Step::Admit(
            gating.clone(),
            c_split,
            envelope(80),
            IdempotencyKey::derive("admit-split"),
        ),
        Step::Admit(
            gating.clone(),
            c_jackknife,
            envelope(80),
            IdempotencyKey::derive("admit-jackknife"),
        ),
        // Undeclared comparison candidate: refused.
        Step::Admit(
            gating.clone(),
            intruder,
            envelope(1),
            IdempotencyKey::derive("admit-intruder"),
        ),
    ];

    // The refutation evidence for the elastic mechanism is itself a
    // REAL ledgered artifact; its content address seals the receipt.
    let refutation_bytes =
        b"manufactured-solution refutation: observed order p on the graded corner family";
    let refutation = design_ledger
        .put_artifact("proof-lane/refutation", refutation_bytes, None)
        .expect("refutation artifact");
    script.push(Step::Finalize(
        FinalizationReceipt::new(m_elastic, TerminalKind::Refuted, None, refutation.hash)
            .expect("sealed receipt"),
        IdempotencyKey::derive("finalize-elastic"),
    ));

    // ---- Run 1: execute the script against the portfolio ledger,
    // persisting every decision row as a structured ledger event.
    let mut portfolio = PortfolioLedger::new(policy());
    let mut expected: Vec<Result<(), LaneError>> = Vec::new();
    for (t, step) in script.iter().enumerate() {
        let outcome = apply(&mut portfolio, step);
        let row = portfolio
            .decisions()
            .last()
            .expect("every request records a decision")
            .to_json();
        design_ledger
            .append_event(&EventRow {
                session: None,
                t: (t + 1) as i64,
                kind: "proof-lane",
                payload: Some(&row),
            })
            .expect("decision event");
        expected.push(outcome);
    }
    assert!(
        expected[0].is_ok() && expected[1].is_ok(),
        "independent lanes admit"
    );
    assert!(
        matches!(expected[2], Err(LaneError::LaneOccupied { .. })),
        "same-lane rival refuses"
    );
    assert!(
        matches!(
            expected[3],
            Err(LaneError::IndependenceClassOccupied { .. })
        ),
        "declared-class split refuses"
    );
    assert!(expected[4].is_ok() && expected[5].is_ok() && expected[6].is_ok());
    assert!(
        matches!(expected[7], Err(LaneError::NotADeclaredCandidate { .. })),
        "undeclared candidate refuses"
    );
    assert!(
        expected[8].is_ok(),
        "refutation releases against the ledgered receipt"
    );
    assert_eq!(
        portfolio.active_count(),
        3,
        "em + two comparison candidates remain"
    );

    // Persist the COMPLETE decision log as the run's content-addressed
    // artifact of record.
    let log_json = portfolio.decisions_json(usize::MAX);
    let log_artifact = design_ledger
        .put_artifact("proof-lane/decision-log", log_json.as_bytes(), None)
        .expect("decision-log artifact");
    assert_eq!(
        design_ledger
            .get_artifact(&log_artifact.hash)
            .expect("artifact read")
            .expect("artifact present"),
        log_json.as_bytes(),
        "the design ledger round-trips the decision log bytes"
    );

    // ---- Package + solver-free checker: the portfolio outcome becomes
    // an independently checkable claim bundle bound to the ledger
    // artifact's content address.
    let pkg = EvidencePackage::new(Provenance::new(
        "frankensim@rjoq.6-e2e",
        "constellation-lock-fixture",
    ))
    .with_claim(Claim::estimated(
        "lane-elasticity",
        format!(
            "elasticity lane terminal: refuted against ledger artifact {}",
            refutation.hash.to_hex()
        ),
        "fs-govern/lanes/v2",
        0.0,
    ))
    .with_claim(Claim::estimated(
        "lane-em",
        "em lane carries one active unproven mechanism",
        "fs-govern/lanes/v2",
        1.0,
    ))
    .with_claim(Claim::estimated(
        "lane-gating",
        format!(
            "gating lane runs a preregistered 2-candidate comparison; decision log {}",
            log_artifact.hash.to_hex()
        ),
        "fs-govern/lanes/v2",
        2.0,
    ));
    let root = pkg.try_merkle_root().expect("bounded package");
    let report = check_against_root(&pkg, root);
    assert!(
        report.passed(),
        "solver-free re-check passes (root binding folds into findings): {report:?}"
    );
    // An adversarial root must NOT pass: the checker really binds it.
    let wrong = check_against_root(&pkg, fs_checker::ContentHash([0x5A; 32]));
    assert!(!wrong.passed(), "a mismatched expected root must refuse");

    // ---- Replay: a FRESH portfolio ledger fed the identical request
    // sequence reproduces the decision log BYTE-FOR-BYTE (G5), and a
    // second pass of byte-identical retries is fully idempotent (G4
    // crash/retry) — same rows, same bytes, no double-charge.
    let mut replayed = PortfolioLedger::new(policy());
    for step in &script {
        let _ = apply(&mut replayed, step);
    }
    assert_eq!(
        replayed.decisions_json(usize::MAX).as_bytes(),
        log_json.as_bytes(),
        "fresh-ledger replay reproduces the retained artifact exactly"
    );
    assert_eq!(replayed.active_count(), portfolio.active_count());
    assert_eq!(replayed.reserved(), portfolio.reserved());
    let rows_before = replayed.decisions().len();
    for step in &script {
        let _ = apply(&mut replayed, step);
    }
    assert_eq!(
        replayed.decisions().len(),
        rows_before,
        "byte-identical retries replay without new rows"
    );
    assert_eq!(
        replayed.decisions_json(usize::MAX).as_bytes(),
        log_json.as_bytes(),
        "idempotent retries leave the log bytes unchanged"
    );
    assert_eq!(
        replayed.reserved(),
        portfolio.reserved(),
        "no double-charge"
    );
}

/// G4 persistence-boundary fault drill (the bead's outstanding storage
/// lane): the process "crashes" mid-sequence — the in-memory portfolio
/// ledger is dropped and the design ledger handle is closed — then the
/// design ledger REOPENS on the same path. Every decision artifact
/// persisted before the fault must survive byte-for-byte, re-persisting
/// a recovered row must dedupe against the surviving artifact (no
/// double write), and re-executing the FULL request sequence with the
/// SAME idempotency keys must converge to a final decision log
/// byte-identical to a never-crashed control run — no divergence, no
/// duplicate effects, no partial admission across the fault.
#[test]
fn lanes_e2e_crash_reopen_recovery() {
    let dir = std::env::temp_dir().join(format!("lanes-crash-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("mkdir");
    let ledger_path = dir.join("crash.led");
    let ledger_path = ledger_path.to_str().expect("utf8");

    let lane_a = charter(
        "claim A survives restarts",
        "restart refutation family",
        "crash-a",
    );
    let lane_b = charter(
        "claim B survives restarts",
        "restart refutation family b",
        "crash-b",
    );
    let m_a = lane_a.mechanism_id("mech a", 1).expect("id");
    let m_b = lane_b.mechanism_id("mech b", 1).expect("id");
    let rival = lane_a.mechanism_id("rival", 1).expect("id");
    let script: Vec<Step> = vec![
        Step::Admit(
            lane_a.clone(),
            m_a,
            envelope(10),
            IdempotencyKey::derive("c-a"),
        ),
        Step::Admit(
            lane_a.clone(),
            rival,
            envelope(5),
            IdempotencyKey::derive("c-rival"),
        ),
        Step::Admit(
            lane_b.clone(),
            m_b,
            envelope(10),
            IdempotencyKey::derive("c-b"),
        ),
    ];

    // Control: the never-crashed run.
    let mut control = PortfolioLedger::new(policy());
    for step in &script {
        let _ = apply(&mut control, step);
    }
    let control_log = control.decisions_json(usize::MAX);

    // Run 1: crash after the FIRST step's decision is persisted.
    let mut persisted_before_crash = Vec::new();
    {
        let ledger = Ledger::open(ledger_path).expect("design ledger");
        let mut portfolio = PortfolioLedger::new(policy());
        let _ = apply(&mut portfolio, &script[0]);
        let row = portfolio.decisions().last().expect("row").to_json();
        let receipt = ledger
            .put_artifact("proof-lane/decision", row.as_bytes(), None)
            .expect("persist decision");
        assert!(!receipt.deduped, "first persistence is a fresh write");
        persisted_before_crash.push((receipt.hash, row));
        // CRASH: both handles drop here; in-memory portfolio state is
        // lost, the design-ledger file survives on disk.
    }

    // Recovery: reopen the SAME path; pre-fault artifacts must survive.
    let reopened = Ledger::open(ledger_path).expect("reopen after fault");
    for (hash, row) in &persisted_before_crash {
        let survived = reopened
            .get_artifact(hash)
            .expect("artifact read")
            .expect("artifact survived the fault");
        assert_eq!(survived, row.as_bytes(), "durable byte-for-byte survival");
    }

    // Re-execute the FULL sequence on a fresh portfolio (deterministic
    // replay makes the pre-fault prefix reproduce identically), then
    // persist every row: the recovered first row must DEDUPE against
    // the surviving pre-fault artifact — no double write.
    let mut recovered = PortfolioLedger::new(policy());
    for step in &script {
        let _ = apply(&mut recovered, step);
    }
    for (index, decision) in recovered.decisions().iter().enumerate() {
        let row = decision.to_json();
        let receipt = reopened
            .put_artifact("proof-lane/decision", row.as_bytes(), None)
            .expect("persist recovered decision");
        if index == 0 {
            assert!(receipt.deduped, "recovered prefix dedupes, not re-writes");
            assert_eq!(receipt.hash, persisted_before_crash[0].0);
        }
    }
    assert_eq!(
        recovered.decisions_json(usize::MAX),
        control_log,
        "post-fault recovery converges to the never-crashed control bytes"
    );

    // Idempotent retries after recovery stay absorbed.
    let rows = recovered.decisions().len();
    for step in &script {
        let _ = apply(&mut recovered, step);
    }
    assert_eq!(recovered.decisions().len(), rows, "no duplicate effects");
}
