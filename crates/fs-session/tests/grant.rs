//! Session-grant authority battery (bead aeq7, increment 1): the G0
//! drills from the bead acceptance — deny-all default, forged/altered,
//! expired, revoked, cross-issuer, ungranted-verb, wildcard-confusion,
//! concurrency-lease, and exact round-trip cases, all failing closed
//! with named typed errors.

use fs_session::{
    CapabilityToken, CoreLeaseBook, IssuerIdentity, IssuerPolicy, NoIssuerPolicy, PolicyDecision,
    SessionError, SessionGrant, SessionId, mint_grant,
};

fn request() -> CapabilityToken {
    CapabilityToken {
        session: SessionId(41),
        ops: vec!["flux.*".to_string(), "ascent.optimize".to_string()],
        core_s: 3600.0,
        mem_bytes: 64 * 1024 * 1024 * 1024,
        wall_s: 7200.0,
        cores: 8,
        ledger_scope: "studies/spout-v3".to_string(),
    }
}

/// A test policy with adjustable expiry and revocation generation.
struct TestPolicy {
    identity: IssuerIdentity,
    expiry_ns: i64,
    generation: std::sync::atomic::AtomicU64,
}

impl TestPolicy {
    fn new(expiry_ns: i64) -> TestPolicy {
        TestPolicy {
            identity: IssuerIdentity::new("ops/test-issuer", "policy-v1").expect("valid identity"),
            expiry_ns,
            generation: std::sync::atomic::AtomicU64::new(1),
        }
    }

    fn revoke(&self) {
        self.generation
            .fetch_add(1, std::sync::atomic::Ordering::SeqCst);
    }
}

impl IssuerPolicy for TestPolicy {
    fn issuer(&self) -> &IssuerIdentity {
        &self.identity
    }

    fn revocation_generation(&self) -> u64 {
        self.generation.load(std::sync::atomic::Ordering::SeqCst)
    }

    fn evaluate(&self, _request: &CapabilityToken, _issuance_ns: i64) -> PolicyDecision {
        PolicyDecision::Granted {
            expiry_ns: self.expiry_ns,
        }
    }
}

fn granted() -> (TestPolicy, SessionGrant) {
    let policy = TestPolicy::new(10_000);
    let grant = mint_grant(&policy, &request(), 1_000).expect("mints");
    (policy, grant)
}

#[test]
fn sg_001_deny_all_default_and_structural_refusals() {
    // Deny-all default: public construction cannot become authority.
    let deny = NoIssuerPolicy::new();
    let refused = mint_grant(&deny, &request(), 1_000);
    assert!(matches!(refused, Err(SessionError::GrantDenied { .. })));

    // Structural refusals fire before the policy ever runs.
    let policy = TestPolicy::new(10_000);
    let mut duplicate = request();
    duplicate.ops.push("flux.*".to_string());
    assert!(matches!(
        mint_grant(&policy, &duplicate, 1_000),
        Err(SessionError::DuplicateOperatorGrant { .. })
    ));
    let mut bad_scope = request();
    bad_scope.ledger_scope = "has whitespace".to_string();
    assert!(matches!(
        mint_grant(&policy, &bad_scope, 1_000),
        Err(SessionError::InvalidLedgerScope { .. })
    ));
    let mut bad_budget = request();
    bad_budget.wall_s = f64::NAN;
    assert!(matches!(
        mint_grant(&policy, &bad_budget, 1_000),
        Err(SessionError::InvalidResource { .. })
    ));
    // Expiry at/before issuance is unrepresentable authority.
    let backwards = TestPolicy::new(1_000);
    assert!(matches!(
        mint_grant(&backwards, &request(), 1_000),
        Err(SessionError::InvalidResource { .. })
    ));
    // Issuer identity fields are bounded canonical ASCII.
    assert!(matches!(
        IssuerIdentity::new("", "p"),
        Err(SessionError::InvalidIssuerField { .. })
    ));
    assert!(matches!(
        IssuerIdentity::new("ops/x", "bad fingerprint"),
        Err(SessionError::InvalidIssuerField { .. })
    ));
}

#[test]
fn sg_002_round_trip_and_admitted_view() {
    let (policy, grant) = granted();
    grant.verify_fresh(&policy, 5_000).expect("fresh grant");
    assert_eq!(grant.session(), SessionId(41));
    assert_eq!(grant.cores(), 8);
    assert!(!grant.digest().is_empty());
    // Admitted view mirrors the ADMITTED (sorted) operator set.
    let admission = grant.to_admission();
    assert_eq!(
        admission.ops,
        vec!["ascent.optimize".to_string(), "flux.*".to_string()],
        "ops are canonically sorted at mint"
    );
    assert_eq!(admission.cores, 8);
    // Verb coverage: exact, namespace, and the confusion cases.
    assert!(grant.grants_op("ascent.optimize"));
    assert!(grant.grants_op("flux.free-surface-lbm"));
    assert!(!grant.grants_op("flux"), "a namespace is not an operator");
    assert!(!grant.grants_op("fluxx.solve"), "prefix confusion refused");
    assert!(!grant.grants_op("ascent.solve-lp"), "exact means exact");
    // Determinism: identical mint inputs give identical digests.
    let again = mint_grant(&policy, &request(), 1_000).expect("mints again");
    assert_eq!(grant.digest(), again.digest());
}

#[test]
fn sg_003_expiry_revocation_and_cross_issuer_fail_closed() {
    let (policy, grant) = granted();
    // Expired: the admitted window is exclusive at expiry.
    assert!(matches!(
        grant.verify_fresh(&policy, 10_000),
        Err(SessionError::GrantExpired { .. })
    ));
    // Revocation: generation advance invalidates without touching the
    // grant bytes.
    grant.verify_fresh(&policy, 5_000).expect("still fresh");
    policy.revoke();
    assert!(matches!(
        grant.verify_fresh(&policy, 5_000),
        Err(SessionError::GrantRevoked { .. })
    ));
    // Cross-issuer: a different issuer (or rotated fingerprint) cannot
    // vouch for this grant.
    let other = TestPolicy {
        identity: IssuerIdentity::new("ops/other-issuer", "policy-v1").expect("valid"),
        expiry_ns: 10_000,
        generation: std::sync::atomic::AtomicU64::new(1),
    };
    assert!(matches!(
        grant.verify_fresh(&other, 5_000),
        Err(SessionError::GrantForged { .. })
    ));
}

#[test]
fn sg_004_core_leases_enforce_verbs_and_concurrency() {
    let (policy, grant) = granted();
    let book = CoreLeaseBook::new();
    // Ungranted verb refuses at acquisition.
    assert!(matches!(
        book.acquire(&grant, &policy, "topo.size", 1, 5_000),
        Err(SessionError::UngrantedVerb { .. })
    ));
    // Concurrency: 5 + 3 fits the 8-core grant; one more core refuses.
    let first = book
        .acquire(&grant, &policy, "flux.free-surface-lbm", 5, 5_000)
        .expect("first lease");
    let second = book
        .acquire(&grant, &policy, "ascent.optimize", 3, 5_000)
        .expect("second lease");
    assert_eq!(book.active_cores(SessionId(41)), 8);
    assert!(matches!(
        book.acquire(&grant, &policy, "ascent.optimize", 1, 5_000),
        Err(SessionError::CoreLeaseExceeded { .. })
    ));
    // Release returns capacity; a revoked grant cannot re-acquire.
    drop(first);
    assert_eq!(book.active_cores(SessionId(41)), 3);
    policy.revoke();
    assert!(matches!(
        book.acquire(&grant, &policy, "ascent.optimize", 1, 5_000),
        Err(SessionError::GrantRevoked { .. })
    ));
    drop(second);
    assert_eq!(book.active_cores(SessionId(41)), 0);
}

// ── Increment 2, slice A: grant-gated session open (bead aeq7) ──────────────

use fs_session::Governor;

/// sg-005: the granted open path registers ONLY the admitted view — the
/// caller supplies no token, and the registered token equals the grant's
/// admitted projection field-for-field.
#[test]
fn sg_005_granted_open_registers_only_the_admitted_view() {
    let policy = TestPolicy::new(1_000_000);
    let grant = mint_grant(&policy, &request(), 10).expect("policy admits");
    let governor = Governor::new();
    let open_id = governor
        .session_open_id(grant.session(), "sg-005-open")
        .expect("open id");
    governor
        .open_session_granted(open_id, &grant, &policy, 20)
        .expect("granted open");
    let registered = governor.token(grant.session()).expect("registered token");
    assert_eq!(registered, grant.admitted_token());
    assert_eq!(registered.ops, vec!["ascent.optimize", "flux.*"]);
    assert_eq!(registered.cores, 8);
}

/// sg-006: stale authority refuses AT the open boundary — expired,
/// rotated (revoked), and cross-issuer grants never reach registration.
#[test]
fn sg_006_stale_or_foreign_grants_refuse_at_the_open_boundary() {
    let governor = Governor::new();

    // Expired at open time.
    let policy = TestPolicy::new(100);
    let grant = mint_grant(&policy, &request(), 10).expect("policy admits");
    let open_id = governor
        .session_open_id(grant.session(), "sg-006-expired")
        .expect("open id");
    assert!(matches!(
        governor.open_session_granted(open_id, &grant, &policy, 100),
        Err(SessionError::GrantExpired { .. })
    ));

    // Policy rotation between mint and open (the rotation drill).
    let rotating = TestPolicy::new(1_000_000);
    let pre_rotation = mint_grant(&rotating, &request(), 10).expect("policy admits");
    rotating.revoke();
    let open_id = governor
        .session_open_id(pre_rotation.session(), "sg-006-rotated")
        .expect("open id");
    assert!(matches!(
        governor.open_session_granted(open_id, &pre_rotation, &rotating, 20),
        Err(SessionError::GrantRevoked { .. })
    ));

    // Cross-issuer presentation: a distinct issuer identity cannot vouch
    // for the grant at the open boundary.
    let issuer_a = TestPolicy::new(1_000_000);
    let issuer_b = TestPolicy {
        identity: IssuerIdentity::new("ops/other-issuer", "policy-v1").expect("valid"),
        expiry_ns: 1_000_000,
        generation: std::sync::atomic::AtomicU64::new(1),
    };
    let foreign = mint_grant(&issuer_a, &request(), 10).expect("policy admits");
    let open_id = governor
        .session_open_id(foreign.session(), "sg-006-foreign")
        .expect("open id");
    assert!(matches!(
        governor.open_session_granted(open_id, &foreign, &issuer_b, 20),
        Err(SessionError::GrantForged { .. })
    ));

    // Nothing registered by any refused path.
    assert!(matches!(
        governor.token(SessionId(41)),
        Err(SessionError::UnknownSession { .. })
    ));
}

/// sg-007: response replay — re-presenting the same grant (same session)
/// cannot open twice; session ids are single-use per governor.
#[test]
fn sg_007_grant_replay_cannot_open_a_second_session() {
    let policy = TestPolicy::new(1_000_000);
    let grant = mint_grant(&policy, &request(), 10).expect("policy admits");
    let governor = Governor::new();
    let open_id = governor
        .session_open_id(grant.session(), "sg-007-first")
        .expect("open id");
    governor
        .open_session_granted(open_id, &grant, &policy, 20)
        .expect("first open");

    // Same grant, fresh open id: replay refused.
    let replay_id = governor
        .session_open_id(grant.session(), "sg-007-replay")
        .expect("open id");
    assert!(matches!(
        governor.open_session_granted(replay_id, &grant, &policy, 30),
        Err(SessionError::SessionAlreadyOpen { .. })
    ));

    // Even a RE-MINTED grant for the same session id is refused: replay
    // protection is the governor's single-use session registry, not
    // grant-value identity.
    let reminted = mint_grant(&policy, &request(), 40).expect("policy admits");
    let reminted_id = governor
        .session_open_id(reminted.session(), "sg-007-reminted")
        .expect("open id");
    assert!(matches!(
        governor.open_session_granted(reminted_id, &reminted, &policy, 50),
        Err(SessionError::SessionAlreadyOpen { .. })
    ));
}

/// sg-008: wildcard confusion at the composed governor+lease boundary — a
/// `flux.*` grant covers namespaced verbs only; the bare namespace, sibling
/// namespaces sharing a prefix, and extended exact verbs all refuse at
/// lease acquisition on a session the governor actually opened.
#[test]
fn sg_008_wildcard_confusion_cannot_escalate_after_granted_open() {
    let policy = TestPolicy::new(1_000_000);
    let grant = mint_grant(&policy, &request(), 10).expect("policy admits");
    let governor = Governor::new();
    let open_id = governor
        .session_open_id(grant.session(), "sg-008-open")
        .expect("open id");
    governor
        .open_session_granted(open_id, &grant, &policy, 20)
        .expect("granted open");

    let book = CoreLeaseBook::new();
    for confused in ["flux", "fluxx.solve", "ascent.optimizer", "ascent"] {
        assert!(
            matches!(
                book.acquire(&grant, &policy, confused, 1, 30),
                Err(SessionError::UngrantedVerb { .. })
            ),
            "verb {confused:?} must not be covered"
        );
    }
    let _solve = book
        .acquire(&grant, &policy, "flux.solve", 4, 30)
        .expect("namespaced verb covered");
    let _exact = book
        .acquire(&grant, &policy, "ascent.optimize", 4, 30)
        .expect("exact verb covered");
    // And the concurrency ceiling still binds across both leases.
    assert!(matches!(
        book.acquire(&grant, &policy, "flux.assemble", 1, 30),
        Err(SessionError::CoreLeaseExceeded { .. })
    ));
}

/// sg-009: the fs-ir bridge — sealed_admission mints GrantBacked evidence
/// whose receipt names this exact grant; stale grants refuse; a receipt
/// naming a different grant digest is not vouched for; an inflated
/// capability projection is not vouched for.
#[test]
fn sg_009_sealed_admission_bridges_grant_authority_into_fs_ir() {
    use fs_ir::admission::{CapabilityEvidenceClass, CapabilityIssuerVerifier as _};
    use fs_session::GrantCapabilityVerifier;

    let policy = TestPolicy::new(1_000_000);
    let grant = mint_grant(&policy, &request(), 10).expect("policy admits");

    let sealed = grant
        .sealed_admission(&policy, 20)
        .expect("fresh grant seals");
    assert_eq!(
        sealed.evidence_class(),
        CapabilityEvidenceClass::GrantBacked
    );
    let receipt = sealed.receipt().expect("grant-backed carries receipt");
    assert_eq!(receipt.grant_digest, grant.digest());
    assert_eq!(receipt.issuer_id, "ops/test-issuer");
    assert_eq!(sealed.capability(), &grant.to_admission());

    // Stale authority refuses at the bridge.
    policy.revoke();
    assert!(matches!(
        grant.sealed_admission(&policy, 30),
        Err(SessionError::GrantRevoked { .. })
    ));

    // The verifier vouches only for its exact grant: a tampered digest or
    // an inflated projection is refused.
    let fresh_policy = TestPolicy::new(1_000_000);
    let fresh_grant = mint_grant(&fresh_policy, &request(), 10).expect("policy admits");
    let bridge = GrantCapabilityVerifier {
        grant: &fresh_grant,
        policy: &fresh_policy,
        now_ns: 20,
    };
    let mut tampered = fresh_grant.admission_receipt();
    tampered.grant_digest = "bb".repeat(32);
    assert!(!bridge.verify(&fresh_grant.to_admission(), &tampered));
    let mut inflated = fresh_grant.to_admission();
    inflated.cores = u64::MAX;
    assert!(!bridge.verify(&inflated, &fresh_grant.admission_receipt()));
}
