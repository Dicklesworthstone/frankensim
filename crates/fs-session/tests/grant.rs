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
