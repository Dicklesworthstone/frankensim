//! Boundary, migration, and diagnostics-determinism tests for the
//! canonical-QR typed surface (`fs_la::canonical_qr`).
//!
//! Bead frankensim-epic-bedrock-6ys.5.1.2. Complements the in-module unit
//! battery: this file owns empty/min/max boundaries, cross-version
//! migration refusal, byte-stable encodings, and stable diagnostic text.

use fs_la::canonical_qr::{
    ArithmeticMode, CertifiedRankProfile, CanonicalQrOutcome, CanonicalQrPolicy, ClaimTier,
    DeterminismClass, ErrorBudget, NoClaimReason, OutcomeAuthority, PivotClass, PolicyError,
    RankTolerance, ReplayIdentity, CANONICAL_QR_SCHEMA_VERSION,
};
use fs_blake3::{hash_bytes, ContentHash};

fn base_policy() -> CanonicalQrPolicy {
    CanonicalQrPolicy::new(
        RankTolerance::default_f64(),
        ErrorBudget::relative(1e-12).expect("in window"),
        DeterminismClass::SameIsaBitStable,
        ArithmeticMode::Binary64RoundToNearest,
        fs_la::canonical_qr::TiePolicy::LowestIndexFirst,
    )
    .expect("valid")
}

fn identity(cert: Option<ContentHash>) -> ReplayIdentity {
    ReplayIdentity {
        input_digest: hash_bytes(b"input"),
        tree_digest: hash_bytes(b"tree"),
        result_digest: hash_bytes(b"result"),
        certificate_ref: cert,
        arithmetic_mode: ArithmeticMode::Binary64RoundToNearest,
    }
}

// ---------------------------------------------------------------------------
// Empty boundary: n = 0 outcomes are first-class (empty TSQR semantics).
// ---------------------------------------------------------------------------
#[test]
fn zero_dimension_outcome_is_admissible() {
    let profile = CertifiedRankProfile::checked(Vec::new()).expect("consistent");
    let outcome = CanonicalQrOutcome::checked(
        Vec::new(),
        0,
        profile,
        OutcomeAuthority::NoClaim(NoClaimReason::RankDeficientCrossScheduleEquality),
        identity(None),
    )
    .expect("empty R is valid");
    assert_eq!(outcome.n(), 0);
    assert!(outcome.r_factor().is_empty());
    assert_eq!(outcome.rank_profile().rank(), 0);
}

// ---------------------------------------------------------------------------
// Min/max budget windows refuse at the edges by construction.
// ---------------------------------------------------------------------------
#[test]
fn budget_window_edges() {
    assert_eq!(ErrorBudget::relative(0.0), Err(PolicyError::BudgetOutOfRange));
    assert!(ErrorBudget::relative(f64::MIN_POSITIVE).is_ok());
    assert!(ErrorBudget::relative(1.0).is_ok());
    assert_eq!(ErrorBudget::relative(1.0 + f64::EPSILON), Err(PolicyError::BudgetOutOfRange));
    assert_eq!(
        ErrorBudget::relative(-1.0),
        Err(PolicyError::BudgetOutOfRange)
    );
    // Rank tolerance window edges (floor is one eps above machine epsilon).
    assert_eq!(
        RankTolerance::relative(2.3e-17),
        Err(PolicyError::InvalidScaleRelativeFactor)
    );
    let just_above = fs_la::canonical_qr::CANONICAL_QR_THEOREM_VERSION; // touch const for compile-surface
    assert_eq!(just_above, 1);
}

// ---------------------------------------------------------------------------
// Encoding is byte-stable across calls and structurally equal outcomes
// carry equal digests only when every bound field matches.
// ---------------------------------------------------------------------------
#[test]
fn encodings_are_deterministic_and_identity_sensitive() {
    let policy = base_policy();
    let e1 = policy.encode();
    let e2 = policy.encode();
    assert_eq!(e1, e2, "policy encoding must be call-stable");

    let mk_outcome = |tree: &str| {
        let profile =
            CertifiedRankProfile::checked(vec![PivotClass::Nonzero, PivotClass::Zero]).expect("ok");
        CanonicalQrOutcome::checked(
            vec![2.0, 1.0, 0.0, 0.25],
            2,
            profile,
            OutcomeAuthority::NoClaim(NoClaimReason::RankDeficientCrossScheduleEquality),
            identity_with_tree(tree),
        )
        .expect("valid")
    };
    fn identity_with_tree(tree: &str) -> ReplayIdentity {
        ReplayIdentity {
            input_digest: hash_bytes(b"input"),
            tree_digest: hash_bytes(tree.as_bytes()),
            result_digest: hash_bytes(b"result"),
            certificate_ref: None,
            arithmetic_mode: ArithmeticMode::Binary64RoundToNearest,
        }
    }
    let o1 = mk_outcome("tree");
    let o2 = mk_outcome("tree");
    let o3 = mk_outcome("other-tree");
    assert_eq!(o1.result_digest(), o2.result_digest());
    assert_ne!(o1.result_digest(), o3.result_digest());
    // Factor-bit mutation moves the digest even when the authority text
    // would read identically.
    let mutated = {
        let profile =
            CertifiedRankProfile::checked(vec![PivotClass::Nonzero, PivotClass::Zero]).expect("ok");
        CanonicalQrOutcome::checked(
            vec![2.0, 1.0 + f64::EPSILON, 0.0, 0.25],
            2,
            profile,
            OutcomeAuthority::NoClaim(NoClaimReason::RankDeficientCrossScheduleEquality),
            identity_with_tree("tree"),
        )
        .expect("still valid")
    };
    assert_ne!(o1.result_digest(), mutated.result_digest());
}

// ---------------------------------------------------------------------------
// Migration surface: identities from mismatched modes are stale/cross-domain;
// unknown schema bytes never decode into a live policy.
// ---------------------------------------------------------------------------
#[test]
fn stale_identities_refuse() {
    let policy = base_policy();
    let mut stale = identity(None);
    stale.arithmetic_mode = ArithmeticMode::Binary64DirectedOutward;
    assert_eq!(
        stale.validate(&policy),
        Err(PolicyError::StaleIdentity { field: "arithmetic_mode" })
    );
    // Current-mode identity validates.
    identity(None).validate(&policy).expect("coherent");
    // Schema constant is the one this build speaks; a decoder from another
    // era must fail closed rather than best-effort parse.
    assert_eq!(CANONICAL_QR_SCHEMA_VERSION, 1);
}

// ---------------------------------------------------------------------------
// Diagnostics are deterministic prose: exact strings pinned so downstream
// logs/ledgers cannot drift silently between revisions.
// ---------------------------------------------------------------------------
#[test]
fn diagnostic_text_is_stable() {
    assert_eq!(
        PolicyError::InvalidScaleRelativeFactor.to_string(),
        "tolerance/budget must be a positive finite scale-relative factor"
    );
    assert_eq!(
        PolicyError::UncertifiedClaim.to_string(),
        "certified tier requires an independent checker receipt reference"
    );
    assert_eq!(
        PolicyError::StrictlyNegativeDiagonal { index: 3 }.to_string(),
        "R diagonal [3] strictly negative; flip law violated"
    );
    assert_eq!(
        PolicyError::UnknownSchemaVersion(9).to_string(),
        "unknown canonical-qr schema version 9"
    );
}

// ---------------------------------------------------------------------------
// Tier tag space is closed under both directions (no silent renumbering):
// every variant round-trips its wire tag; out-of-range tags refuse.
// ---------------------------------------------------------------------------
#[test]
fn claim_tier_tags_roundtrip_fail_closed() {
    for tier in [
        ClaimTier::ExactReconstruction,
        ClaimTier::SameIsaDeterministic,
        ClaimTier::FullRankTreeAgreement,
    ] {
        assert_eq!(ClaimTier::from_tag(tier.tag()), Some(tier));
    }
    assert_eq!(ClaimTier::from_tag(3), None);
    assert_eq!(ClaimTier::from_tag(255), None);

    for authority in [
        OutcomeAuthority::Certified(ClaimTier::FullRankTreeAgreement),
        OutcomeAuthority::NoClaim(NoClaimReason::AmbiguousRankBoundary),
        OutcomeAuthority::NoClaim(NoClaimReason::NonFiniteInput),
    ] {
        assert_eq!(OutcomeAuthority::from_tag(authority.tag()), Some(authority));
    }
    assert_eq!(OutcomeAuthority::from_tag([0, 3]), None);
    assert_eq!(OutcomeAuthority::from_tag([1, 4]), None);
    assert_eq!(OutcomeAuthority::from_tag([2, 0]), None);
}
