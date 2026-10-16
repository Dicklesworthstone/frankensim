//! Canonical replay-identity battery (bead gp3.14): field-mutation
//! coverage (every semantic field moves the root, every documented
//! exclusion does not), delimiter/split collision resistance, type
//! confusion, order sensitivity, float bit-pattern binding, child
//! (dependency) propagation, and fail-closed schema versioning.

use fs_obs::ident::{
    BoundedIdentityBuilder, IDENT_SCHEMA_VERSION, IdentError, IdentityBuildError, IdentityBuilder,
    REPLAY_IDENTITY_DOMAIN, ReplayIdentity, check_version,
};

fn verdict(case: &str, pass: bool, detail: &str) {
    let mut emitter = fs_obs::Emitter::new("fs-obs/ident", case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: "fs-obs/ident".to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed: 0,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("identity verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("identity verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

/// The reference recipe: one field of every type, in a fixed order —
/// the mutation battery flips each in turn.
#[allow(clippy::too_many_arguments)] // one arg per field type, by design
fn recipe(
    kind: &str,
    algo: &str,
    n: u64,
    offset: i64,
    tol: f64,
    payload: &[u8],
    det: bool,
    parent: u64,
) -> u64 {
    IdentityBuilder::new(kind)
        .str("algorithm", algo)
        .u64("n", n)
        .i64("offset", offset)
        .f64_bits("tolerance", tol)
        .bytes("payload", payload)
        .flag("deterministic", det)
        .child_root64("parent", parent)
        .exclude("wall_ns", "timing is evidence, not identity")
        .finish()
        .root()
}

#[test]
fn ident_001_every_semantic_field_moves_the_root() {
    let base = recipe("test-artifact", "gemm", 64, -3, 1e-9, b"pk", true, 7);
    let mutations = [
        (
            "kind",
            recipe("test-artifact2", "gemm", 64, -3, 1e-9, b"pk", true, 7),
        ),
        (
            "str",
            recipe("test-artifact", "gemm2", 64, -3, 1e-9, b"pk", true, 7),
        ),
        (
            "u64",
            recipe("test-artifact", "gemm", 65, -3, 1e-9, b"pk", true, 7),
        ),
        (
            "i64",
            recipe("test-artifact", "gemm", 64, -4, 1e-9, b"pk", true, 7),
        ),
        (
            "f64",
            recipe("test-artifact", "gemm", 64, -3, 2e-9, b"pk", true, 7),
        ),
        (
            "bytes",
            recipe("test-artifact", "gemm", 64, -3, 1e-9, b"pq", true, 7),
        ),
        (
            "flag",
            recipe("test-artifact", "gemm", 64, -3, 1e-9, b"pk", false, 7),
        ),
        (
            "child",
            recipe("test-artifact", "gemm", 64, -3, 1e-9, b"pk", true, 8),
        ),
    ];
    let all_moved = mutations.iter().all(|(_, r)| *r != base);
    let all_distinct = {
        let mut roots: Vec<u64> = mutations.iter().map(|(_, r)| *r).collect();
        roots.push(base);
        roots.sort_unstable();
        roots.windows(2).all(|w| w[0] != w[1])
    };
    // Determinism: the same recipe binds the same root.
    let replay = recipe("test-artifact", "gemm", 64, -3, 1e-9, b"pk", true, 7);
    verdict(
        "ident-001",
        all_moved && all_distinct && replay == base,
        "every semantic field (incl kind and child root) moves the root; all mutants \
         pairwise distinct; replay is bit-stable",
    );
}

#[test]
fn ident_002_documented_exclusions_do_not_move_the_root() {
    let a = IdentityBuilder::new("x")
        .u64("n", 1)
        .exclude("wall_ns", "timing is evidence, not identity")
        .finish();
    let b = IdentityBuilder::new("x").u64("n", 1).finish();
    let c = IdentityBuilder::new("x")
        .u64("n", 1)
        .exclude("hostname", "machine names are provenance-only")
        .exclude("wall_ns", "timing is evidence, not identity")
        .finish();
    let documented = a.exclusions().len() == 1 && c.exclusions().len() == 2;
    verdict(
        "ident-002",
        a.root() == b.root() && a.root() == c.root() && documented,
        "exclude() never moves the root and the exclusion list is the auditable record",
    );
}

#[test]
fn ident_003_delimiter_and_split_collisions_refused() {
    // The classic under-binding: ("ab","c") vs ("a","bc") — identical
    // under concatenation, distinct under length prefixes.
    let ab_c = IdentityBuilder::new("x")
        .str("k1", "ab")
        .str("k2", "c")
        .finish()
        .root();
    let a_bc = IdentityBuilder::new("x")
        .str("k1", "a")
        .str("k2", "bc")
        .finish()
        .root();
    // Adversarial delimiter injection: a value CONTAINING the old "|"
    // separator cannot imitate two fields.
    let joined = IdentityBuilder::new("x").str("k1", "a|b").finish().root();
    let split = IdentityBuilder::new("x")
        .str("k1", "a")
        .str("k1", "b")
        .finish()
        .root();
    // Key/value boundary attack: key "ab" value "" vs key "a" value "b".
    let kv1 = IdentityBuilder::new("x").str("ab", "").finish().root();
    let kv2 = IdentityBuilder::new("x").str("a", "b").finish().root();
    // Empty value vs absent field.
    let empty = IdentityBuilder::new("x").str("k", "").finish().root();
    let absent = IdentityBuilder::new("x").finish().root();
    verdict(
        "ident-003",
        ab_c != a_bc && joined != split && kv1 != kv2 && empty != absent,
        "length prefixes refuse split, delimiter-injection, key/value-boundary, and \
         empty-vs-absent collisions",
    );
}

#[test]
fn ident_004_type_confusion_and_order_are_semantic() {
    // The same logical "1" under four types: four distinct roots.
    let s = IdentityBuilder::new("x").str("k", "1").finish().root();
    let u = IdentityBuilder::new("x").u64("k", 1).finish().root();
    let i = IdentityBuilder::new("x").i64("k", 1).finish().root();
    let b = IdentityBuilder::new("x").bytes("k", b"1").finish().root();
    let mut four = [s, u, i, b];
    four.sort_unstable();
    let types_distinct = four.windows(2).all(|w| w[0] != w[1]);
    // Field order is part of the identity (ordered operation params).
    let ab = IdentityBuilder::new("x")
        .u64("a", 1)
        .u64("b", 2)
        .finish()
        .root();
    let ba = IdentityBuilder::new("x")
        .u64("b", 2)
        .u64("a", 1)
        .finish()
        .root();
    // Floats bind by BIT PATTERN: -0.0 vs 0.0 differ; identical NaNs are
    // stable, while distinct payloads remain distinct even when display text
    // collapses them.
    let z = IdentityBuilder::new("x").f64_bits("v", 0.0).finish().root();
    let nz = IdentityBuilder::new("x")
        .f64_bits("v", -0.0)
        .finish()
        .root();
    let nan1 = IdentityBuilder::new("x")
        .f64_bits("v", f64::NAN)
        .finish()
        .root();
    let nan2 = IdentityBuilder::new("x")
        .f64_bits("v", f64::NAN)
        .finish()
        .root();
    let display_nan_a = f64::from_bits(0x7ff8_0000_0000_0001);
    let display_nan_b = f64::from_bits(0x7ff8_0000_0000_0002);
    let nan_displays_match = display_nan_a.to_string() == display_nan_b.to_string();
    let nan_payload_a = IdentityBuilder::new("x")
        .f64_bits("v", display_nan_a)
        .finish();
    let nan_payload_b = IdentityBuilder::new("x")
        .f64_bits("v", display_nan_b)
        .finish();
    verdict(
        "ident-004",
        types_distinct
            && ab != ba
            && z != nz
            && nan1 == nan2
            && nan_displays_match
            && nan_payload_a.root() != nan_payload_b.root()
            && nan_payload_a.canonical_bytes() != nan_payload_b.canonical_bytes(),
        "type tags, field order, and exact float bits are semantic; same-display NaN payloads remain distinct",
    );
}

#[test]
fn ident_005_dependency_children_propagate() {
    // A downstream identity that binds an upstream root changes when
    // the upstream changes — the dependency-aware edge that names
    // affected goldens.
    let up_v1 = IdentityBuilder::new("kernel").u64("kc", 256).finish();
    let up_v2 = IdentityBuilder::new("kernel").u64("kc", 128).finish();
    let down_v1 = IdentityBuilder::new("golden")
        .str("suite", "gemm")
        .child("kernel", &up_v1)
        .finish();
    let down_v2 = IdentityBuilder::new("golden")
        .str("suite", "gemm")
        .child("kernel", &up_v2)
        .finish();
    verdict(
        "ident-005",
        up_v1.root() != up_v2.root() && down_v1.root() != down_v2.root(),
        "an upstream identity change propagates through child() into every downstream root",
    );
}

#[test]
fn ident_006_unknown_schema_versions_fail_closed() {
    let current_ok = check_version(IDENT_SCHEMA_VERSION).is_ok();
    let future = check_version(IDENT_SCHEMA_VERSION + 1);
    let past = check_version(0);
    let future_refused = matches!(
        future,
        Err(IdentError::UnknownSchemaVersion { declared, .. }) if declared == IDENT_SCHEMA_VERSION + 1
    );
    let past_refused = matches!(past, Err(IdentError::UnknownSchemaVersion { .. }));
    // The version is also part of the hashed frame: an identity's hex
    // form names it explicitly.
    let id = IdentityBuilder::new("x").finish();
    let versioned_display = id
        .hex()
        .starts_with(&format!("{REPLAY_IDENTITY_DOMAIN}-v1:x:"));
    let root_matches_canonical_bytes = id.root() == fs_obs::fnv1a64(id.canonical_bytes());
    let domain_prefix_matches = id
        .canonical_bytes()
        .starts_with(REPLAY_IDENTITY_DOMAIN.as_bytes());
    let version_offset = REPLAY_IDENTITY_DOMAIN.len();
    let mut future_bytes = id.canonical_bytes().to_vec();
    let version_slot_present = if let Some(version_bytes) =
        future_bytes.get_mut(version_offset..version_offset + core::mem::size_of::<u32>())
    {
        version_bytes.copy_from_slice(&(IDENT_SCHEMA_VERSION + 1).to_le_bytes());
        true
    } else {
        false
    };
    let version_bytes_move_root =
        version_slot_present && fs_obs::fnv1a64(&future_bytes) != id.root();
    let mut foreign_domain_bytes = id.canonical_bytes().to_vec();
    let domain_slot_present = if let Some(first) = foreign_domain_bytes.first_mut() {
        *first ^= 1;
        true
    } else {
        false
    };
    let domain_bytes_move_root =
        domain_slot_present && fs_obs::fnv1a64(&foreign_domain_bytes) != id.root();
    verdict(
        "ident-006",
        current_ok
            && future_refused
            && past_refused
            && versioned_display
            && root_matches_canonical_bytes
            && domain_prefix_matches
            && version_bytes_move_root
            && domain_bytes_move_root,
        "declared schema versions other than the supported one are refused (fail closed); \
         the declared domain and version bytes are framed into the root and display form",
    );
}

fn bounded_reference(
    parent: &ReplayIdentity,
    max_canonical_bytes: usize,
) -> Result<ReplayIdentity, IdentityBuildError> {
    let builder = BoundedIdentityBuilder::new("bounded-artifact", max_canonical_bytes)?;
    let builder = builder.str("algorithm", "gemm")?;
    let builder = builder.u64("n", 64)?;
    let builder = builder.i64("offset", -3)?;
    let builder = builder.f64_bits("tolerance", -0.0)?;
    let builder = builder.bytes("payload", b"a|b\0c")?;
    let builder = builder.flag("deterministic", true)?;
    let builder = builder.child("upstream", parent)?;
    let builder = builder.child_root64("legacy", 0x0123_4567_89ab_cdef)?;
    let builder = builder.exclude("wall_ns", "timing is evidence, not identity")?;
    Ok(builder.finish())
}

#[test]
fn ident_007_bounded_builder_is_byte_exact_and_refuses_at_the_cap() {
    let parent = IdentityBuilder::new("kernel").u64("kc", 256).finish();
    let legacy = IdentityBuilder::new("bounded-artifact")
        .str("algorithm", "gemm")
        .u64("n", 64)
        .i64("offset", -3)
        .f64_bits("tolerance", -0.0)
        .bytes("payload", b"a|b\0c")
        .flag("deterministic", true)
        .child("upstream", &parent)
        .child_root64("legacy", 0x0123_4567_89ab_cdef)
        .exclude("wall_ns", "timing is evidence, not identity")
        .finish();
    let exact_limit = legacy.canonical_bytes().len();
    let bounded = bounded_reference(&parent, exact_limit);
    let exact_cap_matches = bounded
        .as_ref()
        .is_ok_and(|identity| identity == &legacy && identity.exclusions() == legacy.exclusions());
    let roomier = bounded_reference(&parent, exact_limit.saturating_add(64));
    let roomier_matches = match (&bounded, &roomier) {
        (Ok(exact), Ok(extra)) => {
            extra == &legacy
                && extra.canonical_bytes() == exact.canonical_bytes()
                && extra.root() == exact.root()
        }
        _ => false,
    };
    let limit_minus_one_refuses = exact_limit.checked_sub(1).is_some_and(|limit| {
        bounded_reference(&parent, limit)
            == Err(IdentityBuildError::CanonicalBytesExceeded {
                requested: exact_limit,
                limit,
            })
    });

    // Derive the header directly from the declared replay domain so this cap
    // cannot preserve an obsolete split-magic assumption.
    let framing_without_payload = REPLAY_IDENTITY_DOMAIN.len()
        + core::mem::size_of::<u32>()
        + core::mem::size_of::<u64>()
        + "x".len()
        + core::mem::size_of::<u8>()
        + core::mem::size_of::<u64>()
        + "k".len()
        + core::mem::size_of::<u64>();
    let payload_limit = framing_without_payload + 4;
    let exact_payload = BoundedIdentityBuilder::new("x", payload_limit)
        .and_then(|builder| builder.bytes("k", b"1234"))
        .map(BoundedIdentityBuilder::finish);
    let exact_payload_fills_cap = exact_payload
        .as_ref()
        .is_ok_and(|identity| identity.canonical_bytes().len() == payload_limit);
    let payload_refuses_limit_plus_one = BoundedIdentityBuilder::new("x", payload_limit)
        .and_then(|builder| builder.bytes("k", b"12345"))
        .map(BoundedIdentityBuilder::finish)
        == Err(IdentityBuildError::CanonicalBytesExceeded {
            requested: payload_limit + 1,
            limit: payload_limit,
        });

    verdict(
        "ident-007",
        exact_cap_matches
            && roomier_matches
            && limit_minus_one_refuses
            && exact_payload_fills_cap
            && payload_refuses_limit_plus_one,
        &format!(
            "bounded builder checks: exact_cap_matches={exact_cap_matches}, \
             roomier_matches={roomier_matches}, limit_minus_one_refuses={limit_minus_one_refuses}, \
             exact_payload_fills_cap={exact_payload_fills_cap}, \
             payload_refuses_limit_plus_one={payload_refuses_limit_plus_one}"
        ),
    );
}
