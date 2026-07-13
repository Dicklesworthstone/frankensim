//! Canonical replay-identity battery (bead gp3.14): field-mutation
//! coverage (every semantic field moves the root, every documented
//! exclusion does not), delimiter/split collision resistance, type
//! confusion, order sensitivity, float bit-pattern binding, child
//! (dependency) propagation, and fail-closed schema versioning.

use fs_obs::ident::{
    BoundedIdentityBuilder, IDENT_SCHEMA_VERSION, IdentError, IdentityBuildError, IdentityBuilder,
    ReplayIdentity, check_version,
};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-obs/ident\",\"case\":\"{case}\",\"verdict\":\"{}\",\"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
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
    // Floats bind by BIT PATTERN: -0.0 vs 0.0 differ; NaN is stable.
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
    verdict(
        "ident-004",
        types_distinct && ab != ba && z != nz && nan1 == nan2,
        "type tags, field order, and float bit patterns are all semantic; NaN is stable",
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
    let versioned_display = id.hex().starts_with("fsid-v1:x:");
    verdict(
        "ident-006",
        current_ok && future_refused && past_refused && versioned_display,
        "declared schema versions other than the supported one are refused (fail closed); \
         the version is framed into the root and the display form",
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
    let bounded = bounded_reference(&parent, exact_limit).expect("exact byte cap is admitted");
    assert_eq!(bounded, legacy);
    assert_eq!(bounded.exclusions(), legacy.exclusions());

    assert_eq!(
        bounded_reference(&parent, exact_limit - 1),
        Err(IdentityBuildError::CanonicalBytesExceeded {
            requested: exact_limit,
            limit: exact_limit - 1,
        })
    );

    // For kind "x" and key "k", the v1 header and field framing consume
    // 35 bytes before payload. This locks exact-limit and limit+1 behavior.
    let payload_limit = 39;
    let exact_payload = BoundedIdentityBuilder::new("x", payload_limit)
        .expect("header fits")
        .bytes("k", b"1234")
        .expect("four-byte payload reaches the exact cap")
        .finish();
    assert_eq!(exact_payload.canonical_bytes().len(), payload_limit);
    let payload_error = BoundedIdentityBuilder::new("x", payload_limit)
        .expect("header fits")
        .bytes("k", b"12345")
        .expect_err("limit+1 payload must refuse before mutation");
    assert_eq!(
        payload_error,
        IdentityBuildError::CanonicalBytesExceeded {
            requested: payload_limit + 1,
            limit: payload_limit,
        }
    );

    verdict(
        "ident-007",
        true,
        "bounded and legacy builders are byte-identical; exact caps admit and limit+1 refuses",
    );
}
