//! Certificate-plugin acceptance (bead 9e8n): two real solver-free
//! families round-trip producer → package bytes → standalone recheck;
//! seeded structurally-valid-but-mathematically-false witnesses fail
//! with localized refutations; bounds, tamper, version substitution,
//! and unknown families all refuse; and the verdict axes stay distinct.

use fs_checker::plugins::{
    MAX_CHAIN_OPS, MAX_RESIDUAL_DIM, PluginRegistry, PluginVerdict, encode_chain_witness,
    encode_residual_witness, witness_identity,
};

fn registry() -> PluginRegistry {
    PluginRegistry::v1()
}

/// Producer side of the chain family: ((a+b)·c) with a claimed enclosure
/// wide enough to contain the outward replay.
fn valid_chain() -> Vec<u8> {
    let inputs = [(1.0f64, 1.0 + 1e-12), (2.0, 2.0), (0.5, 0.5 + 1e-12)];
    // tape indices: 0,1,2 inputs; op0: add(0,1) -> 3; op1: mul(3,2) -> 4.
    let ops = [(1u8, 0u16, 1u16), (2, 3, 2)];
    // Exact value ≈ 1.5·(1+ε); claim a comfortably containing enclosure.
    encode_chain_witness(&inputs, &ops, (1.499_999_999, 1.500_000_001))
}

fn valid_residual() -> Vec<u8> {
    // A = [[2, 0], [0, 4]], x = [1, 0.5], b = [2, 2] → r = 0 exactly.
    encode_residual_witness(2, 2, &[2.0, 0.0, 0.0, 4.0], &[1.0, 0.5], &[2.0, 2.0], 1e-12)
}

#[test]
fn both_families_round_trip_producer_to_standalone_recheck() {
    let registry = registry();
    let chain = valid_chain();
    let verdict = registry.check("interval-enclosure-chain", 1, &chain);
    assert!(verdict.is_verified(), "chain: {verdict}");

    let residual = valid_residual();
    let verdict = registry.check("linear-residual-linf", 1, &residual);
    assert!(verdict.is_verified(), "residual: {verdict}");

    // Content identities are stable and distinct.
    assert_eq!(witness_identity(&chain), witness_identity(&valid_chain()));
    assert_ne!(witness_identity(&chain), witness_identity(&residual));
}

#[test]
fn mathematically_false_witnesses_fail_with_localized_refutations() {
    let registry = registry();

    // Structurally valid chain claiming a TIGHTER enclosure than the
    // outward replay can prove.
    let inputs = [(1.0f64, 1.0), (2.0, 2.0)];
    let ops = [(1u8, 0u16, 1u16)];
    let lying = encode_chain_witness(&inputs, &ops, (3.0, 3.0));
    match registry.check("interval-enclosure-chain", 1, &lying) {
        PluginVerdict::SemanticallyRefuted { location, detail } => {
            assert!(location.contains("op 0"), "localized: {location}");
            assert!(detail.contains("does not contain"));
        }
        other => panic!("tight false claim must refute semantically: {other}"),
    }

    // Structurally valid residual with a bound the true residual breaks:
    // r = b − Ax = [0.5], claimed 1e-6.
    let false_residual = encode_residual_witness(1, 1, &[1.0], &[1.0], &[1.5], 1e-6);
    match registry.check("linear-residual-linf", 1, &false_residual) {
        PluginVerdict::SemanticallyRefuted { location, detail } => {
            assert_eq!(location, "row 0");
            assert!(detail.contains("not provably inside"));
        }
        other => panic!("false bound must refute at the row: {other}"),
    }
}

#[test]
fn unknown_families_and_version_substitution_refuse_explicitly() {
    let registry = registry();
    let chain = valid_chain();

    match registry.check("interval-enclosure-chain", 2, &chain) {
        PluginVerdict::CapabilityRefused {
            family, version, ..
        } => {
            assert_eq!(family, "interval-enclosure-chain");
            assert_eq!(version, 2);
        }
        other => panic!("version substitution must refuse: {other}"),
    }
    match registry.check("spectral-gap-certificate", 1, &chain) {
        PluginVerdict::CapabilityRefused { detail, .. } => {
            assert!(detail.contains("independent semantic verification is unavailable"));
        }
        other => panic!("unknown family must refuse, never Pass: {other}"),
    }
    assert_eq!(registry.families().len(), 2, "v1 registry is closed");
}

#[test]
fn witness_field_tamper_never_verifies_and_changes_identity() {
    let registry = registry();
    let baseline = valid_residual();
    let baseline_id = witness_identity(&baseline);
    let baseline_verdict = registry.check("linear-residual-linf", 1, &baseline);
    assert!(baseline_verdict.is_verified());

    for index in 0..baseline.len() {
        let mut tampered = baseline.clone();
        tampered[index] ^= 0x01;
        let verdict = registry.check("linear-residual-linf", 1, &tampered);
        // A tampered witness may become malformed, refuted, or even
        // still verify a DIFFERENT (weaker) true claim — but it can
        // never verify while pretending to be the original witness.
        if verdict.is_verified() {
            assert_ne!(
                witness_identity(&tampered),
                baseline_id,
                "byte {index}: tamper must change the content identity"
            );
        }
    }
}

#[test]
fn bounds_are_enforced_before_any_arithmetic() {
    let registry = registry();

    let oversized_dim = encode_residual_witness(
        MAX_RESIDUAL_DIM + 1,
        1,
        &vec![0.0; MAX_RESIDUAL_DIM + 1],
        &[0.0],
        &vec![0.0; MAX_RESIDUAL_DIM + 1],
        0.0,
    );
    match registry.check("linear-residual-linf", 1, &oversized_dim) {
        PluginVerdict::WitnessMalformed { detail, .. } => {
            assert!(detail.contains("outside 1..="));
        }
        other => panic!("oversized dimension must be malformed: {other}"),
    }

    let mut runaway_ops = Vec::new();
    for _ in 0..=MAX_CHAIN_OPS {
        runaway_ops.push((1u8, 0u16, 0u16));
    }
    let oversized_chain = encode_chain_witness(&[(0.0, 0.0)], &runaway_ops, (0.0, 0.0));
    match registry.check("interval-enclosure-chain", 1, &oversized_chain) {
        PluginVerdict::WitnessMalformed { detail, .. } => {
            assert!(detail.contains("op count"));
        }
        other => panic!("oversized op tape must be malformed: {other}"),
    }

    let giant = vec![0u8; fs_checker::plugins::MAX_WITNESS_BYTES + 1];
    match registry.check("interval-enclosure-chain", 1, &giant) {
        PluginVerdict::WitnessMalformed { detail, .. } => {
            assert!(detail.contains("byte bound"));
        }
        other => panic!("oversized transport must be malformed: {other}"),
    }

    // Non-finite and inverted inputs refuse structurally.
    let nan_input = encode_chain_witness(&[(f64::NAN, 1.0)], &[(3, 0, 0)], (0.0, 1.0));
    assert!(matches!(
        registry.check("interval-enclosure-chain", 1, &nan_input),
        PluginVerdict::WitnessMalformed { .. }
    ));
    let inverted = encode_chain_witness(&[(2.0, 1.0)], &[(3, 0, 0)], (0.0, 1.0));
    assert!(matches!(
        registry.check("interval-enclosure-chain", 1, &inverted),
        PluginVerdict::WitnessMalformed { .. }
    ));
    let dangling_operand = encode_chain_witness(&[(0.0, 1.0)], &[(1, 0, 7)], (0.0, 2.0));
    assert!(matches!(
        registry.check("interval-enclosure-chain", 1, &dangling_operand),
        PluginVerdict::WitnessMalformed { .. }
    ));
}

#[test]
fn verdict_axes_stay_distinct() {
    // The four outcomes are four types; only one grants verification —
    // integrity (IntegrityStatus) and origin (OriginStatus) authority
    // live on the existing checker surfaces, and the plugin verdict
    // supplies ONLY the independent-semantic-verification axis.
    let verified = PluginVerdict::SemanticallyVerified {
        detail: String::new(),
    };
    let refuted = PluginVerdict::SemanticallyRefuted {
        location: String::new(),
        detail: String::new(),
    };
    let malformed = PluginVerdict::WitnessMalformed {
        offset: 0,
        detail: String::new(),
    };
    let refused = PluginVerdict::CapabilityRefused {
        family: String::new(),
        version: 0,
        detail: String::new(),
    };
    assert!(verified.is_verified());
    for other in [&refuted, &malformed, &refused] {
        assert!(!other.is_verified());
    }
    assert_ne!(verified, refuted);
    assert_ne!(refuted, malformed);
    assert_ne!(malformed, refused);
}
