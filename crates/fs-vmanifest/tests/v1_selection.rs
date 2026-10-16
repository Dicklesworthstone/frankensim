//! V.1.5 deterministic selection conformance (bead i94v.7.1.5): stratum
//! and profile schemas, composition validity, permutation/shard laws,
//! filter algebra, capability routing, empty-selection refusal, legacy
//! alias refusal, and exact semantic diffs under mutation.

use fs_vmanifest::v1::CaseId;
use fs_vmanifest::v1_selection::{
    BuiltinProfile, CompositeProfile, Filter, NamedSkip, ProfileId, SelectableCase, SelectionInput,
    Stratum, expand_selection, semantic_diff,
};
use std::collections::{BTreeMap, BTreeSet};

fn case(id: &str, stratum: Stratum, profiles: &[BuiltinProfile], caps: &[&str]) -> SelectableCase {
    SelectableCase {
        case: CaseId::new(id).expect("case id admits"),
        stratum,
        profiles: profiles.iter().copied().collect(),
        required_capabilities: caps.iter().map(|c| (*c).to_owned()).collect(),
    }
}

fn corpus() -> Vec<SelectableCase> {
    vec![
        case(
            "core/alpha",
            Stratum::Core,
            &[BuiltinProfile::Smoke, BuiltinProfile::Standard],
            &[],
        ),
        case("core/beta", Stratum::Core, &[BuiltinProfile::Standard], &[]),
        case(
            "core/gamma-isa",
            Stratum::Core,
            &[BuiltinProfile::CrossIsa],
            &["x86-64-host"],
        ),
        case(
            "max/delta",
            Stratum::Max,
            &[BuiltinProfile::Standard, BuiltinProfile::Soak],
            &[],
        ),
        case(
            "max/epsilon",
            Stratum::Max,
            &[BuiltinProfile::Adversarial],
            &[],
        ),
        case(
            "core/zeta-slow",
            Stratum::Core,
            &[BuiltinProfile::Standard],
            &[],
        ),
    ]
}

fn input(stratum: Stratum, profile: BuiltinProfile) -> SelectionInput {
    SelectionInput {
        stratum,
        profile: ProfileId::Builtin(profile),
        filters: Vec::new(),
        capabilities: BTreeSet::new(),
        budgets: BTreeMap::from([("time-ns".to_owned(), 60_000_000_000)]),
        skips: Vec::new(),
        shards: 1,
    }
}

#[test]
fn selection_is_invariant_under_enumeration_order_and_presentation() {
    let mut forward = corpus();
    let receipt_forward =
        expand_selection(&forward, &input(Stratum::Max, BuiltinProfile::Standard))
            .expect("selection admits");
    forward.reverse();
    let receipt_reverse =
        expand_selection(&forward, &input(Stratum::Max, BuiltinProfile::Standard))
            .expect("reversed enumeration admits");
    assert_eq!(receipt_forward, receipt_reverse);
    assert_eq!(receipt_forward.digest(), receipt_reverse.digest());
    assert_eq!(
        receipt_forward
            .selected
            .iter()
            .map(CaseId::as_str)
            .collect::<Vec<_>>(),
        ["core/alpha", "core/beta", "core/zeta-slow", "max/delta"]
    );
}

#[test]
fn stratum_selects_scientific_scope_and_profile_selects_intensity() {
    let core = expand_selection(&corpus(), &input(Stratum::Core, BuiltinProfile::Standard))
        .expect("core admits");
    assert!(
        core.selected
            .iter()
            .all(|c| c.as_str().starts_with("core/"))
    );

    let max = expand_selection(&corpus(), &input(Stratum::Max, BuiltinProfile::Standard))
        .expect("max admits");
    assert!(
        max.selected.iter().any(|c| c.as_str() == "max/delta"),
        "max includes the max surface"
    );
    assert!(
        max.selected.iter().any(|c| c.as_str() == "core/beta"),
        "max includes the core surface"
    );

    // Profile change (intensity) cannot smuggle a scope change: the
    // adversarial run at CORE stays empty of max-only cases — and here
    // is EXPLICITLY non-green because nothing at core is adversarial.
    let err = expand_selection(
        &corpus(),
        &input(Stratum::Core, BuiltinProfile::Adversarial),
    )
    .expect_err("empty selection is non-green");
    assert_eq!(err.rule(), "v1-empty-selection");
}

#[test]
fn shards_partition_the_selection_exactly() {
    for shards in [1u32, 2, 3, 5, 8] {
        let mut sel = input(Stratum::Max, BuiltinProfile::Standard);
        sel.shards = shards;
        let receipt = expand_selection(&corpus(), &sel).expect("admits");
        let mut union: Vec<&CaseId> = receipt.shards.iter().flatten().collect();
        union.sort();
        let mut expected: Vec<&CaseId> = receipt.selected.iter().collect();
        expected.sort();
        assert_eq!(union, expected, "shard union = selection at {shards}");
        let total: usize = receipt.shards.iter().map(Vec::len).sum();
        assert_eq!(total, receipt.selected.len(), "disjoint at {shards}");
    }
}

#[test]
fn filter_algebra_applies_in_declared_order() {
    let mut sel = input(Stratum::Max, BuiltinProfile::Standard);
    sel.filters = vec![
        Filter::ExcludePrefix("core/".to_owned()),
        Filter::IncludePrefix("core/zeta".to_owned()),
    ];
    let receipt = expand_selection(&corpus(), &sel).expect("admits");
    assert_eq!(
        receipt
            .selected
            .iter()
            .map(CaseId::as_str)
            .collect::<Vec<_>>(),
        ["core/zeta-slow", "max/delta"],
        "later include re-admits over the earlier exclude"
    );
}

#[test]
fn capability_and_predicate_skips_are_visible_and_named() {
    let mut sel = input(Stratum::Max, BuiltinProfile::CrossIsa);
    sel.skips = vec![NamedSkip {
        name: "quiet-host-required".to_owned(),
        prefix: "core/gamma".to_owned(),
        reason: "perf lane needs a quiet window".to_owned(),
    }];
    // No capabilities supplied: gamma-isa needs x86-64-host.
    let err = expand_selection(&corpus(), &sel).expect_err("only case skipped => empty");
    assert_eq!(err.rule(), "v1-empty-selection");

    // With the capability supplied, the named predicate still skips it —
    // visibly.
    sel.capabilities = BTreeSet::from(["x86-64-host".to_owned()]);
    let err = expand_selection(&corpus(), &sel).expect_err("still empty, but reasoned");
    assert_eq!(err.rule(), "v1-empty-selection");

    // Drop the predicate: the case now runs; remove the capability
    // instead and the skip is NAMED capability:<missing>.
    sel.skips.clear();
    let receipt = expand_selection(&corpus(), &sel).expect("admits");
    assert_eq!(receipt.selected.len(), 1);
    sel.capabilities.clear();
    let err = expand_selection(&corpus(), &sel).expect_err("empty again");
    assert_eq!(err.rule(), "v1-empty-selection");
}

#[test]
fn capability_skips_are_receipted_when_a_selection_survives() {
    let mut sel = input(Stratum::Max, BuiltinProfile::Standard);
    // Add a standard-profile case that needs a missing capability.
    let mut cases = corpus();
    cases.push(case(
        "max/needs-tunnel",
        Stratum::Max,
        &[BuiltinProfile::Standard],
        &["wind-tunnel"],
    ));
    let receipt = expand_selection(&cases, &sel).expect("admits");
    assert!(
        receipt
            .selected
            .iter()
            .all(|c| c.as_str() != "max/needs-tunnel")
    );
    assert_eq!(receipt.skipped.len(), 1);
    assert_eq!(receipt.skipped[0].case.as_str(), "max/needs-tunnel");
    assert_eq!(receipt.skipped[0].reason, "capability:wind-tunnel");
    sel.capabilities.insert("wind-tunnel".to_owned());
    let with_cap = expand_selection(&cases, &sel).expect("admits");
    assert!(with_cap.skipped.is_empty());
    assert_ne!(receipt.digest(), with_cap.digest(), "skips are semantic");
}

#[test]
fn profile_flags_are_exactly_one_atomic_id() {
    assert!(matches!(
        ProfileId::parse_flags(&["standard"]),
        Ok(ProfileId::Builtin(BuiltinProfile::Standard))
    ));
    assert_eq!(
        ProfileId::parse_flags(&[]).expect_err("missing").rule(),
        "v1-profile-missing"
    );
    assert_eq!(
        ProfileId::parse_flags(&["smoke", "soak"])
            .expect_err("repeated flags")
            .rule(),
        "v1-profile-composition"
    );
    assert_eq!(
        ProfileId::parse_atomic("smoke+soak")
            .expect_err("implicit composition")
            .rule(),
        "v1-profile-composition"
    );
    let err = ProfileId::parse_atomic("FULL").expect_err("legacy tier");
    assert_eq!(err.rule(), "v1-legacy-alias");
    assert!(err.ranked_fixes()[0].contains("V.4.6"));
    assert_eq!(
        ProfileId::parse_atomic("MID").expect_err("legacy").rule(),
        "v1-legacy-alias"
    );
    assert_eq!(
        ProfileId::parse_atomic("warp-speed")
            .expect_err("unknown")
            .rule(),
        "v1-profile-unknown"
    );
}

#[test]
fn composite_profiles_are_frozen_ordered_and_duplicate_free() {
    let composite = CompositeProfile {
        id: "release-crossisa".to_owned(),
        version: 1,
        inputs: vec![BuiltinProfile::Release, BuiltinProfile::CrossIsa],
        precedence_rule: "release budgets win conflicts; cross-isa adds the second host".to_owned(),
    };
    composite.validate().expect("composite admits");
    let digest = composite.digest();
    let mut reordered = composite.clone();
    reordered.inputs.reverse();
    assert_ne!(
        digest,
        reordered.digest(),
        "input ORDER is precedence and participates in the frozen digest"
    );

    let mut duplicated = composite.clone();
    duplicated.inputs.push(BuiltinProfile::Release);
    assert_eq!(
        duplicated.validate().expect_err("dup input").rule(),
        "v1-profile-composition"
    );

    let receipt = expand_selection(
        &corpus(),
        &SelectionInput {
            stratum: Stratum::Max,
            profile: ProfileId::Composite(composite),
            filters: Vec::new(),
            capabilities: BTreeSet::from(["x86-64-host".to_owned()]),
            budgets: BTreeMap::new(),
            skips: Vec::new(),
            shards: 2,
        },
    )
    .expect("composite selection admits");
    assert!(
        receipt
            .selected
            .iter()
            .any(|c| c.as_str() == "core/gamma-isa"),
        "the cross-isa member admits its cases"
    );
}

#[test]
fn mutations_produce_exact_diffs_and_nothing_disappears_silently() {
    let baseline = expand_selection(&corpus(), &input(Stratum::Max, BuiltinProfile::Standard))
        .expect("baseline");

    // Core/Max swap: exact case delta, scope flagged.
    let core_only =
        expand_selection(&corpus(), &input(Stratum::Core, BuiltinProfile::Standard)).expect("core");
    let diff = semantic_diff(&baseline, &core_only);
    assert!(diff.scope_changed);
    assert_eq!(
        diff.removed.iter().map(CaseId::as_str).collect::<Vec<_>>(),
        ["max/delta"]
    );
    assert!(diff.added.is_empty());

    // Omit a case from the manifest: the diff names it exactly.
    let mut fewer = corpus();
    fewer.retain(|c| c.case.as_str() != "core/beta");
    let dropped =
        expand_selection(&fewer, &input(Stratum::Max, BuiltinProfile::Standard)).expect("admits");
    let diff = semantic_diff(&baseline, &dropped);
    assert_eq!(
        diff.removed.iter().map(CaseId::as_str).collect::<Vec<_>>(),
        ["core/beta"],
        "no case disappears silently"
    );
    assert_ne!(baseline.digest(), dropped.digest());

    // Budget change: exact row-level diff.
    let mut budgeted = input(Stratum::Max, BuiltinProfile::Standard);
    budgeted
        .budgets
        .insert("time-ns".to_owned(), 30_000_000_000);
    let rebudgeted = expand_selection(&corpus(), &budgeted).expect("admits");
    let diff = semantic_diff(&baseline, &rebudgeted);
    assert_eq!(
        diff.budget_changes,
        vec![(
            "time-ns".to_owned(),
            Some(60_000_000_000),
            Some(30_000_000_000)
        )]
    );
    assert!(!diff.scope_changed);

    // Identity: no mutation, empty diff, equal digests.
    let again =
        expand_selection(&corpus(), &input(Stratum::Max, BuiltinProfile::Standard)).expect("again");
    assert!(semantic_diff(&baseline, &again).is_empty());
    assert_eq!(baseline.digest(), again.digest());
}

#[test]
fn duplicate_case_ids_and_zero_shards_refuse() {
    let mut cases = corpus();
    cases.push(case(
        "core/alpha",
        Stratum::Core,
        &[BuiltinProfile::Smoke],
        &[],
    ));
    assert_eq!(
        expand_selection(&cases, &input(Stratum::Max, BuiltinProfile::Standard))
            .expect_err("duplicate id")
            .rule(),
        "v1-duplicate-case"
    );
    let mut sel = input(Stratum::Max, BuiltinProfile::Standard);
    sel.shards = 0;
    assert_eq!(
        expand_selection(&corpus(), &sel)
            .expect_err("zero shards")
            .rule(),
        "v1-selection-shards"
    );
}
