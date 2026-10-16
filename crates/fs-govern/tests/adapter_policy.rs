//! Battery for the external-adapter policy ratification (bead f85xj.11.1):
//! record completeness, the fail-closed accessor, and the three-way drift
//! gate that keeps the decision record, the AGENTS.md mission text, and the
//! xtask check-deps policy language in literal agreement.

use fs_govern::adapter_policy::{ADAPTER_POLICY_ID, adapter_policy, adapter_policy_json};

/// Whitespace/backtick normalization shared with the drift assertions, so
/// prose reflow cannot fake a policy change and backtick styling cannot hide
/// one.
fn normalized(text: &str) -> String {
    text.chars()
        .filter(|character| *character != '`')
        .collect::<String>()
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
}

#[track_caller]
fn assert_contains(haystack: &str, needle: &str, source: &str) {
    assert!(
        normalized(haystack).contains(&normalized(needle)),
        "policy drift: {source} no longer states {needle:?}"
    );
}

#[test]
fn record_is_complete_and_fail_closed_accessor_admits_it() {
    let record = adapter_policy().expect("the ratified record must validate");
    assert_eq!(record.id, ADAPTER_POLICY_ID);
    assert_eq!(record.chosen_option, "official-quarantined-adapters");
    assert!(
        record
            .options
            .iter()
            .any(|option| option.id == record.chosen_option),
        "chosen option must be among the considered options"
    );
    // Every rejected option carries a reason; a ruling with unreasoned
    // rejections is not a decision record.
    for option in record.options {
        if option.id != record.chosen_option {
            assert!(
                !option.rejection_reason.trim().is_empty(),
                "rejected option {} lacks a reason",
                option.id
            );
        }
    }
    assert!(
        record.options.len() >= 3,
        "all three bead options are recorded"
    );
    assert!(!record.ruling.is_empty());
    assert!(!record.invariants.is_empty());
    assert!(record.falsifiers.len() >= 3, "review triggers are named");
    for falsifier in record.falsifiers {
        assert!(
            falsifier.is_complete(),
            "falsifier {} incomplete",
            falsifier.id
        );
    }
    // The delegation that satisfied the bead's human-sign-off flag is
    // recorded verbatim, not implied.
    assert!(record.authority.contains("delegated"));
    assert!(record.authority.contains("2026-07-23"));
    // The downstream consumers the bead names are gated on this record.
    assert!(
        record
            .downstream_gates
            .contains(&"frankensim-extreal-program-f85xj.6.13")
    );
    assert!(
        record
            .downstream_gates
            .contains(&"frankensim-extreal-program-f85xj.11.5")
    );
}

#[test]
fn json_rendering_is_deterministic_and_carries_the_record() {
    let first = adapter_policy_json().expect("json renders");
    let second = adapter_policy_json().expect("json renders again");
    assert_eq!(first, second);
    for key in [
        "\"adapter_policy\":",
        "\"id\":\"ADPT-2026-07\"",
        "\"chosen_option\":\"official-quarantined-adapters\"",
        "\"options\":[",
        "\"ruling\":[",
        "\"invariants\":[",
        "\"falsifiers\":[",
        "\"downstream_gates\":[",
    ] {
        assert!(first.contains(key), "json lost {key:?}: {first}");
    }
}

#[test]
fn three_way_drift_gate_record_agents_md_and_xtask_agree() {
    let record = adapter_policy().expect("record");
    let agents_md = include_str!("../../../AGENTS.md");
    let xtask_source = include_str!("../../../xtask/src/main.rs");

    // AGENTS.md carries the amendment heading, the record id, and every
    // enforced amendment clause.
    assert_contains(agents_md, record.mission_amendment_heading, "AGENTS.md");
    assert_contains(agents_md, record.id, "AGENTS.md");
    for clause in record.mission_amendment_clauses {
        assert_contains(agents_md, clause, "AGENTS.md");
    }
    // The ruling AMENDS the mission text; it must not have weakened the base
    // Franken-only enumeration it sits beside.
    assert_contains(
        agents_md,
        "No BLAS, LAPACK, C, C++, Fortran, OpenCASCADE, gmsh, FEniCS, MFEM, OpenFOAM",
        "AGENTS.md",
    );
    assert_contains(
        agents_md,
        "check-deps enforcement is unchanged",
        "AGENTS.md",
    );

    // xtask's check-deps policy language cites the same ruling id and the
    // out-of-process boundary, so a policy-text edit in either place that
    // forgets the other fails here.
    assert_contains(xtask_source, record.id, "xtask/src/main.rs");
    assert_contains(
        xtask_source,
        "out-of-process quarantined adapter",
        "xtask/src/main.rs",
    );
    assert_contains(
        xtask_source,
        "never as a dependency edge",
        "xtask/src/main.rs",
    );

    // Controlled negative: the normalization cannot be satisfied vacuously.
    assert!(
        !normalized(agents_md).contains(&normalized("this clause was never ratified")),
        "drift helper must be able to fail"
    );
}
