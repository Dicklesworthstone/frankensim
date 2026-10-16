//! Canonical WFG problem identity conformance (7tv.24.10).
//!
//! This crate-external G0/G3 battery independently reconstructs the public
//! WFG problem identity preimage and mutation-locks every semantic field that
//! study and campaign provenance must bind.
//!
//! This identity names problem semantics only. It does not claim result or
//! campaign identity, optimizer quality, cancellation, cross-ISA stability,
//! external executable parity, or performance.

#![deny(unsafe_code)]

use fs_dfo::wfg::{
    WFG_ARITHMETIC_POLICY, WFG_EQUATION_REFERENCE, WFG_INPUT_DOMAIN_POLICY,
    WFG_PROBLEM_IDENTITY_KIND, WFG_PROBLEM_IDENTITY_SCHEMA_VERSION, WFG_TRACE_SCHEMA, WfgProblem,
    WfgVariant,
};
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};

const DECISION_EXCLUSION: (&str, &str) = (
    "decision-vector",
    "per-evaluation input belongs to the enclosing call identity",
);
const TRACE_EXCLUSION: (&str, &str) = (
    "evaluation-trace",
    "per-evaluation output belongs to the enclosing result identity",
);

#[derive(Debug, Clone, Copy)]
struct IdentityFields<'a> {
    kind: &'a str,
    schema_version: u64,
    variant: &'a str,
    objectives: u64,
    position_parameters: u64,
    distance_parameters: u64,
    equation_reference: &'a str,
    input_domain_policy: &'a str,
    arithmetic_policy: &'a str,
    trace_schema: &'a str,
}

fn independent_identity(fields: IdentityFields<'_>) -> ReplayIdentity {
    IdentityBuilder::new(fields.kind)
        .u64("problem-schema-version", fields.schema_version)
        .str("variant", fields.variant)
        .u64("objectives", fields.objectives)
        .u64("position-parameters", fields.position_parameters)
        .u64("distance-parameters", fields.distance_parameters)
        .str("equation-reference", fields.equation_reference)
        .str("input-domain-policy", fields.input_domain_policy)
        .str("arithmetic-policy", fields.arithmetic_policy)
        .str("trace-schema", fields.trace_schema)
        .exclude(DECISION_EXCLUSION.0, DECISION_EXCLUSION.1)
        .exclude(TRACE_EXCLUSION.0, TRACE_EXCLUSION.1)
        .finish()
}

fn base_fields() -> IdentityFields<'static> {
    IdentityFields {
        kind: "fs-dfo-wfg-problem-v1",
        schema_version: 1,
        variant: "WFG4",
        objectives: 4,
        position_parameters: 12,
        distance_parameters: 6,
        equation_reference: "jmetal-wfg-ea7e882f6b8f94b99535921674e62cda7986f20e",
        input_domain_policy: "normalized-[0,1]-or-canonical-z_i-[0,2(i+1)]-v1",
        arithmetic_policy: "fixed-left-to-right-reductions-fs-math-det-v1",
        trace_schema: "transformed/reduced/positioned/shape/objectives-normalized-equation-space-v1",
    }
}

fn problem_root(
    variant: WfgVariant,
    objectives: usize,
    position_parameters: usize,
    distance_parameters: usize,
) -> u64 {
    WfgProblem::new(
        variant,
        objectives,
        position_parameters,
        distance_parameters,
    )
    .unwrap()
    .replay_identity()
    .root()
}

#[test]
fn public_problem_identity_matches_the_independent_canonical_preimage() {
    assert_eq!(WFG_PROBLEM_IDENTITY_KIND, base_fields().kind);
    assert_eq!(WFG_PROBLEM_IDENTITY_SCHEMA_VERSION, 1);
    assert_eq!(WFG_EQUATION_REFERENCE, base_fields().equation_reference);
    assert_eq!(WFG_INPUT_DOMAIN_POLICY, base_fields().input_domain_policy);
    assert_eq!(WFG_ARITHMETIC_POLICY, base_fields().arithmetic_policy);
    assert_eq!(WFG_TRACE_SCHEMA, base_fields().trace_schema);

    let problem = WfgProblem::new(WfgVariant::Wfg4, 4, 12, 6).unwrap();
    let actual = problem.replay_identity();
    let expected = independent_identity(base_fields());
    assert_eq!(actual, expected);
    assert_eq!(actual.kind(), WFG_PROBLEM_IDENTITY_KIND);
    assert_eq!(actual.exclusions(), &[DECISION_EXCLUSION, TRACE_EXCLUSION]);
}

#[test]
fn every_canonical_identity_field_is_mutation_sensitive() {
    let base = base_fields();
    let base_root = independent_identity(base).root();
    let mutants = [
        IdentityFields {
            kind: "fs-dfo-wfg-problem-mutant",
            ..base
        },
        IdentityFields {
            schema_version: 2,
            ..base
        },
        IdentityFields {
            variant: "WFG5",
            ..base
        },
        IdentityFields {
            objectives: 5,
            ..base
        },
        IdentityFields {
            position_parameters: 6,
            ..base
        },
        IdentityFields {
            distance_parameters: 8,
            ..base
        },
        IdentityFields {
            equation_reference: "jmetal-wfg-mutant",
            ..base
        },
        IdentityFields {
            input_domain_policy: "normalized-only-mutant",
            ..base
        },
        IdentityFields {
            arithmetic_policy: "reassociated-mutant",
            ..base
        },
        IdentityFields {
            trace_schema: "objectives-only-mutant",
            ..base
        },
    ];

    for mutant in mutants {
        assert_ne!(independent_identity(mutant).root(), base_root, "{mutant:?}");
    }
}

#[test]
fn public_problem_mutations_move_roots_and_variants_never_alias() {
    let base_root = problem_root(WfgVariant::Wfg4, 4, 12, 6);
    for mutant_root in [
        problem_root(WfgVariant::Wfg5, 4, 12, 6),
        problem_root(WfgVariant::Wfg4, 5, 12, 6),
        problem_root(WfgVariant::Wfg4, 4, 6, 6),
        problem_root(WfgVariant::Wfg4, 4, 12, 8),
    ] {
        assert_ne!(mutant_root, base_root);
    }

    let roots = WfgVariant::ALL.map(|variant| problem_root(variant, 4, 12, 6));
    for (left_index, left) in roots.iter().enumerate() {
        for right in roots.iter().skip(left_index + 1) {
            assert_ne!(left, right);
        }
    }
}

#[test]
fn evaluation_does_not_mutate_problem_identity() {
    let problem = WfgProblem::new(WfgVariant::Wfg9, 4, 6, 6).unwrap();
    let normalized = [
        0.031_25, 0.093_75, 0.156_25, 0.218_75, 0.281_25, 0.343_75, 0.406_25, 0.468_75, 0.531_25,
        0.593_75, 0.718_75, 0.906_25,
    ];
    let mut canonical = [0.0; 12];
    for (component, (&value, output)) in normalized.iter().zip(&mut canonical).enumerate() {
        *output = value * 2.0 * (component as f64 + 1.0);
    }

    let before = problem.replay_identity();
    problem.evaluate_normalized(&normalized).unwrap();
    problem.evaluate_canonical(&canonical).unwrap();
    assert_eq!(problem.replay_identity(), before);
}
