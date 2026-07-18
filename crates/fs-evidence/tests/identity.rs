//! G0/G3/G4 coverage for typed color-evidence graph identity.
//!
//! A retained cross-ISA known-answer vector is still required for G5.

use fs_blake3::identity::{
    CancellationProbe, CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field,
    FieldSpec, NoClaimState, SchemaId, SourceId, StrongIdentity, TrustState, WireType,
};
use fs_evidence::{
    COLOR_ALGEBRA_VERSION, Color, ColorEvidenceCompositionOpV1, ColorEvidenceIdentityError,
    ColorEvidenceNodeIdV1, ColorEvidenceNodeIdentitySchemaV1, ColorEvidenceNodeKindV1,
    ColorEvidenceNodeV1, ColorEvidenceOperationV1, ColorEvidenceParentSemanticsV1,
    ColorEvidenceSourceIdV1, ColorEvidenceSourceV1, ValidityDomain,
    compose_color_evidence_nodes_v1, identify_color_evidence_source_node_v1,
    identify_color_evidence_source_v1,
};

const LIMITS: CanonicalLimits = CanonicalLimits::new(16_384, 8_192, 32, 64, 256);

enum ForeignSourceSchema {}

impl CanonicalSchema for ForeignSourceSchema {
    const DOMAIN: &'static str = "org.frankensim.tests.foreign-color-source.v1";
    const NAME: &'static str = "foreign-color-source";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "negative child-binding fixture";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("value", WireType::U64)];
}

fn verified(lo: f64, hi: f64) -> Color {
    Color::Verified { lo, hi }
}

fn source_receipt(label: &str) -> ColorEvidenceSourceV1 {
    identify_color_evidence_source_v1(
        "org.frankensim.tests.color-source",
        1,
        label.as_bytes(),
        LIMITS,
        || false,
    )
    .expect("valid retained source")
}

fn source_node(label: &str, output: Color) -> ColorEvidenceNodeV1 {
    identify_color_evidence_source_node_v1(&source_receipt(label), output, LIMITS, || false)
        .expect("valid source node")
}

fn derived(
    operation: ColorEvidenceCompositionOpV1,
    left: &ColorEvidenceNodeV1,
    right: &ColorEvidenceNodeV1,
) -> Result<ColorEvidenceNodeV1, ColorEvidenceIdentityError> {
    compose_color_evidence_nodes_v1(operation, left, right, LIMITS, || false)
}

fn parent_reference_bytes(parent: ColorEvidenceNodeIdV1) -> [u8; 65] {
    let mut output = [0_u8; 65];
    output[0] = ColorEvidenceNodeIdV1::ROLE.tag();
    output[1..33]
        .copy_from_slice(SchemaId::<ColorEvidenceNodeIdentitySchemaV1>::for_schema().as_bytes());
    output[33..].copy_from_slice(parent.as_bytes());
    output
}

#[test]
fn typed_receipts_are_replay_stable_and_explicitly_unanchored() {
    let first = source_node("observation-a", verified(-1.0, 2.0));
    let replay = source_node("observation-a", verified(-1.0, 2.0));

    assert_eq!(first.id(), replay.id());
    assert_eq!(
        first.receipt().canonical_preimage(),
        replay.receipt().canonical_preimage()
    );

    let audit = first.receipt().audit_record();
    assert_eq!(audit.trust(), TrustState::Unanchored);
    assert_eq!(audit.no_claim(), NoClaimState::ExternalTrustRequired);
    assert_eq!(audit.id(), first.id_bytes());
    assert_eq!(first.trust_state(), TrustState::Unanchored);
}

#[test]
fn every_available_semantic_mutation_moves_the_typed_root() {
    let source_a = source_node("observation-a", verified(0.0, 1.0));
    let source_b = source_node("observation-b", verified(1.0, 2.0));
    let source_payload_changed = source_node("observation-a+", verified(0.0, 1.0));
    let output_changed = source_node("observation-a", verified(0.0, 1.5));

    assert_ne!(source_a.id(), source_b.id());
    assert_ne!(source_a.id(), source_payload_changed.id());
    assert_ne!(source_a.id(), output_changed.id());

    let add = derived(ColorEvidenceCompositionOpV1::Add, &source_a, &source_b).expect("valid add");
    let mul =
        derived(ColorEvidenceCompositionOpV1::Mul, &source_a, &source_b).expect("valid multiply");
    assert_ne!(add.id(), mul.id());
}

#[test]
fn commutative_multisets_normalize_order_but_preserve_multiplicity() {
    let a = source_node("a", verified(0.0, 1.0));
    let b = source_node("b", verified(1.0, 2.0));

    for operation in [
        ColorEvidenceCompositionOpV1::Add,
        ColorEvidenceCompositionOpV1::Mul,
        ColorEvidenceCompositionOpV1::Hull,
    ] {
        let ab = derived(operation, &a, &b).expect("a op b");
        let ba = derived(operation, &b, &a).expect("b op a");
        assert_eq!(ab.id(), ba.id());
        assert_eq!(ab.color(), ba.color());
    }

    let ab = derived(ColorEvidenceCompositionOpV1::Add, &a, &b).expect("a+b");
    let aa = derived(ColorEvidenceCompositionOpV1::Add, &a, &a).expect("a+a");
    assert_ne!(ab.id(), aa.id(), "a+a must not collapse into a+b");
}

#[test]
fn composition_recomputes_the_current_color_algebra() {
    let a = source_node("a", verified(-1.0, 2.0));
    let b = source_node("b", verified(3.0, 4.0));

    let add = derived(ColorEvidenceCompositionOpV1::Add, &a, &b).expect("add");
    let mul = derived(ColorEvidenceCompositionOpV1::Mul, &a, &b).expect("multiply");
    let hull = derived(ColorEvidenceCompositionOpV1::Hull, &a, &b).expect("hull");

    assert_eq!(
        add.color(),
        &verified(2.0_f64.next_down(), 6.0_f64.next_up())
    );
    assert_eq!(
        mul.color(),
        &verified((-4.0_f64).next_down(), 8.0_f64.next_up())
    );
    assert_eq!(hull.color(), &verified(-1.0, 4.0));
    assert_eq!(add.operation(), ColorEvidenceOperationV1::Add);
    assert_eq!(add.kind(), ColorEvidenceNodeKindV1::Composition);
    assert_eq!(
        add.parent_semantics(),
        ColorEvidenceParentSemanticsV1::CommutativeMultiset
    );
    assert_ne!(add.id(), mul.id());
    assert_ne!(mul.id(), hull.id());
}

#[test]
fn source_domain_version_and_payload_are_all_identity_bearing() {
    let a = identify_color_evidence_source_v1("domain-a", 1, b"same", LIMITS, || false)
        .expect("domain a v1");
    let b = identify_color_evidence_source_v1("domain-b", 1, b"same", LIMITS, || false)
        .expect("domain b v1");
    let c = identify_color_evidence_source_v1("domain-a", 2, b"same", LIMITS, || false)
        .expect("domain a v2");
    let d = identify_color_evidence_source_v1("domain-a", 1, b"different", LIMITS, || false)
        .expect("different bytes");
    assert_ne!(a.id(), b.id());
    assert_ne!(a.id(), c.id());
    assert_ne!(a.id(), d.id());
    assert_eq!(a.trust_state(), TrustState::Unanchored);

    assert_eq!(
        identify_color_evidence_source_v1("", 1, b"x", LIMITS, || false),
        Err(ColorEvidenceIdentityError::EmptySourceDomain)
    );
    assert_eq!(
        identify_color_evidence_source_v1("domain-a", 0, b"x", LIMITS, || false),
        Err(ColorEvidenceIdentityError::ZeroSourceSchemaVersion)
    );
}

#[test]
fn malformed_output_and_cancellation_publish_no_identity() {
    let source = source_receipt("source");
    let malformed =
        identify_color_evidence_source_node_v1(&source, verified(2.0, 1.0), LIMITS, || false);
    assert!(matches!(
        malformed,
        Err(ColorEvidenceIdentityError::MalformedColor(_))
    ));

    #[derive(Debug)]
    struct CancelNow;
    impl CancellationProbe for CancelNow {
        fn is_cancelled(&mut self) -> bool {
            true
        }
    }

    let cancelled =
        identify_color_evidence_source_node_v1(&source, verified(0.0, 1.0), LIMITS, CancelNow);
    assert!(matches!(
        cancelled,
        Err(ColorEvidenceIdentityError::Canonical(
            CanonicalError::Cancelled { .. }
        ))
    ));

    let left = source_node("cancel-left", verified(0.0, 1.0));
    let right = source_node("cancel-right", verified(1.0, 2.0));
    let cancelled_composition = compose_color_evidence_nodes_v1(
        ColorEvidenceCompositionOpV1::Add,
        &left,
        &right,
        LIMITS,
        CancelNow,
    );
    assert!(matches!(
        cancelled_composition,
        Err(ColorEvidenceIdentityError::Canonical(
            CanonicalError::Cancelled { .. }
        ))
    ));

    #[derive(Debug)]
    struct CancelAfter {
        successful_polls: usize,
    }
    impl CancellationProbe for CancelAfter {
        fn is_cancelled(&mut self) -> bool {
            if self.successful_polls == 0 {
                true
            } else {
                self.successful_polls -= 1;
                false
            }
        }
    }
    let midstream_regime = ValidityDomain::unconstrained()
        .with("axis-c", 0.0, 1.0)
        .with("axis-a", 0.0, 1.0)
        .with("axis-b", 0.0, 1.0);
    let midstream_cancelled = identify_color_evidence_source_node_v1(
        &source,
        Color::Validated {
            regime: midstream_regime,
            dataset: "dataset-v1".to_string(),
        },
        LIMITS,
        CancelAfter {
            successful_polls: 2,
        },
    );
    assert!(matches!(
        midstream_cancelled,
        Err(ColorEvidenceIdentityError::Canonical(
            CanonicalError::Cancelled { .. }
        ))
    ));

    let mut regime = ValidityDomain::unconstrained();
    for index in 0..8 {
        regime = regime.with(format!("axis-{index}"), 0.0, 1.0);
    }
    let tiny = CanonicalLimits::new(512, 256, 32, 2, 64);
    let bounded = identify_color_evidence_source_node_v1(
        &source,
        Color::Validated {
            regime,
            dataset: "dataset-v1".to_string(),
        },
        tiny,
        || false,
    );
    assert!(matches!(
        bounded,
        Err(ColorEvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded { .. }
        ))
    ));

    let field_tiny = CanonicalLimits::new(4096, 33, 32, 64, 64);
    let field_bounded =
        identify_color_evidence_source_node_v1(&source, verified(0.0, 1.0), field_tiny, || false);
    assert!(matches!(
        field_bounded,
        Err(ColorEvidenceIdentityError::Canonical(
            CanonicalError::LimitExceeded {
                kind: fs_blake3::identity::LimitKind::FieldBytes,
                requested: 34,
                limit: 33,
            }
        ))
    ));
}

#[test]
fn source_child_role_and_schema_are_bound_by_the_node_schema() {
    let foreign = CanonicalEncoder::<SourceId<ForeignSourceSchema>, _>::new(LIMITS, || false)
        .expect("foreign source schema")
        .u64(Field::new(0, "value"), 7)
        .expect("foreign source value")
        .finish()
        .expect("foreign source receipt");

    let refusal = CanonicalEncoder::<ColorEvidenceNodeIdV1, _>::new(LIMITS, || false)
        .expect("node schema")
        .variant(Field::new(0, "node-kind"), 1, &[])
        .expect("node kind")
        .variant(Field::new(1, "operation"), 1, &[])
        .expect("operation")
        .variant(Field::new(2, "parent-semantics"), 1, &[])
        .expect("parent semantics")
        .u64(
            Field::new(3, "color-algebra-version"),
            u64::from(COLOR_ALGEBRA_VERSION),
        )
        .expect("algebra version")
        .ordered_children(Field::new(4, "source"), 1, [foreign.id()])
        .expect_err("foreign source role/schema must refuse");
    assert!(matches!(
        refusal,
        CanonicalError::ChildBindingMismatch {
            field: "source",
            ..
        }
    ));
}

#[test]
fn source_node_helper_matches_every_color_canonical_encoding() {
    let source = source_receipt("manual-parity");
    let validated_regime = ValidityDomain::unconstrained()
        .with("z-axis", -3.0, 4.0)
        .with("a-axis", 0.0, 1.0);
    let colors = [
        verified(-2.0, 3.0),
        verified(f64::NEG_INFINITY, f64::INFINITY),
        Color::Validated {
            regime: validated_regime,
            dataset: "dataset-v1".to_string(),
        },
        Color::Estimated {
            estimator: "estimator-v1".to_string(),
            dispersion: -0.0,
        },
        Color::Estimated {
            estimator: "unbounded-estimator-v1".to_string(),
            dispersion: f64::INFINITY,
        },
    ];

    for color in colors {
        let helper =
            identify_color_evidence_source_node_v1(&source, color.clone(), LIMITS, || false)
                .expect("helper source node");
        let color_bytes = color.canonical_bytes();
        let manual = CanonicalEncoder::<ColorEvidenceNodeIdV1, _>::new(LIMITS, || false)
            .expect("node schema")
            .variant(Field::new(0, "node-kind"), 1, &[])
            .expect("node kind")
            .variant(Field::new(1, "operation"), 1, &[])
            .expect("operation")
            .variant(Field::new(2, "parent-semantics"), 1, &[])
            .expect("parent semantics")
            .u64(
                Field::new(3, "color-algebra-version"),
                u64::from(COLOR_ALGEBRA_VERSION),
            )
            .expect("algebra version")
            .ordered_children(Field::new(4, "source"), 1, [source.id()])
            .expect("typed source child")
            .bytes(Field::new(5, "output-color"), &color_bytes)
            .expect("canonical color")
            .ordered_bytes(Field::new(6, "parents"), 0, core::iter::empty::<&[u8]>())
            .expect("empty parent list")
            .finish()
            .expect("manual source node");
        assert_eq!(helper.id(), manual.id());
        assert_eq!(
            helper.receipt().canonical_preimage(),
            manual.canonical_preimage()
        );
    }
}

#[test]
fn composition_helper_matches_descriptor_bound_parent_construction() {
    let a = source_node("composition-parity-a", verified(-1.0, 2.0));
    let b = source_node("composition-parity-b", verified(3.0, 4.0));
    let helper = derived(ColorEvidenceCompositionOpV1::Add, &a, &b).expect("helper add");
    let mut parents = [a.id(), b.id()];
    parents.sort_unstable();
    let parent_rows = parents.map(parent_reference_bytes);
    let output_bytes = verified(2.0_f64.next_down(), 6.0_f64.next_up()).canonical_bytes();

    let manual = CanonicalEncoder::<ColorEvidenceNodeIdV1, _>::new(LIMITS, || false)
        .expect("node schema")
        .variant(Field::new(0, "node-kind"), 2, &[])
        .expect("composition kind")
        .variant(Field::new(1, "operation"), 2, &[])
        .expect("add operation")
        .variant(Field::new(2, "parent-semantics"), 2, &[])
        .expect("commutative multiset")
        .u64(
            Field::new(3, "color-algebra-version"),
            u64::from(COLOR_ALGEBRA_VERSION),
        )
        .expect("algebra version")
        .ordered_children(
            Field::new(4, "source"),
            0,
            core::iter::empty::<ColorEvidenceSourceIdV1>(),
        )
        .expect("no source child")
        .bytes(Field::new(5, "output-color"), &output_bytes)
        .expect("canonical color")
        .ordered_bytes(
            Field::new(6, "parents"),
            2,
            parent_rows.iter().map(|row| row.as_slice()),
        )
        .expect("descriptor-bound parents")
        .finish()
        .expect("manual composition node");

    assert_eq!(helper.id(), manual.id());
    assert_eq!(
        helper.receipt().canonical_preimage(),
        manual.canonical_preimage()
    );

    let mut wrong_role_rows = parent_rows;
    wrong_role_rows[0][0] ^= 0xff;
    let wrong_role = CanonicalEncoder::<ColorEvidenceNodeIdV1, _>::new(LIMITS, || false)
        .expect("node schema")
        .variant(Field::new(0, "node-kind"), 2, &[])
        .expect("composition kind")
        .variant(Field::new(1, "operation"), 2, &[])
        .expect("add operation")
        .variant(Field::new(2, "parent-semantics"), 2, &[])
        .expect("commutative multiset")
        .u64(
            Field::new(3, "color-algebra-version"),
            u64::from(COLOR_ALGEBRA_VERSION),
        )
        .expect("algebra version")
        .ordered_children(
            Field::new(4, "source"),
            0,
            core::iter::empty::<ColorEvidenceSourceIdV1>(),
        )
        .expect("no source child")
        .bytes(Field::new(5, "output-color"), &output_bytes)
        .expect("canonical color")
        .ordered_bytes(
            Field::new(6, "parents"),
            2,
            wrong_role_rows.iter().map(|row| row.as_slice()),
        )
        .expect("schema-shaped but semantically foreign parent rows")
        .finish()
        .expect("foreign manual node");
    assert_ne!(helper.id(), wrong_role.id());
    assert_ne!(manual.id(), wrong_role.id());
}
