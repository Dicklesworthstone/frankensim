//! G0/G3/G5 conformance for bounded hierarchical, nonlinear, and
//! mode-switching tolerance propagation.

use std::num::NonZeroU64;

use fs_toleralloc::{
    AdmittedCorrelationModel, ClampDisposition, ColorRank, CorrelatedStackTerm, InteriorKnotOwner,
    MAX_EXACT_STRUCTURED_WEIGHT_V1, PiecewiseQuadraticLaw, QuadraticResponsePiece,
    StructuredEvaluationStage, StructuredKeyRole, StructuredLawId, StructuredLawIssue,
    StructuredModelIdentity, StructuredMomentQuantity, StructuredMomentScope, StructuredNodeId,
    StructuredNodeSpec, StructuredNumericIssue, StructuredPopulationModel,
    StructuredPropagationError, StructuredResource, StructuredTopologyIssue,
    propagate_correlated_stack, propagate_structured_population,
};

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("fixture weights are positive")
}

fn branch(key: &str, parent: Option<u32>) -> StructuredNodeSpec {
    StructuredNodeSpec::Branch {
        key: key.into(),
        parent: parent.map(StructuredNodeId),
    }
}

fn leaf(key: &str, parent: u32, weight: u64, raw_clearance: f64) -> StructuredNodeSpec {
    StructuredNodeSpec::Leaf {
        key: key.into(),
        parent: StructuredNodeId(parent),
        relative_weight: nz(weight),
        raw_clearance,
        law: StructuredLawId(0),
    }
}

fn backlash_law() -> PiecewiseQuadraticLaw {
    PiecewiseQuadraticLaw {
        key: "backlash-clearance".into(),
        lower_bound: -3.0,
        upper_bound: 3.0,
        knots: vec![-3.0, -1.0, 1.0, 3.0],
        interior_knot_owners: vec![InteriorKnotOwner::UpperPiece, InteriorKnotOwner::LowerPiece],
        pieces: vec![
            QuadraticResponsePiece {
                mode_key: "drive-negative".into(),
                // -(x + 1)^2
                a: -1.0,
                b: -2.0,
                c: -1.0,
            },
            QuadraticResponsePiece {
                mode_key: "lash".into(),
                a: 0.0,
                b: 0.0,
                c: 0.0,
            },
            QuadraticResponsePiece {
                mode_key: "drive-positive".into(),
                // (x - 1)^2
                a: 1.0,
                b: -2.0,
                c: 1.0,
            },
        ],
    }
}

fn gear_model() -> StructuredPopulationModel {
    StructuredPopulationModel {
        identity: StructuredModelIdentity::try_new("gear/backlash-population", nz(1), [0x6d; 32])
            .expect("fixture identity is canonical"),
        nodes: vec![
            branch("gear-population", None),
            branch("process-low", Some(0)),
            branch("lot-low-far", Some(1)),
            leaf("part-low-far-extreme", 2, 1, -4.0),
            leaf("part-low-far-near", 2, 3, -2.0),
            branch("lot-low-lash", Some(1)),
            leaf("part-low-lash-boundary", 5, 3, -1.0),
            leaf("part-low-lash-center", 5, 9, 0.0),
            branch("process-high", Some(0)),
            branch("lot-high-lash", Some(8)),
            leaf("part-high-lash-center", 9, 9, 0.0),
            leaf("part-high-lash-boundary", 9, 27, 1.0),
            branch("lot-high-far", Some(8)),
            leaf("part-high-far-near", 12, 9, 2.0),
            leaf("part-high-far-extreme", 12, 3, 4.0),
        ],
        laws: vec![backlash_law()],
    }
}

fn assert_close(actual: f64, expected: f64) {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    assert!(
        (actual - expected).abs() <= 64.0 * f64::EPSILON * scale,
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

#[test]
fn structured_backlash_matches_the_integer_population_oracle() {
    let model = gear_model();
    let receipt = propagate_structured_population(&model).expect("structured fixture evaluates");

    // Exhaustive integer oracle over the eight compressed observations:
    // sum(w)=64, sum(w*y)=14, sum(w*y^2)=76.
    let expected_mean = 7.0 / 32.0;
    let expected_variance = 1167.0 / 1024.0;
    assert_eq!(receipt.schema_version(), 1);
    assert_eq!(receipt.total_weight(), 64);
    assert_close(receipt.mean(), expected_mean);
    assert_close(receipt.variance(), expected_variance);
    assert_close(receipt.standard_deviation(), expected_variance.sqrt());
    assert_close(receipt.direct_mean(), expected_mean);
    assert_close(receipt.direct_variance(), expected_variance);
    assert_close(receipt.hierarchy_mean_residual(), 0.0);
    assert_close(receipt.hierarchy_variance_residual(), 0.0);
    assert_eq!(receipt.model(), &model);

    let root = &receipt.nodes()[0];
    assert_eq!(root.node(), StructuredNodeId(0));
    assert_eq!(root.relative_weight(), 64);
    assert_eq!(root.descendant_leaf_count(), 8);
    assert_close(root.mean(), expected_mean);
    assert_close(root.within_child_variance(), 255.0 / 256.0);
    assert_close(root.between_child_variance(), 147.0 / 1024.0);
    assert_close(root.total_variance(), expected_variance);
    assert_close(root.decomposition_residual(), 0.0);

    for process_index in [1_usize, 8] {
        let process = &receipt.nodes()[process_index];
        assert_close(process.within_child_variance(), 27.0 / 64.0);
        assert_close(process.between_child_variance(), 147.0 / 256.0);
        assert_close(process.total_variance(), 255.0 / 256.0);
    }
    assert_close(receipt.nodes()[1].mean(), -7.0 / 16.0);
    assert_close(receipt.nodes()[8].mean(), 7.0 / 16.0);

    let modes = receipt.modes();
    assert_eq!(modes.len(), 3);
    assert_eq!(modes[0].relative_weight(), 4);
    assert_eq!(modes[0].leaf_count(), 2);
    assert_close(modes[0].mean().expect("observed"), -7.0 / 4.0);
    assert_close(modes[0].variance().expect("observed"), 27.0 / 16.0);
    assert_eq!(modes[1].relative_weight(), 48);
    assert_eq!(modes[1].leaf_count(), 4);
    assert_close(modes[1].mean().expect("observed"), 0.0);
    assert_close(modes[1].variance().expect("observed"), 0.0);
    assert_eq!(modes[2].relative_weight(), 12);
    assert_eq!(modes[2].leaf_count(), 2);
    assert_close(modes[2].mean().expect("observed"), 7.0 / 4.0);
    assert_close(modes[2].variance().expect("observed"), 27.0 / 16.0);
    assert_close(receipt.within_mode_variance(), 27.0 / 64.0);
    assert_close(receipt.between_mode_variance(), 735.0 / 1024.0);
    assert_close(receipt.mode_decomposed_variance(), expected_variance);
    assert_close(receipt.mode_decomposition_residual(), 0.0);
}

#[test]
fn backlash_modes_are_total_disjoint_boundary_owned_and_saturated() {
    let receipt = propagate_structured_population(&gear_model()).expect("fixture evaluates");
    let leaves = receipt.leaves();
    assert_eq!(leaves.len(), 8);

    assert_eq!(leaves[0].raw_clearance(), -4.0);
    assert_eq!(leaves[0].clamped_clearance(), -3.0);
    assert_eq!(
        leaves[0].clamp_disposition(),
        ClampDisposition::LowerClamped
    );
    assert_eq!(leaves[0].selected_piece(), 0);
    assert_eq!(leaves[0].output(), -4.0);

    assert_eq!(leaves[2].raw_clearance(), -1.0);
    assert_eq!(leaves[2].clamp_disposition(), ClampDisposition::Unchanged);
    assert_eq!(leaves[2].selected_piece(), 1);
    assert_eq!(leaves[2].output().to_bits(), 0);
    assert_eq!(leaves[5].raw_clearance(), 1.0);
    assert_eq!(leaves[5].selected_piece(), 1);
    assert_eq!(leaves[5].output().to_bits(), 0);

    assert_eq!(leaves[7].raw_clearance(), 4.0);
    assert_eq!(leaves[7].clamped_clearance(), 3.0);
    assert_eq!(
        leaves[7].clamp_disposition(),
        ClampDisposition::UpperClamped
    );
    assert_eq!(leaves[7].selected_piece(), 2);
    assert_eq!(leaves[7].output(), 4.0);
}

#[test]
fn nominal_deadband_linearization_honestly_misses_mode_switching() {
    let identity = AdmittedCorrelationModel::try_new(
        "gear/nominal-deadband",
        nz(1),
        [0x31; 32],
        2,
        vec![1.0, 0.0, 0.0, 1.0],
    )
    .expect("identity correlation factor is admitted");
    let zero_nominal_gradient = [
        CorrelatedStackTerm {
            name: "process-offset".into(),
            signed_sensitivity: 0.0,
            sensitivity_color: ColorRank::Verified,
            standard_deviation: 1.0,
        },
        CorrelatedStackTerm {
            name: "part-offset".into(),
            signed_sensitivity: 0.0,
            sensitivity_color: ColorRank::Verified,
            standard_deviation: 1.0,
        },
    ];
    let local = propagate_correlated_stack(&identity, &zero_nominal_gradient)
        .expect("an exact zero local gradient remains explicit");
    assert_eq!(local.independent_variance().to_bits(), 0);
    assert_eq!(local.correlated_variance().to_bits(), 0);

    let structured =
        propagate_structured_population(&gear_model()).expect("finite population evaluates");
    assert_close(structured.mean(), 7.0 / 32.0);
    assert_close(structured.variance(), 1167.0 / 1024.0);
    assert!(structured.variance() > local.correlated_variance());
}

#[test]
fn g3_weight_replication_and_affine_response_preserve_central_semantics() {
    let base_model = gear_model();
    let base = propagate_structured_population(&base_model).expect("base evaluates");

    let mut replicated_model = base_model.clone();
    for node in &mut replicated_model.nodes {
        if let StructuredNodeSpec::Leaf {
            relative_weight, ..
        } = node
        {
            *relative_weight = nz(relative_weight.get() * 3);
        }
    }
    let replicated =
        propagate_structured_population(&replicated_model).expect("replicated weights evaluate");
    assert_eq!(replicated.total_weight(), base.total_weight() * 3);
    assert_close(replicated.mean(), base.mean());
    assert_close(replicated.variance(), base.variance());
    assert_close(
        replicated.within_mode_variance(),
        base.within_mode_variance(),
    );
    assert_close(
        replicated.between_mode_variance(),
        base.between_mode_variance(),
    );

    let mut affine_model = base_model;
    for piece in &mut affine_model.laws[0].pieces {
        piece.a *= 2.0;
        piece.b *= 2.0;
        piece.c = piece.c * 2.0 + 5.0;
    }
    let affine = propagate_structured_population(&affine_model).expect("affine response evaluates");
    assert_close(affine.mean(), 2.0 * base.mean() + 5.0);
    assert_close(affine.variance(), 4.0 * base.variance());
    assert_close(
        affine.within_mode_variance(),
        4.0 * base.within_mode_variance(),
    );
    assert_close(
        affine.between_mode_variance(),
        4.0 * base.between_mode_variance(),
    );
    assert_eq!(
        affine
            .leaves()
            .iter()
            .map(|leaf| leaf.selected_piece())
            .collect::<Vec<_>>(),
        base.leaves()
            .iter()
            .map(|leaf| leaf.selected_piece())
            .collect::<Vec<_>>()
    );
}

#[test]
fn g5_identical_invocation_replays_the_complete_receipt() {
    let model = gear_model();
    let first = propagate_structured_population(&model).expect("first replay evaluates");
    let second = propagate_structured_population(&model).expect("second replay evaluates");
    assert_eq!(first, second);
    assert_eq!(
        first.model().identity.namespace(),
        "gear/backlash-population"
    );
    assert_eq!(first.model().identity.schema_version().get(), 1);
    assert_eq!(first.model().identity.semantic_digest(), [0x6d; 32]);
    assert_eq!(first.model().nodes, model.nodes);
    assert_eq!(first.model().laws, model.laws);
}

#[test]
fn structured_model_admission_is_bounded_and_topology_safe() {
    let mut model = gear_model();
    model.nodes.clear();
    assert_eq!(
        propagate_structured_population(&model),
        Err(StructuredPropagationError::ResourceLimit {
            resource: StructuredResource::Nodes,
            actual: 0,
            max: 8192,
        })
    );

    let mut model = gear_model();
    model.nodes[0] = leaf("not-a-root", 0, 1, 0.0);
    assert_eq!(
        propagate_structured_population(&model),
        Err(StructuredPropagationError::InvalidTopology {
            node: StructuredNodeId(0),
            related: None,
            issue: StructuredTopologyIssue::RootMustBeBranch,
        })
    );

    let mut model = gear_model();
    model.nodes[1] = branch("process-low", Some(8));
    assert_eq!(
        propagate_structured_population(&model),
        Err(StructuredPropagationError::InvalidTopology {
            node: StructuredNodeId(1),
            related: Some(StructuredNodeId(8)),
            issue: StructuredTopologyIssue::ParentNotBeforeChild,
        })
    );

    let mut model = gear_model();
    model.nodes[5] = branch("lot-low-lash", Some(3));
    assert_eq!(
        propagate_structured_population(&model),
        Err(StructuredPropagationError::InvalidTopology {
            node: StructuredNodeId(5),
            related: Some(StructuredNodeId(3)),
            issue: StructuredTopologyIssue::ParentIsLeaf,
        })
    );

    let mut model = gear_model();
    model.nodes.push(branch("unused-branch", Some(0)));
    assert!(matches!(
        propagate_structured_population(&model),
        Err(StructuredPropagationError::InvalidTopology {
            issue: StructuredTopologyIssue::EmptyBranch,
            ..
        })
    ));

    let mut model = gear_model();
    if let StructuredNodeSpec::Leaf {
        relative_weight, ..
    } = &mut model.nodes[3]
    {
        *relative_weight = nz(MAX_EXACT_STRUCTURED_WEIGHT_V1 + 1);
    }
    assert!(matches!(
        propagate_structured_population(&model),
        Err(StructuredPropagationError::ResourceLimit {
            resource: StructuredResource::TotalWeight,
            ..
        })
    ));
}

#[test]
fn structured_laws_refuse_unstable_keys_and_ambiguous_layouts() {
    let mut model = gear_model();
    model.laws[0].pieces[0].mode_key = "Drive-Negative".into();
    assert!(matches!(
        propagate_structured_population(&model),
        Err(StructuredPropagationError::InvalidKey {
            role: StructuredKeyRole::Mode,
            owner_index: Some(0),
            index: 0,
            ..
        })
    ));

    let mut model = gear_model();
    model.laws[0].interior_knot_owners.clear();
    assert_eq!(
        propagate_structured_population(&model),
        Err(StructuredPropagationError::InvalidLawLayout {
            law: StructuredLawId(0),
            issue: StructuredLawIssue::InteriorOwnerCountMismatch,
            knot: None,
        })
    );

    let mut model = gear_model();
    model.laws[0].knots[2] = -1.0;
    assert_eq!(
        propagate_structured_population(&model),
        Err(StructuredPropagationError::InvalidLawLayout {
            law: StructuredLawId(0),
            issue: StructuredLawIssue::KnotsNotIncreasing,
            knot: Some(2),
        })
    );

    let mut model = gear_model();
    model.laws[0].pieces[2].mode_key = "lash".into();
    assert!(matches!(
        propagate_structured_population(&model),
        Err(StructuredPropagationError::DuplicateKey {
            role: StructuredKeyRole::Mode,
            owner_index: Some(0),
            first_index: 1,
            duplicate_index: 2,
            ..
        })
    ));
}

#[test]
fn structured_evaluation_refuses_overflow_and_underflow_without_a_partial_receipt() {
    let identity = StructuredModelIdentity::try_new("gear/numeric-boundary", nz(1), [0x44; 32])
        .expect("identity is admitted");
    let mut overflow = StructuredPopulationModel {
        identity: identity.clone(),
        nodes: vec![branch("root", None), leaf("extreme", 0, 1, f64::MAX)],
        laws: vec![PiecewiseQuadraticLaw {
            key: "overflow-law".into(),
            lower_bound: 0.0,
            upper_bound: f64::MAX,
            knots: vec![0.0, f64::MAX],
            interior_knot_owners: Vec::new(),
            pieces: vec![QuadraticResponsePiece {
                mode_key: "active".into(),
                a: 1.0,
                b: 0.0,
                c: 0.0,
            }],
        }],
    };
    assert_eq!(
        propagate_structured_population(&overflow),
        Err(StructuredPropagationError::InvalidEvaluation {
            node: StructuredNodeId(1),
            law: StructuredLawId(0),
            piece: 0,
            stage: StructuredEvaluationStage::LinearProduct,
            issue: StructuredNumericIssue::NonFinite,
        })
    );

    overflow.nodes[1] = leaf("extreme", 0, 1, f64::MIN_POSITIVE);
    overflow.laws[0].upper_bound = 1.0;
    overflow.laws[0].knots[1] = 1.0;
    overflow.laws[0].pieces[0].a = f64::MIN_POSITIVE;
    assert_eq!(
        propagate_structured_population(&overflow),
        Err(StructuredPropagationError::InvalidEvaluation {
            node: StructuredNodeId(1),
            law: StructuredLawId(0),
            piece: 0,
            stage: StructuredEvaluationStage::QuadraticProduct,
            issue: StructuredNumericIssue::Underflow,
        })
    );
}

#[test]
fn structured_normalization_refuses_distinct_outputs_that_alias_to_one_subnormal() {
    // Dividing by 2^512 maps 7*2^-564 and 9*2^-564 to 1.75 and 2.25
    // minimum-subnormal units. Both round to the same nonzero binary64 value,
    // even though the retained response outputs are distinct.
    let scale = f64::from_bits(1_535_u64 << 52);
    let unit = f64::from_bits(459_u64 << 52);
    let first = 7.0 * unit;
    let second = 9.0 * unit;
    let first_normalized = first / scale;
    let second_normalized = second / scale;
    assert_ne!(first.to_bits(), second.to_bits());
    assert_eq!(first_normalized.to_bits(), 2);
    assert_eq!(second_normalized.to_bits(), 2);

    let model = StructuredPopulationModel {
        identity: StructuredModelIdentity::try_new(
            "gear/normalized-subnormal-alias",
            nz(1),
            [0x58; 32],
        )
        .expect("fixture identity is canonical"),
        nodes: vec![
            branch("root", None),
            leaf("first-tiny", 0, 1, first),
            leaf("second-tiny", 0, 1, second),
            leaf("scale-anchor", 0, 1, 2.0),
        ],
        laws: vec![PiecewiseQuadraticLaw {
            key: "alias-law".into(),
            lower_bound: first,
            upper_bound: 2.0,
            knots: vec![first, 1.0, 2.0],
            interior_knot_owners: vec![InteriorKnotOwner::LowerPiece],
            pieces: vec![
                QuadraticResponsePiece {
                    mode_key: "tiny".into(),
                    a: 0.0,
                    b: 1.0,
                    c: 0.0,
                },
                QuadraticResponsePiece {
                    mode_key: "scale".into(),
                    a: 0.0,
                    b: 0.0,
                    c: scale,
                },
            ],
        }],
    };
    assert_eq!(
        propagate_structured_population(&model),
        Err(StructuredPropagationError::InvalidMoment {
            scope: StructuredMomentScope::Population,
            quantity: StructuredMomentQuantity::NormalizedObservation,
            issue: StructuredNumericIssue::Underflow,
        })
    );
}
