//! G0/G3/G5 conformance for the structured gear-backlash consumer.

use std::num::NonZeroU64;

use fs_toleralloc::{
    GEAR_BACKLASH_CONSUMER_SCHEMA_V1, GearBacklashConsumerDraftV1, GearBacklashConsumptionErrorV1,
    GearBacklashLengthUnitV1, GearBacklashProbabilityErrorV1, GearBacklashProbabilityV1,
    InteriorKnotOwner, MAX_EXACT_STRUCTURED_WEIGHT_V1, MAX_GEAR_BACKLASH_QUANTILES_V1,
    PiecewiseQuadraticLaw, QuadraticResponsePiece, STRUCTURED_PROPAGATION_SCHEMA_V1,
    StructuredLawId, StructuredModelIdentity, StructuredNodeId, StructuredNodeSpec,
    StructuredPopulationModel, propagate_structured_population,
};

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("fixture multiplicities are positive")
}

fn probability(numerator: u64, denominator: u64) -> GearBacklashProbabilityV1 {
    GearBacklashProbabilityV1::try_new(numerator, denominator)
        .expect("fixture probability is inside the closed unit interval")
}

fn branch(key: &str, parent: Option<u32>) -> StructuredNodeSpec {
    StructuredNodeSpec::Branch {
        key: key.into(),
        parent: parent.map(StructuredNodeId),
    }
}

fn leaf(key: &str, parent: u32, weight: u64, raw_clearance: f64) -> StructuredNodeSpec {
    leaf_with_law(key, parent, weight, raw_clearance, StructuredLawId(0))
}

fn leaf_with_law(
    key: &str,
    parent: u32,
    weight: u64,
    raw_clearance: f64,
    law: StructuredLawId,
) -> StructuredNodeSpec {
    StructuredNodeSpec::Leaf {
        key: key.into(),
        parent: StructuredNodeId(parent),
        relative_weight: nz(weight),
        raw_clearance,
        law,
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

fn quantile_requests() -> Vec<GearBacklashProbabilityV1> {
    vec![
        probability(1, 1),
        probability(53, 64),
        probability(0, 7),
        probability(1, 2),
        probability(31, 32),
        probability(1, 64),
        probability(13, 16),
        probability(1, 16),
        probability(61, 64),
        probability(5, 64),
    ]
}

#[test]
#[allow(clippy::too_many_lines)] // One oracle checks every retained support/quantile/mode field.
fn g0_weighted_support_quantiles_units_and_modes_match_the_integer_oracle() {
    let model = gear_model();
    let structured = propagate_structured_population(&model).expect("structured fixture evaluates");
    let report = GearBacklashConsumerDraftV1 {
        output_unit: GearBacklashLengthUnitV1::Micrometre,
        quantiles: quantile_requests(),
    }
    .consume(&structured)
    .expect("admitted receipt is consumable");

    assert_eq!(report.schema_version(), GEAR_BACKLASH_CONSUMER_SCHEMA_V1);
    assert_eq!(
        report.structured_schema_version(),
        STRUCTURED_PROPAGATION_SCHEMA_V1
    );
    assert_eq!(report.model(), &model);
    assert_eq!(report.model_identity(), &model.identity);
    assert_eq!(report.output_unit(), GearBacklashLengthUnitV1::Micrometre);
    assert_eq!(report.output_unit().tag(), 3);
    assert_eq!(report.output_unit().symbol(), "um");
    assert_eq!(report.total_weight(), 64);
    assert_eq!(report.mean().unit(), GearBacklashLengthUnitV1::Micrometre);
    assert_eq!(
        report.variance().unit(),
        GearBacklashLengthUnitV1::Micrometre
    );
    assert_eq!(
        report.standard_deviation().unit(),
        GearBacklashLengthUnitV1::Micrometre
    );

    let scale = GearBacklashLengthUnitV1::Micrometre.metres_per_unit();
    assert_eq!(
        report.mean().source_value().to_bits(),
        structured.mean().to_bits()
    );
    assert_eq!(
        report.mean().metres().to_bits(),
        (structured.mean() * scale).to_bits()
    );
    assert_eq!(
        report.variance().source_variance().to_bits(),
        structured.variance().to_bits()
    );
    assert_eq!(
        report.variance().square_metres().to_bits(),
        (structured.variance() * scale * scale).to_bits()
    );
    assert_eq!(
        report.standard_deviation().source_value().to_bits(),
        structured.standard_deviation().to_bits()
    );
    assert_eq!(
        report.standard_deviation().metres().to_bits(),
        (structured.standard_deviation() * scale).to_bits()
    );

    let expected_support = [
        (-4.0_f64, 1_u64, 0_u64, 1_u64),
        (-1.0, 3, 1, 4),
        (0.0, 48, 4, 52),
        (1.0, 9, 52, 61),
        (4.0, 3, 61, 64),
    ];
    assert_eq!(report.support().len(), expected_support.len());
    for (actual, (source, weight, before, at)) in report.support().iter().zip(expected_support) {
        assert_eq!(actual.value().unit(), GearBacklashLengthUnitV1::Micrometre);
        assert_eq!(actual.value().source_value().to_bits(), source.to_bits());
        assert_eq!(
            actual.value().metres().to_bits(),
            (source * scale).to_bits()
        );
        assert_eq!(actual.relative_weight(), weight);
        assert_eq!(actual.cumulative_before(), before);
        assert_eq!(actual.cumulative_at(), at);
    }

    let expected_quantiles = [
        (probability(0, 1), 0_usize, -4.0_f64, 0_u64, 1_u64),
        (probability(1, 64), 0, -4.0, 0, 1),
        (probability(1, 16), 1, -1.0, 1, 4),
        (probability(5, 64), 2, 0.0, 4, 52),
        (probability(1, 2), 2, 0.0, 4, 52),
        (probability(13, 16), 2, 0.0, 4, 52),
        (probability(53, 64), 3, 1.0, 52, 61),
        (probability(61, 64), 3, 1.0, 52, 61),
        (probability(31, 32), 4, 4.0, 61, 64),
        (probability(1, 1), 4, 4.0, 61, 64),
    ];
    assert_eq!(report.quantiles().len(), expected_quantiles.len());
    for (actual, (level, support_index, source, before, at)) in
        report.quantiles().iter().zip(expected_quantiles)
    {
        assert_eq!(actual.probability(), level);
        assert_eq!(actual.support_index(), support_index);
        assert_eq!(actual.value().source_value().to_bits(), source.to_bits());
        assert_eq!(actual.cumulative_before(), before);
        assert_eq!(actual.cumulative_at(), at);
    }

    let expected_modes = [
        ("drive-negative", 4_u64, 2_usize),
        ("lash", 48, 4),
        ("drive-positive", 12, 2),
    ];
    assert_eq!(report.modes().len(), expected_modes.len());
    for (piece, (actual, (key, weight, leaf_count))) in
        report.modes().iter().zip(expected_modes).enumerate()
    {
        assert_eq!(actual.law(), StructuredLawId(0));
        assert_eq!(actual.law_key(), "backlash-clearance");
        assert_eq!(actual.piece(), piece);
        assert_eq!(actual.mode_key(), key);
        assert_eq!(actual.relative_weight(), weight);
        assert_eq!(actual.leaf_count(), leaf_count);
    }
}

#[test]
fn g5_request_order_and_equivalent_fraction_spellings_are_nonsemantic() {
    let structured =
        propagate_structured_population(&gear_model()).expect("structured fixture evaluates");
    let canonical = GearBacklashConsumerDraftV1 {
        output_unit: GearBacklashLengthUnitV1::Micrometre,
        quantiles: quantile_requests(),
    }
    .consume(&structured)
    .expect("canonical request evaluates");
    let reordered = GearBacklashConsumerDraftV1 {
        output_unit: GearBacklashLengthUnitV1::Micrometre,
        quantiles: vec![
            probability(10, 128),
            probability(62, 64),
            probability(2, 4),
            probability(9, 9),
            probability(2, 32),
            probability(53, 64),
            probability(0, 9),
            probability(2, 128),
            probability(26, 32),
            probability(61, 64),
        ],
    }
    .consume(&structured)
    .expect("equivalent reordered request evaluates");

    assert_eq!(canonical, reordered);
}

#[test]
fn g0_multiple_laws_retain_unobserved_zero_weight_modes() {
    let secondary_law = PiecewiseQuadraticLaw {
        key: "secondary-response".into(),
        lower_bound: -2.0,
        upper_bound: 2.0,
        knots: vec![-2.0, 0.0, 2.0],
        interior_knot_owners: vec![InteriorKnotOwner::UpperPiece],
        pieces: vec![
            QuadraticResponsePiece {
                mode_key: "secondary-negative".into(),
                a: 0.0,
                b: 1.0,
                c: 0.0,
            },
            QuadraticResponsePiece {
                mode_key: "secondary-nonnegative".into(),
                a: 0.0,
                b: 1.0,
                c: 0.0,
            },
        ],
    };
    let model = StructuredPopulationModel {
        identity: StructuredModelIdentity::try_new(
            "gear/backlash-multiple-laws",
            nz(1),
            [0xa7; 32],
        )
        .expect("fixture identity is canonical"),
        nodes: vec![
            branch("population", None),
            leaf_with_law("primary-lash", 0, 3, 0.0, StructuredLawId(0)),
            leaf_with_law("secondary-negative", 0, 5, -1.0, StructuredLawId(1)),
        ],
        laws: vec![backlash_law(), secondary_law],
    };
    let structured = propagate_structured_population(&model).expect("multi-law fixture evaluates");
    let report = GearBacklashConsumerDraftV1 {
        output_unit: GearBacklashLengthUnitV1::Micrometre,
        quantiles: vec![probability(1, 2)],
    }
    .consume(&structured)
    .expect("multi-law receipt is consumable");

    let actual = report
        .modes()
        .iter()
        .map(|mode| {
            (
                mode.law(),
                mode.law_key(),
                mode.piece(),
                mode.mode_key(),
                mode.relative_weight(),
                mode.leaf_count(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        actual,
        vec![
            (
                StructuredLawId(0),
                "backlash-clearance",
                0,
                "drive-negative",
                0,
                0,
            ),
            (StructuredLawId(0), "backlash-clearance", 1, "lash", 3, 1),
            (
                StructuredLawId(0),
                "backlash-clearance",
                2,
                "drive-positive",
                0,
                0,
            ),
            (
                StructuredLawId(1),
                "secondary-response",
                0,
                "secondary-negative",
                5,
                1,
            ),
            (
                StructuredLawId(1),
                "secondary-response",
                1,
                "secondary-nonnegative",
                0,
                0,
            ),
        ]
    );
    assert_eq!(report.total_weight(), 8);
}

#[test]
fn g0_probability_and_request_admission_refuse_ambiguity_before_consumption() {
    assert_eq!(probability(2, 4), probability(1, 2));
    assert_eq!(probability(0, 9), probability(0, 1));
    assert_eq!(
        GearBacklashProbabilityV1::try_new(0, 0),
        Err(GearBacklashProbabilityErrorV1::ZeroDenominator)
    );
    assert_eq!(
        GearBacklashProbabilityV1::try_new(2, 1),
        Err(GearBacklashProbabilityErrorV1::AboveOne)
    );

    let structured =
        propagate_structured_population(&gear_model()).expect("structured fixture evaluates");
    assert!(matches!(
        GearBacklashConsumerDraftV1 {
            output_unit: GearBacklashLengthUnitV1::Millimetre,
            quantiles: Vec::new(),
        }
        .consume(&structured),
        Err(GearBacklashConsumptionErrorV1::NoQuantiles)
    ));
    assert!(matches!(
        GearBacklashConsumerDraftV1 {
            output_unit: GearBacklashLengthUnitV1::Millimetre,
            quantiles: vec![probability(0, 1); MAX_GEAR_BACKLASH_QUANTILES_V1 + 1],
        }
        .consume(&structured),
        Err(GearBacklashConsumptionErrorV1::QuantileLimit {
            actual,
            max: MAX_GEAR_BACKLASH_QUANTILES_V1,
        }) if actual == MAX_GEAR_BACKLASH_QUANTILES_V1 + 1
    ));
    assert_eq!(
        GearBacklashConsumerDraftV1 {
            output_unit: GearBacklashLengthUnitV1::Millimetre,
            quantiles: vec![probability(1, 2), probability(2, 4)],
        }
        .consume(&structured),
        Err(GearBacklashConsumptionErrorV1::DuplicateQuantile {
            probability: probability(1, 2),
        })
    );

    let maximum_request = GearBacklashConsumerDraftV1 {
        output_unit: GearBacklashLengthUnitV1::Millimetre,
        quantiles: (0..MAX_GEAR_BACKLASH_QUANTILES_V1)
            .map(|numerator| {
                probability(
                    u64::try_from(numerator).expect("bounded numerator"),
                    u64::try_from(MAX_GEAR_BACKLASH_QUANTILES_V1 - 1).expect("bounded denominator"),
                )
            })
            .collect(),
    }
    .consume(&structured)
    .expect("the exact request cap is admitted");
    assert_eq!(
        maximum_request.quantiles().len(),
        MAX_GEAR_BACKLASH_QUANTILES_V1
    );
}

#[test]
fn g3_u128_cdf_thresholds_remain_exact_at_the_upstream_weight_cap() {
    let total = MAX_EXACT_STRUCTURED_WEIGHT_V1;
    let model = StructuredPopulationModel {
        identity: StructuredModelIdentity::try_new(
            "gear/backlash-rank-boundary",
            nz(1),
            [0x4b; 32],
        )
        .expect("fixture identity is canonical"),
        nodes: vec![
            branch("population", None),
            leaf("lower", 0, total - 1, -1.0),
            leaf("upper", 0, 1, 1.0),
        ],
        laws: vec![PiecewiseQuadraticLaw {
            key: "identity-response".into(),
            lower_bound: -1.0,
            upper_bound: 1.0,
            knots: vec![-1.0, 1.0],
            interior_knot_owners: Vec::new(),
            pieces: vec![QuadraticResponsePiece {
                mode_key: "identity".into(),
                a: 0.0,
                b: 1.0,
                c: 0.0,
            }],
        }],
    };
    let structured = propagate_structured_population(&model).expect("cap fixture evaluates");
    let twice_total = total.checked_mul(2).expect("2^54 fits u64");
    let report = GearBacklashConsumerDraftV1 {
        output_unit: GearBacklashLengthUnitV1::Metre,
        quantiles: vec![
            probability(0, 1),
            probability(total - 1, total),
            probability(twice_total - 1, twice_total),
            probability(1, 1),
        ],
    }
    .consume(&structured)
    .expect("exact integer thresholding consumes the cap fixture");

    assert_eq!(report.total_weight(), total);
    let selected = report
        .quantiles()
        .iter()
        .map(|quantile| quantile.support_index())
        .collect::<Vec<_>>();
    assert_eq!(selected, vec![0, 0, 1, 1]);
    assert_eq!(report.quantiles()[1].cumulative_at(), total - 1);
    assert_eq!(report.quantiles()[2].cumulative_before(), total - 1);
    assert_eq!(report.quantiles()[2].cumulative_at(), total);
}

#[test]
fn g5_complete_model_and_unit_semantics_move_the_structural_receipt() {
    let model = gear_model();
    let structured = propagate_structured_population(&model).expect("base model evaluates");
    let base = GearBacklashConsumerDraftV1 {
        output_unit: GearBacklashLengthUnitV1::Micrometre,
        quantiles: vec![probability(1, 2)],
    }
    .consume(&structured)
    .expect("base report evaluates");

    let millimetre = GearBacklashConsumerDraftV1 {
        output_unit: GearBacklashLengthUnitV1::Millimetre,
        quantiles: vec![probability(1, 2)],
    }
    .consume(&structured)
    .expect("unit reinterpretation evaluates");
    assert_eq!(
        base.mean().source_value().to_bits(),
        millimetre.mean().source_value().to_bits()
    );
    assert_ne!(
        base.mean().metres().to_bits(),
        millimetre.mean().metres().to_bits()
    );
    assert_ne!(base, millimetre);

    let mut changed_model = gear_model();
    let StructuredNodeSpec::Branch { key, .. } = &mut changed_model.nodes[1] else {
        panic!("fixture node one is a branch");
    };
    *key = "process-low-reidentified".into();
    let changed_structured =
        propagate_structured_population(&changed_model).expect("changed hierarchy evaluates");
    let changed = GearBacklashConsumerDraftV1 {
        output_unit: GearBacklashLengthUnitV1::Micrometre,
        quantiles: vec![probability(1, 2)],
    }
    .consume(&changed_structured)
    .expect("changed hierarchy report evaluates");

    assert_eq!(base.model_identity(), changed.model_identity());
    assert_ne!(base.model(), changed.model());
    assert_ne!(base, changed);
}
