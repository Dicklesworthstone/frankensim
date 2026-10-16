//! I04.1 conservation-defect ontology conformance.
//!
//! The battery locks schema laws, spatial/time recomposition, orientation
//! reversal, unit-scale refusal, duplicate ownership, missing charts,
//! interval containment, migration goldens, and corrupt-receipt mutations.
//! It claims accounting-language semantics only: no physical defect is
//! detected, localized, or attributed here.

use fs_evidence::balance::{
    AccountRole, AccountingChart, AccountingRule, BalanceCodecError, BalanceDefectReceipt,
    BalanceDraft, BalanceError, BalanceRule, BalanceTerm, BoundedId, DefectState, Interval,
    QuantityKind, Semantics, SignConvention, SpatialSupport, SupportKind, TimeSupport,
    UncertaintySplit,
};
use fs_evidence::color::ColorRank;
use std::collections::BTreeSet;
use std::fmt::Write as _;

fn id(text: &str) -> BoundedId {
    BoundedId::new(text).expect("test identifier admits")
}

fn chart() -> AccountingChart {
    AccountingChart::new(
        id("energy-chart"),
        3,
        vec![
            (id("storage"), AccountRole::Storage),
            (id("boundary-flux"), AccountRole::BoundaryFlux),
            (id("production"), AccountRole::Production),
            (id("destruction"), AccountRole::Destruction),
            (id("solver-defect"), AccountRole::SolverDefect),
            (id("jump-work"), AccountRole::JumpWork),
        ],
    )
    .expect("test chart admits")
}

fn cells(ids: &[u64]) -> SpatialSupport {
    SpatialSupport::new(
        id("mesh-a"),
        SupportKind::Cells { dim: 3 },
        ids.iter().copied().collect::<BTreeSet<u64>>(),
    )
    .expect("test support admits")
}

fn window(start: i64, end: i64) -> TimeSupport {
    TimeSupport::Window {
        clock: id("sim-clock"),
        start,
        end,
    }
}

fn term(account: &str, owner: &str, lo: f64, hi: f64, color: ColorRank) -> BalanceTerm {
    BalanceTerm {
        account: id(account),
        owner: id(owner),
        value: Interval::new(lo, hi).expect("test interval admits"),
        color,
    }
}

fn draft(support: SpatialSupport, time: TimeSupport, terms: Vec<BalanceTerm>) -> BalanceDraft {
    BalanceDraft {
        quantity: QuantityKind::Energy,
        semantics: Semantics::ExtensiveAmount,
        sign: SignConvention::InflowPositive,
        scale_pow10: 0,
        chart: chart().reference(),
        support,
        time,
        producer: id("producer-a"),
        terms,
        expected_closure: Interval::point(0.0).expect("zero closure"),
        measured: DefectState::Bounded(Interval::new(-1e-12, 1e-12).expect("defect")),
        uncertainty: UncertaintySplit::new(1e-13, 0.0, 2e-13).expect("uncertainty"),
    }
}

fn admit(
    support: SpatialSupport,
    time: TimeSupport,
    terms: Vec<BalanceTerm>,
) -> BalanceDefectReceipt {
    BalanceDefectReceipt::admit(draft(support, time, terms), &chart()).expect("receipt admits")
}

fn rule_of<T: std::fmt::Debug>(result: Result<T, BalanceError>) -> BalanceRule {
    result.expect_err("must refuse").rule()
}

#[test]
fn identifier_bounds_and_alphabet_refuse_by_name() {
    assert_eq!(rule_of(BoundedId::new("")), BalanceRule::IdBounds);
    assert_eq!(
        rule_of(BoundedId::new(&"x".repeat(97))),
        BalanceRule::IdBounds
    );
    assert_eq!(rule_of(BoundedId::new("Upper")), BalanceRule::IdBounds);
    assert_eq!(rule_of(BoundedId::new("space here")), BalanceRule::IdBounds);
    assert!(BoundedId::new("mesh-a.42:x_1").is_ok());
    assert_eq!(BalanceRule::IdBounds.slug(), "balance-id-bounds");
}

#[test]
fn interval_refuses_nonfinite_inverted_and_signed_zero() {
    assert_eq!(
        rule_of(Interval::new(f64::NAN, 0.0)),
        BalanceRule::NonFiniteValue
    );
    assert_eq!(
        rule_of(Interval::new(0.0, f64::INFINITY)),
        BalanceRule::NonFiniteValue
    );
    assert_eq!(
        rule_of(Interval::new(2.0, 1.0)),
        BalanceRule::InvertedInterval
    );
    assert_eq!(
        rule_of(Interval::new(-0.0, 1.0)),
        BalanceRule::SignedZeroAlias
    );
    assert_eq!(
        rule_of(Interval::new(-1.0, -0.0)),
        BalanceRule::SignedZeroAlias
    );
}

#[test]
fn uncertainty_refuses_negative_and_nonfinite_components() {
    assert_eq!(
        rule_of(UncertaintySplit::new(-1e-9, 0.0, 0.0)),
        BalanceRule::NegativeUncertainty
    );
    assert_eq!(
        rule_of(UncertaintySplit::new(0.0, f64::NAN, 0.0)),
        BalanceRule::NonFiniteValue
    );
    assert_eq!(
        rule_of(UncertaintySplit::new(0.0, 0.0, -0.0)),
        BalanceRule::NegativeUncertainty
    );
}

#[test]
fn quantity_kinds_pin_units_and_accounting_rules() {
    assert_eq!(QuantityKind::Energy.canonical_unit(), "J");
    assert_eq!(QuantityKind::Entropy.canonical_unit(), "J/K");
    assert_eq!(
        QuantityKind::Energy.accounting_rule(),
        AccountingRule::Conserved
    );
    assert_eq!(
        QuantityKind::Entropy.accounting_rule(),
        AccountingRule::ProductionNonNegative
    );
    assert_eq!(
        QuantityKind::Exergy.accounting_rule(),
        AccountingRule::DestructionNonNegative
    );
    assert_eq!(
        QuantityKind::SpeciesMass {
            species: id("iso-octane")
        }
        .accounting_rule(),
        AccountingRule::OpenLedger
    );
    assert_eq!(
        rule_of(QuantityKind::Element { atomic_number: 0 }.validate()),
        BalanceRule::ElementRange
    );
    assert_eq!(
        rule_of(QuantityKind::Element { atomic_number: 119 }.validate()),
        BalanceRule::ElementRange
    );
    assert_eq!(
        rule_of(
            QuantityKind::Momentum {
                frame: id("world"),
                axis: 3
            }
            .validate()
        ),
        BalanceRule::SupportDegreeRange
    );
}

#[test]
fn entropy_refuses_negative_production_terms() {
    let mut d = draft(
        cells(&[1]),
        window(0, 10),
        vec![term(
            "production",
            "combustor",
            -1.0,
            1.0,
            ColorRank::Verified,
        )],
    );
    d.quantity = QuantityKind::Entropy;
    let refused = BalanceDefectReceipt::admit(d, &chart()).expect_err("second law");
    assert_eq!(refused.rule(), BalanceRule::AccountingRuleViolation);
    assert_eq!(refused.rule().slug(), "balance-accounting-rule-violation");

    let mut ok = draft(
        cells(&[1]),
        window(0, 10),
        vec![term(
            "production",
            "combustor",
            0.0,
            1.0,
            ColorRank::Verified,
        )],
    );
    ok.quantity = QuantityKind::Entropy;
    assert!(BalanceDefectReceipt::admit(ok, &chart()).is_ok());

    let energy_negative_production = draft(
        cells(&[1]),
        window(0, 10),
        vec![term(
            "production",
            "model-defect",
            -1.0,
            1.0,
            ColorRank::Verified,
        )],
    );
    assert!(
        BalanceDefectReceipt::admit(energy_negative_production, &chart()).is_ok(),
        "conserved quantities may carry signed model-owned production defects"
    );
}

#[test]
fn missing_chart_account_and_duplicate_owner_refuse() {
    let missing = draft(
        cells(&[1]),
        window(0, 10),
        vec![term("radiation", "wall", 0.0, 1.0, ColorRank::Verified)],
    );
    assert_eq!(
        rule_of(BalanceDefectReceipt::admit(missing, &chart())),
        BalanceRule::MissingChartAccount
    );

    let duplicate = draft(
        cells(&[1]),
        window(0, 10),
        vec![
            term("storage", "cell-owner", 0.0, 1.0, ColorRank::Verified),
            term("storage", "cell-owner", 1.0, 2.0, ColorRank::Verified),
        ],
    );
    assert_eq!(
        rule_of(BalanceDefectReceipt::admit(duplicate, &chart())),
        BalanceRule::DuplicateTermOwner
    );

    assert_eq!(
        rule_of(AccountingChart::new(
            id("c"),
            1,
            vec![
                (id("storage"), AccountRole::Storage),
                (id("storage"), AccountRole::Remap),
            ],
        )),
        BalanceRule::DuplicateChartAccount
    );
}

#[test]
fn chart_reference_binds_exact_account_set() {
    let same = chart();
    let mut receipt_draft = draft(cells(&[1]), window(0, 10), Vec::new());
    let other_chart = AccountingChart::new(
        id("energy-chart"),
        3,
        vec![(id("storage"), AccountRole::Storage)],
    )
    .expect("chart admits");
    assert_eq!(
        rule_of(BalanceDefectReceipt::admit(
            receipt_draft.clone(),
            &other_chart
        )),
        BalanceRule::ChartMismatch,
        "same id and version with a different account set must refuse via digest"
    );
    receipt_draft.chart = same.reference();
    assert!(BalanceDefectReceipt::admit(receipt_draft, &same).is_ok());
}

#[test]
fn extensive_amounts_refuse_instant_supports() {
    let mut d = draft(cells(&[1]), window(0, 10), Vec::new());
    d.time = TimeSupport::Instant {
        clock: id("sim-clock"),
        tick: 5,
    };
    assert_eq!(
        rule_of(BalanceDefectReceipt::admit(d, &chart())),
        BalanceRule::SemanticsMismatch
    );
}

#[test]
fn window_order_and_support_shape_refuse() {
    assert_eq!(rule_of(window(10, 10).validate()), BalanceRule::TimeOrder);
    assert_eq!(rule_of(window(11, 10).validate()), BalanceRule::TimeOrder);
    assert_eq!(
        rule_of(SpatialSupport::new(
            id("mesh-a"),
            SupportKind::Global,
            BTreeSet::from([1u64]),
        )),
        BalanceRule::SupportShape
    );
    assert_eq!(
        rule_of(SpatialSupport::new(
            id("mesh-a"),
            SupportKind::Cells { dim: 3 },
            BTreeSet::new(),
        )),
        BalanceRule::SupportShape
    );
    assert_eq!(
        rule_of(SpatialSupport::new(
            id("mesh-a"),
            SupportKind::Cells { dim: 16 },
            BTreeSet::from([1u64]),
        )),
        BalanceRule::SupportDegreeRange
    );
}

#[test]
fn partition_composition_sums_terms_and_unions_support() {
    let left = admit(
        cells(&[1, 2]),
        window(0, 10),
        vec![
            term("storage", "owner-a", 1.0, 2.0, ColorRank::Verified),
            term("boundary-flux", "owner-a", -1.0, 1.0, ColorRank::Verified),
        ],
    );
    let right = admit(
        cells(&[3, 4]),
        window(0, 10),
        vec![
            term("storage", "owner-a", 10.0, 20.0, ColorRank::Validated),
            term("jump-work", "owner-b", 0.5, 0.5, ColorRank::Verified),
        ],
    );
    let composed = left
        .compose_partition(&right, id("composer"))
        .expect("disjoint partition composes");

    assert_eq!(composed.support().ids(), &BTreeSet::from([1u64, 2, 3, 4]));
    assert_eq!(composed.terms().len(), 3);
    let storage = composed
        .terms()
        .iter()
        .find(|t| t.account.as_str() == "storage")
        .expect("merged storage term");
    assert!(
        storage.value.lo() <= 11.0 && 22.0 <= storage.value.hi(),
        "merged interval [{}, {}] must contain the exact sum [11, 22]",
        storage.value.lo(),
        storage.value.hi()
    );
    assert_eq!(storage.color, ColorRank::Validated, "weakest wins per term");
    assert_eq!(composed.evidence_rank(), ColorRank::Validated);
    assert_eq!(composed.lineage().len(), 2);
    assert_eq!(composed.lineage()[0], left.content_id());
    assert_eq!(composed.lineage()[1], right.content_id());
    assert_eq!(composed.producer().as_str(), "composer");
}

#[test]
fn overlapping_or_foreign_supports_refuse_double_counting() {
    let left = admit(cells(&[1, 2]), window(0, 10), Vec::new());
    let overlap = admit(cells(&[2, 3]), window(0, 10), Vec::new());
    assert_eq!(
        rule_of(left.compose_partition(&overlap, id("composer"))),
        BalanceRule::SupportOverlap
    );

    let other_complex = BalanceDefectReceipt::admit(
        draft(
            SpatialSupport::new(
                id("mesh-b"),
                SupportKind::Cells { dim: 3 },
                BTreeSet::from([9u64]),
            )
            .expect("support admits"),
            window(0, 10),
            Vec::new(),
        ),
        &chart(),
    )
    .expect("receipt admits");
    assert_eq!(
        rule_of(left.compose_partition(&other_complex, id("composer"))),
        BalanceRule::SupportMismatch
    );
}

#[test]
fn window_composition_needs_adjacency_and_extensive_semantics() {
    let a = admit(cells(&[1]), window(0, 10), Vec::new());
    let b = admit(cells(&[1]), window(10, 25), Vec::new());
    let joined = a.compose_windows(&b, id("composer")).expect("adjacent");
    assert_eq!(
        joined.time(),
        &window(0, 25),
        "adjacent windows merge into one"
    );

    let gap = admit(cells(&[1]), window(11, 20), Vec::new());
    assert_eq!(
        rule_of(a.compose_windows(&gap, id("composer"))),
        BalanceRule::WindowsNotAdjacent
    );
    let overlap = admit(cells(&[1]), window(9, 20), Vec::new());
    assert_eq!(
        rule_of(a.compose_windows(&overlap, id("composer"))),
        BalanceRule::WindowsNotAdjacent
    );

    let mut rate = draft(cells(&[1]), window(0, 10), Vec::new());
    rate.semantics = Semantics::Rate;
    let rate_a = BalanceDefectReceipt::admit(rate.clone(), &chart()).expect("rate admits");
    let mut rate_b_draft = draft(cells(&[1]), window(10, 20), Vec::new());
    rate_b_draft.semantics = Semantics::Rate;
    let rate_b = BalanceDefectReceipt::admit(rate_b_draft, &chart()).expect("rate admits");
    assert_eq!(
        rule_of(rate_a.compose_windows(&rate_b, id("composer"))),
        BalanceRule::SemanticsMismatch,
        "rates must not silently average across windows"
    );

    let elsewhere = admit(cells(&[2]), window(10, 20), Vec::new());
    assert_eq!(
        rule_of(a.compose_windows(&elsewhere, id("composer"))),
        BalanceRule::SupportMismatch
    );
}

#[test]
fn header_mismatches_refuse_by_name() {
    let base = admit(cells(&[1]), window(0, 10), Vec::new());

    let mut quantity = draft(cells(&[2]), window(0, 10), Vec::new());
    quantity.quantity = QuantityKind::Charge;
    let quantity = BalanceDefectReceipt::admit(quantity, &chart()).expect("admits");
    assert_eq!(
        rule_of(base.compose_partition(&quantity, id("c"))),
        BalanceRule::QuantityMismatch
    );

    let mut sign = draft(cells(&[2]), window(0, 10), Vec::new());
    sign.sign = SignConvention::OutflowPositive;
    let sign = BalanceDefectReceipt::admit(sign, &chart()).expect("admits");
    assert_eq!(
        rule_of(base.compose_partition(&sign, id("c"))),
        BalanceRule::SignConventionMismatch
    );

    let mut scale = draft(cells(&[2]), window(0, 10), Vec::new());
    scale.scale_pow10 = 3;
    let scale = BalanceDefectReceipt::admit(scale, &chart()).expect("admits");
    assert_eq!(
        rule_of(base.compose_partition(&scale, id("c"))),
        BalanceRule::ScaleMismatch,
        "unit scaling is identity, never silent conversion"
    );

    let mut clock = draft(cells(&[2]), window(0, 10), Vec::new());
    clock.time = TimeSupport::Window {
        clock: id("wall-clock"),
        start: 0,
        end: 10,
    };
    let clock = BalanceDefectReceipt::admit(clock, &chart()).expect("admits");
    assert_eq!(
        rule_of(base.compose_partition(&clock, id("c"))),
        BalanceRule::ClockMismatch
    );
}

#[test]
fn orientation_reversal_is_a_bitwise_involution() {
    let receipt = admit(
        cells(&[1, 2]),
        window(0, 10),
        vec![
            term("storage", "owner-a", 1.0, 2.5, ColorRank::Verified),
            term("boundary-flux", "owner-b", -3.0, 0.0, ColorRank::Estimated),
        ],
    );
    let reversed = receipt.reverse_orientation();
    assert_eq!(reversed.sign(), SignConvention::OutflowPositive);
    let flux = reversed
        .terms()
        .iter()
        .find(|t| t.account.as_str() == "boundary-flux")
        .expect("flux term");
    assert_eq!(
        flux.value.lo().to_bits(),
        0.0f64.to_bits(),
        "-0.0 never appears"
    );
    assert_eq!(flux.value.hi().to_bits(), 3.0f64.to_bits());

    let round_trip = reversed.reverse_orientation();
    assert_eq!(round_trip, receipt, "double reversal restores the receipt");
    assert_eq!(
        round_trip.canonical_bytes(),
        receipt.canonical_bytes(),
        "double reversal restores the exact canonical bytes"
    );
}

#[test]
fn defect_states_stay_distinct_through_composition() {
    let zero = DefectState::ExactZero;
    let bounded = DefectState::Bounded(Interval::new(-1.0, 2.0).expect("i"));
    let unknown = DefectState::Unknown;
    let unowned = DefectState::UnownedRemainder(Interval::new(0.5, 0.75).expect("i"));
    let inapplicable = DefectState::Inapplicable {
        reason: id("rigid-body-has-no-thermal-ledger"),
    };

    assert_eq!(zero.compose(&zero).expect("zz"), DefectState::ExactZero);
    assert_eq!(zero.compose(&bounded).expect("zb"), bounded);
    assert_eq!(bounded.compose(&unknown).expect("bu"), DefectState::Unknown);
    match unowned.compose(&bounded).expect("ub") {
        DefectState::UnownedRemainder(i) => {
            assert!(i.lo() <= -0.5 && 2.75 <= i.hi());
        }
        other => panic!("unowned remainder laundered into {other:?}"),
    }
    assert_eq!(
        unowned.compose(&zero).expect("uz"),
        unowned,
        "missing ownership stays visible"
    );
    assert_eq!(
        inapplicable
            .compose(&inapplicable)
            .expect("same reason composes"),
        inapplicable
    );
    let err = inapplicable.compose(&bounded).expect_err("ib refuses");
    assert_eq!(err.rule(), BalanceRule::DefectStateIncompatible);
    let other_reason = DefectState::Inapplicable {
        reason: id("different-reason"),
    };
    assert_eq!(
        inapplicable
            .compose(&other_reason)
            .expect_err("reasons differ")
            .rule(),
        BalanceRule::DefectStateIncompatible
    );
}

#[test]
fn composed_intervals_contain_the_exact_sums() {
    let a = Interval::new(0.1, 0.2).expect("a");
    let b = Interval::new(0.3, 0.4).expect("b");
    let sum = a.add(&b);
    assert!(sum.lo() <= 0.1 + 0.3 && 0.2 + 0.4 <= sum.hi());
    assert!(sum.lo() < sum.hi());
    let widened = Interval::new(sum.lo(), sum.hi()).expect("still canonical");
    assert!(widened.contains(&Interval::new(0.4, 0.6).expect("exact")));
}

#[test]
fn evidence_rank_is_weakest_wins_and_termless_receipts_stay_estimated() {
    let strong = admit(
        cells(&[1]),
        window(0, 10),
        vec![
            term("storage", "a", 0.0, 1.0, ColorRank::Verified),
            term("jump-work", "b", 0.0, 1.0, ColorRank::Verified),
        ],
    );
    assert_eq!(strong.evidence_rank(), ColorRank::Verified);
    let mixed = admit(
        cells(&[2]),
        window(0, 10),
        vec![
            term("storage", "a", 0.0, 1.0, ColorRank::Verified),
            term("solver-defect", "s", 0.0, 1.0, ColorRank::Estimated),
        ],
    );
    assert_eq!(mixed.evidence_rank(), ColorRank::Estimated);
    let composed = strong
        .compose_partition(&mixed, id("composer"))
        .expect("composes");
    assert_eq!(
        composed.evidence_rank(),
        ColorRank::Estimated,
        "composition never outranks the weakest operand"
    );
    let empty = admit(cells(&[3]), window(0, 10), Vec::new());
    assert_eq!(empty.evidence_rank(), ColorRank::Estimated);
}

#[test]
fn canonical_roundtrip_preserves_identity() {
    let receipt = admit(
        cells(&[1, 2, 7]),
        window(-5, 40),
        vec![
            term("storage", "owner-a", 1.0, 2.0, ColorRank::Verified),
            term("boundary-flux", "owner-b", -0.5, 0.5, ColorRank::Validated),
        ],
    );
    let bytes = receipt.canonical_bytes();
    let decoded = BalanceDefectReceipt::decode(&bytes).expect("roundtrip decodes");
    assert_eq!(decoded, receipt);
    assert_eq!(decoded.content_id(), receipt.content_id());

    let composed = receipt
        .compose_partition(
            &admit(cells(&[9]), window(-5, 40), Vec::new()),
            id("composer"),
        )
        .expect("composes");
    let composed_decoded =
        BalanceDefectReceipt::decode(&composed.canonical_bytes()).expect("decodes");
    assert_eq!(composed_decoded.lineage(), composed.lineage());
}

#[test]
fn transport_refuses_trailing_bytes_truncation_and_bad_magic() {
    let receipt = admit(cells(&[1]), window(0, 10), Vec::new());
    let bytes = receipt.canonical_bytes();

    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(BalanceDefectReceipt::decode(&trailing).is_err());

    for cut in [0usize, 1, 4, bytes.len() / 2, bytes.len() - 1] {
        assert!(
            BalanceDefectReceipt::decode(&bytes[..cut]).is_err(),
            "truncation at {cut} must refuse"
        );
    }

    let mut magic = bytes.clone();
    magic[0] ^= 0xFF;
    let err = BalanceDefectReceipt::decode(&magic).expect_err("bad magic refuses");
    assert_eq!(err.rule_name(), "balance-canonical-identity");
    let _ = BalanceCodecError::to_string(&err);
}

#[test]
fn corrupt_mutations_never_alias_the_original_identity() {
    let receipt = admit(
        cells(&[1, 2]),
        window(0, 10),
        vec![term("storage", "owner-a", 1.0, 2.0, ColorRank::Verified)],
    );
    let bytes = receipt.canonical_bytes();
    let original = receipt.content_id();
    for index in 0..bytes.len() {
        let mut mutated = bytes.clone();
        mutated[index] ^= 0x01;
        match BalanceDefectReceipt::decode(&mutated) {
            Err(_) => {}
            Ok(other) => {
                assert_ne!(
                    other.content_id(),
                    original,
                    "byte {index} flip decoded to the original identity"
                );
            }
        }
    }
}

#[test]
fn migration_golden_v1_bytes_and_identity_stay_decodable() {
    const PINNED_CONTENT_ID: &str =
        "5fe3c97c041e4e8434902bbca3212df85828dd86b8ea5cc709ae296e9e622bef";

    let receipt = admit(
        cells(&[3, 5]),
        window(0, 16),
        vec![
            term("boundary-flux", "port-1", -2.0, -1.0, ColorRank::Validated),
            term("storage", "cell-owner", 0.25, 0.5, ColorRank::Verified),
        ],
    );
    let bytes = receipt.canonical_bytes();
    let mut hex = String::with_capacity(bytes.len() * 2);
    for byte in &bytes {
        write!(&mut hex, "{byte:02x}").expect("writing to a String cannot fail");
    }

    // Schema-v1 golden: the canonical header (magic + version) and the exact
    // content identity of the receipt above. A change here is a WIRE BREAK
    // and needs a version bump plus a migration path, not a golden edit.
    assert!(
        hex.starts_with("4653424400000001"),
        "v1 canonical header drifted: {}",
        &hex[..16.min(hex.len())]
    );
    assert_eq!(
        receipt.content_id().to_hex(),
        PINNED_CONTENT_ID,
        "v1 content identity drifted (full canonical hex: {hex})"
    );

    let decoded = BalanceDefectReceipt::decode(&bytes).expect("v1 golden decodes");
    assert_eq!(decoded.schema_version(), 1);
    assert_eq!(decoded, receipt);
}

#[test]
fn unknown_schema_version_refuses() {
    let receipt = admit(cells(&[1]), window(0, 10), Vec::new());
    let mut bytes = receipt.canonical_bytes();
    bytes[7] = 2;
    let err = BalanceDefectReceipt::decode(&bytes).expect_err("v2 unknown");
    assert!(err.detail().contains("unknown schema version"));
}

#[test]
fn spatial_then_time_recomposition_matches_time_then_spatial_states() {
    let q = |ids: &[u64], t0: i64, t1: i64, lo: f64, hi: f64| {
        let mut d = draft(
            cells(ids),
            window(t0, t1),
            vec![term("storage", "owner-a", lo, hi, ColorRank::Verified)],
        );
        d.measured = DefectState::Bounded(Interval::new(lo, hi).expect("i"));
        BalanceDefectReceipt::admit(d, &chart()).expect("admits")
    };
    let a0 = q(&[1], 0, 10, 0.1, 0.2);
    let b0 = q(&[2], 0, 10, 0.3, 0.4);
    let a1 = q(&[1], 10, 20, 0.5, 0.6);
    let b1 = q(&[2], 10, 20, 0.7, 0.8);

    let space_first = a0
        .compose_partition(&b0, id("c"))
        .expect("s0")
        .compose_windows(&a1.compose_partition(&b1, id("c")).expect("s1"), id("c"))
        .expect("windows");
    let time_first = a0
        .compose_windows(&a1, id("c"))
        .expect("t-a")
        .compose_partition(&b0.compose_windows(&b1, id("c")).expect("t-b"), id("c"))
        .expect("partition");

    assert_eq!(space_first.support(), time_first.support());
    assert_eq!(space_first.time(), time_first.time());
    let exact_sum = 0.1 + 0.3 + 0.5 + 0.7;
    let exact_hi = 0.2 + 0.4 + 0.6 + 0.8;
    for receipt in [&space_first, &time_first] {
        match receipt.measured() {
            DefectState::Bounded(i) => {
                assert!(
                    i.lo() <= exact_sum && exact_hi <= i.hi(),
                    "recomposed bound [{}, {}] must contain [{exact_sum}, {exact_hi}]",
                    i.lo(),
                    i.hi()
                );
            }
            other => panic!("recomposition changed the state family: {other:?}"),
        }
    }
}
