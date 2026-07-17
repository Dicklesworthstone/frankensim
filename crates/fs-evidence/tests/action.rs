//! I08.1 evidence-action vocabulary conformance.
//!
//! The battery locks round trips, correlated/exclusive actions, dependency
//! cycles, unknown costs, expiry, duplicate proposals, unit/currency
//! refusal, malformed budgets, the planning/execution evidence split, and
//! migration goldens. It claims vocabulary semantics only: no portfolio is
//! optimized and no action is executed here.

use fs_evidence::action::{
    ActionDraft, ActionKind, ActionProposal, ActionRule, BoundedId, Budgets, CostState, CostUnit,
    ExecutionReceipt, ExpectedResponse, Expiry, NonNegRange, Portfolio, PortfolioCost, TargetSlice,
    UncertaintyComponent,
};
use fs_evidence::color::ColorRank;
use std::collections::BTreeSet;

fn id(text: &str) -> BoundedId {
    BoundedId::new(text).expect("test identifier admits")
}

fn known(lo: f64, hi: f64, unit: CostUnit) -> CostState {
    CostState::Known {
        range: NonNegRange::new(lo, hi).expect("test range admits"),
        unit,
    }
}

fn usd(lo: f64, hi: f64) -> CostState {
    known(lo, hi, CostUnit::currency("USD").expect("USD admits"))
}

fn budgets() -> Budgets {
    Budgets {
        money: usd(100.0, 200.0),
        compute: known(10.0, 20.0, CostUnit::CpuCoreSeconds),
        memory: known(1e9, 2e9, CostUnit::MemoryBytes),
        lead_time: known(0.0, 1.0, CostUnit::CalendarDays),
        capabilities: BTreeSet::new(),
    }
}

fn draft(kind: ActionKind, claim: &str) -> ActionDraft {
    ActionDraft {
        kind,
        target: TargetSlice {
            claim: id(claim),
            component: UncertaintyComponent::Numerical,
        },
        response: ExpectedResponse::Factor { lo: 0.5, hi: 0.9 },
        budgets: budgets(),
        dependencies: BTreeSet::new(),
        exclusivity_group: None,
        correlation_group: None,
        expiry: None,
        ceiling: ColorRank::Verified,
        proposer: id("planner-a"),
    }
}

fn admit(kind: ActionKind, claim: &str) -> ActionProposal {
    ActionProposal::admit(draft(kind, claim)).expect("proposal admits")
}

fn now() -> Expiry {
    Expiry {
        clock: id("campaign-clock"),
        tick: 1000,
    }
}

#[test]
fn kind_taxonomy_is_complete_and_physicality_is_pinned() {
    let kinds = [
        (ActionKind::SolverTolerance, false),
        (ActionKind::MeshRefinement, false),
        (ActionKind::TimeRefinement, false),
        (ActionKind::RepresentationEscalation, false),
        (ActionKind::UqSamples, false),
        (ActionKind::MaterialCouponTest, true),
        (ActionKind::SensorCampaign, true),
        (ActionKind::Falsification, false),
        (ActionKind::StandardsObligation, false),
        (ActionKind::Refusal, false),
    ];
    for (kind, physical) in kinds {
        assert_eq!(kind.is_physical(), physical, "{kind:?} physicality");
    }
}

#[test]
fn malformed_budgets_and_responses_refuse_by_name() {
    assert_eq!(
        NonNegRange::new(-1.0, 1.0).expect_err("negative").rule(),
        ActionRule::NegativeCost
    );
    assert_eq!(
        NonNegRange::new(2.0, 1.0).expect_err("inverted").rule(),
        ActionRule::InvertedRange
    );
    assert_eq!(
        NonNegRange::new(0.0, f64::NAN).expect_err("nan").rule(),
        ActionRule::NonFiniteValue
    );
    assert_eq!(
        CostUnit::currency("usd").expect_err("lowercase").rule(),
        ActionRule::CurrencyCode
    );
    assert_eq!(
        CostUnit::currency("US").expect_err("short").rule(),
        ActionRule::CurrencyCode
    );

    let mut bad_response = draft(ActionKind::UqSamples, "claim-q");
    bad_response.response = ExpectedResponse::Factor { lo: 0.0, hi: 0.5 };
    assert_eq!(
        ActionProposal::admit(bad_response)
            .expect_err("zero lo")
            .rule(),
        ActionRule::ResponseRange
    );
    let mut above_one = draft(ActionKind::UqSamples, "claim-q");
    above_one.response = ExpectedResponse::Factor { lo: 0.5, hi: 1.5 };
    assert_eq!(
        ActionProposal::admit(above_one).expect_err("hi > 1").rule(),
        ActionRule::ResponseRange
    );
    assert_eq!(ActionRule::ResponseRange.slug(), "action-response-range");
}

#[test]
fn cross_currency_and_cross_unit_sums_refuse() {
    let eur = known(50.0, 60.0, CostUnit::currency("EUR").expect("EUR"));
    let err = usd(1.0, 2.0).add(&eur).expect_err("currency mix");
    assert_eq!(err.rule(), ActionRule::UnitMismatch);
    assert!(err.detail().contains("policy"));

    let wall = known(1.0, 2.0, CostUnit::WallSeconds);
    let cpu = known(1.0, 2.0, CostUnit::CpuCoreSeconds);
    assert_eq!(
        wall.add(&cpu).expect_err("time kinds differ").rule(),
        ActionRule::UnitMismatch
    );

    let same = usd(1.0, 2.0).add(&usd(10.0, 20.0)).expect("same unit adds");
    match same {
        CostState::Known { range, .. } => {
            assert!(range.lo() <= 11.0 && 22.0 <= range.hi());
        }
        other => panic!("known sum became {other:?}"),
    }
}

#[test]
fn unknown_costs_stay_explicit_and_never_read_as_zero() {
    let unknown = CostState::Unknown {
        authority_gap: id("no-coupon-price-authority"),
    };
    match usd(5.0, 6.0).add(&unknown).expect("absorbs explicitly") {
        CostState::Unknown { authority_gap } => {
            assert_eq!(authority_gap.as_str(), "no-coupon-price-authority");
        }
        other => panic!("unknown laundered into {other:?}"),
    }

    let mut priced = draft(ActionKind::MaterialCouponTest, "claim-yield");
    priced.budgets.money = CostState::Unknown {
        authority_gap: id("no-coupon-price-authority"),
    };
    let unpriced = ActionProposal::admit(priced).expect("admits");
    let plain = admit(ActionKind::UqSamples, "claim-q");
    let portfolio = Portfolio::admit(vec![unpriced, plain], &now()).expect("admits");
    match portfolio.total_money() {
        PortfolioCost::Incomparable { gaps, .. } => {
            assert_eq!(gaps.len(), 1);
            assert_eq!(gaps[0].as_str(), "no-coupon-price-authority");
        }
        other => panic!("unknown money total hidden as {other:?}"),
    }
    match portfolio.total_compute() {
        PortfolioCost::Known { range, unit } => {
            assert_eq!(unit, CostUnit::CpuCoreSeconds);
            assert!(range.lo() <= 20.0 && 40.0 <= range.hi());
        }
        other => panic!("compute total should be known, got {other:?}"),
    }
}

#[test]
fn portfolios_with_mixed_currencies_are_incomparable_not_summed() {
    let mut eur_action = draft(ActionKind::SensorCampaign, "claim-flux");
    eur_action.budgets.money = known(9.0, 11.0, CostUnit::currency("EUR").expect("EUR"));
    let eur_action = ActionProposal::admit(eur_action).expect("admits");
    let usd_action = admit(ActionKind::UqSamples, "claim-q");
    let portfolio = Portfolio::admit(vec![eur_action, usd_action], &now()).expect("admits");
    match portfolio.total_money() {
        PortfolioCost::Incomparable { gaps, units_seen } => {
            assert!(gaps.is_empty());
            assert_eq!(units_seen, 2, "both currencies stay visible");
        }
        other => panic!("cross-currency total silently combined: {other:?}"),
    }
}

#[test]
fn duplicate_proposal_identity_is_idempotent_and_refuses() {
    let a = admit(ActionKind::UqSamples, "claim-q");
    let b = admit(ActionKind::UqSamples, "claim-q");
    assert_eq!(a.content_id(), b.content_id(), "identical content, one id");
    let err = Portfolio::admit(vec![a.clone(), b], &now()).expect_err("duplicate");
    assert_eq!(err.rule(), ActionRule::DuplicateProposal);

    let c = admit(ActionKind::UqSamples, "claim-r");
    assert_ne!(a.content_id(), c.content_id());
    assert!(Portfolio::admit(vec![a, c], &now()).is_ok());
}

#[test]
fn dependency_cycles_and_missing_dependencies_refuse() {
    let a = admit(ActionKind::MeshRefinement, "claim-a");
    let mut b_draft = draft(ActionKind::SolverTolerance, "claim-b");
    b_draft.dependencies = BTreeSet::from([a.content_id()]);
    let b = ActionProposal::admit(b_draft).expect("admits");

    let ordered = Portfolio::admit(vec![b.clone(), a.clone()], &now()).expect("dag admits");
    let order = ordered.order();
    let pos_a = order.iter().position(|x| *x == a.content_id()).expect("a");
    let pos_b = order.iter().position(|x| *x == b.content_id()).expect("b");
    assert!(pos_a < pos_b, "dependencies execute first");

    let mut dangling = draft(ActionKind::SolverTolerance, "claim-c");
    dangling.dependencies = BTreeSet::from([a.content_id()]);
    let dangling = ActionProposal::admit(dangling).expect("admits");
    assert_eq!(
        Portfolio::admit(vec![dangling], &now())
            .expect_err("dep absent")
            .rule(),
        ActionRule::MissingDependency
    );

    // Two-node cycle built by fixed-point: impossible to express with
    // content-addressed ids (an id depends on the dependency set), so the
    // canonical cycle test is the self-dependency…
    let mut self_dep = draft(ActionKind::Falsification, "claim-d");
    let probe = ActionProposal::admit(self_dep.clone()).expect("probe");
    self_dep.dependencies = BTreeSet::from([probe.content_id()]);
    let self_like = ActionProposal::admit(self_dep).expect("admits");
    // …whose id CHANGED by declaring the dependency, so it now depends on a
    // missing proposal — content addressing itself forbids true cycles.
    let err = Portfolio::admit(vec![self_like], &now()).expect_err("refuses");
    assert!(
        matches!(
            err.rule(),
            ActionRule::MissingDependency | ActionRule::DependencyCycle
        ),
        "cycle or dangling refusal, never silent admission: {err}"
    );
}

#[test]
fn exclusive_actions_refuse_and_correlated_actions_carry_their_group() {
    let mut coarse = draft(ActionKind::MeshRefinement, "claim-a");
    coarse.exclusivity_group = Some(id("mesh-ladder"));
    let mut fine = draft(ActionKind::MeshRefinement, "claim-b");
    fine.exclusivity_group = Some(id("mesh-ladder"));
    let coarse = ActionProposal::admit(coarse).expect("admits");
    let fine = ActionProposal::admit(fine).expect("admits");
    assert_eq!(
        Portfolio::admit(vec![coarse.clone(), fine], &now())
            .expect_err("exclusive")
            .rule(),
        ActionRule::ExclusivityViolation
    );

    let mut correlated = draft(ActionKind::UqSamples, "claim-c");
    correlated.correlation_group = Some(id("shared-surrogate"));
    let correlated = ActionProposal::admit(correlated).expect("admits");
    let bytes = correlated.canonical_bytes();
    let decoded = ActionProposal::decode(&bytes).expect("roundtrip");
    assert_eq!(
        decoded.correlation_group().map(BoundedId::as_str),
        Some("shared-surrogate")
    );
    assert!(Portfolio::admit(vec![coarse, correlated], &now()).is_ok());
}

#[test]
fn expiry_is_enforced_on_the_admission_clock() {
    let mut expiring = draft(ActionKind::SensorCampaign, "claim-t");
    expiring.expiry = Some(Expiry {
        clock: id("campaign-clock"),
        tick: 1000,
    });
    let expired = ActionProposal::admit(expiring.clone()).expect("admits");
    assert_eq!(
        Portfolio::admit(vec![expired], &now())
            .expect_err("tick 1000 <= now 1000")
            .rule(),
        ActionRule::Expired
    );

    expiring.expiry = Some(Expiry {
        clock: id("campaign-clock"),
        tick: 1001,
    });
    let alive = ActionProposal::admit(expiring.clone()).expect("admits");
    assert!(Portfolio::admit(vec![alive], &now()).is_ok());

    expiring.expiry = Some(Expiry {
        clock: id("wall-clock"),
        tick: 5000,
    });
    let foreign_clock = ActionProposal::admit(expiring).expect("admits");
    assert_eq!(
        Portfolio::admit(vec![foreign_clock], &now())
            .expect_err("clock differs")
            .rule(),
        ActionRule::ClockMismatch
    );
}

#[test]
fn planned_physical_tests_cannot_raise_evidence_until_executed() {
    let coupon = admit(ActionKind::MaterialCouponTest, "claim-yield");
    let portfolio = Portfolio::admit(vec![coupon.clone()], &now()).expect("admits");
    assert!(portfolio.contains_unexecuted_physical());
    assert_eq!(
        portfolio.planned_evidence_effect(),
        None,
        "a plan is never evidence"
    );

    let receipt = ExecutionReceipt::admit(
        &coupon,
        coupon.content_id(),
        ColorRank::Validated,
        now(),
        id("lab-b"),
    )
    .expect("executed receipt admits");
    assert_eq!(receipt.outcome(), ColorRank::Validated);
    assert_eq!(receipt.proposal(), coupon.content_id());

    let mut capped = draft(ActionKind::UqSamples, "claim-q");
    capped.ceiling = ColorRank::Estimated;
    let capped = ActionProposal::admit(capped).expect("admits");
    assert_eq!(
        ExecutionReceipt::admit(
            &capped,
            capped.content_id(),
            ColorRank::Verified,
            now(),
            id("lab-b"),
        )
        .expect_err("outcome above ceiling")
        .rule(),
        ActionRule::CeilingExceeded
    );

    let other = admit(ActionKind::UqSamples, "claim-other");
    assert_eq!(
        ExecutionReceipt::admit(
            &other,
            coupon.content_id(),
            ColorRank::Estimated,
            now(),
            id("lab-b"),
        )
        .expect_err("wrong proposal id")
        .rule(),
        ActionRule::ReceiptProposalMismatch
    );
}

#[test]
fn canonical_roundtrip_preserves_identity_and_every_field() {
    let mut rich = draft(ActionKind::SensorCampaign, "claim-flux");
    rich.budgets.capabilities = BTreeSet::from([id("wind-tunnel-b"), id("pivt-rig")]);
    rich.dependencies = BTreeSet::from([admit(ActionKind::UqSamples, "claim-q").content_id()]);
    rich.exclusivity_group = Some(id("campaign-slot-3"));
    rich.correlation_group = Some(id("shared-instrument"));
    rich.expiry = Some(Expiry {
        clock: id("campaign-clock"),
        tick: 9999,
    });
    rich.ceiling = ColorRank::Validated;
    let proposal = ActionProposal::admit(rich).expect("admits");

    let bytes = proposal.canonical_bytes();
    let decoded = ActionProposal::decode(&bytes).expect("decodes");
    assert_eq!(decoded, proposal);
    assert_eq!(decoded.content_id(), proposal.content_id());
    assert_eq!(decoded.kind(), ActionKind::SensorCampaign);
    assert_eq!(decoded.ceiling(), ColorRank::Validated);
    assert_eq!(decoded.budgets().capabilities.len(), 2);
}

#[test]
fn transport_refuses_trailing_truncation_bad_magic_and_unknown_version() {
    let proposal = admit(ActionKind::UqSamples, "claim-q");
    let bytes = proposal.canonical_bytes();

    let mut trailing = bytes.clone();
    trailing.push(0);
    assert!(ActionProposal::decode(&trailing).is_err());

    for cut in [0usize, 3, 8, bytes.len() / 2, bytes.len() - 1] {
        assert!(
            ActionProposal::decode(&bytes[..cut]).is_err(),
            "truncation at {cut} must refuse"
        );
    }

    let mut magic = bytes.clone();
    magic[0] ^= 0xFF;
    let err = ActionProposal::decode(&magic).expect_err("bad magic");
    assert_eq!(err.rule_name(), "action-canonical-identity");

    let mut version = bytes;
    version[7] = 9;
    let err = ActionProposal::decode(&version).expect_err("unknown version");
    assert!(err.detail().contains("unknown schema version"));
}

#[test]
fn corrupt_mutations_never_alias_the_original_identity() {
    let proposal = admit(ActionKind::MaterialCouponTest, "claim-yield");
    let bytes = proposal.canonical_bytes();
    let original = proposal.content_id();
    for index in 0..bytes.len() {
        let mut mutated = bytes.clone();
        mutated[index] ^= 0x01;
        match ActionProposal::decode(&mutated) {
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
fn migration_golden_v1_header_and_identity_stay_decodable() {
    let proposal = admit(ActionKind::SolverTolerance, "claim-golden");
    let bytes = proposal.canonical_bytes();
    let hex: String = bytes.iter().map(|b| format!("{b:02x}")).collect();

    // Schema-v1 golden: canonical header (magic FSEA + version 1) and the
    // exact content identity of the proposal above. A change here is a WIRE
    // BREAK and needs a version bump plus a migration path, not a golden
    // edit.
    assert!(
        hex.starts_with("4653454100000001"),
        "v1 canonical header drifted: {}",
        &hex[..16.min(hex.len())]
    );
    const PINNED_CONTENT_ID: &str =
        "a2bb9c5ce99a6884347edeb60dbbc0d8df68d1e75b90f15f15b96c215c9af044";
    assert_eq!(
        proposal.content_id().to_hex(),
        PINNED_CONTENT_ID,
        "v1 content identity drifted (full canonical hex: {hex})"
    );

    let decoded = ActionProposal::decode(&bytes).expect("v1 golden decodes");
    assert_eq!(decoded.schema_version(), 1);
    assert_eq!(decoded, proposal);
}
