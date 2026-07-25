//! G0/G3 conformance and evidence-boundary tests for the airflow rung.

use fs_airflow::{
    AirflowError, EnclosureNetwork, FanArrangement, FanBank, FanCurve, FanPoint, LeakageElement,
    LossElement, LossNetwork, LossResistance, SourceProvenance, ToleranceBasis,
    solve_operating_point,
};
use fs_evidence::{NumericalKind, ProvenanceHash, ValidityDomain};
use fs_qty::{Area, Density, DynViscosity, Length, Pressure, VolumetricFlowRate};

fn synthetic_source(id: &str) -> SourceProvenance {
    SourceProvenance::new(
        "Synthetic G0 fixture; not manufacturer performance data",
        id,
    )
}

fn fan_curve(stall_flow: f64) -> FanCurve {
    FanCurve::new(
        "synthetic-reference-fan",
        vec![
            FanPoint::new(VolumetricFlowRate::new(0.00), Pressure::new(160.0)),
            FanPoint::new(VolumetricFlowRate::new(0.04), Pressure::new(130.0)),
            FanPoint::new(VolumetricFlowRate::new(0.08), Pressure::new(70.0)),
            FanPoint::new(VolumetricFlowRate::new(0.12), Pressure::new(0.0)),
        ],
        synthetic_source("synthetic-fan-v1"),
        0.08,
        ToleranceBasis::EngineeringAllowance,
        VolumetricFlowRate::new(stall_flow),
        (0.7, 1.3),
    )
    .expect("valid synthetic fan")
}

fn loss(name: &str, resistance: f64, uncertainty: f64) -> LossElement {
    LossElement::new(
        name,
        LossResistance::new(resistance),
        uncertainty,
        synthetic_source(&format!("synthetic-loss-{name}")),
        ToleranceBasis::EngineeringAllowance,
    )
    .expect("valid synthetic loss")
}

fn network(leakage_resistance: f64) -> EnclosureNetwork {
    let intake = LossNetwork::parallel(vec![
        LossNetwork::Element(loss("left-vent", 45_000.0, 0.10)),
        LossNetwork::Element(loss("right-vent", 45_000.0, 0.10)),
    ])
    .expect("parallel vents");
    let primary = LossNetwork::series(vec![
        intake,
        LossNetwork::Element(loss("heatsink-channel", 30_000.0, 0.12)),
        LossNetwork::Element(loss("outlet", 12_000.0, 0.08)),
    ])
    .expect("series path");
    EnclosureNetwork::new(
        primary,
        LeakageElement::new(loss("case-leakage", leakage_resistance, 0.25)),
    )
}

#[test]
fn interpolation_and_monotone_admission_are_explicit() {
    let bank =
        FanBank::new(fan_curve(0.01), 1, FanArrangement::Series, 1.0).expect("admissible bank");
    let pressure = bank
        .pressure_at(VolumetricFlowRate::new(0.06))
        .expect("inside curve")
        .value();
    assert!((pressure - 100.0).abs() < 1.0e-10, "{pressure}");

    let error = FanCurve::new(
        "bad",
        vec![
            FanPoint::new(VolumetricFlowRate::new(0.02), Pressure::new(80.0)),
            FanPoint::new(VolumetricFlowRate::new(0.01), Pressure::new(70.0)),
        ],
        synthetic_source("bad"),
        0.0,
        ToleranceBasis::Analytic,
        VolumetricFlowRate::new(0.01),
        (1.0, 1.0),
    )
    .expect_err("non-monotone data must refuse");
    assert!(matches!(error, AirflowError::NonMonotoneFlow { .. }));
}

#[test]
fn quadratic_series_and_parallel_composition_obey_g0_identities() {
    let a = LossNetwork::Element(loss("a", 100.0, 0.0));
    let b = LossNetwork::Element(loss("b", 100.0, 0.0));
    let series = LossNetwork::series(vec![a.clone(), b.clone()]).expect("series");
    let parallel = LossNetwork::parallel(vec![a, b]).expect("parallel");
    assert_eq!(series.equivalent_resistance().value(), 200.0);
    assert!((parallel.equivalent_resistance().value() - 25.0).abs() < 1.0e-12);
}

#[test]
fn identical_fans_obey_series_pressure_and_parallel_flow_laws() {
    let single = FanBank::new(fan_curve(0.01), 1, FanArrangement::Series, 1.0).expect("single");
    let series = FanBank::new(fan_curve(0.01), 2, FanArrangement::Series, 1.0).expect("series");
    let parallel =
        FanBank::new(fan_curve(0.01), 2, FanArrangement::Parallel, 1.0).expect("parallel");
    let single_pressure = single
        .pressure_at(VolumetricFlowRate::new(0.06))
        .expect("single point")
        .value();
    let series_pressure = series
        .pressure_at(VolumetricFlowRate::new(0.06))
        .expect("series point")
        .value();
    let parallel_pressure = parallel
        .pressure_at(VolumetricFlowRate::new(0.12))
        .expect("parallel point")
        .value();
    assert!((series_pressure - 2.0 * single_pressure).abs() < 1.0e-9);
    assert!((parallel_pressure - single_pressure).abs() < 1.0e-9);
}

#[test]
fn operating_point_has_unique_sign_changing_nominal_bracket() {
    let fan = FanBank::new(fan_curve(0.01), 1, FanArrangement::Series, 1.0).expect("bank");
    let system = network(180_000.0);
    let point = solve_operating_point(&fan, &system).expect("certified root");
    let bracket = point.nominal_root.flow;
    let resistance = system.equivalent_resistance().value();
    let residual = |q: f64| {
        fan.pressure_at(VolumetricFlowRate::new(q))
            .expect("bracket inside fan curve")
            .value()
            - resistance * q * q
    };
    assert!(residual(bracket.lo()) >= -1.0e-8);
    assert!(residual(bracket.hi()) <= 1.0e-8);
    assert_eq!(point.flow.numerical.kind, NumericalKind::Estimate);
    assert!(point.flow.numerical.lo <= point.flow.value.value());
    assert!(point.flow.value.value() <= point.flow.numerical.hi);
}

#[test]
fn declared_stall_region_refuses() {
    let fan = FanBank::new(fan_curve(0.07), 1, FanArrangement::Series, 1.0).expect("bank");
    let high_resistance_network = network(5.0e8);
    let error = solve_operating_point(&fan, &high_resistance_network)
        .expect_err("intersection below stall boundary must refuse");
    assert!(matches!(error, AirflowError::StallRegion { .. }), "{error}");
}

#[test]
fn three_speed_points_follow_affinity_and_solve_deterministically() {
    let system = network(180_000.0);
    let mut flows = Vec::new();
    for speed in [0.8, 1.0, 1.2] {
        let fan = FanBank::new(fan_curve(0.01), 1, FanArrangement::Series, speed)
            .expect("speed inside declared range");
        flows.push(
            solve_operating_point(&fan, &system)
                .expect("operating point")
                .flow
                .value
                .value(),
        );
    }
    assert!(flows[0] < flows[1] && flows[1] < flows[2], "{flows:?}");
}

#[test]
fn leakage_sensitivity_and_branch_balance_are_visible() {
    let fan = FanBank::new(fan_curve(0.01), 1, FanArrangement::Series, 1.0).expect("bank");
    let leaky = solve_operating_point(&fan, &network(80_000.0)).expect("leaky solve");
    let tight = solve_operating_point(&fan, &network(800_000.0)).expect("tight solve");
    assert!(leaky.leakage_fraction > tight.leakage_fraction);
    let branch_sum: f64 = tight
        .branches
        .iter()
        .filter(|branch| branch.path != "heatsink-channel" && branch.path != "outlet")
        .map(|branch| branch.flow.value.value())
        .sum();
    assert!((branch_sum - tight.flow.value.value()).abs() < 1.0e-10);
    assert!(tight.branches.iter().any(|branch| branch.leakage));
}

#[test]
fn branch_flow_hands_typed_velocity_and_reynolds_to_convection() {
    let fan = FanBank::new(fan_curve(0.01), 1, FanArrangement::Series, 1.0).expect("bank");
    let point = solve_operating_point(&fan, &network(180_000.0)).expect("solve");
    let handoff = point
        .correlation_handoff(
            "heatsink-channel",
            Area::new(0.012),
            Density::new(1.18),
            DynViscosity::new(1.85e-5),
            Length::new(0.008),
            0.71,
        )
        .expect("typed handoff");
    let expected_re = 1.18 * handoff.velocity.value.value() * 0.008 / 1.85e-5;
    assert!((handoff.reynolds - expected_re).abs() < 1.0e-10);
    assert_eq!(handoff.velocity.numerical.kind, NumericalKind::Estimate);
    assert_eq!(handoff.velocity.model, handoff.branch_flow.model);
}

#[test]
fn operating_identity_binds_uncertainty_authority() {
    let fan = FanBank::new(fan_curve(0.01), 1, FanArrangement::Series, 1.0).expect("bank");
    let make_network = |uncertainty| {
        EnclosureNetwork::new(
            LossNetwork::Element(loss("primary", 55_000.0, uncertainty)),
            LeakageElement::new(loss("leakage", 180_000.0, 0.25)),
        )
    };
    let first = solve_operating_point(&fan, &make_network(0.05)).expect("first solve");
    let second = solve_operating_point(&fan, &make_network(0.20)).expect("second solve");

    assert_eq!(first.flow.value, second.flow.value);
    assert_ne!(first.flow.provenance, second.flow.provenance);
}

/// One loss element, optionally carrying a declared regime domain.
fn loss_with_validity(name: &str, axis: Option<(&str, f64, f64)>) -> LossElement {
    let element = loss(name, 55_000.0, 0.10);
    match axis {
        None => element,
        Some((axis, lo, hi)) => element
            .with_regime_validity(ValidityDomain::unconstrained().with(axis, lo, hi))
            .expect("declared domain admits"),
    }
}

fn identity_of(element: LossElement) -> ProvenanceHash {
    let fan = FanBank::new(fan_curve(0.01), 1, FanArrangement::Series, 1.0).expect("bank");
    let network = EnclosureNetwork::new(
        LossNetwork::Element(element),
        LeakageElement::new(loss("leakage", 180_000.0, 0.25)),
    );
    solve_operating_point(&fan, &network)
        .expect("solve")
        .flow
        .provenance
}

#[test]
fn a_declared_regime_domain_changes_the_operating_identity() {
    // The falsifier for bead frankensim-yq435, kept as a regression.
    //
    // `regime_validity` was omitted from the provenance encoding, so two
    // networks differing ONLY in whether a loss element carried a validated
    // operating domain hashed identically. That field is the sole input to
    // regime_audit_cards, so those two runs reach a reviewer as an ADMITTED
    // result and a DEMOTED one — the exact pair a provenance identity exists
    // to separate. A hash that cannot tell them apart is a false certificate,
    // which is worse than no certificate at all.
    let cardless = identity_of(loss_with_validity("heatsink-channel", None));
    let carded = identity_of(loss_with_validity(
        "heatsink-channel",
        Some(("loss_reynolds", 2_000.0, 80_000.0)),
    ));
    assert_ne!(
        cardless, carded,
        "a cardless element and a validity-carded one must not share an identity"
    );
}

#[test]
fn moving_a_validity_bound_changes_the_operating_identity() {
    // Presence alone is not enough: the DOMAIN decides whether an operating
    // point falls inside its card, so an element valid to Re 80 000 and one
    // valid to Re 40 000 are different declarations about the same hardware.
    let wide = identity_of(loss_with_validity(
        "heatsink-channel",
        Some(("loss_reynolds", 2_000.0, 80_000.0)),
    ));
    let narrow = identity_of(loss_with_validity(
        "heatsink-channel",
        Some(("loss_reynolds", 2_000.0, 40_000.0)),
    ));
    let renamed_axis = identity_of(loss_with_validity(
        "heatsink-channel",
        Some(("channel_reynolds", 2_000.0, 80_000.0)),
    ));
    assert_ne!(wide, narrow, "the upper bound must reach the identity");
    assert_ne!(
        wide, renamed_axis,
        "the axis NAME must reach the identity; two axes with equal numbers \
         constrain different physics"
    );
}

#[test]
fn loss_elements_that_compare_unequal_never_share_an_identity() {
    // The general property whose violation was the bug. `LossElement` derives
    // PartialEq INCLUDING regime_validity, so before the fix two elements
    // could compare unequal while hashing equal — and anything deduplicating
    // or caching on the hash would have collapsed them.
    //
    // The property is total over constructible elements: every f64 that
    // reaches the encoding is validated finite (`with_regime_validity` refuses
    // a non-finite bound, `LossElement::new` refuses a non-positive
    // resistance), so there is no NaN case where PartialEq is non-reflexive
    // and the two notions could legitimately disagree.
    let variants = [
        loss_with_validity("heatsink-channel", None),
        loss_with_validity(
            "heatsink-channel",
            Some(("loss_reynolds", 2_000.0, 80_000.0)),
        ),
        loss_with_validity(
            "heatsink-channel",
            Some(("loss_reynolds", 2_000.0, 40_000.0)),
        ),
        loss_with_validity("heatsink-channel", Some(("loss_reynolds", 500.0, 80_000.0))),
        loss_with_validity(
            "heatsink-channel",
            Some(("channel_reynolds", 2_000.0, 80_000.0)),
        ),
    ];
    for (i, left) in variants.iter().enumerate() {
        for (j, right) in variants.iter().enumerate() {
            let same_identity = identity_of(left.clone()) == identity_of(right.clone());
            assert_eq!(
                left == right,
                same_identity,
                "variants {i} and {j}: equality and identity disagree"
            );
        }
    }
}

#[test]
fn an_element_with_a_declared_domain_still_solves_to_the_same_numbers() {
    // The declaration is an EVIDENCE fact, not a physics one: adding a card
    // must move the identity and nothing else. If it also moved the answer,
    // the two runs above would differ for a second reason and the regression
    // test would no longer isolate the provenance defect.
    let fan = FanBank::new(fan_curve(0.01), 1, FanArrangement::Series, 1.0).expect("bank");
    let solve = |element: LossElement| {
        let network = EnclosureNetwork::new(
            LossNetwork::Element(element),
            LeakageElement::new(loss("leakage", 180_000.0, 0.25)),
        );
        solve_operating_point(&fan, &network).expect("solve")
    };
    let cardless = solve(loss_with_validity("heatsink-channel", None));
    let carded = solve(loss_with_validity(
        "heatsink-channel",
        Some(("loss_reynolds", 2_000.0, 80_000.0)),
    ));
    assert_eq!(
        cardless.flow.value.value().to_bits(),
        carded.flow.value.value().to_bits()
    );
    assert_eq!(
        cardless.pressure.value.value().to_bits(),
        carded.pressure.value.value().to_bits()
    );
    assert_ne!(cardless.flow.provenance, carded.flow.provenance);
}
