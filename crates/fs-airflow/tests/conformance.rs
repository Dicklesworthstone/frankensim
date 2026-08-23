//! G0/G3 conformance and evidence-boundary tests for the airflow rung.

use fs_airflow::{
    AirflowError, EnclosureNetwork, FanArrangement, FanBank, FanCurve, FanPoint, LeakageElement,
    LossElement, LossNetwork, LossResistance, ORIFICE_CD_SOURCE_IDENTIFIER,
    ORIFICE_PLATEAU_REYNOLDS, ORIFICE_RESISTANCE_UNCERTAINTY_REL, SHARP_EDGED_ORIFICE_CD,
    SourceProvenance, ToleranceBasis, fan_law_source, sharp_edged_orifice_loss,
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
    fan_curve_with_pressure_scale(stall_flow, 1.0, 0.08)
}

fn fan_curve_with_pressure_scale(
    stall_flow: f64,
    pressure_scale: f64,
    pressure_tolerance_rel: f64,
) -> FanCurve {
    FanCurve::new(
        "synthetic-reference-fan",
        vec![
            FanPoint::new(
                VolumetricFlowRate::new(0.00),
                Pressure::new(160.0 * pressure_scale),
            ),
            FanPoint::new(
                VolumetricFlowRate::new(0.04),
                Pressure::new(130.0 * pressure_scale),
            ),
            FanPoint::new(
                VolumetricFlowRate::new(0.08),
                Pressure::new(70.0 * pressure_scale),
            ),
            FanPoint::new(VolumetricFlowRate::new(0.12), Pressure::new(0.0)),
        ],
        synthetic_source("synthetic-fan-v1"),
        pressure_tolerance_rel,
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
    assert_eq!(
        series.equivalent_resistance().value().to_bits(),
        200.0_f64.to_bits()
    );
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
fn fan_law_source_is_clause_addressed() {
    let source = fan_law_source();
    assert_eq!(source.identifier, "AMCA-201-02-R2011:6.5.1");
    assert_eq!(
        source.citation,
        "AMCA Publication 201-02 (R2011), Fans and Systems, section 6.5.1"
    );
    assert!(
        !source.citation.contains("AMCA 210"),
        "the laboratory test-method standard is not the fan-law authority"
    );
    let card = fan_curve(0.01).model_card();
    assert!(
        card.assumptions.iter().any(|assumption| {
            assumption.contains(&source.citation) && assumption.contains(&source.identifier)
        }),
        "the consumed model card must retain both source fields"
    );
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
fn three_speed_points_follow_affinity_ordering() {
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
fn repeated_operating_point_solve_replays_the_complete_artifact() {
    let fan = FanBank::new(fan_curve(0.01), 1, FanArrangement::Series, 1.0).expect("bank");
    let system = network(180_000.0);
    let first = solve_operating_point(&fan, &system).expect("first solve");

    for replay_index in 0..4 {
        let replay = solve_operating_point(&fan, &system).expect("replay");
        assert_eq!(
            first, replay,
            "complete OperatingPoint changed on replay {replay_index}"
        );
        assert_eq!(
            first.flow.value.value().to_bits(),
            replay.flow.value.value().to_bits()
        );
        assert_eq!(
            first.pressure.value.value().to_bits(),
            replay.pressure.value.value().to_bits()
        );
        assert_eq!(
            first.nominal_root.flow.lo().to_bits(),
            replay.nominal_root.flow.lo().to_bits()
        );
        assert_eq!(
            first.nominal_root.flow.hi().to_bits(),
            replay.nominal_root.flow.hi().to_bits()
        );
    }
}

#[test]
fn pressure_estimate_uses_reachable_paired_tolerance_corners() {
    // Equal parallel branches with R=25_000 and u=0.6 have exact equivalent
    // resistance bounds 2_500 and 10_000. Zero-tolerance replicas of the
    // low-fan/low-resistance and high-fan/high-resistance corners expose their
    // certified roots independently, so the old cross-pairing is falsifiable
    // without hard-coding an output golden.
    let fan = FanBank::new(fan_curve(0.01), 1, FanArrangement::Series, 1.0).expect("bank");
    let system = EnclosureNetwork::new(
        LossNetwork::Element(loss("primary", 25_000.0, 0.6)),
        LeakageElement::new(loss("leakage", 25_000.0, 0.6)),
    );
    let point = solve_operating_point(&fan, &system).expect("solve");

    let low_corner_fan = FanBank::new(
        fan_curve_with_pressure_scale(0.01, 0.92, 0.0),
        1,
        FanArrangement::Series,
        1.0,
    )
    .expect("low-pressure fan corner");
    let low_corner_system = EnclosureNetwork::new(
        LossNetwork::Element(loss("primary", 10_000.0, 0.0)),
        LeakageElement::new(loss("leakage", 10_000.0, 0.0)),
    );
    let low_corner =
        solve_operating_point(&low_corner_fan, &low_corner_system).expect("low corner");
    let low_corner_flow = low_corner.nominal_root.flow.lo();
    let expected_low = 2_500.0 * low_corner_flow * low_corner_flow;

    let high_corner_fan = FanBank::new(
        fan_curve_with_pressure_scale(0.01, 1.08, 0.0),
        1,
        FanArrangement::Series,
        1.0,
    )
    .expect("high-pressure fan corner");
    let high_corner_system = EnclosureNetwork::new(
        LossNetwork::Element(loss("primary", 40_000.0, 0.0)),
        LeakageElement::new(loss("leakage", 40_000.0, 0.0)),
    );
    let high_corner =
        solve_operating_point(&high_corner_fan, &high_corner_system).expect("high corner");
    let high_corner_flow = high_corner.nominal_root.flow.hi();
    let expected_high = 10_000.0 * high_corner_flow * high_corner_flow;

    assert_eq!(
        point.pressure.numerical.lo.to_bits(),
        expected_low.to_bits()
    );
    assert_eq!(
        point.pressure.numerical.hi.to_bits(),
        expected_high.to_bits()
    );

    let impossible_mixed_low = 2_500.0 * point.flow.numerical.lo * point.flow.numerical.lo;
    let impossible_mixed_high = 10_000.0 * point.flow.numerical.hi * point.flow.numerical.hi;
    assert_ne!(
        (
            point.pressure.numerical.lo.to_bits(),
            point.pressure.numerical.hi.to_bits()
        ),
        (
            impossible_mixed_low.to_bits(),
            impossible_mixed_high.to_bits()
        ),
        "pressure bounds must not mix roots and resistance corners"
    );
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

#[test]
fn the_operating_point_is_bit_identical_across_worker_counts() {
    // f85xj.5.5's fourth DONE-WHEN is "network solve deterministic across
    // thread counts". The neighbouring replay test proves REPEAT identity on
    // one thread, which is a different property: it cannot observe a
    // thread-local cache, a lazily-initialised shared table, or a global the
    // solver mutates. This crate spawns no threads of its own, so the claim is
    // structurally true -- but "structurally true" is an argument, and the
    // clause asks for evidence. Running the same solve CONCURRENTLY at several
    // worker counts and requiring bitwise agreement with a single-threaded
    // reference is that evidence, and it is what would actually fail if hidden
    // shared state were ever introduced beneath `solve_operating_point`.
    //
    // Everything compared here is derived from the interval-Newton bracket, so
    // an agreement that held only in the reported scalars while the bracket
    // endpoints drifted would not pass.
    let fan = FanBank::new(fan_curve(0.01), 1, FanArrangement::Series, 1.0).expect("bank");
    let system = network(180_000.0);
    let reference = solve_operating_point(&fan, &system).expect("single-threaded reference");

    for workers in [1usize, 2, 4, 8] {
        let results: Vec<_> = std::thread::scope(|scope| {
            let handles: Vec<_> = (0..workers)
                .map(|_| {
                    scope.spawn(|| {
                        let fan = FanBank::new(fan_curve(0.01), 1, FanArrangement::Series, 1.0)
                            .expect("bank");
                        let system = network(180_000.0);
                        solve_operating_point(&fan, &system).expect("concurrent solve")
                    })
                })
                .collect();
            handles
                .into_iter()
                .map(|handle| handle.join().expect("worker did not panic"))
                .collect()
        });

        assert_eq!(results.len(), workers);
        for (worker, observed) in results.iter().enumerate() {
            assert_eq!(
                &reference, observed,
                "worker {worker} of {workers} produced a different OperatingPoint"
            );
            assert_eq!(
                reference.flow.value.value().to_bits(),
                observed.flow.value.value().to_bits(),
                "worker {worker} of {workers}: flow drifted"
            );
            assert_eq!(
                reference.pressure.value.value().to_bits(),
                observed.pressure.value.value().to_bits(),
                "worker {worker} of {workers}: pressure drifted"
            );
            assert_eq!(
                reference.nominal_root.flow.lo().to_bits(),
                observed.nominal_root.flow.lo().to_bits(),
                "worker {worker} of {workers}: certified bracket lower endpoint drifted"
            );
            assert_eq!(
                reference.nominal_root.flow.hi().to_bits(),
                observed.nominal_root.flow.hi().to_bits(),
                "worker {worker} of {workers}: certified bracket upper endpoint drifted"
            );
            assert_eq!(
                reference.nominal_root.boxes_examined, observed.nominal_root.boxes_examined,
                "worker {worker} of {workers}: interval-Newton search path differed"
            );
            assert_eq!(
                reference.branches.len(),
                observed.branches.len(),
                "worker {worker} of {workers}: branch count differed"
            );
            for (branch, (want, got)) in reference
                .branches
                .iter()
                .zip(observed.branches.iter())
                .enumerate()
            {
                assert_eq!(
                    want.path, got.path,
                    "worker {worker} of {workers}: branch {branch} traversal order differed"
                );
                assert_eq!(
                    want.flow.value.value().to_bits(),
                    got.flow.value.value().to_bits(),
                    "worker {worker} of {workers}: branch {branch} ({}) flow drifted",
                    want.path
                );
                assert_eq!(
                    want.flow.provenance, got.flow.provenance,
                    "worker {worker} of {workers}: branch {branch} provenance drifted"
                );
            }
        }
    }
}

// ---- Sharp-edged-orifice vent lowering (frankensim-frn2i slice 1) ----

#[test]
fn orifice_resistance_matches_the_closed_form_derivation() {
    let element = sharp_edged_orifice_loss("side-vent", Area::new(0.01), Density::new(1.2))
        .expect("physical vent lowers");
    // Independent operation order: zeta = 1/Cd^2 first, then zeta * rho / (2 A^2).
    let zeta = 1.0 / (SHARP_EDGED_ORIFICE_CD * SHARP_EDGED_ORIFICE_CD);
    let expected = zeta * 1.2 / (2.0 * 0.01 * 0.01);
    let got = element.resistance.value();
    assert!(
        ((got - expected) / expected).abs() < 1.0e-14,
        "resistance {got} drifted from the closed form {expected}"
    );
    assert_eq!(
        element.uncertainty_rel.to_bits(),
        ORIFICE_RESISTANCE_UNCERTAINTY_REL.to_bits()
    );
    assert_eq!(element.source.identifier, ORIFICE_CD_SOURCE_IDENTIFIER);
    assert_eq!(
        element.tolerance_basis,
        ToleranceBasis::EngineeringAllowance
    );
}

#[test]
fn orifice_resistance_uncertainty_covers_the_sourced_cd_range() {
    // The declared allowance must bracket both sourced endpoints Cd in
    // {0.60, 0.65} relative to the retained plateau constant.
    for cd in [0.60_f64, 0.65] {
        let scale = SHARP_EDGED_ORIFICE_CD / cd;
        let ratio = scale * scale;
        assert!(
            (ratio - 1.0).abs() < ORIFICE_RESISTANCE_UNCERTAINTY_REL,
            "declared allowance does not cover sourced Cd {cd}: ratio {ratio}"
        );
    }
}

#[test]
fn orifice_lowering_refuses_nonphysical_inputs() {
    for (area, density, field) in [
        (0.0, 1.2, "open area"),
        (-0.01, 1.2, "open area"),
        (f64::NAN, 1.2, "open area"),
        (0.01, 0.0, "air density"),
        (0.01, -1.2, "air density"),
        (0.01, f64::INFINITY, "air density"),
    ] {
        match sharp_edged_orifice_loss("vent", Area::new(area), Density::new(density)) {
            Err(AirflowError::InvalidOrificeInput { field: got, .. }) => {
                assert_eq!(got, field, "wrong refused field for ({area}, {density})");
            }
            other => panic!("({area}, {density}) must refuse as InvalidOrificeInput: {other:?}"),
        }
    }
    match sharp_edged_orifice_loss("", Area::new(0.01), Density::new(1.2)) {
        Err(AirflowError::EmptyElementName) => {}
        other => panic!("empty name must propagate the constructor refusal: {other:?}"),
    }
}

#[test]
fn orifice_element_carries_the_plateau_regime_card() {
    let element = sharp_edged_orifice_loss("intake-vent", Area::new(0.004), Density::new(1.18))
        .expect("physical vent lowers");
    let card = element
        .regime_audit_card()
        .expect("the orifice model must never be cardless");
    assert_eq!(card.name, "airflow.loss.intake-vent");
    assert_eq!(card.version, ORIFICE_CD_SOURCE_IDENTIFIER);
    let bounds = card.validity.bounds();
    let (lo, hi) = bounds
        .get("loss_reynolds")
        .expect("plateau domain is on the loss_reynolds axis");
    assert_eq!((*lo, *hi), ORIFICE_PLATEAU_REYNOLDS);
}

#[test]
fn vent_lowered_network_solves_a_certified_operating_point() {
    let density = Density::new(1.2);
    let primary = LossNetwork::parallel(vec![
        LossNetwork::Element(
            sharp_edged_orifice_loss("intake-vent", Area::new(0.01), density)
                .expect("physical vent"),
        ),
        LossNetwork::Element(
            sharp_edged_orifice_loss("exhaust-vent", Area::new(0.01), density)
                .expect("physical vent"),
        ),
    ])
    .expect("parallel vent network");
    let leakage = LeakageElement::new(
        sharp_edged_orifice_loss("seam-leakage", Area::new(0.002), density)
            .expect("physical leakage orifice"),
    );
    let network = EnclosureNetwork::new(primary, leakage);
    let fan = FanBank::new(fan_curve(0.02), 1, FanArrangement::Series, 1.0).expect("fan bank");
    let point = solve_operating_point(&fan, &network).expect("certified operating point");
    let root = point.nominal_root.flow;
    assert!(root.lo() > 0.0 && root.hi().is_finite() && root.lo() <= root.hi());
    let paths: Vec<&str> = point
        .branches
        .iter()
        .map(|branch| branch.path.as_str())
        .collect();
    assert!(paths.contains(&"intake-vent"));
    assert!(paths.contains(&"exhaust-vent"));
    let leakage_branch = point
        .branches
        .iter()
        .find(|branch| branch.leakage)
        .expect("mandatory leakage branch is reported");
    assert_eq!(leakage_branch.path, "seam-leakage");
    assert!(point.leakage_fraction > 0.0 && point.leakage_fraction < 1.0);
}
