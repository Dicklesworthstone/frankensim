//! G0/G3 request-matrix tests for the executable E09 claim router.

use fs_govern::{
    CLAIM_ROUTER_NO_CLAIM, CLAIM_ROUTER_SCHEMA_VERSION, ChaosBasis, ClaimClass, ClaimExtent,
    ClaimRequest, ClaimRouteDecision, ClaimRouteRefusalCause, ClaimRouterError, DecisionNeed,
    DynamicsProfile, EvidenceRegime, certificate_regime, route_claim,
};

fn decision(tolerance: f64) -> DecisionNeed {
    DecisionNeed::try_new(
        "select evidence for an engineering decision",
        "qoi-unit",
        tolerance,
    )
    .expect("valid decision")
}

fn general() -> DynamicsProfile {
    DynamicsProfile::general()
}

fn chaotic(predictability_horizon: f64, unit: &str) -> DynamicsProfile {
    DynamicsProfile::new(
        false,
        true,
        false,
        ChaosBasis::try_declared(predictability_horizon, unit).expect("chaos basis"),
    )
}

fn canonical_extent(claim: ClaimClass) -> ClaimExtent {
    match claim {
        ClaimClass::RootOrEventTime
        | ClaimClass::ShortHorizonReachability
        | ClaimClass::ConservedQuantity => {
            ClaimExtent::try_finite_horizon(1.0, "seconds").expect("finite extent")
        }
        ClaimClass::LocalStability => ClaimExtent::Local,
        ClaimClass::LongHorizonMeanLoad
        | ClaimClass::BroadbandSpectrum
        | ClaimClass::ExactLongChaoticTrajectory => {
            ClaimExtent::try_long_horizon(100.0, "seconds").expect("long extent")
        }
        ClaimClass::DutyCycleReliability => {
            ClaimExtent::try_population("declared thermal duty cycles").expect("population")
        }
    }
}

fn request(
    claim: ClaimClass,
    extent: ClaimExtent,
    system: DynamicsProfile,
    tolerance: f64,
) -> ClaimRequest {
    ClaimRequest::try_new(
        format!("request/{}", claim.code()),
        claim,
        format!("fixture quantity for {}", claim.code()),
        extent,
        decision(tolerance),
        system,
        vec![
            "model version is fixed for this request".to_string(),
            "units are explicit".to_string(),
        ],
    )
    .expect("valid request")
}

#[test]
fn every_doctrine_row_routes_or_refuses_in_its_canonical_shape() {
    for claim in ClaimClass::ALL {
        let system = if claim == ClaimClass::ExactLongChaoticTrajectory {
            chaotic(10.0, "seconds")
        } else {
            general()
        };
        let route = route_claim(request(claim, canonical_extent(claim), system, 0.5));
        let row = certificate_regime(claim);
        assert_eq!(route.doctrine_row(), row);
        match claim {
            ClaimClass::ExactLongChaoticTrajectory => {
                let refusal = route.refusal().expect("exact chaotic route refuses");
                assert_eq!(
                    refusal.cause(),
                    &ClaimRouteRefusalCause::ExactLongChaoticTrajectoryHasNoUsefulRoute {
                        requested: 100.0,
                        predictability_horizon: 10.0,
                        unit: "seconds".to_string(),
                    }
                );
                assert_eq!(
                    refusal.suggested_reformulation(),
                    ClaimClass::LongHorizonMeanLoad
                );
                assert_eq!(refusal.required_evidence(), EvidenceRegime::NoUsefulBound);
            }
            _ => {
                let routed = route.routed().expect("canonical request routes");
                assert_eq!(routed.row_id(), row.id);
                assert_eq!(routed.evidence(), row.evidence);
            }
        }
        println!("{}", route.render_record());
    }
}

#[test]
fn claim_by_system_request_matrix_is_deterministic_and_fail_closed() {
    for claim in ClaimClass::ALL {
        for (system_name, system) in [
            ("general", general()),
            ("chaotic", chaotic(10.0, "seconds")),
        ] {
            let route = route_claim(request(claim, canonical_extent(claim), system, 0.5));
            if claim == ClaimClass::ExactLongChaoticTrajectory {
                let refusal = route.refusal().expect("trajectory request refuses");
                assert!(matches!(
                    (system_name, refusal.cause()),
                    (
                        "general",
                        ClaimRouteRefusalCause::MissingChaoticClassification
                    ) | (
                        "chaotic",
                        ClaimRouteRefusalCause::ExactLongChaoticTrajectoryHasNoUsefulRoute { .. }
                    )
                ));
            } else {
                assert!(
                    route.routed().is_some(),
                    "{claim:?}/{system_name} should route"
                );
            }
        }
    }
}

#[test]
fn invalid_extent_horizon_and_chaos_pairs_refuse_before_compute() {
    let wrong_extent = route_claim(request(
        ClaimClass::RootOrEventTime,
        ClaimExtent::try_long_horizon(1.0, "seconds").expect("extent"),
        general(),
        0.5,
    ));
    assert!(matches!(
        wrong_extent.refusal().map(|refusal| refusal.cause()),
        Some(ClaimRouteRefusalCause::ExtentMismatch {
            required: "finite-time-or-parameter-domain",
            found: "long-horizon",
        })
    ));

    let parameter_root = route_claim(request(
        ClaimClass::RootOrEventTime,
        ClaimExtent::try_finite_parameter_domain(2.0, "radians").expect("parameter extent"),
        chaotic(10.0, "seconds"),
        0.5,
    ));
    assert!(
        parameter_root.routed().is_some(),
        "a parameter-domain root must not be compared with a time predictability horizon"
    );

    let beyond = route_claim(request(
        ClaimClass::ShortHorizonReachability,
        ClaimExtent::try_finite_horizon(11.0, "seconds").expect("extent"),
        chaotic(10.0, "seconds"),
        0.5,
    ));
    assert!(matches!(
        beyond.refusal().map(|refusal| refusal.cause()),
        Some(ClaimRouteRefusalCause::PredictabilityHorizonExceeded { .. })
    ));

    let inside = route_claim(request(
        ClaimClass::ExactLongChaoticTrajectory,
        ClaimExtent::try_long_horizon(5.0, "seconds").expect("extent"),
        chaotic(10.0, "seconds"),
        0.5,
    ));
    let refusal = inside.refusal().expect("inside-horizon request refuses");
    assert!(matches!(
        refusal.cause(),
        ClaimRouteRefusalCause::InsidePredictabilityHorizon { .. }
    ));
    assert_eq!(
        refusal.suggested_reformulation(),
        ClaimClass::RootOrEventTime
    );

    let unit_mismatch = route_claim(request(
        ClaimClass::ExactLongChaoticTrajectory,
        ClaimExtent::try_long_horizon(100.0, "cycles").expect("extent"),
        chaotic(10.0, "seconds"),
        0.5,
    ));
    assert!(matches!(
        unit_mismatch.refusal().map(|refusal| refusal.cause()),
        Some(ClaimRouteRefusalCause::IncompatibleHorizonUnits { .. })
    ));
}

#[test]
fn assumption_order_and_duplicates_do_not_change_routing() {
    let make = |assumptions: Vec<String>| {
        ClaimRequest::try_new(
            "mean-load",
            ClaimClass::LongHorizonMeanLoad,
            "mean bearing load",
            ClaimExtent::try_long_horizon(1_000.0, "cycles").expect("extent"),
            decision(1.0),
            general(),
            assumptions,
        )
        .expect("request")
    };
    let first = make(vec![
        "stationary operating regime".to_string(),
        "declared sampling law".to_string(),
        "stationary operating regime".to_string(),
    ]);
    let second = make(vec![
        "declared sampling law".to_string(),
        "stationary operating regime".to_string(),
    ]);
    assert_eq!(first.assumptions(), second.assumptions());
    assert_eq!(route_claim(first), route_claim(second));
}

#[test]
fn relaxing_decision_tolerance_never_turns_a_route_into_a_refusal() {
    for tolerance in [0.01, 1.0, 100.0] {
        let accepted = route_claim(request(
            ClaimClass::LongHorizonMeanLoad,
            ClaimExtent::try_long_horizon(100.0, "seconds").expect("extent"),
            general(),
            tolerance,
        ));
        assert!(accepted.routed().is_some());

        let refused = route_claim(request(
            ClaimClass::ExactLongChaoticTrajectory,
            ClaimExtent::try_long_horizon(100.0, "seconds").expect("extent"),
            chaotic(10.0, "seconds"),
            tolerance,
        ));
        assert!(refused.refusal().is_some());
    }
}

#[test]
fn malformed_requests_refuse_construction_and_records_name_no_authority() {
    assert!(matches!(
        DecisionNeed::try_new("decision", "kelvin", 0.0),
        Err(ClaimRouterError::InvalidPositiveFinite {
            field: "decision.tolerance",
            ..
        })
    ));
    assert!(matches!(
        ClaimExtent::try_finite_horizon(f64::INFINITY, "seconds"),
        Err(ClaimRouterError::InvalidPositiveFinite {
            field: "finite_horizon.duration",
            ..
        })
    ));
    assert!(matches!(
        ClaimRequest::try_new(
            "empty-assumptions",
            ClaimClass::LocalStability,
            "spectral abscissa",
            ClaimExtent::Local,
            decision(1.0),
            general(),
            vec![],
        ),
        Err(ClaimRouterError::MissingAssumptions)
    ));
    assert!(matches!(
        ChaosBasis::try_probed(0.0, 10.0, "seconds"),
        Err(ClaimRouterError::InvalidPositiveFinite {
            field: "chaos.local_lyapunov_lower",
            ..
        })
    ));
    assert!(matches!(
        ClaimRequest::try_new(
            "forged-extent",
            ClaimClass::RootOrEventTime,
            "event time",
            ClaimExtent::FiniteHorizon {
                duration: f64::NAN,
                unit: "seconds".to_string(),
            },
            decision(1.0),
            general(),
            vec!["explicit assumption".to_string()],
        ),
        Err(ClaimRouterError::InvalidPositiveFinite {
            field: "finite_horizon.duration",
            ..
        })
    ));
    assert!(matches!(
        ClaimRequest::try_new(
            "forged-chaos",
            ClaimClass::ExactLongChaoticTrajectory,
            "pointwise state",
            ClaimExtent::try_long_horizon(100.0, "seconds").expect("extent"),
            decision(1.0),
            DynamicsProfile::new(
                false,
                false,
                false,
                ChaosBasis::Declared {
                    predictability_horizon: f64::NAN,
                    unit: "seconds".to_string(),
                },
            ),
            vec!["explicit assumption".to_string()],
        ),
        Err(ClaimRouterError::InvalidPositiveFinite {
            field: "chaos.predictability_horizon",
            ..
        })
    ));

    let route = route_claim(request(
        ClaimClass::LocalStability,
        ClaimExtent::Local,
        general(),
        0.5,
    ));
    let record = route.render_record();
    assert!(record.contains(&format!(
        "claim-router-schema={CLAIM_ROUTER_SCHEMA_VERSION}"
    )));
    assert!(record.contains("outcome=routed"));
    assert!(record.contains("doctrine-row=CR-04"));
    assert!(record.contains("capability=fs-spectral/residual-enclosed-spectral-service:available"));
    assert!(record.contains(CLAIM_ROUTER_NO_CLAIM));
    assert!(record.contains("doctrine-no-claim="));
    assert!(matches!(route, ClaimRouteDecision::Routed(_)));
}
