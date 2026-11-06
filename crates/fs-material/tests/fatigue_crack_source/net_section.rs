//! A selected net-ligament limit is separate from the Paris/K validity bounds.
use super::*;
use fs_material::fatigue::crack_growth::{
    CentreCrackTensileError, CentreCrackTensileLimit, FiniteWidthCentreCrack, ParisLaw,
    TabulatedCrackGeometry,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::strength::IsotropicTensileLimit;

fn pressure() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}
fn guard_point(temperature: f64) -> QueryPoint {
    physical_point(temperature)
        .with_quantity("ambient-pressure", pressure(), 1e5)
        .unwrap()
}

fn fixture() -> (MaterialCard, ParisCrackGrowth, CentreCrackTensileLimit) {
    let request = request();
    let mut claims = rows(&request);
    for claim in &mut claims {
        claim.validity =
            claim
                .validity
                .clone()
                .with_quantity("ambient-pressure", pressure(), 1e5, 2e5);
    }
    claims.push(PropertyClaim {
        key: PropertyKey::new("selected_net_section_limit", Pressure::DIMS),
        value: PropertyValue::Scalar {
            value: 12e6,
            dims: Pressure::DIMS,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity("T", absolute_temperature(), 300.0, 400.0)
            .with_quantity("ambient-pressure", pressure(), 1e5, 2e5),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic selected ligament limit, not measured rupture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    });
    let card = card(claims);
    let fit_point = point(&request, 300.0)
        .with_quantity("ambient-pressure", pressure(), 1e5)
        .unwrap();
    let law = Arc::new(resolve_paris(&card, &fit_point, &request, [SINGLE; 4]).unwrap());
    let geometry = FiniteWidthCentreCrack::new(0.01)
        .unwrap()
        .tabulate(vec![0.0005, 0.001, 0.0025, 0.004], 4)
        .unwrap();
    let growth = ParisCrackGrowth::new(
        law,
        geometry,
        0.001,
        CrackGrowthLimits {
            maximum_blocks: 10,
            maximum_cycles: 1e8,
        },
    )
    .unwrap();
    let limit = resolve_limit(&card, &guard_point(300.0));
    let guard = CentreCrackTensileLimit::try_new(&growth, limit).unwrap();
    (card, growth, guard)
}

fn resolve_limit(card: &MaterialCard, point: &QueryPoint) -> IsotropicTensileLimit {
    IsotropicTensileLimit::resolve(
        card,
        point,
        "selected_net_section_limit",
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn net_section_uses_physical_half_length_and_gross_area_stress() {
    let (card, _, guard) = fixture();
    let assessed = guard.assess(0.0025, 6e6).unwrap();
    assert_eq!(assessed.half_length_m(), 0.0025);
    assert_eq!(assessed.gross_stress_pa(), 6e6);
    assert_eq!(assessed.ligament_fraction(), 0.5);
    assert_eq!(assessed.net_stress_pa(), 12e6);
    assert_eq!(assessed.utilization(), 1.0);
    assert!(assessed.reaches_limit());
    assert_eq!(
        guard.check(0.0025, 6e6),
        Err(CentreCrackTensileError::LimitReached(assessed))
    );
    assert_eq!(guard.assess(0.0025, 0.0).unwrap().utilization(), 0.0);
    assert_eq!(guard.plate_width_m(), 0.01);
    assert_eq!(
        guard.limit().resolved().card_identity(),
        card.content_hash()
    );
    assert_eq!(
        guard
            .geometry()
            .finite_width_centre_crack()
            .unwrap()
            .plate_width_m(),
        0.01
    );
    assert!(matches!(
        guard.assess(0.0025, -1.0),
        Err(CentreCrackTensileError::Growth(
            CrackGrowthError::StressOutsideDomain { .. }
        ))
    ));
    for stress in [f64::MAX, f64::from_bits(1)] {
        assert!(matches!(
            guard.assess(0.0025, stress),
            Err(CentreCrackTensileError::Arithmetic)
        ));
    }
}

#[test]
fn guard_requires_exact_typed_geometry_card_and_fixed_operating_query() {
    let (card, growth, guard) = fixture();
    let foreign = MaterialCard::supersede(&card, card.claims().clone(), vec![]).unwrap();
    let extra = guard_point(300.0)
        .with_quantity(
            "extra-condition",
            QuantitySpec::dimensional(Dims::NONE),
            0.2,
        )
        .unwrap();
    for limit in [
        resolve_limit(&foreign, &guard_point(300.0)),
        resolve_limit(&card, &guard_point(310.0)),
        resolve_limit(&card, &extra),
    ] {
        assert!(matches!(
            CentreCrackTensileLimit::try_new(&growth, limit),
            Err(CentreCrackTensileError::Invalid(_))
        ));
    }
    // Even identical Y points and copied prose do not carry typed geometry.
    let anonymous = TabulatedCrackGeometry::new(
        guard.geometry().context().clone(),
        guard.geometry().points().to_vec(),
        4,
    )
    .unwrap();
    let generic = ParisCrackGrowth::new(
        Arc::new(growth.law().clone()),
        anonymous,
        0.001,
        growth.limits(),
    )
    .unwrap();
    assert!(matches!(
        CentreCrackTensileLimit::try_new(&generic, guard.limit().clone()),
        Err(CentreCrackTensileError::Invalid(_))
    ));
    let constant = ParisCrackGrowth::new(
        Arc::new(growth.law().clone()),
        ConstantGeometryCrack {
            geometry_factor: 1.0,
            minimum_length_m: 0.0005,
            maximum_length_m: 0.004,
        },
        0.001,
        growth.limits(),
    )
    .unwrap();
    assert!(matches!(
        CentreCrackTensileLimit::try_new(&constant, guard.limit().clone()),
        Err(CentreCrackTensileError::Invalid(_))
    ));
    let declared = ParisCrackGrowth::new(
        Arc::new(
            ParisLaw::new(
                growth.law().scope().clone(),
                growth.law().parameters(),
                growth.law().domain(),
            )
            .unwrap(),
        ),
        growth.geometry().clone(),
        0.001,
        growth.limits(),
    )
    .unwrap();
    assert!(matches!(
        CentreCrackTensileLimit::try_new(&declared, guard.limit().clone()),
        Err(CentreCrackTensileError::Invalid(_))
    ));
    let mut other_width = ParisCrackGrowth::new(
        Arc::new(growth.law().clone()),
        FiniteWidthCentreCrack::new(0.02)
            .unwrap()
            .tabulate(vec![0.0005, 0.004], 2)
            .unwrap(),
        0.001,
        growth.limits(),
    )
    .unwrap();
    let before = other_width.clone();
    assert!(matches!(
        guard.advance(
            &mut other_width,
            StressCycle {
                minimum_pa: 2e6,
                maximum_pa: 8e6,
                cycles: 1.0
            }
        ),
        Err(CentreCrackTensileError::Invalid(_))
    ));
    assert_eq!(other_width, before);
}

#[test]
fn analytical_limit_forecast_and_guarded_growth_share_the_original_integrator() {
    let (_, mut growth, guard) = fixture();
    let before = growth.clone();
    let forecast = guard.cycles_to_limit(&growth, 2e6, 8e6).unwrap();
    close(forecast.target_length_m, (0.01 / 2.0) * (1.0 - 8e6 / 12e6));
    assert_eq!(growth, before);
    let mut replay = growth.clone();
    replay
        .advance(StressCycle {
            minimum_pa: 2e6,
            maximum_pa: 8e6,
            cycles: forecast.additional_cycles,
        })
        .unwrap();
    close(replay.state().crack_length_m, forecast.target_length_m);
    let refused = guard.advance(
        &mut growth,
        StressCycle {
            minimum_pa: 2e6,
            maximum_pa: 8e6,
            cycles: forecast.additional_cycles * 1.01,
        },
    );
    assert!(matches!(
        refused,
        Err(CentreCrackTensileError::LimitReached(_))
    ));
    assert_eq!(growth, before);
    let accepted = guard
        .advance(
            &mut growth,
            StressCycle {
                minimum_pa: 2e6,
                maximum_pa: 8e6,
                cycles: forecast.additional_cycles * 0.99,
            },
        )
        .unwrap();
    assert!(accepted.after.crack_length_m < forecast.target_length_m);
    assert!(
        !guard
            .assess(accepted.after.crack_length_m, 8e6)
            .unwrap()
            .reaches_limit()
    );
}

#[test]
fn unsupported_targets_static_overloads_and_cancellation_do_not_advance() {
    let (_, mut growth, guard) = fixture();
    let before = growth.clone();
    assert!(matches!(
        guard.cycles_to_limit(&growth, 0.1e6, 0.5e6),
        Err(CentreCrackTensileError::Growth(
            CrackGrowthError::CrackLengthOutsideDomain { .. }
        ))
    ));
    assert!(matches!(
        guard.cycles_to_limit(&growth, 8e6, 8e6),
        Err(CentreCrackTensileError::Growth(CrackGrowthError::Invalid(
            _
        )))
    ));
    assert!(matches!(
        guard.advance(
            &mut growth,
            StressCycle {
                minimum_pa: 10e6,
                maximum_pa: 10e6,
                cycles: 1.0,
            }
        ),
        Err(CentreCrackTensileError::LimitReached(_))
    ));
    assert!(matches!(
        guard.advance_with_cancel(
            &mut growth,
            StressCycle {
                minimum_pa: 2e6,
                maximum_pa: 8e6,
                cycles: 1.0,
            },
            &mut || true
        ),
        Err(CentreCrackTensileError::Growth(CrackGrowthError::Cancelled))
    ));
    assert!(matches!(
        guard.cycles_to_limit_with_cancel(&growth, 2e6, 8e6, &mut || true),
        Err(CentreCrackTensileError::Growth(CrackGrowthError::Cancelled))
    ));
    assert_eq!(growth, before);
}
