//! Original material cards and two-temperature interface claims, with authored
//! analytic values. These fixtures establish source/model behavior, not data quality.
use super::*;
use crate::phase::{EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_qty::Pressure;
use fs_thermal::crystallization::CrystallizationLaw;

fn materials() -> [MaterialCard; 2] {
    ["authored parent solid", "authored melting specimen"].map(|name| {
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "declared parent/crystal pair".into(),
                process: "analytic contact regression".into(),
                revision: 0,
            },
            ClaimSet::new(),
            vec![],
        )
        .unwrap()
    })
}

fn charts(cards: &[MaterialCard; 2]) -> [[EquilibriumEnthalpyPhaseCurve; 2]; 2] {
    std::array::from_fn(|side| {
        std::array::from_fn(|phase| {
            let parent = phase == 0;
            EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
                cards[side].content_hash(),
                if side == 1 && parent {
                    SolidLiquidPhase::Liquid
                } else {
                    SolidLiquidPhase::Solid
                },
                if parent {
                    if side == 0 { 100_000.0 } else { 80_000.0 }
                } else {
                    0.0
                },
                &[280.0, 550.0].map(|temperature_k| HeatCapacityKnot {
                    temperature_k,
                    specific_heat_capacity_j_kg_k: if side == 0 { 2000.0 } else { 1000.0 },
                    bulk_density_kg_m3: if side == 1 && parent { 800.0 } else { 1000.0 },
                }),
                1e-6,
            )
            .unwrap()
        })
    })
}

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn context() -> QueryPoint {
    // Deliberately different from the parcels: preparation must query THEIR T.
    QueryPoint::new()
        .with_quantity("T_a", quantity(QuantityKind::AbsoluteTemperature), 300.0)
        .unwrap()
        .with_quantity("T_b", quantity(QuantityKind::AbsoluteTemperature), 300.0)
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
        .unwrap()
        .with_quantity("gap", QuantitySpec::dimensional(Length::DIMS), 1e-5)
        .unwrap()
}

fn original_resistance(ta: f64, tb: f64) -> f64 {
    1e-4 * (1.0 + 0.001 * (ta - 350.0) + 0.0005 * (tb - 350.0))
}
fn source(cards: &[MaterialCard; 2], zero_resistance: bool) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    // Nonoverlapping claims make source substitution possible to detect: the
    // prepared high-T claim must refuse instead of switching to the other one.
    for (range, original) in [([390.0, 420.0], true), ([280.0, 389.0], false)] {
        let samples = [
            [range[0], 300.0],
            [range[1], 300.0],
            [range[0], 420.0],
            [range[1], 420.0],
        ]
        .map(|ts| {
            (
                ts,
                if original {
                    original_resistance(ts[0], ts[1])
                } else {
                    2e-4
                },
            )
        });
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(
                    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
                ),
                value: if zero_resistance && original {
                    PropertyValue::Scalar {
                        value: 0.0,
                        dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
                    }
                } else {
                    PropertyValue::Surface2D {
                        axes: ["T_a".into(), "T_b".into()],
                        axis_dims: [Temperature::DIMS; 2],
                        samples: samples.to_vec(),
                        triangles: vec![[0, 1, 2], [1, 3, 2]],
                        dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
                    }
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T_a",
                        quantity(QuantityKind::AbsoluteTemperature),
                        range[0],
                        range[1],
                    )
                    .with_quantity(
                        "T_b",
                        quantity(QuantityKind::AbsoluteTemperature),
                        300.0,
                        420.0,
                    )
                    .with_quantity(
                        "pressure",
                        quantity(QuantityKind::Pressure),
                        90_000.0,
                        110_000.0,
                    )
                    .with_quantity("gap", QuantitySpec::dimensional(Length::DIMS), 0.0, 1e-4),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if zero_resistance && original {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                },
                observations: vec![],
                provenance: Provenance {
                    source: "authored two-temperature contact resistance".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |side: usize| SurfaceSpec {
        material: cards[side].id().clone(),
        texture_frame: format!("original ordered side {side}"),
    };
    InterfaceSystemCard::assemble(
        surface(0),
        surface(1),
        SystemContext {
            medium: "explicit contact layer".into(),
            third_body: None,
            environment: "fixed area, pressure and separation".into(),
            history: "no evolving interface history".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn point<'a>(
    interface: &'a InterfaceSystemCard,
    cards: &'a [MaterialCard; 2],
    context: &'a QueryPoint,
) -> ThermalInterfacePoint<'a> {
    ThermalInterfacePoint {
        interface,
        materials: cards.each_ref(),
        context,
        temperature_axes: ["T_a", "T_b"],
        liquid_fraction_axes: [None, None],
        mean_temperature_axis: None,
    }
}
fn prepare<'a>(
    point: ThermalInterfacePoint<'a>,
    pairs: [&CrystallizationCaloricPair<'_>; 2],
    parcels: [&IsobaricCrystallizationParcel<'_, '_>; 2],
) -> Result<CrystallizationInterfaceContact<'a>, CrystallizationInterfaceError> {
    CrystallizationInterfaceContact::new(
        point,
        pairs,
        parcels,
        0.01,
        Some("pressure"),
        Some("gap"),
        MaterialPropertySelection::SingleClaimOnly,
    )
}
fn no_growth(_: f64) -> Result<f64, CrystallizationError> {
    Ok(0.0)
}
fn growth(t: f64) -> Result<f64, CrystallizationError> {
    Ok(0.01 * (350.0 - t).max(0.0))
}
fn side_controls() -> IsobaricExchangeSideControls {
    let mut controls = IsobaricExchangeSideControls::default();
    controls.caloric.energy_tolerance_j = 1e-6;
    controls.adaptivity.maximum_dt_s = 0.05;
    controls
}
fn exchange_controls(upper: f64) -> IsobaricExchangeControls {
    IsobaricExchangeControls {
        heat_bracket_j: [0.0, upper],
        maximum_evaluations: 64,
        residual_tolerance_j: 1e-6,
        energy_tolerance_j: 1e-5,
    }
}
fn parcels<'a, 's>(
    laws: &'a [CrystallizationLaw<'s>; 2],
) -> [IsobaricCrystallizationParcel<'a, 's>; 2] {
    let controls = side_controls().caloric;
    [
        IsobaricCrystallizationParcel::new(
            &laws[0],
            laws[0].state_at_temperature(400.0, 0.0).unwrap(),
            1.0,
            0.001,
            1e6,
            controls,
            || false,
        )
        .unwrap(),
        IsobaricCrystallizationParcel::new(
            &laws[1],
            laws[1].state_at_temperature(350.0, 1.0).unwrap(),
            0.5,
            0.0005,
            2e6,
            controls,
            || false,
        )
        .unwrap()
        .with_equilibrium_remelting(laws[1].equivalent_age_remelting().unwrap())
        .unwrap(),
    ]
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_original_two_temperature_resistance_drives_remelting_and_retains_the_accepted_receipt() {
    let cards = materials();
    let curves = charts(&cards);
    let pairs = curves.each_ref().map(|branches| {
        CrystallizationCaloricPair::new(&branches[0], &branches[1], [280.0, 550.0]).unwrap()
    });
    let laws = [
        pairs[0]
            .law(&no_growth, ContentHash([0xc1; 32]), 1.0, 0.8, 500.0)
            .unwrap(),
        pairs[1]
            .law(&growth, ContentHash([0xc2; 32]), 1.0, 0.8, 350.0)
            .unwrap(),
    ];
    let bodies = parcels(&laws);
    let interface = source(&cards, false);
    let context = context();
    let contact = prepare(
        point(&interface, &cards, &context),
        pairs.each_ref(),
        bodies.each_ref(),
    )
    .unwrap();
    assert_eq!(
        contact.caloric_pair_identities(),
        pairs.each_ref().map(|p| p.identity())
    );
    assert_eq!(
        contact.original_source().state_coordinate("T_a"),
        Some(400.0)
    );
    assert_eq!(
        contact.original_source().state_coordinate("T_b"),
        Some(350.0)
    );
    let original = contact.original_source().clone();
    let controls = side_controls();
    let mut checkpoints = 0;
    let (cooled, melted, report) = contact
        .exchange_heat(
            &bodies[0],
            &bodies[1],
            1.0,
            controls,
            controls,
            exchange_controls(8000.0),
            || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    // R depends on the solved donor T, so Q solves 5e-7 Q² - 1.1 Q + 5000=0.
    // This stable analytic root differs from freezing G at the accepted input.
    let heat = 10_000.0 / (1.1 + 1.2_f64.sqrt());
    close(report.thermal.heat_left_to_right_j, heat, 2e-6);
    close(cooled.state().temperature_k(), 400.0 - heat / 2000.0, 2e-9);
    close(melted.state().temperature_k(), 350.0, 1e-9);
    close(
        melted.state().crystalline_mass_fraction(),
        bodies[1].state().crystalline_mass_fraction() - heat / 40_000.0,
        1e-10,
    );
    let expected_r = original_resistance(
        cooled.state().temperature_k(),
        melted.state().temperature_k(),
    );
    close(report.interface.resistance_m2_k_w, expected_r, 1e-18);
    close(report.interface.conductance_w_k, 0.01 / expected_r, 1e-10);
    assert_eq!(
        report.thermal.conductance_w_k,
        report.interface.conductance_w_k
    );
    let accepted = &report.interface.source;
    assert_eq!(
        accepted.state_coordinate("T_a"),
        Some(cooled.state().temperature_k())
    );
    assert_eq!(
        accepted.state_coordinate("T_b"),
        Some(melted.state().temperature_k())
    );
    assert_eq!(accepted.state_coordinate("pressure"), Some(100_000.0));
    assert_eq!(accepted.state_coordinate("gap"), Some(1e-5));
    assert_eq!(
        accepted.properties()[0].answer().receipt.selected,
        original.properties()[0].answer().receipt.selected
    );
    assert_eq!(accepted.surface_materials(), original.surface_materials());
    assert_eq!(contact.original_source(), &original);
    close(
        report.thermal.combined_internal_energy_change_j + report.thermal.combined_boundary_work_j,
        0.0,
        1e-6,
    );
    assert!(report.thermal.right.boundary_work_j > 0.0);
    // Both temperature coordinates are live, even away from a phase plateau.
    let (right_cold, _) = bodies[1]
        .apply_heat(-5000.0, controls.caloric, || false)
        .unwrap();
    close(
        contact
            .evaluate(&bodies[0], &right_cold)
            .unwrap()
            .resistance_m2_k_w,
        original_resistance(400.0, 340.0),
        1e-18,
    );

    let mut polls = 0;
    let cancelled = contact.exchange_heat(
        &bodies[0],
        &bodies[1],
        1.0,
        controls,
        controls,
        exchange_controls(8000.0),
        || {
            polls += 1;
            polls == checkpoints
        },
    );
    assert!(matches!(
        cancelled,
        Err(CrystallizationInterfaceError::Thermal(
            CrystallizationError::Cancelled
        ))
    ));
    assert_eq!(bodies[0].time_s(), 0.0);
    assert_eq!(bodies[1].time_s(), 0.0);
    assert_eq!(contact.original_source(), &original);
}

#[test]
fn g4_order_phase_and_source_support_refuse_without_substitution_or_partial_publication() {
    let cards = materials();
    let curves = charts(&cards);
    let pairs = curves.each_ref().map(|branches| {
        CrystallizationCaloricPair::new(&branches[0], &branches[1], [280.0, 550.0]).unwrap()
    });
    let laws = [
        pairs[0]
            .law(&no_growth, ContentHash([0xc1; 32]), 1.0, 0.8, 500.0)
            .unwrap(),
        pairs[1]
            .law(&growth, ContentHash([0xc2; 32]), 1.0, 0.8, 350.0)
            .unwrap(),
    ];
    let bodies = parcels(&laws);
    let interface = source(&cards, false);
    let context = context();
    let point = point(&interface, &cards, &context);
    let contact = prepare(point, pairs.each_ref(), bodies.each_ref()).unwrap();
    let original = contact.original_source().clone();
    let original_states = bodies.each_ref().map(|body| body.state());
    assert!(
        prepare(
            ThermalInterfacePoint {
                materials: [&cards[1], &cards[0]],
                ..point
            },
            pairs.each_ref(),
            bodies.each_ref()
        )
        .is_err()
    );
    assert!(contact.evaluate(&bodies[1], &bodies[0]).is_err());
    assert!(matches!(
        prepare(
            ThermalInterfacePoint {
                liquid_fraction_axes: [Some("X"), None],
                ..point
            },
            pairs.each_ref(),
            bodies.each_ref()
        ),
        Err(CrystallizationInterfaceError::Invalid(
            "kinetic contact does not own equilibrium liquid-fraction axes"
        ))
    ));
    let other_pressure = IsobaricCrystallizationParcel::new(
        &laws[0],
        bodies[0].state(),
        1.0,
        0.001,
        1.1e6,
        side_controls().caloric,
        || false,
    )
    .unwrap();
    assert!(contact.evaluate(&other_pressure, &bodies[1]).is_err());
    let different_law = pairs[1]
        .law(&growth, ContentHash([0xcc; 32]), 1.0, 0.8, 350.0)
        .unwrap();
    let different_body = IsobaricCrystallizationParcel::new(
        &different_law,
        different_law.state_at_temperature(350.0, 1.0).unwrap(),
        0.5,
        0.0005,
        2e6,
        side_controls().caloric,
        || false,
    )
    .unwrap();
    assert!(contact.evaluate(&bodies[0], &different_body).is_err());

    // Both material paths remain valid, but the high-T INTERFACE source does
    // not. Another source can answer this endpoint and must not replace the pin.
    let controls = side_controls();
    let (a, _) = bodies[0]
        .advance_heat(
            -30_000.0,
            1.0,
            controls.caloric,
            controls.adaptivity,
            || false,
        )
        .unwrap();
    let (b, _) = bodies[1]
        .advance_heat(30_000.0, 1.0, controls.caloric, controls.adaptivity, || {
            false
        })
        .unwrap();
    let candidate_query = original
        .query()
        .with("T_a", a.state().temperature_k())
        .unwrap()
        .with("T_b", b.state().temperature_k())
        .unwrap();
    let alternative = resolve_interface_state_point(
        &interface,
        &candidate_query,
        &[ScalarPropertyRequirement::try_new(
            AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
            AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_ne!(
        alternative.properties()[0].answer().receipt.selected,
        original.properties()[0].answer().receipt.selected
    );
    let expected = CrystallizationInterfaceError::from(
        original
            .at_query_point(&interface, &candidate_query)
            .unwrap_err(),
    );
    assert_eq!(contact.evaluate(&a, &b).unwrap_err(), expected);
    let refused = contact
        .exchange_heat(
            &bodies[0],
            &bodies[1],
            1.0,
            controls,
            controls,
            exchange_controls(30_000.0),
            || false,
        )
        .err()
        .unwrap();
    assert_eq!(refused, expected);
    assert_eq!(bodies.each_ref().map(|body| body.state()), original_states);
    assert_eq!(contact.original_source(), &original);

    let perfect = source(&cards, true);
    assert!(matches!(
        prepare(
            ThermalInterfacePoint {
                interface: &perfect,
                ..point
            },
            pairs.each_ref(),
            bodies.each_ref()
        ),
        Err(CrystallizationInterfaceError::Invalid(
            "zero contact resistance requires unsupported infinite conductance"
        ))
    ));
    // The pressure coordinate keeps its own meaning, independent of body p.
    assert_eq!(context.axis_quantities()["pressure"].dims(), Pressure::DIMS);
}
