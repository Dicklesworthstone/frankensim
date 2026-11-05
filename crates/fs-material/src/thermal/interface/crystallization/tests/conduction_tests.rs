//! G1/G3: original bulk conductivity and contact resistance share one heat solve.
//! Manufactured laws permit an independent closed form; no material fit is claimed.
use super::*;
use crate::phase::crystallization::properties::CrystallizationMaterialPoint;
use crate::state_point::{THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY};
use fs_qty::Dims;

const DISTANCES: [f64; 2] = [0.001, 0.002];

fn conductivity(axis: &str, range: [f64; 2], values: [f64; 2]) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS),
        value: PropertyValue::Curve {
            abscissa: axis.into(),
            abscissa_dims: Dims::NONE,
            knots: vec![(range[0], values[0]), (range[1], values[1])],
            dims: THERMAL_CONDUCTIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "T",
                quantity(QuantityKind::AbsoluteTemperature),
                280.0,
                550.0,
            )
            .with_quantity(
                axis,
                QuantitySpec::dimensional(Dims::NONE),
                range[0],
                range[1],
            )
            .with_quantity(
                "pressure",
                quantity(QuantityKind::Pressure),
                90_000.0,
                110_000.0,
            ),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "authored conductivity versus absolute kinetic crystal fraction".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn bulk_cards(insulating: bool, restricted: bool) -> [MaterialCard; 2] {
    std::array::from_fn(|side| {
        let mut claims = ClaimSet::new();
        if side == 0 {
            let upper = if restricted { 0.05 } else { 0.8 };
            claims
                .insert_claim(conductivity(
                    "X",
                    [0.0, upper],
                    if insulating {
                        [0.0; 2]
                    } else {
                        [2.0, 2.0 + 2.0 * upper]
                    },
                ))
                .unwrap();
            if restricted {
                claims
                    .insert_claim(conductivity("X", [0.06, 0.8], [10.0; 2]))
                    .unwrap();
            }
        } else {
            claims
                .insert_claim(conductivity("V", [0.0, 0.8], [1.0, 1.8]))
                .unwrap();
        }
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: format!("authored conductive side {side}"),
                phase: "original parent/crystal pair".into(),
                process: "kinetic bulk contact regression".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap()
    })
}

fn bulk_contexts() -> [QueryPoint; 2] {
    ["X", "V"].map(|axis| {
        QueryPoint::new()
            .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 300.0)
            .unwrap()
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), 0.0)
            .unwrap()
            .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
            .unwrap()
    })
}

fn constant_growth(t: f64) -> Result<f64, CrystallizationError> {
    // Constant on every tested trajectory, continuous and zero at/above melting.
    Ok(if t <= 450.0 {
        0.1
    } else {
        (0.002 * (500.0 - t)).max(0.0)
    })
}

fn growing_laws<'a>(pairs: &'a [CrystallizationCaloricPair<'_>; 2]) -> [CrystallizationLaw<'a>; 2] {
    std::array::from_fn(|side| {
        pairs[side]
            .law(
                &constant_growth,
                ContentHash([0xd1 + side as u8; 32]),
                1.0,
                0.8,
                500.0,
            )
            .unwrap()
    })
}

fn growing_parcels<'a, 's>(
    laws: &'a [CrystallizationLaw<'s>; 2],
) -> [IsobaricCrystallizationParcel<'a, 's>; 2] {
    std::array::from_fn(|side| {
        IsobaricCrystallizationParcel::new(
            &laws[side],
            laws[side]
                .state_at_temperature([400.0, 330.0][side], 0.0)
                .unwrap(),
            1.0,
            [0.001, 0.00125][side],
            100_000.0,
            side_controls().caloric,
            || false,
        )
        .unwrap()
    })
}

fn prepare_bulk<'a>(
    cards: &'a [MaterialCard; 2],
    contexts: &'a [QueryPoint; 2],
    pairs: [CrystallizationCaloricPair<'a>; 2],
    bodies: [&IsobaricCrystallizationParcel<'_, 'a>; 2],
    distances: [f64; 2],
) -> Result<CrystallizationBulkConduction<'a>, CrystallizationInterfaceError> {
    CrystallizationBulkConduction::new(
        std::array::from_fn(|side| CrystallizationMaterialPoint {
            material: &cards[side],
            context: &contexts[side],
            temperature_axis: "T",
            crystalline_mass_fraction_axis: (side == 0).then_some("X"),
            crystalline_volume_fraction_axis: (side == 1).then_some("V"),
        }),
        pairs,
        bodies,
        distances,
        [
            MaterialPropertySelection::SingleClaimOnly,
            MaterialPropertySelection::SingleClaimOnly,
        ],
        [side_controls().caloric.caloric; 2],
    )
}

#[test]
fn g1_growth_dependent_bulk_conduction_closes_finite_contact_and_latent_feedback() {
    let cards = bulk_cards(false, false);
    let curves = charts(&cards);
    let pairs = curves
        .each_ref()
        .map(|c| CrystallizationCaloricPair::new(&c[0], &c[1], [280.0, 550.0]).unwrap());
    let laws = growing_laws(&pairs);
    let bodies = growing_parcels(&laws);
    let bulk_context = bulk_contexts();
    let bulk = prepare_bulk(&cards, &bulk_context, pairs, bodies.each_ref(), DISTANCES).unwrap();
    let selected = bulk
        .original_sources()
        .map(|s| s.properties()[0].answer().receipt.selected);
    let interface = source(&cards, true);
    let context = context();
    // Perfect interface contact alone is still singular. The original finite
    // bulk resistance, rather than an artificial contact epsilon, closes G.
    assert!(
        prepare(
            point(&interface, &cards, &context),
            pairs.each_ref(),
            bodies.each_ref()
        )
        .is_err()
    );
    let contact = CrystallizationInterfaceContact::new_with_bulk_conduction(
        point(&interface, &cards, &context),
        pairs.each_ref(),
        bodies.each_ref(),
        0.01,
        Some("pressure"),
        Some("gap"),
        MaterialPropertySelection::SingleClaimOnly,
        bulk,
    )
    .unwrap();
    let run = || {
        contact
            .exchange_heat(
                &bodies[0],
                &bodies[1],
                2.0,
                side_controls(),
                side_controls(),
                exchange_controls(2000.0),
                || false,
            )
            .unwrap()
    };
    let (left, right, report) = run();
    let x = 0.8 * (1.0 - (-0.2_f64).exp());
    let v_right = (x / 1000.0) / ((1.0 - x) / 800.0 + x / 1000.0);
    let k = [2.0 + 2.0 * x, 1.0 + v_right];
    let g = 0.01 / (DISTANCES[0] / k[0] + DISTANCES[1] / k[1]);
    // Original equal-in-phase Cp values give latent recalescence 50X and 80X.
    // With constant K, X is independent of Q: solve backward-Euler heat exactly.
    let q = 2.0 * g * (70.0 - 30.0 * x) / (1.0 + 2.0 * g * (1.0 / 2000.0 + 1.0 / 1000.0));
    close(report.thermal.heat_left_to_right_j, q, 3e-6);
    close(report.interface.conductance_w_k, g, 1e-10);
    close(
        left.state().temperature_k(),
        400.0 + 50.0 * x - q / 2000.0,
        1e-8,
    );
    close(
        right.state().temperature_k(),
        330.0 + 80.0 * x + q / 1000.0,
        1e-8,
    );
    let accepted = report.interface.bulk.as_ref().unwrap();
    for (side, parcel) in [&left, &right].into_iter().enumerate() {
        close(parcel.state().crystalline_mass_fraction(), x, 1e-12);
        assert_eq!(accepted.sources[side].phase(), parcel.state());
        close(accepted.conductivity_w_m_k[side], k[side], 1e-12);
        assert_eq!(
            accepted.sources[side].material().properties()[0]
                .answer()
                .receipt
                .selected,
            selected[side]
        );
        assert_eq!(
            accepted.sources[side]
                .material()
                .state_coordinate("pressure"),
            Some(100_000.0)
        );
    }
    close(
        accepted.sources[0]
            .material()
            .state_coordinate("X")
            .unwrap(),
        x,
        1e-12,
    );
    close(
        accepted.sources[1]
            .material()
            .state_coordinate("V")
            .unwrap(),
        v_right,
        1e-12,
    );
    close(report.thermal.combined_enthalpy_change_j, 0.0, 1e-6);
    close(
        report.thermal.combined_internal_energy_change_j + report.thermal.combined_boundary_work_j,
        0.0,
        1e-6,
    );
    assert!(report.thermal.absolute_energy_residual_j <= 1e-5);
    let frozen_g = 0.01 / (DISTANCES[0] / 2.0 + DISTANCES[1]);
    let frozen_q = 2.0 * frozen_g * (70.0 - 30.0 * x) / (1.0 + 2.0 * frozen_g * 0.0015);
    assert!((q - frozen_q).abs() > 10.0);
    assert_eq!(run().2, report);
    assert_eq!(bodies[0].state().temperature_k(), 400.0);
    assert_eq!(bodies[1].state().temperature_k(), 330.0);
}

#[test]
fn g0_bulk_and_original_interface_resistance_are_both_retained() {
    let cards = bulk_cards(false, false);
    let curves = charts(&cards);
    let pairs = curves
        .each_ref()
        .map(|c| CrystallizationCaloricPair::new(&c[0], &c[1], [280.0, 550.0]).unwrap());
    let laws = growing_laws(&pairs);
    let bodies = growing_parcels(&laws);
    let contexts = bulk_contexts();
    assert!(prepare_bulk(&cards, &contexts, pairs, bodies.each_ref(), [0.0, 0.002]).is_err());
    let wrong_pressure = [
        contexts[0].clone().with("pressure", 100_001.0).unwrap(),
        contexts[1].clone(),
    ];
    assert!(prepare_bulk(&cards, &wrong_pressure, pairs, bodies.each_ref(), DISTANCES).is_err());
    let bulk = prepare_bulk(&cards, &contexts, pairs, bodies.each_ref(), DISTANCES).unwrap();
    let interface = source(&cards, false);
    let context = context();
    let contact = CrystallizationInterfaceContact::new_with_bulk_conduction(
        point(&interface, &cards, &context),
        pairs.each_ref(),
        bodies.each_ref(),
        0.01,
        Some("pressure"),
        Some("gap"),
        MaterialPropertySelection::SingleClaimOnly,
        bulk,
    )
    .unwrap();
    let response = contact.evaluate(&bodies[0], &bodies[1]).unwrap();
    let r = original_resistance(400.0, 330.0);
    close(response.resistance_m2_k_w, r, 1e-15);
    close(
        response.conductance_w_k,
        0.01 / (0.001 / 2.0 + 0.002 + r),
        1e-12,
    );
    assert!(response.conductance_w_k < 0.01 / r);
    assert!(response.conductance_w_k < 0.01 / (0.001 / 2.0 + 0.002));
    assert_eq!(
        response.source.properties()[0].answer().receipt.selected,
        contact.original_source().properties()[0]
            .answer()
            .receipt
            .selected
    );
    assert!(matches!(
        contact.exchange_heat(
            &bodies[0],
            &bodies[1],
            1.0,
            side_controls(),
            side_controls(),
            exchange_controls(2000.0),
            || true
        ),
        Err(CrystallizationInterfaceError::Thermal(
            CrystallizationError::Cancelled
        ))
    ));
}

#[test]
fn g1_supported_zero_conductivity_insulates_without_stopping_phase_history() {
    let cards = bulk_cards(true, false);
    let curves = charts(&cards);
    let pairs = curves
        .each_ref()
        .map(|c| CrystallizationCaloricPair::new(&c[0], &c[1], [280.0, 550.0]).unwrap());
    let laws = growing_laws(&pairs);
    let bodies = growing_parcels(&laws);
    let contexts = bulk_contexts();
    let bulk = prepare_bulk(&cards, &contexts, pairs, bodies.each_ref(), DISTANCES).unwrap();
    let interface = source(&cards, true);
    let context = context();
    let contact = CrystallizationInterfaceContact::new_with_bulk_conduction(
        point(&interface, &cards, &context),
        pairs.each_ref(),
        bodies.each_ref(),
        0.01,
        Some("pressure"),
        Some("gap"),
        MaterialPropertySelection::SingleClaimOnly,
        bulk,
    )
    .unwrap();
    let (left, right, report) = contact
        .exchange_heat(
            &bodies[0],
            &bodies[1],
            0.5,
            side_controls(),
            side_controls(),
            IsobaricExchangeControls {
                heat_bracket_j: [0.0, 0.0],
                ..exchange_controls(0.0)
            },
            || false,
        )
        .unwrap();
    assert_eq!(report.interface.conductance_w_k, 0.0);
    assert_eq!(report.thermal.heat_left_to_right_j, 0.0);
    for (side, parcel) in [&left, &right].into_iter().enumerate() {
        assert!(parcel.state().crystalline_mass_fraction() > 0.0);
        assert!(parcel.state().temperature_k() > bodies[side].state().temperature_k());
        assert_eq!(
            parcel.state().specific_enthalpy_j_kg(),
            bodies[side].state().specific_enthalpy_j_kg()
        );
    }
}

#[test]
fn g4_a_later_bulk_source_failure_cannot_switch_claims_or_publish_either_body() {
    let cards = bulk_cards(false, true);
    let curves = charts(&cards);
    let pairs = curves
        .each_ref()
        .map(|c| CrystallizationCaloricPair::new(&c[0], &c[1], [280.0, 550.0]).unwrap());
    let laws = growing_laws(&pairs);
    let bodies = growing_parcels(&laws);
    let contexts = bulk_contexts();
    let bulk = prepare_bulk(&cards, &contexts, pairs, bodies.each_ref(), DISTANCES).unwrap();
    let selected = bulk.original_sources()[0].properties()[0]
        .answer()
        .receipt
        .selected;
    let interface = source(&cards, true);
    let context = context();
    let contact = CrystallizationInterfaceContact::new_with_bulk_conduction(
        point(&interface, &cards, &context),
        pairs.each_ref(),
        bodies.each_ref(),
        0.01,
        Some("pressure"),
        Some("gap"),
        MaterialPropertySelection::SingleClaimOnly,
        bulk,
    )
    .unwrap();
    let initial = bodies.each_ref().map(|p| p.state());
    let result = contact.exchange_heat(
        &bodies[0],
        &bodies[1],
        2.0,
        side_controls(),
        side_controls(),
        exchange_controls(2000.0),
        || false,
    );
    assert!(matches!(
        result,
        Err(CrystallizationInterfaceError::Property(_))
    ));
    assert_eq!(bodies.each_ref().map(|p| p.state()), initial);
    assert_eq!(
        contact.bulk_conduction().unwrap().original_sources()[0].properties()[0]
            .answer()
            .receipt
            .selected,
        selected
    );
    assert_eq!(
        contact.bulk_conduction().unwrap().original_sources()[0].state_coordinate("X"),
        Some(0.0)
    );
}
