//! Manufactured contact laws test live phase feedback independently of source fitting.
use super::*;
use fs_qty::Dims;

fn axes() -> CrystallizationInterfacePhaseAxes<'static> {
    CrystallizationInterfacePhaseAxes {
        crystalline_mass_fraction: [Some("X_a"), Some("X_b")],
        crystalline_volume_fraction: [Some("V_a"), Some("V_b")],
    }
}

fn phase_context() -> QueryPoint {
    ["X_a", "X_b", "V_a", "V_b"]
        .into_iter()
        .fold(context(), |point, axis| {
            // Each placeholder must be replaced by its actual parcel coordinate.
            point
                .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), 0.7)
                .unwrap()
        })
}

fn phase_source(cards: &[MaterialCard; 2], volume: bool, restricted: bool) -> InterfaceSystemCard {
    let right_axis = if volume { "V_b" } else { "X_b" };
    let mut claims = ClaimSet::new();
    for (range, original) in if restricted {
        vec![([0.45, 0.8], true), ([0.0, 0.449], false)]
    } else {
        vec![([0.0, 0.8], true)]
    } {
        let mut validity = ValidityDomain::unconstrained()
            .with_quantity(
                "T_a",
                quantity(QuantityKind::AbsoluteTemperature),
                280.0,
                550.0,
            )
            .with_quantity(
                "T_b",
                quantity(QuantityKind::AbsoluteTemperature),
                280.0,
                550.0,
            )
            .with_quantity(
                "pressure",
                quantity(QuantityKind::Pressure),
                90_000.0,
                110_000.0,
            )
            .with_quantity("gap", QuantitySpec::dimensional(Length::DIMS), 0.0, 1e-4);
        for axis in ["X_a", "X_b", "V_a", "V_b"] {
            let bounds = if axis == right_axis {
                range
            } else {
                [0.0, 0.8]
            };
            validity = validity.with_quantity(
                axis,
                QuantitySpec::dimensional(Dims::NONE),
                bounds[0],
                bounds[1],
            );
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(
                    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
                ),
                value: PropertyValue::Surface2D {
                    axes: ["X_a".into(), right_axis.into()],
                    axis_dims: [Dims::NONE; 2],
                    samples: [
                        [0.0, range[0]],
                        [0.8, range[0]],
                        [0.0, range[1]],
                        [0.8, range[1]],
                    ]
                    .map(|x| {
                        (
                            x,
                            if original {
                                1e-4 * (1.0 + x[0] + 2.0 * x[1])
                            } else {
                                4e-4
                            },
                        )
                    })
                    .to_vec(),
                    triangles: vec![[0, 1, 2], [1, 3, 2]],
                    dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
                },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored contact resistance versus original crystal fractions".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "side A".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "side B".into(),
        },
        SystemContext {
            medium: "fixed contact layer".into(),
            third_body: None,
            environment: "fixed geometry and pressure".into(),
            history: "uniform bulk kinetic inventories at contact".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn prepare_phase<'a>(
    point: ThermalInterfacePoint<'a>,
    pairs: [&CrystallizationCaloricPair<'a>; 2],
    bodies: [&IsobaricCrystallizationParcel<'_, '_>; 2],
) -> Result<CrystallizationInterfaceContact<'a>, CrystallizationInterfaceError> {
    CrystallizationInterfaceContact::new_with_crystallinity(
        point,
        pairs,
        bodies,
        0.01,
        Some("pressure"),
        Some("gap"),
        axes(),
        MaterialPropertySelection::SingleClaimOnly,
        None,
    )
}

#[test]
fn g1_partial_remelting_changes_contact_at_constant_temperature_and_retains_actual_fraction_receipts()
 {
    let cards = materials();
    let curves = charts(&cards);
    let pairs = curves
        .each_ref()
        .map(|c| CrystallizationCaloricPair::new(&c[0], &c[1], [280.0, 550.0]).unwrap());
    let laws = [
        pairs[0]
            .law(&no_growth, ContentHash([0xe1; 32]), 1.0, 0.8, 500.0)
            .unwrap(),
        pairs[1]
            .law(&growth, ContentHash([0xe2; 32]), 1.0, 0.8, 350.0)
            .unwrap(),
    ];
    let bodies = parcels(&laws);
    let interface = phase_source(&cards, false, false);
    let context = phase_context();
    let contact = prepare_phase(
        point(&interface, &cards, &context),
        pairs.each_ref(),
        bodies.each_ref(),
    )
    .unwrap();
    let original = contact.original_source().clone();
    let initial_x = bodies[1].state().crystalline_mass_fraction();
    assert_eq!(original.state_coordinate("X_a"), Some(0.0));
    assert_eq!(original.state_coordinate("X_b"), Some(initial_x));
    let (left, right, report) = contact
        .exchange_heat(
            &bodies[0],
            &bodies[1],
            1.0,
            side_controls(),
            side_controls(),
            exchange_controls(5000.0),
            || false,
        )
        .unwrap();
    // m_b*L=40000 J; R=1e-4*(1+2*(X0-Q/40000)); T_a=400-Q/2000.
    // The lower quadratic root solves Q=A/R*(T_a-350), independently of the kernel.
    let b = 1.05 + 2.0 * initial_x;
    let heat = 10_000.0 / (b + (b * b - 1.0).sqrt());
    close(report.thermal.heat_left_to_right_j, heat, 3e-6);
    close(right.state().temperature_k(), 350.0, 1e-9);
    close(
        right.state().crystalline_mass_fraction(),
        initial_x - heat / 40_000.0,
        1e-10,
    );
    close(
        report.interface.resistance_m2_k_w,
        1e-4 * (1.0 + 2.0 * right.state().crystalline_mass_fraction()),
        1e-18,
    );
    assert_eq!(
        report.interface.source.state_coordinate("X_b"),
        Some(right.state().crystalline_mass_fraction())
    );
    assert_eq!(
        report.interface.source.state_coordinate("T_a"),
        Some(left.state().temperature_k())
    );
    let x = right.state().crystalline_mass_fraction();
    close(
        report.interface.source.state_coordinate("V_b").unwrap(),
        0.8 * x / (1.0 - 0.2 * x),
        2e-15,
    );
    assert_eq!(
        report.interface.source.properties()[0]
            .answer()
            .receipt
            .selected,
        original.properties()[0].answer().receipt.selected
    );
    assert_eq!(contact.original_source(), &original);
    close(
        report.thermal.combined_internal_energy_change_j + report.thermal.combined_boundary_work_j,
        0.0,
        1e-6,
    );
    let frozen_g = 100.0 / (1.0 + 2.0 * initial_x);
    let frozen_heat = frozen_g * 50.0 / (1.0 + frozen_g / 2000.0);
    assert!(
        heat - frozen_heat > 100.0,
        "stale input crystallinity must change the physical answer"
    );
}

fn constant_growth(t: f64) -> Result<f64, CrystallizationError> {
    Ok(if t <= 450.0 {
        0.1
    } else {
        (0.002 * (500.0 - t)).max(0.0)
    })
}

#[test]
fn g1_live_mass_and_volume_coordinates_close_growth_latent_heat_and_contact_together() {
    let cards = materials();
    let curves = charts(&cards);
    let pairs = curves
        .each_ref()
        .map(|c| CrystallizationCaloricPair::new(&c[0], &c[1], [280.0, 550.0]).unwrap());
    let laws = pairs.each_ref().map(|pair| {
        pair.law(&constant_growth, ContentHash([0xe3; 32]), 1.0, 0.8, 500.0)
            .unwrap()
    });
    let bodies = std::array::from_fn::<_, 2, _>(|side| {
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
    });
    let interface = phase_source(&cards, true, false);
    let context = phase_context();
    let contact = prepare_phase(
        point(&interface, &cards, &context),
        pairs.each_ref(),
        bodies.each_ref(),
    )
    .unwrap();
    let (left, right, report) = contact
        .exchange_heat(
            &bodies[0],
            &bodies[1],
            1.0,
            side_controls(),
            side_controls(),
            exchange_controls(7000.0),
            || false,
        )
        .unwrap();
    let x = 0.8 * -(-0.1_f64).exp_m1();
    let volume = 0.8 * x / (1.0 - 0.2 * x);
    let g = 100.0 / (1.0 + x + 2.0 * volume);
    let heat = g * (70.0 - 30.0 * x) / (1.0 + 0.0015 * g);
    close(report.thermal.heat_left_to_right_j, heat, 3e-6);
    close(left.state().crystalline_mass_fraction(), x, 1e-14);
    close(right.state().crystalline_mass_fraction(), x, 1e-14);
    close(
        left.state().temperature_k(),
        400.0 + 50.0 * x - heat / 2000.0,
        2e-9,
    );
    close(
        right.state().temperature_k(),
        330.0 + 80.0 * x + heat / 1000.0,
        2e-9,
    );
    close(report.interface.conductance_w_k, g, 1e-11);
    close(
        report.interface.source.state_coordinate("V_b").unwrap(),
        volume,
        2e-15,
    );
    assert!(
        (volume - x).abs() > 0.01,
        "unequal phase densities separate volume and mass fractions"
    );
    assert_eq!(
        report.interface.source.state_coordinate("X_a"),
        Some(left.state().crystalline_mass_fraction())
    );
    assert_eq!(
        report.interface.source.state_coordinate("X_b"),
        Some(right.state().crystalline_mass_fraction())
    );
    close(
        report.thermal.combined_internal_energy_change_j + report.thermal.combined_boundary_work_j,
        0.0,
        1e-6,
    );
}

#[test]
fn g4_fraction_domain_failure_is_atomic_and_never_switches_the_pinned_claim() {
    let cards = materials();
    let curves = charts(&cards);
    let pairs = curves
        .each_ref()
        .map(|c| CrystallizationCaloricPair::new(&c[0], &c[1], [280.0, 550.0]).unwrap());
    let laws = [
        pairs[0]
            .law(&no_growth, ContentHash([0xe1; 32]), 1.0, 0.8, 500.0)
            .unwrap(),
        pairs[1]
            .law(&growth, ContentHash([0xe2; 32]), 1.0, 0.8, 350.0)
            .unwrap(),
    ];
    let bodies = parcels(&laws);
    let original_states = bodies.each_ref().map(|body| body.state());
    let interface = phase_source(&cards, false, true);
    let context = phase_context();
    let point = point(&interface, &cards, &context);
    let contact = prepare_phase(point, pairs.each_ref(), bodies.each_ref()).unwrap();
    let original = contact.original_source().clone();
    let (left, _) = bodies[0]
        .apply_heat(-4000.0, side_controls().caloric, || false)
        .unwrap();
    let (right, _) = bodies[1]
        .apply_heat(4000.0, side_controls().caloric, || false)
        .unwrap();
    let query = endpoint_query(
        original.query(),
        point.temperature_axes,
        None,
        [&left, &right],
        contact.phase_coordinates.as_ref(),
    )
    .unwrap();
    let alternative = resolve_interface_state_point(
        &interface,
        &query,
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
        original.at_query_point(&interface, &query).unwrap_err(),
    );
    assert_eq!(contact.evaluate(&left, &right).unwrap_err(), expected);
    assert!(matches!(
        contact.exchange_heat(
            &bodies[0],
            &bodies[1],
            1.0,
            side_controls(),
            side_controls(),
            exchange_controls(5000.0),
            || false,
        ),
        Err(CrystallizationInterfaceError::Material(_))
    ));
    assert_eq!(bodies.each_ref().map(|body| body.state()), original_states);
    assert_eq!(bodies.each_ref().map(|body| body.time_s()), [0.0; 2]);
    assert_eq!(contact.original_source(), &original);
    for invalid_axes in [
        CrystallizationInterfacePhaseAxes {
            crystalline_mass_fraction: [Some("T_a"), Some("X_b")],
            ..axes()
        },
        CrystallizationInterfacePhaseAxes {
            crystalline_volume_fraction: [Some("X_a"), Some("V_b")],
            ..axes()
        },
    ] {
        assert!(
            CrystallizationInterfaceContact::new_with_crystallinity(
                point,
                pairs.each_ref(),
                bodies.each_ref(),
                0.01,
                Some("pressure"),
                Some("gap"),
                invalid_axes,
                MaterialPropertySelection::SingleClaimOnly,
                None,
            )
            .is_err()
        );
    }
    let wrong_units = context
        .clone()
        .with_quantity("X_b", quantity(QuantityKind::Pressure), 0.5)
        .unwrap();
    assert!(matches!(
        prepare_phase(
            ThermalInterfacePoint {
                context: &wrong_units,
                ..point
            },
            pairs.each_ref(),
            bodies.each_ref(),
        ),
        Err(CrystallizationInterfaceError::Invalid(
            "contact crystallinity requires typed dimensional scalar fraction axes"
        ))
    ));
}
