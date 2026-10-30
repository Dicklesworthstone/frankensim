//! Original thermal, speed and path coordinates through the associated ellipse.
use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::{
    EllipticFrictionMode, SourceThermalEllipticCoulomb,
};
use fs_contact::interface_binding::thermal_friction::{FrictionQueryAxis, FrictionQueryAxisKind};

const STATIC_Y: &str = "owned-query-static-y";
const KINETIC_Y: &str = "owned-query-kinetic-y";
const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

#[test]
fn full_matrix_motion_keeps_owned_query_coordinates_and_pins_atomic() {
    use fs_tribo::tangential::elliptic::SymmetricTangentialCompliance;
    let card = source_card();
    let law = law(&card);
    let initial = law.initial_state();
    let query = endpoint();
    let matrix = [[0.01, -0.003], [-0.003, 0.02]];
    let compliance = SymmetricTangentialCompliance::new(matrix).unwrap();
    // At the declared speed=1 and path=0.5 the original claims give
    // capacities [3,2]. Construct half a metre of associated physical slip
    // over half a second, then add the full off-diagonal mechanical recoil.
    let multiplier = 0.5 / fs_math::det::hypot(0.2, 0.4);
    let force = [1.8, 1.6];
    let slip = [0.2 * multiplier, 0.4 * multiplier];
    let mut drive = mechanical();
    drive.dt_s = 0.5;
    drive.slip_increment_m = [
        slip[0] + matrix[0][0] * force[0] + matrix[0][1] * force[1],
        slip[1] + matrix[1][0] * force[0] + matrix[1][1] * force[1],
    ];
    let step = law
        .advance_at_query_with_compliance(
            &initial,
            &drive,
            &query,
            EllipticFrictionMode::Kinetic,
            compliance,
        )
        .unwrap();
    for axis in 0..2 {
        close(step.traction_n[axis], force[axis]);
        close(step.slip_increment_m[axis], slip[axis]);
        let source = &step.state.friction_sources()[axis];
        assert_eq!(source.query(), query);
        for property in source.properties() {
            let original = initial.friction_sources()[axis]
                .property(property.requirement().name())
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                original.answer().receipt.selected
            );
            card.claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
    close(
        fs_math::det::hypot(step.slip_increment_m[0], step.slip_increment_m[1]),
        0.5,
    );
    close(step.dissipated_j, multiplier);
    let accepted = step.state.clone();
    for refused in [
        query.clone().with("path", 0.49).unwrap(),
        query.with("ambient", 296.0).unwrap(),
    ] {
        assert!(
            law.advance_at_query_with_compliance(
                &accepted,
                &drive,
                &refused,
                EllipticFrictionMode::Kinetic,
                compliance,
            )
            .is_err()
        );
        assert_eq!(step.state, accepted);
    }
}

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point() -> QueryPoint {
    let mut p = QueryPoint::new();
    for (name, value) in [
        ("T", 300.0),
        ("Ta", 290.0),
        ("Tb", 310.0),
        ("ambient", 295.0),
    ] {
        p = p
            .with_quantity(name, temperature_quantity(), value)
            .unwrap();
    }
    p.with_quantity("P", pressure_quantity(), 100_000.0)
        .unwrap()
        .with_quantity("speed", QuantitySpec::dimensional(SPEED), 0.5)
        .unwrap()
        .with_quantity("path", QuantitySpec::dimensional(LENGTH), 0.0)
        .unwrap()
}

fn owned() -> Vec<FrictionQueryAxis> {
    use FrictionQueryAxisKind::*;
    [
        ("Ta", Temperature),
        ("Tb", Temperature),
        ("speed", SlidingSpeed),
        ("path", AccumulatedDistance),
    ]
    .into_iter()
    .map(|(name, kind)| FrictionQueryAxis {
        name: name.into(),
        kind,
    })
    .collect()
}

fn mapped() -> ThermalFrictionAxes {
    ThermalFrictionAxes {
        temperature: "T".into(),
        nominal_pressure: Some("P".into()),
    }
}

fn response(name: &str, axis: &str, dims: Dims, endpoints: [(f64, f64); 2]) -> PropertyClaim {
    let mut c = curve(name, [1.0, 1.0]);
    c.value = PropertyValue::Curve {
        abscissa: axis.into(),
        abscissa_dims: dims,
        knots: endpoints.to_vec(),
        dims: Dims::NONE,
    };
    let quantity = if dims == Temperature::DIMS {
        temperature_quantity()
    } else {
        QuantitySpec::dimensional(dims)
    };
    c.validity = ValidityDomain::unconstrained().with_quantity(
        axis,
        quantity,
        endpoints[0].0,
        endpoints[1].0,
    );
    c
}

fn source_card() -> InterfaceSystemCard {
    card(vec![
        response(
            STATIC,
            "Ta",
            Temperature::DIMS,
            [(280.0, 0.9), (360.0, 0.7)],
        ),
        response(KINETIC, "speed", SPEED, [(0.0, 0.4), (2.0, 0.2)]),
        response(
            STATIC_Y,
            "Tb",
            Temperature::DIMS,
            [(280.0, 0.8), (360.0, 0.6)],
        ),
        response(KINETIC_Y, "path", LENGTH, [(0.0, 0.1), (1.0, 0.3)]),
    ])
}

fn principal(
    card: &InterfaceSystemCard,
    names: [&str; 2],
    axes: Vec<FrictionQueryAxis>,
) -> Result<SourceThermalCoulomb, SourceThermalCoulombError> {
    SourceThermalCoulomb::resolve_with_query_axes(
        bind_dry_interface_system_card(card, InputAuthority::SyntheticFixture).unwrap(),
        point(),
        requirement(names[0]),
        requirement(names[1]),
        MaterialPropertySelection::SingleClaimOnly,
        mapped(),
        axes,
    )
}

fn law(card: &InterfaceSystemCard) -> SourceThermalEllipticCoulomb {
    SourceThermalEllipticCoulomb::from_principal_sources(
        [
            principal(card, [STATIC, KINETIC], owned()).unwrap(),
            principal(card, [STATIC_Y, KINETIC_Y], owned()).unwrap(),
        ],
        [1.0, 0.0],
    )
    .unwrap()
}

fn mechanical() -> TangentialStepInput {
    TangentialStepInput {
        nominal_pressure_pa: Some(200_000.0),
        ..input(320.0, [0.03, 0.04])
    }
}

fn endpoint() -> QueryPoint {
    point()
        .with("T", 320.0)
        .unwrap()
        .with("Ta", 320.0)
        .unwrap()
        .with("Tb", 340.0)
        .unwrap()
        .with("P", 200_000.0)
        .unwrap()
        .with("speed", 1.0)
        .unwrap()
        .with("path", 0.5)
        .unwrap()
}

#[test]
fn actual_surface_temperatures_speed_and_path_update_pinned_principal_sources() {
    let card = source_card();
    let law = law(&card);
    let initial = law.initial_state();
    assert!(law.advance_sourced(&initial, &mechanical()).is_err());
    let point = endpoint();
    let step = law
        .advance_at_query(
            &initial,
            &mechanical(),
            &point,
            EllipticFrictionMode::StaticOrSliding,
        )
        .unwrap();
    let work = 0.09_f64.hypot(0.08);
    close(step.traction_n[0], 9.0 * 0.03 / work);
    close(step.traction_n[1], 4.0 * 0.04 / work);
    close(step.dissipated_j, work);
    close(step.state.constitutive_law().static_coefficients()[0], 0.8);
    close(step.state.constitutive_law().static_coefficients()[1], 0.65);
    for (i, source) in step.state.friction_sources().iter().enumerate() {
        assert_eq!(source.query(), point);
        for property in source.properties() {
            assert_eq!(
                property.answer().receipt.selected,
                initial.friction_sources()[i]
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            card.claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
    let continued = law
        .advance_at_query(
            &step.state,
            &mechanical(),
            &point,
            EllipticFrictionMode::Kinetic,
        )
        .unwrap();
    close(continued.state.constitutive_state().slip_m()[0], 0.06);
    close(continued.state.constitutive_state().slip_m()[1], 0.08);
    assert_eq!(initial, law.initial_state());
    let scalar = &law.principal_sources()[0];
    let one = scalar
        .advance_at_query(&scalar.initial_state(), &mechanical(), &point)
        .unwrap();
    close(one.dissipated_j, 3.0 * 0.05);
}

#[test]
fn kinetic_branch_cannot_spend_static_breakaway_capacity_or_replace_its_history() {
    let card = card(vec![curve(STATIC, [0.8, 0.8]), curve(KINETIC, [0.3, 0.3])]);
    let scalar = resolve_at(&card, query(300.0)).unwrap();
    let law =
        SourceThermalEllipticCoulomb::from_principal_sources([scalar.clone(), scalar], [1.0, 0.0])
            .unwrap();
    let initial = law.initial_state();
    let mut drive = input(300.0, [0.006, 0.0]);
    drive.compliance_m_per_n = 0.001;
    let static_step = law
        .advance_at_query(
            &initial,
            &drive,
            &query(300.0),
            EllipticFrictionMode::StaticOrSliding,
        )
        .unwrap();
    assert_eq!(static_step.regime, FrictionRegime::Sticking);
    close(static_step.traction_n[0], 6.0);
    let moving = law
        .advance_at_query(
            &initial,
            &drive,
            &query(300.0),
            EllipticFrictionMode::Kinetic,
        )
        .unwrap();
    assert_eq!(moving.regime, FrictionRegime::Sliding);
    close(moving.traction_n[0], 3.0);
    close(moving.slip_increment_m[0], 0.003);
    close(moving.dissipated_j, 0.009);
    assert_eq!(
        moving.state.constitutive_law().static_coefficients(),
        [0.8; 2]
    );
    let rest = law
        .advance_at_query(
            &moving.state,
            &drive,
            &query(300.0),
            EllipticFrictionMode::StaticOrSliding,
        )
        .unwrap();
    assert_eq!(rest.regime, FrictionRegime::Sticking);
    assert_eq!(
        rest.state.constitutive_state().slip_m(),
        moving.state.constitutive_state().slip_m()
    );
    assert_eq!(initial, law.initial_state());
}

#[test]
fn full_query_refuses_unowned_changes_mismatched_mechanics_and_partial_refresh() {
    let card = source_card();
    let law = law(&card);
    let state = law.initial_state();
    let saved = state.clone();
    let good = endpoint();
    let drive = mechanical();
    for invalid in [
        good.clone().with("ambient", 310.0).unwrap(),
        good.clone().with("new-axis", 1.0).unwrap(),
        good.clone()
            .with_quantity("speed", QuantitySpec::dimensional(LENGTH), 1.0)
            .unwrap(),
        good.clone().with("speed", -1.0).unwrap(),
        good.clone().with("path", -0.1).unwrap(),
        good.clone().with("T", 321.0).unwrap(),
        good.clone().with("P", 200_001.0).unwrap(),
    ] {
        assert!(
            law.advance_at_query(&state, &drive, &invalid, EllipticFrictionMode::Kinetic)
                .is_err()
        );
        assert_eq!(state, saved);
    }
    let outside = good.clone().with("Tb", 400.0).unwrap();
    assert!(matches!(
        law.advance_at_query(&state, &drive, &outside, EllipticFrictionMode::Kinetic),
        Err(SourceThermalCoulombError::Material(_))
    ));
    let mut limited = drive;
    limited.solve.maximum_force_n = 0.1;
    assert!(matches!(
        law.advance_at_query(&state, &limited, &good, EllipticFrictionMode::Kinetic),
        Err(SourceThermalCoulombError::Tribo(_))
    ));
    let retry = law
        .advance_at_query(&state, &drive, &good, EllipticFrictionMode::Kinetic)
        .unwrap();
    let fresh = law
        .advance_at_query(
            &law.initial_state(),
            &drive,
            &good,
            EllipticFrictionMode::Kinetic,
        )
        .unwrap();
    assert_eq!(retry.state, fresh.state);
    assert_eq!(retry.traction_n, fresh.traction_n);
    let retreated = good.clone().with("path", 0.25).unwrap();
    assert!(
        law.advance_at_query(
            &retry.state,
            &drive,
            &retreated,
            EllipticFrictionMode::Kinetic
        )
        .is_err()
    );
    assert_eq!(state, saved);
}

#[test]
fn coordinate_ownership_is_typed_canonical_and_part_of_source_identity() {
    let card = source_card();
    let original = principal(&card, [STATIC, KINETIC], owned()).unwrap();
    let mut reordered = owned();
    reordered.reverse();
    assert_eq!(
        original,
        principal(&card, [STATIC, KINETIC], reordered).unwrap()
    );
    let mut duplicated = owned();
    duplicated.push(duplicated[0].clone());
    assert!(principal(&card, [STATIC, KINETIC], duplicated).is_err());
    let mut wrong = owned();
    wrong[2].kind = FrictionQueryAxisKind::Temperature;
    assert!(principal(&card, [STATIC, KINETIC], wrong).is_err());
    assert!(principal(&card, [STATIC, KINETIC], Vec::new()).is_err());
    let mut broader = owned();
    broader.push(FrictionQueryAxis {
        name: "ambient".into(),
        kind: FrictionQueryAxisKind::Temperature,
    });
    let other = principal(&card, [STATIC_Y, KINETIC_Y], broader).unwrap();
    assert!(
        SourceThermalEllipticCoulomb::from_principal_sources([original, other], [1.0, 0.0])
            .is_err()
    );
}
