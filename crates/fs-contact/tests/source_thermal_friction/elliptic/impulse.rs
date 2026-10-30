//! G1/G3/G4: source roles, preimpact sampling and uninterrupted support history.
use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::SourceThermalEllipticState;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::*;
use fs_tribo::tangential::elliptic::impulse::{EllipticImpulseControls, EllipticImpulseInput};

fn controls() -> EllipticImpulseControls {
    EllipticImpulseControls {
        maximum_iterations: 128,
        maximum_impulse_ns: 100.0,
        impulse_absolute_tolerance_ns: 1e-13,
        impulse_relative_tolerance: 1e-13,
        velocity_tolerance_m_s: 1e-10,
        energy_tolerance_j: 1e-10,
    }
}
fn collision(velocity: [f64; 2]) -> EllipticImpulseInput {
    EllipticImpulseInput {
        normal_impulse_ns: 2.0,
        velocity_before_m_s: velocity,
        mobility_per_kg: [[1.75, 0.0], [0.0, 1.75]],
    }
}
fn boundary(
    query: QueryPoint,
    temperature: f64,
    pressure: Option<f64>,
) -> SourceFrictionConditions {
    SourceFrictionConditions {
        query,
        contact_temperature_k: Some(temperature),
        nominal_pressure_pa: pressure,
    }
}
fn unchanged_history(old: &SourceThermalEllipticState, new: &SourceThermalEllipticState) {
    assert_eq!(
        new.constitutive_state().slip_m(),
        old.constitutive_state().slip_m()
    );
    assert_eq!(
        new.constitutive_state().principal_axis(),
        old.constitutive_state().principal_axis()
    );
    assert_eq!(
        new.constitutive_state().contact_temperature_k(),
        old.constitutive_state().contact_temperature_k()
    );
}
fn original_claims(
    card: &InterfaceSystemCard,
    old: &SourceThermalEllipticState,
    new: &SourceThermalEllipticState,
) {
    for (axis, source) in new.friction_sources().iter().enumerate() {
        for property in source.properties() {
            let previous = old.friction_sources()[axis]
                .property(property.requirement().name())
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                previous.answer().receipt.selected
            );
            card.claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
}

#[test]
fn explicit_roles_keep_original_principal_claims_and_analytic_impact_loss() {
    let card = directional_card();
    let law = law(&card, [1.0, 0.0]);
    let state = law
        .advance_sourced(&law.initial_state(), &input(300.0, [0.01, 0.02]))
        .unwrap()
        .state;
    let original = state.clone();
    for (temperature, role, coefficient, impulse, after, loss) in [
        (
            280.0,
            SourceImpulseCoefficientRole::Kinetic,
            0.6,
            -4.0 / 7.0,
            0.0,
            2.0 / 7.0,
        ),
        (
            320.0,
            SourceImpulseCoefficientRole::Kinetic,
            0.2,
            -0.4,
            0.3,
            0.26,
        ),
        (
            320.0,
            SourceImpulseCoefficientRole::Static,
            0.4,
            -4.0 / 7.0,
            0.0,
            2.0 / 7.0,
        ),
    ] {
        let report = law
            .solve_preimpact_impulse(
                &state,
                &collision([1.0, 0.0]),
                &boundary(query(temperature), temperature, None),
                role,
                controls(),
                &mut || false,
            )
            .unwrap();
        close(report.coefficients_mu[0], coefficient);
        assert_eq!(report.coefficient_role, role);
        close(report.impulse.impulse_ns[0], impulse);
        close(report.impulse.impulse_ns[1], 0.0);
        close(report.impulse.velocity_after_m_s[0], after);
        close(report.impulse.dissipated_j, loss);
        close(report.impulse.kinetic_energy_change_j, -loss);
        close(report.impulse.energy_residual_j, 0.0);
        close(report.queried_speed_m_s, 1.0);
        unchanged_history(&state, &report.state);
        original_claims(&card, &state, &report.state);
        assert!(
            report
                .state
                .friction_sources()
                .iter()
                .all(|source| source.state_coordinate("T") == Some(temperature))
        );
        if temperature == 320.0 {
            for (got, want) in report
                .state
                .constitutive_law()
                .static_coefficients()
                .into_iter()
                .zip([0.4, 0.9])
            {
                close(got, want);
            }
            for (got, want) in report
                .state
                .constitutive_law()
                .kinetic_coefficients()
                .into_iter()
                .zip([0.2, 0.5])
            {
                close(got, want);
            }
        }
        assert_eq!(state, original);
    }
}

#[test]
fn source_bound_off_axis_impulse_retains_texture_full_mobility_and_physical_work() {
    let card = directional_card();
    let law = law(&card, [0.6, 0.8]);
    let state = law.initial_state();
    let world = |v: [f64; 2]| [0.6 * v[0] - 0.8 * v[1], 0.8 * v[0] + 0.6 * v[1]];
    let endpoint = world([0.7, -0.4]);
    let (reaction, _) = expected_force([0.8, 0.6], [0.7, -0.4]);
    let reaction = world(reaction);
    let matrix = [[1.8, 0.35], [0.35, 0.9]];
    let before = [
        endpoint[0] + matrix[0][0] * reaction[0] + matrix[0][1] * reaction[1],
        endpoint[1] + matrix[1][0] * reaction[0] + matrix[1][1] * reaction[1],
    ];
    let request = EllipticImpulseInput {
        mobility_per_kg: matrix,
        ..collision(before)
    };
    let report = law
        .solve_preimpact_impulse(
            &state,
            &request,
            &boundary(query(300.0), 300.0, None),
            SourceImpulseCoefficientRole::Kinetic,
            controls(),
            &mut || false,
        )
        .unwrap();
    for i in 0..2 {
        close(report.impulse.impulse_ns[i], -reaction[i]);
        close(report.impulse.velocity_after_m_s[i], endpoint[i]);
    }
    let loss =
        0.5 * (reaction[0] * (before[0] + endpoint[0]) + reaction[1] * (before[1] + endpoint[1]));
    close(report.impulse.dissipated_j, loss);
    close(report.impulse.kinetic_energy_change_j, -loss);
    assert!(report.impulse.iterations > 0);
    unchanged_history(&state, &report.state);
    original_claims(&card, &state, &report.state);
}

fn multiple_motion_sources(
    card: &InterfaceSystemCard,
) -> (SourceThermalEllipticCoulomb, QueryPoint) {
    let extend = |query: QueryPoint, speed, path| {
        query
            .with_quantity("speed_aux", QuantitySpec::dimensional(SPEED), speed)
            .unwrap()
            .with_quantity("path_aux", QuantitySpec::dimensional(LENGTH), path)
            .unwrap()
    };
    let accepted = extend(point(), 7.0, 1.25);
    let principal = [[STATIC, KINETIC], [STATIC_Y, KINETIC_Y]].map(|names| {
        let owned = [
            ("Ta", FrictionQueryAxisKind::Temperature),
            ("Tb", FrictionQueryAxisKind::Temperature),
            ("speed", FrictionQueryAxisKind::SlidingSpeed),
            ("speed_aux", FrictionQueryAxisKind::SlidingSpeed),
            ("path", FrictionQueryAxisKind::AccumulatedDistance),
            ("path_aux", FrictionQueryAxisKind::AccumulatedDistance),
        ]
        .into_iter()
        .map(|(name, kind)| FrictionQueryAxis {
            name: name.into(),
            kind,
        })
        .collect();
        SourceThermalCoulomb::resolve_with_query_axes(
            bind_dry_interface_system_card(card, InputAuthority::SyntheticFixture).unwrap(),
            accepted.clone(),
            requirement(names[0]),
            requirement(names[1]),
            MaterialPropertySelection::SingleClaimOnly,
            ThermalFrictionAxes {
                temperature: "T".into(),
                nominal_pressure: Some("p".into()),
            },
            owned,
        )
        .unwrap()
    });
    (
        SourceThermalEllipticCoulomb::from_principal_sources(principal, [1.0, 0.0]).unwrap(),
        extend(endpoint(), 99.0, 99.0),
    )
}

#[test]
fn every_speed_uses_actual_preimpact_motion_and_every_path_retains_accepted_bits() {
    let card = motion_card([0.3, 0.5], [0.2, 0.4]);
    let (law, template) = multiple_motion_sources(&card);
    let state = law
        .advance_motion_consistent(
            &law.initial_state(),
            &drive([0.1, 0.0]),
            &template,
            EllipticFrictionMode::Kinetic,
            SymmetricTangentialCompliance::new([[0.0; 2]; 2]).unwrap(),
            solve([0.0, 0.2]),
            || false,
        )
        .unwrap()
        .step
        .state;
    let original = state.clone();
    let conditions = boundary(template, 320.0, Some(2e5));
    let first = law
        .solve_preimpact_impulse(
            &state,
            &collision([0.6, 0.8]),
            &conditions,
            SourceImpulseCoefficientRole::Kinetic,
            controls(),
            &mut || false,
        )
        .unwrap();
    close(first.queried_speed_m_s, 1.0);
    close(first.coefficients_mu[0], 0.4);
    close(first.coefficients_mu[1], 0.235);
    assert!((norm(first.impulse.velocity_after_m_s) - first.queried_speed_m_s).abs() > 0.1);
    for source in first.state.friction_sources() {
        for name in ["speed", "speed_aux"] {
            close(source.state_coordinate(name).unwrap(), 1.0);
        }
        for (name, expected) in [
            ("T", 320.0),
            ("Ta", 310.0),
            ("Tb", 330.0),
            ("p", 2e5),
            ("ambient", 295.0),
        ] {
            assert_eq!(source.state_coordinate(name), Some(expected));
        }
    }
    let second = law
        .solve_preimpact_impulse(
            &first.state,
            &collision(first.impulse.velocity_after_m_s),
            &conditions,
            SourceImpulseCoefficientRole::Kinetic,
            controls(),
            &mut || false,
        )
        .unwrap();
    close(
        second.queried_speed_m_s,
        norm(first.impulse.velocity_after_m_s),
    );
    close(
        second.coefficients_mu[0],
        0.3 + 0.1 * second.queried_speed_m_s,
    );
    for candidate in [&first.state, &second.state] {
        for (axis, source) in candidate.friction_sources().iter().enumerate() {
            for name in ["path", "path_aux"] {
                assert_eq!(
                    source.state_coordinate(name).unwrap().to_bits(),
                    state.friction_sources()[axis]
                        .state_coordinate(name)
                        .unwrap()
                        .to_bits()
                );
            }
        }
        unchanged_history(&state, candidate);
        original_claims(&card, &state, candidate);
    }
    assert_eq!(state, original);
}

#[test]
fn physical_temperature_pressure_and_fixed_query_ownership_remain_mandatory() {
    let card = motion_card([0.3, 0.5], [0.2, 0.4]);
    let law = motion_law(&card, [1.0, 0.0]);
    let state = law.initial_state();
    let original = state.clone();
    let good = boundary(endpoint(), 320.0, Some(2e5));
    let request = collision([0.6, 0.8]);
    let run = |conditions: &SourceFrictionConditions| {
        law.solve_preimpact_impulse(
            &state,
            &request,
            conditions,
            SourceImpulseCoefficientRole::Kinetic,
            controls(),
            &mut || false,
        )
    };
    for conditions in [
        SourceFrictionConditions {
            contact_temperature_k: None,
            ..good.clone()
        },
        SourceFrictionConditions {
            contact_temperature_k: Some(319.0),
            ..good.clone()
        },
        SourceFrictionConditions {
            nominal_pressure_pa: None,
            ..good.clone()
        },
        SourceFrictionConditions {
            nominal_pressure_pa: Some(1e5),
            ..good.clone()
        },
        SourceFrictionConditions {
            query: good.query.clone().with("ambient", 296.0).unwrap(),
            ..good.clone()
        },
        boundary(query(320.0), 320.0, Some(2e5)),
    ] {
        assert!(run(&conditions).is_err());
        assert_eq!(state, original);
    }
    let outside = boundary(
        good.query.clone().with("T", 321.0).unwrap(),
        321.0,
        Some(2e5),
    );
    assert!(matches!(
        run(&outside),
        Err(SourceImpulseError::Source(
            SourceThermalCoulombError::Material(_)
        ))
    ));
    let no_pressure = self::super::super::law(&directional_card(), [1.0, 0.0]);
    assert!(
        no_pressure
            .solve_preimpact_impulse(
                &no_pressure.initial_state(),
                &request,
                &boundary(query(300.0), 300.0, Some(1e5)),
                SourceImpulseCoefficientRole::Kinetic,
                controls(),
                &mut || false,
            )
            .is_err()
    );
    assert_eq!(state, original);
}

#[test]
fn second_source_failure_foreign_history_and_late_cancellation_publish_no_candidate() {
    let mut limited = curve(KINETIC_Y, [0.1, 0.5]);
    limited.validity = ValidityDomain::unconstrained().with("T", 280.0, 310.0);
    let card = card(vec![
        curve(STATIC, [0.8, 0.4]),
        curve(KINETIC, [0.6, 0.2]),
        curve(STATIC_Y, [0.5, 0.9]),
        limited,
    ]);
    let law = law(&card, [1.0, 0.0]);
    let state = law
        .advance_sourced(&law.initial_state(), &input(300.0, [0.01, 0.02]))
        .unwrap()
        .state;
    let saved = state.clone();
    let request = collision([1.0, 0.5]);
    let rejected = law.solve_preimpact_impulse(
        &state,
        &request,
        &boundary(query(320.0), 320.0, None),
        SourceImpulseCoefficientRole::Kinetic,
        controls(),
        &mut || false,
    );
    assert!(matches!(
        rejected,
        Err(SourceImpulseError::Source(
            SourceThermalCoulombError::Material(_)
        ))
    ));
    assert_eq!(state, saved);
    let conditions = boundary(query(305.0), 305.0, None);
    let foreign = self::super::super::law(&card, [0.0, 1.0]);
    assert!(
        foreign
            .solve_preimpact_impulse(
                &state,
                &request,
                &conditions,
                SourceImpulseCoefficientRole::Kinetic,
                controls(),
                &mut || false,
            )
            .is_err()
    );
    for target_poll in [1, 2, 3, 4, 5, 6, 7] {
        let mut polls = 0;
        let result = law.solve_preimpact_impulse(
            &state,
            &request,
            &conditions,
            SourceImpulseCoefficientRole::Kinetic,
            controls(),
            &mut || {
                polls += 1;
                polls == target_poll
            },
        );
        assert!(matches!(result, Err(SourceImpulseError::Cancelled)));
        assert_eq!(state, saved);
    }
    let first = law
        .solve_preimpact_impulse(
            &state,
            &request,
            &conditions,
            SourceImpulseCoefficientRole::Kinetic,
            controls(),
            &mut || false,
        )
        .unwrap();
    let retry = law
        .solve_preimpact_impulse(
            &state,
            &request,
            &conditions,
            SourceImpulseCoefficientRole::Kinetic,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(first, retry);
    assert_eq!(state, saved);
}

#[test]
fn kinetic_impact_preserves_the_full_static_limit_for_subsequent_supported_contact() {
    let card = directional_card();
    let law = law(&card, [1.0, 0.0]);
    let accepted = law
        .advance_sourced(&law.initial_state(), &input(300.0, [0.01, 0.02]))
        .unwrap()
        .state;
    let impact = law
        .solve_preimpact_impulse(
            &accepted,
            &collision([1.0, 0.0]),
            &boundary(query(320.0), 320.0, None),
            SourceImpulseCoefficientRole::Kinetic,
            controls(),
            &mut || false,
        )
        .unwrap();
    close(impact.coefficients_mu[0], 0.2);
    unchanged_history(&accepted, &impact.state);
    // Required support traction is3N: above kinetic2N and below static4N.
    let mut supported = input(320.0, [0.03, 0.0]);
    supported.compliance_m_per_n = 0.01;
    let held = law
        .advance_at_query(
            &impact.state,
            &supported,
            &query(320.0),
            EllipticFrictionMode::StaticOrSliding,
        )
        .unwrap();
    assert_eq!(held.regime, FrictionRegime::Sticking);
    close(held.traction_n[0], 3.0);
    close(held.dissipated_j, 0.0);
    assert_eq!(
        held.state.constitutive_state().slip_m(),
        accepted.constitutive_state().slip_m()
    );
    assert_eq!(
        held.state.constitutive_state().contact_temperature_k(),
        Some(320.0)
    );
    let sliding = law
        .advance_at_query(
            &impact.state,
            &supported,
            &query(320.0),
            EllipticFrictionMode::Kinetic,
        )
        .unwrap();
    assert_eq!(sliding.regime, FrictionRegime::Sliding);
    close(sliding.traction_n[0], 2.0);
    close(sliding.slip_increment_m[0], 0.01);
    close(sliding.dissipated_j, 0.02);
    original_claims(&card, &accepted, &held.state);
    original_claims(&card, &accepted, &sliding.state);
}
