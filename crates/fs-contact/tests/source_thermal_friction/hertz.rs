//! Source-to-Hertz integration using the same original interface fixtures.

use super::*;
use fs_contact::interface_binding::thermal_friction::hertz::SourceHertzWrenchError;
use fs_contact::interface_binding::thermal_friction::{FrictionQueryAxis, FrictionQueryAxisKind};
use fs_tribo::sliding_patch::contact::{
    HertzWrenchConfig, HertzWrenchError, HertzWrenchInput, HertzWrenchMode,
};

fn solve_config() -> HertzWrenchConfig {
    HertzWrenchConfig {
        maximum_static_iterations: 96,
        maximum_iterations: 256,
        maximum_evaluations: 512,
        force_tolerance_n: 2e-9,
        moment_tolerance_nm: 2e-10,
        displacement_tolerance_m: 2e-13,
        energy_tolerance_j: 2e-11,
        relative_matrix_tolerance: 1e-13,
        maximum_force_n: 30.0,
        maximum_spin_moment_nm: 3.0,
    }
}

fn released() -> HertzWrenchInput {
    HertzWrenchInput {
        free_displacement_m: [0.003, -0.004],
        free_spin_angle_rad: 0.25,
        scaled_compliance_m_per_n: [[0.0; 3]; 3],
    }
}

#[test]
fn actual_temperature_drives_joint_force_spin_heat_and_original_receipts() {
    let law = resolve_at(&source_card(), query(300.0)).unwrap();
    let saved = law.clone();
    let cold = law.hertz_at(280.0, None, 10.0, 0.02).unwrap();
    let warm = law.hertz_at(320.0, None, 10.0, 0.02).unwrap();
    for (binding, temperature, static_mu, kinetic_mu) in
        [(&cold, 280.0, 0.8, 0.6), (&warm, 320.0, 0.4, 0.2)]
    {
        assert_eq!(binding.owner_identity(), law.identity());
        let parameters = binding.law().parameters();
        close(parameters.static_mu, static_mu);
        close(parameters.kinetic_mu, kinetic_mu);
        close(parameters.normal_force_n, 10.0);
        close(parameters.contact_radius_m, 0.02);
        assert_eq!(
            binding.friction_source().state_coordinate("T"),
            Some(temperature)
        );
        for property in binding.friction_source().properties() {
            let original = law
                .initial_source()
                .property(property.requirement().name())
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                original.answer().receipt.selected
            );
            assert_eq!(
                property.answer().evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            law.interface()
                .card()
                .claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
    let cold_step = cold
        .law()
        .advance(&released(), solve_config(), || false)
        .unwrap();
    let warm_step = warm
        .law()
        .advance(&released(), solve_config(), || false)
        .unwrap();
    // At epsilon=1 the independent circular-Hertz limits are 9*pi/32 and
    // 9*pi/128 for force/(mu*N) and moment/(mu*N*a), respectively.
    let force_ratio = 9.0 * core::f64::consts::PI / 32.0;
    let moment_ratio = 9.0 * core::f64::consts::PI / 128.0;
    close(cold_step.force_n[0], -6.0 * force_ratio * 0.6);
    close(cold_step.force_n[1], 6.0 * force_ratio * 0.8);
    close(cold_step.spin_moment_nm, -6.0 * 0.02 * moment_ratio);
    close(
        cold_step.dissipated_energy_j,
        6.0 * force_ratio * 0.005 + 6.0 * 0.02 * moment_ratio * 0.25,
    );
    for (c, w) in cold_step.force_n.into_iter().zip(warm_step.force_n) {
        close(c, 3.0 * w);
    }
    close(cold_step.spin_moment_nm, 3.0 * warm_step.spin_moment_nm);
    close(
        cold_step.dissipated_energy_j,
        3.0 * warm_step.dissipated_energy_j,
    );
    assert_eq!(law, saved);
    assert_eq!(
        law.interface().interface().provenance().authority(),
        InputAuthority::SyntheticFixture
    );
}

#[test]
fn sourced_thermal_weakening_breaks_joint_static_contact_with_full_cross_response() {
    let law = resolve_at(&source_card(), query(300.0)).unwrap();
    let cold = law.hertz_at(280.0, None, 10.0, 0.02).unwrap();
    let warm = law.hertz_at(320.0, None, 10.0, 0.02).unwrap();
    let c = [
        [0.001, 0.0002, 0.0001],
        [0.0002, 0.002, 0.0003],
        [0.0001, 0.0003, 0.003],
    ];
    let reaction = [5.0, 0.0, 0.5];
    let free = c.map(|row| {
        row.into_iter()
            .zip(reaction)
            .map(|(a, b)| a * b)
            .sum::<f64>()
    });
    let input = HertzWrenchInput {
        free_displacement_m: [free[0], free[1]],
        free_spin_angle_rad: free[2] / 0.02,
        scaled_compliance_m_per_n: c,
    };
    let held = cold
        .law()
        .advance(&input, solve_config(), || false)
        .unwrap();
    assert_eq!(held.mode, HertzWrenchMode::Sticking);
    close(held.force_n[0], -5.0);
    close(held.force_n[1], 0.0);
    close(held.spin_moment_nm, -0.01);
    assert_eq!(held.dissipated_energy_j, 0.0);
    let slide = warm
        .law()
        .advance(&input, solve_config(), || false)
        .unwrap();
    assert_eq!(slide.mode, HertzWrenchMode::Sliding);
    assert!(slide.dissipated_energy_j > 0.0);
    assert!(slide.force_n[0].hypot(slide.force_n[1]) <= 2.0 + 2e-12);
    let actual = [
        slide.slip_displacement_m[0],
        slide.slip_displacement_m[1],
        0.02 * slide.spin_angle_rad,
    ];
    let force = [
        slide.force_n[0],
        slide.force_n[1],
        slide.spin_moment_nm / 0.02,
    ];
    for i in 0..3 {
        close(
            actual[i],
            free[i] + c[i].into_iter().zip(force).map(|(a, b)| a * b).sum::<f64>(),
        );
    }
    assert!(slide.work_residual_j.abs() <= solve_config().energy_tolerance_j);
    assert_eq!(cold.owner_identity(), warm.owner_identity());
    assert_ne!(
        cold.friction_source().query(),
        warm.friction_source().query()
    );
}

fn pressure_owner() -> SourceThermalCoulomb {
    let card = card(vec![
        pressure_surface(STATIC, 0.0),
        pressure_surface(KINETIC, 0.2),
    ]);
    SourceThermalCoulomb::resolve(
        bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
        query(300.0)
            .with_quantity("p", pressure_quantity(), 1e5)
            .unwrap(),
        requirement(STATIC),
        requirement(KINETIC),
        MaterialPropertySelection::SingleClaimOnly,
        ThermalFrictionAxes {
            temperature: "T".into(),
            nominal_pressure: Some("p".into()),
        },
    )
    .unwrap()
}

#[test]
fn nominal_pressure_is_an_explicit_source_coordinate_independent_of_normal_geometry() {
    let law = pressure_owner();
    let first = law.hertz_at(300.0, Some(1e5), 10.0, 0.02).unwrap();
    let loaded = law.hertz_at(300.0, Some(1e5), 20.0, 0.04).unwrap();
    let pressure_changed = law.hertz_at(300.0, Some(3e5), 10.0, 0.02).unwrap();
    close(first.law().parameters().kinetic_mu, 0.4);
    close(loaded.law().parameters().kinetic_mu, 0.4);
    close(pressure_changed.law().parameters().kinetic_mu, 0.6);
    assert_eq!(first.friction_source(), loaded.friction_source());
    assert_eq!(
        pressure_changed.friction_source().state_coordinate("p"),
        Some(3e5)
    );
    // Double radius and half spin angle keep the same local slip/spin ratio.
    let a = first
        .law()
        .advance(&released(), solve_config(), || false)
        .unwrap();
    let mut changed_motion = released();
    changed_motion.free_spin_angle_rad *= 0.5;
    let b = loaded
        .law()
        .advance(&changed_motion, solve_config(), || false)
        .unwrap();
    close(b.force_n[0], 2.0 * a.force_n[0]);
    close(b.force_n[1], 2.0 * a.force_n[1]);
    close(b.spin_moment_nm, 4.0 * a.spin_moment_nm);
    close(b.dissipated_energy_j, 2.0 * a.dissipated_energy_j);
    for pressure in [None, Some(-1.0), Some(4e5), Some(f64::NAN)] {
        assert!(law.hertz_at(300.0, pressure, 10.0, 0.02).is_err());
    }
    let zero = law.hertz_at(300.0, Some(1e5), 0.0, 0.02).unwrap();
    let step = zero
        .law()
        .advance(&released(), solve_config(), || false)
        .unwrap();
    assert_eq!(step.mode, HertzWrenchMode::Frictionless);
    assert_eq!(step.dissipated_energy_j, 0.0);
    assert!(law.hertz_at(300.0, None, 0.0, 0.02).is_err());
}

#[test]
fn typed_temperature_alias_and_fixed_external_conditions_preserve_original_owner() {
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let properties = [STATIC, KINETIC].map(|name| {
        let mut p = curve(name, [0.6, 0.2]);
        if let PropertyValue::Curve { abscissa, .. } = &mut p.value {
            *abscissa = "contact-temperature".into();
        }
        p.validity = ValidityDomain::unconstrained()
            .with_quantity("contact-temperature", absolute, 280.0, 320.0)
            .with_quantity("ambient", pressure_quantity(), 8e4, 1.2e5);
        p
    });
    let card = card(properties.into());
    let point = QueryPoint::new()
        .with_quantity("contact-temperature", absolute, 300.0)
        .unwrap()
        .with_quantity("ambient", pressure_quantity(), 1e5)
        .unwrap();
    let axes = ThermalFrictionAxes {
        temperature: "contact-temperature".into(),
        nominal_pressure: None,
    };
    let make = |authority| {
        SourceThermalCoulomb::resolve(
            bind_dry_interface_system_card(&card, authority).unwrap(),
            point.clone(),
            requirement(STATIC),
            requirement(KINETIC),
            MaterialPropertySelection::SingleClaimOnly,
            axes.clone(),
        )
        .unwrap()
    };
    let original = make(InputAuthority::SyntheticFixture);
    let binding = original.hertz_at(320.0, None, 10.0, 0.02).unwrap();
    close(binding.law().parameters().static_mu, 0.2);
    close(binding.law().parameters().kinetic_mu, 0.2);
    assert_eq!(
        binding
            .friction_source()
            .state_coordinate("contact-temperature"),
        Some(320.0)
    );
    assert_eq!(
        binding.friction_source().state_coordinate("ambient"),
        Some(1e5)
    );
    assert_eq!(binding.owner_identity(), original.identity());
    assert_ne!(
        binding.owner_identity(),
        make(InputAuthority::Estimated).identity()
    );
    assert!(original.hertz_at(320.0, Some(1e5), 10.0, 0.02).is_err());
}

#[test]
fn pinned_claims_cannot_be_substituted_and_all_refusals_preserve_exact_retry() {
    let mut claims = vec![];
    for name in [STATIC, KINETIC] {
        claims.push(curve(name, [0.6, 0.3]));
        let mut other = curve(name, [0.2, 0.1]);
        other.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(330.0, 0.2), (360.0, 0.1)],
            dims: Dims::NONE,
        };
        other.validity = ValidityDomain::unconstrained().with("T", 330.0, 360.0);
        other.provenance.source = "other high-temperature Hertz coefficient source".into();
        claims.push(other);
    }
    let card = card(claims);
    let pins = |high| {
        MaterialPropertySelection::PinnedByProperty(
            [STATIC, KINETIC]
                .into_iter()
                .map(|name| {
                    (
                        name.into(),
                        card.claims_for(name)
                            .into_iter()
                            .find(|(_, c)| c.provenance.source.starts_with("other") == high)
                            .unwrap()
                            .0,
                    )
                })
                .collect(),
        )
    };
    let law = SourceThermalCoulomb::resolve(
        bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
        query(300.0),
        requirement(STATIC),
        requirement(KINETIC),
        pins(false),
        axes(),
    )
    .unwrap();
    assert!(
        resolve_interface_state_point(
            &card,
            &query(350.0),
            &[requirement(STATIC), requirement(KINETIC)],
            pins(true)
        )
        .is_ok()
    );
    let saved = law.clone();
    let expected = law.hertz_at(310.0, None, 10.0, 0.02).unwrap();
    assert!(matches!(
        law.hertz_at(350.0, None, 10.0, 0.02),
        Err(SourceHertzWrenchError::Source(
            SourceThermalCoulombError::Material(_)
        ))
    ));
    for temperature in [0.0, f64::NAN, f64::INFINITY] {
        assert!(law.hertz_at(temperature, None, 10.0, 0.02).is_err());
    }
    for (normal, radius) in [(-1.0, 0.02), (10.0, 0.0), (10.0, f64::NAN)] {
        assert!(matches!(
            law.hertz_at(310.0, None, normal, radius),
            Err(SourceHertzWrenchError::Hertz(_))
        ));
    }
    let receipt = expected.clone();
    assert_eq!(
        expected.law().advance(&released(), solve_config(), || true),
        Err(HertzWrenchError::Cancelled)
    );
    let mut limited = solve_config();
    limited.maximum_force_n = 0.01;
    assert!(matches!(
        expected.law().advance(&released(), limited, || false),
        Err(HertzWrenchError::EffortBudget { .. })
    ));
    assert_eq!(expected, receipt);
    assert_eq!(law, saved);
    assert_eq!(law.hertz_at(310.0, None, 10.0, 0.02).unwrap(), expected);
    let crossing = resolve_at(
        &super::card(vec![curve(STATIC, [0.8, 0.1]), curve(KINETIC, [0.6, 0.2])]),
        query(280.0),
    )
    .unwrap();
    assert!(matches!(
        crossing.hertz_at(320.0, None, 10.0, 0.02),
        Err(SourceHertzWrenchError::Source(
            SourceThermalCoulombError::Tribo(_)
        ))
    ));
}

#[test]
fn temperature_pressure_entry_refuses_owned_speed_without_freezing_complete_query_behavior() {
    let speed_dims = Dims([1, 0, -1, 0, 0, 0]);
    let speed_quantity = QuantitySpec::dimensional(speed_dims);
    let mut kinetic = curve(KINETIC, [0.6, 0.2]);
    kinetic.value = PropertyValue::Curve {
        abscissa: "speed".into(),
        abscissa_dims: speed_dims,
        knots: vec![(0.0, 0.6), (2.0, 0.2)],
        dims: Dims::NONE,
    };
    kinetic.validity = kinetic
        .validity
        .with_quantity("speed", speed_quantity, 0.0, 2.0);
    let card = card(vec![curve(STATIC, [0.8; 2]), kinetic]);
    let initial_query = query(300.0)
        .with_quantity("speed", speed_quantity, 0.5)
        .unwrap();
    let law = SourceThermalCoulomb::resolve_with_query_axes(
        bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
        initial_query.clone(),
        requirement(STATIC),
        requirement(KINETIC),
        MaterialPropertySelection::SingleClaimOnly,
        axes(),
        vec![FrictionQueryAxis {
            name: "speed".into(),
            kind: FrictionQueryAxisKind::SlidingSpeed,
        }],
    )
    .unwrap();
    let checkpoint = law.initial_state();
    let saved = checkpoint.clone();
    let mechanical = input(300.0, [0.0075, 0.0]);
    let actual_query = initial_query.with("speed", 0.75).unwrap();
    let expected = law
        .advance_at_query(&checkpoint, &mechanical, &actual_query)
        .unwrap();
    close(expected.traction_n[0], 4.5);
    for normal in [0.0, 10.0] {
        assert!(matches!(
            law.hertz_at(300.0, None, normal, 0.02),
            Err(SourceHertzWrenchError::Source(
                SourceThermalCoulombError::Invalid(
                    "owned friction coordinates require the complete physical query"
                )
            ))
        ));
    }
    assert_eq!(checkpoint, saved);
    assert_eq!(
        law.advance_at_query(&checkpoint, &mechanical, &actual_query)
            .unwrap(),
        expected
    );
    assert_eq!(
        expected.state.friction_source().state_coordinate("speed"),
        Some(0.75)
    );
}
