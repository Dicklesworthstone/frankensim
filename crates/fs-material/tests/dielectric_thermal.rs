//! Synthetic Debye temperature laws: potential derivatives, parameter work and
//! complete source refresh. These are constitutive checks, not measured accuracy.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS,
        dielectric::{DebyePoleProperties, ResolvedDebyeCapacitor, thermal::ThermalDebyeCapacitor},
    },
    state_point::{MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement},
};
use fs_qty::{
    Dims, QuantitySpec, Temperature, Time,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            quantity(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with_quantity("p", quantity(QuantityKind::Pressure), 1e5)
        .unwrap()
}
fn path(dims: Dims, values: &[(f64, f64)]) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: values.to_vec(),
        dims,
    }
}
fn constant(dims: Dims, value: f64) -> PropertyValue {
    PropertyValue::Scalar { value, dims }
}
fn card(values: &[(&str, Dims, PropertyValue)]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in values {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(*name, *dims),
                value: value.clone(),
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        280.0,
                        340.0,
                    )
                    .with_quantity("p", quantity(QuantityKind::Pressure), 0.5e5, 2e5),
                interpolation: if matches!(value, PropertyValue::Scalar { .. }) {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic temperature-dependent Debye fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn names() -> Vec<DebyePoleProperties> {
    vec![
        DebyePoleProperties {
            permittivity_increment: "delta1".into(),
            relaxation_time: "tau1".into(),
        },
        DebyePoleProperties {
            permittivity_increment: "delta2".into(),
            relaxation_time: "tau2".into(),
        },
    ]
}
fn resolve(card: &MaterialCard, poles: usize, temperature: f64) -> ResolvedDebyeCapacitor {
    ResolvedDebyeCapacitor::resolve(
        card,
        &point(temperature),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &names()[..poles],
        0.02,
        0.005,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-300),
        "{actual} != {expected}"
    );
}
fn fixture() -> MaterialCard {
    card(&[
        (
            "epsilon_inf",
            ELECTRIC_PERMITTIVITY_DIMS,
            path(
                ELECTRIC_PERMITTIVITY_DIMS,
                &[(280.0, 1.8e-6), (340.0, 2.4e-6)],
            ),
        ),
        (
            "delta1",
            ELECTRIC_PERMITTIVITY_DIMS,
            path(
                ELECTRIC_PERMITTIVITY_DIMS,
                &[(280.0, 3.4e-6), (340.0, 2.2e-6)],
            ),
        ),
        (
            "tau1",
            Time::DIMS,
            path(Time::DIMS, &[(280.0, 0.02), (300.0, 0.01), (340.0, 0.04)]),
        ),
        (
            "delta2",
            ELECTRIC_PERMITTIVITY_DIMS,
            constant(ELECTRIC_PERMITTIVITY_DIMS, 5e-6),
        ),
        ("tau2", Time::DIMS, constant(Time::DIMS, 0.2)),
        (
            "caloric_marker",
            Dims([0; 6]),
            path(Dims([0; 6]), &[(280.0, 100.0), (340.0, 160.0)]),
        ),
    ])
}

#[test]
fn thermal_refresh_preserves_complete_material_and_original_pole_roles() {
    let card = fixture();
    let original = resolve(&card, 2, 300.0);
    let requirement = ScalarPropertyRequirement::try_new(
        "caloric_marker",
        Dims([0; 6]),
        ScalarAdmissibility::StrictlyPositive,
    )
    .unwrap();
    let complete = original
        .material()
        .with_additional_requirements(
            &card,
            &[requirement],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
    let original =
        ResolvedDebyeCapacitor::try_from_material(complete, "epsilon_inf", &names(), 0.02, 0.005)
            .unwrap();
    let thermal = ThermalDebyeCapacitor::new(&card, &original, "T", [280.0, 340.0]).unwrap();
    let refreshed = thermal.at_temperature(320.0).unwrap();
    let source = refreshed.source();
    assert_eq!(source.material().card_identity(), card.content_hash());
    assert_eq!(source.material().properties().len(), 6);
    assert_eq!(source.material().state_coordinate("p"), Some(1e5));
    assert_eq!(source.pole_properties(), original.pole_properties());
    assert_eq!((source.area_m2(), source.gap_m()), (0.02, 0.005));
    close(
        source
            .material()
            .property("caloric_marker")
            .unwrap()
            .value_si(),
        140.0,
        1e-14,
    );
    close(source.poles()[0].relaxation_time_s, 0.025, 1e-14);
    close(
        source.poles()[0].resistance_ohm * source.poles()[0].capacitance_f,
        0.025,
        1e-14,
    );
    for (before, after) in original
        .material()
        .properties()
        .iter()
        .zip(source.material().properties())
    {
        assert_eq!(
            before.answer().receipt.selected,
            after.answer().receipt.selected
        );
        assert_eq!(before.requirement(), after.requirement());
    }
    assert_eq!(thermal.at_source(source).unwrap(), refreshed);
}

#[test]
fn independent_affine_potential_derivatives_cover_every_polarization_charge() {
    let card = fixture();
    let thermal =
        ThermalDebyeCapacitor::new(&card, &resolve(&card, 2, 300.0), "T", [280.0, 340.0]).unwrap();
    let temperature = 315.0;
    let charges = [1e-5, -2e-5, 3e-5];
    let source = thermal.at_temperature(temperature).unwrap();
    let actual = source.state(charges[0], &charges[1..]).unwrap();
    let coefficients = [
        4.0 * (2e-6 + 1e-8 * 15.0),
        4.0 * (3e-6 - 2e-8 * 15.0),
        20e-6,
    ];
    let derivatives = [4e-8, -8e-8, 0.0];
    let mut psi = 0.0;
    let mut entropy = 0.0;
    let mut capacity = 0.0;
    for ((charge, c), slope) in charges.into_iter().zip(coefficients).zip(derivatives) {
        psi += charge * charge / (2.0 * c);
        entropy += charge * charge * slope / (2.0 * c * c);
        capacity -= temperature * charge * charge * slope * slope / (c * c * c);
    }
    close(actual.free_energy_j, psi, 1e-13);
    close(actual.entropy_j_k, entropy, 1e-13);
    close(actual.internal_energy_j, psi + temperature * entropy, 1e-13);
    close(actual.heat_capacity_j_k, capacity, 1e-13);
    let low = thermal
        .at_temperature(temperature - 0.01)
        .unwrap()
        .state(charges[0], &charges[1..])
        .unwrap();
    let high = thermal
        .at_temperature(temperature + 0.01)
        .unwrap()
        .state(charges[0], &charges[1..])
        .unwrap();
    close(
        -(high.free_energy_j - low.free_energy_j) / 0.02,
        entropy,
        2e-8,
    );
    close(
        (high.internal_energy_j - low.internal_energy_j) / 0.02,
        capacity,
        2e-8,
    );
    let zero = source.state(0.0, &[0.0, 0.0]).unwrap();
    assert_eq!(
        [
            zero.free_energy_j,
            zero.entropy_j_k,
            zero.internal_energy_j,
            zero.heat_capacity_j_k
        ],
        [0.0; 4]
    );
    let bound = thermal
        .minimum_heat_capacity_j_k(2e-5, &[3e-5, 4e-5])
        .unwrap();
    assert!(bound < 0.0);
    for i in 0..=30 {
        let at = thermal.at_temperature(280.0 + 2.0 * f64::from(i)).unwrap();
        assert!(at.state(2e-5, &[3e-5, 4e-5]).unwrap().heat_capacity_j_k >= bound);
    }
}

#[test]
fn parameter_work_survives_unchanged_large_instantaneous_energy() {
    let card = card(&[
        (
            "epsilon_inf",
            ELECTRIC_PERMITTIVITY_DIMS,
            constant(ELECTRIC_PERMITTIVITY_DIMS, 0.25),
        ),
        (
            "delta1",
            ELECTRIC_PERMITTIVITY_DIMS,
            path(ELECTRIC_PERMITTIVITY_DIMS, &[(280.0, 0.5), (340.0, 1.0)]),
        ),
        ("tau1", Time::DIMS, constant(Time::DIMS, 0.1)),
    ]);
    let initial = resolve(&card, 1, 280.0);
    let heated = initial.at_query_point(&card, &point(340.0)).unwrap();
    assert_eq!(
        initial.stored_energy_j(1e8, &[1.0]).unwrap(),
        heated.stored_energy_j(1e8, &[1.0]).unwrap()
    );
    let change = initial
        .fixed_charge_free_energy_change_j(&heated, 1e8, &[1.0])
        .unwrap();
    close(change, -0.125, 1e-14);
    assert_eq!(
        heated
            .fixed_charge_free_energy_change_j(&initial, 1e8, &[1.0])
            .unwrap(),
        -change
    );
    assert_eq!(
        initial
            .fixed_charge_free_energy_change_j(&initial, 1e8, &[1.0])
            .unwrap(),
        0.0
    );
}

#[test]
fn thermal_offset_retains_small_charge_work_and_both_temperature_directions() {
    let source_card = fixture();
    let law = ThermalDebyeCapacitor::new(
        &source_card,
        &resolve(&source_card, 2, 300.0),
        "T",
        [280.0, 340.0],
    )
    .unwrap();
    let a = law.at_temperature(300.0).unwrap();
    let b = law.at_temperature(328.0).unwrap();
    let q0 = [1e-5, -2e-5, 3e-5];
    let q1 = [1.1e-5, -1.8e-5, 2.9e-5];
    let before = a.state(q0[0], &q0[1..]).unwrap();
    let after = b.state(q1[0], &q1[1..]).unwrap();
    let expected = 328.0 * after.entropy_j_k - 300.0 * before.entropy_j_k;
    let change = a
        .thermal_offset_change_j(&b, q0[0], &q0[1..], q1[0], &q1[1..])
        .unwrap();
    close(change, expected, 5e-13);
    close(
        b.thermal_offset_change_j(&a, q1[0], &q1[1..], q0[0], &q0[1..])
            .unwrap(),
        -change,
        5e-13,
    );

    let biased_card = card(&[
        (
            "epsilon_inf",
            ELECTRIC_PERMITTIVITY_DIMS,
            path(
                ELECTRIC_PERMITTIVITY_DIMS,
                &[(280.0, 0.245), (340.0, 0.260)],
            ),
        ),
        (
            "delta1",
            ELECTRIC_PERMITTIVITY_DIMS,
            path(ELECTRIC_PERMITTIVITY_DIMS, &[(280.0, 0.45), (340.0, 0.60)]),
        ),
        ("tau1", Time::DIMS, constant(Time::DIMS, 0.1)),
    ]);
    let biased_law = ThermalDebyeCapacitor::new(
        &biased_card,
        &resolve(&biased_card, 1, 300.0),
        "T",
        [280.0, 340.0],
    )
    .unwrap();
    let at = biased_law.at_temperature(300.0).unwrap();
    let old_polarization = 1.0;
    let new_polarization = 1.000001;
    let old = at.state(1e8, &[old_polarization]).unwrap();
    let new = at.state(1e8, &[new_polarization]).unwrap();
    assert_eq!(300.0 * old.entropy_j_k, 300.0 * new.entropy_j_k);
    let resolved = at
        .thermal_offset_change_j(&at, 1e8, &[old_polarization], 1e8, &[new_polarization])
        .unwrap();
    // Independent analytic entropy law S(q)=q² C'/(2 C²) for the small pole.
    let expected = (300.0 * 0.01 / 8.0)
        * (new_polarization - old_polarization)
        * (new_polarization + old_polarization);
    assert!(resolved > 0.0);
    close(resolved, expected, 2e-13);
    assert!(
        at.thermal_offset_change_j(&at, 1e8, &[], 1e8, &[new_polarization])
            .is_err()
    );
}

#[test]
fn relaxation_temperature_dependence_is_friction_not_storage_entropy() {
    let card = card(&[
        (
            "epsilon_inf",
            ELECTRIC_PERMITTIVITY_DIMS,
            constant(ELECTRIC_PERMITTIVITY_DIMS, 2e-6),
        ),
        (
            "delta1",
            ELECTRIC_PERMITTIVITY_DIMS,
            constant(ELECTRIC_PERMITTIVITY_DIMS, 3e-6),
        ),
        (
            "tau1",
            Time::DIMS,
            path(Time::DIMS, &[(280.0, 0.02), (300.0, 0.01), (340.0, 0.2)]),
        ),
    ]);
    for poles in [0, 1] {
        let original = resolve(&card, poles, 300.0);
        let law = ThermalDebyeCapacitor::new(&card, &original, "T", [280.0, 340.0]).unwrap();
        let high = law.at_temperature(340.0).unwrap();
        let q = [2e-5];
        let state = high.state(1e-5, &q[..poles]).unwrap();
        assert_eq!(state.entropy_j_k, 0.0);
        assert_eq!(state.heat_capacity_j_k, 0.0);
        assert_eq!(state.internal_energy_j, state.free_energy_j);
        assert_eq!(
            law.minimum_heat_capacity_j_k(1e-4, &q[..poles]).unwrap(),
            0.0
        );
        assert_eq!(
            original
                .fixed_charge_free_energy_change_j(high.source(), 1e-5, &q[..poles])
                .unwrap(),
            0.0
        );
        if poles == 1 {
            close(
                high.source().poles()[0].resistance_ohm / original.poles()[0].resistance_ohm,
                20.0,
                1e-13,
            );
        }
    }
}

#[test]
fn thermal_source_domain_geometry_pole_history_and_nonthermal_changes_refuse() {
    let card = fixture();
    let original = resolve(&card, 2, 300.0);
    let law = ThermalDebyeCapacitor::new(&card, &original, "T", [280.0, 340.0]).unwrap();
    for bad in [279.0, 341.0, f64::NAN] {
        assert!(law.at_temperature(bad).is_err());
    }
    assert!(
        law.at_temperature(310.0)
            .unwrap()
            .state(0.0, &[0.0])
            .is_err()
    );
    assert!(law.minimum_heat_capacity_j_k(1.0, &[1.0]).is_err());
    assert!(law.minimum_heat_capacity_j_k(-1.0, &[1.0, 1.0]).is_err());
    let changed_geometry = ResolvedDebyeCapacitor::try_from_material(
        original.material().clone(),
        "epsilon_inf",
        &names(),
        0.04,
        0.01,
    )
    .unwrap();
    assert!(law.at_source(&changed_geometry).is_err());
    assert!(
        original
            .fixed_charge_free_energy_change_j(&changed_geometry, 1.0, &[1.0, 1.0])
            .is_err()
    );
    let mut reversed = names();
    reversed.reverse();
    let changed_roles = ResolvedDebyeCapacitor::try_from_material(
        original.material().clone(),
        "epsilon_inf",
        &reversed,
        0.02,
        0.005,
    )
    .unwrap();
    assert!(law.at_source(&changed_roles).is_err());
    let pressure = original
        .at_query_point(&card, &point(300.0).with("p", 1.2e5).unwrap())
        .unwrap();
    assert!(law.at_source(&pressure).is_err());
    assert!(
        original
            .fixed_charge_free_energy_change_j(&pressure, 1.0, &[1.0, 1.0])
            .is_err()
    );
    let wrong_card = self::card(&[(
        "epsilon_inf",
        ELECTRIC_PERMITTIVITY_DIMS,
        constant(ELECTRIC_PERMITTIVITY_DIMS, 2e-6),
    )]);
    assert!(original.at_query_point(&wrong_card, &point(300.0)).is_err());
}

#[test]
fn interior_storage_kinks_and_nonpassive_relaxation_are_not_hidden_by_endpoints() {
    let attempt = |epsilon, tau| {
        let card = card(&[
            ("epsilon_inf", ELECTRIC_PERMITTIVITY_DIMS, epsilon),
            (
                "delta1",
                ELECTRIC_PERMITTIVITY_DIMS,
                constant(ELECTRIC_PERMITTIVITY_DIMS, 1e-6),
            ),
            ("tau1", Time::DIMS, tau),
        ]);
        ThermalDebyeCapacitor::new(&card, &resolve(&card, 1, 300.0), "T", [280.0, 340.0])
    };
    assert!(
        attempt(
            path(
                ELECTRIC_PERMITTIVITY_DIMS,
                &[(280.0, 1e-6), (300.0, 2e-6), (340.0, 3e-6)]
            ),
            constant(Time::DIMS, 0.1)
        )
        .is_err()
    );
    assert!(
        attempt(
            constant(ELECTRIC_PERMITTIVITY_DIMS, 2e-6),
            path(
                Time::DIMS,
                &[(280.0, 0.1), (300.0, 0.1), (320.0, -0.01), (340.0, 0.1)]
            )
        )
        .is_err()
    );
}
