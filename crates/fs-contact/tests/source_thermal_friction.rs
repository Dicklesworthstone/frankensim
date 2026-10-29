//! G1/G3: original interface claims drive temperature/pressure friction and heat.

use fs_contact::interface_binding::bind_dry_interface_system_card;
use fs_contact::interface_binding::thermal_friction::{
    SourceThermalCoulomb, SourceThermalCoulombError, ThermalFrictionAxes,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_interface_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_tribo::{
    FrictionRegime, InputAuthority, TangentialLaw, TangentialSolveConfig, TangentialStepInput,
};

const STATIC: &str = "explicit-static-coefficient";
const KINETIC: &str = "explicit-kinetic-coefficient";

fn curve(name: &str, values: [f64; 2]) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, Dims::NONE),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(280.0, values[0]), (320.0, values[1])],
            dims: Dims::NONE,
        },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: format!("synthetic thermal Coulomb fixture: {name}; not measured data"),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn card(properties: Vec<PropertyClaim>) -> InterfaceSystemCard {
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "synthetic friction fixture".into(),
            revision: 0,
        },
        texture_frame: format!("texture-{name}"),
    };
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic fixture".into(),
            history: "virgin specimen".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn source_card() -> InterfaceSystemCard {
    card(vec![curve(STATIC, [0.8, 0.4]), curve(KINETIC, [0.6, 0.2])])
}

fn requirement(name: &str) -> ScalarPropertyRequirement {
    ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative).unwrap()
}

fn axes() -> ThermalFrictionAxes {
    ThermalFrictionAxes {
        temperature: "T".into(),
        nominal_pressure: None,
    }
}

fn query(temperature: f64) -> QueryPoint {
    QueryPoint::new().with("T", temperature).unwrap()
}

fn resolve_at(
    card: &InterfaceSystemCard,
    point: QueryPoint,
) -> Result<SourceThermalCoulomb, SourceThermalCoulombError> {
    SourceThermalCoulomb::resolve(
        bind_dry_interface_system_card(card, InputAuthority::SyntheticFixture).unwrap(),
        point,
        requirement(STATIC),
        requirement(KINETIC),
        MaterialPropertySelection::SingleClaimOnly,
        axes(),
    )
}

fn input(temperature: f64, slip: [f64; 2]) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: 10.0,
        slip_increment_m: slip,
        dt_s: 0.01,
        contact_temperature_k: Some(temperature),
        nominal_pressure_pa: None,
        compliance_m_per_n: 0.0,
        solve: TangentialSolveConfig {
            max_iterations: 64,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-12,
            force_relative_tolerance: 1e-12,
        },
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
        "actual={actual}, expected={expected}"
    );
}

fn pressure_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}

#[test]
fn original_temperature_curves_change_force_and_heat_without_resetting_slip() {
    let law = resolve_at(&source_card(), query(300.0)).unwrap();
    let cold = law
        .advance(&law.initial_state(), &input(300.0, [0.003, 0.004]))
        .unwrap();
    close(cold.traction_n[0], 2.4);
    close(cold.traction_n[1], 3.2);
    close(cold.dissipated_j, 0.02);
    let warm = law
        .advance(&cold.state, &input(320.0, [-0.006, -0.008]))
        .unwrap();
    close(warm.traction_n[0], -1.2);
    close(warm.traction_n[1], -1.6);
    close(warm.dissipated_j, 0.02);
    close(warm.state.constitutive_state().slip_m()[0], -0.003);
    close(warm.state.constitutive_state().slip_m()[1], -0.004);
    assert_eq!(
        warm.state.constitutive_state().contact_temperature_k(),
        Some(320.0)
    );
    assert_eq!(warm.stored_change_j, 0.0);
    assert_eq!(law.stored_energy_j(&warm.state), 0.0);
    close(warm.closure_residual_j, 0.0);
    close(
        warm.dissipated_j,
        warm.traction_n[0] * warm.slip_increment_m[0]
            + warm.traction_n[1] * warm.slip_increment_m[1],
    );
    for property in warm.state.friction_source().properties() {
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
    assert_eq!(
        law.interface().interface().provenance().authority(),
        InputAuthority::SyntheticFixture
    );
    let stationary = law.advance(&warm.state, &input(280.0, [0.0; 2])).unwrap();
    assert_eq!(stationary.traction_n, [0.0; 2]);
    assert_eq!(stationary.dissipated_j, 0.0);
    assert_eq!(
        stationary.state.constitutive_state().slip_m(),
        warm.state.constitutive_state().slip_m()
    );
}

#[test]
fn temperature_weakening_causes_native_static_breakaway_and_correct_actual_work() {
    let law = resolve_at(&source_card(), query(280.0)).unwrap();
    let mut trial = input(280.0, [0.005, 0.0]);
    trial.compliance_m_per_n = 0.001;
    let held = law.advance(&law.initial_state(), &trial).unwrap();
    assert_eq!(held.regime, FrictionRegime::Sticking);
    close(held.traction_n[0], 5.0);
    assert_eq!(held.slip_increment_m, [0.0; 2]);
    assert_eq!(held.dissipated_j, 0.0);
    trial.contact_temperature_k = Some(320.0);
    let broken = law.advance(&held.state, &trial).unwrap();
    assert_eq!(broken.regime, FrictionRegime::Sliding);
    close(broken.traction_n[0], 2.0);
    close(broken.slip_increment_m[0], 0.003);
    close(broken.dissipated_j, 0.006);
    close(broken.compatibility_residual_m[0], 0.0);
}

#[test]
fn objective_transport_preserves_original_source_and_refreshes_only_on_physical_advance() {
    use fs_tribo::partial_slip::TangentFrame;
    use fs_tribo::tangential::transport::{ObjectiveTangentialLaw, TangentialHistoryRotation};

    let law = resolve_at(&source_card(), query(300.0)).unwrap();
    let loaded = law
        .advance(&law.initial_state(), &input(320.0, [0.003, 0.004]))
        .unwrap();
    let saved = loaded.state.clone();
    let old_frame = TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
    let new_frame = TangentFrame::new([0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap();
    let identity = TangentialHistoryRotation::between(old_frame, old_frame, 0.0).unwrap();
    assert_eq!(law.transport_state(&loaded.state, identity).unwrap(), saved);
    let rotation = TangentialHistoryRotation::between(old_frame, new_frame, 0.0).unwrap();
    let rotated = law.transport_state(&loaded.state, rotation).unwrap();
    close(rotated.constitutive_state().slip_m()[0], 0.004);
    close(rotated.constitutive_state().slip_m()[1], -0.003);
    assert_eq!(rotated.friction_source(), saved.friction_source());
    assert_eq!(
        rotated.constitutive_state().contact_temperature_k(),
        Some(320.0)
    );
    assert_eq!(law.stored_energy_j(&rotated), 0.0);

    // The next physical step, not orientation transport, queries the cold source.
    let continued = law
        .advance(&rotated, &input(280.0, [0.004, -0.003]))
        .unwrap();
    close(continued.traction_n[0], 4.8);
    close(continued.traction_n[1], -3.6);
    close(continued.dissipated_j, 0.03);
    close(continued.state.constitutive_state().slip_m()[0], 0.008);
    close(continued.state.constitutive_state().slip_m()[1], -0.006);
    assert_eq!(
        continued.state.friction_source().state_coordinate("T"),
        Some(280.0)
    );
    for property in continued.state.friction_source().properties() {
        assert_eq!(
            property.answer().receipt.selected,
            saved
                .friction_source()
                .property(property.requirement().name())
                .unwrap()
                .answer()
                .receipt
                .selected
        );
        law.interface()
            .card()
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    let foreign = resolve_at(
        &card(vec![curve(STATIC, [0.9, 0.5]), curve(KINETIC, [0.6, 0.2])]),
        query(300.0),
    )
    .unwrap();
    assert!(foreign.transport_state(&loaded.state, rotation).is_err());
    assert_eq!(loaded.state, saved);
    assert_eq!(
        law.transport_state(&loaded.state, rotation).unwrap(),
        rotated
    );
}

fn pressure_surface(name: &str, offset: f64) -> PropertyClaim {
    let mut property = curve(name, [0.8, 0.4]);
    property.value = PropertyValue::Surface2D {
        axes: ["T".into(), "p".into()],
        axis_dims: [Temperature::DIMS, Pressure::DIMS],
        samples: [[280.0, 1e5], [320.0, 1e5], [280.0, 3e5], [320.0, 3e5]]
            .into_iter()
            .map(|point| {
                (
                    point,
                    0.7 - 0.005 * (point[0] - 280.0) + 1e-6 * (point[1] - 1e5) - offset,
                )
            })
            .collect(),
        triangles: vec![[0, 1, 2], [1, 3, 2]],
        dims: Dims::NONE,
    };
    property.validity = property
        .validity
        .with_quantity("p", pressure_quantity(), 1e5, 3e5);
    property
}

#[test]
fn original_two_dimensional_sources_use_explicit_nominal_pressure_and_temperature() {
    let card = card(vec![
        pressure_surface(STATIC, 0.0),
        pressure_surface(KINETIC, 0.2),
    ]);
    let bound = bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap();
    let point = query(300.0)
        .with_quantity("p", pressure_quantity(), 1e5)
        .unwrap();
    let mapped = ThermalFrictionAxes {
        temperature: "T".into(),
        nominal_pressure: Some("p".into()),
    };
    let law = SourceThermalCoulomb::resolve(
        bound.clone(),
        point,
        requirement(STATIC),
        requirement(KINETIC),
        MaterialPropertySelection::SingleClaimOnly,
        mapped.clone(),
    )
    .unwrap();
    let mut state = law.initial_state();
    for (temperature, pressure, force) in [(300.0, 1e5, 4.0), (300.0, 3e5, 6.0), (320.0, 3e5, 5.0)]
    {
        let mut trial = input(temperature, [0.01, 0.0]);
        trial.nominal_pressure_pa = Some(pressure);
        let step = law.advance(&state, &trial).unwrap();
        close(step.traction_n[0], force);
        close(step.dissipated_j, force * 0.01);
        assert_eq!(
            step.state.friction_source().state_coordinate("p"),
            Some(pressure)
        );
        state = step.state;
    }
    close(state.constitutive_state().slip_m()[0], 0.03);
    let saved = state.clone();
    assert!(law.advance(&state, &input(300.0, [0.01, 0.0])).is_err());
    for pressure in [-1.0, 0.0, 4e5, f64::NAN] {
        let mut trial = input(300.0, [0.01, 0.0]);
        trial.nominal_pressure_pa = Some(pressure);
        assert!(law.advance(&state, &trial).is_err());
    }
    assert_eq!(state, saved);
    let untyped = query(300.0)
        .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), 1e5)
        .unwrap();
    assert!(
        SourceThermalCoulomb::resolve(
            bound,
            untyped,
            requirement(STATIC),
            requirement(KINETIC),
            MaterialPropertySelection::SingleClaimOnly,
            mapped
        )
        .is_err()
    );
}

#[test]
fn runtime_source_and_mechanical_refusals_leave_history_and_receipts_retryable() {
    let law = resolve_at(&source_card(), query(300.0)).unwrap();
    let state = law
        .advance(&law.initial_state(), &input(300.0, [0.01, 0.0]))
        .unwrap()
        .state;
    let saved = state.clone();
    let next = input(310.0, [-0.005, 0.0]);
    let expected = law.advance_sourced(&state, &next).unwrap();
    assert!(matches!(
        law.advance_sourced(&state, &input(350.0, [0.01, 0.0])),
        Err(SourceThermalCoulombError::Material(_))
    ));
    for temperature in [None, Some(0.0), Some(f64::NAN)] {
        let mut bad = next;
        bad.contact_temperature_k = temperature;
        assert!(law.advance_sourced(&state, &bad).is_err());
    }
    let mut bad = next;
    bad.nominal_pressure_pa = Some(1e5);
    assert!(law.advance_sourced(&state, &bad).is_err());
    bad = next;
    bad.dt_s = 0.0;
    assert!(matches!(
        law.advance_sourced(&state, &bad),
        Err(SourceThermalCoulombError::Tribo(_))
    ));
    assert_eq!(state, saved);
    assert_eq!(law.advance_sourced(&state, &next).unwrap(), expected);

    let crossing = resolve_at(
        &card(vec![curve(STATIC, [0.8, 0.1]), curve(KINETIC, [0.6, 0.2])]),
        query(280.0),
    )
    .unwrap();
    let initial = crossing.initial_state();
    assert!(matches!(
        crossing.advance_sourced(&initial, &input(320.0, [0.01, 0.0])),
        Err(SourceThermalCoulombError::Tribo(_))
    ));
    let retry = crossing
        .advance_sourced(&initial, &input(300.0, [0.01, 0.0]))
        .unwrap();
    close(retry.traction_n[0], 4.0);
    close(retry.state.constitutive_state().slip_m()[0], 0.01);
}

#[test]
fn checkpoints_bind_original_order_history_query_roles_selection_and_authority() {
    let card = card(vec![curve(STATIC, [0.6, 0.2]), curve(KINETIC, [0.6, 0.2])]);
    let point = query(300.0)
        .with_quantity("ambient", pressure_quantity(), 1e5)
        .unwrap();
    let law = resolve_at(&card, point.clone()).unwrap();
    let state = law
        .advance(&law.initial_state(), &input(320.0, [0.01, 0.0]))
        .unwrap()
        .state;
    assert_eq!(
        state.friction_source().state_coordinate("ambient"),
        Some(1e5)
    );
    let mut history = card.context().clone();
    history.history = "different specimen history".into();
    let other_history = InterfaceSystemCard::assemble(
        card.surface_a().clone(),
        card.surface_b().clone(),
        history,
        card.claims().clone(),
        vec![],
    )
    .unwrap();
    let reversed = InterfaceSystemCard::assemble(
        card.surface_b().clone(),
        card.surface_a().clone(),
        card.context().clone(),
        card.claims().clone(),
        vec![],
    )
    .unwrap();
    let bound = bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap();
    let pins = law
        .initial_source()
        .properties()
        .iter()
        .map(|property| {
            (
                property.requirement().name().into(),
                property.answer().receipt.selected,
            )
        })
        .collect();
    let variants = vec![
        resolve_at(&other_history, point.clone()).unwrap(),
        resolve_at(&reversed, point.clone()).unwrap(),
        resolve_at(&card, point.clone().with("ambient", 1.01e5).unwrap()).unwrap(),
        SourceThermalCoulomb::resolve(
            bound.clone(),
            point.clone(),
            requirement(KINETIC),
            requirement(STATIC),
            MaterialPropertySelection::SingleClaimOnly,
            axes(),
        )
        .unwrap(),
        SourceThermalCoulomb::resolve(
            bound,
            point.clone(),
            requirement(STATIC),
            requirement(KINETIC),
            MaterialPropertySelection::PinnedByProperty(pins),
            axes(),
        )
        .unwrap(),
        SourceThermalCoulomb::resolve(
            bind_dry_interface_system_card(&card, InputAuthority::Estimated).unwrap(),
            point,
            requirement(STATIC),
            requirement(KINETIC),
            MaterialPropertySelection::SingleClaimOnly,
            axes(),
        )
        .unwrap(),
    ];
    for other in variants {
        assert_ne!(law.identity(), other.identity());
        assert!(other.stored_energy_j(&state).is_nan());
        assert!(other.advance(&state, &input(300.0, [0.01, 0.0])).is_err());
    }
    assert!(law.advance(&state, &input(300.0, [0.01, 0.0])).is_ok());
}

#[test]
fn a_different_valid_claim_cannot_replace_the_original_source_at_high_temperature() {
    let mut properties = vec![];
    for name in [STATIC, KINETIC] {
        properties.push(curve(name, [0.6, 0.3]));
        let mut hot = curve(name, [0.6, 0.3]);
        hot.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(330.0, 0.2), (360.0, 0.1)],
            dims: Dims::NONE,
        };
        hot.validity = ValidityDomain::unconstrained().with("T", 330.0, 360.0);
        hot.provenance.source = "another explicit high-temperature source".into();
        properties.push(hot);
    }
    let card = card(properties);
    let pins = |high: bool| {
        MaterialPropertySelection::PinnedByProperty(
            [STATIC, KINETIC]
                .into_iter()
                .map(|name| {
                    (
                        name.into(),
                        card.claims_for(name)
                            .into_iter()
                            .find(|(_, claim)| {
                                claim.provenance.source.starts_with("another") == high
                            })
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
    let state = law
        .advance(&law.initial_state(), &input(310.0, [0.01, 0.0]))
        .unwrap()
        .state;
    assert!(matches!(
        law.advance_sourced(&state, &input(350.0, [0.01, 0.0])),
        Err(SourceThermalCoulombError::Material(_))
    ));
    assert!(law.advance(&state, &input(300.0, [0.01, 0.0])).is_ok());
}

#[test]
fn kinematic_and_area_bounded_sources_refuse_even_with_constant_coefficients() {
    for axis in ["slip_speed", "slip", "contact_area", "strain", "history"] {
        let mut properties = vec![];
        for name in [STATIC, KINETIC] {
            let mut bounded = curve(name, [0.5; 2]);
            bounded.value = PropertyValue::Scalar {
                value: 0.5,
                dims: Dims::NONE,
            };
            bounded.interpolation = InterpolationPolicy::ConstantWithinValidity;
            bounded.validity = bounded.validity.with(axis, 0.0, 1.0);
            properties.push(bounded);
        }
        assert!(
            resolve_at(&card(properties), query(300.0).with(axis, 0.0).unwrap()).is_err(),
            "{axis}"
        );
    }
}

#[test]
fn missing_wrong_unit_and_unphysical_coefficients_refuse_while_explicit_zero_is_physical() {
    assert!(resolve_at(&card(vec![]), query(300.0)).is_err());
    let mut wrong = curve(KINETIC, [0.4; 2]);
    wrong.key = PropertyKey::new(KINETIC, Pressure::DIMS);
    wrong.value = PropertyValue::Scalar {
        value: 0.4,
        dims: Pressure::DIMS,
    };
    wrong.interpolation = InterpolationPolicy::ConstantWithinValidity;
    assert!(resolve_at(&card(vec![curve(STATIC, [0.6; 2]), wrong]), query(300.0)).is_err());
    assert!(
        resolve_at(
            &card(vec![curve(STATIC, [0.2; 2]), curve(KINETIC, [0.6; 2])]),
            query(300.0)
        )
        .is_err()
    );
    let negative = card(vec![curve(STATIC, [0.2; 2]), curve(KINETIC, [-0.1; 2])]);
    assert!(
        SourceThermalCoulomb::resolve(
            bind_dry_interface_system_card(&negative, InputAuthority::SyntheticFixture).unwrap(),
            query(300.0),
            requirement(STATIC),
            ScalarPropertyRequirement::try_new(KINETIC, Dims::NONE, ScalarAdmissibility::Finite)
                .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            axes()
        )
        .is_err()
    );
    let zero = card(vec![curve("single-explicit-coefficient", [0.0; 2])]);
    let shared = requirement("single-explicit-coefficient");
    let law = SourceThermalCoulomb::resolve(
        bind_dry_interface_system_card(&zero, InputAuthority::SyntheticFixture).unwrap(),
        query(300.0),
        shared.clone(),
        shared,
        MaterialPropertySelection::SingleClaimOnly,
        axes(),
    )
    .unwrap();
    let step = law
        .advance(&law.initial_state(), &input(320.0, [0.01, -0.02]))
        .unwrap();
    assert_eq!(step.traction_n, [0.0; 2]);
    assert_eq!(step.slip_increment_m, [0.01, -0.02]);
    assert_eq!(step.dissipated_j, 0.0);
    assert_eq!(step.stored_change_j, 0.0);
}
