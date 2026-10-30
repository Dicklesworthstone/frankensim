//! Source-temperature feedback on the actual two-sided modal bow and heat owner.
//! All coefficient curves are manufactured conservation fixtures, not rosin data.
use super::*;
use fs_contact::interface_binding::bind_dry_interface_system_card;
use fs_contact::interface_binding::thermal_friction::{
    FrictionQueryAxis, FrictionQueryAxisKind, SourceThermalCoulomb, ThermalFrictionAxes,
};
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use std::cell::Cell;

fn temperature_binding(material: &MaterialCard) -> DryFrictionBinding {
    let mut claims = ClaimSet::new();
    for (name, values) in [("static", [1.0, 0.4]), ("kinetic", [0.5, 0.2])] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, Dims::NONE),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, values[0]), (340.0, values[1])], dims: Dims::NONE },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 340.0)
                .with("p", 0.001, 1e6).with("v", 0.0, 10.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: Provenance {
                source: "manufactured temperature-only dry coefficient; not measured data".into(),
                license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
    }
    let surface = SurfaceSpec { material: material.id().clone(), texture_frame: "fixture tangent".into() };
    let card = InterfaceSystemCard::assemble(surface.clone(), surface, SystemContext {
        medium: "dry".into(), third_body: None, environment: "manufactured".into(),
        history: "unchanged original surfaces".into(),
    }, claims, vec![]).unwrap();
    let semantic = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = QueryPoint::new().with_quantity("T", semantic(QuantityKind::AbsoluteTemperature), 300.0).unwrap()
        .with_quantity("p", semantic(QuantityKind::Pressure), 0.1 / AREA).unwrap()
        .with_quantity("v", QuantitySpec::dimensional(fs_qty::Velocity::DIMS), 0.0).unwrap();
    let requirement = |name| ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative).unwrap();
    let source = SourceThermalCoulomb::resolve_with_query_axes(
        bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(), point,
        requirement("static"), requirement("kinetic"), MaterialPropertySelection::SingleClaimOnly,
        ThermalFrictionAxes { temperature: "T".into(), nominal_pressure: Some("p".into()) },
        vec![FrictionQueryAxis { name: "v".into(), kind: FrictionQueryAxisKind::SlidingSpeed }],
    ).unwrap();
    DryFrictionBinding::new_temperature_source(source,
        DryFrictionApplicability::new(Range::new(280.0, 340.0).unwrap(),
            Range::new(0.001, 1e6).unwrap(), Range::new(0.0, 10.0).unwrap()).unwrap(), AREA).unwrap()
}

fn temperature_bank(cx: &Cx<'_>, gate: &CancelGate, card: &MaterialCard) -> MultiBowStringBank {
    reacting_surface(joint(cx, gate, 7800.0), gate)
        .with_bow_interface(0, temperature_binding(card)).unwrap()
}

#[test]
fn original_temperature_claims_and_finite_capacity_close_same_traction_heat_and_pressure() {
    with_cx(|cx, gate| {
        let mut observations = Vec::new();
        for cp in [20.0, CP] {
            let (card, caloric) = coating_law_with_capacity(false, None, None, cp);
            let bundle = transport(&card);
            let context = point(300.0);
            let depths = [DEPTH; 2];
            let routes = [route(&card, caloric.curve(), &bundle, &context, depths)];
            let mut bank = temperature_bank(cx, gate, &card);
            let source_before = bank.bow_interface(0).unwrap().clone();
            let mut heat = dry_network(caloric.curve(), [300.0; 2], depths);
            let forces = vec![0.0; bank.system().mode_count()];
            let frame = bank.step_with_dry_contact_thermal_motion(&mut heat, &routes, &[],
                solve_controls(), &[0.0; 6], &forces, &[drive(0.05, 0.1), drive(0.0, 0.0)],
                None, None, None, gate, other).unwrap();
            let contact = &frame.dry_contacts[0];
            let temperature = contact.interface_temperature_k.unwrap();
            let traction = frame.mechanics.ports[0].traction_n;
            let slip = frame.mechanics.ports[0].slip_distance_m;
            let dt = bank.system().sample_period_s();
            assert!(temperature > 300.0 && slip > 0.0);
            let coefficient = 0.4 - 0.005 * (temperature - 300.0);
            assert!((traction - 0.1 * coefficient).abs() < 2e-9);
            assert!(traction < 0.04, "a previous-temperature law cannot pass the coupled force equation");
            assert!(contact.interface_energy_residual_j.abs() <= solve_controls().total_energy_tolerance_j);
            assert!((total_heat(&heat) - traction * slip).abs() < 2e-12);
            assert!(frame.total_energy_residual_j.abs() < 2e-12);
            // Independently eliminate the two finite caloric cells in one BE
            // step. The interface has no capacity: G_eff is each body's flux.
            let g = 0.2 * AREA / (0.5 * DEPTH);
            let capacity = DENSITY * AREA * DEPTH * cp;
            let effective = g * capacity / (capacity + dt * g);
            let expected = 300.0 + traction * slip / (dt * 2.0 * effective);
            assert!((temperature - expected).abs() < 1e-6);
            let accepted = bank.bow_interface(0).unwrap();
            assert_eq!(&bank.system().port(0).unwrap().law, accepted.card().friction_law());
            let receipt = accepted.resolved_temperature_source().unwrap();
            assert_eq!(receipt.state_coordinate("T"), Some(temperature));
            assert_eq!(receipt.state_coordinate("p"), Some(0.1 / AREA));
            assert_eq!(receipt.state_coordinate("v"), Some((slip / dt).abs()));
            for property in receipt.properties() {
                assert_eq!(property.answer().receipt.selected,
                    source_before.resolved_temperature_source().unwrap()
                        .property(property.requirement().name()).unwrap().answer().receipt.selected);
            }
            observations.push((temperature, traction, frame.mechanics.observer_pressure_pa));
        }
        assert!(observations[0].0 > observations[1].0);
        assert!(observations[0].1 < observations[1].1);
        assert_ne!(observations[0].2, observations[1].2);
    });
}

#[test]
fn source_feedback_refusal_restores_law_motion_and_cells_then_release_preserves_source() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(false, false);
        let bundle = transport(&card);
        let context = point(300.0);
        let depths = [DEPTH; 2];
        let routes = [route(&card, caloric.curve(), &bundle, &context, depths)];
        let make = || temperature_bank(cx, gate, &card).with_bow_station_motion(motion_limits(), gate).unwrap();
        let mut bank = make();
        let mut reference = make();
        let mut heat = dry_network(caloric.curve(), [300.0; 2], depths);
        let mut reference_heat = dry_network(caloric.curve(), [300.0; 2], depths);
        let forces = vec![0.0; bank.system().mode_count()];
        let drives = [drive(0.05, 0.1), drive(0.0, 0.0)];
        let saved = (states(&bank), bank.bow_interface(0).unwrap().clone(),
            bank.system().port(0).unwrap().clone(), heat.cells().clone(),
            bank.bow_station_fractions().unwrap().to_vec());
        assert!(bank.step(&forces, &drives).is_err(), "evolving law cannot silently take an isothermal step");
        let calls = Cell::new(0);
        let result = bank.step_with_dry_contact_thermal_motion(&mut heat, &routes, &[],
            solve_controls(), &[0.0; 6], &forces, &drives, Some(&[0.202, 0.3]), None, None, gate,
            |endpoint| {
                calls.set(calls.get() + 1);
                if calls.get() == 4 { Err("manufactured late transport refusal") }
                else { Ok(other(endpoint).unwrap()) }
            });
        assert!(result.is_err());
        assert_eq!(calls.get(), 4, "refusal must occur after real constitutive trials");
        assert_eq!(states(&bank), saved.0);
        assert_eq!(bank.bow_interface(0).unwrap(), &saved.1);
        assert_eq!(bank.system().port(0).unwrap(), &saved.2);
        assert_eq!(heat.cells(), &saved.3);
        assert_eq!(bank.bow_station_fractions().unwrap(), saved.4);
        assert_eq!(heat.time_s(), 0.0);
        assert_eq!(bank.system().samples_rendered(), 0);
        let first = bank.step_with_dry_contact_thermal_motion(&mut heat, &routes, &[],
            solve_controls(), &[0.0; 6], &forces, &drives, Some(&[0.202, 0.3]), None, None, gate, other).unwrap();
        let expected = reference.step_with_dry_contact_thermal_motion(&mut reference_heat, &routes, &[],
            solve_controls(), &[0.0; 6], &forces, &drives, Some(&[0.202, 0.3]), None, None, gate, other).unwrap();
        assert_eq!(first.mechanics, expected.mechanics);
        assert_eq!(heat.cells(), reference_heat.cells());
        let source = bank.bow_interface(0).unwrap().clone();
        let cells = heat.cells().clone();
        let guess = guesses(&heat);
        let released = bank.step_with_dry_contact_thermal_motion(&mut heat, &routes, &[],
            solve_controls(), &guess, &forces, &[drive(0.05, 0.0), drive(0.0, 0.0)],
            None, None, None, gate, other).unwrap();
        assert_eq!(released.dry_contacts[0].interface_temperature_k, None);
        assert!(released.dry_contacts[0].cell_heat_j.is_empty());
        assert_eq!(released.routed_dry_heat_j, 0.0);
        assert_eq!(bank.bow_interface(0).unwrap(), &source);
        assert_eq!(heat.cells(), &cells);
    });
}
