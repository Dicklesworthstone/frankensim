//! G0/G1/G4 original caloric sources, adiabatic demixing and atomic source refusal.

use super::*;
use fs_flux::miscibility::material::thermal::BinaryMiscibilityMaterialThermalError;
use fs_flux::miscibility::thermal::BinaryMiscibilityThermalSolve;
use fs_material::liquid::solution::caloric::BinaryLiquidCaloricLaw;
use fs_material::thermal::transport::isobaric::{
    IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery,
};

fn pure_cp(species: &str, molar_mass: f64, molar_cp: [f64; 3]) -> IsobaricHeatCapacityLaw {
    let key = PropertyKey::new("qualified_pure_liquid_cp", SPECIFIC_HEAT_CAPACITY_DIMS);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: [300.0, 450.0, 600.0]
                    .into_iter()
                    .zip(molar_cp.map(|cp| cp / molar_mass))
                    .collect(),
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    quantity(QuantityKind::AbsoluteTemperature),
                    300.0,
                    600.0,
                )
                .with_quantity(
                    "pressure",
                    quantity(QuantityKind::Pressure),
                    100_000.0,
                    100_000.0,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: format!("synthetic qualified pure {species} calorics"),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: species.into(),
            phase: "liquid".into(),
            process: "pure caloric fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = |t| {
        QueryPoint::new()
            .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), t)
            .unwrap()
            .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
            .unwrap()
    };
    IsobaricHeatCapacityLaw::try_new(
        &card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &key,
            temperature_axis: "T",
            pressure_axis: "pressure",
            lower: &point(300.0),
            upper: &point(600.0),
            selection: fs_matdb::ClaimSelection::Policy(fs_matdb::SelectionPolicy::SingleClaimOnly),
        },
        300.0,
    )
    .unwrap()
}

fn caloric<'liquid, 'card>(
    liquid: &'liquid BinaryLiquidSolutionLaw<'card>,
) -> BinaryLiquidCaloricLaw<'liquid, 'card> {
    BinaryLiquidCaloricLaw::try_new(
        liquid,
        [
            (
                SpeciesId::new("A").unwrap(),
                pure_cp("A", 0.018, [40.0, 60.0, 80.0]),
            ),
            (
                SpeciesId::new("B").unwrap(),
                pure_cp("B", 0.10, [80.0, 120.0, 160.0]),
            ),
        ],
        "pressure",
    )
    .unwrap()
}

fn controls() -> BinaryMiscibilityThermalSolve {
    BinaryMiscibilityThermalSolve {
        temperature_bracket: [Temperature::new(300.0), Temperature::new(600.0)],
        temperature_tolerance_k: 1e-10,
        energy_tolerance_j: 1e-8,
        entropy_tolerance_j_per_k: 1e-10,
        maximum_iterations: 128,
    }
}

#[test]
fn g1_adiabatic_demixing_uses_original_pure_cp_and_refreshes_final_materials() {
    let card = temperature_card(1.0, 600.0);
    let liquid = liquid(&card);
    let caloric = caloric(&liquid);
    let source = BinaryLiquidSquareGradientLaw::try_new(&liquid, model(), 1e-12).unwrap();
    let net = network(model(), 1e-5);
    let initial = net
        .state_from_mole_fractions(vec![[0.4, 0.6], [0.6, 0.4]], || false)
        .unwrap();
    let field = net
        .resolve_liquid_field(&initial, &source, 2, || false)
        .unwrap();
    let thermal = field.bind_caloric(&caloric, || false).unwrap();
    assert!(std::ptr::eq(
        field.cells().as_ptr(),
        thermal.cells().as_ptr()
    ));
    assert!(std::ptr::eq(thermal.caloric_law(), &caloric));
    assert_eq!(thermal.state(), &initial);
    let mut polls = 0;
    let candidate = thermal
        .adiabatic_candidate(0.1, BinaryMiscibilitySolve::default(), controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let end = candidate.field();
    let step = candidate.step();
    let temperature = end.state().model().temperature().value();
    let rise = temperature - 300.0;
    assert!(
        rise > 0.0,
        "positive excess enthalpy is released by demixing"
    );
    assert_eq!(end.state().elapsed_s(), 0.1);
    assert_eq!(end.maximum_resolved_cells(), 2);
    assert!(std::ptr::eq(end.caloric_law(), &caloric));
    assert!(std::ptr::eq(end.source().liquid_law(), &liquid));
    assert_eq!(end.source().gradient_model(), end.state().model());
    assert_eq!(end.caloric_law().model(), caloric.model());
    for i in 0..2 {
        assert_eq!(
            end.caloric_law().sources()[i].identity(),
            caloric.sources()[i].identity()
        );
        near(
            end.state().component_moles().iter().map(|n| n[i]).sum(),
            1.0,
        );
    }
    // Exactly one mol of each species: integrated pure Cp is 120*dT+0.2*dT^2.
    // The independent mixture transport Cp would give a different heat account.
    near(
        step.sensible_enthalpy_change_j,
        120.0 * rise + 0.2 * rise * rise,
    );
    assert!(step.energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert!(step.entropy_change_j_per_k >= -controls().entropy_tolerance_j_per_k);
    assert!(step.mixing_enthalpy_change_j < 0.0);
    for cell in end.cells() {
        assert!(cell.thermodynamics().thermodynamic_factor() < 0.0);
        near(cell.phase().volume_m3(), 1e-4);
        let actual = cell.phase().transport_state().unwrap();
        near(
            actual.kinematic_viscosity_m2_s(),
            6e-6 - 5e-6 * rise / 300.0,
        );
        near(actual.specific_heat_j_kg_k(), 3500.0 + 500.0 * rise / 300.0);
        let material = cell.phase().liquid().material();
        assert_eq!(material.state_coordinate("T"), Some(temperature));
        assert_eq!(material.state_coordinate("pressure"), Some(100_000.0));
        assert_eq!(
            material.card_identity(),
            liquid.reference_state().material().card_identity()
        );
    }
    for stop in [1, polls] {
        let mut calls = 0;
        assert!(
            thermal
                .adiabatic_candidate(0.1, BinaryMiscibilitySolve::default(), controls(), || {
                    calls += 1;
                    calls == stop
                })
                .is_err()
        );
        assert_eq!(thermal.state(), &initial);
        assert_eq!(thermal.cells(), field.cells());
    }
    let (accepted, _) = candidate.into_parts();
    let next = accepted
        .adiabatic_candidate(0.1, BinaryMiscibilitySolve::default(), controls(), || false)
        .unwrap();
    assert_eq!(next.field().state().elapsed_s(), 0.2);
    assert!(std::ptr::eq(next.field().caloric_law(), &caloric));
    assert_eq!(field.state(), &initial);
}

#[test]
fn g0_caloric_binding_requires_the_exact_original_liquid_source() {
    let card = temperature_card(1.0, 600.0);
    let liquid = liquid(&card);
    let different_liquid = self::liquid(&card);
    let other_caloric = caloric(&different_liquid);
    let source = BinaryLiquidSquareGradientLaw::try_new(&liquid, model(), 1e-12).unwrap();
    let net = network(model(), 1e-5);
    let initial = net
        .state_from_mole_fractions(vec![[0.4, 0.6]; 2], || false)
        .unwrap();
    let field = net
        .resolve_liquid_field(&initial, &source, 2, || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        field.bind_caloric(&other_caloric, || {
            polls += 1;
            false
        }),
        Err(BinaryMiscibilityMaterialThermalError::Invalid(_))
    ));
    assert_eq!(
        polls, 1,
        "different source owner refuses before numeric preparation"
    );
    assert_eq!(field.state(), &initial);
    assert!(std::ptr::eq(field.source(), &source));
}

#[test]
fn g4_final_adiabatic_density_or_domain_refusal_preserves_the_original_field() {
    for (density_scale, maximum_temperature) in [(0.9, 600.0), (1.0, 300.000001)] {
        let card = temperature_card(density_scale, maximum_temperature);
        let liquid = liquid(&card);
        let caloric = caloric(&liquid);
        let source = BinaryLiquidSquareGradientLaw::try_new(&liquid, model(), 1e-12).unwrap();
        let net = network(model(), 1e-5);
        let initial = net
            .state_from_mole_fractions(vec![[0.4, 0.6], [0.6, 0.4]], || false)
            .unwrap();
        let field = net
            .resolve_liquid_field(&initial, &source, 2, || false)
            .unwrap();
        let thermal = field.bind_caloric(&caloric, || false).unwrap();
        let numerical = net
            .bind_caloric_field(&initial, caloric.model(), || false)
            .unwrap()
            .adiabatic_candidate(0.1, BinaryMiscibilitySolve::default(), controls(), || false)
            .unwrap();
        assert!(numerical.field().state().model().temperature().value() > 300.000001);
        assert!(matches!(
            thermal
                .adiabatic_candidate(0.1, BinaryMiscibilitySolve::default(), controls(), || false),
            Err(BinaryMiscibilityMaterialThermalError::Material(
                BinaryMiscibilityMaterialError::Source { cell: 0, .. }
            ))
        ));
        assert_eq!(thermal.state(), &initial);
        assert_eq!(thermal.cells(), field.cells());
        assert!(std::ptr::eq(thermal.source(), &source));
        assert!(std::ptr::eq(thermal.caloric_law(), &caloric));
    }
}
