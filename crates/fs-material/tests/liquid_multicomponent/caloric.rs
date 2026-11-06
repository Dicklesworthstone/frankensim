//! G1/G3 original-source homogeneous thermal mixing and transactional refresh.

use super::*;
use fs_material::liquid::solution::multicomponent::caloric::{
    MulticomponentLiquidCaloricError, MulticomponentLiquidCaloricLaw,
};
use fs_material::thermal::transport::isobaric::{
    IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery,
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::solution::multicomponent::caloric::MulticomponentCaloricSolveV1;

#[path = "caloric/storage.rs"]
mod storage;

fn liquid(card: &MaterialCard) -> MulticomponentLiquidSolutionLaw<'_> {
    MulticomponentLiquidSolutionLaw::try_new(
        MulticomponentRegularSolutionV1::new(
            (0..3).map(species).collect(),
            vec![0.0, 0.3, -0.1, 0.3, 0.0, 0.6, -0.1, 0.6, 0.0],
            vec![0.0, 120.0, -60.0, 120.0, 0.0, 80.0, -60.0, 80.0, 0.0],
            Temperature::new(250.0),
            Temperature::new(650.0),
        )
        .unwrap(),
        components([2, 0, 1]),
        card,
        reference(card, CompositionBasis::MoleFraction),
        "T",
    )
    .unwrap()
}
fn cp_source(
    index: usize,
    capacities: [f64; 3],
    pressure: f64,
    datum: f64,
    revision: u32,
) -> IsobaricHeatCapacityLaw {
    let key = PropertyKey::new("authored_pure_liquid_cp", SPECIFIC_HEAT_CAPACITY_DIMS);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: [300.0, 400.0, 650.0]
                    .into_iter()
                    .zip(capacities.map(|cp| cp / MASSES[index]))
                    .collect(),
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    300.0,
                    650.0,
                )
                .with_quantity(
                    "source_pressure",
                    quantity(QuantityKind::Pressure),
                    pressure,
                    pressure,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: format!("synthetic pure component {index} Cp source {revision}"),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("component {index}"),
            phase: "liquid".into(),
            process: "authored pure-liquid calorics".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = |temperature| {
        QueryPoint::new()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                temperature,
            )
            .unwrap()
            .with_quantity(
                "source_pressure",
                quantity(QuantityKind::Pressure),
                pressure,
            )
            .unwrap()
    };
    IsobaricHeatCapacityLaw::try_new(
        &card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &key,
            temperature_axis: "temperature",
            pressure_axis: "source_pressure",
            lower: &point(300.0),
            upper: &point(650.0),
            selection: fs_matdb::ClaimSelection::Policy(fs_matdb::SelectionPolicy::SingleClaimOnly),
        },
        datum,
    )
    .unwrap()
}
fn sources() -> Vec<(SpeciesId, IsobaricHeatCapacityLaw)> {
    [2, 0, 1]
        .into_iter()
        .map(|i| {
            (
                species(i),
                cp_source(
                    i,
                    [[80.0; 3], [120.0; 3], [200.0; 3]][i],
                    100_000.0,
                    300.0,
                    0,
                ),
            )
        })
        .collect()
}
fn controls() -> MulticomponentCaloricSolveV1 {
    MulticomponentCaloricSolveV1 {
        temperature_bracket: [Temperature::new(300.0), Temperature::new(650.0)],
        temperature_tolerance_k: 1e-9,
        energy_tolerance_j: 1e-7,
        entropy_tolerance_j_per_k: 1e-10,
        maximum_iterations: 128,
    }
}

#[test]
fn g1_ternary_heat_of_mixing_preserves_pure_sources_and_refreshes_actual_material() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid(&card);
    let offered_sources = sources();
    let law =
        MulticomponentLiquidCaloricLaw::try_new(&liquid, offered_sources.clone(), "p").unwrap();
    for (i, source) in law.sources().iter().enumerate() {
        assert_eq!(
            source,
            &offered_sources
                .iter()
                .find(|(id, _)| *id == species(i))
                .unwrap()
                .1
        );
        assert_eq!(law.model().heat_capacities()[i].species(), &species(i));
    }
    let left = law
        .homogeneous(Temperature::new(320.0), &[1.0, 0.0, 0.0], || false)
        .unwrap();
    let right = law
        .homogeneous(Temperature::new(500.0), &[0.0, 2.0, 1.0], || false)
        .unwrap();
    let result = law.mix(&left, &right, 0.0, controls(), || false).unwrap();
    let expected_t = (80.0 * 320.0 + 440.0 * 500.0 - (95.0 / 3.0) * R) / 520.0;
    near(
        result.final_state().caloric_state().temperature().value(),
        expected_t,
    );
    assert!(result.balance().energy_residual_j().abs() < 1e-7);
    assert!(result.entropy_change_j_per_k() > 0.0);
    assert_eq!(result.liquid().component_moles(), &[1.0, 2.0, 1.0]);
    near(
        result.liquid().mass_kg(),
        MASSES[0] + 2.0 * MASSES[1] + MASSES[2],
    );
    let actual_source = liquid
        .homogeneous(
            result.final_state().caloric_state().temperature(),
            &[1.0, 2.0, 1.0],
            || false,
        )
        .unwrap();
    assert_eq!(result.liquid(), &actual_source);
    assert_ne!(
        result.liquid().transport_state().unwrap(),
        left.liquid().transport_state().unwrap()
    );
    near(
        result
            .liquid()
            .thermal_transport_state()
            .unwrap()
            .temperature_k(),
        expected_t,
    );
    assert_eq!(
        result.liquid().liquid().material().state_coordinate("p"),
        Some(100_000.0)
    );
    assert_eq!(
        result.liquid().liquid().material().card_identity(),
        card.content_hash()
    );
    assert_eq!(result.final_state().sources(), law.sources());
    for property in result.liquid().liquid().material().properties() {
        let original = liquid
            .reference_state()
            .material()
            .property(property.requirement().name())
            .unwrap();
        assert_eq!(
            property.answer().receipt.selected,
            original.answer().receipt.selected
        );
    }
}

#[test]
fn g1_variable_pure_cp_heats_and_cools_without_substituting_transport_cp() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid(&card);
    let mut sources = sources();
    sources[1].1 = cp_source(0, [80.0, 100.0, 150.0], 100_000.0, 300.0, 0);
    let law = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources, "p").unwrap();
    let initial = law
        .homogeneous(Temperature::new(350.0), &[1.0, 2.0, 1.0], || false)
        .unwrap();
    // Cp_A = 20 + 0.2 T; B and C contribute 440 J/K. The 400 K source
    // knot is crossed, preserving the original mass-to-molar conversion.
    let heat = 107.0 * 170.0 + 440.0 * 170.0;
    let result = law.add_heat(&initial, heat, controls(), || false).unwrap();
    near(
        result.final_state().caloric_state().temperature().value(),
        520.0,
    );
    near(
        result
            .final_state()
            .caloric_state()
            .isobaric_heat_capacity_j_per_k(),
        124.0 + 440.0,
    );
    near(
        law.enthalpy_change_j(&initial, result.final_state())
            .unwrap(),
        heat,
    );
    let expected_s = 460.0 * (520.0_f64 / 350.0).ln() + 0.2 * 170.0;
    near(result.entropy_change_j_per_k(), expected_s);
    near(
        law.entropy_change_j_per_k(&initial, result.final_state())
            .unwrap(),
        expected_s,
    );
    let cooled = law
        .add_heat(result.final_state(), -heat, controls(), || false)
        .unwrap();
    near(
        cooled.final_state().caloric_state().temperature().value(),
        350.0,
    );
    near(cooled.entropy_change_j_per_k(), -expected_s);
    assert_eq!(
        cooled.liquid().component_moles(),
        initial.liquid().component_moles()
    );
}

#[test]
fn g3_original_source_binding_and_final_source_support_refuse_atomically() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid(&card);
    let law = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let initial = law
        .homogeneous(Temperature::new(350.0), &[1.0, 2.0, 1.0], || false)
        .unwrap();
    let original = initial.clone();
    let mut wrong = sources();
    wrong[0].0 = species(0);
    assert!(MulticomponentLiquidCaloricLaw::try_new(&liquid, wrong, "p").is_err());
    let mut wrong = sources();
    wrong[0].1 = cp_source(2, [200.0; 3], 200_000.0, 300.0, 0);
    assert!(MulticomponentLiquidCaloricLaw::try_new(&liquid, wrong, "p").is_err());
    let mut wrong = sources();
    wrong[0].1 = cp_source(2, [200.0; 3], 100_000.0, 350.0, 0);
    assert!(MulticomponentLiquidCaloricLaw::try_new(&liquid, wrong, "p").is_err());
    assert!(MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "T").is_err());
    let mut different = sources();
    different[0].1 = cp_source(2, [200.0; 3], 100_000.0, 300.0, 1);
    let other = MulticomponentLiquidCaloricLaw::try_new(&liquid, different, "p").unwrap();
    assert!(other.add_heat(&initial, 0.0, controls(), || false).is_err());
    // Pure sources admit 625 K, but the actual mixture source stops at 600 K.
    let heat = 520.0 * (625.0 - 350.0);
    let numeric = law
        .model()
        .add_heat(initial.caloric_state(), heat, controls(), || false)
        .unwrap();
    near(numeric.final_state().temperature().value(), 625.0);
    assert!(matches!(
        law.add_heat(&initial, heat, controls(), || false),
        Err(MulticomponentLiquidCaloricError::Material(_))
    ));
    let mut calls = 0;
    law.add_heat(&initial, 1000.0, controls(), || {
        calls += 1;
        false
    })
    .unwrap();
    for cutoff in [0, calls / 2, calls - 1] {
        let mut current = 0;
        assert!(
            law.add_heat(&initial, 1000.0, controls(), || {
                let stop = current >= cutoff;
                current += 1;
                stop
            })
            .is_err()
        );
    }
    assert_eq!(initial, original);
}
