//! G1/G3 source-qualified caloric mixing, transport refresh and atomic refusal.

#[path = "caloric/storage.rs"]
mod storage;

use super::*;
use fs_material::liquid::solution::caloric::{BinaryLiquidCaloricError, BinaryLiquidCaloricLaw};
use fs_material::thermal::transport::isobaric::{
    IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery,
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::solution::caloric::BinaryCaloricSolveV1;

fn liquid(card: &MaterialCard, interaction_k: f64) -> BinaryLiquidSolutionLaw<'_> {
    let point = QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 500.0)
        .unwrap()
        .with_quantity("x_a", composition(CompositionBasis::MoleFraction), 0.4)
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
        .unwrap();
    BinaryLiquidSolutionLaw::try_new(
        BinaryRegularSolutionV1::new(
            [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
            0.0,
            interaction_k,
            Temperature::new(290.0),
            Temperature::new(600.0),
        )
        .unwrap(),
        [MolarMass::new(0.018), MolarMass::new(0.10)],
        card,
        resolve_liquid_state_from_kinematic(
            card,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
        "T",
        "x_a",
    )
    .unwrap()
}

fn cp_source(
    species: &str,
    molar_mass: f64,
    molar_cp: [f64; 3],
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
                knots: [290.0, 400.0, 600.0]
                    .into_iter()
                    .zip(molar_cp.map(|cp| cp / molar_mass))
                    .collect(),
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    290.0,
                    600.0,
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
                source: format!("synthetic explicitly pure {species} calorics revision {revision}"),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: species.into(),
            phase: "liquid".into(),
            process: "synthetic caloric fixture".into(),
            revision,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = |t| {
        QueryPoint::new()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                t,
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
            lower: &point(290.0),
            upper: &point(600.0),
            selection: fs_matdb::ClaimSelection::Policy(fs_matdb::SelectionPolicy::SingleClaimOnly),
        },
        datum,
    )
    .unwrap()
}

fn sources() -> [(SpeciesId, IsobaricHeatCapacityLaw); 2] {
    [
        (
            SpeciesId::new("A").unwrap(),
            cp_source("A", 0.018, [100.0; 3], 100_000.0, 300.0, 0),
        ),
        (
            SpeciesId::new("B").unwrap(),
            cp_source("B", 0.10, [100.0; 3], 100_000.0, 300.0, 0),
        ),
    ]
}
fn controls() -> BinaryCaloricSolveV1 {
    BinaryCaloricSolveV1 {
        temperature_bracket: [Temperature::new(290.0), Temperature::new(600.0)],
        temperature_tolerance_k: 1e-10,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 128,
        coexistence: BinaryCoexistenceOptionsV1::default(),
    }
}

#[test]
fn pure_source_calorics_close_mixing_without_replacing_mixture_transport() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid(&card, 100.0);
    let original = sources();
    let thermal = BinaryLiquidCaloricLaw::try_new(&liquid, original.clone(), "pressure").unwrap();
    assert_eq!(
        thermal.temperature_domain(),
        [Temperature::new(290.0), Temperature::new(600.0)]
    );
    let left = thermal
        .equilibrate(
            Temperature::new(300.0),
            [1.0, 0.0],
            controls().coexistence,
            || false,
        )
        .unwrap();
    let right = thermal
        .equilibrate(
            Temperature::new(300.0),
            [0.0, 1.0],
            controls().coexistence,
            || false,
        )
        .unwrap();
    let result = thermal
        .mix(&left, &right, 0.0, controls(), || false)
        .unwrap();
    let expected_t = 300.0 - R / 4.0;
    assert!((result.equilibrium().temperature().value() - expected_t).abs() < 1e-9);
    assert!(result.balance().energy_residual_j().abs() <= controls().energy_tolerance_j);
    assert_eq!(result.equilibrium().total_component_moles(), [1.0, 1.0]);
    assert_eq!(result.equilibrium().phases().len(), 1);
    let phase = &result.equilibrium().phases()[0];
    assert_eq!(phase.liquid().density_kg_m3(), 1100.0);
    assert_eq!(
        phase.transport_state().unwrap().specific_heat_j_kg_k(),
        3500.0
    );
    assert_eq!(
        phase.thermal_transport_state().unwrap().temperature_k(),
        result.equilibrium().temperature().value()
    );
    assert_eq!(phase.liquid().material().state_coordinate("x_a"), Some(0.5));
    for (i, (_, original)) in original.iter().enumerate() {
        assert_eq!(&thermal.sources()[i], original);
        assert_eq!(&result.final_state().sources()[i], original);
        assert_eq!(
            result.final_state().sources()[i].source_line(),
            original.source_line()
        );
    }
    for property in phase.liquid().material().properties() {
        assert_eq!(
            property.answer().receipt.selected,
            liquid
                .reference_state()
                .material()
                .property(property.requirement().name())
                .unwrap()
                .answer()
                .receipt
                .selected
        );
    }
    assert_eq!(left.caloric_state().temperature(), Temperature::new(300.0));
    assert_eq!(right.caloric_state().temperature(), Temperature::new(300.0));
}

#[test]
fn unequal_temperature_parcels_integrate_all_original_cp_segments_and_reheat() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid(&card, 0.0);
    let mut sources = sources();
    // Cp_A = T - 200 J/(mol K); Cp_B = 100 J/(mol K).
    sources[0].1 = cp_source("A", 0.018, [90.0, 200.0, 400.0], 100_000.0, 300.0, 0);
    let thermal = BinaryLiquidCaloricLaw::try_new(&liquid, sources, "pressure").unwrap();
    let left = thermal
        .equilibrate(
            Temperature::new(300.0),
            [1.0, 0.0],
            controls().coexistence,
            || false,
        )
        .unwrap();
    let right = thermal
        .equilibrate(
            Temperature::new(500.0),
            [0.0, 1.0],
            controls().coexistence,
            || false,
        )
        .unwrap();
    let mixed = thermal
        .mix(&left, &right, 5000.0, controls(), || false)
        .unwrap();
    // A absorbs 15 kJ and B releases 10 kJ to reach 400 K.
    assert!((mixed.equilibrium().temperature().value() - 400.0).abs() < 1e-9);
    assert!((mixed.balance().enthalpy_change_j() - 5000.0).abs() <= controls().energy_tolerance_j);
    let reheated = thermal
        .add_heat(mixed.final_state(), 16_250.0, controls(), || false)
        .unwrap();
    assert!((reheated.equilibrium().temperature().value() - 450.0).abs() < 1e-9);
    assert_eq!(reheated.equilibrium().total_component_moles(), [1.0, 1.0]);
    let cooled = thermal
        .add_heat(reheated.final_state(), -16_250.0, controls(), || false)
        .unwrap();
    assert!((cooled.equilibrium().temperature().value() - 400.0).abs() < 1e-9);
    assert_eq!(
        cooled.final_state().sources(),
        mixed.final_state().sources()
    );
}

#[test]
fn source_identity_species_pressure_and_datum_cannot_be_substituted() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid(&card, 0.0);
    let original = sources();
    let thermal = BinaryLiquidCaloricLaw::try_new(&liquid, original.clone(), "pressure").unwrap();
    let left = thermal
        .equilibrate(
            Temperature::new(400.0),
            [1.0, 0.0],
            controls().coexistence,
            || false,
        )
        .unwrap();
    let mut changed = original.clone();
    changed[1].1 = cp_source("B", 0.10, [100.0; 3], 100_000.0, 300.0, 1);
    let other = BinaryLiquidCaloricLaw::try_new(&liquid, changed, "pressure").unwrap();
    let right = other
        .equilibrate(
            Temperature::new(400.0),
            [0.0, 1.0],
            controls().coexistence,
            || false,
        )
        .unwrap();
    // Identical numeric curves do not authorize different source provenance.
    assert!(matches!(
        thermal.mix(&left, &right, 0.0, controls(), || false),
        Err(BinaryLiquidCaloricError::Invalid(_))
    ));
    assert!(matches!(
        thermal.add_heat(&right, 0.0, controls(), || false),
        Err(BinaryLiquidCaloricError::Invalid(_))
    ));
    let rebuilt = BinaryLiquidCaloricLaw::try_new(&liquid, original.clone(), "pressure").unwrap();
    let right = rebuilt
        .equilibrate(
            Temperature::new(400.0),
            [0.0, 1.0],
            controls().coexistence,
            || false,
        )
        .unwrap();
    assert!(
        thermal
            .mix(&left, &right, 0.0, controls(), || false)
            .is_ok()
    );
    let changed_card = super::card(CompositionBasis::MoleFraction, 0.8);
    let changed_liquid = self::liquid(&changed_card, 0.0);
    let changed_material =
        BinaryLiquidCaloricLaw::try_new(&changed_liquid, original.clone(), "pressure").unwrap();
    let other_state = changed_material
        .equilibrate(
            Temperature::new(400.0),
            [0.0, 1.0],
            controls().coexistence,
            || false,
        )
        .unwrap();
    assert!(matches!(
        thermal.mix(&left, &other_state, 0.0, controls(), || false),
        Err(BinaryLiquidCaloricError::Invalid(_))
    ));
    let mut reversed = original.clone();
    reversed.swap(0, 1);
    assert!(BinaryLiquidCaloricLaw::try_new(&liquid, reversed, "pressure").is_err());
    let mut wrong_pressure = original.clone();
    wrong_pressure[1].1 = cp_source("B", 0.10, [100.0; 3], 200_000.0, 300.0, 0);
    assert!(BinaryLiquidCaloricLaw::try_new(&liquid, wrong_pressure, "pressure").is_err());
    let mut wrong_datum = original.clone();
    wrong_datum[1].1 = cp_source("B", 0.10, [100.0; 3], 100_000.0, 350.0, 0);
    assert!(BinaryLiquidCaloricLaw::try_new(&liquid, wrong_datum, "pressure").is_err());
    let untyped = law(&card);
    assert!(BinaryLiquidCaloricLaw::try_new(&untyped, original, "pressure").is_err());
}

#[test]
fn accepted_caloric_phase_selection_is_retained_exactly_on_material_refresh() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid(&card, 900.0);
    let thermal = BinaryLiquidCaloricLaw::try_new(&liquid, sources(), "pressure").unwrap();
    for (temperature, amounts) in [
        (300.0, [2.0, 3.0]),
        (450.0, [1.0, 1.0]),
        (500.0, [1.5e-16, 1.0]),
        (500.0, [1.0, 1.5e-16]),
    ] {
        let state = thermal
            .equilibrate(
                Temperature::new(temperature),
                amounts,
                controls().coexistence,
                || false,
            )
            .unwrap();
        assert_eq!(
            state.caloric_state().equilibrium().unwrap(),
            state.equilibrium().thermodynamic_equilibrium(),
        );
        assert_eq!(state.equilibrium().total_component_moles(), amounts);
    }
    assert!(matches!(
        thermal.equilibrate(
            Temperature::new(500.0),
            [1.0, 1e-20],
            controls().coexistence,
            || false
        ),
        Err(BinaryLiquidCaloricError::Material(
            BinaryLiquidError::Arithmetic(_)
        )),
    ));
}

#[test]
fn thermal_endpoint_source_failure_and_last_cancellation_preserve_inputs() {
    let restricted = card(CompositionBasis::MoleFraction, 0.8);
    let liquid = liquid(&restricted, 900.0);
    let thermal = BinaryLiquidCaloricLaw::try_new(&liquid, sources(), "pressure").unwrap();
    let initial = thermal
        .equilibrate(
            Temperature::new(500.0),
            [2.0, 3.0],
            controls().coexistence,
            || false,
        )
        .unwrap();
    let before = initial.clone();
    let cold = liquid
        .thermodynamics()
        .equilibrate(Temperature::new(300.0), 0.4, controls().coexistence, || {
            false
        })
        .unwrap();
    let heat = 5.0
        * (100.0 * (300.0 - 500.0) + cold.mixing_enthalpy_j_per_mol()
            - initial
                .equilibrium()
                .thermodynamic_equilibrium()
                .mixing_enthalpy_j_per_mol());
    assert!(matches!(
        thermal.add_heat(&initial, heat, controls(), || false),
        Err(BinaryLiquidCaloricError::Material(
            BinaryLiquidError::Source(_)
        ))
    ));
    assert_eq!(initial, before);
    assert_eq!(
        liquid.reference_state().material().state_coordinate("T"),
        Some(500.0)
    );

    let ordinary = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = self::liquid(&ordinary, 900.0);
    let thermal = BinaryLiquidCaloricLaw::try_new(&liquid, sources(), "pressure").unwrap();
    let initial = thermal
        .equilibrate(
            Temperature::new(500.0),
            [2.0, 3.0],
            controls().coexistence,
            || false,
        )
        .unwrap();
    let before = initial.clone();
    let mut polls = 0;
    let completed = thermal
        .add_heat(&initial, -50_000.0, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(completed.equilibrium().temperature().value() < 500.0);
    let mut seen = 0;
    assert_eq!(
        thermal.add_heat(&initial, -50_000.0, controls(), || {
            seen += 1;
            seen == polls
        }),
        Err(BinaryLiquidCaloricError::Cancelled)
    );
    assert_eq!(seen, polls);
    assert_eq!(initial, before);
}
