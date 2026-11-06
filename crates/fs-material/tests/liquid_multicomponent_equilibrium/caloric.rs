//! Original pure-liquid receipts feed the accepted finite phase response.

use super::*;
use fs_material::liquid::solution::multicomponent::caloric::{
    MulticomponentLiquidCaloricError, MulticomponentLiquidCaloricLaw,
};
use fs_material::thermal::transport::isobaric::{
    IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery,
};
use fs_thermochem::solution::multicomponent::caloric::MulticomponentCaloricErrorV1;
use fs_thermochem::solution::multicomponent::equilibrium::caloric::MulticomponentCaloricResponseControlV1;

const MOLAR_CP_AT_400: [f64; 3] = [40.0, 70.0, 110.0];

fn cp_source(index: usize, pressure: f64, datum: f64) -> IsobaricHeatCapacityLaw {
    let key = PropertyKey::new("authored_pure_cp", SPECIFIC_HEAT_CAPACITY_DIMS);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: [300.0, 450.0, 600.0]
                    .into_iter()
                    .zip(
                        [0.6, 1.2, 2.1].map(|scale| scale * MOLAR_CP_AT_400[index] / MASSES[index]),
                    )
                    .collect(),
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    300.0,
                    600.0,
                )
                .with_quantity(
                    "source_pressure",
                    quantity(QuantityKind::Pressure),
                    pressure,
                    pressure,
                ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: format!("authored pure-liquid caloric source for component {index}"),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("authored pure component {index}"),
            phase: "liquid".into(),
            process: "synthetic source-binding fixture".into(),
            revision: 0,
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
            lower: &point(300.0),
            upper: &point(600.0),
            selection: fs_matdb::ClaimSelection::Policy(fs_matdb::SelectionPolicy::SingleClaimOnly),
        },
        datum,
    )
    .unwrap()
}

fn sources() -> Vec<(SpeciesId, IsobaricHeatCapacityLaw)> {
    [2, 0, 1]
        .into_iter()
        .map(|i| (species(i), cp_source(i, 100_000.0, 300.0)))
        .collect()
}

fn enthalpic_model() -> MulticomponentRegularSolutionV1 {
    MulticomponentRegularSolutionV1::new(
        (0..3).map(species).collect(),
        vec![0.0; 9],
        vec![
            0.0, 1600.0, 1600.0, 1600.0, 0.0, 1600.0, 1600.0, 1600.0, 0.0,
        ],
        Temperature::new(250.0),
        Temperature::new(650.0),
    )
    .unwrap()
}

#[test]
fn canonical_sources_preserve_receipts_and_original_molar_mass_conversion() {
    let card = card(1.0);
    let liquid = law(&card, model(0.0));
    let mut original = sources();
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, original.clone(), "p").unwrap();
    original.sort_by(|left, right| left.0.cmp(&right.0));
    assert!(core::ptr::eq(caloric.liquid_law(), &liquid));
    assert_eq!(
        caloric.temperature_domain(),
        [Temperature::new(300.0), Temperature::new(600.0)]
    );
    for (i, (source, pure)) in caloric
        .sources()
        .iter()
        .zip(caloric.model().heat_capacities())
        .enumerate()
    {
        assert_eq!(source, &original[i].1);
        assert_eq!(pure.species(), &species(i));
        for (t, answer) in &source.source_line().knots {
            near(
                pure.heat_capacity_at(Temperature::new(*t)).unwrap().value(),
                answer.evidence.value.value * MASSES[i],
            );
        }
    }
    let amounts = [2.0, 3.0, 5.0];
    let equilibrium = liquid
        .resolve_equilibrium(
            candidate(liquid.thermodynamics(), 400.0, amounts, false),
            || false,
        )
        .unwrap();
    let response = caloric
        .isobaric_caloric_response(
            &equilibrium,
            MulticomponentCaloricResponseControlV1::default(),
            || false,
        )
        .unwrap();
    let expected: f64 = amounts
        .into_iter()
        .zip(MOLAR_CP_AT_400)
        .map(|(n, cp)| n * cp)
        .sum();
    near(
        response.caloric().isobaric_heat_capacity_j_per_k(),
        expected,
    );
    near(
        response.specific_heat_cp_j_kg_k(),
        expected / equilibrium.total_mass_kg(),
    );
    assert_eq!(
        response.caloric().redistribution_heat_capacity_j_per_k(),
        0.0
    );
    assert!(core::ptr::eq(response.equilibrium(), &equilibrium));
    assert_eq!(response.sources(), caloric.sources());
    assert!(core::ptr::eq(
        response.caloric().pure_liquid_heat_capacities(),
        caloric.model().heat_capacities()
    ));
    assert_eq!(response.pressure_axis(), "p");
    assert!(core::ptr::eq(
        response.caloric().candidate(),
        equilibrium.equilibrium_candidate()
    ));
    assert_ne!(
        response.specific_heat_cp_j_kg_k(),
        equilibrium.phases()[0].liquid().specific_heat_j_kg_k()
    );
    let homogeneous = caloric
        .homogeneous(Temperature::new(400.0), &amounts, || false)
        .unwrap();
    near(
        homogeneous.caloric_state().isobaric_heat_capacity_j_per_k(),
        response.caloric().frozen_heat_capacity_j_per_k(),
    );
    assert_eq!(homogeneous.sources(), response.sources());
}

#[test]
fn three_phase_response_tracks_source_enthalpy_and_keeps_mass_specific_scaling() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let amounts = [2.0, 3.0, 5.0];
    let equilibrium = liquid
        .resolve_equilibrium(
            candidate(liquid.thermodynamics(), 400.0, amounts, true),
            || false,
        )
        .unwrap();
    let response = caloric
        .isobaric_caloric_response(
            &equilibrium,
            MulticomponentCaloricResponseControlV1::default(),
            || false,
        )
        .unwrap();
    assert_eq!(response.equilibrium().phases().len(), 3);
    assert!(response.caloric().redistribution_heat_capacity_j_per_k() > 0.0);
    let frozen: f64 = amounts
        .into_iter()
        .zip(MOLAR_CP_AT_400)
        .map(|(n, cp)| n * cp)
        .sum();
    near(response.caloric().frozen_heat_capacity_j_per_k(), frozen);

    // Independent neighboring equilibrium solves and original source enthalpies
    // check the actual fixed-inventory derivative, including phase redistribution.
    let enthalpy = |temperature: f64| {
        let selected = candidate(liquid.thermodynamics(), temperature, amounts, true);
        let sensible: f64 = amounts
            .iter()
            .zip(MASSES)
            .zip(caloric.sources())
            .map(|((n, mass), source)| {
                n * mass * source.specific_enthalpy_j_kg(temperature).unwrap()
            })
            .sum();
        sensible
            + selected
                .phases()
                .iter()
                .map(|phase| {
                    phase.total_moles() * phase.thermodynamics().mixing_enthalpy_j_per_mol()
                })
                .sum::<f64>()
    };
    let slope = (enthalpy(400.001) - enthalpy(399.999)) / 0.002;
    assert!(
        (slope - response.caloric().isobaric_heat_capacity_j_per_k()).abs() < 5e-6 * slope.abs()
    );
    let scaled = liquid
        .resolve_equilibrium(
            candidate(
                liquid.thermodynamics(),
                400.0,
                amounts.map(|n| n * 17.0),
                true,
            ),
            || false,
        )
        .unwrap();
    let scaled_response = caloric
        .isobaric_caloric_response(
            &scaled,
            MulticomponentCaloricResponseControlV1::default(),
            || false,
        )
        .unwrap();
    near(
        scaled_response.specific_heat_cp_j_kg_k(),
        response.specific_heat_cp_j_kg_k(),
    );
    near(
        scaled_response.caloric().isobaric_heat_capacity_j_per_k(),
        response.caloric().isobaric_heat_capacity_j_per_k() * 17.0,
    );
    for (phase, original) in response
        .equilibrium()
        .phases()
        .iter()
        .zip(equilibrium.phases())
    {
        assert!(core::ptr::eq(phase, original));
    }
}

#[test]
fn bindings_refuse_missing_duplicate_species_pressure_datum_and_complete_foreign_owner() {
    let card = card(1.0);
    let liquid = law(&card, model(0.0));
    let mut missing = sources();
    missing.pop();
    let mut duplicate = sources();
    duplicate[0].0 = duplicate[1].0.clone();
    let mut foreign_pressure = sources();
    foreign_pressure[0].1 = cp_source(2, 200_000.0, 300.0);
    let mut foreign_datum = sources();
    foreign_datum[0].1 = cp_source(2, 100_000.0, 310.0);
    for supplied in [missing, duplicate, foreign_pressure] {
        assert!(matches!(
            MulticomponentLiquidCaloricLaw::try_new(&liquid, supplied, "p"),
            Err(MulticomponentLiquidCaloricError::Invalid(_))
        ));
    }
    assert!(matches!(
        MulticomponentLiquidCaloricLaw::try_new(&liquid, foreign_datum, "p"),
        Err(MulticomponentLiquidCaloricError::Caloric(
            MulticomponentCaloricErrorV1::Invalid(_)
        ))
    ));
    for axis in ["missing", "T", "x_a"] {
        assert!(matches!(
            MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), axis),
            Err(MulticomponentLiquidCaloricError::Invalid(_))
        ));
    }
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let mut altered_components = liquid.components().to_vec();
    altered_components[0].molar_mass = MolarMass::new(0.020);
    let altered = MulticomponentLiquidSolutionLaw::try_new(
        liquid.thermodynamics().clone(),
        altered_components,
        &card,
        liquid.reference_state().clone(),
        "T",
    )
    .unwrap();
    let foreign = altered
        .resolve_equilibrium(
            candidate(altered.thermodynamics(), 400.0, [2.0, 3.0, 5.0], false),
            || false,
        )
        .unwrap();
    assert!(matches!(
        caloric.isobaric_caloric_response(
            &foreign,
            MulticomponentCaloricResponseControlV1::default(),
            || false,
        ),
        Err(MulticomponentLiquidCaloricError::Invalid(_))
    ));
    // An identically reconstructed immutable owner remains admissible.
    let identical = law(&card, model(0.0));
    let reconstructed = identical
        .resolve_equilibrium(
            candidate(identical.thermodynamics(), 400.0, [2.0, 3.0, 5.0], false),
            || false,
        )
        .unwrap();
    caloric
        .isobaric_caloric_response(
            &reconstructed,
            MulticomponentCaloricResponseControlV1::default(),
            || false,
        )
        .unwrap();
}

#[test]
fn cancellation_never_changes_original_material_or_source_states() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let original_sources = sources();
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, original_sources, "p").unwrap();
    let equilibrium = liquid
        .resolve_equilibrium(
            candidate(liquid.thermodynamics(), 400.0, [2.0, 3.0, 5.0], true),
            || false,
        )
        .unwrap();
    let original = equilibrium.clone();
    let receipts = caloric.sources().to_vec();
    let mut polls = 0;
    caloric
        .isobaric_caloric_response(
            &equilibrium,
            MulticomponentCaloricResponseControlV1::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut current = 0;
        assert!(
            caloric
                .isobaric_caloric_response(
                    &equilibrium,
                    MulticomponentCaloricResponseControlV1::default(),
                    || {
                        current += 1;
                        current == stop
                    },
                )
                .is_err()
        );
        assert_eq!(equilibrium, original);
        assert_eq!(caloric.sources(), receipts);
    }
}
