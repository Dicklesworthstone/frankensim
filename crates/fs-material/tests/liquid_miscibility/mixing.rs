//! G0/G1/G3 finite isothermal mixing, heat duty and source ownership.

use super::*;
use fs_material::liquid::solution::ResolvedBinaryLiquidEquilibrium;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn pair(card: &MaterialCard, a: f64, b: f64) -> BinaryLiquidSolutionLaw<'_> {
    BinaryLiquidSolutionLaw::try_new(
        BinaryRegularSolutionV1::new(
            [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
            a,
            b,
            Temperature::new(290.0),
            Temperature::new(600.0),
        )
        .unwrap(),
        [MolarMass::new(0.018), MolarMass::new(0.10)],
        card,
        reference(card, CompositionBasis::MoleFraction),
        "T",
        "x_a",
    )
    .unwrap()
}

fn state(
    law: &BinaryLiquidSolutionLaw<'_>,
    temperature: f64,
    components: [f64; 2],
) -> ResolvedBinaryLiquidEquilibrium {
    law.equilibrate(
        Temperature::new(temperature),
        components,
        BinaryCoexistenceOptionsV1::default(),
        || false,
    )
    .unwrap()
}

#[test]
fn pure_component_mixing_has_analytic_heat_and_entropy_signs() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    for (a, b) in [(0.0, 0.0), (0.0, 900.0), (0.0, -900.0), (1.5, -300.0)] {
        let law = pair(&card, a, b);
        let left = state(&law, 500.0, [1.0, 0.0]);
        let right = state(&law, 500.0, [0.0, 1.0]);
        let mixing = law
            .mix_isothermal(&left, &right, BinaryCoexistenceOptionsV1::default(), || {
                false
            })
            .unwrap();
        assert_eq!(mixing.equilibrium().phases().len(), 1);
        assert_eq!(mixing.equilibrium().total_component_moles(), [1.0, 1.0]);
        near(mixing.heat_into_liquid_j(), 0.5 * R * b);
        let entropy = 2.0 * R * (2.0_f64.ln() - a / 4.0);
        near(mixing.liquid_entropy_change_j_per_k(), entropy);
        near(mixing.gibbs_change_j(), 0.5 * R * b - 500.0 * entropy);
        near(mixing.entropy_generated_j_per_k(), entropy - 0.001 * R * b);
        assert!(mixing.gibbs_change_j() < 0.0);
        assert!(mixing.entropy_generated_j_per_k() > 0.0);
        near(
            mixing.liquid_entropy_change_j_per_k() + mixing.reservoir_entropy_change_j_per_k(),
            mixing.entropy_generated_j_per_k(),
        );
        near(mixing.equilibrium().total_mass_kg(), 0.118);
        assert_eq!(
            mixing.equilibrium().phases()[0].liquid().density_kg_m3(),
            1100.0
        );
    }
}

#[test]
fn heat_duty_uses_coexisting_phases_and_retains_species_mass_and_volume() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let law = law(&card);
    let left = state(&law, 300.0, [2.0, 0.0]);
    let right = state(&law, 300.0, [0.0, 3.0]);
    let mixing = law
        .mix_isothermal(&left, &right, BinaryCoexistenceOptionsV1::default(), || {
            false
        })
        .unwrap();
    let result = mixing.equilibrium();
    assert_eq!(result.phases().len(), 2);
    let enthalpy: f64 = result
        .phases()
        .iter()
        .map(|phase| {
            let n: f64 = phase.component_moles().iter().sum();
            let [x, y] = phase.mole_fractions();
            n * R * 900.0 * x * y
        })
        .sum();
    near(mixing.heat_into_liquid_j(), enthalpy);
    let homogeneous_heat = 5.0 * R * 900.0 * 0.4 * 0.6;
    assert!(mixing.heat_into_liquid_j() > 0.0);
    assert!(mixing.heat_into_liquid_j() < 0.5 * homogeneous_heat);
    for i in 0..2 {
        near(
            result
                .phases()
                .iter()
                .map(|phase| phase.component_moles()[i])
                .sum(),
            left.total_component_moles()[i] + right.total_component_moles()[i],
        );
    }
    near(
        result.total_mass_kg(),
        left.total_mass_kg() + right.total_mass_kg(),
    );
    near(
        result.total_volume_m3(),
        result.phases().iter().map(|phase| phase.volume_m3()).sum(),
    );
    // Composition-resolved density admits contraction/expansion; summing the
    // original pure-fluid volumes would erase that source behavior.
    assert!(
        (result.total_volume_m3() - left.total_volume_m3() - right.total_volume_m3()).abs() > 1e-6
    );
    for phase in result.phases() {
        assert_eq!(
            phase.thermal_transport_state().unwrap().temperature_k(),
            300.0
        );
    }
    let reverse = law
        .mix_isothermal(&right, &left, BinaryCoexistenceOptionsV1::default(), || {
            false
        })
        .unwrap();
    assert_eq!(reverse, mixing);
    assert_eq!(reverse.into_equilibrium(), *result);
}

#[test]
fn identical_intensive_states_have_no_fictitious_heat_or_entropy() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let law = law(&card);
    let left = state(&law, 500.0, [2.0, 3.0]);
    let right = state(&law, 500.0, [4.0, 6.0]);
    let mixing = law
        .mix_isothermal(&left, &right, BinaryCoexistenceOptionsV1::default(), || {
            false
        })
        .unwrap();
    assert_eq!(mixing.equilibrium().total_component_moles(), [6.0, 9.0]);
    assert_eq!(mixing.heat_into_liquid_j(), 0.0);
    assert_eq!(mixing.liquid_entropy_change_j_per_k(), 0.0);
    assert_eq!(mixing.gibbs_change_j(), 0.0);
    assert_eq!(mixing.entropy_generated_j_per_k(), 0.0);
}

#[test]
fn same_names_cannot_replace_the_original_thermodynamic_or_material_owner() {
    let original = card(CompositionBasis::MoleFraction, 1.0);
    let changed = card(CompositionBasis::MoleFraction, 0.8);
    let law = law(&original);
    let left = state(&law, 500.0, [1.0, 0.0]);
    let variants = [
        pair(&original, 0.0, 0.0),
        pair(&changed, 0.0, 900.0),
        BinaryLiquidSolutionLaw::try_new(
            BinaryRegularSolutionV1::new(
                [SpeciesId::new("B").unwrap(), SpeciesId::new("A").unwrap()],
                0.0,
                900.0,
                Temperature::new(290.0),
                Temperature::new(600.0),
            )
            .unwrap(),
            [MolarMass::new(0.018), MolarMass::new(0.10)],
            &original,
            reference(&original, CompositionBasis::MoleFraction),
            "T",
            "x_a",
        )
        .unwrap(),
        BinaryLiquidSolutionLaw::try_new(
            model(),
            [MolarMass::new(0.10), MolarMass::new(0.018)],
            &original,
            reference(&original, CompositionBasis::MoleFraction),
            "T",
            "x_a",
        )
        .unwrap(),
    ];
    for other in variants {
        let right = state(&other, 500.0, [0.0, 1.0]);
        assert!(matches!(
            law.mix_isothermal(&left, &right, BinaryCoexistenceOptionsV1::default(), || {
                false
            }),
            Err(BinaryLiquidError::Invalid(_))
        ));
    }
    let right = state(&law, 450.0, [0.0, 1.0]);
    assert!(matches!(
        law.mix_isothermal(&left, &right, BinaryCoexistenceOptionsV1::default(), || {
            false
        }),
        Err(BinaryLiquidError::Invalid(_))
    ));
    // Equality is immutable value identity, not allocation address identity.
    let identical_law = pair(&original, 0.0, 900.0);
    let right = state(&identical_law, 500.0, [0.0, 1.0]);
    assert!(
        law.mix_isothermal(&left, &right, BinaryCoexistenceOptionsV1::default(), || {
            false
        })
        .is_ok()
    );
}

#[test]
fn unsupported_final_source_cannot_consume_either_original_inventory() {
    let ordinary = card(CompositionBasis::MoleFraction, 1.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in ordinary.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == DENSITY_PROPERTY {
            let PropertyValue::Curve { knots, .. } = &mut claim.value else {
                panic!("fixture density must be a curve")
            };
            // Valid endpoints, but the actual final mixture is not an admitted
            // positive-density liquid. Reference x=0.4 is still admissible.
            knots[1].1 = 0.0;
        }
        claims.insert_claim(claim).unwrap();
    }
    let missing = MaterialCard::assemble(ordinary.id().clone(), claims, vec![]).unwrap();
    let law = pair(&missing, 0.0, 0.0);
    let left = state(&law, 500.0, [1.0, 0.0]);
    let right = state(&law, 500.0, [0.0, 1.0]);
    let originals = [left.clone(), right.clone()];
    assert!(matches!(
        law.mix_isothermal(&left, &right, BinaryCoexistenceOptionsV1::default(), || {
            false
        }),
        Err(BinaryLiquidError::Source(_))
    ));
    assert_eq!([left, right], originals);
}

#[test]
fn transport_heat_capacity_does_not_count_the_mixing_enthalpy_twice() {
    let original = card(CompositionBasis::MoleFraction, 1.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == SPECIFIC_HEAT_CAPACITY_PROPERTY {
            let PropertyValue::Curve { knots, .. } = &mut claim.value else {
                panic!("fixture heat capacity must be a curve")
            };
            for (_, value) in knots {
                *value *= 10.0;
            }
        }
        claims.insert_claim(claim).unwrap();
    }
    let changed = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let mut heat = [0.0; 2];
    for (index, card) in [&original, &changed].into_iter().enumerate() {
        let law = law(card);
        let left = state(&law, 500.0, [1.0, 0.0]);
        let right = state(&law, 500.0, [0.0, 1.0]);
        let mixing = law
            .mix_isothermal(&left, &right, BinaryCoexistenceOptionsV1::default(), || {
                false
            })
            .unwrap();
        heat[index] = mixing.heat_into_liquid_j();
        let cp = mixing.equilibrium().phases()[0]
            .transport_state()
            .unwrap()
            .specific_heat_j_kg_k();
        near(cp, if index == 0 { 3500.0 } else { 35000.0 });
    }
    assert_eq!(heat[0], heat[1]);
}

#[test]
fn cancellation_and_inventory_overflow_publish_no_transaction() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let law = law(&card);
    let left = state(&law, 300.0, [2.0, 0.0]);
    let right = state(&law, 300.0, [0.0, 3.0]);
    let originals = [left.clone(), right.clone()];
    let mut count = 0;
    law.mix_isothermal(&left, &right, BinaryCoexistenceOptionsV1::default(), || {
        count += 1;
        false
    })
    .unwrap();
    for stop in [0, 1, count - 1] {
        let mut polls = 0;
        let result =
            law.mix_isothermal(&left, &right, BinaryCoexistenceOptionsV1::default(), || {
                polls += 1;
                polls > stop
            });
        assert!(matches!(
            result,
            Err(BinaryLiquidError::Cancelled)
                | Err(BinaryLiquidError::Thermodynamics(
                    fs_thermochem::solution::BinarySolutionErrorV1::Cancelled
                ))
        ));
        assert_eq!([left.clone(), right.clone()], originals);
    }
    let enormous = state(&law, 500.0, [f64::MAX * 0.75, 0.0]);
    assert!(matches!(
        law.mix_isothermal(
            &enormous,
            &enormous,
            BinaryCoexistenceOptionsV1::default(),
            || false
        ),
        Err(BinaryLiquidError::Arithmetic(_))
    ));
    let large = state(&law, 500.0, [1e20, 0.0]);
    let small = state(&law, 500.0, [1.0, 0.0]);
    assert!(matches!(
        law.mix_isothermal(
            &large,
            &small,
            BinaryCoexistenceOptionsV1::default(),
            || false
        ),
        Err(BinaryLiquidError::Arithmetic(_))
    ));
}
