use super::*;
use crate::solution::caloric::{BinaryCaloricSolveV1, BinaryLiquidCaloricModelV1};
use crate::solution::{BinaryCoexistenceOptionsV1, BinaryRegularSolutionV1};
use crate::{MolarThermalQuantityV1, SpeciesId};

fn ids(n: usize) -> Vec<SpeciesId> {
    (0..n)
        .map(|i| SpeciesId::new(format!("S{i:03}")).unwrap())
        .collect()
}
fn matrix(n: usize, pairs: &[(usize, usize, f64)]) -> Vec<f64> {
    let mut values = vec![0.0; n * n];
    for &(i, j, value) in pairs {
        values[i * n + j] = value;
        values[j * n + i] = value;
    }
    values
}
fn potential(
    n: usize,
    a: &[(usize, usize, f64)],
    b: &[(usize, usize, f64)],
) -> MulticomponentRegularSolutionV1 {
    MulticomponentRegularSolutionV1::new(
        ids(n),
        matrix(n, a),
        matrix(n, b),
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap()
}
fn source(species: SpeciesId, datum: f64, knots: &[(f64, f64)]) -> PureLiquidHeatCapacityV1 {
    PureLiquidHeatCapacityV1::new(
        species,
        Temperature::new(datum),
        knots
            .iter()
            .map(|&(t, cp)| (Temperature::new(t), MolarThermalQuantityV1::new(cp)))
            .collect(),
    )
    .unwrap()
}
fn model(
    a: &[(usize, usize, f64)],
    b: &[(usize, usize, f64)],
    capacities: &[f64],
    datum: f64,
) -> MulticomponentLiquidCaloricModelV1 {
    MulticomponentLiquidCaloricModelV1::new(
        potential(capacities.len(), a, b),
        ids(capacities.len())
            .into_iter()
            .zip(capacities)
            .map(|(id, &cp)| source(id, datum, &[(100.0, cp), (1e15, cp)]))
            .collect(),
    )
    .unwrap()
}
fn controls() -> MulticomponentCaloricSolveV1 {
    MulticomponentCaloricSolveV1 {
        temperature_bracket: [Temperature::new(200.0), Temperature::new(800.0)],
        temperature_tolerance_k: 1e-9,
        energy_tolerance_j: 1e-7,
        entropy_tolerance_j_per_k: 1e-10,
        maximum_iterations: 128,
    }
}
fn parcel(
    model: &MulticomponentLiquidCaloricModelV1,
    t: f64,
    amounts: &[f64],
) -> MulticomponentCaloricStateV1 {
    model
        .homogeneous(Temperature::new(t), amounts, || false)
        .unwrap()
}
fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_ternary_unequal_temperature_mix_closes_all_pair_heat_and_entropy() {
    let a = [(0, 1, 0.3), (0, 2, -0.1), (1, 2, 0.6)];
    let b = [(0, 1, 120.0), (0, 2, -60.0), (1, 2, 80.0)];
    let expected_t = (80.0 * 300.0 + 440.0 * 500.0 - (95.0 / 3.0) * R) / 520.0;
    let x = [0.25_f64, 0.5, 0.25];
    let before_x = [0.0_f64, 2.0 / 3.0, 1.0 / 3.0];
    let ideal = |fractions: &[f64]| -> f64 {
        fractions
            .iter()
            .filter(|&&x| x != 0.0)
            .map(|x| x * x.ln())
            .sum()
    };
    let excess = |fractions: &[f64]| -> f64 {
        a.iter()
            .map(|&(i, j, value)| value * fractions[i] * fractions[j])
            .sum()
    };
    let expected_s = 80.0 * (expected_t / 300.0).ln() + 440.0 * (expected_t / 500.0).ln()
        - R * (4.0 * ideal(&x) - 3.0 * ideal(&before_x))
        - R * (4.0 * excess(&x) - 3.0 * excess(&before_x));
    for datum in [100.0, 1e15] {
        // Huge sensible-reference enthalpies must not erase the mixing heat.
        let model = model(&a, &b, &[80.0, 120.0, 200.0], datum);
        let left = parcel(&model, 300.0, &[1.0, 0.0, 0.0]);
        let right = parcel(&model, 500.0, &[0.0, 2.0, 1.0]);
        let result = model.mix(&left, &right, 0.0, controls(), || false).unwrap();
        near(
            result.final_state().temperature().value(),
            expected_t,
            0.0,
            1e-9,
        );
        near(result.entropy_change_j_per_k(), expected_s, 2e-12, 1e-10);
        near(result.energy_residual_j(), 0.0, 0.0, 1e-7);
        assert_eq!(result.final_state().component_moles(), &[1.0, 2.0, 1.0]);
        assert_eq!(result.final_state().isobaric_heat_capacity_j_per_k(), 520.0);
        assert!(result.entropy_change_j_per_k() > 0.0);
        if datum == 100.0 {
            near(
                result.final_state().enthalpy_j(),
                left.enthalpy_j() + right.enthalpy_j(),
                2e-14,
                1e-7,
            );
        }
        let reversed = model.mix(&right, &left, 0.0, controls(), || false).unwrap();
        near(
            reversed.final_state().temperature().value(),
            expected_t,
            0.0,
            1e-9,
        );
        near(reversed.entropy_change_j_per_k(), expected_s, 2e-12, 1e-10);
    }
}

#[test]
fn g1_variable_cp_heating_and_cooling_cross_original_source_knots() {
    let species = ids(3);
    let model = MulticomponentLiquidCaloricModelV1::new(
        potential(3, &[], &[]),
        vec![
            source(species[2].clone(), 100.0, &[(100.0, 50.0), (700.0, 50.0)]),
            source(
                species[0].clone(),
                100.0,
                &[(100.0, 20.0), (300.0, 40.0), (700.0, 60.0)],
            ),
            source(
                species[1].clone(),
                100.0,
                &[(100.0, 100.0), (400.0, 40.0), (700.0, 70.0)],
            ),
        ],
    )
    .unwrap();
    let initial = parcel(&model, 250.0, &[2.0, 1.0, 0.0]);
    let controls = MulticomponentCaloricSolveV1 {
        temperature_bracket: [Temperature::new(200.0), Temperature::new(700.0)],
        ..controls()
    };
    let result = model
        .add_heat(&initial, 42_250.0, controls, || false)
        .unwrap();
    near(result.final_state().temperature().value(), 550.0, 0.0, 1e-9);
    near(
        result.final_state().isobaric_heat_capacity_j_per_k(),
        160.0,
        0.0,
        1e-9,
    );
    let expected_s = 2.0
        * (0.1 * 50.0
            + 10.0 * (300.0_f64 / 250.0).ln()
            + 0.05 * 250.0
            + 25.0 * (550.0_f64 / 300.0).ln())
        - 0.2 * 150.0
        + 120.0 * (400.0_f64 / 250.0).ln()
        + 0.1 * 150.0;
    near(result.entropy_change_j_per_k(), expected_s, 2e-12, 1e-9);
    near(
        model
            .enthalpy_change_j(&initial, result.final_state())
            .unwrap(),
        42_250.0,
        0.0,
        1e-7,
    );
    near(
        model
            .entropy_change_j_per_k(&initial, result.final_state())
            .unwrap(),
        expected_s,
        2e-12,
        1e-9,
    );
    let cooled = model
        .add_heat(result.final_state(), -42_250.0, controls, || false)
        .unwrap();
    near(cooled.final_state().temperature().value(), 250.0, 0.0, 1e-9);
    near(cooled.entropy_change_j_per_k(), -expected_s, 2e-12, 1e-9);
    assert_eq!(
        cooled.final_state().component_moles(),
        initial.component_moles()
    );
}

#[test]
fn g1_maximum_species_keeps_trace_and_absent_components_without_entropy_floors() {
    let model = model(&[], &[], &[1.0; MAX_REGULAR_SOLUTION_SPECIES_V1], 100.0);
    let mut left = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
    let mut right = left;
    left[0] = 1.0;
    right[1] = 1e-100;
    let left = parcel(&model, 300.0, &left);
    let right = parcel(&model, 300.0, &right);
    let result = model.mix(&left, &right, 0.0, controls(), || false).unwrap();
    assert_eq!(result.final_state().temperature().value(), 300.0);
    assert_eq!(&result.final_state().component_moles()[..2], &[1.0, 1e-100]);
    assert!(
        result.final_state().component_moles()[2..]
            .iter()
            .all(|&n| n == 0.0)
    );
    // Includes the carrier contribution even though its fraction rounds to 1.
    near(
        result.entropy_change_j_per_k(),
        R * 1e-100 * (1.0 - 1e-100_f64.ln()),
        3e-14,
        0.0,
    );
}

#[test]
fn g1_binary_reduction_agrees_with_existing_caloric_owner() {
    for b in [-120.0, 120.0] {
        let model = model(&[(0, 1, 0.2)], &[(0, 1, b)], &[100.0, 180.0], 100.0);
        let pair = BinaryRegularSolutionV1::new(
            [ids(2)[0].clone(), ids(2)[1].clone()],
            0.2,
            b,
            Temperature::new(100.0),
            Temperature::new(1000.0),
        )
        .unwrap();
        let binary = BinaryLiquidCaloricModelV1::new(
            pair,
            [
                model.heat_capacities()[0].clone(),
                model.heat_capacities()[1].clone(),
            ],
        )
        .unwrap();
        let binary_controls = BinaryCaloricSolveV1 {
            temperature_bracket: controls().temperature_bracket,
            temperature_tolerance_k: controls().temperature_tolerance_k,
            energy_tolerance_j: controls().energy_tolerance_j,
            maximum_iterations: controls().maximum_iterations,
            coexistence: BinaryCoexistenceOptionsV1::default(),
        };
        let left = parcel(&model, 320.0, &[1.0, 0.0]);
        let right = parcel(&model, 470.0, &[0.0, 2.0]);
        let binary_left = binary.homogeneous(left.temperature(), [1.0, 0.0]).unwrap();
        let binary_right = binary.homogeneous(right.temperature(), [0.0, 2.0]).unwrap();
        let result = model
            .mix(&left, &right, 500.0, controls(), || false)
            .unwrap();
        let reference = binary
            .mix(&binary_left, &binary_right, 500.0, binary_controls, || {
                false
            })
            .unwrap();
        near(
            result.final_state().temperature().value(),
            reference.final_state().temperature().value(),
            0.0,
            1e-9,
        );
        near(
            result.final_state().enthalpy_j(),
            reference.final_state().enthalpy_j(),
            2e-14,
            1e-7,
        );
        near(
            result.entropy_change_j_per_k(),
            reference.entropy_change_j_per_k(),
            2e-12,
            1e-9,
        );
    }
}

#[test]
fn g1_isolated_homogeneous_mix_refuses_an_entropy_decreasing_endpoint() {
    let model = model(&[(0, 1, 10.0)], &[], &[100.0, 100.0, 100.0], 100.0);
    let left = parcel(&model, 500.0, &[1.0, 0.0, 0.0]);
    let right = parcel(&model, 500.0, &[0.0, 1.0, 0.0]);
    let original = (left.clone(), right.clone());
    match model.mix(&left, &right, 0.0, controls(), || false) {
        Err(MulticomponentCaloricErrorV1::EntropyDecrease { change_j_per_k }) => {
            near(change_j_per_k, 2.0 * R * 2.0_f64.ln() - 5.0 * R, 3e-14, 0.0);
        }
        other => panic!("isolated constrained mixture must refuse: {other:?}"),
    }
    assert_eq!((left, right), original);
}

#[test]
fn g1_source_identity_bounds_and_cancellation_are_transactional() {
    let model = model(&[], &[], &[100.0, 200.0, 300.0], 100.0);
    let initial = parcel(&model, 350.0, &[1.0, 2.0, 0.0]);
    let original = initial.clone();
    assert!(
        MulticomponentLiquidCaloricModelV1::new(
            potential(3, &[], &[]),
            model.heat_capacities()[..2].to_vec()
        )
        .is_err()
    );
    let mut wrong_species = model.heat_capacities().to_vec();
    wrong_species[1] = wrong_species[0].clone();
    assert!(
        MulticomponentLiquidCaloricModelV1::new(potential(3, &[], &[]), wrong_species).is_err()
    );
    let mut wrong_datum = model.heat_capacities().to_vec();
    wrong_datum[1] = source(ids(3)[1].clone(), 300.0, &[(100.0, 200.0), (1000.0, 200.0)]);
    assert!(MulticomponentLiquidCaloricModelV1::new(potential(3, &[], &[]), wrong_datum).is_err());
    assert!(
        model
            .homogeneous(Temperature::new(99.0), &[1.0, 2.0, 0.0], || false)
            .is_err()
    );
    let alien_model = super::tests::model(&[(0, 1, 0.1)], &[], &[100.0, 200.0, 300.0], 100.0);
    let alien = parcel(&alien_model, 350.0, &[1.0, 2.0, 0.0]);
    assert!(
        model
            .mix(&initial, &alien, 0.0, controls(), || false)
            .is_err()
    );
    assert!(matches!(
        model.add_heat(&initial, 1e10, controls(), || false),
        Err(MulticomponentCaloricErrorV1::Caloric(
            BinaryCaloricErrorV1::HeatOutsideBracket { .. }
        ))
    ));
    assert!(matches!(
        model.add_heat(
            &initial,
            12_345.0,
            MulticomponentCaloricSolveV1 {
                maximum_iterations: 1,
                ..controls()
            },
            || false
        ),
        Err(MulticomponentCaloricErrorV1::Caloric(
            BinaryCaloricErrorV1::IterationLimit
        ))
    ));
    assert!(
        model
            .add_heat(&initial, f64::NAN, controls(), || false)
            .is_err()
    );
    assert!(
        model
            .add_heat(
                &initial,
                0.0,
                MulticomponentCaloricSolveV1 {
                    temperature_bracket: [Temperature::new(200.0), Temperature::new(1100.0)],
                    ..controls()
                },
                || false
            )
            .is_err()
    );
    let mut calls = 0;
    model
        .add_heat(&initial, 12_345.0, controls(), || {
            calls += 1;
            false
        })
        .unwrap();
    for cutoff in [0, 1, calls / 2, calls - 1] {
        let mut current = 0;
        assert!(
            model
                .add_heat(&initial, 12_345.0, controls(), || {
                    let stop = current >= cutoff;
                    current += 1;
                    stop
                })
                .is_err()
        );
        assert_eq!(initial, original);
    }
}
