use super::*;
use crate::solution::caloric::{BinaryLiquidCaloricModelV1, PureLiquidHeatCapacityV1};
use crate::solution::{BinaryCoexistenceOptionsV1, BinaryRegularSolutionV1};
use crate::{MolarThermalQuantityV1, SpeciesId};
use fs_qty::Temperature;

fn model(a: f64, b: f64, cp: [[f64; 2]; 2]) -> BinaryLiquidCaloricModelV1 {
    let species = [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()];
    let pair = BinaryRegularSolutionV1::new(
        species.clone(),
        a,
        b,
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap();
    let sources = std::array::from_fn(|i| {
        PureLiquidHeatCapacityV1::new(
            species[i].clone(),
            Temperature::new(100.0),
            vec![
                (
                    Temperature::new(100.0),
                    MolarThermalQuantityV1::new(cp[i][0]),
                ),
                (
                    Temperature::new(1000.0),
                    MolarThermalQuantityV1::new(cp[i][1]),
                ),
            ],
        )
        .unwrap()
    });
    BinaryLiquidCaloricModelV1::new(pair, sources).unwrap()
}

fn state(model: &BinaryLiquidCaloricModelV1, t: f64, n: [f64; 2]) -> BinaryCaloricStateV1 {
    model
        .equilibrate(
            Temperature::new(t),
            n,
            BinaryCoexistenceOptionsV1 {
                relative_tolerance: 4.0 * f64::EPSILON,
                maximum_iterations: 256,
            },
            || false,
        )
        .unwrap()
}

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_homogeneous_capacity_uses_each_original_cp_at_the_actual_temperature() {
    let model = model(0.0, 1000.0, [[20.0, 110.0], [40.0, 220.0]]);
    let homogeneous = model
        .homogeneous(Temperature::new(400.0), [2.0, 3.0])
        .unwrap();
    assert!(homogeneous.equilibrium().is_none());
    // Original source slopes are 0.1 and 0.2 J/(mol K^2).
    assert_eq!(homogeneous.isobaric_heat_capacity_j_per_k().unwrap(), 400.0);
    assert_eq!(
        homogeneous
            .isobaric_heat_capacity_response()
            .unwrap()
            .qualification(),
        BinaryCaloricHeatCapacityQualificationV1::Homogeneous,
    );
    let hot = state(&model, 600.0, [2.0, 3.0]);
    assert_eq!(hot.equilibrium().unwrap().phases().len(), 1);
    assert_eq!(hot.isobaric_heat_capacity_j_per_k().unwrap(), 560.0);
    let pure = state(&model, 400.0, [1.0, 0.0]);
    assert_eq!(pure.isobaric_heat_capacity_j_per_k().unwrap(), 50.0);
    assert_eq!(
        model.heat_capacities()[1]
            .heat_capacity_at(Temperature::new(400.0))
            .unwrap()
            .value(),
        100.0,
    );
    assert!(
        model.heat_capacities()[0]
            .heat_capacity_at(Temperature::new(99.0))
            .is_err()
    );
    assert!(
        state(&model, 400.0, [2.0, 3.0])
            .isobaric_heat_capacity_j_per_k()
            .unwrap()
            > 400.0
    );
}

#[test]
fn g1_equilibrium_cp_matches_independent_enthalpy_and_entropy_derivatives() {
    for (a, b, t) in [
        (0.0, 1000.0, 400.0),
        (0.0, 1000.0, 490.0),
        (5.0, -900.0, 450.0),
    ] {
        let model = model(a, b, [[20.0, 110.0], [40.0, 220.0]]);
        let center = state(&model, t, [2.0, 3.0]);
        let response = center.isobaric_heat_capacity_response().unwrap();
        assert_eq!(
            response.qualification(),
            BinaryCaloricHeatCapacityQualificationV1::Coexisting
        );
        let cp = center.isobaric_heat_capacity_j_per_k().unwrap();
        let sensible = 8.0 * (10.0 + 0.1 * t);
        assert!(cp > sensible + 1.0);
        let step = 0.005;
        let cold = state(&model, t - step, [2.0, 3.0]);
        let hot = state(&model, t + step, [2.0, 3.0]);
        // Differentiate the pre-existing enthalpy and entropy owners. Neither
        // endpoint calculation calls the new heat-capacity response.
        let enthalpy_derivative = (hot.enthalpy_j() - cold.enthalpy_j()) / (2.0 * step);
        let entropy_derivative = model.entropy_change_j_per_k(&cold, &hot).unwrap() / (2.0 * step);
        near(cp, enthalpy_derivative, 3e-6, 0.0);
        near(cp, t * entropy_derivative, 3e-6, 0.0);
    }
}

#[test]
fn g1_capacity_preserves_inventory_scaling_and_species_exchange() {
    let original = model(0.0, 1000.0, [[20.0, 110.0], [40.0, 220.0]]);
    let exchanged = model(0.0, 1000.0, [[40.0, 220.0], [20.0, 110.0]]);
    let cp = state(&original, 400.0, [2.0, 3.0])
        .isobaric_heat_capacity_j_per_k()
        .unwrap();
    let swapped = state(&exchanged, 400.0, [3.0, 2.0])
        .isobaric_heat_capacity_j_per_k()
        .unwrap();
    near(swapped, cp, 0.0, 3e-15);
    for scale in [1e-12, 1e12] {
        let scaled = state(&original, 400.0, [2.0 * scale, 3.0 * scale])
            .isobaric_heat_capacity_j_per_k()
            .unwrap();
        near(scaled, scale * cp, 0.0, 5e-15);
    }
}

#[test]
fn g1_temperature_independent_interactions_have_no_redistribution_capacity() {
    for a in [0.0, 2.0, 3.0] {
        let model = model(a, 0.0, [[20.0, 110.0], [40.0, 220.0]]);
        let center = state(&model, 400.0, [1.0, 1.0]);
        assert_eq!(center.isobaric_heat_capacity_j_per_k().unwrap(), 150.0);
        if let Some(gap) = model
            .thermodynamics()
            .coexistence(
                Temperature::new(400.0),
                BinaryCoexistenceOptionsV1::default(),
                || false,
            )
            .unwrap()
        {
            let fractions = gap.lower_a().mole_fractions();
            let endpoint = model
                .equilibrate(
                    Temperature::new(400.0),
                    fractions,
                    BinaryCoexistenceOptionsV1::default(),
                    || false,
                )
                .unwrap();
            assert_eq!(
                endpoint.equilibrium().unwrap().feed_stability(),
                BinaryMixingStabilityV1::Binodal
            );
            near(
                endpoint.isobaric_heat_capacity_j_per_k().unwrap(),
                50.0 * fractions[0] + 100.0 * fractions[1],
                1e-13,
                0.0,
            );
        }
    }
}

#[test]
fn g1_near_critical_split_capacity_retains_the_finite_one_sided_jump() {
    let model = model(0.0, 1000.0, [[100.0; 2]; 2]);
    let below = state(&model, 500.0_f64.next_down(), [1.0, 1.0]);
    assert!(below.equilibrium().unwrap().feed().chi() > 2.0);
    assert_eq!(below.equilibrium().unwrap().phases().len(), 2);
    // Mean-field split-side limit: 2 R (B/Tc)^2 * (3/8) = 3 R.
    near(
        below.isobaric_heat_capacity_j_per_k().unwrap(),
        200.0 + 3.0 * R,
        0.0,
        3e-15,
    );
    let critical = state(&model, 500.0, [1.0, 1.0]);
    assert!(critical.isobaric_heat_capacity_j_per_k().is_err());
    let tangent = critical.isobaric_heat_capacity_response().unwrap();
    assert_eq!(tangent.heat_capacity_j_per_k(), 200.0);
    assert_eq!(
        tangent.qualification(),
        BinaryCaloricHeatCapacityQualificationV1::HomogeneousPhaseBoundaryTangent
    );
    let above = state(&model, 500.0_f64.next_up(), [1.0, 1.0]);
    assert_eq!(above.isobaric_heat_capacity_j_per_k().unwrap(), 200.0);
    // Criticality is also a composition condition. This off-center feed has
    // an open homogeneous neighborhood at the very same critical temperature.
    let off_center = state(&model, 500.0, [2.0, 3.0]);
    assert_eq!(off_center.isobaric_heat_capacity_j_per_k().unwrap(), 500.0);
    assert_eq!(
        off_center
            .isobaric_heat_capacity_response()
            .unwrap()
            .qualification(),
        BinaryCaloricHeatCapacityQualificationV1::Homogeneous
    );
}

#[test]
fn g3_binodal_and_unresolved_bracket_use_only_an_explicit_newton_tangent() {
    let model = model(0.0, 1000.0, [[100.0; 2]; 2]);
    let controls = BinaryCoexistenceOptionsV1 {
        relative_tolerance: 1e-3,
        maximum_iterations: 128,
    };
    let temperature = Temperature::new(400.0);
    let gap = model
        .thermodynamics()
        .coexistence(temperature, controls, || false)
        .unwrap()
        .unwrap();
    let endpoint = model
        .equilibrate(
            temperature,
            gap.lower_a().mole_fractions(),
            controls,
            || false,
        )
        .unwrap();
    assert_eq!(
        endpoint.equilibrium().unwrap().feed_stability(),
        BinaryMixingStabilityV1::Binodal
    );
    assert!(endpoint.isobaric_heat_capacity_j_per_k().is_err());
    let response = endpoint.isobaric_heat_capacity_response().unwrap();
    assert_eq!(response.heat_capacity_j_per_k(), 100.0);
    assert_eq!(
        response.qualification(),
        BinaryCaloricHeatCapacityQualificationV1::HomogeneousPhaseBoundaryTangent
    );

    let [low, high] = gap.log_odds_bracket();
    let x_lower = dilute_fraction(high);
    let x_upper = dilute_fraction(low);
    let x_unresolved = 0.25 * x_lower + 0.75 * x_upper;
    let unresolved = model
        .equilibrate(
            temperature,
            [x_unresolved, 1.0 - x_unresolved],
            controls,
            || false,
        )
        .unwrap();
    assert!(unresolved.isobaric_heat_capacity_j_per_k().is_err());
    let tangent = unresolved.isobaric_heat_capacity_response().unwrap();
    assert_eq!(tangent.heat_capacity_j_per_k(), 100.0);
    assert_eq!(
        tangent.qualification(),
        BinaryCaloricHeatCapacityQualificationV1::HomogeneousBinodalBracketTangent
    );

    let x_outside = 0.99 * x_lower;
    let outside = model
        .equilibrate(temperature, [x_outside, 1.0 - x_outside], controls, || {
            false
        })
        .unwrap();
    near(
        outside.isobaric_heat_capacity_j_per_k().unwrap(),
        100.0,
        2e-14,
        0.0,
    );
    let x_inside = x_upper + 0.01 * (0.5 - x_upper);
    let inside = model
        .equilibrate(temperature, [x_inside, 1.0 - x_inside], controls, || false)
        .unwrap();
    assert!(inside.isobaric_heat_capacity_j_per_k().unwrap() > 100.0);
}

#[test]
fn g3_capacity_arithmetic_errors_are_not_replaced_by_boundary_tangents() {
    for (cp, amounts) in [(1e-300, [1e-100; 2]), (1e308, [1.0; 2])] {
        let model = model(0.0, 0.0, [[cp; 2]; 2]);
        let state = state(&model, 100.0, amounts);
        assert!(matches!(
            state.isobaric_heat_capacity_j_per_k(),
            Err(BinaryCaloricErrorV1::Arithmetic(_))
        ));
        assert!(matches!(
            state.isobaric_heat_capacity_response(),
            Err(BinaryCaloricErrorV1::Arithmetic(_))
        ));
    }
    let model = model(3.0, 1e-200, [[100.0; 2]; 2]);
    let state = state(&model, 400.0, [1.0; 2]);
    assert!(matches!(
        state.isobaric_heat_capacity_response(),
        Err(BinaryCaloricErrorV1::Arithmetic(_))
    ));
}
