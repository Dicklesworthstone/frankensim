//! Formula and coupling regressions, not experimental validation.

use fs_convection::{
    CorrelationError, CorrelationId, CorrelationInputs, DiscrepancyBasis, ThermalConductivity,
    ThermalDirection, evaluate,
};
use fs_qty::{Length, Temperature};

const CARD: CorrelationId = CorrelationId::SiederTateLaminar;

fn point(re: f64, pr: f64, length: f64, viscosity_ratio: f64) -> CorrelationInputs {
    CorrelationInputs::forced(re, pr)
        .with_length_ratio(length)
        .with_viscosity_ratio(viscosity_ratio)
}

fn nu(values: [f64; 4]) -> f64 {
    evaluate(CARD, point(values[0], values[1], values[2], values[3]))
        .unwrap()
        .evidence()
        .value
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn source_formula_retains_mean_entrance_length_and_viscosity_ratio() {
    // Independent cube-root limits: Gz=27, 1000 and 1e6, with no wall correction.
    for (values, expected) in [
        ([270.0, 10.0, 100.0, 1.0], 5.58),
        ([1000.0, 10.0, 10.0, 1.0], 18.6),
        ([1000.0, 1000.0, 1.0, 1.0], 186.0),
    ] {
        near(nu(values), expected, 3e-13);
    }
    let base = nu([1000.0, 10.0, 10.0, 1.0]);
    near(nu([1000.0, 10.0, 10.0, 4.0]) / base, 1.214_194_884_395_047, 3e-13);
    near(nu([1000.0, 10.0, 10.0, 0.25]) / base, 0.823_591_017_267_573_1, 3e-13);
    // An eightfold heated length lowers the MEAN coefficient by a factor of two.
    near(nu([1000.0, 10.0, 80.0, 1.0]), 0.5 * base, 3e-13);
    let evaluated = evaluate(CARD, point(1000.0, 10.0, 10.0, 1.0)).unwrap();
    assert!(evaluated.card().source.citation.contains("Table 2.4, p. 27"));
    assert_eq!(evaluated.card().discrepancy_basis, DiscrepancyBasis::EngineeringAllowance);
}

#[test]
fn wall_viscosity_controls_heating_and_cooling_without_a_direction_default() {
    for ratio in [0.25, 1.0, 4.0] {
        let inputs = point(1000.0, 10.0, 10.0, ratio);
        let heating = evaluate(CARD, inputs).unwrap();
        let cooling = evaluate(CARD, inputs.with_direction(ThermalDirection::CoolingFluid)).unwrap();
        assert_eq!(heating.evidence().value.to_bits(), cooling.evidence().value.to_bits());
        assert_eq!(heating.groups()["mu_bulk_over_surface"], ratio);
        assert_eq!(
            heating.evidence().value.to_bits(),
            evaluate(CARD, inputs).unwrap().evidence().value.to_bits()
        );
    }
    assert!(matches!(
        evaluate(CARD, CorrelationInputs::forced(1000.0, 10.0).with_length_ratio(10.0)),
        Err(CorrelationError::OutOfDomain { correlation: CARD, .. })
    ));
}

#[test]
fn coupled_entrance_gate_refuses_fully_developed_states_without_clipping() {
    // At the same Re/Pr/L, cold-wall viscosity can leave the entrance envelope.
    assert!(evaluate(CARD, point(100.0, 1.0, 10.0, 1.0)).is_ok());
    let Err(CorrelationError::OutOfDomain { violations, .. }) =
        evaluate(CARD, point(100.0, 1.0, 10.0, 0.1))
    else {
        panic!("the wall correction must participate in the entrance gate");
    };
    assert!(violations.iter().any(|v| v.axis == "Sieder_Tate_entry_factor"));
    assert!(evaluate(CARD, point(79.0, 1.0, 10.0, 1.0)).is_err());
    assert!(evaluate(CARD, point(81.0, 1.0, 10.0, 1.0)).is_ok());
    assert!(evaluate(CARD, point(1000.0, 10.0, 10_000.0, 1.0)).is_err());
}

#[test]
fn source_open_bounds_and_laminar_flow_are_admitted_explicitly() {
    for inputs in [
        point(1000.0, 0.48, 1.0, 1.0),
        point(1000.0, 16_700.0, 100.0, 1.0),
        point(1000.0, 10.0, 10.0, 0.0044),
        point(1000.0, 10.0, 10.0, 9.75),
        point(2000.0_f64.next_up(), 10.0, 10.0, 1.0),
        point(1000.0, 10.0, 1.0_f64.next_down(), 1.0),
    ] {
        assert!(matches!(evaluate(CARD, inputs), Err(CorrelationError::OutOfDomain { .. })));
        assert!(inputs.evaluate_with_sensitivities(CARD).is_err());
    }
    for ratio in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(matches!(
            evaluate(CARD, point(1000.0, 10.0, 10.0, ratio)),
            Err(CorrelationError::InvalidGroup { axis: "mu_bulk_over_surface", .. })
        ));
    }
}

#[test]
fn all_independent_input_partials_match_centered_physical_perturbations() {
    let values = [800.0, 10.0, 40.0, 2.0];
    let result = point(values[0], values[1], values[2], values[3])
        .evaluate_with_sensitivities(CARD)
        .unwrap();
    for (index, axis, elasticity) in [
        (0, "Re", 1.0 / 3.0),
        (1, "Pr", 1.0 / 3.0),
        (2, "L_over_Dh", -1.0 / 3.0),
        (3, "mu_bulk_over_surface", 0.14),
    ] {
        let derivative = result.evidence().sensitivity.d_qoi[axis];
        near(derivative, result.evidence().value * elasticity / values[index], 3e-13);
        let h = 1e-5 * values[index];
        let mut plus = values;
        let mut minus = values;
        plus[index] += h;
        minus[index] -= h;
        near(derivative, (nu(plus) - nu(minus)) / (2.0 * h), 2e-8);
    }
    near(result.reynolds_elasticity().unwrap().unwrap(), 1.0 / 3.0, 1e-14);
}

#[test]
fn mean_coefficient_lowers_into_the_existing_conduction_robin_row() {
    let result = evaluate(CARD, point(1000.0, 10.0, 10.0, 1.0)).unwrap();
    let k = ThermalConductivity::new(0.14);
    let diameter = Length::new(0.01);
    let h = result.heat_transfer_coefficient(k, diameter).unwrap();
    near(h.value.value(), 260.4, 3e-13);
    let boundary = result.robin_boundary(k, diameter, Temperature::new(300.0)).unwrap();
    assert_eq!(boundary.coefficient().value, h.value);
    assert_eq!(boundary.coefficient().model.cards, [CARD.name()]);
}
