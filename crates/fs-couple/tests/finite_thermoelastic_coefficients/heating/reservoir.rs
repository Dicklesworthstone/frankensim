//! Actual prescribed-bath exchange with finite original caloric storage.

use super::*;
use fs_thermal::enthalpy::Reservoir;

fn constant_capacity_card() -> MaterialCard {
    let original = card_with_expansion([0.0; 2], 1.0, [0.0; 3]);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key == PropertyKey::new("specific_heat_capacity", CP_DIMS)
            && let PropertyValue::Curve { knots, .. } = &mut claim.value
        {
            for (_, value) in knots {
                *value = 12.0;
            }
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, original.models().to_vec()).unwrap()
}

#[test]
fn reservoir_matches_constant_capacity_backward_euler_and_separates_temporal_entropy() {
    let card = constant_capacity_card();
    let initial = prepared(&card, &source_curve(&card, 0.0));
    let capacity = 24.0;
    assert_eq!(initial.heat_capacity_j_per_k(), capacity);
    for bath in [340.0, 260.0] {
        // Bath temperatures are deliberately outside the material source path.
        // Area is metadata; G is already total W/K and must not be scaled again.
        let reservoir = Reservoir {
            area_m2: 0.125,
            conductance_w_k: 3.0,
            temperature_k: bath,
        };
        let duration = 2.0;
        let report = initial
            .exchange_with_reservoir(reservoir, duration, controls(), || false)
            .unwrap();
        let expected = (capacity * 300.0 + duration * 3.0 * bath) / (capacity + duration * 3.0);
        let heat = duration * 3.0 * (bath - expected);
        let entropy = capacity * (expected / 300.0).ln();
        near(report.final_state().temperature_k(), expected);
        assert_eq!(report.reservoir(), &reservoir);
        assert_eq!(report.duration_s(), duration);
        near(report.heat_balance().heat_j(), heat);
        near(
            report.heat_balance().internal_energy_change_j(),
            capacity * (expected - 300.0),
        );
        near(report.endpoint_power_into_solid_w(), heat / duration);
        near(report.heat_balance().entropy_change_j_per_k(), entropy);
        near(
            report.endpoint_conduction_entropy_j_per_k(),
            heat * (1.0 / expected - 1.0 / bath),
        );
        near(
            report.temporal_entropy_defect_j_per_k(),
            entropy - heat / expected,
        );
        near(
            report.combined_entropy_change_j_per_k(),
            entropy - heat / bath,
        );
        assert!(report.endpoint_conduction_entropy_j_per_k() > 0.0);
        assert!(report.temporal_entropy_defect_j_per_k() > 0.0);
        assert!(report.heat_balance().energy_residual_j().abs() <= 1e-8);

        let first = initial
            .exchange_with_reservoir(reservoir, duration / 2.0, controls(), || false)
            .unwrap();
        let second = first
            .final_state()
            .exchange_with_reservoir(reservoir, duration / 2.0, controls(), || false)
            .unwrap();
        let exact = bath + (300.0 - bath) * (-3.0 * duration / capacity).exp();
        assert!((second.final_state().temperature_k() - exact).abs() < (expected - exact).abs());
        assert!(
            first.temporal_entropy_defect_j_per_k() + second.temporal_entropy_defect_j_per_k()
                < report.temporal_entropy_defect_j_per_k()
        );
    }
    assert_eq!(initial.temperature_k(), 300.0);
}

#[test]
fn reservoir_crosses_original_cp_knot_updates_mechanics_and_is_datum_independent() {
    let card = card([-12.0, 0.004], 1.0);
    let plain = prepared(&card, &source_curve(&card, 0.0))
        .at_temperature(296.0, || false)
        .unwrap();
    let shifted = prepared(&card, &source_curve(&card, 1e17))
        .at_temperature(296.0, || false)
        .unwrap();
    let (expected_heat, expected_entropy) = expected_change(296.0, 304.0);
    let reservoir = Reservoir {
        area_m2: 0.35,
        conductance_w_k: expected_heat / (2.0 * (330.0 - 304.0)),
        temperature_k: 330.0,
    };
    let heated = plain
        .exchange_with_reservoir(reservoir, 2.0, controls(), || false)
        .unwrap();
    let shifted_heat = shifted
        .exchange_with_reservoir(reservoir, 2.0, controls(), || false)
        .unwrap();
    near(heated.final_state().temperature_k(), 304.0);
    near(heated.heat_balance().heat_j(), expected_heat);
    near(
        heated.heat_balance().entropy_change_j_per_k(),
        expected_entropy,
    );
    assert_eq!(
        heated.final_state().thermodynamics(),
        shifted_heat.final_state().thermodynamics()
    );
    assert_eq!(
        heated.heat_balance().heat_j(),
        shifted_heat.heat_balance().heat_j()
    );
    assert_eq!(
        heated.combined_entropy_change_j_per_k(),
        shifted_heat.combined_entropy_change_j_per_k()
    );
    assert_ne!(
        plain.thermodynamics().first_piola_pa(),
        heated.final_state().thermodynamics().first_piola_pa()
    );
    assert!(
        heated
            .final_state()
            .adiabatic_acoustics()
            .unwrap()
            .acoustic()
            .directional_waves([1.0, 0.0, 0.0])
            .is_ok()
    );
    let cooled = heated
        .final_state()
        .exchange_with_reservoir(
            Reservoir {
                temperature_k: 270.0,
                ..reservoir
            },
            2.0,
            controls(),
            || false,
        )
        .unwrap();
    near(cooled.final_state().temperature_k(), 296.0);
    near(cooled.heat_balance().heat_j(), -expected_heat);
    for report in [&heated, &cooled] {
        let heat = report.heat_balance();
        assert!(heat.energy_residual_j().abs() <= 1e-8);
        assert!(report.endpoint_conduction_entropy_j_per_k() > 0.0);
        assert!(report.temporal_entropy_defect_j_per_k() > 0.0);
        near(
            report.combined_entropy_change_j_per_k(),
            heat.entropy_change_j_per_k() - heat.heat_j() / report.reservoir().temperature_k,
        );
    }
    assert_eq!(plain.temperature_k(), 296.0);
}

#[test]
fn reservoir_noops_and_refusals_preserve_original_state_and_bound_work() {
    let card = card([-12.0, 0.004], 1.0);
    let initial = prepared(&card, &source_curve(&card, 0.0));
    let original = initial.thermodynamics().clone();
    let reservoir = Reservoir {
        area_m2: 0.1,
        conductance_w_k: 2.0,
        temperature_k: 303.0,
    };
    for (bath, duration) in [
        (reservoir, 0.0),
        (
            Reservoir {
                conductance_w_k: 0.0,
                ..reservoir
            },
            2.0,
        ),
        (
            Reservoir {
                temperature_k: 300.0,
                ..reservoir
            },
            2.0,
        ),
    ] {
        let report = initial
            .exchange_with_reservoir(bath, duration, controls(), || false)
            .unwrap();
        assert_eq!(report.final_state().thermodynamics(), &original);
        assert_eq!(report.heat_balance().heat_j(), 0.0);
        assert_eq!(report.heat_balance().state_evaluations(), 0);
        assert_eq!(report.combined_entropy_change_j_per_k(), 0.0);
        let mut invalid = controls();
        invalid.thermal.bracket = [300.0, 300.0];
        assert!(
            initial
                .exchange_with_reservoir(bath, duration, invalid, || false)
                .is_err()
        );
    }
    let mut short = controls();
    short.maximum_state_evaluations = 2;
    assert!(matches!(
        initial.exchange_with_reservoir(reservoir, 2.0, short, || false),
        Err(FiniteThermoelasticHeatError::StateBudget { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        initial.exchange_with_reservoir(reservoir, 2.0, controls(), || {
            polls += 1;
            polls > 8
        }),
        Err(FiniteThermoelasticHeatError::Cancelled)
    ));
    for (bath, duration) in [
        (
            Reservoir {
                area_m2: -1.0,
                ..reservoir
            },
            0.0,
        ),
        (
            Reservoir {
                temperature_k: 0.0,
                ..reservoir
            },
            0.0,
        ),
        (reservoir, -1.0),
        (
            Reservoir {
                temperature_k: 400.0,
                ..reservoir
            },
            1e6,
        ),
        (
            Reservoir {
                conductance_w_k: 1e-300,
                ..reservoir
            },
            1e-300,
        ),
        (
            Reservoir {
                conductance_w_k: 1e-310,
                ..reservoir
            },
            1.0,
        ),
    ] {
        assert!(
            initial
                .exchange_with_reservoir(bath, duration, controls(), || false)
                .is_err()
        );
    }
    assert_eq!(initial.thermodynamics(), &original);
}
