//! G0/G1 adiabatic field balances, entropy admission and thermal time refinement.

use super::{
    BinaryMiscibilityThermalError as Error, BinaryMiscibilityThermalField,
    BinaryMiscibilityThermalSolve,
};
use crate::miscibility::{
    BinaryMiscibilityCell, BinaryMiscibilityError, BinaryMiscibilityFace, BinaryMiscibilityNetwork,
    BinaryMiscibilitySolve, BinaryMiscibilityState,
};
use fs_qty::{AmountConcentration, Temperature};
use fs_thermochem::solution::caloric::{BinaryLiquidCaloricModelV1, PureLiquidHeatCapacityV1};
use fs_thermochem::solution::{BinaryRegularSolutionV1, gradient::BinarySquareGradientV1};
use fs_thermochem::{MolarThermalQuantityV1, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn network(b: f64, kappa: f64) -> BinaryMiscibilityNetwork {
    let pair = BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        0.0,
        b,
        Temperature::new(250.0),
        Temperature::new(600.0),
    )
    .unwrap();
    let model = BinarySquareGradientV1::new(
        pair,
        Temperature::new(300.0),
        AmountConcentration::new(1.0),
        kappa,
    )
    .unwrap();
    BinaryMiscibilityNetwork::new(
        model,
        vec![BinaryMiscibilityCell { volume_m3: 1.0 }; 2],
        vec![BinaryMiscibilityFace {
            cells: [0, 1],
            area_m2: 1.0,
            distance_m: 1.0,
            diffusivity_m2_s: 0.4,
        }],
        16 * 1024 * 1024,
        || false,
    )
    .unwrap()
}

fn state(net: &BinaryMiscibilityNetwork, fractions: [f64; 2]) -> BinaryMiscibilityState {
    net.state_from_mole_fractions(fractions.map(|x| [x, 1.0 - x]).to_vec(), || false)
        .unwrap()
}

fn calorics(
    net: &BinaryMiscibilityNetwork,
    capacities: [f64; 2],
    domain: [f64; 2],
) -> BinaryLiquidCaloricModelV1 {
    let pair = net.model().thermodynamics();
    let sources = [0, 1].map(|species| {
        PureLiquidHeatCapacityV1::new(
            pair.species()[species].clone(),
            Temperature::new(domain[0]),
            domain
                .map(|t| {
                    (
                        Temperature::new(t),
                        MolarThermalQuantityV1::new(capacities[species]),
                    )
                })
                .to_vec(),
        )
        .unwrap()
    });
    BinaryLiquidCaloricModelV1::new(pair.clone(), sources).unwrap()
}

fn field(
    b: f64,
    kappa: f64,
    fractions: [f64; 2],
    capacities: [f64; 2],
) -> BinaryMiscibilityThermalField {
    let net = network(b, kappa);
    let initial = state(&net, fractions);
    let caloric = calorics(&net, capacities, [250.0, 600.0]);
    net.bind_caloric_field(&initial, &caloric, || false)
        .unwrap()
}

fn solve() -> BinaryMiscibilitySolve {
    BinaryMiscibilitySolve {
        absolute_transfer_tolerance_mol: 1e-18,
        relative_transfer_tolerance: 1e-12,
        ..BinaryMiscibilitySolve::default()
    }
}

fn thermal() -> BinaryMiscibilityThermalSolve {
    BinaryMiscibilityThermalSolve {
        temperature_bracket: [Temperature::new(250.0), Temperature::new(600.0)],
        temperature_tolerance_k: 1e-11,
        energy_tolerance_j: 1e-9,
        entropy_tolerance_j_per_k: 1e-10,
        maximum_iterations: 96,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance {tolerance:e}"
    );
}

fn gradient_energy(state: &BinaryMiscibilityState) -> f64 {
    let difference = state.mole_fractions()[0][0] - state.mole_fractions()[1][0];
    0.5 * state.model().gradient_coefficient_j_m() * difference * difference
}

fn mixing_entropy(state: &BinaryMiscibilityState) -> f64 {
    -R * state
        .mole_fractions()
        .iter()
        .flatten()
        .map(|&x| if x == 0.0 { 0.0 } else { x * x.ln() })
        .sum::<f64>()
}

#[test]
fn demixing_heats_the_actual_field_and_closes_energy_species_and_one_physical_clock() {
    let initial = field(900.0, 0.02 * R * 300.0, [0.1, 0.9], [20.0, 30.0]);
    let original = initial.state().clone();
    let candidate = initial
        .adiabatic_candidate(1.0, solve(), thermal(), || false)
        .unwrap();
    let state = candidate.field().state();
    let step = candidate.step();
    let temperature = state.model().temperature().value();
    assert!(temperature > 300.0);
    assert!(state.mole_fractions()[0][0] < original.mole_fractions()[0][0]);
    let bulk = |s: &BinaryMiscibilityState| {
        R * 900.0 * s.mole_fractions().iter().map(|x| x[0] * x[1]).sum::<f64>()
    };
    let mixing_change = bulk(state) - bulk(&original);
    let gradient_change = gradient_energy(state) - gradient_energy(&original);
    near(step.mixing_enthalpy_change_j, mixing_change, 2e-10);
    near(step.interfacial_energy_change_j, gradient_change, 2e-11);
    assert!(step.mixing_enthalpy_change_j < 0.0);
    assert!(step.interfacial_energy_change_j > 0.0);
    near(
        temperature,
        300.0 - (mixing_change + gradient_change) / 50.0,
        3e-11,
    );
    near(
        step.sensible_enthalpy_change_j,
        50.0 * (temperature - 300.0),
        2e-10,
    );
    near(
        step.energy_residual_j,
        step.mixing_enthalpy_change_j
            + step.interfacial_energy_change_j
            + step.sensible_enthalpy_change_j,
        2e-12,
    );
    assert!(step.energy_residual_j.abs() <= thermal().energy_tolerance_j);
    near(
        step.mixing_entropy_change_j_per_k,
        mixing_entropy(state) - mixing_entropy(&original),
        2e-13,
    );
    near(
        step.sensible_entropy_change_j_per_k,
        50.0 * (temperature / 300.0).ln(),
        2e-12,
    );
    near(
        step.entropy_change_j_per_k,
        step.mixing_entropy_change_j_per_k + step.sensible_entropy_change_j_per_k,
        2e-13,
    );
    assert!(step.entropy_change_j_per_k > 0.0);
    let q = step.diffusion.transferred_a_mol[0];
    for cell in 0..2 {
        for species in 0..2 {
            let signed = if cell == species { -q } else { q };
            near(
                state.component_moles()[cell][species] - original.component_moles()[cell][species],
                signed,
                2e-14,
            );
        }
        near(state.component_moles()[cell].iter().sum(), 1.0, 2e-14);
    }
    for species in 0..2 {
        near(
            state.component_moles().iter().map(|n| n[species]).sum(),
            1.0,
            2e-14,
        );
    }
    assert_eq!(state.elapsed_s(), 1.0);
    assert_eq!(step.diffusion.duration_s, 1.0);
    assert_eq!(
        step.temperature_change.temperature_before,
        Temperature::new(300.0)
    );
    assert_eq!(
        step.temperature_change.temperature_after,
        state.model().temperature()
    );
    assert_eq!(candidate.field().caloric_model(), initial.caloric_model());
    let next = candidate
        .field()
        .adiabatic_candidate(0.5, solve(), thermal(), || false)
        .unwrap();
    assert_eq!(next.field().state().elapsed_s(), 1.5);
    assert_eq!(next.field().caloric_model(), initial.caloric_model());
    assert_eq!(initial.state(), &original);
}

#[test]
fn relaxing_an_ideal_interface_heats_the_field_from_gradient_energy_alone() {
    let initial = field(0.0, 0.4 * R * 300.0, [0.2, 0.8], [10.0; 2]);
    let candidate = initial
        .adiabatic_candidate(0.2, solve(), thermal(), || false)
        .unwrap();
    let step = candidate.step();
    let state = candidate.field().state();
    let released = gradient_energy(initial.state()) - gradient_energy(state);
    assert!(released > 0.0);
    assert_eq!(step.mixing_enthalpy_change_j, 0.0);
    near(step.interfacial_energy_change_j, -released, 2e-11);
    near(
        step.sensible_enthalpy_change_j,
        released,
        thermal().energy_tolerance_j,
    );
    near(
        state.model().temperature().value(),
        300.0 + released / 20.0,
        6e-11,
    );
    assert!(step.mixing_entropy_change_j_per_k > 0.0);
    assert!(step.sensible_entropy_change_j_per_k > 0.0);
    assert!(step.entropy_change_j_per_k > 0.0);
}

#[test]
fn ideal_mixing_without_gradient_energy_keeps_temperature_exactly_constant() {
    let initial = field(0.0, 0.0, [1.0, 0.0], [3.0, 70.0]);
    let candidate = initial
        .adiabatic_candidate(0.1, solve(), thermal(), || false)
        .unwrap();
    let step = candidate.step();
    assert_eq!(
        candidate.field().state().model().temperature(),
        Temperature::new(300.0)
    );
    assert_eq!(step.mixing_enthalpy_change_j, 0.0);
    assert_eq!(step.interfacial_energy_change_j, 0.0);
    assert_eq!(step.sensible_enthalpy_change_j, 0.0);
    assert_eq!(step.energy_residual_j, 0.0);
    assert_eq!(step.sensible_entropy_change_j_per_k, 0.0);
    assert!(step.mixing_entropy_change_j_per_k > 0.0);
    for fractions in candidate.field().state().mole_fractions() {
        assert!(fractions.iter().all(|&x| x > 0.0 && x < 1.0));
    }
    assert_eq!(initial.state().mole_fractions(), &[[1.0, 0.0], [0.0, 1.0]]);
}

#[test]
fn energy_closure_alone_cannot_admit_an_entropy_losing_step_but_refinement_can() {
    let initial = field(900.0, 0.001 * R * 300.0, [0.1, 0.9], [0.5; 2]);
    let original = initial.state().clone();
    // Constant total Cp=1 J/K. The old-temperature diffusion step lowers F,
    // and its energy root lies near 446.4 K, yet total entropy loses 0.0568 J/K.
    match initial.adiabatic_candidate(1.0, solve(), thermal(), || false) {
        Err(Error::Entropy {
            change_j_per_k,
            tolerance_j_per_k,
        }) => {
            assert!(change_j_per_k < -0.05 && change_j_per_k > -0.06);
            assert_eq!(tolerance_j_per_k, thermal().entropy_tolerance_j_per_k);
        }
        Err(error) => panic!("expected the physical entropy refusal, got {error:?}"),
        Ok(_) => panic!("energy closure alone admitted an entropy-losing adiabatic step"),
    }
    assert_eq!(initial.state(), &original);
    let mut refined = initial;
    for _ in 0..100 {
        let candidate = refined
            .adiabatic_candidate(0.01, solve(), thermal(), || false)
            .unwrap();
        assert!(candidate.step().entropy_change_j_per_k >= -thermal().entropy_tolerance_j_per_k);
        refined = candidate.into_parts().0;
    }
    near(refined.state().elapsed_s(), 1.0, 2e-14);
    let temperature = refined.state().model().temperature().value();
    assert!(temperature > 320.0 && temperature < 330.0);
    assert!(refined.state().mole_fractions()[0][0] > 0.097);
    assert!(refined.state().mole_fractions()[0][0] < 0.099);
}

#[test]
fn original_caloric_binding_source_domains_and_bounded_thermal_work_are_enforced_atomically() {
    let net = network(900.0, 0.001 * R * 300.0);
    let initial = state(&net, [0.1, 0.9]);
    let foreign = network(901.0, 0.001 * R * 300.0);
    let wrong_pair = calorics(&foreign, [20.0, 30.0], [250.0, 600.0]);
    assert!(matches!(
        net.bind_caloric_field(&initial, &wrong_pair, || false),
        Err(Error::Invalid(_))
    ));
    let outside = calorics(&net, [20.0, 30.0], [350.0, 600.0]);
    assert!(
        net.bind_caloric_field(&initial, &outside, || false)
            .is_err()
    );
    let model = calorics(&net, [20.0, 30.0], [250.0, 600.0]);
    let bound = net.bind_caloric_field(&initial, &model, || false).unwrap();
    let cold_bracket = BinaryMiscibilityThermalSolve {
        temperature_bracket: [Temperature::new(250.0), Temperature::new(300.0)],
        ..thermal()
    };
    assert!(matches!(
        bound.adiabatic_candidate(0.2, solve(), cold_bracket, || false),
        Err(Error::HeatOutsideBracket { .. })
    ));
    let exhausted = BinaryMiscibilityThermalSolve {
        maximum_iterations: 1,
        ..thermal()
    };
    assert!(matches!(
        bound.adiabatic_candidate(0.2, solve(), exhausted, || false),
        Err(Error::IterationLimit)
    ));
    let invalid = BinaryMiscibilityThermalSolve {
        energy_tolerance_j: f64::NAN,
        ..thermal()
    };
    assert!(matches!(
        bound.adiabatic_candidate(0.2, solve(), invalid, || false),
        Err(Error::Invalid(_))
    ));
    let mut calls = 0;
    bound
        .adiabatic_candidate(0.2, solve(), thermal(), || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, calls] {
        let mut polls = 0;
        assert!(matches!(
            bound.adiabatic_candidate(0.2, solve(), thermal(), || {
                polls += 1;
                polls == stop
            }),
            Err(Error::Solver(BinaryMiscibilityError::Cancelled))
        ));
    }
    assert_eq!(bound.state(), &initial);
    assert_eq!(bound.caloric_model(), &model);
}

#[test]
fn refining_the_whole_thermal_step_converges_to_the_conservative_coupled_ode() {
    let kappa = 0.02 * R * 300.0;
    let x0 = 0.2;
    // Independent symmetric two-cell ODE: eliminate temperature using its
    // exact constant-Cp energy invariant, then integrate only composition.
    let temperature = |x: f64| {
        let change_h = 2.0 * R * 900.0 * (x - x0) * (1.0 - x - x0)
            + 0.5 * kappa * (2.0 * (x - x0)) * ((2.0 * x - 1.0) + (2.0 * x0 - 1.0));
        300.0 - change_h / 50.0
    };
    let rate = |x: f64| {
        let t = temperature(x);
        -0.5 * 0.4
            * ((x / (1.0 - x)).ln()
                + (900.0 / t) * (1.0 - 2.0 * x)
                + (kappa / (R * t)) * (2.0 * x - 1.0))
    };
    let mut reference = x0;
    let dt = 1.0 / 1024.0;
    for _ in 0..1024 {
        let a = rate(reference);
        let b = rate(reference + 0.5 * dt * a);
        let c = rate(reference + 0.5 * dt * b);
        let d = rate(reference + dt * c);
        reference += dt * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let mut errors = Vec::new();
    for intervals in [4, 8, 16] {
        let mut current = field(900.0, kappa, [x0, 1.0 - x0], [20.0, 30.0]);
        for _ in 0..intervals {
            current = current
                .adiabatic_candidate(1.0 / intervals as f64, solve(), thermal(), || false)
                .unwrap()
                .into_parts()
                .0;
        }
        let x = current.state().mole_fractions()[0][0];
        let t = current.state().model().temperature().value();
        near(current.state().elapsed_s(), 1.0, 2e-14);
        near(t, temperature(x), 1e-8);
        errors.push((x - reference).abs() + (t - temperature(reference)).abs() / 300.0);
    }
    assert!(errors[1] < errors[0] / 1.7, "{errors:?}");
    assert!(errors[2] < errors[1] / 1.7, "{errors:?}");
    assert!(errors[2] < 0.006, "{errors:?}");
}
