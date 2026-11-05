//! Prescribed quenches preserve field history and change the actual dynamics.

use super::caloric::BinaryMiscibilityCaloricError;
use crate::miscibility::advection::{
    BinaryMiscibilityAdvection, BinaryMiscibilityAdvectionControls, BinaryMiscibilityAdvectionFace,
};
use crate::miscibility::{
    BinaryMiscibilityCell, BinaryMiscibilityError, BinaryMiscibilityFace, BinaryMiscibilityNetwork,
    BinaryMiscibilitySolve, BinaryMiscibilityState,
};
use fs_qty::{AmountConcentration, Temperature};
use fs_thermochem::solution::caloric::{BinaryLiquidCaloricModelV1, PureLiquidHeatCapacityV1};
use fs_thermochem::solution::{BinaryRegularSolutionV1, gradient::BinarySquareGradientV1};
use fs_thermochem::{MolarThermalQuantityV1, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn network(
    a: f64,
    b: f64,
    temperature: f64,
    kappa: f64,
    volumes: &[f64],
) -> BinaryMiscibilityNetwork {
    let pair = BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        a,
        b,
        Temperature::new(250.0),
        Temperature::new(600.0),
    )
    .unwrap();
    let model = BinarySquareGradientV1::new(
        pair,
        Temperature::new(temperature),
        AmountConcentration::new(1.0),
        kappa,
    )
    .unwrap();
    BinaryMiscibilityNetwork::new(
        model,
        volumes
            .iter()
            .map(|&volume_m3| BinaryMiscibilityCell { volume_m3 })
            .collect(),
        (0..volumes.len() - 1)
            .map(|i| BinaryMiscibilityFace {
                cells: [i, i + 1],
                area_m2: 1.0,
                distance_m: 1.0,
                diffusivity_m2_s: 0.4,
            })
            .collect(),
        16 * 1024 * 1024,
        || false,
    )
    .unwrap()
}

fn state(net: &BinaryMiscibilityNetwork, fractions: &[f64]) -> BinaryMiscibilityState {
    net.state_from_mole_fractions(fractions.iter().map(|&x| [x, 1.0 - x]).collect(), || false)
        .unwrap()
}

fn counterflow(net: &BinaryMiscibilityNetwork, rate: f64) -> BinaryMiscibilityAdvection {
    net.prepare_advection(
        vec![BinaryMiscibilityAdvectionFace {
            forward_m3_s: rate,
            reverse_m3_s: rate,
        }],
        0.0,
        || false,
    )
    .unwrap()
}

fn strict() -> BinaryMiscibilitySolve {
    BinaryMiscibilitySolve {
        absolute_transfer_tolerance_mol: 1e-18,
        relative_transfer_tolerance: 1e-12,
        ..BinaryMiscibilitySolve::default()
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance {tolerance:e}"
    );
}

fn gradient_energy(net: &BinaryMiscibilityNetwork, state: &BinaryMiscibilityState) -> f64 {
    net.faces()
        .map(|face| {
            let [i, j] = face.cells;
            let difference = state.mole_fractions()[i][0] - state.mole_fractions()[j][0];
            0.5 * net.model().gradient_coefficient_j_m() * face.area_m2 / face.distance_m
                * difference
                * difference
        })
        .sum()
}

#[test]
fn quench_preserves_clock_interfaces_and_carried_sub_ulp_inventory() {
    let net = network(0.0, 900.0, 500.0, 0.01 * R * 300.0, &[0.7, 1.3]);
    let initial = state(&net, &[0.2, 0.8]);
    let flow = counterflow(&net, 0.1);
    let (advanced, _) = net
        .advect_candidate(
            &initial,
            &flow,
            1e-16,
            BinaryMiscibilityAdvectionControls::default(),
            || false,
        )
        .unwrap()
        .into_parts();
    assert!(
        advanced
            .amount_remainders
            .iter()
            .flatten()
            .any(|&n| n != 0.0)
    );
    assert!(advanced.elapsed_s() > 0.0);
    let original = advanced.clone();
    let candidate = net
        .temperature_candidate(&advanced, Temperature::new(300.0), || false)
        .unwrap();
    assert_eq!(
        candidate.state().component_moles(),
        advanced.component_moles()
    );
    assert_eq!(
        candidate.state().amount_remainders,
        advanced.amount_remainders
    );
    assert_eq!(
        candidate.state().mole_fractions(),
        advanced.mole_fractions()
    );
    assert_eq!(candidate.state().elapsed_s(), advanced.elapsed_s());
    assert_eq!(candidate.network().cells(), net.cells());
    assert_eq!(candidate.network().capacities_mol(), net.capacities_mol());
    assert_eq!(
        candidate.state().model().thermodynamics(),
        net.model().thermodynamics()
    );
    assert_eq!(
        candidate.state().model().gradient_coefficient_j_m(),
        net.model().gradient_coefficient_j_m()
    );
    assert_eq!(
        candidate.state().model().molar_concentration_mol_m3(),
        net.model().molar_concentration_mol_m3()
    );
    assert_eq!(
        gradient_energy(candidate.network(), candidate.state()),
        gradient_energy(&net, &advanced)
    );
    near(candidate.change().chi_before, 1.8, 1e-15);
    near(candidate.change().chi_after, 3.0, 1e-15);
    assert_eq!(advanced, original);
}

#[test]
fn cooling_across_the_spinodal_changes_relaxation_to_demixing_without_a_flash() {
    let warm = network(0.0, 900.0, 500.0, 0.01 * R * 300.0, &[1.0; 2]);
    let initial = state(&warm, &[0.4999, 0.5001]);
    let cold = warm
        .temperature_candidate(&initial, Temperature::new(300.0), || false)
        .unwrap();
    assert_eq!(cold.state().mole_fractions(), initial.mole_fractions());
    let dt = 1e-3;
    let warm_step = warm
        .step_candidate(&initial, dt, strict(), || false)
        .unwrap();
    let cold_step = cold
        .network()
        .step_candidate(cold.state(), dt, strict(), || false)
        .unwrap();
    let amplitude = |s: &BinaryMiscibilityState| 0.5 - s.mole_fractions()[0][0];
    assert!(amplitude(warm_step.state()) < amplitude(&initial));
    assert!(amplitude(cold_step.state()) > amplitude(&initial));
    // Two-cell graph eigenvalue is 2. At x=1/2, the degenerate
    // mobility is D/4 and the gradient coefficient is kappa/(c R T).
    for (net, result) in [(&warm, &warm_step), (cold.network(), &cold_step)] {
        let expected_rate = -0.4
            * 0.25
            * 2.0
            * (4.0 - 2.0 * net.model().chi()
                + 2.0 * net.model().gradient_coefficient_j_m() / net.model().rt_j_per_mol());
        let measured_rate = (amplitude(result.state()) / amplitude(&initial) - 1.0) / dt;
        near(measured_rate, expected_rate, 0.005 * expected_rate.abs());
        assert!(result.step().free_energy_change_j <= strict().free_energy_tolerance_j);
    }
}

#[test]
fn ideal_diffusion_is_temperature_invariant_and_states_require_the_new_owner() {
    let net = network(0.0, 0.0, 300.0, 0.0, &[0.7, 1.3]);
    let initial = state(&net, &[0.17, 0.86]);
    let hot = net
        .temperature_candidate(&initial, Temperature::new(550.0), || false)
        .unwrap();
    let old_step = net
        .step_candidate(&initial, 0.3, strict(), || false)
        .unwrap();
    let hot_step = hot
        .network()
        .step_candidate(hot.state(), 0.3, strict(), || false)
        .unwrap();
    for (old, new) in old_step
        .state()
        .component_moles()
        .iter()
        .zip(hot_step.state().component_moles())
    {
        for species in 0..2 {
            near(old[species], new[species], 2e-13);
        }
    }
    near(
        old_step.step().transferred_a_mol[0],
        hot_step.step().transferred_a_mol[0],
        2e-13,
    );
    assert!(matches!(
        net.step_candidate(hot.state(), 0.3, strict(), || false),
        Err(BinaryMiscibilityError::Invalid(_))
    ));
    assert!(matches!(
        hot.network()
            .step_candidate(&initial, 0.3, strict(), || false),
        Err(BinaryMiscibilityError::Invalid(_))
    ));
}

#[test]
fn unsupported_quenches_and_early_or_final_cancellation_are_atomic() {
    let net = network(0.0, 900.0, 500.0, 0.0, &[1.0; 2]);
    let initial = state(&net, &[0.2, 0.8]);
    let original = initial.clone();
    // A stable model with kappa=0 cannot be quenched into an ill-posed
    // spatial spinodal; unsupported temperatures are never clamped.
    for temperature in [300.0, 249.0, 601.0, f64::NAN] {
        assert!(
            net.temperature_candidate(&initial, Temperature::new(temperature), || false)
                .is_err()
        );
    }
    let mut calls = 0;
    net.temperature_candidate(&initial, Temperature::new(550.0), || {
        calls += 1;
        false
    })
    .unwrap();
    for stop in [1, calls] {
        let mut polls = 0;
        assert!(matches!(
            net.temperature_candidate(&initial, Temperature::new(550.0), || {
                polls += 1;
                polls == stop
            }),
            Err(BinaryMiscibilityError::Cancelled)
        ));
    }
    assert_eq!(initial, original);
}

#[test]
fn repeated_temperature_round_trips_recover_the_original_mobility_without_drift() {
    let original = network(0.5, 500.0, 300.0, 0.01 * R * 300.0, &[0.7, 1.3]);
    let initial = state(&original, &[0.13, 0.82]);
    let unchanged = original
        .temperature_candidate(&initial, Temperature::new(300.0), || false)
        .unwrap();
    assert_eq!(unchanged.network(), &original);
    assert_eq!(unchanged.state(), &initial);
    assert_eq!(unchanged.change().external_free_energy_change_j, 0.0);
    let (mut net, mut current, _) = unchanged.into_parts();
    for _ in 0..4 {
        for temperature in [250.0, 600.0, 310.0, 300.0] {
            (net, current, _) = net
                .temperature_candidate(&current, Temperature::new(temperature), || false)
                .unwrap()
                .into_parts();
        }
        assert_eq!(net, original);
        assert_eq!(current, initial);
    }
}

#[test]
fn direct_temperature_work_retains_trace_entropy_beside_an_unchanged_large_interface() {
    let net = network(0.0, 0.0, 300.0, 1e25, &[1.0; 2]);
    let traces = [1e-30, 3e-30];
    let initial = net
        .state_from_mole_fractions(vec![[1.0, traces[0]], [traces[1], 1.0]], || false)
        .unwrap();
    let next_temperature = f64::from_bits(300.0_f64.to_bits() + 1);
    let candidate = net
        .temperature_candidate(&initial, Temperature::new(next_temperature), || false)
        .unwrap();
    let expected_entropy: f64 = traces
        .into_iter()
        .map(|y| -R * ((-y).ln_1p() + y * y.ln()))
        .sum();
    let expected_change = -(next_temperature - 300.0) * expected_entropy;
    let change = candidate.change();
    // The bulk increment is smaller than one ulp of the unchanged interface.
    assert_eq!(change.free_energy_before_j, change.free_energy_after_j);
    assert_ne!(change.external_free_energy_change_j, 0.0);
    near(
        change.mixing_entropy_j_per_k,
        expected_entropy,
        expected_entropy * 3e-14,
    );
    near(
        change.external_free_energy_change_j,
        expected_change,
        expected_change.abs() * 3e-14,
    );
    assert_eq!(candidate.state().mole_fractions(), initial.mole_fractions());
}

#[test]
fn explicit_flow_rebinding_retains_transport_and_rejects_a_foreign_geometry() {
    let net = network(0.0, 0.0, 300.0, 0.0, &[0.7, 1.3]);
    let initial = state(&net, &[0.2, 0.8]);
    let flow = counterflow(&net, 0.3);
    let hot = net
        .temperature_candidate(&initial, Temperature::new(450.0), || false)
        .unwrap();
    let controls = BinaryMiscibilityAdvectionControls::default();
    assert!(matches!(
        hot.network()
            .advect_candidate(hot.state(), &flow, 0.5, controls, || false),
        Err(BinaryMiscibilityError::Invalid(_))
    ));
    let rebound = hot.network().rebind_advection(&flow, || false).unwrap();
    assert_eq!(rebound.face_rates(), flow.face_rates());
    assert_eq!(
        rebound.maximum_outgoing_rate_s(),
        flow.maximum_outgoing_rate_s()
    );
    let original = net
        .advect_candidate(&initial, &flow, 0.5, controls, || false)
        .unwrap();
    let result = hot
        .network()
        .advect_candidate(hot.state(), &rebound, 0.5, controls, || false)
        .unwrap();
    assert_eq!(
        result.state().component_moles(),
        original.state().component_moles()
    );
    assert_eq!(
        result.state().amount_remainders,
        original.state().amount_remainders
    );
    assert_eq!(
        result.step().transferred_component_moles,
        original.step().transferred_component_moles
    );
    let foreign = network(0.0, 0.0, 450.0, 0.0, &[0.9, 1.1]);
    assert!(matches!(
        foreign.rebind_advection(&flow, || false),
        Err(BinaryMiscibilityError::Invalid(_))
    ));
    assert!(matches!(
        hot.network().rebind_advection(&flow, || true),
        Err(BinaryMiscibilityError::Cancelled)
    ));
}

fn calorics(
    net: &BinaryMiscibilityNetwork,
    reference: f64,
    upper: f64,
    intercept: [f64; 2],
    slope: [f64; 2],
) -> BinaryLiquidCaloricModelV1 {
    let pair = net.model().thermodynamics();
    let capacities = [0, 1].map(|species| {
        PureLiquidHeatCapacityV1::new(
            pair.species()[species].clone(),
            Temperature::new(reference),
            [250.0, upper]
                .map(|t| {
                    (
                        Temperature::new(t),
                        MolarThermalQuantityV1::new(intercept[species] + slope[species] * t),
                    )
                })
                .to_vec(),
        )
        .unwrap()
    });
    BinaryLiquidCaloricModelV1::new(pair.clone(), capacities).unwrap()
}

#[test]
fn caloric_receipt_uses_original_cp_integrals_without_flashing_the_quenched_field() {
    let net = network(0.2, 900.0, 500.0, 0.01 * R * 300.0, &[0.7, 1.3]);
    let initial = state(&net, &[0.15, 0.8]);
    let candidate = net
        .temperature_candidate(&initial, Temperature::new(300.0), || false)
        .unwrap();
    let model = calorics(&net, 250.0, 600.0, [10.0, 40.0], [0.02, 0.01]);
    let receipt = candidate.caloric_change(&model, || false).unwrap();
    let amounts: [f64; 2] =
        [0, 1].map(|species| initial.component_moles().iter().map(|n| n[species]).sum());
    let expected_heat = amounts[0] * -3600.0 + amounts[1] * -8800.0;
    let expected_entropy =
        amounts[0] * (10.0 * 0.6_f64.ln() - 4.0) + amounts[1] * (40.0 * 0.6_f64.ln() - 2.0);
    near(receipt.heat_input_j, expected_heat, 2e-11);
    near(receipt.entropy_change_j_per_k, expected_entropy, 2e-13);
    assert!(receipt.heat_input_j < 0.0);
    assert!(candidate.change().external_free_energy_change_j > 0.0);
    assert_eq!(receipt.temperature_before, Temperature::new(500.0));
    assert_eq!(receipt.temperature_after, Temperature::new(300.0));
    assert_eq!(candidate.state().mole_fractions(), initial.mole_fractions());
    let other_datum = calorics(&net, 600.0, 600.0, [10.0, 40.0], [0.02, 0.01]);
    assert_eq!(
        candidate.caloric_change(&other_datum, || false).unwrap(),
        receipt
    );
    let reverse = candidate
        .network()
        .temperature_candidate(candidate.state(), Temperature::new(500.0), || false)
        .unwrap();
    let reverse_receipt = reverse.caloric_change(&model, || false).unwrap();
    assert_eq!(reverse_receipt.heat_input_j, -receipt.heat_input_j);
    assert_eq!(
        reverse_receipt.entropy_change_j_per_k,
        -receipt.entropy_change_j_per_k
    );
}

#[test]
fn caloric_receipt_retains_independent_trace_inventory_and_a_tiny_temperature_span() {
    let net = network(0.0, 0.0, 300.0, 0.0, &[1.0; 2]);
    let initial = net
        .state_from_mole_fractions(vec![[1.0, 1e-30], [1.0, 3e-30]], || false)
        .unwrap();
    let temperature = f64::from_bits(300.0_f64.to_bits() + 1);
    let candidate = net
        .temperature_candidate(&initial, Temperature::new(temperature), || false)
        .unwrap();
    let model = calorics(&net, 250.0, 600.0, [5.0, 1e30], [0.0; 2]);
    let receipt = candidate.caloric_change(&model, || false).unwrap();
    // The independently represented trace contributes 4 J/K despite the
    // dominant fractions rounding to one. Reconstructing B as 1-A loses it.
    let expected_heat = 14.0 * (temperature - 300.0);
    let expected_entropy = 14.0 * ((temperature - 300.0) / 300.0).ln_1p();
    near(receipt.heat_input_j, expected_heat, expected_heat * 3e-14);
    near(
        receipt.entropy_change_j_per_k,
        expected_entropy,
        expected_entropy * 3e-14,
    );
    assert!(receipt.heat_input_j > 0.0);
}

#[test]
fn caloric_receipt_refuses_changed_pairs_source_loss_and_late_cancellation() {
    let net = network(0.0, 0.0, 300.0, 0.0, &[1.0; 2]);
    let initial = state(&net, &[0.2, 0.8]);
    let candidate = net
        .temperature_candidate(&initial, Temperature::new(500.0), || false)
        .unwrap();
    let original = candidate.state().clone();
    let foreign = network(0.1, 0.0, 300.0, 0.0, &[1.0; 2]);
    let wrong_pair = calorics(&foreign, 250.0, 600.0, [20.0, 30.0], [0.0; 2]);
    assert!(matches!(
        candidate.caloric_change(&wrong_pair, || false),
        Err(BinaryMiscibilityCaloricError::Invalid(_))
    ));
    let unsupported = calorics(&net, 250.0, 400.0, [20.0, 30.0], [0.0; 2]);
    assert!(matches!(
        candidate.caloric_change(&unsupported, || false),
        Err(BinaryMiscibilityCaloricError::Caloric(_))
    ));
    let valid = calorics(&net, 250.0, 600.0, [20.0, 30.0], [0.0; 2]);
    let mut calls = 0;
    candidate
        .caloric_change(&valid, || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, calls] {
        let mut polls = 0;
        assert!(matches!(
            candidate.caloric_change(&valid, || {
                polls += 1;
                polls == stop
            }),
            Err(BinaryMiscibilityCaloricError::Solver(
                BinaryMiscibilityError::Cancelled
            ))
        ));
    }
    assert_eq!(candidate.state(), &original);
}
