//! G0/G1/G3/G4 complete transport: independent dynamics and atomic accounts.

use super::*;
use crate::miscibility::advection::BinaryMiscibilityAdvectionFace;
use crate::miscibility::{BinaryMiscibilityCell, BinaryMiscibilityFace};
use fs_qty::{AmountConcentration, Temperature};
use fs_thermochem::solution::{BinaryRegularSolutionV1, gradient::BinarySquareGradientV1};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn fixture(
    chi: f64,
    kappa: f64,
    volumes: &[f64],
    connections: &[([usize; 2], f64)],
    rates: &[[f64; 2]],
) -> (BinaryMiscibilityNetwork, BinaryMiscibilityAdvection) {
    let pair = BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        chi,
        0.0,
        Temperature::new(250.0),
        Temperature::new(600.0),
    )
    .unwrap();
    let model = BinarySquareGradientV1::new(
        pair,
        Temperature::new(300.0),
        AmountConcentration::new(1.0),
        kappa * R * 300.0,
    )
    .unwrap();
    let network = BinaryMiscibilityNetwork::new(
        model,
        volumes
            .iter()
            .map(|&volume_m3| BinaryMiscibilityCell { volume_m3 })
            .collect(),
        connections
            .iter()
            .map(|&(cells, diffusivity_m2_s)| BinaryMiscibilityFace {
                cells,
                area_m2: 1.0,
                distance_m: 1.0,
                diffusivity_m2_s,
            })
            .collect(),
        1 << 20,
        || false,
    )
    .unwrap();
    let flow = network
        .prepare_advection(
            rates
                .iter()
                .map(
                    |&[forward_m3_s, reverse_m3_s]| BinaryMiscibilityAdvectionFace {
                        forward_m3_s,
                        reverse_m3_s,
                    },
                )
                .collect(),
            0.0,
            || false,
        )
        .unwrap();
    (network, flow)
}

fn state(net: &BinaryMiscibilityNetwork, x: &[f64]) -> BinaryMiscibilityState {
    net.state_from_mole_fractions(x.iter().map(|&x| [x, 1.0 - x]).collect(), || false)
        .unwrap()
}

fn time(tolerance: f64) -> BinaryMiscibilityTimeControl {
    BinaryMiscibilityTimeControl {
        initial_step_s: 0.5,
        minimum_step_s: 1e-9,
        maximum_step_s: 0.5,
        absolute_mole_fraction_tolerance: tolerance,
        relative_mole_fraction_tolerance: 0.0,
        maximum_attempts: 20_000,
    }
}

fn solve() -> BinaryMiscibilitySolve {
    BinaryMiscibilitySolve {
        absolute_transfer_tolerance_mol: 1e-17,
        relative_transfer_tolerance: 1e-12,
        amount_balance_tolerance_mol: 1e-11,
        free_energy_tolerance_j: 1e-8,
        ..BinaryMiscibilitySolve::default()
    }
}

fn advance(
    net: &BinaryMiscibilityNetwork,
    start: &BinaryMiscibilityState,
    flow: &BinaryMiscibilityAdvection,
    duration: f64,
    time: BinaryMiscibilityTimeControl,
) -> BinaryMiscibilityTransportAdvanceCandidate {
    net.advance_transport_candidate(
        start,
        flow,
        duration,
        BinaryMiscibilityAdvectionControls::default(),
        solve(),
        time,
        || false,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance {tolerance:e}"
    );
}

fn check_accounts(
    net: &BinaryMiscibilityNetwork,
    initial: &BinaryMiscibilityState,
    result: &BinaryMiscibilityTransportAdvanceCandidate,
) {
    let report = result.advance();
    let mut expected = initial.component_moles().to_vec();
    for (face, physical) in net.faces().enumerate() {
        for species in 0..2 {
            let transfer = report.transferred_component_moles[face][species];
            expected[physical.cells[0]][species] -= transfer;
            expected[physical.cells[1]][species] += transfer;
            near(
                transfer,
                report.advected_component_moles[face][species]
                    + if species == 0 {
                        report.diffused_a_mol[face]
                    } else {
                        -report.diffused_a_mol[face]
                    },
                2e-13,
            );
        }
    }
    for (cell, (expected, actual)) in expected
        .iter()
        .zip(result.state().component_moles())
        .enumerate()
    {
        for species in 0..2 {
            near(actual[species], expected[species], 3e-12);
        }
        near(actual[0] + actual[1], net.capacities_mol()[cell], 3e-12);
    }
    assert!(
        report
            .component_residual_mol
            .iter()
            .all(|x| x.abs() < 3e-12)
    );
    assert!(report.maximum_capacity_residual_mol < 3e-12);
    near(report.stage_energy_residual_j, 0.0, 2e-8);
    near(
        report.free_energy_change_j,
        report.advection_free_energy_change_j + report.diffusion_free_energy_change_j,
        2e-8,
    );
    assert!(
        report.diffusion_free_energy_change_j + report.diffusion_convex_split_dissipation_j
            <= solve().free_energy_tolerance_j
    );
    assert!(report.maximum_accepted_error_ratio <= 1.0);
    assert_eq!(
        report.completed_advection_substeps,
        report.completed_transport_steps as u64
    );
}

#[test]
fn coupled_refinement_converges_to_independent_three_cell_dynamics() {
    const VOLUMES: [f64; 3] = [1.0, 2.0, 0.8];
    const FACES: [([usize; 2], f64); 3] = [([0, 1], 0.35), ([0, 2], 0.2), ([1, 2], 0.45)];
    const RATES: [[f64; 2]; 3] = [[0.6, 0.0], [0.0, 0.6], [0.6, 0.0]];
    const CHI: f64 = 2.6;
    const KAPPA: f64 = 0.12;
    let (net, flow) = fixture(CHI, KAPPA, &VOLUMES, &FACES, &RATES);
    let initial = state(&net, &[0.12, 0.75, 0.43]);
    // Independent continuous-in-time FV equations. Neither chemical potentials
    // nor exchange rates are evaluated by the production thermodynamic solver.
    let derivative = |x: [f64; 3]| {
        let mut mu: [f64; 3] = x.map(|x| (x / (1.0 - x)).ln() + CHI * (1.0 - 2.0 * x));
        for &([i, j], _) in &FACES {
            let jump = KAPPA * (x[i] - x[j]);
            mu[i] += jump / VOLUMES[i];
            mu[j] -= jump / VOLUMES[j];
        }
        let mut derivative = [0.0; 3];
        for (&([i, j], diffusivity), &[forward, reverse]) in FACES.iter().zip(&RATES) {
            let mean = 0.5 * (x[i] + x[j]);
            let diffusive = diffusivity * mean * (1.0 - mean) * (mu[i] - mu[j]);
            let rate = forward * x[i] - reverse * x[j] + diffusive;
            derivative[i] -= rate / VOLUMES[i];
            derivative[j] += rate / VOLUMES[j];
        }
        derivative
    };
    let mut oracle = [0.12, 0.75, 0.43];
    let dt = 0.7 / 16_384.0;
    for _ in 0..16_384 {
        let a = derivative(oracle);
        let b = derivative(core::array::from_fn(|i| oracle[i] + 0.5 * dt * a[i]));
        let c = derivative(core::array::from_fn(|i| oracle[i] + 0.5 * dt * b[i]));
        let d = derivative(core::array::from_fn(|i| oracle[i] + dt * c[i]));
        for i in 0..3 {
            oracle[i] += dt / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]);
        }
    }
    let coarse = advance(&net, &initial, &flow, 0.7, time(2e-4));
    let fine = advance(&net, &initial, &flow, 0.7, time(2e-7));
    let error = |candidate: &BinaryMiscibilityTransportAdvanceCandidate| {
        candidate
            .state()
            .mole_fractions()
            .iter()
            .zip(oracle)
            .map(|(actual, expected)| (actual[0] - expected).abs())
            .fold(0.0_f64, f64::max)
    };
    let coarse_error = error(&coarse);
    let fine_error = error(&fine);
    assert!(
        fine_error < coarse_error / 15.0,
        "{fine_error:e} vs {coarse_error:e}"
    );
    assert!(fine_error < 1e-4, "{fine_error:e}");
    assert!(fine.advance().accepted_intervals > coarse.advance().accepted_intervals);
    assert!(fine.advance().rejected_intervals > 0);
    assert!(fine.advance().completed_root_iterations > 0);
    assert_eq!(fine.state().elapsed_s(), 0.7);
    check_accounts(&net, &initial, &fine);
}

#[test]
fn cfl_subcycling_cannot_hide_advection_time_error() {
    let (net, flow) = fixture(0.0, 0.0, &[1.0; 2], &[([0, 1], 0.0)], &[[1.0; 2]]);
    let initial = state(&net, &[0.9, 0.1]);
    // A naive full-versus-half driver sees IDENTICAL subcycled Euler paths at
    // dt=1 and reports zero, despite instantaneous artificial homogenization.
    let step = |start, dt| {
        net.transport_candidate(
            start,
            &flow,
            dt,
            BinaryMiscibilityAdvectionControls::default(),
            solve(),
            || false,
        )
        .unwrap()
    };
    let naive_full = step(&initial, 1.0);
    let naive_half = step(&initial, 0.5);
    let naive_fine = step(naive_half.state(), 0.5);
    assert_eq!(naive_full.state(), naive_fine.state());
    near(naive_full.state().mole_fractions()[0][0], 0.5, 1e-15);
    let mut control = time(2e-7);
    control.initial_step_s = 1.0;
    control.maximum_step_s = 1.0;
    let result = advance(&net, &initial, &flow, 1.0, control);
    // Exact zero-net-volume counterflow: x_left' = 1 - 2*x_left.
    let exact = 0.5 + 0.4 * (-2.0_f64).exp();
    near(result.state().mole_fractions()[0][0], exact, 1e-4);
    assert!(result.advance().accepted_intervals > 10);
    assert!(result.advance().maximum_accepted_step_s <= 0.5);
    assert!(result.advance().diffused_a_mol.iter().all(|&q| q == 0.0));
    assert_eq!(result.advance().diffusion_free_energy_change_j, 0.0);
    check_accounts(&net, &initial, &result);
}

#[test]
fn minimum_step_cannot_bypass_the_single_donor_cfl_bound() {
    let (net, flow) = fixture(0.0, 0.0, &[1.0; 2], &[([0, 1], 0.0)], &[[1.0; 2]]);
    let initial = state(&net, &[0.9, 0.1]);
    let original = initial.clone();
    let mut controls = time(1.0);
    controls.minimum_step_s = 0.75;
    controls.initial_step_s = 1.0;
    controls.maximum_step_s = 1.0;
    assert!(matches!(
        net.advance_transport_candidate(
            &initial,
            &flow,
            1.0,
            BinaryMiscibilityAdvectionControls::default(),
            solve(),
            controls,
            || false,
        ),
        Err(BinaryMiscibilityAdvanceError::Solver(
            Error::AdvectionWorkLimit {
                maximum_substeps: 1,
                ..
            }
        ))
    ));
    // A shorter final remainder is still a valid complete requested interval.
    let short = advance(&net, &initial, &flow, 0.2, controls);
    assert_eq!(short.state().elapsed_s(), 0.2);
    assert!(short.advance().suggested_next_step_s <= 0.5);
    assert_eq!(initial, original);
}

#[test]
fn zero_flow_retains_the_original_adaptive_diffusion_limit() {
    let (net, flow) = fixture(0.3, 0.02, &[1.0; 2], &[([0, 1], 0.7)], &[[0.0; 2]]);
    let initial = state(&net, &[0.6, 0.4]);
    let diffusion = net
        .advance_candidate(&initial, 0.8, time(1e-6), solve(), || false)
        .unwrap();
    let complete = advance(&net, &initial, &flow, 0.8, time(1e-6));
    for (actual, expected) in complete
        .state()
        .component_moles()
        .iter()
        .zip(diffusion.state().component_moles())
    {
        for species in 0..2 {
            near(actual[species], expected[species], 2e-12);
        }
    }
    assert_eq!(complete.advance().advection_free_energy_change_j, 0.0);
    assert!(
        complete
            .advance()
            .advected_component_moles
            .iter()
            .flatten()
            .all(|&q| q == 0.0)
    );
    near(
        complete.advance().diffused_a_mol[0],
        diffusion.advance().transferred_a_mol[0],
        2e-12,
    );
    check_accounts(&net, &initial, &complete);
}

#[test]
fn prescribed_stirring_can_raise_free_energy_while_diffusion_dissipates() {
    let (net, flow) = fixture(3.0, 0.001, &[1.0; 2], &[([0, 1], 0.01)], &[[1.0; 2]]);
    let initial = state(&net, &[0.07, 0.93]);
    let result = advance(&net, &initial, &flow, 0.2, time(1e-6));
    let report = result.advance();
    assert!(report.free_energy_change_j > 0.0);
    assert!(report.advection_free_energy_change_j > 0.0);
    assert!(report.diffusion_free_energy_change_j < 0.0);
    assert!(report.diffusion_convex_split_dissipation_j > 0.0);
    assert!(result.state().mole_fractions()[0][0] > 0.1);
    check_accounts(&net, &initial, &result);
}

#[test]
fn unequal_volumes_and_pure_donors_conserve_both_original_inventories() {
    let (net, flow) = fixture(0.0, 0.0, &[1.0, 2.5], &[([0, 1], 0.0)], &[[0.7; 2]]);
    let initial = state(&net, &[0.0, 1.0]);
    let result = advance(&net, &initial, &flow, 0.3, time(2e-7));
    let equilibrium = 2.5 / 3.5;
    let exact = equilibrium * (1.0 - (-0.7_f64 * 1.4 * 0.3).exp());
    near(result.state().mole_fractions()[0][0], exact, 1e-4);
    assert!(
        result
            .state()
            .mole_fractions()
            .iter()
            .flatten()
            .all(|&x| x > 0.0 && x < 1.0)
    );
    check_accounts(&net, &initial, &result);
}

#[test]
fn trace_second_species_controls_the_estimator_without_a_fraction_floor() {
    let (net, flow) = fixture(0.0, 0.0, &[1.0; 2], &[([0, 1], 0.0)], &[[0.4; 2]]);
    let initial = net
        .state_from_mole_fractions(vec![[1.0, 1e-18], [1.0, 2e-18]], || false)
        .unwrap();
    let mut controls = time(1e-30);
    controls.relative_mole_fraction_tolerance = 1e-5;
    let result = advance(&net, &initial, &flow, 0.5, controls);
    let exact = 1.5e-18 - 0.5e-18 * (-0.4_f64).exp();
    near(result.state().component_moles()[0][1], exact, 2e-21);
    assert!(result.advance().accepted_intervals > 1);
    near(
        result.state().component_moles().iter().map(|x| x[1]).sum(),
        3e-18,
        1e-32,
    );
    assert!(result.state().component_moles().iter().all(|x| x[1] > 0.0));
}

#[test]
fn work_accuracy_input_and_late_cancellation_refusals_are_atomic() {
    let (net, flow) = fixture(1.0, 0.1, &[1.0; 2], &[([0, 1], 0.3)], &[[0.2; 2]]);
    let initial = state(&net, &[0.85, 0.15]);
    let original = initial.clone();
    let run = |advection, time, cancelled: &mut dyn FnMut() -> bool| {
        net.advance_transport_candidate(&initial, &flow, 0.8, advection, solve(), time, cancelled)
    };
    let mut control = time(1e-12);
    control.maximum_attempts = 1;
    assert!(matches!(
        run(
            BinaryMiscibilityAdvectionControls::default(),
            control,
            &mut || false
        ),
        Err(BinaryMiscibilityAdvanceError::WorkLimit { .. })
    ));
    control.maximum_attempts = 100;
    control.minimum_step_s = 0.5;
    assert!(matches!(
        run(
            BinaryMiscibilityAdvectionControls::default(),
            control,
            &mut || false
        ),
        Err(BinaryMiscibilityAdvanceError::TemporalAccuracy { .. })
    ));
    let invalid = BinaryMiscibilityAdvectionControls {
        maximum_substeps: 0,
        ..Default::default()
    };
    assert!(matches!(
        run(invalid, time(1e-5), &mut || false),
        Err(BinaryMiscibilityAdvanceError::Solver(Error::Invalid(_)))
    ));
    let mut polls = 0;
    run(
        BinaryMiscibilityAdvectionControls::default(),
        time(1e-5),
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut cancelled_polls = 0;
    assert_eq!(
        run(
            BinaryMiscibilityAdvectionControls::default(),
            time(1e-5),
            &mut || {
                cancelled_polls += 1;
                cancelled_polls == polls
            }
        ),
        Err(BinaryMiscibilityAdvanceError::Solver(Error::Cancelled))
    );
    assert_eq!(initial, original);
}

#[test]
fn continuation_and_a_short_final_remainder_keep_one_physical_clock() {
    let (net, flow) = fixture(0.0, 0.0, &[1.0; 2], &[([0, 1], 0.2)], &[[0.1; 2]]);
    let initial = state(&net, &[0.6, 0.4]);
    let first = advance(&net, &initial, &flow, 0.123, time(1e-5));
    let mut controls = time(1e-5);
    controls.minimum_step_s = 0.001;
    controls.initial_step_s = 0.017;
    controls.maximum_step_s = 0.017;
    let next = advance(&net, first.state(), &flow, 0.1001, controls);
    assert_eq!(next.state().elapsed_s(), 0.123_f64 + 0.1001);
    assert_eq!(
        next.advance().duration_s,
        next.state().elapsed_s() - first.state().elapsed_s()
    );
    let stationary = state(&net, &[0.5, 0.5]);
    let mut controls = time(1e-5);
    controls.minimum_step_s = 0.01;
    let short = advance(&net, &stationary, &flow, 0.5001, controls);
    assert!(short.advance().minimum_accepted_step_s < controls.minimum_step_s);
    assert_eq!(short.state().elapsed_s(), 0.5001);
    assert_eq!(
        short.state().component_moles(),
        stationary.component_moles()
    );
    assert_eq!(short.advance().free_energy_change_j, 0.0);
}
