use super::*;
use crate::miscibility::{BinaryMiscibilityCell, BinaryMiscibilityFace};
use fs_qty::{AmountConcentration, Temperature};
use fs_thermochem::solution::{BinaryRegularSolutionV1, gradient::BinarySquareGradientV1};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn fixture(chi: f64, kappa: f64) -> (BinaryMiscibilityNetwork, BinaryMiscibilityState) {
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
        vec![BinaryMiscibilityCell { volume_m3: 1.0 }; 2],
        vec![BinaryMiscibilityFace {
            cells: [0, 1],
            area_m2: 1.0,
            distance_m: 1.0,
            diffusivity_m2_s: 1.0,
        }],
        1 << 20,
        || false,
    )
    .unwrap();
    let state = network
        .state_from_mole_fractions(vec![[0.6, 0.4], [0.4, 0.6]], || false)
        .unwrap();
    (network, state)
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

#[test]
fn automatic_refinement_improves_physical_time_accuracy_and_retains_flux_balance() {
    let (net, initial) = fixture(0.0, 0.0);
    // Independent exact two-cell ideal-mixture solution. With arithmetic face
    // mobility and x_left+x_right=1, y=2*x_left-1 obeys y'=-2*atanh(y).
    // Integrating that scalar ODE with RK4 avoids reusing the implicit scheme.
    let mut y = 0.2_f64;
    let dt = 1.0 / 8192.0;
    for _ in 0..8192 {
        let derivative = |y: f64| -2.0 * y.atanh();
        let a = derivative(y);
        let b = derivative(y + 0.5 * dt * a);
        let c = derivative(y + 0.5 * dt * b);
        let d = derivative(y + dt * c);
        y += dt / 6.0 * (a + 2.0 * b + 2.0 * c + d);
    }
    let oracle = 0.5 * (1.0 + y);
    let coarse = net
        .advance_candidate(&initial, 1.0, time(2e-4), solve(), || false)
        .unwrap();
    let fine = net
        .advance_candidate(&initial, 1.0, time(2e-7), solve(), || false)
        .unwrap();
    let coarse_error = (coarse.state().mole_fractions()[0][0] - oracle).abs();
    let fine_error = (fine.state().mole_fractions()[0][0] - oracle).abs();
    assert!(
        fine_error < coarse_error / 15.0,
        "{fine_error:e} >= {coarse_error:e}/15"
    );
    assert!(fine_error < 8e-5);
    let report = fine.advance();
    assert!(report.rejected_intervals > 0);
    assert!(report.accepted_intervals > coarse.advance().accepted_intervals);
    assert!(report.maximum_accepted_error_ratio <= 1.0);
    assert_eq!(fine.state().elapsed_s(), 1.0);
    assert!(report.free_energy_change_j < 0.0);
    assert!(
        report.free_energy_change_j + report.convex_split_dissipation_j
            <= solve().free_energy_tolerance_j
    );
    for species in 0..2 {
        let sum: f64 = fine
            .state()
            .component_moles()
            .iter()
            .map(|n| n[species])
            .sum();
        assert!((sum - 1.0).abs() < 1e-12);
        assert!(report.component_residual_mol[species].abs() < 1e-12);
    }
    let actual_transfer = initial.component_moles()[0][0] - fine.state().component_moles()[0][0];
    assert!((report.transferred_a_mol[0] - actual_transfer).abs() < 1e-12);
    assert_eq!(initial.mole_fractions(), &[[0.6, 0.4], [0.4, 0.6]]);
}

#[test]
fn adaptive_spinodal_evolution_remains_homogeneous_local_and_dissipative() {
    let (net, initial) = fixture(3.0, 0.1);
    let candidate = net
        .advance_candidate(&initial, 2.0, time(1e-6), solve(), || false)
        .unwrap();
    let final_x = candidate.state().mole_fractions()[0][0];
    assert!(
        final_x > 0.7 && final_x < 1.0,
        "spinodal amplitude did not grow: {final_x}"
    );
    assert!(candidate.advance().free_energy_change_j < 0.0);
    assert!(candidate.advance().convex_split_dissipation_j > 0.0);
    assert!(
        candidate
            .state()
            .mole_fractions()
            .iter()
            .flatten()
            .all(|&x| x > 0.0 && x < 1.0)
    );
    assert!(candidate.advance().maximum_capacity_residual_mol < 1e-12);
}

#[test]
fn work_accuracy_and_cancellation_refusals_never_publish_a_partial_window() {
    let (net, initial) = fixture(0.0, 0.0);
    let original = initial.clone();
    let mut controls = time(1e-12);
    controls.maximum_attempts = 1;
    assert!(matches!(
        net.advance_candidate(&initial, 1.0, controls, solve(), || false),
        Err(BinaryMiscibilityAdvanceError::WorkLimit { .. })
    ));
    controls.maximum_attempts = 100;
    controls.minimum_step_s = 0.5;
    assert!(matches!(
        net.advance_candidate(&initial, 1.0, controls, solve(), || false),
        Err(BinaryMiscibilityAdvanceError::TemporalAccuracy { .. })
    ));
    let mut polls = 0;
    assert_eq!(
        net.advance_candidate(&initial, 1.0, time(1e-7), solve(), || {
            polls += 1;
            polls > 100
        }),
        Err(BinaryMiscibilityAdvanceError::Solver(Error::Cancelled))
    );
    assert_eq!(initial, original);
}

#[test]
fn continuation_lands_on_requested_time_with_a_short_final_remainder() {
    let (net, initial) = fixture(0.0, 0.0);
    let first = net
        .advance_candidate(&initial, 0.123, time(1e-5), solve(), || false)
        .unwrap();
    let mut controls = time(1e-5);
    controls.minimum_step_s = 0.001;
    controls.initial_step_s = 0.017;
    controls.maximum_step_s = 0.017;
    let next = net
        .advance_candidate(first.state(), 0.1001, controls, solve(), || false)
        .unwrap();
    assert_eq!(next.state().elapsed_s(), 0.123_f64 + 0.1001);
    assert_eq!(
        next.advance().duration_s,
        next.state().elapsed_s() - first.state().elapsed_s()
    );
    let stationary = net
        .state_from_mole_fractions(vec![[0.5, 0.5]; 2], || false)
        .unwrap();
    controls.initial_step_s = 0.5;
    controls.maximum_step_s = 0.5;
    let short = net
        .advance_candidate(&stationary, 0.5001, controls, solve(), || false)
        .unwrap();
    assert!(short.advance().minimum_accepted_step_s < controls.minimum_step_s);
    assert_eq!(
        short.state().component_moles(),
        stationary.component_moles()
    );
    assert_eq!(short.advance().free_energy_change_j, 0.0);
}
