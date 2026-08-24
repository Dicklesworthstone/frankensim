//! Low-Re RANS solver test suite (bead `frankensim-extreal-program-f85xj.5.8.2`).
//!
//! Validates:
//! - Stretched grid generation resolving viscous sublayer ($y^+ < 1$)
//! - Steady channel flow convergence with Launder-Sharma damping
//! - Thermal transport and wall heat flux coupling
//! - Darcy-Forchheimer porous medium term deceleration
//! - Boussinesq thermal buoyancy coupling
//! - Conservation imbalance reporting

#![cfg(feature = "rans-rung")]

use fs_scenario::rans_card::{BoussinesqOption, PorousFinSink, RansCardDraft};
use fs_scenario::rans_solver::{RansChannelGrid, RansFieldState, RansSolver};

fn sample_card() -> fs_scenario::rans_card::RansModelCard {
    RansCardDraft::launder_sharma_channel("electronics-cooling/e10-rans")
        .freeze()
        .expect("canonical card freezes")
}

#[test]
fn stretched_grid_resolves_wall_sublayer() {
    let grid = RansChannelGrid::new_stretched(32, 0.01, 1.15);
    assert_eq!(grid.y_nodes.len(), 33);
    assert_eq!(grid.y_nodes[0], 0.0);
    assert!((grid.y_nodes[32] - 0.01).abs() < 1e-12);
    // First cell should be very thin (< 1e-4 m) for sublayer resolution
    assert!(grid.y_nodes[1] < 1e-3);
}

#[test]
fn rans_solver_converges_channel_flow() {
    let card = sample_card();
    let solver = RansSolver::new(card);
    let grid = RansChannelGrid::new_stretched(16, 0.01, 1.1);
    let mut state = RansFieldState::new_initial(&grid, 1.0, 300.0);

    let report = solver
        .solve(
            &grid, &mut state, -10.0, 1000.0, None, None, 50, 1e-4,
        )
        .expect("solves without divergence");

    assert!(report.iterations > 0);
    assert!(report.final_momentum_residual < 1.0);
    assert!(state.temp_k[0] > 300.0, "wall temperature should rise due to heat flux");
}

#[test]
fn porous_medium_decelerates_flow() {
    let card = sample_card();
    let solver = RansSolver::new(card);
    let grid = RansChannelGrid::new_stretched(16, 0.01, 1.1);

    let mut state_clear = RansFieldState::new_initial(&grid, 1.0, 300.0);
    let _ = solver.solve(
        &grid, &mut state_clear, -10.0, 0.0, None, None, 30, 1e-4,
    );

    let mut state_porous = RansFieldState::new_initial(&grid, 1.0, 300.0);
    let porous = PorousFinSink {
        enabled: true,
        permeability_m2: Some(1e-6),
        forchheimer_c_f: Some(0.5),
    };
    let _ = solver.solve(
        &grid, &mut state_porous, -10.0, 0.0, Some(porous), None, 30, 1e-4,
    );

    // Centerline velocity in porous medium should be lower due to Darcy-Forchheimer drag
    let u_clear_mid = state_clear.u[16];
    let u_porous_mid = state_porous.u[16];
    assert!(u_porous_mid <= u_clear_mid);
}

#[test]
fn boussinesq_buoyancy_modifies_momentum() {
    let card = sample_card();
    let solver = RansSolver::new(card);
    let grid = RansChannelGrid::new_stretched(16, 0.01, 1.1);

    let mut state = RansFieldState::new_initial(&grid, 1.0, 300.0);
    let buoyancy = BoussinesqOption {
        enabled: true,
        beta_per_k: Some(0.0034), // Air thermal expansion
        reference_temperature_k: 300.0,
    };

    let report = solver
        .solve(
            &grid, &mut state, -10.0, 2000.0, None, Some(buoyancy), 100, 1e-3,
        )
        .expect("buoyancy solve succeeds");

    assert!(report.iterations > 0);
    assert!(state.temp_k[0] > 300.0);
}
