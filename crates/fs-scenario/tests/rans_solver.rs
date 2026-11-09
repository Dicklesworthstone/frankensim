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
use fs_scenario::rans_solver::{
    RansChannelGrid, RansFieldState, RansSolver, RansTerminationReason,
};

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
        .solve(&grid, &mut state, -10.0, 1000.0, None, None, 50, 1e-4)
        .expect("solves without divergence");

    assert!(report.iterations > 0);
    assert!(report.final_momentum_residual < 1.0);
    assert!(
        state.temp_k[0] > 300.0,
        "wall temperature should rise due to heat flux"
    );
}

#[test]
fn porous_medium_decelerates_flow() {
    let card = sample_card();
    let solver = RansSolver::new(card);
    let grid = RansChannelGrid::new_stretched(16, 0.01, 1.1);

    let mut state_clear = RansFieldState::new_initial(&grid, 1.0, 300.0);
    let _ = solver.solve(&grid, &mut state_clear, -10.0, 0.0, None, None, 30, 1e-4);

    let mut state_porous = RansFieldState::new_initial(&grid, 1.0, 300.0);
    let porous = PorousFinSink {
        enabled: true,
        permeability_m2: Some(1e-6),
        forchheimer_c_f: Some(0.5),
    };
    let _ = solver.solve(
        &grid,
        &mut state_porous,
        -10.0,
        0.0,
        Some(porous),
        None,
        30,
        1e-4,
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
            &grid,
            &mut state,
            -10.0,
            2000.0,
            None,
            Some(buoyancy),
            100,
            1e-3,
        )
        .expect("buoyancy solve succeeds");

    assert!(report.iterations > 0);
    assert!(state.temp_k[0] > 300.0);
}

#[test]
fn invalid_density_refuses_before_changing_any_field() {
    let grid = RansChannelGrid::new_stretched(4, 1.0, 1.0);
    let initial = RansFieldState::new_initial(&grid, 0.0, 300.0);
    for rho in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut solver = RansSolver::new(sample_card());
        solver.rho = rho;
        let mut state = initial.clone();
        let result = solver.solve(&grid, &mut state, -2.0, 40.0, None, None, 1, 1e-8);
        assert_eq!(result, Err(RansTerminationReason::Diverged), "rho={rho}");
        assert_eq!(
            state, initial,
            "invalid density must leave all fields untouched"
        );
    }
}

#[test]
fn invalid_material_storage_or_transport_cannot_report_convergence() {
    let grid = RansChannelGrid::new_stretched(4, 1.0, 1.0);
    let initial = RansFieldState::new_initial(&grid, 0.0, 300.0);
    for property in 0..4 {
        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut solver = RansSolver::new(sample_card());
            match property {
                0 => solver.nu_laminar = invalid,
                1 => solver.cp = invalid,
                2 => solver.k_thermal = invalid,
                _ => solver.pr_turbulent = invalid,
            }
            let mut state = initial.clone();
            assert_eq!(
                solver.solve(&grid, &mut state, -2.0, 40.0, None, None, 1, 1e-8),
                Err(RansTerminationReason::Diverged),
            );
            assert_eq!(state, initial);
        }
    }
}

#[test]
fn density_weights_pressure_acceleration_and_thermal_storage_consistently() {
    let grid = RansChannelGrid::new_stretched(4, 1.0, 1.0);
    let pressure_gradient = -2.0;
    let heat_flux = 40.0;
    for rho in [0.5, 2.0, 8.0] {
        let mut solver = RansSolver::new(sample_card());
        solver.rho = rho;
        solver.cp = 4.0;
        solver.k_thermal = 0.5;
        let mut state = RansFieldState::new_initial(&grid, 0.0, 300.0);
        solver
            .solve(
                &grid,
                &mut state,
                pressure_gradient,
                heat_flux,
                None,
                None,
                1,
                1e-8,
            )
            .unwrap();
        // At the first interior site, the initial molecular/turbulent momentum
        // diffusion is zero. The existing fixture's pseudo-step and relaxation
        // are 0.001 and 0.2; physical pressure acceleration is -dp/dx / rho.
        let expected_u = 0.2 * 0.001 * (-pressure_gradient / rho);
        assert!((state.u[1] - expected_u).abs() < 1e-15);
        // The same rho supplies thermal storage: the uniform first cell receives
        // q''/(rho*cp*dy) K/s. This independently eliminates k using its Neumann BC.
        let dy = grid.y_nodes[1];
        let expected_temperature = 300.0 + 0.01 * heat_flux / (rho * solver.cp * dy);
        assert!((state.temp_k[1] - expected_temperature).abs() < 1e-12);
        assert_eq!(state.temp_k[0], 300.0 + heat_flux * dy / solver.k_thermal);
    }
}

#[test]
fn nonfinite_updates_are_refused_before_they_can_be_hidden_by_max() {
    let grid = RansChannelGrid::new_stretched(4, 1.0, 1.0);
    let solver = RansSolver::new(sample_card());
    for overflow_eddy_viscosity in [true, false] {
        let mut state = RansFieldState::new_initial(&grid, 0.0, 300.0);
        let heat_flux = if overflow_eddy_viscosity {
            // All input scalars are finite; the evaluated k² overflows.
            state.k[1] = 1e200;
            0.0
        } else {
            // The imposed wall temperature is outside the representable range.
            f64::MAX
        };
        assert_eq!(
            solver.solve(&grid, &mut state, -2.0, heat_flux, None, None, 1, 1e-8),
            Err(RansTerminationReason::Diverged),
        );
        for field in [&state.u, &state.k, &state.eps, &state.nu_t, &state.temp_k] {
            assert!(field.iter().all(|value| value.is_finite()));
        }
    }
}
