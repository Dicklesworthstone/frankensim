//! Battery for the D2Q9 lattice Boltzmann core (fs-lbm). Covers the
//! equilibrium moments, mass conservation, the load-bearing physical check
//! (steady Poiseuille channel flow matches the analytic parabola), and the
//! lattice-scaling assistant's stability bookkeeping.

use fs_lbm::{
    Cell, Color, Grid, Lbm, MACH_LIMIT, Q, equilibrium, plan_scaling, poiseuille_analytic,
};

#[test]
fn the_equilibrium_recovers_its_moments() {
    let (rho, ux, uy) = (1.0, 0.05, -0.02);
    let f = equilibrium(rho, ux, uy);
    let sum: f64 = f.iter().sum();
    assert!((sum - rho).abs() < 1e-12); // density
    // momentum: Σ eₓ fᵢ = ρ uₓ (D2Q9 velocities).
    let ex = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
    let ey = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
    let mx: f64 = f.iter().zip(ex).map(|(fi, e)| fi * e).sum();
    let my: f64 = f.iter().zip(ey).map(|(fi, e)| fi * e).sum();
    assert!((mx - rho * ux).abs() < 1e-12 && (my - rho * uy).abs() < 1e-12);
}

#[test]
fn mass_is_conserved() {
    let mut lbm = Lbm::channel(6, 12, 0.8, 1e-4);
    let m0 = lbm.total_mass();
    lbm.run(200);
    // Mass is conserved BY CONSTRUCTION (collision, forcing, streaming, and
    // bounce-back all preserve Σf), so the only drift is summation roundoff:
    // measured 9.38e-13 over 200 steps on mass 72, BIT-IDENTICAL on aarch64
    // (M4 Pro) and x86-64 (Threadripper 5975WX). Gate at 1e-11 (~10x that
    // roundoff floor) so the CONTRACT's "mass is conserved" claim is verified
    // to roundoff and a future systematic per-step leak is actually caught —
    // the old 1e-9 bound was ~1000x loose and would have passed a real
    // ~5e-12/step leak.
    assert!((lbm.total_mass() - m0).abs() < 1e-11, "mass drifted");
    assert!((m0 - f64::from(6 * 12)).abs() < 1e-9); // unit density
}

#[test]
fn poiseuille_flow_matches_the_analytic_parabola() {
    let (nx, ny, tau, gx) = (4, 25, 0.8, 1e-5);
    let mut lbm = Lbm::channel(nx, ny, tau, gx);
    lbm.run(20_000); // reach steady state
    let profile = lbm.x_velocity_profile();
    let nu = lbm.viscosity();
    // the profile matches the analytic parabola at every row (halfway
    // bounce-back resolves the quadratic exactly).
    let mut max_rel = 0.0_f64;
    for (y, &u) in profile.iter().enumerate() {
        let a = poiseuille_analytic(gx, nu, ny, y);
        max_rel = max_rel.max((u - a).abs() / a.abs());
    }
    assert!(max_rel < 0.03, "profile off by {max_rel:.4}");
    // and it is a parabola: symmetric with its peak at the centre.
    let mid = ny / 2;
    assert!(profile[mid] > profile[0] && profile[mid] > profile[ny - 1]);
    assert!((profile[1] - profile[ny - 2]).abs() / profile[mid] < 0.02);
}

#[test]
fn the_scaling_assistant_derives_tau_and_flags_stability() {
    // Re 100, L 40 lu, u 0.05 -> nu 0.02, tau 0.56, low Mach -> stable.
    let plan = plan_scaling(100.0, 40.0, 0.05);
    assert!((plan.viscosity - 0.02).abs() < 1e-12);
    assert!((plan.tau - 0.56).abs() < 1e-12);
    assert!(plan.stable && plan.tau_margin > 0.0);
    assert!(plan.mach < MACH_LIMIT);
    // a comfortably-stable plan is verified-color.
    assert!(matches!(plan.color(), Color::Verified { .. }));
}

#[test]
fn the_scaling_assistant_rejects_a_high_mach_plan() {
    // too large a lattice velocity breaks the low-Mach (incompressible) regime.
    let plan = plan_scaling(100.0, 20.0, 0.25);
    assert!(plan.mach > MACH_LIMIT);
    assert!(!plan.stable);
    // an unstable plan is not verified-color.
    assert!(matches!(plan.color(), Color::Estimated { .. }));
}

#[test]
#[should_panic(expected = "must be positive")]
fn the_scaling_assistant_rejects_nonsense() {
    let _ = plan_scaling(-1.0, 20.0, 0.1);
}

#[test]
fn the_solver_is_deterministic() {
    let mut a = Lbm::channel(4, 10, 0.7, 1e-4);
    let mut b = Lbm::channel(4, 10, 0.7, 1e-4);
    a.run(100);
    b.run(100);
    assert_eq!(a.x_velocity_profile(), b.x_velocity_profile());
}

#[test]
fn d2q9_wall_momentum_has_exact_link_sign_and_obstacle_selection() {
    let mut grid = Grid::uniform(3, 3, 0.8);
    grid.periodic_x = false;
    grid.periodic_y = false;
    grid.flags.fill(Cell::Gas);
    let fluid = grid.idx(1, 1);
    let left_wall = grid.idx(0, 1);
    let right_wall = grid.idx(2, 1);
    let top_wall = grid.idx(1, 2);
    grid.flags[fluid] = Cell::Fluid;
    grid.flags[left_wall] = Cell::Wall;
    grid.flags[right_wall] = Cell::Wall;
    grid.flags[top_wall] = Cell::Wall;

    let mut post = vec![[0.0; Q]; grid.nx * grid.ny];
    // D2Q9 directions 1 and 3 point east and west. An east-going population
    // of 0.25 transfers +2f = +0.5 x-momentum to the right wall.
    post[fluid][1] = 0.25;
    // The west-going population would transfer -0.75 to the left wall.
    post[fluid][3] = 0.375;
    // A north-going population of 0.125 transfers +0.25 y-momentum to the
    // top wall, independently pinning the lift-axis sign.
    post[fluid][2] = 0.125;

    let mut right_only = vec![false; grid.nx * grid.ny];
    right_only[right_wall] = true;
    let right_receipt = grid.stream_from_with_wall_momentum(&post, &right_only);
    assert_eq!(right_receipt.measured_links, 1);
    assert_eq!(right_receipt.wall_impulse[0].to_bits(), 0.5f64.to_bits());
    assert_eq!(right_receipt.wall_impulse[1].to_bits(), 0.0f64.to_bits());
    assert_eq!(grid.f[fluid][3].to_bits(), 0.25f64.to_bits());

    let mut both = right_only;
    both[left_wall] = true;
    let both_receipt = grid.stream_from_with_wall_momentum(&post, &both);
    assert_eq!(both_receipt.measured_links, 2);
    assert_eq!(both_receipt.wall_impulse[0].to_bits(), (-0.25f64).to_bits());
    assert_eq!(both_receipt.wall_impulse[1].to_bits(), 0.0f64.to_bits());
    assert_eq!(grid.f[fluid][1].to_bits(), 0.375f64.to_bits());

    let mut top_only = vec![false; grid.nx * grid.ny];
    top_only[top_wall] = true;
    let top_receipt = grid.stream_from_with_wall_momentum(&post, &top_only);
    assert_eq!(top_receipt.measured_links, 1);
    assert_eq!(top_receipt.wall_impulse[0].to_bits(), 0.0f64.to_bits());
    assert_eq!(top_receipt.wall_impulse[1].to_bits(), 0.25f64.to_bits());
    assert_eq!(grid.f[fluid][4].to_bits(), 0.125f64.to_bits());
}

#[test]
fn d2q9_wall_momentum_is_zero_at_rest_and_replays_bitwise() {
    let mut resting = Grid::uniform(5, 5, 0.8);
    let resting_wall = resting.idx(2, 2);
    resting.flags[resting_wall] = Cell::Wall;
    let mut resting_mask = vec![false; resting.nx * resting.ny];
    resting_mask[resting_wall] = true;
    let resting_post = resting.f.clone();
    let resting_receipt = resting.stream_from_with_wall_momentum(&resting_post, &resting_mask);
    assert_eq!(resting_receipt.measured_links, 8);
    assert!(resting_receipt.wall_impulse[0].abs() < 1e-15);
    assert!(resting_receipt.wall_impulse[1].abs() < 1e-15);

    let mut first = Grid::uniform(7, 7, 0.8);
    let wall = first.idx(3, 3);
    let upstream = first.idx(2, 3);
    first.flags[wall] = Cell::Wall;
    first.f[upstream] = equilibrium(1.0, 0.04, 0.01);
    let mut second = first.clone();
    let mut legacy = first.clone();
    let mut mask = vec![false; first.nx * first.ny];
    mask[wall] = true;
    let (mut first_scratch, mut second_scratch) = (Vec::new(), Vec::new());
    let mut legacy_scratch = Vec::new();
    let mut observed_nonzero_impulse = false;

    for _ in 0..12 {
        let first_receipt = first.step_with_wall_momentum(&mut first_scratch, &mask);
        let second_receipt = second.step_with_wall_momentum(&mut second_scratch, &mask);
        legacy.step(&mut legacy_scratch);
        observed_nonzero_impulse |= first_receipt.wall_impulse != [0.0, 0.0];
        assert_eq!(first_receipt.measured_links, second_receipt.measured_links);
        assert_eq!(
            first_receipt.wall_impulse.map(f64::to_bits),
            second_receipt.wall_impulse.map(f64::to_bits)
        );
    }
    assert!(observed_nonzero_impulse);
    for (first_cell, (second_cell, legacy_cell)) in
        first.f.iter().zip(second.f.iter().zip(&legacy.f))
    {
        assert_eq!(first_cell.map(f64::to_bits), second_cell.map(f64::to_bits));
        assert_eq!(first_cell.map(f64::to_bits), legacy_cell.map(f64::to_bits));
    }
}

#[test]
fn d2q9_wall_momentum_refuses_invalid_masks_before_advancing() {
    let mut grid = Grid::uniform(3, 3, 0.8);
    let original = grid.f.clone();
    let mut scratch = Vec::new();
    let non_wall_mask = vec![true; grid.nx * grid.ny];
    let refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = grid.step_with_wall_momentum(&mut scratch, &non_wall_mask);
    }));
    assert!(refusal.is_err());
    assert_eq!(grid.f, original);
    assert!(scratch.is_empty());

    let short_mask = vec![false; grid.nx * grid.ny - 1];
    let short_refusal = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = grid.step_with_wall_momentum(&mut scratch, &short_mask);
    }));
    assert!(short_refusal.is_err());
    assert_eq!(grid.f, original);
    assert!(scratch.is_empty());
}
