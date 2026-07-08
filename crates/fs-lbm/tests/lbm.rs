//! Battery for the D2Q9 lattice Boltzmann core (fs-lbm). Covers the
//! equilibrium moments, mass conservation, the load-bearing physical check
//! (steady Poiseuille channel flow matches the analytic parabola), and the
//! lattice-scaling assistant's stability bookkeeping.

use fs_lbm::{Color, Lbm, MACH_LIMIT, equilibrium, plan_scaling, poiseuille_analytic};

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
    assert!((lbm.total_mass() - m0).abs() < 1e-9, "mass drifted");
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
