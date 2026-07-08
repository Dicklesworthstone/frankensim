//! fs-lbm extensions battery (bead tfz.19): non-Newtonian channel
//! profiles, Rayleigh–Bénard onset, level-jump refinement, and the
//! free-surface mass/benchmark/bracketing gates. Verdict-JSON style.

use fs_lbm::rheology::{Rheology, channel_flow, powerlaw_poiseuille_analytic};
use fs_lbm::thermal::{ThermalLbm, gbeta_for_rayleigh};

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

/// lbm-101: power-law Poiseuille — the local shear-rate τ adaptation
/// reproduces the analytic profile u ∝ (H^(1+1/n) − |y|^(1+1/n)) for
/// a shear-thinning n = 0.8 fluid, and the Newtonian limit (n = 1)
/// matches the parabola. τ-floor events are ledgered (zero here —
/// the fixture stays inside the stable window).
#[test]
fn lbm_101_powerlaw_poiseuille() {
    let (nx, ny) = (4usize, 33);
    // Wall τ ≈ 1.1 by design (ν_wall ≈ 0.19 at the wall shear rate);
    // the centerline plug hits TAU_CAP and is ledgered.
    let gx = 1e-5;
    let (k, n) = (0.016, 0.8);
    let (grid, steps, stats) = channel_flow(nx, ny, gx, Rheology::PowerLaw { k, n }, 60_000);
    // Compare the measured profile against the analytic one.
    let mut worst_rel = 0.0f64;
    let mut peak = 0.0f64;
    for y in 1..=ny {
        let got = grid.moments(grid.idx(0, y)).u[0];
        let want = powerlaw_poiseuille_analytic(gx, k, n, ny, y - 1);
        peak = peak.max(want);
        worst_rel = worst_rel.max((got - want).abs());
    }
    let rel = worst_rel / peak;
    verdict(
        "lbm-101-powerlaw-profile",
        rel < 0.03,
        &format!(
            "n=0.8 profile worst dev {rel:.4} of peak {peak:.3e} ({steps} steps, {} floored, {} capped, tau {:.3}..{:.3})",
            stats.floored, stats.capped, stats.tau_range.0, stats.tau_range.1
        ),
    );
    // Newtonian limit through the same machinery.
    let nu = 0.1;
    let (grid_n, _, _) = channel_flow(nx, ny, gx, Rheology::Newtonian { nu }, 40_000);
    let mut worst_n = 0.0f64;
    let mut peak_n = 0.0f64;
    for y in 1..=ny {
        let got = grid_n.moments(grid_n.idx(0, y)).u[0];
        let want = fs_lbm::poiseuille_analytic(gx, nu, ny, y - 1);
        peak_n = peak_n.max(want);
        worst_n = worst_n.max((got - want).abs());
    }
    verdict(
        "lbm-101-newtonian-limit",
        worst_n / peak_n < 0.01,
        &format!("n=1 worst dev {:.4} of peak", worst_n / peak_n),
    );
    // Carreau plateaus: at tiny shear the apparent viscosity is ν0.
    let carreau = Rheology::Carreau {
        nu0: 0.2,
        nu_inf: 0.01,
        lambda: 100.0,
        n: 0.5,
    };
    let lo = carreau.viscosity(1e-9);
    let hi = carreau.viscosity(1e3);
    verdict(
        "lbm-101-carreau-plateaus",
        (lo - 0.2).abs() < 1e-6 && (hi - 0.01).abs() < 1e-3,
        &format!("nu(0)={lo:.4} nu(inf)={hi:.4}"),
    );
}

/// lbm-102: Rayleigh–Bénard onset bracket — a seeded convection mode
/// DECAYS at Ra = 1200 and GROWS at Ra = 2500 (critical Ra ≈ 1708
/// for rigid-rigid), and the convecting state transports heat
/// (Nu > 1). Physics the scheme cannot fake.
#[test]
fn lbm_102_rayleigh_benard_onset() {
    let (nx, ny) = (24usize, 12);
    let (tau_f, tau_g) = (0.7f64, 0.7);
    let run = |ra: f64| -> (f64, f64, f64) {
        let gbeta = gbeta_for_rayleigh(ra, ny, tau_f, tau_g);
        let mut sim = ThermalLbm::slab(nx, ny, tau_f, tau_g, gbeta);
        // Settle the conduction state, then seed the onset mode.
        for _ in 0..500 {
            sim.step();
        }
        sim.perturb(1e-5);
        for _ in 0..1500 {
            sim.step();
        }
        let ke1 = sim.kinetic_energy();
        for _ in 0..3000 {
            sim.step();
        }
        let ke2 = sim.kinetic_energy();
        (ke1, ke2, sim.nusselt())
    };
    let (ke1_lo, ke2_lo, _) = run(1200.0);
    let (ke1_hi, ke2_hi, nu_hi) = run(2500.0);
    verdict(
        "lbm-102-subcritical-decay",
        ke2_lo < ke1_lo,
        &format!("Ra=1200: KE {ke1_lo:.3e} -> {ke2_lo:.3e}"),
    );
    verdict(
        "lbm-102-supercritical-growth",
        ke2_hi > ke1_hi,
        &format!("Ra=2500: KE {ke1_hi:.3e} -> {ke2_hi:.3e}"),
    );
    verdict(
        "lbm-102-nusselt",
        nu_hi > 1.0,
        &format!("convecting Nu = {nu_hi:.3}"),
    );
}
