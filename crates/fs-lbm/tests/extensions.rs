//! fs-lbm extensions battery (bead tfz.19): non-Newtonian channel
//! profiles, Rayleigh–Bénard onset, level-jump refinement, and the
//! free-surface mass/benchmark/bracketing gates. Verdict-JSON style.

use fs_lbm::freesurface::{ContactModel, FreeSurface, dam_break, surge_front};
use fs_lbm::rheology::{Rheology, channel_flow, powerlaw_poiseuille_analytic};
use fs_lbm::thermal::{ThermalLbm, gbeta_for_rayleigh};
use fs_lbm::{Cell, Grid, Q};

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

#[test]
fn lbm_103_gas_cells_do_not_feed_plain_streaming() {
    let mut grid = Grid::uniform(3, 3, 0.8);
    grid.periodic_x = false;
    grid.periodic_y = false;
    let gas = grid.idx(0, 1);
    grid.flags[gas] = Cell::Gas;
    grid.f[gas] = [1000.0; Q];

    let fluid = grid.idx(1, 1);
    let mut scratch = Vec::new();
    grid.step(&mut scratch);

    verdict(
        "lbm-103-gas-neighbor-bounce",
        grid.f[fluid][1] < 1.0,
        &format!(
            "east-moving population after pulling from gas boundary is {:.3e}",
            grid.f[fluid][1]
        ),
    );
}

#[test]
fn lbm_104_thermal_wall_temperatures_match_boundary_values() {
    let sim = ThermalLbm::slab(4, 4, 0.7, 0.7, 0.0);
    let bottom = sim.temperature(0, 0);
    let top = sim.temperature(0, sim.grid.ny - 1);
    verdict(
        "lbm-104-wall-temperatures",
        (bottom - sim.t_bottom).abs() < 1e-12 && (top - sim.t_top).abs() < 1e-12,
        &format!("bottom={bottom:.3}, top={top:.3}"),
    );
}

#[test]
fn lbm_105_invalid_parameters_are_rejected_before_nan_physics() {
    let bad_grid = std::panic::catch_unwind(|| Grid::uniform(0, 2, 0.7));
    let bad_rheology =
        std::panic::catch_unwind(|| Rheology::PowerLaw { k: -1.0, n: 0.8 }.viscosity(1.0));
    let bad_rayleigh = std::panic::catch_unwind(|| gbeta_for_rayleigh(1200.0, 0, 0.7, 0.7));

    verdict(
        "lbm-105-invalid-parameter-rejection",
        bad_grid.is_err() && bad_rheology.is_err() && bad_rayleigh.is_err(),
        "invalid grid, rheology, and Rayleigh setup are rejected before NaNs propagate",
    );
}

/// lbm-104: the free-surface mass ledger is STRICT — over a full dam
/// break with conversions in both directions, total tracked mass
/// (fluid Σf + interface m + carry) stays within 1e-10 relative of
/// its initial value at EVERY step. The make-or-break audit.
#[test]
fn lbm_104_mass_ledger() {
    let mut fs = dam_break(40, 24, 8, 1e-4, 0.0, ContactModel::Neutral);
    let m0 = fs.ledger_mass();
    let mut worst = 0.0f64;
    for _ in 0..600 {
        fs.step();
        worst = worst.max(((fs.ledger_mass() - m0) / m0).abs());
    }
    let c = fs.conversions;
    verdict(
        "lbm-104-mass-ledger",
        worst < 1e-10,
        &format!(
            "worst rel drift {worst:.2e} over 600 steps (conversions: {} to-fluid, {} to-gas, {} gas->int, {} fluid->int)",
            c.to_fluid, c.to_gas, c.gas_to_interface, c.fluid_to_interface
        ),
    );
}

/// lbm-105: dam-break surge front vs the Martin–Moyce-style envelope.
/// COARSE-LATTICE HONESTY BAND: after the initial transient the
/// nondimensional front z = x/a must lie inside [1 + 0.6 t*, 1 + 2.2 t*]
/// (t* = t·sqrt(2g/a)) — the experimental surge sits near 1 + 1.5 t*;
/// the fine-lattice quantitative comparison is perf-lane scope.
#[test]
fn lbm_105_dam_break_front() {
    let a = 10usize;
    let g = 5e-5;
    let mut fs = dam_break(64, 28, a, g, 0.0, ContactModel::Neutral);
    let tstar = |t: usize| (t as f64) * (2.0 * g / a as f64).sqrt();
    let mut ok = true;
    let mut detail = String::new();
    let mut checked = 0;
    for t in 1..=1200 {
        fs.step();
        let ts = tstar(t);
        if ts > 0.5 && ts < 2.0 && t % 150 == 0 {
            let z = surge_front(&fs) as f64 / a as f64;
            let (lo, hi) = (0.6f64.mul_add(ts, 1.0), 2.2f64.mul_add(ts, 1.0));
            use std::fmt::Write as _;
            let _ = write!(detail, "t*={ts:.2}: z={z:.2} in [{lo:.2},{hi:.2}]; ");
            if z < lo || z > hi {
                ok = false;
            }
            checked += 1;
        }
    }
    verdict(
        "lbm-105-dam-break-envelope",
        ok && checked >= 3,
        &format!("HONESTY BAND (coarse lattice): {detail}"),
    );
}

/// lbm-106 (G3): 90° rotation equivariance — a dam break with gravity
/// −ŷ in a W×H box maps cell-for-cell onto the same dam break with
/// gravity −x̂ in the transposed box, to roundoff. Tilt schedules are
/// rotations of g, so this pins the whole moving-frame path.
#[test]
fn lbm_106_rotation_equivariance() {
    let (w, h, a, g) = (36usize, 22usize, 7usize, 1e-4);
    let mut sim1 = dam_break(w, h, a, g, 0.0, ContactModel::Neutral);
    // Build the rotated twin by hand: (x, y) -> (y, x) with gravity
    // along -x. The dam column a wide x 2a tall becomes 2a wide x a
    // tall... which is NOT the same physical problem unless we map
    // the geometry exactly: transpose the box and the column, and
    // point gravity along -x.
    let mut grid = Grid::uniform(h, w, 0.55);
    grid.periodic_x = false;
    grid.periodic_y = false;
    grid.g = [-g, 0.0];
    for i in 0..h * w {
        grid.flags[i] = Cell::Gas;
    }
    for x in 0..h {
        let b = grid.idx(x, 0);
        grid.flags[b] = Cell::Wall;
        let t = grid.idx(x, w - 1);
        grid.flags[t] = Cell::Wall;
    }
    for y in 0..w {
        let l = grid.idx(0, y);
        grid.flags[l] = Cell::Wall;
        let r = grid.idx(h - 1, y);
        grid.flags[r] = Cell::Wall;
    }
    for y in 1..=a.min(w - 2) {
        for x in 1..=(2 * a).min(h - 2) {
            let i = grid.idx(x, y);
            grid.flags[i] = Cell::Fluid;
        }
    }
    let mut sim2 = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    for _ in 0..300 {
        sim1.step();
        sim2.step();
    }
    // Compare fills cell-for-cell under (x, y) -> (y, x).
    let mut worst = 0.0f64;
    for y in 0..h {
        for x in 0..w {
            let f1 = sim1.fill(sim1.grid.idx(x, y));
            let f2 = sim2.fill(sim2.grid.idx(y, x));
            worst = worst.max((f1 - f2).abs());
        }
    }
    let dm = ((sim1.ledger_mass() - sim2.ledger_mass()) / sim1.ledger_mass()).abs();
    verdict(
        "lbm-106-rotation-equivariance",
        worst < 1e-9 && dm < 1e-12,
        &format!("fill dev {worst:.2e}, ledger dev {dm:.2e} after 300 steps"),
    );
}

/// lbm-107: contact-line MODEL BRACKETING + breaking-jet qualitative
/// battery. (a) The same wetting-driven fixture run under the neutral
/// vs wetting contact models produces a REPORTED sensitivity band on
/// the surge front — the honest statement replacing pretended
/// contact-angle certainty. (b) A perturbed liquid strip under
/// surface tension necks and breaks into multiple fragments
/// (Plateau–Rayleigh in the lattice's qualitative regime).
#[test]
fn lbm_107_bracketing_and_jet() {
    // (a) Contact bracket on a surface-tension dam break.
    let run = |cm: ContactModel| -> usize {
        let mut fs = dam_break(48, 24, 8, 8e-5, 0.002, cm);
        for _ in 0..500 {
            fs.step();
        }
        surge_front(&fs)
    };
    let fa = run(ContactModel::Neutral);
    let fb = run(ContactModel::Wetting);
    let band = fa.abs_diff(fb);
    verdict(
        "lbm-107-contact-bracket",
        band > 0,
        &format!(
            "surge front: neutral {fa} vs wetting {fb} cells -> sensitivity band {band} cells (REPORTED, not hidden)"
        ),
    );
    // (b) Plateau–Rayleigh necking: a varicose-perturbed strip breaks.
    let (nx, ny) = (96usize, 24);
    let mut grid = Grid::uniform(nx, ny, 0.55);
    grid.periodic_x = true;
    grid.periodic_y = false;
    for i in 0..nx * ny {
        grid.flags[i] = Cell::Gas;
    }
    for x in 0..nx {
        let b = grid.idx(x, 0);
        grid.flags[b] = Cell::Wall;
        let t = grid.idx(x, ny - 1);
        grid.flags[t] = Cell::Wall;
    }
    let mid = ny / 2;
    for x in 0..nx {
        // Thickness 5 with a 2-cell varicose perturbation, two waves.
        let pert =
            (2.0 * (std::f64::consts::TAU * 2.0 * x as f64 / nx as f64).cos()).round() as i64;
        let half = (5 + pert.max(-4)) / 2;
        let half = usize::try_from(half.max(1)).expect("positive");
        for y in mid.saturating_sub(half)..=(mid + half).min(ny - 2) {
            let i = grid.idx(x, y);
            grid.flags[i] = Cell::Fluid;
        }
    }
    let mut jet = FreeSurface::new(grid, 0.01, ContactModel::Neutral);
    let frags0 = jet.fragment_count();
    let m0 = jet.ledger_mass();
    for _ in 0..800 {
        jet.step();
    }
    let frags = jet.fragment_count();
    let drift = ((jet.ledger_mass() - m0) / m0).abs();
    verdict(
        "lbm-107-jet-breakup",
        frags0 == 1 && frags >= 2 && drift < 1e-10,
        &format!(
            "fragments {frags0} -> {frags} after 800 steps (QUALITATIVE gate), ledger drift {drift:.2e}"
        ),
    );
}
