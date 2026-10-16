//! WS1 3-D feature-preservation battery (bead eg62): thermal
//! Rayleigh–Bénard onset and power-law rheology on D3Q19, certified to
//! the same bars as the 2-D extensions battery. Verdict-JSON style.
//!
//! The determinism hash stays a replay-checked CANDIDATE (not a frozen
//! golden) until the GOLDEN_POLICY four-quadrant ceremony, mirroring the
//! 40p2 precedent.

use fs_lbm::d3q19::{
    PlatesGrid3, ThermalLbm3, gbeta_for_rayleigh3, plate_channel_flow3, update_tau3,
};
use fs_lbm::rheology::{Rheology, powerlaw_poiseuille_analytic};

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

/// lbm3x-001: power-law plate-channel profiles for n ∈ {0.5, 1.0, 1.5}
/// against the SHARED 2-D analytic (the plate geometry is 1-D in z),
/// each to the 2-D battery's 3% bar. n = 1 doubles as the Newtonian
/// consistency check through the same constitutive machinery.
#[test]
fn lbm3x_001_powerlaw_plate_profiles() {
    let (nx, ny, nz) = (4usize, 4, 33);
    let gx = 1e-5;
    for (n, k, max_steps) in [
        (0.5f64, 0.004, 80_000),
        (1.0, 0.05, 60_000),
        (1.5, 0.6, 60_000),
    ] {
        let (grid, steps, stats) =
            plate_channel_flow3(nx, ny, nz, gx, Rheology::PowerLaw { k, n }, max_steps);
        let mut worst = 0.0f64;
        let mut peak = 0.0f64;
        for z in 1..=nz {
            let got = grid.moments3(grid.idx(0, 0, z)).1[0];
            let want = powerlaw_poiseuille_analytic(gx, k, n, nz, z - 1);
            peak = peak.max(want);
            worst = worst.max((got - want).abs());
        }
        let rel = worst / peak;
        verdict(
            &format!("lbm3x-001-powerlaw-n{n:.1}"),
            rel < 0.03,
            &format!(
                "n={n} profile worst dev {rel:.4} of peak {peak:.3e} ({steps} steps, {} floored, {} capped, tau {:.3}..{:.3})",
                stats.floored, stats.capped, stats.tau_range.0, stats.tau_range.1
            ),
        );
    }
}

/// lbm3x-002: Rayleigh–Bénard onset bracket between plates — the seeded
/// roll DECAYS at Ra = 1200 and GROWS at Ra = 2500 (rigid-rigid
/// Ra_c ≈ 1708), and the convecting state transports heat (Nu > 1) —
/// the identical physics gate the 2-D battery pins.
#[test]
fn lbm3x_002_rayleigh_benard_onset_3d() {
    let (nx, ny, nz) = (24usize, 4, 12);
    let (tau_f, tau_g) = (0.7f64, 0.7);
    let run = |ra: f64| -> (f64, f64, f64) {
        let gbeta = gbeta_for_rayleigh3(ra, nz, tau_f, tau_g);
        let mut sim = ThermalLbm3::slab(nx, ny, nz, tau_f, tau_g, gbeta);
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
        "lbm3x-002-subcritical-decay",
        ke2_lo < ke1_lo,
        &format!("Ra=1200: KE {ke1_lo:.3e} -> {ke2_lo:.3e}"),
    );
    verdict(
        "lbm3x-002-supercritical-growth",
        ke2_hi > ke1_hi,
        &format!("Ra=2500: KE {ke1_hi:.3e} -> {ke2_hi:.3e}"),
    );
    verdict(
        "lbm3x-002-nusselt",
        nu_hi > 1.0,
        &format!("convecting Nu = {nu_hi:.3}"),
    );
}

/// lbm3x-003: conservation ledger — total flow mass stays within 1e-11
/// relative of its initial value through a coupled thermal run (Guo
/// forcing is mass-neutral, plates bounce back, x/y wrap).
#[test]
fn lbm3x_003_thermal_mass_conservation() {
    let gbeta = gbeta_for_rayleigh3(2500.0, 8, 0.7, 0.7);
    let mut sim = ThermalLbm3::slab(8, 4, 8, 0.7, 0.7, gbeta);
    sim.perturb(1e-4);
    let m0 = sim.grid.total_mass();
    let mut worst = 0.0f64;
    for _ in 0..300 {
        sim.step();
        worst = worst.max(((sim.grid.total_mass() - m0) / m0).abs());
    }
    verdict(
        "lbm3x-003-mass-ledger",
        worst < 1e-11,
        &format!("worst relative mass drift {worst:.3e} over 300 coupled steps"),
    );
}

/// lbm3x-004: fixed plate temperatures are held exactly by construction
/// and the pure-conduction steady state reproduces the linear profile.
#[test]
fn lbm3x_004_plate_temperatures_and_conduction() {
    let mut sim = ThermalLbm3::slab(4, 4, 8, 0.7, 0.7, 0.0);
    for _ in 0..2000 {
        sim.step();
    }
    let bottom = sim.temperature(0, 0, 0);
    let top = sim.temperature(0, 0, sim.grid.nz - 1);
    verdict(
        "lbm3x-004-plate-temperatures",
        (bottom - sim.t_bottom).abs() < 1e-12 && (top - sim.t_top).abs() < 1e-12,
        &format!("bottom={bottom:.3}, top={top:.3}"),
    );
    let mut worst = 0.0f64;
    for z in 1..sim.grid.nz - 1 {
        let want = sim.t_bottom
            + (sim.t_top - sim.t_bottom) * ((z as f64 - 0.5) / (sim.grid.nz - 2) as f64);
        let got = sim.temperature(1, 2, z);
        worst = worst.max((got - want).abs());
    }
    verdict(
        "lbm3x-004-conduction-profile",
        worst < 1e-6,
        &format!("worst conduction deviation {worst:.3e}"),
    );
}

/// lbm3x-005: seeded determinism CANDIDATE — the coupled thermal state
/// and a rheology plate channel hash replay-stable over exact bits.
/// Freezing into golden-couplings.json is deferred to the four-quadrant
/// ceremony per GOLDEN_POLICY (40p2 precedent).
#[test]
fn lbm3x_005_candidate_hash_is_replay_stable() {
    let run = || -> u64 {
        let mut accumulator = 0xcbf2_9ce4_8422_2325u64;
        let mut feed = |bytes: &[u8]| {
            for byte in bytes {
                accumulator ^= u64::from(*byte);
                accumulator = accumulator.wrapping_mul(0x0000_0100_0000_01b3);
            }
        };
        let gbeta = gbeta_for_rayleigh3(2500.0, 8, 0.7, 0.7);
        let mut sim = ThermalLbm3::slab(8, 4, 8, 0.7, 0.7, gbeta);
        sim.grid.perturb(0x40_02_E6, 1e-4);
        for _ in 0..96 {
            sim.step();
        }
        for (x, y, z) in [(0, 0, 1), (3, 2, 4), (7, 3, 8), (5, 1, 6)] {
            let (rho, u) = sim.grid.moments3(sim.grid.idx(x, y, z));
            feed(&rho.to_bits().to_le_bytes());
            for component in u {
                feed(&component.to_bits().to_le_bytes());
            }
            feed(&sim.temperature(x, y, z).to_bits().to_le_bytes());
        }
        feed(&sim.grid.total_mass().to_bits().to_le_bytes());

        let (grid, _, _) = plate_channel_flow3(
            4,
            4,
            9,
            1e-5,
            Rheology::PowerLaw { k: 0.016, n: 0.8 },
            2_000,
        );
        for z in [1usize, 4, 9] {
            let (rho, u) = grid.moments3(grid.idx(0, 0, z));
            feed(&rho.to_bits().to_le_bytes());
            for component in u {
                feed(&component.to_bits().to_le_bytes());
            }
            feed(&grid.tau[grid.idx(0, 0, z)].to_bits().to_le_bytes());
        }
        accumulator
    };
    let first = run();
    let second = run();
    println!("{{\"test\":\"lbm3x-005-candidate\",\"hash\":\"{first:#018x}\"}}");
    assert_ne!(first, 0);
    assert_eq!(
        first, second,
        "coupled 3-D candidate hash is not replay-stable"
    );
}

/// lbm3x-006: the rheology ledger reports honestly — a Newtonian law
/// through update_tau3 leaves τ exactly uniform with an empty ledger,
/// and invalid parameters refuse before NaN physics.
#[test]
fn lbm3x_006_ledger_honesty_and_refusals() {
    let mut grid = PlatesGrid3::plates(4, 4, 8, 1.0);
    let stats = update_tau3(&mut grid, Rheology::Newtonian { nu: 0.1 });
    let expected_tau = 0.1 / (1.0 / 3.0) + 0.5;
    let uniform = grid
        .tau
        .iter()
        .enumerate()
        .filter(|(i, _)| !grid.is_wall_layer(i / (grid.nx * grid.ny)))
        .all(|(_, t)| (t - expected_tau).abs() < 1e-15);
    verdict(
        "lbm3x-006-newtonian-uniform-tau",
        uniform && stats.floored == 0 && stats.capped == 0,
        &format!(
            "tau -> {expected_tau:.3}, floored {}, capped {}",
            stats.floored, stats.capped
        ),
    );

    let bad_rheology =
        std::panic::catch_unwind(|| Rheology::PowerLaw { k: -1.0, n: 0.8 }.viscosity(1.0));
    let bad_rayleigh = std::panic::catch_unwind(|| gbeta_for_rayleigh3(1200.0, 0, 0.7, 0.7));
    let bad_grid = std::panic::catch_unwind(|| PlatesGrid3::plates(0, 4, 8, 1.0));
    verdict(
        "lbm3x-006-invalid-parameter-rejection",
        bad_rheology.is_err() && bad_rayleigh.is_err() && bad_grid.is_err(),
        "invalid rheology, Rayleigh setup, and grid are rejected before NaNs propagate",
    );
}
