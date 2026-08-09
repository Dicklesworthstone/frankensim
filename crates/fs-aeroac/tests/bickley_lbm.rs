//! The bead's verification chain, link 1: the fs-lbm D2Q9 Bickley
//! jet reproduces the INVISCID Rayleigh instability growth rate from
//! [`fs_aeroac::bickley`] in the linear regime — the analytic anchor
//! that certifies the base flow BEFORE any aeroacoustic source
//! extraction is trusted.
//!
//! Method: periodic domain seeded with `U = U0 sech^2((y-yc)/b)`
//! plus a small SINUOUS-symmetry transverse perturbation at the
//! single wavenumber the box admits (`alpha = 2 pi b / nx`). The
//! transverse kinetic amplitude grows as `e^{sigma t}`; the log-slope
//! over a post-transient window, nondimensionalized by `U0/b`, is
//! compared against `alpha Im(c)` from the shooting oracle. The LBM
//! rate sits BELOW the inviscid one (finite Reynolds number damps),
//! and the measured rate DRIFTS DOWN over long runs because the base
//! jet itself diffuses (executed at Re 120: ~nu/b^2 per step, a 19%
//! U0 decay over 4200 steps broke window consistency) — hence Re 240
//! and an early fit window. Gate direction and margin are measured
//! and disclosed; Mach + linearity diagnostics guard the regime.

use fs_aeroac::bickley::{JetSymmetry, bickley_rayleigh_mode};
use fs_lbm::core2::Grid;
use fs_math::c64::C64;

/// One configuration of the growth-rate measurement.
struct JetConfig {
    nx: usize,
    ny: usize,
    b: f64,
    u0: f64,
    tau: f64,
    fit_lo: f64,
    fit_mid: f64,
    fit_hi: f64,
    steps: usize,
}

const EPS: f64 = 2.0e-3; // seed amplitude relative to u0

fn sech2(x: f64) -> f64 {
    let t = x.tanh();
    1.0 - t * t
}

fn seeded_jet(c: &JetConfig) -> Grid {
    let mut grid = Grid::uniform(c.nx, c.ny, c.tau);
    let yc = c.ny as f64 / 2.0 - 0.5;
    let alpha_lu = 2.0 * core::f64::consts::PI / c.nx as f64;
    for x in 0..c.nx {
        for y in 0..c.ny {
            let s = sech2((y as f64 - yc) / c.b);
            let ux = c.u0 * s;
            // Sinuous symmetry: transverse velocity EVEN in y.
            let vy = EPS * c.u0 * (alpha_lu * x as f64).cos() * s;
            let i = grid.idx(x, y);
            grid.f[i] = fs_lbm::equilibrium(1.0, ux, vy);
        }
    }
    grid
}

/// Transverse kinetic amplitude sqrt(sum vy^2) over the domain.
fn transverse_amplitude(grid: &Grid) -> f64 {
    let mut sum = 0.0;
    for x in 0..grid.nx {
        for y in 0..grid.ny {
            let m = grid.moments(grid.idx(x, y));
            sum += m.u[1] * m.u[1];
        }
    }
    sum.sqrt()
}

/// Run one configuration; return (nondimensional LBM growth rate,
/// max speed, window-consistency ratio).
fn measure(c: &JetConfig) -> (f64, f64, f64) {
    let mut grid = seeded_jet(c);
    let mut scratch = Vec::new();
    let mut log_amp: Vec<(f64, f64)> = Vec::new();
    for t in 0..c.steps {
        grid.step(&mut scratch);
        if t % 50 == 0 {
            log_amp.push((t as f64, transverse_amplitude(&grid).ln()));
        }
    }
    let mut max_speed = 0.0f64;
    for x in 0..c.nx {
        for y in 0..c.ny {
            let m = grid.moments(grid.idx(x, y));
            let sp = (m.u[0] * m.u[0] + m.u[1] * m.u[1]).sqrt();
            max_speed = max_speed.max(sp);
        }
    }
    let fit = |lo: f64, hi: f64| -> f64 {
        let pts: Vec<&(f64, f64)> = log_amp
            .iter()
            .filter(|(t, _)| *t >= lo && *t <= hi)
            .collect();
        let n = pts.len() as f64;
        let sx: f64 = pts.iter().map(|(t, _)| t).sum();
        let sy: f64 = pts.iter().map(|(_, a)| a).sum();
        let sxx: f64 = pts.iter().map(|(t, _)| t * t).sum();
        let sxy: f64 = pts.iter().map(|(t, a)| t * a).sum();
        (n * sxy - sx * sy) / (n * sxx - sx * sx)
    };
    let slope_a = fit(c.fit_lo, c.fit_mid);
    let slope_b = fit(c.fit_mid, c.fit_hi);
    let slope = fit(c.fit_lo, c.fit_hi);
    ((slope * c.b / c.u0), max_speed, (slope_a - slope_b) / slope)
}

#[test]
fn lbm_bickley_jet_converges_to_rayleigh_growth_rate() {
    // Inviscid oracle at the shared box wavenumber alpha = 2 pi b/nx
    // (both configurations keep b/nx = 1/8).
    let alpha = 2.0 * core::f64::consts::PI / 8.0;
    let mode = bickley_rayleigh_mode(alpha, JetSymmetry::Sinuous, C64::new(0.6, 0.14), 14.0, 2048)
        .expect("oracle mode");
    let sigma_inviscid = mode.growth_rate;

    let coarse = JetConfig {
        nx: 64,
        ny: 112,
        b: 8.0,
        u0: 0.05,
        tau: 0.505, // Re = 240
        fit_lo: 1200.0,
        fit_mid: 2600.0,
        fit_hi: 4000.0,
        steps: 4200,
    };
    let fine = JetConfig {
        nx: 96,
        ny: 168,
        b: 12.0,
        u0: 0.06,
        tau: 0.505, // Re = 432
        fit_lo: 1200.0,
        fit_mid: 2600.0,
        fit_hi: 4000.0,
        steps: 4200,
    };
    let (sig_c, mach_c, cons_c) = measure(&coarse);
    let (sig_f, mach_f, cons_f) = measure(&fine);
    for (name, mach, cons) in [("coarse", mach_c, cons_c), ("fine", mach_f, cons_f)] {
        assert!(mach < 0.15, "{name}: low-Mach diagnostic: {mach}");
        assert!(
            cons.abs() < 0.10,
            "{name}: growth not cleanly exponential: consistency {cons:.3}"
        );
    }
    let dev_c = (sig_c - sigma_inviscid) / sigma_inviscid;
    let dev_f = (sig_f - sigma_inviscid) / sigma_inviscid;
    // Physics direction: finite Re + lattice resolution bias the LBM
    // rate LOW; refining resolution AND Reynolds number must move it
    // TOWARD the inviscid oracle (the convergence claim, not a fixed
    // loose envelope). Measured: coarse -21%, fine (Re 432, b = 12)
    // closer — gates authored just outside the measured values.
    assert!(sig_c > 0.0 && sig_f > 0.0, "jet must be unstable");
    assert!(
        dev_c < 0.0 && dev_f < 0.0,
        "viscous bias must be LOW: {dev_c:.3}, {dev_f:.3}"
    );
    assert!(
        dev_f > dev_c + 0.02,
        "refinement must close on the inviscid rate: coarse {dev_c:.3}, fine {dev_f:.3}"
    );
    assert!(
        dev_c.abs() < 0.30,
        "coarse deviation {dev_c:.3} out of envelope"
    );
    assert!(
        dev_f.abs() < 0.20,
        "fine deviation {dev_f:.3} out of envelope"
    );
    println!(
        "{{\"suite\":\"fs-aeroac\",\"case\":\"bickley-lbm\",\"alpha\":{alpha:.4},\"sigma_inviscid\":{sigma_inviscid:.5},\"sigma_coarse\":{sig_c:.5},\"sigma_fine\":{sig_f:.5},\"dev_coarse\":{dev_c:.4},\"dev_fine\":{dev_f:.4},\"verdict\":\"pass\"}}"
    );
}
