//! Qualification of the hybrid compressible Navier-Stokes-Fourier grid (bead pexfp.134, M3.2-W3b).
//!
//! This file grows case by case. Landed so far: the exact linearized-NSF solver's limits (G0), the isothermal-mutation check (G0),
//! the declared fit-window refusals and `quiescent_conduction_fourier` (G1: entropy-mode decay against the exact linearized
//! rate on a three-rung ladder at a fixed physical wavelength, i.e. acoustic scaling with `dt` proportional to `dx`).
//! Not yet landed: the acoustic phase-speed/attenuation ladder at the campaign resolution, the manufactured solution, the
//! Couette case, the Mach-defect map and the typed Mach refusal. Nothing here is tuned; bands are proposals.
use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::dispersion::{linearized_modes, linearized_nsf_modes, GasPoint};
use fs_lbm::nsf3::{NsfGrid3, NsfSpec, WallThermal};
use fs_material::gas::{GasSpec, GasState};

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

fn point_of(g: &GasState, mu_bulk: f64) -> GasPoint {
    GasPoint {
        rho0: g.density,
        t0: g.temperature,
        gas_constant: g.specific_gas_constant,
        gamma: g.gamma,
        mu: g.dynamic_viscosity,
        mu_bulk,
        kappa: g.thermal_conductivity,
    }
}

fn spec(dims: [usize; 3], dx_m: f64) -> NsfSpec {
    NsfSpec {
        dims,
        periodic: [true; 3],
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 293.15,
        pressure_pa: 101_325.0,
        dx_m,
        temperature_window_k: [273.15, 333.15],
        maximum_mach: 0.3,
        maximum_temperature_ratio: 0.2,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "Cramer 2012, Phys. Fluids 24, 066102: mu_b = 0.73 mu for N2 at 293 K",
        maximum_cells: 10_000_000,
        wall: WallThermal::Adiabatic,
        workers: 1, // small grids: one worker is faster than 27 thread spawns per step, and the bits are worker-independent
        transport_points: 33,
    }
}

/// A decay fit window is declared in units of the e-folding time of the mode; shorter than this is refused.
const MIN_DECAY_WINDOW_EFOLDS: f64 = 0.25;
/// An oscillatory fit window is declared in periods; shorter than this is refused (bead pexfp.134).
const MIN_OSCILLATION_WINDOW_PERIODS: f64 = 3.0;

fn admit_decay_window(efolds: f64) -> Result<(), &'static str> {
    if efolds.is_finite() && efolds >= MIN_DECAY_WINDOW_EFOLDS { Ok(()) } else { Err("decay fit window shorter than the declared minimum") }
}
fn admit_oscillation_window(periods: f64) -> Result<(), &'static str> {
    if periods.is_finite() && periods >= MIN_OSCILLATION_WINDOW_PERIODS { Ok(()) } else { Err("oscillation fit window shorter than 3 periods") }
}

/// Least-squares slope of `ln|y|` against `t`.
fn log_slope(t: &[f64], y: &[f64]) -> f64 {
    let n = t.len() as f64;
    let l: Vec<f64> = y.iter().map(|v| v.abs().ln()).collect();
    let (st, sl) = (t.iter().sum::<f64>(), l.iter().sum::<f64>());
    let stl: f64 = t.iter().zip(&l).map(|(a, b)| a * b).sum();
    let stt: f64 = t.iter().map(|a| a * a).sum();
    (n * stl - st * sl) / (n * stt - st * st)
}

#[test]
fn dispersion_solver_limits() {
    // G0. Without transport the five modes are exact: acoustic +-i c k, entropy 0, shear 0.
    let g = air();
    let mut p = point_of(&g, 0.0);
    p.mu = 0.0;
    p.kappa = 0.0;
    let c = p.sound_speed();
    for k in [1.0, 18.3, 5236.0, 1.0e5] {
        let m = linearized_modes(&p, k).unwrap();
        assert!((m[0].im / (c * k) - 1.0).abs() < 1e-13 && m[0].re.abs() < 1e-9 * c * k, "k {k}: {:?}", m[0]);
        assert!((m[1].im / (c * k) + 1.0).abs() < 1e-13);
        assert!(m[2].re.abs() < 1e-9 * c * k && m[3].re == 0.0 && m[4].re == 0.0);
    }
    // Low-k limits with transport: the acoustic decay rate over c is the Stokes-Kirchhoff absorption (bulk viscosity 0 there),
    // the entropy decay is kappa k^2 / (rho c_p) and the shear decay is mu k^2 / rho.
    let p = point_of(&g, 0.0);
    let omega_ref = 2.0 * core::f64::consts::PI * 100.0;
    let k = omega_ref / p.sound_speed();
    let m = linearized_modes(&p, k).unwrap();
    let alpha_sk = g.stokes_kirchhoff_absorption(omega_ref);
    assert!(alpha_sk > 0.0, "an oracle that returns 0.0 must not pass");
    let alpha_exact = -m[0].re / (m[0].im / k);
    assert!((alpha_exact / alpha_sk - 1.0).abs() <= 1e-6, "attenuation {alpha_exact} vs Stokes-Kirchhoff {alpha_sk}");
    let entropy = p.kappa * k * k / (p.rho0 * p.specific_heat_cp());
    assert!((-m[2].re / entropy - 1.0).abs() <= 1e-6, "entropy decay {} vs {entropy}", -m[2].re);
    assert!((-m[3].re / (p.mu * k * k / p.rho0) - 1.0).abs() <= 1e-13);
    // refusals
    assert!(linearized_modes(&p, 0.0).is_err() && linearized_modes(&p, f64::NAN).is_err());
    let mut bad = p;
    bad.gamma = 1.0;
    assert!(linearized_modes(&bad, 1.0).is_err());
}

#[test]
fn isothermal_mutation_fails() {
    // G0. The solver's acoustic phase speed is sqrt(gamma R T), a factor sqrt(gamma) above the isothermal sqrt(R T): an
    // isothermal expectation must be rejected by the same comparison that accepts the adiabatic one.
    let g = air();
    let p = point_of(&g, 0.0);
    let k = 18.3;
    let measured = linearized_modes(&p, k).unwrap()[0].im / k;
    let within = |expected: f64| (measured / expected - 1.0).abs() <= 1e-4;
    assert!(within((p.gamma * p.gas_constant * p.t0).sqrt()));
    assert!(!within((p.gas_constant * p.t0).sqrt()));
}

#[test]
fn fit_window_refusals() {
    // 3 periods admitted, 2.999 refused; the decay window minimum likewise.
    assert!(admit_oscillation_window(3.0).is_ok() && admit_oscillation_window(2.999).is_err());
    assert!(admit_oscillation_window(f64::NAN).is_err());
    assert!(admit_decay_window(MIN_DECAY_WINDOW_EFOLDS).is_ok() && admit_decay_window(MIN_DECAY_WINDOW_EFOLDS * (1.0 - 1e-9)).is_err());
}

/// Entropy-mode decay on a periodic box of `ppw` cells per wavelength (fixed physical wavelength `lambda_m`): relative error of
/// the fitted decay rate against the exact linearized rate.
fn entropy_decay_error(ppw: usize, lambda_m: f64, transverse: usize) -> f64 {
    let a = 5.0e-4;
    let dx = lambda_m / ppw as f64;
    let k_lat = 2.0 * core::f64::consts::PI / ppw as f64;
    let mut g = NsfGrid3::<D3q27>::new(spec([ppw, transverse, transverse], dx), |_| false, move |x| {
        let t_rel = 1.0 + a * (k_lat * x[0] as f64).sin();
        (1.0 / t_rel, [0.0; 3], 293.15 * t_rel) // uniform pressure: the entropy mode
    })
    .unwrap();
    let law = g.law().clone();
    let k_m = 2.0 * core::f64::consts::PI / lambda_m;
    let exact = -linearized_nsf_modes(&law, k_m).unwrap()[2].re * law.dt_s; // decay per step
    let gamma = law.gamma;
    let efold = 1.0 / exact;
    let (t_start, t_end) = (0.02 * efold, 0.35 * efold);
    admit_decay_window((t_end - t_start) / efold).unwrap();
    let sample_every = ((efold / 200.0) as usize).max(1);
    let (mut ts, mut amps) = (Vec::new(), Vec::new());
    let mut step = 0usize;
    while (step as f64) < t_end {
        if (step as f64) >= t_start && step % sample_every == 0 {
            let mut acc = 0.0;
            for x in 0..ppw {
                let (rho, _u, t_k) = g.cell([x, 0, 0]);
                let s = (t_k - law.temperature_ref_k) / law.temperature_ref_k - (gamma - 1.0) * (rho - 1.0);
                acc += s * (k_lat * x as f64).sin();
            }
            ts.push(step as f64);
            amps.push(2.0 * acc / ppw as f64);
        }
        g.step().unwrap();
        step += 1;
    }
    let fitted = -log_slope(&ts, &amps);
    fitted / exact - 1.0
}

#[test]
fn transverse_extent_does_not_matter() {
    // G0. A field uniform across the box gives the same result on a 2 x 2 and a 4 x 4 transverse section (the 2 x 2 box is what the
    // ladders below use to stay affordable): same decay-rate error to 1e-9 relative on the 16-point rung.
    let a = entropy_decay_error(16, 1.2e-3, 2);
    let b = entropy_decay_error(16, 1.2e-3, 4);
    println!("transverse 2x2: {a:+.9e}, 4x4: {b:+.9e}");
    assert!((a - b).abs() <= 1e-9 * a.abs().max(1e-3), "{a} vs {b}");
}

#[test]
fn quiescent_conduction_fourier() {
    // G1. A small isobaric temperature wave in a quiescent periodic box decays as the entropy mode of the exact linearized NSF
    // (Fourier conduction at constant pressure, kappa k^2 / (rho c_p) to leading order), on a 16 / 32 / 64 points-per-wavelength
    // ladder at a fixed physical wavelength of 1.2 mm (2 x 2 transverse section, fit window 0.02 to 0.35 e-folds); observed order >= 1.8 and the finest rung within 1e-3.
    let lambda_m = 1.2e-3;
    let errs: Vec<f64> = [16usize, 32, 64].iter().map(|&n| entropy_decay_error(n, lambda_m, 2)).collect();
    for (n, e) in [16usize, 32, 64].iter().zip(&errs) {
        println!("conduction {n} points per wavelength: decay rate / exact - 1 = {e:+.4e}");
    }
    let order = (errs[1].abs() / errs[2].abs()).log2();
    println!("observed order between 32 and 64: {order:.2}");
    assert!(errs.windows(2).all(|w| w[0].abs() > w[1].abs()), "error must shrink: {errs:?}");
    assert!(order >= 1.8, "observed order {order}");
    assert!(errs[2].abs() <= 1e-3, "finest rung {}", errs[2]);
}
