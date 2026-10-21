//! Qualification of the hybrid compressible Navier-Stokes-Fourier grid (bead pexfp.134, M3.2-W3b).
//!
//! This file grows case by case. Landed so far: the exact linearized-NSF solver's limits and eigenvectors (G0), the
//! isothermal-mutation check (G0), the declared fit-window refusals, `quiescent_conduction_fourier` (G1: entropy-mode decay
//! against the exact linearized rate on a three-rung ladder at a fixed physical wavelength, i.e. acoustic scaling with `dt`
//! proportional to `dx`) and `linear_modes_match_linearized_nsf` (G1: phase speed and attenuation of a travelling acoustic wave at
//! the 10 kHz campaign resolution, and shear-mode decay), `couette_viscous_heating_profile` (G1: compressible Couette flow
//! between isothermal walls, one of them moving) and `finite_mach_and_lattice_defects_map` (G2: Galilean invariance of the
//! thermal and viscous channels at Mach 0.05, 0.1 and 0.2, which sets `MAX_QUALIFIED_MACH`). Not yet landed: the manufactured
//! solution. Nothing here is tuned; bands are proposals.
//!
//! The heavy tests (acoustic, Couette, Mach map) are `#[ignore]`d because a debug build is far too slow for them: run them as
//! `cargo test --release -p fs-lbm --test nsf3_qualification -- --include-ignored`.
use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::dispersion::{
    C64, GasPoint, linearized_modes, linearized_nsf_modes, longitudinal_eigenvector,
};
use fs_lbm::nsf3::{MAX_QUALIFIED_MACH, NsfGrid3, NsfSpec, THETA0, WallThermal};
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
        maximum_mach: MAX_QUALIFIED_MACH,
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
    if efolds.is_finite() && efolds >= MIN_DECAY_WINDOW_EFOLDS {
        Ok(())
    } else {
        Err("decay fit window shorter than the declared minimum")
    }
}
fn admit_oscillation_window(periods: f64) -> Result<(), &'static str> {
    if periods.is_finite() && periods >= MIN_OSCILLATION_WINDOW_PERIODS {
        Ok(())
    } else {
        Err("oscillation fit window shorter than 3 periods")
    }
}

/// Least-squares slope of `y` against `t`.
fn slope(t: &[f64], y: &[f64]) -> f64 {
    let n = t.len() as f64;
    let (st, sy) = (t.iter().sum::<f64>(), y.iter().sum::<f64>());
    let sty: f64 = t.iter().zip(y).map(|(a, b)| a * b).sum();
    let stt: f64 = t.iter().map(|a| a * a).sum();
    (n * sty - st * sy) / (n * stt - st * st)
}

/// Least-squares slope of `ln|y|` against `t`.
fn log_slope(t: &[f64], y: &[f64]) -> f64 {
    slope(t, &y.iter().map(|v| v.abs().ln()).collect::<Vec<_>>())
}

/// Weighted least-squares slope of `y` against `t`.
fn weighted_slope(t: &[f64], y: &[f64], w: &[f64]) -> f64 {
    let sw: f64 = w.iter().sum();
    let tb = t.iter().zip(w).map(|(a, b)| a * b).sum::<f64>() / sw;
    let yb = y.iter().zip(w).map(|(a, b)| a * b).sum::<f64>() / sw;
    let num: f64 = (0..t.len()).map(|i| w[i] * (t[i] - tb) * (y[i] - yb)).sum();
    let den: f64 = (0..t.len()).map(|i| w[i] * (t[i] - tb) * (t[i] - tb)).sum();
    num / den
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
        assert!(
            (m[0].im / (c * k) - 1.0).abs() < 1e-13 && m[0].re.abs() < 1e-9 * c * k,
            "k {k}: {:?}",
            m[0]
        );
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
    assert!(
        (alpha_exact / alpha_sk - 1.0).abs() <= 1e-6,
        "attenuation {alpha_exact} vs Stokes-Kirchhoff {alpha_sk}"
    );
    let entropy = p.kappa * k * k / (p.rho0 * p.specific_heat_cp());
    assert!(
        (-m[2].re / entropy - 1.0).abs() <= 1e-6,
        "entropy decay {} vs {entropy}",
        -m[2].re
    );
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
fn longitudinal_eigenvectors_satisfy_the_momentum_equation() {
    // G0. The eigenvector is built from continuity and energy; the momentum equation is the independent check that it belongs
    // to the root the cubic returned: s u = -i k R T0 (1 + t) - k^2 nu_L u in SI. Residual relative to the sum of the sizes of
    // the three terms, for the right- and left-moving acoustic waves and the entropy mode.
    let g = air();
    let mul = |a: C64, b: C64| C64 {
        re: a.re * b.re - a.im * b.im,
        im: a.re * b.im + a.im * b.re,
    };
    for (f_hz, ratio) in [
        (1.0e3, 0.73),
        (1.0e4, 0.73),
        (1.0e5, 0.0),
        (1.0e5, 0.73),
        (1.0e6, 0.73),
    ] {
        let p = point_of(&g, ratio * g.dynamic_viscosity);
        let c = p.sound_speed();
        let k = 2.0 * core::f64::consts::PI * f_hz / c;
        let modes = linearized_modes(&p, k).unwrap();
        let nu_l = (4.0 * p.mu / 3.0 + p.mu_bulk) / p.rho0;
        for (name, s) in [
            (
                "right",
                C64 {
                    re: modes[0].re,
                    im: -modes[0].im,
                },
            ),
            ("left", modes[0]),
            ("entropy", modes[2]),
        ] {
            let v = longitudinal_eigenvector(&p, k, s).unwrap();
            let su = mul(
                s,
                C64 {
                    re: c * v[1].re,
                    im: c * v[1].im,
                },
            );
            let rhs = C64 {
                re: k * p.gas_constant * p.t0 * v[2].im - k * k * nu_l * c * v[1].re,
                im: -k * p.gas_constant * p.t0 * (1.0 + v[2].re) - k * k * nu_l * c * v[1].im,
            };
            let scale = su.re.hypot(su.im)
                + k * p.gas_constant * p.t0 * (1.0 + v[2].re.hypot(v[2].im))
                + k * k * nu_l * c * v[1].re.hypot(v[1].im);
            let residual = (su.re - rhs.re).hypot(su.im - rhs.im) / scale;
            println!(
                "{name} eigenvector at {f_hz} Hz, bulk ratio {ratio}: momentum residual {residual:.3e}"
            );
            assert!(residual <= 1e-12, "{name} {f_hz} Hz: {residual}");
        }
        let right = longitudinal_eigenvector(
            &p,
            k,
            C64 {
                re: modes[0].re,
                im: -modes[0].im,
            },
        )
        .unwrap();
        assert!(
            right[1].im < 0.0 && right[1].im.abs() > 1e-9,
            "damping must show in the eigenvector at {f_hz} Hz: {:?}",
            right[1]
        );
    }
    // without transport the right-moving mode is the ideal isentropic wave
    let mut p = point_of(&g, 0.0);
    p.mu = 0.0;
    p.kappa = 0.0;
    let m = linearized_modes(&p, 18.3).unwrap()[0];
    let v = longitudinal_eigenvector(
        &p,
        18.3,
        C64 {
            re: m.re,
            im: -m.im,
        },
    )
    .unwrap();
    assert!(
        (v[1].re - 1.0).abs() < 1e-13
            && v[1].im.abs() < 1e-12
            && (v[2].re - (p.gamma - 1.0)).abs() < 1e-13
            && v[2].im.abs() < 1e-12
    );
}

#[test]
fn fit_window_refusals() {
    // 3 periods admitted, 2.999 refused; the decay window minimum likewise.
    assert!(admit_oscillation_window(3.0).is_ok() && admit_oscillation_window(2.999).is_err());
    assert!(admit_oscillation_window(f64::NAN).is_err());
    assert!(
        admit_decay_window(MIN_DECAY_WINDOW_EFOLDS).is_ok()
            && admit_decay_window(MIN_DECAY_WINDOW_EFOLDS * (1.0 - 1e-9)).is_err()
    );
}

/// Steps `g` through the declared decay window (0.02 to 0.35 e-folds of the exact rate `exact` per step), samples the mode
/// amplitude `probe` about 200 times, and returns the fitted decay rate over the exact one, minus 1.
fn fit_decay(g: &mut NsfGrid3<D3q27>, exact: f64, probe: impl Fn(&NsfGrid3<D3q27>) -> f64) -> f64 {
    let efold = 1.0 / exact;
    let (t_start, t_end) = (0.02 * efold, 0.35 * efold);
    admit_decay_window((t_end - t_start) / efold).unwrap();
    let sample_every = ((efold / 200.0) as usize).max(1);
    let (mut ts, mut amps) = (Vec::new(), Vec::new());
    let mut step = 0usize;
    while (step as f64) < t_end {
        if (step as f64) >= t_start && step % sample_every == 0 {
            ts.push(step as f64);
            amps.push(probe(g));
        }
        g.step().unwrap();
        step += 1;
    }
    -log_slope(&ts, &amps) / exact - 1.0
}

/// Entropy-mode decay on a periodic box of `ppw` cells per wavelength (fixed physical wavelength `lambda_m`): relative error of
/// the fitted decay rate against the exact linearized rate.
fn entropy_decay_error(ppw: usize, lambda_m: f64, transverse: usize) -> f64 {
    let a = 5.0e-4;
    let dx = lambda_m / ppw as f64;
    let k_lat = 2.0 * core::f64::consts::PI / ppw as f64;
    let mut g = NsfGrid3::<D3q27>::new(
        spec([ppw, transverse, transverse], dx),
        |_| false,
        move |x| {
            let t_rel = 1.0 + a * (k_lat * x[0] as f64).sin();
            (1.0 / t_rel, [0.0; 3], 293.15 * t_rel) // uniform pressure: the entropy mode
        },
    )
    .unwrap();
    let law = g.law().clone();
    let k_m = 2.0 * core::f64::consts::PI / lambda_m;
    let exact = -linearized_nsf_modes(&law, k_m).unwrap()[2].re * law.dt_s; // decay per step
    let gamma = law.gamma;
    fit_decay(&mut g, exact, |g| {
        let mut acc = 0.0;
        for x in 0..ppw {
            let (rho, _u, t_k) = g.cell([x, 0, 0]);
            let s =
                (t_k - law.temperature_ref_k) / law.temperature_ref_k - (gamma - 1.0) * (rho - 1.0);
            acc += s * (k_lat * x as f64).sin();
        }
        2.0 * acc / ppw as f64
    })
}

/// Shear-mode decay: a transverse velocity wave `u_y = a sin(k x)` in a uniform gas decays at `mu k^2 / rho`; relative error of
/// the fitted rate against the exact linearized one, on the same window as the entropy mode.
fn shear_decay_error(ppw: usize, lambda_m: f64, transverse: usize) -> f64 {
    let a = 1.0e-3; // lattice velocity amplitude: Mach about 1.5e-3
    let dx = lambda_m / ppw as f64;
    let k_lat = 2.0 * core::f64::consts::PI / ppw as f64;
    let mut g = NsfGrid3::<D3q27>::new(
        spec([ppw, transverse, transverse], dx),
        |_| false,
        move |x| (1.0, [0.0, a * (k_lat * x[0] as f64).sin(), 0.0], 293.15),
    )
    .unwrap();
    let k_m = 2.0 * core::f64::consts::PI / lambda_m;
    let exact = -linearized_nsf_modes(g.law(), k_m).unwrap()[3].re * g.law().dt_s;
    fit_decay(&mut g, exact, |g| {
        let acc: f64 = (0..ppw)
            .map(|x| g.cell([x, 0, 0]).1[1] * (k_lat * x as f64).sin())
            .sum();
        2.0 * acc / ppw as f64
    })
}

/// Fitted and exact rates of a travelling acoustic wave, per second.
struct AcousticFit {
    omega_fit: f64,
    sigma_fit: f64,
    omega_exact: f64,
    sigma_exact: f64,
    /// Stokes-Kirchhoff amplitude decay rate `alpha_SK(omega) c`, the scale the excess attenuation is judged against.
    sk_rate: f64,
    /// Mean amplitudes of the left-moving and entropy modes over the window, relative to the right-moving one.
    leak_left: f64,
    leak_entropy: f64,
}

type Cx = (f64, f64);
fn cmul(a: Cx, b: Cx) -> Cx {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}
fn csub(a: Cx, b: Cx) -> Cx {
    (a.0 - b.0, a.1 - b.1)
}
fn cdiv(a: Cx, b: Cx) -> Cx {
    let d = b.0 * b.0 + b.1 * b.1;
    ((a.0 * b.0 + a.1 * b.1) / d, (a.1 * b.0 - a.0 * b.1) / d)
}
fn det3(m: [[Cx; 3]; 3]) -> Cx {
    let minor = |a: Cx, b: Cx, c: Cx, d: Cx| csub(cmul(a, d), cmul(b, c));
    let t0 = cmul(m[0][0], minor(m[1][1], m[1][2], m[2][1], m[2][2]));
    let t1 = cmul(m[0][1], minor(m[1][0], m[1][2], m[2][0], m[2][2]));
    let t2 = cmul(m[0][2], minor(m[1][0], m[1][1], m[2][0], m[2][1]));
    (t0.0 - t1.0 + t2.0, t0.1 - t1.1 + t2.1)
}
/// Solves `V a = x` where `v[c]` is column `c` of `V`, by Cramer's rule (a 3 x 3 complex system, well conditioned here: the three
/// modes are far from parallel).
fn solve3(v: [[Cx; 3]; 3], x: [Cx; 3]) -> [Cx; 3] {
    let matrix = |replace: Option<usize>| -> [[Cx; 3]; 3] {
        let mut m = [[(0.0, 0.0); 3]; 3];
        for (row, m_row) in m.iter_mut().enumerate() {
            for (c, entry) in m_row.iter_mut().enumerate() {
                *entry = if replace == Some(c) {
                    x[row]
                } else {
                    v[c][row]
                };
            }
        }
        m
    };
    let base = det3(matrix(None));
    [
        cdiv(det3(matrix(Some(0))), base),
        cdiv(det3(matrix(Some(1))), base),
        cdiv(det3(matrix(Some(2))), base),
    ]
}

/// A right-moving wave of relative density amplitude `amplitude` (the exact eigenvector of the grid's own gas point) and one
/// wavelength `lambda_m` on `ppw` cells, run for `periods` periods after one settling period. At 16 samples per period the
/// complex amplitudes of the fundamental of density, velocity and temperature are projected on the three exact longitudinal
/// eigenvectors (right-moving, left-moving, entropy), so that the small amounts of the other two modes that the lattice
/// initialisation excites cannot beat against the wave. The frequency is the slope of the phase of the right-moving amplitude
/// against the exact phase advance and the decay rate the slope of its `ln |a|`, both by Hann-weighted least squares over the window.
fn acoustic_fit(
    ppw: usize,
    lambda_m: f64,
    amplitude: f64,
    periods: f64,
    transverse: usize,
) -> AcousticFit {
    let gas = air();
    let gamma = gas.gamma;
    let dx = lambda_m / ppw as f64;
    let k_lat = 2.0 * core::f64::consts::PI / ppw as f64;
    let k_m = 2.0 * core::f64::consts::PI / lambda_m;
    let u_unit = (gamma / 3.0).sqrt(); // c / (dx / dt)
    let law = NsfGrid3::<D3q27>::new(
        spec([ppw, transverse, transverse], dx),
        |_| false,
        |_| (1.0, [0.0; 3], 293.15),
    )
    .unwrap()
    .law()
    .clone();
    let point = GasPoint::from_law(&law);
    let mode = linearized_nsf_modes(&law, k_m).unwrap();
    let rates = [
        C64 {
            re: mode[0].re,
            im: -mode[0].im,
        },
        mode[0],
        mode[2],
    ]; // right, left, entropy
    let vec_of = |s: C64| longitudinal_eigenvector(&point, k_m, s).unwrap();
    let vs = [vec_of(rates[0]), vec_of(rates[1]), vec_of(rates[2])];
    let v = vs[0];
    let mut g = NsfGrid3::<D3q27>::new(
        spec([ppw, transverse, transverse], dx),
        |_| false,
        move |x| {
            let phase = k_lat * x[0] as f64;
            let re = |z: C64| z.re * phase.cos() - z.im * phase.sin();
            (
                1.0 + amplitude * re(v[0]),
                [amplitude * u_unit * re(v[1]), 0.0, 0.0],
                293.15 * (1.0 + amplitude * re(v[2])),
            )
        },
    )
    .unwrap();
    let dt = law.dt_s;
    let omega_step = mode[0].im * dt;
    let period_steps = 2.0 * core::f64::consts::PI / omega_step;
    let (t_start, t_end) = (period_steps, period_steps * (1.0 + periods));
    admit_oscillation_window((t_end - t_start) / period_steps).unwrap();
    let sample_every = ((period_steps / 16.0) as usize).max(1);
    // columns of V are the eigenvectors: V[j][row]
    let col = |w: [C64; 3]| [(w[0].re, w[0].im), (w[1].re, w[1].im), (w[2].re, w[2].im)];
    let basis = [col(vs[0]), col(vs[1]), col(vs[2])];
    let (mut ts, mut phases, mut logs) = (Vec::new(), Vec::new(), Vec::new());
    let (mut left, mut entropy, mut count) = (0.0, 0.0, 0.0);
    let mut step = 0usize;
    while (step as f64) < t_end {
        if (step as f64) >= t_start && step % sample_every == 0 {
            let mut x = [(0.0, 0.0); 3];
            for i in 0..ppw {
                let (rho, u, t_k) = g.cell([i, 0, 0]);
                let f = [rho - 1.0, u[0] / u_unit, t_k / law.temperature_ref_k - 1.0];
                for (r, fr) in f.iter().enumerate() {
                    x[r].0 += fr * (k_lat * i as f64).cos();
                    x[r].1 -= fr * (k_lat * i as f64).sin(); // (2/N) sum f exp(-i k x)
                }
            }
            for xr in &mut x {
                *xr = (2.0 * xr.0 / ppw as f64, 2.0 * xr.1 / ppw as f64);
            }
            let a = solve3(basis, x);
            let theta = omega_step * step as f64;
            ts.push(step as f64);
            phases.push(
                (a[0].1 * theta.cos() + a[0].0 * theta.sin())
                    .atan2(a[0].0 * theta.cos() - a[0].1 * theta.sin()),
            ); // arg(a exp(i theta))
            logs.push(a[0].0.hypot(a[0].1).ln());
            left += a[1].0.hypot(a[1].1) / a[0].0.hypot(a[0].1);
            entropy += a[2].0.hypot(a[2].1) / a[0].0.hypot(a[0].1);
            count += 1.0;
        }
        g.step().unwrap();
        step += 1;
    }
    // Hann weights over the declared window: a residual oscillation of the amplitude cannot bias the slopes through the window ends
    let weights: Vec<f64> = ts
        .iter()
        .map(|&s| {
            (core::f64::consts::PI * (s - t_start) / (t_end - t_start))
                .sin()
                .powi(2)
        })
        .collect();
    let omega_fit = (omega_step - weighted_slope(&ts, &phases, &weights)) / dt;
    let sigma_fit = -weighted_slope(&ts, &logs, &weights) / dt;
    let omega_exact = mode[0].im;
    let alpha_sk = gas.stokes_kirchhoff_absorption(omega_exact);
    assert!(alpha_sk > 0.0, "an oracle that returns 0.0 must not pass");
    AcousticFit {
        omega_fit,
        sigma_fit,
        omega_exact,
        sigma_exact: -mode[0].re,
        sk_rate: alpha_sk * omega_exact / k_m,
        leak_left: left / count,
        leak_entropy: entropy / count,
    }
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
    let errs: Vec<f64> = [16usize, 32, 64]
        .iter()
        .map(|&n| entropy_decay_error(n, lambda_m, 2))
        .collect();
    for (n, e) in [16usize, 32, 64].iter().zip(&errs) {
        println!("conduction {n} points per wavelength: decay rate / exact - 1 = {e:+.4e}");
    }
    let order = (errs[1].abs() / errs[2].abs()).log2();
    println!("observed order between 32 and 64: {order:.2}");
    assert!(
        errs.windows(2).all(|w| w[0].abs() > w[1].abs()),
        "error must shrink: {errs:?}"
    );
    assert!(order >= 1.8, "observed order {order}");
    assert!(errs[2].abs() <= 1e-3, "finest rung {}", errs[2]);
}

/// Wavelength of the acoustic ladder: 34.4 mm = 688 cells of 50 um (the campaign grid), i.e. 9.98 kHz in air at 293.15 K; the
/// declared upper band edge is 10 kHz and the nearest whole-cell wavelength is used.
const ACOUSTIC_LAMBDA_M: f64 = 34.4e-3;
/// Relative density amplitude of the acoustic ladder. The fitted attenuation moves with the square of the amplitude (harmonic
/// generation: +1.5e6 a^2 of the Stokes-Kirchhoff rate, measured over 1e-7 to 1e-5), so the amplitude is chosen to put that
/// shift two orders below the finest-rung excess; round-off sets the lower limit.
const ACOUSTIC_AMPLITUDE: f64 = 4.0e-7;

#[test]
#[ignore = "heavy: run with --release -- --ignored (about 20 s in release)"]
fn acoustic_excess_attenuation_is_amplitude_independent() {
    // G0. The excess attenuation of the ladder must not be a nonlinear or round-off artefact: on the 172-cell rung it agrees
    // between the ladder amplitude and five times it to 2e-5 of the Stokes-Kirchhoff rate, while 25 times it (1e-5) moves it by
    // 1.5e-4, i.e. the check can fail.
    let excess = |a: f64| {
        let f = acoustic_fit(172, ACOUSTIC_LAMBDA_M, a, 8.0, 2);
        (f.sigma_fit - f.sigma_exact) / f.sk_rate
    };
    let (low, five, large) = (
        excess(ACOUSTIC_AMPLITUDE),
        excess(5.0 * ACOUSTIC_AMPLITUDE),
        excess(25.0 * ACOUSTIC_AMPLITUDE),
    );
    println!(
        "excess/SK at 172 cells: a {ACOUSTIC_AMPLITUDE:e}: {low:+.4e}, 5a: {five:+.4e}, 25a: {large:+.4e}"
    );
    assert!((low - five).abs() <= 2e-5, "{low} vs {five}");
    assert!(
        (low - large).abs() > 5e-5,
        "the check must be able to see the nonlinear shift: {low} vs {large}"
    );
}

#[test]
#[ignore = "heavy: campaign-resolution rung, run with --release -- --ignored (about 150 s in release on a loaded host)"]
fn linear_modes_match_linearized_nsf() {
    // G1 (a). Phase speed and attenuation of a travelling acoustic wave at 10 kHz on 172 / 344 / 688 cells per wavelength
    // (dx 200 / 100 / 50 um, acoustic scaling), against the exact linearized-NSF mode with the declared bulk viscosity; the
    // excess attenuation is judged against the Stokes-Kirchhoff rate. Shear-mode decay on a 16 / 32 / 64 ladder at 1.2 mm. The
    // entropy mode is `quiescent_conduction_fourier`.
    let mut phase = Vec::new();
    let mut excess = Vec::new();
    for ppw in [172usize, 344, 688] {
        let f = acoustic_fit(ppw, ACOUSTIC_LAMBDA_M, ACOUSTIC_AMPLITUDE, 8.0, 2);
        let (pe, ex) = (
            f.omega_fit / f.omega_exact - 1.0,
            (f.sigma_fit - f.sigma_exact) / f.sk_rate,
        );
        println!(
            "acoustic {ppw} points per wavelength: phase speed error {pe:+.4e}, attenuation {:.6e} /s vs exact {:.6e} /s (Stokes-Kirchhoff {:.5e} /s), excess / SK {ex:+.4e}; other modes excited: left {:.2e}, entropy {:.2e}",
            f.sigma_fit, f.sigma_exact, f.sk_rate, f.leak_left, f.leak_entropy
        );
        phase.push(pe);
        excess.push(ex);
    }
    let shear: Vec<f64> = [16usize, 32, 64]
        .iter()
        .map(|&n| shear_decay_error(n, 1.2e-3, 2))
        .collect();
    for (n, e) in [16usize, 32, 64].iter().zip(&shear) {
        println!("shear {n} points per wavelength: decay rate / exact - 1 = {e:+.4e}");
    }
    let order = |v: &[f64]| (v[1].abs() / v[2].abs()).log2();
    println!(
        "observed orders (finest pair): phase {:.2}, excess attenuation {:.2}, shear {:.2}",
        order(&phase),
        order(&excess),
        order(&shear)
    );
    assert!(
        phase[2].abs() <= 6e-5,
        "phase speed error at the campaign resolution {}",
        phase[2]
    );
    assert!(
        excess[2].abs() <= 0.10,
        "excess attenuation / Stokes-Kirchhoff {}",
        excess[2]
    );
    assert!(
        phase.windows(2).all(|w| w[0].abs() > w[1].abs()),
        "phase error must shrink: {phase:?}"
    );
    assert!(
        excess.windows(2).all(|w| w[0].abs() > w[1].abs()),
        "attenuation error must shrink: {excess:?}"
    );
    assert!(order(&phase) >= 1.8, "phase order {}", order(&phase));
    assert!(
        order(&excess) >= 1.8,
        "attenuation order {}",
        order(&excess)
    );
    assert!(
        shear.windows(2).all(|w| w[0].abs() > w[1].abs()),
        "shear error must shrink: {shear:?}"
    );
    assert!(
        order(&shear) >= 1.8 && shear[2].abs() <= 1e-3,
        "shear: order {} finest {}",
        order(&shear),
        shear[2]
    );
}

/// The two diffusive modes whose Galilean invariance is mapped.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Diffusive {
    /// Isobaric temperature wave (entropy mode, thermal conduction).
    Entropy,
    /// Transverse velocity wave (shear mode, viscosity).
    Shear,
}

/// Defects of a diffusive wave advected along x at Mach `mach` (uniform velocity, uniform pressure) against the exact solution,
/// which is the quiescent mode translated at the flow speed: the fitted decay rate over the quiescent exact rate, minus 1, and the
/// fitted translation speed over `0.999 * mach * c` (the mean flow actually imposed), minus 1 (zero for `mach = 0`). Same window as
/// the quiescent cases.
fn advected_wave(
    ppw: usize,
    lambda_m: f64,
    mach: f64,
    mode: Diffusive,
    transverse: usize,
) -> (f64, f64) {
    let a = match mode {
        Diffusive::Entropy => 5.0e-4,
        Diffusive::Shear => 1.0e-3,
    };
    let dx = lambda_m / ppw as f64;
    let k_lat = 2.0 * core::f64::consts::PI / ppw as f64;
    let gas = air();
    // c / (dx / dt) = sqrt(gamma / 3). The mean flow is 0.1 % below the nominal Mach number so that the local Mach number, which
    // follows the wave's own temperature and velocity fluctuations, stays within the declared bound at the top of the map.
    let u_lat = mach * (1.0 - 1.0e-3) * (gas.gamma / 3.0).sqrt();
    let mut g = NsfGrid3::<D3q27>::new(
        spec([ppw, transverse, transverse], dx),
        |_| false,
        move |x| {
            let w = (k_lat * x[0] as f64).sin();
            match mode {
                Diffusive::Entropy => (
                    1.0 / (1.0 + a * w),
                    [u_lat, 0.0, 0.0],
                    293.15 * (1.0 + a * w),
                ),
                Diffusive::Shear => (1.0, [u_lat, a * w, 0.0], 293.15),
            }
        },
    )
    .unwrap();
    let law = g.law().clone();
    let k_m = 2.0 * core::f64::consts::PI / lambda_m;
    let index = if mode == Diffusive::Entropy { 2 } else { 3 };
    let exact = -linearized_nsf_modes(&law, k_m).unwrap()[index].re * law.dt_s;
    let efold = 1.0 / exact;
    let (t_start, t_end) = (0.02 * efold, 0.35 * efold);
    admit_decay_window((t_end - t_start) / efold).unwrap();
    let sample_every = ((efold / 200.0) as usize).max(1);
    let (mut ts, mut logs, mut phases) = (Vec::new(), Vec::new(), Vec::new());
    let mut step = 0usize;
    while (step as f64) < t_end {
        if (step as f64) >= t_start && step % sample_every == 0 {
            let (mut re, mut im) = (0.0, 0.0);
            for x in 0..ppw {
                let (rho, u, t_k) = g.cell([x, 0, 0]);
                let s = match mode {
                    Diffusive::Entropy => {
                        (t_k - law.temperature_ref_k) / law.temperature_ref_k
                            - (law.gamma - 1.0) * (rho - 1.0)
                    }
                    Diffusive::Shear => u[1],
                };
                re += s * (k_lat * x as f64).cos();
                im += s * (k_lat * x as f64).sin();
            }
            let (re, im) = (2.0 * re / ppw as f64, 2.0 * im / ppw as f64); // A = re - i im
            let theta = k_lat * u_lat * step as f64;
            ts.push(step as f64);
            logs.push(re.hypot(im));
            // arg(A exp(i theta)): the phase of the translated wave relative to the exact translation
            phases.push(
                (re * theta.sin() - im * theta.cos()).atan2(re * theta.cos() + im * theta.sin()),
            );
        }
        g.step().unwrap();
        step += 1;
    }
    let decay = -log_slope(&ts, &logs) / exact - 1.0;
    let speed = if mach == 0.0 {
        0.0
    } else {
        slope(&ts, &phases) / (k_lat * u_lat)
    };
    (decay, speed)
}

/// Compressible Couette flow between isothermal walls at the reference temperature, the bottom wall at rest and the top wall
/// moving at `mach` times the sound speed: `n` cells across a gap of `gap_m`. The grid starts from the exact steady profiles
/// (linear velocity, `T = T_w + mu U^2 y (h - y) / (2 kappa h^2)`, uniform pressure) and runs `thermal_times` thermal
/// relaxation times `h^2 / (pi^2 alpha)` so that the lattice settles onto its own steady state. Returns the L-infinity error of
/// the temperature over the centre-line rise `mu U^2 / (8 kappa)` and of the velocity over `U`.
fn couette_errors(n: usize, gap_m: f64, mach: f64, thermal_times: f64) -> (f64, f64) {
    let gas = air();
    let dx = gap_m / n as f64;
    let wall_speed = mach * gas.sound_speed;
    let (mu, kappa) = (gas.dynamic_viscosity, gas.thermal_conductivity);
    let heating =
        |y: f64| mu * wall_speed * wall_speed * y * (gap_m - y) / (2.0 * kappa * gap_m * gap_m);
    let mut s = spec([2, n, 2], dx);
    s.periodic = [true, false, true];
    s.wall = WallThermal::Isothermal(293.15);
    let u_unit = (gas.gamma / 3.0).sqrt() / gas.sound_speed; // lattice velocity per m/s
    let mut g = NsfGrid3::<D3q27>::new(
        s,
        |_| false,
        move |x| {
            let y = (x[1] as f64 + 0.5) * dx;
            let t_k = 293.15 + heating(y);
            (
                293.15 / t_k,
                [wall_speed * y / gap_m * u_unit, 0.0, 0.0],
                t_k,
            )
        },
    )
    .unwrap();
    g.set_edge_velocity(1, true, [wall_speed * u_unit, 0.0, 0.0])
        .unwrap();
    let law = g.law().clone();
    let alpha = kappa / (gas.density * gas.specific_heat_cp);
    let tau = gap_m * gap_m / (core::f64::consts::PI.powi(2) * alpha);
    let steps = (thermal_times * tau / law.dt_s).ceil() as usize;
    for _ in 0..steps {
        g.step().unwrap();
    }
    let rise = mu * wall_speed * wall_speed / (8.0 * kappa);
    let (mut t_err, mut u_err) = (0.0f64, 0.0f64);
    for j in 0..n {
        let y = (j as f64 + 0.5) * dx;
        let (_rho, u, t_k) = g.cell([0, j, 0]);
        t_err = t_err.max((t_k - 293.15 - heating(y)).abs() / rise);
        u_err = u_err.max((u[0] / u_unit - wall_speed * y / gap_m).abs() / wall_speed);
    }
    (t_err, u_err)
}

#[test]
#[ignore = "heavy: run with --release -- --ignored (several minutes in release on a loaded host)"]
fn couette_viscous_heating_profile() {
    // G1. Compressible Couette flow between isothermal walls at the reference temperature, the top wall moving at Mach 0.05 (centre
    // line 26 mK warmer), gap 0.2 mm, 20 / 40 / 80 cells across (acoustic scaling), started from the exact profiles and run for six
    // thermal times: the temperature must follow T(y) = T_w + mu U^2 y (h - y) / (2 kappa h^2). The error is the L-infinity norm
    // over the centre-line rise mu U^2 / (8 kappa); a gas that is not heated at all would score 1. The reference holds the gas
    // properties constant, so a floor of order (dln kappa / dln T) * (rise / T) is expected at larger Mach numbers; the Mach
    // dependence on the coarse rung is printed, not asserted.
    let rungs = [20usize, 40, 80];
    let errs: Vec<(f64, f64)> = rungs
        .iter()
        .map(|&n| couette_errors(n, 2.0e-4, 0.05, 6.0))
        .collect();
    for (n, (t, u)) in rungs.iter().zip(&errs) {
        println!("couette {n} cells: temperature L-inf / rise {t:.4e}, velocity L-inf / U {u:.4e}");
    }
    for ma in [0.1, MAX_QUALIFIED_MACH] {
        println!(
            "couette 20 cells at Ma {ma}: temperature L-inf / rise {:.4e}",
            couette_errors(20, 2.0e-4, ma, 6.0).0
        );
    }
    let t: Vec<f64> = errs.iter().map(|e| e.0).collect();
    let order = (t[1] / t[2]).log2();
    println!("observed order between 40 and 80 cells: {order:.2}");
    assert!(
        t.windows(2).all(|w| w[0] > w[1]),
        "error must shrink: {t:?}"
    );
    assert!(order >= 1.8, "observed order {order}");
    assert!(t[2] <= 1e-3, "finest rung {}", t[2]);
    assert!(
        errs.iter().all(|e| e.1 <= 1e-3),
        "velocity profile {errs:?}"
    );
}

/// Mach numbers of the defect map; the translation speed of a wave is undefined at rest, so its defect is taken against the run
/// at `MACH_REFERENCE`.
const MACH_MAP: [f64; 3] = [0.05, 0.1, 0.2];
const MACH_REFERENCE: f64 = 0.025;
/// Proposed bands on the Mach-dependent part of the decay rate (relative to the quiescent rate) and of the translation speed.
const GALILEAN_DECAY_BAND: f64 = 5e-3;
const GALILEAN_SPEED_BAND: f64 = 5e-4;

/// Scaling exponent of a defect with Mach number: the least-squares slope of `ln |d|` against `ln Ma` over `MACH_MAP`.
fn mach_exponent(defects: &[f64]) -> f64 {
    let x: Vec<f64> = MACH_MAP.iter().map(|m| m.ln()).collect();
    slope(
        &x,
        &defects.iter().map(|d| d.abs().ln()).collect::<Vec<_>>(),
    )
}

#[test]
#[ignore = "heavy: run with --release -- --ignored (a few minutes in release on a loaded host)"]
fn finite_mach_and_lattice_defects_map() {
    // G2. Galilean invariance of the two diffusive channels: a thermal (entropy) wave and a shear wave advected at Mach 0.05,
    // 0.1 and 0.2 must decay at the quiescent rate and travel at exactly M c. The defect is the Mach-dependent part of the error
    // at 32 cells per wavelength (the decay error minus the quiescent one; the speed error minus that at Mach 0.025). Every
    // defect must stay within its band over the map; the scaling exponent of each is printed, as is the 64-cell value at the
    // largest Mach number. `MAX_QUALIFIED_MACH` is the largest Mach number of the map and the formulation refuses a larger bound.
    assert_eq!(*MACH_MAP.last().unwrap(), MAX_QUALIFIED_MACH);
    let (ppw, lambda_m) = (32usize, 1.2e-3);
    for mode in [Diffusive::Entropy, Diffusive::Shear] {
        let base = advected_wave(ppw, lambda_m, 0.0, mode, 2).0;
        let reference_speed = advected_wave(ppw, lambda_m, MACH_REFERENCE, mode, 2).1;
        let (mut decay, mut speed) = (Vec::new(), Vec::new());
        for ma in MACH_MAP {
            let (d, s) = advected_wave(ppw, lambda_m, ma, mode, 2);
            println!(
                "{mode:?} Ma {ma}: decay defect {:+.4e}, speed defect {:+.4e}",
                d - base,
                s - reference_speed
            );
            decay.push(d - base);
            speed.push(s - reference_speed);
        }
        println!(
            "{mode:?}: scaling exponent with Mach number: decay {:.2}, speed {:.2}",
            mach_exponent(&decay),
            mach_exponent(&speed)
        );
        let fine_base = advected_wave(2 * ppw, lambda_m, 0.0, mode, 2).0;
        let fine = advected_wave(2 * ppw, lambda_m, MAX_QUALIFIED_MACH, mode, 2).0 - fine_base;
        println!(
            "{mode:?}: decay defect at Ma {MAX_QUALIFIED_MACH}: {:+.4e} on 32 cells, {fine:+.4e} on 64 cells",
            decay[2]
        );
        assert!(
            decay.iter().all(|d| d.abs() <= GALILEAN_DECAY_BAND),
            "{mode:?} decay defects {decay:?}"
        );
        assert!(
            speed.iter().all(|d| d.abs() <= GALILEAN_SPEED_BAND),
            "{mode:?} speed defects {speed:?}"
        );
        assert!(
            decay.windows(2).all(|w| w[0].abs() < w[1].abs()),
            "{mode:?} decay defect must grow with Mach number: {decay:?}"
        );
    }
}

#[test]
fn edge_velocity_refusals() {
    // G0. A wall velocity must sit on a non-periodic axis, be finite and tangential, and respect the Mach bound.
    let mut s = spec([2, 4, 2], 5.0e-5);
    s.periodic = [true, false, true];
    let mut g = NsfGrid3::<D3q27>::new(s, |_| false, |_| (1.0, [0.0; 3], 293.15)).unwrap();
    let c = (g.law().gamma * THETA0).sqrt();
    assert!(
        g.set_edge_velocity(0, true, [0.0; 3]).is_err(),
        "periodic axis"
    );
    assert!(
        g.set_edge_velocity(3, true, [0.0; 3]).is_err(),
        "no such axis"
    );
    assert!(
        g.set_edge_velocity(1, true, [0.0, 0.01, 0.0]).is_err(),
        "wall-normal velocity"
    );
    assert!(
        g.set_edge_velocity(1, true, [f64::NAN, 0.0, 0.0]).is_err(),
        "non-finite velocity"
    );
    assert!(
        g.set_edge_velocity(1, true, [MAX_QUALIFIED_MACH * c * (1.0 + 1e-9), 0.0, 0.0])
            .is_err(),
        "above the Mach bound"
    );
    assert!(
        g.set_edge_velocity(1, true, [MAX_QUALIFIED_MACH * c * (1.0 - 1e-12), 0.0, 0.0])
            .is_ok()
    );
    assert!(
        g.set_edge_velocity(1, false, [0.0, 0.0, 0.1 * c]).is_ok(),
        "the other wall may move along z"
    );
    g.set_edge_velocity(1, true, [0.05 * c, 0.0, 0.0]).unwrap();
    g.step().unwrap();
}
