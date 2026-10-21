//! Battery of the hybrid compressible Navier-Stokes-Fourier grid (bead pexfp.133, M3.2-W3a).
//!
//! Gauntlet tiers are named per test. The heavy rungs (1e4-step closed box, the 192-point acoustic rung) are
//! release-lane runs; every number asserted here is the bead's band or an observed order, never a tuned value.
use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{equilibrium, NsfError, NsfGrid3, NsfSpec, WallThermal, THETA0};
use fs_material::gas::GasSpec;

fn spec(dims: [usize; 3], periodic: [bool; 3], wall: WallThermal, dx_m: f64) -> NsfSpec {
    NsfSpec {
        dims,
        periodic,
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 293.15,
        pressure_pa: 101_325.0,
        dx_m,
        temperature_window_k: [273.15, 333.15],
        maximum_mach: 0.3,
        maximum_temperature_ratio: 0.2,
        bulk_rate: 1.98,
        maximum_cells: 10_000_000,
        wall,
        workers: 4,
        transport_points: 33,
    }
}

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
    }
}

#[test]
fn equilibrium_recovers_conserved_moments() {
    // G0. Mass, momentum, the full second moment and every third moment the D3Q27 lattice can represent
    // (xyz and the six xxy-type components; the cubic xxx components alias onto the first moment, c^3 = c).
    let mut rng = Lcg(11);
    for _ in 0..200 {
        let rho = 1.0 + 0.1 * rng.next();
        let u = [0.05 * rng.next(), 0.05 * rng.next(), 0.05 * rng.next()];
        let theta = THETA0 * (1.0 + 0.1 * rng.next());
        let f = equilibrium::<D3q27>(rho, u, theta);
        let e = <D3q27 as fs_lbm::lattice3::Lattice3<27>>::E;
        let mut m0 = 0.0;
        let mut m1 = [0.0; 3];
        let mut m2 = [[0.0; 3]; 3];
        let mut xyz = 0.0;
        let mut xxy = [[0.0; 3]; 3];
        for q in 0..27 {
            let c = [f64::from(e[q][0]), f64::from(e[q][1]), f64::from(e[q][2])];
            m0 += f[q];
            for a in 0..3 {
                m1[a] += c[a] * f[q];
                for b in 0..3 {
                    m2[a][b] += c[a] * c[b] * f[q];
                    if a != b {
                        xxy[a][b] += c[a] * c[a] * c[b] * f[q];
                    }
                }
            }
            xyz += c[0] * c[1] * c[2] * f[q];
        }
        assert!((m0 - rho).abs() < 1e-13, "mass {m0} vs {rho}");
        for a in 0..3 {
            assert!((m1[a] - rho * u[a]).abs() < 1e-13, "momentum");
            for b in 0..3 {
                let want = rho * (u[a] * u[b] + if a == b { theta } else { 0.0 });
                assert!((m2[a][b] - want).abs() < 1e-13, "second moment {a}{b}: {} vs {want}", m2[a][b]);
                if a != b {
                    let w3 = rho * (u[a] * u[a] * u[b] + theta * u[b]);
                    assert!((xxy[a][b] - w3).abs() < 1e-13, "third moment {a}{a}{b}: {} vs {w3}", xxy[a][b]);
                }
            }
        }
        assert!((xyz - rho * u[0] * u[1] * u[2]).abs() < 1e-13, "xyz moment");
    }
}

#[test]
fn refusals_cap_and_cap_plus_one() {
    let base = || spec([4, 4, 4], [true; 3], WallThermal::Adiabatic, 5.0e-5);
    let rest = |_x: [usize; 3]| (1.0, [0.0; 3], 293.15);
    let free = |_x: [usize; 3]| false;
    // gas receipt window: 50 K admitted, 49.999 K refused
    let mut s = base();
    s.temperature_k = 50.0;
    s.temperature_window_k = [50.0, 60.0];
    assert!(NsfGrid3::<D3q27>::new(s, free, |_| (1.0, [0.0; 3], 50.0)).is_ok());
    let mut s = base();
    s.temperature_k = 49.999;
    s.temperature_window_k = [49.999, 60.0];
    assert!(matches!(NsfGrid3::<D3q27>::new(s, free, |_| (1.0, [0.0; 3], 50.0)), Err(NsfError::Gas(_))));
    // cell cap: exactly the cap admitted, one more refused
    let mut s = base();
    s.maximum_cells = 64;
    assert!(NsfGrid3::<D3q27>::new(s, free, rest).is_ok());
    let mut s = base();
    s.dims = [5, 4, 4];
    s.maximum_cells = 64;
    assert!(matches!(NsfGrid3::<D3q27>::new(s, free, rest), Err(NsfError::Cap { requested: 80, maximum: 64 })));
    // Mach bound (local sound speed sqrt(gamma theta)): at the bound admitted, bound x (1 + 1e-9) refused
    let mach = base().maximum_mach;
    let c = (1.4f64 * THETA0).sqrt();
    let at = |scale: f64| move |_x: [usize; 3]| (1.0, [mach * c * scale, 0.0, 0.0], 293.15);
    assert!(NsfGrid3::<D3q27>::new(base(), free, at(1.0 - 1e-12)).is_ok());
    assert!(matches!(NsfGrid3::<D3q27>::new(base(), free, at(1.0 + 1e-9)), Err(NsfError::Mach { .. })));
    // temperature ratio of the corrected equilibrium and the declared window
    let ratio = base().maximum_temperature_ratio;
    let t_hi = 293.15 * (1.0 + ratio * (1.0 + 1e-9));
    let mut s = base();
    s.temperature_window_k = [273.15, 400.0];
    assert!(matches!(NsfGrid3::<D3q27>::new(s, free, move |_| (1.0, [0.0; 3], t_hi)), Err(NsfError::TemperatureRatio { .. })));
    assert!(matches!(NsfGrid3::<D3q27>::new(base(), free, |_| (1.0, [0.0; 3], 400.0)), Err(NsfError::Window { .. })));
    // non-finite input
    assert!(matches!(NsfGrid3::<D3q27>::new(base(), free, |_| (f64::NAN, [0.0; 3], 293.15)), Err(NsfError::NonFinite { .. })));
    assert!(matches!(NsfGrid3::<D3q27>::new(base(), free, |_| (1.0, [f64::INFINITY, 0.0, 0.0], 293.15)), Err(NsfError::NonFinite { .. })));
}

#[test]
fn closed_box_energy_budget_closes() {
    // G0. Adiabatic rigid box (all six faces halfway bounce-back walls) from a random velocity + temperature state:
    // total energy (internal + kinetic) conserved to <= 1e-11 relative over 1e4 steps (flux telescoping), and the
    // per-step ledger residual closes at the same band.
    let mut rng = Lcg(7);
    let n = 10usize;
    let noise: Vec<(f64, [f64; 3], f64)> = (0..n * n * n)
        .map(|_| (1.0 + 0.02 * rng.next(), [0.03 * rng.next(), 0.03 * rng.next(), 0.03 * rng.next()], 293.15 * (1.0 + 0.02 * rng.next())))
        .collect();
    let s = spec([n, n, n], [false; 3], WallThermal::Adiabatic, 5.0e-5);
    let mut g = NsfGrid3::<D3q27>::new(s, |_| false, |x| noise[x[0] + n * (x[1] + n * x[2])]).unwrap();
    let e0 = g.total_energy_j();
    let mut worst_drift = 0.0f64;
    let mut worst_resid = 0.0f64;
    for step in 0..10_000 {
        let l = g.step().unwrap();
        worst_resid = worst_resid.max(l.residual_rel);
        if step % 100 == 99 {
            worst_drift = worst_drift.max(((g.total_energy_j() - e0) / e0).abs());
        }
    }
    let ledger = g.cumulative_ledger();
    println!("closed box: max drift {worst_drift:.3e}, max step residual {worst_resid:.3e}, kinetic {:.3e} J, viscous heating {:.3e} J", ledger.kinetic_j, ledger.viscous_heating_j);
    assert!(worst_drift <= 1e-11, "total energy drift {worst_drift}");
    assert!(worst_resid <= 1e-11, "ledger residual {worst_resid}");
    assert!(ledger.wall_heat_j.abs() <= 1e-11 * e0, "adiabatic walls must exchange no heat");
}

fn taylor_green(n: usize, dx: f64, workers: usize) -> NsfGrid3<D3q27> {
    let k = 2.0 * core::f64::consts::PI / n as f64;
    let amp = 0.02;
    let mut s = spec([n, n, 4], [true; 3], WallThermal::Adiabatic, dx);
    s.workers = workers;
    NsfGrid3::<D3q27>::new(s, |_| false, move |x| {
        let (xx, yy) = (x[0] as f64 + 0.5, x[1] as f64 + 0.5);
        (1.0, [-amp * (k * xx).cos() * (k * yy).sin(), amp * (k * xx).sin() * (k * yy).cos(), 0.0], 293.15)
    })
    .unwrap()
}

#[test]
fn viscous_heating_equals_kinetic_loss() {
    // G0 term identity: per step, viscous heating = - (kinetic change) - pressure work, to round-off of the TOTAL energy (it follows
    // from total-energy conservation and is bookkeeping; the round-off lives in the large internal energy, not in the small kinetic one).
    // G1 physical content: on a 2-D Taylor-Green vortex, after an initial acoustic adjustment (the start has uniform pressure), the heating
    // accumulated per step equals the continuum dissipation 4 nu k^2 KE, with the observed deviation shrinking as the grid refines.
    let mut deviations = Vec::new();
    let mut rate_deviations = Vec::new();
    for n in [16usize, 32] {
        let mut g = taylor_green(n, 2.0e-5, 4);
        let k = 2.0 * core::f64::consts::PI / n as f64;
        let nu = g.law().mu_lattice(THETA0); // rho = 1
        let mut kin = g.kinetic_energy_j();
        let steps = (0.3 / (4.0 * nu * k * k)) as usize;
        let settle = 60usize;
        let (mut heat, mut analytic) = (0.0, 0.0);
        let mut ke_log: Vec<f64> = Vec::new();
        let mut worst_identity = 0.0f64;
        for s in 0..steps {
            let l = g.step().unwrap();
            let d_kin = l.kinetic_j - kin;
            let id = (l.viscous_heating_j - (-d_kin - l.pressure_work_j)).abs() / g.total_energy_j();
            worst_identity = worst_identity.max(id);
            if s >= settle {
                heat += l.viscous_heating_j;
                analytic += 4.0 * nu * k * k * kin;
                ke_log.push(l.kinetic_j.ln());
            }
            kin = l.kinetic_j;
        }
        let dev = (heat / analytic - 1.0).abs();
        // decay rate of the kinetic energy by a least-squares fit of ln KE after the settle
        let fit = {
            let m = ke_log.len() as f64;
            let xs: Vec<f64> = (0..ke_log.len()).map(|i| i as f64).collect();
            let (sx, sy) = (xs.iter().sum::<f64>(), ke_log.iter().sum::<f64>());
            let sxy: f64 = xs.iter().zip(&ke_log).map(|(x, y)| x * y).sum();
            let sxx: f64 = xs.iter().map(|x| x * x).sum();
            -(m * sxy - sx * sy) / (m * sxx - sx * sx)
        };
        let rate_ratio = fit / (4.0 * nu * k * k);
        println!("TGV n={n}: kinetic-energy decay rate / (4 nu k^2) = {rate_ratio:.5}");
        rate_deviations.push((rate_ratio - 1.0).abs());
        println!("TGV n={n}: steps {steps}, identity residual / total energy {worst_identity:.2e}, heating/analytic - 1 = {:+.3e}", heat / analytic - 1.0);
        assert!(worst_identity <= 1e-13, "term identity {worst_identity}");
        deviations.push(dev);
    }
    assert!(deviations[1] <= 2e-2, "heating vs continuum dissipation at n=32: {}", deviations[1]);
    assert!(deviations[1] < deviations[0], "deviation must shrink with refinement: {deviations:?}");
    assert!(rate_deviations[1] <= 5e-3, "kinetic-energy decay rate at n=32: {rate_deviations:?}");
}

#[test]
fn isothermal_wall_heat_ledger() {
    // G0. Isothermal walls 10 K above a gas at rest: the heat that crossed the walls (time integral of the wall flux) equals
    // the change of total energy to <= 1e-11 relative and it is positive; the semi-infinite-wall estimate Q = 2 kappa dT sqrt(t / (pi alpha))
    // per area is printed for orientation, not asserted (closed box: the heated layer compresses the bulk).
    let n = 24usize;
    let mut s = spec([n, n, n], [false; 3], WallThermal::Isothermal(303.15), 5.0e-6);
    s.temperature_window_k = [273.15, 333.15];
    let mut g = NsfGrid3::<D3q27>::new(s, |_| false, |_| (1.0, [0.0; 3], 293.15)).unwrap();
    let e0 = g.total_energy_j();
    let steps = 600usize;
    let mut wall = 0.0;
    for _ in 0..steps {
        wall += g.step().unwrap().wall_heat_j;
    }
    let d_total = g.total_energy_j() - e0;
    println!("isothermal wall: wall heat {wall:.6e} J, total-energy change {d_total:.6e} J");
    assert!(wall > 0.0);
    assert!(((wall - d_total) / d_total).abs() <= 1e-11, "wall heat {wall} vs energy change {d_total}");
    let law = g.law();
    let kappa = law.kappa_lattice(THETA0);
    let cp = law.gamma / (law.gamma - 1.0);
    let alpha = kappa / cp; // rho = 1
    let dtheta = law.theta_of_kelvin(303.15) - THETA0;
    let area = 6.0 * (n * n) as f64;
    let q_lattice = area * 2.0 * kappa * dtheta * (steps as f64 / (core::f64::consts::PI * alpha)).sqrt();
    let want = q_lattice * law.energy_unit_j();
    println!("semi-infinite-wall estimate {want:.6e} J, ratio {:.3}", wall / want);
    // informational only: the estimate ignores the closed-box piston effect and the boundary layer is about 2 cells thick here
}

fn sound_error(lambda: usize) -> f64 {
    let k = 2.0 * core::f64::consts::PI / lambda as f64;
    let eps = 1.0e-4;
    let s = spec([lambda, 4, 4], [true; 3], WallThermal::Adiabatic, 5.0e-5);
    let gamma = 1.4;
    let mut g = NsfGrid3::<D3q27>::new(s, |_| false, move |x| {
        let rho = 1.0 + eps * (k * x[0] as f64).sin();
        (rho, [0.0; 3], 293.15 * rho.powf(gamma - 1.0))
    })
    .unwrap();
    let c0 = (gamma * THETA0).sqrt();
    let steps = (10.0 * lambda as f64 / c0) as usize;
    let mut amp = Vec::with_capacity(steps);
    for _ in 0..steps {
        let mut a = 0.0;
        for x in 0..lambda {
            a += g.cell([x, 0, 0]).0 * (k * x as f64).sin();
        }
        amp.push(2.0 * a / lambda as f64);
        g.step().unwrap();
    }
    let mut crossings = Vec::new();
    for i in 0..amp.len() - 1 {
        if amp[i] * amp[i + 1] < 0.0 {
            crossings.push(i as f64 + amp[i] / (amp[i] - amp[i + 1]));
        }
    }
    let halves: Vec<f64> = crossings.windows(2).skip(2).map(|w| w[1] - w[0]).collect();
    let period = 2.0 * halves.iter().sum::<f64>() / halves.len() as f64;
    (2.0 * core::f64::consts::PI / period) / k / c0 - 1.0
}

#[test]
fn adiabatic_sound_speed() {
    // G1. Small-amplitude acoustic wave in a periodic box: the measured speed is sqrt(gamma R T), NOT the isothermal
    // sqrt(R T) (ratio sqrt(gamma)), with dispersion error that converges at second order. The bead's band, 1e-4, is
    // asserted at the finest rung of the 48 / 96 / 192 / 288 points-per-wavelength ladder; the 32-point value is printed
    // and not asserted (the athermal D3Q27 LBM alone has about -1.8e-4 there).
    let e32 = sound_error(32);
    println!("32 points per wavelength: relative error {e32:+.3e}");
    let rungs = [48usize, 96, 192, 288];
    let errs: Vec<f64> = rungs.iter().map(|&l| sound_error(l)).collect();
    for (l, e) in rungs.iter().zip(&errs) {
        println!("{l} points per wavelength: relative error {e:+.3e}");
    }
    let order = (errs[2].abs() / errs[3].abs()).log2() / (288.0f64 / 192.0).log2();
    println!("observed order between 192 and 288: {order:.2}");
    assert!(errs[3].abs() <= 1e-4, "finest rung {}", errs[3]);
    assert!(errs.windows(2).all(|w| w[0].abs() > w[1].abs()), "error must shrink: {errs:?}");
    assert!(order >= 1.8, "observed order {order}");
    // mutation: the isothermal sound speed would be a factor sqrt(gamma) lower; the measured speed must NOT be within 1e-4 of it
    let c_ratio = errs[3] + 1.0; // measured / sqrt(gamma theta0)
    let iso = 1.0 / 1.4f64.sqrt();
    assert!((c_ratio / iso - 1.0).abs() > 0.1, "measured speed is isothermal-like: {c_ratio}");
}

#[test]
fn nsf3_bitwise_worker_count_invariance() {
    // G5. The same Taylor-Green run with 1, 2, 3 and 7 workers gives identical state bits and identical ledger bits.
    let mut reference: Option<(u64, [u64; 3])> = None;
    for workers in [1usize, 2, 3, 7] {
        let mut g = taylor_green(16, 2.0e-5, workers);
        let mut last = Default::default();
        for _ in 0..40 {
            last = g.step().unwrap();
        }
        let got = (g.state_hash(), [last.viscous_heating_j.to_bits(), last.pressure_work_j.to_bits(), last.kinetic_j.to_bits()]);
        match &reference {
            None => reference = Some(got),
            Some(r) => assert_eq!(*r, got, "workers {workers}"),
        }
    }
}
