//! A free shell's radiation loss on its own rest modes (`Body::attach_shell_radiation`, MB 14574 (d)).
use super::*;
use crate::shell_radiation::{RayleighGrid, rest_modes};
use std::f64::consts::PI;

fn flat(name: &str, radius: f64, grid: f64, thickness: f64, young: f64, density: f64) -> Body {
    let th = move |_r: f64| thickness;
    Body::plate_without_surface_loss(name, &PlateSpec {
        radius, hole: 0.0, grid, edge: Edge::Free, thickness: &th, rest_height: &|_| 0.0,
        young, poisson: 0.33, density, tension: 0.0, bending: true, nonlinear: Nonlinear::None,
        t60: [(100.0, 1.0), (1000.0, 1.0)], support: None, muffle: None,
    }).unwrap()
}

/// |M^-1/2 (K phi - omega^2 M phi)| for a mode with phi' M phi = 1, K = `rest_stiffness_apply`.
fn rest_residual(body: &mut Body, omega: f64, phi: &[f64]) -> f64 {
    let k_phi = body.rest_stiffness_apply(phi).unwrap();
    k_phi.iter().zip(phi).zip(&body.mass)
        .map(|((k, p), m)| (k - omega * omega * m * p).powi(2) / m)
        .sum::<f64>()
        .sqrt()
}

/// A free disc's rest modes are exact for this body's own operator: each one's residual against
/// `rest_stiffness_apply` is the dense solver's rounding, which scales with the grid's largest
/// eigenvalue omega_max^2, not with the mode's own. The two (2,0) modes, the (0,1) mode and the (3,0)
/// pair are the only ones below lambda^2 = omega a^2 sqrt(rho h / D) = 14, in their parity sectors, and
/// M-orthonormal. How near they lie to Kirchhoff's roots is
/// `a_finer_free_disc_moves_every_rest_mode_toward_the_kirchhoff_root`'s question.
#[test]
fn a_flat_free_disc_has_exact_rest_modes_in_their_parity_sectors() {
    let (r, th, e, nu, rho) = (0.2, 1.0e-3, 1.0e11, 0.33, 8600.0);
    let mut body = flat("disc", r, 8.0e-3, th, e, rho);
    let scale = (rho * th * 12.0 * (1.0 - nu * nu) / (e * th.powi(3))).sqrt() * r * r;
    let rest = rest_modes(&mut body, 14.0 / scale).unwrap();
    assert_eq!(rest.rigid, 3, "an unsupported disc's translation and two tilts");
    assert!(rest.asymmetry < 1e-12 && rest.coupling < 1e-12, "asymmetry {:e}, coupling {:e}", rest.asymmetry, rest.coupling);
    assert_eq!(rest.sectors.iter().sum::<usize>(), body.n);
    let lambda2: Vec<f64> = rest.modes.iter().map(|m| m.0 * scale).collect();
    assert_eq!(lambda2.len(), 5, "{lambda2:?}");
    let top = rest.omega_max * rest.omega_max;
    for (i, (omega, phi, _)) in rest.modes.iter().enumerate() {
        let residual = rest_residual(&mut body, *omega, phi);
        assert!(residual <= 1e-10 * top, "mode {i}: residual {residual:e} against omega_max^2 {top:e}");
    }
    // (2,0): cos 2 theta even in x and y, sin 2 theta odd in both; (0,1) even; (3,0): odd in x or in y
    let pair = |a: usize, b: usize| {
        let mut s = [rest.modes[a].2, rest.modes[b].2];
        s.sort();
        s
    };
    assert_eq!(pair(0, 1), [0, 3]);
    assert_eq!(rest.modes[2].2, 0);
    assert_eq!(pair(3, 4), [1, 2]);
    for (i, (_, phi, _)) in rest.modes.iter().enumerate() {
        let norm: f64 = phi.iter().zip(&body.mass).map(|(p, m)| p * p * m).sum();
        assert!((norm - 1.0).abs() < 1e-12, "mode {i}: phi' M phi = {norm}");
        for (_, other, _) in &rest.modes[..i] {
            let gram: f64 = phi.iter().zip(&body.mass).zip(other).map(|((p, m), q)| p * m * q).sum();
            assert!(gram.abs() < 1e-12, "mode {i}: Gram term {gram:e}");
        }
    }
}

/// Kirchhoff's free circular plate, nu = 0.33: the roots lambda^2 of the frequency equation (Leissa
/// 1969, NASA SP-160, eq. 2.21), computed with scipy (the same script returns Leissa's printed nu = 0.3
/// roots 5.358, 9.003 and 12.439): (2,0) 5.2620, (0,1) 9.0689 and (3,0) 12.2439. On the 8 mm grid (25
/// cells per radius) the (2,0) pair lies 4.57 % and 2.47 % below them, (0,1) 1.03 % and the (3,0) pair
/// 6.67 % (c41b, 11-08): that grid's own staircase free edge, the error growing with the nodal
/// diameters. Halving the step to 4 mm must move every mode toward its root, and the (2,0) and (3,0)
/// deviations to at most 0.6 of the 8 mm ones (a staircase edge is about first order in the step, so
/// they about halve). The factors were set before the 4 mm disc was run. Modes are matched by parity
/// sector and rank within it. Ignored in the lib suite because the 4 mm disc's four dense sectors hold
/// about 1960 vectors each; run by hand with `cargo test --release --lib -- --ignored a_finer_free_disc`.
#[test]
#[ignore = "the 4 mm disc's four dense eigensolves of about 1960 vectors each take minutes in release"]
fn a_finer_free_disc_moves_every_rest_mode_toward_the_kirchhoff_root() {
    let (r, th, e, nu, rho): (f64, f64, f64, f64, f64) = (0.2, 1.0e-3, 1.0e11, 0.33, 8600.0);
    let scale = (rho * th * 12.0 * (1.0 - nu * nu) / (e * th.powi(3))).sqrt() * r * r;
    // each parity sector's modes below lambda^2 = 14, ascending: (2,0) cos and (0,1); (3,0); (3,0); (2,0) sin
    let roots: [&[f64]; 4] = [&[5.2620, 9.0689], &[12.2439], &[12.2439], &[5.2620]];
    let deviations = |grid: f64| -> [Vec<f64>; 4] {
        let mut body = flat("disc", r, grid, th, e, rho);
        let rest = rest_modes(&mut body, 14.0 / scale).unwrap();
        assert_eq!(rest.rigid, 3, "{grid} m grid");
        let top = rest.omega_max * rest.omega_max;
        let mut out: [Vec<f64>; 4] = Default::default();
        for (i, (omega, phi, s)) in rest.modes.iter().enumerate() {
            let residual = rest_residual(&mut body, *omega, phi);
            assert!(residual <= 1e-10 * top, "{grid} m grid, mode {i}: residual {residual:e} against omega_max^2 {top:e}");
            out[*s].push(omega * scale);
        }
        for (s, lambda2) in out.iter_mut().enumerate() {
            assert_eq!(lambda2.len(), roots[s].len(), "{grid} m grid, parity sector {s}: lambda^2 {lambda2:?}");
            for (l, root) in lambda2.iter_mut().zip(roots[s]) {
                *l = *l / root - 1.0;
            }
        }
        println!("{grid} m grid ({} nodes, sectors {:?}): deviation from Kirchhoff by parity sector {out:?}", body.n, rest.sectors);
        out
    };
    let (coarse, fine) = (deviations(8.0e-3), deviations(4.0e-3));
    for (s, (cs, fs)) in coarse.iter().zip(&fine).enumerate() {
        for (k, (c, f)) in cs.iter().zip(fs).enumerate() {
            let factor = if s == 0 && k == 1 { 1.0 } else { 0.6 };
            assert!(f.abs() < factor * c.abs(), "parity sector {s}, mode {k}: 8 mm {c:+.5}, 4 mm {f:+.5}, factor {factor}");
        }
    }
}

/// The FFT convolution is the double sum, and a uniform disc of nodes is the baffled circular
/// piston, R = rho c pi a^2 (1 - J1(2ka)/(ka)), a the equal-area radius; the bracket from scipy, and
/// the same double sum done directly in numpy lies within 2e-4 of it on this disc up to ka = 4.
#[test]
fn the_rayleigh_convolution_is_the_double_sum_and_a_uniform_disc_is_the_baffled_piston() {
    let (a, h, rho, c): (f64, f64, f64, f64) = (0.15, 0.01, 1.2, 343.0);
    let m = (a / h).ceil() as i64 + 1;
    let mut xy = vec![];
    for j in -m..=m {
        for i in -m..=m {
            if ((i * i + j * j) as f64).sqrt() * h <= a + 1e-9 {
                xy.push([i as f64 * h, j as f64 * h]);
            }
        }
    }
    assert_eq!(xy.len(), 709);
    let grid = RayleighGrid::new(&xy, h).unwrap();
    let aeq = (xy.len() as f64 * h * h / PI).sqrt();
    let uniform = vec![1.0; xy.len()];
    for (ka, piston) in [(0.1, 0.004_991_673_6), (0.5, 0.119_898_828_5), (1.0, 0.423_275_192_2), (2.0, 1.033_021_664_0), (4.0, 0.941_340_913_3)] {
        let k = ka / aeq;
        let u = grid.convolve(&uniform, k);
        let r = k * c * rho / (2.0 * PI) * h.powi(4) * u.iter().sum::<f64>();
        let want = rho * c * PI * aeq * aeq * piston;
        assert!((r / want - 1.0).abs() < 1e-3, "ka {ka}: R {r:.6e} vs piston {want:.6e}");
    }
    let shape: Vec<f64> = xy.iter().map(|p| 1.0 + 3.0 * p[0] - 7.0 * p[0] * p[1] + 20.0 * p[1].powi(3)).collect();
    let k = 2.0 / aeq;
    let u = grid.convolve(&shape, k);
    let bound = k * shape.iter().map(|s| s.abs()).sum::<f64>();
    for (x, ux) in xy.iter().zip(&u) {
        let direct: f64 = xy.iter().zip(&shape).map(|(y, s)| {
            let d = ((x[0] - y[0]).powi(2) + (x[1] - y[1]).powi(2)).sqrt();
            s * if d > 0.0 { (k * d).sin() / d } else { k }
        }).sum();
        assert!((ux - direct).abs() <= 1e-12 * bound, "FFT {ux:e} vs direct {direct:e}");
    }
}

/// A light disc (300 kg/m^3, an assumption, so that radiation dominates within milliseconds) started
/// in its most damped loaded mode: the mode's coordinate falls at the attached sigma = R_rad / 2, the
/// energy never rises and the shell's discrete energy plus the radiated work stays the start energy.
#[test]
fn a_rest_mode_decays_at_its_radiation_rate_and_the_energy_closes() {
    let (r, h) = (0.15, 0.01);
    let mut body = flat("light disc", r, h, 1.0e-3, 7.0e10, 300.0);
    let omega_max = rest_modes(&mut body.clone(), 0.0).unwrap().omega_max;
    let lines = body.attach_shell_radiation(&standard_gas(), r, 4000.0, 8000.0).unwrap();
    assert!(lines[0].starts_with("shell radiation (MB 14574 (d), a correction)") && lines[0].contains("upper bounds"), "{lines:?}");
    assert!(lines[1].starts_with("shell radiation modes") && lines[2].starts_with("shell radiation by octave"), "{lines:?}");
    let damping = body.shape_damping.as_ref().unwrap();
    assert_eq!(lines.len(), 3 + damping.two_sigma.len(), "one line per loaded mode");
    let largest = damping.two_sigma.iter().cloned().fold(0.0_f64, f64::max);
    assert!(damping.rate_bound >= largest * (1.0 - 1e-12) && damping.rate_bound <= largest * (1.0 + 1e-6),
        "rate bound {} against 2 max sigma {largest}: exact modes leave the Gram matrix the identity", damping.rate_bound);
    let pick = (0..damping.two_sigma.len()).max_by(|&i, &j| damping.two_sigma[i].total_cmp(&damping.two_sigma[j])).unwrap();
    let sigma = damping.two_sigma[pick] / 2.0;
    assert!(6.907755 / sigma < 0.5, "the test needs a radiation T60 under 0.5 s, got {} s", 6.907755 / sigma);
    let phi: Vec<f64> = damping.mphi[pick].iter().zip(&body.mass).map(|(p, m)| p / m).collect();
    let omega = dot(&phi, &body.rest_stiffness_apply(&phi).unwrap()).sqrt();
    let dt = 1.0 / omega_max;
    for i in 0..body.n {
        body.w[i] = 1e-6 * phi[i];
        body.wp[i] = body.w[i];
    }
    body.sync();
    let initial = body.discrete_energy(dt);
    assert!(initial > 0.0);
    let coordinate_of = |b: &Body| -> f64 { b.w.iter().zip(&b.mass).zip(&phi).map(|((w, m), p)| w * m * p).sum() };
    let mut coordinate = vec![coordinate_of(&body)];
    let steps = (1.5 / (sigma * dt)).ceil() as usize;
    assert!(steps < 400_000, "{steps} steps");
    let mut previous = initial;
    for _ in 0..steps {
        body.linear_forces();
        body.predict(dt);
        body.finish(dt);
        body.rotate();
        coordinate.push(coordinate_of(&body));
        let energy = body.discrete_energy(dt);
        assert!(energy >= 0.0 && energy <= previous + 2e-12 * initial);
        assert!((energy + body.shape_damping_dissipated_j() - initial).abs() < 2e-10 * initial,
            "E {energy:e} + D {:e} vs E0 {initial:e}", body.shape_damping_dissipated_j());
        previous = energy;
    }
    let quarter = (std::f64::consts::FRAC_PI_2 / omega / dt).round() as usize;
    let peaks: Vec<(f64, f64)> = (0..coordinate.len() - quarter)
        .filter(|&i| coordinate[i.saturating_sub(quarter)..=i + quarter].iter().all(|x| x.abs() <= coordinate[i].abs()))
        .map(|i| (i as f64 * dt, coordinate[i].abs().ln()))
        .collect();
    assert!(peaks.len() >= 10, "{} half-cycle peaks", peaks.len());
    let (mean_t, mean_y) = (
        peaks.iter().map(|p| p.0).sum::<f64>() / peaks.len() as f64,
        peaks.iter().map(|p| p.1).sum::<f64>() / peaks.len() as f64,
    );
    let slope = peaks.iter().map(|p| (p.0 - mean_t) * (p.1 - mean_y)).sum::<f64>()
        / peaks.iter().map(|p| (p.0 - mean_t).powi(2)).sum::<f64>();
    assert!((-slope / sigma - 1.0).abs() < 0.01, "the mode falls at {:.3} 1/s, attached {sigma:.3} 1/s", -slope);
}

/// The card crash on its felt, on a coarse 2 cm grid, through the kit: the curvature Hessian passes
/// the symmetry gate, the clock is re-admitted, a second attach refuses, and a struck render stays
/// finite while the load removes energy.
#[test]
fn the_card_crash_takes_the_load_through_the_kit_and_a_struck_render_stays_finite() {
    let c = kit::CymbalSpec { grid: 0.02, ..kit::card_crash() };
    let mut p = kit::cymbal_with_interface(&c, &[(0.0005, 2.0, 0.95, 0.008)], &kit::DECLARED_CYMBAL_INTERFACE).unwrap();
    let lines = kit::attach_shell_radiation(&mut p, &c).unwrap();
    assert!(lines[0].starts_with("shell radiation (MB 14574 (d)") && lines[1].starts_with("shell radiation limits")
        && lines[2].starts_with("shell radiation modes") && lines.len() >= 5, "{lines:?}");
    assert!(kit::attach_shell_radiation(&mut p, &c).is_err());
    let out = kit::render(&mut p, 0.01, |_| {});
    assert!(out[0].iter().all(|v| v.is_finite()));
    let body = p.instrument.bodies.iter().find(|b| b.name == c.name).unwrap();
    assert!(body.shape_damping_dissipated_j() > 0.0);
}

#[test]
fn shell_radiation_refuses_heads_bad_arguments_no_mode_and_a_second_attach() {
    let gas = standard_gas();
    let mut head = Body::plate_without_surface_loss("head", &PlateSpec {
        radius: 0.15, hole: 0.0, grid: 0.01, edge: Edge::Clamped, thickness: &|_| 2.54e-4, rest_height: &|_| 0.0,
        young: 4e9, poisson: 0.38, density: 1390.0, tension: 2000.0, bending: true, nonlinear: Nonlinear::None,
        t60: [(100.0, 1.0), (1000.0, 1.0)], support: None, muffle: None,
    }).unwrap();
    assert!(head.attach_shell_radiation(&gas, 0.15, 1000.0, 2000.0).is_err());
    let mut disc = flat("disc", 0.15, 0.015, 1.0e-3, 7.0e10, 2700.0);
    for (loaded, listed, a) in [(0.0, 100.0, 0.15), (-1.0, 100.0, 0.15), (200.0, 100.0, 0.15), (f64::NAN, 100.0, 0.15),
        (100.0, f64::INFINITY, 0.15), (100.0, 200.0, 0.0), (100.0, 200.0, f64::NAN)] {
        assert!(disc.attach_shell_radiation(&gas, a, loaded, listed).is_err(), "{loaded} {listed} {a}");
        assert!(disc.shape_damping.is_none());
    }
    assert!(disc.attach_shell_radiation(&gas, 0.15, 1.0, 100.0).is_err(), "no elastic mode at or below 1 Hz");
    assert!(disc.shape_damping.is_none());
    disc.attach_shell_radiation(&gas, 0.15, 2000.0, 4000.0).unwrap();
    assert!(disc.attach_shell_radiation(&gas, 0.15, 2000.0, 4000.0).is_err());
    assert!(Body::point("stick", 0.01).attach_shell_radiation(&gas, 0.15, 2000.0, 4000.0).is_err());
}

#[test]
#[should_panic(expected = "shell radiation")]
fn exported_stepper_refuses_shell_radiation() {
    let mut disc = flat("disc", 0.15, 0.015, 1.0e-3, 7.0e10, 2700.0);
    disc.attach_shell_radiation(&standard_gas(), 0.15, 2000.0, 4000.0).unwrap();
    let mut instrument = Instrument::new(96_000.0);
    instrument.add(disc);
    let _ = instrument.image();
}
