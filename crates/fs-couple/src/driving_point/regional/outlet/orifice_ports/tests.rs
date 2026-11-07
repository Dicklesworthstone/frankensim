use super::*;
use fs_duct::{Duct, HoleState, Segment};
use fs_material::gas::GasState;
use fs_material::gas_species::{GasSpecies, resolve_gas_mixture};
use fs_qty::semantic::{Composition, CompositionBasis};

const RATE: u32 = 48_000;
const WALLS: WallModel<'static> = WallModel::Rigid;

fn air() -> GasState {
    resolve_gas_mixture(
        &[GasSpecies::Nitrogen, GasSpecies::Oxygen],
        &Composition::new(CompositionBasis::MoleFraction, vec![0.7, 0.3]).unwrap(),
        293.15,
        100_000.0,
    )
    .unwrap()
}

fn tone_hole(radius: f64, chimney: f64, state: HoleState) -> Segment {
    Segment::ToneHole { hole_radius: radius, chimney_height: chimney, bore_radius: 0.00745, state }
}

fn cylinder(length: f64) -> Segment {
    Segment::Cylinder { radius: 0.00745, length }
}

/// Three open holes of different size (so the port references differ), one closed hole, a 7.45 mm bore.
fn three_hole_duct(ports: usize) -> RegionalDuct {
    let gas = air();
    let segments = vec![
        cylinder(0.12),
        tone_hole(0.00135, 0.0127, HoleState::Open),
        cylinder(0.02),
        tone_hole(0.00235, 0.0069, HoleState::Closed),
        cylinder(0.02),
        tone_hole(0.00235, 0.0069, HoleState::Open),
        cylinder(0.02),
        tone_hole(0.0020, 0.0050, HoleState::Open),
        cylinder(0.25),
    ];
    let count = segments.len();
    RegionalDuct::new(&Duct { segments }, &vec![gas; count])
        .unwrap()
        .with_open_hole_orifice_ports(ports)
        .unwrap()
}

fn one_hole_duct(ports: usize) -> RegionalDuct {
    let gas = air();
    let segments = vec![cylinder(0.15), tone_hole(0.00135, 0.0127, HoleState::Open), cylinder(0.11)];
    RegionalDuct::new(&Duct { segments }, &[gas; 3]).unwrap().with_open_hole_orifice_ports(ports).unwrap()
}

fn lcg(state: &mut u64) -> f64 {
    *state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
    ((*state >> 11) as f64) / ((1_u64 << 53) as f64)
}

#[test]
fn the_card_count_is_bounded() {
    let gas = air();
    let duct = RegionalDuct::new(&Duct { segments: vec![cylinder(0.2)] }, &[gas]).unwrap();
    assert_eq!(duct.open_hole_orifice_ports(), 1);
    assert!(duct.clone().with_open_hole_orifice_ports(0).is_err());
    assert!(duct.clone().with_open_hole_orifice_ports(MAX_ORIFICE_PORTS + 1).is_err());
    assert_eq!(duct.clone().with_open_hole_orifice_ports(3).unwrap().open_hole_orifice_ports(), 3);
    assert_eq!(
        duct.with_open_hole_orifice_ports(MAX_ORIFICE_PORTS).unwrap().open_hole_orifice_ports(),
        MAX_ORIFICE_PORTS,
    );
}

#[test]
fn the_coupled_orifice_solve_converges_to_a_small_residual() {
    let mut seed = 0x9E37_79B9_7F4A_7C15_u64;
    for trial in 0..1200 {
        // n runs over 1..=NP (the cap, 16); the coupling amplitude is 0.3 up to three ports, as before, and
        // 0.9 / n beyond, so every row sum of |D| stays under 0.45 and the system stays passive.
        let n = 1 + trial % NP;
        let amp = 0.3 * 3.0 / n.max(3) as f64;
        let mut d = [[0.0_f64; NP]; NP];
        let mut r = [1.0_f64; NP];
        let mut k = [0.0_f64; NP];
        let mut u_th = [0.0_f64; NP];
        for i in 0..n {
            r[i] = 1.0e5 * (1.0 + lcg(&mut seed));
            k[i] = 10.0_f64.powf(8.0 + 3.0 * lcg(&mut seed));
            u_th[i] = 2.0 * (lcg(&mut seed) - 0.5) * 10.0_f64.powf(-5.0 + 3.0 * lcg(&mut seed));
            for j in 0..n {
                d[i][j] = amp * (lcg(&mut seed) - 0.5);
            }
        }
        let mut i_plus_d = d;
        for i in 0..n {
            i_plus_d[i][i] += 1.0;
        }
        let inverse = invert_real(&i_plus_d, n).unwrap();
        let mut a = [[0.0_f64; NP]; NP];
        for i in 0..n {
            for j in 0..n {
                a[i][j] = (2.0 * inverse[i][j] - if i == j { 1.0 } else { 0.0 }) / r[i];
            }
        }
        let u = solve_orifices(&a, &k, &u_th, n).unwrap();
        let scale = u_th[..n].iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        for i in 0..n {
            let mut f = u[i] - u_th[i];
            for j in 0..n {
                f += a[i][j] * k[j] * u[j].abs() * u[j];
            }
            assert!(f.abs() <= 1.0e-12 * scale, "trial {trial} n {n} port {i}: residual {f:e} vs scale {scale:e}");
        }
    }
}

#[test]
fn zero_loss_leaves_the_linear_flow_untouched() {
    // The solver arrays are NP (the cap, 16) wide; the three active ports sit in the leading block.
    let block = [[1.0e-5, 2.0e-6, 0.0], [2.0e-6, 1.5e-5, 0.0], [0.0, 0.0, 1.0e-5]];
    let mut a = [[0.0_f64; NP]; NP];
    for i in 0..3 {
        a[i][..3].copy_from_slice(&block[i]);
    }
    let mut u_th = [0.0_f64; NP];
    u_th[..3].copy_from_slice(&[1.0e-4, -3.0e-4, 2.0e-4]);
    assert_eq!(solve_orifices(&a, &[0.0; NP], &u_th, 3).unwrap(), u_th);
}

#[test]
fn the_dispatch_keeps_the_first_hole_line_for_one_port_and_for_a_single_open_hole() {
    let walls = WALLS;
    for (duct, label) in [(three_hole_duct(1), "N=1"), (one_hole_duct(3), "one open hole, N=3")] {
        let single = duct.register_hole_shunt_line(RATE, 1024, walls).unwrap().unwrap();
        let chosen = duct.select_register_hole_shunt(Some(single.clone()), RATE, 1024, walls).unwrap().unwrap();
        let RegisterHoleShunt::Single(mut chosen) = chosen else { panic!("{label}: expected the first-hole line") };
        let mut reference = single;
        chosen.prime(0.0).unwrap();
        reference.prime(0.0).unwrap();
        for t in 0..400 {
            let p = 900.0 * (0.07 * f64::from(t)).sin();
            let a = chosen.step(p).unwrap();
            let b = reference.step(p).unwrap();
            assert_eq!(a.to_bits(), b.to_bits(), "{label}: flow differs at sample {t}");
            assert_eq!(chosen.s12_history().to_bits(), reference.s12_history().to_bits(), "{label}: s12 at {t}");
        }
    }
    let ports = three_hole_duct(3);
    let single = ports.register_hole_shunt_line(RATE, 1024, walls).unwrap();
    let chosen = ports.select_register_hole_shunt(single, RATE, 1024, walls).unwrap().unwrap();
    assert!(matches!(&chosen, RegisterHoleShunt::Ports(p) if p.ports() == 3));
    let two = three_hole_duct(2);
    let single = two.register_hole_shunt_line(RATE, 1024, walls).unwrap();
    assert!(matches!(
        two.select_register_hole_shunt(single, RATE, 1024, walls).unwrap().unwrap(),
        RegisterHoleShunt::Ports(p) if p.ports() == 2
    ));
}

#[test]
fn one_realized_port_follows_the_first_hole_line_through_a_strongly_nonlinear_drive() {
    let duct = three_hole_duct(1);
    let walls = WALLS;
    let mut line = duct.register_hole_shunt_line(RATE, 1024, walls).unwrap().unwrap();
    // The first-hole line uses the triangle rule, so the equivalence is exact only under that rule.
    let kernels = duct.realize_orifice_ports_with(RATE, 1024, walls, 1, PassivityRule::Triangle).unwrap();
    let mut ports = RegisterHoleOrificePorts::from_realized(&kernels, 1.0 / f64::from(RATE)).unwrap();
    line.prime(0.0).unwrap();
    ports.prime(0.0).unwrap();
    let mut worst_flow = 0.0_f64;
    let mut worst_history = 0.0_f64;
    let mut peak_flow = 0.0_f64;
    let mut peak_history = 0.0_f64;
    for t in 0..2000 {
        let p = 1500.0 * (0.0654 * f64::from(t)).sin() + 400.0 * (0.211 * f64::from(t)).cos();
        let u_line = line.step(p).unwrap();
        let u_ports = ports.step(p).unwrap()[0];
        worst_flow = worst_flow.max((u_line - u_ports).abs());
        peak_flow = peak_flow.max(u_line.abs());
        let (h_line, h_ports) = (line.s12_history(), ports.s12_history());
        worst_history = worst_history.max((h_line - h_ports).abs());
        peak_history = peak_history.max(h_line.abs());
    }
    assert!(peak_flow > 1.0e-6 && peak_history > 1.0, "the drive must reach the orifice: {peak_flow:e} {peak_history:e}");
    assert!(worst_flow <= 1.0e-6 * peak_flow, "flow: worst {worst_flow:e} of peak {peak_flow:e}");
    assert!(worst_history <= 1.0e-6 * peak_history, "s12: worst {worst_history:e} of peak {peak_history:e}");
}

fn dft_bin(samples: &[f64], cycles: usize) -> C64 {
    // Phasor of the e^{-i w t} convention: (2/N) sum x[t] e^{+i w t}.
    let n = samples.len();
    let mut re = 0.0;
    let mut im = 0.0;
    for (t, &x) in samples.iter().enumerate() {
        let angle = core::f64::consts::TAU * (cycles * t) as f64 / n as f64;
        re += x * angle.cos();
        im += x * angle.sin();
    }
    C64::new(2.0 * re / n as f64, 2.0 * im / n as f64)
}

/// Relative error of the time-domain harmonic steady state for each port against (1) the TMM network,
/// `(I + Y rho)^-1 u_0`, and (2) the same linear network built from the REALIZED kernels' own frequency
/// responses (which isolates the step algebra from the realization's rolloff, truncation and damping);
/// the third value is the tail factor the passivity treatment applied.
fn steady_state_errors(n_grid: usize, rho_scale: [f64; 3], cycles: usize) -> ([f64; 3], [f64; 3], f64) {
    let duct = three_hole_duct(3);
    let kernels = duct.realize_orifice_ports(RATE, n_grid, WALLS, 3).unwrap();
    let rho: Vec<f64> = kernels.r_ref.iter().zip(rho_scale).map(|(r, f)| r * f).collect();
    let window = 4800;
    let f_hz = f64::from(RATE) * cycles as f64 / window as f64;
    let mut ports = RegisterHoleOrificePorts::from_realized(&kernels, 1.0 / f64::from(RATE)).unwrap();
    ports.prime(0.0).unwrap();
    let total = 3 * window;
    let mut flows: Vec<Vec<f64>> = vec![Vec::new(); 3];
    for t in 0..total {
        let p = (core::f64::consts::TAU * f_hz * t as f64 / f64::from(RATE)).cos();
        let u = ports.step_linear(p, &rho).unwrap();
        if t >= total - window {
            for i in 0..3 {
                flows[i].push(u[i]);
            }
        }
    }
    let omega = core::f64::consts::TAU * f_hz;
    let (response, hole_flows, y) = duct
        .response_with_wall_holes_and_ports(omega, LossModel::Bessel, Termination::FlangedOpen, WALLS, 3)
        .unwrap();
    let zc = duct.inlet_impedance();
    let two_over = C64::new(response.impedance.re + zc, response.impedance.im).recip().scale(2.0);
    let u0: Vec<C64> = (0..3).map(|i| hole_flows[i].u_chimney_in_over_u_in * two_over).collect();
    // (I + Y rho) U = u0, solved by Cramer on the 3 x 3 system.
    let mut m = [[C64::ZERO; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = y[i * 3 + j].scale(rho[j]) + if i == j { C64::ONE } else { C64::ZERO };
        }
    }
    let det3 = |a: &[[C64; 3]; 3]| {
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    };
    let d0 = det3(&m);
    // The realized network: S_r(w) = sum_t ir[t] e^{+i w t} and u0_r(w) from the kernels themselves, then
    // [(I + S) - rho R^-1 (S - I)] a = rho u0_r and U = u0_r + R^-1 (S - I) a.
    let theta = omega / f64::from(RATE);
    let transfer = |ir: &[f64]| {
        let (mut re, mut im) = (0.0, 0.0);
        for (t, &x) in ir.iter().enumerate() {
            re += x * (theta * t as f64).cos();
            im += x * (theta * t as f64).sin();
        }
        C64::new(re, im)
    };
    let mut s_real = [[C64::ZERO; NP]; NP];
    let mut system = [[C64::ZERO; NP]; NP];
    for i in 0..3 {
        for j in 0..3 {
            s_real[i][j] = transfer(&kernels.s22[i * 3 + j]);
        }
    }
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { C64::ONE } else { C64::ZERO };
            let s_minus_i = s_real[i][j] - delta;
            system[i][j] = (delta + s_real[i][j]) - s_minus_i.scale(rho[i] / kernels.r_ref[i]);
        }
    }
    let u0_real: Vec<C64> = (0..3).map(|i| transfer(&kernels.s21[i])).collect();
    let inverse = invert_complex(&system, 3).unwrap();
    let mut wave = [C64::ZERO; 3];
    for i in 0..3 {
        for j in 0..3 {
            wave[i] = wave[i] + inverse[i][j] * u0_real[j].scale(rho[j]);
        }
    }
    let mut errors = [0.0_f64; 3];
    let mut errors_real = [0.0_f64; 3];
    for i in 0..3 {
        let mut mi = m;
        for r in 0..3 {
            mi[r][i] = u0[r];
        }
        let predicted = det3(&mi) * d0.recip();
        let measured = dft_bin(&flows[i], cycles);
        errors[i] = (predicted - measured).abs() / predicted.abs();
        let mut b_minus_a = C64::ZERO;
        for j in 0..3 {
            let delta = if i == j { C64::ONE } else { C64::ZERO };
            b_minus_a = b_minus_a + (s_real[i][j] - delta) * wave[j];
        }
        let predicted_real = u0_real[i] + b_minus_a.scale(1.0 / kernels.r_ref[i]);
        errors_real[i] = (predicted_real - measured).abs() / predicted_real.abs();
    }
    (errors, errors_real, kernels.tail_damp)
}

#[test]
fn the_time_domain_ports_reproduce_the_frequency_domain_network_with_linear_drops() {
    // Replace the three orifices by linear resistors Delta p_i = rho_i U_i. The harmonic steady state of the
    // time-domain multiport must then be (I + Y diag(rho))^-1 u_0, with Y and u_0 straight from the TMM chain.
    // The only departures are the realization's own (rolloff, the n/2-tap truncation, the passivity damping),
    // which shrink as the grid grows; rho = 0 isolates the incident-flow kernels from the coupling.
    let mut table = String::new();
    let mut worst_real = 0.0_f64;
    let mut worst_rho0 = 0.0_f64;
    let mut worst_coupled = 0.0_f64;
    for (label, scale) in [("rho 0", [0.0, 0.0, 0.0]), ("rho set", [0.3, 0.6, 1.0])] {
        for cycles in [10_usize, 23] {
            let (e, e_real, damp) = steady_state_errors(4096, scale, cycles);
            table += &format!(
                "{label} cycles {cycles}: vs TMM {:.2e} {:.2e} {:.2e}; vs realized kernels {:.2e} {:.2e} {:.2e}; tail factor {damp:.4}\n",
                e[0], e[1], e[2], e_real[0], e_real[1], e_real[2],
            );
            worst_real = e_real.iter().fold(worst_real, |m, v| m.max(*v));
            if label == "rho 0" {
                worst_rho0 = e.iter().fold(worst_rho0, |m, v| m.max(*v));
            } else {
                worst_coupled = e.iter().fold(worst_coupled, |m, v| m.max(*v));
            }
        }
    }
    for ports in 1..=3 {
        let exact = three_hole_duct(3).realize_orifice_ports(RATE, 4096, WALLS, ports).unwrap();
        let triangle = three_hole_duct(3)
            .realize_orifice_ports_with(RATE, 4096, WALLS, ports, PassivityRule::Triangle)
            .unwrap();
        table += &format!(
            "{ports} port(s): tail factor exact {:.4}, triangle {:.4}\n", exact.tail_damp, triangle.tail_damp,
        );
        // The exact rule keeps every factor the triangle rule keeps, and more (to the 2^-14 bisection step).
        assert!(exact.tail_damp >= triangle.tail_damp - 2.0e-4, "{table}");
    }
    eprintln!("{table}");
    // The step algebra against the realized kernels' own network: only the finite window limits it.
    assert!(worst_real <= 1.5e-2, "step algebra departs from the realized network: {worst_real:e}\n{table}");
    // The incident-flow kernels against the TMM chain: the n/2-tap truncation of a 4096 grid (0.3 % at 16384).
    assert!(worst_rho0 <= 3.0e-2, "incident-flow kernels depart from the chain: {worst_rho0:e}\n{table}");
    // The coupled case also carries the realization's rolloff, truncation and passivity damping. The triangle
    // rule's 0.9675 tail factor cost up to 15 % here (100 Hz); the exact rule leaves the truncation alone.
    assert!(worst_coupled <= 0.06, "coupled network departs from the chain: {worst_coupled:e}\n{table}");
}

#[test]
fn the_realized_matrix_is_passive_in_the_power_wave_norm() {
    let duct = three_hole_duct(3);
    let n = 1024;
    let kernels = duct.realize_orifice_ports(RATE, n, WALLS, 3).unwrap();
    let m = kernels.ports;
    let fft = Fft::new(2 * n);
    let mut scratch = vec![FftC64::new(0.0, 0.0); 2 * n];
    let mut spectra: Vec<Vec<FftC64>> = Vec::new();
    for ir in &kernels.s22 {
        let mut buf = vec![FftC64::new(0.0, 0.0); 2 * n];
        for (dst, &src) in buf.iter_mut().zip(ir.iter()) {
            *dst = FftC64::new(src, 0.0);
        }
        fft.forward(&mut buf, &mut scratch);
        spectra.push(buf);
    }
    let mut worst = 0.0_f64;
    for bin in 0..2 * n {
        let mut t = [[C64::ZERO; NP]; NP];
        for i in 0..m {
            for j in 0..m {
                let c = spectra[i * m + j][bin];
                t[i][j] = C64::new(c.re, c.im);
            }
        }
        worst = worst.max(weighted_spectral_norm(&t, &kernels.r_ref, m));
    }
    assert!(worst < 1.0, "the realized S22 matrix must be strictly passive, worst singular value {worst}");
}

#[test]
fn the_realized_matrix_is_the_admittance_scattering_not_its_transpose() {
    // The hole references differ, so S = (I + R Y)^-1 (I - R Y) and its transpose are different matrices.
    let duct = three_hole_duct(3);
    let n = 4096;
    let kernels = duct.realize_orifice_ports(RATE, n, WALLS, 3).unwrap();
    let m = 3;
    let mut err_s = 0.0;
    let mut err_t = 0.0;
    let mut size = 0.0;
    for k in (8..n / 16).step_by(7) {
        let theta = core::f64::consts::TAU * k as f64 / n as f64;
        let (_, _, y) = duct
            .response_with_wall_holes_and_ports(
                theta * f64::from(RATE), LossModel::Bessel, Termination::FlangedOpen, WALLS, m,
            )
            .unwrap();
        let s = scattering_from_admittance(&y, &kernels.r_ref, m).unwrap();
        for i in 0..m {
            for j in 0..m {
                // Acoustic e^{-iwt}: the transfer is sum ir[t] e^{+i theta t}.
                let mut re = 0.0;
                let mut im = 0.0;
                for (t, &x) in kernels.s22[i * m + j].iter().enumerate().take(n / 2 + 1) {
                    re += x * (theta * t as f64).cos();
                    im += x * (theta * t as f64).sin();
                }
                let realized = C64::new(re, im);
                err_s += (realized - s[i][j]).abs().powi(2);
                err_t += (realized - s[j][i]).abs().powi(2);
                size += s[i][j].abs().powi(2);
            }
        }
    }
    assert!(err_s <= 0.02 * size, "realized S differs from the admittance law: {err_s:e} of {size:e}");
    assert!(err_s < 0.25 * err_t, "the matrix looks transposed: error {err_s:e} vs transposed {err_t:e}");
}

#[test]
fn a_strongly_driven_stiff_three_port_run_stays_finite_and_never_stalls_the_solver() {
    let duct = three_hole_duct(3);
    let kernels = duct.realize_orifice_ports(RATE, 1024, WALLS, 3).unwrap();
    for stiffness in [1.0, 1.0e3] {
        let mut ports = RegisterHoleOrificePorts::from_realized(&kernels, 1.0 / f64::from(RATE)).unwrap();
        for k in &mut ports.k_quad {
            *k *= stiffness;
        }
        ports.prime(0.0).unwrap();
        let mut dissipated = 0.0_f64;
        for t in 0..4000 {
            let phase = 0.0654 * f64::from(t);
            let p = 4000.0 * phase.sin() * (1.0 + 0.5 * (0.0031 * f64::from(t)).sin()) + 700.0 * (3.1 * phase).cos();
            let u = ports.step(p).unwrap_or_else(|e| panic!("stiffness {stiffness}, sample {t}: {e:?}"));
            for i in 0..3 {
                assert!(u[i].is_finite());
                dissipated += ports.k_quad[i] * u[i].abs().powi(3);
            }
            assert!(ports.s12_history().is_finite());
        }
        assert!(dissipated > 0.0, "the orifices must absorb something at stiffness {stiffness}");
    }
}

#[test]
fn a_multi_region_duct_refuses_the_coupled_ports() {
    let air = air();
    let warm = resolve_gas_mixture(
        &[GasSpecies::Nitrogen, GasSpecies::Oxygen],
        &Composition::new(CompositionBasis::MoleFraction, vec![0.7, 0.3]).unwrap(),
        303.15,
        air.pressure,
    )
    .unwrap();
    let segments = vec![
        cylinder(0.12),
        tone_hole(0.00135, 0.0127, HoleState::Open),
        cylinder(0.02),
        tone_hole(0.00235, 0.0069, HoleState::Open),
        cylinder(0.02),
        cylinder(0.25),
    ];
    let gases = [air, air, air, air, air, warm];
    let duct = RegionalDuct::new(&Duct { segments }, &gases).unwrap().with_open_hole_orifice_ports(2).unwrap();
    let single = duct.register_hole_shunt_line(RATE, 1024, WALLS).unwrap();
    assert!(duct.select_register_hole_shunt(single, RATE, 1024, WALLS).is_err());
}
