//! G1/G2/G3 battery for planar distributed-plasticity frames: closed-form
//! elastic stiffness, the corotational pure-bending roll-up, axial-load
//! frequency softening toward Euler buckling, RC pushover equilibrium and
//! capacity, cyclic dissipation, SDOF dynamics against the exact solution,
//! and a multi-story RC frame under base motion with an energy ledger.

use fs_solid::fiber::rc_section;
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel};

const E: f64 = 200e9;
const I: f64 = 8.0e-5;
const A: f64 = 1.0e-2;

fn elastic() -> SectionModel {
    SectionModel::Elastic {
        ea: E * A,
        ei: E * I,
    }
}

/// Vertical cantilever of `n` elements along +x from the origin, base fixed.
fn cantilever(geom: Geometry, l: f64, n: usize) -> (Frame2d, usize) {
    let mut f = Frame2d::new(geom);
    let mut prev = f.add_node(0.0, 0.0);
    f.fix(prev, [true, true, true]);
    for k in 1..=n {
        let node = f.add_node(l * k as f64 / n as f64, 0.0);
        f.add_element(prev, node, 5, 0.0, &elastic).unwrap();
        prev = node;
    }
    (f, prev)
}

#[test]
fn frame_g1_elastic_cantilever_matches_euler_bernoulli_exactly() {
    let (l, p) = (3.0, 1.0e4);
    let (mut f, tip) = cantilever(Geometry::Linear, l, 1);
    let mut load = vec![0.0; f.ndof()];
    load[3 * tip + 1] = p;
    load[3 * tip] = 5.0 * p;
    f.static_load(&load, 1).unwrap();
    let u = f.displacements();
    let v_exact = p * l.powi(3) / (3.0 * E * I);
    let th_exact = p * l * l / (2.0 * E * I);
    let ax_exact = 5.0 * p * l / (E * A);
    assert!(
        (u[3 * tip + 1] - v_exact).abs() < 1e-12 * v_exact.abs().max(1.0),
        "{}",
        u[3 * tip + 1]
    );
    assert!(
        (u[3 * tip + 2] - th_exact).abs() < 1e-12,
        "{}",
        u[3 * tip + 2]
    );
    assert!((u[3 * tip] - ax_exact).abs() < 1e-14, "{}", u[3 * tip]);
}

#[test]
fn frame_g1_portal_lateral_stiffness_matches_closed_form() {
    // Fixed-base portal: K = (24 E I_c / h³)(6ρ + 1)/(6ρ + 4), ρ = (I_b/L)/(I_c/h).
    let (h, span) = (4.0, 6.0);
    let (ic, ib) = (I, 2.0 * I);
    let big = 1e6; // axially rigid members
    let mut f = Frame2d::new(Geometry::Linear);
    let n0 = f.add_node(0.0, 0.0);
    let n1 = f.add_node(0.0, h);
    let n2 = f.add_node(span, h);
    let n3 = f.add_node(span, 0.0);
    f.fix(n0, [true; 3]);
    f.fix(n3, [true; 3]);
    let col = || SectionModel::Elastic {
        ea: E * A * big,
        ei: E * ic,
    };
    let beam = || SectionModel::Elastic {
        ea: E * A * big,
        ei: E * ib,
    };
    f.add_element(n0, n1, 4, 0.0, &col).unwrap();
    f.add_element(n1, n2, 4, 0.0, &beam).unwrap();
    f.add_element(n3, n2, 4, 0.0, &col).unwrap();
    let mut load = vec![0.0; f.ndof()];
    load[3 * n1] = 1.0e5;
    f.static_load(&load, 1).unwrap();
    let drift = f.displacements()[3 * n1];
    let rho = (ib / span) / (ic / h);
    let k_exact = 24.0 * E * ic / h.powi(3) * (6.0 * rho + 1.0) / (6.0 * rho + 4.0);
    let k = 1.0e5 / drift;
    assert!((k - k_exact).abs() / k_exact < 1e-6, "k {k} vs {k_exact}");
}

#[test]
fn frame_g2_corotational_pure_bending_rolls_into_a_half_circle() {
    // Tip moment M = πEI/L bends the cantilever into a half circle of
    // radius L/π: tip at (0, 2L/π), tip rotation π.
    let l = 2.0;
    let (mut f, tip) = cantilever(Geometry::Corotational, l, 20);
    let mut load = vec![0.0; f.ndof()];
    load[3 * tip + 2] = std::f64::consts::PI * E * I / l;
    f.static_load(&load, 40).unwrap();
    let u = f.displacements();
    let (x, y) = (l + u[3 * tip], u[3 * tip + 1]);
    assert!(x.abs() < 5e-3 * l, "tip x {x}");
    assert!(
        (y - 2.0 * l / std::f64::consts::PI).abs() < 5e-3 * l,
        "tip y {y}"
    );
    assert!(
        (u[3 * tip + 2] - std::f64::consts::PI).abs() < 1e-9,
        "tip rotation {}",
        u[3 * tip + 2]
    );
}

#[test]
fn frame_g1_tip_mass_frequency_and_buckling_softening() {
    // Massless cantilever column (along +y) with tip mass m:
    // ω₀ = √(3EI/(mL³)); under axial compression P the lateral frequency
    // softens toward zero at P_cr = π²EI/(4L²) (ω²/ω₀² ≈ 1 − P/P_cr).
    let (l, m) = (3.0, 2000.0);
    let build = || {
        let mut f = Frame2d::new(Geometry::Corotational);
        let mut prev = f.add_node(0.0, 0.0);
        f.fix(prev, [true; 3]);
        for k in 1..=8 {
            let node = f.add_node(0.0, l * f64::from(k) / 8.0);
            f.add_element(prev, node, 5, 0.0, &elastic).unwrap();
            prev = node;
        }
        f.add_mass(prev, m);
        (f, prev)
    };
    let (f0, _) = build();
    let w = f0.natural_frequencies().unwrap();
    let w0 = (3.0 * E * I / (m * l.powi(3))).sqrt();
    // Lowest frequency is lateral bending (axial is far stiffer).
    assert!((w[0] - w0).abs() / w0 < 1e-6, "ω {} vs {w0}", w[0]);
    let pcr = std::f64::consts::PI.powi(2) * E * I / (4.0 * l * l);
    let (mut f1, top) = build();
    let mut load = vec![0.0; f1.ndof()];
    load[3 * top + 1] = -0.5 * pcr;
    f1.static_load(&load, 5).unwrap();
    let w1 = f1.natural_frequencies().unwrap();
    let ratio = (w1[0] / w[0]).powi(2);
    assert!(
        (ratio - 0.5).abs() < 0.03,
        "ω²(P)/ω²(0) = {ratio}, expected ≈ 1 − P/P_cr = 0.5"
    );
}

fn rc() -> SectionModel {
    SectionModel::Fiber(rc_section(0.5, 0.5, 20, 1.5e-3))
}

/// Monotonic moment envelope peak of the RC section at N = 0 over
/// curvatures up to `kappa_max` (section Newton on ε₀ per κ) — an
/// independent oracle for the column's base moment.
fn section_capacity(kappa_max: f64) -> f64 {
    let SectionModel::Fiber(s) = rc() else {
        unreachable!()
    };
    let mut peak = 0.0f64;
    let mut e0 = 0.0;
    for k in 1..=800 {
        let kappa = kappa_max * f64::from(k) / 800.0;
        for _ in 0..60 {
            let r = s.respond(e0, kappa);
            if r.n.abs() < 1.0 {
                break;
            }
            e0 -= r.n / r.tangent[0][0];
        }
        peak = peak.max(s.respond(e0, kappa).m);
    }
    peak
}

#[test]
fn frame_g2_rc_pushover_equilibrium_and_capacity() {
    let l = 3.0;
    let mut f = Frame2d::new(Geometry::Linear);
    let base = f.add_node(0.0, 0.0);
    let top = f.add_node(0.0, l);
    f.fix(base, [true; 3]);
    f.add_element(base, top, 5, 0.0, &rc).unwrap();
    let mut pattern = vec![0.0; f.ndof()];
    pattern[3 * top] = 1.0;
    let targets: Vec<f64> = (1..=60).map(|k| 0.03 * l * f64::from(k) / 60.0).collect();
    let steps = f.displacement_control(&pattern, 3 * top, &targets).unwrap();
    let base_kappa = f.section_deformations(0)[0].1.abs();
    let cap = section_capacity(1.01 * base_kappa);
    let peak_force = steps.iter().fold(0.0f64, |m, s| m.max(s.lambda));
    // Exact equilibrium: base moment = H·L, and it never exceeds the
    // independently computed section capacity.
    let q = f.element_forces(0);
    let h = steps.last().unwrap().lambda;
    assert!(
        (q[1].abs() - h * l).abs() < 1e-6 * h * l,
        "Mᵢ {} vs H·L {}",
        q[1],
        h * l
    );
    assert!(
        peak_force * l <= cap * 1.001,
        "base moment {} above capacity {cap}",
        peak_force * l
    );
    assert!(
        peak_force * l >= 0.85 * cap,
        "pushover peak {} far below capacity {cap}",
        peak_force * l
    );
    // Plasticity concentrates at the base: base curvature ≫ mid-height.
    let d = f.section_deformations(0);
    assert!(d[0].1.abs() > 5.0 * d[2].1.abs(), "{d:?}");
}

#[test]
fn frame_g3_cyclic_pushover_dissipates_energy() {
    let l = 3.0;
    let mut f = Frame2d::new(Geometry::Linear);
    let base = f.add_node(0.0, 0.0);
    let top = f.add_node(0.0, l);
    f.fix(base, [true; 3]);
    f.add_element(base, top, 5, 0.0, &rc).unwrap();
    let mut pattern = vec![0.0; f.ndof()];
    pattern[3 * top] = 1.0;
    let amp = 0.02 * l;
    let mut targets = Vec::new();
    for cycle in [0.5, 1.0] {
        for k in 0..=20 {
            targets.push(cycle * amp * f64::from(k) / 20.0);
        }
        for k in 1..=40 {
            targets.push(cycle * amp * (1.0 - f64::from(k) / 20.0));
        }
        for k in 1..=20 {
            targets.push(cycle * amp * (-1.0 + f64::from(k) / 20.0));
        }
    }
    let steps = f.displacement_control(&pattern, 3 * top, &targets).unwrap();
    let mut work = 0.0;
    let mut prev = (0.0, 0.0);
    for s in &steps {
        let x = s.u[3 * top];
        work += 0.5 * (s.lambda + prev.1) * (x - prev.0);
        prev = (x, s.lambda);
    }
    // Closed cycles back at zero drift: net work = hysteretic dissipation.
    assert!(prev.0.abs() < 1e-12);
    assert!(work > 0.0, "cyclic work {work} must be dissipative");
}

#[test]
fn frame_g1_elastic_sdof_newmark_matches_exact_response() {
    // Cantilever + tip mass under a ground half-sine pulse; compare with the
    // exact undamped SDOF Duhamel solution at the same ω.
    let (l, m) = (3.0, 2000.0);
    let mut f = Frame2d::new(Geometry::Linear);
    let base = f.add_node(0.0, 0.0);
    let top = f.add_node(0.0, l);
    f.fix(base, [true; 3]);
    f.fix(top, [false, true, false]);
    f.add_element(base, top, 5, 0.0, &elastic).unwrap();
    f.add_mass(top, m);
    let w = (3.0 * E * I / (m * l.powi(3))).sqrt();
    let period = 2.0 * std::f64::consts::PI / w;
    let dt = period / 400.0;
    let (ag0, tp) = (2.0, 0.5 * period);
    let nsteps = (2.0 * period / dt) as usize;
    let ground: Vec<f64> = (0..=nsteps)
        .map(|k| {
            let t = k as f64 * dt;
            if t <= tp {
                ag0 * (std::f64::consts::PI * t / tp).sin()
            } else {
                0.0
            }
        })
        .collect();
    let hist = f
        .newmark(&ground, dt, (1.0, 0.0), Rayleigh { a0: 0.0, a1: 0.0 })
        .unwrap();
    // Exact: ü + ω²u = −a_g(t), by numerical Duhamel with fine quadrature.
    let exact = |t: f64| -> f64 {
        let nq = 4000;
        let tau_max = t.min(tp);
        let h = tau_max / f64::from(nq);
        let mut s = 0.0;
        for q in 0..=nq {
            let tau = f64::from(q) * h;
            let wgt = if q == 0 || q == nq { 0.5 } else { 1.0 };
            s += wgt * ag0 * (std::f64::consts::PI * tau / tp).sin() * (w * (t - tau)).sin();
        }
        -s * h / w
    };
    let peak_exact = (0..=nsteps)
        .step_by(10)
        .map(|k| exact(k as f64 * dt).abs())
        .fold(0.0f64, f64::max);
    for k in (0..=nsteps).step_by(25) {
        let t = k as f64 * dt;
        let e = (hist.u[k][3 * top] - exact(t)).abs();
        assert!(
            e < 1e-3 * peak_exact,
            "t {t}: {} vs {}",
            hist.u[k][3 * top],
            exact(t)
        );
    }
    assert!(
        hist.energy_balance_error() < 1e-6,
        "{}",
        hist.energy_balance_error()
    );
}

/// A synthetic broadband ground motion: deterministic sum of sinusoids
/// under a trapezoidal envelope (peak ≈ 0.35 g).
fn ground_motion(dt: f64, duration: f64) -> Vec<f64> {
    let n = (duration / dt) as usize;
    (0..=n)
        .map(|k| {
            let t = k as f64 * dt;
            let env = (t / 1.0).min(1.0) * ((duration - t) / 2.0).clamp(0.0, 1.0);
            let s: f64 = [(1.1, 0.0), (2.3, 1.3), (3.7, 2.1), (5.9, 0.4), (8.3, 2.9)]
                .iter()
                .map(|(fr, ph)| (2.0 * std::f64::consts::PI * fr * t + ph).sin())
                .sum();
            0.35 * 9.81 * env * s / 2.2
        })
        .collect()
}

fn rc_frame(stories: usize, bays: usize) -> (Frame2d, Vec<usize>) {
    let (h, span, floor_mass) = (3.2, 5.0, 40_000.0);
    let mut f = Frame2d::new(Geometry::Corotational);
    let mut grid = vec![vec![0usize; bays + 1]; stories + 1];
    for (s, row) in grid.iter_mut().enumerate() {
        for (b, slot) in row.iter_mut().enumerate() {
            *slot = f.add_node(span * b as f64, h * s as f64);
        }
    }
    for b in 0..=bays {
        f.fix(grid[0][b], [true; 3]);
    }
    let beam = || SectionModel::Fiber(rc_section(0.6, 0.35, 16, 1.2e-3));
    for s in 1..=stories {
        for b in 0..=bays {
            f.add_element(grid[s - 1][b], grid[s][b], 5, 0.0, &rc)
                .unwrap();
            f.add_mass(grid[s][b], floor_mass / (bays + 1) as f64);
        }
        for b in 0..bays {
            f.add_element(grid[s][b], grid[s][b + 1], 4, 0.0, &beam)
                .unwrap();
        }
    }
    let left: Vec<usize> = (0..=stories).map(|s| grid[s][0]).collect();
    (f, left)
}

#[test]
fn frame_g2_multistory_rc_frame_time_history_with_energy_ledger() {
    let (mut f, left) = rc_frame(3, 2);
    // Gravity first (floor weight on the nodes), then the record.
    let mut gravity = vec![0.0; f.ndof()];
    for &n in &left[1..] {
        gravity[3 * n + 1] = -40_000.0 / 3.0 * 9.81;
    }
    f.static_load(&gravity, 4).unwrap();
    let w = f.natural_frequencies().unwrap();
    assert!(w.len() >= 3 && w[0] > 0.0);
    let t1 = 2.0 * std::f64::consts::PI / w[0];
    assert!(
        t1 > 0.1 && t1 < 2.0,
        "fundamental period {t1} s outside the plausible RC band"
    );
    let damping = Rayleigh::from_modes(0.05, w[0], w[2]);
    let dt = 0.01;
    let ground = ground_motion(dt, 8.0);
    let hist = f.newmark(&ground, dt, (1.0, 0.0), damping).unwrap();
    let roof = *left.last().unwrap();
    let peak_roof = hist.u.iter().fold(0.0f64, |m, u| m.max(u[3 * roof].abs()));
    let peak_shear = hist.base_shear.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    assert!(
        peak_roof > 1e-3 && peak_roof < 0.5,
        "roof drift {peak_roof}"
    );
    assert!(peak_shear > 0.0);
    assert!(
        hist.energy_balance_error() < 2e-2,
        "energy balance {}",
        hist.energy_balance_error()
    );
    // Damping and hysteresis both remove energy over the record.
    assert!(*hist.damping_energy.last().unwrap() > 0.0);
}

#[test]
fn frame_g5_time_history_replays_bit_identically() {
    let run = || {
        let (mut f, left) = rc_frame(2, 1);
        let ground = ground_motion(0.01, 2.0);
        let hist = f
            .newmark(&ground, 0.01, (1.0, 0.0), Rayleigh { a0: 0.3, a1: 0.002 })
            .unwrap();
        let roof = *left.last().unwrap();
        hist.u
            .iter()
            .map(|u| u[3 * roof].to_bits())
            .collect::<Vec<u64>>()
    };
    assert_eq!(run(), run());
}

#[test]
fn frame_falsifier_bad_inputs_are_refused() {
    let mut f = Frame2d::new(Geometry::Linear);
    let a = f.add_node(0.0, 0.0);
    let b = f.add_node(0.0, 0.0);
    assert!(
        f.add_element(a, b, 5, 0.0, &elastic).is_err(),
        "zero length"
    );
    let c = f.add_node(1.0, 0.0);
    assert!(
        f.add_element(a, c, 9, 0.0, &elastic).is_err(),
        "unsupported Lobatto count"
    );
    assert!(
        f.add_element(a, 99, 5, 0.0, &elastic).is_err(),
        "unknown node"
    );
    assert!(f.static_load(&[1.0], 1).is_err(), "wrong load length");
}
