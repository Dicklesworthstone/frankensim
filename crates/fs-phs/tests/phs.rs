//! fs-phs conformance battery: admission rejection, discrete-gradient
//! exactness, Gonzalez order-2 on non-quadratic H, interconnection vs
//! monolithic parity + associativity, per-step power balance, mutation
//! visibility through the supply-rate audit, and structure-preserving
//! reduction certified passive by the INDEPENDENT fs-vfit Hamiltonian
//! test.

use fs_math::c64::C64;
use fs_phs::{
    AcousticSection, AcousticTap, MouthFlange, PhsError, PortHamiltonian, QuadraticStorage,
    Storage, ViscothermalPin, acoustic_chain, acoustic_cylinder, acoustic_waveguide,
    common_effort_capacitor, common_effort_dirac, common_effort_star, common_flow_dirac,
    compact_radiation_impedance, discrete_gradient, duffing_oscillator, foster_sqrt_omega_terms,
    helmholtz_resonator, helmholtz_resonator_flow, helmholtz_resonator_radiating, interconnect,
    join_port, kirchhoff_parallel_step, lc_ladder, lc_ladder_terminated, mass_spring_damper,
    modal_bank, modal_bank_ports, moving_end_waveguide, reduce_galerkin, regularized_coulomb,
    series_impedance_ports, slice_linear_taper, spherical_cone, step, step_descriptor, transformer,
    zwikker_kosten_f,
};

fn max_abs(v: &[f64]) -> f64 {
    v.iter().fold(0.0f64, |a, &x| a.max(x.abs()))
}

#[test]
fn admission_rejects_broken_structure() {
    let storage = || Box::new(QuadraticStorage::new(vec![1.0, 0.0, 0.0, 1.0], 2).expect("q"));
    // Non-skew J.
    let bad_j = PortHamiltonian::new(
        2,
        1,
        vec![0.0, 1.0, 1.0, 0.0],
        vec![0.0; 4],
        vec![0.0, 1.0],
        storage(),
    );
    assert!(matches!(bad_j, Err(PhsError::NotSymmetric { what: "J" })));
    // Negative-definite R.
    let bad_r = PortHamiltonian::new(
        2,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, -0.5],
        vec![0.0, 1.0],
        storage(),
    );
    assert!(matches!(bad_r, Err(PhsError::NotPsd { what: "R" })));
    // Asymmetric R.
    let asym_r = PortHamiltonian::new(
        2,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.3, 0.0, 0.5],
        vec![0.0, 1.0],
        storage(),
    );
    assert!(matches!(asym_r, Err(PhsError::NotSymmetric { what: "R" })));
    // Non-PSD Q refuses at the storage.
    assert!(QuadraticStorage::new(vec![1.0, 0.0, 0.0, -1.0], 2).is_err());
}

#[test]
fn discrete_gradient_identity_and_lossless_conservation() {
    // The Gonzalez identity dg.(b-a) = H(b) - H(a) holds exactly.
    let q = QuadraticStorage::new(vec![4.0, 1.0, 1.0, 2.0], 2).expect("q");
    let (a, b) = ([0.3, -1.2], [0.9, 0.4]);
    let dg = discrete_gradient(&q, &a, &b);
    let lhs: f64 = dg
        .iter()
        .zip(b.iter().zip(&a))
        .map(|(&g, (&bi, &ai))| g * (bi - ai))
        .sum();
    let rhs = q.hamiltonian(&b) - q.hamiltonian(&a);
    assert!((lhs - rhs).abs() <= 1.0e-14 * rhs.abs().max(1.0));
    // Undriven lossless LC ladder: H conserved to machine epsilon over
    // thousands of steps (not just small — EXACT up to Newton tol).
    let sys = lc_ladder(4, 0.5, 2.0e-3).expect("lc");
    let mut x = vec![0.0; sys.state_dim()];
    x[0] = 1.0e-3;
    x[3] = 2.0e-4;
    let h0 = sys.hamiltonian(&x);
    let dt = 1.0e-3;
    let mut worst = 0.0f64;
    for _ in 0..2000 {
        let rec = step(&sys, &x, &[0.0], dt).expect("step");
        x = rec.x;
        worst = worst.max((sys.hamiltonian(&x) - h0).abs());
    }
    assert!(
        worst <= 1.0e-10 * h0,
        "lossless H drifted by {worst:.3e} (H0 = {h0:.3e})"
    );
    println!(
        "{{\"suite\":\"fs-phs\",\"case\":\"lossless-conservation\",\"rel_drift\":{:.3e},\"verdict\":\"pass\"}}",
        worst / h0
    );
}

#[test]
fn damped_ledger_matches_exactly() {
    // R > 0: per-step balance residual is solver-zero and the summed
    // ledger reproduces the total H drop.
    let sys = mass_spring_damper(0.02, 800.0, 0.15).expect("msd");
    let mut x = vec![5.0e-3, 0.0];
    let h0 = sys.hamiltonian(&x);
    let dt = 1.0e-4;
    let mut total_dissipated = 0.0;
    let mut worst_residual = 0.0f64;
    let (m, c) = (0.02, 0.15);
    let mut quad_dissipation = 0.0;
    let mut v_prev = x[1] / m;
    for _ in 0..5000 {
        let rec = step(&sys, &x, &[0.0], dt).expect("step");
        worst_residual = worst_residual.max(rec.balance_residual().abs());
        assert!(rec.dissipated >= 0.0, "admitted R cannot un-dissipate");
        total_dissipated += rec.dissipated;
        x = rec.x;
        let v = x[1] / m;
        quad_dissipation += 0.5 * c * (v_prev * v_prev + v * v) * dt;
        v_prev = v;
    }
    let h_end = sys.hamiltonian(&x);
    assert!(worst_residual <= 1.0e-12 * h0);
    // Ledger-vs-H-drop is the accumulated balance residual (same code
    // path — a consistency restatement, kept as a tripwire).
    assert!(
        ((h0 - h_end) - total_dissipated).abs() <= 1.0e-9 * h0,
        "ledger vs H drop mismatch"
    );
    // INDEPENDENT oracle (review finding): the dissipation ledger must
    // match a trapezoidal quadrature of c*v^2 computed from recorded
    // VELOCITIES (states, not StepRecords). Trapezoid error at ~50
    // samples/period is ~0.4%; 2% gate.
    assert!(
        (total_dissipated - quad_dissipation).abs() <= 0.02 * total_dissipated,
        "ledger {total_dissipated:.6e} vs c*v^2 quadrature {quad_dissipation:.6e}"
    );
    println!(
        "{{\"suite\":\"fs-phs\",\"case\":\"damped-ledger\",\"h_drop\":{:.6e},\"dissipated\":{total_dissipated:.6e},\"verdict\":\"pass\"}}",
        h0 - h_end
    );
}

#[test]
fn gonzalez_is_order_two_on_nonquadratic_h() {
    // Richardson on the Duffing oscillator: halving dt must cut the
    // endpoint state error by ~4 (order 2).
    let sys = duffing_oscillator(0.01, 200.0, 5.0e6, 0.0).expect("duffing");
    let x0 = vec![8.0e-3, 0.0];
    let t_end = 0.02;
    let run = |dt: f64| -> Vec<f64> {
        let mut x = x0.clone();
        let steps = (t_end / dt).round() as usize;
        for _ in 0..steps {
            x = step(&sys, &x, &[0.0], dt).expect("step").x;
        }
        x
    };
    let reference = run(t_end / 16384.0);
    let err = |x: &[f64]| -> f64 {
        x.iter()
            .zip(&reference)
            .map(|(&a, &b)| (a - b).abs())
            .fold(0.0f64, f64::max)
    };
    let e1 = err(&run(t_end / 256.0));
    let e2 = err(&run(t_end / 512.0));
    let ratio = e1 / e2;
    assert!(
        (3.5..4.5).contains(&ratio),
        "order-2 Richardson ratio {ratio:.3} (e1 {e1:.3e}, e2 {e2:.3e})"
    );
    // And energy is STILL conserved exactly despite the truncation
    // error in the trajectory (the discrete-gradient point).
    let h0 = sys.hamiltonian(&x0);
    let mut x = x0.clone();
    for _ in 0..200 {
        x = step(&sys, &x, &[0.0], t_end / 200.0).expect("step").x;
    }
    assert!((sys.hamiltonian(&x) - h0).abs() <= 1.0e-10 * h0);
    println!(
        "{{\"suite\":\"fs-phs\",\"case\":\"gonzalez-order2\",\"richardson\":{ratio:.3},\"verdict\":\"pass\"}}"
    );
}

/// Hand-assembled monolithic twin of `interconnect(a, b, [(0, 0)])`
/// for two mass-spring-dampers: block J with the +/- G_a G_b^T
/// coupling, block R, block Q.
fn monolithic_msd_pair(
    (m1, k1, c1): (f64, f64, f64),
    (m2, k2, c2): (f64, f64, f64),
) -> PortHamiltonian {
    let n = 4;
    let q = vec![
        k1,
        0.0,
        0.0,
        0.0, //
        0.0,
        1.0 / m1,
        0.0,
        0.0, //
        0.0,
        0.0,
        k2,
        0.0, //
        0.0,
        0.0,
        0.0,
        1.0 / m2,
    ];
    // J: intra blocks [[0,1],[-1,0]] plus coupling between the two
    // momentum states (g_a = e_p1, g_b = e_p2): -1 at (p1, p2), +1 at
    // (p2, p1).
    let mut j = vec![0.0; n * n];
    j[1] = 1.0;
    j[n] = -1.0;
    j[2 * n + 3] = 1.0;
    j[3 * n + 2] = -1.0;
    j[n + 3] = -1.0;
    j[3 * n + 1] = 1.0;
    let mut r = vec![0.0; n * n];
    r[n + 1] = c1;
    r[3 * n + 3] = c2;
    let g = vec![]; // no external ports after pairing
    let storage = Box::new(QuadraticStorage::new(q, n).expect("q"));
    PortHamiltonian::new(n, 0, j, r, g, storage).expect("monolithic")
}

#[test]
fn helmholtz_resonator_interconnects_with_a_modal_bank() {
    let cavity = helmholtz_resonator(0.0135, 0.043, 0.02, 1.2, 343.0, 0.0).expect("helmholtz");
    let plate = modal_bank(&[2.0e3], &[0.0], &[1.0]).expect("plate");
    let coupled = interconnect(plate, cavity, &[(0, 0)]).expect("plate+cavity");
    assert_eq!(coupled.state_dim(), 4);
    assert_eq!(coupled.port_dim(), 0);
    let mut x = vec![1.0e-4, 0.0, 0.0, 0.0];
    let h0 = coupled.hamiltonian(&x);
    for _ in 0..200 {
        x = step(&coupled, &x, &[], 1.0e-5).expect("step").x;
    }
    let h1 = coupled.hamiltonian(&x);
    assert!(
        (h1 - h0).abs() <= 1.0e-8 * h0.abs().max(1.0e-18),
        "lossless plate+cavity must hold H ({h1} vs {h0})"
    );
    assert!(helmholtz_resonator(0.0, 0.01, 0.01, 1.2, 343.0, 0.0).is_err());
}

#[test]
fn flow_driven_helmholtz_pressure_tracks_injected_volume() {
    let cav = helmholtz_resonator_flow(0.01, 0.02, 0.03, 1.2, 343.0, 0.0).expect("flow");
    let mut x = vec![0.0, 0.0];
    let dt = 1.0e-5;
    for _ in 0..20 {
        x = step(&cav, &x, &[1.0e-4], dt).expect("step").x;
    }
    let p = cav.output(&x)[0];
    assert!(p > 0.0, "injected volume must raise cavity pressure, p={p}");
    assert!(helmholtz_resonator_flow(0.0, 0.02, 0.03, 1.2, 343.0, 0.0).is_err());
}

#[test]
fn interconnection_matches_monolithic_trajectories() {
    let pa = (0.02, 800.0, 0.05);
    let pb = (0.05, 300.0, 0.02);
    let a = mass_spring_damper(pa.0, pa.1, pa.2).expect("a");
    let b = mass_spring_damper(pb.0, pb.1, pb.2).expect("b");
    let coupled = interconnect(a, b, &[(0, 0)]).expect("interconnect");
    assert_eq!(coupled.state_dim(), 4);
    assert_eq!(coupled.port_dim(), 0);
    let mono = monolithic_msd_pair(pa, pb);
    let mut xc = vec![4.0e-3, 0.0, -2.0e-3, 0.0];
    let mut xm = xc.clone();
    let dt = 5.0e-5;
    for _ in 0..4000 {
        xc = step(&coupled, &xc, &[], dt).expect("c").x;
        xm = step(&mono, &xm, &[], dt).expect("m").x;
    }
    let dev = xc
        .iter()
        .zip(&xm)
        .map(|(&p, &q)| (p - q).abs())
        .fold(0.0f64, f64::max);
    let scale = max_abs(&xm).max(1.0e-12);
    assert!(
        dev <= 1.0e-10 * scale,
        "interconnected vs monolithic drifted {dev:.3e}"
    );
    println!(
        "{{\"suite\":\"fs-phs\",\"case\":\"interconnect-monolithic\",\"rel_dev\":{:.3e},\"verdict\":\"pass\"}}",
        dev / scale
    );
}

#[test]
fn interconnection_is_associative_on_structure() {
    // ((A + B) + C) and (A + (B + C)) with the same disjoint pairings
    // must produce IDENTICAL structure matrices (state order [A B C]).
    let mk = |m: f64, k: f64, c: f64, ports: usize| -> PortHamiltonian {
        // msd with `ports` copies of the force port.
        let q = vec![k, 0.0, 0.0, 1.0 / m];
        let mut g = vec![0.0; 2 * ports];
        for p in 0..ports {
            g[ports + p] = 1.0;
        }
        PortHamiltonian::new(
            2,
            ports,
            vec![0.0, 1.0, -1.0, 0.0],
            vec![0.0, 0.0, 0.0, c],
            g,
            Box::new(QuadraticStorage::new(q, 2).expect("q")),
        )
        .expect("sys")
    };
    // A: ports {to B}; B: ports {to A, to C}; C: ports {to B}.
    let left = interconnect(
        interconnect(
            mk(0.02, 800.0, 0.01, 1),
            mk(0.05, 300.0, 0.02, 2),
            &[(0, 0)],
        )
        .expect("ab"),
        mk(0.03, 500.0, 0.03, 1),
        &[(0, 0)],
    )
    .expect("ab_c");
    let right = interconnect(
        mk(0.02, 800.0, 0.01, 1),
        interconnect(
            mk(0.05, 300.0, 0.02, 2),
            mk(0.03, 500.0, 0.03, 1),
            &[(1, 0)],
        )
        .expect("bc"),
        &[(0, 0)],
    )
    .expect("a_bc");
    let (jl, rl, gl) = left.structure();
    let (jr, rr, gr) = right.structure();
    assert_eq!(jl, jr, "J not associative");
    assert_eq!(rl, rr, "R not associative");
    assert_eq!(gl, gr, "G not associative");
}

#[test]
fn driven_power_balance_and_supply_audit() {
    let sys = mass_spring_damper(0.02, 800.0, 0.15).expect("msd");
    let mut x = vec![0.0, 0.0];
    let dt = 1.0e-4;
    let mut worst_defect = f64::NEG_INFINITY;
    for i in 0i32..3000 {
        let u = [0.4 * fs_math::det::sin(2.0 * core::f64::consts::PI * 150.0 * f64::from(i) * dt)];
        let rec = step(&sys, &x, &u, dt).expect("step");
        assert!(rec.balance_residual().abs() <= 1.0e-10);
        // Passive: energy gained never exceeds energy supplied.
        worst_defect = worst_defect.max(rec.supply_defect());
        x = rec.x;
    }
    assert!(
        worst_defect <= 1.0e-12,
        "supply audit violated by {worst_defect:.3e}"
    );
}

#[test]
fn mutations_caught_by_supply_audit() {
    // Symmetric J (energy-pumping "gyrator gone wrong") smuggled past
    // admission via from_raw_parts: the supply audit MUST fire.
    let q = QuadraticStorage::new(vec![800.0, 0.0, 0.0, 50.0], 2).expect("q");
    let sym_j = PortHamiltonian::from_raw_parts(
        2,
        0,
        vec![0.0, 1.0, 1.0, 0.0],
        vec![0.0; 4],
        vec![],
        Box::new(q.clone()),
    );
    let mut x = vec![1.0e-2, 1.0e-2];
    let mut fired = false;
    for _ in 0..2000 {
        let rec = step(&sym_j, &x, &[], 1.0e-4).expect("step");
        if rec.supply_defect() > 1.0e-12 {
            fired = true;
            break;
        }
        x = rec.x;
    }
    assert!(fired, "symmetrized J must violate the supply audit");
    // Sign-flipped R: undriven energy GROWS.
    let neg_r = PortHamiltonian::from_raw_parts(
        2,
        0,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, -0.15],
        vec![],
        Box::new(q.clone()),
    );
    let mut x = vec![1.0e-2, 0.0];
    let mut fired = false;
    for _ in 0..2000 {
        let rec = step(&neg_r, &x, &[], 1.0e-4).expect("step");
        if rec.supply_defect() > 1.0e-12 {
            fired = true;
            break;
        }
        x = rec.x;
    }
    assert!(fired, "sign-flipped R must violate the supply audit");
    // Broken interconnection map (one-sided coupling): the composite J
    // is not skew, so ADMISSION refuses it by name — the same class
    // fs-couple's cross-row residual caught dynamically.
    let mut j = vec![0.0; 16];
    j[1] = 1.0;
    j[4] = -1.0;
    j[11] = 1.0;
    j[14] = -1.0;
    j[4 + 3] = -1.0; // coupling one way only
    let q4 = QuadraticStorage::new(
        vec![
            800.0, 0.0, 0.0, 0.0, //
            0.0, 50.0, 0.0, 0.0, //
            0.0, 0.0, 300.0, 0.0, //
            0.0, 0.0, 0.0, 20.0,
        ],
        4,
    )
    .expect("q4");
    let broken = PortHamiltonian::new(4, 0, j, vec![0.0; 16], vec![], Box::new(q4));
    assert!(matches!(broken, Err(PhsError::NotSymmetric { what: "J" })));
    println!("{{\"suite\":\"fs-phs\",\"case\":\"mutations\",\"verdict\":\"pass\"}}");
}

#[test]
#[allow(clippy::too_many_lines)] // one linear pipeline: reduce -> certify -> H-error
fn reduction_preserves_structure_and_passivity() {
    // 24-mode bank -> keep the 8 modes carrying the drive; the reduced
    // system must re-admit (structural preservation) and its impedance
    // must certify passive under the INDEPENDENT fs-vfit Hamiltonian
    // test (the vector-fitting bead's machinery).
    let nm = 24;
    let omegas: Vec<f64> = (0..nm)
        .map(|i| 2.0 * core::f64::consts::PI * 85.0f64.mul_add(i as f64, 100.0))
        .collect();
    let zetas: Vec<f64> = (0..nm).map(|i| 0.0005f64.mul_add(i as f64, 0.01)).collect();
    let drive: Vec<f64> = (0..nm)
        .map(|i| if i < 8 { 1.0 - 0.05 * i as f64 } else { 0.02 })
        .collect();
    let full = modal_bank(&omegas, &zetas, &drive).expect("bank");
    // Basis: unit vectors of the first 8 modes' (q, p) states.
    let (n, k) = (2 * nm, 16);
    let mut v = vec![0.0; n * k];
    for mode in 0..8 {
        v[(2 * mode) * k + 2 * mode] = 1.0;
        v[(2 * mode + 1) * k + 2 * mode + 1] = 1.0;
    }
    let red = reduce_galerkin(&full, &v, k).expect("reduced admits — structure preserved");
    assert_eq!(red.state_dim(), k);
    // H-error at t = 0 for a state INSIDE the basis is zero (exact
    // certified statement); for a state with energy outside, the
    // deficit equals that energy exactly.
    let mut x_in = vec![0.0; n];
    x_in[0] = 1.0e-3;
    let xr: Vec<f64> = (0..k)
        .map(|c| (0..n).map(|l| v[l * k + c] * x_in[l]).sum())
        .collect();
    assert!((red.hamiltonian(&xr) - full.hamiltonian(&x_in)).abs() <= 1.0e-15);
    // Impedance samples of the reduced system via its linear dynamics,
    // then fs-vfit certification.
    let (jm, rm, gm) = red.structure();
    let a_mat: Vec<f64> = {
        // A = (J - R) Q with Q recovered from unit gradients.
        let mut qm = vec![0.0; k * k];
        for col in 0..k {
            let mut basis = vec![0.0; k];
            basis[col] = 1.0;
            let mut grad = vec![0.0; k];
            // Reduced storage gradient IS linear (quadratic case).
            red_gradient(&red, &basis, &mut grad);
            for row in 0..k {
                qm[row * k + col] = grad[row];
            }
        }
        let mut a = vec![0.0; k * k];
        for i in 0..k {
            for l in 0..k {
                let mut acc = 0.0;
                for t in 0..k {
                    acc += (jm[i * k + t] - rm[i * k + t]) * qm[t * k + l];
                }
                a[i * k + l] = acc;
            }
        }
        // Fold Q into C as well: y = G^T Q x.
        let mut cq = vec![0.0; k];
        for l in 0..k {
            let mut acc = 0.0;
            for t in 0..k {
                acc += gm[t] * qm[t * k + l];
            }
            cq[l] = acc;
        }
        // Sample H(i w) = Cq (i w I - A)^{-1} G on a band and fit.
        let omega_grid: Vec<f64> = (0..300)
            .map(|i| {
                let t = f64::from(i) / 299.0;
                let lo = 300.0f64;
                let hi = 6000.0f64;
                2.0 * core::f64::consts::PI * lo * fs_math::det::exp(t * fs_math::det::ln(hi / lo))
            })
            .collect();
        let mut h = Vec::with_capacity(omega_grid.len());
        for &w in &omega_grid {
            let mut m = vec![C64::ZERO; k * k];
            for i in 0..k {
                for l in 0..k {
                    m[i * k + l] = C64::from_re(-a[i * k + l]);
                }
                m[i * k + i] = m[i * k + i] + C64::new(0.0, w);
            }
            let lu = fs_la::eigen_complex::lu_complex(&m, k).expect("lu");
            let mut xcol: Vec<C64> = gm.iter().map(|&x| C64::from_re(x)).collect();
            lu.solve(&mut xcol);
            let mut acc = C64::ZERO;
            for (ci, xi) in cq.iter().zip(&xcol) {
                acc = acc + xi.scale(*ci);
            }
            h.push(acc);
        }
        let fit = fs_vfit::vector_fit(
            &omega_grid,
            &h,
            &fs_vfit::FitOptions {
                fit_e: false,
                ..fs_vfit::FitOptions::new(16)
            },
        )
        .expect("fit");
        let report = fs_vfit::passivity::check_passivity(
            &fit.model,
            (omega_grid[0], *omega_grid.last().expect("grid")),
        )
        .expect("check");
        assert!(
            report.passive,
            "reduced modal bank must certify passive (worst Re = {:?})",
            report.worst
        );
        println!(
            "{{\"suite\":\"fs-phs\",\"case\":\"reduction-passivity\",\"class\":\"{:?}\",\"worst_re\":{:.3e},\"verdict\":\"pass\"}}",
            report.class, report.worst.0
        );
        a
    };
    let _ = a_mat;
    // Realized H-error on a validation input (drive within the kept
    // modes): logged, and gated by an authored envelope measured at
    // authoring time.
    let dt = 2.0e-5;
    let mut xf = vec![0.0; n];
    let mut xr = vec![0.0; k];
    let mut worst_h_dev = 0.0f64;
    let mut h_scale = 0.0f64;
    for i in 0i32..2000 {
        let u = [0.1 * fs_math::det::sin(2.0 * core::f64::consts::PI * 180.0 * f64::from(i) * dt)];
        xf = step(&full, &xf, &u, dt).expect("full").x;
        xr = step(&red, &xr, &u, dt).expect("red").x;
        let hf = full.hamiltonian(&xf);
        let hr = red.hamiltonian(&xr);
        h_scale = h_scale.max(hf);
        worst_h_dev = worst_h_dev.max((hf - hr).abs());
    }
    // Authored: the truncated modes carry only the residual 0.02 drive
    // weights; measured H deviation was ~2e-3 of peak H at authoring.
    // 2% is the envelope with wide headroom; a-priori trajectory
    // bounds are a recorded no-claim.
    assert!(
        worst_h_dev <= 0.02 * h_scale,
        "reduction H-error {worst_h_dev:.3e} above envelope ({h_scale:.3e} peak H)"
    );
    println!(
        "{{\"suite\":\"fs-phs\",\"case\":\"reduction-h-error\",\"realized\":{worst_h_dev:.3e},\"peak_h\":{h_scale:.3e},\"verdict\":\"pass\"}}"
    );
}

/// Gradient of a (known-quadratic) reduced system via its Storage —
/// helper for the impedance assembly above.
fn red_gradient(sys: &PortHamiltonian, x: &[f64], out: &mut [f64]) {
    // PortHamiltonian doesn't expose storage directly; use output-free
    // probing: gradient = d/dx H at x by central differences would be
    // inexact — instead use the public output() path? The cleanest
    // honest route: H is quadratic, so grad_i = H(x + e_i) - H(x - e_i)
    // over 2 is EXACT for quadratic H with unit steps scaled small
    // enough to avoid roundoff — use symmetric difference with h = 1.0
    // exactness of quadratics: grad(x) . e_i = [H(x + h e_i) - H(x - h
    // e_i)] / (2h) exactly for ANY h on a quadratic.
    let n = x.len();
    for i in 0..n {
        let mut xp = x.to_vec();
        let mut xm = x.to_vec();
        xp[i] += 1.0;
        xm[i] -= 1.0;
        out[i] = 0.5 * (sys.hamiltonian(&xp) - sys.hamiltonian(&xm));
    }
}

/// Reduction linearity probe at TINY gradient scales: a nonlinear
/// (Duffing) storage whose q-row gradients sit far below any absolute
/// floor must still REFUSE quadratic-only reduction, even with an O(1)
/// momentum row alongside (the mixed-scale masking case); a quadratic
/// system at the same tiny scale must still reduce.
#[test]
fn reduction_refuses_tiny_scale_nonlinear_storage_and_admits_tiny_quadratic() {
    let identity = vec![1.0, 0.0, 0.0, 1.0];
    let duffing = fs_phs::duffing_oscillator(1.0, 1.0e-10, 1.0e-9, 0.0).expect("duffing");
    assert!(matches!(
        fs_phs::reduce_galerkin(&duffing, &identity, 2),
        Err(fs_phs::PhsError::Dimension { .. })
    ));
    let msd = fs_phs::mass_spring_damper(1.0, 1.0e-10, 0.0).expect("msd");
    let reduced = fs_phs::reduce_galerkin(&msd, &identity, 2).expect("tiny quadratic reduces");
    assert_eq!(reduced.state_dim(), 2);
}

#[test]
fn compact_radiation_flanged_resists_more_than_unflanged() {
    let (ru, xu) = compact_radiation_impedance(1.2, 343.0, 0.02, 2.0e3, MouthFlange::Unflanged)
        .expect("unflanged");
    let (rf, xf) = compact_radiation_impedance(1.2, 343.0, 0.02, 2.0e3, MouthFlange::Flanged)
        .expect("flanged");
    assert!(ru > 0.0 && rf > ru, "flanged R {rf} vs unflanged {ru}");
    assert!(xu < 0.0 && xf < xu, "flanged mass load {xf} vs {xu}");
    assert!(compact_radiation_impedance(1.2, 343.0, 0.2, 2.0e4, MouthFlange::Flanged).is_err());
}

#[test]
fn radiating_helmholtz_dissipates_stored_energy() {
    let sys = helmholtz_resonator_radiating(0.002, 0.012, 0.03, 1.2, 343.0, MouthFlange::Flanged)
        .expect("radiating");
    let mut x = vec![1.0e-6, 0.0];
    let h0 = sys.hamiltonian(&x);
    for _ in 0..4000 {
        x = step(&sys, &x, &[0.0], 5.0e-5).expect("step").x;
    }
    let h1 = sys.hamiltonian(&x);
    assert!(
        h1 < 0.85 * h0,
        "mouth radiation must drain H ({h1} vs {h0})"
    );
}

#[test]
fn series_flow_cavities_add_pressure() {
    let one = helmholtz_resonator_flow(0.01, 0.02, 0.03, 1.2, 343.0, 0.0).expect("one");
    let series = series_impedance_ports(
        helmholtz_resonator_flow(0.01, 0.02, 0.03, 1.2, 343.0, 0.0).expect("a"),
        helmholtz_resonator_flow(0.01, 0.02, 0.03, 1.2, 343.0, 0.0).expect("b"),
    )
    .expect("series");
    let dt = 1.0e-5;
    let mut x1 = vec![0.0, 0.0];
    let mut xs = vec![0.0; 4];
    for _ in 0..30 {
        x1 = step(&one, &x1, &[1.0e-4], dt).expect("1").x;
        xs = step(&series, &xs, &[1.0e-4], dt).expect("s").x;
    }
    let p1 = one.output(&x1)[0];
    let ps = series.output(&xs)[0];
    assert!(p1 > 0.0);
    assert!(
        (ps / p1 - 2.0).abs() < 0.05,
        "series impedances must add pressures ({ps} vs 2*{p1})"
    );
    assert!(matches!(
        series_impedance_ports(mass_spring_damper(0.01, 10.0, 0.0).expect("a"), {
            let q = vec![10.0, 0.0, 0.0, 100.0];
            PortHamiltonian::new(
                2,
                2,
                vec![0.0, 1.0, -1.0, 0.0],
                vec![0.0; 4],
                vec![0.0, 1.0, 0.0, 1.0],
                Box::new(QuadraticStorage::new(q, 2).expect("q")),
            )
            .expect("two-port")
        },),
        Err(PhsError::BadPortPairing)
    ));
}

#[test]
fn stick_slip_on_a_modal_string_locks_the_bow() {
    assert!(regularized_coulomb(0.4, 1.0, 1.0, 0.01) < 0.0);
    assert!(
        (regularized_coulomb(0.4, 1.0, 0.3, 0.02) + regularized_coulomb(0.4, 1.0, -0.3, 0.02))
            .abs()
            < 1.0e-15
    );

    // 1-DOF slider: stick when |k q| < μN, then slip. Soft enough
    // that a laboratory bow force actually holds.
    let slider = mass_spring_damper(0.02, 80.0, 0.05).expect("slider");
    let v_bow = 0.08;
    let dt = 1.0e-4;
    let mut x = vec![0.0, 0.0];
    let mut stuck = 0usize;
    let mut seen = 0usize;
    let mut slips = 0usize;
    let mut last_stuck = true;
    for i in 0..8_000 {
        let v = slider.output(&x)[0];
        let f = regularized_coulomb(0.5, 4.0, v - v_bow, 0.003);
        x = step(&slider, &x, &[f], dt).expect("step").x;
        if i > 2_000 {
            seen += 1;
            let holding = (v - v_bow).abs() < 0.02;
            if holding {
                stuck += 1;
            } else if last_stuck {
                slips += 1;
            }
            last_stuck = holding;
        }
    }
    let frac = stuck as f64 / seen as f64;
    assert!(
        frac > 0.15 && frac < 0.98,
        "Coulomb port must stick a Helmholtz-like fraction of the cycle (stuck {frac})"
    );
    assert!(slips >= 1, "stick-slip must break at least once");

    // Same port on a modal string bowed at 1/4: even modes are silent
    // in a linear pluck at that point? No — sin(n π/4) is 1 at n=2.
    // The nonlinear force just has to put energy into mode 2.
    let omega1 = 2.0 * core::f64::consts::PI * 40.0;
    let omegas: Vec<f64> = (1..=4).map(|n| omega1 * n as f64).collect();
    let zetas = vec![0.004; 4];
    let drive: Vec<f64> = (1..=4)
        .map(|n| fs_math::det::sin(n as f64 * core::f64::consts::PI * 0.25))
        .collect();
    let string = modal_bank(&omegas, &zetas, &drive).expect("string");
    let mut xs = vec![0.0; 8];
    for _ in 0..4_000 {
        let v = string.output(&xs)[0];
        let f = regularized_coulomb(0.6, 8.0, v - 0.05, 0.004);
        xs = step(&string, &xs, &[f], 5.0e-5).expect("string").x;
    }
    let e1 = 0.5 * (xs[1] * xs[1] + omegas[0] * omegas[0] * xs[0] * xs[0]);
    let e2 = 0.5 * (xs[3] * xs[3] + omegas[1] * omegas[1] * xs[2] * xs[2]);
    assert!(
        e1 > 0.0 && e2 > 0.05 * e1,
        "even string mode from friction ({e2} vs {e1})"
    );
}

#[test]
fn terminated_ladder_radiates_while_lossless_holds() {
    let live = lc_ladder(6, 1.0e-3, 2.0e-6).expect("live");
    let load = lc_ladder_terminated(6, 1.0e-3, 2.0e-6, 8.0).expect("load");
    let mut xl = vec![0.0; 12];
    let mut xd = vec![0.0; 12];
    xl[0] = 1.0e-4;
    xd[0] = 1.0e-4;
    let h0 = live.hamiltonian(&xl);
    let dt = 2.0e-5;
    for _ in 0..2500 {
        xl = step(&live, &xl, &[0.0], dt).expect("live").x;
        xd = step(&load, &xd, &[0.0], dt).expect("load").x;
    }
    let hl = live.hamiltonian(&xl);
    let hd = load.hamiltonian(&xd);
    assert!(
        (hl - h0).abs() <= 1.0e-8 * h0.abs().max(1.0e-18),
        "lossless ladder must hold H ({hl} vs {h0})"
    );
    assert!(
        hd < 0.5 * h0,
        "terminated ladder must radiate ({hd} vs {h0})"
    );
    assert!(lc_ladder_terminated(6, 1.0e-3, 2.0e-6, -1.0).is_err());
}

#[test]
fn common_effort_capacitor_shares_pressure() {
    let c = common_effort_capacitor(2.0e-8).expect("C");
    let mut x = vec![0.0];
    x = step(&c, &x, &[1.0e-4, 2.0e-4], 1.0e-4).expect("step").x;
    let y = c.output(&x);
    assert_eq!(y.len(), 2);
    assert!(
        (y[0] - y[1]).abs() < 1.0e-15,
        "ports must share p ({:?})",
        y
    );
    assert!(y[0] > 0.0, "injected volume must raise pressure");
    assert!(common_effort_capacitor(0.0).is_err());
}

#[test]
fn kirchhoff_parallel_splits_flow_at_common_pressure() {
    let a = helmholtz_resonator_flow(0.01, 0.02, 0.03, 1.2, 343.0, 0.0).expect("a");
    let b = helmholtz_resonator_flow(0.01, 0.02, 0.03, 1.2, 343.0, 0.0).expect("b");
    let xa = vec![0.0, 0.0];
    let xb = vec![0.0, 0.0];
    let u_ext = 2.0e-4;
    let (ra, rb) = kirchhoff_parallel_step(&a, &xa, &b, &xb, u_ext, 1.0e-4).expect("join");
    assert!(
        (ra.y[0] - rb.y[0]).abs() < 1.0e-8 * (1.0 + ra.y[0].abs()),
        "Kirchhoff must share pressure ({:?} vs {:?})",
        ra.y[0],
        rb.y[0]
    );
    // Identical branches take half the flow; one cavity at U_ext/2
    // must reprint the same pressure.
    let solo = step(&a, &xa, &[0.5 * u_ext], 1.0e-4).expect("solo");
    assert!((ra.y[0] - solo.y[0]).abs() < 1.0e-8 * (1.0 + solo.y[0].abs()));
}

#[test]
fn common_effort_dirac_is_the_kirchhoff_dae() {
    let a = helmholtz_resonator_flow(0.01, 0.02, 0.03, 1.2, 343.0, 0.0).expect("a");
    let b = helmholtz_resonator_flow(0.01, 0.02, 0.03, 1.2, 343.0, 0.0).expect("b");
    let sys = common_effort_dirac(
        helmholtz_resonator_flow(0.01, 0.02, 0.03, 1.2, 343.0, 0.0).expect("a2"),
        helmholtz_resonator_flow(0.01, 0.02, 0.03, 1.2, 343.0, 0.0).expect("b2"),
    )
    .expect("dirac");
    assert_eq!(sys.state_dim(), 5);
    assert_eq!(sys.differential_dim(), 4);
    assert_eq!(sys.port_dim(), 1);
    let n = sys.state_dim();
    let j = sys.dirac_j();
    for i in 0..n {
        for k in 0..n {
            assert!(
                (j[i * n + k] + j[k * n + i]).abs() < 1.0e-14,
                "composite J must be skew at ({i},{k})"
            );
        }
    }
    let u_ext = 2.0e-4;
    let dt = 1.0e-4;
    let rec = step_descriptor(&sys, &[0.0; 5], &[u_ext], dt).expect("descriptor step");
    assert_eq!(rec.x.len(), 5);
    let xa = vec![0.0, 0.0];
    let xb = vec![0.0, 0.0];
    let (ra, rb) = kirchhoff_parallel_step(&a, &xa, &b, &xb, u_ext, dt).expect("split");
    assert!(
        (rec.y[0] - ra.y[0]).abs() < 1.0e-6 * (1.0 + ra.y[0].abs()),
        "descriptor pressure {} must reprint the Newton split {}",
        rec.y[0],
        ra.y[0]
    );
    assert!((ra.y[0] - rb.y[0]).abs() < 1.0e-8 * (1.0 + ra.y[0].abs()));
    let e_a = a.output(&rec.x[..2]);
    let e_b = b.output(&rec.x[2..4]);
    assert!(
        (e_a[0] - e_b[0]).abs() < 1.0e-6 * (1.0 + e_a[0].abs()),
        "algebraic row must enforce p_a = p_b ({e_a:?} vs {e_b:?})"
    );
    assert!(common_effort_dirac(a, common_effort_capacitor(1.0e-8).expect("two-port")).is_err());
}

#[test]
fn kirchhoff_star_of_three_flow_cavities_splits_volume() {
    let mk = || helmholtz_resonator_flow(0.01, 0.02, 0.03, 1.2, 343.0, 0.0).expect("cav");
    let sys = common_effort_star(vec![mk(), mk(), mk()]).expect("star");
    assert_eq!(sys.state_dim(), 2 * 3 + 2);
    assert_eq!(sys.differential_dim(), 6);
    let rec = step_descriptor(&sys, &vec![0.0; sys.state_dim()], &[3.0e-4], 1.0e-4).expect("step");
    let solo = step(&mk(), &[0.0, 0.0], &[1.0e-4], 1.0e-4).expect("solo");
    assert!(
        (rec.y[0] - solo.y[0]).abs() < 1.0e-6 * (1.0 + solo.y[0].abs()),
        "identical branches must reprint the 1/3 flow pressure ({} vs {})",
        rec.y[0],
        solo.y[0]
    );
}

#[test]
fn transformer_joins_a_mass_to_a_compliance() {
    let mass = mass_spring_damper(0.02, 0.0, 0.0).expect("mass");
    let cav = helmholtz_resonator_flow(0.01, 0.02, 0.03, 1.2, 343.0, 0.0).expect("cav");
    let sys = transformer(mass, cav, 0, 0, 4.0e-4).expect("area");
    assert_eq!(sys.port_dim(), 0);
    let mut x = vec![0.0; sys.state_dim()];
    x[0] = 1.0e-4;
    let h0 = sys.hamiltonian(&x);
    let rec = step(&sys, &x, &[], 1.0e-4).expect("closed step");
    assert!(
        (rec.delta_h).abs() <= 1.0e-9 * h0.abs().max(1.0e-18),
        "lossless transformer must hold H (ΔH={}, H0={h0})",
        rec.delta_h
    );
}

#[test]
fn moving_end_waveguide_dirac_joins_a_mass() {
    let string = moving_end_waveguide(2, 0.65, 70.0, 5.0e-3, &[0.0, 0.0]).expect("string");
    let y0 = string.output(&vec![0.0; string.state_dim()]);
    assert_eq!(y0.len(), 1);
    let mass = mass_spring_damper(0.05, 2.0e4, 0.0).expect("bridge mass");
    let joined = common_flow_dirac(string, mass).expect("1-junction");
    let rec = step_descriptor(&joined, &vec![0.0; joined.state_dim()], &[1.0], 1.0e-5)
        .expect("join step");
    assert!(rec.y[0].is_finite());
    assert!(rec.x.iter().any(|v| v.abs() > 0.0));
}

#[test]
fn three_phs_string_plate_cavity_is_a_dirac_star_plus_transformer() {
    let string = moving_end_waveguide(1, 0.65, 70.0, 5.0e-3, &[0.0]).expect("string");
    let plate = modal_bank_ports(&[800.0], &[0.0], &[&[2.0], &[0.02]]).expect("plate");
    let cav = helmholtz_resonator_flow(0.002, 0.01, 0.02, 1.2, 343.0, 0.0).expect("cav");
    let plate_cav = transformer(plate, cav, 1, 0, 1.0).expect("area in G");
    assert_eq!(plate_cav.port_dim(), 1);
    let sys = common_flow_dirac(string, plate_cav).expect("three-pHS");
    assert_eq!(sys.port_dim(), 1);
    let rec = step_descriptor(&sys, &vec![0.0; sys.state_dim()], &[0.5], 2.0e-5).expect("step");
    assert!(rec.y[0].is_finite());
    assert!(rec.solver_residual.is_finite());
}

#[test]
fn join_port_keeps_the_leftover_and_holds_energy_when_closed() {
    let mass = mass_spring_damper(0.02, 80.0, 0.0).expect("mass");
    let spring = mass_spring_damper(0.03, 120.0, 0.0).expect("spring");
    let closed = join_port(mass, spring, 0, 0).expect("closed 1-junction");
    assert_eq!(closed.port_dim(), 0);
    let mut x = vec![0.0; closed.state_dim()];
    x[0] = 1.0e-4;
    let h0 = closed.hamiltonian(&x);
    let rec = step_descriptor(&closed, &x, &[], 1.0e-4).expect("closed step");
    assert!(
        (closed.hamiltonian(&rec.x) - h0).abs() <= 1.0e-8 * (1.0 + h0.abs()),
        "closed join_port must hold H"
    );
    let two = modal_bank_ports(&[40.0], &[0.0], &[&[1.0], &[0.5]]).expect("two-port");
    let load = mass_spring_damper(0.01, 50.0, 0.0).expect("load");
    let open = join_port(two, load, 0, 0).expect("leftover");
    assert_eq!(open.port_dim(), 1);
    let rec = step_descriptor(&open, &vec![0.0; open.state_dim()], &[0.2], 1.0e-4).expect("drive");
    assert!(rec.y[0].is_finite());
    assert!(
        join_port(
            mass_spring_damper(0.02, 80.0, 0.0).expect("a"),
            mass_spring_damper(0.03, 120.0, 0.0).expect("b"),
            1,
            0,
        )
        .is_err()
    );
}

#[test]
fn acoustic_cylinder_rings_at_the_quarter_wave() {
    let l = 0.34;
    let c = 343.0;
    let sys = acoustic_cylinder(l, 0.012, 1.2, c, 8, false, 1).expect("closed");
    let dt = 1.0 / 8_000.0;
    let mut x = vec![0.0; sys.state_dim()];
    let mut p = Vec::new();
    for i in 0..640 {
        let u = if i < 16 {
            2.0e-5 * (core::f64::consts::PI * i as f64 / 16.0).sin()
        } else {
            0.0
        };
        let rec = step(&sys, &x, &[u], dt).expect("step");
        p.push(rec.y[0]);
        x = rec.x;
    }
    // Current drive into the first compliance + a lossless last
    // inertance is the quarter-wave (open-mouth / stop) family,
    // period 4L/c. Half-wave closed-closed would need U=0 at both
    // ends, which this Cauer form does not impose.
    let want = 4.0 * l / c;
    let got = dominant_zero_period(&p, dt);
    assert!(
        (got - want).abs() / want < 0.15,
        "cylinder period {got} vs 4L/c {want}"
    );
    assert!(acoustic_cylinder(l, 0.012, 1.2, c, 1, false, 1).is_err());
    let two = acoustic_cylinder(l, 0.012, 1.2, c, 8, false, 2).expect("two inlets");
    assert_eq!(two.port_dim(), 2);
}

#[test]
fn open_tap_raises_the_waveguide_frequency() {
    let l = 0.34;
    let c = 343.0;
    let plain = acoustic_waveguide(l, 0.012, 1.2, c, 8, false, 1, &[]).expect("plain");
    let vented = acoustic_waveguide(
        l,
        0.012,
        1.2,
        c,
        8,
        false,
        1,
        &[AcousticTap {
            station: 0.24,
            neck_length: 0.003,
            neck_radius: 0.003,
        }],
    )
    .expect("vented");
    let dt = 1.0 / 8_000.0;
    let ring = |sys: &PortHamiltonian| {
        let mut x = vec![0.0; sys.state_dim()];
        let mut p = Vec::new();
        for i in 0..640 {
            let u = if i < 16 {
                2.0e-5 * (core::f64::consts::PI * i as f64 / 16.0).sin()
            } else {
                0.0
            };
            let rec = step(sys, &x, &[u], dt).expect("step");
            p.push(rec.y[0]);
            x = rec.x;
        }
        dominant_zero_period(&p, dt)
    };
    let t0 = ring(&plain);
    let t1 = ring(&vented);
    assert!(
        t1 < t0 * 0.98,
        "an open side branch must shorten the period ({t1} vs {t0})"
    );
}

#[test]
fn equal_radius_chain_matches_a_uniform_waveguide() {
    let l = 0.34;
    let c = 343.0;
    let one = acoustic_waveguide(l, 0.012, 1.2, c, 8, false, 1, &[]).expect("one");
    let two = acoustic_chain(
        &[
            AcousticSection {
                length: 0.17,
                radius: 0.012,
                outlet_radius: 0.012,
                cells: 4,
            },
            AcousticSection {
                length: 0.17,
                radius: 0.012,
                outlet_radius: 0.012,
                cells: 4,
            },
        ],
        1.2,
        c,
        false,
        1,
        &[],
        None,
    )
    .expect("two");
    assert_eq!(one.state_dim(), two.state_dim());
    let dt = 1.0 / 8_000.0;
    let mut x1 = vec![0.0; one.state_dim()];
    let mut x2 = vec![0.0; two.state_dim()];
    for i in 0..32 {
        let u = if i < 8 {
            2.0e-5 * (core::f64::consts::PI * i as f64 / 8.0).sin()
        } else {
            0.0
        };
        let a = step(&one, &x1, &[u], dt).expect("a");
        let b = step(&two, &x2, &[u], dt).expect("b");
        assert!(
            (a.y[0] - b.y[0]).abs() <= 1.0e-12 * (1.0 + a.y[0].abs()),
            "equal-radius split must be the same LC line ({} vs {})",
            a.y[0],
            b.y[0]
        );
        x1 = a.x;
        x2 = b.x;
    }
}

#[test]
fn a_constriction_shifts_the_chain_period() {
    let l = 0.34;
    let c = 343.0;
    let plain = acoustic_waveguide(l, 0.012, 1.2, c, 8, false, 1, &[]).expect("plain");
    let stepped = acoustic_chain(
        &[
            AcousticSection {
                length: 0.17,
                radius: 0.012,
                outlet_radius: 0.012,
                cells: 4,
            },
            AcousticSection {
                length: 0.17,
                radius: 0.006,
                outlet_radius: 0.006,
                cells: 4,
            },
        ],
        1.2,
        c,
        false,
        1,
        &[],
        None,
    )
    .expect("stepped");
    let dt = 1.0 / 8_000.0;
    let ring = |sys: &PortHamiltonian| {
        let mut x = vec![0.0; sys.state_dim()];
        let mut p = Vec::new();
        for i in 0..640 {
            let u = if i < 16 {
                2.0e-5 * (core::f64::consts::PI * i as f64 / 16.0).sin()
            } else {
                0.0
            };
            let rec = step(sys, &x, &[u], dt).expect("step");
            p.push(rec.y[0]);
            x = rec.x;
        }
        dominant_zero_period(&p, dt)
    };
    let t0 = ring(&plain);
    let t1 = ring(&stepped);
    assert!(
        (t1 - t0).abs() > 0.02 * t0,
        "an area jump must move the ringing period ({t1} vs {t0})"
    );
    assert!(acoustic_chain(&[], 1.2, c, false, 1, &[], None).is_err());
    assert!(
        acoustic_chain(
            &[AcousticSection {
                length: 0.1,
                radius: 0.01,
                outlet_radius: 0.01,
                cells: 1,
            }],
            1.2,
            c,
            false,
            1,
            &[],
            None,
        )
        .is_err()
    );
}

#[test]
fn viscothermal_pin_damps_the_chain() {
    let l = 0.34;
    let c = 343.0;
    let sections = [AcousticSection {
        length: l,
        radius: 0.006,
        outlet_radius: 0.006,
        cells: 8,
    }];
    let lossless = acoustic_chain(&sections, 1.2, c, false, 1, &[], None).expect("lossless");
    let pin = ViscothermalPin {
        dynamic_viscosity: 1.8e-5,
        gamma: 1.4,
        prandtl: 0.71,
        foster_branches: 0,
    };
    let lossy = acoustic_chain(&sections, 1.2, c, false, 1, &[], Some(&pin)).expect("lossy");
    let zero = ViscothermalPin {
        dynamic_viscosity: 0.0,
        gamma: 1.4,
        prandtl: 0.71,
        foster_branches: 0,
    };
    let muted = acoustic_chain(&sections, 1.2, c, false, 1, &[], Some(&zero)).expect("zero mu");
    let dt = 1.0 / 8_000.0;
    let ring = |sys: &PortHamiltonian| {
        let mut x = vec![0.0; sys.state_dim()];
        let mut p = Vec::new();
        for i in 0..640 {
            let u = if i < 16 {
                2.0e-5 * (core::f64::consts::PI * i as f64 / 16.0).sin()
            } else {
                0.0
            };
            let rec = step(sys, &x, &[u], dt).expect("step");
            p.push(rec.y[0]);
            x = rec.x;
        }
        p
    };
    let p0 = ring(&lossless);
    let p1 = ring(&lossy);
    let p_z = ring(&muted);
    for i in 0..32 {
        assert!(
            (p0[i] - p_z[i]).abs() <= 1.0e-12 * (1.0 + p0[i].abs()),
            "zero viscosity must be the lossless mutation"
        );
    }
    let tail = |p: &[f64]| {
        let t = &p[p.len() / 2..];
        (t.iter().map(|x| x * x).sum::<f64>() / t.len() as f64).sqrt()
    };
    let e0 = tail(&p0);
    let e1 = tail(&p1);
    assert!(
        e1 < e0 * 0.85,
        "wide-tube pin must damp the tail ({e1} vs {e0})"
    );
    assert!(
        acoustic_chain(
            &sections,
            1.2,
            c,
            false,
            1,
            &[],
            Some(&ViscothermalPin {
                dynamic_viscosity: -1.0,
                gamma: 1.4,
                prandtl: 0.71,
                foster_branches: 0,
            }),
        )
        .is_err()
    );
}

#[test]
fn narrow_tube_pin_uses_poiseuille_and_still_damps() {
    let l = 0.05;
    let c = 343.0;
    let sections = [AcousticSection {
        length: l,
        radius: 2.0e-4,
        outlet_radius: 2.0e-4,
        cells: 6,
    }];
    let pin = ViscothermalPin {
        dynamic_viscosity: 1.8e-5,
        gamma: 1.4,
        prandtl: 0.71,
        foster_branches: 0,
    };
    let omega = core::f64::consts::PI * c / (2.0 * l);
    let rv = 2.0e-4 * (omega * 1.2 / 1.8e-5).sqrt();
    assert!(
        rv < 10.0,
        "fixture must sit in the Poiseuille branch (rv={rv})"
    );
    let lossless = acoustic_chain(&sections, 1.2, c, false, 1, &[], None).expect("lossless");
    let lossy = acoustic_chain(&sections, 1.2, c, false, 1, &[], Some(&pin)).expect("poiseuille");
    let dt = 1.0 / 16_000.0;
    let ring = |sys: &PortHamiltonian| {
        let mut x = vec![0.0; sys.state_dim()];
        let mut p = Vec::new();
        for i in 0..400 {
            let u = if i < 8 {
                1.0e-8 * (core::f64::consts::PI * i as f64 / 8.0).sin()
            } else {
                0.0
            };
            let rec = step(sys, &x, &[u], dt).expect("step");
            p.push(rec.y[0]);
            x = rec.x;
        }
        p
    };
    let p0 = ring(&lossless);
    let p1 = ring(&lossy);
    let tail = |p: &[f64]| {
        let t = &p[p.len() / 2..];
        (t.iter().map(|x| x * x).sum::<f64>() / t.len() as f64).sqrt()
    };
    assert!(
        tail(&p1) < tail(&p0) * 0.9,
        "Poiseuille pin must damp a capillary ({} vs {})",
        tail(&p1),
        tail(&p0)
    );
}

#[test]
fn foster_sqrt_omega_matches_the_wall_law_and_adds_states() {
    let gain = 2.0e-3;
    let terms = foster_sqrt_omega_terms(gain, 200.0, 12_800.0, 3).expect("terms");
    assert_eq!(terms.len(), 3);
    assert!(terms.iter().all(|&(g, w)| g > 0.0 && w > 0.0));
    let mid = 1_600.0;
    let re: f64 = terms
        .iter()
        .map(|&(g, w)| g * mid * mid / (mid * mid + w * w))
        .sum();
    let target = gain * mid.sqrt();
    assert!(
        (re - target).abs() < 0.25 * target,
        "Foster Re Z must sit on K√ω ({re} vs {target})"
    );
    let l = 0.34;
    let c = 343.0;
    let sections = [AcousticSection {
        length: l,
        radius: 0.006,
        outlet_radius: 0.006,
        cells: 6,
    }];
    let lumped = ViscothermalPin {
        dynamic_viscosity: 1.8e-5,
        gamma: 1.4,
        prandtl: 0.71,
        foster_branches: 0,
    };
    let foster = ViscothermalPin {
        foster_branches: 3,
        ..lumped
    };
    let lossless = acoustic_chain(&sections, 1.2, c, false, 1, &[], None).expect("lossless");
    let pin = acoustic_chain(&sections, 1.2, c, false, 1, &[], Some(&lumped)).expect("pin");
    let spec = acoustic_chain(&sections, 1.2, c, false, 1, &[], Some(&foster)).expect("foster");
    acoustic_chain(&sections, 1.2, c, true, 1, &[], Some(&foster))
        .expect("open Foster must stay PSD after adding Re Z_rad");
    assert_eq!(spec.state_dim(), pin.state_dim() + 6 * 3 * 2);
    let dt = 1.0 / 8_000.0;
    let ring = |sys: &PortHamiltonian| {
        let mut x = vec![0.0; sys.state_dim()];
        let mut p = Vec::new();
        for i in 0..480 {
            let u = if i < 16 {
                2.0e-5 * (core::f64::consts::PI * i as f64 / 16.0).sin()
            } else {
                0.0
            };
            let rec = step(sys, &x, &[u], dt).expect("step");
            p.push(rec.y[0]);
            x = rec.x;
        }
        p
    };
    let p0 = ring(&pin);
    let p1 = ring(&spec);
    let p_l = ring(&lossless);
    let err: f64 = p0.iter().zip(&p1).map(|(a, b)| (a - b).abs()).sum();
    assert!(
        err > 1.0e-8,
        "Bessel Foster must not reprint the lumped pin"
    );
    let tail = |p: &[f64]| {
        let t = &p[p.len() / 2..];
        (t.iter().map(|x| x * x).sum::<f64>() / t.len() as f64).sqrt()
    };
    assert!(
        tail(&p1) < tail(&p_l) * 0.85,
        "Bessel Foster must damp versus lossless ({} vs {})",
        tail(&p1),
        tail(&p_l)
    );
}

#[test]
fn linear_taper_is_not_the_inlet_cylinder() {
    let c = 343.0;
    let cyl = slice_linear_taper(0.006, 0.006, 0.34, 8).expect("cyl");
    let flare = slice_linear_taper(0.006, 0.018, 0.34, 8).expect("flare");
    assert_eq!(cyl.len(), 1);
    assert_eq!(flare.len(), 1);
    assert_eq!(flare[0].cells, 8);
    assert!(flare[0].outlet_radius > flare[0].radius);
    let a = acoustic_chain(&cyl, 1.2, c, false, 1, &[], None).expect("cyl chain");
    let b = acoustic_chain(&flare, 1.2, c, false, 1, &[], None).expect("flare chain");
    // Cylinder LC is 2 states/cell; the cone adds one near-field
    // shunt inertance per cell. Equal dimensions would mean the
    // flare reprinted the inlet ladder.
    assert_eq!(a.state_dim(), 16);
    assert_eq!(b.state_dim(), 24);
    let dt = 1.0 / 8_000.0;
    let ring = |sys: &PortHamiltonian| {
        let mut x = vec![0.0; sys.state_dim()];
        let mut p = Vec::new();
        for i in 0..640 {
            let u = if i < 16 {
                2.0e-5 * (core::f64::consts::PI * i as f64 / 16.0).sin()
            } else {
                0.0
            };
            let rec = step(sys, &x, &[u], dt).expect("step");
            p.push(rec.y[0]);
            x = rec.x;
        }
        p
    };
    let pa = ring(&a);
    let pb = ring(&b);
    let err: f64 = pa.iter().zip(&pb).map(|(x, y)| (x - y).abs()).sum();
    let peak_a = pa.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
    let peak_b = pb.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
    assert!(
        err > 1.0e-8 * (1.0 + peak_a) * pa.len() as f64
            || (peak_a - peak_b).abs() > 1.0e-8 * (1.0 + peak_a),
        "a linear flare must not reprint the inlet cylinder (err={err}, peaks {peak_a} vs {peak_b})"
    );
    assert!(slice_linear_taper(0.0, 0.01, 0.1, 4).is_err());
}

#[test]
fn spherical_cone_is_not_the_frustum_ladder() {
    let c = 343.0;
    let sph = spherical_cone(0.006, 0.018, 0.34, 1.2, c, 4, false, 1, &[], None).expect("ψ");
    assert_eq!(sph.state_dim(), 12);
    let halves = [
        AcousticSection {
            length: 0.17,
            radius: 0.006,
            outlet_radius: 0.012,
            cells: 2,
        },
        AcousticSection {
            length: 0.17,
            radius: 0.012,
            outlet_radius: 0.018,
            cells: 2,
        },
    ];
    let web = acoustic_chain(&halves, 1.2, c, false, 1, &[], None).expect("frustum");
    let dt = 1.0 / 8_000.0;
    let ring = |sys: &PortHamiltonian| {
        let mut x = vec![0.0; sys.state_dim()];
        let mut p = Vec::new();
        for i in 0..480 {
            let u = if i < 16 {
                2.0e-5 * (core::f64::consts::PI * i as f64 / 16.0).sin()
            } else {
                0.0
            };
            let rec = step(sys, &x, &[u], dt).expect("step");
            p.push(rec.y[0]);
            x = rec.x;
        }
        p
    };
    let pa = ring(&sph);
    let pb = ring(&web);
    let err: f64 = pa.iter().zip(&pb).map(|(x, y)| (x - y).abs()).sum();
    assert!(
        err > 1.0e-6,
        "ψ = x p must not reprint the frustum LC ladder"
    );
    assert!(spherical_cone(0.01, 0.01, 0.2, 1.2, c, 4, false, 1, &[], None).is_err());
}

#[test]
fn zwikker_kosten_f_hits_both_regime_limits() {
    let wide = zwikker_kosten_f(20.0).expect("wide");
    // Wide-tube: F ≈ (1−i) √2 / r_v, Re F ≈ Im(−F) ≈ √2 / 20.
    let eps = core::f64::consts::SQRT_2 / 20.0;
    assert!(
        (wide.re - eps).abs() < 0.4 * eps,
        "wide-tube Re F ({}) should sit near √2/r_v ({eps})",
        wide.re
    );
    let narrow = zwikker_kosten_f(1.0).expect("narrow");
    // r_v → 0 ⇒ F → 1; r_v = 1 is already close.
    assert!(
        (narrow.re - 1.0).abs() < 0.25,
        "narrow-tube Re F should approach 1, got {}",
        narrow.re
    );
    // 1−F for the series law must give Re Z > 0.
    let den = C64::new(1.0 - wide.re, -wide.im);
    assert!(den.im / den.norm_sq() > 0.0);
}

fn dominant_zero_period(x: &[f64], dt: f64) -> f64 {
    let mut prev = x[x.len() / 4];
    let mut times = Vec::new();
    for (i, &s) in x.iter().enumerate().skip(x.len() / 4 + 1) {
        if prev > 0.0 && s <= 0.0 {
            times.push(i as f64 * dt);
        }
        prev = s;
    }
    assert!(times.len() >= 3, "need crossings, got {}", times.len());
    (times[times.len() - 1] - times[0]) / (times.len() - 1) as f64
}
