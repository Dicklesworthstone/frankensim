use super::*;

/// Parallel-RLC impedance `2 zeta w0 g s / (s^2 + 2 zeta w0 s + w0^2)`: a
/// positive-real scalar, returned as its stored pole and residue.
fn rlc(w0: f64, zeta: f64, g: f64) -> (C64, C64) {
    let pole = C64::new(-zeta * w0, w0 * det::sqrt(1.0 - zeta * zeta));
    let residue = (pole.scale(2.0 * zeta * w0 * g)) * C64::new(0.0, 2.0 * pole.im).recip();
    (pole, residue)
}

/// `Z = SUM_k f_k(s) u_k u_k^T + d0 I`: passive by construction, reciprocal.
fn passive_model(d0: f64) -> MimoModel {
    let n = 3;
    let shapes = [
        [1.0, 0.4, -0.2],
        [0.3, 1.0, 0.5],
        [-0.6, 0.2, 1.0],
        [0.5, -0.5, 0.7],
    ];
    let modes = [
        (300.0, 0.05, 2.0),
        (900.0, 0.03, 1.0),
        (2500.0, 0.08, 3.0),
        (6000.0, 0.1, 1.5),
    ];
    let mut poles = Vec::new();
    let mut residues = Vec::new();
    for (u, &(w0, zeta, g)) in shapes.iter().zip(&modes) {
        let (pole, r) = rlc(w0, zeta, g);
        poles.push(PoleTerm::Pair {
            pole,
            residue: C64::ZERO,
        });
        residues.push((0..n * n).map(|i| r.scale(u[i / n] * u[i % n])).collect());
    }
    let mut d = vec![0.0; n * n];
    for i in 0..n {
        d[i * n + i] = d0;
    }
    MimoModel {
        ports: n,
        poles,
        residues,
        d,
        e: vec![0.0; n * n],
        reciprocal: true,
    }
}

fn band() -> Vec<f64> {
    (0..160)
        .map(|i| 100.0 * det::exp(det::ln(100.0) * f64::from(i) / 159.0))
        .collect()
}

#[test]
fn mimo_fit_recovers_a_shared_pole_reciprocal_model_g1() {
    let truth = passive_model(5.0);
    let omega = band();
    let data: Vec<Vec<C64>> = omega.iter().map(|&w| truth.eval_iw(w)).collect();
    let (fit, report) = mimo_vector_fit(&omega, &data, 3, &FitOptions::new(8), true).unwrap();
    assert!(fit.is_stable() && fit.reciprocal);
    assert!(
        report.max_relative_error < 1e-8,
        "known-model recovery: {report:?}"
    );
    assert!(report.data_reciprocity_defect < 1e-14);
    // Every true pole is recovered.
    for t in &truth.poles {
        let PoleTerm::Pair { pole, .. } = *t else {
            unreachable!()
        };
        let best = fit
            .poles
            .iter()
            .map(|f| match *f {
                PoleTerm::Pair { pole: q, .. } => (q - pole).abs(),
                PoleTerm::Real { pole: q, .. } => (C64::from_re(q) - pole).abs(),
            })
            .fold(f64::INFINITY, f64::min);
        assert!(best < 1e-6 * pole.abs(), "pole {pole:?} missing: {best:e}");
    }
    // Planted asymmetry is fitted by its symmetric part and reported.
    let skewed: Vec<Vec<C64>> = data
        .iter()
        .map(|z| {
            let mut z = z.clone();
            z[1] = z[1].scale(1.01);
            z
        })
        .collect();
    let (_, skew_report) = mimo_vector_fit(&omega, &skewed, 3, &FitOptions::new(8), true).unwrap();
    // A symmetric model's misfit is exactly the data's antisymmetric part,
    // half of ||Z - Z^T||, once the symmetric part is fitted to roundoff.
    assert!(skew_report.data_reciprocity_defect > 1e-3);
    assert!(
        (skew_report.max_relative_error - 0.5 * skew_report.data_reciprocity_defect).abs()
            <= 1e-3 * skew_report.data_reciprocity_defect,
        "{skew_report:?}"
    );
}

#[test]
fn mimo_state_space_matches_pole_residue_eval_g0() {
    let mut model = passive_model(2.0);
    model.poles.push(PoleTerm::Real {
        pole: -700.0,
        residue: 0.0,
    });
    model.residues.push(
        (0..9)
            .map(|i| C64::from_re(0.1 * (i as f64 - 4.0)))
            .collect(),
    );
    model.e[4] = 1e-4;
    let ss = model.state_space();
    assert_eq!(ss.n, model.order() * 3);
    for w in [50.0, 333.0, 1234.0, 7000.0] {
        let a = model.eval_iw(w);
        let b = ss.eval(C64::new(0.0, w)).unwrap();
        let gap: Vec<C64> = a.iter().zip(&b).map(|(&x, &y)| x - y).collect();
        assert!(
            frobenius(&gap) <= 1e-12 * frobenius(&a),
            "realization parity at {w}"
        );
    }
}

#[test]
fn mimo_passivity_certifies_passive_and_repairs_a_planted_violation_g1() {
    let passive = passive_model(0.5);
    let report = check_mimo_passivity(&passive, (100.0, 10_000.0)).unwrap();
    assert!(report.passive, "{report:?}");
    assert_eq!(report.class, CertificateClass::HamiltonianExact);
    // A planted off-diagonal residue makes one eigenvalue of Z + Z^H negative
    // near the second mode.
    let mut active = passive.clone();
    for slot in [1, 3] {
        active.residues[1][slot] = active.residues[1][slot] + C64::new(-40.0, 15.0);
    }
    let bad = check_mimo_passivity(&active, (100.0, 10_000.0)).unwrap();
    assert!(!bad.passive && bad.worst.0 < 0.0, "{bad:?}");
    assert!(
        !bad.crossings.is_empty(),
        "the Hamiltonian must see the dip"
    );
    // Each crossing sits where the grid's smallest eigenvalue changes sign.
    for &w in &bad.crossings {
        let below = hermitian_part(&active.eval_iw(w * (1.0 - 1e-4)), 3).0[0];
        let above = hermitian_part(&active.eval_iw(w * (1.0 + 1e-4)), 3).0[0];
        assert!(
            below * above <= 0.0 || below.abs().min(above.abs()) < 1e-6,
            "{w}"
        );
    }
    let (repaired, repair) = repair_mimo_passivity(&active, (100.0, 10_000.0)).unwrap();
    assert!(repair.certificate.passive);
    assert_eq!(repaired.poles, active.poles);
    assert_eq!(
        (repaired.d.clone(), repaired.e.clone()),
        (active.d.clone(), active.e.clone())
    );
    assert!(repair.relative_perturbation < 0.5, "{repair:?}");
    // Reciprocity survives repair.
    for r in &repaired.residues {
        for a in 0..3 {
            for b in 0..3 {
                assert_eq!(r[a * 3 + b], r[b * 3 + a]);
            }
        }
    }
}

#[test]
fn mimo_refusals_g0() {
    let omega = band();
    let truth = passive_model(1.0);
    let data: Vec<Vec<C64>> = omega.iter().map(|&w| truth.eval_iw(w)).collect();
    assert!(mimo_vector_fit(&omega, &data, 3, &FitOptions::new(0), true).is_err());
    assert!(mimo_vector_fit(&omega[..5], &data[..5], 3, &FitOptions::new(8), true).is_err());
    assert!(mimo_vector_fit(&omega, &data, 2, &FitOptions::new(8), true).is_err());
    let mut nan = data.clone();
    nan[7][2] = C64::new(f64::NAN, 0.0);
    assert!(mimo_vector_fit(&omega, &nan, 3, &FitOptions::new(8), true).is_err());
    let mut unstable = truth.clone();
    unstable.poles[0] = PoleTerm::Pair {
        pole: C64::new(10.0, 300.0),
        residue: C64::ZERO,
    };
    assert!(matches!(
        check_mimo_passivity(&unstable, (100.0, 1e4)),
        Err(PassivityError::Unstable)
    ));
    let mut negative = truth.clone();
    negative.d[0] = -3.0;
    assert!(matches!(
        check_mimo_passivity(&negative, (100.0, 1e4)),
        Err(PassivityError::NegativeAsymptote)
    ));
    // A zero D is the named weaker grid-only class, not a refusal.
    let mut lossless = truth.clone();
    lossless.d = vec![0.0; 9];
    assert_eq!(
        check_mimo_passivity(&lossless, (100.0, 1e4)).unwrap().class,
        CertificateClass::GridOnly
    );
}

#[test]
fn delayed_mimo_fit_recovers_known_delays_g1() {
    let truth = passive_model(2.0);
    let delays: Vec<f64> = (0..9)
        .map(|i| {
            let (a, b) = (i / 3, i % 3);
            if a == b { 0.0 } else { 1e-4 * (a + b) as f64 }
        })
        .collect();
    let delayed = DelayedMimoModel {
        delays: delays.clone(),
        model: truth,
    };
    let omega = band();
    let data: Vec<Vec<C64>> = omega.iter().map(|&w| delayed.eval_iw(w)).collect();
    let opts = FitOptions::new(8);
    let (fit, report) = delayed_mimo_vector_fit(&omega, &data, 3, &delays, &opts, true).unwrap();
    assert!(report.max_relative_error < 1e-8, "{report:?}");
    // The same data without delay extraction is far worse at the same order.
    let (_, plain) = mimo_vector_fit(&omega, &data, 3, &opts, true).unwrap();
    assert!(
        plain.max_relative_error > 100.0 * report.max_relative_error,
        "{plain:?}"
    );
    // Delayed (off-diagonal) entries never carry the improper term.
    for a in 0..3 {
        for b in 0..3 {
            if a != b {
                assert_eq!(fit.model.e[a * 3 + b].to_bits(), 0);
            }
        }
    }
    let passive = check_delayed_mimo_passivity(&fit, (100.0, 10_000.0)).unwrap();
    assert_eq!(passive.class, CertificateClass::GridOnly);
    let mut asymmetric = delays.clone();
    asymmetric[1] += 1e-6;
    assert!(delayed_mimo_vector_fit(&omega, &data, 3, &asymmetric, &opts, true).is_err());
    let mut negative = delays;
    negative[1] = -1e-4;
    negative[3] = -1e-4;
    assert!(delayed_mimo_vector_fit(&omega, &data, 3, &negative, &opts, true).is_err());
}

#[test]
fn mimo_zoh_step_response_is_exact_at_the_samples_g1() {
    // A held input makes ZOH exact: x(nT) = (exp(p n T) - 1)/p u.
    let mut model = passive_model(0.7);
    model.poles.push(PoleTerm::Real {
        pole: -450.0,
        residue: 0.0,
    });
    model.residues.push(
        (0..9)
            .map(|i| C64::from_re(0.2 + 0.05 * i as f64))
            .collect(),
    );
    let t_s = 1.0 / 48_000.0;
    let mut zoh = model.zoh(t_s).unwrap();
    let u = [1.0, 0.5, -0.3];
    for step in 0..2000usize {
        let y = zoh.step(&u);
        let t = step as f64 * t_s;
        for a in 0..3 {
            let mut exact = (0..3).map(|b| model.d[a * 3 + b] * u[b]).sum::<f64>();
            for (term, r) in model.poles.iter().zip(&model.residues) {
                let (pole, pair) = match *term {
                    PoleTerm::Real { pole, .. } => (C64::from_re(pole), false),
                    PoleTerm::Pair { pole, .. } => (pole, true),
                };
                let x = pole.scale(t);
                let (s, c) = det::sin_cos(x.im);
                let growth = (C64::new(c, s).scale(det::exp(x.re)) - C64::ONE) * pole.recip();
                let acc = (0..3).fold(C64::ZERO, |acc, b| acc + r[a * 3 + b] * growth.scale(u[b]));
                exact += if pair { 2.0 * acc.re } else { acc.re };
            }
            assert!(
                (y[a] - exact).abs() <= 1e-9 * (1.0 + exact.abs()),
                "step {step} port {a}: {} vs {exact}",
                y[a]
            );
        }
    }
    assert!(model.zoh(0.0).is_err());
    let mut unstable = model.clone();
    unstable.poles[0] = PoleTerm::Pair {
        pole: C64::new(1.0, 300.0),
        residue: C64::ZERO,
    };
    assert!(unstable.zoh(t_s).is_err());
}

#[test]
fn sampled_grid_reads_a_zero_dc_resistance_against_the_direct_term_floor() {
    // Z(s) = 1e8 + r/(s + 1e6) with Z(0) = -delta: a radiation-like load whose
    // DC resistance is zero up to an O(1e8) cancellation, positive above DC.
    let model = |delta: f64| MimoModel {
        ports: 1,
        poles: vec![PoleTerm::Real { pole: -1.0e6, residue: 0.0 }],
        residues: vec![vec![C64::from_re(-(1.0e8 + delta) * 1.0e6)]],
        d: vec![1.0e8],
        e: vec![0.0],
        reciprocal: true,
    };
    // Floor 1e-12 * 1e8 on the eigenvalue of Z + Z^H, which is -2 delta at DC.
    for (delta, passive) in [(0.0, true), (1.0e-5, true), (1.0e-3, false)] {
        let zoh = model(delta).zoh(1.0 / 96_000.0).unwrap();
        for report in [zoh.sampled_passivity(129), zoh.mean_sampled_passivity(129)] {
            assert_eq!(report.passive, passive, "delta {delta}: {report:?}");
            assert!((report.worst.0 + 2.0 * delta).abs() < 1.0e-6, "raw worst kept: {report:?}");
            if !passive {
                assert_eq!(report.violation_freqs, [0.0]);
            }
        }
    }
}

#[test]
fn passivity_certificate_sees_a_negative_dc_resistance_below_its_grid() {
    // Re Z(i w) = r / (w^2 + 1) + 5000 / (w^2 + 1e4): for r = -1 it is -0.5 at
    // DC and crosses zero near 1 rad/s, a decade below the grid's 10 rad/s
    // start. D = 0, so no Hamiltonian test runs: only the DC probe sees it.
    let model = |r: f64| MimoModel {
        ports: 1,
        poles: vec![
            PoleTerm::Real { pole: -1.0, residue: 0.0 },
            PoleTerm::Real { pole: -100.0, residue: 0.0 },
        ],
        residues: vec![vec![C64::from_re(r)], vec![C64::from_re(50.0)]],
        d: vec![0.0],
        e: vec![0.0],
        reciprocal: true,
    };
    let band = (100.0, 10_000.0);
    let fine = check_mimo_passivity(&model(-0.4), band).unwrap();
    assert!(fine.passive && fine.violation_freqs.is_empty(), "{fine:?}");
    let report = check_mimo_passivity(&model(-1.0), band).unwrap();
    assert_eq!(report.class, CertificateClass::GridOnly);
    assert!(!report.passive, "{report:?}");
    assert_eq!(report.violation_freqs, [0.0]);
    assert!((report.worst.0 + 1.0).abs() < 1e-12 && report.worst.1 == 0.0, "{report:?}");
    // The residue repair can act on that DC violation.
    let (repaired, repair) = repair_mimo_passivity(&model(-1.0), band).unwrap();
    assert!(repair.certificate.passive, "{repair:?}");
    assert!(hermitian_part(&repaired.eval_iw(0.0), 1).0[0] >= 0.0);
}
