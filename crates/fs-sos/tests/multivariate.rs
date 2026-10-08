//! G0/G2 battery for the multivariate proof-carrying layer: certified global
//! enclosures against known minima, the Motzkin refusal, Putinar-constrained
//! bounds, tamper falsifiers, and SOS region-of-attraction soundness checked
//! against independent trajectory simulation.

use fs_ivl::Interval;
use fs_sos::{
    DecisionValue, GlobalError, GlobalOptions, Identity, MPoly, RoaOptions, SosProgram,
    VerifyError, certify_roa, minimize, monomials_in_degree_range, verify,
};

fn x(n: usize, k: usize) -> MPoly {
    MPoly::var(n, k)
}

/// Dense grid minimum (independent oracle for low-dimensional fixtures).
fn grid_min_2d(p: &MPoly, lo: f64, hi: f64, steps: usize) -> f64 {
    let mut best = f64::INFINITY;
    for i in 0..=steps {
        for j in 0..=steps {
            let a = lo + (hi - lo) * i as f64 / steps as f64;
            let b = lo + (hi - lo) * j as f64 / steps as f64;
            best = best.min(p.eval(&[a, b]));
        }
    }
    best
}

#[test]
fn sos_g0_quartic_global_minimum_is_enclosed() {
    // p = x⁴ + y⁴ − 3x² + x + y² − 2xy: nonconvex, unique global minimizer.
    let (a, b) = (x(2, 0), x(2, 1));
    let p = a
        .pow(4)
        .add(&b.pow(4))
        .sub(&a.pow(2).scale(3.0))
        .add(&a)
        .add(&b.pow(2))
        .sub(&a.mul(&b).scale(2.0));
    let r = minimize(&p, &[], &GlobalOptions::default()).expect("certified");
    let grid = grid_min_2d(&p, -2.5, 2.5, 500);
    let upper = r.upper.expect("upper bound from moments");
    assert!(r.lower <= upper, "{r:?}");
    // The proved lower bound is below every sampled value...
    assert!(r.lower <= grid + 1e-12, "lower {} > grid {}", r.lower, grid);
    // ...and the enclosure is tight (exact relaxation, small back-off).
    assert!(upper - r.lower < 1e-5, "gap {}", upper - r.lower);
    assert!(upper <= grid + 1e-9, "upper {upper} grid {grid}");
}

#[test]
fn sos_g2_six_hump_camel_matches_literature_value() {
    // f = 4x² − 2.1x⁴ + x⁶/3 + xy − 4y² + 4y⁴; global min −1.031628453489877
    // at (±0.0898, ∓0.7126) — two minimizers (moment matrix of rank 2).
    let (a, b) = (x(2, 0), x(2, 1));
    let p = a
        .pow(2)
        .scale(4.0)
        .sub(&a.pow(4).scale(2.1))
        .add(&a.pow(6).scale(1.0 / 3.0))
        .add(&a.mul(&b))
        .sub(&b.pow(2).scale(4.0))
        .add(&b.pow(4).scale(4.0));
    let r = minimize(&p, &[], &GlobalOptions::default()).expect("certified");
    let fstar = -1.031_628_453_489_877;
    assert!(r.lower <= fstar, "lower {} above the true minimum", r.lower);
    assert!(fstar - r.lower < 1e-5, "lower bound loose: {}", r.lower);
    let upper = r.upper.expect("polished minimizer");
    assert!(
        upper >= fstar - 1e-12 && upper - fstar < 1e-9,
        "upper {upper}"
    );
    let xm = r.minimizer.unwrap();
    assert!((xm[0].abs() - 0.0898).abs() < 1e-3 && (xm[1].abs() - 0.7126).abs() < 1e-3);
}

#[test]
fn sos_falsifier_motzkin_is_refused_not_misclaimed() {
    // Motzkin: x⁴y² + x²y⁴ − 3x²y² + 1 ≥ 0 (min 0) but M − γ is not SOS for
    // ANY γ at the base order. The pipeline must refuse, never claim a bound.
    let (a, b) = (x(2, 0), x(2, 1));
    let m = a
        .pow(4)
        .mul(&b.pow(2))
        .add(&a.pow(2).mul(&b.pow(4)))
        .sub(&a.pow(2).mul(&b.pow(2)).scale(3.0))
        .add(&MPoly::constant(2, 1.0));
    match minimize(&m, &[], &GlobalOptions::default()) {
        Err(GlobalError::RelaxationFailed { .. } | GlobalError::CertificationFailed { .. }) => {}
        Ok(r) => {
            // Any accepted answer must still be SOUND (min is 0).
            assert!(r.lower <= 0.0, "unsound Motzkin bound {}", r.lower);
        }
        Err(e) => panic!("unexpected error kind {e:?}"),
    }
}

#[test]
fn sos_g0_putinar_disk_bound() {
    // min x + y on the unit disk = −√2.
    let (a, b) = (x(2, 0), x(2, 1));
    let p = a.add(&b);
    let g = MPoly::constant(2, 1.0).sub(&a.pow(2)).sub(&b.pow(2));
    let r = minimize(&p, std::slice::from_ref(&g), &GlobalOptions::default()).expect("certified");
    let truth = -std::f64::consts::SQRT_2;
    assert!(r.lower <= truth, "lower {} > −√2", r.lower);
    assert!(truth - r.lower < 1e-5, "loose: {}", r.lower);
    if let (Some(u), Some(xm)) = (r.upper, r.minimizer.as_ref()) {
        assert!(u >= truth - 1e-12);
        let pt: Vec<Interval> = xm.iter().map(|&v| Interval::point(v)).collect();
        assert!(
            g.eval_interval(&pt).lo() >= 0.0,
            "extracted point infeasible"
        );
    }
}

#[test]
fn sos_falsifier_inflated_bound_and_tampered_gram_are_rejected() {
    // p = (x − 1)² + (y + 2)² + 3 has minimum 3. Try to prove γ = 3.1.
    let (a, b) = (x(2, 0), x(2, 1));
    let one = MPoly::constant(2, 1.0);
    let p = a
        .sub(&one)
        .pow(2)
        .add(&b.add(&one.scale(2.0)).pow(2))
        .add(&one.scale(3.0));
    let mut prog = SosProgram::new(2);
    let gamma = prog.scalar("gamma");
    let s0 = prog.sos("s0", monomials_in_degree_range(2, 0, 1)).unwrap();
    prog.add_identity(
        Identity::new(&p)
            .term(one.scale(-1.0), gamma)
            .term(one.scale(-1.0), s0),
    )
    .unwrap();
    // A centred solve at the false level cannot verify.
    let sol = prog
        .solve_centered(&[(gamma, 3.1)], &fs_sos::SdpSettings::default())
        .unwrap();
    assert!(
        verify(&prog, &sol.values).is_err(),
        "false bound 3.1 certified"
    );
    // At a true level it verifies...
    let good = prog
        .solve_centered(&[(gamma, 2.9)], &fs_sos::SdpSettings::default())
        .unwrap();
    verify(&prog, &good.values).expect("γ = 2.9 must certify");
    // ...and a hand-forged Gram claiming γ = 3.5 is refused.
    let mut forged = good.values.clone();
    forged[gamma.index()] = DecisionValue::Scalar(3.5);
    assert!(matches!(
        verify(&prog, &forged),
        Err(VerifyError::NotPositiveDefinite { .. })
    ));
}

#[test]
fn sos_g5_global_bound_replays_bit_identically() {
    let (a, b) = (x(2, 0), x(2, 1));
    let p = a
        .pow(4)
        .add(&b.pow(4))
        .sub(&a.mul(&b).scale(2.0))
        .add(&a.scale(0.25));
    let r1 = minimize(&p, &[], &GlobalOptions::default()).unwrap();
    let r2 = minimize(&p, &[], &GlobalOptions::default()).unwrap();
    assert_eq!(r1.lower.to_bits(), r2.lower.to_bits());
    assert_eq!(r1.upper.map(f64::to_bits), r2.upper.map(f64::to_bits));
    assert_eq!(r1.certificate, r2.certificate);
}

#[test]
fn roa_g1_cubic_scalar_system_is_sound_and_tight() {
    // ẋ = −x + x³: true region of attraction (−1, 1). With P = 1/2,
    // V = x²/2, so any proved level must satisfy 2c < 1 → c < 0.5.
    let xv = x(1, 0);
    let f = vec![xv.scale(-1.0).add(&xv.pow(3))];
    let cert = certify_roa(&f, &RoaOptions::default()).expect("certified");
    assert!((cert.lyapunov[0] - 0.5).abs() < 1e-14);
    assert!(cert.level < 0.5, "UNSOUND level {}", cert.level);
    assert!(cert.level > 0.40, "loose level {}", cert.level);
    assert!(cert.semi_axes[0] < 1.0);
}

/// RK4 integrate and report whether the trajectory converged to 0.
fn converges(f: &[MPoly], mut s: Vec<f64>, dt: f64, steps: usize) -> bool {
    let eval = |s: &[f64]| -> Vec<f64> { f.iter().map(|fi| fi.eval(s)).collect() };
    for _ in 0..steps {
        let k1 = eval(&s);
        let s2: Vec<f64> = s.iter().zip(&k1).map(|(a, k)| a + 0.5 * dt * k).collect();
        let k2 = eval(&s2);
        let s3: Vec<f64> = s.iter().zip(&k2).map(|(a, k)| a + 0.5 * dt * k).collect();
        let k3 = eval(&s3);
        let s4: Vec<f64> = s.iter().zip(&k3).map(|(a, k)| a + dt * k).collect();
        let k4 = eval(&s4);
        for i in 0..s.len() {
            s[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        if s.iter().any(|v| !v.is_finite() || v.abs() > 1e6) {
            return false;
        }
    }
    s.iter().map(|v| v * v).sum::<f64>().sqrt() < 1e-3
}

#[test]
fn roa_g2_reversed_van_der_pol_boundary_trajectories_converge() {
    // Time-reversed Van der Pol: ẋ₁ = −x₂, ẋ₂ = x₁ + (x₁² − 1)x₂. Its ROA is
    // bounded by the (unstable) limit cycle — the classic SOS benchmark.
    let (a, b) = (x(2, 0), x(2, 1));
    let f = vec![b.scale(-1.0), a.add(&a.pow(2).mul(&b)).sub(&b)];
    let cert = certify_roa(&f, &RoaOptions::default()).expect("certified");
    assert!(cert.level > 0.0 && cert.volume > 0.5, "{cert:?}");
    // Independent check: start on the boundary of the certified ellipsoid
    // in 64 directions; every trajectory must converge to the origin, and
    // V̇ < 0 there.
    let p = &cert.lyapunov;
    let v = |s: &[f64]| p[0] * s[0] * s[0] + 2.0 * p[1] * s[0] * s[1] + p[3] * s[1] * s[1];
    for k in 0..64 {
        let th = 2.0 * std::f64::consts::PI * f64::from(k) / 64.0;
        let dir = [th.cos(), th.sin()];
        let scale = (cert.level / v(&dir)).sqrt() * (1.0 - 1e-9);
        let s0 = vec![dir[0] * scale, dir[1] * scale];
        let fx: Vec<f64> = f.iter().map(|fi| fi.eval(&s0)).collect();
        let vdot =
            2.0 * ((p[0] * s0[0] + p[1] * s0[1]) * fx[0] + (p[1] * s0[0] + p[3] * s0[1]) * fx[1]);
        assert!(
            vdot < 0.0,
            "V̇ = {vdot} ≥ 0 on the certified boundary at θ = {th}"
        );
        assert!(
            converges(&f, s0, 1e-2, 4000),
            "trajectory from θ = {th} escaped"
        );
    }
    // Falsifier: just outside the limit cycle (radius ~2 along x₁) the
    // trajectory escapes, and the certificate never reaches it.
    let semi_max = cert.semi_axes.iter().fold(0.0f64, |m, s| m.max(*s));
    assert!(
        semi_max < 2.1,
        "certified set reaches past the limit cycle: {semi_max}"
    );
    assert!(!converges(&f, vec![2.5, 0.0], 1e-2, 4000));
}

#[test]
fn roa_falsifier_unstable_and_non_equilibrium_refused() {
    let (a, b) = (x(2, 0), x(2, 1));
    // Unstable linearization (saddle): P from the Lyapunov equation is not PD.
    let f = vec![a.clone(), b.scale(-1.0)];
    assert!(certify_roa(&f, &RoaOptions::default()).is_err());
    // Origin not an equilibrium.
    let g = vec![a.scale(-1.0).add(&MPoly::constant(2, 0.1)), b.scale(-1.0)];
    assert!(matches!(
        certify_roa(&g, &RoaOptions::default()),
        Err(fs_sos::RoaError::NotEquilibrium { component: 0 })
    ));
}

#[test]
fn roa_linear_system_reaches_the_cap() {
    // Globally stable linear dynamics: every level is provable.
    let (a, b) = (x(2, 0), x(2, 1));
    let f = vec![b.clone(), a.scale(-2.0).sub(&b.scale(0.5))];
    let opts = RoaOptions {
        level_cap: 1e3,
        ..RoaOptions::default()
    };
    let cert = certify_roa(&f, &opts).expect("certified");
    assert_eq!(cert.level.to_bits(), 1e3f64.to_bits());
    assert!(cert.unproved_level.is_none());
}
