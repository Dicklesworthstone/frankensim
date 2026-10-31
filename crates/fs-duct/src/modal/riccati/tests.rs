//! M2.1-W3b continuous multimodal Riccati (IRP pexfp.72). Gauntlet tiers
//! are named per test; each prints one JSON line per rung or check.

use super::*;
use crate::modal::{
    bessel_j0, bessel_j0_certified, bessel_j1, mm_input_impedance, modal_characteristic_impedances,
};
use crate::profile::{MAX_RICCATI_STEPS, profile_impedance};
use crate::{Duct, Segment};
use fs_material::gas::GasSpec;
use fs_mms::{LadderSide, OrderGate, RefinementLadder, fit_order};
use std::time::Instant;

fn air20() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}
fn rel(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs()
}
fn cents(f: f64, reference: f64) -> f64 {
    1200.0 * (f / reference).log2()
}
fn tau_hz(f: f64) -> f64 {
    core::f64::consts::TAU * f
}
fn z00(
    profile: &RadiusProfile,
    f_hz: f64,
    loss: LossModel,
    load: MouthLoad<'_>,
    modes: usize,
    steps: usize,
) -> C64 {
    mm_profile_impedance(
        profile,
        &air20(),
        tau_hz(f_hz),
        loss,
        load,
        ModeCount(modes),
        RiccatiSteps(steps),
    )
    .expect("continuous modal impedance")
    .plane_impedance
}
/// Steps that keep `kappa |h|` just inside the guard for `modes` at the
/// narrowest radius (a conservative bound: cutoff from `gamma/r_min` alone).
fn guarded_steps(length: f64, r_min: f64, modes: usize) -> usize {
    if modes == 1 {
        return 64;
    }
    let gamma = *j1_roots_unchecked(modes - 1).last().unwrap();
    ((gamma / r_min * length / (0.95 * MAX_EVANESCENT_STEP)).ceil() as usize).max(16)
}
fn exp_horn() -> RadiusProfile {
    // 6 mm to 60 mm over 0.5 m: r = r0 exp(a x / 2).
    RadiusProfile::exponential(0.006, 2.0 * det::ln(10.0) / 0.5, 0.5).unwrap()
}
fn complex_exp(z: C64) -> C64 {
    let (s, c) = det::sin_cos(z.im);
    C64::new(c, s).scale(det::exp(z.re))
}

/// `J_n(x) = (1/2 pi) int_0^{2 pi} cos(n t - x sin t) dt` by the trapezoid
/// rule, spectrally accurate on the full period (error ~ 2|J_{2M}(x)|).
fn quad_jn(order: u32, x: f64) -> f64 {
    const M: u32 = 512;
    let mut acc = 0.0;
    for j in 0..M {
        let t = core::f64::consts::TAU * f64::from(j) / f64::from(M);
        let (s, _) = det::sin_cos(t);
        acc += det::sin_cos(f64::from(order) * t - x * s).1;
    }
    acc / f64::from(M)
}
/// `d/dx` of the `J0` quadrature, `-(1/2 pi) int sin t sin(x sin t) dt`:
/// an independent route to `J0'`, so `J0' = -J1` is checked as an identity.
fn quad_j0_prime(x: f64) -> f64 {
    const M: u32 = 512;
    let mut acc = 0.0;
    for j in 0..M {
        let t = core::f64::consts::TAU * f64::from(j) / f64::from(M);
        let (s, _) = det::sin_cos(t);
        acc -= s * det::sin_cos(x * s).0;
    }
    acc / f64::from(M)
}

#[test]
fn mm_bessel_certified_to_gamma16() {
    // G0. The bead also names a J/Y Wronskian. It is omitted on purpose:
    // no Y owner exists here, and the Wronskian is blind to common phase
    // errors (fs-aeroac's own finding). The independent spectral
    // quadrature below bounds J0, J1 and J0' = -J1 pointwise instead.
    let start = Instant::now();
    let (mut e0, mut e1, mut ed, mut legacy) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    let mut grid: Vec<f64> = (0..=1024)
        .map(|i| BESSEL_CERTIFIED_X * f64::from(i) / 1024.0)
        .collect();
    grid.extend([17.999_999_999, 18.0, 18.000_000_001]);
    for &x in &grid {
        e0 = e0.max((bessel_j0_certified(x) - quad_jn(0, x)).abs());
        if x <= 26.0 {
            // Disclosure, not a gate: the staircase's bit-frozen legacy series.
            legacy = legacy.max((bessel_j0(x) - quad_jn(0, x)).abs());
        }
        e1 = e1.max((bessel_j1(x) - quad_jn(1, x)).abs());
        ed = ed.max((quad_j0_prime(x) + bessel_j1(x)).abs());
    }
    let cap = certified_mode_cap();
    let roots = j1_roots_unchecked(cap - 1);
    let worst_root = roots
        .iter()
        .map(|&g| bessel_j1(g).abs())
        .fold(0.0, f64::max);
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"mm-bessel-certified\",\"inputs\":\"x in [0,{BESSEL_CERTIFIED_X}] \
         ({} points incl. the 18 switch)\",\"observed\":{{\"j0_certified\":{e0:e},\"legacy_j0_to_26\":{legacy:e},\"j1\":{e1:e},\"j0_prime_plus_j1\":{ed:e},\
         \"j1_at_roots\":{worst_root:e},\"gamma_max\":{}}},\"tol\":1e-13,\"elapsed_us\":{}}}",
        grid.len(),
        roots.last().unwrap(),
        start.elapsed().as_micros()
    );
    assert!(e0 < 1e-13 && e1 < 1e-13 && ed < 1e-13 && worst_root < 1e-13);
    assert!(
        *roots.last().unwrap() > 51.0,
        "certified range reaches past gamma_16"
    );
}

#[test]
fn mm_mode_cap_refusal() {
    // G0: the cap is computed from the roots; cap admits and cap + 1 refuses.
    let cap = certified_mode_cap();
    let roots = j1_roots_unchecked(cap);
    assert!(roots[cap - 2] < BESSEL_CERTIFIED_X && roots[cap - 1] >= BESSEL_CERTIFIED_X);
    assert!(cap >= 17, "at least 16 modes are admitted (cap {cap})");
    assert!(ModeCount(cap).validate().is_ok());
    assert!(ModeCount(cap + 1).validate().is_err());
    assert!(ModeCount(0).validate().is_err());
    let profile = RadiusProfile::power(0.05, 1.0, 0.0, 0.01).unwrap();
    let call = |modes: usize, steps: usize| {
        mm_profile_impedance(
            &profile,
            &air20(),
            tau_hz(500.0),
            LossModel::Lossless,
            MouthLoad::Matched,
            ModeCount(modes),
            RiccatiSteps(steps),
        )
    };
    assert!(matches!(
        call(cap + 1, 64),
        Err(DuctError::BadParameter { .. })
    ));
    assert!(call(1, 0).is_err() && call(1, MAX_RICCATI_STEPS + 1).is_err());
    // Evanescent guard: a 0.1 m cylinder of radius 10 mm, N = 2, at 100 Hz.
    let cyl = RadiusProfile::power(0.01, 1.0, 0.0, 0.1).unwrap();
    let k = tau_hz(100.0) / air20().sound_speed;
    let gamma1 = j1_roots_unchecked(1)[0];
    let kappa = ((gamma1 / 0.01) * (gamma1 / 0.01) - k * k).sqrt();
    let admitted = (kappa * 0.1 / MAX_EVANESCENT_STEP).ceil() as usize;
    let run = |steps: usize| {
        mm_profile_impedance(
            &cyl,
            &air20(),
            tau_hz(100.0),
            LossModel::Lossless,
            MouthLoad::Matched,
            ModeCount(2),
            RiccatiSteps(steps),
        )
    };
    assert!(
        run(admitted).is_ok(),
        "kappa h = {} admitted",
        kappa * 0.1 / admitted as f64
    );
    assert!(
        run(admitted - 1).is_err(),
        "kappa h = {} refused",
        kappa * 0.1 / (admitted - 1) as f64
    );
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"mm-mode-cap-refusal\",\"observed\":{{\"cap\":{cap},\
         \"guard_steps_admitted\":{admitted}}},\"verdict\":\"pass\"}}"
    );
}

#[test]
fn profile_slope_matches_central_differences() {
    // G0: analytic dr/dx for every profile shape, and C1 at spline knots.
    let stations: Vec<(f64, f64)> = (0..=12)
        .map(|i| {
            let x = 0.05 * f64::from(i);
            (x, 0.01 + 0.004 * x + 0.03 * x * x * x)
        })
        .collect();
    let shapes = [
        exp_horn(),
        RadiusProfile::power(0.007, 0.08, 1.7, 0.4).unwrap(),
        RadiusProfile::spline(stations, Some([0.004, 0.004 + 0.09 * 0.36])).unwrap(),
    ];
    let mut worst = 0.0_f64;
    for p in &shapes {
        let d = 1e-6;
        for i in 1..40 {
            let x = p.length() * f64::from(i) / 40.0;
            let fd = (p.radius(x + d).unwrap() - p.radius(x - d).unwrap()) / (2.0 * d);
            let s = p.slope(x).unwrap();
            worst = worst.max((s - fd).abs() / s.abs().max(1e-3));
        }
        assert!(p.slope(-1e-9).is_err() && p.slope(p.length() + 1e-9).is_err());
    }
    let spline = &shapes[2];
    for i in 1..12 {
        let x = 0.05 * f64::from(i);
        let (l, r) = (
            spline.slope(x - 1e-12).unwrap(),
            spline.slope(x + 1e-12).unwrap(),
        );
        assert!((l - r).abs() < 1e-9, "spline slope continuous at knot {x}");
    }
    assert!(worst < 1e-7, "slope vs central difference {worst:e}");
}

#[test]
fn mm_riccati_n1_equals_plane_riccati() {
    // G2: one mode is the plane-wave Webster/Riccati of M2.1-W3a.
    let profile = exp_horn();
    let state = air20();
    let r_mouth = profile.outlet_radius();
    let mut worst = 0.0_f64;
    for loss in [LossModel::Lossless, LossModel::Bessel] {
        for f in [150.0, 700.0, 2300.0] {
            let omega = tau_hz(f);
            for (name, load, z_load) in [
                (
                    "matched",
                    MouthLoad::Matched,
                    segment_wave(&state, r_mouth, omega, loss)
                        .unwrap()
                        .characteristic_impedance,
                ),
                (
                    "flanged",
                    MouthLoad::Analytic(Termination::FlangedOpen),
                    termination_impedance(Termination::FlangedOpen, &state, r_mouth, omega)
                        .unwrap()
                        .0
                        .unwrap(),
                ),
            ] {
                let mm = z00(&profile, f, loss, load, 1, 400);
                let plane =
                    profile_impedance(&profile, &state, omega, loss, z_load, RiccatiSteps(400))
                        .unwrap()
                        .z_in;
                let e = rel(mm, plane);
                worst = worst.max(e);
                println!(
                    "{{\"suite\":\"fs-duct\",\"case\":\"mm-riccati-n1-plane\",\"inputs\":{{\"f_hz\":{f},\
                     \"loss\":\"{loss:?}\",\"load\":\"{name}\"}},\"modes\":1,\"steps\":400,\
                     \"z00\":[{},{}],\"rel_err\":{e:e},\"tol\":1e-10}}",
                    mm.re, mm.im
                );
            }
        }
    }
    assert!(worst < 1e-10, "N = 1 vs plane Riccati {worst:e}");
}

#[test]
fn mm_riccati_cylinder_decouples() {
    // G2: R' = 0 decouples the modes into the staircase's uniform waves.
    let (r, length, f) = (0.007, 0.4, 500.0);
    let cyl = RadiusProfile::power(r, 1.0, 0.0, length).unwrap();
    let state = air20();
    let omega = tau_hz(f);
    let loss = LossModel::WideTube;
    let steps = guarded_steps(length, r, 4);
    let z = mm_profile_impedance(
        &cyl,
        &state,
        omega,
        loss,
        MouthLoad::Analytic(Termination::Closed),
        ModeCount(4),
        RiccatiSteps(steps),
    )
    .unwrap()
    .impedance_matrix;
    let wave = segment_wave(&state, r, omega, loss).unwrap();
    // Closed end under e^{-i omega t}: Z = i Zc cot(k L).
    let kl = wave.wavenumber.scale(length);
    let (e_pos, e_neg) = (
        complex_exp(C64::new(-kl.im, kl.re)),
        complex_exp(C64::new(kl.im, -kl.re)),
    );
    let cot = C64::new(0.0, 1.0) * (e_pos + e_neg) * (e_pos - e_neg).recip();
    let closed = C64::new(0.0, 1.0) * wave.characteristic_impedance * cot;
    let zc = modal_characteristic_impedances(&state, r, omega, loss, 4).unwrap();
    let mut diag = rel(z[0], closed);
    // These higher modes have kappa*L > 200: their closed-end return has
    // decayed below roundoff, leaving the characteristic impedance locally.
    for n in 1..4 {
        diag = diag.max(rel(z[n * 4 + n], zc[n]));
    }
    let off = (0..16)
        .filter(|i| i / 4 != i % 4)
        .map(|i| z[i].abs())
        .fold(0.0, f64::max);
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"mm-riccati-cylinder\",\"modes\":4,\"steps\":{steps},\
         \"z00\":[{},{}],\"observed\":{{\"diag_rel\":{diag:e},\"off_diag_abs\":{off:e}}},\"tol\":1e-12}}",
        z[0].re, z[0].im
    );
    assert!(diag < 1e-12, "diagonal vs closed form {diag:e}");
    assert!(off <= 1e-14 * z[0].abs(), "off-diagonals {off:e}");
}

/// Staircase `Z00` at slice multipliers 4/8/16/32 for a duct, its observed
/// order through fs-mms, and the first-order Richardson limit. The coarser
/// 1/2/4/8 ladder measured order 1.38 on the cone (pre-asymptotic: first-
/// and second-order terms mixed), so rungs were added, never the band widened.
fn staircase_ladder(duct: &Duct, f: f64, modes: usize) -> (Vec<C64>, f64, C64) {
    let z: Vec<C64> = [4_usize, 8, 16, 32]
        .iter()
        .map(|&e| {
            mm_input_impedance(
                duct,
                &air20(),
                tau_hz(f),
                LossModel::WideTube,
                Termination::FlangedOpen,
                modes,
                e,
            )
            .expect("staircase")
            .plane_impedance
        })
        .collect();
    let diffs: Vec<f64> = z
        .windows(2)
        .map(|w| (w[0] - w[1]).abs() / w[1].abs())
        .collect();
    let ladder = RefinementLadder::new(vec![1.0, 0.5, 0.25], diffs).expect("ladder");
    let verdict = OrderGate { theoretical: 1.0 }
        .check("mm-staircase-order", LadderSide::Primal, &ladder)
        .unwrap_or_else(|e| panic!("staircase is not O(dx): {e}"));
    println!("{}", verdict.json_line(true));
    let limit = z[3].scale(2.0) - z[2];
    (z, verdict.fit.observed, limit)
}

#[test]
fn mm_riccati_matches_extrapolated_staircase() {
    // G1/G3: the staircase converges at O(dx); its Richardson limit is the
    // continuous result (an independent discretization of the same model).
    let modes = 4;
    let f = 900.0;
    // Linear cone 6 -> 30 mm over 0.3 m.
    let cone_duct = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: 0.006,
            outlet_radius: 0.03,
            length: 0.3,
        }],
    };
    let cone = RadiusProfile::power(0.006, 0.006 * 0.3 / 0.024, 1.0, 0.3).unwrap();
    // Exponential flare as 400 exact-station cones; the polygon differs from
    // the curve by ~ dx^2 r''/8 ~ 1e-6 relative, far inside the 1e-4 band.
    let horn = exp_horn();
    let k = 400;
    let horn_duct = Duct {
        segments: (0..k)
            .map(|i| {
                let (xa, xb) = (0.5 * i as f64 / k as f64, 0.5 * (i + 1) as f64 / k as f64);
                Segment::Cone {
                    inlet_radius: horn.radius(xa).unwrap(),
                    outlet_radius: horn.radius(xb).unwrap(),
                    length: xb - xa,
                }
            })
            .collect(),
    };
    for (name, duct, profile) in [
        ("cone", &cone_duct, &cone),
        ("exponential", &horn_duct, &horn),
    ] {
        let (z, order, limit) = staircase_ladder(duct, f, modes);
        let steps = guarded_steps(profile.length(), profile.inlet_radius(), modes).max(1024);
        let cont = z00(
            profile,
            f,
            LossModel::WideTube,
            MouthLoad::Analytic(Termination::FlangedOpen),
            modes,
            steps,
        );
        let cont2 = z00(
            profile,
            f,
            LossModel::WideTube,
            MouthLoad::Analytic(Termination::FlangedOpen),
            modes,
            2 * steps,
        );
        let e = rel(limit, cont);
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"mm-riccati-vs-staircase\",\"inputs\":{{\"fixture\":\"{name}\",\
             \"f_hz\":{f}}},\"modes\":{modes},\"steps\":{steps},\"z00\":[{},{}],\"observed\":{{\
             \"staircase_order\":{order:.3},\"staircase_x1\":[{},{}],\"richardson_rel_err\":{e:e},\
             \"step_refinement_rel\":{:e}}},\"tol\":1e-4}}",
            cont.re,
            cont.im,
            z[0].re,
            z[0].im,
            rel(cont, cont2)
        );
        assert!(
            rel(cont, cont2) < 1e-9,
            "{name}: continuous result is step-converged"
        );
        assert!(e < 1e-4, "{name}: Richardson staircase vs continuous {e:e}");
    }
}

#[test]
fn mm_riccati_step_ladder() {
    // G1: fourth-order Magnus, gated by fs-mms (+-0.2) on a gentle flare.
    let profile = RadiusProfile::exponential(0.02, 2.0 * det::ln(3.0) / 0.5, 0.5).unwrap();
    let (f, modes) = (1500.0, 4);
    let load = MouthLoad::Analytic(Termination::FlangedOpen);
    let reference = z00(&profile, f, LossModel::WideTube, load, modes, 4096);
    let rungs = [32_usize, 64, 128, 256];
    let errors: Vec<f64> = rungs
        .iter()
        .map(|&s| {
            rel(
                z00(&profile, f, LossModel::WideTube, load, modes, s),
                reference,
            )
        })
        .collect();
    for (s, e) in rungs.iter().zip(&errors) {
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"mm-riccati-step-ladder\",\"axis\":\"steps\",\"rung\":{s},\
             \"modes\":{modes},\"rel_err\":{e:e},\"reference_steps\":4096}}"
        );
    }
    assert!(
        errors[0] > 1e-6,
        "liveness: the coarse rung must show discretization error"
    );
    let ladder =
        RefinementLadder::new(rungs.iter().map(|&s| 1.0 / s as f64).collect(), errors).unwrap();
    let verdict = OrderGate { theoretical: 4.0 }
        .check("mm-riccati-step-ladder", LadderSide::Primal, &ladder)
        .unwrap_or_else(|e| panic!("{e}"));
    println!("{}", verdict.json_line(true));
}

/// The trombone bell of Amir, Pagneux & Kergomard (1997, Part II), Table I:
/// the only numeric geometry in either paper. The table header prints both
/// rows as [cm], but the radii are millimetres (the text gives a 25.4 mm
/// throat diameter and Fig. 1 a ~0.115 m mouth); distances are cm. Joined
/// here by the C2 spline (the paper used local cubics through 4 points).
fn trombone_bell() -> RadiusProfile {
    const TABLE_I: [(f64, f64); 22] = [
        (0.0, 12.25),
        (5.5, 12.25),
        (5.6, 12.25),
        (5.7, 12.25),
        (5.8, 12.65),
        (5.9, 12.65),
        (8.3, 13.0),
        (12.4, 14.0),
        (16.2, 15.0),
        (20.3, 16.5),
        (25.3, 18.5),
        (31.2, 21.5),
        (36.0, 24.5),
        (40.0, 28.0),
        (43.2, 32.0),
        (45.7, 37.0),
        (48.2, 44.5),
        (50.0, 52.0),
        (51.5, 62.0),
        (52.7, 74.5),
        (53.6, 89.5),
        (54.4, 114.5),
    ];
    RadiusProfile::spline(
        TABLE_I
            .iter()
            .map(|&(cm, mm)| (cm * 1e-2, mm * 1e-3))
            .collect(),
        None,
    )
    .unwrap()
}

/// Safety factor of the declared algebraic mode tail.
const MODE_TAIL_SAFETY: f64 = 2.0;

/// The declared algebraic mode-tail model: `Z(N) - Z(inf) ~ C N^-p`, with `p`
/// fitted by `fs_mms::fit_order` on the increments at N = 2..16. `deltas`
/// holds `|dZ|` for N = 1->2, 2->4, 4->8 and 8->16. Returns the fit and the
/// tail estimate `|Z(16) - Z(inf)|` before the safety factor.
fn mode_tail(deltas: &[f64]) -> (fs_mms::OrderFit, f64) {
    let fit = fit_order(
        &RefinementLadder::new(vec![0.5, 0.25, 0.125], deltas[1..].to_vec()).expect("ladder"),
    );
    let tail_16 = deltas[3] / (fs_math::det::pow(2.0, fit.observed) - 1.0);
    (fit, tail_16)
}

/// A bracketed |Z00| maximum: walk until the middle sample beats both
/// neighbours, then contract by parabolic vertices kept inside the bracket.
fn refined_peak(eval: &dyn Fn(f64) -> f64, guess: f64, half: f64) -> f64 {
    let mut half = half;
    let mut b = guess;
    let (mut fa, mut fb, mut fc) = (eval(b - half), eval(b), eval(b + half));
    let walk = |b: &mut f64, fa: &mut f64, fb: &mut f64, fc: &mut f64, half: f64| {
        for _ in 0..64 {
            if *fb >= *fa && *fb >= *fc {
                return;
            }
            if *fa > *fc {
                *b -= half;
                (*fc, *fb) = (*fb, *fa);
                *fa = eval(*b - half);
            } else {
                *b += half;
                (*fa, *fb) = (*fb, *fc);
                *fc = eval(*b + half);
            }
        }
        panic!("no |Z00| maximum bracketed near {guess} Hz");
    };
    walk(&mut b, &mut fa, &mut fb, &mut fc, half);
    for _ in 0..3 {
        b += 0.5 * half * (fa - fc) / (fa - 2.0 * fb + fc);
        half *= 0.25;
        (fa, fb, fc) = (eval(b - half), eval(b), eval(b + half));
        walk(&mut b, &mut fa, &mut fb, &mut fc, half);
    }
    b + 0.5 * half * (fa - fc) / (fa - 2.0 * fb + fc)
}

#[test]
fn mm_riccati_mode_ladder() {
    // G1: N = 1, 2, 4, 8, 16 on the published trombone bell; |dZ00| shrinks,
    // and a declared tail bound on Z00(16) is checked at N = 20. Liveness:
    // N = 1 vs N = 8 must move an impedance peak by >= 5 c (the plane image
    // misses mode conversion).
    let bell = trombone_bell();
    let r_min = 0.01225;
    let load = MouthLoad::Analytic(Termination::FlangedOpen);
    let f = 1000.0;
    let mut z = Vec::new();
    for &n in &[1_usize, 2, 4, 8, 16] {
        let steps = guarded_steps(bell.length(), r_min, n).max(256);
        let start = Instant::now();
        let resp = mm_profile_impedance(
            &bell,
            &air20(),
            tau_hz(f),
            LossModel::WideTube,
            load,
            ModeCount(n),
            RiccatiSteps(steps),
        )
        .unwrap();
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"mm-riccati-mode-ladder\",\"inputs\":\"Amir 1997 Table I bell\",\
             \"axis\":\"modes\",\"rung\":{n},\"steps\":{steps},\"z00\":[{},{}],\"tail_bound_estimate\":{:?},\
             \"elapsed_us\":{}}}",
            resp.plane_impedance.re,
            resp.plane_impedance.im,
            resp.tail_bound_estimate,
            start.elapsed().as_micros()
        );
        z.push(resp.plane_impedance);
    }
    let deltas: Vec<f64> = z.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"mm-riccati-mode-ladder\",\"observed\":{{\"abs_dz00\":{deltas:?}}}}}"
    );
    assert!(
        deltas.windows(2).all(|w| w[1] < w[0]),
        "mode ladder must shrink monotonically: {deltas:?}"
    );
    // Declared tail bound. Model the tail as algebraic, Z00(N) − Z00(∞) ≈ C·N^−p,
    // with p fitted on the last three increments (N = 2..16). Then
    // |Z00(16) − Z00(∞)| ≤ MODE_TAIL_SAFETY·|ΔZ(8→16)|/(2^p − 1). The bound is
    // falsified inside the admitted range: N = 20 must move Z00 by no more than
    // the bound, and by the fitted tail's own prediction within a factor 3.
    let (fit, tail_16) = mode_tail(&deltas);
    let order = fit.observed;
    let bound_16 = MODE_TAIL_SAFETY * tail_16;
    let steps_20 = guarded_steps(bell.length(), r_min, 20).max(256);
    let z20 = mm_profile_impedance(
        &bell,
        &air20(),
        tau_hz(f),
        LossModel::WideTube,
        load,
        ModeCount(20),
        RiccatiSteps(steps_20),
    )
    .unwrap()
    .plane_impedance;
    let moved = (z20 - z[4]).abs();
    let predicted = tail_16 * (1.0 - fs_math::det::pow(16.0 / 20.0, order));
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"mm-riccati-mode-tail\",\"inputs\":\"Amir 1997 Table I bell\",\
         \"axis\":\"modes\",\"rungs\":[2,4,8,16,20],\"observed\":{{\"order\":{order:.4},\"rms_residual\":{:.3e},\
         \"tail_bound_16\":{bound_16:.6e},\"dz00_16_to_20\":{moved:.6e},\"predicted\":{predicted:.6e},\
         \"ratio\":{:.4}}},\"tolerance\":{{\"safety\":{MODE_TAIL_SAFETY},\"prediction_band\":[0.3333,3.0]}}}}",
        fit.rms_residual,
        moved / predicted
    );
    assert!(
        order > 0.5,
        "the mode ladder must converge algebraically: order {order}"
    );
    assert!(
        moved <= bound_16,
        "declared tail bound {bound_16:e} is falsified by |Z(20) − Z(16)| = {moved:e}"
    );
    assert!(
        (1.0 / 3.0..=3.0).contains(&(moved / predicted)),
        "N = 20 departs from the fitted algebraic tail: {moved:e} vs {predicted:e}"
    );
    // Liveness: every plane-image peak below 1.5 kHz, refined for N = 1 and
    // N = 8 at one shared step count (so step error cannot pose as physics).
    let steps = guarded_steps(bell.length(), r_min, 8).max(256);
    let mag = |n: usize| {
        move |hz: f64| z00(&trombone_bell(), hz, LossModel::WideTube, load, n, steps).abs()
    };
    let (m1, m8) = (mag(1), mag(8));
    let mut coarse = Vec::new();
    let (mut a, mut b) = (m1(60.0), m1(70.0));
    let mut hz: f64 = 70.0;
    while hz < 1500.0 {
        let c = m1(hz + 10.0);
        if b > a && b > c {
            coarse.push(hz);
        }
        (a, b, hz) = (b, c, hz + 10.0);
    }
    let mut shift = 0.0_f64;
    for &g in &coarse {
        let p1 = refined_peak(&m1, g, 2.0);
        let p8 = refined_peak(&m8, p1, 2.0);
        let c = cents(p8, p1);
        shift = shift.max(c.abs());
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"mm-riccati-mode-liveness\",\"inputs\":{{\"peak_n1_hz\":{p1:.4}}},\
             \"steps\":{steps},\"observed\":{{\"peak_n8_hz\":{p8:.4},\"cents\":{c:.3}}}}}"
        );
    }
    assert!(
        coarse.len() >= 3,
        "expected peaks below 1.5 kHz, found {}",
        coarse.len()
    );
    assert!(
        shift >= 5.0,
        "liveness: N = 1 vs N = 8 must move a peak by >= 5 c (max {shift:.3})"
    );
}

/// Exact one-mode (Webster) state `(P, U)` at `x` of the exponential horn
/// `S = S0 e^{m x}`, lossless, with a mouth impedance `z_l` (Pagneux et al.
/// 1996, Part I, App. E), up to a common factor.
/// `P = e^{-mx/2}(A e^{i b x} + B e^{-i b x})`, `b^2 = k^2 - m^2/4`,
/// `U = S P'/(i omega rho)`.
fn webster_exponential_state(
    r0: f64,
    m: f64,
    length: f64,
    state: &GasState,
    omega: f64,
    z_l: C64,
    x: f64,
) -> (C64, C64) {
    let k = omega / state.sound_speed;
    let beta = C64::from_re(k * k - 0.25 * m * m).sqrt();
    let i = C64::new(0.0, 1.0);
    let lam1 = C64::from_re(-0.5 * m) + i * beta;
    let lam2 = C64::from_re(-0.5 * m) - i * beta;
    let area = |x: f64| core::f64::consts::PI * r0 * r0 * det::exp(m * x);
    let flow = |x: f64| (i.scale(omega * state.density)).recip().scale(area(x));
    let (f1, f2) = (
        complex_exp(lam1.scale(length)),
        complex_exp(lam2.scale(length)),
    );
    let cl = flow(length);
    let det = f1 * f2 * cl * (lam2 - lam1);
    let a = (z_l * cl * lam2 * f2 - f2) * det.recip();
    let b = (f1 - cl * lam1 * f1 * z_l) * det.recip();
    let (e1, e2) = (complex_exp(lam1.scale(x)), complex_exp(lam2.scale(x)));
    (a * e1 + b * e2, flow(x) * (lam1 * a * e1 + lam2 * b * e2))
}

/// Exact one-mode input impedance of the exponential horn, `P(0)/U(0)`.
fn webster_exponential(
    r0: f64,
    m: f64,
    length: f64,
    state: &GasState,
    omega: f64,
    z_l: C64,
) -> C64 {
    let (p, u) = webster_exponential_state(r0, m, length, state, omega, z_l, 0.0);
    p * u.recip()
}

/// Amir, Pagneux & Kergomard (1995), "Wave propagation in acoustic horns
/// through modal decomposition", Computational Acoustics, WIT Trans. Built
/// Environment 10, 37-44 (doi 10.2495/COMAC950051, open access), Fig. 1
/// (p. 42): `|p|` on the axis of an exponential horn, radius 0.5 -> 4 m over
/// 4.5 m, ended by an infinite cylinder and driven by a plane piston at the
/// small end, at 0.3 f_c (top panel) and 0.9 f_c (bottom panel), every 0.5 m
/// from the piston. The figure plots ANSYS (solid line) and a 15-mode modal
/// solution (circles); the paper prints no table. Where the two separate
/// (z >= 3 m at 0.9 f_c, the circles about 0.01 below the line) the reads
/// follow the ANSYS line, a converged finite-element reference.
///
/// Reads (OD17, figure-digitized, lane P8 2026-10-22): a 220 dpi render of
/// the open PDF with a calibrated overlay grid. The stated half-width is
/// 0.005 (top) and 0.015 (bottom) in the plotted units. The paper leaves two
/// inputs implicit; both were identified from its own horn-equation (dashed)
/// curve, which the exact Webster solution reproduces to the read error:
/// f_c is the cutoff of the first axisymmetric mode at the mouth
/// (`k a = gamma_1 = 3.8317`; the 1.8412 reading does not fit), and the
/// plotted quantity is `|p|/(rho c q)` for a piston volume flow `q`.
/// Rows: `(f/f_c, read half-width, |p| at z = 0, 0.5, .., 4.5 m)`.
const WIT_1995_FIG1: [(f64, f64, [f64; 10]); 2] = [
    (
        0.3,
        0.005,
        [
            0.479, 0.331, 0.225, 0.153, 0.109, 0.077, 0.056, 0.045, 0.036, 0.031,
        ],
    ),
    (
        0.9,
        0.015,
        [
            1.159, 0.765, 0.490, 0.395, 0.400, 0.410, 0.390, 0.341, 0.285, 0.209,
        ],
    ),
];

#[test]
fn mm_riccati_pagneux_cases() {
    // G2. Pagneux et al. (1996, Part I; 1997, Part II) publish their
    // benchmark results as figures only (checked 2026-10-20 on the authors'
    // site), so no 1 % comparison against printed numbers is claimed here.
    // What is reproducible: (a) the exact one-mode solution of Part I App. E
    // for the exponential horn, above and below the horn cutoff, (b) the
    // Part I 45-degree truncated cone (0.08 -> 0.30 m, infinite-pipe end,
    // plane piston, 600 Hz, 7 and 17 modes), whose results must be
    // reciprocal and settle with the mode count, and (c) the one open
    // figure by the same authors, [`WIT_1995_FIG1`], against the converged
    // value (reported, not gated).
    let state = air20();
    let horn = exp_horn();
    let m = 2.0 * det::ln(10.0) / 0.5;
    let f_cut = state.sound_speed * m / (4.0 * core::f64::consts::PI);
    let mut worst = 0.0_f64;
    for f in [0.5 * f_cut, 500.0, 1000.0, 2000.0] {
        let omega = tau_hz(f);
        let z_l = segment_wave(&state, horn.outlet_radius(), omega, LossModel::Lossless)
            .unwrap()
            .characteristic_impedance;
        let exact = webster_exponential(0.006, m, 0.5, &state, omega, z_l);
        let mm = z00(&horn, f, LossModel::Lossless, MouthLoad::Matched, 1, 8192);
        let e = rel(mm, exact);
        worst = worst.max(e);
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"mm-riccati-pagneux-webster\",\"inputs\":{{\"fixture\":\
             \"exponential horn, Part I App. E\",\"f_hz\":{f:.2},\"f_cut_hz\":{f_cut:.2}}},\"modes\":1,\
             \"steps\":8192,\"z00\":[{},{}],\"oracle\":[{},{}],\"rel_err\":{e:e},\"tol\":1e-9}}",
            mm.re, mm.im, exact.re, exact.im
        );
    }
    assert!(worst < 1e-9, "one-mode exponential horn vs exact {worst:e}");
    let cone = RadiusProfile::power(0.08, 0.08, 1.0, 0.22).unwrap();
    let mut z_by_n = Vec::new();
    for modes in [7_usize, 17] {
        let steps = guarded_steps(0.22, 0.08, modes).max(64);
        let resp = mm_profile_impedance(
            &cone,
            &state,
            tau_hz(600.0),
            LossModel::Lossless,
            MouthLoad::Matched,
            ModeCount(modes),
            RiccatiSteps(steps),
        )
        .unwrap();
        let z = &resp.impedance_matrix;
        let scale = z.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let asym = (0..modes * modes)
            .map(|idx| (z[idx] - z[(idx % modes) * modes + idx / modes]).abs())
            .fold(0.0, f64::max)
            / scale;
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"mm-riccati-pagneux-cone\",\"inputs\":{{\"fixture\":\
             \"Part I 45-degree cone 0.08->0.30 m, infinite-pipe end\",\"f_hz\":600}},\"modes\":{modes},\
             \"steps\":{steps},\"z00\":[{},{}],\"observed\":{{\"asymmetry\":{asym:e},\"tail_bound\":{:?}}}}}",
            resp.plane_impedance.re, resp.plane_impedance.im, resp.tail_bound_estimate
        );
        assert!(asym < 1e-10, "reciprocity at N = {modes}: {asym:e}");
        z_by_n.push((resp.plane_impedance, resp.tail_bound_estimate.unwrap()));
    }
    // Modal convergence on a 45-degree wall is slow (the rigid-duct modes
    // meet the sloped wall's Neumann condition only in the limit), so the
    // claim is the shrinking tail, not a chosen settling band.
    let settle = rel(z_by_n[0].0, z_by_n[1].0);
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"mm-riccati-pagneux-cone\",\"observed\":{{\"n7_vs_n17_rel\":{settle:e}}}}}"
    );
    assert!(
        z_by_n[1].1 < z_by_n[0].1,
        "the tail bound must shrink from 7 to 17 modes ({:e} -> {:e})",
        z_by_n[0].1,
        z_by_n[1].1
    );
    // (c) The reads are of the ANSYS line, a converged reference, so the
    // comparison is against OUR converged value: |p| at the admitted mode
    // cap. The allowance is the read half-width plus the declared tail of
    // `mode_tail`, fitted per point on the N = 1..16 increments, carried to
    // the cap as tail_16 (16/cap)^p and multiplied by MODE_TAIL_SAFETY. Where
    // the fitted order is not above 0.5 (the mode ladder's own condition)
    // the point has no declared tail. The agreement is reported per point,
    // not gated; the 1 % DONE-WHEN item is the owner's. Gated: step
    // convergence of the compared value, and liveness (the plane image
    // misses both piston reads by more than the half-width).
    const WIT_STEPS: usize = 576; // 9 x 64: stations land on z = 0.5 i (h = 2^-7 m)
    let wit = RadiusProfile::exponential(0.5, 2.0 * det::ln(8.0) / 4.5, 4.5).unwrap();
    let gamma1 = j1_roots_unchecked(1)[0];
    let rho_c = state.density * state.sound_speed;
    let cap = certified_mode_cap();
    let axis = |modes: usize, k: f64, steps: usize| -> Vec<C64> {
        let mut piston = vec![C64::ZERO; modes];
        piston[0] = C64::ONE;
        let field = mm_profile_field(
            &wit,
            &state,
            k * state.sound_speed,
            LossModel::Lossless,
            MouthLoad::Matched,
            ModeCount(modes),
            RiccatiSteps(steps),
            &piston,
        )
        .unwrap();
        (0..10)
            .map(|i| {
                let j = i * steps / 9;
                assert!((field.x[j] - 0.5 * i as f64).abs() < 1e-12, "station on z");
                field.axis_pressure[j].scale(1.0 / rho_c)
            })
            .collect()
    };
    let (mut agree, mut miss, mut untailed, mut plane_outside) =
        (0_usize, 0_usize, 0_usize, 0_usize);
    for &(fraction, half_width, reads) in &WIT_1995_FIG1 {
        let k = fraction * gamma1 / 4.0;
        let start = Instant::now();
        let ladder: Vec<Vec<C64>> = [1_usize, 2, 4, 8, 16]
            .iter()
            .map(|&n| axis(n, k, WIT_STEPS))
            .collect();
        let converged = axis(cap, k, WIT_STEPS);
        let finer = axis(cap, k, 2 * WIT_STEPS);
        for i in 0..10 {
            let deltas: Vec<f64> = ladder
                .windows(2)
                .map(|w| (w[1][i] - w[0][i]).abs())
                .collect();
            let (fit, tail_16) = mode_tail(&deltas);
            let order = fit.observed;
            let tail = (order > 0.5)
                .then(|| MODE_TAIL_SAFETY * tail_16 * det::pow(16.0 / cap as f64, order));
            let value = converged[i].abs();
            let err = (value - reads[i]).abs();
            let verdict = match tail {
                Some(t) if err <= half_width + t => "agree",
                Some(_) => "miss",
                None => "no-tail-model",
            };
            match verdict {
                "agree" => agree += 1,
                "miss" => miss += 1,
                _ => untailed += 1,
            }
            let plane = ladder[0][i].abs();
            plane_outside += usize::from((plane - reads[i]).abs() > half_width);
            let step_move = (finer[i] - converged[i]).abs();
            println!(
                "{{\"suite\":\"fs-duct\",\"case\":\"mm-riccati-pagneux-wit1995\",\"inputs\":{{\"fixture\":\
                 \"Amir-Pagneux-Kergomard 1995 Fig. 1 exp horn 0.5->4 m\",\"f_over_fc\":{fraction},\"z_m\":{}}},\
                 \"modes\":{cap},\"steps\":{WIT_STEPS},\"observed\":{{\"converged\":{value:.5},\"n16\":{:.5},\
                 \"n1\":{plane:.5},\"tail_order\":{order:.3},\"tail_bound\":{},\"step_move\":{step_move:.2e}}},\
                 \"oracle\":{},\"err\":{err:.5},\"tol\":{half_width},\"verdict\":\"{verdict}\",\"elapsed_us\":{}}}",
                0.5 * i as f64,
                ladder[4][i].abs(),
                tail.map_or("null".to_string(), |t| format!("{t:.5}")),
                reads[i],
                start.elapsed().as_micros()
            );
            assert!(
                step_move < 1e-3 * half_width,
                "WIT {fraction} f_c, z = {}: not step-converged ({step_move:e})",
                0.5 * i as f64
            );
        }
        let plane = ladder[0][0].abs();
        assert!(
            (plane - reads[0]).abs() > half_width,
            "liveness: the plane image {plane:.4} must miss the piston read {} +- {half_width}",
            reads[0]
        );
    }
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"mm-riccati-pagneux-wit1995\",\"observed\":{{\"agree\":{agree},\
         \"miss\":{miss},\"no_tail_model\":{untailed},\"n1_outside_half_width\":{plane_outside},\"reads\":20}}}}"
    );
}

#[test]
fn mm_riccati_reciprocity_and_power() {
    // G0: Z = Z^T (lossy and lossless) and, lossless with every retained mode
    // propagating, input power equals the matched mouth's absorbed power.
    let state = air20();
    let cone = RadiusProfile::power(0.08, 0.08, 1.0, 0.22).unwrap();
    let (n, f) = (4_usize, 8000.0);
    for loss in [LossModel::Lossless, LossModel::WideTube] {
        let p = propagate(
            &cone,
            &state,
            tau_hz(f),
            loss,
            MouthLoad::Matched,
            n,
            RiccatiSteps(512),
            true,
        )
        .unwrap();
        let scale = p.zeta.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let asym = (0..n * n)
            .map(|idx| (p.zeta[idx] - p.zeta[(idx % n) * n + idx / n]).abs())
            .fold(0.0, f64::max)
            / scale;
        assert!(asym < 1e-10, "{loss:?}: reciprocity {asym:e}");
        if loss != LossModel::Lossless {
            continue;
        }
        let g = p.flow.unwrap();
        let g_inv = mat_inv(&g, n).unwrap();
        let mouth = modal_wave(
            &state,
            cone.outlet_radius(),
            tau_hz(f),
            loss,
            &j1_roots_unchecked(n - 1),
        )
        .unwrap();
        let mut worst = 0.0_f64;
        for drive in [vec![1.0, 0.0, 0.0, 0.0], vec![0.3, -0.7, 0.2, 0.5]] {
            let u_in: Vec<C64> = drive.iter().map(|&v| C64::from_re(v)).collect();
            let power_in: f64 = (0..n)
                .map(|r| {
                    let pr = (0..n).fold(C64::ZERO, |acc, c| {
                        acc + p.zeta[r * n + c].scale(p.zref) * u_in[c]
                    });
                    (u_in[r].conj() * pr).re
                })
                .sum();
            let u_out: Vec<C64> = (0..n)
                .map(|r| {
                    (0..n).fold(C64::ZERO, |acc, c| {
                        acc + g_inv[r * n + c] * u_in[c].scale(1.0 / p.zref)
                    })
                })
                .map(|v| v.scale(p.zref))
                .collect();
            let power_out: f64 = (0..n)
                .map(|r| (u_out[r].conj() * mouth.zc[r] * u_out[r]).re)
                .sum();
            let e = (power_in - power_out).abs() / power_in.abs();
            worst = worst.max(e);
        }
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"mm-riccati-reciprocity-power\",\"inputs\":{{\"fixture\":\
             \"Pagneux Part I cone\",\"f_hz\":{f}}},\"modes\":{n},\"steps\":512,\"observed\":{{\
             \"asymmetry\":{asym:e},\"power_rel\":{worst:e}}},\"tol\":1e-10}}"
        );
        assert!(worst < 1e-10, "lossless power balance {worst:e}");
    }
}

#[test]
fn mm_riccati_bitwise_repeat() {
    // G5: deterministic fixed-order arithmetic.
    let run = || {
        mm_profile_impedance(
            &trombone_bell(),
            &air20(),
            tau_hz(700.0),
            LossModel::WideTube,
            MouthLoad::Analytic(Termination::FlangedOpen),
            ModeCount(4),
            RiccatiSteps(guarded_steps(0.544, 0.01225, 4).max(256)),
        )
        .unwrap()
        .impedance_matrix
    };
    let (a, b) = (run(), run());
    assert!(
        a.iter()
            .zip(&b)
            .all(|(x, y)| x.re.to_bits() == y.re.to_bits() && x.im.to_bits() == y.im.to_bits())
    );
}

#[test]
fn mm_field_cylinder_modes_closed_form() {
    // G2: on a cylinder each mode travels alone. A matched mouth leaves one
    // outgoing (or decaying) wave per mode, `U_n(x) = U_n(0) e^{i k_n x}`,
    // `P_n = Zc_n U_n`, `k_n = sqrt(k^2 - (gamma_n/R)^2)` with `Im k_n >= 0`;
    // the axis weights `1/|J0(gamma_n)|` are checked against the quadrature.
    let (r, length, f, n) = (0.05, 0.1, 900.0, 4);
    let cyl = RadiusProfile::power(r, 1.0, 0.0, length).unwrap();
    let state = air20();
    let omega = tau_hz(f);
    let zc = modal_characteristic_impedances(&state, r, omega, LossModel::Lossless, n).unwrap();
    let k0 = omega / state.sound_speed;
    let gammas: Vec<f64> = std::iter::once(0.0)
        .chain(j1_roots_unchecked(n - 1))
        .collect();
    let drive: Vec<C64> = [1.0, 0.4, -0.3, 0.2]
        .iter()
        .map(|&v| C64::from_re(v))
        .collect();
    let steps = guarded_steps(length, r, n).max(64);
    let field = mm_profile_field(
        &cyl,
        &state,
        omega,
        LossModel::Lossless,
        MouthLoad::Matched,
        ModeCount(n),
        RiccatiSteps(steps),
        &drive,
    )
    .unwrap();
    assert_eq!(field.n_modes, n);
    assert_eq!(field.x.len(), steps + 1);
    assert!(field.x[0] == 0.0 && (field.x[steps] - length).abs() < 1e-15);
    assert!(field.x.windows(2).all(|w| w[1] > w[0]), "stations ascend");
    let i = C64::new(0.0, 1.0);
    let (mut worst, mut worst_axis) = (0.0_f64, 0.0_f64);
    for (j, &x) in field.x.iter().enumerate() {
        let mut axis = C64::ZERO;
        for m in 0..n {
            let kt = gammas[m] / r;
            let kn = C64::from_re(k0 * k0 - kt * kt).sqrt();
            let u = drive[m] * complex_exp(i * kn.scale(x));
            let p = zc[m] * u;
            worst = worst
                .max((field.flow[j][m] - u).abs() / u.abs())
                .max((field.pressure[j][m] - p).abs() / p.abs());
            axis = axis + p.scale(1.0 / quad_jn(0, gammas[m]).abs());
        }
        worst_axis = worst_axis.max((field.axis_pressure[j] - axis).abs() / axis.abs());
    }
    println!(
        "{{\"suite\":\"fs-duct\",\"case\":\"mm-field-cylinder\",\"modes\":{n},\"steps\":{steps},\
         \"observed\":{{\"modal_rel\":{worst:e},\"axis_rel\":{worst_axis:e}}},\"tol\":1e-12}}"
    );
    assert!(worst < 1e-12, "modal field vs closed form {worst:e}");
    assert!(
        worst_axis < 1e-12,
        "axis pressure vs quadrature weights {worst_axis:e}"
    );
}

#[test]
fn mm_field_n1_equals_exact_webster() {
    // G2: one mode on the exponential horn is the exact Webster field of
    // Part I App. E, station by station, below and above the horn cutoff.
    let state = air20();
    let horn = exp_horn();
    let m = 2.0 * det::ln(10.0) / 0.5;
    let f_cut = state.sound_speed * m / (4.0 * core::f64::consts::PI);
    let steps = 4096;
    let mut worst = 0.0_f64;
    for f in [0.5 * f_cut, 1000.0] {
        let omega = tau_hz(f);
        let z_l = segment_wave(&state, horn.outlet_radius(), omega, LossModel::Lossless)
            .unwrap()
            .characteristic_impedance;
        let field = mm_profile_field(
            &horn,
            &state,
            omega,
            LossModel::Lossless,
            MouthLoad::Matched,
            ModeCount(1),
            RiccatiSteps(steps),
            &[C64::ONE],
        )
        .unwrap();
        let u0 = webster_exponential_state(0.006, m, 0.5, &state, omega, z_l, 0.0).1;
        let exact: Vec<(C64, C64)> = field
            .x
            .iter()
            .map(|&x| {
                let (p, u) = webster_exponential_state(0.006, m, 0.5, &state, omega, z_l, x);
                (p * u0.recip(), u * u0.recip())
            })
            .collect();
        let p_scale = exact.iter().map(|e| e.0.abs()).fold(0.0, f64::max);
        let u_scale = exact.iter().map(|e| e.1.abs()).fold(0.0, f64::max);
        let (mut ep, mut eu) = (0.0_f64, 0.0_f64);
        for (j, &(p, u)) in exact.iter().enumerate() {
            assert!(field.axis_pressure[j] == field.pressure[j][0], "N = 1 axis");
            ep = ep.max((field.pressure[j][0] - p).abs() / p_scale);
            eu = eu.max((field.flow[j][0] - u).abs() / u_scale);
        }
        worst = worst.max(ep).max(eu);
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"mm-field-webster\",\"inputs\":{{\"fixture\":\
             \"exponential horn, Part I App. E\",\"f_hz\":{f:.2},\"f_cut_hz\":{f_cut:.2}}},\"modes\":1,\
             \"steps\":{steps},\"observed\":{{\"p_rel\":{ep:e},\"u_rel\":{eu:e}}},\"tol\":1e-9}}"
        );
    }
    assert!(worst < 1e-9, "one-mode field vs exact Webster {worst:e}");
}

#[test]
fn mm_field_power_flux_conserved() {
    // G0: lossless, the axial power flux `Re sum_n P_n conj(U_n)` (= `Re int
    // p conj(v_x) dS`) is the same at every station, whether the higher
    // modes propagate (8 kHz) or are evanescent (600 Hz). The inlet station
    // is `P = Z U_in` with `Z` from `mm_profile_impedance`, and the matched
    // mouth is `P_n = Zc_n U_n`.
    let state = air20();
    let cone = RadiusProfile::power(0.08, 0.08, 1.0, 0.22).unwrap();
    for (f, n) in [(8000.0, 4_usize), (600.0, 8)] {
        let omega = tau_hz(f);
        let steps = guarded_steps(0.22, 0.08, n).max(512);
        let drive: Vec<C64> = (0..n)
            .map(|m| C64::new(1.0 / (1.0 + m as f64), 0.3 * m as f64 - 0.5))
            .collect();
        let field = mm_profile_field(
            &cone,
            &state,
            omega,
            LossModel::Lossless,
            MouthLoad::Matched,
            ModeCount(n),
            RiccatiSteps(steps),
            &drive,
        )
        .unwrap();
        let flux = |j: usize| -> f64 {
            (0..n)
                .map(|m| (field.pressure[j][m] * field.flow[j][m].conj()).re)
                .sum()
        };
        let w_in = flux(0);
        let drift = (0..=steps)
            .map(|j| (flux(j) - w_in).abs())
            .fold(0.0, f64::max)
            / w_in.abs();
        let z = mm_profile_impedance(
            &cone,
            &state,
            omega,
            LossModel::Lossless,
            MouthLoad::Matched,
            ModeCount(n),
            RiccatiSteps(steps),
        )
        .unwrap()
        .impedance_matrix;
        let p_scale = field.pressure[0]
            .iter()
            .map(|v| v.abs())
            .fold(0.0, f64::max);
        let inlet = (0..n)
            .map(|r| {
                let zu = (0..n).fold(C64::ZERO, |acc, c| acc + z[r * n + c] * drive[c]);
                (field.pressure[0][r] - zu).abs()
            })
            .fold(0.0, f64::max)
            / p_scale;
        let zc = modal_characteristic_impedances(
            &state,
            cone.outlet_radius(),
            omega,
            LossModel::Lossless,
            n,
        )
        .unwrap();
        let mouth = (0..n)
            .map(|m| {
                let p = field.pressure[steps][m];
                (p - zc[m] * field.flow[steps][m]).abs() / p.abs()
            })
            .fold(0.0, f64::max);
        println!(
            "{{\"suite\":\"fs-duct\",\"case\":\"mm-field-power\",\"inputs\":{{\"fixture\":\
             \"Pagneux Part I cone\",\"f_hz\":{f}}},\"modes\":{n},\"steps\":{steps},\"observed\":{{\
             \"flux_drift_rel\":{drift:e},\"inlet_rel\":{inlet:e},\"mouth_rel\":{mouth:e}}},\"tol\":1e-10}}"
        );
        assert!(drift < 1e-10, "{f} Hz: power flux drift {drift:e}");
        assert!(inlet < 1e-12, "{f} Hz: inlet P vs Z U_in {inlet:e}");
        assert!(mouth < 1e-12, "{f} Hz: matched mouth {mouth:e}");
    }
}

#[test]
fn mm_field_refusals() {
    // G0: one finite inlet flow per retained mode, a finite plane-mode mouth
    // load, and the impedance image's own mode and step refusals.
    fn call(
        load: MouthLoad<'_>,
        modes: usize,
        steps: usize,
        drive: &[C64],
    ) -> Result<ModalField, DuctError> {
        mm_profile_field(
            &exp_horn(),
            &air20(),
            tau_hz(500.0),
            LossModel::Lossless,
            load,
            ModeCount(modes),
            RiccatiSteps(steps),
            drive,
        )
    }
    let bad_param =
        |r: Result<ModalField, DuctError>| matches!(r, Err(DuctError::BadParameter { .. }));
    let ok = [C64::ONE, C64::ZERO];
    assert!(call(MouthLoad::Matched, 2, 256, &ok).is_ok());
    assert!(bad_param(call(MouthLoad::Matched, 2, 256, &ok[..1])));
    assert!(bad_param(call(
        MouthLoad::Matched,
        2,
        256,
        &[C64::ONE, C64::ZERO, C64::ZERO]
    )));
    assert!(bad_param(call(
        MouthLoad::Matched,
        2,
        256,
        &[C64::new(f64::NAN, 0.0), C64::ZERO]
    )));
    assert!(bad_param(call(
        MouthLoad::Analytic(Termination::Closed),
        2,
        256,
        &ok
    )));
    let cap = certified_mode_cap();
    assert!(bad_param(call(
        MouthLoad::Matched,
        cap + 1,
        256,
        &vec![C64::ZERO; cap + 1]
    )));
    assert!(call(MouthLoad::Matched, 2, 0, &ok).is_err());
    assert!(call(MouthLoad::Matched, 2, MAX_RICCATI_STEPS + 1, &ok).is_err());
    println!("{{\"suite\":\"fs-duct\",\"case\":\"mm-field-refusals\",\"verdict\":\"pass\"}}");
}

#[test]
fn ideal_mouth_reflects_every_mode_in_both_modal_solvers() {
    // G1/G3: independently separated lossless cylinder modes. In the
    // exp(-i omega t) convention the round trip is +/- exp(2 i k_n L).
    // A short cylinder keeps evanescent returns measurable instead of
    // accidentally testing the semi-infinite (matched) limit.
    let state = air20();
    let (radius, length, modes) = (0.05, 0.015, 3);
    let profile = RadiusProfile::power(radius, 1.0, 0.0, length).unwrap();
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let mut gammas = vec![0.0];
    gammas.extend(j1_roots_unchecked(modes - 1));
    for frequency in [300.0, 5000.0, 10000.0] {
        let omega = tau_hz(frequency);
        let k0 = omega / state.sound_speed;
        let zc = modal_characteristic_impedances(
            &state, radius, omega, LossModel::Lossless, modes,
        ).unwrap();
        for termination in [Termination::Closed, Termination::IdealOpen] {
            let sign = if matches!(termination, Termination::Closed) { 1.0 } else { -1.0 };
            let staircase = mm_input_impedance(
                &duct, &state, omega, LossModel::Lossless, termination, modes, 1,
            ).unwrap();
            let continuous = mm_profile_impedance(
                &profile, &state, omega, LossModel::Lossless,
                MouthLoad::Analytic(termination), ModeCount(modes), RiccatiSteps(32),
            ).unwrap();
            for (name, response) in [("staircase", staircase), ("continuous", continuous)] {
                let reflection = crate::modal::modal_reflection_from_impedance(
                    &response.impedance_matrix, &zc, modes,
                ).unwrap();
                for (mode, &gamma) in gammas.iter().enumerate() {
                    let axial_squared = k0 * k0 - (gamma / radius).powi(2);
                    let expected = if axial_squared >= 0.0 {
                        let phase = 2.0 * axial_squared.sqrt() * length;
                        C64::new(phase.cos(), phase.sin()).scale(sign)
                    } else {
                        C64::from_re(sign * (-2.0 * (-axial_squared).sqrt() * length).exp())
                    };
                    for column in 0..modes {
                        let target = if mode == column { expected } else { C64::ZERO };
                        let actual = reflection[mode * modes + column];
                        assert!(
                            (actual - target).abs() < 5e-11,
                            "{name}, {termination:?}, {frequency} Hz, ({mode},{column}): {actual:?} != {target:?}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn pressure_release_field_annihilates_every_mouth_pressure_mode() {
    // G1: P_n(L)=0 for each independently driven radial mode; the
    // corresponding lossless cylinder flow is U_n(L)=U_n(0)/cos(k_n L).
    let state = air20();
    let (radius, length, modes) = (0.05, 0.015, 3);
    let profile = RadiusProfile::power(radius, 1.0, 0.0, length).unwrap();
    let drive = [C64::ONE, C64::new(0.2, 0.1), C64::new(-0.3, 0.4)];
    let mut gammas = vec![0.0];
    gammas.extend(j1_roots_unchecked(modes - 1));
    for frequency in [300.0, 5000.0, 10000.0] {
        let omega = tau_hz(frequency);
        let field = mm_profile_field(
            &profile, &state, omega, LossModel::Lossless,
            MouthLoad::Analytic(Termination::IdealOpen), ModeCount(modes),
            RiccatiSteps(32), &drive,
        ).unwrap();
        let pressure = field.pressure.last().unwrap();
        let flow = field.flow.last().unwrap();
        for (mode, &gamma) in gammas.iter().enumerate() {
            assert_eq!(pressure[mode], C64::ZERO, "mouth pressure mode {mode}");
            let k0 = omega / state.sound_speed;
            let axial_squared = k0 * k0 - (gamma / radius).powi(2);
            let cosine = if axial_squared >= 0.0 {
                (axial_squared.sqrt() * length).cos()
            } else {
                ((-axial_squared).sqrt() * length).cosh()
            };
            let expected = drive[mode].scale(1.0 / cosine);
            assert!(
                (flow[mode] - expected).abs() < 5e-11 * expected.abs(),
                "{frequency} Hz, mode {mode}: {:?} != {expected:?}", flow[mode]
            );
        }
    }
}
