//! fs-solid structural-elements battery (bead tfz.14).
//!
//! - str-001 G0: rod strain objectivity — a superposed rigid motion
//!   leaves the strain energy EXACTLY invariant.
//! - str-002 G2: elastica — the large-deflection cantilever tip
//!   against an in-test shooting/elliptic-integral oracle.
//! - str-003 G2: the helical family — pure bending (circle), pure
//!   torsion (straight + twist), combined (helix curvature/torsion).
//! - str-004: RC fiber-section moment-curvature hysteresis — closed
//!   loops, positive and growing per-cycle dissipation, peak moment
//!   within a hand capacity band, bitwise section-state determinism.
//! - str-005: batched section updates — bitwise-consistent with the
//!   scalar path, measured throughput ledgered (the §15.2 pairing).
//! - str-006: force-based pushover — post-yield stiffness loss,
//!   hysteretic dissipation under a reversal, and G4 checkpoint/resume
//!   bitwise equality mid-history.

use fs_solid::fiber::rc_section;
use fs_solid::{ForceBasedElement, Rod, RodSection, TipLoad, update_sections_batched};
use std::fmt::Write as _;
use std::time::Instant;

fn verdict(name: &str, pass: bool, details: &str) {
    println!(
        "{{\"test\":\"{name}\",\"verdict\":\"{}\",{details}}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "{name} failed: {details}");
}

// Slenderness EA·L²/EI = 200: stiff enough that axial/shear
// corrections sit at ~0.1% (below every oracle gate), soft enough
// that Newton steps are not rejected by spurious-stretch penalty
// (EA ~ 1e4 puts a quadratic wall around every iterate — measured).
const SECTION: RodSection = RodSection {
    ea: 200.0,
    ga: 200.0,
    gj: 1.0,
    ei: 1.0,
};

// ------------------------------------------------------------------ str-001

#[test]
fn str_001_rod_objectivity() {
    use fs_time::lie::{quat_exp, quat_mul, quat_rotate};
    let mut rod = Rod::straight(1.0, 12, SECTION);
    // A visibly deformed state (smooth bend + twist).
    for (i, q) in rod.quats.iter_mut().enumerate() {
        #[allow(clippy::cast_precision_loss)]
        let s = i as f64 / 12.0;
        *q = quat_exp([0.3 * s, 0.4 * s, 0.1 * s * s]);
    }
    for (i, p) in rod.positions.iter_mut().enumerate() {
        #[allow(clippy::cast_precision_loss)]
        let s = i as f64 / 12.0;
        *p = [s, 0.2 * s * s, 0.1 * s * s * s];
    }
    let e0 = rod.energy();
    // Superpose a rigid motion: rotate everything by R, translate by t.
    let rq = quat_exp([0.4, -0.7, 0.5]);
    let t = [3.0, -2.0, 1.0];
    let mut moved = rod.clone();
    for p in &mut moved.positions {
        let r = quat_rotate(rq, *p);
        *p = [r[0] + t[0], r[1] + t[1], r[2] + t[2]];
    }
    for q in &mut moved.quats {
        *q = quat_mul(rq, *q);
    }
    let e1 = moved.energy();
    let rel = (e1 - e0).abs() / e0.abs().max(1e-30);
    verdict(
        "str-001",
        rel < 1e-12,
        &format!("\"detail\":\"rigid motion leaves rod energy invariant\",\"rel\":{rel:.3e}"),
    );
}

// ------------------------------------------------------------------ str-002

/// Shooting oracle for the cantilever elastica: EIθ″ = −P·cosθ,
/// θ(0) = 0, θ′(L) = 0; tip = (∫cosθ, ∫sinθ).
fn elastica_tip(p_over_ei: f64, l: f64) -> (f64, f64) {
    let rhs = |theta: f64| -p_over_ei * theta.cos();
    let integrate = |th0p: f64| -> (f64, f64, f64) {
        let n = 4000i32;
        let h = l / f64::from(n);
        let (mut th, mut thp) = (0.0f64, th0p);
        let (mut x, mut y) = (0.0f64, 0.0f64);
        for _ in 0..n {
            // RK4 on (θ, θ′); trapezoid on the tip integrals.
            let (t1, p1) = (thp, rhs(th));
            let (t2, p2) = (thp + 0.5 * h * p1, rhs(th + 0.5 * h * t1));
            let (t3, p3) = (thp + 0.5 * h * p2, rhs(th + 0.5 * h * t2));
            let (t4, p4) = (thp + h * p3, rhs(th + h * t3));
            let th_new = th + h / 6.0 * (t1 + 2.0 * t2 + 2.0 * t3 + t4);
            x += 0.5 * h * (th.cos() + th_new.cos());
            y += 0.5 * h * (th.sin() + th_new.sin());
            thp += h / 6.0 * (p1 + 2.0 * p2 + 2.0 * p3 + p4);
            th = th_new;
        }
        (thp, x, y)
    };
    // Bisection on θ′(0) for θ′(L) = 0.
    let (mut a, mut b) = (0.0f64, 3.0 * p_over_ei);
    for _ in 0..80 {
        let m = f64::midpoint(a, b);
        if integrate(m).0 > 0.0 {
            b = m;
        } else {
            a = m;
        }
    }
    let (_, x, y) = integrate(f64::midpoint(a, b));
    (x, y)
}

#[test]
fn str_002_elastica_large_deflection() {
    // PL²/EI = 2 (strongly geometric, pre-loop).
    let p = 2.0;
    let (x_ref, y_ref) = elastica_tip(p, 1.0);
    let mut rod = Rod::straight(1.0, 16, SECTION);
    let load = TipLoad {
        force: [0.0, p * SECTION.ei, 0.0],
        moment: [0.0, 0.0, 0.0],
    };
    rod.solve_static(&load, 8, 1e-7)
        .expect("elastica converges");
    let tip = rod.positions[16];
    let dx = (tip[0] - x_ref).abs();
    let dy = (tip[1] - y_ref).abs();
    let pass = dx < 0.015 && dy < 0.015;
    verdict(
        "str-002",
        pass,
        &format!(
            "\"detail\":\"cantilever elastica vs shooting oracle, PL2/EI=2\",\
             \"tip\":[{:.4},{:.4}],\"oracle\":[{x_ref:.4},{y_ref:.4}]",
            tip[0], tip[1]
        ),
    );
}

// ------------------------------------------------------------------ str-003

#[test]
#[allow(clippy::too_many_lines)]
fn str_003_helical_family() {
    // (a) Pure bending: end moment about z → arc of radius EI/M.
    let m_b = 0.8f64;
    let mut arc = Rod::straight(1.0, 20, SECTION);
    arc.solve_static(
        &TipLoad {
            force: [0.0; 3],
            moment: [0.0, 0.0, m_b],
        },
        6,
        1e-8,
    )
    .expect("arc converges");
    let r_want = SECTION.ei / m_b;
    // Fit: center is at (0, R) for a rod clamped along +x bending in
    // +y; check radii of all nodes.
    let mut arc_dev = 0.0f64;
    for p in &arc.positions {
        let r = (p[0]).hypot(p[1] - r_want);
        arc_dev = arc_dev.max((r - r_want).abs() / r_want);
    }
    // (b) Pure torsion: end moment about x → straight centerline,
    // twist rate M/GJ.
    let m_t = 0.5f64;
    let mut twist = Rod::straight(1.0, 20, SECTION);
    twist
        .solve_static(
            &TipLoad {
                force: [0.0; 3],
                moment: [m_t, 0.0, 0.0],
            },
            4,
            1e-8,
        )
        .expect("twist converges");
    let mut straight_dev = 0.0f64;
    for p in &twist.positions {
        straight_dev = straight_dev.max(p[1].abs().max(p[2].abs()));
    }
    let (_, kappa_end) = twist.strains(10);
    let twist_rate_dev = (kappa_end[0] - m_t / SECTION.gj).abs() / (m_t / SECTION.gj);
    // (c) Helix: combined bending + torsion moments → constant strain
    // state κ = (m_t/GJ, 0, m_b/EI) along the rod.
    let mut helix = Rod::straight(1.0, 24, SECTION);
    helix
        .solve_static(
            &TipLoad {
                force: [0.0; 3],
                moment: [0.4, 0.0, 0.6],
            },
            8,
            1e-8,
        )
        .expect("helix converges");
    let mut kdev = 0.0f64;
    for seg in 2..22 {
        let (_, k) = helix.strains(seg);
        kdev = kdev
            .max((k[0] - 0.4 / SECTION.gj).abs() / (0.4 / SECTION.gj))
            .max((k[2] - 0.6 / SECTION.ei).abs() / (0.6 / SECTION.ei));
    }
    let pass = arc_dev < 0.02 && straight_dev < 1e-6 && twist_rate_dev < 0.02 && kdev < 0.03;
    verdict(
        "str-003",
        pass,
        &format!(
            "\"detail\":\"circle / twist / helix strain states from end moments\",\
             \"arc_radius_dev\":{arc_dev:.3e},\"twist_straightness\":{straight_dev:.3e},\
             \"twist_rate_dev\":{twist_rate_dev:.3e},\"helix_kappa_dev\":{kdev:.3e}"
        ),
    );
}

// ------------------------------------------------------------------ str-004

/// N = 0 section state at prescribed curvature (1D Newton on ε₀).
fn balanced_moment(section: &mut fs_solid::Section, kappa: f64, commit: bool) -> f64 {
    let mut e0 = 0.0f64;
    for _ in 0..80 {
        let r = section.respond(e0, kappa);
        if r.n.abs() < 1e-3 * (1.0 + r.tangent[0][0].abs() * 1e-6) {
            break;
        }
        e0 -= 0.8 * r.n / r.tangent[0][0].max(1e3);
    }
    let m = section.respond(e0, kappa).m;
    if commit {
        section.commit(e0, kappa);
    }
    m
}

#[test]
fn str_004_rc_hysteresis() {
    let run = || {
        let mut section = rc_section(0.5, 0.3, 20, 1e-3);
        let amplitudes = [0.006f64, 0.012, 0.024];
        let mut dissipated = Vec::new();
        let mut peak_m = 0.0f64;
        let mut all_m = Vec::new();
        for &a in &amplitudes {
            // One full cycle 0 → +a → −a → +a-ish closing at 0.
            let mut path = Vec::new();
            let steps = 24;
            for k in 0..=steps {
                path.push(a * f64::from(k) / f64::from(steps));
            }
            for k in (-steps..=steps).rev() {
                path.push(a * f64::from(k) / f64::from(steps));
            }
            for k in -steps..=0 {
                path.push(a * f64::from(k) / f64::from(steps));
            }
            let mut work = 0.0;
            let mut prev: Option<(f64, f64)> = None;
            for &kap in &path {
                let m = balanced_moment(&mut section, kap, true);
                peak_m = peak_m.max(m.abs());
                all_m.push(m);
                if let Some((kp, mp)) = prev {
                    work += f64::midpoint(m, mp) * (kap - kp);
                }
                prev = Some((kap, m));
            }
            dissipated.push(work);
        }
        (dissipated, peak_m, all_m)
    };
    let (dissipated, peak_m, m_a) = run();
    let (_, _, m_b) = run();
    let deterministic = m_a
        .iter()
        .zip(&m_b)
        .all(|(x, y)| x.to_bits() == y.to_bits());
    // Hand capacity: steel couple As·fy·(arm) + concrete block share.
    let steel_couple = 1e-3 * 450e6 * 0.4; // two layers at ±0.2 m
    let growing = dissipated.windows(2).all(|w| w[1] > w[0]);
    let positive = dissipated.iter().all(|&d| d > 0.0);
    let capacity_ok = peak_m > 0.7 * steel_couple && peak_m < 3.0 * steel_couple;
    let pass = positive && growing && capacity_ok && deterministic;
    let mut rows = String::new();
    for (a, d) in [0.006f64, 0.012, 0.024].iter().zip(&dissipated) {
        let _ = write!(rows, "{{\"amplitude\":{a},\"dissipated\":{d:.4e}}},");
    }
    verdict(
        "str-004",
        pass,
        &format!(
            "\"detail\":\"RC section cyclic moment-curvature: dissipation grows, capacity band, \
             bitwise determinism\",\"cycles\":[{}],\"peak_m\":{peak_m:.4e},\
             \"steel_couple\":{steel_couple:.4e},\"deterministic\":{deterministic}",
            rows.trim_end_matches(',')
        ),
    );
}

// ------------------------------------------------------------------ str-005

#[test]
fn str_005_batched_consistency_and_throughput() {
    let m = 4096usize;
    let sections: Vec<_> = (0..m).map(|_| rc_section(0.5, 0.3, 16, 1e-3)).collect();
    #[allow(clippy::cast_precision_loss)]
    let strains: Vec<(f64, f64)> = (0..m)
        .map(|i| {
            let t = i as f64 / m as f64;
            (1e-4 * (t - 0.5), 0.01 * (std::f64::consts::TAU * t).sin())
        })
        .collect();
    let rhs: Vec<(f64, f64)> = (0..m)
        .map(|i| (1.0, if i % 2 == 0 { 0.5 } else { -0.5 }))
        .collect();
    let t0 = Instant::now();
    let (resp, sol) = update_sections_batched(&sections, &strains, &rhs);
    let elapsed = t0.elapsed().as_secs_f64();
    #[allow(clippy::cast_precision_loss)]
    let throughput = m as f64 / elapsed;
    // Scalar-path consistency (response identical; solve to solver tol).
    let mut worst = 0.0f64;
    for i in 0..m {
        let r = sections[i].respond(strains[i].0, strains[i].1);
        worst = worst
            .max((r.n - resp[i].n).abs())
            .max((r.m - resp[i].m).abs());
        // Residual of the 2×2 solve.
        let t = r.tangent;
        let res0 = t[0][0] * sol[i].0 + t[0][1] * sol[i].1 - rhs[i].0;
        let res1 = t[1][0] * sol[i].0 + t[1][1] * sol[i].1 - rhs[i].1;
        let scale = t[0][0].abs().max(t[1][1].abs()).max(1.0);
        worst = worst.max(res0.abs() / scale).max(res1.abs() / scale);
    }
    let pass = worst < 1e-8 && throughput > 1e4;
    verdict(
        "str-005",
        pass,
        &format!(
            "\"detail\":\"batched section updates: scalar-consistent, throughput ledgered\",\
             \"sections\":{m},\"elapsed_s\":{elapsed:.4},\"sections_per_s\":{throughput:.3e},\
             \"worst_dev\":{worst:.3e}"
        ),
    );
}

// ------------------------------------------------------------------ str-006

#[test]
fn str_006_forcebased_pushover_and_resume() {
    let make = || rc_section(0.5, 0.3, 16, 1e-3);
    let mut elem = ForceBasedElement::new(3.0, &make);
    // Pushover: ramp tip shear well past section yield.
    let vy_est = 1e-3 * 450e6 * 0.4 / 3.0; // steel couple / length
    let mut rows = String::new();
    let mut curve: Vec<(f64, f64)> = Vec::new();
    let steps = 14i32;
    let mut checkpoint: Option<(ForceBasedElement, i32)> = None;
    for k in 1..=steps {
        let v = 1.4 * vy_est * f64::from(k) / f64::from(steps);
        let d = elem.tip_deflection_under_shear(v).expect("pushover step");
        curve.push((v, d));
        let (_, kap) = elem.base_committed();
        let _ = write!(
            rows,
            "{{\"v\":{v:.4e},\"d\":{d:.5e},\"base_kappa\":{kap:.4e}}},"
        );
        if k == steps / 2 {
            checkpoint = Some((elem.clone(), k));
        }
    }
    // Post-yield stiffness loss: secant of the last step vs the first.
    let k_first = curve[0].0 / curve[0].1;
    let last = curve.len() - 1;
    let (dv, dd) = (
        curve[last].0 - curve[last - 1].0,
        curve[last].1 - curve[last - 1].1,
    );
    let k_last = dv / dd;
    let softened = k_last < 0.5 * k_first;
    // Reversal dissipation: unload to −0.8 v_max and back.
    let vmax = curve[last].0;
    let mut loop_work = 0.0;
    let mut prev = (vmax, curve[last].1);
    for &v in &[
        0.4 * vmax,
        -0.4 * vmax,
        -0.8 * vmax,
        -0.4 * vmax,
        0.4 * vmax,
        vmax,
    ] {
        let d = elem.tip_deflection_under_shear(v).expect("cycle step");
        loop_work += f64::midpoint(v, prev.0) * (d - prev.1);
        prev = (v, d);
    }
    // G4: resume from the checkpoint replays the remaining pushover
    // bitwise.
    let (mut resumed, from) = checkpoint.expect("checkpoint stored");
    let mut resumed_match = true;
    for k in (from + 1)..=steps {
        let v = 1.4 * vy_est * f64::from(k) / f64::from(steps);
        let d = resumed.tip_deflection_under_shear(v).expect("resume step");
        #[allow(clippy::cast_sign_loss)]
        if d.to_bits() != curve[k as usize - 1].1.to_bits() {
            resumed_match = false;
        }
    }
    let pass = softened && loop_work > 0.0 && resumed_match;
    verdict(
        "str-006",
        pass,
        &format!(
            "\"detail\":\"force-based pushover: softening, hysteretic loop work, G4 resume\",\
             \"k_first\":{k_first:.4e},\"k_last\":{k_last:.4e},\"loop_work\":{loop_work:.4e},\
             \"resume_bitwise\":{resumed_match},\"curve\":[{}]",
            rows.trim_end_matches(',')
        ),
    );
}
