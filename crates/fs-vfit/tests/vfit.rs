//! fs-vfit conformance battery: known-answer identification, front-end
//! cross-check, passivity certification/repair, discretization parity,
//! property tests (idempotence, scaling covariance), determinism, and
//! the mutation-visibility checks the bead demands.

use fs_math::c64::C64;
use fs_vfit::discretize::{bilinear, bilinear_state_space};
use fs_vfit::loewner::{cross_check, loewner_fit};
use fs_vfit::passivity::{CertificateClass, check_passivity, repair_passivity};
use fs_vfit::vf::fit_auto_order;
use fs_vfit::{FitOptions, PoleTerm, RationalModel, WeightPreset, vector_fit};

/// The reference 6-pole passive-style impedance: three weakly damped
/// conjugate pairs with mostly-real residues plus small direct terms.
fn six_pole_model() -> RationalModel {
    RationalModel {
        terms: vec![
            PoleTerm::Pair {
                pole: C64::new(-12.0, 700.0),
                residue: C64::new(90.0, 4.0),
            },
            PoleTerm::Pair {
                pole: C64::new(-25.0, 2100.0),
                residue: C64::new(150.0, -9.0),
            },
            PoleTerm::Pair {
                pole: C64::new(-60.0, 5200.0),
                residue: C64::new(240.0, 15.0),
            },
        ],
        d: 0.02,
        e: 1.5e-6,
    }
}

fn log_grid(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| {
            let t = i as f64 / (n - 1) as f64;
            lo * fs_math::det::exp(t * fs_math::det::ln(hi / lo))
        })
        .collect()
}

fn sample(model: &RationalModel, omega: &[f64]) -> Vec<C64> {
    omega.iter().map(|&w| model.eval_iw(w)).collect()
}

#[test]
fn realization_parity_two_routes() {
    // RationalModel::eval (partial fractions) vs StateSpace::eval
    // (complex LU on the realization) — independent arithmetic routes.
    let model = six_pole_model();
    let ss = model.state_space();
    for &w in &[80.0, 700.0, 1234.5, 5200.0, 9000.0] {
        let s = C64::new(0.0, w);
        let direct = model.eval(s);
        let via_ss = ss.eval(s).expect("nonsingular");
        let denom = direct.abs().max(1e-300);
        assert!(
            (direct - via_ss).abs() / denom < 1.0e-12,
            "realization mismatch at w={w}: {direct:?} vs {via_ss:?}"
        );
    }
    // Off-axis too (general s, not just the imaginary axis).
    let s = C64::new(-30.0, 400.0);
    let direct = model.eval(s);
    let via_ss = ss.eval(s).expect("nonsingular");
    assert!((direct - via_ss).abs() / direct.abs() < 1.0e-12);
}

#[test]
fn known_answer_six_pole_recovery() {
    let truth = six_pole_model();
    let omega = log_grid(50.0, 2.0e4, 400);
    let h = sample(&truth, &omega);
    let outcome = vector_fit(&omega, &h, &FitOptions::new(6)).expect("fit");
    // Response recovery to 1e-9 RELATIVE everywhere on the grid.
    let mut worst = 0.0f64;
    for (&w, hv) in omega.iter().zip(&h) {
        let rel = (outcome.model.eval_iw(w) - *hv).abs() / hv.abs();
        worst = worst.max(rel);
    }
    assert!(
        worst < 1.0e-9,
        "six-pole known answer only recovered to {worst:.3e}"
    );
    // Pole recovery: every truth pole has a fitted pole within 1e-6
    // relative.
    let fitted = outcome.model.poles_expanded();
    for p in truth.poles_expanded() {
        let best = fitted
            .iter()
            .map(|q| (*q - p).abs() / p.abs())
            .fold(f64::INFINITY, f64::min);
        assert!(best < 1.0e-6, "pole {p:?} recovered only to {best:.3e}");
    }
    assert!(outcome.model.is_stable());
    println!(
        "{{\"suite\":\"fs-vfit\",\"case\":\"known-answer-6pole\",\"worst_rel\":{worst:.3e},\"verdict\":\"pass\"}}"
    );
}

#[test]
fn noisy_data_fit_within_noise_floor() {
    let truth = six_pole_model();
    let omega = log_grid(50.0, 2.0e4, 400);
    let mut h = sample(&truth, &omega);
    // Deterministic pseudo-noise at 1e-4 relative (no RNG: phase from
    // a fixed irrational stride).
    let noise_rel = 1.0e-4;
    for (i, v) in h.iter_mut().enumerate() {
        let ph = 2.399963229728653 * i as f64; // golden-angle stride
        let mag = v.abs() * noise_rel;
        *v = *v + C64::new(mag * fs_math::det::cos(ph), mag * fs_math::det::sin(ph));
    }
    let outcome = vector_fit(&omega, &h, &FitOptions::new(6)).expect("fit");
    // The fit should sit AT the noise floor: relative misfit within a
    // small factor of the injected noise, and NOT well below it
    // against the noisy data (that would be overfitting the noise —
    // impossible at order 6, which is the point of the fixture).
    let mut worst = 0.0f64;
    for (&w, hv) in omega.iter().zip(&h) {
        let rel = (outcome.model.eval_iw(w) - *hv).abs() / hv.abs();
        worst = worst.max(rel);
    }
    assert!(
        worst < 8.0 * noise_rel,
        "noisy fit misfit {worst:.3e} far above the {noise_rel:.0e} floor"
    );
    // And against the clean truth it must be at least as good.
    let mut worst_truth = 0.0f64;
    for &w in &omega {
        let rel = (outcome.model.eval_iw(w) - truth.eval_iw(w)).abs() / truth.eval_iw(w).abs();
        worst_truth = worst_truth.max(rel);
    }
    assert!(worst_truth < 8.0 * noise_rel);
    println!(
        "{{\"suite\":\"fs-vfit\",\"case\":\"noise-floor\",\"misfit\":{worst:.3e},\"vs_truth\":{worst_truth:.3e},\"verdict\":\"pass\"}}"
    );
}

#[test]
fn fit_of_a_fit_is_idempotent() {
    let truth = six_pole_model();
    let omega = log_grid(50.0, 2.0e4, 300);
    let h = sample(&truth, &omega);
    let first = vector_fit(&omega, &h, &FitOptions::new(6)).expect("first");
    let h2 = sample(&first.model, &omega);
    let second = vector_fit(&omega, &h2, &FitOptions::new(6)).expect("second");
    for &w in &omega {
        let a = first.model.eval_iw(w);
        let b = second.model.eval_iw(w);
        assert!(
            (a - b).abs() / a.abs() < 1.0e-8,
            "fit-of-a-fit drifted at w={w}"
        );
    }
}

#[test]
fn scaling_covariance() {
    // H -> a*H must scale residues/d/e by a and keep poles fixed.
    let truth = six_pole_model();
    let omega = log_grid(50.0, 2.0e4, 300);
    let h = sample(&truth, &omega);
    let scale = 37.5;
    let h_scaled: Vec<C64> = h.iter().map(|v| v.scale(scale)).collect();
    let base = vector_fit(&omega, &h, &FitOptions::new(6)).expect("base");
    let scaled = vector_fit(&omega, &h_scaled, &FitOptions::new(6)).expect("scaled");
    for &w in &omega {
        let a = base.model.eval_iw(w).scale(scale);
        let b = scaled.model.eval_iw(w);
        assert!(
            (a - b).abs() / a.abs() < 1.0e-9,
            "scaling covariance broken at w={w}"
        );
    }
}

#[test]
fn determinism_bitwise_repeat() {
    let truth = six_pole_model();
    let omega = log_grid(50.0, 2.0e4, 300);
    let h = sample(&truth, &omega);
    let a = vector_fit(&omega, &h, &FitOptions::new(6)).expect("a");
    let b = vector_fit(&omega, &h, &FitOptions::new(6)).expect("b");
    assert_eq!(a.model, b.model, "identical inputs must refit bitwise");
    assert_eq!(a.model.d.to_bits(), b.model.d.to_bits());
    assert_eq!(a.model.e.to_bits(), b.model.e.to_bits());
}

#[test]
fn loewner_cross_check_clean_and_aliased() {
    let truth = six_pole_model();
    let omega = log_grid(50.0, 2.0e4, 200);
    let h = sample(&truth, &omega);
    let opts = FitOptions::new(6);
    let vf_fit = vector_fit(&omega, &h, &opts).expect("vf");
    let (lw_fit, ratios) = loewner_fit(&omega, &h, 6, 1.0e-9, &opts).expect("loewner");
    // Rank reveal: a clean 6-pole system keeps six live directions.
    assert!(
        ratios[5] > 1.0e-6,
        "rank reveal collapsed early: {ratios:?}"
    );
    let agree = cross_check(&omega, &h, &vf_fit.model, &lw_fit.model);
    // Authored agreement gate: the iterated direct-term stripping
    // leaves ~1e-6-relative pole residue on clean improper data, which
    // the sharpest resonance (Q ~ 87) amplifies into the 1e-4 response
    // class; 1e-3 is an order of headroom while staying 30x below the
    // aliased diagnostic asserted after.
    assert!(
        agree.worst_response_mismatch < 1.0e-3,
        "front ends disagree on clean data: {agree:?}"
    );
    // ALIASED data: sample the same system on a grid that undersamples
    // the top resonance region entirely (only 12 points, none near the
    // 5200 rad/s pair). The front ends may each do their best; the
    // cross-check must DIAGNOSE the discrepancy (logged comparison)
    // rather than agree.
    let omega_bad = log_grid(50.0, 900.0, 12);
    let h_bad = sample(&truth, &omega_bad);
    let vf_bad = vector_fit(&omega_bad, &h_bad, &opts).expect("vf aliased");
    let lw_bad = loewner_fit(&omega_bad, &h_bad, 6, 1.0e-9, &opts);
    let diag = match lw_bad {
        Ok((lw_model, _)) => {
            let probe_grid = log_grid(50.0, 2.0e4, 100);
            let probe_h = sample(&truth, &probe_grid);
            cross_check(&probe_grid, &probe_h, &vf_bad.model, &lw_model.model)
                .worst_response_mismatch
        }
        // A typed refusal on degenerate data is ALSO a diagnostic.
        Err(_) => f64::INFINITY,
    };
    assert!(
        diag > 3.0e-2,
        "aliased data should NOT cross-check clean (diag {diag:.3e})"
    );
    println!(
        "{{\"suite\":\"fs-vfit\",\"case\":\"cross-check\",\"clean\":{:.3e},\"aliased_diag\":{diag:.3e},\"verdict\":\"pass\"}}",
        agree.worst_response_mismatch
    );
}

/// A deliberately ACTIVE near-miss: start from a passive-style model
/// and crank one pair's imaginary residue until `Re H` dips negative
/// near the resonance. The fixture asserts its own activity so the
/// repair demonstration cannot go vacuous.
fn near_miss_active_model() -> RationalModel {
    RationalModel {
        terms: vec![
            PoleTerm::Pair {
                pole: C64::new(-12.0, 700.0),
                residue: C64::new(90.0, 4.0),
            },
            PoleTerm::Pair {
                pole: C64::new(-25.0, 2100.0),
                // sigma = 20 makes min Re H ~= -0.017 against d = 0.05
                // (verified numerically at authoring): a NEAR miss, so
                // the repair is a perturbation, not a rewrite.
                residue: C64::new(60.0, 20.0),
            },
        ],
        d: 0.05,
        e: 0.0,
    }
}

#[test]
fn passivity_detects_and_repairs_active_fit() {
    let band = (50.0, 2.0e4);
    // A genuinely passive model certifies green with the EXACT class.
    let good = six_pole_model();
    let good_report = check_passivity(&good, band).expect("check");
    assert!(good_report.passive, "reference model must certify passive");
    assert_eq!(good_report.class, CertificateClass::HamiltonianExact);
    // The near-miss fixture is REALLY active (fixture validity)...
    let bad = near_miss_active_model();
    let bad_report = check_passivity(&bad, band).expect("check");
    assert!(
        !bad_report.passive,
        "fixture must be active; worst Re H = {:?}",
        bad_report.worst
    );
    // ...the Hamiltonian arm sees crossings (MUTATION VISIBILITY: with
    // repair disabled, this is the failing signal)...
    assert!(
        !bad_report.crossings.is_empty(),
        "hamiltonian test must flag the active fixture"
    );
    // ...and the QP repair makes it green with the exact certificate.
    let (repaired, report) = repair_passivity(&bad, band).expect("repair");
    assert!(report.certificate.passive, "repair must reach passivity");
    assert_eq!(report.certificate.class, CertificateClass::HamiltonianExact);
    assert!(report.certificate.crossings.is_empty());
    // The repair is a PERTURBATION, not a rewrite: poles untouched,
    // moderate relative residue movement.
    assert_eq!(
        repaired.poles_expanded(),
        bad.poles_expanded(),
        "repair must not move poles"
    );
    assert!(
        report.relative_perturbation < 1.0,
        "repair blew up the residues: {:.3}",
        report.relative_perturbation
    );
    // KKT stationarity of the final QP: near machine zero relative to
    // the residue scale.
    assert!(
        report.kkt_residual < 1.0e-6,
        "KKT residual {:.3e}",
        report.kkt_residual
    );
    println!(
        "{{\"suite\":\"fs-vfit\",\"case\":\"passivity-repair\",\"rounds\":{},\"rel_pert\":{:.4},\"kkt\":{:.2e},\"verdict\":\"pass\"}}",
        report.rounds, report.relative_perturbation, report.kkt_residual
    );
}

#[test]
fn discretization_matches_continuous_in_band() {
    let model = six_pole_model();
    let fs_hz = 48000.0;
    let t_s = 1.0 / fs_hz;
    let prewarp = 5200.0; // top resonance: where match matters most
    let filt = bilinear(&model, t_s, prewarp).expect("bilinear");
    assert!(filt.is_stable(), "discretized sections must be stable");
    // At the prewarp frequency the match is EXACT (that is what
    // prewarping means); elsewhere in band it is within tolerance.
    let at_pw = filt.eval(prewarp);
    let cont_pw = model.eval_iw(prewarp);
    assert!(
        (at_pw - cont_pw).abs() / cont_pw.abs() < 1.0e-9,
        "prewarp point must match exactly"
    );
    let band = log_grid(50.0, 8000.0, 200);
    let mut worst = 0.0f64;
    for &w in &band {
        let rel = (filt.eval(w) - model.eval_iw(w)).abs() / model.eval_iw(w).abs();
        worst = worst.max(rel);
    }
    // Authored tolerance, derived not guessed: bilinear frequency warp
    // at omega is ~(omega*T)^2/12 relative (0.08% at the 2100 rad/s
    // resonance at fs = 48 kHz with the prewarp pinned at 5200), and a
    // resonance amplifies a frequency shift by Q = omega0/(2|Re p|)
    // ~= 42, so the worst in-band deviation is EXPECTED at the ~5%
    // scale near the sharpest off-prewarp resonance. Measured 5.1e-2
    // at authoring time; 8% keeps headroom without hiding a real
    // regression (a lost prewarp shows up 10x larger — asserted below).
    assert!(
        worst < 0.08,
        "in-band bilinear deviation {worst:.3e} above the authored 8%"
    );
    // MUTATION VISIBILITY (dropped prewarp): the unwarped map is
    // measurably worse AT THE TOP RESONANCE than the prewarped one.
    let nowarp = bilinear(&model, t_s, 0.0).expect("bilinear nowarp");
    let warped_err = (filt.eval(prewarp) - cont_pw).abs() / cont_pw.abs();
    let nowarp_err = (nowarp.eval(prewarp) - cont_pw).abs() / cont_pw.abs();
    assert!(
        nowarp_err > 10.0 * warped_err.max(1.0e-12),
        "dropping prewarp must be visible at the resonance ({nowarp_err:.3e} vs {warped_err:.3e})"
    );
    // State-space route agrees with the section route (independent
    // realizations of the same map; e handled as the leftover section).
    let dss = bilinear_state_space(&model, t_s, prewarp).expect("dss");
    for &w in &[100.0, 700.0, 2100.0, 5200.0] {
        let via_ss = dss.eval(w).expect("eval") + section_e_term(&model, t_s, prewarp, w);
        let via_sections = filt.eval(w);
        assert!(
            (via_ss - via_sections).abs() / via_sections.abs() < 1.0e-9,
            "state-space vs sections mismatch at w={w}"
        );
    }
    // Quantization sensitivity: f32 coefficients stay within a
    // reported (loose) envelope in band — and the report is nonzero
    // (the probe actually does something).
    let mut worst_q = 0.0f64;
    for &w in &band {
        let dq = (filt.eval_f32_quantized(w) - filt.eval(w)).abs() / filt.eval(w).abs();
        worst_q = worst_q.max(dq);
    }
    assert!(worst_q > 0.0, "quantization probe must not be vacuous");
    assert!(
        worst_q < 0.05,
        "f32 quantization deviation {worst_q:.3e} above 5%"
    );
    println!(
        "{{\"suite\":\"fs-vfit\",\"case\":\"discretize\",\"in_band_worst\":{worst:.3e},\"f32_worst\":{worst_q:.3e},\"verdict\":\"pass\"}}"
    );
}

fn section_e_term(model: &RationalModel, t_s: f64, omega_pw: f64, w: f64) -> C64 {
    // The e-section the state-space form leaves to the caller.
    if model.e == 0.0 {
        return C64::ZERO;
    }
    let k = if omega_pw > 0.0 {
        omega_pw / fs_math::det::tan(omega_pw * t_s / 2.0)
    } else {
        2.0 / t_s
    };
    let zi = C64::new(fs_math::det::cos(w * t_s), -fs_math::det::sin(w * t_s));
    let num = C64::from_re(model.e * k) - zi.scale(model.e * k);
    let den = C64::ONE + zi;
    num * den.recip()
}

#[test]
fn auto_order_selection_plateaus_and_refuses_overfit() {
    let truth = six_pole_model();
    let omega = log_grid(50.0, 2.0e4, 300);
    let mut h = sample(&truth, &omega);
    let noise_rel = 1.0e-5;
    for (i, v) in h.iter_mut().enumerate() {
        let ph = 2.399963229728653 * i as f64;
        let mag = v.abs() * noise_rel;
        *v = *v + C64::new(mag * fs_math::det::cos(ph), mag * fs_math::det::sin(ph));
    }
    let base = FitOptions::new(2);
    // plateau_ratio 0.2: an order that fails to improve the weighted
    // error by at least 20% ends the ascent; the floor is 3x the
    // injected noise so the stop is robust to the noise realization.
    let (selected, curve) = fit_auto_order(
        &omega,
        &h,
        &[2, 4, 6, 8, 10, 12],
        &base,
        0.2,
        3.0 * noise_rel,
    )
    .expect("auto order");
    // The curve must show the plateau at the true order: order 6 error
    // is orders of magnitude below order 4.
    let err_at = |o: usize| {
        curve
            .iter()
            .find(|(ord, _)| *ord == o)
            .map(|(_, e)| *e)
            .expect("order in curve")
    };
    assert!(err_at(4) > 100.0 * err_at(6), "no plateau step at order 6");
    // Selection stopped AT the plateau (overfit refusal): no order
    // beyond 8 was even fitted once the floor was hit.
    assert!(
        selected.model.order() <= 8,
        "selected over-rich order {}",
        selected.model.order()
    );
    assert!(
        curve.last().expect("nonempty").0 <= 8,
        "kept fitting past the noise floor: {curve:?}"
    );
    println!(
        "{{\"suite\":\"fs-vfit\",\"case\":\"auto-order\",\"curve\":{curve:?},\"selected\":{},\"verdict\":\"pass\"}}",
        selected.model.order()
    );
}

#[test]
fn weight_preset_changes_the_fit_and_is_logged() {
    // MUTATION VISIBILITY (weight preset swap): under-resolved fits
    // differ measurably between uniform and inverse-magnitude weights,
    // and the report names the preset.
    let truth = six_pole_model();
    let omega = log_grid(50.0, 2.0e4, 300);
    let h = sample(&truth, &omega);
    // Order 4 (under-resolved on purpose: preset choice must matter).
    let uni = vector_fit(
        &omega,
        &h,
        &FitOptions {
            weights: WeightPreset::Uniform,
            ..FitOptions::new(4)
        },
    )
    .expect("uniform");
    let inv = vector_fit(&omega, &h, &FitOptions::new(4)).expect("inverse");
    assert_eq!(uni.report.weights, "uniform");
    assert_eq!(inv.report.weights, "inverse-magnitude");
    // At the antiresonance (response minimum) the inverse-magnitude
    // fit must be relatively better than the uniform fit.
    let w_min = *omega
        .iter()
        .min_by(|a, b| {
            truth
                .eval_iw(**a)
                .abs()
                .total_cmp(&truth.eval_iw(**b).abs())
        })
        .expect("min");
    let rel = |m: &RationalModel| {
        (m.eval_iw(w_min) - truth.eval_iw(w_min)).abs() / truth.eval_iw(w_min).abs()
    };
    assert!(
        rel(&inv.model) < rel(&uni.model),
        "inverse-magnitude weighting must win at the antiresonance ({:.3e} vs {:.3e})",
        rel(&inv.model),
        rel(&uni.model)
    );
}

#[test]
fn nyquist_refusal_on_both_discretization_routes() {
    // Review finding: the Tustin route must refuse a beyond-Nyquist
    // resonance just like the section route (no silent aliasing).
    let model = six_pole_model();
    let t_s = 1.0 / 1000.0; // nyquist = pi*1000 ~= 3141 rad/s < 5200
    assert!(bilinear(&model, t_s, 0.0).is_err());
    assert!(bilinear_state_space(&model, t_s, 0.0).is_err());
}

#[test]
fn refusals_are_typed() {
    let omega = log_grid(50.0, 2.0e4, 40);
    let h: Vec<C64> = omega.iter().map(|&w| C64::new(1.0, w * 1.0e-4)).collect();
    // Order zero.
    assert!(vector_fit(&omega, &h, &FitOptions::new(0)).is_err());
    // Bad sample (negative frequency).
    let mut bad_omega = omega.clone();
    bad_omega[3] = -1.0;
    assert!(vector_fit(&bad_omega, &h, &FitOptions::new(4)).is_err());
    // Too few samples.
    assert!(vector_fit(&omega[..3], &h[..3], &FitOptions::new(8)).is_err());
}
