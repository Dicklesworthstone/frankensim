//! fs-modalid conformance battery: known-answer identification with
//! noise sweeps inside the reported intervals, close-mode cross-check,
//! stabilization automation, MAC gates, window-correction mutation
//! visibility, RFP conditioning demonstration, CSV ingest, and the
//! SNR refusal path.

use fs_math::c64::C64;
use fs_modalid::{
    FrfData, IdentifyOptions, ModalIdError, correct_exponential_window, estimate_snr, identify,
    mac, mac_pairing, rfp_fit,
};
use fs_vfit::vf::FitOptions;

const TWO_PI: f64 = 2.0 * core::f64::consts::PI;

/// Synthetic modal truth: 10 modes with shapes over 4 channels.
struct Truth {
    freqs_hz: Vec<f64>,
    zetas: Vec<f64>,
    /// shapes[mode][channel] (real, sign-varying).
    shapes: Vec<Vec<f64>>,
}

fn truth10() -> Truth {
    let freqs_hz: Vec<f64> = (0..10)
        .map(|i| 97.0f64.mul_add(f64::from(i), 120.0))
        .collect();
    let zetas: Vec<f64> = (0..10)
        .map(|i| 0.001f64.mul_add(f64::from(i), 0.008))
        .collect();
    let shapes: Vec<Vec<f64>> = (0i32..10)
        .map(|k| {
            (0i32..4)
                .map(|c| fs_math::det::sin(f64::from(k + 1) * f64::from(c + 1) * 0.7))
                .collect()
        })
        .collect();
    Truth {
        freqs_hz,
        zetas,
        shapes,
    }
}

/// Engineering-convention (e^{-iwt}) receptance FRF of the truth:
/// H(w) = sum phi_k / (wn^2 - w^2 - 2 i zeta wn w) per channel — the
/// standard measured-FRF form.
fn synth_frf(truth: &Truth, omega: &[f64], channel: usize) -> Vec<C64> {
    omega
        .iter()
        .map(|&w| {
            let mut acc = C64::ZERO;
            for k in 0..truth.freqs_hz.len() {
                let wn = TWO_PI * truth.freqs_hz[k];
                let den = C64::new(wn * wn - w * w, -2.0 * truth.zetas[k] * wn * w);
                acc = acc + C64::from_re(truth.shapes[k][channel]) * den.recip();
            }
            acc
        })
        .collect()
}

fn grid() -> Vec<f64> {
    (0i32..1600)
        .map(|i| TWO_PI * f64::from(i).mul_add(0.68, 60.0))
        .collect()
}

/// Deterministic pseudo-noise (golden-angle phase stride, no RNG).
fn add_noise(h: &mut [C64], rel: f64, seed: usize) {
    for (i, v) in h.iter_mut().enumerate() {
        let ph = 2.399963229728653 * (i + 1000 * seed) as f64;
        let mag = v.abs() * rel;
        *v = *v + C64::new(mag * fs_math::det::cos(ph), mag * fs_math::det::sin(ph));
    }
}

fn data_with_noise(rel: f64, seed: usize) -> FrfData {
    let truth = truth10();
    let omega = grid();
    let channels: Vec<(Vec<C64>, Option<Vec<f64>>)> = (0..4)
        .map(|c| {
            let mut h = synth_frf(&truth, &omega, c);
            add_noise(&mut h, rel, seed);
            (h, None)
        })
        .collect();
    FrfData::new(omega, channels).expect("data")
}

#[test]
fn known_answer_ten_modes_with_noise_inside_intervals() {
    let truth = truth10();
    for (noise, seed) in [(1.0e-6, 1), (1.0e-4, 2), (1.0e-3, 3)] {
        let data = data_with_noise(noise, seed);
        let id = identify(&data, &IdentifyOptions::default()).expect("identify");
        assert_eq!(
            id.modes.len(),
            10,
            "noise {noise:.0e}: expected 10 accepted modes, got {}",
            id.modes.len()
        );
        for (k, mode) in id.modes.iter().enumerate() {
            let df = (mode.frequency_hz - truth.freqs_hz[k]).abs();
            let dz = (mode.damping_ratio - truth.zetas[k]).abs();
            // Within the reported split-sample interval PLUS an
            // authored floor (the CI can legitimately be ~0 on clean
            // data; the floor is the parameterization granularity
            // measured at authoring: 1e-3 Hz / 1e-5 zeta at 1e-6
            // noise, scaling roughly linearly with noise).
            // Absolute caps (review finding: an uncapped CI-based
            // gate is self-referential — a regression that inflates
            // the reported spread stays green).
            let f_gate = mode.frequency_ci_hz.min(0.5).max(2.0e3 * noise + 1.0e-4);
            let z_gate = mode.damping_ci.min(0.01).max(20.0 * noise + 1.0e-7);
            assert!(
                df <= f_gate,
                "noise {noise:.0e} mode {k}: df {df:.3e} Hz above gate {f_gate:.3e}"
            );
            assert!(
                dz <= z_gate,
                "noise {noise:.0e} mode {k}: dz {dz:.3e} above gate {z_gate:.3e}"
            );
        }
        println!(
            "{{\"suite\":\"fs-modalid\",\"case\":\"known-answer\",\"noise\":{noise:.0e},\"modes\":{},\"verdict\":\"pass\"}}",
            id.modes.len()
        );
    }
}

#[test]
fn shapes_recover_truth_by_mac() {
    let truth = truth10();
    let data = data_with_noise(1.0e-5, 4);
    let id = identify(&data, &IdentifyOptions::default()).expect("identify");
    let truth_shapes: Vec<Vec<C64>> = truth
        .shapes
        .iter()
        .map(|s| s.iter().map(|&v| C64::from_re(v)).collect())
        .collect();
    let id_shapes: Vec<Vec<C64>> = id.modes.iter().map(|m| m.shape.clone()).collect();
    let pairing = mac_pairing(&id_shapes, &truth_shapes, 0.95);
    for (i, j, v) in &pairing {
        assert_eq!(
            Some(*i),
            *j,
            "identified mode {i} paired to {j:?} (MAC {v:.4})"
        );
        assert!(*v > 0.99, "mode {i} MAC {v:.4} below 0.99");
    }
}

#[test]
fn close_mode_pair_resolved_and_cross_checked() {
    // Two modes 0.5% apart: a single low-order identifier smears them;
    // the stabilization ladder at sufficient order resolves both, and
    // the RFP cross-check agrees on the pair.
    let omega: Vec<f64> = (0i32..2400)
        .map(|i| TWO_PI * f64::from(i).mul_add(0.025, 180.0))
        .collect();
    let (f1, f2) = (200.0, 201.0); // 0.5% separation
    let truth = Truth {
        freqs_hz: vec![f1, f2],
        zetas: vec![0.004, 0.004],
        shapes: vec![vec![1.0], vec![0.8]],
    };
    let h_eng = synth_frf(&truth, &omega, 0);
    let data = FrfData::new(omega.clone(), vec![(h_eng, None)]).expect("data");
    let opts = IdentifyOptions {
        min_order: 4,
        max_order: 12,
        order_step: 2,
        ..IdentifyOptions::default()
    };
    let id = identify(&data, &opts).expect("identify");
    let freqs: Vec<f64> = id.modes.iter().map(|m| m.frequency_hz).collect();
    assert_eq!(
        freqs.len(),
        2,
        "close pair must resolve into 2 modes: {freqs:?}"
    );
    assert!((freqs[0] - f1).abs() < 0.05 && (freqs[1] - f2).abs() < 0.05);
    // Single-identifier low order (2 = one pair) CANNOT represent both
    // — the failure mode the ladder exists for.
    let h_lap: Vec<C64> = data.channels[0].h.clone();
    let low = fs_vfit::vector_fit(&data.omega, &h_lap, &FitOptions::new(2)).expect("low");
    // Order 2 holds exactly one pair BY CONSTRUCTION; the meaningful
    // assertion (review finding: `<= 1` was a tautology) is that its
    // single frequency cannot serve for both true modes.
    let low_freqs: Vec<f64> = low
        .model
        .terms
        .iter()
        .filter_map(|t| match t {
            fs_vfit::PoleTerm::Pair { pole, .. } => Some(pole.abs() / TWO_PI),
            fs_vfit::PoleTerm::Real { .. } => None,
        })
        .collect();
    assert_eq!(low_freqs.len(), 1, "order-2 fit holds one pair");
    assert!(
        (low_freqs[0] - f1).abs() > 0.1 || (low_freqs[0] - f2).abs() > 0.1,
        "one pair cannot sit on both modes: {low_freqs:?}"
    );
    // RFP cross-check at order 6 finds both frequencies.
    let rfp = rfp_fit(&data.omega, &h_lap, 6, &FitOptions::new(6)).expect("rfp");
    let mut rfp_freqs: Vec<f64> = rfp
        .model
        .terms
        .iter()
        .filter_map(|t| match t {
            fs_vfit::PoleTerm::Pair { pole, .. } => Some(pole.abs() / TWO_PI),
            fs_vfit::PoleTerm::Real { .. } => None,
        })
        .collect();
    rfp_freqs.sort_by(f64::total_cmp);
    let near = |target: f64| rfp_freqs.iter().any(|&f| (f - target).abs() < 0.1);
    assert!(
        near(f1) && near(f2),
        "RFP cross-check missed the close pair: {rfp_freqs:?}"
    );
    println!(
        "{{\"suite\":\"fs-modalid\",\"case\":\"close-modes\",\"vf\":{freqs:?},\"rfp\":{rfp_freqs:?},\"verdict\":\"pass\"}}"
    );
}

#[test]
fn mac_identity_and_orthogonality() {
    let a: Vec<C64> = vec![C64::new(1.0, 0.2), C64::new(-0.5, 0.1), C64::from_re(0.3)];
    assert!((mac(&a, &a) - 1.0).abs() < 1.0e-14);
    // Orthogonal real shapes.
    let e1: Vec<C64> = vec![C64::ONE, C64::ZERO];
    let e2: Vec<C64> = vec![C64::ZERO, C64::ONE];
    assert!(mac(&e1, &e2) < 1.0e-30);
    // Scale and phase invariance.
    let b: Vec<C64> = a.iter().map(|v| *v * C64::new(0.0, 2.5)).collect();
    assert!((mac(&a, &b) - 1.0).abs() < 1.0e-14);
}

#[test]
fn perturbed_plate_modes_flagged_by_mac() {
    // Simulated-vs-perturbed calibration diff: perturb two of six
    // synthetic plate-like shapes; pairing must keep the four intact
    // modes above 0.99 and flag the perturbed ones below the floor.
    let n_pts = 25i32;
    let shape = |kx: u32, ky: u32| -> Vec<C64> {
        (0i32..n_pts)
            .map(|p| {
                let (x, y) = (f64::from(p % 5) / 4.0, f64::from(p / 5) / 4.0);
                C64::from_re(
                    fs_math::det::sin(f64::from(kx) * core::f64::consts::PI * x)
                        * fs_math::det::sin(f64::from(ky) * core::f64::consts::PI * y),
                )
            })
            .collect()
    };
    let sim: Vec<Vec<C64>> = vec![
        shape(1, 1),
        shape(2, 1),
        shape(1, 2),
        shape(2, 2),
        shape(3, 1),
        shape(1, 3),
    ];
    let mut meas = sim.clone();
    // Perturb modes 3 and 5 into different shape families.
    meas[3] = shape(3, 2);
    meas[5] = shape(2, 3);
    let pairing = mac_pairing(&sim, &meas, 0.9);
    for &(i, j, v) in &pairing {
        if i == 3 || i == 5 {
            assert!(
                j.is_none() || v < 0.9,
                "perturbed mode {i} not flagged (MAC {v:.4})"
            );
        } else {
            assert_eq!(j, Some(i), "intact mode {i} mispaired");
            assert!(v > 0.99);
        }
    }
}

#[test]
fn window_correction_mutation_visible() {
    // Synthesize an exponentially WINDOWED measurement: every decay
    // rate gains 1/tau. With the correction the identified zeta
    // matches truth; with the correction DISABLED (the mutation) the
    // bias exceeds the acceptance gate — visibly.
    let tau = 0.4;
    let truth = truth10();
    let omega = grid();
    // Windowed FRF: poles shift re -> re - 1/tau, i.e. zeta_measured
    // = zeta + 1/(tau wn) to first order. Build it exactly by pole
    // shifting.
    let windowed: Vec<C64> = omega
        .iter()
        .map(|&w| {
            let mut acc = C64::ZERO;
            for k in 0..truth.freqs_hz.len() {
                let wn = TWO_PI * truth.freqs_hz[k];
                let zeta = truth.zetas[k];
                let sigma = zeta * wn + 1.0 / tau;
                let wd = wn * fs_math::det::sqrt(1.0 - zeta * zeta);
                // Pair with shifted decay, same wd (engineering conv).
                let p = C64::new(-sigma, wd);
                // Real impulse response needs the residue AND its
                // conjugate on conjugate poles (review finding: a
                // -i factor applied outside the pair put -i*r at both
                // poles and broke conjugate symmetry by ~7%).
                let rr = C64::new(0.0, -truth.shapes[k][0] / (2.0 * wd));
                let s = C64::new(0.0, w);
                let val = rr * (s - p).recip() + rr.conj() * (s - p.conj()).recip();
                acc = acc + val.conj(); // to engineering convention
            }
            acc
        })
        .collect();
    let data = FrfData::new(omega, vec![(windowed, None)]).expect("data");
    let corrected = identify(
        &data,
        &IdentifyOptions {
            window_tau: Some(tau),
            ..IdentifyOptions::default()
        },
    )
    .expect("corrected");
    let uncorrected = identify(&data, &IdentifyOptions::default()).expect("uncorrected");
    assert_eq!(corrected.modes.len(), truth.freqs_hz.len());
    let mut worst_corr = 0.0f64;
    let mut worst_raw = 0.0f64;
    for k in 0..truth.freqs_hz.len() {
        let ec = (corrected.modes[k].damping_ratio - truth.zetas[k]).abs();
        worst_corr = worst_corr.max(ec);
        worst_raw = worst_raw.max((uncorrected.modes[k].damping_ratio - truth.zetas[k]).abs());
    }
    // The windowing bias at the lowest mode is 1/(tau*wn) ~ 3.3e-3 —
    // comparable to the zetas themselves.
    assert!(worst_corr < 2.0e-4, "corrected zeta error {worst_corr:.3e}");
    assert!(
        worst_raw > 10.0 * worst_corr.max(1.0e-5),
        "disabled correction must be visible: raw {worst_raw:.3e} vs corrected {worst_corr:.3e}"
    );
    // Raw values are logged alongside.
    assert_eq!(corrected.damping_raw.len(), corrected.modes.len());
    println!(
        "{{\"suite\":\"fs-modalid\",\"case\":\"window-correction\",\"raw_bias\":{worst_raw:.3e},\"corrected\":{worst_corr:.3e},\"verdict\":\"pass\"}}"
    );
}

#[test]
fn snr_refusal_fires() {
    // Noise at the signal scale: the estimate collapses and identify
    // refuses BY NAME.
    let data = data_with_noise(0.8, 9);
    let err = identify(&data, &IdentifyOptions::default());
    assert!(
        matches!(err, Err(ModalIdError::SnrTooLow { .. })),
        "expected SnrTooLow, got {err:?}"
    );
    // Coherence-driven estimate: a coherence table near 1 reports high
    // SNR, near 0.5 reports ~1.
    let omega = grid();
    let truth = truth10();
    let h = synth_frf(&truth, &omega, 0);
    let good = FrfData::new(
        omega.clone(),
        vec![(h.clone(), Some(vec![0.999; omega.len()]))],
    )
    .expect("good");
    assert!(estimate_snr(&good) > 500.0);
    let bad = FrfData::new(omega.clone(), vec![(h, Some(vec![0.5; omega.len()]))]).expect("bad");
    assert!(estimate_snr(&bad) < 2.0);
}

#[test]
fn rfp_orthogonal_basis_beats_naive_powers() {
    // The conditioning demonstration the bead asks for: identify the
    // 10-mode FRF with RFP at order 20 on the Forsythe basis (works),
    // and show the equivalent monomial-basis normal-equations
    // condition estimate is astronomically worse. The monomial
    // comparison builds the same LS matrix with plain powers and
    // measures column-norm spread (a cheap condition lower bound).
    let truth = truth10();
    let omega = grid();
    let h_eng = synth_frf(&truth, &omega, 0);
    let data = FrfData::new(omega.clone(), vec![(h_eng, None)]).expect("data");
    let h = &data.channels[0].h;
    let rfp = rfp_fit(&data.omega, h, 20, &FitOptions::new(20)).expect("rfp");
    let mut freqs: Vec<f64> = rfp
        .model
        .terms
        .iter()
        .filter_map(|t| match t {
            fs_vfit::PoleTerm::Pair { pole, .. } => Some(pole.abs() / TWO_PI),
            fs_vfit::PoleTerm::Real { .. } => None,
        })
        .collect();
    freqs.sort_by(f64::total_cmp);
    for &f_true in &truth.freqs_hz {
        assert!(
            freqs.iter().any(|&f| (f - f_true).abs() < 0.5),
            "RFP missed {f_true} Hz: {freqs:?}"
        );
    }
    // Conditioning: build the GRAM matrix of unit-normalized scaled
    // monomial columns (degree 0..=20 on the band grid) and take its
    // eigenvalue spread — near-collinearity is the classic RFP
    // catastrophe (column NORMS alone look harmless; the executed
    // first attempt measured norms and saw only 6x). The Forsythe
    // basis is orthonormal under the fit weights BY CONSTRUCTION.
    let wmax = data.omega.iter().fold(0.0f64, |a, &v| a.max(v));
    let deg = 20usize;
    let cols: Vec<Vec<f64>> = (0..=deg)
        .map(|d| {
            let raw: Vec<f64> = data
                .omega
                .iter()
                .map(|&w| fs_math::det::powi(w / wmax, i32::try_from(d).expect("small")))
                .collect();
            let norm = fs_math::det::sqrt(raw.iter().map(|v| v * v).sum::<f64>());
            raw.iter().map(|v| v / norm).collect()
        })
        .collect();
    let m = deg + 1;
    let mut gram = vec![0.0f64; m * m];
    for i in 0..m {
        for j in 0..m {
            gram[i * m + j] = cols[i].iter().zip(&cols[j]).map(|(a, b)| a * b).sum();
        }
    }
    let (values, _) = fs_la::eigen::jacobi_eigh(&gram, m);
    let lam_max = values.iter().fold(0.0f64, |a, &v| a.max(v));
    let lam_min = values.iter().fold(f64::INFINITY, |a, &v| a.min(v.abs()));
    let cond = lam_max / lam_min.max(f64::MIN_POSITIVE);
    // Authored: measured ~1e16-class at degree 20 (numerically
    // singular); require > 1e8 — the point where naive-power RFP
    // normal equations lose all accuracy in f64.
    assert!(cond > 1.0e8, "monomial Gram condition only {cond:.3e}");
    // And the OPERATIVE (review-corrected) claim about the Forsythe
    // basis: after Re/Im row stacking with the fit weights, the REAL
    // Gram is near-identity — condition < 10 (the complex Gram on a
    // positive-frequency grid is NOT orthonormal; the parity
    // structure is what carries the conditioning).
    let ws: Vec<f64> = data.omega.iter().map(|&w| w / wmax).collect();
    let wts: Vec<f64> = h
        .iter()
        .map(|v| {
            let mag = v.abs();
            if mag > 0.0 { 1.0 / mag } else { 1.0 }
        })
        .collect();
    let (basis, _) = fs_modalid::forsythe_basis(&ws, &wts, deg);
    let mf = deg + 1;
    let mut gram_f = vec![0.0f64; mf * mf];
    for i in 0..mf {
        for j in 0..mf {
            let mut acc = 0.0;
            for (k, &wt) in wts.iter().enumerate() {
                acc +=
                    wt * wt * (basis[i][k].re * basis[j][k].re + basis[i][k].im * basis[j][k].im);
            }
            gram_f[i * mf + j] = acc;
        }
    }
    let mut dnorm = vec![0.0f64; mf];
    for i in 0..mf {
        dnorm[i] = fs_math::det::sqrt(gram_f[i * mf + i]).max(f64::MIN_POSITIVE);
    }
    for i in 0..mf {
        for j in 0..mf {
            gram_f[i * mf + j] /= dnorm[i] * dnorm[j];
        }
    }
    let (fvals, _) = fs_la::eigen::jacobi_eigh(&gram_f, mf);
    let f_top = fvals.iter().fold(0.0f64, |a, &v| a.max(v));
    let f_bot = fvals.iter().fold(f64::INFINITY, |a, &v| a.min(v.abs()));
    let cond_f = f_top / f_bot.max(f64::MIN_POSITIVE);
    // Authored: measured 32 on this grid (the wt-vs-wt^2 weighting
    // mismatch between orthogonalization and LS costs a factor); gate
    // at 100 — still 14 orders below the monomial catastrophe, which
    // is the operative claim.
    assert!(
        cond_f < 100.0,
        "real-stacked Forsythe Gram condition {cond_f:.2}"
    );
    println!(
        "{{\"suite\":\"fs-modalid\",\"case\":\"rfp-conditioning\",\"monomial_gram_cond\":{cond:.3e},\"verdict\":\"pass\"}}"
    );
}

#[test]
fn csv_ingest_round_trip() {
    let omega = [100.0, 200.0, 300.0];
    let text = "\
# comment line
100.0, 1.0, -0.5, 0.99, 2.0, 0.25, 0.98
200.0, 0.5, -0.25, 0.97, 1.0, 0.5, 0.96
300.0, 0.25, -0.125, 0.95, 0.5, 0.75, 0.94
";
    let data = FrfData::parse_csv(text, 2, true).expect("csv");
    assert_eq!(data.omega, omega);
    assert_eq!(data.channels.len(), 2);
    // Conjugated on ingest (engineering -> Laplace).
    assert!((data.channels[0].h[0] - C64::new(1.0, 0.5)).abs() < 1.0e-15);
    let coh2 = data.channels[1].coherence.as_ref().expect("coh")[2];
    assert!((coh2 - 0.94).abs() < 1.0e-15);
    // Malformed line refuses with the line number.
    let bad = FrfData::parse_csv("100.0, 1.0\n", 2, true);
    assert!(matches!(bad, Err(ModalIdError::Csv { line: 1 })));
}

#[test]
fn published_benchmark_carcagno_guitar_modes() {
    // Published-parameter benchmark: Carcagno et al., JASA 144(6):3533
    // (2018), CC-BY, Table I (verified verbatim from the saved PDF via
    // pdftotext -layout): Brazilian-rosewood guitar low modes
    // F1 = 97 Hz Q1 = 34, F2 = 177 Hz Q2 = 18, F3 = 336 Hz Q3 = 36.
    // The FRF is SYNTHESIZED from the published modal parameters
    // (mobility magnitudes are figures-only across the CC-BY guitar
    // literature — honest label), so this pins the identification
    // pipeline against published values, not against a measured trace.
    let published = [(97.0, 34.0), (177.0, 18.0), (336.0, 36.0)];
    let truth = Truth {
        freqs_hz: published.iter().map(|&(f, _)| f).collect(),
        zetas: published.iter().map(|&(_, q)| 1.0 / (2.0 * q)).collect(),
        shapes: vec![vec![1.0], vec![0.7], vec![0.5]],
    };
    let omega: Vec<f64> = (0i32..1200)
        .map(|i| TWO_PI * f64::from(i).mul_add(0.35, 60.0))
        .collect();
    let mut h = synth_frf(&truth, &omega, 0);
    add_noise(&mut h, 1.0e-3, 7); // measurement-grade noise
    let data = FrfData::new(omega, vec![(h, None)]).expect("data");
    let opts = IdentifyOptions {
        min_order: 4,
        max_order: 16,
        order_step: 2,
        ..IdentifyOptions::default()
    };
    let id = identify(&data, &opts).expect("identify");
    assert_eq!(id.modes.len(), 3, "expected the 3 published modes");
    for (mode, &(f_pub, q_pub)) in id.modes.iter().zip(&published) {
        let q_id = 1.0 / (2.0 * mode.damping_ratio);
        assert!(
            (mode.frequency_hz - f_pub).abs() <= 0.005 * f_pub,
            "frequency {:.2} vs published {f_pub}",
            mode.frequency_hz
        );
        // Q published as integers: 1-count granularity plus noise.
        assert!(
            (q_id - q_pub).abs() <= 0.05 * q_pub + 1.0,
            "Q {q_id:.1} vs published {q_pub}"
        );
    }
    println!(
        "{{\"suite\":\"fs-modalid\",\"case\":\"carcagno-benchmark\",\"modes\":{:?},\"verdict\":\"pass\"}}",
        id.modes
            .iter()
            .map(|m| (m.frequency_hz, 1.0 / (2.0 * m.damping_ratio)))
            .collect::<Vec<_>>()
    );
}

#[test]
fn window_correction_formula() {
    // zeta_measured = zeta + 1/(tau wn) exactly for pole shifting.
    let (zc, delta) = correct_exponential_window(0.02, 1000.0, 0.1);
    assert!((delta - 0.01).abs() < 1.0e-15);
    assert!((zc - 0.01).abs() < 1.0e-15);
    // Clamp at zero, delta still reported.
    let (zc2, delta2) = correct_exponential_window(0.005, 1000.0, 0.1);
    assert!(zc2.abs() < f64::MIN_POSITIVE, "clamped at zero");
    assert!((delta2 - 0.01).abs() < 1.0e-15);
}
