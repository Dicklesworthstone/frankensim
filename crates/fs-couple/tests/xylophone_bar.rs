//! Xylophone bar: free-free Euler-Bernoulli modal binding + strike
//! gates (music bead `frankensim-music-v8-root-3ez8g.12.1`) — the
//! cheapest FULL filling in the program: analytic bar modes + exact-ZOH
//! runtime + a Hertz strike island, no new physics, only binding and
//! gates.
//!
//! THE BOUNDARY CONDITION IS LOAD-BEARING (the bead's recorded polish
//! catch): a xylophone bar is FREE-FREE — its `beta L` values are the
//! roots of `cosh(x) cos(x) = 1` (4.73004, 7.85320, 10.99561, ...)
//! computed IN-FIXTURE by deterministic Newton (self-verified at
//! machine zero, never transcription-trust), giving the NON-HARMONIC
//! ladder `f2/f1 = (7.85320/4.73004)^2 ~ 2.7565`. The falsifier proves
//! the gate discriminates: the PINNED family (`f_n ~ n^2`, the
//! `fs_nlmodal::prestressed_beam_omega` boundary conditions) measures
//! `f2/f1 = 4` and FAILS the same gate.
//!
//! Material: the Indian-rosewood matdb pack (FPL-GTR-282: E_L = 13530
//! MPa; specific gravity 0.75 — density is the authored Estimate
//! SG*1000 since the pack derives no 12-percent density, disclosed).
//! Bar geometry is an authored chart. The UNDERCUT successor is named:
//! real bars are arched to tune f2/f1 toward 3-4; v1 claims the
//! uniform PRISMATIC bar exactly and only that.

use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_math::det;

const RATE: u32 = 48_000;

/// Authored bar chart + the rosewood pack numbers.
const E_PA: f64 = 13.53e9;
const DENSITY: f64 = 750.0;
const BAR_L: f64 = 0.38;
const BAR_W: f64 = 0.038;
const BAR_T: f64 = 0.018;

/// Roots of `cos(x) = 1/cosh(x)` (the well-conditioned form of
/// `cosh cos = 1`) by deterministic Newton from the asymptotic guesses
/// `(2k+1) pi / 2`.
fn free_free_roots(count: usize) -> Vec<f64> {
    (1..=count)
        .map(|k| {
            let mut x = (2.0 * k as f64 + 1.0) * core::f64::consts::PI / 2.0;
            for _ in 0..60 {
                let g = det::cos(x) - 1.0 / x.cosh();
                let dg = -det::sin(x) + x.sinh() / (x.cosh() * x.cosh());
                let step = g / dg;
                x -= step;
                if step.abs() < 1e-14 {
                    break;
                }
            }
            x
        })
        .collect()
}

/// Free-free mode shape at `xi` in [0, 1] for root `x`.
fn mode_shape(x: f64, xi: f64) -> f64 {
    let sigma = (x.cosh() - det::cos(x)) / (x.sinh() - det::sin(x));
    let a = x * xi;
    a.cosh() + det::cos(a) - sigma * (a.sinh() + det::sin(a))
}

struct BarModes {
    frequencies_hz: Vec<f64>,
    roots: Vec<f64>,
}

/// Bind the bar card to modal frequencies, with refusals.
fn bar_modes(
    e_pa: f64,
    density: f64,
    length: f64,
    width: f64,
    thickness: f64,
    count: usize,
) -> Result<BarModes, &'static str> {
    for (v, what) in [
        (e_pa, "young modulus"),
        (density, "density"),
        (length, "length"),
        (width, "width"),
        (thickness, "thickness"),
    ] {
        if !(v.is_finite() && v > 0.0) {
            return Err(what);
        }
    }
    let inertia = width * thickness * thickness * thickness / 12.0;
    let area = width * thickness;
    let radius = (e_pa * inertia / (density * area)).sqrt();
    let roots = free_free_roots(count);
    let frequencies_hz = roots
        .iter()
        .map(|&x| x * x / (core::f64::consts::TAU * length * length) * radius)
        .collect();
    Ok(BarModes {
        frequencies_hz,
        roots,
    })
}

/// Strike the bar at `xi` with a Hertz mallet; return the audio.
fn strike(modes: &BarModes, xi: f64, seconds: f64) -> Vec<f64> {
    let n_modes = modes.frequencies_hz.len();
    let model_modes: Vec<ModalAcousticMode> = modes
        .frequencies_hz
        .iter()
        .map(|&f| ModalAcousticMode {
            angular_frequency_rad_s: core::f64::consts::TAU * f,
            damping_ratio: 8.0e-4,
            pressure_per_modal_velocity: fs_math::c64::C64::new(1.0, 0.0),
        })
        .collect();
    let mut bar = ModalAcousticTimeModel::try_new(
        RATE,
        model_modes,
        ModalAcousticTimeBudget::audible_reference(),
    )
    .expect("bar modes admit");
    let phi: Vec<f64> = modes.roots.iter().map(|&x| mode_shape(x, xi)).collect();
    let dt = 1.0 / f64::from(RATE);
    let substeps = 64usize;
    let h = dt / substeps as f64;
    let mallet_m = 0.030f64;
    let hertz_k = 2.0e8f64; // rubber-headed mallet, authored
    let mut y_mallet = 0.0f64;
    let mut v_mallet = 1.5f64;
    let n = (seconds * f64::from(RATE)) as usize;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let bar_disp: f64 = bar
            .states()
            .iter()
            .zip(&phi)
            .map(|(s, p)| s.displacement_m_sqrt_kg * p)
            .sum();
        let mut force_avg = 0.0f64;
        if v_mallet != 0.0 || y_mallet > bar_disp {
            for _ in 0..substeps {
                let overlap = y_mallet - bar_disp;
                let f = if overlap > 0.0 {
                    hertz_k * overlap * overlap.sqrt()
                } else {
                    0.0
                };
                force_avg += f;
                v_mallet -= h * f / mallet_m;
                y_mallet += h * v_mallet;
            }
            force_avg /= substeps as f64;
        }
        let generalized: Vec<f64> = phi.iter().map(|p| force_avg * p).collect();
        let frame = bar.step(&generalized).expect("bar step");
        let _ = frame;
        let sample: f64 = bar
            .states()
            .iter()
            .zip(&phi)
            .map(|(s, p)| s.velocity_m_sqrt_kg_per_s * p)
            .sum();
        out.push(sample);
        // Once the mallet has clearly left, freeze it (no re-strike).
        if v_mallet < 0.0 && y_mallet - bar_disp < -1.0e-3 {
            v_mallet = 0.0;
            y_mallet = -1.0;
        }
        let _ = n_modes;
    }
    out
}

fn spectrum_of(signal: &[f64]) -> Vec<f64> {
    use fs_fft::{C64, Fft};
    let n = 1usize << 16;
    let mut buf: Vec<C64> = (0..n)
        .map(|k| {
            let w = 0.5
                - 0.5
                    * det::cos(
                        core::f64::consts::TAU * k as f64 / (signal.len().min(n) - 1) as f64,
                    );
            C64::new(signal.get(k).copied().unwrap_or(0.0) * w, 0.0)
        })
        .collect();
    let mut scratch = vec![C64::new(0.0, 0.0); n];
    Fft::new(n).forward(&mut buf, &mut scratch);
    buf[..n / 2]
        .iter()
        .map(|c| (c.re * c.re + c.im * c.im).sqrt())
        .collect()
}

/// Peak magnitude near `expected_hz` (+-3%), parabolic-refined
/// frequency and the peak magnitude.
fn partial_peak(spectrum: &[f64], expected_hz: f64) -> (f64, f64) {
    let bin_hz = f64::from(RATE) / (2.0 * spectrum.len() as f64);
    let lo = ((expected_hz * 0.97) / bin_hz) as usize;
    let hi = (((expected_hz * 1.03) / bin_hz) as usize).min(spectrum.len() - 2);
    let mut best = lo;
    for k in lo..=hi {
        if spectrum[k] > spectrum[best] {
            best = k;
        }
    }
    let (a, b, c) = (
        spectrum[best - 1].max(1e-30).ln(),
        spectrum[best].max(1e-30).ln(),
        spectrum[best + 1].max(1e-30).ln(),
    );
    let denom = a - 2.0 * b + c;
    let shift = if denom.abs() > 1e-12 {
        0.5 * (a - c) / denom
    } else {
        0.0
    };
    ((best as f64 + shift) * bin_hz, spectrum[best])
}

#[test]
fn xb_001_roots_self_verify_and_match_the_tabulated_values() {
    // Self-verified pin: the Newton roots satisfy the characteristic
    // equation at machine zero AND match the classical tabulated
    // beta*L values at 1e-10 (a cross-check, not transcription-trust),
    // AND approach the (2n+1) pi/2 asymptote from the right side.
    let roots = free_free_roots(5);
    let table = [
        4.730_040_74,
        7.853_204_62,
        10.995_607_84,
        14.137_165_49,
        17.278_759_66,
    ];
    for (i, (&x, &t)) in roots.iter().zip(table.iter()).enumerate() {
        let residual = det::cos(x) - 1.0 / x.cosh();
        assert!(residual.abs() < 1e-13, "root {i} residual {residual:.2e}");
        assert!((x - t).abs() < 1e-8, "root {i}: {x} vs table {t}");
        let asymptote = (2.0 * (i + 1) as f64 + 1.0) * core::f64::consts::PI / 2.0;
        assert!(
            (x - asymptote).abs() < 0.02 * asymptote,
            "root {i} must sit near its asymptote"
        );
    }
    // Refusals by name.
    for (result, expected) in [
        (
            bar_modes(-1.0, DENSITY, BAR_L, BAR_W, BAR_T, 4),
            "young modulus",
        ),
        (
            bar_modes(E_PA, DENSITY, f64::NAN, BAR_W, BAR_T, 4),
            "length",
        ),
        (bar_modes(E_PA, DENSITY, BAR_L, BAR_W, 0.0, 4), "thickness"),
    ] {
        match result {
            Err(what) => assert_eq!(what, expected),
            Ok(_) => panic!("{expected} refusal must fire"),
        }
    }
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"xb-001-roots\",\"verdict\":\"pass\",\
         \"beta_l\":[{:.6},{:.6},{:.6}]}}",
        roots[0], roots[1], roots[2]
    );
}

#[test]
fn xb_002_ratio_gate_holds_and_the_pinned_falsifier_fails_it() {
    // THE GATE THAT MATTERS: the rendered strike's partial ratios must
    // land on the free-free ladder (f2/f1 = 2.7565, f3/f1 = 5.4039),
    // and the PINNED family on the same f1 (f_n = n^2 f1, the
    // prestressed_beam_omega boundary conditions) must FAIL the same
    // gate decisively — the gate provably discriminates the boundary
    // conditions.
    let modes = bar_modes(E_PA, DENSITY, BAR_L, BAR_W, BAR_T, 4).expect("bar");
    let audio = strike(&modes, 0.30, 1.0);
    let spec = spectrum_of(&audio);
    let (f1, _) = partial_peak(&spec, modes.frequencies_hz[0]);
    let (f2, _) = partial_peak(&spec, modes.frequencies_hz[1]);
    let (f3, _) = partial_peak(&spec, modes.frequencies_hz[2]);
    let r21 = f2 / f1;
    let r31 = f3 / f1;
    let law21 = (7.853_204_62f64 / 4.730_040_74).powi(2);
    let law31 = (10.995_607_84f64 / 4.730_040_74).powi(2);
    assert!(
        (r21 - law21).abs() < 0.01 * law21,
        "f2/f1 measured {r21:.4} vs free-free {law21:.4}"
    );
    assert!(
        (r31 - law31).abs() < 0.01 * law31,
        "f3/f1 measured {r31:.4} vs free-free {law31:.4}"
    );
    // FALSIFIER: pinned-family frequencies through the SAME pipeline.
    let pinned = BarModes {
        frequencies_hz: (1..=4)
            .map(|n| f64::from(n * n) * modes.frequencies_hz[0])
            .collect(),
        roots: modes.roots.clone(),
    };
    let pinned_audio = strike(&pinned, 0.30, 1.0);
    let pinned_spec = spectrum_of(&pinned_audio);
    let (p1, _) = partial_peak(&pinned_spec, pinned.frequencies_hz[0]);
    let (p2, _) = partial_peak(&pinned_spec, pinned.frequencies_hz[1]);
    let pinned_r21 = p2 / p1;
    assert!(
        (pinned_r21 - law21).abs() > 0.3 * law21,
        "the pinned family (f2/f1 = {pinned_r21:.3}) must fail the free-free gate"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"xb-002-ratio-gate\",\"verdict\":\"pass\",\
         \"f1\":{f1:.2},\"r21\":{r21:.4},\"r31\":{r31:.4},\"pinned_r21\":{pinned_r21:.3}}}"
    );
}

#[test]
fn xb_003_striking_a_node_suppresses_the_partial() {
    // Strike-position physics, emergent: mode 2 of a free-free bar has
    // a node at the center, and mode 1 at xi ~ 0.224 (the cord-mount
    // point on a real bar) — both computed in-fixture from the mode
    // shapes, both suppressing their partial when struck.
    let modes = bar_modes(E_PA, DENSITY, BAR_L, BAR_W, BAR_T, 4).expect("bar");
    // Mode-1 node by bisection on the in-fixture shape.
    let x1 = modes.roots[0];
    let (mut lo, mut hi) = (0.05f64, 0.45f64);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if mode_shape(x1, mid).signum() == mode_shape(x1, lo).signum() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let node1 = 0.5 * (lo + hi);
    assert!(
        (node1 - 0.224).abs() < 0.01,
        "mode-1 node at {node1:.4} (classical 0.2242)"
    );
    let energy_at = |xi: f64| -> Vec<f64> {
        let audio = strike(&modes, xi, 0.8);
        let spec = spectrum_of(&audio);
        modes
            .frequencies_hz
            .iter()
            .map(|&f| partial_peak(&spec, f).1)
            .collect()
    };
    let generic = energy_at(0.30);
    let center = energy_at(0.5);
    let at_node1 = energy_at(node1);
    for (i, (&g, (&c, &n1))) in generic
        .iter()
        .zip(center.iter().zip(at_node1.iter()))
        .enumerate()
    {
        println!(
            "partial {}: generic {:.3e}, center {:.3e}, node1 {:.3e}",
            i + 1,
            g,
            c,
            n1
        );
    }
    // Mode 2 (antisymmetric) dies at the center strike; mode 1 dies at
    // its own node. Suppression = >20 dB below the generic strike.
    assert!(
        center[1] < 0.1 * generic[1],
        "center strike must suppress partial 2 ({:.3e} vs {:.3e})",
        center[1],
        generic[1]
    );
    assert!(
        at_node1[0] < 0.1 * generic[0],
        "node-1 strike must suppress partial 1 ({:.3e} vs {:.3e})",
        at_node1[0],
        generic[0]
    );
    // And the partials NOT at a node keep speaking.
    assert!(center[0] > 0.3 * generic[0]);
    assert!(at_node1[1] > 0.3 * generic[1]);
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"xb-003-node-suppression\",\"verdict\":\"pass\",\
         \"node1\":{node1:.4}}}"
    );
}

fn repo_root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf()
}

fn sidecar_field(sidecar: &str, key: &str) -> String {
    let tag = format!("\"{key}\":\"");
    let at = sidecar.find(&tag).expect("field") + tag.len();
    let end = sidecar[at..].find('"').expect("end");
    sidecar[at..at + end].to_string()
}

#[test]
#[ignore = "minting run: renders the bar strikes and writes the artifact + receipt"]
fn mint_xylophone_strike_artifact() {
    use fs_psycho::receipt::{ListeningReceipt, ListeningVerdict};
    // Three strikes: generic station (full ladder), center (partial 2
    // suppressed), mode-1 node (fundamental suppressed) — the
    // strike-position physics audible, not asserted.
    let modes = bar_modes(E_PA, DENSITY, BAR_L, BAR_W, BAR_T, 4).expect("bar");
    let mut signal = Vec::new();
    for &xi in &[0.30f64, 0.5, 0.2242] {
        signal.extend(strike(&modes, xi, 0.9));
        signal.extend(std::iter::repeat_n(0.0, 4_800));
    }
    let mean = signal.iter().sum::<f64>() / signal.len() as f64;
    for v in &mut signal {
        *v -= mean;
    }
    let peak = signal.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
    let full_scale = peak * 1.4;
    let (wav, clipped) =
        fs_couple::pcm_wav::encode_pcm16_wav(&signal, RATE, full_scale).expect("wav");
    assert_eq!(clipped, 0, "never clip a listening artifact");
    let hash = fs_blake3::hash_domain("org.frankensim.music-render.wav.v1", &wav);
    let rms = (signal.iter().map(|v| v * v).sum::<f64>() / signal.len() as f64).sqrt();
    let root = repo_root();
    std::fs::write(root.join("data/listening/xylophone-strike.wav"), &wav).expect("wav write");
    let provenance = format!(
        "{{\"schema\":\"frankensim-music-render-provenance-v1\",\"fixture\":\"xylophone-strike \
         (free-free rosewood bar, strikes at xi = 0.30 / 0.50 / 0.2242: full ladder, \
         partial-2-suppressed, fundamental-suppressed; modal-velocity observer at the \
         strike station)\",\
         \"sample_rate_hz\":48000,\"samples\":{},\"block\":43200,\
         \"full_scale_pa\":{full_scale:e},\"clipped_samples\":0,\"peak_pa\":{peak:e},\
         \"rms_pa\":{rms:e},\"wav_blake3\":\"{}\",\"encoder\":\"fs_couple::pcm_wav (mono \
         PCM16, never peak-normalized)\"}}\n",
        signal.len(),
        hash.to_hex()
    );
    std::fs::write(
        root.join("data/listening/xylophone-strike.provenance.json"),
        provenance,
    )
    .expect("sidecar write");
    let lat = fs_psycho::log_attack_time(&signal, f64::from(RATE), 480).expect("attack time");
    let receipt = ListeningReceipt {
        listener: "pending".to_string(),
        session: "2026-08-17".to_string(),
        artifact_hex: hash.to_hex(),
        artifact_ref: "data/listening/xylophone-strike.provenance.json".to_string(),
        question: "does it read as a BAR (inharmonic 2.76/5.40 ladder), not a string — and \
                   does the third strike lose its fundamental?"
            .to_string(),
        verdict: ListeningVerdict::Unadjudicated,
        observations: "three strikes at 0.30/0.50/0.2242 of the free-free rosewood bar; \
                       awaiting the owner's ear"
            .to_string(),
        metrics: fs_psycho::receipt::AttachedMetrics {
            loudness_sone: None,
            sharpness_acum: None,
            log_attack_time: Some(lat),
            spl_db: None,
        },
    };
    let bytes = receipt.to_canonical_bytes().expect("encode");
    std::fs::write(
        root.join("data/listening/xylophone-strike.listening-receipt"),
        &bytes,
    )
    .expect("receipt write");
    println!(
        "minted xylophone-strike: {} samples, peak {peak:.3e}, hash {}",
        signal.len(),
        hash.to_hex()
    );
}

#[test]
fn xb_004_committed_receipt_matches_the_committed_artifact() {
    use fs_psycho::receipt::{ListeningReceipt, ListeningVerdict};
    let root = repo_root();
    let receipt_bytes =
        std::fs::read(root.join("data/listening/xylophone-strike.listening-receipt"))
            .expect("committed xylophone receipt (mint test)");
    let receipt = ListeningReceipt::from_canonical_bytes(&receipt_bytes).expect("decodes");
    let sidecar =
        std::fs::read_to_string(root.join("data/listening/xylophone-strike.provenance.json"))
            .expect("committed sidecar");
    assert_eq!(receipt.artifact_hex, sidecar_field(&sidecar, "wav_blake3"));
    let wav = std::fs::read(root.join("data/listening/xylophone-strike.wav")).expect("wav");
    let hash = fs_blake3::hash_domain("org.frankensim.music-render.wav.v1", &wav);
    assert_eq!(hash.to_hex(), receipt.artifact_hex, "WAV bytes drifted");
    assert_eq!(receipt.verdict, ListeningVerdict::Unadjudicated);
    assert!(!receipt.supports_pass());
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"xb-004-listening-chain\",\"artifact\":\"{}\"}}",
        receipt.artifact_hex
    );
}
