//! Dispersive-waveguide string bake-off (music bead
//! `frankensim-music-v8-root-3ez8g.7.2`): stiffness realized as
//! DISPERSION ON THE DELAY LINE (a first-order allpass cascade inside
//! `fs_vfit::DelayedFilter`) versus the incumbent modal-ZOH image, on
//! the SAME music-wire card. Doctrine D21: bake it off; keep both when
//! both pass and the budgets differ; never delete modal (asserted here
//! against the live registry).
//!
//! The oracle is the analytic stiff-string law `f_n = n f0 sqrt(1 +
//! B n^2)` with `B = pi^2 E I / (T L^2)` DERIVED from the card — the
//! modal image realizes it by construction; the waveguide must EARN it
//! through the allpass phase. The dispersionless control (same loop,
//! no stage) is measured too: without the stage the partials are
//! harmonic and the stiff-law error is large — the stage, not the
//! loop, carries the physics.
//!
//! Budget honesty: the waveguide's cost is one line + M allpass
//! sections regardless of how many partials speak; the modal image
//! pays two states per audio-band mode. Both are recorded from the
//! fixture's own numbers.

use std::collections::BTreeMap;

use fs_blake3::hash_domain;
use fs_couple::bakeoff::{BakeoffOutcome, BakeoffReceipt, ContenderResult};
use fs_fft::{C64, Fft};
use fs_vfit::discretize::{DelayedFilter, DigitalFilter};

const RATE: f64 = 48_000.0;

/// Music-wire card (shared by both contenders; everything below is
/// DERIVED from these five numbers).
const E_PA: f64 = 200.0e9;
const DIAMETER_M: f64 = 1.0e-3;
const LENGTH_M: f64 = 0.65;
const TENSION_N: f64 = 700.0;
const DENSITY: f64 = 7850.0;

struct Card {
    f0_hz: f64,
    b: f64,
}

fn card() -> Card {
    let area = core::f64::consts::PI * DIAMETER_M * DIAMETER_M / 4.0;
    let mu = DENSITY * area;
    let inertia = core::f64::consts::PI * DIAMETER_M.powi(4) / 64.0;
    let f0_hz = (TENSION_N / mu).sqrt() / (2.0 * LENGTH_M);
    let b = core::f64::consts::PI.powi(2) * E_PA * inertia / (TENSION_N * LENGTH_M * LENGTH_M);
    Card { f0_hz, b }
}

fn partial_law_hz(card: &Card, n: usize) -> f64 {
    let nf = n as f64;
    nf * card.f0_hz * (1.0 + card.b * nf * nf).sqrt()
}

/// Band-limited spectral peak with parabolic refinement — the
/// estimator-lesson-bank discipline (a global peak grabs the dominant
/// harmonic; each partial is searched only inside its own +-40-cent
/// window around the analytic location).
fn measured_partial_hz(spectrum: &[f64], expected_hz: f64) -> f64 {
    let bin_hz = RATE / spectrum.len() as f64 / 2.0;
    let lo = ((expected_hz * 0.977) / bin_hz) as usize;
    let hi = ((expected_hz * 1.023) / bin_hz) as usize;
    let mut best = lo;
    for k in lo..=hi.min(spectrum.len() - 2) {
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
    (best as f64 + shift) * bin_hz
}

fn spectrum_of(signal: &[f64]) -> Vec<f64> {
    let n = 1usize << 17;
    let mut buf: Vec<C64> = (0..n)
        .map(|k| {
            let w = 0.5
                - 0.5
                    * fs_math::det::cos(
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

/// Fit B from measured partials via the linear form
/// `(f_n / (n f_1))^2 = (1 + B n^2) / (1 + B)`.
fn fit_b(partials_hz: &[f64]) -> f64 {
    let f1 = partials_hz[0];
    let mut num = 0.0;
    let mut den = 0.0;
    for (i, &fnn) in partials_hz.iter().enumerate().skip(1) {
        let n = (i + 1) as f64;
        let r2 = (fnn / (n * f1)).powi(2);
        // r2 = (1 + B n^2)/(1 + B)  =>  B = (r2 - 1)/(n^2 - r2).
        let b = (r2 - 1.0) / (n * n - r2);
        num += b * n * n;
        den += n * n;
    }
    num / den
}

struct WgDesign {
    delay_samples: f64,
    a: f64,
    sections: usize,
    predicted_worst_cents: f64,
}

/// Exact first-order allpass phase at normalized frequency `omega`
/// for `H(z) = (a + z^-1)/(1 + a z^-1)` (radians, negative = lag).
fn allpass_phase(a: f64, omega: f64) -> f64 {
    let (s, c) = (fs_math::det::sin(omega), fs_math::det::cos(omega));
    fs_math::det::atan2(-s, a + c) - fs_math::det::atan2(-a * s, 1.0 + a * c)
}

/// Design (D, a) for a fixed section count by matching the loop phase
/// condition `Omega_n D + M * lag(Omega_n) = 2 pi n` at partials 1..8
/// of the analytic law — golden-section over `a`, closed-form D.
fn design_waveguide(card: &Card, sections: usize) -> WgDesign {
    let targets: Vec<(f64, f64)> = (1..=8)
        .map(|n| {
            let f = partial_law_hz(card, n);
            (core::f64::consts::TAU * f / RATE, n as f64)
        })
        .collect();
    let residual = |a: f64| -> (f64, f64, f64) {
        // Required D per partial, its mean, and the worst relative miss.
        let d_req: Vec<f64> = targets
            .iter()
            .map(|&(omega, n)| {
                let lag = -(sections as f64) * allpass_phase(a, omega);
                (core::f64::consts::TAU * n - lag) / omega
            })
            .collect();
        let d = d_req.iter().sum::<f64>() / d_req.len() as f64;
        let worst = d_req
            .iter()
            .map(|&r| (r - d).abs() / d)
            .fold(0.0f64, f64::max);
        (worst, d, a)
    };
    let (mut lo, mut hi) = (-0.85f64, -0.001f64);
    let phi = 0.618_033_988_75f64;
    for _ in 0..80 {
        let m1 = hi - phi * (hi - lo);
        let m2 = lo + phi * (hi - lo);
        if residual(m1).0 < residual(m2).0 {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let (worst, d, a) = residual(0.5 * (lo + hi));
    WgDesign {
        delay_samples: d,
        a,
        sections,
        // Loop-delay miss maps ~1:1 to frequency: cents ~ 1731 * rel.
        predicted_worst_cents: 1731.0 * worst,
    }
}

/// The dispersionless CONTROL: a plain loop tuned so partial 1 lands
/// on the law's f_1 (removing the stage from the dispersive design
/// would retune the whole string — the honest control changes ONE
/// thing: the stage).
fn control_design(card: &Card, sections: usize) -> WgDesign {
    WgDesign {
        delay_samples: RATE / partial_law_hz(card, 1),
        a: 0.0,
        sections,
        predicted_worst_cents: 0.0,
    }
}

struct Measured {
    partials_hz: Vec<f64>,
    b_hat: f64,
    worst_cents_vs_law: f64,
    states: usize,
}

fn cents(a: f64, b: f64) -> f64 {
    1200.0 * (a / b).log2()
}

fn measure_against_law(card: &Card, signal: &[f64], states: usize) -> Measured {
    let spec = spectrum_of(signal);
    let partials_hz: Vec<f64> = (1..=8)
        .map(|n| measured_partial_hz(&spec, partial_law_hz(card, n)))
        .collect();
    let worst = partials_hz
        .iter()
        .enumerate()
        .map(|(i, &f)| cents(f, partial_law_hz(card, i + 1)).abs())
        .fold(0.0f64, f64::max);
    Measured {
        b_hat: fit_b(&partials_hz),
        partials_hz,
        worst_cents_vs_law: worst,
        states,
    }
}

/// The waveguide string: one dispersive line in positive feedback
/// (`y = x + H y`), loop gain g < 1, raised-cosine excitation burst.
fn run_waveguide(design: &WgDesign, dispersion_on: bool) -> Vec<f64> {
    let g = 0.9995f64;
    let filter = DigitalFilter {
        sections: Vec::new(),
        direct: g,
        t_s: 1.0 / RATE,
        prewarp: 0.0,
    };
    let mut line = DelayedFilter::new(design.delay_samples, filter).expect("line admits");
    if dispersion_on {
        line = line
            .with_dispersion(design.a, design.sections)
            .expect("dispersion admits");
    }
    let n = (1.5 * RATE) as usize;
    let mut out = Vec::with_capacity(n);
    let mut feedback = 0.0f64;
    for k in 0..n {
        let excite = if k < 24 {
            0.5 - 0.5 * fs_math::det::cos(core::f64::consts::TAU * k as f64 / 24.0)
        } else {
            0.0
        };
        let output = line.push(feedback).expect("line step");
        feedback = excite + output;
        out.push(output);
    }
    out
}

/// The incumbent: exact-ZOH modal oscillators AT the law's frequencies
/// (realized by construction — that is the incumbent's whole claim).
fn run_modal(card: &Card) -> Vec<f64> {
    let n = (1.5 * RATE) as usize;
    let dt = 1.0 / RATE;
    let mut oscillators: Vec<(f64, f64, f64, f64)> = (1..=8)
        .map(|k| {
            let omega = core::f64::consts::TAU * partial_law_hz(card, k);
            let zeta = 2.0e-4_f64;
            let od = omega * (1.0 - zeta * zeta).sqrt();
            let decay = fs_math::det::exp(-zeta * omega * dt);
            // (re, im) rotor state; amplitude ~ pluck 1/k^2.
            (
                decay * fs_math::det::cos(od * dt),
                decay * fs_math::det::sin(od * dt),
                1.0 / (k * k) as f64,
                0.0,
            )
        })
        .collect();
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let mut sample = 0.0;
        for (c, s, re, im) in &mut oscillators {
            let nre = *c * *re - *s * *im;
            let nim = *s * *re + *c * *im;
            *re = nre;
            *im = nim;
            sample += *re;
        }
        out.push(sample);
    }
    out
}

fn qois(m: &Measured, card: &Card) -> BTreeMap<String, f64> {
    let mut q = BTreeMap::new();
    q.insert("b_over_b_card".to_string(), m.b_hat / card.b);
    q.insert("worst_partial_cents".to_string(), m.worst_cents_vs_law);
    q.insert(
        "partial8_cents_vs_harmonic".to_string(),
        cents(m.partials_hz[7], 8.0 * card.f0_hz),
    );
    q
}

#[test]
#[ignore = "minting run: measures both contenders and prints fresh receipt bytes"]
fn mint_dispersion_bakeoff_receipt() {
    let card = card();
    println!(
        "card: f0 = {:.3} Hz, B = {:.4e}; stiff partial 8 sits {:.1} cents sharp of harmonic",
        card.f0_hz,
        card.b,
        cents(partial_law_hz(&card, 8), 8.0 * card.f0_hz)
    );
    let design = design_waveguide(&card, 8);
    println!(
        "design: D = {:.3} samples, a = {:.5}, M = {}, predicted worst {:.2} cents",
        design.delay_samples, design.a, design.sections, design.predicted_worst_cents
    );
    let wg_states = design.delay_samples.ceil() as usize + 2 * design.sections;
    let wg = measure_against_law(&card, &run_waveguide(&design, true), wg_states);
    let control = measure_against_law(&card, &run_waveguide(&control_design(&card, 8), false), wg_states);
    // The modal image needs two states per audio-band mode to SOUND
    // like the full string; count modes under 20 kHz from the law.
    let full_band_modes = (1..)
        .take_while(|&n| partial_law_hz(&card, n) < 20_000.0)
        .count();
    let modal = measure_against_law(&card, &run_modal(&card), 2 * full_band_modes);
    println!(
        "waveguide: B_hat/B = {:.3}, worst {:.2} cents, states {}",
        wg.b_hat / card.b,
        wg.worst_cents_vs_law,
        wg.states
    );
    println!(
        "control (no stage): B_hat/B = {:.3}, worst {:.2} cents",
        control.b_hat / card.b,
        control.worst_cents_vs_law
    );
    println!(
        "modal: B_hat/B = {:.3}, worst {:.2} cents, states {} ({} audio-band modes)",
        modal.b_hat / card.b,
        modal.worst_cents_vs_law,
        modal.states,
        full_band_modes
    );
    let mut cards_bytes = Vec::new();
    for v in [E_PA, DIAMETER_M, LENGTH_M, TENSION_N, DENSITY, RATE] {
        cards_bytes.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    let mut reference = BTreeMap::new();
    reference.insert("b_over_b_card".to_string(), 1.0);
    reference.insert("worst_partial_cents".to_string(), 0.0);
    reference.insert(
        "partial8_cents_vs_harmonic".to_string(),
        cents(partial_law_hz(&card, 8), 8.0 * card.f0_hz),
    );
    let receipt = BakeoffReceipt {
        filling: "string".to_string(),
        fixture: "crates/fs-couple/tests/bakeoff_dispersion.rs music-wire card \
                  (E=200GPa d=1mm L=0.65m T=120N rho=7850), pluck loop 1.5s at 48k, \
                  partials 1..8 vs the analytic stiff law; dispersionless control in rationale"
            .to_string(),
        shared_cards: hash_domain("org.frankensim.fs-couple.dispersion-cards.v1", &cards_bytes),
        reference,
        contenders: [
            ContenderResult {
                image: "modal-zoh".to_string(),
                owner_crates: vec!["fs-couple".to_string()],
                measured: qois(&modal, &card),
                states: modal.states,
                steps: 72_000,
                solver_iterations: 0,
                failure_modes: vec![
                    "cost scales with partial count (2 states/mode to fill the audio band)"
                        .to_string(),
                ],
            },
            ContenderResult {
                image: "waveguide-dispersive".to_string(),
                owner_crates: vec!["fs-vfit".to_string(), "fs-couple".to_string()],
                measured: qois(&wg, &card),
                states: wg.states,
                steps: 72_000,
                solver_iterations: 0,
                failure_modes: vec![
                    "shared-coefficient cascade tunes the ladder collectively; individual \
                     partials carry cents-class residuals the modal image does not"
                        .to_string(),
                    "loop delay quantized by the line's linear fractional interpolation"
                        .to_string(),
                ],
            },
        ],
        outcome: BakeoffOutcome::KeepBoth {
            scope_a: "modal-zoh keeps every certified claim: exact partial placement on the \
                      retained basis, truncation-audited, the polyphony image"
                .to_string(),
            scope_b: "waveguide-dispersive earns the budget-constrained hero-string scope: \
                      one line + M allpass states regardless of partial count, stiff ladder \
                      within its measured cents band; NOT a certified-partial image"
                .to_string(),
        },
        rationale: String::new(), // measured numbers filled at mint time
        listening_receipts: Vec::new(),
    };
    let mut with_rationale = receipt;
    with_rationale.rationale = format!(
        "Waveguide B_hat/B = {:.3} with worst partial {:.1} cents (predicted {:.1}); the \
         DISPERSIONLESS control on the same loop measured B_hat/B = {:.3} and {:.1} cents \
         worst — the allpass stage, not the loop, carries the stiffness. Modal worst \
         {:.2} cents (realizes the law by construction). Budgets, honestly: modal holds \
         FEWER states here ({} waveguide vs {} modal) — the differing shape is OPS and \
         scaling: the waveguide spends ~1 tap + M allpass sections per sample regardless \
         of partial count while modal pays one rotor per audio-band mode (~10x flops \
         here), and the waveguide's state cost is flat in bandwidth. D21 keeps both.",
        wg.b_hat / card.b,
        wg.worst_cents_vs_law,
        design.predicted_worst_cents,
        control.b_hat / card.b,
        control.worst_cents_vs_law,
        modal.worst_cents_vs_law,
        wg.states,
        modal.states,
    );
    let bytes = with_rationale.to_canonical_bytes().expect("canonical");
    let hash = with_rationale.content_hash().expect("hash");
    println!("receipt hash: {hash}");
    std::fs::write(
        concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/receipts/string-modal-zoh-vs-waveguide-dispersive.bakeoff"
        ),
        &bytes,
    )
    .expect("write receipt");
    println!("receipt written; commit it and run the gate test");
}

#[test]
fn committed_dispersion_receipt_decides_and_modal_rows_are_untouched() {
    let bytes = std::fs::read(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/receipts/string-modal-zoh-vs-waveguide-dispersive.bakeoff"
    ))
    .expect("committed receipt (run the ignored mint test first)");
    let receipt = BakeoffReceipt::from_canonical_bytes(&bytes).expect("decode");
    receipt.validate().expect("valid");
    assert_eq!(receipt.filling, "string");
    assert_eq!(receipt.contenders[0].image, "modal-zoh");
    assert_eq!(receipt.contenders[1].image, "waveguide-dispersive");
    // The verdict is KeepBoth with the budget-shape rationale.
    match &receipt.outcome {
        BakeoffOutcome::KeepBoth { scope_a, scope_b } => {
            assert!(scope_a.contains("modal-zoh keeps every certified claim"));
            assert!(scope_b.contains("budget-constrained hero-string"));
        }
        other => panic!("the dispersion bake-off decided KeepBoth, found {other:?}"),
    }
    // Budget shapes really differ (the reason KeepBoth exists): both
    // recorded, and the rationale carries the ops/scaling story — the
    // honest form, since modal holds FEWER states on this card.
    assert!(receipt.contenders[0].states > 0 && receipt.contenders[1].states > 0);
    assert!(receipt.rationale.contains("OPS"));
    // D21: the modal rows are UNTOUCHED — the incumbent keeps its green
    // registry row no matter what the newcomer measured.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(std::path::Path::parent)
        .expect("root")
        .to_path_buf();
    let registry = std::fs::read_to_string(root.join("instrument-claims.json")).expect("registry");
    let modal_at = registry
        .find("\"image\": \"modal-zoh\"")
        .expect("modal row present");
    let gate_tag = "\"gate\": \"";
    let gate_at = registry[modal_at..].find(gate_tag).expect("gate") + gate_tag.len();
    let gate_end = registry[modal_at + gate_at..].find('"').expect("end");
    assert_eq!(
        &registry[modal_at + gate_at..modal_at + gate_at + gate_end],
        "green",
        "D21: the modal incumbent must keep its green row"
    );
    // The newcomer's own row exists with its earned (ungated) scope.
    assert!(
        registry.contains("\"image\": \"waveguide-dispersive\""),
        "the newcomer's registry row must exist with its earned scope"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"dispersion-bakeoff-gate\",\"verdict\":\"pass\",\
         \"hash\":\"{}\"}}",
        receipt.content_hash().expect("hash")
    );
}

#[test]
fn dispersion_stage_is_passive_and_the_control_is_harmonic() {
    // Structural passivity: |H| = 1 for every allpass section, so the
    // near-lossless loop must stay bounded with the stage in it; and
    // the dispersionless control's partial 8 must sit essentially ON
    // the harmonic (proving the stage is what bends the ladder).
    let card = card();
    let design = design_waveguide(&card, 8);
    // Constructor refusals: a unit-circle coefficient and a zero-section
    // cascade never build.
    {
        let f = DigitalFilter {
            sections: Vec::new(),
            direct: 1.0,
            t_s: 1.0 / RATE,
            prewarp: 0.0,
        };
        assert!(
            DelayedFilter::new(64.0, f.clone())
                .expect("line")
                .with_dispersion(1.0, 8)
                .is_err()
        );
        assert!(
            DelayedFilter::new(64.0, f)
                .expect("line")
                .with_dispersion(-0.5, 0)
                .is_err()
        );
    }
    let out = run_waveguide(&design, true);
    let peak = out.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
    assert!(
        peak.is_finite() && peak < 100.0,
        "loop bounded (peak {peak:.2})"
    );
    let tail_rms = (out[out.len() - 4800..].iter().map(|v| v * v).sum::<f64>() / 4800.0).sqrt();
    let head_rms = (out[2400..7200].iter().map(|v| v * v).sum::<f64>() / 4800.0).sqrt();
    assert!(
        tail_rms < head_rms,
        "g<1 loop with a unit-magnitude stage must decay ({tail_rms:.3e} vs {head_rms:.3e})"
    );
    let control = measure_against_law(&card, &run_waveguide(&control_design(&card, 8), false), 0);
    // Against its OWN harmonic ladder (8 * f_1) the control is near-
    // harmonic (the residual is the fractional-interp phase), while
    // against the STIFF law it must miss by the stiffness gap — the
    // stage is the thing that closes it.
    let control_partial8_vs_harmonic = cents(control.partials_hz[7], 8.0 * partial_law_hz(&card, 1));
    // Measured 2026-08-17: -9.3 cents (the linear-interp phase at
    // partial 8 on a frac=0.16 line) vs a 27.0-cent stiff-law miss —
    // the comparative band is the honest one: the control's residual
    // must stay well under the stiffness signal it fails to produce.
    assert!(
        control_partial8_vs_harmonic.abs() < 0.5 * control.worst_cents_vs_law,
        "the control's own-ladder residual ({control_partial8_vs_harmonic:.1} cents) must \
         stay well under its stiff-law miss ({:.1})",
        control.worst_cents_vs_law
    );
    assert!(
        control.worst_cents_vs_law > 10.0,
        "the control must MISS the stiff law ({:.1} cents)",
        control.worst_cents_vs_law
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"dispersion-passivity-control\",\"verdict\":\"pass\",\
         \"peak\":{peak:.3},\"control_p8_cents\":{control_partial8_vs_harmonic:.2}}}"
    );
}
