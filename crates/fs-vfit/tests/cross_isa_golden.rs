//! Cross-ISA determinism goldens for the `fs-vfit` fit + discretize +
//! filter-step pipeline (bead `frankensim-music-v8-root-3ez8g.13.4`).
//!
//! Vector fitting (pole relocation, residue pass), Tustin bilinear
//! discretization with prewarp, and DF-II filter stepping must produce
//! bit-identical results on both reference ISA families in both build
//! modes. The authored truth model is sampled through the crate's own
//! `eval_iw` (det-arithmetic), so any platform-libm leak on the pipeline
//! shows up as a digest mismatch — a golden event: bisect stage-wise,
//! name the hazard, route it through `det::`, same commit.

use fs_math::c64::C64;
use fs_math::det;
use fs_vfit::discretize::bilinear;
use fs_vfit::{FitOptions, PoleTerm, RationalModel, WeightPreset, vector_fit};

const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

fn fold(acc: u64, v: f64) -> u64 {
    v.to_bits()
        .to_le_bytes()
        .iter()
        .fold(acc, |a, &b| (a ^ u64::from(b)).wrapping_mul(FNV_PRIME))
}

fn fold_u64(acc: u64, v: u64) -> u64 {
    (acc ^ v).wrapping_mul(FNV_PRIME)
}

/// Authored passive-ish three-pole-pair truth model.
fn truth_model() -> RationalModel {
    let pair = |re: f64, im: f64, cre: f64, cim: f64| PoleTerm::Pair {
        pole: C64::new(re, im),
        residue: C64::new(cre, cim),
    };
    RationalModel {
        terms: vec![
            pair(-30.0, 2.0 * core::f64::consts::PI * 180.0, 80.0, 5.0),
            pair(-45.0, 2.0 * core::f64::consts::PI * 900.0, -55.0, 12.0),
            pair(-60.0, 2.0 * core::f64::consts::PI * 3_400.0, 34.0, -7.0),
        ],
        d: 0.02,
        e: 0.0,
    }
}

/// Log-spaced grid via det-arithmetic only.
fn log_grid(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    let lr = det::ln(hi / lo);
    (0..n)
        .map(|k| lo * det::exp(lr * k as f64 / (n - 1) as f64))
        .collect()
}

/// Verified bit-identical aarch64-apple and x86_64-linux (debug) on
/// 2026-08-23, bead frankensim-music-v8-root-3ez8g.13.4.
const GOLDEN_HASH: u64 = 0xd00d_69e2_b740_e56b;

#[test]
fn fit_and_filter_step_digest_is_cross_isa_golden() {
    let omega = log_grid(
        2.0 * core::f64::consts::PI * 20.0,
        2.0 * core::f64::consts::PI * 20.0e3,
        240,
    );
    let truth = truth_model();
    let h: Vec<C64> = omega.iter().map(|&w| truth.eval_iw(w)).collect();

    // Stage 1: vector fit of the sampled truth.
    let mut opts = FitOptions::new(6);
    opts.weights = WeightPreset::InverseMagnitude;
    opts.fit_e = false;
    opts.fit_d = true;
    let outcome = vector_fit(&omega, &h, &opts).expect("fit");
    let mut acc = FNV_OFFSET;
    acc = fold(acc, outcome.model.d);
    acc = fold(acc, outcome.model.e);
    for term in &outcome.model.terms {
        match *term {
            PoleTerm::Real { pole, residue } => {
                acc = fold(acc, pole);
                acc = fold(acc, residue);
            }
            PoleTerm::Pair { pole, residue } => {
                acc = fold(acc, pole.re);
                acc = fold(acc, pole.im);
                acc = fold(acc, residue.re);
                acc = fold(acc, residue.im);
            }
        }
    }
    acc = fold(acc, outcome.report.weighted_rms);
    acc = fold(acc, outcome.report.max_abs_error);
    acc = fold_u64(
        acc,
        u64::try_from(outcome.report.iterations_run).expect("iters"),
    );

    // Stage 2: Tustin discretize with prewarp at the top resonance.
    let t_s = 1.0 / 48_000.0;
    let prewarp = 2.0 * core::f64::consts::PI * 3_400.0;
    let filt = bilinear(&outcome.model, t_s, prewarp).expect("bilinear");
    assert!(filt.is_stable(), "discretized sections must be stable");
    let mut state = filt.zero_state();
    for k in 0..480 {
        let input = 0.5 * det::sin(core::f64::consts::TAU * 220.0 * k as f64 * t_s);
        let out = filt.step(&mut state, input).expect("step");
        acc = fold(acc, out);
    }

    println!(
        "{{\"suite\":\"fs-vfit\",\"case\":\"cross-isa-fit-filter\",\"arch\":\"{}\",\
         \"profile\":\"{}\",\"digest\":\"{acc:#018x}\",\"verdict\":\"golden-check\"}}",
        std::env::consts::ARCH,
        if cfg!(debug_assertions) {
            "debug"
        } else {
            "release"
        },
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "fit/filter bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — cross-ISA golden \
         event: bisect stage-wise, name the hazard, route through det:: in the same commit"
    );
}
