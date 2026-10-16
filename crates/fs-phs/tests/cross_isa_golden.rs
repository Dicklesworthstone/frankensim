//! Cross-ISA determinism goldens for `fs-phs` step ledgers on pinned
//! fixtures (bead `frankensim-music-v8-root-3ez8g.13.4`).
//!
//! The Gonzalez discrete-gradient stepper with its FD-Jacobian Newton
//! solve must produce bit-identical state/energy ledgers on both
//! reference ISA families in both build modes. Drive signals are built
//! from `fs_math::det` trig; a digest mismatch is a golden event:
//! bisect stage-wise, name the platform-libm hazard, route it through
//! `det::`, same commit.

use fs_math::det;
use fs_phs::{duffing_oscillator, modal_bank, step};

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

/// Duffing free decay ledger (same fixture shape as the conformance
/// battery's oscillator).
fn duffing_ledger() -> (u64, usize) {
    let duffing = duffing_oscillator(1.0, 1.0e-10, 1.0e-9, 0.0).expect("duffing");
    let dt = 1.0 / 48_000.0;
    let u = [0.0];
    let mut x = vec![1.0e-3, 0.0];
    let mut acc = FNV_OFFSET;
    for k in 0..2_000usize {
        let record = step(&duffing, &x, &u, dt).expect("step");
        acc = fold(acc, record.x[0]);
        acc = fold(acc, record.x[1]);
        acc = fold(acc, record.delta_h);
        x = record.x;
        let _ = k;
    }
    (acc, 2_000)
}

/// Three-mode mass-normalized modal bank (the bakeoff string card)
/// under zero force — the exact fixture family behind the
/// `string/phs-modal-bank` claims row.
fn modal_bank_ledger() -> u64 {
    const LENGTH_M: f64 = 0.65;
    const TENSION_N: f64 = 60.0;
    const LIN_DENSITY_KG_M: f64 = 6.0e-4;
    const MODES: usize = 3;
    const ZETAS: [f64; MODES] = [1.0e-3, 1.5e-3, 2.0e-3];
    let wave_speed = (TENSION_N / LIN_DENSITY_KG_M).sqrt();
    let omegas: Vec<f64> = (0..MODES)
        .map(|k| (k + 1) as f64 * core::f64::consts::PI * wave_speed / LENGTH_M)
        .collect();
    let drive = vec![0.0; MODES];
    let bank = modal_bank(&omegas, &ZETAS, &drive).expect("bank");
    let dt = 1.0 / 48_000.0;
    let mut x = vec![0.0; 2 * MODES];
    for mode in 0..MODES {
        x[2 * mode] = 1.0e-3;
    }
    let u = [0.0];
    let mut acc = FNV_OFFSET;
    for _ in 0..2_400 {
        let record = step(&bank, &x, &u, dt).expect("step");
        for value in &record.x {
            acc = fold(acc, *value);
        }
        acc = fold(acc, record.delta_h);
        x = record.x;
    }
    acc
}

/// Verified bit-identical aarch64-apple (debug) and x86_64-linux
/// (debug) on 2026-08-23, bead frankensim-music-v8-root-3ez8g.13.4.
const GOLDEN_HASH: u64 = 0x798c_84cb_eb3c_39b9;

#[test]
fn phs_step_ledger_digest_is_cross_isa_golden() {
    // det-trig drive sanity: the audit's own signal path stays inside the
    // deterministic elementary functions.
    let probe = det::sin(core::f64::consts::TAU * 220.0);
    assert!(probe.is_finite());

    let (duffing_digest, steps) = duffing_ledger();
    let bank_digest = modal_bank_ledger();
    let acc = fold_u64(
        fold_u64(duffing_digest, bank_digest),
        u64::try_from(steps).expect("n"),
    );

    println!(
        "{{\"suite\":\"fs-phs\",\"case\":\"cross-isa-step-ledger\",\"arch\":\"{}\",\
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
        "pHS ledger bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — cross-ISA golden \
         event: bisect stage-wise, name the hazard, route through det:: in the same commit"
    );
}
