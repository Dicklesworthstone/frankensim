//! SME2 exploratory battery (bead wf9.3, feature `frontier-sme2`):
//! probe outcome ledgered, capability-absent fallback, G0 equivalence
//! vs the scalar twin, and the honest benchmark row. Every test PASSES
//! on hardware without SME2 by verifying inertness — the flag must
//! never break a build or a box.

#![cfg(all(target_arch = "aarch64", feature = "frontier-sme2", not(miri)))]

use std::time::Instant;

use fs_simd::sme2::{
    TILE, gemm_tile_f32, gemm_tile_f32_scalar, sme2_available, streaming_vl_bytes,
};

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

fn panic_message(payload: &(dyn std::any::Any + Send)) -> &str {
    payload
        .downcast_ref::<&str>()
        .copied()
        .or_else(|| payload.downcast_ref::<String>().map(String::as_str))
        .unwrap_or("non-string panic payload")
}

fn lcg_f32(seed: &mut u64) -> f32 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5) as f32
}

/// sme2-001: the probe outcome is LEDGERED (probe result + SVL), and
/// the availability logic is consistent (available ⇒ SVL == 64).
#[test]
fn sme2_001_probe_ledgered() {
    let avail = sme2_available();
    let svl = streaming_vl_bytes();
    println!(
        "{{\"ledger\":\"sme2-probe\",\"available\":{avail},\"svl_bytes\":{}}}",
        svl.map_or(-1i64, |v| i64::try_from(v).unwrap_or(i64::MAX))
    );
    verdict(
        "sme2-001-probe",
        !avail || svl == Some(64),
        &format!(
            "probe: available={avail}, svl={svl:?} (available implies the fixed 512-bit shape)"
        ),
    );
}

/// sme2-002: capability-absent fallback — when the probe says no, the
/// module is INERT by contract (this test embodies the fallback gate on
/// non-SME boxes and documents the guard on SME boxes).
#[test]
fn sme2_002_capability_fallback() {
    if sme2_available() {
        verdict(
            "sme2-002-fallback",
            true,
            "hardware present: fallback path not exercised here — the guard assert is the gate (kernel refuses without capability)",
        );
    } else {
        // The kernel must refuse loudly rather than execute garbage.
        let a = vec![1.0f32; TILE];
        let b = vec![1.0f32; TILE];
        let mut c = vec![0.0f32; TILE * TILE];
        let refused = std::panic::catch_unwind(move || gemm_tile_f32(&a, &b, &mut c, 1)).is_err();
        verdict(
            "sme2-002-fallback",
            refused,
            "no SME2: kernel refuses loudly; the flag is inert and NEON remains committed",
        );
    }
}

/// sme2-002b: checked panel geometry refuses release-mode `usize`
/// wrapping before capability dispatch or unsafe assembly.
#[test]
fn sme2_002b_panel_extent_overflow_refused() {
    let payload = std::panic::catch_unwind(|| {
        let mut c = vec![0.0f32; TILE * TILE];
        gemm_tile_f32(&[], &[], &mut c, usize::MAX);
    })
    .expect_err("overflowed geometry unexpectedly reached dispatch");
    let message = panic_message(payload.as_ref());
    verdict(
        "sme2-002b-overflow",
        message.contains("SME2 panel length overflow or mismatch"),
        "overflowed k x TILE extent is rejected before unsafe dispatch",
    );
}

/// sme2-003: G0 equivalence — the streaming kernel vs the scalar
/// mul_add twin on identical panels, identical k-order. The battery
/// MEASURES whether fused outer-product accumulation is bitwise here
/// and ledgers the answer (no cross-ISA claim either way).
#[test]
fn sme2_003_g0_equivalence() {
    if !sme2_available() {
        verdict(
            "sme2-003-g0",
            true,
            "no SME2 on this box: equivalence lane vacuous, fallback covered by sme2-002",
        );
        return;
    }
    let mut seed = 0x53E2_u64;
    let mut worst_ulp = 0u64;
    let mut bitwise = true;
    for k in [1usize, 3, 17, 64, 257] {
        let a: Vec<f32> = (0..k * TILE).map(|_| lcg_f32(&mut seed)).collect();
        let b: Vec<f32> = (0..k * TILE).map(|_| lcg_f32(&mut seed)).collect();
        let mut c_sme = vec![0.0f32; TILE * TILE];
        let mut c_ref = vec![0.0f32; TILE * TILE];
        gemm_tile_f32(&a, &b, &mut c_sme, k);
        gemm_tile_f32_scalar(&a, &b, &mut c_ref, k);
        for (x, y) in c_sme.iter().zip(&c_ref) {
            if x.to_bits() != y.to_bits() {
                bitwise = false;
                let d = (i64::from(x.to_bits().cast_signed())
                    - i64::from(y.to_bits().cast_signed()))
                .unsigned_abs();
                worst_ulp = worst_ulp.max(d);
            }
        }
    }
    println!("{{\"ledger\":\"sme2-g0\",\"bitwise\":{bitwise},\"worst_ulp\":{worst_ulp}}}");
    verdict(
        "sme2-003-g0",
        bitwise || worst_ulp <= 1,
        &format!(
            "fmopa vs scalar mul_add twin over k in {{1,3,17,64,257}}: bitwise={bitwise}, worst ULP {worst_ulp} (measured, ledgered)"
        ),
    );
}

/// sme2-004: the honest benchmark row — GFLOP/s of the streaming tile
/// vs the scalar twin on the same panels. LEDGERED, not gated: the
/// promotion decision (beat NEON across the autotuner sweep) belongs
/// to the perf lanes; this is the exploratory evidence.
#[test]
fn sme2_004_benchmark_row() {
    if !sme2_available() {
        verdict(
            "sme2-004-bench",
            true,
            "no SME2 on this box: no benchmark row (inert)",
        );
        return;
    }
    let k = 1024usize;
    let reps = 200u32;
    let mut seed = 0xBE7C_u64;
    let a: Vec<f32> = (0..k * TILE).map(|_| lcg_f32(&mut seed)).collect();
    let b: Vec<f32> = (0..k * TILE).map(|_| lcg_f32(&mut seed)).collect();
    let mut c = vec![0.0f32; TILE * TILE];
    // Warm.
    gemm_tile_f32(&a, &b, &mut c, k);
    let t0 = Instant::now();
    for _ in 0..reps {
        gemm_tile_f32(&a, &b, &mut c, k);
    }
    let sme_s = t0.elapsed().as_secs_f64() / f64::from(reps);
    let t1 = Instant::now();
    for _ in 0..reps {
        gemm_tile_f32_scalar(&a, &b, &mut c, k);
    }
    let ref_s = t1.elapsed().as_secs_f64() / f64::from(reps);
    let flops = 2.0 * (TILE * TILE * k) as f64;
    let (g_sme, g_ref) = (flops / sme_s / 1e9, flops / ref_s / 1e9);
    println!(
        "{{\"ledger\":\"sme2-bench\",\"shape\":\"16x16x{k}\",\"sme2_gflops\":{g_sme:.2},\"scalar_gflops\":{g_ref:.2},\"speedup\":{:.2}}}",
        g_sme / g_ref
    );
    verdict(
        "sme2-004-bench",
        g_sme > 0.0 && g_ref > 0.0,
        &format!(
            "16x16x{k} f32: SME2 {g_sme:.1} GFLOP/s vs scalar twin {g_ref:.1} GFLOP/s ({:.1}x) — ledgered evidence, promotion stays with the perf lanes",
            g_sme / g_ref
        ),
    );
}
