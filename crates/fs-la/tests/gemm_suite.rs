//! GEMM integration gates (bead xlvx): row-band parallel GEMM bitwise
//! equality with the serial kernel across thread counts.

use std::any::Any;

fn panic_text(payload: &(dyn Any + Send)) -> String {
    if let Some(text) = payload.downcast_ref::<String>() {
        text.clone()
    } else if let Some(text) = payload.downcast_ref::<&str>() {
        (*text).to_string()
    } else {
        "non-string panic".to_string()
    }
}

fn assert_extent_overflow(f: impl FnOnce() + std::panic::UnwindSafe, label: &str) {
    let panic = std::panic::catch_unwind(f).expect_err("overflow must panic");
    let text = panic_text(panic.as_ref());
    assert!(
        text.contains("extent overflow"),
        "{label}: unexpected panic: {text}"
    );
}

/// G0: safe GEMM facades reject impossible extents identically in debug and
/// release. In particular, validation happens before beta can mutate C.
#[test]
fn public_facades_reject_extent_overflow_before_mutation() {
    let mut c64 = [7.0f64];
    assert_extent_overflow(
        std::panic::AssertUnwindSafe(|| {
            fs_la::gemm_f64(usize::MAX, 0, 2, 1.0, &[], &[], 0.0, &mut c64);
        }),
        "gemm_f64",
    );
    assert_eq!(c64[0].to_bits(), 7.0f64.to_bits());

    let mut c32 = [7.0f32];
    assert_extent_overflow(
        std::panic::AssertUnwindSafe(|| {
            fs_la::gemm_f32(usize::MAX, 0, 2, 1.0, &[], &[], 0.0, &mut c32);
        }),
        "gemm_f32",
    );
    assert_eq!(c32[0].to_bits(), 7.0f32.to_bits());

    let mut mixed = [7.0f64];
    assert_extent_overflow(
        std::panic::AssertUnwindSafe(|| {
            fs_la::gemm_mixed(usize::MAX, 0, 2, 1.0, &[], &[], 0.0, &mut mixed);
        }),
        "gemm_mixed",
    );
    assert_eq!(mixed[0].to_bits(), 7.0f64.to_bits());

    let mut parallel = [7.0f64];
    assert_extent_overflow(
        std::panic::AssertUnwindSafe(|| {
            fs_la::gemm_f64_parallel(usize::MAX, 0, 2, 1.0, &[], &[], 0.0, &mut parallel, 4);
        }),
        "gemm_f64_parallel",
    );
    assert_eq!(parallel[0].to_bits(), 7.0f64.to_bits());

    let mut op = [7.0f64];
    assert_extent_overflow(
        std::panic::AssertUnwindSafe(|| {
            fs_la::gemm_f64_op(
                2,
                0,
                1,
                1.0,
                &[],
                usize::MAX,
                fs_la::Trans::N,
                &[],
                1,
                fs_la::Trans::N,
                0.0,
                &mut op,
                1,
            );
        }),
        "gemm_f64_op",
    );
    assert_eq!(op[0].to_bits(), 7.0f64.to_bits());
}

/// G0: tuner quanta need not align to the microkernel. Packed storage must
/// account for the padded MR/NR tails even when a full KC panel is active.
#[test]
fn parallel_unaligned_tune_quanta_size_padded_packs() {
    let (m, n, k) = (257usize, 5usize, 256usize);
    let a: Vec<f64> = (0..m * k).map(|i| (i % 31) as f64 - 15.0).collect();
    let b: Vec<f64> = (0..k * n).map(|i| (i % 17) as f64 - 8.0).collect();
    let mut serial = vec![0.25; m * n];
    let mut parallel = serial.clone();
    fs_la::gemm_f64(m, n, k, 0.75, &a, &b, 0.5, &mut serial);
    fs_la::gemm_f64_parallel_with(m, n, k, 0.75, &a, &b, 0.5, &mut parallel, 2, 60, 5);
    assert!(
        serial
            .iter()
            .zip(&parallel)
            .all(|(lhs, rhs)| lhs.to_bits() == rhs.to_bits())
    );
}

/// G0: the producer-owned routing query is the single source of truth for
/// whether a session may publish MC/NC timing evidence.
#[test]
fn tuning_effectiveness_refuses_noop_serial_and_small_routes() {
    assert!(!fs_la::gemm_tuning_is_effective(512, 512, 512, 1.0, 1));
    assert!(!fs_la::gemm_tuning_is_effective(255, 512, 512, 1.0, 8));
    assert!(!fs_la::gemm_tuning_is_effective(512, 0, 512, 1.0, 8));
    assert!(!fs_la::gemm_tuning_is_effective(512, 512, 0, 1.0, 8));
    assert!(!fs_la::gemm_tuning_is_effective(512, 512, 512, -0.0, 8));
    assert!(fs_la::gemm_tuning_is_effective(256, 1, 1, 1.0, 2));
    assert!(!fs_la::gemm_execution_tier().is_empty());
    assert_eq!(fs_la::GEMM_IMPLEMENTATION_VERSION, 1);
    assert_eq!(
        fs_la::GEMM_PARALLEL_IMPLEMENTATION,
        "std-thread-scope-work-stealing-unpinned-v1"
    );
}

/// xlvx item 3: row-band parallel GEMM is BITWISE equal to serial at
/// every thread count (per-element accumulation order is independent
/// of m — xdgf's recorded fact (b), now gated).
#[test]
fn parallel_gemm_bitwise_across_thread_counts() {
    // m >= 2*MC so the THREADED path runs (below that the facade
    // falls back to serial and the gate would test nothing); all three
    // dims deliberately unaligned to MR/NR/KC/MC.
    let (m, n, k) = (391usize, 173, 83);
    let a: Vec<f64> = (0..m * k).map(|i| ((i as f64) * 0.7).sin()).collect();
    let b: Vec<f64> = (0..k * n).map(|i| ((i as f64) * 1.3).cos()).collect();
    let mut c_ref: Vec<f64> = (0..m * n).map(|i| (i as f64) * 0.01 - 3.0).collect();
    let c0 = c_ref.clone();
    fs_la::gemm_f64(m, n, k, 1.25, &a, &b, 0.5, &mut c_ref);
    for t in [1usize, 2, 3, 5, 8, 16] {
        let mut c_par = c0.clone();
        fs_la::gemm_f64_parallel(m, n, k, 1.25, &a, &b, 0.5, &mut c_par, t);
        assert!(
            c_ref
                .iter()
                .zip(&c_par)
                .all(|(x, y)| x.to_bits() == y.to_bits()),
            "parallel gemm (t={t}) != serial bitwise"
        );
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"xlvx-parallel-bitwise\",\"verdict\":\"pass\",\"detail\":\"row-band parallel GEMM bitwise == serial for t in 1/2/3/5/8/16 on unaligned 391x173x83 (threaded path)\"}}"
    );
}

/// xlvx segment 5: MC/NC are BIT-NEUTRAL (the module-doc contract the
/// adaptive-MC parallel path and the autotune sweep both lean on) —
/// gate it empirically over an (mc, nc) grid, including MR/NR-unaligned
/// quanta and nc wider than n.
#[test]
fn parallel_gemm_bitwise_across_blockings() {
    let (m, n, k) = (391usize, 173, 83);
    let a: Vec<f64> = (0..m * k).map(|i| ((i as f64) * 0.7).sin()).collect();
    let b: Vec<f64> = (0..k * n).map(|i| ((i as f64) * 1.3).cos()).collect();
    let mut c_ref: Vec<f64> = (0..m * n).map(|i| (i as f64) * 0.01 - 3.0).collect();
    let c0 = c_ref.clone();
    fs_la::gemm_f64(m, n, k, 1.25, &a, &b, 0.5, &mut c_ref);
    for (mc, nc) in [
        (8usize, 4usize),
        (16, 64),
        (40, 128),
        (60, 100), // deliberately MR/NR-unaligned quanta
        (128, 512),
        (128, 1024), // nc > n
    ] {
        let mut c_par = c0.clone();
        fs_la::gemm_f64_parallel_with(m, n, k, 1.25, &a, &b, 0.5, &mut c_par, 7, mc, nc);
        assert!(
            c_ref
                .iter()
                .zip(&c_par)
                .all(|(x, y)| x.to_bits() == y.to_bits()),
            "parallel gemm (mc={mc}, nc={nc}) != serial bitwise"
        );
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"xlvx-blocking-bitwise\",\"verdict\":\"pass\",\"detail\":\"parallel GEMM bitwise == serial over (mc,nc) grid incl unaligned quanta on 391x173x83\"}}"
    );
}
