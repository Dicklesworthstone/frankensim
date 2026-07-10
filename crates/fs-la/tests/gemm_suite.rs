//! GEMM integration gates (bead xlvx): row-band parallel GEMM bitwise
//! equality with the serial kernel across thread counts.

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
