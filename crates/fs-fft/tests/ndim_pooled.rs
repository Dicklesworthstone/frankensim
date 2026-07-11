//! Executor-tiled N-D FFT gates (bead 27d3): pooled == serial BITWISE
//! at every worker count (the P2 law — parallel placement changes
//! timing, never bits), plus the cancellation contract.

use fs_exec::{CancelGate, PoolConfig, TilePool};
use fs_fft::{C64, FftNd};

fn fixture(total: usize) -> Vec<C64> {
    let mut seed = 0xD1D5_u64;
    (0..total)
        .map(|_| {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let re = ((seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5;
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let im = ((seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5;
            C64::new(re, im)
        })
        .collect()
}

fn bits_equal(a: &[C64], b: &[C64]) -> bool {
    a.iter()
        .zip(b)
        .all(|(x, y)| x.re.to_bits() == y.re.to_bits() && x.im.to_bits() == y.im.to_bits())
}

#[test]
fn pooled_ndim_is_bitwise_across_worker_counts() {
    for dims in [
        vec![8usize, 16],
        vec![4, 8, 2],
        vec![2, 2, 2, 4],
        vec![1, 16, 1, 4],
        vec![32],
    ] {
        let plan = FftNd::new(&dims);
        let x0 = fixture(plan.total());
        // Serial reference: forward then inverse.
        let mut fwd_ref = x0.clone();
        plan.forward(&mut fwd_ref);
        let mut inv_ref = fwd_ref.clone();
        plan.inverse(&mut inv_ref);
        for workers in [1usize, 2, 3, 7] {
            let pool = TilePool::new(PoolConfig::for_host(workers, 0xFD1D));
            let gate = CancelGate::new();
            let mut fwd = x0.clone();
            plan.forward_pooled(&mut fwd, &pool, &gate)
                .expect("pooled forward runs");
            assert!(
                bits_equal(&fwd, &fwd_ref),
                "forward bits drift: dims {dims:?} workers {workers}"
            );
            let mut inv = fwd.clone();
            plan.inverse_pooled(&mut inv, &pool, &gate)
                .expect("pooled inverse runs");
            assert!(
                bits_equal(&inv, &inv_ref),
                "inverse bits drift: dims {dims:?} workers {workers}"
            );
        }
    }
    println!(
        "{{\"suite\":\"fs-fft\",\"case\":\"ndim-pooled-bitwise\",\"verdict\":\"pass\",\
         \"detail\":\"pooled == serial bitwise over 5 shapes x 4 worker counts, forward+inverse\"}}"
    );
}

#[test]
fn pooled_ndim_cancellation_is_structured() {
    let plan = FftNd::new(&[8, 16]);
    let mut data = fixture(plan.total());
    let pool = TilePool::new(PoolConfig::for_host(2, 0xFD1D));
    let gate = CancelGate::new();
    gate.request(); // pre-cancelled: trips at the first bounded check
    let err = plan
        .forward_pooled(&mut data, &pool, &gate)
        .expect_err("pre-requested gate cancels");
    assert!(
        matches!(err, fs_exec::RunError::Cancelled { .. }),
        "expected structured cancellation, got {err:?}"
    );
    println!(
        "{{\"suite\":\"fs-fft\",\"case\":\"ndim-pooled-cancel\",\"verdict\":\"pass\",\
         \"detail\":\"pre-requested gate yields structured RunError::Cancelled; buffer contract documented\"}}"
    );
}
