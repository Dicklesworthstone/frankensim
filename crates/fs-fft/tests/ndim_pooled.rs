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
fn observed_pass_reports_bind_geometry_and_do_not_move_bits() {
    let dims = [4usize, 8, 2];
    let plan = FftNd::new(&dims);
    let x0 = fixture(plan.total());
    let pool = TilePool::new(PoolConfig::for_host(3, 0xFD1D));
    let gate = CancelGate::new();
    let mut plain = x0.clone();
    plan.forward_pooled(&mut plain, &pool, &gate)
        .expect("plain forward");
    let mut observed = x0.clone();
    let mut passes = Vec::new();
    plan.forward_pooled_observed(&mut observed, &pool, &gate, &mut |report| {
        passes.push(report)
    })
    .expect("observed forward");
    assert!(
        bits_equal(&observed, &plain),
        "the observer must not move bits"
    );
    // dims [4,8,2]: axis 0 is the column kernel (outer 1, stride 16;
    // group = ceil(16/3) = 6 -> 3 tiles); axis 1 is the block kernel
    // over 4 outer blocks; axis 2 (stride 1, outer 32) over 32 blocks.
    assert_eq!(passes.len(), 3, "one report per non-unit axis");
    assert_eq!(passes[0].kernel, "fs-fft/ndim-pencil-column-v1");
    assert_eq!(
        (
            passes[0].axis,
            passes[0].n,
            passes[0].stride,
            passes[0].outer
        ),
        (0, 4, 16, 1)
    );
    assert_eq!((passes[0].tiles, passes[0].completed), (3, 3));
    assert_eq!(passes[1].kernel, "fs-fft/ndim-pencil-block-v1");
    assert_eq!(
        (
            passes[1].axis,
            passes[1].n,
            passes[1].stride,
            passes[1].outer
        ),
        (1, 8, 2, 4)
    );
    assert_eq!((passes[1].tiles, passes[1].completed), (4, 4));
    assert_eq!(passes[2].kernel, "fs-fft/ndim-pencil-block-v1");
    assert_eq!(
        (
            passes[2].axis,
            passes[2].n,
            passes[2].stride,
            passes[2].outer
        ),
        (2, 2, 1, 32)
    );
    assert_eq!((passes[2].tiles, passes[2].completed), (32, 32));
    assert!(passes.iter().all(|p| p.workers == 3 && !p.inverse));

    // Inverse direction tags its passes; bits match the SERIAL inverse
    // (the roundtrip is approximate, the P2 law is serial == pooled).
    let mut inv_ref = plain.clone();
    plan.inverse(&mut inv_ref);
    let mut inv = observed;
    let mut inv_passes = Vec::new();
    plan.inverse_pooled_observed(&mut inv, &pool, &gate, &mut |report| {
        inv_passes.push(report)
    })
    .expect("observed inverse");
    assert!(
        bits_equal(&inv, &inv_ref),
        "observed inverse matches the serial inverse bitwise"
    );
    assert_eq!(inv_passes.len(), 3);
    assert!(inv_passes.iter().all(|p| p.inverse));

    // A pre-cancelled run still observes the interrupted pass, with
    // nothing completed.
    let cancelled_gate = CancelGate::new();
    cancelled_gate.request();
    let mut cancelled = x0.clone();
    let mut cancelled_passes = Vec::new();
    let err = plan
        .forward_pooled_observed(&mut cancelled, &pool, &cancelled_gate, &mut |report| {
            cancelled_passes.push(report);
        })
        .expect_err("pre-requested gate cancels");
    assert!(matches!(err, fs_exec::RunError::Cancelled { .. }));
    assert_eq!(
        cancelled_passes.len(),
        1,
        "the interrupted pass is observed"
    );
    assert_eq!(cancelled_passes[0].completed, 0);
    println!(
        "{{\"suite\":\"fs-fft\",\"case\":\"ndim-pass-observation\",\"verdict\":\"pass\",\
         \"detail\":\"per-pass reports bind kernel/geometry/tiles, bits identical, \
         interrupted pass observed with completed 0\"}}"
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
