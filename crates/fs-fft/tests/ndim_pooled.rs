//! Executor-tiled N-D FFT gates (bead 27d3): pooled == serial BITWISE
//! at every worker count (the P2 law — parallel placement changes
//! timing, never bits), plus the cancellation contract.

use fs_exec::{CancelGate, LatencyLane, PoolConfig, TilePool};
use fs_fft::{C64, FftNd};

fn in_task<R>(f: impl FnOnce(asupersync::Cx) -> R + Send + 'static) -> R
where
    R: Send + 'static,
{
    let lane = LatencyLane::new(1).expect("one-worker latency lane");
    let root = lane
        .runtime()
        .handle()
        .spawn(async move { f(asupersync::Cx::current().expect("live task context")) });
    lane.block_on(root)
}

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
fn task_scoped_pooled_ndim_preserves_bits_and_refuses_cancelled_task() {
    let dims = [64_usize, 256];
    let plan = FftNd::new(&dims);
    let input = fixture(plan.total());
    let mut serial = input.clone();
    plan.forward(&mut serial);

    let scoped_input = input.clone();
    let scoped = in_task(move |task_cx| {
        let plan = FftNd::new(&dims);
        let pool = TilePool::new(PoolConfig::for_host(3, 0x5C0E_D));
        let gate = CancelGate::new();
        let mut output = scoped_input;
        plan.forward_pooled_scoped(&mut output, &pool, &task_cx, &gate)
            .expect("task-scoped forward runs");
        output
    });
    assert!(
        bits_equal(&scoped, &serial),
        "task-scoped pooled transform must preserve serial bits"
    );

    let serial_plan = FftNd::new(&[32]);
    let refused_input = fixture(serial_plan.total());
    let refused_expected = refused_input.clone();
    let (error, after) = in_task(move |task_cx| {
        let plan = FftNd::new(&[32]);
        let pool = TilePool::new(PoolConfig::for_host(2, 0x5C0E_D));
        let gate = CancelGate::new();
        task_cx.set_cancel_requested(true);
        let mut output = refused_input.clone();
        let error = plan
            .forward_pooled_scoped(&mut output, &pool, &task_cx, &gate)
            .expect_err("cancelled task refuses before the first pencil pass");
        (error, output)
    });
    assert!(matches!(
        error,
        fs_exec::RunError::Cancelled { completed: 0, .. }
    ));
    assert!(
        bits_equal(&after, &refused_expected),
        "task cancellation at a serial pencil boundary must leave the scratch buffer untouched"
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
    // dims [4,8,2] under the 3f6c work floor (4096 elems/tile): every
    // pass collapses to ONE tile — a 64-element problem is serial-sized
    // and micro-tiling it is exactly the measured pathology.
    assert_eq!(passes.len(), 3, "one report per non-unit axis");
    assert_eq!(passes[0].kernel, "fs-fft/ndim-pencil-column-v2");
    assert_eq!(
        (
            passes[0].axis,
            passes[0].n,
            passes[0].stride,
            passes[0].outer
        ),
        (0, 4, 16, 1)
    );
    assert_eq!((passes[0].tiles, passes[0].completed), (1, 1));
    assert_eq!(passes[1].kernel, "fs-fft/ndim-pencil-block-v2");
    assert_eq!(
        (
            passes[1].axis,
            passes[1].n,
            passes[1].stride,
            passes[1].outer
        ),
        (1, 8, 2, 4)
    );
    assert_eq!((passes[1].tiles, passes[1].completed), (1, 1));
    assert_eq!(passes[2].kernel, "fs-fft/ndim-pencil-block-v2");
    assert_eq!(
        (
            passes[2].axis,
            passes[2].n,
            passes[2].stride,
            passes[2].outer
        ),
        (2, 2, 1, 32)
    );
    assert_eq!((passes[2].tiles, passes[2].completed), (1, 1));
    assert!(passes.iter().all(|p| p.workers == 3 && !p.inverse));

    // A work-floor-clearing shape splits into real tiles: [64,256] at 3
    // workers gives the column pass groups of 64 pencils (4 tiles) and
    // the block pass groups of 16 blocks x 256 elems (4 tiles).
    let big = FftNd::new(&[64, 256]);
    let big0 = fixture(big.total());
    let mut big_serial = big0.clone();
    big.forward(&mut big_serial);
    let mut big_pooled = big0.clone();
    let mut big_passes = Vec::new();
    big.forward_pooled_observed(&mut big_pooled, &pool, &gate, &mut |report| {
        big_passes.push(report)
    })
    .expect("big observed forward");
    assert!(bits_equal(&big_pooled, &big_serial), "grouped tiling bits");
    assert_eq!(big_passes.len(), 2);
    assert_eq!(big_passes[0].kernel, "fs-fft/ndim-pencil-column-v2");
    assert_eq!((big_passes[0].tiles, big_passes[0].completed), (4, 4));
    assert_eq!(big_passes[1].kernel, "fs-fft/ndim-pencil-block-v2");
    assert_eq!((big_passes[1].tiles, big_passes[1].completed), (4, 4));

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
