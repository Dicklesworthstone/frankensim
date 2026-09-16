//! Tests for `brownian_frames` / `brownian_frames_window`.
//! Bead: am-fs-export-brownian-frames-nhm.

use fs_rand::{STREAM_CHECKPOINT_CANONICAL_LEN, StreamCheckpoint, StreamKey};
use fs_wasm::brownian::{
    derived_stream_index, draws_per_step, quantity_binding, reconstruct_particle,
};
use fs_wasm::{
    BROWNIAN_MAX_OUTPUT_LEN, BROWNIAN_STREAM_KERNEL_ID, admit_brownian_checkpoint, brownian_frames,
    brownian_frames_window,
};

fn last_column(buf: &[f64], n: usize, steps: usize) -> Vec<f64> {
    let cols = steps + 1;
    (0..n).map(|p| buf[p * cols + steps]).collect()
}

fn concat_drop_first(
    first: &[f64],
    next: &[f64],
    n: usize,
    first_steps: usize,
    next_steps: usize,
) -> Vec<f64> {
    let total_steps = first_steps + next_steps;
    let mut out = vec![0.0; n * (total_steps + 1)];
    let fcols = first_steps + 1;
    let ncols = next_steps + 1;
    for p in 0..n {
        for s in 0..=first_steps {
            out[p * (total_steps + 1) + s] = first[p * fcols + s];
        }
        for s in 1..=next_steps {
            out[p * (total_steps + 1) + first_steps + s] = next[p * ncols + s];
        }
    }
    out
}

fn verdict(case: &str, step_kernel: u32, seed: u64, outcome: &str) {
    eprintln!(
        "{{\"suite\":\"brownian_frames\",\"beadId\":\"am-fs-export-brownian-frames-nhm\",\"caseId\":\"{case}\",\"stepKernel\":{step_kernel},\"seed\":{seed},\"outcome\":\"{outcome}\"}}"
    );
}

#[test]
fn layout_and_initial_positions() {
    for &(n, steps) in &[(1usize, 1usize), (3, 4), (7, 8)] {
        let buf = brownian_frames(n, steps, 0, 1, 1.0, 1.0).expect("valid");
        assert_eq!(buf.len(), n * (steps + 1));
        for p in 0..n {
            assert_eq!(buf[p * (steps + 1)], 0.0);
        }
    }
    verdict("layout_and_initial_positions", 0, 1, "pass");
}

#[test]
fn per_kernel_step_rules() {
    let d = 1.25e-13;
    let dt = 0.1;
    let s = (2.0 * d * dt).sqrt();
    let h = (6.0 * d * dt).sqrt();
    let n = 400usize;
    let buf0 = brownian_frames(n, 1, 0, 7, d, dt).unwrap();
    let mut m2 = 0.0;
    let mut m4 = 0.0;
    for p in 0..n {
        let step = buf0[p * 2 + 1];
        assert!(step == s || step == -s, "coin step must be ±s");
        m2 += step * step;
        m4 += step.powi(4);
    }
    m2 /= n as f64;
    m4 /= n as f64;
    let sigma2 = 2.0 * d * dt;
    assert!(
        (m2 - sigma2).abs()
            <= 4.0 * f64::EPSILON * s * s.max(1.0) + s * s / (n as f64).sqrt() * 6.0
    );
    assert!((m4 / sigma2.powi(2) - 1.0).abs() < 0.25);

    let buf1 = brownian_frames(n, 1, 1, 7, d, dt).unwrap();
    m2 = 0.0;
    m4 = 0.0;
    for p in 0..n {
        let step = buf1[p * 2 + 1];
        assert!(step >= -h && step < h, "uniform support [-h, h)");
        m2 += step * step;
        m4 += step.powi(4);
    }
    m2 /= n as f64;
    m4 /= n as f64;
    let se = ((1.8 * sigma2.powi(2) - sigma2.powi(2)) / n as f64).sqrt();
    assert!((m2 - sigma2).abs() < 6.0 * se.max(1e-30));
    assert!((m4 / sigma2.powi(2) - 1.8).abs() < 0.6);

    let buf3 = brownian_frames(n, 1, 3, 7, d, dt).unwrap();
    m2 = 0.0;
    m4 = 0.0;
    for p in 0..n {
        let step = buf3[p * 2 + 1];
        m2 += step * step;
        m4 += step.powi(4);
    }
    m2 /= n as f64;
    m4 /= n as f64;
    let se_g = ((3.0 * sigma2.powi(2) - sigma2.powi(2)) / n as f64).sqrt();
    assert!((m2 - sigma2).abs() < 6.0 * se_g.max(1e-30));
    assert!((m4 / sigma2.powi(2) - 3.0).abs() < 1.5);
    verdict("per_kernel_step_rules", 3, 7, "pass");
}

#[test]
fn draw_consumption() {
    assert_eq!(draws_per_step(0), Some(1));
    assert_eq!(draws_per_step(1), Some(1));
    assert_eq!(draws_per_step(2), Some(2));
    assert_eq!(draws_per_step(3), Some(2));
    assert_eq!(derived_stream_index(0, 0, 10), Some(10));
    assert_eq!(derived_stream_index(3, 5, 7), Some(24));
    verdict("draw_consumption", 3, 0, "pass");
}

#[test]
fn arithmetic_reconstruction() {
    let d = 2.0;
    let dt = 0.25;
    let steps = 16;
    let seed = 99u64;
    for kernel in [0u32, 1] {
        let buf = brownian_frames(5, steps, kernel, seed, d, dt).unwrap();
        for p in 0..5u32 {
            let recon = reconstruct_particle(p, steps, kernel, seed, d, dt).unwrap();
            let cols = steps + 1;
            for s in 0..=steps {
                assert_eq!(
                    buf[p as usize * cols + s].to_bits(),
                    recon[s].to_bits(),
                    "kernel {kernel} particle {p} step {s}"
                );
            }
        }
    }
    verdict("arithmetic_reconstruction", 0, 99, "pass");
}

#[test]
fn kernel3_mean_square_growth() {
    let d = 1.0;
    let dt = 1.0;
    let n = 400usize;
    let steps = 8usize;
    let buf = brownian_frames(n, steps, 3, 123, d, dt).unwrap();
    let cols = steps + 1;
    for s in [1usize, 4, 8] {
        let mut ss = 0.0;
        for p in 0..n {
            let x = buf[p * cols + s];
            ss += x * x;
        }
        let expected = 2.0 * d * (s as f64) * dt;
        let ratio = (ss / n as f64) / expected;
        assert!(
            (0.7..=1.3).contains(&ratio),
            "step {s} mean-square ratio {ratio}"
        );
    }
    verdict("kernel3_mean_square_growth", 3, 123, "pass");
}

#[test]
fn kernel2_resolution() {
    assert_eq!(quantity_binding(2), Some(("walkStepCoordinate1d", "step")));
    assert_eq!(quantity_binding(3), Some(("latentPosition1d", "metre")));
    let a = brownian_frames(4, 8, 2, 5, 1.0, 1.0).unwrap();
    let b = brownian_frames(4, 8, 2, 5, 9.0, 0.01).unwrap();
    assert_eq!(a, b, "kernel 2 must ignore D and dt");
    let c = brownian_frames(4, 8, 3, 5, 1.0, 1.0).unwrap();
    assert_ne!(a, c, "kernel 2 is not an alias of kernel 3");
    verdict("kernel2_resolution", 2, 5, "pass");
}

#[test]
fn prefix_stability_and_regeneration() {
    let small = brownian_frames(3, 12, 1, 42, 1.0, 0.5).unwrap();
    let large = brownian_frames(100, 12, 1, 42, 1.0, 0.5).unwrap();
    let cols = 13;
    for p in 0..3 {
        for s in 0..cols {
            assert_eq!(small[p * cols + s].to_bits(), large[p * cols + s].to_bits());
        }
    }
    let one = brownian_frames(1, 12, 1, 42, 1.0, 0.5).unwrap();
    for s in 0..cols {
        assert_eq!(one[s].to_bits(), large[s].to_bits());
    }
    verdict("prefix_stability_and_regeneration", 1, 42, "pass");
}

#[test]
fn seed_boundaries() {
    let a = brownian_frames(2, 4, 0, 0, 1.0, 1.0).unwrap();
    let b = brownian_frames(2, 4, 0, 9_007_199_254_740_992, 1.0, 1.0).unwrap();
    let c = brownian_frames(2, 4, 0, 9_007_199_254_740_993, 1.0, 1.0).unwrap();
    let d = brownian_frames(2, 4, 0, u64::MAX, 1.0, 1.0).unwrap();
    assert_ne!(b, c);
    assert_ne!(a, b);
    assert_ne!(a, d);
    verdict("seed_boundaries", 0, 0, "pass");
}

#[test]
fn refusals_not_clamps() {
    let bad_k = brownian_frames(2, 2, 4, 1, 1.0, 1.0).unwrap_err();
    assert_eq!(bad_k.code, "unsupported-kernel");
    let nan = brownian_frames(2, 2, 0, 1, f64::NAN, 1.0).unwrap_err();
    assert_eq!(nan.code, "nonfinite-input");
    let neg = brownian_frames(2, 2, 0, 1, -1.0, 1.0).unwrap_err();
    assert_eq!(neg.code, "invalid-parameter");
    let zero_n = brownian_frames(0, 2, 0, 1, 1.0, 1.0).unwrap_err();
    assert_eq!(zero_n.code, "invalid-parameter");
    let zero_s = brownian_frames(2, 0, 0, 1, 1.0, 1.0).unwrap_err();
    assert_eq!(zero_s.code, "invalid-parameter");
    let zeros = brownian_frames(2, 3, 0, 1, 0.0, 1.0).unwrap();
    assert!(zeros.iter().all(|&x| x == 0.0));
    let over = brownian_frames(BROWNIAN_MAX_OUTPUT_LEN, 1, 0, 1, 1.0, 1.0).unwrap_err();
    assert_eq!(over.code, "budget-exhausted");
    verdict("refusals_not_clamps", 0, 1, "pass");
}

#[test]
fn window_concatenation_is_bitwise() {
    let n = 7usize;
    let seed = 11u64;
    let d = 0.5;
    let dt = 0.25;
    for kernel in [0u32, 1, 3] {
        let one = brownian_frames(n, 64, kernel, seed, d, dt).unwrap();
        let splits: &[&[usize]] = &[&[64], &[32, 32], &[1, 63], &[16, 16, 16, 16]];
        for split in splits {
            let mut start_step = 0usize;
            let mut starts = vec![0.0; n];
            let mut pieces: Vec<(usize, Vec<f64>)> = Vec::new();
            for &w in *split {
                let buf =
                    brownian_frames_window(n, start_step, w, kernel, seed, d, dt, &starts).unwrap();
                starts = last_column(&buf, n, w);
                pieces.push((w, buf));
                start_step += w;
            }
            let mut acc = pieces[0].1.clone();
            let mut acc_steps = pieces[0].0;
            for (w, buf) in pieces.iter().skip(1) {
                acc = concat_drop_first(&acc, buf, n, acc_steps, *w);
                acc_steps += *w;
            }
            assert_eq!(acc.len(), one.len());
            for i in 0..acc.len() {
                assert_eq!(
                    acc[i].to_bits(),
                    one[i].to_bits(),
                    "kernel {kernel} idx {i}"
                );
            }
        }
    }
    verdict("window_concatenation_is_bitwise", 3, 11, "pass");
}

#[test]
fn window_start_position_validation() {
    let err_len = brownian_frames_window(2, 1, 2, 0, 1, 1.0, 1.0, &[0.0]).unwrap_err();
    assert_eq!(err_len.code, "invalid-parameter");
    let err_nan = brownian_frames_window(1, 1, 2, 0, 1, 1.0, 1.0, &[f64::NAN]).unwrap_err();
    assert_eq!(err_nan.code, "invalid-parameter");
    let err_inf = brownian_frames_window(1, 1, 2, 0, 1, 1.0, 1.0, &[f64::INFINITY]).unwrap_err();
    assert_eq!(err_inf.code, "invalid-parameter");
    let err_zero = brownian_frames_window(1, 0, 2, 0, 1, 1.0, 1.0, &[1.0]).unwrap_err();
    assert_eq!(err_zero.code, "invalid-parameter");
    verdict("window_start_position_validation", 0, 1, "pass");
}

#[test]
fn window_stream_overflow() {
    let start = (u64::MAX / 2) as usize + 1;
    let err = brownian_frames_window(1, start, 1, 3, 1, 1.0, 1.0, &[0.0]).unwrap_err();
    assert_eq!(err.code, "stream-index-overflow");
    assert!(err.details.contains("startIndex"));
    assert!(err.details.contains("draws"));
    assert!(err.details.contains("maxIndex"));
    verdict("window_stream_overflow", 3, 1, "pass");
}

#[test]
fn window_checkpoint_disagreement() {
    let good = StreamCheckpoint::current(
        StreamKey {
            seed: 1,
            kernel: BROWNIAN_STREAM_KERNEL_ID,
            tile: 0,
        },
        10,
    )
    .to_canonical_le_bytes();
    admit_brownian_checkpoint(1, 0, 10, 0, &good).expect("matching checkpoint");

    let wrong_tile = StreamCheckpoint::current(
        StreamKey {
            seed: 1,
            kernel: BROWNIAN_STREAM_KERNEL_ID,
            tile: 7,
        },
        10,
    )
    .to_canonical_le_bytes();
    let e = admit_brownian_checkpoint(1, 0, 10, 0, &wrong_tile).unwrap_err();
    assert_eq!(e.code, "invalid-parameter");
    assert!(e.details.contains("tile"));

    let mut bad_sem = StreamCheckpoint::current(
        StreamKey {
            seed: 1,
            kernel: BROWNIAN_STREAM_KERNEL_ID,
            tile: 0,
        },
        10,
    );
    bad_sem.stream_semantics_version = 9;
    let e = admit_brownian_checkpoint(1, 0, 10, 0, &bad_sem.to_canonical_le_bytes()).unwrap_err();
    assert_eq!(e.code, "invalid-parameter");
    assert!(e.details.contains("stream_semantics_version"));

    let wrong_idx = StreamCheckpoint::current(
        StreamKey {
            seed: 1,
            kernel: BROWNIAN_STREAM_KERNEL_ID,
            tile: 0,
        },
        11,
    )
    .to_canonical_le_bytes();
    let e = admit_brownian_checkpoint(1, 0, 10, 0, &wrong_idx).unwrap_err();
    assert!(e.details.contains("index"));
    assert_eq!(good.len(), STREAM_CHECKPOINT_CANONICAL_LEN);
    verdict("window_checkpoint_disagreement", 0, 1, "pass");
}

#[test]
fn one_shot_is_a_window_from_zeros() {
    let a = brownian_frames(5, 9, 3, 8, 0.4, 0.2).unwrap();
    let b = brownian_frames_window(5, 0, 9, 3, 8, 0.4, 0.2, &[0.0; 5]).unwrap();
    assert_eq!(a, b);
    verdict("one_shot_is_a_window_from_zeros", 3, 8, "pass");
}
