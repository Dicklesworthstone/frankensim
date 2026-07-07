//! fs-tilelang battery (wf9.11): the three reference kernels of the
//! acceptance criteria (batched axpy, stencil apply, SDF-style
//! trilinear grid sample) plus a deterministic-sum reduction kernel —
//! each written ONCE, lowered to scalar + lane variants, checked
//! against hand-written oracles, tier-equivalent bitwise (G0),
//! deterministic on repeat (G5), with intensity metadata logged for
//! the roofline harness. The macro's auto-generated twin tests run
//! alongside (visible as `__twin_tests` in the test list).

use fs_rand::StreamKey;
use fs_tilelang::{DeterminismClass, ReductionKind, kernel};

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-tilelang\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

fn rand_vec(n: usize, tile: u32) -> Vec<f64> {
    let mut s = StreamKey {
        seed: 77,
        kernel: 0x71E5,
        tile,
    }
    .stream();
    (0..n).map(|_| 2.0f64.mul_add(s.next_f64(), -1.0)).collect()
}

// Reference kernel 1: batched axpy-like (pure elementwise map).
kernel! {
    name: axpy_k,
    reads: [x, y],
    params: [alpha],
    writes: [out],
    reduction: none,
    body: {
        out = alpha.mul_add(x, y);
    },
}

// Reference kernel 2a: 1D 3-point stencil (literal shifts, halo 1) —
// gets the auto twin tests.
kernel! {
    name: stencil3_k,
    reads: [u],
    params: [c0, c1],
    writes: [out],
    halo: 1,
    reduction: none,
    body: {
        out = c0.mul_add(u, c1 * (shift_sub(u, 1) + shift_add(u, 1)));
    },
}

// Reference kernel 2b: 3D 7-point stencil via stride uparams (halo =
// one xy-plane); uparam kernels drive their own twin checks here.
kernel! {
    name: stencil7_k,
    reads: [u],
    uparams: [nx, nxy],
    params: [c0, c1],
    writes: [out],
    halo: nxy,
    reduction: none,
    body: {
        let ring = shift_sub(u, 1) + shift_add(u, 1)
            + shift_sub(u, nx) + shift_add(u, nx)
            + shift_sub(u, nxy) + shift_add(u, nxy);
        out = c0.mul_add(u, c1 * ring);
    },
}

// Reference kernel 3: SDF-style trilinear grid sample (gather form).
kernel! {
    name: trilinear_k,
    reads: [g, fx, fy, fz],
    index_reads: [ix, iy, iz],
    uparams: [nx, nxy],
    writes: [out],
    reduction: none,
    body: {
        let base = ix + nx * iy + nxy * iz;
        let c00 = gather(g, base) * (1.0 - fx) + gather(g, base + 1) * fx;
        let c10 = gather(g, base + nx) * (1.0 - fx) + gather(g, base + nx + 1) * fx;
        let c01 = gather(g, base + nxy) * (1.0 - fx) + gather(g, base + nxy + 1) * fx;
        let c11 = gather(g, base + nxy + nx) * (1.0 - fx) + gather(g, base + nxy + nx + 1) * fx;
        let c0 = c00 * (1.0 - fy) + c10 * fy;
        let c1 = c01 * (1.0 - fy) + c11 * fy;
        out = c0 * (1.0 - fz) + c1 * fz;
    },
}

// Reduction kernel: deterministic dot product.
kernel! {
    name: dot_k,
    reads: [x, y],
    writes: [],
    reduction: deterministic_sum,
    body: {
        acc = x * y;
    },
}

#[test]
fn axpy_matches_oracle_and_meta() {
    let n = 1537;
    let (x, y) = (rand_vec(n, 1), rand_vec(n, 2));
    let mut out = vec![0.0f64; n];
    axpy_k::run(&x, &y, 1.75, &mut out);
    for i in [0usize, n / 2, n - 1] {
        assert_eq!(
            out[i].to_bits(),
            1.75f64.mul_add(x[i], y[i]).to_bits(),
            "axpy oracle mismatch at {i}"
        );
    }
    assert_eq!(axpy_k::META.flops_per_elem, 2);
    assert_eq!(axpy_k::META.bytes_per_elem, 24);
    assert_eq!(axpy_k::META.reduction, ReductionKind::None);
    assert_eq!(axpy_k::META.determinism, DeterminismClass::BitwiseAllTiers);
    log("axpy", "pass", &axpy_k::META.descr());
}

#[test]
fn stencil3_matches_oracle_and_halo() {
    let n = 513;
    let u = rand_vec(n, 3);
    let mut out = vec![f64::NAN; n];
    stencil3_k::run(&u, 0.5, 0.25, &mut out);
    // Halo untouched (still NaN), interior matches the oracle.
    assert!(
        out[0].is_nan() && out[n - 1].is_nan(),
        "halo must be untouched"
    );
    for i in 1..n - 1 {
        let expect = 0.5f64.mul_add(u[i], 0.25 * (u[i - 1] + u[i + 1]));
        assert_eq!(
            out[i].to_bits(),
            expect.to_bits(),
            "stencil3 mismatch at {i}"
        );
    }
    assert_eq!(stencil3_k::META.halo, 1);
    log("stencil3", "pass", &stencil3_k::META.descr());
}

#[test]
fn stencil7_matches_oracle_and_tier_twins() {
    // 3D grid flattened: nx=12, ny=11, nz=9. Halo of one plane means
    // some non-interior cells get written too (wrap-reads within
    // bounds) — the ORACLE uses identical index arithmetic, so
    // equality is exact everywhere the kernel writes.
    let (nx, ny, nz) = (12usize, 11, 9);
    let nxy = nx * ny;
    let n = nxy * nz;
    let u = rand_vec(n, 4);
    let mut out = vec![0.0f64; n];
    stencil7_k::run(&u, nx, nxy, 0.4, 0.1, &mut out);
    let mut worst = 0u64;
    for i in nxy..n - nxy {
        let ring = u[i - 1] + u[i + 1] + u[i - nx] + u[i + nx] + u[i - nxy] + u[i + nxy];
        let expect = 0.4f64.mul_add(u[i], 0.1 * ring);
        worst = worst.max(out[i].to_bits() ^ expect.to_bits());
    }
    assert_eq!(worst, 0, "stencil7 oracle mismatch");
    // uparam kernels drive their own tier twins (macro can't guess
    // strides): all lane widths bitwise-equal to scalar.
    let mut out_s = vec![0.0f64; n];
    stencil7_k::run_scalar(&u, nx, nxy, 0.4, 0.1, &mut out_s);
    for w in [2usize, 4, 8] {
        let mut out_w = vec![0.0f64; n];
        match w {
            2 => stencil7_k::run_lanes::<2>(&u, nx, nxy, 0.4, 0.1, &mut out_w),
            4 => stencil7_k::run_lanes::<4>(&u, nx, nxy, 0.4, 0.1, &mut out_w),
            _ => stencil7_k::run_lanes::<8>(&u, nx, nxy, 0.4, 0.1, &mut out_w),
        }
        assert!(
            out_s
                .iter()
                .zip(&out_w)
                .all(|(a, b)| a.to_bits() == b.to_bits()),
            "stencil7 lane width {w} diverges from scalar"
        );
    }
    log("stencil7", "pass", &stencil7_k::META.descr());
}

#[test]
fn trilinear_matches_oracle_and_tier_twins() {
    // Grid 8×7×6, 500 query points with in-range bases.
    let (nx, ny, nz) = (8usize, 7, 6);
    let nxy = nx * ny;
    let g = rand_vec(nxy * nz, 5);
    let m = 500usize;
    let mut s = StreamKey {
        seed: 78,
        kernel: 0x71E5,
        tile: 6,
    }
    .stream();
    let ix: Vec<u32> = (0..m)
        .map(|_| u32::try_from(s.next_below(nx as u64 - 1)).expect("small"))
        .collect();
    let iy: Vec<u32> = (0..m)
        .map(|_| u32::try_from(s.next_below(ny as u64 - 1)).expect("small"))
        .collect();
    let iz: Vec<u32> = (0..m)
        .map(|_| u32::try_from(s.next_below(nz as u64 - 1)).expect("small"))
        .collect();
    let fx: Vec<f64> = (0..m).map(|_| s.next_f64()).collect();
    let fy: Vec<f64> = (0..m).map(|_| s.next_f64()).collect();
    let fz: Vec<f64> = (0..m).map(|_| s.next_f64()).collect();
    let mut out = vec![0.0f64; m];
    trilinear_k::run(&g, &fx, &fy, &fz, &ix, &iy, &iz, nx, nxy, &mut out);
    // Hand-written oracle with identical arithmetic order.
    for q in [0usize, 250, 499] {
        let base = ix[q] as usize + nx * (iy[q] as usize) + nxy * (iz[q] as usize);
        let c00 = g[base] * (1.0 - fx[q]) + g[base + 1] * fx[q];
        let c10 = g[base + nx] * (1.0 - fx[q]) + g[base + nx + 1] * fx[q];
        let c01 = g[base + nxy] * (1.0 - fx[q]) + g[base + nxy + 1] * fx[q];
        let c11 = g[base + nxy + nx] * (1.0 - fx[q]) + g[base + nxy + nx + 1] * fx[q];
        let c0 = c00 * (1.0 - fy[q]) + c10 * fy[q];
        let c1 = c01 * (1.0 - fy[q]) + c11 * fy[q];
        let expect = c0 * (1.0 - fz[q]) + c1 * fz[q];
        assert_eq!(
            out[q].to_bits(),
            expect.to_bits(),
            "trilinear mismatch at {q}"
        );
    }
    // Gather kernels drive their own tier twins.
    let mut out_s = vec![0.0f64; m];
    trilinear_k::run_scalar(&g, &fx, &fy, &fz, &ix, &iy, &iz, nx, nxy, &mut out_s);
    let mut out_w = vec![0.0f64; m];
    trilinear_k::run_lanes::<4>(&g, &fx, &fy, &fz, &ix, &iy, &iz, nx, nxy, &mut out_w);
    assert!(
        out_s
            .iter()
            .zip(&out_w)
            .all(|(a, b)| a.to_bits() == b.to_bits()),
        "trilinear lanes diverge from scalar"
    );
    log("trilinear", "pass", &trilinear_k::META.descr());
}

#[test]
fn dot_reduction_deterministic_and_tier_equal() {
    let n = 10_007;
    let (x, y) = (rand_vec(n, 7), rand_vec(n, 8));
    let d1 = dot_k::run(&x, &y);
    let d2 = dot_k::run_scalar(&x, &y);
    let d4 = dot_k::run_lanes::<4>(&x, &y);
    let d8 = dot_k::run_lanes::<8>(&x, &y);
    assert_eq!(d1.to_bits(), d2.to_bits(), "resolved-tier vs scalar");
    assert_eq!(d2.to_bits(), d4.to_bits(), "lanes 4 vs scalar");
    assert_eq!(d2.to_bits(), d8.to_bits(), "lanes 8 vs scalar");
    // Against the runtime's fixed-shape reference combiner.
    let prods: Vec<f64> = x.iter().zip(&y).map(|(a, b)| a * b).collect();
    assert_eq!(
        d2.to_bits(),
        fs_tilelang::deterministic_sum(&prods).to_bits(),
        "kernel reduction must equal the fixed-shape reference"
    );
    // Sanity vs a naive fold (envelope, not bitwise).
    let naive: f64 = prods.iter().sum();
    assert!((d2 - naive).abs() < 1e-9 * naive.abs().max(1.0));
    assert_eq!(dot_k::META.reduction, ReductionKind::DeterministicSum);
    log("dot", "pass", &dot_k::META.descr());
}

#[test]
fn metadata_feeds_the_roofline_table() {
    // The per-kernel variant/intensity table (P6 evidence, logged).
    for meta in [
        axpy_k::META,
        stencil3_k::META,
        stencil7_k::META,
        trilinear_k::META,
        dot_k::META,
    ] {
        assert!(meta.flops_per_elem > 0, "{}: zero flops counted", meta.name);
        assert!(meta.bytes_per_elem > 0);
        assert!(meta.intensity() > 0.0);
        log("roofline-meta", "info", &meta.descr());
    }
}
