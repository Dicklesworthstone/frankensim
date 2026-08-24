//! WASM SIMD128 Tier-1w capsule (bead `frankensim-wf-root-guzez.1.5`, E0.5).
//!
//! WebAssembly 128-bit vector (`v128`) SIMD intrinsics behind safe façades.
//! Scalar referee handles tails and fallback.
#![allow(unsafe_code)] // registered capsule — see SAFETY.md beside this file

#[cfg(all(target_arch = "wasm32", target_feature = "simd128", not(miri)))]
use core::arch::wasm32::{
    f64x2_add, f64x2_extract_lane, f64x2_mul, f64x2_splat, v128, v128_load, v128_store,
};

const LANES: usize = 2; // v128 contains 2 f64 lanes

/// y[i] = a * x[i] + y[i]
pub fn axpy(a: f64, x: &[f64], y: &mut [f64]) {
    assert_eq!(x.len(), y.len(), "axpy length mismatch");
    let (xc, xt) = x.as_chunks::<LANES>();
    let (yc, yt) = y.as_chunks_mut::<LANES>();

    #[cfg(all(target_arch = "wasm32", target_feature = "simd128", not(miri)))]
    unsafe {
        let va = f64x2_splat(a);
        for (xk, yk) in xc.iter().zip(yc) {
            let vx = v128_load(xk.as_ptr() as *const v128);
            let vy = v128_load(yk.as_ptr() as *const v128);
            let vprod = f64x2_mul(va, vx);
            let vres = f64x2_add(vprod, vy);
            v128_store(yk.as_mut_ptr() as *mut v128, vres);
        }
    }

    #[cfg(not(all(target_arch = "wasm32", target_feature = "simd128", not(miri))))]
    {
        for (xk, yk) in xc.iter().zip(yc) {
            yk[0] += a * xk[0];
            yk[1] += a * xk[1];
        }
    }

    crate::scalar::axpy(a, xt, yt);
}

/// x[i] *= a
pub fn scale(a: f64, x: &mut [f64]) {
    let (xc, xt) = x.as_chunks_mut::<LANES>();

    #[cfg(all(target_arch = "wasm32", target_feature = "simd128", not(miri)))]
    unsafe {
        let va = f64x2_splat(a);
        for xk in xc {
            let vx = v128_load(xk.as_ptr() as *const v128);
            let vres = f64x2_mul(va, vx);
            v128_store(xk.as_mut_ptr() as *mut v128, vres);
        }
    }

    #[cfg(not(all(target_arch = "wasm32", target_feature = "simd128", not(miri))))]
    {
        for xk in xc {
            xk[0] *= a;
            xk[1] *= a;
        }
    }

    crate::scalar::scale(a, xt);
}

/// out[i] = a[i] * b[i]
pub fn mul_elem(a: &[f64], b: &[f64], out: &mut [f64]) {
    assert_eq!(a.len(), b.len(), "mul_elem length mismatch");
    assert_eq!(a.len(), out.len(), "mul_elem length mismatch");
    let (ac, at) = a.as_chunks::<LANES>();
    let (bc, bt) = b.as_chunks::<LANES>();
    let (oc, ot) = out.as_chunks_mut::<LANES>();

    #[cfg(all(target_arch = "wasm32", target_feature = "simd128", not(miri)))]
    unsafe {
        for ((ak, bk), ok) in ac.iter().zip(bc).zip(oc) {
            let va = v128_load(ak.as_ptr() as *const v128);
            let vb = v128_load(bk.as_ptr() as *const v128);
            let vres = f64x2_mul(va, vb);
            v128_store(ok.as_mut_ptr() as *mut v128, vres);
        }
    }

    #[cfg(not(all(target_arch = "wasm32", target_feature = "simd128", not(miri))))]
    {
        for ((ak, bk), ok) in ac.iter().zip(bc).zip(oc) {
            ok[0] = ak[0] * bk[0];
            ok[1] = ak[1] * bk[1];
        }
    }

    crate::scalar::mul_elem(at, bt, ot);
}

/// out[i] = a[i] * b[i] + c[i]
pub fn fma3(a: &[f64], b: &[f64], c: &[f64], out: &mut [f64]) {
    assert_eq!(a.len(), b.len(), "fma3 length mismatch");
    assert_eq!(a.len(), c.len(), "fma3 length mismatch");
    assert_eq!(a.len(), out.len(), "fma3 length mismatch");
    let (ac, at) = a.as_chunks::<LANES>();
    let (bc, bt) = b.as_chunks::<LANES>();
    let (cc, ct) = c.as_chunks::<LANES>();
    let (oc, ot) = out.as_chunks_mut::<LANES>();

    #[cfg(all(target_arch = "wasm32", target_feature = "simd128", not(miri)))]
    unsafe {
        for (((ak, bk), ck), ok) in ac.iter().zip(bc).zip(cc).zip(oc) {
            let va = v128_load(ak.as_ptr() as *const v128);
            let vb = v128_load(bk.as_ptr() as *const v128);
            let vc = v128_load(ck.as_ptr() as *const v128);
            let vprod = f64x2_mul(va, vb);
            let vres = f64x2_add(vprod, vc);
            v128_store(ok.as_mut_ptr() as *mut v128, vres);
        }
    }

    #[cfg(not(all(target_arch = "wasm32", target_feature = "simd128", not(miri))))]
    {
        for (((ak, bk), ck), ok) in ac.iter().zip(bc).zip(cc).zip(oc) {
            ok[0] = ak[0].mul_add(bk[0], ck[0]);
            ok[1] = ak[1].mul_add(bk[1], ck[1]);
        }
    }

    crate::scalar::fma3(at, bt, ct, ot);
}

/// acc[i] = a[i] * b[i] + acc[i]
pub fn fmacc(a: &[f64], b: &[f64], acc: &mut [f64]) {
    assert_eq!(a.len(), b.len(), "fmacc length mismatch");
    assert_eq!(a.len(), acc.len(), "fmacc length mismatch");
    let (ac, at) = a.as_chunks::<LANES>();
    let (bc, bt) = b.as_chunks::<LANES>();
    let (oc, ot) = acc.as_chunks_mut::<LANES>();

    #[cfg(all(target_arch = "wasm32", target_feature = "simd128", not(miri)))]
    unsafe {
        for ((ak, bk), ok) in ac.iter().zip(bc).zip(oc) {
            let va = v128_load(ak.as_ptr() as *const v128);
            let vb = v128_load(bk.as_ptr() as *const v128);
            let vc = v128_load(ok.as_ptr() as *const v128);
            let vprod = f64x2_mul(va, vb);
            let vres = f64x2_add(vprod, vc);
            v128_store(ok.as_mut_ptr() as *mut v128, vres);
        }
    }

    #[cfg(not(all(target_arch = "wasm32", target_feature = "simd128", not(miri))))]
    {
        for ((ak, bk), ok) in ac.iter().zip(bc).zip(oc) {
            ok[0] = ak[0].mul_add(bk[0], ok[0]);
            ok[1] = ak[1].mul_add(bk[1], ok[1]);
        }
    }

    crate::scalar::fmacc(at, bt, ot);
}

/// dot product
pub fn dot(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(x.len(), y.len(), "dot length mismatch");
    let (xc, xt) = x.as_chunks::<LANES>();
    let (yc, yt) = y.as_chunks::<LANES>();

    let mut sum = 0.0f64;

    #[cfg(all(target_arch = "wasm32", target_feature = "simd128", not(miri)))]
    unsafe {
        let mut vsum = f64x2_splat(0.0);
        for (xk, yk) in xc.iter().zip(yc) {
            let vx = v128_load(xk.as_ptr() as *const v128);
            let vy = v128_load(yk.as_ptr() as *const v128);
            vsum = f64x2_add(vsum, f64x2_mul(vx, vy));
        }
        sum = f64x2_extract_lane::<0>(vsum) + f64x2_extract_lane::<1>(vsum);
    }

    #[cfg(not(all(target_arch = "wasm32", target_feature = "simd128", not(miri))))]
    {
        for (xk, yk) in xc.iter().zip(yc) {
            sum += xk[0] * yk[0] + xk[1] * yk[1];
        }
    }

    sum + crate::scalar::dot(xt, yt)
}

/// Sum elements
pub fn sum(x: &[f64]) -> f64 {
    let (xc, xt) = x.as_chunks::<LANES>();
    let mut total = 0.0f64;

    #[cfg(all(target_arch = "wasm32", target_feature = "simd128", not(miri)))]
    unsafe {
        let mut vsum = f64x2_splat(0.0);
        for xk in xc {
            let vx = v128_load(xk.as_ptr() as *const v128);
            vsum = f64x2_add(vsum, vx);
        }
        total = f64x2_extract_lane::<0>(vsum) + f64x2_extract_lane::<1>(vsum);
    }

    #[cfg(not(all(target_arch = "wasm32", target_feature = "simd128", not(miri))))]
    {
        for xk in xc {
            total += xk[0] + xk[1];
        }
    }

    total + crate::scalar::sum(xt)
}
