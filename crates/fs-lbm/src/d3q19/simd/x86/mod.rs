//! AVX2 lane capsule for the frozen Duct axial-z BGK collision.
//!
//! Four vector lanes are four independent cells. Every per-cell expression
//! uses explicit multiply/add/subtract/divide operations in the scalar twin's
//! order; FMA is admitted as part of the house x86 tier but is not used to
//! refreeze the existing BGK arithmetic.
#![allow(unsafe_code)] // registered capsule; SAFETY.md is adjacent

use core::arch::x86_64::{
    _mm256_add_pd, _mm256_div_pd, _mm256_loadu_pd, _mm256_mul_pd, _mm256_set1_pd, _mm256_storeu_pd,
    _mm256_sub_pd,
};

use super::{
    D3q19BgkSimdTier, MacroscopicTile, TileInput, TileKernel, TileOutput, scalar_kernel,
    selected_x86_tier_for, validate_output,
};
use crate::d3q19::{CollisionError3, E3, Q3, TILE_CELLS, W3};

pub(super) fn select_kernel() -> (TileKernel, D3q19BgkSimdTier) {
    let avx2_available = std::arch::is_x86_feature_detected!("avx2");
    let fma_available = std::arch::is_x86_feature_detected!("fma");
    let tier = selected_x86_tier_for(avx2_available, fma_available);
    if tier == D3q19BgkSimdTier::Avx2 {
        (selected_kernel, tier)
    } else {
        (scalar_kernel, D3q19BgkSimdTier::Scalar)
    }
}

fn selected_kernel(
    input: &TileInput<'_>,
    output: &mut TileOutput<'_>,
    macros: &MacroscopicTile,
    tau: f64,
    gz: f64,
) -> Result<(), CollisionError3> {
    assert!(
        TILE_CELLS.is_multiple_of(4),
        "D3Q19 AVX2 tile extent must be a multiple of four"
    );
    // SAFETY: this private thunk can enter the process-wide table only through
    // `select_kernel`, immediately after AVX2+FMA detection.
    unsafe { kernel_256(input, output, macros, tau, gz) };
    validate_output(output)
}

/// # Safety
/// The caller must establish AVX2+FMA availability. Tile references are
/// fixed-size, non-aliasing Rust borrows; unaligned loads/stores stay within
/// each four-element chunk.
#[target_feature(enable = "avx2,fma")]
unsafe fn kernel_256(
    input: &TileInput<'_>,
    output: &mut TileOutput<'_>,
    macros: &MacroscopicTile,
    tau: f64,
    gz: f64,
) {
    const LANES: usize = 4;
    // SAFETY: every offset is a multiple of four below TILE_CELLS=64. Each
    // load/store therefore covers one complete in-bounds chunk.
    unsafe {
        let one = _mm256_set1_pd(1.0);
        let three = _mm256_set1_pd(3.0);
        let four_point_five = _mm256_set1_pd(4.5);
        let one_point_five = _mm256_set1_pd(1.5);
        let nine = _mm256_set1_pd(9.0);
        let tau_v = _mm256_set1_pd(tau);
        let gz_v = _mm256_set1_pd(gz);
        let coefficient = _mm256_set1_pd(1.0 - 0.5 / tau);

        for offset in (0..TILE_CELLS).step_by(LANES) {
            let rho = _mm256_loadu_pd(macros.rho.as_ptr().add(offset));
            let ux = _mm256_loadu_pd(macros.velocity[0].as_ptr().add(offset));
            let uy = _mm256_loadu_pd(macros.velocity[1].as_ptr().add(offset));
            let uz = _mm256_loadu_pd(macros.velocity[2].as_ptr().add(offset));
            let usq = _mm256_add_pd(
                _mm256_add_pd(_mm256_mul_pd(ux, ux), _mm256_mul_pd(uy, uy)),
                _mm256_mul_pd(uz, uz),
            );

            for direction in 0..Q3 {
                let population = _mm256_loadu_pd(input[direction].as_ptr().add(offset));
                let ex = _mm256_set1_pd(f64::from(E3[direction].0));
                let ey = _mm256_set1_pd(f64::from(E3[direction].1));
                let ez = _mm256_set1_pd(f64::from(E3[direction].2));
                let weight = _mm256_set1_pd(W3[direction]);
                let eu = _mm256_add_pd(
                    _mm256_add_pd(_mm256_mul_pd(ex, ux), _mm256_mul_pd(ey, uy)),
                    _mm256_mul_pd(ez, uz),
                );
                let polynomial = _mm256_sub_pd(
                    _mm256_add_pd(
                        _mm256_add_pd(one, _mm256_mul_pd(three, eu)),
                        _mm256_mul_pd(_mm256_mul_pd(four_point_five, eu), eu),
                    ),
                    _mm256_mul_pd(one_point_five, usq),
                );
                let equilibrium = _mm256_mul_pd(_mm256_mul_pd(weight, rho), polynomial);
                let force_projection = _mm256_add_pd(
                    _mm256_mul_pd(three, _mm256_sub_pd(ez, uz)),
                    _mm256_mul_pd(_mm256_mul_pd(nine, eu), ez),
                );
                let forcing = _mm256_mul_pd(
                    _mm256_mul_pd(_mm256_mul_pd(coefficient, weight), force_projection),
                    gz_v,
                );
                let relaxed = _mm256_div_pd(_mm256_sub_pd(equilibrium, population), tau_v);
                let value = _mm256_add_pd(_mm256_add_pd(population, relaxed), forcing);
                _mm256_storeu_pd(output[direction].as_mut_ptr().add(offset), value);
            }
        }
    }
}
