//! AVX2 row-move capsule for the Duct pull-stream stencil.
//!
//! One vector is one four-cell x-row. X-moving directions permute the current
//! source row and blend one neighboring-tile or bounce-back edge value; all
//! other directions are a direct load/store. Population bits never enter an
//! arithmetic instruction.
#![allow(unsafe_code)] // registered capsule; SAFETY.md is adjacent

use core::arch::x86_64::{
    _mm256_blend_pd, _mm256_loadu_pd, _mm256_permute4x64_pd, _mm256_set1_pd, _mm256_storeu_pd,
};

use super::{
    D3q19StreamSimdTier, StreamKernel, pull_source, row_lane, scalar_kernel, selected_x86_tier_for,
    tile_index,
};
use crate::d3q19::{E3, OPP3, Q3, TILE, TILE_CELLS, Tile};

pub(super) fn select_kernel() -> (StreamKernel, D3q19StreamSimdTier) {
    let tier = selected_x86_tier_for(std::arch::is_x86_feature_detected!("avx2"));
    if tier == D3q19StreamSimdTier::Avx2 {
        (selected_kernel, tier)
    } else {
        (scalar_kernel, D3q19StreamSimdTier::Scalar)
    }
}

fn selected_kernel(
    post: &[Vec<Tile>; Q3],
    populations: &mut [Vec<Tile>; Q3],
    nx: usize,
    ny: usize,
    nz: usize,
) {
    assert_eq!(TILE, 4, "D3Q19 AVX2 stream requires four-cell tile rows");
    assert!(
        TILE_CELLS.is_multiple_of(4),
        "D3Q19 AVX2 stream tile extent must be a multiple of four"
    );
    // SAFETY: this private thunk enters the table only after AVX2 detection;
    // the safe outer façade already validated dimensions and field lengths.
    unsafe { kernel_256(post, populations, nx, ny, nz) };
}

/// # Safety
/// Requires AVX2 plus validated positive, tile-aligned dimensions and matching
/// input/output field lengths. The private selected thunk establishes both.
#[target_feature(enable = "avx2")]
unsafe fn kernel_256(
    post: &[Vec<Tile>; Q3],
    populations: &mut [Vec<Tile>; Q3],
    nx: usize,
    ny: usize,
    nz: usize,
) {
    let (ntx, nty) = (nx / TILE, ny / TILE);
    // SAFETY: the outer façade proved every field has exactly ntx*nty*ntz
    // tiles. Mapping helpers keep tile/lane indices in those extents; each
    // intrinsic accesses exactly one four-element row.
    unsafe {
        for direction in 0..Q3 {
            let (ex, ey, ez) = E3[direction];
            for z in 0..nz {
                let sz = pull_source(z, nz, ez, true).expect("periodic z source");
                for y in 0..ny {
                    let sy = pull_source(y, ny, ey, false);
                    for tx in 0..ntx {
                        let dtile = tile_index(ntx, nty, tx, y / TILE, z / TILE);
                        let drow = row_lane(y, z);
                        let value = if let Some(sy) = sy {
                            let stile = tile_index(ntx, nty, tx, sy / TILE, sz / TILE);
                            let srow = row_lane(sy, sz);
                            let current =
                                _mm256_loadu_pd(post[direction][stile].0.as_ptr().add(srow));
                            match ex {
                                0 => current,
                                1 => {
                                    let edge = if tx == 0 {
                                        post[OPP3[direction]][dtile].0[drow]
                                    } else {
                                        let left =
                                            tile_index(ntx, nty, tx - 1, sy / TILE, sz / TILE);
                                        post[direction][left].0[srow + TILE - 1]
                                    };
                                    let shifted = _mm256_permute4x64_pd::<0x90>(current);
                                    _mm256_blend_pd::<0x1>(shifted, _mm256_set1_pd(edge))
                                }
                                -1 => {
                                    let edge = if tx + 1 == ntx {
                                        post[OPP3[direction]][dtile].0[drow + TILE - 1]
                                    } else {
                                        let right =
                                            tile_index(ntx, nty, tx + 1, sy / TILE, sz / TILE);
                                        post[direction][right].0[srow]
                                    };
                                    let shifted = _mm256_permute4x64_pd::<0xf9>(current);
                                    _mm256_blend_pd::<0x8>(shifted, _mm256_set1_pd(edge))
                                }
                                _ => unreachable!(),
                            }
                        } else {
                            _mm256_loadu_pd(post[OPP3[direction]][dtile].0.as_ptr().add(drow))
                        };
                        _mm256_storeu_pd(
                            populations[direction][dtile].0.as_mut_ptr().add(drow),
                            value,
                        );
                    }
                }
            }
        }
    }
}
