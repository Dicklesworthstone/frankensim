//! NEON row-move capsule for the Duct pull-stream stencil.
//!
//! Two registers cover one four-cell x-row. `vext` shifts x-moving rows and
//! inserts one neighboring-tile or bounce-back edge value; other directions
//! are direct loads/stores. Population bits are moved without arithmetic.
#![allow(unsafe_code)] // registered capsule; SAFETY.md is adjacent

use core::arch::aarch64::{vdupq_n_f64, vextq_f64, vld1q_f64, vst1q_f64};

use super::{pull_source, row_lane, tile_index};
use crate::d3q19::{E3, OPP3, Q3, TILE, TILE_CELLS, Tile};

pub(super) fn selected_kernel(
    post: &[Vec<Tile>; Q3],
    populations: &mut [Vec<Tile>; Q3],
    nx: usize,
    ny: usize,
    nz: usize,
) {
    assert_eq!(TILE, 4, "D3Q19 NEON stream requires four-cell tile rows");
    assert!(
        TILE_CELLS.is_multiple_of(2),
        "D3Q19 NEON stream tile extent must be a multiple of two"
    );
    // SAFETY: NEON is architectural on aarch64; the safe outer façade already
    // validated dimensions and field lengths.
    unsafe { kernel_neon(post, populations, nx, ny, nz) };
}

unsafe fn kernel_neon(
    post: &[Vec<Tile>; Q3],
    populations: &mut [Vec<Tile>; Q3],
    nx: usize,
    ny: usize,
    nz: usize,
) {
    let (ntx, nty) = (nx / TILE, ny / TILE);
    // SAFETY: the outer façade proved every field has exactly ntx*nty*ntz
    // tiles. Mapping helpers keep tile/lane indices in those extents; each
    // intrinsic accesses exactly one two-element half-row.
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
                        let (low, high) = if let Some(sy) = sy {
                            let stile = tile_index(ntx, nty, tx, sy / TILE, sz / TILE);
                            let srow = row_lane(sy, sz);
                            let current_low =
                                vld1q_f64(post[direction][stile].0.as_ptr().add(srow));
                            let current_high =
                                vld1q_f64(post[direction][stile].0.as_ptr().add(srow + 2));
                            match ex {
                                0 => (current_low, current_high),
                                1 => {
                                    let edge = if tx == 0 {
                                        post[OPP3[direction]][dtile].0[drow]
                                    } else {
                                        let left =
                                            tile_index(ntx, nty, tx - 1, sy / TILE, sz / TILE);
                                        post[direction][left].0[srow + TILE - 1]
                                    };
                                    (
                                        vextq_f64::<1>(vdupq_n_f64(edge), current_low),
                                        vextq_f64::<1>(current_low, current_high),
                                    )
                                }
                                -1 => {
                                    let edge = if tx + 1 == ntx {
                                        post[OPP3[direction]][dtile].0[drow + TILE - 1]
                                    } else {
                                        let right =
                                            tile_index(ntx, nty, tx + 1, sy / TILE, sz / TILE);
                                        post[direction][right].0[srow]
                                    };
                                    (
                                        vextq_f64::<1>(current_low, current_high),
                                        vextq_f64::<1>(current_high, vdupq_n_f64(edge)),
                                    )
                                }
                                _ => unreachable!(),
                            }
                        } else {
                            (
                                vld1q_f64(post[OPP3[direction]][dtile].0.as_ptr().add(drow)),
                                vld1q_f64(post[OPP3[direction]][dtile].0.as_ptr().add(drow + 2)),
                            )
                        };
                        let destination = populations[direction][dtile].0.as_mut_ptr().add(drow);
                        vst1q_f64(destination, low);
                        vst1q_f64(destination.add(2), high);
                    }
                }
            }
        }
    }
}
