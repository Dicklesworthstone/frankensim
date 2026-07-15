//! NEON lane capsule for the frozen Duct axial-z BGK collision.
//!
//! Two vector lanes are two independent cells. Separate multiply/add/subtract
//! and divide intrinsics preserve the scalar twin's existing operation tree.
#![allow(unsafe_code)] // registered capsule; SAFETY.md is adjacent

use core::arch::aarch64::{
    vaddq_f64, vdivq_f64, vdupq_n_f64, vld1q_f64, vmulq_f64, vst1q_f64, vsubq_f64,
};

use super::{MacroscopicTile, TileInput, TileOutput, validate_output};
use crate::d3q19::{CollisionError3, E3, Q3, TILE_CELLS, W3};

pub(super) fn selected_kernel(
    input: &TileInput<'_>,
    output: &mut TileOutput<'_>,
    macros: &MacroscopicTile,
    tau: f64,
    gz: f64,
) -> Result<(), CollisionError3> {
    assert!(
        TILE_CELLS.is_multiple_of(2),
        "D3Q19 NEON tile extent must be a multiple of two"
    );
    // SAFETY: NEON is architectural on aarch64; fixed tile references prove
    // every two-lane load/store extent.
    unsafe { kernel_neon(input, output, macros, tau, gz) };
    validate_output(output)
}

unsafe fn kernel_neon(
    input: &TileInput<'_>,
    output: &mut TileOutput<'_>,
    macros: &MacroscopicTile,
    tau: f64,
    gz: f64,
) {
    const LANES: usize = 2;
    // SAFETY: every offset is a multiple of two below TILE_CELLS=64. Each
    // load/store therefore covers one complete in-bounds chunk.
    unsafe {
        let one = vdupq_n_f64(1.0);
        let three = vdupq_n_f64(3.0);
        let four_point_five = vdupq_n_f64(4.5);
        let one_point_five = vdupq_n_f64(1.5);
        let nine = vdupq_n_f64(9.0);
        let tau_v = vdupq_n_f64(tau);
        let gz_v = vdupq_n_f64(gz);
        let coefficient = vdupq_n_f64(1.0 - 0.5 / tau);

        for offset in (0..TILE_CELLS).step_by(LANES) {
            let rho = vld1q_f64(macros.rho.as_ptr().add(offset));
            let ux = vld1q_f64(macros.velocity[0].as_ptr().add(offset));
            let uy = vld1q_f64(macros.velocity[1].as_ptr().add(offset));
            let uz = vld1q_f64(macros.velocity[2].as_ptr().add(offset));
            let usq = vaddq_f64(
                vaddq_f64(vmulq_f64(ux, ux), vmulq_f64(uy, uy)),
                vmulq_f64(uz, uz),
            );

            for direction in 0..Q3 {
                let population = vld1q_f64(input[direction].as_ptr().add(offset));
                let ex = vdupq_n_f64(f64::from(E3[direction].0));
                let ey = vdupq_n_f64(f64::from(E3[direction].1));
                let ez = vdupq_n_f64(f64::from(E3[direction].2));
                let weight = vdupq_n_f64(W3[direction]);
                let eu = vaddq_f64(
                    vaddq_f64(vmulq_f64(ex, ux), vmulq_f64(ey, uy)),
                    vmulq_f64(ez, uz),
                );
                let polynomial = vsubq_f64(
                    vaddq_f64(
                        vaddq_f64(one, vmulq_f64(three, eu)),
                        vmulq_f64(vmulq_f64(four_point_five, eu), eu),
                    ),
                    vmulq_f64(one_point_five, usq),
                );
                let equilibrium = vmulq_f64(vmulq_f64(weight, rho), polynomial);
                let force_projection = vaddq_f64(
                    vmulq_f64(three, vsubq_f64(ez, uz)),
                    vmulq_f64(vmulq_f64(nine, eu), ez),
                );
                let forcing = vmulq_f64(
                    vmulq_f64(vmulq_f64(coefficient, weight), force_projection),
                    gz_v,
                );
                let relaxed = vdivq_f64(vsubq_f64(equilibrium, population), tau_v);
                let value = vaddq_f64(vaddq_f64(population, relaxed), forcing);
                vst1q_f64(output[direction].as_mut_ptr().add(offset), value);
            }
        }
    }
}
