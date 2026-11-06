//! Isolated rectangular-grid convolution for the symmetric free operator.
//!
//! Padding each axis to at least `2*n-1` retains every possible cell difference
//! without periodic images in the original grid. This applies the canonical
//! lower-triangle kernel used by dense assembly, including its rounded samples.
//! Original physical opening and acceptance still belong to the direct kernel.

use std::sync::Arc;

use fs_exec::Cx;
use fs_fft::{C64, FftNd, FftNdScratch};

use super::{ComplianceKernel, FiniteGapPlasticError, checkpoint, invalid};

const TILE: usize = 64;
const MAX_PADDED_CELLS: usize = 16 * super::super::MAX_FINITE_GAP_CELLS;

/// Immutable cached spectrum and plans; clones share their numeric payload.
#[derive(Debug, Clone)]
pub(super) struct FftCompliance {
    cells_x: usize,
    cells_y: usize,
    padded_x: usize,
    padded_y: usize,
    self_compliance: f64,
    plan: Arc<FftNd>,
    spectrum: Arc<[C64]>,
}

impl PartialEq for FftCompliance {
    fn eq(&self, other: &Self) -> bool {
        // Plans are deterministic functions of the retained padded dimensions.
        self.cells_x == other.cells_x
            && self.cells_y == other.cells_y
            && self.padded_x == other.padded_x
            && self.padded_y == other.padded_y
            && self.self_compliance == other.self_compliance
            && self.spectrum == other.spectrum
    }
}

/// Per-solve reusable work; no allocation occurs in an operator application.
pub(super) struct FftComplianceScratch {
    work: Vec<C64>,
    pencils: FftNdScratch,
}

impl FftCompliance {
    pub(super) fn build(
        kernel: &ComplianceKernel,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteGapPlasticError> {
        checkpoint(cx)?;
        let count = kernel
            .cells_x
            .checked_mul(kernel.cells_y)
            .filter(|count| *count > 0 && *count <= super::super::MAX_FINITE_GAP_CELLS)
            .ok_or_else(|| invalid("bounded FFT contact grid"))?;
        let padded_axis = |n: usize| {
            n.checked_mul(2)
                .and_then(|value| value.checked_sub(1))
                .and_then(usize::checked_next_power_of_two)
        };
        let padded_x = padded_axis(kernel.cells_x)
            .ok_or_else(|| invalid("representable FFT contact width"))?;
        let padded_y = padded_axis(kernel.cells_y)
            .ok_or_else(|| invalid("representable FFT contact depth"))?;
        let padded_count = padded_x
            .checked_mul(padded_y)
            .filter(|total| *total <= MAX_PADDED_CELLS)
            .ok_or_else(|| invalid("bounded FFT contact padding"))?;
        let rows = 2 * kernel.cells_y - 1;
        if kernel.stride != 2 * kernel.cells_x - 1 || kernel.values.len() != kernel.stride * rows {
            return Err(invalid("rectangular FFT compliance cache"));
        }
        let self_compliance = kernel.compliance(0, 0);
        if !self_compliance.is_finite() || self_compliance <= 0.0 {
            return Err(invalid("finite positive FFT self compliance"));
        }
        // Every original axis is <= count, so every padded axis is <=32768
        // under the contact grid cap. Existing Fft therefore stays below its
        // optional 65536-element six-step route and allocates no hidden work.
        debug_assert!(padded_x < 4 * count && padded_y < 4 * count);
        debug_assert!(
            padded_x <= 2 * super::super::MAX_FINITE_GAP_CELLS
                && padded_y <= 2 * super::super::MAX_FINITE_GAP_CELLS
        );
        let plan = Arc::new(FftNd::new(&[padded_y, padded_x]));
        checkpoint(cx)?;
        let mut spectrum = vec![C64::default(); padded_count];
        for row in 0..rows {
            checkpoint(cx)?;
            let dy = row as isize - (kernel.cells_y - 1) as isize;
            let y = if dy < 0 {
                padded_y - dy.unsigned_abs()
            } else {
                dy as usize
            };
            for column in 0..kernel.stride {
                if column % TILE == 0 {
                    checkpoint(cx)?;
                }
                let dx = column as isize - (kernel.cells_x - 1) as isize;
                let x = if dx < 0 {
                    padded_x - dx.unsigned_abs()
                } else {
                    dx as usize
                };
                // Dense assembly chooses the larger row-major cell index as
                // receiver. Its source-minus-receiver offset is always in
                // this negative half-plane, not an independent abs(x),abs(y).
                let (kx, ky) = if dy < 0 || (dy == 0 && dx <= 0) {
                    (dx, dy)
                } else {
                    (-dx, -dy)
                };
                let source_x = (kx + (kernel.cells_x - 1) as isize) as usize;
                let source_y = (ky + (kernel.cells_y - 1) as isize) as usize;
                let value = kernel.values[source_y * kernel.stride + source_x] / self_compliance;
                if !value.is_finite() {
                    return Err(invalid("finite normalized FFT compliance"));
                }
                spectrum[y * padded_x + x].re = value;
            }
        }
        let mut pencils = plan.scratch();
        plan.forward_with_scratch_cancellable(&mut spectrum, &mut pencils, || checkpoint(cx))?;
        if spectrum
            .iter()
            .any(|value| !value.re.is_finite() || !value.im.is_finite())
        {
            return Err(invalid("finite FFT compliance spectrum"));
        }
        checkpoint(cx)?;
        Ok(Self {
            cells_x: kernel.cells_x,
            cells_y: kernel.cells_y,
            padded_x,
            padded_y,
            self_compliance,
            plan,
            spectrum: spectrum.into(),
        })
    }

    pub(super) fn scratch(&self) -> FftComplianceScratch {
        FftComplianceScratch {
            work: vec![C64::default(); self.spectrum.len()],
            pencils: self.plan.scratch(),
        }
    }

    /// Apply `(C / divisor) input` on the original full grid.
    ///
    /// The input may have either sign, including CG search directions. The
    /// canonical kernel and input are normalized before spectral arithmetic;
    /// caller-owned free gather/scatter and spring terms remain separate.
    /// Scratch may be overwritten on refusal, while output publishes only
    /// after all finite-value and cancellation checks have passed.
    pub(super) fn apply_scaled(
        &self,
        input: &[f64],
        output: &mut [f64],
        divisor: f64,
        scratch: &mut FftComplianceScratch,
        cx: &Cx<'_>,
    ) -> Result<(), FiniteGapPlasticError> {
        checkpoint(cx)?;
        let count = self.cells_x * self.cells_y;
        if input.len() != count
            || output.len() != count
            || scratch.work.len() != self.spectrum.len()
        {
            return Err(invalid("FFT compliance vector dimensions"));
        }
        let relative_scale = self.self_compliance / divisor;
        if !divisor.is_finite()
            || divisor <= 0.0
            || !relative_scale.is_finite()
            || relative_scale <= 0.0
        {
            return Err(invalid("representable positive FFT compliance scale"));
        }
        let mut input_scale = 0.0_f64;
        for chunk in input.chunks(TILE) {
            checkpoint(cx)?;
            for value in chunk {
                if !value.is_finite() {
                    return Err(invalid("finite FFT compliance input"));
                }
                input_scale = input_scale.max(value.abs());
            }
        }
        if input_scale == 0.0 {
            checkpoint(cx)?;
            output.fill(0.0);
            return Ok(());
        }
        for chunk in scratch.work.chunks_mut(TILE) {
            checkpoint(cx)?;
            chunk.fill(C64::default());
        }
        for (row, source) in input.chunks_exact(self.cells_x).enumerate() {
            checkpoint(cx)?;
            for (column, value) in source.iter().enumerate() {
                scratch.work[row * self.padded_x + column].re = value / input_scale;
            }
        }
        self.plan.forward_with_scratch_cancellable(
            &mut scratch.work,
            &mut scratch.pencils,
            || checkpoint(cx),
        )?;
        for (values, coefficients) in scratch
            .work
            .chunks_mut(TILE)
            .zip(self.spectrum.chunks(TILE))
        {
            checkpoint(cx)?;
            for (value, coefficient) in values.iter_mut().zip(coefficients) {
                let input = *value;
                value.re = input
                    .re
                    .mul_add(coefficient.re, -(input.im * coefficient.im));
                value.im = input.re.mul_add(coefficient.im, input.im * coefficient.re);
                if !value.re.is_finite() || !value.im.is_finite() {
                    return Err(invalid("finite FFT compliance product"));
                }
            }
        }
        self.plan.inverse_with_scratch_cancellable(
            &mut scratch.work,
            &mut scratch.pencils,
            || checkpoint(cx),
        )?;
        for row in 0..self.cells_y {
            checkpoint(cx)?;
            for column in 0..self.cells_x {
                let value = &mut scratch.work[row * self.padded_x + column];
                value.re = (value.re * relative_scale) * input_scale;
                if !value.re.is_finite() || !value.im.is_finite() {
                    return Err(invalid("finite FFT compliance opening"));
                }
            }
        }
        checkpoint(cx)?;
        // A bounded original-grid copy is the only publication. There are no
        // remaining fallible operations or incomplete outputs on this path.
        for (row, destination) in output.chunks_exact_mut(self.cells_x).enumerate() {
            for (column, value) in destination.iter_mut().enumerate() {
                *value = scratch.work[row * self.padded_x + column].re;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    use super::super::super::{FiniteGapContactError, FiniteGapGrid};

    fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
        let gate = CancelGate::new_clock_free();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 47,
                    kernel_id: 131,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx, &gate)
        });
    }

    fn kernel(nx: usize, ny: usize) -> ComplianceKernel {
        ComplianceKernel::build(
            FiniteGapGrid {
                cells_x: nx,
                cells_y: ny,
                cell_width_m: 7.0e-4,
                cell_depth_m: 1.1e-3,
            },
            1.3e9,
        )
    }

    fn direct(kernel: &ComplianceKernel, input: &[f64], divisor: f64) -> Vec<f64> {
        (0..input.len())
            .map(|i| {
                input.iter().enumerate().fold(0.0, |sum, (j, value)| {
                    (kernel.compliance(i.max(j), i.min(j)) / divisor).mul_add(*value, sum)
                })
            })
            .collect()
    }

    fn assert_close(actual: &[f64], expected: &[f64], scale: f64) {
        for (i, (a, b)) in actual.iter().zip(expected).enumerate() {
            assert!(
                (a - b).abs() <= 5.0e-12 * scale,
                "cell {i}: {a:.16e} versus {b:.16e}"
            );
        }
    }

    #[test]
    fn isolated_edge_impulses_match_canonical_kernel_on_rectangular_grids() {
        with_cx(|cx, _| {
            for (nx, ny) in [(1, 1), (1, 7), (9, 1), (5, 7), (10, 3)] {
                let kernel = kernel(nx, ny);
                let operator = FftCompliance::build(&kernel, cx).unwrap();
                let mut scratch = operator.scratch();
                let mut input = vec![0.0; nx * ny];
                let mut output = vec![0.0; input.len()];
                let divisor = kernel.compliance(0, 0);
                for source in [0, nx - 1, (ny - 1) * nx, nx * ny - 1] {
                    input.fill(0.0);
                    input[source] = 1.0;
                    operator
                        .apply_scaled(&input, &mut output, divisor, &mut scratch, cx)
                        .unwrap();
                    assert_close(&output, &direct(&kernel, &input, divisor), 1.0);
                }
            }
        });
    }

    #[test]
    fn signed_krylov_vectors_preserve_rectangular_coupling_scaling_and_replay() {
        with_cx(|cx, _| {
            let kernel = kernel(7, 5);
            let operator = FftCompliance::build(&kernel, cx).unwrap();
            let mut scratch = operator.scratch();
            let base: Vec<f64> = (0..35).map(|i| ((i * 7 % 11) as f64 - 5.0) / 5.0).collect();
            let divisor = 2.5 * kernel.compliance(0, 0);
            for scale in [1.0e-120, 1.0, 1.0e120] {
                let input: Vec<f64> = base.iter().map(|value| value * scale).collect();
                let mut output = vec![0.0; input.len()];
                operator
                    .apply_scaled(&input, &mut output, divisor, &mut scratch, cx)
                    .unwrap();
                assert_close(&output, &direct(&kernel, &input, divisor), scale);
                let first = output.clone();
                operator
                    .apply_scaled(&input, &mut output, divisor, &mut scratch, cx)
                    .unwrap();
                assert_eq!(output, first);
            }
            let mut output = vec![1.0; 35];
            operator
                .apply_scaled(&[0.0; 35], &mut output, divisor, &mut scratch, cx)
                .unwrap();
            assert_eq!(output, vec![0.0; 35]);
            assert_eq!(operator, operator.clone());
        });
    }

    #[test]
    fn numerical_and_cancellation_refusals_preserve_output() {
        with_cx(|cx, gate| {
            let kernel = kernel(5, 3);
            let operator = FftCompliance::build(&kernel, cx).unwrap();
            let mut scratch = operator.scratch();
            let mut input = vec![1.0; 15];
            let mut output = vec![7.0; 15];
            let saved = output.clone();
            input[14] = f64::NAN;
            assert!(
                operator
                    .apply_scaled(&input, &mut output, 1.0, &mut scratch, cx)
                    .is_err()
            );
            assert_eq!(output, saved);
            input[14] = 1.0;
            assert!(
                operator
                    .apply_scaled(&input, &mut output, 0.0, &mut scratch, cx)
                    .is_err()
            );
            assert_eq!(output, saved);
            let invalid_shape = ComplianceKernel {
                cells_x: usize::MAX,
                cells_y: 2,
                stride: 0,
                values: vec![],
            };
            assert!(FftCompliance::build(&invalid_shape, cx).is_err());
            gate.request();
            assert!(matches!(
                operator.apply_scaled(&input, &mut output, 1.0, &mut scratch, cx),
                Err(FiniteGapPlasticError::Contact(
                    FiniteGapContactError::Cancelled
                ))
            ));
            assert_eq!(output, saved);
            assert!(matches!(
                FftCompliance::build(&kernel, cx),
                Err(FiniteGapPlasticError::Contact(
                    FiniteGapContactError::Cancelled
                ))
            ));
        });
    }
}
