//! Reusable, cancellable pencils for the original separable N-D transform.

use crate::{C64, FftNd};

/// Reusable pencil and transform work for an [`FftNd`] plan.
///
/// The two buffers each hold the longest axis. They contain no accepted state
/// and may be reused after cancellation. The separable driver allocates nothing
/// during execution; the existing optional `frontier-sixstep` implementation
/// can still allocate inside a sufficiently large individual 1D transform.
#[derive(Debug, Clone)]
pub struct FftNdScratch {
    line: Vec<C64>,
    transform: Vec<C64>,
}

impl FftNd {
    /// Prepare both pencil buffers once for repeated transforms.
    #[must_use]
    pub fn scratch(&self) -> FftNdScratch {
        let longest = self.dims.iter().copied().max().unwrap_or(0);
        FftNdScratch {
            line: vec![C64::default(); longest],
            transform: vec![C64::default(); longest],
        }
    }

    /// Forward transform with reusable scratch and caller-owned interruption.
    ///
    /// Arithmetic, axis order and pencil order are identical to [`Self::forward`].
    /// The callback runs before and after every 1D pencil, and at entry/exit.
    /// A single 1D transform is the smallest cancellable unit; callers bound
    /// their axis lengths when they require bounded cancellation latency.
    ///
    /// # Errors
    /// Returns the callback's original error. On refusal, `data` and scratch
    /// may contain partial transforms; stage output when publication is atomic.
    ///
    /// # Panics
    /// If the data length differs from [`Self::total`] or scratch cannot hold
    /// the longest axis. Scratch from a larger plan can also be reused.
    pub fn forward_with_scratch_cancellable<E>(
        &self,
        data: &mut [C64],
        scratch: &mut FftNdScratch,
        checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<(), E> {
        self.run_with_scratch(data, scratch, false, checkpoint)
    }

    /// Inverse transform, normalized by `1/total`, with reusable pencil work.
    ///
    /// Cancellation and buffer semantics match
    /// [`Self::forward_with_scratch_cancellable`]. The original 1D inverse
    /// supplies normalization along each axis; no extra scale is applied.
    ///
    /// # Errors
    /// Returns the unchanged callback error without claiming complete output.
    ///
    /// # Panics
    /// For mismatched data or insufficient scratch, as in the forward method.
    pub fn inverse_with_scratch_cancellable<E>(
        &self,
        data: &mut [C64],
        scratch: &mut FftNdScratch,
        checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<(), E> {
        self.run_with_scratch(data, scratch, true, checkpoint)
    }

    pub(super) fn run_with_scratch<E>(
        &self,
        data: &mut [C64],
        workspace: &mut FftNdScratch,
        inverse: bool,
        mut checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<(), E> {
        let total = self.total();
        assert_eq!(
            data.len(),
            total,
            "buffer length {} must equal the product of dims {total}",
            data.len()
        );
        let longest = self.dims.iter().copied().max().unwrap_or(0);
        assert!(
            workspace.line.len() >= longest && workspace.transform.len() >= longest,
            "FFT pencil scratch must hold the longest axis {longest}"
        );
        checkpoint()?;
        for (ax, plan) in self.plans.iter().enumerate() {
            let n = self.dims[ax];
            if n == 1 {
                continue;
            }
            let stride: usize = self.dims[ax + 1..].iter().product();
            let outer: usize = self.dims[..ax].iter().product();
            let line = &mut workspace.line[..n];
            let scratch = &mut workspace.transform[..n];
            for o in 0..outer {
                for i in 0..stride {
                    checkpoint()?;
                    let base = o * n * stride + i;
                    for (t, slot) in line.iter_mut().enumerate() {
                        *slot = data[base + t * stride];
                    }
                    if inverse {
                        plan.inverse(line, scratch);
                    } else {
                        plan.forward(line, scratch);
                    }
                    checkpoint()?;
                    for (t, &value) in line.iter().enumerate() {
                        data[base + t * stride] = value;
                    }
                }
            }
        }
        checkpoint()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn input(rows: usize, columns: usize) -> Vec<C64> {
        (0..rows * columns)
            .map(|i| C64::new((i % 7) as f64 - 3.0, (i % 5) as f64 * 0.2 - 0.3))
            .collect()
    }

    fn direct_dft(input: &[C64], rows: usize, columns: usize) -> Vec<C64> {
        (0..input.len())
            .map(|k| {
                let mut sum = C64::default();
                for (j, value) in input.iter().enumerate() {
                    let turns = (k / columns * (j / columns)) as f64 / rows as f64
                        + (k % columns * (j % columns)) as f64 / columns as f64;
                    let (sine, cosine) = (-2.0 * core::f64::consts::PI * turns).sin_cos();
                    sum.re += value.re * cosine - value.im * sine;
                    sum.im += value.re * sine + value.im * cosine;
                }
                sum
            })
            .collect()
    }

    #[test]
    fn reusable_rectangular_pencils_match_direct_dft_and_round_trip() {
        for (rows, columns) in [(1, 1), (1, 16), (4, 8), (8, 2)] {
            let plan = FftNd::new(&[rows, columns]);
            let mut scratch = plan.scratch();
            let original = input(rows, columns);
            let expected = direct_dft(&original, rows, columns);
            let mut actual = original.clone();
            plan.forward_with_scratch_cancellable(&mut actual, &mut scratch, || Ok::<_, ()>(()))
                .unwrap();
            for (a, b) in actual.iter().zip(&expected) {
                assert!((a.re - b.re).abs() < 1.0e-10);
                assert!((a.im - b.im).abs() < 1.0e-10);
            }
            let first = actual.clone();
            plan.inverse_with_scratch_cancellable(&mut actual, &mut scratch, || Ok::<_, ()>(()))
                .unwrap();
            for (a, b) in actual.iter().zip(&original) {
                assert!((a.re - b.re).abs() < 2.0e-13);
                assert!((a.im - b.im).abs() < 2.0e-13);
            }
            actual.copy_from_slice(&original);
            plan.forward_with_scratch_cancellable(&mut actual, &mut scratch, || Ok::<_, ()>(()))
                .unwrap();
            assert_eq!(actual, first);
        }
    }

    #[test]
    fn interruption_keeps_its_reason_and_scratch_supports_retry() {
        let plan = FftNd::new(&[4, 8]);
        let mut scratch = plan.scratch();
        let original = input(4, 8);
        let mut actual = original.clone();
        let mut calls = 0;
        let result = plan.forward_with_scratch_cancellable(&mut actual, &mut scratch, || {
            calls += 1;
            if calls == 4 {
                Err("pencil budget")
            } else {
                Ok(())
            }
        });
        assert_eq!(result, Err("pencil budget"));
        assert_eq!(calls, 4);
        actual.copy_from_slice(&original);
        plan.forward_with_scratch_cancellable(&mut actual, &mut scratch, || Ok::<_, ()>(()))
            .unwrap();
        for (a, b) in actual.iter().zip(direct_dft(&original, 4, 8)) {
            assert!((a.re - b.re).abs() < 1.0e-10);
            assert!((a.im - b.im).abs() < 1.0e-10);
        }
    }
}
