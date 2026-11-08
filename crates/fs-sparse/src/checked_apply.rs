//! Bounded cooperative sparse applies with staged output publication.
//!
//! These are local numerical operations. The checkpoint may poll a `Cx`, a
//! deadline, or a caller's work budget; it does not mint executor, allocation,
//! or invocation-finalization authority. The existing unchecked kernels keep
//! their allocation-free API and arithmetic.

use crate::Csr;
use std::fmt;

/// Maximum sparse traversal/arithmetic units between cooperative checkpoints.
///
/// A unit is a visited row or storage slot, one scalar fused multiply-add, or
/// one staged output write. BSR padding is charged even though it never reads
/// the RHS. This is a work bound, not a wall-clock latency guarantee.
pub const APPLY_CHECKPOINT_WORK: usize = 256;

/// A sparse apply refused before publishing any part of its output.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SparseApplyError<E> {
    /// A dense shape cannot be represented by `usize`.
    DimensionOverflow,
    /// An input or output slice has the wrong length.
    DimensionMismatch {
        /// `"input"` or `"output"`.
        buffer: &'static str,
        /// Required scalar length.
        expected: usize,
        /// Supplied scalar length.
        actual: usize,
    },
    /// Storage for the staged output could not be reserved.
    Allocation,
    /// The caller's unchanged cancellation, deadline, or budget error.
    Interrupted(E),
}

impl<E: fmt::Display> fmt::Display for SparseApplyError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DimensionOverflow => f.write_str("sparse apply dimensions overflow usize"),
            Self::DimensionMismatch {
                buffer,
                expected,
                actual,
            } => {
                write!(
                    f,
                    "sparse apply {buffer} length: expected {expected}, got {actual}"
                )
            }
            Self::Allocation => f.write_str("sparse apply could not reserve its staged output"),
            Self::Interrupted(error) => write!(f, "sparse apply interrupted: {error}"),
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for SparseApplyError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Interrupted(error) => Some(error),
            _ => None,
        }
    }
}

pub(crate) struct WorkPoll<F> {
    checkpoint: F,
    remaining: usize,
}

impl<E, F: FnMut() -> Result<(), E>> WorkPoll<F> {
    pub(crate) fn checkpoint(&mut self) -> Result<(), E> {
        (self.checkpoint)()?;
        self.remaining = APPLY_CHECKPOINT_WORK;
        Ok(())
    }

    #[inline(always)]
    pub(crate) fn work(&mut self) -> Result<(), E> {
        if self.remaining == 0 {
            self.checkpoint()?;
        }
        self.remaining -= 1;
        Ok(())
    }
}

pub(crate) fn apply<E, F, K>(
    nrows: usize,
    ncols: usize,
    input_len: usize,
    nrhs: usize,
    output: &mut [f64],
    checkpoint: F,
    kernel: K,
) -> Result<(), SparseApplyError<E>>
where
    F: FnMut() -> Result<(), E>,
    K: FnOnce(&mut [f64], &mut WorkPoll<F>) -> Result<(), E>,
{
    let input_expected = ncols
        .checked_mul(nrhs)
        .ok_or(SparseApplyError::DimensionOverflow)?;
    let output_expected = nrows
        .checked_mul(nrhs)
        .ok_or(SparseApplyError::DimensionOverflow)?;
    for (buffer, expected, actual) in [
        ("input", input_expected, input_len),
        ("output", output_expected, output.len()),
    ] {
        if expected != actual {
            return Err(SparseApplyError::DimensionMismatch {
                buffer,
                expected,
                actual,
            });
        }
    }
    let mut poll = WorkPoll {
        checkpoint,
        remaining: 0,
    };
    poll.checkpoint().map_err(SparseApplyError::Interrupted)?;
    let mut staged = Vec::new();
    staged
        .try_reserve_exact(output_expected)
        .map_err(|_| SparseApplyError::Allocation)?;
    // Initialization must also be interruptible for a mostly empty matrix with
    // a large output. A single resize would move the unbounded gap elsewhere.
    while staged.len() < output_expected {
        poll.checkpoint().map_err(SparseApplyError::Interrupted)?;
        let count = APPLY_CHECKPOINT_WORK.min(output_expected - staged.len());
        staged.resize(staged.len() + count, 0.0);
    }
    poll.checkpoint().map_err(SparseApplyError::Interrupted)?;
    // In particular, a 0x0 matrix with an enormous RHS count has no work.
    // Do not walk its empty RHS tiles merely because the count is nonzero.
    if !staged.is_empty() {
        kernel(&mut staged, &mut poll).map_err(SparseApplyError::Interrupted)?;
    }
    // This is the publication boundary. A request observed here leaves the
    // destination intact. Once accepted, the copy completes without polling:
    // returning a refusal midway through it would expose a torn destination.
    poll.checkpoint().map_err(SparseApplyError::Interrupted)?;
    output.copy_from_slice(&staged);
    Ok(())
}

impl Csr {
    /// Compute `y = A x`, with bounded checkpoints and all-or-nothing output.
    ///
    /// The successful result has the same ascending-column FMA sequence as
    /// [`Self::spmv`], including explicit zeros and non-finite IEEE arithmetic.
    /// A refusal or unwinding callback panic leaves `y` unchanged. One extra
    /// output-sized allocation is used. Allocation itself and the final
    /// publication copy are not interruptible; cancellation observed before
    /// that copy starts refuses, and a later request applies to subsequent work.
    /// This is local output atomicity, not a root-invocation transaction.
    ///
    /// # Errors
    /// Returns a shape/allocation refusal or the checkpoint's original error.
    pub fn spmv_with_checkpoint<E>(
        &self,
        x: &[f64],
        y: &mut [f64],
        checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<(), SparseApplyError<E>> {
        self.spmm_with_checkpoint(x, 1, y, checkpoint)
    }

    /// Compute row-major `Y = A X` with the guarantees of
    /// [`Self::spmv_with_checkpoint`]. Zero RHS count accepts empty buffers.
    /// Checkpoints cover wide rows, RHS tails and zero-entry rows, not just
    /// row boundaries. The existing runtime FMA dispatch is retained.
    ///
    /// # Errors
    /// Returns a shape/overflow/allocation refusal or the checkpoint's error.
    pub fn spmm_with_checkpoint<E>(
        &self,
        x: &[f64],
        nrhs: usize,
        y: &mut [f64],
        checkpoint: impl FnMut() -> Result<(), E>,
    ) -> Result<(), SparseApplyError<E>> {
        apply(
            self.nrows(),
            self.ncols(),
            x.len(),
            nrhs,
            y,
            checkpoint,
            |out, poll| crate::fma::csr_checkpointed_dispatch(self, x, nrhs, out, poll),
        )
    }

    #[allow(clippy::inline_always)] // inherit the registered FMA capsule
    #[inline(always)]
    pub(crate) fn spmm_checkpointed_body<E, F: FnMut() -> Result<(), E>>(
        &self,
        x: &[f64],
        nrhs: usize,
        y: &mut [f64],
        poll: &mut WorkPoll<F>,
    ) -> Result<(), E> {
        const WIDTH: usize = 8;
        for first in (0..nrhs).step_by(WIDTH) {
            let width = WIDTH.min(nrhs - first);
            for row in 0..self.nrows() {
                poll.work()?;
                let mut sums = [0.0; WIDTH];
                let (columns, values) = self.row(row);
                for (&column, &value) in columns.iter().zip(values) {
                    poll.work()?;
                    let offset = column * nrhs + first;
                    for (sum, &rhs) in sums[..width].iter_mut().zip(&x[offset..offset + width]) {
                        poll.work()?;
                        *sum = value.mul_add(rhs, *sum);
                    }
                }
                let offset = row * nrhs + first;
                for (out, &sum) in y[offset..offset + width].iter_mut().zip(&sums[..width]) {
                    poll.work()?;
                    *out = sum;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Bsr, Coo, CsrCompact, Sell};
    use std::convert::Infallible;

    enum Format<'a> {
        Csr(&'a Csr),
        Bsr(Bsr),
        Sell(Sell),
        Compact(CsrCompact),
    }

    impl Format<'_> {
        #[inline(never)]
        fn portable(&self, x: &[f64], nrhs: usize, y: &mut [f64]) {
            let mut poll = WorkPoll {
                checkpoint: || Ok::<_, Infallible>(()),
                remaining: APPLY_CHECKPOINT_WORK,
            };
            match self {
                Self::Csr(a) => a.spmm_checkpointed_body(x, nrhs, y, &mut poll),
                Self::Bsr(a) => a.spmm_checkpointed_body(x, nrhs, y, &mut poll),
                Self::Sell(a) => a.spmm_checkpointed_body(x, nrhs, y, &mut poll),
                Self::Compact(a) => a.spmm_checkpointed_body(x, nrhs, y, &mut poll),
            }
            .unwrap();
        }

        fn apply<E>(
            &self,
            x: &[f64],
            nrhs: usize,
            y: &mut [f64],
            checkpoint: impl FnMut() -> Result<(), E>,
        ) -> Result<(), SparseApplyError<E>> {
            match self {
                Self::Csr(a) => a.spmm_with_checkpoint(x, nrhs, y, checkpoint),
                Self::Bsr(a) => a.spmm_with_checkpoint(x, nrhs, y, checkpoint),
                Self::Sell(a) => a.spmm_with_checkpoint(x, nrhs, y, checkpoint),
                Self::Compact(a) => a.spmm_with_checkpoint(x, nrhs, y, checkpoint),
            }
        }
    }

    fn formats(a: &Csr, block: (usize, usize)) -> [Format<'_>; 4] {
        [
            Format::Csr(a),
            Format::Bsr(Bsr::from_csr(a, block.0, block.1)),
            Format::Sell(Sell::from_csr(a, 4, 8).unwrap()),
            Format::Compact(CsrCompact::from_csr(a)),
        ]
    }

    fn bits(values: &[f64]) -> Vec<u64> {
        values.iter().map(|v| v.to_bits()).collect()
    }

    #[test]
    fn g0_g5_checkpointed_formats_preserve_fma_bits_and_rhs_tails() {
        let mut coo = Coo::new(6, 9);
        for (r, c, v) in [
            (0, 0, 1e100),
            (0, 1, -1e100),
            (0, 7, 3.0),
            (1, 1, -0.0),
            (1, 5, 1.25),
            (3, 2, -1e-300),
            (4, 6, -2.0),
            (5, 0, 0.0),
            (5, 7, 5.0),
        ] {
            coo.push(r, c, v);
        }
        let a = coo.assemble();
        let formats = formats(&a, (2, 3));
        for nrhs in [0, 1, 2, 7, 8, 9, 17] {
            let x: Vec<f64> = (0..9 * nrhs)
                .map(|i| {
                    if i / nrhs == 8 {
                        f64::NAN
                    } else if i / nrhs == 2 {
                        1e-100
                    } else {
                        0.25 * (i % 7 + 1) as f64
                    }
                })
                .collect();
            let mut expected = vec![9.0; 6 * nrhs];
            a.spmm(&x, nrhs, &mut expected);
            for format in &formats {
                let mut output = vec![-9.0; expected.len()];
                format
                    .apply(&x, nrhs, &mut output, || Ok::<_, Infallible>(()))
                    .unwrap();
                assert_eq!(bits(&output), bits(&expected));
                let mut portable = vec![-19.0; expected.len()];
                format.portable(&x, nrhs, &mut portable);
                assert_eq!(bits(&output), bits(&portable));
                assert!(
                    output.iter().all(|v| v.is_finite()),
                    "padding must not read NaN"
                );
            }
        }
    }

    #[test]
    fn g1_checkpointed_products_match_an_exact_integer_oracle() {
        let mut coo = Coo::new(4, 6);
        let entries = [
            (0, 0, 3),
            (0, 5, -2),
            (1, 2, 7),
            (2, 3, -4),
            (2, 4, 1),
            (3, 0, -1),
            (3, 1, 2),
        ];
        for &(r, c, v) in &entries {
            coo.push(r, c, f64::from(v));
        }
        let a = coo.assemble();
        let nrhs = 11;
        let x: Vec<f64> = (0..6 * nrhs).map(|i| (i as i32 - 20) as f64).collect();
        let mut oracle = vec![0_i64; 4 * nrhs];
        for &(r, c, v) in &entries {
            for j in 0..nrhs {
                oracle[r * nrhs + j] += i64::from(v) * (c * nrhs + j) as i64 - i64::from(v) * 20;
            }
        }
        for format in formats(&a, (2, 3)) {
            let mut actual = vec![0.0; oracle.len()];
            format
                .apply(&x, nrhs, &mut actual, || Ok::<_, Infallible>(()))
                .unwrap();
            assert_eq!(actual, oracle.iter().map(|&v| v as f64).collect::<Vec<_>>());
        }
    }

    #[test]
    fn g4_every_observed_wide_row_checkpoint_preserves_the_whole_destination() {
        let columns = APPLY_CHECKPOINT_WORK * 5;
        let a = Csr::from_parts(
            1,
            columns,
            vec![0, columns],
            (0..columns).collect(),
            vec![0.5; columns],
        );
        let x = vec![2.0; columns * 9];
        for format in formats(&a, (1, columns)) {
            let original = vec![
                -0.0,
                21.0,
                f64::from_bits(0x7ff8_0000_0000_0042),
                4.0,
                5.0,
                6.0,
                7.0,
                8.0,
                9.0,
            ];
            let mut complete = original.clone();
            let mut polls = 0;
            format
                .apply(&x, 9, &mut complete, || {
                    polls += 1;
                    Ok::<_, usize>(())
                })
                .unwrap();
            assert!(
                polls > columns / APPLY_CHECKPOINT_WORK,
                "a single wide row must poll internally"
            );
            assert_eq!(complete, vec![columns as f64; 9]);
            for stop_at in 1..=polls {
                let mut output = original.clone();
                let mut seen = 0;
                let result = format.apply(&x, 9, &mut output, || {
                    seen += 1;
                    if seen == stop_at {
                        Err(stop_at)
                    } else {
                        Ok(())
                    }
                });
                assert_eq!(result, Err(SparseApplyError::Interrupted(stop_at)));
                assert_eq!(seen, stop_at, "no further callbacks after refusal");
                assert_eq!(bits(&output), bits(&original));
            }
        }
    }

    #[test]
    fn g4_block_padding_and_empty_output_rows_are_checkpointed() {
        let columns = APPLY_CHECKPOINT_WORK * 8;
        let a = Csr::from_parts(1, columns, vec![0, 1], vec![0], vec![2.0]);
        let bsr = Bsr::from_csr(&a, 1, columns);
        let mut x = vec![f64::NAN; columns];
        x[0] = 3.0;
        let mut output = [17.0];
        let mut polls = 0;
        bsr.spmv_with_checkpoint(&x, &mut output, || {
            polls += 1;
            Ok::<_, ()>(())
        })
        .unwrap();
        assert_eq!(output, [6.0]);
        assert!(polls >= columns / APPLY_CHECKPOINT_WORK);

        let rows = APPLY_CHECKPOINT_WORK * 5;
        let empty = Csr::from_parts(rows, 0, vec![0; rows + 1], vec![], vec![]);
        for format in formats(&empty, (1, 1)) {
            let mut output = vec![7.0; rows];
            let mut seen = 0;
            let result = format.apply(&[], 1, &mut output, || {
                seen += 1;
                if seen == 3 {
                    Err("during staging")
                } else {
                    Ok(())
                }
            });
            assert_eq!(result, Err(SparseApplyError::Interrupted("during staging")));
            assert_eq!(output, vec![7.0; rows]);
            let mut polls = 0;
            format
                .apply(&[], 1, &mut output, || {
                    polls += 1;
                    Ok::<_, ()>(())
                })
                .unwrap();
            assert!(polls >= 2 * rows / APPLY_CHECKPOINT_WORK);
            assert!(output.iter().all(|v| v.to_bits() == 0));
        }
    }

    #[test]
    fn g0_invalid_shapes_and_overflow_refuse_before_callbacks_or_writes() {
        let a = Csr::identity(2);
        for format in formats(&a, (1, 1)) {
            let mut output = [11.0, 12.0];
            let mut called = false;
            let result = format.apply(&[1.0], 1, &mut output, || {
                called = true;
                Ok::<_, ()>(())
            });
            assert_eq!(
                result,
                Err(SparseApplyError::DimensionMismatch {
                    buffer: "input",
                    expected: 2,
                    actual: 1,
                })
            );
            assert!(!called);
            assert_eq!(output, [11.0, 12.0]);
            let result = format.apply(&[1.0, 2.0], 1, &mut output[..1], || {
                called = true;
                Ok::<_, ()>(())
            });
            assert_eq!(
                result,
                Err(SparseApplyError::DimensionMismatch {
                    buffer: "output",
                    expected: 2,
                    actual: 1,
                })
            );
            assert!(!called);
        }
        let enormous = Csr::from_parts(0, usize::MAX, vec![0], vec![], vec![]);
        assert_eq!(
            enormous.spmm_with_checkpoint(&[], 2, &mut [], || Ok::<_, ()>(())),
            Err(SparseApplyError::DimensionOverflow)
        );
        let empty = Csr::from_parts(0, 0, vec![0], vec![], vec![]);
        for format in formats(&empty, (1, 1)) {
            let mut polls = 0;
            format
                .apply(&[], usize::MAX, &mut [], || {
                    polls += 1;
                    Ok::<_, ()>(())
                })
                .unwrap();
            assert_eq!(polls, 3, "zero output must not traverse empty RHS tiles");
        }
    }

    #[test]
    fn g4_unwinding_callback_panic_never_publishes_staged_rows() {
        let a = Csr::identity(300);
        for format in formats(&a, (1, 1)) {
            let mut output = vec![19.0; 300];
            let original = output.clone();
            let mut polls = 0;
            let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                let _ = format.apply(&vec![2.0; 300], 1, &mut output, || {
                    polls += 1;
                    assert!(polls != 6, "injected callback panic during arithmetic");
                    Ok::<_, ()>(())
                });
            }));
            assert!(outcome.is_err());
            assert_eq!(output, original);
        }
    }

    #[test]
    fn g0_poll_work_bound_includes_the_exact_limit() {
        let mut count = 0;
        let mut poll = WorkPoll {
            checkpoint: || {
                count += 1;
                Ok::<_, ()>(())
            },
            remaining: 0,
        };
        for _ in 0..APPLY_CHECKPOINT_WORK {
            poll.work().unwrap();
        }
        assert_eq!(poll.remaining, 0);
        poll.work().unwrap();
        assert_eq!(poll.remaining, APPLY_CHECKPOINT_WORK - 1);
        drop(poll);
        assert_eq!(count, 2);
    }
}
