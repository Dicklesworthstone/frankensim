//! Increments of a fixed symmetric quadratic form on actual endpoint vectors.
//!
//! The caller supplies the original linear action. This module retains the low
//! words of endpoint differences and means before that action and the reduction;
//! it never substitutes a diagonal, lumped or reconstructed operator. Symmetry,
//! operator accuracy and any positivity requirement remain with the caller.

use crate::dd::Dd;

/// Two independent terms of a fixed symmetric quadratic form.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuadraticIncrement {
    /// `(after + before)^T A (after - before) / 2`, the storage change.
    pub change: Dd,
    /// `(after - before)^T A (after - before) / 2`.
    ///
    /// For positive mass and backward Euler this is the numerical kinetic loss.
    /// It is not material dissipation or heat. Centered integration has no such
    /// loss, even though this independent quadratic is generally nonzero.
    pub increment_quadratic: Dd,
}

/// Refusal before publishing a complete quadratic increment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QuadraticError<E> {
    /// Endpoint vectors have different dimensions.
    Dimension,
    /// The caller requested cancellation.
    Cancelled,
    /// An endpoint, operator output or accumulated result is not finite.
    NonFinite,
    /// The original operator refused its action.
    Operator(E),
}

/// Evaluate a symmetric quadratic increment without subtracting rounded totals.
///
/// `apply` must overwrite every output entry with the same unchanged linear
/// operator applied to its input. It receives the caller's cancellation probe
/// so long operator actions can retain their own bounded checkpoints and errors.
/// The surrounding vector passes poll every 256 entries and around each action.
///
/// Two vector buffers are reused. One action suffices when all differences are
/// exact in binary64; a second applies the low difference words otherwise.
/// Double-double leading, cross and trailing reductions retain small exchanges
/// when large coordinate transfers cancel. The returned words can be combined
/// with independent work before rounding. This is not an interval enclosure or
/// a bound for an inaccurate operator, overflow or arbitrary ill-conditioning.
pub fn increment<E, C: FnMut() -> bool>(
    before: &[f64],
    after: &[f64],
    mut apply: impl FnMut(&[f64], &mut [f64], &mut C) -> Result<(), E>,
    cancelled: &mut C,
) -> Result<QuadraticIncrement, QuadraticError<E>> {
    poll(cancelled)?;
    if before.len() != after.len() {
        return Err(QuadraticError::Dimension);
    }
    let mut delta = Vec::with_capacity(before.len());
    let mut has_tail = false;
    for (i, (&old, &new)) in before.iter().zip(after).enumerate() {
        if i % 256 == 0 {
            poll(cancelled)?;
        }
        let difference = Dd::from_f64(finite(new)?) - Dd::from_f64(finite(old)?);
        delta.push(finite(difference.hi)?);
        has_tail |= finite(difference.lo)? != 0.0;
    }
    let mut action = vec![f64::NAN; delta.len()];
    poll(cancelled)?;
    apply(&delta, &mut action, cancelled).map_err(QuadraticError::Operator)?;
    poll(cancelled)?;
    let mut change = Layers::default();
    let mut increment_quadratic = Layers::default();
    for (i, ((&old, &new), &output)) in before.iter().zip(after).zip(&action).enumerate() {
        if i % 256 == 0 {
            poll(cancelled)?;
        }
        let mean = Dd::from_f64(0.5 * old) + Dd::from_f64(0.5 * new);
        let difference = Dd::from_f64(new) - Dd::from_f64(old);
        let output = Dd::from_f64(finite(output)?);
        change.leading = change.leading + Dd::from_f64(mean.hi) * output;
        change.cross = change.cross + Dd::from_f64(mean.lo) * output;
        increment_quadratic.leading =
            increment_quadratic.leading + Dd::from_f64(0.5 * difference.hi) * output;
        increment_quadratic.cross =
            increment_quadratic.cross + Dd::from_f64(0.5 * difference.lo) * output;
    }
    if has_tail {
        for (i, ((value, &old), &new)) in delta.iter_mut().zip(before).zip(after).enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            *value = (Dd::from_f64(new) - Dd::from_f64(old)).lo;
        }
        for chunk in action.chunks_mut(256) {
            poll(cancelled)?;
            chunk.fill(f64::NAN);
        }
        poll(cancelled)?;
        apply(&delta, &mut action, cancelled).map_err(QuadraticError::Operator)?;
        poll(cancelled)?;
        for (i, ((&old, &new), &output)) in before.iter().zip(after).zip(&action).enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            let mean = Dd::from_f64(0.5 * old) + Dd::from_f64(0.5 * new);
            let difference = Dd::from_f64(new) - Dd::from_f64(old);
            let output = Dd::from_f64(finite(output)?);
            change.cross = change.cross + Dd::from_f64(mean.hi) * output;
            change.trailing = change.trailing + Dd::from_f64(mean.lo) * output;
            increment_quadratic.cross =
                increment_quadratic.cross + Dd::from_f64(0.5 * difference.hi) * output;
            increment_quadratic.trailing =
                increment_quadratic.trailing + Dd::from_f64(0.5 * difference.lo) * output;
        }
    }
    poll(cancelled)?;
    Ok(QuadraticIncrement {
        change: change.finish()?,
        increment_quadratic: increment_quadratic.finish()?,
    })
}

#[derive(Default)]
struct Layers {
    leading: Dd,
    cross: Dd,
    trailing: Dd,
}

impl Layers {
    fn finish<E>(self) -> Result<Dd, QuadraticError<E>> {
        let value = (self.leading + self.cross) + self.trailing;
        finite(value.hi)?;
        finite(value.lo)?;
        Ok(value)
    }
}

fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), QuadraticError<E>> {
    if cancelled() {
        Err(QuadraticError::Cancelled)
    } else {
        Ok(())
    }
}

fn finite<E>(value: f64) -> Result<f64, QuadraticError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(QuadraticError::NonFinite)
    }
}

#[cfg(test)]
mod tests {
    use super::{QuadraticError, increment};
    use crate::dd::Dd;
    use std::convert::Infallible;

    #[test]
    fn coupled_countermotion_retains_storage_and_separate_euler_loss() {
        let common = 2.0_f64.powi(54);
        let before = [common, common];
        let after = [common + 4.0, common - 4.0];
        let mut calls = 0;
        let account = increment(
            &before,
            &after,
            |x, y, _| {
                calls += 1;
                y[0] = 2.0 * x[0] + x[1];
                y[1] = x[0] + 2.0 * x[1];
                Ok::<_, Infallible>(())
            },
            &mut || false,
        )
        .unwrap();
        assert_eq!(calls, 1);
        assert_eq!(account.change, Dd::from_f64(16.0));
        assert_eq!(account.increment_quadratic, Dd::from_f64(16.0));
        // The same physical endpoints have BE endpoint work 32 J and midpoint
        // work 16 J. The integration owner decides which numerical loss applies.
        assert_eq!(
            (account.change + account.increment_quadratic).to_f64(),
            32.0
        );
    }

    #[test]
    fn exchanged_inventory_retains_difference_tails_before_mass_action() {
        let large = 2.0_f64.powi(54);
        let before = [large, 1.0];
        let after = [2.0, large];
        for (old, new, expected) in [(&before, &after, 1.5), (&after, &before, -1.5)] {
            let mut calls = 0;
            let account = increment(
                old,
                new,
                |x, y, _| {
                    calls += 1;
                    y.copy_from_slice(x);
                    Ok::<_, Infallible>(())
                },
                &mut || false,
            )
            .unwrap();
            assert_eq!(calls, 2);
            assert_eq!(account.change, Dd::from_f64(expected));
            assert!(account.increment_quadratic.hi > 0.0);
        }
    }

    #[test]
    fn retained_low_word_closes_independent_work_before_rounding() {
        let before = [2.0_f64.powi(54), 1.0];
        let after = [before[0] + 4.0, 2.0];
        let account = increment(
            &before,
            &after,
            |x, y, _| {
                y.copy_from_slice(x);
                Ok::<_, Infallible>(())
            },
            &mut || false,
        )
        .unwrap();
        let large_work = Dd::from_f64(4.0 * before[0]) + Dd::from_f64(8.0);
        assert_eq!((account.change - large_work).to_f64(), 1.5);
        assert_ne!(account.change.lo, 0.0);
        assert_eq!(account.increment_quadratic, Dd::from_f64(8.5));
    }

    #[test]
    fn operator_failure_and_mid_action_cancellation_remain_distinct() {
        let failure = increment(&[1.0], &[2.0], |_, _, _| Err("mass refused"), &mut || false);
        assert_eq!(failure, Err(QuadraticError::Operator("mass refused")));
        let requested = std::cell::Cell::new(false);
        let mut callback_observed_cancel = false;
        let result = increment(
            &[1.0],
            &[2.0],
            |x, y, cancelled| {
                requested.set(true);
                callback_observed_cancel = cancelled();
                y.copy_from_slice(x);
                Ok::<_, Infallible>(())
            },
            &mut || requested.get(),
        );
        assert!(callback_observed_cancel);
        assert_eq!(result, Err(QuadraticError::Cancelled));
    }

    #[test]
    fn incomplete_operator_output_and_nonfinite_endpoints_refuse() {
        let incomplete = increment(
            &[1.0],
            &[2.0],
            |_, _, _| Ok::<_, Infallible>(()),
            &mut || false,
        );
        assert_eq!(incomplete, Err(QuadraticError::NonFinite));
        let nonfinite = increment(
            &[f64::NAN],
            &[2.0],
            |_, _, _| panic!("invalid endpoint reached the operator"),
            &mut || false,
        );
        assert_eq!(nonfinite, Err(QuadraticError::<Infallible>::NonFinite));
    }
}
