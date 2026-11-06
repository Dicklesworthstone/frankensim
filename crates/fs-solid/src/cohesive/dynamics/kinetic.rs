//! Storage change on the unchanged symmetric reduced mass.
//!
//! The mechanical owner supplies its original mass action; this leaf owns only
//! compensated endpoint arithmetic and reduction, never a replacement operator.

use fs_math::dd::Dd;

#[cfg(test)]
#[path = "kinetic_tests.rs"]
mod tests;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum KineticError {
    Dimension,
    Cancelled,
    NonFinite,
}

pub(super) fn change(
    before: &[f64],
    after: &[f64],
    mut apply_mass: impl FnMut(&[f64], &mut [f64]),
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, KineticError> {
    let mut poll = || {
        if cancelled() {
            Err(KineticError::Cancelled)
        } else {
            Ok(())
        }
    };
    poll()?;
    if before.len() != after.len() {
        return Err(KineticError::Dimension);
    }
    let mut delta = Vec::with_capacity(before.len());
    let mut has_tail = false;
    for (i, (&old, &new)) in before.iter().zip(after).enumerate() {
        if i % 256 == 0 {
            poll()?;
        }
        let difference = Dd::from_f64(finite(new)?) - Dd::from_f64(finite(old)?);
        delta.push(finite(difference.hi)?);
        has_tail |= finite(difference.lo)? != 0.0;
    }
    let mut action = vec![f64::NAN; delta.len()];
    poll()?;
    apply_mass(&delta, &mut action);
    poll()?;
    let mut leading = Dd::ZERO;
    let mut cross = Dd::ZERO;
    let mut trailing = Dd::ZERO;
    for (i, ((&old, &new), &momentum)) in before.iter().zip(after).zip(&action).enumerate() {
        if i % 256 == 0 {
            poll()?;
        }
        let mean = Dd::from_f64(0.5 * old) + Dd::from_f64(0.5 * new);
        let momentum = Dd::from_f64(finite(momentum)?);
        leading = leading + Dd::from_f64(mean.hi) * momentum;
        cross = cross + Dd::from_f64(mean.lo) * momentum;
    }
    if has_tail {
        for (i, ((value, &old), &new)) in delta.iter_mut().zip(before).zip(after).enumerate() {
            if i % 256 == 0 {
                poll()?;
            }
            *value = (Dd::from_f64(new) - Dd::from_f64(old)).lo;
        }
        for chunk in action.chunks_mut(256) {
            poll()?;
            chunk.fill(f64::NAN);
        }
        poll()?;
        apply_mass(&delta, &mut action);
        poll()?;
        for (i, ((&old, &new), &momentum)) in before.iter().zip(after).zip(&action).enumerate() {
            if i % 256 == 0 {
                poll()?;
            }
            let mean = Dd::from_f64(0.5 * old) + Dd::from_f64(0.5 * new);
            let momentum = Dd::from_f64(finite(momentum)?);
            cross = cross + Dd::from_f64(mean.hi) * momentum;
            trailing = trailing + Dd::from_f64(mean.lo) * momentum;
        }
    }
    poll()?;
    finite(((leading + cross) + trailing).to_f64())
}

fn finite(value: f64) -> Result<f64, KineticError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(KineticError::NonFinite)
    }
}
